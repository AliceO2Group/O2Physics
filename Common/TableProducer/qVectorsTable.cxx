// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   qVectorsTable.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Task calculating the Q-vectors for each collision in a bunch crossing
///         (with or without corrections) and save the results in a dedicated table.
///

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DetectorsCommonDataFormats/AlignParam.h>
#include <FT0Base/Geometry.h>
#include <FV0Base/Geometry.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <TComplex.h>
#include <TH3.h>
#include <TMath.h>
#include <TProfile3D.h>
#include <TString.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected,
                               aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>;

struct qVectorsTable {
  enum Detectors {
    kFT0C = 0,
    kFT0A = 1,
    kFT0M,
    kFV0A,
    kTPCPos,
    kTPCNeg,
    kTPCAll,
    kNDetectors
  };
  enum Corrections {
    kNoCorr = 0,
    kRecenter,
    kTwist,
    kRescale,
    kNCorrections
  };
  enum MultNorms {
    kNoNorm = 0,
    kScalarProd,
    kEsE,
    kMultNormTypes
  };

  // Configurables.
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL",
                                     "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than",
                                      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
                                      "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  Configurable<int> cfgCentEsti{"cfgCentEsti",
                                2, "Centrality estimator (Run3): 0 = FT0M, 1 = FT0A, 2 = FT0C, 3 = FV0A"};

  Configurable<float> cfgMinPtOnTPC{"cfgMinPtOnTPC", 0.15, "minimum transverse momentum selection for TPC tracks participating in Q-vector reconstruction"};
  Configurable<float> cfgMaxPtOnTPC{"cfgMaxPtOnTPC", 5., "maximum transverse momentum selection for TPC tracks participating in Q-vector reconstruction"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgEtaMin{"cfgEtaMin", -0.8, "Minimum pseudorapidiy for charged track"};

  Configurable<int> cfgCorrLevel{"cfgCorrLevel", 4, "calibration step: 0 = no corr, 1 = gain corr, 2 = rectr, 3 = twist, 4 = full"};
  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2, 3}, "Modulation of interest"};
  Configurable<float> cfgMaxCentrality{"cfgMaxCentrality", 100.f, "max. centrality for Q vector calibration"};

  Configurable<bool> useCorrectionForRun{"useCorrectionForRun", true, "Get Qvector corrections based on run number instead of timestamp"};
  Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Users/j/junlee/Qvector/GainEq", "CCDB path for gain equalization constants"};
  Configurable<std::string> cfgQvecCalibPath{"cfgQvecCalibPath", "Analysis/EventPlane/QVecCorrections", "CCDB path for Q-vector calibration constants"};

  Configurable<bool> cfgShiftCorr{"cfgShiftCorr", false, "configurable flag for shift correction"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "", "CCDB path for shift correction"};

  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, ""};

  Configurable<bool> cfgUseFT0C{"cfgUseFT0C", false, "Initial value for using FT0C. By default obtained from DataModel."};
  Configurable<bool> cfgUseFT0A{"cfgUseFT0A", false, "Initial value for using FT0A. By default obtained from DataModel."};
  Configurable<bool> cfgUseFT0M{"cfgUseFT0M", false, "Initial value for using FT0M. By default obtained from DataModel."};
  Configurable<bool> cfgUseFV0A{"cfgUseFV0A", false, "Initial value for using FV0A. By default obtained from DataModel."};
  Configurable<bool> cfgUseTPCpos{"cfgUseTPCpos", false, "Initial value for using TPCPos. By default obtained from DataModel."};
  Configurable<bool> cfgUseTPCneg{"cfgUseTPCneg", false, "Initial value for using TPCNeg. By default obtained from DataModel."};
  Configurable<bool> cfgUseTPCall{"cfgUseTPCall", false, "Initial value for using TPCAll. By default obtained from DataModel."};
  Configurable<bool> cfgProduceRedQVecs{"cfgProduceRedQVecs", false, "Produce reduced Q-vectors for Event-Shape Engineering"};

  // Table.
  Produces<aod::Qvectors> qVector;
  Produces<aod::QvectorFT0Cs> qVectorFT0C;
  Produces<aod::QvectorFT0As> qVectorFT0A;
  Produces<aod::QvectorFT0Ms> qVectorFT0M;
  Produces<aod::QvectorFV0As> qVectorFV0A;
  Produces<aod::QvectorTPCposs> qVectorTPCPos;
  Produces<aod::QvectorTPCnegs> qVectorTPCNeg;
  Produces<aod::QvectorTPCalls> qVectorTPCAll;

  Produces<aod::QvectorFT0CVecs> qVectorFT0CVec;
  Produces<aod::QvectorFT0AVecs> qVectorFT0AVec;
  Produces<aod::QvectorFT0MVecs> qVectorFT0MVec;
  Produces<aod::QvectorFV0AVecs> qVectorFV0AVec;
  Produces<aod::QvectorTPCposVecs> qVectorTPCPosVec;
  Produces<aod::QvectorTPCnegVecs> qVectorTPCNegVec;
  Produces<aod::QvectorTPCallVecs> qVectorTPCAllVec;

  Produces<aod::EseQvectors> eseQVector;
  Produces<aod::EseQvecFT0Cs> eseQVectorFT0C;
  Produces<aod::EseQvecFT0As> eseQVectorFT0A;
  Produces<aod::EseQvecFT0Ms> eseQVectorFT0M;
  Produces<aod::EseQvecFV0As> eseQVectorFV0A;
  Produces<aod::EseQvecTPCposs> eseQVectorTPCPos;
  Produces<aod::EseQvecTPCnegs> eseQVectorTPCNeg;
  Produces<aod::EseQvecTPCalls> eseQVectorTPCAll;

  Produces<aod::EseQvecFT0CVecs> eseQVectorFT0CVec;
  Produces<aod::EseQvecFT0AVecs> eseQVectorFT0AVec;
  Produces<aod::EseQvecFT0MVecs> eseQVectorFT0MVec;
  Produces<aod::EseQvecFV0AVecs> eseQVectorFV0AVec;
  Produces<aod::EseQvecTPCposVecs> eseQVectorTPCPosVec;
  Produces<aod::EseQvecTPCnegVecs> eseQVectorTPCNegVec;
  Produces<aod::EseQvecTPCallVecs> eseQVectorTPCAllVec;
  Produces<aod::EseQvecPercs> eseQVectorPerc;

  std::vector<float> FT0RelGainConst{};
  std::vector<float> FV0RelGainConst{};

  // Enable access to the CCDB for the offset and correction constants and save them
  // in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // geometry instances for V0 and T0
  o2::fv0::Geometry* fv0geom;
  o2::ft0::Geometry ft0geom;

  // Variables for other classes.
  EventPlaneHelper helperEP;

  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  const float minAmplitude = 1.0e-8f;
  int runNumber{-1};
  float cent;

  std::vector<TH3F*> corrsQvecSp{};
  std::vector<TH3F*> corrsQvecEse{};
  std::vector<TProfile3D*> shiftProfileSp{};
  std::vector<TProfile3D*> shiftProfileEse{};

  // Deprecated, will be removed in future after transition time //
  Configurable<bool> cfgUseBPos{"cfgUseBPos", false, "Initial value for using BPos. By default obtained from DataModel."};
  Configurable<bool> cfgUseBNeg{"cfgUseBNeg", false, "Initial value for using BNeg. By default obtained from DataModel."};
  Configurable<bool> cfgUseBTot{"cfgUseBTot", false, "Initial value for using BTot. By default obtained from DataModel."};

  Produces<aod::QvectorBPoss> qVectorBPos;
  Produces<aod::QvectorBNegs> qVectorBNeg;
  Produces<aod::QvectorBTots> qVectorBTot;

  Produces<aod::QvectorBPosVecs> qVectorBPosVec;
  Produces<aod::QvectorBNegVecs> qVectorBNegVec;
  Produces<aod::QvectorBTotVecs> qVectorBTotVec;
  /////////////////////////////////////////////////////////////////

  std::unordered_map<std::string, bool> useDetector = {
    {"QvectorBTots", cfgUseBTot},
    {"QvectorBNegs", cfgUseBNeg},
    {"QvectorBPoss", cfgUseBPos},
    {"QvectorTPCalls", cfgUseTPCall},
    {"QvectorTPCnegs", cfgUseTPCneg},
    {"QvectorTPCposs", cfgUseTPCpos},
    {"QvectorFV0As", cfgUseFV0A},
    {"QvectorFT0Ms", cfgUseFT0M},
    {"QvectorFT0As", cfgUseFT0A},
    {"QvectorFT0Cs", cfgUseFT0C}};

  void init(InitContext& initContext)
  {
    // Check the sub-detector used
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      for (auto const& input : device.inputs) {
        if (input.matcher.binding == "Qvectors") {
          for (auto const& det : useDetector) {
            useDetector[det.first.data()] = true;
          }
          LOGF(info, "Using all detectors.");
          goto allDetectorsInUse; // Added to break from nested loop if all detectors are in use.
        }
        for (auto const& det : useDetector) {
          std::string tableNameWithVector = det.first; // for replacing s with Vecs at the end.
          if (input.matcher.binding == det.first || input.matcher.binding == tableNameWithVector.replace(tableNameWithVector.size() - 1, 1, "Vecs")) {
            useDetector[det.first.data()] = true;
            LOGF(info, Form("Using detector: %s.", det.first.data()));
          }
        }
      }
    }

  // Exit point in case all detectors are being used.
  allDetectorsInUse:
    // Setup the access to the CCDB objects of interest.
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(cfgCcdbParam.nolaterthan.value);
    ccdb->setFatalWhenNull(false);

    AxisSpec axisPt = {40, 0.0, 4.0};
    AxisSpec axisEta = {32, -0.8, 0.8};
    AxisSpec axisPhi = {32, 0, constants::math::TwoPI};
    AxisSpec axixCent = {20, 0, 100};

    AxisSpec axisFITamp{cfgaxisFITamp, "FIT amp"};
    AxisSpec axisChID = {220, 0, 220};

    fv0geom = o2::fv0::Geometry::instance(o2::fv0::Geometry::eUninitialized);

    histosQA.add("ChTracks", "", {HistType::kTHnSparseF, {axisPt, axisEta, axisPhi, axixCent}});
    histosQA.add("FT0Amp", "", {HistType::kTH2F, {axisFITamp, axisChID}});
    histosQA.add("FT0AmpCor", "", {HistType::kTH2F, {axisFITamp, axisChID}});
    histosQA.add("FV0Amp", "", {HistType::kTH2F, {axisFITamp, axisChID}});
    histosQA.add("FV0AmpCor", "", {HistType::kTH2F, {axisFITamp, axisChID}});
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    FT0RelGainConst.clear();
    FV0RelGainConst.clear();
    FT0RelGainConst = {};
    FV0RelGainConst = {};

    std::string fullPath;

    auto timestamp = bc.timestamp();
    auto runnumber = bc.runNumber();

    auto offsetFT0 = getForTsOrRun<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp, runnumber);
    auto offsetFV0 = getForTsOrRun<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", timestamp, runnumber);

    if (offsetFT0 != nullptr) {
      helperEP.SetOffsetFT0A((*offsetFT0)[0].getX(), (*offsetFT0)[0].getY());
      helperEP.SetOffsetFT0C((*offsetFT0)[1].getX(), (*offsetFT0)[1].getY());
    } else {
      LOGF(fatal, "Could not get the alignment parameters for FT0.");
    }

    if (offsetFV0 != nullptr) {
      helperEP.SetOffsetFV0left((*offsetFV0)[0].getX(), (*offsetFV0)[0].getY());
      helperEP.SetOffsetFV0right((*offsetFV0)[1].getX(), (*offsetFV0)[1].getY());
    } else {
      LOGF(fatal, "Could not get the alignment parameters for FV0.");
    }

    corrsQvecSp.clear();
    for (std::size_t i = 0; i < cfgnMods->size(); i++) {
      int ind = cfgnMods->at(i);
      fullPath = cfgQvecCalibPath;
      fullPath += "/v";
      fullPath += std::to_string(ind);
      auto modeCorrQvecSp = getForTsOrRun<TH3F>(fullPath, timestamp, runnumber);
      if (!modeCorrQvecSp) {
        fullPath = cfgQvecCalibPath;
        fullPath += "/v2";
        modeCorrQvecSp = getForTsOrRun<TH3F>(fullPath, timestamp, runnumber);
      }
      corrsQvecSp.push_back(modeCorrQvecSp);
    }

    corrsQvecEse.clear();
    for (std::size_t i = 0; i < cfgnMods->size(); i++) {
      int ind = cfgnMods->at(i);
      fullPath = cfgQvecCalibPath;
      fullPath += "/eseq";
      fullPath += std::to_string(ind);
      auto modeCorrQvecEse = getForTsOrRun<TH3F>(fullPath, timestamp, runnumber);
      if (!modeCorrQvecEse) {
        fullPath = cfgQvecCalibPath;
        fullPath += "/eseq2";
        modeCorrQvecEse = getForTsOrRun<TH3F>(fullPath, timestamp, runnumber);
      }
      corrsQvecEse.push_back(modeCorrQvecEse);
    }

    if (cfgShiftCorr) {
      shiftProfileSp.clear();
      for (std::size_t i = 0; i < cfgnMods->size(); i++) {
        int ind = cfgnMods->at(i);
        fullPath = cfgShiftPath;
        fullPath += "/v";
        fullPath += std::to_string(ind);
        auto objshift = getForTsOrRun<TProfile3D>(fullPath, timestamp, runnumber);
        shiftProfileSp.push_back(objshift);
      }

      if (cfgProduceRedQVecs) {
        shiftProfileEse.clear();
        for (std::size_t i = 0; i < cfgnMods->size(); i++) {
          int ind = cfgnMods->at(i);
          fullPath = cfgShiftPath;
          fullPath += "/eseq";
          fullPath += std::to_string(ind);
          auto objshift = getForTsOrRun<TProfile3D>(fullPath, timestamp, runnumber);
          shiftProfileEse.push_back(objshift);
        }
      }
    }

    fullPath = cfgGainEqPath;
    fullPath += "/FT0";
    const int nPixelsFT0 = 208;
    const auto objft0Gain = getForTsOrRun<std::vector<float>>(fullPath, timestamp, runnumber);
    if (!objft0Gain || cfgCorrLevel == 0) {
      for (auto i{0u}; i < nPixelsFT0; i++) {
        FT0RelGainConst.push_back(1.);
      }
    } else {
      FT0RelGainConst = *(objft0Gain);
    }

    fullPath = cfgGainEqPath;
    fullPath += "/FV0";
    const int nChannelsFV0 = 48;
    const auto objfv0Gain = getForTsOrRun<std::vector<float>>(fullPath, timestamp, runnumber);
    if (!objfv0Gain || cfgCorrLevel == 0) {
      for (auto i{0u}; i < nChannelsFV0; i++) {
        FV0RelGainConst.push_back(1.);
      }
    } else {
      FV0RelGainConst = *(objfv0Gain);
    }
  }

  template <typename TrackType>
  bool selTrack(const TrackType track)
  {
    if (track.pt() < cfgMinPtOnTPC)
      return false;
    if (track.pt() > cfgMaxPtOnTPC)
      return false;
    if (!track.passedITSNCls())
      return false;
    if (!track.passedITSChi2NDF())
      return false;
    if (!track.passedITSHits())
      return false;
    if (!track.passedTPCCrossedRowsOverNCls())
      return false;
    if (!track.passedTPCChi2NDF())
      return false;
    if (!track.passedDCAxy())
      return false;
    if (!track.passedDCAz())
      return false;

    return true;
  }

  /// Function to get corrections from CCDB eithr using the timestamp or the runnumber
  /// \param fullPath is the path to correction in CCDB
  /// \param timestamp is the collision timestamp
  /// \param runNumber is the collision run number
  /// \return CCDB correction
  template <typename CorrectionType>
  CorrectionType* getForTsOrRun(std::string const& fullPath, int64_t timestamp, int runNumber)
  {
    if (useCorrectionForRun) {
      return ccdb->getForRun<CorrectionType>(fullPath, runNumber);
    } else {
      return ccdb->getForTimeStamp<CorrectionType>(fullPath, timestamp);
    }
  }

  /// Function to normalize the q-vectors
  /// \param qVecReNorm is the vector with the normalized real part of the q-vector for each detector and correction step
  /// \param qVecImNorm is the vector with the normalized imaginary part of the q-vector for each detector and correction step
  /// \param qVecReRaw is the vector with the raw real part of the q-vector for each detector and correction step
  /// \param qVecImRaw is the vector with the raw imaginary part of the q-vector for each detector and correction step
  /// \param qVecAmp is the vector with the amplitude of the q-vector for each detector and correction step
  /// \param normType is the type of normalization to apply to the q-vectors
  void normalizeQVec(std::vector<float>& qVecReNorm,
                     std::vector<float>& qVecImNorm,
                     std::vector<float> const& qVecReRaw,
                     std::vector<float> const& qVecImRaw,
                     std::vector<float> const& qVecAmp,
                     MultNorms normType)
  {
    for (std::size_t i = 0; i < kNDetectors; i++) {
      float qVecDetReNorm{-999.}, qVecDetImNorm{-999.};
      if (qVecAmp[i] > minAmplitude) {
        switch (normType) {
          case MultNorms::kScalarProd:
            qVecDetReNorm = qVecReRaw[i] / qVecAmp[i];
            qVecDetImNorm = qVecImRaw[i] / qVecAmp[i];
            break;
          case MultNorms::kEsE:
            qVecDetReNorm = qVecReRaw[i] / std::sqrt(qVecAmp[i]);
            qVecDetImNorm = qVecImRaw[i] / std::sqrt(qVecAmp[i]);
            break;
          default:
            LOGP(fatal, "Undefined normalization type for Q-vector amplitude. Check the configuration.");
            break;
        }
      }
      for (int iCorr = 0; iCorr < kNCorrections; iCorr++) {
        qVecReNorm.push_back(qVecDetReNorm);
        qVecImNorm.push_back(qVecDetImNorm);
      }
    }
  }

  /// Function to calculate the un-normalized q-vectors
  /// \param cent is the collision centrality
  /// \param qVecRe is the vector with the real part of the q-vector for each detector and correction step
  /// \param qVecIm is the vector with the imaginary part of the q-vector for each detector and correction step
  /// \param histsCorrs is the vector with the histograms with the correction constants for each detector and correction step
  /// \param nMode is the modulation of interest
  void correctQVec(float cent, std::vector<float>& qVecRe, std::vector<float>& qVecIm, TH3F* histsCorrs, std::vector<TProfile3D*>& shiftProfile, int nMode)
  {
    int nCorrections = static_cast<int>(kNCorrections);
    if (cent < cfgMaxCentrality) {
      for (auto i{0u}; i < kTPCAll + 1; i++) {
        int idxDet = i * kNCorrections;
        helperEP.DoRecenter(qVecRe[idxDet + kRecenter], qVecIm[idxDet + kRecenter],
                            histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));

        helperEP.DoRecenter(qVecRe[idxDet + kTwist], qVecIm[idxDet + kTwist],
                            histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));
        helperEP.DoTwist(qVecRe[idxDet + kTwist], qVecIm[idxDet + kTwist],
                         histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 3, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 4, i + 1));

        helperEP.DoRecenter(qVecRe[idxDet + kRescale], qVecIm[idxDet + kRescale],
                            histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));
        helperEP.DoTwist(qVecRe[idxDet + kRescale], qVecIm[idxDet + kRescale],
                         histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 3, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 4, i + 1));
        helperEP.DoRescale(qVecRe[idxDet + kRescale], qVecIm[idxDet + kRescale],
                           histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 5, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 6, i + 1));
      }
      if (cfgShiftCorr) {
        auto deltaPsiFT0C = 0.0;
        auto deltaPsiFT0A = 0.0;
        auto deltaPsiFT0M = 0.0;
        auto deltaPsiFV0A = 0.0;
        auto deltaPsiTPCPos = 0.0;
        auto deltaPsiTPCNeg = 0.0;
        auto deltaPsiTPCAll = 0.0;

        auto psiDefFT0C = std::atan2(qVecIm[kFT0C * nCorrections + kRescale], qVecRe[kFT0C * nCorrections + kRescale]) / static_cast<float>(nMode);
        auto psiDefFT0A = std::atan2(qVecIm[kFT0A * nCorrections + kRescale], qVecRe[kFT0A * nCorrections + kRescale]) / static_cast<float>(nMode);
        auto psiDefFT0M = std::atan2(qVecIm[kFT0M * nCorrections + kRescale], qVecRe[kFT0M * nCorrections + kRescale]) / static_cast<float>(nMode);
        auto psiDefFV0A = std::atan2(qVecIm[kFV0A * nCorrections + kRescale], qVecRe[kFV0A * nCorrections + kRescale]) / static_cast<float>(nMode);
        auto psiDefTPCPos = std::atan2(qVecIm[kTPCPos * nCorrections + kRescale], qVecRe[kTPCPos * nCorrections + kRescale]) / static_cast<float>(nMode);
        auto psiDefTPCNeg = std::atan2(qVecIm[kTPCNeg * nCorrections + kRescale], qVecRe[kTPCNeg * nCorrections + kRescale]) / static_cast<float>(nMode);
        auto psiDefTPCAll = std::atan2(qVecIm[kTPCAll * nCorrections + kRescale], qVecRe[kTPCAll * nCorrections + kRescale]) / static_cast<float>(nMode);

        for (int iShift = 1; iShift <= 10; iShift++) {
          auto coeffShiftXFT0C = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kFT0C, iShift - 0.5));
          auto coeffShiftYFT0C = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kFT0C + 1, iShift - 0.5));
          auto coeffShiftXFT0A = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kFT0A, iShift - 0.5));
          auto coeffShiftYFT0A = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kFT0A + 1, iShift - 0.5));
          auto coeffShiftXFT0M = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kFT0M, iShift - 0.5));
          auto coeffShiftYFT0M = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kFT0M + 1, iShift - 0.5));
          auto coeffShiftXFV0A = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kFV0A, iShift - 0.5));
          auto coeffShiftYFV0A = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kFV0A + 1, iShift - 0.5));
          auto coeffShiftXTPCPos = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kTPCPos, iShift - 0.5));
          auto coeffShiftYTPCPos = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kTPCPos + 1, iShift - 0.5));
          auto coeffShiftXTPCNeg = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kTPCNeg, iShift - 0.5));
          auto coeffShiftYTPCNeg = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kTPCNeg + 1, iShift - 0.5));
          auto coeffShiftXTPCAll = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kTPCAll, iShift - 0.5));
          auto coeffShiftYTPCAll = shiftProfile.at(nMode - 2)->GetBinContent(shiftProfile.at(nMode - 2)->FindBin(cent, 2 * kTPCAll + 1, iShift - 0.5));

          deltaPsiFT0C += ((2. / (1.0 * iShift)) * (-coeffShiftXFT0C * std::cos(iShift * static_cast<float>(nMode) * psiDefFT0C) + coeffShiftYFT0C * std::sin(iShift * static_cast<float>(nMode) * psiDefFT0C))) / static_cast<float>(nMode);
          deltaPsiFT0A += ((2. / (1.0 * iShift)) * (-coeffShiftXFT0A * std::cos(iShift * static_cast<float>(nMode) * psiDefFT0A) + coeffShiftYFT0A * std::sin(iShift * static_cast<float>(nMode) * psiDefFT0A))) / static_cast<float>(nMode);
          deltaPsiFT0M += ((2. / (1.0 * iShift)) * (-coeffShiftXFT0M * std::cos(iShift * static_cast<float>(nMode) * psiDefFT0M) + coeffShiftYFT0M * std::sin(iShift * static_cast<float>(nMode) * psiDefFT0M))) / static_cast<float>(nMode);
          deltaPsiFV0A += ((2. / (1.0 * iShift)) * (-coeffShiftXFV0A * std::cos(iShift * static_cast<float>(nMode) * psiDefFV0A) + coeffShiftYFV0A * std::sin(iShift * static_cast<float>(nMode) * psiDefFV0A))) / static_cast<float>(nMode);
          deltaPsiTPCPos += ((2. / (1.0 * iShift)) * (-coeffShiftXTPCPos * std::cos(iShift * static_cast<float>(nMode) * psiDefTPCPos) + coeffShiftYTPCPos * std::sin(iShift * static_cast<float>(nMode) * psiDefTPCPos))) / static_cast<float>(nMode);
          deltaPsiTPCNeg += ((2. / (1.0 * iShift)) * (-coeffShiftXTPCNeg * std::cos(iShift * static_cast<float>(nMode) * psiDefTPCNeg) + coeffShiftYTPCNeg * std::sin(iShift * static_cast<float>(nMode) * psiDefTPCNeg))) / static_cast<float>(nMode);
          deltaPsiTPCAll += ((2. / (1.0 * iShift)) * (-coeffShiftXTPCAll * std::cos(iShift * static_cast<float>(nMode) * psiDefTPCAll) + coeffShiftYTPCAll * std::sin(iShift * static_cast<float>(nMode) * psiDefTPCAll))) / static_cast<float>(nMode);
        }

        deltaPsiFT0C *= static_cast<float>(nMode);
        deltaPsiFT0A *= static_cast<float>(nMode);
        deltaPsiFT0M *= static_cast<float>(nMode);
        deltaPsiFV0A *= static_cast<float>(nMode);
        deltaPsiTPCPos *= static_cast<float>(nMode);
        deltaPsiTPCNeg *= static_cast<float>(nMode);
        deltaPsiTPCAll *= static_cast<float>(nMode);

        float qVecReShiftedFT0C = qVecRe[kFT0C * nCorrections + kRescale] * std::cos(deltaPsiFT0C) - qVecIm[kFT0C * nCorrections + kRescale] * std::sin(deltaPsiFT0C);
        float qVecImShiftedFT0C = qVecRe[kFT0C * nCorrections + kRescale] * std::sin(deltaPsiFT0C) + qVecIm[kFT0C * nCorrections + kRescale] * std::cos(deltaPsiFT0C);
        float qVecReShiftedFT0A = qVecRe[kFT0A * nCorrections + kRescale] * std::cos(deltaPsiFT0A) - qVecIm[kFT0A * nCorrections + kRescale] * std::sin(deltaPsiFT0A);
        float qVecImShiftedFT0A = qVecRe[kFT0A * nCorrections + kRescale] * std::sin(deltaPsiFT0A) + qVecIm[kFT0A * nCorrections + kRescale] * std::cos(deltaPsiFT0A);
        float qVecReShiftedFT0M = qVecRe[kFT0M * nCorrections + kRescale] * std::cos(deltaPsiFT0M) - qVecIm[kFT0M * nCorrections + kRescale] * std::sin(deltaPsiFT0M);
        float qVecImShiftedFT0M = qVecRe[kFT0M * nCorrections + kRescale] * std::sin(deltaPsiFT0M) + qVecIm[kFT0M * nCorrections + kRescale] * std::cos(deltaPsiFT0M);
        float qVecReShiftedFV0A = qVecRe[kFV0A * nCorrections + kRescale] * std::cos(deltaPsiFV0A) - qVecIm[kFV0A * nCorrections + kRescale] * std::sin(deltaPsiFV0A);
        float qVecImShiftedFV0A = qVecRe[kFV0A * nCorrections + kRescale] * std::sin(deltaPsiFV0A) + qVecIm[kFV0A * nCorrections + kRescale] * std::cos(deltaPsiFV0A);
        float qVecReShiftedTPCPos = qVecRe[kTPCPos * nCorrections + kRescale] * std::cos(deltaPsiTPCPos) - qVecIm[kTPCPos * nCorrections + kRescale] * std::sin(deltaPsiTPCPos);
        float qVecImShiftedTPCPos = qVecRe[kTPCPos * nCorrections + kRescale] * std::sin(deltaPsiTPCPos) + qVecIm[kTPCPos * nCorrections + kRescale] * std::cos(deltaPsiTPCPos);
        float qVecReShiftedTPCNeg = qVecRe[kTPCNeg * nCorrections + kRescale] * std::cos(deltaPsiTPCNeg) - qVecIm[kTPCNeg * nCorrections + kRescale] * std::sin(deltaPsiTPCNeg);
        float qVecImShiftedTPCNeg = qVecRe[kTPCNeg * nCorrections + kRescale] * std::sin(deltaPsiTPCNeg) + qVecIm[kTPCNeg * nCorrections + kRescale] * std::cos(deltaPsiTPCNeg);
        float qVecReShiftedTPCAll = qVecRe[kTPCAll * nCorrections + kRescale] * std::cos(deltaPsiTPCAll) - qVecIm[kTPCAll * nCorrections + kRescale] * std::sin(deltaPsiTPCAll);
        float qVecImShiftedTPCAll = qVecRe[kTPCAll * nCorrections + kRescale] * std::sin(deltaPsiTPCAll) + qVecIm[kTPCAll * nCorrections + kRescale] * std::cos(deltaPsiTPCAll);

        qVecRe[kFT0C * nCorrections + kRescale] = qVecReShiftedFT0C;
        qVecIm[kFT0C * nCorrections + kRescale] = qVecImShiftedFT0C;
        qVecRe[kFT0A * nCorrections + kRescale] = qVecReShiftedFT0A;
        qVecIm[kFT0A * nCorrections + kRescale] = qVecImShiftedFT0A;
        qVecRe[kFT0M * nCorrections + kRescale] = qVecReShiftedFT0M;
        qVecIm[kFT0M * nCorrections + kRescale] = qVecImShiftedFT0M;
        qVecRe[kFV0A * nCorrections + kRescale] = qVecReShiftedFV0A;
        qVecIm[kFV0A * nCorrections + kRescale] = qVecImShiftedFV0A;
        qVecRe[kTPCPos * nCorrections + kRescale] = qVecReShiftedTPCPos;
        qVecIm[kTPCPos * nCorrections + kRescale] = qVecImShiftedTPCPos;
        qVecRe[kTPCNeg * nCorrections + kRescale] = qVecReShiftedTPCNeg;
        qVecIm[kTPCNeg * nCorrections + kRescale] = qVecImShiftedTPCNeg;
        qVecRe[kTPCAll * nCorrections + kRescale] = qVecReShiftedTPCAll;
        qVecIm[kTPCAll * nCorrections + kRescale] = qVecImShiftedTPCAll;
      }
    }
  }

  /// Function to calculate the un-normalized q-vectors
  /// \param nMode is the harmonic number of the q-vector
  /// \param coll is the collision object
  /// \param tracks are the tracks associated to the collision
  /// \param qVecRe is the vector with the real part of the q-vector for each detector
  /// \param qVecIm is the vector with the imaginary part of the q-vector for each detector
  /// \param qVecAmp is the vector with the amplitude of the signal in each detector
  /// \param trkTPCPosLabel is the vector with the number of TPC tracks with positive eta
  /// \param trkTPCNegLabel is the vector with the number of TPC tracks with negative eta
  /// \param trkTPCAllLabel is the vector with the number of TPC tracks with any eta
  template <typename Nmode, typename CollType, typename TrackType>
  void calcQVec(const Nmode nMode, const CollType& coll, const TrackType& tracks, std::vector<float>& qVecRe, std::vector<float>& qVecIm, std::vector<float>& qVecAmp, std::vector<int>& trkTPCPosLabel, std::vector<int>& trkTPCNegLabel, std::vector<int>& trkTPCAllLabel)
  {
    float qVectFT0A[2] = {-999., -999.};
    float qVectFT0C[2] = {-999., -999.};
    float qVectFT0M[2] = {-999., -999.};
    float qVectFV0A[2] = {-999., -999.};
    float qVectTPCPos[2] = {0., 0.}; // Always computed
    float qVectTPCNeg[2] = {0., 0.}; // Always computed
    float qVectTPCAll[2] = {0., 0.}; // Always computed

    TComplex qVecDet(0);
    TComplex qVecFT0M(0);
    float sumAmplFT0A = 0.;
    float sumAmplFT0C = 0.;
    float sumAmplFT0M = 0.;
    float sumAmplFV0A = 0.;

    if (coll.has_foundFT0() && (useDetector["QvectorFT0As"] || useDetector["QvectorFT0Cs"] || useDetector["QvectorFT0Ms"])) {
      auto ft0 = coll.foundFT0();

      if (useDetector["QvectorFT0As"]) {
        for (std::size_t iChA = 0; iChA < ft0.channelA().size(); iChA++) {
          float ampl = ft0.amplitudeA()[iChA];
          int FT0AchId = ft0.channelA()[iChA];

          histosQA.fill(HIST("FT0Amp"), ampl, FT0AchId);
          histosQA.fill(HIST("FT0AmpCor"), ampl / FT0RelGainConst[FT0AchId], FT0AchId);

          helperEP.SumQvectors(0, FT0AchId, ampl / FT0RelGainConst[FT0AchId], nMode, qVecDet, sumAmplFT0A, ft0geom, fv0geom);
          helperEP.SumQvectors(0, FT0AchId, ampl / FT0RelGainConst[FT0AchId], nMode, qVecFT0M, sumAmplFT0M, ft0geom, fv0geom);
        }
        if (sumAmplFT0A > minAmplitude) {
          qVectFT0A[0] = qVecDet.Re();
          qVectFT0A[1] = qVecDet.Im();
        }
      }

      if (useDetector["QvectorFT0Cs"]) {
        qVecDet = TComplex(0., 0.);
        for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
          float ampl = ft0.amplitudeC()[iChC];
          int FT0CchId = ft0.channelC()[iChC] + 96;

          histosQA.fill(HIST("FT0Amp"), ampl, FT0CchId);
          histosQA.fill(HIST("FT0AmpCor"), ampl / FT0RelGainConst[FT0CchId], FT0CchId);

          helperEP.SumQvectors(0, FT0CchId, ampl / FT0RelGainConst[FT0CchId], nMode, qVecDet, sumAmplFT0C, ft0geom, fv0geom);
          helperEP.SumQvectors(0, FT0CchId, ampl / FT0RelGainConst[FT0CchId], nMode, qVecFT0M, sumAmplFT0M, ft0geom, fv0geom);
        }

        if (sumAmplFT0C > minAmplitude) {
          qVectFT0C[0] = qVecDet.Re();
          qVectFT0C[1] = qVecDet.Im();
        }
        if (sumAmplFT0M > minAmplitude && useDetector["QvectorFT0Ms"]) {
          qVectFT0M[0] = qVecFT0M.Re();
          qVectFT0M[1] = qVecFT0M.Im();
        }
      }

      qVecDet = TComplex(0., 0.);
      sumAmplFV0A = 0;
      if (coll.has_foundFV0() && useDetector["QvectorFV0As"]) {
        auto fv0 = coll.foundFV0();

        for (std::size_t iCh = 0; iCh < fv0.channel().size(); iCh++) {
          float ampl = fv0.amplitude()[iCh];
          int FV0AchId = fv0.channel()[iCh];
          histosQA.fill(HIST("FV0Amp"), ampl, FV0AchId);
          histosQA.fill(HIST("FV0AmpCor"), ampl / FV0RelGainConst[FV0AchId], FV0AchId);

          helperEP.SumQvectors(1, FV0AchId, ampl / FV0RelGainConst[FV0AchId], nMode, qVecDet, sumAmplFV0A, ft0geom, fv0geom);
        }

        if (sumAmplFV0A > minAmplitude) {
          qVectFV0A[0] = qVecDet.Re();
          qVectFV0A[1] = qVecDet.Im();
        }
      }
    }

    int nTrkTPCPos = 0;
    int nTrkTPCNeg = 0;
    int nTrkTPCAll = 0;

    for (auto const& trk : tracks) {
      if (!selTrack(trk)) {
        continue;
      }
      histosQA.fill(HIST("ChTracks"), trk.pt(), trk.eta(), trk.phi(), cent);
      if (trk.eta() > cfgEtaMax) {
        continue;
      }
      if (trk.eta() < cfgEtaMin) {
        continue;
      }
      qVectTPCAll[0] += trk.pt() * std::cos(trk.phi() * nMode);
      qVectTPCAll[1] += trk.pt() * std::sin(trk.phi() * nMode);
      trkTPCAllLabel.push_back(trk.globalIndex());
      nTrkTPCAll++;
      if (std::abs(trk.eta()) < 0.1) {
        continue;
      }
      if (trk.eta() > 0 && (useDetector["QvectorTPCposs"] || useDetector["QvectorBPoss"])) {
        qVectTPCPos[0] += trk.pt() * std::cos(trk.phi() * nMode);
        qVectTPCPos[1] += trk.pt() * std::sin(trk.phi() * nMode);
        trkTPCPosLabel.push_back(trk.globalIndex());
        nTrkTPCPos++;
      } else if (trk.eta() < 0 && (useDetector["QvectorTPCnegs"] || useDetector["QvectorBNegs"])) {
        qVectTPCNeg[0] += trk.pt() * std::cos(trk.phi() * nMode);
        qVectTPCNeg[1] += trk.pt() * std::sin(trk.phi() * nMode);
        trkTPCNegLabel.push_back(trk.globalIndex());
        nTrkTPCNeg++;
      }
    }

    qVecRe.push_back(qVectFT0C[0]);
    qVecIm.push_back(qVectFT0C[1]);
    qVecRe.push_back(qVectFT0A[0]);
    qVecIm.push_back(qVectFT0A[1]);
    qVecRe.push_back(qVectFT0M[0]);
    qVecIm.push_back(qVectFT0M[1]);
    qVecRe.push_back(qVectFV0A[0]);
    qVecIm.push_back(qVectFV0A[1]);
    qVecRe.push_back(qVectTPCPos[0]);
    qVecIm.push_back(qVectTPCPos[1]);
    qVecRe.push_back(qVectTPCNeg[0]);
    qVecIm.push_back(qVectTPCNeg[1]);
    qVecRe.push_back(qVectTPCAll[0]);
    qVecIm.push_back(qVectTPCAll[1]);

    qVecAmp.push_back(sumAmplFT0C);
    qVecAmp.push_back(sumAmplFT0A);
    qVecAmp.push_back(sumAmplFT0M);
    qVecAmp.push_back(sumAmplFV0A);
    qVecAmp.push_back(static_cast<float>(nTrkTPCPos));
    qVecAmp.push_back(static_cast<float>(nTrkTPCNeg));
    qVecAmp.push_back(static_cast<float>(nTrkTPCAll));
  }

  void process(MyCollisions::iterator const& coll, aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&, MyTracks const& tracks)
  {
    std::vector<int> trkTPCPosLabel{};
    std::vector<int> trkTPCNegLabel{};
    std::vector<int> trkTPCAllLabel{};
    std::vector<float> qVecAmp{};

    std::vector<float> qVecReSp{};
    std::vector<float> qVecImSp{};
    std::vector<float> qVecReFT0CSp{};
    std::vector<float> qVecImFT0CSp{};
    std::vector<float> qVecReFT0ASp{};
    std::vector<float> qVecImFT0ASp{};
    std::vector<float> qVecReFT0MSp{};
    std::vector<float> qVecImFT0MSp{};
    std::vector<float> qVecReFV0ASp{};
    std::vector<float> qVecImFV0ASp{};
    std::vector<float> qVecReTPCPosSp{};
    std::vector<float> qVecImTPCPosSp{};
    std::vector<float> qVecReTPCNegSp{};
    std::vector<float> qVecImTPCNegSp{};
    std::vector<float> qVecReTPCAllSp{};
    std::vector<float> qVecImTPCAllSp{};

    std::vector<float> qVecReEse{};
    std::vector<float> qVecImEse{};
    std::vector<float> qVecReFT0CEse{};
    std::vector<float> qVecImFT0CEse{};
    std::vector<float> qVecReFT0AEse{};
    std::vector<float> qVecImFT0AEse{};
    std::vector<float> qVecReFT0MEse{};
    std::vector<float> qVecImFT0MEse{};
    std::vector<float> qVecReFV0AEse{};
    std::vector<float> qVecImFV0AEse{};
    std::vector<float> qVecReTPCPosEse{};
    std::vector<float> qVecImTPCPosEse{};
    std::vector<float> qVecReTPCNegEse{};
    std::vector<float> qVecImTPCNegEse{};
    std::vector<float> qVecReTPCAllEse{};
    std::vector<float> qVecImTPCAllEse{};

    auto bc = coll.bc_as<aod::BCsWithTimestamps>();
    int currentRun = bc.runNumber();
    if (runNumber != currentRun) {
      initCCDB(bc);
      runNumber = currentRun;
    }

    const float centAllEstim[4] = {
      coll.centFT0M(), coll.centFT0A(), coll.centFT0C(),
      coll.centFV0A()};
    cent = centAllEstim[cfgCentEsti];
    bool isCalibrated = true;
    if (cent < 0. || cent > cfgMaxCentrality) {
      cent = 110.;
      isCalibrated = false;
    }

    for (std::size_t id = 0; id < cfgnMods->size(); id++) {
      int nMode = cfgnMods->at(id);

      // Raw Q-vectors, no multiplicity normalization and no corrections
      std::vector<float> qVecReRaw{};
      std::vector<float> qVecImRaw{};
      calcQVec(nMode, coll, tracks, qVecReRaw, qVecImRaw, qVecAmp, trkTPCPosLabel, trkTPCNegLabel, trkTPCAllLabel);

      // Scalar Product Q-vectors, normalization by multiplicity/amplitude
      std::vector<float> nModeQVecReSp{};
      std::vector<float> nModeQVecImSp{};
      normalizeQVec(nModeQVecReSp, nModeQVecImSp, qVecReRaw, qVecImRaw, qVecAmp, MultNorms::kScalarProd);
      correctQVec(cent, nModeQVecReSp, nModeQVecImSp, corrsQvecSp[id], shiftProfileSp, nMode);
      // Add to summary vector
      qVecReSp.insert(qVecReSp.end(), nModeQVecReSp.begin(), nModeQVecReSp.end());
      qVecImSp.insert(qVecImSp.end(), nModeQVecImSp.begin(), nModeQVecImSp.end());

      // Pick the desired correction level for the Q-vectors to be stored in the analysis table
      int corrLevel = cfgCorrLevel == 0 ? 0 : cfgCorrLevel - 1;
      int nCorrections = static_cast<int>(kNCorrections);

      qVecReFT0CSp.push_back(nModeQVecReSp[kFT0C * nCorrections + corrLevel]);
      qVecImFT0CSp.push_back(nModeQVecImSp[kFT0C * nCorrections + corrLevel]);
      qVecReFT0ASp.push_back(nModeQVecReSp[kFT0A * nCorrections + corrLevel]);
      qVecImFT0ASp.push_back(nModeQVecImSp[kFT0A * nCorrections + corrLevel]);
      qVecReFT0MSp.push_back(nModeQVecReSp[kFT0M * nCorrections + corrLevel]);
      qVecImFT0MSp.push_back(nModeQVecImSp[kFT0M * nCorrections + corrLevel]);
      qVecReFV0ASp.push_back(nModeQVecReSp[kFV0A * nCorrections + corrLevel]);
      qVecImFV0ASp.push_back(nModeQVecImSp[kFV0A * nCorrections + corrLevel]);
      qVecReTPCPosSp.push_back(nModeQVecReSp[kTPCPos * nCorrections + corrLevel]);
      qVecImTPCPosSp.push_back(nModeQVecImSp[kTPCPos * nCorrections + corrLevel]);
      qVecReTPCNegSp.push_back(nModeQVecReSp[kTPCNeg * nCorrections + corrLevel]);
      qVecImTPCNegSp.push_back(nModeQVecImSp[kTPCNeg * nCorrections + corrLevel]);
      qVecReTPCAllSp.push_back(nModeQVecReSp[kTPCAll * nCorrections + corrLevel]);
      qVecImTPCAllSp.push_back(nModeQVecImSp[kTPCAll * nCorrections + corrLevel]);

      if (cfgProduceRedQVecs) {

        // Ese Q-vectors, normalization by sqrt(multiplicity/amplitude)
        std::vector<float> nModeQVecReEse{};
        std::vector<float> nModeQVecImEse{};
        normalizeQVec(nModeQVecReEse, nModeQVecImEse, qVecReRaw, qVecImRaw, qVecAmp, MultNorms::kEsE);
        correctQVec(cent, nModeQVecReEse, nModeQVecImEse, corrsQvecEse[id], shiftProfileEse, nMode);
        // Add to summary vector
        qVecReEse.insert(qVecReEse.end(), nModeQVecReEse.begin(), nModeQVecReEse.end());
        qVecImEse.insert(qVecImEse.end(), nModeQVecImEse.begin(), nModeQVecImEse.end());

        qVecReFT0CEse.push_back(nModeQVecReEse[kFT0C * nCorrections + corrLevel]);
        qVecImFT0CEse.push_back(nModeQVecImEse[kFT0C * nCorrections + corrLevel]);
        qVecReFT0AEse.push_back(nModeQVecReEse[kFT0A * nCorrections + corrLevel]);
        qVecImFT0AEse.push_back(nModeQVecImEse[kFT0A * nCorrections + corrLevel]);
        qVecReFT0MEse.push_back(nModeQVecReEse[kFT0M * nCorrections + corrLevel]);
        qVecImFT0MEse.push_back(nModeQVecImEse[kFT0M * nCorrections + corrLevel]);
        qVecReFV0AEse.push_back(nModeQVecReEse[kFV0A * nCorrections + corrLevel]);
        qVecImFV0AEse.push_back(nModeQVecImEse[kFV0A * nCorrections + corrLevel]);
        qVecReTPCPosEse.push_back(nModeQVecReEse[kTPCPos * nCorrections + corrLevel]);
        qVecImTPCPosEse.push_back(nModeQVecImEse[kTPCPos * nCorrections + corrLevel]);
        qVecReTPCNegEse.push_back(nModeQVecReEse[kTPCNeg * nCorrections + corrLevel]);
        qVecImTPCNegEse.push_back(nModeQVecImEse[kTPCNeg * nCorrections + corrLevel]);
        qVecReTPCAllEse.push_back(nModeQVecReEse[kTPCAll * nCorrections + corrLevel]);
        qVecImTPCAllEse.push_back(nModeQVecImEse[kTPCAll * nCorrections + corrLevel]);
      }
    }

    // Fill the columns of the Qvectors table.
    qVector(cent, isCalibrated, qVecReSp, qVecImSp, qVecAmp);
    if (useDetector["QvectorFT0Cs"])
      qVectorFT0C(isCalibrated, qVecReFT0CSp.at(0), qVecImFT0CSp.at(0), qVecAmp[kFT0C]);
    if (useDetector["QvectorFT0As"])
      qVectorFT0A(isCalibrated, qVecReFT0ASp.at(0), qVecImFT0ASp.at(0), qVecAmp[kFT0A]);
    if (useDetector["QvectorFT0Ms"])
      qVectorFT0M(isCalibrated, qVecReFT0MSp.at(0), qVecImFT0MSp.at(0), qVecAmp[kFT0M]);
    if (useDetector["QvectorFV0As"])
      qVectorFV0A(isCalibrated, qVecReFV0ASp.at(0), qVecImFV0ASp.at(0), qVecAmp[kFV0A]);
    if (useDetector["QvectorTPCposs"])
      qVectorTPCPos(isCalibrated, qVecReTPCPosSp.at(0), qVecImTPCPosSp.at(0), qVecAmp[kTPCPos], trkTPCPosLabel);
    if (useDetector["QvectorTPCnegs"])
      qVectorTPCNeg(isCalibrated, qVecReTPCNegSp.at(0), qVecImTPCNegSp.at(0), qVecAmp[kTPCNeg], trkTPCNegLabel);
    if (useDetector["QvectorTPCalls"])
      qVectorTPCAll(isCalibrated, qVecReTPCAllSp.at(0), qVecImTPCAllSp.at(0), qVecAmp[kTPCAll], trkTPCAllLabel);

    qVectorFT0CVec(isCalibrated, qVecReFT0CSp, qVecImFT0CSp, qVecAmp[kFT0C]);
    qVectorFT0AVec(isCalibrated, qVecReFT0ASp, qVecImFT0ASp, qVecAmp[kFT0A]);
    qVectorFT0MVec(isCalibrated, qVecReFT0MSp, qVecImFT0MSp, qVecAmp[kFT0M]);
    qVectorFV0AVec(isCalibrated, qVecReFV0ASp, qVecImFV0ASp, qVecAmp[kFV0A]);
    qVectorTPCPosVec(isCalibrated, qVecReTPCPosSp, qVecImTPCPosSp, qVecAmp[kTPCPos], trkTPCPosLabel);
    qVectorTPCNegVec(isCalibrated, qVecReTPCNegSp, qVecImTPCNegSp, qVecAmp[kTPCNeg], trkTPCNegLabel);
    qVectorTPCAllVec(isCalibrated, qVecReTPCAllSp, qVecImTPCAllSp, qVecAmp[kTPCAll], trkTPCAllLabel);

    // Deprecated, will be removed in future after transition time //
    if (useDetector["QvectorBPoss"])
      qVectorBPos(isCalibrated, qVecReTPCPosSp.at(0), qVecImTPCPosSp.at(0), qVecAmp[kTPCPos], trkTPCPosLabel);
    if (useDetector["QvectorBNegs"])
      qVectorBNeg(isCalibrated, qVecReTPCNegSp.at(0), qVecImTPCNegSp.at(0), qVecAmp[kTPCNeg], trkTPCNegLabel);
    if (useDetector["QvectorBTots"])
      qVectorBTot(isCalibrated, qVecReTPCAllSp.at(0), qVecImTPCAllSp.at(0), qVecAmp[kTPCAll], trkTPCAllLabel);

    qVectorBPosVec(isCalibrated, qVecReTPCPosSp, qVecImTPCPosSp, qVecAmp[kTPCPos], trkTPCPosLabel);
    qVectorBNegVec(isCalibrated, qVecReTPCNegSp, qVecImTPCNegSp, qVecAmp[kTPCNeg], trkTPCNegLabel);
    qVectorBTotVec(isCalibrated, qVecReTPCAllSp, qVecImTPCAllSp, qVecAmp[kTPCAll], trkTPCAllLabel);
    /////////////////////////////////////////////////////////////////

    if (cfgProduceRedQVecs) {
      eseQVector(cent, isCalibrated, qVecReEse, qVecImEse, qVecAmp);
      eseQVectorFT0C(isCalibrated, qVecReFT0CEse.at(0), qVecImFT0CEse.at(0), qVecAmp[kFT0C]);
      eseQVectorFT0A(isCalibrated, qVecReFT0AEse.at(0), qVecImFT0AEse.at(0), qVecAmp[kFT0A]);
      eseQVectorFT0M(isCalibrated, qVecReFT0MEse.at(0), qVecImFT0MEse.at(0), qVecAmp[kFT0M]);
      eseQVectorFV0A(isCalibrated, qVecReFV0AEse.at(0), qVecImFV0AEse.at(0), qVecAmp[kFV0A]);
      eseQVectorTPCPos(isCalibrated, qVecReTPCPosEse.at(0), qVecImTPCPosEse.at(0), qVecAmp[kTPCPos], trkTPCPosLabel);
      eseQVectorTPCNeg(isCalibrated, qVecReTPCNegEse.at(0), qVecImTPCNegEse.at(0), qVecAmp[kTPCNeg], trkTPCNegLabel);
      eseQVectorTPCAll(isCalibrated, qVecReTPCAllEse.at(0), qVecImTPCAllEse.at(0), qVecAmp[kTPCAll], trkTPCAllLabel);
      eseQVectorFT0CVec(isCalibrated, qVecReFT0CEse, qVecImFT0CEse, qVecAmp[kFT0C]);
      eseQVectorFT0AVec(isCalibrated, qVecReFT0AEse, qVecImFT0AEse, qVecAmp[kFT0A]);
      eseQVectorFT0MVec(isCalibrated, qVecReFT0MEse, qVecImFT0MEse, qVecAmp[kFT0M]);
      eseQVectorFV0AVec(isCalibrated, qVecReFV0AEse, qVecImFV0AEse, qVecAmp[kFV0A]);
      eseQVectorTPCPosVec(isCalibrated, qVecReTPCPosEse, qVecImTPCPosEse, qVecAmp[kTPCPos], trkTPCPosLabel);
      eseQVectorTPCNegVec(isCalibrated, qVecReTPCNegEse, qVecImTPCNegEse, qVecAmp[kTPCNeg], trkTPCNegLabel);
      eseQVectorTPCAllVec(isCalibrated, qVecReTPCAllEse, qVecImTPCAllEse, qVecAmp[kTPCAll], trkTPCAllLabel);
      eseQVectorPerc(qVecReFT0CEse.at(0), qVecImFT0CEse.at(0), qVecAmp[kFT0C],
                     qVecReFT0AEse.at(0), qVecImFT0AEse.at(0), qVecAmp[kFT0A],
                     qVecReFT0MEse.at(0), qVecImFT0MEse.at(0), qVecAmp[kFT0M],
                     qVecReFV0AEse.at(0), qVecImFV0AEse.at(0), qVecAmp[kFV0A],
                     qVecReTPCPosEse.at(0), qVecImTPCPosEse.at(0), qVecAmp[kTPCPos],
                     qVecReTPCNegEse.at(0), qVecImTPCNegEse.at(0), qVecAmp[kTPCNeg],
                     qVecReTPCAllEse.at(0), qVecImTPCAllEse.at(0), qVecAmp[kTPCAll]);
    }
  } // End process.
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsTable>(cfgc)};
}
