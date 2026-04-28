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
  enum {
    kFT0C = 0,
    kFT0A = 1,
    kFT0M,
    kFV0A,
    kTPCpos,
    kTPCneg,
    kTPCall,
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
  Configurable<std::string> cfgQvecEseCalibPath{"cfgQvecEseCalibPath", "Analysis/EventPlane/QVecEseCorrections", "CCDB path for EsE Q-vector calibration constants"};

  Configurable<bool> cfgShiftCorr{"cfgShiftCorr", false, "configurable flag for shift correction"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "", "CCDB path for shift correction"};

  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, ""};

  Configurable<bool> cfgUseFT0C{"cfgUseFT0C", false, "Initial value for using FT0C. By default obtained from DataModel."};
  Configurable<bool> cfgUseFT0A{"cfgUseFT0A", false, "Initial value for using FT0A. By default obtained from DataModel."};
  Configurable<bool> cfgUseFT0M{"cfgUseFT0M", false, "Initial value for using FT0M. By default obtained from DataModel."};
  Configurable<bool> cfgUseFV0A{"cfgUseFV0A", false, "Initial value for using FV0A. By default obtained from DataModel."};
  Configurable<bool> cfgUseTPCpos{"cfgUseTPCpos", false, "Initial value for using TPCpos. By default obtained from DataModel."};
  Configurable<bool> cfgUseTPCneg{"cfgUseTPCneg", false, "Initial value for using TPCneg. By default obtained from DataModel."};
  Configurable<bool> cfgUseTPCall{"cfgUseTPCall", false, "Initial value for using TPCall. By default obtained from DataModel."};
  Configurable<bool> cfgProduceRedQVecs{"cfgProduceRedQVecs", false, "Produce reduced Q-vectors for Event-Shape Engineering"};

  // Table.
  Produces<aod::Qvectors> qVector;
  Produces<aod::QvectorFT0Cs> qVectorFT0C;
  Produces<aod::QvectorFT0As> qVectorFT0A;
  Produces<aod::QvectorFT0Ms> qVectorFT0M;
  Produces<aod::QvectorFV0As> qVectorFV0A;
  Produces<aod::QvectorTPCposs> qVectorTPCpos;
  Produces<aod::QvectorTPCnegs> qVectorTPCneg;
  Produces<aod::QvectorTPCalls> qVectorTPCall;

  Produces<aod::QvectorFT0CVecs> qVectorFT0CVec;
  Produces<aod::QvectorFT0AVecs> qVectorFT0AVec;
  Produces<aod::QvectorFT0MVecs> qVectorFT0MVec;
  Produces<aod::QvectorFV0AVecs> qVectorFV0AVec;
  Produces<aod::QvectorTPCposVecs> qVectorTPCposVec;
  Produces<aod::QvectorTPCnegVecs> qVectorTPCnegVec;
  Produces<aod::QvectorTPCallVecs> qVectorTPCallVec;

  Produces<aod::EseQvectors> eseQVector;
  Produces<aod::EseQvecFT0Cs> eseQVectorFT0C;
  Produces<aod::EseQvecFT0As> eseQVectorFT0A;
  Produces<aod::EseQvecFT0Ms> eseQVectorFT0M;
  Produces<aod::EseQvecFV0As> eseQVectorFV0A;
  Produces<aod::EseQvecTPCposs> eseQVectorTPCpos;
  Produces<aod::EseQvecTPCnegs> eseQVectorTPCneg;
  Produces<aod::EseQvecTPCalls> eseQVectorTPCall;

  Produces<aod::EseQvecFT0CVecs> eseQVectorFT0CVec;
  Produces<aod::EseQvecFT0AVecs> eseQVectorFT0AVec;
  Produces<aod::EseQvecFT0MVecs> eseQVectorFT0MVec;
  Produces<aod::EseQvecFV0AVecs> eseQVectorFV0AVec;
  Produces<aod::EseQvecTPCposVecs> eseQVectorTPCposVec;
  Produces<aod::EseQvecTPCnegVecs> eseQVectorTPCnegVec;
  Produces<aod::EseQvecTPCallVecs> eseQVectorTPCallVec;
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

  int runNumber{-1};
  float cent;

  std::vector<TH3F*> corrsQvecSp{};
  std::vector<TH3F*> corrsQvecEse{};
  std::vector<TProfile3D*> shiftprofile{};

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
          std::string table_name_with_vector = det.first; // for replacing s with Vecs at the end.
          if (input.matcher.binding == det.first || input.matcher.binding == table_name_with_vector.replace(table_name_with_vector.size() - 1, 1, "Vecs")) {
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
      fullPath = cfgQvecEseCalibPath;
      fullPath += "/eseq";
      fullPath += std::to_string(ind);
      auto modeCorrQvecEse = getForTsOrRun<TH3F>(fullPath, timestamp, runnumber);
      if (!modeCorrQvecEse) {
        fullPath = cfgQvecCalibPath; // cfgQvecEseCalibPath;
        fullPath += "/v2"; // "/eseq2";
        modeCorrQvecEse = getForTsOrRun<TH3F>(fullPath, timestamp, runnumber);
      }
      corrsQvecEse.push_back(modeCorrQvecEse);
    }

    if (cfgShiftCorr) {
      shiftprofile.clear();
      for (std::size_t i = 0; i < cfgnMods->size(); i++) {
        int ind = cfgnMods->at(i);
        fullPath = cfgShiftPath;
        fullPath += "/v";
        fullPath += std::to_string(ind);
        auto objshift = getForTsOrRun<TProfile3D>(fullPath, timestamp, runnumber);
        shiftprofile.push_back(objshift);
      }
    }

    fullPath = cfgGainEqPath;
    fullPath += "/FT0";
    const auto objft0Gain = getForTsOrRun<std::vector<float>>(fullPath, timestamp, runnumber);
    if (!objft0Gain || cfgCorrLevel == 0) {
      for (auto i{0u}; i < 208; i++) {
        FT0RelGainConst.push_back(1.);
      }
    } else {
      FT0RelGainConst = *(objft0Gain);
    }

    fullPath = cfgGainEqPath;
    fullPath += "/FV0";
    const auto objfv0Gain = getForTsOrRun<std::vector<float>>(fullPath, timestamp, runnumber);
    if (!objfv0Gain || cfgCorrLevel == 0) {
      for (auto i{0u}; i < 48; i++) {
        FV0RelGainConst.push_back(1.);
      }
    } else {
      FV0RelGainConst = *(objfv0Gain);
    }
  }

  template <typename TrackType>
  bool SelTrack(const TrackType track)
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
  /// \param QvecReNorm is the vector with the normalized real part of the q-vector for each detector and correction step
  /// \param QvecImNorm is the vector with the normalized imaginary part of the q-vector for each detector and correction step
  /// \param QvecReRaw is the vector with the raw real part of the q-vector for each detector and correction step
  /// \param QvecImRaw is the vector with the raw imaginary part of the q-vector for each detector and correction step
  /// \param QvecAmp is the vector with the amplitude of the q-vector for each detector and correction step
  /// \param normType is the type of normalization to apply to the q-vectors
  void NormalizeQvec(std::vector<float>& QvecReNorm, 
                     std::vector<float>& QvecImNorm, 
                     std::vector<float> QvecReRaw, 
                     std::vector<float> QvecImRaw, 
                     std::vector<float>& QvecAmp, 
                     MultNorms normType) 
  {
    for (std::size_t i = 0; i < kNDetectors; i++) {
      float qVecDetReNorm{999.}, qVecDetImNorm{999.};
      if (QvecAmp[i] > 1e-8) {
        switch (normType) {
          case MultNorms::kScalarProd:
            qVecDetReNorm = QvecReRaw[i] / QvecAmp[i];
            qVecDetImNorm = QvecImRaw[i] / QvecAmp[i];
            break;
          case MultNorms::kEsE:
            qVecDetReNorm = QvecReRaw[i] / std::sqrt(QvecAmp[i]);
            qVecDetImNorm = QvecImRaw[i] / std::sqrt(QvecAmp[i]);
            break;
          default:
            LOGP(fatal, "Undefined normalization type for Q-vector amplitude. Check the configuration.");
            break;
        }
        std::cout << "[NORMALIZED] " << i << " Re: " << qVecDetReNorm << ", Im: " << qVecDetImNorm << ", amp: " << QvecAmp[i] << std::endl;
      }
      for (int iCorr=0; iCorr < Corrections::kNCorrections; iCorr++) {
        QvecReNorm.push_back(qVecDetReNorm);
        QvecImNorm.push_back(qVecDetImNorm);
      }
    }
  }

  /// Function to calculate the un-normalized q-vectors
  /// \param cent is the collision centrality
  /// \param qvecRe is the vector with the real part of the q-vector for each detector and correction step
  /// \param qvecIm is the vector with the imaginary part of the q-vector for each detector and correction step
  /// \param histsCorrs is the vector with the histograms with the correction constants for each detector and correction step
  /// \param nMode is the modulation of interest
  void CorrectQvec(float cent, std::vector<float>& qvecRe, std::vector<float>& qvecIm, TH3F* histsCorrs, int nMode) {
    int nCorrections = static_cast<int>(Corrections::kNCorrections);
    if (cent < cfgMaxCentrality) {
      for (auto i{0u}; i < kTPCall + 1; i++) {
        int idxDet = i * Corrections::kNCorrections;
        helperEP.DoRecenter(qvecRe[idxDet + Corrections::kRecenter], qvecIm[idxDet + Corrections::kRecenter],
                            histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));

        helperEP.DoRecenter(qvecRe[idxDet + Corrections::kTwist], qvecIm[idxDet + Corrections::kTwist],
                            histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));
        helperEP.DoTwist(qvecRe[idxDet + Corrections::kTwist], qvecIm[idxDet + Corrections::kTwist],
                         histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 3, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 4, i + 1));

        helperEP.DoRecenter(qvecRe[idxDet + Corrections::kRescale], qvecIm[idxDet + Corrections::kRescale],
                            histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));
        helperEP.DoTwist(qvecRe[idxDet + Corrections::kRescale], qvecIm[idxDet + Corrections::kRescale],
                         histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 3, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 4, i + 1));
        helperEP.DoRescale(qvecRe[idxDet + Corrections::kRescale], qvecIm[idxDet + Corrections::kRescale],
                           histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 5, i + 1), histsCorrs->GetBinContent(static_cast<int>(cent) + 1, 6, i + 1));
      }
      if (cfgShiftCorr) {
        auto deltapsiFT0C = 0.0;
        auto deltapsiFT0A = 0.0;
        auto deltapsiFT0M = 0.0;
        auto deltapsiFV0A = 0.0;
        auto deltapsiTPCpos = 0.0;
        auto deltapsiTPCneg = 0.0;
        auto deltapsiTPCall = 0.0;

        auto psidefFT0C = TMath::ATan2(qvecIm[kFT0C * nCorrections + Corrections::kRescale], qvecRe[kFT0C * nCorrections + Corrections::kRescale]) / static_cast<float>(nMode);
        auto psidefFT0A = TMath::ATan2(qvecIm[kFT0A * nCorrections + Corrections::kRescale], qvecRe[kFT0A * nCorrections + Corrections::kRescale]) / static_cast<float>(nMode);
        auto psidefFT0M = TMath::ATan2(qvecIm[kFT0M * nCorrections + Corrections::kRescale], qvecRe[kFT0M * nCorrections + Corrections::kRescale]) / static_cast<float>(nMode);
        auto psidefFV0A = TMath::ATan2(qvecIm[kFV0A * nCorrections + Corrections::kRescale], qvecRe[kFV0A * nCorrections + Corrections::kRescale]) / static_cast<float>(nMode);
        auto psidefTPCpos = TMath::ATan2(qvecIm[kTPCpos * nCorrections + Corrections::kRescale], qvecRe[kTPCpos * nCorrections + Corrections::kRescale]) / static_cast<float>(nMode);
        auto psidefTPCneg = TMath::ATan2(qvecIm[kTPCneg * nCorrections + Corrections::kRescale], qvecRe[kTPCneg * nCorrections + Corrections::kRescale]) / static_cast<float>(nMode);
        auto psidefTPCall = TMath::ATan2(qvecIm[kTPCall * nCorrections + Corrections::kRescale], qvecRe[kTPCall * nCorrections + Corrections::kRescale]) / static_cast<float>(nMode);

        for (int ishift = 1; ishift <= 10; ishift++) {
          auto coeffshiftxFT0C = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kFT0C, ishift - 0.5));
          auto coeffshiftyFT0C = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kFT0C + 1, ishift - 0.5));
          auto coeffshiftxFT0A = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kFT0A, ishift - 0.5));
          auto coeffshiftyFT0A = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kFT0A + 1, ishift - 0.5));
          auto coeffshiftxFT0M = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kFT0M, ishift - 0.5));
          auto coeffshiftyFT0M = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kFT0M + 1, ishift - 0.5));
          auto coeffshiftxFV0A = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kFV0A, ishift - 0.5));
          auto coeffshiftyFV0A = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kFV0A + 1, ishift - 0.5));
          auto coeffshiftxTPCpos = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kTPCpos, ishift - 0.5));
          auto coeffshiftyTPCpos = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kTPCpos + 1, ishift - 0.5));
          auto coeffshiftxTPCneg = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kTPCneg, ishift - 0.5));
          auto coeffshiftyTPCneg = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kTPCneg + 1, ishift - 0.5));
          auto coeffshiftxTPCall = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kTPCall, ishift - 0.5));
          auto coeffshiftyTPCall = shiftprofile.at(nMode - 2)->GetBinContent(shiftprofile.at(nMode - 2)->FindBin(cent, 2 * kTPCall + 1, ishift - 0.5));

          deltapsiFT0C += ((2. / (1.0 * ishift)) * (-coeffshiftxFT0C * TMath::Cos(ishift * static_cast<float>(nMode) * psidefFT0C) + coeffshiftyFT0C * TMath::Sin(ishift * static_cast<float>(nMode) * psidefFT0C))) / static_cast<float>(nMode);
          deltapsiFT0A += ((2. / (1.0 * ishift)) * (-coeffshiftxFT0A * TMath::Cos(ishift * static_cast<float>(nMode) * psidefFT0A) + coeffshiftyFT0A * TMath::Sin(ishift * static_cast<float>(nMode) * psidefFT0A))) / static_cast<float>(nMode);
          deltapsiFT0M += ((2. / (1.0 * ishift)) * (-coeffshiftxFT0M * TMath::Cos(ishift * static_cast<float>(nMode) * psidefFT0M) + coeffshiftyFT0M * TMath::Sin(ishift * static_cast<float>(nMode) * psidefFT0M))) / static_cast<float>(nMode);
          deltapsiFV0A += ((2. / (1.0 * ishift)) * (-coeffshiftxFV0A * TMath::Cos(ishift * static_cast<float>(nMode) * psidefFV0A) + coeffshiftyFV0A * TMath::Sin(ishift * static_cast<float>(nMode) * psidefFV0A))) / static_cast<float>(nMode);
          deltapsiTPCpos += ((2. / (1.0 * ishift)) * (-coeffshiftxTPCpos * TMath::Cos(ishift * static_cast<float>(nMode) * psidefTPCpos) + coeffshiftyTPCpos * TMath::Sin(ishift * static_cast<float>(nMode) * psidefTPCpos))) / static_cast<float>(nMode);
          deltapsiTPCneg += ((2. / (1.0 * ishift)) * (-coeffshiftxTPCneg * TMath::Cos(ishift * static_cast<float>(nMode) * psidefTPCneg) + coeffshiftyTPCneg * TMath::Sin(ishift * static_cast<float>(nMode) * psidefTPCneg))) / static_cast<float>(nMode);
          deltapsiTPCall += ((2. / (1.0 * ishift)) * (-coeffshiftxTPCall * TMath::Cos(ishift * static_cast<float>(nMode) * psidefTPCall) + coeffshiftyTPCall * TMath::Sin(ishift * static_cast<float>(nMode) * psidefTPCall))) / static_cast<float>(nMode);
        }

        deltapsiFT0C *= static_cast<float>(nMode);
        deltapsiFT0A *= static_cast<float>(nMode);
        deltapsiFT0M *= static_cast<float>(nMode);
        deltapsiFV0A *= static_cast<float>(nMode);
        deltapsiTPCpos *= static_cast<float>(nMode);
        deltapsiTPCneg *= static_cast<float>(nMode);
        deltapsiTPCall *= static_cast<float>(nMode);

        float qvecReShiftedFT0C = qvecRe[kFT0C * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiFT0C) - qvecIm[kFT0C * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiFT0C);
        float qvecImShiftedFT0C = qvecRe[kFT0C * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiFT0C) + qvecIm[kFT0C * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiFT0C);
        float qvecReShiftedFT0A = qvecRe[kFT0A * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiFT0A) - qvecIm[kFT0A * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiFT0A);
        float qvecImShiftedFT0A = qvecRe[kFT0A * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiFT0A) + qvecIm[kFT0A * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiFT0A);
        float qvecReShiftedFT0M = qvecRe[kFT0M * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiFT0M) - qvecIm[kFT0M * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiFT0M);
        float qvecImShiftedFT0M = qvecRe[kFT0M * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiFT0M) + qvecIm[kFT0M * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiFT0M);
        float qvecReShiftedFV0A = qvecRe[kFV0A * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiFV0A) - qvecIm[kFV0A * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiFV0A);
        float qvecImShiftedFV0A = qvecRe[kFV0A * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiFV0A) + qvecIm[kFV0A * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiFV0A);
        float qvecReShiftedTPCpos = qvecRe[kTPCpos * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiTPCpos) - qvecIm[kTPCpos * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiTPCpos);
        float qvecImShiftedTPCpos = qvecRe[kTPCpos * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiTPCpos) + qvecIm[kTPCpos * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiTPCpos);
        float qvecReShiftedTPCneg = qvecRe[kTPCneg * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiTPCneg) - qvecIm[kTPCneg * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiTPCneg);
        float qvecImShiftedTPCneg = qvecRe[kTPCneg * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiTPCneg) + qvecIm[kTPCneg * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiTPCneg);
        float qvecReShiftedTPCall = qvecRe[kTPCall * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiTPCall) - qvecIm[kTPCall * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiTPCall);
        float qvecImShiftedTPCall = qvecRe[kTPCall * nCorrections + Corrections::kRescale] * TMath::Sin(deltapsiTPCall) + qvecIm[kTPCall * nCorrections + Corrections::kRescale] * TMath::Cos(deltapsiTPCall);

        qvecRe[kFT0C * nCorrections + Corrections::kRescale] = qvecReShiftedFT0C;
        qvecIm[kFT0C * nCorrections + Corrections::kRescale] = qvecImShiftedFT0C;
        qvecRe[kFT0A * nCorrections + Corrections::kRescale] = qvecReShiftedFT0A;
        qvecIm[kFT0A * nCorrections + Corrections::kRescale] = qvecImShiftedFT0A;
        qvecRe[kFT0M * nCorrections + Corrections::kRescale] = qvecReShiftedFT0M;
        qvecIm[kFT0M * nCorrections + Corrections::kRescale] = qvecImShiftedFT0M;
        qvecRe[kFV0A * nCorrections + Corrections::kRescale] = qvecReShiftedFV0A;
        qvecIm[kFV0A * nCorrections + Corrections::kRescale] = qvecImShiftedFV0A;
        qvecRe[kTPCpos * nCorrections + Corrections::kRescale] = qvecReShiftedTPCpos;
        qvecIm[kTPCpos * nCorrections + Corrections::kRescale] = qvecImShiftedTPCpos;
        qvecRe[kTPCneg * nCorrections + Corrections::kRescale] = qvecReShiftedTPCneg;
        qvecIm[kTPCneg * nCorrections + Corrections::kRescale] = qvecImShiftedTPCneg;
        qvecRe[kTPCall * nCorrections + Corrections::kRescale] = qvecReShiftedTPCall;
        qvecIm[kTPCall * nCorrections + Corrections::kRescale] = qvecImShiftedTPCall;
      }
    }
  }

  /// Function to calculate the un-normalized q-vectors
  /// \param nMode is the harmonic number of the q-vector
  /// \param coll is the collision object
  /// \param track are the tracks associated to the collision
  /// \param QvecRe is the vector with the real part of the q-vector for each detector
  /// \param QvecIm is the vector with the imaginary part of the q-vector for each detector
  /// \param QvecAmp is the vector with the amplitude of the signal in each detector
  /// \param TrkTPCposLabel is the vector with the number of TPC tracks with positive eta
  /// \param TrkTPCnegLabel is the vector with the number of TPC tracks with negative eta
  /// \param TrkTPCallLabel is the vector with the number of TPC tracks with any eta
  template <typename Nmode, typename CollType, typename TrackType>
  void CalQvec(const Nmode nMode, const CollType& coll, const TrackType& track, std::vector<float>& QvecRe, std::vector<float>& QvecIm, std::vector<float>& QvecAmp, std::vector<int>& TrkTPCposLabel, std::vector<int>& TrkTPCnegLabel, std::vector<int>& TrkTPCallLabel)
  {
    float qVectFT0A[2] = {-999., -999.};
    float qVectFT0C[2] = {-999., -999.};
    float qVectFT0M[2] = {-999., -999.};
    float qVectFV0A[2] = {-999., -999.};
    float qVectTPCpos[2] = {0., 0.};  // Always computed
    float qVectTPCneg[2] = {0., 0.};  // Always computed
    float qVectTPCall[2] = {0., 0.};  // Always computed

    TComplex QvecDet(0);
    TComplex QvecFT0M(0);
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

          helperEP.SumQvectors(0, FT0AchId, ampl / FT0RelGainConst[FT0AchId], nMode, QvecDet, sumAmplFT0A, ft0geom, fv0geom);
          helperEP.SumQvectors(0, FT0AchId, ampl / FT0RelGainConst[FT0AchId], nMode, QvecFT0M, sumAmplFT0M, ft0geom, fv0geom);
        }
        if (sumAmplFT0A > 1e-8) {
          qVectFT0A[0] = QvecDet.Re();
          qVectFT0A[1] = QvecDet.Im();
        }
      }

      if (useDetector["QvectorFT0Cs"]) {
        QvecDet = TComplex(0., 0.);
        for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
          float ampl = ft0.amplitudeC()[iChC];
          int FT0CchId = ft0.channelC()[iChC] + 96;

          histosQA.fill(HIST("FT0Amp"), ampl, FT0CchId);
          histosQA.fill(HIST("FT0AmpCor"), ampl / FT0RelGainConst[FT0CchId], FT0CchId);

          helperEP.SumQvectors(0, FT0CchId, ampl / FT0RelGainConst[FT0CchId], nMode, QvecDet, sumAmplFT0C, ft0geom, fv0geom);
          helperEP.SumQvectors(0, FT0CchId, ampl / FT0RelGainConst[FT0CchId], nMode, QvecFT0M, sumAmplFT0M, ft0geom, fv0geom);
        }

        if (sumAmplFT0C > 1e-8) {
          qVectFT0C[0] = QvecDet.Re();
          qVectFT0C[1] = QvecDet.Im();
        }
        if (sumAmplFT0M > 1e-8 && useDetector["QvectorFT0Ms"]) {
          qVectFT0M[0] = QvecFT0M.Re();
          qVectFT0M[1] = QvecFT0M.Im();
        }
      }

      QvecDet = TComplex(0., 0.);
      sumAmplFV0A = 0;
      if (coll.has_foundFV0() && useDetector["QvectorFV0As"]) {
        auto fv0 = coll.foundFV0();

        for (std::size_t iCh = 0; iCh < fv0.channel().size(); iCh++) {
          float ampl = fv0.amplitude()[iCh];
          int FV0AchId = fv0.channel()[iCh];
          histosQA.fill(HIST("FV0Amp"), ampl, FV0AchId);
          histosQA.fill(HIST("FV0AmpCor"), ampl / FV0RelGainConst[FV0AchId], FV0AchId);

          helperEP.SumQvectors(1, FV0AchId, ampl / FV0RelGainConst[FV0AchId], nMode, QvecDet, sumAmplFV0A, ft0geom, fv0geom);
        }

        if (sumAmplFV0A > 1e-8) {
          qVectFV0A[0] = QvecDet.Re();
          qVectFV0A[1] = QvecDet.Im();
        }
      }
    }

    int nTrkTPCpos = 0;
    int nTrkTPCneg = 0;
    int nTrkTPCall = 0;

    for (auto const& trk : track) {
      if (!SelTrack(trk)) {
        continue;
      }
      histosQA.fill(HIST("ChTracks"), trk.pt(), trk.eta(), trk.phi(), cent);
      if (trk.eta() > cfgEtaMax) {
        continue;
      }
      if (trk.eta() < cfgEtaMin) {
        continue;
      }
      qVectTPCall[0] += trk.pt() * std::cos(trk.phi() * nMode);
      qVectTPCall[1] += trk.pt() * std::sin(trk.phi() * nMode);
      TrkTPCallLabel.push_back(trk.globalIndex());
      nTrkTPCall++;
      if (std::abs(trk.eta()) < 0.1) {
        continue;
      }
      if (trk.eta() > 0 && (useDetector["QvectorTPCposs"] || useDetector["QvectorBPoss"])) {
        qVectTPCpos[0] += trk.pt() * std::cos(trk.phi() * nMode);
        qVectTPCpos[1] += trk.pt() * std::sin(trk.phi() * nMode);
        TrkTPCposLabel.push_back(trk.globalIndex());
        nTrkTPCpos++;
      } else if (trk.eta() < 0 && (useDetector["QvectorTPCnegs"] || useDetector["QvectorBNegs"])) {
        qVectTPCneg[0] += trk.pt() * std::cos(trk.phi() * nMode);
        qVectTPCneg[1] += trk.pt() * std::sin(trk.phi() * nMode);
        TrkTPCnegLabel.push_back(trk.globalIndex());
        nTrkTPCneg++;
      }
    }

    QvecRe.push_back(qVectFT0C[0]);
    QvecIm.push_back(qVectFT0C[1]);
    QvecRe.push_back(qVectFT0A[0]);
    QvecIm.push_back(qVectFT0A[1]);
    QvecRe.push_back(qVectFT0M[0]);
    QvecIm.push_back(qVectFT0M[1]);
    QvecRe.push_back(qVectFV0A[0]);
    QvecIm.push_back(qVectFV0A[1]);
    QvecRe.push_back(qVectTPCpos[0]);
    QvecIm.push_back(qVectTPCpos[1]);
    QvecRe.push_back(qVectTPCneg[0]);
    QvecIm.push_back(qVectTPCneg[1]);
    QvecRe.push_back(qVectTPCall[0]);
    QvecIm.push_back(qVectTPCall[1]);

    QvecAmp.push_back(sumAmplFT0C);
    QvecAmp.push_back(sumAmplFT0A);
    QvecAmp.push_back(sumAmplFT0M);
    QvecAmp.push_back(sumAmplFV0A);
    QvecAmp.push_back(static_cast<float>(nTrkTPCpos));
    QvecAmp.push_back(static_cast<float>(nTrkTPCneg));
    QvecAmp.push_back(static_cast<float>(nTrkTPCall));

    LOG(info) << "[RAW] qVectFT0A: " << qVectFT0A[0] << ", " << qVectFT0A[1] << ", ampl: " << sumAmplFT0A;
    LOG(info) << "[RAW] qVectFT0C: " << qVectFT0C[0] << ", " << qVectFT0C[1] << ", ampl: " << sumAmplFT0C;
    LOG(info) << "[RAW] qVectFT0M: " << qVectFT0M[0] << ", " << qVectFT0M[1] << ", ampl: " << sumAmplFT0M;
    LOG(info) << "[RAW] qVectFV0A: " << qVectFV0A[0] << ", " << qVectFV0A[1] << ", ampl: " << sumAmplFV0A;
    LOG(info) << "[RAW] qVectTPCpos: " << qVectTPCpos[0] << ", " << qVectTPCpos[1] << ", nTrk: " << nTrkTPCpos;
    LOG(info) << "[RAW] qVectTPCneg: " << qVectTPCneg[0] << ", " << qVectTPCneg[1] << ", nTrk: " << nTrkTPCneg;
    LOG(info) << "[RAW] qVectTPCall: " << qVectTPCall[0] << ", " << qVectTPCall[1] << ", nTrk: " << nTrkTPCall;

  }

  void process(MyCollisions::iterator const& coll, aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&, MyTracks const& tracks)
  {
    LOG(info) << "---------------------------- Processing Event ---------------------------";
    std::vector<int> TrkTPCposLabel{};
    std::vector<int> TrkTPCnegLabel{};
    std::vector<int> TrkTPCallLabel{};
    std::vector<float> qvecAmp{};

    std::vector<float> qvecReSp{};
    std::vector<float> qvecImSp{};
    std::vector<float> qvecReFT0CSp{};
    std::vector<float> qvecImFT0CSp{};
    std::vector<float> qvecReFT0ASp{};
    std::vector<float> qvecImFT0ASp{};
    std::vector<float> qvecReFT0MSp{};
    std::vector<float> qvecImFT0MSp{};
    std::vector<float> qvecReFV0ASp{};
    std::vector<float> qvecImFV0ASp{};
    std::vector<float> qvecReTPCposSp{};
    std::vector<float> qvecImTPCposSp{};
    std::vector<float> qvecReTPCnegSp{};
    std::vector<float> qvecImTPCnegSp{};
    std::vector<float> qvecReTPCallSp{};
    std::vector<float> qvecImTPCallSp{};

    std::vector<float> qvecReEse{};
    std::vector<float> qvecImEse{};
    std::vector<float> qvecReFT0CEse{};
    std::vector<float> qvecImFT0CEse{};
    std::vector<float> qvecReFT0AEse{};
    std::vector<float> qvecImFT0AEse{};
    std::vector<float> qvecReFT0MEse{};
    std::vector<float> qvecImFT0MEse{};
    std::vector<float> qvecReFV0AEse{};
    std::vector<float> qvecImFV0AEse{};
    std::vector<float> qvecReTPCposEse{};
    std::vector<float> qvecImTPCposEse{};
    std::vector<float> qvecReTPCnegEse{};
    std::vector<float> qvecImTPCnegEse{};
    std::vector<float> qvecReTPCallEse{};
    std::vector<float> qvecImTPCallEse{};

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
    bool IsCalibrated = true;
    if (cent < 0. || cent > cfgMaxCentrality) {
      cent = 110.;
      IsCalibrated = false;
    }

    for (std::size_t id = 0; id < cfgnMods->size(); id++) {
      int nMode = cfgnMods->at(id);

      // Raw Q-vectors, no multiplicity normalization and no corrections
      std::vector<float> qvecReRaw{};
      std::vector<float> qvecImRaw{};
      CalQvec(nMode, coll, tracks, qvecReRaw, qvecImRaw, qvecAmp, TrkTPCposLabel, TrkTPCnegLabel, TrkTPCallLabel);

      // Scalar Product Q-vectors, normalization by multiplicity/amplitude
      std::vector<float> nModeQvecReSp{};
      std::vector<float> nModeQvecImSp{};
      NormalizeQvec(nModeQvecReSp, nModeQvecImSp, qvecReRaw, qvecImRaw, qvecAmp, MultNorms::kScalarProd);
      CorrectQvec(cent, nModeQvecReSp, nModeQvecImSp, corrsQvecSp[id], nMode);
      // Add to summary vector
      qvecReSp.insert(qvecReSp.end(), nModeQvecReSp.begin(), nModeQvecReSp.end());
      qvecImSp.insert(qvecImSp.end(), nModeQvecImSp.begin(), nModeQvecImSp.end());

      // Ese Q-vectors, normalization by sqrt(multiplicity/amplitude)
      std::vector<float> nModeQvecReEse{};
      std::vector<float> nModeQvecImEse{};
      NormalizeQvec(nModeQvecReEse, nModeQvecImEse, qvecReRaw, qvecImRaw, qvecAmp, MultNorms::kEsE);
      CorrectQvec(cent, nModeQvecReEse, nModeQvecImEse, corrsQvecEse[id], nMode);
      // Add to summary vector
      qvecReEse.insert(qvecReEse.end(), nModeQvecReEse.begin(), nModeQvecReEse.end());
      qvecImEse.insert(qvecImEse.end(), nModeQvecImEse.begin(), nModeQvecImEse.end());

      // Pick the desired correction level for the Q-vectors to be stored in the analysis table
      int CorrLevel = cfgCorrLevel == 0 ? 0 : cfgCorrLevel - 1;
      int nCorrections = static_cast<int>(Corrections::kNCorrections);

      qvecReFT0CSp.push_back(nModeQvecReSp[kFT0C * nCorrections + CorrLevel]);
      qvecImFT0CSp.push_back(nModeQvecImSp[kFT0C * nCorrections + CorrLevel]);
      qvecReFT0ASp.push_back(nModeQvecReSp[kFT0A * nCorrections + CorrLevel]);
      qvecImFT0ASp.push_back(nModeQvecImSp[kFT0A * nCorrections + CorrLevel]);
      qvecReFT0MSp.push_back(nModeQvecReSp[kFT0M * nCorrections + CorrLevel]);
      qvecImFT0MSp.push_back(nModeQvecImSp[kFT0M * nCorrections + CorrLevel]);
      qvecReFV0ASp.push_back(nModeQvecReSp[kFV0A * nCorrections + CorrLevel]);
      qvecImFV0ASp.push_back(nModeQvecImSp[kFV0A * nCorrections + CorrLevel]);
      qvecReTPCposSp.push_back(nModeQvecReSp[kTPCpos * nCorrections + CorrLevel]);
      qvecImTPCposSp.push_back(nModeQvecImSp[kTPCpos * nCorrections + CorrLevel]);
      qvecReTPCnegSp.push_back(nModeQvecReSp[kTPCneg * nCorrections + CorrLevel]);
      qvecImTPCnegSp.push_back(nModeQvecImSp[kTPCneg * nCorrections + CorrLevel]);
      qvecReTPCallSp.push_back(nModeQvecReSp[kTPCall * nCorrections + CorrLevel]);
      qvecImTPCallSp.push_back(nModeQvecImSp[kTPCall * nCorrections + CorrLevel]);

      qvecReFT0CEse.push_back(nModeQvecReEse[kFT0C * nCorrections + CorrLevel]);
      qvecImFT0CEse.push_back(nModeQvecImEse[kFT0C * nCorrections + CorrLevel]);
      qvecReFT0AEse.push_back(nModeQvecReEse[kFT0A * nCorrections + CorrLevel]);
      qvecImFT0AEse.push_back(nModeQvecImEse[kFT0A * nCorrections + CorrLevel]);
      qvecReFT0MEse.push_back(nModeQvecReEse[kFT0M * nCorrections + CorrLevel]);
      qvecImFT0MEse.push_back(nModeQvecImEse[kFT0M * nCorrections + CorrLevel]);
      qvecReFV0AEse.push_back(nModeQvecReEse[kFV0A * nCorrections + CorrLevel]);
      qvecImFV0AEse.push_back(nModeQvecImEse[kFV0A * nCorrections + CorrLevel]);
      qvecReTPCposEse.push_back(nModeQvecReEse[kTPCpos * nCorrections + CorrLevel]);
      qvecImTPCposEse.push_back(nModeQvecImEse[kTPCpos * nCorrections + CorrLevel]);
      qvecReTPCnegEse.push_back(nModeQvecReEse[kTPCneg * nCorrections + CorrLevel]);
      qvecImTPCnegEse.push_back(nModeQvecImEse[kTPCneg * nCorrections + CorrLevel]);
      qvecReTPCallEse.push_back(nModeQvecReEse[kTPCall * nCorrections + CorrLevel]);
      qvecImTPCallEse.push_back(nModeQvecImEse[kTPCall * nCorrections + CorrLevel]);
    }

    // Fill the columns of the Qvectors table.
    qVector(cent, IsCalibrated, qvecReSp, qvecImSp, qvecAmp);
    if (useDetector["QvectorFT0Cs"])
      qVectorFT0C(IsCalibrated, qvecReFT0CSp.at(0), qvecImFT0CSp.at(0), qvecAmp[kFT0C]);
    if (useDetector["QvectorFT0As"])
      qVectorFT0A(IsCalibrated, qvecReFT0ASp.at(0), qvecImFT0ASp.at(0), qvecAmp[kFT0A]);
    if (useDetector["QvectorFT0Ms"])
      qVectorFT0M(IsCalibrated, qvecReFT0MSp.at(0), qvecImFT0MSp.at(0), qvecAmp[kFT0M]);
    if (useDetector["QvectorFV0As"])
      qVectorFV0A(IsCalibrated, qvecReFV0ASp.at(0), qvecImFV0ASp.at(0), qvecAmp[kFV0A]);
    if (useDetector["QvectorTPCposs"])
      qVectorTPCpos(IsCalibrated, qvecReTPCposSp.at(0), qvecImTPCposSp.at(0), qvecAmp[kTPCpos], TrkTPCposLabel);
    if (useDetector["QvectorTPCnegs"])
      qVectorTPCneg(IsCalibrated, qvecReTPCnegSp.at(0), qvecImTPCnegSp.at(0), qvecAmp[kTPCneg], TrkTPCnegLabel);
    if (useDetector["QvectorTPCalls"])
      qVectorTPCall(IsCalibrated, qvecReTPCallSp.at(0), qvecImTPCallSp.at(0), qvecAmp[kTPCall], TrkTPCallLabel);

    // Debug prints of values after corrections
    std::cout << "[CORRECTED] FT0C, Re: " << qvecReFT0CSp.at(0) << ", Im: " << qvecImFT0CSp.at(0) << std::endl;
    std::cout << "[CORRECTED] FT0A, Re: " << qvecReFT0ASp.at(0) << ", Im: " << qvecImFT0ASp.at(0) << std::endl;
    std::cout << "[CORRECTED] FT0M, Re: " << qvecReFT0MSp.at(0) << ", Im: " << qvecImFT0MSp.at(0) << std::endl;
    std::cout << "[CORRECTED] FV0A, Re: " << qvecReFV0ASp.at(0) << ", Im: " << qvecImFV0ASp.at(0) << std::endl;
    std::cout << "[CORRECTED] TPCpos, Re: " << qvecReTPCposSp.at(0) << ", Im: " << qvecImTPCposSp.at(0) << std::endl;
    std::cout << "[CORRECTED] TPCneg, Re: " << qvecReTPCnegSp.at(0) << ", Im: " << qvecImTPCnegSp.at(0) << std::endl;
    std::cout << "[CORRECTED] TPCall, Re: " << qvecReTPCallSp.at(0) << ", Im: " << qvecImTPCallSp.at(0) << std::endl;


    qVectorFT0CVec(IsCalibrated, qvecReFT0CSp, qvecImFT0CSp, qvecAmp[kFT0C]);
    qVectorFT0AVec(IsCalibrated, qvecReFT0ASp, qvecImFT0ASp, qvecAmp[kFT0A]);
    qVectorFT0MVec(IsCalibrated, qvecReFT0MSp, qvecImFT0MSp, qvecAmp[kFT0M]);
    qVectorFV0AVec(IsCalibrated, qvecReFV0ASp, qvecImFV0ASp, qvecAmp[kFV0A]);
    qVectorTPCposVec(IsCalibrated, qvecReTPCposSp, qvecImTPCposSp, qvecAmp[kTPCpos], TrkTPCposLabel);
    qVectorTPCnegVec(IsCalibrated, qvecReTPCnegSp, qvecImTPCnegSp, qvecAmp[kTPCneg], TrkTPCnegLabel);
    qVectorTPCallVec(IsCalibrated, qvecReTPCallSp, qvecImTPCallSp, qvecAmp[kTPCall], TrkTPCallLabel);

    // Deprecated, will be removed in future after transition time //
    if (useDetector["QvectorBPoss"])
      qVectorBPos(IsCalibrated, qvecReTPCposSp.at(0), qvecImTPCposSp.at(0), qvecAmp[kTPCpos], TrkTPCposLabel);
    if (useDetector["QvectorBNegs"])
      qVectorBNeg(IsCalibrated, qvecReTPCnegSp.at(0), qvecImTPCnegSp.at(0), qvecAmp[kTPCneg], TrkTPCnegLabel);
    if (useDetector["QvectorBTots"])
      qVectorBTot(IsCalibrated, qvecReTPCallSp.at(0), qvecImTPCallSp.at(0), qvecAmp[kTPCall], TrkTPCallLabel);

    qVectorBPosVec(IsCalibrated, qvecReTPCposSp, qvecImTPCposSp, qvecAmp[kTPCpos], TrkTPCposLabel);
    qVectorBNegVec(IsCalibrated, qvecReTPCnegSp, qvecImTPCnegSp, qvecAmp[kTPCneg], TrkTPCnegLabel);
    qVectorBTotVec(IsCalibrated, qvecReTPCallSp, qvecImTPCallSp, qvecAmp[kTPCall], TrkTPCallLabel);
    /////////////////////////////////////////////////////////////////

    if (cfgProduceRedQVecs) {
      eseQVector(cent, IsCalibrated, qvecReEse, qvecImEse, qvecAmp);
      eseQVectorFT0C(IsCalibrated, qvecReFT0CEse.at(0), qvecImFT0CEse.at(0), qvecAmp[kFT0C]);
      eseQVectorFT0A(IsCalibrated, qvecReFT0AEse.at(0), qvecImFT0AEse.at(0), qvecAmp[kFT0A]);
      eseQVectorFT0M(IsCalibrated, qvecReFT0MEse.at(0), qvecImFT0MEse.at(0), qvecAmp[kFT0M]);
      eseQVectorFV0A(IsCalibrated, qvecReFV0AEse.at(0), qvecImFV0AEse.at(0), qvecAmp[kFV0A]);
      eseQVectorTPCpos(IsCalibrated, qvecReTPCposEse.at(0), qvecImTPCposEse.at(0), qvecAmp[kTPCpos], TrkTPCposLabel);
      eseQVectorTPCneg(IsCalibrated, qvecReTPCnegEse.at(0), qvecImTPCnegEse.at(0), qvecAmp[kTPCneg], TrkTPCnegLabel);
      eseQVectorTPCall(IsCalibrated, qvecReTPCallEse.at(0), qvecImTPCallEse.at(0), qvecAmp[kTPCall], TrkTPCallLabel);
      eseQVectorFT0CVec(IsCalibrated, qvecReFT0CEse, qvecImFT0CEse, qvecAmp[kFT0C]);
      eseQVectorFT0AVec(IsCalibrated, qvecReFT0AEse, qvecImFT0AEse, qvecAmp[kFT0A]);
      eseQVectorFT0MVec(IsCalibrated, qvecReFT0MEse, qvecImFT0MEse, qvecAmp[kFT0M]);
      eseQVectorFV0AVec(IsCalibrated, qvecReFV0AEse, qvecImFV0AEse, qvecAmp[kFV0A]);
      eseQVectorTPCposVec(IsCalibrated, qvecReTPCposEse, qvecImTPCposEse, qvecAmp[kTPCpos], TrkTPCposLabel);
      eseQVectorTPCnegVec(IsCalibrated, qvecReTPCnegEse, qvecImTPCnegEse, qvecAmp[kTPCneg], TrkTPCnegLabel);
      eseQVectorTPCallVec(IsCalibrated, qvecReTPCallEse, qvecImTPCallEse, qvecAmp[kTPCall], TrkTPCallLabel);
      eseQVectorPerc(qvecReFT0CEse.at(0), qvecImFT0CEse.at(0), qvecAmp[kFT0C],
                     qvecReFT0AEse.at(0), qvecImFT0AEse.at(0), qvecAmp[kFT0A],
                     qvecReFT0MEse.at(0), qvecImFT0MEse.at(0), qvecAmp[kFT0M],
                     qvecReFV0AEse.at(0), qvecImFV0AEse.at(0), qvecAmp[kFV0A],
                     qvecReTPCposEse.at(0), qvecImTPCposEse.at(0), qvecAmp[kTPCpos],
                     qvecReTPCnegEse.at(0), qvecImTPCnegEse.at(0), qvecAmp[kTPCneg],
                     qvecReTPCallEse.at(0), qvecImTPCallEse.at(0), qvecAmp[kTPCall]);
    }
  } // End process.
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsTable>(cfgc)};
}
