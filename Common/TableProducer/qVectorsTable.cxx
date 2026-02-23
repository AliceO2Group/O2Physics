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
#include <TString.h>

#include <chrono>
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
    kTPCall
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
  Configurable<std::string> cfgQvecCalibPath{"cfgQvecCalibPath", "Analysis/EventPlane/QVecCorrections", "CCDB pasth for Q-vecteor calibration constants"};

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

  std::vector<TH3F*> objQvec{};
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

    objQvec.clear();
    for (std::size_t i = 0; i < cfgnMods->size(); i++) {
      int ind = cfgnMods->at(i);
      fullPath = cfgQvecCalibPath;
      fullPath += "/v";
      fullPath += std::to_string(ind);
      auto objqvec = getForTsOrRun<TH3F>(fullPath, timestamp, runnumber);
      if (!objqvec) {
        fullPath = cfgQvecCalibPath;
        fullPath += "/v2";
        objqvec = getForTsOrRun<TH3F>(fullPath, timestamp, runnumber);
      }
      objQvec.push_back(objqvec);
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

  template <typename Nmode, typename CollType, typename TrackType>
  void CalQvec(const Nmode nmode, const CollType& coll, const TrackType& track, std::vector<float>& QvecRe, std::vector<float>& QvecIm, std::vector<float>& QvecAmp, std::vector<int>& TrkTPCposLabel, std::vector<int>& TrkTPCnegLabel, std::vector<int>& TrkTPCallLabel)
  {
    float qVectFT0A[2] = {0.};
    float qVectFT0C[2] = {0.};
    float qVectFT0M[2] = {0.};
    float qVectFV0A[2] = {0.};
    float qVectTPCpos[2] = {0.};
    float qVectTPCneg[2] = {0.};
    float qVectTPCall[2] = {0.};

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

          helperEP.SumQvectors(0, FT0AchId, ampl / FT0RelGainConst[FT0AchId], nmode, QvecDet, sumAmplFT0A, ft0geom, fv0geom);
          helperEP.SumQvectors(0, FT0AchId, ampl / FT0RelGainConst[FT0AchId], nmode, QvecFT0M, sumAmplFT0M, ft0geom, fv0geom);
        }
        if (sumAmplFT0A > 1e-8) {
          QvecDet /= sumAmplFT0A;
          qVectFT0A[0] = QvecDet.Re();
          qVectFT0A[1] = QvecDet.Im();
        }
      } else {
        qVectFT0A[0] = 999.;
        qVectFT0A[1] = 999.;
      }

      if (useDetector["QvectorFT0Cs"]) {
        QvecDet = TComplex(0., 0.);
        for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
          float ampl = ft0.amplitudeC()[iChC];
          int FT0CchId = ft0.channelC()[iChC] + 96;

          histosQA.fill(HIST("FT0Amp"), ampl, FT0CchId);
          histosQA.fill(HIST("FT0AmpCor"), ampl / FT0RelGainConst[FT0CchId], FT0CchId);

          helperEP.SumQvectors(0, FT0CchId, ampl / FT0RelGainConst[FT0CchId], nmode, QvecDet, sumAmplFT0C, ft0geom, fv0geom);
          helperEP.SumQvectors(0, FT0CchId, ampl / FT0RelGainConst[FT0CchId], nmode, QvecFT0M, sumAmplFT0M, ft0geom, fv0geom);
        }

        if (sumAmplFT0C > 1e-8) {
          QvecDet /= sumAmplFT0C;
          qVectFT0C[0] = QvecDet.Re();
          qVectFT0C[1] = QvecDet.Im();
        } else {
          qVectFT0C[0] = 999.;
          qVectFT0C[1] = 999.;
        }
      } else {
        qVectFT0C[0] = -999.;
        qVectFT0C[1] = -999.;
      }

      if (sumAmplFT0M > 1e-8 && useDetector["QvectorFT0Ms"]) {
        QvecFT0M /= sumAmplFT0M;
        qVectFT0M[0] = QvecFT0M.Re();
        qVectFT0M[1] = QvecFT0M.Im();
      } else {
        qVectFT0M[0] = 999.;
        qVectFT0M[1] = 999.;
      }
    } else {
      qVectFT0A[0] = -999.;
      qVectFT0A[1] = -999.;
      qVectFT0C[0] = -999.;
      qVectFT0C[1] = -999.;
      qVectFT0M[0] = -999.;
      qVectFT0M[1] = -999.;
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

        helperEP.SumQvectors(1, FV0AchId, ampl / FV0RelGainConst[FV0AchId], nmode, QvecDet, sumAmplFV0A, ft0geom, fv0geom);
      }

      if (sumAmplFV0A > 1e-8) {
        QvecDet /= sumAmplFV0A;
        qVectFV0A[0] = QvecDet.Re();
        qVectFV0A[1] = QvecDet.Im();
      } else {
        qVectFV0A[0] = 999.;
        qVectFV0A[1] = 999.;
      }
    } else {
      qVectFV0A[0] = -999.;
      qVectFV0A[1] = -999.;
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
      qVectTPCall[0] += trk.pt() * std::cos(trk.phi() * nmode);
      qVectTPCall[1] += trk.pt() * std::sin(trk.phi() * nmode);
      TrkTPCallLabel.push_back(trk.globalIndex());
      nTrkTPCall++;
      if (std::abs(trk.eta()) < 0.1) {
        continue;
      }
      if (trk.eta() > 0 && (useDetector["QvectorTPCposs"] || useDetector["QvectorBPoss"])) {
        qVectTPCpos[0] += trk.pt() * std::cos(trk.phi() * nmode);
        qVectTPCpos[1] += trk.pt() * std::sin(trk.phi() * nmode);
        TrkTPCposLabel.push_back(trk.globalIndex());
        nTrkTPCpos++;
      } else if (trk.eta() < 0 && (useDetector["QvectorTPCnegs"] || useDetector["QvectorBNegs"])) {
        qVectTPCneg[0] += trk.pt() * std::cos(trk.phi() * nmode);
        qVectTPCneg[1] += trk.pt() * std::sin(trk.phi() * nmode);
        TrkTPCnegLabel.push_back(trk.globalIndex());
        nTrkTPCneg++;
      }
    }
    if (nTrkTPCpos > 0) {
      qVectTPCpos[0] /= nTrkTPCpos;
      qVectTPCpos[1] /= nTrkTPCpos;
    } else {
      qVectTPCpos[0] = 999.;
      qVectTPCpos[1] = 999.;
    }

    if (nTrkTPCneg > 0) {
      qVectTPCneg[0] /= nTrkTPCneg;
      qVectTPCneg[1] /= nTrkTPCneg;
    } else {
      qVectTPCneg[0] = 999.;
      qVectTPCneg[1] = 999.;
    }

    if (nTrkTPCall > 0) {
      qVectTPCall[0] /= nTrkTPCall;
      qVectTPCall[1] /= nTrkTPCall;
    } else {
      qVectTPCall[0] = 999.;
      qVectTPCall[1] = 999.;
    }

    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectFT0C[0]);
      QvecIm.push_back(qVectFT0C[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectFT0A[0]);
      QvecIm.push_back(qVectFT0A[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectFT0M[0]);
      QvecIm.push_back(qVectFT0M[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectFV0A[0]);
      QvecIm.push_back(qVectFV0A[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectTPCpos[0]);
      QvecIm.push_back(qVectTPCpos[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectTPCneg[0]);
      QvecIm.push_back(qVectTPCneg[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectTPCall[0]);
      QvecIm.push_back(qVectTPCall[1]);
    }

    QvecAmp.push_back(sumAmplFT0C);
    QvecAmp.push_back(sumAmplFT0A);
    QvecAmp.push_back(sumAmplFT0M);
    QvecAmp.push_back(sumAmplFV0A);
    QvecAmp.push_back(static_cast<float>(nTrkTPCpos));
    QvecAmp.push_back(static_cast<float>(nTrkTPCneg));
    QvecAmp.push_back(static_cast<float>(nTrkTPCall));
  }

  void process(MyCollisions::iterator const& coll, aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&, MyTracks const& tracks)
  {
    std::vector<int> TrkTPCposLabel{};
    std::vector<int> TrkTPCnegLabel{};
    std::vector<int> TrkTPCallLabel{};
    std::vector<float> qvecRe{};
    std::vector<float> qvecIm{};
    std::vector<float> qvecAmp{};

    std::vector<float> qvecReFT0C{};
    std::vector<float> qvecImFT0C{};
    std::vector<float> qvecReFT0A{};
    std::vector<float> qvecImFT0A{};
    std::vector<float> qvecReFT0M{};
    std::vector<float> qvecImFT0M{};
    std::vector<float> qvecReFV0A{};
    std::vector<float> qvecImFV0A{};
    std::vector<float> qvecReTPCpos{};
    std::vector<float> qvecImTPCpos{};
    std::vector<float> qvecReTPCneg{};
    std::vector<float> qvecImTPCneg{};
    std::vector<float> qvecReTPCall{};
    std::vector<float> qvecImTPCall{};

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
      int nmode = cfgnMods->at(id);
      CalQvec(nmode, coll, tracks, qvecRe, qvecIm, qvecAmp, TrkTPCposLabel, TrkTPCnegLabel, TrkTPCallLabel);
      if (cent < cfgMaxCentrality) {
        for (auto i{0u}; i < kTPCall + 1; i++) {
          helperEP.DoRecenter(qvecRe[(kTPCall + 1) * 4 * id + i * 4 + 1], qvecIm[(kTPCall + 1) * 4 * id + i * 4 + 1],
                              objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));

          helperEP.DoRecenter(qvecRe[(kTPCall + 1) * 4 * id + i * 4 + 2], qvecIm[(kTPCall + 1) * 4 * id + i * 4 + 2],
                              objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));
          helperEP.DoTwist(qvecRe[(kTPCall + 1) * 4 * id + i * 4 + 2], qvecIm[(kTPCall + 1) * 4 * id + i * 4 + 2],
                           objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 3, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 4, i + 1));

          helperEP.DoRecenter(qvecRe[(kTPCall + 1) * 4 * id + i * 4 + 3], qvecIm[(kTPCall + 1) * 4 * id + i * 4 + 3],
                              objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));
          helperEP.DoTwist(qvecRe[(kTPCall + 1) * 4 * id + i * 4 + 3], qvecIm[(kTPCall + 1) * 4 * id + i * 4 + 3],
                           objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 3, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 4, i + 1));
          helperEP.DoRescale(qvecRe[(kTPCall + 1) * 4 * id + i * 4 + 3], qvecIm[(kTPCall + 1) * 4 * id + i * 4 + 3],
                             objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 5, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 6, i + 1));
        }
        if (cfgShiftCorr) {
          auto deltapsiFT0C = 0.0;
          auto deltapsiFT0A = 0.0;
          auto deltapsiFT0M = 0.0;
          auto deltapsiFV0A = 0.0;
          auto deltapsiTPCpos = 0.0;
          auto deltapsiTPCneg = 0.0;
          auto deltapsiTPCall = 0.0;

          auto psidefFT0C = TMath::ATan2(qvecIm[(kTPCall + 1) * 4 * id + kFT0C * 4 + 3], qvecRe[(kTPCall + 1) * 4 * id + kFT0C * 4 + 3]) / static_cast<float>(nmode);
          auto psidefFT0A = TMath::ATan2(qvecIm[(kTPCall + 1) * 4 * id + kFT0A * 4 + 3], qvecRe[(kTPCall + 1) * 4 * id + kFT0A * 4 + 3]) / static_cast<float>(nmode);
          auto psidefFT0M = TMath::ATan2(qvecIm[(kTPCall + 1) * 4 * id + kFT0M * 4 + 3], qvecRe[(kTPCall + 1) * 4 * id + kFT0M * 4 + 3]) / static_cast<float>(nmode);
          auto psidefFV0A = TMath::ATan2(qvecIm[(kTPCall + 1) * 4 * id + kFV0A * 4 + 3], qvecRe[(kTPCall + 1) * 4 * id + kFV0A * 4 + 3]) / static_cast<float>(nmode);
          auto psidefTPCpos = TMath::ATan2(qvecIm[(kTPCall + 1) * 4 * id + kTPCpos * 4 + 3], qvecRe[(kTPCall + 1) * 4 * id + kTPCpos * 4 + 3]) / static_cast<float>(nmode);
          auto psidefTPCneg = TMath::ATan2(qvecIm[(kTPCall + 1) * 4 * id + kTPCneg * 4 + 3], qvecRe[(kTPCall + 1) * 4 * id + kTPCneg * 4 + 3]) / static_cast<float>(nmode);
          auto psidefTPCall = TMath::ATan2(qvecIm[(kTPCall + 1) * 4 * id + kTPCall * 4 + 3], qvecRe[(kTPCall + 1) * 4 * id + kTPCall * 4 + 3]) / static_cast<float>(nmode);

          for (int ishift = 1; ishift <= 10; ishift++) {
            auto coeffshiftxFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kFT0C, ishift - 0.5));
            auto coeffshiftyFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kFT0C + 1, ishift - 0.5));
            auto coeffshiftxFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kFT0A, ishift - 0.5));
            auto coeffshiftyFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kFT0A + 1, ishift - 0.5));
            auto coeffshiftxFT0M = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kFT0M, ishift - 0.5));
            auto coeffshiftyFT0M = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kFT0M + 1, ishift - 0.5));
            auto coeffshiftxFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kFV0A, ishift - 0.5));
            auto coeffshiftyFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kFV0A + 1, ishift - 0.5));
            auto coeffshiftxTPCpos = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kTPCpos, ishift - 0.5));
            auto coeffshiftyTPCpos = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kTPCpos + 1, ishift - 0.5));
            auto coeffshiftxTPCneg = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kTPCneg, ishift - 0.5));
            auto coeffshiftyTPCneg = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kTPCneg + 1, ishift - 0.5));
            auto coeffshiftxTPCall = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kTPCall, ishift - 0.5));
            auto coeffshiftyTPCall = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(cent, 2 * kTPCall + 1, ishift - 0.5));

            deltapsiFT0C += ((2. / (1.0 * ishift)) * (-coeffshiftxFT0C * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0C) + coeffshiftyFT0C * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0C))) / static_cast<float>(nmode);
            deltapsiFT0A += ((2. / (1.0 * ishift)) * (-coeffshiftxFT0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0A) + coeffshiftyFT0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0A))) / static_cast<float>(nmode);
            deltapsiFT0M += ((2. / (1.0 * ishift)) * (-coeffshiftxFT0M * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0M) + coeffshiftyFT0M * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0M))) / static_cast<float>(nmode);
            deltapsiFV0A += ((2. / (1.0 * ishift)) * (-coeffshiftxFV0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFV0A) + coeffshiftyFV0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFV0A))) / static_cast<float>(nmode);
            deltapsiTPCpos += ((2. / (1.0 * ishift)) * (-coeffshiftxTPCpos * TMath::Cos(ishift * static_cast<float>(nmode) * psidefTPCpos) + coeffshiftyTPCpos * TMath::Sin(ishift * static_cast<float>(nmode) * psidefTPCpos))) / static_cast<float>(nmode);
            deltapsiTPCneg += ((2. / (1.0 * ishift)) * (-coeffshiftxTPCneg * TMath::Cos(ishift * static_cast<float>(nmode) * psidefTPCneg) + coeffshiftyTPCneg * TMath::Sin(ishift * static_cast<float>(nmode) * psidefTPCneg))) / static_cast<float>(nmode);
            deltapsiTPCall += ((2. / (1.0 * ishift)) * (-coeffshiftxTPCall * TMath::Cos(ishift * static_cast<float>(nmode) * psidefTPCall) + coeffshiftyTPCall * TMath::Sin(ishift * static_cast<float>(nmode) * psidefTPCall))) / static_cast<float>(nmode);
          }

          deltapsiFT0C *= static_cast<float>(nmode);
          deltapsiFT0A *= static_cast<float>(nmode);
          deltapsiFT0M *= static_cast<float>(nmode);
          deltapsiFV0A *= static_cast<float>(nmode);
          deltapsiTPCpos *= static_cast<float>(nmode);
          deltapsiTPCneg *= static_cast<float>(nmode);
          deltapsiTPCall *= static_cast<float>(nmode);

          float qvecReShiftedFT0C = qvecRe[(kTPCall + 1) * 4 * id + kFT0C * 4 + 3] * TMath::Cos(deltapsiFT0C) - qvecIm[(kTPCall + 1) * 4 * id + kFT0C * 4 + 3] * TMath::Sin(deltapsiFT0C);
          float qvecImShiftedFT0C = qvecRe[(kTPCall + 1) * 4 * id + kFT0C * 4 + 3] * TMath::Sin(deltapsiFT0C) + qvecIm[(kTPCall + 1) * 4 * id + kFT0C * 4 + 3] * TMath::Cos(deltapsiFT0C);
          float qvecReShiftedFT0A = qvecRe[(kTPCall + 1) * 4 * id + kFT0A * 4 + 3] * TMath::Cos(deltapsiFT0A) - qvecIm[(kTPCall + 1) * 4 * id + kFT0A * 4 + 3] * TMath::Sin(deltapsiFT0A);
          float qvecImShiftedFT0A = qvecRe[(kTPCall + 1) * 4 * id + kFT0A * 4 + 3] * TMath::Sin(deltapsiFT0A) + qvecIm[(kTPCall + 1) * 4 * id + kFT0A * 4 + 3] * TMath::Cos(deltapsiFT0A);
          float qvecReShiftedFT0M = qvecRe[(kTPCall + 1) * 4 * id + kFT0M * 4 + 3] * TMath::Cos(deltapsiFT0M) - qvecIm[(kTPCall + 1) * 4 * id + kFT0M * 4 + 3] * TMath::Sin(deltapsiFT0M);
          float qvecImShiftedFT0M = qvecRe[(kTPCall + 1) * 4 * id + kFT0M * 4 + 3] * TMath::Sin(deltapsiFT0M) + qvecIm[(kTPCall + 1) * 4 * id + kFT0M * 4 + 3] * TMath::Cos(deltapsiFT0M);
          float qvecReShiftedFV0A = qvecRe[(kTPCall + 1) * 4 * id + kFV0A * 4 + 3] * TMath::Cos(deltapsiFV0A) - qvecIm[(kTPCall + 1) * 4 * id + kFV0A * 4 + 3] * TMath::Sin(deltapsiFV0A);
          float qvecImShiftedFV0A = qvecRe[(kTPCall + 1) * 4 * id + kFV0A * 4 + 3] * TMath::Sin(deltapsiFV0A) + qvecIm[(kTPCall + 1) * 4 * id + kFV0A * 4 + 3] * TMath::Cos(deltapsiFV0A);
          float qvecReShiftedTPCpos = qvecRe[(kTPCall + 1) * 4 * id + kTPCpos * 4 + 3] * TMath::Cos(deltapsiTPCpos) - qvecIm[(kTPCall + 1) * 4 * id + kTPCpos * 4 + 3] * TMath::Sin(deltapsiTPCpos);
          float qvecImShiftedTPCpos = qvecRe[(kTPCall + 1) * 4 * id + kTPCpos * 4 + 3] * TMath::Sin(deltapsiTPCpos) + qvecIm[(kTPCall + 1) * 4 * id + kTPCpos * 4 + 3] * TMath::Cos(deltapsiTPCpos);
          float qvecReShiftedTPCneg = qvecRe[(kTPCall + 1) * 4 * id + kTPCneg * 4 + 3] * TMath::Cos(deltapsiTPCneg) - qvecIm[(kTPCall + 1) * 4 * id + kTPCneg * 4 + 3] * TMath::Sin(deltapsiTPCneg);
          float qvecImShiftedTPCneg = qvecRe[(kTPCall + 1) * 4 * id + kTPCneg * 4 + 3] * TMath::Sin(deltapsiTPCneg) + qvecIm[(kTPCall + 1) * 4 * id + kTPCneg * 4 + 3] * TMath::Cos(deltapsiTPCneg);
          float qvecReShiftedTPCall = qvecRe[(kTPCall + 1) * 4 * id + kTPCall * 4 + 3] * TMath::Cos(deltapsiTPCall) - qvecIm[(kTPCall + 1) * 4 * id + kTPCall * 4 + 3] * TMath::Sin(deltapsiTPCall);
          float qvecImShiftedTPCall = qvecRe[(kTPCall + 1) * 4 * id + kTPCall * 4 + 3] * TMath::Sin(deltapsiTPCall) + qvecIm[(kTPCall + 1) * 4 * id + kTPCall * 4 + 3] * TMath::Cos(deltapsiTPCall);

          qvecRe[(kTPCall + 1) * 4 * id + kFT0C * 4 + 3] = qvecReShiftedFT0C;
          qvecIm[(kTPCall + 1) * 4 * id + kFT0C * 4 + 3] = qvecImShiftedFT0C;
          qvecRe[(kTPCall + 1) * 4 * id + kFT0A * 4 + 3] = qvecReShiftedFT0A;
          qvecIm[(kTPCall + 1) * 4 * id + kFT0A * 4 + 3] = qvecImShiftedFT0A;
          qvecRe[(kTPCall + 1) * 4 * id + kFT0M * 4 + 3] = qvecReShiftedFT0M;
          qvecIm[(kTPCall + 1) * 4 * id + kFT0M * 4 + 3] = qvecImShiftedFT0M;
          qvecRe[(kTPCall + 1) * 4 * id + kFV0A * 4 + 3] = qvecReShiftedFV0A;
          qvecIm[(kTPCall + 1) * 4 * id + kFV0A * 4 + 3] = qvecImShiftedFV0A;
          qvecRe[(kTPCall + 1) * 4 * id + kTPCpos * 4 + 3] = qvecReShiftedTPCpos;
          qvecIm[(kTPCall + 1) * 4 * id + kTPCpos * 4 + 3] = qvecImShiftedTPCpos;
          qvecRe[(kTPCall + 1) * 4 * id + kTPCneg * 4 + 3] = qvecReShiftedTPCneg;
          qvecIm[(kTPCall + 1) * 4 * id + kTPCneg * 4 + 3] = qvecImShiftedTPCneg;
          qvecRe[(kTPCall + 1) * 4 * id + kTPCall * 4 + 3] = qvecReShiftedTPCall;
          qvecIm[(kTPCall + 1) * 4 * id + kTPCall * 4 + 3] = qvecImShiftedTPCall;
        }
      }
      int CorrLevel = cfgCorrLevel == 0 ? 0 : cfgCorrLevel - 1;
      qvecReFT0C.push_back(qvecRe[(kTPCall + 1) * 4 * id + kFT0C * 4 + CorrLevel]);
      qvecImFT0C.push_back(qvecIm[(kTPCall + 1) * 4 * id + kFT0C * 4 + CorrLevel]);
      qvecReFT0A.push_back(qvecRe[(kTPCall + 1) * 4 * id + kFT0A * 4 + CorrLevel]);
      qvecImFT0A.push_back(qvecIm[(kTPCall + 1) * 4 * id + kFT0A * 4 + CorrLevel]);
      qvecReFT0M.push_back(qvecRe[(kTPCall + 1) * 4 * id + kFT0M * 4 + CorrLevel]);
      qvecImFT0M.push_back(qvecIm[(kTPCall + 1) * 4 * id + kFT0M * 4 + CorrLevel]);
      qvecReFV0A.push_back(qvecRe[(kTPCall + 1) * 4 * id + kFV0A * 4 + CorrLevel]);
      qvecImFV0A.push_back(qvecIm[(kTPCall + 1) * 4 * id + kFV0A * 4 + CorrLevel]);
      qvecReTPCpos.push_back(qvecRe[(kTPCall + 1) * 4 * id + kTPCpos * 4 + CorrLevel]);
      qvecImTPCpos.push_back(qvecIm[(kTPCall + 1) * 4 * id + kTPCpos * 4 + CorrLevel]);
      qvecReTPCneg.push_back(qvecRe[(kTPCall + 1) * 4 * id + kTPCneg * 4 + CorrLevel]);
      qvecImTPCneg.push_back(qvecIm[(kTPCall + 1) * 4 * id + kTPCneg * 4 + CorrLevel]);
      qvecReTPCall.push_back(qvecRe[(kTPCall + 1) * 4 * id + kTPCall * 4 + CorrLevel]);
      qvecImTPCall.push_back(qvecIm[(kTPCall + 1) * 4 * id + kTPCall * 4 + CorrLevel]);
    }

    // Fill the columns of the Qvectors table.
    qVector(cent, IsCalibrated, qvecRe, qvecIm, qvecAmp);
    if (useDetector["QvectorFT0Cs"])
      qVectorFT0C(IsCalibrated, qvecReFT0C.at(0), qvecImFT0C.at(0), qvecAmp[kFT0C]);
    if (useDetector["QvectorFT0As"])
      qVectorFT0A(IsCalibrated, qvecReFT0A.at(0), qvecImFT0A.at(0), qvecAmp[kFT0A]);
    if (useDetector["QvectorFT0Ms"])
      qVectorFT0M(IsCalibrated, qvecReFT0M.at(0), qvecImFT0M.at(0), qvecAmp[kFT0M]);
    if (useDetector["QvectorFV0As"])
      qVectorFV0A(IsCalibrated, qvecReFV0A.at(0), qvecImFV0A.at(0), qvecAmp[kFV0A]);
    if (useDetector["QvectorTPCposs"])
      qVectorTPCpos(IsCalibrated, qvecReTPCpos.at(0), qvecImTPCpos.at(0), qvecAmp[kTPCpos], TrkTPCposLabel);
    if (useDetector["QvectorTPCnegs"])
      qVectorTPCneg(IsCalibrated, qvecReTPCneg.at(0), qvecImTPCneg.at(0), qvecAmp[kTPCneg], TrkTPCnegLabel);
    if (useDetector["QvectorTPCalls"])
      qVectorTPCall(IsCalibrated, qvecReTPCall.at(0), qvecImTPCall.at(0), qvecAmp[kTPCall], TrkTPCallLabel);

    qVectorFT0CVec(IsCalibrated, qvecReFT0C, qvecImFT0C, qvecAmp[kFT0C]);
    qVectorFT0AVec(IsCalibrated, qvecReFT0A, qvecImFT0A, qvecAmp[kFT0A]);
    qVectorFT0MVec(IsCalibrated, qvecReFT0M, qvecImFT0M, qvecAmp[kFT0M]);
    qVectorFV0AVec(IsCalibrated, qvecReFV0A, qvecImFV0A, qvecAmp[kFV0A]);
    qVectorTPCposVec(IsCalibrated, qvecReTPCpos, qvecImTPCpos, qvecAmp[kTPCpos], TrkTPCposLabel);
    qVectorTPCnegVec(IsCalibrated, qvecReTPCneg, qvecImTPCneg, qvecAmp[kTPCneg], TrkTPCnegLabel);
    qVectorTPCallVec(IsCalibrated, qvecReTPCall, qvecImTPCall, qvecAmp[kTPCall], TrkTPCallLabel);

    // Deprecated, will be removed in future after transition time //
    if (useDetector["QvectorBPoss"])
      qVectorBPos(IsCalibrated, qvecReTPCpos.at(0), qvecImTPCpos.at(0), qvecAmp[kTPCpos], TrkTPCposLabel);
    if (useDetector["QvectorBNegs"])
      qVectorBNeg(IsCalibrated, qvecReTPCneg.at(0), qvecImTPCneg.at(0), qvecAmp[kTPCneg], TrkTPCnegLabel);
    if (useDetector["QvectorBTots"])
      qVectorBTot(IsCalibrated, qvecReTPCall.at(0), qvecImTPCall.at(0), qvecAmp[kTPCall], TrkTPCallLabel);

    qVectorBPosVec(IsCalibrated, qvecReTPCpos, qvecImTPCpos, qvecAmp[kTPCpos], TrkTPCposLabel);
    qVectorBNegVec(IsCalibrated, qvecReTPCneg, qvecImTPCneg, qvecAmp[kTPCneg], TrkTPCnegLabel);
    qVectorBTotVec(IsCalibrated, qvecReTPCall, qvecImTPCall, qvecAmp[kTPCall], TrkTPCallLabel);
    /////////////////////////////////////////////////////////////////

  } // End process.
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsTable>(cfgc)};
}
