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

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TH3F.h>

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/Qvectors.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
// o2 includes.
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

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
    kBPos,
    kBNeg,
    kBTot
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
  Configurable<int> cfgCorrLevel{"cfgCorrLevel", 4, "calibration step: 0 = no corr, 1 = gain corr, 2 = rectr, 3 = twist, 4 = full"};
  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2, 3}, "Modulation of interest"};

  Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Users/j/junlee/Qvector/GainEq", "CCDB path for gain equalization constants"};
  Configurable<std::string> cfgQvecCalibPath{"cfgQvecCalibPath", "Analysis/EventPlane/QVecCorrections", "CCDB pasth for Q-vecteor calibration constants"};

  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, ""};

  Configurable<bool> cfgUseFT0C{"cfgUseFT0C", false, "Initial value for using FT0C. By default obtained from DataModel."};
  Configurable<bool> cfgUseFT0A{"cfgUseFT0A", false, "Initial value for using FT0A. By default obtained from DataModel."};
  Configurable<bool> cfgUseFT0M{"cfgUseFT0M", false, "Initial value for using FT0M. By default obtained from DataModel."};
  Configurable<bool> cfgUseFV0A{"cfgUseFV0A", false, "Initial value for using FV0A. By default obtained from DataModel."};
  Configurable<bool> cfgUseBPos{"cfgUseBPos", false, "Initial value for using BPos. By default obtained from DataModel."};
  Configurable<bool> cfgUseBNeg{"cfgUseBNeg", false, "Initial value for using BNeg. By default obtained from DataModel."};
  Configurable<bool> cfgUseBTot{"cfgUseBTot", false, "Initial value for using BTot. By default obtained from DataModel."};

  // Table.
  Produces<aod::Qvectors> qVector;
  Produces<aod::QvectorFT0Cs> qVectorFT0C;
  Produces<aod::QvectorFT0As> qVectorFT0A;
  Produces<aod::QvectorFT0Ms> qVectorFT0M;
  Produces<aod::QvectorFV0As> qVectorFV0A;
  Produces<aod::QvectorBPoss> qVectorBPos;
  Produces<aod::QvectorBNegs> qVectorBNeg;
  Produces<aod::QvectorBTots> qVectorBTot;

  Produces<aod::QvectorFT0CVecs> qVectorFT0CVec;
  Produces<aod::QvectorFT0AVecs> qVectorFT0AVec;
  Produces<aod::QvectorFT0MVecs> qVectorFT0MVec;
  Produces<aod::QvectorFV0AVecs> qVectorFV0AVec;
  Produces<aod::QvectorBPosVecs> qVectorBPosVec;
  Produces<aod::QvectorBNegVecs> qVectorBNegVec;
  Produces<aod::QvectorBTotVecs> qVectorBTotVec;

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

  std::unordered_map<string, bool> useDetector = {
    {"QvectorBTots", cfgUseBTot},
    {"QvectorBNegs", cfgUseBNeg},
    {"QvectorBPoss", cfgUseBPos},
    {"QvectorFV0As", cfgUseFV0A},
    {"QvectorFT0Ms", cfgUseFT0M},
    {"QvectorFT0As", cfgUseFT0A},
    {"QvectorFT0Cs", cfgUseFT0C}};

  void init(InitContext& initContext)
  {
    // Check the sub-detector used
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      for (auto const& input : device.inputs) {
        if (input.matcher.binding == "Qvectors") {
          for (auto det : useDetector) {
            useDetector[det.first.data()] = true;
          }
          LOGF(info, "Using all detectors.");
          goto allDetectorsInUse; // Added to break from nested loop if all detectors are in use.
        }
        for (auto det : useDetector) {
          if (input.matcher.binding == det.first) {
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

    auto offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    auto offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", timestamp);

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
    for (auto i = 0; i < cfgnMods->size(); i++) {
      int ind = cfgnMods->at(i);
      fullPath = cfgQvecCalibPath;
      fullPath += "/v";
      fullPath += std::to_string(ind);
      auto objqvec = ccdb->getForTimeStamp<TH3F>(fullPath, timestamp);
      if (!objqvec) {
        fullPath = cfgQvecCalibPath;
        fullPath += "/v2";
        objqvec = ccdb->getForTimeStamp<TH3F>(fullPath, timestamp);
      }
      objQvec.push_back(objqvec);
    }
    fullPath = cfgGainEqPath;
    fullPath += "/FT0";
    auto objft0Gain = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
    if (!objft0Gain || cfgCorrLevel == 0) {
      for (auto i{0u}; i < 208; i++) {
        FT0RelGainConst.push_back(1.);
      }
    } else {
      FT0RelGainConst = *(objft0Gain);
    }

    fullPath = cfgGainEqPath;
    fullPath += "/FV0";
    auto objfv0Gain = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
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

  template <typename Nmode, typename CollType, typename TrackType>
  void CalQvec(const Nmode nmode, const CollType& coll, const TrackType& track, std::vector<float>& QvecRe, std::vector<float>& QvecIm, std::vector<float>& QvecAmp, std::vector<int>& TrkBPosLabel, std::vector<int>& TrkBNegLabel, std::vector<int>& TrkBTotLabel)
  {
    float qVectFT0A[2] = {0.};
    float qVectFT0C[2] = {0.};
    float qVectFT0M[2] = {0.};
    float qVectFV0A[2] = {0.};
    float qVectBPos[2] = {0.};
    float qVectBNeg[2] = {0.};
    float qVectBTot[2] = {0.};

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

    int nTrkBPos = 0;
    int nTrkBNeg = 0;
    int nTrkBTot = 0;

    for (auto& trk : track) {
      if (!SelTrack(trk)) {
        continue;
      }
      histosQA.fill(HIST("ChTracks"), trk.pt(), trk.eta(), trk.phi(), cent);
      if (std::abs(trk.eta()) < 0.1 || std::abs(trk.eta()) > 0.8) {
        continue;
      }
      if (trk.eta() > 0 && useDetector["QvectorBPoss"]) {
        qVectBPos[0] += trk.pt() * std::cos(trk.phi() * nmode);
        qVectBPos[1] += trk.pt() * std::sin(trk.phi() * nmode);
        TrkBPosLabel.push_back(trk.globalIndex());
        nTrkBPos++;
      } else if (trk.eta() < 0 && useDetector["QvectorBNegs"]) {
        qVectBNeg[0] += trk.pt() * std::cos(trk.phi() * nmode);
        qVectBNeg[1] += trk.pt() * std::sin(trk.phi() * nmode);
        TrkBNegLabel.push_back(trk.globalIndex());
        nTrkBNeg++;
      }
      qVectBTot[0] += trk.pt() * std::cos(trk.phi() * nmode);
      qVectBTot[1] += trk.pt() * std::sin(trk.phi() * nmode);
      TrkBTotLabel.push_back(trk.globalIndex());
      nTrkBTot++;
    }
    if (nTrkBPos > 0) {
      qVectBPos[0] /= nTrkBPos;
      qVectBPos[1] /= nTrkBPos;
    } else {
      qVectBPos[0] = 999.;
      qVectBPos[1] = 999.;
    }

    if (nTrkBNeg > 0) {
      qVectBNeg[0] /= nTrkBNeg;
      qVectBNeg[1] /= nTrkBNeg;
    } else {
      qVectBNeg[0] = 999.;
      qVectBNeg[1] = 999.;
    }

    if (nTrkBTot > 0) {
      qVectBTot[0] /= nTrkBTot;
      qVectBTot[1] /= nTrkBTot;
    } else {
      qVectBTot[0] = 999.;
      qVectBTot[1] = 999.;
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
      QvecRe.push_back(qVectBPos[0]);
      QvecIm.push_back(qVectBPos[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectBNeg[0]);
      QvecIm.push_back(qVectBNeg[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      QvecRe.push_back(qVectBTot[0]);
      QvecIm.push_back(qVectBTot[1]);
    }

    QvecAmp.push_back(sumAmplFT0C);
    QvecAmp.push_back(sumAmplFT0A);
    QvecAmp.push_back(sumAmplFT0M);
    QvecAmp.push_back(sumAmplFV0A);
    QvecAmp.push_back(static_cast<float>(nTrkBPos));
    QvecAmp.push_back(static_cast<float>(nTrkBNeg));
    QvecAmp.push_back(static_cast<float>(nTrkBTot));
  }

  void process(MyCollisions::iterator const& coll, aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&, MyTracks const& tracks)
  {
    std::vector<int> TrkBPosLabel{};
    std::vector<int> TrkBNegLabel{};
    std::vector<int> TrkBTotLabel{};
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
    std::vector<float> qvecReBPos{};
    std::vector<float> qvecImBPos{};
    std::vector<float> qvecReBNeg{};
    std::vector<float> qvecImBNeg{};
    std::vector<float> qvecReBTot{};
    std::vector<float> qvecImBTot{};

    auto bc = coll.bc_as<aod::BCsWithTimestamps>();
    int currentRun = bc.runNumber();
    if (runNumber != currentRun) {
      initCCDB(bc);
      runNumber = currentRun;
    }

    float centAllEstim[4] = {
      coll.centFT0M(), coll.centFT0A(), coll.centFT0C(),
      coll.centFV0A()};
    cent = centAllEstim[cfgCentEsti];
    bool IsCalibrated = true;
    if (cent < 0. || cent > 80.) {
      cent = 110.;
      IsCalibrated = false;
    }
    for (auto id = 0; id < cfgnMods->size(); id++) {
      int ind = cfgnMods->at(id);
      CalQvec(ind, coll, tracks, qvecRe, qvecIm, qvecAmp, TrkBPosLabel, TrkBNegLabel, TrkBTotLabel);
      if (cent < 80) {
        for (auto i{0u}; i < 6; i++) {
          helperEP.DoRecenter(qvecRe[(kBTot + 1) * 4 * id + i * 4 + 1], qvecIm[(kBTot + 1) * 4 * id + i * 4 + 1],
                              objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));

          helperEP.DoRecenter(qvecRe[(kBTot + 1) * 4 * id + i * 4 + 2], qvecIm[(kBTot + 1) * 4 * id + i * 4 + 2],
                              objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));
          helperEP.DoTwist(qvecRe[(kBTot + 1) * 4 * id + i * 4 + 2], qvecIm[(kBTot + 1) * 4 * id + i * 4 + 2],
                           objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 3, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 4, i + 1));

          helperEP.DoRecenter(qvecRe[(kBTot + 1) * 4 * id + i * 4 + 3], qvecIm[(kBTot + 1) * 4 * id + i * 4 + 3],
                              objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 1, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 2, i + 1));
          helperEP.DoTwist(qvecRe[(kBTot + 1) * 4 * id + i * 4 + 3], qvecIm[(kBTot + 1) * 4 * id + i * 4 + 3],
                           objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 3, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 4, i + 1));
          helperEP.DoRescale(qvecRe[(kBTot + 1) * 4 * id + i * 4 + 3], qvecIm[(kBTot + 1) * 4 * id + i * 4 + 3],
                             objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 5, i + 1), objQvec.at(id)->GetBinContent(static_cast<int>(cent) + 1, 6, i + 1));
        }
      }
      int CorrLevel = cfgCorrLevel == 0 ? 0 : cfgCorrLevel - 1;
      qvecReFT0C.push_back(qvecRe[(kBTot + 1) * 4 * id + kFT0C * 4 + CorrLevel]);
      qvecImFT0C.push_back(qvecIm[(kBTot + 1) * 4 * id + kFT0C * 4 + CorrLevel]);
      qvecReFT0A.push_back(qvecRe[(kBTot + 1) * 4 * id + kFT0A * 4 + CorrLevel]);
      qvecImFT0A.push_back(qvecIm[(kBTot + 1) * 4 * id + kFT0A * 4 + CorrLevel]);
      qvecReFT0M.push_back(qvecRe[(kBTot + 1) * 4 * id + kFT0M * 4 + CorrLevel]);
      qvecImFT0M.push_back(qvecIm[(kBTot + 1) * 4 * id + kFT0M * 4 + CorrLevel]);
      qvecReFV0A.push_back(qvecRe[(kBTot + 1) * 4 * id + kFV0A * 4 + CorrLevel]);
      qvecImFV0A.push_back(qvecIm[(kBTot + 1) * 4 * id + kFV0A * 4 + CorrLevel]);
      qvecReBPos.push_back(qvecRe[(kBTot + 1) * 4 * id + kBPos * 4 + CorrLevel]);
      qvecImBPos.push_back(qvecIm[(kBTot + 1) * 4 * id + kBPos * 4 + CorrLevel]);
      qvecReBNeg.push_back(qvecRe[(kBTot + 1) * 4 * id + kBNeg * 4 + CorrLevel]);
      qvecImBNeg.push_back(qvecIm[(kBTot + 1) * 4 * id + kBNeg * 4 + CorrLevel]);
      qvecReBTot.push_back(qvecRe[(kBTot + 1) * 4 * id + kBTot * 4 + CorrLevel]);
      qvecImBTot.push_back(qvecIm[(kBTot + 1) * 4 * id + kBTot * 4 + CorrLevel]);
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
    if (useDetector["QvectorBPoss"])
      qVectorBPos(IsCalibrated, qvecReBPos.at(0), qvecImBPos.at(0), qvecAmp[kBPos], TrkBPosLabel);
    if (useDetector["QvectorBNegs"])
      qVectorBNeg(IsCalibrated, qvecReBNeg.at(0), qvecImBNeg.at(0), qvecAmp[kBNeg], TrkBNegLabel);
    if (useDetector["QvectorBTots"])
      qVectorBTot(IsCalibrated, qvecReBTot.at(0), qvecImBTot.at(0), qvecAmp[kBTot], TrkBTotLabel);

    qVectorFT0CVec(IsCalibrated, qvecReFT0C, qvecImFT0C, qvecAmp[kFT0C]);
    qVectorFT0AVec(IsCalibrated, qvecReFT0A, qvecImFT0A, qvecAmp[kFT0A]);
    qVectorFT0MVec(IsCalibrated, qvecReFT0M, qvecImFT0M, qvecAmp[kFT0M]);
    qVectorFV0AVec(IsCalibrated, qvecReFV0A, qvecImFV0A, qvecAmp[kFV0A]);
    qVectorBPosVec(IsCalibrated, qvecReBPos, qvecImBPos, qvecAmp[kBPos], TrkBPosLabel);
    qVectorBNegVec(IsCalibrated, qvecReBNeg, qvecImBNeg, qvecAmp[kBNeg], TrkBNegLabel);
    qVectorBTotVec(IsCalibrated, qvecReBTot, qvecImBTot, qvecAmp[kBTot], TrkBTotLabel);
  } // End process.
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsTable>(cfgc)};
}
