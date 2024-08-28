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

// q vector framework with ESE (20/08/2024)
//
/// \author Joachim Hansen <joachim.hansen@cern.ch>
//

#include <chrono>
#include <string>
#include <TComplex.h>
#include <algorithm>
#include <numeric>
#include <vector>

#include "Framework/ASoA.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DataFormatsParameters/GRPMagField.h>

#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/FT0Corrected.h"

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FT0Base/Geometry.h"

#include "PWGCF/Flow/DataModel/FlowESE.h"

#include "FFitWeights.h"
#include "TSpline.h"

using namespace o2;
using namespace o2::framework;
// using namespace o2::framework::expressions;

// using CollWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;
using CollWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

struct ProduceFlowESE {
  Produces<o2::aod::qVecFV0As> qVectorFV0A;
  Produces<o2::aod::qVecFT0Cs> qVectorFT0C;
  Produces<o2::aod::qPercentileFT0Cs> qPercs;

  OutputObj<FFitWeights> FFitObj{FFitWeights("weights")};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<int> cfgCorrLevel{"cfgCorrLevel", 3, "calibration step: 0 = no corr, 1 = gain corr, 2 = rectr, 4 = full"};
  Configurable<bool> cfgESE{"cfgESE", 1, "ese actovation step: false = no ese, true = evaluate splines and fill table"};

  Configurable<std::string> cfgCalibrationPath{"cfgCalibrationPath", "Users/j/joachiha/Calibration", "CCDB path for gain equalization constants"};
  Configurable<std::string> cfgEsePath{"cfgEsePath", "Users/j/joachiha/ESE/local/splines", "CCDB path for ese splines"};

  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, "Fit amplitude range"};
  ConfigurableAxis cfgaxisq2FT0C{"cfgaxisq2FT0C", {200, 0, 25}, "q2 amplitude range"};
  ConfigurableAxis cfgaxisQFT0C{"cfgaxisQFT0C", {250, -3500, 3500}, "q2 amplitude range"};

  AxisSpec axisFITamp{cfgaxisFITamp, "FIT amp"};
  AxisSpec axisq2FT0C{cfgaxisq2FT0C, "q2 FIT amp"};
  AxisSpec axisQFT0C{cfgaxisQFT0C, "Q xy Range"};

  AxisSpec axisCentralityR = {100, 0, 100};
  AxisSpec axisChID = {220, 0, 220};
  std::vector<float> FV0RelGainConst{};
  std::vector<float> FT0RelGainConst{};
  int runNumber{-1};

  o2::ft0::Geometry ft0geom;
  double mOffsetFT0AX = 0.; // X-coordinate of the offset of FT0-A.
  double mOffsetFT0AY = 0.; // Y-coordinate of the offset of FT0-A.
  double mOffsetFT0CX = 0.; // X-coordinate of the offset of FT0-C.
  double mOffsetFT0CY = 0.; // Y-coordinate of the offset of FT0-C.

  FFitWeights* weights{nullptr};
  TList* eselist{nullptr};
  std::vector<std::vector<TSpline3*>> spl{2, std::vector<TSpline3*>(90)}; // 90x2

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  void init(o2::framework::InitContext&)
  {

    LOGF(info, "ESETable::init()");

    registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});

    registry.add("h_qx2VecFT0C", "qxVecFT0C;qxVecFT0C;entries", {HistType::kTH1F, {{1000, 0, 1000}}});
    registry.add("h_qy2VecFT0C", "qyVecFT0C;qyVecFT0C;entries", {HistType::kTH1F, {{1000, 0, 1000}}});
    registry.add("h_qx3VecFT0C", "qxVecFT0C;qxVecFT0C;entries", {HistType::kTH1F, {{1000, 0, 1000}}});
    registry.add("h_qy3VecFT0C", "qyVecFT0C;qyVecFT0C;entries", {HistType::kTH1F, {{1000, 0, 1000}}});

    registry.add("h_Cent_q2FT0C", "q^{FT0C}_{2};Centrality (FT0A);", {HistType::kTH2F, {axisCentralityR, axisq2FT0C}});
    registry.add("h_Cent_q3FT0C", "q^{FT0C}_{3};Centrality (FT0A);", {HistType::kTH2F, {axisCentralityR, axisq2FT0C}});

    registry.add("h_Psi2", "#Psi_{2}^{FT0C};Centrality;", {HistType::kTH2F, {axisCentralityR, {100, -5, 5}}});
    registry.add("h_Psi3", "#Psi_{3}^{FT0C};Centrality;", {HistType::kTH2F, {axisCentralityR, {100, -5, 5}}});

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    FFitObj->Init();
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    FV0RelGainConst.clear();
    FT0RelGainConst.clear();
    FV0RelGainConst = {};
    FT0RelGainConst = {};

    auto timestamp = bc.timestamp();

    // auto offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    // auto offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", timestamp);

    // if (!objfv0Gain || cfgCorrLevel == 0) {
    //   for (auto i{0u}; i < 48; i++) {
    //     FV0RelGainConst.push_back(1.);
    //   }
    // } else {
    //   FV0RelGainConst = *(objfv0Gain);
    // }
    for (auto i{0u}; i < 48; i++) {
      FV0RelGainConst.push_back(1.);
    }

    std::string fullPath = cfgCalibrationPath;

    weights = ccdb->getForTimeStamp<FFitWeights>(fullPath, timestamp);
    if (!weights || cfgCorrLevel == 0) {
      for (auto i{0u}; i < 208; i++) {
        FT0RelGainConst.push_back(1.);
      }
    } else {
      weights->CreateGain();
      FT0RelGainConst = weights->GetGain();
    }
    if (cfgCorrLevel > 1) {
      weights->CreateRecenter("x");
      weights->CreateRMS("x");
      weights->CreateRecenter("y");
      weights->CreateRMS("y");
    }

    if (cfgESE) {
      eselist = ccdb->getForTimeStamp<TList>(cfgEsePath, timestamp);
      if (!LoadSplines())
        LOGF(fatal, "failed loading splines with ese flag");
      LOGF(info, "successfully loaded splines");
    }
  }

  // void processQVecFV0 (CollWithFV0 const& collision, aod::BCsWithTimestamps const&) {
  // template <typename CollType>
  // void QVecFV0 (CollType const& collision) {
  //   TComplex Qvec(0);
  //   float qVecFV0A[2] = {0.};
  //   float sumAmplFV0A{0.0f};

  //   if (collision.has_foundFV0()){
  //     auto fv0 = collision.foundFV0();

  //     for (std::size_t iCh = 0; iCh < fv0.channel().size(); iCh++) {
  //       float ampl = fv0.amplitude()[iCh];
  //       int FV0AchId = fv0.channel()[iCh];
  //       registry.fill(HIST("FV0Amplitude"), ampl, FV0AchId);
  //       registry.fill(HIST("FV0AmplitudeCorr"), ampl / FV0RelGainConst[FV0AchId], FV0AchId);

  //       // helperEP.SumQvectors(1, FV0AchId, ampl / FV0RelGainConst[FV0AchId], nmode, QvecDet, sumAmplFV0A, ft0geom, fv0geom);
  //     }

  //     if (sumAmplFV0A > 1e-8) {
  //       Qvec /= sumAmplFV0A;
  //       qVecFV0A[0] = Qvec.Re();
  //       qVecFV0A[1] = Qvec.Im();
  //     } else {
  //       qVecFV0A[0] = 999.;
  //       qVecFV0A[1] = 999.;
  //     }
  //   }
  //   else {
  //     qVecFV0A[0] = -999.;
  //     qVecFV0A[1] = -999.;
  //   }

  //   //// qVector(Qvec.Re(),Qvec.Im());
  //   registry.fill(HIST("h_collisions"), 1.5);
  //   qVectorFV0A(1,-2);
  // }
  bool LoadSplines()
  {
    for (int i{0}; i < 90; i++) {
      for (int j{0}; j < 2; j++) {
        spl[j][i] = static_cast<TSpline3*>(eselist->FindObject(Form("sp_q%iFT0C_%i", j + 2, i)));
      }
    }
    return true;
  }

  float Calcqn(const float& Qx, const float& Qy, const float& Mult)
  {
    float dqn{0.0f};
    float qn{0.0f};

    dqn = Qx * Qx + Qy * Qy;
    qn = TMath::Sqrt(dqn) / TMath::Sqrt(Mult);
    return qn;
  }

  double GetPhiFT0(const int& chno, o2::ft0::Geometry ft0geom)
  {
    /* Calculate the azimuthal angle in FT0 for the channel number 'chno'. The offset
      of FT0-A is taken into account if chno is between 0 and 95. */

    float offsetX = 0.;
    float offsetY = 0.; // No offset for FT0-C (default case).

    if (chno < 96) { // Channel in FT0-A, non-zero offset must be applied. // LOKI: make general.
      offsetX = mOffsetFT0AX;
      offsetY = mOffsetFT0AY;
    }

    ft0geom.calculateChannelCenter();
    auto chPos = ft0geom.getChannelCenter(chno);
    /// printf("Channel id: %d X: %.3f Y: %.3f\n", chno, chPos.X(), chPos.Y());

    return TMath::ATan2(chPos.Y() + offsetY, chPos.X() + offsetX);
  }

  float EventPlane(const float& x, const float& y, const float& nHarm)
  {
    return 1 / nHarm * TMath::ATan2(y, x);
  }

  void SumQvectors(const int& chno, const float& ampl, const int& nHarm, TComplex& Qvec, float& sum, o2::ft0::Geometry ft0geom)
  {
    /* Calculate the complex Q-vector for the provided detector and channel number,
      before adding it to the total Q-vector given as argument. */
    double phi = -999.;

    phi = GetPhiFT0(chno, ft0geom); // already be given in the right range.

    if (phi < -900) {
      printf("Error on phi. Skip\n");
      return;
    }
    Qvec += TComplex(ampl * TMath::Cos(phi * nHarm), ampl * TMath::Sin(phi * nHarm));
    sum += ampl;
  }

  template <typename CollType>
  void QVecFT0C(CollType const& collision, const int& nHarm, std::vector<float>& qx, std::vector<float>& qy, std::vector<float>& qnP)
  {
    TComplex Qvec(0);
    float qVecFT0C[2] = {0.};
    float sumAmplFT0C{0.0f};
    float qn{0.0f};
    bool fCalc{0};

    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();
      // auto ft0 = collision.foundFT0();
      for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
        float ampl = ft0.amplitudeC()[iChC];
        int FT0CchId = ft0.channelC()[iChC] + 96;

        // registry.fill(HIST("FT0Amplitude"), FT0CchId, ampl);
        // registry.fill(HIST("FT0AmplitudeCorr"), FT0CchId, ampl / FT0RelGainConst[FT0CchId]);
        if (nHarm == 2)
          FFitObj->FillFT0(FT0CchId, ampl, FT0RelGainConst[FT0CchId]);

        SumQvectors(FT0CchId, ampl / FT0RelGainConst[FT0CchId], nHarm, Qvec, sumAmplFT0C, ft0geom);
      }

      if (sumAmplFT0C > 1e-8) {
        // Qvec /= sumAmplFT0C;
        qVecFT0C[0] = Qvec.Re();
        qVecFT0C[1] = Qvec.Im();
        fCalc = true;
      } else {
        qVecFT0C[0] = 999.;
        qVecFT0C[1] = 999.;
        fCalc = false;
      }
    } else {
      qVecFT0C[0] = 999.;
      qVecFT0C[1] = 999.;
    }

    if (!fCalc) {
      qn = 0;
    } else {
      qn = Calcqn(qVecFT0C[0], qVecFT0C[1], sumAmplFT0C);
    }

    FFitObj->FillQ(collision.centFT0C(), qVecFT0C[0], nHarm, "x", "");
    FFitObj->FillQ(collision.centFT0C(), qVecFT0C[1], nHarm, "y", "");

    if (cfgCorrLevel > 1) {
      int centr = static_cast<int>(collision.centFT0C());
      if (fCalc) {
        qVecFT0C[0] = qVecFT0C[0] - weights->GetRecVal(centr, "x", nHarm);
        qVecFT0C[1] = qVecFT0C[1] - weights->GetRecVal(centr, "y", nHarm);
        FFitObj->FillQ(collision.centFT0C(), qVecFT0C[0], nHarm, "x", "_Rec");
        FFitObj->FillQ(collision.centFT0C(), qVecFT0C[1], nHarm, "y", "_Rec");

        qVecFT0C[0] = qVecFT0C[0] / weights->GetRMSVal(centr, "x", nHarm);
        qVecFT0C[1] = qVecFT0C[1] / weights->GetRMSVal(centr, "y", nHarm);
        FFitObj->FillQ(collision.centFT0C(), qVecFT0C[0], nHarm, "x", "_RecTot");
        FFitObj->FillQ(collision.centFT0C(), qVecFT0C[1], nHarm, "y", "_RecTot");
      }
    }

    float Psi = EventPlane(qVecFT0C[0], qVecFT0C[1], nHarm);

    if (nHarm == 2) {
      registry.fill(HIST("h_qx2VecFT0C"), qVecFT0C[0]);
      registry.fill(HIST("h_qy2VecFT0C"), qVecFT0C[1]);
      if (fCalc) {
        registry.fill(HIST("h_Cent_q2FT0C"), collision.centFT0C(), qn);
        registry.fill(HIST("h_Psi2"), collision.centFT0C(), Psi);
      }
    } else if (nHarm == 3) {
      registry.fill(HIST("h_qx3VecFT0C"), qVecFT0C[0]);
      registry.fill(HIST("h_qy3VecFT0C"), qVecFT0C[1]);
      if (fCalc) {
        registry.fill(HIST("h_Cent_q3FT0C"), collision.centFT0C(), qn);
        registry.fill(HIST("h_Psi3"), collision.centFT0C(), Psi);
      }
    }

    if (cfgESE) {
      int qSpCent = static_cast<int>(collision.centFT0C());
      float qnCent{-1};
      if (qSpCent > 0 && qSpCent < 90)
        qnCent = 100. * spl[nHarm - 2][qSpCent]->Eval(qn);

      qnP.push_back(qnCent);
    } else {
      qnP.push_back(-1);
    }

    qx.push_back(qVecFT0C[0]);
    qy.push_back(qVecFT0C[1]);
    registry.fill(HIST("h_collisions"), 1.5);
  }

  void processQVecs(CollWithMults::iterator const& collision, aod::BCsWithTimestamps const&, aod::FV0As const&, aod::FT0s const&)
  {

    std::vector<float> qvecRe{};
    std::vector<float> qvecIm{};
    std::vector<float> qnp{};

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRun = bc.runNumber();
    if (runNumber != currentRun) {
      runNumber = currentRun;
      initCCDB(bc);
    }
    registry.fill(HIST("h_collisions"), 0.5);

    // QVecFV0(collision);
    QVecFT0C(collision, 2, qvecRe, qvecIm, qnp);
    QVecFT0C(collision, 3, qvecRe, qvecIm, qnp);

    qVectorFT0C(qvecRe, qvecIm);
    qPercs(qnp);
  }
  PROCESS_SWITCH(ProduceFlowESE, processQVecs, "procc q vectors ", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<ProduceFlowESE>(cfgc, TaskName{"flow-qv-ese"})}; }
