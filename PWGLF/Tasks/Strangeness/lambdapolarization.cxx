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

/// \author Junlee Kim (jikim1290@gmail.com)

#include <cmath>
#include <array>
#include <cstdlib>
#include <chrono>
#include <string>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TVector2.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"

#include "CommonConstants/PhysicsConstants.h"

#include "ReconstructionDataFormats/Track.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct lambdapolarization {
  //  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::Qvectors>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using V0TrackCandidate = aod::V0Datas;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL",
                                     "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than",
                                      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
                                      "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<float> cfgCentSel{"cfgCentSel", 80., "Centrality selection"};
  Configurable<int> cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPosToPVMin{"cfgDCAPosToPVMin", 0.05, "minimum DCA to PV for positive track"};
  Configurable<float> cfgDCANegToPVMin{"cfgDCANegToPVMin", 0.2, "minimum DCA to PV for negative track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<bool> cfgQAv0{"cfgQAv0", false, "QA plot"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 70, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 5, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 5, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.5, "minimum daughter pion pt"};

  Configurable<int> cfgnMods{"cfgnMods", 1, "The number of modulations of interest starting from 2"};
  Configurable<int> cfgNQvec{"cfgNQvec", 7, "The number of total Qvectors for looping over the task"};

  Configurable<std::string> cfgQvecDetName{"cfgQvecDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgQvecRefAName{"cfgQvecRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgQvecRefBName{"cfgQvecRefBName", "TPCneg", "The name of detector for reference B"};

  Configurable<bool> cfgShiftCorr{"cfgShiftCorr", false, "additional shift correction"};
  Configurable<bool> cfgShiftCorrDef{"cfgShiftCorrDef", false, "additional shift correction definition"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "Users/j/junlee/Qvector/QvecCalib/Shift", "Path for Shift"};

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  int DetId;
  int RefAId;
  int RefBId;

  int QvecDetInd;
  int QvecRefAInd;
  int QvecRefBInd;

  float centrality;

  double angle;
  double relphi;

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  std::vector<TProfile3D*> shiftprofile{};

  std::string fullCCDBShiftCorrPath;

  template <typename T>
  int GetDetId(const T& name)
  {
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos") {
      return 4;
    } else if (name.value == "TPCneg") {
      return 5;
    } else {
      return 0;
    }
  }

  void init(o2::framework::InitContext&)
  {
    AxisSpec massAxis = {100, 1.065, 1.165};
    AxisSpec ptAxis = {100, 0.0, 10.0};
    AxisSpec cosAxis = {110, -1.05, 1.05};
    AxisSpec centAxis = {8, 0.0, 80.0};
    AxisSpec centQaAxis = {80, 0.0, 80.0};
    AxisSpec epAxis = {6, 0.0, 2.0 * constants::math::PI};
    AxisSpec epQaAxis = {100, -1.0 * constants::math::PI, constants::math::PI};

    AxisSpec pidAxis = {100, -10, 10};

    AxisSpec shiftAxis = {10, 0, 10, "shift"};
    AxisSpec basisAxis = {20, 0, 20, "basis"};

    for (auto i = 2; i < cfgnMods + 2; i++) {
      histos.add(Form("psi%d/h_lambda_cos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_alambda_cos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_lambda_cos2", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_alambda_cos2", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});

      histos.add(Form("psi%d/h_lambda_cossin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_alambda_cossin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
    }

    if (cfgQAv0) {
      histos.add("QA/CentDist", "", {HistType::kTH1F, {centQaAxis}});

      histos.add("QA/nsigma_tpc_pt_ppr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_ppi", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpi", "", {HistType::kTH2F, {ptAxis, pidAxis}});

      for (auto i = 2; i < cfgnMods + 2; i++) {
        histos.add(Form("psi%d/QA/EP_Det", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_RefA", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_RefB", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});

        histos.add(Form("psi%d/QA/EPRes_Det_RefA", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_Det_RefB", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_RefA_RefB", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});

        histos.add(Form("psi%d/QA/EP_FT0C_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_FT0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});

        histos.add(Form("psi%d/QA/EPRes_FT0C_FT0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_FT0C_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_FT0A_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
      }
    }

    if (cfgShiftCorrDef) {
      for (auto i = 2; i < cfgnMods + 2; i++) {
        histos.add(Form("psi%d/ShiftFIT", i), "", kTProfile3D, {centQaAxis, basisAxis, shiftAxis});
      }
    }

    DetId = GetDetId(cfgQvecDetName);
    RefAId = GetDetId(cfgQvecRefAName);
    RefBId = GetDetId(cfgQvecRefBName);

    if (DetId == RefAId || DetId == RefBId || RefAId == RefBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      DetId = 0;
      RefAId = 4;
      RefBId = 5;
    }

    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  ROOT::Math::PxPyPzMVector ProtonVec, PionVec, LambdaVec, ProtonBoostedVec, PionBoostedVec;

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return 0;
    }

    if (cfgCentSel < centrality) {
      return 0;
    }
    /*
        auto multNTracksPV = collision.multNTracksPV();
        if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
          return 0;
        }
        if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
          return 0;
        }
    */
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }

    return 1;
  } // event selection

  template <typename TCollision, typename V0>
  bool SelectionV0(TCollision const& collision, V0 const& candidate)
  {
    if (candidate.v0radius() < cfgv0radiusMin)
      return false;
    if (std::abs(candidate.dcapostopv()) < cfgDCAPosToPVMin)
      return false;
    if (std::abs(candidate.dcanegtopv()) < cfgDCANegToPVMin)
      return false;
    if (candidate.v0cosPA() < cfgv0CosPA)
      return false;
    if (candidate.dcaV0daughters() > cfgDCAV0Dau)
      return false;
    if (candidate.pt() < cfgV0PtMin)
      return false;
    if (candidate.yLambda() < cfgV0EtaMin)
      return false;
    if (candidate.yLambda() > cfgV0EtaMax)
      return false;
    if (candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > cfgV0LifeTime)
      return false;

    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, int pid) // pid 0: proton, pid 1: pion
  {
    if (track.tpcNClsFound() < cfgDaughTPCnclsMin)
      return false;
    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr)
      return false;
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi)
      return false;
    if (track.eta() > cfgDaughEtaMax)
      return false;
    if (track.eta() < cfgDaughEtaMin)
      return false;
    if (pid == 0 && track.pt() < cfgDaughPrPt)
      return false;
    if (pid == 1 && track.pt() < cfgDaughPiPt)
      return false;

    return true;
  }

  template <typename TCollision>
  void FillShiftCorrection(TCollision const& collision, int nmode)
  {
    QvecDetInd = DetId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefAInd = RefAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefBInd = RefBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    for (int ishift = 1; ishift <= 10; ishift++) {
      if (nmode == 2) {
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi2/ShiftFIT"), centrality, 2.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 3.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi2/ShiftFIT"), centrality, 4.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 5.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
      } else if (nmode == 3) {
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi3/ShiftFIT"), centrality, 2.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 3.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi3/ShiftFIT"), centrality, 4.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 5.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
      } else if (nmode == 4) {
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi4/ShiftFIT"), centrality, 2.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 3.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi4/ShiftFIT"), centrality, 4.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 5.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
      }
    }
  }

  template <typename TCollision>
  void FillEPQA(TCollision const& collision, int nmode)
  {
    QvecDetInd = DetId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefAInd = RefAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefBInd = RefBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    if (collision.qvecAmp()[DetId] < 1e-5 || collision.qvecAmp()[RefAId] < 1e-5 || collision.qvecAmp()[RefBId] < 1e-5)
      return;

    if (nmode == 2) {
      histos.fill(HIST("psi2/QA/EP_Det"), centrality, TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi2/QA/EP_RefA"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi2/QA/EP_RefB"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi2/QA/EPRes_Det_RefA"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd])));
      histos.fill(HIST("psi2/QA/EPRes_Det_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
      histos.fill(HIST("psi2/QA/EPRes_RefA_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
    } else if (nmode == 3) {
      histos.fill(HIST("psi3/QA/EP_Det"), centrality, TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi3/QA/EP_RefA"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi3/QA/EP_RefB"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi3/QA/EPRes_Det_RefA"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd])));
      histos.fill(HIST("psi3/QA/EPRes_Det_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
      histos.fill(HIST("psi3/QA/EPRes_RefA_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
    } else if (nmode == 4) {
      histos.fill(HIST("psi4/QA/EP_Det"), centrality, TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi4/QA/EP_RefA"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi4/QA/EP_RefB"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi4/QA/EPRes_Det_RefA"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd])));
      histos.fill(HIST("psi4/QA/EPRes_Det_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
      histos.fill(HIST("psi4/QA/EPRes_RefA_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
    }

    if (cfgShiftCorr) {
      auto deltapsiFT0C = 0.0;
      auto deltapsiFT0A = 0.0;
      auto deltapsiFV0A = 0.0;

      auto psidefFT0C = TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode);
      auto psidefFT0A = TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode);
      auto psidefFV0A = TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode);
      for (int ishift = 1; ishift <= 10; ishift++) {
        auto coeffshiftxFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 0.5, ishift - 0.5));
        auto coeffshiftyFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 1.5, ishift - 0.5));
        auto coeffshiftxFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 2.5, ishift - 0.5));
        auto coeffshiftyFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 3.5, ishift - 0.5));
        auto coeffshiftxFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 4.5, ishift - 0.5));
        auto coeffshiftyFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 5.5, ishift - 0.5));

        deltapsiFT0C += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0C) + coeffshiftyFT0C * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0C)));
        deltapsiFT0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0A) + coeffshiftyFT0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0A)));
        deltapsiFV0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFV0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFV0A) + coeffshiftyFV0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFV0A)));
      }
      if (nmode == 2) {
        histos.fill(HIST("psi2/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi2/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi2/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi2/QA/EPRes_FT0C_FT0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi2/QA/EPRes_FT0C_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi2/QA/EPRes_FT0A_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      } else if (nmode == 3) {
        histos.fill(HIST("psi3/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi3/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi3/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi3/QA/EPRes_FT0C_FT0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi3/QA/EPRes_FT0C_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi3/QA/EPRes_FT0A_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      } else if (nmode == 4) {
        histos.fill(HIST("psi4/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi4/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi4/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi4/QA/EPRes_FT0C_FT0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi4/QA/EPRes_FT0C_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi4/QA/EPRes_FT0A_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      }
    }
  }

  template <typename TCollision, typename V0>
  void FillHistograms(TCollision const& collision, V0 const& V0s, int nmode)
  {
    QvecDetInd = DetId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefAInd = RefAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefBInd = RefBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    for (auto& v0 : V0s) {
      auto postrack = v0.template posTrack_as<TrackCandidates>();
      auto negtrack = v0.template negTrack_as<TrackCandidates>();

      double nTPCSigmaPosPr = postrack.tpcNSigmaPr();
      double nTPCSigmaNegPi = negtrack.tpcNSigmaPi();

      double nTPCSigmaNegPr = negtrack.tpcNSigmaPr();
      double nTPCSigmaPosPi = postrack.tpcNSigmaPi();

      if (cfgQAv0 && nmode == 2) {
        histos.fill(HIST("QA/nsigma_tpc_pt_ppr"), postrack.pt(), nTPCSigmaPosPr);
        histos.fill(HIST("QA/nsigma_tpc_pt_ppi"), postrack.pt(), nTPCSigmaPosPi);

        histos.fill(HIST("QA/nsigma_tpc_pt_mpr"), negtrack.pt(), nTPCSigmaNegPr);
        histos.fill(HIST("QA/nsigma_tpc_pt_mpi"), negtrack.pt(), nTPCSigmaNegPi);
      }

      int LambdaTag = 0;
      int aLambdaTag = 0;

      if (isSelectedV0Daughter(postrack, 0) && isSelectedV0Daughter(negtrack, 1)) {
        LambdaTag = 1;
      }
      if (isSelectedV0Daughter(negtrack, 0) && isSelectedV0Daughter(postrack, 1)) {
        aLambdaTag = 1;
      }

      if (LambdaTag == aLambdaTag)
        continue;

      if (!SelectionV0(collision, v0))
        continue;

      if (LambdaTag) {
        ProtonVec = ROOT::Math::PxPyPzMVector(postrack.px(), postrack.py(), postrack.pz(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(negtrack.px(), negtrack.py(), negtrack.pz(), massPi);
      }
      if (aLambdaTag) {
        ProtonVec = ROOT::Math::PxPyPzMVector(negtrack.px(), negtrack.py(), negtrack.pz(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(postrack.px(), postrack.py(), postrack.pz(), massPi);
      }
      LambdaVec = ProtonVec + PionVec;
      LambdaVec.SetM(massLambda);

      ROOT::Math::Boost boost{LambdaVec.BoostToCM()};
      ProtonBoostedVec = boost(ProtonVec);

      angle = ProtonBoostedVec.Pz() / ProtonBoostedVec.P();
      relphi = TVector2::Phi_0_2pi(static_cast<float>(nmode) * (LambdaVec.Phi() - TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));

      if (cfgShiftCorr) {
        auto deltapsiFT0C = 0.0;
        auto deltapsiFT0A = 0.0;
        auto deltapsiFV0A = 0.0;

        auto psidefFT0C = TMath::ATan2(collision.qvecIm()[3 + (nmode - 2) * 28], collision.qvecRe()[3 + (nmode - 2) * 28]) / static_cast<float>(nmode);
        auto psidefFT0A = TMath::ATan2(collision.qvecIm()[3 + 4 + (nmode - 2) * 28], collision.qvecRe()[3 + 4 + (nmode - 2) * 28]) / static_cast<float>(nmode);
        auto psidefFV0A = TMath::ATan2(collision.qvecIm()[3 + 12 + (nmode - 2) * 28], collision.qvecRe()[3 + 12 + (nmode - 2) * 28]) / static_cast<float>(nmode);
        for (int ishift = 1; ishift <= 10; ishift++) {
          auto coeffshiftxFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 1.5, ishift - 0.5));
          auto coeffshiftxFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 2.5, ishift - 0.5));
          auto coeffshiftyFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 3.5, ishift - 0.5));
          auto coeffshiftxFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 4.5, ishift - 0.5));
          auto coeffshiftyFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 5.5, ishift - 0.5));

          deltapsiFT0C += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0C) + coeffshiftyFT0C * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0C)));
          deltapsiFT0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0A) + coeffshiftyFT0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0A)));
          deltapsiFV0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFV0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFV0A) + coeffshiftyFV0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFV0A)));
        }
        relphi = TVector2::Phi_0_2pi(static_cast<float>(nmode) * (LambdaVec.Phi() - psidefFT0C - deltapsiFT0C));
      }

      if (nmode == 2) { ////////////
        if (LambdaTag) {
          histos.fill(HIST("psi2/h_lambda_cos"), v0.mLambda(), v0.pt(), angle, centrality, relphi);
          histos.fill(HIST("psi2/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi2/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * TMath::Sin(relphi), centrality);
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi2/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle, centrality, relphi);
          histos.fill(HIST("psi2/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi2/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(relphi), centrality);
        }
      } else if (nmode == 3) {
        if (LambdaTag) {
          histos.fill(HIST("psi3/h_lambda_cos"), v0.mLambda(), v0.pt(), angle, centrality, relphi);
          histos.fill(HIST("psi3/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi3/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * TMath::Sin(relphi), centrality);
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi3/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle, centrality, relphi);
          histos.fill(HIST("psi3/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi3/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(relphi), centrality);
        }
      } else if (nmode == 4) {
        if (LambdaTag) {
          histos.fill(HIST("psi4/h_lambda_cos"), v0.mLambda(), v0.pt(), angle, centrality, relphi);
          histos.fill(HIST("psi4/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi4/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * TMath::Sin(relphi), centrality);
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi4/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle, centrality, relphi);
          histos.fill(HIST("psi4/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi4/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(relphi), centrality);
        }
      } ////////// FIXME: not possible to get histograms using nmode
    }
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s,
                   aod::BCsWithTimestamps const&)
  {
    if (cfgCentEst == 1) {
      centrality = collision.centFT0C();
    } else if (cfgCentEst == 2) {
      centrality = collision.centFT0M();
    }
    if (!eventSelected(collision)) {
      return;
    }
    histos.fill(HIST("QA/CentDist"), centrality, 1.0);
    if (cfgShiftCorr) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      currentRunNumber = bc.runNumber();
      if (currentRunNumber != lastRunNumber) {
        shiftprofile.clear();
        for (int i = 2; i < cfgnMods + 2; i++) {
          fullCCDBShiftCorrPath = cfgShiftPath;
          fullCCDBShiftCorrPath += "/v";
          fullCCDBShiftCorrPath += std::to_string(i);
          auto objshift = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPath, bc.timestamp());
          shiftprofile.push_back(objshift);
        }
        lastRunNumber = currentRunNumber;
      }
    }
    for (int i = 2; i < cfgnMods + 2; i++) {
      if (cfgShiftCorrDef) {
        FillShiftCorrection(collision, i);
      }
      if (cfgQAv0) {
        FillEPQA(collision, i);
      }
      FillHistograms(collision, V0s, i);
    } // FIXME: need to fill different histograms for different harmonic
  }
  PROCESS_SWITCH(lambdapolarization, processData, "Process Event for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdapolarization>(cfgc, TaskName{"lf-lambdapolarization"})};
}
