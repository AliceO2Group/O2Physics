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
/// \file flowDirectedFlowTask.cxx
/// \brief Task for q1-dependent directed flow and global polarization
/// \since May 2026

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <TProfile2D.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;
using namespace o2::constants::physics;

struct flowDirectedFlowTask {
  Configurable<std::string> cfgUrl{"cfgUrl", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<int64_t> nolaterthan{"nolaterthan", -1, "Latest acceptable CCDB timestamp; -1 = current"};

  Configurable<bool> cfgPVSel{"cfgPVSel", true, "flag for the primary vertex selection"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<bool> cfgAddEvtSelPileup{"cfgAddEvtSelPileup", true, "flag for the pilup rejection"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  Configurable<int> cfgq1mode{"cfgq1mode", 0, "0: |QC-QA|, 1: |QA|, 2: |QC|, 3: 0.5*(|QA|+|QC|)"};

  Configurable<bool> cfgUseAccCorr{"cfgUseAccCorr", true, "Use acceptance correction profile"};
  Configurable<std::string> cfgAccPathL{"cfgAccPathL", "Users/p/prottay/My/Object/acccorrL", "CCDB path to Lambda acceptance correction"};
  Configurable<std::string> cfgAccPathAL{"cfgAccPathAL", "Users/p/prottay/My/Object/acccorrAL", "CCDB path to AntiLambda acceptance correction"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0.f, "Minimum V0 pT"};
  Configurable<float> cfgV0Rap{"cfgV0Rap", 0.8f, "V0 rapidity range"};
  Configurable<float> cfgV0DCADaughMax{"cfgV0DCADaughMax", 0.2, "Maximum DCA between V0 daughters"};
  Configurable<float> cfgV0CPAMin{"cfgV0CPAMin", 0.9998, "Minimum V0 CPA"};
  Configurable<float> cfgV0TranRadV0Min{"cfgV0TranRadV0Min", 1.5f, "Minimum V0 transverse radius"};
  Configurable<float> cfgV0TranRadV0Max{"cfgV0TranRadV0Max", 100.f, "Maximum V0 transverse radius"};
  Configurable<float> cfgMaxV0DCA{"cfgMaxV0DCA", 1.2, "Maximum V0 DCA to PV"};
  Configurable<float> cfgMinV0DCAPr{"cfgMinV0DCAPr", 0.05, "Minimum proton DCA to PV"};
  Configurable<float> cfgMinV0DCAPi{"cfgMinV0DCAPi", 0.05, "Minimum pion DCA to PV"};
  Configurable<float> cfgMaxV0LifeTime{"cfgMaxV0LifeTime", 20, "Maximum V0 lifetime"};

  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};

  Configurable<float> cfgDaughEta{"cfgDaughEta", 0.8f, "Daughter eta range"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.4f, "Minimum daughter proton pT"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2f, "Minimum daughter pion pT"};
  Configurable<float> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 50.f, "Minimum TPC found clusters"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Minimum TPC crossed rows"};
  Configurable<float> cfgRcrfc{"cfgRcrfc", 0.8f, "TPC crossed rows / findable clusters"};
  Configurable<float> cfgDaughPIDCuts{"cfgDaughPIDCuts", 3.f, "TPC PID nsigma cut"};

  Configurable<float> cfgTrackPtMin{"cfgTrackPtMin", 0.2f, "Minimum track pT for v1"};
  Configurable<float> cfgTrackPtMax{"cfgTrackPtMax", 5.0f, "Maximum track pT for v1"};
  Configurable<float> cfgTrackEtaMax{"cfgTrackEtaMax", 0.8f, "Maximum |eta| for v1"};
  Configurable<int> cfgMinTPCClustersTrack{"cfgMinTPCClustersTrack", 70, "Minimum TPC clusters for v1 tracks"};

  struct : ConfigurableGroup {
    Configurable<bool> cfgRequireRCTFlagChecker{"cfgRequireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", true, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  float massLambda = o2::constants::physics::MassLambda0;
  float massPr = o2::constants::physics::MassProton;
  float massPi = o2::constants::physics::MassPionCharged;

  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 80.}, "Centrality (%)"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.2, 0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 6.5, 8., 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis etaAxis{"etaAxis", {VARIABLE_WIDTH, -0.8, -0.4, -0.2, 0., 0.2, 0.4, 0.8}, "#eta"};
  ConfigurableAxis massAxis{"massAxis", {100, 1.0, 1.2}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis polAxis{"polAxis", {200, -1.0, 1.0}, "polarization observable"};
  ConfigurableAxis spAxis{"spAxis", {400, -10.0, 10.0}, "SP observable"};
  ConfigurableAxis q1Axis{"q1Axis", {VARIABLE_WIDTH, 0., 0.5, 1., 1.5, 2., 3., 5., 10.}, "q_{1}^{ZDC}"};
  ConfigurableAxis qAxis{"qAxis", {200, -10., 10.}, "Q"};
  ConfigurableAxis resAxis{"resAxis", {200, -1., 1.}, "resolution / correlation"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  TProfile2D* accprofileL = nullptr;
  TProfile2D* accprofileAL = nullptr;
  int currentRunNumber = -999;
  int lastRunNumber = -999;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);

  RCTFlagsChecker rctChecker;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::SPCalibrationTables, aod::Mults>>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa>>;
  using ResoV0s = aod::V0Datas;

  void init(InitContext&)
  {
    ccdb->setURL(cfgUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    histos.add("hCentrality", "Centrality", kTH1F, {{centAxis}});
    histos.add("hQ1ZDCvscent", "q1 vs centrality", kTH2F, {{centAxis}, {q1Axis}});
    histos.add("hQxAvscent", "Qx A vs centrality", kTH2F, {{centAxis}, {qAxis}});
    histos.add("hQyAvscent", "Qy A vs centrality", kTH2F, {{centAxis}, {qAxis}});
    histos.add("hQxCvscent", "Qx C vs centrality", kTH2F, {{centAxis}, {qAxis}});
    histos.add("hQyCvscent", "Qy C vs centrality", kTH2F, {{centAxis}, {qAxis}});

    histos.add("hpResCosAC", "cos(#Psi_{A}-#Psi_{C}) vs centrality", kTH3F, {{centAxis}, {resAxis}, {q1Axis}});
    histos.add("hpDotAC", "Q_{A}#upoint Q_{C} vs centrality", kTH3F, {{centAxis}, {qAxis}, {q1Axis}});
    histos.add("hpResDotAC", "Q_{A}#upoint Q_{C} vs centrality", kTH3F, {{centAxis}, {qAxis}, {q1Axis}});
    histos.add("hpQxAQxC", "QxA QxC", kTH2F, {{centAxis}, {resAxis}});
    histos.add("hpQyAQyC", "QyA QyC", kTH2F, {{centAxis}, {resAxis}});
    histos.add("hpQxAQyC", "QxA QyC", kTH2F, {{centAxis}, {resAxis}});
    histos.add("hpQxCQyA", "QxC QyA", kTH2F, {{centAxis}, {resAxis}});

    std::vector<AxisSpec> axesV1SP = {centAxis, ptAxis, etaAxis, q1Axis, spAxis};
    std::vector<AxisSpec> axesV1EP = {centAxis, ptAxis, etaAxis, q1Axis, polAxis};

    histos.add("hV1SPFullQ1", "charged v1 SP numerator u#upoint(QC-QA) vs q1", HistType::kTHnSparseF, axesV1SP, true);
    histos.add("hV1SPAQ1", "charged v1 SP numerator u#upointQA vs q1", HistType::kTHnSparseF, axesV1SP, true);
    histos.add("hV1SPCQ1", "charged v1 SP numerator u#upointQC vs q1", HistType::kTHnSparseF, axesV1SP, true);
    histos.add("hV1EPFullQ1", "charged v1 EP QA cos(phi-PsiZDC) vs q1", HistType::kTHnSparseF, axesV1EP, true);
    histos.add("hV1EPAQ1", "charged v1 EP QA cos(phi-PsiA) vs q1", HistType::kTHnSparseF, axesV1EP, true);
    histos.add("hV1EPCQ1", "charged v1 EP QA cos(phi-PsiC) vs q1", HistType::kTHnSparseF, axesV1EP, true);

    std::vector<AxisSpec> axesPolSPQ1 = {massAxis, ptAxis, spAxis, centAxis, q1Axis};

    histos.add("hSparseLambdaPolSPQ1", "Lambda SP polarization numerator vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaPolSPwgtQ1", "Lambda SP polarization numerator / acceptance vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaPolSPAQ1", "Lambda SP polarization numerator A vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaPolSPCQ1", "Lambda SP polarization numerator C vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaPolSPAwgtQ1", "Lambda SP polarization numerator A / acceptance vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaPolSPCwgtQ1", "Lambda SP polarization numerator C / acceptance vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);

    histos.add("hSparseAntiLambdaPolSPQ1", "AntiLambda SP polarization numerator vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaPolSPwgtQ1", "AntiLambda SP polarization numerator / acceptance vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaPolSPAQ1", "AntiLambda SP polarization numerator A vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaPolSPCQ1", "AntiLambda SP polarization numerator C vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaPolSPAwgtQ1", "AntiLambda SP polarization numerator A / acceptance vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaPolSPCwgtQ1", "AntiLambda SP polarization numerator C / acceptance vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);

    histos.add("hSparseLambdaPolEPQ1", "Lambda EP QA sin(phi*-PsiZDC) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaPolEPAQ1", "Lambda EP QA sin(phi*-PsiA) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaPolEPCQ1", "Lambda EP QA sin(phi*-PsiC) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaPolEPQ1", "AntiLambda EP QA sin(phi*-PsiZDC) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaPolEPAQ1", "AntiLambda EP QA sin(phi*-PsiA) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaPolEPCQ1", "AntiLambda EP QA sin(phi*-PsiC) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);

    histos.add("hSparseLambdaCosPsiQ1", "Lambda cos(PsiZDC) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaSinPsiQ1", "Lambda sin(PsiZDC) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaCosPsiQ1", "AntiLambda cos(PsiZDC) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaSinPsiQ1", "AntiLambda sin(PsiZDC) vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);

    histos.add("hSparseLambdaCorrSinPhiStarQ1", "Lambda sin(phi*) acceptance correction vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaCorrCosPhiStarQ1", "Lambda cos(phi*) acceptance correction vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaCorrSinThetaStarQ1", "Lambda sin(theta*) acceptance correction vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaAvgUxQ1", "Lambda <u_{x}> vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseLambdaAvgUyQ1", "Lambda <u_{y}> vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);

    histos.add("hSparseAntiLambdaCorrSinPhiStarQ1", "AntiLambda sin(phi*) acceptance correction vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaCorrCosPhiStarQ1", "AntiLambda cos(phi*) acceptance correction vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaCorrSinThetaStarQ1", "AntiLambda sin(theta*) acceptance correction vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaAvgUxQ1", "AntiLambda <u_{x}> vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
    histos.add("hSparseAntiLambdaAvgUyQ1", "AntiLambda <u_{y}> vs q1", HistType::kTHnSparseF, axesPolSPQ1, true);
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (cfgPVSel && std::abs(collision.posZ()) > cfgCutVertex) {
      return 0;
    }
    if (cfgAddEvtSelPileup && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy) {
      return 0;
    }
    if (rctCut.cfgRequireRCTFlagChecker && !rctChecker(collision)) {
      return 0;
    }

    return 1;
  }

  template <typename Collision>
  bool getZDCQ1(Collision const& collision, float& qxA, float& qyA, float& qxC, float& qyC,
                float& psiA, float& psiC, float& psiFull, float& q1)
  {
    const float rawQxA = collision.qxZDCA();
    const float rawQyA = collision.qyZDCA();
    const float rawQxC = collision.qxZDCC();
    const float rawQyC = collision.qyZDCC();

    psiA = collision.psiZDCA();
    psiC = collision.psiZDCC();

    const float magA = std::sqrt(rawQxA * rawQxA + rawQyA * rawQyA);
    const float magC = std::sqrt(rawQxC * rawQxC + rawQyC * rawQyC);
    qxA = magA * std::cos(psiA);
    qyA = magA * std::sin(psiA);
    qxC = magC * std::cos(psiC);
    qyC = magC * std::sin(psiC);

    const float qxFull = qxC - qxA;
    const float qyFull = qyC - qyA;
    psiFull = std::atan2(qyFull, qxFull);

    const float q1A = std::sqrt(qxA * qxA + qyA * qyA);
    const float q1C = std::sqrt(qxC * qxC + qyC * qyC);
    const float q1Full = std::sqrt(qxFull * qxFull + qyFull * qyFull);
    const float q1Mean = 0.5 * (q1A + q1C);

    if (cfgq1mode == 1) {
      q1 = q1A;
    } else if (cfgq1mode == 2) {
      q1 = q1C;
    } else if (cfgq1mode == 3) {
      q1 = q1Mean;
    } else {
      q1 = q1Full;
    }

    return std::isfinite(qxA) && std::isfinite(qyA) && std::isfinite(qxC) && std::isfinite(qyC) && std::isfinite(psiA) && std::isfinite(psiC) && std::isfinite(psiFull) && std::isfinite(q1);
  }

  template <typename Track>
  bool selectV1Track(Track const& track)
  {
    if (!track.isGlobalTrack()) {
      return false;
    }
    if (track.pt() < cfgTrackPtMin || track.pt() > cfgTrackPtMax) {
      return false;
    }
    if (std::abs(track.eta()) > cfgTrackEtaMax) {
      return false;
    }
    if (track.tpcNClsFound() < cfgMinTPCClustersTrack) {
      return false;
    }
    return true;
  }

  template <typename V0>
  bool isCompatible(V0 const& v0, int pid)
  {
    // pid = 0: Lambda, pid = 1: AntiLambda
    auto posTrack = v0.template posTrack_as<AllTrackCandidates>();
    auto negTrack = v0.template negTrack_as<AllTrackCandidates>();

    if (pid == 0 && (v0.positivept() < cfgDaughPrPt || v0.negativept() < cfgDaughPiPt)) {
      return false;
    }
    if (pid == 1 && (v0.positivept() < cfgDaughPiPt || v0.negativept() < cfgDaughPrPt)) {
      return false;
    }
    if (std::abs(v0.positiveeta()) > cfgDaughEta || std::abs(v0.negativeeta()) > cfgDaughEta) {
      return false;
    }

    if (posTrack.tpcNClsCrossedRows() < cfgTPCcluster || negTrack.tpcNClsCrossedRows() < cfgTPCcluster) {
      return false;
    }
    if (posTrack.tpcNClsFound() < cfgDaughTPCnclsMin || negTrack.tpcNClsFound() < cfgDaughTPCnclsMin) {
      return false;
    }
    if (posTrack.tpcCrossedRowsOverFindableCls() < cfgRcrfc || negTrack.tpcCrossedRowsOverFindableCls() < cfgRcrfc) {
      return false;
    }

    if (pid == 0 && (std::abs(posTrack.tpcNSigmaPr()) > cfgDaughPIDCuts || std::abs(negTrack.tpcNSigmaPi()) > cfgDaughPIDCuts)) {
      return false;
    }
    if (pid == 1 && (std::abs(posTrack.tpcNSigmaPi()) > cfgDaughPIDCuts || std::abs(negTrack.tpcNSigmaPr()) > cfgDaughPIDCuts)) {
      return false;
    }

    if (pid == 0 && (std::abs(v0.dcapostopv()) < cfgMinV0DCAPr || std::abs(v0.dcanegtopv()) < cfgMinV0DCAPi)) {
      return false;
    }
    if (pid == 1 && (std::abs(v0.dcapostopv()) < cfgMinV0DCAPi || std::abs(v0.dcanegtopv()) < cfgMinV0DCAPr)) {
      return false;
    }
    return true;
  }

  template <typename Collision, typename V0>
  bool selectV0(Collision const& collision, V0 const& v0)
  {
    const float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda;
    if (std::abs(v0.dcav0topv()) > cfgMaxV0DCA) {
      return false;
    }
    if (v0.pt() < cfgV0PtMin) {
      return false;
    }
    if (std::abs(v0.dcaV0daughters()) > cfgV0DCADaughMax) {
      return false;
    }
    if (v0.v0cosPA() < cfgV0CPAMin) {
      return false;
    }
    if (v0.v0radius() < cfgV0TranRadV0Min || v0.v0radius() > cfgV0TranRadV0Max) {
      return false;
    }
    if (std::abs(ctauLambda) > cfgMaxV0LifeTime) {
      return false;
    }
    if (std::abs(v0.yLambda()) > cfgV0Rap) {
      return false;
    }
    return true;
  }

  float acceptanceWeight(bool isLambda, float eta, float pt)
  {
    if (!cfgUseAccCorr) {
      return 1.0;
    }
    TProfile2D* prof = isLambda ? accprofileL : accprofileAL;
    if (!prof) {
      return 1.0;
    }
    const int binx = prof->GetXaxis()->FindBin(eta);
    const int biny = prof->GetYaxis()->FindBin(pt);
    const float acc = prof->GetBinContent(binx, biny);
    if (acc <= 0.0 || !std::isfinite(acc)) {
      return 1.0;
    }
    return acc;
  }

  template <typename V0>
  void fillPolarization(bool isLambda, V0 const& v0, float centrality, float q1, float psiA, float psiC, float psiFull, float qxA, float qyA, float qxC, float qyC)
  {
    ROOT::Math::PxPyPzMVector daughter;
    ROOT::Math::PxPyPzMVector pion;
    float mass = 0.0;

    if (isLambda) {
      daughter = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
      pion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
      mass = v0.mLambda();
    } else {
      daughter = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
      pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
      mass = v0.mAntiLambda();
    }

    const auto hyperon = daughter + pion;
    ROOT::Math::Boost boostToRest(hyperon.BoostToCM());
    const auto daughterStar = boostToRest(daughter);

    float phiStar = std::atan2(daughterStar.Py(), daughterStar.Px());
    const float ux = std::cos(phiStar);
    const float uy = std::sin(phiStar);

    const float cosThetaStar = daughterStar.Pz() / daughterStar.P();
    const float sinThetaStar = std::sqrt(std::max(0.0, 1.0 - cosThetaStar * cosThetaStar));
    const float sinPhiStar = std::sin(phiStar);
    const float cosPhiStar = std::cos(phiStar);

    float polEP_A = std::sin(phiStar - psiA);
    float polEP_C = std::sin(phiStar - psiC);
    float polEP = std::sin(phiStar - psiFull);

    const float qxFull = qxC - qxA;
    const float qyFull = qyC - qyA;

    float polSP = uy * qxFull - ux * qyFull;
    float polSP_A = uy * qxA - ux * qyA;
    float polSP_C = uy * qxC - ux * qyC;

    float accDen = acceptanceWeight(isLambda, v0.eta(), v0.pt());
    float wgt = 1.0; // efficiency

    polSP /= accDen;
    polSP_A /= accDen;
    polSP_C /= accDen;

    polEP /= accDen;
    polEP_A /= accDen;
    polEP_C /= accDen;

    const float cosPsi = std::cos(psiFull);
    const float sinPsi = std::sin(psiFull);

    if (isLambda) {
      histos.fill(HIST("hSparseLambdaPolSPQ1"), mass, v0.pt(), polSP, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaPolSPAQ1"), mass, v0.pt(), polSP_A, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaPolSPCQ1"), mass, v0.pt(), polSP_C, centrality, q1, wgt);

      histos.fill(HIST("hSparseLambdaPolEPQ1"), mass, v0.pt(), polEP, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaPolEPAQ1"), mass, v0.pt(), polEP_A, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaPolEPCQ1"), mass, v0.pt(), polEP_C, centrality, q1, wgt);

      histos.fill(HIST("hSparseLambdaCosPsiQ1"), mass, v0.pt(), cosPsi, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaSinPsiQ1"), mass, v0.pt(), sinPsi, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaCorrSinPhiStarQ1"), mass, v0.pt(), sinPhiStar, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaCorrCosPhiStarQ1"), mass, v0.pt(), cosPhiStar, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaCorrSinThetaStarQ1"), mass, v0.pt(), sinThetaStar, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaAvgUxQ1"), mass, v0.pt(), ux, centrality, q1, wgt);
      histos.fill(HIST("hSparseLambdaAvgUyQ1"), mass, v0.pt(), uy, centrality, q1, wgt);
    } else {
      histos.fill(HIST("hSparseAntiLambdaPolSPQ1"), mass, v0.pt(), polSP, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaPolSPAQ1"), mass, v0.pt(), polSP_A, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaPolSPCQ1"), mass, v0.pt(), polSP_C, centrality, q1, wgt);

      histos.fill(HIST("hSparseAntiLambdaPolEPQ1"), mass, v0.pt(), polEP, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaPolEPAQ1"), mass, v0.pt(), polEP_A, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaPolEPCQ1"), mass, v0.pt(), polEP_C, centrality, q1, wgt);

      histos.fill(HIST("hSparseAntiLambdaCosPsiQ1"), mass, v0.pt(), cosPsi, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaSinPsiQ1"), mass, v0.pt(), sinPsi, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaCorrSinPhiStarQ1"), mass, v0.pt(), sinPhiStar, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaCorrCosPhiStarQ1"), mass, v0.pt(), cosPhiStar, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaCorrSinThetaStarQ1"), mass, v0.pt(), sinThetaStar, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaAvgUxQ1"), mass, v0.pt(), ux, centrality, q1, wgt);
      histos.fill(HIST("hSparseAntiLambdaAvgUyQ1"), mass, v0.pt(), uy, centrality, q1, wgt);
    }
  }

  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const& tracks, ResoV0s const& v0s, aod::BCsWithTimestamps const&)
  {
    if (!eventSelected(collision)) {
      return;
    }
    const float centrality = collision.centFT0C();

    if (cfgUseAccCorr) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      currentRunNumber = bc.runNumber();
      if (currentRunNumber != lastRunNumber) {
        accprofileL = ccdb->getForTimeStamp<TProfile2D>(cfgAccPathL.value, bc.timestamp());
        accprofileAL = ccdb->getForTimeStamp<TProfile2D>(cfgAccPathAL.value, bc.timestamp());
        lastRunNumber = currentRunNumber;
      }
    }

    float qxA = 0.0, qyA = 0.0, qxC = 0.0, qyC = 0.0;
    float psiA = 0.0, psiC = 0.0, psiFull = 0.0, q1 = 0.0;
    if (!getZDCQ1(collision, qxA, qyA, qxC, qyC, psiA, psiC, psiFull, q1)) {
      return;
    }

    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hQ1ZDCvscent"), centrality, q1);
    histos.fill(HIST("hQxAvscent"), centrality, qxA);
    histos.fill(HIST("hQyAvscent"), centrality, qyA);
    histos.fill(HIST("hQxCvscent"), centrality, qxC);
    histos.fill(HIST("hQyCvscent"), centrality, qyC);

    float magA = std::sqrt(qxA * qxA + qyA * qyA);
    float magC = std::sqrt(qxC * qxC + qyC * qyC);
    float qxFull = qxC - qxA;
    float qyFull = qyC - qyA;
    float dotAC = qxA * qxC + qyA * qyC;
    float resDot = dotAC / (magA * magC);

    histos.fill(HIST("hpResCosAC"), centrality, std::cos(psiA - psiC), q1);
    histos.fill(HIST("hpDotAC"), centrality, dotAC, q1);
    histos.fill(HIST("hpResDotAC"), centrality, resDot, q1);
    histos.fill(HIST("hpQxAQxC"), centrality, qxA * qxC);
    histos.fill(HIST("hpQyAQyC"), centrality, qyA * qyC);
    histos.fill(HIST("hpQxAQyC"), centrality, qxA * qyC);
    histos.fill(HIST("hpQxCQyA"), centrality, qxC * qyA);

    for (const auto& track : tracks) {
      if (!selectV1Track(track)) {
        continue;
      }

      float phi = track.phi();
      float ux = std::cos(phi);
      float uy = std::sin(phi);

      float v1SPFull = ux * qxFull + uy * qyFull;
      float v1SPA = ux * qxA + uy * qyA;
      float v1SPC = ux * qxC + uy * qyC;

      histos.fill(HIST("hV1SPFullQ1"), centrality, track.pt(), track.eta(), q1, v1SPFull);
      histos.fill(HIST("hV1SPAQ1"), centrality, track.pt(), track.eta(), q1, v1SPA);
      histos.fill(HIST("hV1SPCQ1"), centrality, track.pt(), track.eta(), q1, v1SPC);

      histos.fill(HIST("hV1EPFullQ1"), centrality, track.pt(), track.eta(), q1, std::cos(phi - psiFull));
      histos.fill(HIST("hV1EPAQ1"), centrality, track.pt(), track.eta(), q1, std::cos(phi - psiA));
      histos.fill(HIST("hV1EPCQ1"), centrality, track.pt(), track.eta(), q1, std::cos(phi - psiC));
    }

    for (const auto& v0 : v0s) {
      bool lambdaTag = isCompatible(v0, 0);
      bool antiLambdaTag = isCompatible(v0, 1);
      if (!lambdaTag && !antiLambdaTag) {
        continue;
      }
      if (!selectV0(collision, v0)) {
        continue;
      }
      if (lambdaTag) {
        fillPolarization(true, v0, centrality, q1, psiA, psiC, psiFull, qxA, qyA, qxC, qyC);
      }
      if (antiLambdaTag) {
        fillPolarization(false, v0, centrality, q1, psiA, psiC, psiFull, qxA, qyA, qxC, qyC);
      }
    }
  }
  PROCESS_SWITCH(flowDirectedFlowTask, processData, "Process for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<flowDirectedFlowTask>(cfgc)};
}
