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

/// \file   SpectatorPlaneTutorial.cxx
/// \author Noor Koster
/// \since  11/2025
/// \brief  This is a tutorial task to show how to use the ZDC q-vectors and the spectator plane resolution.

#include "PWGCF/DataModel/SPTableZDC.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include "TF1.h"
#include "TPDGCode.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;
// using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct SpectatorPlaneTutorial {
  RCTFlagsChecker rctChecker;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgEvtUseRCTFlagChecker, bool, false, "Evt sel: use RCT flag checker");
    O2_DEFINE_CONFIGURABLE(cfgEvtRCTFlagCheckerLabel, std::string, "CBT_hadronPID", "Evt sel: RCT flag checker label (CBT, CBT_hadronPID)"); // all Labels can be found in Common/CCDB/RCTSelectionFlags.h
    O2_DEFINE_CONFIGURABLE(cfgEvtRCTFlagCheckerZDCCheck, bool, false, "Evt sel: RCT flag checker ZDC check");
    O2_DEFINE_CONFIGURABLE(cfgEvtRCTFlagCheckerLimitAcceptAsBad, bool, false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad");
  } rctFlags;

  // settings for CCDB data
  O2_DEFINE_CONFIGURABLE(cfgCCDBdir_QQ, std::string, "Users/c/ckoster/ZDC/LHC23_PbPb_pass5/meanQQ/Default", "ccdb dir for average QQ values in 1% centrality bins");
  O2_DEFINE_CONFIGURABLE(cfgCCDBdir_SP, std::string, "Users/c/ckoster/ZDC/LHC23_zzh_pass4/SPPlaneRes", "ccdb dir for average event plane resolution in 1% centrality bins");
  // Confogirable axis
  ConfigurableAxis axisCentrality{"axisCentrality", {20, 0, 100}, "Centrality bins for vn "};

  // Configurables containing vector
  float vtxZ = 10.0;
  float etaMax = 0.8;
  float ptMin = 0.2;
  float ptMax = 10.0;
  float dcaXYMax = 0.2;
  float dcaZMax = 2.0;

  Filter collisionFilter = nabs(aod::collision::posZ) < vtxZ;
  Filter trackFilter = nabs(aod::track::eta) < etaMax && aod::track::pt > ptMin&& aod::track::pt < ptMax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < dcaXYMax&& nabs(aod::track::dcaZ) < dcaZMax;
  using GeneralCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>;
  using UnfilteredTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;

  using UsedTracks = soa::Filtered<UnfilteredTracks>;
  using ZDCCollisions = soa::Filtered<soa::Join<GeneralCollisions, aod::SPTableZDC>>; // IMPORTANT: ZDCCollisions must be used to get access to ZDC q-vectors --> PWGCF/DataModel/SPTableZDC.h & PWGCF/Flow/TableProducer/zdcQVectors.cxx

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    AxisSpec axisPhi = {60, 0, constants::math::TwoPI, "#varphi"};
    AxisSpec axisEta = {64, -1.6, 1.6, "#eta"};
    AxisSpec axisEtaVn = {8, -.8, .8, "#eta"};
    AxisSpec axisCent = {90, 0, 90, "Centrality(%)"};
    AxisSpec axisPhiPlane = {100, -constants::math::PI, constants::math::PI, "#Psi"};
    AxisSpec axisQQ = {100, -0.2, 0.2, "#LT Q_{X}^{A}Q_{Y}^{C} #GT"};

    std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
    AxisSpec axisPt = {ptbinning, "#it{p}_{T} GeV/#it{c}"};

    rctChecker.init(rctFlags.cfgEvtRCTFlagCheckerLabel, rctFlags.cfgEvtRCTFlagCheckerZDCCheck, rctFlags.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    registry.add("hCentrality", "Centrality; #Events; Centrality (%)", {HistType::kTH1D, {axisCent}});
    registry.add("hSPplaneA", "Spectator Plane Angle A; #Events; #Psi_{A}", {HistType::kTH1D, {axisPhiPlane}});
    registry.add("qAqCXY", "Qx ZDCA x Qy ZDCC vs cent ; #ltQ_{X}^{A}Q_{Y}^{C}#GT; Centrality (%)", {HistType::kTH2D, {axisCent, axisQQ}});
    // Add here the histograms for the spectator plane angle C
    // registry.add(.........)
    registry.add("hSPplaneAvsC", "Spectator Plane Angle A vs C; #Events; #Psi_{A}; #Psi_{C}", {HistType::kTH2D, {axisPhiPlane, axisPhiPlane}});
    registry.add("hSPplaneFull", "Spectator Plane Angle Full; #Events; #Psi_{Full}", {HistType::kTH1D, {axisPhiPlane}});

    // Note: we will fill these with data from the CCDB, just to take a look!
    registry.add("CalibHistos/hQQx", "QQx; #Events; QQx", {HistType::kTProfile, {axisCent}});
    registry.add("CalibHistos/hQQy", "QQy; #Events; QQy", {HistType::kTProfile, {axisCent}});
    registry.add("CalibHistos/hQQ", "QQ; #Events; QQx + QQy", {HistType::kTProfile, {axisCent}});
    // Add here the histograms for the cross-terms of the Q-vectors from the ZDC

    registry.add("CalibHistos/hEvPlaneRes", "Event Plane Resolution; #Events; Event Plane Resolution", {HistType::kTProfile, {axisCent}});

    // Flow Histograms
    registry.add("flow/v1A", "", {HistType::kTProfile, {axisPt}});
    registry.add("flow/v1C", "", {HistType::kTProfile, {axisPt}});

    registry.add("flow/vnAxCxUxMH", "", {HistType::kTProfile, {axisCent}});
    registry.add("flow/vnAyCyUxMH", "", {HistType::kTProfile, {axisCent}});
    registry.add("flow/vnAxCyUyMH", "", {HistType::kTProfile, {axisCent}});
    registry.add("flow/vnAyCxUyMH", "", {HistType::kTProfile, {axisCent}});
  }

  void process(ZDCCollisions::iterator const& collision, aod::BCsWithTimestamps const&, UsedTracks const& tracks)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    float centrality = collision.centFT0C();
    registry.fill(HIST("hCentrality"), centrality);

    float centMin = 0;
    float centMax = 80;

    if (centrality > centMax || centrality < centMin)
      return;

    if (collision.isSelected() == false)
      return;

    // Get the ZDC q-vectors
    double qxA = collision.qxA();
    double qyA = collision.qyA();
    double qxC = collision.qxC();
    double qyC = collision.qyC();

    // Calculate the spectator plane angles
    double psiA = 1.0 * std::atan2(qyA, qxA);
    registry.fill(HIST("hSPplaneA"), psiA);

    // Add the PsiA vs PsiC as a TH2D

    double psiC = 1.0 * std::atan2(qyC, qxC);
    double psiFull = 1.0 * std::atan2(qyA + qyC, qxA + qxC);
    registry.fill(HIST("hSPplaneFull"), psiFull);

    registry.fill(HIST("hSPplaneAvsC"), psiA, psiC);

    // Fill the q-vector correlations
    registry.fill(HIST("qAqCXY"), centrality, qxA * qxC + qyA * qyC);

    double corrQQx = 1;
    double corrQQy = 1;
    double corrQQ = 1;
    double evPlaneRes = 1;

    // Get QQ-correlations from CCDB
    if (cfgCCDBdir_QQ.value.empty() == false) {
      TList* list = ccdb->getForTimeStamp<TList>(cfgCCDBdir_QQ.value, bc.timestamp());
      TProfile* qAqCX = reinterpret_cast<TProfile*>(list->FindObject("qAqCX"));
      TProfile* qAqCY = reinterpret_cast<TProfile*>(list->FindObject("qAqCY"));
      TProfile* qAqCXY = reinterpret_cast<TProfile*>(list->FindObject("qAqCXY"));
      // The sum is qAqCXY
      corrQQx = qAqCX->GetBinContent(centrality);
      corrQQy = qAqCY->GetBinContent(centrality);
      corrQQ = qAqCXY->GetBinContent(centrality);
      registry.fill(HIST("CalibHistos/hQQx"), centrality, corrQQx);
      registry.fill(HIST("CalibHistos/hQQy"), centrality, corrQQy);
      registry.fill(HIST("CalibHistos/hQQ"), centrality, corrQQ);
    }
    // Get event plane resolution from CCDB
    if (cfgCCDBdir_SP.value.empty() == false) {
      evPlaneRes = ccdb->getForTimeStamp<TProfile>(cfgCCDBdir_SP.value, bc.timestamp())->GetBinContent(centrality);
      registry.fill(HIST("CalibHistos/hEvPlaneRes"), centrality, evPlaneRes);
    }

    for (const auto& track : tracks) {

      // constrain angle to 0 -> [0,0+2pi]
      auto phi = RecoDecay::constrainAngle(track.phi(), 0);

      double ux = std::cos(phi);
      double uy = std::sin(phi);

      double uxMH = std::cos(2 * phi);
      double uyMH = std::sin(2 * phi);

      double v1A = (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ));
      double v1C = (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ));

      double v2AxCxUxMH = (uxMH * qxA * qxC) / corrQQx;
      double v2AyCyUxMH = (uxMH * qyA * qyC) / corrQQy;
      double v2AxCyUyMH = (uyMH * qxA * qyC) / corrQQx;
      double v2AyCxUyMH = (uyMH * qyA * qxC) / corrQQy;

      registry.fill(HIST("flow/v1A"), track.eta(), v1A);
      registry.fill(HIST("flow/v1C"), track.eta(), v1C);

      registry.fill(HIST("flow/v2AxCxUxMH"), centrality, v2AxCxUxMH);
      registry.fill(HIST("flow/v2AyCyUxMH"), centrality, v2AyCyUxMH);
      registry.fill(HIST("flow/v2AxCyUyMH"), centrality, v2AxCyUyMH);
      registry.fill(HIST("flow/v2AyCxUyMH"), centrality, v2AyCxUyMH);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SpectatorPlaneTutorial>(cfgc)};
}
