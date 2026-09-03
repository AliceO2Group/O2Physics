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
/// \file ptmultCorr.cxx
///
/// \brief task for analysis of charged-particle RAA at midrapidity in light-ion collisions
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since October 01, 2025

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;
using namespace o2::aod::rctsel;

using ColDataTablePbPb = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Cs>;
using ColDataTablepp = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults>;
using ColMCRecTablePbPb = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Cs>>;
using ColMCRecTablepp = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults>>;
using ColMCTrueTable = soa::Join<aod::McCollisions, aod::MultMCExtras>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>;
using TrackMCTrueTable = aod::McParticles;

enum {
  kTrackTypebegin = 0,
  kGlobalplusITS,
  kGlobalonly,
  kITSonly,
  kTrackTypeend
};

enum {
  kGenpTbegin = 0,
  kNoGenpTVar,
  kGenpTup,
  kGenpTdown,
  kGenpTend
};

// ── CHANGE 1 ────────────────────────────────────────────────────────────────
// Added kGenAllCharged as the new last bin before kGenTrkTypeend.
// Numerator  : kGenAll        → MC truth physical primaries
// Denominator: kGenAllCharged → MC truth all charged (primaries + secondaries)
// Primary fraction from MC = kGenAll / kGenAllCharged per pT bin.
// axisGenTrkType (built from kGenTrkTypeend) auto-expands from 5 → 6 bins.
// ────────────────────────────────────────────────────────────────────────────
enum {
  kGenTrkTypebegin = 0,
  kGenAll,
  kGenPion,
  kGenKaon,
  kGenProton,
  kGenOther,
  kGenAllCharged, // all charged particles (primaries + secondaries) — denominator for primary fraction
  kGenTrkTypeend
};

enum {
  kRecTrkTypebegin = 0,
  kRecoAll,
  kRecoHasmc,
  kRecoPrimary,
  kRecoPion,
  kRecoKaon,
  kRecoProton,
  kRecoOther,
  kRecoSecondary,
  kRecoWeakDecay,
  kRecoMaterial,
  kRecoFake,
  kRecoBkg,
  kRecTrkTypeend
};

const AxisSpec axisEvent{15, 0.5, 15.5, "#Event", "EventAxis"};
const AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};
const AxisSpec axisEta{40, -2, 2, "#eta", "EtaAxis"};
const AxisSpec axisPhi{{0, o2::constants::math::PIQuarter, o2::constants::math::PIHalf, o2::constants::math::PIQuarter * 3., o2::constants::math::PI, o2::constants::math::PIQuarter * 5., o2::constants::math::PIHalf * 3., o2::constants::math::PIQuarter * 7., o2::constants::math::TwoPI}, "#phi", "PhiAxis"};
const AxisSpec axisPhi2{629, 0, o2::constants::math::TwoPI, "#phi"};
const AxisSpec axisCent{100, 0, 100, "#Cent"};
const AxisSpec axisDeltaPt{50, -1.0, +1.0, "#Delta(pT)"};
const AxisSpec axisDCAxy{600, -3.0, 3.0, "DCA_{xy} (cm)", "DCAxyAxis"}; // range ±3 cm for TPC-only case
const AxisSpec axisGenPtVary = {kGenpTend - 1, +kGenpTbegin + 0.5, +kGenpTend - 0.5, "", "GenpTVaryAxis"};
const AxisSpec axisGenTrkType = {kGenTrkTypeend - 1, +kGenTrkTypebegin + 0.5, +kGenTrkTypeend - 0.5, "", "GenTrackTypeAxis"};
const AxisSpec axisRecTrkType = {kRecTrkTypeend - 1, +kRecTrkTypebegin + 0.5, +kRecTrkTypeend - 0.5, "", "RecTrackTypeAxis"};
const AxisSpec axisTrackType = {kTrackTypeend - 1, +kTrackTypebegin + 0.5, +kTrackTypeend - 0.5, "", "TrackTypeAxis"};
auto static constexpr KminCharge = 3.f;

struct PtmultCorr {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg{};
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;
  Configurable<float> etaRange{"etaRange", 0.8f, "Eta range to consider"};
  Configurable<float> vtxRange{"vtxRange", 10.0f, "Vertex Z range to consider"};
  Configurable<float> occuRange{"occuRange", 500.0f, "Occupancy range to consider"};
  Configurable<float> cfgPtCutMin{"cfgPtCutMin", 0.1f, "minimum accepted track pT"};
  Configurable<float> cfgPtCutMax{"cfgPtCutMax", 1e9f, "maximum accepted track pT"};
  Configurable<float> extraphicut1{"extraphicut1", 3.07666f, "Extra Phi cut 1"};
  Configurable<float> extraphicut2{"extraphicut2", 3.12661f, "Extra Phi cut 2"};
  Configurable<float> extraphicut3{"extraphicut3", 0.03f, "Extra Phi cut 3"};
  Configurable<float> extraphicut4{"extraphicut4", 6.253f, "Extra Phi cut 4"};
  Configurable<bool> isZvtxPosSelMC{"isZvtxPosSelMC", true, "Zvtx position selection for MC events"};

  Configurable<int16_t> cfgMinNCrossedRows{"cfgMinNCrossedRows", 70, "Minimum TPC crossed rows"};
  Configurable<bool> cfgUseNclsPID{"cfgUseNclsPID", false, "Use NclsPID instead of NclsFound for the Ncls cut"};
  Configurable<float> cfgMinChi2ClsTPC{"cfgMinChi2ClsTPC", 0.5f, "Minimum TPC chi2/cls"};
  Configurable<float> cfgMaxChi2ClsTPC{"cfgMaxChi2ClsTPC", 4.0f, "Maximum TPC chi2/cls"};
  Configurable<float> cfgMaxChi2ClsITS{"cfgMaxChi2ClsITS", 36.0f, "Maximum ITS chi2/cls"};
  Configurable<float> cfgDCAz{"cfgDCAz", 0.04, "DCAz cut in cm"};

  ConfigurableAxis ptHistBin{"ptHistBin", {200, 0., 20.}, ""};
  ConfigurableAxis centralityBinning{"centralityBinning", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, ""};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {VARIABLE_WIDTH, 0.0, 3.00065, 4.28798, 6.14552, 7.6196, 8.90942, 10.0897, 11.2002, 12.2709, 13.3167, 14.4173, 23.2518}, "Binning of the impact parameter axis"};

  Configurable<bool> isSameBunchPileup{"isSameBunchPileup", false, "Enable SameBunchPileup cut"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", false, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> applyExtraPhiCut{"applyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "Enable NoCollInTimeRangeStandard cut"};
  Configurable<bool> isNoCollInTimeRangeStrict{"isNoCollInTimeRangeStrict", false, "use isNoCollInTimeRangeStrict?"};

  Configurable<bool> selHasBC{"selHasBC", true, "Require has_foundBC"};
  Configurable<bool> selHasFT0{"selHasFT0", true, "Require has_foundFT0"};

  Configurable<bool> isNoCollInRofStandard{"isNoCollInRofStandard", false, "Enable NoCollInRofStandard cut"};
  Configurable<bool> isNoCollInRofStrict{"isNoCollInRofStrict", false, "use isNoCollInRofStrict?"};

  Configurable<bool> isNoHighMultCollInPrevRof{"isNoHighMultCollInPrevRof", false, "use isNoHighMultCollInPrevRof?"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", false, "use isNoCollInTimeRangeNarrow?"};

  Configurable<bool> applyFT0CbasedOccupancy{"applyFT0CbasedOccupancy", false, "Enable FT0CbasedOccupancy cut"};
  Configurable<bool> applyInelgt0{"applyInelgt0", true, "Enable INEL > 0 condition"};
  Configurable<bool> applyOccuCut{"applyOccuCut", true, "Enable occupancy selection"};

  Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "RCT selection flag"};
  Configurable<bool> rctCheckZDC{"rctCheckZDC", false, "Check ZDC in RCT"};
  Configurable<bool> rctTreatLimitedAcceptanceAsBad{"rctTreatLimitedAcceptanceAsBad", false, "Treat limited acceptance as bad"};
  Configurable<bool> requireGoodRct{"requireGoodRct", true, "Apply RCT selection"};

  // RCT checker instance
  RCTFlagsChecker rctChecker;

  static constexpr int kZeroInt{0};
  static constexpr float kOne{1.0f};

  void init(InitContext const&)
  {
    if (requireGoodRct) {
      rctChecker.init(rctLabel.value, rctCheckZDC.value, rctTreatLimitedAcceptanceAsBad.value);
    }
    AxisSpec centAxis = {centralityBinning, "Centrality", "CentralityAxis"};
    AxisSpec axisPt = {ptHistBin, "pT", "pTAxis"};
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter"};

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);
    histos.add("CentPercentileHist", "CentPercentileHist", kTH1D, {axisCent}, false);
    histos.add("PhiVsEtaHistNoCut", "PhiVsEtaHistNoCut", kTH2D, {axisPhi2, axisEta}, false);
    histos.add("PhiVsEtaHistWithCut", "PhiVsEtaHistWithCut", kTH2D, {axisPhi2, axisEta}, false);

    if (doprocessDataPbPb) {
      histos.add("hdatazvtxcent", "hdatazvtxcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hdatahistPbPb", "hdatahistPbPb", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisTrackType, axisDCAxy}, false);
    }

    if (doprocessDatapp) {
      histos.add("hdatahistpp", "hdatahistpp", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisTrackType, axisDCAxy}, false);
    }

    if (doprocessMCeffPbPb) {
      histos.add("hPbPbGenMCvtxz", "hPbPbGenMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hPbPbGenMCvtxzcent", "hPbPbGenMCvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hPbPbGenMCAssoRecvtxz", "hPbPbGenMCAssoRecvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hPbPbGenMCAssoRecvtxzcent", "hPbPbGenMCAssoRecvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hPbPbGenMCdndpt", "hPbPbGenMCdndpt", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi}, false);
      // axisGenTrkType now has 6 bins: kGenAll(1)…kGenOther(5), kGenAllCharged(6)
      histos.add("hPbPbGenMCAssoRecdndpt", "hPbPbGenMCAssoRecdndpt", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisGenTrkType, axisGenPtVary}, false);

      histos.add("hPbPbRecMCvtxz", "hPbPbRecMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hPbPbRecMCvtxzcent", "hPbPbRecMCvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hPbPbRecMCcent", "hPbPbRecMCcent", kTH1D, {axisCent}, false);
      histos.add("hPbPbRecMCdndpt", "hPbPbRecMCdndpt", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisRecTrkType, axisTrackType, axisDCAxy}, false);
      histos.add("hPbPbEtaReso", "hPbPbEtaReso", kTH2D, {axisPt, axisDeltaPt});
    }

    if (doprocessMCeffpp) {
      histos.add("hppGenMCvtxz", "hppGenMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppGenMCAssoRecvtxz", "hppGenMCAssoRecvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppGenMCdndpt", "hppGenMCdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi}, false);
      // axisGenTrkType now has 6 bins: kGenAll(1)…kGenOther(5), kGenAllCharged(6)
      histos.add("hppGenMCAssoRecdndpt", "hppGenMCAssoRecdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisGenTrkType, axisGenPtVary}, false);

      histos.add("hppRecMCvtxz", "hppRecMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppRecMCdndpt", "hppRecMCdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisRecTrkType, axisTrackType, axisDCAxy}, false);
      histos.add("hppEtaReso", "hppEtaReso", kTH2D, {axisPt, axisDeltaPt});
    }

    if (doprocessEvtLossSigLossMCpp || doprocessEvtLossSigLossMCPbPb) {
      histos.add("MCEventHist", "MCEventHist", kTH1F, {axisEvent}, false);
      auto hstat = histos.get<TH1>(HIST("MCEventHist"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All MC events");
      x->SetBinLabel(2, "MC events with atleast one reco event");
      histos.add("hgenptBeforeEvtSel", "hgenptBeforeEvtSel", kTH1F, {axisPt}, false);
      histos.add("hgenptAfterEvtSel", "hgenptAfterEvtSel", kTH1F, {axisPt}, false);
      histos.add("hgenptBeforeEvtSelPbPb", "hgenptBeforeEvtSelPbPb", kTH2F, {axisPt, impactParAxis}, false);
      histos.add("hgenptAfterEvtSelPbPb", "hgenptAfterEvtSelPbPb", kTH2F, {axisPt, impactParAxis}, false);
      histos.add("hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {axisCent, impactParAxis});

      // CO map: centrality vs generated Nch (for N_rec_at_least_once events)
      histos.add("hNchVsCent_WithRecoEvt", "hNchVsCent_WithRecoEvt",
                 kTH2F, {axisCent, {500, 0, 500, "Gen N_{ch} |#eta|<0.8"}});

      // Event loss denominator: Nch distribution of all generated events
      histos.add("hNchMC_AllGen", "hNchMC_AllGen",
                 kTH1F, {{500, 0, 500, "Gen N_{ch} |#eta|<0.8"}});

      // Event loss numerator: Nch distribution of N_rec_at_least_once events
      histos.add("hNchMC_WithRecoEvt", "hNchMC_WithRecoEvt",
                 kTH1F, {{500, 0, 500, "Gen N_{ch} |#eta|<0.8"}});

      // Event splitting: centrality distributions
      histos.add("hCent_WRecoEvtWSelCri", "N_rec_at_least_once (with sel)", kTH1F, {axisCent});
      histos.add("hCent_AllRecoEvt", "N_rec (with sel)", kTH1F, {axisCent});

      // Signal loss: pT vs Nch (2D), for all gen events and N_rec_at_least_once
      histos.add("hPtVsNchMC_AllGen", "hPtVsNchMC_AllGen",
                 kTH2F, {{axisPt}, {500, 0, 500, "Gen N_{ch} |#eta|<0.8"}});
      histos.add("hPtVsNchMC_WithRecoEvt", "hPtVsNchMC_WithRecoEvt",
                 kTH2F, {{axisPt}, {500, 0, 500, "Gen N_{ch} |#eta|<0.8"}});
    }

    if (doprocessEvtLossSigLossMCpp) {
      histos.add("hNch_AllRecoEvt",
                 "All reco collisions passing partial cuts (event split denom pp)",
                 kTH1F, {{500, 0, 500, "N_{ch}  |#eta|<0.8"}});
      histos.add("hNch_WRecoEvtWSelCri",
                 "Best reco collision passing full selection (event split num pp)",
                 kTH1F, {{500, 0, 500, "N_{ch}  |#eta|<0.8"}});
    }
    if (doprocessDatapp || doprocessDataPbPb) {
      histos.add("RCTSel", "All=1 | RCT passed=2", kTH1F, {{2, 0.5, 2.5}});
      auto hrct = histos.get<TH1>(HIST("RCTSel"));
      hrct->GetXaxis()->SetBinLabel(1, "All");
      hrct->GetXaxis()->SetBinLabel(2, "RCT passed");
    }
    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "has_foundBC");
    x->SetBinLabel(3, "has_foundFT0");
    x->SetBinLabel(4, "kIsTriggerTVX");
    x->SetBinLabel(5, "kNoITSROFrameBorder");
    x->SetBinLabel(6, "kNoTimeFrameBorder");
    x->SetBinLabel(7, "|vz|<vtxRange");
    x->SetBinLabel(8, "kIsGoodZvtxFT0vsPV");
    x->SetBinLabel(9, "kNoSameBunchPileup");
    x->SetBinLabel(10, "kNoCollInTimeRangeStrict");
    x->SetBinLabel(11, "kNoCollInRofStrict");
    x->SetBinLabel(12, "kNoHighMultCollInPrevRof");
    x->SetBinLabel(13, "INEL>0");
    x->SetBinLabel(14, "Occupancy");
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    histos.fill(HIST("EventHist"), 1);

    if (!col.sel8()) {
      return false;
    }
    // TVX trigger (replaces sel8 as primary trigger requirement)
    // if (!col.selection_bit(o2::aod::evsel::kIsTriggerTVX))
    // return false;
    // histos.fill(HIST("EventHist"), 4);

    // ITS ROF border
    if (!col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 5);

    // Time frame border
    if (!col.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 6);

    // vtxZ
    if (std::abs(col.posZ()) > vtxRange) {
      return false;
    }
    histos.fill(HIST("EventHist"), 7);

    // Good ZvtxFT0vsPV
    if (isGoodZvtxFT0vsPV &&
        !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 8);

    // No same bunch pileup
    if (isSameBunchPileup &&
        !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 9);

    // Time range isolation — Strict (replaces Standard)
    if (isNoCollInTimeRangeStrict &&
        !col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 10);

    // ROF isolation — Strict (replaces Standard)
    if (isNoCollInRofStrict &&
        !col.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 11);

    // No high mult collision in previous ROF
    if (isNoHighMultCollInPrevRof &&
        !col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 12);

    // INEL > 0
    if (applyInelgt0 && !col.isInelGt0()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 13);

    // Occupancy
    auto occu = applyFT0CbasedOccupancy
                  ? col.ft0cOccupancyInTimeRange()
                  : col.trackOccupancyInTimeRange();
    if (applyOccuCut && occu > occuRange) {
      return false;
    }
    histos.fill(HIST("EventHist"), 14);

    return true;
  }

  template <typename CheckTrack>
  bool isTrackSelected(CheckTrack const& track)
  {
    // --- eta (applied to all track types) ---
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }

    // --- pT (applied to all track types) ---
    if (track.pt() < cfgPtCutMin || track.pt() > cfgPtCutMax) {
      return false;
    }

    // ITS inner-barrel hit: must fire layer 0 or layer 1
    if (!(track.itsClusterMap() & 0x01) && !(track.itsClusterMap() & 0x02)) {
      return false;
    }

    // ITS chi2 cut (also part of ITS cuts)
    if (track.itsChi2NCl() > cfgMaxChi2ClsITS) {
      return false;
    }

    if (track.hasTPC()) {

      if (track.tpcNClsCrossedRows() < cfgMinNCrossedRows) {
        return false;
      }

      if (track.tpcChi2NCl() < cfgMinChi2ClsTPC || track.tpcChi2NCl() > cfgMaxChi2ClsTPC) {
        return false;
      }
    }

    if (std::abs(track.dcaZ()) > cfgDCAz) {
      return false;
    }

    // --- optional phi cut (applied to all track types) ---
    histos.fill(HIST("PhiVsEtaHistNoCut"), track.phi(), track.eta());
    if (applyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    histos.fill(HIST("PhiVsEtaHistWithCut"), track.phi(), track.eta());

    return true;
  }

  // Selects MC truth physical primaries (numerator for primary fraction).
  // Requires: isPhysicalPrimary + producedByGenerator + charged + |eta| < etaRange
  template <typename CheckGenTrack>
  bool isGenTrackSelected(CheckGenTrack const& track)
  {
    if (!track.isPhysicalPrimary()) {
      return false;
    }
    if (track.pt() < cfgPtCutMin || track.pt() > cfgPtCutMax) {
      return false;
    }
    auto pdgTrack = pdg->GetParticle(track.pdgCode());
    if (pdgTrack == nullptr) {
      return false;
    }
    if (std::abs(pdgTrack->Charge()) < KminCharge) {
      return false;
    }
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (applyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    return true;
  }

  // ── CHANGE 2 ──────────────────────────────────────────────────────────────
  // New helper: selects ALL charged particles in acceptance (denominator).
  // Identical to isGenTrackSelected but WITHOUT the isPhysicalPrimary() check,
  // so secondaries (weak-decay daughters, material interactions) also pass.
  // Used to fill kGenAllCharged → primary fraction = kGenAll / kGenAllCharged.
  // ──────────────────────────────────────────────────────────────────────────
  template <typename CheckGenTrack>
  bool isGenChargedTrackSelected(CheckGenTrack const& track)
  {
    if (track.pt() < cfgPtCutMin || track.pt() > cfgPtCutMax) {
      return false;
    }
    auto pdgTrack = pdg->GetParticle(track.pdgCode());
    if (pdgTrack == nullptr) {
      return false;
    }
    if (std::abs(pdgTrack->Charge()) < KminCharge) {
      return false;
    }
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (applyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    return true;
  }

  void processDataPbPb(ColDataTablePbPb::iterator const& cols, TrackDataTable const& tracks)
  {
    if (requireGoodRct) {
      histos.fill(HIST("RCTSel"), 1); // all events
      if (!rctChecker(cols)) {
        return;
      }
      histos.fill(HIST("RCTSel"), 2); // passed RCT
    }
    if (!isEventSelected(cols)) {
      return;
    }
    histos.fill(HIST("VtxZHist"), cols.posZ());
    histos.fill(HIST("CentPercentileHist"), cols.centFT0C());
    histos.fill(HIST("hdatazvtxcent"), cols.posZ(), cols.centFT0C());

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      histos.fill(HIST("hdatahistPbPb"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), kGlobalplusITS, track.dcaXY());
      if (track.hasTPC()) {
        histos.fill(HIST("hdatahistPbPb"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), kGlobalonly, track.dcaXY());
      } else {
        histos.fill(HIST("hdatahistPbPb"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), kITSonly, track.dcaXY());
      }
    }
  }

  void processDatapp(ColDataTablepp::iterator const& cols, TrackDataTable const& tracks)
  {
    if (requireGoodRct) {
      histos.fill(HIST("RCTSel"), 1); // all events
      if (!rctChecker(cols)) {
        return;
      }
      histos.fill(HIST("RCTSel"), 2); // passed RCT
    }
    if (!isEventSelected(cols)) {
      return;
    }
    histos.fill(HIST("VtxZHist"), cols.posZ());

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      histos.fill(HIST("hdatahistpp"), cols.posZ(), track.pt(), track.phi(), kGlobalplusITS, track.dcaXY());
      if (track.hasTPC()) {
        histos.fill(HIST("hdatahistpp"), cols.posZ(), track.pt(), track.phi(), kGlobalonly, track.dcaXY());
      } else {
        histos.fill(HIST("hdatahistpp"), cols.posZ(), track.pt(), track.phi(), kITSonly, track.dcaXY());
      }
    }
  }

  void processMCeffPbPb(ColMCTrueTable::iterator const& mcCollision, ColMCRecTablePbPb const& RecCols, TrackMCTrueTable const& GenParticles, TrackMCRecTable const& RecTracks)
  {
    float gencent = -999.f;
    bool atLeastOne = false;
    int bestCollisionIndex = -1;
    int biggestNContribs = -1;
    for (const auto& RecCol : RecCols) {
      if (biggestNContribs < RecCol.numContrib()) {
        biggestNContribs = RecCol.numContrib();
        bestCollisionIndex = RecCol.globalIndex();
      }
    }
    for (const auto& RecCol : RecCols) {

      if (!isEventSelected(RecCol)) {
        continue;
      }

      if (RecCol.globalIndex() != bestCollisionIndex) {
        continue;
      }

      atLeastOne = true;
      gencent = RecCol.centFT0C();
    }
    histos.fill(HIST("hPbPbGenMCvtxz"), mcCollision.posZ());
    histos.fill(HIST("hPbPbGenMCvtxzcent"), mcCollision.posZ(), gencent);
    if (atLeastOne) {
      histos.fill(HIST("hPbPbGenMCAssoRecvtxz"), mcCollision.posZ());
      histos.fill(HIST("hPbPbGenMCAssoRecvtxzcent"), mcCollision.posZ(), gencent);
    }

    // ── CHANGE 3 ────────────────────────────────────────────────────────────
    // Restructured PbPb gen particle loop for primary fraction from MC truth.
    //
    // Pass 1 (loose): isGenChargedTrackSelected — all charged in acceptance.
    //   → fills kGenAllCharged (denominator) when associated with a reco event.
    //
    // Pass 2 (tight): additionally require isPhysicalPrimary (via continue).
    //   → fills kGenAll and species bins (numerator) — identical to before.
    //
    // Primary fraction (PbPb) = hPbPbGenMCAssoRecdndpt[kGenAll]
    //                         / hPbPbGenMCAssoRecdndpt[kGenAllCharged]  per pT bin
    // ────────────────────────────────────────────────────────────────────────
    for (const auto& particle : GenParticles) {
      // Loose check: charged + |eta| + phi (no isPhysicalPrimary requirement)
      if (!isGenChargedTrackSelected(particle)) {
        continue;
      }
      // Fill denominator: all charged particles (primaries + secondaries)
      if (atLeastOne) {
        histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAllCharged), kNoGenpTVar);
      }
      // Tight check: physical primaries only — skip secondaries from this point on
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      // Fill primary-only histograms (numerator + species breakdown)
      histos.fill(HIST("hPbPbGenMCdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi());
      if (atLeastOne) {
        histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kNoGenpTVar);
        if (particle.pt() < cfgPtCutMin) {
          histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTup);
          histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown);
        }
        int pid = 0;
        switch (std::abs(particle.pdgCode())) {
          case PDG_t::kPiPlus:
            pid = kGenPion;
            break;
          case PDG_t::kKPlus:
            pid = kGenKaon;
            break;
          case PDG_t::kProton:
            pid = kGenProton;
            break;
          default:
            pid = kGenOther;
            break;
        }
        histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(pid), kNoGenpTVar);
      } // Associated with reco col
    } // track (mcgen) loop

    for (const auto& RecCol : RecCols) {

      if (!isEventSelected(RecCol)) {
        continue;
      }

      if (RecCol.globalIndex() != bestCollisionIndex) {
        continue;
      }

      histos.fill(HIST("hPbPbRecMCvtxz"), RecCol.posZ());
      histos.fill(HIST("hPbPbRecMCcent"), RecCol.centFT0C());
      histos.fill(HIST("hPbPbRecMCvtxzcent"), RecCol.posZ(), RecCol.centFT0C());
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      std::vector<int> mclabels;
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        auto trkType = Rectrack.hasTPC() ? kGlobalonly : kITSonly;

        histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoAll), kGlobalplusITS, Rectrack.dcaXY());
        histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoAll), trkType, Rectrack.dcaXY());

        if (Rectrack.has_mcParticle()) {
          int pid = 0;
          auto mcpart = Rectrack.mcParticle();
          histos.fill(HIST("hPbPbEtaReso"), Rectrack.pt(), Rectrack.pt() - mcpart.pt());
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoHasmc), kGlobalplusITS, Rectrack.dcaXY());
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoHasmc), trkType, Rectrack.dcaXY());
          if (mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoPrimary), kGlobalplusITS, Rectrack.dcaXY());
            histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoPrimary), trkType, Rectrack.dcaXY());
            switch (std::abs(mcpart.pdgCode())) {
              case PDG_t::kPiPlus:
                pid = kRecoPion;
                break;
              case PDG_t::kKPlus:
                pid = kRecoKaon;
                break;
              case PDG_t::kProton:
                pid = kRecoProton;
                break;
              default:
                pid = kRecoOther;
                break;
            }
          } else {
            // non-primary → split weak-decay vs material by production mechanism
            if (mcpart.getProcess() == TMCProcess::kPDecay) {
              pid = kRecoWeakDecay;
            } else {
              pid = kRecoMaterial;
            }
            // also fill the union bin as a cross-check (weak + material)
            histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoSecondary), kGlobalplusITS, Rectrack.dcaXY());
            histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoSecondary), trkType, Rectrack.dcaXY());
          }

          // duplicate-track (split) relabel — overrides origin, as before
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kRecoFake;
          }
          mclabels.push_back(Rectrack.mcParticleId());

          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(pid), kGlobalplusITS, Rectrack.dcaXY());
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(pid), trkType, Rectrack.dcaXY());
        } else {
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoBkg), kGlobalplusITS, Rectrack.dcaXY());
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoBkg), trkType, Rectrack.dcaXY());
        }
      }
    } // collision loop
  }

  void processMCeffpp(ColMCTrueTable::iterator const& mcCollision, ColMCRecTablepp const& RecCols, TrackMCTrueTable const& GenParticles, TrackMCRecTable const& RecTracks)
  {
    bool atLeastOne = false;
    int bestCollisionIndex = -1;
    int biggestNContribs = -1;
    for (const auto& RecCol : RecCols) {
      if (biggestNContribs < RecCol.numContrib()) {
        biggestNContribs = RecCol.numContrib();
        bestCollisionIndex = RecCol.globalIndex();
      }
    }
    for (const auto& RecCol : RecCols) {

      if (!isEventSelected(RecCol)) {
        continue;
      }

      if (RecCol.globalIndex() != bestCollisionIndex) {
        continue;
      }

      atLeastOne = true;
    }
    histos.fill(HIST("hppGenMCvtxz"), mcCollision.posZ());
    if (atLeastOne) {
      histos.fill(HIST("hppGenMCAssoRecvtxz"), mcCollision.posZ());
    }

    // ── CHANGE 4 ────────────────────────────────────────────────────────────
    // Restructured pp gen particle loop for primary fraction from MC truth.
    //   Pass 1 (loose): isGenChargedTrackSelected → fills kGenAllCharged (denom).
    //   Pass 2 (tight):  + isPhysicalPrimary      → fills kGenAll + species (num).
    //   Primary fraction (pp) = kGenAll / kGenAllCharged per pT bin.
    // ────────────────────────────────────────────────────────────────────────
    for (const auto& particle : GenParticles) {
      if (!isGenChargedTrackSelected(particle)) {
        continue;
      }
      if (atLeastOne) {
        histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAllCharged), kNoGenpTVar);
      }
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      histos.fill(HIST("hppGenMCdndpt"), mcCollision.posZ(), particle.pt(), particle.phi());
      if (atLeastOne) {
        histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kNoGenpTVar);
        if (particle.pt() < cfgPtCutMin) {
          histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTup);
          histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown);
        }
        int pid = 0;
        switch (std::abs(particle.pdgCode())) {
          case PDG_t::kPiPlus:
            pid = kGenPion;
            break;
          case PDG_t::kKPlus:
            pid = kGenKaon;
            break;
          case PDG_t::kProton:
            pid = kGenProton;
            break;
          default:
            pid = kGenOther;
            break;
        }
        histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(pid), kNoGenpTVar);
      } // Associated with reco col
    } // track (mcgen) loop

    for (const auto& RecCol : RecCols) {

      if (!isEventSelected(RecCol)) {
        continue;
      }

      if (RecCol.globalIndex() != bestCollisionIndex) {
        continue;
      }

      histos.fill(HIST("hppRecMCvtxz"), RecCol.posZ());
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      std::vector<int> mclabels;
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        auto trkType = Rectrack.hasTPC() ? kGlobalonly : kITSonly;
        histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoAll), kGlobalplusITS, Rectrack.dcaXY());
        histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoAll), trkType, Rectrack.dcaXY());
        if (Rectrack.has_mcParticle()) {
          int pid = 0;
          auto mcpart = Rectrack.mcParticle();
          histos.fill(HIST("hppEtaReso"), Rectrack.pt(), Rectrack.pt() - mcpart.pt());
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoHasmc), kGlobalplusITS, Rectrack.dcaXY());
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoHasmc), trkType, Rectrack.dcaXY());
          if (mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoPrimary), kGlobalplusITS, Rectrack.dcaXY());
            histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoPrimary), trkType, Rectrack.dcaXY());
            switch (std::abs(mcpart.pdgCode())) {
              case PDG_t::kPiPlus:
                pid = kRecoPion;
                break;
              case PDG_t::kKPlus:
                pid = kRecoKaon;
                break;
              case PDG_t::kProton:
                pid = kRecoProton;
                break;
              default:
                pid = kRecoOther;
                break;
            }
          } else {
            // non-primary → split weak-decay vs material by production mechanism
            if (mcpart.getProcess() == TMCProcess::kPDecay) {
              pid = kRecoWeakDecay;
            } else {
              pid = kRecoMaterial;
            }
            // union cross-check bin (weak + material)
            histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoSecondary), kGlobalplusITS, Rectrack.dcaXY());
            histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoSecondary), trkType, Rectrack.dcaXY());
          }
          // duplicate-track (split) relabel — overrides origin, as before
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kRecoFake;
          }
          mclabels.push_back(Rectrack.mcParticleId());
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(pid), kGlobalplusITS, Rectrack.dcaXY());
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(pid), trkType, Rectrack.dcaXY());
        } else {
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoBkg), kGlobalplusITS, Rectrack.dcaXY());
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoBkg), trkType, Rectrack.dcaXY());
        }
      } // track (mcrec) loop
    } // collision loop
  }
  void processEvtLossSigLossMCpp(ColMCTrueTable::iterator const& mcCollision, ColMCRecTablepp const& RecCols, TrackMCTrueTable const& GenParticles)
  {

    // ── Count generated Nch in |eta| < etaRange for event loss ──────────────
    int nChMC = 0;
    int nChMCEta1 = 0;
    for (const auto& particle : GenParticles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle) {
        continue;
      }
      if (std::abs(pdgParticle->Charge()) < KminCharge) {
        continue;
      }
      if (std::abs(particle.eta()) < etaRange) {
        nChMC++;
      }
      if (std::abs(particle.eta()) > kOne) {
        continue;
      }
      nChMCEta1++;
    }

    if (isZvtxPosSelMC && (std::fabs(mcCollision.posZ()) > vtxRange)) {
      return;
    }

    //---------------------------
    // Select only INEL > 0 generated events?
    //---------------------------
    if (applyInelgt0) {
      if (!(nChMCEta1 > kZeroInt)) {
        return;
      }
    }

    bool atLeastOne = false;
    int bestCollisionIndex = -1;
    int biggestNContribs = -1;

    for (const auto& RecCol : RecCols) {

      if (biggestNContribs < RecCol.numContrib()) {
        biggestNContribs = RecCol.numContrib();
        bestCollisionIndex = RecCol.globalIndex();
      }

      if (!RecCol.sel8()) {
        continue;
      }

      if (!RecCol.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }
      if (!RecCol.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }
      if (std::fabs(RecCol.posZ()) > vtxRange) {
        continue;
      }
      if (isGoodZvtxFT0vsPV &&
          !RecCol.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        continue;
      }
      if (isSameBunchPileup &&
          !RecCol.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        continue;
      }

      histos.fill(HIST("hNch_AllRecoEvt"), nChMC); // denominator for event splitting
    }

    // ── Find best collision, apply full selection ────────────────────────────

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }

      if (RecCol.globalIndex() != bestCollisionIndex) {
        continue;
      }
      atLeastOne = true;
      histos.fill(HIST("hNch_WRecoEvtWSelCri"), nChMC); // numerator for event splitting
    }

    // ── Event loss histograms ────────────────────────────────────────────────
    histos.fill(HIST("MCEventHist"), 1);
    histos.fill(HIST("hNchMC_AllGen"), nChMC); // denominator for event loss

    if (atLeastOne) {
      histos.fill(HIST("MCEventHist"), 2);
      histos.fill(HIST("hNchMC_WithRecoEvt"), nChMC); // numerator for event loss
    }

    // ── Signal loss particle loop ────────────────────────────────────────────
    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }

      histos.fill(HIST("hPtVsNchMC_AllGen"), particle.pt(), nChMC); // denominator: all generated events
      histos.fill(HIST("hgenptBeforeEvtSel"), particle.pt());

      if (atLeastOne) {
        histos.fill(HIST("hPtVsNchMC_WithRecoEvt"), particle.pt(), nChMC); // numerator: ≥1 reco collision passing selection
        histos.fill(HIST("hgenptAfterEvtSel"), particle.pt());
      }
    }
  }

  void processEvtLossSigLossMCPbPb(ColMCTrueTable::iterator const& mcCollision, ColMCRecTablePbPb const& RecCols, TrackMCTrueTable const& GenParticles)
  {
    // ── Count generated Nch in |eta| < etaRange for event loss ──────────────
    int nChMC = 0;

    int nChMCEta1 = 0;
    for (const auto& particle : GenParticles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle) {
        continue;
      }
      if (std::abs(pdgParticle->Charge()) < KminCharge) {
        continue;
      }
      if (std::abs(particle.eta()) < etaRange) {
        nChMC++;
      }

      if (std::abs(particle.eta()) > kOne) {
        continue;
      }
      nChMCEta1++;
    }

    if (isZvtxPosSelMC && (std::fabs(mcCollision.posZ()) > vtxRange)) {
      return;
    }

    //---------------------------
    // Select only INEL > 0 generated events?
    //---------------------------
    if (applyInelgt0) {
      if (!(nChMCEta1 > kZeroInt)) {
        return;
      }
    }

    // ── Event splitting denominator ──────────────────────────────────────────
    bool atLeastOne = false;
    int bestCollisionIndex = -1;
    int biggestNContribs = -1;

    for (const auto& RecCol : RecCols) {

      if (biggestNContribs < RecCol.numContrib()) {
        biggestNContribs = RecCol.numContrib();
        bestCollisionIndex = RecCol.globalIndex();
      }

      if (!RecCol.sel8()) {
        continue;
      }

      if (!RecCol.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }
      if (!RecCol.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }
      if (std::fabs(RecCol.posZ()) > vtxRange) {
        continue;
      }
      if (isGoodZvtxFT0vsPV &&
          !RecCol.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        continue;
      }
      if (isSameBunchPileup &&
          !RecCol.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        continue;
      }

      histos.fill(HIST("hCent_AllRecoEvt"), RecCol.centFT0C()); // denominator for event splitting
    }

    // ── Find best collision, apply full selection ────────────────────────────

    auto centrality = -999.f;
    for (const auto& RecCol : RecCols) {
      if (RecCol.globalIndex() != bestCollisionIndex) {
        continue;
      }

      if (!isEventSelected(RecCol)) {
        continue;
      }

      atLeastOne = true;
      centrality = RecCol.centFT0C();
      histos.fill(HIST("hCent_WRecoEvtWSelCri"), RecCol.centFT0C()); // numerator for event splitting
    }

    // ── Event loss histograms ────────────────────────────────────────────────
    histos.fill(HIST("MCEventHist"), 1);
    histos.fill(HIST("hImpactParameterGen"), mcCollision.impactParameter());
    histos.fill(HIST("hNchMC_AllGen"), nChMC); // denominator for event loss

    if (atLeastOne) {
      histos.fill(HIST("MCEventHist"), 2);
      histos.fill(HIST("hImpactParameterRec"), mcCollision.impactParameter());
      histos.fill(HIST("hImpactParvsCentrRec"), centrality, mcCollision.impactParameter());
      histos.fill(HIST("hNchMC_WithRecoEvt"), nChMC);                 // numerator for event loss
      histos.fill(HIST("hNchVsCent_WithRecoEvt"), centrality, nChMC); // CO map
    }

    // ── Signal loss particle loop ────────────────────────────────────────────
    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }

      histos.fill(HIST("hPtVsNchMC_AllGen"), particle.pt(), nChMC); // denominator: all generated particles from all events
      histos.fill(HIST("hgenptBeforeEvtSelPbPb"), particle.pt(), mcCollision.impactParameter());

      if (atLeastOne) {
        histos.fill(HIST("hPtVsNchMC_WithRecoEvt"), particle.pt(), nChMC); // numerator: particles from events with ≥1 reco collision passing selection
        histos.fill(HIST("hgenptAfterEvtSelPbPb"), particle.pt(), mcCollision.impactParameter());
      }
    }
  }

  PROCESS_SWITCH(PtmultCorr, processDataPbPb, "process data heavy-ion", false);
  PROCESS_SWITCH(PtmultCorr, processDatapp, "process data pp", false);
  PROCESS_SWITCH(PtmultCorr, processMCeffPbPb, "process MC heavy-ion", false);
  PROCESS_SWITCH(PtmultCorr, processMCeffpp, "process MC pp", false);
  PROCESS_SWITCH(PtmultCorr, processEvtLossSigLossMCpp, "process Signal Loss, Event Loss", false);
  PROCESS_SWITCH(PtmultCorr, processEvtLossSigLossMCPbPb, "process Signal Loss, Event Loss", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PtmultCorr>(cfgc)};
}
