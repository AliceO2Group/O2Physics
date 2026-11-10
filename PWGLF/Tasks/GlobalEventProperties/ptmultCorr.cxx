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
const AxisSpec axisGenPtVary = {kGenpTend - 1, +kGenpTbegin + 0.5, +kGenpTend - 0.5, "", "GenpTVaryAxis"};
const AxisSpec axisGenTrkType = {kGenTrkTypeend - 1, +kGenTrkTypebegin + 0.5, +kGenTrkTypeend - 0.5, "", "GenTrackTypeAxis"};
const AxisSpec axisRecTrkType = {kRecTrkTypeend - 1, +kRecTrkTypebegin + 0.5, +kRecTrkTypeend - 0.5, "", "RecTrackTypeAxis"};
const AxisSpec axisTrackType = {kTrackTypeend - 1, +kTrackTypebegin + 0.5, +kTrackTypeend - 0.5, "", "TrackTypeAxis"};
auto static constexpr KminCharge = 3.f;
static constexpr float kZeroFloat{0.f};

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

  // ── pT-dependent DCAxy cut ──────────────────────────────────────────────────
  // The THnSparses keep the FULL DCAxy distribution (needed for the template fit).
  // When the flag is on, additional *DcaCut histograms (no DCAxy axis) are
  // filled only with tracks passing |DCAxy| < N * (A + B/pT^C), using reco pT.
  Configurable<bool> cfgApplyPtDepDcaXY{"cfgApplyPtDepDcaXY", false, "fill *DcaCut spectra with |DCAxy| < N*(A + B/pT^C)"};
  Configurable<float> cfgDcaXYSigmaA{"cfgDcaXYSigmaA", 0.0015f, "sigma_DCAxy(pT) = A + B/pT^C : A (cm)"};
  Configurable<float> cfgDcaXYSigmaB{"cfgDcaXYSigmaB", 0.0050f, "sigma_DCAxy(pT) = A + B/pT^C : B (cm (GeV/c)^C)"};
  Configurable<float> cfgDcaXYSigmaC{"cfgDcaXYSigmaC", 1.1f, "sigma_DCAxy(pT) = A + B/pT^C : C"};
  Configurable<float> cfgDcaXYNSigma{"cfgDcaXYNSigma", 7.f, "number of sigmas for the DCAxy cut"};

  // DCAxy binning of the THnSparses, by regions in |DCAxy| from 0 outwards:
  // region i spans (edge[i-1], edge[i]] with bin width width[i].
  // Default = previous uniform binning (20 um over +-5 cm).
  // Finer core, e.g.: edges {0.05, 0.5, 5.0}, widths {0.0005, 0.005, 0.1}
  Configurable<std::vector<double>> cfgDcaXYRegionEdges{"cfgDcaXYRegionEdges", {5.0}, "DCAxy binning: outer |DCAxy| edge of each region (cm), increasing"};
  Configurable<std::vector<double>> cfgDcaXYRegionWidths{"cfgDcaXYRegionWidths", {0.002}, "DCAxy binning: bin width in each region (cm)"};

  // Fine DCAxy vs reco pT of primaries (sigma(pT) extraction, MC only)
  ConfigurableAxis binsDcaXYFine{"binsDcaXYFine", {4000, -0.4, 0.4}, "fine DCAxy binning for sigma(pT) (cm)"};
  ConfigurableAxis binsPtDcaFine{"binsPtDcaFine", {VARIABLE_WIDTH, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.6, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 20.0, 30.0, 50.0}, "reco pT binning of the fine DCAxy histograms"};

  bool passDcaXY(float pt, float dcaXY) const
  {
    if (pt <= kZeroFloat) {
      return false;
    }
    const float sigma = cfgDcaXYSigmaA.value + cfgDcaXYSigmaB.value / std::pow(pt, cfgDcaXYSigmaC.value);
    return std::fabs(dcaXY) < cfgDcaXYNSigma.value * sigma;
  }

  static std::vector<double> makeDcaXYEdges(const std::vector<double>& regionEdges, const std::vector<double>& regionWidths)
  {
    if (regionEdges.empty() || regionEdges.size() != regionWidths.size()) {
      LOGF(fatal, "DCAxy binning: %zu region edges vs %zu widths", regionEdges.size(), regionWidths.size());
    }
    std::vector<double> pos{0.};
    double x0 = 0.;
    for (size_t i = 0; i < regionEdges.size(); ++i) {
      const double xMax = regionEdges[i];
      const double w = regionWidths[i];
      if (xMax <= x0 || w <= kZeroFloat) {
        LOGF(fatal, "DCAxy binning: region %zu invalid (edge %g after %g, width %g)", i, xMax, x0, w);
      }
      const int n = std::max(1, static_cast<int>(std::lround((xMax - x0) / w)));
      for (int j = 1; j <= n; ++j) {
        pos.push_back(x0 + j * (xMax - x0) / n); // integer steps: region edges are exact
      }
      if (std::fabs((xMax - x0) / n - w) > kWidthRelTol * w) {
        LOGF(warning, "DCAxy binning: region %zu width adjusted from %g to %g cm to fit the range", i, w, (xMax - x0) / n);
      }
      x0 = xMax;
    }
    std::vector<double> edges;
    for (auto it = pos.rbegin(); it != pos.rend(); ++it) {
      if (*it > kZeroFloat) {
        edges.push_back(-*it);
      }
    }
    edges.insert(edges.end(), pos.begin(), pos.end());
    LOGF(info, "DCAxy axis: %zu bins in [%g, %g] cm", edges.size() - 1, edges.front(), edges.back());
    return edges;
  }

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
  static constexpr double kWidthRelTol{1e-3};

  void init(InitContext const&)
  {
    if (requireGoodRct) {
      rctChecker.init(rctLabel.value, rctCheckZDC.value, rctTreatLimitedAcceptanceAsBad.value);
    }
    AxisSpec centAxis = {centralityBinning, "Centrality", "CentralityAxis"};
    AxisSpec axisPt = {ptHistBin, "pT", "pTAxis"};
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter"};
    const AxisSpec axisDCAxy{makeDcaXYEdges(cfgDcaXYRegionEdges.value, cfgDcaXYRegionWidths.value), "DCA_{xy} (cm)", "DCAxyAxis"};
    const AxisSpec axisDcaXYFine{binsDcaXYFine, "DCA_{xy} (cm)", "DCAxyFineAxis"};
    const AxisSpec axisPtDcaFine{binsPtDcaFine, "#it{p}_{T}^{rec} (GeV/#it{c})", "pTDcaFineAxis"};

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);
    histos.add("CentPercentileHist", "CentPercentileHist", kTH1D, {axisCent}, false);
    histos.add("PhiVsEtaHistNoCut", "PhiVsEtaHistNoCut", kTH2D, {axisPhi2, axisEta}, false);
    histos.add("PhiVsEtaHistWithCut", "PhiVsEtaHistWithCut", kTH2D, {axisPhi2, axisEta}, false);

    if (doprocessDataPbPb) {
      histos.add("hdatazvtxcent", "hdatazvtxcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hdatahistPbPb", "hdatahistPbPb", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisTrackType, axisDCAxy}, false);
      if (cfgApplyPtDepDcaXY) {
        histos.add("hdatahistPbPbDcaCut", "hdatahistPbPbDcaCut", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisTrackType}, false);
      }
    }

    if (doprocessDatapp) {
      histos.add("hdatahistpp", "hdatahistpp", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisTrackType, axisDCAxy}, false);
      if (cfgApplyPtDepDcaXY) {
        histos.add("hdatahistppDcaCut", "hdatahistppDcaCut", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisTrackType}, false);
      }
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
      if (cfgApplyPtDepDcaXY) {
        histos.add("hPbPbRecMCdndptDcaCut", "hPbPbRecMCdndptDcaCut", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisRecTrkType, axisTrackType}, false);
      }
      histos.add("hPbPbEtaReso", "hPbPbEtaReso", kTH2D, {axisPt, axisDeltaPt});
      // DCAxy of reco primaries (global tracks, same MC collision) vs RECO pT and centrality
      histos.add("hPbPbDcaXYPrimFine", "hPbPbDcaXYPrimFine", kTH3F, {axisPtDcaFine, axisDcaXYFine, centAxis}, false);
      // primaries whose MC particle belongs to another MC collision (mis-association)
      histos.add("hPbPbDcaXYPrimWrongColl", "hPbPbDcaXYPrimWrongColl", kTH2F, {axisPtDcaFine, axisDcaXYFine}, false);
    }

    if (doprocessMCeffpp) {
      histos.add("hppGenMCvtxz", "hppGenMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppGenMCAssoRecvtxz", "hppGenMCAssoRecvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppGenMCdndpt", "hppGenMCdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi}, false);
      // axisGenTrkType now has 6 bins: kGenAll(1)…kGenOther(5), kGenAllCharged(6)
      histos.add("hppGenMCAssoRecdndpt", "hppGenMCAssoRecdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisGenTrkType, axisGenPtVary}, false);

      histos.add("hppRecMCvtxz", "hppRecMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppRecMCdndpt", "hppRecMCdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisRecTrkType, axisTrackType, axisDCAxy}, false);
      if (cfgApplyPtDepDcaXY) {
        histos.add("hppRecMCdndptDcaCut", "hppRecMCdndptDcaCut", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisRecTrkType, axisTrackType}, false);
      }
      histos.add("hppEtaReso", "hppEtaReso", kTH2D, {axisPt, axisDeltaPt});
      histos.add("hppDcaXYPrimFine", "hppDcaXYPrimFine", kTH2F, {axisPtDcaFine, axisDcaXYFine}, false);
      histos.add("hppDcaXYPrimWrongColl", "hppDcaXYPrimWrongColl", kTH2F, {axisPtDcaFine, axisDcaXYFine}, false);
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

  // All track cuts EXCEPT the pT-dependent DCAxy cut (the flat DCAz cut is included).
  // The DCAxy cut is applied at the fill points via passDcaXY(), so that the
  // THnSparses keep the full DCAxy distribution for the template fit.
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
      auto trkType = track.hasTPC() ? kGlobalonly : kITSonly;
      // full DCAxy distribution (template fit)
      histos.fill(HIST("hdatahistPbPb"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), kGlobalplusITS, track.dcaXY());
      histos.fill(HIST("hdatahistPbPb"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), trkType, track.dcaXY());
      // spectrum with the pT-dependent DCAxy cut
      if (cfgApplyPtDepDcaXY && passDcaXY(track.pt(), track.dcaXY())) {
        histos.fill(HIST("hdatahistPbPbDcaCut"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), kGlobalplusITS);
        histos.fill(HIST("hdatahistPbPbDcaCut"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), trkType);
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
      auto trkType = track.hasTPC() ? kGlobalonly : kITSonly;
      // full DCAxy distribution (template fit)
      histos.fill(HIST("hdatahistpp"), cols.posZ(), track.pt(), track.phi(), kGlobalplusITS, track.dcaXY());
      histos.fill(HIST("hdatahistpp"), cols.posZ(), track.pt(), track.phi(), trkType, track.dcaXY());
      // spectrum with the pT-dependent DCAxy cut
      if (cfgApplyPtDepDcaXY && passDcaXY(track.pt(), track.dcaXY())) {
        histos.fill(HIST("hdatahistppDcaCut"), cols.posZ(), track.pt(), track.phi(), kGlobalplusITS);
        histos.fill(HIST("hdatahistppDcaCut"), cols.posZ(), track.pt(), track.phi(), trkType);
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
        if (!isTrackSelected(Rectrack)) { // all cuts except the pT-dependent DCAxy cut
          continue;
        }
        auto trkType = Rectrack.hasTPC() ? kGlobalonly : kITSonly;
        // DCAxy cut always evaluated with RECO pT (as in data)
        const bool passCut = cfgApplyPtDepDcaXY && passDcaXY(Rectrack.pt(), Rectrack.dcaXY());

        // fills the full-DCAxy THnSparse (templates) and, for tracks passing
        // the DCAxy cut, the *DcaCut copy (spectrum / efficiency numerator)
        auto fillRec = [&](float pt, float phi, int cls) {
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), pt, phi, static_cast<double>(cls), kGlobalplusITS, Rectrack.dcaXY());
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), pt, phi, static_cast<double>(cls), trkType, Rectrack.dcaXY());
          if (passCut) {
            histos.fill(HIST("hPbPbRecMCdndptDcaCut"), RecCol.posZ(), RecCol.centFT0C(), pt, phi, static_cast<double>(cls), kGlobalplusITS);
            histos.fill(HIST("hPbPbRecMCdndptDcaCut"), RecCol.posZ(), RecCol.centFT0C(), pt, phi, static_cast<double>(cls), trkType);
          }
        };

        fillRec(Rectrack.pt(), Rectrack.phi(), kRecoAll);

        if (Rectrack.has_mcParticle()) {
          int pid = 0;
          auto mcpart = Rectrack.mcParticle();
          histos.fill(HIST("hPbPbEtaReso"), Rectrack.pt(), Rectrack.pt() - mcpart.pt());
          fillRec(mcpart.pt(), mcpart.phi(), kRecoHasmc);
          if (mcpart.isPhysicalPrimary()) {
            // fine DCAxy for sigma(pT): global tracks, reco pT, no DCAxy cut
            if (Rectrack.hasTPC()) {
              if (mcpart.mcCollisionId() == mcCollision.globalIndex()) {
                histos.fill(HIST("hPbPbDcaXYPrimFine"), Rectrack.pt(), Rectrack.dcaXY(), RecCol.centFT0C());
              } else {
                histos.fill(HIST("hPbPbDcaXYPrimWrongColl"), Rectrack.pt(), Rectrack.dcaXY());
              }
            }
            fillRec(mcpart.pt(), mcpart.phi(), kRecoPrimary);
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
            fillRec(mcpart.pt(), mcpart.phi(), kRecoSecondary);
          }

          // duplicate-track (split) relabel — overrides origin, as before
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kRecoFake;
          }
          mclabels.push_back(Rectrack.mcParticleId());

          fillRec(mcpart.pt(), mcpart.phi(), pid);
        } else {
          fillRec(Rectrack.pt(), Rectrack.phi(), kRecoBkg);
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
        if (!isTrackSelected(Rectrack)) { // all cuts except the pT-dependent DCAxy cut
          continue;
        }
        auto trkType = Rectrack.hasTPC() ? kGlobalonly : kITSonly;
        // DCAxy cut always evaluated with RECO pT (as in data)
        const bool passCut = cfgApplyPtDepDcaXY && passDcaXY(Rectrack.pt(), Rectrack.dcaXY());

        auto fillRec = [&](float pt, float phi, int cls) {
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), pt, phi, static_cast<double>(cls), kGlobalplusITS, Rectrack.dcaXY());
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), pt, phi, static_cast<double>(cls), trkType, Rectrack.dcaXY());
          if (passCut) {
            histos.fill(HIST("hppRecMCdndptDcaCut"), RecCol.posZ(), pt, phi, static_cast<double>(cls), kGlobalplusITS);
            histos.fill(HIST("hppRecMCdndptDcaCut"), RecCol.posZ(), pt, phi, static_cast<double>(cls), trkType);
          }
        };

        fillRec(Rectrack.pt(), Rectrack.phi(), kRecoAll);

        if (Rectrack.has_mcParticle()) {
          int pid = 0;
          auto mcpart = Rectrack.mcParticle();
          histos.fill(HIST("hppEtaReso"), Rectrack.pt(), Rectrack.pt() - mcpart.pt());
          fillRec(mcpart.pt(), mcpart.phi(), kRecoHasmc);
          if (mcpart.isPhysicalPrimary()) {
            // fine DCAxy for sigma(pT): global tracks, reco pT, no DCAxy cut
            if (Rectrack.hasTPC()) {
              if (mcpart.mcCollisionId() == mcCollision.globalIndex()) {
                histos.fill(HIST("hppDcaXYPrimFine"), Rectrack.pt(), Rectrack.dcaXY());
              } else {
                histos.fill(HIST("hppDcaXYPrimWrongColl"), Rectrack.pt(), Rectrack.dcaXY());
              }
            }
            fillRec(mcpart.pt(), mcpart.phi(), kRecoPrimary);
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
            fillRec(mcpart.pt(), mcpart.phi(), kRecoSecondary);
          }
          // duplicate-track (split) relabel — overrides origin, as before
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kRecoFake;
          }
          mclabels.push_back(Rectrack.mcParticleId());
          fillRec(mcpart.pt(), mcpart.phi(), pid);
        } else {
          fillRec(Rectrack.pt(), Rectrack.phi(), kRecoBkg);
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
