// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file foxWolframCorrelation.cxx
/// \author Andreicovici Iulian Florin
/// \brief RECO-level event shape analysis using Fox-Wolfram Moments for O+O collisions.
///        Classifies events by topology (triangular), fills correlation histograms,
///        and produces custom SOA tables for mixed-event background estimation.
/// \since 10 Nov 2026

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TGraph.h>

#include <fmt/format.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::soa;

// Define compact event and track tables for mixed-event processing.
namespace o2::aod
{
namespace collision
{
DECLARE_SOA_COLUMN(Multip, multip, int);
}

DECLARE_SOA_TABLE(TriangleCollisions, "AOD", "TRIANGLECOLLISION",
                  o2::soa::Index<>,
                  collision::PosZ,
                  collision::Multip);

using TriangleCollision = TriangleCollisions::iterator;

namespace track_tri
{
DECLARE_SOA_INDEX_COLUMN(TriangleCollision, triangleCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace track_tri

DECLARE_SOA_TABLE(TriangleTracks, "AOD", "TRITRACK",
                  o2::soa::Index<>,
                  track_tri::TriangleCollisionId,
                  track_tri::Pt, track_tri::Eta, track_tri::Phi);

} // namespace o2::aod

// Wrap the azimuthal difference to the correlation interval.
double computeDeltaPhi(double phi1, double phi2)
{
  return RecoDecay::constrainAngle(phi1 - phi2, -o2::constants::math::PIHalf);
}

// Name the available pair-selection modes.
enum LeadAssocMode : int { // 3 distinct modes
  kAllPairs = 0,
  kpTLeadConstrain = 1,
  kpTLeadFree = 2
};

// Name the available event-topology selections.
enum TopologyMode : int {
  kAllEvSetNom = 0,
  kTriEvSetNom = 1,
  kIsoEvSetNom = 2,
  kTriEvSet1p = 3,
  kTriEvSet2p = 4,
  kTriEvSet3p = 5,
  kTriEvSet1m = 6,
  kTriEvSet2m = 7,
  kTriEvSet3m = 8
};

template <typename Track, typename TLead>
// Apply the same pair-selection rules to signal and background.
bool passesLeadAssoc(int mode, bool foundPart, Track const& track1, Track const& track2,
                     TLead const& triggTrk, TLead const& triggTrkPart)
{
  auto ptUpperLim1 = 3.f, ptUpperLim2 = 2.f; // low and high limits in pt

  if (mode == kAllPairs) { // all i-j pairs
    return (track1.pt() < ptUpperLim1) && (track2.pt() < ptUpperLim1);
  }
  if (mode == kpTLeadConstrain) { // L-A: Leading pt in (2,3), Assoc pt < 2
    return foundPart && (track1.globalIndex() == triggTrkPart.globalIndex()) && (track2.pt() < ptUpperLim2);
  }
  if (mode == kpTLeadFree) { // L-A
    return (track1.globalIndex() == triggTrk.globalIndex()) && (track2.pt() < track1.pt());
  }

  return false;
}

template <typename T>
// Find the overall leading track and the leading track in the selected momentum window.
auto findLeadingPart(T const& tracks)
{
  int leadingAllId = 0, leadingPartId = 0;
  int i = 0;

  double ptMaxAll = -1., ptMaxPart = -1.;
  bool foundPart = false;

  auto ptLim1 = 2.f, ptLim2 = 3.f; // low and high limits in pt

  for (const auto& track : tracks) {
    // pT Leading (no constraint)
    if (track.pt() > ptMaxAll) {
      ptMaxAll = track.pt();
      leadingAllId = i;
    }

    // pT Leading in (2,3)
    if (track.pt() > ptLim1 && track.pt() < ptLim2 && track.pt() > ptMaxPart) {
      ptMaxPart = track.pt();
      leadingPartId = i;
      foundPart = true;
    }
    i++;
  }

  return std::make_tuple(tracks.iteratorAt(leadingAllId), tracks.iteratorAt(leadingPartId), foundPart);
}

// Store selected tracks and collisions across data frames.
using TrackTuple = std::tuple<float, float, float>;
using CollisionTuple = std::tuple<float, int, std::vector<TrackTuple>>;

// Define the reconstructed collision and track views.
using RecCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
using CollisionRecTable = RecCollisions::iterator;
using RecTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

// Select events, compute moments, and fill same-event observables.
struct FoxWolframCorrelation {

  SliceCache cache;
  Preslice<RecTracks> perCollision = aod::track::collisionId;

  Produces<aod::TriangleCollisions> triCollisions;
  Produces<aod::TriangleTracks> triTracks;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  //// Separate the 120 folded spectra to stay below the registry capacity ////
  HistogramRegistry foldedHistos{"foldedHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // Configure histogram axes, event selection, and topology.
  Configurable<float> vtxRange{"vtxRange", 10.0f, "Vertex Z range to consider"};
  Configurable<float> etaRange{"etaRange", 0.8f, "eta range to consider"};
  Configurable<int> nSelEv{"nSelEv", 200, "no. of mixing Ev."};
  Configurable<bool> doSameEvent{"doSameEvent", true, "Fill same-event TRI correlations"};
  Configurable<bool> doFWMCorrelations{"doFWMCorrelations", false, "Fill H_i-H_j corr. hist"};
  Configurable<int> topologyMode{"topologyMode", kTriEvSetNom,
                                 "SelOpt: 0 AllEv, 1 TriNom, 2 IsoNom, 3 TriSet1p, 4 TriSet2p, 5 TriSet3p, 6 TriSet1m, 7 TriSet2m, 8 TriSet3m"};

  //// Minimum associated pT included in mean-pT profiles ////
  Configurable<double> meanPtMin{"meanPtMin", 0.5, "Associated pT threshold for mean-pT profiles (GeV/c)"};

  // Choose the pair selection for same-event correlations.
  Configurable<int> leadAssocMode{"leadAssocMode", 2, "0: ALL i-j pairs (pT<3); 1: pTLead in (2,3) & pTAssoc < 2; 2: pTAssoc < pTLead"};

  Configurable<bool> isVtxRange{"isVtxRange", true, "isVtxRange"};
  Configurable<bool> isSel8{"isSel8", true, "isSel8"};
  Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", true, "isNoSameBunchPileup"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", true, "isGoodZvtxFT0vsPV"};
  Configurable<bool> isVertexITSTPC{"isVertexITSTPC", false, "isVertexITSTPC"};
  Configurable<bool> isVertexTOFmatched{"isVertexTOFmatched", false, "isVertexTOFmatched"};
  Configurable<bool> kIsGoodITSLayer0123{"kIsGoodITSLayer0123", false, "kIsGoodITSLayer0123"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", false, "isNoCollInTimeRangeNarrow"};

  //// Configurable FWM limits: H1, H2, H3, H4, H5, H6, H7, H8 ////
  Configurable<std::vector<double>> kTriFWMSetNom{"kTriFWMSetNom", {0.05, 0.30, 0.040, 0.190, 0.020, 0.150, 0.016, 0.120}, "kTriFWMSetNom"};

  Configurable<std::vector<double>> kTriFWMSet1p{"kTriFWMSet1p", {0.0525, 0.315, 0.0380, 0.1995, 0.0190, 0.1575, 0.0152, 0.1260}, "kTriFWMSet1p"};
  Configurable<std::vector<double>> kTriFWMSet2p{"kTriFWMSet2p", {0.0550, 0.330, 0.0360, 0.2090, 0.0180, 0.1650, 0.0144, 0.1320}, "kTriFWMSet2p"};
  Configurable<std::vector<double>> kTriFWMSet3p{"kTriFWMSet3p", {0.0575, 0.345, 0.0340, 0.2185, 0.0170, 0.1725, 0.0136, 0.1380}, "kTriFWMSet3p"};

  Configurable<std::vector<double>> kTriFWMSet1m{"kTriFWMSet1m", {0.0475, 0.2850, 0.0420, 0.1805, 0.0210, 0.1425, 0.0168, 0.1140}, "kTriFWMSet1m"};
  Configurable<std::vector<double>> kTriFWMSet2m{"kTriFWMSet2m", {0.0450, 0.2700, 0.0440, 0.1710, 0.0220, 0.1350, 0.0176, 0.1080}, "kTriFWMSet2m"};
  Configurable<std::vector<double>> kTriFWMSet3m{"kTriFWMSet3m", {0.0425, 0.2550, 0.0460, 0.1615, 0.0230, 0.1275, 0.0184, 0.1020}, "kTriFWMSet3m"};

  //// Named TRI selection indices ////
  enum TriSelection {
    kTriNominal,
    kTriSet1p,
    kTriSet2p,
    kTriSet3p,
    kTriSet1m,
    kTriSet2m,
    kTriSet3m,
    kNTriSel
  };
  using TriFlags = std::array<bool, kNTriSel>;

  ConfigurableAxis cfgAxisPt{"cfgAxisPt", {VARIABLE_WIDTH, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0}, "p_{T} Bins"};
  ConfigurableAxis cfgAxisMultip{"cfgAxisMultip", {351, -0.5, 350.5}, "Multiplicity Bins"};

  Configurable<std::vector<double>> isoFWMnominal{"isoFWMnominal", {0.02, 0.27, 0.02, 0.16, 0.02, 0.12, 0.02, 0.10}, "isoFWMnominal"};

  // Multiplicity bins: bin i is [multLower[i], multUpper[i]] (inclusive limits).
  // The first three entries are the broad ranges (RANGE1-3: FWM_BIN1-3, lead-assoc, mixing buffers),
  // the remaining ones are the finer intervals (BIN4...). Bins may overlap. Max 10 bins.
  // The same limits must be given to FoxWolframCorrelationMixing.
  Configurable<std::vector<int>> multLower{"multLower", {40, 86, 141, 40, 61, 81, 101, 121, 141, 161}, "lower multiplicity limit of each bin (inclusive)"};
  Configurable<std::vector<int>> multUpper{"multUpper", {85, 140, 350, 60, 80, 100, 120, 140, 160, 350}, "upper multiplicity limit of each bin (inclusive)"};

  static constexpr int MaxMultBins = 10;
  int nMultBins = 0; // set in init() from multLower/multUpper

  static constexpr std::array<std::string_view, 10> BinNames = {
    "BIN1", "BIN2", "BIN3", "BIN4", "BIN5",
    "BIN6", "BIN7", "BIN8", "BIN9", "BIN10"};

  //// Class-level names for compile-time HIST lookups in nested loops ////
  static constexpr std::array<std::string_view, 6> NearNames = {
    "nearPtD0_", "nearPtD1_", "nearPtD2_",
    "nearPtD3_", "nearPtD4_", "nearPtD5_"};

  static constexpr std::array<std::string_view, 6> AwayNames = {
    "awayPtD0_", "awayPtD1_", "awayPtD2_",
    "awayPtD3_", "awayPtD4_", "awayPtD5_"};

  //// Profile names match the folded LR declarations ////
  static constexpr std::array<std::string_view, 6> NearMeanNames = {
    "nearPtDMean0_", "nearPtDMean1_", "nearPtDMean2_",
    "nearPtDMean3_", "nearPtDMean4_", "nearPtDMean5_"};

  static constexpr std::array<std::string_view, 6> AwayMeanNames = {
    "awayPtDMean0_", "awayPtDMean1_", "awayPtDMean2_",
    "awayPtDMean3_", "awayPtDMean4_", "awayPtDMean5_"};

  static constexpr std::array<std::string_view, 3> AllLrNames = {
    "allLRNear", "allLRValley", "allLRAway"};

  static constexpr std::array<std::string_view, 3> AllSrNames = {
    "allSRNear", "allSRValley", "allSRAway"};

  static constexpr std::array<std::string_view, 4> TriLrNames = {
    "triLRPeak1", "triLRPeak2", "triLRPeak3", "triLRValleys"};

  static constexpr std::array<std::string_view, 4> TriSrNames = {
    "triSRPeak1", "triSRPeak2", "triSRPeak3", "triSRValleys"};

  // Map each unique pair of Fox-Wolfram moments to its histogram.
  static constexpr std::array<std::string_view, 3> FwmRangeNames = {"FWM_BIN1/", "FWM_BIN2/", "FWM_BIN3/"};
  static constexpr std::array<std::string_view, 8> FwmHlNames = {"H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8"};
  static constexpr std::array<std::string_view, 8> FwmMulthNames = {"Multip_H1", "Multip_H2", "Multip_H3", "Multip_H4", "Multip_H5", "Multip_H6", "Multip_H7", "Multip_H8"};

  static constexpr std::array<std::string_view, 28> FwmCorrelationNames = {"H1_H2", "H1_H3", "H1_H4", "H1_H5", "H1_H6", "H1_H7", "H1_H8", "H2_H3", "H2_H4", "H2_H5", "H2_H6", "H2_H7", "H2_H8", "H3_H4", "H3_H5", "H3_H6", "H3_H7", "H3_H8", "H4_H5", "H4_H6", "H4_H7", "H4_H8", "H5_H6", "H5_H7", "H5_H8", "H6_H7", "H6_H8", "H7_H8"};
  static constexpr std::array<int, 28> FwmCorrelationFirst = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6};
  static constexpr std::array<int, 28> FwmCorrelationSecond = {1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7};

  // Register event, track, moment, and correlation histograms.
  void init(InitContext const&)
  {
    // Reject unsupported topology modes at startup.
    if (topologyMode.value < kAllEvSetNom || topologyMode.value > kTriEvSet3m) {
      LOGF(fatal, "Invalid topologyMode=%d. Allowed values: 0 ALL, 1 TRI nominal, 2 ISO, 3 TRI 1p, 4 TRI 2p, 5 TRI 3p, 6 TRI 1m, 7 TRI 2m, 8 TRI 3m", topologyMode.value);
    }

    // Validate the configurable multiplicity bins.
    if (multLower.value.size() != multUpper.value.size() || multLower.value.empty() || multLower.value.size() > static_cast<size_t>(MaxMultBins)) {
      LOGF(fatal, "multLower/multUpper must have the same size, between 1 and %d entries", MaxMultBins);
    }
    nMultBins = static_cast<int>(multLower.value.size());

    for (int i = 0; i < nMultBins; ++i) {
      if (multLower.value[i] > multUpper.value[i]) {
        LOGF(fatal, "Multiplicity bin %d has multLower (%d) > multUpper (%d)", i + 1, multLower.value[i], multUpper.value[i]);
      }
    }

    // Define axes shared by the registered histograms.
    const AxisSpec axispT{cfgAxisPt, "p_{T} (GeV/c)"};
    const AxisSpec axisMultip{cfgAxisMultip, "Multiplicity"};
    const AxisSpec axisDeltaEta{32, -1.6, +1.6, "#Delta#eta"};
    const AxisSpec axisDeltaPhi{72, -constants::math::PIHalf, +3 * constants::math::PIHalf, "#Delta#phi"};
    const AxisSpec axisPhi{100, 0, constants::math::TwoPI, "#phi"};
    const AxisSpec axisEta{40, -0.8, 0.8, "#eta"};
    const AxisSpec axisHlrec{50, 0, 1, "H_{l}"};
    const AxisSpec axisHlfine{100, 0, 1, "H_{l}"};
    const AxisSpec axisEvent{10, 0.5, 10.5, "", "EventAxis"};
    const AxisSpec axisRecoil{100, 0., 1, "Recoil"};

    // Label same-event pair counts by multiplicity interval.
    histos.add("SamePairs/SamePairsNoLA", "SamePairsNoLA", kTH1D, {axisEvent}, false);
    auto* samePairsNoLA = histos.get<TH1>(HIST("SamePairs/SamePairsNoLA"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      samePairsNoLA->SetBinLabel(i + 1, i < 3 ? fmt::format("samePairsLA_RANGE{}", i + 1).c_str() : fmt::format("samePairsLA_BIN{}", i + 1).c_str());
    }

    // ALL Ev per Sel. BIN
    histos.add("Events_Details/allEvSetNom", "SetAllNominal", kTH1D, {axisEvent}, false);
    auto* allEvSetNom = histos.get<TH1>(HIST("Events_Details/allEvSetNom"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      allEvSetNom->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // ISO Ev per Sel. BIN
    histos.add("Events_Details/eventsNoISO", "SetIsoNominal", kTH1D, {axisEvent}, false);
    auto* eventsNoISO = histos.get<TH1>(HIST("Events_Details/eventsNoISO"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      eventsNoISO->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // TRI Ev per Sel. BIN (NOMINAL)
    histos.add("Events_Details/triEvSetNom", "SetTriNominal", kTH1D, {axisEvent}, false);
    auto* triEvSetNom = histos.get<TH1>(HIST("Events_Details/triEvSetNom"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      triEvSetNom->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // TRI Ev per Sel. BIN (SET1PLUS)
    histos.add("Events_Details/triEvSet1p", "Set1TriPlus", kTH1D, {axisEvent}, false);
    auto* triEvSet1p = histos.get<TH1>(HIST("Events_Details/triEvSet1p"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      triEvSet1p->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // TRI Ev per Sel. BIN (SET2PLUS)
    histos.add("Events_Details/triEvSet2p", "Set2TriPlus", kTH1D, {axisEvent}, false);
    auto* triEvSet2p = histos.get<TH1>(HIST("Events_Details/triEvSet2p"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      triEvSet2p->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // TRI Ev per Sel. BIN (SET3PLUS)
    histos.add("Events_Details/triEvSet3p", "Set3TriPlus", kTH1D, {axisEvent}, false);
    auto* triEvSet3p = histos.get<TH1>(HIST("Events_Details/triEvSet3p"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      triEvSet3p->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // TRI Ev per Sel. BIN (SET1MINUS)
    histos.add("Events_Details/triEvSet1m", "Set1TriMinus", kTH1D, {axisEvent}, false);
    auto* triEvSet1m = histos.get<TH1>(HIST("Events_Details/triEvSet1m"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      triEvSet1m->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // TRI Ev per Sel. BIN (SET2MINUS)
    histos.add("Events_Details/triEvSet2m", "Set2TriMinus", kTH1D, {axisEvent}, false);
    auto* triEvSet2m = histos.get<TH1>(HIST("Events_Details/triEvSet2m"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      triEvSet2m->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // TRI Ev per Sel. BIN (SET3MINUS)
    histos.add("Events_Details/triEvSet3m", "Set3TriMinus", kTH1D, {axisEvent}, false);
    auto* triEvSet3m = histos.get<TH1>(HIST("Events_Details/triEvSet3m"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      triEvSet3m->SetBinLabel(i + 1, i < 3 ? fmt::format("RANGE{}", i + 1).c_str() : fmt::format("BIN{}", i + 1).c_str());
    }

    // Register inclusive multiplicity, track spectra, and recoil histograms.
    histos.add("Multip/multipAllEv", "multipAllEv; trks; Events", kTH1D, {axisMultip});
    histos.add("Multip/multipTriEv", "multipTriEv; trks; Events", kTH1D, {axisMultip});

    histos.add("Multip/nTrkTriEvS1p", "nTrkTriEvS1p", kTH1D, {axisMultip});
    histos.add("Multip/nTrkTriEvS2p", "nTrkTriEvS2p", kTH1D, {axisMultip});
    histos.add("Multip/nTrkTriEvS3p", "nTrkTriEvS3p", kTH1D, {axisMultip});
    histos.add("Multip/nTrkTriEvS1m", "nTrkTriEvS1m", kTH1D, {axisMultip});
    histos.add("Multip/nTrkTriEvS2m", "nTrkTriEvS2m", kTH1D, {axisMultip});
    histos.add("Multip/nTrkTriEvS3m", "nTrkTriEvS3m", kTH1D, {axisMultip});

    // Assoc. pT spectra
    histos.add("pTSpectra/pTAssocAllEv", "pTAssocAllEv", kTH1D, {axispT});
    histos.add("pTSpectra/pTAssocTriEv", "pTAssocTriEv", kTH1D, {axispT});
    histos.add("pTSpectra/meanpTAssocMultipAllEv", "meanpTAssocMultipAllEv", kTProfile, {axisMultip});
    histos.add("pTSpectra/meanpTAssocMultipTriEv", "meanpTAssocMultipTriEv", kTProfile, {axisMultip});

    // Lead. pT spectra
    histos.add("pTSpectra/pTLeadAllEv", "pTLeadAllEv", kTH1D, {axispT});
    histos.add("pTSpectra/pTLeadTriEv", "pTLeadTriEv", kTH1D, {axispT});
    histos.add("pTSpectra/ptLeadMultipAllEv", "ptLeadMultipAllEv", kTH2D, {axisMultip, axispT});
    histos.add("pTSpectra/ptLeadMultipTriEv", "ptLeadMultipTriEv", kTH2D, {axisMultip, axispT});
    histos.add("pTSpectra/meanptLeadMultipAllEv", "meanptLeadMultipAllEv", kTProfile, {axisMultip});
    histos.add("pTSpectra/meanptLeadMultipTriEv", "meanptLeadMultipTriEv", kTProfile, {axisMultip});

    histos.add("pTSpectra/RecoilAllEv", "RecoilAllEv", kTH2D, {axisMultip, axisRecoil});
    histos.add("pTSpectra/RecoilTriEv", "RecoilTriEv", kTH2D, {axisMultip, axisRecoil});

    // Register track spectra for each multiplicity interval.
    for (int i = 1; i <= nMultBins; i++) {

      // eta-phi and pt-phi
      histos.add(fmt::format("Spectra/etaphiAllEv_BIN{}", i).c_str(), fmt::format("etaphiAllEv_BIN{}", i).c_str(), kTH2D, {axisPhi, axisEta});
      histos.add(fmt::format("Spectra/ptphiAllEv_BIN{}", i).c_str(), fmt::format("ptphiAllEv_BIN{}", i).c_str(), kTH2D, {axisPhi, axispT});

      histos.add(fmt::format("Spectra/etaphiTriEv_BIN{}", i).c_str(), fmt::format("etaphiTriEv_BIN{}", i).c_str(), kTH2D, {axisPhi, axisEta});
      histos.add(fmt::format("Spectra/ptphiTriEv_BIN{}", i).c_str(), fmt::format("ptphiTriEv_BIN{}", i).c_str(), kTH2D, {axisPhi, axispT});

      //// Folded LR pT spectra and mean pT versus multiplicity ////
      for (int AngBin = 0; AngBin < 6; ++AngBin) {

        foldedHistos.add(fmt::format("triEvFoldedLR/nearPtD{}_BIN{}", AngBin, i).c_str(), fmt::format("LR near: {}-{} deg from 0", 10 * AngBin, 10 * (AngBin + 1)).c_str(), kTH1D, {axispT});
        foldedHistos.add(fmt::format("triEvFoldedLR/nearPtDMean{}_BIN{}", AngBin, i).c_str(), fmt::format("LR near: {}-{} deg from 0", 10 * AngBin, 10 * (AngBin + 1)).c_str(), kTProfile, {axisMultip});

        foldedHistos.add(fmt::format("triEvFoldedLR/awayPtD{}_BIN{}", AngBin, i).c_str(), fmt::format("LR away: {}-{} deg from 120 and 240", 10 * AngBin, 10 * (AngBin + 1)).c_str(), kTH1D, {axispT});
        foldedHistos.add(fmt::format("triEvFoldedLR/awayPtDMean{}_BIN{}", AngBin, i).c_str(), fmt::format("LR away: {}-{} deg from 120 and 240", 10 * AngBin, 10 * (AngBin + 1)).c_str(), kTProfile, {axisMultip});
      }

      //// ALL: pT spectra (TH1D) and mean associated pT versus multiplicity (TProfile) ////
      histos.add(fmt::format("allEvPkVall/allLRNearPt_BIN{}", i).c_str(), fmt::format("allLRNearPt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("allEvPkVall/allLRNearPtMultip_BIN{}", i).c_str(), fmt::format("allLRNearPtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});

      histos.add(fmt::format("allEvPkVall/allSRNearPt_BIN{}", i).c_str(), fmt::format("allSRNearPt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("allEvPkVall/allSRNearPtMultip_BIN{}", i).c_str(), fmt::format("allSRNearPtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});

      histos.add(fmt::format("allEvPkVall/allLRValleyPt_BIN{}", i).c_str(), fmt::format("allLRValleyPt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("allEvPkVall/allLRValleyPtMultip_BIN{}", i).c_str(), fmt::format("allLRValleyPtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});
      histos.add(fmt::format("allEvPkVall/allSRValleyPt_BIN{}", i).c_str(), fmt::format("allSRValleyPt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("allEvPkVall/allSRValleyPtMultip_BIN{}", i).c_str(), fmt::format("allSRValleyPtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});

      histos.add(fmt::format("allEvPkVall/allLRAwayPt_BIN{}", i).c_str(), fmt::format("allLRAwayPt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("allEvPkVall/allLRAwayPtMultip_BIN{}", i).c_str(), fmt::format("allLRAwayPtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});
      histos.add(fmt::format("allEvPkVall/allSRAwayPt_BIN{}", i).c_str(), fmt::format("allSRAwayPt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("allEvPkVall/allSRAwayPtMultip_BIN{}", i).c_str(), fmt::format("allSRAwayPtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});

      //// TRI: three peaks and two valleys combined into one spectrum/profile ////
      histos.add(fmt::format("triEvPkVall/triLRPeak1Pt_BIN{}", i).c_str(), fmt::format("triLRPeak1Pt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("triEvPkVall/triLRPeak1PtMultip_BIN{}", i).c_str(), fmt::format("triLRPeak1PtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});
      histos.add(fmt::format("triEvPkVall/triSRPeak1Pt_BIN{}", i).c_str(), fmt::format("triSRPeak1Pt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("triEvPkVall/triSRPeak1PtMultip_BIN{}", i).c_str(), fmt::format("triSRPeak1PtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});

      histos.add(fmt::format("triEvPkVall/triLRPeak2Pt_BIN{}", i).c_str(), fmt::format("triLRPeak2Pt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("triEvPkVall/triLRPeak2PtMultip_BIN{}", i).c_str(), fmt::format("triLRPeak2PtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});
      histos.add(fmt::format("triEvPkVall/triSRPeak2Pt_BIN{}", i).c_str(), fmt::format("triSRPeak2Pt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("triEvPkVall/triSRPeak2PtMultip_BIN{}", i).c_str(), fmt::format("triSRPeak2PtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});

      histos.add(fmt::format("triEvPkVall/triLRPeak3Pt_BIN{}", i).c_str(), fmt::format("triLRPeak3Pt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("triEvPkVall/triLRPeak3PtMultip_BIN{}", i).c_str(), fmt::format("triLRPeak3PtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});
      histos.add(fmt::format("triEvPkVall/triSRPeak3Pt_BIN{}", i).c_str(), fmt::format("triSRPeak3Pt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("triEvPkVall/triSRPeak3PtMultip_BIN{}", i).c_str(), fmt::format("triSRPeak3PtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});

      histos.add(fmt::format("triEvPkVall/triLRValleysPt_BIN{}", i).c_str(), fmt::format("triLRValleysPt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("triEvPkVall/triLRValleysPtMultip_BIN{}", i).c_str(), fmt::format("triLRValleysPtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});
      histos.add(fmt::format("triEvPkVall/triSRValleysPt_BIN{}", i).c_str(), fmt::format("triSRValleysPt_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("triEvPkVall/triSRValleysPtMultip_BIN{}", i).c_str(), fmt::format("triSRValleysPtMultip_BIN{}", i).c_str(), kTProfile, {axisMultip});

      // Same Even Fct
      histos.add(fmt::format("Corelations/sTRI_BIN{}", i).c_str(), fmt::format("sTRI_BIN{}", i).c_str(), kTH2D, {axisDeltaPhi, axisDeltaEta});
      histos.add(fmt::format("Corelations/sProjTRI_BIN{}", i).c_str(), fmt::format("sProjTRI_BIN{}", i).c_str(), kTH1D, {axisDeltaPhi});
    }

    // Register moment distributions in the broad ranges.
    const std::array<int, 3> bins = {1, 2, 3};
    for (int j = 1; j <= 8; j++) {
      for (const auto& i : bins) {
        histos.add(fmt::format("FWM_BIN{}/H{}", i, j).c_str(), fmt::format("H_{{{}}}^{{rec}}", j).c_str(), kTH1F, {axisHlrec});
        histos.add(fmt::format("FWM_BIN{}/Multip_H{}", i, j).c_str(), fmt::format("H_{{{}}} vs Multip", j).c_str(), kTH2D, {axisMultip, axisHlfine});
      }
    }

    // Register optional pairwise moment correlations.
    if (doFWMCorrelations) {

      for (size_t k = 0; k < FwmCorrelationNames.size(); ++k) {
        const int first = FwmCorrelationFirst[k] + 1;
        const int second = FwmCorrelationSecond[k] + 1;

        histos.add(fmt::format("FWM_BIN1/{}", FwmCorrelationNames[k]).c_str(), fmt::format("H_{{{}}} vs H_{{{}}}", first, second).c_str(), kTH2D, {axisHlfine, axisHlfine});
        histos.add(fmt::format("FWM_BIN2/{}", FwmCorrelationNames[k]).c_str(), fmt::format("H_{{{}}} vs H_{{{}}}", first, second).c_str(), kTH2D, {axisHlfine, axisHlfine});
        histos.add(fmt::format("FWM_BIN3/{}", FwmCorrelationNames[k]).c_str(), fmt::format("H_{{{}}} vs H_{{{}}}", first, second).c_str(), kTH2D, {axisHlfine, axisHlfine});
      }
    }
  }

  // Apply reconstructed collision quality requirements.
  template <typename TCollision>
  bool isEventSelected(const TCollision& RecCollision)
  {
    if (isVtxRange && std::abs(RecCollision.posZ()) >= vtxRange) {
      return false;
    }

    if (isSel8 && !RecCollision.sel8()) {
      return false;
    }

    if (isNoSameBunchPileup && !RecCollision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (isGoodZvtxFT0vsPV && !RecCollision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (isVertexITSTPC && !RecCollision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }

    if (isVertexTOFmatched && !RecCollision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }

    if (kIsGoodITSLayer0123 && !RecCollision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      return false;
    }

    if (isNoCollInTimeRangeNarrow && !RecCollision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }

    return true;
  }

  // Compute normalized Fox-Wolfram moments through order eight.
  template <typename T>
  std::array<double, 9> computeFWM(T const& tracks)
  {
    std::array<double, 9> hSN{};
    std::array<double, 9> hS{};

    // Accumulate transverse-momentum-weighted contributions from track pairs.
    for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {

      double c = std::cos(t1.phi() - t2.phi());
      double w = t1.pt() * t2.pt();
      double p0 = 1.0, p1 = c;

      hSN[0] += w;
      hSN[1] += w * p1;

      // Evaluate higher Legendre orders with the recurrence relation.
      for (int l = 2; l < 9; ++l) {

        double p = ((2 * l - 1) * c * p1 - (l - 1) * p0) / l;

        hSN[l] += w * p;
        p0 = p1;
        p1 = p;
      }
    }

    // Normalize each moment by the zeroth moment.
    for (int l = 0; l < 9; ++l) {
      hS[l] = hSN[l] / hSN[0];
    }

    return hS;
  }

  // Calculate normalized transverse recoil.
  template <typename Tracks>
  double computeRecoil(const Tracks& tracks)
  {
    double sumPx = 0.0;
    double sumPy = 0.0;
    double sumPt = 0.0;

    // Sum track momenta to form the transverse recoil.
    for (const auto& track : tracks) {
      sumPx += track.px();
      sumPy += track.py();
      sumPt += track.pt();
    }

    if (sumPt == 0.0) {
      return 0.0;
    }

    return std::sqrt(sumPx * sumPx + sumPy * sumPy) / sumPt;
  }

  // Fill inclusive and topology-selected observables in one multiplicity interval.
  template <typename T>
  void fillHistoMult(int binIndex, int minMult, int maxMult,
                     bool isSelectedTopology, const TriFlags& triFlags, bool isIsoNominal,
                     const T& tracks)
  {
    if (tracks.size() < minMult || tracks.size() > maxMult) {
      return;
    }

    // Count all events and each topology selection in this interval.
    histos.fill(HIST("Events_Details/allEvSetNom"), binIndex);

    if (isIsoNominal) {
      histos.fill(HIST("Events_Details/eventsNoISO"), binIndex);
    }

    if (triFlags[kTriNominal]) {
      histos.fill(HIST("Events_Details/triEvSetNom"), binIndex);
    }

    // The remaining filters
    if (triFlags[kTriSet1p]) {
      histos.fill(HIST("Events_Details/triEvSet1p"), binIndex);
    }

    if (triFlags[kTriSet2p]) {
      histos.fill(HIST("Events_Details/triEvSet2p"), binIndex);
    }

    if (triFlags[kTriSet3p]) {
      histos.fill(HIST("Events_Details/triEvSet3p"), binIndex);
    }

    if (triFlags[kTriSet1m]) {
      histos.fill(HIST("Events_Details/triEvSet1m"), binIndex);
    }

    if (triFlags[kTriSet2m]) {
      histos.fill(HIST("Events_Details/triEvSet2m"), binIndex);
    }

    if (triFlags[kTriSet3m]) {
      histos.fill(HIST("Events_Details/triEvSet3m"), binIndex);
    }

    // Fill inclusive spectra and selected-topology spectra for this interval.
    static_for<0, 9>([&](auto i) {
      constexpr int IdxVal = i.value;

      if (binIndex == IdxVal + 1) {

        for (const auto& track : tracks) {

          histos.fill(HIST("Spectra/etaphiAllEv_") + HIST(BinNames[IdxVal]), track.phi(), track.eta());
          histos.fill(HIST("Spectra/ptphiAllEv_") + HIST(BinNames[IdxVal]), track.phi(), track.pt());

          if (!isSelectedTopology) {
            continue;
          } // TRI Ev.

          histos.fill(HIST("Spectra/etaphiTriEv_") + HIST(BinNames[IdxVal]), track.phi(), track.eta());
          histos.fill(HIST("Spectra/ptphiTriEv_") + HIST(BinNames[IdxVal]), track.phi(), track.pt());
        }
      }
    }); // per Multip Bin
  }

  //// Helper for one multiplicity bin (IdxVal is a template parameter, not a captured lambda variable) ////
  //// Avoids the GCC 14 ICE in tsubst_expr caused by nested generic static_for lambdas ////
  template <int IdxVal>
  void fillFoldedBin(int nTracks, double deltaPhi, double pt)
  {
    const double degToRad = constants::math::PI / 180.;

    //// Six angular intervals with the existing strict boundary cuts ////
    static_for<0, 5>([&](auto j) {
      constexpr int AngBin = j.value;
      constexpr double LowLim = 10. * AngBin;
      constexpr double HighLim = 10. * (AngBin + 1);

      //// Combine the two sides of the near peak at zero ////
      const bool inNear = (deltaPhi > -HighLim * degToRad && deltaPhi < -LowLim * degToRad) || (deltaPhi > LowLim * degToRad && deltaPhi < HighLim * degToRad);

      if (inNear) {
        if (pt > meanPtMin.value) {
          foldedHistos.fill(HIST("triEvFoldedLR/") + HIST(NearNames[AngBin]) + HIST(BinNames[IdxVal]), pt);
          foldedHistos.fill(HIST("triEvFoldedLR/") + HIST(NearMeanNames[AngBin]) + HIST(BinNames[IdxVal]), nTracks, pt);
        } // for pt > meanPtMin
      } // PEAK 1

      //// Combine both sides of 120 and 240 degrees ////
      //// The final window maps angles above 270 degrees to negative deltaPhi ////
      const bool inAway = (deltaPhi > (120. - HighLim) * degToRad && deltaPhi < (120. - LowLim) * degToRad) || (deltaPhi > (120. + LowLim) * degToRad && deltaPhi < (120. + HighLim) * degToRad) ||
                          (deltaPhi > (240. - HighLim) * degToRad && deltaPhi < (240. - LowLim) * degToRad) || (deltaPhi > (240. + LowLim) * degToRad && deltaPhi < (240. + HighLim) * degToRad) ||
                          (deltaPhi > (240. + LowLim - 360.) * degToRad && deltaPhi < (240. + HighLim - 360.) * degToRad);

      if (inAway) {
        if (pt > meanPtMin.value) {
          foldedHistos.fill(HIST("triEvFoldedLR/") + HIST(AwayNames[AngBin]) + HIST(BinNames[IdxVal]), pt);
          foldedHistos.fill(HIST("triEvFoldedLR/") + HIST(AwayMeanNames[AngBin]) + HIST(BinNames[IdxVal]), nTracks, pt);
        } // for pt > meanPtMin
      } // PEAKS 2 AND 3
    });
  }

  //// TRI LR pT spectra in six symmetric angular intervals ////
  template <typename Track>
  void fillSameEventFoldedLR(int binIndex,
                             int nTracks,
                             Track const& associatedTrack,
                             double deltaPhi,
                             double deltaEta,
                             bool isLeadAssocPair)
  {
    //// Event topology is already selected by the caller ////
    const bool isAnyTriMode = topologyMode.value == kTriEvSetNom ||
                              topologyMode.value == kTriEvSet1p ||
                              topologyMode.value == kTriEvSet2p ||
                              topologyMode.value == kTriEvSet3p ||
                              topologyMode.value == kTriEvSet1m ||
                              topologyMode.value == kTriEvSet2m ||
                              topologyMode.value == kTriEvSet3m;

    if (!isAnyTriMode || leadAssocMode.value != kpTLeadFree || !isLeadAssocPair) {
      return;
    }

    //// Keep only 1.25 < |deltaEta| < 1.60 ////
    const double absEta = std::abs(deltaEta);
    if (absEta <= 1.25 || absEta >= 1.60) {
      return;
    }

    //// Fill only the requested multiplicity interval ////
    static_for<0, 9>([&](auto i) {
      constexpr int IdxVal = i.value;

      if (binIndex != IdxVal + 1) {
        return;
      }

      fillFoldedBin<IdxVal>(nTracks, deltaPhi, associatedTrack.pt());
    });
  }

  //// Helper for one multiplicity bin (IdxVal is a template parameter, not a captured lambda variable) ////
  template <int IdxVal>
  void fillPtBin(bool isAll, bool isTri,
                 const std::array<bool, 3>& allWindows, const std::array<bool, 4>& triWindows,
                 bool inLR, bool inSR, int nTracks, double pt)
  {
    //// Fill ALL spectra and profiles with one entry per accepted pair ////
    if (isAll) {
      static_for<0, 2>([&](auto j) {
        constexpr int RegionVal = j.value;

        if (!allWindows[RegionVal]) {
          return;
        }

        if (inLR) {
          if (pt > meanPtMin.value) {
            histos.fill(HIST("allEvPkVall/") + HIST(AllLrNames[RegionVal]) + HIST("Pt_") + HIST(BinNames[IdxVal]), pt);
            histos.fill(HIST("allEvPkVall/") + HIST(AllLrNames[RegionVal]) + HIST("PtMultip_") + HIST(BinNames[IdxVal]), nTracks, pt);
          } // for pt > meanPtMin
        }

        if (inSR) {

          if (pt > meanPtMin.value) {
            histos.fill(HIST("allEvPkVall/") + HIST(AllSrNames[RegionVal]) + HIST("Pt_") + HIST(BinNames[IdxVal]), pt);
            histos.fill(HIST("allEvPkVall/") + HIST(AllSrNames[RegionVal]) + HIST("PtMultip_") + HIST(BinNames[IdxVal]), nTracks, pt);
          } // for pt > meanPtMin
        }
      });
    }

    //// Fill TRI spectra and profiles with one entry per accepted pair ////
    if (isTri) {
      static_for<0, 3>([&](auto j) {
        constexpr int RegionVal = j.value;

        if (!triWindows[RegionVal]) {
          return;
        }

        if (inLR) {
          if (pt > meanPtMin.value) {
            histos.fill(HIST("triEvPkVall/") + HIST(TriLrNames[RegionVal]) + HIST("Pt_") + HIST(BinNames[IdxVal]), pt);
            histos.fill(HIST("triEvPkVall/") + HIST(TriLrNames[RegionVal]) + HIST("PtMultip_") + HIST(BinNames[IdxVal]), nTracks, pt);
          } // for pt > meanPtMin
        }

        if (inSR) {
          if (pt > meanPtMin.value) {
            histos.fill(HIST("triEvPkVall/") + HIST(TriSrNames[RegionVal]) + HIST("Pt_") + HIST(BinNames[IdxVal]), pt);
            histos.fill(HIST("triEvPkVall/") + HIST(TriSrNames[RegionVal]) + HIST("PtMultip_") + HIST(BinNames[IdxVal]), nTracks, pt);
          } // for pt > meanPtMin
        }
      });
    }
  }

  //// Fixed-window pT spectra and mean-pT profiles for ALL and TRI ////
  template <typename Track>
  void fillSameEventPtHist(int binIndex,
                           int nTracks,
                           Track const& associatedTrack,
                           double deltaPhi,
                           double deltaEta,
                           bool isLeadAssocPair)
  {
    if (!isLeadAssocPair || leadAssocMode.value != kpTLeadFree) {
      return;
    }

    //// Associated pT and LR/SR pair selections ////
    const double pt = associatedTrack.pt();
    const double absEta = std::abs(deltaEta);
    const double degToRad = constants::math::PI / 180.;

    const bool inLR = absEta > 1.25 && absEta < 1.60;
    const bool inSR = absEta < 0.35;

    //// Preserve the original open angular windows ////
    const auto inPhi = [&](double lowLim, double highLim) {
      return deltaPhi > lowLim * degToRad && deltaPhi < highLim * degToRad;
    };

    //// Event topology ////
    const bool isAll = topologyMode.value == kAllEvSetNom;

    const bool isTri = topologyMode.value == kTriEvSetNom ||
                       topologyMode.value == kTriEvSet1p ||
                       topologyMode.value == kTriEvSet2p ||
                       topologyMode.value == kTriEvSet3p ||
                       topologyMode.value == kTriEvSet1m ||
                       topologyMode.value == kTriEvSet2m ||
                       topologyMode.value == kTriEvSet3m;

    //// ALL: near, between, away ////
    const std::array<bool, 3> allWindows = {
      inPhi(-15., 15.), inPhi(75., 105.), inPhi(165., 195.)};

    //// TRI: peaks at 0, 120, 240 degrees and the two combined valleys ////
    const std::array<bool, 4> triWindows = {
      inPhi(-15., 15.), inPhi(105., 135.), inPhi(225., 255.), inPhi(45., 75.) || inPhi(165., 195.)};

    //// Fill only the requested multiplicity interval ////
    static_for<0, 9>([&](auto i) {
      constexpr int IdxVal = i.value;

      if (binIndex != IdxVal + 1) {
        return;
      }

      fillPtBin<IdxVal>(isAll, isTri, allWindows, triWindows, inLR, inSR, nTracks, pt);
    });
  }

  //// Same-event pair counts and angular correlations only ////
  template <typename Track>
  void fillSameEventHist(int binIndex,
                         Track const& /*associatedTrack*/,
                         double deltaPhi,
                         double deltaEta,
                         bool isLeadAssocPair)
  {

    static_for<0, 9>([&](auto i) {
      constexpr int IdxVal = i.value;

      if (binIndex != IdxVal + 1) {
        return;
      }

      if (!isLeadAssocPair) {
        return;
      }

      //// Count the pair and fill its original angular coordinates ////
      histos.fill(HIST("SamePairs/SamePairsNoLA"), binIndex);
      histos.fill(HIST("Corelations/sProjTRI_") + HIST(BinNames[IdxVal]), deltaPhi);
      histos.fill(HIST("Corelations/sTRI_") + HIST(BinNames[IdxVal]), deltaPhi, deltaEta);
    });
  }

  /// Helper template for FWMs, FWMs vs Multip, FWMs i-j corr.
  template <int R>
  void fillFWMRange(const std::array<double, 9>& H, int nTrk)
  {
    // H[0] = H0 (normalizare), H[1..8] = H1..H8
    static_for<0, 7>([&](auto j) {
      constexpr int JjIt = j.value;
      histos.fill(HIST(FwmRangeNames[R]) + HIST(FwmHlNames[JjIt]), H[JjIt + 1]);
      histos.fill(HIST(FwmRangeNames[R]) + HIST(FwmMulthNames[JjIt]), nTrk, H[JjIt + 1]);
    });

    if (doFWMCorrelations) {
      static_for<0, 27>([&](auto k) {
        constexpr int KkIt = k.value;
        histos.fill(HIST(FwmRangeNames[R]) + HIST(FwmCorrelationNames[KkIt]), H[FwmCorrelationFirst[KkIt] + 1], H[FwmCorrelationSecond[KkIt] + 1]);
      });
    }
  }

  std::array<std::vector<CollisionTuple>, 3> triTuples;

  Filter itsTpcTracks = requireGlobalTrackInFilter();
  using GlobalTracks = soa::Filtered<RecTracks>;

  // Process reconstructed events and buffer selected tracks for mixing.
  void process(CollisionRecTable const& collision,
               GlobalTracks const& tracks)
  {

    // Stop before analysis when the collision fails quality cuts.
    if (!isEventSelected(collision)) {
      return;
    }

    if (tracks.size() == 0) {
      return;
    }

    // Determine track multiplicity, recoil, and leading tracks.
    int nGlobalTracks = tracks.size();
    auto nRecoil = computeRecoil(tracks);

    auto [triggTrk, triggTrkPart, foundPart] = findLeadingPart(tracks);

    // Assign the event to broad and fine multiplicity intervals.
    std::array<bool, MaxMultBins> inBin{};
    for (int i = 0; i < nMultBins; ++i) {
      inBin[i] = (nGlobalTracks >= multLower.value[i] && nGlobalTracks <= multUpper.value[i]);
    }
    const bool multRANGE1 = inBin[0];
    const bool multRANGE2 = inBin[1];
    const bool multRANGE3 = inBin[2];

    // Compute event moments and evaluate topology selections.
    auto fwH = computeFWM(tracks);
    auto [H0, H1, H2, H3, H4, H5, H6, H7, H8] = fwH;

    // FWMs, FWMs vs Multip, FWMs i-j corr.
    if (multRANGE1) {
      fillFWMRange<0>(fwH, nGlobalTracks);
    }

    if (multRANGE2) {
      fillFWMRange<1>(fwH, nGlobalTracks);
    }

    if (multRANGE3) {
      fillFWMRange<2>(fwH, nGlobalTracks);
    }

    // Apply the nominal triangular moment limits
    const bool isTriSetNom = (H1 < kTriFWMSetNom.value[0]) && (H2 < kTriFWMSetNom.value[1]) && (H3 > kTriFWMSetNom.value[2]) && (H4 < kTriFWMSetNom.value[3]) && (H5 > kTriFWMSetNom.value[4]) && (H6 < kTriFWMSetNom.value[5]) && (H7 > kTriFWMSetNom.value[6]) && (H8 < kTriFWMSetNom.value[7]);

    // Apply three progressively looser triangular selections
    const bool isTriSet1p = (H1 < kTriFWMSet1p.value[0]) && (H2 < kTriFWMSet1p.value[1]) && (H3 > kTriFWMSet1p.value[2]) && (H4 < kTriFWMSet1p.value[3]) && (H5 > kTriFWMSet1p.value[4]) && (H6 < kTriFWMSet1p.value[5]) && (H7 > kTriFWMSet1p.value[6]) && (H8 < kTriFWMSet1p.value[7]);
    const bool isTriSet2p = (H1 < kTriFWMSet2p.value[0]) && (H2 < kTriFWMSet2p.value[1]) && (H3 > kTriFWMSet2p.value[2]) && (H4 < kTriFWMSet2p.value[3]) && (H5 > kTriFWMSet2p.value[4]) && (H6 < kTriFWMSet2p.value[5]) && (H7 > kTriFWMSet2p.value[6]) && (H8 < kTriFWMSet2p.value[7]);
    const bool isTriSet3p = (H1 < kTriFWMSet3p.value[0]) && (H2 < kTriFWMSet3p.value[1]) && (H3 > kTriFWMSet3p.value[2]) && (H4 < kTriFWMSet3p.value[3]) && (H5 > kTriFWMSet3p.value[4]) && (H6 < kTriFWMSet3p.value[5]) && (H7 > kTriFWMSet3p.value[6]) && (H8 < kTriFWMSet3p.value[7]);

    // Apply three progressively tighter triangular selections
    const bool isTriSet1m = (H1 < kTriFWMSet1m.value[0]) && (H2 < kTriFWMSet1m.value[1]) && (H3 > kTriFWMSet1m.value[2]) && (H4 < kTriFWMSet1m.value[3]) && (H5 > kTriFWMSet1m.value[4]) && (H6 < kTriFWMSet1m.value[5]) && (H7 > kTriFWMSet1m.value[6]) && (H8 < kTriFWMSet1m.value[7]);
    const bool isTriSet2m = (H1 < kTriFWMSet2m.value[0]) && (H2 < kTriFWMSet2m.value[1]) && (H3 > kTriFWMSet2m.value[2]) && (H4 < kTriFWMSet2m.value[3]) && (H5 > kTriFWMSet2m.value[4]) && (H6 < kTriFWMSet2m.value[5]) && (H7 > kTriFWMSet2m.value[6]) && (H8 < kTriFWMSet2m.value[7]);
    const bool isTriSet3m = (H1 < kTriFWMSet3m.value[0]) && (H2 < kTriFWMSet3m.value[1]) && (H3 > kTriFWMSet3m.value[2]) && (H4 < kTriFWMSet3m.value[3]) && (H5 > kTriFWMSet3m.value[4]) && (H6 < kTriFWMSet3m.value[5]) && (H7 > kTriFWMSet3m.value[6]) && (H8 < kTriFWMSet3m.value[7]);

    const TriFlags triFlags = {isTriSetNom, isTriSet1p, isTriSet2p, isTriSet3p, isTriSet1m, isTriSet2m, isTriSet3m};

    // Apply the configurable isotropic moment limits.
    const bool isIsoNominal = (H1 < isoFWMnominal.value[0]) && (H2 < isoFWMnominal.value[1]) && (H3 < isoFWMnominal.value[2]) && (H4 < isoFWMnominal.value[3]) && (H5 < isoFWMnominal.value[4]) && (H6 < isoFWMnominal.value[5]) && (H7 < isoFWMnominal.value[6]) && (H8 < isoFWMnominal.value[7]);

    // Select all, triangular, or isotropic events for correlations.
    const bool isSelectedTopology =
      (topologyMode.value == kAllEvSetNom) ||
      (topologyMode.value == kTriEvSetNom && isTriSetNom) ||
      (topologyMode.value == kIsoEvSetNom && isIsoNominal) ||
      (topologyMode.value == kTriEvSet1p && isTriSet1p) ||
      (topologyMode.value == kTriEvSet2p && isTriSet2p) ||
      (topologyMode.value == kTriEvSet3p && isTriSet3p) ||
      (topologyMode.value == kTriEvSet1m && isTriSet1m) ||
      (topologyMode.value == kTriEvSet2m && isTriSet2m) ||
      (topologyMode.value == kTriEvSet3m && isTriSet3m);

    // TRI Ev. Multip (the remaining filters) - p == permissive and m == restrictive
    if (triFlags[kTriSet1p]) {
      histos.fill(HIST("Multip/nTrkTriEvS1p"), nGlobalTracks);
    }

    if (triFlags[kTriSet2p]) {
      histos.fill(HIST("Multip/nTrkTriEvS2p"), nGlobalTracks);
    }

    if (triFlags[kTriSet3p]) {
      histos.fill(HIST("Multip/nTrkTriEvS3p"), nGlobalTracks);
    }

    if (triFlags[kTriSet1m]) {
      histos.fill(HIST("Multip/nTrkTriEvS1m"), nGlobalTracks);
    }
    if (triFlags[kTriSet2m]) {
      histos.fill(HIST("Multip/nTrkTriEvS2m"), nGlobalTracks);
    }
    if (triFlags[kTriSet3m]) {
      histos.fill(HIST("Multip/nTrkTriEvS3m"), nGlobalTracks);
    }

    // Fill inclusive event multiplicity, recoil and average pt lead.

    histos.fill(HIST("Multip/multipAllEv"), nGlobalTracks);
    histos.fill(HIST("pTSpectra/RecoilAllEv"), nGlobalTracks, nRecoil);
    histos.fill(HIST("pTSpectra/pTLeadAllEv"), triggTrk.pt());
    histos.fill(HIST("pTSpectra/ptLeadMultipAllEv"), nGlobalTracks, triggTrk.pt());
    histos.fill(HIST("pTSpectra/meanptLeadMultipAllEv"), nGlobalTracks, triggTrk.pt());

    //// Inclusive spectra for the configured selected topology ////
    if (isSelectedTopology) {
      histos.fill(HIST("Multip/multipTriEv"), nGlobalTracks);
      histos.fill(HIST("pTSpectra/RecoilTriEv"), nGlobalTracks, nRecoil);
      histos.fill(HIST("pTSpectra/pTLeadTriEv"), triggTrk.pt());
      histos.fill(HIST("pTSpectra/ptLeadMultipTriEv"), nGlobalTracks, triggTrk.pt());
      histos.fill(HIST("pTSpectra/meanptLeadMultipTriEv"), nGlobalTracks, triggTrk.pt());
    }

    // Fill inclusive and topology-selected leading-track spectra.
    for (const auto& track : tracks) {

      if (track.globalIndex() == triggTrk.globalIndex()) {
        continue;
      } // just the assoc. particles

      histos.fill(HIST("pTSpectra/pTAssocAllEv"), track.pt());
      histos.fill(HIST("pTSpectra/meanpTAssocMultipAllEv"), nGlobalTracks, track.pt());

      if (!isSelectedTopology) {
        continue;
      }

      histos.fill(HIST("pTSpectra/pTAssocTriEv"), track.pt());
      histos.fill(HIST("pTSpectra/meanpTAssocMultipTriEv"), nGlobalTracks, track.pt());
    }

    // Build same-event angular correlations for the selected topology.
    if (doSameEvent && isSelectedTopology) {

      // Form distinct track pairs within the selected event.
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {

        if (track1.globalIndex() == track2.globalIndex()) {
          continue;
        }

        // Calculate pair separations and classify azimuthal windows.
        const double deltaPhi = computeDeltaPhi(track1.phi(), track2.phi());
        const double deltaEta = track1.eta() - track2.eta();

        // Apply the configured pair selection to each interval.
        const bool isLeadAssocPair = passesLeadAssoc(leadAssocMode.value, foundPart, track1, track2, triggTrk, triggTrkPart);

        for (int i = 0; i < nMultBins; ++i) {
          if (inBin[i]) {
            //// Angular correlations for the selected topology ////
            fillSameEventHist(i + 1, track2, deltaPhi, deltaEta, isLeadAssocPair);
            //// Fixed-window spectra and mean-pT profiles ////
            fillSameEventPtHist(i + 1, nGlobalTracks, track2, deltaPhi, deltaEta, isLeadAssocPair);
            //// TRI LR spectra in symmetric 10-degree windows ////
            fillSameEventFoldedLR(i + 1, nGlobalTracks, track2, deltaPhi, deltaEta, isLeadAssocPair);
          }
        }
      }
    }

    // Map the event to a broad range for mixed-event buffering.
    int range = -1;

    if (multRANGE1) {
      range = 0;
    } else if (multRANGE2) {
      range = 1;
    } else if (multRANGE3) {
      range = 2;
    }

    // Buffer up to the requested number of selected events per range.
    if (isSelectedTopology && range >= 0 && triTuples[range].size() < static_cast<size_t>(nSelEv)) {

      std::vector<TrackTuple> tracksTRI;
      tracksTRI.reserve(nGlobalTracks);

      for (const auto& track : tracks) {
        tracksTRI.emplace_back(track.pt(), track.eta(), track.phi());
      }

      triTuples[range].emplace_back(collision.posZ(), nGlobalTracks, std::move(tracksTRI));
    }

    // Fill event counts and spectra for all matching multiplicity intervals.
    for (int i = 0; i < nMultBins; ++i) {
      fillHistoMult(i + 1, multLower.value[i], multUpper.value[i], isSelectedTopology, triFlags, isIsoNominal, tracks);
    }
  }

  int countDF = 0;

  // Write buffered events to the tables consumed by the mixing task.
  void processFILL(aod::Collisions const& /*collisions*/, aod::Tracks const& /*tracks*/)
  {
    // Track data frames and report buffered event counts.
    countDF++;

    // LOGF(info, "Input data - ALL Collisions %d, Tracks %d -> DF no %d", collisions.size(), tracks.size(), countDF);

    for (int range = 0; range < 3; ++range) {

      auto& events = triTuples[range];
      // LOGF(info, "Selected events buffered in RANGE%d: %zu/%d", range + 1, events.size(), static_cast<int>(nSelEv));

      // Wait until a range reaches its event quota.
      if (events.size() < static_cast<size_t>(nSelEv)) {
        continue;
      }

      // LOGF(info, "-> FILL TABLES FOR RANGE%d <-", range + 1);

      // Write each buffered collision and its tracks to the output tables.
      for (const auto& col : events) {

        float posZ = std::get<0>(col);
        int multip = std::get<1>(col);
        triCollisions(posZ, multip);

        for (const auto& track : std::get<2>(col)) {
          triTracks(triCollisions.lastIndex(), std::get<0>(track), std::get<1>(track), std::get<2>(track));
        }
      }

      // Release events after their tables are filled.
      events.clear();
    }
  }

  PROCESS_SWITCH(FoxWolframCorrelation, processFILL, "Fill the triTables", true);
};

// Build mixed-event backgrounds from the selected event tables.
struct FoxWolframCorrelationMixing {

  SliceCache cacheBin1, cacheBin2, cacheBin3, cacheBin4, cacheBin5, cacheBin6, cacheBin7, cacheBin8, cacheBin9, cacheBin10;
  Preslice<aod::TriangleTracks> triCol = aod::track_tri::triangleCollisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr std::array<std::string_view, 10> MixedBinNames = {
    "BIN1", "BIN2", "BIN3", "BIN4", "BIN5", "BIN6", "BIN7", "BIN8", "BIN9", "BIN10"};

  Configurable<int> nBin{"nBin", 100, "No of bis per axis"};
  ConfigurableAxis axisVtx{"axisVtx", {VARIABLE_WIDTH, -10, 10}, "vertex axis for mixed event histograms"};

  // Use half-integer edges for inclusive integer multiplicity bins.
  ConfigurableAxis axisMult1{"axisMult1", {VARIABLE_WIDTH, 39.5, 85.5}, "axisMult1"};
  ConfigurableAxis axisMult2{"axisMult2", {VARIABLE_WIDTH, 85.5, 140.5}, "axisMult2"};
  ConfigurableAxis axisMult3{"axisMult3", {VARIABLE_WIDTH, 140.5, 350.5}, "axisMult3"};
  ConfigurableAxis axisMult4{"axisMult4", {VARIABLE_WIDTH, 39.5, 60.5}, "axisMult4"};
  ConfigurableAxis axisMult5{"axisMult5", {VARIABLE_WIDTH, 60.5, 80.5}, "axisMult5"};
  ConfigurableAxis axisMult6{"axisMult6", {VARIABLE_WIDTH, 80.5, 100.5}, "axisMult6"};
  ConfigurableAxis axisMult7{"axisMult7", {VARIABLE_WIDTH, 100.5, 120.5}, "axisMult7"};
  ConfigurableAxis axisMult8{"axisMult8", {VARIABLE_WIDTH, 120.5, 140.5}, "axisMult8"};
  ConfigurableAxis axisMult9{"axisMult9", {VARIABLE_WIDTH, 140.5, 160.5}, "axisMult9"};
  ConfigurableAxis axisMult10{"axisMult10", {VARIABLE_WIDTH, 160.5, 350.5}, "axisMult10"};

  Configurable<bool> doMixedEvent{"doMixedEvent", true, "Fill mixed-event TRI correlations"};
  Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of mixed events per trigger event"};

  // Same multiplicity limits as in FoxWolframCorrelation (inclusive, same values in both tasks).
  // They set how many mixing pools are active and are checked against axisMult1..10.
  Configurable<std::vector<int>> multLower{"multLower", {40, 86, 141, 40, 61, 81, 101, 121, 141, 161}, "lower multip. limits (inclusive)"};
  Configurable<std::vector<int>> multUpper{"multUpper", {85, 140, 350, 60, 80, 100, 120, 140, 160, 350}, "upper multip. limits (inclusive)"};

  static constexpr int MaxMultBins = 10;
  int nMultBins = 0; // set in init() from multLower/multUpper

  // Choose the pair selection for mixed-event correlations.
  Configurable<int> leadAssocMode{"leadAssocMode", 2, "0: ALL i-j pairs (pT<3); 1: pTLead in (2,3) & pTAssoc < 2; 2: pTAssoc < pTLead"};

  // Register mixed-event pair counts and angular distributions.
  void init(InitContext const&)
  {
    // Validate the configurable multiplicity bins.
    if (multLower.value.size() != multUpper.value.size() || multLower.value.empty() || multLower.value.size() > static_cast<size_t>(MaxMultBins)) {
      LOGF(fatal, "multLower/multUpper must have the same size, between 1 and %d entries", MaxMultBins);
    }
    nMultBins = static_cast<int>(multLower.value.size());

    // The binning policies are fixed objects, so check that each axisMultK matches [multLower-0.5, multUpper+0.5].
    const std::array<std::vector<double>, MaxMultBins> axisEdges = {
      axisMult1.value, axisMult2.value, axisMult3.value, axisMult4.value, axisMult5.value,
      axisMult6.value, axisMult7.value, axisMult8.value, axisMult9.value, axisMult10.value};

    //// Each mixing axis must cover its configured integer multiplicity interval ////
    for (int i = 0; i < nMultBins; ++i) {
      const auto& edges = axisEdges[i];
      const double lowLim = multLower.value[i] - 0.5;
      const double highLim = multUpper.value[i] + 0.5;
      if (multLower.value[i] > multUpper.value[i] || edges.size() < 2 ||
          std::abs(edges[edges.size() - 2] - lowLim) > 1.e-6 ||
          std::abs(edges.back() - highLim) > 1.e-6) {
        LOGF(fatal, "axisMult%d must match multiplicity limits [%.1f, %.1f]", i + 1, lowLim, highLim);
      }
    }

    pools = {&pairBIN1, &pairBIN2, &pairBIN3, &pairBIN4, &pairBIN5,
             &pairBIN6, &pairBIN7, &pairBIN8, &pairBIN9, &pairBIN10};

    const AxisSpec axisDeltaEta{32, -1.6, +1.6, "#Delta#eta"};
    const AxisSpec axisDeltaPhi{72, -constants::math::PIHalf, +3 * constants::math::PIHalf, "#Delta#phi"};
    const AxisSpec axisEvent{10, 0.5, 10.5, "#Multip. Bin", "EventAxis"};

    // Label mixed-event pair counts by multiplicity interval.
    histos.add("MixingPairs/MixingPairsNoLA", "MixingPairsNoLA", kTH1D, {axisEvent}, false);
    auto* mixingPairsNoLA = histos.get<TH1>(HIST("MixingPairs/MixingPairsNoLA"))->GetXaxis();
    for (int i = 0; i < MaxMultBins; ++i) {
      mixingPairsNoLA->SetBinLabel(i + 1, i < 3 ? fmt::format("mixPairsLA_RANGE{}", i + 1).c_str() : fmt::format("mixPairsLA_BIN{}", i + 1).c_str());
    }

    // Register mixed-event correlation histograms for all intervals.
    for (int i = 1; i <= nMultBins; i++) {
      histos.add(fmt::format("Corelations/bTRI_BIN{}", i).c_str(), fmt::format("bTRI_BIN{}", i).c_str(), kTH2D, {axisDeltaPhi, axisDeltaEta});
      histos.add(fmt::format("Corelations/bProjTRI_BIN{}", i).c_str(), fmt::format("bProjTRI_BIN{}", i).c_str(), kTH1D, {axisDeltaPhi});
    }
  }

  // Mix events in vertex and multiplicity intervals.
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::Multip>;

  BinningType binningOnPosMultip1{{axisVtx, axisMult1}};
  BinningType binningOnPosMultip2{{axisVtx, axisMult2}};
  BinningType binningOnPosMultip3{{axisVtx, axisMult3}};
  BinningType binningOnPosMultip4{{axisVtx, axisMult4}};
  BinningType binningOnPosMultip5{{axisVtx, axisMult5}};
  BinningType binningOnPosMultip6{{axisVtx, axisMult6}};
  BinningType binningOnPosMultip7{{axisVtx, axisMult7}};
  BinningType binningOnPosMultip8{{axisVtx, axisMult8}};
  BinningType binningOnPosMultip9{{axisVtx, axisMult9}};
  BinningType binningOnPosMultip10{{axisVtx, axisMult10}};

  // Create independent mixing pools for the ten multiplicity intervals.
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN1{binningOnPosMultip1, nMixedEvents, -1, &cacheBin1};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN2{binningOnPosMultip2, nMixedEvents, -1, &cacheBin2};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN3{binningOnPosMultip3, nMixedEvents, -1, &cacheBin3};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN4{binningOnPosMultip4, nMixedEvents, -1, &cacheBin4};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN5{binningOnPosMultip5, nMixedEvents, -1, &cacheBin5};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN6{binningOnPosMultip6, nMixedEvents, -1, &cacheBin6};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN7{binningOnPosMultip7, nMixedEvents, -1, &cacheBin7};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN8{binningOnPosMultip8, nMixedEvents, -1, &cacheBin8};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN9{binningOnPosMultip9, nMixedEvents, -1, &cacheBin9};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN10{binningOnPosMultip10, nMixedEvents, -1, &cacheBin10};

  // Pointers to the pools so they can be processed in a loop (filled in init()).
  using MixPair = SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType>;
  std::array<MixPair*, MaxMultBins> pools{};

  // Fill mixed-event distributions for a selected multiplicity interval.
  template <typename PairType>
  void fillMixedEventBin(PairType& pair, int binIndex)
  {
    // Mix populated event pairs within the selected pool.
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {
      if (tracks1.size() == 0 || tracks2.size() == 0) {
        continue;
      }

      // Find leading tracks in the first mixed event.
      auto [triggTrk, triggTrkPart, foundPart] = findLeadingPart(tracks1);

      // LOGF(debug, "Mix TRI Ev (BIN%d): (%d, %d), z=(%.3f, %.3f), mult=(%d, %d)",
      //      binIndex, c1.globalIndex(), c2.globalIndex(), c1.posZ(), c2.posZ(), tracks1.size(), tracks2.size());

      // Select tracks from different events using the configured pair rules.
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!passesLeadAssoc(leadAssocMode.value, foundPart, track1, track2, triggTrk, triggTrkPart)) {
          continue;
        }

        const double deltaPhi = computeDeltaPhi(track1.phi(), track2.phi());
        const double deltaEta = track1.eta() - track2.eta();

        // Count accepted pairs in the corresponding multiplicity interval.
        histos.fill(HIST("MixingPairs/MixingPairsNoLA"), binIndex);
        static_for<0, 9>([&](auto i) {
          constexpr int IdxVal = i.value;
          if (binIndex == IdxVal + 1) {
            histos.fill(HIST("Corelations/bProjTRI_") + HIST(MixedBinNames[IdxVal]), deltaPhi);
            histos.fill(HIST("Corelations/bTRI_") + HIST(MixedBinNames[IdxVal]), deltaPhi, deltaEta);
          }
        });
      }
    }
  }

  // Run mixed-event pairing for all multiplicity intervals.
  void process(aod::TriangleCollisions const& collisions,
               aod::TriangleTracks const& tracks)
  {
    // Skip mixing when disabled or when the input tables are empty.
    if (!doMixedEvent || collisions.size() == 0 || tracks.size() == 0) {
      return;
    }

    // Process the active independent mixing pools.
    for (int i = 0; i < nMultBins; ++i) {
      fillMixedEventBin(*pools[i], i + 1);
    }
  }
};

// Expose the event analysis and mixed-event tasks as one workflow.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FoxWolframCorrelation>(cfgc),
    adaptAnalysisTask<FoxWolframCorrelationMixing>(cfgc),
  };
}
