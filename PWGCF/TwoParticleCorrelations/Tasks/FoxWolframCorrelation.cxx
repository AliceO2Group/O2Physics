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

/// \file FoxWolframCorrelation.cxx
/// \author Andreicovici Iulian Florin
/// \brief RECO-level event shape analysis using Fox-Wolfram Moments for O+O collisions.
///        Classifies events by topology (triangular), fills correlation histograms,
///        and produces custom SOA tables for mixed-event background estimation.
/// \since 25 Sept 2026

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
enum LeadAssocMode : int {
  kAllPairs = 0,
  kAllPairsEtaConstrain = 1,
  kpTLeadConstrain = 2,
  kpTLeadFree = 3
};

// Name the available event-topology selections.
enum TopologyMode : int {
  kAllEvents = 0,
  kTriEvents = 1,
  kIsoEvents = 2
};

template <typename Track, typename TLead>
// Apply the same pair-selection rules to signal and background.
bool passesLeadAssoc(int mode, bool foundPart, Track const& track1, Track const& track2,
                     TLead const& triggTrk, TLead const& triggTrkPart)
{
  if (mode == kAllPairs) {
    return (track1.pt() < 3.f) && (track2.pt() < 3.f);
  }
  if (mode == kAllPairsEtaConstrain) {

    const bool firstPositiveSecondNegative =
      track1.eta() > 0.1f && track1.eta() < 0.8f &&
      track2.eta() > -0.8f && track2.eta() < -0.1f;

    const bool firstNegativeSecondPositive =
      track1.eta() > -0.8f && track1.eta() < -0.1f &&
      track2.eta() > 0.1f && track2.eta() < 0.8f;

    return track1.pt() < 3.f && track2.pt() < 3.f &&
           (firstPositiveSecondNegative || firstNegativeSecondPositive);
  }
  if (mode == kpTLeadConstrain) {
    return foundPart && (track1.globalIndex() == triggTrkPart.globalIndex()) && (track2.pt() < 2.f);
  }
  if (mode == kpTLeadFree) {
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

  for (const auto& track : tracks) {
    if (track.pt() > ptMaxAll) {
      ptMaxAll = track.pt();
      leadingAllId = i;
    }
    if (track.pt() > 2.f && track.pt() < 3.f && track.pt() > ptMaxPart) {
      ptMaxPart = track.pt();
      leadingPartId = i;
      foundPart = true;
    }
    i++;
  }

  return std::make_tuple(tracks.iteratorAt(leadingAllId),
                         tracks.iteratorAt(leadingPartId),
                         foundPart);
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

  // Configure histogram axes, event selection, and topology.
  Configurable<int> nBins{"nBins", 100, "no. bins in all histos"};
  Configurable<int> nSelEv{"nSelEv", 200, "no. of mixing Ev."};
  Configurable<float> vtxRange{"vtxRange", 10.0f, "Vertex Z range to consider"};
  Configurable<float> etaRange{"etaRange", 0.8f, "eta range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "custom DCA Z cut (ignored if negative)"};
  Configurable<bool> doSameEvent{"doSameEvent", true, "Fill same-event TRI correlations"};
  Configurable<bool> doFWMCorrelations{"doFWMCorrelations", false,
                                       "Fill H_i-H_j correlation histograms for RANGE1-RANGE3"};
  Configurable<int> topologyMode{"topologyMode", kTriEvents,
                                 "Topology used for selected-event outputs: 0 ALL, 1 TRI, 2 ISO"};

  // Choose the pair selection for same-event correlations.
  Configurable<int> leadAssocMode{"leadAssocMode", 0,
                                  "0: ALL i-j pairs (pT<3); 1: ALL i-j pairs (pT < 3, Eta constrain) 2: pTLead in (2,3) & pTAssoc < 2; 3: pTAssoc < pTLead"};

  Configurable<bool> isVtxRange{"isVtxRange", true, "isVtxRange"};
  Configurable<bool> isSel8{"isSel8", true, "isSel8"};
  Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", true, "isNoSameBunchPileup"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", true, "isGoodZvtxFT0vsPV"};
  Configurable<bool> isVertexITSTPC{"isVertexITSTPC", false, "isVertexITSTPC"};
  Configurable<bool> isVertexTOFmatched{"isVertexTOFmatched", false, "isVertexTOFmatched"};
  Configurable<bool> kIsGoodITSLayer0123{"kIsGoodITSLayer0123", false, "kIsGoodITSLayer0123"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", false, "isNoCollInTimeRangeNarrow"};

  // Set the isotropic event selection thresholds.
  Configurable<double> h1IsoMax{"h1IsoMax", 0.02, "h1IsoMax"};
  Configurable<double> h2IsoMax{"h2IsoMax", 0.27, "h2IsoMax"};
  Configurable<double> h3IsoMax{"h3IsoMax", 0.02, "h3IsoMax"};
  Configurable<double> h4IsoMax{"h4IsoMax", 0.16, "h4IsoMax"};
  Configurable<double> h5IsoMax{"h5IsoMax", 0.02, "h5IsoMax"};
  Configurable<double> h6IsoMax{"h6IsoMax", 0.12, "h6IsoMax"};
  Configurable<double> h7IsoMax{"h7IsoMax", 0.02, "h7IsoMax"};
  Configurable<double> h8IsoMax{"h8IsoMax", 0.10, "h8IsoMax"};

  // Define the three broad multiplicity ranges.
  Configurable<int> range1Min{"range1Min", 40, "minRANGE1 "};
  Configurable<int> range1Max{"range1Max", 85, "maxRANGE1"};
  Configurable<int> range2Min{"range2Min", 86, "minRANGE2"};
  Configurable<int> range2Max{"range2Max", 140, "maxRANGE2"};
  Configurable<int> range3Min{"range3Min", 141, "minRANGE3"};
  Configurable<int> range3Max{"range3Max", 280, "maxRANGE3"};

  // Define the seven finer multiplicity intervals.
  Configurable<int> bin4Min{"bin4Min", 40, "minBIN4"};
  Configurable<int> bin4Max{"bin4Max", 60, "maxBIN4"};
  Configurable<int> bin5Min{"bin5Min", 61, "minBIN5"};
  Configurable<int> bin5Max{"bin5Max", 80, "maxBIN5"};
  Configurable<int> bin6Min{"bin6Min", 81, "minBIN6"};
  Configurable<int> bin6Max{"bin6Max", 100, "maxBIN6"};
  Configurable<int> bin7Min{"bin7Min", 101, "minBIN7"};
  Configurable<int> bin7Max{"bin7Max", 120, "maxBIN7"};
  Configurable<int> bin8Min{"bin8Min", 121, "minBIN8"};
  Configurable<int> bin8Max{"bin8Max", 140, "maxBIN8"};
  Configurable<int> bin9Min{"bin9Min", 141, "minBIN9"};
  Configurable<int> bin9Max{"bin9Max", 160, "maxBIN9"};
  Configurable<int> bin10Min{"bin10Min", 161, "minBIN10"};
  Configurable<int> bin10Max{"bin10Max", 280, "maxBIN10"};

  Configurable<int> noOfMultipBins{"noOfMultipBins", 10, "noOfMultipBins"};

  static constexpr std::array<std::string_view, 10> binNames = {
    "BIN1", "BIN2", "BIN3", "BIN4", "BIN5",
    "BIN6", "BIN7", "BIN8", "BIN9", "BIN10"};

  // Map each unique pair of Fox-Wolfram moments to its histogram.
  static constexpr std::array<std::string_view, 28> fwmCorrelationNames = {"H1_H2", "H1_H3", "H1_H4", "H1_H5", "H1_H6", "H1_H7", "H1_H8", "H2_H3", "H2_H4", "H2_H5", "H2_H6", "H2_H7", "H2_H8", "H3_H4", "H3_H5", "H3_H6", "H3_H7", "H3_H8", "H4_H5", "H4_H6", "H4_H7", "H4_H8", "H5_H6", "H5_H7", "H5_H8", "H6_H7", "H6_H8", "H7_H8"};
  static constexpr std::array<int, 28> fwmCorrelationFirst = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6};
  static constexpr std::array<int, 28> fwmCorrelationSecond = {1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7};

  // Register event, track, moment, and correlation histograms.
  void init(InitContext const&)
  {
    // Reject unsupported topology modes at startup.
    if (topologyMode.value < kAllEvents || topologyMode.value > kIsoEvents) {
      LOGF(fatal, "Invalid topologyMode=%d. Allowed values: 0 ALL, 1 TRI, 2 ISO", topologyMode.value);
    }

    // Define axes shared by the registered histograms.
    const AxisSpec axispT{100, 0.15, 10, "p_{T}"};
    const AxisSpec axispTL{nBins, 2, 3, "p_{T}^{Lead}"};
    const AxisSpec axispTA{nBins, 0.15, 2, "p_{T}^{Assoc}"};
    const AxisSpec axisDeltaEta{32, -1.6, +1.6, "#Delta#eta"};
    const AxisSpec axisDeltaPhi{72, -constants::math::PIHalf, +3 * constants::math::PIHalf, "#Delta#phi"};
    const AxisSpec axisPhi{100, 0, constants::math::TwoPI, "#phi"};
    const AxisSpec axisEta{40, -0.8, 0.8, "#eta"};
    const AxisSpec axisMultip{350, 0, 350, "N_{ch}"};
    const AxisSpec axisHlrec{nBins / 2, 0, 1, "H_{l}^{rec}"};
    const AxisSpec axisHlfine{nBins, 0, 1, "H_{l}^{rec}"};
    const AxisSpec axisEvent{10, 0.5, 10.5, "", "EventAxis"};
    const AxisSpec axisRecoil{nBins, 0., 1, "Recoil"};

    // Label same-event pair counts by multiplicity interval.
    histos.add("SamePairs/SamePairsNoLA", "SamePairsNoLA", kTH1D, {axisEvent}, false);
    auto hstatLA = histos.get<TH1>(HIST("SamePairs/SamePairsNoLA"));
    auto* hLA = hstatLA->GetXaxis();
    hLA->SetBinLabel(1, "sPairsLA_RANGE1");
    hLA->SetBinLabel(2, "sPairsLA_RANGE2");
    hLA->SetBinLabel(3, "sPairsLA_RANGE3");
    hLA->SetBinLabel(4, "sPairsLA_BIN4");
    hLA->SetBinLabel(5, "sPairsLA_BIN5");
    hLA->SetBinLabel(6, "sPairsLA_BIN6");
    hLA->SetBinLabel(7, "sPairsLA_BIN7");
    hLA->SetBinLabel(8, "sPairsLA_BIN8");
    hLA->SetBinLabel(9, "sPairsLA_BIN9");
    hLA->SetBinLabel(10, "sPairsLA_BIN10");

    // Label event counts for all multiplicity intervals.
    histos.add("Events_Details/eventsNo", "No of EVENTS", kTH1D, {axisEvent}, false);
    auto hstat2 = histos.get<TH1>(HIST("Events_Details/eventsNo"));
    auto* h2 = hstat2->GetXaxis();
    h2->SetBinLabel(1, "RANGE1");
    h2->SetBinLabel(2, "RANGE2");
    h2->SetBinLabel(3, "RANGE3");
    h2->SetBinLabel(4, "BIN4");
    h2->SetBinLabel(5, "BIN5");
    h2->SetBinLabel(6, "BIN6");
    h2->SetBinLabel(7, "BIN7");
    h2->SetBinLabel(8, "BIN8");
    h2->SetBinLabel(9, "BIN9");
    h2->SetBinLabel(10, "BIN10");

    // Count events passing the nominal triangular selection.
    histos.add("Events_Details/eventsNoTRI", "No of TRIANGLE EVENTS", kTH1D, {axisEvent}, false);
    auto hstat3 = histos.get<TH1>(HIST("Events_Details/eventsNoTRI"));
    auto* h3 = hstat3->GetXaxis();
    h3->SetBinLabel(1, "RANGE1");
    h3->SetBinLabel(2, "RANGE2");
    h3->SetBinLabel(3, "RANGE3");
    h3->SetBinLabel(4, "BIN4");
    h3->SetBinLabel(5, "BIN5");
    h3->SetBinLabel(6, "BIN6");
    h3->SetBinLabel(7, "BIN7");
    h3->SetBinLabel(8, "BIN8");
    h3->SetBinLabel(9, "BIN9");
    h3->SetBinLabel(10, "BIN10");

    // Count events passing the first alternate triangular selection.
    histos.add("Events_Details/eventsNoTRI1", "No of TRI events selected by SET1", kTH1D, {axisEvent}, false);
    auto hstatTRI1 = histos.get<TH1>(HIST("Events_Details/eventsNoTRI1"));
    auto* hTRI1 = hstatTRI1->GetXaxis();
    hTRI1->SetBinLabel(1, "RANGE1");
    hTRI1->SetBinLabel(2, "RANGE2");
    hTRI1->SetBinLabel(3, "RANGE3");
    hTRI1->SetBinLabel(4, "BIN4");
    hTRI1->SetBinLabel(5, "BIN5");
    hTRI1->SetBinLabel(6, "BIN6");
    hTRI1->SetBinLabel(7, "BIN7");
    hTRI1->SetBinLabel(8, "BIN8");
    hTRI1->SetBinLabel(9, "BIN9");
    hTRI1->SetBinLabel(10, "BIN10");

    // Count events passing the second alternate triangular selection.
    histos.add("Events_Details/eventsNoTRI2", "No of TRI events selected by SET2", kTH1D, {axisEvent}, false);
    auto hstatTRI2 = histos.get<TH1>(HIST("Events_Details/eventsNoTRI2"));
    auto* hTRI2 = hstatTRI2->GetXaxis();
    hTRI2->SetBinLabel(1, "RANGE1");
    hTRI2->SetBinLabel(2, "RANGE2");
    hTRI2->SetBinLabel(3, "RANGE3");
    hTRI2->SetBinLabel(4, "BIN4");
    hTRI2->SetBinLabel(5, "BIN5");
    hTRI2->SetBinLabel(6, "BIN6");
    hTRI2->SetBinLabel(7, "BIN7");
    hTRI2->SetBinLabel(8, "BIN8");
    hTRI2->SetBinLabel(9, "BIN9");
    hTRI2->SetBinLabel(10, "BIN10");

    // Count events passing the third alternate triangular selection.
    histos.add("Events_Details/eventsNoTRI3", "No of TRI events selected by SET3", kTH1D, {axisEvent}, false);
    auto hstatTRI3 = histos.get<TH1>(HIST("Events_Details/eventsNoTRI3"));
    auto* hTRI3 = hstatTRI3->GetXaxis();
    hTRI3->SetBinLabel(1, "RANGE1");
    hTRI3->SetBinLabel(2, "RANGE2");
    hTRI3->SetBinLabel(3, "RANGE3");
    hTRI3->SetBinLabel(4, "BIN4");
    hTRI3->SetBinLabel(5, "BIN5");
    hTRI3->SetBinLabel(6, "BIN6");
    hTRI3->SetBinLabel(7, "BIN7");
    hTRI3->SetBinLabel(8, "BIN8");
    hTRI3->SetBinLabel(9, "BIN9");
    hTRI3->SetBinLabel(10, "BIN10");

    // Count events passing the tighter auxiliary triangular selection.
    histos.add("Events_Details/eventsNoTRIAux", "No of TRI events selected by auxiliary SET", kTH1D, {axisEvent}, false);
    auto hstatTRIAux = histos.get<TH1>(HIST("Events_Details/eventsNoTRIAux"));
    auto* hTRIAux = hstatTRIAux->GetXaxis();
    hTRIAux->SetBinLabel(1, "RANGE1");
    hTRIAux->SetBinLabel(2, "RANGE2");
    hTRIAux->SetBinLabel(3, "RANGE3");
    hTRIAux->SetBinLabel(4, "BIN4");
    hTRIAux->SetBinLabel(5, "BIN5");
    hTRIAux->SetBinLabel(6, "BIN6");
    hTRIAux->SetBinLabel(7, "BIN7");
    hTRIAux->SetBinLabel(8, "BIN8");
    hTRIAux->SetBinLabel(9, "BIN9");
    hTRIAux->SetBinLabel(10, "BIN10");

    // Count events passing the isotropic selection.
    histos.add("Events_Details/eventsNoISO", "No of ISO EVENTS", kTH1D, {axisEvent}, false);
    auto hstat4 = histos.get<TH1>(HIST("Events_Details/eventsNoISO"));
    auto* h4 = hstat4->GetXaxis();
    h4->SetBinLabel(1, "RANGE1");
    h4->SetBinLabel(2, "RANGE2");
    h4->SetBinLabel(3, "RANGE3");
    h4->SetBinLabel(4, "BIN4");
    h4->SetBinLabel(5, "BIN5");
    h4->SetBinLabel(6, "BIN6");
    h4->SetBinLabel(7, "BIN7");
    h4->SetBinLabel(8, "BIN8");
    h4->SetBinLabel(9, "BIN9");
    h4->SetBinLabel(10, "BIN10");

    // Register inclusive multiplicity, track spectra, and recoil histograms.
    histos.add("Multip/ITSTPC_Multiplicity", "Reconstructed ITS-TPC Multiplicity; trks; Events", kTH1D, {axisMultip});

    histos.add("pTSpectra/pT", "pT", kTH1D, {axispT});
    histos.add("Spectra/phi", "phi", kTH1D, {axisPhi});

    histos.add("pTSpectra/RecoilALL", "RecoilALL", kTH2D, {axisMultip, axisRecoil});
    histos.add("pTSpectra/RecoilTRI", "RecoilTRI", kTH2D, {axisMultip, axisRecoil});

    // Register track spectra for each multiplicity interval.
    for (int i = 1; i <= noOfMultipBins; i++) {
      histos.add(fmt::format("Multip/GlobalTrk_BIN{}", i).c_str(), fmt::format("GlobalTrk_BIN{}", i).c_str(), kTH1D, {axisMultip});
      histos.add(fmt::format("Multip/GlobalTrkTri_BIN{}", i).c_str(), fmt::format("GlobalTrkTri_BIN{}", i).c_str(), kTH1D, {axisMultip});

      histos.add(fmt::format("pTSpectra/pT_BIN{}", i).c_str(), fmt::format("pT_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("pTSpectra/pTTri_BIN{}", i).c_str(), fmt::format("pTTri_BIN{}", i).c_str(), kTH1D, {axispT});

      histos.add(fmt::format("Spectra/phi_BIN{}", i).c_str(), fmt::format("phi_BIN{}", i).c_str(), kTH1D, {axisPhi});
      histos.add(fmt::format("Spectra/phiTri_BIN{}", i).c_str(), fmt::format("phiTri_BIN{}", i).c_str(), kTH1D, {axisPhi});

      histos.add(fmt::format("Spectra/etaphi_BIN{}", i).c_str(), fmt::format("etaphi_BIN{}", i).c_str(), kTH2D, {axisPhi, axisEta});
      histos.add(fmt::format("Spectra/etaphiTri_BIN{}", i).c_str(), fmt::format("etaphiTri_BIN{}", i).c_str(), kTH2D, {axisPhi, axisEta});

      histos.add(fmt::format("Spectra/ptphi_BIN{}", i).c_str(), fmt::format("ptphi_BIN{}", i).c_str(), kTH2D, {axisPhi, axispT});
      histos.add(fmt::format("Spectra/ptphiTri_BIN{}", i).c_str(), fmt::format("ptphiTri_BIN{}", i).c_str(), kTH2D, {axisPhi, axispT});
    }

    // Register leading and associated track spectra in the broad ranges.
    for (int i = 1; i <= 3; i++) {
      histos.add(fmt::format("pTSpectra/pT_LA_BIN{}", i).c_str(), fmt::format("pT_LA_BIN{}", i).c_str(), kTH2D, {axispT, axispT});
      histos.add(fmt::format("pTSpectra/pT_LA_TriBIN{}", i).c_str(), fmt::format("pT_LA_BIN{}", i).c_str(), kTH2D, {axispT, axispT});

      histos.add(fmt::format("pTSpectra/fpT_LA_BIN{}", i).c_str(), fmt::format("fpT_LA_BIN{}", i).c_str(), kTH2D, {axispTL, axispTA});
      histos.add(fmt::format("pTSpectra/fpT_LA_TriBIN{}", i).c_str(), fmt::format("fpT_LA_BIN{}", i).c_str(), kTH2D, {axispTL, axispTA});

      histos.add(fmt::format("pTSpectra/peakPT_BIN{}", i).c_str(), fmt::format("peakPT_BIN{}", i).c_str(), kTH1D, {axispT});
      histos.add(fmt::format("pTSpectra/valleyPT_BIN{}", i).c_str(), fmt::format("valleyPT_BIN{}", i).c_str(), kTH1D, {axispT});
    }

    // Register same-event angular correlations.
    for (int i = 1; i <= noOfMultipBins; ++i) {
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

      for (size_t k = 0; k < fwmCorrelationNames.size(); ++k) {

        const int first = fwmCorrelationFirst[k] + 1;
        const int second = fwmCorrelationSecond[k] + 1;

        histos.add(fmt::format("FWM_BIN1/{}", fwmCorrelationNames[k]).c_str(), fmt::format("H_{{{}}} vs H_{{{}}}", first, second).c_str(), kTH2D, {axisHlfine, axisHlfine});
        histos.add(fmt::format("FWM_BIN2/{}", fwmCorrelationNames[k]).c_str(), fmt::format("H_{{{}}} vs H_{{{}}}", first, second).c_str(), kTH2D, {axisHlfine, axisHlfine});
        histos.add(fmt::format("FWM_BIN3/{}", fwmCorrelationNames[k]).c_str(), fmt::format("H_{{{}}} vs H_{{{}}}", first, second).c_str(), kTH2D, {axisHlfine, axisHlfine});
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
                     bool isSelectedTopology, bool isTRI, bool isTRI1,
                     bool isTRI2, bool isTRI3, bool isTRIAux, bool isISO,
                     const T& tracks)
  {
    if (tracks.size() < minMult || tracks.size() > maxMult) {
      return;
    }

    // Count all events and each topology selection in this interval.
    histos.fill(HIST("Events_Details/eventsNo"), binIndex);

    if (isISO) {
      histos.fill(HIST("Events_Details/eventsNoISO"), binIndex);
    }

    if (isTRI) {
      histos.fill(HIST("Events_Details/eventsNoTRI"), binIndex);
    }

    if (isTRI1) {
      histos.fill(HIST("Events_Details/eventsNoTRI1"), binIndex);
    }
    if (isTRI2) {
      histos.fill(HIST("Events_Details/eventsNoTRI2"), binIndex);
    }
    if (isTRI3) {
      histos.fill(HIST("Events_Details/eventsNoTRI3"), binIndex);
    }
    if (isTRIAux) {
      histos.fill(HIST("Events_Details/eventsNoTRIAux"), binIndex);
    }

    // Fill inclusive spectra and selected-topology spectra for this interval.
    static_for<0, 9>([&](auto i) {
      constexpr int idx = i.value;

      if (binIndex == idx + 1) {

        histos.fill(HIST("Multip/GlobalTrk_") + HIST(binNames[idx]), tracks.size());

        for (const auto& track : tracks) {

          histos.fill(HIST("pTSpectra/pT_") + HIST(binNames[idx]), track.pt());
          histos.fill(HIST("Spectra/phi_") + HIST(binNames[idx]), track.phi());

          histos.fill(HIST("Spectra/etaphi_") + HIST(binNames[idx]), track.phi(), track.eta());
          histos.fill(HIST("Spectra/ptphi_") + HIST(binNames[idx]), track.phi(), track.pt());
        }

        // Fill spectra for the configured topology selection.
        if (isSelectedTopology) {

          histos.fill(HIST("Multip/GlobalTrkTri_") + HIST(binNames[idx]), tracks.size());

          for (const auto& track : tracks) {
            histos.fill(HIST("pTSpectra/pTTri_") + HIST(binNames[idx]), track.pt());
            histos.fill(HIST("Spectra/phiTri_") + HIST(binNames[idx]), track.phi());

            histos.fill(HIST("Spectra/etaphiTri_") + HIST(binNames[idx]), track.phi(), track.eta());
            histos.fill(HIST("Spectra/ptphiTri_") + HIST(binNames[idx]), track.phi(), track.pt());
          }
        }
      }
    });
  }

  // Fill same-event pair distributions for one multiplicity interval.
  template <typename Track>
  void fillSameEventHist(int binIndex,
                         Track const& associatedTrack,
                         double deltaPhi,
                         double deltaEta,
                         bool onPeaks,
                         bool onValleys,
                         bool isLeadAssocPair)
  {

    static_for<0, 9>([&](auto i) {
      constexpr int idx = i.value;

      if (binIndex != idx + 1) {
        return;
      }

      if (!isLeadAssocPair) {
        return;
      }

      // Record peak and valley track spectra in the broad ranges.
      if constexpr (idx < 3) {
        if (leadAssocMode.value == 3) {
          if (onPeaks) {
            histos.fill(HIST("pTSpectra/peakPT_") + HIST(binNames[idx]), associatedTrack.pt());
          }
          if (onValleys) {
            histos.fill(HIST("pTSpectra/valleyPT_") + HIST(binNames[idx]), associatedTrack.pt());
          }
        }
      }

      // Count accepted pairs and fill their angular distributions.
      histos.fill(HIST("SamePairs/SamePairsNoLA"), binIndex);
      histos.fill(HIST("Corelations/sProjTRI_") + HIST(binNames[idx]), deltaPhi);
      histos.fill(HIST("Corelations/sTRI_") + HIST(binNames[idx]), deltaPhi, deltaEta);
    });
  }

  std::array<std::vector<CollisionTuple>, 3> triTuples;

  Filter ITSTPCTracks = requireGlobalTrackInFilter();
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
    bool multRANGE1 = (nGlobalTracks >= range1Min && nGlobalTracks <= range1Max);
    bool multRANGE2 = (nGlobalTracks >= range2Min && nGlobalTracks <= range2Max);
    bool multRANGE3 = (nGlobalTracks >= range3Min && nGlobalTracks <= range3Max);
    bool multBIN4 = (nGlobalTracks >= bin4Min && nGlobalTracks <= bin4Max);
    bool multBIN5 = (nGlobalTracks >= bin5Min && nGlobalTracks <= bin5Max);
    bool multBIN6 = (nGlobalTracks >= bin6Min && nGlobalTracks <= bin6Max);
    bool multBIN7 = (nGlobalTracks >= bin7Min && nGlobalTracks <= bin7Max);
    bool multBIN8 = (nGlobalTracks >= bin8Min && nGlobalTracks <= bin8Max);
    bool multBIN9 = (nGlobalTracks >= bin9Min && nGlobalTracks <= bin9Max);
    bool multBIN10 = (nGlobalTracks >= bin10Min && nGlobalTracks <= bin10Max);

    // Compute event moments and evaluate topology selections.
    auto fwH = computeFWM(tracks);
    auto [H0, H1, H2, H3, H4, H5, H6, H7, H8] = fwH;

    // Use tighter moment limits for the auxiliary triangular sample.
    const bool isTRIAux = (H1 < 0.0475) && (H2 < 0.2850) && (H3 > 0.0420) && (H4 < 0.1805) && (H5 > 0.0210) && (H6 < 0.1425) && (H7 > 0.0168) && (H8 < 0.1140);

    // Apply the nominal triangular moment limits.
    const bool isTRI = (H1 < 0.05) && (H2 < 0.30) && (H3 > 0.040) && (H4 < 0.190) && (H5 > 0.020) && (H6 < 0.150) && (H7 > 0.016) && (H8 < 0.120);

    // Apply three progressively looser triangular selections.
    const bool isTRI1 = (H1 < 0.0525) && (H2 < 0.315) && (H3 > 0.0380) && (H4 < 0.1995) && (H5 > 0.0190) && (H6 < 0.1575) && (H7 > 0.0152) && (H8 < 0.1260);

    const bool isTRI2 = (H1 < 0.0550) && (H2 < 0.330) && (H3 > 0.0360) && (H4 < 0.2090) && (H5 > 0.0180) && (H6 < 0.1650) && (H7 > 0.0144) && (H8 < 0.1320);

    const bool isTRI3 = (H1 < 0.0575) && (H2 < 0.345) && (H3 > 0.0340) && (H4 < 0.2185) && (H5 > 0.0170) && (H6 < 0.1725) && (H7 > 0.0136) && (H8 < 0.1380);

    // Apply the configurable isotropic moment limits.
    const bool isISO = (H1 < h1IsoMax) && (H2 < h2IsoMax) && (H3 < h3IsoMax) && (H4 < h4IsoMax) && (H5 < h5IsoMax) && (H6 < h6IsoMax) && (H7 < h7IsoMax) && (H8 < h8IsoMax);

    // Select all, triangular, or isotropic events for correlations.
    const bool isSelectedTopology =
      (topologyMode.value == kAllEvents) ||
      (topologyMode.value == kTriEvents && isTRI) ||
      (topologyMode.value == kIsoEvents && isISO);

    // Fill inclusive event multiplicity and recoil.
    histos.fill(HIST("Multip/ITSTPC_Multiplicity"), nGlobalTracks);
    histos.fill(HIST("pTSpectra/RecoilALL"), nGlobalTracks, nRecoil);

    if (isSelectedTopology) {
      histos.fill(HIST("pTSpectra/RecoilTRI"), nGlobalTracks, nRecoil);
    }

    // Fill moments and optional moment correlations in the first broad range.
    if (multRANGE1) {

      histos.fill(HIST("FWM_BIN1/H1"), H1);
      histos.fill(HIST("FWM_BIN1/H2"), H2);
      histos.fill(HIST("FWM_BIN1/H3"), H3);
      histos.fill(HIST("FWM_BIN1/H4"), H4);
      histos.fill(HIST("FWM_BIN1/H5"), H5);
      histos.fill(HIST("FWM_BIN1/H6"), H6);
      histos.fill(HIST("FWM_BIN1/H7"), H7);
      histos.fill(HIST("FWM_BIN1/H8"), H8);
      histos.fill(HIST("FWM_BIN1/Multip_H1"), nGlobalTracks, H1);
      histos.fill(HIST("FWM_BIN1/Multip_H2"), nGlobalTracks, H2);
      histos.fill(HIST("FWM_BIN1/Multip_H3"), nGlobalTracks, H3);
      histos.fill(HIST("FWM_BIN1/Multip_H4"), nGlobalTracks, H4);
      histos.fill(HIST("FWM_BIN1/Multip_H5"), nGlobalTracks, H5);
      histos.fill(HIST("FWM_BIN1/Multip_H6"), nGlobalTracks, H6);
      histos.fill(HIST("FWM_BIN1/Multip_H7"), nGlobalTracks, H7);
      histos.fill(HIST("FWM_BIN1/Multip_H8"), nGlobalTracks, H8);
      if (doFWMCorrelations) {
        histos.fill(HIST("FWM_BIN1/H1_H2"), H1, H2);
        histos.fill(HIST("FWM_BIN1/H1_H3"), H1, H3);
        histos.fill(HIST("FWM_BIN1/H1_H4"), H1, H4);
        histos.fill(HIST("FWM_BIN1/H1_H5"), H1, H5);
        histos.fill(HIST("FWM_BIN1/H1_H6"), H1, H6);
        histos.fill(HIST("FWM_BIN1/H1_H7"), H1, H7);
        histos.fill(HIST("FWM_BIN1/H1_H8"), H1, H8);
        histos.fill(HIST("FWM_BIN1/H2_H3"), H2, H3);
        histos.fill(HIST("FWM_BIN1/H2_H4"), H2, H4);
        histos.fill(HIST("FWM_BIN1/H2_H5"), H2, H5);
        histos.fill(HIST("FWM_BIN1/H2_H6"), H2, H6);
        histos.fill(HIST("FWM_BIN1/H2_H7"), H2, H7);
        histos.fill(HIST("FWM_BIN1/H2_H8"), H2, H8);
        histos.fill(HIST("FWM_BIN1/H3_H4"), H3, H4);
        histos.fill(HIST("FWM_BIN1/H3_H5"), H3, H5);
        histos.fill(HIST("FWM_BIN1/H3_H6"), H3, H6);
        histos.fill(HIST("FWM_BIN1/H3_H7"), H3, H7);
        histos.fill(HIST("FWM_BIN1/H3_H8"), H3, H8);
        histos.fill(HIST("FWM_BIN1/H4_H5"), H4, H5);
        histos.fill(HIST("FWM_BIN1/H4_H6"), H4, H6);
        histos.fill(HIST("FWM_BIN1/H4_H7"), H4, H7);
        histos.fill(HIST("FWM_BIN1/H4_H8"), H4, H8);
        histos.fill(HIST("FWM_BIN1/H5_H6"), H5, H6);
        histos.fill(HIST("FWM_BIN1/H5_H7"), H5, H7);
        histos.fill(HIST("FWM_BIN1/H5_H8"), H5, H8);
        histos.fill(HIST("FWM_BIN1/H6_H7"), H6, H7);
        histos.fill(HIST("FWM_BIN1/H6_H8"), H6, H8);
        histos.fill(HIST("FWM_BIN1/H7_H8"), H7, H8);
      }
    }

    // Fill moments and optional moment correlations in the second broad range.
    if (multRANGE2) {

      histos.fill(HIST("FWM_BIN2/H1"), H1);
      histos.fill(HIST("FWM_BIN2/H2"), H2);
      histos.fill(HIST("FWM_BIN2/H3"), H3);
      histos.fill(HIST("FWM_BIN2/H4"), H4);
      histos.fill(HIST("FWM_BIN2/H5"), H5);
      histos.fill(HIST("FWM_BIN2/H6"), H6);
      histos.fill(HIST("FWM_BIN2/H7"), H7);
      histos.fill(HIST("FWM_BIN2/H8"), H8);
      histos.fill(HIST("FWM_BIN2/Multip_H1"), nGlobalTracks, H1);
      histos.fill(HIST("FWM_BIN2/Multip_H2"), nGlobalTracks, H2);
      histos.fill(HIST("FWM_BIN2/Multip_H3"), nGlobalTracks, H3);
      histos.fill(HIST("FWM_BIN2/Multip_H4"), nGlobalTracks, H4);
      histos.fill(HIST("FWM_BIN2/Multip_H5"), nGlobalTracks, H5);
      histos.fill(HIST("FWM_BIN2/Multip_H6"), nGlobalTracks, H6);
      histos.fill(HIST("FWM_BIN2/Multip_H7"), nGlobalTracks, H7);
      histos.fill(HIST("FWM_BIN2/Multip_H8"), nGlobalTracks, H8);
      if (doFWMCorrelations) {
        histos.fill(HIST("FWM_BIN2/H1_H2"), H1, H2);
        histos.fill(HIST("FWM_BIN2/H1_H3"), H1, H3);
        histos.fill(HIST("FWM_BIN2/H1_H4"), H1, H4);
        histos.fill(HIST("FWM_BIN2/H1_H5"), H1, H5);
        histos.fill(HIST("FWM_BIN2/H1_H6"), H1, H6);
        histos.fill(HIST("FWM_BIN2/H1_H7"), H1, H7);
        histos.fill(HIST("FWM_BIN2/H1_H8"), H1, H8);
        histos.fill(HIST("FWM_BIN2/H2_H3"), H2, H3);
        histos.fill(HIST("FWM_BIN2/H2_H4"), H2, H4);
        histos.fill(HIST("FWM_BIN2/H2_H5"), H2, H5);
        histos.fill(HIST("FWM_BIN2/H2_H6"), H2, H6);
        histos.fill(HIST("FWM_BIN2/H2_H7"), H2, H7);
        histos.fill(HIST("FWM_BIN2/H2_H8"), H2, H8);
        histos.fill(HIST("FWM_BIN2/H3_H4"), H3, H4);
        histos.fill(HIST("FWM_BIN2/H3_H5"), H3, H5);
        histos.fill(HIST("FWM_BIN2/H3_H6"), H3, H6);
        histos.fill(HIST("FWM_BIN2/H3_H7"), H3, H7);
        histos.fill(HIST("FWM_BIN2/H3_H8"), H3, H8);
        histos.fill(HIST("FWM_BIN2/H4_H5"), H4, H5);
        histos.fill(HIST("FWM_BIN2/H4_H6"), H4, H6);
        histos.fill(HIST("FWM_BIN2/H4_H7"), H4, H7);
        histos.fill(HIST("FWM_BIN2/H4_H8"), H4, H8);
        histos.fill(HIST("FWM_BIN2/H5_H6"), H5, H6);
        histos.fill(HIST("FWM_BIN2/H5_H7"), H5, H7);
        histos.fill(HIST("FWM_BIN2/H5_H8"), H5, H8);
        histos.fill(HIST("FWM_BIN2/H6_H7"), H6, H7);
        histos.fill(HIST("FWM_BIN2/H6_H8"), H6, H8);
        histos.fill(HIST("FWM_BIN2/H7_H8"), H7, H8);
      }
    }

    // Fill moments and optional moment correlations in the third broad range.
    if (multRANGE3) {

      histos.fill(HIST("FWM_BIN3/H1"), H1);
      histos.fill(HIST("FWM_BIN3/H2"), H2);
      histos.fill(HIST("FWM_BIN3/H3"), H3);
      histos.fill(HIST("FWM_BIN3/H4"), H4);
      histos.fill(HIST("FWM_BIN3/H5"), H5);
      histos.fill(HIST("FWM_BIN3/H6"), H6);
      histos.fill(HIST("FWM_BIN3/H7"), H7);
      histos.fill(HIST("FWM_BIN3/H8"), H8);
      histos.fill(HIST("FWM_BIN3/Multip_H1"), nGlobalTracks, H1);
      histos.fill(HIST("FWM_BIN3/Multip_H2"), nGlobalTracks, H2);
      histos.fill(HIST("FWM_BIN3/Multip_H3"), nGlobalTracks, H3);
      histos.fill(HIST("FWM_BIN3/Multip_H4"), nGlobalTracks, H4);
      histos.fill(HIST("FWM_BIN3/Multip_H5"), nGlobalTracks, H5);
      histos.fill(HIST("FWM_BIN3/Multip_H6"), nGlobalTracks, H6);
      histos.fill(HIST("FWM_BIN3/Multip_H7"), nGlobalTracks, H7);
      histos.fill(HIST("FWM_BIN3/Multip_H8"), nGlobalTracks, H8);
      if (doFWMCorrelations) {
        histos.fill(HIST("FWM_BIN3/H1_H2"), H1, H2);
        histos.fill(HIST("FWM_BIN3/H1_H3"), H1, H3);
        histos.fill(HIST("FWM_BIN3/H1_H4"), H1, H4);
        histos.fill(HIST("FWM_BIN3/H1_H5"), H1, H5);
        histos.fill(HIST("FWM_BIN3/H1_H6"), H1, H6);
        histos.fill(HIST("FWM_BIN3/H1_H7"), H1, H7);
        histos.fill(HIST("FWM_BIN3/H1_H8"), H1, H8);
        histos.fill(HIST("FWM_BIN3/H2_H3"), H2, H3);
        histos.fill(HIST("FWM_BIN3/H2_H4"), H2, H4);
        histos.fill(HIST("FWM_BIN3/H2_H5"), H2, H5);
        histos.fill(HIST("FWM_BIN3/H2_H6"), H2, H6);
        histos.fill(HIST("FWM_BIN3/H2_H7"), H2, H7);
        histos.fill(HIST("FWM_BIN3/H2_H8"), H2, H8);
        histos.fill(HIST("FWM_BIN3/H3_H4"), H3, H4);
        histos.fill(HIST("FWM_BIN3/H3_H5"), H3, H5);
        histos.fill(HIST("FWM_BIN3/H3_H6"), H3, H6);
        histos.fill(HIST("FWM_BIN3/H3_H7"), H3, H7);
        histos.fill(HIST("FWM_BIN3/H3_H8"), H3, H8);
        histos.fill(HIST("FWM_BIN3/H4_H5"), H4, H5);
        histos.fill(HIST("FWM_BIN3/H4_H6"), H4, H6);
        histos.fill(HIST("FWM_BIN3/H4_H7"), H4, H7);
        histos.fill(HIST("FWM_BIN3/H4_H8"), H4, H8);
        histos.fill(HIST("FWM_BIN3/H5_H6"), H5, H6);
        histos.fill(HIST("FWM_BIN3/H5_H7"), H5, H7);
        histos.fill(HIST("FWM_BIN3/H5_H8"), H5, H8);
        histos.fill(HIST("FWM_BIN3/H6_H7"), H6, H7);
        histos.fill(HIST("FWM_BIN3/H6_H8"), H6, H8);
        histos.fill(HIST("FWM_BIN3/H7_H8"), H7, H8);
      }
    }

    // Fill inclusive and topology-selected leading-track spectra.
    for (const auto& track : tracks) {

      histos.fill(HIST("pTSpectra/pT"), track.pt());
      histos.fill(HIST("Spectra/phi"), track.phi());

      if (track.globalIndex() == triggTrk.globalIndex()) {
        continue;
      }

      if (multRANGE1) {
        histos.fill(HIST("pTSpectra/pT_LA_BIN1"), triggTrk.pt(), track.pt());
        histos.fill(HIST("pTSpectra/fpT_LA_BIN1"), triggTrk.pt(), track.pt());
      }
      if (multRANGE2) {
        histos.fill(HIST("pTSpectra/pT_LA_BIN2"), triggTrk.pt(), track.pt());
        histos.fill(HIST("pTSpectra/fpT_LA_BIN2"), triggTrk.pt(), track.pt());
      }
      if (multRANGE3) {
        histos.fill(HIST("pTSpectra/pT_LA_BIN3"), triggTrk.pt(), track.pt());
        histos.fill(HIST("pTSpectra/fpT_LA_BIN3"), triggTrk.pt(), track.pt());
      }

      if (!isSelectedTopology) {
        continue;
      }

      if (multRANGE1) {
        histos.fill(HIST("pTSpectra/pT_LA_TriBIN1"), triggTrk.pt(), track.pt());
        histos.fill(HIST("pTSpectra/fpT_LA_TriBIN1"), triggTrk.pt(), track.pt());
      }
      if (multRANGE2) {
        histos.fill(HIST("pTSpectra/pT_LA_TriBIN2"), triggTrk.pt(), track.pt());
        histos.fill(HIST("pTSpectra/fpT_LA_TriBIN2"), triggTrk.pt(), track.pt());
      }
      if (multRANGE3) {
        histos.fill(HIST("pTSpectra/pT_LA_TriBIN3"), triggTrk.pt(), track.pt());
        histos.fill(HIST("pTSpectra/fpT_LA_TriBIN3"), triggTrk.pt(), track.pt());
      }
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

        const bool onPeaks = (deltaPhi > -0.4 && deltaPhi < 0.4) ||
                             (deltaPhi > 1.6 && deltaPhi < 2.4) ||
                             (deltaPhi > 3.6 && deltaPhi < 4.4);

        const bool onValleys = (deltaPhi > -2.6 && deltaPhi < -1.4) ||
                               (deltaPhi > 0.6 && deltaPhi < 1.4) ||
                               (deltaPhi > 2.6 && deltaPhi < 3.4);

        // Apply the configured pair selection to each interval.
        const bool isLeadAssocPair = passesLeadAssoc(leadAssocMode.value, foundPart, track1, track2, triggTrk, triggTrkPart);

        if (multRANGE1) {
          fillSameEventHist(1, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }

        if (multRANGE2) {
          fillSameEventHist(2, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }

        if (multRANGE3) {
          fillSameEventHist(3, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }

        if (multBIN4) {
          fillSameEventHist(4, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }
        if (multBIN5) {
          fillSameEventHist(5, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }
        if (multBIN6) {
          fillSameEventHist(6, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }
        if (multBIN7) {
          fillSameEventHist(7, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }
        if (multBIN8) {
          fillSameEventHist(8, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }
        if (multBIN9) {
          fillSameEventHist(9, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
        }
        if (multBIN10) {
          fillSameEventHist(10, track2, deltaPhi, deltaEta, onPeaks, onValleys, isLeadAssocPair);
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
    fillHistoMult(1, range1Min, range1Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(2, range2Min, range2Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(3, range3Min, range3Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(4, bin4Min, bin4Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(5, bin5Min, bin5Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(6, bin6Min, bin6Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(7, bin7Min, bin7Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(8, bin8Min, bin8Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(9, bin9Min, bin9Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
    fillHistoMult(10, bin10Min, bin10Max, isSelectedTopology, isTRI, isTRI1, isTRI2, isTRI3, isTRIAux, isISO, tracks);
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

  SliceCache cache_bin1, cache_bin2, cache_bin3, cache_bin4, cache_bin5, cache_bin6;
  SliceCache cache_bin7, cache_bin8, cache_bin9, cache_bin10;
  Preslice<aod::TriangleTracks> triCol = aod::track_tri::triangleCollisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr std::array<std::string_view, 10> mixedBinNames = {
    "BIN1", "BIN2", "BIN3", "BIN4", "BIN5", "BIN6",
    "BIN7", "BIN8", "BIN9", "BIN10"};

  Configurable<int> nBin{"nBin", 100, "No of bis per axis"};
  ConfigurableAxis axisVtx{"axisVtx", {VARIABLE_WIDTH, -10, 10}, "vertex axis for mixed event histograms"};

  // Use half-integer edges for inclusive integer multiplicity bins.
  ConfigurableAxis axisMult1{"axisMult1", {VARIABLE_WIDTH, 39.5, 85.5}, "axisMult1"};
  ConfigurableAxis axisMult2{"axisMult2", {VARIABLE_WIDTH, 85.5, 140.5}, "axisMult2"};
  ConfigurableAxis axisMult3{"axisMult3", {VARIABLE_WIDTH, 140.5, 280.5}, "axisMult3"};
  ConfigurableAxis axisMult4{"axisMult4", {VARIABLE_WIDTH, 39.5, 60.5}, "axisMult4"};
  ConfigurableAxis axisMult5{"axisMult5", {VARIABLE_WIDTH, 60.5, 80.5}, "axisMult5"};
  ConfigurableAxis axisMult6{"axisMult6", {VARIABLE_WIDTH, 80.5, 100.5}, "axisMult6"};
  ConfigurableAxis axisMult7{"axisMult7", {VARIABLE_WIDTH, 100.5, 120.5}, "axisMult7"};
  ConfigurableAxis axisMult8{"axisMult8", {VARIABLE_WIDTH, 120.5, 140.5}, "axisMult8"};
  ConfigurableAxis axisMult9{"axisMult9", {VARIABLE_WIDTH, 140.5, 160.5}, "axisMult9"};
  ConfigurableAxis axisMult10{"axisMult10", {VARIABLE_WIDTH, 160.5, 280.5}, "axisMult10"};

  Configurable<bool> doMixedEvent{"doMixedEvent", true, "Fill mixed-event TRI correlations"};
  Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of mixed events per trigger event"};

  // Choose the pair selection for mixed-event correlations.
  Configurable<int> leadAssocMode{"leadAssocMode", 0,
                                  "0: ALL i-j pairs (pT<3); 1: ALL i-j pairs (pT < 3, Eta constrain) 2: pTLead in (2,3) & pTAssoc < 2; 3: pTAssoc < pTLead"};

  // Register mixed-event pair counts and angular distributions.
  void init(InitContext const&)
  {
    const AxisSpec axisDeltaEta{32, -1.6, +1.6, "#Delta#eta"};
    const AxisSpec axisDeltaPhi{72, -constants::math::PIHalf, +3 * constants::math::PIHalf, "#Delta#phi"};
    const AxisSpec axisEvent{10, 0.5, 10.5, "#Multip. Bin", "EventAxis"};

    // Label mixed-event pair counts by multiplicity interval.
    histos.add("MixingPairs/MixingPairsNoLA", "MixingPairsNoLA", kTH1D, {axisEvent}, false);
    auto hstatLA = histos.get<TH1>(HIST("MixingPairs/MixingPairsNoLA"));
    auto* hLA = hstatLA->GetXaxis();
    hLA->SetBinLabel(1, "bPairsLA_BIN1");
    hLA->SetBinLabel(2, "bPairsLA_BIN2");
    hLA->SetBinLabel(3, "bPairsLA_BIN3");
    hLA->SetBinLabel(4, "bPairsLA_BIN4");
    hLA->SetBinLabel(5, "bPairsLA_BIN5");
    hLA->SetBinLabel(6, "bPairsLA_BIN6");
    hLA->SetBinLabel(7, "bPairsLA_BIN7");
    hLA->SetBinLabel(8, "bPairsLA_BIN8");
    hLA->SetBinLabel(9, "bPairsLA_BIN9");
    hLA->SetBinLabel(10, "bPairsLA_BIN10");

    // Register mixed-event correlation histograms for all intervals.
    for (int i = 1; i <= 10; i++) {
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
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN1{binningOnPosMultip1, nMixedEvents, -1, &cache_bin1};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN2{binningOnPosMultip2, nMixedEvents, -1, &cache_bin2};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN3{binningOnPosMultip3, nMixedEvents, -1, &cache_bin3};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN4{binningOnPosMultip4, nMixedEvents, -1, &cache_bin4};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN5{binningOnPosMultip5, nMixedEvents, -1, &cache_bin5};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN6{binningOnPosMultip6, nMixedEvents, -1, &cache_bin6};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN7{binningOnPosMultip7, nMixedEvents, -1, &cache_bin7};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN8{binningOnPosMultip8, nMixedEvents, -1, &cache_bin8};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN9{binningOnPosMultip9, nMixedEvents, -1, &cache_bin9};
  SameKindPair<aod::TriangleCollisions, aod::TriangleTracks, BinningType> pairBIN10{binningOnPosMultip10, nMixedEvents, -1, &cache_bin10};

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
          constexpr int idx = i.value;
          if (binIndex == idx + 1) {
            histos.fill(HIST("Corelations/bProjTRI_") + HIST(mixedBinNames[idx]), deltaPhi);
            histos.fill(HIST("Corelations/bTRI_") + HIST(mixedBinNames[idx]), deltaPhi, deltaEta);
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

    // Process all ten independent mixing pools.
    fillMixedEventBin(pairBIN1, 1);
    fillMixedEventBin(pairBIN2, 2);
    fillMixedEventBin(pairBIN3, 3);
    fillMixedEventBin(pairBIN4, 4);
    fillMixedEventBin(pairBIN5, 5);
    fillMixedEventBin(pairBIN6, 6);
    fillMixedEventBin(pairBIN7, 7);
    fillMixedEventBin(pairBIN8, 8);
    fillMixedEventBin(pairBIN9, 9);
    fillMixedEventBin(pairBIN10, 10);
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
