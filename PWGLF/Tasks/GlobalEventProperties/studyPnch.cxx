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
/// \file studyPnch.cxx
///
/// \brief task for analysis of charged-particle multiplicity distributions
/// \author Abhi Modak (abhi.modak@cern.ch), Lucas José (lucas.jose.franco.da.silva@cern.ch)
/// \since September 10, 2025

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

using ColDataTable = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using FilTrackDataTable = soa::Filtered<TrackDataTable>;
using ColMCTrueTable = aod::McCollisions;
using TrackMCTrueTable = aod::McParticles;
using ColMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults>>;
using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>;
using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;

static constexpr TrackSelectionFlags::flagtype TrackSelectionIts =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;
static constexpr TrackSelectionFlags::flagtype TrackSelectionTpc =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDca =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDcaxyOnly =
  TrackSelectionFlags::kDCAxy;

AxisSpec axisEvent{10, 0.5, 10.5, "#Event", "EventAxis"};
AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};
AxisSpec axisEta{40, -2, 2, "#eta", "EtaAxis"};
AxisSpec axisPhi{629, 0, o2::constants::math::TwoPI, "#phi"};
AxisSpec axisCollSel{5, 0.5, 5.5, "#Event", "CollSelAxis"};
auto static constexpr kMinCharge = 3.f;
auto static constexpr kMinPtCut = 0.1f;

struct StudyPnch {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  Configurable<float> etaRange{"etaRange", 1.0f, "Eta range to consider"};
  Configurable<float> vtxRange{"vtxRange", 10.0f, "Vertex Z range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  Configurable<float> extraphicut1{"extraphicut1", 3.07666f, "Extra Phi cut 1"};
  Configurable<float> extraphicut2{"extraphicut2", 3.12661f, "Extra Phi cut 2"};
  Configurable<float> extraphicut3{"extraphicut3", 0.03f, "Extra Phi cut 3"};
  Configurable<float> extraphicut4{"extraphicut4", 6.253f, "Extra Phi cut 4"};
  ConfigurableAxis multHistBin{"multHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis pvHistBin{"pvHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis fv0aMultHistBin{"fv0aMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ft0aMultHistBin{"ft0aMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ft0cMultHistBin{"ft0cMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ptHistBin{"ptHistBin", {200, 0., 20.}, ""};
  ConfigurableAxis binsDCA{"binsDCA", {500, -10.0f, 10.0f}, ""};
  ConfigurableAxis countNumberTracks{"countNumberTracks", {10, -0.5, 9.5}, ""};

  Configurable<bool> isApplyTFcut{"isApplyTFcut", true, "Enable TimeFrameBorder cut"};
  Configurable<bool> isApplyITSROcut{"isApplyITSROcut", true, "Enable ITS ReadOutFrameBorder cut"};
  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyInelgt0{"isApplyInelgt0", false, "Enable INEL > 0 condition"};
  Configurable<bool> isApplyExtraPhiCut{"isApplyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<bool> isApplyTVX{"isApplyTVX", false, "Enable TVX trigger sel"};
  Configurable<bool> isApplyCheckID{"isApplyCheckID", true, "Select Tracks evaluating Collision ID"};
  Configurable<bool> isApplyDuplicatedTrack{"isApplyDuplicatedTrack", true, "Select tracks that are not duplicated"};
  Configurable<bool> isApplyPhiSelection{"isApplyPhiSelection", false, "Select tracks in specific phi range"};
  Configurable<float> minPhi{"minPhi", 0.f, "Minimum phi value for track selection"};
  Configurable<float> maxPhi{"maxPhi", 6.283185f, "Maximum phi value for track selection"};
  Configurable<bool> ispTincrease{"ispTincrease", false, "Varies low pT particles by a conservative amount of +100%"};
  Configurable<bool> ispTdecrease{"ispTdecrease", false, "Varies low pT particles by a conservative amount of -50%"};
  Configurable<bool> isApplyStrangenessSysUncert{"isApplyStrangenessSysUncert", false, "Enable the evaluation of systematics due to strange particle contribution"};

  void init(InitContext const&)
  {
    AxisSpec axisMult = {multHistBin, "Mult", "MultAxis"};
    AxisSpec axisPV = {pvHistBin, "PV", "PVAxis"};
    AxisSpec axisFv0aMult = {fv0aMultHistBin, "fv0a", "FV0AMultAxis"};
    AxisSpec axisFt0aMult = {ft0aMultHistBin, "ft0a", "FT0AMultAxis"};
    AxisSpec axisFt0cMult = {ft0cMultHistBin, "ft0c", "FT0CMultAxis"};
    AxisSpec axisPt = {ptHistBin, "pT", "pTAxis"};
    AxisSpec axisCountNumberTracks = {countNumberTracks, "Count", "CountAxis"};
    AxisSpec dcaAxis = {binsDCA, "DCA vs PV"};

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);

    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "kIsTriggerTVX");
    x->SetBinLabel(3, "kNoTimeFrameBorder");
    x->SetBinLabel(4, "kNoITSROFrameBorder");
    x->SetBinLabel(5, "kNoSameBunchPileup"); // reject collisions in case of pileup with another collision in the same foundBC
    x->SetBinLabel(6, "INEL > 0");
    x->SetBinLabel(7, "|vz| < 10");

    histos.add("SelCollsHist", "SelCollsHist", kTH1D, {axisCollSel}, false);
    auto hstat_colls = histos.get<TH1>(HIST("SelCollsHist"));
    auto* xColls = hstat_colls->GetXaxis();
    xColls->SetBinLabel(1, "All collisions");
    xColls->SetBinLabel(2, "Best Collision Selection");
    xColls->SetBinLabel(3, "Has MC Collision Selection");

    if (doprocessData || doprocessCorrelation || doprocessMonteCarlo) {
      histos.add("PhiVsEtaHist", "PhiVsEtaHist", kTH2F, {axisPhi, axisEta}, false);
      histos.add("EtaHist", "EtaHist", kTH1D, {axisEta}, false);
      histos.add("PhiHist", "PhiHist", kTH1D, {axisPhi}, false);
      histos.add("hdcaxy", "dca to pv in the xy plane", kTH1D, {dcaAxis}, false);
      histos.add("hdcaz", "dca to pv in the z axis", kTH1D, {dcaAxis}, false);
    }
    if (doprocessData) {
      histos.add("hMultiplicityData", "hMultiplicityData", kTH1F, {axisMult}, true);
    }
    if (doprocessCorrelation) {
      histos.add("GlobalMult_vs_FT0A", "GlobalMult_vs_FT0A", kTH2F, {axisMult, axisFt0aMult}, true);
      histos.add("GlobalMult_vs_FT0C", "GlobalMult_vs_FT0C", kTH2F, {axisMult, axisFt0cMult}, true);
      histos.add("GlobalMult_vs_FV0A", "GlobalMult_vs_FV0A", kTH2F, {axisMult, axisFv0aMult}, true);
      histos.add("NPVtracks_vs_FT0C", "NPVtracks_vs_FT0C", kTH2F, {axisPV, axisFt0cMult}, true);
      histos.add("NPVtracks_vs_GlobalMult", "NPVtracks_vs_GlobalMult", kTH2F, {axisPV, axisMult}, true);
    }
    if (doprocessMonteCarlo) {
      histos.add("PhiVsEtaGenHist", "PhiVsEtaGenHist", kTH2F, {axisPhi, axisEta}, false);
      histos.add("EtaGenHist", "EtaGenHist", kTH1D, {axisEta}, false);
      histos.add("PhiGenHist", "PhiGenHist", kTH1D, {axisPhi}, false);
      histos.add("hMultiplicityMCrec", "hMultiplicityMCrec", kTH1F, {axisMult}, true);
      histos.add("hMultiplicityMCgen", "hMultiplicityMCgen", kTH1F, {axisMult}, true);
      histos.add("hResponseMatrix", "hResponseMatrix", kTH2F, {axisMult, axisMult}, true);
      histos.add("hCountNTracks", "hCountNTracks", kTH1F, {axisCountNumberTracks}, true);
    }
    if (ispTincrease || ispTdecrease) {
      histos.add("hMultiplicityMCgenPtCut", "hMultiplicityMCgenPtCut", kTH1F, {axisMult}, true);
      histos.add("hResponseMatrixPtCut", "hResponseMatrixPtCut", kTH2F, {axisMult, axisMult}, true);
    }
    if (isApplyStrangenessSysUncert) {
      histos.add("hMultiplicityMCStangeDecay", "hMultiplicityMCStangeDecay", kTH1F, {axisMult}, true);
      histos.add("hMultiplicityMCSubtractionSDecay", "hMultiplicityMCSubtractionSDecay", kTH1F, {axisMult}, true);
      histos.add("hResponseMatrixStrangeDecay", "hResponseMatrixStrangeDecay", kTH2F, {axisMult, axisMult}, true);
      histos.add("hResponseMatrixSubtractionSDecay", "hResponseMatrixSubtractionSDecay", kTH2F, {axisMult, axisMult}, true);
    }
    if (doprocessEvtLossSigLossMC) {
      histos.add("MCEventHist", "MCEventHist", kTH1F, {axisEvent}, false);
      auto hstat = histos.get<TH1>(HIST("MCEventHist"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All MC events");
      x->SetBinLabel(2, "MC events with atleast one reco event");
      histos.add("hMultiplicityMCgenAll", "hMultiplicityMCgenAll", kTH1F, {axisMult}, true);
      histos.add("hMultiplicityMCgenSel", "hMultiplicityMCgenSel", kTH1F, {axisMult}, true);
    }
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    histos.fill(HIST("EventHist"), 1);
    if (!col.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 2);
    if (isApplyTFcut && !col.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 3);
    if (isApplyITSROcut && !col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 4);
    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 5);
    if (isApplyInelgt0 && !col.isInelGt0()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 6);
    if (std::abs(col.posZ()) >= vtxRange) {
      return false;
    }
    histos.fill(HIST("EventHist"), 7);
    histos.fill(HIST("VtxZHist"), col.posZ());
    return true;
  }

  template <typename CheckTrack>
  bool isTrackSelected(CheckTrack const& track)
  {
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (isApplyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    return true;
  }

  template <typename CheckGenTrack>
  bool isGenTrackSelected(CheckGenTrack const& track)
  {
    if (!track.isPhysicalPrimary()) {
      return false;
    }
    if (!track.producedByGenerator()) {
      return false;
    }
    auto pdgTrack = pdg->GetParticle(track.pdgCode());
    if (pdgTrack == nullptr) {
      return false;
    }
    if (std::abs(pdgTrack->Charge()) < kMinCharge) {
      return false;
    }
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (isApplyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    return true;
  }

  template <typename countTrk>
  int countNTracks(countTrk const& tracks)
  {
    auto nTrk = 0;
    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      if (isApplyPhiSelection && (track.phi() < minPhi || track.phi() > maxPhi)) {
        continue;
      }
      histos.fill(HIST("hdcaxy"), track.dcaXY());
      histos.fill(HIST("hdcaz"), track.dcaZ());
      histos.fill(HIST("EtaHist"), track.eta());
      histos.fill(HIST("PhiHist"), track.phi());
      histos.fill(HIST("PhiVsEtaHist"), track.phi(), track.eta());
      nTrk++;
    }
    return nTrk;
  }

  template <typename countTrk, typename McColType>
  int countGenTracks(countTrk const& tracks, McColType const& McCol)
  {
    auto nTrk = 0;
    for (const auto& track : tracks) {
      if (!isGenTrackSelected(track)) {
        continue;
      }
      if (track.mcCollisionId() != McCol.mcCollisionId()) {
        continue;
      }
      if (isApplyPhiSelection && (track.phi() < minPhi || track.phi() > maxPhi)) {
        continue;
      }
      histos.fill(HIST("EtaGenHist"), track.eta());
      histos.fill(HIST("PhiGenHist"), track.phi());
      histos.fill(HIST("PhiVsEtaGenHist"), track.phi(), track.eta());
      nTrk++;
    }
    return nTrk;
  }

  template <typename countTrk, typename McColType>
  int countNTracksMcCol(countTrk const& tracks, McColType const& McCol)
  {
    auto nTrk = 0;
    std::vector<int> mcRecIDs;
    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      if (track.has_mcParticle()) {
        auto particle = track.mcParticle();
        if (isApplyCheckID && particle.mcCollisionId() != McCol.mcCollisionId()) {
          continue;
        }
        if (isApplyDuplicatedTrack && find(mcRecIDs.begin(), mcRecIDs.end(), particle.globalIndex()) != mcRecIDs.end()) {
          continue;
        }
        mcRecIDs.push_back(particle.globalIndex());
        if (isApplyPhiSelection && (track.phi() < minPhi || track.phi() > maxPhi)) {
          continue;
        }
        nTrk++;
      }
      histos.fill(HIST("hdcaxy"), track.dcaXY());
      histos.fill(HIST("hdcaz"), track.dcaZ());
      histos.fill(HIST("EtaHist"), track.eta());
      histos.fill(HIST("PhiHist"), track.phi());
      histos.fill(HIST("PhiVsEtaHist"), track.phi(), track.eta());
    }
    return nTrk;
  }

  template <typename countTrk, typename McColType>
  int countStrangeTracksMcCol(countTrk const& tracks, McColType const& McCol)
  {
    auto nTrk_strange = 0;
    std::vector<int> mcRecIDs;
    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      if (track.has_mcParticle()) {
        auto particle = track.mcParticle();
        if (isApplyCheckID && particle.mcCollisionId() != McCol.mcCollisionId()) {
          continue;
        }
        if (isApplyDuplicatedTrack && find(mcRecIDs.begin(), mcRecIDs.end(), particle.globalIndex()) != mcRecIDs.end()) {
          continue;
        }
        mcRecIDs.push_back(particle.globalIndex());
        if (particle.has_mothers()) {
          auto mcMother = particle.template mothers_as<aod::McParticles>().front();
          if (mcMother.pdgCode() == PDG_t::kK0Short || std::abs(mcMother.pdgCode()) == PDG_t::kLambda0) {
            nTrk_strange++;
          }
        }
      }
    }
    return nTrk_strange;
  }

  template <typename countTrk, typename McColType>
  int countTracksPtCut(countTrk const& tracks, McColType const& McCol)
  {
    auto nTrk_lowpT = 0;
    auto nTrk_highpT = 0;
    auto nTrk = 0;
    for (const auto& track : tracks) {
      if (!isGenTrackSelected(track)) {
        continue;
      }
      if (track.mcCollisionId() != McCol.mcCollisionId()) {
        continue;
      }
      // Evaluate low pT extrapolation
      if (track.pt() < kMinPtCut) {
        // nTrk_lowpT++;
        if (ispTincrease) {
          nTrk_lowpT += 2 - 10 * track.pt();
        }
        if (ispTdecrease) {
          nTrk_lowpT += 0.5 + 5 * track.pt();
        }
      } else {
        nTrk_highpT++;
      }
    }
    nTrk = nTrk_lowpT + nTrk_highpT;
    return nTrk;
  }

  Filter fTrackSelectionITS = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                              ncheckbit(aod::track::trackCutFlag, TrackSelectionIts);
  Filter fTrackSelectionTPC = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true);
  Filter fTrackSelectionDCA = ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));

  void processData(ColDataTable::iterator const& cols, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(cols)) {
      return;
    }
    auto mult = countNTracks(tracks);
    if (mult > 0) {
      histos.fill(HIST("hMultiplicityData"), mult);
    }
  }

  void processCorrelation(ColDataTable::iterator const& cols, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(cols)) {
      return;
    }
    auto mult = countNTracks(tracks);
    histos.fill(HIST("GlobalMult_vs_FT0A"), mult, cols.multFT0A());
    histos.fill(HIST("GlobalMult_vs_FT0C"), mult, cols.multFT0C());
    histos.fill(HIST("GlobalMult_vs_FV0A"), mult, cols.multFV0A());
    histos.fill(HIST("NPVtracks_vs_FT0C"), cols.multNTracksPV(), cols.multFT0C());
    histos.fill(HIST("NPVtracks_vs_GlobalMult"), cols.multNTracksPV(), mult);
  }

  void processMonteCarlo(soa::Join<aod::McCollisions, aod::McCollsExtra>::iterator const& mcCollision, ColMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      histos.fill(HIST("SelCollsHist"), 1);
      // Evaluation of reconstructed collisions with more than 1 contributor
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      histos.fill(HIST("SelCollsHist"), 2);
      if (!RecCol.has_mcCollision()) {
        continue;
      }
      histos.fill(HIST("SelCollsHist"), 3);
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      auto multrec = countNTracksMcCol(recTracksPart, RecCol);
      if (multrec > 0) {
        histos.fill(HIST("hMultiplicityMCrec"), multrec);
      }
      auto multgen = countGenTracks(GenParticles, RecCol);
      if (multgen > 0 && multrec > 0) {
        histos.fill(HIST("hMultiplicityMCgen"), multgen);
        histos.fill(HIST("hResponseMatrix"), multrec, multgen);
      }
      if (ispTincrease || ispTdecrease) {
        auto nTrkPtCut = countTracksPtCut(GenParticles, RecCol);
        if (nTrkPtCut > 0) {
          histos.fill(HIST("hMultiplicityMCgenPtCut"), nTrkPtCut);
          histos.fill(HIST("hResponseMatrixPtCut"), multrec, nTrkPtCut);
        }
      }
      if (isApplyStrangenessSysUncert) {
        auto nTrk_strange = countStrangeTracksMcCol(recTracksPart, RecCol);
        auto nSubtract_strange = multrec - nTrk_strange;
        if (multrec > 0) {
          histos.fill(HIST("hMultiplicityMCStangeDecay"), nTrk_strange);
          histos.fill(HIST("hMultiplicityMCSubtractionSDecay"), nSubtract_strange);
          histos.fill(HIST("hResponseMatrixStrangeDecay"), nTrk_strange, multgen);
          histos.fill(HIST("hResponseMatrixSubtractionSDecay"), nSubtract_strange, multgen);
        }
      }
    }
  }

  void processEvtLossSigLossMC(soa::Join<ColMCTrueTable, aod::MultMCExtras>::iterator const& mcCollision, ColMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles)
  {
    if (isApplyInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }
    if (isApplyTVX && !(mcCollision.multMCFT0C() > 0 && mcCollision.multMCFT0A() > 0)) {
      return;
    }
    if (std::abs(mcCollision.posZ()) >= vtxRange) {
      return;
    }
    // All generated events
    histos.fill(HIST("MCEventHist"), 1);
    auto nTrk_multAll = 0;
    for (const auto& GenParticle : GenParticles) {
      if (!isGenTrackSelected(GenParticle)) {
        continue;
      }
      nTrk_multAll++;
    }
    if (nTrk_multAll > 0) {
      histos.fill(HIST("hMultiplicityMCgenAll"), nTrk_multAll);
    }

    bool atLeastOne = false;
    auto numcontributors = -999;
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.numContrib() <= numcontributors) {
        continue;
      } else {
        numcontributors = RecCol.numContrib();
      }
      atLeastOne = true;
    }

    if (atLeastOne) {
      histos.fill(HIST("MCEventHist"), 2);
      auto nTrk_multSel = 0;
      for (const auto& GenParticle : GenParticles) {
        if (!isGenTrackSelected(GenParticle)) {
          continue;
        }
        nTrk_multSel++;
      }
      if (nTrk_multSel > 0) {
        histos.fill(HIST("hMultiplicityMCgenSel"), nTrk_multSel);
      }
    }
  }

  PROCESS_SWITCH(StudyPnch, processData, "process data CentFT0C", false);
  PROCESS_SWITCH(StudyPnch, processCorrelation, "do correlation study in data", false);
  PROCESS_SWITCH(StudyPnch, processMonteCarlo, "process MC CentFT0C", false);
  PROCESS_SWITCH(StudyPnch, processEvtLossSigLossMC, "process Signal Loss, Event Loss", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<StudyPnch>(cfgc)};
}
