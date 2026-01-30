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
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since September 10, 2025

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TPDGCode.h>

#include <cmath>
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
auto static constexpr kMinCharge = 3.f;

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

  Configurable<bool> isApplyTFcut{"isApplyTFcut", true, "Enable TimeFrameBorder cut"};
  Configurable<bool> isApplyITSROcut{"isApplyITSROcut", true, "Enable ITS ReadOutFrameBorder cut"};
  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyInelgt0{"isApplyInelgt0", false, "Enable INEL > 0 condition"};
  Configurable<bool> isApplyExtraPhiCut{"isApplyExtraPhiCut", false, "Enable extra phi cut"};

  void init(InitContext const&)
  {
    AxisSpec axisMult = {multHistBin, "Mult", "MultAxis"};
    AxisSpec axisPV = {pvHistBin, "PV", "PVAxis"};
    AxisSpec axisFv0aMult = {fv0aMultHistBin, "fv0a", "FV0AMultAxis"};
    AxisSpec axisFt0aMult = {ft0aMultHistBin, "ft0a", "FT0AMultAxis"};
    AxisSpec axisFt0cMult = {ft0cMultHistBin, "ft0c", "FT0CMultAxis"};
    AxisSpec axisPt = {ptHistBin, "pT", "pTAxis"};

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

    if (doprocessData || doprocessCorrelation || doprocessMonteCarlo || doprocessTreatedMonteCarlo) {
      histos.add("PhiVsEtaHist", "PhiVsEtaHist", kTH2F, {axisPhi, axisEta}, false);
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
      histos.add("hMultiplicityMCrec", "hMultiplicityMCrec", kTH1F, {axisMult}, true);
      histos.add("hMultiplicityMCgen", "hMultiplicityMCgen", kTH1F, {axisMult}, true);
      histos.add("hResponseMatrix", "hResponseMatrix", kTH2F, {axisMult, axisMult}, true);
    }
    if (doprocessTreatedMonteCarlo) {
      histos.add("hMultiplicityTreatMCrec", "hMultiplicityTreatMCrec", kTH1F, {axisMult}, true);
      histos.add("hMultiplicityTreatMCgen", "hMultiplicityTreatMCgen", kTH1F, {axisMult}, true);
      histos.add("hResponseMatrixTreat", "hResponseMatrixTreat", kTH2F, {axisMult, axisMult}, true);
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
      // Verify that the track belongs to the given MC collision
      if (track.mcCollisionId() != McCol.globalIndex()) {
        continue;
      }
    histos.fill(HIST("PhiVsEtaHist"), track.phi(), track.eta());
    nTrk++;
    }
    return nTrk;
  }

  template <typename countTrk, typename McColType>
  int countNTracksMcCol(countTrk const& tracks, McColType const& McCol)
  {
    auto nTrk = 0;
    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      // Verify that the track belongs to the given MC collision
      if (track.has_mcParticle()) {
        auto particle = track.mcParticle();
        if (particle.mcCollisionId() != McCol.mcCollisionId()) {
          continue;
        }
      }
      histos.fill(HIST("PhiVsEtaHist"), track.phi(), track.eta());
      nTrk++;
    }
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
    histos.fill(HIST("hMultiplicityData"), mult);
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

  void processMonteCarlo(ColMCTrueTable::iterator const& mcCollision, ColMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      auto multrec = countNTracksMcCol(recTracksPart, RecCol);
      histos.fill(HIST("hMultiplicityMCrec"), multrec);
      auto multgen = countGenTracks(GenParticles, mcCollision);
      histos.fill(HIST("hMultiplicityMCgen"), multgen);
      histos.fill(HIST("hResponseMatrix"), multrec, multgen);
    }
  }

  void processTreatedMonteCarlo(ColMCTrueTable::iterator const& mcCollision, ColMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    // Count generated tracks at each iterator
    auto multgen = countGenTracks(GenParticles, mcCollision);
    histos.fill(HIST("hMultiplicityTreatMCgen"), multgen);
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      // Verify that the reconstructed collision corresponds to the given MC collision
      if (RecCol.mcCollisionId() != mcCollision.globalIndex()) {
        continue;
      }
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      auto multrec = countNTracksMcCol(recTracksPart, RecCol);
      histos.fill(HIST("hMultiplicityTreatMCrec"), multrec);
      histos.fill(HIST("hResponseMatrixTreat"), multrec, multgen);
    }
  }

  void processEvtLossSigLossMC(soa::Join<ColMCTrueTable, aod::MultMCExtras>::iterator const& mcCollision, ColMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles)
  {
    if (isApplyInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }
    if (std::abs(mcCollision.posZ()) >= vtxRange) {
      return;
    }
    // All generated events
    histos.fill(HIST("MCEventHist"), 1);
    auto multAll = countGenTracks(GenParticles, mcCollision);
    histos.fill(HIST("hMultiplicityMCgenAll"), multAll);

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
      auto multSel = countGenTracks(GenParticles, mcCollision);
      histos.fill(HIST("hMultiplicityMCgenSel"), multSel);
    }
  }

  PROCESS_SWITCH(StudyPnch, processData, "process data CentFT0C", false);
  PROCESS_SWITCH(StudyPnch, processCorrelation, "do correlation study in data", false);
  PROCESS_SWITCH(StudyPnch, processMonteCarlo, "process MC CentFT0C", false);
  PROCESS_SWITCH(StudyPnch, processTreatedMonteCarlo, "process Treated MC CentFT0C", false);
  PROCESS_SWITCH(StudyPnch, processEvtLossSigLossMC, "process Signal Loss, Event Loss", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<StudyPnch>(cfgc)};
}
