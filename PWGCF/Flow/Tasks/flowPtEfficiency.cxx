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
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include <iostream>
#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct flowPtEfficiency {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgTrkSelRun3ITSMatch, bool, false, "GlobalTrackRun3ITSMatching::Run3ITSall7Layers selection")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 70.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCcrossedrows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "DCAxy cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxyppPass3Enabled, bool, false, "switch of ppPass3 DCAxy pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAzPtDepEnabled, bool, false, "switch of DCAz pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgSelRunNumberEnabled, bool, false, "switch of run number selection")
  Configurable<std::vector<int>> cfgRunNumberList{"cfgRunNumberList", std::vector<int>{-1}, "runnumber list in consideration for analysis"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};

  // Filter the tracks
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>>;

  // Filter for collisions
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  // Filter for MCParticle
  Filter particleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPtMin) && (aod::mcparticle::pt < cfgCutPtMax);
  using myMcParticles = soa::Filtered<aod::McParticles>;

  // Filter for MCcollisions
  Filter mccollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  using myMcCollisions = soa::Filtered<aod::McCollisions>;

  // Additional filters for tracks
  TrackSelection myTrackSel;

  // Define the output
  HistogramRegistry registry{"registry"};

  bool isStable(int pdg)
  {
    if (abs(pdg) == 211)
      return true;
    if (abs(pdg) == 321)
      return true;
    if (abs(pdg) == 2212)
      return true;
    if (abs(pdg) == 11)
      return true;
    if (abs(pdg) == 13)
      return true;
    return false;
  }

  void init(InitContext const&)
  {
    const AxisSpec axisCounter{1, 0, +1, ""};
    // create histograms
    registry.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCRec", "Monte Carlo Reco", {HistType::kTH1D, {axisPt}});

    registry.add("mcEventCounter", "Monte Carlo Truth EventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCGen", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});

    if (cfgTrkSelRun3ITSMatch) {
      myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    } else {
      myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    }
    if (cfgCutDCAxyppPass3Enabled) {
      myTrackSel.SetMaxDcaXYPtDep([](float pt) { return 0.004f + 0.013f / pt; });
    } else {
      myTrackSel.SetMaxDcaXY(cfgCutDCAxy);
    }
    myTrackSel.SetMinNClustersTPC(cfgCutTPCclu);
    myTrackSel.SetMinNCrossedRowsTPC(cfgCutTPCcrossedrows);
    if (!cfgCutDCAzPtDepEnabled)
      myTrackSel.SetMaxDcaZ(cfgCutDCAz);
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    if (cfgCutDCAzPtDepEnabled && (track.dcaZ() > (0.004f + 0.013f / track.pt())))
      return false;
    return myTrackSel.IsSelected(track);
  }

  void processReco(myCollisions::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks, aod::McParticles const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    if (!collision.sel8())
      return;
    if (tracks.size() < 1)
      return;
    if (cfgSelRunNumberEnabled) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      int RunNumber = bc.runNumber();
      if (!std::count(cfgRunNumberList.value.begin(), cfgRunNumberList.value.end(), RunNumber))
        return;
    }
    for (const auto& track : tracks) {
      if (!trackSelected(track))
        continue;
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        if (isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("hPtMCRec"), track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(flowPtEfficiency, processReco, "process reconstructed information", true);

  void processSim(myMcCollisions::iterator const& collision, aod::BCsWithTimestamps const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, myMcParticles const& mcParticles)
  {
    if (cfgSelRunNumberEnabled) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      int RunNumber = bc.runNumber();
      if (!std::count(cfgRunNumberList.value.begin(), cfgRunNumberList.value.end(), RunNumber))
        return;
    }
    if (collisions.size() > -1) {
      registry.fill(HIST("mcEventCounter"), 0.5);
      for (const auto& mcParticle : mcParticles) {
        if (mcParticle.isPhysicalPrimary() && isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("hPtMCGen"), mcParticle.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(flowPtEfficiency, processSim, "process pure simulation information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowPtEfficiency>(cfgc)};
}
