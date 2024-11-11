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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct GFWEfficiency {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "DCAxy cut for tracks")

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};

  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "centrality axis for histograms"};


  // Filter the tracks
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (nabs(aod::track::dcaXY) < cfgCutDCAxy);
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>>;
  using Colls = soa::Join<aod::Collisions, aod::CentFT0Cs>;

  // Filter for MCParticle
  Filter particleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPtMin) && (aod::mcparticle::pt < cfgCutPtMax);
  using myMcParticles = soa::Filtered<aod::McParticles>;
  using myMcCollisionsFT0Cs = soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::CentFT0Cs>;

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
    registry.add("hCenMCRec", "Monte Carlo Reco", {HistType::kTH1D, {axisCentrality}});
    registry.add("hPtMCRec05", "Monte Carlo Reco", {HistType::kTH1D, {axisPt}});
    registry.add("hPtMCRec5060", "Monte Carlo Reco", {HistType::kTH1D, {axisPt}});

    registry.add("mcEventCounter", "Monte Carlo Truth EventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCGen", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});
    registry.add("hCenMCGen", "Monte Carlo Truth", {HistType::kTH1D, {axisCentrality}});
    registry.add("hPtMCGen05", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});
    registry.add("hPtMCGen5060", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});
  }

  void processReco(Colls::iterator const& collision, myTracks const& tracks, aod::McParticles const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);
      const auto centrality = collision.centFT0C();
      registry.fill(HIST("hCenMCRec"), centrality);
    for (const auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < 70)
        continue;

      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        if (isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("hPtMCRec"), track.pt());
            if (centrality >0 && centrality<=5){
                registry.fill(HIST("hPtMCRec05"), track.pt());
            }
            if (centrality >=50 && centrality<=60){
                registry.fill(HIST("hPtMCRec5060"), track.pt());
            }
        }
      }
    }
  }
  PROCESS_SWITCH(GFWEfficiency, processReco, "process reconstructed information", true);

  void processSim(aod::McCollision const&, soa::SmallGroups<soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>> const& collisions, myMcParticles const& mcParticles, myMcCollisionsFT0Cs const& mcCollisionsFT0Cs )
  {
    if (collisions.size() > -1) {
      registry.fill(HIST("mcEventCounter"), 0.5);
        for (const auto& mcCollisionsFT0C : mcCollisionsFT0Cs) {
          registry.fill(HIST("hCenMCGen"), mcCollisionsFT0C.centFT0C());
        }

    for (const auto& mcCollisionsFT0C : mcCollisionsFT0Cs) {
        const auto centrality = mcCollisionsFT0C.centFT0C();
      for (const auto& mcParticle : mcParticles) {
          if (mcParticle.isPhysicalPrimary() && isStable(mcParticle.pdgCode())) {
              registry.fill(HIST("hPtMCGen"), mcParticle.pt());
              if (centrality >0 && centrality<=5){
                  registry.fill(HIST("hPtMCGen05"), mcParticle.pt());
              }
              if (centrality >=50 && centrality<=60){
                  registry.fill(HIST("hPtMCGen5060"), mcParticle.pt());
              }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(GFWEfficiency, processSim, "process pure simulation information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GFWEfficiency>(cfgc)};
}
