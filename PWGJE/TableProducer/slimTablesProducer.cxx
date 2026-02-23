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

/// \file slimTablesProducer.cxx
/// \brief Task to produce a reduced version of Tables for tracks, collisions, mcparticles and mccollisions.
/// \author Millot Louise <louise.millot@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/SlimTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> checkCentFT0M{"checkCentFT0M", false, "0: centFT0C as default, 1: use centFT0M estimator"};
  Configurable<float> centralityMin{"centralityMin", -999, ""};
  Configurable<float> centralityMax{"centralityMax", 999, ""};
  Configurable<float> minPt{"minPt", 0.15, "min pT to save"};
  Configurable<float> maxPt{"maxPt", 200.0, "max pT to save"};
  Configurable<float> minEta{"minEta", -0.9, "min eta to save"};
  Configurable<float> maxEta{"maxEta", 0.9, "max eta to save"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "Event selection"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events; jet-level rejection can also be applied at the jet finder level for jets only, here rejection is applied for collision and track process functions for the first time, and on jets in case it was set to false at the jet finder level"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};

  std::vector<int> eventSelectionBits;
  bool doSumw2 = false;

  void init(InitContext&)
  {
    doSumw2 = skipMBGapEvents; // true or false : storage of square erros when jet-jet

    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};

    histos.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    histos.add("h2_centrality_collisions", "event status vs. centrality;entries;centrality", {HistType::kTH2F, {centralityAxis, {4, 0.0, 4.0}}}, doSumw2);
    auto hColl = histos.get<TH1>(HIST("h_collisions"));
    hColl->GetXaxis()->SetBinLabel(1, "All");
    hColl->GetXaxis()->SetBinLabel(2, "eventSelection");

    histos.add("h_mcCollMCD_counts_weight", "MC event status;event status;weighted entries", {HistType::kTH1F, {{5, 0.0, 5.0}}});
    histos.add("h2_centrality_MCD", "mc event status vs. centrality;entries;centrality", {HistType::kTH2F, {centralityAxis, {4, 0.0, 4.0}}}, doSumw2);
    auto hMCD = histos.get<TH1>(HIST("h_mcCollMCD_counts_weight"));
    hMCD->GetXaxis()->SetBinLabel(1, "All");
    hMCD->GetXaxis()->SetBinLabel(2, "Has MC coll + eventSelection ");

    histos.add("h_mcCollMCP_counts_weight", "MC event status;event status;weighted entries", {HistType::kTH1F, {{7, 0.0, 7.0}}});
    histos.add("h2_centrality_MCP", "mc event status vs. centrality;entries;centrality", {HistType::kTH2F, {centralityAxis, {4, 0.0, 4.0}}}, doSumw2);
    auto hMCP = histos.get<TH1>(HIST("h_mcCollMCP_counts_weight"));
    hMCP->GetXaxis()->SetBinLabel(1, "All");
    hMCP->GetXaxis()->SetBinLabel(2, "ZVertex");
    hMCP->GetXaxis()->SetBinLabel(3, "Collision size");
    hMCP->GetXaxis()->SetBinLabel(4, "eventSelection");
    hMCP->GetXaxis()->SetBinLabel(5, "eventSelectionMC");

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
  }

  Produces<o2::aod::SlimCollisions> slimCollisions;
  Produces<o2::aod::SlMcCollisions> slimMcCollisions;
  Produces<o2::aod::SlimTracks> slimTracks;
  Produces<o2::aod::SlimParticles> slimParticles;

  Filter trackFilter = (aod::jtrack::pt >= minPt && aod::jtrack::pt < maxPt && aod::jtrack::eta > minEta && aod::jtrack::eta < maxEta);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut &&
                      (checkCentFT0M ? aod::jcollision::centFT0M : aod::jcollision::centFT0C) >= centralityMin &&
                      (checkCentFT0M ? aod::jcollision::centFT0M : aod::jcollision::centFT0C) < centralityMax);
  Filter mcCollisionFilter = (nabs(aod::jmccollision::posZ) < vertexZCut && aod::jmccollision::centFT0M >= centralityMin && aod::jmccollision::centFT0M < centralityMax); // no centFT0C for mccollisions, using centFT0M for both
  Filter particleCuts = (aod::jmcparticle::pt >= minPt && aod::jmcparticle::pt < maxPt && aod::jmcparticle::eta > minEta && aod::jmcparticle::eta < maxEta);

  void processData(soa::Filtered<o2::aod::JetCollisions>::iterator const& collision,
                   soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>> const& tracks)
  {
    histos.fill(HIST("h_collisions"), 0.5);
    float centrality = -1.0;
    checkCentFT0M ? centrality = collision.centFT0M() : centrality = collision.centFT0C();
    histos.fill(HIST("h2_centrality_collisions"), centrality, 0.5, 1.0);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, false, applyRCTSelections)) {
      return;
    }
    histos.fill(HIST("h_collisions"), 1.5);

    slimCollisions(collision.posZ());
    auto slimCollIndex = slimCollisions.lastIndex();
    for (const auto& track : tracks) {
      float mass = jetderiveddatautilities::mPion;
      float p = track.pt() * std::cosh(track.eta());
      float energy = std::sqrt(p * p + mass * mass);
      slimTracks(slimCollIndex, track.pt(), track.eta(), track.phi(), track.px(), track.py(), track.pz(), energy);
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processData, "process collisions and tracks for Data and MCD", false);

  void processMCD(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                  soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const&, // join the weight
                  soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>> const& tracks)
  {
    float eventWeight = collision.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>().weight();
    histos.fill(HIST("h_mcCollMCD_counts_weight"), 0.5, eventWeight);

    float centrality = -1.0;
    checkCentFT0M ? centrality = collision.centFT0M() : centrality = collision.centFT0C();
    histos.fill(HIST("h2_centrality_MCD"), centrality, 0.5, eventWeight);

    if (!collision.has_mcCollision()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) {
      return;
    }
    histos.fill(HIST("h_mcCollMCD_counts_weight"), 1.5, eventWeight);
    auto slimCollIndex = slimCollisions.lastIndex();
    slimCollisions(collision.posZ());
    for (const auto& track : tracks) {
      float mass = jetderiveddatautilities::mPion;
      float p = track.pt() * std::cosh(track.eta());
      float energy = std::sqrt(p * p + mass * mass);
      slimTracks(slimCollIndex, track.pt(), track.eta(), track.phi(), track.px(), track.py(), track.pz(), energy);
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processMCD, "process collisions and tracks for MCD", false);

  void processMCP(soa::Filtered<o2::aod::JetMcCollisions>::iterator const& mcCollision,
                  soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                  soa::Filtered<aod::JetParticles> const& particles)
  {
    float eventWeight = mcCollision.weight();
    float centrality = mcCollision.centFT0M(); // checkCentFT0M ? centrality = mccollision.centFT0M() : centrality = mccollision.centFT0C();
    histos.fill(HIST("h_mcCollMCP_counts_weight"), 0.5, eventWeight);
    histos.fill(HIST("h2_centrality_MCP"), centrality, 0.5, eventWeight);
    if (std::abs(mcCollision.posZ()) > vertexZCut) {
      return;
    }
    histos.fill(HIST("h_mcCollMCP_counts_weight"), 1.5, eventWeight);
    if (collisions.size() < 1) {
      return;
    }
    histos.fill(HIST("h_mcCollMCP_counts_weight"), 2.5, eventWeight);
    bool hasSel8Coll = false;
    for (auto const& collision : collisions) {
      if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) { // look if the rec collision associated to the mc collision passes the event selection
        hasSel8Coll = true;
      }
    }
    if (!hasSel8Coll) {
      return;
    }
    histos.fill(HIST("h_mcCollMCP_counts_weight"), 3.5, eventWeight);
    if (!jetderiveddatautilities::selectMcCollision(mcCollision, skipMBGapEvents, applyRCTSelections)) {
      return;
    }
    histos.fill(HIST("h_mcCollMCP_counts_weight"), 4.5, eventWeight);
    auto slimMcCollIndex = slimMcCollisions.lastIndex();
    slimMcCollisions(mcCollision.posZ());
    for (const auto& particle : particles) {
      slimParticles(slimMcCollIndex, particle.pt(), particle.eta(), particle.phi(), particle.px(), particle.py(), particle.pz());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processMCP, "process mccollisions and mcparticles for MCD", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc)};
}
