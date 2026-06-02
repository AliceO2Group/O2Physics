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

#include <Framework/ASoA.h>
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

#include <cstdlib>
#include <string>
#include <vector>

namespace o2::aod
{
namespace slimcollision
{
DECLARE_SOA_COLUMN(Weight, weight, float);
}
DECLARE_SOA_TABLE(SlimCollisions, "AOD", "SLIMCOLLISION",
                  o2::soa::Index<>,
                  o2::aod::collision::PosZ,
                  o2::aod::collision::CollisionTime,
                  slimcollision::Weight);
using SlimCollision = SlimCollisions::iterator;
namespace slmccollision
{
DECLARE_SOA_COLUMN(McWeight, mcWeight, float);
}
DECLARE_SOA_TABLE(SlMcCollisions, "AOD", "SLMCCOLLISION",
                  o2::soa::Index<>,
                  o2::aod::mccollision::PosZ,
                  slmccollision::McWeight);
using SlMcCollision = SlMcCollisions::iterator;
namespace slimtracks
{
DECLARE_SOA_INDEX_COLUMN(SlimCollision, slimCollision);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(E, e, float);
} // namespace slimtracks
DECLARE_SOA_TABLE(SlimTracks, "AOD", "SLIMTRACK",
                  o2::soa::Index<>,
                  slimtracks::SlimCollisionId,
                  slimtracks::Px,
                  slimtracks::Py,
                  slimtracks::Pz);
using SlimTrack = SlimTracks::iterator;
namespace slimparticles
{
DECLARE_SOA_INDEX_COLUMN(SlMcCollision, slMcCollision);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(E, e, float);
} // namespace slimparticles
DECLARE_SOA_TABLE(SlimParticles, "AOD", "SLIMPARTICLE",
                  o2::soa::Index<>,
                  slimparticles::SlMcCollisionId,
                  slimparticles::Px,
                  slimparticles::Py,
                  slimparticles::Pz,
                  slimparticles::E);
using SlimParticle = SlimParticles::iterator;
} // namespace o2::aod

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
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections; other option: uniformTracks"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events; jet-level rejection can also be applied at the jet finder level for jets only, here rejection is applied for collision and track process functions for the first time, and on jets in case it was set to false at the jet finder level"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};

  std::vector<int> eventSelectionBits;
  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
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
    auto hMCD = histos.get<TH1>(HIST("h_mcCollMCD_counts_weight"));
    hMCD->GetXaxis()->SetBinLabel(1, "All");
    hMCD->GetXaxis()->SetBinLabel(2, "hasMcCollision");
    hMCD->GetXaxis()->SetBinLabel(3, "selectCollision");

    histos.add("h_mcCollMCP_counts_weight", "MC event status;event status;weighted entries", {HistType::kTH1F, {{7, 0.0, 7.0}}});
    auto hMCP = histos.get<TH1>(HIST("h_mcCollMCP_counts_weight"));
    hMCP->GetXaxis()->SetBinLabel(1, "All");
    hMCP->GetXaxis()->SetBinLabel(2, "Zvertex");
    hMCP->GetXaxis()->SetBinLabel(3, "selectMcCollision");

    histos.add("Ntracks_pT", "track pT distribution;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, -0.5, 199.5}}}, doSumw2);
    histos.add("Nparticles_pT", "particle pT distribution;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, -0.5, 199.5}}}, doSumw2);

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
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

  Preslice<aod::JetTracksMCD> perCollisionTracks = aod::jtrack::collisionId;
  Preslice<aod::JetParticles> perMcCollisionParticles = aod::jmcparticle::mcCollisionId;

  void processData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                   soa::Filtered<aod::JetTracks> const& tracks)
  {
    histos.fill(HIST("h_collisions"), 0.5);
    float centrality = -1.0;
    checkCentFT0M ? centrality = collision.centFT0M() : centrality = collision.centFT0C();
    histos.fill(HIST("h2_centrality_collisions"), centrality, 0.5, 1.0);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, false, applyRCTSelections)) {
      return;
    }
    histos.fill(HIST("h_collisions"), 1.5);
    slimCollisions(collision.posZ(), collision.collisionTime(), 1.0);
    auto slimCollIndex = slimCollisions.lastIndex();
    for (const auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      slimTracks(slimCollIndex, track.px(), track.py(), track.pz());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processData, "process collisions and tracks for data", false);

  void processMC(soa::Filtered<aod::JetMcCollisions>::iterator const& mccollision,
                 soa::SmallGroups<aod::JetCollisionsMCD> const& collisions, // join the weight
                 soa::Filtered<aod::JetTracksMCD> const& tracks,
                 soa::Filtered<aod::JetParticles> const& particles)
  {
    float eventWeightMC = mccollision.weight();
    if (collisions.size() != 1) { // skip the mccollision if it has mre than 1 associated rec collision
      return;
    }
    histos.fill(HIST("h_mcCollMCP_counts_weight"), 0.5, eventWeightMC);
    if (std::abs(mccollision.posZ()) > vertexZCut) {
      return;
    }
    histos.fill(HIST("h_mcCollMCP_counts_weight"), 1.5, eventWeightMC);
    if (!jetderiveddatautilities::selectCollision(mccollision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) {
      return;
    }
    histos.fill(HIST("h_mcCollMCP_counts_weight"), 2.5, eventWeightMC);
    for (auto const& collision : collisions) {
      float eventWeight = collision.weight();
      histos.fill(HIST("h_mcCollMCD_counts_weight"), 0.5, eventWeightMC);
      if (!collision.has_mcCollision()) {
        continue;
      }
      histos.fill(HIST("h_mcCollMCD_counts_weight"), 1.5, eventWeightMC);
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) {
        continue;
      }
      histos.fill(HIST("h_mcCollMCD_counts_weight"), 2.5, eventWeightMC);
      slimCollisions(collision.posZ(), collision.collisionTime(), eventWeight);
      auto slimCollIndex = slimCollisions.lastIndex();
      auto slicedTracks = tracks.sliceBy(perCollisionTracks, collision.globalIndex()); // tracks associated to the rec collision
      for (const auto& track : slicedTracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection))
          continue;
        histos.fill(HIST("Ntracks_pT"), track.pt(), eventWeight);
        slimTracks(slimCollIndex, track.px(), track.py(), track.pz());
      }
      slimMcCollisions(mccollision.posZ(), eventWeightMC);
      auto slimMcCollIndex = slimMcCollisions.lastIndex();
      for (const auto& particle : particles) {
        if (!particle.isPhysicalPrimary())
          continue;
        auto pdgParticle = pdgDatabase->GetParticle(particle.pdgCode());
        if (!pdgParticle)
          continue;
        if (pdgParticle->Charge() == 0) // keep charged particles, exclude neutrals
          continue;
        histos.fill(HIST("Nparticles_pT"), particle.pt(), eventWeightMC);
        slimParticles(slimMcCollIndex, particle.px(), particle.py(), particle.pz(), particle.energy());
      }
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processMC, "process collisions & tracks, MCcollisions & particles for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc)};
}
