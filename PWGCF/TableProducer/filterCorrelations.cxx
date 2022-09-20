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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "MathUtils/detail/TypeTruncation.h"

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

#include <TH3F.h>
#include <TDatabasePDG.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::math_utils::detail;

#define FLOAT_PRECISION 0xFFFFFFF0

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FilterCF {
  Service<TDatabasePDG> pdg;

  // Configuration
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 7.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutMCPt, float, 0.5f, "Minimal pT for particles (WARNING only for multiplicity estimate)")
  O2_DEFINE_CONFIGURABLE(cfgCutMCEta, float, 0.8f, "Eta range for particles (WARNING only for multiplicity estimate)")

  // Filters and input definitions
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = aod::cent::centRun2V0M >= 0.0f && aod::cent::centRun2V0M <= 100.0f;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt);
  Filter collisionVertexTypeFilter = (aod::collision::flags & (uint16_t)aod::collision::CollisionFlagsRun2::Run2VertexerTracks) == (uint16_t)aod::collision::CollisionFlagsRun2::Run2VertexerTracks;
  Filter trackSelection = (requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true);

  // TODO cannot be used yet as there is no index rewriting below for the tracks
  // Filter mcParticleFilter = (nabs(aod::mcparticle::eta) < cfgCutMCEta) && (aod::mcparticle::pt > cfgCutMCPt);

  OutputObj<TH3F> yields{TH3F("yields", "centrality vs pT vs eta", 100, 0, 100, 40, 0, 20, 100, -2, 2)};
  OutputObj<TH3F> etaphi{TH3F("etaphi", "centrality vs eta vs phi", 100, 0, 100, 100, -2, 2, 200, 0, 2 * M_PI)};

  Produces<aod::CFCollisions> outputCollisions;
  Produces<aod::CFTracks> outputTracks;

  Produces<aod::CFMcCollisions> outputMcCollisions;
  Produces<aod::CFMcParticles> outputMcParticles;

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
  {
    LOGF(info, "processData: Tracks for collision: %d | Vertex: %.1f (%d) | INT7: %d | V0M: %.1f", tracks.size(), collision.posZ(), collision.flags(), collision.sel7(), collision.centRun2V0M());

    if (!collision.alias()[kINT7] || !collision.sel7()) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    outputCollisions(-1, bc.runNumber(), collision.posZ(), collision.centRun2V0M(), bc.timestamp());

    for (auto& track : tracks) {
      uint8_t trackType = 0;
      if (track.isGlobalTrack()) {
        trackType = 1;
      } else if (track.isGlobalTrackSDD()) {
        trackType = 2;
      }

      outputTracks(outputCollisions.lastIndex(), -1, track.pt(), track.eta(), track.phi(), track.sign(), trackType);

      yields->Fill(collision.centRun2V0M(), track.pt(), track.eta());
      etaphi->Fill(collision.centRun2V0M(), track.eta(), track.phi());
    }
  }
  PROCESS_SWITCH(FilterCF, processData, "Process data", true);

  void processMC1(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TrackSelection>> const& tracks)
  {
    LOGF(info, "processMC1: Tracks for collision: %d | Vertex: %.1f (%d) | INT7: %d", tracks.size(), collision.posZ(), collision.flags(), collision.sel7());

    if (!collision.alias()[kINT7] || !collision.sel7()) {
      return;
    }

    float multiplicity = tracks.size();

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // NOTE only works if we save all MC collisions...
    outputCollisions(collision.mcCollisionId(), bc.runNumber(), collision.posZ(), multiplicity, bc.timestamp());

    for (auto& track : tracks) {
      uint8_t trackType = 0;
      if (track.isGlobalTrack()) {
        trackType = 1;
      } else if (track.isGlobalTrackSDD()) {
        trackType = 2;
      }

      // NOTE only works if we save all MC tracks...
      outputTracks(outputCollisions.lastIndex(), track.mcParticleId(), track.pt(), track.eta(), track.phi(), track.sign(), trackType);

      yields->Fill(multiplicity, track.pt(), track.eta());
      etaphi->Fill(multiplicity, track.eta(), track.phi());
    }
  }
  PROCESS_SWITCH(FilterCF, processMC1, "Process MC: data part", false);

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  void processMC2(aod::McCollision const& collision, aod::McParticles const& particles, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TrackSelection>> const& tracks)
  {
    LOGF(info, "processMC2: Particles for MC collision: %d | Vertex: %.1f", particles.size(), collision.posZ());

    bool* reconstructed = new bool[particles.size()];
    for (int i = 0; i < particles.size(); i++) {
      reconstructed[i] = false;
    }
    if (collisions.size() > 0) {
      // TODO deal with case of more than 1 collision
      auto collision = collisions.begin();
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      LOGF(info, "  Reconstructed collision at vtx-z = %f which has %d tracks", collision.posZ(), groupedTracks.size());

      for (auto& track : groupedTracks) {
        if (track.has_mcParticle()) {
          reconstructed[track.mcParticleId() - particles.begin().globalIndex()] = true;
        }
      }
    }

    int multiplicity = 0;
    for (auto& particle : particles) {
      int8_t sign = 0;
      TParticlePDG* pdgparticle = pdg->GetParticle(particle.pdgCode());
      if (pdgparticle != nullptr) {
        sign = (pdgparticle->Charge() > 0) ? 1.0 : ((pdgparticle->Charge() < 0) ? -1.0 : 0.0);
      }
      if (particle.isPhysicalPrimary() && sign != 0 && std::abs(particle.eta()) < cfgCutMCEta && particle.pt() > cfgCutMCPt) {
        multiplicity++;
      }
      // use highest bit to flag if it is reconstructed
      uint8_t flags = particle.flags() & ~aod::cfmcparticle::kReconstructed; // clear bit in case of clashes in the future
      if (reconstructed[particle.index()]) {
        flags |= aod::cfmcparticle::kReconstructed;
      }

      // NOTE using "outputMcCollisions.lastIndex()+1" here to allow filling of outputMcCollisions *after* the loop
      outputMcParticles(outputMcCollisions.lastIndex() + 1, truncateFloatFraction(particle.pt(), FLOAT_PRECISION), truncateFloatFraction(particle.eta(), FLOAT_PRECISION),
                        truncateFloatFraction(particle.phi(), FLOAT_PRECISION), sign, particle.pdgCode(), flags);
    }

    outputMcCollisions(collision.posZ(), multiplicity);

    delete[] reconstructed;
  }
  PROCESS_SWITCH(FilterCF, processMC2, "Process MC: MC part", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FilterCF>(cfgc)};
}
