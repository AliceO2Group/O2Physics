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
#include "Framework/O2DatabasePDGPlugin.h"

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

namespace o2::aod
{
namespace cfmultiplicity
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Centrality/multiplicity value
} // namespace cfmultiplicity
DECLARE_SOA_TABLE(CFMultiplicities, "AOD", "CFMULTIPLICITY", cfmultiplicity::Multiplicity); //! Transient multiplicity table

using CFMultiplicity = CFMultiplicities::iterator;
} // namespace o2::aod

struct FilterCF {
  Service<o2::framework::O2DatabasePDG> pdg;

  // Configuration
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 7.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutMCPt, float, 0.5f, "Minimal pT for particles")
  O2_DEFINE_CONFIGURABLE(cfgCutMCEta, float, 0.8f, "Eta range for particles")
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 1, "Verbosity level (0 = major, 1 = per collision)")
  O2_DEFINE_CONFIGURABLE(cfgTrigger, int, 7, "Trigger choice: (0 = none, 7 = sel7, 8 = sel8)")
  O2_DEFINE_CONFIGURABLE(cfgCollisionFlags, uint16_t, aod::collision::CollisionFlagsRun2::Run2VertexerTracks, "Request collision flags if non-zero (0 = off, 1 = Run2VertexerTracks)")

  // Filters and input definitions
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter collisionVertexTypeFilter = (cfgCollisionFlags == 0) || ((aod::collision::flags & cfgCollisionFlags) == cfgCollisionFlags);

  // TODO how to have this in the second task? For now they are copied
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt);
  Filter trackSelection = (requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true);

  Filter mcCollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;

  OutputObj<TH3F> yields{TH3F("yields", "centrality vs pT vs eta", 100, 0, 100, 40, 0, 20, 100, -2, 2)};
  OutputObj<TH3F> etaphi{TH3F("etaphi", "centrality vs eta vs phi", 100, 0, 100, 100, -2, 2, 200, 0, 2 * M_PI)};

  Produces<aod::CFCollisions> outputCollisions;
  Produces<aod::CFTracks> outputTracks;

  Produces<aod::CFCollLabels> outputMcCollisionLabels;
  Produces<aod::CFTrackLabels> outputTrackLabels;

  Produces<aod::CFMcCollisions> outputMcCollisions;
  Produces<aod::CFMcParticles> outputMcParticles;

  template <typename TCollision>
  bool keepCollision(TCollision& collision)
  {
    if (cfgTrigger == 0) {
      return true;
    } else if (cfgTrigger == 7) {
      return collision.alias_bit(kINT7) && collision.sel7();
    } else if (cfgTrigger == 8) {
      return collision.sel8();
    }
    return false;
  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CFMultiplicities>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "processData: Tracks for collision: %d | Vertex: %.1f (%d) | INT7: %d | Multiplicity: %.1f", tracks.size(), collision.posZ(), collision.flags(), collision.sel7(), collision.multiplicity());
    }

    if (!keepCollision(collision)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    outputCollisions(bc.runNumber(), collision.posZ(), collision.multiplicity(), bc.timestamp());

    for (auto& track : tracks) {
      uint8_t trackType = 0;
      if (track.isGlobalTrack()) {
        trackType = 1;
      } else if (track.isGlobalTrackSDD()) {
        trackType = 2;
      }

      outputTracks(outputCollisions.lastIndex(), track.pt(), track.eta(), track.phi(), track.sign(), trackType);

      yields->Fill(collision.multiplicity(), track.pt(), track.eta());
      etaphi->Fill(collision.multiplicity(), track.eta(), track.phi());
    }
  }
  PROCESS_SWITCH(FilterCF, processData, "Process data", true);

  // NOTE not filtering collisions here because in that case there can be tracks referring to MC particles which are not part of the selected MC collisions
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  void processMC(aod::McCollisions const& mcCollisions, aod::McParticles const& allParticles,
                 soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CFMultiplicities> const& allCollisions,
                 soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TrackSelection>> const& tracks,
                 aod::BCsWithTimestamps const&)
  {
    bool* reconstructed = new bool[allParticles.size()];
    int* mcParticleLabels = new int[allParticles.size()];
    for (int i = 0; i < allParticles.size(); i++) {
      reconstructed[i] = false;
      mcParticleLabels[i] = -1;
    }

    // PASS 1 on collisions: check which particles are kept
    for (auto& collision : allCollisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "processMC:   Tracks for collision %d: %d | Vertex: %.1f (%d) | INT7: %d", collision.globalIndex(), groupedTracks.size(), collision.posZ(), collision.flags(), collision.sel7());
      }

      if (!keepCollision(collision)) {
        continue;
      }

      for (auto& track : groupedTracks) {
        if (track.has_mcParticle()) {
          reconstructed[track.mcParticleId()] = true;
        }
      }
    }

    for (auto& mcCollision : mcCollisions) {
      auto particles = allParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

      if (cfgVerbosity > 0) {
        LOGF(info, "processMC: Particles for MC collision %d: %d | Vertex: %.1f", mcCollision.globalIndex(), particles.size(), mcCollision.posZ());
      }

      // Store selected MC particles and MC collisions
      int multiplicity = 0;
      for (auto& particle : particles) {
        int8_t sign = 0;
        TParticlePDG* pdgparticle = pdg->GetParticle(particle.pdgCode());
        if (pdgparticle != nullptr) {
          sign = (pdgparticle->Charge() > 0) ? 1.0 : ((pdgparticle->Charge() < 0) ? -1.0 : 0.0);
        }
        bool primary = particle.isPhysicalPrimary() && sign != 0 && std::abs(particle.eta()) < cfgCutMCEta && particle.pt() > cfgCutMCPt;
        if (primary) {
          multiplicity++;
        }
        if (reconstructed[particle.globalIndex()] || primary) {
          // keep particle

          // use highest bit to flag if it is reconstructed
          uint8_t flags = particle.flags() & ~aod::cfmcparticle::kReconstructed; // clear bit in case of clashes in the future
          if (reconstructed[particle.globalIndex()]) {
            flags |= aod::cfmcparticle::kReconstructed;
          }

          // NOTE using "outputMcCollisions.lastIndex()+1" here to allow filling of outputMcCollisions *after* the loop
          outputMcParticles(outputMcCollisions.lastIndex() + 1, truncateFloatFraction(particle.pt(), FLOAT_PRECISION), truncateFloatFraction(particle.eta(), FLOAT_PRECISION),
                            truncateFloatFraction(particle.phi(), FLOAT_PRECISION), sign, particle.pdgCode(), flags);

          // relabeling array
          mcParticleLabels[particle.globalIndex()] = outputMcParticles.lastIndex();
        }
      }
      outputMcCollisions(mcCollision.posZ(), multiplicity);
    }

    // PASS 2 on collisions: store collisions and tracks
    for (auto& collision : allCollisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "processMC:   Tracks for collision %d: %d | Vertex: %.1f (%d) | INT7: %d", collision.globalIndex(), groupedTracks.size(), collision.posZ(), collision.flags(), collision.sel7());
      }

      if (!keepCollision(collision)) {
        continue;
      }

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      // NOTE works only when we store all MC collisions (as we do here)
      outputCollisions(bc.runNumber(), collision.posZ(), collision.multiplicity(), bc.timestamp());
      outputMcCollisionLabels(collision.mcCollisionId());

      for (auto& track : groupedTracks) {
        uint8_t trackType = 0;
        if (track.isGlobalTrack()) {
          trackType = 1;
        } else if (track.isGlobalTrackSDD()) {
          trackType = 2;
        }

        int mcParticleId = track.mcParticleId();
        if (mcParticleId >= 0) {
          mcParticleId = mcParticleLabels[track.mcParticleId()];
          if (mcParticleId < 0) {
            LOGP(fatal, "processMC:     Track {} is referring to a MC particle which we do not store {} {} (reco flag {})", track.index(), track.mcParticleId(), mcParticleId, reconstructed[track.mcParticleId()]);
          }
        }
        outputTracks(outputCollisions.lastIndex(),
                     truncateFloatFraction(track.pt()), truncateFloatFraction(track.eta()), truncateFloatFraction(track.phi()), track.sign(), trackType);
        outputTrackLabels(mcParticleId);

        yields->Fill(collision.multiplicity(), track.pt(), track.eta());
        etaphi->Fill(collision.multiplicity(), track.eta(), track.phi());
      }
    }

    delete[] reconstructed;
    delete[] mcParticleLabels;
  }
  PROCESS_SWITCH(FilterCF, processMC, "Process MC", false);
};

struct MultiplicitySelector {
  Produces<aod::CFMultiplicities> output;

  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")

  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt);
  Filter trackSelection = (requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true);

  void init(InitContext& context)
  {
    int enabledFunctions = 0;
    if (doprocessRun2V0M) {
      enabledFunctions++;
    }
    if (doprocessTracks) {
      enabledFunctions++;
    }

    if (enabledFunctions != 1) {
      LOGP(fatal, "{} multiplicity selectors enabled but we need exactly 1.", enabledFunctions);
    }
  }

  void processTracks(aod::Collision const&, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
  {
    output(tracks.size());
  }
  PROCESS_SWITCH(MultiplicitySelector, processTracks, "Select track count as multiplicity", false);

  void processRun2V0M(aod::CentRun2V0Ms const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centRun2V0M());
    }
  }
  PROCESS_SWITCH(MultiplicitySelector, processRun2V0M, "Select V0M centrality as multiplicity", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FilterCF>(cfgc),
    adaptAnalysisTask<MultiplicitySelector>(cfgc)};
}
