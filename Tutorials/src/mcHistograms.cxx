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
/// \brief Accessing MC data and the related MC truth.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace constants::math;

// Simple access to collision
struct VertexDistribution {
  OutputObj<TH1F> vertex{TH1F("vertex", "vertex", 100, -10, 10)};

  Configurable<int> reduceOutput{"reduce-output", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};

  // loop over MC truth McCollisions
  void process(aod::McCollision const& mcCollision)
  {
    if (reduceOutput < 2) {
      LOGF(info, "MC. vtx-z = %f", mcCollision.posZ());
    }
    vertex->Fill(mcCollision.posZ());
  }
};

// Simple analysis of PhysicalPrimary particles
struct PhysicalPrimaryCharge {
  OutputObj<TH1F> charge{TH1F("charge_prim", "charge_prim", 100, -5, 5)};
  Service<o2::framework::O2DatabasePDG> pdg;

  void process(aod::McParticles const& mcParticles)
  {
    for (auto& particle : mcParticles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle) {
        continue;
      }
      charge->Fill(pdgParticle->Charge() / 3.); // note that charge comes in units of 1/3
    }
  }
};

// Grouping between MC particles and collisions
struct AccessMcData {
  OutputObj<TH1F> phiH{TH1F("phi", "phi", 100, 0., TwoPI)};
  OutputObj<TH1F> etaH{TH1F("eta", "eta", 102, -2.01, 2.01)};

  Configurable<int> reduceOutput{"reduce-output", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};

  // group according to McCollisions
  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    // access MC truth information with mcCollision() and mcParticle() methods
    if (reduceOutput < 2) {
      LOGF(info, "MC. vtx-z = %f", mcCollision.posZ());
      LOGF(info, "First: %d | Length: %d", mcParticles.begin().index(), mcParticles.size());
    }
    int count = 0;
    for (auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        phiH->Fill(mcParticle.phi());
        etaH->Fill(mcParticle.eta());
        count++;
        // Loop over mothers and daughters
        if (mcParticle.has_mothers()) {
          // Check first mother
          auto const& mother = mcParticle.mothers_first_as<aod::McParticles>();
          if (reduceOutput == 0) {
            LOGF(info, "First mother: %d has pdg code %d", mother.globalIndex(), mother.pdgCode());
          }
          // Loop over all mothers (needed for some MCs with junctions etc.)
          for (auto& m : mcParticle.mothers_as<aod::McParticles>()) {
            LOGF(debug, "M2 %d %d", mcParticle.globalIndex(), m.globalIndex());
          }
        }
        if (mcParticle.has_daughters()) {
          for (auto& d : mcParticle.daughters_as<aod::McParticles>()) {
            LOGF(debug, "D2 %d %d", mcParticle.globalIndex(), d.globalIndex());
          }
        }
      }
    }
    if (reduceOutput < 2) {
      LOGF(info, "Primaries for this collision: %d", count);
    }
  }
};

// Access from tracks to MC particle
struct AccessMcTruth {
  OutputObj<TH1F> etaDiff{TH1F("etaDiff", ";eta_{MC} - eta_{Rec}", 100, -2, 2)};
  OutputObj<TH1F> phiDiff{TH1F("phiDiff", ";phi_{MC} - phi_{Rec}", 100, -PI, PI)};

  Configurable<int> reduceOutput{"reduce-output", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};

  // group according to reconstructed Collisions
  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks,
               aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    // access MC truth information with mcCollision() and mcParticle() methods
    if (!collision.has_mcCollision()) {
      LOGF(warning, "No MC collision for this collision, skip...");
      return;
    }
    if (reduceOutput < 2) {
      LOGF(info, "vtx-z (data) = %f | vtx-z (MC) = %f", collision.posZ(), collision.mcCollision().posZ());
    }
    for (auto& track : tracks) {
      // if (track.trackType() != 0)
      //   continue;
      // if (track.labelMask() != 0)
      //   continue;
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      if (particle.isPhysicalPrimary()) {
        etaDiff->Fill(particle.eta() - track.eta());
        auto delta = particle.phi() - track.phi();
        if (delta > PI) {
          delta -= TwoPI;
        }
        if (delta < -PI) {
          delta += TwoPI;
        }
        phiDiff->Fill(delta);
      }
      // LOGF(info, "eta: %.2f %.2f \t phi: %.2f %.2f | %d", track.mcParticle().eta(), track.eta(), track.mcParticle().phi(), track.phi(), track.mcParticle().index());
    }
  }
};

// Loop over MCCollisions and get corresponding collisions (there can be more than one)
// For each of them get the corresponding tracks
// Note the use of "SmallGroups" template, that allows to handle both Run 2, where
// we have exactly 1-to-1 correspondence between collisions and mc collisions, and
// Run 3, where we can have 0, 1, or more collisions for a given mc collision
struct LoopOverMcMatched {
  OutputObj<TH1F> etaDiff{TH1F("etaDiff", ";eta_{MC} - eta_{Rec}", 100, -2, 2)};

  Configurable<int> reduceOutput{"reduce-output", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};

  using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void process(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions,
               LabeledTracks const& tracks, aod::McParticles const& mcParticles)
  {
    // access MC truth information with mcCollision() and mcParticle() methods
    if (reduceOutput < 2) {
      LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisions.size());
    }
    for (auto& collision : collisions) {
      if (reduceOutput < 2) {
        LOGF(info, "  Reconstructed collision at vtx-z = %f", collision.posZ());
      }

      // NOTE this will be replaced by a improved grouping in the future
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      if (reduceOutput < 2) {
        LOGF(info, "  which has %d tracks", groupedTracks.size());
      }
      for (auto& track : groupedTracks) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        }
        etaDiff->Fill(track.mcParticle().eta() - track.eta());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<VertexDistribution>(cfgc),
    adaptAnalysisTask<PhysicalPrimaryCharge>(cfgc),
    adaptAnalysisTask<AccessMcData>(cfgc),
    adaptAnalysisTask<AccessMcTruth>(cfgc),
    adaptAnalysisTask<LoopOverMcMatched>(cfgc)};
}
