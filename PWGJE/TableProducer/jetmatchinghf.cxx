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

/// \file jetmatchinghf.cxx
/// \brief Match jets containing the same D0s
///
/// \author Jochen Klein <jochen.klein@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetMatchingHF {
  using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using Tracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using HfCandidates = soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate, aod::HfCandProng2MCRec>;
  using McParticles = soa::Join<aod::McParticles, aod::HfCandProng2MCGen>;
  using DetectorLevelJets = soa::Join<aod::MCDetectorLevelHFJets, aod::MCDetectorLevelHFJetConstituents>;
  using ParticleLevelJets = soa::Join<aod::MCParticleLevelHFJets, aod::MCParticleLevelHFJetConstituents>;

  Produces<aod::MatchedMCParticleDetectorLevelHFJets> jetsPartToDetMatching;
  Produces<aod::MatchedMCDetectorParticleLevelHFJets> jetsDetToPartMatching;

  Preslice<ParticleLevelJets> ParticleLevelJetsPerMcCollision = aod::jet::mcCollisionId;

  void init(InitContext const&)
  {
  }

  void process(Collisions::iterator const& collision, aod::McCollisions const& mcCollisions,
               ParticleLevelJets const& jetsMC, DetectorLevelJets const& jetsRec,
               HfCandidates const& hfcandidates, Tracks const& tracks,
               McParticles const& particlesMC)
  {
    const auto jetsPL = jetsMC.sliceBy(ParticleLevelJetsPerMcCollision, collision.mcCollisionId());

    // match rec to MC
    for (const auto& jet : jetsRec) {
      LOGF(info, "jet index: %d (coll %d, pt %g, phi %g) with %d tracks, %d HF candidates",
           jet.index(), jet.collisionId(), jet.pt(), jet.phi(), jet.tracks().size(), jet.hfcandidates().size());

      const auto& cands = jet.hfcandidates_as<HfCandidates>();
      int matchedIdx = -1;
      if ((cands.front().flagMCMatchRec() & (1 << aod::hf_cand_prong2::DecayType::D0ToPiK)) == 0) {
        jetsDetToPartMatching(matchedIdx);
        continue;
      }
      for (const auto& cand : cands) {
        const auto& daughter0 = cand.index0_as<Tracks>();
        const auto& daughter1 = cand.index1_as<Tracks>();
        if (!daughter0.has_mcParticle() || !daughter1.has_mcParticle()) {
          LOGF(warning, "Encountered candidate daughter (%d or %d) without MC particle", daughter0.globalIndex(), daughter1.globalIndex());
          continue;
        }
        const auto mother0Id = daughter0.mcParticle_as<McParticles>().mothers_as<McParticles>().front().globalIndex();
        const auto mother1Id = daughter1.mcParticle_as<McParticles>().mothers_as<McParticles>().front().globalIndex();
        LOGF(debug, "MC candidate %d with prongs: %d (MC %d), %d (MC %d)", cand.globalIndex(),
             daughter0.globalIndex(), daughter0.mcParticle_as<McParticles>().globalIndex(),
             daughter1.globalIndex(), daughter1.mcParticle_as<McParticles>().globalIndex());
        LOGF(info, "MC ids of mothers: %d - %d", mother0Id, mother1Id);
        if ((mother0Id == mother1Id) &&
            std::abs(daughter0.mcParticle_as<McParticles>().mothers_as<McParticles>().front().flagMCMatchGen()) & (1 << aod::hf_cand_prong2::DecayType::D0ToPiK)) {
          LOGF(info, "D0 - looking for jet");
          for (const auto& pjet : jetsPL) {
            for (const auto& cand : pjet.hfcandidates_as<McParticles>()) {
              if (mother0Id == cand.globalIndex()) {
                matchedIdx = pjet.globalIndex();
                LOGF(info, "Found match det to part: %d (pt %g) -> %d (pt %g)",
                     jet.globalIndex(), jet.pt(), matchedIdx, pjet.pt());
              }
            }
          }
        }
      }
      jetsDetToPartMatching(matchedIdx);
    }

    // match MC to rec
    for (const auto& jet : jetsPL) {
      LOGF(info, "MC jet index: %d (coll %d) with %d tracks, %d HF candidates",
           jet.index(), jet.mcCollisionId(), jet.tracks().size(), jet.hfcandidates().size());

      int matchedIdx = -1;
      for (const auto& cand : jet.hfcandidates_as<McParticles>()) {
        const auto& daughters = cand.daughters_as<McParticles>();
        LOGF(info, "MC candidate %d with daughters %d, %d", cand.globalIndex(), daughters.iteratorAt(0).globalIndex(), daughters.iteratorAt(1).globalIndex());
        int index0 = -1, index1 = -1;
        for (const auto& track : tracks) {
          if (!track.has_mcParticle())
            continue;
          if (track.mcParticle_as<McParticles>().globalIndex() == daughters.iteratorAt(0).globalIndex() &&
              index0 < 0) {
            index0 = track.globalIndex();
            LOGF(info, "Found track for daughter 0: %d", index0);
          }
          if (track.mcParticle_as<McParticles>().globalIndex() == daughters.iteratorAt(1).globalIndex() &&
              index1 < 0) {
            index1 = track.globalIndex();
            LOGF(info, "Found track for daughter 1: %d", index1);
          }
        }
        if (index0 < 0 || index1 < 0)
          continue;
        int candIdx = 0;
        for (const auto& prong : hfcandidates) {
          LOGF(info, "checking prong %d with daughters %d-%d, %d-%d",
               prong.globalIndex(), prong.index0Id(), prong.index0_as<Tracks>().globalIndex(), prong.index1Id(), prong.index1_as<Tracks>().globalIndex());
          if ((prong.index0_as<Tracks>().globalIndex() == index0 && prong.index1_as<Tracks>().globalIndex() == index1) ||
              (prong.index0Id() == index1 && prong.index1Id() == index0)) {
            candIdx = prong.globalIndex();
            LOGF(info, "Found matching 2prong candidate: %d", candIdx);
          }
        }
        for (const auto& djet : jetsRec) {
          if (djet.hfcandidates_as<HfCandidates>().front().globalIndex() == candIdx) {
            matchedIdx = djet.globalIndex();
            LOGF(info, "Found match part to det: %d -> %d", jet.globalIndex(), matchedIdx);
          }
        }
      }
      jetsPartToDetMatching(matchedIdx);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetMatchingHF>(cfgc, TaskName{"jet-matching-hf"})};
}
