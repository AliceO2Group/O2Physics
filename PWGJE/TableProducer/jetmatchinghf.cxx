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
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename BaseJetCollection, typename TagJetCollection,
          typename BaseToTagMatchingTable, typename TagToBaseMatchingTable, typename HfCandidates>
struct JetMatchingHF {
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.4f, "Max matching distance"};

  Produces<BaseToTagMatchingTable> jetsBaseToTag;
  Produces<TagToBaseMatchingTable> jetsTagToBase;

  // preslicing jet collections, only for MC-based collection
  // Preslice<BaseJetCollection> baseJetsPerCollision = aod::jet::mcCollisionId;
  Preslice<TagJetCollection> tagJetsPerCollision = aod::jet::mcCollisionId;

  using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using Tracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using McParticles = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

  void init(InitContext const&)
  {
  }

  void process(Collisions::iterator const& collision, aod::McCollisions const& mcCollisions,
               BaseJetCollection const& jetsBase, TagJetCollection const& jetsTag,
               Tracks const& tracks, McParticles const& particlesMC,
               HfCandidates const& hfcandidates)
  {
    // slicing base jets, only for MC-based collection
    auto jetsBasePerColl = jetsBase;
    constexpr bool jetsBaseIsMC = o2::soa::is_binding_compatible_v<BaseJetCollection, std::decay_t<aod::McCollisions>>();
    if constexpr (jetsBaseIsMC) {
      LOGF(info, "slicing base jet collection by MC collision");
    } else {
      LOGF(info, "not slicing base jet collection by MC collision");
    }
    // jetsBasePerColl = jetsBase.sliceBy(baseJetsPerCollision, collision.mcCollisionId());

    // slicing tag jets, only for MC-based collection
    constexpr bool jetsTagIsMC = o2::soa::is_binding_compatible_v<TagJetCollection, std::decay_t<aod::McCollisions>>();
    // decltype(jetsTag.sliceBy(tagJetsPerCollision, collision.mcCollisionId())) jetsTagPerColl;
    if (jetsTagIsMC) {
      LOGF(info, "slicing tag jet collection by MC collision");
    } else {
      LOGF(info, "not slicing tag jet collection by MC collision");
    }
    const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, collision.mcCollisionId());

    // geometric matching
    std::vector<double> jetsBasePhi(jetsBase.size());
    std::vector<double> jetsBaseEta(jetsBase.size());
    for (auto jet : jetsBasePerColl) {
      jetsBasePhi.emplace_back(jet.phi());
      jetsBaseEta.emplace_back(jet.eta());
    }
    std::vector<double> jetsTagPhi(jetsTag.size());
    std::vector<double> jetsTagEta(jetsTag.size());
    for (auto& jet : jetsTagPerColl) {
      jetsTagPhi.emplace_back(jet.phi());
      jetsTagEta.emplace_back(jet.eta());
    }
    auto&& [baseToTagIndexMap, tagToBaseIndexMap] = JetUtilities::MatchJetsGeometrically(jetsBasePhi, jetsBaseEta, jetsTagPhi, jetsTagEta, maxMatchingDistance);

    // forward matching
    for (const auto& jet : jetsBasePerColl) {
      LOGF(info, "jet index: %d (coll %d, pt %g, phi %g) with %d tracks, %d HF candidates",
           jet.index(), jet.collisionId(), jet.pt(), jet.phi(), jet.tracks().size(), jet.hfcandidates().size());

      const auto& cands = jet.template hfcandidates_as<HfCandidates>();
      int matchedIdx = -1;
      if ((cands.front().flagMcMatchRec() & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) != 0) {
        for (const auto& cand : cands) {
          const auto& daughter0 = cand.template prong0_as<Tracks>();
          const auto& daughter1 = cand.template prong1_as<Tracks>();
          if (!daughter0.has_mcParticle() || !daughter1.has_mcParticle()) {
            LOGF(warning, "Encountered candidate daughter (%d or %d) without MC particle", daughter0.globalIndex(), daughter1.globalIndex());
            continue;
          }
          const auto mother0Id = daughter0.template mcParticle_as<McParticles>().template mothers_as<McParticles>().front().globalIndex();
          const auto mother1Id = daughter1.template mcParticle_as<McParticles>().template mothers_as<McParticles>().front().globalIndex();
          LOGF(debug, "MC candidate %d with prongs: %d (MC %d), %d (MC %d)", cand.globalIndex(),
               daughter0.globalIndex(), daughter0.template mcParticle_as<McParticles>().globalIndex(),
               daughter1.globalIndex(), daughter1.template mcParticle_as<McParticles>().globalIndex());
          LOGF(info, "MC ids of mothers: %d - %d", mother0Id, mother1Id);
          if ((mother0Id == mother1Id) &&
              std::abs(daughter0.template mcParticle_as<McParticles>().template mothers_as<McParticles>().front().flagMcMatchGen()) & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            LOGF(info, "D0 - looking for jet");
            for (const auto& pjet : jetsTagPerColl) {
              for (const auto& cand : pjet.template hfcandidates_as<McParticles>()) {
                if (mother0Id == cand.globalIndex()) {
                  matchedIdx = pjet.globalIndex();
                  LOGF(info, "Found match det to part: %d (pt %g) -> %d (pt %g)",
                       jet.globalIndex(), jet.pt(), matchedIdx, pjet.pt());
                }
              }
            }
          }
        }
      }
      jetsBaseToTag(matchedIdx, baseToTagIndexMap[jet.index()]); // TODO: check usage of index
    }

    // backward matching
    for (const auto& jet : jetsTagPerColl) {
      LOGF(info, "MC jet index: %d (coll %d) with %d tracks, %d HF candidates",
           jet.index(), jet.mcCollisionId(), jet.tracks().size(), jet.hfcandidates().size());

      int matchedIdx = -1;
      for (const auto& cand : jet.template hfcandidates_as<McParticles>()) {
        const auto& daughters = cand.template daughters_as<McParticles>();
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
               prong.globalIndex(), prong.prong0Id(), prong.template prong0_as<Tracks>().globalIndex(), prong.prong1Id(), prong.template prong1_as<Tracks>().globalIndex());
          if ((prong.template prong0_as<Tracks>().globalIndex() == index0 && prong.template prong1_as<Tracks>().globalIndex() == index1) ||
              (prong.prong0Id() == index1 && prong.prong1Id() == index0)) {
            candIdx = prong.globalIndex();
            LOGF(info, "Found matching 2prong candidate: %d", candIdx);
          }
        }
        for (const auto& djet : jetsBasePerColl) {
          if (djet.template hfcandidates_as<HfCandidates>().front().globalIndex() == candIdx) {
            matchedIdx = djet.globalIndex();
            LOGF(info, "Found match part to det: %d -> %d", jet.globalIndex(), matchedIdx);
          }
        }
      }
      jetsTagToBase(matchedIdx, tagToBaseIndexMap[jet.index()]); // TODO: check usage of index
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetMatchingHF<soa::Join<aod::D0MCDJets, aod::D0MCDJetConstituents>,
                                    soa::Join<aod::D0MCPJets, aod::D0MCPJetConstituents>,
                                    aod::D0MCDJetsMatchedToD0MCPJets, aod::D0MCPJetsMatchedToD0MCDJets,
                                    soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>>(cfgc, TaskName{"jet-matching-hf"})};
}
