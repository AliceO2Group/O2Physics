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

/// \file jetmatching.cxx
/// \brief Unified implementation of jet matching based on different criteria
/// expanding on previously separate implementations of geometric matching
/// (by Raymond Ehlers) and heavy-flavour matching
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
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
  static constexpr bool jetsBaseIsMC = o2::soa::relatedByIndex<aod::McCollisions, BaseJetCollection>();
  static constexpr bool jetsTagIsMC = o2::soa::relatedByIndex<aod::McCollisions, TagJetCollection>();
  Preslice<BaseJetCollection> baseJetsPerCollision = jetsBaseIsMC ? aod::jet::mcCollisionId : aod::jet::collisionId;
  Preslice<TagJetCollection> tagJetsPerCollision = jetsTagIsMC ? aod::jet::mcCollisionId : aod::jet::collisionId;

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
    // const auto jetsBasePerColl = jetsBase.sliceBy(baseJetsPerCollision, jetsBaseIsMC ? collision.mcCollisionId() : collision.globalIndex());
    const auto jetsBasePerColl = jetsBase;
    const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, jetsTagIsMC ? collision.mcCollisionId() : collision.globalIndex());

    // geometric matching
    std::vector<double> jetsBasePhi(jetsBasePerColl.size());
    std::vector<double> jetsBaseEta(jetsBasePerColl.size());
    for (auto jet : jetsBasePerColl) {
      jetsBasePhi.emplace_back(jet.phi());
      jetsBaseEta.emplace_back(jet.eta());
    }
    std::vector<double> jetsTagPhi(jetsTagPerColl.size());
    std::vector<double> jetsTagEta(jetsTagPerColl.size());
    for (auto& jet : jetsTagPerColl) {
      jetsTagPhi.emplace_back(jet.phi());
      jetsTagEta.emplace_back(jet.eta());
    }
    auto&& [baseToTagGeo, tagToBaseGeo] = JetUtilities::MatchJetsGeometrically(jetsBasePhi, jetsBaseEta, jetsTagPhi, jetsTagEta, maxMatchingDistance);

    std::vector<std::size_t> baseToTagHF(jetsBase.size(), -1);
    std::vector<std::size_t> tagToBaseHF(jetsTag.size(), -1);
    // HF matching
    for (const auto& bjet : jetsBasePerColl) {
      LOGF(info, "jet index: %d (coll %d, pt %g, phi %g) with %d tracks, %d HF candidates",
           bjet.index(), bjet.collisionId(), bjet.pt(), bjet.phi(), bjet.tracks().size(), bjet.hfcandidates().size());

      const auto& hfcand = bjet.template hfcandidates_as<HfCandidates>().front();
      if (hfcand.flagMcMatchRec() & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        const auto hfCandMcId = hfcand.template prong0_as<Tracks>().template mcParticle_as<McParticles>().template mothers_as<McParticles>().front().globalIndex();
        for (const auto& tjet : jetsTagPerColl) {
          const auto& cand = tjet.template hfcandidates_as<McParticles>().front();
          if (hfCandMcId == cand.globalIndex()) {
            LOGF(info, "Found match: %d (pt %g) <-> %d (pt %g)",
                 bjet.globalIndex(), bjet.pt(), tjet.globalIndex(), tjet.pt());
            baseToTagHF[bjet.index()] = tjet.globalIndex();
            tagToBaseHF[tjet.index()] = bjet.globalIndex();
          }
        }
      }
    }

    for (const auto& jet : jetsBasePerColl) {
      jetsBaseToTag(baseToTagGeo[jet.index()], baseToTagHF[jet.index()]);
    }

    for (const auto& jet : jetsTagPerColl) {
      jetsTagToBase(tagToBaseGeo[jet.index()], tagToBaseHF[jet.index()]);
    }
  }
};

using JetMatchingD0 = JetMatchingHF<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>,
                                    soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>,
                                    aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets,
                                    aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets,
                                    soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetMatchingD0>(cfgc, TaskName{"jet-matching-hf"})};
}
