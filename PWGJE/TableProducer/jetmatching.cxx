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
struct JetMatching {
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

  static constexpr int8_t getHfFlag()
  {
    if (std::is_same<BaseToTagMatchingTable, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>::value &&
        std::is_same<TagToBaseMatchingTable, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>::value)
      return 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;

    if (std::is_same<BaseToTagMatchingTable, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>::value &&
        std::is_same<TagToBaseMatchingTable, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>::value)
      return 1 << aod::hf_cand_3prong::DecayType::LcToPKPi;

    if (std::is_same<BaseToTagMatchingTable, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>::value &&
        std::is_same<TagToBaseMatchingTable, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>::value)
      return 1 << aod::hf_cand_bplus::DecayType::BplusToD0Pi;

    return -1;
  }

  void init(InitContext const&)
  {
  }

  // for now:
  // BaseJetCollection must contain detector level jets
  // TagJetCollection must contain particle level jets
  void processJets(Collisions::iterator const& collision, aod::McCollisions const& mcCollisions,
                   BaseJetCollection const& jetsBase, TagJetCollection const& jetsTag,
                   Tracks const& tracks, McParticles const& particlesMC,
                   HfCandidates const& hfcandidates)
  {
    // TODO: check whether we need to handle jets in MC collisions without a reconstructed collision

    // waiting for framework fix to make sliced collection of same type as original collection:
    // const auto jetsBasePerColl = jetsBaseIsMC ? jetsBase.sliceBy(baseJetsPerCollision, collision.mcCollisionId());
    const auto jetsBasePerColl = jetsBase;
    const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, jetsTagIsMC ? collision.mcCollisionId() : collision.globalIndex());

    // geometric matching
    LOGF(info, "performing geometric matching for collision %d (%d / %d jets)",
         collision.globalIndex(), jetsBasePerColl.size(), jetsTagPerColl.size());
    std::vector<double> jetsBasePhi;
    std::vector<double> jetsBaseEta;
    for (const auto& jet : jetsBasePerColl) {
      LOGF(info, "base jet index: %d-%d (coll %d, pt %g, phi %g) with %d tracks, %d HF candidates",
           jet.index(), jet.globalIndex(), jet.collisionId(), jet.pt(), jet.phi(), jet.tracks().size(), jet.hfcandidates().size());
      jetsBasePhi.emplace_back(jet.phi());
      jetsBaseEta.emplace_back(jet.eta());
    }
    std::vector<double> jetsTagPhi;
    std::vector<double> jetsTagEta;
    for (const auto& jet : jetsTagPerColl) {
      LOGF(info, "tag jet index: %d-%d (coll %d, pt %g, phi %g) with %d tracks, %d HF candidates",
           jet.index(), jet.globalIndex(), jet.mcCollisionId(), jet.pt(), jet.phi(), jet.tracks().size(), jet.hfcandidates().size());
      jetsTagPhi.emplace_back(jet.phi());
      jetsTagEta.emplace_back(jet.eta());
    }
    auto&& [baseToTagGeo, tagToBaseGeo] = JetUtilities::MatchJetsGeometrically(jetsBasePhi, jetsBaseEta, jetsTagPhi, jetsTagEta, maxMatchingDistance);
    LOGF(debug, "geometric matching: %d - %d jets", baseToTagGeo.size(), tagToBaseGeo.size());
    for (std::size_t i = 0; i < baseToTagGeo.size(); ++i) {
      LOGF(debug, "bjet %i -> %i", i, baseToTagGeo[i]);
    }
    for (std::size_t i = 0; i < tagToBaseGeo.size(); ++i) {
      LOGF(debug, "tjet %i -> %i", i, tagToBaseGeo[i]);
    }

    // HF matching
    std::vector<int> baseToTagHF(jetsBasePerColl.size(), -1);
    std::vector<int> tagToBaseHF(jetsTagPerColl.size(), -1);
    if constexpr (getHfFlag() > 0) {
      LOGF(info, "performing HF matching for collision %d", collision.globalIndex());
      for (const auto& bjet : jetsBasePerColl) {
        LOGF(info, "jet index: %d (coll %d, pt %g, phi %g) with %d tracks, %d HF candidates",
             bjet.index(), bjet.collisionId(), bjet.pt(), bjet.phi(), bjet.tracks().size(), bjet.hfcandidates().size());

        const auto hfcand = bjet.template hfcandidates_first_as<HfCandidates>();
        if (hfcand.flagMcMatchRec() & getHfFlag()) {
          const auto hfCandMC = hfcand.template prong0_as<Tracks>().template mcParticle_as<McParticles>();
          const auto hfCandMcId = hfCandMC.template mothers_first_as<McParticles>().globalIndex();
          for (const auto& tjet : jetsTagPerColl) {
            const auto cand = tjet.template hfcandidates_first_as<McParticles>();
            if (cand.globalIndex() == hfCandMcId) {
              LOGF(info, "Found HF match: %d (pt %g) <-> %d (pt %g)",
                   bjet.globalIndex(), bjet.pt(), tjet.globalIndex(), tjet.pt());
              baseToTagHF[bjet.index()] = tjet.globalIndex();
              tagToBaseHF[tjet.index()] = bjet.globalIndex();
            }
          }
        }
      }
    }

    for (const auto& jet : jetsBasePerColl) {
      int geojetid = baseToTagGeo[jet.index()];
      if (geojetid > -1 && geojetid < jetsTagPerColl.size())
        geojetid = jetsTagPerColl.iteratorAt(geojetid).globalIndex();
      else
        geojetid = -1;
      LOGF(info, "registering matches for base jet %d (%d): geo -> %d (%d), HF -> %d",
           jet.index(), jet.globalIndex(), geojetid, baseToTagGeo[jet.index()], baseToTagHF[jet.index()]);
      jetsBaseToTag(geojetid, baseToTagHF[jet.index()]);
    }

    for (const auto& jet : jetsTagPerColl) {
      int geojetid = tagToBaseGeo[jet.index()];
      if (geojetid > -1 && geojetid < jetsBasePerColl.size())
        geojetid = jetsBasePerColl.iteratorAt(geojetid).globalIndex();
      else
        geojetid = -1;
      LOGF(info, "registering matches for tag jet %d (%d): geo -> %d (%d), HF -> %d",
           jet.index(), jet.globalIndex(), geojetid, tagToBaseGeo[jet.index()], tagToBaseHF[jet.index()]);
      jetsTagToBase(geojetid, tagToBaseHF[jet.index()]);
    }
  }
  PROCESS_SWITCH(JetMatching, processJets, "Perform jet matching", true);
};

using ChargedJetMatching = JetMatching<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>,
                                       soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>,
                                       aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets,
                                       aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets,
                                       aod::Tracks>;
using D0ChargedJetMatching = JetMatching<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>,
                                         soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>,
                                         aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets,
                                         aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets,
                                         soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;
using LcChargedJetMatching = JetMatching<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>,
                                         soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>,
                                         aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets,
                                         aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets,
                                         soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>;
using BplusChargedJetMatching = JetMatching<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>,
                                            soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>,
                                            aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets,
                                            aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets,
                                            soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<ChargedJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-ch"}));
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-d0-ch"}));
  tasks.emplace_back(adaptAnalysisTask<LcChargedJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-lc-ch"}));
  tasks.emplace_back(adaptAnalysisTask<BplusChargedJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-bplus-ch"}));

  return WorkflowSpec{tasks};
}
