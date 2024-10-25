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

/// \file jetmatchingmc.cxx
/// \brief matching detector level and generator level jets
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Aimeric Lanodu <aimeric.landou@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetMatchingUtilities.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetsBase, typename JetsTag, typename JetsBasetoTagMatchingTable, typename JetsTagtoBaseMatchingTable, typename CandidatesBase, typename CandidatesTag, typename ClustersBase>
struct JetMatchingMc {

  Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};

  Produces<JetsBasetoTagMatchingTable> jetsBasetoTagMatchingTable;
  Produces<JetsTagtoBaseMatchingTable> jetsTagtoBaseMatchingTable;

  // preslicing jet collections, only for Mc-based collection
  static constexpr bool jetsBaseIsMc = o2::soa::relatedByIndex<aod::JetMcCollisions, JetsBase>();
  static constexpr bool jetsTagIsMc = o2::soa::relatedByIndex<aod::JetMcCollisions, JetsTag>();

  Preslice<JetsBase> baseJetsPerCollision = jetsBaseIsMc ? aod::jet::mcCollisionId : aod::jet::collisionId;
  Preslice<JetsTag> tagJetsPerCollision = jetsTagIsMc ? aod::jet::mcCollisionId : aod::jet::collisionId;

  PresliceUnsorted<aod::JetCollisionsMCD> CollisionsPerMcCollision = aod::jmccollisionlb::mcCollisionId;

  void init(InitContext const&)
  {
  }

  void processDummy(aod::JetMcCollisions const&)
  {
  }
  PROCESS_SWITCH(JetMatchingMc, processDummy, "Dummy process", true);

  void processJets(aod::JetMcCollisions const& mcCollisions, aod::JetCollisionsMCD const& collisions,
                   JetsBase const& jetsBase, JetsTag const& jetsTag,
                   aod::JetTracksMCD const& tracks,
                   ClustersBase const& clusters,
                   aod::JetParticles const& particles,
                   CandidatesBase const& candidatesBase,
                   CandidatesTag const& candidatesTag)
  {

    // initialise objects used to store the matching index arrays (array in case a mcCollision is split) before filling the matching tables
    std::vector<std::vector<int>> jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF;
    std::vector<std::vector<int>> jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF;
    //  waiting for framework fix to make sliced collection of same type as original collection:
    jetsBasetoTagMatchingGeo.assign(jetsBase.size(), {});
    jetsBasetoTagMatchingPt.assign(jetsBase.size(), {});
    jetsBasetoTagMatchingHF.assign(jetsBase.size(), {});
    jetsTagtoBaseMatchingGeo.assign(jetsTag.size(), {});
    jetsTagtoBaseMatchingPt.assign(jetsTag.size(), {});
    jetsTagtoBaseMatchingHF.assign(jetsTag.size(), {});

    for (const auto& mcCollision : mcCollisions) {

      const auto collisionsPerMcColl = collisions.sliceBy(CollisionsPerMcCollision, mcCollision.globalIndex());

      for (const auto& collision : collisionsPerMcColl) {

        const auto jetsBasePerColl = jetsBase.sliceBy(baseJetsPerCollision, jetsBaseIsMc ? mcCollision.globalIndex() : collision.globalIndex());
        const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, jetsTagIsMc ? mcCollision.globalIndex() : collision.globalIndex());

        jetmatchingutilities::doAllMatching<jetsBaseIsMc, jetsTagIsMc>(jetsBasePerColl, jetsTagPerColl, jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF, jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF, candidatesBase, candidatesTag, tracks, clusters, particles, particles, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
      }
    }
    for (auto i = 0; i < jetsBase.size(); ++i) {
      jetsBasetoTagMatchingTable(jetsBasetoTagMatchingGeo[i], jetsBasetoTagMatchingPt[i], jetsBasetoTagMatchingHF[i]); // is (and needs to) be filled in order
    }
    for (auto i = 0; i < jetsTag.size(); i++) {
      jetsTagtoBaseMatchingTable(jetsTagtoBaseMatchingGeo[i], jetsTagtoBaseMatchingPt[i], jetsTagtoBaseMatchingHF[i]); // is (and needs to) be filled in order
    }
  }
  PROCESS_SWITCH(JetMatchingMc, processJets, "Perform jet matching", false);
};

using ChargedJetMatching = JetMatchingMc<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>,
                                         soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>,
                                         aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets,
                                         aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets,
                                         aod::JCollisions,
                                         aod::JMcCollisions,
                                         aod::JDummys>;
using FullJetMatching = JetMatchingMc<soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>,
                                      soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>,
                                      aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets,
                                      aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets,
                                      aod::JCollisions,
                                      aod::JMcCollisions,
                                      aod::JetClustersMCD>;
using NeutralJetMatching = JetMatchingMc<soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents>,
                                         soa::Join<aod::NeutralMCParticleLevelJets, aod::NeutralMCParticleLevelJetConstituents>,
                                         aod::NeutralMCDetectorLevelJetsMatchedToNeutralMCParticleLevelJets,
                                         aod::NeutralMCParticleLevelJetsMatchedToNeutralMCDetectorLevelJets,
                                         aod::JCollisions,
                                         aod::JMcCollisions,
                                         aod::JetClustersMCD>;
using D0ChargedJetMatching = JetMatchingMc<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>,
                                           soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>,
                                           aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets,
                                           aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets,
                                           aod::CandidatesD0MCD,
                                           aod::CandidatesD0MCP,
                                           aod::JDummys>;
using LcChargedJetMatching = JetMatchingMc<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>,
                                           soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>,
                                           aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets,
                                           aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets,
                                           aod::CandidatesLcMCD,
                                           aod::CandidatesLcMCP,
                                           aod::JDummys>;
/*using BplusChargedJetMatching = JetMatchingMc<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>,
                                              soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>,
                                              aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets,
                                              aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets,
                                              aod::CandidatesBplusMCD,
                                              aod::CandidatesBplusMCP,
                                              aod::JDummys>>;*/
using V0ChargedJetMatching = JetMatchingMc<soa::Join<aod::V0ChargedMCDetectorLevelJets, aod::V0ChargedMCDetectorLevelJetConstituents>,
                                           soa::Join<aod::V0ChargedMCParticleLevelJets, aod::V0ChargedMCParticleLevelJetConstituents>,
                                           aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets,
                                           aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets,
                                           aod::CandidatesV0MCD,
                                           aod::CandidatesV0MCP,
                                           aod::JDummys>;

using DielectronChargedJetMatching = JetMatchingMc<soa::Join<aod::DielectronChargedMCDetectorLevelJets, aod::DielectronChargedMCDetectorLevelJetConstituents>,
                                                   soa::Join<aod::DielectronChargedMCParticleLevelJets, aod::DielectronChargedMCParticleLevelJetConstituents>,
                                                   aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets,
                                                   aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets,
                                                   aod::CandidatesDielectronMCD,
                                                   aod::CandidatesDielectronMCP,
                                                   aod::JDummys>;
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<ChargedJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-mc-ch"}));
  tasks.emplace_back(adaptAnalysisTask<FullJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-mc-full"}));
  tasks.emplace_back(adaptAnalysisTask<NeutralJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-mc-neutral"}));
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetMatching>(cfgc, TaskName{"jet-matching-mc-d0-ch"}));
  tasks.emplace_back(adaptAnalysisTask<LcChargedJetMatching>(cfgc, TaskName{"jet-matching-mc-lc-ch"}));
  // tasks.emplace_back(adaptAnalysisTask<BplusChargedJetMatching>(cfgc, TaskName{"jet-matching-mc-bplus-ch"}));
  tasks.emplace_back(adaptAnalysisTask<V0ChargedJetMatching>(cfgc, TaskName{"jet-matching-mc-v0-ch"}));
  tasks.emplace_back(adaptAnalysisTask<DielectronChargedJetMatching>(cfgc, TaskName{"jet-matching-mc-dielectron-ch"}));

  return WorkflowSpec{tasks};
}
