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
/// \brief matching duplicate jets
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

template <typename JetsBase, typename JetsTag, typename JetsBasetoTagMatchingTable, typename JetsTagtoBaseMatchingTable, typename Tracks, typename Candidates>
struct JetMatchingDuplicates {

  Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};

  Produces<JetsBasetoTagMatchingTable> jetsBasetoTagMatchingTable;
  Produces<JetsTagtoBaseMatchingTable> jetsTagtoBaseMatchingTable;

  // preslicing jet collections, only for Mc-based collection
  static constexpr bool jetsBaseIsMc = o2::soa::relatedByIndex<aod::JMcCollisions, JetsBase>();
  static constexpr bool jetsTagIsMc = o2::soa::relatedByIndex<aod::JMcCollisions, JetsTag>();

  Preslice<JetsBase> baseJetsPerCollision = jetsBaseIsMc ? aod::jet::mcCollisionId : aod::jet::collisionId;
  Preslice<JetsTag> tagJetsPerCollision = jetsTagIsMc ? aod::jet::mcCollisionId : aod::jet::collisionId;

  void init(InitContext const&)
  {
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetMatchingDuplicates, processDummy, "Dummy process", true);

  void processJets(aod::JetCollisions const& collisions,
                   JetsBase const& jetsBase, JetsTag const& jetsTag,
                   Tracks const& tracks, Candidates const& candidates)
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

    for (const auto& collision : collisions) {

      const auto jetsBasePerColl = jetsBase.sliceBy(baseJetsPerCollision, collision.globalIndex());
      const auto jetsTagPerColl = jetsTag.sliceBy(tagJetsPerCollision, collision.globalIndex());
      // initialise template parameters as false since even if they are Mc we are not matching between detector and particle level
      jetmatchingutilities::doAllMatching<false, false>(jetsBasePerColl, jetsTagPerColl, jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF, jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF, candidates, candidates, tracks, tracks, tracks, tracks, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
    }

    for (auto i = 0; i < jetsBase.size(); ++i) {
      jetsBasetoTagMatchingTable(jetsBasetoTagMatchingGeo[i], jetsBasetoTagMatchingPt[i], jetsBasetoTagMatchingHF[i]); // is (and needs to) be filled in order
    }
    for (auto i = 0; i < jetsTag.size(); i++) {
      jetsTagtoBaseMatchingTable(jetsTagtoBaseMatchingGeo[i], jetsTagtoBaseMatchingPt[i], jetsTagtoBaseMatchingHF[i]); // is (and needs to) be filled in order
    }
  }
  PROCESS_SWITCH(JetMatchingDuplicates, processJets, "Perform jet matching", false);
};

using Charged1JetDataMatching = JetMatchingDuplicates<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>,
                                                      soa::Join<aod::Charged1Jets, aod::Charged1JetConstituents>,
                                                      aod::ChargedJetsMatchedToCharged1Jets,
                                                      aod::Charged1JetsMatchedToChargedJets,
                                                      aod::JTracks,
                                                      aod::JDummys>;

using Charged1JetMCDMatching = JetMatchingDuplicates<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>,
                                                     soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>,
                                                     aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets,
                                                     aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets,
                                                     aod::JTracks,
                                                     aod::JDummys>;

using Charged1JetMCPMatching = JetMatchingDuplicates<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>,
                                                     soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>,
                                                     aod::ChargedMCParticleLevelJetsMatchedToCharged1MCParticleLevelJets,
                                                     aod::Charged1MCParticleLevelJetsMatchedToChargedMCParticleLevelJets,
                                                     aod::JMcParticles,
                                                     aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<Charged1JetDataMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-data-ch-1"}));
  tasks.emplace_back(adaptAnalysisTask<Charged1JetMCDMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-mcd-ch-1"}));
  tasks.emplace_back(adaptAnalysisTask<Charged1JetMCPMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-mcp-ch-1"}));

  return WorkflowSpec{tasks};
}
