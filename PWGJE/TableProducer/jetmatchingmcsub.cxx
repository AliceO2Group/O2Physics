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

/// \file jetmatchingmcsub.cxx
/// \brief matching event-wise constituent subtracted detector level and unsubtracted generated level jets (this is usseful as a template for embedding  matching)
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

template <typename JetsBase, typename JetsTag, typename JetsBasetoTagMatchingTable, typename JetsTagtoBaseMatchingTable, typename Candidates>
struct JetMatchingMcSub {

  Configurable<bool> doMatchingGeo{"doMatchingGeo", true, "Enable geometric matching"};
  Configurable<bool> doMatchingPt{"doMatchingPt", true, "Enable pt matching"};
  Configurable<bool> doMatchingHf{"doMatchingHf", false, "Enable HF matching"};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.24f, "Max matching distance"};
  Configurable<float> minPtFraction{"minPtFraction", 0.5f, "Minimum pt fraction for pt matching"};

  Produces<JetsBasetoTagMatchingTable> jetsBasetoTagMatchingTable;
  Produces<JetsTagtoBaseMatchingTable> jetsTagtoBaseMatchingTable;

  // preslicing jet collections, only for Mc-based collection
  static constexpr bool jetsBaseIsMc = false;
  static constexpr bool jetsTagIsMc = false;

  Preslice<JetsBase> baseJetsPerCollision = aod::jet::collisionId;
  Preslice<JetsTag> tagJetsPerCollision = aod::jet::collisionId;

  void init(InitContext const&)
  {
  }

  void processDummy(JetCollisions const& mcCollisions)
  {
  }
  PROCESS_SWITCH(JetMatchingMcSub, processDummy, "Dummy process", true);

  void processJets(JetCollisions const& collisions,
                   JetsBase const& jetsBase, JetsTag const& jetsTag,
                   JetTracks const& tracks,
                   JetTracksSub const& tracksSub,
                   Candidates const& candidates)
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

      jetmatchingutilities::doAllMatching<jetsBaseIsMc, jetsTagIsMc>(jetsBasePerColl, jetsTagPerColl, jetsBasetoTagMatchingGeo, jetsBasetoTagMatchingPt, jetsBasetoTagMatchingHF, jetsTagtoBaseMatchingGeo, jetsTagtoBaseMatchingPt, jetsTagtoBaseMatchingHF, candidates, candidates, tracks, tracksSub, doMatchingGeo, doMatchingHf, doMatchingPt, maxMatchingDistance, minPtFraction);
    }

    for (auto i = 0; i < jetsBase.size(); ++i) {
      jetsBasetoTagMatchingTable(jetsBasetoTagMatchingGeo[i], jetsBasetoTagMatchingPt[i], jetsBasetoTagMatchingHF[i]); // is (and needs to) be filled in order
    }
    for (auto i = 0; i < jetsTag.size(); i++) {
      jetsTagtoBaseMatchingTable(jetsTagtoBaseMatchingGeo[i], jetsTagtoBaseMatchingPt[i], jetsTagtoBaseMatchingHF[i]); // is (and needs to) be filled in order
    }
  }
  PROCESS_SWITCH(JetMatchingMcSub, processJets, "Perform jet matching", false);
};

using ChargedJetMatching = JetMatchingMcSub<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>,
                                            soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                            aod::ChargedMCDetectorLevelJetsMatchedToChargedMCDetectorLevelEventWiseSubtractedJets,
                                            aod::ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToChargedMCDetectorLevelJets,
                                            aod::JDummys>;
using D0ChargedJetMatching = JetMatchingMcSub<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>,
                                              soa::Join<aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                              aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCDetectorLevelEventWiseSubtractedJets,
                                              aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToD0ChargedMCDetectorLevelJets,
                                              CandidatesD0MCD>;
/*using LcChargedJetMatching = JetMatchingMcSub<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>,
                                              soa::Join<aod::LcChargedMCDetectorLevelEventWiseSubtractedJets, aod::LcChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                              aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCDetectorLevelEventWiseSubtractedJets,
                                              aod::LcChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToLcChargedMCDetectorLevelJets,
                                              CandidatesLcMCD>;
using BplusChargedJetMatching = JetMatchingMcSub<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>,
                                                 soa::Join<aod::BplusChargedMCDetectorLevelEventWiseSubtractedJets, aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                                 aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCDetectorLevelEventWiseSubtractedJets,
                                                 aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToBplusChargedMCDetectorLevelJets,
                                                 CandidatesBplusMCD>;*/

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<ChargedJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-mc-sub-ch"}));
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetMatching>(cfgc, TaskName{"jet-matching-mc-sub-d0-ch"}));
  // tasks.emplace_back(adaptAnalysisTask<LcChargedJetMatching>(cfgc, TaskName{"jet-matching-mc-sub-lc-ch"}));
  // tasks.emplace_back(adaptAnalysisTask<BplusChargedJetMatching>(cfgc, TaskName{"jet-matching-mc-sub-bplus-ch"}));

  return WorkflowSpec{tasks};
}
