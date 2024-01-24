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
/// \brief matching event-wise constituent subtracted data jets and unsubtracted data jets
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

template <typename JetsBase, typename JetsTag, typename JetsBasetoTagMatchingTable, typename JetsTagtoBaseMatchingTable, typename TracksTag, typename Candidates>
struct JetMatchingSub {

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

  void processDummy(JetCollisions const& collisions)
  {
  }
  PROCESS_SWITCH(JetMatchingSub, processDummy, "Dummy process", true);

  void processJets(JetCollisions const& collisions,
                   JetsBase const& jetsBase, JetsTag const& jetsTag,
                   JetTracks const& tracks, TracksTag const& tracksSub, Candidates const& candidates)
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
  PROCESS_SWITCH(JetMatchingSub, processJets, "Perform jet matching", false);
};

using ChargedJetMatching = JetMatchingSub<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>,
                                          soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>,
                                          aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets,
                                          aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets,
                                          aod::JTrackSubs,
                                          aod::JDummys>;
using D0ChargedJetMatching = JetMatchingSub<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>,
                                            soa::Join<aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents>,
                                            aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets,
                                            aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets,
                                            aod::JTrackD0Subs,
                                            CandidatesD0Data>;
/*using LcChargedJetMatching = JetMatchingSub<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents>,
                                            soa::Join<aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents>,
                                            aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets,
                                            aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets,
                                            aod::JTrackLcSubs,
                                            CandidatesLcData>;
using BplusChargedJetMatching = JetMatchingSub<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents>,
                                               soa::Join<aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents>,
                                               aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets,
                                               aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets,
                                               aod::JTrackBplusSubs,
                                               CandidatesBplusData>;*/

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<ChargedJetMatching>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-sub-ch"}));
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetMatching>(cfgc, TaskName{"jet-matching-sub-d0-ch"}));
  // tasks.emplace_back(adaptAnalysisTask<LcChargedJetMatching>(cfgc, TaskName{"jet-matching-sub-lc-ch"}));
  // tasks.emplace_back(adaptAnalysisTask<BplusChargedJetMatching>(cfgc, TaskName{"jet-matching-sub-bplus-ch"}));

  return WorkflowSpec{tasks};
}
