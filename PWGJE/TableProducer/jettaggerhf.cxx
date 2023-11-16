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

// Task to produce a table joinable to the jet tables for hf jet tagging
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename JetTableMCD, typename JetTaggingTableData, typename JetTaggingTableMCD>
struct JetTaggerHFTask {

  Produces<JetTaggingTableData> taggingTableData;
  Produces<JetTaggingTableMCD> taggingTableMCD;

  Configurable<bool> doAlgorithm1{"doAlgorithm1", false, "fill table for algoithm 1"};
  Configurable<bool> doAlgorithm2{"doAlgorithm2", false, "fill table for algoithm 2"};
  Configurable<bool> doAlgorithm3{"doAlgorithm3", false, "fill table for algoithm 3"};
  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};

  void processDummy(aod::Collision const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDummy, "Dummy process", true);

  void processData(aod::JCollision const& collision, JetTableData const& jets, aod::JTracks const& tracks)
  {
    for (auto& jet : jets) {

      int algorithm1 = jet.globalIndex(); // This needs to be changed. It is only done because O2Physics compilation breaks if jet is unused
      int algorithm2 = 0;
      int algorithm3 = 0;
      // if (doAlgorithm1) algorithm1 = JetTaggingUtilities::Algorithm1((mcdjet, tracks);
      // if (doAlgorithm2) algorithm2 = JetTaggingUtilities::Algorithm2((mcdjet, tracks);
      // if (doAlgorithm3) algorithm3 = JetTaggingUtilities::Algorithm3((mcdjet, tracks);
      taggingTableData(0, algorithm1, algorithm2, algorithm3);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processData, "Fill tagging decision for data jets", false);

  void processMCD(aod::JCollision const& collision, JetTableMCD const& mcdjets, soa::Join<aod::JTracks, aod::McTrackLabels> const& tracks, aod::JMcParticles const& particles)
  {
    for (auto& mcdjet : mcdjets) {

      int origin = JetTaggingUtilities::mcdJetFromHFShower(mcdjet, tracks, particles, maxDeltaR);
      int algorithm1 = 0;
      int algorithm2 = 0;
      int algorithm3 = 0;
      // if (doAlgorithm1) algorithm1 = JetTaggingUtilities::Algorithm1((mcdjet, tracks);
      // if (doAlgorithm2) algorithm2 = JetTaggingUtilities::Algorithm2((mcdjet, tracks);
      // if (doAlgorithm3) algorithm3 = JetTaggingUtilities::Algorithm3((mcdjet, tracks);
      taggingTableMCD(origin, algorithm1, algorithm2, algorithm3);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processMCD, "Fill tagging decision for mcd jets", false);
};

using JetTaggerChargedJets = JetTaggerHFTask<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedJetTags, aod::ChargedMCDetectorLevelJetTags>;
using JetTaggerFullJets = JetTaggerHFTask<soa::Join<aod::FullJets, aod::FullJetConstituents>, soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>, aod::FullJetTags, aod::FullMCDetectorLevelJetTags>;
// using JetTaggerNeutralJets = JetTaggerHFTask<soa::Join<aod::NeutralJets, aod::NeutralJetConstituents>,soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents>, aod::NeutralJetTags, aod::NeutralMCDetectorLevelJetTags>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerChargedJets>(cfgc,
                                            SetDefaultProcesses{}, TaskName{"jet-taggerhf-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerFullJets>(cfgc,
                                         SetDefaultProcesses{}, TaskName{"jet-taggerhf-full"}));
  /*
    tasks.emplace_back(
      adaptAnalysisTask<JetTaggerNeutralJets>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-taggerhf-neutral"}));
  */
  return WorkflowSpec{tasks};
}
