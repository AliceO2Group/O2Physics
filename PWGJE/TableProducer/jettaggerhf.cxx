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
/// \author Hanseo Park <hanseo.park@cern.ch>

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

  Configurable<bool> doWShower{"doWShower", false, "find jet origin included gluon spliting"}; // true:: remove gluon spliting
  Configurable<bool> doTC{"doTC", false, "fill table for track counting algorithm"};
  Configurable<bool> doSV{"doSV", false, "fill table for secondary vertex algorithm"};
  Configurable<bool> doML{"doML", false, "fill table for machine learning"};
  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};

  void processDummy(JetCollisions const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDummy, "Dummy process", true);

  void processData(JetCollision const& collision, JetTableData const& jets, aod::JTracks const& tracks)
  {
    for (auto& jet : jets) {

      int algorithm1 = jet.globalIndex(); // This needs to be changed. It is only done because O2Physics compilation breaks if jet is unused
      int algorithm2 = 0;
      int algorithm3 = 0;
      // if (doTC) algorithm1 = jettaggingutilities::Algorithm1((mcdjet, tracks);
      // if (doSV) algorithm2 = jettaggingutilities::Algorithm2((mcdjet, tracks);
      // if (doML) algorithm3 = jettaggingutilities::Algorithm3((mcdjet, tracks);
      taggingTableData(0, algorithm1, algorithm2, algorithm3);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processData, "Fill tagging decision for data jets", false);

  void processMCD(JetCollision const& collision, JetTableMCD const& mcdjets, JetTracksMCD const& tracks, JetParticles const& particles)
  {
    for (auto& mcdjet : mcdjets) {
      typename JetTracksMCD::iterator hftrack;
      int origin = 0;
      if (!doWShower)
        origin = jettaggingutilities::mcdJetFromHFShower(mcdjet, tracks, particles, maxDeltaR); // TODO: it is not working due to getOriginalMotherIndex
      else
        origin = jettaggingutilities::jetTrackFromHFShower(mcdjet, tracks, particles, hftrack);
      int algorithm1 = 0;
      int algorithm2 = 0;
      int algorithm3 = 0;
      // if (doTC) algorithm1 = jettaggingutilities::Algorithm1((mcdjet, tracks);
      // if (doSV) algorithm2 = jettaggingutilities::Algorithm2((mcdjet, tracks);
      // if (doML) algorithm3 = jettaggingutilities::Algorithm3((mcdjet, tracks);
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
