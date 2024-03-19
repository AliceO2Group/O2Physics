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

#include <TF1.h>
#include <TH1.h>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/trackUtilities.h"

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
  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};

  using JetTagTracksData = soa::Join<JetTracks, aod::JTrackPIs, aod::JTracksTag>;
  using JetTagTracksMCD = soa::Join<JetTracksMCD, aod::JTrackPIs, aod::JTracksTag>;
  using OriTracksData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection>;
  using OriTracksMCD = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::McTrackLabels>;

  void processDummy(JetCollisions const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDummy, "Dummy process", true);

  void processData(JetCollision const& collision, JetTableData const& jets, JetTagTracksData const& jtracks, OriTracksData const& tracks)
  {
    for (auto& jet : jets) {
      float jetProb = 0;
      int algorithm2 = 0;
      int algorithm3 = 0;
      std::unique_ptr<TF1> fSignImpXYSig(jettaggingutilities::getResolutionFunction(0));
      if (doTC) {
        jetProb = jettaggingutilities::getJetProbability(fSignImpXYSig, collision, jet, jtracks, tracks);
      }
      // if (doSV) algorithm2 = jettaggingutilities::Algorithm2((mcdjet, tracks);
      taggingTableData(0, jetProb, algorithm2, algorithm3);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processData, "Fill tagging decision for data jets", false);

  void processMCD(JetCollision const& collision, JetTableMCD const& mcdjets, JetTagTracksMCD const& jtracks, OriTracksMCD const& tracks, JetParticles const& particles)
  {
    for (auto& mcdjet : mcdjets) {
      typename JetTagTracksMCD::iterator hftrack;
      int origin = 0;
      if (!doWShower)
        origin = jettaggingutilities::mcdJetFromHFShower(mcdjet, jtracks, particles, maxDeltaR);
      else
        origin = jettaggingutilities::jetTrackFromHFShower(mcdjet, jtracks, particles, hftrack);
      float jetProb = 0;
      int algorithm2 = 0;
      int algorithm3 = 0;
      std::unique_ptr<TF1> fSignImpXYSig(jettaggingutilities::getResolutionFunction(origin));
      if (doTC) {
        jetProb = jettaggingutilities::getJetProbability(fSignImpXYSig, collision, mcdjet, jtracks, tracks);
      }
      // if (doSV) algorithm2 = jettaggingutilities::Algorithm2((mcdjet, tracks);
      taggingTableMCD(origin, jetProb, algorithm2, algorithm3);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processMCD, "Fill tagging decision for mcd jets", false);

  void processTraining(JetCollision const& collision, JetTableMCD const& mcdjets, JetTagTracksMCD const& tracks)
  {
    // To create table for ML
  }
  PROCESS_SWITCH(JetTaggerHFTask, processTraining, "Fill tagging decision for mcd jets", false);
};

struct JetTaggerHFExtTask {

  Produces<aod::JTracksTag> jTracksTagTable;

  void init(InitContext const&)
  {
  }

  void processTracks(soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov>::iterator const& track)
  {
    float dcaXYZ = 0;
    float sigmaDcaXYZ2 = 0;

    jettaggingutilities::calculateDcaXYZ(dcaXYZ, sigmaDcaXYZ2, track.dcaXY(), track.dcaZ(), track.cYY(), track.cZY(), track.cZZ(), track.sigmaDcaXY2(), track.sigmaDcaZ2());
    jTracksTagTable(dcaXYZ, sigmaDcaXYZ2);
  }
  PROCESS_SWITCH(JetTaggerHFExtTask, processTracks, "produces derived track table for tagging", true);
};

using JetTaggerChargedJets = JetTaggerHFTask<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedJetTags, aod::ChargedMCDetectorLevelJetTags>;
using JetTaggerFullJets = JetTaggerHFTask<soa::Join<aod::FullJets, aod::FullJetConstituents>, soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>, aod::FullJetTags, aod::FullMCDetectorLevelJetTags>;
// using JetTaggerNeutralJets = JetTaggerHFTask<soa::Join<aod::NeutralJets, aod::NeutralJetConstituents>,soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents>, aod::NeutralJetTags, aod::NeutralMCDetectorLevelJetTags>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerHFExtTask>(cfgc,
                                          SetDefaultProcesses{}, TaskName{"jet-taggerhf-extension"}));

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
