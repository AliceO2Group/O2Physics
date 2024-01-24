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
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename JetTableMCD, typename JetTaggingTableData, typename ConstTaggingTableData, typename ConstTCTaggingTableData, typename JetTaggingTableMCD, typename ConstTaggingTableMCD, typename ConstTCTaggingTableMCD>
struct JetTaggerHFTask {

  Produces<JetTaggingTableData> taggingTableData;
  Produces<ConstTaggingTableData> constTableData;
  Produces<ConstTCTaggingTableData> constTCTableData;
  Produces<JetTaggingTableMCD> taggingTableMCD;
  Produces<ConstTaggingTableMCD> constTableMCD;
  Produces<ConstTCTaggingTableMCD> constTCTableMCD;

  Configurable<bool> doWShower{"doWShower", false, "find jet origin included gluon spliting"}; // true:: remove gluon spliting
  Configurable<bool> doAlgorithm1{"doAlgorithm1", false, "fill table for algoithm 1"};
  Configurable<bool> doAlgorithm2{"doAlgorithm2", false, "fill table for algoithm 2"};
  Configurable<bool> doAlgorithm3{"doAlgorithm3", false, "fill table for algoithm 3"};
  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};

  using JetTagTracksData = soa::Join<aod::JTracks, aod::JTrackTagDcas, aod::JTrackTagDcaCovs>;
  using JetTagTracksMC = soa::Join<aod::JTracks, aod::JMcTrackLbs, aod::JTrackTagDcas, aod::JTrackTagDcaCovs>;


  void processDummy(aod::Collision const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDummy, "Dummy process", true);

  void processData(soa::Join<aod::JCollisions, aod::JCollisionPIs>::iterator const& jcollision, aod::Collisions&, JetTableData const& jets, aod::JTracks&, JetTagTracksData const& tracks)
  {
    auto collision = jcollision.template collision_as<aod::Collisions>();
    for (auto& jet : jets) {
      int algorithm1 = jet.globalIndex(); // This needs to be changed. It is only done because O2Physics compilation breaks if jet is unused
      int algorithm2 = 0;
      int algorithm3 = 0;
      // if (doAlgorithm1) algorithm1 = JetTaggingUtilities::Algorithm1((mcdjet, tracks);
      // if (doAlgorithm2) algorithm2 = JetTaggingUtilities::Algorithm2((mcdjet, tracks);
      // if (doAlgorithm3) algorithm3 = JetTaggingUtilities::Algorithm3((mcdjet, tracks);
      taggingTableData(0, algorithm1, algorithm2, algorithm3);
      std::vector<int> trackconst;
      std::vector<int> vecGeoSign;
      std::vector<float> vecTrackPt, vecTrackEta, vecTrackPhi, vecSignedIP2D, vecSignedIP2Ds, vecSignedIP3D, vecSignedIP3Ds;
      for (auto& track : jet.template tracks_as<JetTagTracksData>()) {
        trackconst.push_back(track.globalIndex());
        int geoSign = JetTaggingUtilities::getGeoSign(collision, jet, track);
        trackconst.push_back(track.globalIndex());
        vecGeoSign.push_back(geoSign);
        vecTrackPt.push_back(track.pt());
        vecTrackEta.push_back(track.eta());
        vecTrackPhi.push_back(track.phi());
        vecSignedIP2D.push_back(track.dcaXY() / TMath::Sqrt(track.sigmaDcaXY2()) * JetTaggingUtilities::cmTomum);
        vecSignedIP2Ds.push_back(geoSign * TMath::Abs(track.dcaXY()) / TMath::Sqrt(track.sigmaDcaXY2()) * JetTaggingUtilities::cmTomum);
        vecSignedIP3D.push_back(track.dcaXYZ() / TMath::Sqrt(track.sigmaDcaXYZ2()) * JetTaggingUtilities::cmTomum);
        vecSignedIP3Ds.push_back(geoSign * TMath::Abs(track.dcaXYZ()) / TMath::Sqrt(track.sigmaDcaXYZ2()) * JetTaggingUtilities::cmTomum);
      }
      constTableData(taggingTableData.lastIndex(), trackconst);
      constTCTableData(taggingTableData.lastIndex(), 0, jet.pt(), jet.eta(), jet.phi(), vecGeoSign, vecTrackPt, vecTrackEta, vecTrackPhi, vecSignedIP2D, vecSignedIP2Ds, vecSignedIP3D, vecSignedIP3Ds);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processData, "Fill tagging decision for data jets", false);

  void processMCD(soa::Join<aod::JCollisions, aod::JCollisionPIs>::iterator const& jcollision, aod::Collisions&, JetTableMCD const& mcdjets, aod::JTracks&, JetTagTracksMC const& tracks, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles)
  {
    auto collision = jcollision.template collision_as<aod::Collisions>();

    for (auto& mcdjet : mcdjets) {
      typename JetTagTracksMC::iterator hftrack;
      int origin = 0;
      if (!doWShower)
        origin = JetTaggingUtilities::mcdJetFromHFShower(mcdjet, tracks, particles, maxDeltaR);
      else
        origin = JetTaggingUtilities::jetTrackFromHFShower(mcdjet, tracks, particles, hftrack);
      int algorithm1 = 0;
      int algorithm2 = 0;
      int algorithm3 = 0;
      // if (doAlgorithm1) algorithm1 = JetTaggingUtilities::Algorithm1((mcdjet, tracks);
      // if (doAlgorithm2) algorithm2 = JetTaggingUtilities::Algorithm2((mcdjet, tracks);
      // if (doAlgorithm3) algorithm3 = JetTaggingUtilities::Algorithm3((mcdjet, tracks);
      taggingTableMCD(origin, algorithm1, algorithm2, algorithm3);
      std::vector<int> trackconst;
      std::vector<int> vecGeoSign;
      std::vector<float> vecTrackPt, vecTrackEta, vecTrackPhi, vecSignedIP2D, vecSignedIP2Ds, vecSignedIP3D, vecSignedIP3Ds;
      for (auto& track : mcdjet.template tracks_as<JetTagTracksMC>()) {
        int geoSign = JetTaggingUtilities::getGeoSign(collision, mcdjet, track);
        trackconst.push_back(track.globalIndex());
        vecGeoSign.push_back(geoSign);
        vecTrackPt.push_back(track.pt());
        vecTrackEta.push_back(track.eta());
        vecTrackPhi.push_back(track.phi());
        vecSignedIP2D.push_back(track.dcaXY() / TMath::Sqrt(track.sigmaDcaXY2()));
        vecSignedIP2Ds.push_back(geoSign * TMath::Abs(track.dcaXY()) / TMath::Sqrt(track.sigmaDcaXY2()));
        vecSignedIP3D.push_back(track.dcaXYZ() / TMath::Sqrt(track.sigmaDcaXYZ2()));
        vecSignedIP3Ds.push_back(geoSign * TMath::Abs(track.dcaXYZ()) / TMath::Sqrt(track.sigmaDcaXYZ2()));
      }
      constTableMCD(taggingTableMCD.lastIndex(), trackconst);
      constTCTableMCD(taggingTableMCD.lastIndex(), origin, mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), vecGeoSign, vecTrackPt, vecTrackEta, vecTrackPhi, vecSignedIP2D, vecSignedIP2Ds, vecSignedIP3D, vecSignedIP3Ds);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processMCD, "Fill tagging decision for mcd jets", false);
};

using JetTaggerChargedJets = JetTaggerHFTask<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedJetTags, aod::ChargedJetTagConstituents, aod::ChargedJetTagConstituentTCs, aod::ChargedMCDetectorLevelJetTags, aod::ChargedMCDetectorLevelJetTagConstituents, aod::ChargedMCDetectorLevelJetTagConstituentTCs>;
using JetTaggerFullJets = JetTaggerHFTask<soa::Join<aod::FullJets, aod::FullJetConstituents>, soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>, aod::FullJetTags, aod::FullJetTagConstituents, aod::FullJetTagConstituentTCs, aod::FullMCDetectorLevelJetTags, aod::FullMCDetectorLevelJetTagConstituents, aod::FullMCDetectorLevelJetTagConstituentTCs>;
// using JetTaggerNeutralJets = JetTaggerHFTask<soa::Join<aod::NeutralJets, aod::NeutralJetConstituents>,soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents>, aod::NeutralJetTags, aod::NeutralJetTagConstituents, aod::NeutralMCDetectorLevelJetTags, aod::NeutralMCDetectorLevelJetTagConstituents>;

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
