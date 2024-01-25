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
#include "Common/Core/trackUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename JetTableMCD, typename JetTaggingTableData, typename StoredJetTaggingTableData, typename JetTaggingTableMCD, typename StoredJetTaggingTableMCD>
struct JetTaggerHFTask {

  Produces<JetTaggingTableData> taggingTableData;
  Produces<StoredJetTaggingTableData> storedTaggingTableData;
  Produces<JetTaggingTableMCD> taggingTableMCD;
  Produces<StoredJetTaggingTableMCD> storedTaggingTableMCD;

  Configurable<bool> doWShower{"doWShower", false, "find jet origin included gluon spliting"}; // true:: remove gluon spliting
  Configurable<bool> doTC{"doTC", false, "fill table for track counting algorithm"};
  Configurable<bool> doSV{"doSV", false, "fill table for secondary vertex algorithm"};
  Configurable<bool> doML{"doML", false, "fill table for machine learning"};
  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};

  using OriTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov>;
  using OriTracksMC = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::McTrackLabels>;

  void processDummy(aod::Collision const& collision)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDummy, "Dummy process", true);

  void processData(soa::Join<aod::JCollisions, aod::JCollisionPIs>::iterator const& jcollision, aod::Collisions&, JetTableData const& jets, soa::Join<JetTracks, aod::JTrackPIs> const& jtracks, OriTracks&)
  {
    auto collision = jcollision.template collision_as<aod::Collisions>();
    for (auto& jet : jets) {
      int algorithm1 = jet.globalIndex(); // This needs to be changed. It is only done because O2Physics compilation breaks if jet is unused
      int algorithm2 = 0;
      int algorithm3 = 0;
      // if (doTC) algorithm1 = jettaggingutilities::Algorithm1((mcdjet, tracks);
      // if (doSV) algorithm2 = jettaggingutilities::Algorithm2((mcdjet, tracks);
      // if (doML) algorithm3 = jettaggingutilities::Algorithm3((mcdjet, tracks);
      taggingTableData(0, algorithm1, algorithm2, algorithm3);

      if (doML) {
        std::vector<int> vecGeoSign;
        std::vector<float> vecTrackPt, vecTrackEta, vecTrackPhi, vecSignedIP2D, vecSignedIP2Ds, vecSignedIP3D, vecSignedIP3Ds;
        for (auto& jtrack : jet.template tracks_as<soa::Join<JetTracks, aod::JTrackPIs>>()) {
          auto track = jtrack.template track_as<OriTracks>();
          int geoSign = jettaggingutilities::getGeoSign(collision, jet, track);
          float dcaXYZ = 0.;
          float sigmaDcaXYZ2 = 0.;
          jettaggingutilities::calculateDcaXYZ(dcaXYZ, sigmaDcaXYZ2, track.dcaXY(), track.dcaZ(), track.cYY(), track.cZY(), track.cZZ(), track.sigmaDcaXY2(), track.sigmaDcaZ2());
          vecGeoSign.push_back(geoSign);
          vecTrackPt.push_back(track.pt());
          vecTrackEta.push_back(track.eta());
          vecTrackPhi.push_back(track.phi());
          vecSignedIP2D.push_back(track.dcaXY() / TMath::Sqrt(track.sigmaDcaXY2()) * jettaggingutilities::cmTomum);
          vecSignedIP2Ds.push_back(geoSign * TMath::Abs(track.dcaXY()) / TMath::Sqrt(track.sigmaDcaXY2()) * jettaggingutilities::cmTomum);
          vecSignedIP3D.push_back(dcaXYZ / TMath::Sqrt(sigmaDcaXYZ2) * jettaggingutilities::cmTomum);
          vecSignedIP3Ds.push_back(geoSign * TMath::Abs(dcaXYZ) / TMath::Sqrt(sigmaDcaXYZ2) * jettaggingutilities::cmTomum);
        }
        storedTaggingTableData(0, jet.pt(), jet.eta(), jet.phi(), vecGeoSign, vecTrackPt, vecTrackEta, vecTrackPhi, vecSignedIP2D, vecSignedIP2Ds, vecSignedIP3D, vecSignedIP3Ds);
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processData, "Fill tagging decision for data jets", false);

  void processMCD(soa::Join<aod::JCollisions, aod::JCollisionPIs>::iterator const& jcollision, aod::Collisions&, JetTableMCD const& mcdjets, JetTracksMCD const& jtracks, OriTracksMC&, JetParticles const& jparticles)
  {
    auto collision = jcollision.template collision_as<aod::Collisions>();

    for (auto& mcdjet : mcdjets) {
      typename JetTracksMCD::iterator hfjtrack;
      int origin = 0;
      if (!doWShower)
        origin = jettaggingutilities::mcdJetFromHFShower(mcdjet, jtracks, jparticles, maxDeltaR);
      else
        origin = jettaggingutilities::jetTrackFromHFShower(mcdjet, jtracks, jparticles, hfjtrack);
      int algorithm1 = 0;
      int algorithm2 = 0;
      int algorithm3 = 0;
      // if (doTC) algorithm1 = jettaggingutilities::Algorithm1((mcdjet, jtracks);
      // if (doSV) algorithm2 = jettaggingutilities::Algorithm2((mcdjet, jtracks);
      // if (doML) algorithm3 = jettaggingutilities::Algorithm3((mcdjet, jtracks);
      taggingTableMCD(origin, algorithm1, algorithm2, algorithm3);

      if (doML) {
        std::vector<int> vecGeoSign;
        std::vector<float> vecTrackPt, vecTrackEta, vecTrackPhi, vecSignedIP2D, vecSignedIP2Ds, vecSignedIP3D, vecSignedIP3Ds;
        for (auto& jtrack : mcdjet.template tracks_as<soa::Join<JetTracksMCD, aod::JTrackPIs>>()) {
          auto track = jtrack.template track_as<OriTracks>();
          int geoSign = jettaggingutilities::getGeoSign(collision, mcdjet, track);
          float dcaXYZ = 0.;
          float sigmaDcaXYZ2 = 0.;
          jettaggingutilities::calculateDcaXYZ(dcaXYZ, sigmaDcaXYZ2, track.dcaXY(), track.dcaZ(), track.cYY(), track.cZY(), track.cZZ(), track.sigmaDcaXY2(), track.sigmaDcaZ2());
          vecTrackPt.push_back(track.pt());
          vecTrackEta.push_back(track.eta());
          vecTrackPhi.push_back(track.phi());
          vecSignedIP2D.push_back(track.dcaXY() / TMath::Sqrt(track.sigmaDcaXY2()) * jettaggingutilities::cmTomum);
          vecSignedIP2Ds.push_back(geoSign * TMath::Abs(track.dcaXY()) / TMath::Sqrt(track.sigmaDcaXY2()) * jettaggingutilities::cmTomum);
          vecSignedIP3D.push_back(dcaXYZ / TMath::Sqrt(sigmaDcaXYZ2) * jettaggingutilities::cmTomum);
          vecSignedIP3Ds.push_back(geoSign * TMath::Abs(dcaXYZ) / TMath::Sqrt(sigmaDcaXYZ2) * jettaggingutilities::cmTomum);
        }
        storedTaggingTableMCD(origin, mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), vecGeoSign, vecTrackPt, vecTrackEta, vecTrackPhi, vecSignedIP2D, vecSignedIP2Ds, vecSignedIP3D, vecSignedIP3Ds);
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processMCD, "Fill tagging decision for mcd jets", false);
};

using JetTaggerChargedJets = JetTaggerHFTask<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedJetTags, aod::StoredChargedJetTags, aod::ChargedMCDetectorLevelJetTags, aod::StoredChargedMCDetectorLevelJetTags>;
using JetTaggerFullJets = JetTaggerHFTask<soa::Join<aod::FullJets, aod::FullJetConstituents>, soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>, aod::FullJetTags, aod::StoredFullJetTags, aod::FullMCDetectorLevelJetTags, aod::StoredFullMCDetectorLevelJetTags>;
// using JetTaggerNeutralJets = JetTaggerHFTask<soa::Join<aod::NeutralJets, aod::NeutralJetConstituents>, soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents>, aod::NeutralJetTags, aod::StoredNeutralJetTags, aod::NeutralMCDetectorLevelJetTags, aod::StoredNeutralMCDetectorLevelJetTags>;

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
