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

// FullJet Spectra in pp
//
/// \author Archita Rani Dash <archita.rani.dash@cern.ch>
#include <vector>
#include <iostream>
#include <utility>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"
// #include "PWGJE/Core/FastJetUtilitites.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace std;
using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename JetConstituentTableData>
struct FullJetSpectraTask {

  HistogramRegistry registry;

  // Event configurables
  Configurable<float> VertexZCut{"VertexZCut", 10.0f, "Accepted z-vertex range"};

  // Jet configurables
  Configurable<double> JetRadii{"JetRadii", 0.4, "jet resolution parameters"};
  Configurable<float> jetpTMin{"jetpTMin", 0., "minimum jet pT"};
  Configurable<float> jetpTMax{"jetpTMax", 350., "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -1.0, "minimum jet eta"};
  Configurable<float> jetEtaMax{"jetEtaMax", 1.0, "maximum jet eta"};
  Configurable<float> jetPhiMin{"jetPhiMin", 0., "minimum jet phi"};
  Configurable<float> jetPhiMax{"jetPhiMax", 7., "maximum jet phi"};

  // Track configurables
  Configurable<float> trackpTMin{"trackpTMin", 0., "minimum track pT"};
  Configurable<float> trackpTMax{"trackpTMax", 200., "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -1.0, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 1.0, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", 0., "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 7., "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // Cluster configurables
  Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"};
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};
  Configurable<float> clusterPhiMin{"clusterPhiMin", -999., "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 999., "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -999., "minimum cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 999., "maximum cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  int trackSelection = -1;
  int eventSelection = -1;
  std::string particleSelection;

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    // JetTrack QA histograms
    if (doprocessTracks) {
      registry.add("h_collisions", "event status; event status;entries", {HistType::kTH1F, {{4, 0., 4.0}}});
      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, 0., 7.}}});

      // Cluster QA histograms
      registry.add("h_cluster_pt", "cluster pT;#it{p}_{T_cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_cluster_eta", "cluster #eta;#eta_{cluster};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_cluster_phi", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_cluster_energy", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}});
    }

    // Jet QA histograms
    if (doprocessJetsData) {
      registry.add("h_full_jet_pt", "jet pT;#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
    }
  } // init

  using JetTableDataJoined = soa::Join<JetTableData, JetConstituentTableData>;

  // Applying some cuts(filters) on collisions, tracks, clusters
  Filter eventCuts = (nabs(aod::jcollision::posZ) < VertexZCut);
  Filter trackCuts = (aod::jtrack::pt >= trackpTMin && aod::jtrack::pt < trackpTMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));

  template <typename T>
  void fillJetHistograms(T const& jet, float weight = 1.0)
  {
    if (jet.r() == round(JetRadii * 100.0)) {
      registry.fill(HIST("h_full_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_full_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_full_jet_phi"), jet.phi(), weight);
    }
  }

  template <typename T, typename U>
  void fillTrackHistograms(T const& tracks, U const& clusters, float weight = 1.0)
  {
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt(), weight);
      registry.fill(HIST("h_track_eta"), track.eta(), weight);
      registry.fill(HIST("h_track_phi"), track.phi(), weight);
    }
    for (auto const& cluster : clusters) {
      double clusterpt = cluster.energy() / std::cosh(cluster.eta());
      registry.fill(HIST("h_cluster_pt"), clusterpt, weight);
      registry.fill(HIST("h_cluster_eta"), cluster.eta(), weight);
      registry.fill(HIST("h_cluster_phi"), cluster.phi(), weight);
      registry.fill(HIST("h_cluster_energy"), cluster.energy(), weight);
    }
  }

  void processDummy(JetCollisions const& collisions)
  {
  }
  PROCESS_SWITCH(FullJetSpectraTask, processDummy, "dummy task", true);

  void processJetsData(soa::Filtered<JetCollisions>::iterator const& collision, JetTableDataJoined const& jets, JetTracks const& tracks, JetClusters const& clusters)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectraTask, processJetsData, "Full Jets Data", false);

  void processTracks(JetCollision const& collision, soa::Filtered<JetTracks> const& tracks, soa::Filtered<JetClusters> const& clusters)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::eventEMCAL(collision)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    fillTrackHistograms(tracks, clusters);
  }
  PROCESS_SWITCH(FullJetSpectraTask, processTracks, "QA for fulljet tracks", false);
}; // struct

using FullJetsSpectraTask = FullJetSpectraTask<aod::FullJets, aod::FullJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<FullJetsSpectraTask>(cfgc, TaskName{"full-jet-spectra-pp"})}; }
