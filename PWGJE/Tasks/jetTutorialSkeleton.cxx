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

// jet tutorial skeleton task for hands on tutorial session (09/11/2023)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetTutorialSkeletonTask {
  HistogramRegistry registry{"registry",
                             {{"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}}},
                              {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_track_chi2PerCluster", "track #chi^{2} per cluste ;#chi^{2};entries", {HistType::kTH1F, {{100, 0, 40}}}},
                              {"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_jet_pt_bkgsub", "jet pT bkg sub;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_part_jet_pt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_part_jet_eta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_part_jet_phi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_full_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_full_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_full_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_full_jet_nclusters", "jet N clusters;N_{jet clusters};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_full_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_part_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_recoil_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_recoil_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_recoil_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_recoil_jet_dphi", "hadron-jet #Delta#phi;#Delta#phi_{jet,trigger hadron};entries", {HistType::kTH1F, {{40, -2.0, 2.0}}}},
                              {"h_matched_jets_pt", "#it{p}_{T,jet part}; #it{p}_{T,jet det}", {HistType::kTH2F, {{100, 0., 20.}, {100, 0., 20.0}}}},
                              {"h_matched_jets_eta", "#eta_{jet part}; #eta_{jet det}", {HistType::kTH2F, {{100, -1.0, 1.0}, {100, -1.0, 1.0}}}},
                              {"h_matched_jets_phi", "#phi_{jet part}; #phi_{jet det}", {HistType::kTH2F, {{80, -1.0, 7.}, {80, -1.0, 7.}}}}}};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = JetDerivedDataUtilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = JetDerivedDataUtilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);

  void processCollisions(aod::JCollision const& collision, aod::JTracks const& tracks)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processCollisions, "process self contained collisions", true);

  void processCollisionsWithExternalTracks(aod::JCollision const& collision, soa::Join<aod::JTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection> const& originalTracks)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processCollisionsWithExternalTracks, "process non self contained collisions", true);

  void processDataCharged(soa::Filtered<aod::ChargedJets>::iterator const& jet)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataCharged, "jets data", true);

  void processMCDetectorLevelCharged(soa::Filtered<aod::ChargedMCDetectorLevelJets>::iterator const& jet)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCDetectorLevelCharged, "jets on detector level MC", false);

  void processMCParticleLevel(soa::Filtered<aod::ChargedMCParticleLevelJets>::iterator const& jet)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCParticleLevel, "jets on particle level MC", false);

  void processMCCharged(aod::JCollisions const& collision, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCCharged, "jets on detector and particle level MC", false);

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processMCChargedMatched(aod::JCollision const& collision,
                               soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                               JetMCPTable const& mcpjets,
                               aod::JTracks const& tracks,
                               aod::JMcParticles const& particles)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCChargedMatched, "jet finder QA matched mcp and mcd", false);

  void processDataChargedSubstructure(soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>::iterator const& jet, aod::JTracks const& tracks)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataChargedSubstructure, "jet substructure charged jets", false);

  void processMCParticleSubstructure(soa::Filtered<soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>>::iterator const& jet, aod::JMcParticles const& particles)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCParticleSubstructure, "jet substructure particle level full jets", false);

  void processDataFull(soa::Filtered<aod::FullJets>::iterator const& jet)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataFull, "jets data", true);

  void processDataFullSubstructure(soa::Filtered<soa::Join<aod::FullJets, aod::FullJetConstituents>>::iterator const& jet, aod::JTracks const& tracks, aod::JClusters const& clusters)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataFullSubstructure, "jet substructure full jets", false);

  void processDataRecoil(aod::JCollision const& collision, soa::Filtered<aod::ChargedJets> const& jets, aod::JTracks const& tracks)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataRecoil, "hadron-recoil jets", false);

  void processDataBackgroundSubtracted(soa::Join<aod::JCollisions, aod::JCollisionRhos>::iterator const& collision, soa::Filtered<aod::ChargedJets> const& jets)
  {
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataBackgroundSubtracted, "baackground subtracted jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetTutorialSkeletonTask>(cfgc, TaskName{"jet-tutorial-skeleton"})}; }
