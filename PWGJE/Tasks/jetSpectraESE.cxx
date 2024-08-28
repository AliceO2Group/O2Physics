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

// jet analysis framework with ESE (19/08/2024)
//
/// \author Joachim Hansen <joachim.hansen@cern.ch>
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

#include "PWGCF/Flow/DataModel/FlowESE.h"
#include "PWGCF/Flow/Core/FFitWeights.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

using ColWqVecFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::qVecFT0Cs, aod::qPercentileFT0Cs>;
using JColwESE = soa::Join<aod::JCollisions, aod::CentFT0Cs, aod::qVecFT0Cs, aod::qPercentileFT0Cs>;

struct JetSpectraEseTask {

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
                              {"h_matched_jets_phi", "#phi_{jet part}; #phi_{jet det}", {HistType::kTH2F, {{80, -1.0, 7.}, {80, -1.0, 7.}}}},

                              {"h_Cent_Psi2_FV0A", "#Psi_{2};Centrality", {HistType::kTH2F, {{90, 0, 90}, {150, -3, 3}}}},
                              {"h_Cent_Qx_FT0C", "Qx_{2};Centrality", {HistType::kTH2F, {{90, 0, 90}, {250, -250, 250}}}},
                              {"h_Cent_Qy_FV0A", "Qy_{2};Centrality;", {HistType::kTH2F, {{90, 0, 90}, {250, -250, 250}}}},
                              {"h_Cent_q2Percs", "qPercs;Centrality;", {HistType::kTH2F, {{100, 0, 100}, {250, 0, 5000}}}},

                              {"h_jet_pt_q2_0", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_pt_q2_1", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_pt_in", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_pt_out", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}}

                             }};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<std::vector<float>> qLimits{"qLimits", {30, 70}, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);

  void processCollisions(JetCollision const& collision, JetTracks const& tracks)
  {

    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processCollisions, "process self contained collisions", true);

  int spCode(const int qPerc)
  {
    if (qPerc < 0)
      return -1;

    if (qPerc <= qLimits->at(0)) {
      return 0;
    } else if (qPerc >= qLimits->at(1)) {
      return 1;
    } else {
      return -1;
    }
  }
  std::string jetplane(float deltaPhi)
  {
    std::string s = "";
    if (deltaPhi < TMath::Pi() / 6.) {
      s = "in";
    } else if (deltaPhi > TMath::Pi() / 3.) {
      s = "out";
    }

    return s;
  }

  void processESEDataCharged(JColwESE::iterator const& collision, soa::Filtered<aod::ChargedJets> const& jets, aod::BCsWithTimestamps const&)
  {

    float vPsi2 = FFitWeights::EventPlane(collision, 2);

    auto qPerc = collision.qPERCFT0C();
    if (qPerc[0] > 0)
      registry.fill(HIST("h_Cent_q2Percs"), collision.centFT0C(), qPerc[0]);

    int code = spCode(qPerc[0]);

    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    registry.fill(HIST("h_collisions"), 1.5);
    float plane{0};
    std::string pl{""};
    for (auto const& jet : jets) {

      plane = vPsi2 - jet.phi();
      pl = jetplane(plane);
      registry.fill(HIST("h_jet_pt"), jet.pt());

      if (code < 0)
        continue;

      if (code == 0)
        registry.fill(HIST("h_jet_pt_q2_0"), jet.pt());
      if (code == 1)
        registry.fill(HIST("h_jet_pt_q2_1"), jet.pt());
      if (pl == "in")
        registry.fill(HIST("h_jet_pt_in"), jet.pt());
      if (pl == "out")
        registry.fill(HIST("h_jet_pt_out"), jet.pt());
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEDataCharged, "process self contained collisions", true);

  void processCollisionsWithExternalTracks(JetCollision const& collision, soa::Join<JetTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection> const&)
  {

    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
      auto originalTrack = track.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>();
      registry.fill(HIST("h_track_chi2PerCluster"), originalTrack.tpcChi2NCl());
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processCollisionsWithExternalTracks, "process non self contained collisions", true);

  void processDataCharged(soa::Filtered<aod::ChargedJets>::iterator const& jet)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processDataCharged, "jets data", true);

  void processMCDetectorLevelCharged(soa::Filtered<aod::ChargedMCDetectorLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCDetectorLevelCharged, "jets on detector level MC", false);

  void processMCParticleLevel(soa::Filtered<aod::ChargedMCParticleLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_part_jet_pt"), jet.pt());
    registry.fill(HIST("h_part_jet_eta"), jet.eta());
    registry.fill(HIST("h_part_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCParticleLevel, "jets on particle level MC", false);

  void processMCCharged(JetCollision const&, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets)
  {
    for (auto& mcdjet : mcdjets) {
      registry.fill(HIST("h_jet_pt"), mcdjet.pt());
      registry.fill(HIST("h_jet_eta"), mcdjet.eta());
      registry.fill(HIST("h_jet_phi"), mcdjet.phi());
    }
    for (auto& mcpjet : mcpjets) {
      registry.fill(HIST("h_part_jet_pt"), mcpjet.pt());
      registry.fill(HIST("h_part_jet_eta"), mcpjet.eta());
      registry.fill(HIST("h_part_jet_phi"), mcpjet.phi());
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCCharged, "jets on detector and particle level MC", false);

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processMCChargedMatched(JetCollision const&,
                               soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                               JetMCPTable const&,
                               JetTracks const&,
                               JetParticles const&)
  {
    for (const auto& mcdjet : mcdjets) {

      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetMCPTable>()) {
        // for (auto& mcpjet : mcdjet.template matchedJetPt_as<JetMCPTable>()) {

        registry.fill(HIST("h_matched_jets_pt"), mcpjet.pt(), mcdjet.pt());
        registry.fill(HIST("h_matched_jets_pt"), mcpjet.phi(), mcdjet.phi());
        registry.fill(HIST("h_matched_jets_pt"), mcpjet.eta(), mcdjet.eta());
      }
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCChargedMatched, "jet finder QA matched mcp and mcd", false);

  void processDataChargedSubstructure(soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>::iterator const& jet, JetTracks const&)
  {
    // add aditional selection on jet eta
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
    registry.fill(HIST("h_jet_ntracks"), jet.tracksIds().size());
    double angularity = 0.0;
    for (auto& jetConstituent : jet.tracks_as<JetTracks>()) {
      angularity += jetConstituent.pt() * TMath::Sqrt(TMath::Power(jet.phi() - jetConstituent.phi(), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0));
    }
    registry.fill(HIST("h_jet_angularity"), angularity / (jet.pt() * round(jet.r() / 100.0f)));
  }
  PROCESS_SWITCH(JetSpectraEseTask, processDataChargedSubstructure, "jet substructure charged jets", false);

  void processMCParticleSubstructure(soa::Filtered<soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>>::iterator const& jet, JetParticles const&)
  {
    double angularity = 0.0;
    for (auto& jetConstituents : jet.tracks_as<aod::JMcParticles>()) {
      angularity += jetConstituents.pt() * TMath::Sqrt(TMath::Power(jet.phi() - jetConstituents.phi(), 2.0) + TMath::Power(jet.eta() - jetConstituents.eta(), 2.0));
    }
    registry.fill(HIST("h_part_jet_angularity"), angularity / (jet.pt() * round(jet.r() / 100.0f)));
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCParticleSubstructure, "jet substructure particle level full jets", false);

  void processDataFull(soa::Filtered<aod::FullJets>::iterator const& jet)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processDataFull, "jets data", true);

  void processDataFullSubstructure(soa::Filtered<soa::Join<aod::FullJets, aod::FullJetConstituents>>::iterator const& jet, JetTracks const&, JetClusters const&)
  {
    // add aditional selection on jet eta
    registry.fill(HIST("h_full_jet_pt"), jet.pt());
    registry.fill(HIST("h_full_jet_eta"), jet.eta());
    registry.fill(HIST("h_full_jet_phi"), jet.phi());
    registry.fill(HIST("h_full_jet_ntracks"), jet.tracksIds().size());
    registry.fill(HIST("h_full_jet_nclusters"), jet.clustersIds().size());
    double angularity = 0.0;
    for (auto& jetTrack : jet.tracks_as<JetTracks>()) {
      angularity += jetTrack.pt() * TMath::Sqrt(TMath::Power(jet.phi() - jetTrack.phi(), 2.0) + TMath::Power(jet.eta() - jetTrack.eta(), 2.0));
    }

    for (auto& jetCluster : jet.clusters_as<JetClusters>()) {
      angularity += jetCluster.energy() * TMath::Sqrt(TMath::Power(jet.phi() - jetCluster.phi(), 2.0) + TMath::Power(jet.eta() - jetCluster.eta(), 2.0));
    }

    registry.fill(HIST("h_full_jet_angularity"), angularity / (jet.pt() * round(jet.r() * 100.0f)));
  }
  PROCESS_SWITCH(JetSpectraEseTask, processDataFullSubstructure, "jet substructure full jets", false);

  void processDataRecoil(JetCollision const& collision, soa::Filtered<aod::ChargedJets> const& jets, JetTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    double leadingTrackpT = 0.0;
    double leadingTrackPhi = 0.0;
    for (auto& track : tracks) {
      if (track.pt() > 6.0 && track.pt() < 10.0) {
        if (track.pt() > leadingTrackpT) {
          leadingTrackpT = track.pt();
          leadingTrackPhi = track.phi();
        }
      }
    }
    if (leadingTrackpT == 0.0)
      return;
    for (auto& jet : jets) {
      if (TMath::Abs(RecoDecay::constrainAngle(RecoDecay::constrainAngle(jet.phi(), -o2::constants::math::PIHalf) - RecoDecay::constrainAngle(leadingTrackPhi, -o2::constants::math::PIHalf), -o2::constants::math::PIHalf) > 0.6)) {
        registry.fill(HIST("h_recoil_jet_pt"), jet.pt());
        registry.fill(HIST("h_recoil_jet_eta"), jet.eta());
        registry.fill(HIST("h_recoil_jet_phi"), jet.phi());
        registry.fill(HIST("h_recoil_jet_dphi"), jet.phi() - leadingTrackPhi);
      }
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processDataRecoil, "hadron-recoil jets", false);

  // void processDataBackgroundSubtracted(soa::Join<JetCollisions, aod::JCollisionRhos>::iterator const& collision, soa::Filtered<aod::ChargedJets> const& jets)
  // {
  //   for (auto jet : jets) {
  //     registry.fill(HIST("h_jet_pt"), jet.pt());
  //     registry.fill(HIST("h_jet_pt_bkgsub"), jet.pt() - (collision.rho() * jet.area()));
  //     registry.fill(HIST("h_jet_eta"), jet.eta());
  //     registry.fill(HIST("h_jet_phi"), jet.phi());
  //   }
  // }
  // PROCESS_SWITCH(JetSpectraEseTask, processDataBackgroundSubtracted, "baackground subtracted jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetSpectraEseTask>(cfgc, TaskName{"jet-spectra-ese"})}; }
