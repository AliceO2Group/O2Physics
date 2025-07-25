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

// jet finder v0 QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TMathBase.h>

#include <cmath>
#include <string>
#include <type_traits>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename JetConstituentTableData, typename JetMatchingTableData, typename CandidateTableData, typename JetTableMCD, typename JetConstituentTableMCD, typename JetMatchingTableMCDMCP, typename JetTableMCDWeighted, typename CandidateTableMCD, typename JetTableMCP, typename JetConstituentTableMCP, typename JetMatchingTableMCPMCD, typename JetTableMCPWeighted, typename CandidateTableMCP>
struct JetFinderV0QATask {
  HistogramRegistry registry;

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    auto jetRadiiBins = (std::vector<double>)jetRadii;
    if (jetRadiiBins.size() > 1) {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (TMath::Abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
    } else {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
    }

    if (doprocessJetsData || doprocessJetsMCD || doprocessJetsMCDWeighted) {
      registry.add("h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h2_centrality_jet_pt", "centrality vs #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, 0., 200.}}});
      registry.add("h2_centrality_jet_eta", "centrality vs #eta_{jet}; centrality; #eta_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_jet_phi", "centrality vs #varphi_{jet}; centrality; #varphi_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});
      registry.add("h2_centrality_jet_ntracks", "centrality vs N_{jet tracks}; centrality; N_{jet tracks}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_centrality", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});centrality", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {1200, -10.0, 110.0}}});
      registry.add("h3_jet_r_jet_pt_jet_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_jet_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_eta_jet_phi", "#it{R}_{jet};#eta_{jet};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_jet_ntracks", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_jet_area", "#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c}); #it{area}_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {300, 0., 3.}}});
      registry.add("h3_jet_r_jet_pt_track_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_candidate_invmass_jet_pt_candidate_pt", ";#it{m}_{inv, candidate} (GeV/#it{c}^{2}); #it{p}_{T,jet} (GeV/#it{c}) ;#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{500, 0.0, 5.0}, {200, 0.0, 200.0}, {200, 0.0, 200.0}}});
    }

    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      registry.add("h_jet_pt_part", "jet pT;#it{p}_{T,jet}^{part}(GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta_part", "jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_jet_phi_part", "jet #varphi;#varphi_{jet}^{part};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_jet_ntracks_part", "jet N tracks;N_{jet tracks}^{part};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_eta_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_phi_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_eta_part_jet_phi_part", ";#it{R}_{jet}^{part};#eta_{jet}^{part};#varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_ntracks_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_pt_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet tracks}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_eta_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_phi_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_pt_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,candidate}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_eta_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{candidate}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_phi_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi{candidate}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_y_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});y_{candidate}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    }

    if (doprocessTracks || doprocessTracksWeighted) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_collisions", "centrality vs collisions; centrality; collisions", {HistType::kTH2F, {{1200, -10.0, 110.0}, {4, 0.0, 4.0}}});
      registry.add("h2_centrality_track_pt", "centrality vs track pT; centrality; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, 0., 200.}}});
      registry.add("h2_centrality_track_eta", "centrality vs track #eta; centrality; #eta_{track}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_track_phi", "centrality vs track #varphi; centrality; #varphi_{track}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});
      registry.add("h2_centrality_track_energy", "centrality vs track energy; centrality; Energy GeV", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, 0.0, 100.0}}});
      if (doprocessTracksWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }
    if (doprocessMCCollisionsWeighted) {
      AxisSpec weightAxis = {{VARIABLE_WIDTH, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}, "weights"};
      registry.add("h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {weightAxis}});
    }

    if (doprocessCandidates) {
      registry.add("h_candidate_invmass", ";#it{m}_{inv, candidate} (GeV/#it{c}^{2});counts", {HistType::kTH1F, {{500, 0.0, 5.0}}});
      registry.add("h_candidate_pt", ";#it{p}_{T,candidate} (GeV/#it{c});counts", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      registry.add("h_candidate_y", ";y_{candidate};counts", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h2_centrality_ncandidates", "centrality vs N_{candidates};centrality;N_{candidates};", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, 0.0, 100.0}}});
    }
  }

  using JetTableDataJoined = soa::Join<JetTableData, JetConstituentTableData>;
  using JetTableDataMatchedJoined = soa::Join<JetTableData, JetConstituentTableData, JetMatchingTableData>;

  using JetTableMCDJoined = soa::Join<JetTableMCD, JetConstituentTableMCD>;
  using JetTableMCDWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetTableMCDWeighted>;
  using JetTableMCDMatchedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP>;
  using JetTableMCDMatchedWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP, JetTableMCDWeighted>;

  using JetTableMCPJoined = soa::Join<JetTableMCP, JetConstituentTableMCP>;
  using JetTableMCPWeightedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetTableMCPWeighted>;
  using JetTableMCPMatchedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetMatchingTableMCPMCD>;
  using JetTableMCPMatchedWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCP, JetMatchingTableMCPMCD, JetTableMCPWeighted>;

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centFT0M >= centralityMin && aod::jcollision::centFT0M < centralityMax);

  template <typename T, typename U, typename V>
  bool isAcceptedJet(V const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    if (leadingConstituentPtMin > -98.0) {
      bool isMinleadingConstituent = false;
      for (auto& constituent : jet.template tracks_as<T>()) {
        if (constituent.pt() >= leadingConstituentPtMin) {
          isMinleadingConstituent = true;
          break;
        }
      }
      for (auto& hfconstituent : jet.template candidates_as<U>()) {
        if (hfconstituent.pt() >= leadingConstituentPtMin) {
          isMinleadingConstituent = true;
          break;
        }
      }
      if (!isMinleadingConstituent) {
        return false;
      }
    }
    return true;
  }

  template <typename T, typename U>
  void fillHistograms(T const& jet, float centrality, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }

    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracksIds().size() + jet.candidatesIds().size(), weight);
      registry.fill(HIST("h2_centrality_jet_pt"), centrality, jet.pt(), weight);
      registry.fill(HIST("h2_centrality_jet_eta"), centrality, jet.eta(), weight);
      registry.fill(HIST("h2_centrality_jet_phi"), centrality, jet.phi(), weight);
      registry.fill(HIST("h2_centrality_jet_ntracks"), centrality, jet.tracksIds().size() + jet.candidatesIds().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_jet_pt_centrality"), jet.r() / 100.0, jet.pt(), centrality, weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_eta"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_phi"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_eta_jet_phi"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_ntracks"), jet.r() / 100.0, jet.pt(), jet.tracksIds().size() + jet.candidatesIds().size(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_area"), jet.r() / 100.0, jet.pt(), jet.area(), weight);

    for (auto& constituent : jet.template tracks_as<aod::JetTracks>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }

    for (auto& v0candidate : jet.template candidates_as<std::decay_t<U>>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt"), jet.r() / 100.0, jet.pt(), v0candidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta"), jet.r() / 100.0, jet.pt(), v0candidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi"), jet.r() / 100.0, jet.pt(), v0candidate.phi(), weight);

      registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt"), jet.r() / 100.0, jet.pt(), v0candidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta"), jet.r() / 100.0, jet.pt(), v0candidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi"), jet.r() / 100.0, jet.pt(), v0candidate.phi(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_y"), jet.r() / 100.0, jet.pt(), v0candidate.y(), weight);

      if (jet.r() == round(selectedJetsRadius * 100.0f)) {
        registry.fill(HIST("h3_candidate_invmass_jet_pt_candidate_pt"), v0candidate.mLambda(), jet.pt(), v0candidate.pt(), weight);
        registry.fill(HIST("h3_candidate_invmass_jet_pt_candidate_pt"), v0candidate.mAntiLambda(), jet.pt(), v0candidate.pt(), weight);
        registry.fill(HIST("h3_candidate_invmass_jet_pt_candidate_pt"), v0candidate.mK0Short(), jet.pt(), v0candidate.pt(), weight);
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat) {
      return;
    }

    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks_part"), jet.tracksIds().size() + jet.candidatesIds().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_eta_part"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_phi_part"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_eta_part_jet_phi_part"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_ntracks_part"), jet.r() / 100.0, jet.pt(), jet.tracksIds().size() + jet.candidatesIds().size(), weight);

    for (auto& constituent : jet.template tracks_as<std::decay_t<U>>()) {

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }

    for (auto& v0candidate : jet.template candidates_as<std::decay_t<V>>()) {

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), v0candidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), v0candidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), v0candidate.phi(), weight);

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_pt_part"), jet.r() / 100.0, jet.pt(), v0candidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_eta_part"), jet.r() / 100.0, jet.pt(), v0candidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_phi_part"), jet.r() / 100.0, jet.pt(), v0candidate.phi(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_y_part"), jet.r() / 100.0, jet.pt(), v0candidate.y(), weight);
    }
  }

  template <typename T, typename U>
  void fillTrackHistograms(T const& collision, U const& tracks, float weight = 1.0)
  {
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h2_centrality_track_pt"), collision.centFT0M(), track.pt(), weight);
      registry.fill(HIST("h2_centrality_track_eta"), collision.centFT0M(), track.eta(), weight);
      registry.fill(HIST("h2_centrality_track_phi"), collision.centFT0M(), track.phi(), weight);
      registry.fill(HIST("h2_centrality_track_energy"), collision.centFT0M(), track.energy(), weight);
    }
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetFinderV0QATask, processDummy, "dummy task", true);

  void processJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision, JetTableDataJoined const& jets, CandidateTableData const&, aod::JetTracks const&)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks, CandidateTableData>(jet)) {
        continue;
      }
      fillHistograms<typename JetTableDataJoined::iterator, CandidateTableData>(jet, collision.centFT0M());
    }
  }
  PROCESS_SWITCH(JetFinderV0QATask, processJetsData, "jet finder HF QA data", false);

  void processJetsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, JetTableMCDJoined const& jets, CandidateTableMCD const&, aod::JetTracks const&)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks, CandidateTableMCD>(jet)) {
        continue;
      }
      fillHistograms<typename JetTableMCDJoined::iterator, CandidateTableMCD>(jet, collision.centFT0M());
    }
  }
  PROCESS_SWITCH(JetFinderV0QATask, processJetsMCD, "jet finder HF QA mcd", false);

  void processJetsMCDWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision, JetTableMCDWeightedJoined const& jets, CandidateTableMCD const&, aod::JetTracks const&)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks, CandidateTableMCD>(jet)) {
        continue;
      }
      fillHistograms<typename JetTableMCDWeightedJoined::iterator, CandidateTableMCD>(jet, collision.centFT0M(), jet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderV0QATask, processJetsMCDWeighted, "jet finder HF QA mcd on weighted events", false);

  void processJetsMCP(typename JetTableMCPJoined::iterator const& jet, aod::JetParticles const&, CandidateTableMCP const&)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    if (!isAcceptedJet<aod::JetParticles, CandidateTableMCP>(jet)) {
      return;
    }
    fillMCPHistograms<typename JetTableMCPJoined::iterator, aod::JetParticles, CandidateTableMCP>(jet);
  }
  PROCESS_SWITCH(JetFinderV0QATask, processJetsMCP, "jet finder HF QA mcp", false);

  void processJetsMCPWeighted(typename JetTableMCPWeightedJoined::iterator const& jet, aod::JetParticles const&, CandidateTableMCP const&)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    if (!isAcceptedJet<aod::JetParticles, CandidateTableMCP>(jet)) {
      return;
    }
    fillMCPHistograms<typename JetTableMCPWeightedJoined::iterator, aod::JetParticles, CandidateTableMCP>(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderV0QATask, processJetsMCPWeighted, "jet finder HF QA mcp on weighted events", false);

  void processMCCollisionsWeighted(aod::JetMcCollision const& collision)
  {
    registry.fill(HIST("h_collision_eventweight_part"), collision.weight());
  }
  PROCESS_SWITCH(JetFinderV0QATask, processMCCollisionsWeighted, "collision QA for weighted events", false);

  void processTracks(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                     soa::Filtered<aod::JetTracks> const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centFT0M(), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centFT0M(), 1.5);
    fillTrackHistograms(collision, tracks);
  }
  PROCESS_SWITCH(JetFinderV0QATask, processTracks, "QA for charged tracks", false);

  void processTracksWeighted(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                             aod::JetMcCollisions const&,
                             soa::Filtered<aod::JetTracks> const& tracks)
  {
    float eventWeight = collision.mcCollision().weight();
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    fillTrackHistograms(collision, tracks, eventWeight);
  }
  PROCESS_SWITCH(JetFinderV0QATask, processTracksWeighted, "QA for charged tracks weighted", false);

  void processCandidates(soa::Filtered<aod::JetCollisions>::iterator const& collision, CandidateTableData const& candidates)
  {

    for (auto const& candidate : candidates) {
      registry.fill(HIST("h_candidate_invmass"), candidate.mLambda());
      registry.fill(HIST("h_candidate_invmass"), candidate.mAntiLambda());
      registry.fill(HIST("h_candidate_invmass"), candidate.mK0Short());
      registry.fill(HIST("h_candidate_pt"), candidate.pt());
      registry.fill(HIST("h_candidate_y"), candidate.y());
    }
    registry.fill(HIST("h2_centrality_ncandidates"), collision.centFT0M(), candidates.size());
  }
  PROCESS_SWITCH(JetFinderV0QATask, processCandidates, "HF candidate QA", false);
};

using JetFinderV0QA = JetFinderV0QATask<aod::V0ChargedJets, aod::V0ChargedJetConstituents, aod::V0ChargedJetsMatchedToV0ChargedEventWiseSubtractedJets, aod::CandidatesV0Data, aod::V0ChargedMCDetectorLevelJets, aod::V0ChargedMCDetectorLevelJetConstituents, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets, aod::V0ChargedMCDetectorLevelJetEventWeights, aod::CandidatesV0MCD, aod::V0ChargedMCParticleLevelJets, aod::V0ChargedMCParticleLevelJetConstituents, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets, aod::V0ChargedMCParticleLevelJetEventWeights, aod::CandidatesV0MCP>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderV0QA>(cfgc,
                                                      SetDefaultProcesses{},
                                                      TaskName{"jet-finder-charged-v0-qa"}));

  return WorkflowSpec{tasks};
}
