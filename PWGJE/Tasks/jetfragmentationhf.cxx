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

// \file jetfragmentationhf.cxx
// \note Extended from jetfinderhfQA.cxx and jetfragmentation.cxx

// HF jets fragmentation function task
//
// Authors: Hanseo Park

#include <string>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;

#include "Framework/runDataProcessing.h"

struct JetFragmentationHFTask {

  Configurable<double> yCandMax{"yCandMax", -1, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};

  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {200, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binNtracks{"binNtracks", {100, -0.5, 99.5}, ""};
  ConfigurableAxis binZ{"binZ", {100, -5e-3f, 1.f + 5e-3f}, ""};
  ConfigurableAxis binJetR{"binJetR", {6, 0.05f, 0.65f}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binD0CandidatePt{"binD0CandidatePt", {360, 0.f, 36.f}, ""};
  ConfigurableAxis binD0CandidateMass{"binD0CandidatMass", {500, 0.f, 5.f}, ""};

  // Axis
  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T, jet}"};
  AxisSpec etaAxis = {binEta, "#eta"};
  AxisSpec phiAxis = {binPhi, "#phi"};
  AxisSpec ntracksAxis = {binNtracks, "#it{N}_{tracks}"};
  AxisSpec zAxis = {binZ, "#it{z}"};
  AxisSpec jetRAxis = {binJetR, "#it{R}"};
  AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{track}"};
  AxisSpec d0candidatePtAxis = {binD0CandidatePt, "#it{p}_{T}^{D0 candidate}"};
  AxisSpec d0candidateMassAxis = {binD0CandidateMass, "inv. mass(#pi K)"};

  HistogramRegistry registry{"registry",
                             {{"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxis}}, true},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {etaAxis}}, true},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {phiAxis}}, true},
                              {"h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {ntracksAxis}}, true},
                              {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {trackPtAxis}}, true},
                              {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {etaAxis}}, true},
                              {"h_track_phi", "track #phi;#phi_{track};entries", {HistType::kTH1F, {phiAxis}}, true},
                              {"h2_jet_pt_track_pt", ";#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, trackPtAxis}}, true},
                              {"h2_jet_pt_track_eta", ";#it{p}_{T,jet} (GeV/#it{c}); #eta_{track}", {HistType::kTH2F, {jetPtAxis, etaAxis}}, true},
                              {"h2_jet_pt_track_phi", ";#it{p}_{T,jet} (GeV/#it{c}); #phi_{track}", {HistType::kTH2F, {jetPtAxis, phiAxis}}, true},
                              {"h_part_jet_pt", "Particle Level jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxis}}, true},
                              {"h_part_jet_eta", "Particle Level jet #eta;#eta_{jet};entries", {HistType::kTH1F, {etaAxis}}, true},
                              {"h_part_jet_phi", "Particle Level jet #phi;#phi_{jet};entries", {HistType::kTH1F, {phiAxis}}, true},
                              {"h_part_jet_ntracks", "Particle Level jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {ntracksAxis}}, true},
                              {"h_part_track_pt", "Particle Level track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {trackPtAxis}}, true},
                              {"h_part_track_eta", "Particle Level track #eta;#eta_{track};entries", {HistType::kTH1F, {etaAxis}}, true},
                              {"h_part_track_phi", "Particle Level track #phi;#phi_{track};entries", {HistType::kTH1F, {phiAxis}}, true},
                              {"h2_part_jet_pt_part_track_pt", "Particle Level;#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, trackPtAxis}}, true},
                              {"h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {{100, 0.0, 10.0}}}, true}}};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    // Data and MC (Detector level, Rec)
    registry.add("h_d0candidate_pt", "D0 candidates #it{p}_{T};#it{p}_{T,candidate} (GeV/#it{c});entries", {HistType::kTH1F, {d0candidatePtAxis}});
    registry.add("h_d0candidate_eta", "D0 candidates #eta;#eta_{candidate};entries", {HistType::kTH1F, {etaAxis}});
    registry.add("h_d0candidate_phi", "D0 candidates #phi;#phi_{candidate};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("h2_jet_pt_d0candidate_pt", "#it{p}_{T, jet} vs #it{t}_{T, D0 candidates};#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, d0candidatePtAxis}});
    registry.add("h2_jet_pt_d0candidate_eta", "#it{p}_{T,jet} vs #eta_{D0 candidates};#it{p}_{T,jet} (GeV/#it{c}); #eta_{candidate}", {HistType::kTH2F, {jetPtAxis, etaAxis}});
    registry.add("h2_jet_pt_d0candidate_phi", "#it{p}_{T,jet} vs #phi_{D0 candidates};#it{p}_{T,jet} (GeV/#it{c}); #phi_{candidate}", {HistType::kTH2F, {jetPtAxis, phiAxis}});
    registry.add("h_d0candidate_mass", "Inv. mass distribution of D0 candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {d0candidateMassAxis}});
    registry.add("h2_d0candidate_mass_d0candidate_pt", "Inv. mass of D0 candidates vs #it{p}_{T, D0 candidates};inv. mass (#pi K) (GeV/#it{c}^{2});entries ", {HistType::kTH2F, {d0candidateMassAxis, {vbins, "#it{p}_{T,candidate} (GeV/#it{c})"}}});
    registry.add("h2_d0candidate_mass_d0candidate_phi", "Inv. mass of D0 candidates vs #phi_{D0 candidates};inv. mass (#pi K) (GeV/#it{c}^{2});entries ", {HistType::kTH2F, {d0candidateMassAxis, phiAxis}});
    registry.add("h_jet_chargefrag", "#it{p}_{T,track}/#it{p}_{T,jet};entries", {HistType::kTH1F, {zAxis}});
    registry.add("h_jet_d0frag", "#it{p}_{T,D0 candidates}/#it{p}_{T,jet};entries", {HistType::kTH1F, {zAxis}});
    // prong
    registry.add("h_d0candidate_prong0_pt", "D0 candidates prong0#it{p}_{T};#it{p}_{T,candidate} (GeV/#it{c});entries", {HistType::kTH1F, {d0candidatePtAxis}});
    registry.add("h_d0candidate_prong1_pt", "D0 candidates prong1#it{p}_{T};#it{p}_{T,candidate} (GeV/#it{c});entries", {HistType::kTH1F, {d0candidatePtAxis}});
    registry.add("h_decay_length", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} GeV/#it{c}"}}});
    registry.add("h_decay_length_xy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} GeV/#it{c}"}}});
    registry.add("h_decay_length_error", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_decay_length_xy_error", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_impact_parameter_d0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_impact_parameter_d0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_impact_parameter_d0Prong0_error", "2-prong candidates;prong 0 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_impact_parameter_d0Prong1_error", "2-prong candidates;prong 1 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_jet_d0Prong0frag", "#it{p}_{T,D0 candidates Prong0}/#it{p}_{T,jet};entries", {HistType::kTH1F, {zAxis}});
    registry.add("h_jet_d0Prong1frag", "#it{p}_{T,D0 candidates Prong1}/#it{p}_{T,jet};entries", {HistType::kTH1F, {zAxis}});

    // MC (Partcle level, Gen)
    registry.add("h_part_d0candidate_pt", "Particle Level D0 candidates #it{p}_{T};#it{p}_{T,candidate} (GeV/#it{c}) ();entries", {HistType::kTH1F, {d0candidatePtAxis}});
    registry.add("h_part_d0candidate_eta", "Particle Level D0 candidates #eta;#eta_{candidate};entries", {HistType::kTH1F, {etaAxis}});
    registry.add("h_part_d0candidate_phi", "Particle Level D0 candidates #phi;#phi_{candidate};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("h2_part_jet_pt_part_d0candidate_pt", "#it{p}_{T, jet}^{Part} vs #it{t}_{T, D0 candidates}^{Part};#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, d0candidatePtAxis}});
    registry.add("h2_part_jet_pt_part_d0candidate_eta", "#it{p}_{T,jet}^{Part} vs #eta_{D0 candidates}^{Part};#it{p}_{T,jet} (GeV/#it{c}); #eta_{candidate}", {HistType::kTH2F, {jetPtAxis, etaAxis}});
    registry.add("h2_part_jet_pt_part_d0candidate_phi", "#it{p}_{T,jet}^{Part} vs #phi_{D0 candidates}^{Part};#it{p}_{T,jet} (GeV/#it{c}); #phi_{candidate}", {HistType::kTH2F, {jetPtAxis, phiAxis}});
    registry.add("h_part_jet_part_chargefrag", "#it{p}_{T,track}^{Part}/#it{p}_{T,jet}^{Part};entries", {HistType::kTH1F, {zAxis}});
    registry.add("h_part_jet_part_d0frag", "#it{p}_{T,D0 candidates}^{Part}/#it{p}_{T,jet}^{Part};entries", {HistType::kTH1F, {zAxis}});
    registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);D^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);D^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    // Matched
    registry.add("h2_jet_pt_part_jet_pt", "#it{p}_{T,jet}^{part} vs #it{p}_{T,jet}^{Dec} ;#it{p}_{T,jet}^{DeC} (GeV/#it{c}); #it{p}_{T,jet}^{part}", {HistType::kTH2F, {binJetPt, binJetPt}});

    // Comp. hf jets and inclusive jets
  }

  using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using CandidateD0Data = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using CandidateD0MC = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>;
  using JetParticles2Prong = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  template <typename T>
  void fillDataHistograms(T const& jets, float weight = 1.0)
  {
    for (const auto& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);
      for (auto& track : jet.template tracks_as<JetTracks>()) {
        auto chargeFrag = track.pt() / jet.pt();
        registry.fill(HIST("h_jet_chargefrag"), chargeFrag, weight);
        registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), track.pt(), weight);
        registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), track.eta(), weight);
        registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), track.phi(), weight);
        registry.fill(HIST("h_track_pt"), track.pt(), weight);
        registry.fill(HIST("h_track_eta"), track.eta(), weight);
        registry.fill(HIST("h_track_phi"), track.phi(), weight);
      }
      for (auto& d0candidate : jet.template hfcandidates_as<CandidateD0Data>()) {
        auto massD0 = invMassD0ToPiK(d0candidate);
        auto D0Frag = d0candidate.pt() / jet.pt();
        auto D0Prong0Frag = d0candidate.ptProng0() / jet.pt();
        auto D0Prong1Frag = d0candidate.ptProng1() / jet.pt();
        registry.fill(HIST("h_jet_d0frag"), D0Frag, weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_pt"), jet.pt(), d0candidate.pt(), weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_eta"), jet.pt(), d0candidate.eta(), weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_phi"), jet.pt(), d0candidate.phi(), weight);
        registry.fill(HIST("h_d0candidate_pt"), d0candidate.pt(), weight);
        registry.fill(HIST("h_d0candidate_eta"), d0candidate.eta(), weight);
        registry.fill(HIST("h_d0candidate_phi"), d0candidate.phi(), weight);

        registry.fill(HIST("h_jet_d0Prong0frag"), D0Prong0Frag, weight);
        registry.fill(HIST("h_jet_d0Prong1frag"), D0Prong1Frag, weight);
        registry.fill(HIST("h_d0candidate_prong0_pt"), d0candidate.ptProng0(), weight);
        registry.fill(HIST("h_d0candidate_prong1_pt"), d0candidate.ptProng1(), weight);
        registry.fill(HIST("h_decay_length"), d0candidate.decayLength(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_xy"), d0candidate.decayLengthXY(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_error"), d0candidate.errorDecayLength(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_xy_error"), d0candidate.errorDecayLengthXY(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong0"), d0candidate.impactParameter0(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong1"), d0candidate.impactParameter1(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong0_error"), d0candidate.errorImpactParameter0(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong1_error"), d0candidate.errorImpactParameter1(), d0candidate.pt(), weight);

        if (d0candidate.isSelD0() >= 1 || d0candidate.isSelD0bar() >= 1) {
          registry.fill(HIST("h_d0candidate_mass"), massD0);
          registry.fill(HIST("h2_d0candidate_mass_d0candidate_pt"), massD0, d0candidate.pt());
          registry.fill(HIST("h2_d0candidate_mass_d0candidate_phi"), massD0, d0candidate.phi());
        }
      }
    }
  }

  template <typename T>
  void fillMCDHistograms(T const& jets, float weight = 1.0)
  {
    for (const auto& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);
      for (auto& track : jet.template tracks_as<JetTracks>()) {
        auto chargeFrag = track.pt() / jet.pt();
        registry.fill(HIST("h_jet_chargefrag"), chargeFrag, weight);
        registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), track.pt(), weight);
        registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), track.eta(), weight);
        registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), track.phi(), weight);
        registry.fill(HIST("h_track_pt"), track.pt(), weight);
        registry.fill(HIST("h_track_eta"), track.eta(), weight);
        registry.fill(HIST("h_track_phi"), track.phi(), weight);
      }
      for (auto& d0candidate : jet.template hfcandidates_as<CandidateD0MC>()) {
        auto massD0 = invMassD0ToPiK(d0candidate);
        auto D0Frag = d0candidate.pt() / jet.pt();
        auto D0Prong0Frag = d0candidate.ptProng0() / jet.pt();
        auto D0Prong1Frag = d0candidate.ptProng1() / jet.pt();
        registry.fill(HIST("h_jet_d0frag"), D0Frag, weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_pt"), jet.pt(), d0candidate.pt(), weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_eta"), jet.pt(), d0candidate.eta(), weight);
        registry.fill(HIST("h2_jet_pt_d0candidate_phi"), jet.pt(), d0candidate.phi(), weight);
        registry.fill(HIST("h_d0candidate_pt"), d0candidate.pt(), weight);
        registry.fill(HIST("h_d0candidate_eta"), d0candidate.eta(), weight);
        registry.fill(HIST("h_d0candidate_phi"), d0candidate.phi(), weight);

        registry.fill(HIST("h_jet_d0Prong0frag"), D0Prong0Frag, weight);
        registry.fill(HIST("h_jet_d0Prong1frag"), D0Prong1Frag, weight);
        registry.fill(HIST("h_d0candidate_prong0_pt"), d0candidate.ptProng0(), weight);
        registry.fill(HIST("h_d0candidate_prong1_pt"), d0candidate.ptProng1(), weight);
        registry.fill(HIST("h_decay_length"), d0candidate.decayLength(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_xy"), d0candidate.decayLengthXY(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_error"), d0candidate.errorDecayLength(), d0candidate.pt(), weight);
        registry.fill(HIST("h_decay_length_xy_error"), d0candidate.errorDecayLengthXY(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong0"), d0candidate.impactParameter0(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong1"), d0candidate.impactParameter1(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong0_error"), d0candidate.errorImpactParameter0(), d0candidate.pt(), weight);
        registry.fill(HIST("h_impact_parameter_d0Prong1_error"), d0candidate.errorImpactParameter1(), d0candidate.pt(), weight);

        if (d0candidate.isSelD0() >= 1 || d0candidate.isSelD0bar() >= 1) {
          registry.fill(HIST("h_d0candidate_mass"), massD0);
          registry.fill(HIST("h2_d0candidate_mass_d0candidate_pt"), massD0, d0candidate.pt());
          registry.fill(HIST("h2_d0candidate_mass_d0candidate_phi"), massD0, d0candidate.phi());
        }
      }
    }
  }

  template <typename T>
  void fillMCPHistograms(T const& partjets, float weight = 1.0)
  {
    for (const auto& jet : partjets) {
      registry.fill(HIST("h_part_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_part_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_part_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_part_jet_ntracks"), jet.tracks().size(), weight);
      for (auto& track : jet.template tracks_as<JetParticles2Prong>()) {
        auto chargeFrag = track.pt() / jet.pt();
        registry.fill(HIST("h_part_jet_part_chargefrag"), chargeFrag, weight);
        registry.fill(HIST("h2_part_jet_pt_part_track_pt"), jet.pt(), track.pt(), weight);
        registry.fill(HIST("h2_part_jet_pt_part_track_eta"), jet.pt(), track.eta(), weight);
        registry.fill(HIST("h2_part_jet_pt_part_track_phi"), jet.pt(), track.phi(), weight);
        registry.fill(HIST("h_part_track_pt"), track.pt(), weight);
        registry.fill(HIST("h_part_track_eta"), track.eta(), weight);
        registry.fill(HIST("h_part_track_phi"), track.phi(), weight);
      }
      for (auto& d0candidate : jet.template hfcandidates_as<JetParticles2Prong>()) {
        auto D0Frag = d0candidate.pt() / jet.pt();
        registry.fill(HIST("h_part_jet_part_d0frag"), D0Frag, weight);
        registry.fill(HIST("h2_part_jet_pt_part_d0candidate_pt"), jet.pt(), d0candidate.pt(), weight);
        registry.fill(HIST("h2_part_jet_pt_part_d0candidate_eta"), jet.pt(), d0candidate.eta(), weight);
        registry.fill(HIST("h2_part_jet_pt_part_d0candidate_phi"), jet.pt(), d0candidate.phi(), weight);
        registry.fill(HIST("h_part_d0candidate_pt"), d0candidate.pt(), weight);
        registry.fill(HIST("h_part_d0candidate_eta"), d0candidate.eta(), weight);
        registry.fill(HIST("h_part_d0candidate_phi"), d0candidate.phi(), weight);

        int counter = 0;
        std::array<float, 2> ptProngs;
        std::array<float, 2> etaProngs;
        std::array<float, 2> yProngs;
        for (auto const& daught : d0candidate.template daughters_as<JetParticles2Prong>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(array{daught.px(), daught.py(), daught.pz()}, RecoDecay::getMassPDG(daught.pdgCode()));
          counter++;
        }

        if (isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1])) {
          registry.fill(HIST("hPtGenWithProngsInAcceptance"), d0candidate.pt());
          registry.fill(HIST("hYGenWithProngsInAcceptance"), d0candidate.y(), d0candidate.pt());
          registry.fill(HIST("hEtaGenWithProngsInAcceptance"), d0candidate.eta(), d0candidate.pt());
        }
      }
    }
  }

  template <typename T>
  void fillMCMatchedHistograms(T const& mcdjet, float weight = 1.0)
  {
    if (mcdjet.has_matchedJetCand() && mcdjet.matchedJetCandId() >= 0) {
      const auto& mcpjet = mcdjet.template matchedJetCand_as<soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>>();
      registry.fill(HIST("h2_jet_pt_part_jet_pt"), mcpjet.pt(), mcdjet.pt(), weight);
    }
  }

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(JetFragmentationHFTask, processDummy, "Dummy process function turned on by default", true);

  void processJetsData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                       soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& hfjets,
                       CandidateD0Data const& candidates,
                       JetTracks const& tracks)
  {
    fillDataHistograms(hfjets);
  }
  PROCESS_SWITCH(JetFragmentationHFTask, processJetsData, "Task of jet fragmentation for heavy flavor (Data)", false);

  void processJetsMCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                      soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents> const& hfjets,
                      CandidateD0MC const& candidates,
                      JetTracks const& tracks)
  {
    fillMCDHistograms(hfjets);
  }
  PROCESS_SWITCH(JetFragmentationHFTask, processJetsMCD, "Task of jet fragmentation for heavy flavor (MCD)", false);

  void processJetsMCP(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                      soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents> const& hfjets,
                      JetParticles2Prong const& particles,
                      JetTracks const& tracks)
  {
    fillMCPHistograms(hfjets);
  }
  PROCESS_SWITCH(JetFragmentationHFTask, processJetsMCP, "Task of jet fragmentaion for heavy flavor (MCP)", false);

  void processJetsMCPMCDMatched(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                                soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets> const& mcdjets,
                                soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets> const& mcpjets,
                                JetParticles2Prong const& particles,
                                JetTracks const& tracks)
  {
    for (const auto& mcdjet : mcdjets) {
      fillMCMatchedHistograms(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFragmentationHFTask, processJetsMCPMCDMatched, "Matching of detector level jets and particle level jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetFragmentationHFTask>(cfgc, TaskName{"jet-fragmentation-hf"})}; }
