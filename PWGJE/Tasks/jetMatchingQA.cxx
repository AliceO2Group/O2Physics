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

/// \file jetmatchinghfqa.cxx
/// \brief Basic QA of HF jet matching
///
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Aimeric Lanodu <aimeric.landou@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename BaseJetCollection, typename TagJetCollection>
struct JetMatchingQA {

  HistogramRegistry registry{
    "registry",
    {
      {"h_jet_match_hf_pt", "hf-matched jets", {HistType::kTH2F, {{1000, 0.0f, 100.0f, "#it{p}_{T}^{gen} (GeV/#it{c})"}, {1000, 0.0f, 100.0f, "#it{p}_{T}^{det} (GeV/#it{c})"}}}},
      {"h_jet_match_hf_deta_dphi", "hf-matched jets", {HistType::kTH2F, {{100, -2. * TMath::Pi(), 2. * TMath::Pi(), "jet #Delta#phi"}, {100, -2.0f, 2.0f, "#Delta#eta"}}}},
      {"h_jet_match_geo_pt", "geo-matched jets", {HistType::kTH2F, {{1000, 0.0f, 100.0f, "#it{p}_{T}^{gen} (GeV/#it{c})"}, {1000, 0.0f, 100.0f, "#it{p}_{T}^{det} (GeV/#it{c})"}}}},
      {"h_jet_match_geo_deta_dphi", "geo-matched jets", {HistType::kTH2F, {{100, -2. * TMath::Pi(), 2. * TMath::Pi(), "jet #Delta#phi"}, {100, -2.0f, 2.0f, "#Delta#eta"}}}},
      {"h_jet_match_pt_pt", "pt-matched jets", {HistType::kTH2F, {{1000, 0.0f, 100.0f, "#it{p}_{T}^{gen} (GeV/#it{c})"}, {1000, 0.0f, 100.0f, "#it{p}_{T}^{det} (GeV/#it{c})"}}}},
      {"h_jet_match_pt_deta_dphi", "pt-matched jets", {HistType::kTH2F, {{100, -2. * TMath::Pi(), 2. * TMath::Pi(), "jet #Delta#phi"}, {100, -2.0f, 2.0f, "#Delta#eta"}}}},

      {"h_jet_match_geo_pt_zoom", "geo-matched jets", {HistType::kTH2F, {{1000, 0.0f, 10.0f, "#it{p}_{T} (particle level, GeV/#it{c})"}, {1000, 0.0f, 10.0f, "#it{p}_{T} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_geo_dpt", "geo-matched jets", {HistType::kTH1F, {{100, -10.0f, 10.0f, "#Delta#it{p}_{T} (particle level vs detector level, GeV/#it{c})"}}}},
      {"h2_jet_pt_jet_match_geo_dptoverpt", "geo-matched jets", {HistType::kTH2F, {{2000, 0.0f, 200.0f, "#it{p}_{T} (detector level, GeV/#it{c})"}, {700, -5.0f, 2.0f, "(#it{p}_{T, part}-#it{p}_{T, det})/#it{p}_{T,part} (particle level vs detector level, GeV/#it{c})"}}}},
      {"h_jet_match_geo_PtLeadingPart", "geo-matched jets", {HistType::kTH2F, {{1000, 0.0f, 10.0f, "#it{p}_{T}^{leading} (particle level, GeV/#it{c})"}, {1000, 0.0f, 10.0f, "#it{p}_{T}^{leading} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_geo_phi", "geo-matched jets", {HistType::kTH2F, {{80, -1.0f, 7.0f, "#phi_{jet} (particle level, rad)"}, {80, -1.0f, 7.0f, "#phi_{jet} (detector level, rad)"}}}},
      {"h_jet_match_geo_eta", "geo-matched jets", {HistType::kTH2F, {{70, -0.7f, 0.7f, "#it{p}_{T}^{particle level} (GeV/#it{c})"}, {70, -0.7f, 0.7f, "#it{p}_{T} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_geo_Nconst", "geo-matched jets", {HistType::kTH2F, {{30, 0.0f, 30.0f, "n constituents (particle level)"}, {30, 0.0f, 30.0f, "n constituents (detector level)"}}}},

      {"h_jet_match_hf_pt_zoom", "hf-matched jets", {HistType::kTH2F, {{1000, 0.0f, 10.0f, "#it{p}_{T} (particle level, GeV/#it{c})"}, {1000, 0.0f, 10.0f, "#it{p}_{T} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_hf_dpt", "hf-matched jets", {HistType::kTH1F, {{100, -10.0f, 10.0f, "#Delta#it{p}_{T} (particle level vs detector level, GeV/#it{c})"}}}},
      {"h_jet_match_hf_PtLeadingPart", "hf-matched jets", {HistType::kTH2F, {{1000, 0.0f, 10.0f, "#it{p}_{T}^{leading} (particle level, GeV/#it{c})"}, {1000, 0.0f, 10.0f, "#it{p}_{T}^{leading} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_hf_phi", "hf-matched jets", {HistType::kTH2F, {{80, -1.0f, 7.0f, "#phi_{jet} (particle level, rad)"}, {80, -1.0f, 7.0f, "#phi_{jet} (detector level, rad)"}}}},
      {"h_jet_match_hf_eta", "hf-matched jets", {HistType::kTH2F, {{70, -0.7f, 0.7f, "#it{p}_{T}^{particle level} (GeV/#it{c})"}, {70, -0.7f, 0.7f, "#it{p}_{T} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_hf_Nconst", "hf-matched jets", {HistType::kTH2F, {{30, 0.0f, 30.0f, "n constituents (particle level)"}, {30, 0.0f, 30.0f, "n constituents (detector level)"}}}},

      {"h_jet_match_pt_pt_zoom", "pt-matched jets", {HistType::kTH2F, {{1000, 0.0f, 10.0f, "#it{p}_{T} (particle level, GeV/#it{c})"}, {1000, 0.0f, 10.0f, "#it{p}_{T} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_pt_dpt", "pt-matched jets", {HistType::kTH1F, {{100, -10.0f, 10.0f, "#Delta#it{p}_{T} (particle level vs detector level, GeV/#it{c})"}}}},
      {"h_jet_match_pt_PtLeadingPart", "pt-matched jets", {HistType::kTH2F, {{1000, 0.0f, 10.0f, "#it{p}_{T}^{leading} (particle level, GeV/#it{c})"}, {1000, 0.0f, 10.0f, "#it{p}_{T}^{leading} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_pt_phi", "pt-matched jets", {HistType::kTH2F, {{80, -1.0f, 7.0f, "#phi_{jet} (particle level, rad)"}, {80, -1.0f, 7.0f, "#phi_{jet} (detector level, rad)"}}}},
      {"h_jet_match_pt_eta", "pt-matched jets", {HistType::kTH2F, {{70, -0.7f, 0.7f, "#it{p}_{T}^{particle level} (GeV/#it{c})"}, {70, -0.7f, 0.7f, "#it{p}_{T} (detector level, GeV/#it{c})"}}}},
      {"h_jet_match_pt_Nconst", "pt-matched jets", {HistType::kTH2F, {{30, 0.0f, 30.0f, "n constituents (particle level)"}, {30, 0.0f, 30.0f, "n constituents (detector level)"}}}},

      {"h_jet_det_pt", "detector level jets", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "jet p_{T}^{det} (GeV/#it{c})"}}}},
      {"h_jet_gen_pt", "particle level jets", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "jet p_{T}^{gen} (GeV/#it{c})"}}}},
      {"h_jet_det_phi", "detector level jet #phi", {HistType::kTH1F, {{140, -7.0f, 7.0f, "#phi"}}}},
      {"h_jet_gen_phi", "particle level jets #phi", {HistType::kTH1F, {{140, -7.0f, 7.0f, "#phi"}}}},
      {"h_jet_det_eta", "detector level jets #eta", {HistType::kTH1F, {{30, -1.5f, 1.5f, "#eta"}}}},
      {"h_jet_gen_eta", "particle level jets #eta", {HistType::kTH1F, {{30, -1.5f, 1.5f, "#eta"}}}},
      {"h_jet_det_ntracks", "detector level jets N tracks", {HistType::kTH1F, {{150, -0.5f, 99.5f, "N tracks"}}}},
      {"h_jet_gen_ntracks", "particle level jets N tracks", {HistType::kTH1F, {{150, -0.5f, 99.5f, "N tracks"}}}},
    },
  };

  void init(InitContext const&)
  {
  }

  void processDummy(aod::JetMcCollisions const&)
  {
  }
  PROCESS_SWITCH(JetMatchingQA, processDummy, "Dummy process", true);

  void processMCD(aod::JetCollision const&, aod::JetParticles const&, aod::JetTracksMCD const&,
                  BaseJetCollection const& djets, TagJetCollection const&)
  {
    for (const auto& djet : djets) {
      if (djet.has_matchedJetCand() || djet.has_matchedJetGeo()) {
        registry.fill(HIST("h_jet_det_pt"), djet.pt());
        registry.fill(HIST("h_jet_det_phi"), djet.phi());
        registry.fill(HIST("h_jet_det_eta"), djet.eta());
        registry.fill(HIST("h_jet_det_ntracks"), djet.tracksIds().size() + 1); // adding HF candidate
      }

      // HF matching QA
      for (auto& pjet : djet.template matchedJetCand_as<TagJetCollection>()) {
        registry.fill(HIST("h_jet_match_hf_pt"), pjet.pt(), djet.pt());
        const auto dphi = -TMath::Pi() + fmod(2 * TMath::Pi() + fmod(djet.phi() - pjet.phi() + TMath::Pi(), 2 * TMath::Pi()), 2 * TMath::Pi());
        registry.fill(HIST("h_jet_match_hf_deta_dphi"), dphi, djet.eta() - pjet.eta());

        registry.fill(HIST("h_jet_match_hf_pt_zoom"), pjet.pt(), djet.pt());
        registry.fill(HIST("h_jet_match_hf_dpt"), pjet.pt() - djet.pt());
        registry.fill(HIST("h_jet_match_hf_phi"), pjet.phi(), djet.phi());
        registry.fill(HIST("h_jet_match_hf_eta"), pjet.eta(), djet.eta());
        registry.fill(HIST("h_jet_match_hf_Nconst"), pjet.tracksIds().size(), djet.tracksIds().size());

        double pjet_pt_lead = 0.;
        for (auto& mcparticle : pjet.template tracks_as<aod::JetParticles>()) {
          if (mcparticle.pt() > pjet_pt_lead) {
            pjet_pt_lead = mcparticle.pt();
          }
        }
        double djet_pt_lead = 0.;
        for (auto& track : djet.template tracks_as<aod::JetTracksMCD>()) {
          if (track.pt() > djet_pt_lead) {
            djet_pt_lead = track.pt();
          }
        }
        registry.fill(HIST("h_jet_match_hf_PtLeadingPart"), djet_pt_lead, pjet_pt_lead);
      }

      // geo matching QA
      for (auto& pjet : djet.template matchedJetGeo_as<TagJetCollection>()) {
        registry.fill(HIST("h_jet_match_geo_pt"), pjet.pt(), djet.pt());
        const auto dphi = -TMath::Pi() + fmod(2 * TMath::Pi() + fmod(djet.phi() - pjet.phi() + TMath::Pi(), 2 * TMath::Pi()), 2 * TMath::Pi());
        registry.fill(HIST("h_jet_match_geo_deta_dphi"), dphi, djet.eta() - pjet.eta());

        registry.fill(HIST("h_jet_match_geo_pt_zoom"), pjet.pt(), djet.pt());
        registry.fill(HIST("h_jet_match_geo_dpt"), pjet.pt() - djet.pt());
        registry.fill(HIST("h2_jet_pt_jet_match_geo_dptoverpt"), pjet.pt(), (pjet.pt() - djet.pt()) * 1. / pjet.pt());
        registry.fill(HIST("h_jet_match_geo_phi"), pjet.phi(), djet.phi());
        registry.fill(HIST("h_jet_match_geo_eta"), pjet.eta(), djet.eta());
        registry.fill(HIST("h_jet_match_geo_Nconst"), pjet.tracksIds().size(), djet.tracksIds().size());

        double pjet_pt_lead = 0.;
        for (auto& mcparticle : pjet.template tracks_as<aod::JetParticles>()) {
          if (mcparticle.pt() > pjet_pt_lead) {
            pjet_pt_lead = mcparticle.pt();
          }
        }
        double djet_pt_lead = 0.;
        for (auto& track : djet.template tracks_as<aod::JetTracksMCD>()) {
          if (track.pt() > djet_pt_lead) {
            djet_pt_lead = track.pt();
          }
        }
        registry.fill(HIST("h_jet_match_geo_PtLeadingPart"), djet_pt_lead, pjet_pt_lead);
      }

      // pT matching QA
      for (auto& pjet : djet.template matchedJetPt_as<TagJetCollection>()) {
        registry.fill(HIST("h_jet_match_pt_pt"), pjet.pt(), djet.pt());
        const auto dphi = -TMath::Pi() + fmod(2 * TMath::Pi() + fmod(djet.phi() - pjet.phi() + TMath::Pi(), 2 * TMath::Pi()), 2 * TMath::Pi());
        registry.fill(HIST("h_jet_match_pt_deta_dphi"), dphi, djet.eta() - pjet.eta());

        registry.fill(HIST("h_jet_match_pt_pt_zoom"), pjet.pt(), djet.pt());
        registry.fill(HIST("h_jet_match_pt_dpt"), pjet.pt() - djet.pt());
        registry.fill(HIST("h_jet_match_pt_phi"), pjet.phi(), djet.phi());
        registry.fill(HIST("h_jet_match_pt_eta"), pjet.eta(), djet.eta());
        registry.fill(HIST("h_jet_match_pt_Nconst"), pjet.tracksIds().size(), djet.tracksIds().size());

        double pjet_pt_lead = 0.;
        for (auto& mcparticle : pjet.template tracks_as<aod::JetParticles>()) {
          if (mcparticle.pt() > pjet_pt_lead) {
            pjet_pt_lead = mcparticle.pt();
          }
        }
        double djet_pt_lead = 0.;
        for (auto& track : djet.template tracks_as<aod::JetTracksMCD>()) {
          if (track.pt() > djet_pt_lead) {
            djet_pt_lead = track.pt();
          }
        }
        registry.fill(HIST("h_jet_match_pt_PtLeadingPart"), djet_pt_lead, pjet_pt_lead);
      }
    }
  }
  PROCESS_SWITCH(JetMatchingQA, processMCD, "QA on detector-level jets", false);

  void processMCP(aod::JetMcCollision const&,
                  TagJetCollection const& pjets, BaseJetCollection const&)
  {
    for (const auto& pjet : pjets) {
      if (pjet.has_matchedJetCand() || pjet.has_matchedJetGeo()) {
        registry.fill(HIST("h_jet_gen_pt"), pjet.pt());
        registry.fill(HIST("h_jet_gen_phi"), pjet.phi());
        registry.fill(HIST("h_jet_gen_eta"), pjet.eta());
        registry.fill(HIST("h_jet_gen_ntracks"), pjet.tracksIds().size() + 1); // adding HF candidate
      }
    }
  }
  PROCESS_SWITCH(JetMatchingQA, processMCP, "QA on generator-level jets", false);
};

using ChargedDetectorLevelJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;
using ChargedParticleLevelJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;
using ChargedJetMatchingQA = JetMatchingQA<ChargedDetectorLevelJets, ChargedParticleLevelJets>;

using D0ChargedDetectorLevelJets = soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;
using D0ChargedParticleLevelJets = soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;
using D0ChargedJetMatchingQA = JetMatchingQA<D0ChargedDetectorLevelJets, D0ChargedParticleLevelJets>;

using LcChargedDetectorLevelJets = soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>;
using LcChargedParticleLevelJets = soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>;
using LcChargedJetMatchingQA = JetMatchingQA<LcChargedDetectorLevelJets, LcChargedParticleLevelJets>;

using BplusChargedDetectorLevelJets = soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>;
using BplusChargedParticleLevelJets = soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>;
using BplusChargedJetMatchingQA = JetMatchingQA<BplusChargedDetectorLevelJets, BplusChargedParticleLevelJets>;

using DielectronChargedDetectorLevelJets = soa::Join<aod::DielectronChargedMCDetectorLevelJets, aod::DielectronChargedMCDetectorLevelJetConstituents, aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets>;
using DielectronChargedParticleLevelJets = soa::Join<aod::DielectronChargedMCParticleLevelJets, aod::DielectronChargedMCParticleLevelJetConstituents, aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets>;
using DielectronChargedJetMatchingQA = JetMatchingQA<DielectronChargedDetectorLevelJets, DielectronChargedParticleLevelJets>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<ChargedJetMatchingQA>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-qa-ch"}));
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetMatchingQA>(cfgc, TaskName{"jet-matching-qa-d0-ch"}));
  tasks.emplace_back(adaptAnalysisTask<LcChargedJetMatchingQA>(cfgc, TaskName{"jet-matching-qa-lc-ch"}));
  // tasks.emplace_back(adaptAnalysisTask<BplusChargedJetMatchingQA>(cfgc, TaskName{"jet-matching-qa-bplus-ch"}));
  tasks.emplace_back(adaptAnalysisTask<DielectronChargedJetMatchingQA>(cfgc, TaskName{"jet-matching-qa-dielectron-ch"}));

  return WorkflowSpec{tasks};
}
