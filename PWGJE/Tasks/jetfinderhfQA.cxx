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

// jet finder hf QA task
//
// Authors: Nima Zardoshti

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;

#include "Framework/runDataProcessing.h"

struct JetFinderHFQATask {
  OutputObj<TH2F> h2JetPt{"h2_jet_pt"};
  OutputObj<TH2F> h2JetPhi{"h2_jet_phi"};
  OutputObj<TH2F> h2JetEta{"h2_jet_eta"};
  OutputObj<TH2F> h2JetNTracks{"h2_jet_ntracks"};
  OutputObj<TH1F> hJetPt{"h_jet_pt"};
  OutputObj<TH1F> hJetPhi{"h_jet_phi"};
  OutputObj<TH1F> hJetEta{"h_jet_eta"};
  OutputObj<TH1F> hJetNTracks{"h_jet_ntracks"};
  OutputObj<TH1F> hCandPt{"h_cand_pt"};

  OutputObj<TH2F> h2JetPt_Part{"h2_jet_pt_part"};
  OutputObj<TH2F> h2JetPhi_Part{"h2_jet_phi_part"};
  OutputObj<TH2F> h2JetEta_Part{"h2_jet_eta_part"};
  OutputObj<TH2F> h2JetNTracks_Part{"h2_jet_ntracks_part"};
  OutputObj<TH1F> hJetPt_Part{"h_jet_pt_part"};
  OutputObj<TH1F> hJetPhi_Part{"h_jet_phi_part"};
  OutputObj<TH1F> hJetEta_Part{"h_jet_eta_part"};
  OutputObj<TH1F> hJetNTracks_Part{"h_jet_ntracks_part"};
  OutputObj<TH1F> hCandPt_Part{"h_cand_pt_part"};

  void init(o2::framework::InitContext&)
  {
    h2JetPt.setObject(new TH2F("h2_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})", 100, 0., 100., 10, 0.05, 1.05));
    h2JetPhi.setObject(new TH2F("h2_jet_phi", "jet #phi;#phi", 80, -1., 7., 10, 0.05, 1.05));
    h2JetEta.setObject(new TH2F("h2_jet_eta", "jet #eta;#eta", 70, -0.7, 0.7, 10, 0.05, 1.05));
    h2JetNTracks.setObject(new TH2F("h2_jet_ntracks", "jet n;n constituents", 30, 0., 30., 10, 0.05, 1.05));

    hJetPt.setObject(new TH1F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})", 100, 0., 100.));
    hJetPhi.setObject(new TH1F("h_jet_phi", "jet #phi; #phi", 140, -7.0, 7.0));
    hJetEta.setObject(new TH1F("h_jet_eta", "jet #eta; #eta", 30, -1.5, 1.5));
    hJetNTracks.setObject(new TH1F("h_jet_ntracks", "jet N tracks ; N tracks", 150, -0.5, 99.5));
    hCandPt.setObject(new TH1F("h_cand_pt", "jet p_{T,cand};p_{T,cand} (GeV/#it{c})", 100, 0., 100.));

    h2JetPt_Part.setObject(new TH2F("h2_jet_pt_part", "jet p_{T};p_{T} (GeV/#it{c})", 100, 0., 100., 10, 0.05, 1.05));
    h2JetPhi_Part.setObject(new TH2F("h2_jet_phi_part", "jet #phi;#phi", 80, -1., 7., 10, 0.05, 1.05));
    h2JetEta_Part.setObject(new TH2F("h2_jet_eta_part", "jet #eta;#eta", 70, -0.7, 0.7, 10, 0.05, 1.05));
    h2JetNTracks_Part.setObject(new TH2F("h2_jet_ntracks_part", "jet n;n constituents", 30, 0., 30., 10, 0.05, 1.05));

    hJetPt_Part.setObject(new TH1F("h_jet_pt_part", "jet p_{T};p_{T} (GeV/#it{c})", 100, 0., 100.));
    hJetPhi_Part.setObject(new TH1F("h_jet_phi_part", "jet #phi; #phi", 140, -7.0, 7.0));
    hJetEta_Part.setObject(new TH1F("h_jet_eta_part", "jet #eta; #eta", 30, -1.5, 1.5));
    hJetNTracks_Part.setObject(new TH1F("h_jet_ntracks_part", "jet N tracks ; N tracks", 150, -0.5, 99.5));
    hCandPt_Part.setObject(new TH1F("h_cand_pt_part", "jet p_{T,cand};p_{T,cand} (GeV/#it{c})", 100, 0., 100.));
  }

  using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using CandidateD0Data = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using CandidateD0MC = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>;
  using JetParticles2Prong = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
  using JetParticles3Prong = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  using JetParticlesBPlus = soa::Join<aod::McParticles, aod::HfCandBplusMcGen>;

  void processDummy(aod::Tracks const& track) {}
  PROCESS_SWITCH(JetFinderHFQATask, processDummy, "Dummy process function turned on by default", true);

  void processJetsData(soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>::iterator const& jet, CandidateD0Data const& candidates, JetTracks const& tracks)
  {
    h2JetPt->Fill(jet.pt(), jet.r() / 100.0);
    h2JetPhi->Fill(jet.phi(), jet.r() / 100.0);
    h2JetEta->Fill(jet.eta(), jet.r() / 100.0);
    h2JetNTracks->Fill(jet.tracks().size() + jet.hfcandidates().size(), jet.r() / 100.0);
    hJetPt->Fill(jet.pt());
    hJetPhi->Fill(jet.phi());
    hJetEta->Fill(jet.eta());
    hJetNTracks->Fill(jet.tracks().size() + jet.hfcandidates().size());
    for (auto& candidate : jet.hfcandidates_as<CandidateD0Data>()) {
      hCandPt->Fill(candidate.pt());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsData, "jet finder HF QA data", false);

  void processJetsMCD(soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>::iterator const& jet, CandidateD0MC const& candidates, JetTracks const& tracks)
  {
    h2JetPt->Fill(jet.pt(), jet.r() / 100.0);
    h2JetPhi->Fill(jet.phi(), jet.r() / 100.0);
    h2JetEta->Fill(jet.eta(), jet.r() / 100.0);
    h2JetNTracks->Fill(jet.tracks().size() + jet.hfcandidates().size(), jet.r() / 100.0);
    hJetPt->Fill(jet.pt());
    hJetPhi->Fill(jet.phi());
    hJetEta->Fill(jet.eta());
    hJetNTracks->Fill(jet.tracks().size() + jet.hfcandidates().size());
    for (auto& candidate : jet.hfcandidates_as<CandidateD0MC>()) {
      hCandPt->Fill(candidate.pt());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCD, "jet finder HF QA mcd", false);

  void processJetsMCP(soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>::iterator const& jet, JetParticles2Prong const& particles)
  {
    h2JetPt_Part->Fill(jet.pt(), jet.r() / 100.0);
    h2JetPhi_Part->Fill(jet.phi(), jet.r() / 100.0);
    h2JetEta_Part->Fill(jet.eta(), jet.r() / 100.0);
    h2JetNTracks_Part->Fill(jet.tracks().size() + jet.hfcandidates().size(), jet.r() / 100.0);
    hJetPt_Part->Fill(jet.pt());
    hJetPhi_Part->Fill(jet.phi());
    hJetEta_Part->Fill(jet.eta());
    hJetNTracks_Part->Fill(jet.tracks().size() + jet.hfcandidates().size());
    for (auto& candidate : jet.hfcandidates_as<JetParticles2Prong>()) {
      hCandPt_Part->Fill(candidate.pt());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCP, "jet finder HF QA mcp", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetFinderHFQATask>(cfgc, TaskName{"jet-finder-hf-qa"})}; }
