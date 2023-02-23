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

struct JetMatchingHFQA {
  using DetectorLevelJets = soa::Join<aod::MCDetectorLevelHFJets, aod::MCDetectorLevelHFJetConstituents, aod::MatchedMCDetectorParticleLevelHFJets>;
  using ParticleLevelJets = soa::Join<aod::MCParticleLevelHFJets, aod::MCParticleLevelHFJetConstituents, aod::MatchedMCParticleDetectorLevelHFJets>;

  OutputObj<TH2F> hJetPt{"h_jet_pt"};
  OutputObj<TH2F> hJetDetaDphi{"h_jet_deta_dphi"};
  OutputObj<TH1F> hJetDetPt{"h_jet_det_pt"};
  OutputObj<TH1F> hJetGenPt{"h_jet_gen_pt"};
  OutputObj<TH1F> hJetDetPhi{"h_jet_det_phi"};
  OutputObj<TH1F> hJetGenPhi{"h_jet_gen_phi"};
  OutputObj<TH1F> hJetDetEta{"h_jet_det_eta"};
  OutputObj<TH1F> hJetGenEta{"h_jet_gen_eta"};
  OutputObj<TH1F> hJetDetNTracks{"h_jet_det_ntracks"};
  OutputObj<TH1F> hJetGenNTracks{"h_jet_gen_ntracks"};

  void init(InitContext const&)
  {
    hJetPt.setObject(new TH2F("h_jet_pt", "HF-matched jets;jet p_{T}^{gen} (GeV/#it{c});jet p_{T}^{det} (GeV/#it{c})",
                              100, 0., 100., 100, 0., 100.));
    hJetDetaDphi.setObject(new TH2F("h_jet_deta_dphi", "HF-matched jets;jet #Delta#phi;#Delta#eta",
                                    100, -2. * TMath::Pi(), 2. * TMath::Pi(), 100, -2., 2.));
    hJetDetPt.setObject(new TH1F("h_jet_det_pt", "HF-matched jets;jet p_{T}^{det} (GeV/#it{c})", 100, 0., 100.));
    hJetGenPt.setObject(new TH1F("h_jet_gen_pt", "HF-matched jets;jet p_{T}^{gen} (GeV/#it{c})", 100, 0., 100.));
    hJetDetPhi.setObject(new TH1F("h_jet_det_phi", "jet #phi; #phi", 140, -7.0, 7.0));
    hJetGenPhi.setObject(new TH1F("h_jet_gen_phi", "jet #phi; #phi", 140, -7.0, 7.0));
    hJetDetEta.setObject(new TH1F("h_jet_det_eta", "jet #eta; #eta", 30, -1.5, 1.5));
    hJetGenEta.setObject(new TH1F("h_jet_gen_eta", "jet #eta; #eta", 30, -1.5, 1.5));
    hJetDetNTracks.setObject(new TH1F("h_jet_det_ntracks", "jet N tracks ; N tracks", 150, -0.5, 99.5));
    hJetGenNTracks.setObject(new TH1F("h_jet_gen_ntracks", "jet N tracks ; N tracks", 150, -0.5, 99.5));
  }

  void process(aod::Collisions::iterator const& collision,
               DetectorLevelJets const& djets, ParticleLevelJets const& pjets)
  {
    for (const auto& djet : djets) {
      if (djet.has_matchedJet() && djet.matchedJetId() >= 0) {
        const auto& pjet = djet.matchedJet_as<ParticleLevelJets>();
        LOGF(info, "djet %d (pt of %g GeV/c) is matched to %d (pt of %g GeV/c)",
             djet.globalIndex(), djet.pt(), djet.matchedJetId(), pjet.pt());
        hJetPt->Fill(pjet.pt(), djet.pt());
        hJetDetPt->Fill(djet.pt());
        hJetDetPhi->Fill(djet.phi());
        hJetDetEta->Fill(djet.eta());
        hJetDetNTracks->Fill(djet.tracksIds().size() + 1); // adding HF candidate
        const auto dphi = -TMath::Pi() + fmod(2 * TMath::Pi() + fmod(djet.phi() - pjet.phi() + TMath::Pi(), 2 * TMath::Pi()), 2 * TMath::Pi());
        hJetDetaDphi->Fill(dphi, djet.eta() - pjet.eta());
      }
    }
  }

  void processMCP(aod::McCollision const& collision,
                  ParticleLevelJets const& pjets, DetectorLevelJets const& djets)
  {
    LOGF(info, "analysing MC collision %d", collision.globalIndex());
    for (const auto& pjet : pjets) {
      if (pjet.has_matchedJet() && pjet.matchedJetId() >= 0) {
        const auto& djet = pjet.matchedJet_as<DetectorLevelJets>();
        LOGF(info, "pjet %d (pt of %g GeV/c) is matched to %d (pt of %g GeV/c)",
             pjet.globalIndex(), pjet.pt(), pjet.matchedJetId(), djet.pt());
        hJetGenPt->Fill(pjet.pt());
        hJetGenPhi->Fill(pjet.phi());
        hJetGenEta->Fill(pjet.eta());
        hJetGenNTracks->Fill(pjet.tracksIds().size() + 1); // adding HF candidate
      }
    }
  }
  PROCESS_SWITCH(JetMatchingHFQA, processMCP, "QA on generator-level jets", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetMatchingHFQA>(cfgc, TaskName{"jet-matching-hf-qa"})};
}
