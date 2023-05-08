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

template <typename BaseJetCollection, typename TagJetCollection>
struct JetMatchingQA {
  OutputObj<TH2F> hJetPt{"h_jet_pt"};
  OutputObj<TH2F> hJetDetaDphi{"h_jet_deta_dphi"};
  OutputObj<TH2F> hJetGeoPt{"h_jet_geo_pt"};
  OutputObj<TH2F> hJetGeoDetaDphi{"h_jet_geo_deta_dphi"};
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
    hJetGeoPt.setObject(new TH2F("h_jet_geo_pt", "geo-matched jets;jet p_{T}^{gen} (GeV/#it{c});jet p_{T}^{det} (GeV/#it{c})",
                                 100, 0., 100., 100, 0., 100.));
    hJetGeoDetaDphi.setObject(new TH2F("h_jet_geo_deta_dphi", "geo-matched jets;jet #Delta#phi;#Delta#eta",
                                       100, -2. * TMath::Pi(), 2. * TMath::Pi(), 100, -2., 2.));
    hJetPt.setObject(new TH2F("h_jet_pt", "HF-matched jets;jet p_{T}^{gen} (GeV/#it{c});jet p_{T}^{det} (GeV/#it{c})",
                              100, 0., 100., 100, 0., 100.));
    hJetDetaDphi.setObject(new TH2F("h_jet_deta_dphi", "HF-matched jets;jet #Delta#phi;#Delta#eta",
                                    100, -2. * TMath::Pi(), 2. * TMath::Pi(), 100, -2., 2.));

    hJetDetPt.setObject(new TH1F("h_jet_det_pt", "detector level jets;jet p_{T}^{det} (GeV/#it{c})", 100, 0., 100.));
    hJetGenPt.setObject(new TH1F("h_jet_gen_pt", "particle level jets;jet p_{T}^{gen} (GeV/#it{c})", 100, 0., 100.));
    hJetDetPhi.setObject(new TH1F("h_jet_det_phi", "jet #phi; #phi", 140, -7.0, 7.0));
    hJetGenPhi.setObject(new TH1F("h_jet_gen_phi", "jet #phi; #phi", 140, -7.0, 7.0));
    hJetDetEta.setObject(new TH1F("h_jet_det_eta", "jet #eta; #eta", 30, -1.5, 1.5));
    hJetGenEta.setObject(new TH1F("h_jet_gen_eta", "jet #eta; #eta", 30, -1.5, 1.5));
    hJetDetNTracks.setObject(new TH1F("h_jet_det_ntracks", "jet N tracks ; N tracks", 150, -0.5, 99.5));
    hJetGenNTracks.setObject(new TH1F("h_jet_gen_ntracks", "jet N tracks ; N tracks", 150, -0.5, 99.5));
  }

  void processMCD(aod::Collisions::iterator const& collision,
                  BaseJetCollection const& djets, TagJetCollection const& pjets)
  {
    for (const auto& djet : djets) {
      if (djet.has_matchedJetCand() || djet.has_matchedJetGeo()) {
        hJetDetPhi->Fill(djet.phi());
        hJetDetEta->Fill(djet.eta());
        hJetDetNTracks->Fill(djet.tracksIds().size() + 1); // adding HF candidate
      }

      if (djet.has_matchedJetCand() && djet.matchedJetCandId() >= 0) {
        const auto& pjet = djet.template matchedJetCand_as<TagJetCollection>();
        LOGF(info, "djet %d (pt of %g GeV/c) is HF-matched to %d (pt of %g GeV/c)",
             djet.globalIndex(), djet.pt(), djet.matchedJetCandId(), pjet.pt());
        hJetPt->Fill(pjet.pt(), djet.pt());
        const auto dphi = -TMath::Pi() + fmod(2 * TMath::Pi() + fmod(djet.phi() - pjet.phi() + TMath::Pi(), 2 * TMath::Pi()), 2 * TMath::Pi());
        hJetDetaDphi->Fill(dphi, djet.eta() - pjet.eta());
      }

      if (djet.has_matchedJetGeo()) {
        const auto& pjet = djet.template matchedJetGeo_as<TagJetCollection>();
        LOGF(info, "djet %d (pt of %g GeV/c) is geo-matched to %d (pt of %g GeV/c)",
             djet.globalIndex(), djet.pt(), djet.matchedJetGeoId(), pjet.pt());
        hJetGeoPt->Fill(pjet.pt(), djet.pt());
        const auto dphi = -TMath::Pi() + fmod(2 * TMath::Pi() + fmod(djet.phi() - pjet.phi() + TMath::Pi(), 2 * TMath::Pi()), 2 * TMath::Pi());
        hJetGeoDetaDphi->Fill(dphi, djet.eta() - pjet.eta());
      }
    }
  }
  PROCESS_SWITCH(JetMatchingQA, processMCD, "QA on detector-level jets", true);

  void processMCP(aod::McCollision const& collision,
                  TagJetCollection const& pjets, BaseJetCollection const& djets)
  {
    for (const auto& pjet : pjets) {
      if (pjet.has_matchedJetCand() || pjet.has_matchedJetGeo()) {
        hJetGenPt->Fill(pjet.pt());
        hJetGenPhi->Fill(pjet.phi());
        hJetGenEta->Fill(pjet.eta());
        hJetGenNTracks->Fill(pjet.tracksIds().size() + 1); // adding HF candidate
      }

      if (pjet.has_matchedJetCand() && pjet.matchedJetCandId() >= 0) {
        const auto& djet = pjet.template matchedJetCand_as<BaseJetCollection>();
        LOGF(info, "pjet %d (pt of %g GeV/c) is HF-matched to %d (pt of %g GeV/c)",
             pjet.globalIndex(), pjet.pt(), pjet.matchedJetCandId(), djet.pt());
      }

      if (pjet.has_matchedJetGeo()) {
        const auto& djet = pjet.template matchedJetGeo_as<BaseJetCollection>();
        LOGF(info, "pjet %d (pt of %g GeV/c) is geo-matched to %d (pt of %g GeV/c)",
             pjet.globalIndex(), pjet.pt(), pjet.matchedJetGeoId(), djet.pt());
      }
    }
  }
  PROCESS_SWITCH(JetMatchingQA, processMCP, "QA on generator-level jets", true);
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

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<ChargedJetMatchingQA>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-qa-ch"}));
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetMatchingQA>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-qa-d0-ch"}));
  tasks.emplace_back(adaptAnalysisTask<LcChargedJetMatchingQA>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-qa-lc-ch"}));
  tasks.emplace_back(adaptAnalysisTask<BplusChargedJetMatchingQA>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-qa-bplus-ch"}));

  return WorkflowSpec{tasks};
}
