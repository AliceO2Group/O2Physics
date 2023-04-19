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

// jet finder task
//
// Authors: Nima Zardoshti, Jochen Klein

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec hfjetMode = {
    "hfjetMode",
    VariantType::String,
    "",
    {"HF jet finder mode."},
  };
  workflowOptions.push_back(hfjetMode);
}

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename JetTable, typename TrackConstituentTable>
struct JetFinderHFTask {
  Produces<JetTable> jetsTable;
  Produces<TrackConstituentTable> trackConstituents;
  OutputObj<TH1F> hJetPt{"h_jet_pt"};
  OutputObj<TH1F> hJetPhi{"h_jet_phi"};
  OutputObj<TH1F> hJetEta{"h_jet_eta"};
  OutputObj<TH1F> hJetNTracks{"h_jet_ntracks"};
  OutputObj<TH1F> hD0Pt{"h_D0_pt"};

  Service<O2DatabasePDG> pdg;

  TrackSelection globalTracks;

  std::vector<fastjet::PseudoJet> jets;
  std::vector<fastjet::PseudoJet> inputParticles;
  JetFinder jetFinder;

  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.8, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.8, "maximum track eta"};
  Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};

  void init(InitContext const&)
  {
    // set up global tracks and adjust as necessary
    globalTracks = getGlobalTrackSelection();
    globalTracks.SetEtaRange(trackEtaMin, trackEtaMax);

    hJetPt.setObject(new TH1F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                              100, 0., 100.));
    hJetPhi.setObject(new TH1F("h_jet_phi", "jet #phi; #phi",
                               140, -7.0, 7.0));
    hJetEta.setObject(new TH1F("h_jet_eta", "jet #eta; #eta",
                               30, -1.5, 1.5));
    hJetNTracks.setObject(new TH1F("h_jet_ntracks", "jet N tracks ; N tracks",
                                   150, -0.5, 99.5));
    hD0Pt.setObject(new TH1F("h_D0_pt", "jet p_{T,D};p_{T,D} (GeV/#it{c})",
                             100, 0., 100.));

    jetFinder.etaMin = trackEtaMin;
    jetFinder.etaMax = trackEtaMax;
    jetFinder.jetPtMin = jetPtMin;
    jetFinder.jetPtMax = jetPtMax;
    jetFinder.jetR = jetR;
  }

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  //need enum as configurable
  enum pdgCode { pdgD0 = 421 };

  Filter trackCuts = (aod::track::pt >= trackPtMin && aod::track::pt < trackPtMax && aod::track::eta > trackEtaMin && aod::track::eta < trackEtaMax);
  Filter partCuts = (aod::mcparticle::pt >= trackPtMin && aod::mcparticle::pt < trackPtMax);
  Filter candCuts = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);

  using JetTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;

  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                   JetTracks const& tracks,
                   soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> const& candidates)
  {
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;

    if (!collision.sel8())
      return;

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      if (yD0(candidate) < candYMin || yD0(candidate) > candYMax) {
        continue;
      }
      if (candidate.pt() < candPtMin || candidate.pt() >= candPtMax) {
        continue;
      }
      jets.clear();
      inputParticles.clear();
      for (auto& track : tracks) {
        if (!globalTracks.IsSelected(track)) {
          continue;
        }
        if (candidate.prong0_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
        fillConstituents(track, inputParticles, track.globalIndex());
      }
      fillConstituents(candidate, inputParticles, -1, RecoDecay::getMassPDG(pdgD0));
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
      for (const auto& jet : jets) {
        isHFJet = false;
        if (jet.eta() < trackEtaMin + jetR || jet.eta() > trackEtaMax - jetR) {
          continue;
        }
        if (jet.perp() < jetPtMin || jet.perp() >= jetPtMax) {
          continue;
        }
        std::vector<int> trackconst;
        std::vector<int> candconst;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == -1) {
            isHFJet = true;
            break;
          }
        }

        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          const auto& constituents = sorted_by_pt(jet.constituents());
          for (const auto& constituent : constituents) {
            if (constituent.user_index() != -1) {
              // auto track = tracks.rawIteratorAt(constituent.user_index() - tracks.offset());
              trackconst.push_back(constituent.user_index());
            }
          }
          candconst.push_back(candidate.globalIndex()); // is this grouped per collision too?
          trackConstituents(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
          hJetPhi->Fill(jet.phi());
          hJetEta->Fill(jet.rap());
          hJetNTracks->Fill(jet.constituents().size());
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processData, "HF jet finding on data", true);

  void processMCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                  JetTracks const& tracks,
                  soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> const& candidates)
  {
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;
    //this loop should be made more efficient
    // TODO: should probably refine the candidate selection
    for (auto& candidate : candidates) {
      if (yD0(candidate) < candYMin || yD0(candidate) > candYMax) {
        continue;
      }
      if (candidate.pt() < candPtMin || candidate.pt() >= candPtMax) {
        continue;
      }
      // are the next two ifs needed?
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (!(std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      jets.clear();
      inputParticles.clear();
      for (auto& track : tracks) {
        if (!globalTracks.IsSelected(track)) {
          continue;
        }
        if (candidate.prong0_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
        fillConstituents(track, inputParticles, track.globalIndex());
      }
      fillConstituents(candidate, inputParticles, -1, RecoDecay::getMassPDG(pdgD0));
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
      for (const auto& jet : jets) {
        isHFJet = false;
        if (jet.eta() < trackEtaMin + jetR || jet.eta() > trackEtaMax - jetR) { // is this needed?
          continue;
        }
        if (jet.perp() < jetPtMin || jet.perp() >= jetPtMax) {
          continue;
        }
        std::vector<int> trackconst;
        std::vector<int> candconst;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == -1) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          const auto& constituents = sorted_by_pt(jet.constituents());
          for (const auto& constituent : constituents) {
            if (constituent.user_index() != -1) {
              // auto track = tracks.rawIteratorAt(constituent.user_index() - tracks.offset());
              trackconst.push_back(constituent.user_index());
            }
          }
          candconst.push_back(candidate.globalIndex()); // check if its correct
          trackConstituents(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
          hJetPhi->Fill(jet.phi());
          hJetEta->Fill(jet.rap());
          hJetNTracks->Fill(jet.constituents().size());
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processMCD, "HF jet finding on MC detector level", false);

  void processMCP(aod::McCollision const& collision,
                  soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>> const& particles)
  {
    LOG(debug) << "Per Event MCP";
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;

    // TODO: probably should do this as a filter
    std::vector<soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>::iterator> candidates;
    for (auto const& part : particles) {
      // TODO: generalise to any D0
      if (std::abs(part.flagMcMatchGen()) & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        auto partY = RecoDecay::y(array{part.px(), part.py(), part.pz()}, RecoDecay::getMassPDG(part.pdgCode()));
        if (partY < candYMin || partY > candYMax) {
          continue;
        }
        if (part.pt() < candPtMin || part.pt() >= candPtMax) {
          continue;
        }
        candidates.push_back(part);
      }
    }

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      jets.clear();
      inputParticles.clear();
      for (auto& track : particles) {
        // exclude neutral particles
        // TODO: can we do this through the filter?
        if (track.eta() < trackEtaMin || track.eta() > trackEtaMax) {
          continue;
        }
        auto p = pdg->GetParticle(track.pdgCode());
        //   track.globalIndex(), track.getGenStatusCode(), p ? std::abs(p->Charge()) : -999.);
        if ((track.getGenStatusCode() != 1) || (p ? std::abs(p->Charge()) : 0.) < 3.) {
          continue;
        }
        // TODO: check what mass to use?
        const auto daughters = candidate.daughtersIds();
        if (std::find(std::begin(daughters), std::end(daughters), track.globalIndex()) != std::end(daughters)) {
          continue;
        }
        fillConstituents(track, inputParticles, track.globalIndex(), RecoDecay::getMassPDG(track.pdgCode()));
      }
      fillConstituents(candidate, inputParticles, -1, RecoDecay::getMassPDG(candidate.pdgCode()));

      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        if (jet.eta() < trackEtaMin + jetR || jet.eta() > trackEtaMax - jetR) {
          continue;
        }
        if (jet.perp() < jetPtMin || jet.perp() >= jetPtMax) {
          continue;
        }
        isHFJet = false;
        std::vector<int> trackconst;
        std::vector<int> candconst;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == -1) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          for (const auto& constituent : jet.constituents()) {
            if (constituent.user_index() == -1)
              continue;
            trackconst.push_back(constituent.user_index());
          }

          candconst.push_back(candidate.globalIndex());
          trackConstituents(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
          hJetPhi->Fill(jet.phi());
          hJetEta->Fill(jet.rap());
          hJetNTracks->Fill(jet.constituents().size());
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processMCP, "HF jet finding on MC particle level", false);
};

using JetFinderHF = JetFinderHFTask<o2::aod::HFJets, o2::aod::HFJetConstituents>;
using MCParticleLevelJetFinderHF = JetFinderHFTask<o2::aod::MCParticleLevelHFJets, o2::aod::MCParticleLevelHFJetConstituents>;
using MCDetectorLevelJetFinderHF = JetFinderHFTask<o2::aod::MCDetectorLevelHFJets, o2::aod::MCDetectorLevelHFJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  auto hfjetMode = cfgc.options().get<std::string>("hfjetMode");

  if (hfjetMode.find("data") != std::string::npos || hfjetMode.empty())
    tasks.emplace_back(adaptAnalysisTask<JetFinderHF>(cfgc,
                                                      SetDefaultProcesses{{{"processData", true}, {"processMCP", false}, {"processMCD", false}}},
                                                      TaskName{"jet-finder-hf-data"}));

  if (hfjetMode.find("mcp") != std::string::npos || hfjetMode.empty())
    tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetFinderHF>(cfgc,
                                                                     SetDefaultProcesses{{{"processData", false}, {"processMCP", true}, {"processMCD", false}}},
                                                                     TaskName{"jet-finder-hf-mcp"}));

  if (hfjetMode.find("mcd") != std::string::npos || hfjetMode.empty())
    tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetFinderHF>(cfgc,
                                                                     SetDefaultProcesses{{{"processData", false}, {"processMCP", false}, {"processMCD", true}}},
                                                                     TaskName{"jet-finder-hf-mcd"}));

  return WorkflowSpec{tasks};
}
