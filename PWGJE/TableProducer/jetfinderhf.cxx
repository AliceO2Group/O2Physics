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
// Author: Nima Zardoshti
/// \author Jochen Klein <jochen.klein@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

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
  OutputObj<TH1F> hJetPtTrue{"h_jet_pt_true"};
  OutputObj<TH1F> hD0Pt{"h_D0_pt"};

  Service<TDatabasePDG> pdg;

  TrackSelection globalTracks;

  std::vector<fastjet::PseudoJet> jets;
  std::vector<fastjet::PseudoJet> inputParticles;
  JetFinder jetFinder;

  void init(InitContext const&)
  {
    // set up global tracks and adjust as necessary
    globalTracks = getGlobalTrackSelection();
    globalTracks.SetEtaRange(-.9, .9);

    hJetPt.setObject(new TH1F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                              100, 0., 100.));
    hJetPtTrue.setObject(new TH1F("h_jet_pt_true", "jet p_{T};p_{T} (GeV/#it{c})",
                                  100, 0., 100.));
    hD0Pt.setObject(new TH1F("h_D0_pt", "jet p_{T,D};p_{T,D} (GeV/#it{c})",
                             60, 0., 60.));
  }

  Configurable<int> d_selectionFlagD0{"d_selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> d_selectionFlagD0bar{"d_selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  //need enum as configurable
  enum pdgCode { pdgD0 = 421 };

  Filter trackCuts = (aod::track::pt > 0.15f && aod::track::eta > -0.9f && aod::track::eta < 0.9f);
  Filter partCuts = (aod::mcparticle::pt > 0.15f && aod::mcparticle::eta > -0.9f && aod::mcparticle::eta < 0.9f);
  Filter seltrack = (aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar);

  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                   soa::Filtered<aod::Tracks> const& tracks,
                   soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>> const& candidates)
  {
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;

    if (!collision.sel8())
      return;

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      jets.clear();
      inputParticles.clear();
      for (auto& track : tracks) {
        auto energy = std::sqrt(track.p() * track.p() + JetFinder::mPion * JetFinder::mPion);
        if (candidate.index0().globalIndex() == track.globalIndex() || candidate.index1().globalIndex() == track.globalIndex()) { //is it global index?
          continue;
        }
        inputParticles.emplace_back(track.px(), track.py(), track.pz(), energy);
        inputParticles.back().set_user_index(track.globalIndex());
      }
      inputParticles.emplace_back(candidate.px(), candidate.py(), candidate.pz(), candidate.e(RecoDecay::getMassPDG(pdgD0)));
      inputParticles.back().set_user_index(1);

      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        isHFJet = false;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == 1 && (candidate.isSelD0() == 1 || candidate.isSelD0bar() == 1)) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          // for (const auto& constituent : jet.constituents()) {
          // trackConstituents(jetsTable.lastIndex(), constituent.user_index());
          // }
          hJetPt->Fill(jet.pt());
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processData, "HF jet finding on data", true);

  void processMCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                  soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>> const& tracks,
                  soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate, aod::HfCandProng2MCRec>> const& candidates)
  {
    LOG(debug) << "Per Event MCP";
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;

    //this loop should be made more efficient
    // TODO: should probably refine the candidate selection
    for (auto& candidate : candidates) {
      jets.clear();
      inputParticles.clear();
      for (auto& track : tracks) {
        if (!globalTracks.IsSelected(track)) {
          LOGF(info, "Rejecting track %d with track cuts", track.globalIndex());
          continue;
        }
        if (candidate.index0().globalIndex() == track.globalIndex() || candidate.index1().globalIndex() == track.globalIndex()) {
          LOGF(info, "Rejecting track %d as daughter of candidate %d", track.globalIndex(), candidate.globalIndex());
          continue;
        }
        // LOGF(info, "Adding track %d with pt %g", track.globalIndex(), track.pt());
        auto energy = std::sqrt(track.p() * track.p() + JetFinder::mPion * JetFinder::mPion);
        inputParticles.emplace_back(track.px(), track.py(), track.pz(), energy);
        inputParticles.back().set_user_index(track.globalIndex());
      }
      inputParticles.emplace_back(candidate.px(), candidate.py(), candidate.pz(), candidate.e(RecoDecay::getMassPDG(pdgD0)));
      inputParticles.back().set_user_index(1);

      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        isHFJet = false;
        std::vector<int> trackconst;
        std::vector<int> candconst;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == 1 && (candidate.isSelD0() == 1 || candidate.isSelD0bar() == 1)) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          const auto& constituents = sorted_by_pt(jet.constituents());
          for (const auto& constituent : constituents) {
            if (constituent.user_index() != 1) {
              LOGF(info, "jet %d (coll %d) has constituent %d", jetsTable.lastIndex(), collision.globalIndex(), constituent.user_index());
              auto track = tracks.rawIteratorAt(constituent.user_index());
              LOGF(info, "constituent %d points to track %d (coll %d)", constituent.user_index(), track.globalIndex(), 0); // , track.collisionId()); // .globalIndex());
              trackconst.push_back(constituent.user_index());
            }
          }
          LOGF(info, "jet %d (coll %d) has candidate %d-%d", jetsTable.lastIndex(), collision.globalIndex(), candidate.index(), candidate.globalIndex());
          candconst.push_back(candidate.globalIndex());
          trackConstituents(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
          if (candidate.flagMCMatchRec() & (1 << aod::hf_cand_prong2::DecayType::D0ToPiK))
            hJetPtTrue->Fill(jet.pt());
          hD0Pt->Fill(candidate.pt());
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processMCD, "HF jet finding on MC detector level", false);

  void processMCP(aod::McCollision const& collision,
                  soa::Filtered<soa::Join<aod::McParticles, aod::HfCandProng2MCGen>> const& particles)
  {
    LOG(debug) << "Per Event MCP";
    // TODO: retrieve pion mass from somewhere
    bool isHFJet;

    // TODO: probably should do this as a filter
    std::vector<soa::Filtered<soa::Join<aod::McParticles, aod::HfCandProng2MCGen>>::iterator> candidates;
    for (auto const& part : particles) {
      // TODO: generalise to any D0
      if (std::abs(part.flagMCMatchGen()) & (1 << aod::hf_cand_prong2::DecayType::D0ToPiK)) {
        candidates.push_back(part);
        LOGF(info, "MC candidate %d -> %d", part.globalIndex(), candidates.back().globalIndex());
      }
    }

    //this loop should be made more efficient
    for (auto& candidate : candidates) {
      jets.clear();
      inputParticles.clear();
      for (auto& track : particles) {
        // exclude neutral particles
        // TODO: can we do this through the filter?
        auto p = pdg->GetParticle(track.pdgCode());
        // LOGF(info, "Checking particle %i with status %d and charge %g",
        //   track.globalIndex(), track.getGenStatusCode(), p ? std::abs(p->Charge()) : -999.);
        if ((track.getGenStatusCode() != 1) || (p ? std::abs(p->Charge()) : 0.) < 3.) {
          LOGF(info, "Rejecting particle %d with status %d and charge %g", track.globalIndex(), track.getGenStatusCode(), p ? std::abs(p->Charge()) : -999.);
          continue;
        }

        // TODO: check what mass to use?
        auto energy = std::sqrt(track.p() * track.p() + JetFinder::mPion * JetFinder::mPion);
        const auto daughters = candidate.daughtersIds();
        if (std::find(std::begin(daughters), std::end(daughters), track.globalIndex()) != std::end(daughters)) {
          LOGF(info, "Rejecting particle %d as daughter of candidate %d", track.globalIndex(), candidate.globalIndex());
          continue;
        }
        // LOGF(info, "Adding particle %d with pt %g", track.globalIndex(), track.pt());
        inputParticles.emplace_back(track.px(), track.py(), track.pz(), energy);
        inputParticles.back().set_user_index(track.globalIndex());
      }
      inputParticles.emplace_back(candidate.px(), candidate.py(), candidate.pz(), candidate.e());
      inputParticles.back().set_user_index(1);

      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        isHFJet = false;
        std::vector<int> trackconst;
        std::vector<int> candconst;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.user_index() == 1) {
            isHFJet = true;
            break;
          }
        }
        if (isHFJet) {
          jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                    jet.E(), jet.m(), jet.area(), jetFinder.jetR);
          for (const auto& constituent : jet.constituents()) {
            if (constituent.user_index() == 1)
              continue;
            LOGF(info, "MC jet %d (MC coll %d) has constituent %d", jetsTable.lastIndex(), collision.globalIndex(), constituent.user_index());
            trackconst.push_back(constituent.user_index());
          }
          LOGF(info, "MC jet %d (MC coll %d) has candidate %d",
               jetsTable.lastIndex(), collision.globalIndex(), candidates.back().globalIndex());
          candconst.push_back(candidate.globalIndex());
          LOGF(info, "MC jet %d has %d track and %d candidate constituents", jetsTable.lastIndex(), trackconst.size(), candconst.size());
          trackConstituents(jetsTable.lastIndex(), trackconst, std::vector<int>(), candconst);
          hJetPt->Fill(jet.pt());
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
