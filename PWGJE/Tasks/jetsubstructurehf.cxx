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

// heavy-flavour jet substructure task (subscribing to jet finder hf task)
//
// Author: Nima Zardoshti
//

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec hfjetsubstructureMode = {
    "hfjetsubstructureMode",
    VariantType::String,
    "",
    {"HF jet substrcture mode."},
  };
  workflowOptions.push_back(hfjetsubstructureMode);
}

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

namespace o2::aod
{
namespace jetsubstructurehf
{
// add an index to the jet table
// add a coloumn to make it unique for data, MC det and MC gen
DECLARE_SOA_COLUMN(Zg, zg, float);
DECLARE_SOA_COLUMN(Rg, rg, float);
DECLARE_SOA_COLUMN(Nsd, nsd, float);
} // namespace jetsubstructurehf
DECLARE_SOA_TABLE(JetSubtructureHF, "AOD", "JETSUBSTRUCTHF", jetsubstructurehf::Zg, jetsubstructurehf::Rg, jetsubstructurehf::Nsd);
} // namespace o2::aod

template <typename JetTable, typename TrackConstituentTable> // add one more for output table type
struct JetSubstructureHFTask {
  // Produces<aod::JetSubtructureHF> jetSubstructurehf;
  OutputObj<TH1F> hZg{"h_jet_zg"};
  OutputObj<TH1F> hRg{"h_jet_rg"};
  OutputObj<TH1F> hNsd{"h_jet_nsd"};

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<bool> doConstSub{"doConstSub", false, "do constituent subtraction"};

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  void init(InitContext const&)
  {
    hZg.setObject(new TH1F("h_jet_zg", "zg ;zg",
                           8, 0.1, 0.5));
    hRg.setObject(new TH1F("h_jet_rg", "rg ;rg",
                           10, 0.0, 0.5));
    hNsd.setObject(new TH1F("h_jet_nsd", "nsd ;nsd",
                            7, -0.5, 6.5));
    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
  }
  enum pdgCode { pdgD0 = 421 };

  // Filter jetCuts = aod::jet::pt > f_jetPtMin;

  void processData(soa::Join<aod::HFJets, aod::HFJetConstituents>::iterator const& jet, // add template back
                   soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
                   aod::Tracks const& tracks)
  {
    jetConstituents.clear();
    jetReclustered.clear();
    // if (b_DoConstSub) {
    // for (const auto& constituent : constituentsSub) {
    //  fillConstituents(constituent, jetConstituents);
    //  }
    // } else {
    //  for (auto& jetConstituentIndex : jet.trackIds()) {
    //  auto jetConstituent = tracks.rawIteratorAt(jetConstituentIndex - tracks.offset());
    for (auto& jetConstituent : jet.tracks_as<aod::Tracks>()) {
      fillConstituents(jetConstituent, jetConstituents);
    }
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>()) { // should only be one at the moment
      fillConstituents(jetHFCandidate, jetConstituents, -1, RecoDecay::getMassPDG(pdgD0));
    }
    //}
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    auto nsd = 0.0;
    auto zg = -1.0;
    auto rg = -1.0;
    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);
      }
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      if (z >= zCut * TMath::Power(theta / jetR, beta)) {
        if (!softDropped) {
          zg = z;
          rg = theta;
          hZg->Fill(zg);
          hRg->Fill(rg);
          softDropped = true;
        }
        nsd++;
      }
      bool isHFInSubjet1 = false;
      for (auto& subjet1Constituent : parentSubJet1.constituents()) {
        if (subjet1Constituent.user_index() == -1) {
          isHFInSubjet1 = true;
          break;
        }
      }
      if (isHFInSubjet1) {
        daughterSubJet = parentSubJet1;
      } else {
        daughterSubJet = parentSubJet2;
      }
    }
    hNsd->Fill(nsd);
    // jetSubstructurehf(zg, rg, nsd);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processData, "HF jet substructure on data", true);

  void processMCD(soa::Join<aod::MCDetectorLevelHFJets, aod::MCDetectorLevelHFJetConstituents>::iterator const& jet,
                  soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const& candidates,
                  // aod::JetConstituentsSub const& constituentsSub,
                  aod::Tracks const& tracks)
  {
    jetConstituents.clear();
    jetReclustered.clear();
    // if (b_DoConstSub) {
    // for (const auto& constituent : constituentsSub) {
    //  fillConstituents(constituent, jetConstituents);
    //  }
    // } else {
    for (auto& jetConstituent : jet.tracks_as<aod::Tracks>()) {
      fillConstituents(jetConstituent, jetConstituents);
    }
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>()) { // should only be one at the moment
      fillConstituents(jetHFCandidate, jetConstituents, -1, RecoDecay::getMassPDG(pdgD0));
    }
    //}
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    auto nsd = 0.0;
    auto zg = -1.0;
    auto rg = -1.0;
    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);
      }
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      if (z >= zCut * TMath::Power(theta / jetR, beta)) {
        if (!softDropped) {
          zg = z;
          rg = theta;
          hZg->Fill(zg);
          hRg->Fill(rg);
          softDropped = true;
        }
        nsd++;
      }
      bool isHFInSubjet1 = false;
      for (auto& subjet1Constituent : parentSubJet1.constituents()) {
        if (subjet1Constituent.user_index() == -1) {
          isHFInSubjet1 = true;
          break;
        }
      }
      if (isHFInSubjet1) {
        daughterSubJet = parentSubJet1;
      } else {
        daughterSubJet = parentSubJet2;
      }
    }
    hNsd->Fill(nsd);
    // jetSubstructurehf(zg, rg, nsd);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processMCD, "HF jet substructure on MC detector level", true);

  void processMCP(soa::Join<aod::MCParticleLevelHFJets, aod::MCParticleLevelHFJetConstituents>::iterator const& jet,
                  // aod::JetConstituentsSub const& constituentsSub,
                  soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& particles)
  {
    jetConstituents.clear();
    jetReclustered.clear();
    // if (b_DoConstSub) {
    // for (const auto& constituent : constituentsSub) {
    //  fillConstituents(constituent, jetConstituents);
    //  }
    // } else {
    for (auto& jetConstituent : jet.tracks_as<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>()) {
      fillConstituents(jetConstituent, jetConstituents, -2, RecoDecay::getMassPDG(jetConstituent.pdgCode()));
    }
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>()) { // should only be one at the moment
      fillConstituents(jetHFCandidate, jetConstituents, -1, RecoDecay::getMassPDG(jetHFCandidate.pdgCode()));
    }
    //}
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    auto nsd = 0.0;
    auto zg = -1.0;
    auto rg = -1.0;
    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);
      }
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      if (z >= zCut * TMath::Power(theta / jetR, beta)) {
        if (!softDropped) {
          zg = z;
          rg = theta;
          hZg->Fill(zg);
          hRg->Fill(rg);
          softDropped = true;
        }
        nsd++;
      }
      bool isHFInSubjet1 = false;
      for (auto& subjet1Constituent : parentSubJet1.constituents()) {
        if (subjet1Constituent.user_index() == -1) {
          isHFInSubjet1 = true;
          break;
        }
      }
      if (isHFInSubjet1) {
        daughterSubJet = parentSubJet1;
      } else {
        daughterSubJet = parentSubJet2;
      }
    }
    hNsd->Fill(nsd);
    // jetSubstructurehf(zg, rg, nsd);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processMCP, "HF jet substructure on MC particle level", true);
};
using JetSubstructureHFData = JetSubstructureHFTask<aod::HFJets, aod::HFJetConstituents>;
using JetSubstructureHFMCParticleLevel = JetSubstructureHFTask<aod::MCParticleLevelHFJets, aod::MCParticleLevelHFJetConstituents>;
using JetSubstructureHFMCDetectorLevel = JetSubstructureHFTask<aod::MCDetectorLevelHFJets, aod::MCDetectorLevelHFJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  auto hfjetsubstructureMode = cfgc.options().get<std::string>("hfjetsubstructureMode");

  if (hfjetsubstructureMode.find("data") != std::string::npos || hfjetsubstructureMode.empty())
    tasks.emplace_back(adaptAnalysisTask<JetSubstructureHFData>(cfgc,
                                                                SetDefaultProcesses{{{"processData", true}, {"processMCP", false}, {"processMCD", false}}},
                                                                TaskName{"jet-substructure-hf-data"}));

  if (hfjetsubstructureMode.find("mcp") != std::string::npos || hfjetsubstructureMode.empty())
    tasks.emplace_back(adaptAnalysisTask<JetSubstructureHFMCParticleLevel>(cfgc,
                                                                           SetDefaultProcesses{{{"processData", false}, {"processMCP", true}, {"processMCD", false}}},
                                                                           TaskName{"jet-substructure-hf-mcp"}));

  if (hfjetsubstructureMode.find("mcd") != std::string::npos || hfjetsubstructureMode.empty())
    tasks.emplace_back(adaptAnalysisTask<JetSubstructureHFMCDetectorLevel>(cfgc,
                                                                           SetDefaultProcesses{{{"processData", false}, {"processMCP", false}, {"processMCD", true}}},
                                                                           TaskName{"jet-substructure-hf-mcd"}));

  return WorkflowSpec{tasks};
}
