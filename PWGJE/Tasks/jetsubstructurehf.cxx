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
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename SubstructureTable>
struct JetSubstructureHFTask {
  Produces<SubstructureTable> jetSubstructurehfTable;
  OutputObj<TH1F> hZg{"h_jet_zg"};
  OutputObj<TH1F> hRg{"h_jet_rg"};
  OutputObj<TH1F> hNsd{"h_jet_nsd"};

  // HF candidate level configurables
  Configurable<std::string> candSpecie_s{"candSpecie_s", "D0", "options are D0, Lc, Bplus"};

  // Jet level configurables
  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<bool> doConstSub{"doConstSub", false, "do constituent subtraction"};

  Service<O2DatabasePDG> pdg;
  int candPDG;

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

    auto candSpecie = static_cast<std::string>(candSpecie_s);
    if (candSpecie == "D0") {
      candPDG = static_cast<int>(pdg::Code::kD0);
    }
    if (candSpecie == "Bplus") {
      candPDG = static_cast<int>(pdg::Code::kBPlus);
    }
    if (candSpecie == "Lc") {
      candPDG = static_cast<int>(pdg::Code::kLambdaCPlus);
    }
    if (candSpecie == "JPsi") {
      candPDG = static_cast<int>(pdg::Code::kJPsi);
    }
  }

  template <typename T>
  void processConstituents(T const& jet)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<aod::Tracks>()) {
      FastJetUtilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
  }

  template <typename T>
  void processGenLevel(T const& jet)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<aod::McParticles>()) {
      FastJetUtilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), RecoDecay::getMassPDG(jetConstituent.pdgCode()));
    }
    for (auto& jetHFCandidate : jet.template hfcandidates_as<aod::McParticles>()) { // should only be one at the moment
      FastJetUtilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(jetHFCandidate.pdgCode()));
    }
    jetReclustering(jet);
  }

  template <typename T>
  void jetReclustering(T const& jet)
  {
    jetReclustered.clear();
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
        if (subjet1Constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
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
    jetSubstructurehfTable(jet.globalIndex(), zg, rg, nsd);
  }

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processDummy, "Dummy process function turned on by default", true);

  void processD0Data(soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>::iterator const& jet,
                     soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
                     aod::Tracks const& tracks)
  {
    processConstituents(jet);
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>()) { // should only be one at the moment
      FastJetUtilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processD0Data, "D0 jet substructure on data", false);

  void processD0MCD(soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                    soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const& candidates,
                    aod::Tracks const& tracks)
  {
    processConstituents(jet);
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>()) { // should only be one at the moment
      FastJetUtilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processD0MCD, "D0 jet substructure on MC detector level", false);

  void processD0MCP(soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>::iterator const& jet,
                    aod::McParticles const& particles)
  {
    processGenLevel(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processD0MCP, "D0 jet substructure on MC particle level", false);

  void processLcData(soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents>::iterator const& jet,
                     soa::Join<aod::HfCand2Prong, aod::HfSelLc> const& candidates,
                     aod::Tracks const& tracks)
  {
    processConstituents(jet);
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelLc>>()) { // should only be one at the moment
      FastJetUtilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processLcData, "Lc jet substructure on data", false);

  void processLcMCD(soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                    soa::Join<aod::HfCand2Prong, aod::HfSelLc, aod::HfCand2ProngMcRec> const& candidates,
                    aod::Tracks const& tracks)
  {
    processConstituents(jet);
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelLc, aod::HfCand2ProngMcRec>>()) { // should only be one at the moment
      FastJetUtilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processLcMCD, "Lc jet substructure on MC detector level", false);

  void processLcMCP(soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>::iterator const& jet,
                    aod::McParticles const& particles)
  {
    processGenLevel(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processLcMCP, "Lc jet substructure on MC particle level", false);

  void processBplusData(soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents>::iterator const& jet,
                        soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi> const& candidates,
                        aod::Tracks const& tracks)
  {
    processConstituents(jet);
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>>()) { // should only be one at the moment
      FastJetUtilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processBplusData, "Bplus jet substructure on data", false);

  void processBplusMCD(soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                       soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec> const& candidates,
                       aod::Tracks const& tracks)
  {
    processConstituents(jet);
    for (auto& jetHFCandidate : jet.hfcandidates_as<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>>()) { // should only be one at the moment
      FastJetUtilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processBplusMCD, "Bplus jet substructure on MC detector level", false);

  void processBplusMCP(soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>::iterator const& jet,
                       aod::McParticles const& particles)
  {
    processGenLevel(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processBplusMCP, "Bplus jet substructure on MC particle level", false);
};
using JetSubstructureD0 = JetSubstructureHFTask<o2::aod::D0ChargedJetSubstructure>;
using MCDetectorLevelJetSubstructureD0 = JetSubstructureHFTask<o2::aod::D0ChargedMCDetectorLevelJetSubstructure>;
using MCParticleLevelJetSubstructureD0 = JetSubstructureHFTask<o2::aod::D0ChargedMCParticleLevelJetSubstructure>;
using JetSubstructureLc = JetSubstructureHFTask<o2::aod::LcChargedJetSubstructure>;
using MCDetectorLevelJetSubstructureLc = JetSubstructureHFTask<o2::aod::LcChargedMCDetectorLevelJetSubstructure>;
using MCParticleLevelJetSubstructureLc = JetSubstructureHFTask<o2::aod::LcChargedMCParticleLevelJetSubstructure>;
using JetSubstructureBplus = JetSubstructureHFTask<o2::aod::BplusChargedJetSubstructure>;
using MCDetectorLevelJetSubstructureBplus = JetSubstructureHFTask<o2::aod::BplusChargedMCDetectorLevelJetSubstructure>;
using MCParticleLevelJetSubstructureBplus = JetSubstructureHFTask<o2::aod::BplusChargedMCParticleLevelJetSubstructure>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureD0>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-substructure-D0-data"}));

  tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetSubstructureD0>(cfgc,
                                                                         SetDefaultProcesses{},
                                                                         TaskName{"jet-substructure-D0-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetSubstructureD0>(cfgc,
                                                                         SetDefaultProcesses{},
                                                                         TaskName{"jet-substructure-D0-mcp"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureLc>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-substructure-Lc-data"}));

  tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetSubstructureLc>(cfgc,
                                                                         SetDefaultProcesses{},
                                                                         TaskName{"jet-substructure-Lc-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetSubstructureLc>(cfgc,
                                                                         SetDefaultProcesses{},
                                                                         TaskName{"jet-substructure-Lc-mcp"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureBplus>(cfgc,
                                                             SetDefaultProcesses{},
                                                             TaskName{"jet-substructure-Bplus-data"}));

  tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetSubstructureBplus>(cfgc,
                                                                            SetDefaultProcesses{},
                                                                            TaskName{"jet-substructure-Bplus-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetSubstructureBplus>(cfgc,
                                                                            SetDefaultProcesses{},
                                                                            TaskName{"jet-substructure-Bplus-mcp"}));

  return WorkflowSpec{tasks};
}
