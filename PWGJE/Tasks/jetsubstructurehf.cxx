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

template <typename JetTable, typename CandidateTable, typename JetTableMCP, typename SubstructureTable>
struct JetSubstructureHFTask {
  Produces<SubstructureTable> jetSubstructurehfTable;
  OutputObj<TH2F> hZg{"h_jet_zg_jet_pt"};
  OutputObj<TH2F> hRg{"h_jet_rg_jet_pt"};
  OutputObj<TH2F> hNsd{"h_jet_nsd_jet_pt"};

  // Jet level configurables
  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};

  Service<o2::framework::O2DatabasePDG> pdg;
  int candPDG;

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  void init(InitContext const&)
  {
    hZg.setObject(new TH2F("h_jet_zg_jet_pt", ";z_{g}; #it{p}_{T,jet} (GeV/#it{c})",
                           10, 0.0, 0.5, 200, 0.0, 200.0));
    hRg.setObject(new TH2F("h_jet_rg_jet_pt", ";R_{g}; #it{p}_{T,jet} (GeV/#it{c})",
                           10, 0.0, 0.5, 200, 0.0, 200.0));
    hNsd.setObject(new TH2F("h_jet_nsd_jet_pt", ";n_{SD}; #it{p}_{T,jet} (GeV/#it{c})",
                            7, -0.5, 6.5, 200, 0.0, 200.0));

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;

    if constexpr (std::is_same_v<std::decay_t<JetTableMCP>, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>>) {
      candPDG = static_cast<int>(pdg::Code::kD0);
    }
    if constexpr (std::is_same_v<std::decay_t<JetTableMCP>, soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>>) {
      candPDG = static_cast<int>(pdg::Code::kBPlus);
    }
    if constexpr (std::is_same_v<std::decay_t<JetTableMCP>, soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>>) {
      candPDG = static_cast<int>(pdg::Code::kLambdaCPlus);
    }
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
      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        if (!softDropped) {
          zg = z;
          rg = theta;
          hZg->Fill(zg, jet.pt());
          hRg->Fill(rg, jet.pt());
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
    hNsd->Fill(nsd, jet.pt());
    jetSubstructurehfTable(zg, rg, nsd);
  }

  void processDummy(aod::Tracks const& tracks)
  {
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsHF(typename JetTable::iterator const& jet,
                            CandidateTable const& candidates,
                            aod::Tracks const& tracks)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<aod::Tracks>()) {
      FastJetUtilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    for (auto& jetHFCandidate : jet.template hfcandidates_as<CandidateTable>()) { // should only be one at the moment
      FastJetUtilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsHF, "HF jet substructure", false);

  void processChargedJetsHFMCP(typename JetTableMCP::iterator const& jet,
                               aod::McParticles const& particles)
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
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsHFMCP, "HF jet substructure on MC particle level", false);
};
using JetSubstructureD0 = JetSubstructureHFTask<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>, soa::Join<aod::HfCand2Prong, aod::HfSelD0>, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>, o2::aod::D0ChargedJetSubstructures>;
using MCDetectorLevelJetSubstructureD0 = JetSubstructureHFTask<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>, o2::aod::D0ChargedMCDetectorLevelJetSubstructures>;
using MCParticleLevelJetSubstructureD0 = JetSubstructureHFTask<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>, o2::aod::D0ChargedMCParticleLevelJetSubstructures>;

// using JetSubstructureBplus = JetSubstructureHFTask<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents>,soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>,soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>,o2::aod::BplusChargedJetSubstructures>;
// using MCDetectorLevelJetSubstructureBplus = JetSubstructureHFTask<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>,soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>,soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>,o2::aod::BplusChargedMCDetectorLevelJetSubstructures>;
// using MCParticleLevelJetSubstructureBplus = JetSubstructureHFTask<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>,soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>,soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>,o2::aod::BplusChargedMCParticleLevelJetSubstructures>;

// using JetSubstructureLc = JetSubstructureHFTask<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents>, soa::Join<aod::HfCand3Prong, aod::HfSelLc>, soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>, o2::aod::LcChargedJetSubstructures>;
//  using MCDetectorLevelJetSubstructureLc = JetSubstructureHFTask<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>,soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>,soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>,o2::aod::LcChargedMCDetectorLevelJetSubstructures>;
//  using MCParticleLevelJetSubstructureLc = JetSubstructureHFTask<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>,soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>,soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>,o2::aod::LcChargedMCParticleLevelJetSubstructures>;

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
    /*
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
    */
    return WorkflowSpec{tasks};
}
