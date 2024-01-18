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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

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
#include "PWGJE/Core/JetHFUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename JetTableData, typename JetTableMCD, typename JetTableMCP, typename JetTableDataSub, typename CandidateTable, typename CandidateTableMCP, typename SubstructureTableData, typename SubstructureTableMCD, typename SubstructureTableMCP, typename SubstructureTableDataSub, typename TracksSub>
struct JetSubstructureHFTask {
  Produces<SubstructureTableData> jetSubstructureDataTable;
  Produces<SubstructureTableMCD> jetSubstructureMCDTable;
  Produces<SubstructureTableMCP> jetSubstructureMCPTable;
  Produces<SubstructureTableDataSub> jetSubstructureDataSubTable;
  OutputObj<TH2F> hZg{"h_jet_zg_jet_pt"};
  OutputObj<TH2F> hRg{"h_jet_rg_jet_pt"};
  OutputObj<TH2F> hNsd{"h_jet_nsd_jet_pt"};

  // Jet level configurables
  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};

  Service<o2::framework::O2DatabasePDG> pdg;
  int candMass;

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

    candMass = jethfutilities::getTablePDGMass<CandidateTable>();
  }

  template <typename T, typename U>
  void jetReclustering(T const& jet, U& outputTable)
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
        if (subjet1Constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
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
    outputTable(zg, rg, nsd);
  }

  template <typename T, typename U, typename V, typename M>
  void analyseCharged(T const& jet, U const& tracks, V const& candidates, M& outputTable)
  {

    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    for (auto& jetHFCandidate : jet.template hfcandidates_as<V>()) { // should only be one at the moment
      fastjetutilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), candMass);
    }
    jetReclustering(jet, outputTable);
  }

  void processDummy(JetTracks const& tracks)
  {
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(typename JetTableData::iterator const& jet,
                              CandidateTable const& candidates,
                              JetTracks const& tracks)
  {
    analyseCharged(jet, tracks, candidates, jetSubstructureDataTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsData, "HF jet substructure on data", false);

  void processChargedJetsDataSub(typename JetTableDataSub::iterator const& jet,
                                 CandidateTable const& candidates,
                                 TracksSub const& tracks)
  {
    analyseCharged(jet, tracks, candidates, jetSubstructureDataSubTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsDataSub, "HF jet substructure on data", false);

  void processChargedJetsMCD(typename JetTableMCD::iterator const& jet,
                             CandidateTable const& candidates,
                             JetTracks const& tracks)
  {
    analyseCharged(jet, tracks, candidates, jetSubstructureMCDTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsMCD, "HF jet substructure on data", false);

  void processChargedJetsMCP(typename JetTableMCP::iterator const& jet,
                             JetParticles const& particles,
                             CandidateTableMCP const& candidates)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<JetParticles>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
    }
    for (auto& jetHFCandidate : jet.template hfcandidates_as<CandidateTableMCP>()) {
      fastjetutilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), candMass);
    }
    jetReclustering(jet, jetSubstructureMCPTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsMCP, "HF jet substructure on MC particle level", false);
};
using JetSubstructureD0 = JetSubstructureHFTask<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>, soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>, soa::Join<aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents>, CandidatesD0Data, CandidatesD0MCP, aod::D0ChargedJetSubstructures, aod::D0ChargedMCDetectorLevelJetSubstructures, aod::D0ChargedMCParticleLevelJetSubstructures, aod::D0ChargedEventWiseSubtractedJetSubstructures, aod::JTrackD0Subs>;
// using JetSubstructureLc = JetSubstructureHFTask<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents>,soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>,soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>,soa::Join<aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents>, CandidatesLcData, CandidatesLcMCP, aod::LcChargedJetSubstructures,aod::LcChargedMCDetectorLevelJetSubstructures,aod::LcChargedMCParticleLevelJetSubstructures, aod::LcChargedEventWiseSubtractedJetSubstructures, aod::JTrackLcSubs>;
// using JetSubstructureBplus = JetSubstructureHFTask<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents>,soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>,soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>,soa::Join<aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents>, CandidatesBplusData, CandidatesBplusMCP, aod::BplusChargedJetSubstructures,aod::BplusChargedMCDetectorLevelJetSubstructures,aod::BplusChargedMCParticleLevelJetSubstructures, aod::BplusChargedEventWiseSubtractedJetSubstructures, aod::JTrackBplusSubs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureD0>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-substructure-D0"}));
  /*
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureLc>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-substructure-Lc"}));

    tasks.emplace_back(adaptAnalysisTask<JetSubstructureBplus>(cfgc,
                                                               SetDefaultProcesses{},
                                                               TaskName{"jet-substructure-Bplus"}));
  */
  return WorkflowSpec{tasks};
}
