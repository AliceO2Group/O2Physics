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

// heavy-flavour jet substructure tree filling task (subscribing to jet finder hf and jet substructure hf tasks)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename ParticleTable, typename CandidateTableData, typename CandidateTableMCD, typename JetTableData, typename OutputTableData, typename SubstructureOutputTableData, typename JetTableMCD, typename OutputTableMCD, typename SubstructureOutputTableMCD, typename JetTableMCP, typename OutputTableMCP, typename SubstructureOutputTableMCP>
struct JetSubstructureHFOutputTask {
  Produces<OutputTableData> jetOutputTableData;
  Produces<SubstructureOutputTableData> jetSubstructureOutputTableData;
  Produces<OutputTableMCD> jetOutputTableMCD;
  Produces<SubstructureOutputTableMCD> jetSubstructureOutputTableMCD;
  Produces<OutputTableMCP> jetOutputTableMCP;
  Produces<SubstructureOutputTableMCP> jetSubstructureOutputTableMCP;

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;
  }

  Filter jetSelection = aod::jet::pt >= jetPtMin;

  template <typename T, typename U, typename V, typename M, typename N>
  void fillTables(T const& collision, U const& jet, V const& cand, M& jetOutputTable, N& jetSubstructureOutputTable, std::vector<int> geoMatching, std::vector<int> ptMatching, std::vector<int> candMatching)
  {
    jetOutputTable(collision.globalIndex(), jet.globalIndex(), cand.globalIndex(), geoMatching, ptMatching, candMatching, jet.pt(), jet.phi(), jet.eta(), jet.r(), jet.tracks().size() + jet.hfcandidates().size());
    jetSubstructureOutputTable(jet.globalIndex(), jet.zg(), jet.rg(), jet.nsd());
  }

  void processDummy(aod::JCollision const& collision) {}
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processOutputData(aod::JCollision const& collision,
                         JetTableData const& jets,
                         aod::JTracks const& tracks,
                         CandidateTableData const& candidates)
  {
    std::vector<int> geoMatching{-1};
    std::vector<int> ptMatching{-1};
    std::vector<int> candMatching{-1};
    for (const auto& jet : jets) {
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          auto cands = jet.template hfcandidates_as<CandidateTableData>();
          auto cand = cands[0];
          fillTables(collision, jet, cand, jetOutputTableData, jetSubstructureOutputTableData, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputData, "hf jet substructure output Data", false);

  void processOutputMCD(aod::JCollision const& collision,
                        JetTableMCD const& mcdjets,
                        JetTableMCP const& mcpjets,
                        aod::JTracks const& tracks,
                        CandidateTableMCD const& candidates)
  {
    for (const auto& mcdjet : mcdjets) {
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (mcdjet.r() == round(jetRadiiValue * 100.0f)) {
          auto cands = mcdjet.template hfcandidates_as<CandidateTableMCD>();
          auto cand = cands[0];
          std::vector<int> geoMatching;
          std::vector<int> ptMatching;
          std::vector<int> candMatching;
          for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCP>()) {
            geoMatching.push_back(mcpjet.globalIndex());
          }
          for (auto& mcpjet : mcdjet.template matchedJetPt_as<JetTableMCP>()) {
            ptMatching.push_back(mcpjet.globalIndex());
          }
          for (auto& mcpjet : mcdjet.template matchedJetCand_as<JetTableMCP>()) {
            candMatching.push_back(mcpjet.globalIndex());
          }
          fillTables(collision, mcdjet, cand, jetOutputTableMCD, jetSubstructureOutputTableMCD, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputMCD, "hf jet substructure output MCD", false);

  void processOutputMCP(aod::JMcCollision const& collision,
                        JetTableMCP const& mcpjets,
                        JetTableMCD const& mcdjets,
                        ParticleTable const& particles)
  {
    for (const auto& mcpjet : mcpjets) {
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (mcpjet.r() == round(jetRadiiValue * 100.0f)) {
          auto cands = mcpjet.template hfcandidates_as<ParticleTable>();
          auto cand = cands[0];
          std::vector<int> geoMatching;
          std::vector<int> ptMatching;
          std::vector<int> candMatching;
          for (auto& mcdjet : mcpjet.template matchedJetGeo_as<JetTableMCD>()) {
            geoMatching.push_back(mcdjet.globalIndex());
          }
          for (auto& mcdjet : mcpjet.template matchedJetPt_as<JetTableMCD>()) {
            ptMatching.push_back(mcdjet.globalIndex());
          }
          for (auto& mcdjet : mcpjet.template matchedJetCand_as<JetTableMCD>()) {
            candMatching.push_back(mcdjet.globalIndex());
          }
          fillTables(collision, mcpjet, cand, jetOutputTableMCP, jetSubstructureOutputTableMCP, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputMCP, "hf jet substructure output MCP", false);
};
using JetSubstructureOutputD0 = JetSubstructureHFOutputTask<soa::Join<aod::JMcParticles, aod::HfCand2ProngMcGen>, soa::Join<aod::HfCand2Prong, aod::HfSelD0>, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>, soa::Filtered<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0ChargedJetSubstructures>>, aod::D0ChargedJetOutput, aod::D0ChargedJetSubstructureOutput, soa::Filtered<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetSubstructures, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>>, aod::D0ChargedMCDetectorLevelJetOutput, aod::D0ChargedMCDetectorLevelJetSubstructureOutput, soa::Filtered<soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetSubstructures, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>>, aod::D0ChargedMCParticleLevelJetOutput, aod::D0ChargedMCParticleLevelJetSubstructureOutput>;
// using JetSubstructureOutputLc = JetSubstructureHFOutputTask<soa::Join<aod::JMcParticles, aod::HfCand3ProngMcGen>, soa::Join<aod::HfCand3Prong, aod::HfSelLc>, soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>, soa::Filtered<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcChargedJetSubstructures>>, aod::LcChargedJetOutput, aod::LcChargedJetSubstructureOutput,soa::Filtered<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetSubstructures,aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>>, aod::LcChargedMCDetectorLevelJetOutput, aod::LcChargedMCDetectorLevelJetSubstructureOutput, soa::Filtered<soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetSubstructures,aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>>, aod::LcChargedMCParticleLevelJetOutput, aod::LcChargedMCParticleLevelJetSubstructureOutput>;
// using JetSubstructureOutputBplus = JetSubstructureHFOutputTask<soa::Join<aod::JMcParticles, aod::HfCandBplusMcGen>, soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>, soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>, soa::Filtered<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusChargedJetSubstructures>>, aod::BplusChargedJetOutput, aod::BplusChargedJetSubstructureOutput,soa::Filtered<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelJetSubstructures,aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>>, aod::BplusChargedMCDetectorLevelJetOutput, aod::BplusChargedMCDetectorLevelJetSubstructureOutput, soa::Filtered<soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelJetSubstructures,aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>>, aod::BplusChargedMCParticleLevelJetOutput, aod::BplusChargedMCParticleLevelJetSubstructureOutput>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputD0>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructured0-output"}));
  // tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputLc>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructurelc-output"}));
  // tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputBplus>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructurebplus-output"}));
  return WorkflowSpec{tasks};
}
