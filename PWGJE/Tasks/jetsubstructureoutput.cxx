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

// jet substructure tree filling task (subscribing to jet finder hf and jet substructure tasks)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename JetTableData, typename OutputTableData, typename SubstructureOutputTableData, typename JetTableMCD, typename OutputTableMCD, typename SubstructureOutputTableMCD, typename JetTableMCP, typename OutputTableMCP, typename SubstructureOutputTableMCP>
struct JetSubstructureOutputTask {
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

  template <typename T, typename U, typename V, typename M>
  void fillTables(T const& collision, U const& jet, V& jetOutputTable, M& jetSubstructureOutputTable, std::vector<int> geoMatching, std::vector<int> ptMatching, std::vector<int> candMatching)
  {
    jetOutputTable(collision.globalIndex(), jet.globalIndex(), -1, geoMatching, ptMatching, candMatching, jet.pt(), jet.phi(), jet.eta(), jet.r(), jet.tracks().size());
    jetSubstructureOutputTable(jet.globalIndex(), jet.zg(), jet.rg(), jet.nsd());
  }

  void processDummy(aod::JCollision const& collision) {}
  PROCESS_SWITCH(JetSubstructureOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processOutputData(aod::JCollision const& collision,
                         JetTableData const& jets,
                         aod::JTracks const& tracks)
  {
    std::vector<int> geoMatching{-1};
    std::vector<int> ptMatching{-1};
    std::vector<int> candMatching{-1};
    for (const auto& jet : jets) {
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          fillTables(collision, jet, jetOutputTableData, jetSubstructureOutputTableData, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputData, "jet substructure output Data", false);

  void processOutputMCD(aod::JCollision const& collision,
                        JetTableMCD const& mcdjets,
                        JetTableMCP const& mcpjets,
                        aod::JTracks const& tracks)
  {
    std::vector<int> candMatching{-1};
    for (const auto& mcdjet : mcdjets) {
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (mcdjet.r() == round(jetRadiiValue * 100.0f)) {
          std::vector<int> geoMatching;
          std::vector<int> ptMatching;
          for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCP>()) {
            geoMatching.push_back(mcpjet.globalIndex());
          }
          for (auto& mcpjet : mcdjet.template matchedJetPt_as<JetTableMCP>()) {
            ptMatching.push_back(mcpjet.globalIndex());
          }
          fillTables(collision, mcdjet, jetOutputTableMCD, jetSubstructureOutputTableMCD, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMCD, "jet substructure output MCD", false);

  void processOutputMCP(aod::JMcCollision const& collision,
                        JetTableMCP const& mcpjets,
                        JetTableMCD const& mcdjets,
                        aod::JMcParticles const& particles)
  {
    std::vector<int> candMatching{-1};
    for (const auto& mcpjet : mcpjets) {
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (mcpjet.r() == round(jetRadiiValue * 100.0f)) {
          std::vector<int> geoMatching;
          std::vector<int> ptMatching;
          for (auto& mcdjet : mcpjet.template matchedJetGeo_as<JetTableMCD>()) {
            geoMatching.push_back(mcdjet.globalIndex());
          }
          for (auto& mcdjet : mcpjet.template matchedJetPt_as<JetTableMCD>()) {
            ptMatching.push_back(mcdjet.globalIndex());
          }
          fillTables(collision, mcpjet, jetOutputTableMCP, jetSubstructureOutputTableMCP, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMCP, "jet substructure output MCP", false);
};
using JetSubstructureOutput = JetSubstructureOutputTask<soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetSubstructures>>, aod::ChargedJetOutput, aod::ChargedJetSubstructureOutput, soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetSubstructures, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>, aod::ChargedMCDetectorLevelJetOutput, aod::ChargedMCDetectorLevelJetSubstructureOutput, soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetSubstructures, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>, aod::ChargedMCParticleLevelJetOutput, aod::ChargedMCParticleLevelJetSubstructureOutput>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutput>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-output"}));

  return WorkflowSpec{tasks};
}
