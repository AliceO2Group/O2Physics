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
// Author: Nima Zardoshti
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

template <typename CollisionTable, typename JetTable, typename OutputTable, typename SubstructureOutputTable>
struct JetSubstructureOutputTask {
  Produces<OutputTable> jetOutputTable;
  Produces<SubstructureOutputTable> jetSubstructureOutputTable;

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;
  }

  Filter jetSelection = aod::jet::pt >= jetPtMin;

  template <typename T, typename U>
  void fillTables(T const& collision, U const& jets)
  {
    for (const auto& jet : jets) {
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          jetOutputTable(collision.globalIndex(), jet.globalIndex(), -1, jet.pt(), jet.phi(), jet.eta(), jet.r(), jet.tracks().size());
          jetSubstructureOutputTable(jet.globalIndex(), jet.zg(), jet.rg(), jet.nsd());
        }
      }
    }
  }

  void processDummy(typename CollisionTable::iterator const& collision) {}
  PROCESS_SWITCH(JetSubstructureOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processOutput(typename CollisionTable::iterator const& collision,
                     JetTable const& jets)
  {
    fillTables(collision, jets);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutput, "jet substructure output", false);
};
using JetSubstructureOutputData = JetSubstructureOutputTask<aod::Collisions, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetSubstructures>>, aod::ChargedJetOutput, aod::ChargedJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevel = JetSubstructureOutputTask<aod::Collisions, soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetSubstructures>>, aod::ChargedMCDetectorLevelJetOutput, aod::ChargedMCDetectorLevelJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevel = JetSubstructureOutputTask<aod::McCollisions, soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetSubstructures>>, aod::ChargedMCParticleLevelJetOutput, aod::ChargedMCParticleLevelJetSubstructureOutput>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputData>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-output-data"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCDetectorLevel>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-output-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCParticleLevel>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-output-mcp"}));

  return WorkflowSpec{tasks};
}
