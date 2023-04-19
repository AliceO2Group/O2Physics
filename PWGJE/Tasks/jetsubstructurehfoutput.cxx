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
// Author: Nima Zardoshti
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

#include "PWGJE/DataModel/JetHF.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename OutputTable>
struct JetSubstructureHFOutputTask {
  Produces<OutputTable> jetSubstructurehfoutputTable;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for Lc->PKPi"};
  Configurable<int> selectionFlagLcToPiKP{"selectionFlagLcToPiKP", 1, "Selection Flag for Lc->PiKP"};
  Configurable<int> selectionFlagBPlus{"selectionFlagBPlus", 1, "Selection Flag for B+"};

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processD0Data(soa::Join<aod::D0Jets, aod::D0JetConstituents, aod::D0JetSubstructure>::iterator const& jet, // add template back
                     soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
                     aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>();
    auto cand = cands[0];

    auto invMassCand = -1;
    auto invMassCandBar = -1;
    if (cand.isSelD0() >= selectionFlagD0) {
      invMassCand = invMassD0ToPiK(cand);
    }
    if (cand.isSelD0bar() >= selectionFlagD0bar) {
      invMassCandBar = invMassD0barToKPi(cand);
    }

    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), cand.pt(), cand.phi(), cand.eta(), yD0(cand), invMassCand, invMassCandBar, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processD0Data, "D0 jet substructure output on data", false);

  void processD0MCD(soa::Join<aod::MCDetectorLevelD0Jets, aod::MCDetectorLevelD0JetConstituents, aod::MCDetectorLevelD0JetSubstructure>::iterator const& jet,
                    soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const& candidates,
                    aod::Tracks const& tracks)
  {

    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>();
    auto cand = cands[0];
    auto invMassCand = -1;
    auto invMassCandBar = -1;
    if (cand.isSelD0() >= selectionFlagD0) {
      invMassCand = invMassD0ToPiK(cand);
    }
    if (cand.isSelD0bar() >= selectionFlagD0bar) {
      invMassCandBar = invMassD0barToKPi(cand);
    }
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), cand.pt(), cand.phi(), cand.eta(), yD0(cand), invMassCand, invMassCandBar, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processD0MCD, "D0 jet substructure output on MC detector level", false);

  void processD0MCP(soa::Join<aod::MCParticleLevelD0Jets, aod::MCParticleLevelD0JetConstituents, aod::MCParticleLevelD0JetSubstructure>::iterator const& jet,
                    aod::McParticles const& particles)
  {

    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    auto hfparticle = hfparticles[0];
    auto Y = RecoDecay::y(array{hfparticle.px(), hfparticle.py(), hfparticle.pz()}, RecoDecay::getMassPDG(hfparticle.pdgCode()));
    auto M = RecoDecay::getMassPDG(hfparticle.pdgCode());
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), hfparticle.pt(), hfparticle.phi(), hfparticle.eta(), Y, M, M, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processD0MCP, "D0 jet substructure output on MC particle level", false);

  void processLcData(soa::Join<aod::LcJets, aod::LcJetConstituents, aod::LcJetSubstructure>::iterator const& jet, // add template back
                     soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidates,
                     aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
    auto cand = cands[0];

    auto invMassCand = -1;
    auto invMassCandBar = -1;
    if (cand.isSelLcToPKPi() >= selectionFlagLcToPKPi) {
      invMassCand = invMassLcToPKPi(cand);
    }
    if (cand.isSelLcToPiKP() >= selectionFlagLcToPiKP) {
      invMassCandBar = invMassLcToPiKP(cand);
    }

    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), cand.pt(), cand.phi(), cand.eta(), yLc(cand), invMassCand, invMassCandBar, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processLcData, "Lc jet substructure output on data", false);

  void processLcMCD(soa::Join<aod::MCDetectorLevelLcJets, aod::MCDetectorLevelLcJetConstituents, aod::MCDetectorLevelLcJetSubstructure>::iterator const& jet,
                    soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const& candidates,
                    aod::Tracks const& tracks)
  {

    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>();
    auto cand = cands[0];
    auto invMassCand = -1;
    auto invMassCandBar = -1;
    if (cand.isSelLcToPKPi() >= selectionFlagLcToPKPi) {
      invMassCand = invMassLcToPKPi(cand);
    }
    if (cand.isSelLcToPiKP() >= selectionFlagLcToPiKP) {
      invMassCandBar = invMassLcToPiKP(cand);
    }
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), cand.pt(), cand.phi(), cand.eta(), yLc(cand), invMassCand, invMassCandBar, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processLcMCD, "Lc jet substructure output on MC detector level", false);

  void processLcMCP(soa::Join<aod::MCParticleLevelLcJets, aod::MCParticleLevelLcJetConstituents, aod::MCParticleLevelLcJetSubstructure>::iterator const& jet,
                    aod::McParticles const& particles)
  {

    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    auto hfparticle = hfparticles[0];
    auto Y = RecoDecay::y(array{hfparticle.px(), hfparticle.py(), hfparticle.pz()}, RecoDecay::getMassPDG(hfparticle.pdgCode()));
    auto M = RecoDecay::getMassPDG(hfparticle.pdgCode());
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), hfparticle.pt(), hfparticle.phi(), hfparticle.eta(), Y, M, M, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processLcMCP, "Lc jet substructure output on MC particle level", false);

  void processBPlusData(soa::Join<aod::BPlusJets, aod::BPlusJetConstituents, aod::BPlusJetSubstructure>::iterator const& jet, // add template back
                        soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi> const& candidates,
                        aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>>();
    auto cand = cands[0];

    auto invMassCand = -1;
    auto invMassCandBar = -1;
    invMassCand = invMassBplusToD0Pi(cand);
    invMassCandBar = invMassCand;

    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), cand.pt(), cand.phi(), cand.eta(), yBplus(cand), invMassCand, invMassCandBar, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBPlusData, "B+ jet substructure output on data", false);

  void processBPlusMCD(soa::Join<aod::MCDetectorLevelBPlusJets, aod::MCDetectorLevelBPlusJetConstituents, aod::MCDetectorLevelBPlusJetSubstructure>::iterator const& jet,
                       soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec> const& candidates,
                       aod::Tracks const& tracks)
  {

    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>>();
    auto cand = cands[0];
    auto invMassCand = -1;
    auto invMassCandBar = -1;
    invMassCand = invMassBplusToD0Pi(cand);
    invMassCandBar = invMassCand;
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), cand.pt(), cand.phi(), cand.eta(), yBplus(cand), invMassCand, invMassCandBar, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBPlusMCD, "B+ jet substructure output on MC detector level", false);

  void processBPlusMCP(soa::Join<aod::MCParticleLevelBPlusJets, aod::MCParticleLevelBPlusJetConstituents, aod::MCParticleLevelBPlusJetSubstructure>::iterator const& jet,
                       aod::McParticles const& particles)
  {

    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    auto hfparticle = hfparticles[0];
    auto Y = RecoDecay::y(array{hfparticle.px(), hfparticle.py(), hfparticle.pz()}, RecoDecay::getMassPDG(hfparticle.pdgCode()));
    auto M = RecoDecay::getMassPDG(hfparticle.pdgCode());
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), hfparticle.pt(), hfparticle.phi(), hfparticle.eta(), Y, M, M, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBPlusMCP, "B+ jet substructure output on MC particle level", false);
};
using JetSubstructureOutputDataD0 = JetSubstructureHFOutputTask<aod::D0JetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelD0 = JetSubstructureHFOutputTask<aod::MCDetectorLevelD0JetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelD0 = JetSubstructureHFOutputTask<aod::MCParticleLevelD0JetSubstructureOutput>;
using JetSubstructureOutputDataLc = JetSubstructureHFOutputTask<aod::LcJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelLc = JetSubstructureHFOutputTask<aod::MCDetectorLevelLcJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelLc = JetSubstructureHFOutputTask<aod::MCParticleLevelLcJetSubstructureOutput>;
using JetSubstructureOutputDataBPlus = JetSubstructureHFOutputTask<aod::BPlusJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelBPlus = JetSubstructureHFOutputTask<aod::MCDetectorLevelBPlusJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelBPlus = JetSubstructureHFOutputTask<aod::MCParticleLevelBPlusJetSubstructureOutput>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDataD0>(cfgc,
                                                                    SetDefaultProcesses{},
                                                                    TaskName{"jet-substructure-output-D0-data"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCDetectorLevelD0>(cfgc,
                                                                               SetDefaultProcesses{},
                                                                               TaskName{"jet-substructure-output-D0-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCParticleLevelD0>(cfgc,
                                                                               SetDefaultProcesses{},
                                                                               TaskName{"jet-substructure-output-D0-mcp"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDataLc>(cfgc,
                                                                    SetDefaultProcesses{},
                                                                    TaskName{"jet-substructure-output-Lc-data"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCDetectorLevelLc>(cfgc,
                                                                               SetDefaultProcesses{},
                                                                               TaskName{"jet-substructure-output-Lc-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCParticleLevelLc>(cfgc,
                                                                               SetDefaultProcesses{},
                                                                               TaskName{"jet-substructure-output-Lc-mcp"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDataBPlus>(cfgc,
                                                                       SetDefaultProcesses{},
                                                                       TaskName{"jet-substructure-output-BPlus-data"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCDetectorLevelBPlus>(cfgc,
                                                                                  SetDefaultProcesses{},
                                                                                  TaskName{"jet-substructure-output-BPlus-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCParticleLevelBPlus>(cfgc,
                                                                                  SetDefaultProcesses{},
                                                                                  TaskName{"jet-substructure-output-BPlus-mcp"}));

  return WorkflowSpec{tasks};
}
