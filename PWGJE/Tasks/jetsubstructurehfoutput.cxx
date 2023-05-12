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

#include "PWGJE/DataModel/Jet.h"
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
  Configurable<int> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for B+"};

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processD0Data(soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0ChargedJetSubstructure>::iterator const& jet, // add template back
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

  void processD0MCD(soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetSubstructure>::iterator const& jet,
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

  void processD0MCP(soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetSubstructure>::iterator const& jet,
                    aod::McParticles const& particles)
  {

    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    auto hfparticle = hfparticles[0];
    auto Y = RecoDecay::y(array{hfparticle.px(), hfparticle.py(), hfparticle.pz()}, RecoDecay::getMassPDG(hfparticle.pdgCode()));
    auto M = RecoDecay::getMassPDG(hfparticle.pdgCode());
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), hfparticle.pt(), hfparticle.phi(), hfparticle.eta(), Y, M, M, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processD0MCP, "D0 jet substructure output on MC particle level", false);

  void processLcData(soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcChargedJetSubstructure>::iterator const& jet, // add template back
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

  void processLcMCD(soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetSubstructure>::iterator const& jet,
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

  void processLcMCP(soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetSubstructure>::iterator const& jet,
                    aod::McParticles const& particles)
  {

    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    auto hfparticle = hfparticles[0];
    auto Y = RecoDecay::y(array{hfparticle.px(), hfparticle.py(), hfparticle.pz()}, RecoDecay::getMassPDG(hfparticle.pdgCode()));
    auto M = RecoDecay::getMassPDG(hfparticle.pdgCode());
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), hfparticle.pt(), hfparticle.phi(), hfparticle.eta(), Y, M, M, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processLcMCP, "Lc jet substructure output on MC particle level", false);

  void processBplusData(soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusChargedJetSubstructure>::iterator const& jet, // add template back
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
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBplusData, "B+ jet substructure output on data", false);

  void processBplusMCD(soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelJetSubstructure>::iterator const& jet,
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
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBplusMCD, "B+ jet substructure output on MC detector level", false);

  void processBplusMCP(soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelJetSubstructure>::iterator const& jet,
                       aod::McParticles const& particles)
  {

    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    auto hfparticle = hfparticles[0];
    auto Y = RecoDecay::y(array{hfparticle.px(), hfparticle.py(), hfparticle.pz()}, RecoDecay::getMassPDG(hfparticle.pdgCode()));
    auto M = RecoDecay::getMassPDG(hfparticle.pdgCode());
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), hfparticle.pt(), hfparticle.phi(), hfparticle.eta(), Y, M, M, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBplusMCP, "B+ jet substructure output on MC particle level", false);
};
using JetSubstructureOutputDataD0 = JetSubstructureHFOutputTask<aod::D0ChargedJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelD0 = JetSubstructureHFOutputTask<aod::D0ChargedMCDetectorLevelJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelD0 = JetSubstructureHFOutputTask<aod::D0ChargedMCParticleLevelJetSubstructureOutput>;
using JetSubstructureOutputDataLc = JetSubstructureHFOutputTask<aod::LcChargedJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelLc = JetSubstructureHFOutputTask<aod::LcChargedMCDetectorLevelJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelLc = JetSubstructureHFOutputTask<aod::LcChargedMCParticleLevelJetSubstructureOutput>;
using JetSubstructureOutputDataBplus = JetSubstructureHFOutputTask<aod::BplusChargedJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelBplus = JetSubstructureHFOutputTask<aod::BplusChargedMCDetectorLevelJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelBplus = JetSubstructureHFOutputTask<aod::BplusChargedMCParticleLevelJetSubstructureOutput>;

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

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDataBplus>(cfgc,
                                                                       SetDefaultProcesses{},
                                                                       TaskName{"jet-substructure-output-Bplus-data"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCDetectorLevelBplus>(cfgc,
                                                                                  SetDefaultProcesses{},
                                                                                  TaskName{"jet-substructure-output-Bplus-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCParticleLevelBplus>(cfgc,
                                                                                  SetDefaultProcesses{},
                                                                                  TaskName{"jet-substructure-output-Bplus-mcp"}));

  return WorkflowSpec{tasks};
}
