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

// jet finder hf task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Jochen Klein <jochen.klein@cern.ch>

#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/Core/JetFindingUtilities.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/*
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> hfjetworkflows{{"d0-data-charged", VariantType::Int, 1, {"D0 jets charged data"}},
                                       {"d0-mcd-charged", VariantType::Int, 0, {"D0 jets charged MCD"}},
                                       {"d0-mcp-charged", VariantType::Int, 0, {"D0 jets charged MCD"}},
                                       {"bplus-data-charged", VariantType::Int, 0, {"B+ jets charged MCD"}},
                                       {"bplus-mcd-charged", VariantType::Int, 0, {"B+ jets charged MCD"}},
                                       {"bplus-mcp-charged", VariantType::Int, 0, {"B+ jets charged MCD"}},
                                       {"lc-data-charged", VariantType::Int, 0, {"Lc jets charged MCD"}},
                                       {"lc-mcd-charged", VariantType::Int, 0, {"Lc jets charged MCD"}},
                                       {"lc-mcp-charged", VariantType::Int, 0, {"Lc jets charged MCD"}}};
  std::swap(workflowOptions, hfjetworkflows);
}
*/

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename CandidateTableData, typename CandidateTableMCD, typename CandidateTableMCP, typename JetTracksSubTable, typename JetTable, typename ConstituentTable, typename JetEvtWiseSubTable, typename ConstituentEvtWiseSubTable>
struct JetFinderHFTask {
  Produces<JetTable> jetsTable;
  Produces<ConstituentTable> constituentsTable;
  Produces<JetEvtWiseSubTable> jetsEvtWiseSubTable;
  Produces<ConstituentEvtWiseSubTable> constituentsEvtWiseSubTable;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // cluster level configurables
  Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.71, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.71, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterPhiMin{"clusterPhiMin", 1.39, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 3.27, "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -25., "minimum Cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 25., "maximum Cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  // HF candidate level configurables
  Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};
  // HF candidiate selection configurables
  Configurable<bool> rejectBackgroundMCDCandidates{"rejectBackgroundMCDCandidates", false, "reject background HF candidates at MC detector level"};
  Configurable<bool> rejectIncorrectDecaysMCP{"rejectIncorrectDecaysMCP", true, "reject HF paticles decaying to the non-analysed decay channels at MC generator level"};

  // jet level configurables
  Configurable<std::vector<double>> jetRadius{"jetRadius", {0.4}, "jet resolution parameters"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<int> jetTypeParticleLevel{"jetTypeParticleLevel", 1, "Type of stored jets. 0 = full, 1 = charged, 2 = neutral"};
  Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
  Configurable<int> ghostRepeat{"ghostRepeat", 1, "set to 0 to gain speed if you dont need area calculation"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};

  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
  int eventSelection = -1;
  std::string particleSelection;

  JetFinder jetFinder;
  std::vector<fastjet::PseudoJet> inputParticles;

  double candMass;

  void init(InitContext const&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    jetFinder.etaMin = trackEtaMin;
    jetFinder.etaMax = trackEtaMax;
    jetFinder.jetPtMin = jetPtMin;
    jetFinder.jetPtMax = jetPtMax;
    jetFinder.jetEtaMin = jetEtaMin;
    jetFinder.jetEtaMax = jetEtaMax;
    if (jetEtaMin < -98.0) {
      jetFinder.jetEtaDefault = true;
    }
    jetFinder.algorithm = static_cast<fastjet::JetAlgorithm>(static_cast<int>(jetAlgorithm));
    jetFinder.recombScheme = static_cast<fastjet::RecombinationScheme>(static_cast<int>(jetRecombScheme));
    jetFinder.ghostArea = jetGhostArea;
    jetFinder.ghostRepeatN = ghostRepeat;

    candMass = jethfutilities::getTablePDGMass<CandidateTableData>();
  }

  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  Filter trackSubCuts = (aod::jtracksub::pt >= trackPtMin && aod::jtracksub::pt < trackPtMax && aod::jtracksub::eta > trackEtaMin && aod::jtracksub::eta < trackEtaMax && aod::jtracksub::phi >= trackPhiMin && aod::jtracksub::phi <= trackPhiMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));
  // Filter candidateCuts = (aod::hfcand::pt >= candPtMin && aod::hfcand::pt < candPtMax && aod::hfcand::y >= candYMin && aod::hfcand::y < candYMax);

  PresliceOptional<soa::Filtered<JetTracksSubTable>> perD0Candidate = aod::bkgd0::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perLcCandidate = aod::bkglc::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perBplusCandidate = aod::bkgbplus::candidateId;

  // function that generalically processes Data and reco level events
  template <bool isEvtWiseSub, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& collision, U const& tracks, V const& candidate, M& jetsTableInput, N& constituentsTableInput, O& originalTracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    inputParticles.clear();

    if constexpr (jethfutilities::isHFCandidate<V>()) {
      if (!jetfindingutilities::analyseCandidate(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax, candMass)) {
        return;
      }
    }

    if constexpr (jethfutilities::isHFMcCandidate<V>()) {
      if (!jetfindingutilities::analyseCandidateMC(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax, candMass, rejectBackgroundMCDCandidates)) {
        return;
      }
    }
    if constexpr (isEvtWiseSub) {
      jetfindingutilities::analyseTracks<U, typename U::iterator>(inputParticles, tracks, trackSelection);
    } else {
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, std::optional{candidate});
    }
    jetfindingutilities::findJets(jetFinder, inputParticles, jetRadius, jetAreaFractionMin, collision, jetsTableInput, constituentsTableInput, true);
  }

  // function that generalically processes gen level events
  template <typename T, typename U, typename V>
  void analyseMCP(T const& collision, U const& particles, V const& candidate, int jetTypeParticleLevel)
  {
    if (rejectIncorrectDecaysMCP && !jethfutilities::isMatchedHFCandidate(candidate)) { // is this even needed in the new derived format? it means any simulations run have to force the decay channel
      return;
    }

    inputParticles.clear();
    if (!jetfindingutilities::analyseCandidate(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax, candMass)) {
      return;
    }
    jetfindingutilities::analyseParticles(inputParticles, particleSelection, jetTypeParticleLevel, particles, pdgDatabase, std::optional{candidate});
    jetfindingutilities::findJets(jetFinder, inputParticles, jetRadius, jetAreaFractionMin, collision, jetsTable, constituentsTable, true);
  }

  void processDummy(JetCollisions const& collisions)
  {
  }
  PROCESS_SWITCH(JetFinderHFTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<JetTracks> const& tracks, CandidateTableData const& candidates)
  {
    for (typename CandidateTableData::iterator const& candidate : candidates) { // why can the type not be auto?  try const auto
      analyseCharged<false>(collision, tracks, candidate, jetsTable, constituentsTable, tracks);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsData, "charged hf jet finding on data", false);

  void processChargedEvtWiseSubJetsData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<JetTracksSubTable> const& tracks, CandidateTableData const& candidates)
  {
    for (typename CandidateTableData::iterator const& candidate : candidates) {
      analyseCharged<true>(collision, jethfutilities::slicedPerCandidate(tracks, candidate, perD0Candidate, perLcCandidate, perBplusCandidate), candidate, jetsEvtWiseSubTable, constituentsEvtWiseSubTable, tracks);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedEvtWiseSubJetsData, "charged hf jet finding on data with event-wise constituent subtraction", false);

  void processChargedJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<JetTracks> const& tracks, CandidateTableMCD const& candidates)
  {
    for (typename CandidateTableMCD::iterator const& candidate : candidates) {
      analyseCharged<false>(collision, tracks, candidate, jetsTable, constituentsTable, tracks);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsMCD, "charged hf jet finding on MC detector level", false);

  void processChargedEvtWiseSubJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<JetTracksSubTable> const& tracks, CandidateTableMCD const& candidates)
  {
    for (typename CandidateTableMCD::iterator const& candidate : candidates) {
      analyseCharged<true>(collision, jethfutilities::slicedPerCandidate(tracks, candidate, perD0Candidate, perLcCandidate, perBplusCandidate), candidate, jetsEvtWiseSubTable, constituentsEvtWiseSubTable, tracks);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedEvtWiseSubJetsMCD, "charged hf jet finding on MC detector level with event-wise constituent subtraction", false);

  void processChargedJetsMCP(JetMcCollision const& collision,
                             soa::Filtered<JetParticles> const& particles,
                             CandidateTableMCP const& candidates)
  {
    for (typename CandidateTableMCP::iterator const& candidate : candidates) {
      analyseMCP(collision, particles, candidate, 1);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsMCP, "hf jet finding on MC particle level", false);
};
/*

using JetFinderD0DataCharged = JetFinderHFTask<CandidatesD0Data,CandidatesD0MCD,CandidatesD0MCP,JetTracksSubD0, aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents>;
using JetFinderD0MCDetectorLevelCharged = JetFinderHFTask<CandidatesD0Data,CandidatesD0MCD,CandidatesD0MCP,JetTracksSubD0,aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
using JetFinderD0MCParticleLevelCharged = JetFinderHFTask<CandidatesD0Data,CandidatesD0MCD,CandidatesD0MCP,JetTracksSubD0,aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelEventWiseSubtractedJets, aod::D0ChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

using JetFinderBplusDataCharged = JetFinderHFTask<CandidatesBplusData,CandidatesBplusMCD,CandidatesBplusMCP,JetTracksSubBplus,aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents>;
using JetFinderBplusMCDetectorLevelCharged = JetFinderHFTask<CandidatesBplusData,CandidatesBplusMCD,CandidatesBplusMCP,JetTracksSubBplus,aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelEventWiseSubtractedJets, aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
using JetFinderBplusMCParticleLevelCharged = JetFinderHFTask<CandidatesBplusData,CandidatesBplusMCD,CandidatesBplusMCP,JetTracksSubBplus,aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelEventWiseSubtractedJets, aod::BplusChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

using JetFinderLcDataCharged = JetFinderHFTask<CandidatesLcData,CandidatesLcMCD,CandidatesLcMCP,JetTracksSubLc,aod::LcChargedJets, aod::LcChargedJetConstituents,aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents>;
using JetFinderLcMCDetectorLevelCharged = JetFinderHFTask<CandidatesLcData,CandidatesLcMCD,CandidatesLcMCP,JetTracksSubLc,aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelEventWiseSubtractedJets, aod::LcChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;
using JetFinderLcMCParticleLevelCharged = JetFinderHFTask<CandidatesLcData,CandidatesLcMCD,CandidatesLcMCP,JetTracksSubLc,aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelEventWiseSubtractedJets, aod::LcChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderD0DataCharged>(cfgc,
                                                               SetDefaultProcesses{},
                                                               TaskName{"jet-finder-d0-data-charged"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderD0MCDetectorLevelCharged>(cfgc,
                                                                          SetDefaultProcesses{},
                                                                          TaskName{"jet-finder-d0-mcd-charged"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderD0MCParticleLevelCharged>(cfgc,
                                                                          SetDefaultProcesses{},
                                                                          TaskName{"jet-finder-d0-mcp-charged"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderBplusDataCharged>(cfgc,
                                                                  SetDefaultProcesses{},
                                                                  TaskName{"jet-finder-bplus-data-charged"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderBplusMCDetectorLevelCharged>(cfgc,
                                                                             SetDefaultProcesses{},
                                                                             TaskName{"jet-finder-bplus-mcd-charged"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderBplusMCParticleLevelCharged>(cfgc,
                                                                             SetDefaultProcesses{},
                                                                             TaskName{"jet-finder-bplus-mcp-charged"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderLcDataCharged>(cfgc,
                                                               SetDefaultProcesses{},
                                                               TaskName{"jet-finder-lc-data-charged"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderLcMCDetectorLevelCharged>(cfgc,
                                                                          SetDefaultProcesses{},
                                                                          TaskName{"jet-finder-lc-mcd-charged"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderLcMCParticleLevelCharged>(cfgc,
                                                                          SetDefaultProcesses{},
                                                                          TaskName{"jet-finder-lc-mcp-charged"}));

  return WorkflowSpec{tasks};
}
*/
/*
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<JetFinderD0DataCharged>(cfgc)};
  workflow.clear();
  if (cfgc.hfjetworkflows().get<int>("d0-data-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderD0DataCharged>(cfgc));
  }
  if (cfgc.hfjetworkflows().get<int>("d0-mcd-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderD0MCDetectorLevelCharged>(cfgc));
  }
  if (cfgc.hfjetworkflows().get<int>("d0-mcp-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderD0MCParticleLevelCharged>(cfgc));
  }
  if (cfgc.hfjetworkflows().get<int>("bplus-data-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderBplusDataCharged>(cfgc));
  }
  if (cfgc.hfjetworkflows().get<int>("bplus-mcd-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderBplusMCDetectorLevelCharged>(cfgc));
  }
  if (cfgc.hfjetworkflows().get<int>("bplus-mcp-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderBplusMCParticleLevelCharged>(cfgc));
  }
  if (cfgc.hfjetworkflows().get<int>("lc-data-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderLcDataCharged>(cfgc));
  }
  if (cfgc.hfjetworkflows().get<int>("lc-mcd-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderLcMCDetectorLevelCharged>(cfgc));
  }
  if (cfgc.hfjetworkflows().get<int>("lc-mcp-charged")) {
    workflow.push_back(adaptAnalysisTask<JetFinderLcMCParticleLevelCharged>(cfgc));
  }
  return workflow;
}
*/
