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
// Authors: Nima Zardoshti, Jochen Klein

#include "PWGJE/TableProducer/jetfinder.h"

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

template <typename CandidateTableData, typename CandidateTableMCD, typename ParticleTable, typename JetTable, typename ConstituentTable, typename ConstituentSubTable>
struct JetFinderHFTask {
  Produces<JetTable> jetsTable;
  Produces<ConstituentTable> constituentsTable;
  Produces<ConstituentSubTable> constituentsSubTable;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

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
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterPhiMin{"clusterPhiMin", -999, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 999, "maximum cluster phi"};

  // HF candidate level configurables
  Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};
  // HF candidiate selection configurables
  Configurable<bool> rejectBackgroundMCCandidates{"rejectBackgroundMCCandidates", true, "reject background HF candidates at MC detector level"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for Lc->PKPi"};
  Configurable<int> selectionFlagLcToPiPK{"selectionFlagLcToPiPK", 1, "Selection Flag for Lc->PiPK"};
  Configurable<int> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for B+"};

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
  Configurable<bool> DoRhoAreaSub{"DoRhoAreaSub", false, "do rho area subtraction"};
  Configurable<bool> DoConstSub{"DoConstSub", false, "do constituent subtraction"};

  Service<o2::framework::O2DatabasePDG> pdg;
  std::string trackSelection;
  std::string eventSelection;
  std::string particleSelection;

  JetFinder jetFinder;
  std::vector<fastjet::PseudoJet> inputParticles;

  int candPDG;
  int candDecay;

  void init(InitContext const&)
  {
    trackSelection = static_cast<std::string>(trackSelections);
    eventSelection = static_cast<std::string>(eventSelections);
    particleSelection = static_cast<std::string>(particleSelections);

    if (DoRhoAreaSub) {
      jetFinder.setBkgSubMode(JetFinder::BkgSubMode::rhoAreaSub);
    }
    if (DoConstSub) {
      jetFinder.setBkgSubMode(JetFinder::BkgSubMode::constSub);
    }

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

    if constexpr (std::is_same_v<std::decay_t<CandidateTableData>, CandidatesD0Data>) { // Note : need to be careful if configurable workflow options are added later
      candPDG = static_cast<int>(pdg::Code::kD0);
      candDecay = static_cast<int>(aod::hf_cand_2prong::DecayType::D0ToPiK);
    }
    if constexpr (std::is_same_v<std::decay_t<CandidateTableData>, CandidatesBplusData>) {
      candPDG = static_cast<int>(pdg::Code::kBPlus);
      candDecay = static_cast<int>(aod::hf_cand_bplus::DecayType::BplusToD0Pi);
    }
    if constexpr (std::is_same_v<std::decay_t<CandidateTableData>, CandidatesLcData>) {
      candPDG = static_cast<int>(pdg::Code::kLambdaCPlus);
      candDecay = static_cast<int>(aod::hf_cand_3prong::DecayType::LcToPKPi);
    }
  }

  o2::aod::EMCALClusterDefinition clusterDefinition = o2::aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter collisionFilter = (nabs(aod::collision::posZ) < vertexZCut);
  Filter trackCuts = (aod::track::pt >= trackPtMin && aod::track::pt < trackPtMax && aod::track::eta > trackEtaMin && aod::track::eta < trackEtaMax && aod::track::phi >= trackPhiMin && aod::track::phi <= trackPhiMax);
  Filter partCuts = (aod::mcparticle::pt >= trackPtMin && aod::mcparticle::pt < trackPtMax);
  Filter clusterFilter = (o2::aod::emcalcluster::definition == static_cast<int>(clusterDefinition) && aod::emcalcluster::eta > clusterEtaMin && aod::emcalcluster::eta < clusterEtaMax && aod::emcalcluster::phi >= clusterPhiMin && aod::emcalcluster::phi <= clusterPhiMax);
  // Filter candidateCuts = (aod::hfcand::pt >= candPtMin && aod::hfcand::pt < candPtMax && aod::hfcand::y >= candYMin && aod::hfcand::y < candYMax);
  Filter candidateCutsD0 = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);
  Filter candidateCutsLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLcToPKPi || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLcToPiPK);
  Filter candidateCutsBplus = (aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBplus);

  // function that processes data for all 2Prong candidates
  template <typename T, typename U, typename M>
  void analyseData(T const& collision, U const& tracks, M const& candidates)
  {
    if (!selectCollision(collision, eventSelection)) {
      return;
    }

    for (auto& candidate : candidates) {
      inputParticles.clear();
      if (!analyseCandidate(inputParticles, candPDG, candPtMin, candPtMax, candYMin, candYMax, candidate)) {
        continue;
      }
      analyseTracks(inputParticles, tracks, trackSelection, std::optional{candidate});
      findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, DoConstSub, true);
    }
  }

  // function that processes MC det for all 2Prong candidates
  template <typename T, typename U, typename M>
  void analyseMCD(T const& collision, U const& tracks, M const& candidates)
  {
    if (!selectCollision(collision, eventSelection)) {
      return;
    }

    for (auto& candidate : candidates) {
      inputParticles.clear();
      if (!analyseCandidateMC(inputParticles, candPDG, candDecay, candPtMin, candPtMax, candYMin, candYMax, candidate, rejectBackgroundMCCandidates)) {
        continue;
      }
      analyseTracks(inputParticles, tracks, trackSelection, std::optional{candidate});
      findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, DoConstSub, true);
    }
  }

  // function that generalically processes gen level events
  template <typename T, typename U>
  void analyseMCP(T const& collision, U const& particles)
  {
    inputParticles.clear();
    std::vector<typename U::iterator> candidates;
    candidates.clear();

    for (auto const& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) & (1 << candDecay)) {
        auto particleY = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        if (particleY < candYMin || particleY > candYMax) {
          continue;
        }
        if (particle.pt() < candPtMin || particle.pt() >= candPtMax) {
          continue;
        }
        candidates.push_back(particle);
      }
    }
    for (auto& candidate : candidates) {
      analyseParticles(inputParticles, particleSelection, jetTypeParticleLevel, particles, pdg->Instance(), std::optional{candidate});
      FastJetUtilities::fillTracks(candidate, inputParticles, candidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candidate.pdgCode()));
      findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, DoConstSub, true);
    }
  }

  void processDummy(aod::Collisions const& collision)
  {
  }
  PROCESS_SWITCH(JetFinderHFTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, JetTracks const& tracks, CandidateTableData const& candidates) { analyseData(collision, tracks, candidates); }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsData, "charged hf jet finding on data", false);

  void processChargedJetsMCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, JetTracks const& tracks, CandidateTableMCD const& candidates) { analyseMCD(collision, tracks, candidates); }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsMCD, "charged hf jet finding on MC detector level", false);

  void processChargedJetsMCP(aod::McCollision const& collision,
                             ParticleTable const& particles)
  {
    analyseMCP(collision, particles);
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsMCP, "hf jet finding on MC particle level", false);
};
/*

using JetFinderD0DataCharged = JetFinderHFTask<CandidatesD0Data,CandidatesD0MCD,ParticlesD0,o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents, o2::aod::D0ChargedJetConstituentsSub>;
using JetFinderD0MCDetectorLevelCharged = JetFinderHFTask<CandidatesD0Data,CandidatesD0MCD,ParticlesD0,o2::aod::D0ChargedMCDetectorLevelJets, o2::aod::D0ChargedMCDetectorLevelJetConstituents, o2::aod::D0ChargedMCDetectorLevelJetConstituentsSub>;
using JetFinderD0MCParticleLevelCharged = JetFinderHFTask<CandidatesD0Data,CandidatesD0MCD,ParticlesD0,o2::aod::D0ChargedMCParticleLevelJets, o2::aod::D0ChargedMCParticleLevelJetConstituents, o2::aod::D0ChargedMCParticleLevelJetConstituentsSub>;

using JetFinderBplusDataCharged = JetFinderHFTask<CandidatesBplusData,CandidatesBplusMCD,ParticlesBplus,o2::aod::BplusChargedJets, o2::aod::BplusChargedJetConstituents, o2::aod::BplusChargedJetConstituentsSub>;
using JetFinderBplusMCDetectorLevelCharged = JetFinderHFTask<CandidatesBplusData,CandidatesBplusMCD,ParticlesBplus,o2::aod::BplusChargedMCDetectorLevelJets, o2::aod::BplusChargedMCDetectorLevelJetConstituents, o2::aod::BplusChargedMCDetectorLevelJetConstituentsSub>;
using JetFinderBplusMCParticleLevelCharged = JetFinderHFTask<CandidatesBplusData,CandidatesBplusMCD,ParticlesBplus,o2::aod::BplusChargedMCParticleLevelJets, o2::aod::BplusChargedMCParticleLevelJetConstituents, o2::aod::BplusChargedMCParticleLevelJetConstituentsSub>;

using JetFinderLcDataCharged = JetFinderHFTask<CandidatesLcData,CandidatesLcMCD,ParticlesLc,o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents, o2::aod::LcChargedJetConstituentsSub>;
using JetFinderLcMCDetectorLevelCharged = JetFinderHFTask<CandidatesLcData,CandidatesLcMCD,ParticlesLc,o2::aod::LcChargedMCDetectorLevelJets, o2::aod::LcChargedMCDetectorLevelJetConstituents, o2::aod::LcChargedMCDetectorLevelJetConstituentsSub>;
using JetFinderLcMCParticleLevelCharged = JetFinderHFTask<CandidatesLcData,CandidatesLcMCD,ParticlesLc,o2::aod::LcChargedMCParticleLevelJets, o2::aod::LcChargedMCParticleLevelJetConstituents, o2::aod::LcChargedMCParticleLevelJetConstituentsSub>;

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
