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

// jet finder task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL

#include "PWGJE/TableProducer/jetfinder.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTable, typename ConstituentTable, typename ConstituentSubTable>
struct JetFinderTask {
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
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -999., "minimum Cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 999., "maximum Cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  // jet level configurables
  Configurable<std::vector<double>> jetRadius{"jetRadius", {0.4}, "jet resolution parameters"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
  Configurable<int> ghostRepeat{"ghostRepeat", 1, "set to 0 to gain speed if you dont need area calculation"};
  Configurable<bool> DoTriggering{"DoTriggering", false, "used for the charged jet trigger to remove the eta constraint on the jet axis"};

  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
  int eventSelection = -1;
  std::string particleSelection;

  JetFinder jetFinder;
  std::vector<fastjet::PseudoJet> inputParticles;

  bool doConstSub = false;

  void init(InitContext const&)
  {
    trackSelection = JetDerivedDataUtilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelection = JetDerivedDataUtilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
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
    if (DoTriggering) {
      jetFinder.isTriggering = true;
    }
  }

  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax); // do we need eta cut both here and in globalselection?
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta > trackEtaMin && aod::jmcparticle::eta < trackEtaMax);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));

  void processDummy(aod::JCollision const& collision)
  {
  }

  PROCESS_SWITCH(JetFinderTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJets(soa::Filtered<aod::JCollisions>::iterator const& collision,
                          JetTracks const& tracks)
  {
    if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection)) {
      return;
    }
    inputParticles.clear();
    analyseTracks<JetTracks, JetTracks::iterator>(inputParticles, tracks, trackSelection);
    findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, doConstSub);
  }

  PROCESS_SWITCH(JetFinderTask, processChargedJets, "Data jet finding for charged jets", false);

  void processNeutralJets(soa::Filtered<aod::JCollisions>::iterator const& collision,
                          JetClusters const& clusters)
  {
    if (!JetDerivedDataUtilities::eventEMCAL(collision)) {
      return;
    }
    inputParticles.clear();
    analyseClusters(inputParticles, &clusters);
    findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, doConstSub);
  }
  PROCESS_SWITCH(JetFinderTask, processNeutralJets, "Data jet finding for neutral jets", false);

  void processFullJets(soa::Filtered<aod::JCollisions>::iterator const& collision,
                       JetTracks const& tracks,
                       JetClusters const& clusters)
  {
    if (!JetDerivedDataUtilities::eventEMCAL(collision)) {
      return;
    }
    inputParticles.clear();
    analyseTracks<JetTracks, JetTracks::iterator>(inputParticles, tracks, trackSelection);
    analyseClusters(inputParticles, &clusters);
    findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, doConstSub);
  }
  PROCESS_SWITCH(JetFinderTask, processFullJets, "Data jet finding for full and neutral jets", false);

  void processParticleLevelChargedJets(aod::JMcCollision const& collision, soa::Filtered<aod::JMcParticles> const& particles)
  {
    // TODO: MC event selection?
    analyseParticles<soa::Filtered<aod::JMcParticles>, soa::Filtered<aod::JMcParticles>::iterator>(inputParticles, particleSelection, 1, particles, pdgDatabase);
    findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, doConstSub);
  }
  PROCESS_SWITCH(JetFinderTask, processParticleLevelChargedJets, "Particle level charged jet finding", false);

  void processParticleLevelNeutralJets(aod::JMcCollision const& collision, soa::Filtered<aod::JMcParticles> const& particles)
  {
    // TODO: MC event selection?
    analyseParticles<soa::Filtered<aod::JMcParticles>, soa::Filtered<aod::JMcParticles>::iterator>(inputParticles, particleSelection, 2, particles, pdgDatabase);
    findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, doConstSub);
  }
  PROCESS_SWITCH(JetFinderTask, processParticleLevelNeutralJets, "Particle level neutral jet finding", false);

  void processParticleLevelFullJets(aod::JMcCollision const& collision, soa::Filtered<aod::JMcParticles> const& particles)
  {
    // TODO: MC event selection?
    analyseParticles<soa::Filtered<aod::JMcParticles>, soa::Filtered<aod::JMcParticles>::iterator>(inputParticles, particleSelection, 0, particles, pdgDatabase);
    findJets(jetFinder, inputParticles, jetRadius, collision, jetsTable, constituentsTable, constituentsSubTable, doConstSub);
  }

  PROCESS_SWITCH(JetFinderTask, processParticleLevelFullJets, "Particle level full jet finding", false);
};

using JetFinderDataCharged = JetFinderTask<o2::aod::ChargedJets, o2::aod::ChargedJetConstituents, o2::aod::ChargedJetConstituentsSub>;
using JetFinderDataFull = JetFinderTask<o2::aod::FullJets, o2::aod::FullJetConstituents, o2::aod::FullJetConstituentsSub>;
using JetFinderDataNeutral = JetFinderTask<o2::aod::NeutralJets, o2::aod::NeutralJetConstituents, o2::aod::NeutralJetConstituentsSub>;
using JetFinderMCDetectorLevelCharged = JetFinderTask<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents, o2::aod::ChargedMCDetectorLevelJetConstituentsSub>;
using JetFinderMCDetectorLevelFull = JetFinderTask<o2::aod::FullMCDetectorLevelJets, o2::aod::FullMCDetectorLevelJetConstituents, o2::aod::FullMCDetectorLevelJetConstituentsSub>;
using JetFinderMCDetectorLevelNeutral = JetFinderTask<o2::aod::NeutralMCDetectorLevelJets, o2::aod::NeutralMCDetectorLevelJetConstituents, o2::aod::NeutralMCDetectorLevelJetConstituentsSub>;
using JetFinderMCParticleLevelCharged = JetFinderTask<o2::aod::ChargedMCParticleLevelJets, o2::aod::ChargedMCParticleLevelJetConstituents, o2::aod::ChargedMCParticleLevelJetConstituentsSub>;
using JetFinderMCParticleLevelFull = JetFinderTask<o2::aod::FullMCParticleLevelJets, o2::aod::FullMCParticleLevelJetConstituents, o2::aod::FullMCParticleLevelJetConstituentsSub>;
using JetFinderMCParticleLevelNeutral = JetFinderTask<o2::aod::NeutralMCParticleLevelJets, o2::aod::NeutralMCParticleLevelJetConstituents, o2::aod::NeutralMCParticleLevelJetConstituentsSub>;
// using JetFinderHybridIntermediate = JetFinderTask<o2::aod::HybridIntermediateJets, o2::aod::HybridIntermediateJetConstituents, o2::aod::HybridIntermediateJetConstituentsSub>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderDataCharged>(cfgc,
                                            SetDefaultProcesses{}, TaskName{"jet-finder-data-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderDataFull>(cfgc,
                                         SetDefaultProcesses{}, TaskName{"jet-finder-data-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderDataNeutral>(cfgc,
                                            SetDefaultProcesses{}, TaskName{"jet-finder-data-neutral"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCDetectorLevelCharged>(cfgc,
                                                       SetDefaultProcesses{}, TaskName{"jet-finder-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCDetectorLevelFull>(cfgc,
                                                    SetDefaultProcesses{}, TaskName{"jet-finder-mcd-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCDetectorLevelNeutral>(cfgc,
                                                       SetDefaultProcesses{}, TaskName{"jet-finder-mcd-neutral"}));

  // tasks.emplace_back(
  // adaptAnalysisTask<JetFinderHybridIntermediate>(cfgc,
  // SetDefaultProcesses{}, TaskName{"jet-finder-hybrid-intermedaite-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCParticleLevelCharged>(cfgc,
                                                       SetDefaultProcesses{}, TaskName{"jet-finder-mcp-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCParticleLevelFull>(cfgc,
                                                    SetDefaultProcesses{}, TaskName{"jet-finder-mcp-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderMCParticleLevelNeutral>(cfgc,
                                                       SetDefaultProcesses{}, TaskName{"jet-finder-mcp-neutral"}));

  return WorkflowSpec{tasks};
}
