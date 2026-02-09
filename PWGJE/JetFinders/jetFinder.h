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

#ifndef PWGJE_JETFINDERS_JETFINDER_H_
#define PWGJE_JETFINDERS_JETFINDER_H_

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h> // IWYU pragma: export

#include <THn.h>
#include <TMathBase.h>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <memory>
#include <string>
#include <vector>

template <typename JetTable, typename ConstituentTable, typename JetEvtWiseSubTable, typename ConstituentEvtWiseSubTable>
struct JetFinderTask {
  o2::framework::Produces<JetTable> jetsTable;
  o2::framework::Produces<ConstituentTable> constituentsTable;
  o2::framework::Produces<JetEvtWiseSubTable> jetsEvtWiseSubTable;
  o2::framework::Produces<ConstituentEvtWiseSubTable> constituentsEvtWiseSubTable;

  o2::framework::HistogramRegistry registry;

  // event level configurables
  o2::framework::Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  o2::framework::Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  o2::framework::Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  o2::framework::Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  o2::framework::Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  o2::framework::Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  o2::framework::Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};
  o2::framework::Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};

  // track level configurables
  o2::framework::Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  o2::framework::Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  o2::framework::Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  o2::framework::Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  o2::framework::Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  o2::framework::Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  o2::framework::Configurable<bool> applyTrackingEfficiency{"applyTrackingEfficiency", {false}, "configurable to decide whether to apply artificial tracking efficiency (discarding tracks) in jet finding"};
  o2::framework::Configurable<std::vector<double>> trackingEfficiencyPtBinning{"trackingEfficiencyPtBinning", {0., 10, 999.}, "pt binning of tracking efficiency array if applyTrackingEfficiency is true"};
  o2::framework::Configurable<std::vector<double>> trackingEfficiency{"trackingEfficiency", {1.0, 1.0}, "tracking efficiency array applied to jet finding if applyTrackingEfficiency is true"};
  o2::framework::Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  o2::framework::Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // cluster level configurables
  o2::framework::Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  o2::framework::Configurable<float> clusterEtaMin{"clusterEtaMin", -0.71, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  o2::framework::Configurable<float> clusterEtaMax{"clusterEtaMax", 0.71, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  o2::framework::Configurable<float> clusterPhiMin{"clusterPhiMin", 1.39, "minimum cluster phi"};
  o2::framework::Configurable<float> clusterPhiMax{"clusterPhiMax", 3.27, "maximum cluster phi"};
  o2::framework::Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  o2::framework::Configurable<float> clusterTimeMin{"clusterTimeMin", -25., "minimum Cluster time (ns)"};
  o2::framework::Configurable<float> clusterTimeMax{"clusterTimeMax", 25., "maximum Cluster time (ns)"};
  o2::framework::Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};
  o2::framework::Configurable<int> hadronicCorrectionType{"hadronicCorrectionType", 0, "0 = no correction, 1 = CorrectedOneTrack1, 2 = CorrectedOneTrack2, 3 = CorrectedAllTracks1, 4 = CorrectedAllTracks2"};
  o2::framework::Configurable<bool> doEMCALEventSelection{"doEMCALEventSelection", true, "apply the selection to the event alias_bit for full and neutral jets"};
  o2::framework::Configurable<bool> doEMCALEventSelectionChargedJets{"doEMCALEventSelectionChargedJets", false, "apply the selection to the event alias_bit for charged jets"};

  // jet level configurables
  o2::framework::Configurable<std::vector<double>> jetRadius{"jetRadius", {0.4}, "jet resolution parameters"};
  o2::framework::Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  o2::framework::Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  o2::framework::Configurable<float> jetEWSPtMin{"jetEWSPtMin", 0.0, "minimum event-wise subtracted jet pT"};
  o2::framework::Configurable<float> jetEWSPtMax{"jetEWSPtMax", 1000.0, "maximum event-wise subtracted jet pT"};
  o2::framework::Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  o2::framework::Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  o2::framework::Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  o2::framework::Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  o2::framework::Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
  o2::framework::Configurable<int> ghostRepeat{"ghostRepeat", 1, "set to 0 to gain speed if you dont need area calculation"};
  o2::framework::Configurable<bool> DoTriggering{"DoTriggering", false, "used for the charged jet trigger to remove the eta constraint on the jet axis"};
  o2::framework::Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  o2::framework::Configurable<int> jetPtBinWidth{"jetPtBinWidth", 5, "used to define the width of the jetPt bins for the THnSparse"};
  o2::framework::Configurable<bool> fillTHnSparse{"fillTHnSparse", false, "switch to fill the THnSparse"};
  o2::framework::Configurable<double> jetExtraParam{"jetExtraParam", -99.0, "sets the _extra_param in fastjet"};

  o2::framework::Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
  std::vector<int> eventSelectionBits;
  std::string particleSelection;

  JetFinder jetFinder;
  std::vector<fastjet::PseudoJet> inputParticles;

  std::vector<int> triggerMaskBits;

  o2::aod::EMCALClusterDefinition clusterDefinition;

  void init(o2::framework::InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    clusterDefinition = o2::aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);

    jetFinder.etaMin = trackEtaMin;
    jetFinder.etaMax = trackEtaMax;
    jetFinder.phiMin = trackPhiMin;
    jetFinder.phiMax = trackPhiMax;
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
    jetFinder.fastjetExtraParam = jetExtraParam;

    auto jetRadiiBins = (std::vector<double>)jetRadius;
    if (jetRadiiBins.size() > 1) {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (TMath::Abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
    } else {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
    }
    int jetPtMaxInt = static_cast<int>(jetPtMax);
    int jetPtMinInt = static_cast<int>(jetPtMin);
    jetPtMinInt = (jetPtMinInt / jetPtBinWidth) * jetPtBinWidth;
    jetPtMaxInt = ((jetPtMaxInt + jetPtBinWidth - 1) / jetPtBinWidth) * jetPtBinWidth;
    int jetPtBinNumber = (jetPtMaxInt - jetPtMinInt) / jetPtBinWidth;
    double jetPtMinDouble = static_cast<double>(jetPtMinInt);
    double jetPtMaxDouble = static_cast<double>(jetPtMaxInt);

    if (fillTHnSparse) {
      registry.add("hJet", "sparse for data or mcd jets", {o2::framework::HistType::kTHnD, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
      registry.add("hJetEWS", "sparse for data or mcd event-wise subtracted jets", {o2::framework::HistType::kTHnD, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
      registry.add("hJetMCP", "sparse for mcp jets", {o2::framework::HistType::kTHnD, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
      registry.add("hJetEWSMCP", "sparse for mcp event-wise subtracted jets", {o2::framework::HistType::kTHnD, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
    }

    if (applyTrackingEfficiency) {
      if (trackingEfficiencyPtBinning->size() < 2) {
        LOGP(fatal, "jetFinder workflow: trackingEfficiencyPtBinning configurable should have at least two bin edges");
      }
      if (trackingEfficiency->size() + 1 != trackingEfficiencyPtBinning->size()) {
        LOGP(fatal, "jetFinder workflow: trackingEfficiency configurable should have exactly one less entry than the number of bin edges set in trackingEfficiencyPtBinning configurable");
      }
    }
  }

  o2::framework::expressions::Filter collisionFilter = (nabs(o2::aod::jcollision::posZ) < vertexZCut && o2::aod::jcollision::centFT0M >= centralityMin && o2::aod::jcollision::centFT0M < centralityMax && o2::aod::jcollision::trackOccupancyInTimeRange <= trackOccupancyInTimeRangeMax); // should we add a posZ vtx cut here or leave it to analysers?
  o2::framework::expressions::Filter mcCollisionFilter = (nabs(o2::aod::jmccollision::posZ) < vertexZCut);
  o2::framework::expressions::Filter trackCuts = (o2::aod::jtrack::pt >= trackPtMin && o2::aod::jtrack::pt < trackPtMax && o2::aod::jtrack::eta >= trackEtaMin && o2::aod::jtrack::eta <= trackEtaMax && o2::aod::jtrack::phi >= trackPhiMin && o2::aod::jtrack::phi <= trackPhiMax); // do we need eta cut both here and in globalselection?
  o2::framework::expressions::Filter partCuts = (o2::aod::jmcparticle::pt >= trackPtMin && o2::aod::jmcparticle::pt < trackPtMax && o2::aod::jmcparticle::eta >= trackEtaMin && o2::aod::jmcparticle::eta <= trackEtaMax && o2::aod::jmcparticle::phi >= trackPhiMin && o2::aod::jmcparticle::phi <= trackPhiMax);
  o2::framework::expressions::Filter clusterFilter = (o2::aod::jcluster::definition == static_cast<int>(clusterDefinition) && o2::aod::jcluster::eta >= clusterEtaMin && o2::aod::jcluster::eta <= clusterEtaMax && o2::aod::jcluster::phi >= clusterPhiMin && o2::aod::jcluster::phi <= clusterPhiMax && o2::aod::jcluster::energy >= clusterEnergyMin && o2::aod::jcluster::time > clusterTimeMin && o2::aod::jcluster::time < clusterTimeMax && (!clusterRejectExotics || o2::aod::jcluster::isExotic != true));

  void processChargedJets(o2::soa::Filtered<o2::aod::JetCollisions>::iterator const& collision,
                          o2::soa::Filtered<o2::aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections) || !jetderiveddatautilities::selectTrigger(collision, triggerMaskBits) || (doEMCALEventSelectionChargedJets && !jetderiveddatautilities::eventEMCAL(collision))) {
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseTracks<o2::soa::Filtered<o2::aod::JetTracks>, o2::soa::Filtered<o2::aod::JetTracks>::iterator>(inputParticles, tracks, trackSelection, applyTrackingEfficiency, trackingEfficiency, trackingEfficiencyPtBinning);
    jetfindingutilities::findJets(jetFinder, inputParticles, jetPtMin, jetPtMax, jetRadius, jetAreaFractionMin, collision, jetsTable, constituentsTable, fillTHnSparse ? registry.get<THn>(HIST("hJet")) : std::shared_ptr<THn>(nullptr), fillTHnSparse);
  }

  PROCESS_SWITCH(JetFinderTask, processChargedJets, "Data and reco level jet finding for charged jets", false);

  void processChargedEvtWiseSubJets(o2::soa::Filtered<o2::aod::JetCollisions>::iterator const& collision,
                                    o2::soa::Filtered<o2::aod::JetTracksSub> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections) || !jetderiveddatautilities::selectTrigger(collision, triggerMaskBits) || (doEMCALEventSelectionChargedJets && !jetderiveddatautilities::eventEMCAL(collision))) {
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseTracks<o2::soa::Filtered<o2::aod::JetTracksSub>, o2::soa::Filtered<o2::aod::JetTracksSub>::iterator>(inputParticles, tracks, trackSelection, applyTrackingEfficiency, trackingEfficiency, trackingEfficiencyPtBinning);
    jetfindingutilities::findJets(jetFinder, inputParticles, jetEWSPtMin, jetEWSPtMax, jetRadius, jetAreaFractionMin, collision, jetsEvtWiseSubTable, constituentsEvtWiseSubTable, fillTHnSparse ? registry.get<THn>(HIST("hJetEWS")) : std::shared_ptr<THn>(nullptr), fillTHnSparse);
  }

  PROCESS_SWITCH(JetFinderTask, processChargedEvtWiseSubJets, "Data and reco level jet finding for charged jets with event-wise constituent subtraction", false);

  void processNeutralJets(o2::soa::Filtered<o2::aod::JetCollisions>::iterator const& collision,
                          o2::soa::Filtered<o2::aod::JetClusters> const& clusters)
  {
    if ((doEMCALEventSelection && !jetderiveddatautilities::eventEMCAL(collision)) || !jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections, "CBT_calo") || !jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseClusters(inputParticles, &clusters);
    jetfindingutilities::findJets(jetFinder, inputParticles, jetPtMin, jetPtMax, jetRadius, jetAreaFractionMin, collision, jetsTable, constituentsTable, fillTHnSparse ? registry.get<THn>(HIST("hJet")) : std::shared_ptr<THn>(nullptr), fillTHnSparse);
  }
  PROCESS_SWITCH(JetFinderTask, processNeutralJets, "Data and reco level jet finding for neutral jets", false);

  void processFullJets(o2::soa::Filtered<o2::aod::JetCollisions>::iterator const& collision,
                       o2::soa::Filtered<o2::aod::JetTracks> const& tracks,
                       o2::soa::Filtered<o2::aod::JetClusters> const& clusters)
  {
    if ((doEMCALEventSelection && !jetderiveddatautilities::eventEMCAL(collision)) || !jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections, "CBT_calo") || !jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseTracks<o2::soa::Filtered<o2::aod::JetTracks>, o2::soa::Filtered<o2::aod::JetTracks>::iterator>(inputParticles, tracks, trackSelection, applyTrackingEfficiency, trackingEfficiency, trackingEfficiencyPtBinning);
    jetfindingutilities::analyseClusters(inputParticles, &clusters, hadronicCorrectionType);
    jetfindingutilities::findJets(jetFinder, inputParticles, jetPtMin, jetPtMax, jetRadius, jetAreaFractionMin, collision, jetsTable, constituentsTable, fillTHnSparse ? registry.get<THn>(HIST("hJet")) : std::shared_ptr<THn>(nullptr), fillTHnSparse);
  }
  PROCESS_SWITCH(JetFinderTask, processFullJets, "Data and reco level jet finding for full and neutral jets", false);

  void processParticleLevelChargedJets(o2::soa::Filtered<o2::aod::JetMcCollisions>::iterator const& mcCollision, o2::soa::Filtered<o2::aod::JetParticles> const& particles)
  {
    if (!jetderiveddatautilities::selectMcCollision(mcCollision, skipMBGapEvents, applyRCTSelections)) {
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseParticles<false, o2::soa::Filtered<o2::aod::JetParticles>, o2::soa::Filtered<o2::aod::JetParticles>::iterator>(inputParticles, particleSelection, 1, particles, pdgDatabase);
    jetfindingutilities::findJets(jetFinder, inputParticles, jetPtMin, jetPtMax, jetRadius, jetAreaFractionMin, mcCollision, jetsTable, constituentsTable, fillTHnSparse ? registry.get<THn>(HIST("hJetMCP")) : std::shared_ptr<THn>(nullptr), fillTHnSparse);
  }
  PROCESS_SWITCH(JetFinderTask, processParticleLevelChargedJets, "Particle level charged jet finding", false);

  void processParticleLevelChargedEvtWiseSubJets(o2::soa::Filtered<o2::aod::JetMcCollisions>::iterator const& mcCollision, o2::soa::Filtered<o2::aod::JetParticlesSub> const& particles)
  {
    if (!jetderiveddatautilities::selectMcCollision(mcCollision, skipMBGapEvents, applyRCTSelections)) {
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseParticles<false, o2::soa::Filtered<o2::aod::JetParticlesSub>, o2::soa::Filtered<o2::aod::JetParticlesSub>::iterator>(inputParticles, particleSelection, 1, particles, pdgDatabase);
    jetfindingutilities::findJets(jetFinder, inputParticles, jetEWSPtMin, jetEWSPtMax, jetRadius, jetAreaFractionMin, mcCollision, jetsEvtWiseSubTable, constituentsEvtWiseSubTable, fillTHnSparse ? registry.get<THn>(HIST("hJetEWSMCP")) : std::shared_ptr<THn>(nullptr), fillTHnSparse);
  }
  PROCESS_SWITCH(JetFinderTask, processParticleLevelChargedEvtWiseSubJets, "Particle level charged with event-wise constituent subtraction jet finding", false);

  void processParticleLevelNeutralJets(o2::soa::Filtered<o2::aod::JetMcCollisions>::iterator const& mcCollision, o2::soa::Filtered<o2::aod::JetParticles> const& particles)
  {
    if (!jetderiveddatautilities::selectMcCollision(mcCollision, skipMBGapEvents, applyRCTSelections, "CBT_calo")) {
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseParticles<false, o2::soa::Filtered<o2::aod::JetParticles>, o2::soa::Filtered<o2::aod::JetParticles>::iterator>(inputParticles, particleSelection, 2, particles, pdgDatabase);
    jetfindingutilities::findJets(jetFinder, inputParticles, jetPtMin, jetPtMax, jetRadius, jetAreaFractionMin, mcCollision, jetsTable, constituentsTable, fillTHnSparse ? registry.get<THn>(HIST("hJetMCP")) : std::shared_ptr<THn>(nullptr), fillTHnSparse);
  }
  PROCESS_SWITCH(JetFinderTask, processParticleLevelNeutralJets, "Particle level neutral jet finding", false);

  void processParticleLevelFullJets(o2::soa::Filtered<o2::aod::JetMcCollisions>::iterator const& mcCollision, o2::soa::Filtered<o2::aod::JetParticles> const& particles)
  {
    if (!jetderiveddatautilities::selectMcCollision(mcCollision, skipMBGapEvents, applyRCTSelections, "CBT_calo")) {
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseParticles<false, o2::soa::Filtered<o2::aod::JetParticles>, o2::soa::Filtered<o2::aod::JetParticles>::iterator>(inputParticles, particleSelection, 0, particles, pdgDatabase);
    jetfindingutilities::findJets(jetFinder, inputParticles, jetPtMin, jetPtMax, jetRadius, jetAreaFractionMin, mcCollision, jetsTable, constituentsTable, fillTHnSparse ? registry.get<THn>(HIST("hJetMCP")) : std::shared_ptr<THn>(nullptr), fillTHnSparse);
  }

  PROCESS_SWITCH(JetFinderTask, processParticleLevelFullJets, "Particle level full jet finding", false);
};

#endif // PWGJE_JETFINDERS_JETFINDER_H_
