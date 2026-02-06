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

// Task to produce a table joinable to the jcollision or candidate table with the mean background pT density
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <string>
#include <type_traits>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct RhoEstimatorTask {
  Produces<aod::BkgChargedRhos> rhoChargedTable;
  Produces<aod::BkgChargedMcRhos> rhoChargedMcTable;
  Produces<aod::BkgD0Rhos> rhoD0Table;
  Produces<aod::BkgD0McRhos> rhoD0McTable;
  Produces<aod::BkgDplusRhos> rhoDplusTable;
  Produces<aod::BkgDplusMcRhos> rhoDplusMcTable;
  Produces<aod::BkgDsRhos> rhoDsTable;
  Produces<aod::BkgDsMcRhos> rhoDsMcTable;
  Produces<aod::BkgDstarRhos> rhoDstarTable;
  Produces<aod::BkgDstarMcRhos> rhoDstarMcTable;
  Produces<aod::BkgLcRhos> rhoLcTable;
  Produces<aod::BkgLcMcRhos> rhoLcMcTable;
  Produces<aod::BkgB0Rhos> rhoB0Table;
  Produces<aod::BkgB0McRhos> rhoB0McTable;
  Produces<aod::BkgBplusRhos> rhoBplusTable;
  Produces<aod::BkgBplusMcRhos> rhoBplusMcTable;
  Produces<aod::BkgXicToXiPiPiRhos> rhoXicToXiPiPiTable;
  Produces<aod::BkgXicToXiPiPiMcRhos> rhoXicToXiPiPiMcTable;
  Produces<aod::BkgDielectronRhos> rhoDielectronTable;
  Produces<aod::BkgDielectronMcRhos> rhoDielectronMcTable;

  struct : ConfigurableGroup {

    Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
    Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
    Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
    Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
    Configurable<float> trackPhiMin{"trackPhiMin", -99.0, "minimum track phi"};
    Configurable<float> trackPhiMax{"trackPhiMax", 99.0, "maximum track phi"};
    Configurable<bool> applyTrackingEfficiency{"applyTrackingEfficiency", {false}, "configurable to decide whether to apply artificial tracking efficiency (discarding tracks) in the collision analysed by this task"};
    Configurable<std::vector<double>> trackingEfficiencyPtBinning{"trackingEfficiencyPtBinning", {0., 10, 999.}, "pt binning of tracking efficiency array if applyTrackingEfficiency is true"};
    Configurable<std::vector<double>> trackingEfficiency{"trackingEfficiency", {1.0, 1.0}, "tracking efficiency array applied if applyTrackingEfficiency is true"};
    Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

    Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

    Configurable<int> jetAlgorithm{"jetAlgorithm", 0, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
    Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
    Configurable<float> bkgjetR{"bkgjetR", 0.2, "jet resolution parameter for determining background density"};
    Configurable<float> bkgEtaMin{"bkgEtaMin", -0.7, "minimim pseudorapidity for determining background density"};
    Configurable<float> bkgEtaMax{"bkgEtaMax", 0.7, "maximum pseudorapidity for determining background density"};
    Configurable<float> bkgPhiMin{"bkgPhiMin", -6.283, "minimim phi for determining background density"};
    Configurable<float> bkgPhiMax{"bkgPhiMax", 6.283, "maximum phi for determining background density"};
    Configurable<bool> doSparse{"doSparse", false, "perfom sparse estimation"};
    Configurable<double> ghostRapMax{"ghostRapMax", 0.9, "Ghost rapidity max"};
    Configurable<int> ghostRepeat{"ghostRepeat", 1, "Ghost tiling repeats"};
    Configurable<double> ghostArea{"ghostArea", 0.005, "Area per ghost"};
    Configurable<double> ghostGridScatter{"ghostGridScatter", 1.0, "Grid scatter"};
    Configurable<double> ghostKtScatter{"ghostKtScatter", 0.1, "kT scatter"};
    Configurable<double> ghostMeanPt{"ghostMeanPt", 1e-100, "Mean ghost pT"};

    Configurable<float> thresholdTriggerTrackPtMin{"thresholdTriggerTrackPtMin", 0.0, "Minimum trigger track pt to accept event"};
    Configurable<float> thresholdClusterEnergyMin{"thresholdClusterEnergyMin", 0.0, "Minimum cluster energy to accept event"};
    Configurable<bool> performTriggerTrackSelection{"performTriggerTrackSelection", false, "only accept trigger tracks that pass one of the track selections"};
    Configurable<float> triggerTrackPtSelectionMin{"triggerTrackPtSelectionMin", 0.15, "only accept trigger tracks that have a pT larger than this pT"};
    Configurable<float> triggerTrackEtaSelectionMax{"triggerTrackEtaSelectionMax", 0.9, "only accept trigger tracks that have an eta smaller than this eta"};

    Configurable<float> vertexZCut{"vertexZCut", 10.0, "z-vertex cut on event"};
    Configurable<std::string> eventSelections{"eventSelections", "", "choose event selection"};
    Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
    Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
    Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
    Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};
    Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};

    Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  } config;

  JetBkgSubUtils bkgSub;
  float bkgPhiMax_;
  float bkgPhiMin_;
  std::vector<fastjet::PseudoJet> inputParticles;
  int trackSelection = -1;
  std::string particleSelection;

  std::vector<bool> collisionFlag;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  std::vector<int> eventSelectionBits;
  std::vector<int> triggerMaskBits;
  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(config.trackSelections));
    particleSelection = static_cast<std::string>(config.particleSelections);

    bkgSub.setJetAlgorithmAndScheme(static_cast<fastjet::JetAlgorithm>(static_cast<int>(config.jetAlgorithm)), static_cast<fastjet::RecombinationScheme>(static_cast<int>(config.jetRecombScheme)));
    bkgSub.setJetBkgR(config.bkgjetR);
    bkgSub.setEtaMinMax(config.bkgEtaMin, config.bkgEtaMax);
    bkgPhiMax_ = config.bkgPhiMax;
    bkgPhiMin_ = config.bkgPhiMin;
    if (config.bkgPhiMax > 98.0) {
      bkgPhiMax_ = 2.0 * M_PI;
    }
    if (config.bkgPhiMin < -98.0) {
      bkgPhiMin_ = -2.0 * M_PI;
    }
    bkgSub.setPhiMinMax(bkgPhiMin_, bkgPhiMax_);

    fastjet::GhostedAreaSpec ghostAreaSpec(config.ghostRapMax, config.ghostRepeat, config.ghostArea,
                                           config.ghostGridScatter, config.ghostKtScatter, config.ghostMeanPt);
    bkgSub.setGhostAreaSpec(ghostAreaSpec);

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(config.eventSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(config.triggerMasks);

    if (config.applyTrackingEfficiency) {
      if (config.trackingEfficiencyPtBinning->size() < 2) {
        LOGP(fatal, "rhoEstimator workflow: trackingEfficiencyPtBinning configurable should have at least two bin edges");
      }
      if (config.trackingEfficiency->size() + 1 != config.trackingEfficiencyPtBinning->size()) {
        LOGP(fatal, "rhoEstimator workflow: trackingEfficiency configurable should have exactly one less entry than the number of bin edges set in trackingEfficiencyPtBinning configurable");
      }
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= config.trackPtMin && aod::jtrack::pt < config.trackPtMax && aod::jtrack::eta > config.trackEtaMin && aod::jtrack::eta < config.trackEtaMax && aod::jtrack::phi >= config.trackPhiMin && aod::jtrack::phi <= config.trackPhiMax);
  Filter partCuts = (aod::jmcparticle::pt >= config.trackPtMin && aod::jmcparticle::pt < config.trackPtMax && aod::jmcparticle::eta >= config.trackEtaMin && aod::jmcparticle::eta <= config.trackEtaMax && aod::jmcparticle::phi >= config.trackPhiMin && aod::jmcparticle::phi <= config.trackPhiMax);

  void processSetupCollisionSelection(aod::JCollisions const& collisions)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size(), false);
  }

  void processSetupEventTriggering(aod::JCollision const& collision)
  {
    if (jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      collisionFlag[collision.globalIndex()] = true;
    }
  }

  template <typename T>
  void processSelectionObjects(T& selectionObjects)
  {
    float selectionObjectPtMin = 0.0;
    if constexpr (std::is_same_v<std::decay_t<T>, aod::JTracks>) {
      selectionObjectPtMin = config.thresholdTriggerTrackPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::JClusters>) {
      selectionObjectPtMin = config.thresholdClusterEnergyMin;
    } else {
      selectionObjectPtMin = 0.0;
    }
    for (const auto& selectionObject : selectionObjects) {
      bool isTriggerObject = false;
      if constexpr (std::is_same_v<std::decay_t<T>, aod::JClusters>) {
        if (selectionObject.energy() >= selectionObjectPtMin) {
          isTriggerObject = true;
        }
      } else {
        if constexpr (std::is_same_v<std::decay_t<T>, aod::JTracks>) {
          if (config.performTriggerTrackSelection && !(selectionObject.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) {
            continue;
          }
          if (selectionObject.pt() < config.triggerTrackPtSelectionMin || std::abs(selectionObject.eta()) > config.triggerTrackEtaSelectionMax) {
            continue;
          }
        }
        if (selectionObject.pt() >= selectionObjectPtMin) {
          isTriggerObject = true;
        }
      }
      if (isTriggerObject) {
        if (selectionObject.collisionId() >= 0) {
          collisionFlag[selectionObject.collisionId()] = true;
        }
      }
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processSetupCollisionSelection, "setup the writing for data based on collisions", false);
  PROCESS_SWITCH(RhoEstimatorTask, processSetupEventTriggering, "process software triggers", false);
  PROCESS_SWITCH_FULL(RhoEstimatorTask, processSelectionObjects<aod::JClusters>, processSelectingClusters, "process EMCal clusters", false);
  PROCESS_SWITCH_FULL(RhoEstimatorTask, processSelectionObjects<aod::JTracks>, processSelectingTracks, "process high pt tracks", false);

  void processChargedCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
      rhoChargedTable(0.0, 0.0);
      return;
    }
    if (collisionFlag.size() != 0 && !collisionFlag[collision.globalIndex()]) {
      rhoChargedTable(0.0, 0.0);
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseTracks<soa::Filtered<aod::JetTracks>, soa::Filtered<aod::JetTracks>::iterator>(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning);
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
    rhoChargedTable(rho, rhoM);
  }
  PROCESS_SWITCH(RhoEstimatorTask, processChargedCollisions, "Fill rho tables for collisions using charged tracks", true);

  void processChargedMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles)
  {
    if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
      rhoChargedMcTable(0.0, 0.0);
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseParticles<false, soa::Filtered<aod::JetParticles>, soa::Filtered<aod::JetParticles>::iterator>(inputParticles, particleSelection, 1, particles, pdgDatabase);
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
    rhoChargedMcTable(rho, rhoM);
  }
  PROCESS_SWITCH(RhoEstimatorTask, processChargedMcCollisions, "Fill rho tables for MC collisions using charged tracks", false);

  void processD0Collisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesD0Data const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoD0Table(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoD0Table(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processD0Collisions, "Fill rho tables for collisions with D0 candidates", false);

  void processD0McCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesD0MCP const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoD0McTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoD0McTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processD0McCollisions, "Fill rho tables for collisions with D0 MCP candidates", false);

  void processDplusCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesDplusData const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoDplusTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoDplusTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDplusCollisions, "Fill rho tables for collisions with Dplus candidates", false);

  void processDplusMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesDplusMCP const& candidates)
  {
    for (auto& candidate : candidates) {

      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoDplusMcTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoDplusMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDplusMcCollisions, "Fill rho tables for collisions with Dplus MCP candidates", false);

  void processDsCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesDsData const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoDsTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoDsTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDsCollisions, "Fill rho tables for collisions with Ds candidates", false);

  void processDsMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesDsMCP const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoDsMcTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoDsMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDsMcCollisions, "Fill rho tables for collisions with Ds MCP candidates", false);

  void processDstarCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesDstarData const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoDstarTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoDstarTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDstarCollisions, "Fill rho tables for collisions with Dstar candidates", false);

  void processDstarMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesDstarMCP const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoDstarMcTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoDstarMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDstarMcCollisions, "Fill rho tables for collisions with Dstar MCP candidates", false);

  void processLcCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesLcData const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoLcTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoLcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processLcCollisions, "Fill rho tables for collisions with Lc candidates", false);

  void processLcMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesLcMCP const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoLcMcTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoLcMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processLcMcCollisions, "Fill rho tables for collisions with Lc MCP candidates", false);

  void processB0Collisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesB0Data const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoB0Table(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoB0Table(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processB0Collisions, "Fill rho tables for collisions with B0 candidates", false);

  void processB0McCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesB0MCP const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoB0McTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoB0McTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processB0McCollisions, "Fill rho tables for collisions with B0 MCP candidates", false);

  void processBplusCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesBplusData const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoBplusTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoBplusTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processBplusCollisions, "Fill rho tables for collisions with Bplus candidates", false);

  void processBplusMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesBplusMCP const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoBplusMcTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoBplusMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processBplusMcCollisions, "Fill rho tables for collisions with Bplus MCP candidates", false);

  void processXicToXiPiPiCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesXicToXiPiPiData const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoXicToXiPiPiTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoXicToXiPiPiTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processXicToXiPiPiCollisions, "Fill rho tables for collisions with XicToXiPiPi candidates", false);

  void processXicToXiPiPiMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesXicToXiPiPiMCP const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoXicToXiPiPiMcTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoXicToXiPiPiMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processXicToXiPiPiMcCollisions, "Fill rho tables for collisions with XicToXiPiPi MCP candidates", false);

  void processDielectronCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesDielectronData const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections) || collision.centFT0M() < config.centralityMin || collision.centFT0M() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
        rhoDielectronTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, config.applyTrackingEfficiency, config.trackingEfficiency, config.trackingEfficiencyPtBinning, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoDielectronTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDielectronCollisions, "Fill rho tables for collisions with Dielectron candidates", false);

  void processDielectronMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesDielectronMCP const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!jetderiveddatautilities::selectMcCollision(mcCollision, config.skipMBGapEvents, config.applyRCTSelections) || std::abs(mcCollision.posZ()) > config.vertexZCut) {
        rhoDielectronMcTable(0.0, 0.0);
        continue;
      }
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate);

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, config.doSparse);
      rhoDielectronMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDielectronMcCollisions, "Fill rho tables for collisions with Dielectron MCP candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<RhoEstimatorTask>(cfgc, TaskName{"estimator-rho"})}; }
