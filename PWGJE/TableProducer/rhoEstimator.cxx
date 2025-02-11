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

// Task to produce a table joinable to the jcollision table with the mean background pT density
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "Framework/runDataProcessing.h"

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
  Produces<aod::BkgLcRhos> rhoLcTable;
  Produces<aod::BkgLcMcRhos> rhoLcMcTable;
  Produces<aod::BkgBplusRhos> rhoBplusTable;
  Produces<aod::BkgBplusMcRhos> rhoBplusMcTable;
  Produces<aod::BkgDielectronRhos> rhoDielectronTable;
  Produces<aod::BkgDielectronMcRhos> rhoDielectronMcTable;

  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};

  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -99.0, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 99.0, "maximum track phi"};
  Configurable<double> trackingEfficiency{"trackingEfficiency", 1.0, "tracking efficiency applied to jet finding"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  Configurable<int> jetAlgorithm{"jetAlgorithm", 0, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  Configurable<float> bkgjetR{"bkgjetR", 0.2, "jet resolution parameter for determining background density"};
  Configurable<float> bkgEtaMin{"bkgEtaMin", -0.7, "minimim pseudorapidity for determining background density"};
  Configurable<float> bkgEtaMax{"bkgEtaMax", 0.7, "maximum pseudorapidity for determining background density"};
  Configurable<float> bkgPhiMin{"bkgPhiMin", -99.0, "minimim phi for determining background density"};
  Configurable<float> bkgPhiMax{"bkgPhiMax", 99.0, "maximum phi for determining background density"};
  Configurable<bool> doSparse{"doSparse", false, "perfom sparse estimation"};

  JetBkgSubUtils bkgSub;
  float bkgPhiMax_;
  std::vector<fastjet::PseudoJet> inputParticles;
  int trackSelection = -1;
  std::string particleSelection;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    bkgSub.setJetAlgorithmAndScheme(static_cast<fastjet::JetAlgorithm>(static_cast<int>(jetAlgorithm)), static_cast<fastjet::RecombinationScheme>(static_cast<int>(jetRecombScheme)));
    bkgSub.setJetBkgR(bkgjetR);
    bkgSub.setEtaMinMax(bkgEtaMin, bkgEtaMax);
    if (bkgPhiMax > 98.0) {
      bkgPhiMax_ = 2.0 * M_PI;
    }
    bkgSub.setPhiMinMax(bkgPhiMin, bkgPhiMax_);
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta >= trackEtaMin && aod::jmcparticle::eta <= trackEtaMax && aod::jmcparticle::phi >= trackPhiMin && aod::jmcparticle::phi <= trackPhiMax);

  void processChargedCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (skipMBGapEvents && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      rhoChargedTable(-1., -1.);
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseTracks<soa::Filtered<aod::JetTracks>, soa::Filtered<aod::JetTracks>::iterator>(inputParticles, tracks, trackSelection, trackingEfficiency);
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
    rhoChargedTable(rho, rhoM);
  }
  PROCESS_SWITCH(RhoEstimatorTask, processChargedCollisions, "Fill rho tables for collisions using charged tracks", true);

  void processChargedMcCollisions(aod::JetMcCollision const& mcCollision, soa::Filtered<aod::JetParticles> const& particles)
  {
    if (skipMBGapEvents && mcCollision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      rhoChargedTable(-1., -1.);
      return;
    }
    inputParticles.clear();
    jetfindingutilities::analyseParticles<true, soa::Filtered<aod::JetParticles>, soa::Filtered<aod::JetParticles>::iterator>(inputParticles, particleSelection, 1, particles, pdgDatabase);
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
    rhoChargedMcTable(rho, rhoM);
  }
  PROCESS_SWITCH(RhoEstimatorTask, processChargedMcCollisions, "Fill rho tables for MC collisions using charged tracks", false);

  void processD0Collisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesD0Data const& candidates)
  {
    inputParticles.clear();
    for (auto& candidate : candidates) {
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, trackingEfficiency, std::optional{candidate});

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
      rhoD0Table(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processD0Collisions, "Fill rho tables for collisions with D0 candidates", false);

  void processD0McCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesD0MCP const& candidates)
  {
    inputParticles.clear();
    for (auto& candidate : candidates) {
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, std::optional{candidate});

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
      rhoD0McTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processD0McCollisions, "Fill rho tables for collisions with D0 MCP candidates", false);

  void processDplusCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesDplusData const& candidates)
  {
    inputParticles.clear();
    for (auto& candidate : candidates) {
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, trackingEfficiency, std::optional{candidate});

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
      rhoDplusTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDplusCollisions, "Fill rho tables for collisions with Dplus candidates", false);

  void processDplusMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesDplusMCP const& candidates)
  {
    inputParticles.clear();
    for (auto& candidate : candidates) {
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, std::optional{candidate});

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
      rhoDplusMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDplusMcCollisions, "Fill rho tables for collisions with Dplus MCP candidates", false);

  void processLcCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesLcData const& candidates)
  {
    inputParticles.clear();
    for (auto& candidate : candidates) {
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, trackingEfficiency, std::optional{candidate});

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
      rhoLcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processLcCollisions, "Fill rho tables for collisions with Lc candidates", false);

  void processLcMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesLcMCP const& candidates)
  {
    inputParticles.clear();
    for (auto& candidate : candidates) {
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, std::optional{candidate});

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
      rhoLcMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processLcMcCollisions, "Fill rho tables for collisions with Lc MCP candidates", false);

    void processBplusCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesBplusData const& candidates)
    {
      inputParticles.clear();
      for (auto& candidate : candidates) {
        inputParticles.clear();
        jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, trackingEfficiency, std::optional{candidate});

        auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
        rhoBplusTable(rho, rhoM);
      }
    }
    PROCESS_SWITCH(RhoEstimatorTask, processBplusCollisions, "Fill rho tables for collisions with Bplus candidates", false);

    void processBplusMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesBplusMCP const& candidates)
    {
      inputParticles.clear();
      for (auto& candidate : candidates) {
        inputParticles.clear();
        jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, std::optional{candidate});

        auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
        rhoBplusMcTable(rho, rhoM);
      }
    }
    PROCESS_SWITCH(RhoEstimatorTask, processBplusMcCollisions, "Fill rho tables for collisions with Bplus MCP candidates", false);

  void processDielectronCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, aod::CandidatesDielectronData const& candidates)
  {
    inputParticles.clear();
    for (auto& candidate : candidates) {
      inputParticles.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, trackingEfficiency, std::optional{candidate});

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
      rhoDielectronTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDielectronCollisions, "Fill rho tables for collisions with Dielectron candidates", false);

  void processDielectronMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& particles, aod::CandidatesDielectronMCP const& candidates)
  {
    inputParticles.clear();
    for (auto& candidate : candidates) {
      inputParticles.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, std::optional{candidate});

      auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
      rhoDielectronMcTable(rho, rhoM);
    }
  }
  PROCESS_SWITCH(RhoEstimatorTask, processDielectronMcCollisions, "Fill rho tables for collisions with Dielectron MCP candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<RhoEstimatorTask>(cfgc, TaskName{"estimator-rho"})}; }
