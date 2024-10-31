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
  Produces<aod::BkgD0Rhos> rhoD0Table;
  Produces<aod::BkgLcRhos> rhoLcTable;
  Produces<aod::BkgBplusRhos> rhoBplusTable;
  Produces<aod::BkgDielectronRhos> rhoDielectronTable;

  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<double> trackingEfficiency{"trackingEfficiency", 1.0, "tracking efficiency applied to jet finding"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<float> bkgjetR{"bkgjetR", 0.2, "jet resolution parameter for determining background density"};
  Configurable<float> bkgEtaMin{"bkgEtaMin", -0.9, "minimim pseudorapidity for determining background density"};
  Configurable<float> bkgEtaMax{"bkgEtaMax", 0.9, "maximum pseudorapidity for determining background density"};
  Configurable<float> bkgPhiMin{"bkgPhiMin", 0., "minimim phi for determining background density"};
  Configurable<float> bkgPhiMax{"bkgPhiMax", 99.0, "maximum phi for determining background density"};
  Configurable<bool> doSparse{"doSparse", false, "perfom sparse estimation"};

  JetBkgSubUtils bkgSub;
  float bkgPhiMax_;
  std::vector<fastjet::PseudoJet> inputParticles;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    bkgSub.setJetBkgR(bkgjetR);
    bkgSub.setEtaMinMax(bkgEtaMin, bkgEtaMax);
    if (bkgPhiMax > 98.0) {
      bkgPhiMax_ = 2.0 * M_PI;
    }
    bkgSub.setPhiMinMax(bkgPhiMin, bkgPhiMax_);
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);

  void processChargedCollisions(aod::JetCollision const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    inputParticles.clear();
    jetfindingutilities::analyseTracks<soa::Filtered<aod::JetTracks>, soa::Filtered<aod::JetTracks>::iterator>(inputParticles, tracks, trackSelection, trackingEfficiency);
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
    rhoChargedTable(rho, rhoM);
  }
  PROCESS_SWITCH(RhoEstimatorTask, processChargedCollisions, "Fill rho tables for collisions using charged tracks", true);

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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<RhoEstimatorTask>(cfgc, TaskName{"estimator-rho"})}; }
