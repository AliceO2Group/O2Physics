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

// Task to produce a table joinable to the collision table with the mean background pT density
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct RhoEstimatorTask {
  Produces<aod::JCollisionRhos> rhoTable;

  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> bkgjetR{"bkgjetR", 0.2, "jet resolution parameter for determining background density"};
  Configurable<float> bkgEtaMin{"bkgEtaMin", -0.9, "minimim pseudorapidity for determining background density"};
  Configurable<float> bkgEtaMax{"bkgEtaMax", 0.9, "maximum pseudorapidity for determining background density"};
  Configurable<float> bkgPhiMin{"bkgPhiMin", -99., "minimim phi for determining background density"};
  Configurable<float> bkgPhiMax{"bkgPhiMax", 99., "maximum phi for determining background density"};
  Configurable<bool> doSparse{"doSparse", false, "perfom sparse estimation"};

  JetBkgSubUtils bkgSub;
  std::vector<fastjet::PseudoJet> inputParticles;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    trackSelection = JetDerivedDataUtilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    bkgSub.setJetBkgR(bkgjetR);
    bkgSub.setEtaMinMax(bkgEtaMin, bkgEtaMax);
    bkgSub.setPhiMinMax(bkgPhiMin, bkgPhiMax);
  }

  void processCollisions(aod::JCollisions const& collision, aod::JTracks const& tracks)
  {
    inputParticles.clear();
    for (auto& track : tracks) {
      if (!JetDerivedDataUtilities::selectTrack(track, trackSelection)) {
        continue;
      }
      FastJetUtilities::fillTracks(track, inputParticles, track.globalIndex());
    }

    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(inputParticles, doSparse);
    rhoTable(rho, rhoM);
  }
  PROCESS_SWITCH(RhoEstimatorTask, processCollisions, "Fill rho tables for collisions", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<RhoEstimatorTask>(cfgc, TaskName{"rho-estimator"})}; }
