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
///
/// \file collisionCuts.cxx
/// \brief Tutorial: Using CollisionCutsGroup for automatic configurable registration
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
/// \since 30/12/2024

#include "PWGLF/Utils/collisionCuts.h"

#include "PWGLF/Utils/collisionCutsGroup.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Tutorial: Collision Cuts with Automatic Configurable Registration
// This tutorial demonstrates the new CollisionCutsGroup pattern that enables
// automatic configurable registration, reducing boilerplate from ~30 lines to ~5 lines.
struct CollisionCutsTutorial { // o2-linter: disable=name/workflow-file (Tutorial file naming convention)
  // NEW PATTERN: Single line configurable declaration!
  // All collision selection parameters are automatically registered by DPL framework
  Configurable<CollisionCutsGroup> collCuts{"collCuts", {}, "Collision event selection cuts"};

  // Collision cuts checker object
  o2::analysis::CollisonCuts colCutsChecker;

  // Histogram registry
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    LOGF(info, "==========================================================");
    LOGF(info, "Tutorial: collisionCutsGroup Automatic Configurable Registration");
    LOGF(info, "==========================================================");
    LOGF(info, "");

    // Print all settings from CollisionCutsGroup
    LOGF(info, "Step 1: Printing CollisionCutsGroup settings...");
    collCuts->printSelections();

    // Apply settings from CollisionCutsGroup to CollisonCuts FIRST
    LOGF(info, "Step 2: Applying CollisionCutsGroup settings to CollisonCuts...");
    colCutsChecker.setCuts(collCuts, /*isRun3*/ true);

    // Then initialize collision cuts checker (creates histograms based on settings)
    LOGF(info, "Step 3: Initializing CollisonCuts object...");
    colCutsChecker.init(&histos);

    // Print CollisonCuts settings for verification
    LOGF(info, "Step 4: Printing CollisonCuts settings after application...");
    colCutsChecker.printCuts();

    LOGF(info, "");
    LOGF(info, "==========================================================");
    LOGF(info, "Tutorial: Initialization Complete!");
    LOGF(info, "==========================================================");
    LOGF(info, "");
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As>::iterator const& collision)
  {
    // Test collision selection
    bool isSelected = colCutsChecker.isSelected(collision);

    if (isSelected) {
      LOGF(debug, "Collision passed selection: z-vertex = %.2f cm", collision.posZ());
      colCutsChecker.fillQA(collision);
    } else {
      LOGF(debug, "Collision rejected by selection");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CollisionCutsTutorial>(cfgc)};
}
