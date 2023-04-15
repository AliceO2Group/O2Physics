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
/// \brief Femtodream Tutorial 0
/// \author Luca Barioglio, Anton Riedel

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Example task illustrating how to create histograms and fill them with basic
// information. A basic event selection is applied.

struct CFTutorialTask0 {
  HistogramRegistry histos{
    "Histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // TODO:
    // Define your axes
    // Constant bin width axis:

    // Variable bin width axis:

    // TODO:
    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("dummy", "dummy-x ; dummy-y", kTH1F, {{100, -10, 10}});
  }

  // Equivalent of the AliRoot task UserExec
  void process(aod::Collision const& coll, aod::Tracks const& inputTracks)
  {

    // TODO:
    // Fill histograms before and after event selection
    // Performing the event selection, i.e. z vertex cut

    // TODO:
    // Loop over tracks and fill histograms
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<CFTutorialTask0>(cfgc)};
  return workflow;
}
