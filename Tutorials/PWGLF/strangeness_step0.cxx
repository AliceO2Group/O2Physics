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
/// \brief this is a starting point for the Strangeness tutorial
/// \author
/// \since 12/05/2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all V0s and fill invariant mass histogram
struct strangeness_tutorial {

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{"registry",
                             {{"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
                              {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}}}};

  // Defining filters for events (event selection)
  // Processed events will be already fulfulling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               aod::V0Datas const& V0s)
  {
    // Fill the event counter
    registry.fill(HIST("hVertexZ"), collision.posZ());

    for (auto& v0 : V0s) {
      registry.fill(HIST("hMassK0Short"), v0.mK0Short());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<strangeness_tutorial>(cfgc)}; }
