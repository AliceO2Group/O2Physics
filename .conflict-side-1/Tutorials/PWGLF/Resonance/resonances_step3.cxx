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
/// \brief this is a starting point for the Resonances tutorial for MC part
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>
/// \since 11/10/2024

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 3
// Starting point for MC: loop over all Generated MC particles
struct resonances_tutorial {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  // Configurable for min pT cut
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minimum pt cut"};

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
  }

  // MC particle selection
  template <typename ParticleType>
  bool ptCut(const ParticleType resoParents)
  {
    // basic pt cuts
    if (std::abs(resoParents.pt()) < cMinPtcut)
      return false;

    return true;
  }

  // Fill histograms (main function)
  template <typename CollisionType, typename ParticleType>
  void fillHistograms(const CollisionType& /*collision*/, const ParticleType& resoParents)
  {
    for (auto part : resoParents) { // loop over all resoParents
      if (!ptCut(part))
        continue; // pt selection

      // QA plots
      histos.fill(HIST("hEta"), part.eta());
    }
  }

  // Process the MC
  void process(aod::ResoCollision& collision, aod::ResoMCParents& resoParents)
  {
    // Fill the event counter
    histos.fill(HIST("hVertexZ"), collision.posZ());
    fillHistograms(collision, resoParents); // Fill histograms, MC
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<resonances_tutorial>(cfgc)}; }
