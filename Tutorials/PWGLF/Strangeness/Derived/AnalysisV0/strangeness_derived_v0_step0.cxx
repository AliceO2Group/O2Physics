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
/// \brief Step0 of the V0 strangeness tutorial on derived data (K0S, Lambda, AntiLambda)
/// \author Sara Pucillo
/// based on the O2AT strangeness derived-data tutorial from:
/// \author Romain Schotter
/// and on the original codes from:
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all V0s and fill invariant mass histograms

struct StrangenessDerivedV0Step0 {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rK0Short{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda{"lambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0ShortMassAxis = {100, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec LambdaMassAxis = {100, 1.08f, 1.15f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // K0S/Lambda/AntiLambda reconstruction
    rK0Short.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {K0ShortMassAxis}});

    rLambda.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {LambdaMassAxis}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  void process(soa::Filtered<soa::Join<aod::StraCollisions, aod::StraEvSels>>::iterator const& collision,
               aod::V0Cores const& V0s)
  {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    // V0s
    for (const auto& v0 : V0s) {
      rK0Short.fill(HIST("hMassK0Short"), v0.mK0Short());
      rLambda.fill(HIST("hMassLambda"), v0.mLambda());
      rLambda.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangenessDerivedV0Step0>(cfgc, TaskName{"strangeness-derived-v0-tutorial"})};
}
