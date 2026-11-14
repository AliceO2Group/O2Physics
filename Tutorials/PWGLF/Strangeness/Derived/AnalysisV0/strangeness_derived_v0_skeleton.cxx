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
/// \brief Skeleton of the V0 strangeness tutorial on derived data (K0S, Lambda, AntiLambda)
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

// SKELETON
// Starting point of the V0 tutorial: loop over the selected events (sel8, |z| < 10 cm)
// and fill a histogram with the z position of the primary vertex.
// Goal of the hands-on: improve this file step by step up to step 4
// (checkpoints: strangeness_derived_v0_stepX.cxx).

struct StrangenessDerivedV0Skeleton {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  void process(soa::Filtered<soa::Join<aod::StraCollisions, aod::StraEvSels>>::iterator const& collision)
  {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangenessDerivedV0Skeleton>(cfgc, TaskName{"strangeness-derived-v0-tutorial"})};
}
