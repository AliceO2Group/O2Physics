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

/// \file upcTrueGapService.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Apr/7/2026
/// \brief task providing truegap

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UpcService.h"
#include "PWGUD/DataModel/UDTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct UpcTrueGapService {
  SGSelector sgSelector;
  Configurable<float> cfgCutFV0{"cfgCutFV0", 50., "FV0A threshold"};
  Configurable<float> cfgCutFT0A{"cfgCutFT0A", 150., "FT0A threshold"};
  Configurable<float> cfgCutFT0C{"cfgCutFT0C", 50., "FT0C threshold"};
  Configurable<float> cfgCutZDC{"cfgCutZDC", 10., "ZDC threshold"};

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDCollisionSelExtras>;

  Produces<aod::Truegapside> truegapside;
  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::vector<double> binEdges = {-1.5, -0.5, 0.5, 1.5, 2.5, 3.5};
    AxisSpec axisgap = {binEdges, "true gap side"};
    registry.add("truegap", "truegap", {HistType::kTH1D, {axisgap}});
  }

  void process(UDCollisionsFull::iterator const& collision)
  {
    auto truegap = sgSelector.trueGap(collision, cfgCutFV0, cfgCutFT0A, cfgCutFT0C, cfgCutZDC);
    truegapside(truegap);
    registry.fill(HIST("truegap"), truegap);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTrueGapService>(cfgc),
  };
}
