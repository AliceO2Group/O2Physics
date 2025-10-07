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
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;

// Converts bc_000 into bc_001
struct bcConverter {
  Produces<aod::BCs_001> bc_001;

  void process(aod::BCs_000 const& bcTable)
  {
    for (auto& bc : bcTable) {
      constexpr uint64_t lEmptyTriggerInputs = 0;
      bc_001(bc.runNumber(), bc.globalBC(), bc.triggerMask(), lEmptyTriggerInputs);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<bcConverter>(cfgc),
  };
}
