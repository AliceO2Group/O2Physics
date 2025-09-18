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

using namespace o2;
using namespace o2::framework;

struct mcCollisionConverter {
  Produces<aod::McCollisions_001> mcCollisions_001;

  void process(aod::McCollisions_000 const& mcCollisionTable)
  {
    for (auto& mcCollision : mcCollisionTable) {

      // Repopulate new table
      mcCollisions_001(
        mcCollision.bcId(),
        mcCollision.generatorsID(),
        mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
        mcCollision.t(), mcCollision.weight(),
        mcCollision.impactParameter(),
        0.0f); // dummy event plane, not available in _000
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mcCollisionConverter>(cfgc),
  };
}
