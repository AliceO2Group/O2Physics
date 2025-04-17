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
/// \file   udMcCollisions2udCollisions.cxx
/// \author Roman Lavička
/// \since  2025-04-15
/// \brief  A task to create a reverse index from UDMcCollisions to UDCollisions
///

#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/DataModel/UDIndex.h"

using namespace o2;
using namespace o2::framework;

struct UDMcCollisions2UDCollisions {
  using LabeledCollisions = soa::Join<aod::UDCollisions, aod::UDMcCollsLabels>;
  Produces<aod::UDMcCollisionsToUDCollisions> udmcc2udc;

  std::vector<int> collisionIds;

  void init(InitContext&)
  {
  }

  void process(aod::UDMcCollisions const& mcCollisions)
  {
    if (doprocessIndexingCentral || doprocessIndexingCentralFast) {
      udmcc2udc.reserve(mcCollisions.size());
    }
  }

  void processIndexingCentralFast(aod::UDMcCollisions const& mcCollisions, LabeledCollisions const& collisions)
  {
    // faster version, but will use more memory due to pre-allocation
    std::vector<std::vector<int>> mccoll2coll(mcCollisions.size());
    for (const auto& collision : collisions) {
      if (collision.has_udMcCollision())
        mccoll2coll[collision.udMcCollisionId()].push_back(collision.globalIndex());
    }
    for (const auto& mcCollision : mcCollisions) {
      udmcc2udc(mccoll2coll[mcCollision.globalIndex()]);
    }
  }
  PROCESS_SWITCH(UDMcCollisions2UDCollisions, processIndexingCentralFast, "Create reverse index from mccollisions to collision: more memory use but potentially faster", true);

  void processIndexingCentral(aod::UDMcCollisions const&, soa::SmallGroups<LabeledCollisions> const& collisions)
  {
    collisionIds.clear();
    for (const auto& collision : collisions) {
      collisionIds.push_back(collision.globalIndex());
    }
    udmcc2udc(collisionIds);
  }
  PROCESS_SWITCH(UDMcCollisions2UDCollisions, processIndexingCentral, "Create reverse index from mccollisions to collision", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UDMcCollisions2UDCollisions>(cfgc)};
}
