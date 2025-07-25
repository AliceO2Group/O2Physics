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
/// \brief Dynamic columns are computed on-the-fly when attached to an existing table
/// \author
/// \since

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

namespace o2::aod
{
namespace extension
{
DECLARE_SOA_EXPRESSION_COLUMN(P2, p2, float, track::p* track::p);
} // namespace extension
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;

struct ExtendTable {
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void process(aod::Collisions const& collisions, aod::Tracks const& tracks)
  {
    // note that this needs to be done only once, as it is done for the whole table
    // the extension is temporary and is lost when the variable it is assigned to
    // goes out of scope
    auto table_extension = soa::Extend<aod::Tracks, aod::extension::P2>(tracks);
    for (auto& collision : collisions) {
      auto trackSlice = table_extension.sliceBy(perCollision, collision.globalIndex());
      LOGP(info, "Collision {}", collision.globalIndex());
      for (auto& row : trackSlice) {
        if (row.trackType() != 3) {
          if (row.index() % 100 == 0) {
            LOGP(info, "P^2 = {}", row.p2());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExtendTable>(cfgc),
  };
}
