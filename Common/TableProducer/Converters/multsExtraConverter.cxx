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

#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct MultsExtraConverter {
  Produces<aod::MultsExtra_001> multsExtra_001;
  void process(aod::MultsExtra_000 const& multsExtra_000)
  {
    for (const auto& r : multsExtra_000) {
      multsExtra_001(r.multPVTotalContributors(), r.multPVChi2(),
                     r.multCollisionTimeRes(), r.multRunNumber(), r.multPVz(), r.multSel8(),
                     r.multNTracksHasITS(), r.multNTracksHasTPC(), r.multNTracksHasTOF(),
                     r.multNTracksHasTRD(), r.multNTracksITSOnly(),
                     r.multNTracksTPCOnly(), r.multNTracksITSTPC(),
                     r.multAllTracksTPCOnly(), r.multAllTracksITSTPC(),
                     r.trackOccupancyInTimeRange(),
                     0.0f,
                     r.flags());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultsExtraConverter>(cfgc)};
}
