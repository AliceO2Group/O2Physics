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
/// \file   pidTPCBase.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TPC PID tasks.
///

#include <utility>
#include <vector>
#include <string>

// O2 includes
#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/FT0Corrected.h"
#include "TableHelper.h"
#include "pidTPCBase.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

struct PidMultiplicity {
  Produces<aod::PIDMults> mult;

  // For vertex-Z corrections in calibration
  Partition<soa::Join<aod::Tracks, aod::TracksExtra>> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  bool enableTable = false;

  void init(InitContext& initContext)
  {
    LOG(info) << "Initializing PID Mult Task";
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "PIDMults");
    if (enableTable) {
      LOG(info) << "Table TPC PID Multiplicity enabled!";
    }
  }

  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra> const& tracksExtra)
  {
    if (!enableTable) {
      return;
    }
    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex());
    mult(tracksGrouped.size());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidMultiplicity>(cfgc)};
}
