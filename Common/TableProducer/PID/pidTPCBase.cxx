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
  SliceCache cache;
  Produces<aod::PIDMults> mult;

  bool enableTable = false;
  void init(InitContext& initContext)
  {
    LOG(info) << "Initializing PID Mult Task";
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "PIDMults");
    if (enableTable) {
      LOG(info) << "Table TPC PID Multiplicity enabled!";
    }
    if (doprocessStandard == true && doprocessIU == true) {
      LOG(fatal) << "Both processStandard and processIU are enabled, pick one!";
    }
  }

  using TrksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
  Partition<TrksIU> tracksWithTPCIU = (aod::track::tpcNClsFindable > (uint8_t)0);
  void processIU(aod::Collision const& collision, TrksIU const&)
  {
    if (!enableTable) {
      return;
    }
    auto tracksGrouped = tracksWithTPCIU->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    mult(tracksGrouped.size());
  }
  PROCESS_SWITCH(PidMultiplicity, processIU, "Process with IU tracks, faster but works on Run3 only", false);

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra>;
  Partition<Trks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  void processStandard(aod::Collision const& collision, Trks const&)
  {
    if (!enableTable) {
      return;
    }
    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    mult(tracksGrouped.size());
  }
  PROCESS_SWITCH(PidMultiplicity, processStandard, "Process with tracks, needs propagated tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidMultiplicity>(cfgc)};
}
