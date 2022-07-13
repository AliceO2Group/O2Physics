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
/// \author Christian Sonnabend christian.sonnabend@cern.ch
/// \author Annalena Kalteyer annalena.sophie.kalteyer@cern.ch
/// \brief  Base to build tasks for TPC PID tasks.
///

#include "Framework/Configurable.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include "TableHelper.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "pidTPCBase.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the TPC multiplicity table
struct tpcMultiplicity {
  // Tables to produce
  Produces<o2::aod::TPCMult> tableMult;

  bool enableTable = false;

  void init(o2::framework::InitContext& initContext)
  {
    if (doprocessRun2 == true && doprocessRun3 == true) {
      LOGF(fatal, "Cannot enable processRun2 and processRun3 at the same time. Please choose one.");
    }

    enableTable = isTableRequiredInWorkflow(initContext, "TPCMult");
    if (enableTable) {
      LOG(info) << "Table TPCMult enabled!";
    }
  }
  Partition<soa::Join<aod::TracksIU, aod::TracksExtra>> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);

  ///
  /// Process function to prepare the event for each track on Run 2 data
  void processRun2(aod::Run2MatchedSparse::iterator const& collision, soa::Join<aod::TracksIU, aod::TracksExtra> const& tracksExtra)
  {
    if (!enableTable) {
      return;
    }

    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex());
    tableMult(tracksGrouped.size());
  }
  PROCESS_SWITCH(tpcMultiplicity, processRun2, "Process with Run2 data", false);

  ///
  /// Process function to prepare the event for each track on Run 3 data
  void processRun3(aod::Collisions::iterator const& collision, soa::Join<aod::TracksIU, aod::TracksExtra> const& tracksExtra)
  {
    if (!enableTable) {
      return;
    }
    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex());
    tableMult(tracksGrouped.size());
  }

  PROCESS_SWITCH(tpcMultiplicity, processRun3, "Process with Run3 data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tpcMultiplicity>(cfgc)};
}
