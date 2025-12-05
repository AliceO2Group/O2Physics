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

/// \file slimTablesProducer.cxx
/// \brief Task to produce a reduced version of Tables for tracks, collisions, mcparticles and mccollisions.
/// \author Millot Louise <louise.millot@cern.ch>

#include "PWGJE/DataModel/SlimTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SlimTablesProducer {

  void init(InitContext&)
  {
  }

  Produces<o2::aod::SlimCollisions> slimCollisions;
  Produces<o2::aod::SlimMcCollisions> slimMcCollisions;
  Produces<o2::aod::SlimTracks> slimTracks;
  Produces<o2::aod::SlimParticles> slimParticles;

  void processCollision(aod::SlimCollisions const& collisions)
  {
    for (const auto& coll : collisions) {
      slimCollisions(coll.posZ(), coll.centFT0C(), coll.centFT0M(), coll.weight(), coll.eventSel(), coll.trackOccupancyInTimeRange());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processCollision, "Produce slim collision table", true);

  void processMcCollision(aod::SlimMcCollisions const& mccollisions)
  {
    for (const auto& mccoll : mccollisions) {
      slimMcCollisions(mccoll.posZ(), mccoll.centFT0M(), mccoll.weight(), mccoll.accepted(), mccoll.ptHard());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processMcCollision, "Produce slim mc collision table", true);

  void processTracks(aod::SlimTracks const& tracks)
  {
    for (const auto& trk : tracks) {
      slimTracks(trk.collision(), trk.pt(), trk.eta(), trk.phi(), trk.dcaXY());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processTracks, "Produce slim track table", true);

  void processParticles(aod::McParticles const& parts)
  {
    for (const auto& p : parts) {
      slimParticles(p.mcCollision(), p.pt(), p.eta(), p.phi());
    }
  }
  PROCESS_SWITCH(SlimTablesProducer, processParticles, "Produce slim particles", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SlimTablesProducer>(cfgc)};
}
