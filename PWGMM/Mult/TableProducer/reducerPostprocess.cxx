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

#include "Framework/AnalysisTask.h"

#include "ReducedTables.h"

#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ReducerPostprocess {
  Produces<aod::RFeatMins> features;

  using MCFull = soa::Join<aod::RMCCollisions, aod::RHepMCinfos>;
  void processFull(MCFull const&,
                   soa::Join<aod::RCollisions, aod::RMCColLabels> const& cols)
  {
    for (auto& col : cols) {
      auto mccol = col.rmccollision_as<MCFull>();
      features(mccol.multMCNParticlesEta10(), col.multNTracksPVeta1(), mccol.processId(), mccol.impactParameter(), col.posX(), col.posY(), col.posZ(), col.collisionTimeRes(), col.multFT0A(), col.multFT0C());
    }
  }

  PROCESS_SWITCH(ReducerPostprocess, processFull, "Process with HepMC", false);

  using MCLite = aod::RMCCollisions;
  void processLite(MCLite const&,
                   soa::Join<aod::RCollisions, aod::RMCColLabels> const& cols)
  {
    for (auto& col : cols) {
      auto mccol = col.rmccollision_as<MCLite>();
      features(mccol.multMCNParticlesEta10(), col.multNTracksPVeta1(), -1, -1.f, col.posX(), col.posY(), col.posZ(), col.collisionTimeRes(), col.multFT0A(), col.multFT0C());
    }
  }

  PROCESS_SWITCH(ReducerPostprocess, processLite, "Process without HepMC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<ReducerPostprocess>(cfgc)};
}
