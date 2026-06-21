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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
// Task used to convert the data model from the old format to the new format. To avoid
// the conflict with the old data model.

#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct reducedMCeventConverter002 {
  Produces<aod::ReducedMCEvents_002> reducedMCevent_002;

  void init(InitContext const&)
  {
    if (doprocessV000ToV002 == false && doprocessV001ToV002 == false) {
      LOGF(fatal, "Neither processV000ToV002 nor processV001ToV002 is enabled. Please choose one!");
    }
    if (doprocessV000ToV002 == true && doprocessV001ToV002 == true) {
      LOGF(fatal, "Both processV000ToV002 and processV001ToV002 are enabled. Please choose only one!");
    }
  }

  void processV000ToV002(aod::ReducedMCEvents_000 const& events)
  {
    for (const auto& event : events) {
      uint64_t globalBc = 0;
      reducedMCevent_002(globalBc, event.generatorsID(), event.mcPosX(), event.mcPosY(), event.mcPosZ(),
                         event.t(), event.weight(), event.impactParameter(),
                         -1.0f, -1.0f, -1.0f, -1.0f);
    }
  }
  PROCESS_SWITCH(reducedMCeventConverter002, processV000ToV002, "process v000-to-v002 conversion", false);

  void processV001ToV002(aod::ReducedMCEvents_001 const& events)
  {
    for (const auto& event : events) {
      uint64_t globalBc = 0;
      reducedMCevent_002(globalBc, event.generatorsID(), event.mcPosX(), event.mcPosY(), event.mcPosZ(),
                         event.t(), event.weight(), event.impactParameter(),
                         event.centFT0C(), event.multMCNParticlesEta05(), event.multMCNParticlesEta08(), event.multMCNParticlesEta10());
    }
  }
  PROCESS_SWITCH(reducedMCeventConverter002, processV001ToV002, "process v001-to-v002 conversion", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<reducedMCeventConverter002>(cfgc)};
}
