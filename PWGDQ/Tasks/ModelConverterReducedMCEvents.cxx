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

// other includes
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <THashList.h>
#include <TList.h>
#include <TString.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct reducedMCeventConverter000_001 {
  Produces<aod::ReducedMCEvents_001> reducedMCevent_001;

  void process(aod::ReducedMCEvents_000 const& events)
  {
    for (const auto& event : events) {
      reducedMCevent_001(event.generatorsID(), event.mcPosX(), event.mcPosY(), event.mcPosZ(),
                         event.t(), event.weight(), event.impactParameter(),
                         -1.0f, -1.0f, -1.0f, -1.0f);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<reducedMCeventConverter000_001>(cfgc)};
}
