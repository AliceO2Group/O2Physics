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

/// \file UDCollisionSelExtrasV002Converter.cxx
/// \brief Converts UDCollisionSelExtras table from version 000 to 002 and 001 to 002
/// \author Roman Lavicka <roman.lavicka@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2;
using namespace o2::framework;

// Converts UDCollisions for version 000 to 002 and 001 to 002
struct UDCollisionSelExtrasV002Converter {
  Produces<o2::aod::UDCollisionSelExtras_002> udCollisionSelExtras_002;

  void init(InitContext const&)
  {
    if (!doprocessV000ToV002 && !doprocessV001ToV002) {
      LOGF(fatal, "Neither processV000ToV002 nor processV001ToV002 is enabled. Please choose one!");
    }
    if (doprocessV000ToV002 && doprocessV001ToV002) {
      LOGF(fatal, "Both processV000ToV002 and processV001ToV002 are enabled. Please choose only one!");
    }
  }

  void processV000ToV002(o2::aod::UDCollisionSelExtras_000 const& collisions)
  {

    for (const auto& collision : collisions) {

      udCollisionSelExtras_002(collision.chFT0A(),
                               collision.chFT0C(),
                               collision.chFDDA(),
                               collision.chFDDC(),
                               collision.chFV0A(),
                               0,    // dummy occupancy
                               0.0f, // dummy rate
                               0,    // dummy trs
                               0,    // dummy trofs
                               0,    // dummy hmpr
                               0,    // dummy tfb
                               0,    // dummy itsROFb
                               0,    // dummy sbp
                               0,    // dummy zVtxFT0vPV
                               0);   // dummy vtxITSTPC
    }
  }
  PROCESS_SWITCH(UDCollisionSelExtrasV002Converter, processV000ToV002, "process v000-to-v002 conversion", false);

  void processV001ToV002(o2::aod::UDCollisionSelExtras_001 const& collisions)
  {

    for (const auto& collision : collisions) {

      udCollisionSelExtras_002(collision.chFT0A(),
                               collision.chFT0C(),
                               collision.chFDDA(),
                               collision.chFDDC(),
                               collision.chFV0A(),
                               collision.occupancyInTime(),
                               collision.hadronicRate(),
                               collision.trs(),
                               collision.trofs(),
                               collision.hmpr(),
                               0,  // dummy tfb
                               0,  // dummy itsROFb
                               0,  // dummy sbp
                               0,  // dummy zVtxFT0vPV
                               0); // dummy vtxITSTPC
    }
  }
  PROCESS_SWITCH(UDCollisionSelExtrasV002Converter, processV001ToV002, "process v001-to-v002 conversion", true);
};

/// Spawn the extended table for UDCollisionSelExtras002 to avoid the call to the internal spawner and a consequent circular dependency
// struct UDCollisionSelExtrasSpawner {
//   Spawns<aod::UDCollisionSelExtras_002> udCollisionSelExtras_002;
// };

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDCollisionSelExtrasV002Converter>(cfgc),
    //    adaptAnalysisTask<UDCollisionSelExtrasSpawner>(cfgc),
  };
}
