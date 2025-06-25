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

/// \file UDCollisionSelExtrasV003Converter.cxx
/// \brief Converts UDCollisionSelExtras table from version 000 to 003 and 001 to 003 and 002 to 003
/// \author Adam Matyja <adam.tomasz.matyja@cern.ch>

#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

// Converts UDCollisions for version 000 to 003 and 001 to 003 and 002 to 003
struct UDCollisionSelExtrasV003Converter {
  Produces<o2::aod::UDCollisionSelExtras_003> udCollisionSelExtras_003;

  void init(InitContext const&)
  {
    if (!doprocessV000ToV003 && !doprocessV001ToV003 && !doprocessV002ToV003) {
      LOGF(fatal, "Neither processV000ToV003 nor processV001ToV003 nor processV002ToV003 is enabled. Please choose one!");
    }
    if (static_cast<int>(doprocessV000ToV003) + static_cast<int>(doprocessV001ToV003) + static_cast<int>(doprocessV002ToV003) > 1) {
      LOGF(fatal, "More than one among processV000ToV003, processV001ToV003, processV002ToV003 is enabled. Please choose only one!");
    }
  }

  void processV000ToV003(o2::aod::UDCollisionSelExtras_000 const& collisions)
  {

    for (const auto& collision : collisions) {

      udCollisionSelExtras_003(collision.chFT0A(),
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
                               0,    // dummy vtxITSTPC
                               0);   // dummy rct
    }
  }
  PROCESS_SWITCH(UDCollisionSelExtrasV003Converter, processV000ToV003, "process v000-to-v003 conversion", false);

  void processV001ToV003(o2::aod::UDCollisionSelExtras_001 const& collisions)
  {

    for (const auto& collision : collisions) {

      udCollisionSelExtras_003(collision.chFT0A(),
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
                               0,  // dummy vtxITSTPC
                               0); // dummy rct
    }
  }
  PROCESS_SWITCH(UDCollisionSelExtrasV003Converter, processV001ToV003, "process v001-to-v003 conversion", false);

  void processV002ToV003(o2::aod::UDCollisionSelExtras_002 const& collisions)
  {

    for (const auto& collision : collisions) {

      udCollisionSelExtras_003(collision.chFT0A(),
                               collision.chFT0C(),
                               collision.chFDDA(),
                               collision.chFDDC(),
                               collision.chFV0A(),
                               collision.occupancyInTime(),
                               collision.hadronicRate(),
                               collision.trs(),
                               collision.trofs(),
                               collision.hmpr(),
                               collision.tfb(),
                               collision.itsROFb(),
                               collision.sbp(),
                               collision.zVtxFT0vPV(),
                               collision.vtxITSTPC(),
                               0); // dummy rct
    }
  }
  PROCESS_SWITCH(UDCollisionSelExtrasV003Converter, processV002ToV003, "process v002-to-v003 conversion", true);
};

/// Spawn the extended table for UDCollisionSelExtras003 to avoid the call to the internal spawner and a consequent circular dependency
// struct UDCollisionSelExtrasSpawner {
//   Spawns<aod::UDCollisionSelExtras_003> udCollisionSelExtras_003;
// };

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDCollisionSelExtrasV003Converter>(cfgc),
    //    adaptAnalysisTask<UDCollisionSelExtrasSpawner>(cfgc),
  };
}
