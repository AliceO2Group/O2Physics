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

/// \file UDCollisionSelExtrasConverter.cxx
/// \brief Converts UDCollisionSelExtras table from version 000 to 001

/// This task allows for the conversion of the UDCollisionSelExtras table from the version 000,
/// to include occupancy and interaction rate
/// to the version 001, that includes it

/// executable name o2-analysis-ud-collisionselectras-converter

/// \author Sasha Bylinkin <alexandr.bylinkin@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2;
using namespace o2::framework;

// Converts UDCollisions for version 000 to 001
struct UDCollisionSelExtrasConverter {
  Produces<o2::aod::UDCollisionSelExtras_001> udCollisionSelExtras_001;

  void process(o2::aod::UDCollisionSelExtras_000 const& collisions)
  {

    for (const auto& collision : collisions) {

      udCollisionSelExtras_001(collision.chFT0A(),
                               collision.chFT0C(),
                               collision.chFDDA(),
                               collision.chFDDC(),
                               collision.chFV0A(),
                               0,    // dummy occupancy
                               0.0f, // dummy rate
                               0,    // dummy trs
                               0,    // dummy trofs
                               0);   // dummy hmpr
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDCollisionSelExtrasConverter>(cfgc),
  };
}
