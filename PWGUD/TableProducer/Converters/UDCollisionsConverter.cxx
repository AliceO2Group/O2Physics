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

/// \file UDCollisionsConverter.cxx
/// \brief Converts UDCollisions table from version 000 to 001

/// This task allows for the conversion of the UDCollisions table from the version 000,
/// that is before the introduction of Flags for UPC reconstruction
/// to the version 001, that includes it

/// executable name o2-analysis-ud-collisions-converter

/// \author Sasha Bylinkin <alexandr.bylinkin@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2;
using namespace o2::framework;

// Converts UDCollisions for version 000 to 001
struct UDCollisionsConverter {
  Produces<o2::aod::UDCollisions_001> udCollisions_001;

  void process(o2::aod::UDCollisions_000 const& collisions)
  {

    for (const auto& collision : collisions) {

      udCollisions_001(collision.globalBC(),
                       collision.runNumber(),
                       collision.posX(),
                       collision.posY(),
                       collision.posZ(),
                       0.0f, // dummy UPC reco flag, not available in version 000
                       collision.numContrib(),
                       collision.netCharge(),
                       collision.rgtrwTOF());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDCollisionsConverter>(cfgc),
  };
}
