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

/// \file UDCollisionFITBitsConverter.cxx
/// \brief Converts UDCollisionFITBits table from version 000 to 001
///
/// This task converts the UDCollisionFITBits table from version 000,
/// without an explicit UDCollisionId link, to version 001, which includes
/// the UDCollisionId column.
///
/// The conversion assumes the version-000 table has one row per UDCollision
/// and follows the same ordering as the UDCollisions table.
///
/// executable name o2-analysis-ud-fit-bits-converter
///
/// \author Sandor Lokos <sandor.lokos@cern.ch>

#include "PWGUD/DataModel/UDTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

// Converts UDCollisionFITBits from version 000 to 001
struct UDCollisionFITBitsConverter {
  Produces<o2::aod::UDCollisionFITBits_001> udCollisionFITBits_001;

  void process(o2::aod::UDCollisionFITBits_000 const& fitBits)
  {
    udCollisionFITBits_001.reserve(fitBits.size());

    for (const auto& fitBit : fitBits) {
      udCollisionFITBits_001(fitBit.globalIndex(), // UDCollisionId; version 000 has no explicit link
                             fitBit.thr1W0(),
                             fitBit.thr1W1(),
                             fitBit.thr1W2(),
                             fitBit.thr1W3(),
                             fitBit.thr2W0(),
                             fitBit.thr2W1(),
                             fitBit.thr2W2(),
                             fitBit.thr2W3());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDCollisionFITBitsConverter>(cfgc),
  };
}