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

/// \file UDCollisionFITBitsConverter.cxx
/// \brief Converts UDCollisionFITBits from version 000 to 001
///
/// The legacy table was written with exactly one FIT-bit row for every
/// UDCollision, in identical row order. Therefore row i of
/// UDCollisionFITBits_000 belongs to row i of UDCollisions. Because the
/// legacy table does not store this association explicitly, the converter
/// materialises it by writing the FIT-bit row's global index as
/// UDCollisionId. The eight packed FIT-bit words are copied unchanged.
///
/// Executable: o2-analysis-ud-collision-fit-bits-converter
///
/// \author Sandor Lokos, sandor.lokos@cern.ch

#include "PWGUD/DataModel/UDTables.h"

#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct UDCollisionFITBitsConverter {
  Produces<aod::UDCollisionFITBits_001> outputFITBits;

  void process(aod::UDCollisionFITBits_000 const& legacyFITBits)
  {
    outputFITBits.reserve(legacyFITBits.size());

    for (const auto& bits : legacyFITBits) {
      outputFITBits(bits.globalIndex(),
                    bits.thr1W0(),
                    bits.thr1W1(),
                    bits.thr1W2(),
                    bits.thr1W3(),
                    bits.thr2W0(),
                    bits.thr2W1(),
                    bits.thr2W2(),
                    bits.thr2W3());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDCollisionFITBitsConverter>(cfgc)};
}
