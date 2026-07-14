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

/// \file multBCsConverter.cxx
/// \brief Converts MultBCs and MultBcSel table from version 000 to 001
/// \author Jesper Karlsson Gumrpecht

#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct MultBCsConverter {
  Produces<aod::MultBCs_001> multBC;
  Produces<aod::MultBcSel_001> multBcSel;

  static constexpr float DummyValue = -1.f;
  static constexpr int DummyRct = 0;
  void process(soa::Join<aod::MultBCs_000, aod::MultBcSel_000> const& multBCs)
  {
    multBC.reserve(multBCs.size());
    multBcSel.reserve(multBCs.size());
    for (const auto& multbc : multBCs) {
      multBC(
        multbc.multFT0A(),
        multbc.multFT0C(),
        multbc.multFV0A(),
        multbc.multFDDA(),
        multbc.multFDDC(),
        multbc.multZNA(),
        multbc.multZNC(),
        multbc.multZEM1(),
        multbc.multZEM2(),
        multbc.multZPA(),
        multbc.multZPC(),
        DummyValue,  // dummy amplitude for FV0A Outer
        DummyValue); // dummy amplitude for FT0A Outer

      multBcSel(
        multbc.selection_raw(),
        DummyRct, // all flags to false
        multbc.flags(),
        multbc.timestamp(),
        multbc.multFT0PosZ(),
        multbc.multFT0PosZValid(),
        multbc.multV0triggerBits(),
        multbc.multT0triggerBits(),
        multbc.multFDDtriggerBits(),
        multbc.multTriggerMask(),
        multbc.multCollidingBC(),
        multbc.multTVX(),
        multbc.multFV0OrA());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultBCsConverter>(cfgc)};
}
