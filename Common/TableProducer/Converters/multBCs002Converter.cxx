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

/// \file multBCs002Converter.cxx
/// \brief Converts MultBCs table from version 001 to 002
/// \author Jesper Karlsson Gumrpecht

#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

struct MultBCs002Converter {
  o2::framework::Produces<o2::aod::MultBCs_002> multBC;
  static constexpr float DummyValue = -1.f;
  void process(o2::aod::MultBCs_001 const& multBCs)
  {
    multBC.reserve(multBCs.size());
    for (const auto& multbc : multBCs) {
      multBC(multbc.multFT0A(),
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
             multbc.multFV0AOuter(),
             multbc.multFT0AOuter(),
             DummyValue /* dummy amplitude for FT0C Outer */);
    }
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{adaptAnalysisTask<MultBCs002Converter>(cfgc)};
}
