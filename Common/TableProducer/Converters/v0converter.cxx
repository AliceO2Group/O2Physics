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
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002

struct V0Converter {
  Produces<aod::V0s_002> v0s_002;

  void process(aod::V0s_001 const& v0s)
  {
    for (auto& v0 : v0s) {
      uint8_t bitMask = static_cast<uint8_t>(1); // first bit on
      v0s_002(v0.collisionId(), v0.posTrackId(), v0.negTrackId(), bitMask);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0Converter>(cfgc)};
}
