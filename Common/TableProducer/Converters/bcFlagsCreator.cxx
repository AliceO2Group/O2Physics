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

// Creates an empty BCFlags for data that doesn't have it to be used seamlessly
// n.b. this will overwrite existing BCFlags, to be discussed if data in mixed condition
struct bcFlagsCreator {
  Produces<aod::BCFlags> bcFlags;

  void process(aod::BCs const& bcTable)
  {
    for (int64_t i = 0; i < bcTable.size(); ++i) {
      bcFlags(0);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<bcFlagsCreator>(cfgc),
  };
}
