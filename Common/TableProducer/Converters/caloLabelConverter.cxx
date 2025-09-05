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
#include <vector>

using namespace o2;
using namespace o2::framework;

// Converts the old McCaloLabels_000 table to the new McCaloLabels_001 table where we have a variable size array for associated MCParticles for each calo cell
struct caloLabelConverter {
  Produces<aod::McCaloLabels_001> McCaloLabels_001;

  void process(aod::McCaloLabels_000 const& mccalolabelTable)
  {
    std::vector<float> amplitude = {0};
    std::vector<int32_t> particleId = {0};
    for (auto& mccalolabel : mccalolabelTable) {
      particleId[0] = mccalolabel.mcParticleId();
      // Repopulate new table
      McCaloLabels_001(
        particleId,
        amplitude);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<caloLabelConverter>(cfgc),
  };
}
