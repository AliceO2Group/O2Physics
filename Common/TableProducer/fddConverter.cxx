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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;

// Converts FDD table from version 000 to 001

struct FddConverter {
  Produces<aod::FDDs_001> fdd_001;

  void process(aod::FDDs_000 const& fdd_000)
  {
    for (auto& p : fdd_000) {
      int16_t amplitudeA_001[8] = {0u};
      int16_t amplitudeC_001[8] = {0u};

      for (int i = 0; i < 4; i++) {
        amplitudeA_001[i] = p.amplitudeA_000()[i];
        amplitudeA_001[i+4] = p.amplitudeA_000()[i];
	
	amplitudeC_001[i] = p.amplitudeC_000()[i];
        amplitudeC_001[i+4] = p.amplitudeC_000()[i];
        }

      fdd_001(p.bcId(), amplitudeA_001, amplitudeC_001,
              p.timeA(), p.timeC(), p.triggerMask());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FddConverter>(cfgc),
  };
}
