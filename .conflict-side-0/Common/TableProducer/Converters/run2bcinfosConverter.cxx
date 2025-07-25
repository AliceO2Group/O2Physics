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

struct Run2BCInfosConverter {
  Produces<aod::Run2BCInfos_001> Run2BCInfos_001;
  void process(aod::Run2BCInfos_000 const& Run2BCInfos_000)
  {

    for (const auto& entry : Run2BCInfos_000) {
      Run2BCInfos_001(entry.eventCuts(),
                      entry.triggerMaskNext50(), entry.l0TriggerInputMask(),
                      entry.spdClustersL0(), entry.spdClustersL1(),
                      entry.spdFiredChipsL0(), entry.spdFiredChipsL1(),
                      entry.spdFiredFastOrL0(), entry.spdFiredFastOrL1(),
                      entry.v0TriggerChargeA(), entry.v0TriggerChargeC(),
                      0, 0);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Run2BCInfosConverter>(cfgc)};
}
