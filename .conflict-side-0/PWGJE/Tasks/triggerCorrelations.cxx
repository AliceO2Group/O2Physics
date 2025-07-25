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

// QA correlation task for jet trigger
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TriggerCorrelationsTask {

  HistogramRegistry registry;

  std::vector<int> triggerMaskBits;
  void init(o2::framework::InitContext&)
  {
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(jetderiveddatautilities::JTriggerMasks);

    std::vector<std::string> trigSelLabels = {"JetChLowPt", "JetChHighPt", "TrackLowPt", "TrackHighPt", "JetD0ChLowPt", "JetD0ChHighPt", "JetLcChLowPt", "JetLcChHighPt", "EMCALReadout", "JetFullHighPt", "JetFullLowPt", "JetNeutralHighPt", "JetNeutralLowPt", "GammaVeryHighPtEMCAL", "GammaVeryHighPtDCAL", "GammaHighPtEMCAL", "GammaHighPtDCAL", "GammaLowPtEMCAL", "GammaLowPtDCAL", "GammaVeryLowPtEMCAL", "GammaVeryLowPtDCAL"};
    registry.add("triggerCorrelations", "Correlation between jet triggers", HistType::kTH2D, {{static_cast<int>(trigSelLabels.size()), -0.5, static_cast<double>(trigSelLabels.size()) - 0.5, "primary trigger"}, {static_cast<int>(trigSelLabels.size()), -0.5, static_cast<double>(trigSelLabels.size()) - 0.5, "secondary trigger"}});
    auto triggerCorrelation = registry.get<TH2>(HIST("triggerCorrelations"));
    for (std::vector<std::string>::size_type iTrigs = 0; iTrigs < trigSelLabels.size(); iTrigs++) {
      triggerCorrelation->GetXaxis()->SetBinLabel(iTrigs + 1, trigSelLabels[iTrigs].data());
      triggerCorrelation->GetYaxis()->SetBinLabel(iTrigs + 1, trigSelLabels[iTrigs].data());
    }
  }

  template <typename T>
  void fillCorrelationsHistogram(T const& collision, bool fill = false, int iCurrentTrig = -1)
  {
    for (std::vector<int>::size_type iTrig = 0; iTrig < triggerMaskBits.size(); iTrig++) {
      if (fill) {
        if (jetderiveddatautilities::selectTrigger(collision, triggerMaskBits[iTrig])) {
          registry.fill(HIST("triggerCorrelations"), iCurrentTrig, iTrig);
        }
      } else {
        if (jetderiveddatautilities::selectTrigger(collision, triggerMaskBits[iTrig])) {
          fillCorrelationsHistogram(collision, true, iTrig);
        }
      }
    }
  }

  void processTriggeredCorrelations(soa::Join<aod::JCollisions, aod::JChTrigSels, aod::JFullTrigSels, aod::JChHFTrigSels>::iterator const& collision)
  {
    fillCorrelationsHistogram(collision);
  }
  PROCESS_SWITCH(TriggerCorrelationsTask, processTriggeredCorrelations, "QA for trigger correlations", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TriggerCorrelationsTask>(cfgc, TaskName{"trigger-correlations"})};
}
