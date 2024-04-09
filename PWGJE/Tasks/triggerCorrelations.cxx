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

#include <string>

#include "EMCALBase/Geometry.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/CCDB/TriggerAliases.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct TriggerCorrelationsTask {

  HistogramRegistry registry;

  int nChTrigs;
  int nFullTrigs;
  int nChHFTrigs;
  int nAllTrigs;

  int chTrigOffset;
  int fullTrigOffset;
  int chHFTrigOffset;

  void init(o2::framework::InitContext&)
  {
    std::vector<std::string> trigSelChLabels = {"chargedLow", "chargedHigh", "trackPt"};
    std::vector<std::string> trigSelFullLabels = {"fullHigh", "fullLow", "neutralHigh", "neutralLow", "gammaVeryHighEMCAL", "gammaHighEMCAL", "gammaLowEMCAL", "gammaVeryLowEMCAL", "gammaVeryHighDCAL", "gammaHighDCAL", "gammaLowDCAL", "gammaVeryLowDCAL"};
    std::vector<std::string> trigSelChHFLabels = {"chargedD0Low", "chargedD0High", "chargedLcLow", "chargedLcHigh"};
    nChTrigs = trigSelChLabels.size();
    nFullTrigs = trigSelFullLabels.size();
    nChHFTrigs = trigSelChHFLabels.size();
    nAllTrigs = nChTrigs + nFullTrigs + nChHFTrigs;
    chTrigOffset = 0;
    fullTrigOffset = chTrigOffset + nChTrigs;
    chHFTrigOffset = fullTrigOffset + nFullTrigs;
    registry.add("triggerCorrelations", "Correlation between jet triggers", HistType::kTH2D, {{nAllTrigs, -0.5, static_cast<double>(nAllTrigs) - 0.5, "primary trigger"}, {nAllTrigs, -0.5, static_cast<double>(nAllTrigs) - 0.5, "secondary trigger"}});
    auto triggerCorrelation = registry.get<TH2>(HIST("triggerCorrelations"));
    for (auto iChTrigs = 0; iChTrigs < nChTrigs; iChTrigs++) {
      triggerCorrelation->GetXaxis()->SetBinLabel(iChTrigs + chTrigOffset + 1, trigSelChLabels[iChTrigs].data());
      triggerCorrelation->GetYaxis()->SetBinLabel(iChTrigs + chTrigOffset + 1, trigSelChLabels[iChTrigs].data());
    }
    for (auto iFullTrigs = 0; iFullTrigs < nFullTrigs; iFullTrigs++) {
      triggerCorrelation->GetXaxis()->SetBinLabel(iFullTrigs + fullTrigOffset + 1, trigSelFullLabels[iFullTrigs].data());
      triggerCorrelation->GetYaxis()->SetBinLabel(iFullTrigs + fullTrigOffset + 1, trigSelFullLabels[iFullTrigs].data());
    }
    for (auto iChHFTrigs = 0; iChHFTrigs < nChHFTrigs; iChHFTrigs++) {
      triggerCorrelation->GetXaxis()->SetBinLabel(iChHFTrigs + chHFTrigOffset + 1, trigSelChHFLabels[iChHFTrigs].data());
      triggerCorrelation->GetYaxis()->SetBinLabel(iChHFTrigs + chHFTrigOffset + 1, trigSelChHFLabels[iChHFTrigs].data());
    }
  }

  template <typename T>
  void fillCorrelationsHistogram(T const& collision, bool fill = false, int iTrig = -1)
  {

    for (auto iChTrigs = 0; iChTrigs < nChTrigs; iChTrigs++) {
      if (fill) {
        if (jetderiveddatautilities::selectChargedTrigger(collision, iChTrigs + 1)) {
          registry.fill(HIST("triggerCorrelations"), iTrig, iChTrigs + chTrigOffset);
        }
      } else {
        if (jetderiveddatautilities::selectChargedTrigger(collision, iChTrigs + 1)) {
          fillCorrelationsHistogram(collision, true, iChTrigs + chTrigOffset);
        }
      }
    }

    for (auto iFullTrigs = 0; iFullTrigs < nFullTrigs; iFullTrigs++) {
      if (fill) {
        if (jetderiveddatautilities::selectFullTrigger(collision, iFullTrigs + 1)) {
          registry.fill(HIST("triggerCorrelations"), iTrig, iFullTrigs + fullTrigOffset);
        }
      } else {
        if (jetderiveddatautilities::selectFullTrigger(collision, iFullTrigs + 1)) {
          fillCorrelationsHistogram(collision, true, iFullTrigs + fullTrigOffset);
        }
      }
    }

    for (auto iChHFTrigs = 0; iChHFTrigs < nFullTrigs; iChHFTrigs++) {
      if (fill) {
        if (jetderiveddatautilities::selectChargedHFTrigger(collision, iChHFTrigs + 1)) {
          registry.fill(HIST("triggerCorrelations"), iTrig, iChHFTrigs + chHFTrigOffset);
        }
      } else {
        if (jetderiveddatautilities::selectChargedHFTrigger(collision, iChHFTrigs + 1)) {
          fillCorrelationsHistogram(collision, true, iChHFTrigs + chHFTrigOffset);
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
