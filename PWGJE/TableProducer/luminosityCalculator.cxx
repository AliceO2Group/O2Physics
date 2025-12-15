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

/// \file luminositycalculator.cxx
/// \brief Task to calculate luminosity of measured data. The luminosity per TVX trigger can be obtained from the eventselectiontask (at the time of writing the variable csTVX)
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/DataModel/JetReducedData.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct LuminosityCalculator {

  HistogramRegistry registry;
  void init(InitContext&)
  {
    std::vector<std::string> histLabels = {"BC", "BC+TVX", "BC+TVX+NoTFB", "BC+TVX+NoTFB+NoITSROFB", "Coll", "Coll+TVX", "Coll+TVX+VtxZ+Sel8", "Coll+TVX+VtxZ+Sel8Full", "Coll+TVX+VtxZ+Sel8FullPbPb", "Coll+TVX+VtxZ+SelMC", "Coll+TVX+VtxZ+SelMCFull", "Coll+TVX+VtxZ+SelMCFullPbPb", "Coll+TVX+VtxZ+SelUnanchoredMC", "Coll+TVX+VtxZ+SelTVX", "Coll+TVX+VtxZ+Sel7", "Coll+TVX+VtxZ+Sel7KINT7", "custom"};
    registry.add("counter", "BCs and Collisions", HistType::kTH1D, {{static_cast<int>(histLabels.size()), -0.5, static_cast<double>(histLabels.size()) - 0.5}});
    auto counter = registry.get<TH1>(HIST("counter"));
    for (std::vector<std::string>::size_type iCounter = 0; iCounter < histLabels.size(); iCounter++) {
      counter->GetXaxis()->SetBinLabel(iCounter + 1, histLabels[iCounter].data());
    }
  }

  void processCalculateLuminosity(aod::StoredBCCounts const& bcCounts, aod::StoredCollisionCounts const& collisionCounts)
  {
    int readBC = 0;
    int readBCWithTVXCounter = 0;
    int readBCWithTVXAndNoTFBCounter = 0;
    int readBCWithTVXAndNoTFBAndNoITSROFB = 0;

    for (const auto& bcCount : bcCounts) {
      readBC += bcCount.counts();
      readBCWithTVXCounter += bcCount.countsWithTVX();
      readBCWithTVXAndNoTFBCounter += bcCount.countsWithTVXAndNoTFB();
      readBCWithTVXAndNoTFBAndNoITSROFB += bcCount.countsWithTVXAndNoTFBAndNoITSROFB();
    }

    int readCollision = 0;
    int readCollisionWithTVXCounter = 0;
    int readCollisionWithTVXAndZVertexAndSel8Counter = 0;
    int readCollisionWithTVXAndZVertexAndSel8FullCounter = 0;
    int readCollisionWithTVXAndZVertexAndSel8FullPbPbCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelMCCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelMCFullCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelTVXCounter = 0;
    int readCollisionWithTVXAndZVertexAndSel7Counter = 0;
    int readCollisionWithTVXAndZVertexAndSel7KINT7Counter = 0;
    int readCollisionWithCustomCounter = 0;

    for (const auto& collisionCount : collisionCounts) {
      readCollision += collisionCount.counts();
      readCollisionWithTVXCounter += collisionCount.countsWithTVX();
      readCollisionWithTVXAndZVertexAndSel8Counter += collisionCount.countsWithTVXAndZVertexAndSel8();
      readCollisionWithTVXAndZVertexAndSel8FullCounter += collisionCount.countsWithTVXAndZVertexAndSel8Full();
      readCollisionWithTVXAndZVertexAndSel8FullPbPbCounter += collisionCount.countsWithTVXAndZVertexAndSel8FullPbPb();
      readCollisionWithTVXAndZVertexAndSelMCCounter += collisionCount.countsWithTVXAndZVertexAndSelMC();
      readCollisionWithTVXAndZVertexAndSelMCFullCounter += collisionCount.countsWithTVXAndZVertexAndSelMCFull();
      readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounter += collisionCount.countsWithTVXAndZVertexAndSelMCFullPbPb();
      readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounter += collisionCount.countsWithTVXAndZVertexAndSelUnanchoredMC();
      readCollisionWithTVXAndZVertexAndSelTVXCounter += collisionCount.countsWithTVXAndZVertexAndSelTVX();
      readCollisionWithTVXAndZVertexAndSel7Counter += collisionCount.countsWithTVXAndZVertexAndSel7();
      readCollisionWithTVXAndZVertexAndSel7KINT7Counter += collisionCount.countsWithTVXAndZVertexAndSel7KINT7();
      readCollisionWithCustomCounter += collisionCount.countsWithCustomSelection();
    }

    registry.get<TH1>(HIST("counter"))->SetBinContent(1, registry.get<TH1>(HIST("counter"))->GetBinContent(1) + readBC);
    registry.get<TH1>(HIST("counter"))->SetBinContent(2, registry.get<TH1>(HIST("counter"))->GetBinContent(2) + readBCWithTVXCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(3, registry.get<TH1>(HIST("counter"))->GetBinContent(3) + readBCWithTVXAndNoTFBCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(4, registry.get<TH1>(HIST("counter"))->GetBinContent(4) + readBCWithTVXAndNoTFBAndNoITSROFB);
    registry.get<TH1>(HIST("counter"))->SetBinContent(5, registry.get<TH1>(HIST("counter"))->GetBinContent(5) + readCollision);
    registry.get<TH1>(HIST("counter"))->SetBinContent(6, registry.get<TH1>(HIST("counter"))->GetBinContent(6) + readCollisionWithTVXCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(7, registry.get<TH1>(HIST("counter"))->GetBinContent(7) + readCollisionWithTVXAndZVertexAndSel8Counter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(8, registry.get<TH1>(HIST("counter"))->GetBinContent(8) + readCollisionWithTVXAndZVertexAndSel8FullCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(9, registry.get<TH1>(HIST("counter"))->GetBinContent(9) + readCollisionWithTVXAndZVertexAndSel8FullPbPbCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(10, registry.get<TH1>(HIST("counter"))->GetBinContent(10) + readCollisionWithTVXAndZVertexAndSelMCCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(11, registry.get<TH1>(HIST("counter"))->GetBinContent(11) + readCollisionWithTVXAndZVertexAndSelMCFullCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(12, registry.get<TH1>(HIST("counter"))->GetBinContent(12) + readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(13, registry.get<TH1>(HIST("counter"))->GetBinContent(13) + readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(14, registry.get<TH1>(HIST("counter"))->GetBinContent(14) + readCollisionWithTVXAndZVertexAndSelTVXCounter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(15, registry.get<TH1>(HIST("counter"))->GetBinContent(15) + readCollisionWithTVXAndZVertexAndSel7Counter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(16, registry.get<TH1>(HIST("counter"))->GetBinContent(16) + readCollisionWithTVXAndZVertexAndSel7KINT7Counter);
    registry.get<TH1>(HIST("counter"))->SetBinContent(17, registry.get<TH1>(HIST("counter"))->GetBinContent(17) + readCollisionWithCustomCounter);
  }
  PROCESS_SWITCH(LuminosityCalculator, processCalculateLuminosity, "calculate ingredients for luminosity and fill a histogram", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<LuminosityCalculator>(cfgc, TaskName{"jet-luminosity-calculator"}));

  return WorkflowSpec{tasks};
}
