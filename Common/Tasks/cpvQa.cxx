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
///
/// \brief The task produces basic QA histograms for CPV clusters found in AO2Ds
/// \author Sergey Evdokimov
/// \since 25.10.2022

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <string_view>

using namespace o2;
using namespace o2::framework;

struct cpvQa {
  HistogramRegistry histos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};

  void init(InitContext const&)
  {
    const AxisSpec axisPosX{1000, -100., 100., "x (cm)"};
    const AxisSpec axisPosZ{1000, -100., 100., "z (cm)"};
    const AxisSpec axisAmplitude{20000, 0., 20000., "amplitude"};
    const AxisSpec axisPadMultiplicity{32, 0., 32., "number of pads"};
    const AxisSpec axisNclusters{100, 0., 100., "clusters per event"};

    histos.add("hClustersPerEventTotal", "Total number of clusters per event", kTH1F, {axisNclusters});
    histos.add("hAmplitude", "Cluster total amplitude (all modules)", kTH1F, {axisAmplitude});
    histos.add("hPadMultiplicity", "Number of pads in clusters (all modules)", kTH1F, {axisPadMultiplicity});

    for (int m = 0; m < 3; m++) {
      histos.add(Form("hPosXZ_M%d", m + 2), Form("Cluster position M%d", m + 2), kTH2F, {axisPosX, axisPosZ});
      histos.add(Form("hAmplitude_M%d", m + 2), Form("Cluster total amplitude M%d", m + 2), kTH1F, {axisAmplitude});
      histos.add(Form("hPadMultiplicity_M%d", m + 2), Form("Number of pads in cluster M%d", m + 2), kTH1F, {axisPadMultiplicity});
      histos.add(Form("hClustersPerEvent_M%d", m + 2), Form("Clusters per event M%d", m + 2), kTH1F, {axisNclusters});
    }
  }

  void process(aod::BC const&, aod::CPVClusters const& clusters)
  {
    histos.get<TH1>(HIST("hClustersPerEventTotal"))->Fill(clusters.size());
    if (clusters.size() == 0) {
      return;
    }

    int cluPerEvent[] = {0, 0, 0};
    for (const auto& clu : clusters) {
      static constexpr std::string_view modules[] = {"_M2", "_M3", "_M4"};
      static_for<0, 2>([&](auto mod) {
        constexpr int index = mod.value;
        if (index != clu.moduleNumber() - 2) {
          return;
        }
        cluPerEvent[index]++;
        histos.get<TH2>(HIST("hPosXZ") + HIST(modules[index]))->Fill(clu.posX(), clu.posZ());
        histos.get<TH1>(HIST("hAmplitude") + HIST(modules[index]))->Fill(clu.amplitude());
        histos.get<TH1>(HIST("hPadMultiplicity") + HIST(modules[index]))->Fill(clu.padMult());
      });
      histos.get<TH1>(HIST("hAmplitude"))->Fill(clu.amplitude());
      histos.get<TH1>(HIST("hPadMultiplicity"))->Fill(clu.padMult());
    }

    histos.get<TH1>(HIST("hClustersPerEvent_M2"))->Fill(cluPerEvent[0]);
    histos.get<TH1>(HIST("hClustersPerEvent_M3"))->Fill(cluPerEvent[1]);
    histos.get<TH1>(HIST("hClustersPerEvent_M4"))->Fill(cluPerEvent[2]);
  }
}; // struct cpvQa

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cpvQa>(cfgc, TaskName{"qa-cpv"})};
}
