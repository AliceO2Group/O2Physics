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
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/GroupSlicer.h"
#include "Framework/StepTHn.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "CommonConstants/MathConstants.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "DataFormatsParameters/GRPObject.h"

#include <TH1F.h>
#include <cmath>
#include <TDirectory.h>
#include <THn.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct TaskHfCorrelations {

  HistogramRegistry registry{"registry"};
  OutputObj<CorrelationContainer> same{"sameEvent"};

  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.5f, "Minimum pT for tracks"};

  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  // Track filters
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt) && ((aod::track::isGlobalTrack == (uint8_t) true) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtended, aod::TrackSelection>>;

  void init(o2::framework::InitContext&)
  {
    registry.add("eventCounter", "eventCounter", {HistType::kTH1F, {{2, 0.5, 2.5}}});
    registry.add("yields", "centrality vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("etaphi", "centrality vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "centrality"}, {100, -2, 2, "#eta"}, {200, 0, 2 * M_PI, "#varphi"}}});
    registry.add("pT", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("eta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("phi", "phi", {HistType::kTH1F, {{100, 0, 2 * M_PI, "#varphi"}}});

    //  set axes of the event counter histogram
    const int nBins = 2;
    std::string labels[nBins];
    labels[0] = "all";
    labels[1] = "after sel8";
    for (int iBin = 0; iBin < nBins; iBin++) {
      registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    
    //  set axes of the correlation container
    std::vector<AxisSpec> axisList = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / centrality"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisEtaEfficiency, "#eta"},
                                      {axisPtEfficiency, "p_{T} (GeV/c)"},
                                      {axisVertexEfficiency, "z-vtx (cm)"}};
    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", axisList));
  }

  template <typename TCollision>
  bool fillCollision(TCollision collision, float centrality)
  {
    registry.fill(HIST("eventCounter"), 1);

    if (!collision.sel8()) {
      return false;
    }

    registry.fill(HIST("eventCounter"), 2);

    return true;
  }

  template <typename TCollision, typename TTracks>
  void fillQA(TCollision collision, float centrality, TTracks tracks)
  {
    for (auto& track1 : tracks) {
      registry.fill(HIST("pT"), track1.pt());
      registry.fill(HIST("eta"), track1.eta());
      registry.fill(HIST("phi"), track1.phi());
      registry.fill(HIST("yields"), centrality, track1.pt(), track1.eta());
      registry.fill(HIST("etaphi"), centrality, track1.eta(), track1.phi());
    }
  }

  template <typename TTarget, typename TTracks>
  void fillCorrelations(TTarget target, TTracks tracks1, TTracks tracks2, float centrality, float posZ)
  {

    auto triggerWeight = 1;
    auto associatedWeight = 1;

    for (auto& track1 : tracks1) {

      //  add getter for NUE trigger efficiency here

      target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), centrality, posZ, triggerWeight);

      for (auto& track2 : tracks2) {

        if (track1 == track2) {
          continue;
        }

        //  add getter for NUE associated efficiency here
    
        //  add pair cuts on phi*

        float deltaPhi = track1.phi() - track2.phi();
        if (deltaPhi > 1.5f * PI) {
          deltaPhi -= TwoPI;
        }
        if (deltaPhi < -PIHalf) {
          deltaPhi += TwoPI;
        }

        target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                    track1.eta() - track2.eta(), track2.pt(), track1.pt(), centrality, deltaPhi, posZ,
                                    triggerWeight * associatedWeight);
      }
    }

  }

  // =====================================
  //    process same event correlations
  // =====================================
  void processSameAOD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aodTracks const& tracks)
  {
    const auto centrality = 50;
    //const auto centrality = collision.centV0M();

    if (fillCollision(collision, centrality) == false) {
      return;
    }
    fillQA(collision, centrality, tracks);
    fillCorrelations(same, tracks, tracks, centrality, collision.posZ());
  }
  PROCESS_SWITCH(TaskHfCorrelations, processSameAOD, "Process same event on AOD", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TaskHfCorrelations>(cfgc)};
}
