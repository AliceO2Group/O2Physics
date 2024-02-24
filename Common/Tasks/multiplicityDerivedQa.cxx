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
//
// This code does QA based on a saved derived dataset using the
// tables provided by multiplicityTable.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;

struct MultiplicityDerivedQa {
  // Raw multiplicities
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<int> selection{"sel", 8, "trigger: 7 - sel7, 8 - sel8"};
  Configurable<float> vtxZsel{"vtxZsel", 10, "max vertex Z (cm)"};
  Configurable<bool> INELgtZERO{"INELgtZERO", true, "0 - no, 1 - yes"};
  Configurable<bool> do2Dplots{"do2Dplots", false, "0 - no, 1 - yes"};

  ConfigurableAxis axisMultFV0{"axisMultFV0", {10000, 0, 500000}, "FV0 amplitude"};
  ConfigurableAxis axisMultFT0{"axisMultFT0", {10000, 0, 40000}, "FT0 amplitude"};
  ConfigurableAxis axisMultFT0A{"axisMultFT0A", {10000, 0, 30000}, "FT0A amplitude"};
  ConfigurableAxis axisMultFT0C{"axisMultFT0C", {10000, 0, 10000}, "FT0C amplitude"};
  ConfigurableAxis axisMultFDD{"axisMultFDD", {1000, 0, 4000}, "FDD amplitude"};
  ConfigurableAxis axisMultNTracks{"axisMultNTracks", {500, 0, 500}, "N_{tracks}"};

  ConfigurableAxis axisMultZNA{"axisMultZNA", {1000, 0, 4000}, "ZNA amplitude"};
  ConfigurableAxis axisMultZNC{"axisMultZNC", {1000, 0, 4000}, "ZNC amplitude"};
  ConfigurableAxis axisMultZEM1{"axisMultZEM1", {1000, 0, 4000}, "ZEM1 amplitude"};
  ConfigurableAxis axisMultZEM2{"axisMultZEM2", {1000, 0, 4000}, "ZEM2 amplitude"};
  ConfigurableAxis axisMultZPA{"axisMultZPA", {1000, 0, 4000}, "ZPA amplitude"};
  ConfigurableAxis axisMultZPC{"axisMultZPC", {1000, 0, 4000}, "ZPC amplitude"};

  ConfigurableAxis axisVertexZ{"axisVertexZ", {60, -15, 15}, "Vertex z (cm)"};

  ConfigurableAxis axisContributors{"axisContributors", {100, -0.5f, 99.5f}, "Vertex z (cm)"};
  ConfigurableAxis axisNumberOfPVs{"axisNumberOfPVs", {10, -0.5f, 9.5f}, "Number of reconstructed PVs"};
  ConfigurableAxis axisNchFT0{"axisNchFT0", {500, -0.5f, 499.5f}, "Number of charged particles in FT0 acceptance"};

  void init(InitContext&)
  {
    const AxisSpec axisEvent{10, 0, 10, "Event counter"};

    // Base histograms
    histos.add("multiplicityQa/hEventCounter", "Event counter", kTH1D, {axisEvent});

    histos.add("multiplicityQa/hRawFV0", "Raw FV0", kTH1D, {axisMultFV0});
    histos.add("multiplicityQa/hRawFT0", "Raw FT0", kTH1D, {axisMultFT0});
    histos.add("multiplicityQa/hRawFT0A", "Raw FT0A", kTH1D, {axisMultFT0A});
    histos.add("multiplicityQa/hRawFT0C", "Raw FT0C", kTH1D, {axisMultFT0C});
    histos.add("multiplicityQa/hRawFDD", "Raw FDD", kTH1D, {axisMultFDD});
    histos.add("multiplicityQa/hRawNTracksPV", "Raw NTracks", kTH1D, {axisMultNTracks});

    // ZDC information
    histos.add("multiplicityQa/hRawZNA", "Raw ZNA", kTH1D, {axisMultZNA});
    histos.add("multiplicityQa/hRawZNC", "Raw ZNC", kTH1D, {axisMultZNC});
    histos.add("multiplicityQa/hRawZEM1", "Raw ZEM1", kTH1D, {axisMultZEM1});
    histos.add("multiplicityQa/hRawZEM2", "Raw ZEM2", kTH1D, {axisMultZEM2});
    histos.add("multiplicityQa/hRawZPA", "Raw ZPA", kTH1D, {axisMultZPA});
    histos.add("multiplicityQa/hRawZPC", "Raw ZPC", kTH1D, {axisMultZPC});

    // Vertex-Z profiles for vertex-Z dependency estimate
    histos.add("multiplicityQa/hVtxZFV0A", "Av FV0A vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFT0A", "Av FT0A vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFT0C", "Av FT0C vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFDDA", "Av FDDA vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFDDC", "Av FDDC vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZNTracksPV", "Av NTracks vs vertex Z", kTProfile, {axisVertexZ});

    // two-dimensional histograms
    if (do2Dplots) {
      histos.add("multiplicityQa/h2dNchVsFV0", "FV0", kTH2F, {axisMultFV0, axisMultNTracks});
      histos.add("multiplicityQa/h2dNchVsFT0", "FT0", kTH2F, {axisMultFT0, axisMultNTracks});
      histos.add("multiplicityQa/h2dNchVsFT0A", "FT0A", kTH2F, {axisMultFT0A, axisMultNTracks});
      histos.add("multiplicityQa/h2dNchVsFT0C", "FT0C", kTH2F, {axisMultFT0C, axisMultNTracks});
      histos.add("multiplicityQa/h2dNchVsFDD", "FDD", kTH2F, {axisMultFDD, axisMultNTracks});

      // correlate T0 and V0
      histos.add("multiplicityQa/h2dFT0VsFV0", "FDD", kTH2F, {axisMultFV0, axisMultFT0});

      // correlate ZDC and FT0C
      histos.add("multiplicityQa/h2dZNAVsFT0C", "ZNAvsFT0C", kTH2F, {axisMultFT0C, axisMultZNA});
      histos.add("multiplicityQa/h2dZNCVsFT0C", "ZNPvsFT0C", kTH2F, {axisMultFT0C, axisMultZNC});
      histos.add("multiplicityQa/h2dZEM1VsFT0C", "ZEM1vsFT0C", kTH2F, {axisMultFT0C, axisMultZEM1});
      histos.add("multiplicityQa/h2dZEM2VsFT0C", "ZEM2vsFT0C", kTH2F, {axisMultFT0C, axisMultZEM2});
      histos.add("multiplicityQa/h2dZPAVsFT0C", "ZPAvsFT0C", kTH2F, {axisMultFT0C, axisMultZPA});
      histos.add("multiplicityQa/h2dZPCVsFT0C", "ZPCvsFT0C", kTH2F, {axisMultFT0C, axisMultZPC});
    }
  }

  void process(soa::Join<aod::Mults, aod::MultsExtra>::iterator const& col)
  {
    histos.fill(HIST("multiplicityQa/hEventCounter"), 0.5);
    if (selection == 8 && !col.multSel8()) {
      return;
    }
    if (selection != 8 && selection >= 0) {
      LOGF(fatal, "Unknown selection type!");
    }
    histos.fill(HIST("multiplicityQa/hEventCounter"), 1.5);
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    histos.fill(HIST("multiplicityQa/hEventCounter"), 2.5);

    // Vertex-Z dependencies, necessary for CCDB objects
    histos.fill(HIST("multiplicityQa/hVtxZFV0A"), col.multPVz(), col.multFV0A());
    histos.fill(HIST("multiplicityQa/hVtxZFT0A"), col.multPVz(), col.multFT0A());
    histos.fill(HIST("multiplicityQa/hVtxZFT0C"), col.multPVz(), col.multFT0C());
    histos.fill(HIST("multiplicityQa/hVtxZFDDA"), col.multPVz(), col.multFDDA());
    histos.fill(HIST("multiplicityQa/hVtxZFDDC"), col.multPVz(), col.multFDDC());
    histos.fill(HIST("multiplicityQa/hVtxZNTracksPV"), col.multPVz(), col.multNTracksPV());

    if (fabs(col.multPVz()) > vtxZsel) {
      return;
    }

    histos.fill(HIST("multiplicityQa/hEventCounter"), 3.5);

    LOGF(debug, "multFV0A=%5.0f multFV0C=%5.0f multFV0M=%5.0f multFT0A=%5.0f multFT0C=%5.0f multFT0M=%5.0f multFDDA=%5.0f multFDDC=%5.0f", col.multFV0A(), col.multFV0C(), col.multFV0M(), col.multFT0A(), col.multFT0C(), col.multFT0M(), col.multFDDA(), col.multFDDC());

    // Raw multiplicities
    histos.fill(HIST("multiplicityQa/hRawFV0"), col.multFV0A());
    histos.fill(HIST("multiplicityQa/hRawFT0"), col.multFT0M());
    histos.fill(HIST("multiplicityQa/hRawFT0A"), col.multFT0A());
    histos.fill(HIST("multiplicityQa/hRawFT0C"), col.multFT0C());
    histos.fill(HIST("multiplicityQa/hRawFDD"), col.multFDDM());
    histos.fill(HIST("multiplicityQa/hRawNTracksPV"), col.multNTracksPV());

    histos.fill(HIST("multiplicityQa/hRawZNA"), col.multZNA());
    histos.fill(HIST("multiplicityQa/hRawZNC"), col.multZNC());
    histos.fill(HIST("multiplicityQa/hRawZEM1"), col.multZEM1());
    histos.fill(HIST("multiplicityQa/hRawZEM2"), col.multZEM2());
    histos.fill(HIST("multiplicityQa/hRawZPA"), col.multZPA());
    histos.fill(HIST("multiplicityQa/hRawZPC"), col.multZPC());

    // Profiles
    if (do2Dplots) {
      histos.fill(HIST("multiplicityQa/h2dNchVsFV0"), col.multFV0A(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0"), col.multFT0A() + col.multFT0C(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0A"), col.multFT0A(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0C"), col.multFT0C(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFDD"), col.multFDDA() + col.multFDDC(), col.multNTracksPV());

      histos.fill(HIST("multiplicityQa/h2dZNAVsFT0C"), col.multFT0C(), col.multZNA());
      histos.fill(HIST("multiplicityQa/h2dZNCVsFT0C"), col.multFT0C(), col.multZNC());
      histos.fill(HIST("multiplicityQa/h2dZEM1VsFT0C"), col.multFT0C(), col.multZEM1());
      histos.fill(HIST("multiplicityQa/h2dZEM2VsFT0C"), col.multFT0C(), col.multZEM2());
      histos.fill(HIST("multiplicityQa/h2dZPAVsFT0C"), col.multFT0C(), col.multZPA());
      histos.fill(HIST("multiplicityQa/h2dZPCVsFT0C"), col.multFT0C(), col.multZPC());

      // 2d FT0 vs FV0 fill
      histos.fill(HIST("multiplicityQa/h2dFT0VsFV0"), col.multFV0A(), col.multFT0A() + col.multFT0C());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiplicityDerivedQa>(cfgc)};
}
