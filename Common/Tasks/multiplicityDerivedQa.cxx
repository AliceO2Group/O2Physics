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

#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cmath>

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
  ConfigurableAxis axisMultZEM{"axisMultZEM", {1000, 0, 4000}, "ZEM amplitude"};
  ConfigurableAxis axisMultZPA{"axisMultZPA", {1000, 0, 4000}, "ZPA amplitude"};
  ConfigurableAxis axisMultZPC{"axisMultZPC", {1000, 0, 4000}, "ZPC amplitude"};

  ConfigurableAxis axisVertexZ{"axisVertexZ", {60, -15, 15}, "Vertex z (cm)"};

  ConfigurableAxis axisContributors{"axisContributors", {100, -0.5f, 99.5f}, "Vertex z (cm)"};
  ConfigurableAxis axisNumberOfPVs{"axisNumberOfPVs", {10, -0.5f, 9.5f}, "Number of reconstructed PVs"};
  ConfigurableAxis axisNchFT0{"axisNchFT0", {500, -0.5f, 499.5f}, "Number of charged particles in FT0 acceptance"};

  // artificially select run if desired
  Configurable<int> selectRun{"selectRun", -1, "explicit run selection (-1: no sel)"};

  // Selection criteria for QC studies in 2D plots
  // parameters:
  // --- maxFT0C -> max FT0C value for which this cut will be applied. Nothing done if maxFT0C < 0.0f
  // --- A, B, C  -> A*<variable> + B*FT0C + C > 0, arbitrary linear selection
  static constexpr float default2dCuts[1][4] = {{-1., -1., -1., -1.}};
  Configurable<LabeledArray<float>> selZNA{"selZNA", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selZNA"};
  Configurable<LabeledArray<float>> selZNC{"selZNC", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selZNC"};
  Configurable<LabeledArray<float>> selZPA{"selZPA", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selZPA"};
  Configurable<LabeledArray<float>> selZPC{"selZPC", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selZPC"};
  Configurable<LabeledArray<float>> selZEM1{"selZEM1", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selZEM1"};
  Configurable<LabeledArray<float>> selZEM2{"selZEM2", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selZEM2"};
  Configurable<LabeledArray<float>> selZEM{"selZEM", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selZEM"};
  Configurable<LabeledArray<float>> selFV0A{"selFV0A", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selFV0A"};
  Configurable<LabeledArray<float>> selFT0A{"selFT0A", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selFT0A"};
  Configurable<LabeledArray<float>> selNch{"selNch", {default2dCuts[0], 4, {"maxFT0C", "A", "B", "C"}}, "selNch"};

  ConfigurableAxis axisRawCentrality{"axisRawCentrality", {VARIABLE_WIDTH, 0.000f, 52.320f, 75.400f, 95.719f, 115.364f, 135.211f, 155.791f, 177.504f, 200.686f, 225.641f, 252.645f, 281.906f, 313.850f, 348.302f, 385.732f, 426.307f, 470.146f, 517.555f, 568.899f, 624.177f, 684.021f, 748.734f, 818.078f, 892.577f, 973.087f, 1058.789f, 1150.915f, 1249.319f, 1354.279f, 1465.979f, 1584.790f, 1710.778f, 1844.863f, 1985.746f, 2134.643f, 2291.610f, 2456.943f, 2630.653f, 2813.959f, 3006.631f, 3207.229f, 3417.641f, 3637.318f, 3865.785f, 4104.997f, 4354.938f, 4615.786f, 4885.335f, 5166.555f, 5458.021f, 5762.584f, 6077.881f, 6406.834f, 6746.435f, 7097.958f, 7462.579f, 7839.165f, 8231.629f, 8635.640f, 9052.000f, 9484.268f, 9929.111f, 10389.350f, 10862.059f, 11352.185f, 11856.823f, 12380.371f, 12920.401f, 13476.971f, 14053.087f, 14646.190f, 15258.426f, 15890.617f, 16544.433f, 17218.024f, 17913.465f, 18631.374f, 19374.983f, 20136.700f, 20927.783f, 21746.796f, 22590.880f, 23465.734f, 24372.274f, 25314.351f, 26290.488f, 27300.899f, 28347.512f, 29436.133f, 30567.840f, 31746.818f, 32982.664f, 34276.329f, 35624.859f, 37042.588f, 38546.609f, 40139.742f, 41837.980f, 43679.429f, 45892.130f, 400000.000f}, "raw centrality signal"}; // for QA

  void init(InitContext&)
  {
    const AxisSpec axisCentrality{100, 0, 100, "Centrality"};

    // Base histograms
    histos.add("multiplicityQa/hCentralityAll", "Centrality", kTH1D, {axisCentrality});
    histos.add("multiplicityQa/hCentralitySel8", "Centrality", kTH1D, {axisCentrality});
    histos.add("multiplicityQa/hCentralityINELgtZERO", "Centrality", kTH1D, {axisCentrality});
    histos.add("multiplicityQa/hCentralityVertexZ", "Centrality", kTH1D, {axisCentrality});
    histos.add("multiplicityQa/hCentralitySelected", "Centrality", kTH1D, {axisCentrality});

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
      histos.add("multiplicityQa/h2dNchVsCentrality", "Centrality", kTH2F, {axisCentrality, axisMultNTracks});

      // correlate T0 and V0
      histos.add("multiplicityQa/h2dFT0MVsFV0M", "FDD", kTH2F, {axisMultFV0, axisMultFT0});

      // correlate FIT signals and FT0C
      histos.add("multiplicityQa/h2dFV0AVsFT0C", "FV0AvsFT0C", kTH2F, {axisMultFT0C, axisMultFV0});
      histos.add("multiplicityQa/h2dFT0AVsFT0C", "FT0AvsFT0C", kTH2F, {axisMultFT0C, axisMultFT0A});

      // correlate ZDC and FT0C
      histos.add("multiplicityQa/h2dZNAVsFT0C", "ZNAvsFT0C", kTH2F, {axisMultFT0C, axisMultZNA});
      histos.add("multiplicityQa/h2dZNCVsFT0C", "ZNPvsFT0C", kTH2F, {axisMultFT0C, axisMultZNC});
      histos.add("multiplicityQa/h2dZEM1VsFT0C", "ZEM1vsFT0C", kTH2F, {axisMultFT0C, axisMultZEM1});
      histos.add("multiplicityQa/h2dZEM2VsFT0C", "ZEM2vsFT0C", kTH2F, {axisMultFT0C, axisMultZEM2});
      histos.add("multiplicityQa/h2dZEMVsFT0C", "ZEMvsFT0C", kTH2F, {axisMultFT0C, axisMultZEM});
      histos.add("multiplicityQa/h2dZPAVsFT0C", "ZPAvsFT0C", kTH2F, {axisMultFT0C, axisMultZPA});
      histos.add("multiplicityQa/h2dZPCVsFT0C", "ZPCvsFT0C", kTH2F, {axisMultFT0C, axisMultZPC});

      // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
      // vs centrality
      histos.add("multiplicityQa/h2dFV0AVsCentrality", "FV0AVsCentrality", kTH2F, {axisCentrality, axisMultFV0});
      histos.add("multiplicityQa/h2dFT0AVsCentrality", "FT0AVsCentrality", kTH2F, {axisCentrality, axisMultFT0A});
      histos.add("multiplicityQa/h2dZNAVsCentrality", "ZNAVsCentrality", kTH2F, {axisCentrality, axisMultZNA});
      histos.add("multiplicityQa/h2dZNCVsCentrality", "ZNPVsCentrality", kTH2F, {axisCentrality, axisMultZNC});
      histos.add("multiplicityQa/h2dZEM1VsCentrality", "ZEM1VsCentrality", kTH2F, {axisCentrality, axisMultZEM1});
      histos.add("multiplicityQa/h2dZEM2VsCentrality", "ZEM2VsCentrality", kTH2F, {axisCentrality, axisMultZEM2});
      histos.add("multiplicityQa/h2dZEMVsCentrality", "ZEMVsCentrality", kTH2F, {axisCentrality, axisMultZEM});
      histos.add("multiplicityQa/h2dZPAVsCentrality", "ZPAVsCentrality", kTH2F, {axisCentrality, axisMultZPA});
      histos.add("multiplicityQa/h2dZPCVsCentrality", "ZPCVsCentrality", kTH2F, {axisCentrality, axisMultZPC});
      // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
    }

    // for QA and test purposes
    auto hRawCentrality = histos.add<TH1>("hRawCentrality", "hRawCentrality", kTH1F, {axisRawCentrality});

    for (int ii = 1; ii < 101; ii++) {
      float value = 100.5f - static_cast<float>(ii);
      hRawCentrality->SetBinContent(ii, value);
    }
  }

  void process(soa::Join<aod::Mults, aod::MultsExtra>::iterator const& col)
  {
    // for analysis of slim derived data (MC, etc) and specific run selection
    if (selectRun > -1 && col.multRunNumber() != selectRun)
      return;

    auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
    float centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(col.multFT0C()));
    histos.fill(HIST("multiplicityQa/hCentralityAll"), centrality);

    if (selection == 8 && !col.multSel8()) {
      return;
    }
    if (selection != 8 && selection >= 0) {
      LOGF(fatal, "Unknown selection type!");
    }
    histos.fill(HIST("multiplicityQa/hCentralitySel8"), centrality);
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    histos.fill(HIST("multiplicityQa/hCentralityINELgtZERO"), centrality);

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
    histos.fill(HIST("multiplicityQa/hCentralityVertexZ"), centrality);

    // apply special event selections
    if (selZNA->get("maxFT0C") > -0.5f && col.multFT0C() < selZNA->get("maxFT0C") && (selZNA->get("A") * col.multZNA() + selZNA->get("B") * col.multFT0C() + selZNA->get("C") < 0.0f))
      return;
    if (selZNC->get("maxFT0C") > -0.5f && col.multFT0C() < selZNC->get("maxFT0C") && (selZNC->get("A") * col.multZNC() + selZNC->get("B") * col.multFT0C() + selZNC->get("C") < 0.0f))
      return;
    if (selZPA->get("maxFT0C") > -0.5f && col.multFT0C() < selZPA->get("maxFT0C") && (selZPA->get("A") * col.multZPA() + selZPA->get("B") * col.multFT0C() + selZPA->get("C") < 0.0f))
      return;
    if (selZPC->get("maxFT0C") > -0.5f && col.multFT0C() < selZPC->get("maxFT0C") && (selZPC->get("A") * col.multZPC() + selZPC->get("B") * col.multFT0C() + selZPC->get("C") < 0.0f))
      return;
    if (selZEM1->get("maxFT0C") > -0.5f && col.multFT0C() < selZEM1->get("maxFT0C") && (selZEM1->get("A") * col.multZEM1() + selZEM1->get("B") * col.multFT0C() + selZEM1->get("C") < 0.0f))
      return;
    if (selZEM2->get("maxFT0C") > -0.5f && col.multFT0C() < selZEM2->get("maxFT0C") && (selZEM2->get("A") * col.multZEM2() + selZEM2->get("B") * col.multFT0C() + selZEM2->get("C") < 0.0f))
      return;
    if (selZEM->get("maxFT0C") > -0.5f && col.multFT0C() < selZEM->get("maxFT0C") && (selZEM->get("A") * col.multZEM2() + selZEM->get("B") * col.multFT0C() + selZEM->get("C") < 0.0f))
      return;
    if (selFV0A->get("maxFT0C") > -0.5f && col.multFT0C() < selFV0A->get("maxFT0C") && (selFV0A->get("A") * col.multFV0A() + selFV0A->get("B") * col.multFT0C() + selFV0A->get("C") < 0.0f))
      return;
    if (selFT0A->get("maxFT0C") > -0.5f && col.multFT0C() < selFT0A->get("maxFT0C") && (selFT0A->get("A") * col.multFT0A() + selFT0A->get("B") * col.multFT0C() + selFT0A->get("C") < 0.0f))
      return;
    if (selNch->get("maxFT0C") > -0.5f && col.multFT0C() < selNch->get("maxFT0C") && (selNch->get("A") * col.multNTracksPV() + selNch->get("B") * col.multFT0C() + selNch->get("C") < 0.0f))
      return;

    histos.fill(HIST("multiplicityQa/hCentralitySelected"), centrality);

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
      histos.fill(HIST("multiplicityQa/h2dNchVsCentrality"), centrality, col.multNTracksPV());

      // correlate FIT signals and FT0C
      histos.fill(HIST("multiplicityQa/h2dFV0AVsFT0C"), col.multFT0C(), col.multFV0A());
      histos.fill(HIST("multiplicityQa/h2dFT0AVsFT0C"), col.multFT0C(), col.multFT0A());

      // correlate ZDC and FT0C
      histos.fill(HIST("multiplicityQa/h2dZNAVsFT0C"), col.multFT0C(), col.multZNA());
      histos.fill(HIST("multiplicityQa/h2dZNCVsFT0C"), col.multFT0C(), col.multZNC());
      histos.fill(HIST("multiplicityQa/h2dZEM1VsFT0C"), col.multFT0C(), col.multZEM1());
      histos.fill(HIST("multiplicityQa/h2dZEM2VsFT0C"), col.multFT0C(), col.multZEM2());
      histos.fill(HIST("multiplicityQa/h2dZEMVsFT0C"), col.multFT0C(), col.multZEM1() + col.multZEM2());
      histos.fill(HIST("multiplicityQa/h2dZPAVsFT0C"), col.multFT0C(), col.multZPA());
      histos.fill(HIST("multiplicityQa/h2dZPCVsFT0C"), col.multFT0C(), col.multZPC());

      // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
      // vs centrality
      histos.fill(HIST("multiplicityQa/h2dFV0AVsCentrality"), centrality, col.multFV0A());
      histos.fill(HIST("multiplicityQa/h2dFT0AVsCentrality"), centrality, col.multFT0A());
      histos.fill(HIST("multiplicityQa/h2dZNAVsCentrality"), centrality, col.multZNA());
      histos.fill(HIST("multiplicityQa/h2dZNCVsCentrality"), centrality, col.multZNC());
      histos.fill(HIST("multiplicityQa/h2dZEM1VsCentrality"), centrality, col.multZEM1());
      histos.fill(HIST("multiplicityQa/h2dZEM2VsCentrality"), centrality, col.multZEM2());
      histos.fill(HIST("multiplicityQa/h2dZEMVsCentrality"), centrality, col.multZEM1() + col.multZEM2());
      histos.fill(HIST("multiplicityQa/h2dZPAVsCentrality"), centrality, col.multZPA());
      histos.fill(HIST("multiplicityQa/h2dZPCVsCentrality"), centrality, col.multZPC());
      // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+

      // 2d FT0 vs FV0 fill
      histos.fill(HIST("multiplicityQa/h2dFT0MVsFV0M"), col.multFV0A(), col.multFT0A() + col.multFT0C());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiplicityDerivedQa>(cfgc)};
}
