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

/// \file longRangeDihadronCor.cxx
/// \brief long range di-hadron correlation for O-O, Pb-Pb collisions
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Sep/10/2025

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

#include "TF1.h"
#include "TRandom3.h"
#include <TPDGCode.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct LongRangeDihadronCor {
  Service<ccdb::BasicCCDBManager> ccdb;

  O2_DEFINE_CONFIGURABLE(cfgCutVtxZ, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")
  O2_DEFINE_CONFIGURABLE(cfgSelCollByNch, bool, true, "Select collisions by Nch or centrality")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMin, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMax, int, 10, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMin, float, 60.0f, "Minimum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMax, float, 80.0f, "Maximum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgMixEventNumMin, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgCentTableUnavailable, bool, false, "if a dataset does not provide centrality information")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoSameBunchPileup, bool, false, "rejects collisions which are associated with the same found-by-T0 bunch crossing")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoITSROFrameBorder, bool, false, "reject events at ITS ROF border")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoTimeFrameBorder, bool, false, "reject events at TF border")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodZvtxFT0vsPV, bool, false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInTimeRangeStandard, bool, false, "no collisions in specified time range")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayersAll, bool, true, "cut time intervals with dead ITS staves")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInRofStandard, bool, false, "no other collisions in this Readout Frame with per-collision multiplicity above threshold")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoHighMultCollInPrevRof, bool, false, "veto an event if FT0C amplitude in previous ITS ROF is above threshold")
  O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelV0AT0ACut, bool, true, "V0A T0A 5 sigma cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelOccupancy, bool, true, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 2000, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgSwitchCor, int, 0, "0:TPC-FT0A; 1:TPC-FT0C")
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgMultCentHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 10.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultCentLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultT0CCutEnabled, bool, false, "Enable Global multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultT0CCutPars{"cfgMultT0CCutPars", std::vector<double>{143.04, -4.58368, 0.0766055, -0.000727796, 2.86153e-06, 23.3108, -0.36304, 0.00437706, -4.717e-05, 1.98332e-07}, "Global multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultPVT0CCutEnabled, bool, false, "Enable PV multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultPVT0CCutPars{"cfgMultPVT0CCutPars", std::vector<double>{195.357, -6.15194, 0.101313, -0.000955828, 3.74793e-06, 30.0326, -0.43322, 0.00476265, -5.11206e-05, 2.13613e-07}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVHighCutFunction, std::string, "[0]+[1]*x + 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVLowCutFunction, std::string, "[0]+[1]*x - 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCutEnabled, bool, false, "Enable global multiplicity vs PV multiplicity cut")
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.140809, 0.734344, 2.77495, 0.0165935}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0AHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 4.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ALowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ACutEnabled, bool, false, "Enable global multiplicity vs V0A multiplicity cut")
    Configurable<std::vector<double>> cfgMultMultV0ACutPars{"cfgMultMultV0ACutPars", std::vector<double>{534.893, 184.344, 0.423539, -0.00331436, 5.34622e-06, 871.239, 53.3735, -0.203528, 0.000122758, 5.41027e-07}, "Global multiplicity vs V0A multiplicity cut parameter values"};
    std::vector<double> multT0CCutPars;
    std::vector<double> multPVT0CCutPars;
    std::vector<double> multGlobalPVCutPars;
    std::vector<double> multMultV0ACutPars;
    TF1* fMultPVT0CCutLow = nullptr;
    TF1* fMultPVT0CCutHigh = nullptr;
    TF1* fMultT0CCutLow = nullptr;
    TF1* fMultT0CCutHigh = nullptr;
    TF1* fMultGlobalPVCutLow = nullptr;
    TF1* fMultGlobalPVCutHigh = nullptr;
    TF1* fMultMultV0ACutLow = nullptr;
    TF1* fMultMultV0ACutHigh = nullptr;
    TF1* fT0AV0AMean = nullptr;
    TF1* fT0AV0ASigma = nullptr;
  } cfgFuncParas;

  SliceCache cache;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {32, -5.8, -2.6}, "delta eta axis, -5.8~-2.6 for TPC-FT0A, 1.2~4.2 for TPC-FT0C"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt associated axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisAmplitudeFt0a{"axisAmplitudeFt0a", {5000, 0, 1000}, "FT0A amplitude"};
  ConfigurableAxis axisChannelFt0aAxis{"axisChannelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};

  // make the filters and cuts.
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVtxZ);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  // FT0 geometry
  o2::ft0::Geometry ft0Det;
  std::vector<o2::detectors::AlignParam>* offsetFT0;

  // Corrections
  TH3D* mEfficiency = nullptr;
  TH1D* mCentralityWeight = nullptr;
  bool correctionsLoaded = false;

  // Define the outputs
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};
  HistogramRegistry registry{"registry"};

  // define global variables
  TRandom3* gRandom = new TRandom3();
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    // Count the total number of enum
    kCount_CentEstimators
  };
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };
  enum FITIndex {
    kFT0A = 0,
    kFT0C = 1
  };

  void init(InitContext&)
  {
    if (cfgCentTableUnavailable && !cfgSelCollByNch) {
      LOGF(fatal, "Centrality table is unavailable, cannot select collisions by centrality");
    }
    const AxisSpec axisPhi{72, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    // Event Counter
    if (doprocessSame && cfgUseAdditionalEventCut) {
      registry.add("hEventCountSpecific", "Number of Event;; Count", {HistType::kTH1D, {{12, 0, 12}}});
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(1, "after sel8");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(2, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(3, "kNoITSROFrameBorder");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(4, "kNoTimeFrameBorder");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(5, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(7, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(8, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(9, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(10, "occupancy");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(11, "MultCorrelation");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(12, "cfgEvSelV0AT0ACut");
    }

    if (cfgEvSelMultCorrelation) {
      cfgFuncParas.multT0CCutPars = cfgFuncParas.cfgMultT0CCutPars;
      cfgFuncParas.multPVT0CCutPars = cfgFuncParas.cfgMultPVT0CCutPars;
      cfgFuncParas.multGlobalPVCutPars = cfgFuncParas.cfgMultGlobalPVCutPars;
      cfgFuncParas.multMultV0ACutPars = cfgFuncParas.cfgMultMultV0ACutPars;
      cfgFuncParas.fMultPVT0CCutLow = new TF1("fMultPVT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutLow->SetParameters(&(cfgFuncParas.multPVT0CCutPars[0]));
      cfgFuncParas.fMultPVT0CCutHigh = new TF1("fMultPVT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutHigh->SetParameters(&(cfgFuncParas.multPVT0CCutPars[0]));

      cfgFuncParas.fMultT0CCutLow = new TF1("fMultT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutLow->SetParameters(&(cfgFuncParas.multT0CCutPars[0]));
      cfgFuncParas.fMultT0CCutHigh = new TF1("fMultT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutHigh->SetParameters(&(cfgFuncParas.multT0CCutPars[0]));

      cfgFuncParas.fMultGlobalPVCutLow = new TF1("fMultGlobalPVCutLow", cfgFuncParas.cfgMultMultPVLowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutLow->SetParameters(&(cfgFuncParas.multGlobalPVCutPars[0]));
      cfgFuncParas.fMultGlobalPVCutHigh = new TF1("fMultGlobalPVCutHigh", cfgFuncParas.cfgMultMultPVHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutHigh->SetParameters(&(cfgFuncParas.multGlobalPVCutPars[0]));

      cfgFuncParas.fMultMultV0ACutLow = new TF1("fMultMultV0ACutLow", cfgFuncParas.cfgMultMultV0ALowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutLow->SetParameters(&(cfgFuncParas.multMultV0ACutPars[0]));
      cfgFuncParas.fMultMultV0ACutHigh = new TF1("fMultMultV0ACutHigh", cfgFuncParas.cfgMultMultV0AHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutHigh->SetParameters(&(cfgFuncParas.multMultV0ACutPars[0]));

      cfgFuncParas.fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      cfgFuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      cfgFuncParas.fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      cfgFuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(cfgCentEstimator);
    // Make histograms to check the distributions after cuts
    if (doprocessSame) {
      registry.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
      registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
      registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
      registry.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
      registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
      registry.add("pTCorrected", "pTCorrected", {HistType::kTH1D, {axisPtTrigger}});
      registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      registry.add("Nch_used", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}}); // histogram to see how many events are in the same and mixed event
      registry.add("Centrality", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}});
      registry.add("CentralityWeighted", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}});
      registry.add("Centrality_used", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}}); // histogram to see how many events are in the same and mixed event
      registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      registry.add("zVtx_used", "zVtx_used", {HistType::kTH1D, {axisVertex}});
      registry.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      registry.add("Assoc_amp_same", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Assoc_amp_mixed", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
    }

    registry.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");
    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {
      {axisEtaEfficiency, "#eta"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisVertexEfficiency, "z-vtx (cm)"},
    };
    std::vector<AxisSpec> userAxis;

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));

    LOGF(info, "End of init");
  }

  double getPhiFT0(uint64_t chno, int i)
  {
    // offsetFT0[0]: FT0A, offsetFT0[1]: FT0C
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return RecoDecay::phi(chPos.X() + (*offsetFT0)[i].getX(), chPos.Y() + (*offsetFT0)[i].getY());
  }

  double getEtaFT0(uint64_t chno, int i)
  {
    // offsetFT0[0]: FT0A, offsetFT0[1]: FT0C
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  template <typename TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent;
    switch (cfgCentEstimator) {
      case kCentFT0C:
        cent = collision.centFT0C();
        break;
      case kCentFT0CVariant1:
        cent = collision.centFT0CVariant1();
        break;
      case kCentFT0M:
        cent = collision.centFT0M();
        break;
      case kCentFV0A:
        cent = collision.centFV0A();
        break;
      default:
        cent = collision.centFT0C();
    }
    return cent;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    return ((track.tpcNClsFound() >= cfgCutTPCclu) && (track.tpcNClsCrossedRows() >= cfgCutTPCCrossedRows) && (track.itsNCls() >= cfgCutITSclu));
  }

  void loadAlignParam(uint64_t timestamp)
  {
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    if (offsetFT0 == nullptr) {
      LOGF(fatal, "Could not load FT0/Calib/Align for timestamp %d", timestamp);
    }
  }

  void loadCorrection(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (cfgEfficiency.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiency.value.c_str(), "READ");
        mEfficiency = reinterpret_cast<TH3D*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        mEfficiency = ccdb->getForTimeStamp<TH3D>(cfgEfficiency, timestamp);
      }
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    if (cfgCentralityWeight.value.empty() == false) {
      mCentralityWeight = ccdb->getForTimeStamp<TH1D>(cfgCentralityWeight, timestamp);
      if (mCentralityWeight == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgCentralityWeight.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgCentralityWeight.value.c_str(), (void*)mCentralityWeight);
    }
    correctionsLoaded = true;
  }

  bool getEfficiencyCorrection(float& weight_nue, float eta, float pt, float posZ)
  {
    float eff = 1.;
    if (mEfficiency) {
      int etaBin = mEfficiency->GetXaxis()->FindBin(eta);
      int ptBin = mEfficiency->GetYaxis()->FindBin(pt);
      int zBin = mEfficiency->GetZaxis()->FindBin(posZ);
      eff = mEfficiency->GetBinContent(etaBin, ptBin, zBin);
    } else {
      eff = 1.0;
    }
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    return true;
  }

  bool getCentralityWeight(float& weightCent, const float centrality)
  {
    float weight = 1.;
    if (mCentralityWeight)
      weight = mCentralityWeight->GetBinContent(mCentralityWeight->FindBin(centrality));
    else
      weight = 1.0;
    if (weight == 0)
      return false;
    weightCent = weight;
    return true;
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    float weff1 = 1;
    float vtxz = collision.posZ();
    for (auto const& track1 : tracks) {
      if (!trackSelected(track1))
        continue;
      if (!getEfficiencyCorrection(weff1, track1.eta(), track1.pt(), vtxz))
        continue;
      registry.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      registry.fill(HIST("Eta"), track1.eta());
      registry.fill(HIST("EtaCorrected"), track1.eta(), weff1);
      registry.fill(HIST("pT"), track1.pt());
      registry.fill(HIST("pTCorrected"), track1.pt(), weff1);
    }
  }

  template <typename TFT0s>
  void getChannel(TFT0s const& ft0, std::size_t const& iCh, int& id, float& ampl)
  {
    int switchCor = cfgSwitchCor;
    if (switchCor == kFT0C) {
      id = ft0.channelC()[iCh];
      ampl = ft0.amplitudeC()[iCh];
    } else if (switchCor == kFT0A) {
      id = ft0.channelA()[iCh];
      ampl = ft0.amplitudeA()[iCh];
    } else {
      LOGF(fatal, "Cor Index %d out of range", switchCor);
    }
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TFT0s>
  void fillCorrelations(TTracks tracks1, TFT0s const& ft0, float posZ, int system, float cent, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    if (system == SameEvent) {
      if (!cfgCentTableUnavailable)
        registry.fill(HIST("Centrality_used"), cent);
      registry.fill(HIST("Nch_used"), tracks1.size());
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;
      if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
        continue;
      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt(), eventWeight * triggerWeight);
      }

      std::size_t channelSize = 0;
      int switchCor = cfgSwitchCor;
      if (switchCor == kFT0C) {
        channelSize = ft0.channelC().size();
      } else if (switchCor == kFT0A) {
        channelSize = ft0.channelA().size();
      } else {
        LOGF(fatal, "Cor Index %d out of range", switchCor);
      }
      for (std::size_t iCh = 0; iCh < channelSize; iCh++) {
        int chanelid = 0;
        float ampl = 0.;
        getChannel(ft0, iCh, chanelid, ampl);
        auto phi = getPhiFT0(chanelid, switchCor);
        auto eta = getEtaFT0(chanelid, switchCor);
        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - phi, -PIHalf);
        float deltaEta = track1.eta() - eta;
        // fill the right sparse and histograms
        if (system == SameEvent) {
          registry.fill(HIST("Assoc_amp_same"), chanelid, ampl);
          same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track1.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight);
          registry.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight);
        } else if (system == MixedEvent) {
          registry.fill(HIST("Assoc_amp_mixed"), chanelid, ampl);
          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track1.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight);
          registry.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight);
        }
      }
    }
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int multTrk, const float centrality, const bool fillCounter)
  {
    registry.fill(HIST("hEventCountSpecific"), 0.5);
    if (cfgEvSelkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (fillCounter && cfgEvSelkNoSameBunchPileup)
      registry.fill(HIST("hEventCountSpecific"), 1.5);
    if (cfgEvSelkNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgEvSelkNoITSROFrameBorder)
      registry.fill(HIST("hEventCountSpecific"), 2.5);
    if (cfgEvSelkNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgEvSelkNoTimeFrameBorder)
      registry.fill(HIST("hEventCountSpecific"), 3.5);
    if (cfgEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (fillCounter && cfgEvSelkIsGoodZvtxFT0vsPV)
      registry.fill(HIST("hEventCountSpecific"), 4.5);
    if (cfgEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (fillCounter && cfgEvSelkNoCollInTimeRangeStandard)
      registry.fill(HIST("hEventCountSpecific"), 5.5);
    if (cfgEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (fillCounter && cfgEvSelkIsGoodITSLayersAll)
      registry.fill(HIST("hEventCountSpecific"), 6.5);
    if (cfgEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (fillCounter && cfgEvSelkNoCollInRofStandard)
      registry.fill(HIST("hEventCountSpecific"), 7.5);
    if (cfgEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (fillCounter && cfgEvSelkNoHighMultCollInPrevRof)
      registry.fill(HIST("hEventCountSpecific"), 8.5);
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgEvSelOccupancy && (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh))
      return 0;
    if (fillCounter && cfgEvSelOccupancy)
      registry.fill(HIST("hEventCountSpecific"), 9.5);

    auto multNTracksPV = collision.multNTracksPV();
    if (cfgEvSelMultCorrelation) {
      if (cfgFuncParas.cfgMultPVT0CCutEnabled && !cfgCentTableUnavailable) {
        if (multNTracksPV < cfgFuncParas.fMultPVT0CCutLow->Eval(centrality))
          return 0;
        if (multNTracksPV > cfgFuncParas.fMultPVT0CCutHigh->Eval(centrality))
          return 0;
      }
      if (cfgFuncParas.cfgMultT0CCutEnabled && !cfgCentTableUnavailable) {
        if (multTrk < cfgFuncParas.fMultT0CCutLow->Eval(centrality))
          return 0;
        if (multTrk > cfgFuncParas.fMultT0CCutHigh->Eval(centrality))
          return 0;
      }
      if (cfgFuncParas.cfgMultGlobalPVCutEnabled) {
        if (multTrk < cfgFuncParas.fMultGlobalPVCutLow->Eval(multNTracksPV))
          return 0;
        if (multTrk > cfgFuncParas.fMultGlobalPVCutHigh->Eval(multNTracksPV))
          return 0;
      }
      if (cfgFuncParas.cfgMultMultV0ACutEnabled) {
        if (collision.multFV0A() < cfgFuncParas.fMultMultV0ACutLow->Eval(multTrk))
          return 0;
        if (collision.multFV0A() > cfgFuncParas.fMultMultV0ACutHigh->Eval(multTrk))
          return 0;
      }
    }
    if (fillCounter && cfgEvSelMultCorrelation)
      registry.fill(HIST("hEventCountSpecific"), 10.5);

    // V0A T0A 5 sigma cut
    float sigma = 5.0;
    if (cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > sigma * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (fillCounter && cfgEvSelV0AT0ACut)
      registry.fill(HIST("hEventCountSpecific"), 11.5);

    return 1;
  }

  void processSame(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = -1.;
    float weightCent = 1.0f;
    if (!cfgCentTableUnavailable) {
      cent = getCentrality(collision);
    }
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent, true))
      return;
    if (!collision.has_foundFT0())
      return;
    loadAlignParam(bc.timestamp());
    loadCorrection(bc.timestamp());
    if (!cfgCentTableUnavailable) {
      getCentralityWeight(weightCent, cent);
      registry.fill(HIST("Centrality"), cent);
      registry.fill(HIST("CentralityWeighted"), cent, weightCent);
    }
    registry.fill(HIST("Nch"), tracks.size());
    registry.fill(HIST("zVtx"), collision.posZ());

    if (cfgSelCollByNch && (tracks.size() < cfgCutMultMin || tracks.size() >= cfgCutMultMax)) {
      return;
    }
    if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent < cfgCutCentMin || cent >= cfgCutCentMax)) {
      return;
    }

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    fillYield(collision, tracks);

    same->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);
    const auto& ft0 = collision.foundFT0();
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, ft0, collision.posZ(), SameEvent, cent, weightCent);
  }
  PROCESS_SWITCH(LongRangeDihadronCor, processSame, "Process same event", true);

  // the process for filling the mixed events
  void processMixed(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (cfgSelCollByNch && (tracks1.size() < cfgCutMultMin || tracks1.size() >= cfgCutMultMax))
        continue;

      if (cfgSelCollByNch && (tracks2.size() < cfgCutMultMin || tracks2.size() >= cfgCutMultMax))
        continue;

      float cent1 = -1;
      float cent2 = -1;
      if (!cfgCentTableUnavailable) {
        cent1 = getCentrality(collision1);
        cent2 = getCentrality(collision2);
      }
      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), cent1, false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), cent2, false))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent1 < cfgCutCentMin || cent1 >= cfgCutCentMax))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent2 < cfgCutCentMin || cent2 >= cfgCutCentMax))
        continue;

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0()))
        continue;

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;
      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }
      float weightCent = 1.0f;
      if (!cfgCentTableUnavailable)
        getCentralityWeight(weightCent, cent1);
      const auto& ft0 = collision2.foundFT0();
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, collision1.posZ(), MixedEvent, cent1, eventWeight * weightCent);
    }
  }

  PROCESS_SWITCH(LongRangeDihadronCor, processMixed, "Process mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LongRangeDihadronCor>(cfgc),
  };
}
