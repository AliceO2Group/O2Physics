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
// \author: sourav kundu
// \email: sourav.kundu@cern.ch

#include "PWGLF/DataModel/ZDC2StageCalibrationTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile3D.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

using BCsRun3 = o2::soa::Join<o2::aod::BCsWithTimestamps, o2::aod::Run3MatchedToBCSparse>;

struct zdc2stagecalibration {

  o2::framework::Produces<o2::aod::ZDC2StageCalibs> zdccaltable;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  } cfgCcdbParam;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 100.0f, "Centrality cut Max"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Centrality cut Min"};
  Configurable<int> calibrationStage{"calibrationStage", 2, "1: derive detector calibration QA, 2: apply detector calibration and produce corrected Q"};
  Configurable<float> cfgTimeSliceMinutes{"cfgTimeSliceMinutes", 30.0f, "Calibration time-slice width in minutes"};
  Configurable<float> cfgMaxRunHours{"cfgMaxRunHours", 48.0f, "Maximum run duration shown in QA"};
  Configurable<bool> followpub{"followpub", false, "flag to use alphaZDC"};
  Configurable<bool> usecfactor{"usecfactor", false, "use c factor"};
  Configurable<bool> useGainCallib{"useGainCallib", true, "use stage-1 tower gain calibration"};
  Configurable<std::string> confGainPath{"confGainPath", "Users/s/skundu/ZDC/Stage1Gain", "CCDB path to the run-wise stage-1 gain map"};
  Configurable<bool> derive1nTowerResponse{"derive1nTowerResponse", true, "save joint relative tower fractions for run-wise 1n response fit"};
  // Select one spatial correction: 0=gain only, 1=November24 1n Q-map, 2=hadronic cross-coordinate x centrality.
  // For backward compatibility, use1nTowerCorrection=true selects mode 1 if this is 0.
  Configurable<int> spatialCorrectionMode{"spatialCorrectionMode", 0, "0: no spatial, 1: November24 Q-space tower map, 2: centrality-dependent cross-space tower map"};
  Configurable<bool> deriveCentralitySpatial{"deriveCentralitySpatial", false, "write gain-only hadronic cross-space x centrality regression moments (train/validation)"};
  Configurable<std::string> confCentralitySpatialPath{"confCentralitySpatialPath", "", "CCDB TH3D (crossX,crossY,centClass*8+arm*4+tower)"};
  Configurable<std::string> confCentralitySpatialLocalFile{"confCentralitySpatialLocalFile", "", "Optional local ROOT file with CentralitySpatialTowerFactors for dry-run; overrides CCDB"};
  Configurable<int> cfgCentralitySpatialLocalRun{"cfgCentralitySpatialLocalRun", 0, "Expected run for local spatial map (0: allow any; single-run test only)"};
  Configurable<float> cfgCentralitySpatialMaxChange{"cfgCentralitySpatialMaxChange", 0.15f, "reject CCDB tower factors exceeding this absolute deviation from unity"};
  Configurable<bool> use1nTowerCorrection{"use1nTowerCorrection", false, "apply Q-dependent relative tower-energy factors derived from 1n"};
  Configurable<std::string> conf1nTowerPath{"conf1nTowerPath", "", "CCDB TH3D: gain Qx,Qy, A1-4/C1-4 relative response factors"};
  Configurable<bool> deriveSingleNeutronResponse{"deriveSingleNeutronResponse", true, "compact 1n/2n/3n neutron response QA"};
  Configurable<float> cfgSingleNeutronRawEnergyMax{"cfgSingleNeutronRawEnergyMax", 12.f, "maximum common energy for compact neutron QA"};
  // Generous, disjoint common-energy windows (native units). Each side is selected independently.
  // Defaults are approximate 1n, 2n, 3n regions; adjust from measured Ecommon peaks.
  Configurable<float> cfg1nMin{"cfg1nMin", 1.3f, "1n window lower Ecommon"};
  Configurable<float> cfg1nMax{"cfg1nMax", 4.3f, "1n window upper Ecommon"};
  Configurable<float> cfg2nMin{"cfg2nMin", 4.3f, "2n window lower Ecommon"};
  Configurable<float> cfg2nMax{"cfg2nMax", 7.0f, "2n window upper Ecommon"};
  Configurable<float> cfg3nMin{"cfg3nMin", 7.0f, "3n window lower Ecommon"};
  Configurable<float> cfg3nMax{"cfg3nMax", 10.5f, "3n window upper Ecommon"};
  Configurable<bool> cfgTPCV1Correlation{"cfgTPCV1Correlation", false, "enable TPC track loop and TPC-ZDC first-harmonic correlations versus eta and centrality (default off)"};
  Configurable<bool> deriveQRecentering{"deriveQRecentering", false, "store structured71 recentering moments after enabled tower-energy corrections"};
  Configurable<bool> useQRecentering{"useQRecentering", false, "apply Q-vector recentering calibration"};
  Configurable<int> qRecenteringMethod{"qRecenteringMethod", 4, "Only Structured71 recentering is supported (must be 4)"};
  Configurable<std::string> confQRecenteringStructured71Path{"confQRecenteringStructured71Path", "", "CCDB path to the 71-feature structured Q-vector recentering coefficients"};
  Configurable<float> cfgTPCResponsePtMin{"cfgTPCResponsePtMin", 0.5f, "TPC response QA minimum pT (strict)"};
  Configurable<float> cfgTPCResponsePtMax{"cfgTPCResponsePtMax", 3.0f, "TPC response QA maximum pT (strict)"};
  Configurable<float> cfgTPCResponseEtaMax{"cfgTPCResponseEtaMax", 0.8f, "TPC response QA maximum absolute eta (strict)"};
  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", true, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  RCTFlagsChecker rctChecker;
  int currentRunNumber = -999;
  int64_t bcSOR = -1;
  uint64_t sorTimestamp = 0;
  uint64_t eorTimestamp = 0;
  TH2D* gainprofile = nullptr;
  TH3D* qRecenteringProfile = nullptr;
  // Eight position-dependent relative tower factors, Z bins A1-A4,C1-C4.
  // Map lookup coordinates are the gain-only reconstructed Qx/Qy (cm).
  TH3D* oneNeutronTowerMap = nullptr;

  // Fixed centrality classes used both by the collection task and fitting macro.
  // 0:[0,5), 1:[5,10), 2:[10,20), 3:[20,30).
  static constexpr int kCentSpatialClasses = 4;
  static constexpr int kCentSpatialMomentCount = 15; // 10 zij products, 4 zi sums, 1 count
  int getSpatialMode() const
  {
    return (spatialCorrectionMode.value == 0 && use1nTowerCorrection.value) ? 1 : spatialCorrectionMode.value;
  }
  int getCentSpatialClass(double cent) const
  {
    if (cent < 0. || cent >= 30. || !std::isfinite(cent))
      return -1;
    return cent < 5. ? 0 : cent < 10. ? 1
                         : cent < 20. ? 2
                                      : 3;
  }
  TH3D* centralitySpatialTowerMap = nullptr;
  std::unique_ptr<TH3D> localCentralitySpatialMap;
  static constexpr int kQRecenteringNComponents = 4;
  static constexpr int kQRecenteringMultidimNFeatures = 35;
  // Structured regression diagnostic basis. The existing 35 cubic features are
  // kept unchanged and augmented with truncated-cubic 1D basis functions for
  // time, Vx, Vy and Vz. The additional basis is intentionally overcomplete so
  // that the final subset and regularisation can be chosen offline without
  // another derivation pass.
  static constexpr int kQRecenteringStructuredNKnots = 9;
  static constexpr int kQRecenteringStructuredNFeatures = kQRecenteringMultidimNFeatures + 4 * kQRecenteringStructuredNKnots;
  static constexpr int kQRecenteringStructuredNMatrixMoments = kQRecenteringStructuredNFeatures * (kQRecenteringStructuredNFeatures + 1) / 2;
  static constexpr int kQRecenteringStructuredNMoments = kQRecenteringStructuredNMatrixMoments + kQRecenteringNComponents * kQRecenteringStructuredNFeatures + kQRecenteringNComponents;
  static constexpr int kQRecenteringStructuredFlatnessNMoments = kQRecenteringStructuredNFeatures + 2 * kQRecenteringNComponents;

  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    if (spatialCorrectionMode.value < 0 || spatialCorrectionMode.value > 2) {
      LOGF(fatal, "spatialCorrectionMode must be 0 (gain only), 1 (November24), or 2 (centrality cross-space)");
    }
    if (spatialCorrectionMode.value == 2 && use1nTowerCorrection.value) {
      LOGF(fatal, "Do not enable old use1nTowerCorrection together with spatialCorrectionMode=2");
    }
    if (getSpatialMode() != 0 && calibrationStage.value != 2) {
      LOGF(fatal, "Spatial correction application requires calibrationStage=2");
    }
    if (getSpatialMode() == 2 && (!useGainCallib.value || (confCentralitySpatialPath.value.empty() && confCentralitySpatialLocalFile.value.empty()))) {
      LOGF(fatal, "Mode 2 requires useGainCallib=true and either a CCDB path or a local ROOT map file");
    }
    if (!(cfgCentralitySpatialMaxChange.value > 0.f && cfgCentralitySpatialMaxChange.value < 1.f)) {
      LOGF(fatal, "cfgCentralitySpatialMaxChange must be between 0 and 1");
    }
    if (deriveCentralitySpatial.value && (calibrationStage.value != 2 || !useGainCallib.value || getSpatialMode() != 0)) {
      LOGF(fatal, "deriveCentralitySpatial requires stage 2, valid stage-1 gain, and spatialCorrectionMode=0");
    }
    if (qRecenteringMethod.value != 4) {
      LOGF(fatal, "Only qRecenteringMethod=4 (structured71) is supported");
    }
    if (deriveQRecentering.value && useQRecentering.value) {
      LOGF(fatal, "Choose either deriveQRecentering or useQRecentering, not both");
    }
    if (cfgTPCV1Correlation.value &&
        !(cfgTPCResponsePtMax.value > cfgTPCResponsePtMin.value && cfgTPCResponseEtaMax.value > 0.f)) {
      LOGF(fatal, "Invalid TPC v1 QA track selection");
    }
    if ((deriveSingleNeutronResponse.value || derive1nTowerResponse.value) && (!(cfgSingleNeutronRawEnergyMax.value > 0.f) || cfgSingleNeutronRawEnergyMax.value > 12.f)) {
      LOGF(fatal, "Single-neutron common-energy cutoff must be >0 and <=12 (native units)");
    }
    if ((deriveSingleNeutronResponse.value || derive1nTowerResponse.value) &&
        !(0.f < cfg1nMin.value && cfg1nMin.value < cfg1nMax.value &&
          cfg1nMax.value <= cfg2nMin.value && cfg2nMin.value < cfg2nMax.value &&
          cfg2nMax.value <= cfg3nMin.value && cfg3nMin.value < cfg3nMax.value &&
          cfg3nMax.value <= cfgSingleNeutronRawEnergyMax.value)) {
      LOGF(fatal, "Invalid/overlapping neutron energy windows; require ordered 1n,2n,3n windows within Emax");
    }
    if (derive1nTowerResponse.value && !useGainCallib.value) {
      LOGF(fatal, "1n response derivation requires useGainCallib=true");
    }
    if (derive1nTowerResponse.value && getSpatialMode() == 1) {
      LOGF(fatal, "Fit-input collection requires gain-only towers: disable use1nTowerCorrection");
    }
    if (cfgTPCV1Correlation.value && calibrationStage.value != 2) {
      LOGF(fatal, "TPC v1(eta) numerator QA requires stage 2");
    }
    const int nTimeBins = static_cast<int>(std::ceil(cfgMaxRunHours.value * 60.f / cfgTimeSliceMinutes.value));
    AxisSpec timeAxis = {nTimeBins, 0.0, cfgMaxRunHours.value, "time from SOR (h)"};
    AxisSpec towerAxis = {4, 0.0, 4.0, "tower"};
    AxisSpec crossAxis = {120, -1.0, 1.0, "cross asymmetry"};
    AxisSpec momentAxis = {15, 0.0, 15.0, "regression moment"};
    AxisSpec phiAxis = {72, -3.14159265358979323846, 3.14159265358979323846, "#phi"};
    AxisSpec crossCoordAxis = {2, 0.0, 2.0, "coordinate"};
    AxisSpec ratioAxis = {240, 0.0, 2.4, "#Sigma tower/common"};
    AxisSpec qRecenteringStructuredMomentAxis = {kQRecenteringStructuredNMoments, 0.0, static_cast<double>(kQRecenteringStructuredNMoments), "structured regression moment"};
    AxisSpec qRecenteringStructuredFlatnessMomentAxis = {kQRecenteringStructuredFlatnessNMoments, 0.0, static_cast<double>(kQRecenteringStructuredFlatnessNMoments), "structured flatness sufficient-statistic index"};
    AxisSpec centralityAxis = {100, 0.0, 100.0, "centrality (%)"};
    if (deriveSingleNeutronResponse.value) {
      const AxisSpec energyAxis{120, 0., static_cast<double>(cfgSingleNeutronRawEnergyMax.value), "Ecommon (native units)"};
      const AxisSpec fractionAxis{80, 0., 2., "gain tower / common"};
      const AxisSpec qAxis{100, -2.5, 2.5, "Q (cm)"};
      const AxisSpec snMomentAxis{14, 0., 14., "Q moment index"};
      const AxisSpec statusAxis{4, 0., 4., "0=pre tower cut;1=both valid;2=ZNA valid;3=ZNC valid"};
      const AxisSpec crossQA{30, -1., 1., "gain-stage cross coordinate"};
      histos.add("NeutronQA/hSelectionStatus", "Energy-class selection before and after tower cut;side/class (A1,A2,A3,C1,C2,C3);status", kTH2F, {{6, 0., 6.}, statusAxis});
      histos.add("NeutronQA/hEnergyVsCentA", "ZNA energy spectrum;centrality (%);Ecommon", kTH2F, {centralityAxis, energyAxis});
      histos.add("NeutronQA/hEnergyVsCentC", "ZNC energy spectrum;centrality (%);Ecommon", kTH2F, {centralityAxis, energyAxis});
      // Only 1n for the response constraints; total around 12 small 30x30 maps.
      histos.add("NeutronQA/hCross1nA", "ZNA 1n gain-stage cross;crossX;crossY", kTH2F, {crossQA, crossQA});
      histos.add("NeutronQA/hCross1nC", "ZNC 1n gain-stage cross;crossX;crossY", kTH2F, {crossQA, crossQA});
      histos.add("NeutronQA/pRecenteredQxVsCross1nA", "ZNA 1n Qx after recenter vs gain cross;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pRecenteredQyVsCross1nA", "ZNA 1n Qy after recenter vs gain cross;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pRecenteredQxVsCross1nC", "ZNC 1n Qx after recenter vs gain cross;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pRecenteredQyVsCross1nC", "ZNC 1n Qy after recenter vs gain cross;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pGainTowerFractionVsCross1nA1", "A 1n gain tower1/common;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pGainTowerFractionVsCross1nA2", "A 1n gain tower2/common;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pGainTowerFractionVsCross1nA3", "A 1n gain tower3/common;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pGainTowerFractionVsCross1nA4", "A 1n gain tower4/common;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pGainTowerFractionVsCross1nC1", "C 1n gain tower1/common;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pGainTowerFractionVsCross1nC2", "C 1n gain tower2/common;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pGainTowerFractionVsCross1nC3", "C 1n gain tower3/common;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/pGainTowerFractionVsCross1nC4", "C 1n gain tower4/common;crossX;crossY", kTProfile2D, {crossQA, crossQA});
      histos.add("NeutronQA/hTowerFractionA1", "A tower1/common;Ecommon;fraction", kTH2F, {energyAxis, fractionAxis});
      histos.add("NeutronQA/hTowerFractionA2", "A tower2/common;Ecommon;fraction", kTH2F, {energyAxis, fractionAxis});
      histos.add("NeutronQA/hTowerFractionA3", "A tower3/common;Ecommon;fraction", kTH2F, {energyAxis, fractionAxis});
      histos.add("NeutronQA/hTowerFractionA4", "A tower4/common;Ecommon;fraction", kTH2F, {energyAxis, fractionAxis});
      histos.add("NeutronQA/hTowerFractionC1", "C tower1/common;Ecommon;fraction", kTH2F, {energyAxis, fractionAxis});
      histos.add("NeutronQA/hTowerFractionC2", "C tower2/common;Ecommon;fraction", kTH2F, {energyAxis, fractionAxis});
      histos.add("NeutronQA/hTowerFractionC3", "C tower3/common;Ecommon;fraction", kTH2F, {energyAxis, fractionAxis});
      histos.add("NeutronQA/hTowerFractionC4", "C tower4/common;Ecommon;fraction", kTH2F, {energyAxis, fractionAxis});
      histos.add("NeutronQA/hQxyRecentered1nA", "A 1n Q after Recentered;Qx (cm);Qy (cm)", kTH2F, {qAxis, qAxis});
      histos.add("NeutronQA/pMomentsRecentered1nA", "A 1n Recentered Q moments;centrality (%);moment index;mean", kTProfile2D, {centralityAxis, snMomentAxis});
      histos.add("NeutronQA/hQxyRecentered2nA", "A 2n Q after Recentered;Qx (cm);Qy (cm)", kTH2F, {qAxis, qAxis});
      histos.add("NeutronQA/pMomentsRecentered2nA", "A 2n Recentered Q moments;centrality (%);moment index;mean", kTProfile2D, {centralityAxis, snMomentAxis});
      histos.add("NeutronQA/hQxyRecentered3nA", "A 3n Q after Recentered;Qx (cm);Qy (cm)", kTH2F, {qAxis, qAxis});
      histos.add("NeutronQA/pMomentsRecentered3nA", "A 3n Recentered Q moments;centrality (%);moment index;mean", kTProfile2D, {centralityAxis, snMomentAxis});
      histos.add("NeutronQA/hQxyRecentered1nC", "C 1n Q after Recentered;Qx (cm);Qy (cm)", kTH2F, {qAxis, qAxis});
      histos.add("NeutronQA/pMomentsRecentered1nC", "C 1n Recentered Q moments;centrality (%);moment index;mean", kTProfile2D, {centralityAxis, snMomentAxis});
      histos.add("NeutronQA/hQxyRecentered2nC", "C 2n Q after Recentered;Qx (cm);Qy (cm)", kTH2F, {qAxis, qAxis});
      histos.add("NeutronQA/pMomentsRecentered2nC", "C 2n Recentered Q moments;centrality (%);moment index;mean", kTProfile2D, {centralityAxis, snMomentAxis});
      histos.add("NeutronQA/hQxyRecentered3nC", "C 3n Q after Recentered;Qx (cm);Qy (cm)", kTH2F, {qAxis, qAxis});
      histos.add("NeutronQA/pMomentsRecentered3nC", "C 3n Recentered Q moments;centrality (%);moment index;mean", kTProfile2D, {centralityAxis, snMomentAxis});
    }
    if (cfgTPCV1Correlation.value) {
      // Track-by-track numerators, 16 eta bins and 40 centrality bins; no THnSparse.
      // indices 0:Q1xTPC*QxA, 1:Q1yTPC*QyA, 2:Q1xTPC*QxC, 3:Q1yTPC*QyC.
      histos.add("TPCV1EtaQA/pNumerators", "TPC-ZDC first harmonic;centrality (%);track #eta;component", kTProfile3D,
                 {{100, 0., 100.}, {16, -0.8, 0.8}, {4, 0., 4.}});
      histos.add("TPCV1EtaQA/hTrackCounts", "Selected tracks;centrality (%);track #eta", kTH2F,
                 {{100, 0., 100.}, {16, -0.8, 0.8}});
    }

    AxisSpec resolutionCentralityAxis = {100, 0.0, 100.0, "centrality (%)"};
    AxisSpec resolutionVzAxis = {20, -10.0, 10.0, "v_{z} (cm)"};
    AxisSpec vertexXYAxis = {100, -0.1, 0.1, "vertex x/y (cm)"};
    AxisSpec vertexZAxis = {100, -10.0, 10.0, "vertex z (cm)"};
    AxisSpec qComponentAxis = {4, 0.0, 4.0, "Q component"};
    AxisSpec qFinalMomentAxis = {10, 0.0, 10.0, "Q moments (A: 0-4, C: 5-9)"};
    AxisSpec qCorrelationComponentAxis = {4, 0.0, 4.0, "A-C correlation component"};
    AxisSpec phiHarmonicAxis = {8, 0.0, 8.0, "cos1,sin1,...,cos4,sin4"};
    AxisSpec qValueAxis = {200, -4.0, 4.0, "Q"};

    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{10, 0.0, 10.0}});
    auto hEvtSelInfo = histos.get<TH1>(HIST("hEvtSelInfo"));
    hEvtSelInfo->GetXaxis()->SetBinLabel(1, "All");
    hEvtSelInfo->GetXaxis()->SetBinLabel(2, "has ZDC");
    hEvtSelInfo->GetXaxis()->SetBinLabel(3, "valid common");
    hEvtSelInfo->GetXaxis()->SetBinLabel(4, "valid towers");
    hEvtSelInfo->GetXaxis()->SetBinLabel(5, "RCT CBT+ZDC");
    hEvtSelInfo->GetXaxis()->SetBinLabel(6, "event sel.");
    hEvtSelInfo->GetXaxis()->SetBinLabel(7, "gain valid");
    hEvtSelInfo->GetXaxis()->SetBinLabel(8, "valid Q");
    hEvtSelInfo->GetXaxis()->SetBinLabel(9, "accepted");
    if (deriveCentralitySpatial.value || getSpatialMode() == 2) {
      const AxisSpec centCrossAxis{40, -1., 1., "gain-only cross coordinate"};
      const AxisSpec centClassAxis{kCentSpatialClasses, 0., static_cast<double>(kCentSpatialClasses), "centrality class (0:0-5,1:5-10,2:10-20,3:20-30)"};
      const AxisSpec spatialQAxis{80, -2.5, 2.5, "Q (cm), linear tower weighting"};
      const AxisSpec deltaQAxis{100, 0., 50., "centrality (%)"};
      // TProfile3D visualises the mean energy response, not a forced-flat occupancy.
      for (const char* side : {"ZA", "ZC"}) {
        for (int i = 0; i < 5; ++i) {
          const std::string title = i < 4 ? "gain tower/common" : "gain sum/common";
          const std::string name = std::string("CentSpatial/pGainResponse") + std::to_string(i) + side;
          histos.add(name.c_str(), (title + ";Xcross;Ycross;centrality class;response").c_str(), kTProfile3D, {centCrossAxis, centCrossAxis, centClassAxis});
        }
        const std::string after = std::string("CentSpatial/pCorrectedTotal") + side;
        histos.add(after.c_str(), "corrected tower sum/common in GAIN cross coordinates;Xcross;Ycross;centrality class;response", kTProfile3D, {centCrossAxis, centCrossAxis, centClassAxis});
      }
      histos.add("CentSpatial/hQGainZA", "Gain-only ZNA Q;Qx;Qy;centrality class", kTH3D, {spatialQAxis, spatialQAxis, centClassAxis});
      histos.add("CentSpatial/hQGainZC", "Gain-only ZNC Q;Qx;Qy;centrality class", kTH3D, {spatialQAxis, spatialQAxis, centClassAxis});
      histos.add("CentSpatial/hQSpatialZA", "After selected spatial correction ZNA Q;Qx;Qy;centrality class", kTH3D, {spatialQAxis, spatialQAxis, centClassAxis});
      histos.add("CentSpatial/hQSpatialZC", "After selected spatial correction ZNC Q;Qx;Qy;centrality class", kTH3D, {spatialQAxis, spatialQAxis, centClassAxis});
      histos.add("CentSpatial/pDeltaQxZA", "ZNA <QxSpatial-QxGain>;centrality (%);<DeltaQx>", kTProfile, {deltaQAxis});
      histos.add("CentSpatial/pDeltaQyZA", "ZNA <QySpatial-QyGain>;centrality (%);<DeltaQy>", kTProfile, {deltaQAxis});
      histos.add("CentSpatial/pDeltaQMagZA", "ZNA <|QSpatial-QGain|>;centrality (%);<|DeltaQ|>", kTProfile, {deltaQAxis});
      histos.add("CentSpatial/pDeltaQxZC", "ZNC <QxSpatial-QxGain>;centrality (%);<DeltaQx>", kTProfile, {deltaQAxis});
      histos.add("CentSpatial/pDeltaQyZC", "ZNC <QySpatial-QyGain>;centrality (%);<DeltaQy>", kTProfile, {deltaQAxis});
      histos.add("CentSpatial/pDeltaQMagZC", "ZNC <|QSpatial-QGain|>;centrality (%);<|DeltaQ|>", kTProfile, {deltaQAxis});
      histos.add("CentSpatial/hCorrectionStatus", "Mode2 status per side: 0=no map/class,1=ZNA only,2=ZNC only,3=both;status;events", kTH1F, {{4, 0., 4.}});
      if (deriveCentralitySpatial.value) {
        // Per (Xcross,Ycross,centrality-class) cell, retain the EXACT sums
        // Σzi*zj (i<=j), Σzi and N, zi=E_i/Ecommon. Independent train and
        // validation splits. This is sufficient for a regularized global fit.
        // 15 histograms x 2 arms x 2 splits; no 1n selection.
        for (const char* side : {"ZA", "ZC"}) {
          for (const char* split : {"Train", "Valid"}) {
            for (int m = 0; m < kCentSpatialMomentCount; ++m) {
              const std::string name = std::string("CentSpatial/hMom") + split + side + "_" + (m < 10 ? "0" : "") + std::to_string(m);
              histos.add(name.c_str(), "Sum of normalized energy regression moment;gain Xcross;gain Ycross;centrality class", kTH3D, {centCrossAxis, centCrossAxis, centClassAxis});
            }
          }
        }
      }
    }

    histos.add("GainQA/hSumOverCommonVsTimeZNA", "ZNA tower sum/common vs time;time from SOR (h);#Sigma T_{i}/C", kTH2F, {timeAxis, ratioAxis});
    histos.add("GainQA/hSumOverCommonVsTimeZNC", "ZNC tower sum/common vs time;time from SOR (h);#Sigma T_{i}/C", kTH2F, {timeAxis, ratioAxis});
    histos.add("GainQA/pTowerOverCommonVsTimeZNA", "ZNA individual tower/common vs time;time from SOR (h);tower;<T_{i}/C>", kTProfile2D, {timeAxis, towerAxis});
    histos.add("GainQA/pTowerOverCommonVsTimeZNC", "ZNC individual tower/common vs time;time from SOR (h);tower;<T_{i}/C>", kTProfile2D, {timeAxis, towerAxis});
    histos.add("TowerResponseQA/pCrossVsTimeZNA", "ZNA cross asymmetry vs time;time from SOR (h);0:X cross, 1:Y cross;<cross>", kTProfile2D, {timeAxis, crossCoordAxis});
    histos.add("TowerResponseQA/pCrossVsTimeZNC", "ZNC cross asymmetry vs time;time from SOR (h);0:X cross, 1:Y cross;<cross>", kTProfile2D, {timeAxis, crossCoordAxis});
    histos.add("TowerResponseQA/pResponseVsCrossXYZNA", "ZNA response map;X cross;Y cross;<#Sigma T_{i}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pResponseVsCrossXYZNC", "ZNC response map;X cross;Y cross;<#Sigma T_{i}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pCentroidVsTimeZNA", "ZNA centroid vs time;time from SOR (h);0:x, 1:y;<centroid (cm)>", kTProfile2D, {timeAxis, crossCoordAxis});
    histos.add("TowerResponseQA/pCentroidVsTimeZNC", "ZNC centroid vs time;time from SOR (h);0:x, 1:y;<centroid (cm)>", kTProfile2D, {timeAxis, crossCoordAxis});
    histos.add("TowerResponseQA/hCrossXYZNA", "ZNA cross occupancy;X cross;Y cross", kTH2F, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/hCrossXYZNC", "ZNC cross occupancy;X cross;Y cross", kTH2F, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pTower1OverCommonVsCrossXYZNA", "ZNA tower 1 response;X cross;Y cross;<T_{1}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pTower2OverCommonVsCrossXYZNA", "ZNA tower 2 response;X cross;Y cross;<T_{2}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pTower3OverCommonVsCrossXYZNA", "ZNA tower 3 response;X cross;Y cross;<T_{3}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pTower4OverCommonVsCrossXYZNA", "ZNA tower 4 response;X cross;Y cross;<T_{4}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pTower1OverCommonVsCrossXYZNC", "ZNC tower 1 response;X cross;Y cross;<T_{1}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pTower2OverCommonVsCrossXYZNC", "ZNC tower 2 response;X cross;Y cross;<T_{2}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pTower3OverCommonVsCrossXYZNC", "ZNC tower 3 response;X cross;Y cross;<T_{3}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("TowerResponseQA/pTower4OverCommonVsCrossXYZNC", "ZNC tower 4 response;X cross;Y cross;<T_{4}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("PhiQA/hPhiRawZNA", "ZNA #phi raw;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiRawZNC", "ZNC #phi raw;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterGainZNA", "ZNA #phi after gain calibration;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterGainZNC", "ZNC #phi after gain calibration;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterTowerCorrectionZNA", "ZNA #phi after 1n tower-energy correction (when enabled);#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterTowerCorrectionZNC", "ZNC #phi after 1n tower-energy correction (when enabled);#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiRawVsCentralityZNA", "ZNA #phi raw vs centrality;centrality (%);#phi", kTH2F, {centralityAxis, phiAxis});
    histos.add("PhiQA/hPhiRawVsCentralityZNC", "ZNC #phi raw vs centrality;centrality (%);#phi", kTH2F, {centralityAxis, phiAxis});
    histos.add("PhiQA/hPhiAfterGainVsCentralityZNA", "ZNA #phi after gain vs centrality;centrality (%);#phi", kTH2F, {centralityAxis, phiAxis});
    histos.add("PhiQA/hPhiAfterGainVsCentralityZNC", "ZNC #phi after gain vs centrality;centrality (%);#phi", kTH2F, {centralityAxis, phiAxis});
    histos.add("PhiQA/hPhiAfterTowerCorrectionVsCentralityZNA", "ZNA #phi after 1n tower-energy correction (when enabled) vs centrality;centrality (%);#phi", kTH2F, {centralityAxis, phiAxis});
    histos.add("PhiQA/hPhiAfterTowerCorrectionVsCentralityZNC", "ZNC #phi after 1n tower-energy correction (when enabled) vs centrality;centrality (%);#phi", kTH2F, {centralityAxis, phiAxis});

    // Stage-1 correction input. For each time bin the first 10 y bins contain sum(T_i*T_j), the next 4 contain sum(C*T_i), and the last contains the event count.
    histos.add("GainCalibration/hGainMomentsZNA", "ZNA linear-regression moments;time from SOR (h);moment index", kTH2D, {timeAxis, momentAxis});
    histos.add("GainCalibration/hGainMomentsZNC", "ZNC linear-regression moments;time from SOR (h);moment index", kTH2D, {timeAxis, momentAxis});
    if (deriveQRecentering.value) {
      histos.add("QRecenteringCalibration/hRegressionMomentsStructuredTrain", "Structured Q recentering moments (train);centrality (%);moment index", kTH2D, {centralityAxis, qRecenteringStructuredMomentAxis});
      histos.add("QRecenteringCalibration/hRegressionMomentsStructuredValidation", "Structured Q recentering moments (validation);centrality (%);moment index", kTH2D, {centralityAxis, qRecenteringStructuredMomentAxis});
      histos.add("QRecenteringCalibration/hStructuredFlatnessMomentsTimeTrain", "Structured flatness statistics vs time (train);centrality (%);time from SOR (h);statistic index", kTH3D, {centralityAxis, timeAxis, qRecenteringStructuredFlatnessMomentAxis});
      histos.add("QRecenteringCalibration/hStructuredFlatnessMomentsVxTrain", "Structured flatness statistics vs v_{x} (train);centrality (%);v_{x} (cm);statistic index", kTH3D, {centralityAxis, vertexXYAxis, qRecenteringStructuredFlatnessMomentAxis});
      histos.add("QRecenteringCalibration/hStructuredFlatnessMomentsVyTrain", "Structured flatness statistics vs v_{y} (train);centrality (%);v_{y} (cm);statistic index", kTH3D, {centralityAxis, vertexXYAxis, qRecenteringStructuredFlatnessMomentAxis});
      histos.add("QRecenteringCalibration/hStructuredFlatnessMomentsVzTrain", "Structured flatness statistics vs v_{z} (train);centrality (%);v_{z} (cm);statistic index", kTH3D, {centralityAxis, vertexZAxis, qRecenteringStructuredFlatnessMomentAxis});
      histos.add("QRecenteringCalibration/hStructuredFlatnessMomentsTimeValidation", "Structured flatness statistics vs time (validation);centrality (%);time from SOR (h);statistic index", kTH3D, {centralityAxis, timeAxis, qRecenteringStructuredFlatnessMomentAxis});
      histos.add("QRecenteringCalibration/hStructuredFlatnessMomentsVxValidation", "Structured flatness statistics vs v_{x} (validation);centrality (%);v_{x} (cm);statistic index", kTH3D, {centralityAxis, vertexXYAxis, qRecenteringStructuredFlatnessMomentAxis});
      histos.add("QRecenteringCalibration/hStructuredFlatnessMomentsVyValidation", "Structured flatness statistics vs v_{y} (validation);centrality (%);v_{y} (cm);statistic index", kTH3D, {centralityAxis, vertexXYAxis, qRecenteringStructuredFlatnessMomentAxis});
      histos.add("QRecenteringCalibration/hStructuredFlatnessMomentsVzValidation", "Structured flatness statistics vs v_{z} (validation);centrality (%);v_{z} (cm);statistic index", kTH3D, {centralityAxis, vertexZAxis, qRecenteringStructuredFlatnessMomentAxis});
    }

    histos.add("QRecenteringQA/pQBeforeVsCentrality", "Q before recentering vs centrality;centrality (%);Q component;<Q>", kTProfile2D, {centralityAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQBeforeVsTime", "Q before recentering vs time;time from SOR (h);Q component;<Q>", kTProfile2D, {timeAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQBeforeVsVx", "Q before recentering vs v_{x};v_{x} (cm);Q component;<Q>", kTProfile2D, {vertexXYAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQBeforeVsVy", "Q before recentering vs v_{y};v_{y} (cm);Q component;<Q>", kTProfile2D, {vertexXYAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQBeforeVsVz", "Q before recentering vs v_{z};v_{z} (cm);Q component;<Q>", kTProfile2D, {vertexZAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQAfterVsCentrality", "Q after recentering vs centrality;centrality (%);Q component;<Q>", kTProfile2D, {centralityAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQAfterVsTime", "Q after recentering vs time;time from SOR (h);Q component;<Q>", kTProfile2D, {timeAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQAfterVsVx", "Q after recentering vs v_{x};v_{x} (cm);Q component;<Q>", kTProfile2D, {vertexXYAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQAfterVsVy", "Q after recentering vs v_{y};v_{y} (cm);Q component;<Q>", kTProfile2D, {vertexXYAxis, qComponentAxis});
    histos.add("QRecenteringQA/pQAfterVsVz", "Q after recentering vs v_{z};v_{z} (cm);Q component;<Q>", kTProfile2D, {vertexZAxis, qComponentAxis});
    histos.add("ResolutionQA/pCosPsiAPsiCVsCentrality", "ZDC A-C spectator-plane correlation;centrality (%);<#cos(#Psi_{A}-#Psi_{C})>", kTProfile, {resolutionCentralityAxis});
    histos.add("ResolutionQA/pCosPsiAPsiCVsCentralityVz", "ZDC A-C correlation;centrality (%);v_{z} (cm); <#cos(#Psi_{A}-#Psi_{C})>", kTProfile2D, {resolutionCentralityAxis, resolutionVzAxis});
    // Monitor Q azimuth after structured71 recentering.

    // A-C Q-vector correlation QA. These histograms are diagnostic only; no Qx/Qy mixing is applied.
    histos.add("PhiQA/hPhiAfterRecenteringZNA", "ZNA #phi after Q recentering;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterRecenteringZNC", "ZNC #phi after Q recentering;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterRecenteringVsCentralityZNA", "ZNA #phi after Q recentering vs centrality;centrality (%);#phi", kTH2F, {centralityAxis, phiAxis});
    histos.add("PhiQA/hPhiAfterRecenteringVsCentralityZNC", "ZNC #phi after Q recentering vs centrality;centrality (%);#phi", kTH2F, {centralityAxis, phiAxis});

    // A-C Q correlations versus centrality. Component convention:
    // 0: QxA*QxC, 1: QyA*QyC, 2: QxA*QyC, 3: QyA*QxC.
    // These are QA only; they are never used in the recentering correction.
    histos.add("QCorrelationQA/pQACorrelationBeforeVsCentrality", "A-C Q correlations before recentering;centrality (%);correlation component;<Q_{A}Q_{C}>", kTProfile2D, {centralityAxis, qCorrelationComponentAxis});
    histos.add("QCorrelationQA/pQACorrelationAfterVsCentrality", "A-C Q correlations after recentering;centrality (%);correlation component;<Q_{A}Q_{C}>", kTProfile2D, {centralityAxis, qCorrelationComponentAxis});

    // Final-calibration QA: after structured71 recentering (or un-recentered in derivation pass).
    histos.add("final_calibration/hQxAVsCentrality", "Final QxA vs centrality;centrality (%);Q_{x}^{A}", kTH2F, {centralityAxis, qValueAxis});
    histos.add("final_calibration/hQyAVsCentrality", "Final QyA vs centrality;centrality (%);Q_{y}^{A}", kTH2F, {centralityAxis, qValueAxis});
    histos.add("final_calibration/hQxCVsCentrality", "Final QxC vs centrality;centrality (%);Q_{x}^{C}", kTH2F, {centralityAxis, qValueAxis});
    histos.add("final_calibration/hQyCVsCentrality", "Final QyC vs centrality;centrality (%);Q_{y}^{C}", kTH2F, {centralityAxis, qValueAxis});
    histos.add("final_calibration/pQxAVsCentrality", "Final <QxA> vs centrality;centrality (%);<Q_{x}^{A}>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pQyAVsCentrality", "Final <QyA> vs centrality;centrality (%);<Q_{y}^{A}>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pQxCVsCentrality", "Final <QxC> vs centrality;centrality (%);<Q_{x}^{C}>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pQyCVsCentrality", "Final <QyC> vs centrality;centrality (%);<Q_{y}^{C}>", kTProfile, {centralityAxis});
    histos.add("final_calibration/hPhiAVsCentrality", "Final ZNA #phi vs centrality;centrality (%);#phi_{A}", kTH2F, {centralityAxis, phiAxis});
    histos.add("final_calibration/hPhiCVsCentrality", "Final ZNC #phi vs centrality;centrality (%);#phi_{C}", kTH2F, {centralityAxis, phiAxis});
    histos.add("final_calibration/pCosPsiAPsiCVsCentrality", "Final A-C plane correlation;centrality (%);<cos(#Psi_{A}-#Psi_{C})>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pCosPsiAPsiCVsCentralityVz", "Final A-C plane correlation vs centrality and v_{z};centrality (%);v_{z} (cm);<cos(#Psi_{A}-#Psi_{C})>", kTProfile2D, {centralityAxis, resolutionVzAxis});
    histos.add("final_calibration/pQxAQxCVsCentrality", "Final <QxA QxC>;centrality (%);<Q_{x}^{A}Q_{x}^{C}>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pQyAQyCVsCentrality", "Final <QyA QyC>;centrality (%);<Q_{y}^{A}Q_{y}^{C}>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pQxAQyCVsCentrality", "Final <QxA QyC>;centrality (%);<Q_{x}^{A}Q_{y}^{C}>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pQyAQxCVsCentrality", "Final <QyA QxC>;centrality (%);<Q_{y}^{A}Q_{x}^{C}>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pQxAVsCentralityVz", "Final <QxA> vs centrality and v_{z};centrality (%);v_{z} (cm);<Q_{x}^{A}>", kTProfile2D, {centralityAxis, resolutionVzAxis});
    histos.add("final_calibration/pQyAVsCentralityVz", "Final <QyA> vs centrality and v_{z};centrality (%);v_{z} (cm);<Q_{y}^{A}>", kTProfile2D, {centralityAxis, resolutionVzAxis});
    histos.add("final_calibration/pQxCVsCentralityVz", "Final <QxC> vs centrality and v_{z};centrality (%);v_{z} (cm);<Q_{x}^{C}>", kTProfile2D, {centralityAxis, resolutionVzAxis});
    histos.add("final_calibration/pQyCVsCentralityVz", "Final <QyC> vs centrality and v_{z};centrality (%);v_{z} (cm);<Q_{y}^{C}>", kTProfile2D, {centralityAxis, resolutionVzAxis});
    histos.add("final_calibration/pQMomentsVsCentrality", "Final Q moments: Qx,Qy,Qx^{2},Qy^{2},QxQy for A then C;centrality (%);moment index;<value>", kTProfile2D, {centralityAxis, qFinalMomentAxis});
    histos.add("final_calibration/pPhiHarmonicsAVsCentrality", "Final ZNA #phi harmonics: cos(n#phi),sin(n#phi) for n=1..4;centrality (%);harmonic index;<value>", kTProfile2D, {centralityAxis, phiHarmonicAxis});
    histos.add("final_calibration/pPhiHarmonicsCVsCentrality", "Final ZNC #phi harmonics: cos(n#phi),sin(n#phi) for n=1..4;centrality (%);harmonic index;<value>", kTProfile2D, {centralityAxis, phiHarmonicAxis});
    histos.add("final_calibration/pQMagAVsCentrality", "Final |Q_{A}| vs centrality;centrality (%);<|Q_{A}|>", kTProfile, {centralityAxis});
    histos.add("final_calibration/pQMagCVsCentrality", "Final |Q_{C}| vs centrality;centrality (%);<|Q_{C}|>", kTProfile, {centralityAxis});

    if (derive1nTowerResponse.value) {
      // Three independent normalized tower fractions f1,f2,f3; f4=1-f1-f2-f3.
      // Retains the joint four-tower energy-sharing information needed by the fit.
      // Run-wise 1n train / held-out validation, plus 2n and 3n transfer controls.
      const AxisSpec frac3Axis{64, 0., 1., "E_{i}/#Sigma E_{j}"};
      histos.add("OneNeutronTowerCalibration/hRelative1nA_Train", "ZNA 1n gain-only tower fractions (train);f1;f2;f3", kTH3F, {frac3Axis, frac3Axis, frac3Axis});
      histos.add("OneNeutronTowerCalibration/hRelative1nC_Train", "ZNC 1n gain-only tower fractions (train);f1;f2;f3", kTH3F, {frac3Axis, frac3Axis, frac3Axis});
      histos.add("OneNeutronTowerCalibration/hRelative1nA_Validation", "ZNA 1n gain-only tower fractions (validation);f1;f2;f3", kTH3F, {frac3Axis, frac3Axis, frac3Axis});
      histos.add("OneNeutronTowerCalibration/hRelative1nC_Validation", "ZNC 1n gain-only tower fractions (validation);f1;f2;f3", kTH3F, {frac3Axis, frac3Axis, frac3Axis});
      histos.add("OneNeutronTowerCalibration/hRelative2nA", "ZNA 2n gain-only fractions;f1;f2;f3", kTH3F, {frac3Axis, frac3Axis, frac3Axis});
      histos.add("OneNeutronTowerCalibration/hRelative2nC", "ZNC 2n gain-only fractions;f1;f2;f3", kTH3F, {frac3Axis, frac3Axis, frac3Axis});
      histos.add("OneNeutronTowerCalibration/hRelative3nA", "ZNA 3n gain-only fractions;f1;f2;f3", kTH3F, {frac3Axis, frac3Axis, frac3Axis});
      histos.add("OneNeutronTowerCalibration/hRelative3nC", "ZNC 3n gain-only fractions;f1;f2;f3", kTH3F, {frac3Axis, frac3Axis, frac3Axis});
    }
    if (deriveSingleNeutronResponse.value || derive1nTowerResponse.value || getSpatialMode() == 1) {
      histos.add("OneNeutronTowerQA/hCorrectionStatus", "Tower-map application status (0=disabled/missing, 1=A only, 2=C only, 3=both);status;events", kTH1F, {{4, 0., 4.}});
      const AxisSpec qAxisTower{100, -2.5, 2.5, "Q (cm)"};
      // 1n and 2n/3n before/after the tower correction, independent of recentering.
      histos.add("OneNeutronTowerQA/hBefore1nA", "Before Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hBefore1nC", "Before Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hBefore2nA", "Before Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hBefore2nC", "Before Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hBefore3nA", "Before Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hBefore3nC", "Before Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hAfter1nA", "After Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hAfter1nC", "After Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hAfter2nA", "After Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hAfter2nC", "After Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hAfter3nA", "After Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
      histos.add("OneNeutronTowerQA/hAfter3nC", "After Q before recentering;Qx (cm);Qy (cm)", kTH2F, {qAxisTower, qAxisTower});
    }

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  void updateRun(int runNumber, uint64_t timestamp)
  {
    if (runNumber == currentRunNumber) {
      return;
    }

    currentRunNumber = runNumber;
    auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), runNumber);
    sorTimestamp = runInfo.sor;
    eorTimestamp = runInfo.eor;
    bcSOR = runInfo.orbitSOR * static_cast<int64_t>(o2::constants::lhc::LHCMaxBunches);
    gainprofile = nullptr;
    qRecenteringProfile = nullptr;
    oneNeutronTowerMap = nullptr;
    centralitySpatialTowerMap = nullptr;
    localCentralitySpatialMap.reset();
    if (calibrationStage.value == 2 && useGainCallib.value && !confGainPath.value.empty()) {
      gainprofile = ccdb->getForTimeStamp<TH2D>(confGainPath.value, timestamp);
      if (!gainprofile) {
        LOGF(warn, "No ZDC stage-1 gain map found for run %d at timestamp %llu", runNumber, timestamp);
      }
    }
    if (calibrationStage.value == 2 && getSpatialMode() == 1 && conf1nTowerPath.value.empty()) {
      LOGF(warn, "1n relative tower correction is enabled but CCDB path is empty for run %d", runNumber);
    }
    if (calibrationStage.value == 2 && getSpatialMode() == 1 && !conf1nTowerPath.value.empty()) {
      oneNeutronTowerMap = ccdb->getForTimeStamp<TH3D>(conf1nTowerPath.value, timestamp);
      if (!oneNeutronTowerMap) {
        LOGF(warn, "Missing 1n relative-tower map for run %d; skipping only this step", runNumber);
      } else if (oneNeutronTowerMap->GetNbinsZ() != 8 ||
                 oneNeutronTowerMap->GetXaxis()->GetXmin() > -1.75 || oneNeutronTowerMap->GetXaxis()->GetXmax() < 1.75 ||
                 oneNeutronTowerMap->GetYaxis()->GetXmin() > -1.75 || oneNeutronTowerMap->GetYaxis()->GetXmax() < 1.75 ||
                 oneNeutronTowerMap->GetZaxis()->GetXmin() != 0. || oneNeutronTowerMap->GetZaxis()->GetXmax() != 8.) {
        LOGF(warn, "Incompatible 1n tower map dimensions for run %d; step skipped", runNumber);
        oneNeutronTowerMap = nullptr;
      }
    }
    if (calibrationStage.value == 2 && getSpatialMode() == 2) {
      if (!confCentralitySpatialLocalFile.value.empty()) {
        if (cfgCentralitySpatialLocalRun.value > 0 && cfgCentralitySpatialLocalRun.value != runNumber) {
          LOGF(fatal, "Local spatial map run mismatch: configured %d, current %d", cfgCentralitySpatialLocalRun.value, runNumber);
        }
        TFile localFile(confCentralitySpatialLocalFile.value.c_str(), "READ");
        if (!localFile.IsZombie()) {
          TH3D* localHist = nullptr;
          localFile.GetObject("CentralitySpatialTowerFactors", localHist);
          if (localHist) {
            localCentralitySpatialMap.reset(static_cast<TH3D*>(localHist->Clone("LocalCentSpatialTowerFactors")));
            localCentralitySpatialMap->SetDirectory(nullptr);
            centralitySpatialTowerMap = localCentralitySpatialMap.get();
          }
        }
      } else {
        centralitySpatialTowerMap = ccdb->getForTimeStamp<TH3D>(confCentralitySpatialPath.value, timestamp);
      }
      if (!centralitySpatialTowerMap || centralitySpatialTowerMap->GetNbinsZ() != 8 * kCentSpatialClasses ||
          centralitySpatialTowerMap->GetZaxis()->GetXmin() != 0. ||
          centralitySpatialTowerMap->GetZaxis()->GetXmax() != 8 * kCentSpatialClasses ||
          centralitySpatialTowerMap->GetXaxis()->GetXmin() > -1. || centralitySpatialTowerMap->GetXaxis()->GetXmax() < 1. ||
          centralitySpatialTowerMap->GetYaxis()->GetXmin() > -1. || centralitySpatialTowerMap->GetYaxis()->GetXmax() < 1.) {
        LOGF(fatal, "Missing/incompatible centrality-dependent cross-coordinate map for run %d at %llu", runNumber, timestamp);
      }
    }
    if (calibrationStage.value == 2 && useQRecentering.value) {
      const std::string qRecenteringPath = confQRecenteringStructured71Path.value;
      if (!qRecenteringPath.empty()) {
        qRecenteringProfile = ccdb->getForTimeStamp<TH3D>(qRecenteringPath, timestamp);
      }
      if (!qRecenteringProfile) {
        LOGF(warn, "No ZDC Q-recentering calibration found for run %d at timestamp %llu", runNumber, timestamp);
      } else {
        const bool structuredMapOK = qRecenteringProfile->GetYaxis()->GetNbins() == kQRecenteringStructuredNFeatures &&
                                     std::string(qRecenteringProfile->GetYaxis()->GetBinLabel(kQRecenteringMultidimNFeatures + 1)) == "(T--0.8)_+^3" &&
                                     std::string(qRecenteringProfile->GetYaxis()->GetBinLabel(kQRecenteringStructuredNFeatures)) == "(Vz-+0.8)_+^3";
        if (!structuredMapOK) {
          LOGF(warn, "Invalid structured71 map %s run %d: recentering skipped", qRecenteringPath.c_str(), runNumber);
          qRecenteringProfile = nullptr;
        }
      }
    }
  }

  float getTimeFromSOR(uint64_t globalBC) const
  {
    if (bcSOR < 0 || globalBC < static_cast<uint64_t>(bcSOR)) {
      return -1.f;
    }
    const auto deltaBC = static_cast<double>(globalBC - static_cast<uint64_t>(bcSOR));
    return static_cast<float>(deltaBC * o2::constants::lhc::LHCBunchSpacingNS * 1.e-9 / 3600.);
  }

  double getGain(float timeFromSOR, int channel, bool& gainOK) const
  {
    if (calibrationStage.value != 2 || !useGainCallib.value) {
      return 1.0;
    }
    if (!gainprofile || timeFromSOR < 0.f) {
      gainOK = false;
      return 1.0;
    }
    const double gain = gainprofile->GetBinContent(gainprofile->FindBin(timeFromSOR + 1.e-7, channel + 0.5));
    if (!std::isfinite(gain) || gain <= 0.0) {
      gainOK = false;
      return 1.0;
    }
    return gain;
  }

  std::array<double, kQRecenteringMultidimNFeatures> makeQRecenteringMultidimFeatures(float timeFromSOR, float vx, float vy, float vz) const
  {
    const double runHours = (eorTimestamp > sorTimestamp) ? static_cast<double>(eorTimestamp - sorTimestamp) * 1.e-3 / 3600.0 : static_cast<double>(cfgMaxRunHours.value);
    const double t = (runHours > 0.0) ? (2.0 * static_cast<double>(timeFromSOR) / runHours - 1.0) : 0.0;
    const double x = static_cast<double>(vx) / 0.1;
    const double y = static_cast<double>(vy) / 0.1;
    const double z = static_cast<double>(vz) / 10.0;
    const double t2 = t * t;
    const double x2 = x * x;
    const double y2 = y * y;
    const double z2 = z * z;

    return {1.0, t, x, y, z,
            t2, x2, y2, z2, t * x, t * y, t * z, x * y, x * z, y * z,
            t2 * t, x2 * x, y2 * y, z2 * z,
            t2 * x, t2 * y, t2 * z, x2 * t, x2 * y, x2 * z,
            y2 * t, y2 * x, y2 * z, z2 * t, z2 * x, z2 * y,
            t * x * y, t * x * z, t * y * z, x * y * z};
  }

  std::array<double, kQRecenteringStructuredNFeatures> makeQRecenteringStructuredFeatures(float timeFromSOR, float vx, float vy, float vz) const
  {
    const auto baseFeatures = makeQRecenteringMultidimFeatures(timeFromSOR, vx, vy, vz);
    std::array<double, kQRecenteringStructuredNFeatures> features{};

    for (int i = 0; i < kQRecenteringMultidimNFeatures; ++i) {
      features[i] = baseFeatures[i];
    }

    const double runHours = (eorTimestamp > sorTimestamp) ? static_cast<double>(eorTimestamp - sorTimestamp) * 1.e-3 / 3600.0 : static_cast<double>(cfgMaxRunHours.value);
    const double t = (runHours > 0.0) ? (2.0 * static_cast<double>(timeFromSOR) / runHours - 1.0) : 0.0;
    const double x = static_cast<double>(vx) / 0.1;
    const double y = static_cast<double>(vy) / 0.1;
    const double z = static_cast<double>(vz) / 10.0;

    // Internal knots in the normalized [-1,1] coordinate. Together with the
    // existing 1,u,u^2,u^3 terms this is a cubic regression-spline basis.
    constexpr std::array<double, kQRecenteringStructuredNKnots> knots = {
      -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};

    int index = kQRecenteringMultidimNFeatures;
    auto appendSplineTerms = [&](double value) {
      for (int i = 0; i < kQRecenteringStructuredNKnots; ++i) {
        const double delta = value - knots[i];
        features[index++] = delta > 0.0 ? delta * delta * delta : 0.0;
      }
    };

    appendSplineTerms(t);
    appendSplineTerms(x);
    appendSplineTerms(y);
    appendSplineTerms(z);

    return features;
  }

  using MyCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::Mults, o2::aod::FT0sCorrected, o2::aod::CentFT0Cs>;

  using ResponseTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection>;

  void process(MyCollisions::iterator const& collision, ResponseTracks const& tracks, o2::aod::FT0s const&, o2::aod::FV0As const&, BCsRun3 const&, o2::aod::Zdcs const&)
  {
    histos.fill(HIST("hEvtSelInfo"), 0.5);

    const float centrality = collision.centFT0C();
    const float vx = collision.posX();
    const float vy = collision.posY();
    const float vz = collision.posZ();
    auto bc = collision.foundBC_as<BCsRun3>();
    const int runNumber = bc.runNumber();
    const uint64_t globalBC = bc.globalBC();
    const uint64_t timestamp = bc.timestamp();

    updateRun(runNumber, timestamp);
    const float timeFromSOR = getTimeFromSOR(globalBC);

    auto fillTable = [&](bool trigger, float qxA, float qxC, float qyA, float qyC) {
      zdccaltable(collision.globalIndex(), trigger, runNumber, globalBC, centrality, vx, vy, vz, qxA, qxC, qyA, qyC);
    };

    if (!bc.has_zdc()) {
      fillTable(false, -999.f, -999.f, -999.f, -999.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 1.5);

    auto zdc = bc.zdc();
    auto znaEnergy = zdc.energySectorZNA();
    auto zncEnergy = zdc.energySectorZNC();
    const float znaEnergycommon = zdc.energyCommonZNA();
    const float zncEnergycommon = zdc.energyCommonZNC();

    if (!std::isfinite(znaEnergycommon) || !std::isfinite(zncEnergycommon) || znaEnergycommon <= 0.f || zncEnergycommon <= 0.f) {
      fillTable(false, -999.f, -999.f, -999.f, -999.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 2.5);

    // Classification uses common-channel energy only; no Q-angle or opposite-arm condition.
    const auto neutronCategory = [&](double e) -> int {
      if (e >= cfg1nMin.value && e < cfg1nMax.value)
        return 1;
      if (e >= cfg2nMin.value && e < cfg2nMax.value)
        return 2;
      if (e >= cfg3nMin.value && e < cfg3nMax.value)
        return 3;
      return 0;
    };
    const bool neutronQAEnabled = deriveSingleNeutronResponse.value || derive1nTowerResponse.value || getSpatialMode() == 1;
    const int nClassA = neutronQAEnabled ? neutronCategory(znaEnergycommon) : 0;
    const int nClassC = neutronQAEnabled ? neutronCategory(zncEnergycommon) : 0;
    bool validTowerA = true, validTowerC = true;
    for (int i = 0; i < 4; ++i) {
      validTowerA &= std::isfinite(znaEnergy[i]) && znaEnergy[i] > 0.f;
      validTowerC &= std::isfinite(zncEnergy[i]) && zncEnergy[i] > 0.f;
    }
    if (deriveSingleNeutronResponse.value) {
      // Fill before valid-tower return, so rejected neutron candidates are counted.
      if (nClassA) {
        histos.fill(HIST("NeutronQA/hSelectionStatus"), nClassA - .5, .5);
        if (validTowerA)
          histos.fill(HIST("NeutronQA/hSelectionStatus"), nClassA - .5, 2.5);
        if (validTowerA && validTowerC)
          histos.fill(HIST("NeutronQA/hSelectionStatus"), nClassA - .5, 1.5);
      }
      if (nClassC) {
        histos.fill(HIST("NeutronQA/hSelectionStatus"), nClassC + 2.5, .5);
        if (validTowerC)
          histos.fill(HIST("NeutronQA/hSelectionStatus"), nClassC + 2.5, 3.5);
        if (validTowerA && validTowerC)
          histos.fill(HIST("NeutronQA/hSelectionStatus"), nClassC + 2.5, 1.5);
      }
    }
    bool validTowers = true;
    for (int i = 0; i < 4; ++i) {
      if (!std::isfinite(znaEnergy[i]) || !std::isfinite(zncEnergy[i]) || znaEnergy[i] <= 0.f || zncEnergy[i] <= 0.f) {
        validTowers = false;
      }
    }
    if (!validTowers) {
      fillTable(false, -999.f, -999.f, -999.f, -999.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 3.5);

    if (rctCut.requireRCTFlagChecker && !rctChecker(collision)) {
      fillTable(false, -999.f, -999.f, -999.f, -999.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 4.5);

    const bool acceptedCollision = collision.sel8() && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && std::isfinite(centrality) && centrality >= cfgCutCentralityMin && centrality < cfgCutCentralityMax && std::abs(vz) < cfgCutVertex && collision.has_foundFT0();
    if (!acceptedCollision) {
      fillTable(false, -999.f, -999.f, -999.f, -999.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 5.5);

    if (timeFromSOR < 0.f || timeFromSOR >= cfgMaxRunHours.value) {
      fillTable(false, -999.f, -999.f, -999.f, -999.f);
      return;
    }

    constexpr std::array<double, 4> x = {-1.75, 1.75, -1.75, 1.75};
    constexpr std::array<double, 4> y = {-1.75, -1.75, 1.75, 1.75};

    double rawXA = 0.0;
    double rawYA = 0.0;
    double rawXC = 0.0;
    double rawYC = 0.0;
    double rawSumA = 0.0;
    double rawSumC = 0.0;
    for (int i = 0; i < 4; ++i) {
      rawXA -= znaEnergy[i] * x[i];
      rawYA += znaEnergy[i] * y[i];
      rawXC += zncEnergy[i] * x[i];
      rawYC += zncEnergy[i] * y[i];
      rawSumA += znaEnergy[i];
      rawSumC += zncEnergy[i];
    }
    rawXA /= rawSumA;
    rawYA /= rawSumA;
    rawXC /= rawSumC;
    rawYC /= rawSumC;
    const double phiRawA = std::atan2(rawYA, rawXA);
    const double phiRawC = std::atan2(rawYC, rawXC);
    histos.fill(HIST("PhiQA/hPhiRawZNA"), phiRawA);
    histos.fill(HIST("PhiQA/hPhiRawZNC"), phiRawC);
    histos.fill(HIST("PhiQA/hPhiRawVsCentralityZNA"), centrality, phiRawA);
    histos.fill(HIST("PhiQA/hPhiRawVsCentralityZNC"), centrality, phiRawC);

    if (calibrationStage.value == 1) {
      std::array<double, 15> momentsA{};
      std::array<double, 15> momentsC{};
      int moment = 0;
      for (int i = 0; i < 4; ++i) {
        for (int j = i; j < 4; ++j) {
          momentsA[moment] = static_cast<double>(znaEnergy[i]) * static_cast<double>(znaEnergy[j]);
          momentsC[moment] = static_cast<double>(zncEnergy[i]) * static_cast<double>(zncEnergy[j]);
          ++moment;
        }
      }
      for (int i = 0; i < 4; ++i) {
        momentsA[10 + i] = static_cast<double>(znaEnergycommon) * static_cast<double>(znaEnergy[i]);
        momentsC[10 + i] = static_cast<double>(zncEnergycommon) * static_cast<double>(zncEnergy[i]);
      }
      momentsA[14] = 1.0;
      momentsC[14] = 1.0;
      for (int i = 0; i < 15; ++i) {
        histos.fill(HIST("GainCalibration/hGainMomentsZNA"), timeFromSOR, i + 0.5, momentsA[i]);
        histos.fill(HIST("GainCalibration/hGainMomentsZNC"), timeFromSOR, i + 0.5, momentsC[i]);
      }
    }

    bool gainOK = true;
    std::array<double, 4> znaCorr{};
    std::array<double, 4> zncCorr{};
    std::array<double, 8> gains{};
    for (int i = 0; i < 8; ++i) {
      gains[i] = getGain(timeFromSOR, i, gainOK);
    }
    if (!gainOK) {
      // The whole gain step is skipped, never only a subset of channels.
      // Recentering stays independently enabled; the new tower fit requires valid gain.
      LOGF(debug, "Gain coefficients unavailable/invalid for run %d; using unity gains for this event", runNumber);
      if (getSpatialMode() == 2 || deriveCentralitySpatial.value) {
        LOGF(fatal, "Centrality spatial correction/derivation requires a valid stage-1 gain map for run %d", runNumber);
      }
      gains.fill(1.0);
    }
    for (int i = 0; i < 4; ++i) {
      znaCorr[i] = gains[i] * znaEnergy[i];
      zncCorr[i] = gains[i + 4] * zncEnergy[i];
    }
    histos.fill(HIST("hEvtSelInfo"), 6.5);
    // Independent 1n/2n/3n selections based only on each arm's common energy.
    // The existing event Q reconstruction still requires all eight towers valid.
    if (deriveSingleNeutronResponse.value) {
      histos.fill(HIST("NeutronQA/hEnergyVsCentA"), centrality, znaEnergycommon);
      histos.fill(HIST("NeutronQA/hEnergyVsCentC"), centrality, zncEnergycommon);
    }

    double gainXA = 0.0;
    double gainYA = 0.0;
    double gainXC = 0.0;
    double gainYC = 0.0;
    double gainSumA = 0.0;
    double gainSumC = 0.0;
    for (int i = 0; i < 4; ++i) {
      gainXA -= znaCorr[i] * x[i];
      gainYA += znaCorr[i] * y[i];
      gainXC += zncCorr[i] * x[i];
      gainYC += zncCorr[i] * y[i];
      gainSumA += znaCorr[i];
      gainSumC += zncCorr[i];
    }
    gainXA /= gainSumA;
    gainYA /= gainSumA;
    gainXC /= gainSumC;
    gainYC /= gainSumC;
    const double phiGainA = std::atan2(gainYA, gainXA);
    const double phiGainC = std::atan2(gainYC, gainXC);
    histos.fill(HIST("PhiQA/hPhiAfterGainZNA"), phiGainA);
    histos.fill(HIST("PhiQA/hPhiAfterGainZNC"), phiGainC);
    histos.fill(HIST("PhiQA/hPhiAfterGainVsCentralityZNA"), centrality, phiGainA);
    histos.fill(HIST("PhiQA/hPhiAfterGainVsCentralityZNC"), centrality, phiGainC);

    const double crossLookupXA = (znaCorr[3] - znaCorr[0]) / (znaCorr[3] + znaCorr[0]);
    const double crossLookupYA = (znaCorr[2] - znaCorr[1]) / (znaCorr[2] + znaCorr[1]);
    const double crossLookupXC = (zncCorr[3] - zncCorr[0]) / (zncCorr[3] + zncCorr[0]);
    const double crossLookupYC = (zncCorr[2] - zncCorr[1]) / (zncCorr[2] + zncCorr[1]);
    // Hadronic, centrality-resolved regression input. All accepted events,
    // with gain-only cross coordinates and the common channel reference.
    if (deriveCentralitySpatial.value || getSpatialMode() == 2) {
      const int centClass = getCentSpatialClass(centrality);
      if (centClass >= 0) {
        const double cc = centClass + 0.5;
        histos.fill(HIST("CentSpatial/pGainResponse0ZA"), crossLookupXA, crossLookupYA, cc, znaCorr[0] / znaEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse1ZA"), crossLookupXA, crossLookupYA, cc, znaCorr[1] / znaEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse2ZA"), crossLookupXA, crossLookupYA, cc, znaCorr[2] / znaEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse3ZA"), crossLookupXA, crossLookupYA, cc, znaCorr[3] / znaEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse4ZA"), crossLookupXA, crossLookupYA, cc, (znaCorr[0] + znaCorr[1] + znaCorr[2] + znaCorr[3]) / znaEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse0ZC"), crossLookupXC, crossLookupYC, cc, zncCorr[0] / zncEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse1ZC"), crossLookupXC, crossLookupYC, cc, zncCorr[1] / zncEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse2ZC"), crossLookupXC, crossLookupYC, cc, zncCorr[2] / zncEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse3ZC"), crossLookupXC, crossLookupYC, cc, zncCorr[3] / zncEnergycommon);
        histos.fill(HIST("CentSpatial/pGainResponse4ZC"), crossLookupXC, crossLookupYC, cc, (zncCorr[0] + zncCorr[1] + zncCorr[2] + zncCorr[3]) / zncEnergycommon);

        if (deriveCentralitySpatial.value && gainOK) {
          uint64_t splitKey = globalBC + 0xd6e8feb86659fd93ULL + static_cast<uint64_t>(collision.globalIndex());
          splitKey = (splitKey ^ (splitKey >> 30)) * 0xbf58476d1ce4e5b9ULL;
          splitKey = (splitKey ^ (splitKey >> 27)) * 0x94d049bb133111ebULL;
          splitKey ^= (splitKey >> 31);
          const bool isValidation = (splitKey & 1ULL) != 0ULL;
          const auto fillMoments = [&](bool sideA, bool valid, double cx, double cy, const std::array<double, 4>& energy, double common) {
            std::array<double, kCentSpatialMomentCount> moments{};
            std::array<double, 4> z{};
            for (int i = 0; i < 4; ++i)
              z[i] = energy[i] / common;
            int m = 0;
            for (int i = 0; i < 4; ++i)
              for (int j = i; j < 4; ++j)
                moments[m++] = z[i] * z[j];
            for (int i = 0; i < 4; ++i)
              moments[10 + i] = z[i];
            moments[14] = 1.;
            if (sideA) {
              if (valid) {
                histos.fill(HIST("CentSpatial/hMomValidZA_00"), cx, cy, cc, moments[0]);
                histos.fill(HIST("CentSpatial/hMomValidZA_01"), cx, cy, cc, moments[1]);
                histos.fill(HIST("CentSpatial/hMomValidZA_02"), cx, cy, cc, moments[2]);
                histos.fill(HIST("CentSpatial/hMomValidZA_03"), cx, cy, cc, moments[3]);
                histos.fill(HIST("CentSpatial/hMomValidZA_04"), cx, cy, cc, moments[4]);
                histos.fill(HIST("CentSpatial/hMomValidZA_05"), cx, cy, cc, moments[5]);
                histos.fill(HIST("CentSpatial/hMomValidZA_06"), cx, cy, cc, moments[6]);
                histos.fill(HIST("CentSpatial/hMomValidZA_07"), cx, cy, cc, moments[7]);
                histos.fill(HIST("CentSpatial/hMomValidZA_08"), cx, cy, cc, moments[8]);
                histos.fill(HIST("CentSpatial/hMomValidZA_09"), cx, cy, cc, moments[9]);
                histos.fill(HIST("CentSpatial/hMomValidZA_10"), cx, cy, cc, moments[10]);
                histos.fill(HIST("CentSpatial/hMomValidZA_11"), cx, cy, cc, moments[11]);
                histos.fill(HIST("CentSpatial/hMomValidZA_12"), cx, cy, cc, moments[12]);
                histos.fill(HIST("CentSpatial/hMomValidZA_13"), cx, cy, cc, moments[13]);
                histos.fill(HIST("CentSpatial/hMomValidZA_14"), cx, cy, cc, moments[14]);
              } else {
                histos.fill(HIST("CentSpatial/hMomTrainZA_00"), cx, cy, cc, moments[0]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_01"), cx, cy, cc, moments[1]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_02"), cx, cy, cc, moments[2]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_03"), cx, cy, cc, moments[3]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_04"), cx, cy, cc, moments[4]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_05"), cx, cy, cc, moments[5]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_06"), cx, cy, cc, moments[6]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_07"), cx, cy, cc, moments[7]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_08"), cx, cy, cc, moments[8]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_09"), cx, cy, cc, moments[9]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_10"), cx, cy, cc, moments[10]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_11"), cx, cy, cc, moments[11]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_12"), cx, cy, cc, moments[12]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_13"), cx, cy, cc, moments[13]);
                histos.fill(HIST("CentSpatial/hMomTrainZA_14"), cx, cy, cc, moments[14]);
              }
            } else {
              if (valid) {
                histos.fill(HIST("CentSpatial/hMomValidZC_00"), cx, cy, cc, moments[0]);
                histos.fill(HIST("CentSpatial/hMomValidZC_01"), cx, cy, cc, moments[1]);
                histos.fill(HIST("CentSpatial/hMomValidZC_02"), cx, cy, cc, moments[2]);
                histos.fill(HIST("CentSpatial/hMomValidZC_03"), cx, cy, cc, moments[3]);
                histos.fill(HIST("CentSpatial/hMomValidZC_04"), cx, cy, cc, moments[4]);
                histos.fill(HIST("CentSpatial/hMomValidZC_05"), cx, cy, cc, moments[5]);
                histos.fill(HIST("CentSpatial/hMomValidZC_06"), cx, cy, cc, moments[6]);
                histos.fill(HIST("CentSpatial/hMomValidZC_07"), cx, cy, cc, moments[7]);
                histos.fill(HIST("CentSpatial/hMomValidZC_08"), cx, cy, cc, moments[8]);
                histos.fill(HIST("CentSpatial/hMomValidZC_09"), cx, cy, cc, moments[9]);
                histos.fill(HIST("CentSpatial/hMomValidZC_10"), cx, cy, cc, moments[10]);
                histos.fill(HIST("CentSpatial/hMomValidZC_11"), cx, cy, cc, moments[11]);
                histos.fill(HIST("CentSpatial/hMomValidZC_12"), cx, cy, cc, moments[12]);
                histos.fill(HIST("CentSpatial/hMomValidZC_13"), cx, cy, cc, moments[13]);
                histos.fill(HIST("CentSpatial/hMomValidZC_14"), cx, cy, cc, moments[14]);
              } else {
                histos.fill(HIST("CentSpatial/hMomTrainZC_00"), cx, cy, cc, moments[0]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_01"), cx, cy, cc, moments[1]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_02"), cx, cy, cc, moments[2]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_03"), cx, cy, cc, moments[3]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_04"), cx, cy, cc, moments[4]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_05"), cx, cy, cc, moments[5]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_06"), cx, cy, cc, moments[6]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_07"), cx, cy, cc, moments[7]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_08"), cx, cy, cc, moments[8]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_09"), cx, cy, cc, moments[9]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_10"), cx, cy, cc, moments[10]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_11"), cx, cy, cc, moments[11]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_12"), cx, cy, cc, moments[12]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_13"), cx, cy, cc, moments[13]);
                histos.fill(HIST("CentSpatial/hMomTrainZC_14"), cx, cy, cc, moments[14]);
              }
            }
          };
          fillMoments(true, isValidation, crossLookupXA, crossLookupYA, znaCorr, znaEnergycommon);
          fillMoments(false, isValidation, crossLookupXC, crossLookupYC, zncCorr, zncEnergycommon);
        }
      }
    }

    if (deriveSingleNeutronResponse.value) {
      if (nClassA) {
        histos.fill(HIST("NeutronQA/hTowerFractionA1"), znaEnergycommon, znaCorr[0] / znaEnergycommon);
        histos.fill(HIST("NeutronQA/hTowerFractionA2"), znaEnergycommon, znaCorr[1] / znaEnergycommon);
        histos.fill(HIST("NeutronQA/hTowerFractionA3"), znaEnergycommon, znaCorr[2] / znaEnergycommon);
        histos.fill(HIST("NeutronQA/hTowerFractionA4"), znaEnergycommon, znaCorr[3] / znaEnergycommon);
      }
      if (nClassC) {
        histos.fill(HIST("NeutronQA/hTowerFractionC1"), zncEnergycommon, zncCorr[0] / zncEnergycommon);
        histos.fill(HIST("NeutronQA/hTowerFractionC2"), zncEnergycommon, zncCorr[1] / zncEnergycommon);
        histos.fill(HIST("NeutronQA/hTowerFractionC3"), zncEnergycommon, zncCorr[2] / zncEnergycommon);
        histos.fill(HIST("NeutronQA/hTowerFractionC4"), zncEnergycommon, zncCorr[3] / zncEnergycommon);
      }
    }

    if (deriveSingleNeutronResponse.value) {
      if (nClassA == 1) {
        histos.fill(HIST("NeutronQA/hCross1nA"), crossLookupXA, crossLookupYA);
        histos.fill(HIST("NeutronQA/pGainTowerFractionVsCross1nA1"), crossLookupXA, crossLookupYA, znaCorr[0] / znaEnergycommon);
        histos.fill(HIST("NeutronQA/pGainTowerFractionVsCross1nA2"), crossLookupXA, crossLookupYA, znaCorr[1] / znaEnergycommon);
        histos.fill(HIST("NeutronQA/pGainTowerFractionVsCross1nA3"), crossLookupXA, crossLookupYA, znaCorr[2] / znaEnergycommon);
        histos.fill(HIST("NeutronQA/pGainTowerFractionVsCross1nA4"), crossLookupXA, crossLookupYA, znaCorr[3] / znaEnergycommon);
      }
      if (nClassC == 1) {
        histos.fill(HIST("NeutronQA/hCross1nC"), crossLookupXC, crossLookupYC);
        histos.fill(HIST("NeutronQA/pGainTowerFractionVsCross1nC1"), crossLookupXC, crossLookupYC, zncCorr[0] / zncEnergycommon);
        histos.fill(HIST("NeutronQA/pGainTowerFractionVsCross1nC2"), crossLookupXC, crossLookupYC, zncCorr[1] / zncEnergycommon);
        histos.fill(HIST("NeutronQA/pGainTowerFractionVsCross1nC3"), crossLookupXC, crossLookupYC, zncCorr[2] / zncEnergycommon);
        histos.fill(HIST("NeutronQA/pGainTowerFractionVsCross1nC4"), crossLookupXC, crossLookupYC, zncCorr[3] / zncEnergycommon);
      }
    }

    // Derive 1n response from JOINT gain-corrected tower energy fractions.
    // No absolute common-channel normalization enters this fit; common energy
    // only selects the 1n/2n/3n classes. Fit each arm independently.
    if (derive1nTowerResponse.value && calibrationStage.value == 2 && gainOK) {
      uint64_t splitKey = globalBC + 0x517cc1b727220a95ULL + static_cast<uint64_t>(collision.globalIndex());
      splitKey = (splitKey ^ (splitKey >> 30)) * 0xbf58476d1ce4e5b9ULL;
      splitKey = (splitKey ^ (splitKey >> 27)) * 0x94d049bb133111ebULL;
      splitKey ^= splitKey >> 31;
      const bool isValidation = (splitKey & 1ULL) != 0ULL;
      const auto fillRelative = [&](bool sideA, int neutronClass, const std::array<double, 4>& energy, double sum) {
        if (neutronClass == 0 || !(sum > 0.)) {
          return;
        }
        const double f1 = energy[0] / sum, f2 = energy[1] / sum, f3 = energy[2] / sum;
        if (sideA) {
          if (neutronClass == 1) {
            if (isValidation) {
              histos.fill(HIST("OneNeutronTowerCalibration/hRelative1nA_Validation"), f1, f2, f3);
            } else {
              histos.fill(HIST("OneNeutronTowerCalibration/hRelative1nA_Train"), f1, f2, f3);
            }
          } else if (neutronClass == 2) {
            histos.fill(HIST("OneNeutronTowerCalibration/hRelative2nA"), f1, f2, f3);
          } else if (neutronClass == 3) {
            histos.fill(HIST("OneNeutronTowerCalibration/hRelative3nA"), f1, f2, f3);
          }
        } else {
          if (neutronClass == 1) {
            if (isValidation) {
              histos.fill(HIST("OneNeutronTowerCalibration/hRelative1nC_Validation"), f1, f2, f3);
            } else {
              histos.fill(HIST("OneNeutronTowerCalibration/hRelative1nC_Train"), f1, f2, f3);
            }
          } else if (neutronClass == 2) {
            histos.fill(HIST("OneNeutronTowerCalibration/hRelative2nC"), f1, f2, f3);
          } else if (neutronClass == 3) {
            histos.fill(HIST("OneNeutronTowerCalibration/hRelative3nC"), f1, f2, f3);
          }
        }
      };
      fillRelative(true, nClassA, znaCorr, gainSumA);
      fillRelative(false, nClassC, zncCorr, gainSumC);
    }

    // Apply relative tower corrections atomically for each arm at the GAIN-only
    // reconstructed Q position (not recentered Q; do not iterate the lookup).
    // Each factor is positive; a common multiplicative factor cancels in Q.
    bool applied1nA = false, applied1nC = false;
    if (calibrationStage.value == 2 && getSpatialMode() == 1 && gainOK && oneNeutronTowerMap) {
      const auto applyTowerResponse = [&](std::array<double, 4>& energies, double qx, double qy, int zOffset) -> bool {
        if (!std::isfinite(qx) || !std::isfinite(qy) ||
            qx < oneNeutronTowerMap->GetXaxis()->GetXmin() || qx >= oneNeutronTowerMap->GetXaxis()->GetXmax() ||
            qy < oneNeutronTowerMap->GetYaxis()->GetXmin() || qy >= oneNeutronTowerMap->GetYaxis()->GetXmax())
          return false;
        std::array<double, 4> factor{};
        for (int i = 0; i < 4; ++i) {
          factor[i] = oneNeutronTowerMap->GetBinContent(oneNeutronTowerMap->FindBin(qx, qy, zOffset + i + 0.5));
          if (!std::isfinite(factor[i]) || factor[i] <= 0. || factor[i] > 2.)
            return false;
        }
        for (int i = 0; i < 4; ++i)
          energies[i] *= factor[i];
        return true;
      };
      applied1nA = applyTowerResponse(znaCorr, gainXA, gainYA, 0);
      applied1nC = applyTowerResponse(zncCorr, gainXC, gainYC, 4);
    }
    if (calibrationStage.value == 2 &&
        (deriveSingleNeutronResponse.value || derive1nTowerResponse.value || getSpatialMode() == 1)) {
      histos.fill(HIST("OneNeutronTowerQA/hCorrectionStatus"),
                  0.5 + static_cast<int>(applied1nA) + 2 * static_cast<int>(applied1nC));
    }

    // Independent mode 2: centrality-dependent correction in DETECTOR cross space.
    // Lookup once using gain-only E; no Q-dependent fitting or recentering here.
    bool appliedCentA = false, appliedCentC = false;
    if (calibrationStage.value == 2 && getSpatialMode() == 2 && gainOK && centralitySpatialTowerMap) {
      const int centClass = getCentSpatialClass(centrality);
      if (centClass >= 0) {
        const auto applyCentFactors = [&](std::array<double, 4>& energy, double cx, double cy, int sideOffset) -> bool {
          if (!std::isfinite(cx) || !std::isfinite(cy) ||
              cx < centralitySpatialTowerMap->GetXaxis()->GetXmin() || cx >= centralitySpatialTowerMap->GetXaxis()->GetXmax() ||
              cy < centralitySpatialTowerMap->GetYaxis()->GetXmin() || cy >= centralitySpatialTowerMap->GetYaxis()->GetXmax())
            return false;
          std::array<double, 4> factors{};
          for (int i = 0; i < 4; ++i) {
            const double packZ = centClass * 8 + sideOffset + i + 0.5;
            factors[i] = centralitySpatialTowerMap->GetBinContent(centralitySpatialTowerMap->FindBin(cx, cy, packZ));
            if (!std::isfinite(factors[i]) || !(factors[i] > 0.) ||
                std::abs(factors[i] - 1.) > cfgCentralitySpatialMaxChange.value + 1.e-5) {
              LOGF(fatal, "Invalid centrality spatial factor %.6f: run %d centrality %.3f X %.3f Y %.3f tower %d", factors[i], currentRunNumber, centrality, cx, cy, i);
            }
          }
          for (int i = 0; i < 4; ++i)
            energy[i] *= factors[i];
          return true;
        };
        appliedCentA = applyCentFactors(znaCorr, crossLookupXA, crossLookupYA, 0);
        appliedCentC = applyCentFactors(zncCorr, crossLookupXC, crossLookupYC, 4);
      }
    }
    if (getSpatialMode() == 2) {
      histos.fill(HIST("CentSpatial/hCorrectionStatus"), 0.5 + static_cast<int>(appliedCentA) + 2 * static_cast<int>(appliedCentC));
      if (getCentSpatialClass(centrality) >= 0 && (!appliedCentA || !appliedCentC)) {
        LOGF(fatal, "Mode2 spatial map failed for selected event at centrality %.3f run %d", centrality, runNumber);
      }
    }

    double sumEnergyA = 0.0;
    double sumEnergyC = 0.0;
    for (int i = 0; i < 4; ++i) {
      sumEnergyA += znaCorr[i];
      sumEnergyC += zncCorr[i];
    }

    const double ratioA = sumEnergyA / znaEnergycommon;
    const double ratioC = sumEnergyC / zncEnergycommon;
    const double crossXA = (znaCorr[3] - znaCorr[0]) / (znaCorr[3] + znaCorr[0]);
    const double crossYA = (znaCorr[2] - znaCorr[1]) / (znaCorr[2] + znaCorr[1]);
    const double crossXC = (zncCorr[3] - zncCorr[0]) / (zncCorr[3] + zncCorr[0]);
    const double crossYC = (zncCorr[2] - zncCorr[1]) / (zncCorr[2] + zncCorr[1]);

    double centroidXA = 0.0;
    double centroidYA = 0.0;
    double centroidXC = 0.0;
    double centroidYC = 0.0;
    for (int i = 0; i < 4; ++i) {
      centroidXA -= znaCorr[i] * x[i];
      centroidYA += znaCorr[i] * y[i];
      centroidXC += zncCorr[i] * x[i];
      centroidYC += zncCorr[i] * y[i];
    }
    centroidXA /= sumEnergyA;
    centroidYA /= sumEnergyA;
    centroidXC /= sumEnergyC;
    centroidYC /= sumEnergyC;
    const double phiTowerCorrectedA = std::atan2(centroidYA, centroidXA);
    const double phiTowerCorrectedC = std::atan2(centroidYC, centroidXC);
    if (deriveCentralitySpatial.value || getSpatialMode() == 2) {
      const int cclass = getCentSpatialClass(centrality);
      if (cclass >= 0) {
        const double zz = cclass + 0.5;
        histos.fill(HIST("CentSpatial/hQGainZA"), gainXA, gainYA, zz);
        histos.fill(HIST("CentSpatial/hQGainZC"), gainXC, gainYC, zz);
        histos.fill(HIST("CentSpatial/hQSpatialZA"), centroidXA, centroidYA, zz);
        histos.fill(HIST("CentSpatial/hQSpatialZC"), centroidXC, centroidYC, zz);
        const double dxA = centroidXA - gainXA, dyA = centroidYA - gainYA;
        const double dxC = centroidXC - gainXC, dyC = centroidYC - gainYC;
        histos.fill(HIST("CentSpatial/pDeltaQxZA"), centrality, dxA);
        histos.fill(HIST("CentSpatial/pDeltaQyZA"), centrality, dyA);
        histos.fill(HIST("CentSpatial/pDeltaQMagZA"), centrality, std::hypot(dxA, dyA));
        histos.fill(HIST("CentSpatial/pDeltaQxZC"), centrality, dxC);
        histos.fill(HIST("CentSpatial/pDeltaQyZC"), centrality, dyC);
        histos.fill(HIST("CentSpatial/pDeltaQMagZC"), centrality, std::hypot(dxC, dyC));
        // Response ratio after application, binned by the ORIGINAL gain-only cross coordinates.
        histos.fill(HIST("CentSpatial/pCorrectedTotalZA"), crossLookupXA, crossLookupYA, zz, sumEnergyA / znaEnergycommon);
        histos.fill(HIST("CentSpatial/pCorrectedTotalZC"), crossLookupXC, crossLookupYC, zz, sumEnergyC / zncEnergycommon);
      }
    }
    histos.fill(HIST("PhiQA/hPhiAfterTowerCorrectionZNA"), phiTowerCorrectedA);
    histos.fill(HIST("PhiQA/hPhiAfterTowerCorrectionZNC"), phiTowerCorrectedC);
    histos.fill(HIST("PhiQA/hPhiAfterTowerCorrectionVsCentralityZNA"), centrality, phiTowerCorrectedA);
    histos.fill(HIST("PhiQA/hPhiAfterTowerCorrectionVsCentralityZNC"), centrality, phiTowerCorrectedC);

    histos.fill(HIST("GainQA/hSumOverCommonVsTimeZNA"), timeFromSOR, ratioA);
    histos.fill(HIST("GainQA/hSumOverCommonVsTimeZNC"), timeFromSOR, ratioC);
    for (int i = 0; i < 4; ++i) {
      histos.fill(HIST("GainQA/pTowerOverCommonVsTimeZNA"), timeFromSOR, i + 0.5, znaCorr[i] / znaEnergycommon);
      histos.fill(HIST("GainQA/pTowerOverCommonVsTimeZNC"), timeFromSOR, i + 0.5, zncCorr[i] / zncEnergycommon);
    }
    histos.fill(HIST("TowerResponseQA/pCrossVsTimeZNA"), timeFromSOR, 0.5, crossXA);
    histos.fill(HIST("TowerResponseQA/pCrossVsTimeZNA"), timeFromSOR, 1.5, crossYA);
    histos.fill(HIST("TowerResponseQA/pCrossVsTimeZNC"), timeFromSOR, 0.5, crossXC);
    histos.fill(HIST("TowerResponseQA/pCrossVsTimeZNC"), timeFromSOR, 1.5, crossYC);
    histos.fill(HIST("TowerResponseQA/pResponseVsCrossXYZNA"), crossLookupXA, crossLookupYA, ratioA);
    histos.fill(HIST("TowerResponseQA/pResponseVsCrossXYZNC"), crossLookupXC, crossLookupYC, ratioC);
    histos.fill(HIST("TowerResponseQA/pCentroidVsTimeZNA"), timeFromSOR, 0.5, centroidXA);
    histos.fill(HIST("TowerResponseQA/pCentroidVsTimeZNA"), timeFromSOR, 1.5, centroidYA);
    histos.fill(HIST("TowerResponseQA/pCentroidVsTimeZNC"), timeFromSOR, 0.5, centroidXC);
    histos.fill(HIST("TowerResponseQA/pCentroidVsTimeZNC"), timeFromSOR, 1.5, centroidYC);
    histos.fill(HIST("TowerResponseQA/hCrossXYZNA"), crossLookupXA, crossLookupYA);
    histos.fill(HIST("TowerResponseQA/hCrossXYZNC"), crossLookupXC, crossLookupYC);

    histos.fill(HIST("TowerResponseQA/pTower1OverCommonVsCrossXYZNA"), crossLookupXA, crossLookupYA, znaCorr[0] / znaEnergycommon);
    histos.fill(HIST("TowerResponseQA/pTower2OverCommonVsCrossXYZNA"), crossLookupXA, crossLookupYA, znaCorr[1] / znaEnergycommon);
    histos.fill(HIST("TowerResponseQA/pTower3OverCommonVsCrossXYZNA"), crossLookupXA, crossLookupYA, znaCorr[2] / znaEnergycommon);
    histos.fill(HIST("TowerResponseQA/pTower4OverCommonVsCrossXYZNA"), crossLookupXA, crossLookupYA, znaCorr[3] / znaEnergycommon);

    histos.fill(HIST("TowerResponseQA/pTower1OverCommonVsCrossXYZNC"), crossLookupXC, crossLookupYC, zncCorr[0] / zncEnergycommon);
    histos.fill(HIST("TowerResponseQA/pTower2OverCommonVsCrossXYZNC"), crossLookupXC, crossLookupYC, zncCorr[1] / zncEnergycommon);
    histos.fill(HIST("TowerResponseQA/pTower3OverCommonVsCrossXYZNC"), crossLookupXC, crossLookupYC, zncCorr[2] / zncEnergycommon);
    histos.fill(HIST("TowerResponseQA/pTower4OverCommonVsCrossXYZNC"), crossLookupXC, crossLookupYC, zncCorr[3] / zncEnergycommon);
    const double alphaZDC = 0.395;
    double qxZDCA = 0.0;
    double qxZDCC = 0.0;
    double qyZDCA = 0.0;
    double qyZDCC = 0.0;
    double sumA = 0.0;
    double sumC = 0.0;

    for (int i = 0; i < 4; ++i) {
      double amplA = znaCorr[i];
      double amplC = zncCorr[i];
      if (followpub) {
        amplA = std::pow(amplA, alphaZDC);
        amplC = std::pow(amplC, alphaZDC);
      }
      qxZDCA -= amplA * x[i];
      qyZDCA += amplA * y[i];
      qxZDCC += amplC * x[i];
      qyZDCC += amplC * y[i];
      sumA += amplA;
      sumC += amplC;
    }

    if (sumA <= 1.e-4 || sumC <= 1.e-4) {
      fillTable(false, -999.f, -999.f, -999.f, -999.f);
      return;
    }

    double cZNA = 1.0;
    double cZNC = 1.0;
    if (usecfactor) {
      const double beamEne = 5.36 * 0.5;
      const double nSpecnA = sumEnergyA / beamEne;
      const double nSpecnC = sumEnergyC / beamEne;
      cZNA = 1.89358 - 0.71262 / (nSpecnA + 0.71789);
      cZNC = 1.89358 - 0.71262 / (nSpecnC + 0.71789);
    }

    qxZDCA = cZNA * qxZDCA / sumA;
    qyZDCA = cZNA * qyZDCA / sumA;
    qxZDCC = cZNC * qxZDCC / sumC;
    qyZDCC = cZNC * qyZDCC / sumC;

    if (!std::isfinite(qxZDCA) || !std::isfinite(qyZDCA) || !std::isfinite(qxZDCC) || !std::isfinite(qyZDCC)) {
      fillTable(false, -999.f, -999.f, -999.f, -999.f);
      return;
    }

    // Fixed before/after distributions for each neutron class, independent of
    // whether structured71 is enabled. No class-specific recentering is used.
    if (calibrationStage.value == 2 &&
        (deriveSingleNeutronResponse.value || derive1nTowerResponse.value || getSpatialMode() == 1)) {
      const auto fillTowerQ = [&](int neutronClass, bool sideA, double qxBefore, double qyBefore,
                                  double qxAfter, double qyAfter) {
        if (sideA) {
          if (neutronClass == 1) {
            histos.fill(HIST("OneNeutronTowerQA/hBefore1nA"), qxBefore, qyBefore);
            histos.fill(HIST("OneNeutronTowerQA/hAfter1nA"), qxAfter, qyAfter);
          }
          if (neutronClass == 2) {
            histos.fill(HIST("OneNeutronTowerQA/hBefore2nA"), qxBefore, qyBefore);
            histos.fill(HIST("OneNeutronTowerQA/hAfter2nA"), qxAfter, qyAfter);
          }
          if (neutronClass == 3) {
            histos.fill(HIST("OneNeutronTowerQA/hBefore3nA"), qxBefore, qyBefore);
            histos.fill(HIST("OneNeutronTowerQA/hAfter3nA"), qxAfter, qyAfter);
          }
        } else {
          if (neutronClass == 1) {
            histos.fill(HIST("OneNeutronTowerQA/hBefore1nC"), qxBefore, qyBefore);
            histos.fill(HIST("OneNeutronTowerQA/hAfter1nC"), qxAfter, qyAfter);
          }
          if (neutronClass == 2) {
            histos.fill(HIST("OneNeutronTowerQA/hBefore2nC"), qxBefore, qyBefore);
            histos.fill(HIST("OneNeutronTowerQA/hAfter2nC"), qxAfter, qyAfter);
          }
          if (neutronClass == 3) {
            histos.fill(HIST("OneNeutronTowerQA/hBefore3nC"), qxBefore, qyBefore);
            histos.fill(HIST("OneNeutronTowerQA/hAfter3nC"), qxAfter, qyAfter);
          }
        }
      };
      fillTowerQ(nClassA, true, gainXA, gainYA, qxZDCA, qyZDCA);
      fillTowerQ(nClassC, false, gainXC, gainYC, qxZDCC, qyZDCC);
    }

    std::array<double, 4> qValues = {qxZDCA, qyZDCA, qxZDCC, qyZDCC};

    for (int i = 0; i < 4; ++i) {
      const double component = i + 0.5;
      histos.fill(HIST("QRecenteringQA/pQBeforeVsCentrality"), centrality, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQBeforeVsTime"), timeFromSOR, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQBeforeVsVx"), vx, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQBeforeVsVy"), vy, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQBeforeVsVz"), vz, component, qValues[i]);
    }

    histos.fill(HIST("QCorrelationQA/pQACorrelationBeforeVsCentrality"), centrality, 0.5, qValues[0] * qValues[2]);
    histos.fill(HIST("QCorrelationQA/pQACorrelationBeforeVsCentrality"), centrality, 1.5, qValues[1] * qValues[3]);
    histos.fill(HIST("QCorrelationQA/pQACorrelationBeforeVsCentrality"), centrality, 2.5, qValues[0] * qValues[3]);
    histos.fill(HIST("QCorrelationQA/pQACorrelationBeforeVsCentrality"), centrality, 3.5, qValues[1] * qValues[2]);

    if (calibrationStage.value == 2 && deriveQRecentering.value && !useQRecentering.value) {
      {
        // Structured71 regression moments. Events are deterministically split
        // 50/50 using a mixed global-BC key, avoiding any dependence on the
        // regression variables themselves.
        const auto structuredFeatures = makeQRecenteringStructuredFeatures(timeFromSOR, vx, vy, vz);
        uint64_t splitKey = globalBC + 0x9e3779b97f4a7c15ULL;
        splitKey = (splitKey ^ (splitKey >> 30)) * 0xbf58476d1ce4e5b9ULL;
        splitKey = (splitKey ^ (splitKey >> 27)) * 0x94d049bb133111ebULL;
        splitKey = splitKey ^ (splitKey >> 31);
        const bool isValidation = (splitKey & 1ULL) != 0ULL;

        int structuredMoment = 0;
        for (int i = 0; i < kQRecenteringStructuredNFeatures; ++i) {
          for (int j = i; j < kQRecenteringStructuredNFeatures; ++j) {
            const double value = structuredFeatures[i] * structuredFeatures[j];
            if (isValidation) {
              histos.fill(HIST("QRecenteringCalibration/hRegressionMomentsStructuredValidation"), centrality, structuredMoment + 0.5, value);
            } else {
              histos.fill(HIST("QRecenteringCalibration/hRegressionMomentsStructuredTrain"), centrality, structuredMoment + 0.5, value);
            }
            ++structuredMoment;
          }
        }

        for (int component = 0; component < kQRecenteringNComponents; ++component) {
          for (int i = 0; i < kQRecenteringStructuredNFeatures; ++i) {
            const int index = kQRecenteringStructuredNMatrixMoments + component * kQRecenteringStructuredNFeatures + i;
            const double value = structuredFeatures[i] * qValues[component];
            if (isValidation) {
              histos.fill(HIST("QRecenteringCalibration/hRegressionMomentsStructuredValidation"), centrality, index + 0.5, value);
            } else {
              histos.fill(HIST("QRecenteringCalibration/hRegressionMomentsStructuredTrain"), centrality, index + 0.5, value);
            }
          }

          const int q2Index = kQRecenteringStructuredNMatrixMoments + kQRecenteringNComponents * kQRecenteringStructuredNFeatures + component;
          const double q2Value = qValues[component] * qValues[component];
          if (isValidation) {
            histos.fill(HIST("QRecenteringCalibration/hRegressionMomentsStructuredValidation"), centrality, q2Index + 0.5, q2Value);
          } else {
            histos.fill(HIST("QRecenteringCalibration/hRegressionMomentsStructuredTrain"), centrality, q2Index + 0.5, q2Value);
          }
        }

        for (int i = 0; i < kQRecenteringStructuredNFeatures; ++i) {
          const double value = structuredFeatures[i];
          if (isValidation) {
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsTimeValidation"), centrality, timeFromSOR, i + 0.5, value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVxValidation"), centrality, vx, i + 0.5, value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVyValidation"), centrality, vy, i + 0.5, value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVzValidation"), centrality, vz, i + 0.5, value);
          } else {
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsTimeTrain"), centrality, timeFromSOR, i + 0.5, value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVxTrain"), centrality, vx, i + 0.5, value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVyTrain"), centrality, vy, i + 0.5, value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVzTrain"), centrality, vz, i + 0.5, value);
          }
        }

        for (int component = 0; component < kQRecenteringNComponents; ++component) {
          const int qIndex = kQRecenteringStructuredNFeatures + component;
          const int q2Index = kQRecenteringStructuredNFeatures + kQRecenteringNComponents + component;
          const double qValue = qValues[component];
          const double q2Value = qValue * qValue;

          if (isValidation) {
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsTimeValidation"), centrality, timeFromSOR, qIndex + 0.5, qValue);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVxValidation"), centrality, vx, qIndex + 0.5, qValue);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVyValidation"), centrality, vy, qIndex + 0.5, qValue);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVzValidation"), centrality, vz, qIndex + 0.5, qValue);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsTimeValidation"), centrality, timeFromSOR, q2Index + 0.5, q2Value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVxValidation"), centrality, vx, q2Index + 0.5, q2Value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVyValidation"), centrality, vy, q2Index + 0.5, q2Value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVzValidation"), centrality, vz, q2Index + 0.5, q2Value);
          } else {
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsTimeTrain"), centrality, timeFromSOR, qIndex + 0.5, qValue);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVxTrain"), centrality, vx, qIndex + 0.5, qValue);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVyTrain"), centrality, vy, qIndex + 0.5, qValue);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVzTrain"), centrality, vz, qIndex + 0.5, qValue);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsTimeTrain"), centrality, timeFromSOR, q2Index + 0.5, q2Value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVxTrain"), centrality, vx, q2Index + 0.5, q2Value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVyTrain"), centrality, vy, q2Index + 0.5, q2Value);
            histos.fill(HIST("QRecenteringCalibration/hStructuredFlatnessMomentsVzTrain"), centrality, vz, q2Index + 0.5, q2Value);
          }
        }
      }
    }

    if (calibrationStage.value == 2 && useQRecentering.value && qRecenteringProfile &&
        centrality >= qRecenteringProfile->GetXaxis()->GetXmin() &&
        centrality < qRecenteringProfile->GetXaxis()->GetXmax()) {
      const auto features = makeQRecenteringStructuredFeatures(timeFromSOR, vx, vy, vz);
      std::array<double, 4> bias{};
      bool valid = true;
      for (int component = 0; component < kQRecenteringNComponents; ++component) {
        for (std::size_t i = 0; i < features.size(); ++i) {
          const double coeff = qRecenteringProfile->GetBinContent(qRecenteringProfile->FindBin(centrality, i + 0.5, component + 0.5));
          if (!std::isfinite(coeff)) {
            valid = false;
            break;
          }
          bias[component] += coeff * features[i];
        }
        if (!std::isfinite(bias[component]))
          valid = false;
      }
      if (valid) {
        for (int i = 0; i < 4; ++i)
          qValues[i] -= bias[i];
      } else {
        LOGF(debug, "Invalid recentering map values run %d: recentering skipped without altering Q", runNumber);
      }
    }

    // Compact Q moments: 0 Qx, 1 Qy, 2 Qx^2, 3 Qy^2, 4 QxQy,
    // 5 cos(phi), 6 sin(phi), 7 cos(2phi), 8 sin(2phi),
    // 9 cos(4phi), 10 sin(4phi), 11 |Q|, 12 |Q|^2, 13 entries.
    // 1n means are not forced to zero: subtract <Q>_1n OFFLINE when studying shapes.
    const auto fillNeutronQ = [&](const std::array<double, 4>& q) {
      if (!deriveSingleNeutronResponse.value || calibrationStage.value != 2 || !useQRecentering.value)
        return;
      const auto fillSide = [&](int side, int neutronClass, double qx, double qy) {
        if (neutronClass == 0)
          return;
        const double mag = std::hypot(qx, qy);
        const double phi = std::atan2(qy, qx);
        const double vals[12] = {qx, qy, qx * qx, qy * qy, qx * qy,
                                 mag > 1.e-12 ? std::cos(phi) : 0., mag > 1.e-12 ? std::sin(phi) : 0.,
                                 mag > 1.e-12 ? std::cos(2. * phi) : 0., mag > 1.e-12 ? std::sin(2. * phi) : 0.,
                                 mag > 1.e-12 ? std::cos(4. * phi) : 0., mag > 1.e-12 ? std::sin(4. * phi) : 0., mag};
        // Keep fixed histogram names for the existing 1n/2n/3n plotting macro.
        for (int k = 0; k < 14; ++k) {
          const double value = k < 12 ? vals[k] : (k == 12 ? mag * mag : 1.);
          if (side == 0) {
            if (neutronClass == 1)
              histos.fill(HIST("NeutronQA/pMomentsRecentered1nA"), centrality, k + .5, value);
            else if (neutronClass == 2)
              histos.fill(HIST("NeutronQA/pMomentsRecentered2nA"), centrality, k + .5, value);
            else if (neutronClass == 3)
              histos.fill(HIST("NeutronQA/pMomentsRecentered3nA"), centrality, k + .5, value);
          } else {
            if (neutronClass == 1)
              histos.fill(HIST("NeutronQA/pMomentsRecentered1nC"), centrality, k + .5, value);
            else if (neutronClass == 2)
              histos.fill(HIST("NeutronQA/pMomentsRecentered2nC"), centrality, k + .5, value);
            else if (neutronClass == 3)
              histos.fill(HIST("NeutronQA/pMomentsRecentered3nC"), centrality, k + .5, value);
          }
        }
        if (side == 0) {
          if (neutronClass == 1)
            histos.fill(HIST("NeutronQA/hQxyRecentered1nA"), qx, qy);
          else if (neutronClass == 2)
            histos.fill(HIST("NeutronQA/hQxyRecentered2nA"), qx, qy);
          else if (neutronClass == 3)
            histos.fill(HIST("NeutronQA/hQxyRecentered3nA"), qx, qy);
        } else {
          if (neutronClass == 1)
            histos.fill(HIST("NeutronQA/hQxyRecentered1nC"), qx, qy);
          else if (neutronClass == 2)
            histos.fill(HIST("NeutronQA/hQxyRecentered2nC"), qx, qy);
          else if (neutronClass == 3)
            histos.fill(HIST("NeutronQA/hQxyRecentered3nC"), qx, qy);
        }
      };
      fillSide(0, nClassA, q[0], q[1]);
      fillSide(1, nClassC, q[2], q[3]);
    };
    if (deriveSingleNeutronResponse.value && calibrationStage.value == 2 && useQRecentering.value) {
      if (nClassA == 1) {
        histos.fill(HIST("NeutronQA/pRecenteredQxVsCross1nA"), crossLookupXA, crossLookupYA, qValues[0]);
        histos.fill(HIST("NeutronQA/pRecenteredQyVsCross1nA"), crossLookupXA, crossLookupYA, qValues[1]);
      }
      if (nClassC == 1) {
        histos.fill(HIST("NeutronQA/pRecenteredQxVsCross1nC"), crossLookupXC, crossLookupYC, qValues[2]);
        histos.fill(HIST("NeutronQA/pRecenteredQyVsCross1nC"), crossLookupXC, crossLookupYC, qValues[3]);
      }
      fillNeutronQ(qValues);
    }

    for (int i = 0; i < 4; ++i) {
      const double component = i + 0.5;
      histos.fill(HIST("QRecenteringQA/pQAfterVsCentrality"), centrality, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQAfterVsTime"), timeFromSOR, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQAfterVsVx"), vx, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQAfterVsVy"), vy, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQAfterVsVz"), vz, component, qValues[i]);
    }

    histos.fill(HIST("QCorrelationQA/pQACorrelationAfterVsCentrality"), centrality, 0.5, qValues[0] * qValues[2]);
    histos.fill(HIST("QCorrelationQA/pQACorrelationAfterVsCentrality"), centrality, 1.5, qValues[1] * qValues[3]);
    histos.fill(HIST("QCorrelationQA/pQACorrelationAfterVsCentrality"), centrality, 2.5, qValues[0] * qValues[3]);
    histos.fill(HIST("QCorrelationQA/pQACorrelationAfterVsCentrality"), centrality, 3.5, qValues[1] * qValues[2]);

    const double phiFinalA = std::atan2(qValues[1], qValues[0]);
    const double phiFinalC = std::atan2(qValues[3], qValues[2]);
    histos.fill(HIST("PhiQA/hPhiAfterRecenteringZNA"), phiFinalA);
    histos.fill(HIST("PhiQA/hPhiAfterRecenteringZNC"), phiFinalC);
    histos.fill(HIST("PhiQA/hPhiAfterRecenteringVsCentralityZNA"), centrality, phiFinalA);
    histos.fill(HIST("PhiQA/hPhiAfterRecenteringVsCentralityZNC"), centrality, phiFinalC);

    // Final Q vectors are from enabled tower-energy and structured71 steps.
    const double phiOutputA = std::atan2(qValues[1], qValues[0]);
    const double phiOutputC = std::atan2(qValues[3], qValues[2]);

    // Direction-sensitive TPC first-harmonic numerators for each ZNA/C side and x/y.
    // No track iteration when cfgTPCV1Correlation is false.
    // Same-side Q after the enabled corrections, per-track unit weights.
    // This is NOT a fitted v1 nor resolution-corrected; use a separate plotting macro.
    if (cfgTPCV1Correlation.value && calibrationStage.value == 2) {
      for (const auto& track : tracks) {
        if (!track.isGlobalTrack() || !std::isfinite(track.pt()) || !std::isfinite(track.eta()) || !std::isfinite(track.phi()) ||
            track.pt() <= cfgTPCResponsePtMin.value || track.pt() >= cfgTPCResponsePtMax.value ||
            std::abs(track.eta()) >= cfgTPCResponseEtaMax.value)
          continue;
        const double tx = std::cos(track.phi()), ty = std::sin(track.phi());
        const double eta = track.eta();
        histos.fill(HIST("TPCV1EtaQA/hTrackCounts"), centrality, eta);
        histos.fill(HIST("TPCV1EtaQA/pNumerators"), centrality, eta, .5, tx * qValues[0]);
        histos.fill(HIST("TPCV1EtaQA/pNumerators"), centrality, eta, 1.5, ty * qValues[1]);
        histos.fill(HIST("TPCV1EtaQA/pNumerators"), centrality, eta, 2.5, tx * qValues[2]);
        histos.fill(HIST("TPCV1EtaQA/pNumerators"), centrality, eta, 3.5, ty * qValues[3]);
      }
    }
    // ZDC resolution after the enabled corrections
    const double qMagA = std::hypot(qValues[0], qValues[1]);
    const double qMagC = std::hypot(qValues[2], qValues[3]);
    if (centrality >= 0.f && centrality < cfgCutCentralityMax.value && qMagA > 0. && qMagC > 0.) {
      const double correlation = std::cos(phiOutputA - phiOutputC);
      histos.fill(HIST("ResolutionQA/pCosPsiAPsiCVsCentrality"), centrality, correlation);
      histos.fill(HIST("ResolutionQA/pCosPsiAPsiCVsCentralityVz"), centrality, vz, correlation);
    }

    // Always describe the Q vectors actually written to the stage-2 table, regardless of enabled corrections.
    if (calibrationStage.value == 2) {
      histos.fill(HIST("final_calibration/hQxAVsCentrality"), centrality, qValues[0]);
      histos.fill(HIST("final_calibration/hQyAVsCentrality"), centrality, qValues[1]);
      histos.fill(HIST("final_calibration/hQxCVsCentrality"), centrality, qValues[2]);
      histos.fill(HIST("final_calibration/hQyCVsCentrality"), centrality, qValues[3]);
      histos.fill(HIST("final_calibration/pQxAVsCentrality"), centrality, qValues[0]);
      histos.fill(HIST("final_calibration/pQyAVsCentrality"), centrality, qValues[1]);
      histos.fill(HIST("final_calibration/pQxCVsCentrality"), centrality, qValues[2]);
      histos.fill(HIST("final_calibration/pQyCVsCentrality"), centrality, qValues[3]);
      if (qMagA > 0.) {
        histos.fill(HIST("final_calibration/hPhiAVsCentrality"), centrality, phiOutputA);
      }
      if (qMagC > 0.) {
        histos.fill(HIST("final_calibration/hPhiCVsCentrality"), centrality, phiOutputC);
      }
      if (qMagA > 0. && qMagC > 0.) {
        const double finalCorrelation = std::cos(phiOutputA - phiOutputC);
        histos.fill(HIST("final_calibration/pCosPsiAPsiCVsCentrality"), centrality, finalCorrelation);
        histos.fill(HIST("final_calibration/pCosPsiAPsiCVsCentralityVz"), centrality, vz, finalCorrelation);
      }
      histos.fill(HIST("final_calibration/pQxAQxCVsCentrality"), centrality, qValues[0] * qValues[2]);
      histos.fill(HIST("final_calibration/pQyAQyCVsCentrality"), centrality, qValues[1] * qValues[3]);
      histos.fill(HIST("final_calibration/pQxAQyCVsCentrality"), centrality, qValues[0] * qValues[3]);
      histos.fill(HIST("final_calibration/pQyAQxCVsCentrality"), centrality, qValues[1] * qValues[2]);
      histos.fill(HIST("final_calibration/pQxAVsCentralityVz"), centrality, vz, qValues[0]);
      histos.fill(HIST("final_calibration/pQyAVsCentralityVz"), centrality, vz, qValues[1]);
      histos.fill(HIST("final_calibration/pQxCVsCentralityVz"), centrality, vz, qValues[2]);
      histos.fill(HIST("final_calibration/pQyCVsCentralityVz"), centrality, vz, qValues[3]);
      // Ten basic Q moments retained for final calibration QA (no twist applied).
      const std::array<double, 10> finalMoments = {
        qValues[0], qValues[1], qValues[0] * qValues[0],
        qValues[1] * qValues[1], qValues[0] * qValues[1],
        qValues[2], qValues[3], qValues[2] * qValues[2],
        qValues[3] * qValues[3], qValues[2] * qValues[3]};
      for (int i = 0; i < 10; ++i) {
        histos.fill(HIST("final_calibration/pQMomentsVsCentrality"), centrality, i + 0.5, finalMoments[i]);
      }
      for (int n = 1; n <= 4; ++n) {
        if (qMagA > 0.) {
          histos.fill(HIST("final_calibration/pPhiHarmonicsAVsCentrality"), centrality, 2 * n - 1.5, std::cos(n * phiOutputA));
          histos.fill(HIST("final_calibration/pPhiHarmonicsAVsCentrality"), centrality, 2 * n - 0.5, std::sin(n * phiOutputA));
        }
        if (qMagC > 0.) {
          histos.fill(HIST("final_calibration/pPhiHarmonicsCVsCentrality"), centrality, 2 * n - 1.5, std::cos(n * phiOutputC));
          histos.fill(HIST("final_calibration/pPhiHarmonicsCVsCentrality"), centrality, 2 * n - 0.5, std::sin(n * phiOutputC));
        }
      }
      histos.fill(HIST("final_calibration/pQMagAVsCentrality"), centrality, qMagA);
      histos.fill(HIST("final_calibration/pQMagCVsCentrality"), centrality, qMagC);
    }

    histos.fill(HIST("hEvtSelInfo"), 7.5);
    fillTable(true, static_cast<float>(qValues[0]), static_cast<float>(qValues[2]), static_cast<float>(qValues[1]), static_cast<float>(qValues[3]));
    histos.fill(HIST("hEvtSelInfo"), 8.5);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<zdc2stagecalibration>(context)};
}
