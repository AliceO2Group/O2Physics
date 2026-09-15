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

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
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
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Centrality cut Max"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Centrality cut Min"};
  Configurable<int> calibrationStage{"calibrationStage", 1, "1: derive detector calibration QA, 2: apply detector calibration and produce corrected Q"};
  Configurable<float> cfgTimeSliceMinutes{"cfgTimeSliceMinutes", 30.0f, "Calibration time-slice width in minutes"};
  Configurable<float> cfgMaxRunHours{"cfgMaxRunHours", 48.0f, "Maximum run duration shown in QA"};
  Configurable<bool> followpub{"followpub", true, "flag to use alphaZDC"};
  Configurable<bool> usecfactor{"usecfactor", false, "use c factor"};
  Configurable<bool> useGainCallib{"useGainCallib", false, "use stage-1 tower gain calibration"};
  Configurable<std::string> confGainPath{"confGainPath", "", "CCDB path to the run-wise stage-1 gain map"};
  Configurable<bool> useSpatialCalib{"useSpatialCalib", false, "use ZDC spatial calibration"};
  Configurable<std::string> confSpatialPath{"confSpatialPath", "", "CCDB path to ZDC spatial calibration"};
  Configurable<bool> deriveSpatialCalib{"deriveSpatialCalib", false, "store spatial regression moments after gain correction"};
  Configurable<bool> deriveQRecentering{"deriveQRecentering", false, "store Q-vector recentering regression moments after gain and spatial correction"};
  Configurable<bool> useQRecentering{"useQRecentering", false, "apply Q-vector recentering calibration"};
  Configurable<std::string> confQRecenteringPath{"confQRecenteringPath", "", "CCDB path to the run-wise Q-vector recentering coefficients"};

  Configurable<int> cfgSpatialCalibBins{
    "cfgSpatialCalibBins", 30,
    "number of bins per cross coordinate for spatial calibration"};
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
  TH3D* spatialprofile = nullptr;
  TH2D* qRecenteringProfile = nullptr;

  static constexpr int kQRecenteringNFeatures = 21;
  static constexpr int kQRecenteringNMatrixMoments = kQRecenteringNFeatures * (kQRecenteringNFeatures + 1) / 2;
  static constexpr int kQRecenteringNComponents = 4;
  static constexpr int kQRecenteringNMoments = kQRecenteringNMatrixMoments + kQRecenteringNComponents * kQRecenteringNFeatures;

  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    const int nTimeBins = static_cast<int>(std::ceil(cfgMaxRunHours.value * 60.f / cfgTimeSliceMinutes.value));
    AxisSpec timeAxis = {nTimeBins, 0.0, cfgMaxRunHours.value, "time from SOR (h)"};
    AxisSpec towerAxis = {4, 0.0, 4.0, "tower"};
    AxisSpec crossAxis = {120, -1.0, 1.0, "cross asymmetry"};
    AxisSpec momentAxis = {15, 0.0, 15.0, "regression moment"};
    AxisSpec spatialMomentAxis = {6, 0.0, 6.0, "normalized pair regression moment"};
    AxisSpec spatialCalibAxis = {cfgSpatialCalibBins.value, -1.0, 1.0, "cross asymmetry"};
    AxisSpec phiAxis = {72, -3.14159265358979323846, 3.14159265358979323846, "#phi"};
    AxisSpec crossCoordAxis = {2, 0.0, 2.0, "coordinate"};
    AxisSpec ratioAxis = {240, 0.0, 2.4, "#Sigma tower/common"};
    AxisSpec qRecenteringMomentAxis = {kQRecenteringNMoments, 0.0, static_cast<double>(kQRecenteringNMoments), "regression moment"};
    AxisSpec centralityAxis = {80, 0.0, 80.0, "centrality (%)"};
    AxisSpec vertexXYAxis = {100, -0.5, 0.5, "vertex x/y (cm)"};
    AxisSpec vertexZAxis = {100, -10.0, 10.0, "vertex z (cm)"};
    AxisSpec qComponentAxis = {4, 0.0, 4.0, "Q component"};

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
    histos.add("GainQA/hSumOverCommonVsTimeZNA", "ZNA tower sum/common vs time;time from SOR (h);#Sigma T_{i}/C", kTH2F, {timeAxis, ratioAxis});
    histos.add("GainQA/hSumOverCommonVsTimeZNC", "ZNC tower sum/common vs time;time from SOR (h);#Sigma T_{i}/C", kTH2F, {timeAxis, ratioAxis});
    histos.add("GainQA/pTowerOverCommonVsTimeZNA", "ZNA individual tower/common vs time;time from SOR (h);tower;<T_{i}/C>", kTProfile2D, {timeAxis, towerAxis});
    histos.add("GainQA/pTowerOverCommonVsTimeZNC", "ZNC individual tower/common vs time;time from SOR (h);tower;<T_{i}/C>", kTProfile2D, {timeAxis, towerAxis});
    histos.add("SpatialQA/pCrossVsTimeZNA", "ZNA cross asymmetry vs time;time from SOR (h);0:X cross, 1:Y cross;<cross>", kTProfile2D, {timeAxis, crossCoordAxis});
    histos.add("SpatialQA/pCrossVsTimeZNC", "ZNC cross asymmetry vs time;time from SOR (h);0:X cross, 1:Y cross;<cross>", kTProfile2D, {timeAxis, crossCoordAxis});
    histos.add("SpatialQA/pResponseVsCrossXYZNA", "ZNA response map;X cross;Y cross;<#Sigma T_{i}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pResponseVsCrossXYZNC", "ZNC response map;X cross;Y cross;<#Sigma T_{i}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pCentroidVsTimeZNA", "ZNA centroid vs time;time from SOR (h);0:x, 1:y;<centroid (cm)>", kTProfile2D, {timeAxis, crossCoordAxis});
    histos.add("SpatialQA/pCentroidVsTimeZNC", "ZNC centroid vs time;time from SOR (h);0:x, 1:y;<centroid (cm)>", kTProfile2D, {timeAxis, crossCoordAxis});
    histos.add("SpatialQA/hCrossXYZNA", "ZNA cross occupancy;X cross;Y cross", kTH2F, {crossAxis, crossAxis});
    histos.add("SpatialQA/hCrossXYZNC", "ZNC cross occupancy;X cross;Y cross", kTH2F, {crossAxis, crossAxis});
    histos.add("SpatialQA/pTower1OverCommonVsCrossXYZNA", "ZNA tower 1 response;X cross;Y cross;<T_{1}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pTower2OverCommonVsCrossXYZNA", "ZNA tower 2 response;X cross;Y cross;<T_{2}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pTower3OverCommonVsCrossXYZNA", "ZNA tower 3 response;X cross;Y cross;<T_{3}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pTower4OverCommonVsCrossXYZNA", "ZNA tower 4 response;X cross;Y cross;<T_{4}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pTower1OverCommonVsCrossXYZNC", "ZNC tower 1 response;X cross;Y cross;<T_{1}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pTower2OverCommonVsCrossXYZNC", "ZNC tower 2 response;X cross;Y cross;<T_{2}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pTower3OverCommonVsCrossXYZNC", "ZNC tower 3 response;X cross;Y cross;<T_{3}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("SpatialQA/pTower4OverCommonVsCrossXYZNC", "ZNC tower 4 response;X cross;Y cross;<T_{4}/C>", kTProfile2D, {crossAxis, crossAxis});
    histos.add("PhiQA/hPhiRawZNA", "ZNA #phi raw;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiRawZNC", "ZNC #phi raw;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterGainZNA", "ZNA #phi after gain calibration;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterGainZNC", "ZNC #phi after gain calibration;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterSpatialZNA", "ZNA #phi after spatial calibration;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterSpatialZNC", "ZNC #phi after spatial calibration;#phi;events", kTH1F, {phiAxis});

    // Stage-1 correction input. For each time bin the first 10 y bins contain sum(T_i*T_j), the next 4 contain sum(C*T_i), and the last contains the event count.
    histos.add("GainCalibration/hGainMomentsZNA", "ZNA linear-regression moments;time from SOR (h);moment index", kTH2D, {timeAxis, momentAxis});
    histos.add("GainCalibration/hGainMomentsZNC", "ZNC linear-regression moments;time from SOR (h);moment index", kTH2D, {timeAxis, momentAxis});
    // Normalized two-pair spatial regression moments. For each (Xcross,Ycross) cell:
    // u=(T1+T4)/C, v=(T2+T3)/C, and moments are u^2, u*v, v^2, u, v, N.
    histos.add("SpatialCalibration/hSpatialMomentsZNA", "ZNA normalized pair spatial regression moments;X cross;Y cross;moment index", kTH3D, {spatialCalibAxis, spatialCalibAxis, spatialMomentAxis});
    histos.add("SpatialCalibration/hSpatialMomentsZNC", "ZNC normalized pair spatial regression moments;X cross;Y cross;moment index", kTH3D, {spatialCalibAxis, spatialCalibAxis, spatialMomentAxis});

    // Q-recentering regression moments. Feature order:
    // 0:1, 1:C, 2:T, 3:Vx, 4:Vy, 5:Vz, 6:C2, 7:T2, 8:Vx2, 9:Vy2, 10:Vz2,
    // 11:C*T, 12:C*Vx, 13:C*Vy, 14:C*Vz, 15:T*Vx, 16:T*Vy, 17:T*Vz,
    // 18:Vx*Vy, 19:Vx*Vz, 20:Vy*Vz.
    // Moments 0..230 are the upper triangle of sum(X_i X_j).
    // Then 4 blocks of 21 sum(X_i Q): QxA, QyA, QxC, QyC.
    histos.add("QRecenteringCalibration/hRegressionMoments", "Q recentering regression moments;moment index;sum", kTH1D, {qRecenteringMomentAxis});

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
    histos.add("PhiQA/hPhiAfterRecenteringZNA", "ZNA #phi after Q recentering;#phi;events", kTH1F, {phiAxis});
    histos.add("PhiQA/hPhiAfterRecenteringZNC", "ZNC #phi after Q recentering;#phi;events", kTH1F, {phiAxis});

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
    bcSOR = static_cast<int64_t>(runInfo.orbitSOR) * static_cast<int64_t>(o2::constants::lhc::LHCMaxBunches);
    gainprofile = nullptr;
    spatialprofile = nullptr;
    qRecenteringProfile = nullptr;
    if (calibrationStage.value == 2 && useGainCallib.value && !confGainPath.value.empty()) {
      gainprofile = ccdb->getForTimeStamp<TH2D>(confGainPath.value, timestamp);
      if (!gainprofile) {
        LOGF(warn, "No ZDC stage-1 gain map found for run %d at timestamp %llu", runNumber, static_cast<unsigned long long>(timestamp));
      }
    }
    if (calibrationStage.value == 2 && useSpatialCalib.value && !confSpatialPath.value.empty()) {
      spatialprofile = ccdb->getForTimeStamp<TH3D>(confSpatialPath.value, timestamp);

      if (!spatialprofile) {
        LOGF(warn, "No ZDC spatial calibration found for run %d at timestamp %llu", runNumber, static_cast<unsigned long long>(timestamp));
      }
    }
    if (calibrationStage.value == 2 && useQRecentering.value && !confQRecenteringPath.value.empty()) {
      qRecenteringProfile = ccdb->getForTimeStamp<TH2D>(confQRecenteringPath.value, timestamp);
      if (!qRecenteringProfile) {
        LOGF(warn, "No ZDC Q-recentering calibration found for run %d at timestamp %llu", runNumber, static_cast<unsigned long long>(timestamp));
      }
    }
  }

  float getTimeFromSOR(uint64_t globalBC) const
  {
    if (bcSOR < 0 || globalBC < static_cast<uint64_t>(bcSOR)) {
      return -1.f;
    }
    const double deltaBC = static_cast<double>(globalBC - static_cast<uint64_t>(bcSOR));
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

  std::array<double, kQRecenteringNFeatures> makeQRecenteringFeatures(float centrality, float timeFromSOR, float vx, float vy, float vz) const
  {
    // Fixed scaling keeps the polynomial basis numerically well behaved and must be
    // used identically when solving/applying the CCDB coefficients.
    const double c = (static_cast<double>(centrality) - 40.0) / 40.0;
    const double runHours = (eorTimestamp > sorTimestamp) ? static_cast<double>(eorTimestamp - sorTimestamp) * 1.e-3 / 3600.0 : static_cast<double>(cfgMaxRunHours.value);
    const double t = (runHours > 0.0) ? (2.0 * static_cast<double>(timeFromSOR) / runHours - 1.0) : 0.0;
    const double x = static_cast<double>(vx) / 0.1;
    const double y = static_cast<double>(vy) / 0.1;
    const double z = static_cast<double>(vz) / 10.0;

    return {1.0, c, t, x, y, z,
            c * c, t * t, x * x, y * y, z * z,
            c * t, c * x, c * y, c * z,
            t * x, t * y, t * z,
            x * y, x * z, y * z};
  }

  using MyCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::Mults, o2::aod::FT0sCorrected, o2::aod::CentFT0Cs>;

  void process(MyCollisions::iterator const& collision, o2::aod::FT0s const& /*ft0s*/, o2::aod::FV0As const& /*fv0s*/, BCsRun3 const& /*bcs*/, o2::aod::Zdcs const& /*zdcs*/)
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
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 1.5);

    auto zdc = bc.zdc();
    auto znaEnergy = zdc.energySectorZNA();
    auto zncEnergy = zdc.energySectorZNC();
    const float znaEnergycommon = zdc.energyCommonZNA();
    const float zncEnergycommon = zdc.energyCommonZNC();

    if (!std::isfinite(znaEnergycommon) || !std::isfinite(zncEnergycommon) || znaEnergycommon <= 0.f || zncEnergycommon <= 0.f) {
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 2.5);

    bool validTowers = true;
    for (int i = 0; i < 4; ++i) {
      if (!std::isfinite(znaEnergy[i]) || !std::isfinite(zncEnergy[i]) || znaEnergy[i] <= 0.f || zncEnergy[i] <= 0.f) {
        validTowers = false;
      }
    }
    if (!validTowers) {
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 3.5);

    if (rctCut.requireRCTFlagChecker && !rctChecker(collision)) {
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 4.5);

    if (!(collision.sel8() && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && std::isfinite(centrality) && centrality > cfgCutCentralityMin && centrality < cfgCutCentralityMax && std::abs(vz) < cfgCutVertex && collision.has_foundFT0())) {
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 5.5);

    if (timeFromSOR < 0.f || timeFromSOR >= cfgMaxRunHours.value) {
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }

    constexpr double x[4] = {-1.75, 1.75, -1.75, 1.75};
    constexpr double y[4] = {-1.75, -1.75, 1.75, 1.75};

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
    for (int i = 0; i < 4; ++i) {
      znaCorr[i] = getGain(timeFromSOR, i, gainOK) * znaEnergy[i];
      zncCorr[i] = getGain(timeFromSOR, i + 4, gainOK) * zncEnergy[i];
    }
    if (!gainOK) {
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 6.5);

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

    const double crossLookupXA = (znaCorr[3] - znaCorr[0]) / (znaCorr[3] + znaCorr[0]);
    const double crossLookupYA = (znaCorr[2] - znaCorr[1]) / (znaCorr[2] + znaCorr[1]);
    const double crossLookupXC = (zncCorr[3] - zncCorr[0]) / (zncCorr[3] + zncCorr[0]);
    const double crossLookupYC = (zncCorr[2] - zncCorr[1]) / (zncCorr[2] + zncCorr[1]);
    if (calibrationStage.value == 2 && useGainCallib.value && deriveSpatialCalib.value && !useSpatialCalib.value) {

      const double uA = (znaCorr[0] + znaCorr[3]) / static_cast<double>(znaEnergycommon);
      const double vA = (znaCorr[1] + znaCorr[2]) / static_cast<double>(znaEnergycommon);
      const double uC = (zncCorr[0] + zncCorr[3]) / static_cast<double>(zncEnergycommon);
      const double vC = (zncCorr[1] + zncCorr[2]) / static_cast<double>(zncEnergycommon);

      if (std::isfinite(uA) && std::isfinite(vA) && std::isfinite(uC) && std::isfinite(vC)) {
        const std::array<double, 6> spatialMomentsA = {
          uA * uA,
          uA * vA,
          vA * vA,
          uA,
          vA,
          1.0};
        const std::array<double, 6> spatialMomentsC = {
          uC * uC,
          uC * vC,
          vC * vC,
          uC,
          vC,
          1.0};

        for (int i = 0; i < 6; ++i) {
          histos.fill(HIST("SpatialCalibration/hSpatialMomentsZNA"), crossLookupXA, crossLookupYA, i + 0.5, spatialMomentsA[i]);
          histos.fill(HIST("SpatialCalibration/hSpatialMomentsZNC"), crossLookupXC, crossLookupYC, i + 0.5, spatialMomentsC[i]);
        }
      }
    }

    bool spatialOK = true;
    if (calibrationStage.value == 2 && useSpatialCalib.value) {
      if (!spatialprofile) {
        spatialOK = false;
      } else {
        for (int i = 0; i < 4; ++i) {
          const double corrA = spatialprofile->GetBinContent(spatialprofile->FindBin(crossLookupXA, crossLookupYA, i + 0.5));
          const double corrC = spatialprofile->GetBinContent(spatialprofile->FindBin(crossLookupXC, crossLookupYC, i + 4.5));
          if (!std::isfinite(corrA) || corrA <= 0.0 || !std::isfinite(corrC) || corrC <= 0.0) {
            spatialOK = false;
            break;
          }
          znaCorr[i] *= corrA;
          zncCorr[i] *= corrC;
        }
      }
    }
    if (!spatialOK) {
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
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
    const double phiSpatialA = std::atan2(centroidYA, centroidXA);
    const double phiSpatialC = std::atan2(centroidYC, centroidXC);
    histos.fill(HIST("PhiQA/hPhiAfterSpatialZNA"), phiSpatialA);
    histos.fill(HIST("PhiQA/hPhiAfterSpatialZNC"), phiSpatialC);

    histos.fill(HIST("GainQA/hSumOverCommonVsTimeZNA"), timeFromSOR, ratioA);
    histos.fill(HIST("GainQA/hSumOverCommonVsTimeZNC"), timeFromSOR, ratioC);
    for (int i = 0; i < 4; ++i) {
      histos.fill(HIST("GainQA/pTowerOverCommonVsTimeZNA"), timeFromSOR, i + 0.5, znaCorr[i] / znaEnergycommon);
      histos.fill(HIST("GainQA/pTowerOverCommonVsTimeZNC"), timeFromSOR, i + 0.5, zncCorr[i] / zncEnergycommon);
    }
    histos.fill(HIST("SpatialQA/pCrossVsTimeZNA"), timeFromSOR, 0.5, crossXA);
    histos.fill(HIST("SpatialQA/pCrossVsTimeZNA"), timeFromSOR, 1.5, crossYA);
    histos.fill(HIST("SpatialQA/pCrossVsTimeZNC"), timeFromSOR, 0.5, crossXC);
    histos.fill(HIST("SpatialQA/pCrossVsTimeZNC"), timeFromSOR, 1.5, crossYC);
    histos.fill(HIST("SpatialQA/pResponseVsCrossXYZNA"), crossLookupXA, crossLookupYA, ratioA);
    histos.fill(HIST("SpatialQA/pResponseVsCrossXYZNC"), crossLookupXC, crossLookupYC, ratioC);
    histos.fill(HIST("SpatialQA/pCentroidVsTimeZNA"), timeFromSOR, 0.5, centroidXA);
    histos.fill(HIST("SpatialQA/pCentroidVsTimeZNA"), timeFromSOR, 1.5, centroidYA);
    histos.fill(HIST("SpatialQA/pCentroidVsTimeZNC"), timeFromSOR, 0.5, centroidXC);
    histos.fill(HIST("SpatialQA/pCentroidVsTimeZNC"), timeFromSOR, 1.5, centroidYC);
    histos.fill(HIST("SpatialQA/hCrossXYZNA"), crossLookupXA, crossLookupYA);
    histos.fill(HIST("SpatialQA/hCrossXYZNC"), crossLookupXC, crossLookupYC);

    histos.fill(HIST("SpatialQA/pTower1OverCommonVsCrossXYZNA"), crossLookupXA, crossLookupYA, znaCorr[0] / znaEnergycommon);
    histos.fill(HIST("SpatialQA/pTower2OverCommonVsCrossXYZNA"), crossLookupXA, crossLookupYA, znaCorr[1] / znaEnergycommon);
    histos.fill(HIST("SpatialQA/pTower3OverCommonVsCrossXYZNA"), crossLookupXA, crossLookupYA, znaCorr[2] / znaEnergycommon);
    histos.fill(HIST("SpatialQA/pTower4OverCommonVsCrossXYZNA"), crossLookupXA, crossLookupYA, znaCorr[3] / znaEnergycommon);

    histos.fill(HIST("SpatialQA/pTower1OverCommonVsCrossXYZNC"), crossLookupXC, crossLookupYC, zncCorr[0] / zncEnergycommon);
    histos.fill(HIST("SpatialQA/pTower2OverCommonVsCrossXYZNC"), crossLookupXC, crossLookupYC, zncCorr[1] / zncEnergycommon);
    histos.fill(HIST("SpatialQA/pTower3OverCommonVsCrossXYZNC"), crossLookupXC, crossLookupYC, zncCorr[2] / zncEnergycommon);
    histos.fill(HIST("SpatialQA/pTower4OverCommonVsCrossXYZNC"), crossLookupXC, crossLookupYC, zncCorr[3] / zncEnergycommon);
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
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
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
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }

    const auto qFeatures = makeQRecenteringFeatures(centrality, timeFromSOR, vx, vy, vz);
    std::array<double, 4> qValues = {qxZDCA, qyZDCA, qxZDCC, qyZDCC};

    for (int i = 0; i < 4; ++i) {
      const double component = i + 0.5;
      histos.fill(HIST("QRecenteringQA/pQBeforeVsCentrality"), centrality, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQBeforeVsTime"), timeFromSOR, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQBeforeVsVx"), vx, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQBeforeVsVy"), vy, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQBeforeVsVz"), vz, component, qValues[i]);
    }

    if (calibrationStage.value == 2 && useGainCallib.value && useSpatialCalib.value && deriveQRecentering.value && !useQRecentering.value) {
      int moment = 0;
      for (int i = 0; i < kQRecenteringNFeatures; ++i) {
        for (int j = i; j < kQRecenteringNFeatures; ++j) {
          histos.fill(HIST("QRecenteringCalibration/hRegressionMoments"), moment + 0.5, qFeatures[i] * qFeatures[j]);
          ++moment;
        }
      }
      for (int component = 0; component < kQRecenteringNComponents; ++component) {
        for (int i = 0; i < kQRecenteringNFeatures; ++i) {
          const int index = kQRecenteringNMatrixMoments + component * kQRecenteringNFeatures + i;
          histos.fill(HIST("QRecenteringCalibration/hRegressionMoments"), index + 0.5, qFeatures[i] * qValues[component]);
        }
      }
    }

    bool qRecenteringOK = true;
    if (calibrationStage.value == 2 && useQRecentering.value) {
      if (!qRecenteringProfile) {
        qRecenteringOK = false;
      } else {
        for (int component = 0; component < kQRecenteringNComponents; ++component) {
          double predictedBias = 0.0;
          for (int i = 0; i < kQRecenteringNFeatures; ++i) {
            const double coefficient = qRecenteringProfile->GetBinContent(qRecenteringProfile->FindBin(i + 0.5, component + 0.5));
            if (!std::isfinite(coefficient)) {
              qRecenteringOK = false;
              break;
            }
            predictedBias += coefficient * qFeatures[i];
          }
          if (!qRecenteringOK || !std::isfinite(predictedBias)) {
            qRecenteringOK = false;
            break;
          }
          qValues[component] -= predictedBias;
        }
      }
    }
    if (!qRecenteringOK) {
      fillTable(false, 0.f, 0.f, 0.f, 0.f);
      return;
    }

    for (int i = 0; i < 4; ++i) {
      const double component = i + 0.5;
      histos.fill(HIST("QRecenteringQA/pQAfterVsCentrality"), centrality, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQAfterVsTime"), timeFromSOR, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQAfterVsVx"), vx, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQAfterVsVy"), vy, component, qValues[i]);
      histos.fill(HIST("QRecenteringQA/pQAfterVsVz"), vz, component, qValues[i]);
    }

    const double phiRecenteringA = std::atan2(qValues[1], qValues[0]);
    const double phiRecenteringC = std::atan2(qValues[3], qValues[2]);
    histos.fill(HIST("PhiQA/hPhiAfterRecenteringZNA"), phiRecenteringA);
    histos.fill(HIST("PhiQA/hPhiAfterRecenteringZNC"), phiRecenteringC);

    histos.fill(HIST("hEvtSelInfo"), 7.5);
    fillTable(true, static_cast<float>(qValues[0]), static_cast<float>(qValues[2]), static_cast<float>(qValues[1]), static_cast<float>(qValues[3]));
    histos.fill(HIST("hEvtSelInfo"), 8.5);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<zdc2stagecalibration>(cfgc)};
}
