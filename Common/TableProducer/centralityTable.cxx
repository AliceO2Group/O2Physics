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
/// \file centralityTable.cxx
/// \brief Task to produce the centrality tables associated to each of the required centrality estimators
///
/// \author ALICE
//

#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include <CCDB/BasicCCDBManager.h>
#include <TH1F.h>
#include <TFormula.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "MetadataHelper.h"
#include "TableHelper.h"
#include "TList.h"

using namespace o2;
using namespace o2::framework;

MetadataHelper metadataInfo; // Metadata helper

static constexpr int kCentRun2V0Ms = 0;
static constexpr int kCentRun2V0As = 1;
static constexpr int kCentRun2SPDTrks = 2;
static constexpr int kCentRun2SPDClss = 3;
static constexpr int kCentRun2CL0s = 4;
static constexpr int kCentRun2CL1s = 5;
static constexpr int kCentFV0As = 6;
static constexpr int kCentFT0Ms = 7;
static constexpr int kCentFT0As = 8;
static constexpr int kCentFT0Cs = 9;
static constexpr int kCentFT0CVariant1s = 10;
static constexpr int kCentFDDMs = 11;
static constexpr int kCentNTPVs = 12;
static constexpr int kCentNGlobals = 13;
static constexpr int kCentMFTs = 14;
static constexpr int NTables = 15;
static constexpr int NParameters = 1;
static const std::vector<std::string> tableNames{"CentRun2V0Ms",
                                                 "CentRun2V0As",
                                                 "CentRun2SPDTrks",
                                                 "CentRun2SPDClss",
                                                 "CentRun2CL0s",
                                                 "CentRun2CL1s",
                                                 "CentFV0As",
                                                 "CentFT0Ms",
                                                 "CentFT0As",
                                                 "CentFT0Cs",
                                                 "CentFT0CVariant1s",
                                                 "CentFDDMs",
                                                 "CentNTPVs",
                                                 "CentNGlobals",
                                                 "CentMFTs"};
static const std::vector<std::string> parameterNames{"Enable"};
static const int defaultParameters[NTables][NParameters]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};

struct CentralityTable {
  Produces<aod::CentRun2V0Ms> centRun2V0M;
  Produces<aod::CentRun2V0As> centRun2V0A;
  Produces<aod::CentRun2SPDTrks> centRun2SPDTracklets;
  Produces<aod::CentRun2SPDClss> centRun2SPDClusters;
  Produces<aod::CentRun2CL0s> centRun2CL0;
  Produces<aod::CentRun2CL1s> centRun2CL1;
  Produces<aod::CentFV0As> centFV0A;
  Produces<aod::CentFT0Ms> centFT0M;
  Produces<aod::CentFT0As> centFT0A;
  Produces<aod::CentFT0Cs> centFT0C;
  Produces<aod::CentFT0CVariant1s> centFT0CVariant1;
  Produces<aod::CentFDDMs> centFDDM;
  Produces<aod::CentNTPVs> centNTPV;
  Produces<aod::CentNGlobals> centNGlobals;
  Produces<aod::CentMFTs> centMFTs;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], NTables, NParameters, tableNames, parameterNames},
                                                "Produce tables depending on needs. Values different than -1 override the automatic setup: the corresponding table can be set off (0) or on (1)"};
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
    Configurable<std::string> ccdbPath{"ccdbPath", "Centrality/Estimators", "The CCDB path for centrality/multiplicity information"};
    Configurable<std::string> genName{"genName", "", "Genearator name: HIJING, PYTHIA8, ... Default: \"\""};
    Configurable<bool> doNotCrashOnNull{"doNotCrashOnNull", false, {"Option to not crash on null and instead fill required tables with dummy info"}};
    Configurable<std::string> reconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};
  } ccdbConfig;

  Configurable<bool> embedINELgtZEROselection{"embedINELgtZEROselection", false, {"Option to do percentile 100.5 if not INELgtZERO"}};
  Configurable<bool> produceHistograms{"produceHistograms", false, {"Option to produce debug histograms"}};
  ConfigurableAxis binsPercentile{"binsPercentile", {VARIABLE_WIDTH, 0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0}, "Binning of the percentile axis"};

  int mRunNumber;
  struct TagRun2V0MCalibration {
    bool mCalibrationStored = false;
    TFormula* mMCScale = nullptr;
    float mMCScalePars[6] = {0.0};
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhVtxAmpCorrV0C = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2V0MInfo;
  struct TagRun2V0ACalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2V0AInfo;
  struct TagRun2SPDTrackletsCalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2SPDTksInfo;
  struct TagRun2SPDClustersCalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorrCL0 = nullptr;
    TH1* mhVtxAmpCorrCL1 = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2SPDClsInfo;
  struct TagRun2CL0Calibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2CL0Info;
  struct TagRun2CL1Calibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2CL1Info;
  struct CalibrationInfo {
    std::string name = "";
    bool mCalibrationStored = false;
    TH1* mhMultSelCalib = nullptr;
    float mMCScalePars[6] = {0.0};
    TFormula* mMCScale = nullptr;
    explicit CalibrationInfo(std::string name)
      : name(name),
        mCalibrationStored(false),
        mhMultSelCalib(nullptr),
        mMCScalePars{0.0},
        mMCScale(nullptr)
    {
    }
    bool isSane(bool fatalize = false)
    {
      if (!mhMultSelCalib) {
        return true;
      }
      for (int i = 1; i < mhMultSelCalib->GetNbinsX() + 1; i++) {
        if (mhMultSelCalib->GetXaxis()->GetBinLowEdge(i) > mhMultSelCalib->GetXaxis()->GetBinUpEdge(i)) {
          if (fatalize) {
            LOG(fatal) << "Centrality calibration table " << name << " has bins with low edge > up edge";
          }
          LOG(warning) << "Centrality calibration table " << name << " has bins with low edge > up edge";
          return false;
        }
      }
      return true;
    }
  };
  CalibrationInfo fv0aInfo = CalibrationInfo("FV0");
  CalibrationInfo ft0mInfo = CalibrationInfo("FT0");
  CalibrationInfo ft0aInfo = CalibrationInfo("FT0A");
  CalibrationInfo ft0cInfo = CalibrationInfo("FT0C");
  CalibrationInfo ft0cVariant1Info = CalibrationInfo("FT0Cvar1");
  CalibrationInfo fddmInfo = CalibrationInfo("FDD");
  CalibrationInfo ntpvInfo = CalibrationInfo("NTracksPV");
  CalibrationInfo nGlobalInfo = CalibrationInfo("NGlobal");
  CalibrationInfo mftInfo = CalibrationInfo("MFT");
  std::vector<int> mEnabledTables; // Vector of enabled tables
  std::array<bool, NTables> isTableEnabled;

  // Debug output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TList> listCalib{"calib-list", OutputObjHandlingPolicy::QAObject};

  void init(InitContext& context)
  {
    LOG(info) << "Initializing centrality table producer";
    if (doprocessRun3FT0 == true) {
      LOG(fatal) << "FT0 only mode is automatically enabled in Run3 mode. Please disable it and enable processRun3.";
    }
    if (doprocessRun2 == false && doprocessRun3 == false && doprocessRun3Complete == false) {
      LOGF(fatal, "Neither processRun2 nor processRun3 nor processRun3Complete enabled. Please choose one.");
    }
    if (doprocessRun2 == true && doprocessRun3 == true) {
      LOGF(fatal, "Cannot enable processRun2 and processRun3 at the same time. Please choose one.");
    }

    /* Checking the tables which are requested in the workflow and enabling them */
    for (int i = 0; i < NTables; i++) {
      int f = enabledTables->get(tableNames[i].c_str(), "Enable");
      enableFlagIfTableRequired(context, tableNames[i], f);
      if (f == 1) {
        if (tableNames[i].find("Run2") != std::string::npos) {
          if (doprocessRun3) {
            LOG(fatal) << "Cannot enable Run2 table `" << tableNames[i] << "` while running in Run3 mode. Please check and disable them.";
          }
        } else {
          if (doprocessRun2) {
            LOG(fatal) << "Cannot enable Run3 table `" << tableNames[i] << "` while running in Run2 mode. Please check and disable them.";
          }
        }
        isTableEnabled[i] = true;
        mEnabledTables.push_back(i);
      }
    }

    if (mEnabledTables.size() == 0) {
      LOGF(fatal, "No table enabled. Please enable at least one table.");
    }
    std::sort(mEnabledTables.begin(), mEnabledTables.end());

    // Check if FT0 is the only centrality needed
    if (mEnabledTables.size() == 1 && isTableEnabled[kCentFT0Ms] == true) {
      LOG(info) << "FT0 only mode is enabled";
      doprocessRun3FT0.value = true;
      doprocessRun3.value = false;
    }

    ccdb->setURL(ccdbConfig.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    mRunNumber = 0;
    if (!produceHistograms.value) {
      return;
    }

    histos.add("FT0M/percentile", "FT0M percentile.", HistType::kTH1D, {{binsPercentile, "FT0M percentile"}});
    histos.add("FT0M/percentilevsPV", "percentile vs PV mult.", HistType::kTH2D, {{binsPercentile, "FT0M percentile"}, {100, 0, 100, "PV mult."}});
    histos.add("FT0M/MultvsPV", "FT0M mult. vs PV mult.", HistType::kTH2D, {{1000, 0, 5000, "FT0M mult."}, {100, 0, 100, "PV mult."}});

    histos.add("FT0A/percentile", "FT0A percentile.", HistType::kTH1D, {{binsPercentile, "FT0A percentile"}});
    histos.add("FT0A/percentilevsPV", "percentile vs PV mult.", HistType::kTH2D, {{binsPercentile, "FT0A percentile"}, {100, 0, 100, "PV mult."}});
    histos.add("FT0A/MultvsPV", "FT0A mult. vs PV mult.", HistType::kTH2D, {{1000, 0, 5000, "FT0A mult."}, {100, 0, 100, "PV mult."}});

    histos.add("FT0C/percentile", "FT0C percentile.", HistType::kTH1D, {{binsPercentile, "FT0C percentile"}});
    histos.add("FT0C/percentilevsPV", "percentile vs PV mult.", HistType::kTH2D, {{binsPercentile, "FT0C percentile"}, {100, 0, 100, "PV mult."}});
    histos.add("FT0C/MultvsPV", "FT0C mult. vs PV mult.", HistType::kTH2D, {{1000, 0, 5000, "FT0C mult."}, {100, 0, 100, "PV mult."}});

    histos.addClone("FT0M/", "sel8FT0M/");
    histos.addClone("FT0C/", "sel8FT0C/");
    histos.addClone("FT0A/", "sel8FT0A/");

    histos.print();
    listCalib.setObject(new TList);
  }

  using BCsWithTimestampsAndRun2Infos = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

  void processRun2(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, BCsWithTimestampsAndRun2Infos const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<BCsWithTimestampsAndRun2Infos>();
    if (bc.runNumber() != mRunNumber) {
      mRunNumber = bc.runNumber(); // mark that this run has been attempted already regardless of outcome
      LOGF(debug, "timestamp=%llu", bc.timestamp());
      TList* callst = nullptr;
      if (ccdbConfig.reconstructionPass.value == "") {
        callst = ccdb->getForRun<TList>(ccdbConfig.ccdbPath, bc.runNumber());
      } else if (ccdbConfig.reconstructionPass.value == "metadata") {
        std::map<std::string, std::string> metadata;
        metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
        LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
        callst = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, bc.runNumber(), metadata);
      } else {
        std::map<std::string, std::string> metadata;
        metadata["RecoPassName"] = ccdbConfig.reconstructionPass.value;
        LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", ccdbConfig.reconstructionPass.value);
        callst = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, bc.runNumber(), metadata);
      }

      Run2V0MInfo.mCalibrationStored = false;
      Run2V0AInfo.mCalibrationStored = false;
      Run2SPDTksInfo.mCalibrationStored = false;
      Run2SPDClsInfo.mCalibrationStored = false;
      Run2CL0Info.mCalibrationStored = false;
      Run2CL1Info.mCalibrationStored = false;
      if (callst != nullptr) {
        auto getccdb = [callst](const char* ccdbhname) {
          TH1* h = reinterpret_cast<TH1*>(callst->FindObject(ccdbhname));
          return h;
        };
        auto getformulaccdb = [callst](const char* ccdbhname) {
          TFormula* f = reinterpret_cast<TFormula*>(callst->FindObject(ccdbhname));
          return f;
        };
        if (isTableEnabled[kCentRun2V0Ms]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2V0MInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          Run2V0MInfo.mhVtxAmpCorrV0C = getccdb("hVtx_fAmplitude_V0C_Normalized");
          Run2V0MInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0M");
          Run2V0MInfo.mMCScale = getformulaccdb(TString::Format("%s-V0M", ccdbConfig.genName->c_str()).Data());
          if ((Run2V0MInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0MInfo.mhVtxAmpCorrV0C != nullptr) && (Run2V0MInfo.mhMultSelCalib != nullptr)) {
            if (ccdbConfig.genName->length() != 0) {
              if (Run2V0MInfo.mMCScale != nullptr) {
                for (int ixpar = 0; ixpar < 6; ++ixpar) {
                  Run2V0MInfo.mMCScalePars[ixpar] = Run2V0MInfo.mMCScale->GetParameter(ixpar);
                }
              } else {
                if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
                  LOGF(fatal, "MC Scale information from V0M for run %d not available", bc.runNumber());
                } else { // only if asked: continue filling with non-valid values (105)
                  LOGF(info, "MC Scale information from V0M for run %d not available", bc.runNumber());
                }
              }
            }
            Run2V0MInfo.mCalibrationStored = true;
          } else {
            if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
              LOGF(fatal, "Calibration information from V0M for run %d corrupted", bc.runNumber());
            } else { // only if asked: continue filling with non-valid values (105)
              LOGF(info, "Calibration information from V0M for run %d corrupted, will fill V0M tables with dummy values", bc.runNumber());
            }
          }
        }
        if (isTableEnabled[kCentRun2V0As]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2V0AInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          Run2V0AInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0A");
          if ((Run2V0AInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0AInfo.mhMultSelCalib != nullptr)) {
            Run2V0AInfo.mCalibrationStored = true;
          } else {
            if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
              LOGF(fatal, "Calibration information from V0A for run %d corrupted", bc.runNumber());
            } else { // only if asked: continue filling with non-valid values (105)
              LOGF(info, "Calibration information from V0A for run %d corrupted, will fill V0A tables with dummy values", bc.runNumber());
            }
          }
        }
        if (isTableEnabled[kCentRun2SPDTrks]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2SPDTksInfo.mhVtxAmpCorr = getccdb("hVtx_fnTracklets_Normalized");
          Run2SPDTksInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDTracklets");
          if ((Run2SPDTksInfo.mhVtxAmpCorr != nullptr) && (Run2SPDTksInfo.mhMultSelCalib != nullptr)) {
            Run2SPDTksInfo.mCalibrationStored = true;
          } else {
            if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
              LOGF(fatal, "Calibration information from SPD tracklets for run %d corrupted", bc.runNumber());
            } else { // only if asked: continue filling with non-valid values (105)
              LOGF(info, "Calibration information from SPD tracklets for run %d corrupted, will fill SPD tracklets tables with dummy values", bc.runNumber());
            }
          }
        }
        if (isTableEnabled[kCentRun2SPDClss]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2SPDClsInfo.mhVtxAmpCorrCL0 = getccdb("hVtx_fnSPDClusters0_Normalized");
          Run2SPDClsInfo.mhVtxAmpCorrCL1 = getccdb("hVtx_fnSPDClusters1_Normalized");
          Run2SPDClsInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDClusters");
          if ((Run2SPDClsInfo.mhVtxAmpCorrCL0 != nullptr) && (Run2SPDClsInfo.mhVtxAmpCorrCL1 != nullptr) && (Run2SPDClsInfo.mhMultSelCalib != nullptr)) {
            Run2SPDClsInfo.mCalibrationStored = true;
          } else {
            if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
              LOGF(fatal, "Calibration information from SPD clusters for run %d corrupted", bc.runNumber());
            } else { // only if asked: continue filling with non-valid values (105)
              LOGF(info, "Calibration information from SPD clusters for run %d corrupted, will fill SPD clusters tables with dummy values", bc.runNumber());
            }
          }
        }
        if (isTableEnabled[kCentRun2CL0s]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2CL0Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters0_Normalized");
          Run2CL0Info.mhMultSelCalib = getccdb("hMultSelCalib_CL0");
          if ((Run2CL0Info.mhVtxAmpCorr != nullptr) && (Run2CL0Info.mhMultSelCalib != nullptr)) {
            Run2CL0Info.mCalibrationStored = true;
          } else {
            if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
              LOGF(fatal, "Calibration information from CL0 multiplicity for run %d corrupted", bc.runNumber());
            } else { // only if asked: continue filling with non-valid values (105)
              LOGF(info, "Calibration information from CL0 multiplicity for run %d corrupted, will fill CL0 multiplicity tables with dummy values", bc.runNumber());
            }
          }
        }
        if (isTableEnabled[kCentRun2CL1s]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2CL1Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters1_Normalized");
          Run2CL1Info.mhMultSelCalib = getccdb("hMultSelCalib_CL1");
          if ((Run2CL1Info.mhVtxAmpCorr != nullptr) && (Run2CL1Info.mhMultSelCalib != nullptr)) {
            Run2CL1Info.mCalibrationStored = true;
          } else {
            if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
              LOGF(fatal, "Calibration information from CL1 multiplicity for run %d corrupted", bc.runNumber());
            } else { // only if asked: continue filling with non-valid values (105)
              LOGF(info, "Calibration information from CL1 multiplicity for run %d corrupted, will fill CL1 multiplicity tables with dummy values", bc.runNumber());
            }
          }
        }
      } else {
        if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
          LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        } else { // only if asked: continue filling with non-valid values (105)
          LOGF(info, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", bc.runNumber(), bc.timestamp());
        }
      }
    }

    auto scaleMC = [](float x, float pars[6]) {
      return std::pow(((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
    };

    if (isTableEnabled[kCentRun2V0Ms]) {
      float cV0M = 105.0f;
      if (Run2V0MInfo.mCalibrationStored) {
        float v0m;
        if (Run2V0MInfo.mMCScale != nullptr) {
          v0m = scaleMC(collision.multFV0M(), Run2V0MInfo.mMCScalePars);
          LOGF(debug, "Unscaled v0m: %f, scaled v0m: %f", collision.multFV0M(), v0m);
        } else {
          v0m = collision.multFV0A() * Run2V0MInfo.mhVtxAmpCorrV0A->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0A->FindFixBin(collision.posZ())) +
                collision.multFV0C() * Run2V0MInfo.mhVtxAmpCorrV0C->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0C->FindFixBin(collision.posZ()));
        }
        cV0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));
      }
      LOGF(debug, "centRun2V0M=%.0f", cV0M);
      // fill centrality columns
      centRun2V0M(cV0M);
    }
    if (isTableEnabled[kCentRun2V0As]) {
      float cV0A = 105.0f;
      if (Run2V0AInfo.mCalibrationStored) {
        float v0a = collision.multFV0A() * Run2V0AInfo.mhVtxAmpCorrV0A->GetBinContent(Run2V0AInfo.mhVtxAmpCorrV0A->FindFixBin(collision.posZ()));
        cV0A = Run2V0AInfo.mhMultSelCalib->GetBinContent(Run2V0AInfo.mhMultSelCalib->FindFixBin(v0a));
      }
      LOGF(debug, "centRun2V0A=%.0f", cV0A);
      // fill centrality columns
      centRun2V0A(cV0A);
    }
    if (isTableEnabled[kCentRun2SPDTrks]) {
      float cSPD = 105.0f;
      if (Run2SPDTksInfo.mCalibrationStored) {
        float spdm = collision.multTracklets() * Run2SPDTksInfo.mhVtxAmpCorr->GetBinContent(Run2SPDTksInfo.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cSPD = Run2SPDTksInfo.mhMultSelCalib->GetBinContent(Run2SPDTksInfo.mhMultSelCalib->FindFixBin(spdm));
      }
      LOGF(debug, "centSPDTracklets=%.0f", cSPD);
      centRun2SPDTracklets(cSPD);
    }
    if (isTableEnabled[kCentRun2SPDClss]) {
      float cSPD = 105.0f;
      if (Run2SPDClsInfo.mCalibrationStored) {
        float spdm = bc.spdClustersL0() * Run2SPDClsInfo.mhVtxAmpCorrCL0->GetBinContent(Run2SPDClsInfo.mhVtxAmpCorrCL0->FindFixBin(collision.posZ())) +
                     bc.spdClustersL1() * Run2SPDClsInfo.mhVtxAmpCorrCL1->GetBinContent(Run2SPDClsInfo.mhVtxAmpCorrCL1->FindFixBin(collision.posZ()));
        cSPD = Run2SPDClsInfo.mhMultSelCalib->GetBinContent(Run2SPDClsInfo.mhMultSelCalib->FindFixBin(spdm));
      }
      LOGF(debug, "centSPDClusters=%.0f", cSPD);
      centRun2SPDClusters(cSPD);
    }
    if (isTableEnabled[kCentRun2CL0s]) {
      float cCL0 = 105.0f;
      if (Run2CL0Info.mCalibrationStored) {
        float cl0m = bc.spdClustersL0() * Run2CL0Info.mhVtxAmpCorr->GetBinContent(Run2CL0Info.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cCL0 = Run2CL0Info.mhMultSelCalib->GetBinContent(Run2CL0Info.mhMultSelCalib->FindFixBin(cl0m));
      }
      LOGF(debug, "centCL0=%.0f", cCL0);
      centRun2CL0(cCL0);
    }
    if (isTableEnabled[kCentRun2CL1s]) {
      float cCL1 = 105.0f;
      if (Run2CL1Info.mCalibrationStored) {
        float cl1m = bc.spdClustersL1() * Run2CL1Info.mhVtxAmpCorr->GetBinContent(Run2CL1Info.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cCL1 = Run2CL1Info.mhMultSelCalib->GetBinContent(Run2CL1Info.mhMultSelCalib->FindFixBin(cl1m));
      }
      LOGF(debug, "centCL1=%.0f", cCL1);
      centRun2CL1(cCL1);
    }
  }

  using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

  template <bool enableCentFV0 = true,
            bool enableCentFT0 = true,
            bool enableCentFDD = true,
            bool enableCentNTPV = true,
            bool enableCentNGlobal = false,
            bool enableCentMFT = false,
            typename CollisionType>
  void produceRun3Tables(CollisionType const& collisions)
  {
    // do memory reservation for the relevant tables only, please
    for (auto const& table : mEnabledTables) {
      switch (table) {
        case kCentFV0As:
          centFV0A.reserve(collisions.size());
          break;
        case kCentFT0Ms:
          centFT0M.reserve(collisions.size());
          break;
        case kCentFT0As:
          centFT0A.reserve(collisions.size());
          break;
        case kCentFT0Cs:
          centFT0C.reserve(collisions.size());
          break;
        case kCentFT0CVariant1s:
          centFT0CVariant1.reserve(collisions.size());
          break;
        case kCentFDDMs:
          centFDDM.reserve(collisions.size());
          break;
        case kCentNTPVs:
          centNTPV.reserve(collisions.size());
          break;
        case kCentNGlobals:
          centNGlobals.reserve(collisions.size());
          break;
        case kCentMFTs:
          centMFTs.reserve(collisions.size());
          break;
        default:
          LOGF(fatal, "Table %d not supported in Run3", table);
          break;
      }
    }

    for (auto const& collision : collisions) {
      /* check the previous run number */
      auto bc = collision.template bc_as<BCsWithTimestamps>();
      if (bc.runNumber() != mRunNumber) {
        mRunNumber = bc.runNumber(); // mark that this run has been attempted already regardless of outcome
        LOGF(info, "timestamp=%llu, run number=%d", bc.timestamp(), bc.runNumber());
        TList* callst = nullptr;
        // Check if the ccdb path is a root file
        if (ccdbConfig.ccdbPath.value.find(".root") != std::string::npos) {
          TFile f(ccdbConfig.ccdbPath.value.c_str(), "READ");
          f.GetObject(ccdbConfig.reconstructionPass.value.c_str(), callst);
          if (!callst) {
            f.ls();
            LOG(fatal) << "No calibration list " << ccdbConfig.reconstructionPass.value << " found in the file " << ccdbConfig.ccdbPath.value;
          }
        } else {
          if (ccdbConfig.reconstructionPass.value == "") {
            callst = ccdb->getForRun<TList>(ccdbConfig.ccdbPath, bc.runNumber());
          } else if (ccdbConfig.reconstructionPass.value == "metadata") {
            std::map<std::string, std::string> metadata;
            metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
            LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
            callst = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, bc.runNumber(), metadata);
          } else {
            std::map<std::string, std::string> metadata;
            metadata["RecoPassName"] = ccdbConfig.reconstructionPass.value;
            LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", ccdbConfig.reconstructionPass.value);
            callst = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, bc.runNumber(), metadata);
          }
        }

        fv0aInfo.mCalibrationStored = false;
        ft0mInfo.mCalibrationStored = false;
        ft0aInfo.mCalibrationStored = false;
        ft0cInfo.mCalibrationStored = false;
        ft0cVariant1Info.mCalibrationStored = false;
        fddmInfo.mCalibrationStored = false;
        ntpvInfo.mCalibrationStored = false;
        nGlobalInfo.mCalibrationStored = false;
        mftInfo.mCalibrationStored = false;
        if (callst != nullptr) {
          if (produceHistograms) {
            listCalib->Add(callst->Clone(Form("%i", bc.runNumber())));
          }
          LOGF(info, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          auto getccdb = [callst, bc](struct CalibrationInfo& estimator, const Configurable<std::string> generatorName) { // TODO: to consider the name inside the estimator structure
            estimator.mhMultSelCalib = reinterpret_cast<TH1*>(callst->FindObject(TString::Format("hCalibZeq%s", estimator.name.c_str()).Data()));
            estimator.mMCScale = reinterpret_cast<TFormula*>(callst->FindObject(TString::Format("%s-%s", generatorName->c_str(), estimator.name.c_str()).Data()));
            if (estimator.mhMultSelCalib != nullptr) {
              if (generatorName->length() != 0) {
                LOGF(info, "Retrieving MC calibration for %d, generator name: %s", bc.runNumber(), generatorName->c_str());
                if (estimator.mMCScale != nullptr) {
                  for (int ixpar = 0; ixpar < 6; ++ixpar) {
                    estimator.mMCScalePars[ixpar] = estimator.mMCScale->GetParameter(ixpar);
                    LOGF(info, "Parameter index %i value %.5f", ixpar, estimator.mMCScalePars[ixpar]);
                  }
                } else {
                  LOGF(warning, "MC Scale information from %s for run %d not available", estimator.name.c_str(), bc.runNumber());
                }
              }
              estimator.mCalibrationStored = true;
              estimator.isSane();
            } else {
              LOGF(error, "Calibration information from %s for run %d not available", estimator.name.c_str(), bc.runNumber());
            }
          };

          for (auto const& table : mEnabledTables) {
            switch (table) {
              case kCentFV0As:
                getccdb(fv0aInfo, ccdbConfig.genName);
                break;
              case kCentFT0Ms:
                getccdb(ft0mInfo, ccdbConfig.genName);
                break;
              case kCentFT0As:
                getccdb(ft0aInfo, ccdbConfig.genName);
                break;
              case kCentFT0Cs:
                getccdb(ft0cInfo, ccdbConfig.genName);
                break;
              case kCentFT0CVariant1s:
                getccdb(ft0cVariant1Info, ccdbConfig.genName);
                break;
              case kCentFDDMs:
                getccdb(fddmInfo, ccdbConfig.genName);
                break;
              case kCentNTPVs:
                getccdb(ntpvInfo, ccdbConfig.genName);
                break;
              case kCentNGlobals:
                getccdb(nGlobalInfo, ccdbConfig.genName);
                break;
              case kCentMFTs:
                getccdb(mftInfo, ccdbConfig.genName);
                break;
              default:
                LOGF(fatal, "Table %d not supported in Run3", table);
                break;
            }
          }
        } else {
          if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
            LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
          } else { // only if asked: continue filling with non-valid values (105)
            LOGF(info, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", bc.runNumber(), bc.timestamp());
          }
        }
      }

      /**
       * @brief Populates a table with data based on the given calibration information and multiplicity.
       *
       * @param table The table to populate.
       * @param estimator The calibration information.
       * @param multiplicity The multiplicity value.
       */

      auto populateTable = [&](auto& table, struct CalibrationInfo& estimator, float multiplicity) {
        const bool assignOutOfRange = embedINELgtZEROselection && !collision.isInelGt0();
        auto scaleMC = [](float x, float pars[6]) {
          return std::pow(((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
        };

        float percentile = 105.0f;
        float scaledMultiplicity = multiplicity;
        if (estimator.mCalibrationStored) {
          if (estimator.mMCScale != nullptr) {
            scaledMultiplicity = scaleMC(multiplicity, estimator.mMCScalePars);
            LOGF(debug, "Unscaled %s multiplicity: %f, scaled %s multiplicity: %f", estimator.name.c_str(), multiplicity, estimator.name.c_str(), scaledMultiplicity);
          }
          percentile = estimator.mhMultSelCalib->GetBinContent(estimator.mhMultSelCalib->FindFixBin(scaledMultiplicity));
          if (assignOutOfRange)
            percentile = 100.5f;
        }
        LOGF(debug, "%s centrality/multiplicity percentile = %.0f for a zvtx eq %s value %.0f", estimator.name.c_str(), percentile, estimator.name.c_str(), scaledMultiplicity);
        table(percentile);
        return percentile;
      };

      for (auto const& table : mEnabledTables) {
        switch (table) {
          case kCentFV0As:
            if constexpr (enableCentFV0) {
              populateTable(centFV0A, fv0aInfo, collision.multZeqFV0A());
            }
            break;
          case kCentFT0Ms:
            if constexpr (enableCentFT0) {
              const float perC = populateTable(centFT0M, ft0mInfo, collision.multZeqFT0A() + collision.multZeqFT0C());
              if (produceHistograms.value) {
                histos.fill(HIST("FT0M/percentile"), perC);
                histos.fill(HIST("FT0M/percentilevsPV"), perC, collision.multNTracksPV());
                histos.fill(HIST("FT0M/MultvsPV"), collision.multZeqFT0A() + collision.multZeqFT0C(), collision.multNTracksPV());
                if (collision.sel8()) {
                  histos.fill(HIST("sel8FT0M/percentile"), perC);
                  histos.fill(HIST("sel8FT0M/percentilevsPV"), perC, collision.multNTracksPV());
                  histos.fill(HIST("sel8FT0M/MultvsPV"), collision.multZeqFT0A() + collision.multZeqFT0C(), collision.multNTracksPV());
                }
              }
            }
            break;
          case kCentFT0As:
            if constexpr (enableCentFT0) {
              const float perC = populateTable(centFT0A, ft0aInfo, collision.multZeqFT0A());
              if (produceHistograms.value) {
                histos.fill(HIST("FT0A/percentile"), perC);
                histos.fill(HIST("FT0A/percentilevsPV"), perC, collision.multNTracksPV());
                histos.fill(HIST("FT0A/MultvsPV"), collision.multZeqFT0A() + collision.multZeqFT0C(), collision.multNTracksPV());
                if (collision.sel8()) {
                  histos.fill(HIST("sel8FT0A/percentile"), perC);
                  histos.fill(HIST("sel8FT0A/percentilevsPV"), perC, collision.multNTracksPV());
                  histos.fill(HIST("sel8FT0A/MultvsPV"), collision.multZeqFT0A() + collision.multZeqFT0C(), collision.multNTracksPV());
                }
              }
            }
            break;
          case kCentFT0Cs:
            if constexpr (enableCentFT0) {
              const float perC = populateTable(centFT0C, ft0cInfo, collision.multZeqFT0C());
              if (produceHistograms.value) {
                histos.fill(HIST("FT0C/percentile"), perC);
                histos.fill(HIST("FT0C/percentilevsPV"), perC, collision.multNTracksPV());
                histos.fill(HIST("FT0C/MultvsPV"), collision.multZeqFT0A() + collision.multZeqFT0C(), collision.multNTracksPV());
                if (collision.sel8()) {
                  histos.fill(HIST("sel8FT0C/percentile"), perC);
                  histos.fill(HIST("sel8FT0C/percentilevsPV"), perC, collision.multNTracksPV());
                  histos.fill(HIST("sel8FT0C/MultvsPV"), collision.multZeqFT0A() + collision.multZeqFT0C(), collision.multNTracksPV());
                }
              }
            }
            break;
          case kCentFT0CVariant1s:
            if constexpr (enableCentFT0) {
              populateTable(centFT0CVariant1, ft0cVariant1Info, collision.multZeqFT0C());
            }
            break;
          case kCentFDDMs:
            if constexpr (enableCentFDD) {
              populateTable(centFDDM, fddmInfo, collision.multZeqFDDA() + collision.multZeqFDDC());
            }
            break;
          case kCentNTPVs:
            if constexpr (enableCentNTPV) {
              populateTable(centNTPV, ntpvInfo, collision.multZeqNTracksPV());
            }
            break;
          case kCentNGlobals:
            if constexpr (enableCentNGlobal) {
              populateTable(centNGlobals, nGlobalInfo, collision.multNTracksGlobal());
            }
            break;
          case kCentMFTs:
            if constexpr (enableCentMFT) {
              populateTable(centMFTs, mftInfo, collision.mftNtracks());
            }
            break;
          default:
            LOGF(fatal, "Table %d not supported in Run3", table);
            break;
        }
      }
    }
  }

  void processRun3Complete(soa::Join<aod::Collisions, aod::PVMults, aod::MultZeqs, aod::EvSels, aod::MultsGlobal, aod::MFTMults> const& collisions, BCsWithTimestamps const&)
  {
    produceRun3Tables<true, true, true, true, true, true>(collisions);
  }

  void processRun3(soa::Join<aod::Collisions, aod::PVMults, aod::MultZeqs, aod::EvSels> const& collisions, BCsWithTimestamps const&)
  {
    produceRun3Tables(collisions);
  }

  void processRun3FT0(soa::Join<aod::Collisions, aod::PVMults, aod::FT0MultZeqs, aod::PVMultZeqs, aod::EvSels> const& collisions, BCsWithTimestamps const&)
  {
    produceRun3Tables<false, // FV0
                      true,  // FT0
                      false, // PV
                      false>(collisions);
  }

  // Process switches
  PROCESS_SWITCH(CentralityTable, processRun2, "Provide Run2 calibrated centrality/multiplicity percentiles tables", true);
  PROCESS_SWITCH(CentralityTable, processRun3Complete, "Provide Run3 calibrated centrality/multiplicity percentiles tables using MFT and global tracks (requires extra subscriptions)", false);
  PROCESS_SWITCH(CentralityTable, processRun3, "Provide Run3 calibrated centrality/multiplicity percentiles tables", false);
  PROCESS_SWITCH(CentralityTable, processRun3FT0, "Provide Run3 calibrated centrality/multiplicity percentiles tables for FT0 only", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<CentralityTable>(cfgc)};
}
