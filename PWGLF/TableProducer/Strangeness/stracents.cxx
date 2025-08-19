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
/// \file straCents.cxx
/// \brief Task to re-produce the strangeness centrality tables, repeating what has been done the centralityTable.cxx and multiplicityTable.cxx tasks
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
#include "Framework/RunningWorkflowInfo.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "MetadataHelper.h"
#include "TableHelper.h"
#include "TList.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

struct straCents {
  Produces<aod::StraCents> strangeCents;
  Produces<aod::StraCentsRun2> strangeCentsRun2;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> doVertexZeq{"doVertexZeq", 1, "if 1: do vertex Z eq mult table"};
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
    Configurable<std::string> ccdbPath_Centrality{"ccdbPath_Centrality", "Centrality/Estimators", "The CCDB path for centrality/multiplicity information"};
    Configurable<std::string> ccdbPath_Multiplicity{"ccdbPath_Multiplicity", "Centrality/Calibration", "The CCDB path for centrality/multiplicity information"};
    Configurable<std::string> genName{"genName", "", "Genearator name: HIJING, PYTHIA8, ... Default: \"\""};
    Configurable<bool> doNotCrashOnNull{"doNotCrashOnNull", false, {"Option to not crash on null and instead fill required tables with dummy info"}};
    Configurable<std::string> reconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};
  } ccdbConfig;

  Configurable<bool> embedINELgtZEROselection{"embedINELgtZEROselection", false, {"Option to do percentile 100.5 if not INELgtZERO"}};
  Configurable<bool> produceHistograms{"produceHistograms", false, {"Option to produce debug histograms"}};
  ConfigurableAxis binsPercentile{"binsPercentile", {VARIABLE_WIDTH, 0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0}, "Binning of the percentile axis"};

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

  int mRunNumber;
  bool lCalibLoaded;
  TProfile* hVtxZFV0A;
  TProfile* hVtxZFT0A;
  TProfile* hVtxZFT0C;
  TProfile* hVtxZFDDA;

  TProfile* hVtxZFDDC;
  TProfile* hVtxZNTracks;

  // Debug output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TList> listCalib{"calib-list", OutputObjHandlingPolicy::QAObject};

  void init(InitContext&)
  {
    LOG(info) << "Initializing centrality table producer";

    if (metadataInfo.isFullyDefined() && !doprocessRun2 && !doprocessRun3) { // Check if the metadata is initialized (only if not forced from the workflow configuration)
      if (metadataInfo.isRun3()) {
        doprocessRun3.value = true;
      } else {
        doprocessRun2.value = false;
      }
    }

    if (doprocessRun2 == false && doprocessRun3 == false) {
      LOGF(fatal, "Neither processRun2 nor processRun3 nor processRun3Complete enabled. Please choose one.");
    }
    if (doprocessRun2 == true && doprocessRun3 == true) {
      LOGF(fatal, "Cannot enable processRun2 and processRun3 at the same time. Please choose one.");
    }

    mRunNumber = 0;
    lCalibLoaded = false;

    hVtxZFV0A = nullptr;
    hVtxZFT0A = nullptr;
    hVtxZFT0C = nullptr;
    hVtxZFDDA = nullptr;
    hVtxZFDDC = nullptr;
    hVtxZNTracks = nullptr;

    ccdb->setURL(ccdbConfig.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false); // don't fatal, please - exception is caught explicitly (as it should)

    if (produceHistograms) {
      histos.add("FT0M/Mult", "FT0M mult.", HistType::kTH1D, {{1000, 0, 5000, "FT0M mult."}});
      histos.add("FT0M/percentile", "FT0M percentile.", HistType::kTH1D, {{binsPercentile, "FT0M percentile"}});
      histos.add("FT0M/percentilevsPV", "percentile vs PV mult.", HistType::kTH2D, {{binsPercentile, "FT0M percentile"}, {100, 0, 100, "PV mult."}});
      histos.add("FT0M/MultvsPV", "FT0M mult. vs PV mult.", HistType::kTH2D, {{1000, 0, 5000, "FT0M mult."}, {100, 0, 100, "PV mult."}});

      histos.add("FT0A/Mult", "FT0A mult", HistType::kTH1D, {{1000, 0, 1000, "FT0A multiplicity"}});
      histos.add("FT0A/percentile", "FT0A percentile.", HistType::kTH1D, {{binsPercentile, "FT0A percentile"}});
      histos.add("FT0A/percentilevsPV", "percentile vs PV mult.", HistType::kTH2D, {{binsPercentile, "FT0A percentile"}, {100, 0, 100, "PV mult."}});
      histos.add("FT0A/MultvsPV", "FT0A mult. vs PV mult.", HistType::kTH2D, {{1000, 0, 5000, "FT0A mult."}, {100, 0, 100, "PV mult."}});

      histos.add("FT0C/Mult", "FT0C mult", HistType::kTH1D, {{1000, 0, 5000, "FT0C multiplicity"}});
      histos.add("FT0C/percentile", "FT0C percentile.", HistType::kTH1D, {{binsPercentile, "FT0C percentile"}});
      histos.add("FT0C/percentilevsPV", "percentile vs PV mult.", HistType::kTH2D, {{binsPercentile, "FT0C percentile"}, {100, 0, 100, "PV mult."}});
      histos.add("FT0C/MultvsPV", "FT0C mult. vs PV mult.", HistType::kTH2D, {{1000, 0, 5000, "FT0C mult."}, {100, 0, 100, "PV mult."}});

      histos.addClone("FT0M/", "sel8FT0M/");
      histos.addClone("FT0C/", "sel8FT0C/");
      histos.addClone("FT0A/", "sel8FT0A/");

      histos.print();
    }
    listCalib.setObject(new TList);
  }

  template <typename TCollision>
  void initCCDB(TCollision collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    mRunNumber = collision.runNumber();
    LOGF(info, "timestamp=%llu, run number=%d", collision.timestamp(), collision.runNumber());

    TList* lCalibObjects_Centrality = nullptr;
    TList* lCalibObjects_Multiplicity = nullptr;
    if (ccdbConfig.reconstructionPass.value == "") {
      lCalibObjects_Centrality = ccdb->getForRun<TList>(ccdbConfig.ccdbPath_Centrality, mRunNumber);
      lCalibObjects_Multiplicity = ccdb->getForRun<TList>(ccdbConfig.ccdbPath_Multiplicity, mRunNumber);
    } else if (ccdbConfig.reconstructionPass.value == "metadata") {
      std::map<std::string, std::string> metadata;
      metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
      LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
      lCalibObjects_Centrality = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath_Centrality, mRunNumber, metadata);
      lCalibObjects_Multiplicity = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath_Multiplicity, mRunNumber, metadata);
    } else {
      std::map<std::string, std::string> metadata;
      metadata["RecoPassName"] = ccdbConfig.reconstructionPass.value;
      LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", ccdbConfig.reconstructionPass.value);
      lCalibObjects_Centrality = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath_Centrality, mRunNumber, metadata);
      lCalibObjects_Multiplicity = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath_Multiplicity, mRunNumber, metadata);
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

    Run2V0MInfo.mCalibrationStored = false;
    Run2V0AInfo.mCalibrationStored = false;
    Run2SPDTksInfo.mCalibrationStored = false;
    Run2SPDClsInfo.mCalibrationStored = false;
    Run2CL0Info.mCalibrationStored = false;
    Run2CL1Info.mCalibrationStored = false;

    if (lCalibObjects_Centrality) {
      if (produceHistograms) {
        listCalib->Add(lCalibObjects_Centrality->Clone(Form("Centrality_%i", collision.runNumber())));
      }

      LOGF(info, "Getting new histograms with %d run number for %d run number", mRunNumber, collision.runNumber());

      if constexpr (requires { collision.sel7(); }) { // check if we are in Run 2
        auto getccdb = [lCalibObjects_Centrality](const char* ccdbhname) {
          TH1* h = reinterpret_cast<TH1*>(lCalibObjects_Centrality->FindObject(ccdbhname));
          return h;
        };
        auto getformulaccdb = [lCalibObjects_Centrality](const char* ccdbhname) {
          TFormula* f = reinterpret_cast<TFormula*>(lCalibObjects_Centrality->FindObject(ccdbhname));
          return f;
        };

        // V0M
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
                LOGF(fatal, "MC Scale information from V0M for run %d not available", collision.runNumber());
              } else { // only if asked: continue filling with non-valid values (105)
                LOGF(warning, "MC Scale information from V0M for run %d not available", collision.runNumber());
              }
            }
          }
          Run2V0MInfo.mCalibrationStored = true;
        } else {
          if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
            LOGF(fatal, "Calibration information from V0M for run %d corrupted", collision.runNumber());
          } else { // only if asked: continue filling with non-valid values (105)
            LOGF(warning, "Calibration information from V0M for run %d corrupted, will fill V0M tables with dummy values", collision.runNumber());
          }
        }

        // V0A
        Run2V0AInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
        Run2V0AInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0A");
        if ((Run2V0AInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0AInfo.mhMultSelCalib != nullptr)) {
          Run2V0AInfo.mCalibrationStored = true;
        } else {
          if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
            LOGF(fatal, "Calibration information from V0A for run %d corrupted", collision.runNumber());
          } else { // only if asked: continue filling with non-valid values (105)
            LOGF(warning, "Calibration information from V0A for run %d corrupted, will fill V0A tables with dummy values", collision.runNumber());
          }
        }

        // SPD tracklets
        Run2SPDTksInfo.mhVtxAmpCorr = getccdb("hVtx_fnTracklets_Normalized");
        Run2SPDTksInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDTracklets");
        if ((Run2SPDTksInfo.mhVtxAmpCorr != nullptr) && (Run2SPDTksInfo.mhMultSelCalib != nullptr)) {
          Run2SPDTksInfo.mCalibrationStored = true;
        } else {
          if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
            LOGF(fatal, "Calibration information from SPD tracklets for run %d corrupted", collision.runNumber());
          } else { // only if asked: continue filling with non-valid values (105)
            LOGF(warning, "Calibration information from SPD tracklets for run %d corrupted, will fill SPD tracklets tables with dummy values", collision.runNumber());
          }
        }

        // SPD clusters
        Run2SPDClsInfo.mhVtxAmpCorrCL0 = getccdb("hVtx_fnSPDClusters0_Normalized");
        Run2SPDClsInfo.mhVtxAmpCorrCL1 = getccdb("hVtx_fnSPDClusters1_Normalized");
        Run2SPDClsInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDClusters");
        if ((Run2SPDClsInfo.mhVtxAmpCorrCL0 != nullptr) && (Run2SPDClsInfo.mhVtxAmpCorrCL1 != nullptr) && (Run2SPDClsInfo.mhMultSelCalib != nullptr)) {
          Run2SPDClsInfo.mCalibrationStored = true;
        } else {
          if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
            LOGF(fatal, "Calibration information from SPD clusters for run %d corrupted", collision.runNumber());
          } else { // only if asked: continue filling with non-valid values (105)
            LOGF(warning, "Calibration information from SPD clusters for run %d corrupted, will fill SPD clusters tables with dummy values", collision.runNumber());
          }
        }
      } else {
        // we are in Run 3
        auto getccdb = [lCalibObjects_Centrality, collision](struct CalibrationInfo& estimator, const Configurable<std::string> generatorName, const Configurable<bool> notCrashOnNull) { // TODO: to consider the name inside the estimator structure
          estimator.mhMultSelCalib = reinterpret_cast<TH1*>(lCalibObjects_Centrality->FindObject(TString::Format("hCalibZeq%s", estimator.name.c_str()).Data()));
          estimator.mMCScale = reinterpret_cast<TFormula*>(lCalibObjects_Centrality->FindObject(TString::Format("%s-%s", generatorName->c_str(), estimator.name.c_str()).Data()));
          if (estimator.mhMultSelCalib != nullptr) {
            if (generatorName->length() != 0) {
              LOGF(info, "Retrieving MC calibration for %d, generator name: %s", collision.runNumber(), generatorName->c_str());
              if (estimator.mMCScale != nullptr) {
                for (int ixpar = 0; ixpar < 6; ++ixpar) {
                  estimator.mMCScalePars[ixpar] = estimator.mMCScale->GetParameter(ixpar);
                  LOGF(info, "Parameter index %i value %.5f", ixpar, estimator.mMCScalePars[ixpar]);
                }
              } else {
                LOGF(warning, "MC Scale information from %s for run %d not available", estimator.name.c_str(), collision.runNumber());
              }
            }
            estimator.mCalibrationStored = true;
            estimator.isSane();
          } else {
            if (!notCrashOnNull.value) { // default behaviour: crash
              LOGF(error, "Calibration information from %s for run %d not available", estimator.name.c_str(), collision.runNumber());
            } else {
              LOGF(warning, "Calibration information from %s for run %d not available", estimator.name.c_str(), collision.runNumber());
            }
          }
        };
        getccdb(fv0aInfo, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);
        getccdb(ft0mInfo, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);
        getccdb(ft0aInfo, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);
        getccdb(ft0cInfo, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);
        getccdb(ft0cVariant1Info, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);
        getccdb(fddmInfo, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);
        getccdb(ntpvInfo, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);
        getccdb(nGlobalInfo, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);
        getccdb(mftInfo, ccdbConfig.genName, ccdbConfig.doNotCrashOnNull);

        if (lCalibObjects_Multiplicity) {
          if (produceHistograms) {
            listCalib->Add(lCalibObjects_Multiplicity->Clone(Form("Multiplicity_%i", collision.runNumber())));
          }

          hVtxZFV0A = static_cast<TProfile*>(lCalibObjects_Multiplicity->FindObject("hVtxZFV0A"));
          hVtxZFT0A = static_cast<TProfile*>(lCalibObjects_Multiplicity->FindObject("hVtxZFT0A"));
          hVtxZFT0C = static_cast<TProfile*>(lCalibObjects_Multiplicity->FindObject("hVtxZFT0C"));
          hVtxZFDDA = static_cast<TProfile*>(lCalibObjects_Multiplicity->FindObject("hVtxZFDDA"));
          hVtxZFDDC = static_cast<TProfile*>(lCalibObjects_Multiplicity->FindObject("hVtxZFDDC"));
          hVtxZNTracks = static_cast<TProfile*>(lCalibObjects_Multiplicity->FindObject("hVtxZNTracksPV"));
          lCalibLoaded = true;
          // Capture error
          if (!hVtxZFV0A || !hVtxZFT0A || !hVtxZFT0C || !hVtxZFDDA || !hVtxZFDDC || !hVtxZNTracks) {
            LOGF(error, "Problem loading CCDB objects! Please check");
            lCalibLoaded = false;
          }
        } else {
          if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
            LOGF(fatal, "Multiplicity calibration is not available in CCDB for run=%d at timestamp=%llu", collision.runNumber(), collision.timestamp());
          } else { // only if asked: continue filling with non-valid values (105)
            LOGF(warning, "Multiplicity calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", collision.runNumber(), collision.timestamp());
          }
          lCalibLoaded = false;
        }
      } // end we are in Run 3
    } else {
      if (!ccdbConfig.doNotCrashOnNull) { // default behaviour: crash
        LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", collision.runNumber(), collision.timestamp());
      } else { // only if asked: continue filling with non-valid values (105)
        LOGF(warning, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", collision.runNumber(), collision.timestamp());
      }
      lCalibLoaded = false;
    }
  }

  void processRun2(soa::Join<aod::StraCollisions, aod::StraEvSelsRun2, aod::StraStamps> const& collisions)
  {
    for (auto const& collision : collisions) {
      if (doVertexZeq > 0) {
        initCCDB(collision);
      }

      float centRun2V0M = 105.0f;
      float centRun2V0A = 105.0f;
      float centRun2SPDTrks = 105.0f;
      float centRun2SPDClss = 105.0f;

      auto scaleMC = [](float x, float pars[6]) {
        return std::pow(((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
      };

      // Run 2 V0M
      if (Run2V0MInfo.mCalibrationStored) {
        float v0m;
        if (Run2V0MInfo.mMCScale != nullptr) {
          v0m = scaleMC(collision.multFV0A() + collision.multFV0C(), Run2V0MInfo.mMCScalePars);
          LOGF(debug, "Unscaled v0m: %f, scaled v0m: %f", collision.multFV0A() + collision.multFV0C(), v0m);
        } else {
          v0m = collision.multFV0A() * Run2V0MInfo.mhVtxAmpCorrV0A->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0A->FindFixBin(collision.posZ())) +
                collision.multFV0C() * Run2V0MInfo.mhVtxAmpCorrV0C->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0C->FindFixBin(collision.posZ()));
        }
        centRun2V0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));
      }
      LOGF(debug, "centRun2V0M=%.0f", centRun2V0M);

      // Run 2 V0A
      if (Run2V0AInfo.mCalibrationStored) {
        float v0a = collision.multFV0A() * Run2V0AInfo.mhVtxAmpCorrV0A->GetBinContent(Run2V0AInfo.mhVtxAmpCorrV0A->FindFixBin(collision.posZ()));
        centRun2V0A = Run2V0AInfo.mhMultSelCalib->GetBinContent(Run2V0AInfo.mhMultSelCalib->FindFixBin(v0a));
      }
      LOGF(debug, "centRun2V0A=%.0f", centRun2V0A);

      // Run 2 SPD tracklets
      if (Run2SPDTksInfo.mCalibrationStored) {
        float spdm = collision.multTracklets() * Run2SPDTksInfo.mhVtxAmpCorr->GetBinContent(Run2SPDTksInfo.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        centRun2SPDTrks = Run2SPDTksInfo.mhMultSelCalib->GetBinContent(Run2SPDTksInfo.mhMultSelCalib->FindFixBin(spdm));
      }
      LOGF(debug, "centSPDTracklets=%.0f", centRun2SPDTrks);

      // Run 2 SPD Cls
      if (Run2SPDClsInfo.mCalibrationStored) {
        // spdClustersL0 and spdClustersL1 not available in strangeness data model
        float spdm = collision.spdClustersL0() * Run2SPDClsInfo.mhVtxAmpCorrCL0->GetBinContent(Run2SPDClsInfo.mhVtxAmpCorrCL0->FindFixBin(collision.posZ())) +
                     collision.spdClustersL1() * Run2SPDClsInfo.mhVtxAmpCorrCL1->GetBinContent(Run2SPDClsInfo.mhVtxAmpCorrCL1->FindFixBin(collision.posZ()));
        centRun2SPDClss = Run2SPDClsInfo.mhMultSelCalib->GetBinContent(Run2SPDClsInfo.mhMultSelCalib->FindFixBin(spdm));
      }
      LOGF(debug, "centSPDClusters=%.0f", centRun2SPDClss);

      strangeCentsRun2(centRun2V0M, centRun2V0A,
                       centRun2SPDTrks, centRun2SPDClss);
    }
  }

  void processRun3(soa::Join<aod::StraCollisions, aod::StraEvSels, aod::StraStamps> const& collisions)
  {
    for (auto const& collision : collisions) {
      if (doVertexZeq > 0) {
        initCCDB(collision);
      }

      float multZeqFV0A = 0.f;
      float multZeqFT0A = 0.f;
      float multZeqFT0C = 0.f;
      // float multZeqFDDA = 0.f;
      // float multZeqFDDC = 0.f;

      if (std::fabs(collision.posZ()) < 15.0f && hVtxZFV0A) {
        multZeqFV0A = hVtxZFV0A->Interpolate(0.0) * collision.multFV0A() / hVtxZFV0A->Interpolate(collision.posZ());
      }
      if (std::fabs(collision.posZ()) < 15.0f && hVtxZFT0A) {
        multZeqFT0A = hVtxZFT0A->Interpolate(0.0) * collision.multFT0A() / hVtxZFT0A->Interpolate(collision.posZ());
      }
      if (std::fabs(collision.posZ()) < 15.0f && hVtxZFT0C) {
        multZeqFT0C = hVtxZFT0C->Interpolate(0.0) * collision.multFT0C() / hVtxZFT0C->Interpolate(collision.posZ());
      }
      // if (std::fabs(collision.posZ()) < 15.0f && hVtxZFDDA) {
      //   multZeqFDDA = hVtxZFDDA->Interpolate(0.0) * collision.multFDDA() / hVtxZFDDA->Interpolate(collision.posZ());
      // }
      // if (std::fabs(collision.posZ()) < 15.0f && hVtxZFDDC) {
      //   multZeqFDDC = hVtxZFDDC->Interpolate(0.0) * collision.multFDDC() / hVtxZFDDC->Interpolate(collision.posZ());
      // }

      /**
       * @brief Get centrality value based on the given calibration information and multiplicity.
       *        Modified version of populateTable (https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/centralityTable.cxx#L648)
       *
       * @param estimator The calibration information.
       * @param multiplicity The multiplicity value.
       */

      auto getCentrality = [&](struct CalibrationInfo& estimator, float multiplicity) {
        const bool assignOutOfRange = embedINELgtZEROselection && !(collision.multNTracksPVeta1() > 0);
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
        return percentile;
      };

      float centFT0M = getCentrality(ft0mInfo, multZeqFT0A + multZeqFT0C);
      float centFT0A = getCentrality(ft0aInfo, multZeqFT0A);
      float centFT0C = getCentrality(ft0cInfo, multZeqFT0C);
      float centFV0A = getCentrality(fv0aInfo, multZeqFV0A);
      float centFT0CVariant1 = getCentrality(ft0cVariant1Info, multZeqFT0C);
      float centMFT = 100.5f; // missing mftNtracks in strangeness data model
      float centNGlobal = getCentrality(nGlobalInfo, collision.multNTracksGlobal());

      strangeCents(centFT0M, centFT0A,
                   centFT0C, centFV0A, centFT0CVariant1,
                   centMFT, centNGlobal);

      if (produceHistograms.value) {
        histos.fill(HIST("FT0M/Mult"), multZeqFT0A + multZeqFT0C);
        histos.fill(HIST("FT0M/percentile"), centFT0M);
        histos.fill(HIST("FT0M/percentilevsPV"), centFT0M, collision.multNTracksPVeta1());
        histos.fill(HIST("FT0M/MultvsPV"), multZeqFT0A + multZeqFT0C, collision.multNTracksPVeta1());

        histos.fill(HIST("FT0A/Mult"), multZeqFT0A);
        histos.fill(HIST("FT0A/percentile"), centFT0A);
        histos.fill(HIST("FT0A/percentilevsPV"), centFT0A, collision.multNTracksPVeta1());
        histos.fill(HIST("FT0A/MultvsPV"), multZeqFT0A, collision.multNTracksPVeta1());

        histos.fill(HIST("FT0C/Mult"), multZeqFT0C);
        histos.fill(HIST("FT0C/percentile"), centFT0C);
        histos.fill(HIST("FT0C/percentilevsPV"), centFT0C, collision.multNTracksPVeta1());
        histos.fill(HIST("FT0C/MultvsPV"), multZeqFT0C, collision.multNTracksPVeta1());
        if (collision.sel8()) {
          histos.fill(HIST("sel8FT0M/Mult"), multZeqFT0A + multZeqFT0C);
          histos.fill(HIST("sel8FT0M/percentile"), centFT0M);
          histos.fill(HIST("sel8FT0M/percentilevsPV"), centFT0M, collision.multNTracksPVeta1());
          histos.fill(HIST("sel8FT0M/MultvsPV"), multZeqFT0A + multZeqFT0C, collision.multNTracksPVeta1());

          histos.fill(HIST("sel8FT0A/Mult"), multZeqFT0A);
          histos.fill(HIST("sel8FT0A/percentile"), centFT0A);
          histos.fill(HIST("sel8FT0A/percentilevsPV"), centFT0A, collision.multNTracksPVeta1());
          histos.fill(HIST("sel8FT0A/MultvsPV"), multZeqFT0A, collision.multNTracksPVeta1());

          histos.fill(HIST("sel8FT0C/Mult"), multZeqFT0C);
          histos.fill(HIST("sel8FT0C/percentile"), centFT0C);
          histos.fill(HIST("sel8FT0C/percentilevsPV"), centFT0C, collision.multNTracksPVeta1());
          histos.fill(HIST("sel8FT0C/MultvsPV"), multZeqFT0C, collision.multNTracksPVeta1());
        }
      }
    }
  }

  // Process switches
  PROCESS_SWITCH(straCents, processRun2, "Provide Run2 calibrated centrality/multiplicity percentiles tables", false);
  PROCESS_SWITCH(straCents, processRun3, "Provide Run3 calibrated centrality/multiplicity percentiles tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<straCents>(cfgc)};
}
