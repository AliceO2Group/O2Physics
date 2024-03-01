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

/// \file centrality.cxx
/// \brief Task to produce the centrality tables associated to each of the required centrality estimators

#include <CCDB/BasicCCDBManager.h>
#include <TH1F.h>
#include <TFormula.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "TableHelper.h"

using namespace o2;
using namespace o2::framework;

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
static constexpr int kCentFDDMs = 10;
static constexpr int kCentNTPVs = 11;
static constexpr int nTables = 12;
static constexpr int nParameters = 1;
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
                                                 "CentFDDMs",
                                                 "CentNTPVs"};
static const std::vector<std::string> parameterNames{"Enable"};
static const int defaultParameters[nTables][nParameters]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};

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
  Produces<aod::CentFDDMs> centFDDM;
  Produces<aod::CentNTPVs> centNTPV;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], nTables, nParameters, tableNames, parameterNames},
                                                "Produce tables depending on needs. Values different than -1 override the automatic setup: the corresponding table can be set off (0) or on (1)"};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
  Configurable<std::string> ccdbPath{"ccdbpath", "Centrality/Estimators", "The CCDB path for centrality/multiplicity information"};
  Configurable<std::string> genName{"genname", "", "Genearator name: HIJING, PYTHIA8, ... Default: \"\""};
  Configurable<bool> doNotCrashOnNull{"doNotCrashOnNull", false, {"Option to not crash on null and instead fill required tables with dummy info"}};
  Configurable<bool> embedINELgtZEROselection{"embedINELgtZEROselection", false, {"Option to do percentile 100.5 if not INELgtZERO"}};

  int mRunNumber;
  struct tagRun2V0MCalibration {
    bool mCalibrationStored = false;
    TFormula* mMCScale = nullptr;
    float mMCScalePars[6] = {0.0};
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhVtxAmpCorrV0C = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2V0MInfo;
  struct tagRun2V0ACalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2V0AInfo;
  struct tagRun2SPDTrackletsCalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2SPDTksInfo;
  struct tagRun2SPDClustersCalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorrCL0 = nullptr;
    TH1* mhVtxAmpCorrCL1 = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2SPDClsInfo;
  struct tagRun2CL0Calibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2CL0Info;
  struct tagRun2CL1Calibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2CL1Info;
  struct calibrationInfo {
    std::string name = "";
    bool mCalibrationStored = false;
    TH1* mhMultSelCalib = nullptr;
    float mMCScalePars[6] = {0.0};
    TFormula* mMCScale = nullptr;
    explicit calibrationInfo(std::string name)
      : name(name),
        mCalibrationStored(false),
        mhMultSelCalib(nullptr),
        mMCScalePars{0.0},
        mMCScale(nullptr)
    {
    }
  };
  calibrationInfo FV0AInfo = calibrationInfo("FV0");
  calibrationInfo FT0MInfo = calibrationInfo("FT0");
  calibrationInfo FT0AInfo = calibrationInfo("FT0A");
  calibrationInfo FT0CInfo = calibrationInfo("FT0C");
  calibrationInfo FDDMInfo = calibrationInfo("FDD");
  calibrationInfo NTPVInfo = calibrationInfo("NTracksPV");
  std::vector<int> mEnabledTables; // Vector of enabled tables
  std::array<bool, nTables> isTableEnabled;

  void init(InitContext& context)
  {
    LOG(info) << "Initializing centrality table producer";
    if (doprocessRun2 == false && doprocessRun3 == false) {
      LOGF(fatal, "Neither processRun2 nor processRun3 enabled. Please choose one.");
    }
    if (doprocessRun2 == true && doprocessRun3 == true) {
      LOGF(fatal, "Cannot enable processRun2 and processRun3 at the same time. Please choose one.");
    }

    /* Checking the tables which are requested in the workflow and enabling them */
    for (int i = 0; i < nTables; i++) {
      int f = enabledTables->get(tableNames[i].c_str(), "Enable");
      enableFlagIfTableRequired(context, tableNames[i], f);
      if (f == 1) {
        if (tableNames[i].find("Run2") != std::string::npos) {
          if (doprocessRun3) {
            LOGF(fatal, "Cannot enable Run2 tables in Run3 mode. Please check and disable them.");
          }
        } else {
          if (doprocessRun2) {
            LOGF(fatal, "Cannot enable Run3 tables in Run2 mode. Please check and disable them.");
          }
        }
        isTableEnabled[i] = true;
        mEnabledTables.push_back(i);
      }
    }

    if (mEnabledTables.size() == 0) {
      LOGF(fatal, "No table enabled. Please enable at least one table.");
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    mRunNumber = 0;
  }

  using BCsWithTimestampsAndRun2Infos = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

  void processRun2(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, BCsWithTimestampsAndRun2Infos const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<BCsWithTimestampsAndRun2Infos>();
    if (bc.runNumber() != mRunNumber) {
      LOGF(debug, "timestamp=%llu", bc.timestamp());
      TList* callst = ccdb->getForTimeStamp<TList>(ccdbPath, bc.timestamp());

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
          Run2V0MInfo.mMCScale = getformulaccdb(TString::Format("%s-V0M", genName->c_str()).Data());
          if ((Run2V0MInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0MInfo.mhVtxAmpCorrV0C != nullptr) && (Run2V0MInfo.mhMultSelCalib != nullptr)) {
            if (genName->length() != 0) {
              if (Run2V0MInfo.mMCScale != nullptr) {
                for (int ixpar = 0; ixpar < 6; ++ixpar) {
                  Run2V0MInfo.mMCScalePars[ixpar] = Run2V0MInfo.mMCScale->GetParameter(ixpar);
                }
              } else {
                LOGF(fatal, "MC Scale information from V0M for run %d not available", bc.runNumber());
              }
            }
            Run2V0MInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from V0M for run %d corrupted", bc.runNumber());
          }
        }
        if (isTableEnabled[kCentRun2V0As]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2V0AInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          Run2V0AInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0A");
          if ((Run2V0AInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0AInfo.mhMultSelCalib != nullptr)) {
            Run2V0AInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from V0A for run %d corrupted", bc.runNumber());
          }
        }
        if (isTableEnabled[kCentRun2SPDTrks]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2SPDTksInfo.mhVtxAmpCorr = getccdb("hVtx_fnTracklets_Normalized");
          Run2SPDTksInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDTracklets");
          if ((Run2SPDTksInfo.mhVtxAmpCorr != nullptr) && (Run2SPDTksInfo.mhMultSelCalib != nullptr)) {
            Run2SPDTksInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from SPD tracklets for run %d corrupted", bc.runNumber());
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
            LOGF(fatal, "Calibration information from SPD clusters for run %d corrupted", bc.runNumber());
          }
        }
        if (isTableEnabled[kCentRun2CL0s]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2CL0Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters0_Normalized");
          Run2CL0Info.mhMultSelCalib = getccdb("hMultSelCalib_CL0");
          if ((Run2CL0Info.mhVtxAmpCorr != nullptr) && (Run2CL0Info.mhMultSelCalib != nullptr)) {
            Run2CL0Info.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from CL0 multiplicity for run %d corrupted", bc.runNumber());
          }
        }
        if (isTableEnabled[kCentRun2CL1s]) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2CL1Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters1_Normalized");
          Run2CL1Info.mhMultSelCalib = getccdb("hMultSelCalib_CL1");
          if ((Run2CL1Info.mhVtxAmpCorr != nullptr) && (Run2CL1Info.mhMultSelCalib != nullptr)) {
            Run2CL1Info.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from CL1 multiplicity for run %d corrupted", bc.runNumber());
          }
        }
        if (Run2V0MInfo.mCalibrationStored || Run2V0AInfo.mCalibrationStored || Run2SPDTksInfo.mCalibrationStored || Run2SPDClsInfo.mCalibrationStored || Run2CL0Info.mCalibrationStored || Run2CL1Info.mCalibrationStored) {
          mRunNumber = bc.runNumber();
        }
      } else {
        if (!doNotCrashOnNull) { // default behaviour: crash
          LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        } else { // only if asked: continue filling with non-valid values (105)
          LOGF(info, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", bc.runNumber(), bc.timestamp());
          mRunNumber = bc.runNumber();
        }
      }
    }

    auto scaleMC = [](float x, float pars[6]) {
      return pow(((pars[0] + pars[1] * pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
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
  PROCESS_SWITCH(CentralityTable, processRun2, "Provide Run2 calibrated centrality/multiplicity percentiles tables", true);

  using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

  void processRun3(soa::Join<aod::Collisions, aod::Mults, aod::MultZeqs> const& collisions, BCsWithTimestamps const&)
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
        case kCentFDDMs:
          centFDDM.reserve(collisions.size());
          break;
        case kCentNTPVs:
          centNTPV.reserve(collisions.size());
          break;
        default:
          LOGF(fatal, "Table %d not supported in Run3", table);
          break;
      }
    }

    for (auto const& collision : collisions) {
      /* check the previous run number */
      auto bc = collision.bc_as<BCsWithTimestamps>();
      if (bc.runNumber() != mRunNumber) {
        LOGF(info, "timestamp=%llu, run number=%d", bc.timestamp(), bc.runNumber());
        TList* callst = ccdb->getForTimeStamp<TList>(ccdbPath, bc.timestamp());

        FV0AInfo.mCalibrationStored = false;
        FT0MInfo.mCalibrationStored = false;
        FT0AInfo.mCalibrationStored = false;
        FT0CInfo.mCalibrationStored = false;
        FDDMInfo.mCalibrationStored = false;
        NTPVInfo.mCalibrationStored = false;
        if (callst != nullptr) {
          LOGF(info, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          auto getccdb = [callst, bc](struct calibrationInfo& estimator, const Configurable<std::string> generatorName) { // TODO: to consider the name inside the estimator structure
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
            } else {
              LOGF(error, "Calibration information from %s for run %d not available", estimator.name.c_str(), bc.runNumber());
            }
          };

          for (auto const& table : mEnabledTables) {
            switch (table) {
              case kCentFV0As:
                getccdb(FV0AInfo, genName);
                break;
              case kCentFT0Ms:
                getccdb(FT0MInfo, genName);
                break;
              case kCentFT0As:
                getccdb(FT0AInfo, genName);
                break;
              case kCentFT0Cs:
                getccdb(FT0CInfo, genName);
                break;
              case kCentFDDMs:
                getccdb(FDDMInfo, genName);
                break;
              case kCentNTPVs:
                getccdb(NTPVInfo, genName);
                break;
              default:
                LOGF(fatal, "Table %d not supported in Run3", table);
                break;
            }
          }
          mRunNumber = bc.runNumber();
        } else {
          if (!doNotCrashOnNull) { // default behaviour: crash
            LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
          } else { // only if asked: continue filling with non-valid values (105)
            LOGF(info, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu, will fill tables with dummy values", bc.runNumber(), bc.timestamp());
            mRunNumber = bc.runNumber();
          }
        }
      }

      auto populateTable = [](auto& table, struct calibrationInfo& estimator, float multiplicity, bool assignOutOfRange) {
        auto scaleMC = [](float x, float pars[6]) {
          return pow(((pars[0] + pars[1] * pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
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
      };

      for (auto const& table : mEnabledTables) {
        switch (table) {
          case kCentFV0As:
            populateTable(centFV0A, FV0AInfo, collision.multZeqFV0A(), collision.multNTracksPVeta1() < 1 && embedINELgtZEROselection);
            break;
          case kCentFT0Ms:
            populateTable(centFT0M, FT0MInfo, collision.multZeqFT0A() + collision.multZeqFT0C(), collision.multNTracksPVeta1() < 1 && embedINELgtZEROselection);
            break;
          case kCentFT0As:
            populateTable(centFT0A, FT0AInfo, collision.multZeqFT0A(), collision.multNTracksPVeta1() < 1 && embedINELgtZEROselection);
            break;
          case kCentFT0Cs:
            populateTable(centFT0C, FT0CInfo, collision.multZeqFT0C(), collision.multNTracksPVeta1() < 1 && embedINELgtZEROselection);
            break;
          case kCentFDDMs:
            populateTable(centFDDM, FDDMInfo, collision.multZeqFDDA() + collision.multZeqFDDC(), collision.multNTracksPVeta1() < 1 && embedINELgtZEROselection);
            break;
          case kCentNTPVs:
            populateTable(centNTPV, NTPVInfo, collision.multZeqNTracksPV(), collision.multNTracksPVeta1() < 1 && embedINELgtZEROselection);
            break;
          default:
            LOGF(fatal, "Table %d not supported in Run3", table);
            break;
        }
      }
    }
  }
  PROCESS_SWITCH(CentralityTable, processRun3, "Provide Run3 calibrated centrality/multiplicity percentiles tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityTable>(cfgc)};
}
