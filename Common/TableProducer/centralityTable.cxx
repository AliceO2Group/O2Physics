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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include <CCDB/BasicCCDBManager.h>
#include <TH1F.h>
#include <TFormula.h>

using namespace o2;
using namespace o2::framework;

struct CentralityTable {
  Produces<aod::CentRun2V0Ms> centRun2V0M;
  Produces<aod::CentRun2SPDTrks> centRun2SPDTracklets;
  Produces<aod::CentRun2SPDClss> centRun2SPDClusters;
  Produces<aod::CentRun2CL0s> centRun2CL0;
  Produces<aod::CentRun2CL1s> centRun2CL1;
  Produces<aod::CentFV0As> centFV0A;
  Produces<aod::CentFT0Ms> centFT0M;
  Produces<aod::CentFDDMs> centFDDM;
  Produces<aod::CentNTPVs> centNTPV;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> estRun2V0M{"estRun2V0M", -1, {"Produces centrality percentiles using V0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2SPDTrklets{"estRun2SPDtks", -1, {"Produces Run2 centrality percentiles using SPD tracklets multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2SPDClusters{"estRun2SPDcls", -1, {"Produces Run2 centrality percentiles using SPD clusters multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2CL0{"estRun2CL0", -1, {"Produces Run2 centrality percentiles using CL0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2CL1{"estRun2CL1", -1, {"Produces Run2 centrality percentiles using CL1 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estFV0A{"estFV0A", -1, {"Produces centrality percentiles using FV0A multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estFT0M{"estFT0M", -1, {"Produces centrality percentiles using FT0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estFDDM{"estFDDM", -1, {"Produces centrality percentiles using FDD multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estNTPV{"estNTPV", -1, {"Produces centrality percentiles using number of tracks contributing to the PV. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
  Configurable<std::string> ccdbPath{"ccdbpath", "Centrality/Estimators", "The CCDB path for centrality/multiplicity information"};
  Configurable<std::string> genName{"genname", "", "Genearator name: HIJING, PYTHIA8, ... Default: \"\""};

  int mRunNumber;
  struct tagRun2V0MCalibration {
    bool mCalibrationStored = false;
    TFormula* mMCScale = nullptr;
    float mMCScalePars[6] = {0.0};
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhVtxAmpCorrV0C = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } Run2V0MInfo;
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
    calibrationInfo(std::string name)
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
  calibrationInfo FDDMInfo = calibrationInfo("FDD");
  calibrationInfo NTPVInfo = calibrationInfo("NTracksPV");

  void init(InitContext& context)
  {
    /* Checking the tables which are requested in the workflow and enabling them */
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec device : workflows.devices) {
      for (auto input : device.inputs) {
        auto enable = [&input](const std::string detector, Configurable<int>& flag) {
          const std::string table = "Cent" + detector + "s";
          if (input.matcher.binding == table) {
            if (flag < 0) {
              flag.value = 1;
              LOGF(info, "Auto-enabling table: %s", table.c_str());
            } else if (flag > 0) {
              flag.value = 1;
              LOGF(info, "Table %s already enabled", table.c_str());
            } else {
              LOGF(info, "Table %s disabled", table.c_str());
            }
          }
        };
        enable("Run2V0M", estRun2V0M);
        enable("Run2SPDTrk", estRun2SPDTrklets);
        enable("Run2SPDCls", estRun2SPDClusters);
        enable("Run2CL0", estRun2CL0);
        enable("Run2CL1", estRun2CL1);
        enable("FV0A", estFV0A);
        enable("FT0M", estFT0M);
        enable("FDDM", estFDDM);
        enable("NTPV", estNTPV);
      }
    }
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
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
      Run2SPDTksInfo.mCalibrationStored = false;
      Run2SPDClsInfo.mCalibrationStored = false;
      Run2CL0Info.mCalibrationStored = false;
      Run2CL1Info.mCalibrationStored = false;
      if (callst != nullptr) {
        auto getccdb = [callst](const char* ccdbhname) {
          TH1* h = (TH1*)callst->FindObject(ccdbhname);
          return h;
        };
        auto getformulaccdb = [callst](const char* ccdbhname) {
          TFormula* f = (TFormula*)callst->FindObject(ccdbhname);
          return f;
        };
        if (estRun2V0M == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2V0MInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          Run2V0MInfo.mhVtxAmpCorrV0C = getccdb("hVtx_fAmplitude_V0C_Normalized");
          Run2V0MInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0M");
          Run2V0MInfo.mMCScale = getformulaccdb(TString::Format("%s-V0M", genName->c_str()).Data());
          if ((Run2V0MInfo.mhVtxAmpCorrV0A != nullptr) and (Run2V0MInfo.mhVtxAmpCorrV0C != nullptr) and (Run2V0MInfo.mhMultSelCalib != nullptr)) {
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
        if (estRun2SPDTrklets == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2SPDTksInfo.mhVtxAmpCorr = getccdb("hVtx_fnTracklets_Normalized");
          Run2SPDTksInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDTracklets");
          if ((Run2SPDTksInfo.mhVtxAmpCorr != nullptr) and (Run2SPDTksInfo.mhMultSelCalib != nullptr)) {
            Run2SPDTksInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from SPD tracklets for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2SPDClusters == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2SPDClsInfo.mhVtxAmpCorrCL0 = getccdb("hVtx_fnSPDClusters0_Normalized");
          Run2SPDClsInfo.mhVtxAmpCorrCL1 = getccdb("hVtx_fnSPDClusters1_Normalized");
          Run2SPDClsInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDClusters");
          if ((Run2SPDClsInfo.mhVtxAmpCorrCL0 != nullptr) and (Run2SPDClsInfo.mhVtxAmpCorrCL1 != nullptr) and (Run2SPDClsInfo.mhMultSelCalib != nullptr)) {
            Run2SPDClsInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from SPD clusters for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2CL0 == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2CL0Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters0_Normalized");
          Run2CL0Info.mhMultSelCalib = getccdb("hMultSelCalib_CL0");
          if ((Run2CL0Info.mhVtxAmpCorr != nullptr) and (Run2CL0Info.mhMultSelCalib != nullptr)) {
            Run2CL0Info.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from CL0 multiplicity for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2CL1 == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run2CL1Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters1_Normalized");
          Run2CL1Info.mhMultSelCalib = getccdb("hMultSelCalib_CL1");
          if ((Run2CL1Info.mhVtxAmpCorr != nullptr) and (Run2CL1Info.mhMultSelCalib != nullptr)) {
            Run2CL1Info.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from CL1 multiplicity for run %d corrupted", bc.runNumber());
          }
        }
        if (Run2V0MInfo.mCalibrationStored or Run2SPDTksInfo.mCalibrationStored or Run2SPDClsInfo.mCalibrationStored or Run2CL0Info.mCalibrationStored or Run2CL1Info.mCalibrationStored) {
          mRunNumber = bc.runNumber();
        }
      } else {
        LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
    }

    auto scaleMC = [](float x, float pars[6]) {
      return pow(((pars[0] + pars[1] * pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
    };

    if (estRun2V0M == 1) {
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
    if (estRun2SPDTrklets == 1) {
      float cSPD = 105.0f;
      if (Run2SPDTksInfo.mCalibrationStored) {
        float spdm = collision.multTracklets() * Run2SPDTksInfo.mhVtxAmpCorr->GetBinContent(Run2SPDTksInfo.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cSPD = Run2SPDTksInfo.mhMultSelCalib->GetBinContent(Run2SPDTksInfo.mhMultSelCalib->FindFixBin(spdm));
      }
      LOGF(debug, "centSPDTracklets=%.0f", cSPD);
      centRun2SPDTracklets(cSPD);
    }
    if (estRun2SPDClusters == 1) {
      float cSPD = 105.0f;
      if (Run2SPDClsInfo.mCalibrationStored) {
        float spdm = bc.spdClustersL0() * Run2SPDClsInfo.mhVtxAmpCorrCL0->GetBinContent(Run2SPDClsInfo.mhVtxAmpCorrCL0->FindFixBin(collision.posZ())) +
                     bc.spdClustersL1() * Run2SPDClsInfo.mhVtxAmpCorrCL1->GetBinContent(Run2SPDClsInfo.mhVtxAmpCorrCL1->FindFixBin(collision.posZ()));
        cSPD = Run2SPDClsInfo.mhMultSelCalib->GetBinContent(Run2SPDClsInfo.mhMultSelCalib->FindFixBin(spdm));
      }
      LOGF(debug, "centSPDClusters=%.0f", cSPD);
      centRun2SPDClusters(cSPD);
    }
    if (estRun2CL0 == 1) {
      float cCL0 = 105.0f;
      if (Run2CL0Info.mCalibrationStored) {
        float cl0m = bc.spdClustersL0() * Run2CL0Info.mhVtxAmpCorr->GetBinContent(Run2CL0Info.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cCL0 = Run2CL0Info.mhMultSelCalib->GetBinContent(Run2CL0Info.mhMultSelCalib->FindFixBin(cl0m));
      }
      LOGF(debug, "centCL0=%.0f", cCL0);
      centRun2CL0(cCL0);
    }
    if (estRun2CL1 == 1) {
      float cCL1 = 105.0f;
      if (Run2CL1Info.mCalibrationStored) {
        float cl1m = bc.spdClustersL1() * Run2CL1Info.mhVtxAmpCorr->GetBinContent(Run2CL1Info.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cCL1 = Run2CL1Info.mhMultSelCalib->GetBinContent(Run2CL1Info.mhMultSelCalib->FindFixBin(cl1m));
      }
      LOGF(debug, "centCL1=%.0f", cCL1);
      centRun2CL1(cCL1);
    }
  }
  PROCESS_SWITCH(CentralityTable, processRun2, "Provide Run2 calibrated centrality/multiplicity percentiles tables", false);

  using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

  void processRun3(soa::Join<aod::Collisions, aod::Mults, aod::MultZeqs>::iterator const& collision, BCsWithTimestamps const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<BCsWithTimestamps>();
    if (bc.runNumber() != mRunNumber) {
      LOGF(info, "timestamp=%llu, run number=%d", bc.timestamp(), bc.runNumber());
      TList* callst = ccdb->getForTimeStamp<TList>(ccdbPath, bc.timestamp());

      FV0AInfo.mCalibrationStored = false;
      FT0MInfo.mCalibrationStored = false;
      FDDMInfo.mCalibrationStored = false;
      NTPVInfo.mCalibrationStored = false;
      if (callst != nullptr) {
        LOGF(info, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
        auto getccdb = [callst, bc](struct calibrationInfo& estimator, const Configurable<std::string> generatorName) { // TODO: to consider the name inside the estimator structure
          estimator.mhMultSelCalib = (TH1*)callst->FindObject(TString::Format("hCalibZeq%s", estimator.name.c_str()).Data());
          estimator.mMCScale = (TFormula*)callst->FindObject(TString::Format("%s-%s", generatorName->c_str(), estimator.name.c_str()).Data());
          if (estimator.mhMultSelCalib != nullptr) {
            if (generatorName->length() != 0) {
              if (estimator.mMCScale != nullptr) {
                for (int ixpar = 0; ixpar < 6; ++ixpar) {
                  estimator.mMCScalePars[ixpar] = estimator.mMCScale->GetParameter(ixpar);
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
        if (estFV0A == 1) {
          getccdb(FV0AInfo, genName);
        }
        if (estFT0M == 1) {
          getccdb(FT0MInfo, genName);
        }
        if (estFDDM == 1) {
          getccdb(FDDMInfo, genName);
        }
        if (estNTPV == 1) {
          getccdb(NTPVInfo, genName);
        }
        mRunNumber = bc.runNumber();
      } else {
        LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
    }

    auto populateTable = [](auto& table, struct calibrationInfo& estimator, float multiplicity) {
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
      }
      LOGF(debug, "%s centrality/multiplicity percentile = %.0f for a zvtx eq %s value %.0f", estimator.name.c_str(), percentile, estimator.name.c_str(), scaledMultiplicity);
      table(percentile);
    };

    if (estFV0A == 1) {
      populateTable(centFV0A, FV0AInfo, collision.multZeqFV0A());
    }
    if (estFT0M == 1) {
      populateTable(centFT0M, FT0MInfo, collision.multZeqFT0A() + collision.multZeqFT0C());
    }
    if (estFDDM == 1) {
      populateTable(centFDDM, FDDMInfo, collision.multZeqFDDA() + collision.multZeqFDDC());
    }
    if (estNTPV == 1) {
      populateTable(centNTPV, NTPVInfo, collision.multZeqNTracksPV());
    }
  }
  PROCESS_SWITCH(CentralityTable, processRun3, "Provide Run3 calibrated centrality/multiplicity percentiles tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityTable>(cfgc)};
}
