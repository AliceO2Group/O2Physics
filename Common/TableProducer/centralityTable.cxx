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
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> estRun2V0M{"estRun2V0M", -1, {"Produces centrality percentiles using V0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2SPDTrklets{"estRun2SPDtks", -1, {"Produces Run2 centrality percentiles using SPD tracklets multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2SPDClusters{"estRun2SPDcls", -1, {"Produces Run2 centrality percentiles using SPD clusters multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2CL0{"estRun2CL0", -1, {"Produces Run2 centrality percentiles using CL0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2CL1{"estRun2CL1", -1, {"Produces Run2 centrality percentiles using CL1 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
  Configurable<std::string> ccdbPath{"ccdbpath", "Centrality/Estimators", "The CCDB path for centrality/multiplicity information"};
  Configurable<std::string> genName{"genname", "", "Genearator name: HIJING, PYTHIA8, ... Default: \"\""};

  int mRunNumber;
  struct tagV0MCalibration {
    bool mCalibrationStored = false;
    TFormula* mMCScale = nullptr;
    float mMCScalePars[6] = {0.0};
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhVtxAmpCorrV0C = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } V0MInfo;
  struct tagSPDTrackletsCalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } SPDTksInfo;
  struct tagSPDClustersCalibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorrCL0 = nullptr;
    TH1* mhVtxAmpCorrCL1 = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } SPDClsInfo;
  struct tagCL0Calibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } CL0Info;
  struct tagCL1Calibration {
    bool mCalibrationStored = false;
    TH1* mhVtxAmpCorr = nullptr;
    TH1* mhMultSelCalib = nullptr;
  } CL1Info;

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
      }
    }
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    mRunNumber = 0;
  }

  using BCsWithTimestampsAndRun2Infos = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

  void process(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, BCsWithTimestampsAndRun2Infos const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<BCsWithTimestampsAndRun2Infos>();
    if (bc.runNumber() != mRunNumber) {
      LOGF(debug, "timestamp=%llu", bc.timestamp());
      TList* callst = ccdb->getForTimeStamp<TList>(ccdbPath, bc.timestamp());

      V0MInfo.mCalibrationStored = false;
      SPDTksInfo.mCalibrationStored = false;
      SPDClsInfo.mCalibrationStored = false;
      CL0Info.mCalibrationStored = false;
      CL1Info.mCalibrationStored = false;
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
          V0MInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          V0MInfo.mhVtxAmpCorrV0C = getccdb("hVtx_fAmplitude_V0C_Normalized");
          V0MInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0M");
          V0MInfo.mMCScale = getformulaccdb(TString::Format("%s-V0M", genName->c_str()).Data());
          if ((V0MInfo.mhVtxAmpCorrV0A != nullptr) and (V0MInfo.mhVtxAmpCorrV0C != nullptr) and (V0MInfo.mhMultSelCalib != nullptr)) {
            if (genName->length() != 0) {
              if (V0MInfo.mMCScale != nullptr) {
                for (int ixpar = 0; ixpar < 6; ++ixpar) {
                  V0MInfo.mMCScalePars[ixpar] = V0MInfo.mMCScale->GetParameter(ixpar);
                }
              } else {
                LOGF(fatal, "MC Scale information from V0M for run %d not available", bc.runNumber());
              }
            }
            V0MInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from V0M for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2SPDTrklets == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          SPDTksInfo.mhVtxAmpCorr = getccdb("hVtx_fnTracklets_Normalized");
          SPDTksInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDTracklets");
          if ((SPDTksInfo.mhVtxAmpCorr != nullptr) and (SPDTksInfo.mhMultSelCalib != nullptr)) {
            SPDTksInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from SPD tracklets for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2SPDClusters == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          SPDClsInfo.mhVtxAmpCorrCL0 = getccdb("hVtx_fnSPDClusters0_Normalized");
          SPDClsInfo.mhVtxAmpCorrCL1 = getccdb("hVtx_fnSPDClusters1_Normalized");
          SPDClsInfo.mhMultSelCalib = getccdb("hMultSelCalib_SPDClusters");
          if ((SPDClsInfo.mhVtxAmpCorrCL0 != nullptr) and (SPDClsInfo.mhVtxAmpCorrCL1 != nullptr) and (SPDClsInfo.mhMultSelCalib != nullptr)) {
            SPDClsInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from SPD clusters for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2CL0 == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          CL0Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters0_Normalized");
          CL0Info.mhMultSelCalib = getccdb("hMultSelCalib_CL0");
          if ((CL0Info.mhVtxAmpCorr != nullptr) and (CL0Info.mhMultSelCalib != nullptr)) {
            CL0Info.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from CL0 multiplicity for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2CL1 == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          CL1Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters1_Normalized");
          CL1Info.mhMultSelCalib = getccdb("hMultSelCalib_CL1");
          if ((CL1Info.mhVtxAmpCorr != nullptr) and (CL1Info.mhMultSelCalib != nullptr)) {
            CL1Info.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from CL1 multiplicity for run %d corrupted", bc.runNumber());
          }
        }
        if (V0MInfo.mCalibrationStored or SPDTksInfo.mCalibrationStored or SPDClsInfo.mCalibrationStored or CL0Info.mCalibrationStored or CL1Info.mCalibrationStored) {
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
      if (V0MInfo.mCalibrationStored) {
        float v0m;
        if (V0MInfo.mMCScale != nullptr) {
          v0m = scaleMC(collision.multFV0M(), V0MInfo.mMCScalePars);
          LOGF(debug, "Unscaled v0m: %f, scaled v0m: %f", collision.multFV0M(), v0m);
        } else {
          v0m = collision.multFV0A() * V0MInfo.mhVtxAmpCorrV0A->GetBinContent(V0MInfo.mhVtxAmpCorrV0A->FindFixBin(collision.posZ())) +
                collision.multFV0C() * V0MInfo.mhVtxAmpCorrV0C->GetBinContent(V0MInfo.mhVtxAmpCorrV0C->FindFixBin(collision.posZ()));
        }
        cV0M = V0MInfo.mhMultSelCalib->GetBinContent(V0MInfo.mhMultSelCalib->FindFixBin(v0m));
      }
      LOGF(debug, "centRun2V0M=%.0f", cV0M);
      // fill centrality columns
      centRun2V0M(cV0M);
    }
    if (estRun2SPDTrklets == 1) {
      float cSPD = 105.0f;
      if (SPDTksInfo.mCalibrationStored) {
        float spdm = collision.multTracklets() * SPDTksInfo.mhVtxAmpCorr->GetBinContent(SPDTksInfo.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cSPD = SPDTksInfo.mhMultSelCalib->GetBinContent(SPDTksInfo.mhMultSelCalib->FindFixBin(spdm));
      }
      LOGF(debug, "centSPDTracklets=%.0f", cSPD);
      centRun2SPDTracklets(cSPD);
    }
    if (estRun2SPDClusters == 1) {
      float cSPD = 105.0f;
      if (SPDClsInfo.mCalibrationStored) {
        float spdm = bc.spdClustersL0() * SPDClsInfo.mhVtxAmpCorrCL0->GetBinContent(SPDClsInfo.mhVtxAmpCorrCL0->FindFixBin(collision.posZ())) +
                     bc.spdClustersL1() * SPDClsInfo.mhVtxAmpCorrCL1->GetBinContent(SPDClsInfo.mhVtxAmpCorrCL1->FindFixBin(collision.posZ()));
        cSPD = SPDClsInfo.mhMultSelCalib->GetBinContent(SPDClsInfo.mhMultSelCalib->FindFixBin(spdm));
      }
      LOGF(debug, "centSPDClusters=%.0f", cSPD);
      centRun2SPDClusters(cSPD);
    }
    if (estRun2CL0 == 1) {
      float cCL0 = 105.0f;
      if (CL0Info.mCalibrationStored) {
        float cl0m = bc.spdClustersL0() * CL0Info.mhVtxAmpCorr->GetBinContent(CL0Info.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cCL0 = CL0Info.mhMultSelCalib->GetBinContent(CL0Info.mhMultSelCalib->FindFixBin(cl0m));
      }
      LOGF(debug, "centCL0=%.0f", cCL0);
      centRun2CL0(cCL0);
    }
    if (estRun2CL1 == 1) {
      float cCL1 = 105.0f;
      if (CL1Info.mCalibrationStored) {
        float cl1m = bc.spdClustersL1() * CL1Info.mhVtxAmpCorr->GetBinContent(CL1Info.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        cCL1 = CL1Info.mhMultSelCalib->GetBinContent(CL1Info.mhMultSelCalib->FindFixBin(cl1m));
      }
      LOGF(debug, "centCL1=%.0f", cCL1);
      centRun2CL1(cCL1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityTable>(cfgc)};
}
