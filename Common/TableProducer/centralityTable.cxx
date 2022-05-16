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
  Produces<aod::CentRun3V0Ms> centRun3V0M;
  Produces<aod::CentRun3FT0Ms> centRun3FT0M;
  Produces<aod::CentRun3FDDMs> centRun3FDDM;
  Produces<aod::CentRun3NTPVs> centRun3NTPVs;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> estRun2V0M{"estRun2V0M", -1, {"Produces centrality percentiles using V0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2SPDTrklets{"estRun2SPDtks", -1, {"Produces Run2 centrality percentiles using SPD tracklets multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2SPDClusters{"estRun2SPDcls", -1, {"Produces Run2 centrality percentiles using SPD clusters multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2CL0{"estRun2CL0", -1, {"Produces Run2 centrality percentiles using CL0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2CL1{"estRun2CL1", -1, {"Produces Run2 centrality percentiles using CL1 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun3V0M{"estRun3V0M", -1, {"Produces centrality percentiles using V0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun3FT0M{"estRun3FT0M", -1, {"Produces centrality percentiles using FT0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun3FDDM{"estRun3FDDM", -1, {"Produces centrality percentiles using FDD multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun3NTPV{"estRun3NTPV", -1, {"Produces centrality percentiles using number of tracks contributing to the PV. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
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
  struct tagRun3Calibration {
    bool mCalibrationStored = false;
    TH1* mhMultSelCalib = nullptr;
    float mMCScalePars[6] = {0.0};
    TFormula* mMCScale = nullptr;
  } Run3V0MInfo, Run3FT0MInfo, Run3FDDMInfo, Run3NTPVInfo;

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
        enable("Run3V0M", estRun3V0M);
        enable("Run3FT0M", estRun3FT0M);
        enable("Run3FDDM", estRun3FDDM);
        enable("Run3NTPV", estRun3NTPV);
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
  PROCESS_SWITCH(CentralityTable, processRun2, "Provide Run2 calibrated centrality/multiplicity percentiles tables", true);

  using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

  void processRun3(soa::Join<aod::Collisions, aod::Mults, aod::MultZeqs>::iterator const& collision, BCsWithTimestamps const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<BCsWithTimestamps>();
    if (bc.runNumber() != mRunNumber) {
      LOGF(debug, "timestamp=%llu", bc.timestamp());
      TList* callst = ccdb->getForTimeStamp<TList>(ccdbPath, bc.timestamp());

      Run3V0MInfo.mCalibrationStored = false;
      Run3FT0MInfo.mCalibrationStored = false;
      Run3FDDMInfo.mCalibrationStored = false;
      Run3NTPVInfo.mCalibrationStored = false;
      if (callst != nullptr) {
        auto getccdb = [callst](const char* ccdbhname) {
          TH1* h = (TH1*)callst->FindObject(ccdbhname);
          return h;
        };
        auto getformulaccdb = [callst](const char* ccdbhname) {
          TFormula* f = (TFormula*)callst->FindObject(ccdbhname);
          return f;
        };
        if (estRun3V0M == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run3V0MInfo.mhMultSelCalib = getccdb("hCalibZeqV0");
          Run3V0MInfo.mMCScale = getformulaccdb(TString::Format("%s-V0", genName->c_str()).Data());
          if (Run3V0MInfo.mhMultSelCalib != nullptr) {
            if (genName->length() != 0) {
              if (Run3V0MInfo.mMCScale != nullptr) {
                for (int ixpar = 0; ixpar < 6; ++ixpar) {
                  Run3V0MInfo.mMCScalePars[ixpar] = Run3V0MInfo.mMCScale->GetParameter(ixpar);
                }
              } else {
                LOGF(fatal, "MC Scale information from V0M for run %d not available", bc.runNumber());
              }
            }
            Run3V0MInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from V0M for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun3FT0M == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run3FT0MInfo.mhMultSelCalib = getccdb("hCalibZeqT0");
          if (Run3FT0MInfo.mhMultSelCalib != nullptr) {
            Run3FT0MInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from FT0 for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun3FDDM == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run3FDDMInfo.mhMultSelCalib = getccdb("hCalibZeqFDD");
          if (Run3FDDMInfo.mhMultSelCalib != nullptr) {
            Run3FDDMInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from FDD for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun3NTPV == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          Run3NTPVInfo.mhMultSelCalib = getccdb("hCalibZeqNTracks");
          if (Run3NTPVInfo.mhMultSelCalib != nullptr) {
            Run3NTPVInfo.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from NTPV for run %d corrupted", bc.runNumber());
          }
        }
        if (Run3V0MInfo.mCalibrationStored or Run3FT0MInfo.mCalibrationStored or Run3FDDMInfo.mCalibrationStored or Run3NTPVInfo.mCalibrationStored) {
          mRunNumber = bc.runNumber();
        }
      } else {
        LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
    }

    auto scaleMC = [](float x, float pars[6]) {
      return pow(((pars[0] + pars[1] * pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
    };

    if (estRun3V0M == 1) {
      float cV0M = 105.0f;
      if (Run3V0MInfo.mCalibrationStored) {
        float v0m;
        if (Run3V0MInfo.mMCScale != nullptr) {
          v0m = scaleMC(collision.multFV0M(), Run2V0MInfo.mMCScalePars);
          LOGF(debug, "Unscaled v0m: %f, scaled v0m: %f", collision.multFV0M(), v0m);
        } else {
          v0m = collision.multZeqFV0A();
        }
        cV0M = Run3V0MInfo.mhMultSelCalib->GetBinContent(Run3V0MInfo.mhMultSelCalib->FindFixBin(v0m));
      }
      LOGF(debug, "centRun2V0M=%.0f", cV0M);
      // fill centrality columns
      centRun3V0M(cV0M);
    }
    if (estRun3FT0M == 1) {
      float cFT0 = 105.0f;
      if (Run3FT0MInfo.mCalibrationStored) {
        float cft0m = collision.multZeqFT0A() + collision.multZeqFT0C();
        cFT0 = Run3FT0MInfo.mhMultSelCalib->GetBinContent(Run3FT0MInfo.mhMultSelCalib->FindFixBin(cft0m));
      }
      LOGF(debug, "centRun3FT0M=%.0f", cFT0);
      centRun3FT0M(cFT0);
    }
    if (estRun3FDDM == 1) {
      float cFDD = 105.0f;
      if (Run3FDDMInfo.mCalibrationStored) {
        float cfddm = collision.multZeqFDDA() + collision.multZeqFDDC();
        cFDD = Run3FDDMInfo.mhMultSelCalib->GetBinContent(Run3FDDMInfo.mhMultSelCalib->FindFixBin(cfddm));
      }
      LOGF(debug, "centRun3FDDM=%.0f", cFDD);
      centRun3FDDM(cFDD);
    }
    if (estRun3NTPV == 1) {
      float cNTPV = 105.0f;
      if (Run3NTPVInfo.mCalibrationStored) {
        float cntv = collision.multZeqNTracksPV();
        cNTPV = Run3NTPVInfo.mhMultSelCalib->GetBinContent(Run3NTPVInfo.mhMultSelCalib->FindFixBin(cntv));
      }
      LOGF(debug, "centRun3NTPV=%.0f", cNTPV);
      centRun3FT0M(cNTPV);
    }
  }
  PROCESS_SWITCH(CentralityTable, processRun3, "Provide Run3 calibrated centrality/multiplicity percentiles tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityTable>(cfgc)};
}
