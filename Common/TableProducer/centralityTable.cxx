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

using namespace o2;
using namespace o2::framework;

struct CentralityTable {
  Produces<aod::CentV0Ms> centVOM;
  Produces<aod::CentRun2SPDs> centRun2SPD;
  Produces<aod::CentRun2CL0s> centRun2CL1;
  Produces<aod::CentRun2CL1s> centRun2CL2;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> estV0M{"estV0M", -1, {"Produces centrality percentiles using V0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2SPD{"estRun2SPD", -1, {"Produces Run2 centrality percentiles using SPD tracklets multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2CL0{"estRun2CL0", -1, {"Produces Run2 centrality percentiles using CL0 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};
  Configurable<int> estRun2CL1{"estRun2CL1", -1, {"Produces Run2 centrality percentiles using CL1 multiplicity. -1: auto, 0: don't, 1: yes. Default: auto (-1)"}};

  int mRunNumber;
  struct tagV0MCalibration {
    bool mV0MCalibrationStored = false;
    TH1* mhVtxAmpCorrV0A = nullptr;
    TH1* mhVtxAmpCorrV0C = nullptr;
    TH1* mhMultSelCalibV0M = nullptr;
  } V0MInfo;
  struct tagSPDTrackletsCalibration {
    bool mSPDCalibrationStored = false;
    TH1* mhVtxAmpCorrSPD = nullptr;
    TH1* mhMultSelCalibSPD = nullptr;
  } SPDInfo;

  void
    init(InitContext& context)
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
        enable("V0M", estV0M);
        enable("Run2SPD", estRun2SPD);
        enable("Run2CL0", estRun2CL0);
        enable("Run2CL1", estRun2CL1);
      }
    }
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    mRunNumber = 0;
  }

  void process(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, aod::BCsWithTimestamps const&, aod::Tracks const& tracks)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (bc.runNumber() != mRunNumber) {
      LOGF(debug, "timestamp=%llu", bc.timestamp());
      TList* callst = ccdb->getForTimeStamp<TList>("Centrality/Estimators", bc.timestamp());

      V0MInfo.mV0MCalibrationStored = false;
      SPDInfo.mSPDCalibrationStored = false;
      if (callst != nullptr) {
        auto getccdb = [callst](const char* ccdbhname) {
          TH1* h = (TH1*)callst->FindObject(ccdbhname);
          return h;
        };
        if (estV0M == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          V0MInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          V0MInfo.mhVtxAmpCorrV0C = getccdb("hVtx_fAmplitude_V0C_Normalized");
          V0MInfo.mhMultSelCalibV0M = getccdb("hMultSelCalib_V0M");
          if ((V0MInfo.mhVtxAmpCorrV0A != nullptr) and (V0MInfo.mhVtxAmpCorrV0C != nullptr) and (V0MInfo.mhMultSelCalibV0M != nullptr)) {
            V0MInfo.mV0MCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from V0M for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2SPD == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          SPDInfo.mhVtxAmpCorrSPD = getccdb("hVtx_fnTracklets_Normalized");
          SPDInfo.mhMultSelCalibSPD = getccdb("hMultSelCalib_SPDTracklets");
          if ((SPDInfo.mhVtxAmpCorrSPD != nullptr) and (SPDInfo.mhMultSelCalibSPD != nullptr)) {
            SPDInfo.mSPDCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from SPD tracklets for run %d corrupted", bc.runNumber());
          }
        }
        if (estRun2CL0 == 1) {
          LOGF(fatal, "Run2 calibration information estimated from CL0 still not available");
        }
        if (estRun2CL1 == 1) {
          LOGF(fatal, "Run2 calibration information estimated from CL1 still not available");
        }
        if (V0MInfo.mV0MCalibrationStored or SPDInfo.mSPDCalibrationStored) {
          mRunNumber = bc.runNumber();
        }
      } else {
        LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
    }
    if (estV0M == 1) {
      float cV0M = 105.0f;
      if (V0MInfo.mV0MCalibrationStored) {
        float v0m = collision.multV0A() * V0MInfo.mhVtxAmpCorrV0A->GetBinContent(V0MInfo.mhVtxAmpCorrV0A->FindFixBin(collision.posZ())) +
                    collision.multV0C() * V0MInfo.mhVtxAmpCorrV0C->GetBinContent(V0MInfo.mhVtxAmpCorrV0C->FindFixBin(collision.posZ()));
        cV0M = V0MInfo.mhMultSelCalibV0M->GetBinContent(V0MInfo.mhMultSelCalibV0M->FindFixBin(v0m));
      }
      LOGF(debug, "centV0M=%.0f", cV0M);
      // fill centrality columns
      centVOM(cV0M);
    }
    if (estRun2SPD == 1) {
      float cSPD = 105.0f;
      if (SPDInfo.mSPDCalibrationStored) {
        float spdm = collision.multTracklets() * SPDInfo.mhVtxAmpCorrSPD->GetBinContent(SPDInfo.mhVtxAmpCorrSPD->FindFixBin(collision.posZ()));
        cSPD = SPDInfo.mhMultSelCalibSPD->GetBinContent(SPDInfo.mhMultSelCalibSPD->FindFixBin(spdm));
      }
      LOGF(debug, "centSPD=%.0f", cSPD);
      centRun2SPD(cSPD);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityTable>(cfgc)};
}
