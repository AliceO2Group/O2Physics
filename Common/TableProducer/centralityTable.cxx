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
  bool mV0MCalibrationStored;
  TH1* mhVtxAmpCorrV0A;
  TH1* mhVtxAmpCorrV0C;
  TH1* mhMultSelCalibV0M;

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
    mV0MCalibrationStored = false;
    mhVtxAmpCorrV0A = nullptr;
    mhVtxAmpCorrV0C = nullptr;
    mhMultSelCalibV0M = nullptr;
  }

  void process(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, aod::BCsWithTimestamps const&, aod::Tracks const& tracks)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (bc.runNumber() != mRunNumber) {
      LOGF(debug, "timestamp=%llu", bc.timestamp());
      TList* callst = ccdb->getForTimeStamp<TList>("Centrality/Estimators", bc.timestamp());

      if (callst != nullptr) {
        auto getccdb = [callst](const char* ccdbhname) {
          TH1* h = (TH1*)callst->FindObject(ccdbhname);
          return h;
        };
        if (estV0M == 1) {
          LOGF(debug, "Getting new histograms with %d run number for %d run number", mRunNumber, bc.runNumber());
          mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
          mhVtxAmpCorrV0C = getccdb("hVtx_fAmplitude_V0C_Normalized");
          mhMultSelCalibV0M = getccdb("hMultSelCalib_V0M");
          if ((mhVtxAmpCorrV0A != nullptr) and (mhVtxAmpCorrV0C != nullptr) and (mhMultSelCalibV0M != nullptr)) {
            mV0MCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from V0M for run %d corrupted");
          }
        }
        if (estRun2SPD == 1) {
          LOGF(fatal, "Run2 calibration information estimated from SPD tracklets still not available");
        }
        if (estRun2CL0 == 1) {
          LOGF(fatal, "Run2 calibration information estimated from CL0 still not available");
        }
        if (estRun2CL1 == 1) {
          LOGF(fatal, "Run2 calibration information estimated from CL1 still not available");
        }
        if (mV0MCalibrationStored) {
          mRunNumber = bc.runNumber();
        }
      } else {
        /* we dont change the run number to keep trying */
        mV0MCalibrationStored = false;
        LOGF(fatal, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
    }
    if (estV0M == 1) {
      float centV0M = 105.0f;
      if (mV0MCalibrationStored) {
        float v0m = collision.multV0A() * mhVtxAmpCorrV0A->GetBinContent(mhVtxAmpCorrV0A->FindFixBin(collision.posZ())) +
                    collision.multV0C() * mhVtxAmpCorrV0C->GetBinContent(mhVtxAmpCorrV0C->FindFixBin(collision.posZ()));
        centV0M = mhMultSelCalibV0M->GetBinContent(mhMultSelCalibV0M->FindFixBin(v0m));
      }
      LOGF(debug, "centV0M=%.0f", centV0M);
      // fill centrality columns
      centVOM(centV0M);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityTable>(cfgc)};
}
