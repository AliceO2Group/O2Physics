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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "TH1F.h"

using namespace o2;
using namespace o2::framework;

struct CentralityQa {
  OutputObj<TH1F> hCentV0M{TH1F("hCentV0M", "V0M", 21, 0, 105.)};
  OutputObj<TH1F> hCentSPDTks{TH1F("hCentSPDTks", "SPD Tracklets", 21, 0, 105.)};
  OutputObj<TH1F> hCentSPDCls{TH1F("hCentSPDCls", "SPD Clusters", 21, 0, 105.)};
  OutputObj<TH1F> hCentCL0{TH1F("hCentCL0", "CL0", 21, 0, 105.)};
  OutputObj<TH1F> hCentCL1{TH1F("hCentCL1", "CL1", 21, 0, 105.)};
  void processPP(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2SPDClss>::iterator const& col)
  {
    if (!col.alias()[kINT7]) {
      return;
    }
    if (!col.sel7()) {
      return;
    }

    LOGF(debug, "centV0M=%.0f", col.centRun2V0M());
    LOGF(debug, "centSPDTracklets=%.0f", col.centRun2SPDTracklets());
    LOGF(debug, "centSPDClusters=%.0f", col.centRun2SPDClusters());
    // fill centrality histos
    hCentV0M->Fill(col.centRun2V0M());
    hCentSPDTks->Fill(col.centRun2SPDTracklets());
    hCentSPDCls->Fill(col.centRun2SPDClusters());
  }
  PROCESS_SWITCH(CentralityQa, processPP, "Process with SPD clusters centrality/multiplicity estimation", false);

  void processPbPb(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2CL0s, aod::CentRun2CL1s>::iterator const& col)
  {
    if (!col.alias()[kINT7]) {
      return;
    }
    if (!col.sel7()) {
      return;
    }

    LOGF(debug, "centV0M=%.0f", col.centRun2V0M());
    LOGF(debug, "centSPDTracklets=%.0f", col.centRun2SPDTracklets());
    LOGF(debug, "centCL0=%.0f", col.centRun2CL0());
    LOGF(debug, "centCL1=%.0f", col.centRun2CL1());
    // fill centrality histos
    hCentV0M->Fill(col.centRun2V0M());
    hCentSPDTks->Fill(col.centRun2SPDTracklets());
    hCentCL0->Fill(col.centRun2CL0());
    hCentCL1->Fill(col.centRun2CL1());
  }
  PROCESS_SWITCH(CentralityQa, processPbPb, "Process with CL0 and CL1 multiplicities centrality/multiplicity  estimation", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityQa>(cfgc)};
}
