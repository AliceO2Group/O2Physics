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
  OutputObj<TH1F> hCentRun2V0M{TH1F("hCentRun2V0M", "V0M", 21, 0, 105.)};
  OutputObj<TH1F> hCentRun2SPDTks{TH1F("hCentRun2SPDTks", "SPD Tracklets", 21, 0, 105.)};
  OutputObj<TH1F> hCentRun2SPDCls{TH1F("hCentRun2SPDCls", "SPD Clusters", 21, 0, 105.)};
  OutputObj<TH1F> hCentRun2CL0{TH1F("hCentRun2CL0", "CL0", 21, 0, 105.)};
  OutputObj<TH1F> hCentRun2CL1{TH1F("hCentRun2CL1", "CL1", 21, 0, 105.)};
  OutputObj<TH1F> hCentRun3V0M{TH1F("hCentRun3V0M", "V0M", 21, 0, 105.)};
  OutputObj<TH1F> hCentRun3FT0M{TH1F("hCentRun3FT0M", "FT0M", 21, 0, 105.)};
  OutputObj<TH1F> hCentRun3FDDM{TH1F("hCentRun3FDDM", "FDDM", 21, 0, 105.)};
  OutputObj<TH1F> hCentRun3NTPV{TH1F("hCentRun3NTPV", "NTPV", 21, 0, 105.)};
  void processRun2PP(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2SPDClss>::iterator const& col)
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
    hCentRun2V0M->Fill(col.centRun2V0M());
    hCentRun2SPDTks->Fill(col.centRun2SPDTracklets());
    hCentRun2SPDCls->Fill(col.centRun2SPDClusters());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PP, "Process with Run2 SPD clusters centrality/multiplicity estimation", false);

  void processRun2PbPb(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2CL0s, aod::CentRun2CL1s>::iterator const& col)
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
    hCentRun2V0M->Fill(col.centRun2V0M());
    hCentRun2SPDTks->Fill(col.centRun2SPDTracklets());
    hCentRun2CL0->Fill(col.centRun2CL0());
    hCentRun2CL1->Fill(col.centRun2CL1());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PbPb, "Process with Run2 CL0 and CL1 multiplicities centrality/multiplicity  estimation", false);

  void processRun3(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun3V0Ms, aod::CentRun3FT0Ms, aod::CentRun3FDDMs, aod::CentRun3NTPVs>::iterator const& col)
  {
    if (!col.sel8()) {
      return;
    }

    LOGF(debug, "centV0M=%.0f", col.centRun3V0M());
    LOGF(debug, "centFT0M=%.0f", col.centRun3FT0M());
    LOGF(debug, "centFDDM=%.0f", col.centRun3FDDM());
    LOGF(debug, "centNTPV=%.0f", col.centRun3NTPV());
    // fill centrality histos
    hCentRun3V0M->Fill(col.centRun3V0M());
    hCentRun3FT0M->Fill(col.centRun3FT0M());
    hCentRun3FDDM->Fill(col.centRun3FDDM());
    hCentRun3NTPV->Fill(col.centRun3NTPV());
  }
  PROCESS_SWITCH(CentralityQa, processRun3, "Process with Run3 centrality/multiplicity estimators", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityQa>(cfgc)};
}
