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
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TProfile.h>

using namespace o2;
using namespace o2::framework;

struct CentralityQa {
  Configurable<int> nBins{"nBins", 1050, "number of bins"};
  Configurable<bool> INELgtZERO{"INELgtZERO", 1, "0 - no, 1 - yes"};
  OutputObj<TH1F> hCentRun2V0M{TH1F("hCentRun2V0M", "V0M", nBins, 0, 105.)};
  OutputObj<TH1F> hCentRun2V0A{TH1F("hCentRun2V0A", "V0A", nBins, 0, 105.)};
  OutputObj<TH1F> hCentRun2SPDTks{TH1F("hCentRun2SPDTks", "SPD Tracklets", nBins, 0, 105.)};
  OutputObj<TH1F> hCentRun2SPDCls{TH1F("hCentRun2SPDCls", "SPD Clusters", nBins, 0, 105.)};
  OutputObj<TH1F> hCentRun2CL0{TH1F("hCentRun2CL0", "CL0", nBins, 0, 105.)};
  OutputObj<TH1F> hCentRun2CL1{TH1F("hCentRun2CL1", "CL1", nBins, 0, 105.)};
  OutputObj<TH1F> hCentFV0A{TH1F("hCentFV0A", "FV0A", nBins, 0, 105.)};
  OutputObj<TH1F> hCentFT0M{TH1F("hCentFT0M", "FT0M", nBins, 0, 105.)};
  OutputObj<TH1F> hCentFT0A{TH1F("hCentFT0A", "FT0A", nBins, 0, 105.)};
  OutputObj<TH1F> hCentFT0C{TH1F("hCentFT0C", "FT0C", nBins, 0, 105.)};
  OutputObj<TH1F> hCentFDDM{TH1F("hCentFDDM", "FDDM", nBins, 0, 105.)};
  OutputObj<TH1F> hCentNTPV{TH1F("hCentNTPV", "NTPV", nBins, 0, 105.)};

  // profiles of midrapidity multiplicity density
  OutputObj<TProfile> hCentProfileFV0A{TProfile("hCentProfileFV0A", "FV0A", nBins, 0, 105.)};
  OutputObj<TProfile> hCentProfileFT0M{TProfile("hCentProfileFT0M", "FT0M", nBins, 0, 105.)};
  OutputObj<TProfile> hCentProfileFT0A{TProfile("hCentProfileFT0A", "FT0A", nBins, 0, 105.)};
  OutputObj<TProfile> hCentProfileFT0C{TProfile("hCentProfileFT0C", "FT0C", nBins, 0, 105.)};
  OutputObj<TProfile> hCentProfileFDDM{TProfile("hCentProfileFDDM", "FDDM", nBins, 0, 105.)};
  OutputObj<TProfile> hCentProfileNTPV{TProfile("hCentProfileNTPV", "NTPV", nBins, 0, 105.)};

  void processRun2PP(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2SPDClss, aod::Mults>::iterator const& col)
  {
    if (!col.alias_bit(kINT7)) {
      return;
    }
    if (!col.sel7()) {
      return;
    }
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
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

  void processRun2PbPb(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2CL0s, aod::CentRun2CL1s, aod::Mults>::iterator const& col)
  {
    if (!col.alias_bit(kINT7)) {
      return;
    }
    if (!col.sel7()) {
      return;
    }
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
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

  void processRun2PPb(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0As, aod::Mults>::iterator const& col)
  {
    if (!col.alias_bit(kINT7)) {
      return;
    }
    if (!col.sel7()) {
      return;
    }
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    LOGF(debug, "centV0A=%.0f", col.centRun2V0A());
    // fill centrality histos
    hCentRun2V0A->Fill(col.centRun2V0A());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PPb, "Process with Run2 V0A multiplicitY centrality/multiplicity  estimation", false);

  void processRun3_FV0A(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As>::iterator const& col)
  {
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    if (!col.sel8())
      return;
    LOGF(debug, "centFV0A=%.0f", col.centFV0A());
    hCentFV0A->Fill(col.centFV0A());
    hCentProfileFV0A->Fill(col.centFV0A(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FV0A, "Process with Run 3 FV0A estimator", false);

  void processRun3_FT0M(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator const& col)
  {
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    if (!col.sel8())
      return;
    LOGF(debug, "centFT0M=%.0f", col.centFT0M());
    hCentFT0M->Fill(col.centFT0M());
    hCentProfileFT0M->Fill(col.centFT0M(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0M, "Process with Run 3 FT0M estimator", false);

  void processRun3_FT0A(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0As>::iterator const& col)
  {
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    if (!col.sel8())
      return;
    hCentFT0A->Fill(col.centFT0A());
    hCentProfileFT0A->Fill(col.centFT0A(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0A, "Process with Run 3 FT0A estimator", false);

  void processRun3_FT0C(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>::iterator const& col)
  {
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    if (!col.sel8())
      return;
    hCentFT0C->Fill(col.centFT0C());
    hCentProfileFT0C->Fill(col.centFT0C(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0C, "Process with Run 3 FT0A estimator", false);

  void processRun3_FDDM(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFDDMs>::iterator const& col)
  {
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    if (!col.sel8())
      return;
    hCentFDDM->Fill(col.centFDDM());
    hCentProfileFDDM->Fill(col.centFDDM(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FDDM, "Process with Run 3 FDDM estimator", false);

  void processRun3_NTPV(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentNTPVs>::iterator const& col)
  {
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    if (!col.sel8())
      return;
    hCentNTPV->Fill(col.centNTPV());
    hCentProfileNTPV->Fill(col.centNTPV(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_NTPV, "Process with Run 3 NTPV estimator", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityQa>(cfgc)};
}
