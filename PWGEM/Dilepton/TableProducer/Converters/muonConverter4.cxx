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
// ========================
//
// This code produces muon table 004 from 003.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct muonConverter4 {
  Produces<o2::aod::EMPrimaryMuons_004> muon_004;
  Produces<o2::aod::EMPrimaryMuonsCov_003> muoncov_003;

  void processMuon(o2::aod::EMPrimaryMuons_003 const& muons)
  {
    for (const auto& muon : muons) {
      muon_004(
        muon.collisionId(),
        muon.fwdtrackId(), muon.mfttrackId(), muon.mchtrackId(), muon.trackType(),
        muon.pt(), muon.eta(), muon.phi(), muon.sign(),
        muon.fwdDcaX(), muon.fwdDcaY(), muon.cXX(), muon.cYY(), muon.cXY(),
        muon.ptMatchedMCHMID(), muon.etaMatchedMCHMID(), muon.phiMatchedMCHMID(),
        muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(),
        muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(), muon.diffChi2MatchingMCHMFT(),
        muon.mchBitMap(), muon.midBitMap(), muon.midBoards(),
        muon.mftClusterSizesAndTrackFlags(), muon.chi2MFT(), muon.isAssociatedToMPC(), muon.isAmbiguous());
    } // end of muon loop
  }
  PROCESS_SWITCH(muonConverter4, processMuon, "convert primary muon 003 into 004", true);

  void processMuonCov(soa::Join<o2::aod::EMPrimaryMuons_003, o2::aod::EMPrimaryMuonsCov_002> const& muons)
  {
    for (const auto& muon : muons) {
      muoncov_003(muon.x(), muon.y(), muon.y(),
                  muon.cXX(),
                  muon.cYY(), muon.cXY(),
                  muon.cPhiX(), muon.cPhiY(), muon.cPhiPhi(),
                  muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                  muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2());
    } // end of muon loop
  }
  PROCESS_SWITCH(muonConverter4, processMuonCov, "convert primary muon 002 into 003", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<muonConverter4>(cfgc, TaskName{"muon-converter4"})};
}
