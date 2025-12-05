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
// This code produces muon table 001 from 000.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct muonConverter1 {
  Produces<aod::EMPrimaryMuons_001> muon_001;

  void process(aod::EMPrimaryMuons_000 const& muons)
  {
    for (const auto& muon : muons) {
      muon_001(
        muon.collisionId(),
        muon.fwdtrackId(), muon.mfttrackId(), muon.mchtrackId(), muon.trackType(),
        muon.pt(), muon.eta(), muon.phi(), muon.sign(),
        muon.fwdDcaX(), muon.fwdDcaY(), muon.cXXatDCA(), muon.cYYatDCA(), muon.cXYatDCA(),
        muon.ptMatchedMCHMID(), muon.etaMatchedMCHMID(), muon.phiMatchedMCHMID(),
        0, 0, 0, 0,
        muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(),
        muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
        muon.mchBitMap(), muon.midBitMap(), muon.midBoards(),
        muon.mftClusterSizesAndTrackFlags(), muon.chi2MFT(), muon.isAssociatedToMPC(), muon.isAmbiguous());
    } // end of muon loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<muonConverter1>(cfgc, TaskName{"muon-converter1"})};
}
