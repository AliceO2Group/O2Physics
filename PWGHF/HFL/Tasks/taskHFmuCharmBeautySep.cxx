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
///
/// \Task to estimate HF->mu in the forward direction and use DCA observable to separate b-> mu, c-> mu.
/// \Command to run on terminal : o2-analysis-hfmu-charm-beauty-separation | o2-analysis-fwdtrackextension --aod-file AO2D.root
/// \Author : Shreyasi Acharya <shreyasi.acharya@cern.ch>, LPC, France

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "TDatabasePDG.h"
#include "PWGDQ/Core/VarManager.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
// Iterate on muon using the collision iterator in the dq-analysis style
struct IterateFwdTracksMuons {

  HistogramRegistry spectra{"ForwardTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  void init(o2::framework::InitContext&)
  {
    AxisSpec trackTypeAxis = {6, -0.5, 5.5, "Track Type"};
    AxisSpec ptRecoAxis = {1500, 0, 15, "#it{p}_{T}_{Reco}"};
    AxisSpec dcaxAxis = {1000, -5.0, 5.0, "DCA{x or y} (cm)"};
    AxisSpec etaRecoAxis = {150, -5, -2, "#eta_{Reco}"};
    AxisSpec rAbsAxis = {100, 0, 100, "R_{abs}"};
    AxisSpec pdcaAxis = {450, 0, 450, "p_{DCA}"};
    AxisSpec chi2Axis = {170, -1.5, 150.5, "#chi^{2} global"};
    AxisSpec chi2MCHMFTAxis = {170, -1.5, 150.5, "#chi^{2} MCH-MFT"};
    AxisSpec chi2MCHMIDAxis = {170, -1.5, 150.5, "#chi^{2} MCH-MID"};

    HistogramConfigSpec HistVariable({HistType::kTHnSparseF, {ptRecoAxis, dcaxAxis, etaRecoAxis, chi2MCHMFTAxis, chi2Axis, chi2MCHMIDAxis, rAbsAxis, pdcaAxis}});
    spectra.add("BasicDist", "", HistVariable);
    HistogramConfigSpec HistTrackType({HistType::kTH1F, {trackTypeAxis}});
    spectra.add("TrackType", "", HistTrackType);

    spectra.add("hDCAxMuons_Type0", " dca", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
    spectra.add("hDCAyMuons_Type0", " dca", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
    spectra.add("hDCAxyMuons_Type0", " dca", {HistType::kTH1F, {{1000, 0.0, 10.0}}});
  }
  void process(aod::Collisions::iterator const& collision, soa::Join<aod::FwdTracks, aod::FwdTracksDCA> const& muons)
  {
    auto pt = 0.;
    auto dcax = 0.;
    auto dcay = 0.;
    auto eta = 0.;
    auto chi2MatchMCHMFT = 0.;
    auto chi2MatchMCHMID = 0.;
    auto chi2 = 0.;
    auto rAbs = 0.;
    auto pDca = 0.;

    for (auto& muon : muons) {
      spectra.fill(HIST("TrackType"), muon.trackType());
      if (muon.has_collision()) {
        if (muon.trackType() == 0) {

          pt = muon.pt();
          dcax = muon.fwdDcaX();
          dcay = muon.fwdDcaY();
          eta = muon.eta();
          chi2MatchMCHMFT = muon.chi2MatchMCHMFT();
          chi2MatchMCHMID = muon.chi2MatchMCHMID();
          chi2 = muon.chi2();
          rAbs = muon.rAtAbsorberEnd();
          pDca = muon.pDca();

          spectra.fill(HIST("BasicDist"), pt, dcax, eta, chi2MatchMCHMFT, chi2, chi2MatchMCHMID, rAbs, pDca);
          spectra.fill(HIST("hDCAxMuons_Type0"), dcax);
          spectra.fill(HIST("hDCAyMuons_Type0"), dcay);
          spectra.fill(HIST("hDCAxyMuons_Type0"), std::sqrt(dcax * dcax + dcay * dcay));
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<IterateFwdTracksMuons>(cfgc),
  };
}
