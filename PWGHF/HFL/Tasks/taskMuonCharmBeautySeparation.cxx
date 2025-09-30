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
/// \file taskMuonCharmBeautySeparation.cxx
/// \note This workflow requires o2-analysis-fwdtrackextension as a dependency.
/// \brief Task to estimate HF->mu in the forward direction and use DCA observable to separate b-> mu, c-> mu.
/// \author Shreyasi Acharya <shreyasi.acharya@cern.ch>, LPC, France

#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskMuonCharmBeautySeparation {
  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    AxisSpec const trackTypeAxis = {6, -0.5, 5.5, "Track Type"};
    AxisSpec const ptRecoAxis = {1500, 0, 15, "#it{p}_{T}_{Reco}"};
    AxisSpec const dcaxAxis = {1000, -5.0, 5.0, "DCA {x or y} (cm)"};
    AxisSpec const dcaAxis = {1000, 0.0, 100.0, "DCA {xy} (cm)"};
    AxisSpec const zvtxAxis = {400, -20.0, 20.0, "zvtx (cm)"};
    AxisSpec const etaRecoAxis = {150, -5, -2, "#eta_{Reco}"};
    AxisSpec const rAbsAxis = {100, 0, 100, "R_{abs}"};
    AxisSpec const pdcaAxis = {450, 0, 450, "p_{DCA}"};
    AxisSpec const chi2GlobalAxis = {170, -1.5, 150.5, "#chi^{2} global"};
    AxisSpec const chi2MCHMFTAxis = {170, -1.5, 150.5, "#chi^{2} MCH-MFT"};
    AxisSpec const chi2MCHMIDAxis = {170, -1.5, 150.5, "#chi^{2} MCH-MID"};

    HistogramConfigSpec const histVariable({HistType::kTHnSparseF, {ptRecoAxis, dcaxAxis, dcaxAxis, dcaAxis, zvtxAxis}});
    registry.add("hBasicDist", "", histVariable);
    registry.add("hTrackType", "hTrackType", {HistType::kTH1F, {trackTypeAxis}});
    registry.add("hZvtx", "Zvtx in cm", {HistType::kTH1F, {zvtxAxis}});
    registry.add("hZvtx_WithMuons", "Zvtx with muons", {HistType::kTH1F, {zvtxAxis}});
    registry.add("hSign", "Sign of the track eletric charge", {HistType::kTH1F, {{5, -2.5, 2.5}}});
    registry.add("hForwardMultiplicity", "Multiplicity in forward direction", {HistType::kTH1F, {{20, 0, 20}}});
  }

  void process(aod::Collisions::iterator const& collision,
               soa::Join<aod::FwdTracks, aod::FwdTracksDCA> const& tracks)
  {
    auto pt = 0.;
    auto dcax = 0.;
    auto dcay = 0.;
    auto dca = 0.;
    auto chargeSign = 0.;
    auto zvtx = 0.;
    auto nFwdTracks = 0.;
    zvtx = collision.posZ();

    registry.fill(HIST("hZvtx"), zvtx);

    for (const auto& muon : tracks) {
      registry.fill(HIST("hTrackType"), muon.trackType());
      if (muon.has_collision()) {
        if (muon.trackType() == 0) {
          nFwdTracks++;
          pt = muon.pt();
          dcax = muon.fwdDcaX();
          dcay = muon.fwdDcaY();
          dca = std::sqrt(dcax * dcax + dcay * dcay);
          chargeSign = muon.sign();
          registry.fill(HIST("hBasicDist"), pt, dcax, dcay, dca, zvtx);
          registry.fill(HIST("hSign"), chargeSign);
        }
      }
    }
    registry.fill(HIST("hForwardMultiplicity"), nFwdTracks);
    if (nFwdTracks > 0) {
      registry.fill(HIST("hZvtx_WithMuons"), zvtx);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskMuonCharmBeautySeparation>(cfgc),
  };
}
