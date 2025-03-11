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

/// \file taskSingleMuonMult.cxx
/// \brief Task used to study the Open heavy flavour decay muon production as a function of multiplicity.
/// \author Md Samsul Islam <md.samsul.islam@cern.ch>, IITB

#include <cmath>
#include <vector>
#include <TPDGCode.h>
#include <TString.h>
#include "Framework/StaticFor.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/FT0Corrected.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"

#include "Common/CCDB/EventSelectionParams.h"

#include "ReconstructionDataFormats/TrackFwd.h"

#include "TableHelper.h"

#include "PWGLF/DataModel/LFResonanceTables.h"
#include "CommonConstants/MathConstants.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;
using MyTracks = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksIU, aod::TracksDCA, aod::TrackSelection>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
using MyMcMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MFTTracksExtra = soa::Join<aod::MFTTracks>;

struct HfTaskSingleMuonMult {
  Configurable<float> etaMin{"etaMin", -3.6, "eta minimum value"};
  Configurable<float> etaMax{"etaMax", -2.5, "eta maximum value"};
  Configurable<float> pDcaMin{"pDcaMin", 324., "p*DCA maximum value for small Rabs"};
  Configurable<float> pDcaMax{"pDcaMax", 594., "p*DCA maximum value for large Rabs"};
  Configurable<float> rAbsMin{"rAbsMin", 17.6, "R at absorber end minimum value"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5, "R at absorber end maximum value"};
  Configurable<float> rAbsMid{"rAbsMid", 26.5, "R at absorber end split point for different p*DCA selections"};
  Configurable<float> zVtx{"zVtx", 10., "Z edge of primary vertex [cm]"};
  Configurable<bool> reduceOrphMft{"reduceOrphMft", true, "reduce orphan MFT tracks"};

  o2::framework::HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  static constexpr std::string_view kTrackType[] = {"TrackType0", "TrackType1", "TrackType2", "TrackType3", "TrackType4"};
  uint8_t globalMuonTrack = o2::aod::fwdtrack::GlobalMuonTrack;

  void init(InitContext&)
  {
    AxisSpec axisCent = {101, -0.5, 100.5, "centrality"};
    AxisSpec axisEvent{3, 0.5, 3.5, "Event Selection"};
    AxisSpec axisEventSize{500, 0.5, 500.5, "Event Size"};
    AxisSpec axisVtxZ{80, -20., 20., "#it{z}_{vtx} (cm)"};
    AxisSpec axisMuTrk{5, 0.5, 5.5, "Muon Selection"};
    AxisSpec axisNch{500, 0.5, 500.5, "N_{ch}"};
    AxisSpec axisNmu{20, -0.5, 19.5, "N_{#mu}"};
    AxisSpec axisPt{1000, 0., 500., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{250, -5., 5., "#it{#eta}"};
    AxisSpec axisTheta{500, 170., 180., "#it{#theta}"};
    AxisSpec axisRabs{1000, 0., 100., "R_{abs} (cm)"};
    AxisSpec axisDCA{500, 0., 5., "#it{DCA}_{xy} (cm)"};
    AxisSpec axisChi2MatchMCHMFT{1000, 0., 1000., "MCH-MFT matching #chi^{2}"};
    AxisSpec axisSign{5, -2.5, 2.5, "Charge"};
    AxisSpec axisPDca{100000, 0, 100000, "#it{p} #times DCA (GeV/#it{c} * cm)"};
    AxisSpec axisDCAx{1000, -5., 5., "#it{DCA}_{x or y} (cm)"};
    AxisSpec axisEtaDif{200, -2., 2., "#it{#eta} diff"};
    AxisSpec axisDeltaPt{10000, -50, 50, "#Delta #it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisTrackType{8, -1.5, 6.5, "TrackType"};
    AxisSpec axisPtDif{200, -2., 2., "#it{p}_{T} diff (GeV/#it{c})"};

    HistogramConfigSpec hCentrality{HistType::kTH1F, {axisCent}};
    HistogramConfigSpec hEvent{HistType::kTH1F, {axisEvent}};
    HistogramConfigSpec hEventSize{HistType::kTH1F, {axisEventSize}};
    HistogramConfigSpec hVtxZ{HistType::kTH1F, {axisVtxZ}};

    HistogramConfigSpec hMuTrkSel{HistType::kTH1F, {axisMuTrk}};
    HistogramConfigSpec hTHnMu{HistType::kTHnSparseF, {axisCent, axisNch, axisPt, axisEta, axisTheta, axisRabs, axisDCA, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10};
    HistogramConfigSpec hTHnMuDeltaPt{HistType::kTHnSparseF, {axisCent, axisNch, axisPt, axisEta, axisTheta, axisRabs, axisDCA, axisPDca, axisChi2MatchMCHMFT, axisDeltaPt}, 10};
    HistogramConfigSpec h3DCA{HistType::kTH3F, {axisDCAx, axisDCAx, axisTrackType}};
    HistogramConfigSpec hTHnCh{HistType::kTHnSparseF, {axisCent, axisNch, axisPt, axisEta, axisSign}, 5};
    HistogramConfigSpec h3MultNchNmu{HistType::kTH3F, {axisCent, axisNch, axisNmu}};

    HistogramConfigSpec h2PtMc{HistType::kTH2F, {axisPt, axisPtDif}};
    HistogramConfigSpec h2EtaMc{HistType::kTH2F, {axisEta, axisEtaDif}};

    registry.add("hCentrality", "", hCentrality);
    registry.add("hEvent", "", hEvent);
    registry.add("hEventSize", "", hEventSize);
    registry.add("hVtxZBeforeSel", "", hVtxZ);
    registry.add("hVtxZAfterSel", "", hVtxZ);

    registry.add("hMuTrkSel", "", hMuTrkSel);
    registry.add("hMuBeforeMatchMFT", "", hTHnMu);
    registry.add("hMuBeforeAccCuts", "", hTHnMu);
    registry.add("h3DCABeforeAccCuts", "", h3DCA);
    registry.add("hMuDeltaPtBeforeAccCuts", "", hTHnMuDeltaPt);
    registry.add("hMuAfterEtaCuts", "", hTHnMu);
    registry.add("hMuAfterRabsCuts", "", hTHnMu);
    registry.add("hMuAfterPdcaCuts", "", hTHnMu);
    registry.add("hMuAfterAccCuts", "", hTHnMu);
    registry.add("h3DCAAfterAccCuts", "", h3DCA);
    registry.add("hMuDeltaPtAfterAccCuts", "", hTHnMuDeltaPt);

    registry.add("hTHnTrk", "", hTHnCh);
    registry.add("h3MultNchNmu", "", h3MultNchNmu);

    auto hEvstat = registry.get<TH1>(HIST("hEvent"));
    auto* xEv = hEvstat->GetXaxis();
    xEv->SetBinLabel(1, "All events");
    xEv->SetBinLabel(2, "sel8");
    xEv->SetBinLabel(3, "VtxZAfterSel");

    auto hMustat = registry.get<TH1>(HIST("hMuTrkSel"));
    auto* xMu = hMustat->GetXaxis();
    xMu->SetBinLabel(1, "noCut");
    xMu->SetBinLabel(2, "etaCut");
    xMu->SetBinLabel(3, "RabsCut");
    xMu->SetBinLabel(4, "pDcaCut");
    xMu->SetBinLabel(5, "chi2Cut");

    const uint8_t muonTrackType[]{0, 1, 2, 3, 4};
    for (const auto& trktype : muonTrackType) {
      registry.add(Form("h3MultNchNmu_TrackType%d", trktype), "", h3MultNchNmu);
    }
  }

  template <typename TCollision, typename TTracks, typename TMuons>
  void runMuonSel(TCollision const& collision, TTracks const& tracks, TMuons const& muons)
  {
    registry.fill(HIST("hEvent"), 1);

    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEvent"), 2);
    registry.fill(HIST("hVtxZBeforeSel"), collision.posZ());
    LOGP(debug, "SAMSUL_{} collisions", collision.globalIndex());

    if (std::abs(collision.posZ()) > zVtx) {
      return;
    }
    registry.fill(HIST("hEvent"), 3);
    registry.fill(HIST("hVtxZAfterSel"), collision.posZ());

    const auto cent = collision.centFT0M();
    registry.fill(HIST("hCentrality"), cent);

    int nCh = 0.;
    int nMu = 0.;
    const int nTypes = 5;
    int nMuTrackType[nTypes] = {0};

    std::vector<typename std::decay_t<decltype(tracks)>::iterator> chTracks;
    for (const auto& track : tracks) {
      if (track.isGlobalTrack() && std::abs(track.eta()) < 0.8 && track.pt() > 0.15) {
        chTracks.push_back(track);
      }
    }
    nCh = chTracks.size();
    if (nCh < 1) {
      return;
    }
    registry.fill(HIST("hEventSize"), nCh);

    for (int isize = 0; isize < nCh; isize++) {
      auto chTrack = chTracks[isize];
      registry.fill(HIST("hTHnTrk"), cent, nCh, chTrack.pt(), chTrack.eta(), chTrack.sign());
    }

    // muons
    for (const auto& muon : muons) {
      const auto pt(muon.pt()), eta(muon.eta()), theta(90 - ((std::atan(muon.tgl())) * (180. / constants::math::PI))), pDca(muon.pDca()), rAbs(muon.rAtAbsorberEnd()), chi2(muon.chi2MatchMCHMFT());
      const auto dcaXY(RecoDecay::sqrtSumOfSquares(muon.fwdDcaX(), muon.fwdDcaY()));
      const auto muTrackType(muon.trackType());

      registry.fill(HIST("hMuBeforeMatchMFT"), cent, nCh, pt, eta, theta, rAbs, dcaXY, pDca, chi2, muTrackType);

      // histograms before the acceptance cuts
      registry.fill(HIST("hMuTrkSel"), 1);
      registry.fill(HIST("hMuBeforeAccCuts"), cent, nCh, pt, eta, theta, rAbs, dcaXY, pDca, chi2, muTrackType);
      registry.fill(HIST("h3DCABeforeAccCuts"), muon.fwdDcaX(), muon.fwdDcaY(), muTrackType);

      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<MyMuons>();
        auto dpt(muonType3.pt() - pt);
        if (muTrackType == globalMuonTrack) {
          registry.fill(HIST("hMuDeltaPtBeforeAccCuts"), cent, nCh, pt, eta, theta, rAbs, dcaXY, pDca, chi2, dpt);
        }
      }

      // Apply various standard muon acceptance cuts
      // eta cuts
      if ((eta >= etaMax) || (eta < etaMin)) {
        continue;
      }
      registry.fill(HIST("hMuTrkSel"), 2);
      registry.fill(HIST("hMuAfterEtaCuts"), cent, nCh, pt, eta, theta, rAbs, dcaXY, pDca, chi2, muTrackType);

      // Rabs cuts
      if ((rAbs < rAbsMin) || (rAbs >= rAbsMax)) {
        continue;
      }
      registry.fill(HIST("hMuTrkSel"), 3);
      registry.fill(HIST("hMuAfterRabsCuts"), cent, nCh, pt, eta, theta, rAbs, dcaXY, pDca, chi2, muTrackType);

      if ((rAbs < rAbsMid) && (pDca >= pDcaMin)) {
        continue;
      }
      if ((rAbs >= rAbsMid) && (pDca >= pDcaMax)) {
        continue;
      }
      registry.fill(HIST("hMuTrkSel"), 4);
      registry.fill(HIST("hMuAfterPdcaCuts"), cent, nCh, pt, eta, theta, rAbs, dcaXY, pDca, chi2, muTrackType);

      //  MCH-MFT matching chi2
      if ((muon.chi2() >= 1e6) || (muon.chi2() < 0)) {
        continue;
      }
      registry.fill(HIST("hMuTrkSel"), 5);

      // histograms after acceptance cuts
      registry.fill(HIST("hMuAfterAccCuts"), cent, nCh, pt, eta, theta, rAbs, dcaXY, pDca, chi2, muTrackType);
      registry.fill(HIST("h3DCAAfterAccCuts"), muon.fwdDcaX(), muon.fwdDcaY(), muTrackType);
      nMu++;
      nMuTrackType[muTrackType]++;

      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<MyMuons>();
        auto dpt(muonType3.pt() - pt);

        if (muTrackType == globalMuonTrack) {
          registry.fill(HIST("hMuDeltaPtAfterAccCuts"), cent, nCh, pt, eta, theta, rAbs, dcaXY, pDca, chi2, dpt);
        }
      }
    }

    registry.fill(HIST("h3MultNchNmu"), cent, nCh, nMu);

    static_for<0, 4>([&](auto i) {
      constexpr int kIndex = i.value;
      if (nMuTrackType[kIndex] > 0)
        registry.fill(HIST("h3MultNchNmu_") + HIST(kTrackType[kIndex]), cent, nCh, nMuTrackType[kIndex]);
    });
    chTracks.clear();
  }
  void process(MyCollisions::iterator const& collision,
               MyTracks const& tracks,
               MyMuons const& muons)
  {
    runMuonSel(collision, tracks, muons);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskSingleMuonMult>(cfgc)};
}
