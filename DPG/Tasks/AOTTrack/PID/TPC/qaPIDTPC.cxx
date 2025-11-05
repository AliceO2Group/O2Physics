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
/// \file   qaPIDTPC.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Implementation for QA tasks of the TPC PID quantities
///

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the TPC QA plots
struct tpcPidQa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};
  static constexpr std::string_view hdelta[Np] = {"delta/El", "delta/Mu", "delta/Pi",
                                                  "delta/Ka", "delta/Pr", "delta/De",
                                                  "delta/Tr", "delta/He", "delta/Al"};
  static constexpr std::string_view hdelta_pt[Np] = {"delta/pt/El", "delta/pt/Mu", "delta/pt/Pi",
                                                     "delta/pt/Ka", "delta/pt/Pr", "delta/pt/De",
                                                     "delta/pt/Tr", "delta/pt/He", "delta/pt/Al"};
  static constexpr std::string_view hexpsigma[Np] = {"expsigma/El", "expsigma/Mu", "expsigma/Pi",
                                                     "expsigma/Ka", "expsigma/Pr", "expsigma/De",
                                                     "expsigma/Tr", "expsigma/He", "expsigma/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hnsigma_pt[Np] = {"nsigma/pt/El", "nsigma/pt/Mu", "nsigma/pt/Pi",
                                                      "nsigma/pt/Ka", "nsigma/pt/Pr", "nsigma/pt/De",
                                                      "nsigma/pt/Tr", "nsigma/pt/He", "nsigma/pt/Al"};
  static constexpr std::string_view hnsigma_p_eta_Ncl[Np] = {"nsigma/sparsePinEtaNcl/El", "nsigma/sparsePinEtaNcl/Mu", "nsigma/sparsePinEtaNcl/Pi",
                                                             "nsigma/sparsePinEtaNcl/Ka", "nsigma/sparsePinEtaNcl/Pr", "nsigma/sparsePinEtaNcl/De",
                                                             "nsigma/sparsePinEtaNcl/Tr", "nsigma/sparsePinEtaNcl/He", "nsigma/sparsePinEtaNcl/Al"};

  // With TOF
  static constexpr std::string_view hexpected_wTOF[Np] = {"wTOF/expected/El", "wTOF/expected/Mu", "wTOF/expected/Pi",
                                                          "wTOF/expected/Ka", "wTOF/expected/Pr", "wTOF/expected/De",
                                                          "wTOF/expected/Tr", "wTOF/expected/He", "wTOF/expected/Al"};
  static constexpr std::string_view hdelta_wTOF[Np] = {"wTOF/delta/El", "wTOF/delta/Mu", "wTOF/delta/Pi",
                                                       "wTOF/delta/Ka", "wTOF/delta/Pr", "wTOF/delta/De",
                                                       "wTOF/delta/Tr", "wTOF/delta/He", "wTOF/delta/Al"};
  static constexpr std::string_view hdelta_pt_pos_wTOF[Np] = {"wTOF/delta/pt/pos/El", "wTOF/delta/pt/pos/Mu", "wTOF/delta/pt/pos/Pi",
                                                              "wTOF/delta/pt/pos/Ka", "wTOF/delta/pt/pos/Pr", "wTOF/delta/pt/pos/De",
                                                              "wTOF/delta/pt/pos/Tr", "wTOF/delta/pt/pos/He", "wTOF/delta/pt/pos/Al"};
  static constexpr std::string_view hdelta_pt_neg_wTOF[Np] = {"wTOF/delta/pt/neg/El", "wTOF/delta/pt/neg/Mu", "wTOF/delta/pt/neg/Pi",
                                                              "wTOF/delta/pt/neg/Ka", "wTOF/delta/pt/neg/Pr", "wTOF/delta/pt/neg/De",
                                                              "wTOF/delta/pt/neg/Tr", "wTOF/delta/pt/neg/He", "wTOF/delta/pt/neg/Al"};
  static constexpr std::string_view hexpsigma_wTOF[Np] = {"wTOF/expsigma/El", "wTOF/expsigma/Mu", "wTOF/expsigma/Pi",
                                                          "wTOF/expsigma/Ka", "wTOF/expsigma/Pr", "wTOF/expsigma/De",
                                                          "wTOF/expsigma/Tr", "wTOF/expsigma/He", "wTOF/expsigma/Al"};
  static constexpr std::string_view hnsigma_wTOF[Np] = {"wTOF/nsigma/El", "wTOF/nsigma/Mu", "wTOF/nsigma/Pi",
                                                        "wTOF/nsigma/Ka", "wTOF/nsigma/Pr", "wTOF/nsigma/De",
                                                        "wTOF/nsigma/Tr", "wTOF/nsigma/He", "wTOF/nsigma/Al"};
  static constexpr std::string_view hnsigma_pt_wTOF[Np] = {"wTOF/nsigma/pt/El", "wTOF/nsigma/pt/Mu", "wTOF/nsigma/pt/Pi",
                                                           "wTOF/nsigma/pt/Ka", "wTOF/nsigma/pt/Pr", "wTOF/nsigma/pt/De",
                                                           "wTOF/nsigma/pt/Tr", "wTOF/nsigma/pt/He", "wTOF/nsigma/pt/Al"};
  static constexpr std::string_view hsignal_wTOF[Np] = {"wTOF/signal/El", "wTOF/signal/Mu", "wTOF/signal/Pi",
                                                        "wTOF/signal/Ka", "wTOF/signal/Pr", "wTOF/signal/De",
                                                        "wTOF/signal/Tr", "wTOF/signal/He", "wTOF/signal/Al"};
  static constexpr std::string_view hnsigma_p_eta_Ncl_wTOF[Np] = {"wTOF/nsigma/sparsePinEtaNcl/El", "wTOF/nsigma/sparsePinEtaNcl/Mu", "wTOF/nsigma/sparsePinEtaNcl/Pi",
                                                                  "wTOF/nsigma/sparsePinEtaNcl/Ka", "wTOF/nsigma/sparsePinEtaNcl/Pr", "wTOF/nsigma/sparsePinEtaNcl/De",
                                                                  "wTOF/nsigma/sparsePinEtaNcl/Tr", "wTOF/nsigma/sparsePinEtaNcl/He", "wTOF/nsigma/sparsePinEtaNcl/Al"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 1, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  ConfigurableAxis etaBins{"etaBins", {100, -1.f, 1.f}, "Binning in eta"};
  ConfigurableAxis phiBins{"phiBins", {100, 0, TMath::TwoPi()}, "Binning in phi"};
  ConfigurableAxis trackLengthBins{"trackLengthBins", {100, 0, 1000.f}, "Binning in track length plot"};
  ConfigurableAxis deltaBins{"deltaBins", {200, -1000.f, 1000.f}, "Binning in Delta (dEdx - expected dEdx)"};
  ConfigurableAxis expSigmaBins{"expSigmaBins", {200, 0.f, 200.f}, "Binning in expected Sigma"};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {401, -10.025f, 10.025f}, "Binning in NSigma"};
  ConfigurableAxis dEdxBins{"dEdxBins", {5000, 0.f, 5000.f}, "Binning in dE/dx"};
  // Axes for optional THnSparse
  ConfigurableAxis binsPForSparse{"binsPForSparse", {200, 0.f, 20.f}, "Binning in momentum for optional THnSparse"};
  ConfigurableAxis binsEtaForSparse{"binsEtaForSparse", {50, -1.f, 1.f}, "Binning in eta for optional THnSparse"};
  ConfigurableAxis binsnSigmaForSparse{"binsnSigmaForSparse", {101, -7.575, 7.575}, "Binning in nsigma for optional THnSparse"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> applyRapidityCut{"applyRapidityCut", false, "Flag to apply rapidity cut"};
  Configurable<bool> splitSignalPerCharge{"splitSignalPerCharge", true, "Split the signal per charge (reduces memory footprint if off)"};
  Configurable<bool> enableDeDxPlot{"enableDeDxPlot", true, "Enables the dEdx plot (reduces memory footprint if off)"};
  Configurable<int16_t> minTPCNcls{"minTPCNcls", 0, "Minimum number or TPC Clusters for tracks"};
  ConfigurableAxis tpcNclsBins{"tpcNclsBins", {16, 0, 160}, "Binning in number of clusters in TPC"};
  Configurable<bool> fillTHnSparses{"fillTHnSparses", false, "Flag to fill multidimensional histograms for nsigma vs pt, eta, Ncls"};

  template <o2::track::PID::ID id>
  void initPerParticle(const AxisSpec& pAxis,
                       const AxisSpec& ptAxis,
                       const AxisSpec& dedxAxis,
                       const AxisSpec& chargeAxis)
  {
    static_assert(id >= 0 && id <= PID::Alpha && "Particle index outside limits");
    bool enableTOFHistos = false;
    bool enableFullHistos = false;
    int enabledProcesses = 0;
    switch (id) { // Skipping disabled particles
#define particleCase(particleId)                                                                     \
  case PID::particleId:                                                                              \
    if (!doprocess##particleId && !doprocessFull##particleId && !doprocessFullWithTOF##particleId) { \
      return;                                                                                        \
    }                                                                                                \
    if (doprocess##particleId) {                                                                     \
      enabledProcesses++;                                                                            \
    }                                                                                                \
    if (doprocessFull##particleId) {                                                                 \
      enableFullHistos = true;                                                                       \
      enabledProcesses++;                                                                            \
    }                                                                                                \
    if (doprocessFullWithTOF##particleId) {                                                          \
      enableFullHistos = true;                                                                       \
      enableTOFHistos = true;                                                                        \
      enabledProcesses++;                                                                            \
    }                                                                                                \
    LOGF(info, "Enabled TPC QA for %s %s", #particleId, pT[id]);                                     \
    break;

      particleCase(Electron);
      particleCase(Muon);
      particleCase(Pion);
      particleCase(Kaon);
      particleCase(Proton);
      particleCase(Deuteron);
      particleCase(Triton);
      particleCase(Helium3);
      particleCase(Alpha);
#undef particleCase
    }
    if (enabledProcesses != 1) {
      LOG(fatal) << "Cannot enable more than one process function per particle, check and retry!";
    }

    // NSigma
    const char* axisTitle = Form("N_{#sigma}^{TPC}(%s)", pT[id]);
    const AxisSpec nSigmaAxis{nSigmaBins, axisTitle};
    histos.add(hnsigma[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    if (splitSignalPerCharge) {
      histos.add(hnsigma_pt[id].data(), axisTitle, kTH3F, {ptAxis, nSigmaAxis, chargeAxis});
    } else {
      histos.add(hnsigma_pt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    }

    if (!enableFullHistos) { // Enabling only NSigma for tiny tables
      return;
    }

    // Exp signal
    const AxisSpec expAxis{1000, 0, 1000, Form("d#it{E}/d#it{x}(%s) A.U.", pT[id])};
    histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{deltaBins, Form("d#it{E}/d#it{x} - d#it{E}/d#it{x}(%s)", pT[id])};
    axisTitle = Form("#Delta^{TPC}(%s)", pT[id]);
    histos.add(hdelta[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
    if (splitSignalPerCharge) {
      histos.add(hdelta_pt[id].data(), axisTitle, kTH3F, {ptAxis, deltaAxis, chargeAxis});
    } else {
      histos.add(hdelta_pt[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
    }

    // Exp Sigma
    const AxisSpec expSigmaAxis{expSigmaBins, Form("Exp_{#sigma}^{TPC}(%s)", pT[id])};
    histos.add(hexpsigma[id].data(), "", kTH2F, {pAxis, expSigmaAxis});

    const AxisSpec etaAxis{etaBins, "#it{#eta}"};
    const AxisSpec tpcnclsAxis{tpcNclsBins, "TPC #cls"};
    const AxisSpec sparseMomentumAxis{binsPForSparse, "#it{p} (GeV/#it{c})"};
    const AxisSpec sparseEtaAxis{binsEtaForSparse, "#eta"};
    const AxisSpec sparseNSigmaAxis{binsnSigmaForSparse, "#n_{#sigma}^{TPC}"};
    HistogramConfigSpec particleSparseHists{HistType::kTHnSparseF, {sparseMomentumAxis, sparseEtaAxis, sparseNSigmaAxis, tpcnclsAxis}};
    if (fillTHnSparses) {
      histos.add(hnsigma_p_eta_Ncl[id].data(), axisTitle, particleSparseHists);
    }

    if (!enableTOFHistos) { // Returning if the plots with TOF are not requested
      return;
    }

    histos.add(hexpected_wTOF[id].data(), "With TOF", kTH2F, {pAxis, expAxis});
    histos.add(hdelta_wTOF[id].data(), "With TOF", kTH2F, {pAxis, deltaAxis});
    histos.add(hdelta_pt_pos_wTOF[id].data(), "With TOF Positive", kTH2F, {ptAxis, deltaAxis});
    histos.add(hdelta_pt_neg_wTOF[id].data(), "With TOF Negative", kTH2F, {ptAxis, deltaAxis});
    histos.add(hexpsigma_wTOF[id].data(), "With TOF", kTH2F, {pAxis, expSigmaAxis});
    histos.add(hnsigma_wTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {pAxis, nSigmaAxis});

    HistogramConfigSpec particleSparseHists_wTOF{HistType::kTHnSparseF, {sparseMomentumAxis, sparseEtaAxis, sparseNSigmaAxis, tpcnclsAxis}};
    if (fillTHnSparses) {
      histos.add(hnsigma_p_eta_Ncl_wTOF[id].data(), Form("With TOF %s", axisTitle), particleSparseHists_wTOF);
    }

    if (splitSignalPerCharge) {
      histos.add(hnsigma_pt_wTOF[id].data(), Form("With TOF %s", axisTitle), kTH3F, {ptAxis, nSigmaAxis, chargeAxis});
      histos.add(hsignal_wTOF[id].data(), "With TOF", kTH3F, {pAxis, dedxAxis, chargeAxis});
    } else {
      histos.add(hnsigma_pt_wTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hsignal_wTOF[id].data(), "With TOF", kTH2F, {pAxis, dedxAxis});
    }
  }

  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec etaAxis{etaBins, "#it{#eta}"};
    const AxisSpec phiAxis{phiBins, "#it{#phi}"};
    const AxisSpec lAxis{trackLengthBins, "Track length (cm)"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T}/|Z| (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p}/|Z| (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
    }
    const AxisSpec dedxAxis{dEdxBins, "d#it{E}/d#it{x} Arb. units"};
    const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};

    // Event properties
    auto h = histos.add<TH1>("event/evsel", "", kTH1D, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed vtx Z");

    h = histos.add<TH1>("event/trackselection", "", kTH1D, {{10, 0.5, 10.5, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "isGlobalTrack");
    h->GetXaxis()->SetBinLabel(3, "hasITS");
    h->GetXaxis()->SetBinLabel(4, "hasTPC");
    h->GetXaxis()->SetBinLabel(5, Form("tpcNClsFound > %i", minTPCNcls.value));

    histos.add("event/vertexz", "", kTH1D, {vtxZAxis});
    h = histos.add<TH1>("event/particlehypo", "", kTH1D, {{10, 0, 10, "PID in tracking"}});
    for (int i = 0; i < 9; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, PID::getName(i));
    }
    if (enableDeDxPlot) {
      if (splitSignalPerCharge) {
        histos.add("event/tpcsignal", "", kTH3F, {pAxis, dedxAxis, chargeAxis});
        histos.add("event/tpcsignalvspt", "", kTH3F, {ptAxis, dedxAxis, chargeAxis});
      } else {
        histos.add("event/tpcsignal", "", kTH2F, {pAxis, dedxAxis});
        histos.add("event/tpcsignalvspt", "", kTH2F, {ptAxis, dedxAxis});
      }
    }
    histos.add("event/eta", "", kTH1D, {etaAxis});
    histos.add("event/phi", "", kTH1D, {phiAxis});
    histos.add("event/etaphi", "", kTH2F, {etaAxis, phiAxis});
    histos.add("event/length", "", kTH1D, {lAxis});
    histos.add("event/pt", "", kTH1D, {ptAxis});
    histos.add("event/p", "", kTH1D, {pAxis});

    static_for<0, 8>([&](auto i) {
      initPerParticle<i>(pAxis, ptAxis, dedxAxis, chargeAxis);
    });
    LOG(info) << "QA PID TPC histograms:";
    histos.print();
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& /*tracks*/)
  {

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 1);
    }
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 2);
    }

    if (std::abs(collision.posZ()) > 10.f) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 3);
      histos.fill(HIST("event/vertexz"), collision.posZ());
    }
    return true;
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isTrackSelected(const CollisionType& /*collision*/, const TrackType& track)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 1.f);
    }
    if (!track.isGlobalTrack()) { // Skipping non global tracks
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 2.f);
    }
    if (!track.hasITS()) { // Skipping tracks without ITS
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 3.f);
    }
    if (!track.hasTPC()) { // Skipping tracks without TPC
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 4.f);
    }
    if (track.tpcNClsFound() < minTPCNcls) { // Skipping tracks without enough TPC clusters
      return false;
    }

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 5.f);
      histos.fill(HIST("event/particlehypo"), track.pidForTracking());
      if (enableDeDxPlot) {
        if (splitSignalPerCharge) {
          histos.fill(HIST("event/tpcsignal"), track.tpcInnerParam(), track.tpcSignal(), track.sign());
          histos.fill(HIST("event/tpcsignalvspt"), track.pt(), track.tpcSignal(), track.sign());
        } else {
          histos.fill(HIST("event/tpcsignal"), track.tpcInnerParam(), track.tpcSignal());
          histos.fill(HIST("event/tpcsignalvspt"), track.pt(), track.tpcSignal());
        }
      }
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/phi"), track.phi());
      histos.fill(HIST("event/etaphi"), track.eta(), track.phi());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
      // histos.fill(HIST("event/ptreso"), track.p(), track.sigma1Pt() * track.pt() * track.pt());
    }
    return true;
  }

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (o2::aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (o2::aod::evsel::sel8 == true));
  Filter trackFilter = ((trackSelection.node() == 0) ||
                        ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                        ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                        ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                        ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                        ((trackSelection.node() == 5) && requireInAcceptanceTracksInFilter()));
  using CollisionCandidate = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  void process(CollisionCandidate const& collision,
               soa::Filtered<TrackCandidates> const& tracks)
  {
    isEventSelected<true>(collision, tracks);
    for (const auto& t : tracks) {
      isTrackSelected<true>(collision, t);
    }
  }

  template <o2::track::PID::ID id, bool fillFullHistograms,
            bool fillWithTOFHistograms,
            typename TrackType>
  void processSingleParticle(CollisionCandidate const& collision,
                             TrackType const& tracks)
  {
    if (!isEventSelected<false>(collision, tracks)) {
      return;
    }

    for (const auto& t : tracks) {
      if (!isTrackSelected<false>(collision, t)) {
        continue;
      }

      if (applyRapidityCut) {
        if (std::abs(t.rapidity(PID::getMass(id))) > 0.5) {
          continue;
        }
      }

      const auto nsigma = o2::aod::pidutils::tpcNSigma<id>(t);
      histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
      if (splitSignalPerCharge) {
        histos.fill(HIST(hnsigma_pt[id]), t.pt(), nsigma, t.sign());
      } else {
        histos.fill(HIST(hnsigma_pt[id]), t.pt(), nsigma);
      }

      if constexpr (fillFullHistograms) {
        const auto& diff = o2::aod::pidutils::tpcExpSignalDiff<id>(t);
        // Fill histograms
        histos.fill(HIST(hexpected[id]), t.tpcInnerParam(), t.tpcSignal() - diff);
        histos.fill(HIST(hdelta[id]), t.tpcInnerParam(), diff);

        if (fillTHnSparses) {
          histos.fill(HIST(hnsigma_p_eta_Ncl[id]), t.p(), t.eta(), nsigma, t.tpcNClsFindable());
        }

        if (splitSignalPerCharge) {
          histos.fill(HIST(hdelta_pt[id]), t.pt(), diff, t.sign());
        } else {
          histos.fill(HIST(hdelta_pt[id]), t.pt(), diff);
        }
        histos.fill(HIST(hexpsigma[id]), t.tpcInnerParam(), o2::aod::pidutils::tpcExpSigma<id>(t));
        if constexpr (fillWithTOFHistograms) {
          if (std::abs(o2::aod::pidutils::tofNSigma<id>(t)) < 3.f) {
            histos.fill(HIST(hexpected_wTOF[id]), t.tpcInnerParam(), t.tpcSignal() - diff);
            histos.fill(HIST(hdelta_wTOF[id]), t.tpcInnerParam(), diff);
            histos.fill(HIST(hexpsigma_wTOF[id]), t.p(), o2::aod::pidutils::tpcExpSigma<id>(t));
          }
        }
      }
      if constexpr (fillWithTOFHistograms) { // Filling nsigma (common to full and tiny)
        const auto& nsigmatof = o2::aod::pidutils::tofNSigma<id>(t);
        if (std::abs(nsigmatof) < 3.f) {
          histos.fill(HIST(hnsigma_wTOF[id]), t.p(), nsigma);
          if (fillTHnSparses) {
            histos.fill(HIST(hnsigma_p_eta_Ncl_wTOF[id]), t.p(), t.eta(), nsigma, t.tpcNClsFindable());
          }
          if (splitSignalPerCharge) {
            histos.fill(HIST(hnsigma_pt_wTOF[id]), t.pt(), nsigma, t.sign());
            histos.fill(HIST(hsignal_wTOF[id]), t.tpcInnerParam(), t.tpcSignal(), t.sign());
          } else {
            histos.fill(HIST(hnsigma_pt_wTOF[id]), t.pt(), nsigma);
            histos.fill(HIST(hsignal_wTOF[id]), t.tpcInnerParam(), t.tpcSignal());
          }
          // histos.fill(HIST("event/signedtpcsignal"), t.tpcInnerParam() * t.sign(), t.tpcSignal());
        }
      }
    }
  }

  // QA of nsigma only tables
#define makeProcessFunction(inputPid, particleId)                                             \
  void process##particleId(CollisionCandidate const& collision,                               \
                           soa::Filtered<soa::Join<TrackCandidates, inputPid>> const& tracks) \
  {                                                                                           \
    processSingleParticle<PID::particleId, false, false>(collision, tracks);                  \
  }                                                                                           \
  PROCESS_SWITCH(tpcPidQa, process##particleId, Form("Process for the %s hypothesis for TPC NSigma QA", #particleId), false);

  makeProcessFunction(aod::pidTPCEl, Electron);
  makeProcessFunction(aod::pidTPCMu, Muon);
  makeProcessFunction(aod::pidTPCPi, Pion);
  makeProcessFunction(aod::pidTPCKa, Kaon);
  makeProcessFunction(aod::pidTPCPr, Proton);
  makeProcessFunction(aod::pidTPCDe, Deuteron);
  makeProcessFunction(aod::pidTPCTr, Triton);
  makeProcessFunction(aod::pidTPCHe, Helium3);
  makeProcessFunction(aod::pidTPCAl, Alpha);
#undef makeProcessFunction

// QA of full tables
#define makeProcessFunction(inputPid, particleId)                                                 \
  void processFull##particleId(CollisionCandidate const& collision,                               \
                               soa::Filtered<soa::Join<TrackCandidates, inputPid>> const& tracks) \
  {                                                                                               \
    processSingleParticle<PID::particleId, true, false>(collision, tracks);                       \
  }                                                                                               \
  PROCESS_SWITCH(tpcPidQa, processFull##particleId, Form("Process for the %s hypothesis for full TPC PID QA", #particleId), false);

  makeProcessFunction(aod::pidTPCFullEl, Electron);
  makeProcessFunction(aod::pidTPCFullMu, Muon);
  makeProcessFunction(aod::pidTPCFullPi, Pion);
  makeProcessFunction(aod::pidTPCFullKa, Kaon);
  makeProcessFunction(aod::pidTPCFullPr, Proton);
  makeProcessFunction(aod::pidTPCFullDe, Deuteron);
  makeProcessFunction(aod::pidTPCFullTr, Triton);
  makeProcessFunction(aod::pidTPCFullHe, Helium3);
  makeProcessFunction(aod::pidTPCFullAl, Alpha);
#undef makeProcessFunction

  // QA of full tables with TOF information
#define makeProcessFunction(inputPid, inputPidTOF, particleId)                                                        \
  void processFullWithTOF##particleId(CollisionCandidate const& collision,                                            \
                                      soa::Filtered<soa::Join<TrackCandidates, inputPid, inputPidTOF>> const& tracks) \
  {                                                                                                                   \
    processSingleParticle<PID::particleId, true, true>(collision, tracks);                                            \
  }                                                                                                                   \
  PROCESS_SWITCH(tpcPidQa, processFullWithTOF##particleId, Form("Process for the %s hypothesis for full TPC PID QA with the TOF info added", #particleId), false);

  makeProcessFunction(aod::pidTPCFullEl, aod::pidTOFFullEl, Electron);
  makeProcessFunction(aod::pidTPCFullMu, aod::pidTOFFullMu, Muon);
  makeProcessFunction(aod::pidTPCFullPi, aod::pidTOFFullPi, Pion);
  makeProcessFunction(aod::pidTPCFullKa, aod::pidTOFFullKa, Kaon);
  makeProcessFunction(aod::pidTPCFullPr, aod::pidTOFFullPr, Proton);
  makeProcessFunction(aod::pidTPCFullDe, aod::pidTOFFullDe, Deuteron);
  makeProcessFunction(aod::pidTPCFullTr, aod::pidTOFFullTr, Triton);
  makeProcessFunction(aod::pidTPCFullHe, aod::pidTOFFullHe, Helium3);
  makeProcessFunction(aod::pidTPCFullAl, aod::pidTOFFullAl, Alpha);
#undef makeProcessFunction
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tpcPidQa>(cfgc)};
}
