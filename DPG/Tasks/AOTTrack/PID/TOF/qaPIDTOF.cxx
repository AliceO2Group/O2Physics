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
/// \file   qaPIDTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Implementation for QA tasks of the TOF PID quantities
///

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the TOF QA plots
struct tofPidQa {
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

  static constexpr std::string_view hdelta_etaphi[Np] = {"delta/etaphi/El", "delta/etaphi/Mu", "delta/etaphi/Pi",
                                                         "delta/etaphi/Ka", "delta/etaphi/Pr", "delta/etaphi/De",
                                                         "delta/etaphi/Tr", "delta/etaphi/He", "delta/etaphi/Al"};

  // Ev. Time fill
  static constexpr std::string_view hdelta_evtime_fill[Np] = {"delta/evtime/fill/El", "delta/evtime/fill/Mu", "delta/evtime/fill/Pi",
                                                              "delta/evtime/fill/Ka", "delta/evtime/fill/Pr", "delta/evtime/fill/De",
                                                              "delta/evtime/fill/Tr", "delta/evtime/fill/He", "delta/evtime/fill/Al"};
  static constexpr std::string_view hdelta_pt_evtime_fill[Np] = {"delta/pt/evtime/fill/El", "delta/pt/evtime/fill/Mu", "delta/pt/evtime/fill/Pi",
                                                                 "delta/pt/evtime/fill/Ka", "delta/pt/evtime/fill/Pr", "delta/pt/evtime/fill/De",
                                                                 "delta/pt/evtime/fill/Tr", "delta/pt/evtime/fill/He", "delta/pt/evtime/fill/Al"};
  // Ev. Time TOF
  static constexpr std::string_view hdelta_evtime_tof[Np] = {"delta/evtime/tof/El", "delta/evtime/tof/Mu", "delta/evtime/tof/Pi",
                                                             "delta/evtime/tof/Ka", "delta/evtime/tof/Pr", "delta/evtime/tof/De",
                                                             "delta/evtime/tof/Tr", "delta/evtime/tof/He", "delta/evtime/tof/Al"};
  static constexpr std::string_view hdelta_pt_evtime_tof[Np] = {"delta/pt/evtime/tof/El", "delta/pt/evtime/tof/Mu", "delta/pt/evtime/tof/Pi",
                                                                "delta/pt/evtime/tof/Ka", "delta/pt/evtime/tof/Pr", "delta/pt/evtime/tof/De",
                                                                "delta/pt/evtime/tof/Tr", "delta/pt/evtime/tof/He", "delta/pt/evtime/tof/Al"};
  // Ev. Time FT0
  static constexpr std::string_view hdelta_evtime_ft0[Np] = {"delta/evtime/ft0/El", "delta/evtime/ft0/Mu", "delta/evtime/ft0/Pi",
                                                             "delta/evtime/ft0/Ka", "delta/evtime/ft0/Pr", "delta/evtime/ft0/De",
                                                             "delta/evtime/ft0/Tr", "delta/evtime/ft0/He", "delta/evtime/ft0/Al"};
  static constexpr std::string_view hdelta_pt_evtime_ft0[Np] = {"delta/pt/evtime/ft0/El", "delta/pt/evtime/ft0/Mu", "delta/pt/evtime/ft0/Pi",
                                                                "delta/pt/evtime/ft0/Ka", "delta/pt/evtime/ft0/Pr", "delta/pt/evtime/ft0/De",
                                                                "delta/pt/evtime/ft0/Tr", "delta/pt/evtime/ft0/He", "delta/pt/evtime/ft0/Al"};
  // Ev. Time TOF+FT0
  static constexpr std::string_view hdelta_evtime_tofft0[Np] = {"delta/evtime/tofft0/El", "delta/evtime/tofft0/Mu", "delta/evtime/tofft0/Pi",
                                                                "delta/evtime/tofft0/Ka", "delta/evtime/tofft0/Pr", "delta/evtime/tofft0/De",
                                                                "delta/evtime/tofft0/Tr", "delta/evtime/tofft0/He", "delta/evtime/tofft0/Al"};
  static constexpr std::string_view hdelta_pt_evtime_tofft0[Np] = {"delta/pt/evtime/tofft0/El", "delta/pt/evtime/tofft0/Mu", "delta/pt/evtime/tofft0/Pi",
                                                                   "delta/pt/evtime/tofft0/Ka", "delta/pt/evtime/tofft0/Pr", "delta/pt/evtime/tofft0/De",
                                                                   "delta/pt/evtime/tofft0/Tr", "delta/pt/evtime/tofft0/He", "delta/pt/evtime/tofft0/Al"};
  static constexpr std::string_view hexpsigma[Np] = {"expsigma/El", "expsigma/Mu", "expsigma/Pi",
                                                     "expsigma/Ka", "expsigma/Pr", "expsigma/De",
                                                     "expsigma/Tr", "expsigma/He", "expsigma/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hnsigma_pt[Np] = {"nsigma/pt/El", "nsigma/pt/Mu", "nsigma/pt/Pi",
                                                      "nsigma/pt/Ka", "nsigma/pt/Pr", "nsigma/pt/De",
                                                      "nsigma/pt/Tr", "nsigma/pt/He", "nsigma/pt/Al"};

  // Ev. Time fill
  static constexpr std::string_view hnsigma_evtime_fill[Np] = {"nsigma/evtime/fill/El", "nsigma/evtime/fill/Mu", "nsigma/evtime/fill/Pi",
                                                               "nsigma/evtime/fill/Ka", "nsigma/evtime/fill/Pr", "nsigma/evtime/fill/De",
                                                               "nsigma/evtime/fill/Tr", "nsigma/evtime/fill/He", "nsigma/evtime/fill/Al"};
  static constexpr std::string_view hnsigma_pt_evtime_fill[Np] = {"nsigma/pt/evtime/fill/El", "nsigma/pt/evtime/fill/Mu", "nsigma/pt/evtime/fill/Pi",
                                                                  "nsigma/pt/evtime/fill/Ka", "nsigma/pt/evtime/fill/Pr", "nsigma/pt/evtime/fill/De",
                                                                  "nsigma/pt/evtime/fill/Tr", "nsigma/pt/evtime/fill/He", "nsigma/pt/evtime/fill/Al"};
  // Ev. Time TOF
  static constexpr std::string_view hnsigma_evtime_tof[Np] = {"nsigma/evtime/tof/El", "nsigma/evtime/tof/Mu", "nsigma/evtime/tof/Pi",
                                                              "nsigma/evtime/tof/Ka", "nsigma/evtime/tof/Pr", "nsigma/evtime/tof/De",
                                                              "nsigma/evtime/tof/Tr", "nsigma/evtime/tof/He", "nsigma/evtime/tof/Al"};
  static constexpr std::string_view hnsigma_pt_evtime_tof[Np] = {"nsigma/pt/evtime/tof/El", "nsigma/pt/evtime/tof/Mu", "nsigma/pt/evtime/tof/Pi",
                                                                 "nsigma/pt/evtime/tof/Ka", "nsigma/pt/evtime/tof/Pr", "nsigma/pt/evtime/tof/De",
                                                                 "nsigma/pt/evtime/tof/Tr", "nsigma/pt/evtime/tof/He", "nsigma/pt/evtime/tof/Al"};
  // Ev. Time FT0
  static constexpr std::string_view hnsigma_evtime_ft0[Np] = {"nsigma/evtime/ft0/El", "nsigma/evtime/ft0/Mu", "nsigma/evtime/ft0/Pi",
                                                              "nsigma/evtime/ft0/Ka", "nsigma/evtime/ft0/Pr", "nsigma/evtime/ft0/De",
                                                              "nsigma/evtime/ft0/Tr", "nsigma/evtime/ft0/He", "nsigma/evtime/ft0/Al"};
  static constexpr std::string_view hnsigma_pt_evtime_ft0[Np] = {"nsigma/pt/evtime/ft0/El", "nsigma/pt/evtime/ft0/Mu", "nsigma/pt/evtime/ft0/Pi",
                                                                 "nsigma/pt/evtime/ft0/Ka", "nsigma/pt/evtime/ft0/Pr", "nsigma/pt/evtime/ft0/De",
                                                                 "nsigma/pt/evtime/ft0/Tr", "nsigma/pt/evtime/ft0/He", "nsigma/pt/evtime/ft0/Al"};
  // Ev. Time TOF+FT0
  static constexpr std::string_view hnsigma_evtime_tofft0[Np] = {"nsigma/evtime/tofft0/El", "nsigma/evtime/tofft0/Mu", "nsigma/evtime/tofft0/Pi",
                                                                 "nsigma/evtime/tofft0/Ka", "nsigma/evtime/tofft0/Pr", "nsigma/evtime/tofft0/De",
                                                                 "nsigma/evtime/tofft0/Tr", "nsigma/evtime/tofft0/He", "nsigma/evtime/tofft0/Al"};
  static constexpr std::string_view hnsigma_pt_evtime_tofft0[Np] = {"nsigma/pt/evtime/tofft0/El", "nsigma/pt/evtime/tofft0/Mu", "nsigma/pt/evtime/tofft0/Pi",
                                                                    "nsigma/pt/evtime/tofft0/Ka", "nsigma/pt/evtime/tofft0/Pr", "nsigma/pt/evtime/tofft0/De",
                                                                    "nsigma/pt/evtime/tofft0/Tr", "nsigma/pt/evtime/tofft0/He", "nsigma/pt/evtime/tofft0/Al"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.1f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 5.f, "Maximum momentum in range"};
  ConfigurableAxis etaBins{"etaBins", {100, -1.f, 1.f}, "Binning in eta"};
  ConfigurableAxis phiBins{"phiBins", {100, 0, TMath::TwoPi()}, "Binning in eta"};
  ConfigurableAxis trackLengthBins{"trackLengthBins", {100, 0, 1000.f}, "Binning in track length plot"};
  ConfigurableAxis deltaBins{"deltaBins", {200, -1000.f, 1000.f}, "Binning in Delta (T-Texp-T0)"};
  ConfigurableAxis expSigmaBins{"expSigmaBins", {200, 0.f, 200.f}, "Binning in expected Sigma"};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {401, -10.025f, 10.025f}, "Binning in NSigma"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> applyRapidityCut{"applyRapidityCut", false, "Flag to apply rapidity cut"};
  Configurable<bool> enableEvTimeSplitting{"enableEvTimeSplitting", false, "Flag to enable histograms splitting depending on the Event Time used"};
  Configurable<bool> produceDeltaTEtaPhiMap{"produceDeltaTEtaPhiMap", false, "Produces the map of the delta time as a function of eta and phi"};
  Configurable<float> ptDeltaTEtaPhiMapMin{"ptDeltaTEtaPhiMapMin", 1.45f, "Threshold in pT to build the map of the delta time as a function of eta and phi"};
  Configurable<float> ptDeltaTEtaPhiMapMax{"ptDeltaTEtaPhiMapMax", 1.55f, "Threshold in pT to build the map of the delta time as a function of eta and phi"};
  Configurable<bool> splitSignalPerCharge{"splitSignalPerCharge", true, "Split the signal per charge (reduces memory footprint if off)"};
  Configurable<int> enableVsMomentumHistograms{"enableVsMomentumHistograms", 0, "1: Enables plots vs momentum instead of just pT 2: Enables plots vs momentum vs eta instead of just pT (reduces memory footprint if off)"};
  Configurable<bool> requireGoodMatchTracks{"requireGoodMatchTracks", false, "Require good match tracks"};
  Configurable<float> pvContributorsMin{"pvContributorsMin", -10, "Minimum pvContributors"};
  Configurable<float> pvContributorsMax{"pvContributorsMax", 10000, "Maximum pvContributors"};

  template <o2::track::PID::ID id>
  void initPerParticle(const AxisSpec& pAxis,
                       const AxisSpec& ptAxis,
                       const AxisSpec& etaAxis,
                       const AxisSpec& phiAxis,
                       const AxisSpec& chargeAxis)
  {
    static_assert(id >= 0 && id <= PID::Alpha && "Particle index outside limits");
    bool enableFullHistos = false;
    int enabledProcesses = 0;
    switch (id) { // Skipping disabled particles
#define particleCase(particleId)                                 \
  case PID::particleId:                                          \
    if (!doprocess##particleId && !doprocessFull##particleId) {  \
      return;                                                    \
    }                                                            \
    if (doprocess##particleId) {                                 \
      enabledProcesses++;                                        \
    }                                                            \
    if (doprocessFull##particleId) {                             \
      enableFullHistos = true;                                   \
      enabledProcesses++;                                        \
    }                                                            \
    LOGF(info, "Enabled TOF QA for %s %s", #particleId, pT[id]); \
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
    const char* axisTitle = Form("N_{#sigma}^{TOF}(%s)", pT[id]);
    const AxisSpec nSigmaAxis{nSigmaBins, axisTitle};
    histos.add(hnsigma[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    if (splitSignalPerCharge) {
      histos.add(hnsigma_pt[id].data(), axisTitle, kTH3F, {ptAxis, nSigmaAxis, chargeAxis});
    } else {
      histos.add(hnsigma_pt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    }
    if (enableEvTimeSplitting) {
      histos.add(hnsigma_evtime_fill[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      histos.add(hnsigma_evtime_tof[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      histos.add(hnsigma_evtime_ft0[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      histos.add(hnsigma_evtime_tofft0[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});

      if (splitSignalPerCharge) {
        histos.add(hnsigma_pt_evtime_fill[id].data(), axisTitle, kTH3F, {ptAxis, nSigmaAxis, chargeAxis});
        histos.add(hnsigma_pt_evtime_tof[id].data(), axisTitle, kTH3F, {ptAxis, nSigmaAxis, chargeAxis});
        histos.add(hnsigma_pt_evtime_ft0[id].data(), axisTitle, kTH3F, {ptAxis, nSigmaAxis, chargeAxis});
        histos.add(hnsigma_pt_evtime_tofft0[id].data(), axisTitle, kTH3F, {ptAxis, nSigmaAxis, chargeAxis});
      } else {
        histos.add(hnsigma_pt_evtime_fill[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
        histos.add(hnsigma_pt_evtime_tof[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
        histos.add(hnsigma_pt_evtime_ft0[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
        histos.add(hnsigma_pt_evtime_tofft0[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
      }
    }

    if (!enableFullHistos) { // Enabling only NSigma for tiny tables
      return;
    }

    // Exp signal
    const AxisSpec expAxis{1000, 0, 2e6, Form("t_{exp}(%s) (ps)", pT[id])};
    histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{deltaBins, Form("t-t_{ev}-t_{exp}(%s) (ps)", pT[id])};
    axisTitle = Form("#Delta^{TOF}(%s)", pT[id]);
    histos.add(hdelta[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
    if (splitSignalPerCharge) {
      histos.add(hdelta_pt[id].data(), axisTitle, kTH3F, {ptAxis, deltaAxis, chargeAxis});
    } else {
      histos.add(hdelta_pt[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
    }
    if (produceDeltaTEtaPhiMap) {
      histos.add(hdelta_etaphi[id].data(), Form("%s, %.2f < #it{p}_{T} < %.2f", axisTitle, ptDeltaTEtaPhiMapMin.value, ptDeltaTEtaPhiMapMax.value), kTH3F, {etaAxis, phiAxis, deltaAxis});
    }

    // Exp Sigma
    const AxisSpec expSigmaAxis{expSigmaBins, Form("Exp_{#sigma}^{TOF}(%s) (ps)", pT[id])};
    histos.add(hexpsigma[id].data(), "", kTH2F, {pAxis, expSigmaAxis});

    if (!enableEvTimeSplitting) { // Returning if the plots with the different event time are not reqested
      return;
    }

    if (enableVsMomentumHistograms == 1) {
      histos.add(hdelta_evtime_fill[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
      histos.add(hdelta_evtime_tof[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
      histos.add(hdelta_evtime_ft0[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
      histos.add(hdelta_evtime_tofft0[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
    } else if (enableVsMomentumHistograms == 2) {
      histos.add(hdelta_evtime_fill[id].data(), axisTitle, kTH3F, {pAxis, etaAxis, deltaAxis});
      histos.add(hdelta_evtime_tof[id].data(), axisTitle, kTH3F, {pAxis, etaAxis, deltaAxis});
      histos.add(hdelta_evtime_ft0[id].data(), axisTitle, kTH3F, {pAxis, etaAxis, deltaAxis});
      histos.add(hdelta_evtime_tofft0[id].data(), axisTitle, kTH3F, {pAxis, etaAxis, deltaAxis});
    }

    if (splitSignalPerCharge) {
      histos.add(hdelta_pt_evtime_fill[id].data(), axisTitle, kTH3F, {ptAxis, deltaAxis, chargeAxis});
      histos.add(hdelta_pt_evtime_tof[id].data(), axisTitle, kTH3F, {ptAxis, deltaAxis, chargeAxis});
      histos.add(hdelta_pt_evtime_ft0[id].data(), axisTitle, kTH3F, {ptAxis, deltaAxis, chargeAxis});
      histos.add(hdelta_pt_evtime_tofft0[id].data(), axisTitle, kTH3F, {ptAxis, deltaAxis, chargeAxis});
    } else {
      histos.add(hdelta_pt_evtime_fill[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
      histos.add(hdelta_pt_evtime_tof[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
      histos.add(hdelta_pt_evtime_ft0[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
      histos.add(hdelta_pt_evtime_tofft0[id].data(), axisTitle, kTH2F, {ptAxis, deltaAxis});
    }
  }

  void init(o2::framework::InitContext&)
  {
    const AxisSpec multAxis{100, 0, 100, "TOF multiplicity"};
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec contributorsAxis{100, 0, 1000, "PV contributors"};
    const AxisSpec etaAxis{etaBins, "#it{#eta}"};
    const AxisSpec phiAxis{phiBins, "#it{#phi}"};
    const AxisSpec colTimeAxis{100, -2000, 2000, "Collision time (ps)"};
    const AxisSpec colTimeResoAxis{100, 0, 1000, "#sigma_{Collision time} (ps)"};
    const AxisSpec lAxis{trackLengthBins, "Track length (cm)"};
    const AxisSpec ptResoAxis{100, 0, 0.1, "#sigma_{#it{p}_{T}}"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T}/|Z| (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p}/|Z| (GeV/#it{c})"};
    AxisSpec pExpAxis{nBinsP, minP, maxP, "#it{p}_{Exp. TOF} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
      pExpAxis.makeLogarithmic();
    }
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal (ps)"};
    const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};

    // Event properties
    auto h = histos.add<TH1>("event/evsel", "", kTH1D, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed vtx Z");
    h->GetXaxis()->SetBinLabel(4, Form("Passed pvContributorsMin %f", pvContributorsMin.value));
    h->GetXaxis()->SetBinLabel(5, Form("Passed pvContributorsMax %f", pvContributorsMax.value));

    h = histos.add<TH1>("event/trackselection", "", kTH1D, {{10, 0.5, 10.5, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "isGlobalTrack");
    h->GetXaxis()->SetBinLabel(3, "hasITS");
    h->GetXaxis()->SetBinLabel(4, "hasTPC");
    h->GetXaxis()->SetBinLabel(5, "hasTOF");
    h->GetXaxis()->SetBinLabel(6, "goodTOFMatch");

    histos.add("event/pvcontributors", "", kTH1D, {contributorsAxis});
    histos.add("event/vertexz", "", kTH1D, {vtxZAxis});
    h = histos.add<TH1>("event/particlehypo", "", kTH1D, {{10, 0, 10, "PID in tracking"}});
    for (int i = 0; i < 9; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, PID::getName(i));
    }

    histos.add("event/evtime/colltime", "collisionTime()", kTH1D, {colTimeAxis});
    histos.add("event/evtime/colltimereso", "collisionTimeRes()", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/undef", "Undefined event time", kTH1D, {colTimeAxis});
    histos.add("event/evtime/undefreso", "Undefined event time reso.", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/avail", "Available event time", kTH1D, {colTimeAxis});
    histos.add("event/evtime/availreso", "Available event time reso.", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/ft0tof", "FT0+TOF event time", kTH1D, {colTimeAxis});
    histos.add("event/evtime/ft0tofreso", "FT0+TOF event time reso.", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/tof", "TOF event time", kTH1D, {colTimeAxis});
    histos.add("event/evtime/tofreso", "TOF event time reso.", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/evtime/ft0", "FT0 event time", kTH1D, {colTimeAxis});
    histos.add("event/evtime/ft0reso", "FT0 event time reso.", kTH2F, {multAxis, colTimeResoAxis});

    histos.add("event/tofsignal", "TOF signal", kTH2F, {pAxis, tofAxis});
    histos.add("event/tofsignalunassigned", "TOF signal (unassigned tracks)", kTH2F, {pAxis, tofAxis});
    histos.add("event/pexp", "", kTH2F, {pAxis, pExpAxis});
    histos.add("event/eta", "", kTH1D, {etaAxis});
    histos.add("event/phi", "", kTH1D, {phiAxis});
    histos.add("event/etaphi", "", kTH2F, {etaAxis, phiAxis});
    histos.add("event/length", "", kTH1D, {lAxis});
    histos.add("event/pt", "", kTH1D, {ptAxis});
    histos.add("event/p", "", kTH1D, {pAxis});
    // histos.add("event/ptreso", "", kTH2F, {pAxis, ptResoAxis});

    static_for<0, 8>([&](auto i) {
      initPerParticle<i>(pAxis, ptAxis, etaAxis, phiAxis, chargeAxis);
    });
    LOG(info) << "QA PID TOF histograms:";
    histos.print();
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& tracks)
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

    int tofmult = 0;
    float evtime = 0.f;
    float evtimereso = 0.f;
    int evtimeflag = 0;

    if constexpr (fillHistograms) {
      for (auto t : tracks) {
        if (!t.hasTOF()) { // Skipping tracks without TOF
          continue;
        }
        tofmult++;
        evtime = t.tofEvTime();
        evtimereso = t.tofEvTimeErr();
        evtimeflag = 0;
        if (t.isEvTimeDefined()) {
          evtimeflag = 1;
        }
        if (t.isEvTimeTOF() && t.isEvTimeT0AC()) {
          evtimeflag = 2;
        } else if (t.isEvTimeTOF()) {
          evtimeflag = 3;
        } else if (t.isEvTimeT0AC()) {
          evtimeflag = 4;
        }
      }
    }
    if (std::abs(collision.posZ()) > 10.f) {
      return false;
    }
    // Count the number of contributors
    int pvContributors = 0;
    for (const auto& trk : tracks) {
      if (trk.isPVContributor()) {
        pvContributors++;
      }
    }
    histos.fill(HIST("event/pvcontributors"), pvContributors);
    if (pvContributors < pvContributorsMin) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 4);
    }
    if (pvContributors > pvContributorsMax) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 5);
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 6);
      histos.fill(HIST("event/vertexz"), collision.posZ());

      histos.fill(HIST("event/evtime/colltime"), collision.collisionTime() * 1000.f);
      histos.fill(HIST("event/evtime/colltimereso"), tofmult, collision.collisionTimeRes() * 1000.f);

      switch (evtimeflag) {
        case 0:
          histos.fill(HIST("event/evtime/undef"), evtime);
          histos.fill(HIST("event/evtime/undefreso"), tofmult, evtimereso);
          break;
        case 1:
          histos.fill(HIST("event/evtime/avail"), evtime);
          histos.fill(HIST("event/evtime/availreso"), tofmult, evtimereso);
          break;
        case 2:
          histos.fill(HIST("event/evtime/ft0tof"), evtime);
          histos.fill(HIST("event/evtime/ft0tofreso"), tofmult, evtimereso);
          break;
        case 3:
          histos.fill(HIST("event/evtime/tof"), evtime);
          histos.fill(HIST("event/evtime/tofreso"), tofmult, evtimereso);
          break;
        case 4:
          histos.fill(HIST("event/evtime/tof"), evtime);
          histos.fill(HIST("event/evtime/tofreso"), tofmult, evtimereso);
          break;
        default:
          LOG(fatal) << "Unrecognized Event time flag";
          break;
      }
    }
    return true;
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isTrackSelected(const CollisionType&, const TrackType& track)
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
    if (!track.hasTOF()) { // Skipping tracks without TOF
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 5.f);
    }
    if (requireGoodMatchTracks.value && !track.goodTOFMatch()) { // Skipping tracks without good match
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 6.f);
      histos.fill(HIST("event/particlehypo"), track.pidForTracking());
      if (track.has_collision()) {
        histos.fill(HIST("event/tofsignal"), track.p(), track.tofSignal());
      } else {
        histos.fill(HIST("event/tofsignalunassigned"), track.p(), track.tofSignal());
      }
      histos.fill(HIST("event/pexp"), track.p(), track.tofExpMom());
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
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireInAcceptanceTracksInFilter());
  using CollisionCandidate = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                    aod::pidEvTimeFlags, aod::TOFSignal, aod::TOFEvTime,
                                    aod::pidTOFFlags>;

  void process(CollisionCandidate const& collision,
               soa::Filtered<TrackCandidates> const& tracks)
  {
    isEventSelected<true>(collision, tracks);
    for (auto t : tracks) {
      isTrackSelected<true>(collision, t);
    }
  }

  template <o2::track::PID::ID id, bool fillFullHistograms,
            typename TrackType>
  void processSingleParticle(CollisionCandidate const& collision,
                             TrackType const& tracks)
  {
    if (!isEventSelected<false>(collision, tracks)) {
      return;
    }

    for (auto t : tracks) {
      if (!isTrackSelected<false>(collision, t)) {
        continue;
      }

      if (applyRapidityCut) {
        if (std::abs(t.rapidity(PID::getMass(id))) > 0.5) {
          continue;
        }
      }

      const auto nsigma = o2::aod::pidutils::tofNSigma<id>(t);
      histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
      if (splitSignalPerCharge) {
        histos.fill(HIST(hnsigma_pt[id]), t.pt(), nsigma, t.sign());
      } else {
        histos.fill(HIST(hnsigma_pt[id]), t.pt(), nsigma);
      }
      // Filling info split per ev. time
      if (enableEvTimeSplitting) {
        if (t.isEvTimeTOF() && t.isEvTimeT0AC()) { // TOF + FT0 Ev. Time
          histos.fill(HIST(hnsigma_evtime_tofft0[id]), t.p(), nsigma);
          if (splitSignalPerCharge) {
            histos.fill(HIST(hnsigma_pt_evtime_tofft0[id]), t.pt(), nsigma, t.sign());
          } else {
            histos.fill(HIST(hnsigma_pt_evtime_tofft0[id]), t.pt(), nsigma);
          }
        } else if (t.isEvTimeT0AC()) { // FT0 Ev. Time
          histos.fill(HIST(hnsigma_evtime_ft0[id]), t.p(), nsigma);
          if (splitSignalPerCharge) {
            histos.fill(HIST(hnsigma_pt_evtime_ft0[id]), t.pt(), nsigma, t.sign());
          } else {
            histos.fill(HIST(hnsigma_pt_evtime_ft0[id]), t.pt(), nsigma);
          }
        } else if (t.isEvTimeTOF()) { // TOF Ev. Time
          histos.fill(HIST(hnsigma_evtime_tof[id]), t.p(), nsigma);
          if (splitSignalPerCharge) {
            histos.fill(HIST(hnsigma_pt_evtime_tof[id]), t.pt(), nsigma, t.sign());
          } else {
            histos.fill(HIST(hnsigma_pt_evtime_tof[id]), t.pt(), nsigma);
          }
        } else { // No Ev. Time -> Fill Ev. Time
          histos.fill(HIST(hnsigma_evtime_fill[id]), t.p(), nsigma);
          if (splitSignalPerCharge) {
            histos.fill(HIST(hnsigma_pt_evtime_fill[id]), t.pt(), nsigma, t.sign());
          } else {
            histos.fill(HIST(hnsigma_pt_evtime_fill[id]), t.pt(), nsigma);
          }
        }
      }

      if constexpr (fillFullHistograms) {
        const float& tof = t.tofSignal() - t.tofEvTime();
        const auto& diff = o2::aod::pidutils::tofExpSignalDiff<id>(t);
        // Fill histograms
        histos.fill(HIST(hexpected[id]), t.p(), tof - diff);
        histos.fill(HIST(hdelta[id]), t.p(), diff);
        if (splitSignalPerCharge) {
          histos.fill(HIST(hdelta_pt[id]), t.pt(), diff, t.sign());
        } else {
          histos.fill(HIST(hdelta_pt[id]), t.pt(), diff);
        }

        if (produceDeltaTEtaPhiMap) {
          if (t.pt() > ptDeltaTEtaPhiMapMin && t.pt() < ptDeltaTEtaPhiMapMax) {
            histos.fill(HIST(hdelta_etaphi[id]), t.eta(), t.phi(), diff);
          }
        }
        histos.fill(HIST(hexpsigma[id]), t.p(), o2::aod::pidutils::tofExpSigma<id>(t));

        // Filling info split per ev. time
        if (enableEvTimeSplitting) {
          if (t.isEvTimeTOF() && t.isEvTimeT0AC()) { // TOF + FT0 Ev. Time
            if (enableVsMomentumHistograms == 1) {
              histos.fill(HIST(hdelta_evtime_tofft0[id]), t.p(), diff);
            } else if (enableVsMomentumHistograms == 2) {
              histos.fill(HIST(hdelta_evtime_tofft0[id]), t.p(), t.eta(), diff);
            }
            if (splitSignalPerCharge) {
              histos.fill(HIST(hdelta_pt_evtime_tofft0[id]), t.pt(), diff, t.sign());
            } else {
              histos.fill(HIST(hdelta_pt_evtime_tofft0[id]), t.pt(), diff);
            }
          } else if (t.isEvTimeT0AC()) { // FT0 Ev. Time
            if (enableVsMomentumHistograms == 1) {
              histos.fill(HIST(hdelta_evtime_ft0[id]), t.p(), diff);
            } else if (enableVsMomentumHistograms == 2) {
              histos.fill(HIST(hdelta_evtime_ft0[id]), t.p(), t.eta(), diff);
            }
            if (splitSignalPerCharge) {
              histos.fill(HIST(hdelta_pt_evtime_ft0[id]), t.pt(), diff, t.sign());
            } else {
              histos.fill(HIST(hdelta_pt_evtime_ft0[id]), t.pt(), diff);
            }
          } else if (t.isEvTimeTOF()) { // TOF Ev. Time
            if (enableVsMomentumHistograms == 1) {
              histos.fill(HIST(hdelta_evtime_tof[id]), t.p(), diff);
            } else if (enableVsMomentumHistograms == 2) {
              histos.fill(HIST(hdelta_evtime_tof[id]), t.p(), t.eta(), diff);
            }
            if (splitSignalPerCharge) {
              histos.fill(HIST(hdelta_pt_evtime_tof[id]), t.pt(), diff, t.sign());
            } else {
              histos.fill(HIST(hdelta_pt_evtime_tof[id]), t.pt(), diff);
            }
          } else { // No Ev. Time -> Fill Ev. Time
            if (enableVsMomentumHistograms == 1) {
              histos.fill(HIST(hdelta_evtime_fill[id]), t.p(), diff);
            } else if (enableVsMomentumHistograms == 2) {
              histos.fill(HIST(hdelta_evtime_fill[id]), t.p(), t.eta(), diff);
            }
            if (splitSignalPerCharge) {
              histos.fill(HIST(hdelta_pt_evtime_fill[id]), t.pt(), diff, t.sign());
            } else {
              histos.fill(HIST(hdelta_pt_evtime_fill[id]), t.pt(), diff);
            }
          }
        }
      }
    }
  }

  // QA of nsigma only tables
#define makeProcessFunction(inputPid, particleId)                                             \
  void process##particleId(CollisionCandidate const& collision,                               \
                           soa::Filtered<soa::Join<TrackCandidates, inputPid>> const& tracks) \
  {                                                                                           \
    processSingleParticle<PID::particleId, false>(collision, tracks);                         \
  }                                                                                           \
  PROCESS_SWITCH(tofPidQa, process##particleId, Form("Process for the %s hypothesis for TOF NSigma QA", #particleId), false);

  makeProcessFunction(aod::pidTOFEl, Electron);
  makeProcessFunction(aod::pidTOFMu, Muon);
  makeProcessFunction(aod::pidTOFPi, Pion);
  makeProcessFunction(aod::pidTOFKa, Kaon);
  makeProcessFunction(aod::pidTOFPr, Proton);
  makeProcessFunction(aod::pidTOFDe, Deuteron);
  makeProcessFunction(aod::pidTOFTr, Triton);
  makeProcessFunction(aod::pidTOFHe, Helium3);
  makeProcessFunction(aod::pidTOFAl, Alpha);
#undef makeProcessFunction

// QA of full tables
#define makeProcessFunction(inputPid, particleId)                                                 \
  void processFull##particleId(CollisionCandidate const& collision,                               \
                               soa::Filtered<soa::Join<TrackCandidates, inputPid>> const& tracks) \
  {                                                                                               \
    processSingleParticle<PID::particleId, true>(collision, tracks);                              \
  }                                                                                               \
  PROCESS_SWITCH(tofPidQa, processFull##particleId, Form("Process for the %s hypothesis for full TOF PID QA", #particleId), false);

  makeProcessFunction(aod::pidTOFFullEl, Electron);
  makeProcessFunction(aod::pidTOFFullMu, Muon);
  makeProcessFunction(aod::pidTOFFullPi, Pion);
  makeProcessFunction(aod::pidTOFFullKa, Kaon);
  makeProcessFunction(aod::pidTOFFullPr, Proton);
  makeProcessFunction(aod::pidTOFFullDe, Deuteron);
  makeProcessFunction(aod::pidTOFFullTr, Triton);
  makeProcessFunction(aod::pidTOFFullHe, Helium3);
  makeProcessFunction(aod::pidTOFFullAl, Alpha);
#undef makeProcessFunction
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofPidQa>(cfgc)};
}
