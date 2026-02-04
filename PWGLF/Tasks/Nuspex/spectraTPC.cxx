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
/// \file   spectraTPC.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
///
/// \brief Task for the analysis of the spectra with the TPC detector
///        In addition the task makes histograms of the TPC signal with TOF selections.
///

// O2 includes
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{
    {"add-tof-histos", VariantType::Int, 0, {"Generate TPC with TOF histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

// Spectra task
struct tpcSpectra {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  Configurable<float> cfgNSigmaCut{"cfgNSigmaCut", 3, "Value of the Nsigma cut"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<bool> isRun2{"isRun2", false, "Flag to process Run 2 data"};

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  static constexpr std::string_view hp[Np] = {"p/El", "p/Mu", "p/Pi",
                                              "p/Ka", "p/Pr", "p/De",
                                              "p/Tr", "p/He", "p/Al"};
  static constexpr std::string_view hpt[Np] = {"pt/El", "pt/Mu", "pt/Pi",
                                               "pt/Ka", "pt/Pr", "pt/De",
                                               "pt/Tr", "pt/He", "pt/Al"};
  static constexpr std::string_view hdcaxy[Np] = {"dcaxy/El", "dcaxy/Mu", "dcaxy/Pi",
                                                  "dcaxy/Ka", "dcaxy/Pr", "dcaxy/De",
                                                  "dcaxy/Tr", "dcaxy/He", "dcaxy/Al"};
  static constexpr std::string_view hdcaz[Np] = {"dcaz/El", "dcaz/Mu", "dcaz/Pi",
                                                 "dcaz/Ka", "dcaz/Pr", "dcaz/De",
                                                 "dcaz/Tr", "dcaz/He", "dcaz/Al"};
  static constexpr std::string_view hdcaxyphi[Np] = {"dcaxyphi/El", "dcaxyphi/Mu", "dcaxyphi/Pi",
                                                     "dcaxyphi/Ka", "dcaxyphi/Pr", "dcaxyphi/De",
                                                     "dcaxyphi/Tr", "dcaxyphi/He", "dcaxyphi/Al"};

  TrackSelection globalTrackswoPrim; // Track without cut for primaries

  void init(o2::framework::InitContext&)
  {
    globalTrackswoPrim = getGlobalTrackSelection();
    globalTrackswoPrim.SetMaxDcaXYPtDep([](float pt) { return 3.f + pt; });
    globalTrackswoPrim.SetRequireGoldenChi2(false);
    if (!isRun2) {
      globalTrackswoPrim.SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    }

    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("event/vertexz", "", HistType::kTH1F, {vtxZAxis});
    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Ev. sel. passed");
    h->GetXaxis()->SetBinLabel(3, "posZ passed");
    h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Eta passed");
    h->GetXaxis()->SetBinLabel(3, "Quality passed");
    histos.add("p/Unselected", "Unselected", kTH1F, {pAxis});
    histos.add("pt/Unselected", "Unselected", kTH1F, {ptAxis});
    for (int i = 0; i < Np; i++) {
      histos.add(hp[i].data(), pT[i], kTH1F, {pAxis});
      histos.add(hpt[i].data(), pT[i], kTH1F, {ptAxis});
    }

    // DCAxy
    const AxisSpec dcaXyAxis{600, -3.005, 2.995, "DCA_{xy} (cm)"};
    const AxisSpec phiAxis{200, 0, 7, "#it{#varphi} (rad)"};
    const AxisSpec dcaZAxis{600, -3.005, 2.995, "DCA_{z} (cm)"};
    for (int i = 0; i < Np; i++) {
      histos.add(hdcaxy[i].data(), pT[i], kTH2F, {ptAxis, dcaXyAxis});
      histos.add(hdcaz[i].data(), pT[i], kTH2F, {ptAxis, dcaZAxis});
      histos.add(hdcaxyphi[i].data(), Form("%s -- 0.9 < #it{p}_{T} < 1.1 GeV/#it{c}", pT[i]), kTH2F, {phiAxis, dcaXyAxis});
    }
  }

  template <PID::ID id, typename T>
  void fillParticleHistos(const T& track)
  {
    const float y = TMath::ASinH(track.pt() / TMath::Sqrt(PID::getMass2(id) + track.pt() * track.pt()) * TMath::SinH(track.eta()));
    if (std::abs(y) > 0.5) {
      return;
    }
    const auto& nsigma = o2::aod::pidutils::tpcNSigma<id>(track);
    if (std::abs(nsigma) < 2) {
      histos.fill(HIST(hdcaxy[id]), track.pt(), track.dcaXY());
      histos.fill(HIST(hdcaz[id]), track.pt(), track.dcaZ());
      if (track.pt() < 1.1 && track.pt() > 0.9) {
        histos.fill(HIST(hdcaxyphi[id]), track.phi(), track.dcaXY());
      }
    }
    if (!track.isGlobalTrack()) {
      return;
    }
    if (std::abs(nsigma) > cfgNSigmaCut) {
      return;
    }
    histos.fill(HIST(hp[id]), track.p());
    histos.fill(HIST(hpt[id]), track.pt());
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                    aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                    aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                    aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                                    aod::TrackSelection>;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               TrackCandidates const& tracks)
  {
    histos.fill(HIST("evsel"), 1);
    if (isRun2 && !collision.sel7()) {
      return;

    } else if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("evsel"), 2);
    if (std::abs(collision.posZ()) > cfgCutVertex) {
      return;
    }
    histos.fill(HIST("evsel"), 3);
    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (const auto& track : tracks) {
      histos.fill(HIST("tracksel"), 1);
      if (std::abs(track.eta()) > cfgCutEta) {
        continue;
      }
      histos.fill(HIST("tracksel"), 2);
      if (!globalTrackswoPrim.IsSelected(track)) {
        continue;
      }
      histos.fill(HIST("tracksel"), 3);

      histos.fill(HIST("p/Unselected"), track.p());
      histos.fill(HIST("pt/Unselected"), track.pt());

      fillParticleHistos<PID::Electron>(track);
      fillParticleHistos<PID::Muon>(track);
      fillParticleHistos<PID::Pion>(track);
      fillParticleHistos<PID::Kaon>(track);
      fillParticleHistos<PID::Proton>(track);
      fillParticleHistos<PID::Deuteron>(track);
      fillParticleHistos<PID::Triton>(track);
      fillParticleHistos<PID::Helium3>(track);
      fillParticleHistos<PID::Alpha>(track);
    }
  } // end of the process function
}; // end of spectra task

struct tpcPidQaSignalwTof {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 200, "Number of bins for the NSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Minimum NSigma in range"};
  Configurable<float> maxNSigma{"maxNSigma", 10.f, "Maximum NSigma in range"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  static constexpr std::string_view htpcsignal[Np] = {"tpcsignal/El", "tpcsignal/Mu", "tpcsignal/Pi",
                                                      "tpcsignal/Ka", "tpcsignal/Pr", "tpcsignal/De",
                                                      "tpcsignal/Tr", "tpcsignal/He", "tpcsignal/Al"};
  static constexpr std::string_view hnsigmatpctof[Np] = {"nsigmatpctof/El", "nsigmatpctof/Mu", "nsigmatpctof/Pi",
                                                         "nsigmatpctof/Ka", "nsigmatpctof/Pr", "nsigmatpctof/De",
                                                         "nsigmatpctof/Tr", "nsigmatpctof/He", "nsigmatpctof/Al"};
  static constexpr std::string_view hdcaxy[Np] = {"dcaxy/El", "dcaxy/Mu", "dcaxy/Pi",
                                                  "dcaxy/Ka", "dcaxy/Pr", "dcaxy/De",
                                                  "dcaxy/Tr", "dcaxy/He", "dcaxy/Al"};
  static constexpr std::string_view hdcaz[Np] = {"dcaz/El", "dcaz/Mu", "dcaz/Pi",
                                                 "dcaz/Ka", "dcaz/Pr", "dcaz/De",
                                                 "dcaz/Tr", "dcaz/He", "dcaz/Al"};
  static constexpr std::string_view hdcaxyphi[Np] = {"dcaxyphi/El", "dcaxyphi/Mu", "dcaxyphi/Pi",
                                                     "dcaxyphi/Ka", "dcaxyphi/Pr", "dcaxyphi/De",
                                                     "dcaxyphi/Tr", "dcaxyphi/He", "dcaxyphi/Al"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    histos.add("event/vertexz", "", HistType::kTH1F, {vtxZAxis});
    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "posZ passed");
    h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Eta passed");
    h->GetXaxis()->SetBinLabel(3, "Quality passed");
    h->GetXaxis()->SetBinLabel(4, "TOF passed");

    AxisSpec pAxis{1000, 0.001, 20, "#it{p} (GeV/#it{c})"};
    pAxis.makeLogarithmic();
    AxisSpec ptAxis{1000, 0.001, 20, "#it{p}_{T} (GeV/#it{c})"};
    ptAxis.makeLogarithmic();
    const AxisSpec axisSignal{1000, 0, 1000, "TPC Signal"};
    const AxisSpec dcaXyAxis{600, -3.005, 2.995, "DCA_{xy} (cm)"};
    const AxisSpec phiAxis{200, 0, 7, "#it{#varphi} (rad)"};
    const AxisSpec dcaZAxis{600, -3.005, 2.995, "DCA_{z} (cm)"};

    for (int i = 0; i < Np; i++) {
      const AxisSpec nSigmaTPCAxis{nBinsNSigma, minNSigma, maxNSigma, Form("N_{#sigma}^{TPC}(%s)", pT[i])};
      const AxisSpec nSigmaTOFAxis{nBinsNSigma, minNSigma, maxNSigma, Form("N_{#sigma}^{TOF}(%s)", pT[i])};
      // TPC Signal
      histos.add(htpcsignal[i].data(), pT[i], kTH3D, {pAxis, axisSignal, nSigmaTOFAxis});
      histos.add(hnsigmatpctof[i].data(), pT[i], kTH3F, {ptAxis, nSigmaTPCAxis, nSigmaTOFAxis});
      // DCAxy
      histos.add(hdcaxy[i].data(), pT[i], kTH2F, {ptAxis, dcaXyAxis});
      histos.add(hdcaz[i].data(), pT[i], kTH2F, {ptAxis, dcaZAxis});
      histos.add(hdcaxyphi[i].data(), Form("%s -- 0.9 < #it{p}_{T} < 1.1 GeV/#it{c}", pT[i]), kTH2F, {phiAxis, dcaXyAxis});
    }
  }

  template <PID::ID id, typename T>
  void fillParticleHistos(const T& track)
  {
    const auto& nsigmaTOF = o2::aod::pidutils::tofNSigma<id>(track);
    const auto& nsigmaTPC = o2::aod::pidutils::tpcNSigma<id>(track);
    if (std::sqrt(nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) < 2) {
      histos.fill(HIST(hdcaxy[id]), track.pt(), track.dcaXY());
      if (track.pt() < 1.1 && track.pt() > 0.9) {
        histos.fill(HIST(hdcaxyphi[id]), track.phi(), track.dcaXY());
      }
      histos.fill(HIST(hdcaz[id]), track.pt(), track.dcaZ());
    }
    if (!track.isGlobalTrack()) {
      return;
    }
    histos.fill(HIST(htpcsignal[id]), track.tpcInnerParam(), track.tpcSignal(), nsigmaTOF);
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                    aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                    aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                    aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                                    aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                    aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                                    aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                                    aod::TrackSelection>;
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               TrackCandidates const& tracks)
  {
    histos.fill(HIST("evsel"), 1);
    if (std::abs(collision.posZ()) > cfgCutVertex) {
      return;
    }
    histos.fill(HIST("evsel"), 2);
    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (const auto& track : tracks) {
      histos.fill(HIST("tracksel"), 1);
      if (std::abs(track.eta()) > cfgCutEta) {
        continue;
      }
      histos.fill(HIST("tracksel"), 2);
      if (!track.isGlobalTrack()) {
        continue;
      }
      histos.fill(HIST("tracksel"), 3);
      if (!track.hasTOF()) {
        continue;
      }
      histos.fill(HIST("tracksel"), 4);

      fillParticleHistos<PID::Electron>(track);
      fillParticleHistos<PID::Muon>(track);
      fillParticleHistos<PID::Pion>(track);
      fillParticleHistos<PID::Kaon>(track);
      fillParticleHistos<PID::Proton>(track);
      fillParticleHistos<PID::Deuteron>(track);
      fillParticleHistos<PID::Triton>(track);
      fillParticleHistos<PID::Helium3>(track);
      fillParticleHistos<PID::Alpha>(track);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<tpcSpectra>(cfgc)};
  if (cfgc.options().get<int>("add-tof-histos")) {
    workflow.push_back(adaptAnalysisTask<tpcPidQaSignalwTof>(cfgc));
  }
  return workflow;
}
