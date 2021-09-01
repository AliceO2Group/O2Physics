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

// O2 includes
#include "Framework/AnalysisTask.h"
#include "ALICE3/DataModel/RICH.h"
#include "Common/Core/MC.h"
#include "Common/Core/PID/PIDResponse.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{
    {"qa-el", VariantType::Int, 1, {"Produce PID information for the electron mass hypothesis"}},
    {"qa-mu", VariantType::Int, 1, {"Produce PID information for the muon mass hypothesis"}},
    {"qa-pikapr", VariantType::Int, 1, {"Produce PID information for the Pion, Kaon, Proton mass hypothesis"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

namespace o2::aod
{

namespace indices
{
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_INDEX_COLUMN(RICH, rich);
DECLARE_SOA_INDEX_COLUMN(FRICH, frich);
} // namespace indices

DECLARE_SOA_INDEX_TABLE_USER(RICHTracksIndex, Tracks, "RICHTRK", indices::TrackId, indices::RICHId);
DECLARE_SOA_INDEX_TABLE_USER(FRICHTracksIndex, Tracks, "FRICHTRK", indices::TrackId, indices::FRICHId);
} // namespace o2::aod

struct richIndexBuilder { // Builder of the RICH-track index linkage
  Builds<o2::aod::RICHTracksIndex> indB;
  Builds<o2::aod::FRICHTracksIndex> indF;
  void init(o2::framework::InitContext&)
  {
  }
};

template <o2::track::PID::ID pid_type>
struct richPidQaMc {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};
  Configurable<int> pdgCode{"pdgCode", 0, "pdg code of the particles to accept"};
  Configurable<int> useOnlyPhysicsPrimary{"useOnlyPhysicsPrimary", 1,
                                          "Whether to use only physical primary particles."};
  Configurable<int> useTOF{"useTOF", 0,
                           "Whether to use the TOF information"};
  Configurable<float> minLength{"minLength", 0, "Minimum length of accepted tracks (cm)"};
  Configurable<float> maxLength{"maxLength", 1000, "Maximum length of accepted tracks (cm)"};
  Configurable<float> minEta{"minEta", -1.4, "Minimum eta of accepted tracks"};
  Configurable<float> maxEta{"maxEta", 1.4, "Maximum eta of accepted tracks"};
  Configurable<int> nBinsP{"nBinsP", 500, "Number of momentum bins"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum plotted (GeV/c)"};
  Configurable<float> maxP{"maxP", 100, "Maximum momentum plotted (GeV/c)"};
  Configurable<int> nBinsNsigma{"nBinsNsigma", 600, "Number of Nsigma bins"};
  Configurable<float> minNsigma{"minNsigma", -100.f, "Minimum Nsigma plotted"};
  Configurable<float> maxNsigma{"maxNsigma", 100.f, "Maximum Nsigma plotted"};
  Configurable<int> nBinsDelta{"nBinsDelta", 600, "Number of delta bins"};
  Configurable<float> minDelta{"minDelta", -0.4f, "Minimum delta plotted (rad)"};
  Configurable<float> maxDelta{"maxDelta", 0.4f, "Maximum delta plotted (rad)"};
  Configurable<int> logAxis{"logAxis", 1, "Flag to use a log momentum axis"};

  static constexpr int Np = 5;
  static constexpr std::string_view hdelta[Np] = {"delta/El", "delta/Mu", "delta/Pi", "delta/Ka", "delta/Pr"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi", "nsigma/Ka", "nsigma/Pr"};
  static constexpr std::string_view hnsigmaprm[Np] = {"nsigmaprm/El", "nsigmaprm/Mu", "nsigmaprm/Pi", "nsigmaprm/Ka", "nsigmaprm/Pr"};
  static constexpr std::string_view hnsigmasec[Np] = {"nsigmasec/El", "nsigmasec/Mu", "nsigmasec/Pi", "nsigmasec/Ka", "nsigmasec/Pr"};
  static constexpr std::string_view hnsigmaMC[Np] = {"nsigmaMC/El", "nsigmaMC/Mu", "nsigmaMC/Pi", "nsigmaMC/Ka", "nsigmaMC/Pr"};
  static constexpr std::string_view hnsigmaMCsec[Np] = {"nsigmaMCsec/El", "nsigmaMCsec/Mu", "nsigmaMCsec/Pi", "nsigmaMCsec/Ka", "nsigmaMCsec/Pr"};
  static constexpr std::string_view hnsigmaMCprm[Np] = {"nsigmaMCprm/El", "nsigmaMCprm/Mu", "nsigmaMCprm/Pi", "nsigmaMCprm/Ka", "nsigmaMCprm/Pr"};

  static constexpr std::string_view hfrichnsigma[Np] = {"fRICH/nsigma/El", "fRICH/nsigma/Mu", "fRICH/nsigma/Pi", "fRICH/nsigma/Ka", "fRICH/nsigma/Pr"};
  static constexpr std::string_view hfrichnsigmaprm[Np] = {"fRICH/nsigmaprm/El", "fRICH/nsigmaprm/Mu", "fRICH/nsigmaprm/Pi", "fRICH/nsigmaprm/Ka", "fRICH/nsigmaprm/Pr"};
  static constexpr std::string_view hfrichnsigmasec[Np] = {"fRICH/nsigmasec/El", "fRICH/nsigmasec/Mu", "fRICH/nsigmasec/Pi", "fRICH/nsigmasec/Ka", "fRICH/nsigmasec/Pr"};

  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p"};
  static constexpr int PDGs[Np] = {11, 13, 211, 321, 2212};
  template <uint8_t i>
  void addParticleHistos()
  {
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogaritmic();
    }
    const AxisSpec nsigmaAxis{nBinsNsigma, minNsigma, maxNsigma, Form("N_{#sigma}^{RICH}(%s)", pT[pid_type])};

    TString tit = Form("%s", pT[i]);
    if (useTOF) {
      tit = Form("TOF Selected %s", pT[i]);
    }
    // NSigma
    histos.add(hnsigmaMC[i].data(), "True " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hnsigmaMCprm[i].data(), "True Primary " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hnsigmaMCsec[i].data(), "True Secondary " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});
  }
  void init(o2::framework::InitContext&)
  {
    AxisSpec momAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    if (logAxis) {
      momAxis.makeLogaritmic();
      ptAxis.makeLogaritmic();
    }

    const AxisSpec sigAxis{1000, 0, 0.3, "Cherenkov angle (rad)"};
    const AxisSpec nsigmaAxis{nBinsNsigma, minNsigma, maxNsigma, Form("N_{#sigma}^{RICH}(%s)", pT[pid_type])};
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("#Delta(%s) (rad)", pT[pid_type])};

    histos.add("event/vertexz", ";Vtx_{z} (cm);Entries", kTH1F, {{100, -20, 20}});
    histos.add("p/Unselected", "Unselected", kTH1F, {momAxis});
    histos.add("p/Prim", "Primaries", kTH1F, {momAxis});
    histos.add("p/Sec", "Secondaries", kTH1F, {momAxis});
    histos.add("pt/Unselected", "Unselected", kTH1F, {momAxis});
    histos.add("qa/signal", "", kTH1F, {sigAxis});
    histos.add("qa/eta", ";#it{#eta}", kTH1F, {{100, -4, 4}});
    histos.add("qa/signalerror", ";Cherenkov angle error (rad)", kTH1F, {{100, 0, 1}});
    histos.add("qa/signalvsP", "", kTH2F, {momAxis, sigAxis});
    histos.add("qa/signalvsPPrim", "", kTH2F, {momAxis, sigAxis});
    histos.add("qa/signalvsPSec", "", kTH2F, {momAxis, sigAxis});
    histos.add(hdelta[pid_type].data(), "", kTH2F, {momAxis, deltaAxis});
    histos.add(hnsigma[pid_type].data(), "", HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hnsigmaprm[pid_type].data(), "Primary", HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hnsigmasec[pid_type].data(), "Secondary", HistType::kTH2F, {ptAxis, nsigmaAxis});

    histos.add(hfrichnsigma[pid_type].data(), "", HistType::kTH2F, {ptAxis, nsigmaAxis});

    addParticleHistos<0>();
    addParticleHistos<1>();
    addParticleHistos<2>();
    addParticleHistos<3>();
    addParticleHistos<4>();
  }

  template <uint8_t pidIndex, typename T, typename TTT, typename TT>
  void fillNsigma(const T& track, const TTT& particle, const TT& mcParticles, const float& nsigma)
  {
    if (abs(particle.pdgCode()) == PDGs[pidIndex]) {
      histos.fill(HIST(hnsigmaMC[pidIndex]), track.pt(), nsigma);

      if (MC::isPhysicalPrimary(particle)) { // Selecting primaries
        histos.fill(HIST(hnsigmaMCprm[pidIndex]), track.pt(), nsigma);
      } else {
        histos.fill(HIST(hnsigmaMCsec[pidIndex]), track.pt(), nsigma);
      }
    }
  }

  using Trks = soa::Join<aod::Tracks, aod::RICHTracksIndex, aod::TracksExtra,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                         aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using TrksfRICH = soa::Join<aod::Tracks, aod::FRICHTracksIndex, aod::TracksExtra>;
  void process(const aod::McParticles& mcParticles,
               const Trks& tracks,
               const TrksfRICH& tracksfrich,
               const aod::McTrackLabels& labels,
               const aod::RICHs&,
               const aod::FRICHs&,
               const aod::Collisions& colls)
  {
    for (const auto& col : colls) {
      histos.fill(HIST("event/vertexz"), col.posZ());
    }

    for (const auto& track : tracks) { // Barrel RICH
      if (!track.has_rich()) {
        continue;
      }
      if (track.length() < minLength) {
        continue;
      }
      if (track.length() > maxLength) {
        continue;
      }
      if (track.eta() > maxEta || track.eta() < minEta) {
        continue;
      }
      const auto mcParticle = labels.iteratorAt(track.globalIndex()).mcParticle();
      if (pdgCode != 0 && abs(mcParticle.pdgCode()) != pdgCode) {
        continue;
      }
      if (useTOF && !track.hasTOF()) {
        continue;
      }

      histos.fill(HIST("p/Unselected"), track.p());
      histos.fill(HIST("pt/Unselected"), track.pt());
      histos.fill(HIST("qa/eta"), track.eta());
      histos.fill(HIST("qa/signal"), track.rich().richSignal());
      histos.fill(HIST("qa/signalerror"), track.rich().richSignalError());
      histos.fill(HIST("qa/signalvsP"), track.p(), track.rich().richSignal());

      float delta = -999.f;
      float nsigma = -999.f;
      if constexpr (pid_type == 0) {
        delta = track.rich().richDeltaEl();
        nsigma = track.rich().richNsigmaEl();
        if (useTOF && abs(track.tofNSigmaEl()) > 3.f) {
          continue;
        }
      } else if constexpr (pid_type == 1) {
        delta = track.rich().richDeltaMu();
        nsigma = track.rich().richNsigmaMu();
        if (useTOF && abs(track.tofNSigmaMu()) > 3.f) {
          continue;
        }
      } else if constexpr (pid_type == 2) {
        delta = track.rich().richDeltaPi();
        nsigma = track.rich().richNsigmaPi();
        if (useTOF && abs(track.tofNSigmaPi()) > 3.f) {
          continue;
        }
      } else if constexpr (pid_type == 3) {
        delta = track.rich().richDeltaKa();
        nsigma = track.rich().richNsigmaKa();
        if (useTOF && abs(track.tofNSigmaKa()) > 3.f) {
          continue;
        }
      } else if constexpr (pid_type == 4) {
        delta = track.rich().richDeltaPr();
        nsigma = track.rich().richNsigmaPr();
        if (useTOF && abs(track.tofNSigmaPr()) > 3.f) {
          continue;
        }
      }
      histos.fill(HIST(hnsigma[pid_type]), track.pt(), nsigma);
      histos.fill(HIST(hdelta[pid_type]), track.p(), delta);
      if (MC::isPhysicalPrimary(mcParticle)) { // Selecting primaries
        histos.fill(HIST(hnsigmaprm[pid_type]), track.pt(), nsigma);
        histos.fill(HIST("p/Prim"), track.p());
      } else {
        histos.fill(HIST(hnsigmasec[pid_type]), track.pt(), nsigma);
        histos.fill(HIST("p/Sec"), track.p());
      }

      fillNsigma<0>(track, mcParticle, mcParticles, nsigma);
      fillNsigma<1>(track, mcParticle, mcParticles, nsigma);
      fillNsigma<2>(track, mcParticle, mcParticles, nsigma);
      fillNsigma<3>(track, mcParticle, mcParticles, nsigma);
      fillNsigma<4>(track, mcParticle, mcParticles, nsigma);
    }

    for (const auto& track : tracksfrich) { // Forward RICH
      if (!track.has_frich()) {
        continue;
      }
      if (track.length() < minLength) {
        continue;
      }
      if (track.length() > maxLength) {
        continue;
      }
      if (track.eta() > maxEta || track.eta() < minEta) {
        continue;
      }
      const auto mcParticle = labels.iteratorAt(track.globalIndex()).mcParticle();
      if (pdgCode != 0 && abs(mcParticle.pdgCode()) != pdgCode) {
        continue;
      }

      float delta = -999.f;
      float nsigma = -999.f;
      if constexpr (pid_type == 0) {
        delta = track.frich().frichDeltaEl();
        nsigma = track.frich().frichNsigmaEl();
      } else if constexpr (pid_type == 1) {
        delta = track.frich().frichDeltaMu();
        nsigma = track.frich().frichNsigmaMu();
      } else if constexpr (pid_type == 2) {
        delta = track.frich().frichDeltaPi();
        nsigma = track.frich().frichNsigmaPi();
      } else if constexpr (pid_type == 3) {
        delta = track.frich().frichDeltaKa();
        nsigma = track.frich().frichNsigmaKa();
      } else if constexpr (pid_type == 4) {
        delta = track.frich().frichDeltaPr();
        nsigma = track.frich().frichNsigmaPr();
      }
      histos.fill(HIST(hfrichnsigma[pid_type]), track.pt(), nsigma);
      if (MC::isPhysicalPrimary(mcParticle)) { // Selecting primaries
        histos.fill(HIST(hfrichnsigmaprm[pid_type]), track.pt(), nsigma);
      } else {
        histos.fill(HIST(hfrichnsigmasec[pid_type]), track.pt(), nsigma);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<richIndexBuilder>(cfg)};
  if (cfg.options().get<int>("qa-el")) {
    workflow.push_back(adaptAnalysisTask<richPidQaMc<PID::Electron>>(cfg, TaskName{"pidRICH-qa-El"}));
  }
  if (cfg.options().get<int>("qa-mu")) {
    workflow.push_back(adaptAnalysisTask<richPidQaMc<PID::Muon>>(cfg, TaskName{"pidRICH-qa-Mu"}));
  }
  if (cfg.options().get<int>("qa-pikapr")) {
    workflow.push_back(adaptAnalysisTask<richPidQaMc<PID::Pion>>(cfg, TaskName{"pidRICH-qa-Pi"}));
    workflow.push_back(adaptAnalysisTask<richPidQaMc<PID::Kaon>>(cfg, TaskName{"pidRICH-qa-Ka"}));
    workflow.push_back(adaptAnalysisTask<richPidQaMc<PID::Proton>>(cfg, TaskName{"pidRICH-qa-Pr"}));
  }
  return workflow;
}
