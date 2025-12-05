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
/// \file   qaTOFMC.cxx
/// \author Nicolò Jacazio
/// \brief  Task to produce QA output of the PID with ALICE3 RICH running on the MC.
///

// O2 includes
#include "ALICE3/DataModel/RICH.h"
#include "Common/DataModel/PIDResponseTOF.h"

#include "Framework/AnalysisTask.h"
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
  Configurable<float> minP{"minP", -2, "Minimum momentum plotted (GeV/c)"};
  Configurable<float> maxP{"maxP", 2, "Maximum momentum plotted (GeV/c)"};
  Configurable<int> nBinsNsigma{"nBinsNsigma", 600, "Number of Nsigma bins"};
  Configurable<float> minNsigma{"minNsigma", -100.f, "Minimum Nsigma plotted"};
  Configurable<float> maxNsigma{"maxNsigma", 100.f, "Maximum Nsigma plotted"};
  Configurable<int> nBinsDelta{"nBinsDelta", 600, "Number of delta bins"};
  Configurable<float> minDelta{"minDelta", -0.4f, "Minimum delta plotted (rad)"};
  Configurable<float> maxDelta{"maxDelta", 0.4f, "Maximum delta plotted (rad)"};
  Configurable<int> logAxis{"logAxis", 1, "Flag to use a log momentum axis"};
  Configurable<float> nSigmaNorm{"nSigmaNorm", 0.7071067811865475f, "Normalization for the combined Nsigma"};

  static constexpr int Np = 5;
  static constexpr std::string_view hbRICHDelta[Np] = {"bRICH/delta/El", "bRICH/delta/Mu", "bRICH/delta/Pi", "bRICH/delta/Ka", "bRICH/delta/Pr"};
  static constexpr std::string_view hbRICHNSigma[Np] = {"bRICH/nsigma/El", "bRICH/nsigma/Mu", "bRICH/nsigma/Pi", "bRICH/nsigma/Ka", "bRICH/nsigma/Pr"};
  static constexpr std::string_view hbRICHNSigmaPrm[Np] = {"bRICH/nsigmaprm/El", "bRICH/nsigmaprm/Mu", "bRICH/nsigmaprm/Pi", "bRICH/nsigmaprm/Ka", "bRICH/nsigmaprm/Pr"};
  static constexpr std::string_view hbRICHNSigmaSec[Np] = {"bRICH/nsigmasec/El", "bRICH/nsigmasec/Mu", "bRICH/nsigmasec/Pi", "bRICH/nsigmasec/Ka", "bRICH/nsigmasec/Pr"};
  static constexpr std::string_view hbRICHNSigmaMC[Np] = {"bRICH/nsigmaMC/El", "bRICH/nsigmaMC/Mu", "bRICH/nsigmaMC/Pi", "bRICH/nsigmaMC/Ka", "bRICH/nsigmaMC/Pr"};
  static constexpr std::string_view hbRICHNSigmaMCSec[Np] = {"bRICH/nsigmaMCsec/El", "bRICH/nsigmaMCsec/Mu", "bRICH/nsigmaMCsec/Pi", "bRICH/nsigmaMCsec/Ka", "bRICH/nsigmaMCsec/Pr"};
  static constexpr std::string_view hbRICHNSigmaMCPrm[Np] = {"bRICH/nsigmaMCprm/El", "bRICH/nsigmaMCprm/Mu", "bRICH/nsigmaMCprm/Pi", "bRICH/nsigmaMCprm/Ka", "bRICH/nsigmaMCprm/Pr"};

  static constexpr std::string_view hfRICHDelta[Np] = {"fRICH/delta/El", "fRICH/delta/Mu", "fRICH/delta/Pi", "fRICH/delta/Ka", "fRICH/delta/Pr"};
  static constexpr std::string_view hfRICHNSigma[Np] = {"fRICH/nsigma/El", "fRICH/nsigma/Mu", "fRICH/nsigma/Pi", "fRICH/nsigma/Ka", "fRICH/nsigma/Pr"};
  static constexpr std::string_view hfRICHNSigmaPrm[Np] = {"fRICH/nsigmaprm/El", "fRICH/nsigmaprm/Mu", "fRICH/nsigmaprm/Pi", "fRICH/nsigmaprm/Ka", "fRICH/nsigmaprm/Pr"};
  static constexpr std::string_view hfRICHNSigmaSec[Np] = {"fRICH/nsigmasec/El", "fRICH/nsigmasec/Mu", "fRICH/nsigmasec/Pi", "fRICH/nsigmasec/Ka", "fRICH/nsigmasec/Pr"};
  static constexpr std::string_view hfRICHNSigmaMC[Np] = {"fRICH/nsigmaMC/El", "fRICH/nsigmaMC/Mu", "fRICH/nsigmaMC/Pi", "fRICH/nsigmaMC/Ka", "fRICH/nsigmaMC/Pr"};
  static constexpr std::string_view hfRICHNSigmaMCSec[Np] = {"fRICH/nsigmaMCsec/El", "fRICH/nsigmaMCsec/Mu", "fRICH/nsigmaMCsec/Pi", "fRICH/nsigmaMCsec/Ka", "fRICH/nsigmaMCsec/Pr"};
  static constexpr std::string_view hfRICHNSigmaMCPrm[Np] = {"fRICH/nsigmaMCprm/El", "fRICH/nsigmaMCprm/Mu", "fRICH/nsigmaMCprm/Pi", "fRICH/nsigmaMCprm/Ka", "fRICH/nsigmaMCprm/Pr"};
  static constexpr std::string_view hfRICHNSigmaVsp[Np] = {"fRICH/nsigmavsp/El", "fRICH/nsigmavsp/Mu", "fRICH/nsigmavsp/Pi", "fRICH/nsigmavsp/Ka", "fRICH/nsigmavsp/Pr"};
  static constexpr std::string_view hfRICHNSigmaMCVsp[Np] = {"fRICH/nsigmaMCvsp/El", "fRICH/nsigmaMCvsp/Mu", "fRICH/nsigmaMCvsp/Pi", "fRICH/nsigmaMCvsp/Ka", "fRICH/nsigmaMCvsp/Pr"};
  static constexpr std::string_view hfRICHNSigmaMCSecVsp[Np] = {"fRICH/nsigmaMCsecvsp/El", "fRICH/nsigmaMCsecvsp/Mu", "fRICH/nsigmaMCsecvsp/Pi", "fRICH/nsigmaMCsecvsp/Ka", "fRICH/nsigmaMCsecvsp/Pr"};
  static constexpr std::string_view hfRICHNSigmaMCPrmVsp[Np] = {"fRICH/nsigmaMCprmvsp/El", "fRICH/nsigmaMCprmvsp/Mu", "fRICH/nsigmaMCprmvsp/Pi", "fRICH/nsigmaMCprmvsp/Ka", "fRICH/nsigmaMCprmvsp/Pr"};

  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p"};
  static constexpr int PDGs[Np] = {11, 13, 211, 321, 2212};
  template <uint8_t i>
  void addParticleHistos()
  {
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    if (logAxis) {
      pAxis.makeLogarithmic();
      ptAxis.makeLogarithmic();
    }
    const AxisSpec nsigmaAxis{nBinsNsigma, minNsigma, maxNsigma, Form("N_{#sigma}^{RICH}(%s)", pT[pid_type])};

    TString tit = Form("%s", pT[i]);
    if (useTOF) {
      tit = Form("TOF Selected %s", pT[i]);
    }
    // NSigma
    histos.add(hbRICHNSigmaMC[i].data(), "True " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hbRICHNSigmaMCPrm[i].data(), "True Primary " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hbRICHNSigmaMCSec[i].data(), "True Secondary " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});

    histos.add(hfRICHNSigmaMC[i].data(), "True " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hfRICHNSigmaMCPrm[i].data(), "True Primary " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hfRICHNSigmaMCSec[i].data(), "True Secondary " + tit, HistType::kTH2F, {ptAxis, nsigmaAxis});

    histos.add(hfRICHNSigmaMCVsp[i].data(), "True " + tit, HistType::kTH2F, {pAxis, nsigmaAxis});
    histos.add(hfRICHNSigmaMCPrmVsp[i].data(), "True Primary " + tit, HistType::kTH2F, {pAxis, nsigmaAxis});
    histos.add(hfRICHNSigmaMCSecVsp[i].data(), "True Secondary " + tit, HistType::kTH2F, {pAxis, nsigmaAxis});
  }

  void init(o2::framework::InitContext&)
  {
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    if (logAxis) {
      pAxis.makeLogarithmic();
      ptAxis.makeLogarithmic();
    }

    const AxisSpec sigAxis{1000, 0, 0.3, "Cherenkov angle (rad)"};
    const AxisSpec lengthAxis{1000, 0, 3000, "Track length (cm)"};
    const AxisSpec sigErrAxis{100, 0, 1, "Cherenkov angle error (rad)"};
    const char* detName = useTOF ? "TOF-RICH" : "RICH";
    const AxisSpec nsigmaAxis{nBinsNsigma, minNsigma, maxNsigma, Form("N_{#sigma}^{%s}(%s)", detName, pT[pid_type])};
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("#Delta(%s) (rad)", pT[pid_type])};
    const AxisSpec etaAxis{100, -4, 4, "#it{#eta}"};

    histos.add("event/vertexz", ";Vtx_{z} (cm);Entries", kTH1F, {{100, -20, 20}});
    histos.add("particle/p", "", kTH1F, {pAxis});
    histos.add("particle/pt", "", kTH1F, {ptAxis});
    histos.add("particle/eta", "", kTH1F, {etaAxis});
    histos.add("tracks/p", "", kTH1F, {pAxis});
    histos.add("tracks/pt", "", kTH1F, {ptAxis});
    histos.add("tracks/eta", "", kTH1F, {etaAxis});
    histos.add("tracks/length", "", kTH1F, {lengthAxis});

#define MakeRICHHistos(rich)                                                                                                 \
  histos.add(rich "/TrackSelection", "", kTH1F, {{10, 0.5, 10.5}});                                                          \
  histos.get<TH1>(HIST(rich "/TrackSelection"))->GetXaxis()->SetBinLabel(1, "Tracks read");                                  \
  histos.get<TH1>(HIST(rich "/TrackSelection"))->GetXaxis()->SetBinLabel(2, "Passed RICH");                                  \
  histos.get<TH1>(HIST(rich "/TrackSelection"))->GetXaxis()->SetBinLabel(3, "Passed TOF");                                   \
  histos.get<TH1>(HIST(rich "/TrackSelection"))->GetXaxis()->SetBinLabel(4, Form("Passed minLength %.2f", minLength.value)); \
  histos.get<TH1>(HIST(rich "/TrackSelection"))->GetXaxis()->SetBinLabel(5, Form("Passed maxLength %.2f", maxLength.value)); \
  histos.get<TH1>(HIST(rich "/TrackSelection"))->GetXaxis()->SetBinLabel(6, Form("Passed minEta %.2f", minEta.value));       \
  histos.get<TH1>(HIST(rich "/TrackSelection"))->GetXaxis()->SetBinLabel(7, Form("Passed maxEta %.2f", maxEta.value));       \
  histos.get<TH1>(HIST(rich "/TrackSelection"))->GetXaxis()->SetBinLabel(8, "Passed PDG");                                   \
  histos.add(rich "/pt", "Unselected", kTH1F, {ptAxis});                                                                     \
  histos.add(rich "/ptPrm", "Primaries", kTH1F, {ptAxis});                                                                   \
  histos.add(rich "/ptSec", "Secondaries", kTH1F, {ptAxis});                                                                 \
  histos.add(rich "/p", "Unselected", kTH1F, {pAxis});                                                                       \
  histos.add(rich "/pPrm", "Primaries", kTH1F, {pAxis});                                                                     \
  histos.add(rich "/pSec", "Secondaries", kTH1F, {pAxis});                                                                   \
  histos.add(rich "/signal", "", kTH1F, {sigAxis});                                                                          \
  histos.add(rich "/eta", "", kTH1F, {etaAxis});                                                                             \
  histos.add(rich "/signalerror", "", kTH1F, {sigErrAxis});                                                                  \
  histos.add(rich "/signalvsP", "Unselected", kTH2F, {pAxis, sigAxis});                                                      \
  histos.add(rich "/signalvsPPrm", "Primaries", kTH2F, {pAxis, sigAxis});                                                    \
  histos.add(rich "/signalvsPSec", "Secondaries", kTH2F, {pAxis, sigAxis});

    MakeRICHHistos("bRICH");
    MakeRICHHistos("fRICH");

#undef MakeRICHHistos

    histos.add(hbRICHDelta[pid_type].data(), "", kTH2F, {pAxis, deltaAxis});
    histos.add(hbRICHNSigma[pid_type].data(), "", HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hbRICHNSigmaPrm[pid_type].data(), "Primary", HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hbRICHNSigmaSec[pid_type].data(), "Secondary", HistType::kTH2F, {ptAxis, nsigmaAxis});

    histos.add(hfRICHDelta[pid_type].data(), "", kTH2F, {pAxis, deltaAxis});
    histos.add(hfRICHNSigma[pid_type].data(), "", HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hfRICHNSigmaPrm[pid_type].data(), "Primary", HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hfRICHNSigmaSec[pid_type].data(), "Secondary", HistType::kTH2F, {ptAxis, nsigmaAxis});
    histos.add(hfRICHNSigmaVsp[pid_type].data(), "", HistType::kTH2F, {pAxis, nsigmaAxis});

    addParticleHistos<0>();
    addParticleHistos<1>();
    addParticleHistos<2>();
    addParticleHistos<3>();
    addParticleHistos<4>();
  }

  template <uint8_t pidIndex, typename T, typename TT>
  void fillNsigma(const T& track, const TT& particle, const float& nsigma)
  {
    if (abs(particle.pdgCode()) != PDGs[pidIndex]) {
      return;
    }
    histos.fill(HIST(hbRICHNSigmaMC[pidIndex]), track.pt(), nsigma);

    if (particle.isPhysicalPrimary()) { // Selecting primaries
      histos.fill(HIST(hbRICHNSigmaMCPrm[pidIndex]), track.pt(), nsigma);
    } else {
      histos.fill(HIST(hbRICHNSigmaMCSec[pidIndex]), track.pt(), nsigma);
    }
  }

  template <uint8_t pidIndex, typename T, typename TT>
  void fillNsigmafRICH(const T& track, const TT& particle, const float& nsigma)
  {
    if (abs(particle.pdgCode()) != PDGs[pidIndex]) {
      return;
    }
    histos.fill(HIST(hfRICHNSigmaMC[pidIndex]), track.pt(), nsigma);
    histos.fill(HIST(hfRICHNSigmaMCVsp[pidIndex]), track.p(), nsigma);

    if (particle.isPhysicalPrimary()) { // Selecting primaries
      histos.fill(HIST(hfRICHNSigmaMCPrm[pidIndex]), track.pt(), nsigma);
      histos.fill(HIST(hfRICHNSigmaMCPrmVsp[pidIndex]), track.p(), nsigma);
    } else {
      histos.fill(HIST(hfRICHNSigmaMCSec[pidIndex]), track.pt(), nsigma);
      histos.fill(HIST(hfRICHNSigmaMCSecVsp[pidIndex]), track.p(), nsigma);
    }
  }

  using Trks = soa::Join<aod::Tracks, aod::RICHTracksIndex, aod::TracksExtra,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                         aod::pidTOFFullKa, aod::pidTOFFullPr, aod::McTrackLabels>;
  using TrksfRICH = soa::Join<aod::Tracks, aod::FRICHTracksIndex, aod::TracksExtra, aod::McTrackLabels>;
  void process(const Trks& tracks,
               const aod::McParticles_000& mcParticles,
               const TrksfRICH& tracksfrich,
               const aod::RICHs&,
               const aod::FRICHs&,
               const aod::Collisions& colls)
  {
    for (const auto& col : colls) {
      histos.fill(HIST("event/vertexz"), col.posZ());
    }
    for (const auto& p : mcParticles) {
      if (pdgCode != 0 && abs(p.pdgCode()) != pdgCode) {
        continue;
      }
      histos.fill(HIST("particle/p"), p.p());
      histos.fill(HIST("particle/pt"), p.pt());
      histos.fill(HIST("particle/eta"), p.eta());
    }

    for (const auto& t : tracks) {
      histos.fill(HIST("tracks/p"), t.p());
      histos.fill(HIST("tracks/pt"), t.pt());
      histos.fill(HIST("tracks/eta"), t.eta());
      histos.fill(HIST("tracks/length"), t.length());
    }

    auto rejectTrack = [&](const auto& t, auto h) {
      if (t.length() < minLength) {
        return true;
      }
      histos.fill(h, 4);
      if (t.length() > maxLength) {
        return true;
      }
      histos.fill(h, 5);
      if (t.eta() < minEta) {
        return true;
      }
      histos.fill(h, 6);
      if (t.eta() > maxEta) {
        return true;
      }
      histos.fill(h, 7);
      const auto mcParticle = t.template mcParticle_as<aod::McParticles_000>();
      if (pdgCode != 0 && abs(mcParticle.pdgCode()) != pdgCode) {
        return true;
      }
      histos.fill(h, 8);
      return false;
    };

    for (const auto& track : tracks) { // Barrel RICH
      histos.fill(HIST("bRICH/TrackSelection"), 1);
      if (!track.has_rich()) {
        continue;
      }
      histos.fill(HIST("bRICH/TrackSelection"), 2);
      if (useTOF && !track.hasTOF()) {
        continue;
      }
      histos.fill(HIST("bRICH/TrackSelection"), 3);
      if (rejectTrack(track, HIST("bRICH/TrackSelection"))) {
        continue;
      }

      const float delta = track.rich().richDelta(pid_type);
      float nsigma = track.rich().richNsigma(pid_type);
      if (useTOF) {
        if constexpr (pid_type == 0) { // Electron
          nsigma += track.tofNSigmaEl();
        } else if constexpr (pid_type == 1) { // Muon
          nsigma += track.tofNSigmaMu();
        } else if constexpr (pid_type == 2) { // Pion
          nsigma += track.tofNSigmaPi();
        } else if constexpr (pid_type == 3) { // Kaon
          nsigma += track.tofNSigmaKa();
        } else if constexpr (pid_type == 4) { // Proton
          nsigma += track.tofNSigmaPr();
        }
        nsigma *= nSigmaNorm; // Normalize to 1
      }

      histos.fill(HIST("bRICH/p"), track.p());
      histos.fill(HIST("bRICH/pt"), track.pt());
      histos.fill(HIST("bRICH/eta"), track.eta());
      histos.fill(HIST("bRICH/signal"), track.rich().richSignal());
      histos.fill(HIST("bRICH/signalerror"), track.rich().richSignalError());
      histos.fill(HIST("bRICH/signalvsP"), track.p(), track.rich().richSignal());

      histos.fill(HIST(hbRICHNSigma[pid_type]), track.pt(), nsigma);
      histos.fill(HIST(hbRICHDelta[pid_type]), track.p(), delta);
      // const auto mcParticle = labels.iteratorAt(track.globalIndex()).mcParticle_as<aod::McParticles_000>();
      const auto mcParticle = track.mcParticle_as<aod::McParticles_000>();
      if (mcParticle.isPhysicalPrimary()) { // Selecting primaries
        histos.fill(HIST(hbRICHNSigmaPrm[pid_type]), track.pt(), nsigma);
        histos.fill(HIST("bRICH/signalvsPPrm"), track.p(), track.rich().richSignal());
        histos.fill(HIST("bRICH/pPrm"), track.p());
        histos.fill(HIST("bRICH/ptPrm"), track.pt());
      } else {
        histos.fill(HIST(hbRICHNSigmaSec[pid_type]), track.pt(), nsigma);
        histos.fill(HIST("bRICH/signalvsPSec"), track.p(), track.rich().richSignal());
        histos.fill(HIST("bRICH/pSec"), track.p());
        histos.fill(HIST("bRICH/ptSec"), track.pt());
      }
      fillNsigma<0>(track, mcParticle, nsigma);
      fillNsigma<1>(track, mcParticle, nsigma);
      fillNsigma<2>(track, mcParticle, nsigma);
      fillNsigma<3>(track, mcParticle, nsigma);
      fillNsigma<4>(track, mcParticle, nsigma);
    }

    for (const auto& track : tracksfrich) { // Forward RICH
      histos.fill(HIST("fRICH/TrackSelection"), 1);
      if (!track.has_frich()) {
        continue;
      }
      histos.fill(HIST("fRICH/TrackSelection"), 2);
      if (rejectTrack(track, HIST("fRICH/TrackSelection"))) {
        continue;
      }

      const float delta = track.frich().frichDelta(pid_type);
      const float nsigma = track.frich().frichNsigma(pid_type);

      histos.fill(HIST("fRICH/p"), track.p());
      histos.fill(HIST("fRICH/pt"), track.pt());
      histos.fill(HIST("fRICH/eta"), track.eta());
      histos.fill(HIST("fRICH/signal"), track.frich().frichSignal());
      histos.fill(HIST("fRICH/signalerror"), track.frich().frichSignalError());
      histos.fill(HIST("fRICH/signalvsP"), track.p(), track.frich().frichSignal());

      histos.fill(HIST(hfRICHNSigma[pid_type]), track.pt(), nsigma);
      histos.fill(HIST(hfRICHDelta[pid_type]), track.p(), delta);
      const auto mcParticle = track.mcParticle_as<aod::McParticles_000>();
      if (mcParticle.isPhysicalPrimary()) { // Selecting primaries
        histos.fill(HIST(hfRICHNSigmaPrm[pid_type]), track.pt(), nsigma);
        histos.fill(HIST("fRICH/signalvsPPrm"), track.p(), track.frich().frichSignal());
        histos.fill(HIST("fRICH/pPrm"), track.p());
        histos.fill(HIST("fRICH/ptPrm"), track.pt());
      } else {
        histos.fill(HIST(hfRICHNSigmaSec[pid_type]), track.pt(), nsigma);
        histos.fill(HIST("fRICH/signalvsPSec"), track.p(), track.frich().frichSignal());
        histos.fill(HIST("fRICH/pSec"), track.p());
        histos.fill(HIST("fRICH/ptSec"), track.pt());
      }
      fillNsigmafRICH<0>(track, mcParticle, nsigma);
      fillNsigmafRICH<1>(track, mcParticle, nsigma);
      fillNsigmafRICH<2>(track, mcParticle, nsigma);
      fillNsigmafRICH<3>(track, mcParticle, nsigma);
      fillNsigmafRICH<4>(track, mcParticle, nsigma);
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
