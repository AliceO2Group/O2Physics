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
/// \file   spectraTOFtiny.h
/// \author Nicolo' Jacazio
///
/// \brief Task for the analysis of the spectra with the TOF detector using the tiny tables
///

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Spectra task
struct tofSpectraTiny {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hp[Np] = {"p/El", "p/Mu", "p/Pi", "p/Ka", "p/Pr", "p/De", "p/Tr", "p/He", "p/Al"};
  static constexpr std::string_view hpt[Np] = {"pt/El", "pt/Mu", "pt/Pi", "pt/Ka", "pt/Pr", "pt/De", "pt/Tr", "pt/He", "pt/Al"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

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
    histos.add("p/Unselected", "Unselected;#it{p} (GeV/#it{c})", kTH1F, {{100, 0, 20}});
    histos.add("pt/Unselected", "Unselected;#it{p}_{T} (GeV/#it{c})", kTH1F, {{100, 0, 20}});
    for (int i = 0; i < Np; i++) {
      histos.add(hp[i].data(), Form("%s;#it{p} (GeV/#it{c})", pT[i]), kTH1F, {{100, 0, 20}});
      histos.add(hpt[i].data(), Form("%s;#it{p}_{T} (GeV/#it{c})", pT[i]), kTH1F, {{100, 0, 20}});
    }
  }

  template <std::size_t i, typename T>
  void fillParticleHistos(const T& track, const float& nsigma)
  {
    if (abs(nsigma) > cfgNSigmaCut) {
      return;
    }
    histos.fill(HIST(hp[i]), track.p());
    histos.fill(HIST(hpt[i]), track.pt());
  }

  //Defining filters and input
  Configurable<float> cfgNSigmaCut{"cfgNSigmaCut", 3, "Value of the Nsigma cut"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra,
                                    aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                    aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe,
                                    aod::pidTOFTr, aod::pidTOFHe, aod::pidTOFAl,
                                    aod::pidEvTimeFlags, aod::TOFSignal, aod::TrackSelection>;
  void process(aod::Collision const& collision,
               TrackCandidates const& tracks)
  {
    histos.fill(HIST("evsel"), 1);
    if (abs(collision.posZ()) > cfgCutVertex) {
      return;
    }
    histos.fill(HIST("evsel"), 2);
    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (const auto& track : tracks) {
      histos.fill(HIST("tracksel"), 1);
      if (abs(track.eta()) > cfgCutEta) {
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

      histos.fill(HIST("p/Unselected"), track.p());
      histos.fill(HIST("pt/Unselected"), track.pt());

      fillParticleHistos<0>(track, track.tofNSigmaEl());
      fillParticleHistos<1>(track, track.tofNSigmaMu());
      fillParticleHistos<2>(track, track.tofNSigmaPi());
      fillParticleHistos<3>(track, track.tofNSigmaKa());
      fillParticleHistos<4>(track, track.tofNSigmaPr());
      fillParticleHistos<5>(track, track.tofNSigmaDe());
      fillParticleHistos<6>(track, track.tofNSigmaTr());
      fillParticleHistos<7>(track, track.tofNSigmaHe());
      fillParticleHistos<8>(track, track.tofNSigmaAl());
    }
  } // end of the process function
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofSpectraTiny>(cfgc)};
}
