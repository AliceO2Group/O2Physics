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
/// \author Dushmanta Sahu  (dushmanta.sahus@cern.ch)
/// \since September 1, 2026
/// \file flattenicityPtSpectra.cxx
/// \brief  Analysis to do flattenicity pt spectra

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>

#include <array>
#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

using MyCollisions = soa::Join<aod::Collisions,
                               aod::Mults,
                               aod::FT0sCorrected,
                               aod::EvSels>;

struct FlattenicityPtSpectra {

  Configurable<float> etaMin{"etaMin", -0.8f, "Min eta for mid-rapidity tracks"};
  Configurable<float> etaMax{"etaMax", 0.8f, "Max eta for mid-rapidity tracks"};
  Configurable<float> ptMin{"ptMin", 0.15f, "Min pT (GeV/c)"};
  Configurable<float> ptMax{"ptMax", 20.0f, "Max pT (GeV/c)"};
  Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "Min TPC found clusters"};
  Configurable<float> maxDCAz{"maxDCAz", 2.0f, "Max |DCAz| (cm)"};
  Configurable<float> maxVtxZ{"maxVtxZ", 10.0f, "Max |vtx z| (cm)"};

  Configurable<int> minActiveFT0Ch{"minActiveFT0Ch", 4, "Min active FT0 channels for flattenicity"};

  static constexpr int NChA = 96;
  static constexpr int NChC = 96;
  static constexpr int NCells = NChA + NChC; // 192 total

  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  template <typename TrackType>
  bool isGoodTrack(TrackType const& track) const
  {
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound.value)
      return false;
    if (track.eta() < etaMin.value || track.eta() > etaMax.value)
      return false;
    if (std::abs(track.dcaZ()) > maxDCAz.value)
      return false;
    if (track.pt() < ptMin.value || track.pt() > ptMax.value)
      return false;
    return true;
  }

  float calculateFlattenicityFT0(aod::FT0 const& ft0, int minActive) const
  {
    std::array<float, NCells> signals{};
    signals.fill(0.f);

    int nActive = 0;

    int chIdx = 0;
    for (const auto& a : ft0.amplitudeA()) {
      if (chIdx < NChA) {
        float amp = static_cast<float>(a);
        if (amp > 0.f) {
          signals[chIdx] = amp;
          ++nActive;
        }
      }
      ++chIdx;
    }

    chIdx = 0;
    for (const auto& a : ft0.amplitudeC()) {
      if (chIdx < NChC) {
        float amp = static_cast<float>(a);
        if (amp > 0.f) {
          signals[NChA + chIdx] = amp;
          ++nActive;
        }
      }
      ++chIdx;
    }

    if (nActive < minActive)
      return -1.f;

    float mRho = 0.f;
    for (int i = 0; i < NCells; ++i) {
      mRho += signals[i];
    }
    mRho /= static_cast<float>(NCells);
    if (mRho <= 0.f)
      return -1.f;

    float sRhoTmp = 0.f;
    for (int i = 0; i < NCells; ++i) {
      if (signals[i] > 0.f) {
        float d = signals[i] - mRho;
        sRhoTmp += d * d;
      }
    }
    sRhoTmp /= static_cast<float>(NCells) * static_cast<float>(NCells);
    float sRho = std::sqrt(sRhoTmp);

    // flattenicity = 1 - sigma/mean, clamped to [0, 1]
    float cv = sRho / mRho;
    return std::max(0.f, std::min(1.f, 1.f - cv));
  }

  void init(InitContext const&)
  {
    AxisSpec ptAxis{200, 0.0f, 20.0f, "#it{p}_{T} (GeV/c)"};
    AxisSpec flatAxis{102, -0.01f, 1.01f, "Flattenicity"};
    AxisSpec multAxis{300, 0.f, 3000.f, "FT0M amplitude (a.u.)"};
    AxisSpec nchAxis{500, -0.5f, 499.5f, "Raw N_{ch} (|#eta|<0.8, p_{T}>0.15)"};
    AxisSpec vtxzAxis{200, -20.f, 20.f, "v_{z} (cm)"};

    registry.add("hEventCounter",
                 "Event counter",
                 HistType::kTH1F, {{5, 0.5, 5.5}});
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(1, "All");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(2, "vtxZ");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(3, "hasFT0");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(4, "FT0signal");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(5, "goodFlat");

    registry.add("hVtxZ",
                 "Primary vertex z",
                 HistType::kTH1F, {vtxzAxis});

    registry.add("hFT0MMultiplicity",
                 "FT0M amplitude sum",
                 HistType::kTH1F, {multAxis});

    registry.add("hFlattenicityDistribution",
                 "Flattenicity (Gyula method, 192 FT0 ch)",
                 HistType::kTH1F, {flatAxis});

    registry.add("hFT0AAmplitude",
                 "Total FT0A amplitude",
                 HistType::kTH1F, {{200, 0.f, 2000.f, "Amplitude (a.u.)"}});

    registry.add("hFT0CAmplitude",
                 "Total FT0C amplitude",
                 HistType::kTH1F, {{200, 0.f, 2000.f, "Amplitude (a.u.)"}});

    // QA: active channel count — use this to tune minActiveFT0Ch
    registry.add("hNActiveFT0Channels",
                 "Active FT0 channels per event (amp > 0)",
                 HistType::kTH1F, {{193, -0.5f, 192.5f, "N active channels"}});

    registry.add("hEtaDistribution",
                 "Eta distribution (selected tracks)",
                 HistType::kTH1F, {{100, -1.0f, 1.0f, "#eta"}});

    registry.add("hPhiDistribution",
                 "Phi distribution (selected tracks)",
                 HistType::kTH1F, {{100, 0.f, o2::constants::math::TwoPI, "#phi"}});

    registry.add("hPtSpectrum_All",
                 "pT spectrum (all selected events)",
                 HistType::kTH1F, {ptAxis});

    registry.add("hFT0MMultVsFlattenicity",
                 "FT0M multiplicity vs Flattenicity;"
                 "FT0M amplitude (a.u.);Flattenicity",
                 HistType::kTH2F, {multAxis, flatAxis});

    registry.add("hFT0MMultVsFlattenicityVsPt",
                 "FT0M multiplicity vs Flattenicity vs #it{p}_{T};"
                 "FT0M amplitude (a.u.);Flattenicity;#it{p}_{T} (GeV/c)",
                 HistType::kTH3F, {multAxis, flatAxis, ptAxis});

    registry.add("hNchDistribution",
                 "Raw charged particle multiplicity",
                 HistType::kTH1F, {nchAxis});

    registry.add("hNchVsFlattenicity",
                 "Raw N_{ch} vs Flattenicity;"
                 "Raw N_{ch};Flattenicity",
                 HistType::kTH2F, {nchAxis, flatAxis});

    registry.add("hNchVsFlattenicityVsPt",
                 "Raw N_{ch} vs Flattenicity vs #it{p}_{T};"
                 "Raw N_{ch};Flattenicity;#it{p}_{T} (GeV/c)",
                 HistType::kTH3F, {nchAxis, flatAxis, ptAxis});

    LOG(info) << "FlattenicityPtSpectra::init — done"
              << "  NCells=" << NCells
              << "  minActiveFT0Ch=" << minActiveFT0Ch.value
              << "  Flattenicity method: Gyula (fixed denominator = total channels)";
  }

  void process(MyCollisions::iterator const& collision,
               MyTracks const& tracks,
               aod::FT0s const& /*ft0s*/)
  {
    registry.fill(HIST("hEventCounter"), 1.f);

    if (std::abs(collision.posZ()) > maxVtxZ.value)
      return;
    registry.fill(HIST("hEventCounter"), 2.f);
    registry.fill(HIST("hVtxZ"), collision.posZ());

    if (!collision.has_foundFT0())
      return;
    registry.fill(HIST("hEventCounter"), 3.f);

    auto ft0 = collision.foundFT0();

    float sumA = 0.f, sumC = 0.f;
    int nActive = 0;
    for (const auto& a : ft0.amplitudeA()) {
      if (a > 0.f) {
        sumA += a;
        ++nActive;
      }
    }
    for (const auto& a : ft0.amplitudeC()) {
      if (a > 0.f) {
        sumC += a;
        ++nActive;
      }
    }

    registry.fill(HIST("hFT0AAmplitude"), sumA);
    registry.fill(HIST("hFT0CAmplitude"), sumC);
    registry.fill(HIST("hNActiveFT0Channels"), static_cast<float>(nActive));

    if (sumA + sumC <= 0.f)
      return;
    registry.fill(HIST("hEventCounter"), 4.f);

    float flattenicity = calculateFlattenicityFT0(ft0, minActiveFT0Ch.value);
    if (flattenicity < 0.f)
      return;
    registry.fill(HIST("hEventCounter"), 5.f);

    float multFT0M = collision.multFT0M();

    registry.fill(HIST("hFT0MMultiplicity"), multFT0M);
    registry.fill(HIST("hFlattenicityDistribution"), flattenicity);
    registry.fill(HIST("hFT0MMultVsFlattenicity"), multFT0M, flattenicity);

    std::vector<typename MyTracks::iterator> goodTracks;
    goodTracks.reserve(1000);

    for (const auto& track : tracks) {
      if (!isGoodTrack(track))
        continue;
      goodTracks.push_back(track);
    }

    int rawNch = goodTracks.size();

    // Fill per-event histograms with Nch
    registry.fill(HIST("hNchDistribution"), rawNch);
    registry.fill(HIST("hNchVsFlattenicity"), rawNch, flattenicity);

    // Fill per-track histograms using the collected good tracks
    for (const auto& track : goodTracks) {
      registry.fill(HIST("hEtaDistribution"), track.eta());
      registry.fill(HIST("hPhiDistribution"), track.phi());
      registry.fill(HIST("hPtSpectrum_All"), track.pt());
      registry.fill(HIST("hFT0MMultVsFlattenicityVsPt"),
                    multFT0M, flattenicity, track.pt());

      // NEW: Fill Nch-based 3D histogram
      registry.fill(HIST("hNchVsFlattenicityVsPt"),
                    rawNch, flattenicity, track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlattenicityPtSpectra>(cfgc)};
}
