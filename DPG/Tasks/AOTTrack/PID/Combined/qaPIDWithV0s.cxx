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
/// \file   qaPIDWithV0s.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2022-11-14
/// \brief  Task to monitor the PID performance making use of V0s
///

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;

using PIDTracks = soa::Join<aod::Tracks, aod::TracksExtra,
                            aod::pidTOFbeta, aod::pidTOFmass, aod::TrackSelection,
                            aod::pidTPCFullPi, aod::pidTPCFullPr,
                            aod::pidTOFFullPi, aod::pidTOFFullPr>;
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

struct pidQaWithV0s {
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, ""};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {1000, -100.f, 100.f}, "Binning for the nSigma histograms"};
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5, "minimumV0Radius (cm)"};
  Configurable<float> nSigTPC{"nSigTPC", 10., "nSigTPC"};
  Configurable<float> minMassK0s{"minMassK0s", 0.47f, "Minimum K0s mass"};
  Configurable<float> maxMassK0s{"maxMassK0s", 0.52f, "Maximum K0s mass"};
  Configurable<float> minMassLambda{"minMassLambda", 1.105f, "Minimum Lambda mass"};
  Configurable<float> maxMassLambda{"maxMassLambda", 1.125f, "Maximum Lambda mass"};
  Configurable<float> minMassAntiLambda{"minMassAntiLambda", 1.105f, "Minimum AntiLambda mass"};
  Configurable<float> maxMassAntiLambda{"maxMassAntiLambda", 1.125f, "Maximum AntiLambda mass"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  HistogramRegistry histos{"K0sTrackingEfficiency"};
#define fillHistogram(name, ...) histos.fill(HIST(name), __VA_ARGS__)

  void init(InitContext const&)
  {
    const AxisSpec cospaAxis{200, 0.9f, 1.01f, "Cos. PA"};
    const AxisSpec decayRadiusAxis{400, 0.f, 400.f, "Decay radius (cm)"};
    const AxisSpec mAxisK0s{200, 0.4f, 0.6f, "#it{m} (GeV/#it{c}^{2})"};
    const AxisSpec mAxisLambda{400, 1.0f, 1.250f, "#it{m} (GeV/#it{c}^{2})"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec betaAxis{100, 0, 2, "TOF #it{#beta}"};
    const AxisSpec dedxAxis{500, 0, 1000, "d#it{E}/d#it{x} A.U."};
    const AxisSpec nsigmaTOFAxis{nSigmaBins, "N#sigma(TOF)"};
    const AxisSpec nsigmaTPCAxis{nSigmaBins, "N#sigma(TPC)"};

    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1F, {{2, 0.5, 2.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Ev. sel. passed");

    h = histos.add<TH1>("tracks", "tracks", HistType::kTH1F, {{6, 0.5, 6.5}});
    h->GetXaxis()->SetBinLabel(1, "tracks read");
    h->GetXaxis()->SetBinLabel(2, "pos tracks read");
    h->GetXaxis()->SetBinLabel(3, "neg tracks read");
    h->GetXaxis()->SetBinLabel(4, "tracks passed cuts");
    h->GetXaxis()->SetBinLabel(5, "pos tracks passed cuts");
    h->GetXaxis()->SetBinLabel(6, "neg tracks passed cuts");

    h = histos.add<TH1>("v0s", "v0s", HistType::kTH1F, {{2, 0.5, 2.5}});
    h->GetXaxis()->SetBinLabel(1, "V0s read");
    h->GetXaxis()->SetBinLabel(2, "V0s accepted");
    histos.add("beforeselections/cospa", "cospa", HistType::kTH1F, {cospaAxis});
    histos.add("afterselections/cospa", "cospa", HistType::kTH1F, {cospaAxis});
    histos.add("beforeselections/decayradius", "decayradius", HistType::kTH1F, {decayRadiusAxis});
    histos.add("afterselections/decayradius", "decayradius", HistType::kTH1F, {decayRadiusAxis});

    histos.add("K0s/mass", "mass", HistType::kTH1F, {mAxisK0s});
    histos.add("K0s/masscut", "masscut", HistType::kTH1F, {mAxisK0s});
    histos.add<TH2>("K0s/TOF/pos", "TOF pos", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(#pi^{+})");
    histos.add<TH2>("K0s/TOF/neg", "TOF neg", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(#pi^{-})");
    histos.add<TH2>("K0s/TPC/pos", "TPC pos", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(#pi^{+})");
    histos.add<TH2>("K0s/TPC/neg", "TPC neg", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(#pi^{-})");
    histos.add<TH2>("K0s/TOF/beta/pos", "TOF beta pos", HistType::kTH2F, {ptAxis, betaAxis});
    histos.add<TH2>("K0s/TOF/beta/neg", "TOF beta neg", HistType::kTH2F, {ptAxis, betaAxis});
    histos.add<TH2>("K0s/TPC/dedx/pos", "TPC dedx pos", HistType::kTH2F, {ptAxis, dedxAxis});
    histos.add<TH2>("K0s/TPC/dedx/neg", "TPC dedx neg", HistType::kTH2F, {ptAxis, dedxAxis});

    histos.add("Lambda/mass", "mass", HistType::kTH1F, {mAxisLambda});
    histos.add("Lambda/masscut", "masscut", HistType::kTH1F, {mAxisLambda});
    histos.add<TH2>("Lambda/TOF/pos", "TOF pos", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(p)");
    histos.add<TH2>("Lambda/TOF/neg", "TOF neg", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(#pi^{-})");
    histos.add<TH2>("Lambda/TPC/pos", "TPC pos", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(p)");
    histos.add<TH2>("Lambda/TPC/neg", "TPC neg", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(#pi^{-})");
    histos.add<TH2>("Lambda/TOF/beta/pos", "TOF beta pos", HistType::kTH2F, {ptAxis, betaAxis});
    histos.add<TH2>("Lambda/TOF/beta/neg", "TOF beta neg", HistType::kTH2F, {ptAxis, betaAxis});
    histos.add<TH2>("Lambda/TPC/dedx/pos", "TPC dedx pos", HistType::kTH2F, {ptAxis, dedxAxis});
    histos.add<TH2>("Lambda/TPC/dedx/neg", "TPC dedx neg", HistType::kTH2F, {ptAxis, dedxAxis});

    histos.add("AntiLambda/mass", "mass", HistType::kTH1F, {mAxisLambda});
    histos.add("AntiLambda/masscut", "mass", HistType::kTH1F, {mAxisLambda});
    histos.add<TH2>("AntiLambda/TOF/pos", "TOF pos", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(#pi^{+})");
    histos.add<TH2>("AntiLambda/TOF/neg", "TOF neg", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(#bar{p})");
    histos.add<TH2>("AntiLambda/TPC/pos", "TPC pos", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(#pi^{+})");
    histos.add<TH2>("AntiLambda/TPC/neg", "TPC neg", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(#bar{p})");
    histos.add<TH2>("AntiLambda/TOF/beta/pos", "TOF beta pos", HistType::kTH2F, {ptAxis, betaAxis});
    histos.add<TH2>("AntiLambda/TOF/beta/neg", "TOF beta neg", HistType::kTH2F, {ptAxis, betaAxis});
    histos.add<TH2>("AntiLambda/TPC/dedx/pos", "TPC dedx pos", HistType::kTH2F, {ptAxis, dedxAxis});
    histos.add<TH2>("AntiLambda/TPC/dedx/neg", "TPC dedx neg", HistType::kTH2F, {ptAxis, dedxAxis});
  }

  template <typename T1, typename C>
  bool acceptV0s(const T1& v0, const C& /*collision*/) // Apply general selections on V0
  {
    float cospa = v0.v0cosPA();
    fillHistogram("beforeselections/cospa", cospa);
    fillHistogram("beforeselections/decayradius", v0.v0radius());
    if (cospa < v0cospa) {
      return false;
    }
    if (v0.v0radius() < minimumV0Radius) {
      return false;
    }
    fillHistogram("afterselections/cospa", cospa);
    fillHistogram("afterselections/decayradius", v0.v0radius());
    return true;
  }

  template <typename T1, typename T2, typename C>
  bool acceptK0s(const T1& v0, const T2& ntrack, const T2& ptrack, const C& /*collision*/) // Apply selection for K0s
  {
    if (TMath::Abs(v0.yK0Short()) > rapidity) {
      return false;
    }

    // Apply selections on V0 daughters
    if (!ntrack.hasTPC() || !ptrack.hasTPC())
      return false;
    if (ntrack.tpcNSigmaPi() > nSigTPC || ptrack.tpcNSigmaPi() > nSigTPC)
      return false;
    return true;
  }

  template <typename T1, typename T2, typename C>
  bool acceptLambda(const T1& v0, const T2& ntrack, const T2& ptrack, const C& /*collision*/) // Apply selection for Lambdas
  {
    if (TMath::Abs(v0.yLambda()) > rapidity) {
      return false;
    }

    // Apply selections on V0 daughters
    if (!ntrack.hasTPC() || !ptrack.hasTPC())
      return false;
    if (ntrack.tpcNSigmaPi() > nSigTPC && ptrack.tpcNSigmaPr() > nSigTPC)
      return false;
    if (ntrack.tpcNSigmaPr() > nSigTPC && ptrack.tpcNSigmaPi() > nSigTPC)
      return false;
    return true;
  }

  void process(SelectedCollisions::iterator const& collision,
               aod::V0Datas const& fullV0s,
               PIDTracks const& tracks)
  {
    fillHistogram("evsel", 1.);
    if (eventSelection && !collision.sel8()) {
      return;
    }
    fillHistogram("evsel", 2.);
    for (auto& trk : tracks) {
      fillHistogram("tracks", 1.f);
      if (trk.sign() > 0) {
        fillHistogram("tracks", 2.f);
      } else {
        fillHistogram("tracks", 3.f);
      }
      if (!trk.isGlobalTrack()) {
        continue;
      }
      fillHistogram("tracks", 4.f);
      if (trk.sign() > 0) {
        fillHistogram("tracks", 5.f);
      } else {
        fillHistogram("tracks", 6.f);
      }
    }

    for (auto& v0 : fullV0s) {
      fillHistogram("v0s", 1.);
      if (!acceptV0s(v0, collision)) {
        continue;
      }
      fillHistogram("v0s", 2.);

      const auto& posTrack = v0.posTrack_as<PIDTracks>();
      const auto& negTrack = v0.negTrack_as<PIDTracks>();
      if (acceptK0s(v0, negTrack, posTrack, collision)) {
        fillHistogram("K0s/mass", v0.mK0Short());
        if (v0.mK0Short() < minMassK0s || v0.mK0Short() > maxMassK0s) {
          continue;
        }
        fillHistogram("K0s/masscut", v0.mK0Short());

        fillHistogram("K0s/TOF/pos", posTrack.pt(), posTrack.tofNSigmaPi());
        fillHistogram("K0s/TOF/neg", negTrack.pt(), negTrack.tofNSigmaPi());
        fillHistogram("K0s/TPC/pos", posTrack.pt(), posTrack.tpcNSigmaPi());
        fillHistogram("K0s/TPC/neg", negTrack.pt(), negTrack.tpcNSigmaPi());

        fillHistogram("K0s/TOF/beta/pos", posTrack.pt(), posTrack.beta());
        fillHistogram("K0s/TOF/beta/neg", negTrack.pt(), negTrack.beta());
        fillHistogram("K0s/TPC/dedx/pos", posTrack.pt(), posTrack.tpcSignal());
        fillHistogram("K0s/TPC/dedx/neg", negTrack.pt(), negTrack.tpcSignal());
      }

      if (acceptLambda(v0, negTrack, posTrack, collision)) {
        fillHistogram("Lambda/mass", v0.mLambda());
        if (v0.mLambda() > minMassLambda && v0.mLambda() < maxMassLambda) {
          fillHistogram("Lambda/masscut", v0.mLambda());

          fillHistogram("Lambda/TOF/pos", posTrack.pt(), posTrack.tofNSigmaPr());
          fillHistogram("Lambda/TOF/neg", negTrack.pt(), negTrack.tofNSigmaPi());
          fillHistogram("Lambda/TPC/pos", posTrack.pt(), posTrack.tpcNSigmaPr());
          fillHistogram("Lambda/TPC/neg", negTrack.pt(), negTrack.tpcNSigmaPi());

          fillHistogram("Lambda/TOF/beta/pos", posTrack.pt(), posTrack.beta());
          fillHistogram("Lambda/TOF/beta/neg", negTrack.pt(), negTrack.beta());
          fillHistogram("Lambda/TPC/dedx/pos", posTrack.pt(), posTrack.tpcSignal());
          fillHistogram("Lambda/TPC/dedx/neg", negTrack.pt(), negTrack.tpcSignal());
        }

        fillHistogram("AntiLambda/mass", v0.mAntiLambda());
        if (v0.mAntiLambda() > minMassAntiLambda && v0.mAntiLambda() < maxMassAntiLambda) {
          fillHistogram("AntiLambda/masscut", v0.mAntiLambda());

          fillHistogram("AntiLambda/TOF/pos", posTrack.pt(), posTrack.tofNSigmaPi());
          fillHistogram("AntiLambda/TOF/neg", negTrack.pt(), negTrack.tofNSigmaPr());
          fillHistogram("AntiLambda/TPC/pos", posTrack.pt(), posTrack.tpcNSigmaPi());
          fillHistogram("AntiLambda/TPC/neg", negTrack.pt(), negTrack.tpcNSigmaPr());

          fillHistogram("AntiLambda/TOF/beta/pos", posTrack.pt(), posTrack.beta());
          fillHistogram("AntiLambda/TOF/beta/neg", negTrack.pt(), negTrack.beta());
          fillHistogram("AntiLambda/TPC/dedx/pos", posTrack.pt(), posTrack.tpcSignal());
          fillHistogram("AntiLambda/TPC/dedx/neg", negTrack.pt(), negTrack.tpcSignal());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pidQaWithV0s>(cfgc)};
}
