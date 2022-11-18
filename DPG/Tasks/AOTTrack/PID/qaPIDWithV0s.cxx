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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;

using PIDTracks = soa::Join<aod::Tracks, aod::TracksExtra,
                            aod::pidTPCFullPi, aod::pidTPCFullPr,
                            aod::pidTOFFullPi, aod::pidTOFFullPr>;
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

struct pidQaWithV0s {
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, ""};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {1000, -100.f, 100.f}, "Binning for the nSigma histograms"};
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> nSigTPC{"nSigTPC", 10., "nSigTPC"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  HistogramRegistry histos{"K0sTrackingEfficiency"};
#define fillHistogram(name, ...) histos.fill(HIST(name), __VA_ARGS__)

  void init(InitContext const&)
  {
    const AxisSpec mAxisK0s{200, 0.4f, 0.6f, "#it{m} (GeV/#it{c}^{2})"};
    const AxisSpec mAxisLambda{400, 1.0f, 1.250f, "#it{m} (GeV/#it{c}^{2})"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nsigmaTOFAxis{nSigmaBins, "N#sigma(TOF)"};
    const AxisSpec nsigmaTPCAxis{nSigmaBins, "N#sigma(TPC)"};

    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1F, {{2, 0.5, 2.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Ev. sel. passed");

    histos.add("K0s/mass", "mass", HistType::kTH1F, {mAxisK0s});
    histos.add("K0s/masscut", "masscut", HistType::kTH1F, {mAxisK0s});
    histos.add("K0s/TOF/pos", "TOF pos", HistType::kTH2F, {ptAxis, nsigmaTOFAxis});
    histos.add("K0s/TOF/neg", "TOF neg", HistType::kTH2F, {ptAxis, nsigmaTOFAxis});
    histos.add("K0s/TPC/pos", "TPC pos", HistType::kTH2F, {ptAxis, nsigmaTPCAxis});
    histos.add("K0s/TPC/neg", "TPC neg", HistType::kTH2F, {ptAxis, nsigmaTPCAxis});

    histos.add("Lambda/mass", "mass", HistType::kTH1F, {mAxisLambda});
    histos.add("Lambda/masscut", "masscut", HistType::kTH1F, {mAxisLambda});
    histos.add<TH2>("Lambda/TOF/pos", "TOF pos", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(p)");
    histos.add<TH2>("Lambda/TOF/neg", "TOF neg", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(#pi^{-})");
    histos.add<TH2>("Lambda/TPC/pos", "TPC pos", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(p)");
    histos.add<TH2>("Lambda/TPC/neg", "TPC neg", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(#pi^{-})");

    histos.add("AntiLambda/mass", "mass", HistType::kTH1F, {mAxisLambda});
    histos.add("AntiLambda/masscut", "mass", HistType::kTH1F, {mAxisLambda});
    histos.add<TH2>("AntiLambda/TOF/pos", "TOF pos", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(#pi^{+})");
    histos.add<TH2>("AntiLambda/TOF/neg", "TOF neg", HistType::kTH2F, {ptAxis, nsigmaTOFAxis})->GetYaxis()->SetTitle("N#sigma_{TOF}(#bar{p})");
    histos.add<TH2>("AntiLambda/TPC/pos", "TPC pos", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(#pi^{+})");
    histos.add<TH2>("AntiLambda/TPC/neg", "TPC neg", HistType::kTH2F, {ptAxis, nsigmaTPCAxis})->GetYaxis()->SetTitle("N#sigma_{TPC}(#bar{p})");
  }

  template <typename T1, typename T2, typename C>
  bool acceptK0s(const T1& v0, const T2& ntrack, const T2& ptrack, const C& collision)
  {
    // Apply selections on V0
    if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa)
      return kFALSE;
    if (TMath::Abs(v0.yK0Short()) > rapidity)
      return kFALSE;

    // Apply selections on V0 daughters
    if (!ntrack.hasTPC() || !ptrack.hasTPC())
      return kFALSE;
    if (ntrack.tpcNSigmaPi() > nSigTPC || ptrack.tpcNSigmaPi() > nSigTPC)
      return kFALSE;
    return kTRUE;
  }

  template <typename T1, typename T2, typename C>
  bool acceptLambda(const T1& v0, const T2& ntrack, const T2& ptrack, const C& collision)
  {
    // Apply selections on V0
    if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa)
      return kFALSE;
    if (TMath::Abs(v0.yLambda()) > rapidity)
      return kFALSE;

    // Apply selections on V0 daughters
    if (!ntrack.hasTPC() || !ptrack.hasTPC())
      return kFALSE;
    if (ntrack.tpcNSigmaPi() > nSigTPC && ptrack.tpcNSigmaPr() > nSigTPC)
      return kFALSE;
    if (ntrack.tpcNSigmaPr() > nSigTPC && ptrack.tpcNSigmaPi() > nSigTPC)
      return kFALSE;
    return kTRUE;
  }

  void process(SelectedCollisions::iterator const& collision,
               aod::V0Datas const& fullV0s,
               PIDTracks const& tracks)
  {
    histos.fill(HIST("evsel"), 0.);
    if (eventSelection && !collision.sel8()) {
      return;
    }
    histos.fill(HIST("evsel"), 1.);

    for (auto& v0 : fullV0s) {
      const auto& posTrack = v0.posTrack_as<PIDTracks>();
      const auto& negTrack = v0.negTrack_as<PIDTracks>();
      if (acceptK0s(v0, negTrack, posTrack, collision)) {
        fillHistogram("K0s/mass", v0.mK0Short());
        if (v0.mK0Short() < 0.47f || v0.mK0Short() > 0.52f) {
          continue;
        }
        fillHistogram("K0s/masscut", v0.mK0Short());
        fillHistogram("K0s/TOF/pos", posTrack.pt(), posTrack.tofNSigmaPi());

        fillHistogram("K0s/TOF/neg", negTrack.pt(), negTrack.tofNSigmaPi());
        fillHistogram("K0s/TPC/pos", posTrack.pt(), posTrack.tpcNSigmaPi());
        fillHistogram("K0s/TPC/neg", negTrack.pt(), negTrack.tpcNSigmaPi());
      }

      if (acceptLambda(v0, negTrack, posTrack, collision)) {
        fillHistogram("Lambda/mass", v0.mLambda());
        if (v0.mLambda() > 1.105f && v0.mLambda() < 1.125f) {
          fillHistogram("Lambda/masscut", v0.mLambda());
          fillHistogram("Lambda/TOF/pos", posTrack.pt(), posTrack.tofNSigmaPr());
          fillHistogram("Lambda/TOF/neg", negTrack.pt(), negTrack.tofNSigmaPi());
          fillHistogram("Lambda/TPC/pos", posTrack.pt(), posTrack.tpcNSigmaPr());
          fillHistogram("Lambda/TPC/neg", negTrack.pt(), negTrack.tpcNSigmaPi());
        }

        fillHistogram("AntiLambda/mass", v0.mAntiLambda());
        if (v0.mAntiLambda() > 1.105f && v0.mAntiLambda() < 1.125f) {
          fillHistogram("AntiLambda/masscut", v0.mAntiLambda());
          fillHistogram("AntiLambda/TOF/pos", posTrack.pt(), posTrack.tofNSigmaPi());
          fillHistogram("AntiLambda/TOF/neg", negTrack.pt(), negTrack.tofNSigmaPr());
          fillHistogram("AntiLambda/TPC/pos", posTrack.pt(), posTrack.tpcNSigmaPi());
          fillHistogram("AntiLambda/TPC/neg", negTrack.pt(), negTrack.tpcNSigmaPr());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pidQaWithV0s>(cfgc)};
}
