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

using PIDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi>;
using PIDTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi>;
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

struct qaK0sTrackingEfficiency {
  ConfigurableAxis nSigmaBins{"nSigmaBins", {1000, -100.f, 100.f}, "Binning for the nSigma histograms"};

  HistogramRegistry registry{"K0sTrackingEfficiency"};
#define fillHistogram(name, ...) registry.fill(HIST(name), __VA_ARGS__)

  void init(InitContext const&)
  {
    const AxisSpec RAxis{500, 0.f, 50.f, "#it{R} (cm)"};
    const AxisSpec pTAxis{200, 0.f, 10.f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec mAxis{200, 0.4f, 0.6f, "#it{m} (GeV/#it{c}^{2})"};
    const AxisSpec etaAxis{200, -1.f, 1.0f, "#it{#eta}"};
    const AxisSpec phiAxis{200, -1.f, 1.0f, "#it{#phi}"};
    const AxisSpec nsigmaAxis{200, -1.f, 1.0f, "#it{#phi}"};
    const AxisSpec statusAxis{2, -0.5f, 1.5f, ""};
    const AxisSpec hitMapAxis{128, -0.5f, 127.5f, ""};

    registry.add("h_EventCounter", "", kTH1D, {{2, -0.5f, 1.5f, ""}});
    registry.get<TH1>(HIST("h_EventCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("h_EventCounter"))->GetXaxis()->SetBinLabel(2, "Selected");

    registry.add("h5_RpTmassITSHitMap", "h5_RpTmassITSHitMap", {HistType::kTHnSparseD, {RAxis, pTAxis, mAxis, hitMapAxis, hitMapAxis}});

    registry.add("Test/h_R", "h_R", {HistType::kTH1D, {RAxis}});
    registry.add("Test/h_pT", "h_pT", {HistType::kTH1D, {pTAxis}});
    registry.add("Test/h_mass", "h_mass", {HistType::kTH1D, {mAxis}});
    registry.add("Test/h_negITSStatus", "h_negITSStatus", {HistType::kTH1I, {statusAxis}});
    registry.add("Test/h_posITSStatus", "h_posITSStatus", {HistType::kTH1I, {statusAxis}});
    registry.add("Test/h_negITSHitMap", "h_negITSHitMap", {HistType::kTH1I, {hitMapAxis}});
    registry.add("Test/h_posITSHitMap", "h_posITSHitMap", {HistType::kTH1I, {hitMapAxis}});
    registry.add("VsPt/h_mass", "h_mass", HistType::kTH2F, {pTAxis, mAxis});
    registry.add("VsPt/h_tofnsigma", "h_tofnsigma", HistType::kTH2F, {pTAxis, nSigmaBins});
    registry.add("VsPt/h_tpcnsigma", "h_tpcnsigma", HistType::kTH2F, {pTAxis, nSigmaBins});
    registry.add("VsEta/h_mass", "h_mass", HistType::kTH2F, {etaAxis, mAxis});
    registry.add("VsPhi/h_mass", "h_mass", HistType::kTH2F, {phiAxis, mAxis});
    registry.add("VsRadius/h_mass", "h_mass", HistType::kTH2F, {RAxis, mAxis});

    if (!doprocessIU) {
      return;
    }
    // IU info
    registry.add("VsRadiusPosIU/h_mass", "h_mass", HistType::kTH2F, {RAxis, mAxis});
    registry.add("VsRadiusPosIU/h_radius", "IU Radius of the negative track", HistType::kTH2F, {RAxis, RAxis});
    registry.add("VsRadiusPosIU/h_phi", "Phi", HistType::kTH2F, {RAxis, phiAxis});

    registry.add("VsRadiusNegIU/h_mass", "h_mass", HistType::kTH2F, {RAxis, mAxis});
    registry.add("VsRadiusNegIU/h_radius", "IU Radius of the positive track", HistType::kTH2F, {RAxis, RAxis});
    registry.add("VsRadiusNegIU/h_phi", "Phi", HistType::kTH2F, {RAxis, phiAxis});
  }

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> nSigTPC{"nSigTPC", 10., "nSigTPC"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  template <typename T1, typename T2, typename C>
  bool acceptV0(const T1& v0, const T2& ntrack, const T2& ptrack, const C& /*collision*/)
  {
    // Apply selections on V0
    if (v0.v0cosPA() < v0cospa)
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

  void process(SelectedCollisions::iterator const& collision, aod::V0Datas const& fullV0s, PIDTracks const&)
  // TODO: add centrality
  {
    registry.fill(HIST("h_EventCounter"), 0.);
    if (eventSelection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("h_EventCounter"), 1.);

    for (auto& v0 : fullV0s) {

      const auto& recopostrack = v0.posTrack_as<PIDTracks>();
      const auto& reconegtrack = v0.negTrack_as<PIDTracks>();

      if (acceptV0(v0, reconegtrack, recopostrack, collision)) {
        registry.fill(HIST("Test/h_R"), v0.v0radius());
        registry.fill(HIST("Test/h_pT"), v0.pt());
        registry.fill(HIST("Test/h_mass"), v0.mK0Short());
        registry.fill(HIST("VsPt/h_mass"), v0.pt(), v0.mK0Short());
        registry.fill(HIST("VsPt/h_tpcnsigma"), v0.pt(), recopostrack.tpcNSigmaPi());
        registry.fill(HIST("VsPt/h_tofnsigma"), v0.pt(), recopostrack.tofNSigmaPi());
        registry.fill(HIST("VsEta/h_mass"), v0.eta(), v0.mK0Short());
        registry.fill(HIST("VsPhi/h_mass"), v0.phi(), v0.mK0Short());
        registry.fill(HIST("VsRadius/h_mass"), v0.v0radius(), v0.mK0Short());

        registry.fill(HIST("Test/h_negITSStatus"), reconegtrack.hasITS());
        registry.fill(HIST("Test/h_posITSStatus"), recopostrack.hasITS());

        uint8_t negITSHitMap = reconegtrack.itsClusterMap();
        uint8_t posITSHitMap = recopostrack.itsClusterMap();
        registry.fill(HIST("Test/h_negITSHitMap"), negITSHitMap);
        registry.fill(HIST("Test/h_posITSHitMap"), posITSHitMap);
        registry.fill(HIST("h5_RpTmassITSHitMap"), v0.v0radius(), v0.pt(), v0.mK0Short(), negITSHitMap, posITSHitMap);
      }
    }
  }

  void processIU(SelectedCollisions::iterator const& collision,
                 aod::V0Datas const& fullV0s,
                 PIDTracksIU const&)
  // TODO: add centrality
  {
    if (eventSelection && !collision.sel8()) {
      return;
    }

    float posRadius = 0.f;
    float negRadius = 0.f;
    for (auto& v0 : fullV0s) {
      const auto& posTrackIU = v0.posTrack_as<PIDTracksIU>();
      const auto& negTrackIU = v0.negTrack_as<PIDTracksIU>();
      if (!acceptV0(v0, negTrackIU, posTrackIU, collision)) {
        continue;
      }
      posRadius = sqrt(posTrackIU.x() * posTrackIU.x() + posTrackIU.y() * posTrackIU.y());
      negRadius = sqrt(negTrackIU.x() * negTrackIU.x() + negTrackIU.y() * negTrackIU.y());

      fillHistogram("VsRadiusPosIU/h_mass", posRadius, v0.mK0Short());
      fillHistogram("VsRadiusPosIU/h_radius", posRadius, negRadius);
      fillHistogram("VsRadiusPosIU/h_phi", posRadius, posTrackIU.phi());

      fillHistogram("VsRadiusNegIU/h_mass", negRadius, v0.mK0Short());
      fillHistogram("VsRadiusNegIU/h_radius", negRadius, posRadius);
      fillHistogram("VsRadiusNegIU/h_phi", negRadius, posTrackIU.phi());
    }
  }
  PROCESS_SWITCH(qaK0sTrackingEfficiency, processIU, "Run also on IU", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaK0sTrackingEfficiency>(cfgc)};
}
