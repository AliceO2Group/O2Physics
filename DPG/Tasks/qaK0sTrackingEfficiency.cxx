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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;

using PIDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi>;
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

struct qaK0sTrackingEfficiency {

  HistogramRegistry registry{"K0sTrackingEfficiency"};

  void init(InitContext const&)
  {
    const AxisSpec RAxis{100, 0.f, 10.f, "#it{R} (cm)"};
    const AxisSpec pTAxis{200, 0.f, 10.f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec mAxis{200, 0.4f, 0.6f, "#it{m} (GeV/#it{c}^{2})"};
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
  }

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> nSigTPC{"nSigTPC", 10., "nSigTPC"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  template <typename T1, typename T2, typename C>
  bool acceptV0(T1 v0, T2 ntrack, T2 ptrack, C collision)
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

  void process(SelectedCollisions::iterator const& collision, aod::V0Datas const& fullV0s, PIDTracks const& tracks)
  // TODO: add centrality
  {
    registry.fill(HIST("h_EventCounter"), 0.);
    if (eventSelection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("h_EventCounter"), 1.);

    for (auto& v0 : fullV0s) {

      auto reconegtrack = v0.negTrack_as<PIDTracks>();
      auto recopostrack = v0.posTrack_as<PIDTracks>();

      if (acceptV0(v0, reconegtrack, recopostrack, collision)) {
        registry.fill(HIST("Test/h_R"), v0.v0radius());
        registry.fill(HIST("Test/h_pT"), v0.pt());
        registry.fill(HIST("Test/h_mass"), v0.mK0Short());

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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaK0sTrackingEfficiency>(cfgc)};
}
