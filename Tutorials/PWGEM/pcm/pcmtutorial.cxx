// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \brief tutorial for pcm analysis.
/// \author daiki.sekihata@cern.ch

#include <TMath.h>
#include <Math/Vector4D.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct PCMTutorial {

  Configurable<int> mincrossedrows{"mincrossedrows", 40, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> minpt{"minpt", 0.05, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"Event/hVertexZ", "z vtx; z vtx (cm);Number of Events", {HistType::kTH1F, {{100, -50.f, +50.f}}}},
      {"V0/hPt", "p_{T} of V0;V0 p_{T} (GeV/c)", {HistType::kTH1F, {{1000, 0.f, 10.f}}}},
      {"V0/hMGamma", "mass #gamma;m_{ee} (GeV/c^{2})", {HistType::kTH1F, {{100, 0.f, 0.1f}}}},
      {"V0/hEtaPhi", "#eta vs. #varphi of V0;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.f, TMath::TwoPi()}, {200, -1.f, +1.f}}}},
      {"V0/hAP", "Armenteros Podolanski;#alpha;q_{T} (GeV/c)", {HistType::kTH2F, {{200, -1, +1}, {250, 0.0, 0.25f}}}},
      {"V0Leg/hdEdx_Pin", "dE/dx in TPC;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"Diphoton/hMgg", "2-photon invariant mass;m_{#gamma#gamma} (GeV/c^{2})", {HistType::kTH1F, {{200, 0.f, 0.8f}}}},
    },
  };

  void init(InitContext const&)
  {
    auto hEventCounter = fRegistry.add<TH1>("hEventCounter", "hEventCounter", kTH1F, {{5, 0.5f, 5.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "FT0AND");
    hEventCounter->GetXaxis()->SetBinLabel(3, "|Z_{vtx}| < 10 cm");
  }

  template <typename TTrack, typename TV0>
  bool checkV0(TV0 const& v0)
  {
    if (!checkAP(v0.alpha(), v0.qtarm(), 0.95, 0.03)) { // select photon conversions
      return false;
    }
    auto pos = v0.template posTrack_as<TTrack>(); // positive daughter
    auto ele = v0.template negTrack_as<TTrack>(); // negative daughter
    if (!checkV0leg(pos) || !checkV0leg(ele)) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool checkV0leg(TTrack const& track)
  {
    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }

    if (track.pt() < minpt || abs(track.eta()) > maxeta) {
      return false;
    }
    if (abs(track.dcaXY()) < dcamin || dcamax < abs(track.dcaXY())) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows || track.tpcChi2NCl() > maxchi2tpc) {
      return false;
    }
    if (abs(track.tpcNSigmaEl()) > maxTPCNsigmaEl) {
      return false;
    }
    return true;
  }

  using MyTracks = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullEl>;
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  Preslice<aod::V0Datas> perCollision = aod::v0data::collisionId;

  void process(MyCollisions const& collisions, aod::BCsWithTimestamps const&, aod::V0Datas const& v0s, MyTracks const& tracks)
  {
    for (auto& collision : collisions) {
      fRegistry.fill(HIST("hEventCounter"), 1);

      if (!collision.sel8()) {
        continue;
      }
      fRegistry.fill(HIST("hEventCounter"), 2);

      if (abs(collision.posZ()) > 10.f) {
        continue;
      }
      fRegistry.fill(HIST("hEventCounter"), 3);
      fRegistry.fill(HIST("Event/hVertexZ"), collision.posZ());

      auto v0s_per_coll = v0s.sliceBy(perCollision, collision.globalIndex());
      for (auto& v0 : v0s_per_coll) {
        fRegistry.fill(HIST("V0/hAP"), v0.alpha(), v0.qtarm());

        if (!checkV0<MyTracks>(v0)) {
          continue;
        }

        fRegistry.fill(HIST("V0/hPt"), v0.pt());
        fRegistry.fill(HIST("V0/hEtaPhi"), v0.phi(), v0.eta());
        fRegistry.fill(HIST("V0/hMGamma"), v0.mGamma());

        auto pos = v0.template posTrack_as<MyTracks>(); // positive daughter
        auto ele = v0.template negTrack_as<MyTracks>(); // negative daughter
        for (auto& leg : {pos, ele}) {
          fRegistry.fill(HIST("V0Leg/hdEdx_Pin"), leg.tpcInnerParam(), leg.tpcSignal());

        } // end of leg loop
      }   // end of v0 loop

      for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(v0s_per_coll, v0s_per_coll))) {
        if (!checkV0<MyTracks>(g1) || !checkV0<MyTracks>(g2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        fRegistry.fill(HIST("Diphoton/hMgg"), v12.M());
      } // end of pairing loop
    }   // end of collision loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PCMTutorial>(cfgc, TaskName{"pcm-tutorial"})};
}
