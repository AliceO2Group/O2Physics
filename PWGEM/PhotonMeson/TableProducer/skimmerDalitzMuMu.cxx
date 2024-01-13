// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \brief write relevant information for dalitz mumu analysis to an AO2D.root file. This file is then the only necessary input to perform pcm analysis.
/// \author daiki.sekihata@cern.ch

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyCollisions = soa::Join<aod::EMReducedEvents, aod::EMReducedEventsMult, aod::EMReducedEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMReducedEventIds>;
using MyTrack = MyTracks::iterator;

struct skimmerDalitzMuMu {
  enum class EM_MuMuPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSmm = -1,
  };

  SliceCache cache;
  Preslice<MyTracks> perCol = o2::aod::emprimarymuon::emreducedeventId;
  Produces<aod::DalitzMuMus> dalitzmumus;
  Produces<o2::aod::DalitzMuMuEMReducedEventIds> dalitz_mumu_eventid;
  Produces<o2::aod::EMReducedEventsNmumu> event_nmumu;

  // Configurables
  Configurable<float> maxMmumu{"maxMmumu", 1.1, "max. mmumu to store mumu pairs"};
  Configurable<bool> storeLS{"storeLS", false, "flag to store LS pairs"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hNpairs", "hNpairs;pair type;Number of Pairs", {HistType::kTH1F, {{3, -1.5f, +1.5f}}}},
    },
  };

  void init(InitContext const&) {}

  template <EM_MuMuPairType pairtype, typename TCollision, typename TTracks1, typename TTracks2>
  int fillPairTable(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    int npair = 0;
    if constexpr (pairtype == EM_MuMuPairType::kULS) { // ULS
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMmumu) { // don't store
          continue;
        }
        float phiv = 0.f;
        float opangle = 0.f;
        float dcaxy1 = t1.dcaXY() / sqrt(t1.cYY());
        float dcaxy2 = t2.dcaXY() / sqrt(t2.cYY());
        float dcamumuxy = sqrt((pow(dcaxy1, 2) + pow(dcaxy2, 2)) / 2.);
        float dcaz1 = t1.dcaZ() / sqrt(t1.cZZ());
        float dcaz2 = t2.dcaZ() / sqrt(t2.cZZ());
        float dcamumuz = sqrt((pow(dcaz1, 2) + pow(dcaz2, 2)) / 2.);

        dalitzmumus(collision.collisionId(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), phiv, opangle, dcamumuxy, dcamumuz, static_cast<int>(pairtype));
        dalitz_mumu_eventid(collision.globalIndex());
        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        npair++;
      }      // end of pairing loop
    } else { // LS
      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMmumu) { // don't store
          continue;
        }
        float phiv = 0.f;
        float opangle = 0.f;
        float dcaxy1 = t1.dcaXY() / sqrt(t1.cYY());
        float dcaxy2 = t2.dcaXY() / sqrt(t2.cYY());
        float dcamumuxy = sqrt((pow(dcaxy1, 2) + pow(dcaxy2, 2)) / 2.);
        float dcaz1 = t1.dcaZ() / sqrt(t1.cZZ());
        float dcaz2 = t2.dcaZ() / sqrt(t2.cZZ());
        float dcamumuz = sqrt((pow(dcaz1, 2) + pow(dcaz2, 2)) / 2.);

        dalitzmumus(collision.collisionId(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), phiv, opangle, dcamumuxy, dcamumuz, static_cast<int>(pairtype));
        dalitz_mumu_eventid(collision.globalIndex());
        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        npair++;
      } // end of pairing loop
    }
    return npair;
  }

  Partition<MyTracks> posTracks = o2::aod::emprimarymuon::sign > 0;
  Partition<MyTracks> negTracks = o2::aod::emprimarymuon::sign < 0;

  void process(MyCollisions const& collisions, MyTracks const& tracks)
  {
    for (auto& collision : collisions) {
      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimarymuon::emreducedeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimarymuon::emreducedeventId, collision.globalIndex(), cache);

      int npair_uls = 0, npair_lspp = 0, npair_lsmm = 0;
      npair_uls = fillPairTable<EM_MuMuPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        npair_lspp = fillPairTable<EM_MuMuPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        npair_lsmm = fillPairTable<EM_MuMuPairType::kLSmm>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
      event_nmumu(npair_uls, npair_lspp, npair_lsmm);
    } // end of collision loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerDalitzMuMu>(cfgc, TaskName{"skimmer-dalitz-mumu"})};
}
