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

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds>;
using MyTrack = MyTracks::iterator;

struct skimmerDalitzMuMu {
  enum class EM_MuMuPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSmm = -1,
  };

  SliceCache cache;
  Preslice<MyTracks> perCol = o2::aod::emprimarymuon::emeventId;
  Produces<aod::DalitzMuMus> dalitzmumus;
  Produces<o2::aod::DalitzMuMuEMEventIds> dalitz_mumu_eventid;
  Produces<o2::aod::EMEventsNmumu> event_nmumu;

  // Configurables
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<float> maxMmumu{"maxMmumu", 1.1, "max. mmumu to store mumu pairs"};
  Configurable<bool> storeLS{"storeLS", false, "flag to store LS pairs"};
  Configurable<float> minpt{"minpt", 0.05, "min pt for track"};
  Configurable<float> maxpt{"maxpt", 0.6, "max pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 10, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<int> minitsncls{"minitsncls", 5, "min. number of ITS clusters"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 5.0, "max. chi2/NclsITS"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0f, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0f, "max DCAz in cm"};
  Configurable<float> dca_3d_sigma_max{"dca_3d_sigma_max", 1.0f, "max DCA 3D in sigma"};
  Configurable<float> minTPCNsigmaMu{"minTPCNsigmaMu", -3.0, "min. TPC n sigma for muon inclusion"};
  Configurable<float> maxTPCNsigmaMu{"maxTPCNsigmaMu", 3.0, "max. TPC n sigma for muon inclusion"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hNpairs", "hNpairs;pair type;Number of Pairs", {HistType::kTH1F, {{3, -1.5f, +1.5f}}}},
    },
  };

  void init(InitContext const&) {}

  std::pair<int8_t, std::set<uint8_t>> itsRequirement = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.

  template <typename TTrack>
  bool checkTrack(TTrack const& track)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsNCls() < minitsncls) {
      return false;
    }

    auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement.first) {
      return false;
    }

    if (track.tpcNClsFound() < min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < min_tpc_cr_findable_ratio) {
      return false;
    }

    // float rel_diff = (track.tpcInnerParam() - track.p())/track.p();
    // if (rel_diff < -0.0156432+-0.00154524/pow(track.p(),2.29244) - 0.02 || -0.0156432+-0.00154524/pow(track.p(),2.29244) + 0.02 < rel_diff) {
    //   return false;
    // }

    float dca_3d = 999.f;
    float det = track.cYY() * track.cZZ() - track.cZY() * track.cZY();
    if (det < 0) {
      dca_3d = 999.f;
    } else {
      float chi2 = (track.dcaXY() * track.dcaXY() * track.cZZ() + track.dcaZ() * track.dcaZ() * track.cYY() - 2. * track.dcaXY() * track.dcaZ() * track.cZY()) / det;
      dca_3d = std::sqrt(std::abs(chi2) / 2.);
    }
    if (dca_3d > dca_3d_sigma_max) {
      return false;
    }

    return true;
  }

  template <EM_MuMuPairType pairtype, typename TCollision, typename TTracks1, typename TTracks2>
  int fillPairTable(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    int npair = 0;
    const float phiv = 0.f;
    const float opangle = 0.f;
    if constexpr (pairtype == EM_MuMuPairType::kULS) { // ULS
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack(t1) || !checkTrack(t2)) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMmumu) { // don't store
          continue;
        }
        dalitzmumus(collision.collisionId(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), v12.Rapidity(), phiv, opangle, static_cast<int>(pairtype));
        dalitz_mumu_eventid(collision.globalIndex());
        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        npair++;
      }      // end of pairing loop
    } else { // LS
      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack(t1) || !checkTrack(t2)) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMmumu) { // don't store
          continue;
        }
        dalitzmumus(collision.collisionId(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), v12.Rapidity(), phiv, opangle, static_cast<int>(pairtype));
        dalitz_mumu_eventid(collision.globalIndex());
        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        npair++;
      } // end of pairing loop
    }
    return npair;
  }

  Partition<MyTracks> posTracks = minpt < o2::aod::track::pt && o2::aod::track::pt < maxpt && nabs(o2::aod::track::eta) < maxeta && o2::aod::emprimarymuon::sign > int8_t(0) && minTPCNsigmaMu < o2::aod::pidtpc::tpcNSigmaMu&& o2::aod::pidtpc::tpcNSigmaMu < maxTPCNsigmaMu;
  Partition<MyTracks> negTracks = minpt < o2::aod::track::pt && o2::aod::track::pt < maxpt && nabs(o2::aod::track::eta) < maxeta && o2::aod::emprimarymuon::sign < int8_t(0) && minTPCNsigmaMu < o2::aod::pidtpc::tpcNSigmaMu && o2::aod::pidtpc::tpcNSigmaMu < maxTPCNsigmaMu;
  void processAnalysis(MyCollisions const& collisions, MyTracks const&)
  {
    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        event_nmumu(0, 0, 0);
        continue;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);

      int npair_uls = 0, npair_lspp = 0, npair_lsmm = 0;
      npair_uls = fillPairTable<EM_MuMuPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        npair_lspp = fillPairTable<EM_MuMuPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        npair_lsmm = fillPairTable<EM_MuMuPairType::kLSmm>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
      event_nmumu(npair_uls, npair_lspp, npair_lsmm);
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerDalitzMuMu, processAnalysis, "Process dalitz mumu for analysis", true);

  void processOnlyNmumu(soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent> const& collisions)
  {
    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        event_nmumu(0, 0, 0);
        continue;
      }
      event_nmumu(0, 0, 0);
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerDalitzMuMu, processOnlyNmumu, "Process only nmumu", false); // for central event filter processing
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerDalitzMuMu>(cfgc, TaskName{"skimmer-dalitz-mumu"})};
}
