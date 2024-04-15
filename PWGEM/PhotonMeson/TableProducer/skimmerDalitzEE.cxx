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

/// \brief write relevant information for dalitz ee analysis to an AO2D.root file. This file is then the only necessary input to perform pcm analysis.
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

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsBz>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds>;
using MyTrack = MyTracks::iterator;

struct skimmerDalitzEE {
  enum class EM_EEPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSmm = -1,
  };

  SliceCache cache;
  Preslice<MyTracks> perCol = o2::aod::emprimaryelectron::emeventId;

  SliceCache cache_cefp;
  PresliceUnsorted<aod::EMPrimaryElectrons> perCol_cefp = o2::aod::emprimaryelectron::collisionId;

  Produces<aod::DalitzEEs> dalitzees;
  Produces<o2::aod::DalitzEEEMEventIds> dalitz_ee_eventid;
  Produces<o2::aod::EMEventsNee> event_nee;

  // Configurables
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<float> maxMee{"maxMee", 0.5, "max. mee to store ee pairs"};
  Configurable<bool> storeLS{"storeLS", false, "flag to store LS pairs"};
  Configurable<float> minpt{"minpt", 0.2, "min pt for track for loose track sample"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance for loose track sample"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 3.0, "max. TPC n sigma for electron inclusion"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 10, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 100, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<int> minitsncls{"minitsncls", 5, "min. number of ITS clusters"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 5.0, "max. chi2/NclsITS"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0f, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0f, "max DCAz in cm"};
  Configurable<float> dca_3d_sigma_max{"dca_3d_sigma_max", 1.f, "max DCA 3D in sigma"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hNpairs", "hNpairs;pair type;Number of Pairs", {HistType::kTH1F, {{3, -1.5f, +1.5f}}}},
      {"hNele", "hNele;centrality FT0C;Number of electrons", {HistType::kTH2F, {{110, 0, 110}, {101, -0.5f, +100.5f}}}},
      {"hNpos", "hNpos;centrality FT0C;Number of positrons", {HistType::kTH2F, {{110, 0, 110}, {101, -0.5f, +100.5f}}}},
    },
  };
  std::pair<int8_t, std::set<uint8_t>> itsRequirement = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.

  void init(InitContext const&) {}

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

  template <EM_EEPairType pairtype, bool isCEFP, typename TCollision, typename TTracks1, typename TTracks2>
  int fillPairTable(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    int npair = 0;
    if constexpr (pairtype == EM_EEPairType::kULS) { // ULS
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack(t1) || !checkTrack(t2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMee) { // don't store
          continue;
        }
        float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), collision.bz());
        float opangle = getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());

        // if (!std::isfinite(phiv)) {
        //   LOGF(info, "t1.px() = %f, t1.py() = %f, t1.pz() = %f, t2.px() = %f, t2.py() = %f, t2.pz() = %f", t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());
        // }

        if constexpr (isCEFP) {
          dalitzees(collision.globalIndex(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), v12.Rapidity(), phiv, opangle, static_cast<int>(pairtype));
          dalitz_ee_eventid(collision.globalIndex());
        } else { // for analysis
          dalitzees(collision.collisionId(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), v12.Rapidity(), phiv, opangle, static_cast<int>(pairtype));
          dalitz_ee_eventid(collision.globalIndex());
        }

        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        npair++;
      }      // end of pairing loop
    } else { // LS
      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack(t1) || !checkTrack(t2)) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMee) { // don't store
          continue;
        }
        float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), collision.bz());
        float opangle = getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());

        if constexpr (isCEFP) {
          dalitzees(collision.globalIndex(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), v12.M(), phiv, opangle, static_cast<int>(pairtype));
          dalitz_ee_eventid(collision.globalIndex());
        } else { // for analysis
          dalitzees(collision.collisionId(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), v12.Rapidity(), phiv, opangle, static_cast<int>(pairtype));
          dalitz_ee_eventid(collision.globalIndex());
        }

        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        npair++;
      } // end of pairing loop
    }
    return npair;
  }

  Partition<MyTracks> posTracks = o2::aod::emprimaryelectron::sign > 0 && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  Partition<MyTracks> negTracks = o2::aod::emprimaryelectron::sign < 0 && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  void processAnalysis(MyCollisions const& collisions, MyTracks const& tracks)
  {
    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        event_nee(0, 0, 0);
        continue;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      fRegistry.fill(HIST("hNpos"), collision.centFT0C(), posTracks_per_coll.size());
      fRegistry.fill(HIST("hNele"), collision.centFT0C(), negTracks_per_coll.size());
      // LOGF(info, "collision.centFT0C() = %f, posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", collision.centFT0C() , posTracks_per_coll.size(), negTracks_per_coll.size());

      int npair_uls = 0, npair_lspp = 0, npair_lsmm = 0;
      npair_uls = fillPairTable<EM_EEPairType::kULS, false>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        npair_lspp = fillPairTable<EM_EEPairType::kLSpp, false>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        npair_lsmm = fillPairTable<EM_EEPairType::kLSmm, false>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
      event_nee(npair_uls, npair_lspp, npair_lsmm);
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerDalitzEE, processAnalysis, "Process dalitz ee for analysis", true);

  Partition<aod::EMPrimaryElectrons> posTracks_cefp = o2::aod::emprimaryelectron::sign > 0 && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  Partition<aod::EMPrimaryElectrons> negTracks_cefp = o2::aod::emprimaryelectron::sign < 0 && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  void processCEFP(soa::Join<aod::Collisions, aod::EMEventsBz> const& collisions, aod::EMPrimaryElectrons const& tracks)
  {
    for (auto& collision : collisions) {
      auto posTracks_per_coll = posTracks_cefp->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache_cefp);
      auto negTracks_per_coll = negTracks_cefp->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache_cefp);

      int npair_uls = 0, npair_lspp = 0, npair_lsmm = 0;
      npair_uls = fillPairTable<EM_EEPairType::kULS, true>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        npair_lspp = fillPairTable<EM_EEPairType::kLSpp, true>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        npair_lsmm = fillPairTable<EM_EEPairType::kLSmm, true>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
      event_nee(npair_uls, npair_lspp, npair_lsmm);
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerDalitzEE, processCEFP, "Process dalitz ee for CEFP", false); // for central event filter processing
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerDalitzEE>(cfgc, TaskName{"skimmer-dalitz-ee"})};
}
