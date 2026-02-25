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

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/GeometryManager.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep

#include <cmath>
#include <set>
#include <string>
#include <utility>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyCollisions = soa::Join<aod::EMEvents_004, aod::EMEventsMult_000, aod::EMEventsCent_000>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronDaEMEventIds>;
using MyTrack = MyTracks::iterator;

using MyTracksCEFP = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov>;
using MyTrackCEFP = MyTracksCEFP::iterator;

struct skimmerDalitzEE {
  enum class EM_EEPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSmm = -1,
  };

  SliceCache cache;
  Preslice<MyTracks> perCol = o2::aod::emprimaryelectronda::emphotoneventId;

  SliceCache cache_cefp;
  PresliceUnsorted<MyTracksCEFP> perCol_cefp = o2::aod::emprimaryelectron::collisionId;

  Produces<aod::DalitzEEs> dalitzees;
  Produces<o2::aod::DalitzEEEMEventIds> dalitz_ee_eventid;
  Produces<o2::aod::EMEventsNee> event_nee;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<float> maxMee{"maxMee", 1e+10, "max. mee to store ee pairs"};
  Configurable<bool> storeLS{"storeLS", false, "flag to store LS pairs"};
  Configurable<float> minpt{"minpt", 0.1, "min pt for track for loose track sample"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance for loose track sample"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.5, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 3.5, "max. TPC n sigma for electron inclusion"};
  Configurable<float> maxTOFNsigmaEl{"maxTOFNsigmaEl", 4.0, "max. TOF n sigma for electron inclusion"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 10, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 40, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<int> minitsncls{"minitsncls", 4, "min. number of ITS clusters"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0f, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0f, "max DCAz in cm"};
  Configurable<float> dca_3d_sigma_max{"dca_3d_sigma_max", 1e+10, "max DCA 3D in sigma"};
  Configurable<float> max_mean_itsob_cluster_size{"max_mean_itsob_cluster_size", 16.f, "max. <ITSob cluster size> x cos(lambda)"}; // this is to suppress random combination. default 4 + 1 for skimming.

  Configurable<bool> applyTPChadrejORTOFreq{"applyTPChadrejORTOFreq", false, "flag to apply TPChadrej-or-TOFreq at the skimming level"};
  Configurable<bool> applyPiRej_TPC{"applyPiRej_TPC", false, "flag to apply Pion rejection in TPC at the skimming level"};
  Configurable<bool> applyKaRej_TPC{"applyKaRej_TPC", false, "flag to apply Kaon rejection in TPC at the skimming level"};
  Configurable<bool> applyPrRej_TPC{"applyPrRej_TPC", false, "flag to apply Proton rejection in TPC at the skimming level"};
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 2.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> maxTPCNsigmaKa{"maxTPCNsigmaKa", 2.0, "max. TPC n sigma for kaon exclusion"};
  Configurable<float> maxTPCNsigmaPr{"maxTPCNsigmaPr", 2.0, "max. TPC n sigma for proton exclusion"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hNpairs", "hNpairs;pair type;Number of Pairs", {HistType::kTH1F, {{3, -1.5f, +1.5f}}}},
      {"hNele", "hNele;centrality FT0C;Number of electrons", {HistType::kTH2F, {{110, 0, 110}, {101, -0.5f, +100.5f}}}},
      {"hNpos", "hNpos;centrality FT0C;Number of positrons", {HistType::kTH2F, {{110, 0, 110}, {101, -0.5f, +100.5f}}}},
    },
  };

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  void init(InitContext const&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = collision.runNumber();
  }

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

    uint32_t itsClusterSizes = track.itsClusterSizes();
    int total_cluster_size = 0, nl = 0;
    for (unsigned int layer = 3; layer < 7; layer++) {
      int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
      if (cluster_size_per_layer > 0) {
        nl++;
      }
      total_cluster_size += cluster_size_per_layer;
    }
    if (static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(track.tgl())) > max_mean_itsob_cluster_size) {
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

    if (!isElectron(track)) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    if (applyTPChadrejORTOFreq) {
      return isElectron_TPChadrej(track) || isElectron_TOFrequire(track);
    } else {
      return true;
    }
    return true;
  }

  template <typename TTrack>
  bool isElectron_TPChadrej(TTrack const& track)
  {
    if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
      return false;
    }
    if (applyPiRej_TPC && std::abs(track.tpcNSigmaPi()) < maxTPCNsigmaPi) {
      return false;
    }
    if (applyKaRej_TPC && std::abs(track.tpcNSigmaKa()) < maxTPCNsigmaKa) {
      return false;
    }
    if (applyPrRej_TPC && std::abs(track.tpcNSigmaPr()) < maxTPCNsigmaPr) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron_TOFrequire(TTrack const& track)
  {
    if (applyPiRej_TPC && std::abs(track.tpcNSigmaPi()) < maxTPCNsigmaPi) {
      return false;
    }
    return minTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < maxTPCNsigmaEl && std::abs(track.tofNSigmaEl()) < maxTOFNsigmaEl;
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
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);
        float opangle = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());

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
      } // end of pairing loop
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
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);
        float opangle = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());

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

  Partition<MyTracks> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0) && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  Partition<MyTracks> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0) && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  void processAnalysis(MyCollisions const& collisions, MyTracks const&)
  {
    for (auto& collision : collisions) {
      initCCDB(collision);
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        event_nee(0, 0, 0);
        continue;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectronda::emphotoneventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectronda::emphotoneventId, collision.globalIndex(), cache);
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

  Partition<MyTracksCEFP> posTracks_cefp = o2::aod::emprimaryelectron::sign > int8_t(0) && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  Partition<MyTracksCEFP> negTracks_cefp = o2::aod::emprimaryelectron::sign < int8_t(0) && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  void processCEFP(aod::Collisions const& collisions, MyTracksCEFP const&, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
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

  void processOnlyNee(soa::Join<aod::EMEvents_004, aod::EMEventsMult_000, aod::EMEventsCent_000> const& collisions)
  {
    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        event_nee(0, 0, 0);
        continue;
      }
      event_nee(0, 0, 0);
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerDalitzEE, processOnlyNee, "Process only nee", false); // for central event filter processing
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerDalitzEE>(cfgc, TaskName{"skimmer-dalitz-ee"})};
}
