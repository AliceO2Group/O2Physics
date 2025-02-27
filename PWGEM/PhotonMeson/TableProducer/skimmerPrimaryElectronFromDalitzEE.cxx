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

/// \brief write relevant information about primary electrons.
/// \author daiki.sekihata@cern.ch

#include <map>
#include <string>
#include <vector>
#include <utility>

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::pwgem::photonmeson;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEventsNgPCM>;
using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov, aod::pidTPCFullEl, aod::pidTPCFullPi>;
using MyTrack = MyTracks::iterator;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;
using MyTrackMC = MyTracksMC::iterator;

struct skimmerPrimaryElectronFromDalitzEE {

  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  Produces<aod::EMPrimaryElectronsFromDalitz> emprimaryelectrons;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<bool> applyEveSel_at_skimming{"applyEveSel_at_skimming", false, "flag to apply minimal event selection at the skimming level"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
  Configurable<int> min_ncluster_its{"min_ncluster_its", 4, "min ncluster its"};
  Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 1, "min ncluster itsib"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};
  Configurable<float> minpt{"minpt", 0.1, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 0.05, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 0.05, "max DCAz in cm"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.5, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 3.5, "max. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 0.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", 0.0, "min. TPC n sigma for pion exclusion"};
  Configurable<float> maxMee{"maxMee", 0.06, "max. mee to store dalitz ee pairs"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {20, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/hDCAxy_Pt", "DCA_{xy} vs. pT;p_{T} (GeV/c);DCA_{xy} (cm)", kTH2F, {{1000, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/hDCAz_Pt", "DCA_{z} vs. pT;p_{T} (GeV/c);DCA_{z} (cm)", kTH2F, {{1000, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Pair/hMeePtee_ULS", "mee vs. pTee for dalitz ee ULS;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{100, 0, 0.1}, {100, 0, 10}}, false);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool checkTrack(TCollision const&, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (track.tpcChi2NCl() > maxchi2tpc) {
      return false;
    }

    if (track.itsChi2NCl() > maxchi2its) {
      return false;
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    if (track.itsNCls() < min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < min_ncluster_itsib) {
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

    if (track.tpcFractionSharedCls() > max_frac_shared_clusters_tpc) {
      return false;
    }

    if (std::fabs(track.dcaXY()) > dca_xy_max || std::fabs(track.dcaZ()) > dca_z_max) {
      return false;
    }

    if (track.pt() < minpt || std::fabs(track.eta()) > maxeta) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
      return false;
    }
    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi) {
      return false;
    }
    return true;
  }

  template <typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track)
  {
    if (std::find(stored_trackIds.begin(), stored_trackIds.end(), std::make_pair(collision.globalIndex(), track.globalIndex())) == stored_trackIds.end()) {
      emprimaryelectrons(collision.globalIndex(), track.globalIndex(), track.sign(),
                         track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(), track.cYY(), track.cZY(), track.cZZ(),
                         track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                         track.tpcChi2NCl(), track.tpcInnerParam(),
                         track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaPi(),
                         track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap(), track.tgl());

      fRegistry.fill(HIST("Track/hPt"), track.pt());
      fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/hEtaPhi"), track.phi(), track.eta());
      fRegistry.fill(HIST("Track/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/hDCAxy_Pt"), track.pt(), track.dcaXY());
      fRegistry.fill(HIST("Track/hDCAz_Pt"), track.pt(), track.dcaZ());
      fRegistry.fill(HIST("Track/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());

      stored_trackIds.emplace_back(std::make_pair(collision.globalIndex(), track.globalIndex()));
    }
  }

  template <bool isMC, typename TCollision, typename TTracks1, typename TTracks2>
  void fillPairInfo(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    for (const auto& t1 : tracks1) {
      for (const auto& t2 : tracks2) {
        if (!checkTrack<isMC>(collision, t1) || !checkTrack<isMC>(collision, t2)) {
          continue;
        }
        if (!isElectron(t1) || !isElectron(t2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        if (v12.M() > maxMee) { // don't store
          continue;
        }
        fRegistry.fill(HIST("Pair/hMeePtee_ULS"), v12.M(), v12.Pt());
        fillTrackTable(collision, t1);
        fillTrackTable(collision, t2);
      } // end of t2
    } // end of t1
  }

  std::vector<std::pair<int64_t, int64_t>> stored_trackIds;
  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& o2::aod::track::itsChi2NCl < maxchi2its&& ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true && nabs(o2::aod::track::dcaXY) < dca_xy_max&& nabs(o2::aod::track::dcaZ) < dca_z_max;
  Filter pidFilter = minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  using MyFilteredTracks = soa::Filtered<MyTracks>;
  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  // ---------- for data ----------
  void processRec(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    stored_trackIds.reserve(tracks.size());

    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (applyEveSel_at_skimming && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (collision.ngpcm() < 1) {
        continue;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      fillPairInfo<false>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectronFromDalitzEE, processRec, "process reconstructed info only", true); // standalone

  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyFilteredTracksMC> posTracksMC = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracksMC = o2::aod::track::signed1Pt < 0.f;
  // ---------- for MC ----------
  void processMC(MyCollisionsMC const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks)
  {
    stored_trackIds.reserve(tracks.size());

    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (applyEveSel_at_skimming && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (collision.ngpcm() < 1) {
        continue;
      }

      auto posTracks_per_coll = posTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      fillPairInfo<true>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectronFromDalitzEE, processMC, "process reconstructed and MC info ", false);
};

// struct associateAmbiguousElectron {
//   Produces<aod::EMAmbiguousElectronSelfIds> em_amb_ele_ids;
//
//   SliceCache cache;
//   PresliceUnsorted<aod::EMPrimaryElectrons> perTrack = o2::aod::emprimaryelectron::trackId;
//   std::vector<int> ambele_self_Ids;
//
//   void process(aod::EMPrimaryElectrons const& electrons)
//   {
//     for (const auto& electron : electrons) {
//       auto electrons_with_same_trackId = electrons.sliceBy(perTrack, electron.trackId());
//       ambele_self_Ids.reserve(electrons_with_same_trackId.size());
//       for (const auto& amp_ele : electrons_with_same_trackId) {
//         if (amp_ele.globalIndex() == electron.globalIndex()) { // don't store myself.
//           continue;
//         }
//         ambele_self_Ids.emplace_back(amp_ele.globalIndex());
//       }
//       em_amb_ele_ids(ambele_self_Ids);
//       ambele_self_Ids.clear();
//       ambele_self_Ids.shrink_to_fit();
//     }
//   }
// };
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryElectronFromDalitzEE>(cfgc, TaskName{"skimmer-primary-electron-from-dalitzee"}),
    //    adaptAnalysisTask<associateAmbiguousElectron>(cfgc, TaskName{"associate-ambiguous-electron"})
  };
}
