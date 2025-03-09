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
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::pwgem::photonmeson;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEventsNgPCM>;
using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFbeta>;
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
  Configurable<float> minpt{"minpt", 0.05, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 2.0, "max eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 0.05, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 0.05, "max DCAz in cm"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.5, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", +3.5, "max. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 0.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", 0.0, "min. TPC n sigma for pion exclusion"};
  Configurable<float> minTOFNsigmaEl{"minTOFNsigmaEl", -3.5, "min. TOF n sigma for electron inclusion"};
  Configurable<float> maxTOFNsigmaEl{"maxTOFNsigmaEl", +3.5, "max. TOF n sigma for electron inclusion"};
  Configurable<float> maxMee{"maxMee", 0.03, "max. mee to store dalitz ee pairs"};
  Configurable<bool> fillLS{"fillLS", true, "flag to fill LS histograms for QA"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view dileptonSigns[3] = {"uls/", "lspp/", "lsmm/"};

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
    fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {400, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/hDCAxy_Pt", "DCA_{xy} vs. pT;p_{T} (GeV/c);DCA_{xy} (cm)", kTH2F, {{200, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/hDCAz_Pt", "DCA_{z} vs. pT;p_{T} (GeV/c);DCA_{z} (cm)", kTH2F, {{200, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);

    // TPC
    fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
    fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);

    // ITS
    fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);

    // TOF
    fRegistry.add("Track/hChi2TOF", "chi2 of TOF", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTOFbeta", "TOF beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
    fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);

    // pair
    fRegistry.add("Pair/uls/hMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{100, 0, 0.1}, {200, 0, 2}}, false);
    fRegistry.add("Pair/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{180, 0, M_PI}, {100, 0, 0.1}}, false);
    fRegistry.addClone("Pair/uls/", "Pair/lspp/");
    fRegistry.addClone("Pair/uls/", "Pair/lsmm/");
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
    if (track.hasTOF() && (track.tofNSigmaEl() < minTOFNsigmaEl || maxTOFNsigmaEl < track.tofNSigmaEl())) { // TOFif
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
                         track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
                         track.tpcChi2NCl(), track.tpcInnerParam(),
                         track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaPi(),
                         track.beta(), track.tofNSigmaEl(), track.tofNSigmaPi(),
                         track.itsClusterSizes(), track.itsChi2NCl(), track.tofChi2(), track.detectorMap(), track.tgl());

      fRegistry.fill(HIST("Track/hPt"), track.pt());
      fRegistry.fill(HIST("Track/hEtaPhi"), track.phi(), track.eta());
      fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/hDCAxy_Pt"), track.pt(), track.dcaXY());
      fRegistry.fill(HIST("Track/hDCAz_Pt"), track.pt(), track.dcaZ());
      fRegistry.fill(HIST("Track/hDCAxyzSigma"), track.dcaXY() / std::sqrt(track.cYY()), track.dcaZ() / std::sqrt(track.cZZ()));
      fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/hDCAzRes_Pt"), track.pt(), std::sqrt(track.cZZ()) * 1e+4);  // convert cm to um

      fRegistry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
      fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());

      fRegistry.fill(HIST("Track/hChi2TOF"), track.tofChi2());
      fRegistry.fill(HIST("Track/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());

      fRegistry.fill(HIST("Track/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());

      int nsize = 0;
      for (int il = 0; il < 7; il++) {
        nsize += track.itsClsSizeInLayer(il);
      }
      fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), track.p(), static_cast<float>(nsize) / static_cast<float>(track.itsNCls()) * std::cos(std::atan(track.tgl())));

      stored_trackIds.emplace_back(std::make_pair(collision.globalIndex(), track.globalIndex()));
    }
  }

  template <bool isMC, int pairtype, typename TCollision, typename TTracks1, typename TTracks2>
  void fillPairInfo(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    if constexpr (pairtype == 0) { // ULS
      for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack<isMC>(collision, t1) || !checkTrack<isMC>(collision, t2)) {
          continue;
        }
        if (!isElectron(t1) || !isElectron(t2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);
        fRegistry.fill(HIST("Pair/") + HIST(dileptonSigns[pairtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/") + HIST(dileptonSigns[pairtype]) + HIST("hMvsPhiV"), phiv, v12.M());

        if (v12.M() > maxMee) { // don't store
          continue;
        }
        fillTrackTable(collision, t1);
        fillTrackTable(collision, t2);
      } // end of pairing
    } else { // LS
      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack<isMC>(collision, t1) || !checkTrack<isMC>(collision, t2)) {
          continue;
        }
        if (!isElectron(t1) || !isElectron(t2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);
        fRegistry.fill(HIST("Pair/") + HIST(dileptonSigns[pairtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/") + HIST(dileptonSigns[pairtype]) + HIST("hMvsPhiV"), phiv, v12.M());
      } // end of pairing
    }
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

      fillPairInfo<false, 0>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (fillLS) {
        fillPairInfo<false, 1>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        fillPairInfo<false, 2>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
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

      fillPairInfo<true, 0>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (fillLS) {
        fillPairInfo<true, 1>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        fillPairInfo<true, 2>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectronFromDalitzEE, processMC, "process reconstructed and MC info ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerPrimaryElectronFromDalitzEE>(cfgc, TaskName{"skimmer-primary-electron-from-dalitzee"})};
}
