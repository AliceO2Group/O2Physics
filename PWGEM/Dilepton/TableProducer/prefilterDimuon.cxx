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
//
// ========================
//
// This code produces information on prefilter for dimuons.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TString.h"

#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

struct prefilterDimuon {
  using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
  using MyCollision = MyCollisions::iterator;

  using MyTracks = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds, aod::EMGlobalMuonSelfIds>;
  using MyTrack = MyTracks::iterator;

  Produces<aod::EMPrimaryMuonsPrefilterBitDerived> pfb_derived;

  // // Configurables
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  // Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
  } eventcuts;

  DimuonCut fDimuonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dimuoncut_group";

    // for deta-dphi prefilter
    Configurable<bool> cfg_apply_detadphi_uls{"cfg_apply_detadphi_uls", false, "flag to apply generator deta-dphi elliptic cut in ULS"}; // region to be rejected
    Configurable<bool> cfg_apply_detadphi_ls{"cfg_apply_detadphi_ls", false, "flag to apply generator deta-dphi elliptic cut in LS"};    // region to be rejected
    Configurable<float> cfg_min_deta_ls{"cfg_min_deta_ls", 0.04, "deta between 2 electrons (elliptic cut)"};                             // region to be rejected
    Configurable<float> cfg_min_dphi_ls{"cfg_min_dphi_ls", 0.2, "dphi between 2 electrons (elliptic cut)"};                              // region to be rejected
    Configurable<float> cfg_min_deta_uls{"cfg_min_deta_uls", 0.04, "deta between 2 electrons (elliptic cut)"};                           // region to be rejected
    Configurable<float> cfg_min_dphi_uls{"cfg_min_dphi_uls", 0.2, "dphi between 2 electrons (elliptic cut)"};                            // region to be rejected

    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+6, "max chi2/ndf"};
    Configurable<float> cfg_max_chi2mft{"cfg_max_chi2mft", 1e+6, "max chi2/ndf"};
    // Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 40, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_border_pt_for_chi2mchmft{"cfg_border_pt_for_chi2mchmft", 0, "border pt for different max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mftmch_lowPt{"cfg_max_matching_chi2_mftmch_lowPt", 8, "max chi2 for MFT-MCH matching for low pT"};
    Configurable<float> cfg_max_matching_chi2_mftmch_highPt{"cfg_max_matching_chi2_mftmch_highPt", 40, "max chi2 for MFT-MCH matching for high pT"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
    Configurable<float> cfg_max_relDPt_wrt_matchedMCHMID{"cfg_max_relDPt_wrt_matchedMCHMID", 1e+10f, "max. relative dpt between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_DEta_wrt_matchedMCHMID{"cfg_max_DEta_wrt_matchedMCHMID", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_DPhi_wrt_matchedMCHMID{"cfg_max_DPhi_wrt_matchedMCHMID", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
    Configurable<bool> requireMFTHitMap{"requireMFTHitMap", false, "flag to apply MFT hit map"};
    Configurable<std::vector<int>> requiredMFTDisks{"requiredMFTDisks", std::vector<int>{0}, "hit map on MFT disks [0,1,2,3,4]. logical-OR of each double-sided disk"};
  } dimuoncuts;

  // o2::ccdb::CcdbApi ccdbApi;
  // Service<o2::ccdb::BasicCCDBManager> ccdb;
  // int mRunNumber;
  // float d_bz;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext& /*context*/)
  {
    DefineEMEventCut();
    DefineDimuonCut();
    addhistograms();

    // mRunNumber = 0;
    // d_bz = 0;

    // ccdb->setURL(ccdburl);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setFatalWhenNull(false);
  }

  // template <typename TCollision>
  // void initCCDB(TCollision const& collision)
  // {
  //   if (mRunNumber == collision.runNumber()) {
  //     return;
  //   }

  //   // In case override, don't proceed, please - no CCDB access required
  //   if (d_bz_input > -990) {
  //     d_bz = d_bz_input;
  //     o2::parameters::GRPMagField grpmag;
  //     if (fabs(d_bz) > 1e-5) {
  //       grpmag.setL3Current(30000.f / (d_bz / 5.0f));
  //     }
  //     o2::base::Propagator::initFieldFromGRP(&grpmag);
  //     mRunNumber = collision.runNumber();
  //     return;
  //   }

  //   auto run3grp_timestamp = collision.timestamp();
  //   o2::parameters::GRPObject* grpo = 0x0;
  //   o2::parameters::GRPMagField* grpmag = 0x0;
  //   if (!skipGRPOquery)
  //     grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
  //   if (grpo) {
  //     o2::base::Propagator::initFieldFromGRP(grpo);
  //     // Fetch magnetic field from ccdb for current collision
  //     d_bz = grpo->getNominalL3Field();
  //     LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
  //   } else {
  //     grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
  //     if (!grpmag) {
  //       LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
  //     }
  //     o2::base::Propagator::initFieldFromGRP(grpmag);
  //     // Fetch magnetic field from ccdb for current collision
  //     d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
  //     LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
  //   }
  //   mRunNumber = collision.runNumber();
  // }

  ~prefilterDimuon() {}

  void addhistograms()
  {
    const AxisSpec axis_mass{380, 0.2, 4, "m_{#mu#mu} (GeV/c^{2})"};
    const AxisSpec axis_pair_pt{100, 0, 10, "p_{T,#mu#mu} (GeV/c)"};

    // for pair
    fRegistry.add("Pair/before/uls/hMvsPt", "m_{#mu#mu} vs. p_{T,#mu#mu}", kTH2D, {axis_mass, axis_pair_pt}, true);
    fRegistry.add("Pair/before/uls/hDeltaEtaDeltaPhi", "#Delta#eta-#Delta#varphi between 2 tracks;#Delta#varphi (rad.);#Delta#eta;", kTH2D, {{180, -M_PI, M_PI}, {400, -2, +2}}, true);
    fRegistry.addClone("Pair/before/uls/", "Pair/before/lspp/");
    fRegistry.addClone("Pair/before/uls/", "Pair/before/lsmm/");
    fRegistry.addClone("Pair/before/", "Pair/after/");
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(eventcuts.cfgZvtxMin, eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
  }

  void DefineDimuonCut()
  {
    fDimuonCut = DimuonCut("fDimuonCut", "fDimuonCut");

    // don't apply pair cut in prefilter!

    // for track
    fDimuonCut.SetTrackType(dimuoncuts.cfg_track_type);
    fDimuonCut.SetTrackPtRange(dimuoncuts.cfg_min_pt_track, dimuoncuts.cfg_max_pt_track);
    fDimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_eta_track, dimuoncuts.cfg_max_eta_track);
    fDimuonCut.SetTrackPhiRange(dimuoncuts.cfg_min_phi_track, dimuoncuts.cfg_max_phi_track);
    fDimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 20);
    fDimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    fDimuonCut.SetChi2MFT(0.f, dimuoncuts.cfg_max_chi2mft);
    // fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMaxMatchingChi2MCHMFTPtDep([&](float pt) { return (pt < dimuoncuts.cfg_border_pt_for_chi2mchmft ? dimuoncuts.cfg_max_matching_chi2_mftmch_lowPt : dimuoncuts.cfg_max_matching_chi2_mftmch_highPt); });
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
    fDimuonCut.SetMaxdPtdEtadPhiwrtMCHMID(dimuoncuts.cfg_max_relDPt_wrt_matchedMCHMID, dimuoncuts.cfg_max_DEta_wrt_matchedMCHMID, dimuoncuts.cfg_max_DPhi_wrt_matchedMCHMID); // this is relevant for global muons
    fDimuonCut.SetMFTHitMap(dimuoncuts.requireMFTHitMap, dimuoncuts.requiredMFTDisks);
    fDimuonCut.EnableTTCA(dimuoncuts.enableTTCA);
  }

  std::unordered_map<int, bool> map_best_match_globalmuon;
  std::unordered_map<int, uint16_t> map_pfb; // map track.globalIndex -> prefilter bit

  SliceCache cache;
  Preslice<MyTracks> perCollision_track = aod::emprimarymuon::emeventId;
  Partition<MyTracks> posTracks = o2::aod::emprimarymuon::sign > int8_t(0);
  Partition<MyTracks> negTracks = o2::aod::emprimarymuon::sign < int8_t(0);

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processPFB(FilteredMyCollisions const& collisions, MyTracks const& tracks)
  {
    map_best_match_globalmuon = findBestMatchMap(tracks, fDimuonCut);

    for (const auto& track : tracks) {
      map_pfb[track.globalIndex()] = 0;
    } // end of track loop

    for (const auto& collision : collisions) {
      // initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      bool is_cent_ok = true;
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        is_cent_ok = false;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);

      if (!fEMEventCut.IsSelected(collision) || !is_cent_ok) {
        for (const auto& pos : posTracks_per_coll) {
          map_pfb[pos.globalIndex()] = 0;
        }
        for (const auto& neg : negTracks_per_coll) {
          map_pfb[neg.globalIndex()] = 0;
        }
        continue;
      }

      // LOGF(info, "centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        if (!fDimuonCut.IsSelectedTrack(pos) || !fDimuonCut.IsSelectedTrack(neg)) {
          continue;
        }
        if (!map_best_match_globalmuon[pos.globalIndex()] || !map_best_match_globalmuon[neg.globalIndex()]) {
          continue;
        }

        // don't apply pair cut when you produce prefilter bit.

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(neg.pt(), neg.eta(), neg.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float deta = pos.sign() * v1.Pt() > neg.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = pos.sign() * v1.Pt() > neg.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/before/uls/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/before/uls/hDeltaEtaDeltaPhi"), dphi, deta);

        if (dimuoncuts.cfg_apply_detadphi_uls && std::pow(deta / dimuoncuts.cfg_min_deta_uls, 2) + std::pow(dphi / dimuoncuts.cfg_min_dphi_uls, 2) < 1.f) {
          map_pfb[pos.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS);
          map_pfb[neg.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS);
        }
      } // end of ULS pairing

      for (const auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        if (!fDimuonCut.IsSelectedTrack(pos1) || !fDimuonCut.IsSelectedTrack(pos2)) {
          continue;
        }
        if (!map_best_match_globalmuon[pos1.globalIndex()] || !map_best_match_globalmuon[pos2.globalIndex()]) {
          continue;
        }
        // don't apply pair cut when you produce prefilter bit.

        ROOT::Math::PtEtaPhiMVector v1(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float deta = pos1.sign() * v1.Pt() > pos2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = pos1.sign() * v1.Pt() > pos2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/before/lspp/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/before/lspp/hDeltaEtaDeltaPhi"), dphi, deta);

        if (dimuoncuts.cfg_apply_detadphi_ls && std::pow(deta / dimuoncuts.cfg_min_deta_ls, 2) + std::pow(dphi / dimuoncuts.cfg_min_dphi_ls, 2) < 1.f) {
          map_pfb[pos1.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS);
          map_pfb[pos2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS);
        }
      } // end of LS++ pairing

      for (const auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        if (!fDimuonCut.IsSelectedTrack(neg1) || !fDimuonCut.IsSelectedTrack(neg2)) {
          continue;
        }
        if (!map_best_match_globalmuon[neg1.globalIndex()] || !map_best_match_globalmuon[neg2.globalIndex()]) {
          continue;
        }
        // don't apply pair cut when you produce prefilter bit.

        ROOT::Math::PtEtaPhiMVector v1(neg1.pt(), neg1.eta(), neg1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(neg2.pt(), neg2.eta(), neg2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float deta = neg1.sign() * v1.Pt() > neg2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = neg1.sign() * v1.Pt() > neg2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/before/lsmm/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/before/lsmm/hDeltaEtaDeltaPhi"), dphi, deta);

        if (dimuoncuts.cfg_apply_detadphi_ls && std::pow(deta / dimuoncuts.cfg_min_deta_ls, 2) + std::pow(dphi / dimuoncuts.cfg_min_dphi_ls, 2) < 1.f) {
          map_pfb[neg1.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS);
          map_pfb[neg2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS);
        }
      } // end of LS-- pairing

    } // end of collision loop

    for (const auto& track : tracks) {
      // LOGF(info, "map_pfb[%d] = %d", track.globalIndex(), map_pfb[track.globalIndex()]);
      pfb_derived(map_pfb[track.globalIndex()]);
    } // end of track loop

    // check pfb.
    for (const auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        if (!fDimuonCut.IsSelectedTrack(pos) || !fDimuonCut.IsSelectedTrack(neg)) {
          continue;
        }
        if (!map_best_match_globalmuon[pos.globalIndex()] || !map_best_match_globalmuon[neg.globalIndex()]) {
          continue;
        }
        if (map_pfb[pos.globalIndex()] != 0 || map_pfb[neg.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(neg.pt(), neg.eta(), neg.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float deta = pos.sign() * v1.Pt() > neg.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = pos.sign() * v1.Pt() > neg.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/after/uls/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/after/uls/hDeltaEtaDeltaPhi"), dphi, deta);
      }

      for (const auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        if (!fDimuonCut.IsSelectedTrack(pos1) || !fDimuonCut.IsSelectedTrack(pos2)) {
          continue;
        }
        if (!map_best_match_globalmuon[pos1.globalIndex()] || !map_best_match_globalmuon[pos2.globalIndex()]) {
          continue;
        }
        if (map_pfb[pos1.globalIndex()] != 0 || map_pfb[pos2.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float deta = pos1.sign() * v1.Pt() > pos2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = pos1.sign() * v1.Pt() > pos2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/after/lspp/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/after/lspp/hDeltaEtaDeltaPhi"), dphi, deta);
      }

      for (const auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        if (!fDimuonCut.IsSelectedTrack(neg1) || !fDimuonCut.IsSelectedTrack(neg2)) {
          continue;
        }
        if (!map_best_match_globalmuon[neg1.globalIndex()] || !map_best_match_globalmuon[neg2.globalIndex()]) {
          continue;
        }
        if (map_pfb[neg1.globalIndex()] != 0 || map_pfb[neg2.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(neg1.pt(), neg1.eta(), neg1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(neg2.pt(), neg2.eta(), neg2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float deta = neg1.sign() * v1.Pt() > neg2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = neg1.sign() * v1.Pt() > neg2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/after/lsmm/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/after/lsmm/hDeltaEtaDeltaPhi"), dphi, deta);
      }

    } // end of collision loop
    map_pfb.clear();
    map_best_match_globalmuon.clear();
  } // end of process
  PROCESS_SWITCH(prefilterDimuon, processPFB, "produce prefilter bit", false);

  void processDummy(MyTracks const& tracks)
  {
    for (int i = 0; i < tracks.size(); i++) {
      pfb_derived(0);
    }
  }
  PROCESS_SWITCH(prefilterDimuon, processDummy, "dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<prefilterDimuon>(cfgc, TaskName{"prefilter-dimuon"})};
}
