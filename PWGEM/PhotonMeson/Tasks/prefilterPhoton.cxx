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
// This code produces information on prefilter for photon.
//    Please write to: daiki.sekihata@cern.ch

#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

// #include "TString.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
// #include "Common/Core/RecoDecay.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"

#include "Common/Core/trackUtilities.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::photonpair;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyPrimaryElectrons = soa::Join<aod::EMPrimaryElectronsFromDalitz, aod::EMPrimaryElectronEMEventIds>;
using MyPrimaryElectron = MyPrimaryElectrons::iterator;

struct prefilterPhoton {
  Produces<aod::V0PhotonsKFPrefilterBitDerived> pfb_v0_derived;
  Produces<aod::EMPrimaryElectronsPrefilterBitDerived> pfb_ele_derived;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
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

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";

    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_pt_v0{"cfg_max_pt_v0", 1e+10, "max pT for v0 photons at PV"};
    Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.9, "min eta for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", +0.9, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.99, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 1.5, "max distance btween 2 legs"};
    Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};
  } pcmcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";

    Configurable<float> cfg_min_mee{"cfg_min_mee", 0.0, "min mass"};
    Configurable<float> cfg_max_mee{"cfg_max_mee", 0.02, "max mass"};
    // Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", false, "flag to apply phiv cut"};
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.05, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -2.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 2.0, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 40, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 0, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.f, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.f, "max dca Z for single track in cm"};
    Configurable<float> cfg_max_dca3dsigma_track{"cfg_max_dca3dsigma_track", 1.5, "max DCA 3D in sigma"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTOFif), "pid scheme [kTOFif : 0, kTPConly : 1]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", 0.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", 0.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
  } dileptoncuts;

  struct : ConfigurableGroup {
    std::string prefix = "ggcut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.10, "min mass for prefilter"}; // region to be rejected
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.15, "max mass for prefilter"}; // region to be rejected
  } ggcuts;

  struct : ConfigurableGroup {
    std::string prefix = "eegcut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.10, "min mass for prefilter"}; // region to be rejected
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.15, "max mass for prefilter"}; // region to be rejected
  } eegcuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;

  void init(InitContext& /*context*/)
  {
    DefineEMEventCut();
    DefinePCMCut();
    addhistograms();

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  ~prefilterPhoton() {}

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

  void addhistograms()
  {
    const AxisSpec axis_mass{200, 0, 0.8, "m_{#gamma#gamma} (GeV/c^{2})"};
    const AxisSpec axis_pair_pt{100, 0, 10, "p_{T,#gamma#gamma} (GeV/c)"};
    const AxisSpec axis_phiv{180, 0, M_PI, "#varphi_{V} (rad.)"};

    // for pair
    fRegistry.add("Pair/PCMPCM/before/hMvsPt", "m_{#gamma#gamma} vs. p_{T,#gamma#gamma}", kTH2D, {axis_mass, axis_pair_pt}, true);
    fRegistry.add("Pair/PCMDalitzEE/before/hMvsPt", "m_{ee#gamma} vs. p_{T,ee#gamma}", kTH2D, {axis_mass, axis_pair_pt}, true);
    fRegistry.add("Pair/PCMDalitzEE/before/hMvsPhiV", "m_{ee} vs. #varphi_{V}", kTH2D, {{180, 0, M_PI}, {100, 0, 0.1}}, true);
    fRegistry.addClone("Pair/PCMPCM/before/", "Pair/PCMPCM/after/");
    fRegistry.addClone("Pair/PCMDalitzEE/before/", "Pair/PCMDalitzEE/after/");
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, pcmcuts.cfg_max_pt_v0);
    fV0PhotonCut.SetV0EtaRange(pcmcuts.cfg_min_eta_v0, pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfg_max_pca);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.cfg_max_chi2kf);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfg_min_v0radius, pcmcuts.cfg_max_v0radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfg_max_alpha_ap, pcmcuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(pcmcuts.cfg_reject_v0_on_itsib);

    // for track
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfg_min_ncluster_tpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfg_max_frac_shared_clusters_tpc);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfg_disable_itsonly_track);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.cfg_disable_tpconly_track);

    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
  }

  void DefineDileptonCut()
  {
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    // for pair
    fDileptonCut.SetMeeRange(dileptoncuts.cfg_min_mee, dileptoncuts.cfg_max_mee);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfg_phiv_intercept) / dileptoncuts.cfg_phiv_slope; });
    fDileptonCut.ApplyPhiV(false);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfg_require_itsib_any);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfg_require_itsib_1st);

    // for track
    fDileptonCut.SetTrackPtRange(dileptoncuts.cfg_min_pt_track, dileptoncuts.cfg_max_pt_track);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfg_min_eta_track, +dileptoncuts.cfg_max_eta_track);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfg_min_ncluster_tpc);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfg_min_ncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetMaxFracSharedClustersTPC(dileptoncuts.cfg_max_frac_shared_clusters_tpc);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfg_max_chi2tpc);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfg_max_chi2its);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfg_min_ncluster_its, 7);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfg_max_dcaxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfg_max_dcaz);
    fDileptonCut.SetTrackDca3DRange(0.f, dileptoncuts.cfg_max_dca3dsigma_track); // in sigma

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);
  }

  template <PairType pairtype, typename TCollisions, typename TPhotons1, typename TPhotons2, typename TSubInfos1, typename TSubInfos2, typename TPreslice1, typename TPreslice2, typename TCut1, typename TCut2>
  void runPairing(TCollisions const& collisions,
                  TPhotons1 const& photons1, TPhotons2 const& photons2,
                  TSubInfos1 const&, TSubInfos2 const&,
                  TPreslice1 const& perCollision1, TPreslice2 const& perCollision2,
                  TCut1 const& cut1, TCut2 const& cut2)
  {
    if constexpr (pairtype == PairType::kPCMPCM) {
      for (const auto& photon1 : photons1) {
        map_pfb_v0[photon1.globalIndex()] = 0;
      } // end of v0 loop
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      for (const auto& photon1 : photons1) {
        map_pfb_v0[photon1.globalIndex()] = 0;
      } // end of v0 loop
      for (const auto& photon2 : photons2) {
        map_pfb_ele[photon2.globalIndex()] = 0;
      } // end of electron loop
    }

    if constexpr (pairtype == PairType::kPCMPCM) {
      for (const auto& collision : collisions) {
        initCCDB(collision);
        const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
        bool is_cent_ok = true;
        if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
          is_cent_ok = false;
        }

        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex());

        if (!fEMEventCut.IsSelected(collision) || !is_cent_ok) {
          for (const auto& photon1 : photons1_per_collision) {
            map_pfb_v0[photon1.globalIndex()] = 0;
          }
          continue;
        }
        for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }
          // don't apply pair cut when you produce prefilter bit.

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.f);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.f);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          fRegistry.fill(HIST("Pair/PCMPCM/before/hMvsPt"), v12.M(), v12.Pt());

          if (ggcuts.cfg_min_mass < v12.M() && v12.M() < ggcuts.cfg_max_mass) {
            map_pfb_v0[g1.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::PhotonPrefilterBitDerived::kPhotonFromPi0gg);
            map_pfb_v0[g2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::PhotonPrefilterBitDerived::kPhotonFromPi0gg);
          }
        } // end of 2photon pairing loop
      } // end of collision loop
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      for (const auto& collision : collisions) {
        initCCDB(collision);
        const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
        bool is_cent_ok = true;
        if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
          is_cent_ok = false;
        }

        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto positrons_per_collision = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto electrons_per_collision = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

        if (!fEMEventCut.IsSelected(collision) || !is_cent_ok) {
          for (const auto& photon1 : photons1_per_collision) {
            map_pfb_v0[photon1.globalIndex()] = 0;
          }
          for (const auto& pos : positrons_per_collision) {
            map_pfb_ele[pos.globalIndex()] = 0;
          }
          for (const auto& ele : electrons_per_collision) {
            map_pfb_ele[ele.globalIndex()] = 0;
          }
          continue;
        }

        for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_collision, photons1_per_collision))) { // PCM-PCM // cut, and subinfo is different from kPCMPCM
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut1.template IsSelected<TSubInfos1>(g2)) {
            continue;
          }
          // don't apply pair cut when you produce prefilter bit.

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.f);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.f);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          fRegistry.fill(HIST("Pair/PCMPCM/before/hMvsPt"), v12.M(), v12.Pt());

          if (ggcuts.cfg_min_mass < v12.M() && v12.M() < ggcuts.cfg_max_mass) {
            map_pfb_v0[g1.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::PhotonPrefilterBitDerived::kPhotonFromPi0gg);
            map_pfb_v0[g2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::PhotonPrefilterBitDerived::kPhotonFromPi0gg);
          }
        } // end of 2photon pairing loop

        for (const auto& g1 : photons1_per_collision) { // PCM-DalitzEE
          if (!cut1.template IsSelected<TSubInfos1>(g1)) {
            continue;
          }
          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

          for (const auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
            if (pos2.trackId() == ele2.trackId()) { // this is protection against pairing identical 2 tracks.
              continue;
            }
            if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) {
              continue;
            }

            if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
              continue;
            }

            ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;
            if (!(dileptoncuts.cfg_min_mee < v_ee.M() && v_ee.M() < dileptoncuts.cfg_max_mee)) {
              continue;
            }
            ROOT::Math::PtEtaPhiMVector veeg = v_gamma + v_pos + v_ele;
            fRegistry.fill(HIST("Pair/PCMDalitzEE/before/hMvsPt"), veeg.M(), veeg.Pt());

            if (eegcuts.cfg_min_mass < veeg.M() && veeg.M() < eegcuts.cfg_max_mass) {
              map_pfb_v0[g1.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::PhotonPrefilterBitDerived::kPhotonFromPi0eeg);
              map_pfb_ele[pos2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::ElectronPrefilterBitDerived::kElectronFromPi0eeg);
              map_pfb_ele[ele2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::ElectronPrefilterBitDerived::kElectronFromPi0eeg);
            }
          } // end of dielectron loop
        } // end of g1 loop

        for (const auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
          if (pos2.trackId() == ele2.trackId()) { // this is protection against pairing identical 2 tracks.
            continue;
          }

          if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos2.px(), pos2.py(), pos2.pz(), ele2.px(), ele2.py(), ele2.pz(), pos2.sign(), ele2.sign(), d_bz);
          fRegistry.fill(HIST("Pair/PCMDalitzEE/before/hMvsPhiV"), phiv, v_ee.M());

          if (v_ee.M() < phiv * dileptoncuts.cfg_phiv_slope + dileptoncuts.cfg_phiv_intercept) {
            map_pfb_ele[pos2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::ElectronPrefilterBitDerived::kElectronFromFakePC);
            map_pfb_ele[ele2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::ElectronPrefilterBitDerived::kElectronFromFakePC);
          }
        } // end of dielectron loop to reject photon conversion
      } // end of collision loop
    }

    if constexpr (pairtype == PairType::kPCMPCM) {
      for (const auto& photon1 : photons1) {
        pfb_v0_derived(map_pfb_v0[photon1.globalIndex()]);
      } // end of v0 loop
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      for (const auto& photon1 : photons1) {
        pfb_v0_derived(map_pfb_v0[photon1.globalIndex()]);
      } // end of v0 loop
      for (const auto& photon2 : photons2) {
        pfb_ele_derived(map_pfb_ele[photon2.globalIndex()]);
      } // end of electron loop
    }

    // check pfb.
    if constexpr (pairtype == PairType::kPCMPCM) {
      for (auto& collision : collisions) {
        const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
        if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
          continue;
        }

        if (!fEMEventCut.IsSelected(collision)) {
          continue;
        }

        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex());

        for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }
          if (map_pfb_v0[g1.globalIndex()] != 0 || map_pfb_v0[g2.globalIndex()] != 0) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.f);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.f);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          fRegistry.fill(HIST("Pair/PCMPCM/after/hMvsPt"), v12.M(), v12.Pt());
        }
      } // end of collision loop
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      for (auto& collision : collisions) {
        initCCDB(collision);
        const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
        if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
          continue;
        }

        if (!fEMEventCut.IsSelected(collision)) {
          continue;
        }

        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto positrons_per_collision = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto electrons_per_collision = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

        for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_collision, photons1_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut1.template IsSelected<TSubInfos1>(g2)) {
            continue;
          }
          if (map_pfb_v0[g1.globalIndex()] != 0 || map_pfb_v0[g2.globalIndex()] != 0) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.f);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.f);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          fRegistry.fill(HIST("Pair/PCMPCM/after/hMvsPt"), v12.M(), v12.Pt());
        }

        for (const auto& g1 : photons1_per_collision) {
          if (!cut1.template IsSelected<TSubInfos1>(g1)) {
            continue;
          }
          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

          for (const auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
            if (pos2.trackId() == ele2.trackId()) { // this is protection against pairing identical 2 tracks.
              continue;
            }
            if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) {
              continue;
            }

            if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
              continue;
            }
            if (map_pfb_v0[g1.globalIndex()] != 0 || map_pfb_ele[pos2.globalIndex()] != 0 || map_pfb_ele[ele2.globalIndex()] != 0) {
              continue;
            }

            ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;
            float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos2.px(), pos2.py(), pos2.pz(), ele2.px(), ele2.py(), ele2.pz(), pos2.sign(), ele2.sign(), d_bz);
            if (!(dileptoncuts.cfg_min_mee < v_ee.M() && v_ee.M() < dileptoncuts.cfg_max_mee)) {
              continue;
            }
            ROOT::Math::PtEtaPhiMVector veeg = v_gamma + v_pos + v_ele;
            fRegistry.fill(HIST("Pair/PCMDalitzEE/after/hMvsPt"), veeg.M(), veeg.Pt());
            fRegistry.fill(HIST("Pair/PCMDalitzEE/after/hMvsPhiV"), phiv, v_ee.M());
          } // end of dielectron loop
        } // end of g1 loop
      } // end of collision loop
    }

    map_pfb_v0.clear();
    map_pfb_ele.clear();
  }

  std::unordered_map<int, uint16_t> map_pfb_v0;  // map v0.globalIndex -> prefilter bit
  std::unordered_map<int, uint16_t> map_pfb_ele; // map ele.globalIndex -> prefilter bit

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_v0 = aod::v0photonkf::emeventId;
  Preslice<MyPrimaryElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Partition<MyPrimaryElectrons> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<MyPrimaryElectrons> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0);

  void processPCMPCM(FilteredMyCollisions const& collisions, MyV0Photons const& v0s, aod::V0Legs const& v0legs)
  {
    runPairing<PairType::kPCMPCM>(collisions, v0s, v0s, v0legs, v0legs, perCollision_v0, perCollision_v0, fV0PhotonCut, fV0PhotonCut); // produces filter bit for both photons
  }
  PROCESS_SWITCH(prefilterPhoton, processPCMPCM, "produce prefilter bit for PCM-PCM", false);

  void processPCMDalitzEE(FilteredMyCollisions const& collisions, MyV0Photons const& v0s, aod::V0Legs const& v0legs, MyPrimaryElectrons const& primaryelectrons)
  {
    runPairing<PairType::kPCMDalitzEE>(collisions, v0s, primaryelectrons, v0legs, primaryelectrons, perCollision_v0, perCollision_electron, fV0PhotonCut, fDileptonCut); // produces filter bit for both photons and electrons
  }
  PROCESS_SWITCH(prefilterPhoton, processPCMDalitzEE, "produce prefilter bit for PCM-DalitzEE", false);

  void processDummyV0(MyV0Photons const& v0s)
  {
    for (int i = 0; i < v0s.size(); i++) {
      pfb_v0_derived(0);
    }
  }
  PROCESS_SWITCH(prefilterPhoton, processDummyV0, "dummy for v0s", true);

  void processDummyElectron(MyPrimaryElectrons const& primaryelectrons)
  {
    for (int i = 0; i < primaryelectrons.size(); i++) {
      pfb_ele_derived(0);
    }
  }
  PROCESS_SWITCH(prefilterPhoton, processDummyElectron, "dummy for electrons", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<prefilterPhoton>(cfgc, TaskName{"prefilter-photon"})};
}
