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

#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <tuple>

#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0PhotonsKFCov, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

struct prefilterPhoton {
  Produces<aod::V0PhotonsKFPrefilterBitDerived> pfb_derived;

  // Configurables
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

    // for mgg prefilter
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.12, "min mass for prefilter"}; // region to be rejected
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.15, "max mass for prefilter"}; // region to be rejected

    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.9, "min eta for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", +0.9, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 1.5, "max distance btween 2 legs"};
    Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfg_require_v0_with_correct_xz{"cfg_require_v0_with_correct_xz", true, "flag to select V0s with correct xz"};
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

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext& /*context*/)
  {
    DefineEMEventCut();
    DefinePCMCut();
    addhistograms();
  }

  ~prefilterPhoton() {}

  void addhistograms()
  {
    const AxisSpec axis_mass{200, 0, 0.8, "m_{#gamma#gamma} (GeV/c^{2})"};
    const AxisSpec axis_pair_pt{100, 0, 10, "p_{T,#gamma#gamma} (GeV/c)"};

    // for pair
    fRegistry.add("Pair/before/hMvsPt", "m_{#gamma#gamma} vs. p_{T,#gamma#gamma}", kTH2D, {axis_mass, axis_pair_pt}, true);
    fRegistry.addClone("Pair/before/", "Pair/after/");
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
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, 1e10f);
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
    fV0PhotonCut.SetIsWithinBeamPipe(pcmcuts.cfg_require_v0_with_correct_xz);
  }

  std::unordered_map<int, uint16_t> map_pfb; // map v0.globalIndex -> prefilter bit

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_v0 = aod::v0photonkf::emeventId;

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processPFB(FilteredMyCollisions const& collisions, MyV0Photons const& v0s, aod::V0Legs const&)
  {

    for (const auto& v0 : v0s) {
      map_pfb[v0.globalIndex()] = 0;
    } // end of v0 loop

    for (const auto& collision : collisions) {
      // initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      bool is_cent_ok = true;
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        is_cent_ok = false;
      }

      auto v0s_per_collision = v0s.sliceBy(perCollision_v0, collision.globalIndex());

      if (!fEMEventCut.IsSelected(collision) || !is_cent_ok) {
        for (const auto& v0 : v0s_per_collision) {
          map_pfb[v0.globalIndex()] = 0;
        }
        continue;
      }

      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(v0s_per_collision, v0s_per_collision))) {
        if (!fV0PhotonCut.template IsSelected<aod::V0Legs>(g1) || !fV0PhotonCut.template IsSelected<aod::V0Legs>(g2)) {
          continue;
        }
        // don't apply pair cut when you produce prefilter bit.

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.f);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.f);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        fRegistry.fill(HIST("Pair/before/hMvsPt"), v12.M(), v12.Pt());

        if (pcmcuts.cfg_min_mass < v12.M() && v12.M() < pcmcuts.cfg_max_mass) {
          map_pfb[g1.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::PhotonPrefilterBitDerived::kPhotonFromPi0Default);
          map_pfb[g2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::photonmeson::utils::pairutil::PhotonPrefilterBitDerived::kPhotonFromPi0Default);
        }
      }
    } // end of collision loop

    for (auto& v0 : v0s) {
      // LOGF(info, "map_pfb[%d] = %d", v0.globalIndex(), map_pfb[v0.globalIndex()]);
      pfb_derived(map_pfb[v0.globalIndex()]);
    } // end of v0 loop

    // check pfb.
    for (auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      auto v0s_per_collision = v0s.sliceBy(perCollision_v0, collision.globalIndex());

      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(v0s_per_collision, v0s_per_collision))) {
        if (!fV0PhotonCut.template IsSelected<aod::V0Legs>(g1) || !fV0PhotonCut.template IsSelected<aod::V0Legs>(g2)) {
          continue;
        }
        if (map_pfb[g1.globalIndex()] != 0 || map_pfb[g2.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.f);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.f);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/after/hMvsPt"), v12.M(), v12.Pt());
      }
    } // end of collision loop

    map_pfb.clear();

  } // end of process
  PROCESS_SWITCH(prefilterPhoton, processPFB, "produce prefilter bit", false);

  void processDummy(MyV0Photons const& v0s)
  {
    for (int i = 0; i < v0s.size(); i++) {
      pfb_derived(0);
    }
  }
  PROCESS_SWITCH(prefilterPhoton, processDummy, "dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<prefilterPhoton>(cfgc, TaskName{"prefilter-photon"})};
}
