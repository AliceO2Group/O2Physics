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

#include <cstdlib>
#include <tuple>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Math/Boost.h"
#include "Math/Vector4D.h"
#include "TVector3.h"

// \struct RealPhotonHBTTask
/// \brief Simple real photon HBT Task to exract necessary V0 and cluster information
/// \author Stefanie Mrozinski <stefanie.mrozinski@cern.ch>, Goethe University Frankfurt
/// \since 03.02.2025
///
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyEMEvents = soa::Join<aod::EMEvents, aod::EMEventsCent>;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::BinnedGenPts>;
using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

struct TaskRealPhotonHBT {

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<bool> cfg_require_v0_on_wwire_ib{"cfg_require_v0_on_wwire_ib", false, "flag to select V0s on W wires ITSib"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.8, "min eta for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", +0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.997, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 3.0, "max distance btween 2 legs"};
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
  } pcmcuts;

  struct : ConfigurableGroup {
    std::string prefix = "mixingConfig";
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgCentBins{"cfgCentBins", {VARIABLE_WIDTH, 0.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Mixing bins - centrality"};
    Configurable<int> cfgMixingDepth{"cfgMixingDepth", 2, "Mixing depth"};
  } mixingConfig;

  HistogramRegistry mHistManager{"mHistManager", {}, OutputObjHandlingPolicy::AnalysisObject, false};

  ConfigurableAxis KBinning{"K", {800, 0.0f, 1.0f}, "Binning used along K axis for average pair momentum"};
  ConfigurableAxis QinvBinning{"QInv", {800, 0.0f, 2.0f}, "Binning used along QInv axis"};

  // creating histograms

  using o2HistType = HistType;
  using o2Axis = AxisSpec;

  // photon related histograms

  void init(InitContext const&)
  {

    mHistManager.add("Same/QinvVsK", "Qinv and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
    mHistManager.add("Same/qlcmsVsK", "qout and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
    mHistManager.add("Same/qoutVsK", "qout and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
    mHistManager.add("Same/qlongVsK", "qlong and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
    mHistManager.add("Same/qsideVsK", "qside and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});

    mHistManager.add("Mixed/QinvVsK", "Qinv and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
    mHistManager.add("Mixed/qlcmsVsK", "qout and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
    mHistManager.add("Mixed/qoutVsK", "qout and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
    mHistManager.add("Mixed/qlongVsK", "qlong and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
    mHistManager.add("Mixed/qsideVsK", "qside and average pair momentum of V0 candidates", o2HistType::kTH2D, {QinvBinning, KBinning});
  }

  tuple<float, float, float, float, float, float> getPairProperties(auto const& v01, auto const& v02)
  {
    ROOT::Math::PxPyPzEVector v1(v01.px(), v01.py(), v01.pz(), v01.e());
    ROOT::Math::PxPyPzEVector v2(v02.px(), v02.py(), v02.pz(), v02.e());

    ROOT::Math::PxPyPzEVector photonpair = v1 + v2;

    float Qinv = std::sqrt(std::max(0.0, photonpair.M2()));
    float kt = 0.5 * (v01.pt() + v02.pt());

    double betaLong = photonpair.pz() / photonpair.e();

    // Create the boost transformation along the Z-axis
    ROOT::Math::BoostZ boostZ(-betaLong);

    // Apply the boost to transform the four-momentum vectors
    ROOT::Math::PxPyPzEVector p1LCMS(v01.px(), v01.py(), v01.pz(), v01.e());
    ;
    p1LCMS = boostZ(p1LCMS);
    ROOT::Math::PxPyPzEVector p2LCMS(v02.px(), v02.py(), v02.pz(), v02.e());
    ;
    p1LCMS = boostZ(p2LCMS);

    // Compute relative momentum q in LCMS
    ROOT::Math::PxPyPzEVector q4 = p1LCMS - p2LCMS;

    float qlong = q4.Pz();
    float q_lcms = q4.P();
    float qTmag = std::sqrt(q4.Px() * q4.Px() + q4.Py() * q4.Py());
    float pTmag = std::sqrt(photonpair.Px() * photonpair.Px() + photonpair.Py() * photonpair.Py());
    float qout = (pTmag > 0) ? (q4.Px() * photonpair.Px() + q4.Py() * photonpair.Py()) / pTmag : 0;
    float qside = (pTmag > 0) ? std::sqrt(qTmag * qTmag - qout * qout) : qTmag;

    return {kt, q_lcms, Qinv, qside, qlong, qout};
  }

  void processSameEvents(MyEMEvents const& collisions, MyV0Photons const& v0photons)
  {

    for (const auto& collision : collisions) {

      auto v0photons_collision = v0photons.sliceBy(perCollision, collision.globalIndex());

      for (auto& [v01, v02] : combinations(CombinationsStrictlyUpperIndexPolicy(v0photons_collision, v0photons_collision))) {

        auto [kt, q_lcms, Qinv, qside, qlong, qout] = getPairProperties(v01, v02);

        mHistManager.fill(HIST("Same/QinvVsK"), Qinv, kt);
        mHistManager.fill(HIST("Same/qlcmsVsK"), q_lcms, kt);
        mHistManager.fill(HIST("Same/qsideVsK"), qside, kt);
        mHistManager.fill(HIST("Same/qlongVsK"), qlong, kt);
        mHistManager.fill(HIST("Same/qoutVsK"), qout, kt);
      }
    }
  }

  PROCESS_SWITCH(TaskRealPhotonHBT, processSameEvents, "Process Same Events", false);

  SliceCache cache;

  void processMixedEvents(MyEMEvents const& collisions, MyV0Photons const& v0photons)
  {
    auto getPhotonCount =
      [&v0photons, this](MyEMEvents::iterator const& col) {
        auto associatedPhotons = v0photons.sliceByCached(v0photonkf::emeventId, col.globalIndex(), this->cache);
        return associatedPhotons.size();
      };

    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getPhotonCount)>, aod::collision::PosZ>;
    BinningType binningWithLambda{{getPhotonCount}, {mixingConfig.cfgVtxBins}, true};

    auto photonsTuple = std::make_tuple(v0photons);

    SameKindPair<MyEMEvents, MyV0Photons, BinningType> pair{binningWithLambda, 5, -1, collisions, photonsTuple, &cache};

    for (auto& [c1, v0photons1, c2, v0photons2] : pair) {

      for (auto& [v01, v02] : combinations(CombinationsFullIndexPolicy(v0photons1, v0photons2))) {
        auto [kt, q_lcms, Qinv, qside, qlong, qout] = getPairProperties(v01, v02);
        mHistManager.fill(HIST("Mixed/QinvVsK"), Qinv, kt);
        mHistManager.fill(HIST("Mixed/qlcmsVsK"), q_lcms, kt);
        mHistManager.fill(HIST("Mixed/qsideVsK"), qside, kt);
        mHistManager.fill(HIST("Mixed/qlongVsK"), qlong, kt);
        mHistManager.fill(HIST("Mixed/qoutVsK"), qout, kt);
      }
    }
  }

  PROCESS_SWITCH(TaskRealPhotonHBT, processMixedEvents, "Process Mixed Events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskRealPhotonHBT>(cfgc)};
}
