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

/// \file PhotonHBT.h
/// \brief This code loops over v0 photons and makes pairs for photon HBT analysis.
/// \author Daiki Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_PHOTONHBT_H_
#define PWGEM_PHOTONMESON_CORE_PHOTONHBT_H_

#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
//
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h> // IWYU pragma: keep
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace o2::aod::pwgem::photon::core::photonhbt
{
enum class ggHBTPairType : int {
  kPCMPCM = 0,
};
} // namespace o2::aod::pwgem::photon::core::photonhbt

using namespace o2;                                      // o2-linter: disable=using-directive (required by O2 framework, inherited from upstream PWGEM)
using namespace o2::aod;                                 // o2-linter: disable=using-directive (required by O2 framework, inherited from upstream PWGEM)
using namespace o2::framework;                           // o2-linter: disable=using-directive (required by O2 framework, inherited from upstream PWGEM)
using namespace o2::framework::expressions;              // o2-linter: disable=using-directive (required by O2 framework, inherited from upstream PWGEM)
using namespace o2::soa;                                 // o2-linter: disable=using-directive (required by O2 framework, inherited from upstream PWGEM)
using namespace o2::aod::pwgem::dilepton::utils;         // o2-linter: disable=using-directive (required by O2 framework, inherited from upstream PWGEM)
using namespace o2::aod::pwgem::photon::core::photonhbt; // o2-linter: disable=using-directive (required by O2 framework, inherited from upstream PWGEM)

using MyCollisions = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsPhiVPsi>;
using MyV0Photon = MyV0Photons::iterator;

template <ggHBTPairType pairtype, typename... Types>
struct PhotonHBT {

  // Single-photon:
  //   0 = Inclusive  (all photons)
  //   1 = ITSTPC_ITSTPC   — both legs have ITS+TPC
  //   2 = ITSTPC_TPCOnly  — one ITS+TPC, one TPC-only
  //   3 = TPCOnly_TPCOnly — both legs TPC-only
  //
  // Pair combo.
  //   0 = Inclusive  (all recognised pairs)
  //   1 = ITSTPC_ITSTPC  × ITSTPC_ITSTPC
  //   2 = ITSTPC_ITSTPC  × ITSTPC_TPCOnly
  //   3 = ITSTPC_ITSTPC  × TPCOnly_TPCOnly
  //   4 = ITSTPC_TPCOnly × ITSTPC_TPCOnly
  //   5 = ITSTPC_TPCOnly × TPCOnly_TPCOnly
  //   6 = TPCOnly_TPCOnly × TPCOnly_TPCOnly

  template <typename TGamma, typename TSubInfos>
  static inline int classifyV0ComboIdx(TGamma const& g)
  {
    auto pos = g.template posTrack_as<TSubInfos>();
    auto neg = g.template negTrack_as<TSubInfos>();
    const bool posII = pos.hasITS() && pos.hasTPC();
    const bool posTPC = !pos.hasITS() && pos.hasTPC();
    const bool negII = neg.hasITS() && neg.hasTPC();
    const bool negTPC = !neg.hasITS() && neg.hasTPC();
    if (posII && negII)
      return 1;
    if ((posII && negTPC) || (posTPC && negII))
      return 2;
    if (posTPC && negTPC)
      return 3;
    return 0;
  }

  static inline int pairComboBin(int c1, int c2)
  {
    if (c1 <= 0 || c2 <= 0)
      return 0;
    if (c1 > c2)
      std::swap(c1, c2);
    static constexpr int kTable[4][4] = {
      {0, 0, 0, 0},
      {0, 1, 2, 3},
      {0, 2, 4, 5},
      {0, 3, 5, 6}};
    return kTable[c1][c2];
  }

  Configurable<bool> cfgDo3D{"cfgDo3D", false, "enable 3D analysis"};
  Configurable<int> cfgEp2EstimatorForMix{"cfgEp2EstimatorForMix", 3, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<int> cfgEp3EstimatorForMix{"cfgEp3EstimatorForMix", 3, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5 (for psi3 mixing)"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<float> maxY{"maxY", 0.8, "maximum rapidity for reconstructed particles"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  Configurable<uint64_t> ndiffBcMix{"ndiffBcMix", 594, "difference in global BC required in mixed events"};
  Configurable<bool> cfgUseLcms{"cfgUseLcms", true, "measure relative momentum in LCMS for 1D"};

  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis confCentBins{"confCentBins", {VARIABLE_WIDTH, 0.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis confEpBins{"confEpBins", {8, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}, "Mixing bins - psi2 event plane"};

  ConfigurableAxis confEp3Bins{"confEp3Bins", {1, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}, "Mixing bins - psi3 event plane (set to 1 bin to disable)"};
  ConfigurableAxis confOccupancyBins{"confOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  ConfigurableAxis confQBins{"confQBins", {60, 0, +0.3f}, "q bins for output histograms"}; // o2-linter: disable=name/configurable (Q is a physics symbol for momentum transfer; cannot be lowercased without losing meaning)
  ConfigurableAxis confKtBins{"confKtBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0}, "kT bins for output histograms"};

  // QA axis configurables
  ConfigurableAxis confPtBins{"confPtBins", {100, 0.f, 2.f}, "pT bins (GeV/c)"};
  ConfigurableAxis confEtaBins{"confEtaBins", {80, -0.8f, 0.8f}, "eta bins"};
  ConfigurableAxis confPhiBins{"confPhiBins", {90, -o2::constants::math::PI, o2::constants::math::PI}, "phi bins (rad)"};
  ConfigurableAxis confDeltaEtaBins{"confDeltaEtaBins", {100, -0.5f, +0.5f}, "Delta-eta bins"};
  ConfigurableAxis confDeltaPhiBins{"confDeltaPhiBins", {100, -0.5f, +0.5f}, "Delta-phi bins (rad)"};
  ConfigurableAxis confEllipseValBins{"confEllipseValBins", {200, 0.f, 10.f}, "ellipse value bins"};
  ConfigurableAxis confCosThetaBins{"confCosThetaBins", {100, 0.f, 1.f}, "cos(theta*) bins"};
  ConfigurableAxis confOpeningAngleBins{"confOpeningAngleBins", {100, 0.f, o2::constants::math::PI}, "opening angle bins (rad)"};

  ConfigurableAxis confRxyBins{"confRxyBins", {90, 0.f, 90.f}, "Rxy bins (cm)"};

  ConfigurableAxis confPsiPairBins{"confPsiPairBins", {100, -0.2f, 0.2f}, "psi_pair bins (rad)"};

  // ===========================================================================
  // QA flags
  // ===========================================================================

  Configurable<bool> doPairQa{"doPairQa", true, "fill pair QA histograms (Before/After ellipse cut, with combo axis)"};
  Configurable<bool> doSinglePhotonQa{"doSinglePhotonQa", true, "fill single-photon QA histograms (pT, eta, phi with V0 combo axis)"};

  // ===========================================================================
  // Pair cuts
  // ===========================================================================

  struct : ConfigurableGroup {
    std::string prefix = "ggpaircut_group";
    Configurable<float> cfgMinDrCosOa{"cfgMinDrCosOa", -1, "min. dr/cosOA for kPCMPCM"};
    Configurable<bool> cfgApplyEllipseCut{"cfgApplyEllipseCut", false, "reject pairs inside ellipse in DeltaEta-DeltaPhi"};
    Configurable<float> cfgEllipseSigEta{"cfgEllipseSigEta", 0.02f, "sigma_eta for ellipse cut"};
    Configurable<float> cfgEllipseSigPhi{"cfgEllipseSigPhi", 0.02f, "sigma_phi for ellipse cut"};
    Configurable<float> cfgEllipseR2{"cfgEllipseR2", 1.0f, "R^2 threshold: reject if value < R^2"};
  } ggpaircuts;

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND"};                  // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require no TF border"};              // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS ROF border"}; // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC"};             // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx FT0 vs PV"}; // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2.f, "min. FT0C occupancy"};         // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000.f, "max. FT0C occupancy"}; // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};   // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};         // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"}; // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "ITS layer 3 chips OK"};                                            // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "ITS layers 0-3 chips OK"};                                   // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "all ITS layers chips OK"};                                   // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfgRequireV0WithItstpc{"cfgRequireV0WithItstpc", false, "flag to select V0s with ITS-TPC matched tracks"};    // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireV0WithItsonly{"cfgRequireV0WithItsonly", false, "flag to select V0s with ITSonly tracks"};          // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRequireV0WithTpconly{"cfgRequireV0WithTpconly", false, "flag to select V0s with TPConly tracks"};          // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMinPtV0{"cfgMinPtV0", 0.1, "min pT for v0 photons at PV"};                                                // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxEtaV0{"cfgMaxEtaV0", 0.8, "max eta for v0 photons at PV"};                                             // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMinV0Radius{"cfgMinV0Radius", 16.0, "min v0 radius"};                                                     // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxV0Radius{"cfgMaxV0Radius", 90.0, "max v0 radius"};                                                     // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxAlphaAp{"cfgMaxAlphaAp", 0.95, "max alpha for AP cut"};                                                // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxQtAp{"cfgMaxQtAp", 0.01, "max qT for AP cut"};                                                         // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMinCospa{"cfgMinCospa", 0.997, "min V0 CosPA"};                                                           // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxPca{"cfgMaxPca", 3.0, "max distance btween 2 legs"};                                                   // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxChi2Kf{"cfgMaxChi2Kf", 1e+10, "max chi2/ndf with KF"};                                                 // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgRejectV0OnItsib{"cfgRejectV0OnItsib", true, "flag to reject V0s on ITSib"};                                // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgDisableItsonlyTrack{"cfgDisableItsonlyTrack", false, "flag to disable ITSonly tracks"};                    // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<bool> cfgDisableTpconlyTrack{"cfgDisableTpconlyTrack", false, "flag to disable TPConly tracks"};                    // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<int> cfgMinNclusterTpc{"cfgMinNclusterTpc", 0, "min ncluster tpc"};                                                 // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<int> cfgMinNcrossedrows{"cfgMinNcrossedrows", 40, "min ncrossed rows"};                                             // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxFracSharedClustersTpc{"cfgMaxFracSharedClustersTpc", 999.f, "max fraction of shared clusters in TPC"}; // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxChi2Tpc{"cfgMaxChi2Tpc", 4.0, "max chi2/NclsTPC"};                                                     // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxChi2Its{"cfgMaxChi2Its", 36.0, "max chi2/NclsITS"};                                                    // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMinTpcNsigmaEl{"cfgMinTpcNsigmaEl", -3.0, "min. TPC n sigma for electron"};                               // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
    Configurable<float> cfgMaxTpcNsigmaEl{"cfgMaxTpcNsigmaEl", +3.0, "max. TPC n sigma for electron"};                               // o2-linter: disable=name/configurable (upstream PWGEM name; JSON string cannot be changed without breaking existing configurations)
  } pcmcuts;

  ~PhotonHBT()
  {
    delete emh1;
    emh1 = 0x0;
    delete emh2;
    emh2 = 0x0;
    map_mixed_eventId_to_globalBC.clear();
    used_photonIds_per_col.clear();
    used_photonIds_per_col.shrink_to_fit();
  }

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::mt19937 engine;
  std::uniform_int_distribution<int> dist01;
  int mRunNumber;

  std::vector<float> zvtx_bin_edges;
  std::vector<float> cent_bin_edges;
  std::vector<float> ep2_bin_edges;
  std::vector<float> ep3_bin_edges;
  std::vector<float> occ_bin_edges;

  static constexpr float MinFloatTol = 1e-9f;
  static constexpr float CosineMin = -1.f;
  static constexpr float CosineMax = 1.f;
  static constexpr float HalfAngle = 2.f;
  static constexpr float PairEtaWeight = 0.5f;

  inline bool isInsideEllipse(float deta, float dphi) const
  {
    if (!ggpaircuts.cfgApplyEllipseCut.value)
      return false;
    const float sE = ggpaircuts.cfgEllipseSigEta.value;
    const float sP = ggpaircuts.cfgEllipseSigPhi.value;
    if (sE < MinFloatTol || sP < MinFloatTol)
      return false;
    return (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP) < ggpaircuts.cfgEllipseR2.value;
  }

  static inline float normDphi(float dphi)
  {
    while (dphi > o2::constants::math::PI)
      dphi -= o2::constants::math::TwoPI; // o2-linter: disable=two-pi-add-subtract (PWGEM has no PWGHF dependency; manual wrap is established pattern)
    while (dphi < -o2::constants::math::PI)
      dphi += o2::constants::math::TwoPI; // o2-linter: disable=two-pi-add-subtract (PWGEM has no PWGHF dependency; manual wrap is established pattern)
    return dphi;
  }

  static inline float computeCosTheta(const ROOT::Math::PtEtaPhiMVector& v1,
                                      const ROOT::Math::PtEtaPhiMVector& v2)
  {
    ROOT::Math::PxPyPzEVector p1(v1), p2(v2);
    ROOT::Math::PxPyPzEVector pair = p1 + p2;
    ROOT::Math::Boost boost(-pair.BoostToCM());
    ROOT::Math::PxPyPzEVector p1cm = boost(p1);
    ROOT::Math::XYZVector pairDir(pair.Px(), pair.Py(), pair.Pz());
    ROOT::Math::XYZVector p1cmDir(p1cm.Px(), p1cm.Py(), p1cm.Pz());
    if (pairDir.R() < MinFloatTol || p1cmDir.R() < MinFloatTol)
      return -1.f;
    return static_cast<float>(pairDir.Unit().Dot(p1cmDir.Unit()));
  }

  void init(InitContext& /*context*/)
  {
    mRunNumber = 0;

    auto parseBinsVerbatim = [](const ConfigurableAxis& cfg, std::vector<float>& edges) {
      if (cfg.value[0] == VARIABLE_WIDTH) {
        edges = std::vector<float>(cfg.value.begin(), cfg.value.end());
        edges.erase(edges.begin());
      } else {
        int nbins = static_cast<int>(cfg.value[0]);
        float xmin = static_cast<float>(cfg.value[1]);
        float xmax = static_cast<float>(cfg.value[2]);
        edges.resize(nbins + 1);
        for (int i = 0; i < nbins + 1; i++)
          edges[i] = (xmax - xmin) / nbins * i + xmin;
      }
    };

    parseBinsVerbatim(confVtxBins, zvtx_bin_edges);
    parseBinsVerbatim(confCentBins, cent_bin_edges);
    parseBinsVerbatim(confEpBins, ep2_bin_edges);
    parseBinsVerbatim(confEp3Bins, ep3_bin_edges);
    LOGF(info, "cfgOccupancyEstimator = %d", cfgOccupancyEstimator.value);
    parseBinsVerbatim(confOccupancyBins, occ_bin_edges);

    emh1 = new MyEMH(ndepth);
    emh2 = new MyEMH(ndepth);

    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    DefineEMEventCut();
    DefinePCMCut();
    addhistograms();

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_int_distribution<int>(0, 1);

    fRegistry.add("Pair/mix/hDiffBC",
                  "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|",
                  kTH1D, {{10001, -0.5, 10000.5}}, true);
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber())
      return;
    mRunNumber = collision.runNumber();
  }

  void addhistograms()
  {
    static constexpr std::string_view qvec_det_names[6] = {"FT0M", "FT0A", "FT0C", "BTot", "BPos", "BNeg"};
    fRegistry.add("Event/before/hEP2_CentFT0C_forMix",
                  Form("2nd harmonics event plane for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", qvec_det_names[cfgEp2EstimatorForMix].data()),
                  kTH2F, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);
    fRegistry.add("Event/after/hEP2_CentFT0C_forMix",
                  Form("2nd harmonics event plane for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", qvec_det_names[cfgEp2EstimatorForMix].data()),
                  kTH2F, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);

    const AxisSpec axis_kt{confKtBins, "k_{T} (GeV/c)"};
    const AxisSpec axis_qinv{confQBins, "q_{inv} (GeV/c)"};
    const AxisSpec axis_qabs_lcms{confQBins, "|#bf{q}|^{LCMS} (GeV/c)"};
    const AxisSpec axis_qout{confQBins, "q_{out} (GeV/c)"};
    const AxisSpec axis_qside{confQBins, "q_{side} (GeV/c)"};
    const AxisSpec axis_qlong{confQBins, "q_{long} (GeV/c)"};

    const AxisSpec axisPt{confPtBins, "p_{T} (GeV/c)"};
    const AxisSpec axisEta{confEtaBins, "#eta"};
    const AxisSpec axisPhi{confPhiBins, "#phi (rad)"};
    const AxisSpec axisDeltaEta{confDeltaEtaBins, "#Delta#eta"};
    const AxisSpec axisDeltaPhi{confDeltaPhiBins, "#Delta#phi (rad)"};
    const AxisSpec axisEllipseVal{confEllipseValBins, "Ellipse val"};
    const AxisSpec axisCosTheta{confCosThetaBins, "cos(#theta*)"};
    const AxisSpec axisOpeningAngle{confOpeningAngleBins, "#alpha (rad)"};

    const AxisSpec axisV0Combo{4, -0.5f, 3.5f, "V0 combo (0=Incl,1=II,2=IT,3=TT)"};
    const AxisSpec axisPairCombo{7, -0.5f, 6.5f, "Pair combo (0=Incl,1=II-II,2=II-IT,3=II-TT,4=IT-IT,5=IT-TT,6=TT-TT)"};
    const AxisSpec axisRxy{confRxyBins, "R_{xy} (cm)"};
    const AxisSpec axisPsiPair{confPsiPairBins, "#psi_{pair} (rad)"};

    // ── Single-photon QA ─────────────────────────────────────────────────────

    fRegistry.add("SinglePhoton/hPt", "V0 photon p_{T};p_{T} (GeV/c);V0 combo", kTH2F, {axisPt, axisV0Combo}, true);
    fRegistry.add("SinglePhoton/hEta", "V0 photon #eta;#eta;V0 combo", kTH2F, {axisEta, axisV0Combo}, true);
    fRegistry.add("SinglePhoton/hPhi", "V0 photon #phi;#phi (rad);V0 combo", kTH2F, {axisPhi, axisV0Combo}, true);
    fRegistry.add("SinglePhoton/hEtaVsPhi", "V0 photon acceptance;#phi (rad);#eta;V0 combo", kTH3F, {axisPhi, axisEta, axisV0Combo}, true);
    fRegistry.add("SinglePhoton/hRxy", "Conversion R_{xy};R_{xy} (cm);V0 combo", kTH2F, {axisRxy, axisV0Combo}, true);
    fRegistry.add("SinglePhoton/hPsiPair", "#psi_{pair};#psi_{pair} (rad);V0 combo", kTH2F, {axisPsiPair, axisV0Combo}, true);

    if (cfgDo3D) {
      fRegistry.add("Pair/same/hs_3d", "diphoton correlation 3D LCMS", kTHnSparseD, {axis_qout, axis_qside, axis_qlong, axis_kt}, true);
    } else {
      if (cfgUseLcms)
        fRegistry.add("Pair/same/hs_1d", "diphoton correlation 1D LCMS", kTHnSparseD, {axis_qabs_lcms, axis_kt}, true);
      else
        fRegistry.add("Pair/same/hs_1d", "diphoton correlation 1D", kTHnSparseD, {axis_qinv, axis_kt}, true);
    }
    if constexpr (pairtype == ggHBTPairType::kPCMPCM)
      fRegistry.add("Pair/same/hDeltaRCosOA", "distance between 2 conversion points;#Deltar/cos(#theta_{op}/2) (cm)", kTH1D, {{100, 0, 100}}, true);

    fRegistry.add("Pair/same/hKt", "k_{T};k_{T} (GeV/c);counts", kTH1F, {axis_kt}, true);

    for (const auto& step : {"Before", "After"}) {
      const std::string s = std::string("Pair/same/QA/") + step + "/";

      fRegistry.add((s + "hDeltaEta").c_str(), "#Delta#eta;#Delta#eta;Pair combo", kTH2F, {axisDeltaEta, axisPairCombo}, true);
      fRegistry.add((s + "hDeltaPhi").c_str(), "#Delta#phi;#Delta#phi (rad);Pair combo", kTH2F, {axisDeltaPhi, axisPairCombo}, true);
      fRegistry.add((s + "hDEtaDPhi").c_str(), "#Delta#eta vs #Delta#phi;#Delta#eta;#Delta#phi (rad);Pair combo", kTH3F, {axisDeltaEta, axisDeltaPhi, axisPairCombo}, true);
      fRegistry.add((s + "hDeltaEtaVsPairEta").c_str(), "#Delta#eta vs #LT#eta#GT_{pair};#LT#eta#GT_{pair};#Delta#eta;Pair combo", kTH3F, {axisEta, axisDeltaEta, axisPairCombo}, true);
      fRegistry.add((s + "hCosTheta").c_str(), "cos(#theta*);cos(#theta*);Pair combo", kTH2F, {axisCosTheta, axisPairCombo}, true);
      fRegistry.add((s + "hOpeningAngle").c_str(), "Opening angle;#alpha (rad);Pair combo", kTH2F, {axisOpeningAngle, axisPairCombo}, true);
      fRegistry.add((s + "hEllipseVal").c_str(), "Ellipse value;value;Pair combo", kTH2F, {axisEllipseVal, axisPairCombo}, true);
    }

    fRegistry.addClone("Pair/same/", "Pair/mix/");
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(eventcuts.cfgZvtxMin, eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
    fEMEventCut.SetRequireGoodITSLayer3(eventcuts.cfgRequireGoodITSLayer3);
    fEMEventCut.SetRequireGoodITSLayer0123(eventcuts.cfgRequireGoodITSLayer0123);
    fEMEventCut.SetRequireGoodITSLayersAll(eventcuts.cfgRequireGoodITSLayersAll);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfgMinPtV0, 1e10f);
    fV0PhotonCut.SetV0EtaRange(-pcmcuts.cfgMaxEtaV0, +pcmcuts.cfgMaxEtaV0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfgMinCospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfgMaxPca);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.cfgMaxChi2Kf);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfgMinV0Radius, pcmcuts.cfgMaxV0Radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfgMaxAlphaAp, pcmcuts.cfgMaxQtAp);
    fV0PhotonCut.RejectITSib(pcmcuts.cfgRejectV0OnItsib);
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfgMinNclusterTpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfgMinNcrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfgMaxFracSharedClustersTpc);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfgMaxChi2Tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfgMinTpcNsigmaEl, pcmcuts.cfgMaxTpcNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfgMaxChi2Its);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfgDisableItsonlyTrack);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.cfgDisableTpconlyTrack);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetRequireITSTPC(pcmcuts.cfgRequireV0WithItstpc);
    fV0PhotonCut.SetRequireITSonly(pcmcuts.cfgRequireV0WithItsonly);
    fV0PhotonCut.SetRequireTPConly(pcmcuts.cfgRequireV0WithTpconly);
  }

  template <int ev_id, typename TCollision>
  void fillPairHistogram(TCollision const&,
                         const ROOT::Math::PtEtaPhiMVector v1,
                         const ROOT::Math::PtEtaPhiMVector v2,
                         const float weight = 1.f)
  {
    float rndm = std::pow(-1, dist01(engine) % 2);
    ROOT::Math::PtEtaPhiMVector q12 = (v1 - v2) * rndm;
    ROOT::Math::PtEtaPhiMVector k12 = 0.5 * (v1 + v2);
    float qinv = -q12.M();
    float kt = k12.Pt();
    ROOT::Math::XYZVector uv_out(k12.Px() / k12.Pt(), k12.Py() / k12.Pt(), 0);
    ROOT::Math::XYZVector uv_long(0, 0, 1);
    ROOT::Math::XYZVector uv_side = uv_out.Cross(uv_long);
    ROOT::Math::PxPyPzEVector v1c(v1), v2c(v2);
    ROOT::Math::PxPyPzEVector q12c = (v1c - v2c) * rndm;
    float beta_z = (v1 + v2).Beta() * std::cos((v1 + v2).Theta());
    ROOT::Math::Boost bst_z(0, 0, -beta_z);
    ROOT::Math::PxPyPzEVector q12_lcms = bst_z(q12c);
    ROOT::Math::XYZVector q_3d_lcms = q12_lcms.Vect();
    float qabs_lcms = q_3d_lcms.R();
    float qout_lcms = q_3d_lcms.Dot(uv_out);
    float qside_lcms = q_3d_lcms.Dot(uv_side);
    float qlong_lcms = q_3d_lcms.Dot(uv_long);
    if (cfgDo3D)
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_3d"),
                     std::fabs(qout_lcms), std::fabs(qside_lcms), std::fabs(qlong_lcms), kt, weight);
    else {
      if (cfgUseLcms)
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_1d"), qabs_lcms, kt, weight);
      else
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_1d"), qinv, kt, weight);
    }

    fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hKt"), kt, weight);
  }

  template <int ev_id, bool IsBefore>
  inline void fillPairQAStep(float deta, float dphi, float pairEta,
                             float cosTheta, float openingAngle, float ellVal,
                             int pairComboIdx)
  {
    if (!doPairQa)
      return;

    if constexpr (IsBefore) {
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hDeltaEta"), deta, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hDeltaPhi"), dphi, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hDEtaDPhi"), deta, dphi, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hDeltaEtaVsPairEta"), pairEta, deta, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hCosTheta"), cosTheta, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hOpeningAngle"), openingAngle, 0.f);
      if (ellVal >= 0.f)
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hEllipseVal"), ellVal, 0.f);
      if (pairComboIdx > 0) {
        const float pcb = float(pairComboIdx);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hDeltaEta"), deta, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hDeltaPhi"), dphi, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hDEtaDPhi"), deta, dphi, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hDeltaEtaVsPairEta"), pairEta, deta, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hCosTheta"), cosTheta, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hOpeningAngle"), openingAngle, pcb);
        if (ellVal >= 0.f)
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/Before/hEllipseVal"), ellVal, pcb);
      }
    } else {
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hDeltaEta"), deta, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hDeltaPhi"), dphi, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hDEtaDPhi"), deta, dphi, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hDeltaEtaVsPairEta"), pairEta, deta, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hCosTheta"), cosTheta, 0.f);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hOpeningAngle"), openingAngle, 0.f);
      if (ellVal >= 0.f)
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hEllipseVal"), ellVal, 0.f);
      if (pairComboIdx > 0) {
        const float pcb = float(pairComboIdx);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hDeltaEta"), deta, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hDeltaPhi"), dphi, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hDEtaDPhi"), deta, dphi, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hDeltaEtaVsPairEta"), pairEta, deta, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hCosTheta"), cosTheta, pcb);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hOpeningAngle"), openingAngle, pcb);
        if (ellVal >= 0.f)
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("QA/After/hEllipseVal"), ellVal, pcb);
      }
    }
  }

  // ===========================================================================
  // runPairing
  // ===========================================================================

  template <typename TCollisions, typename TPhotons1, typename TPhotons2,
            typename TSubInfos1, typename TSubInfos2,
            typename TPreslice1, typename TPreslice2,
            typename TCut1, typename TCut2>
  void runPairing(TCollisions const& collisions,
                  TPhotons1 const& photons1, TPhotons2 const& photons2,
                  TSubInfos1 const&, TSubInfos2 const&,
                  TPreslice1 const& perCollision1, TPreslice2 const& perCollision2,
                  TCut1 const& cut1, TCut2 const& cut2)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      int ndiphoton = 0;

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator])
        continue;

      const float ep2_arr[6] = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(),
                                collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      const float ep3_arr[6] = {collision.ep3ft0m(), collision.ep3ft0a(), collision.ep3ft0c(),
                                collision.ep3btot(), collision.ep3bpos(), collision.ep3bneg()};
      float ep2 = ep2_arr[cfgEp2EstimatorForMix];
      float ep3 = ep3_arr[cfgEp3EstimatorForMix];

      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision))
        continue;
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      auto clampBin = [](int b, int nmax) {
        return (b < 0) ? 0 : (b > nmax ? nmax : b);
      };
      auto binOf = [&](const std::vector<float>& edges, float val) {
        int b = static_cast<int>(std::lower_bound(edges.begin(), edges.end(), val) - edges.begin()) - 1;
        return clampBin(b, static_cast<int>(edges.size()) - 2);
      };

      int zbin = binOf(zvtx_bin_edges, collision.posZ());
      int centbin = binOf(cent_bin_edges, centralities[cfgCentEstimator]);
      int ep2bin = binOf(ep2_bin_edges, ep2);
      int ep3bin = binOf(ep3_bin_edges, ep3);
      int occbin = binOf(occ_bin_edges,
                         cfgOccupancyEstimator == 1
                           ? static_cast<float>(collision.trackOccupancyInTimeRange())
                           : collision.ft0cOccupancyInTimeRange());

      auto key_bin = std::make_tuple(zbin, centbin, ep2bin, ep3bin, occbin);
      auto key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if constexpr (pairtype == ggHBTPairType::kPCMPCM) {

        auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_coll = photons2.sliceBy(perCollision2, collision.globalIndex());

        // ── Single-photon QA ─────────────────────────────────────────────────
        if (doSinglePhotonQa) {
          for (const auto& g : photons1_coll) {
            if (!cut1.template IsSelected<decltype(g), TSubInfos1>(g))
              continue;
            const int ci = classifyV0ComboIdx<decltype(g), TSubInfos1>(g); // 0=Other/1/2/3
            fRegistry.fill(HIST("SinglePhoton/hPt"), g.pt(), 0.f);
            fRegistry.fill(HIST("SinglePhoton/hEta"), g.eta(), 0.f);
            fRegistry.fill(HIST("SinglePhoton/hPhi"), g.phi(), 0.f);
            fRegistry.fill(HIST("SinglePhoton/hEtaVsPhi"), g.phi(), g.eta(), 0.f);
            fRegistry.fill(HIST("SinglePhoton/hRxy"), std::hypot(g.vx(), g.vy()), 0.f);
            if constexpr (requires { g.psipair(); }) {
              fRegistry.fill(HIST("SinglePhoton/hPsiPair"), g.psipair(), 0.f);
            }
            if (ci > 0) {
              const float cb = float(ci);
              fRegistry.fill(HIST("SinglePhoton/hPt"), g.pt(), cb);
              fRegistry.fill(HIST("SinglePhoton/hEta"), g.eta(), cb);
              fRegistry.fill(HIST("SinglePhoton/hPhi"), g.phi(), cb);
              fRegistry.fill(HIST("SinglePhoton/hEtaVsPhi"), g.phi(), g.eta(), cb);
              fRegistry.fill(HIST("SinglePhoton/hRxy"), std::hypot(g.vx(), g.vy()), cb);
              if constexpr (requires { g.psipair(); }) {
                fRegistry.fill(HIST("SinglePhoton/hPsiPair"), g.psipair(), cb);
              }
            }
          }
        }

        // ── Same-event pair loop ──────────────────────────────────────────────
        for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_coll, photons2_coll))) {
          if (!cut1.template IsSelected<decltype(g1), TSubInfos1>(g1) ||
              !cut2.template IsSelected<decltype(g2), TSubInfos2>(g2))
            continue;

          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          auto pos2 = g2.template posTrack_as<TSubInfos2>();
          auto ele2 = g2.template negTrack_as<TSubInfos2>();
          if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId())
            continue;

          const int c1 = classifyV0ComboIdx<decltype(g1), TSubInfos1>(g1);
          const int c2 = classifyV0ComboIdx<decltype(g2), TSubInfos2>(g2);
          const int pcb = pairComboBin(c1, c2); // 0 if either is Other

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);

          float dr = std::sqrt(std::pow(g1.vx() - g2.vx(), 2) + std::pow(g1.vy() - g2.vy(), 2) + std::pow(g1.vz() - g2.vz(), 2));
          ROOT::Math::XYZVector cp1(g1.vx(), g1.vy(), g1.vz());
          ROOT::Math::XYZVector cp2(g2.vx(), g2.vy(), g2.vz());
          float opa = std::acos(std::clamp(static_cast<float>(cp1.Dot(cp2) / (std::sqrt(cp1.Mag2()) * std::sqrt(cp2.Mag2()))), CosineMin, CosineMax));
          o2::math_utils::bringTo02Pi(opa);
          if (opa > o2::constants::math::PI)
            opa -= o2::constants::math::PI;
          float cosOA = std::cos(opa / HalfAngle);
          if (dr / cosOA < ggpaircuts.cfgMinDrCosOa)
            continue;
          fRegistry.fill(HIST("Pair/same/hDeltaRCosOA"), dr / cosOA);

          float deta = g1.eta() - g2.eta();
          float dphi = normDphi(g1.phi() - g2.phi());
          float pairEta = PairEtaWeight * (g1.eta() + g2.eta());
          float cosTheta = std::fabs(computeCosTheta(v1, v2));
          const float sE = ggpaircuts.cfgEllipseSigEta.value;
          const float sP = ggpaircuts.cfgEllipseSigPhi.value;
          float ellVal = (sE > MinFloatTol && sP > MinFloatTol)
                           ? (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP)
                           : -1.f;

          // QA Before ellipse cut
          fillPairQAStep<0, true>(deta, dphi, pairEta, cosTheta, opa, ellVal, pcb);

          // Ellipse cut
          if (isInsideEllipse(deta, dphi))
            continue;

          // QA After ellipse cut
          fillPairQAStep<0, false>(deta, dphi, pairEta, cosTheta, opa, ellVal, pcb);

          fillPairHistogram<0>(collision, v1, v2, 1.f);
          ndiphoton++;

          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g1.globalIndex()) == used_photonIds_per_col.end()) {
            EMPair g1tmp = EMPair(g1.pt(), g1.eta(), g1.phi(), 0);
            g1tmp.setConversionPointXYZ(g1.vx(), g1.vy(), g1.vz());
            emh1->AddTrackToEventPool(key_df_collision, g1tmp);
            used_photonIds_per_col.emplace_back(g1.globalIndex());
          }
          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g2.globalIndex()) == used_photonIds_per_col.end()) {
            EMPair g2tmp = EMPair(g2.pt(), g2.eta(), g2.phi(), 0);
            g2tmp.setConversionPointXYZ(g2.vx(), g2.vy(), g2.vz());
            emh1->AddTrackToEventPool(key_df_collision, g2tmp);
            used_photonIds_per_col.emplace_back(g2.globalIndex());
          }
        } // end same-event pair loop
      }

      used_photonIds_per_col.clear();
      used_photonIds_per_col.shrink_to_fit();

      // ── Mixed-event loop ────────────────────────────────────────────────────
      if (!cfgDoMix || !(ndiphoton > 0))
        continue;

      auto selected_photons1_in_this_event = emh1->GetTracksPerCollision(key_df_collision);
      auto collisionIds1_in_mixing_pool = emh1->GetCollisionIdsFromEventPool(key_bin);

      if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
        for (const auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;
          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId)
            continue;

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
          if (diffBC < ndiffBcMix)
            continue;

          auto photons1_from_event_pool = emh1->GetTracksPerCollision(mix_dfId_collisionId);

          for (const auto& g1 : selected_photons1_in_this_event) {
            for (const auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);

              float dr = std::sqrt(std::pow(g1.vx() - g2.vx(), 2) + std::pow(g1.vy() - g2.vy(), 2) + std::pow(g1.vz() - g2.vz(), 2));
              ROOT::Math::XYZVector cp1(g1.vx(), g1.vy(), g1.vz());
              ROOT::Math::XYZVector cp2(g2.vx(), g2.vy(), g2.vz());
              float opa = std::acos(std::clamp(static_cast<float>(cp1.Dot(cp2) / (std::sqrt(cp1.Mag2()) * std::sqrt(cp2.Mag2()))), CosineMin, CosineMax));
              o2::math_utils::bringTo02Pi(opa);
              if (opa > o2::constants::math::PI)
                opa -= o2::constants::math::PI;
              float cosOA = std::cos(opa / HalfAngle);
              if (dr / cosOA < ggpaircuts.cfgMinDrCosOa)
                continue;
              fRegistry.fill(HIST("Pair/mix/hDeltaRCosOA"), dr / cosOA);

              float deta = g1.eta() - g2.eta();
              float dphi = normDphi(g1.phi() - g2.phi());
              float pairEta = PairEtaWeight * (g1.eta() + g2.eta());
              float cosTheta = std::fabs(computeCosTheta(v1, v2));
              const float sE = ggpaircuts.cfgEllipseSigEta.value;
              const float sP = ggpaircuts.cfgEllipseSigPhi.value;
              float ellVal = (sE > MinFloatTol && sP > MinFloatTol)
                               ? (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP)
                               : -1.f;

              fillPairQAStep<1, true>(deta, dphi, pairEta, cosTheta, opa, ellVal, 0);

              if (isInsideEllipse(deta, dphi))
                continue;

              fillPairQAStep<1, false>(deta, dphi, pairEta, cosTheta, opa, ellVal, 0);

              fillPairHistogram<1>(collision, v1, v2, 1.f);
            }
          }
        } // end mixed event pool loop
      }

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh2->AddCollisionIdAtLast(key_bin, key_df_collision);
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
      }
    } // end collision loop
  }

  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<
    std::tuple<int, int, int, int, int>,
    std::pair<int, int>, EMPair>;
  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;
  std::vector<int> used_photonIds_per_col;
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::pmeventId;

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) ||
                                      (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) ||
                                      (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange &&
                                           o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange &&
                                          o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  int ndf = 0;
  void processAnalysis(FilteredMyCollisions const& collisions, Types const&... args)
  {
    if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      runPairing(collisions, v0photons, v0photons, v0legs, v0legs,
                 perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut);
    }
    ndf++;
  }
  PROCESS_SWITCH(PhotonHBT, processAnalysis, "pairing for analysis", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(PhotonHBT, processDummy, "Dummy function", true);
};

#endif // PWGEM_PHOTONMESON_CORE_PHOTONHBT_H_
