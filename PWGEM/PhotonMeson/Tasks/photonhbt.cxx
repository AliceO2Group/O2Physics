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

/// \file photonhbt.cxx
/// \brief V0 photon HBT analysis.
/// \author Daiki Sekihata, daiki.sekihata@cern.ch
///         Stefanie Mrozinski, stefanie.mrozinski@cern.ch

#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

#include "Common/Core/RecoDecay.h"
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
#include <Math/Vector3D.h> // IWYU pragma: keep (do not replace with Math/Vector3Dfwd.h)
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <initializer_list>
#include <map>
#include <random>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

/// Single-photon track-type combo.
enum class V0Combo : int {
  Inclusive = 0,
  ItstpcItstpc = 1,
  ItstpcTpconly = 2,
  TpconlyTpconly = 3,
};

/// Photon-pair track-type combo.
enum class PairCombo : int {
  Inclusive = 0,
  IiXIi = 1,
  IiXIt = 2,
  IiXTt = 3,
  ItXIt = 4,
  ItXTt = 5,
  TtXTt = 6,
};

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils;

// EMMCEventLabels needed for processMC truth-efficiency loop
using MyCollisions = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000,
                               aod::EMEventsCent_000, aod::EMEventsQvec_001,
                               aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsPhiVPsi>;
using MyV0Photon = MyV0Photons::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

// ─── MC truth classification types ────────────────────────────────────────────

struct PhotonMCInfo {
  bool hasMC = false;
  bool sameMother = false;
  bool isTruePhoton = false;
  int mcPosId = -1;
  int mcNegId = -1;
  int motherId = -1;
  int motherPdg = 0;
  bool isPhysicalPrimary = false;
};

enum class PairTruthType : uint8_t {
  Unknown = 0,
  TrueTrueDistinct,
  TrueTrueSamePhoton,
  SharedMcLeg,
  TrueFake,
  FakeFake,
  Pi0Daughters,
};

static constexpr float kMinMagnitude = 1e-12f;
static constexpr float kMinCosine = 1e-12f;
static constexpr float kMinSigma = 1e-9;

struct Photonhbt {

  template <is_iterator TGamma, is_table TSubInfos>
  static inline V0Combo classifyV0Combo(TGamma const& g)
  {
    const auto pos = g.template posTrack_as<TSubInfos>();
    const auto neg = g.template negTrack_as<TSubInfos>();
    const bool posII = pos.hasITS() && pos.hasTPC();
    const bool posTPC = !pos.hasITS() && pos.hasTPC();
    const bool negII = neg.hasITS() && neg.hasTPC();
    const bool negTPC = !neg.hasITS() && neg.hasTPC();
    if (posII && negII)
      return V0Combo::ItstpcItstpc;
    if ((posII && negTPC) || (posTPC && negII))
      return V0Combo::ItstpcTpconly;
    if (posTPC && negTPC)
      return V0Combo::TpconlyTpconly;
    return V0Combo::Inclusive;
  }

  static inline PairCombo classifyPairCombo(V0Combo c1, V0Combo c2)
  {
    const int i1 = static_cast<int>(c1);
    const int i2 = static_cast<int>(c2);
    if (i1 <= 0 || i2 <= 0)
      return PairCombo::Inclusive;
    const int lo = std::min(i1, i2);
    const int hi = std::max(i1, i2);
    static constexpr std::array<std::array<int, 4>, 4> kTable = {
      {{0, 0, 0, 0}, {0, 1, 2, 3}, {0, 2, 4, 5}, {0, 3, 5, 6}}};
    return static_cast<PairCombo>(kTable[lo][hi]);
  }

  // ─── Configurables ────────────────────────────────────────────────────────

  ConfigurableAxis confQBins{"confQBins", {60, 0, +0.3f}, "q bins for output histograms"};
  ConfigurableAxis confKtBins{"confKtBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75}, "kT bins"};
  ConfigurableAxis confPtBins{"confPtBins", {100, 0.f, 2.f}, "pT bins (GeV/c)"};
  ConfigurableAxis confEtaBins{"confEtaBins", {80, -0.8f, 0.8f}, "eta bins"};
  ConfigurableAxis confPhiBins{"confPhiBins", {90, 0.f, o2::constants::math::TwoPI}, "phi bins (rad)"};
  ConfigurableAxis confDeltaEtaBins{"confDeltaEtaBins", {180, -1.6f, +1.6f}, "Delta-eta bins"};
  ConfigurableAxis confDeltaPhiBins{"confDeltaPhiBins", {180, -o2::constants::math::PI, o2::constants::math::PI}, "Delta-phi bins (rad)"};
  ConfigurableAxis confEllipseValBins{"confEllipseValBins", {200, 0.f, 10.f}, "ellipse value bins"};
  ConfigurableAxis confCosThetaBins{"confCosThetaBins", {100, 0.f, 1.f}, "cos(theta*) bins"};
  ConfigurableAxis confOpeningAngleBins{"confOpeningAngleBins", {100, 0.f, o2::constants::math::PI}, "opening angle bins (rad)"};
  ConfigurableAxis confRBins{"confRBins", {100, 0.f, 100.f}, "conversion radius bins (cm)"};
  ConfigurableAxis confDeltaRBins{"confDeltaRBins", {120, 0.f, 30.f}, "|R1-R2| bins (cm)"};
  ConfigurableAxis confDeltaR3DBins{"confDeltaR3DBins", {100, 0.f, 100.f}, "3D distance between conversion points (cm)"};
  ConfigurableAxis confDeltaRxyBins{"confDeltaRxyBins", {100, 0.f, 100.f}, "xy distance between conversion points (cm)"};
  ConfigurableAxis confZConvBins{"confZConvBins", {200, -100.f, 100.f}, "conversion z (cm)"};
  ConfigurableAxis confDeltaZBins{"confDeltaZBins", {200, -100.f, 100.f}, "#Deltaz bins (cm)"};
  ConfigurableAxis confOccupancyQA{"confOccupancyQA", {100, 0.f, 50000.f}, "occupancy"};
  ConfigurableAxis confCentQABins{"confCentQABins", {110, 0.f, 110.f}, "centrality (%)"};

  // ─── Axis specs ────────────────────────────────────────────────────────────

  const AxisSpec axisKt{confKtBins, "k_{T} (GeV/c)"};
  const AxisSpec axisQinv{confQBins, "q_{inv} (GeV/c)"};
  const AxisSpec axisQabsLcms{confQBins, "|#bf{q}|^{LCMS} (GeV/c)"};
  const AxisSpec axisQout{confQBins, "q_{out} (GeV/c)"};
  const AxisSpec axisQside{confQBins, "q_{side} (GeV/c)"};
  const AxisSpec axisQlong{confQBins, "q_{long} (GeV/c)"};
  const AxisSpec axisPt{confPtBins, "p_{T} (GeV/c)"};
  const AxisSpec axisEta{confEtaBins, "#eta"};
  const AxisSpec axisPhi{confPhiBins, "#phi (rad)"};
  const AxisSpec axisDeltaEta{confDeltaEtaBins, "#Delta#eta"};
  const AxisSpec axisDeltaPhi{confDeltaPhiBins, "#Delta#phi (rad)"};
  const AxisSpec axisEllipseVal{confEllipseValBins, "(#Delta#eta/#sigma_{#eta})^{2}+(#Delta#phi/#sigma_{#phi})^{2}"};
  const AxisSpec axisCosTheta{confCosThetaBins, "cos(#theta*)"};
  const AxisSpec axisOpeningAngle{confOpeningAngleBins, "Opening angle (rad)"};
  const AxisSpec axisR{confRBins, "R_{conv} (cm)"};
  const AxisSpec axisDeltaR{confDeltaRBins, "|R_{1}-R_{2}| (cm)"};
  const AxisSpec axisDeltaR3D{confDeltaR3DBins, "|#vec{r}_{1}-#vec{r}_{2}| (cm)"};
  const AxisSpec axisDeltaRxy{confDeltaRxyBins, "#Delta r_{xy} (cm)"};
  const AxisSpec axisZConv{confZConvBins, "z_{conv} (cm)"};
  const AxisSpec axisDeltaZ{confDeltaZBins, "#Delta z (cm)"};
  const AxisSpec axisOccupancy{confOccupancyQA, "occupancy"};
  const AxisSpec axisCentQA{confCentQABins, "centrality (%)"};

  // ─── Configurables: QA flags ───────────────────────────────────────────────

  struct : ConfigurableGroup {
    std::string prefix = "qaflags_group";
    Configurable<bool> doPairQa{"doPairQa", true, "fill pair QA histograms at each cut step"};
    Configurable<bool> doSinglePhotonQa{"doSinglePhotonQa", true, "fill single-photon QA histograms (pT, eta, phi)"};
    Configurable<float> cfgMaxQinvForQA{"cfgMaxQinvForQA", 0.1f, "fill per-step pair QA histograms only when q_inv < this value. Set <= 0 to disable."};
    Configurable<float> cfgMaxQinvForFullRange{"cfgMaxQinvForFullRange", 0.3f, "fill full-range histograms only when q_inv < this value. Set <= 0 to disable."};
    Configurable<float> cfgMaxQinvForMCQA{"cfgMaxQinvForMCQA", 0.3f,
                                          "fill MC truth 1D histograms (hQinv, hKt, hDeltaEta, ...) only when q_inv < this value. "
                                          "hDEtaDPhi is always filled (needs full sample). Set <= 0 to disable. Default 0.6 cuts "
                                          "most combinatorics while covering well beyond the CF range for systematics."};
  } qaflags;

  // ─── HBT analysis mode ───────────────────────────────────────────────────────────
  struct : ConfigurableGroup {
    std::string prefix = "hbtanalysis_group";
    Configurable<bool> cfgDo3D{"cfgDo3D", false, "enable 3D (qout,qside,qlong) analysis"};
    Configurable<bool> cfgDo2D{"cfgDo2D", false, "enable 2D (qout,qinv) projection (requires cfgDo3D)"};
    Configurable<bool> cfgUseLCMS{"cfgUseLCMS", false, "measure 1D relative momentum in LCMS"};
  } hbtanalysis;

  // ─── Event mixing ─────────────────────────────────────────────────────────────
  struct : ConfigurableGroup {
    std::string prefix = "mixing_group";
    Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
    Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
    Configurable<uint64_t> ndiffBCMix{"ndiffBCMix", 594, "difference in global BC required for mixed events"};
    Configurable<int> cfgEP2EstimatorForMix{"cfgEP2EstimatorForMix", 3, "FT0M:0, FT0A:1, FT0C:2, FV0A:3, BTot:4, BPos:5, BNeg:6"};
    Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
    Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};

    ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis confCentBins{"confCentBins", {VARIABLE_WIDTH, 0.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f, 999.f}, "Mixing bins - centrality"};
    ConfigurableAxis confEPBinsBins{"confEPBinsBins", {16, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}, "Mixing bins - EP angle"};
    ConfigurableAxis confOccupancyBins{"confOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};
  } mixing;

  // ─── Centrality slection ─────────────────────────────────────────────────
  struct : ConfigurableGroup {
    std::string prefix = "centralitySelection_group";
    Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  } centralitySelection;

  struct : ConfigurableGroup {
    std::string prefix = "mctruth_group";
    Configurable<float> cfgMCMaxQinv{"cfgMCMaxQinv", 0.3f, "..."};
    Configurable<float> cfgMCMinKt{"cfgMCMinKt", 0.0f, "..."};
    Configurable<float> cfgMCMaxKt{"cfgMCMaxKt", 0.7f, "..."};
    Configurable<bool> cfgDoTruthMix{"cfgDoTruthMix", false, "..."};
    Configurable<int> cfgTruthMixDepth{"cfgTruthMixDepth", 10, "..."};
    Configurable<float> cfgMCMinV0Pt{"cfgMCMinV0Pt", 0.1f,
                                     "min pT for true photons in truth-efficiency loop (GeV/c); "
                                     "0 = fall back to pcmcuts.cfgMinPtV0"};
    Configurable<float> cfgMCMinLegPt{"cfgMCMinLegPt", 0.0f, "min pT for true e^{+}/e^{-} legs in truth-efficiency loop (GeV/c);"};
  } mctruth;

  struct : ConfigurableGroup {
    std::string prefix = "mctruthSparse_group";
    Configurable<bool> cfgFillDEtaDPhiVsQinvTrueTrueDistinct{"cfgFillDEtaDPhiVsQinvTrueTrueDistinct", true, "fill hDEtaDPhiVsQinv for TrueTrueDistinct pairs"};
    Configurable<bool> cfgFillDEtaDPhiVsQinvTrueTrueSamePhoton{"cfgFillDEtaDPhiVsQinvTrueTrueSamePhoton", false, "fill hDEtaDPhiVsQinv for TrueTrueSamePhoton pairs"};
    Configurable<bool> cfgFillDEtaDPhiVsQinvSharedMcLeg{"cfgFillDEtaDPhiVsQinvSharedMcLeg", false, "fill hDEtaDPhiVsQinv for SharedMcLeg pairs"};
    Configurable<bool> cfgFillDEtaDPhiVsQinvTrueFake{"cfgFillDEtaDPhiVsQinvTrueFake", false, "fill hDEtaDPhiVsQinv for TrueFake pairs"};
    Configurable<bool> cfgFillDEtaDPhiVsQinvFakeFake{"cfgFillDEtaDPhiVsQinvFakeFake", true, "fill hDEtaDPhiVsQinv for FakeFake pairs"};
    Configurable<bool> cfgFillDEtaDPhiVsQinvPi0Daughters{"cfgFillDEtaDPhiVsQinvPi0Daughters", false, "fill hDEtaDPhiVsQinv for Pi0Daughters pairs"};
    Configurable<bool> cfgFillDRDZQinvTrueTrueDistinct{"cfgFillDRDZQinvTrueTrueDistinct", true, "fill hSparseDeltaRDeltaZQinv for TrueTrueDistinct pairs"};
    Configurable<bool> cfgFillDRDZQinvTrueTrueSamePhoton{"cfgFillDRDZQinvTrueTrueSamePhoton", false, "fill hSparseDeltaRDeltaZQinv for TrueTrueSamePhoton pairs"};
    Configurable<bool> cfgFillDRDZQinvSharedMcLeg{"cfgFillDRDZQinvSharedMcLeg", false, "fill hSparseDeltaRDeltaZQinv for SharedMcLeg pairs"};
    Configurable<bool> cfgFillDRDZQinvTrueFake{"cfgFillDRDZQinvTrueFake", false, "fill hSparseDeltaRDeltaZQinv for TrueFake pairs"};
    Configurable<bool> cfgFillDRDZQinvFakeFake{"cfgFillDRDZQinvFakeFake", true, "fill hSparseDeltaRDeltaZQinv for FakeFake pairs"};
    Configurable<bool> cfgFillDRDZQinvPi0Daughters{"cfgFillDRDZQinvPi0Daughters", false, "fill hSparseDeltaRDeltaZQinv for Pi0Daughters pairs"};
  } mctruthSparse;

  struct : ConfigurableGroup {
    std::string prefix = "ggpaircut_group";
    Configurable<float> cfgMinDRCosOA{"cfgMinDRCosOA", -1.f, "min. dr/cosOA; <0 = disabled"};
    Configurable<bool> cfgDoRCut{"cfgDoRCut", false, "apply |R1-R2| > cfgMinDeltaR cut"};
    Configurable<float> cfgMinDeltaR{"cfgMinDeltaR", 0.f, "minimum |R1-R2| (cm)"};
    Configurable<bool> cfgDoZCut{"cfgDoZCut", false, "apply |DeltaZ| > cfgMinDeltaZ cut"};
    Configurable<float> cfgMinDeltaZ{"cfgMinDeltaZ", 0.f, "minimum |DeltaZ| (cm)"};
    Configurable<bool> cfgDoEllipseCut{"cfgDoEllipseCut", false, "reject pairs inside ellipse in DeltaEta-DeltaPhi"};
    Configurable<float> cfgEllipseSigEta{"cfgEllipseSigEta", 0.02f, "sigma_eta for ellipse cut"};
    Configurable<float> cfgEllipseSigPhi{"cfgEllipseSigPhi", 0.02f, "sigma_phi for ellipse cut"};
    Configurable<float> cfgEllipseR2{"cfgEllipseR2", 1.0f, "R^2 threshold: reject if ellipse value < R^2"};
    Configurable<float> cfgMaxAsymmetry{"cfgMaxAsymmetry", -1.f, "max |p_{T, 1} - p_{T, 2}|/(p_{T, 1} + p_{T, 2}) asymmetry cut"};
  } ggpaircuts;

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require no TF border"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS ROF border"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx FT0 vs PV"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2.f, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000.f, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "no coll in time range std"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "no coll in time range strict"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "no coll in ITS ROF std"};
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "no coll in ITS ROF strict"};
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "no HM coll in prev ROF"};
    Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "ITS layer 3 OK"};
    Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "ITS layers 0-3 OK"};
    Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "all ITS layers OK"};
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfgRequireV0WithITSTPC{"cfgRequireV0WithITSTPC", false, "select V0s with ITS-TPC tracks"};
    Configurable<bool> cfgRequireV0WithITSOnly{"cfgRequireV0WithITSOnly", false, "select V0s with ITS-only tracks"};
    Configurable<bool> cfgRequireV0WithTPCOnly{"cfgRequireV0WithTPCOnly", false, "select V0s with TPC-only tracks"};
    Configurable<float> cfgMinPtV0{"cfgMinPtV0", 0.1, "min pT for V0 photons at PV"};
    Configurable<float> cfgMaxEtaV0{"cfgMaxEtaV0", 0.8, "max eta for V0 photons at PV"};
    Configurable<float> cfgMinV0Radius{"cfgMinV0Radius", 16.0, "min V0 radius"};
    Configurable<float> cfgMaxV0Radius{"cfgMaxV0Radius", 90.0, "max V0 radius"};
    Configurable<float> cfgMaxAlphaAP{"cfgMaxAlphaAP", 0.95, "max alpha for AP cut"};
    Configurable<float> cfgMaxQtAP{"cfgMaxQtAP", 0.01, "max qT for AP cut"};
    Configurable<float> cfgMinCosPA{"cfgMinCosPA", 0.997, "min V0 CosPA"};
    Configurable<float> cfgMaxPCA{"cfgMaxPCA", 3.0, "max distance between 2 legs"};
    Configurable<float> cfgMaxChi2KF{"cfgMaxChi2KF", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfgRejectV0OnITSIB{"cfgRejectV0OnITSIB", true, "reject V0s on ITSib"};
    Configurable<bool> cfgDisableITSOnlyTrack{"cfgDisableITSOnlyTrack", false, "disable ITS-only tracks"};
    Configurable<bool> cfgDisableTPCOnlyTrack{"cfgDisableTPCOnlyTrack", false, "disable TPC-only tracks"};
    Configurable<int> cfgMinNClusterTPC{"cfgMinNClusterTPC", 70, "min ncluster TPC"};
    Configurable<int> cfgMinNCrossedRows{"cfgMinNCrossedRows", 70, "min crossed rows"};
    Configurable<float> cfgMaxFracSharedClustersTPC{"cfgMaxFracSharedClustersTPC", 999.f, "max fraction of shared TPC clusters"};
    Configurable<float> cfgMaxChi2TPC{"cfgMaxChi2TPC", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfgMaxChi2ITS{"cfgMaxChi2ITS", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfgMinTPCNsigmaEl{"cfgMinTPCNsigmaEl", -3.5, "min TPC nsigma electron"};
    Configurable<float> cfgMaxTPCNsigmaEl{"cfgMaxTPCNsigmaEl", +3.5, "max TPC nsigma electron"};
  } pcmcuts;

  ~Photonhbt()
  {
    delete emh1;
    emh1 = nullptr;
    delete emh2;
    emh2 = nullptr;
    mapMixedEventIdToGlobalBC.clear();
    usedPhotonIdsPerCol.clear();
    truthGammaPool.clear();
  }

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry fRegistryPairQA{"outputPairQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry fRegistryPairMC{"outputPairMC", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry fRegistryMC{"outputMC", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::mt19937 engine;
  std::uniform_int_distribution<int> dist01;
  int mRunNumber{0};

  std::vector<float> ztxBinEdges;
  std::vector<float> centBinEdges;
  std::vector<float> epBinEgdes;
  std::vector<float> occBinEdges;

  inline bool isInsideEllipse(float deta, float dphi) const
  {
    if (!ggpaircuts.cfgDoEllipseCut.value)
      return false;
    const float sE = ggpaircuts.cfgEllipseSigEta.value;
    const float sP = ggpaircuts.cfgEllipseSigPhi.value;
    if (sE < kMinSigma || sP < kMinSigma)
      return false;
    return (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP) < ggpaircuts.cfgEllipseR2.value;
  }

  inline bool passRZCut(float deltaR, float deltaZ) const
  {
    if (ggpaircuts.cfgDoRCut.value && deltaR < ggpaircuts.cfgMinDeltaR.value)
      return false;
    if (ggpaircuts.cfgDoZCut.value && std::fabs(deltaZ) < ggpaircuts.cfgMinDeltaZ.value)
      return false;
    return true;
  }

  inline bool passAsymmetryCut(float pt1, float pt2) const
  {
    if (ggpaircuts.cfgMaxAsymmetry.value < 0.f) {
      return true;
    }

    const float sum = pt1 + pt2;
    if (sum < kMinSigma) {
      return false;
    }
    return std::fabs(pt1 - pt2) / sum < ggpaircuts.cfgMaxAsymmetry.value;
  }

  inline bool passQinvQAGate(float qinv) const
  {
    const float limit = qaflags.cfgMaxQinvForQA.value;
    return (limit <= 0.f) || (qinv < limit);
  }

  inline bool passQinvFullRangeGate(float qinv) const
  {
    const float limit = qaflags.cfgMaxQinvForFullRange.value;
    return (limit <= 0.f) || (qinv < limit);
  }

  inline bool passQinvMCQAGate(float qinv) const
  {
    const float limit = qaflags.cfgMaxQinvForMCQA.value;
    return (limit <= 0.f) || (qinv < limit);
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
    if (pairDir.R() < kMinSigma || p1cmDir.R() < kMinSigma)
      return -1.f;
    return static_cast<float>(pairDir.Unit().Dot(p1cmDir.Unit()));
  }

  static void parseBins(const ConfigurableAxis& cfg, std::vector<float>& edges)
  {
    if (cfg.value[0] == VARIABLE_WIDTH) {
      edges = std::vector<float>(cfg.value.begin(), cfg.value.end());
      edges.erase(edges.begin());
    } else {
      const int n = static_cast<int>(cfg.value[0]);
      const float xmin = static_cast<float>(cfg.value[1]);
      const float xmax = static_cast<float>(cfg.value[2]);
      edges.resize(n + 1);
      for (int i = 0; i <= n; ++i)
        edges[i] = xmin + (xmax - xmin) / n * i;
    }
  }

  static int clampBin(int b, int nmax) { return std::clamp(b, 0, nmax); }

  static int binOf(const std::vector<float>& edges, float val)
  {
    const int b = static_cast<int>(
                    std::lower_bound(edges.begin(), edges.end(), val) - edges.begin()) -
                  1;
    return clampBin(b, static_cast<int>(edges.size()) - 2); //
  }

  template <int ev_id, int step_id>
  static constexpr const char* qaPrefix()
  {
    if constexpr (ev_id == 0) {
      if constexpr (step_id == 0)
        return "Pair/same/QA/Before/";
      if constexpr (step_id == 1)
        return "Pair/same/QA/AfterDRCosOA/";
      if constexpr (step_id == 2) // o2-linter: disable=magic-number (just counting the step of a cut)
        return "Pair/same/QA/AfterRZ/";
      return "Pair/same/QA/AfterEllipse/";
    } else {
      if constexpr (step_id == 0)
        return "Pair/mix/QA/Before/";
      if constexpr (step_id == 1)
        return "Pair/mix/QA/AfterDRCosOA/";
      if constexpr (step_id == 2) // o2-linter: disable=magic-number (just counting the step of a cut)
        return "Pair/mix/QA/AfterRZ/";
      return "Pair/mix/QA/AfterEllipse/";
    }
  }

  template <int ev_id>
  static constexpr const char* fullRangePrefix()
  {
    if constexpr (ev_id == 0)
      return "Pair/same/FullRange/";
    return "Pair/mix/FullRange/";
  }

  void init(InitContext& /*context*/)
  {
    mRunNumber = 0;
    parseBins(mixing.confVtxBins, ztxBinEdges);
    parseBins(mixing.confCentBins, centBinEdges);
    parseBins(mixing.confEPBinsBins, epBinEgdes);
    parseBins(mixing.confOccupancyBins, occBinEdges);
    emh1 = new MyEMH(mixing.ndepth);
    emh2 = new MyEMH(mixing.ndepth);
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    DefineEMEventCut();
    DefinePCMCut();
    addhistograms();
    std::random_device seedGen;
    engine = std::mt19937(seedGen());
    dist01 = std::uniform_int_distribution<int>(0, 1);
    fRegistry.add("Pair/mix/hDiffBC",
                  "diff. global BC in mixed event;|BC_{current}-BC_{mixed}|",
                  kTH1D, {{10001, -0.5, 10000.5}}, true);

    // Print histogram counts and memory estimates for all registries
    // LOGF(info, "=== photonhbt histogram summary ===");
    // fRegistry.print();
    // fRegistryPairQA.print();
    // fRegistryPairMC.print();
    // fRegistryMC.print();
    // LOGF(info, "===================================");
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber())
      return;
    mRunNumber = collision.runNumber();
  }

  struct PairQAObservables {
    ROOT::Math::PtEtaPhiMVector v1, v2, k12;
    float x1 = 0.f, y1 = 0.f, z1 = 0.f, x2 = 0.f, y2 = 0.f, z2 = 0.f;
    float r1 = 0.f, r2 = 0.f, dx = 0.f, dy = 0.f, dz = 0.f;
    float deltaR = 0.f, deltaZ = 0.f, deltaRxy = 0.f, deltaR3D = 0.f;
    float opa = 0.f, cosOA = 0.f, drOverCosOA = 0.f;
    float deta = 0.f, dphi = 0.f, pairEta = 0.f, pairPhi = 0.f;
    float kt = 0.f, qinv = 0.f, cosTheta = 0.f, openingAngle = 0.f;
    bool valid = true;
  };

  struct TruthGamma {
    int id = -1, posId = -1, negId = -1;
    float eta = 0.f, phi = 0.f, pt = 0.f;
    float rTrue = -1.f;
    float legDRtrue = -1.f;
    float legDEta = 0.f; // ← neu
    float legDPhi = 0.f; // ← neu
    float alphaTrue = 0.f;
  };

  std::map<std::tuple<int, int, int, int>, std::deque<std::vector<TruthGamma>>> truthGammaPool;

  void addSinglePhotonQAHistogramsForStep(const std::string& path)
  {
    fRegistryPairQA.add((path + "hEtaVsPhiPt").c_str(), "acceptance;#phi (rad);#eta", kTH3D, {axisPhi, axisEta, axisPt}, true);
    fRegistryPairQA.add((path + "hRVsZConvPt").c_str(), "R_{conv} vs z_{conv};z_{conv} (cm);R_{conv} (cm)", kTH3D, {axisZConv, axisR, axisPt}, true);
  }

  void addFullRangeHistograms(const std::string& path)
  {
    fRegistry.add((path + "hDeltaRVsQinv").c_str(), "|R_{1}-R_{2}| vs q_{inv};q_{inv} (GeV/c);|R_{1}-R_{2}| (cm)", kTH2D, {axisQinv, axisDeltaR}, true);
    fRegistry.add((path + "hDeltaZVsQinv").c_str(), "#Delta z vs q_{inv};q_{inv} (GeV/c);#Delta z (cm)", kTH2D, {axisQinv, axisDeltaZ}, true);
    fRegistry.add((path + "hDeltaR3DVsQinv").c_str(), "#Delta r_{3D} vs q_{inv};q_{inv} (GeV/c);#Delta r_{3D} (cm)", kTH2D, {axisQinv, axisDeltaR3D}, true);
    fRegistry.add((path + "hQinvVsCent").c_str(), "q_{inv} vs centrality;centrality (%);q_{inv} (GeV/c)", kTH2D, {axisCentQA, axisQinv}, true);
    fRegistry.add((path + "hQinvVsOccupancy").c_str(), "q_{inv} vs occupancy;occupancy;q_{inv} (GeV/c)", kTH2D, {axisOccupancy, axisQinv}, true);
    fRegistry.add((path + "hSparseDeltaRDeltaZQinv").c_str(), "|R_{1}-R_{2}|,#Delta z,q_{inv}", kTHnSparseD, {axisDeltaR, axisDeltaZ, axisQinv}, true);
    fRegistry.add((path + "hDeltaRCosOAVsQinv").c_str(), "#Delta r/cos(#theta_{op}/2) vs q_{inv};q_{inv} (GeV/c);#Delta r/cos(#theta_{op}/2) (cm)", kTH2D, {axisQinv, {100, 0, 100}}, true);
  }

  template <int ev_id>
  inline void fillFullRangeQA(PairQAObservables const& obs, float cent, float occupancy)
  {
    constexpr auto base = fullRangePrefix<ev_id>();
    fRegistry.fill(HIST(base) + HIST("hDeltaRVsQinv"), obs.qinv, obs.deltaR);
    fRegistry.fill(HIST(base) + HIST("hDeltaZVsQinv"), obs.qinv, obs.deltaZ);
    fRegistry.fill(HIST(base) + HIST("hDeltaR3DVsQinv"), obs.qinv, obs.deltaR3D);
    fRegistry.fill(HIST(base) + HIST("hQinvVsCent"), cent, obs.qinv);
    fRegistry.fill(HIST(base) + HIST("hQinvVsOccupancy"), occupancy, obs.qinv);
    fRegistry.fill(HIST(base) + HIST("hSparseDeltaRDeltaZQinv"), obs.deltaR, obs.deltaZ, obs.qinv);
  }

  template <int ev_id>
  inline void fillFullRangeDeltaRCosOA(float qinv, float drOverCosOA)
  {
    constexpr auto base = fullRangePrefix<ev_id>();
    fRegistry.fill(HIST(base) + HIST("hDeltaRCosOAVsQinv"), qinv, drOverCosOA);
  }
  void addQAHistogramsForStep(const std::string& path)
  {
    //  Ellipse
    fRegistryPairQA.add((path + "hEllipseVal").c_str(), "(#Delta#eta/#sigma)^{2}+(#Delta#phi/#sigma)^{2};value;counts", kTH1D, {axisEllipseVal}, true);

    //  Conversion point
    fRegistryPairQA.add((path + "hR1VsR2").c_str(), "R_{1} vs R_{2};R_{1} (cm);R_{2} (cm)", kTH2D, {axisR, axisR}, true);
    fRegistryPairQA.add((path + "hDeltaRxyKt").c_str(), "#Delta r_{xy} vs k_{T};#Delta r_{xy} (cm);k_{T} (GeV/c)", kTH2D, {axisDeltaRxy, axisKt}, true);
    fRegistryPairQA.add((path + "hDeltaR3DKt").c_str(), "#Delta r_{3D} vs k_{T};#Delta r_{3D} (cm);k_{T} (GeV/c)", kTH2D, {axisDeltaR3D, axisKt}, true);

    //  Delta Eta QA
    fRegistryPairQA.add((path + "hDeltaEtaDeltaRKt").c_str(), "#Delta#eta,|R_{1}-R_{2}|,k_{T}", kTHnSparseD, {axisDeltaEta, axisDeltaR, axisKt}, true);
    fRegistryPairQA.add((path + "hDeltaEtaDeltaZKt").c_str(), "#Delta#eta,#Delta z,k_{T}", kTHnSparseD, {axisDeltaEta, axisDeltaZ, axisKt}, true);
    fRegistryPairQA.add((path + "hDeltaEtaEtaKt").c_str(), "#Delta#eta,#eta_{pair},k_{T}", kTHnSparseD, {axisDeltaEta, axisEta, axisKt}, true);
    fRegistryPairQA.add((path + "hDeltaEtaPhiKt").c_str(), "#Delta#eta,#phi_{pair},k_{T}", kTHnSparseD, {axisDeltaEta, axisPhi, axisKt}, true);

    // Delta Phi QA
    fRegistryPairQA.add((path + "hSparseDEtaDPhiKt").c_str(), "#Delta#eta,#Delta#phi,k_{T}", kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisKt}, true);
    fRegistryPairQA.add((path + "hDeltaPhiDeltaRKt").c_str(), "#Delta#phi,|R_{1}-R_{2}|,k_{T}", kTHnSparseD, {axisDeltaPhi, axisDeltaR, axisKt}, true);
    fRegistryPairQA.add((path + "hDeltaPhiDeltaZKt").c_str(), "#Delta#phi,#Delta z,k_{T}", kTHnSparseD, {axisDeltaPhi, axisDeltaZ, axisKt}, true);
    fRegistryPairQA.add((path + "hDeltaPhiPhiKt").c_str(), "#Delta#phi,#phi_{pair},k_{T}", kTHnSparseD, {axisDeltaPhi, axisPhi, axisKt}, true);
    fRegistryPairQA.add((path + "hDeltaPhiEtaKt").c_str(), "#Delta#phi,#eta_{pair},k_{T}", kTHnSparseD, {axisDeltaPhi, axisEta, axisKt}, true);

    // Delta Eta Delta Phi Diagnostics
    fRegistryPairQA.add((path + "hPhiVsEtaKt").c_str(), "#phi_{pair},#eta_{pair},k_{T}", kTHnSparseD, {axisPhi, axisEta, axisKt}, true);
    fRegistryPairQA.add((path + "hSparseDeltaRDeltaZKt").c_str(), "|R_{1}-R_{2}|,#Delta z,k_{T}", kTHnSparseD, {axisDeltaR, axisDeltaZ, axisKt}, true);
  }

  void addhistograms()
  {
    static constexpr std::string_view det[6] = {"FT0M", "FT0A", "FT0C", "BTot", "BPos", "BNeg"};
    fRegistry.add("Event/before/hEP2_CentFT0C_forMix",
                  Form("2nd harmonics EP for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", det[mixing.cfgEP2EstimatorForMix].data()),
                  kTH2D, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);
    fRegistry.add("Event/after/hEP2_CentFT0C_forMix",
                  Form("2nd harmonics EP for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", det[mixing.cfgEP2EstimatorForMix].data()),
                  kTH2D, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);

    addSinglePhotonQAHistogramsForStep("SinglePhoton/Before/");
    addSinglePhotonQAHistogramsForStep("SinglePhoton/AfterDRCosOA/");
    addSinglePhotonQAHistogramsForStep("SinglePhoton/AfterRZ/");
    addSinglePhotonQAHistogramsForStep("SinglePhoton/AfterEllipse/");

    if (hbtanalysis.cfgDo3D) {
      fRegistry.add("Pair/same/CF_3D", "diphoton correlation 3D LCMS", kTHnSparseD, {axisQout, axisQside, axisQlong, axisKt}, true);
      if (hbtanalysis.cfgDo2D)
        fRegistry.add("Pair/same/CF_2D", "diphoton correlation 2D (qout,qinv)", kTHnSparseD, {axisQout, axisQinv, axisKt}, true);
    } else {
      if (hbtanalysis.cfgUseLCMS)
        fRegistry.add("Pair/same/CF_1D", "diphoton correlation 1D LCMS", kTH2D, {axisQabsLcms, axisKt}, true);
      else
        fRegistry.add("Pair/same/CF_1D", "diphoton correlation 1D (qinv)", kTH2D, {axisQinv, axisKt}, true);
    }

    fRegistry.add("Pair/same/hDeltaRCosOA", "distance between 2 conversion points / cos(#theta_{op}/2);#Delta r / cos(#theta_{op}/2) (cm);counts", kTH1D, {{100, 0, 100}}, true);
    fRegistry.add("Pair/same/hSparse_DEtaDPhi_kT",
                  "same-event (#Delta#eta,#Delta#phi,q_{inv},k_{T}) for efficiency reweighting;"
                  "#Delta#eta;#Delta#phi (rad);q_{inv} (GeV/c);k_{T} (GeV/c)",
                  kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisKt}, true);

    fRegistry.add("Pair/same/hPhi_lowerPtV0", "azimuthal angle of lower-p_{T} V0 in pair;#phi (rad);counts", kTH1D, {axisPhi}, true);
    fRegistry.add("Pair/same/hSparse_DEtaDPhi_qinv_kT", "azimuthal angle of lower-p_{T} V0 in pair;#phi (rad);counts", kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisQinv, axisKt}, true);

    addQAHistogramsForStep("Pair/same/QA/Before/");
    addQAHistogramsForStep("Pair/same/QA/AfterDRCosOA/");
    addQAHistogramsForStep("Pair/same/QA/AfterRZ/");
    addQAHistogramsForStep("Pair/same/QA/AfterEllipse/");

    addMCHistograms();
    fRegistryPairQA.addClone("Pair/same/QA/", "Pair/mix/QA/");
    fRegistryPairMC.addClone("Pair/same/MC/", "Pair/mix/MC/");
    addFullRangeHistograms("Pair/same/FullRange/");
    fRegistry.addClone("Pair/same/", "Pair/mix/");

    fRegistry.add("Pair/same/EtaTopology/hSparse_DEtaDPhi_kT_sameSideV0_sameSideLegs",
                  "both V0 same #eta-side, all legs same side;"
                  "#Delta#eta;#Delta#phi (rad);k_{T} (GeV/c)",
                  kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisKt}, true);
    fRegistry.add("Pair/same/EtaTopology/hSparse_DEtaDPhi_kT_sameSideV0_mixedLegs",
                  "both V0 same #eta-side, legs mixed sides;"
                  "#Delta#eta;#Delta#phi (rad);k_{T} (GeV/c)",
                  kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisKt}, true);
    fRegistry.add("Pair/same/EtaTopology/hSparse_DEtaDPhi_kT_diffSideV0_matchingLegs",
                  "V0 on opposite #eta-sides, legs match their V0 side;"
                  "#Delta#eta;#Delta#phi (rad);k_{T} (GeV/c)",
                  kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisKt}, true);
    fRegistry.add("Pair/same/EtaTopology/hSparse_DEtaDPhi_kT_diffSideV0_crossedLegs",
                  "V0 on opposite #eta-sides, legs cross their V0 side;"
                  "#Delta#eta;#Delta#phi (rad);k_{T} (GeV/c)",
                  kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisKt}, true);
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
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfgMinCosPA);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfgMaxPCA);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.cfgMaxChi2KF);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfgMinV0Radius, pcmcuts.cfgMaxV0Radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfgMaxAlphaAP, pcmcuts.cfgMaxQtAP);
    fV0PhotonCut.RejectITSib(pcmcuts.cfgRejectV0OnITSIB);
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfgMinNClusterTPC);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfgMinNCrossedRows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfgMaxFracSharedClustersTPC);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfgMaxChi2TPC);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfgMinTPCNsigmaEl, pcmcuts.cfgMaxTPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfgMaxChi2ITS);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfgDisableITSOnlyTrack);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.cfgDisableTPCOnlyTrack);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetRequireITSTPC(pcmcuts.cfgRequireV0WithITSTPC);
    fV0PhotonCut.SetRequireITSonly(pcmcuts.cfgRequireV0WithITSOnly);
    fV0PhotonCut.SetRequireTPConly(pcmcuts.cfgRequireV0WithTPCOnly);
  }

  template <int step_id>
  static constexpr const char* singlePhotonQAPrefix()
  {
    if constexpr (step_id == 0)
      return "SinglePhoton/Before/";
    if constexpr (step_id == 1)
      return "SinglePhoton/AfterDRCosOA/";
    if constexpr (step_id == 2) // o2-linter: disable=magic-number (just counting the step of a cut)
      return "SinglePhoton/AfterRZ/";
    return "SinglePhoton/AfterEllipse/";
  }

  template <int step_id, typename TPhoton>
  inline void fillSinglePhotonQAStep(TPhoton const& g)
  {
    if (!qaflags.doSinglePhotonQa)
      return;
    constexpr auto base = singlePhotonQAPrefix<step_id>();
    const float r = std::sqrt(g.vx() * g.vx() + g.vy() * g.vy());
    fRegistryPairQA.fill(HIST(base) + HIST("hEtaVsPhiPt"), g.phi(), g.eta(), g.pt());
    fRegistryPairQA.fill(HIST(base) + HIST("hRVsZConvPt"), g.vz(), r, g.pt());
  }

  template <int ev_id, typename TCollision>
  void fillPairHistogram(TCollision const& /*collision*/,
                         ROOT::Math::PtEtaPhiMVector v1,
                         ROOT::Math::PtEtaPhiMVector v2,
                         float weight = 1.f)
  {
    float rndm = std::pow(-1, dist01(engine) % 2);
    auto k12 = 0.5 * (v1 + v2);
    float kt = k12.Pt();
    float qinv = -(((v1 - v2) * rndm).M());
    ROOT::Math::XYZVector uv_out(k12.Px() / k12.Pt(), k12.Py() / k12.Pt(), 0);
    ROOT::Math::XYZVector uv_long(0, 0, 1);
    ROOT::Math::XYZVector uv_side = uv_out.Cross(uv_long);
    ROOT::Math::PxPyPzEVector v1c(v1), v2c(v2);
    float beta_z = (v1 + v2).Beta() * std::cos((v1 + v2).Theta());
    ROOT::Math::Boost bst_z(0, 0, -beta_z);
    auto q12_lcms = bst_z((v1c - v2c) * rndm);
    auto q3_lcms = q12_lcms.Vect();
    float qabs_lcms = q3_lcms.R();
    float qout_lcms = q3_lcms.Dot(uv_out);
    float qside_lcms = q3_lcms.Dot(uv_side);
    float qlong_lcms = q3_lcms.Dot(uv_long);
    if (hbtanalysis.cfgDo3D) {
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("CF_3D"),
                     std::fabs(qout_lcms), std::fabs(qside_lcms), std::fabs(qlong_lcms), kt, weight);
      if (hbtanalysis.cfgDo2D)
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("CF_2D"),
                       std::fabs(qout_lcms), std::fabs(qinv), kt, weight);
    } else {
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("CF_1D"),
                     hbtanalysis.cfgUseLCMS ? qabs_lcms : qinv, kt, weight);
    }
    float deta_pair = v1.Eta() - v2.Eta();
    float dphi_pair = v1.Phi() - v2.Phi();
    dphi_pair = RecoDecay::constrainAngle(dphi_pair, -o2::constants::math::PI);
    if constexpr (ev_id == 0) {
      fRegistry.fill(HIST("Pair/same/hSparse_DEtaDPhi_qinv_kT"), deta_pair, dphi_pair, qinv, kt, weight);
    } else {
      fRegistry.fill(HIST("Pair/mix/hSparse_DEtaDPhi_qinv_kT"), deta_pair, dphi_pair, qinv, kt, weight);
    }
  }

  template <int ev_id, PairTruthType TruthT, typename TCollision>
  void fillPairHistogramMC(TCollision const& /*collision*/,
                           ROOT::Math::PtEtaPhiMVector v1,
                           ROOT::Math::PtEtaPhiMVector v2,
                           float weight = 1.f)
  {
    float rndm = std::pow(-1, dist01(engine) % 2);
    auto k12 = 0.5 * (v1 + v2);
    float kt = k12.Pt();
    float qinv = -(((v1 - v2) * rndm).M());
    ROOT::Math::XYZVector uv_out(k12.Px() / k12.Pt(), k12.Py() / k12.Pt(), 0);
    ROOT::Math::XYZVector uv_long(0, 0, 1);
    ROOT::Math::XYZVector uv_side = uv_out.Cross(uv_long);
    ROOT::Math::PxPyPzEVector v1c(v1), v2c(v2);
    float beta_z = (v1 + v2).Beta() * std::cos((v1 + v2).Theta());
    ROOT::Math::Boost bst_z(0, 0, -beta_z);
    auto q12_lcms = bst_z((v1c - v2c) * rndm);
    auto q3_lcms = q12_lcms.Vect();
    float qabs_lcms = q3_lcms.R();
    float qout_lcms = q3_lcms.Dot(uv_out);
    float qside_lcms = q3_lcms.Dot(uv_side);
    float qlong_lcms = q3_lcms.Dot(uv_long);
    constexpr auto mcDir = []() constexpr -> const char* {
      if constexpr (ev_id == 0) {
        if constexpr (TruthT == PairTruthType::TrueTrueDistinct)
          return "Pair/same/MC/TrueTrueDistinct/";
        if constexpr (TruthT == PairTruthType::TrueTrueSamePhoton)
          return "Pair/same/MC/TrueTrueSamePhoton/";
        if constexpr (TruthT == PairTruthType::SharedMcLeg)
          return "Pair/same/MC/SharedMcLeg/";
        if constexpr (TruthT == PairTruthType::TrueFake)
          return "Pair/same/MC/TrueFake/";
        if constexpr (TruthT == PairTruthType::FakeFake)
          return "Pair/same/MC/FakeFake/";
        return "Pair/same/MC/Pi0Daughters/";
      } else {
        if constexpr (TruthT == PairTruthType::TrueTrueDistinct)
          return "Pair/mix/MC/TrueTrueDistinct/";
        if constexpr (TruthT == PairTruthType::TrueTrueSamePhoton)
          return "Pair/mix/MC/TrueTrueSamePhoton/";
        if constexpr (TruthT == PairTruthType::SharedMcLeg)
          return "Pair/mix/MC/SharedMcLeg/";
        if constexpr (TruthT == PairTruthType::TrueFake)
          return "Pair/mix/MC/TrueFake/";
        if constexpr (TruthT == PairTruthType::FakeFake)
          return "Pair/mix/MC/FakeFake/";
        return "Pair/mix/MC/Pi0Daughters/";
      }
    }();
    if (hbtanalysis.cfgDo3D) {
      fRegistryPairMC.fill(HIST(mcDir) + HIST("CF_3D"),
                           std::fabs(qout_lcms), std::fabs(qside_lcms), std::fabs(qlong_lcms), kt, weight);
      if (hbtanalysis.cfgDo2D)
        fRegistryPairMC.fill(HIST(mcDir) + HIST("CF_2D"), std::fabs(qout_lcms), std::fabs(qinv), kt, weight);
    } else {
      fRegistryPairMC.fill(HIST(mcDir) + HIST("CF_1D"), hbtanalysis.cfgUseLCMS ? qabs_lcms : qinv, kt, weight);
    }
    float deta_pair = v1.Eta() - v2.Eta();
    float dphi_pair = v1.Phi() - v2.Phi();
    dphi_pair = RecoDecay::constrainAngle(dphi_pair, -o2::constants::math::PI);
    if constexpr (ev_id == 0) {
      fRegistry.fill(HIST("Pair/same/hSparse_DEtaDPhi_qinv_kT"), deta_pair, dphi_pair, qinv, kt, weight);
    } else {
      fRegistry.fill(HIST("Pair/mix/hSparse_DEtaDPhi_qinv_kT"), deta_pair, dphi_pair, qinv, kt, weight);
    }
  }

  template <typename TG1, typename TG2>
  PairQAObservables buildPairQAObservables(TG1 const& g1, TG2 const& g2)
  {
    PairQAObservables o{};
    o.x1 = g1.vx();
    o.y1 = g1.vy();
    o.z1 = g1.vz();
    o.x2 = g2.vx();
    o.y2 = g2.vy();
    o.z2 = g2.vz();
    o.r1 = std::sqrt(o.x1 * o.x1 + o.y1 * o.y1);
    o.r2 = std::sqrt(o.x2 * o.x2 + o.y2 * o.y2);
    o.dx = o.x1 - o.x2;
    o.dy = o.y1 - o.y2;
    o.dz = o.z1 - o.z2;
    o.deltaR = std::fabs(o.r1 - o.r2);
    o.deltaZ = o.dz;
    o.deltaRxy = std::sqrt(o.dx * o.dx + o.dy * o.dy);
    o.deltaR3D = std::sqrt(o.dx * o.dx + o.dy * o.dy + o.dz * o.dz);
    ROOT::Math::XYZVector cp1(o.x1, o.y1, o.z1), cp2(o.x2, o.y2, o.z2);
    const float mag1 = std::sqrt(cp1.Mag2()), mag2 = std::sqrt(cp2.Mag2());
    if (mag1 < kMinMagnitude || mag2 < kMinMagnitude) {
      o.valid = false;
      return o;
    }
    float cosPA = static_cast<float>(cp1.Dot(cp2) / (mag1 * mag2));
    cosPA = std::clamp(cosPA, -1.f, 1.f);
    o.opa = std::acos(cosPA);
    o2::math_utils::bringTo02Pi(o.opa);
    if (o.opa > o2::constants::math::PI)
      o.opa -= o2::constants::math::PI;
    o.cosOA = std::cos(o.opa / 2.f);
    o.drOverCosOA = (std::fabs(o.cosOA) < kMinCosine) ? 1e12f : (o.deltaR3D / o.cosOA);
    o.v1 = ROOT::Math::PtEtaPhiMVector(g1.pt(), g1.eta(), g1.phi(), 0.f);
    o.v2 = ROOT::Math::PtEtaPhiMVector(g2.pt(), g2.eta(), g2.phi(), 0.f);
    o.k12 = 0.5f * (o.v1 + o.v2);
    o.deta = g1.eta() - g2.eta();
    o.dphi = RecoDecay::constrainAngle(g1.phi() - g2.phi(), -o2::constants::math::PI);
    o.pairEta = 0.5f * (g1.eta() + g2.eta());
    o.pairPhi = RecoDecay::constrainAngle(o.k12.Phi(), 0.f);
    o.kt = o.k12.Pt();
    o.qinv = std::fabs((o.v1 - o.v2).M());
    o.cosTheta = std::fabs(computeCosTheta(o.v1, o.v2));
    o.openingAngle = o.opa;
    return o;
  }

  template <int ev_id, int step_id>
  inline void fillPairQAStep(PairQAObservables const& o, float /*cent*/, float /*occupancy*/)
  {
    if (!qaflags.doPairQa)
      return;
    constexpr auto base = qaPrefix<ev_id, step_id>();

    ///// Delta Eta QA

    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaEtaDeltaRKt"), o.deta, o.deltaR, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaEtaEtaKt"), o.deta, o.pairEta, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaEtaPhiKt"), o.deta, o.pairPhi, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaEtaDeltaZKt"), o.deta, o.deltaZ, o.kt);

    ///// Delta Phi QA

    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaPhiDeltaRKt"), o.dphi, o.deltaR, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaPhiDeltaZKt"), o.dphi, o.deltaZ, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hSparseDEtaDPhiKt"), o.deta, o.dphi, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaPhiPhiKt"), o.dphi, o.pairPhi, o.kt);

    // Delta Eta Dleta Phi Stuff

    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaPhiEtaKt"), o.dphi, o.pairEta, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hPhiVsEtaKt"), o.pairPhi, o.pairEta, o.kt);

    //// Delta R (Conversion point) QA

    fRegistryPairQA.fill(HIST(base) + HIST("hR1VsR2"), o.r1, o.r2);
    fRegistryPairQA.fill(HIST(base) + HIST("hSparseDeltaRDeltaZKt"), o.deltaR, o.deltaZ, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaRxyKt"), o.deltaRxy, o.kt);
    fRegistryPairQA.fill(HIST(base) + HIST("hDeltaR3DKt"), o.deltaR3D, o.kt);

    const float sE = ggpaircuts.cfgEllipseSigEta.value, sP = ggpaircuts.cfgEllipseSigPhi.value;
    if (sE > kMinSigma && sP > kMinSigma)
      fRegistryPairQA.fill(HIST(base) + HIST("hEllipseVal"), (o.deta / sE) * (o.deta / sE) + (o.dphi / sP) * (o.dphi / sP));
  }

  template <typename TPhoton, typename TLegs, typename TMCParticles>
  static PhotonMCInfo buildPhotonMCInfo(TPhoton const& g, TMCParticles const& mcParticles)
  {
    PhotonMCInfo info{};
    const auto pos = g.template posTrack_as<TLegs>();
    const auto neg = g.template negTrack_as<TLegs>();
    if (!pos.has_emmcparticle() || !neg.has_emmcparticle())
      return info;
    info.hasMC = true;
    info.mcPosId = pos.emmcparticleId();
    info.mcNegId = neg.emmcparticleId();
    const auto mcPos = pos.template emmcparticle_as<TMCParticles>();
    const auto mcNeg = neg.template emmcparticle_as<TMCParticles>();
    if (!mcPos.has_mothers() || !mcNeg.has_mothers())
      return info;
    const int mothIdPos = mcPos.mothersIds()[0], mothIdNeg = mcNeg.mothersIds()[0];
    if (mothIdPos != mothIdNeg)
      return info;
    info.sameMother = true;
    info.motherId = mothIdPos;
    const auto mother = mcParticles.iteratorAt(mothIdPos);
    info.motherPdg = mother.pdgCode();
    info.isTruePhoton = (info.motherPdg == kGamma);
    info.isPhysicalPrimary = mother.isPhysicalPrimary();
    return info;
  }

  static PairTruthType classifyPairTruth(PhotonMCInfo const& m1, PhotonMCInfo const& m2)
  {
    const bool t1 = m1.hasMC && m1.sameMother && m1.isTruePhoton;
    const bool t2 = m2.hasMC && m2.sameMother && m2.isTruePhoton;
    if (m1.hasMC && m2.hasMC) {
      if ((m1.mcPosId >= 0 && (m1.mcPosId == m2.mcPosId || m1.mcPosId == m2.mcNegId)) ||
          (m1.mcNegId >= 0 && (m1.mcNegId == m2.mcPosId || m1.mcNegId == m2.mcNegId)))
        return PairTruthType::SharedMcLeg;
    }
    if (!t1 && !t2)
      return PairTruthType::FakeFake;
    if (t1 != t2)
      return PairTruthType::TrueFake;
    if (m1.motherId >= 0 && m1.motherId == m2.motherId)
      return PairTruthType::TrueTrueSamePhoton;
    return PairTruthType::TrueTrueDistinct;
  }

  enum class EtaTopology : uint8_t {
    SameSideV0SameSideLegs = 0, ///< both V0 eta same sign; all 4 legs same sign
    SameSideV0MixedLegs = 1,    ///< both V0 eta same sign; legs not all on that sign
    DiffSideV0MatchingLegs = 2, ///< V0s on opposite sides; each V0's legs match its own side
    DiffSideV0CrossedLegs = 3,  ///< V0s on opposite sides; legs do NOT match their V0's side
  };

  /// Classify the eta-side topology of a photon pair from the 6 raw eta values.
  static EtaTopology classifyEtaTopology(float v1eta, float v2eta,
                                         float pos1eta, float neg1eta,
                                         float pos2eta, float neg2eta)
  {
    const bool v1pos = v1eta >= 0.f;
    const bool v2pos = v2eta >= 0.f;
    if (v1pos == v2pos) { // same-side V0s
      const bool allSame =
        ((pos1eta >= 0.f) == v1pos) && ((neg1eta >= 0.f) == v1pos) &&
        ((pos2eta >= 0.f) == v1pos) && ((neg2eta >= 0.f) == v1pos);
      return allSame ? EtaTopology::SameSideV0SameSideLegs
                     : EtaTopology::SameSideV0MixedLegs;
    }
    // different-side V0s
    const bool v1match = ((pos1eta >= 0.f) == v1pos) && ((neg1eta >= 0.f) == v1pos);
    const bool v2match = ((pos2eta >= 0.f) == v2pos) && ((neg2eta >= 0.f) == v2pos);
    return (v1match && v2match) ? EtaTopology::DiffSideV0MatchingLegs
                                : EtaTopology::DiffSideV0CrossedLegs;
  }

  inline void fillEtaTopologyHisto(EtaTopology topo, float deta, float dphi, float kt)
  {
    switch (topo) {
      case EtaTopology::SameSideV0SameSideLegs:
        fRegistry.fill(HIST("Pair/same/EtaTopology/hSparse_DEtaDPhi_kT_sameSideV0_sameSideLegs"),
                       deta, dphi, kt);
        break;
      case EtaTopology::SameSideV0MixedLegs:
        fRegistry.fill(HIST("Pair/same/EtaTopology/hSparse_DEtaDPhi_kT_sameSideV0_mixedLegs"),
                       deta, dphi, kt);
        break;
      case EtaTopology::DiffSideV0MatchingLegs:
        fRegistry.fill(HIST("Pair/same/EtaTopology/hSparse_DEtaDPhi_kT_diffSideV0_matchingLegs"),
                       deta, dphi, kt);
        break;
      case EtaTopology::DiffSideV0CrossedLegs:
        fRegistry.fill(HIST("Pair/same/EtaTopology/hSparse_DEtaDPhi_kT_diffSideV0_crossedLegs"),
                       deta, dphi, kt);
        break;
    }
  }

  template <typename TMCParticles>
  static bool isPi0DaughterPair(PhotonMCInfo const& m1, PhotonMCInfo const& m2,
                                TMCParticles const& mcParticles)
  {
    if (!m1.isTruePhoton || !m2.isTruePhoton || m1.motherId < 0 || m2.motherId < 0)
      return false;
    const auto ph1 = mcParticles.iteratorAt(m1.motherId);
    const auto ph2 = mcParticles.iteratorAt(m2.motherId);
    if (!ph1.has_mothers() || !ph2.has_mothers())
      return false;
    const int gm1 = ph1.mothersIds()[0], gm2 = ph2.mothersIds()[0];
    if (gm1 != gm2)
      return false;
    return (std::abs(mcParticles.iteratorAt(gm1).pdgCode()) == kPi0);
  }

  static constexpr std::string_view pairTruthLabel(PairTruthType t)
  {
    switch (t) {
      case PairTruthType::TrueTrueDistinct:
        return "TrueTrueDistinct/";
      case PairTruthType::TrueTrueSamePhoton:
        return "TrueTrueSamePhoton/";
      case PairTruthType::SharedMcLeg:
        return "SharedMcLeg/";
      case PairTruthType::TrueFake:
        return "TrueFake/";
      case PairTruthType::FakeFake:
        return "FakeFake/";
      case PairTruthType::Pi0Daughters:
        return "Pi0Daughters/";
      default:
        return "Unknown/";
    }
  }
  void addMCHistograms()
  {
    const AxisSpec axisTruthType{{0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5},
                                 "truth type (1=TrueTrueDistinct,2=TrueTrueSamePhoton,3=SharedMcLeg,"
                                 "4=TrueFake,5=FakeFake,6=Pi0Daughters)"};
    const AxisSpec axisDeltaEtaMC{90, -1.6f, +1.6f, "#Delta#eta"};
    const AxisSpec axisDeltaPhiMC{90, -o2::constants::math::PI, +o2::constants::math::PI, "#Delta#phi (rad)"};
    const AxisSpec axQinvMC{60, 0.f, 0.3f, "q_{inv}^{true} (GeV/c)"};
    const AxisSpec axRconv{180, 0.f, 90.f, "R_{conv}^{true} (cm)"};
    const AxisSpec axAlpha{100, -1.f, 1.f, "#alpha^{true}"};
    const AxisSpec axLegDR{100, 0.f, 0.3f, "leg #Delta R^{true}"};

    // fRegistryPairMC  — reco-level pair histograms, per MC truth type

    //  Per-type CF + observables
    static constexpr std::array<std::string_view, 6> kTypes = {
      "TrueTrueDistinct/", "TrueTrueSamePhoton/", "SharedMcLeg/",
      "TrueFake/", "FakeFake/", "Pi0Daughters/"};

    for (const auto& label : kTypes) {
      const std::string base = std::string("Pair/same/MC/") + std::string(label);

      // CF
      if (hbtanalysis.cfgDo3D) {
        fRegistryPairMC.add((base + "CF_3D").c_str(), "MC CF 3D LCMS",
                            kTHnSparseD, {axisQout, axisQside, axisQlong, axisKt}, true);
        if (hbtanalysis.cfgDo2D)
          fRegistryPairMC.add((base + "CF_2D").c_str(), "MC CF 2D",
                              kTHnSparseD, {axisQout, axisQinv, axisKt}, true);
      } else {
        fRegistryPairMC.add((base + "CF_1D").c_str(),
                            hbtanalysis.cfgUseLCMS ? "MC CF 1D LCMS" : "MC CF 1D (qinv)",
                            kTH2D, {hbtanalysis.cfgUseLCMS ? axisQabsLcms : axisQinv, axisKt}, true);
      }

      // 1D observables
      fRegistryPairMC.add((base + "hQinv").c_str(), "q_{inv};q_{inv} (GeV/c);counts", kTH1D, {axisQinv}, true);
      fRegistryPairMC.add((base + "hDeltaR").c_str(), "|R_{1}-R_{2}|;|R_{1}-R_{2}| (cm);counts", kTH1D, {axisDeltaR}, true);
      fRegistryPairMC.add((base + "hDeltaZ").c_str(), "#Delta z;#Delta z (cm);counts", kTH1D, {axisDeltaZ}, true);
      fRegistryPairMC.add((base + "hDeltaR3D").c_str(), "#Delta r_{3D};#Delta r_{3D} (cm);counts", kTH1D, {axisDeltaR3D}, true);

      // 2D observables
      fRegistryPairMC.add((base + "hDEtaDPhi").c_str(), "#Delta#eta vs #Delta#phi", kTH2D, {axisDeltaEta, axisDeltaPhi}, true);
      fRegistryPairMC.add((base + "hDeltaRVsQinv").c_str(), "|R_{1}-R_{2}| vs q_{inv}", kTH2D, {axisQinv, axisDeltaR}, true);
      fRegistryPairMC.add((base + "hDeltaZVsQinv").c_str(), "#Delta z vs q_{inv}", kTH2D, {axisQinv, axisDeltaZ}, true);
      fRegistryPairMC.add((base + "hDeltaR3DVsQinv").c_str(), "#Delta r_{3D} vs q_{inv}", kTH2D, {axisQinv, axisDeltaR3D}, true);

      // Sparse (conditional)
      fRegistryPairMC.add((base + "hSparse_DEtaDPhi_kT").c_str(),
                          "#Delta#eta,#Delta#phi,k_{T};#Delta#eta;#Delta#phi (rad);k_{T} (GeV/c)",
                          kTHnSparseD, {axisDeltaEtaMC, axisDeltaPhiMC, axisKt}, true);

      const bool addDEtaDPhiVsQinv =
        (label == "TrueTrueDistinct/")     ? mctruthSparse.cfgFillDEtaDPhiVsQinvTrueTrueDistinct.value
        : (label == "TrueTrueSamePhoton/") ? mctruthSparse.cfgFillDEtaDPhiVsQinvTrueTrueSamePhoton.value
        : (label == "SharedMcLeg/")        ? mctruthSparse.cfgFillDEtaDPhiVsQinvSharedMcLeg.value
        : (label == "TrueFake/")           ? mctruthSparse.cfgFillDEtaDPhiVsQinvTrueFake.value
        : (label == "FakeFake/")           ? mctruthSparse.cfgFillDEtaDPhiVsQinvFakeFake.value
                                           : mctruthSparse.cfgFillDEtaDPhiVsQinvPi0Daughters.value;
      if (addDEtaDPhiVsQinv)
        fRegistryPairMC.add((base + "hDEtaDPhiVsQinv").c_str(),
                            "#Delta#eta vs #Delta#phi vs q_{inv}", kTHnSparseD,
                            {axisDeltaEtaMC, axisDeltaPhiMC, axisQinv}, true);

      const bool addDRDZQinv =
        (label == "TrueTrueDistinct/")     ? mctruthSparse.cfgFillDRDZQinvTrueTrueDistinct.value
        : (label == "TrueTrueSamePhoton/") ? mctruthSparse.cfgFillDRDZQinvTrueTrueSamePhoton.value
        : (label == "SharedMcLeg/")        ? mctruthSparse.cfgFillDRDZQinvSharedMcLeg.value
        : (label == "TrueFake/")           ? mctruthSparse.cfgFillDRDZQinvTrueFake.value
        : (label == "FakeFake/")           ? mctruthSparse.cfgFillDRDZQinvFakeFake.value
                                           : mctruthSparse.cfgFillDRDZQinvPi0Daughters.value;
      if (addDRDZQinv)
        fRegistryPairMC.add((base + "hSparseDeltaRDeltaZQinv").c_str(),
                            "|R_{1}-R_{2}|,#Delta z,q_{inv}", kTHnSparseD,
                            {axisDeltaR, axisDeltaZ, axisQinv}, true);
    }

    fRegistryPairMC.add("Pair/same/MC/hTruthTypeVsQinv",
                        "truth type vs q_{inv};q_{inv} (GeV/c);truth type", kTH2D, {axisQinv, axisTruthType}, true);
    fRegistryPairMC.add("Pair/same/MC/hTruthTypeVsKt",
                        "truth type vs k_{T};k_{T} (GeV/c);truth type", kTH2D, {axisKt, axisTruthType}, true);

    fRegistryPairMC.add("Pair/same/MC/hDEtaDPhi_truePairs",
                        "true reco pairs (TrueTrueDistinct+SamePhoton+Pi0);"
                        "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad)",
                        kTH2D, {axisDeltaEta, axisDeltaPhi}, true);
    fRegistryPairMC.add("Pair/same/MC/hDEtaDPhi_fakePairs",
                        "fake reco pairs (FakeFake+TrueFake+SharedMcLeg);"
                        "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad)",
                        kTH2D, {axisDeltaEta, axisDeltaPhi}, true);
    fRegistryPairMC.add("Pair/same/MC/hSparse_DEtaDPhi_kT_truePairs",
                        "true pairs: #Delta#eta,#Delta#phi,k_{T};"
                        "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad);k_{T} (GeV/c)",
                        kTHnSparseD, {axisDeltaEtaMC, axisDeltaPhiMC, axisKt}, true);
    fRegistryPairMC.add("Pair/same/MC/hSparse_DEtaDPhi_kT_fakePairs",
                        "fake pairs: #Delta#eta,#Delta#phi,k_{T};"
                        "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad);k_{T} (GeV/c)",
                        kTHnSparseD, {axisDeltaEtaMC, axisDeltaPhiMC, axisKt}, true);
    fRegistryPairMC.add("Pair/same/MC/hSparse_DEtaDPhi_qinv_truePairs",
                        "true pairs: #Delta#eta,#Delta#phi,q_{inv};"
                        "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad);q_{inv} (GeV/c)",
                        kTHnSparseD, {axisDeltaEtaMC, axisDeltaPhiMC, axisQinv}, true);
    fRegistryPairMC.add("Pair/same/MC/hSparse_DEtaDPhi_qinv_fakePairs",
                        "fake pairs: #Delta#eta,#Delta#phi,q_{inv};"
                        "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad);q_{inv} (GeV/c)",
                        kTHnSparseD, {axisDeltaEtaMC, axisDeltaPhiMC, axisQinv}, true);

    // ─── Pairs with missing MC label ─────────────────────────────────────────
    if (hbtanalysis.cfgDo3D) {
      fRegistryPairMC.add("Pair/same/MC/NoLabel/CF_3D",
                          "missing MC label — CF 3D LCMS", kTHnSparseD, {axisQout, axisQside, axisQlong, axisKt}, true);
      if (hbtanalysis.cfgDo2D)
        fRegistryPairMC.add("Pair/same/MC/NoLabel/CF_2D",
                            "missing MC label — CF 2D", kTHnSparseD, {axisQout, axisQinv, axisKt}, true);
    } else {
      fRegistryPairMC.add("Pair/same/MC/NoLabel/CF_1D",
                          hbtanalysis.cfgUseLCMS ? "missing MC label — CF 1D LCMS" : "missing MC label — CF 1D (qinv)",
                          kTH2D, {hbtanalysis.cfgUseLCMS ? axisQabsLcms : axisQinv, axisKt}, true);
    }
    fRegistryPairMC.add("Pair/same/MC/NoLabel/hDEtaDPhi",
                        "missing MC label: #Delta#eta vs #Delta#phi;"
                        "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad)",
                        kTH2D, {axisDeltaEta, axisDeltaPhi}, true);
    fRegistryPairMC.add("Pair/same/MC/NoLabel/hQinv",
                        "missing MC label: q_{inv};q_{inv} (GeV/c);counts", kTH1D, {axisQinv}, true);
    fRegistryPairMC.add("Pair/same/MC/NoLabel/hKt",
                        "missing MC label: k_{T};k_{T} (GeV/c);counts", kTH1D, {axisKt}, true);

    // fRegistryMC  — truth-level histograms

    // ─── Truth-level CF
    fRegistryMC.add("MC/TruthCF/hQinvVsKt_same", "truth-level same-event CF;k_{T} (GeV/c);q_{inv}^{true} (GeV/c)", kTH2D, {axisKt, axQinvMC}, true);
    fRegistryMC.add("MC/TruthCF/hDEtaDPhi_same",
                    "truth-level same-event #Delta#eta vs #Delta#phi;"
                    "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad)",
                    kTH2D, {axisDeltaEta, axisDeltaPhi}, true);
    fRegistryMC.add("MC/TruthCF/hQinvVsKt_mix", "truth-level mixed-event CF;k_{T} (GeV/c);q_{inv}^{true} (GeV/c)", kTH2D, {axisKt, axQinvMC}, true);
    fRegistryMC.add("MC/TruthCF/hDEtaDPhi_mix",
                    "truth-level mixed-event #Delta#eta vs #Delta#phi;"
                    "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad)",
                    kTH2D, {axisDeltaEta, axisDeltaPhi}, true);

    auto addStageHistos = [&](const char* suffix, const char* title, auto... axes) {
      for (const auto& stage : {"truthConverted", "all4LegsThisColl",
                                "bothPhotonsBuilt", "bothPhotonsSelected"}) {
        const std::string name = std::string("MC/TruthAO2D/") + suffix + std::string("_") + stage;
        const std::string ttl = std::string(title) + std::string(" [") + stage + "]";
        fRegistryMC.add(name.c_str(), ttl.c_str(), kTHnD, {axes...}, true);
      }
    };
    auto addStageHistos2D = [&](const char* suffix, const char* title, auto ax1, auto ax2) {
      for (const auto& stage : {"truthConverted", "all4LegsThisColl",
                                "bothPhotonsBuilt", "bothPhotonsSelected"}) {
        const std::string name = std::string("MC/TruthAO2D/") + suffix + std::string("_") + stage;
        fRegistryMC.add(name.c_str(), title, kTH2D, {ax1, ax2}, true);
      }
    };

    addStageHistos("hSparse_DEtaDPhi_qinv",
                   "#Delta#eta,#Delta#phi,q_{inv}^{true};"
                   "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad);q_{inv}^{true} (GeV/c)",
                   axisDeltaEta, axisDeltaPhi, axQinvMC);

    addStageHistos("hSparse_DEtaDPhi_kT",
                   "#Delta#eta,#Delta#phi,k_{T};"
                   "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad);k_{T} (GeV/c)",
                   axisDeltaEta, axisDeltaPhi, axisKt);

    addStageHistos2D("hQinvVsKt",
                     "q_{inv}^{true} vs k_{T};k_{T} (GeV/c);q_{inv}^{true} (GeV/c)",
                     axisKt, axQinvMC);

    addStageHistos2D("hDEtaDPhi",
                     "#Delta#eta vs #Delta#phi;#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad)",
                     axisDeltaEta, axisDeltaPhi);

    // Rconv waterfall (nur converted + selected)
    fRegistryMC.add("MC/TruthAO2D/hRconv1_vs_Rconv2_truthConverted",
                    "denominator: R_{conv,1} vs R_{conv,2};R_{conv,1}^{true} (cm);R_{conv,2}^{true} (cm)",
                    kTH2D, {axRconv, axRconv}, true);
    fRegistryMC.add("MC/TruthAO2D/hRconv1_vs_Rconv2_bothPhotonsSelected",
                    "numerator: R_{conv,1} vs R_{conv,2};R_{conv,1}^{true} (cm);R_{conv,2}^{true} (cm)",
                    kTH2D, {axRconv, axRconv}, true);
    fRegistryMC.add("MC/TruthAO2D/hMinRconv_vs_kT_truthConverted",
                    "denominator: min(R_{conv}) vs k_{T};k_{T} (GeV/c);min(R_{conv}^{true}) (cm)",
                    kTH2D, {axisKt, axRconv}, true);
    fRegistryMC.add("MC/TruthAO2D/hMinRconv_vs_kT_bothPhotonsSelected",
                    "numerator: min(R_{conv}) vs k_{T};k_{T} (GeV/c);min(R_{conv}^{true}) (cm)",
                    kTH2D, {axisKt, axRconv}, true);

    // Stage waterfall summary + consistency
    fRegistryMC.add("MC/TruthAO2D/hStage_vs_kT",
                    "efficiency waterfall;k_{T} (GeV/c);stage (0=converted,1=all4legs,2=bothBuilt,3=bothSel)",
                    kTH2D, {axisKt, AxisSpec{4, -0.5f, 3.5f, "stage"}}, true);
    fRegistryMC.add("MC/TruthAO2D/hStageConsistency",
                    "stage consistency (expect all at 0);N(V0 built but legs not found);counts",
                    kTH1D, {AxisSpec{20, -0.5f, 19.5f, "N_{bad}"}}, true);

    // ─── Single-leg diagnostics ───────────────────────────────────────────────
    fRegistryMC.add("MC/LegDiag/hLegDRtrue_vs_pt_legFound", "leg found: #Delta R^{true} vs p_{T};p_{T,#gamma}^{true} (GeV/c);#Delta R_{e^{+}e^{-}}^{true}", kTH2D, {axisPt, axLegDR}, true);
    fRegistryMC.add("MC/LegDiag/hLegDRtrue_vs_pt_legMissing", "leg missing: #Delta R^{true} vs p_{T};p_{T,#gamma}^{true} (GeV/c);#Delta R_{e^{+}e^{-}}^{true}", kTH2D, {axisPt, axLegDR}, true);
    fRegistryMC.add("MC/LegDiag/hLegDEta_legFound_vs_pt", "leg found: |#Delta#eta| vs p_{T};p_{T,#gamma}^{true} (GeV/c);|#Delta#eta_{e^{+}e^{-}}|", kTH2D, {axisPt, AxisSpec{100, 0.f, 0.5f, "|#Delta#eta_{legs}|"}}, true);
    fRegistryMC.add("MC/LegDiag/hLegDEta_legMissing_vs_pt", "leg missing: |#Delta#eta| vs p_{T};p_{T,#gamma}^{true} (GeV/c);|#Delta#eta_{e^{+}e^{-}}|", kTH2D, {axisPt, AxisSpec{100, 0.f, 0.5f, "|#Delta#eta_{legs}|"}}, true);
    fRegistryMC.add("MC/LegDiag/hLegDPhi_legFound_vs_pt", "leg found: |#Delta#phi| vs p_{T};p_{T,#gamma}^{true} (GeV/c);|#Delta#phi_{e^{+}e^{-}}| (rad)", kTH2D, {axisPt, AxisSpec{100, 0.f, 0.5f, "|#Delta#phi_{legs}|"}}, true);
    fRegistryMC.add("MC/LegDiag/hLegDPhi_legMissing_vs_pt", "leg missing: |#Delta#phi| vs p_{T};p_{T,#gamma}^{true} (GeV/c);|#Delta#phi_{e^{+}e^{-}}| (rad)", kTH2D, {axisPt, AxisSpec{100, 0.f, 0.5f, "|#Delta#phi_{legs}|"}}, true);
    fRegistryMC.add("MC/LegDiag/hAlphaTrue_legFound_vs_pt", "leg found: #alpha^{true} vs p_{T};p_{T,#gamma}^{true} (GeV/c);#alpha^{true}", kTH2D, {axisPt, axAlpha}, true);
    fRegistryMC.add("MC/LegDiag/hAlphaTrue_legMissing_vs_pt", "leg missing: #alpha^{true} vs p_{T};p_{T,#gamma}^{true} (GeV/c);#alpha^{true}", kTH2D, {axisPt, axAlpha}, true);
    fRegistryMC.add("MC/LegDiag/hAlpha_vs_legDR_legMissing", "leg missing: #alpha^{true} vs #Delta R^{true};#Delta R_{e^{+}e^{-}}^{true};#alpha^{true}", kTH2D, {axLegDR, axAlpha}, true);

    // ─── Pair-level leg diagnostics ───────────────────────────────────────────
    fRegistryMC.add("MC/LegDiag/hNLegsPair_vs_kT", "N legs found per pair vs k_{T};k_{T} (GeV/c);N_{legs found} (0-4)", kTH2D, {axisKt, AxisSpec{5, -0.5f, 4.5f, "N_{legs found}"}}, true);
    fRegistryMC.add("MC/LegDiag/hMissingLegPt_vs_kT", "missing leg p_{T}^{true} vs pair k_{T};k_{T} (GeV/c);p_{T,leg}^{true} (GeV/c)", kTH2D, {axisKt, AxisSpec{100, 0.f, 0.5f, "p_{T,leg}^{true} (GeV/c)"}}, true);
    fRegistryMC.add("MC/LegDiag/hMissingLegRconv_vs_kT", "missing leg R_{conv}^{true} vs pair k_{T};k_{T} (GeV/c);R_{conv}^{true} (cm)", kTH2D, {axisKt, axisR}, true);

    // ─── Cross-built V0 pairs ─────────────────────────────────────────────────
    fRegistryMC.add("MC/PairCrossBuild/hSparse_DEtaDPhi_kT",
                    "cross-built V0 pairs: #Delta#eta,#Delta#phi,k_{T};"
                    "#Delta#eta_{#gamma#gamma};#Delta#phi_{#gamma#gamma} (rad);k_{T} (GeV/c)",
                    kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisKt}, true);
    fRegistryMC.add("MC/PairCrossBuild/hStageOut_vs_kT",
                    "cross-built pairs: N correctly built vs k_{T};"
                    "k_{T} (GeV/c);N photons correctly built (0/1/2)",
                    kTH2D, {axisKt, AxisSpec{3, -0.5f, 2.5f, "N correctly built"}}, true);
  }

  template <PairTruthType TruthT, bool IsMix>
  inline void fillMCPairQATyped(PairQAObservables const& obs, bool doSparse, bool doMCQA)
  {
    constexpr auto base = []() constexpr -> const char* {
      if constexpr (!IsMix) {
        if constexpr (TruthT == PairTruthType::TrueTrueDistinct)
          return "Pair/same/MC/TrueTrueDistinct/";
        if constexpr (TruthT == PairTruthType::TrueTrueSamePhoton)
          return "Pair/same/MC/TrueTrueSamePhoton/";
        if constexpr (TruthT == PairTruthType::SharedMcLeg)
          return "Pair/same/MC/SharedMcLeg/";
        if constexpr (TruthT == PairTruthType::TrueFake)
          return "Pair/same/MC/TrueFake/";
        if constexpr (TruthT == PairTruthType::FakeFake)
          return "Pair/same/MC/FakeFake/";
        return "Pair/same/MC/Pi0Daughters/";
      } else {
        if constexpr (TruthT == PairTruthType::TrueTrueDistinct)
          return "Pair/mix/MC/TrueTrueDistinct/";
        if constexpr (TruthT == PairTruthType::TrueTrueSamePhoton)
          return "Pair/mix/MC/TrueTrueSamePhoton/";
        if constexpr (TruthT == PairTruthType::SharedMcLeg)
          return "Pair/mix/MC/SharedMcLeg/";
        if constexpr (TruthT == PairTruthType::TrueFake)
          return "Pair/mix/MC/TrueFake/";
        if constexpr (TruthT == PairTruthType::FakeFake)
          return "Pair/mix/MC/FakeFake/";
        return "Pair/mix/MC/Pi0Daughters/";
      }
    }();
    if (doMCQA) {
      fRegistryPairMC.fill(HIST(base) + HIST("hDEtaDPhi"), obs.deta, obs.dphi);
      fRegistryPairMC.fill(HIST(base) + HIST("hQinv"), obs.qinv);
      fRegistryPairMC.fill(HIST(base) + HIST("hDeltaR"), obs.deltaR);
      fRegistryPairMC.fill(HIST(base) + HIST("hDeltaZ"), obs.deltaZ);
      fRegistryPairMC.fill(HIST(base) + HIST("hDeltaR3D"), obs.deltaR3D);
    }
    if (doSparse)
      fRegistryPairMC.fill(HIST(base) + HIST("hSparse_DEtaDPhi_kT"), obs.deta, obs.dphi, obs.kt);
    constexpr auto summaryDir = IsMix ? "Pair/mix/MC/" : "Pair/same/MC/";
    if (doMCQA) {
      fRegistryPairMC.fill(HIST(summaryDir) + HIST("hTruthTypeVsQinv"), obs.qinv, static_cast<int>(TruthT));
      fRegistryPairMC.fill(HIST(summaryDir) + HIST("hTruthTypeVsKt"), obs.kt, static_cast<int>(TruthT));
    }
  }

  template <bool IsMix>
  inline void fillMCPairQA(PairTruthType truthType, PairQAObservables const& obs, bool doSparse, bool doMCQA)
  {
    switch (truthType) {
      case PairTruthType::TrueTrueDistinct:
        fillMCPairQATyped<PairTruthType::TrueTrueDistinct, IsMix>(obs, doSparse, doMCQA);
        break;
      case PairTruthType::TrueTrueSamePhoton:
        fillMCPairQATyped<PairTruthType::TrueTrueSamePhoton, IsMix>(obs, doSparse, doMCQA);
        break;
      case PairTruthType::SharedMcLeg:
        fillMCPairQATyped<PairTruthType::SharedMcLeg, IsMix>(obs, doSparse, doMCQA);
        break;
      case PairTruthType::TrueFake:
        fillMCPairQATyped<PairTruthType::TrueFake, IsMix>(obs, doSparse, doMCQA);
        break;
      case PairTruthType::FakeFake:
        fillMCPairQATyped<PairTruthType::FakeFake, IsMix>(obs, doSparse, doMCQA);
        break;
      case PairTruthType::Pi0Daughters:
        fillMCPairQATyped<PairTruthType::Pi0Daughters, IsMix>(obs, doSparse, doMCQA);
        break;
      default:
        break;
    }
  }

  template <PairTruthType TruthT, bool IsMix>
  inline void fillMCPairQAFullRangeTyped(PairQAObservables const& obs)
  {
    constexpr auto base = []() constexpr -> const char* {
      if constexpr (!IsMix) {
        if constexpr (TruthT == PairTruthType::TrueTrueDistinct)
          return "Pair/same/MC/TrueTrueDistinct/";
        if constexpr (TruthT == PairTruthType::TrueTrueSamePhoton)
          return "Pair/same/MC/TrueTrueSamePhoton/";
        if constexpr (TruthT == PairTruthType::SharedMcLeg)
          return "Pair/same/MC/SharedMcLeg/";
        if constexpr (TruthT == PairTruthType::TrueFake)
          return "Pair/same/MC/TrueFake/";
        if constexpr (TruthT == PairTruthType::FakeFake)
          return "Pair/same/MC/FakeFake/";
        return "Pair/same/MC/Pi0Daughters/";
      } else {
        if constexpr (TruthT == PairTruthType::TrueTrueDistinct)
          return "Pair/mix/MC/TrueTrueDistinct/";
        if constexpr (TruthT == PairTruthType::TrueTrueSamePhoton)
          return "Pair/mix/MC/TrueTrueSamePhoton/";
        if constexpr (TruthT == PairTruthType::SharedMcLeg)
          return "Pair/mix/MC/SharedMcLeg/";
        if constexpr (TruthT == PairTruthType::TrueFake)
          return "Pair/mix/MC/TrueFake/";
        if constexpr (TruthT == PairTruthType::FakeFake)
          return "Pair/mix/MC/FakeFake/";
        return "Pair/mix/MC/Pi0Daughters/";
      }
    }();
    fRegistryPairMC.fill(HIST(base) + HIST("hDeltaRVsQinv"), obs.qinv, obs.deltaR);
    fRegistryPairMC.fill(HIST(base) + HIST("hDeltaZVsQinv"), obs.qinv, obs.deltaZ);
    fRegistryPairMC.fill(HIST(base) + HIST("hDeltaR3DVsQinv"), obs.qinv, obs.deltaR3D);
    const bool fillDRDZ = ((TruthT == PairTruthType::TrueTrueDistinct) ? mctruthSparse.cfgFillDRDZQinvTrueTrueDistinct.value : (TruthT == PairTruthType::TrueTrueSamePhoton) ? mctruthSparse.cfgFillDRDZQinvTrueTrueSamePhoton.value
                                                                                                                             : (TruthT == PairTruthType::SharedMcLeg)        ? mctruthSparse.cfgFillDRDZQinvSharedMcLeg.value
                                                                                                                             : (TruthT == PairTruthType::TrueFake)           ? mctruthSparse.cfgFillDRDZQinvTrueFake.value
                                                                                                                             : (TruthT == PairTruthType::FakeFake)           ? mctruthSparse.cfgFillDRDZQinvFakeFake.value
                                                                                                                                                                             : mctruthSparse.cfgFillDRDZQinvPi0Daughters.value);
    if (fillDRDZ)
      fRegistryPairMC.fill(HIST(base) + HIST("hSparseDeltaRDeltaZQinv"), obs.deltaR, obs.deltaZ, obs.qinv);
    const bool enabled = ((TruthT == PairTruthType::TrueTrueDistinct) ? mctruthSparse.cfgFillDEtaDPhiVsQinvTrueTrueDistinct.value : (TruthT == PairTruthType::TrueTrueSamePhoton) ? mctruthSparse.cfgFillDEtaDPhiVsQinvTrueTrueSamePhoton.value
                                                                                                                                  : (TruthT == PairTruthType::SharedMcLeg)        ? mctruthSparse.cfgFillDEtaDPhiVsQinvSharedMcLeg.value
                                                                                                                                  : (TruthT == PairTruthType::TrueFake)           ? mctruthSparse.cfgFillDEtaDPhiVsQinvTrueFake.value
                                                                                                                                  : (TruthT == PairTruthType::FakeFake)           ? mctruthSparse.cfgFillDEtaDPhiVsQinvFakeFake.value
                                                                                                                                                                                  : mctruthSparse.cfgFillDEtaDPhiVsQinvPi0Daughters.value);
    if (enabled)
      fRegistryPairMC.fill(HIST(base) + HIST("hDEtaDPhiVsQinv"), obs.deta, obs.dphi, obs.qinv);
  }

  template <bool IsMix>
  inline void fillMCPairQAFullRange(PairTruthType truthType, PairQAObservables const& obs)
  {
    switch (truthType) {
      case PairTruthType::TrueTrueDistinct:
        fillMCPairQAFullRangeTyped<PairTruthType::TrueTrueDistinct, IsMix>(obs);
        break;
      case PairTruthType::TrueTrueSamePhoton:
        fillMCPairQAFullRangeTyped<PairTruthType::TrueTrueSamePhoton, IsMix>(obs);
        break;
      case PairTruthType::SharedMcLeg:
        fillMCPairQAFullRangeTyped<PairTruthType::SharedMcLeg, IsMix>(obs);
        break;
      case PairTruthType::TrueFake:
        fillMCPairQAFullRangeTyped<PairTruthType::TrueFake, IsMix>(obs);
        break;
      case PairTruthType::FakeFake:
        fillMCPairQAFullRangeTyped<PairTruthType::FakeFake, IsMix>(obs);
        break;
      case PairTruthType::Pi0Daughters:
        fillMCPairQAFullRangeTyped<PairTruthType::Pi0Daughters, IsMix>(obs);
        break;
      default:
        break;
    }
  }
  template <typename TCollision>
  void fillPairHistogramNoLabel(TCollision const& /*collision*/,
                                ROOT::Math::PtEtaPhiMVector v1,
                                ROOT::Math::PtEtaPhiMVector v2)
  {
    float rndm = std::pow(-1, dist01(engine) % 2);
    auto k12 = 0.5 * (v1 + v2);
    float kt = k12.Pt();
    float qinv = -(((v1 - v2) * rndm).M());
    ROOT::Math::XYZVector uv_out(k12.Px() / k12.Pt(), k12.Py() / k12.Pt(), 0);
    ROOT::Math::XYZVector uv_long(0, 0, 1);
    ROOT::Math::XYZVector uv_side = uv_out.Cross(uv_long);
    ROOT::Math::PxPyPzEVector v1c(v1), v2c(v2);
    float beta_z = (v1 + v2).Beta() * std::cos((v1 + v2).Theta());
    ROOT::Math::Boost bst_z(0, 0, -beta_z);
    auto q12_lcms = bst_z((v1c - v2c) * rndm);
    auto q3_lcms = q12_lcms.Vect();
    float qabs_lcms = q3_lcms.R();
    float qout_lcms = q3_lcms.Dot(uv_out);
    float qside_lcms = q3_lcms.Dot(uv_side);
    float qlong_lcms = q3_lcms.Dot(uv_long);
    if (hbtanalysis.cfgDo3D) {
      fRegistryPairMC.fill(HIST("Pair/same/MC/NoLabel/CF_3D"),
                           std::fabs(qout_lcms), std::fabs(qside_lcms), std::fabs(qlong_lcms), kt);
      if (hbtanalysis.cfgDo2D)
        fRegistryPairMC.fill(HIST("Pair/same/MC/NoLabel/CF_2D"), std::fabs(qout_lcms), std::fabs(qinv), kt);
    } else {
      fRegistryPairMC.fill(HIST("Pair/same/MC/NoLabel/CF_1D"), hbtanalysis.cfgUseLCMS ? qabs_lcms : qinv, kt);
    }
  }

  template <typename TCollisions,
            typename TPhotons1, typename TPhotons2,
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
      const float cent[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (cent[mixing.cfgCentEstimator] < centralitySelection.cfgCentMin || centralitySelection.cfgCentMax < cent[mixing.cfgCentEstimator])
        continue;
      const std::array<float, 7> epArr = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(),
                                          collision.ep2fv0a(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      const float ep2 = epArr[mixing.cfgEP2EstimatorForMix];
      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision))
        continue;
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0);
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0);
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);
      const float occupancy = (mixing.cfgOccupancyEstimator == 1)
                                ? static_cast<float>(collision.trackOccupancyInTimeRange())
                                : collision.ft0cOccupancyInTimeRange();
      const float centForQA = cent[mixing.cfgCentEstimator];
      const int zbin = binOf(ztxBinEdges, collision.posZ()), centbin = binOf(centBinEdges, centForQA);
      const int epbin = binOf(epBinEgdes, ep2), occbin = binOf(occBinEdges, occupancy);
      auto keyBin = std::make_tuple(zbin, centbin, epbin, occbin);
      auto keyDFCollision = std::make_pair(ndf, collision.globalIndex());
      auto photons1Coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      auto photons2Coll = photons2.sliceBy(perCollision2, collision.globalIndex());
      if (qaflags.doSinglePhotonQa)
        for (const auto& g : photons1Coll)
          if (cut1.template IsSelected<decltype(g), TSubInfos1>(g))
            fillSinglePhotonQAStep<0>(g);
      std::unordered_set<int> idsAfterDR, idsAfterRZ, idsAfterEllipse;
      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1Coll, photons2Coll))) {
        if (!cut1.template IsSelected<decltype(g1), TSubInfos1>(g1) ||
            !cut2.template IsSelected<decltype(g2), TSubInfos2>(g2))
          continue;
        const auto pos1 = g1.template posTrack_as<TSubInfos1>(), ele1 = g1.template negTrack_as<TSubInfos1>();
        const auto pos2 = g2.template posTrack_as<TSubInfos2>(), ele2 = g2.template negTrack_as<TSubInfos2>();
        if (pos1.trackId() == pos2.trackId() || pos1.trackId() == ele2.trackId() ||
            ele1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId())
          continue;
        if (!passAsymmetryCut(g1.pt(), g2.pt()))
          continue;
        auto obs = buildPairQAObservables(g1, g2);

        if (!obs.valid)
          continue;
        const bool doQA = passQinvQAGate(obs.qinv), doFR = passQinvFullRangeGate(obs.qinv);
        if (doQA)
          fillPairQAStep<0, 0>(obs, centForQA, occupancy);
        if (doFR)
          fillFullRangeDeltaRCosOA<0>(obs.qinv, obs.drOverCosOA);
        fRegistry.fill(HIST("Pair/same/hDeltaRCosOA"), obs.drOverCosOA);
        if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA)
          continue;
        idsAfterDR.insert(g1.globalIndex());
        idsAfterDR.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 1>(obs, centForQA, occupancy);
        if (!passRZCut(obs.deltaR, obs.deltaZ))
          continue;
        idsAfterRZ.insert(g1.globalIndex());
        idsAfterRZ.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 2>(obs, centForQA, occupancy);
        if (isInsideEllipse(obs.deta, obs.dphi))
          continue;
        idsAfterEllipse.insert(g1.globalIndex());
        idsAfterEllipse.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 3>(obs, centForQA, occupancy);
        if (doFR)
          fillFullRangeQA<0>(obs, centForQA, occupancy);
        fillPairHistogram<0>(collision, obs.v1, obs.v2, 1.f);
        ndiphoton++;
        fRegistry.fill(HIST("Pair/same/hPhi_lowerPtV0"),
                       (g1.pt() < g2.pt()) ? g1.phi() : g2.phi());

        fillEtaTopologyHisto(
          classifyEtaTopology(g1.eta(), g2.eta(),
                              pos1.eta(), ele1.eta(),
                              pos2.eta(), ele2.eta()),
          obs.deta, obs.dphi, obs.kt);
        auto addToPool = [&](auto const& g) {
          if (usedPhotonIdsPerCol.insert(g.globalIndex()).second) {
            EMPair gtmp(g.pt(), g.eta(), g.phi(), 0.f);
            gtmp.setConversionPointXYZ(g.vx(), g.vy(), g.vz());
            emh1->AddTrackToEventPool(keyDFCollision, gtmp);
          }
        };
        addToPool(g1);
        addToPool(g2);
      }
      if (qaflags.doSinglePhotonQa) {
        for (const auto& g : photons1Coll) {
          if (cut1.template IsSelected<decltype(g), TSubInfos1>(g)) {
            const int gid = g.globalIndex();
            if (idsAfterDR.count(gid)) {
              fillSinglePhotonQAStep<1>(g);
            }
            if (idsAfterRZ.count(gid)) {
              fillSinglePhotonQAStep<2>(g);
            }
            if (idsAfterEllipse.count(gid)) {
              fillSinglePhotonQAStep<3>(g);
            }
          }
        }
      }
      usedPhotonIdsPerCol.clear();
      if (!mixing.cfgDoMix || ndiphoton == 0)
        continue;
      auto selectedPhotons = emh1->GetTracksPerCollision(keyDFCollision);
      auto poolIDs = emh1->GetCollisionIdsFromEventPool(keyBin);
      for (const auto& mixID : poolIDs) {
        if (mixID.second == collision.globalIndex() && mixID.first == ndf)
          continue;
        const uint64_t bcMix = mapMixedEventIdToGlobalBC[mixID];
        const uint64_t diffBC = std::max(collision.globalBC(), bcMix) - std::min(collision.globalBC(), bcMix);
        fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
        if (diffBC < mixing.ndiffBCMix)
          continue;
        auto poolPhotons = emh1->GetTracksPerCollision(mixID);
        for (const auto& g1 : selectedPhotons)
          for (const auto& g2 : poolPhotons) {
            if (!passAsymmetryCut(g1.pt(), g2.pt()))
              continue;
            auto obs = buildPairQAObservables(g1, g2);
            if (!obs.valid)
              continue;
            const bool doQA = passQinvQAGate(obs.qinv), doFR = passQinvFullRangeGate(obs.qinv);
            if (doQA)
              fillPairQAStep<1, 0>(obs, centForQA, occupancy);
            if (doFR)
              fillFullRangeDeltaRCosOA<1>(obs.qinv, obs.drOverCosOA);
            fRegistry.fill(HIST("Pair/mix/hDeltaRCosOA"), obs.drOverCosOA);
            if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA)
              continue;
            if (doQA)
              fillPairQAStep<1, 1>(obs, centForQA, occupancy);
            if (!passRZCut(obs.deltaR, obs.deltaZ))
              continue;
            if (doQA)
              fillPairQAStep<1, 2>(obs, centForQA, occupancy);
            if (isInsideEllipse(obs.deta, obs.dphi))
              continue;
            if (doQA)
              fillPairQAStep<1, 3>(obs, centForQA, occupancy);
            if (doFR)
              fillFullRangeQA<1>(obs, centForQA, occupancy);
            fillPairHistogram<1>(collision, obs.v1, obs.v2, 1.f);
            fRegistry.fill(HIST("Pair/mix/hPhi_lowerPtV0"),
                           (g1.pt() < g2.pt()) ? g1.phi() : g2.phi());
          }
      }
      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(keyBin, keyDFCollision);
        emh2->AddCollisionIdAtLast(keyBin, keyDFCollision);
        mapMixedEventIdToGlobalBC[keyDFCollision] = collision.globalBC();
      }
    }
  }

  template <soa::is_table TCollisions, soa::is_table TPhotons,
            soa::is_table TLegs, soa::is_table TMCParticles,
            typename TPreslice, typename TCut>
  void runPairingMC(TCollisions const& collisions, TPhotons const& photons,
                    TLegs const& /*legs*/, TMCParticles const& mcParticles,
                    TPreslice const& perCollision, TCut const& cut)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      int ndiphoton = 0;
      const float cent[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (cent[mixing.cfgCentEstimator] < centralitySelection.cfgCentMin || centralitySelection.cfgCentMax < cent[mixing.cfgCentEstimator])
        continue;
      const std::array<float, 7> epArr = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(),
                                          collision.ep2fv0a(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      const float ep2 = epArr[mixing.cfgEP2EstimatorForMix];
      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision))
        continue;
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0);
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0);
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);
      const float occupancy = (mixing.cfgOccupancyEstimator == 1)
                                ? static_cast<float>(collision.trackOccupancyInTimeRange())
                                : collision.ft0cOccupancyInTimeRange();
      const float centForQA = cent[mixing.cfgCentEstimator];
      const int zbin = binOf(ztxBinEdges, collision.posZ()), centbin = binOf(centBinEdges, centForQA);
      const int epbin = binOf(epBinEgdes, ep2), occbin = binOf(occBinEdges, occupancy);
      auto keyBin = std::make_tuple(zbin, centbin, epbin, occbin);
      auto keyDFCollision = std::make_pair(ndf, collision.globalIndex());
      auto photonsColl = photons.sliceBy(perCollision, collision.globalIndex());
      if (qaflags.doSinglePhotonQa)
        for (const auto& g : photonsColl)
          if (cut.template IsSelected<decltype(g), TLegs>(g))
            fillSinglePhotonQAStep<0>(g);
      std::unordered_set<int> idsAfterDR, idsAfterRZ, idsAfterEllipse;
      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsColl, photonsColl))) {
        if (!cut.template IsSelected<decltype(g1), TLegs>(g1) || !cut.template IsSelected<decltype(g2), TLegs>(g2))
          continue;
        const auto pos1 = g1.template posTrack_as<TLegs>(), ele1 = g1.template negTrack_as<TLegs>();
        const auto pos2 = g2.template posTrack_as<TLegs>(), ele2 = g2.template negTrack_as<TLegs>();
        if (pos1.trackId() == pos2.trackId() || pos1.trackId() == ele2.trackId() ||
            ele1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId())
          continue;
        const auto mc1 = buildPhotonMCInfo<decltype(g1), TLegs>(g1, mcParticles);
        const auto mc2 = buildPhotonMCInfo<decltype(g2), TLegs>(g2, mcParticles);
        auto truthType = classifyPairTruth(mc1, mc2);
        if (truthType == PairTruthType::TrueTrueDistinct && isPi0DaughterPair(mc1, mc2, mcParticles))
          truthType = PairTruthType::Pi0Daughters;
        if (!passAsymmetryCut(g1.pt(), g2.pt()))
          continue;
        auto obs = buildPairQAObservables(g1, g2);
        if (!obs.valid)
          continue;
        const bool doQA = passQinvQAGate(obs.qinv), doFR = passQinvFullRangeGate(obs.qinv);
        if (doQA)
          fillPairQAStep<0, 0>(obs, centForQA, occupancy);
        if (doFR)
          fillFullRangeDeltaRCosOA<0>(obs.qinv, obs.drOverCosOA);
        fRegistry.fill(HIST("Pair/same/hDeltaRCosOA"), obs.drOverCosOA);
        if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA)
          continue;
        idsAfterDR.insert(g1.globalIndex());
        idsAfterDR.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 1>(obs, centForQA, occupancy);
        if (!passRZCut(obs.deltaR, obs.deltaZ))
          continue;
        idsAfterRZ.insert(g1.globalIndex());
        idsAfterRZ.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 2>(obs, centForQA, occupancy);
        if (isInsideEllipse(obs.deta, obs.dphi))
          continue;
        idsAfterEllipse.insert(g1.globalIndex());
        idsAfterEllipse.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 3>(obs, centForQA, occupancy);
        if (doFR)
          fillFullRangeQA<0>(obs, centForQA, occupancy);
        fillPairHistogram<0>(collision, obs.v1, obs.v2, 1.f);
        fRegistry.fill(HIST("Pair/same/hPhi_lowerPtV0"),
                       (g1.pt() < g2.pt()) ? g1.phi() : g2.phi());

        fillEtaTopologyHisto(
          classifyEtaTopology(g1.eta(), g2.eta(),
                              pos1.eta(), ele1.eta(),
                              pos2.eta(), ele2.eta()),
          obs.deta, obs.dphi, obs.kt);
        ndiphoton++;
        if (!mc1.hasMC || !mc2.hasMC) {
          fillPairHistogramNoLabel(collision, obs.v1, obs.v2);
          fRegistryPairMC.fill(HIST("Pair/same/MC/NoLabel/hDEtaDPhi"), obs.deta, obs.dphi);
          fRegistryPairMC.fill(HIST("Pair/same/MC/NoLabel/hKt"), obs.kt);
          fRegistryPairMC.fill(HIST("Pair/same/MC/NoLabel/hQinv"), obs.qinv);
        } else {
          const bool doMCQA = passQinvMCQAGate(obs.qinv);
          fillMCPairQA<false>(truthType, obs, doQA, doMCQA);
          if (doFR)
            fillMCPairQAFullRange<false>(truthType, obs);
          const bool isTruePair = (truthType == PairTruthType::TrueTrueDistinct ||
                                   truthType == PairTruthType::TrueTrueSamePhoton ||
                                   truthType == PairTruthType::Pi0Daughters);
          if (isTruePair) {
            fRegistryPairMC.fill(HIST("Pair/same/MC/hDEtaDPhi_truePairs"), obs.deta, obs.dphi);
            if (doMCQA) {
              fRegistryPairMC.fill(HIST("Pair/same/MC/hSparse_DEtaDPhi_kT_truePairs"), obs.deta, obs.dphi, obs.kt);
              fRegistryPairMC.fill(HIST("Pair/same/MC/hSparse_DEtaDPhi_qinv_truePairs"), obs.deta, obs.dphi, obs.qinv);
            }
          } else {
            fRegistryPairMC.fill(HIST("Pair/same/MC/hDEtaDPhi_fakePairs"), obs.deta, obs.dphi);
            if (doMCQA) {
              fRegistryPairMC.fill(HIST("Pair/same/MC/hSparse_DEtaDPhi_kT_fakePairs"), obs.deta, obs.dphi, obs.kt);
              fRegistryPairMC.fill(HIST("Pair/same/MC/hSparse_DEtaDPhi_qinv_fakePairs"), obs.deta, obs.dphi, obs.qinv);
            }
          }

          switch (truthType) {
            case PairTruthType::TrueTrueDistinct:
              fillPairHistogramMC<0, PairTruthType::TrueTrueDistinct>(collision, obs.v1, obs.v2);
              break;
            case PairTruthType::TrueTrueSamePhoton:
              fillPairHistogramMC<0, PairTruthType::TrueTrueSamePhoton>(collision, obs.v1, obs.v2);
              break;
            case PairTruthType::SharedMcLeg:
              fillPairHistogramMC<0, PairTruthType::SharedMcLeg>(collision, obs.v1, obs.v2);
              break;
            case PairTruthType::TrueFake:
              fillPairHistogramMC<0, PairTruthType::TrueFake>(collision, obs.v1, obs.v2);
              break;
            case PairTruthType::FakeFake:
              fillPairHistogramMC<0, PairTruthType::FakeFake>(collision, obs.v1, obs.v2);
              break;
            case PairTruthType::Pi0Daughters:
              fillPairHistogramMC<0, PairTruthType::Pi0Daughters>(collision, obs.v1, obs.v2);
              break;
            default:
              break;
          }
        }

        auto addToPool = [&](auto const& g) {
          if (usedPhotonIdsPerCol.insert(g.globalIndex()).second) {
            EMPair gtmp(g.pt(), g.eta(), g.phi(), 0.f);
            gtmp.setConversionPointXYZ(g.vx(), g.vy(), g.vz());
            emh1->AddTrackToEventPool(keyDFCollision, gtmp);
          }
        };
        addToPool(g1);
        addToPool(g2);
      }
      if (qaflags.doSinglePhotonQa) {
        for (const auto& g : photonsColl) {
          if (cut.template IsSelected<decltype(g), TLegs>(g)) {
            const int gid = g.globalIndex();
            if (idsAfterDR.count(gid)) {
              fillSinglePhotonQAStep<1>(g);
            }
            if (idsAfterRZ.count(gid)) {
              fillSinglePhotonQAStep<2>(g);
            }
            if (idsAfterEllipse.count(gid)) {
              fillSinglePhotonQAStep<3>(g);
            }
          }
        }
      }
      usedPhotonIdsPerCol.clear();
      if (!mixing.cfgDoMix || ndiphoton == 0)
        continue;
      auto selectedPhotons = emh1->GetTracksPerCollision(keyDFCollision);
      auto poolIDs = emh1->GetCollisionIdsFromEventPool(keyBin);
      for (const auto& mixID : poolIDs) {
        if (mixID.second == collision.globalIndex() && mixID.first == ndf)
          continue;
        const uint64_t bcMix = mapMixedEventIdToGlobalBC[mixID];
        const uint64_t diffBC = std::max(collision.globalBC(), bcMix) - std::min(collision.globalBC(), bcMix);
        fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
        if (diffBC < mixing.ndiffBCMix)
          continue;
        auto poolPhotons = emh1->GetTracksPerCollision(mixID);
        for (const auto& g1 : selectedPhotons)
          for (const auto& g2 : poolPhotons) {
            if (!passAsymmetryCut(g1.pt(), g2.pt()))
              continue;
            auto obs = buildPairQAObservables(g1, g2);

            if (!obs.valid)
              continue;
            const bool doQA = passQinvQAGate(obs.qinv), doFR = passQinvFullRangeGate(obs.qinv);
            if (doQA)
              fillPairQAStep<1, 0>(obs, centForQA, occupancy);
            if (doFR)
              fillFullRangeDeltaRCosOA<1>(obs.qinv, obs.drOverCosOA);
            fRegistry.fill(HIST("Pair/mix/hDeltaRCosOA"), obs.drOverCosOA);
            if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA)
              continue;
            if (doQA)
              fillPairQAStep<1, 1>(obs, centForQA, occupancy);
            if (!passRZCut(obs.deltaR, obs.deltaZ))
              continue;
            if (doQA)
              fillPairQAStep<1, 2>(obs, centForQA, occupancy);
            if (isInsideEllipse(obs.deta, obs.dphi))
              continue;
            if (doQA)
              fillPairQAStep<1, 3>(obs, centForQA, occupancy);
            if (doFR)
              fillFullRangeQA<1>(obs, centForQA, occupancy);
            fillPairHistogram<1>(collision, obs.v1, obs.v2, 1.f);
            fRegistry.fill(HIST("Pair/mix/hPhi_lowerPtV0"),
                           (g1.pt() < g2.pt()) ? g1.phi() : g2.phi());
          }
      }
      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(keyBin, keyDFCollision);
        emh2->AddCollisionIdAtLast(keyBin, keyDFCollision);
        mapMixedEventIdToGlobalBC[keyDFCollision] = collision.globalBC();
      }
    }
  }
  template <soa::is_table TCollisions, soa::is_table TPhotons,
            soa::is_table TLegs, soa::is_table TMCParticles,
            soa::is_table TMCEvents,
            typename TPresliceMCParts, typename TPresliceLegs, typename TCut>
  void runTruthEfficiency(TCollisions const& collisions,
                          TPhotons const& v0photons,
                          TLegs const& v0legs,
                          TMCParticles const& emmcParticles,
                          TMCEvents const& /*mcEvents*/,
                          TPresliceMCParts const& perMCCollision,
                          TPresliceLegs const& perCollisionLegs,
                          TCut const& cut)
  {
    auto wrapPhi = [](float dphi) -> float {
      return RecoDecay::constrainAngle(dphi, -o2::constants::math::PI);
    };

    std::unordered_set<int> mcIdsWithAnyV0Leg;
    mcIdsWithAnyV0Leg.reserve(v0legs.size() * 2);
    for (const auto& leg : v0legs) {
      if (leg.has_emmcparticle())
        mcIdsWithAnyV0Leg.insert(leg.emmcparticleId());
    }

    for (const auto& collision : collisions) {
      if (!fEMEventCut.IsSelected(collision))
        continue;
      const float cent[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (cent[mixing.cfgCentEstimator] < centralitySelection.cfgCentMin ||
          centralitySelection.cfgCentMax < cent[mixing.cfgCentEstimator])
        continue;
      if (!collision.has_emmcevent())
        continue;

      const int64_t thisCollisionId = collision.globalIndex();
      const int mcEventId = collision.template emmcevent_as<TMCEvents>().globalIndex();

      auto recoPhotonsColl = v0photons.sliceBy(perCollisionPCM, thisCollisionId);
      auto emmcPartsColl = emmcParticles.sliceBy(perMCCollision, mcEventId);
      auto legsColl = v0legs.sliceBy(perCollisionLegs, thisCollisionId);

      std::unordered_set<int> legIdsThisCollision;
      legIdsThisCollision.reserve(legsColl.size());
      for (const auto& leg : legsColl)
        if (leg.has_emmcparticle())
          legIdsThisCollision.insert(leg.emmcparticleId());

      std::unordered_map<int, std::unordered_set<int>> crossBuildMap;

      struct PhotonRecoInfo {
        bool hasV0 = false, passesCut = false;
      };
      std::unordered_map<int, PhotonRecoInfo> gammaRecoMap;
      gammaRecoMap.reserve(recoPhotonsColl.size());

      for (const auto& g : recoPhotonsColl) {
        const auto pos = g.template posTrack_as<TLegs>();
        const auto neg = g.template negTrack_as<TLegs>();
        if (pos.collisionId() != thisCollisionId || neg.collisionId() != thisCollisionId)
          continue;
        if (!pos.has_emmcparticle() || !neg.has_emmcparticle())
          continue;
        const auto mcPos = pos.template emmcparticle_as<aod::EMMCParticles>();
        const auto mcNeg = neg.template emmcparticle_as<aod::EMMCParticles>();
        if (!mcPos.has_mothers() || !mcNeg.has_mothers())
          continue;
        const int posMotherId = mcPos.mothersIds()[0], negMotherId = mcNeg.mothersIds()[0];
        if (posMotherId != negMotherId) {
          const auto posMother = emmcParticles.iteratorAt(posMotherId);
          const auto negMother = emmcParticles.iteratorAt(negMotherId);
          if (posMother.pdgCode() == kGamma && negMother.pdgCode() == kGamma) {
            crossBuildMap[posMotherId].insert(negMotherId);
            crossBuildMap[negMotherId].insert(posMotherId);
          }
          continue;
        }
        const int gammaId = posMotherId;
        if (emmcParticles.iteratorAt(gammaId).pdgCode() != kGamma)
          continue;
        const bool passes = cut.template IsSelected<std::decay_t<decltype(g)>, TLegs>(g);
        auto& info = gammaRecoMap[gammaId];
        info.hasV0 = true;
        info.passesCut = info.passesCut || passes;
      }

      // ─── Build true gamma list ────────────────────────────────────────────────
      std::vector<TruthGamma> trueGammas;
      trueGammas.reserve(32);

      for (const auto& g : emmcPartsColl) {
        if (g.pdgCode() != kGamma)
          continue;
        if (!g.isPhysicalPrimary() && !g.producedByGenerator())
          continue;
        if (std::fabs(g.eta()) > pcmcuts.cfgMaxEtaV0.value)
          continue;
        const float mcV0PtMin = (mctruth.cfgMCMinV0Pt.value > 0.f)
                                  ? mctruth.cfgMCMinV0Pt.value
                                  : pcmcuts.cfgMinPtV0.value;
        if (g.pt() < mcV0PtMin)
          continue;
        if (!g.has_daughters())
          continue;

        int posId = -1, negId = -1;
        float rTrue = -1.f;
        for (const auto& dId : g.daughtersIds()) {
          if (dId < 0) {
            continue;
          }
          const auto d = emmcParticles.iteratorAt(dId);
          if (d.pdgCode() == kElectron) {
            posId = dId;
            rTrue = std::sqrt(d.vx() * d.vx() + d.vy() * d.vy());
          } else if (d.pdgCode() == kPositron) {
            negId = dId;
          }
        }
        if (posId < 0 || negId < 0) {
          continue;
        }

        const auto mcPosE = emmcParticles.iteratorAt(posId);
        const auto mcNegE = emmcParticles.iteratorAt(negId);

        if (mctruth.cfgMCMinLegPt.value > 0.f &&
            (static_cast<float>(mcPosE.pt()) < mctruth.cfgMCMinLegPt.value ||
             static_cast<float>(mcNegE.pt()) < mctruth.cfgMCMinLegPt.value))
          continue;

        const float deTrE = static_cast<float>(mcPosE.eta() - mcNegE.eta());
        const float dpTrE = wrapPhi(static_cast<float>(mcPosE.phi() - mcNegE.phi()));
        const float legDRt = std::sqrt(deTrE * deTrE + dpTrE * dpTrE);

        const float pxG = static_cast<float>(g.px()), pyG = static_cast<float>(g.py()),
                    pzG = static_cast<float>(g.pz());
        const float magG = std::sqrt(pxG * pxG + pyG * pyG + pzG * pzG);
        float alphaTrue = 0.f;
        if (magG > kMinSigma) {
          const float ux = pxG / magG, uy = pyG / magG, uz = pzG / magG;
          const float pLpos = static_cast<float>(mcPosE.px()) * ux +
                              static_cast<float>(mcPosE.py()) * uy +
                              static_cast<float>(mcPosE.pz()) * uz;
          const float pLneg = static_cast<float>(mcNegE.px()) * ux +
                              static_cast<float>(mcNegE.py()) * uy +
                              static_cast<float>(mcNegE.pz()) * uz;
          const float sumPL = pLpos + pLneg;
          if (std::fabs(sumPL) > kMinSigma)
            alphaTrue = (pLpos - pLneg) / sumPL;
        }

        trueGammas.push_back({static_cast<int>(g.globalIndex()), posId, negId,
                              static_cast<float>(g.eta()), static_cast<float>(g.phi()),
                              static_cast<float>(g.pt()), rTrue, legDRt,
                              deTrE,
                              dpTrE,
                              alphaTrue});
      }

      {
        int nBad = 0;
        for (const auto& tg : trueGammas) {
          const auto it = gammaRecoMap.find(tg.id);
          if (it != gammaRecoMap.end() && it->second.hasV0 &&
              (mcIdsWithAnyV0Leg.count(tg.posId) == 0 || mcIdsWithAnyV0Leg.count(tg.negId) == 0))
            ++nBad;
        }
        fRegistryMC.fill(HIST("MC/TruthAO2D/hStageConsistency"), static_cast<float>(nBad));
      }

      for (const auto& tg : trueGammas) {
        const bool posFound = legIdsThisCollision.count(tg.posId) > 0;
        const bool negFound = legIdsThisCollision.count(tg.negId) > 0;
        const bool bothFound = posFound && negFound;

        for (const auto& [legId, legFound] :
             std::initializer_list<std::pair<int, bool>>{{tg.posId, posFound}, {tg.negId, negFound}}) {
          if (legId < 0)
            continue;
          if (legFound) {
            fRegistryMC.fill(HIST("MC/LegDiag/hLegDRtrue_vs_pt_legFound"), tg.pt, tg.legDRtrue);
            fRegistryMC.fill(HIST("MC/LegDiag/hLegDEta_legFound_vs_pt"), tg.pt, std::fabs(tg.legDEta));
            fRegistryMC.fill(HIST("MC/LegDiag/hLegDPhi_legFound_vs_pt"), tg.pt, std::fabs(tg.legDPhi));
          } else {
            fRegistryMC.fill(HIST("MC/LegDiag/hLegDRtrue_vs_pt_legMissing"), tg.pt, tg.legDRtrue);
            fRegistryMC.fill(HIST("MC/LegDiag/hLegDEta_legMissing_vs_pt"), tg.pt, std::fabs(tg.legDEta));
            fRegistryMC.fill(HIST("MC/LegDiag/hLegDPhi_legMissing_vs_pt"), tg.pt, std::fabs(tg.legDPhi));
          }
        }

        // ─── Armenteros-α diagnostics per photon ─────────────────────────────
        if (bothFound) {
          fRegistryMC.fill(HIST("MC/LegDiag/hAlphaTrue_legFound_vs_pt"), tg.pt, tg.alphaTrue);
        } else {
          // At least one leg missing
          fRegistryMC.fill(HIST("MC/LegDiag/hAlphaTrue_legMissing_vs_pt"), tg.pt, tg.alphaTrue);
          fRegistryMC.fill(HIST("MC/LegDiag/hAlpha_vs_legDR_legMissing"), tg.legDRtrue, tg.alphaTrue);
        }
      }

      // ─── Pair loop: efficiency  ─────────────────────────────────────
      for (size_t i = 0; i < trueGammas.size(); ++i) {
        for (size_t j = i + 1; j < trueGammas.size(); ++j) {
          const auto& g1 = trueGammas[i];
          const auto& g2 = trueGammas[j];
          const float deta = g1.eta - g2.eta;
          const float dphi = wrapPhi(g1.phi - g2.phi);
          const float px1 = g1.pt * std::cos(g1.phi), py1 = g1.pt * std::sin(g1.phi);
          const float px2 = g2.pt * std::cos(g2.phi), py2 = g2.pt * std::sin(g2.phi);
          const float kt = 0.5f * std::sqrt((px1 + px2) * (px1 + px2) + (py1 + py2) * (py1 + py2));

          if (mctruth.cfgMCMinKt > 0.f && kt < mctruth.cfgMCMinKt)
            continue;
          if (mctruth.cfgMCMaxKt > 0.f && kt > mctruth.cfgMCMaxKt)
            continue;

          const float e1 = g1.pt * std::cosh(g1.eta), e2 = g2.pt * std::cosh(g2.eta);
          const float dot = e1 * e2 - (px1 * px2 + py1 * py2 +
                                       g1.pt * std::sinh(g1.eta) * g2.pt * std::sinh(g2.eta));
          const float qinv_true = std::sqrt(std::max(0.f, 2.f * dot));

          if (mctruth.cfgMCMaxQinv > 0.f && qinv_true > mctruth.cfgMCMaxQinv)
            continue;

          auto it1 = gammaRecoMap.find(g1.id), it2 = gammaRecoMap.find(g2.id);
          const bool g1Built = (it1 != gammaRecoMap.end()) && it1->second.hasV0;
          const bool g2Built = (it2 != gammaRecoMap.end()) && it2->second.hasV0;
          const bool g1Sel = (it1 != gammaRecoMap.end()) && it1->second.passesCut;
          const bool g2Sel = (it2 != gammaRecoMap.end()) && it2->second.passesCut;

          const bool pairAll4LegsThisColl =
            legIdsThisCollision.count(g1.posId) > 0 && legIdsThisCollision.count(g1.negId) > 0 &&
            legIdsThisCollision.count(g2.posId) > 0 && legIdsThisCollision.count(g2.negId) > 0;

          fRegistryMC.fill(HIST("MC/TruthAO2D/hSparse_DEtaDPhi_qinv_truthConverted"), deta, dphi, qinv_true);
          fRegistryMC.fill(HIST("MC/TruthAO2D/hSparse_DEtaDPhi_kT_truthConverted"), deta, dphi, kt);
          fRegistryMC.fill(HIST("MC/TruthAO2D/hQinvVsKt_truthConverted"), kt, qinv_true);
          fRegistryMC.fill(HIST("MC/TruthAO2D/hDEtaDPhi_truthConverted"), deta, dphi);

          if (pairAll4LegsThisColl) {
            fRegistryMC.fill(HIST("MC/TruthAO2D/hSparse_DEtaDPhi_qinv_all4LegsThisColl"), deta, dphi, qinv_true);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hSparse_DEtaDPhi_kT_all4LegsThisColl"), deta, dphi, kt);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hQinvVsKt_all4LegsThisColl"), kt, qinv_true);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hDEtaDPhi_all4LegsThisColl"), deta, dphi);
          }
          if (g1Built && g2Built) {
            fRegistryMC.fill(HIST("MC/TruthAO2D/hSparse_DEtaDPhi_qinv_bothPhotonsBuilt"), deta, dphi, qinv_true);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hSparse_DEtaDPhi_kT_bothPhotonsBuilt"), deta, dphi, kt);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hQinvVsKt_bothPhotonsBuilt"), kt, qinv_true);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hDEtaDPhi_bothPhotonsBuilt"), deta, dphi);
          }
          if (g1Sel && g2Sel) {
            fRegistryMC.fill(HIST("MC/TruthAO2D/hSparse_DEtaDPhi_qinv_bothPhotonsSelected"), deta, dphi, qinv_true);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hSparse_DEtaDPhi_kT_bothPhotonsSelected"), deta, dphi, kt);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hQinvVsKt_bothPhotonsSelected"), kt, qinv_true);
            fRegistryMC.fill(HIST("MC/TruthAO2D/hDEtaDPhi_bothPhotonsSelected"), deta, dphi);
          }

          if (g1.rTrue >= 0.f && g2.rTrue >= 0.f) {
            fRegistryMC.fill(HIST("MC/TruthAO2D/hRconv1_vs_Rconv2_truthConverted"), g1.rTrue, g2.rTrue);
            if (g1Sel && g2Sel)
              fRegistryMC.fill(HIST("MC/TruthAO2D/hRconv1_vs_Rconv2_bothPhotonsSelected"), g1.rTrue, g2.rTrue);
          }
          const float minRconv = (g1.rTrue >= 0.f && g2.rTrue >= 0.f)
                                   ? std::min(g1.rTrue, g2.rTrue)
                                   : (g1.rTrue >= 0.f ? g1.rTrue : g2.rTrue);
          if (minRconv >= 0.f) {
            fRegistryMC.fill(HIST("MC/TruthAO2D/hMinRconv_vs_kT_truthConverted"), kt, minRconv);
            if (g1Sel && g2Sel)
              fRegistryMC.fill(HIST("MC/TruthAO2D/hMinRconv_vs_kT_bothPhotonsSelected"), kt, minRconv);
          }

          fRegistryMC.fill(HIST("MC/TruthAO2D/hStage_vs_kT"), kt, 0.f);
          if (pairAll4LegsThisColl)
            fRegistryMC.fill(HIST("MC/TruthAO2D/hStage_vs_kT"), kt, 1.f);
          if (g1Built && g2Built)
            fRegistryMC.fill(HIST("MC/TruthAO2D/hStage_vs_kT"), kt, 2.f);
          if (g1Sel && g2Sel)
            fRegistryMC.fill(HIST("MC/TruthAO2D/hStage_vs_kT"), kt, 3.f);

          const auto itCB = crossBuildMap.find(g1.id);
          if (itCB != crossBuildMap.end() && itCB->second.count(g2.id) > 0) {
            fRegistryMC.fill(HIST("MC/PairCrossBuild/hSparse_DEtaDPhi_kT"), deta, dphi, kt);
            fRegistryMC.fill(HIST("MC/PairCrossBuild/hStageOut_vs_kT"),
                             kt, static_cast<float>((g1Built ? 1 : 0) + (g2Built ? 1 : 0)));
          }

          const int nLegsThisColl =
            (legIdsThisCollision.count(g1.posId) > 0 ? 1 : 0) +
            (legIdsThisCollision.count(g1.negId) > 0 ? 1 : 0) +
            (legIdsThisCollision.count(g2.posId) > 0 ? 1 : 0) +
            (legIdsThisCollision.count(g2.negId) > 0 ? 1 : 0);
          fRegistryMC.fill(HIST("MC/LegDiag/hNLegsPair_vs_kT"), kt, static_cast<float>(nLegsThisColl));

          for (const auto& [tgRef, legId] :
               std::initializer_list<std::pair<std::reference_wrapper<const TruthGamma>, int>>{
                 {std::cref(g1), g1.posId}, {std::cref(g1), g1.negId}, {std::cref(g2), g2.posId}, {std::cref(g2), g2.negId}}) {
            if (legId < 0 || legIdsThisCollision.count(legId) > 0)
              continue;
            const auto& tg_parent = tgRef.get();
            const float legPtTrue = static_cast<float>(emmcParticles.iteratorAt(legId).pt());
            fRegistryMC.fill(HIST("MC/LegDiag/hMissingLegPt_vs_kT"), kt, legPtTrue);
            if (tg_parent.rTrue >= 0.f)
              fRegistryMC.fill(HIST("MC/LegDiag/hMissingLegRconv_vs_kT"), kt, tg_parent.rTrue);
          }
        }
      }

      // ─── Truth-level CF mixing ────────────────────────────────────────────────
      if (mctruth.cfgDoTruthMix.value) {
        for (size_t i = 0; i < trueGammas.size(); ++i) {
          for (size_t j = i + 1; j < trueGammas.size(); ++j) {
            const auto& g1 = trueGammas[i];
            const auto& g2 = trueGammas[j];
            if (!passAsymmetryCut(g1.pt, g2.pt))
              continue;
            const float deta = g1.eta - g2.eta;
            const float dphi = wrapPhi(g1.phi - g2.phi);
            const float px1 = g1.pt * std::cos(g1.phi), py1 = g1.pt * std::sin(g1.phi);
            const float px2 = g2.pt * std::cos(g2.phi), py2 = g2.pt * std::sin(g2.phi);
            const float kt = 0.5f * std::sqrt((px1 + px2) * (px1 + px2) + (py1 + py2) * (py1 + py2));
            const float e1 = g1.pt * std::cosh(g1.eta), e2 = g2.pt * std::cosh(g2.eta);
            const float dot = e1 * e2 - (px1 * px2 + py1 * py2 +
                                         g1.pt * std::sinh(g1.eta) * g2.pt * std::sinh(g2.eta));
            const float qinv_true = std::sqrt(std::max(0.f, 2.f * dot));
            fRegistryMC.fill(HIST("MC/TruthCF/hQinvVsKt_same"), kt, qinv_true);
            fRegistryMC.fill(HIST("MC/TruthCF/hDEtaDPhi_same"), deta, dphi);
          }
        }

        const float centForBin = cent[mixing.cfgCentEstimator.value];
        const std::array<float, 7> epArr = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(),
                                            collision.ep2fv0a(), collision.ep2btot(), collision.ep2bpos(),
                                            collision.ep2bneg()};
        const float ep2 = epArr[mixing.cfgEP2EstimatorForMix.value];
        const float occupancy = (mixing.cfgOccupancyEstimator.value == 1)
                                  ? static_cast<float>(collision.trackOccupancyInTimeRange())
                                  : collision.ft0cOccupancyInTimeRange();
        auto keyBin = std::make_tuple(binOf(ztxBinEdges, collision.posZ()),
                                      binOf(centBinEdges, centForBin),
                                      binOf(epBinEgdes, ep2),
                                      binOf(occBinEdges, occupancy));

        if (truthGammaPool.count(keyBin)) {
          for (const auto& poolEvent : truthGammaPool[keyBin]) {
            for (const auto& g1 : trueGammas) {
              for (const auto& g2 : poolEvent) {
                if (!passAsymmetryCut(g1.pt, g2.pt))
                  continue;
                const float deta = g1.eta - g2.eta;
                const float dphi = wrapPhi(g1.phi - g2.phi);
                const float px1 = g1.pt * std::cos(g1.phi), py1 = g1.pt * std::sin(g1.phi);
                const float px2 = g2.pt * std::cos(g2.phi), py2 = g2.pt * std::sin(g2.phi);
                const float kt = 0.5f * std::sqrt((px1 + px2) * (px1 + px2) + (py1 + py2) * (py1 + py2));
                const float e1 = g1.pt * std::cosh(g1.eta), e2 = g2.pt * std::cosh(g2.eta);
                const float dot = e1 * e2 - (px1 * px2 + py1 * py2 +
                                             g1.pt * std::sinh(g1.eta) * g2.pt * std::sinh(g2.eta));
                const float qinv_true = std::sqrt(std::max(0.f, 2.f * dot));
                fRegistryMC.fill(HIST("MC/TruthCF/hQinvVsKt_mix"), kt, qinv_true);
                fRegistryMC.fill(HIST("MC/TruthCF/hDEtaDPhi_mix"), deta, dphi);
              }
            }
          }
        }

        if (!trueGammas.empty()) {
          auto& poolBin = truthGammaPool[keyBin];
          poolBin.push_back(trueGammas);
          if (static_cast<int>(poolBin.size()) > mctruth.cfgTruthMixDepth.value)
            poolBin.pop_front();
        }
      } // end cfgDoTruthMix
    } // end collision loop
  } // end runTruthEfficiency

  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<
    std::tuple<int, int, int, int>, std::pair<int, int>, EMPair>;

  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;
  std::unordered_set<int> usedPhotonIdsPerCol;
  std::map<std::pair<int, int>, uint64_t> mapMixedEventIdToGlobalBC;

  SliceCache cache;
  Preslice<MyV0Photons> perCollisionPCM = aod::v0photonkf::pmeventId;
  PresliceUnsorted<MyMCV0Legs> perCollisionV0Legs = aod::v0leg::collisionId;
  PresliceUnsorted<aod::EMMCParticles> perMCCollisionEMMCParts = aod::emmcparticle::emmceventId;

  Filter collisionFilterCentrality =
    (centralitySelection.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < centralitySelection.cfgCentMax) ||
    (centralitySelection.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < centralitySelection.cfgCentMax) ||
    (centralitySelection.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < centralitySelection.cfgCentMax);
  Filter collisionFilterOccupancyTrack =
    eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange &&
    o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilterOccupancyFT0c =
    eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange &&
    o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;

  using FilteredMyCollisions = soa::Filtered<MyCollisions>;
  int ndf = 0;

  void processAnalysis(FilteredMyCollisions const& collisions,
                       MyV0Photons const& v0photons,
                       aod::V0Legs const& v0legs)
  {
    runPairing(collisions, v0photons, v0photons, v0legs, v0legs,
               perCollisionPCM, perCollisionPCM, fV0PhotonCut, fV0PhotonCut);
    ndf++;
  }
  PROCESS_SWITCH(Photonhbt, processAnalysis, "pairing for analysis", true);

  void processMC(FilteredMyCollisions const& collisions,
                 MyV0Photons const& v0photons,
                 MyMCV0Legs const& v0legs,
                 aod::EMMCParticles const& mcParticles,
                 aod::EMMCEvents const& mcEvents)
  {

    runPairingMC(collisions, v0photons, v0legs, mcParticles,
                 perCollisionPCM, fV0PhotonCut);
    runTruthEfficiency(collisions, v0photons, v0legs, mcParticles, mcEvents,
                       perMCCollisionEMMCParts, perCollisionV0Legs, fV0PhotonCut);

    ndf++;
  }
  PROCESS_SWITCH(Photonhbt, processMC, "MC CF + truth efficiency maps for CF correction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Photonhbt>(cfgc)};
}
