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
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
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
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h>
#include <Math/Vector4Dfwd.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

/// Single-photon track-type combo.
enum class V0Combo : int {
  Inclusive = 0,
  ItstpcItstpc = 1,   ///< both legs ITS+TPC
  ItstpcTpconly = 2,  ///< one ITS+TPC leg, one TPC-only
  TpconlyTpconly = 3, ///< both legs TPC-only
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

using MyCollisions = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000,
                               aod::EMEventsCent_000, aod::EMEventsQvec_001>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsPhiVPsi>;
using MyV0Photon = MyV0Photons::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

// ─── MC truth classification types ────────────────────────────────────────────

/// Per-photon MC truth information built from the two V0 legs.
struct PhotonMCInfo {
  bool hasMC = false;        // both legs have a valid MC label
  bool sameMother = false;   // both legs share the same MC mother
  bool isTruePhoton = false; // mother PDG == 22

  int mcPosId = -1;  // MC particle index of the positive leg
  int mcNegId = -1;  // MC particle index of the negative leg
  int motherId = -1; //  MC particle index of the common mother
  int motherPdg = 0;

  bool isPhysicalPrimary = false;
};

/// Classification of a photon pair at the MC-truth level.
enum class PairTruthType : uint8_t {
  Unknown = 0,
  TrueTrueDistinct,   // both photons are true, from different MC photons
  TrueTrueSamePhoton, // both photons are true, same MC photon (clone/split)
  SharedMcLeg,        // different reco tracks but same MC-level leg
  TrueFake,           // one photon is true, one is fake
  FakeFake,           // both photons are fake
  Pi0Daughters,       // both photons come from the same MC pi0
};

struct photonhbt {

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

  // ─── Configurables: histogram axis bins ───────────────────────────────────

  // HBT physics
  ConfigurableAxis confQBins{"confQBins", {60, 0, +0.3f}, "q bins for output histograms"};
  ConfigurableAxis confKtBins{"confKtBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}, "kT bins"};

  // Single-photon QA
  ConfigurableAxis confPtBins{"confPtBins", {100, 0.f, 2.f}, "pT bins (GeV/c)"};
  ConfigurableAxis confEtaBins{"confEtaBins", {80, -0.8f, 0.8f}, "eta bins"};
  ConfigurableAxis confPhiBins{"confPhiBins", {90, 0.f, o2::constants::math::TwoPI}, "phi bins (rad) — O2 track phi is in [0, 2pi]"};

  // Pair angular
  ConfigurableAxis confDeltaEtaBins{"confDeltaEtaBins", {100, -0.9f, +0.9f}, "Delta-eta bins"};
  ConfigurableAxis confDeltaPhiBins{"confDeltaPhiBins", {100, -o2::constants::math::PI, o2::constants::math::PI}, "Delta-phi bins (rad)"};
  ConfigurableAxis confEllipseValBins{"confEllipseValBins", {200, 0.f, 10.f}, "ellipse value bins"};
  ConfigurableAxis confCosThetaBins{"confCosThetaBins", {100, 0.f, 1.f}, "cos(theta*) bins"};
  ConfigurableAxis confOpeningAngleBins{"confOpeningAngleBins", {100, 0.f, o2::constants::math::PI}, "opening angle bins (rad)"};

  // Pair geometry
  ConfigurableAxis confRBins{"confRBins", {100, 0.f, 100.f}, "conversion radius bins (cm)"};
  ConfigurableAxis confDeltaRBins{"confDeltaRBins", {120, 0.f, 30.f}, "|R1-R2| bins (cm)"};
  ConfigurableAxis confDeltaR3DBins{"confDeltaR3DBins", {100, 0.f, 100.f}, "3D distance between conversion points (cm)"};
  ConfigurableAxis confDeltaRxyBins{"confDeltaRxyBins", {100, 0.f, 100.f}, "xy distance between conversion points (cm)"};
  ConfigurableAxis confZConvBins{"confZConvBins", {200, -100.f, 100.f}, "conversion z (cm)"};
  ConfigurableAxis confDeltaZBins{"confDeltaZBins", {200, -100.f, 100.f}, "#Deltaz bins (cm)"}; ///< FIX: was missing

  // Event QA
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
  const AxisSpec axisDeltaZ{confDeltaZBins, "#Delta z (cm)"}; ///< FIX: was missing
  const AxisSpec axisOccupancy{confOccupancyQA, "occupancy"};
  const AxisSpec axisCentQA{confCentQABins, "centrality (%)"};

  // ─── Configurables: QA flags ───────────────────────────────────────────────

  struct : ConfigurableGroup {
    std::string prefix = "qaflags_group";
    Configurable<bool> doPairQa{"doPairQa", true, "fill pair QA histograms at each cut step"};
    Configurable<bool> doSinglePhotonQa{"doSinglePhotonQa", true, "fill single-photon QA histograms (pT, eta, phi)"};
    Configurable<float> cfgMaxQinvForQA{"cfgMaxQinvForQA", 0.1f,
                                        "fill per-step pair QA histograms (hDeltaEta, hDeltaPhi, THnSparses, ...) "
                                        "only when q_inv < this value (GeV/c). "
                                        "Set to the HBT signal region, typically 0.1. "
                                        "Set <= 0 to disable the gate. The CF is always filled regardless."};
    Configurable<float> cfgMaxQinvForFullRange{"cfgMaxQinvForFullRange", 0.3f,
                                               "fill full-range histograms (hDeltaRVsQinv, hSparseDeltaRDeltaZQinv, ...) "
                                               "only when q_inv < this value (GeV/c). "
                                               "Should match the upper edge of confQBins (default 0.3) — "
                                               "fills beyond the axis range only go into the overflow bin. "
                                               "Set <= 0 to disable the gate."};
  } qaflags;

  // ─── Configurables: HBT kind ───────────────────────────────────────────────

  Configurable<bool> cfgDo3D{"cfgDo3D", false, "enable 3D (qout,qside,qlong) analysis"};
  Configurable<bool> cfgDo2D{"cfgDo2D", false, "enable 2D (qout,qinv) projection (requires cfgDo3D)"};
  Configurable<bool> cfgUseLCMS{"cfgUseLCMS", true, "measure 1D relative momentum in LCMS"};

  // ─── Configurables: events ─────────────────────────────────────────────────

  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity"};

  // ─── Configurables: mixed event ────────────────────────────────────────────

  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  Configurable<uint64_t> ndiffBCMix{"ndiffBCMix", 594, "difference in global BC required for mixed events"};
  Configurable<int> cfgEP2EstimatorForMix{"cfgEP2EstimatorForMix", 3, "FT0M:0, FT0A:1, FT0C:2, FV0A:3, BTot:4, BPos:5, BNeg:6"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};

  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {16, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}, "Mixing bins - EP angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  // ─── Configurables: pair cuts ──────────────────────────────────────────────

  struct : ConfigurableGroup {
    std::string prefix = "ggpaircut_group";
    // dr/cosOA cut
    Configurable<float> cfgMinDRCosOA{"cfgMinDRCosOA", -1.f, "min. dr/cosOA; <0 = disabled"};
    // R/Z geometry cuts
    Configurable<bool> cfgDoRCut{"cfgDoRCut", false, "apply |R1-R2| > cfgMinDeltaR cut"};
    Configurable<float> cfgMinDeltaR{"cfgMinDeltaR", 0.f, "minimum |R1-R2| (cm)"};
    Configurable<bool> cfgDoZCut{"cfgDoZCut", false, "apply |DeltaZ| > cfgMinDeltaZ cut"};
    Configurable<float> cfgMinDeltaZ{"cfgMinDeltaZ", 0.f, "minimum |DeltaZ| (cm)"};
    // Ellipse cut in (DeltaEta, DeltaPhi)
    Configurable<bool> cfgDoEllipseCut{"cfgDoEllipseCut", false, "reject pairs inside ellipse in DeltaEta-DeltaPhi"};
    Configurable<float> cfgEllipseSigEta{"cfgEllipseSigEta", 0.02f, "sigma_eta for ellipse cut"};
    Configurable<float> cfgEllipseSigPhi{"cfgEllipseSigPhi", 0.02f, "sigma_phi for ellipse cut"};
    Configurable<float> cfgEllipseR2{"cfgEllipseR2", 1.0f, "R^2 threshold: reject if ellipse value < R^2"};
  } ggpaircuts;

  // ─── Event cut ─────────────────────────────────────────────────────────────

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

  // ─── PCM cut ───────────────────────────────────────────────────────────────

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

  ~photonhbt()
  {
    delete emh1;
    emh1 = nullptr;
    delete emh2;
    emh2 = nullptr;
    mapMixedEventIdToGlobalBC.clear();
    usedPhotonIdsPerCol.clear();
  }

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::mt19937 engine;
  std::uniform_int_distribution<int> dist01;
  int mRunNumber{0};

  std::vector<float> ztxBinEdges;
  std::vector<float> centBinEdges;
  std::vector<float> epBinEgdes;
  std::vector<float> occBinEdges;

  // ─── Pair-cut helpers ──────────────────────────────────────────────────────

  inline bool isInsideEllipse(float deta, float dphi) const
  {
    if (!ggpaircuts.cfgDoEllipseCut.value) // .value needed: operator T() is non-const in O2
      return false;
    const float sE = ggpaircuts.cfgEllipseSigEta.value;
    const float sP = ggpaircuts.cfgEllipseSigPhi.value;
    if (sE < 1e-9f || sP < 1e-9f)
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

  static inline float computeCosTheta(const ROOT::Math::PtEtaPhiMVector& v1,
                                      const ROOT::Math::PtEtaPhiMVector& v2)
  {
    ROOT::Math::PxPyPzEVector p1(v1), p2(v2);
    ROOT::Math::PxPyPzEVector pair = p1 + p2;
    ROOT::Math::Boost boost(-pair.BoostToCM());
    ROOT::Math::PxPyPzEVector p1cm = boost(p1);
    ROOT::Math::XYZVector pairDir(pair.Px(), pair.Py(), pair.Pz());
    ROOT::Math::XYZVector p1cmDir(p1cm.Px(), p1cm.Py(), p1cm.Pz());
    if (pairDir.R() < 1e-9 || p1cmDir.R() < 1e-9)
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
    return clampBin(b, static_cast<int>(edges.size()) - 2);
  }

  /// ev_id  : 0 = same-event, 1 = mixed-event
  /// step_id: 0 = Before, 1 = AfterDRCosOA, 2 = AfterRZ, 3 = AfterEllipse
  template <int ev_id, int step_id>
  static constexpr const char* qaPrefix()
  {
    if constexpr (ev_id == 0) {
      if constexpr (step_id == 0)
        return "Pair/same/QA/Before/";
      if constexpr (step_id == 1)
        return "Pair/same/QA/AfterDRCosOA/";
      if constexpr (step_id == 2)
        return "Pair/same/QA/AfterRZ/";
      return "Pair/same/QA/AfterEllipse/";
    } else {
      if constexpr (step_id == 0)
        return "Pair/mix/QA/Before/";
      if constexpr (step_id == 1)
        return "Pair/mix/QA/AfterDRCosOA/";
      if constexpr (step_id == 2)
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

    parseBins(ConfVtxBins, ztxBinEdges);
    parseBins(ConfCentBins, centBinEdges);
    parseBins(ConfEPBins, epBinEgdes);
    parseBins(ConfOccupancyBins, occBinEdges);

    emh1 = new MyEMH(ndepth);
    emh2 = new MyEMH(ndepth);

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
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber())
      return;
    mRunNumber = collision.runNumber();
  }

  // ─── PairQAObservables ─────────────────────────────────────────────────────
  /// Plain data struct holding all observables computed from a photon pair.
  /// Defined early so all histogram-booking and fill functions can use it.

  struct PairQAObservables {
    // conversion-point coordinates
    float x1 = 0.f, y1 = 0.f, z1 = 0.f;
    float x2 = 0.f, y2 = 0.f, z2 = 0.f;
    // conversion-point radii and distances
    float r1 = 0.f, r2 = 0.f;
    float dx = 0.f, dy = 0.f, dz = 0.f;
    float deltaR = 0.f;   ///< |R1-R2|
    float deltaZ = 0.f;   ///< z1-z2
    float deltaRxy = 0.f; ///< sqrt(dx^2+dy^2)
    float deltaR3D = 0.f; ///< sqrt(dx^2+dy^2+dz^2)
    // opening angle of conversion-point vectors
    float opa = 0.f;
    float cosOA = 0.f;
    float drOverCosOA = 0.f;
    // photon four-vectors and pair kinematics
    ROOT::Math::PtEtaPhiMVector v1;
    ROOT::Math::PtEtaPhiMVector v2;
    ROOT::Math::PtEtaPhiMVector k12;
    float deta = 0.f, dphi = 0.f;
    float pairEta = 0.f, pairPhi = 0.f;
    float kt = 0.f, qinv = 0.f;
    float cosTheta = 0.f;
    float openingAngle = 0.f;
    // validity flag
    bool valid = true;
  };

  void addSinglePhotonQAHistogramsForStep(const std::string& path)
  {
    fRegistry.add((path + "hPt").c_str(), "p_{T};p_{T} (GeV/c);counts", kTH1D, {axisPt}, true);
    fRegistry.add((path + "hEta").c_str(), "#eta;#eta;counts", kTH1D, {axisEta}, true);
    fRegistry.add((path + "hPhi").c_str(), "#phi;#phi (rad);counts", kTH1D, {axisPhi}, true);
    fRegistry.add((path + "hEtaVsPhi").c_str(), "acceptance;#phi (rad);#eta", kTH2D, {axisPhi, axisEta}, true);
    fRegistry.add((path + "hR").c_str(), "R_{conv};R_{conv} (cm);counts", kTH1D, {axisR}, true);
    fRegistry.add((path + "hZConv").c_str(), "z_{conv};z_{conv} (cm);counts", kTH1D, {axisZConv}, true);
    fRegistry.add((path + "hRVsZConv").c_str(), "R_{conv} vs z_{conv};z_{conv} (cm);R_{conv} (cm)", kTH2D, {axisZConv, axisR}, true);
  }

  void addFullRangeHistograms(const std::string& path)
  {
    fRegistry.add((path + "hDeltaRVsQinv").c_str(), "|R_{1}-R_{2}| vs q_{inv} (full range);q_{inv} (GeV/c);|R_{1}-R_{2}| (cm)", kTH2D, {axisQinv, axisDeltaR}, true);
    fRegistry.add((path + "hDeltaZVsQinv").c_str(), "#Delta z vs q_{inv} (full range);q_{inv} (GeV/c);#Delta z (cm)", kTH2D, {axisQinv, axisDeltaZ}, true);
    fRegistry.add((path + "hDeltaR3DVsQinv").c_str(), "#Delta r_{3D} vs q_{inv} (full range);q_{inv} (GeV/c);#Delta r_{3D} (cm)", kTH2D, {axisQinv, axisDeltaR3D}, true);
    fRegistry.add((path + "hQinvVsCent").c_str(), "q_{inv} vs centrality (full range);centrality (%);q_{inv} (GeV/c)", kTH2D, {axisCentQA, axisQinv}, true);
    fRegistry.add((path + "hQinvVsOccupancy").c_str(), "q_{inv} vs occupancy (full range);occupancy;q_{inv} (GeV/c)", kTH2D, {axisOccupancy, axisQinv}, true);
    fRegistry.add((path + "hSparseDeltaRDeltaZQinv").c_str(), "|R_{1}-R_{2}|,#Delta z,q_{inv} (full range)", kTHnSparseD, {axisDeltaR, axisDeltaZ, axisQinv}, true);
    fRegistry.add((path + "hDeltaRCosOAVsQinv").c_str(), "#Delta r/cos(#theta_{op}/2) vs q_{inv};q_{inv} (GeV/c);#Delta r/cos(#theta_{op}/2) (cm)", kTH2D, {axisQinv, {100, 0, 100}}, true);
  }

  /// ev_id : 0 = same-event, 1 = mixed-event
  template <int ev_id>
  inline void fillFullRangeQA(PairQAObservables const& obs, float cent, float occupancy)
  {
    constexpr auto base = fullRangePrefix<ev_id>();
    fRegistry.fill(HIST(base) + HIST("hDeltaRVsQinv"), obs.qinv, obs.deltaR);
    fRegistry.fill(HIST(base) + HIST("hDeltaZVsQinv"), obs.qinv, obs.deltaZ);
    fRegistry.fill(HIST(base) + HIST("hDeltaR3DVsQinv"), obs.qinv, obs.deltaR3D);
    fRegistry.fill(HIST(base) + HIST("hQinvVsCent"), cent, obs.qinv);
    fRegistry.fill(HIST(base) + HIST("hQinvVsOccupancy"), occupancy, obs.qinv);
    fRegistry.fill(HIST(base) + HIST("hSparseDeltaRDeltaZQinv"),
                   obs.deltaR, obs.deltaZ, obs.qinv);
  }

  template <int ev_id>
  inline void fillFullRangeDeltaRCosOA(float qinv, float drOverCosOA)
  {
    constexpr auto base = fullRangePrefix<ev_id>();
    fRegistry.fill(HIST(base) + HIST("hDeltaRCosOAVsQinv"), qinv, drOverCosOA);
  }

  void addQAHistogramsForStep(const std::string& path)
  {
    // ── 1D: photon kinematics ────────────────────────────────────────────────
    fRegistry.add((path + "hPairEta").c_str(), "pair #eta;#eta_{pair};counts", kTH1D, {axisEta}, true);
    fRegistry.add((path + "hPairPhi").c_str(), "pair #phi;#phi_{pair} (rad);counts", kTH1D, {axisPhi}, true);
    fRegistry.add((path + "hPairKt").c_str(), "pair k_{T};k_{T} (GeV/c);counts", kTH1D, {axisKt}, true);
    fRegistry.add((path + "hQinv").c_str(), "q_{inv};q_{inv} (GeV/c);counts", kTH1D, {axisQinv}, true);

    // ── 1D: angular ─────────────────────────────────────────────────────────
    fRegistry.add((path + "hDeltaEta").c_str(), "#Delta#eta;#Delta#eta;counts", kTH1D, {axisDeltaEta}, true);
    fRegistry.add((path + "hDeltaPhi").c_str(), "#Delta#phi;#Delta#phi (rad);counts", kTH1D, {axisDeltaPhi}, true);
    fRegistry.add((path + "hCosTheta").c_str(), "cos(#theta*) in pair rest frame;cos(#theta*);counts", kTH1D, {axisCosTheta}, true);
    fRegistry.add((path + "hOpeningAngle").c_str(), "Opening angle;#alpha (rad);counts", kTH1D, {axisOpeningAngle}, true);
    fRegistry.add((path + "hEllipseVal").c_str(), "(#Delta#eta/#sigma_{#eta})^{2}+(#Delta#phi/#sigma_{#phi})^{2};value;counts", kTH1D, {axisEllipseVal}, true);

    // ── 1D: geometry ────────────────────────────────────────────────────────
    fRegistry.add((path + "hR1").c_str(), "R_{conv,1};R_{1} (cm);counts", kTH1D, {axisR}, true);
    fRegistry.add((path + "hR2").c_str(), "R_{conv,2};R_{2} (cm);counts", kTH1D, {axisR}, true);
    fRegistry.add((path + "hDeltaR").c_str(), "|R_{1}-R_{2}|;|R_{1}-R_{2}| (cm);counts", kTH1D, {axisDeltaR}, true);
    fRegistry.add((path + "hDeltaZ").c_str(), "#Delta z;#Delta z (cm);counts", kTH1D, {axisDeltaZ}, true);
    fRegistry.add((path + "hDeltaRxy").c_str(), "#Delta r_{xy};#Delta r_{xy} (cm);counts", kTH1D, {axisDeltaRxy}, true);
    fRegistry.add((path + "hDeltaR3D").c_str(), "|#vec{r}_{1}-#vec{r}_{2}|;#Delta r_{3D} (cm);counts", kTH1D, {axisDeltaR3D}, true);

    // ── 1D: event-level ─────────────────────────────────────────────────────
    fRegistry.add((path + "hCent").c_str(), "centrality;centrality (%);counts", kTH1D, {axisCentQA}, true);
    fRegistry.add((path + "hOccupancy").c_str(), "occupancy;occupancy;counts", kTH1D, {axisOccupancy}, true);

    // ── 2D: angular ─────────────────────────────────────────────────────────
    fRegistry.add((path + "hDEtaDPhi").c_str(), "#Delta#eta vs #Delta#phi;#Delta#eta;#Delta#phi (rad)", kTH2D, {axisDeltaEta, axisDeltaPhi}, true);
    fRegistry.add((path + "hDeltaEtaVsPairEta").c_str(), "#Delta#eta vs #LT#eta#GT_{pair};#LT#eta#GT_{pair};#Delta#eta", kTH2D, {axisEta, axisDeltaEta}, true);

    // ── 2D: geometry ────────────────────────────────────────────────────────
    fRegistry.add((path + "hR1VsR2").c_str(), "R_{1} vs R_{2};R_{1} (cm);R_{2} (cm)", kTH2D, {axisR, axisR}, true);
    fRegistry.add((path + "hDeltaRVsDeltaZ").c_str(), "|R_{1}-R_{2}| vs #Delta z;|R_{1}-R_{2}| (cm);#Delta z (cm)", kTH2D, {axisDeltaR, axisDeltaZ}, true);

    // ── 2D: geometry vs kT ──────────────────────────────────────────────────
    // Note: hDeltaRVsQinv, hDeltaZVsQinv live in FullRange/ (always filled, full q range)
    fRegistry.add((path + "hDeltaRVsKt").c_str(), "|R_{1}-R_{2}| vs k_{T};k_{T} (GeV/c);|R_{1}-R_{2}| (cm)", kTH2D, {axisKt, axisDeltaR}, true);
    fRegistry.add((path + "hDeltaZVsKt").c_str(), "#Delta z vs k_{T};k_{T} (GeV/c);#Delta z (cm)", kTH2D, {axisKt, axisDeltaZ}, true);

    // ── 2D: angular vs geometry ─────────────────────────────────────────────
    fRegistry.add((path + "hDeltaPhiVsDeltaR").c_str(), "#Delta#phi vs |R_{1}-R_{2}|;|R_{1}-R_{2}| (cm);#Delta#phi (rad)", kTH2D, {axisDeltaR, axisDeltaPhi}, true);
    fRegistry.add((path + "hDeltaEtaVsDeltaR").c_str(), "#Delta#eta vs |R_{1}-R_{2}|;|R_{1}-R_{2}| (cm);#Delta#eta", kTH2D, {axisDeltaR, axisDeltaEta}, true);
    fRegistry.add((path + "hDeltaPhiVsDeltaZ").c_str(), "#Delta#phi vs #Delta z;#Delta z (cm);#Delta#phi (rad)", kTH2D, {axisDeltaZ, axisDeltaPhi}, true);
    fRegistry.add((path + "hDeltaEtaVsDeltaZ").c_str(), "#Delta#eta vs #Delta z;#Delta z (cm);#Delta#eta", kTH2D, {axisDeltaZ, axisDeltaEta}, true);

    // ── 2D: vs event properties ─────────────────────────────────────────────
    fRegistry.add((path + "hDeltaRVsCent").c_str(), "|R_{1}-R_{2}| vs centrality;centrality (%);|R_{1}-R_{2}| (cm)", kTH2D, {axisCentQA, axisDeltaR}, true);
    fRegistry.add((path + "hDeltaRVsOccupancy").c_str(), "|R_{1}-R_{2}| vs occupancy;occupancy;|R_{1}-R_{2}| (cm)", kTH2D, {axisOccupancy, axisDeltaR}, true);

    fRegistry.add((path + "hSparseDEtaDPhiCent").c_str(),
                  "#Delta#eta,#Delta#phi,centrality",
                  kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisCentQA}, true);
    fRegistry.add((path + "hSparseDEtaDPhiOcc").c_str(),
                  "#Delta#eta,#Delta#phi,occupancy",
                  kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisOccupancy}, true);
    fRegistry.add((path + "hSparseDEtaDPhiKt").c_str(),
                  "#Delta#eta,#Delta#phi,k_{T}",
                  kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisKt}, true);
    fRegistry.add((path + "hSparseDeltaRDeltaZKt").c_str(),
                  "|R_{1}-R_{2}|,#Delta z,k_{T}",
                  kTHnSparseD, {axisDeltaR, axisDeltaZ, axisKt}, true);
  }

  void addhistograms()
  {
    static constexpr std::string_view det[6] = {"FT0M", "FT0A", "FT0C", "BTot", "BPos", "BNeg"};
    fRegistry.add("Event/before/hEP2_CentFT0C_forMix",
                  Form("2nd harmonics EP for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", det[cfgEP2EstimatorForMix].data()),
                  kTH2D, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);
    fRegistry.add("Event/after/hEP2_CentFT0C_forMix",
                  Form("2nd harmonics EP for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", det[cfgEP2EstimatorForMix].data()),
                  kTH2D, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);

    addSinglePhotonQAHistogramsForStep("SinglePhoton/Before/");
    addSinglePhotonQAHistogramsForStep("SinglePhoton/AfterDRCosOA/");
    addSinglePhotonQAHistogramsForStep("SinglePhoton/AfterRZ/");
    addSinglePhotonQAHistogramsForStep("SinglePhoton/AfterEllipse/");

    // ── HBT correlation functions ─────────────────────────────────────────────
    if (cfgDo3D) {
      fRegistry.add("Pair/same/CF_3D", "diphoton correlation 3D LCMS",
                    kTHnSparseD, {axisQout, axisQside, axisQlong, axisKt}, true);
      if (cfgDo2D) {
        fRegistry.add("Pair/same/CF_2D", "diphoton correlation 2D (qout,qinv)",
                      kTHnSparseD, {axisQout, axisQinv, axisKt}, true);
      }
    } else {
      if (cfgUseLCMS) {
        fRegistry.add("Pair/same/CF_1D", "diphoton correlation 1D LCMS", kTH2D, {axisQabsLcms, axisKt}, true);
      } else {
        fRegistry.add("Pair/same/CF_1D", "diphoton correlation 1D (qinv)", kTH2D, {axisQinv, axisKt}, true);
      }
    }

    fRegistry.add("Pair/same/hDeltaRCosOA",
                  "distance between 2 conversion points / cos(#theta_{op}/2);#Delta r / cos(#theta_{op}/2) (cm);counts",
                  kTH1D, {{100, 0, 100}}, true);

    // ── QA steps (same-event; mix-event cloned below) ─────────────────────────
    addQAHistogramsForStep("Pair/same/QA/Before/");
    addQAHistogramsForStep("Pair/same/QA/AfterDRCosOA/");
    addQAHistogramsForStep("Pair/same/QA/AfterRZ/");
    addQAHistogramsForStep("Pair/same/QA/AfterEllipse/");

    // ── MC truth histograms (same-event; mix-event cloned below) ─────────────
    addMCHistograms();

    // ── Full-range histograms: always filled, qinv as axis ───────────────────
    addFullRangeHistograms("Pair/same/FullRange/");

    // Clone all Pair/same/ histograms to Pair/mix/
    fRegistry.addClone("Pair/same/", "Pair/mix/");
  }

  // ─── DefineEMEventCut ──────────────────────────────────────────────────────

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

  // ─── DefinePCMCut ──────────────────────────────────────────────────────────

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

  /// step_id: 0 = Before, 1 = AfterDRCosOA, 2 = AfterRZ, 3 = AfterEllipse
  template <int step_id>
  static constexpr const char* singlePhotonQAPrefix()
  {
    if constexpr (step_id == 0)
      return "SinglePhoton/Before/";
    if constexpr (step_id == 1)
      return "SinglePhoton/AfterDRCosOA/";
    if constexpr (step_id == 2)
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
    fRegistry.fill(HIST(base) + HIST("hPt"), g.pt());
    fRegistry.fill(HIST(base) + HIST("hEta"), g.eta());
    fRegistry.fill(HIST(base) + HIST("hPhi"), g.phi());
    fRegistry.fill(HIST(base) + HIST("hEtaVsPhi"), g.phi(), g.eta());
    fRegistry.fill(HIST(base) + HIST("hR"), r);
    fRegistry.fill(HIST(base) + HIST("hZConv"), g.vz());
    fRegistry.fill(HIST(base) + HIST("hRVsZConv"), g.vz(), r);
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

    if (cfgDo3D) {
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("CF_3D"),
                     std::fabs(qout_lcms), std::fabs(qside_lcms), std::fabs(qlong_lcms), kt, weight);
      if (cfgDo2D) {
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("CF_2D"),
                       std::fabs(qout_lcms), std::fabs(qinv), kt, weight);
      }
    } else {
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("CF_1D"),
                     cfgUseLCMS ? qabs_lcms : qinv, kt, weight);
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

    if (cfgDo3D) {
      fRegistry.fill(HIST(mcDir) + HIST("CF_3D"),
                     std::fabs(qout_lcms), std::fabs(qside_lcms), std::fabs(qlong_lcms), kt, weight);
      if (cfgDo2D) {
        fRegistry.fill(HIST(mcDir) + HIST("CF_2D"),
                       std::fabs(qout_lcms), std::fabs(qinv), kt, weight);
      }
    } else {
      fRegistry.fill(HIST(mcDir) + HIST("CF_1D"),
                     cfgUseLCMS ? qabs_lcms : qinv, kt, weight);
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

    ROOT::Math::XYZVector cp1(o.x1, o.y1, o.z1);
    ROOT::Math::XYZVector cp2(o.x2, o.y2, o.z2);
    const float mag1 = std::sqrt(cp1.Mag2());
    const float mag2 = std::sqrt(cp2.Mag2());
    if (mag1 < 1e-12f || mag2 < 1e-12f) {
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
    o.drOverCosOA = (std::fabs(o.cosOA) < 1e-12f) ? 1e12f : (o.deltaR3D / o.cosOA);

    o.v1 = ROOT::Math::PtEtaPhiMVector(g1.pt(), g1.eta(), g1.phi(), 0.f);
    o.v2 = ROOT::Math::PtEtaPhiMVector(g2.pt(), g2.eta(), g2.phi(), 0.f);
    o.k12 = 0.5f * (o.v1 + o.v2);

    o.deta = g1.eta() - g2.eta();
    o.dphi = RecoDecay::constrainAngle(g1.phi() - g2.phi(), -o2::constants::math::PI); // dphi in [-pi, pi]
    o.pairEta = 0.5f * (g1.eta() + g2.eta());
    o.pairPhi = RecoDecay::constrainAngle(o.k12.Phi(), 0.f); // pair phi in [0, 2pi] — matches axisPhi
    o.kt = o.k12.Pt();
    o.qinv = std::fabs((o.v1 - o.v2).M());
    o.cosTheta = std::fabs(computeCosTheta(o.v1, o.v2));
    o.openingAngle = o.opa;

    return o;
  }

  template <int ev_id, int step_id>
  inline void fillPairQAStep(PairQAObservables const& o, float cent, float occupancy)
  {
    if (!qaflags.doPairQa)
      return;

    constexpr auto base = qaPrefix<ev_id, step_id>();

    // 1D: kinematics
    fRegistry.fill(HIST(base) + HIST("hPairEta"), o.pairEta);
    fRegistry.fill(HIST(base) + HIST("hPairPhi"), o.pairPhi);
    fRegistry.fill(HIST(base) + HIST("hPairKt"), o.kt);
    fRegistry.fill(HIST(base) + HIST("hQinv"), o.qinv);

    // 1D: angular
    fRegistry.fill(HIST(base) + HIST("hDeltaEta"), o.deta);
    fRegistry.fill(HIST(base) + HIST("hDeltaPhi"), o.dphi);
    fRegistry.fill(HIST(base) + HIST("hCosTheta"), o.cosTheta);
    fRegistry.fill(HIST(base) + HIST("hOpeningAngle"), o.openingAngle);

    // 1D: geometry
    fRegistry.fill(HIST(base) + HIST("hR1"), o.r1);
    fRegistry.fill(HIST(base) + HIST("hR2"), o.r2);
    fRegistry.fill(HIST(base) + HIST("hDeltaR"), o.deltaR);
    fRegistry.fill(HIST(base) + HIST("hDeltaZ"), o.deltaZ);
    fRegistry.fill(HIST(base) + HIST("hDeltaRxy"), o.deltaRxy);
    fRegistry.fill(HIST(base) + HIST("hDeltaR3D"), o.deltaR3D);

    // 1D: event
    fRegistry.fill(HIST(base) + HIST("hCent"), cent);
    fRegistry.fill(HIST(base) + HIST("hOccupancy"), occupancy);

    // 1D: ellipse value (diagnostic, conditional on cut being configured)
    const float sE = ggpaircuts.cfgEllipseSigEta.value;
    const float sP = ggpaircuts.cfgEllipseSigPhi.value;
    if (sE > 1e-9f && sP > 1e-9f) {
      const float ellipseVal = (o.deta / sE) * (o.deta / sE) + (o.dphi / sP) * (o.dphi / sP);
      fRegistry.fill(HIST(base) + HIST("hEllipseVal"), ellipseVal);
    }

    // 2D: angular
    fRegistry.fill(HIST(base) + HIST("hDEtaDPhi"), o.deta, o.dphi);
    fRegistry.fill(HIST(base) + HIST("hDeltaEtaVsPairEta"), o.pairEta, o.deta);

    // 2D: geometry
    fRegistry.fill(HIST(base) + HIST("hR1VsR2"), o.r1, o.r2);
    fRegistry.fill(HIST(base) + HIST("hDeltaRVsDeltaZ"), o.deltaR, o.deltaZ);

    // 2D: geometry vs kT (qinv variants live in FullRange/)
    fRegistry.fill(HIST(base) + HIST("hDeltaRVsKt"), o.kt, o.deltaR);
    fRegistry.fill(HIST(base) + HIST("hDeltaZVsKt"), o.kt, o.deltaZ);

    // 2D: angular vs geometry
    fRegistry.fill(HIST(base) + HIST("hDeltaPhiVsDeltaR"), o.deltaR, o.dphi);
    fRegistry.fill(HIST(base) + HIST("hDeltaEtaVsDeltaR"), o.deltaR, o.deta);
    fRegistry.fill(HIST(base) + HIST("hDeltaPhiVsDeltaZ"), o.deltaZ, o.dphi);
    fRegistry.fill(HIST(base) + HIST("hDeltaEtaVsDeltaZ"), o.deltaZ, o.deta);

    // 2D: vs event properties (qinv variants live in FullRange/)
    fRegistry.fill(HIST(base) + HIST("hDeltaRVsCent"), cent, o.deltaR);
    fRegistry.fill(HIST(base) + HIST("hDeltaRVsOccupancy"), occupancy, o.deltaR);

    // THnSparse (hSparseDeltaRDeltaZQinv lives in FullRange/)
    fRegistry.fill(HIST(base) + HIST("hSparseDEtaDPhiCent"), o.deta, o.dphi, cent);
    fRegistry.fill(HIST(base) + HIST("hSparseDEtaDPhiOcc"), o.deta, o.dphi, occupancy);
    fRegistry.fill(HIST(base) + HIST("hSparseDEtaDPhiKt"), o.deta, o.dphi, o.kt);
    fRegistry.fill(HIST(base) + HIST("hSparseDeltaRDeltaZKt"), o.deltaR, o.deltaZ, o.kt);
  }

  template <typename TPhoton, typename TLegs, typename TMCParticles>
  static PhotonMCInfo buildPhotonMCInfo(TPhoton const& g,
                                        TMCParticles const& mcParticles)
  {
    PhotonMCInfo info{};

    const auto pos = g.template posTrack_as<TLegs>();
    const auto neg = g.template negTrack_as<TLegs>();

    // PWGEM uses emmcparticle, not the standard mcParticle accessor
    if (!pos.has_emmcparticle() || !neg.has_emmcparticle())
      return info;

    info.hasMC = true;
    info.mcPosId = pos.emmcparticleId();
    info.mcNegId = neg.emmcparticleId();

    const auto mcPos = pos.template emmcparticle_as<TMCParticles>();
    const auto mcNeg = neg.template emmcparticle_as<TMCParticles>();

    if (!mcPos.has_mothers() || !mcNeg.has_mothers())
      return info;

    const int mothIdPos = mcPos.mothersIds()[0];
    const int mothIdNeg = mcNeg.mothersIds()[0];
    if (mothIdPos != mothIdNeg)
      return info;

    info.sameMother = true;
    info.motherId = mothIdPos;

    const auto mother = mcParticles.iteratorAt(mothIdPos);
    info.motherPdg = mother.pdgCode();
    info.isTruePhoton = (info.motherPdg == 22);
    info.isPhysicalPrimary = mother.isPhysicalPrimary();

    return info;
  }

  static PairTruthType classifyPairTruth(PhotonMCInfo const& m1,
                                         PhotonMCInfo const& m2)
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

    // Both are true photons — same or different MC photon?
    if (m1.motherId >= 0 && m1.motherId == m2.motherId)
      return PairTruthType::TrueTrueSamePhoton;

    return PairTruthType::TrueTrueDistinct;
  }

  template <typename TMCParticles>
  static bool isPi0DaughterPair(PhotonMCInfo const& m1,
                                PhotonMCInfo const& m2,
                                TMCParticles const& mcParticles)
  {
    if (!m1.isTruePhoton || !m2.isTruePhoton)
      return false;
    if (m1.motherId < 0 || m2.motherId < 0)
      return false;
    // The photons themselves must have the same grandmother = pi0
    const auto ph1 = mcParticles.iteratorAt(m1.motherId);
    const auto ph2 = mcParticles.iteratorAt(m2.motherId);
    if (!ph1.has_mothers() || !ph2.has_mothers())
      return false;
    const int gm1 = ph1.mothersIds()[0];
    const int gm2 = ph2.mothersIds()[0];
    if (gm1 != gm2)
      return false;
    return (std::abs(mcParticles.iteratorAt(gm1).pdgCode()) == 111);
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
    const AxisSpec axisTruthType{{0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5}, "truth type (1=TrueTrueDistinct,2=TrueTrueSamePhoton,3=SharedMcLeg,4=TrueFake,5=FakeFake,6=Pi0Daughters)"};

    static constexpr std::array<std::string_view, 6> kTypes = {
      "TrueTrueDistinct/",
      "TrueTrueSamePhoton/",
      "SharedMcLeg/",
      "TrueFake/",
      "FakeFake/",
      "Pi0Daughters/"};

    for (const auto& label : kTypes) {
      const std::string base = std::string("Pair/same/MC/") + std::string(label);

      if (cfgDo3D) {
        fRegistry.add((base + "CF_3D").c_str(), "MC CF 3D LCMS", kTHnSparseD, {axisQout, axisQside, axisQlong, axisKt}, true);
        if (cfgDo2D) {
          fRegistry.add((base + "CF_2D").c_str(), "MC CF 2D", kTHnSparseD, {axisQout, axisQinv, axisKt}, true);
        }
      } else {
        if (cfgUseLCMS) {
          fRegistry.add((base + "CF_1D").c_str(), "MC CF 1D LCMS", kTH2D, {axisQabsLcms, axisKt}, true);
        } else {
          fRegistry.add((base + "CF_1D").c_str(), "MC CF 1D (qinv)", kTH2D, {axisQinv, axisKt}, true);
        }
      }

      fRegistry.add((base + "hQinv").c_str(), "q_{inv};q_{inv} (GeV/c);counts", kTH1D, {axisQinv}, true);
      fRegistry.add((base + "hDeltaEta").c_str(), "#Delta#eta;#Delta#eta;counts", kTH1D, {axisDeltaEta}, true);
      fRegistry.add((base + "hDeltaPhi").c_str(), "#Delta#phi;#Delta#phi (rad);counts", kTH1D, {axisDeltaPhi}, true);
      fRegistry.add((base + "hDEtaDPhi").c_str(), "#Delta#eta vs #Delta#phi;#Delta#eta;#Delta#phi", kTH2D, {axisDeltaEta, axisDeltaPhi}, true);
      fRegistry.add((base + "hDeltaR").c_str(), "|R_{1}-R_{2}|;|R_{1}-R_{2}| (cm);counts", kTH1D, {axisDeltaR}, true);
      fRegistry.add((base + "hDeltaZ").c_str(), "#Delta z;#Delta z (cm);counts", kTH1D, {axisDeltaZ}, true);
      fRegistry.add((base + "hDeltaR3D").c_str(), "#Delta r_{3D};#Delta r_{3D} (cm);counts", kTH1D, {axisDeltaR3D}, true);
      fRegistry.add((base + "hKt").c_str(), "k_{T};k_{T} (GeV/c);counts", kTH1D, {axisKt}, true);
      fRegistry.add((base + "hDeltaRVsQinv").c_str(), "|R_{1}-R_{2}| vs q_{inv};q_{inv} (GeV/c);|R_{1}-R_{2}| (cm)", kTH2D, {axisQinv, axisDeltaR}, true);
      fRegistry.add((base + "hDeltaZVsQinv").c_str(), "#Delta z vs q_{inv};q_{inv} (GeV/c);#Delta z (cm)", kTH2D, {axisQinv, axisDeltaZ}, true);
      fRegistry.add((base + "hDeltaR3DVsQinv").c_str(), "#Delta r_{3D} vs q_{inv};q_{inv} (GeV/c);#Delta r_{3D} (cm)", kTH2D, {axisQinv, axisDeltaR3D}, true);
      fRegistry.add((base + "hDEtaDPhiVsQinv").c_str(), "#Delta#eta vs #Delta#phi vs q_{inv};#Delta#eta;#Delta#phi;q_{inv}", kTHnSparseD, {axisDeltaEta, axisDeltaPhi, axisQinv}, true);
      fRegistry.add((base + "hSparseDeltaRDeltaZQinv").c_str(), "|R_{1}-R_{2}|,#Delta z,q_{inv};|R_{1}-R_{2}| (cm);#Delta z (cm);q_{inv} (GeV/c)", kTHnSparseD, {axisDeltaR, axisDeltaZ, axisQinv}, true);
    }

    fRegistry.add("Pair/same/MC/hTruthTypeVsQinv", "truth type vs q_{inv};q_{inv} (GeV/c);truth type", kTH2D, {axisQinv, axisTruthType}, true);
    fRegistry.add("Pair/same/MC/hTruthTypeVsKt", "truth type vs k_{T};k_{T} (GeV/c);truth type", kTH2D, {axisKt, axisTruthType}, true);
  }

  template <PairTruthType TruthT, bool IsMix>
  inline void fillMCPairQATyped(PairQAObservables const& obs)
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

    fRegistry.fill(HIST(base) + HIST("hQinv"), obs.qinv);
    fRegistry.fill(HIST(base) + HIST("hDeltaEta"), obs.deta);
    fRegistry.fill(HIST(base) + HIST("hDeltaPhi"), obs.dphi);
    fRegistry.fill(HIST(base) + HIST("hDEtaDPhi"), obs.deta, obs.dphi);
    fRegistry.fill(HIST(base) + HIST("hDeltaR"), obs.deltaR);
    fRegistry.fill(HIST(base) + HIST("hDeltaZ"), obs.deltaZ);
    fRegistry.fill(HIST(base) + HIST("hDeltaR3D"), obs.deltaR3D);
    fRegistry.fill(HIST(base) + HIST("hKt"), obs.kt);

    constexpr auto summaryDir = IsMix ? "Pair/mix/MC/" : "Pair/same/MC/";
    const int typeIdx = static_cast<int>(TruthT);
    fRegistry.fill(HIST(summaryDir) + HIST("hTruthTypeVsQinv"), obs.qinv, typeIdx);
    fRegistry.fill(HIST(summaryDir) + HIST("hTruthTypeVsKt"), obs.kt, typeIdx);
  }

  template <bool IsMix>
  inline void fillMCPairQA(PairTruthType truthType, PairQAObservables const& obs)
  {
    switch (truthType) {
      case PairTruthType::TrueTrueDistinct:
        fillMCPairQATyped<PairTruthType::TrueTrueDistinct, IsMix>(obs);
        break;
      case PairTruthType::TrueTrueSamePhoton:
        fillMCPairQATyped<PairTruthType::TrueTrueSamePhoton, IsMix>(obs);
        break;
      case PairTruthType::SharedMcLeg:
        fillMCPairQATyped<PairTruthType::SharedMcLeg, IsMix>(obs);
        break;
      case PairTruthType::TrueFake:
        fillMCPairQATyped<PairTruthType::TrueFake, IsMix>(obs);
        break;
      case PairTruthType::FakeFake:
        fillMCPairQATyped<PairTruthType::FakeFake, IsMix>(obs);
        break;
      case PairTruthType::Pi0Daughters:
        fillMCPairQATyped<PairTruthType::Pi0Daughters, IsMix>(obs);
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

    fRegistry.fill(HIST(base) + HIST("hDeltaRVsQinv"), obs.qinv, obs.deltaR);
    fRegistry.fill(HIST(base) + HIST("hDeltaZVsQinv"), obs.qinv, obs.deltaZ);
    fRegistry.fill(HIST(base) + HIST("hDeltaR3DVsQinv"), obs.qinv, obs.deltaR3D);
    fRegistry.fill(HIST(base) + HIST("hDEtaDPhiVsQinv"), obs.deta, obs.dphi, obs.qinv);
    fRegistry.fill(HIST(base) + HIST("hSparseDeltaRDeltaZQinv"), obs.deltaR, obs.deltaZ, obs.qinv);
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

      // ── Centrality selection ──────────────────────────────────────────────
      const float cent[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (cent[cfgCentEstimator] < cfgCentMin || cfgCentMax < cent[cfgCentEstimator])
        continue;

      const std::array<float, 7> epArr = {
        collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(),
        collision.ep2fv0a(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      const float ep2 = epArr[cfgEP2EstimatorForMix];

      // ── Event QA and event cut ────────────────────────────────────────────
      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision))
        continue;
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      // ── Event mixing bins ─────────────────────────────────────────────────
      const float occupancy = (cfgOccupancyEstimator == 1)
                                ? static_cast<float>(collision.trackOccupancyInTimeRange())
                                : collision.ft0cOccupancyInTimeRange();
      const float centForQA = cent[cfgCentEstimator];

      const int zbin = binOf(ztxBinEdges, collision.posZ());
      const int centbin = binOf(centBinEdges, centForQA);
      const int epbin = binOf(epBinEgdes, ep2);
      const int occbin = binOf(occBinEdges, occupancy);

      auto keyBin = std::make_tuple(zbin, centbin, epbin, occbin);
      auto keyDFCollision = std::make_pair(ndf, collision.globalIndex());

      // ── Slice photons for this collision ──────────────────────────────────
      auto photons1Coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      auto photons2Coll = photons2.sliceBy(perCollision2, collision.globalIndex());

      // ── Single-photon QA
      if (qaflags.doSinglePhotonQa) {
        for (const auto& g : photons1Coll) {
          if (!cut1.template IsSelected<decltype(g), TSubInfos1>(g))
            continue;
          fillSinglePhotonQAStep<0>(g);
        }
      }

      std::unordered_set<int> photonIdsAfterDRCosOA;
      std::unordered_set<int> photonIdsAfterRZ;
      std::unordered_set<int> photonIdsAfterEllipse;

      // ── Same-event pair loop ──────────────────────────────────────────────
      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1Coll, photons2Coll))) {
        if (!cut1.template IsSelected<decltype(g1), TSubInfos1>(g1) ||
            !cut2.template IsSelected<decltype(g2), TSubInfos2>(g2))
          continue;

        const auto pos1 = g1.template posTrack_as<TSubInfos1>();
        const auto ele1 = g1.template negTrack_as<TSubInfos1>();
        const auto pos2 = g2.template posTrack_as<TSubInfos2>();
        const auto ele2 = g2.template negTrack_as<TSubInfos2>();
        if (pos1.trackId() == pos2.trackId() ||
            pos1.trackId() == ele2.trackId() ||
            ele1.trackId() == pos2.trackId() ||
            ele1.trackId() == ele2.trackId())
          continue;

        auto obs = buildPairQAObservables(g1, g2);
        if (!obs.valid)
          continue;

        const bool doQA = passQinvQAGate(obs.qinv);
        const bool doFullRange = passQinvFullRangeGate(obs.qinv);

        // ── QA: Before any pair cut ───────────────────────────────────────
        if (doQA)
          fillPairQAStep<0, 0>(obs, centForQA, occupancy);

        // ── Cut 1: dr/cosOA ───────────────────────────────────────────────
        if (doFullRange)
          fillFullRangeDeltaRCosOA<0>(obs.qinv, obs.drOverCosOA);
        fRegistry.fill(HIST("Pair/same/hDeltaRCosOA"), obs.drOverCosOA);
        if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA)
          continue;

        photonIdsAfterDRCosOA.insert(g1.globalIndex());
        photonIdsAfterDRCosOA.insert(g2.globalIndex());

        // ── QA: After dr/cosOA cut ────────────────────────────────────────
        if (doQA)
          fillPairQAStep<0, 1>(obs, centForQA, occupancy);

        // ── Cut 2: R/Z geometry ───────────────────────────────────────────
        if (!passRZCut(obs.deltaR, obs.deltaZ))
          continue;

        photonIdsAfterRZ.insert(g1.globalIndex());
        photonIdsAfterRZ.insert(g2.globalIndex());

        // ── QA: After R/Z cut ─────────────────────────────────────────────
        if (doQA)
          fillPairQAStep<0, 2>(obs, centForQA, occupancy);

        // ── Cut 3: Ellipse in (DeltaEta, DeltaPhi) ────────────────────────
        if (isInsideEllipse(obs.deta, obs.dphi))
          continue;

        photonIdsAfterEllipse.insert(g1.globalIndex());
        photonIdsAfterEllipse.insert(g2.globalIndex());

        // ── QA: After ellipse cut = final accepted pairs ──────────────────
        if (doQA)
          fillPairQAStep<0, 3>(obs, centForQA, occupancy);

        if (doFullRange)
          fillFullRangeQA<0>(obs, centForQA, occupancy);

        fillPairHistogram<0>(collision, obs.v1, obs.v2, 1.f);
        ndiphoton++;

        auto addToPool = [&](auto const& g) {
          if (usedPhotonIdsPerCol.insert(g.globalIndex()).second) {
            EMPair gtmp(g.pt(), g.eta(), g.phi(), 0.f);
            gtmp.setConversionPointXYZ(g.vx(), g.vy(), g.vz());
            emh1->AddTrackToEventPool(keyDFCollision, gtmp);
          }
        };
        addToPool(g1);
        addToPool(g2);
      } // end same-event pair loop

      if (qaflags.doSinglePhotonQa) {
        for (const auto& g : photons1Coll) {
          if (!cut1.template IsSelected<decltype(g), TSubInfos1>(g))
            continue;
          const int gid = g.globalIndex();
          if (photonIdsAfterDRCosOA.count(gid))
            fillSinglePhotonQAStep<1>(g);
          if (photonIdsAfterRZ.count(gid))
            fillSinglePhotonQAStep<2>(g);
          if (photonIdsAfterEllipse.count(gid))
            fillSinglePhotonQAStep<3>(g);
        }
      }

      usedPhotonIdsPerCol.clear();

      if (!cfgDoMix || ndiphoton == 0)
        continue;

      auto selectedPhotons = emh1->GetTracksPerCollision(keyDFCollision);
      auto poolIDs = emh1->GetCollisionIdsFromEventPool(keyBin);

      for (const auto& mixID : poolIDs) {
        // skip same event
        if (mixID.second == collision.globalIndex() && mixID.first == ndf)
          continue;

        const uint64_t bcMix = mapMixedEventIdToGlobalBC[mixID];
        const uint64_t diffBC = std::max(collision.globalBC(), bcMix) -
                                std::min(collision.globalBC(), bcMix);
        fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
        if (diffBC < ndiffBCMix)
          continue;

        auto poolPhotons = emh1->GetTracksPerCollision(mixID);

        for (const auto& g1 : selectedPhotons) {
          for (const auto& g2 : poolPhotons) {

            auto obs = buildPairQAObservables(g1, g2);
            if (!obs.valid)
              continue;

            const bool doQA = passQinvQAGate(obs.qinv);
            const bool doFullRange = passQinvFullRangeGate(obs.qinv);

            // ── QA: Before any pair cut ─────────────────────────────────
            if (doQA)
              fillPairQAStep<1, 0>(obs, centForQA, occupancy);

            // ── Cut 1: dr/cosOA ─────────────────────────────────────────
            if (doFullRange)
              fillFullRangeDeltaRCosOA<1>(obs.qinv, obs.drOverCosOA);
            fRegistry.fill(HIST("Pair/mix/hDeltaRCosOA"), obs.drOverCosOA);
            if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA)
              continue;

            // ── QA: After dr/cosOA cut ──────────────────────────────────
            if (doQA)
              fillPairQAStep<1, 1>(obs, centForQA, occupancy);

            // ── Cut 2: R/Z geometry ─────────────────────────────────────
            if (!passRZCut(obs.deltaR, obs.deltaZ))
              continue;

            // ── QA: After R/Z cut ───────────────────────────────────────
            if (doQA)
              fillPairQAStep<1, 2>(obs, centForQA, occupancy);

            // ── Cut 3: Ellipse ──────────────────────────────────────────
            if (isInsideEllipse(obs.deta, obs.dphi))
              continue;

            // ── QA: After ellipse cut ───────────────────────────────────
            if (doQA)
              fillPairQAStep<1, 3>(obs, centForQA, occupancy);

            // ── Full-range fills ────────────────────────────────────────
            if (doFullRange)
              fillFullRangeQA<1>(obs, centForQA, occupancy);

            // ── Fill CF histogram — always ──────────────────────────────
            fillPairHistogram<1>(collision, obs.v1, obs.v2, 1.f);
          }
        }
      } // end mixed-event loop

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(keyBin, keyDFCollision);
        emh2->AddCollisionIdAtLast(keyBin, keyDFCollision);
        mapMixedEventIdToGlobalBC[keyDFCollision] = collision.globalBC();
      }
    } // end collision loop
  }

  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<
    std::tuple<int, int, int, int>,
    std::pair<int, int>,
    EMPair>;

  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;

  std::unordered_set<int> usedPhotonIdsPerCol;
  std::map<std::pair<int, int>, uint64_t> mapMixedEventIdToGlobalBC;

  SliceCache cache;
  Preslice<MyV0Photons> perCollisionPCM = aod::v0photonkf::pmeventId;

  Filter collisionFilterCentrality =
    (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) ||
    (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) ||
    (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
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
    runPairing(collisions,
               v0photons, v0photons,
               v0legs, v0legs,
               perCollisionPCM, perCollisionPCM,
               fV0PhotonCut, fV0PhotonCut);
    ndf++;
  }

  PROCESS_SWITCH(photonhbt, processAnalysis, "pairing for analysis", true);

  template <typename TCollisions,
            typename TPhotons,
            typename TLegs,
            typename TMCParticles,
            typename TPreslice,
            typename TCut>
  void runPairingMC(TCollisions const& collisions,
                    TPhotons const& photons,
                    TLegs const& /*legs*/,
                    TMCParticles const& mcParticles,
                    TPreslice const& perCollision,
                    TCut const& cut)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      int ndiphoton = 0;

      const float cent[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (cent[cfgCentEstimator] < cfgCentMin || cfgCentMax < cent[cfgCentEstimator])
        continue;

      const std::array<float, 7> epArr = {
        collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(),
        collision.ep2fv0a(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      const float ep2 = epArr[cfgEP2EstimatorForMix];

      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision))
        continue;

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      const float occupancy = (cfgOccupancyEstimator == 1)
                                ? static_cast<float>(collision.trackOccupancyInTimeRange())
                                : collision.ft0cOccupancyInTimeRange();
      const float centForQA = cent[cfgCentEstimator];

      const int zbin = binOf(ztxBinEdges, collision.posZ());
      const int centbin = binOf(centBinEdges, centForQA);
      const int epbin = binOf(epBinEgdes, ep2);
      const int occbin = binOf(occBinEdges, occupancy);

      auto keyBin = std::make_tuple(zbin, centbin, epbin, occbin);
      auto keyDFCollision = std::make_pair(ndf, collision.globalIndex());

      auto photonsColl = photons.sliceBy(perCollision, collision.globalIndex());

      if (qaflags.doSinglePhotonQa) {
        for (const auto& g : photonsColl) {
          if (!cut.template IsSelected<decltype(g), TLegs>(g))
            continue;
          fillSinglePhotonQAStep<0>(g);
        }
      }

      std::unordered_set<int> photonIdsAfterDRCosOA;
      std::unordered_set<int> photonIdsAfterRZ;
      std::unordered_set<int> photonIdsAfterEllipse;

      // ── Same-event pair loop ──────────────────────────────────────────────
      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsColl, photonsColl))) {
        if (!cut.template IsSelected<decltype(g1), TLegs>(g1) ||
            !cut.template IsSelected<decltype(g2), TLegs>(g2))
          continue;

        const auto pos1 = g1.template posTrack_as<TLegs>();
        const auto ele1 = g1.template negTrack_as<TLegs>();
        const auto pos2 = g2.template posTrack_as<TLegs>();
        const auto ele2 = g2.template negTrack_as<TLegs>();
        if (pos1.trackId() == pos2.trackId() ||
            pos1.trackId() == ele2.trackId() ||
            ele1.trackId() == pos2.trackId() ||
            ele1.trackId() == ele2.trackId())
          continue;

        // ── MC truth classification ───────────────────────────────────────
        const auto mc1 = buildPhotonMCInfo<decltype(g1), TLegs>(g1, mcParticles);
        const auto mc2 = buildPhotonMCInfo<decltype(g2), TLegs>(g2, mcParticles);
        auto truthType = classifyPairTruth(mc1, mc2);
        if (truthType == PairTruthType::TrueTrueDistinct &&
            isPi0DaughterPair(mc1, mc2, mcParticles))
          truthType = PairTruthType::Pi0Daughters;

        auto obs = buildPairQAObservables(g1, g2);
        if (!obs.valid)
          continue;

        const bool doQA = passQinvQAGate(obs.qinv);
        const bool doFullRange = passQinvFullRangeGate(obs.qinv);

        // ── Pair QA: Before ───────────────────────────────────────────────
        if (doQA)
          fillPairQAStep<0, 0>(obs, centForQA, occupancy);

        // ── Cut 1: dr/cosOA ───────────────────────────────────────────────
        if (doFullRange)
          fillFullRangeDeltaRCosOA<0>(obs.qinv, obs.drOverCosOA);
        fRegistry.fill(HIST("Pair/same/hDeltaRCosOA"), obs.drOverCosOA);
        if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA)
          continue;

        photonIdsAfterDRCosOA.insert(g1.globalIndex());
        photonIdsAfterDRCosOA.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 1>(obs, centForQA, occupancy);

        // ── Cut 2: R/Z geometry ───────────────────────────────────────────
        if (!passRZCut(obs.deltaR, obs.deltaZ))
          continue;

        photonIdsAfterRZ.insert(g1.globalIndex());
        photonIdsAfterRZ.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 2>(obs, centForQA, occupancy);

        // ── Cut 3: Ellipse ────────────────────────────────────────────────
        if (isInsideEllipse(obs.deta, obs.dphi))
          continue;

        photonIdsAfterEllipse.insert(g1.globalIndex());
        photonIdsAfterEllipse.insert(g2.globalIndex());
        if (doQA)
          fillPairQAStep<0, 3>(obs, centForQA, occupancy);

        // ── Full-range fills ──────────────────────────────────────────────
        if (doFullRange)
          fillFullRangeQA<0>(obs, centForQA, occupancy);

        // ── Fill inclusive CF — always ────────────────────────────────────
        fillPairHistogram<0>(collision, obs.v1, obs.v2, 1.f);
        ndiphoton++;

        if (doQA)
          fillMCPairQA</*IsMix=*/false>(truthType, obs);
        if (doFullRange)
          fillMCPairQAFullRange</*IsMix=*/false>(truthType, obs);
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

        auto addToPool = [&](auto const& g) {
          if (usedPhotonIdsPerCol.insert(g.globalIndex()).second) {
            EMPair gtmp(g.pt(), g.eta(), g.phi(), 0.f);
            gtmp.setConversionPointXYZ(g.vx(), g.vy(), g.vz());
            emh1->AddTrackToEventPool(keyDFCollision, gtmp);
          }
        };
        addToPool(g1);
        addToPool(g2);
      } // end same-event pair loop

      if (qaflags.doSinglePhotonQa) {
        for (const auto& g : photonsColl) {
          if (!cut.template IsSelected<decltype(g), TLegs>(g))
            continue;
          const int gid = g.globalIndex();
          if (photonIdsAfterDRCosOA.count(gid))
            fillSinglePhotonQAStep<1>(g);
          if (photonIdsAfterRZ.count(gid))
            fillSinglePhotonQAStep<2>(g);
          if (photonIdsAfterEllipse.count(gid))
            fillSinglePhotonQAStep<3>(g);
        }
      }

      usedPhotonIdsPerCol.clear();

      if (!cfgDoMix || ndiphoton == 0)
        continue;

      auto selectedPhotons = emh1->GetTracksPerCollision(keyDFCollision);
      auto poolIDs = emh1->GetCollisionIdsFromEventPool(keyBin);

      for (const auto& mixID : poolIDs) {
        if (mixID.second == collision.globalIndex() && mixID.first == ndf)
          continue;
        const uint64_t bcMix = mapMixedEventIdToGlobalBC[mixID];
        const uint64_t diffBC = std::max(collision.globalBC(), bcMix) -
                                std::min(collision.globalBC(), bcMix);
        fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
        if (diffBC < ndiffBCMix)
          continue;

        auto poolPhotons = emh1->GetTracksPerCollision(mixID);
        for (const auto& g1 : selectedPhotons) {
          for (const auto& g2 : poolPhotons) {
            auto obs = buildPairQAObservables(g1, g2);
            if (!obs.valid)
              continue;
            const bool doQA = passQinvQAGate(obs.qinv);
            const bool doFullRange = passQinvFullRangeGate(obs.qinv);
            if (doQA)
              fillPairQAStep<1, 0>(obs, centForQA, occupancy);
            if (doFullRange)
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
            if (doFullRange)
              fillFullRangeQA<1>(obs, centForQA, occupancy);
            fillPairHistogram<1>(collision, obs.v1, obs.v2, 1.f);
          }
        }
      }

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(keyBin, keyDFCollision);
        emh2->AddCollisionIdAtLast(keyBin, keyDFCollision);
        mapMixedEventIdToGlobalBC[keyDFCollision] = collision.globalBC();
      }
    } // end collision loop
  }

  void processMC(FilteredMyCollisions const& collisions,
                 MyV0Photons const& v0photons,
                 MyMCV0Legs const& v0legs,
                 aod::EMMCParticles const& mcParticles,
                 aod::EMMCEvents const& /*mcEvents*/)
  {
    runPairingMC(collisions, v0photons, v0legs, mcParticles,
                 perCollisionPCM, fV0PhotonCut);
    ndf++;
  }

  PROCESS_SWITCH(photonhbt, processMC, "pairing with MC truth classification", false);
};

// ============================================================================
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<photonhbt>(cfgc, TaskName{"photonhbt"})};
}
