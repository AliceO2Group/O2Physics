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
/// \brief This code loops over v0 photons and makes pairs for photon HBT analysis.
/// \author Daiki Sekihata, daiki.sekihata@cern.ch and Stefanie Mrozinski stefanie.mrozinski@cern.ch

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
#include <utility>
#include <vector>

/// Single-photon track-type combo.
enum class V0Combo : int {
  Inclusive = 0,
  ItstpcItstpc = 1,   ///< both legs ITS+TPC
  ItstpcTpconly = 2,  /// one ITS+TPC leg, one TPC-only
  TpconlyTpconly = 3, ///  both legs TPC-only
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

using MyCollisions = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsPhiVPsi>;
using MyV0Photon = MyV0Photons::iterator;

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
    static constexpr std::array<std::array<int, 4>, 4> kTable = {{{0, 0, 0, 0}, {0, 1, 2, 3}, {0, 2, 4, 5}, {0, 3, 5, 6}}};
    return static_cast<PairCombo>(kTable[lo][hi]);
  }

  // ---------------------------------------------------------------------------
  // Configurables: histogram axes
  // ---------------------------------------------------------------------------

  // HBT physics
  ConfigurableAxis confQBins{"confQBins", {60, 0, +0.3f}, "q bins for output histograms"};
  ConfigurableAxis confKtBins{"confKtBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}, "kT bins"};

  // Single-photon QA
  ConfigurableAxis confPtBins{"confPtBins", {100, 0.f, 2.f}, "pT bins (GeV/c)"};
  ConfigurableAxis confEtaBins{"confEtaBins", {80, -0.8f, 0.8f}, "eta bins"};
  ConfigurableAxis confPhiBins{"confPhiBins", {90, -o2::constants::math::PI, o2::constants::math::PI}, "phi bins (rad)"};

  // Pair QA
  ConfigurableAxis confDeltaEtaBins{"confDeltaEtaBins", {100, -0.9f, +0.9f}, "Delta-eta bins"};
  ConfigurableAxis confDeltaPhiBins{"confDeltaPhiBins", {100, -o2::constants::math::PI, o2::constants::math::PI}, "Delta-phi bins (rad)"};
  ConfigurableAxis confEllipseValBins{"confEllipseValBins", {200, 0.f, 10.f}, "ellipse value bins"};
  ConfigurableAxis confCosThetaBins{"confCosThetaBins", {100, 0.f, 1.f}, "cos(theta*) bins"};
  ConfigurableAxis confOpeningAngleBins{"confOpeningAngleBins", {100, 0.f, o2::constants::math::PI}, "opening angle bins (rad)"};

  // Axis specs
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

  // ---------------------------------------------------------------------------
  // Configurables: QA flags
  // ---------------------------------------------------------------------------

  struct : ConfigurableGroup {
    std::string prefix = "qaflags_group";
    Configurable<bool> doPairQa{"doPairQa", true, "fill pair QA histograms (Before/After ellipse cut)"};
    Configurable<bool> doSinglePhotonQa{"doSinglePhotonQa", true, "fill single-photon QA histograms (pT, eta, phi)"};
  } qaflags;

  // ---------------------------------------------------------------------------
  // Configurables: HBT kind
  // ---------------------------------------------------------------------------

  Configurable<bool> cfgDo3D{"cfgDo3D", false, "enable 3D analysis"};
  Configurable<bool> cfgUseLCMS{"cfgUseLCMS", true, "measure relative momentum in LCMS for 1D"};

  // ---------------------------------------------------------------------------
  // Configurables: events
  // ---------------------------------------------------------------------------

  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity"};

  // ---------------------------------------------------------------------------
  // Configurables: mixed event
  // ---------------------------------------------------------------------------

  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  Configurable<uint64_t> ndiffBCMix{"ndiffBCMix", 594, "difference in global BC required in mixed events"};
  Configurable<int> cfgEP2EstimatorForMix{"cfgEP2EstimatorForMix", 3, "FT0M:0, FT0A:1, FT0C:2, FV0A:3, BTot:4, BPos:5, BNeg:6"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};

  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {16, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}, "Mixing bins - EP angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  // ---------------------------------------------------------------------------
  // Configurables: pair cuts
  // ---------------------------------------------------------------------------

  struct : ConfigurableGroup {
    std::string prefix = "ggpaircut_group";
    Configurable<float> cfgMinDR_CosOA{"cfgMinDR_CosOA", -1, "min. dr/cosOA for kPCMPCM"};
    Configurable<bool> cfgApplyEllipseCut{"cfgApplyEllipseCut", false, "reject pairs inside ellipse in DeltaEta-DeltaPhi"};
    Configurable<float> cfgEllipseSigEta{"cfgEllipseSigEta", 0.02f, "sigma_eta for ellipse cut"};
    Configurable<float> cfgEllipseSigPhi{"cfgEllipseSigPhi", 0.02f, "sigma_phi for ellipse cut"};
    Configurable<float> cfgEllipseR2{"cfgEllipseR2", 1.0f, "R^2 threshold: reject if value < R^2"};
  } ggpaircuts;

  // ---------------------------------------------------------------------------
  // Event cut
  // ---------------------------------------------------------------------------

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

  // ---------------------------------------------------------------------------
  // PCM cut
  // ---------------------------------------------------------------------------

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
    usedPhotonIdsPerCol.shrink_to_fit();
  }

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::mt19937 engine;
  std::uniform_int_distribution<int> dist01;
  int mRunNumber;

  std::vector<float> ztxBinEdges;
  std::vector<float> centBinEdges;
  std::vector<float> epBinEgdes;
  std::vector<float> occBinEdges;

  inline bool isInsideEllipse(float deta, float dphi) const
  {
    if (!ggpaircuts.cfgApplyEllipseCut.value)
      return false;
    const float sE = ggpaircuts.cfgEllipseSigEta.value;
    const float sP = ggpaircuts.cfgEllipseSigPhi.value;
    if (sE < 1e-9f || sP < 1e-9f)
      return false;
    return (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP) < ggpaircuts.cfgEllipseR2.value;
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

  /// Clamp bin index to valid range [0, nmax].
  static int clampBin(int b, int nmax) { return std::clamp(b, 0, nmax); }

  /// Find the bin index for val in a sorted edge vector.
  static int binOf(const std::vector<float>& edges, float val)
  {
    const int b = static_cast<int>(std::lower_bound(edges.begin(), edges.end(), val) - edges.begin()) - 1;
    return clampBin(b, static_cast<int>(edges.size()) - 2);
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
  void addQAHistogramsForStep(const std::string& path)
  {
    fRegistry.add((path + "hDeltaEta").c_str(), "#Delta#eta;#Delta#eta;counts", kTH1F, {axisDeltaEta}, true);
    fRegistry.add((path + "hDeltaPhi").c_str(), "#Delta#phi;#Delta#phi (rad);counts", kTH1F, {axisDeltaPhi}, true);
    fRegistry.add((path + "hDEtaDPhi").c_str(), "#Delta#eta vs #Delta#phi;#Delta#eta;#Delta#phi (rad)", kTH2F, {axisDeltaEta, axisDeltaPhi}, true);
    fRegistry.add((path + "hDeltaEtaVsPairEta").c_str(), "#Delta#eta vs #LT#eta#GT_{pair};#LT#eta#GT_{pair};#Delta#eta", kTH2F, {axisEta, axisDeltaEta}, true);
    fRegistry.add((path + "hCosTheta").c_str(), "cos(#theta*) in pair rest frame;cos(#theta*);counts", kTH1F, {axisCosTheta}, true);
    fRegistry.add((path + "hOpeningAngle").c_str(), "Opening angle between conversion points;#alpha (rad);counts", kTH1F, {axisOpeningAngle}, true);
    fRegistry.add((path + "hEllipseVal").c_str(), "(#Delta#eta/#sigma_{#eta})^{2}+(#Delta#phi/#sigma_{#phi})^{2};value;counts", kTH1D, {axisEllipseVal}, true);
  }

  void addhistograms()
  {
    static constexpr std::string_view det[6] = {"FT0M", "FT0A", "FT0C", "BTot", "BPos", "BNeg"};
    fRegistry.add("Event/before/hEP2_CentFT0C_forMix", Form("2nd harmonics EP for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", det[cfgEP2EstimatorForMix].data()), kTH2F, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);
    fRegistry.add("Event/after/hEP2_CentFT0C_forMix", Form("2nd harmonics EP for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", det[cfgEP2EstimatorForMix].data()), kTH2F, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);

    // ── Single-photon QA ─────────────────────────────────────────────────────
    fRegistry.add("SinglePhoton/hPt", "V0 photon p_{T};p_{T} (GeV/c);counts", kTH1F, {axisPt}, true);
    fRegistry.add("SinglePhoton/hEta", "V0 photon #eta;#eta;counts", kTH1F, {axisEta}, true);
    fRegistry.add("SinglePhoton/hPhi", "V0 photon #phi;#phi (rad);counts", kTH1F, {axisPhi}, true);
    fRegistry.add("SinglePhoton/hEtaVsPhi", "V0 photon acceptance;#phi (rad);#eta", kTH2F, {axisPhi, axisEta}, true);

    // ── HBT physics ──────────────────────────────────────────────────────────
    if (cfgDo3D) {
      fRegistry.add("Pair/same/hs_3d", "diphoton correlation 3D LCMS", kTHnSparseD, {axisQout, axisQside, axisQlong, axisKt}, true);
    } else {
      if (cfgUseLCMS) {
        fRegistry.add("Pair/same/hs_1d", "diphoton correlation 1D LCMS", kTHnSparseD, {axisQabsLcms, axisKt}, true);
      } else {
        fRegistry.add("Pair/same/hs_1d", "diphoton correlation 1D", kTHnSparseD, {axisQinv, axisKt}, true);
      }
    }

    fRegistry.add("Pair/same/hDeltaRCosOA", "distance between 2 conversion points;#Deltar/cos(#theta_{op}/2) (cm)", kTH1D, {{100, 0, 100}}, true);

    addQAHistogramsForStep("Pair/same/QA/Before/");
    addQAHistogramsForStep("Pair/same/QA/After/");

    fRegistry.addClone("Pair/same/", "Pair/mix/");
  }

  // ---------------------------------------------------------------------------
  // DefineEMEventCut / DefinePCMCut
  // ---------------------------------------------------------------------------

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

  template <int ev_id, typename TCollision>
  void fillPairHistogram(TCollision const&,
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
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_3d"),
                     std::fabs(qout_lcms), std::fabs(qside_lcms), std::fabs(qlong_lcms), kt, weight);
    } else {
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_1d"),
                     cfgUseLCMS ? qabs_lcms : qinv, kt, weight);
    }
  }

  template <int ev_id, bool IsBefore>
  inline void fillPairQAStep(float deta, float dphi, float pairEta,
                             float cosTheta, float openingAngle)
  {
    if (!qaflags.doPairQa)
      return;

    const float sE = ggpaircuts.cfgEllipseSigEta.value;
    const float sP = ggpaircuts.cfgEllipseSigPhi.value;

    if constexpr (ev_id == 0 && IsBefore) {
      fRegistry.fill(HIST("Pair/same/QA/Before/hDeltaEta"), deta);
      fRegistry.fill(HIST("Pair/same/QA/Before/hDeltaPhi"), dphi);
      fRegistry.fill(HIST("Pair/same/QA/Before/hDEtaDPhi"), deta, dphi);
      fRegistry.fill(HIST("Pair/same/QA/Before/hDeltaEtaVsPairEta"), pairEta, deta);
      fRegistry.fill(HIST("Pair/same/QA/Before/hCosTheta"), cosTheta);
      fRegistry.fill(HIST("Pair/same/QA/Before/hOpeningAngle"), openingAngle);
      if (sE > 1e-9f && sP > 1e-9f)
        fRegistry.fill(HIST("Pair/same/QA/Before/hEllipseVal"), (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP));
    } else if constexpr (ev_id == 0 && !IsBefore) {
      fRegistry.fill(HIST("Pair/same/QA/After/hDeltaEta"), deta);
      fRegistry.fill(HIST("Pair/same/QA/After/hDeltaPhi"), dphi);
      fRegistry.fill(HIST("Pair/same/QA/After/hDEtaDPhi"), deta, dphi);
      fRegistry.fill(HIST("Pair/same/QA/After/hDeltaEtaVsPairEta"), pairEta, deta);
      fRegistry.fill(HIST("Pair/same/QA/After/hCosTheta"), cosTheta);
      fRegistry.fill(HIST("Pair/same/QA/After/hOpeningAngle"), openingAngle);
      if (sE > 1e-9f && sP > 1e-9f)
        fRegistry.fill(HIST("Pair/same/QA/After/hEllipseVal"), (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP));
    } else if constexpr (ev_id == 1 && IsBefore) {
      fRegistry.fill(HIST("Pair/mix/QA/Before/hDeltaEta"), deta);
      fRegistry.fill(HIST("Pair/mix/QA/Before/hDeltaPhi"), dphi);
      fRegistry.fill(HIST("Pair/mix/QA/Before/hDEtaDPhi"), deta, dphi);
      fRegistry.fill(HIST("Pair/mix/QA/Before/hDeltaEtaVsPairEta"), pairEta, deta);
      fRegistry.fill(HIST("Pair/mix/QA/Before/hCosTheta"), cosTheta);
      fRegistry.fill(HIST("Pair/mix/QA/Before/hOpeningAngle"), openingAngle);
      if (sE > 1e-9f && sP > 1e-9f)
        fRegistry.fill(HIST("Pair/mix/QA/Before/hEllipseVal"), (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP));
    } else {
      fRegistry.fill(HIST("Pair/mix/QA/After/hDeltaEta"), deta);
      fRegistry.fill(HIST("Pair/mix/QA/After/hDeltaPhi"), dphi);
      fRegistry.fill(HIST("Pair/mix/QA/After/hDEtaDPhi"), deta, dphi);
      fRegistry.fill(HIST("Pair/mix/QA/After/hDeltaEtaVsPairEta"), pairEta, deta);
      fRegistry.fill(HIST("Pair/mix/QA/After/hCosTheta"), cosTheta);
      fRegistry.fill(HIST("Pair/mix/QA/After/hOpeningAngle"), openingAngle);
      if (sE > 1e-9f && sP > 1e-9f)
        fRegistry.fill(HIST("Pair/mix/QA/After/hEllipseVal"), (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP));
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
      if (cent[cfgCentEstimator] < cfgCentMin || cfgCentMax < cent[cfgCentEstimator])
        continue;

      const std::array<float, 7> epArr = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(),
                                          collision.ep2fv0a(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      float ep2 = epArr[cfgEP2EstimatorForMix];

      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision))
        continue;
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      // Mixing-bin indices — uses static binOf helper:
      int zbin = binOf(ztxBinEdges, collision.posZ());
      int centbin = binOf(centBinEdges, cent[cfgCentEstimator]);
      int epbin = binOf(epBinEgdes, ep2);
      int occbin = binOf(occBinEdges,
                         cfgOccupancyEstimator == 1
                           ? static_cast<float>(collision.trackOccupancyInTimeRange())
                           : collision.ft0cOccupancyInTimeRange());

      auto keyBin = std::make_tuple(zbin, centbin, epbin, occbin);
      auto keyDFCollision = std::make_pair(ndf, collision.globalIndex());

      auto photons1Coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      auto photons2Coll = photons2.sliceBy(perCollision2, collision.globalIndex());

      // ── Single-photon QA ─────────────────────────────────────────────────
      if (qaflags.doSinglePhotonQa) {
        for (const auto& g : photons1Coll) {
          if (!cut1.template IsSelected<decltype(g), TSubInfos1>(g))
            continue;
          fRegistry.fill(HIST("SinglePhoton/hPt"), g.pt());
          fRegistry.fill(HIST("SinglePhoton/hEta"), g.eta());
          fRegistry.fill(HIST("SinglePhoton/hPhi"), g.phi());
          fRegistry.fill(HIST("SinglePhoton/hEtaVsPhi"), g.phi(), g.eta());
        }
      }

      // ── Same-event pair loop ──────────────────────────────────────────────
      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1Coll, photons2Coll))) {
        if (!cut1.template IsSelected<decltype(g1), TSubInfos1>(g1) ||
            !cut2.template IsSelected<decltype(g2), TSubInfos2>(g2))
          continue;

        auto pos1 = g1.template posTrack_as<TSubInfos1>();
        auto ele1 = g1.template negTrack_as<TSubInfos1>();
        auto pos2 = g2.template posTrack_as<TSubInfos2>();
        auto ele2 = g2.template negTrack_as<TSubInfos2>();
        if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId())
          continue;

        ROOT::Math::XYZVector cp1(g1.vx(), g1.vy(), g1.vz());
        ROOT::Math::XYZVector cp2(g2.vx(), g2.vy(), g2.vz());
        float dr = std::sqrt(std::pow(g1.vx() - g2.vx(), 2) + std::pow(g1.vy() - g2.vy(), 2) + std::pow(g1.vz() - g2.vz(), 2));
        float opa = std::acos(std::clamp(static_cast<float>(cp1.Dot(cp2) / (std::sqrt(cp1.Mag2()) * std::sqrt(cp2.Mag2()))), -1.f, 1.f));
        o2::math_utils::bringTo02Pi(opa);
        if (opa > o2::constants::math::PI)
          opa -= o2::constants::math::PI;
        float cosOA = std::cos(opa / 2.f);
        if (dr / cosOA < ggpaircuts.cfgMinDR_CosOA)
          continue;
        fRegistry.fill(HIST("Pair/same/hDeltaRCosOA"), dr / cosOA);

        // Kinematic variables for QA
        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        float deta = g1.eta() - g2.eta();
        float dphi = RecoDecay::constrainAngle(g1.phi() - g2.phi(), -o2::constants::math::PI);
        float pairEta = 0.5f * (g1.eta() + g2.eta());
        float cosTheta = std::fabs(computeCosTheta(v1, v2));
        float openingAngle = opa;

        // ── QA: Before ellipse cut ──────────────────────────────────────
        fillPairQAStep<0, true>(deta, dphi, pairEta, cosTheta, openingAngle);

        // ── Ellipse cut ─────────────────────────────────────────────────
        if (isInsideEllipse(deta, dphi))
          continue;

        // ── QA: After ellipse cut ───────────────────────────────────────
        fillPairQAStep<0, false>(deta, dphi, pairEta, cosTheta, openingAngle);

        // ── Physics ─────────────────────────────────────────────────
        fillPairHistogram<0>(collision, v1, v2, 1.f);
        ndiphoton++;

        auto addToPool = [&](auto const& g) {
          if (std::find(usedPhotonIdsPerCol.begin(), usedPhotonIdsPerCol.end(),
                        g.globalIndex()) == usedPhotonIdsPerCol.end()) {
            EMPair gtmp(g.pt(), g.eta(), g.phi(), 0);
            gtmp.setConversionPointXYZ(g.vx(), g.vy(), g.vz());
            emh1->AddTrackToEventPool(keyDFCollision, gtmp);
            usedPhotonIdsPerCol.emplace_back(g.globalIndex());
          }
        };
        addToPool(g1);
        addToPool(g2);

        // end same-event pair loop
      }

      usedPhotonIdsPerCol.clear();
      usedPhotonIdsPerCol.shrink_to_fit();

      // ── Mixed-event loop ────────────────────────────────────────────────────
      if (!cfgDoMix || !(ndiphoton > 0))
        continue;

      auto selectedPhotons = emh1->GetTracksPerCollision(keyDFCollision);
      auto poolIDs = emh1->GetCollisionIdsFromEventPool(keyBin);

      for (const auto& mixID : poolIDs) {
        if (mixID.second == collision.globalIndex() && mixID.first == ndf)
          continue;

        uint64_t bcMix = mapMixedEventIdToGlobalBC[mixID];
        uint64_t diffBC = std::max(collision.globalBC(), bcMix) - std::min(collision.globalBC(), bcMix);
        fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
        if (diffBC < ndiffBCMix)
          continue;

        auto poolPhotons = emh1->GetTracksPerCollision(mixID);

        for (const auto& g1 : selectedPhotons) {
          for (const auto& g2 : poolPhotons) {

            ROOT::Math::XYZVector cp1(g1.vx(), g1.vy(), g1.vz());
            ROOT::Math::XYZVector cp2(g2.vx(), g2.vy(), g2.vz());
            float dr = std::sqrt(std::pow(g1.vx() - g2.vx(), 2) + std::pow(g1.vy() - g2.vy(), 2) + std::pow(g1.vz() - g2.vz(), 2));
            float opa = std::acos(std::clamp(static_cast<float>(cp1.Dot(cp2) / (std::sqrt(cp1.Mag2()) * std::sqrt(cp2.Mag2()))), -1.f, 1.f));
            o2::math_utils::bringTo02Pi(opa);
            if (opa > o2::constants::math::PI)
              opa -= o2::constants::math::PI;
            float cosOA = std::cos(opa / 2.f);
            if (dr / cosOA < ggpaircuts.cfgMinDR_CosOA)
              continue;
            fRegistry.fill(HIST("Pair/mix/hDeltaRCosOA"), dr / cosOA);

            ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
            ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
            float deta = g1.eta() - g2.eta();
            float dphi = RecoDecay::constrainAngle(g1.phi() - g2.phi(), -o2::constants::math::PI);
            float pairEta = 0.5f * (g1.eta() + g2.eta());
            float cosTheta = std::fabs(computeCosTheta(v1, v2));
            float openingAngle = opa;

            // QA Before cut — mix/QA/Before/ histograms.
            fillPairQAStep<1, true>(deta, dphi, pairEta, cosTheta, openingAngle);

            // Apply ellipse cut
            if (isInsideEllipse(deta, dphi))
              continue;

            // QA After cut
            fillPairQAStep<1, false>(deta, dphi, pairEta, cosTheta, openingAngle);

            fillPairHistogram<1>(collision, v1, v2, 1.f);
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

  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMPair>;
  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;
  std::vector<int> usedPhotonIdsPerCol;
  std::map<std::pair<int, int>, uint64_t> mapMixedEventIdToGlobalBC;

  SliceCache cache;
  Preslice<MyV0Photons> perCollisionPCM = aod::v0photonkf::pmeventId;

  Filter collisionFilterCentrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) ||
                                     (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) ||
                                     (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilterOccupancyTrack = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange &&
                                         o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilterOccupancyFT0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange &&
                                        o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  int ndf = 0;
  void processAnalysis(FilteredMyCollisions const& collisions, MyV0Photons const& v0photons, aod::V0Legs const& v0legs)
  {
    runPairing(collisions, v0photons, v0photons, v0legs, v0legs,
               perCollisionPCM, perCollisionPCM, fV0PhotonCut, fV0PhotonCut);
    ndf++;
  }
  PROCESS_SWITCH(photonhbt, processAnalysis, "pairing for analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<photonhbt>(cfgc, TaskName{"photonhbt"})};
}
