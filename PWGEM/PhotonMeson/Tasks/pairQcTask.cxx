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

/// \file pairQCTask.cxx
/// \brief This task looks at the pair of V0 candidates.
/// \author Stefanie Mrozinski, stefanie.mrozinski@cern.ch

#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/AnalyticV0.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector3D.h> // IWYU pragma: keep
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::pwgem::photonmeson;
namespace mcutil = o2::aod::pwgem::photonmeson::utils::mcutil;
namespace pairutil = o2::aod::pwgem::photonmeson::utils::pairutil;

// ─── Event Information Tables ────────────────────────────────────────────

using MyCollisions = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001>;
using MyCollisionsMC = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001, aod::EMMCEventLabels>;

// ─── Photon Tables ────────────────────────────────────────────

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsPhiVPsi>;

// ─── Assoc. Track Tables ────────────────────────────────────────────

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

struct PairQCTask {

  static constexpr float KMinMagnitude = 1e-12f;
  static constexpr float KMinCosine = 1e-12f;
  static constexpr float KMinSigma = 1e-9f;
  static constexpr float KMe = 0.000510999f; // electron mass, GeV/c^2

  struct LegLite {
    float pt{0.f}, eta{0.f}, phi{0.f};
    float nsigEl{0.f}, dcaZ{0.f};
    float nClsFindable{0.f}, nClsITS{0.f};
    bool timed{false}; // ITS, TRD or TOF: the track has its own time
    bool foreign{false};
    int64_t trackId{-1};
  };

  struct PhotonCand {
    int64_t gi{-1};
    float fPt{0.f}, fEta{0.f}, fPhi{0.f};
    float fVx{0.f}, fVy{0.f}, fVz{0.f}; // conversion point, global (cm)
    float fVtxZ{0.f};                   // primary vertex z of the collision
    float fDcaXYToPV{0.f}, fDcaZToPV{0.f};
    float psipair{0.f}, phiv{0.f}, chi2{0.f}, cospa{0.f}, pca{0.f};
    std::array<LegLite, 2> leg{}; // 0 = e+, 1 = e-
    pairutil::V0PhotonLegCounts legCounts{};
    pairutil::DedupCand dedupCand{};
    bool dedupRejected{false};
    AnalyticV0 own{}; // the photon rebuilt from its own legs (WP-B, WP-F)
    mcutil::PhotonMCInfo mc{};
    [[nodiscard]] float pt() const { return fPt; }
    [[nodiscard]] float eta() const { return fEta; }
    [[nodiscard]] float phi() const { return fPhi; }
    [[nodiscard]] float vx() const { return fVx; }
    [[nodiscard]] float vy() const { return fVy; }
    [[nodiscard]] float vz() const { return fVz; }
    [[nodiscard]] float rConv() const { return std::hypot(fVx, fVy); }
    [[nodiscard]] int nTimedLegs() const { return (leg[0].timed ? 1 : 0) + (leg[1].timed ? 1 : 0); }
    [[nodiscard]] int nForeignLegs() const { return (leg[0].foreign ? 1 : 0) + (leg[1].foreign ? 1 : 0); }
    [[nodiscard]] bool sharesTrackWith(PhotonCand const& o) const
    {
      return leg[0].trackId == o.leg[0].trackId || leg[0].trackId == o.leg[1].trackId ||
             leg[1].trackId == o.leg[0].trackId || leg[1].trackId == o.leg[1].trackId;
    }
    [[nodiscard]] LegParam legParam(int il) const
    {
      LegParam l;
      l.xyz = {fVx, fVy, fVz - fVtxZ};
      const float ptL = leg[static_cast<size_t>(il)].pt;
      l.mom = {ptL * std::cos(leg[static_cast<size_t>(il)].phi), ptL * std::sin(leg[static_cast<size_t>(il)].phi), ptL * std::sinh(leg[static_cast<size_t>(il)].eta)};
      l.charge = (il == 0) ? +1 : -1;
      return l;
    }
  };

  // Pair observables from the two conversion points and momenta.
  struct PairObs {
    ROOT::Math::PtEtaPhiMVector v1, v2;
    float deltaR{0.f}, deltaZ{0.f}, deltaR3D{0.f};
    float opa{0.f}, drOverCosOA{0.f};
    float deta{0.f}, dphi{0.f};
    pairutil::PairQ q{};
    bool valid{true};
  };

  struct CrossObs {
    std::array<AltPairing, 2> fit{};
    std::array<AnalyticV0, 2> ana{};
    int nAltValidFit{0}, nAltValidAna{0};
    float maxDcaFit{999.f}, maxDcaXYFit{999.f}, maxDcaZFit{999.f};
    float maxDcaAna{999.f}, maxDcaXYAna{999.f}, maxDcaZAna{999.f};
    std::array<float, 2> mee{999.f, 999.f};
    float meeOverQ{999.f};
    std::array<float, 2> ownMee{999.f, 999.f};
    float ownMeeMax{999.f}, deltaMee{-999.f}, m4{999.f};
    bool evaluated{false};
  };

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  Configurable<std::string> cfgCcdbUrl{"cfgCcdbUrl", "http://alice-ccdb.cern.ch", "CCDB url"};
  Configurable<float> cfgBzOverrideT{"cfgBzOverrideT", -999.f, "Bz in Tesla; used instead of CCDB if > -100"};

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

  struct : ConfigurableGroup {
    std::string prefix = "centralitySelection_group";
    Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
    Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  } centralitySelection;

  struct : ConfigurableGroup {
    std::string prefix = "ggpaircut_group";
    Configurable<float> cfgMinDRCosOA{"cfgMinDRCosOA", -1.f, "min. dr/cosOA; <0 = disabled"};
    Configurable<bool> cfgDoRCut{"cfgDoRCut", false, "apply |R1-R2| > cfgMinDeltaR cut"};
    Configurable<float> cfgMinDeltaR{"cfgMinDeltaR", 0.f, "minimum |R1-R2| (cm)"};
    Configurable<bool> cfgDoZCut{"cfgDoZCut", false, "apply |DeltaZ| > cfgMinDeltaZ cut"};
    Configurable<float> cfgMinDeltaZ{"cfgMinDeltaZ", 0.f, "minimum |DeltaZ| (cm)"};
    Configurable<bool> cfgDoEllipseCut{"cfgDoEllipseCut", false, "reject pairs inside ellipse in DeltaEta-DeltaPhi"};
    Configurable<float> cfgEllipseSigEta{"cfgEllipseSigEta", 0.1f, "sigma_eta for ellipse cut"};
    Configurable<float> cfgEllipseSigPhi{"cfgEllipseSigPhi", 0.1f, "sigma_phi for ellipse cut"};
    Configurable<float> cfgEllipseR2{"cfgEllipseR2", 1.0f, "R^2 threshold: reject if ellipse value < R^2"};
    Configurable<float> cfgMaxAsymmetry{"cfgMaxAsymmetry", -1.f, "max |p_{T, 1} - p_{T, 2}|/(p_{T, 1} + p_{T, 2}) asymmetry cut"};
    Configurable<float> cfgMaxDcaZToPV{"cfgMaxDcaZToPV", 999.f, "max |DCAz| of both V0 photons to the primary vertex (cm); 999 = off"};
    Configurable<float> cfgMaxDcaXYToPV{"cfgMaxDcaXYToPV", 999.f, "max |DCAxy| of both V0 photons to the primary vertex (cm); 999 = off"};
  } ggpaircuts;

  struct : ConfigurableGroup {
    std::string prefix = "dedup_group";
    Configurable<bool> cfgDoDedup{"cfgDoDedup", false, "REMOVE duplicate-like candidates before pairing (as photonhbt would); the QA is filled either way"};
    Configurable<float> cfgDupMaxLegDEta{"cfgDupMaxLegDEta", 0.01f, "leg identical if |dEta| below this"};
    Configurable<float> cfgDupMaxLegDPhi{"cfgDupMaxLegDPhi", 0.01f, "leg identical if |dPhi| below this (rad)"};
    Configurable<float> cfgDupMaxLegPtAsym{"cfgDupMaxLegPtAsym", 0.05f, "leg identical if |pt1-pt2|/(pt1+pt2) below this"};
    Configurable<float> cfgDupMaxLegDeDxAsym{"cfgDupMaxLegDeDxAsym", 0.05f, "leg identical if |dEdx1-dEdx2|/(dEdx1+dEdx2) below this; skipped if a leg has no TPC"};
    Configurable<bool> cfgDupRequireBothLegs{"cfgDupRequireBothLegs", false, "true: both same-charge leg pairs must be identical; false: one is enough"};
    Configurable<float> cfgDupMaxDVtx3D{"cfgDupMaxDVtx3D", 999.f, "additionally require the 3D conversion-point distance below this (cm); 999 = off"};
  } dedup;

  struct : ConfigurableGroup {
    std::string prefix = "crosspair_group";
    Configurable<bool> cfgDoCrossPairCut{"cfgDoCrossPairCut", false, "apply the crossed-pair veto (DCAFitter alternatives) as the last pair cut, as photonhbt does"};
    Configurable<int> cfgCrossVetoMode{"cfgCrossVetoMode", 1, "veto criterion: 1 = geometric AND-chain (cfgAlt*), 2 = maxAltDca alone (cfgCrossMinMaxAltDca), 3 = 1 OR 2"};
    Configurable<float> cfgCrossMinMaxAltDca{"cfgCrossMinMaxAltDca", 4.9f, "veto mode 2/3: reject when maxAltDca is BELOW this (cm); 0 disables"};
    Configurable<float> cfgAltMaxQinv{"cfgAltMaxQinv", 0.15f, "evaluate the alternative pairings only below this q_inv"};
    Configurable<float> cfgAltMaxDca{"cfgAltMaxDca", 1.5f, "alternative pairing photon-like if its legs approach closer than this at their PCA (cm) ..."};
    Configurable<float> cfgAltMinR{"cfgAltMinR", 1.f, "... and the PCA lies at a radius above this (cm) ..."};
    Configurable<float> cfgAltMaxR{"cfgAltMaxR", 90.f, "... and below this (cm) ..."};
    Configurable<float> cfgAltMinCosPA{"cfgAltMinCosPA", 0.98f, "... and cosPA above this ..."};
    Configurable<float> cfgAltMaxDcaXY{"cfgAltMaxDcaXY", 3.f, "... and DCAxy of the photon line to the vertex below this (cm) ..."};
    Configurable<float> cfgAltMaxDcaZ{"cfgAltMaxDcaZ", 999.f, "... and DCAz below this (cm); 999 = not used"};
    Configurable<float> cfgAltMaxMee{"cfgAltMaxMee", 0.1f, "... and its m_ee below this (GeV/c^2)"};
    Configurable<bool> cfgAltRequireBoth{"cfgAltRequireBoth", true, "true: veto only if BOTH alternative pairings are photon-like; false: one is enough"};
  } crosspair;

  struct : ConfigurableGroup {
    std::string prefix = "qaflags_group";
    Configurable<float> cfgMaxQinvForQA{"cfgMaxQinvForQA", 0.3f, "pair QA only below this q_inv (GeV/c); the cheap gate runs before any boost"};
    Configurable<float> cfgMaxQinvForMCQA{"cfgMaxQinvForMCQA", 0.3f, "MC pair diagnostics (ancestry, true conversions) only below this q_inv"};
    Configurable<bool> cfgDoCrossPairQA{"cfgDoCrossPairQA", true, "WP-F: evaluate the alternative pairings (two DCAFitter fits AND two analytic builds per pair below cfgAltMaxQinv) and fill Pair/CrossPair/*"};
    Configurable<bool> cfgDoAnalyticV0QA{"cfgDoAnalyticV0QA", true, "WP-B: rebuild every selected photon from its own legs and fill the residuals against the builder (Photon/AnalyticV0/*)"};
    Configurable<bool> cfgDoAncestry{"cfgDoAncestry", true, "WP-E (MC): census of the four legs, common ancestor, true conversion points (Pair/MC/Ancestry/*)"};
    Configurable<int> cfgAncestorMaxGen{"cfgAncestorMaxGen", 4, "how many MC generations to walk up when looking for a common ancestor"};
    Configurable<bool> cfgFillLegSimilarity{"cfgFillLegSimilarity", true, "like-sign leg similarity sparse (duplicate diagnostic)"};
    Configurable<bool> cfgFillTypeSparses{"cfgFillTypeSparses", true, "MC: per truth type also the (dEta, dPhi, q_inv) and (dR, dZ, q_inv) sparses (large)"};
    Configurable<float> cfgScoreWeight{"cfgScoreWeight", 0.5f, "builder score weight w for the analytic V0 (as mixing_group.cfgScoreWeight of photonhbt)"};
    Configurable<float> cfgMaxQinvForSparses{"cfgMaxQinvForSparses", 0.1f, "the sparses with 4 or more axes carry the COARSE q axis (confQBinsCoarse) and are filled only below this q_inv; keep it equal to the upper edge of confQBinsCoarse"};
  } qaflags;

  ConfigurableAxis confQBins{"confQBins", {60, 0, +0.3f}, "q bins"};
  ConfigurableAxis confQBinsCoarse{"confQBinsCoarse", {VARIABLE_WIDTH, 0.0, 0.01, 0.02, 0.03, 0.05, 0.1}, "coarse q bins for the sparses with 4 or more axes (merging size)"};
  ConfigurableAxis confKtBins{"confKtBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75}, "kT bins"};

  AxisSpec axisQinv{confQBins, "q_{inv} (GeV/c)"};
  AxisSpec axisQinvCoarse{confQBinsCoarse, "q_{inv} (GeV/c), coarse"};
  AxisSpec axisKt{confKtBins, "k_{T} (GeV/c)"};
  AxisSpec axisDeltaEta{90, -1.6f, +1.6f, "#Delta#eta"};
  AxisSpec axisDeltaPhi{90, -o2::constants::math::PI, +o2::constants::math::PI, "#Delta#phi (rad)"};
  AxisSpec axisDeltaR{90, 0.f, 30.f, "|R_{1}-R_{2}| (cm)"};
  AxisSpec axisDeltaZ{200, -100.f, 100.f, "#Delta z (cm)"};
  AxisSpec axisDeltaR3D{100, 0.f, 100.f, "|#vec{r}_{1}-#vec{r}_{2}| (cm)"};
  AxisSpec axisTruthType{{0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5}, "truth type (1=TrueTrueDistinct,2=TrueTrueSamePhoton,3=SharedMcLeg,4=TrueFake,5=FakeFake,6=Pi0Daughters)"};
  AxisSpec axisMethod{2, -0.5f, 1.5f, "0 = DCAFitter, 1 = analytic"};
  [[nodiscard]] static AxisSpec makeAxisEnum(int nClasses, const char* title)
  {
    return AxisSpec{nClasses, -0.5f, static_cast<float>(nClasses) - 0.5f, title};
  }

  static constexpr int kNPhotonSources = 6;
  static constexpr const char* kPhotonSourceTitle = "parent photon source (0 direct, 1 #pi^{0}, 2 #eta, 3 other hadron, 4 electron = bremsstrahlung, 5 no photon parent)";
  template <typename TMCParticles>
  [[nodiscard]] static int photonSource(TMCParticles const& mcParticles, int photonId)
  {
    if (photonId < 0) {
      return 5; // o2-linter: disable=magic-number (source code)
    }
    const auto photon = mcParticles.iteratorAt(photonId);
    if (!photon.has_mothers()) {
      return 0;
    }
    const int motherPdg = std::abs(mcParticles.iteratorAt(photon.mothersIds()[0]).pdgCode());
    if (motherPdg == PDG_t::kPi0) {
      return 1;
    }
    if (motherPdg == o2::constants::physics::Pdg::kEta) {
      return 2; // o2-linter: disable=magic-number (source code)
    }
    if (motherPdg == PDG_t::kElectron) {
      return 4; // o2-linter: disable=magic-number (source code)
    }
    return 3; // o2-linter: disable=magic-number (source code)
  }

  // Track-class combination of the two legs of a V0, unordered:
  // 0 ITSTPC+ITSTPC, 1 ITSTPC+ITSonly, 2 ITSTPC+TPConly, 3 ITSonly+ITSonly, 4 ITSonly+TPConly, 5 TPConly+TPConly
  static constexpr int kNV0Classes = 6;
  static constexpr const char* kV0ClassTitle = "V0 leg classes (0 ITSTPC+ITSTPC, 1 ITSTPC+ITSonly, 2 ITSTPC+TPConly, 3 ITSonly+ITSonly, 4 ITSonly+TPConly, 5 TPConly+TPConly)";
  [[nodiscard]] static int v0LegClass(pairutil::V0PhotonLegCounts const& c)
  {
    if (c.nITSTPC == 2) { // o2-linter: disable=magic-number (two legs)
      return 0;
    }
    if (c.nITSTPC == 1) {
      return (c.nITSOnly == 1) ? 1 : 2; // o2-linter: disable=magic-number (class code)
    }
    if (c.nITSOnly == 2) { // o2-linter: disable=magic-number (two legs)
      return 3;            // o2-linter: disable=magic-number (class code)
    }
    return (c.nITSOnly == 1) ? 4 : 5; // o2-linter: disable=magic-number (class code)
  }

  HistogramRegistry fRegistry{"pairQC", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  EMPhotonEventCut fEMEventCut;
  V0PhotonCut fV0PhotonCut;
  bool isMC = false;
  int mRunNumber{0};
  float mBzT{0.f};
  PhotonLikeCuts mAltCuts{};
  pairutil::DedupConfig mDedupCfg{};

  SliceCache cache;
  Preslice<MyV0Photons> perCollisionPCM = aod::v0photonkf::pmeventId;

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
  using FilteredMyMCCollisions = soa::Filtered<MyCollisionsMC>;

  void init(InitContext& context)
  {
    isMC = context.mOptions.get<bool>("processMC");
    defineEMEventCut();
    definePCMCut();
    mAltCuts.maxDca = crosspair.cfgAltMaxDca.value;
    mAltCuts.minR = crosspair.cfgAltMinR.value;
    mAltCuts.maxR = crosspair.cfgAltMaxR.value;
    mAltCuts.minCosPA = crosspair.cfgAltMinCosPA.value;
    mAltCuts.maxMee = crosspair.cfgAltMaxMee.value;
    mAltCuts.maxDcaXY = crosspair.cfgAltMaxDcaXY.value;
    mAltCuts.maxDcaZ = crosspair.cfgAltMaxDcaZ.value;
    mDedupCfg.maxLegDEta = dedup.cfgDupMaxLegDEta.value;
    mDedupCfg.maxLegDPhi = dedup.cfgDupMaxLegDPhi.value;
    mDedupCfg.maxLegPtAsym = dedup.cfgDupMaxLegPtAsym.value;
    mDedupCfg.maxLegDeDxAsym = dedup.cfgDupMaxLegDeDxAsym.value;
    mDedupCfg.maxDVtx3D = dedup.cfgDupMaxDVtx3D.value;
    mDedupCfg.requireBothLegs = dedup.cfgDupRequireBothLegs.value;
    ccdb->setURL(cfgCcdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    fRegistry.add("Photon/hPtEtaPhi", "selected V0 photons;p_{T} (GeV/c);#eta;#varphi (rad)", kTH3F, {{100, 0.f, 2.f}, {40, -0.8f, 0.8f}, {72, 0.f, 6.2832f}}, true); // o2-linter: disable=magic-number (axis definition)
    fRegistry.add("Photon/hNPhotonsPerEvent", "selected V0 photons per event;N_{#gamma};events", kTH1F, {{21, -0.5f, 20.5f}}, true);
    fRegistry.add("Photon/hV0LegClass", "leg classes of the selected V0s;;V0s", kTH1F, {makeAxisEnum(kNV0Classes, kV0ClassTitle)}, true);
    fRegistry.add("Photon/hV0LegClass_vs_Pt", "leg classes of the selected V0s vs p_{T};V0 leg class;p_{T} (GeV/c)", kTH2F, {makeAxisEnum(kNV0Classes, kV0ClassTitle), {40, 0.f, 2.f}}, true);
    fRegistry.add("Photon/hV0LegClass_vs_R", "leg classes of the selected V0s vs R_{conv};V0 leg class;R_{conv} (cm)", kTH2F, {makeAxisEnum(kNV0Classes, kV0ClassTitle), {50, 0.f, 100.f}}, true);
    addAnalyticV0Histograms();
    addPairQAHistograms();
    addDedupHistograms();
    addCrossPairHistograms();
    if (isMC) {
      addPairTruthHistograms(); // WP-C
      addAncestryHistograms();  // WP-E
    }
    LOGF(info, "pairQC: isMC = %d, pair QA below q_inv = %.2f, cross-pair QA %s (below %.3f), analytic V0 QA %s, ancestry %s, dedup %s",
         isMC, qaflags.cfgMaxQinvForQA.value, qaflags.cfgDoCrossPairQA.value ? "on" : "off", crosspair.cfgAltMaxQinv.value, qaflags.cfgDoAnalyticV0QA.value ? "on" : "off", qaflags.cfgDoAncestry.value ? "on" : "off", dedup.cfgDoDedup.value ? "ON" : "off");
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }
    mRunNumber = collision.runNumber();
    if (cfgBzOverrideT.value > -100.f) { // o2-linter: disable=magic-number (override sentinel)
      mBzT = cfgBzOverrideT.value;
      return;
    }
    auto grpmag = ccdb->getForRun<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", mRunNumber);
    mBzT = 0.1f * static_cast<float>(grpmag->getNominalL3Field()); // o2-linter: disable=magic-number (kGauss -> Tesla)
    LOGF(info, "pairQC: run %d, Bz = %.2f T", mRunNumber, mBzT);
  }

  void defineEMEventCut()
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

  void definePCMCut()
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

  void addAnalyticV0Histograms()
  {
    if (!qaflags.cfgDoAnalyticV0QA.value) {
      return;
    }
    const AxisSpec axR{50, 0.f, 100.f, "R_{conv}^{builder} (cm)"};
    const AxisSpec axTruth{3, -1.5f, 1.5f, "-1 = data / no label, 0 = fake V0, 1 = true photon"};
    fRegistry.add("Photon/AnalyticV0/hOk", "analytic rebuild succeeded;ok;photons", kTH1F, {{2, -0.5f, 1.5f}}, true);
    fRegistry.add("Photon/AnalyticV0/hdR_R_Truth", "R_{ana} - R_{builder};R_{ana} - R_{builder} (cm);R_{conv}^{builder} (cm);truth", kTH3F, {{200, -10.f, 10.f}, axR, axTruth}, true);
    fRegistry.add("Photon/AnalyticV0/hdZ_R_Truth", "z_{ana} - z_{builder};z_{ana} - z_{builder} (cm);R_{conv}^{builder} (cm);truth", kTH3F, {{200, -10.f, 10.f}, axR, axTruth}, true);
    fRegistry.add("Photon/AnalyticV0/hdPhiConv_R_Truth", "#varphi_{ana} - #varphi_{builder} of the conversion point;#Delta#varphi (rad);R_{conv}^{builder} (cm);truth", kTH3F, {{200, -0.2f, 0.2f}, axR, axTruth}, true);
    fRegistry.add("Photon/AnalyticV0/hdPtRel_Pt_Truth", "(p_{T,ana} - p_{T,builder}) / p_{T,builder};relative p_{T} difference;p_{T}^{builder} (GeV/c);truth", kTH3F, {{200, -0.5f, 0.5f}, {40, 0.f, 2.f}, axTruth}, true);
    fRegistry.add("Photon/AnalyticV0/hdCosPA_Truth", "cosPA_{ana} - cosPA_{builder};#Delta cosPA;truth", kTH2F, {{200, -0.02f, 0.02f}, axTruth}, true);
    fRegistry.add("Photon/AnalyticV0/hPcaAna_PcaBuilder_Truth", "PCA of the two legs;PCA_{ana} (cm);PCA_{builder} (cm);truth", kTH3F, {{60, 0.f, 6.f}, {60, 0.f, 6.f}, axTruth}, true);
    fRegistry.add("Photon/AnalyticV0/hScore_R_Truth", "builder score of the own analytic V0;score;R_{conv}^{builder} (cm);truth", kTH3F, {{200, 0.f, 30.f}, axR, axTruth}, true);
    fRegistry.add("Photon/AnalyticV0/hPhotonLike_Truth", "own analytic V0 passes the cfgAlt* cuts;photon-like;truth", kTH2F, {{2, -0.5f, 1.5f}, axTruth}, true);
  }

  void addPairQAHistograms()
  {
    const AxisSpec axDrCos{100, 0.f, 200.f, "|#vec{r}_{1}-#vec{r}_{2}| / cos(#alpha/2) (cm)"};
    const AxisSpec axOpa{100, 0.f, 0.5f, "opening angle of the conversion points (rad)"};
    const AxisSpec axNTimed{5, -0.5f, 4.5f, "legs with a time (ITS, TRD or TOF) among the four"};
    const AxisSpec axNForeign{5, -0.5f, 4.5f, "legs whose collisionId differs from their V0"};
    const AxisSpec axMaxDcaZ{50, 0.f, 10.f, "max |DCAz| of the two V0 photons to the PV (cm)"};
    const AxisSpec axMaxDcaXY{30, 0.f, 3.f, "max |DCAxy| of the two V0 photons to the PV (cm)"};
    const std::string base = "Pair/QA/Before/";
    fRegistry.add((base + "hDEtaDPhi_Qinv").c_str(), "pair (#Delta#eta,#Delta#phi,q_{inv});#Delta#eta;#Delta#phi (rad);q_{inv} (GeV/c)", kTH3F, {axisDeltaEta, axisDeltaPhi, axisQinv}, true);
    fRegistry.add((base + "hDeltaRVsQinv").c_str(), "|R_{1}-R_{2}| vs q_{inv};q_{inv} (GeV/c);|R_{1}-R_{2}| (cm)", kTH2F, {axisQinv, axisDeltaR}, true);
    fRegistry.add((base + "hDeltaZVsQinv").c_str(), "#Delta z vs q_{inv};q_{inv} (GeV/c);#Delta z (cm)", kTH2F, {axisQinv, axisDeltaZ}, true);
    fRegistry.add((base + "hDeltaR3DVsQinv").c_str(), "|#vec{r}_{1}-#vec{r}_{2}| vs q_{inv};q_{inv} (GeV/c);|#vec{r}_{1}-#vec{r}_{2}| (cm)", kTH2F, {axisQinv, axisDeltaR3D}, true);
    fRegistry.add((base + "hDrOverCosOAVsQinv").c_str(), "dr/cosOA vs q_{inv};q_{inv} (GeV/c);|#vec{r}_{1}-#vec{r}_{2}| / cos(#alpha/2) (cm)", kTH2F, {axisQinv, axDrCos}, true);
    fRegistry.add((base + "hOpeningAngleVsQinv").c_str(), "opening angle of the conversion points vs q_{inv};q_{inv} (GeV/c);#alpha (rad)", kTH2F, {axisQinv, axOpa}, true);
    fRegistry.add((base + "hQinv_NTimed_NForeign").c_str(), "pairs by leg collision association;q_{inv} (GeV/c);legs with a time;legs with foreign collisionId", kTH3F, {axisQinv, axNTimed, axNForeign}, true);
    fRegistry.add((base + "hQinv_MaxDcaZ_NTimed").c_str(), "pairs by V0 pointing to the primary vertex, z;q_{inv} (GeV/c);max |DCAz| (cm);legs with a time", kTH3F, {axisQinv, axMaxDcaZ, axNTimed}, true);
    fRegistry.add((base + "hQinv_MaxDcaXY_NTimed").c_str(), "pairs by V0 pointing to the primary vertex, xy;q_{inv} (GeV/c);max |DCAxy| (cm);legs with a time", kTH3F, {axisQinv, axMaxDcaXY, axNTimed}, true);
    fRegistry.add((base + "hRin_Rout_OuterSofter_Qinv").c_str(), "ordered conversion radii;R_{in} (cm);R_{out} (cm);outer photon softer (0/1);q_{inv} (GeV/c)", kTHnSparseF, {{50, 0.f, 100.f}, {50, 0.f, 100.f}, {2, -0.5f, 1.5f}, axisQinvCoarse}, true);
    fRegistry.add((base + "hLegStart_dNFindable_dNits_Qinv").c_str(), "leg start-point asymmetry per V0;|#Delta N_{findable}^{TPC}| (pos-neg);|#Delta N_{cls}^{ITS}|;q_{inv} (GeV/c)", kTH3F, {{40, -0.5f, 39.5f}, {8, -0.5f, 7.5f}, axisQinv}, true);
    fRegistry.add((base + "hV0_Eta_Phi_R_Qinv").c_str(), "eta, phi and conversion radius of each V0 of a pair;#eta;#varphi (rad);R_{conv} (cm);q_{inv} (GeV/c)", kTHnSparseF, {{20, -1.f, 1.f}, {72, 0.f, o2::constants::math::TwoPI}, {25, 0.f, 100.f}, axisQinvCoarse}, true);
    fRegistry.add((base + "hPhi_lowerPtV0").c_str(), "azimuthal angle of the lower-p_{T} V0;#phi (rad);pairs", kTH1F, {{90, 0.f, o2::constants::math::TwoPI}}, true);
    const AxisSpec axPairEta{20, -1.f, 1.f, "#eta_{pair}"};
    const AxisSpec axPairPhi{72, 0.f, o2::constants::math::TwoPI, "#varphi_{pair} (rad)"};
    const AxisSpec axDeltaRAng{100, 0.f, 0.5f, "#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#varphi^{2}}"};
    fRegistry.add((base + "hDEta_DPhi_PairEta_Qinv").c_str(), "#Delta#eta, #Delta#varphi vs pair #eta;#Delta#eta;#Delta#phi (rad);#eta_{pair};q_{inv} (GeV/c)", kTHnSparseF, {axisDeltaEta, axisDeltaPhi, axPairEta, axisQinvCoarse}, true);
    fRegistry.add((base + "hDEta_DPhi_PairPhi_Qinv").c_str(), "#Delta#eta, #Delta#varphi vs pair #varphi;#Delta#eta;#Delta#phi (rad);#varphi_{pair} (rad);q_{inv} (GeV/c)", kTHnSparseF, {axisDeltaEta, axisDeltaPhi, axPairPhi, axisQinvCoarse}, true);
    fRegistry.add((base + "hDPhi_Phi1_Phi2_Qinv").c_str(), "#Delta#varphi vs the azimuth of each V0;#Delta#phi (rad);#varphi_{1} (rad);#varphi_{2} (rad);q_{inv} (GeV/c)", kTHnSparseF, {axisDeltaPhi, axPairPhi, axPairPhi, axisQinvCoarse}, true);
    fRegistry.add((base + "hDEta_Eta1_Eta2_Qinv").c_str(), "#Delta#eta vs the pseudorapidity of each V0;#Delta#eta;#eta_{1};#eta_{2};q_{inv} (GeV/c)", kTHnSparseF, {axisDeltaEta, axPairEta, axPairEta, axisQinvCoarse}, true);
    fRegistry.add((base + "hDeltaRAngVsQinv").c_str(), "#DeltaR of the pair vs q_{inv};q_{inv} (GeV/c);#DeltaR", kTH2F, {axisQinv, axDeltaRAng}, true);
    fRegistry.add((base + "hDeltaRAng_PairEta_PairPhi").c_str(), "#DeltaR of the pair vs where it sits;#DeltaR;#eta_{pair};#varphi_{pair} (rad)", kTH3F, {axDeltaRAng, axPairEta, axPairPhi}, true);
    fRegistry.add((base + "hV0Class1_V0Class2_Qinv").c_str(), "leg classes of the two V0s of a pair (lower class first);V0 leg class 1;V0 leg class 2;q_{inv} (GeV/c)", kTH3F, {makeAxisEnum(kNV0Classes, kV0ClassTitle), makeAxisEnum(kNV0Classes, kV0ClassTitle), axisQinv}, true);
    fRegistry.add((base + "hNLegsITSTPC_NLegsTPConly_Qinv").c_str(), "legs of the four by class;N_{ITS-TPC legs} (0-4);N_{TPC-only legs} (0-4);q_{inv} (GeV/c)", kTH3F, {{5, -0.5f, 4.5f}, {5, -0.5f, 4.5f}, axisQinv}, true);
    fRegistry.addClone("Pair/QA/Before/", "Pair/QA/AfterPairCuts/");

    if (qaflags.cfgFillLegSimilarity.value) {
      fRegistry.add("Pair/LegSimilarity/hdEta_dPhi_ptRatio_dNsig_Qinv", "like-sign leg similarity;#Delta#eta(LS legs);#Delta#varphi(LS legs) (rad);|p_{T,1}-p_{T,2}|/(p_{T,1}+p_{T,2});|#Delta n#sigma_{e}^{TPC}|;q_{inv} (GeV/c)",
                    kTHnSparseF, {{100, -0.1f, 0.1f}, {100, -0.1f, 0.1f}, {50, 0.f, 1.f}, {40, 0.f, 8.f}, axisQinvCoarse}, true); // o2-linter: disable=magic-number (axis definition)
    }
  }

  void addDedupHistograms()
  {
    const AxisSpec axDEta{100, 0.f, 0.5f, "|#Delta#eta| (candidate pair)"};
    const AxisSpec axDPhi{100, 0.f, 0.5f, "|#Delta#varphi| (rad)"};
    const AxisSpec axPtAsym{50, 0.f, 1.f, "|p_{T,1}-p_{T,2}|/(p_{T,1}+p_{T,2})"};
    const AxisSpec axDR{100, 0.f, 20.f, "|R_{1}-R_{2}| (cm)"};
    fRegistry.add("Dedup/hdEta_dPhi_ptAsym", "distance between two V0 candidates of the same collision;|#Delta#eta|;|#Delta#varphi| (rad);p_{T} asymmetry", kTH3F, {axDEta, axDPhi, axPtAsym}, true);
    fRegistry.add("Dedup/hdR_ptAsym", "distance between two V0 candidates of the same collision;|R_{1}-R_{2}| (cm);p_{T} asymmetry", kTH2F, {axDR, axPtAsym}, true);
    fRegistry.add("Dedup/hNDupPartners_vs_Pt", "duplicate partners per candidate (configured tolerances);N_{partners};p_{T,#gamma} (GeV/c)", kTH2F, {{6, -0.5f, 5.5f}, {50, 0.f, 2.f}}, true);
    fRegistry.add("Dedup/hNCandBeforeAfter", "selected candidates per collision;0 = before dedup, 1 = after dedup;N_{cand}", kTH2F, {{2, -0.5f, 1.5f}, {21, -0.5f, 20.5f}}, true);
    const AxisSpec axLegTruth{5, -0.5f, 4.5f, "0=diffReco/diffMC, 1=diffReco/sameMC, 2=sameReco/sameMC, 3=sameReco/diffMC, 4=noMC"};
    fRegistry.add("Dedup/hLeg_dEta_dPhi_Truth", "same-charge legs of two candidates;|#Delta#eta|;|#Delta#varphi| (rad);track-level truth", kTH3F, {{50, 0.f, 0.05f}, {50, 0.f, 0.05f}, axLegTruth}, true);
    fRegistry.add("Dedup/hLeg_ptAsym_dEdxAsym_Truth", "same-charge legs of two candidates;p_{T} asym.;dE/dx asym.;track-level truth", kTH3F, {{50, 0.f, 0.5f}, {50, 0.f, 0.5f}, axLegTruth}, true);
    fRegistry.add("Dedup/hLeg_fShared_Truth", "same-charge legs of two candidates;max f_{shared};track-level truth", kTH2F, {{20, 0.f, 1.f}, axLegTruth}, true);
    if (!isMC) {
      return;
    }
    const AxisSpec axDup = makeAxisEnum(static_cast<int>(mcutil::DupClass::NClasses), "dup class (0 unknown, 1 same MC photon, 2 shared MC leg, 3 cross-built, 4 distinct true, 5 involves fake)");
    const AxisSpec axFlag{2, -0.5f, 1.5f, "flagged by the criterion"};
    fRegistry.add("Dedup/MC/hNCandPerMcPhoton", "V0 candidates per MC photon;N_{cand};MC photons", kTH1F, {{6, 0.5f, 6.5f}}, true);
    fRegistry.add("Dedup/MC/hNCandPerMcPhoton_vs_Pt", "V0 candidates per MC photon vs p_{T};N_{cand};p_{T,#gamma}^{reco} (GeV/c)", kTH2F, {{6, 0.5f, 6.5f}, {50, 0.f, 2.f}}, true);
    fRegistry.add("Dedup/MC/hDupClass", "MC class of each candidate pair of a collision;dup class;pairs", kTH1F, {axDup}, true);
    fRegistry.add("Dedup/MC/hDupClass_vs_Flagged", "MC class vs dedup decision;dup class;flagged", kTH2F, {axDup, axFlag}, true);
    fRegistry.add("Dedup/MC/hDupClass_dEta_dPhi", "candidate distance by MC class;dup class;|#Delta#eta|;|#Delta#varphi| (rad)", kTH3F, {axDup, axDEta, axDPhi}, true);
    fRegistry.add("Dedup/MC/hDupClass_ptAsym_dR", "candidate distance by MC class;dup class;p_{T} asymmetry;|R_{1}-R_{2}| (cm)", kTH3F, {axDup, axPtAsym, axDR}, true);
    fRegistry.add("Dedup/MC/hSamePhotonDuplicateType", "same MC photon: split pattern;0 = no split reco leg, 1 = split e^{+}, 2 = split e^{-}, 3 = both;pairs", kTH1F, {{4, -0.5f, 3.5f}}, true);
    fRegistry.add("Dedup/MC/hSamePhotonDuplicateType_vs_Flagged", "same MC photon: split pattern vs dedup decision;split pattern;flagged", kTH2F, {{4, -0.5f, 3.5f}, axFlag}, true);
    fRegistry.add("Dedup/MC/hSplitLegPairClass", "pair class for pairs containing a split reco leg;dup class;pairs", kTH1F, {axDup}, true);
    fRegistry.add("Dedup/MC/hSplitLegV0Type", "V0 truth of both candidates with a split leg;0 = same true photon, 1 = one true + one fake, 2 = both fake, 3 = other;pairs", kTH1F, {{4, -0.5f, 3.5f}}, true);
    fRegistry.add("Dedup/MC/hSplitLeg_ptAsym_dEdxAsym_DupClass", "split legs: distance vs pair class;p_{T} asym.;dE/dx asym.;dup class", kTH3F, {{50, 0.f, 0.5f}, {50, 0.f, 0.5f}, axDup}, true);
    fRegistry.add("Dedup/MC/hSplitLeg_fShared_DupClass", "split legs: shared clusters vs pair class;max f_{shared};dup class", kTH2F, {{20, 0.f, 1.f}, axDup}, true);
    fRegistry.add("Dedup/MC/hDupTrackClassMatrix", "true duplicates: ITS-TPC legs A vs B;n_{ITS-TPC}(A);n_{ITS-TPC}(B)", kTH2F, {{3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}}, true);
    fRegistry.add("Dedup/MC/hDupDeltaRconv", "true duplicates: radius difference;|R_{1}-R_{2}| (cm);pairs", kTH1F, {{200, 0.f, 20.f}}, true);
    fRegistry.add("Dedup/MC/hChoiceOutcome", "dedup choice for flagged pairs;0 true kept/fake dropped, 1 fake kept/true dropped, 2 both true same photon, 3 both fake, 4 both true distinct;flagged pairs", kTH1F, {{5, -0.5f, 4.5f}}, true);
    fRegistry.add("Dedup/MC/hLegTrackTruthClass", "track id vs MC id for same-charge legs;0 diff reco/diff MC, 1 diff reco/same MC, 2 same reco/same MC, 3 same reco/diff MC, 4 no MC;leg pairs", kTH1F, {{5, -0.5f, 4.5f}}, true);
    fRegistry.add("Dedup/MC/hLegTrackTruthClass_vs_Flagged", "track-level truth vs leg flagged;class;leg flagged", kTH2F, {{5, -0.5f, 4.5f}, axFlag}, true);
    fRegistry.add("Pair/MC/hTruthType_Qinv_DedupFlag", "would the dedup remove this pair?;truth type;q_{inv} (GeV/c);0 = stays, 1 = removed", kTH3F, {axisTruthType, axisQinv, axFlag}, true);
  }

  void addCrossPairHistograms()
  {
    auto h = fRegistry.add<TH1>("Pair/CrossPair/hVetoCounter", "crossed-pair veto (DCAFitter, cfgCrossVetoMode);;pairs", kTH1D, {{2, -0.5f, 1.5f}}, true);
    h->GetXaxis()->SetBinLabel(1, "evaluated");
    h->GetXaxis()->SetBinLabel(2, "vetoed");
    if (!qaflags.cfgDoCrossPairQA.value) {
      return;
    }
    const AxisSpec axNAlt{3, -0.5f, 2.5f, "N photon-like alternative pairings"};
    const AxisSpec axDca{60, 0.f, 12.f, "max DCA of the two alternative pairings (cm)"};
    const AxisSpec axDcaXY{60, 0.f, 12.f, "max DCAxy to the vertex of the two alternative pairings (cm)"};
    const AxisSpec axDcaZ{60, 0.f, 12.f, "max DCAz to the vertex of the two alternative pairings (cm)"};
    for (const auto& m : {std::string("Fit"), std::string("Ana")}) {
      fRegistry.add(("Pair/CrossPair/hQinv_NAltValid_MaxAltDca_" + m).c_str(), ("alternative pairings of the four legs (" + m + ");q_{inv} (GeV/c);N photon-like;max DCA (cm)").c_str(), kTH3F, {axisQinv, axNAlt, axDca}, true);
      fRegistry.add(("Pair/CrossPair/hQinv_NAltValid_MaxAltDcaXY_" + m).c_str(), ("alternative pairings, xy (" + m + ");q_{inv} (GeV/c);N photon-like;max DCAxy (cm)").c_str(), kTH3F, {axisQinv, axNAlt, axDcaXY}, true);
      fRegistry.add(("Pair/CrossPair/hQinv_NAltValid_MaxAltDcaZ_" + m).c_str(), ("alternative pairings, z (" + m + ");q_{inv} (GeV/c);N photon-like;max DCAz (cm)").c_str(), kTH3F, {axisQinv, axNAlt, axDcaZ}, true);
    }
    fRegistry.add("Pair/CrossPair/hAltDca_Fit_vs_Ana", "alternative pairing: DCA of the legs, DCAFitter vs analytic;DCA_{fit} (cm);DCA_{ana} (cm)", kTH2F, {{120, 0.f, 12.f}, {120, 0.f, 12.f}}, true);
    fRegistry.add("Pair/CrossPair/hAltCosPA_Fit_vs_Ana", "alternative pairing: cosPA, DCAFitter vs analytic;cosPA_{fit};cosPA_{ana}", kTH2F, {{100, 0.9f, 1.f}, {100, 0.9f, 1.f}}, true);
    fRegistry.add("Pair/CrossPair/hNAltValid_Fit_vs_Ana", "N photon-like alternatives, DCAFitter vs analytic;N_{fit};N_{ana}", kTH2F, {axNAlt, axNAlt}, true);
    fRegistry.add("Pair/CrossPair/hMeeRatio_dR_Qinv", "crossed-hypothesis mass;min(m_{ee}^{cross})/q_{inv};|R_{1}-R_{2}| (cm);q_{inv} (GeV/c)", kTH3F, {{100, 0.f, 2.f}, {80, 0.f, 80.f}, axisQinv}, true);
    fRegistry.add("Pair/CrossPair/hRepairDeltaMee_OwnMee_Qinv", "re-pairing test;max(m_{ee}^{own}) - max(m_{ee}^{re-paired}) (GeV/c^{2});max(m_{ee}^{own}) (GeV/c^{2});q_{inv} (GeV/c)", kTH3F, {{120, -0.06f, 0.06f}, {60, 0.f, 0.06f}, axisQinv}, true);
    fRegistry.add("Pair/CrossPair/hMee1_Mee2_M4_Qinv", "pair Dalitz plane;m_{ee}(#gamma_{1}) (GeV/c^{2});m_{ee}(#gamma_{2}) (GeV/c^{2});M(4 legs) (GeV/c^{2});q_{inv} (GeV/c)", kTHnSparseF, {{30, 0.f, 0.06f}, {30, 0.f, 0.06f}, {75, 0.f, 0.3f}, axisQinvCoarse}, true);
    if (!isMC) {
      return;
    }
    for (const auto& m : {std::string("Fit"), std::string("Ana")}) {
      fRegistry.add(("Pair/MC/CrossPair/hQinv_NAltValid_MaxAltDca_" + m + "_Type").c_str(), ("alternative pairings per truth type (" + m + ");q_{inv} (GeV/c);N photon-like;max DCA (cm);truth type").c_str(), kTHnSparseF, {axisQinv, axNAlt, axDca, axisTruthType}, true);
      fRegistry.add(("Pair/MC/CrossPair/hQinv_NAltValid_MaxAltDcaZ_" + m + "_Type").c_str(), ("alternative pairings per truth type, z (" + m + ");q_{inv} (GeV/c);N photon-like;max DCAz (cm);truth type").c_str(), kTHnSparseF, {axisQinv, axNAlt, axDcaZ, axisTruthType}, true);
    }
    fRegistry.add("Pair/MC/CrossPair/hMeeRatio_Qinv_Type", "crossed-pair hypothesis per truth type;min(m_{ee}^{cross})/q_{inv};q_{inv} (GeV/c);truth type", kTH3F, {{100, 0.f, 2.f}, axisQinv, axisTruthType}, true);
    fRegistry.add("Pair/MC/CrossPair/hMee1_Mee2_M4_Qinv_Type", "pair Dalitz plane per truth type;m_{ee}(#gamma_{1});m_{ee}(#gamma_{2});M(4 legs);q_{inv} (GeV/c);truth type", kTHnSparseF, {{30, 0.f, 0.06f}, {30, 0.f, 0.06f}, {75, 0.f, 0.3f}, axisQinvCoarse, axisTruthType}, true);
    fRegistry.add("Pair/MC/CrossPair/hQinv_MaxDcaZ_Type", "V0 pointing per truth type, z;q_{inv} (GeV/c);max |DCAz| (cm);truth type", kTH3F, {axisQinv, {50, 0.f, 10.f}, axisTruthType}, true);
    fRegistry.add("Pair/MC/CrossPair/hQinv_MaxDcaXY_Type", "V0 pointing per truth type, xy;q_{inv} (GeV/c);max |DCAxy| (cm);truth type", kTH3F, {axisQinv, {30, 0.f, 3.f}, axisTruthType}, true);
    const AxisSpec axEmu{4, -0.5f, 3.5f, "0 = not evaluated, 1 = kept, 2 = swapped, 3 = cross not photon-like"};
    fRegistry.add("Pair/MC/CrossPair/hEmulation_vs_Truth", "cross-pairing emulation (analytic score) vs truth;q_{inv} (GeV/c);truth type;outcome", kTH3F, {axisQinv, axisTruthType, axEmu}, true);
    fRegistry.add("Pair/MC/CrossPair/hEmuDeltaScore_vs_Truth", "S_{cross} - S_{true} vs truth;q_{inv} (GeV/c);truth type;#DeltaS", kTH3F, {axisQinv, axisTruthType, {200, -20.f, 20.f}}, true);
    fRegistry.add("Pair/MC/CrossPair/hFakeSubtype_Emu_Qinv", "fake subtype vs emulation outcome;fake subtype (0 ordinary, 1 full cross, 2 half cross);outcome;q_{inv} (GeV/c)", kTH3F, {{3, -0.5f, 2.5f}, axEmu, axisQinv}, true);
  }

  void addPairTruthHistograms()
  {
    fRegistry.add("Pair/MC/hTruthTypeVsQinv", "truth type vs q_{inv};q_{inv} (GeV/c);truth type", kTH2F, {axisQinv, axisTruthType}, true);
    fRegistry.add("Pair/MC/hTruthTypeVsKt", "truth type vs k_{T};k_{T} (GeV/c);truth type", kTH2F, {axisKt, axisTruthType}, true);
    fRegistry.add("Pair/MC/hFakeSubtypeVsQinv", "fake pair subtype vs q_{inv};q_{inv} (GeV/c);0=ordinary, 1=full cross, 2=half cross", kTH2F, {axisQinv, {3, -0.5f, 2.5f}}, true);
    fRegistry.add("Pair/MC/hV0Class1_V0Class2_Qinv_Type", "leg classes of the two V0s of a pair vs truth type;V0 leg class 1;V0 leg class 2;q_{inv} (GeV/c);truth type", kTHnSparseF, {makeAxisEnum(kNV0Classes, kV0ClassTitle), makeAxisEnum(kNV0Classes, kV0ClassTitle), axisQinv, axisTruthType}, true);
    fRegistry.add("Pair/MC/hV0Class_IsTrue_Qinv", "per V0: leg class vs V0 truth;V0 leg class;V0 is a true photon;q_{inv} (GeV/c)", kTH3F, {makeAxisEnum(kNV0Classes, kV0ClassTitle), {2, -0.5f, 1.5f}, axisQinv}, true);
    static constexpr std::array<std::string_view, 6> kTypes = {"TrueTrueDistinct/", "TrueTrueSamePhoton/", "SharedMcLeg/", "TrueFake/", "FakeFake/", "Pi0Daughters/"};
    for (const auto& label : kTypes) {
      const std::string base = std::string("Pair/MC/") + std::string(label);
      fRegistry.add((base + "hDeltaRVsQinv").c_str(), "|R_{1}-R_{2}| vs q_{inv};q_{inv} (GeV/c);|R_{1}-R_{2}| (cm)", kTH2F, {axisQinv, axisDeltaR}, true);
      fRegistry.add((base + "hDeltaZVsQinv").c_str(), "#Delta z vs q_{inv};q_{inv} (GeV/c);#Delta z (cm)", kTH2F, {axisQinv, axisDeltaZ}, true);
      fRegistry.add((base + "hDeltaR3DVsQinv").c_str(), "|#vec{r}_{1}-#vec{r}_{2}| vs q_{inv};q_{inv} (GeV/c);|#vec{r}_{1}-#vec{r}_{2}| (cm)", kTH2F, {axisQinv, axisDeltaR3D}, true);
      fRegistry.add((base + "hDEtaDPhi").c_str(), "#Delta#eta vs #Delta#phi;#Delta#eta;#Delta#phi (rad)", kTH2F, {axisDeltaEta, axisDeltaPhi}, true);
      if (qaflags.cfgFillTypeSparses.value) {
        fRegistry.add((base + "hDEtaDPhi_Qinv").c_str(), "#Delta#eta, #Delta#phi, q_{inv};#Delta#eta;#Delta#phi (rad);q_{inv} (GeV/c)", kTH3F, {axisDeltaEta, axisDeltaPhi, axisQinv}, true);
        fRegistry.add((base + "hDeltaR_DeltaZ_Qinv").c_str(), "|R_{1}-R_{2}|, #Delta z, q_{inv};|R_{1}-R_{2}| (cm);#Delta z (cm);q_{inv} (GeV/c)", kTH3F, {axisDeltaR, axisDeltaZ, axisQinv}, true);
      }
    }
    fRegistry.add("Pair/MC/hRreco_minus_Rmc", "conversion radius resolution by V0 truth;R_{reco} - R_{MC}(leg) (cm);R_{MC}(leg) (cm);V0 is a true photon", kTH3F, {{120, -30.f, 30.f}, {50, 0.f, 100.f}, {2, -0.5f, 1.5f}}, true);
    fRegistry.add("Pair/MC/hPsiPair_PhiV_Qinv_Type_IsTrue", "per-photon conversion topology;|#psi_{pair}| (rad);#varphi_{V} (rad);q_{inv} (GeV/c);truth type;is true photon", kTHnSparseF, {{90, 0.f, o2::constants::math::PIHalf}, {90, 0.f, o2::constants::math::PI}, axisQinvCoarse, axisTruthType, {2, -0.5f, 1.5f}}, true);
    fRegistry.add("Pair/MC/hMcEventPattern_NOwn_Qinv_Type", "MC event of the four legs;0 all same, 1 two consistent V0 from different events, 2 a V0 mixes events;N legs from own event;q_{inv} (GeV/c);truth type", kTHnSparseF, {{3, -0.5f, 2.5f}, {5, -0.5f, 4.5f}, axisQinv, axisTruthType}, true);
  }

  void addAncestryHistograms()
  {
    if (!qaflags.cfgDoAncestry.value) {
      return;
    }
    const AxisSpec axAnc = makeAxisEnum(static_cast<int>(mcutil::AncestorKind::NClasses), "ancestor kind (0 none, 1 string/parton, 2 beam/nucleus, 3 #pi^{0}, 4 #eta, 5 other meson, 6 baryon, 7 photon, 8 electron, 9 other)");
    const AxisSpec axReco = makeAxisEnum(static_cast<int>(mcutil::V0RecoStatus::NClasses), "V0 reco status (0 no label, 1 no conv leg, 2 one conv leg, 3 cross-built, 4 true photon)");
    const AxisSpec axSplit{6, -0.5f, 5.5f, "split code (2 n_{consistent} + shared)"};
    const AxisSpec axSame{2, -0.5f, 1.5f, "0 = different origin particles, 1 = the same"};
    const AxisSpec axSrc = makeAxisEnum(kNPhotonSources, kPhotonSourceTitle);
    const AxisSpec axQtrue{60, 0.f, 0.3f, "q_{inv}^{true} (GeV/c)"};
    const std::string base = "Pair/MC/Ancestry/";

    fRegistry.add((base + "hCensus_Qinv_Type").c_str(), "census of the four legs: code = (legs from a conversion) + 5 (distinct parent photons);census code;q_{inv} (GeV/c);truth type", kTH3F, {{25, -0.5f, 24.5f}, axisQinvCoarse, axisTruthType}, true);
    fRegistry.add((base + "hSplit_Qinv_Type").c_str(), "split code of the four legs;split code;q_{inv} (GeV/c);truth type", kTH3F, {axSplit, axisQinvCoarse, axisTruthType}, true);

    fRegistry.add((base + "hQtrue_Qreco_Split").c_str(), "pairs with exactly two parent photons: true q_{inv} of those two vs reconstructed, by pair class;q_{inv}^{true} (GeV/c);q_{inv}^{reco} (GeV/c);split code", kTH3F, {axQtrue, axisQinvCoarse, axSplit}, true);
    fRegistry.add((base + "hAncKind_Qinv_Type").c_str(), "common ancestor of the two parent photons;common ancestor;q_{inv} (GeV/c);truth type", kTH3F, {axAnc, axisQinvCoarse, axisTruthType}, true);
    fRegistry.add((base + "hAnc1_Anc2_Same").c_str(), "origin particle of each parent photon (walk stops before beam/string);origin of photon A;origin of photon B;same origin particle?", kTH3F, {axAnc, axAnc, axSame}, true);
    fRegistry.add((base + "hSame_Qinv_Type").c_str(), "the two parent photons come from the same origin particle;same origin?;q_{inv} (GeV/c);truth type", kTH3F, {axSame, axisQinvCoarse, axisTruthType}, true);
    fRegistry.add((base + "hSrc1_Src2_Type").c_str(), "pairs with exactly two parent photons: source of each (lower code first);source of photon A;source of photon B;truth type", kTH3F, {axSrc, axSrc, axisTruthType}, true);
    fRegistry.add((base + "hTrueConv_d3D_Qinv_Type").c_str(), "true conversion points of the two parent photons;|#vec{r}_{A}-#vec{r}_{B}|^{true} (cm);q_{inv} (GeV/c);truth type", kTH3F, {{50, 0.f, 50.f}, axisQinvCoarse, axisTruthType}, true);
    fRegistry.add((base + "hTrueConv_dR_dZ_Type").c_str(), "true conversion points of the two parent photons;|R_{A}-R_{B}|^{true} (cm);|z_{A}-z_{B}|^{true} (cm);truth type", kTH3F, {{50, 0.f, 25.f}, {50, 0.f, 25.f}, axisTruthType}, true);
    fRegistry.add((base + "hRMeanTrue_DRTrue_Split").c_str(), "per pair: mean and difference of the two true conversion radii, by pair class;(R_{A}+R_{B})/2 (cm);|R_{A}-R_{B}| (cm);split code", kTH3F, {{100, 0.f, 100.f}, {50, 0.f, 25.f}, axSplit}, true);
    fRegistry.add((base + "hTrueConv_dTheta_Qtrue_Type").c_str(), "true opening angle and q of the two parent photons;#Delta#theta_{#gamma#gamma}^{true} (mrad);q_{inv}^{true} (GeV/c);truth type", kTH3F, {{50, 0.f, 50.f}, axQtrue, axisTruthType}, true);
    fRegistry.add((base + "hTrueConv_DeltaQ_Qtrue_Type").c_str(), "q_{inv}^{reco} - q_{inv}^{true} per pair;q_{inv}^{reco} - q_{inv}^{true} (GeV/c);q_{inv}^{true} (GeV/c);truth type", kTH3F, {{140, -0.05f, 0.30f}, axQtrue, axisTruthType}, true);

    fRegistry.add((base + "hEta_Phi_Reco").c_str(), "per V0 candidate: eta, phi and reco status;#eta;#varphi (rad);V0 reco status", kTH3F, {{20, -1.f, 1.f}, {72, 0.f, o2::constants::math::TwoPI}, axReco}, true);
    fRegistry.add((base + "hR_Reco_Qinv").c_str(), "per V0 candidate: conversion radius and reco status;R_{conv} (cm);V0 reco status;q_{inv} (GeV/c)", kTH3F, {{25, 0.f, 100.f}, axReco, axisQinvCoarse}, true);
    fRegistry.add((base + "hVx_Vy_Split").c_str(), "reconstructed conversion point of each V0, labelled by the class of its PAIR;V_{x} (cm);V_{y} (cm);split code", kTH3F, {{200, -100.f, 100.f}, {200, -100.f, 100.f}, axSplit}, true);
    fRegistry.add((base + "hVxTrue_VyTrue_Split").c_str(), "true conversion point (production vertex of the e^{+} leg), by pair class;V_{x}^{true} (cm);V_{y}^{true} (cm);split code", kTH3F, {{200, -100.f, 100.f}, {200, -100.f, 100.f}, axSplit}, true);
    fRegistry.add((base + "hRTruePos_RTrueNeg_Reco").c_str(), "true conversion radius of each leg, by reco status;R^{true}(e^{+}) (cm);R^{true}(e^{-}) (cm);V0 reco status", kTH3F, {{100, 0.f, 100.f}, {100, 0.f, 100.f}, axReco}, true);
    fRegistry.add((base + "hRTrueMean_RReco_Reco").c_str(), "mean true conversion radius of the two legs vs reconstructed V0 radius;(R^{true}(e^{+}) + R^{true}(e^{-}))/2 (cm);R^{reco} (cm);V0 reco status", kTH3F, {{100, 0.f, 100.f}, {100, 0.f, 100.f}, axReco}, true);
    fRegistry.add((base + "hSrcPos_SrcNeg_Reco").c_str(), "per V0 candidate: parent photon sources and reco status;source of the e^{+} parent;source of the e^{-} parent;V0 reco status", kTH3F, {axSrc, axSrc, axReco}, true);
    fRegistry.add((base + "hSrcPos_SrcNeg_Type").c_str(), "per V0 candidate: parent photon sources and the truth type of its pair;source of the e^{+} parent;source of the e^{-} parent;truth type", kTH3F, {axSrc, axSrc, axisTruthType}, true);
    fRegistry.add((base + "hAncPos_AncNeg_Reco").c_str(), "per V0: origin particle of each leg's parent photon;origin of the e^{+} parent;origin of the e^{-} parent;V0 reco status", kTH3F, {axAnc, axAnc, axReco}, true);
    fRegistry.add((base + "hAncSame_Reco").c_str(), "per V0: the two legs' parent photons have the same origin particle;same origin?;V0 reco status", kTH2F, {axSame, axReco}, true);
    std::vector<double> psiEdges;
    constexpr int kPsiFine = 60, kPsiCoarse = 13;
    for (int i = 0; i <= kPsiFine; ++i) {
      psiEdges.push_back(0.005 * i); // o2-linter: disable=magic-number (5 mrad steps)
    }
    for (int i = 1; i <= kPsiCoarse; ++i) {
      psiEdges.push_back(0.30 + 0.1 * i); // o2-linter: disable=magic-number (coarse part up to pi/2)
    }
    std::vector<double> chi2Edges;
    constexpr int kChi2Fine = 50, kChi2Coarse = 20;
    for (int i = 0; i <= kChi2Fine; ++i) {
      chi2Edges.push_back(0.2 * i); // o2-linter: disable=magic-number (0.2 steps to 10)
    }
    for (int i = 1; i <= kChi2Coarse; ++i) {
      chi2Edges.push_back(10. + 2.0 * i); // o2-linter: disable=magic-number (coarse part to 50)
    }
    fRegistry.add((base + "hPsiPair_Chi2_Split").c_str(), "per V0: the Run-2 cut plane, by pair class;|#psi_{pair}| (rad);#chi^{2}/ndf (KF);split code", kTH3F, {{psiEdges, "|#psi_{pair}| (rad)"}, {chi2Edges, "#chi^{2}/ndf (KF)"}, axSplit}, true);
    fRegistry.add((base + "hPsiPair_R_Split").c_str(), "per V0: |psi_pair| vs conversion radius, by pair class;|#psi_{pair}| (rad);R_{conv}^{reco} (cm);split code", kTH3F, {{psiEdges, "|#psi_{pair}| (rad)"}, {10, 0.f, 100.f}, axSplit}, true);
  }

  [[nodiscard]] inline bool passAsymmetryCut(float pt1, float pt2) const
  {
    if (ggpaircuts.cfgMaxAsymmetry.value < 0.f) {
      return true;
    }
    const float sum = pt1 + pt2;
    if (sum < KMinSigma) {
      return false;
    }
    return std::fabs(pt1 - pt2) / sum < ggpaircuts.cfgMaxAsymmetry.value;
  }

  [[nodiscard]] inline bool isInsideEllipse(float deta, float dphi) const
  {
    if (!ggpaircuts.cfgDoEllipseCut.value) {
      return false;
    }
    const float sE = ggpaircuts.cfgEllipseSigEta.value;
    const float sP = ggpaircuts.cfgEllipseSigPhi.value;
    if (sE < KMinSigma || sP < KMinSigma) {
      return false;
    }
    return (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP) < ggpaircuts.cfgEllipseR2.value;
  }

  [[nodiscard]] inline bool passRZCut(float deltaR, float deltaZ) const
  {
    if (ggpaircuts.cfgDoRCut.value && deltaR < ggpaircuts.cfgMinDeltaR.value) {
      return false;
    }
    if (ggpaircuts.cfgDoZCut.value && std::fabs(deltaZ) < ggpaircuts.cfgMinDeltaZ.value) {
      return false;
    }
    return true;
  }

  [[nodiscard]] bool passPairCuts(PairObs const& obs, PhotonCand const& a, PhotonCand const& b) const
  {
    if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA.value) {
      return false;
    }
    if (!passRZCut(obs.deltaR, obs.deltaZ)) {
      return false;
    }
    if (isInsideEllipse(obs.deta, obs.dphi)) {
      return false;
    }
    return std::max(std::fabs(a.fDcaZToPV), std::fabs(b.fDcaZToPV)) < ggpaircuts.cfgMaxDcaZToPV.value &&
           std::max(std::fabs(a.fDcaXYToPV), std::fabs(b.fDcaXYToPV)) < ggpaircuts.cfgMaxDcaXYToPV.value;
  }

  [[nodiscard]] PairObs buildPairObs(PhotonCand const& g1, PhotonCand const& g2) const
  {
    PairObs o;
    const float dx = g1.fVx - g2.fVx, dy = g1.fVy - g2.fVy, dz = g1.fVz - g2.fVz;
    o.deltaR = std::fabs(g1.rConv() - g2.rConv());
    o.deltaZ = dz;
    o.deltaR3D = std::sqrt(dx * dx + dy * dy + dz * dz);
    const ROOT::Math::XYZVector cp1(g1.fVx, g1.fVy, g1.fVz), cp2(g2.fVx, g2.fVy, g2.fVz);
    const float mag1 = std::sqrt(cp1.Mag2()), mag2 = std::sqrt(cp2.Mag2());
    if (mag1 < KMinMagnitude || mag2 < KMinMagnitude) {
      o.valid = false;
      return o;
    }
    const float cosPA = std::clamp(static_cast<float>(cp1.Dot(cp2) / (mag1 * mag2)), -1.f, 1.f);
    o.opa = std::acos(cosPA);
    const float cosHalf = std::cos(0.5f * o.opa);
    o.drOverCosOA = (std::fabs(cosHalf) < KMinCosine) ? 1e12f : (o.deltaR3D / cosHalf); // o2-linter: disable=magic-number (sentinel for a degenerate opening angle)
    o.v1 = ROOT::Math::PtEtaPhiMVector(g1.fPt, g1.fEta, g1.fPhi, 0.f);
    o.v2 = ROOT::Math::PtEtaPhiMVector(g2.fPt, g2.fEta, g2.fPhi, 0.f);
    o.deta = g1.fEta - g2.fEta;
    o.dphi = RecoDecay::constrainAngle(g1.fPhi - g2.fPhi, -o2::constants::math::PI);
    o.q = pairutil::computePairQ(o.v1, o.v2);
    return o;
  }

  [[nodiscard]] CrossObs computeCrossObs(PhotonCand const& a, PhotonCand const& b, float qinv) const
  {
    CrossObs c;
    if (qinv < crosspair.cfgAltMaxQinv.value) {
      const float bzkG = 10.f * mBzT; // o2-linter: disable=magic-number (Tesla -> kGauss)
      c.fit[0] = evaluateAltPairingFitter(a.legParam(0), b.legParam(1), mBzT);
      c.fit[1] = evaluateAltPairingFitter(b.legParam(0), a.legParam(1), mBzT);
      c.ana[0] = buildAnalyticV0(a.legParam(0), b.legParam(1), bzkG, qaflags.cfgScoreWeight.value);
      c.ana[1] = buildAnalyticV0(b.legParam(0), a.legParam(1), bzkG, qaflags.cfgScoreWeight.value);
      c.nAltValidFit = (fitIsPhotonLike(c.fit[0]) ? 1 : 0) + (fitIsPhotonLike(c.fit[1]) ? 1 : 0);
      c.nAltValidAna = (isPhotonLike(c.ana[0], mAltCuts) ? 1 : 0) + (isPhotonLike(c.ana[1], mAltCuts) ? 1 : 0);
      c.maxDcaFit = std::max(c.fit[0].dca, c.fit[1].dca);
      c.maxDcaXYFit = std::max(c.fit[0].dcaxy, c.fit[1].dcaxy);
      c.maxDcaZFit = std::max(c.fit[0].dcaz, c.fit[1].dcaz);
      c.maxDcaAna = std::max(c.ana[0].pca, c.ana[1].pca);
      c.maxDcaXYAna = std::max(c.ana[0].dcaxy, c.ana[1].dcaxy);
      c.maxDcaZAna = std::max(c.ana[0].dcaz, c.ana[1].dcaz);
      c.evaluated = true;
    }
    auto legVec = [](PhotonCand const& p, int i) {
      const auto& l = p.leg[static_cast<size_t>(i)];
      return ROOT::Math::PtEtaPhiMVector(l.pt, l.eta, l.phi, KMe);
    };
    c.mee[0] = static_cast<float>((legVec(a, 0) + legVec(b, 1)).M());
    c.mee[1] = static_cast<float>((legVec(b, 0) + legVec(a, 1)).M());
    if (qinv > 0.f) {
      c.meeOverQ = std::min(c.mee[0], c.mee[1]) / qinv;
    }
    c.ownMee[0] = static_cast<float>((legVec(a, 0) + legVec(a, 1)).M());
    c.ownMee[1] = static_cast<float>((legVec(b, 0) + legVec(b, 1)).M());
    c.ownMeeMax = std::max(c.ownMee[0], c.ownMee[1]);
    c.deltaMee = c.ownMeeMax - std::max(c.mee[0], c.mee[1]);
    c.m4 = static_cast<float>((legVec(a, 0) + legVec(a, 1) + legVec(b, 0) + legVec(b, 1)).M());
    return c;
  }

  [[nodiscard]] inline bool fitIsPhotonLike(AltPairing const& r) const
  {
    return r.fitted && r.dca < mAltCuts.maxDca && r.rxy > mAltCuts.minR && r.rxy < mAltCuts.maxR &&
           r.cospa > mAltCuts.minCosPA && r.mee < mAltCuts.maxMee && r.dcaxy < mAltCuts.maxDcaXY && r.dcaz < mAltCuts.maxDcaZ;
  }

  [[nodiscard]] inline bool passCrossPairVeto(CrossObs const& c) const
  {
    const int need = crosspair.cfgAltRequireBoth.value ? 2 : 1;
    const bool passGeom = c.nAltValidFit < need;
    const bool passDca = (crosspair.cfgCrossMinMaxAltDca.value <= 0.f) || (c.maxDcaFit > crosspair.cfgCrossMinMaxAltDca.value);
    switch (crosspair.cfgCrossVetoMode.value) {
      case 2: // o2-linter: disable=magic-number (veto mode)
        return passDca;
      case 3: // o2-linter: disable=magic-number (veto mode)
        return passGeom && passDca;
      default:
        return passGeom;
    }
  }
  struct Emulation {
    int outcome{0};
    float deltaScore{0.f};
  };
  [[nodiscard]] Emulation emulate(CrossObs const& c, PhotonCand const& a, PhotonCand const& b) const
  {
    Emulation e;
    if (!c.evaluated || !a.own.ok || !b.own.ok) {
      return e;
    }
    if (c.nAltValidAna < 2) { // o2-linter: disable=magic-number (both cross candidates must be photon-like)
      e.outcome = 3;          // o2-linter: disable=magic-number (outcome code)
      return e;
    }
    e.deltaScore = (c.ana[0].score + c.ana[1].score) - (a.own.score + b.own.score);
    e.outcome = (e.deltaScore < 0.f) ? 2 : 1; // o2-linter: disable=magic-number (outcome code)
    return e;
  }

  template <int step>
  void fillPairQA(PairObs const& obs, PhotonCand const& a, PhotonCand const& b)
  {
    constexpr auto dir = (step == 0) ? "Pair/QA/Before/" : "Pair/QA/AfterPairCuts/";
    const float qinv = obs.q.qinv;
    fRegistry.fill(HIST(dir) + HIST("hDEtaDPhi_Qinv"), obs.deta, obs.dphi, qinv);
    fRegistry.fill(HIST(dir) + HIST("hDeltaRVsQinv"), qinv, obs.deltaR);
    fRegistry.fill(HIST(dir) + HIST("hDeltaZVsQinv"), qinv, obs.deltaZ);
    fRegistry.fill(HIST(dir) + HIST("hDeltaR3DVsQinv"), qinv, obs.deltaR3D);
    fRegistry.fill(HIST(dir) + HIST("hDrOverCosOAVsQinv"), qinv, std::min(obs.drOverCosOA, 199.9f)); // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST(dir) + HIST("hOpeningAngleVsQinv"), qinv, obs.opa);
    const int nTimed = a.nTimedLegs() + b.nTimedLegs();
    fRegistry.fill(HIST(dir) + HIST("hQinv_NTimed_NForeign"), qinv, static_cast<float>(nTimed), static_cast<float>(a.nForeignLegs() + b.nForeignLegs()));
    fRegistry.fill(HIST(dir) + HIST("hQinv_MaxDcaZ_NTimed"), qinv, std::min(std::max(std::fabs(a.fDcaZToPV), std::fabs(b.fDcaZToPV)), 9.9f), static_cast<float>(nTimed));    // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST(dir) + HIST("hQinv_MaxDcaXY_NTimed"), qinv, std::min(std::max(std::fabs(a.fDcaXYToPV), std::fabs(b.fDcaXYToPV)), 2.9f), static_cast<float>(nTimed)); // o2-linter: disable=magic-number (clamp to the axis)
    const float ra = a.rConv(), rb = b.rConv();
    const bool aInner = ra <= rb;
    const float rIn = aInner ? ra : rb, rOut = aInner ? rb : ra;
    const bool outerSofter = aInner ? (b.fPt < a.fPt) : (a.fPt < b.fPt);
    const bool inSparseWindow = qinv < qaflags.cfgMaxQinvForSparses.value;
    if (inSparseWindow) {
      fRegistry.fill(HIST(dir) + HIST("hRin_Rout_OuterSofter_Qinv"), rIn, rOut, outerSofter ? 1.f : 0.f, qinv);
    }
    for (const PhotonCand& p : {std::cref(a), std::cref(b)}) {
      fRegistry.fill(HIST(dir) + HIST("hLegStart_dNFindable_dNits_Qinv"), std::fabs(p.leg[0].nClsFindable - p.leg[1].nClsFindable), std::fabs(p.leg[0].nClsITS - p.leg[1].nClsITS), qinv);
      if (inSparseWindow) {
        fRegistry.fill(HIST(dir) + HIST("hV0_Eta_Phi_R_Qinv"), p.fEta, p.fPhi, std::min(p.rConv(), 99.9f), qinv); // o2-linter: disable=magic-number (clamp to the axis)
      }
    }
    fRegistry.fill(HIST(dir) + HIST("hPhi_lowerPtV0"), (a.fPt < b.fPt) ? a.fPhi : b.fPhi);
    const float pairEta = 0.5f * (a.fEta + b.fEta);
    const float pairPhi = RecoDecay::constrainAngle(static_cast<float>((obs.v1 + obs.v2).Phi()), 0.f);
    const float deltaRAng = std::hypot(obs.deta, obs.dphi);
    if (inSparseWindow) {
      fRegistry.fill(HIST(dir) + HIST("hDEta_DPhi_PairEta_Qinv"), obs.deta, obs.dphi, pairEta, qinv);
      fRegistry.fill(HIST(dir) + HIST("hDEta_DPhi_PairPhi_Qinv"), obs.deta, obs.dphi, pairPhi, qinv);
      fRegistry.fill(HIST(dir) + HIST("hDPhi_Phi1_Phi2_Qinv"), obs.dphi, a.fPhi, b.fPhi, qinv);
      fRegistry.fill(HIST(dir) + HIST("hDEta_Eta1_Eta2_Qinv"), obs.deta, a.fEta, b.fEta, qinv);
    }
    fRegistry.fill(HIST(dir) + HIST("hDeltaRAngVsQinv"), qinv, std::min(deltaRAng, 0.499f));                       // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST(dir) + HIST("hDeltaRAng_PairEta_PairPhi"), std::min(deltaRAng, 0.499f), pairEta, pairPhi); // o2-linter: disable=magic-number (clamp to the axis)
    const int c1 = v0LegClass(a.legCounts), c2 = v0LegClass(b.legCounts);
    fRegistry.fill(HIST(dir) + HIST("hV0Class1_V0Class2_Qinv"), static_cast<float>(std::min(c1, c2)), static_cast<float>(std::max(c1, c2)), qinv);
    fRegistry.fill(HIST(dir) + HIST("hNLegsITSTPC_NLegsTPConly_Qinv"), static_cast<float>(a.legCounts.nITSTPC + b.legCounts.nITSTPC), static_cast<float>(a.legCounts.nTPCOnly + b.legCounts.nTPCOnly), qinv);
  }

  void fillLegSimilarity(PhotonCand const& a, PhotonCand const& b, float qinv)
  {
    if (!qaflags.cfgFillLegSimilarity.value || qinv >= qaflags.cfgMaxQinvForSparses.value) {
      return;
    }
    for (int il = 0; il < 2; ++il) { // o2-linter: disable=magic-number (e+ with e+, e- with e-)
      const auto& la = a.leg[static_cast<size_t>(il)];
      const auto& lb = b.leg[static_cast<size_t>(il)];
      const float ptSum = la.pt + lb.pt;
      const float ptRatio = (ptSum > KMinSigma) ? std::fabs(la.pt - lb.pt) / ptSum : 1.f;
      fRegistry.fill(HIST("Pair/LegSimilarity/hdEta_dPhi_ptRatio_dNsig_Qinv"), la.eta - lb.eta,
                     RecoDecay::constrainAngle(la.phi - lb.phi, -o2::constants::math::PI), ptRatio,
                     std::min(std::fabs(la.nsigEl - lb.nsigEl), 7.99f), qinv); // o2-linter: disable=magic-number (clamp to the axis)
    }
  }

  template <bool IsMC>
  void fillAnalyticV0QA(PhotonCand const& p)
  {
    if (!qaflags.cfgDoAnalyticV0QA.value) {
      return;
    }
    fRegistry.fill(HIST("Photon/AnalyticV0/hOk"), p.own.ok ? 1.f : 0.f);
    if (!p.own.ok) {
      return;
    }
    float truth = -1.f;
    if constexpr (IsMC) {
      if (p.mc.hasBothLabels()) {
        truth = p.mc.isTruePhoton ? 1.f : 0.f;
      }
    }
    const float rB = p.rConv();
    const float xA = p.own.vtx[0], yA = p.own.vtx[1], zA = p.own.vtx[2] + p.fVtxZ;
    fRegistry.fill(HIST("Photon/AnalyticV0/hdR_R_Truth"), std::hypot(xA, yA) - rB, rB, truth);
    fRegistry.fill(HIST("Photon/AnalyticV0/hdZ_R_Truth"), zA - p.fVz, rB, truth);
    fRegistry.fill(HIST("Photon/AnalyticV0/hdPhiConv_R_Truth"), RecoDecay::constrainAngle(std::atan2(yA, xA) - std::atan2(p.fVy, p.fVx), -o2::constants::math::PI), rB, truth);
    fRegistry.fill(HIST("Photon/AnalyticV0/hdPtRel_Pt_Truth"), (p.fPt > KMinSigma) ? (p.own.pt() - p.fPt) / p.fPt : 0.f, p.fPt, truth);
    fRegistry.fill(HIST("Photon/AnalyticV0/hdCosPA_Truth"), p.own.cospa - p.cospa, truth);
    fRegistry.fill(HIST("Photon/AnalyticV0/hPcaAna_PcaBuilder_Truth"), std::min(p.own.pca, 5.99f), std::min(p.pca, 5.99f), truth); // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST("Photon/AnalyticV0/hScore_R_Truth"), std::min(p.own.score, 29.9f), rB, truth);                             // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST("Photon/AnalyticV0/hPhotonLike_Truth"), isPhotonLike(p.own, mAltCuts) ? 1.f : 0.f, truth);
  }
  template <bool IsMC>
  void runDedupQA(std::vector<PhotonCand>& cands)
  {
    const auto n = static_cast<int>(cands.size());
    std::vector<int> nDupPartners(cands.size(), 0);
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        auto& a = cands[static_cast<size_t>(i)];
        auto& b = cands[static_cast<size_t>(j)];
        const float dEta = std::fabs(a.fEta - b.fEta);
        const float dPhi = std::fabs(RecoDecay::constrainAngle(a.fPhi - b.fPhi, -o2::constants::math::PI));
        const float ptSum = a.fPt + b.fPt;
        const float ptAsym = (ptSum > KMinSigma) ? std::fabs(a.fPt - b.fPt) / ptSum : 0.f;
        const float dR = std::fabs(a.rConv() - b.rConv());
        const bool flagged = pairutil::isDuplicatePair(a.dedupCand, b.dedupCand, mDedupCfg);
        fRegistry.fill(HIST("Dedup/hdEta_dPhi_ptAsym"), dEta, dPhi, ptAsym);
        fRegistry.fill(HIST("Dedup/hdR_ptAsym"), dR, ptAsym);
        for (int il = 0; il < 2; ++il) { // o2-linter: disable=magic-number (e+ and e- legs)
          const auto ld = pairutil::legDistance(a.dedupCand, b.dedupCand, il);
          float truth = 4.f; // o2-linter: disable=magic-number (no / invalid MC info)
          if constexpr (IsMC) {
            const int ida = (il == 0) ? a.mc.mcPosId : a.mc.mcNegId;
            const int idb = (il == 0) ? b.mc.mcPosId : b.mc.mcNegId;
            const int idbX = (il == 0) ? b.mc.mcNegId : b.mc.mcPosId;
            if (ida >= 0 && (idb >= 0 || idbX >= 0)) {
              const bool sameRecoTrack = (a.leg[static_cast<size_t>(il)].trackId == b.leg[static_cast<size_t>(il)].trackId) ||
                                         (a.leg[static_cast<size_t>(il)].trackId == b.leg[static_cast<size_t>(1 - il)].trackId);
              const bool sameMC = (ida == idb) || (ida == idbX);
              // 0 diff reco/diff MC, 1 diff reco/same MC, 2 same reco/same MC, 3 same reco/diff MC
              if (!sameRecoTrack) {
                truth = sameMC ? 1.f : 0.f;
              } else {
                truth = sameMC ? 2.f : 3.f; // o2-linter: disable=magic-number (class codes, see booking)
              }
            }
            const bool legFlagged = pairutil::legsIdentical(a.dedupCand, b.dedupCand, il, mDedupCfg);
            fRegistry.fill(HIST("Dedup/MC/hLegTrackTruthClass"), truth);
            fRegistry.fill(HIST("Dedup/MC/hLegTrackTruthClass_vs_Flagged"), truth, legFlagged ? 1.f : 0.f);
          }
          fRegistry.fill(HIST("Dedup/hLeg_dEta_dPhi_Truth"), ld.dEta, ld.dPhi, truth);
          fRegistry.fill(HIST("Dedup/hLeg_ptAsym_dEdxAsym_Truth"), ld.ptAsym, ld.dedxAsym, truth);
          fRegistry.fill(HIST("Dedup/hLeg_fShared_Truth"), ld.fShared, truth);
        }
        if constexpr (IsMC) {
          const auto dupClass = static_cast<float>(static_cast<int>(mcutil::classifyDupPair(a.mc, b.mc)));
          fRegistry.fill(HIST("Dedup/MC/hDupClass"), dupClass);
          fRegistry.fill(HIST("Dedup/MC/hDupClass_vs_Flagged"), dupClass, flagged ? 1.f : 0.f);
          fRegistry.fill(HIST("Dedup/MC/hDupClass_dEta_dPhi"), dupClass, dEta, dPhi);
          fRegistry.fill(HIST("Dedup/MC/hDupClass_ptAsym_dR"), dupClass, ptAsym, dR);
          const std::array<bool, 2> split = {
            a.leg[0].trackId != b.leg[0].trackId && a.mc.mcPosId >= 0 && a.mc.mcPosId == b.mc.mcPosId,
            a.leg[1].trackId != b.leg[1].trackId && a.mc.mcNegId >= 0 && a.mc.mcNegId == b.mc.mcNegId};
          if (split[0] || split[1]) {
            fRegistry.fill(HIST("Dedup/MC/hSplitLegPairClass"), dupClass);
            float v0Type = 3.f; // o2-linter: disable=magic-number (other)
            if (a.mc.isTruePhoton && b.mc.isTruePhoton && a.mc.posPhotonId == b.mc.posPhotonId) {
              v0Type = 0.f;
            } else if (a.mc.isTruePhoton != b.mc.isTruePhoton) {
              v0Type = 1.f;
            } else if (!a.mc.isTruePhoton && !b.mc.isTruePhoton) {
              v0Type = 2.f; // o2-linter: disable=magic-number (both fake)
            }
            fRegistry.fill(HIST("Dedup/MC/hSplitLegV0Type"), v0Type);
            for (int il = 0; il < 2; ++il) { // o2-linter: disable=magic-number (two legs)
              if (!split[static_cast<size_t>(il)]) {
                continue;
              }
              const auto ldSplit = pairutil::legDistance(a.dedupCand, b.dedupCand, il);
              fRegistry.fill(HIST("Dedup/MC/hSplitLeg_ptAsym_dEdxAsym_DupClass"), ldSplit.ptAsym, ldSplit.dedxAsym, dupClass);
              fRegistry.fill(HIST("Dedup/MC/hSplitLeg_fShared_DupClass"), ldSplit.fShared, dupClass);
            }
          }
          if (mcutil::classifyDupPair(a.mc, b.mc) == mcutil::DupClass::SameMcPhoton) {
            // 0 no split leg, 1 split e+, 2 split e-, 3 both
            const float splitType = static_cast<float>((split[0] ? 1 : 0) + (split[1] ? 2 : 0)); // o2-linter: disable=magic-number (bit code, see booking)
            fRegistry.fill(HIST("Dedup/MC/hSamePhotonDuplicateType"), splitType);
            fRegistry.fill(HIST("Dedup/MC/hSamePhotonDuplicateType_vs_Flagged"), splitType, flagged ? 1.f : 0.f);
            fRegistry.fill(HIST("Dedup/MC/hDupTrackClassMatrix"), static_cast<float>(a.dedupCand.nITSTPC), static_cast<float>(b.dedupCand.nITSTPC));
            fRegistry.fill(HIST("Dedup/MC/hDupDeltaRconv"), dR);
          }
        }
        if (flagged) {
          ++nDupPartners[static_cast<size_t>(i)];
          ++nDupPartners[static_cast<size_t>(j)];
          const bool keepA = pairutil::isBetterCand(a.dedupCand, b.dedupCand);
          (keepA ? b : a).dedupRejected = true;
          if constexpr (IsMC) {
            const auto& kept = keepA ? a : b;
            const auto& dropped = keepA ? b : a;
            float outcome = 4.f; // o2-linter: disable=magic-number (both true, distinct photons)
            if (kept.mc.isTruePhoton && !dropped.mc.isTruePhoton) {
              outcome = 0.f;
            } else if (!kept.mc.isTruePhoton && dropped.mc.isTruePhoton) {
              outcome = 1.f;
            } else if (kept.mc.isTruePhoton && dropped.mc.isTruePhoton && kept.mc.posPhotonId == dropped.mc.posPhotonId) {
              outcome = 2.f; // o2-linter: disable=magic-number (both true, same photon)
            } else if (!kept.mc.isTruePhoton && !dropped.mc.isTruePhoton) {
              outcome = 3.f; // o2-linter: disable=magic-number (both fake)
            }
            fRegistry.fill(HIST("Dedup/MC/hChoiceOutcome"), outcome);
          }
        }
      }
    }
    if constexpr (IsMC) {
      std::unordered_map<int, std::pair<int, float>> perMcPhoton; // photon id -> (count, reco pT of the last candidate)
      for (const auto& c : cands) {
        if (c.mc.isTruePhoton) {
          auto& entry = perMcPhoton[c.mc.posPhotonId];
          ++entry.first;
          entry.second = c.fPt;
        }
      }
      for (const auto& [mcId, entry] : perMcPhoton) {
        const auto nCand = static_cast<float>(std::min(entry.first, 6)); // o2-linter: disable=magic-number (last bin)
        fRegistry.fill(HIST("Dedup/MC/hNCandPerMcPhoton"), nCand);
        fRegistry.fill(HIST("Dedup/MC/hNCandPerMcPhoton_vs_Pt"), nCand, entry.second);
      }
    }
    int nRejected = 0;
    for (int i = 0; i < n; ++i) {
      fRegistry.fill(HIST("Dedup/hNDupPartners_vs_Pt"), static_cast<float>(std::min(nDupPartners[static_cast<size_t>(i)], 5)), cands[static_cast<size_t>(i)].fPt); // o2-linter: disable=magic-number (last bin)
      nRejected += cands[static_cast<size_t>(i)].dedupRejected ? 1 : 0;
    }
    fRegistry.fill(HIST("Dedup/hNCandBeforeAfter"), 0.f, static_cast<float>(n));
    fRegistry.fill(HIST("Dedup/hNCandBeforeAfter"), 1.f, static_cast<float>(n - nRejected));
  }

  template <bool IsMC>
  void fillCrossPairQA(CrossObs const& c, PairObs const& obs, float truthAxis, PhotonCand const& a, PhotonCand const& b)
  {
    const float qinv = obs.q.qinv;
    fRegistry.fill(HIST("Pair/CrossPair/hQinv_NAltValid_MaxAltDca_Fit"), qinv, static_cast<float>(c.nAltValidFit), std::min(c.maxDcaFit, 11.9f));     // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST("Pair/CrossPair/hQinv_NAltValid_MaxAltDca_Ana"), qinv, static_cast<float>(c.nAltValidAna), std::min(c.maxDcaAna, 11.9f));     // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST("Pair/CrossPair/hQinv_NAltValid_MaxAltDcaXY_Fit"), qinv, static_cast<float>(c.nAltValidFit), std::min(c.maxDcaXYFit, 11.9f)); // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST("Pair/CrossPair/hQinv_NAltValid_MaxAltDcaXY_Ana"), qinv, static_cast<float>(c.nAltValidAna), std::min(c.maxDcaXYAna, 11.9f)); // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST("Pair/CrossPair/hQinv_NAltValid_MaxAltDcaZ_Fit"), qinv, static_cast<float>(c.nAltValidFit), std::min(c.maxDcaZFit, 11.9f));   // o2-linter: disable=magic-number (clamp to the axis)
    fRegistry.fill(HIST("Pair/CrossPair/hQinv_NAltValid_MaxAltDcaZ_Ana"), qinv, static_cast<float>(c.nAltValidAna), std::min(c.maxDcaZAna, 11.9f));   // o2-linter: disable=magic-number (clamp to the axis)
    if (c.evaluated) {
      for (int k = 0; k < 2; ++k) { // o2-linter: disable=magic-number (the two alternative pairings)
        if (c.fit[static_cast<size_t>(k)].fitted && c.ana[static_cast<size_t>(k)].ok) {
          fRegistry.fill(HIST("Pair/CrossPair/hAltDca_Fit_vs_Ana"), std::min(c.fit[static_cast<size_t>(k)].dca, 11.9f), std::min(c.ana[static_cast<size_t>(k)].pca, 11.9f)); // o2-linter: disable=magic-number (clamp to the axis)
          fRegistry.fill(HIST("Pair/CrossPair/hAltCosPA_Fit_vs_Ana"), c.fit[static_cast<size_t>(k)].cospa, c.ana[static_cast<size_t>(k)].cospa);
        }
      }
      fRegistry.fill(HIST("Pair/CrossPair/hNAltValid_Fit_vs_Ana"), static_cast<float>(c.nAltValidFit), static_cast<float>(c.nAltValidAna));
    }
    if (c.meeOverQ < 900.f) { // o2-linter: disable=magic-number (qinv was zero, skip)
      fRegistry.fill(HIST("Pair/CrossPair/hMeeRatio_dR_Qinv"), c.meeOverQ, obs.deltaR, qinv);
    }
    fRegistry.fill(HIST("Pair/CrossPair/hRepairDeltaMee_OwnMee_Qinv"), c.deltaMee, c.ownMeeMax, qinv);
    if (qinv < qaflags.cfgMaxQinvForSparses.value) {
      fRegistry.fill(HIST("Pair/CrossPair/hMee1_Mee2_M4_Qinv"), c.ownMee[0], c.ownMee[1], c.m4, qinv);
    }
    if constexpr (IsMC) {
      if (truthAxis < 0.f) {
        return;
      }
      fRegistry.fill(HIST("Pair/MC/CrossPair/hQinv_NAltValid_MaxAltDca_Fit_Type"), qinv, static_cast<float>(c.nAltValidFit), std::min(c.maxDcaFit, 11.9f), truthAxis);   // o2-linter: disable=magic-number (clamp to the axis)
      fRegistry.fill(HIST("Pair/MC/CrossPair/hQinv_NAltValid_MaxAltDca_Ana_Type"), qinv, static_cast<float>(c.nAltValidAna), std::min(c.maxDcaAna, 11.9f), truthAxis);   // o2-linter: disable=magic-number (clamp to the axis)
      fRegistry.fill(HIST("Pair/MC/CrossPair/hQinv_NAltValid_MaxAltDcaZ_Fit_Type"), qinv, static_cast<float>(c.nAltValidFit), std::min(c.maxDcaZFit, 11.9f), truthAxis); // o2-linter: disable=magic-number (clamp to the axis)
      fRegistry.fill(HIST("Pair/MC/CrossPair/hQinv_NAltValid_MaxAltDcaZ_Ana_Type"), qinv, static_cast<float>(c.nAltValidAna), std::min(c.maxDcaZAna, 11.9f), truthAxis); // o2-linter: disable=magic-number (clamp to the axis)
      if (c.meeOverQ < 900.f) {                                                                                                                                          // o2-linter: disable=magic-number (qinv was zero, skip)
        fRegistry.fill(HIST("Pair/MC/CrossPair/hMeeRatio_Qinv_Type"), c.meeOverQ, qinv, truthAxis);
      }
      if (qinv < qaflags.cfgMaxQinvForSparses.value) {
        fRegistry.fill(HIST("Pair/MC/CrossPair/hMee1_Mee2_M4_Qinv_Type"), c.ownMee[0], c.ownMee[1], c.m4, qinv, truthAxis);
      }
      fRegistry.fill(HIST("Pair/MC/CrossPair/hQinv_MaxDcaZ_Type"), qinv, std::min(std::max(std::fabs(a.fDcaZToPV), std::fabs(b.fDcaZToPV)), 9.9f), truthAxis);    // o2-linter: disable=magic-number (clamp to the axis)
      fRegistry.fill(HIST("Pair/MC/CrossPair/hQinv_MaxDcaXY_Type"), qinv, std::min(std::max(std::fabs(a.fDcaXYToPV), std::fabs(b.fDcaXYToPV)), 2.9f), truthAxis); // o2-linter: disable=magic-number (clamp to the axis)
      const Emulation emu = emulate(c, a, b);
      fRegistry.fill(HIST("Pair/MC/CrossPair/hEmulation_vs_Truth"), qinv, truthAxis, static_cast<float>(emu.outcome));
      if (emu.outcome >= 1) {
        fRegistry.fill(HIST("Pair/MC/CrossPair/hEmuDeltaScore_vs_Truth"), qinv, truthAxis, std::clamp(emu.deltaScore, -19.99f, 19.99f)); // o2-linter: disable=magic-number (clamp to the axis)
      }
      fRegistry.fill(HIST("Pair/MC/CrossPair/hFakeSubtype_Emu_Qinv"), static_cast<float>(mcutil::classifyFakeSubtype(a.mc, b.mc)), static_cast<float>(emu.outcome), qinv);
    }
  }

  template <mcutil::PairTruthType T>
  void fillTypeObservables(PairObs const& obs)
  {
    constexpr auto dir =
      (T == mcutil::PairTruthType::TrueTrueDistinct)     ? "Pair/MC/TrueTrueDistinct/"
      : (T == mcutil::PairTruthType::TrueTrueSamePhoton) ? "Pair/MC/TrueTrueSamePhoton/"
      : (T == mcutil::PairTruthType::SharedMcLeg)        ? "Pair/MC/SharedMcLeg/"
      : (T == mcutil::PairTruthType::TrueFake)           ? "Pair/MC/TrueFake/"
      : (T == mcutil::PairTruthType::FakeFake)           ? "Pair/MC/FakeFake/"
                                                         : "Pair/MC/Pi0Daughters/";
    const float qinv = obs.q.qinv;
    fRegistry.fill(HIST(dir) + HIST("hDeltaRVsQinv"), qinv, obs.deltaR);
    fRegistry.fill(HIST(dir) + HIST("hDeltaZVsQinv"), qinv, obs.deltaZ);
    fRegistry.fill(HIST(dir) + HIST("hDeltaR3DVsQinv"), qinv, obs.deltaR3D);
    fRegistry.fill(HIST(dir) + HIST("hDEtaDPhi"), obs.deta, obs.dphi);
    if (qaflags.cfgFillTypeSparses.value) {
      fRegistry.fill(HIST(dir) + HIST("hDEtaDPhi_Qinv"), obs.deta, obs.dphi, qinv);
      fRegistry.fill(HIST(dir) + HIST("hDeltaR_DeltaZ_Qinv"), obs.deltaR, obs.deltaZ, qinv);
    }
  }

  void fillTypeObservables(mcutil::PairTruthType t, PairObs const& obs)
  {
    switch (t) {
      case mcutil::PairTruthType::TrueTrueDistinct:
        fillTypeObservables<mcutil::PairTruthType::TrueTrueDistinct>(obs);
        break;
      case mcutil::PairTruthType::TrueTrueSamePhoton:
        fillTypeObservables<mcutil::PairTruthType::TrueTrueSamePhoton>(obs);
        break;
      case mcutil::PairTruthType::SharedMcLeg:
        fillTypeObservables<mcutil::PairTruthType::SharedMcLeg>(obs);
        break;
      case mcutil::PairTruthType::TrueFake:
        fillTypeObservables<mcutil::PairTruthType::TrueFake>(obs);
        break;
      case mcutil::PairTruthType::FakeFake:
        fillTypeObservables<mcutil::PairTruthType::FakeFake>(obs);
        break;
      case mcutil::PairTruthType::Pi0Daughters:
        fillTypeObservables<mcutil::PairTruthType::Pi0Daughters>(obs);
        break;
      default:
        break;
    }
  }

  template <typename TMCParticles>
  void fillAncestry(TMCParticles const& mcParticles, PhotonCand const& a, PhotonCand const& b, PairObs const& obs, float truthAxis)
  {
    if (!qaflags.cfgDoAncestry.value) {
      return;
    }
    const float qinv = obs.q.qinv;
    if (qinv >= qaflags.cfgMaxQinvForSparses.value) {
      return;
    }
    const int maxGen = qaflags.cfgAncestorMaxGen.value;
    const auto census = mcutil::censusMothers(a.mc, b.mc);
    const auto splitAxis = static_cast<float>(census.splitCode());
    const auto censusCode = static_cast<float>(census.nConvLegs + 5 * census.nMothers); // o2-linter: disable=magic-number (code = legs + 5 mothers)
    fRegistry.fill(HIST("Pair/MC/Ancestry/hCensus_Qinv_Type"), censusCode, qinv, truthAxis);
    fRegistry.fill(HIST("Pair/MC/Ancestry/hSplit_Qinv_Type"), splitAxis, qinv, truthAxis);

    if (census.nMothers == 2) { // o2-linter: disable=magic-number (two parent photons)
      const int idA = census.distinct[0], idB = census.distinct[1];
      const float qTrue = mcutil::trueQinvOfPhotons(mcParticles, idA, idB);
      const auto ancKind = static_cast<float>(static_cast<int>(mcutil::commonAncestorKind(mcParticles, idA, idB, maxGen)));
      fRegistry.fill(HIST("Pair/MC/Ancestry/hQtrue_Qreco_Split"), qTrue, qinv, splitAxis);
      fRegistry.fill(HIST("Pair/MC/Ancestry/hAncKind_Qinv_Type"), ancKind, qinv, truthAxis);
      const int origA = mcutil::originParticleId(mcParticles, idA), origB = mcutil::originParticleId(mcParticles, idB);
      const auto kindA = static_cast<float>(static_cast<int>(mcutil::ancestorKindOfPdg(mcParticles.iteratorAt(origA).pdgCode())));
      const auto kindB = static_cast<float>(static_cast<int>(mcutil::ancestorKindOfPdg(mcParticles.iteratorAt(origB).pdgCode())));
      const float same = (origA == origB) ? 1.f : 0.f;
      fRegistry.fill(HIST("Pair/MC/Ancestry/hAnc1_Anc2_Same"), std::min(kindA, kindB), std::max(kindA, kindB), same);
      fRegistry.fill(HIST("Pair/MC/Ancestry/hSame_Qinv_Type"), same, qinv, truthAxis);
      const int srcA = photonSource(mcParticles, idA), srcB = photonSource(mcParticles, idB);
      fRegistry.fill(HIST("Pair/MC/Ancestry/hSrc1_Src2_Type"), static_cast<float>(std::min(srcA, srcB)), static_cast<float>(std::max(srcA, srcB)), truthAxis);
      std::array<float, 3> cA{}, cB{};
      if (mcutil::conversionPoint(mcParticles, idA, cA) && mcutil::conversionPoint(mcParticles, idB, cB)) {
        const float rA = std::hypot(cA[0], cA[1]), rB = std::hypot(cB[0], cB[1]);
        const float d3 = std::hypot(cA[0] - cB[0], cA[1] - cB[1], cA[2] - cB[2]);
        fRegistry.fill(HIST("Pair/MC/Ancestry/hTrueConv_d3D_Qinv_Type"), std::min(d3, 49.9f), qinv, truthAxis);                                                   // o2-linter: disable=magic-number (clamp to the axis)
        fRegistry.fill(HIST("Pair/MC/Ancestry/hTrueConv_dR_dZ_Type"), std::min(std::fabs(rA - rB), 24.9f), std::min(std::fabs(cA[2] - cB[2]), 24.9f), truthAxis); // o2-linter: disable=magic-number (clamp to the axes)
        fRegistry.fill(HIST("Pair/MC/Ancestry/hRMeanTrue_DRTrue_Split"), 0.5f * (rA + rB), std::min(std::fabs(rA - rB), 24.9f), splitAxis);                       // o2-linter: disable=magic-number (clamp to the axis)
      }
      const float theta = mcutil::trueOpeningAngle(mcParticles, idA, idB);
      if (theta >= 0.f) {
        fRegistry.fill(HIST("Pair/MC/Ancestry/hTrueConv_dTheta_Qtrue_Type"), std::min(theta * 1000.f, 49.9f), qTrue, truthAxis); // o2-linter: disable=magic-number (rad -> mrad, clamp)
      }
      fRegistry.fill(HIST("Pair/MC/Ancestry/hTrueConv_DeltaQ_Qtrue_Type"), qinv - qTrue, qTrue, truthAxis);
    }

    // per V0 candidate
    for (const PhotonCand& p : {std::cref(a), std::cref(b)}) {
      const auto reco = static_cast<float>(static_cast<int>(mcutil::classifyV0Reco(p.mc)));
      const float rReco = p.rConv();
      {
        const auto srcPos = static_cast<float>(photonSource(mcParticles, p.mc.posPhotonId));
        const auto srcNeg = static_cast<float>(photonSource(mcParticles, p.mc.negPhotonId));
        fRegistry.fill(HIST("Pair/MC/Ancestry/hSrcPos_SrcNeg_Reco"), srcPos, srcNeg, reco);
        fRegistry.fill(HIST("Pair/MC/Ancestry/hSrcPos_SrcNeg_Type"), srcPos, srcNeg, truthAxis);
      }
      fRegistry.fill(HIST("Pair/MC/Ancestry/hEta_Phi_Reco"), p.fEta, p.fPhi, reco);
      fRegistry.fill(HIST("Pair/MC/Ancestry/hR_Reco_Qinv"), std::min(rReco, 99.9f), reco, qinv); // o2-linter: disable=magic-number (clamp to the axis)
      fRegistry.fill(HIST("Pair/MC/Ancestry/hVx_Vy_Split"), p.fVx, p.fVy, splitAxis);
      fRegistry.fill(HIST("Pair/MC/Ancestry/hPsiPair_Chi2_Split"), std::fabs(p.psipair), std::min(p.chi2, 49.9f), splitAxis); // o2-linter: disable=magic-number (clamp to the axis)
      fRegistry.fill(HIST("Pair/MC/Ancestry/hPsiPair_R_Split"), std::fabs(p.psipair), std::min(rReco, 99.9f), splitAxis);     // o2-linter: disable=magic-number (clamp to the axis)
      if (p.mc.mcPosId >= 0 && p.mc.mcNegId >= 0) {
        const auto ePos = mcParticles.iteratorAt(p.mc.mcPosId);
        const auto eNeg = mcParticles.iteratorAt(p.mc.mcNegId);
        const float rPos = std::hypot(static_cast<float>(ePos.vx()), static_cast<float>(ePos.vy()));
        const float rNeg = std::hypot(static_cast<float>(eNeg.vx()), static_cast<float>(eNeg.vy()));
        fRegistry.fill(HIST("Pair/MC/Ancestry/hVxTrue_VyTrue_Split"), static_cast<float>(ePos.vx()), static_cast<float>(ePos.vy()), splitAxis);
        fRegistry.fill(HIST("Pair/MC/Ancestry/hRTruePos_RTrueNeg_Reco"), std::min(rPos, 99.9f), std::min(rNeg, 99.9f), reco);                // o2-linter: disable=magic-number (clamp to the axes)
        fRegistry.fill(HIST("Pair/MC/Ancestry/hRTrueMean_RReco_Reco"), std::min(0.5f * (rPos + rNeg), 99.9f), std::min(rReco, 99.9f), reco); // o2-linter: disable=magic-number (clamp to the axes)
      }
      if (p.mc.posPhotonId >= 0 && p.mc.negPhotonId >= 0) {
        const int oPos = mcutil::originParticleId(mcParticles, p.mc.posPhotonId);
        const int oNeg = mcutil::originParticleId(mcParticles, p.mc.negPhotonId);
        const auto kPos = static_cast<float>(static_cast<int>(mcutil::ancestorKindOfPdg(mcParticles.iteratorAt(oPos).pdgCode())));
        const auto kNeg = static_cast<float>(static_cast<int>(mcutil::ancestorKindOfPdg(mcParticles.iteratorAt(oNeg).pdgCode())));
        fRegistry.fill(HIST("Pair/MC/Ancestry/hAncPos_AncNeg_Reco"), kPos, kNeg, reco);
        fRegistry.fill(HIST("Pair/MC/Ancestry/hAncSame_Reco"), (oPos == oNeg) ? 1.f : 0.f, reco);
      }
    }
  }

  template <bool IsMC, typename TPhotons, typename TLegs, typename TMCParticles>
  void collectPhotons(TPhotons const& photonsColl, TMCParticles const& mcParticles, float vtxZ, std::vector<PhotonCand>& cands)
  {
    cands.clear();
    cands.reserve(photonsColl.size());
    for (const auto& g : photonsColl) {
      if (!fV0PhotonCut.template IsSelected<decltype(g), TLegs>(g)) {
        continue;
      }
      const auto pos = g.template posTrack_as<TLegs>();
      const auto ele = g.template negTrack_as<TLegs>();
      PhotonCand c;
      c.gi = g.globalIndex();
      c.fPt = g.pt();
      c.fEta = g.eta();
      c.fPhi = g.phi();
      c.fVx = g.vx();
      c.fVy = g.vy();
      c.fVz = g.vz();
      c.fVtxZ = vtxZ;
      c.fDcaXYToPV = static_cast<float>(g.dcaXYtopv());
      c.fDcaZToPV = static_cast<float>(g.dcaZtopv());
      c.psipair = static_cast<float>(g.psipair());
      c.phiv = static_cast<float>(g.phiv());
      c.chi2 = static_cast<float>(g.chiSquareNDF());
      c.cospa = static_cast<float>(g.cospa());
      c.pca = static_cast<float>(g.pca());
      auto fillLeg = [&](LegLite& l, auto const& t) {
        l.pt = static_cast<float>(t.pt());
        l.eta = static_cast<float>(t.eta());
        l.phi = static_cast<float>(t.phi());
        l.nsigEl = static_cast<float>(t.tpcNSigmaEl());
        l.dcaZ = static_cast<float>(t.dcaZ());
        l.nClsFindable = static_cast<float>(t.tpcNClsFindable());
        l.nClsITS = static_cast<float>(t.itsNCls());
        l.timed = t.hasITS() || t.hasTRD() || t.hasTOF();
        l.foreign = (t.collisionId() != g.collisionId());
        l.trackId = t.trackId();
      };
      fillLeg(c.leg[0], pos);
      fillLeg(c.leg[1], ele);
      c.legCounts = pairutil::getV0PhotonLegCounts(pos, ele);
      c.dedupCand = pairutil::makeDedupCand(g, pos, ele);
      if (qaflags.cfgDoAnalyticV0QA.value || qaflags.cfgDoCrossPairQA.value) {
        c.own = buildAnalyticV0(c.legParam(0), c.legParam(1), 10.f * mBzT, qaflags.cfgScoreWeight.value); // o2-linter: disable=magic-number (Tesla -> kGauss)
      }
      if constexpr (IsMC) {
        c.mc = mcutil::makePhotonMCInfo<TLegs>(g, mcParticles);
      }
      cands.push_back(c);
    }
  }

  template <bool IsMC, typename TCollision, typename TMCParticles>
  void runPairQa(TCollision const& collision, std::vector<PhotonCand> const& cands, TMCParticles const& mcParticles)
  {
    const int evOwn = [&]() {
      if constexpr (IsMC) {
        return collision.has_emmcevent() ? static_cast<int>(collision.emmceventId()) : -1;
      } else {
        return -1;
      }
    }();
    for (size_t i = 0; i < cands.size(); ++i) {
      for (size_t j = i + 1; j < cands.size(); ++j) {
        const auto& g1 = cands[i];
        const auto& g2 = cands[j];
        if (g1.sharesTrackWith(g2)) {
          continue;
        }
        const bool dedupWouldRemove = g1.dedupRejected || g2.dedupRejected;
        if (dedup.cfgDoDedup.value && dedupWouldRemove) {
          continue;
        }
        if (!passAsymmetryCut(g1.fPt, g2.fPt)) {
          continue;
        }
        if (!pairutil::pairBelowQmax(g1.fPt, g1.fEta, g1.fPhi, g2.fPt, g2.fEta, g2.fPhi, qaflags.cfgMaxQinvForQA.value)) {
          continue;
        }
        const PairObs obs = buildPairObs(g1, g2);
        if (!obs.valid) {
          continue;
        }
        const float qinv = obs.q.qinv;

        float truthAxis = -1.f;
        mcutil::PairTruthType truthType = mcutil::PairTruthType::Unknown;
        if constexpr (IsMC) {
          if (g1.mc.hasBothLabels() && g2.mc.hasBothLabels()) {
            truthType = mcutil::pairTruthType(g1.mc, g2.mc, mcParticles);
            truthAxis = static_cast<float>(static_cast<int>(truthType));
          }
        }

        fillPairQA<0>(obs, g1, g2);
        fillLegSimilarity(g1, g2, qinv);

        if (!passPairCuts(obs, g1, g2)) {
          continue;
        }
        CrossObs cross;
        if (qaflags.cfgDoCrossPairQA.value || crosspair.cfgDoCrossPairCut.value) {
          cross = computeCrossObs(g1, g2, qinv);
          if (qaflags.cfgDoCrossPairQA.value) {
            fillCrossPairQA<IsMC>(cross, obs, truthAxis, g1, g2);
          }
          const bool vetoed = !passCrossPairVeto(cross);
          fRegistry.fill(HIST("Pair/CrossPair/hVetoCounter"), vetoed ? 1.f : 0.f);
          if (crosspair.cfgDoCrossPairCut.value && vetoed) {
            continue;
          }
        }

        fillPairQA<1>(obs, g1, g2);

        if constexpr (IsMC) {
          if (truthAxis < 0.f) {
            continue;
          }
          fRegistry.fill(HIST("Pair/MC/hTruthTypeVsQinv"), qinv, truthAxis);
          fRegistry.fill(HIST("Pair/MC/hTruthTypeVsKt"), obs.q.kt, truthAxis);
          fRegistry.fill(HIST("Pair/MC/hTruthType_Qinv_DedupFlag"), truthAxis, qinv, dedupWouldRemove ? 1.f : 0.f);
          if (truthType == mcutil::PairTruthType::FakeFake || truthType == mcutil::PairTruthType::TrueFake) {
            fRegistry.fill(HIST("Pair/MC/hFakeSubtypeVsQinv"), qinv, static_cast<float>(mcutil::classifyFakeSubtype(g1.mc, g2.mc)));
          }
          fillTypeObservables(truthType, obs);
          {
            const int c1 = v0LegClass(g1.legCounts), c2 = v0LegClass(g2.legCounts);
            fRegistry.fill(HIST("Pair/MC/hV0Class1_V0Class2_Qinv_Type"), static_cast<float>(std::min(c1, c2)), static_cast<float>(std::max(c1, c2)), qinv, truthAxis);
            fRegistry.fill(HIST("Pair/MC/hV0Class_IsTrue_Qinv"), static_cast<float>(c1), g1.mc.isTruePhoton ? 1.f : 0.f, qinv);
            fRegistry.fill(HIST("Pair/MC/hV0Class_IsTrue_Qinv"), static_cast<float>(c2), g2.mc.isTruePhoton ? 1.f : 0.f, qinv);
          }
          for (const PhotonCand& p : {std::cref(g1), std::cref(g2)}) {
            if (qinv < qaflags.cfgMaxQinvForSparses.value) {
              fRegistry.fill(HIST("Pair/MC/hPsiPair_PhiV_Qinv_Type_IsTrue"), std::fabs(p.psipair), p.phiv, qinv, truthAxis, p.mc.isTruePhoton ? 1.f : 0.f);
            }
            const float rReco = p.rConv();
            const std::array<int, 2> mcLegIds{p.mc.mcPosId, p.mc.mcNegId};
            for (const auto& mcLegId : mcLegIds) {
              const auto legPart = mcParticles.iteratorAt(mcLegId);
              const float rMc = std::hypot(static_cast<float>(legPart.vx()), static_cast<float>(legPart.vy()));
              fRegistry.fill(HIST("Pair/MC/hRreco_minus_Rmc"), rReco - rMc, rMc, p.mc.isTruePhoton ? 1.f : 0.f);
            }
          }
          {
            const std::array<int, 4> legIds = {g1.mc.mcPosId, g1.mc.mcNegId, g2.mc.mcPosId, g2.mc.mcNegId};
            std::array<int, 4> ev{};
            int nOwn = 0;
            for (size_t il = 0; il < legIds.size(); ++il) {
              ev[il] = static_cast<int>(mcParticles.iteratorAt(legIds[il]).emmceventId());
              nOwn += (ev[il] == evOwn) ? 1 : 0;
            }
            const float pattern = (ev[0] != ev[1] || ev[2] != ev[3]) ? 2.f : (ev[0] != ev[2]) ? 1.f
                                                                                              : 0.f; // o2-linter: disable=magic-number (pattern codes, see booking)
            fRegistry.fill(HIST("Pair/MC/hMcEventPattern_NOwn_Qinv_Type"), pattern, static_cast<float>(nOwn), qinv, truthAxis);
          }
          if (qinv < qaflags.cfgMaxQinvForMCQA.value) {
            fillAncestry(mcParticles, g1, g2, obs, truthAxis);
          }
        }
      }
    }
  }

  template <bool IsMC, typename TCollisions, typename TPhotons, typename TLegs, typename TMCParticles>
  void runCollisions(TCollisions const& collisions, TPhotons const& photons, TLegs const& /*legs*/, TMCParticles const& mcParticles)
  {
    std::vector<PhotonCand> cands;
    for (const auto& collision : collisions) {
      initCCDB(collision);
      const std::array<float, 3> cent = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (cent[centralitySelection.cfgCentEstimator] < centralitySelection.cfgCentMin ||
          centralitySelection.cfgCentMax < cent[centralitySelection.cfgCentEstimator]) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // o2-linter: disable=magic-number (counter bin, as in photonhbt)
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // o2-linter: disable=magic-number (counter bin, as in photonhbt)
      auto photonsColl = photons.sliceBy(perCollisionPCM, collision.globalIndex());

      collectPhotons<IsMC, decltype(photonsColl), TLegs>(photonsColl, mcParticles, collision.posZ(), cands);
      fRegistry.fill(HIST("Photon/hNPhotonsPerEvent"), static_cast<float>(std::min(cands.size(), static_cast<size_t>(20)))); // o2-linter: disable=magic-number (last bin)
      for (const auto& c : cands) {
        fRegistry.fill(HIST("Photon/hPtEtaPhi"), c.fPt, c.fEta, c.fPhi);
        const auto cls = static_cast<float>(v0LegClass(c.legCounts));
        fRegistry.fill(HIST("Photon/hV0LegClass"), cls);
        fRegistry.fill(HIST("Photon/hV0LegClass_vs_Pt"), cls, c.fPt);
        fRegistry.fill(HIST("Photon/hV0LegClass_vs_R"), cls, std::min(c.rConv(), 99.9f)); // o2-linter: disable=magic-number (clamp to the axis)
        fillAnalyticV0QA<IsMC>(c);
      }
      runDedupQA<IsMC>(cands);
      if (cands.size() < 2) { // o2-linter: disable=magic-number (suplicate number)
        continue;
      }
      runPairQa<IsMC>(collision, cands, mcParticles);
    }
  }

  void processReco(FilteredMyCollisions const& collisions, MyV0Photons const& v0photons, aod::V0Legs const& v0legs)
  {
    runCollisions<false>(collisions, v0photons, v0legs, v0legs /* placeholder, unused in data */);
  }
  PROCESS_SWITCH(PairQCTask, processReco, "pair QC on reconstructed data", true);

  void processMC(FilteredMyMCCollisions const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs,
                 aod::EMMCParticles const& mcParticles, aod::EMMCEvents const& /*mcEvents*/)
  {
    runCollisions<true>(collisions, v0photons, v0legs, mcParticles);
  }
  PROCESS_SWITCH(PairQCTask, processMC, "pair QC on MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PairQCTask>(cfgc)};
}
