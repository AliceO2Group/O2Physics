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
/// \brief V0 photon HBT analysis on the artificial events of builder run 2 (correlation functions, pair-class QA)
///        and the generated-level MC. Helpers in PWGEM/PhotonMeson/Utils/PhotonHbtUtilities.h.
/// \author Daiki Sekihata, daiki.sekihata@cern.ch
///         Stefanie Mrozinski, stefanie.mrozinski@cern.ch

// =====================================================================================================================
// INCLUDES AND ALIASES
// =====================================================================================================================

#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PhotonHbtUtilities.h"
#include "PWGEM/PhotonMeson/Utils/V0DedupUtilities.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsTPC/VDriftCorrFact.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h> // IWYU pragma: keep
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TPDGCode.h>
#include <TString.h>

#include <GPUROOTCartesianFwd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <initializer_list>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils;
namespace pairutil = o2::aod::pwgem::photonmeson::utils::pairutil;
using namespace o2::pwgem::photonmeson;
using namespace o2::aod::pwgem::photonmeson::hbt;

// generated-level MC (processGenerated): collisions, photons and legs of the classic derived tables
using MyCollisionsMC = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001, aod::EMMCEventLabels>;
using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsPhiVPsi>;
using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

// output of builder run 2 (processPairs): photons of the artificial events of PCMMixCollisions
using MyV0PhotonsMix = soa::Join<aod::V0PhotonsKFMix, aod::V0PhotonsPhiVPsiMix>;

struct Photonhbt {

  // ===================================================================================================================
  // CONFIGURABLES
  // - event: event cuts of the generated-level part, centrality window
  // - photon: photon cut (fV0PhotonCut), photon-class pair selection
  // - pair cuts: pair geometry, same-sign leg separation (pair merging), crossed (leg-swapped) combinations
  // - mixing: pools, co-reconstructed mixed event, partner pairing
  // - output and QA: blocks on/off, correlation-function axes, QA switches, generated-level MC
  // ===================================================================================================================

  // ─── Event ────────────────────────────────────────────────────────────────
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
    std::string prefix = "centralitySelection_group";
    Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  } centralitySelection;

  // ─── Photon ───────────────────────────────────────────────────────────────
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
    Configurable<float> cfgMinNCrossedRowsOverFindable{"cfgMinNCrossedRowsOverFindable", 0.8f, "min. TPC crossed rows / findable clusters of the legs"};
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

  Configurable<bool> cfgDoPhotonClassPairCut{"cfgDoPhotonClassPairCut", false, "apply photon-class pair selection A x B (leg track composition)"};

  struct : ConfigurableGroup {
    std::string prefix = "photonclassA_group";
    Configurable<int> cfgMinNLegsITSTPC{"cfgMinNLegsITSTPC", 0, "min number of ITS-TPC legs (0-2)"};
    Configurable<int> cfgMaxNLegsITSTPC{"cfgMaxNLegsITSTPC", 2, "max number of ITS-TPC legs (0-2)"};
    Configurable<int> cfgMinNLegsITSOnly{"cfgMinNLegsITSOnly", 0, "min number of ITS-only legs (0-2)"};
    Configurable<int> cfgMaxNLegsITSOnly{"cfgMaxNLegsITSOnly", 2, "max number of ITS-only legs (0-2)"};
    Configurable<int> cfgMinNLegsTPCOnly{"cfgMinNLegsTPCOnly", 0, "min number of TPC-only legs (0-2)"};
    Configurable<int> cfgMaxNLegsTPCOnly{"cfgMaxNLegsTPCOnly", 2, "max number of TPC-only legs (0-2)"};
    Configurable<int> cfgMinNLegsTRD{"cfgMinNLegsTRD", 0, "min number of legs with TRD (0-2)"};
    Configurable<int> cfgMaxNLegsTRD{"cfgMaxNLegsTRD", 2, "max number of legs with TRD (0-2)"};
    Configurable<int> cfgMinNLegsTOF{"cfgMinNLegsTOF", 0, "min number of legs with TOF (0-2)"};
    Configurable<int> cfgMaxNLegsTOF{"cfgMaxNLegsTOF", 2, "max number of legs with TOF (0-2)"};
  } photonclassA;

  struct : ConfigurableGroup {
    std::string prefix = "photonclassB_group";
    Configurable<int> cfgMinNLegsITSTPC{"cfgMinNLegsITSTPC", 0, "min number of ITS-TPC legs (0-2)"};
    Configurable<int> cfgMaxNLegsITSTPC{"cfgMaxNLegsITSTPC", 2, "max number of ITS-TPC legs (0-2)"};
    Configurable<int> cfgMinNLegsITSOnly{"cfgMinNLegsITSOnly", 0, "min number of ITS-only legs (0-2)"};
    Configurable<int> cfgMaxNLegsITSOnly{"cfgMaxNLegsITSOnly", 2, "max number of ITS-only legs (0-2)"};
    Configurable<int> cfgMinNLegsTPCOnly{"cfgMinNLegsTPCOnly", 0, "min number of TPC-only legs (0-2)"};
    Configurable<int> cfgMaxNLegsTPCOnly{"cfgMaxNLegsTPCOnly", 2, "max number of TPC-only legs (0-2)"};
    Configurable<int> cfgMinNLegsTRD{"cfgMinNLegsTRD", 0, "min number of legs with TRD (0-2)"};
    Configurable<int> cfgMaxNLegsTRD{"cfgMaxNLegsTRD", 2, "max number of legs with TRD (0-2)"};
    Configurable<int> cfgMinNLegsTOF{"cfgMinNLegsTOF", 0, "min number of legs with TOF (0-2)"};
    Configurable<int> cfgMaxNLegsTOF{"cfgMaxNLegsTOF", 2, "max number of legs with TOF (0-2)"};
  } photonclassB;

  // ─── Pair cuts ────────────────────────────────────────────────────────────
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
  } ggpaircuts;

  struct : ConfigurableGroup {
    std::string prefix = "pairsep_group";
    Configurable<bool> cfgDoPairSepQA{"cfgDoPairSepQA", true, "fill SE/ME same-sign leg separation histograms"};
    Configurable<float> cfgPairSepMaxQinv{"cfgPairSepMaxQinv", 0.5f, "fill PairSep histograms only for qinv < this value (<=0: no gate)"};
    Configurable<float> cfgSigRPhi{"cfgSigRPhi", 1.f, "rphi scale of the pair-merging metric (cm)"};
    Configurable<float> cfgSigZ{"cfgSigZ", 1.f, "z scale of the pair-merging metric (cm)"};
    Configurable<float> cfgSigZNs{"cfgSigZNs", -1.f, "z scale in ns of drift time; <=0 keeps cfgSigZ in cm"};
    Configurable<float> cfgDPairClose{"cfgDPairClose", 1.f, "a sampled radius counts as close when D_pair < this"};
    Configurable<std::vector<float>> cfgCloseThrCm{"cfgCloseThrCm", {0.5f, 1.0f, 1.5f, 2.0f}, "four d_rphi thresholds (cm) for N_close/f_close"};
    Configurable<bool> cfgDoPairMergeCut{"cfgDoPairMergeCut", false, "apply the pair-merging cut to SE AND ME"};
    Configurable<float> cfgDPairCut{"cfgDPairCut", 1.f, "reject pairs with min D_pair below this (needs cfgDoPairMergeCut)"};
    Configurable<int> cfgNCloseCut{"cfgNCloseCut", 0, "reject pairs with at least this many close radii; 0 disables"};
  } pairsep;

  struct : ConfigurableGroup {
    std::string prefix = "crosspair_group";
    Configurable<bool> cfgDoCrossPairCut{"cfgDoCrossPairCut", false, "reject pairs with a photon-like crossed combination; applied to SE AND ME"};
    Configurable<float> cfgCrossMaxMeeRatio{"cfgCrossMaxMeeRatio", 0.2f, "veto: min(m_ee^cross)/q_inv below this; genuine pairs sit near 0.5, leg swaps near 0"};
    Configurable<bool> cfgDoCrossGeomCut{"cfgDoCrossGeomCut", false, "reject pairs whose leg-swapped combination is geometrically a photon: helix circles of the swapped legs touching (|gap| < cfgCrossMaxGap), |tgl+ - tgl-| < cfgCrossMaxDTgl and the legs' z at the touching point within cfgCrossMaxDz (vertex-subtracted); applied to SE AND ME. A stored leg-swapped pair (two crossed fakes of two collinear photons) always has such a combination (the true photons), a genuine pair only when the two photons are collinear and convert on one line"};
    Configurable<float> cfgCrossMaxGap{"cfgCrossMaxGap", 1.5f, "geometric crossed veto: max |circle gap| of the swapped legs (cm)"};
    Configurable<float> cfgCrossMaxDTgl{"cfgCrossMaxDTgl", 0.1f, "geometric crossed veto: max |tgl+ - tgl-| of the swapped legs"};
    Configurable<float> cfgCrossMaxDz{"cfgCrossMaxDz", 2.f, "geometric crossed veto: max |z+ - z-| of the swapped legs at the touching point (cm)"};
  } crosspair;

  // ─── Mixing ───────────────────────────────────────────────────────────────
  struct : ConfigurableGroup {
    std::string prefix = "mixing_group";
    Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
    Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
    Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
    ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis confCentBins{"confCentBins", {VARIABLE_WIDTH, 0.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f, 999.f}, "Mixing bins - centrality"};
    ConfigurableAxis confOccupancyBins{"confOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};
    // correlation block: which mixed event fills Pair/mix
    Configurable<int> cfgCoRecoMixMode{"cfgCoRecoMixMode", 0, "correlation block: 0 co-reconstructed mixed event (pairs inside the mixed artificial events), 1 standard mixed event built here from the same-event artificial events (pools in z, centrality and occupancy, no BC gap), for old vs new on the same photons"};
    Configurable<bool> cfgCoRecoSkipOneSource{"cfgCoRecoSkipOneSource", true, "correlation block: skip pairs of a mixed artificial event whose four legs all come from one source event (same-event pairs by construction)"};
    Configurable<bool> cfgCoRecoMC{"cfgCoRecoMC", false, "correlation block: MC pair classes from the leg labels carried in PCMMixLegs (CoReco/MC, and CoReco/MC/Note with cfgFillNoteQA)"};
    Configurable<int> cfgCoRecoPartnerPairing{"cfgCoRecoPartnerPairing", 0, "partner pairing of the run-2 crossed photons (correlation and pair-class blocks, every artificial event): a pair containing a crossed photon X (built in run 2 from legs of two run-1 photons) is kept only if the other photon is a partner of X. 0 off (all pairs kept), 1 partner = the two photons share a run-1 source photon through their legs (PCMMixLegs sourcePhotonUid: the other half of a run-1 photon that X took apart, e.g. the complementary crossed photon of a full swap or the second piece of a half swap), 2 as 1 plus: both photons hold a pooled free leg of the same source event with opposite charge (possible legs of one photon that run 1 did not pair). The same event has no crossed photons, so the rule only acts on the mixed event; diagnostic PairClass/Pairs/hPartner"};
  } mixing;

  // ─── Output and QA ────────────────────────────────────────────────────────
  struct : ConfigurableGroup {
    std::string prefix = "hbtanalysis_group";
    Configurable<bool> cfgDoCorrelation{"cfgDoCorrelation", true, "processPairs: run the correlation block (Pair/, PairSep/, SinglePhoton/, CoReco/ histograms)"};
    Configurable<bool> cfgFillPipe{"cfgFillPipe", false, "correlation block: book/fill the pipeline histograms CoReco/Pipe (data and MC) and CoReco/MC/Pipe (MC, with cfgCoRecoMC), for the pairs of the correlation functions after all pair cuts; the pipeline maps hold only pairs with Q_inv inside the Q axis (confQBins)"};
    Configurable<float> cfgPipeTagMaxQinv{"cfgPipeTagMaxQinv", 0.3f, "pairs with Q_inv below this value (GeV/c; keep it on a bin edge of confQBins) get the pointing tag; above it the tag is not computed and the pair is not filled into the tag maps (CoReco/Pipe/hTagMap, CoReco/MC/Pipe/hClassTagMap)"};
    Configurable<bool> cfgDo3D{"cfgDo3D", false, "book Pair/{same,mix}/CF_2D (with cfgDo2D) in place of Pair/{same,mix}/CF_1D"};
    Configurable<bool> cfgDo2D{"cfgDo2D", false, "enable 2D (qout,qinv) projection (requires cfgDo3D)"};
    Configurable<bool> cfgUseLCMS{"cfgUseLCMS", false, "measure 1D relative momentum in LCMS"};
    ConfigurableAxis confMultNTracksBins{"confMultNTracksBins", {VARIABLE_WIDTH, 0., 100., 200., 300., 400., 500., 600, 700, 800, 900, 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200., 2300., 2400., 2500., 2600., 2700., 2800., 2900., 3000., 4000., 5000.}, "N_{tracks}^{PV} |#eta|<1 bins for CF sparses"};
    ConfigurableAxis confMultFT0MBins{"confMultFT0MBins", {VARIABLE_WIDTH, 0., 25., 50., 100., 200., 400., 800., 2000., 2500., 3000., 3500., 4000., 5000., 6000., 7000., 8000., 9000., 10000., 15000., 40000., 100000., 250000.}, "FT0M amplitude bins for CF sparses"};
  } hbtanalysis;

  struct : ConfigurableGroup {
    std::string prefix = "qaflags_group";
    Configurable<bool> doPairQa{"doPairQa", true, "fill the pair QA Pair/same/QA/Before/hR1VsR2 (same event, before the pair cuts)"};
    Configurable<bool> doSinglePhotonQa{"doSinglePhotonQa", true, "fill single-photon QA histograms (SinglePhoton/Before, SinglePhotonMC/Response)"};
    Configurable<float> cfgMaxQinvForQA{"cfgMaxQinvForQA", 0.1f, "fill the pair QA only when q_inv < this value (<= 0: no limit)"};
    Configurable<float> cfgMaxQinvForFullRange{"cfgMaxQinvForFullRange", 0.3f, "fill full-range histograms only when q_inv < this value"};
    Configurable<int> cfgQaLevel{"cfgQaLevel", 2, "QA level: 0 none, 1 single-photon and pair QA, 2 also the same-sign leg separation PairSep"};
    Configurable<float> cfgMaxQinvForProcessing{"cfgMaxQinvForProcessing", 0.5, "standard mixed event (cfgCoRecoMixMode 1): skip pairs with q_inv above this before building observables"};
    // pair-class QA (processPairs with cfgDoPairClassQA)
    Configurable<bool> cfgDoPairClassQA{"cfgDoPairClassQA", false, "processPairs: run the pair-class QA block (PairClass/ histograms; same event, its own standard mixed event, co-reconstructed mixed event)"};
    Configurable<float> cfgPairClassMaxQ{"cfgPairClassMaxQ", 0.5f, "pair-class QA: max q_inv of the pairs (GeV/c)"};
    Configurable<int> cfgPairClassMixDepth{"cfgPairClassMixDepth", 10, "pair-class QA, standard mixed event: same-event artificial events kept per pool"};
    Configurable<float> cfgPairClassMixZBinWidth{"cfgPairClassMixZBinWidth", 2.f, "pair-class QA, standard mixed event: pool z-vertex bin width (cm)"};
    Configurable<float> cfgPairClassMixCentBinWidth{"cfgPairClassMixCentBinWidth", 10.f, "pair-class QA, standard mixed event: pool centrality (FT0C) bin width (%)"};
    Configurable<float> cfgPairClassMaxSharedFrac{"cfgPairClassMaxSharedFrac", 0.4f, "pair-class QA, clean flag: all four legs with TPC shared-cluster fraction below this"};
    Configurable<float> cfgNoteMaxQinv{"cfgNoteMaxQinv", 0.3f, "correlation block with cfgFillNoteQA: pairs with reconstructed q_inv above this value (GeV/c) are skipped before anything is computed for them (CPU); 0.3 is the range of the CoReco/MC/Note histograms. The true-pair histogram is in true q_inv: its last bins lose the pairs reconstructed above the limit"};
    Configurable<bool> cfgFillNoteQA{"cfgFillNoteQA", false, "book/fill the analysis-note QA, every histogram with a before/after-cut axis: PairClass/Note (pair-class QA: classes vs centrality and occupancy, duplicate cluster counts, shared-cluster scan, photon kinematics per class), CoReco/Note and CoReco/MC/Note (correlation block: classes before/after the pair cuts, q-component response, true pairs vs true q), MC/TruthAO2D/Note (processGenerated: fiducial true pairs per stage)"};
    Configurable<bool> cfgFillHalfSwapQA{"cfgFillHalfSwapQA", false, "book/fill the half-swap QA PairClass/HalfSwap (pair-class QA, MC): hStray = MC relation of the two stray legs of a half swap (class 4) to the shared true photon; hCrossBoth = photon-likeness of BOTH leg-swapped combinations of a pair per pair class (design of an asymmetric crossed veto: one combination a good photon, the other not)"};
  } qaflags;

  struct : ConfigurableGroup {
    std::string prefix = "mctruth_group";
    Configurable<float> cfgMCMaxQinv{"cfgMCMaxQinv", 0.3f, "processGenerated: max true q_inv of the true-photon pairs (GeV/c); <= 0 no limit"};
    Configurable<float> cfgMCMinKt{"cfgMCMinKt", 0.0f, "processGenerated: min true k_T of the true-photon pairs (GeV/c); <= 0 no limit"};
    Configurable<float> cfgMCMaxKt{"cfgMCMaxKt", 0.7f, "processGenerated: max true k_T of the true-photon pairs (GeV/c); <= 0 no limit"};
    Configurable<float> cfgMCMinV0Pt{"cfgMCMinV0Pt", 0.1f,
                                     "min pT for true photons in truth-efficiency loop (GeV/c); "
                                     "0 = fall back to pcmcuts.cfgMinPtV0"};
    Configurable<float> cfgMCMinLegPt{"cfgMCMinLegPt", 0.0f, "min pT for true e^{+}/e^{-} legs in truth-efficiency loop (GeV/c);"};
  } mctruth;

  ConfigurableAxis confQBins{"confQBins", {60, 0, +0.3f}, "q bins for output histograms"};
  ConfigurableAxis confKtBins{"confKtBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75}, "kT bins"};
  ConfigurableAxis confPtBins{"confPtBins", {100, 0.f, 2.f}, "pT bins (GeV/c)"};
  ConfigurableAxis confEtaBins{"confEtaBins", {80, -0.8f, 0.8f}, "eta bins"};
  ConfigurableAxis confPhiBins{"confPhiBins", {90, 0.f, o2::constants::math::TwoPI}, "phi bins (rad)"};
  ConfigurableAxis confRBins{"confRBins", {100, 0.f, 100.f}, "conversion radius bins (cm)"};
  ConfigurableAxis confZConvBins{"confZConvBins", {200, -100.f, 100.f}, "conversion z (cm)"};
  ConfigurableAxis confOccupancyQA{"confOccupancyQA", {100, 0.f, 50000.f}, "occupancy"};

  // ─── CCDB ─────────────────────────────────────────────────────────────────
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  Configurable<std::string> cfgCcdbUrl{"cfgCcdbUrl", "http://alice-ccdb.cern.ch", "CCDB url"};
  Configurable<float> cfgBzOverrideT{"cfgBzOverrideT", -999.f, "Bz in Tesla; used instead of CCDB if > -100"};

  AxisSpec axisKt{std::vector<double>{}};
  AxisSpec axisQinv{std::vector<double>{}};
  AxisSpec axisQabsLcms{std::vector<double>{}};
  AxisSpec axisQout{std::vector<double>{}};
  AxisSpec axisQside{std::vector<double>{}};
  AxisSpec axisQlong{std::vector<double>{}};
  AxisSpec axisPt{std::vector<double>{}};
  AxisSpec axisEta{std::vector<double>{}};
  AxisSpec axisPhi{std::vector<double>{}};
  AxisSpec axisR{std::vector<double>{}};
  AxisSpec axisZConv{std::vector<double>{}};
  AxisSpec axisOccupancy{std::vector<double>{}};
  [[nodiscard]] AxisSpec makeAxisMultNTracks() const { return AxisSpec{hbtanalysis.confMultNTracksBins, "N_{tracks}^{PV}, |#eta| < 1"}; }
  [[nodiscard]] AxisSpec makeAxisMultFT0M() const { return AxisSpec{hbtanalysis.confMultFT0MBins, "mult. FT0M (amplitude)"}; }

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry fRegistryCF{"cf", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry fRegistryPairQA{"pairQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry fRegistryPairMC{"pairMC", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry fRegistryTruthMC{"truthMC", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  EMPhotonEventCut fEMEventCut;
  V0PhotonCut fV0PhotonCut;

  float mBzT{0.f};           // signed, Tesla
  float mVDriftCmPerNs{0.f}; // TPC drift velocity from CCDB; 0 = unavailable
  int mRunNumber{0};
  bool isMC = false;
  double mPipeQinvMax = 0.;    // upper edge of the Q_inv axis (axisQinv): the pipeline maps hold pairs below it only
  bool mDoCorrelation = false; // processPairs and cfgDoCorrelation
  bool mDoPairClassQA = false; // processPairs and cfgDoPairClassQA
  int ndf = 0;                 // DF counter: makes the keys of the standard mixed event unique over DFs
  pairutil::V0PhotonClassSelection mPhotonClassSelA{};
  pairutil::V0PhotonClassSelection mPhotonClassSelB{};

  struct QaGates {
    bool pairQa{true};
    bool singlePhotonQa{true};
    bool pairSepQa{true};
    bool legPropagation{true}; // something reads the propagated legs: PairSep QA or the pair-merging cut
  } mQa;

  // standard mixed event of the correlation block (cfgCoRecoMixMode 1)
  std::vector<float> ztxBinEdges;
  std::vector<float> centBinEdges;
  std::vector<float> occBinEdges;
  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, PhotonWithLegs>;
  std::shared_ptr<MyEMH> emh1;
  std::unordered_set<int> usedPhotonIdsPerCol;

  // pair-class QA: standard mixed event pools
  std::map<int, std::deque<std::vector<PcPhoton>>> mPcPools;
  float mPcCent = -1.f;               // centrality FT0C
  float mPcNLegs = -1.f;              // legs (PCMMixLegs rows), for hPairDens
  float mPcOcc = -1.f, mPcVtxZ = 0.f; // occupancy and z vertex (analysis-note QA)
  std::vector<int8_t> mPhotonPass;    // photon cut per photon of the DF: -1 not yet evaluated, 0 rejected, 1 selected

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

  using FilteredMyMCCollisions = soa::Filtered<MyCollisionsMC>;

  void init(InitContext& context)
  {

    buildConfigurableAxes();
    mPipeQinvMax = axisQinv.binEdges.back();
    const bool doPairs = context.mOptions.get<bool>("processPairs");
    if (doPairs && !hbtanalysis.cfgDoCorrelation.value && !qaflags.cfgDoPairClassQA.value) {
      LOGF(fatal, "processPairs needs cfgDoCorrelation or cfgDoPairClassQA: with both off it fills nothing");
    }
    mDoCorrelation = doPairs && hbtanalysis.cfgDoCorrelation.value;
    mDoPairClassQA = doPairs && qaflags.cfgDoPairClassQA.value;
    isMC = context.mOptions.get<bool>("processGenerated");
    if (mDoCorrelation && isMC) {
      LOGF(fatal, "processPairs with cfgDoCorrelation sets isMC from cfgCoRecoMC, processGenerated needs it on: run processGenerated in its own photonhbt instance");
    }
    if (mDoCorrelation) {
      isMC = mixing.cfgCoRecoMC.value;
    }
    if (mDoPairClassQA) {
      addPairClassMixHistograms();
    }
    if (pairsep.cfgCloseThrCm.value.size() != 4) { // o2-linter: disable=magic-number (four thresholds)
      LOGF(fatal, "cfgCloseThrCm must contain exactly four thresholds, got %zu",
           pairsep.cfgCloseThrCm.value.size());
    }
    if (pairsep.cfgCloseThrCm.value.front() <= 0.f ||
        !std::is_sorted(pairsep.cfgCloseThrCm.value.begin(), pairsep.cfgCloseThrCm.value.end())) {
      LOGF(fatal, "cfgCloseThrCm must be positive and sorted in ascending order");
    }
    const int qaLevel = qaflags.cfgQaLevel.value;
    mQa.pairQa = qaflags.doPairQa.value && qaLevel >= 1;                 // o2-linter: disable=magic-number (QA set-up)
    mQa.singlePhotonQa = qaflags.doSinglePhotonQa.value && qaLevel >= 1; // o2-linter: disable=magic-number (QA set-up)
    mQa.pairSepQa = pairsep.cfgDoPairSepQA.value && qaLevel >= 2;        // o2-linter: disable=magic-number (QA set-up)
    mQa.legPropagation = mQa.pairSepQa || pairsep.cfgDoPairMergeCut.value;
    LOGF(info, "photonhbt QA level %d -> pairQA %d, singlePhotonQA %d, pairSep %d",
         qaLevel, mQa.pairQa, mQa.singlePhotonQa, mQa.pairSepQa);
    LOGF(info, "photonhbt: processPairs blocks -> correlation %d, pair-class QA %d", mDoCorrelation, mDoPairClassQA);
    mRunNumber = 0;
    parseBins(mixing.confVtxBins, ztxBinEdges);
    parseBins(mixing.confCentBins, centBinEdges);
    parseBins(mixing.confOccupancyBins, occBinEdges);
    emh1 = std::make_shared<MyEMH>(mixing.ndepth);
    DefineEMEventCut();
    DefinePCMCut();
    mPhotonClassSelA = pairutil::buildV0PhotonClassSelection(photonclassA);
    mPhotonClassSelB = pairutil::buildV0PhotonClassSelection(photonclassB);
    addhistograms();
    if (mDoCorrelation) {
      addCoRecoHistograms();
    }
    if (mDoCorrelation && hbtanalysis.cfgFillPipe.value) {
      addPipeHistograms();
    }
    ccdb->setURL(cfgCcdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    const float worst = selfTestLegHelix();
    constexpr float kSelfTestTol = 1.e-4f;
    if (worst < kSelfTestTol) {
      LOGF(info, "photonhbt: leg-helix origin-limit self test PASSED (max |dphi| = %.2e rad)", worst);
    } else {
      LOGF(fatal,
           "photonhbt: leg-helix origin-limit self test FAILED, max |dphi| = %.3e rad. "
           "The helix-centre sign convention or the branch choice is wrong; every "
           "leg-separation variable would be mirrored. Refusing to run.",
           worst);
    }
    LOGF(info, "photonhbt: leg propagation sampled at %zu radii, %.0f-%.0f cm in %.0f cm steps", kPhiStarRadiiM.size(), kPhiStarRadiiM.front() * kCmPerM, kPhiStarRadiiM.back() * kCmPerM, (kPhiStarRadiiM.size() > 1) ? (kPhiStarRadiiM[1] - kPhiStarRadiiM[0]) * kCmPerM : 0.f);
    LOGF(info, "photonhbt: PairSep QA %s", pairsep.cfgDoPairSepQA.value ? "on" : "off");
    LOGF(info, "photonhbt r-phi closeness thresholds: %.2f, %.2f, %.2f, %.2f cm", pairsep.cfgCloseThrCm.value[0], pairsep.cfgCloseThrCm.value[1], pairsep.cfgCloseThrCm.value[2], pairsep.cfgCloseThrCm.value[3]);
    LOGF(info, "photonhbt pair merging: metric sigma_rphi = %.2f cm, sigma_z = %.2f cm (or %.1f ns if > 0), close at D_pair < %.2f", pairsep.cfgSigRPhi.value, pairsep.cfgSigZ.value, pairsep.cfgSigZNs.value, pairsep.cfgDPairClose.value);
    LOGF(info,
         "photonhbt pair separation: all merging observables (d_rphi^min, D_pair^min, "
         "N_close, f_close) use SAME-SIDE radii only; the pure r-phi quantities are "
         "selected by d_rphi, not by d_3D");
    LOGF(info, "photonhbt pair merging cut: %s, D_pair > %.2f, N_close < %d (0 = off)", pairsep.cfgDoPairMergeCut.value ? "ON" : "off", pairsep.cfgDPairCut.value, pairsep.cfgNCloseCut.value);
    LOGF(info, "photonhbt: crossed-pair cut %s (min(m_ee)/q_inv > %.2f), geometric cut %s (|gap| < %.2f cm, |dtgl| < %.3f, |dz| < %.2f cm)", crosspair.cfgDoCrossPairCut.value ? "ON" : "off", crosspair.cfgCrossMaxMeeRatio.value, crosspair.cfgDoCrossGeomCut.value ? "ON" : "off", crosspair.cfgCrossMaxGap.value, crosspair.cfgCrossMaxDTgl.value, crosspair.cfgCrossMaxDz.value);
    LOGF(info, "photonhbt: q_inv gates -- QA < %.3f, FullRange < %.3f, PairSep < %.3f", qaflags.cfgMaxQinvForQA.value, qaflags.cfgMaxQinvForFullRange.value, pairsep.cfgPairSepMaxQinv.value);
  }

  void buildConfigurableAxes()
  {
    axisKt = AxisSpec{confKtBins, "k_{T} (GeV/c)"};
    axisQinv = AxisSpec{confQBins, "q_{inv} (GeV/c)"};
    axisQabsLcms = AxisSpec{confQBins, "|#bf{q}|^{LCMS} (GeV/c)"};
    axisQout = AxisSpec{confQBins, "q_{out} (GeV/c)"};
    axisQside = AxisSpec{confQBins, "q_{side} (GeV/c)"};
    axisQlong = AxisSpec{confQBins, "q_{long} (GeV/c)"};
    axisPt = AxisSpec{confPtBins, "p_{T} (GeV/c)"};
    axisEta = AxisSpec{confEtaBins, "#eta"};
    axisPhi = AxisSpec{confPhiBins, "#phi (rad)"};
    axisR = AxisSpec{confRBins, "R_{conv} (cm)"};
    axisZConv = AxisSpec{confZConvBins, "z_{conv} (cm)"};
    axisOccupancy = AxisSpec{confOccupancyQA, "occupancy"};
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }
    mRunNumber = collision.runNumber();

    auto vd = ccdb->getForRun<o2::tpc::VDriftCorrFact>("TPC/Calib/VDriftTgl", mRunNumber);
    if (vd != nullptr) {
      mVDriftCmPerNs = vd->refVDrift * vd->corrFact * 1e-3f; // o2-linter: disable=magic-number (cm/us -> cm/ns)
      LOGF(info, "photonhbt: run %d, TPC vdrift = %.6f cm/ns (%.4f cm/us); 1 cm in z = %.0f ns",
           mRunNumber, mVDriftCmPerNs, 1e3f * mVDriftCmPerNs,
           (mVDriftCmPerNs > 0.f) ? 1.f / mVDriftCmPerNs : -1.f);
    } else {
      mVDriftCmPerNs = 0.f;
      LOGF(warn, "photonhbt: no TPC VDrift object for run %d -- the z scale stays in cm", mRunNumber);
    }

    if (cfgBzOverrideT.value > -100.f) { // o2-linter: disable=magic-number (number in case B-field is overridden in case not fetched from CCDB)
      mBzT = cfgBzOverrideT.value;
      return;
    }
    auto grpmag = ccdb->getForRun<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", mRunNumber);
    mBzT = 0.1f * static_cast<float>(grpmag->getNominalL3Field());
    LOGF(info, "photonhbt: run %d, Bz = %.2f T", mRunNumber, mBzT);
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
    fEMEventCut.SetTrackOccupancyRange(eventcuts.cfgTrackOccupancyMin, eventcuts.cfgTrackOccupancyMax);
    fEMEventCut.SetFT0COccupancyRange(eventcuts.cfgFT0COccupancyMin, eventcuts.cfgFT0COccupancyMax);
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
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(pcmcuts.cfgMinNCrossedRowsOverFindable);
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

  void addhistograms()
  {
    addPairCFHistograms();
    addPairSepHistograms();
    addSinglePhotonQAHistograms();
    addPairQAHistograms();

    if (isMC) {
      addSinglePhotonMCHistograms();
      addTruthMCHistograms();
    }
  }

  // ─── Pair/: correlation functions and FullRange (fRegistryCF) ─────────────
  void addPairCFHistograms()
  {
    const AxisSpec axisMultNTracks = makeAxisMultNTracks();
    const AxisSpec axisMultFT0M = makeAxisMultFT0M();
    if (hbtanalysis.cfgDo3D) {
      if (hbtanalysis.cfgDo2D) {
        fRegistryCF.add("Pair/same/CF_2D", "diphoton correlation 2D (qout,qinv)", kTHnSparseD, {axisQout, axisQinv, axisMultNTracks, axisMultFT0M, axisKt}, true);
      }
    } else {
      fRegistryCF.add("Pair/same/CF_1D", hbtanalysis.cfgUseLCMS ? "diphoton correlation 1D LCMS" : "diphoton correlation 1D (qinv)", kTH2D, {hbtanalysis.cfgUseLCMS ? axisQabsLcms : axisQinv, axisKt}, true);
    }
    fRegistryCF.add("Pair/same/FullRange/hQinvVsOccupancy", "q_{inv} vs occupancy;occupancy;q_{inv} (GeV/c)", kTH2D, {axisOccupancy, axisQinv}, true);

    fRegistryCF.add("Pair/same/CF_QLcms_Qinv", "diphoton CF |q|_{LCMS} vs. q_{inv}", kTHnSparseD, {axisQabsLcms, axisQinv, axisKt}, true);

    fRegistryCF.addClone("Pair/same/", "Pair/mix/");
  }

  // ─── PairSep/: same-sign leg separation (fRegistryCF) ─────────────────────
  void addPairSepHistograms()
  {
    const AxisSpec axisQcal4{{0.00, 0.01, 0.02, 0.04, 0.06, 0.10, 0.15, 0.20, 0.30, 0.50}};

    std::vector<std::string> dirs = {"PairSep/same/", "PairSep/mix/"};
    for (const auto& sm : dirs) {
      fRegistryCF.add((sm + "hSparse_dRPhiMin_fClose_Thr_Qinv").c_str(), "pure r#varphi;d_{r#varphi}^{min} (cm);f_{close}^{r#varphi};threshold index;q_{inv} (GeV/c)", kTHnSparseF, {{40, 0.f, 10.f}, {18, -0.028f, 1.028f}, {4, -0.5f, 3.5f}, axisQcal4}, true);
      fRegistryCF.add((sm + "hSparse_DPair_NClose_Qinv").c_str(), "scaled pair distance;D_{pair}^{min};N_{close};q_{inv} (GeV/c)", kTHnSparseF, {{40, 0.f, 10.f}, {18, -0.5f, 17.5f}, {50, 0.f, 0.5f}}, true);
    }
  }

  // ─── SinglePhoton/: selected photons of the same event (fRegistryPairQA) ──
  void addSinglePhotonQAHistograms()
  {
    fRegistryPairQA.add("SinglePhoton/Before/hEtaVsPhiPt", "acceptance;#phi (rad);#eta", kTH3D, {axisPhi, axisEta, axisPt}, true);
    fRegistryPairQA.add("SinglePhoton/Before/hRVsZConvPt", "R_{conv} vs z_{conv};z_{conv} (cm);R_{conv} (cm)", kTH3D, {axisZConv, axisR, axisPt}, true);
  }

  // ─── Pair/same/QA/: same-event pairs before the pair cuts (fRegistryPairQA)
  void addPairQAHistograms()
  {
    //  Conversion point
    fRegistryPairQA.add("Pair/same/QA/Before/hR1VsR2", "R_{1} vs R_{2};R_{1} (cm);R_{2} (cm)", kTH2D, {axisR, axisR}, true);
  }

  // ─── SinglePhotonMC/: single-photon response (fRegistryPairMC) ────────────
  void addSinglePhotonMCHistograms()
  {
    fRegistryPairMC.add("SinglePhotonMC/Response/hResponse", "true reconstructed photons: response (MC)", kTHnSparseF, {{40, 0.f, 4.f, "#it{p}_{T}^{true} (GeV/#it{c})"}, {60, 0.3f, 1.5f, "#it{E}^{rec} / #it{E}^{true}"}, {40, 0.f, 0.04f, "angle(rec, true) (rad)"}, {18, 0.f, 90.f, "#it{R}_{conv} (cm)"}, {21, -0.05f, 1.f, "hardest photon radiated by a leg: #it{E} / #it{E}^{true} (first bin: none)"}, {2, -0.5f, 1.5f, "that radiated photon reconstructed in this event"}, {3, -0.5f, 2.5f, "this photon radiated by a conversion e^{#pm}: 0 no, 1 yes and parent #gamma lost, 2 yes and parent reconstructed"}, {21, -0.05f, 1.f, "#it{E}^{true} / #it{E}(parent #gamma) (first bin: n/a)"}, {51, -0.005f, 0.25f, "angle(this, parent #gamma) MC (rad) (first bin: n/a)"}}, true);
    fRegistryPairMC.add("SinglePhotonMC/Response/hGenPt", "generated primary photons of the selected collisions, |#eta| < cfgMaxEtaV0;#it{p}_{T} (GeV/#it{c})", kTH1F, {{80, 0.f, 4.f}}, true);
  }

  // ─── MC/TruthAO2D/: generated true-photon pairs (fRegistryTruthMC) ────────
  void addTruthMCHistograms()
  {
    const AxisSpec axRconv{180, 0.f, 90.f, "R_{conv}^{true} (cm)"};

    // ─── Rconv waterfall
    fRegistryTruthMC.add("MC/TruthAO2D/hRconv1_vs_Rconv2_truthConverted", "denominator: R_{conv,1} vs R_{conv,2};R_{conv,1}^{true} (cm);R_{conv,2}^{true} (cm)", kTH2D, {axRconv, axRconv}, true);
    fRegistryTruthMC.add("MC/TruthAO2D/hRconv1_vs_Rconv2_bothPhotonsSelected", "numerator: R_{conv,1} vs R_{conv,2};R_{conv,1}^{true} (cm);R_{conv,2}^{true} (cm)", kTH2D, {axRconv, axRconv}, true);
    if (qaflags.cfgFillNoteQA.value) {
      fRegistryTruthMC.add("MC/TruthAO2D/Note/hTruePairs", "fiducial true pairs of converted photons per stage", kTHnSparseF, {{4, -0.5f, 3.5f, "stage: 0 converted, 1 all four legs tracked in this collision, 2 both photons built, 3 both photons selected"}, {60, 0.f, 0.3f, "#it{q}_{inv}^{true} (GeV/#it{c})"}, {std::vector<double>{0., 0.2, 0.4, 0.6, 1.0, 2.0}, "#it{k}_{T}^{true} (GeV/#it{c})"}, {std::vector<double>{0., 10., 20., 30., 40., 50., 60., 80., 100.}, "centrality (%)"}, {std::vector<double>{0., 0.5, 1., 2., 3., 5., 7.5, 10., 15., 20., 30., 50., 100., 1000.}, "distance of the two true conversion points (cm)"}}, true);
    }
  }

  // ─── PairClass/: pair-class QA block (fRegistryPairMC) ────────────────────
  void addPairClassMixHistograms()
  {
    if (qaflags.cfgFillNoteQA.value) {
      addNotePairClassHistograms();
    }
    const int nq = std::max(1, static_cast<int>(std::lround(qaflags.cfgPairClassMaxQ.value * 1000.f)));
    const AxisSpec axSample{4, -0.5f, 3.5f, "0 same event, 1 standard mixed event, 2 co-reconstructed mixed event, 3 co-reco ME, one source event"};
    const AxisSpec axClass{21, -0.5f, 20.5f, "pair class: 13 true distinct, 14 full swap, 4 half swap, 1-3 duplicates, 16 bremsstrahlung (both true), 5 bremsstrahlung with a fake, 18 photons of the same mother (pi0/eta, both true), 19 swap of two such photons, 20 same mother with a fake, 17 related swap, 15 no label, 6-12 other fakes"};
    const AxisSpec axQ{nq, 0.f, qaflags.cfgPairClassMaxQ.value, "#it{q}_{inv} (GeV/#it{c})"};
    const AxisSpec axOrigin{4, -0.5f, 3.5f, "stored photons: 0 both SVertexer originals, 1 one finder original, 2 both finder, 3 a crossed pairing"};
    const AxisSpec axClean{2, -0.5f, 1.5f, "all four legs below the TPC shared-cluster cut"};
    fRegistryPairMC.add("PairClass/Pairs/hPairClass", "photon pairs by MC class", kTHnSparseF, {axSample, axClass, axQ, axOrigin, axClean});
    const AxisSpec axPass{3, -0.5f, 2.5f, "photons of the pair passing the analysis photon cut (0, 1, 2)"};
    fRegistryPairMC.add("PairClass/Pairs/hPairTag", "photon pairs with q_{inv} < 100 MeV: bremsstrahlung pointing tag vs the analysis cuts", kTHnSparseF, {axSample, axClass, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {std::vector<double>{0., 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 10.}, "|#Delta#it{E}| (GeV)"}, {std::vector<double>{0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1., 10.}, "#it{k}_{T} (GeV/#it{c})"}, {std::vector<double>{0., 10., 20., 40., 60., 80., 100.}, "centrality FT0C (%)"}, {5, -0.5f, 4.5f, "pointing d_{min}: 0 < 0.1, 1 0.1-0.3, 2 0.3-1, 3 > 1 cm, 4 no allowed point"}, axPass, {2, -0.5f, 1.5f, "same-sign leg separation d_{min} > cfgDPairCut"}, {std::vector<double>{0., 0.04, 0.1, 0.2, 0.3, 10.}, "|#it{q}_{out}| (LCMS) (GeV/#it{c})"}});
    fRegistryPairMC.add("PairClass/Pairs/hBremsFake", "bremsstrahlung with a fake photon (class 5): the fake photon and its foreign leg", kTHnSparseF, {axSample, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {4, -0.5f, 3.5f, "fake photon: 0 the radiated photon (one of its legs + a foreign leg), 1 the parent photon (radiating e^{#pm} + a foreign leg), 2 both, 3 no radiated-photon leg (leg is the daughter of a leg)"}, {6, -0.5f, 5.5f, "foreign leg of the fake: 0 conversion e^{#pm} of an unrelated photon, 1 from the same cascade, 2 e^{#pm} not from a conversion, 3 not an electron, 4 no MC label, 5 n/a"}, {std::vector<double>{0., 10., 20., 40., 60., 80., 100.}, "centrality FT0C (%)"}, {5, -0.5f, 4.5f, "pointing d_{min}: 0 < 0.1, 1 0.1-0.3, 2 0.3-1, 3 > 1 cm, 4 no allowed point"}, axPass, {2, -0.5f, 1.5f, "same-sign leg separation d_{min} > cfgDPairCut"}, {3, -0.5f, 2.5f, "the fake photon is: 0 an SVertexer photon, 1 a finder photon, 2 a crossed pairing (run 2)"}});
    fRegistryPairMC.add("PairClass/Pairs/hPartner", "pairs with a run-2 crossed photon: partner relation (photonhbt cfgCoRecoPartnerPairing) before the rule", kTHnSparseF, {axSample, axClass, axQ, {3, 0.5f, 3.5f, "1 shared run-1 source photon, 2 free legs of one source event with opposite charge, 3 unrelated"}, {2, -0.5f, 1.5f, "0 removed, 1 kept"}});
    fRegistryPairMC.add("PairClass/Pairs/hDEQinv", "photon pairs: |#Delta E| vs q_{inv} by MC class", kTHnSparseF, {axSample, axClass, {100, 0.f, 0.5f, "#it{q}_{inv} (GeV/#it{c})"}, {100, 0.f, 0.5f, "|#it{E}_{1} #minus #it{E}_{2}| (GeV)"}, {std::vector<double>{0., 0.2, 0.4, 0.6, 0.9, 1.5, 3.}, "#it{k}_{T} (GeV/#it{c})"}});
    fRegistryPairMC.add("PairClass/Pairs/hFreeLegPair", "photon pairs by MC class: photons built with a pooled free leg", kTHnSparseF, {axSample, axClass, axQ, {3, -0.5f, 2.5f, "photons of the pair with a free leg: 0, 1, 2"}, {3, -0.5f, 2.5f, "crossed photons (run-2 pairings) in the pair: 0, 1, 2"}});
    fRegistryPairMC.add("PairClass/Pairs/hHalfSwapLegs", "half swaps: run-1 origin of the two legs of the shared true photon", kTHnSparseF, {axSample, axQ, {4, -0.5f, 3.5f, "legs of the shared photon: 0 both free legs, 1 one stored photon (taken apart), 2 two different stored photons, 3 one free one stored"}, {3, -0.5f, 2.5f, "photons of the pair with a free leg: 0, 1, 2"}});
    if (qaflags.cfgFillHalfSwapQA.value) {
      const AxisSpec axQ100{20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"};
      const AxisSpec axStrayA{7, -0.5f, 6.5f, "stray leg of photon 1: 0 conversion e^{#pm} of an unrelated photon, 1 other unrelated e^{#pm}, 2 same cascade as the shared photon, 3 same mother as the shared photon, 4 same physical primary otherwise, 5 not an electron, 6 no MC label"};
      const AxisSpec axStrayB{7, -0.5f, 6.5f, "stray leg of photon 2 (categories as photon 1)"};
      const AxisSpec axCent6{std::vector<double>{0., 10., 20., 40., 60., 80., 100.}, "centrality FT0C (%)"};
      fRegistryPairMC.add("PairClass/HalfSwap/hStray", "half swaps (class 4): MC relation of the two stray legs to the shared true photon", kTHnSparseF, {axSample, axQ100, axStrayA, axStrayB, {3, -0.5f, 2.5f, "the two stray legs: 0 unrelated, 1 same mother, 2 same physical primary"}, axCent6, axPass, {3, -0.5f, 2.5f, "stray legs that are pooled free legs (0, 1, 2)"}});
      const std::vector<double> sBins{0., 0.25, 0.5, 0.75, 1., 1.5, 2., 3., 5., 10., 1000.};
      fRegistryPairMC.add("PairClass/HalfSwap/hCrossBoth", "photon pairs with q_{inv} < 100 MeV: photon-likeness s of both leg-swapped combinations by MC class (s = max of |circle gap|, |#Delta z|, |#Delta tan#lambda| in units of cfgCrossMaxGap / Dz / DTgl; s < 1 = inside the geometric crossed veto window)", kTHnSparseF, {axSample, axClass, axQ100, {sBins, "s of the more photon-like leg-swapped combination"}, {sBins, "s of the other leg-swapped combination"}, axPass, {2, -0.5f, 1.5f, "same-sign leg separation d_{min} > cfgDPairCut"}});
    }
    fRegistryPairMC.add("PairClass/Pairs/hKtQ", "photon pairs: k_{T} vs q_{inv} by MC class", kTHnSparseF, {axSample, axClass, axQ, {100, 0.f, 2.f, "#it{k}_{T} (GeV/#it{c})"}});
    fRegistryPairMC.add("PairClass/Pairs/hQinvQout", "photon pairs: q_{inv} vs |q_{out}| (LCMS) by MC class", kTHnSparseF, {axSample, axClass, axQ, {60, 0.f, 0.3f, "|#it{q}_{out}| (LCMS) (GeV/#it{c})"}, {std::vector<double>{0., 0.2, 0.4, 0.6, 0.9, 1.5, 3.}, "#it{k}_{T} (GeV/#it{c})"}});
    fRegistryPairMC.add("PairClass/Pairs/hPairOpening", "photon pairs with q_{inv} < 50 MeV: orientation of the opening by MC class", kTHnSparseF, {axSample, axClass, {10, 0.f, 0.05f, "#it{q}_{inv} (GeV/#it{c})"}, {30, 0.f, 0.15f, "|#it{q}_{side}| (LCMS) (GeV/#it{c})"}, {30, 0.f, 0.15f, "|#it{q}_{long}| (LCMS) (GeV/#it{c})"}, {18, 0.f, 1.5708f, "#psi = atan(|#Delta#eta| / |#Delta#varphi|) (rad)"}});
    fRegistryPairMC.add("PairClass/Pairs/hSoftHard", "photon pairs with q_{inv} < 100 MeV: soft vs hard photon (reconstructed) by MC class", kTHnSparseF, {axSample, axClass, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {20, 0.f, 1.f, "#it{E}_{soft} / #it{E}_{hard}"}, {60, -60.f, 60.f, "#it{R}_{xy}(soft) #minus #it{R}_{xy}(hard) (cm)"}, {18, 0.f, 90.f, "#it{R}_{xy}(hard) (cm)"}, {18, 0.f, 1.5708f, "#psi (rad)"}});
    fRegistryPairMC.add("PairClass/Pairs/hPairRes", "pairs of two true photons, q_{inv} < 100 MeV: resolution by MC class", kTHnSparseF, {axSample, axClass, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {40, -0.02f, 0.02f, "#it{q}_{inv}^{rec} #minus #it{q}_{inv}^{true} (GeV/#it{c})"}, {40, -0.1f, 0.1f, "#Delta#it{E}^{rec} #minus #Delta#it{E}^{true} (GeV)"}, {std::vector<double>{0., 10., 20., 40., 80., 100.}, "centrality FT0C (%)"}});
    fRegistryPairMC.add("PairClass/Pairs/hPairDens", "photon pairs by MC class vs the number of legs in the artificial event", kTHnSparseF, {axSample, axClass, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {std::vector<double>{0., 25., 50., 75., 100., 125., 150., 200., 250., 300., 400., 500., 600., 800., 1000., 1500.}, "legs in the artificial event"}, {std::vector<double>{0., 10., 20., 30., 40., 50., 60., 80., 100.}, "centrality FT0C (%)"}});
    fRegistryPairMC.add("PairClass/Pairs/hPairDist", "photon pairs: conversion-point distances by MC class", kTHnSparseF, {axSample, axClass, {50, 0.f, 0.5f, "#it{q}_{inv} (GeV/#it{c})"}, {100, 0.f, 50.f, "reconstructed: distance of the two conversion points (cm)"}, {104, -2.f, 50.f, "true: largest distance of the legs' production vertices (cm)"}});
    fRegistryPairMC.add("PairClass/Pairs/hPairQ", "photon pairs: fitted q, four-leg invariant mass, true q by MC class", kTHnSparseF, {axSample, axClass, {100, 0.f, 0.2f, "#it{q}_{inv} of the fitted photons (GeV/#it{c})"}, {100, 0.f, 0.2f, "#it{M} of the four legs (GeV/#it{c}^{2})"}, {100, 0.f, 0.2f, "#it{q}_{inv} of the two true photons (GeV/#it{c})"}, {6, 0.f, 1.5f, "#it{k}_{T} (GeV/#it{c})"}});
    fRegistryPairMC.add("PairClass/Pairs/hDupPieces", "duplicates (classes 1, 2, 3, 12): the two pieces of the shared MC particle", kTHnSparseF, {axSample, {4, -0.5f, 3.5f, "class: 0 = 1 (same photon twice), 1 = 2, 2 = 3, 3 = 12 (same track)"}, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {40, 0.f, 0.2f, "|tan#lambda_{1} #minus tan#lambda_{2}| of the two pieces"}, {20, 0.f, 1.f, "|#it{p}_{T,1} #minus #it{p}_{T,2}| / max(#it{p}_{T}) of the two pieces"}, {3, -0.5f, 2.5f, "charge slots: 0 same (both +/both #minus), 1 opposite, 2 no shared particle found"}, axOrigin, {3, -0.5f, 2.5f, "pieces with ITS (of 2)"}, {3, -0.5f, 2.5f, "photons of the pair passing the analysis photon cut (0, 1, 2)"}});
    fRegistryPairMC.add("PairClass/Pairs/hTrueOrigin", "true distinct pairs (class 13, MC), same event, both photons primary: which collision the photons come from (true origin z on the beam line vs reconstructed vertex; needs cfgFillNoteQA)", kTHnSparseF, {{15, 0.f, 0.3f, "q_{inv} (GeV/c)"}, {std::vector<double>{0., 1000., 3000., 7000., 12000., 20000., 35000., 1.e6}, "occupancy"}, {4, -0.5f, 3.5f, "0 both this collision, 1 one other, 2 both the same other, 3 two others"}, {2, -0.5f, 1.5f, "0 opposite eta sides, 1 same side"}});
    fRegistryPairMC.add("PairClass/Pairs/hTrueOriginDz", "photons of true distinct pairs (class 13, MC), same event, both primary: true origin z minus reconstructed vertex z (needs cfgFillNoteQA)", kTHnSparseF, {{120, -30.f, 30.f, "z_{origin}^{true} - z_{vtx} (cm)"}, {8, -0.8f, 0.8f, "#eta"}, {std::vector<double>{0., 1000., 3000., 7000., 12000., 20000., 35000., 1.e6}, "occupancy"}});
    fRegistryPairMC.add("PairClass/Pairs/hTrueEta12", "true distinct pairs (class 13, MC), q_inv < 0.25 GeV/c: eta of the two photons (random order) vs occupancy (needs cfgFillNoteQA)", kTHnSparseF, {axSample, {8, -0.8f, 0.8f, "#eta of one photon"}, {8, -0.8f, 0.8f, "#eta of the other photon"}, {std::vector<double>{0., 1000., 3000., 7000., 12000., 20000., 35000., 1.e6}, "occupancy"}});
    fRegistryPairMC.add("PairClass/Pairs/hTrueAngleEvt", "true distinct pairs (class 13, MC), q_inv < 0.25 GeV/c: |d eta| vs event properties (needs cfgFillNoteQA)", kTHnSparseF, {axSample, {8, 0.f, 0.8f, "|#Delta#eta|"}, {5, 0.f, 10.f, "|z_{vtx}| (cm)"}, {std::vector<double>{0., 1000., 3000., 7000., 12000., 1.e6}, "occupancy"}, {4, -0.8f, 0.8f, "mean #eta of the pair"}});
    fRegistryPairMC.add("PairClass/Pairs/hTrueAngle", "true distinct pairs (class 13, MC): angular distance of the two photons", kTHnSparseF, {axSample, {5, 0.f, 0.25f, "#it{q}_{inv} (GeV/#it{c})"}, {16, 0.f, 0.8f, "|#Delta#eta|"}, {16, 0.f, 0.8f, "|#Delta#varphi| (rad)"}, {3, -0.5f, 2.5f, "photons of the pair that are physical primaries"}});
    fRegistryPairMC.add("PairClass/Pairs/hTrueRelation", "true distinct pairs (class 13, MC): decay relation of the two photons", kTHnSparseF, {axSample, {60, 0.f, 0.3f, "#it{q}_{inv} (GeV/#it{c})"}, {5, -0.5f, 4.5f, "0 unrelated, 1 mothers are siblings (e.g. #eta #rightarrow 3#pi^{0}, K^{0}_{S} #rightarrow 2#pi^{0}), 2 one photon is a sibling of the other's mother (e.g. #omega #rightarrow #pi^{0}#gamma), 3 same first physical primary otherwise, 4 same root of the decay chain otherwise"}, {7, -0.5f, 6.5f, "common particle: 0 K^{0}_{S}, 1 #eta, 2 #omega, 3 #eta', 4 baryon, 5 other, 6 none"}, {6, 0.f, 3.f, "#it{k}_{T} (GeV/#it{c})"}});
    fRegistryPairMC.add("PairClass/Pairs/hPairR", "photon pairs: conversion radii by MC class", kTHnSparseF, {axSample, axClass, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {90, 0.f, 90.f, "#it{R}_{xy} of the inner photon (cm)"}, {90, 0.f, 90.f, "#it{R}_{xy} of the outer photon (cm)"}, {91, -1.f, 90.f, "true #it{R}_{xy} of the inner photon (cm), -1 without MC"}});
    fRegistryPairMC.add("PairClass/Pairs/hBrems", "bremsstrahlung pairs (classes 16 and 5, MC)", kTHnSparseF, {axSample, {2, -0.5f, 1.5f, "0 class 16 (both true), 1 class 5 (with a fake)"}, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {5, -0.5f, 4.5f, "radiating electron from: 0 photon conversion, 1 #pi^{0} Dalitz, 2 #eta Dalitz, 3 other, 4 not found"}, {6, -0.5f, 5.5f, "mother of the parent photon: 0 #pi^{0}, 1 #eta, 2 #omega, 3 e^{#pm} (bremsstrahlung), 4 other, 5 none / not a photon"}, {50, 0.f, 0.25f, "#theta(#gamma_{brems}, radiating e^{#pm} rec.) (rad)"}, {50, 0.f, 0.25f, "#theta(#gamma_{brems}, parent #gamma) (rad), MC"}, {60, -40.f, 80.f, "true #it{R}_{xy}(#gamma_{brems}) #minus #it{R}_{xy}(parent) (cm)"}, {20, 0.f, 1.f, "#it{E}(#gamma_{brems}) / #it{E}(parent), MC"}});
    fRegistryPairMC.add("PairClass/Pairs/hPhotonQuality", "photons of the pairs by MC class (each photon once per pair)", kTHnSparseF, {axSample, axClass, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {2, -0.5f, 1.5f, "0 photon radiated by an e^{#pm} (MC mother e^{#pm}), 1 other"}, {50, 0.f, 50.f, "KF #chi^{2}/ndf"}, {60, 0.f, 0.3f, "pointing angle to the PV (rad)"}, {50, 0.f, 5.f, "|DCA_{xy}| to the PV (cm)"}, {50, 0.f, 5.f, "|DCA_{z}| to the PV (cm)"}});
    {
      const char* anc = "0 #pi^{0}, 1 #eta, 2 #omega, 3 #eta', 4 #rho/K*, 5 K^{0}, 6 baryon, 7 primary (from the collision), 8 other";
      fRegistryPairMC.add("PairClass/Pairs/hBremsAncestry", "bremsstrahlung pairs (classes 16 and 5, MC): ancestry of the parent photon (the photon whose electron radiated)", kTHnSparseF, {axSample, {2, -0.5f, 1.5f, "0 class 16 (both true), 1 class 5 (with a fake)"}, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {9, -0.5f, 8.5f, Form("mother of the parent photon: %s", anc)}, {9, -0.5f, 8.5f, Form("mother of that particle: %s", anc)}, {9, -0.5f, 8.5f, Form("root of the decay chain: %s", anc)}});
    }
    fRegistryPairMC.add("PairClass/Pairs/hPhotonKF", "photons of the pairs by MC class (each photon once per pair)", kTHnSparseF, {axSample, axClass, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {2, -0.5f, 1.5f, "0 photon radiated by an e^{#pm} (MC mother e^{#pm}), 1 other"}, {50, 0.f, 25.f, "KF #chi^{2}/ndf of the conversion-vertex fit"}, {60, 0.f, 300.f, "KF #chi^{2} to the primary vertex"}, {40, -0.2f, 0.2f, "#psi_{pair} (rad)"}});
    fRegistryPairMC.add("PairClass/Pairs/hHiddenBrems", "true distinct pairs (class 13) with a photon radiated by an e^{#pm} (MC)", kTHnSparseF, {axSample, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {3, -0.5f, 2.5f, "photons of the pair radiated by an e^{#pm}"}, {4, -0.5f, 3.5f, "partner vs the radiating e^{#pm}: 0 unrelated, 1 same mother (e.g. Dalitz #gamma), 2 same root ancestor, 3 partner is itself radiated"}, {6, -0.5f, 5.5f, "radiating e^{#pm} from: 0 photon conversion, 1 #pi^{0} Dalitz, 2 #eta Dalitz, 3 heavy flavour (D, B), 4 other hadron, 5 none / unknown"}, {3, -0.5f, 2.5f, "partner vs the radiating e^{#pm}: 0 different first physical primary, 1 same first physical primary, 2 n/a"}});
    fRegistryPairMC.add("PairClass/Pairs/hPairSep", "photon pairs with q_{inv} < 100 MeV: minimum same-sign leg distance in the TPC by MC class", kTHnSparseF, {axSample, axClass, {20, 0.f, 0.1f, "#it{q}_{inv} (GeV/#it{c})"}, {50, 0.f, 10.f, "min #sqrt{#Delta(r#varphi)^{2} + #Delta z^{2}} of the same-sign legs (cm)"}, {25, 0.f, 5.f, "min #Delta(r#varphi) of the same-sign legs (cm)"}, {10, 0.f, 0.5f, "|#it{E}_{1} #minus #it{E}_{2}| (GeV)"}, {3, -0.5f, 2.5f, "conversion radius R_{xy}: 0 both < 70 cm, 1 one #geq 70 cm, 2 both #geq 70 cm"}, {2, -0.5f, 1.5f, "all four legs below the TPC shared-cluster cut"}, {3, -0.5f, 2.5f, "photons of the pair passing the analysis photon cut (0, 1, 2)"}});
    fRegistryPairMC.add("PairClass/Photons/hLegTrueDist", "stored photons: distance of the true production vertices of the two legs;0 same event, 1 mixed event;0 SVertexer original, 1 finder original, 2 crossed pairing;distance (cm)", kTH3F, {{2, -0.5f, 1.5f}, {3, -0.5f, 2.5f}, {104, -2.f, 50.f}});
    auto hEv = fRegistryPairMC.add<TH1>("PairClass/Events/hArtificialEvents", "artificial events;;events", kTH1F, {{4, -0.5f, 3.5f}});
    hEv->GetXaxis()->SetBinLabel(1, "same event");
    hEv->GetXaxis()->SetBinLabel(2, "mixed event");
    hEv->GetXaxis()->SetBinLabel(3, "standard mixing: partner events used");
    hEv->GetXaxis()->SetBinLabel(4, "standard mixing: events without a partner");
    fRegistryPairMC.add("PairClass/Events/hNPhotons", "photons per artificial event;0 same event, 1 mixed event;photons", kTH2F, {{2, -0.5f, 1.5f}, {100, -0.5f, 99.5f}});
    fRegistryPairMC.add("PairClass/Photons/hOrigin", "stored photons;0 same event, 1 mixed event;0 SVertexer original, 1 finder original, 2 crossed pairing", kTH2F, {{2, -0.5f, 1.5f}, {3, -0.5f, 2.5f}});
  }

  static constexpr float kNoteMaxQ = 0.1f;
  static std::vector<double> noteCentEdges() { return {0., 10., 20., 30., 40., 50., 60., 80., 100.}; }

  static std::vector<double> noteOccEdges()
  {
    std::vector<double> edges;
    for (int i = 0; i <= 30; ++i) { // o2-linter: disable=magic-number (30 bins of 1000)
      edges.push_back(1000. * i);
    }
    edges.insert(edges.end(), {40000., 50000., 1.e6});
    return edges;
  }
  static std::vector<double> noteConvDistEdges() { return {0., 0.5, 1., 2., 3., 5., 7.5, 10., 15., 20., 30., 50., 100., 1000.}; }

  void addNotePairClassHistograms()
  {
    const AxisSpec axSample{4, -0.5f, 3.5f, "0 same event, 1 standard mixed event, 2 co-reconstructed mixed event, 3 co-reco ME, one source event"};
    const AxisSpec axClass{21, -0.5f, 20.5f, "pair class (as PairClass/Pairs/hPairClass)"};
    const AxisSpec axQ{20, 0.f, kNoteMaxQ, "#it{q}_{inv} (GeV/#it{c})"};
    const AxisSpec axQCoarse{std::vector<double>{0., 0.02, 0.05, 0.1}, "#it{q}_{inv} (GeV/#it{c})"};
    const AxisSpec axCent{noteCentEdges(), "centrality FT0C (%)"};
    const AxisSpec axOcc{noteOccEdges(), "occupancy"};
    const AxisSpec axPass{3, -0.5f, 2.5f, "photons of the pair passing the analysis photon cut (0, 1, 2)"};
    const AxisSpec axSep{2, -0.5f, 1.5f, "same-sign leg separation d_{min} > cfgDPairCut"};
    fRegistryPairMC.add("PairClass/Note/hEvents", "artificial events;0 same event, 1 mixed event;centrality FT0C (%);occupancy", kTH3F, {{2, -0.5f, 1.5f}, axCent, axOcc});
    fRegistryPairMC.add("PairClass/Note/hClassCentOcc", "photon pairs by MC class vs centrality and occupancy", kTHnSparseF, {axSample, axClass, axQ, axCent, axOcc, axPass, axSep});
    fRegistryPairMC.add("PairClass/Note/hDupCls", "duplicates (classes 1, 2, 3, 12): TPC clusters of the two pieces of the shared MC particle", kTHnSparseF, {axSample, {4, -0.5f, 3.5f, "class: 0 = 1 (same photon twice), 1 = 2, 2 = 3, 3 = 12 (same track)"}, axQ, {32, 0.f, 160.f, "TPC clusters of the shorter piece"}, {32, 0.f, 160.f, "TPC clusters of the longer piece"}, {20, 0.f, 1.f, "larger TPC shared-cluster fraction of the two pieces"}, {3, -0.5f, 2.5f, "charge slots: 0 same, 1 opposite, 2 no shared particle found"}, axPass, axSep});
    fRegistryPairMC.add("PairClass/Note/hLegCls", "photon pairs by MC class: TPC clusters of the four legs (threshold scan)", kTHnSparseF, {axSample, axClass, axQ, {20, 0.f, 1.f, "largest TPC shared-cluster fraction of the four legs"}, {16, 0.f, 160.f, "smallest TPC cluster count of the four legs"}, axPass, axSep});
    fRegistryPairMC.add("PairClass/Note/hPhotonKine", "photons of the pairs by MC class (each photon once per pair)", kTHnSparseF, {axSample, axClass, axQCoarse, {20, -1.f, 1.f, "#eta"}, {36, 0.f, o2::constants::math::TwoPI, "#varphi (rad)"}, {std::vector<double>{0., 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.5, 2., 3., 10.}, "#it{p}_{T} (GeV/#it{c})"}, axPass, axSep});
    fRegistryPairMC.add("PairClass/Note/hPhotonConv", "photons of the pairs by MC class (each photon once per pair): conversion point and vertex", kTHnSparseF, {axSample, axClass, axQCoarse, {18, 0.f, 90.f, "#it{R}_{conv} (cm)"}, {24, -120.f, 120.f, "#it{z}_{conv} (cm)"}, {10, -10.f, 10.f, "#it{z}_{vtx} (cm)"}, axPass, axSep});
  }

  void addCoRecoHistograms()
  {
    if (qaflags.cfgFillNoteQA.value) {
      addNoteCoRecoHistograms();
    }
    const AxisSpec axSample{2, -0.5f, 1.5f, "0 same event, 1 mixed event"};
    const AxisSpec axOrigin{4, -0.5f, 3.5f, "stored photons: 0 both SVertexer originals, 1 one finder original, 2 both finder, 3 a crossed pairing"};
    const AxisSpec axSign{2, -0.5f, 1.5f, "#alpha_{1}#alpha_{2} > 0 (same charge carries more momentum)"};
    auto hEv = fRegistry.add<TH1>("CoReco/hArtificialEvents", "artificial events (PCMMixCollisions);;events", kTH1F, {{4, -0.5f, 3.5f}});
    hEv->GetXaxis()->SetBinLabel(1, "same event");
    hEv->GetXaxis()->SetBinLabel(2, "same event, centrality");
    hEv->GetXaxis()->SetBinLabel(3, "mixed event");
    hEv->GetXaxis()->SetBinLabel(4, "mixed event, centrality");
    fRegistry.add("CoReco/hNPhotons", "photons per artificial event (before the photon cut);0 same event, 1 mixed event;photons", kTH2F, {axSample, {100, -0.5f, 99.5f}});
    fRegistry.add("CoReco/hOneSourcePairs", "mixed artificial events: pairs with all four legs from one source event;q_{inv} (GeV/c);pairs", kTH1F, {axisQinv});
    fRegistry.add("CoReco/hPairOrigin", "accepted pairs by origin of the stored photons", kTHnSparseF, {axSample, axOrigin, axisQinv, axisKt, axSign});
    fRegistry.add("CoReco/hQPartner", "accepted pairs by partner-event distance (same event: 0)", kTHnSparseF, {axSample, {100, 0.f, 0.5f, "#it{q}_{inv} (GeV/#it{c})"}, {std::vector<double>{0., 0.2, 0.4, 0.6, 1.0, 2.0}, "#it{k}_{T} (GeV/#it{c})"}, {std::vector<double>{0., 0.25, 0.5, 1., 2., 20.}, "|#Delta z_{vtx}| partner (cm)"}, {std::vector<double>{0., 1., 2., 5., 10., 100.}, "|#Delta centrality| partner (%)"}, {std::vector<double>{0., 250., 500., 1000., 2000., 1.e6}, "|#Delta occupancy| partner"}, {std::vector<double>{0., 1., 10., 60., 600., 1.e7}, "|#Delta t| partner (s)"}});
    if (!isMC) {
      return;
    }
    const AxisSpec axClass{21, -0.5f, 20.5f, "pair class: 13 true distinct, 14 full swap, 4 half swap, 1-3 duplicates, 16 bremsstrahlung (both true), 5 bremsstrahlung with a fake, 18 photons of the same mother (pi0/eta, both true), 19 swap of two such photons, 20 same mother with a fake, 17 related swap, 15 no label, 6-12 other fakes"};
    const AxisSpec axQFine{300, 0.f, 0.3f, "#it{q}_{inv} (GeV/#it{c})"};
    const AxisSpec axQTrue{301, -0.001f, 0.3f, "#it{q}_{inv} of the two MC photons (GeV/#it{c}), -1 if not two photons"};
    const AxisSpec axKtCoarse{std::vector<double>{0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 2.0}, "#it{k}_{T} (GeV/#it{c})"};
    fRegistryPairMC.add("CoReco/MC/hClassQ", "accepted pairs by MC class: reconstructed vs true q", kTHnSparseF, {axSample, axClass, axQFine, axQTrue, axKtCoarse});
    AxisSpec axKtCls = axisKt; // edges of confKtBins, the title the readers of the class maps expect
    axKtCls.title = "#it{k}_{T} (GeV/#it{c})";
    fRegistryPairMC.add("CoReco/MC/hClassQinvQout", "accepted pairs by MC class: q_{inv} vs |q_{out}| (LCMS); standard mixed event: class 13 both true, class 0 with a fake", kTHnSparseF, {axSample, axClass, axisQinv, axisQout, axKtCls});
  }

  void addNoteCoRecoHistograms()
  {
    const AxisSpec axSample{2, -0.5f, 1.5f, "0 same event, 1 mixed event"};
    const AxisSpec axCent{noteCentEdges(), "centrality (%)"};
    const AxisSpec axOcc{noteOccEdges(), "occupancy"};
    fRegistry.add("CoReco/Note/hEvents", "artificial events in the centrality window;0 same event, 1 mixed event;centrality (%);occupancy", kTH3F, {axSample, axCent, axOcc});
    if (!isMC) {
      return;
    }
    const AxisSpec axClass{21, -0.5f, 20.5f, "pair class (as PairClass/Pairs/hPairClass)"};
    const AxisSpec axStage{2, -0.5f, 1.5f, "0 two selected photons, before the pair cuts, 1 after all pair cuts"};
    const AxisSpec axKtCoarse{std::vector<double>{0., 0.2, 0.4, 0.6, 1.0, 2.0}, "#it{k}_{T} (GeV/#it{c})"};
    fRegistryPairMC.add("CoReco/MC/Note/hClassStage", "pairs by MC class vs centrality and occupancy, before / after the pair cuts", kTHnSparseF, {axSample, axClass, {60, 0.f, 0.3f, "#it{q}_{inv} (GeV/#it{c})"}, axCent, axOcc, axStage});
    fRegistryPairMC.add("CoReco/MC/Note/hClassQinvQout", "pairs by MC class: q_{inv} vs |q_{out}| (LCMS)", kTHnSparseF, {axSample, axClass, axisQinv, axisQout, axKtCoarse, axStage});
    fRegistryPairMC.add("CoReco/MC/Note/hClassQinvQside", "pairs by MC class: q_{inv} vs |q_{side}| (LCMS)", kTHnSparseF, {axSample, axClass, axisQinv, axisQside, axKtCoarse, axStage});
    fRegistryPairMC.add("CoReco/MC/Note/hClassQinvQlong", "pairs by MC class: q_{inv} vs |q_{long}| (LCMS)", kTHnSparseF, {axSample, axClass, axisQinv, axisQlong, axKtCoarse, axStage});
    const AxisSpec axRec{150, 0.f, 0.3f, "reconstructed (GeV/#it{c})"};
    const AxisSpec axTrue{151, -0.002f, 0.3f, "of the two MC photons (GeV/#it{c}), < 0 if not two photons"};
    fRegistryPairMC.add("CoReco/MC/Note/hRespQout", "response |q_{out}| (LCMS) by MC class", kTHnSparseF, {axSample, axClass, axRec, axTrue, axKtCoarse, axStage});
    fRegistryPairMC.add("CoReco/MC/Note/hRespQside", "response |q_{side}| (LCMS) by MC class", kTHnSparseF, {axSample, axClass, axRec, axTrue, axKtCoarse, axStage});
    fRegistryPairMC.add("CoReco/MC/Note/hRespQlong", "response |q_{long}| (LCMS) by MC class", kTHnSparseF, {axSample, axClass, axRec, axTrue, axKtCoarse, axStage});
    fRegistryPairMC.add("CoReco/MC/Note/hRespDE", "response |#Delta E| by MC class", kTHnSparseF, {axSample, axClass, {100, 0.f, 0.5f, "|#Delta#it{E}| reconstructed (GeV)"}, {101, -0.005f, 0.5f, "|#Delta#it{E}| of the two MC photons (GeV), < 0 if not two photons"}, {std::vector<double>{0., 0.01, 0.02, 0.03, 0.05, 0.1, 0.3}, "#it{q}_{inv} (GeV/#it{c})"}, axStage});
    fRegistryPairMC.add("CoReco/MC/Note/hTruePairs", "reconstructed pairs of two MC photons vs the true pair variables", kTHnSparseF, {axSample, {3, -0.5f, 2.5f, "0 true distinct (class 13), 1 full swap (class 14), 2 other class"}, {60, 0.f, 0.3f, "#it{q}_{inv}^{true} (GeV/#it{c})"}, {std::vector<double>{0., 0.2, 0.4, 0.6, 1.0, 2.0}, "#it{k}_{T}^{true} (GeV/#it{c})"}, axCent, {noteConvDistEdges(), "distance of the two true conversion points (cm)"}, axStage});
  }

  static std::vector<double> pipeOccEdges() { return {0., 1000., 3000., 7000., 12000., 20000., 35000., 1.e6}; } // as PairClass/Pairs/hTrueOrigin
  void addPipeHistograms()
  {
    const AxisSpec axSample{2, -0.5f, 1.5f, "0 same event, 1 mixed event"};
    const AxisSpec axTag{5, -0.5f, 4.5f, "pointing d_{min}: 0 < 0.1, 1 0.1-0.3, 2 0.3-1, 3 > 1 cm, 4 no allowed point"};
    AxisSpec axQoutAbs = axisQout; // binning of CF_2D
    axQoutAbs.title = "|#it{q}_{out}| (LCMS) (GeV/#it{c})";
    const AxisSpec axOcc{pipeOccEdges(), "occupancy"};
    const AxisSpec axQCoarse{std::vector<double>{0., 0.02, 0.15, 0.3}, "#it{q}_{inv} (GeV/#it{c})"};
    fRegistry.add("CoReco/Pipe/hTagMap", "accepted pairs with q_{inv} < cfgPipeTagMaxQinv: bremsstrahlung pointing tag", kTHnSparseF, {axSample, axTag, axQoutAbs, axisQinv, axisKt});
    fRegistry.add("CoReco/Pipe/hOriginMap", "accepted pairs by origin of the stored photons and partner relation", kTHnSparseF, {axSample, {4, -0.5f, 3.5f, "stored photons: 0 both SVertexer originals, 1 one finder original, 2 both finder, 3 a crossed pairing"}, {4, -0.5f, 3.5f, "partner relation: 0 both run-1 pairings, 1 shared run-1 source photon, 2 free legs of one source event with opposite charge, 3 unrelated"}, axQoutAbs, axisQinv, axisKt});
    fRegistry.add("CoReco/Pipe/hPhotonDcaz", "photons of the accepted pairs (each photon once per pair)", kTHnSparseF, {axSample, {100, -5.f, 5.f, "DCA_{z} of the photon to the PV (cm)"}, {3, -0.5f, 2.5f, "legs with ITS hits"}, axOcc, {std::vector<double>{0., 0.2, 0.4, 0.6, 1., 2., 5.}, "#it{p}_{T} (GeV/#it{c})"}});
    fRegistry.add("CoReco/Pipe/hDEtaDPhiOcc", "accepted pairs: #Delta#eta vs #Delta#varphi vs occupancy", kTHnSparseF, {axSample, {64, -1.6f, 1.6f, "#Delta#eta"}, {36, -o2::constants::math::PI, o2::constants::math::PI, "#Delta#varphi (rad)"}, axOcc, axQCoarse});
    fRegistry.add("CoReco/Pipe/hEta1Eta2Occ", "accepted pairs: #eta of the two photons vs occupancy", kTHnSparseF, {axSample, {32, -0.8f, 0.8f, "#eta of the lower-#it{p}_{T} photon"}, {32, -0.8f, 0.8f, "#eta of the higher-#it{p}_{T} photon"}, axOcc, axQCoarse});
    fRegistry.add("CoReco/Pipe/hPtKt", "accepted pairs: #it{p}_{T} of the two photons vs #it{k}_{T}", kTHnSparseF, {axSample, {40, 0.f, 2.f, "#it{p}_{T} of the lower-#it{p}_{T} photon (GeV/#it{c})"}, {40, 0.f, 2.f, "#it{p}_{T} of the higher-#it{p}_{T} photon (GeV/#it{c})"}, axisKt, axQCoarse});
    if (!isMC) {
      return;
    }
    const AxisSpec axClass{21, -0.5f, 20.5f, "pair class (as PairClass/Pairs/hPairClass)"};
    AxisSpec axKtCls = axisKt; // as the k_T axis of CoReco/MC/hClassQinvQout
    axKtCls.title = "#it{k}_{T} (GeV/#it{c})";
    AxisSpec axQoutRec = axisQout, axQinvRec = axisQinv; // binning of CoReco/MC/hClassQinvQout
    axQoutRec.title = "|#it{q}_{out}| reconstructed (GeV/#it{c})";
    axQinvRec.title = "#it{q}_{inv} reconstructed (GeV/#it{c})";
    fRegistryPairMC.add("CoReco/MC/Pipe/hClassTagMap", "accepted pairs with q_{inv} < cfgPipeTagMaxQinv by MC class: bremsstrahlung pointing tag", kTHnSparseF, {axSample, axClass, axTag, axQoutAbs, axisQinv, axKtCls});
    fRegistryPairMC.add("CoReco/MC/Pipe/hResp", "accepted true distinct pairs (class 13), same event: joint response of |q_{out}| (LCMS) and q_{inv}", kTHnSparseF, {axQoutRec, axQinvRec, {60, 0.f, 0.3f, "|#it{q}_{out}| of the two MC photons (GeV/#it{c})"}, {60, 0.f, 0.3f, "#it{q}_{inv} of the two MC photons (GeV/#it{c})"}, {std::vector<double>{0.2, 0.4, 0.6, 1., 2.}, "#it{k}_{T} reconstructed (GeV/#it{c})"}});
  }

  /// centrality of the estimator cfgCentEstimator
  template <typename TCollision>
  float eventCentrality(TCollision const& collision) const
  {
    const std::array<float, 3> cent = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
    return cent[mixing.cfgCentEstimator.value];
  }

  [[nodiscard]] bool isSelectedCentrality(const float cent) const
  {
    return !(cent < centralitySelection.cfgCentMin.value || centralitySelection.cfgCentMax.value < cent); // NOLINT(readability-simplify-boolean-expr) a NaN centrality passes, as the old continue did
  }

  bool isSelectedPhoton(MyV0PhotonsMix const& photons, const int64_t id)
  {
    auto& pass = mPhotonPass[id];
    if (pass < 0) {
      const auto g = photons.rawIteratorAt(id);
      pass = fV0PhotonCut.IsSelected<decltype(g), aod::V0LegsMix>(g) ? 1 : 0;
    }
    return pass == 1;
  }

  inline bool isInsideEllipse(float deta, float dphi) const
  {
    if (!ggpaircuts.cfgDoEllipseCut.value) {
      return false;
    }
    const float sE = ggpaircuts.cfgEllipseSigEta.value;
    const float sP = ggpaircuts.cfgEllipseSigPhi.value;
    if (sE < kMinSigma || sP < kMinSigma) {
      return false;
    }
    return (deta / sE) * (deta / sE) + (dphi / sP) * (dphi / sP) < ggpaircuts.cfgEllipseR2.value;
  }

  inline bool passRZCut(float deltaR, float deltaZ) const
  {
    if (ggpaircuts.cfgDoRCut.value && deltaR < ggpaircuts.cfgMinDeltaR.value) {
      return false;
    }
    if (ggpaircuts.cfgDoZCut.value && std::fabs(deltaZ) < ggpaircuts.cfgMinDeltaZ.value) {
      return false;
    }
    return true;
  }

  [[nodiscard]] inline bool passPhotonClassPairCut(pairutil::V0PhotonLegCounts const& c1,
                                                   pairutil::V0PhotonLegCounts const& c2) const
  {
    if (!cfgDoPhotonClassPairCut.value) {
      return true;
    }
    return pairutil::isPairPhotonClassSelected(c1, c2, mPhotonClassSelA, mPhotonClassSelB);
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

  [[nodiscard]] bool fillsPairSep(float qinv) const
  {
    const float limit = pairsep.cfgPairSepMaxQinv.value;
    return mQa.pairSepQa && !(limit > 0.f && qinv >= limit); // NOLINT(readability-simplify-boolean-expr) exact negation of the early return in fillPairSep, also for NaN
  }

  /// computePairSep is needed for this pair: its only readers are the PairSep QA and the pair-merging cut
  [[nodiscard]] bool needsPairSep(float qinv) const
  {
    return fillsPairSep(qinv) || pairsep.cfgDoPairMergeCut.value;
  }

  [[nodiscard]] PairSepConfig pairSepConfig() const
  {
    PairSepConfig c;
    c.sigR = std::max(pairsep.cfgSigRPhi.value, kMinSigma);
    const float sigZcm = (pairsep.cfgSigZNs.value > 0.f && mVDriftCmPerNs > 0.f)
                           ? pairsep.cfgSigZNs.value * mVDriftCmPerNs
                           : pairsep.cfgSigZ.value;
    c.sigZ = std::max(sigZcm, kMinSigma);
    c.dClose = pairsep.cfgDPairClose.value;
    for (std::size_t t = 0; t < c.closeThr.size(); ++t) {
      c.closeThr[t] = pairsep.cfgCloseThrCm.value[t];
    }
    return c;
  }

  template <typename TG1, typename TG2>
  [[nodiscard]] inline bool passFastQinvGate(TG1 const& g1, TG2 const& g2) const
  {
    const float qMax = qaflags.cfgMaxQinvForProcessing.value;
    if (qMax > 1e9f) { // o2-linter: disable=magic-number (skip if non-sensical value is chosen)
      return true;
    }
    const float dEta = g1.eta() - g2.eta();
    const float dPhi = RecoDecay::constrainAngle(g1.phi() - g2.phi(), -o2::constants::math::PI);
    const float q2 = 2.f * g1.pt() * g2.pt() * (std::cosh(dEta) - std::cos(dPhi));
    return q2 <= qMax * qMax;
  }

  [[nodiscard]] inline bool passPairMergeCut(PairSep const& s) const
  {
    if (!pairsep.cfgDoPairMergeCut.value) {
      return true;
    }
    if (s.nCommon <= 0) { // nothing was propagated: no statement possible
      return true;
    }
    if (s.dPairMin < pairsep.cfgDPairCut.value) {
      return false;
    }
    const int nCut = pairsep.cfgNCloseCut.value;
    // dipping close once is survivable, running alongside is not
    return nCut <= 0 || s.nCloseScaled < nCut;
  }

  [[nodiscard]] inline bool passCrossPairVeto(CrossObs const& c) const
  {
    if (crosspair.cfgDoCrossPairCut.value && c.meeOverQ <= crosspair.cfgCrossMaxMeeRatio.value) {
      return false;
    }
    if (crosspair.cfgDoCrossGeomCut.value && c.geomPhotonLike) {
      return false;
    }
    return true;
  }

  [[nodiscard]] bool passPairCuts(PhotonWithLegs const& g1, PhotonWithLegs const& g2, PairQAObservables const& obs, PairSep const& sep) const
  {
    if (!passPhotonClassPairCut(g1.legCounts(), g2.legCounts())) {
      return false;
    }
    if (obs.drOverCosOA < ggpaircuts.cfgMinDRCosOA.value) {
      return false;
    }
    if (!passPairMergeCut(sep)) {
      return false;
    }
    if (!passRZCut(obs.deltaR, obs.deltaZ)) {
      return false;
    }
    if (isInsideEllipse(obs.deta, obs.dphi)) {
      return false;
    }
    if (crosspair.cfgDoCrossPairCut.value || crosspair.cfgDoCrossGeomCut.value) {
      const auto cross = computeCrossObs(g1, g2, obs.qinv, mBzT, crosspair.cfgCrossMaxGap.value, crosspair.cfgCrossMaxDTgl.value, crosspair.cfgCrossMaxDz.value);
      if (!passCrossPairVeto(cross)) {
        return false;
      }
    }
    return true;
  }

  static void parseBins(const ConfigurableAxis& cfg, std::vector<float>& edges)
  {
    if (cfg.value[0] == VARIABLE_WIDTH) {
      edges = std::vector<float>(cfg.value.begin(), cfg.value.end());
      edges.erase(edges.begin());
    } else {
      const auto n = static_cast<int>(cfg.value[0]);
      const auto xmin = static_cast<float>(cfg.value[1]);
      const auto xmax = static_cast<float>(cfg.value[2]);
      edges.resize(n + 1);
      for (int i = 0; i <= n; ++i) {
        edges[i] = xmin + (xmax - xmin) / n * i;
      }
    }
  }

  static int clampBin(int b, int nmax) { return std::clamp(b, 0, nmax); }

  static int binOf(const std::vector<float>& edges, float val)
  {
    const auto pos = static_cast<int>(std::lower_bound(edges.begin(), edges.end(), val) - edges.begin());
    return clampBin(pos - 1, static_cast<int>(edges.size()) - 2);
  }

  template <int ev_id>
  static constexpr const char* pairSepPrefix()
  {
    return (ev_id == 0) ? "PairSep/same/" : "PairSep/mix/";
  }

  template <int ev_id>
  static constexpr const char* fullRangePrefix()
  {
    return (ev_id == 0) ? "Pair/same/FullRange/" : "Pair/mix/FullRange/";
  }

  template <int ev_id>
  inline void fillFullRangeQA(PairQAObservables const& obs, float occupancy)
  {
    constexpr auto base = fullRangePrefix<ev_id>();
    fRegistryCF.fill(HIST(base) + HIST("hQinvVsOccupancy"), occupancy, obs.qinv);
  }

  template <int ev_id>
  inline void fillPairSep(PairSep const& s, PairQAObservables const& obs)
  {
    if (!mQa.pairSepQa) {
      return;
    }
    const float limit = pairsep.cfgPairSepMaxQinv.value;
    if (limit > 0.f && obs.qinv >= limit) {
      return;
    }
    constexpr auto dir = pairSepPrefix<ev_id>();
    if (s.nCommon > 0) {
      for (std::size_t t = 0; t < s.fCloseRPhiThr.size(); ++t) {
        fRegistryCF.fill(HIST(dir) + HIST("hSparse_dRPhiMin_fClose_Thr_Qinv"), s.dMinRPhi,
                         s.fCloseRPhiThr[t], static_cast<float>(t), obs.qinv);
      }
      fRegistryCF.fill(HIST(dir) + HIST("hSparse_DPair_NClose_Qinv"), s.dPairMin, static_cast<float>(s.nCloseScaled), obs.qinv);
    }
  }

  template <typename TPhoton>
  inline void fillSinglePhotonQA(TPhoton const& g)
  {
    if (!mQa.singlePhotonQa) {
      return;
    }
    const float r = std::sqrt(g.vx() * g.vx() + g.vy() * g.vy());
    fRegistryPairQA.fill(HIST("SinglePhoton/Before/hEtaVsPhiPt"), g.phi(), g.eta(), g.pt());
    fRegistryPairQA.fill(HIST("SinglePhoton/Before/hRVsZConvPt"), g.vz(), r, g.pt());
  }

  inline void fillPairQA(PairQAObservables const& o)
  {
    if (!mQa.pairQa) {
      return;
    }
    //// Delta R (Conversion point) QA
    fRegistryPairQA.fill(HIST("Pair/same/QA/Before/hR1VsR2"), o.r1, o.r2);
  }

  template <int ev_id, typename TCollision>
  void fillPairHistogram(TCollision const& collision,
                         ROOT::Math::PtEtaPhiMVector const& v1,
                         ROOT::Math::PtEtaPhiMVector const& v2,
                         float weight = 1.f)
  {
    const float multNTracks = collision.multNTracksPVeta1();
    const float multFT0M = collision.multFT0M();
    const auto p12 = v1 + v2;
    auto k12 = 0.5 * p12;
    float kt = k12.Pt();
    float qinv = -((v1 - v2).M()); // M() < 0 for the space-like difference of two photons
    ROOT::Math::PxPyPzEVector v1c(v1), v2c(v2);
    float beta_z = p12.Beta() * std::cos(p12.Theta());
    ROOT::Math::Boost bst_z(0, 0, -beta_z);
    auto q3_lcms = bst_z(v1c - v2c).Vect();
    float qabs_lcms = q3_lcms.R();
    if (hbtanalysis.cfgDo3D) {
      if (hbtanalysis.cfgDo2D) {
        ROOT::Math::XYZVector uv_out(k12.Px() / k12.Pt(), k12.Py() / k12.Pt(), 0);
        float qout_lcms = q3_lcms.Dot(uv_out);
        if constexpr (ev_id == 0) {
          fRegistryCF.fill(HIST("Pair/same/CF_2D"), std::fabs(qout_lcms), std::fabs(qinv), multNTracks, multFT0M, kt, weight);
        } else {
          fRegistryCF.fill(HIST("Pair/mix/CF_2D"), std::fabs(qout_lcms), std::fabs(qinv), multNTracks, multFT0M, kt, weight);
        }
      }
    } else {
      if constexpr (ev_id == 0) {
        fRegistryCF.fill(HIST("Pair/same/CF_1D"), hbtanalysis.cfgUseLCMS ? qabs_lcms : qinv, kt, weight);
      } else {
        fRegistryCF.fill(HIST("Pair/mix/CF_1D"), hbtanalysis.cfgUseLCMS ? qabs_lcms : qinv, kt, weight);
      }
    }
    if constexpr (ev_id == 0) {
      fRegistryCF.fill(HIST("Pair/same/CF_QLcms_Qinv"), qabs_lcms, std::fabs(qinv), kt, weight);
    } else {
      fRegistryCF.fill(HIST("Pair/mix/CF_QLcms_Qinv"), qabs_lcms, std::fabs(qinv), kt, weight);
    }
  }

  void fillCoRecoNote(const float sample, const float stage, PcPhoton const& a, PcPhoton const& b, PairQAObservables const& obs, const float cent, const float occupancy)
  {
    if (obs.qinv >= qaflags.cfgNoteMaxQinv.value) {
      return;
    }
    const int clsId = pcPairClass(a, b);
    const auto cls = static_cast<float>(clsId);
    constexpr float kMaxQ = 0.3f, kMaxDE = 0.5f;
    const auto q3 = lcmsQ(obs.v1, obs.v2);
    const float dERec = std::fabs(static_cast<float>(obs.v1.E() - obs.v2.E()));
    if (obs.qinv < kMaxQ) {
      fRegistryPairMC.fill(HIST("CoReco/MC/Note/hClassStage"), sample, cls, obs.qinv, cent, occupancy, stage);
    }
    fRegistryPairMC.fill(HIST("CoReco/MC/Note/hClassQinvQout"), sample, cls, obs.qinv, q3[0], obs.kt, stage);
    fRegistryPairMC.fill(HIST("CoReco/MC/Note/hClassQinvQside"), sample, cls, obs.qinv, q3[1], obs.kt, stage);
    fRegistryPairMC.fill(HIST("CoReco/MC/Note/hClassQinvQlong"), sample, cls, obs.qinv, q3[2], obs.kt, stage);
    std::array<float, 3> q3True = {-0.001f, -0.001f, -0.001f};
    float dETrue = -0.001f;
    const auto mothers = pcTrueMothers(a, b);
    if (mothers[0] != nullptr) {
      const ROOT::Math::PxPyPzMVector m1(mothers[0]->mpx, mothers[0]->mpy, mothers[0]->mpz, 0.f), m2(mothers[1]->mpx, mothers[1]->mpy, mothers[1]->mpz, 0.f);
      q3True = lcmsQ(ROOT::Math::PtEtaPhiMVector(m1), ROOT::Math::PtEtaPhiMVector(m2));
      dETrue = std::fabs(static_cast<float>(m1.E() - m2.E()));
      const float qTrue = o2::pwgem::photonmeson::qInv(mothers[0]->mpx, mothers[0]->mpy, mothers[0]->mpz, mothers[1]->mpx, mothers[1]->mpy, mothers[1]->mpz);
      // fiducial true photons as in runTruthEfficiency: physical primary (the leg's first physical-primary ancestor is its mother), |eta|, pT
      const float ptMin = mctruth.cfgMCMinV0Pt.value > 0.f ? mctruth.cfgMCMinV0Pt.value : pcmcuts.cfgMinPtV0.value;
      auto fiducial = [&](PcLeg const& l, ROOT::Math::PxPyPzMVector const& m) { return l.anc >= 0 && l.anc == l.mid && std::fabs(m.Eta()) <= pcmcuts.cfgMaxEtaV0.value && m.Pt() >= ptMin; };
      if (qTrue < kMaxQ && fiducial(*mothers[0], m1) && fiducial(*mothers[1], m2)) {
        const float ktTrue = 0.5f * static_cast<float>((m1 + m2).Pt());
        const float dConv = std::hypot(mothers[0]->vx - mothers[1]->vx, mothers[0]->vy - mothers[1]->vy, mothers[0]->vz - mothers[1]->vz);
        const float kind = clsId == kClsTrueDistinct ? 0.f : (clsId == kClsFullSwap ? 1.f : 2.f);
        fRegistryPairMC.fill(HIST("CoReco/MC/Note/hTruePairs"), sample, kind, qTrue, std::min(ktTrue, 1.99f), cent, std::min(dConv, 999.f), stage);
      }
    }
    if (q3[0] < kMaxQ) {
      fRegistryPairMC.fill(HIST("CoReco/MC/Note/hRespQout"), sample, cls, q3[0], std::min(q3True[0], 0.2999f), obs.kt, stage);
    }
    if (q3[1] < kMaxQ) {
      fRegistryPairMC.fill(HIST("CoReco/MC/Note/hRespQside"), sample, cls, q3[1], std::min(q3True[1], 0.2999f), obs.kt, stage);
    }
    if (q3[2] < kMaxQ) {
      fRegistryPairMC.fill(HIST("CoReco/MC/Note/hRespQlong"), sample, cls, q3[2], std::min(q3True[2], 0.2999f), obs.kt, stage);
    }
    if (dERec < kMaxDE && obs.qinv < kMaxQ) {
      fRegistryPairMC.fill(HIST("CoReco/MC/Note/hRespDE"), sample, cls, dERec, std::min(dETrue, 0.4999f), obs.qinv, stage);
    }
  }

  void fillPipePair(const float sample, PcPhoton const* a, PcPhoton const* b, PairQAObservables const& obs, const float qOut, const int clsId, const float origin, const float occupancy, aod::PCMMixLegs const* mixLegs)
  {
    // only pairs with Q_inv inside the Q axis: the Q_inv overflow is read by nothing
    if (!(obs.qinv < mPipeQinvMax)) {
      return;
    }
    const float occ = std::max(occupancy, 0.f); // as PairClass/Pairs/hTrueOrigin: a negative occupancy (none) goes to the first bin
    // ─── kinematic maps: every pair, also the standard mixed event ───
    const bool firstSofter = obs.v1.Pt() <= obs.v2.Pt(); // equal p_T: the first photon counts as the lower-p_T one
    const auto& soft = firstSofter ? obs.v1 : obs.v2;
    const auto& hard = firstSofter ? obs.v2 : obs.v1;
    fRegistry.fill(HIST("CoReco/Pipe/hDEtaDPhiOcc"), sample, obs.deta, obs.dphi, occ, obs.qinv);
    fRegistry.fill(HIST("CoReco/Pipe/hEta1Eta2Occ"), sample, soft.Eta(), hard.Eta(), occ, obs.qinv);
    fRegistry.fill(HIST("CoReco/Pipe/hPtKt"), sample, soft.Pt(), hard.Pt(), obs.kt, obs.qinv);
    // standard mixed event: the pooled photons carry no leg records, so no origin, DCA_z or tag
    if (a == nullptr || b == nullptr) {
      return;
    }
    // ─── origin and partner relation, photon DCA_z ───
    fRegistry.fill(HIST("CoReco/Pipe/hOriginMap"), sample, origin, static_cast<float>(pcPartner(*a, *b)), qOut, obs.qinv, obs.kt);
    const std::array<std::pair<PcPhoton const*, float>, 2> photonsPt{{{a, static_cast<float>(obs.v1.Pt())}, {b, static_cast<float>(obs.v2.Pt())}}};
    for (const auto& [p, pt] : photonsPt) {
      fRegistry.fill(HIST("CoReco/Pipe/hPhotonDcaz"), sample, p->dcaz, static_cast<float>(static_cast<int>(p->pos.its) + static_cast<int>(p->ele.its)), occ, pt);
    }
    // ─── bremsstrahlung pointing tag: only below cfgPipeTagMaxQinv (helix scan of the four legs) ───
    const bool tagged = obs.qinv < hbtanalysis.cfgPipeTagMaxQinv.value;
    float tag = 0.f;
    if (tagged) {
      tag = pcPointingCategory(pcBremsPointingDMin(*a, *b, mixLegs, mBzT));
      fRegistry.fill(HIST("CoReco/Pipe/hTagMap"), sample, tag, qOut, obs.qinv, obs.kt);
    }
    if (isMC) {
      fillPipePairMC(sample, *a, *b, obs, qOut, clsId, tagged, tag);
    }
  }

  /// MC part of fillPipePair: tag map by class (tagged pairs), joint response of the true distinct pairs of the same event
  void fillPipePairMC(const float sample, PcPhoton const& a, PcPhoton const& b, PairQAObservables const& obs, const float qOut, const int clsId, const bool tagged, const float tag)
  {
    if (tagged) {
      fRegistryPairMC.fill(HIST("CoReco/MC/Pipe/hClassTagMap"), sample, static_cast<float>(clsId), tag, qOut, obs.qinv, obs.kt);
    }
    if (sample > 0.5f || clsId != kClsTrueDistinct) { // o2-linter: disable=magic-number (same event only)
      return;
    }
    // true |q_out| and q_inv of the two MC photons, as in fillCoRecoNote; class 13 always has exactly two
    const auto mothers = pcTrueMothers(a, b);
    if (mothers[0] == nullptr) {
      return;
    }
    const ROOT::Math::PxPyPzMVector m1(mothers[0]->mpx, mothers[0]->mpy, mothers[0]->mpz, 0.f), m2(mothers[1]->mpx, mothers[1]->mpy, mothers[1]->mpz, 0.f);
    const auto q3True = lcmsQ(ROOT::Math::PtEtaPhiMVector(m1), ROOT::Math::PtEtaPhiMVector(m2));
    const float qTrue = o2::pwgem::photonmeson::qInv(mothers[0]->mpx, mothers[0]->mpy, mothers[0]->mpz, mothers[1]->mpx, mothers[1]->mpy, mothers[1]->mpz);
    fRegistryPairMC.fill(HIST("CoReco/MC/Pipe/hResp"), qOut, obs.qinv, q3True[0], qTrue, obs.kt);
  }

  template <int ev_id, typename TCollision>
  int pairCoRecoEvent(TCollision const& collision, std::vector<int64_t> const& ids, MyV0PhotonsMix const& photons,
                      std::vector<PcPhoton> const& pc, aod::PCMMixLegs const& mixLegs, const float centForQA, const float occupancy, const bool pool)
  {
    constexpr bool kSame = ev_id == 0;
    constexpr float kSample = kSame ? 0.f : 1.f;
    // ─── photon selection ─────────────────────────────────────
    std::vector<bool> selected(ids.size(), false);
    for (size_t i = 0; i < ids.size(); ++i) {
      selected[i] = isSelectedPhoton(photons, ids[i]);
      if (kSame && selected[i] && mQa.singlePhotonQa) {
        fillSinglePhotonQA(photons.rawIteratorAt(ids[i]));
      }
    }
    const auto keyDFCollision = std::make_pair(ndf, static_cast<int>(collision.globalIndex()));
    const PairSepConfig sepConfig = pairSepConfig();
    // leg-propagated photons, built at most once per photon of this event and only for photons that reach a pair
    std::vector<PhotonWithLegs> withLegs(ids.size());
    std::vector<bool> withLegsBuilt(ids.size(), false);
    auto photonWithLegs = [&](const size_t k) -> PhotonWithLegs const& {
      if (!withLegsBuilt[k]) {
        const auto g = photons.rawIteratorAt(ids[k]);
        const auto pos = g.template posTrack_as<aod::V0LegsMix>(), ele = g.template negTrack_as<aod::V0LegsMix>();
        withLegs[k] = makePhotonWithLegs(g, pos, ele, collision.posZ(), mBzT, mQa.legPropagation);
        if (isMC) {
          withLegs[k].fIsTruePhoton = pcHasLabels(pc[ids[k]]) ? (pcTruePhoton(pc[ids[k]]) ? 1 : 0) : -1;
        }
        withLegsBuilt[k] = true;
      }
      return withLegs[k];
    };
    int ndiphoton = 0;
    for (size_t i = 0; i < ids.size(); ++i) {
      if (!selected[i]) {
        continue;
      }
      for (size_t j = i + 1; j < ids.size(); ++j) {
        if (!selected[j]) {
          continue;
        }
        const auto g1 = photons.rawIteratorAt(ids[i]);
        const auto g2 = photons.rawIteratorAt(ids[j]);
        const auto& a = pc[ids[i]];
        const auto& b = pc[ids[j]];
        // ─── pair cuts before the pair observables ───────────────
        if (pcSharedTrack(a, b)) {
          continue;
        }
        if constexpr (!kSame) {
          if (a.pos.evt == a.ele.evt && a.pos.evt == b.pos.evt && a.pos.evt == b.ele.evt) {
            fRegistry.fill(HIST("CoReco/hOneSourcePairs"), o2::pwgem::photonmeson::qInv(a.px, a.py, a.pz, b.px, b.py, b.pz));
            if (mixing.cfgCoRecoSkipOneSource.value) {
              continue;
            }
          }
        }
        if (!passAsymmetryCut(g1.pt(), g2.pt())) {
          continue;
        }
        auto obs = buildPairQAObservables(g1, g2);
        if (!obs.valid) {
          continue;
        }
        const bool doQA = passQinvQAGate(obs.qinv), doFR = passQinvFullRangeGate(obs.qinv);
        const auto& pwl1 = photonWithLegs(i);
        const auto& pwl2 = photonWithLegs(j);
        if (isMC && qaflags.cfgFillNoteQA.value) {
          fillCoRecoNote(kSample, 0.f, a, b, obs, centForQA, occupancy);
        }
        PairSep sep;
        if (needsPairSep(obs.qinv)) {
          sep = computePairSep(pwl1, pwl2, sepConfig);
        }

        // ─── QA before the pair cuts ──────────────────────────────
        fillPairSep<ev_id>(sep, obs);
        if (doQA) {
          if constexpr (kSame) {
            fillPairQA(obs);
          }
        }

        // ─── pair cuts ────────────────────────────────────────────
        if (!passPairCuts(pwl1, pwl2, obs, sep)) {
          continue;
        }

        // ─── partner pairing of the run-2 crossed photons (cfgCoRecoPartnerPairing) ───
        if (const int partner = pcPartner(a, b); partner > 0) {
          const bool kept = pcPartnerKept(partner, mixing.cfgCoRecoPartnerPairing.value);
          if (!kept) {
            continue;
          }
        }

        // ─── fills after all pair cuts ────────────────────────────
        if (doFR) {
          fillFullRangeQA<ev_id>(obs, occupancy);
        }
        fillPairHistogram<ev_id>(collision, obs.v1, obs.v2, 1.f);
        // |q_out| (LCMS) and MC class of the pair, computed once for the MC maps and the pipeline histograms
        const bool pipePair = hbtanalysis.cfgFillPipe.value && obs.qinv < mPipeQinvMax; // pipeline maps: Q_inv inside the Q axis
        const auto q3 = (isMC || pipePair) ? lcmsQ(obs.v1, obs.v2) : std::array<float, 3>{};
        const int clsId = isMC ? pcPairClass(a, b) : 0;
        const auto origin = static_cast<float>(pcPairOrigin(a, b)); // for CoReco/hPairOrigin and the pipeline maps
        if (pipePair) {
          fillPipePair(kSample, &a, &b, obs, q3[0], clsId, origin, occupancy, &mixLegs);
        }
        ndiphoton++;
        fRegistry.fill(HIST("CoReco/hQPartner"), kSample, obs.qinv, obs.kt, std::fabs(collision.partnerDz()), std::fabs(collision.partnerDCent()),
                       std::fabs(collision.partnerDOcc()), collision.partnerDt());
        const float sameSign = pcAlpha(a) * pcAlpha(b) > 0.f ? 1.f : 0.f;
        fRegistry.fill(HIST("CoReco/hPairOrigin"), kSample, origin, obs.qinv, obs.kt, sameSign);

        if (isMC) {
          const auto cls = static_cast<float>(clsId);
          fRegistryPairMC.fill(HIST("CoReco/MC/hClassQ"), kSample, cls, obs.qinv, pcTrueQ(a, b), obs.kt);
          fRegistryPairMC.fill(HIST("CoReco/MC/hClassQinvQout"), kSample, cls, obs.qinv, std::fabs(q3[0]), obs.kt);
          if (qaflags.cfgFillNoteQA.value) {
            fillCoRecoNote(kSample, 1.f, a, b, obs, centForQA, occupancy);
          }
        }
        if (pool) {
          if (usedPhotonIdsPerCol.insert(static_cast<int>(ids[i])).second) {
            emh1->AddTrackToEventPool(keyDFCollision, pwl1);
          }
          if (usedPhotonIdsPerCol.insert(static_cast<int>(ids[j])).second) {
            emh1->AddTrackToEventPool(keyDFCollision, pwl2);
          }
        }
      }
    }
    usedPhotonIdsPerCol.clear();
    return ndiphoton;
  }

  template <typename TCollision>
  void mixStandardCoReco(TCollision const& collision, std::tuple<int, int, int, int> const& keyBin, const float occupancy)
  {
    const auto keyDFCollision = std::make_pair(ndf, static_cast<int>(collision.globalIndex()));
    const auto& selectedPhotons = emh1->GetTracksPerCollision(keyDFCollision);
    const auto& poolIDs = emh1->GetCollisionIdsFromEventPool(keyBin);
    const PairSepConfig sepConfig = pairSepConfig();
    for (const auto& mixID : poolIDs) {
      if (mixID == keyDFCollision) {
        continue;
      }
      const auto& poolPhotons = emh1->GetTracksPerCollision(mixID);
      for (const auto& g1 : selectedPhotons) {
        for (const auto& g2 : poolPhotons) {
          // ─── kinematic gates ──────────────────────────────────
          if (!passAsymmetryCut(g1.pt(), g2.pt())) {
            continue;
          }
          if (!passFastQinvGate(g1, g2)) {
            continue;
          }
          auto obs = buildPairQAObservables(g1, g2);
          if (!obs.valid) {
            continue;
          }
          const bool doFR = passQinvFullRangeGate(obs.qinv);
          // ─── QA before the pair cuts ──────────────────────────
          PairSep sep;
          if (needsPairSep(obs.qinv)) {
            sep = computePairSep(g1, g2, sepConfig);
          }
          fillPairSep<1>(sep, obs);
          // ─── pair cuts ────────────────────────────────────────
          if (!passPairCuts(g1, g2, obs, sep)) {
            continue;
          }
          // ─── fills after all pair cuts ────────────────────────
          if (doFR) {
            fillFullRangeQA<1>(obs, occupancy);
          }
          fillPairHistogram<1>(collision, obs.v1, obs.v2, 1.f);
          if (hbtanalysis.cfgFillPipe.value) {
            fillPipePair(1.f, nullptr, nullptr, obs, 0.f, 0, 0.f, occupancy, nullptr); // pooled photons: no leg records, kinematic maps only
          }
          if (isMC) {
            fRegistryPairMC.fill(HIST("CoReco/MC/hClassQinvQout"), 1.f, (g1.fIsTruePhoton == 1 && g2.fIsTruePhoton == 1) ? static_cast<float>(kClsTrueDistinct) : 0.f, obs.qinv, std::fabs(lcmsQ(obs.v1, obs.v2)[0]), obs.kt);
          }
        }
      }
    }
    emh1->AddCollisionIdAtLast(keyBin, keyDFCollision);
  }

  template <typename TCollision>
  void fillCorrelationEvent(TCollision const& collision, std::vector<int64_t> const& ids, MyV0PhotonsMix const& photons, std::vector<PcPhoton> const& pc, aod::PCMMixLegs const& mixLegs)
  {
    const bool standardMix = mixing.cfgCoRecoMixMode.value == 1;
    const bool mixedEvent = collision.isMixed();
    fRegistry.fill(HIST("CoReco/hArtificialEvents"), mixedEvent ? 2.f : 0.f);
    fRegistry.fill(HIST("CoReco/hNPhotons"), mixedEvent ? 1.f : 0.f, static_cast<float>(ids.size()));
    const float centForQA = eventCentrality(collision);
    if (!isSelectedCentrality(centForQA)) {
      return;
    }
    fRegistry.fill(HIST("CoReco/hArtificialEvents"), mixedEvent ? 3.f : 1.f);
    if (qaflags.cfgFillNoteQA.value) {
      fRegistry.fill(HIST("CoReco/Note/hEvents"), mixedEvent ? 1.f : 0.f, centForQA, collision.occupancy());
    }
    if (mixedEvent && standardMix) {
      return;
    }
    const float occupancy = collision.occupancy();
    if (mixedEvent) {
      pairCoRecoEvent<1>(collision, ids, photons, pc, mixLegs, centForQA, occupancy, false);
      return;
    }
    const int ndiphoton = pairCoRecoEvent<0>(collision, ids, photons, pc, mixLegs, centForQA, occupancy, standardMix && mixing.cfgDoMix.value);
    if (standardMix && mixing.cfgDoMix.value && ndiphoton > 0) {
      const auto keyBin = std::make_tuple(binOf(ztxBinEdges, collision.posZ()), binOf(centBinEdges, centForQA), 0, binOf(occBinEdges, occupancy));
      mixStandardCoReco(collision, keyBin, occupancy);
    }
  }

  bool pcApplyPartnerRule(const float sample, PcPhoton const& a, PcPhoton const& b)
  {
    const int partner = pcPartner(a, b);
    const bool kept = pcPartnerKept(partner, mixing.cfgCoRecoPartnerPairing.value);
    if (partner > 0) {
      const float q = o2::pwgem::photonmeson::qInv(a.px, a.py, a.pz, b.px, b.py, b.pz);
      if (q < qaflags.cfgPairClassMaxQ.value) {
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hPartner"), sample, static_cast<float>(pcPairClass(a, b)), q, static_cast<float>(partner), kept ? 1.f : 0.f);
      }
    }
    return kept;
  }

  void fillPcHalfSwapStray(const float sample, const float q, PcPhoton const& a, PcPhoton const& b, const float nPass)
  {
    const std::array<std::pair<const PcLeg*, const PcLeg*>, 2> legsA{{{&a.pos, &a.ele}, {&a.ele, &a.pos}}};
    const std::array<std::pair<const PcLeg*, const PcLeg*>, 2> legsB{{{&b.pos, &b.ele}, {&b.ele, &b.pos}}};
    for (const auto& [x, xo] : legsA) {
      for (const auto& [y, yo] : legsB) {
        if (x->mpdg == kPcGamma && pcSameMother(*x, *y)) { // x, y = the two legs of the shared photon; xo, yo = the stray legs
          const bool sameEvt = xo->evt == yo->evt && xo->pid >= 0 && yo->pid >= 0;
          const float strays = (sameEvt && pcSameMother(*xo, *yo)) ? 1.f : ((sameEvt && xo->anc >= 0 && xo->anc == yo->anc) ? 2.f : 0.f);
          fRegistryPairMC.fill(HIST("PairClass/HalfSwap/hStray"), sample, q, static_cast<float>(pcStrayRelation(*xo, *x, *y)), static_cast<float>(pcStrayRelation(*yo, *x, *y)), strays, mPcCent, nPass, static_cast<float>(static_cast<int>(xo->free) + static_cast<int>(yo->free)));
          return;
        }
      }
    }
  }

  void fillPcNote(const float sample, const float cls, const float q, PcPhoton const& a, PcPhoton const& b, const float nPass, const float sepPass)
  {
    if (!qaflags.cfgFillNoteQA.value || q >= kNoteMaxQ) {
      return;
    }
    fRegistryPairMC.fill(HIST("PairClass/Note/hClassCentOcc"), sample, cls, q, mPcCent, mPcOcc, nPass, sepPass);
    const float sharedMax = std::max({a.pos.sharedFrac, a.ele.sharedFrac, b.pos.sharedFrac, b.ele.sharedFrac});
    const float nClsMin = std::min({a.pos.nCls, a.ele.nCls, b.pos.nCls, b.ele.nCls});
    fRegistryPairMC.fill(HIST("PairClass/Note/hLegCls"), sample, cls, q, std::min(sharedMax, 0.999f), std::min(nClsMin, 159.f), nPass, sepPass);
    const int icls = static_cast<int>(cls);
    if (icls == kClsSamePhoton || icls == kClsSharedLegPhotons || icls == kClsSharedLeg || icls == kClsSharedTrack) { // duplicate classes
      const auto [x, y, slots] = pcSharedPieces(a, b, icls == kClsSharedTrack);
      if (x != nullptr) {
        const float clsBin = icls == kClsSharedTrack ? 3.f : static_cast<float>(icls - 1); // o2-linter: disable=magic-number (bin of class 12)
        fRegistryPairMC.fill(HIST("PairClass/Note/hDupCls"), sample, clsBin, q, std::min(std::min(x->nCls, y->nCls), 159.f), std::min(std::max(x->nCls, y->nCls), 159.f),
                             std::min(std::max(x->sharedFrac, y->sharedFrac), 0.999f), static_cast<float>(slots), nPass, sepPass);
      }
    }
    for (const auto* p : {&a, &b}) {
      const float pt = std::hypot(p->px, p->py);
      const float eta = std::asinh(p->pz / std::max(pt, 1.e-6f));
      const float phi = RecoDecay::constrainAngle(std::atan2(p->py, p->px), 0.f);
      fRegistryPairMC.fill(HIST("PairClass/Note/hPhotonKine"), sample, cls, q, std::clamp(eta, -0.999f, 0.999f), phi, std::min(pt, 9.99f), nPass, sepPass);
      fRegistryPairMC.fill(HIST("PairClass/Note/hPhotonConv"), sample, cls, q, std::min(std::hypot(p->vx, p->vy), 89.9f), std::clamp(p->vz, -119.9f, 119.9f), std::clamp(mPcVtxZ, -9.99f, 9.99f), nPass, sepPass);
    }
  }

  void fillPcPair(const float sample, PcPhoton const& a, PcPhoton const& b, aod::PCMMixLegs const& mixLegs)
  {
    const float q = o2::pwgem::photonmeson::qInv(a.px, a.py, a.pz, b.px, b.py, b.pz);
    if (q >= qaflags.cfgPairClassMaxQ.value) {
      return;
    }
    const float clean = std::max({a.pos.sharedFrac, a.ele.sharedFrac, b.pos.sharedFrac, b.ele.sharedFrac}) < qaflags.cfgPairClassMaxSharedFrac.value ? 1.f : 0.f;
    const int clsId = pcPairClass(a, b);
    const auto cls = static_cast<float>(clsId);
    fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairClass"), sample, cls, q, static_cast<float>(pcPairOrigin(a, b)), clean);
    const auto nPass = static_cast<float>(static_cast<int>(a.pass) + static_cast<int>(b.pass));
    if (q < 0.1f) { // o2-linter: disable=magic-number (diagnostic range)
      fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairDens"), sample, cls, q, mPcNLegs, mPcCent);
    }
    fRegistryPairMC.fill(HIST("PairClass/Pairs/hKtQ"), sample, cls, q, 0.5f * std::hypot(a.px + b.px, a.py + b.py));
    if (q < 0.1f) { // o2-linter: disable=magic-number (diagnostic range)
      const auto sep = pcPairSep(a, b, mBzT);
      const float eAs = std::hypot(a.px, a.py, a.pz), eBs = std::hypot(b.px, b.py, b.pz);
      constexpr float kRTpcEntranceCm = 70.f; // conversions in front of the TPC (inner vessel / field cage): swaps the overlay misses
      const auto rRegion = static_cast<float>((std::hypot(a.vx, a.vy) >= kRTpcEntranceCm ? 1 : 0) + (std::hypot(b.vx, b.vy) >= kRTpcEntranceCm ? 1 : 0));
      fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairSep"), sample, cls, q, std::min(sep.dMin, 9.99f), std::min(sep.dRPhiMin, 4.99f), std::min(std::fabs(eAs - eBs), 0.499f), rRegion, clean, nPass);
      fillPcNote(sample, cls, q, a, b, nPass, sep.dMin > pairsep.cfgDPairCut.value ? 1.f : 0.f);
      if (sample != 1.f) { // not the standard ME: its pool photons can come from earlier DFs (leg rows of another PCMMixLegs table)
        const float dPoint = pcBremsPointingDMin(a, b, &mixLegs, mBzT);
        const float dCat = pcPointingCategory(dPoint);
        const float kTt = 0.5f * std::hypot(a.px + b.px, a.py + b.py);
        const float sepPass = sep.dMin > pairsep.cfgDPairCut.value ? 1.f : 0.f;
        // |q_out| in the LCMS (as PairClass/Pairs/hQinvQout): the tagged sample in the fit plane (q_out, Q_inv)
        const ROOT::Math::PxPyPzMVector vTagA(a.px, a.py, a.pz, 0.f), vTagB(b.px, b.py, b.pz, 0.f);
        const auto sumTag = vTagA + vTagB;
        const ROOT::Math::Boost bstTag(0., 0., -sumTag.Pz() / sumTag.E());
        const float qOutTag = std::fabs(bstTag(ROOT::Math::PxPyPzEVector(vTagA - vTagB)).Vect().Dot(ROOT::Math::XYZVector(sumTag.Px(), sumTag.Py(), 0.).Unit()));
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairTag"), sample, cls, q, std::min(std::fabs(eAs - eBs), 9.99f), std::min(kTt, 9.99f), mPcCent, dCat, nPass, sepPass, std::min(qOutTag, 9.99f));
        if (qaflags.cfgFillHalfSwapQA.value) {
          const auto sCross = pcCrossScores(a, b, &mixLegs, mBzT, crosspair.cfgCrossMaxGap.value, crosspair.cfgCrossMaxDz.value, crosspair.cfgCrossMaxDTgl.value);
          fRegistryPairMC.fill(HIST("PairClass/HalfSwap/hCrossBoth"), sample, cls, q, std::min(sCross[0], 999.f), std::min(sCross[1], 999.f), nPass, sepPass);
        }
        if (clsId == kClsBremsFake) {
          // x = leg of the radiated photon, y = the radiating e+- (x's grandmother); xo, yo = the other legs of the two photons
          const PcLeg *x = nullptr, *y = nullptr, *xo = nullptr, *yo = nullptr;
          const PcPhoton *pB = nullptr, *pA = nullptr; // pB holds x (radiated side), pA holds y (parent side)
          const std::array<std::pair<const PcPhoton*, const PcPhoton*>, 2> orders{{{&a, &b}, {&b, &a}}};
          for (const auto& [p1, p2] : orders) {
            const std::array<std::pair<const PcLeg*, const PcLeg*>, 2> legs1{{{&p1->pos, &p1->ele}, {&p1->ele, &p1->pos}}};
            const std::array<std::pair<const PcLeg*, const PcLeg*>, 2> legs2{{{&p2->pos, &p2->ele}, {&p2->ele, &p2->pos}}};
            for (const auto& [l1, l1o] : legs1) {
              for (const auto& [l2, l2o] : legs2) {
                if (!x && pcRadiatedBy(*l1, *l2)) {
                  x = l1;
                  xo = l1o;
                  y = l2;
                  yo = l2o;
                  pB = p1;
                  pA = p2;
                }
              }
            }
          }
          float side = 3.f, foreign = 5.f, fakeOrigin = 0.f;
          if (x) {
            const bool bTrue = pcTruePhoton(*pB), aTrue = pcTruePhoton(*pA);
            side = (bTrue == aTrue) ? 2.f : (bTrue ? 1.f : 0.f);
            const PcLeg* f = !bTrue ? xo : yo; // foreign leg of the fake photon (of the radiated side if both are fakes)
            const PcPhoton* pf = !bTrue ? pB : pA;
            fakeOrigin = pf->original ? (pf->svertexer ? 0.f : 1.f) : 2.f;
            const std::array<const PcLeg*, 3> others{!bTrue ? x : y, !bTrue ? y : x, !bTrue ? yo : xo};
            bool cascade = false;
            for (const auto* l : others) {
              if (f->evt != l->evt || f->pid < 0 || l->pid < 0) {
                continue;
              }
              cascade = cascade || (f->gmid >= 0 && f->gmid == l->pid) || (f->mid >= 0 && f->mid == l->pid) || f->pid == l->gmid || f->pid == l->mid || (f->mid >= 0 && f->mid == l->mid);
            }
            if (f->pid < 0) {
              foreign = 4.f;
            } else if (std::abs(f->pdg) != kPcElectron) {
              foreign = 3.f;
            } else if (cascade) {
              foreign = 1.f;
            } else if (f->mpdg == kPcGamma) {
              foreign = 0.f;
            } else {
              foreign = 2.f;
            }
          } else {
            const PcPhoton& pf = !pcTruePhoton(a) ? a : b;
            fakeOrigin = pf.original ? (pf.svertexer ? 0.f : 1.f) : 2.f;
          }
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hBremsFake"), sample, q, side, foreign, mPcCent, dCat, nPass, sepPass, fakeOrigin);
        }
      }
    }
    if (clsId == kClsHalfSwap && q < 0.1f && qaflags.cfgFillHalfSwapQA.value) { // o2-linter: disable=magic-number (diagnostic range)
      fillPcHalfSwapStray(sample, q, a, b, nPass);
    }
    if (clsId == kClsHalfSwap) {
      for (const auto* x : {&a.pos, &a.ele}) {
        for (const auto* y : {&b.pos, &b.ele}) {
          if (x->mpdg == kPcGamma && pcSameMother(*x, *y)) {
            const float src = (x->free && y->free) ? 0.f : (x->free || y->free) ? 3.f
                                                                                : (x->srcPhoton == y->srcPhoton ? 1.f : 2.f);
            fRegistryPairMC.fill(HIST("PairClass/Pairs/hHalfSwapLegs"), sample, q, src, static_cast<float>((a.nFree > 0 ? 1 : 0) + (b.nFree > 0 ? 1 : 0)));
          }
        }
      }
    }
    {
      const float eA = std::sqrt(a.px * a.px + a.py * a.py + a.pz * a.pz), eB = std::sqrt(b.px * b.px + b.py * b.py + b.pz * b.pz);
      fRegistryPairMC.fill(HIST("PairClass/Pairs/hDEQinv"), sample, cls, q, std::fabs(eA - eB), 0.5f * std::hypot(a.px + b.px, a.py + b.py));
    }
    fRegistryPairMC.fill(HIST("PairClass/Pairs/hFreeLegPair"), sample, cls, q, static_cast<float>((a.nFree > 0 ? 1 : 0) + (b.nFree > 0 ? 1 : 0)), static_cast<float>((a.original ? 0 : 1) + (b.original ? 0 : 1)));
    {
      const ROOT::Math::PxPyPzMVector v1(a.px, a.py, a.pz, 0.f), v2(b.px, b.py, b.pz, 0.f);
      const auto sum = v1 + v2;
      const ROOT::Math::Boost bstZ(0., 0., -sum.Pz() / sum.E());
      const auto q3 = bstZ(ROOT::Math::PxPyPzEVector(v1 - v2)).Vect();
      const float qout = std::fabs(q3.Dot(ROOT::Math::XYZVector(sum.Px(), sum.Py(), 0.).Unit()));
      const float kt = 0.5f * std::hypot(a.px + b.px, a.py + b.py);
      fRegistryPairMC.fill(HIST("PairClass/Pairs/hQinvQout"), sample, cls, q, qout, std::min(kt, 2.99f));
      const auto ktUnit = ROOT::Math::XYZVector(sum.Px(), sum.Py(), 0.).Unit();
      const float qside = std::fabs(q3.Dot(ROOT::Math::XYZVector(-ktUnit.Y(), ktUnit.X(), 0.)));
      const float qlong = std::fabs(q3.Z());
      auto etaOf = [](PcPhoton const& p) { return std::asinh(p.pz / std::max(std::hypot(p.px, p.py), 1.e-6f)); };
      const float dEta = std::fabs(etaOf(a) - etaOf(b));
      const float dPhi = std::fabs(RecoDecay::constrainAngle(std::atan2(a.py, a.px) - std::atan2(b.py, b.px), -o2::constants::math::PI));
      const float psi = std::atan2(dEta, dPhi);
      if (q < 0.05f) { // o2-linter: disable=magic-number (diagnostic range)
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairOpening"), sample, cls, q, std::min(qside, 0.149f), std::min(qlong, 0.149f), std::min(psi, 1.57f));
      }
      if (q < 0.1f) { // o2-linter: disable=magic-number (diagnostic range)
        const float eA = std::hypot(a.px, a.py, a.pz), eB = std::hypot(b.px, b.py, b.pz);
        PcPhoton const& soft = eA <= eB ? a : b;
        PcPhoton const& hard = eA <= eB ? b : a;
        const float rSoft = std::hypot(soft.vx, soft.vy), rHard = std::hypot(hard.vx, hard.vy);
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hSoftHard"), sample, cls, q, std::min(std::min(eA, eB) / std::max(std::max(eA, eB), 1.e-6f), 0.999f), std::clamp(rSoft - rHard, -59.9f, 59.9f), std::min(rHard, 89.9f), std::min(psi, 1.57f));
        if (pcTruePhoton(a) && pcTruePhoton(b)) {
          const float tA = std::hypot(a.pos.mpx, a.pos.mpy, a.pos.mpz), tB = std::hypot(b.pos.mpx, b.pos.mpy, b.pos.mpz);
          const float qT = o2::pwgem::photonmeson::qInv(a.pos.mpx, a.pos.mpy, a.pos.mpz, b.pos.mpx, b.pos.mpy, b.pos.mpz);
          const float dERec = eA >= eB ? eA - eB : eB - eA;
          const float dETrue = eA >= eB ? tA - tB : tB - tA; // same photon order as the reconstructed one
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairRes"), sample, cls, q, std::clamp(q - qT, -0.0199f, 0.0199f), std::clamp(dERec - dETrue, -0.0999f, 0.0999f), mPcCent);
        }
      }
    }
    fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairDist"), sample, cls, q, std::hypot(a.vx - b.vx, a.vy - b.vy, a.vz - b.vz), pcTrueDist({&a.pos, &a.ele, &b.pos, &b.ele}));
    // four-leg invariant mass (pairing invariant) and the true q of the two mother photons
    ROOT::Math::PxPyPzMVector sum4;
    const std::array<const PcLeg*, 4> legs4{&a.pos, &a.ele, &b.pos, &b.ele};
    for (const auto* l : legs4) {
      sum4 += ROOT::Math::PxPyPzMVector(l->px, l->py, l->pz, o2::constants::physics::MassElectron);
    }
    float qTrue = -1.f;
    {
      std::vector<const PcLeg*> mothers;
      for (const auto* l : legs4) {
        if (l->mid < 0 || l->mpdg != kPcGamma) {
          continue;
        }
        bool seen = false;
        for (const auto* m : mothers) {
          seen = seen || pcSameMother(*l, *m);
        }
        if (!seen) {
          mothers.push_back(l);
        }
      }
      if (mothers.size() == 2) {
        qTrue = o2::pwgem::photonmeson::qInv(mothers[0]->mpx, mothers[0]->mpy, mothers[0]->mpz, mothers[1]->mpx, mothers[1]->mpy, mothers[1]->mpz);
      }
    }
    fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairQ"), sample, cls, q, static_cast<float>(sum4.M()), qTrue, 0.5f * std::hypot(a.px + b.px, a.py + b.py));
    if (clsId == kClsTrueDistinct && q < 0.3f) { // o2-linter: disable=magic-number (histogram range)
      // decay relation of two different true photons: leg -> photon (mid) -> photon's mother (gmid) -> its mother (ggid)
      float rel = 0.f;
      int common = 0;
      const PcLeg& la = a.pos;
      const PcLeg& lb = b.pos;
      if (la.evt == lb.evt) {
        if (la.ggid >= 0 && la.ggid == lb.ggid) {
          rel = 1.f;
          common = la.ggpdg;
        } else if (la.gmid >= 0 && la.gmid == lb.ggid) {
          rel = 2.f;
          common = la.gpdg;
        } else if (lb.gmid >= 0 && lb.gmid == la.ggid) {
          rel = 2.f;
          common = lb.gpdg;
        } else if (la.anc >= 0 && la.anc == lb.anc) {
          rel = 3.f;
        } else if (la.root >= 0 && la.root == lb.root) {
          rel = 4.f;
          common = la.rootPdg;
        }
      }
      float cp = 6.f; // o2-linter: disable=magic-number (category codes)
      if (rel > 0.f) {
        const int ac = std::abs(common);
        cp = 5.f;                                  // o2-linter: disable=magic-number (category codes)
        if (ac == 310 || ac == 130 || ac == 311) { // o2-linter: disable=pdg/explicit-code (K0)
          cp = 0.f;
        } else if (ac == 221) { // o2-linter: disable=pdg/explicit-code (eta)
          cp = 1.f;
        } else if (ac == 223) { // o2-linter: disable=pdg/explicit-code (omega)
          cp = 2.f;
        } else if (ac == 331) { // o2-linter: disable=pdg/explicit-code (eta')
          cp = 3.f;
        } else if (ac > 1000 && ac < 10000) { // o2-linter: disable=magic-number (baryons)
          cp = 4.f;
        }
      }
      fRegistryPairMC.fill(HIST("PairClass/Pairs/hTrueRelation"), sample, q, rel, cp, std::min(0.5f * std::hypot(a.px + b.px, a.py + b.py), 2.99f));
      if (q < 0.25f) { // o2-linter: disable=magic-number (histogram range)
        const float etaA = std::asinh(a.pz / std::max(std::hypot(a.px, a.py), 1e-6f)), etaB = std::asinh(b.pz / std::max(std::hypot(b.px, b.py), 1e-6f));
        float dphi = std::fabs(std::atan2(a.py, a.px) - std::atan2(b.py, b.px));
        if (dphi > o2::constants::math::PI) {
          dphi = o2::constants::math::TwoPI - dphi;
        }
        const int nPrim = ((la.anc >= 0 && la.anc == la.mid) ? 1 : 0) + ((lb.anc >= 0 && lb.anc == lb.mid) ? 1 : 0);
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hTrueAngle"), sample, q, std::min(std::fabs(etaA - etaB), 0.799f), std::min(dphi, 0.799f), static_cast<float>(nPrim));
        {
          const bool swapOrder = (static_cast<int>(std::fabs(etaA + etaB) * 1.e6f) & 1) != 0; // pseudo-random order of the two photons
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hTrueEta12"), sample, std::clamp(swapOrder ? etaB : etaA, -0.799f, 0.799f), std::clamp(swapOrder ? etaA : etaB, -0.799f, 0.799f), std::max(mPcOcc, 0.f));
        }
        if (sample < 0.5f && nPrim == 2) { // o2-linter: disable=magic-number (same event, both photons primary)
          auto originZ = [](PcLeg const& l) { return l.vz - (l.vx * l.mpx + l.vy * l.mpy) / std::max(l.mpx * l.mpx + l.mpy * l.mpy, 1e-9f) * l.mpz; };
          const float zA = originZ(la), zB = originZ(lb);
          const bool forA = std::fabs(zA - mPcVtxZ) > 0.5f, forB = std::fabs(zB - mPcVtxZ) > 0.5f;         // o2-linter: disable=magic-number (cm)
          const int nFar = static_cast<int>(forA) + static_cast<int>(forB);                                // photons whose true origin is more than 0.5 cm from the vertex
          const float cat = nFar == 0 ? 0.f : (nFar == 1 ? 1.f : (std::fabs(zA - zB) < 0.1f ? 2.f : 3.f)); // o2-linter: disable=magic-number (category codes, cm)
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hTrueOrigin"), q, std::max(mPcOcc, 0.f), cat, etaA * etaB > 0.f ? 1.f : 0.f);
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hTrueOriginDz"), std::clamp(zA - mPcVtxZ, -29.99f, 29.99f), std::clamp(etaA, -0.799f, 0.799f), std::max(mPcOcc, 0.f));
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hTrueOriginDz"), std::clamp(zB - mPcVtxZ, -29.99f, 29.99f), std::clamp(etaB, -0.799f, 0.799f), std::max(mPcOcc, 0.f));
        }
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hTrueAngleEvt"), sample, std::min(std::fabs(etaA - etaB), 0.799f), std::min(std::fabs(mPcVtxZ), 9.99f), std::max(mPcOcc, 0.f), std::clamp(0.5f * (etaA + etaB), -0.799f, 0.799f));
      }
    }
    if (q < 0.1f) {
      const int icls = clsId;
      if (icls == kClsSamePhoton || icls == kClsSharedLegPhotons || icls == kClsSharedLeg || icls == kClsSharedTrack) { // duplicate classes
        const std::array<const PcLeg*, 2> la{&a.pos, &a.ele}, lb{&b.pos, &b.ele};
        const PcLeg* x = nullptr;
        const PcLeg* y = nullptr;
        int slots = 2;
        auto shares = [icls](PcLeg const& u, PcLeg const& v) {
          return icls == kClsSharedTrack ? (u.evt == v.evt && u.trk >= 0 && u.trk == v.trk) : pcSamePart(u, v);
        };
        for (int i = 0; i < 2 && x == nullptr; ++i) { // o2-linter: disable=magic-number (two legs)
          for (int j = 0; j < 2; ++j) {               // o2-linter: disable=magic-number (two legs)
            if (shares(*la[i], *lb[j])) {
              x = la[i];
              y = lb[j];
              slots = i == j ? 0 : 1;
              break;
            }
          }
        }
        float dTglPieces = 0.199f, dPtPieces = 0.999f, nItsPieces = 0.f;
        if (x != nullptr) {
          const float pt1 = std::hypot(x->px, x->py), pt2 = std::hypot(y->px, y->py);
          if (pt1 > 0.f && pt2 > 0.f) {
            dTglPieces = std::min(std::fabs(x->pz / pt1 - y->pz / pt2), 0.199f);
            dPtPieces = std::min(std::fabs(pt1 - pt2) / std::max(pt1, pt2), 0.999f);
          }
          nItsPieces = static_cast<float>((x->its ? 1 : 0) + (y->its ? 1 : 0));
        }
        const float clsBin = icls == kClsSharedTrack ? 3.f : static_cast<float>(icls - 1); // o2-linter: disable=magic-number (bin of class 12)
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hDupPieces"), sample, clsBin, q, dTglPieces, dPtPieces, static_cast<float>(slots), static_cast<float>(pcPairOrigin(a, b)), nItsPieces, nPass);
      }
      {
        const float rA = std::hypot(a.vx, a.vy), rB = std::hypot(b.vx, b.vy);
        PcPhoton const& inner = rA <= rB ? a : b;
        const float rTrue = inner.pos.pid >= 0 ? std::hypot(inner.pos.vx, inner.pos.vy) : (inner.ele.pid >= 0 ? std::hypot(inner.ele.vx, inner.ele.vy) : -0.5f);
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hPairR"), sample, cls, q, std::min(std::min(rA, rB), 89.9f), std::min(std::max(rA, rB), 89.9f), std::min(rTrue, 89.9f));
      }
      {
        auto radiated = [](PcPhoton const& p) { return p.pos.mpdg == kPcGamma && p.pos.mid >= 0 && std::abs(p.pos.gpdg) == kPcElectron; };
        for (const auto* p : {&a, &b}) {
          const float ang = std::acos(std::clamp(p->cospa, -1.f, 1.f));
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hPhotonQuality"), sample, cls, q, radiated(*p) ? 0.f : 1.f, std::min(p->chi2, 49.9f), std::min(ang, 0.299f), std::min(std::fabs(p->dcaxy), 4.99f), std::min(std::fabs(p->dcaz), 4.99f));
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hPhotonKF"), sample, cls, q, radiated(*p) ? 0.f : 1.f, std::clamp(p->chi2Vtx, 0.f, 24.99f), std::clamp(p->chi2PV, 0.f, 299.f), std::clamp(p->psipair, -0.199f, 0.199f));
        }
        if (clsId == kClsTrueDistinct && (radiated(a) || radiated(b))) {
          const PcPhoton& r = radiated(a) ? a : b;
          const PcPhoton& o = radiated(a) ? b : a;
          float rel = 0.f;
          if (radiated(o)) {
            rel = 3.f;
          } else if (o.pos.evt == r.pos.evt && o.pos.gmid >= 0 && o.pos.gmid == r.pos.ggid) {
            rel = 1.f; // partner photon and the radiating electron have the same mother
          } else if (o.pos.evt == r.pos.evt && o.pos.root >= 0 && o.pos.root == r.pos.root) {
            rel = 2.f;
          }
          // where the radiating electron comes from: the mother of the photon's mother (the electron) = ggpdg of the photon's leg
          const int em = std::abs(r.pos.ggpdg);
          float radFrom = 5.f; // o2-linter: disable=magic-number (category codes)
          if (em == kPcGamma) {
            radFrom = 0.f;
          } else if (em == 111) { // o2-linter: disable=pdg/explicit-code (pi0)
            radFrom = 1.f;
          } else if (em == 221) { // o2-linter: disable=pdg/explicit-code (eta)
            radFrom = 2.f;
          } else if ((em >= 400 && em < 600) || (em >= 4000 && em < 6000)) { // o2-linter: disable=magic-number (charm and beauty hadrons)
            radFrom = 3.f;
          } else if (em > 100) { // o2-linter: disable=magic-number (hadrons)
            radFrom = 4.f;
          }
          // first physical primary of the brems photon (via its mother electron) vs that of the partner photon
          const float sameAnc = (r.pos.anc < 0 || o.pos.anc < 0 || r.pos.evt != o.pos.evt) ? 2.f : (r.pos.anc == o.pos.anc ? 1.f : 0.f);
          fRegistryPairMC.fill(HIST("PairClass/Pairs/hHiddenBrems"), sample, q, static_cast<float>((radiated(a) ? 1 : 0) + (radiated(b) ? 1 : 0)), rel, radFrom, sameAnc);
        }
      }
      if (clsId == kClsBrems || clsId == kClsBremsFake) {
        // x = leg of the bremsstrahlung photon, y = the radiating electron (a leg of the partner photon): x's grandmother is y
        const PcLeg *x = nullptr, *y = nullptr;
        const std::array<std::pair<const PcPhoton*, const PcPhoton*>, 2> orders{{{&a, &b}, {&b, &a}}}; // both photon orders
        for (const auto& [p1, p2] : orders) {
          for (const auto* l1 : {&p1->pos, &p1->ele}) {
            for (const auto* l2 : {&p2->pos, &p2->ele}) {
              if (!x && pcRadiatedBy(*l1, *l2)) {
                x = l1;
                y = l2;
              }
            }
          }
        }
        auto pdgClass = [](int pdg) { pdg = std::abs(pdg); return pdg == 111 ? 0 : (pdg == 221 ? 1 : (pdg == 223 ? 2 : (pdg == 11 ? 3 : (pdg == 0 ? 5 : 4)))); }; // o2-linter: disable=magic-number (PDG codes)
        auto angle = [](float ax, float ay, float az, float bx, float by, float bz) {
          const float n = std::sqrt((ax * ax + ay * ay + az * az) * (bx * bx + by * by + bz * bz));
          return n > 0.f ? std::acos(std::clamp((ax * bx + ay * by + az * bz) / n, -1.f, 1.f)) : 0.249f;
        };
        float origin = 4.f, parentMother = 5.f, thE = 0.249f, thG = 0.249f, dR = -39.9f, eFrac = 0.999f;
        if (x && y) {
          const int mp = std::abs(y->mpdg);
          origin = mp == kPcGamma ? 0.f : (mp == 111 ? 1.f : (mp == 221 ? 2.f : 3.f)); // o2-linter: disable=magic-number (PDG codes)
          parentMother = mp == kPcGamma ? static_cast<float>(pdgClass(y->gpdg)) : 5.f;
          thE = std::min(angle(x->mpx, x->mpy, x->mpz, y->px, y->py, y->pz), 0.249f);
          if (mp == kPcGamma) {
            thG = std::min(angle(x->mpx, x->mpy, x->mpz, y->mpx, y->mpy, y->mpz), 0.249f);
            const float eB = std::hypot(x->mpx, x->mpy, x->mpz), eA = std::hypot(y->mpx, y->mpy, y->mpz);
            eFrac = eA > 0.f ? std::min(eB / eA, 0.999f) : 0.999f;
          }
          dR = std::clamp(std::hypot(x->vx, x->vy) - std::hypot(y->vx, y->vy), -39.9f, 79.9f);
        }
        fRegistryPairMC.fill(HIST("PairClass/Pairs/hBrems"), sample, clsId == kClsBrems ? 0.f : 1.f, q, origin, parentMother, thE, thG, dR, eFrac);
        if (x && y && std::abs(y->mpdg) == kPcGamma) {
          auto ancClass = [](int pdg) { // o2-linter: disable=pdg/explicit-code (diagnostic classes)
            pdg = std::abs(pdg);
            switch (pdg) {
              case 111:
                return 0.f;
              case 221:
                return 1.f;
              case 223:
                return 2.f;
              case 331:
                return 3.f;
              case 113:
              case 213:
              case 313:
              case 323:
                return 4.f;
              case 310:
              case 130:
              case 311:
                return 5.f;
              default:
                break;
            }
            if (pdg == 11 || pdg == 13 || pdg == 22) {
              return 8.f; // o2-linter: disable=magic-number (leptons, photon: other)
            }
            if (pdg > 1000) {
              return 6.f; // o2-linter: disable=magic-number (baryons)
            }
            if (pdg < 100) {
              return 7.f; // o2-linter: disable=magic-number (no hadron mother: quark, gluon, string, cluster, none)
            }
            return 8.f;
          };
          // a particle whose mother is the root of the decay chain (the beam record) comes directly from the collision: "primary" (7)
          const float ancMother = (y->gmid >= 0 && y->gmid == y->root) ? 7.f : ancClass(y->gpdg); // o2-linter: disable=magic-number (class)
          const float ancGrand = (y->ggid >= 0 && y->ggid == y->root) ? 7.f : ancClass(y->ggpdg); // o2-linter: disable=magic-number (class)

          fRegistryPairMC.fill(HIST("PairClass/Pairs/hBremsAncestry"), sample, clsId == kClsBrems ? 0.f : 1.f, q, ancMother, ancGrand, ancClass(y->rootPdg));
        }
      }
    }
  }

  template <typename TCollision>
  void fillPairClassQAEvent(TCollision const& collision, std::vector<int64_t> const& ids, MyV0PhotonsMix const& photons, std::vector<PcPhoton> const& pc, aod::PCMMixLegs const& mixLegs)
  {
    mPcCent = collision.centFT0C(); // standard ME: centrality of the current event (pool partners are in the same centrality bin)
    mPcNLegs = static_cast<float>(collision.nLegs());
    mPcOcc = collision.occupancy();
    mPcVtxZ = collision.posZ();
    if (qaflags.cfgFillNoteQA.value) {
      fRegistryPairMC.fill(HIST("PairClass/Note/hEvents"), collision.isMixed() ? 1.f : 0.f, mPcCent, mPcOcc);
    }
    std::vector<PcPhoton> list;
    list.reserve(ids.size());
    for (const auto id : ids) {
      list.push_back(pc[id]);
      list.back().pass = isSelectedPhoton(photons, id); // the analysis photon cut, as in the correlation block
    }
    const float type = collision.isMixed() ? 1.f : 0.f;
    fRegistryPairMC.fill(HIST("PairClass/Events/hArtificialEvents"), type);
    fRegistryPairMC.fill(HIST("PairClass/Events/hNPhotons"), type, static_cast<float>(list.size()));
    for (const auto& p : list) {
      const float origin = p.original ? (p.svertexer ? 0.f : 1.f) : 2.f;
      fRegistryPairMC.fill(HIST("PairClass/Photons/hOrigin"), type, origin);
      fRegistryPairMC.fill(HIST("PairClass/Photons/hLegTrueDist"), type, origin, pcTrueDist({&p.pos, &p.ele}));
    }
    for (size_t i = 0; i < list.size(); ++i) {
      for (size_t j = i + 1; j < list.size(); ++j) {
        const auto& a = list[i];
        const auto& b = list[j];
        float sample = 0.f;
        if (collision.isMixed()) {
          const bool oneSource = a.pos.evt == a.ele.evt && a.pos.evt == b.pos.evt && a.pos.evt == b.ele.evt;
          sample = oneSource ? 3.f : 2.f;
        }
        if (!pcApplyPartnerRule(sample, a, b)) {
          continue;
        }
        fillPcPair(sample, a, b, mixLegs);
      }
    }
    if (collision.isMixed()) {
      return;
    }
    const int zBin = static_cast<int>(std::floor((collision.posZ() + 10.f) / qaflags.cfgPairClassMixZBinWidth.value));
    const int cBin = static_cast<int>(std::floor(collision.centFT0C() / qaflags.cfgPairClassMixCentBinWidth.value));
    auto& pool = mPcPools[zBin * 1000 + cBin];
    if (pool.empty()) {
      fRegistryPairMC.fill(HIST("PairClass/Events/hArtificialEvents"), 3.f);
    }
    for (const auto& partner : pool) {
      fRegistryPairMC.fill(HIST("PairClass/Events/hArtificialEvents"), 2.f);
      for (const auto& a : list) {
        for (const auto& b : partner) {
          fillPcPair(1.f, a, b, mixLegs);
        }
      }
    }
    pool.push_back(list);
    if (static_cast<int>(pool.size()) > qaflags.cfgPairClassMixDepth.value) {
      pool.pop_front();
    }
  }

  template <typename TLegs, typename TCollision, typename TPhotons, typename TMCParticles, typename TCut>
  void fillSinglePhotonResponse(TCollision const& collision, TPhotons const& photonsColl, TMCParticles const& mcParticles, TCut const& cut)
  {
    if (collision.has_emmcevent()) {
      auto gen = mcParticles.sliceBy(perMCCollisionEMMCParts, collision.emmceventId());
      for (const auto& m : gen) {
        if (m.pdgCode() != kGamma || !(m.isPhysicalPrimary() || m.producedByGenerator()) || std::fabs(m.eta()) > pcmcuts.cfgMaxEtaV0.value) {
          continue;
        }
        fRegistryPairMC.fill(HIST("SinglePhotonMC/Response/hGenPt"), m.pt());
      }
    }
    auto angleOf = [](float ax, float ay, float az, float bx, float by, float bz) {
      const float n = std::sqrt((ax * ax + ay * ay + az * az) * (bx * bx + by * by + bz * bz));
      return n > 0.f ? std::acos(std::clamp((ax * bx + ay * by + az * bz) / n, -1.f, 1.f)) : 0.f;
    };
    std::unordered_set<int> recoMothers;
    for (const auto& g : photonsColl) {
      if (!cut.template IsSelected<decltype(g), TLegs>(g)) {
        continue;
      }
      const auto mcg = buildPhotonMCInfo<decltype(g), TLegs>(g, mcParticles);
      if (mcg.hasMC && mcg.sameMother && mcg.isTruePhoton) {
        recoMothers.insert(mcg.motherId);
      }
    }
    for (const auto& g : photonsColl) {
      if (!cut.template IsSelected<decltype(g), TLegs>(g)) {
        continue;
      }
      const auto mcg = buildPhotonMCInfo<decltype(g), TLegs>(g, mcParticles);
      if (!(mcg.hasMC && mcg.sameMother && mcg.isTruePhoton)) {
        continue;
      }
      const auto mother = mcParticles.iteratorAt(mcg.motherId);
      const float eTrue = std::hypot(mother.px(), mother.py(), mother.pz());
      if (eTrue <= 0.f) {
        continue;
      }
      const float eRec = std::hypot(g.px(), g.py(), g.pz());
      const float angRec = std::min(angleOf(g.px(), g.py(), g.pz(), mother.px(), mother.py(), mother.pz()), 0.0399f);
      // the hardest photon radiated by one of the two legs, and whether it was reconstructed here
      float eRadFrac = -0.025f;
      float radReco = 0.f;
      for (const auto& legId : {mcg.mcPosId, mcg.mcNegId}) {
        const auto leg = mcParticles.iteratorAt(legId);
        if (!leg.has_daughters()) {
          continue;
        }
        for (const auto& dId : leg.daughtersIds()) {
          const auto dtr = mcParticles.iteratorAt(dId);
          if (dtr.pdgCode() != kGamma) {
            continue;
          }
          const float frac = std::hypot(dtr.px(), dtr.py(), dtr.pz()) / eTrue;
          if (frac > eRadFrac) {
            eRadFrac = std::min(frac, 0.999f);
            radReco = recoMothers.contains(dId) ? 1.f : 0.f;
          }
        }
      }
      // this photon itself radiated by a conversion electron: mother photon <- electron <- parent photon
      float isRad = 0.f, fracParent = -0.025f, angParent = -0.0025f;
      if (mother.has_mothers()) {
        const auto ele = mcParticles.iteratorAt(mother.mothersIds()[0]);
        if (std::abs(ele.pdgCode()) == PDG_t::kElectron && ele.has_mothers()) {
          const auto parent = mcParticles.iteratorAt(ele.mothersIds()[0]);
          if (parent.pdgCode() == kGamma) {
            isRad = recoMothers.contains(parent.globalIndex()) ? 2.f : 1.f;
            const float eParent = std::hypot(parent.px(), parent.py(), parent.pz());
            fracParent = eParent > 0.f ? std::min(eTrue / eParent, 0.999f) : 0.999f;
            angParent = std::min(angleOf(mother.px(), mother.py(), mother.pz(), parent.px(), parent.py(), parent.pz()), 0.249f);
          }
        }
      }
      fRegistryPairMC.fill(HIST("SinglePhotonMC/Response/hResponse"), mother.pt(), std::clamp(eRec / eTrue, 0.301f, 1.499f), angRec, std::min(std::hypot(g.vx(), g.vy()), 89.9f),
                           eRadFrac, radReco, isRad, fracParent, angParent);
    }
  }

  template <soa::is_table TLegs, soa::is_table TCollisions, soa::is_table TPhotons, soa::is_table TMCParticles, typename TPreslice, typename TCut>
  void runSinglePhotonResponse(TCollisions const& collisions, TPhotons const& photons, TMCParticles const& mcParticles,
                               TPreslice const& perCollision, TCut const& cut)
  {
    for (const auto& collision : collisions) {
      if (!isSelectedCentrality(eventCentrality(collision))) {
        continue;
      }
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      auto photonsColl = photons.sliceBy(perCollision, collision.globalIndex());
      if (mQa.singlePhotonQa) {
        fillSinglePhotonResponse<TLegs>(collision, photonsColl, mcParticles, cut);
      }
    }
  }

  template <soa::is_table TMCEvents, soa::is_table TCollisions, soa::is_table TPhotons,
            soa::is_table TLegs, soa::is_table TMCParticles,
            typename TPresliceMCParts, typename TPresliceLegs, typename TCut>
  void runTruthEfficiency(TCollisions const& collisions,
                          TPhotons const& v0photons,
                          TLegs const& v0legs,
                          TMCParticles const& emmcParticles,
                          TPresliceMCParts const& perMCCollision,
                          TPresliceLegs const& perCollisionLegs,
                          TCut const& cut)
  {
    for (const auto& collision : collisions) {
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      const float centNote = eventCentrality(collision);
      if (!isSelectedCentrality(centNote)) {
        continue;
      }
      if (!collision.has_emmcevent()) {
        continue;
      }

      const int64_t thisCollisionId = collision.globalIndex();
      const int mcEventId = collision.template emmcevent_as<TMCEvents>().globalIndex();

      auto recoPhotonsColl = v0photons.sliceBy(perCollisionPCM, thisCollisionId);
      auto emmcPartsColl = emmcParticles.sliceBy(perMCCollision, mcEventId);
      auto legsColl = v0legs.sliceBy(perCollisionLegs, thisCollisionId);

      std::unordered_set<int> legIdsThisCollision;
      legIdsThisCollision.reserve(legsColl.size());
      for (const auto& leg : legsColl) {
        if (leg.has_emmcparticle()) {
          legIdsThisCollision.insert(leg.emmcparticleId());
        }
      }

      struct PhotonRecoInfo {
        bool hasV0 = false, passesCut = false;
        bool hasCounts = false;
        pairutil::V0PhotonLegCounts legCounts{};
      };
      std::unordered_map<int, PhotonRecoInfo> gammaRecoMap;
      gammaRecoMap.reserve(recoPhotonsColl.size());

      for (const auto& g : recoPhotonsColl) {
        const auto pos = g.template posTrack_as<TLegs>();
        const auto neg = g.template negTrack_as<TLegs>();
        if (pos.collisionId() != thisCollisionId || neg.collisionId() != thisCollisionId) {
          continue;
        }
        if (!pos.has_emmcparticle() || !neg.has_emmcparticle()) {
          continue;
        }
        const auto mcPos = pos.template emmcparticle_as<aod::EMMCParticles>();
        const auto mcNeg = neg.template emmcparticle_as<aod::EMMCParticles>();
        if (!mcPos.has_mothers() || !mcNeg.has_mothers()) {
          continue;
        }
        const int posMotherId = mcPos.mothersIds()[0], negMotherId = mcNeg.mothersIds()[0];
        if (posMotherId != negMotherId) {
          continue;
        }
        const int gammaId = posMotherId;
        if (emmcParticles.iteratorAt(gammaId).pdgCode() != kGamma) {
          continue;
        }
        const bool passes = cut.template IsSelected<std::decay_t<decltype(g)>, TLegs>(g);
        auto& info = gammaRecoMap[gammaId];
        info.hasV0 = true;
        info.passesCut = info.passesCut || passes;
        if (passes && !info.hasCounts) {
          info.legCounts = pairutil::getV0PhotonLegCounts(pos, neg);
          info.hasCounts = true;
        }
      }

      // ─── true photons that converted ──────────────────────────
      std::vector<TruthGamma> trueGammas;
      trueGammas.reserve(32);

      for (const auto& g : emmcPartsColl) {
        if (g.pdgCode() != kGamma) {
          continue;
        }
        if (!g.isPhysicalPrimary() && !g.producedByGenerator()) {
          continue;
        }
        if (std::fabs(g.eta()) > pcmcuts.cfgMaxEtaV0.value) {
          continue;
        }
        const float mcV0PtMin = (mctruth.cfgMCMinV0Pt.value > 0.f)
                                  ? mctruth.cfgMCMinV0Pt.value
                                  : pcmcuts.cfgMinPtV0.value;
        if (g.pt() < mcV0PtMin) {
          continue;
        }
        if (!g.has_daughters()) {
          continue;
        }

        int posId = -1, negId = -1;
        float rTrue = -1.f;
        for (const auto& dId : g.daughtersIds()) {
          if (dId < 0) {
            continue;
          }
          const auto d = emmcParticles.iteratorAt(dId);
          if (d.pdgCode() == kPositron) {
            posId = dId;
            rTrue = std::sqrt(d.vx() * d.vx() + d.vy() * d.vy());
          } else if (d.pdgCode() == kElectron) {
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
             static_cast<float>(mcNegE.pt()) < mctruth.cfgMCMinLegPt.value)) {
          continue;
        }

        TruthGamma tg;
        tg.id = static_cast<int>(g.globalIndex());
        tg.posId = posId;
        tg.negId = negId;
        tg.eta = static_cast<float>(g.eta());
        tg.phi = static_cast<float>(g.phi());
        tg.pt = static_cast<float>(g.pt());
        tg.rTrue = rTrue;
        tg.conv = {static_cast<float>(mcPosE.vx()), static_cast<float>(mcPosE.vy()), static_cast<float>(mcPosE.vz())};
        trueGammas.push_back(tg);
      }

      // ─── pairs of true photons: stages ────────────────────────
      for (size_t i = 0; i < trueGammas.size(); ++i) {
        for (size_t j = i + 1; j < trueGammas.size(); ++j) {
          const auto& g1 = trueGammas[i];
          const auto& g2 = trueGammas[j];
          const float px1 = g1.pt * std::cos(g1.phi), py1 = g1.pt * std::sin(g1.phi);
          const float px2 = g2.pt * std::cos(g2.phi), py2 = g2.pt * std::sin(g2.phi);
          const float kt = 0.5f * std::sqrt((px1 + px2) * (px1 + px2) + (py1 + py2) * (py1 + py2));

          if (mctruth.cfgMCMinKt > 0.f && kt < mctruth.cfgMCMinKt) {
            continue;
          }
          if (mctruth.cfgMCMaxKt > 0.f && kt > mctruth.cfgMCMaxKt) {
            continue;
          }

          const float e1 = g1.pt * std::cosh(g1.eta), e2 = g2.pt * std::cosh(g2.eta);
          const float dot = e1 * e2 - (px1 * px2 + py1 * py2 + g1.pt * std::sinh(g1.eta) * g2.pt * std::sinh(g2.eta));
          const float qinv_true = std::sqrt(std::max(0.f, 2.f * dot));

          if (mctruth.cfgMCMaxQinv > 0.f && qinv_true > mctruth.cfgMCMaxQinv) {
            continue;
          }

          auto it1 = gammaRecoMap.find(g1.id), it2 = gammaRecoMap.find(g2.id);
          const bool g1Built = (it1 != gammaRecoMap.end()) && it1->second.hasV0;
          const bool g2Built = (it2 != gammaRecoMap.end()) && it2->second.hasV0;
          const bool g1Sel = (it1 != gammaRecoMap.end()) && it1->second.passesCut;
          const bool g2Sel = (it2 != gammaRecoMap.end()) && it2->second.passesCut;
          const bool pairClassSel = !cfgDoPhotonClassPairCut.value || (g1Sel && g2Sel && pairutil::isPairPhotonClassSelected(it1->second.legCounts, it2->second.legCounts, mPhotonClassSelA, mPhotonClassSelB));

          const bool pairAll4LegsThisColl =
            legIdsThisCollision.contains(g1.posId) > 0 && legIdsThisCollision.contains(g1.negId) > 0 &&
            legIdsThisCollision.contains(g2.posId) > 0 && legIdsThisCollision.contains(g2.negId) > 0;

          if (g1.rTrue >= 0.f && g2.rTrue >= 0.f) {
            fRegistryTruthMC.fill(HIST("MC/TruthAO2D/hRconv1_vs_Rconv2_truthConverted"), g1.rTrue, g2.rTrue);
            if (g1Sel && g2Sel && pairClassSel) {
              fRegistryTruthMC.fill(HIST("MC/TruthAO2D/hRconv1_vs_Rconv2_bothPhotonsSelected"), g1.rTrue, g2.rTrue);
            }
          }
          if (qaflags.cfgFillNoteQA.value) {
            const float dConv = std::min(std::hypot(g1.conv[0] - g2.conv[0], g1.conv[1] - g2.conv[1], g1.conv[2] - g2.conv[2]), 999.f);
            const float ktNote = std::min(kt, 1.99f);
            fRegistryTruthMC.fill(HIST("MC/TruthAO2D/Note/hTruePairs"), 0.f, qinv_true, ktNote, centNote, dConv);
            if (pairAll4LegsThisColl) {
              fRegistryTruthMC.fill(HIST("MC/TruthAO2D/Note/hTruePairs"), 1.f, qinv_true, ktNote, centNote, dConv);
            }
            if (g1Built && g2Built) {
              fRegistryTruthMC.fill(HIST("MC/TruthAO2D/Note/hTruePairs"), 2.f, qinv_true, ktNote, centNote, dConv);
            }
            if (g1Sel && g2Sel && pairClassSel) {
              fRegistryTruthMC.fill(HIST("MC/TruthAO2D/Note/hTruePairs"), 3.f, qinv_true, ktNote, centNote, dConv);
            }
          }
        }
      }
    } // end collision loop
  } // end runTruthEfficiency

  void processPairs(aod::PCMMixCollisions const& collisions, MyV0PhotonsMix const& photons, aod::V0LegsMix const& legs, aod::PCMMixLegs const& mixLegs,
                    aod::V0PhotonsKFChi2Mix const& kfChi2)
  {

    std::vector<std::vector<int64_t>> photonsOf(collisions.size());
    std::vector<PcPhoton> pc(photons.size());
    for (const auto& photon : photons) {
      if (photon.collisionId() >= 0 && photon.collisionId() < static_cast<int>(collisions.size())) {
        photonsOf[photon.collisionId()].push_back(photon.globalIndex());
      }
      pc[photon.globalIndex()] = pcPhotonFrom(photon, legs, mixLegs, kfChi2);
    }
    mPhotonPass.assign(photons.size(), -1);
    for (const auto& collision : collisions) {
      initCCDB(collision); // magnetic field and drift velocity of the run (cached per run)
      const auto& ids = photonsOf[collision.globalIndex()];
      if (mDoCorrelation) {
        fillCorrelationEvent(collision, ids, photons, pc, mixLegs);
      }
      if (mDoPairClassQA) {
        fillPairClassQAEvent(collision, ids, photons, pc, mixLegs);
      }
    }
    ndf++;
  }
  PROCESS_SWITCH(Photonhbt, processPairs, "photon pairs of the artificial events of builder run 2 (PCMMixCollisions): correlation block (cfgDoCorrelation; same event, co-reconstructed or standard mixed event) and pair-class QA block (cfgDoPairClassQA)", false);

  void processGenerated(FilteredMyMCCollisions const& mccollisions,
                        MyV0Photons const& v0photons,
                        MyMCV0Legs const& v0legs,
                        aod::EMMCParticles const& mcParticles,
                        aod::EMMCEvents const& /*mcEvents*/)
  {
    runSinglePhotonResponse<MyMCV0Legs>(mccollisions, v0photons, mcParticles,
                                        perCollisionPCM, fV0PhotonCut);
    runTruthEfficiency<aod::EMMCEvents>(mccollisions, v0photons, v0legs, mcParticles,
                                        perMCCollisionEMMCParts, perCollisionV0Legs, fV0PhotonCut);
  }
  PROCESS_SWITCH(Photonhbt, processGenerated, "generated-level MC: true-pair stages (MC/TruthAO2D) and single-photon response (SinglePhotonMC/Response)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<Photonhbt>(context)};
}
