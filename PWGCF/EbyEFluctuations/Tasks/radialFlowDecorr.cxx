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

/// \file radialFlowDecorr.cxx
/// \brief Analysis task for event-by-event radial-flow decorrelation measurement.
/// \author Somadutta Bhatta

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DetectorsCommonDataFormats/AlignParam.h>
#include <FT0Base/Geometry.h>
#include <FV0Base/Geometry.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TArrayD.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <THnSparse.h>
#include <TList.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom3.h>
#include <TString.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

struct RadialFlowDecorr {

  // --- fixed constants ---------------------------------------------------------
  static constexpr int KIntM = 4; // pT-moment order used in the sums (need m up to 3 for c3)
  static constexpr int KIntK = 4; // weight-power order used in the sums (need k up to 3 for c3)

  static constexpr int KNEtaHalfBinsMax = 8;                // 0.8 / 0.1: finest (DATA) bins per side
  static constexpr int KNEtaMax = 2 * KNEtaHalfBinsMax + 1; // + index-0 full-range reference bin

  // Bootstrap: at most KMaxBoot replicas of 8 profiles each, in one registry of at most
  // 512 objects. nBoot = min(cfgNBootstrap, KMaxBoot) replicas are filled.
  static constexpr int KMaxBoot = 60;

  static constexpr float KFloatEpsilon = 1e-6f;
  static constexpr float KEtaEdgeTolerance = 1e-3f; // slack for cfgCutEta being an integer multiple of the bin width
  static constexpr float KBinOffset = 0.5f;
  static constexpr float KPhiMin = 0.f;
  static constexpr int KNbinsZvtx = 240;
  static constexpr float KZvtxMin = -12.f;
  static constexpr float KZvtxMax = 12.f;
  static constexpr float KPMin = 0.f;
  static constexpr float KPMax = 10.f;
  static constexpr int KNbinsPt = 200;
  static constexpr float KPtMin = 0.15f;
  static constexpr float KPtMax = 10.f;
  static constexpr float KEtaMin = -1.2f;
  static constexpr float KEtaMax = 1.2f;
  static constexpr int KNbinsPhi = 64;
  static constexpr int KNbinsPtRes = 50;
  static constexpr int KNbinsEtaRes = 100;
  static constexpr int KNbinsVz = 80;
  static constexpr float KVzMin = -40.f;
  static constexpr float KVzMax = 40.f;
  static constexpr int KNbinsEtaFine = 20;
  static constexpr float KEtaFineMax = 1.f;
  // occupancy QA axis
  static constexpr int KNbinsOccQA = 150;
  static constexpr float KOccQAMax = 15000.f;

  // Forward detectors: FT0 channels 0-95 are A side, 96-207 C side; FV0-A has 48 cells.
  static constexpr int KNChFT0A = 96;
  static constexpr int KNChFT0 = 208;
  static constexpr int KNChFV0 = 48;
  static constexpr int KNFwdMax = 32; // detector totals + eta sub-bins

  enum ECentralityEstimator {
    kCentFT0C = 1,
    kCentFT0M = 2,
    kCentFDDM = 3,
    kCentFV0A = 4
  };
  enum SystemType {
    kPbPb = 1,
    kNeNe = 2,
    kOO = 3,
    kpp = 4
  };
  enum EITSLayerCut {
    kITSLayerAuto = -1, // nominal: pp kIsGoodITSLayer0123, other systems kIsGoodITSLayersAll
    kITSLayerNone = 0,
    kITSLayer3 = 1,
    kITSLayer0123 = 2,
    kITSLayersAll = 3
  };
  enum EFwdDet {
    kFT0A = 0,
    kFT0C,
    kFV0A,
    kNFwdDet
  };
  // Histogram level: DATA uses level 0 only; MC fills truth, reco and efficiency-corrected reco.
  enum ELevel {
    kGen = 0,
    kReco,
    kRecoEffCorr,
    kNLevels
  };
  // Last axis of the multi-dimensional objects.
  enum EMeanObs { // Prof_Mean_Nch_etabin_obs
    kMeanPt = 0,
    kMeanN,
    kNMeanObs
  };
  enum EFlucObs { // Prof_Obs_Nch_etabin_obs
    kObsN = 0,    // N(eta), efficiency-corrected
    kObsMeanPt,   // [pT](eta)
    kObsC2,       // c2 within the bin
    kObsC3,       // c3 within the bin
    kObsC2Sub,    // <d[pT]_A d[pT]_C>, mirror pair, at the A bin
    kObsC3Sub,    // c2_A <d[pT]>_C at A (AAB) and c2_C <d[pT]>_A at C (ABB), mirror pair
    kNFlucObs
  };
  enum EEffCat {         // efficiency maps h3_* (centrality, pT, eta)
    kEffAllPrim = 0,     // charged primaries of the MC collision
    kEffRecoMatchedPrim, // reco tracks of primaries of the same MC collision, at the MC pT and eta
    kEffAllReco,
    kEffRecoSec,       // secondaries of the same MC collision
    kEffRecoOtherColl, // particles of another MC collision (pile-up, wrong association)
    kEffRecoFake,      // reco tracks without MC label
    kNEffCat
  };
  inline static const std::array<std::string, kNEffCat> effMapName = {
    "h3_AllPrimary", "h3_RecoMatchedToPrimary", "h3_AllReco", "h3_RecoUnMatchedToPrimary_Secondary",
    "h3_RecoUnMatchedToPrimary_OtherCollision", "h3_RecoUnMatchedToPrimary_Fake"};
  enum ETrackOrigin { // MCQA/Hist3D_Nch_pT_origin
    kOriginPrimary = 0,
    kOriginSecondary,
    kOriginOtherColl,
    kOriginNoLabel,
    kNOrigin
  };
  enum EVtxQA {     // MCQA/Hist2D_dVz_type
    kVtxSingle = 0, // the only reco collision of the MC collision
    kVtxSplitBest,  // the one used among several
    kVtxSplitOther, // the others (skipped)
    kNVtxQA
  };
  enum EFlatType { // hFlatMaps
    kFlatRaw = 0,
    kFlatEffWtd, // (1 - fake) / eff (Base only)
    kFlatWtd,    // all weights of the pass
    kNFlatType
  };
  enum EChanQA { // FwdQA/Prof_Channel_qa
    kChAmpRaw = 0,
    kChAmpCorr,
    kChEta,
    kNChQA
  };
  // Event cut flow (hEvtCount, x axis); the last two bins count selected events only
  enum EEvtCut {
    kEvAll = 0,
    kEvSel8,
    kEvVtxZ,
    kEvGoodZvtx,
    kEvTFBorder,
    kEvITSROFBorder,
    kEvVertexITSTPC,
    kEvSameBunch,
    kEvTimeRange,
    kEvITSLayer,
    kEvAddPileup,
    kEvVariationCut,
    kEvCent,
    kEvInfoPrevRof,    // selected events passing kNoHighMultCollInPrevRof
    kEvInfoOccUnknown, // selected events with occupancy < 0 (time-frame border)
    kNEvCut
  };
  // Track cut flow (hTrkCount, x axis)
  enum ETrkCut {
    kTrkAll = 0,
    kTrkSign,
    kTrkHasColl,
    kTrkPVContrib,
    kTrkPtEta,
    kTrkGlobal,
    kTrkITS,
    kTrkTPC,
    kTrkTPCChi2,
    kTrkDCAz,
    kTrkVariationCut,
    kTrkWeight, // valid efficiency (and flattening) weight: used
    kNTrkCut
  };

  // Cut QA (processCutQA): each quantity the Base selection cuts on, before and after it.
  enum ECutQAStage {
    kQABefore = 0, // events: all collisions; tracks: all tracks of the selected events
    kQAAfter,      // Base-selected events / tracks
    kNQAStage
  };
  enum ECutQAEvt1D {
    kQAVtxZ = 0,
    kQACent,
    kQAOccupancy,
    kNQAEvt1D
  };
  enum ECutQAEvt2D {
    kQANchNPV = 0,
    kQANchCent,
    kNQAEvt2D
  };
  enum ECutQATrk1D {
    kQAPt = 0,
    kQAEta,
    kQAItsNCls,
    kQAItsChi2,
    kQATpcNClsFound,
    kQATpcCrossedRows,
    kQATpcCrossedRowsOverFindable,
    kQATpcChi2,
    kQAPVContributor,
    kQAGlobalTrack,
    kNQATrk1D
  };
  enum ECutQATrk2D {
    kQADcaXYPt = 0,
    kQADcaZPt,
    kNQATrk2D
  };

  // Systematic variations, all evaluated in every job. Base is the nominal measurement.
  // The cut values are the K* constants below; only the occupancy limits are configurable.
  //   Cent      : FT0M instead of the nominal estimator (cfgCentralityChoice, FT0C)
  //   Chi2TPC   : chi2 per TPC cluster < 2 instead of cfgCutTpcChi2NCl
  //   DCAxy     : |DCAxy| < 5/7 x (globalTrack pT-dependent DCAxy cut), i.e. 5 instead of 7 sigma
  //   DCAz      : |DCAz| < 1 cm instead of cfgCutTracKDcaMaxZ
  //               DCAxy and DCAz are null checks: with the nominal isPVContributor requirement even
  //               5 sigma DCAxy or |DCAz| < 0.05 cm removes only ~0.2% of the selected tracks
  //   Eff     : efficiency scaled by -10% at all pT (cancels in [pT]: a check)
  //   EffpTDep  : efficiency scaled by -10% at cfgPtMin, linearly to 0 at 10 GeV/c
  //   Flat      : no acceptance (flattening) weights
  //   ITSLayer  : ITS good-layer requirement one step looser (All -> 0123, 0123 -> Layer3)
  //   NITS      : itsNCls >= 4 instead of cfgITScluster
  //   NTPC      : tpcNClsFound >= 80 instead of cfgTPCcluster
  //   Pileup    : additionally kNoCollInRofStandard
  //   Vz        : |z_vtx| < 8 cm instead of cfgVtxZCut
  //   NchBin2   : Base selection and inputs, observables (and its own mean pass) in N_ch bins
  //               of 2 (PbPb, OO) or 1 (pp) tracks: cross-check of the N_ch binning of the mean subtraction
  //   Occupancy : 0 <= trackOccupancyInTimeRange < cfgSystOccMax* of the system
  //   PileupNSigma : n-sigma pile-up bands narrower by 1 (mean and fluctuation passes)
  // Since the observables are stored versus N_ch, Cent differs from Base only in the events
  // kept (cent < cfgCentMax) and in its (centrality, N_ch) map, which defines its classes.
  enum ESystType {
    kSystBase = 0,
    kSystCent,
    kSystChi2TPC,
    kSystDCAxy,
    kSystDCAz,
    kSystEff,
    kSystEffPtDep,
    kSystFlat,
    kSystITSLayer,
    kSystNITS,
    kSystNTPC,
    kSystPileup,
    kSystVz,
    kSystNchBin,
    kSystOccupancy,
    kSystPileupNSigma,
    kNSyst
  };
  inline static const std::array<std::string, kNSyst> systName = {
    "Base", "systCent", "systChi2TPC", "systDCAxy", "systDCAz", "systEff", "systEffpTDep",
    "systFlat", "systITSLayer", "systNITS", "systNTPC", "systPileup", "systVz", "systNchBin2",
    "systOccupancy", "systPileupNSigma"};
  static constexpr uint32_t KAllSystMask = (1 << kNSyst) - 1;
  static constexpr int KNSystPerRegistry = 4; // variations 1..15 are spread over 4 registries
  // Cut of each variation (fixed; the occupancy limits are configurables, cfgSystOccMax*).
  // Pileup adds kNoCollInRofStandard; ITSLayer loosens the nominal ITS good-layer requirement by one step.
  static constexpr int KSystCentralityChoice = 2;     // Cent: FT0M
  static constexpr float KSystChi2TPC = 2.f;          // Chi2TPC: chi2 per TPC cluster < 2 (nominal 2.5)
  static constexpr float KSystDcaXyScale = 5.f / 7.f; // DCAxy: 5 sigma (nominal 7 sigma)
  static constexpr float KSystDcaMaxZ = 1.f;          // DCAz: |DCAz| < 1 cm (nominal 2 cm)
  static constexpr float KEffVarPercent = -10.f;      // Eff, EffpTDep: efficiency scaled by -10%
  static constexpr float KEffVarPtZero = 10.f;        // EffpTDep: variation falls linearly to 0 at 10 GeV/c
  static constexpr int KSystITScluster = 4;           // NITS: itsNCls >= 4 (nominal 5)
  static constexpr int KSystTPCcluster = 80;          // NTPC: tpcNClsFound >= 80 (nominal 70)
  static constexpr float KSystVtxZCut = 8.f;          // Vz: |z_vtx| < 8 cm (nominal 10 cm)
  static constexpr int KSystNchBinWidthPbPb = 2;      // NchBin2: N_ch bin width, PbPb (nominal 10)
  static constexpr int KSystNchBinWidthOO = 2;        // NchBin2: OO and NeNe (nominal 5)
  static constexpr int KSystNchBinWidthPP = 1;        // NchBin2: pp (nominal 2)
  static constexpr float KSystPupnSigStep = 1.f;      // PileupNSigma: n-sigma bands narrower by 1
  static constexpr double KMinBandEntries = 50.;      // n-sigma band bins with fewer flattening-pass events apply no band cut

  // Output folder of a variation (Base keeps the top level).
  static std::string systPrefix(int v) { return v == kSystBase ? std::string{} : systName[v] + "/"; }
  // Suffix appended to every CCDB object path of a variation.
  static std::string systSuffix(int v) { return v == kSystBase ? std::string{} : "_" + systName[v]; }
  // Variation whose efficiency map v uses. Cent, Eff, EffpTDep, Flat, NchBin2 and PileupNSigma
  // have the Base selection in the efficiency pass and use the Base map; Eff and EffpTDep scale
  // it when applied. (The n-sigma bands are applied only in the mean and fluctuation passes.)
  static constexpr int effSource(int v)
  {
    return (v == kSystCent || v == kSystEff || v == kSystEffPtDep || v == kSystFlat || v == kSystNchBin || v == kSystPileupNSigma) ? static_cast<int>(kSystBase) : v;
  }
  // Variation whose flattening map v uses. Cent, NchBin2 and PileupNSigma use the Base map;
  // Flat uses none (-1).
  static constexpr int flatSource(int v)
  {
    if (v == kSystFlat) {
      return -1;
    }
    return (v == kSystCent || v == kSystNchBin || v == kSystPileupNSigma) ? static_cast<int>(kSystBase) : v;
  }

  static constexpr float KinvalidCentrality = -1.0f;

  // Observable eta binning, filled at init(): index 0 = full-range reference bin.
  std::vector<float> etaLw;
  std::vector<float> etaUp;
  int nEta = 9;        // active eta-bin count (set in init())
  int nBoot = 0;       // active bootstrap samples (set in init())
  bool doBoot = false; // bootstrap active for this run (base data fluc only)

  // Variations filled by each pass (set in init()): maskEff in the MC efficiency pass,
  // maskFlat in the flattening passes, maskAll in the mean and fluctuation passes.
  uint32_t maskAll = 1;
  uint32_t maskEff = 1;
  uint32_t maskFlat = 1;
  int itsLayerNom = kITSLayersAll;
  int itsLayerVar = kITSLayer0123;
  int occupancyMax = 0;    // Occupancy variation cut of this system
  float pupNSigma = 5.f;   // width of the n-sigma pile-up bands of this system
  int nchWidth = 10;       // N_ch bin width of the observables of this system
  int nchWidthNchBin2 = 2; // N_ch bin width of NchBin2 for this system

  // Forward bins: 0-2 detector totals (FT0A, FT0C, FV0A), then the eta sub-bins of each detector.
  int nFwd = kNFwdDet;
  std::array<std::vector<float>, kNFwdDet> fwdEdges;
  std::array<int, kNFwdDet> fwdOffset{};
  std::vector<std::string> fwdLabels;

  // --- configurables -----------------------------------------------------------
  // Event selection (nominal)
  struct : ConfigurableGroup {
    Configurable<float> cfgVtxZCut{"cfgVtxZCut", 10.f, "|z_{vtx}| acceptance (cm)"};
    Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 1, "Which centrality estimator? 1-->FT0C, 2-->FT0M, 3-->FDDM, 4-->FV0A"};
    Configurable<float> cfgCentMax{"cfgCentMax", 99.f, "Events with centrality below this value (%) are kept"};
    Configurable<bool> cfgIsGoodZvtxFT0VsPV{"cfgIsGoodZvtxFT0VsPV", true, "Good Vertexing cut"};
    Configurable<bool> cfgEvSelNoTimeFrameBorder{"cfgEvSelNoTimeFrameBorder", true, "kNoTimeFrameBorder"};
    Configurable<bool> cfgEvSelNoITSROFrameBorder{"cfgEvSelNoITSROFrameBorder", true, "kNoITSROFrameBorder"};
    Configurable<bool> cfgEvSelIsVertexITSTPC{"cfgEvSelIsVertexITSTPC", true, "kIsVertexITSTPC"};
    Configurable<bool> cfgEvSelNoSameBunchPileup{"cfgEvSelNoSameBunchPileup", true, "Pileup removal: kNoSameBunchPileup"};
    Configurable<bool> cfgEvSelNoCollInTimeRangeStd{"cfgEvSelNoCollInTimeRangeStd", true, "Pileup removal: kNoCollInTimeRangeStandard"};
    Configurable<int> cfgITSLayerCut{"cfgITSLayerCut", -1, "ITS good-layer cut: -1 auto (pp 0123, else All), 0 none, 1 Layer3, 2 Layer0123, 3 LayersAll"};
  } cfgEvSel;

  // Pile-up cuts on the correlations of N_ch (global tracks) with N_PV and centrality.
  // The linear cut is applied in every pass. The flattening passes (MC and DATA) fill
  // Hist2D_globalTracks_cent and Hist2D_globalTracks_PVTracks after it; the mean +- n sigma
  // bands of N_ch from these maps are then applied in the mean and fluctuation passes.
  struct : ConfigurableGroup {
    Configurable<float> cfgPupnSigPbPb{"cfgPupnSigPbPb", 5.f, "Width of the n-sigma bands (number of RMS), PbPb"};
    Configurable<float> cfgPupnSigOO{"cfgPupnSigOO", 5.f, "Width of the n-sigma bands (number of RMS), OO and NeNe"};
    Configurable<float> cfgPupnSigPP{"cfgPupnSigPP", 6.f, "Width of the n-sigma bands (number of RMS), pp"};
    Configurable<bool> cfgApplySigPupCut{"cfgApplySigPupCut", 0, "Mean/Fluc: n-sigma band of N_ch versus centrality"};
    Configurable<bool> cfgApplySigPupCutPV{"cfgApplySigPupCutPV", 0, "Mean/Fluc: n-sigma band of N_ch versus N_PV"};
    Configurable<bool> cfgApplyLinPupCut{"cfgApplyLinPupCut", 0, "Linear cut of N_ch versus N_PV (every pass)"};
    Configurable<float> cfgLinPupParam0{"cfgLinPupParam0", 3.0f, "(Upper) Linear Pileup Cut Const"};
    Configurable<float> cfgLinPupParam1{"cfgLinPupParam1", 3.0f, "(Upper) Linear Pileup Slope"};
    Configurable<float> cfgLinPupParam2{"cfgLinPupParam2", 3.0f, "(Lower) Linear Pileup Cut Const"};
    Configurable<float> cfgLinPupParam3{"cfgLinPupParam3", 3.0f, "(Lower) Linear Pileup Slope"};
  } cfgPileup;

  // Track selection (nominal), on top of isGlobalTrack (Run 3 globalTrack, incl. its 7 sigma pT-dependent DCAxy cut)
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "min pT for observables"};
  Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "max pT for observables"};
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut (track selection)"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 10.0f, "Higher pT cut (track selection)"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute Eta cut"};
  Configurable<float> cfgEtaBinWidth{"cfgEtaBinWidth", 0.1f, "DATA only: width of each narrow observable eta bin (eta units); cfgCutEta must be an integer multiple so eta=0 is a bin edge. MC is pinned to 0.2."};
  Configurable<float> cfgCutTracKDcaMaxZ{"cfgCutTracKDcaMaxZ", 2.0f, "Maximum |DCAz| (cm)"};
  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Nominal: itsNCls >= value"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Nominal: tpcNClsFound >= value"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Nominal: tpcChi2NCl < value"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};

  // Systematic variations (names stay flat in the JSON)
  struct : ConfigurableGroup {
    Configurable<int> cfgSystMask{"cfgSystMask", static_cast<int>(KAllSystMask), "Bitmask of active variations, bit i = ESystType i: 0 Base (always on), 1 Cent, 2 Chi2TPC, 3 DCAxy, 4 DCAz, 5 Eff, 6 EffpTDep, 7 Flat, 8 ITSLayer, 9 NITS, 10 NTPC, 11 Pileup, 12 Vz, 13 NchBin2, 14 Occupancy, 15 PileupNSigma"};
    Configurable<float> cfgDcaXyP0{"cfgDcaXyP0", 0.0105f, "globalTrack DCAxy cut p0 + p1 / pT^p2 (trackselection dcaSetup=0; dcaSetup=1: 0.004, 0.013, 1)"};
    Configurable<float> cfgDcaXyP1{"cfgDcaXyP1", 0.035f, "globalTrack DCAxy cut p1"};
    Configurable<float> cfgDcaXyP2{"cfgDcaXyP2", 1.1f, "globalTrack DCAxy cut p2"};
    Configurable<int> cfgSystOccMaxPbPb{"cfgSystOccMaxPbPb", 2000, "Occupancy: maximum trackOccupancyInTimeRange, PbPb"};
    Configurable<int> cfgSystOccMaxOO{"cfgSystOccMaxOO", 1000, "Occupancy: maximum trackOccupancyInTimeRange, OO and NeNe"};
    Configurable<int> cfgSystOccMaxPP{"cfgSystOccMaxPP", 500, "Occupancy: maximum trackOccupancyInTimeRange, pp"};
  } cfgSyst;

  // N_ch axes. Largest N_ch (global tracks) seen in the data: 3542 (PbPb), 328 (OO), 72 (pp)
  struct : ConfigurableGroup {
    Configurable<int> cfgNchPbMax{"cfgNchPbMax", 4000, "Max N_ch (global tracks) and N_PV of the axes, PbPb"};
    Configurable<int> cfgNchOMax{"cfgNchOMax", 400, "Max N_ch (global tracks) and N_PV of the axes, OO and NeNe"};
    Configurable<int> cfgNchPPMax{"cfgNchPPMax", 100, "Max N_ch (global tracks) and N_PV of the axes, pp"};
    Configurable<int> cfgNchBinWidthPbPb{"cfgNchBinWidthPbPb", 10, "N_ch bin width (tracks) of the observables, PbPb; Mean and Fluc passes must use the same value"};
    Configurable<int> cfgNchBinWidthOO{"cfgNchBinWidthOO", 5, "N_ch bin width (tracks) of the observables, OO and NeNe"};
    Configurable<int> cfgNchBinWidthPP{"cfgNchBinWidthPP", 2, "N_ch bin width (tracks) of the observables, pp"};
  } cfgNch;

  Configurable<int> cfgSys{"cfgSys", 1, "Which collision system? 1-->PbPb, 2-->NeNe, 3-->OO, 4-->pp"};
  Configurable<int> cfgNBootstrap{"cfgNBootstrap", 20, "Number of Poisson bootstrap samples (Base, DataFluc only)"};
  Configurable<int> cfgBootstrapSeed{"cfgBootstrapSeed", 0, "TRandom3 seed for bootstrap (0 = machine-random per job)"};

  Configurable<bool> cfgFlat{"cfgFlat", false, "Whether to use flattening weights"};
  Configurable<bool> cfgEff{"cfgEff", false, "Whether to use Efficiency weights"};

  // Forward multiplicity (FT0A, FT0C, FV0A amplitudes) correlated with mid-rapidity [pT] and N
  struct : ConfigurableGroup {
    Configurable<bool> cfgFwd{"cfgFwd", true, "Mean/Fluc: forward amplitudes and their covariance with [pT](eta) and N(eta)"};
    Configurable<bool> cfgFwdGainEq{"cfgFwdGainEq", true, "DATA: apply the FT0/FV0 gain equalisation from cfgGainEqPath"};
    Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Analysis/EventPlane/GainEq", "CCDB path of the FT0/FV0 gain equalisation (+ /FT0, /FV0); central objects per run for LHC23zzh Pb-Pb and 2025 O-O/Ne-Ne"};
    Configurable<bool> cfgFwdVzCorr{"cfgFwdVzCorr", true, "Channel eta relative to the primary vertex (else to z = 0)"};
    Configurable<std::vector<float>> cfgFwdEtaEdgesFT0A{"cfgFwdEtaEdgesFT0A", {3.5f, 4.2f, 4.9f}, "FT0A eta sub-bin edges"};
    Configurable<std::vector<float>> cfgFwdEtaEdgesFT0C{"cfgFwdEtaEdgesFT0C", {-3.3f, -2.7f, -2.1f}, "FT0C eta sub-bin edges"};
    Configurable<std::vector<float>> cfgFwdEtaEdgesFV0A{"cfgFwdEtaEdgesFV0A", {2.2f, 2.7f, 3.2f, 3.9f, 4.5f, 5.1f}, "FV0A eta sub-bin edges"};
  } cfgFwdOpt;

  Configurable<std::string> cfgCCDBurl{"cfgCCDBurl", "https://alice-ccdb.cern.ch", "ccdb url"};
  Configurable<std::string> cfgCCDBUserPath{"cfgCCDBUserPath", "/Users/s/somadutt", "Base CCDB path"};

  // --- axes --------------------------------------------------------------------
  AxisSpec centAxis1Per{100, 0.0, 100.0, "Centrality (%)"};
  AxisSpec centAxisBoot{20, 0.0, 100.0, "Centrality (%)"};                                             // bootstrap replicas: final 5% classes
  AxisSpec multQAAxis{1, 0., 1., "Nch", "Nch"};                                                        // QA (N_PV correlations)
  AxisSpec nchAxisObs{1, 0., 1., "N_{ch}"};                                                            // every observable (bins of nchWidth)
  AxisSpec nchAxisNchBin2{1, 0., 1., "N_{ch}"};                                                        // NchBin2 (bins of nchWidthNchBin2)
  AxisSpec centAxisEff{{0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "Centrality (%)"}; // efficiency maps

  AxisSpec vzAxis{5, -12.5, 12.5, "Vz"};
  AxisSpec chgAxis{3, -1.5, 1.5};
  AxisSpec pTAxis{{0.0, 0.2, 0.4, 0.6, 0.8, 1, 3, 5, 7, 10}, "pT Axis"};
  AxisSpec phiAxis{KNbinsPhi, KPhiMin, TwoPI, "#phi"};

  // etaFlatAxis is the *fixed* granularity used only for the flattening map, so
  // that a flattening map is reusable across the 0.1/0.2 observable binnings.
  AxisSpec etaFlatAxis{{-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8}, "#eta"};

  // etaAxis (physical) and etaBinAxis (integer index) follow the observable
  // binning and are rebuilt in init(). Placeholders here.
  AxisSpec etaAxis{9, -0.9, 0.9, "#eta"};
  AxisSpec etaBinAxis{10, -0.5, 9.5, "#eta bin Number"};
  AxisSpec gapAxis{15, -1.5, 1.5, "#Delta#eta"};
  AxisSpec sumAxis{15, -1.5, 1.5, "#Sigma#eta"};
  AxisSpec fwdAxis{kNFwdDet, -0.5, static_cast<double>(kNFwdDet) - 0.5, "forward bin"};

  // --- process switches --------------------------------------------------------
  Configurable<bool> cfgRunGetEff{"cfgRunGetEff", false, "Run MC pass to build efficiency/fake maps"};
  Configurable<bool> cfgRunGetMCFlat{"cfgRunGetMCFlat", false, "Run MC to get flattening weights"};
  Configurable<bool> cfgRunMCMean{"cfgRunMCMean", false, "Run MC mean(pT)"};
  Configurable<bool> cfgRunMCFluc{"cfgRunMCFluc", false, "Run MC fluctuations (C2, subevent)"};

  Configurable<bool> cfgRunGetDataFlat{"cfgRunGetDataFlat", false, "Run data get flattening weights"};
  Configurable<bool> cfgRunDataMean{"cfgRunDataMean", false, "Run DATA mean(pT)"};
  Configurable<bool> cfgRunDataFluc{"cfgRunDataFluc", false, "Run DATA fluctuations (C2, subevent)"};

  Service<ccdb::BasicCCDBManager> ccdb{};
  Service<o2::framework::O2DatabasePDG> pdg{};

  // histos: Base and QA; histosBoot: bootstrap replicas; histosSystA-D: the variations,
  // 4 per registry (a registry holds at most 512 objects).
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosBoot{"HistosBoot", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosSystA{"HistosSystA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosSystB{"HistosSystB", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosSystC{"HistosSystC", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosSystD{"HistosSystD", {}, OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry& reg(int v)
  {
    if (v == kSystBase) {
      return histos;
    }
    switch ((v - 1) / KNSystPerRegistry) {
      case 0:
        return histosSystA;
      case 1:
        return histosSystB;
      case 2:
        return histosSystC;
      default:
        return histosSystD;
    }
  }

  const AxisSpec& nchAxisOf(int v) const { return v == kSystNchBin ? nchAxisNchBin2 : nchAxisObs; }
  int nchBinWidthOf(int v) const { return v == kSystNchBin ? nchWidthNchBin2 : nchWidth; }

  TRandom3 rng; // bootstrap Poisson weights

  // --- per-variation inputs ----------------------------------------------------
  struct SystInputs {
    // efficiency and fake rate versus (centrality, pT, eta), and integrated over centrality (one
    // centrality bin: used where the cell is empty or not a valid efficiency)
    TH3* hEff = nullptr;
    TH3* hFake = nullptr;
    TH3* hEffPtEta = nullptr;
    TH3* hFakePtEta = nullptr;
    THnSparseF* hFlat = nullptr;
    // mean pass, per level: <[pT]> and <N> vs (N_ch, eta bin), <a> vs (N_ch, forward bin)
    std::array<TProfile3D*, kNLevels> pmean{};
    std::array<TProfile2D*, kNLevels> pmeanFwd{};
  };
  std::array<SystInputs, kNSyst> inputs{};

  // mean +- n sigma of N_ch per bin of x (centrality or N_PV), from the flattening pass
  struct NSigmaBand {
    std::vector<std::pair<float, float>> meanRms; // per bin of x; rms < 0 for an empty bin
    float xMin = 0.f, xMax = 0.f;
    bool applied = false;
  };
  enum EPileupBand {
    kBandCent = 0,
    kBandPV,
    kNBands
  };
  std::array<NSigmaBand, kNBands> pileupBands{};

  // --- per-variation histograms, filled through pointers (HIST() needs names known at compile time) ---
  struct EffHists {
    std::array<std::shared_ptr<TH3>, kNEffCat> maps{}; // (centrality, pT, eta), one per EEffCat
  };
  struct FlatHists {
    std::shared_ptr<THnSparse> maps; // (vz, sign, pT, eta, phi, EFlatType)
  };
  struct MeanHists {
    std::shared_ptr<TProfile3D> mean;    // (N_ch, eta bin, EMeanObs)
    std::shared_ptr<TProfile2D> fwdMean; // (N_ch, forward bin): <a>
  };
  struct FlucHists {
    std::shared_ptr<TProfile3D> obs;        // (N_ch, eta bin, EFlucObs)
    std::shared_ptr<TProfile3D> c2Sub2D;    // (N_ch, eta_A, eta_C)
    std::shared_ptr<TProfile3D> c3Sub2D;    // (N_ch, eta_A pair, eta_C single)
    std::shared_ptr<TProfile3D> gapSum2D;   // (N_ch, Delta eta, Sigma eta) of the C2 map
    std::shared_ptr<TProfile3D> c3GapSum2D; // (N_ch, eta_pair - eta_single, eta_pair + eta_single)
    std::shared_ptr<TProfile2D> fwdVar;     // (N_ch, forward bin): <da^2>
    std::shared_ptr<TProfile3D> covFwdPt;   // (N_ch, eta bin, forward bin): <d[pT] da>
    std::shared_ptr<TProfile3D> covFwdN;    // (N_ch, eta bin, forward bin): <dN da>
  };
  struct {
    std::shared_ptr<TH2> evtCount;                         // (EEvtCut, variation)
    std::shared_ptr<TH2> trkCount;                         // (ETrkCut, variation)
    std::array<std::shared_ptr<TH1>, kNSyst> evtCount1D{}; // EEvtCut, one per variation
    std::array<std::shared_ptr<TH1>, kNSyst> trkCount1D{}; // ETrkCut, one per variation
    std::shared_ptr<TH2> nchOccupancy;                     // (N_ch, occupancy) of the Base events
    std::array<std::array<std::shared_ptr<TH1>, kNQAStage>, kNQAEvt1D> qaEvt1D{};
    std::array<std::array<std::shared_ptr<TH2>, kNQAStage>, kNQAEvt2D> qaEvt2D{};
    std::array<std::array<std::shared_ptr<TH1>, kNQAStage>, kNQATrk1D> qaTrk1D{};
    std::array<std::array<std::shared_ptr<TH2>, kNQAStage>, kNQATrk2D> qaTrk2D{};
    std::array<std::shared_ptr<TH2>, kNSyst> nchCent{}; // (centrality, N_ch) per variation
    std::array<EffHists, kNSyst> eff{};
    std::array<FlatHists, kNSyst> flat{};
    std::array<std::array<MeanHists, kNLevels>, kNSyst> mean{};
    std::array<std::array<FlucHists, kNLevels>, kNSyst> fluc{};
    std::array<FlucHists, KMaxBoot> boot{};
    std::shared_ptr<TProfile2D> fwdChannel; // (channel: FT0 0-207, FV0 208-255; EChanQA)
    std::shared_ptr<TProfile> fwdEta;       // amplitude-weighted eta per forward bin
    std::shared_ptr<TH3> mcTrackOrigin;     // (N_ch, pT, ETrackOrigin) of the Base tracks
    std::shared_ptr<TH2> mcVtx;             // (z_reco - z_true, EVtxQA)
  } hists;

  // --- forward detector calibration (per run) ----------------------------------
  struct {
    o2::ft0::Geometry ft0Det;
    o2::fv0::Geometry* fv0Det = nullptr;
    int run = -1;
    std::array<std::array<float, 3>, KNChFT0> posFT0{};
    std::array<std::array<float, 3>, KNChFV0> posFV0{};
    std::array<float, KNChFT0> gainFT0{};
    std::array<float, KNChFV0> gainFV0{};
  } fwdCalib;

  // --- per-event work arrays ---------------------------------------------------
  // Power sums S_{m,k} = sum w^k pT^m per eta bin.
  struct PowerSums {
    std::array<std::array<std::array<double, KIntK>, KIntM>, KNEtaMax> pmwk{};
    std::array<std::array<double, KIntK>, KNEtaMax> wk{};
    void reset()
    {
      for (auto& a : pmwk) {
        for (auto& b : a) {
          b.fill(0.);
        }
      }
      for (auto& a : wk) {
        a.fill(0.);
      }
    }
    void add(int ieta, double pt, double w)
    {
      double wPow = 1.;
      for (int k = 0; k < KIntK; ++k) {
        double ptPow = 1.;
        for (int m = 0; m < KIntM; ++m) {
          pmwk[ieta][m][k] += wPow * ptPow;
          ptPow *= pt;
        }
        wk[ieta][k] += wPow;
        wPow *= w;
      }
    }
  };
  std::array<std::array<PowerSums, kNLevels>, kNSyst> sums{};
  // track cut flow of the current event, added to hTrkCount once per event
  std::array<std::array<double, kNTrkCut>, kNSyst> trkCounts{};

  struct FwdSignal {
    std::array<double, KNFwdMax> ampl{}; // amplitude (data, MC reco) or number of particles (MC truth)
    std::array<bool, KNFwdMax> valid{};  // detector read out in this event
  };

  struct EventSelection {
    uint32_t mask = 0;
    std::array<float, kNSyst> cent{};
    int occupancy = -1;
  };

  // Per-event observables of one (variation, level).
  struct EventObs {
    float nCh = 0.f;
    std::array<double, KNEtaMax> nW{}, mean{}, c2{}, c3{}, p1kBar{}, dN{};
    std::array<double, KNFwdMax> fwd{}, dFwd{};
    std::array<bool, KNFwdMax> fwdValid{};
  };

  // ===========================================================================
  // Selection helpers
  // ===========================================================================
  template <typename T>
  static bool passItsLayerCut(int cut, const T& col)
  {
    switch (cut) {
      case kITSLayer3:
        return col.selection_bit(o2::aod::evsel::kIsGoodITSLayer3);
      case kITSLayer0123:
        return col.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123);
      case kITSLayersAll:
        return col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);
      default:
        return true;
    }
  }

  // N_ch within mean +- nSigma RMS at x (events outside the x range of the band are rejected,
  // empty bins of the band accept every event).
  static bool insideBand(const NSigmaBand& band, float x, float nCh, float nSigma)
  {
    const auto& mr = band.meanRms;
    int bin = 1 + static_cast<int>((x - band.xMin) / (band.xMax - band.xMin) * (mr.size() - 2));
    if (bin < 1 || bin >= static_cast<int>(mr.size() - 1)) {
      return false;
    }
    if (mr[bin].second < 0.f) {
      return true;
    }
    return std::abs(nCh - mr[bin].first) <= nSigma * mr[bin].second;
  }

  // Linear cut and n-sigma bands of width nSigma.
  bool isPassAddPileup(float multPV, int trksize, float cent, float nSigma)
  {
    if (pileupBands[kBandCent].applied && !insideBand(pileupBands[kBandCent], cent, trksize, nSigma)) {
      return false;
    }
    if (pileupBands[kBandPV].applied && !insideBand(pileupBands[kBandPV], multPV, trksize, nSigma)) {
      return false;
    }
    if (cfgPileup.cfgApplyLinPupCut) {
      if (trksize > (cfgPileup.cfgLinPupParam0 + cfgPileup.cfgLinPupParam1 * multPV)) {
        return false;
      }
      if (trksize < (cfgPileup.cfgLinPupParam2 + cfgPileup.cfgLinPupParam3 * multPV)) {
        return false;
      }
    }
    return true;
  }

  // One event at cut-flow stage s of variation v: hEvtCountSyst and that variation's hEvtCount.
  void countEvent(int s, int v)
  {
    hists.evtCount->Fill(s, v);
    if (hists.evtCount1D[v]) {
      hists.evtCount1D[v]->Fill(s);
    }
  }

  // Event selection of all active variations: returns the variations the event passes and its
  // centrality for each of them. The cuts common to all variations come first. The event cut
  // flow of every active variation is filled unless count is false (cut QA).
  template <typename T>
  EventSelection isEventSelected(const T& col, int nCh, uint32_t active, bool count = true)
  {
    EventSelection ev;
    auto countAll = [&](int s) {
      for (int v = 0; v < kNSyst && count; ++v) {
        if (TESTBIT(active, v)) {
          countEvent(s, v);
        }
      }
    };
    countAll(kEvAll);
    if (!col.sel8()) {
      return ev;
    }
    countAll(kEvSel8);
    if (std::abs(col.posZ()) > cfgEvSel.cfgVtxZCut) {
      return ev;
    }
    countAll(kEvVtxZ);
    if (cfgEvSel.cfgIsGoodZvtxFT0VsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return ev;
    }
    countAll(kEvGoodZvtx);
    if (cfgEvSel.cfgEvSelNoTimeFrameBorder && !col.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return ev;
    }
    countAll(kEvTFBorder);
    if (cfgEvSel.cfgEvSelNoITSROFrameBorder && !col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return ev;
    }
    countAll(kEvITSROFBorder);
    if (cfgEvSel.cfgEvSelIsVertexITSTPC && !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return ev;
    }
    countAll(kEvVertexITSTPC);

    const float centBase = getCentrality(col, cfgEvSel.cfgCentralityChoice);
    const float multPV = col.multNTracksPV();
    ev.occupancy = col.trackOccupancyInTimeRange();
    const bool passAddPileup = isPassAddPileup(multPV, nCh, centBase, pupNSigma);
    const bool passAddPileupTight = isPassAddPileup(multPV, nCh, centBase, pupNSigma - KSystPupnSigStep);

    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(active, v)) {
        continue;
      }
      auto countOne = [&](int s) {
        if (count) {
          countEvent(s, v);
        }
      };
      const float cent = (v == kSystCent) ? getCentrality(col, KSystCentralityChoice) : centBase;
      ev.cent[v] = cent;
      if (cfgEvSel.cfgEvSelNoSameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      countOne(kEvSameBunch);
      if (cfgEvSel.cfgEvSelNoCollInTimeRangeStd && !col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        continue;
      }
      countOne(kEvTimeRange);
      if (!passItsLayerCut(v == kSystITSLayer ? itsLayerVar : itsLayerNom, col)) {
        continue;
      }
      countOne(kEvITSLayer);
      if (!(v == kSystPileupNSigma ? passAddPileupTight : passAddPileup)) {
        continue;
      }
      countOne(kEvAddPileup);
      // the cut of the event-selection variations
      if (v == kSystPileup && !col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        continue;
      }
      if (v == kSystVz && !(std::abs(col.posZ()) < KSystVtxZCut)) {
        continue;
      }
      if (v == kSystOccupancy && (ev.occupancy < 0 || ev.occupancy >= occupancyMax)) {
        continue;
      }
      countOne(kEvVariationCut);
      if (!(cent < cfgEvSel.cfgCentMax)) {
        continue;
      }
      countOne(kEvCent);
      SETBIT(ev.mask, v);
      if (col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        countOne(kEvInfoPrevRof);
      }
      if (ev.occupancy < 0) {
        countOne(kEvInfoOccUnknown);
      }
    }
    return ev;
  }

  // Track selection of all variations the event passed (evMask): returns the variations the
  // track passes. The cuts common to all variations come first. The cut flow is collected in
  // trkCounts (unless count is false: cut QA) and added to the track counters by flushTrackCounts().
  template <typename T>
  uint32_t isTrackSelected(const T& trk, uint32_t evMask, bool count = true)
  {
    const double inc = count ? 1. : 0.;
    auto countAll = [&](int s) {
      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(evMask, v)) {
          trkCounts[v][s] += inc;
        }
      }
    };
    countAll(kTrkAll);
    if (trk.sign() == 0) {
      return 0;
    }
    countAll(kTrkSign);
    if (!trk.has_collision()) {
      return 0;
    }
    countAll(kTrkHasColl);
    if (!trk.isPVContributor()) {
      return 0;
    }
    countAll(kTrkPVContrib);
    if (trk.pt() < cfgCutPtLower || trk.pt() > cfgCutPtUpper || std::abs(trk.eta()) > cfgCutEta) {
      return 0;
    }
    countAll(kTrkPtEta);
    if (!trk.isGlobalTrack() || trk.tpcNClsCrossedRows() < cfgTPCnCrossedRows) {
      return 0;
    }
    countAll(kTrkGlobal);

    const float maxDcaXY = KSystDcaXyScale * (cfgSyst.cfgDcaXyP0 + cfgSyst.cfgDcaXyP1 / std::pow(trk.pt(), cfgSyst.cfgDcaXyP2.value));
    uint32_t mask = 0;
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(evMask, v)) {
        continue;
      }
      if (!(trk.itsNCls() >= (v == kSystNITS ? KSystITScluster : cfgITScluster.value))) {
        continue;
      }
      trkCounts[v][kTrkITS] += inc;
      if (!(trk.tpcNClsFound() >= (v == kSystNTPC ? KSystTPCcluster : cfgTPCcluster.value))) {
        continue;
      }
      trkCounts[v][kTrkTPC] += inc;
      if (!(trk.tpcChi2NCl() < (v == kSystChi2TPC ? KSystChi2TPC : cfgCutTpcChi2NCl.value))) {
        continue;
      }
      trkCounts[v][kTrkTPCChi2] += inc;
      if (!(std::abs(trk.dcaZ()) < (v == kSystDCAz ? KSystDcaMaxZ : cfgCutTracKDcaMaxZ.value))) {
        continue;
      }
      trkCounts[v][kTrkDCAz] += inc;
      // the cut of the DCAxy variation
      if (v == kSystDCAxy && !(std::abs(trk.dcaXY()) < maxDcaXY)) {
        continue;
      }
      trkCounts[v][kTrkVariationCut] += inc;
      SETBIT(mask, v);
    }
    return mask;
  }

  // Adds the track cut flow of the event to hTrkCountSyst and hTrkCount.
  void flushTrackCounts()
  {
    for (int v = 0; v < kNSyst; ++v) {
      for (int s = 0; s < kNTrkCut; ++s) {
        if (trkCounts[v][s] > 0.) {
          hists.trkCount->Fill(s, v, trkCounts[v][s]);
          if (hists.trkCount1D[v]) {
            hists.trkCount1D[v]->Fill(s, trkCounts[v][s]);
          }
          trkCounts[v][s] = 0.;
        }
      }
    }
  }

  // Charged primary in the observable acceptance (event vertex cuts are applied to the reco collision).
  template <typename T>
  bool isParticleSelected(const T& particle)
  {
    auto* pd = pdg->GetParticle(particle.pdgCode());
    if (!pd) {
      return false;
    }
    if (std::abs(pd->Charge()) == 0) {
      return false;
    }
    if (particle.pt() < cfgCutPtLower || particle.pt() > cfgCutPtUpper || std::abs(particle.eta()) > cfgCutEta) {
      return false;
    }
    return true;
  }

  // Origin of a reco MC track (ETrackOrigin) with respect to the MC collision of coll.
  template <typename T, typename C>
  int trackOrigin(const T& track, const C& coll) const
  {
    if (!track.has_mcParticle()) {
      return kOriginNoLabel;
    }
    const auto& particle = track.mcParticle();
    if (particle.mcCollisionId() != coll.mcCollisionId()) {
      return kOriginOtherColl;
    }
    return particle.isPhysicalPrimary() ? kOriginPrimary : kOriginSecondary;
  }

  float getCentrality(const auto& col, int choice) const
  {
    if (choice == kCentFT0C) {
      return col.centFT0C();
    }
    if (choice == kCentFT0M) {
      return col.centFT0M();
    }
    if (choice == kCentFDDM) {
      return col.centFDDM();
    }
    if (choice == kCentFV0A) {
      return col.centFV0A();
    }
    return KinvalidCentrality;
  }

  // Observable eta bins containing eta: the full-range bin 0 and one narrow bin.
  int etaBinsOf(float eta, std::array<int, 2>& bins) const
  {
    int n = 0;
    for (int ieta = 0; ieta < nEta && n < 2; ++ieta) {
      if (eta > etaLw[ieta] && eta <= etaUp[ieta]) {
        bins[n++] = ieta;
      }
    }
    return n;
  }

  // ===========================================================================
  // Corrections
  // ===========================================================================
  // Inclusive efficiency/fake lookup (no species dependence).
  static float getEfficiency(const TH3* h, float cent, float pt, float eta, bool isEff)
  {
    const float invalid = std::numeric_limits<float>::quiet_NaN();
    if (!h || !std::isfinite(cent) || !std::isfinite(pt) || !std::isfinite(eta)) {
      return invalid;
    }
    int ibx = h->GetXaxis()->FindFixBin(cent);
    // Reject centrality underflow, clamp overflow.
    if (ibx < 1) {
      return invalid;
    }
    ibx = std::min(ibx, h->GetNbinsX());
    int iby = h->GetYaxis()->FindFixBin(pt);
    int ibz = h->GetZaxis()->FindFixBin(eta);
    // Do not extrapolate in pT or eta.
    if (iby < 1 || iby > h->GetNbinsY() || ibz < 1 || ibz > h->GetNbinsZ()) {
      return invalid;
    }
    const double val = h->GetBinContent(ibx, iby, ibz);
    const double err = h->GetBinError(ibx, iby, ibz);
    // A bin without a positive error is not a valid measurement, except a fake rate of 0
    // (no fake track in the bin), which has no error.
    if (!std::isfinite(val) || !std::isfinite(err) || (err <= 0.0 && (isEff || val != 0.0))) {
      return invalid;
    }
    // Physical-domain checks.
    if (val < 0.0 || val > 1.0 || (isEff && val == 0.0)) {
      return invalid;
    }
    return static_cast<float>(val);
  }

  // Nominal efficiency and fake rate of the map of variation v (no variation factor).
  bool getValidEffFake(int v, float cent, float pt, float eta, float& eff, float& fake)
  {
    if (!cfgEff) {
      eff = 1.f;
      fake = 0.f;
      return true;
    }
    const auto& in = inputs[v];
    auto valid = [&]() { return std::isfinite(eff) && std::isfinite(fake) && eff > KFloatEpsilon && eff <= 1.f && fake >= 0.f && fake < 1.f; };
    eff = getEfficiency(in.hEff, cent, pt, eta, true);
    fake = getEfficiency(in.hFake, cent, pt, eta, false);
    if (valid()) {
      return true;
    }
    // empty or statistically invalid cell, or no valid centrality: the map integrated over centrality
    const float centAll = in.hEffPtEta ? in.hEffPtEta->GetXaxis()->GetBinCenter(1) : 0.f;
    eff = getEfficiency(in.hEffPtEta, centAll, pt, eta, true);
    fake = getEfficiency(in.hFakePtEta, centAll, pt, eta, false);
    return valid();
  }

  float efficiencyVariationFactor(int v, float pt)
  {
    if (v == kSystEff) {
      return 1.f + 0.01f * KEffVarPercent;
    }
    if (v == kSystEffPtDep) {
      // Linear in pT: full variation at cfgPtMin, none at KEffVarPtZero.
      const float shape = std::clamp((KEffVarPtZero - pt) / (KEffVarPtZero - cfgPtMin), 0.f, 1.f);
      return 1.f + 0.01f * KEffVarPercent * shape;
    }
    return 1.f;
  }

  static float getFlatteningWeight(const THnSparseF* h, float vz, float chg, float pt, float eta, float phi)
  {
    if (!h) {
      return 0.0;
    }
    std::array<int, 5> bins{};
    bins[0] = h->GetAxis(0)->FindBin(vz);
    bins[1] = h->GetAxis(1)->FindBin(chg);
    bins[2] = h->GetAxis(2)->FindBin(pt);
    bins[3] = h->GetAxis(3)->FindBin(eta);
    bins[4] = h->GetAxis(4)->FindBin(phi);
    return h->GetBinContent(bins.data());
  }

  // Weight w = w_flat (1 - fake) / eff of the track for every variation in mask; a map shared
  // by several variations is looked up once. Returns the variations with a usable weight and
  // the Base efficiency, fake rate and flattening weight for the QA profiles.
  uint32_t getTrackWeights(uint32_t mask, float cent, float vz, float sign, float pt, float eta, float phi, bool useFlat,
                           std::array<float, kNSyst>& weights, float& effBase, float& fakeBase, float& flatBase)
  {
    std::array<int8_t, kNSyst> effDone{};
    std::array<float, kNSyst> effS{}, fakeS{}, flatS{};
    std::array<bool, kNSyst> flatDone{};
    uint32_t weightMask = 0;
    effBase = std::numeric_limits<float>::quiet_NaN();
    fakeBase = std::numeric_limits<float>::quiet_NaN();
    flatBase = 1.f;
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(mask, v)) {
        continue;
      }
      const int es = effSource(v);
      if (effDone[es] == 0) {
        effDone[es] = getValidEffFake(es, cent, pt, eta, effS[es], fakeS[es]) ? 1 : -1;
      }
      if (effDone[es] < 0) {
        continue;
      }
      const float eff = effS[es] * efficiencyVariationFactor(v, pt);
      if (!(eff > KFloatEpsilon)) {
        continue;
      }
      float flat = 1.f;
      const int fs = flatSource(v);
      if (useFlat && cfgFlat && fs >= 0) {
        if (!flatDone[fs]) {
          flatS[fs] = getFlatteningWeight(inputs[fs].hFlat, vz, sign, pt, eta, phi);
          flatDone[fs] = true;
        }
        flat = flatS[fs];
      }
      const float weight = flat * (1.0f - fakeS[es]) / eff;
      if (!std::isfinite(weight) || weight <= 0.f) {
        continue;
      }
      weights[v] = weight;
      SETBIT(weightMask, v);
      trkCounts[v][kTrkWeight] += 1.;
      if (v == kSystBase) {
        effBase = eff;
        fakeBase = fakeS[es];
        flatBase = flat;
      }
    }
    return weightMask;
  }

  // Per-event two- and three-particle pT correlators (standard method) from the
  // power sums, following arXiv:2112.03397, Eqs. 2-3.
  struct C2C3Result {
    double c2 = std::numeric_limits<double>::quiet_NaN();
    double c3 = std::numeric_limits<double>::quiet_NaN();
  };

  template <int M, int K>
  C2C3Result calculateC2C3FromSums(const std::array<std::array<double, K>, M>& sumpmwk, const std::array<double, K>& sumwk, float referenceMeanPt) const
  {
    C2C3Result r;
    const double s1 = sumwk[1];
    if (s1 <= 0.) {
      return r;
    }
    const double refpt = referenceMeanPt;
    const double p11 = sumpmwk[1][1] / s1;
    const double p1Bar1 = p11 - refpt;

    // --- c2: needs >=2 particles (D2 = 1 - tau1 > eps) ---
    const double tau1 = sumwk[2] / (s1 * s1);
    const double denom2 = 1. - tau1;
    if (denom2 > KFloatEpsilon) {
      const double p12 = sumpmwk[1][2] / sumwk[2];
      const double p22 = sumpmwk[2][2] / sumwk[2];
      const double p2Bar2 = p22 - 2.0 * p12 * refpt + refpt * refpt;
      r.c2 = (p1Bar1 * p1Bar1 - tau1 * p2Bar2) / denom2;
    }

    // --- c3: needs >=3 particles (D3 = 1 - 3 tau1 + 2 tau2 > eps) ---
    const double tau2 = sumwk[3] / (s1 * s1 * s1);
    const double denom3 = 1. - 3.0 * tau1 + 2.0 * tau2;
    if (denom3 > KFloatEpsilon) {
      const double p12 = sumpmwk[1][2] / sumwk[2];
      const double p22 = sumpmwk[2][2] / sumwk[2];
      const double p13 = sumpmwk[1][3] / sumwk[3];
      const double p23 = sumpmwk[2][3] / sumwk[3];
      const double p33 = sumpmwk[3][3] / sumwk[3];
      const double p2Bar2 = p22 - 2.0 * p12 * refpt + refpt * refpt;
      const double p3Bar3 = p33 - 3.0 * p23 * refpt + 3.0 * p13 * refpt * refpt - refpt * refpt * refpt;
      r.c3 = (p1Bar1 * p1Bar1 * p1Bar1 - 3.0 * tau1 * p2Bar2 * p1Bar1 + 2.0 * tau2 * p3Bar3) / denom3;
    }
    return r;
  }

  // ===========================================================================
  // Forward detectors
  // ===========================================================================
  // Eta sub-bin of detector det containing eta, -1 outside its edges.
  int fwdSubBin(int det, float eta) const
  {
    const auto& e = fwdEdges[det];
    if (e.size() < 2 || eta < e.front() || eta >= e.back()) {
      return -1;
    }
    const int i = static_cast<int>(std::upper_bound(e.begin(), e.end(), eta) - e.begin()) - 1;
    return fwdOffset[det] + i;
  }

  void addFwd(int det, float eta, float ampl, FwdSignal& fwd)
  {
    fwd.ampl[det] += ampl;
    const int sub = fwdSubBin(det, eta);
    if (sub >= 0) {
      fwd.ampl[sub] += ampl;
    }
    if (hists.fwdEta) {
      hists.fwdEta->Fill(det, eta, ampl);
      if (sub >= 0) {
        hists.fwdEta->Fill(sub, eta, ampl);
      }
    }
  }

  void markFwdValid(int det, FwdSignal& fwd) const
  {
    fwd.valid[det] = true;
    for (int i = 0; i + 1 < static_cast<int>(fwdEdges[det].size()); ++i) {
      fwd.valid[fwdOffset[det] + i] = true;
    }
  }

  float channelEta(const std::array<float, 3>& p, float vz)
  {
    const float r = std::hypot(p[0], p[1]);
    const float dz = p[2] - (cfgFwdOpt.cfgFwdVzCorr ? vz : 0.f);
    return r > 0.f ? std::asinh(dz / r) : std::numeric_limits<float>::quiet_NaN();
  }

  void fillChannelQA(int ch, float raw, float corr, float eta)
  {
    if (hists.fwdChannel) {
      hists.fwdChannel->Fill(ch, kChAmpRaw, raw);
      hists.fwdChannel->Fill(ch, kChAmpCorr, corr);
      hists.fwdChannel->Fill(ch, kChEta, eta);
    }
  }

  // Channel positions with the FT0/FV0 alignment (as in PWGCF/flowDecorrelation) and the
  // gain equalisation, loaded at each new run.
  template <typename B>
  void loadFwdCalib(const B& bc, bool applyGain)
  {
    if (bc.runNumber() == fwdCalib.run) {
      return;
    }
    fwdCalib.run = bc.runNumber();
    const uint64_t ts = bc.timestamp();

    auto* alFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", ts);
    auto* alFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", ts);
    if (!alFT0 || alFT0->size() < 2 || !alFV0 || alFV0->size() < 2) {
      LOGF(fatal, "Could not load FT0/FV0 alignment (FT0/Calib/Align, FV0/Calib/Align) for timestamp %llu", ts);
      return;
    }

    // offsetFT0[0]: FT0A, offsetFT0[1]: FT0C
    fwdCalib.ft0Det.calculateChannelCenter();
    for (int ch = 0; ch < KNChFT0; ++ch) {
      const int side = ch < KNChFT0A ? 0 : 1;
      const auto c = fwdCalib.ft0Det.getChannelCenter(ch);
      float z = c.Z() + (*alFT0)[side].getZ();
      if (ch >= KNChFT0A) {
        z = -z;
      }
      fwdCalib.posFT0[ch] = {static_cast<float>(c.X() + (*alFT0)[side].getX()), static_cast<float>(c.Y() + (*alFT0)[side].getY()), z};
    }
    // offsetFV0[0]: left half, offsetFV0[1]: right half
    static constexpr std::array<int, 24> KCellsInLeft = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    for (int ch = 0; ch < KNChFV0; ++ch) {
      const bool left = std::find(KCellsInLeft.begin(), KCellsInLeft.end(), ch) != KCellsInLeft.end();
      const auto& off = (*alFV0)[left ? 0 : 1];
      const auto c = fwdCalib.fv0Det->getReadoutCenter(ch);
      fwdCalib.posFV0[ch] = {static_cast<float>(c.x + off.getX()), static_cast<float>(c.y + off.getY()), static_cast<float>(c.z + off.getZ())};
    }

    fwdCalib.gainFT0.fill(1.f);
    fwdCalib.gainFV0.fill(1.f);
    if (applyGain && cfgFwdOpt.cfgFwdGainEq) {
      auto* gFT0 = ccdb->getForTimeStamp<std::vector<float>>(cfgFwdOpt.cfgGainEqPath.value + "/FT0", ts);
      auto* gFV0 = ccdb->getForTimeStamp<std::vector<float>>(cfgFwdOpt.cfgGainEqPath.value + "/FV0", ts);
      if (gFT0 && gFT0->size() >= static_cast<size_t>(KNChFT0)) {
        std::copy_n(gFT0->begin(), KNChFT0, fwdCalib.gainFT0.begin());
      } else {
        LOGF(warning, "No FT0 gain equalisation at %s/FT0 for run %d: unit gains", cfgFwdOpt.cfgGainEqPath.value.c_str(), fwdCalib.run);
      }
      if (gFV0 && gFV0->size() >= static_cast<size_t>(KNChFV0)) {
        std::copy_n(gFV0->begin(), KNChFV0, fwdCalib.gainFV0.begin());
      } else {
        LOGF(warning, "No FV0 gain equalisation at %s/FV0 for run %d: unit gains", cfgFwdOpt.cfgGainEqPath.value.c_str(), fwdCalib.run);
      }
    }
    LOGF(info, "Forward-detector calibration loaded for run %d", fwdCalib.run);
  }

  // Gain-equalised amplitudes of FT0A, FT0C and FV0A, totals and eta sub-bins.
  template <typename C>
  void computeFwdReco(const C& col, float vz, FwdSignal& fwd, bool qa)
  {
    fwd = {};
    auto corr = [](float a, float g) { return g > 0.f ? a / g : a; };
    if (col.has_foundFT0()) {
      const auto& ft0 = col.foundFT0();
      markFwdValid(kFT0A, fwd);
      markFwdValid(kFT0C, fwd);
      for (std::size_t i = 0; i < ft0.channelA().size(); ++i) {
        const int id = ft0.channelA()[i];
        if (id < 0 || id >= KNChFT0A) {
          continue;
        }
        const float a = corr(ft0.amplitudeA()[i], fwdCalib.gainFT0[id]);
        const float eta = channelEta(fwdCalib.posFT0[id], vz);
        addFwd(kFT0A, eta, a, fwd);
        if (qa) {
          fillChannelQA(id, ft0.amplitudeA()[i], a, eta);
        }
      }
      for (std::size_t i = 0; i < ft0.channelC().size(); ++i) {
        const int id = ft0.channelC()[i] + KNChFT0A;
        if (id < KNChFT0A || id >= KNChFT0) {
          continue;
        }
        const float a = corr(ft0.amplitudeC()[i], fwdCalib.gainFT0[id]);
        const float eta = channelEta(fwdCalib.posFT0[id], vz);
        addFwd(kFT0C, eta, a, fwd);
        if (qa) {
          fillChannelQA(id, ft0.amplitudeC()[i], a, eta);
        }
      }
    }
    if (col.has_foundFV0()) {
      const auto& fv0 = col.foundFV0();
      markFwdValid(kFV0A, fwd);
      for (std::size_t i = 0; i < fv0.channel().size(); ++i) {
        const int id = fv0.channel()[i];
        if (id < 0 || id >= KNChFV0) {
          continue;
        }
        const float a = corr(fv0.amplitude()[i], fwdCalib.gainFV0[id]);
        const float eta = channelEta(fwdCalib.posFV0[id], vz);
        addFwd(kFV0A, eta, a, fwd);
        if (qa) {
          fillChannelQA(KNChFT0 + id, fv0.amplitude()[i], a, eta);
        }
      }
    }
  }

  // Number of charged primaries in the acceptance of each forward detector (MC truth).
  template <typename P>
  void computeFwdTruth(const P& particles, FwdSignal& fwd)
  {
    fwd = {};
    for (int det = 0; det < kNFwdDet; ++det) {
      markFwdValid(det, fwd);
    }
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      const float eta = particle.eta();
      if (std::abs(eta) < 1.5f) {
        continue;
      }
      auto* pd = pdg->GetParticle(particle.pdgCode());
      if (!pd || std::abs(pd->Charge()) == 0) {
        continue;
      }
      for (int det = 0; det < kNFwdDet; ++det) {
        const int sub = fwdSubBin(det, eta);
        if (sub >= 0) {
          fwd.ampl[det] += 1.;
          fwd.ampl[sub] += 1.;
        }
      }
    }
  }

  // ===========================================================================
  // Observables
  // ===========================================================================
  // Bin content of p at (x, ybin[, zbin]), NaN for an empty bin.
  static double profileValue(const TProfile2D* p, double x, int ybin)
  {
    if (!p) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const int gb = p->GetBin(p->GetXaxis()->FindBin(x), ybin);
    if (p->GetBinEntries(gb) <= 0.) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return p->GetBinContent(gb);
  }
  static double profileValue(const TProfile3D* p, double x, int ybin, int zbin)
  {
    if (!p) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const int gb = p->GetBin(p->GetXaxis()->FindBin(x), ybin, zbin);
    if (p->GetBinEntries(gb) <= 0.) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return p->GetBinContent(gb);
  }

  void computeEventObs(const PowerSums& s, int v, int level, const FwdSignal* fwd, float nCh, EventObs& evObs)
  {
    const auto& input = inputs[v];
    const double nan = std::numeric_limits<double>::quiet_NaN();
    evObs.nCh = nCh;
    for (int ieta = 0; ieta < nEta; ++ieta) {
      evObs.nW[ieta] = s.wk[ieta][1];
      evObs.mean[ieta] = (s.wk[ieta][1] > KFloatEpsilon) ? (s.pmwk[ieta][1][1] / s.wk[ieta][1]) : nan;
      evObs.c2[ieta] = nan;
      evObs.c3[ieta] = nan;
      evObs.p1kBar[ieta] = nan;
      // c2 (>=2 tracks), c3 (>=3 tracks) and the pT deviation need a usable reference
      const double mmpt = profileValue(input.pmean[level], nCh, ieta + 1, kMeanPt + 1);
      if (std::isfinite(mmpt) && mmpt != 0) {
        const auto corr = calculateC2C3FromSums<KIntM, KIntK>(s.pmwk[ieta], s.wk[ieta], mmpt);
        evObs.c2[ieta] = corr.c2;
        evObs.c3[ieta] = corr.c3;
        if (std::isfinite(evObs.mean[ieta])) {
          evObs.p1kBar[ieta] = evObs.mean[ieta] - mmpt;
        }
      }
      evObs.dN[ieta] = evObs.nW[ieta] - profileValue(input.pmean[level], nCh, ieta + 1, kMeanN + 1);
    }
    evObs.fwdValid.fill(false);
    if (!fwd || !input.pmeanFwd[level]) {
      return;
    }
    for (int k = 0; k < nFwd; ++k) {
      if (!fwd->valid[k]) {
        continue;
      }
      const double m = profileValue(input.pmeanFwd[level], nCh, k + 1);
      if (!std::isfinite(m)) {
        continue;
      }
      evObs.fwd[k] = fwd->ampl[k];
      evObs.dFwd[k] = fwd->ampl[k] - m;
      evObs.fwdValid[k] = true;
    }
  }

  // Fill the fluctuation profiles of one (variation, level) at x = N_ch, or of one bootstrap
  // replica with weight w at x = centrality.
  void fillFluc(const FlucHists& h, const EventObs& evObs, double w, float x)
  {
    auto fill = [w](const std::shared_ptr<TProfile3D>& p, double a, double b, double c, double val) {
      if (p) {
        p->Fill(a, b, c, val, w);
      }
    };

    // multiplicity, meanpT, C2, C3 vs eta bin
    for (int ieta = 0; ieta < nEta; ++ieta) {
      fill(h.obs, x, ieta, kObsN, evObs.nW[ieta]);
      if (std::isfinite(evObs.mean[ieta])) {
        fill(h.obs, x, ieta, kObsMeanPt, evObs.mean[ieta]);
      }
      if (std::isfinite(evObs.c2[ieta])) {
        fill(h.obs, x, ieta, kObsC2, evObs.c2[ieta]);
      }
      if (std::isfinite(evObs.c3[ieta])) {
        fill(h.obs, x, ieta, kObsC3, evObs.c3[ieta]);
      }
    }

    // mirror-pair subevent (C2Sub, C3Sub) vs eta bin
    for (int ietaA = 1; ietaA <= (nEta - 1) / 2; ++ietaA) {
      int ietaC = nEta - ietaA;
      float c2Sub = evObs.p1kBar[ietaA] * evObs.p1kBar[ietaC];
      // Mirror pairs (A, C = -A): both orderings are filled, each at the bin of its pair,
      // so bins A < 0 hold AAB and bins C > 0 hold ABB of the same mirror pair.  For a
      // symmetric collision system the two agree in expectation: eta -> -eta maps
      // AAB(a, -a) onto AAB(-a, a) = ABB(a, -a).  Comparing bin A with bin nEta - A is
      // therefore a forward/backward check, not two independent observables.
      float c3SubA = evObs.c2[ietaA] * evObs.p1kBar[ietaC]; // AAB: 2 particles from A, 1 from C
      float c3SubC = evObs.c2[ietaC] * evObs.p1kBar[ietaA]; // ABB: 2 particles from C, 1 from A

      if (std::isfinite(c2Sub)) {
        fill(h.obs, x, ietaA, kObsC2Sub, c2Sub);
      }
      if (std::isfinite(c3SubA)) {
        fill(h.obs, x, ietaA, kObsC3Sub, c3SubA);
      }
      if (std::isfinite(c3SubC)) {
        fill(h.obs, x, ietaC, kObsC3Sub, c3SubC);
      }
    }

    // full 2D subevent maps, in (eta_A, eta_C) and in (Delta eta, Sigma eta)
    for (int ietaA = 1; ietaA < nEta; ++ietaA) {
      for (int ietaC = 1; ietaC < nEta; ++ietaC) {
        float etaValA = (etaLw[ietaA] + etaUp[ietaA]) / 2.0f;
        float etaValB = (etaLw[ietaC] + etaUp[ietaC]) / 2.0f;
        const float bw = etaUp[ietaA] - etaLw[ietaA];
        const float gap = (ietaA - ietaC) * bw;
        const float sum = 2.f * etaLw[1] + (ietaA + ietaC - 1) * bw; // = cA + cC on the lattice

        float c2Sub = (ietaA == ietaC) ? static_cast<float>(evObs.c2[ietaA]) : evObs.p1kBar[ietaA] * evObs.p1kBar[ietaC];
        // C3 sub-event map, ORDERED: cell (A,C) = c2[A] p1kBar[C] = <dpT dpT>_A <dpT>_C,
        // two particles from A and one from C ("AAB").  The same pair of bins read the
        // other way round, cell (C,A) = c2[C] p1kBar[A], is "ABB" (one from A, two from C).
        // Both are filled, so the map carries AAB and ABB of every bin pair: the row at
        // eta_A = eta_ref is AAB with the pair at eta_ref, the column is ABB.  Unlike C2,
        // this map is NOT symmetric under A <-> C, and the two triangles must not be
        // averaged or treated as duplicates.  Diagonal: the within-bin c3 ("AAA").
        float c3Sub2D = (ietaA == ietaC) ? static_cast<float>(evObs.c3[ietaA]) : evObs.c2[ietaA] * evObs.p1kBar[ietaC];

        if (std::isfinite(c2Sub)) {
          fill(h.c2Sub2D, x, etaValA, etaValB, c2Sub);
          fill(h.gapSum2D, x, gap, sum, c2Sub);
        }
        if (std::isfinite(c3Sub2D)) {
          fill(h.c3Sub2D, x, etaValA, etaValB, c3Sub2D);
          fill(h.c3GapSum2D, x, gap, sum, c3Sub2D);
        }
      }
    }

    // forward amplitude a_k: <da^2>, <d[pT](eta) da> and <dN(eta) da>, with
    // da = a - <a>, d[pT] = [pT] - <[pT]>, dN = N - <N>, the means at N_ch from the mean pass
    for (int k = 0; k < nFwd; ++k) {
      if (!evObs.fwdValid[k]) {
        continue;
      }
      const double da = evObs.dFwd[k];
      if (h.fwdVar) {
        h.fwdVar->Fill(x, k, da * da, w);
      }
      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (std::isfinite(evObs.p1kBar[ieta])) {
          fill(h.covFwdPt, x, ieta, k, evObs.p1kBar[ieta] * da);
        }
        if (std::isfinite(evObs.dN[ieta])) {
          fill(h.covFwdN, x, ieta, k, evObs.dN[ieta] * da);
        }
      }
    }
  }

  // ===========================================================================
  // Table joins
  // ===========================================================================
  using GeneralCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                      aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFDDMs, aod::CentFV0As,
                                      aod::CentNTPVs>;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgEvSel.cfgVtxZCut;
  using AodCollisionsSel = soa::Filtered<GeneralCollisions>;

  using UnfilteredTracks = soa::Join<
    aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  Filter trackFilter = aod::track::pt > KPtMin&& aod::track::pt < KPtMax&& requireGlobalTrackInFilter();
  using AodTracksSel = soa::Filtered<UnfilteredTracks>;
  using TCs = soa::Join<UnfilteredTracks, aod::McTrackLabels>;
  using FilteredTCs = soa::Filtered<TCs>;

  using MyRun3MCCollisions = soa::Join<
    aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra,
    aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFDDMs, aod::CentFV0As,
    aod::CentNGlobals, aod::McCollisionLabels>;
  // reconstructed collisions of one MC collision
  using McRecoCollisions = soa::SmallGroups<MyRun3MCCollisions>;

  Preslice<FilteredTCs> tracksPerCollision = aod::track::collisionId;

  // ===========================================================================
  // Histogram declarations
  // ===========================================================================
  // Bin labels only on TH1/TH2 and THnSparse axes: ROOT merges labelled axes by label, which for
  // profiles and TH3 extends the axis and sums all bins into a new one, or drops the entries of
  // bins with zero content. Their bins are named in the axis title instead.
  static void labelAxis(TAxis* ax, const std::vector<const char*>& labels)
  {
    for (std::size_t i = 0; i < labels.size() && static_cast<int>(i) < ax->GetNbins(); ++i) {
      ax->SetBinLabel(i + 1, labels[i]);
    }
  }

  void labelSystAxis(TAxis* ax) const
  {
    for (int v = 0; v < kNSyst && v < ax->GetNbins(); ++v) {
      ax->SetBinLabel(v + 1, systName[v].c_str());
    }
  }

  void declareCommonQA()
  {
    histos.add("hVtxZ_after_sel", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
    histos.add("hVtxZ", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{centAxis1Per}});
    histos.add("Hist2D_globalTracks_PVTracks", ";N_{PV};N_{ch} (global tracks)", kTH2F, {{multQAAxis}, {multQAAxis}});
    histos.add("Hist2D_PVTracks_cent", "cent (%);N_{PV}", kTH2F, {{centAxis1Per}, {multQAAxis}});
    hists.nchOccupancy = histos.add<TH2>("Hist2D_Nch_occupancy", ";N_{ch} (global tracks);trackOccupancyInTimeRange", kTH2F, {{nchAxisObs}, {KNbinsOccQA, 0., KOccQAMax}});

    histos.add("hP", ";p (GeV/c)", kTH1F, {{KNbinsPt, KPMin, KPMax}});
    histos.add("hPt", ";p_{T} (GeV/c)", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
    histos.add("hEta", ";#eta", kTH1F, {{KNbinsEtaFine, KEtaMin, KEtaMax}});
    histos.add("hPhi", ";#phi", kTH1F, {{KNbinsPhi, KPhiMin, TwoPI}});

    // cut flows of every variation (y bin 1 = Base); 1D copies per variation in declareCounts()
    hists.evtCount = histos.add<TH2>("hEvtCountSyst", "Number of events;;variation", kTH2D, {{kNEvCut, -0.5, static_cast<double>(kNEvCut) - 0.5}, {kNSyst, -0.5, static_cast<double>(kNSyst) - 0.5}});
    labelAxis(hists.evtCount->GetXaxis(), evtCutLabels());
    labelSystAxis(hists.evtCount->GetYaxis());

    hists.trkCount = histos.add<TH2>("hTrkCountSyst", "Number of tracks;;variation", kTH2D, {{kNTrkCut, -0.5, static_cast<double>(kNTrkCut) - 0.5}, {kNSyst, -0.5, static_cast<double>(kNSyst) - 0.5}});
    labelAxis(hists.trkCount->GetXaxis(), trkCutLabels());
    labelSystAxis(hists.trkCount->GetYaxis());
  }

  static std::vector<const char*> evtCutLabels()
  {
    return {"all", "sel8", "|z_{vtx}|", "kIsGoodZvtxFT0vsPV", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kIsVertexITSTPC",
            "kNoSameBunchPileup", "kNoCollInTimeRangeStandard",
            "ITS good layers", "n#sigma/linear pile-up", "variation cut", "cent < cfgCentMax (selected)",
            "selected: kNoHighMultCollInPrevRof", "selected: occupancy < 0"};
  }
  static std::vector<const char*> trkCutLabels()
  {
    return {"all", "sign#neq0", "has_collision", "isPVContributor", "p_{T}, #eta", "isGlobalTrack, crossed rows",
            "itsNCls", "tpcNClsFound", "tpcChi2NCl", "DCAz", "variation cut (selected)", "valid weight (used)"};
  }

  // hEvtCount and hTrkCount of each variation of this pass (Base at the top level, the others
  // in their folder), so that every variation's output carries its own cut flow.
  void declareCounts(uint32_t mask)
  {
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(mask, v)) {
        continue;
      }
      auto& r = reg(v);
      hists.evtCount1D[v] = r.add<TH1>((systPrefix(v) + "hEvtCount").c_str(), ("Number of events, " + systName[v]).c_str(), kTH1D, {{kNEvCut, -0.5, static_cast<double>(kNEvCut) - 0.5}});
      labelAxis(hists.evtCount1D[v]->GetXaxis(), evtCutLabels());
      hists.trkCount1D[v] = r.add<TH1>((systPrefix(v) + "hTrkCount").c_str(), ("Number of tracks, " + systName[v]).c_str(), kTH1D, {{kNTrkCut, -0.5, static_cast<double>(kNTrkCut) - 0.5}});
      labelAxis(hists.trkCount1D[v]->GetXaxis(), trkCutLabels());
    }
  }

  void declareCutQA()
  {
    const AxisSpec ptAxis{100, 0., 10., "p_{T} (GeV/c)"};
    for (int st = 0; st < kNQAStage; ++st) {
      const std::string d = st == kQABefore ? "CutQA/before/" : "CutQA/after/";
      auto h1 = [&](const char* name, const char* title, const AxisSpec& ax) { return histos.add<TH1>((d + name).c_str(), title, kTH1F, {ax}); };
      auto h2 = [&](const char* name, const char* title, const AxisSpec& ax, const AxisSpec& ay) { return histos.add<TH2>((d + name).c_str(), title, kTH2F, {ax, ay}); };
      hists.qaEvt1D[kQAVtxZ][st] = h1("hVtxZ", ";z_{vtx} (cm)", {200, -20., 20.});
      hists.qaEvt1D[kQACent][st] = h1("hCent", ";centrality (%)", {centAxis1Per});
      hists.qaEvt1D[kQAOccupancy][st] = h1("hOccupancy", ";trackOccupancyInTimeRange", {KNbinsOccQA, 0., KOccQAMax});
      hists.qaEvt2D[kQANchNPV][st] = h2("hNch_NPV", ";N_{PV};N_{ch} (global tracks)", multQAAxis, multQAAxis);
      hists.qaEvt2D[kQANchCent][st] = h2("hNch_cent", ";centrality (%);N_{ch} (global tracks)", centAxis1Per, multQAAxis);
      hists.qaTrk1D[kQAPt][st] = h1("hTrkPt", "", ptAxis);
      hists.qaTrk1D[kQAEta][st] = h1("hTrkEta", ";#eta", {60, -1.5, 1.5});
      hists.qaTrk1D[kQAItsNCls][st] = h1("hTrkItsNCls", ";ITS clusters", {8, -0.5, 7.5});
      hists.qaTrk1D[kQAItsChi2][st] = h1("hTrkItsChi2NCl", ";ITS #chi^{2}/cluster", {80, 0., 40.});
      hists.qaTrk1D[kQATpcNClsFound][st] = h1("hTrkTpcNClsFound", ";TPC found clusters", {161, -0.5, 160.5});
      hists.qaTrk1D[kQATpcCrossedRows][st] = h1("hTrkTpcCrossedRows", ";TPC crossed rows", {161, -0.5, 160.5});
      hists.qaTrk1D[kQATpcCrossedRowsOverFindable][st] = h1("hTrkTpcCrossedRowsOverFindable", ";TPC crossed rows / findable clusters", {75, 0., 1.5});
      hists.qaTrk1D[kQATpcChi2][st] = h1("hTrkTpcChi2NCl", ";TPC #chi^{2}/cluster", {100, 0., 10.});
      hists.qaTrk1D[kQAPVContributor][st] = h1("hTrkPVContributor", ";isPVContributor", {2, -0.5, 1.5});
      hists.qaTrk1D[kQAGlobalTrack][st] = h1("hTrkGlobalTrack", ";isGlobalTrack", {2, -0.5, 1.5});
      hists.qaTrk2D[kQADcaXYPt][st] = h2("hTrkDcaXY_pt", ";p_{T} (GeV/c);DCA_{xy} (cm)", ptAxis, {400, -1., 1.});
      hists.qaTrk2D[kQADcaZPt][st] = h2("hTrkDcaZ_pt", ";p_{T} (GeV/c);DCA_{z} (cm)", ptAxis, {300, -3., 3.});
    }
  }

  template <typename C>
  void fillCutQAEvent(int st, const C& col, int nCh, float cent)
  {
    hists.qaEvt1D[kQAVtxZ][st]->Fill(col.posZ());
    hists.qaEvt1D[kQACent][st]->Fill(cent);
    hists.qaEvt1D[kQAOccupancy][st]->Fill(col.trackOccupancyInTimeRange());
    hists.qaEvt2D[kQANchNPV][st]->Fill(col.multNTracksPV(), nCh);
    hists.qaEvt2D[kQANchCent][st]->Fill(cent, nCh);
  }

  template <typename T>
  void fillCutQATrack(int st, const T& trk)
  {
    hists.qaTrk1D[kQAPt][st]->Fill(trk.pt());
    hists.qaTrk1D[kQAEta][st]->Fill(trk.eta());
    hists.qaTrk1D[kQAItsNCls][st]->Fill(trk.itsNCls());
    hists.qaTrk1D[kQAItsChi2][st]->Fill(trk.itsChi2NCl());
    hists.qaTrk1D[kQATpcNClsFound][st]->Fill(trk.tpcNClsFound());
    hists.qaTrk1D[kQATpcCrossedRows][st]->Fill(trk.tpcNClsCrossedRows());
    hists.qaTrk1D[kQATpcCrossedRowsOverFindable][st]->Fill(trk.tpcCrossedRowsOverFindableCls());
    hists.qaTrk1D[kQATpcChi2][st]->Fill(trk.tpcChi2NCl());
    hists.qaTrk1D[kQAPVContributor][st]->Fill(trk.isPVContributor());
    hists.qaTrk1D[kQAGlobalTrack][st]->Fill(trk.isGlobalTrack());
    hists.qaTrk2D[kQADcaXYPt][st]->Fill(trk.pt(), trk.dcaXY());
    hists.qaTrk2D[kQADcaZPt][st]->Fill(trk.pt(), trk.dcaZ());
  }

  // (centrality, N_ch) of the selected events of every variation: its event count and the
  // centrality of each N_ch bin with the estimator of the variation.
  void declareNchCent(uint32_t mask)
  {
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(mask, v) || hists.nchCent[v]) {
        continue;
      }
      hists.nchCent[v] = reg(v).add<TH2>((systPrefix(v) + "Hist2D_globalTracks_cent").c_str(), "cent (%);N_{ch} (global tracks)", kTH2F, {{centAxis1Per}, {nchAxisOf(v)}});
    }
  }

  void declareFwdQA()
  {
    if (hists.fwdChannel) {
      return;
    }
    hists.fwdChannel = histos.add<TProfile2D>("FwdQA/Prof_Channel_qa", ";channel (FT0A 0-95, FT0C 96-207, FV0A 208-255);quantity (0 amplitude, 1 gain-equalised amplitude, 2 #eta aligned vs z_{vtx})", kTProfile2D, {{KNChFT0 + KNChFV0, -0.5, KNChFT0 + KNChFV0 - 0.5}, {kNChQA, -0.5, static_cast<double>(kNChQA) - 0.5}});
    hists.fwdEta = histos.add<TProfile>("FwdQA/Prof_Eta_fwdbin", ";forward bin;amplitude-weighted #LT#eta#GT", kTProfile, {fwdAxis});
  }

  void declareEff(uint32_t mask)
  {
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(mask, v)) {
        continue;
      }
      for (int c = 0; c < kNEffCat; ++c) {
        hists.eff[v].maps[c] = reg(v).add<TH3>((systPrefix(v) + effMapName[c]).c_str(), ";centrality (%);p_{T} (GeV/c);#eta", kTH3F,
                                               {{centAxisEff}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      }
    }
  }

  void declareMCQA()
  {
    histos.add("ptResolution", ";p_{T}^{MC};(p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}", kTH2F, {{KNbinsPtRes, KPtMin, KPtMax}, {100, -0.2, 0.2}});
    histos.add("etaResolution", ";#eta^{MC};#eta^{reco}-#eta^{MC}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {100, -0.02, 0.02}});
    histos.add("etaTruthReco", ";#eta^{MC};#eta^{reco}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}});
    histos.add("TruthTracKVz", ";Vz^{MC};Vz^{Reco}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {KNbinsVz, KVzMin, KVzMax}});
    histos.add("vzResolution", ";Vz^{MC};(Vz^{reco}-Vz^{MC})/Vz^{MC}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {100, -0.1, 0.1}});
    hists.mcTrackOrigin = histos.add<TH3>("MCQA/Hist3D_Nch_pT_origin", ";N_{ch};p_{T} (GeV/c);origin (0 primary, 1 secondary, 2 other MC collision, 3 no label)", kTH3F, {{nchAxisObs}, {KNbinsPtRes, KPtMin, KPtMax}, {kNOrigin, -0.5, static_cast<double>(kNOrigin) - 0.5}});
  }

  // vertex QA of all MC passes
  void declareMCVtxQA()
  {
    if (hists.mcVtx) {
      return;
    }
    hists.mcVtx = histos.add<TH2>("MCQA/Hist2D_dVz_type", ";z_{vtx}^{reco} - z_{vtx}^{MC} (cm);reco collision", kTH2F, {{200, -0.5, 0.5}, {kNVtxQA, -0.5, static_cast<double>(kNVtxQA) - 0.5}});
    labelAxis(hists.mcVtx->GetYaxis(), {"single", "split: used", "split: skipped"});
  }

  // (vz, sign, pT, eta, phi) maps, unweighted, efficiency-weighted and with all weights;
  // dir = "" (DATA) or "MCReco/".
  void declareFlat(uint32_t mask, const std::string& dir)
  {
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(mask, v) || hists.flat[v].maps) {
        continue;
      }
      hists.flat[v].maps = reg(v).add<THnSparse>((systPrefix(v) + dir + "hFlatMaps").c_str(), ";vz;sign;pt;eta;phi;weight", kTHnSparseF,
                                                 {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}, {kNFlatType, -0.5, static_cast<double>(kNFlatType) - 0.5}});
      labelAxis(hists.flat[v].maps->GetAxis(5), {"raw", "(1-fake)/eff", "all weights"});
    }
  }

  // Mean-pass profiles; MC: in MCGen/, MCReco/, MCRecoEffCorr/.
  void declareMean(uint32_t mask, bool isMC)
  {
    static const std::array<std::string, kNLevels> mcDir = {"MCGen/", "MCReco/", "MCRecoEffCorr/"};
    const int nLev = isMC ? kNLevels : 1;
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(mask, v)) {
        continue;
      }
      auto& r = reg(v);
      for (int l = 0; l < nLev; ++l) {
        const std::string p = systPrefix(v) + (isMC ? mcDir[l] : std::string{});
        auto& h = hists.mean[v][l];
        h.mean = r.add<TProfile3D>((p + "Prof_Mean_Nch_etabin_obs").c_str(), ";N_{ch};#eta-bin;observable (0 #LT[p_{T}]#GT, 1 #LT N#GT)", kTProfile3D, {{nchAxisOf(v)}, {etaBinAxis}, {kNMeanObs, -0.5, static_cast<double>(kNMeanObs) - 0.5}});
        if (cfgFwdOpt.cfgFwd) {
          h.fwdMean = r.add<TProfile2D>((p + "Prof_FwdAmpMean_Nch_fwdbin").c_str(), ";N_{ch};forward bin;#LT a#GT", kTProfile2D, {{nchAxisOf(v)}, {fwdAxis}});
        }
      }
    }
  }

  // Fluctuation-pass profiles named prefix + name + suffix, with x axis ax.
  // withGapSum: also the maps in (Delta eta, Sigma eta).
  FlucHists declareFluc(HistogramRegistry& r, const std::string& p, const std::string& s, const AxisSpec& ax, bool withGapSum)
  {
    FlucHists h;
    auto n = [&](const char* name) { return p + name + s; };
    h.obs = r.add<TProfile3D>(n("Prof_Obs_Nch_etabin_obs").c_str(), ";N_{ch};#eta-bin;observable (0 N, 1 [p_{T}], 2 c_{2}, 3 c_{3}, 4 C_{2}^{sub} mirror, 5 C_{3}^{sub} mirror)", kTProfile3D, {{ax}, {etaBinAxis}, {kNFlucObs, -0.5, static_cast<double>(kNFlucObs) - 0.5}});
    h.c2Sub2D = r.add<TProfile3D>(n("Prof_C2Sub2D_Nch_etaA_etaC").c_str(), ";N_{ch};#eta_{A};#eta_{C}", kTProfile3D, {{ax}, {etaAxis}, {etaAxis}});
    h.c3Sub2D = r.add<TProfile3D>(n("Prof_C3Sub2D_Nch_etaA_etaC").c_str(), ";N_{ch};#eta_{A} (pair);#eta_{C} (single)", kTProfile3D, {{ax}, {etaAxis}, {etaAxis}});
    if (withGapSum) {
      h.gapSum2D = r.add<TProfile3D>(n("Prof_GapSum2D_Nch_gap_sum").c_str(), ";N_{ch};#eta_{A} #minus #eta_{C} (Gap);#eta_{A} + #eta_{C} (Sum)", kTProfile3D, {{ax}, {gapAxis}, {sumAxis}});
      h.c3GapSum2D = r.add<TProfile3D>(n("Prof_C3GapSum2D_Nch_gap_sum").c_str(), ";N_{ch};#eta_{pair} #minus #eta_{single};#eta_{pair} + #eta_{single}", kTProfile3D, {{ax}, {gapAxis}, {sumAxis}});
    }
    if (cfgFwdOpt.cfgFwd) {
      h.fwdVar = r.add<TProfile2D>(n("Prof_FwdAmpVar_Nch_fwdbin").c_str(), ";N_{ch};forward bin;#LT #delta a^{2}#GT", kTProfile2D, {{ax}, {fwdAxis}});
      h.covFwdPt = r.add<TProfile3D>(n("Prof_CovFwdAmpMeanpT_Nch_etabin_fwdbin").c_str(), ";N_{ch};#eta-bin;forward bin;#LT #delta[p_{T}] #delta a#GT", kTProfile3D, {{ax}, {etaBinAxis}, {fwdAxis}});
      h.covFwdN = r.add<TProfile3D>(n("Prof_CovFwdAmpNch_Nch_etabin_fwdbin").c_str(), ";N_{ch};#eta-bin;forward bin;#LT #delta N #delta a#GT", kTProfile3D, {{ax}, {etaBinAxis}, {fwdAxis}});
    }
    return h;
  }

  // Gap/sum maps for every DATA variation; in MC only for Base, to limit the memory.
  void declareFlucAll(uint32_t mask, bool isMC)
  {
    static const std::array<std::string, kNLevels> mcDir = {"MCGen/", "MCReco/", "MCRecoEffCorr/"};
    const int nLev = isMC ? kNLevels : 1;
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(mask, v)) {
        continue;
      }
      for (int l = 0; l < nLev; ++l) {
        hists.fluc[v][l] = declareFluc(reg(v), systPrefix(v) + (isMC ? mcDir[l] : std::string{}), "", nchAxisOf(v), !isMC || v == kSystBase);
      }
    }
  }

  // Poisson-bootstrap replicas of all fluctuation profiles (Base, DataFluc only). Their x axis
  // is the Base centrality in 5% classes, the binning in which the errors are evaluated.
  void declareBootstrapHists()
  {
    for (int s = 0; s < nBoot; ++s) {
      hists.boot[s] = declareFluc(histosBoot, "Bootstrap/", Form("_sample%d", s), centAxisBoot, true);
    }
  }

  // Histograms of each pass, booked after the inputs so that switched-off variations get none.
  void declareMCGetEffHists()
  {
    declareNchCent(maskEff);
    declareEff(maskEff);
    declareMCQA();
  }

  void declareMCGetFlatHists()
  {
    declareNchCent(maskFlat);
    declareFlat(maskFlat, "MCReco/");
  }

  void declareMCMeanHists()
  {
    declareNchCent(maskAll);
    declareFlat(1, "MCReco/"); // Base map with all weights, as QA
    if (cfgFwdOpt.cfgFwd) {
      declareFwdQA();
    }
    declareMean(maskAll, true);
    histos.add("Prof_Weights_pT_type", ";p_{T};weight (0 eff, 1 fake, 2 w)", kTProfile2D, {{KNbinsPtRes, KPtMin, KPtMax}, {3, -0.5, 2.5}});
  }

  void declareMCFlucHists()
  {
    declareNchCent(maskAll);
    declareFlat(1, "MCReco/");
    if (cfgFwdOpt.cfgFwd) {
      declareFwdQA();
    }
    declareFlucAll(maskAll, true);
  }

  void declareDataGetFlatHists()
  {
    declareNchCent(maskFlat);
    declareFlat(maskFlat, "");
  }

  void declareDataMeanHists()
  {
    declareNchCent(maskAll);
    declareFlat(1, "");
    if (cfgFwdOpt.cfgFwd) {
      declareFwdQA();
    }
    declareMean(maskAll, false);
    histos.add("Prof_Weights_pt_eta_type", ";p_{T} (GeV/c);#eta;weight (0 eff, 1 fake, 2 w_{#phi})", kTProfile3D, {{KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}, {3, -0.5, 2.5}});
  }

  void declareDataFlucHists()
  {
    declareNchCent(maskAll);
    declareFlat(1, "");
    if (cfgFwdOpt.cfgFwd) {
      declareFwdQA();
    }
    declareFlucAll(maskAll, false);
    if (doBoot) {
      declareBootstrapHists();
    }
  }

  // ===========================================================================
  // CCDB helpers
  // ===========================================================================
  static void copyAxis(TAxis* dst, const TAxis* src)
  {
    const TArrayD* xb = src->GetXbins();
    if (xb->GetSize() > 0) {
      dst->Set(src->GetNbins(), xb->GetArray());
    } else {
      dst->Set(src->GetNbins(), src->GetXmin(), src->GetXmax());
    }
  }

  // Flattening weights w_phi = <N_phi> / N_phi from bin `type` of the last axis of hFlatMaps.
  THnSparseF* buildWeightMapFromRaw(const THnSparse* hRaw, int type, const char* mapName)
  {
    if (!hRaw || hRaw->GetNdimensions() != 6) {
      LOGF(error, "Raw eta-phi map for '%s' is null or not 6D; no flattening will be applied.", mapName);
      return nullptr;
    }
    std::array<int, 5> nb{};
    for (int i = 0; i < 5; ++i) {
      nb[i] = hRaw->GetAxis(i)->GetNbins();
    }
    auto* hWMap = new THnSparseF(mapName, Form("Flattening Weight Map %s (w_{#phi} = <N_{#phi}> / N_{#phi})", mapName), 5, nb.data());
    for (int i = 0; i < 5; ++i) {
      copyAxis(hWMap->GetAxis(i), hRaw->GetAxis(i));
    }

    std::array<int, 6> bins{};
    std::array<int, 5> wbins{};
    bins[5] = type + 1;
    for (int iv = 1; iv <= nb[0]; ++iv) {
      bins[0] = wbins[0] = iv;
      for (int ichg = 1; ichg <= nb[1]; ++ichg) {
        bins[1] = wbins[1] = ichg;
        for (int ipt = 1; ipt <= nb[2]; ++ipt) {
          bins[2] = wbins[2] = ipt;
          for (int ie = 1; ie <= nb[3]; ++ie) {
            bins[3] = wbins[3] = ie;
            double sum = 0.0;
            const int nphi = nb[4];
            for (int ip = 1; ip <= nphi; ++ip) {
              bins[4] = ip;
              sum += hRaw->GetBinContent(bins.data());
            }
            const double avg = (nphi > 0 ? sum / nphi : 0.0);
            for (int ip = 1; ip <= nphi; ++ip) {
              bins[4] = wbins[4] = ip;
              const double raw = hRaw->GetBinContent(bins.data());
              const double w = (avg > 0.0 && raw > 0.0) ? (avg / raw) : 1.0;
              hWMap->SetBinContent(wbins.data(), w);
            }
          }
        }
      }
    }
    LOGF(info, "Flattening weight map '%s' built.", mapName);
    return hWMap;
  }

  template <typename TP>
  void loadProfileFromList(TList* src, const char* name, TP*& target)
  {
    if (!src) {
      return;
    }
    auto* obj = src->FindObject(name);
    if (!obj) {
      LOGF(error, "Profile %s missing in CCDB TList", name);
      return;
    }
    auto* tp = dynamic_cast<TP*>(obj);
    if (!tp) {
      LOGF(error, "%s is not the expected profile type (it is %s)", name, obj->ClassName());
      return;
    }
    target = dynamic_cast<TP*>(tp->Clone());
    target->SetDirectory(nullptr);
    LOGF(info, "Loaded %s from list", name);
  }

  // Switch variation v off because an input is missing; a missing Base input is fatal.
  void dropSyst(int v, const std::string& why)
  {
    if (v == kSystBase) {
      LOGF(fatal, "Base input missing: %s", why.c_str());
      return;
    }
    LOGF(error, "Variation %s disabled: %s", systName[v].c_str(), why.c_str());
    CLRBIT(maskAll, v);
    CLRBIT(maskEff, v);
    CLRBIT(maskFlat, v);
  }

  // Efficiency and fake-rate maps of variation s versus (centrality, pT, eta):
  // eff = matched primaries / primaries, fake = (secondaries + other collision + no label) / all reco.
  bool loadEffMaps(int s, const std::string& path, int64_t now)
  {
    auto* lst = ccdb->getForTimeStamp<TList>(path, now);
    if (!lst) {
      return false;
    }
    std::array<TH3*, kNEffCat> h{};
    for (int c = 0; c < kNEffCat; ++c) {
      h[c] = dynamic_cast<TH3*>(lst->FindObject(effMapName[c].c_str()));
      if (!h[c]) {
        LOGF(error, "Missing efficiency map %s in %s", effMapName[c].c_str(), path.c_str());
        return false;
      }
    }
    auto* eff = dynamic_cast<TH3*>(h[kEffRecoMatchedPrim]->Clone(Form("hEff_%d", s)));
    auto* fake = dynamic_cast<TH3*>(h[kEffRecoSec]->Clone(Form("hFake_%d", s)));
    eff->SetDirectory(nullptr);
    fake->SetDirectory(nullptr);
    fake->Add(h[kEffRecoOtherColl]);
    fake->Add(h[kEffRecoFake]);
    // integrated over centrality (one centrality bin), for empty or invalid cells
    const int nCent = eff->GetNbinsX();
    auto ratio = [&](TH3* num, TH3* den, const char* name) {
      auto* r = num->RebinX(nCent, name);
      std::unique_ptr<TH3> d(den->RebinX(nCent, Form("%s_den", name)));
      r->SetDirectory(nullptr);
      r->Divide(d.get());
      return r;
    };
    inputs[s].hEffPtEta = ratio(eff, h[kEffAllPrim], Form("hEff_%d_ptEta", s));
    inputs[s].hFakePtEta = ratio(fake, h[kEffAllReco], Form("hFake_%d_ptEta", s));
    eff->Divide(h[kEffAllPrim]);
    fake->Divide(h[kEffAllReco]);
    inputs[s].hEff = eff;
    inputs[s].hFake = fake;
    return true;
  }

  THnSparseF* loadFlatMap(int s, const std::string& path, const std::string& dir, int64_t now)
  {
    auto* lst = ccdb->getForTimeStamp<TList>(path, now);
    if (!lst) {
      return nullptr;
    }
    const std::string hName = dir + "hFlatMaps";
    const auto* hRaw = dynamic_cast<const THnSparse*>(lst->FindObject(hName.c_str()));
    if (!hRaw) {
      LOGF(error, "Flattening map '%s' not found in %s", hName.c_str(), path.c_str());
      return nullptr;
    }
    return buildWeightMapFromRaw(hRaw, cfgEff ? kFlatWtd : kFlatRaw, Form("hFlatWeight_%d", s));
  }

  // Mean-pass objects of variation v for the fluctuation pass.
  void loadMeans(int v, const std::string& path, bool isMC, int64_t now)
  {
    static const std::array<std::string, kNLevels> mcDir = {"MCGen/", "MCReco/", "MCRecoEffCorr/"};
    auto* lst = ccdb->getForTimeStamp<TList>(path, now);
    if (!lst) {
      dropSyst(v, "no mean-pT list at " + path);
      return;
    }
    const int nLev = isMC ? kNLevels : 1;
    for (int l = 0; l < nLev; ++l) {
      const std::string d = isMC ? mcDir[l] : std::string{};
      loadProfileFromList(lst, (d + "Prof_Mean_Nch_etabin_obs").c_str(), inputs[v].pmean[l]);
      if (!inputs[v].pmean[l]) {
        dropSyst(v, "mean-pT profile missing in " + path);
        return;
      }
      // The fluctuation pass reads its reference at (N_ch bin, eta bin), so the binning of the
      // mean pass has to be the same.
      const auto* p = inputs[v].pmean[l];
      if (p->GetYaxis()->GetNbins() != nEta + 1 ||
          p->GetXaxis()->GetNbins() != nchAxisOf(v).getNbins() ||
          std::abs(p->GetXaxis()->GetBinWidth(1) - nchBinWidthOf(v)) > KFloatEpsilon) {
        LOGF(fatal,
             "Mean pass %s has %d N_ch bins of width %.1f and %d eta bins; this job uses %d bins of "
             "width %d and %d eta bins. Rerun the mean pass with the same binning.",
             path.c_str(), p->GetXaxis()->GetNbins(), p->GetXaxis()->GetBinWidth(1), p->GetYaxis()->GetNbins(),
             nchAxisOf(v).getNbins(), nchBinWidthOf(v), nEta + 1);
      }
      if (cfgFwdOpt.cfgFwd) {
        loadProfileFromList(lst, (d + "Prof_FwdAmpMean_Nch_fwdbin").c_str(), inputs[v].pmeanFwd[l]);
        if (inputs[v].pmeanFwd[l] && inputs[v].pmeanFwd[l]->GetYaxis()->GetNbins() != nFwd) {
          LOGF(error, "%s: mean pass has %d forward bins, this job %d; forward covariances not filled",
               path.c_str(), inputs[v].pmeanFwd[l]->GetYaxis()->GetNbins(), nFwd);
          inputs[v].pmeanFwd[l] = nullptr;
        }
      }
    }
  }

  // ===========================================================================
  // init
  // ===========================================================================
  void init(InitContext&)
  {
    TH1::SetDefaultSumw2(true);
    // Multiplicity axes: multQAAxis for the QA, nchAxisObs (nchAxisNchBin2) for the observables.
    const int nMax = (cfgSys == kPbPb) ? cfgNch.cfgNchPbMax.value : (cfgSys == kpp ? cfgNch.cfgNchPPMax.value : cfgNch.cfgNchOMax.value);
    if (cfgSys == kPbPb) {
      multQAAxis = {nMax / 2, KBinOffset, nMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    } else {
      multQAAxis = {nMax, KBinOffset, nMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    }
    {
      // bin k holds k*width ... (k+1)*width - 1 tracks
      auto makeNchAxis = [&](int width) {
        if (width < 1) {
          LOGF(fatal, "N_ch bin width must be >= 1 (got %d)", width);
        }
        const int nb = (nMax + width - 1) / width;
        return AxisSpec{nb, -KBinOffset, nb * width - KBinOffset, "N_{ch}"};
      };
      const bool pp = (cfgSys == kpp);
      nchWidth = (cfgSys == kPbPb) ? cfgNch.cfgNchBinWidthPbPb.value : (pp ? cfgNch.cfgNchBinWidthPP.value : cfgNch.cfgNchBinWidthOO.value);
      nchWidthNchBin2 = (cfgSys == kPbPb) ? KSystNchBinWidthPbPb : (pp ? KSystNchBinWidthPP : KSystNchBinWidthOO);
      nchAxisObs = makeNchAxis(nchWidth);
      nchAxisNchBin2 = makeNchAxis(nchWidthNchBin2);
    }

    // ---- run type: MC is pinned to the reference binning, DATA is configurable ----
    // cfgEtaBinWidth exists only to vary the DATA measurement. MC (efficiency,
    // closure) always uses 0.2-wide bins so its mean and fluctuation passes stay
    // mutually consistent and reproducible.
    const bool isMcRun = (cfgRunGetEff || cfgRunGetMCFlat || cfgRunMCMean || cfgRunMCFluc);
    const bool isDataRun = (cfgRunGetDataFlat || cfgRunDataMean || cfgRunDataFluc);
    if (isMcRun && isDataRun) {
      LOGF(fatal,
           "MC and DATA process switches are both enabled in one job: MC binning is "
           "pinned to 0.2 while DATA uses cfgEtaBinWidth. Run them as separate jobs.");
    }
    const float mcEtaBinWidth = 0.2f;
    const float effEtaBinWidth = isMcRun ? mcEtaBinWidth : static_cast<float>(cfgEtaBinWidth);

    // ---- observable eta binning (width = effEtaBinWidth) --------------------------
    {
      const float halfEta = cfgCutEta;
      const float width = effEtaBinWidth;
      if (width <= 0.f) {
        LOGF(fatal, "eta bin width must be > 0 (got %.3f)", width);
      }
      const int nHalf = static_cast<int>(std::lround(halfEta / width));
      if (nHalf < 1 || std::abs(nHalf * width - halfEta) > KEtaEdgeTolerance) {
        LOGF(fatal,
             "cfgCutEta=%.3f is not an integer multiple of the eta bin width=%.3f; "
             "eta=0 must be a bin edge for mirror pairing.",
             halfEta, width);
      }
      const int nbins = 2 * nHalf; // even -> symmetric about eta = 0
      const float lo = -halfEta;
      const float hi = halfEta;

      etaLw.clear();
      etaUp.clear();
      etaLw.push_back(lo); // index 0: full-range reference bin
      etaUp.push_back(hi);
      for (int i = 0; i < nbins; ++i) {
        etaLw.push_back(lo + i * width);
        etaUp.push_back(lo + (i + 1) * width);
      }
      nEta = nbins + 1;
      if (nEta > KNEtaMax) {
        LOGF(fatal, "nEta=%d exceeds KNEtaMax=%d (raise the cap or widen the eta bin width)", nEta, KNEtaMax);
      }

      std::vector<double> obsEdges;
      obsEdges.reserve(nEta);
      obsEdges.push_back(etaLw[1]);
      for (int i = 1; i < nEta; ++i) {
        obsEdges.push_back(etaUp[i]);
      }
      etaAxis = AxisSpec{obsEdges, "#eta"};
      etaBinAxis = AxisSpec{nEta + 1, -0.5, static_cast<double>(nEta) + 0.5, "#eta bin Number"};
      // (Delta eta, Sigma eta): fills are integer multiples of w, so 2N-1 bins
      // of width w centred on them, out to (N-1/2) w.
      const int nGS = 2 * nbins - 1;
      const double gsMax = (nbins - 0.5) * width;
      gapAxis = AxisSpec{nGS, -gsMax, gsMax, "#Delta#eta"};
      sumAxis = AxisSpec{nGS, -gsMax, gsMax, "#Sigma#eta"};
    }

    // ---- forward bins: detector totals, then the eta sub-bins of each detector ----
    {
      static const std::array<const char*, kNFwdDet> detName = {"FT0A", "FT0C", "FV0A"};
      fwdEdges = {cfgFwdOpt.cfgFwdEtaEdgesFT0A.value, cfgFwdOpt.cfgFwdEtaEdgesFT0C.value, cfgFwdOpt.cfgFwdEtaEdgesFV0A.value};
      fwdLabels.assign(detName.begin(), detName.end());
      nFwd = kNFwdDet;
      for (int det = 0; det < kNFwdDet; ++det) {
        auto& e = fwdEdges[det];
        if (e.size() < 2 || !std::is_sorted(e.begin(), e.end())) {
          LOGF(fatal, "%s eta edges must be >= 2 increasing values", detName[det]);
        }
        fwdOffset[det] = nFwd;
        for (std::size_t i = 0; i + 1 < e.size(); ++i) {
          fwdLabels.push_back(Form("%s %.1f<#eta<%.1f", detName[det], e[i], e[i + 1]));
          ++nFwd;
        }
      }
      if (nFwd > KNFwdMax) {
        LOGF(fatal, "%d forward bins exceed KNFwdMax=%d", nFwd, KNFwdMax);
      }
      fwdAxis = AxisSpec{nFwd, -0.5, nFwd - 0.5, "forward bin"};
      std::string list;
      for (int k = 0; k < nFwd; ++k) {
        list += Form("%s%d %s", k > 0 ? ", " : "", k, fwdLabels[k].c_str());
      }
      LOGF(info, "RadialFlowDecorr forward bins: %s", list.c_str());
    }

    // ---- variations ----
    const uint32_t active = (static_cast<uint32_t>(cfgSyst.cfgSystMask.value) | 1) & KAllSystMask;
    maskAll = maskEff = maskFlat = 0;
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(active, v)) {
        continue;
      }
      SETBIT(maskAll, v);
      if (effSource(v) == v) {
        SETBIT(maskEff, v);
      }
      if (flatSource(v) == v) {
        SETBIT(maskFlat, v);
      }
    }
    const bool isPP = (cfgSys == kpp);
    itsLayerNom = cfgEvSel.cfgITSLayerCut >= 0 ? cfgEvSel.cfgITSLayerCut.value : (isPP ? kITSLayer0123 : kITSLayersAll);
    if (itsLayerNom == kITSLayersAll) {
      itsLayerVar = kITSLayer0123;
    } else if (itsLayerNom == kITSLayer0123) {
      itsLayerVar = kITSLayer3;
    } else {
      itsLayerVar = kITSLayersAll;
    }
    occupancyMax = (cfgSys == kPbPb) ? cfgSyst.cfgSystOccMaxPbPb.value : (isPP ? cfgSyst.cfgSystOccMaxPP.value : cfgSyst.cfgSystOccMaxOO.value);
    pupNSigma = (cfgSys == kPbPb) ? cfgPileup.cfgPupnSigPbPb.value : (isPP ? cfgPileup.cfgPupnSigPP.value : cfgPileup.cfgPupnSigOO.value);
    if (cfgEvSel.cfgCentralityChoice == KSystCentralityChoice && TESTBIT(maskAll, kSystCent)) {
      LOGF(warning, "systCent uses the same centrality estimator as Base (%d).", cfgEvSel.cfgCentralityChoice.value);
    }

    // bootstrap: Base data fluctuation pass only
    doBoot = cfgRunDataFluc && (cfgNBootstrap > 0);
    nBoot = std::min<int>(cfgNBootstrap, KMaxBoot);
    if (cfgNBootstrap > KMaxBoot) {
      LOGF(warning, "cfgNBootstrap=%d capped at %d", cfgNBootstrap.value, KMaxBoot);
    }
    rng.SetSeed(cfgBootstrapSeed);

    // ---- CCDB ----
    ccdb->setURL(cfgCCDBurl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
    // a missing object returns nullptr: every lookup is checked (a missing variation input switches
    // the variation off, a missing Base input or alignment is fatal, a missing gain gives unit gains)
    ccdb->setFatalWhenNull(false);
    fwdCalib.fv0Det = o2::fv0::Geometry::instance(o2::fv0::Geometry::eUninitialized);

    std::string sysDir;
    switch (cfgSys) {
      case kPbPb:
        sysDir = "PbPbTest";
        break;
      case kNeNe:
        sysDir = "NeNeTest";
        break;
      case kOO:
        sysDir = "OOTest";
        break;
      case kpp:
        sysDir = "ppTest";
        break;
      default:
        LOGF(fatal, "Invalid cfgSys value: %d", cfgSys.value);
    }
    const std::string pathEff = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_EffMaps";
    const std::string pathMCFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_MCFlatMaps";
    const std::string pathMCMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_MCMean";
    const std::string pathDataFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_DataFlatMaps";
    const std::string pathDataMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_DataMean";

    // ---- inputs of each variation; a variation without its inputs is switched off ----
    // variations used by the enabled passes
    uint32_t used = 0;
    if (cfgRunGetEff) {
      used |= maskEff;
    }
    if (cfgRunGetMCFlat || cfgRunGetDataFlat) {
      used |= maskFlat;
    }
    if (cfgRunMCMean || cfgRunMCFluc || cfgRunDataMean || cfgRunDataFluc) {
      used |= maskAll;
    }

    // efficiency/fake maps
    if (cfgEff && (cfgRunGetMCFlat || cfgRunMCMean || cfgRunMCFluc || isDataRun)) {
      std::array<int8_t, kNSyst> loaded{};
      for (int v = 0; v < kNSyst; ++v) {
        if (!TESTBIT(used, v)) {
          continue;
        }
        const int s = effSource(v);
        if (loaded[s] == 0) {
          loaded[s] = loadEffMaps(s, pathEff + systSuffix(s), now) ? 1 : -1;
        }
        if (loaded[s] < 0) {
          dropSyst(v, "no efficiency maps at " + pathEff + systSuffix(s));
          continue;
        }
        inputs[v].hEff = inputs[s].hEff;
        inputs[v].hFake = inputs[s].hFake;
        inputs[v].hEffPtEta = inputs[s].hEffPtEta;
        inputs[v].hFakePtEta = inputs[s].hFakePtEta;
      }
    }

    // flattening maps
    if (cfgFlat && (cfgRunMCMean || cfgRunMCFluc || cfgRunDataMean || cfgRunDataFluc)) {
      const bool mc = cfgRunMCMean || cfgRunMCFluc;
      const std::string& path = mc ? pathMCFlat : pathDataFlat;
      const std::string dir = mc ? "MCReco/" : "";
      std::array<int8_t, kNSyst> loaded{};
      for (int v = 0; v < kNSyst; ++v) {
        const int s = flatSource(v);
        if (!TESTBIT(maskAll, v) || s < 0) {
          continue;
        }
        if (loaded[s] == 0) {
          inputs[s].hFlat = loadFlatMap(s, path + systSuffix(s), dir, now);
          loaded[s] = inputs[s].hFlat ? 1 : -1;
        }
        if (loaded[s] < 0) {
          dropSyst(v, "no flattening map at " + path + systSuffix(s));
          continue;
        }
        inputs[v].hFlat = inputs[s].hFlat;
      }
    }

    // ---- n-sigma pile-up bands, from the Base flattening output of the same production ----
    // The flattening and efficiency passes only apply the linear cut: their maps define the bands.
    {
      const bool meanOrFluc = cfgRunMCMean || cfgRunMCFluc || cfgRunDataMean || cfgRunDataFluc;
      const bool wantCent = cfgPileup.cfgApplySigPupCut && meanOrFluc;
      const bool wantPV = cfgPileup.cfgApplySigPupCutPV && meanOrFluc;
      if ((wantCent || wantPV) && (cfgRunGetEff || cfgRunGetMCFlat || cfgRunGetDataFlat)) {
        LOGF(fatal, "The n-sigma pile-up bands come from the flattening pass: run it in a separate job before the mean/fluctuation passes.");
      }
      if (wantCent || wantPV) {
        const std::string& bandPath = isMcRun ? pathMCFlat : pathDataFlat;
        auto* bandList = ccdb->getForTimeStamp<TList>(bandPath, now);
        if (!bandList) {
          LOGF(fatal, "n-sigma pile-up bands requested but no flattening output at %s", bandPath.c_str());
        }
        auto makeBand = [&](const char* name, NSigmaBand& band) {
          auto* h2 = dynamic_cast<TH2*>(bandList->FindObject(name));
          if (!h2) {
            LOGF(fatal, "n-sigma pile-up band: %s missing in %s", name, bandPath.c_str());
            return;
          }
          std::unique_ptr<TProfile> prof(h2->ProfileX(Form("%s_band", name), 1, -1, "S"));
          const int nBins = prof->GetNbinsX();
          band.xMin = prof->GetXaxis()->GetXmin();
          band.xMax = prof->GetXaxis()->GetXmax();
          band.meanRms.assign(nBins + 2, {0.f, -1.f});
          for (int i = 1; i <= nBins; ++i) {
            // a sparse bin (high-N_PV tail, peripheral edge) gets no band: its RMS is not measured
            if (prof->GetBinEntries(i) >= KMinBandEntries) {
              band.meanRms[i] = {static_cast<float>(prof->GetBinContent(i)), static_cast<float>(prof->GetBinError(i))};
            }
          }
          band.applied = true;
          LOGF(info, "n-sigma pile-up band (%.1f sigma, PileupNSigma %.1f) from %s in %s", pupNSigma,
               pupNSigma - KSystPupnSigStep, name, bandPath.c_str());
        };
        if (wantCent) {
          makeBand("Hist2D_globalTracks_cent", pileupBands[kBandCent]);
        }
        if (wantPV) {
          makeBand("Hist2D_globalTracks_PVTracks", pileupBands[kBandPV]);
        }
      }
    }

    // ---- mean profiles for the fluctuation passes (one per variation) ----
    if (cfgRunMCFluc || cfgRunDataFluc) {
      const std::string& path = cfgRunMCFluc ? pathMCMean : pathDataMean;
      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(maskAll, v)) {
          LOGF(info, "Loading mean profiles of %s from %s", systName[v].c_str(), (path + systSuffix(v)).c_str());
          loadMeans(v, path + systSuffix(v), cfgRunMCFluc, now);
        }
      }
    }

    // ---- declarations ----
    declareCommonQA();
    if (isMcRun) {
      declareMCVtxQA();
    }
    declareCounts(cfgRunGetEff ? maskEff : ((cfgRunGetMCFlat || cfgRunGetDataFlat) ? maskFlat : maskAll));
    if (doprocessCutQA) {
      declareCutQA();
    }
    if (cfgRunGetEff) {
      declareMCGetEffHists();
    }
    if (cfgRunGetMCFlat) {
      declareMCGetFlatHists();
    }
    if (cfgRunMCMean) {
      declareMCMeanHists();
    }
    if (cfgRunMCFluc) {
      declareMCFlucHists();
    }
    if (cfgRunGetDataFlat) {
      declareDataGetFlatHists();
    }
    if (cfgRunDataMean) {
      declareDataMeanHists();
    }
    if (cfgRunDataFluc) {
      declareDataFlucHists();
    }

    std::string list;
    for (int v = 0; v < kNSyst; ++v) {
      if (TESTBIT(used & maskAll, v) || (cfgRunGetEff && TESTBIT(maskEff, v)) || ((cfgRunGetMCFlat || cfgRunGetDataFlat) && TESTBIT(maskFlat, v))) {
        list += " " + systName[v];
      }
    }
    LOGF(info, "RadialFlowDecorr variations in this job:%s (ITS layer cut nominal %d, variation %d; N_ch bins of %d, NchBin2 %d; occupancy cut %d)",
         list.c_str(), itsLayerNom, itsLayerVar, nchWidth, nchWidthNchBin2, occupancyMax);
    LOGF(info, "CCDB initialization complete for RadialFlowDecorr.");
  }

  // ===========================================================================
  // Helpers of the process functions
  // ===========================================================================
  // (vz, sign, pT, eta, phi) map of variation v: raw, with all weights w and, for Base, with
  // the efficiency weight wEff = (1 - fake) / eff.
  void fillFlat(int v, float vz, float sign, float pt, float eta, float phi, float w, float wEff)
  {
    const auto& h = hists.flat[v].maps;
    if (!h) {
      return;
    }
    std::array<double, 6> x = {vz, sign, pt, eta, phi, static_cast<double>(kFlatRaw)};
    h->Fill(x.data(), 1.);
    x[5] = kFlatWtd;
    h->Fill(x.data(), w);
    if (v == kSystBase) {
      x[5] = kFlatEffWtd;
      h->Fill(x.data(), wEff);
    }
  }

  // The reco collision with the most PV contributors; with split vertices the others are
  // skipped, so the truth of an MC collision enters once. Fills the vertex QA.
  template <typename R>
  int64_t bestRecoCollision(const R& recoCollisions, float mcVz)
  {
    int64_t best = -1;
    int maxContrib = -1;
    int nReco = 0;
    for (const auto& rc : recoCollisions) {
      ++nReco;
      if (rc.numContrib() > maxContrib) {
        maxContrib = rc.numContrib();
        best = rc.globalIndex();
      }
    }
    for (const auto& rc : recoCollisions) {
      const int type = (nReco == 1) ? kVtxSingle : (rc.globalIndex() == best ? kVtxSplitBest : kVtxSplitOther);
      hists.mcVtx->Fill(rc.posZ() - mcVz, type);
    }
    return best;
  }

  // ===========================================================================
  // MC: build efficiency / fake maps (inclusive)
  // ===========================================================================
  void processGetEffHists(aod::McCollisions::iterator const& mcCollision, McRecoCollisions const& recoCollisions,
                          FilteredTCs const& mcTracks, aod::McParticles const& mcParticles)
  {
    const int64_t best = bestRecoCollision(recoCollisions, mcCollision.posZ());
    for (const auto& coll : recoCollisions) {
      if (coll.globalIndex() != best) {
        continue;
      }
      const auto tracks = mcTracks.sliceBy(tracksPerCollision, coll.globalIndex());
      histos.fill(HIST("hVtxZ"), coll.posZ());
      const auto ev = isEventSelected(coll, tracks.size(), maskEff);
      if (ev.mask == 0) {
        return;
      }
      const float nCh = tracks.size();
      const float multPV = coll.multNTracksPV();
      const float vz = coll.posZ();
      const float cent = ev.cent[kSystBase];

      if (TESTBIT(ev.mask, kSystBase)) {
        histos.fill(HIST("hVtxZ_after_sel"), vz);
        histos.fill(HIST("hCentrality"), ev.cent[kSystBase]);
        histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, nCh);
        histos.fill(HIST("Hist2D_PVTracks_cent"), ev.cent[kSystBase], multPV);
        hists.nchOccupancy->Fill(nCh, ev.occupancy);
      }
      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(ev.mask, v) && hists.nchCent[v]) {
          hists.nchCent[v]->Fill(ev.cent[v], nCh);
        }
      }

      // one entry in the efficiency map of category cat of every variation in mask
      auto fillEffMaps = [&](uint32_t mask, float pt, float eta, int cat) {
        for (int v = 0; v < kNSyst; ++v) {
          if (TESTBIT(mask, v)) {
            hists.eff[v].maps[cat]->Fill(cent, pt, eta);
          }
        }
      };

      // --- truth: charged primaries of the MC collision ---
      for (const auto& particle : mcParticles) {
        if (!isParticleSelected(particle) || !particle.isPhysicalPrimary()) {
          continue;
        }
        fillEffMaps(ev.mask, particle.pt(), particle.eta(), kEffAllPrim);
      }

      // --- reco ---
      for (const auto& track : tracks) {
        const uint32_t trackMask = isTrackSelected(track, ev.mask);
        if (trackMask == 0) {
          continue;
        }
        float pt = track.pt(), eta = track.eta();
        const bool isBase = TESTBIT(trackMask, kSystBase);
        const int origin = trackOrigin(track, coll);
        if (isBase) {
          histos.fill(HIST("hPt"), pt);
          histos.fill(HIST("hEta"), eta);
          hists.mcTrackOrigin->Fill(nCh, pt, origin);
        }
        fillEffMaps(trackMask, pt, eta, kEffAllReco);
        if (origin == kOriginPrimary) {
          const auto& mcP = track.mcParticle();
          if (isBase) {
            histos.fill(HIST("ptResolution"), mcP.pt(), (pt - mcP.pt()) / mcP.pt());
            histos.fill(HIST("etaResolution"), mcP.eta(), eta - mcP.eta());
            histos.fill(HIST("etaTruthReco"), mcP.eta(), eta);
            histos.fill(HIST("vzResolution"), mcP.vz(), (vz - mcP.vz()) / mcP.vz());
            histos.fill(HIST("TruthTracKVz"), mcP.vz(), vz);
          }
          fillEffMaps(trackMask, mcP.pt(), mcP.eta(), kEffRecoMatchedPrim);
        } else if (origin == kOriginSecondary) {
          fillEffMaps(trackMask, pt, eta, kEffRecoSec);
        } else if (origin == kOriginOtherColl) {
          fillEffMaps(trackMask, pt, eta, kEffRecoOtherColl);
        } else {
          fillEffMaps(trackMask, pt, eta, kEffRecoFake);
        }
      }
      flushTrackCounts();
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetEffHists, "process MC to calculate EffWeights", cfgRunGetEff);

  // ===========================================================================
  // MC: build flattening maps (inclusive)
  // ===========================================================================
  void processMCFlat(aod::McCollisions::iterator const& mcCollision, McRecoCollisions const& recoCollisions, FilteredTCs const& mcTracks)
  {
    const int64_t best = bestRecoCollision(recoCollisions, mcCollision.posZ());
    for (const auto& coll : recoCollisions) {
      if (coll.globalIndex() != best) {
        continue;
      }
      const auto tracks = mcTracks.sliceBy(tracksPerCollision, coll.globalIndex());
      histos.fill(HIST("hVtxZ"), coll.posZ());
      const auto ev = isEventSelected(coll, tracks.size(), maskFlat);
      if (ev.mask == 0) {
        return;
      }
      const float nCh = tracks.size();
      const float multPV = coll.multNTracksPV();
      const float vz = coll.posZ();

      if (TESTBIT(ev.mask, kSystBase)) {
        histos.fill(HIST("hVtxZ_after_sel"), vz);
        histos.fill(HIST("hCentrality"), ev.cent[kSystBase]);
        histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, nCh);
        histos.fill(HIST("Hist2D_PVTracks_cent"), ev.cent[kSystBase], multPV);
        hists.nchOccupancy->Fill(nCh, ev.occupancy);
      }
      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(ev.mask, v) && hists.nchCent[v]) {
          hists.nchCent[v]->Fill(ev.cent[v], nCh);
        }
      }

      std::array<float, kNSyst> w{};
      for (const auto& track : tracks) {
        const uint32_t trackMask = isTrackSelected(track, ev.mask);
        if (trackMask == 0) {
          continue;
        }
        float pt = track.pt();
        if (pt <= cfgPtMin || pt > cfgPtMax) {
          continue;
        }
        float eta = track.eta(), phi = track.phi();
        float sign = track.sign();
        if (TESTBIT(trackMask, kSystBase)) {
          histos.fill(HIST("hPt"), pt);
          histos.fill(HIST("hEta"), eta);
          histos.fill(HIST("hPhi"), phi);
        }

        // efficiency weights only: the flattening map is what this pass measures
        float eff = 1.f, fake = 0.f, flatWeight = 1.f;
        const uint32_t weightMask = getTrackWeights(trackMask, ev.cent[kSystBase], vz, sign, pt, eta, phi, false, w, eff, fake, flatWeight);
        for (int v = 0; v < kNSyst; ++v) {
          if (TESTBIT(weightMask, v)) {
            fillFlat(v, vz, sign, pt, eta, phi, w[v], w[v]);
          }
        }
      }
      flushTrackCounts();
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFlat, "process MC to calculate FlatWeights", cfgRunGetMCFlat);

  // ===========================================================================
  // MC: mean pT (truth / reco / reco-eff-corrected)
  // ===========================================================================
  void processMCMean(aod::McCollisions::iterator const& mcCollision, McRecoCollisions const& recoCollisions,
                     aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&,
                     FilteredTCs const& mcTracks, aod::McParticles const& mcParticles)
  {
    const int64_t best = bestRecoCollision(recoCollisions, mcCollision.posZ());
    for (const auto& coll : recoCollisions) {
      if (coll.globalIndex() != best) {
        continue;
      }
      const auto tracks = mcTracks.sliceBy(tracksPerCollision, coll.globalIndex());
      histos.fill(HIST("hVtxZ"), coll.posZ());
      const auto ev = isEventSelected(coll, tracks.size(), maskAll);
      if (ev.mask == 0) {
        return;
      }
      const float nCh = tracks.size();
      const float multPV = coll.multNTracksPV();
      const float vz = coll.posZ();

      if (TESTBIT(ev.mask, kSystBase)) {
        histos.fill(HIST("hVtxZ_after_sel"), vz);
        histos.fill(HIST("hCentrality"), ev.cent[kSystBase]);
        histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, nCh);
        histos.fill(HIST("Hist2D_PVTracks_cent"), ev.cent[kSystBase], multPV);
        hists.nchOccupancy->Fill(nCh, ev.occupancy);
      }
      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(ev.mask, v) && hists.nchCent[v]) {
          hists.nchCent[v]->Fill(ev.cent[v], nCh);
        }
      }

      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(ev.mask, v)) {
          for (auto& s : sums[v]) {
            s.reset();
          }
        }
      }

      // --- truth: independent of the track cuts, copied to every selected variation ---
      PowerSums& sumsTru = sums[kSystBase][kGen];
      if (!TESTBIT(ev.mask, kSystBase)) {
        sumsTru.reset();
      }
      for (const auto& particle : mcParticles) {
        if (!isParticleSelected(particle) || !particle.isPhysicalPrimary()) {
          continue;
        }
        float pt = particle.pt(), eta = particle.eta();
        if (pt <= cfgPtMin || pt > cfgPtMax) {
          continue;
        }
        for (int ieta = 0; ieta < nEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
            continue;
          }
          sumsTru.add(ieta, pt, 1.);
        }
      }
      for (int v = 1; v < kNSyst; ++v) {
        if (TESTBIT(ev.mask, v)) {
          sums[v][kGen] = sumsTru;
        }
      }

      // --- reco ---
      std::array<float, kNSyst> w{};
      for (const auto& track : tracks) {
        const uint32_t trackMask = isTrackSelected(track, ev.mask);
        if (trackMask == 0) {
          continue;
        }
        float pt = track.pt();
        if (pt <= cfgPtMin || pt > cfgPtMax) {
          continue;
        }
        float eta = track.eta(), phi = track.phi();
        float sign = track.sign();
        if (TESTBIT(trackMask, kSystBase)) {
          histos.fill(HIST("hPt"), pt);
          histos.fill(HIST("hEta"), eta);
          histos.fill(HIST("hPhi"), phi);
        }

        // reco level: every selected track, independent of the corrections
        for (int ieta = 0; ieta < nEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
            continue;
          }
          for (int v = 0; v < kNSyst; ++v) {
            if (TESTBIT(trackMask, v)) {
              sums[v][kReco].add(ieta, pt, 1.);
            }
          }
        }

        float eff = 1.f, fake = 0.f, flatWeight = 1.f;
        const uint32_t weightMask = getTrackWeights(trackMask, ev.cent[kSystBase], vz, sign, pt, eta, phi, true, w, eff, fake, flatWeight);
        if (weightMask == 0) {
          continue;
        }
        if (TESTBIT(weightMask, kSystBase)) {
          histos.fill(HIST("Prof_Weights_pT_type"), pt, 0, eff);
          histos.fill(HIST("Prof_Weights_pT_type"), pt, 1, fake);
          histos.fill(HIST("Prof_Weights_pT_type"), pt, 2, w[kSystBase]);
          fillFlat(kSystBase, vz, sign, pt, eta, phi, w[kSystBase], (1.f - fake) / eff);
        }

        // reco, efficiency-corrected
        for (int ieta = 0; ieta < nEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
            continue;
          }
          for (int v = 0; v < kNSyst; ++v) {
            if (TESTBIT(weightMask, v)) {
              sums[v][kRecoEffCorr].add(ieta, pt, w[v]);
            }
          }
        }
      }
      flushTrackCounts();

      // --- forward detectors: charged primaries (truth) and amplitudes (reco) ---
      FwdSignal fwdTru, fwdReco;
      if (cfgFwdOpt.cfgFwd) {
        loadFwdCalib(coll.template bc_as<aod::BCsWithTimestamps>(), false);
        computeFwdTruth(mcParticles, fwdTru);
        computeFwdReco(coll, vz, fwdReco, TESTBIT(ev.mask, kSystBase));
      }

      // --- <[pT]> and <N> per eta bin, <a> per forward bin ---
      for (int v = 0; v < kNSyst; ++v) {
        if (!TESTBIT(ev.mask, v)) {
          continue;
        }
        for (int l = 0; l < kNLevels; ++l) {
          const auto& s = sums[v][l];
          const auto& h = hists.mean[v][l];
          for (int ieta = 0; ieta < nEta; ++ieta) {
            const double n = s.wk[ieta][1];
            h.mean->Fill(nCh, ieta, kMeanN, n);
            if (n <= KFloatEpsilon) {
              continue;
            }
            const double mpt = s.pmwk[ieta][1][1] / n;
            if (std::isfinite(mpt)) {
              h.mean->Fill(nCh, ieta, kMeanPt, mpt);
            }
          }
          const FwdSignal& fwd = (l == kGen) ? fwdTru : fwdReco;
          if (cfgFwdOpt.cfgFwd && h.fwdMean) {
            for (int k = 0; k < nFwd; ++k) {
              if (fwd.valid[k]) {
                h.fwdMean->Fill(nCh, k, fwd.ampl[k]);
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCMean, "process MC to calculate mean pt", cfgRunMCMean);

  // ===========================================================================
  // MC: fluctuations (C2, subevent, forward covariances) at three levels
  // ===========================================================================
  void processMCFluc(aod::McCollisions::iterator const& mcCollision, McRecoCollisions const& recoCollisions,
                     aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&,
                     FilteredTCs const& mcTracks, aod::McParticles const& mcParticles)
  {
    const int64_t best = bestRecoCollision(recoCollisions, mcCollision.posZ());
    for (const auto& coll : recoCollisions) {
      if (coll.globalIndex() != best) {
        continue;
      }
      const auto tracks = mcTracks.sliceBy(tracksPerCollision, coll.globalIndex());
      histos.fill(HIST("hVtxZ"), coll.posZ());
      const auto ev = isEventSelected(coll, tracks.size(), maskAll);
      if (ev.mask == 0) {
        return;
      }
      const float nCh = tracks.size();
      const float multPV = coll.multNTracksPV();
      const float vz = coll.posZ();

      if (TESTBIT(ev.mask, kSystBase)) {
        histos.fill(HIST("hVtxZ_after_sel"), vz);
        histos.fill(HIST("hCentrality"), ev.cent[kSystBase]);
        histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, nCh);
        histos.fill(HIST("Hist2D_PVTracks_cent"), ev.cent[kSystBase], multPV);
        hists.nchOccupancy->Fill(nCh, ev.occupancy);
      }
      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(ev.mask, v) && hists.nchCent[v]) {
          hists.nchCent[v]->Fill(ev.cent[v], nCh);
        }
      }

      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(ev.mask, v)) {
          for (auto& s : sums[v]) {
            s.reset();
          }
        }
      }

      // --- truth sums: independent of the track cuts, copied to every selected variation ---
      PowerSums& sumsTru = sums[kSystBase][kGen];
      if (!TESTBIT(ev.mask, kSystBase)) {
        sumsTru.reset();
      }
      for (const auto& particle : mcParticles) {
        if (!isParticleSelected(particle) || !particle.isPhysicalPrimary()) {
          continue;
        }
        float pt = particle.pt(), eta = particle.eta();
        if (pt <= cfgPtMin || pt > cfgPtMax) {
          continue;
        }
        for (int ieta = 0; ieta < nEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
            continue;
          }
          sumsTru.add(ieta, pt, 1.);
        }
      }
      for (int v = 1; v < kNSyst; ++v) {
        if (TESTBIT(ev.mask, v)) {
          sums[v][kGen] = sumsTru;
        }
      }

      // --- reco sums ---
      std::array<float, kNSyst> w{};
      for (const auto& track : tracks) {
        const uint32_t trackMask = isTrackSelected(track, ev.mask);
        if (trackMask == 0) {
          continue;
        }
        float pt = track.pt();
        if (pt <= cfgPtMin || pt > cfgPtMax) {
          continue;
        }
        float eta = track.eta(), phi = track.phi();
        float sign = track.sign();
        if (TESTBIT(trackMask, kSystBase)) {
          histos.fill(HIST("hPt"), pt);
          histos.fill(HIST("hEta"), eta);
          histos.fill(HIST("hPhi"), phi);
        }

        // reco level: every selected track, independent of the corrections
        for (int ieta = 0; ieta < nEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
            continue;
          }
          for (int v = 0; v < kNSyst; ++v) {
            if (TESTBIT(trackMask, v)) {
              sums[v][kReco].add(ieta, pt, 1.);
            }
          }
        }

        float eff = 1.f, fake = 0.f, flatWeight = 1.f;
        const uint32_t weightMask = getTrackWeights(trackMask, ev.cent[kSystBase], vz, sign, pt, eta, phi, true, w, eff, fake, flatWeight);
        if (weightMask == 0) {
          continue;
        }
        if (TESTBIT(weightMask, kSystBase)) {
          fillFlat(kSystBase, vz, sign, pt, eta, phi, w[kSystBase], (1.f - fake) / eff);
        }

        // reco, efficiency-corrected
        for (int ieta = 0; ieta < nEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
            continue;
          }
          for (int v = 0; v < kNSyst; ++v) {
            if (TESTBIT(weightMask, v)) {
              sums[v][kRecoEffCorr].add(ieta, pt, w[v]);
            }
          }
        }
      }
      flushTrackCounts();

      // --- forward detectors: charged primaries (truth) and amplitudes (reco) ---
      FwdSignal fwdTru, fwdReco;
      if (cfgFwdOpt.cfgFwd) {
        loadFwdCalib(coll.template bc_as<aod::BCsWithTimestamps>(), false);
        computeFwdTruth(mcParticles, fwdTru);
        computeFwdReco(coll, vz, fwdReco, TESTBIT(ev.mask, kSystBase));
      }

      // --- c2, c3, subevent correlators and forward covariances of every variation and level ---
      EventObs evObs;
      for (int v = 0; v < kNSyst; ++v) {
        if (!TESTBIT(ev.mask, v)) {
          continue;
        }
        for (int l = 0; l < kNLevels; ++l) {
          const FwdSignal* fwd = cfgFwdOpt.cfgFwd ? (l == kGen ? &fwdTru : &fwdReco) : nullptr;
          computeEventObs(sums[v][l], v, l, fwd, nCh, evObs);
          fillFluc(hists.fluc[v][l], evObs, 1., nCh);
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFluc, "process MC to calculate pt fluc", cfgRunMCFluc);

  // ===========================================================================
  // DATA: build flattening maps (inclusive)
  // ===========================================================================
  void processGetDataFlat(AodCollisionsSel::iterator const& coll, AodTracksSel const& tracks)
  {
    histos.fill(HIST("hVtxZ"), coll.posZ());
    const auto ev = isEventSelected(coll, tracks.size(), maskFlat);
    if (ev.mask == 0) {
      return;
    }
    const float nCh = tracks.size();
    const float multPV = coll.multNTracksPV();
    const float vz = coll.posZ();

    if (TESTBIT(ev.mask, kSystBase)) {
      histos.fill(HIST("hVtxZ_after_sel"), vz);
      histos.fill(HIST("hCentrality"), ev.cent[kSystBase]);
      histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, nCh);
      histos.fill(HIST("Hist2D_PVTracks_cent"), ev.cent[kSystBase], multPV);
      hists.nchOccupancy->Fill(nCh, ev.occupancy);
    }
    for (int v = 0; v < kNSyst; ++v) {
      if (TESTBIT(ev.mask, v) && hists.nchCent[v]) {
        hists.nchCent[v]->Fill(ev.cent[v], nCh);
      }
    }

    std::array<float, kNSyst> w{};
    for (const auto& track : tracks) {
      const uint32_t trackMask = isTrackSelected(track, ev.mask);
      if (trackMask == 0) {
        continue;
      }
      float pt = track.pt();
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }
      float eta = track.eta(), phi = track.phi();
      float sign = track.sign();
      if (TESTBIT(trackMask, kSystBase)) {
        histos.fill(HIST("hPt"), pt);
        histos.fill(HIST("hEta"), eta);
        histos.fill(HIST("hPhi"), phi);
      }

      // efficiency weights only: the flattening map is what this pass measures
      float eff = 1.f, fake = 0.f, flatWeight = 1.f;
      const uint32_t weightMask = getTrackWeights(trackMask, ev.cent[kSystBase], vz, sign, pt, eta, phi, false, w, eff, fake, flatWeight);
      for (int v = 0; v < kNSyst; ++v) {
        if (TESTBIT(weightMask, v)) {
          fillFlat(v, vz, sign, pt, eta, phi, w[v], w[v]);
        }
      }
    }
    flushTrackCounts();
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetDataFlat, "process data to calculate Flattening maps", cfgRunGetDataFlat);

  // ===========================================================================
  // DATA: mean pT
  // ===========================================================================
  void processDataMean(AodCollisionsSel::iterator const& coll, aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&, AodTracksSel const& tracks)
  {
    const auto ev = isEventSelected(coll, tracks.size(), maskAll);
    if (ev.mask == 0) {
      return;
    }
    const float nCh = tracks.size();
    const float multPV = coll.multNTracksPV();
    const float vz = coll.posZ();

    if (TESTBIT(ev.mask, kSystBase)) {
      histos.fill(HIST("hVtxZ_after_sel"), vz);
      histos.fill(HIST("hCentrality"), ev.cent[kSystBase]);
      histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, nCh);
      histos.fill(HIST("Hist2D_PVTracks_cent"), ev.cent[kSystBase], multPV);
      hists.nchOccupancy->Fill(nCh, ev.occupancy);
    }
    for (int v = 0; v < kNSyst; ++v) {
      if (TESTBIT(ev.mask, v) && hists.nchCent[v]) {
        hists.nchCent[v]->Fill(ev.cent[v], nCh);
      }
    }

    for (int v = 0; v < kNSyst; ++v) {
      if (TESTBIT(ev.mask, v)) {
        sums[v][0].reset();
      }
    }
    std::array<float, kNSyst> w{};
    for (const auto& track : tracks) {
      const uint32_t trackMask = isTrackSelected(track, ev.mask);
      if (trackMask == 0) {
        continue;
      }
      float p = track.p();
      float pt = track.pt();
      if (p < KFloatEpsilon) {
        continue;
      }
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }
      float eta = track.eta(), phi = track.phi();
      float sign = track.sign();
      if (TESTBIT(trackMask, kSystBase)) {
        histos.fill(HIST("hP"), p);
        histos.fill(HIST("hPt"), pt);
        histos.fill(HIST("hEta"), eta);
        histos.fill(HIST("hPhi"), phi);
      }

      float eff = 1.f, fake = 0.f, flatWeight = 1.f;
      const uint32_t weightMask = getTrackWeights(trackMask, ev.cent[kSystBase], vz, sign, pt, eta, phi, true, w, eff, fake, flatWeight);
      if (weightMask == 0) {
        continue;
      }
      if (TESTBIT(weightMask, kSystBase)) {
        histos.fill(HIST("Prof_Weights_pt_eta_type"), pt, eta, 0, eff);
        histos.fill(HIST("Prof_Weights_pt_eta_type"), pt, eta, 1, fake);
        histos.fill(HIST("Prof_Weights_pt_eta_type"), pt, eta, 2, flatWeight);
        fillFlat(kSystBase, vz, sign, pt, eta, phi, w[kSystBase], (1.f - fake) / eff);
      }

      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
          continue;
        }
        for (int v = 0; v < kNSyst; ++v) {
          if (TESTBIT(weightMask, v)) {
            sums[v][0].add(ieta, pt, w[v]);
          }
        }
      }
    }
    flushTrackCounts();

    // --- forward detectors: gain-equalised amplitudes ---
    FwdSignal fwd;
    if (cfgFwdOpt.cfgFwd) {
      loadFwdCalib(coll.template bc_as<aod::BCsWithTimestamps>(), true);
      computeFwdReco(coll, vz, fwd, TESTBIT(ev.mask, kSystBase));
    }

    // --- <[pT]> and <N> per eta bin, <a> per forward bin ---
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(ev.mask, v)) {
        continue;
      }
      const auto& s = sums[v][0];
      const auto& h = hists.mean[v][0];
      for (int ieta = 0; ieta < nEta; ++ieta) {
        const double n = s.wk[ieta][1];
        h.mean->Fill(nCh, ieta, kMeanN, n);
        if (n <= KFloatEpsilon) {
          continue;
        }
        const double mpt = s.pmwk[ieta][1][1] / n;
        if (std::isfinite(mpt)) {
          h.mean->Fill(nCh, ieta, kMeanPt, mpt);
        }
      }
      if (cfgFwdOpt.cfgFwd && h.fwdMean) {
        for (int k = 0; k < nFwd; ++k) {
          if (fwd.valid[k]) {
            h.fwdMean->Fill(nCh, k, fwd.ampl[k]);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataMean, "process data to calculate mean pT", cfgRunDataMean);

  // ===========================================================================
  // DATA: fluctuations (C2, subevent, forward covariances) + Poisson bootstrap for Base
  // ===========================================================================
  void processDataFluc(AodCollisionsSel::iterator const& coll, aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&, AodTracksSel const& tracks)
  {
    const auto ev = isEventSelected(coll, tracks.size(), maskAll);
    if (ev.mask == 0) {
      return;
    }
    const float nCh = tracks.size();
    const float multPV = coll.multNTracksPV();
    const float vz = coll.posZ();

    if (TESTBIT(ev.mask, kSystBase)) {
      histos.fill(HIST("hVtxZ_after_sel"), vz);
      histos.fill(HIST("hCentrality"), ev.cent[kSystBase]);
      histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, nCh);
      histos.fill(HIST("Hist2D_PVTracks_cent"), ev.cent[kSystBase], multPV);
      hists.nchOccupancy->Fill(nCh, ev.occupancy);
    }
    for (int v = 0; v < kNSyst; ++v) {
      if (TESTBIT(ev.mask, v) && hists.nchCent[v]) {
        hists.nchCent[v]->Fill(ev.cent[v], nCh);
      }
    }

    for (int v = 0; v < kNSyst; ++v) {
      if (TESTBIT(ev.mask, v)) {
        sums[v][0].reset();
      }
    }
    std::array<float, kNSyst> w{};
    for (const auto& track : tracks) {
      const uint32_t trackMask = isTrackSelected(track, ev.mask);
      if (trackMask == 0) {
        continue;
      }
      float p = track.p();
      float pt = track.pt();
      if (p < KFloatEpsilon) {
        continue;
      }
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }
      float eta = track.eta(), phi = track.phi();
      float sign = track.sign();
      if (TESTBIT(trackMask, kSystBase)) {
        histos.fill(HIST("hP"), p);
        histos.fill(HIST("hPt"), pt);
        histos.fill(HIST("hEta"), eta);
        histos.fill(HIST("hPhi"), phi);
      }

      float eff = 1.f, fake = 0.f, flatWeight = 1.f;
      const uint32_t weightMask = getTrackWeights(trackMask, ev.cent[kSystBase], vz, sign, pt, eta, phi, true, w, eff, fake, flatWeight);
      if (weightMask == 0) {
        continue;
      }
      if (TESTBIT(weightMask, kSystBase)) {
        fillFlat(kSystBase, vz, sign, pt, eta, phi, w[kSystBase], (1.f - fake) / eff);
      }

      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
          continue;
        }
        for (int v = 0; v < kNSyst; ++v) {
          if (TESTBIT(weightMask, v)) {
            sums[v][0].add(ieta, pt, w[v]);
          }
        }
      }
    }
    flushTrackCounts();

    // --- forward detectors: gain-equalised amplitudes ---
    FwdSignal fwd;
    if (cfgFwdOpt.cfgFwd) {
      loadFwdCalib(coll.template bc_as<aod::BCsWithTimestamps>(), true);
      computeFwdReco(coll, vz, fwd, TESTBIT(ev.mask, kSystBase));
    }

    // --- Poisson bootstrap: one weight per sample per event ---
    std::array<double, KMaxBoot> poisW{};
    if (doBoot) {
      for (int s = 0; s < nBoot; ++s) {
        poisW[s] = rng.Poisson(1.0);
      }
    }

    // --- c2, c3, subevent correlators and forward covariances of every variation ---
    EventObs evObs;
    for (int v = 0; v < kNSyst; ++v) {
      if (!TESTBIT(ev.mask, v)) {
        continue;
      }
      computeEventObs(sums[v][0], v, 0, cfgFwdOpt.cfgFwd ? &fwd : nullptr, nCh, evObs);
      fillFluc(hists.fluc[v][0], evObs, 1., nCh);
      // bootstrap replicas of Base, in the 5% centrality classes
      if (v == kSystBase && doBoot) {
        for (int s = 0; s < nBoot; ++s) {
          if (poisW[s] > 0.) {
            fillFluc(hists.boot[s], evObs, poisW[s], ev.cent[kSystBase]);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataFluc, "process data to calculate fluc pT", cfgRunDataFluc);

  // ===========================================================================
  // Cut QA
  // ===========================================================================
  // On the unfiltered collisions and tracks (the analysis passes see only the filtered ones):
  // every quantity the Base selection cuts on, before and after it. The counters are not
  // filled here. N_ch is the number of tracks of the analysis track filter.
  void processCutQA(GeneralCollisions::iterator const& col, UnfilteredTracks const& tracks)
  {
    int nCh = 0;
    for (const auto& trk : tracks) {
      if (trk.isGlobalTrack() && trk.pt() > KPtMin && trk.pt() < KPtMax) {
        ++nCh;
      }
    }
    const float cent = getCentrality(col, cfgEvSel.cfgCentralityChoice);
    fillCutQAEvent(kQABefore, col, nCh, cent);
    const auto ev = isEventSelected(col, nCh, 1 << kSystBase, false);
    if (!TESTBIT(ev.mask, kSystBase)) {
      return;
    }
    fillCutQAEvent(kQAAfter, col, nCh, cent);
    for (const auto& trk : tracks) {
      fillCutQATrack(kQABefore, trk);
      if (TESTBIT(isTrackSelected(trk, 1 << kSystBase, false), kSystBase)) {
        fillCutQATrack(kQAAfter, trk);
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processCutQA, "QA of the event and track cut quantities before and after the Base selection", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<RadialFlowDecorr>(cfgc)};
  return workflow;
}
