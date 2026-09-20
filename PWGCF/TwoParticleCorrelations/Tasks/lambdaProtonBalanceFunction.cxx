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

/// \file lambdaProtonBalanceFunction.cxx
/// \brief Single- and two-particle density correlations (rho1, rho2) and R2 correlations for Lambda-proton, proton-proton and Lambda-Lambda pairs (plus antiparticle combinations) to measure the balance function.
/// \author Anoop Poruthiyil <anoop.poruthiyil@cern.ch>

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/EndOfStreamContext.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA,
                           aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPr,
                           aod::pidTOFbeta, aod::TrackSelection>;

// Per-centrality-bin histogram set. Defined at file scope (NOT inside the task struct),
// otherwise O2's struct reflection in adaptAnalysisTask fails to compile.
struct CentRho2Set {
  // eta space (only booked if cFillEtaSpace && cFillCentEtaSpace)
  std::shared_ptr<TH2> Lp, LAp, ALp, ALAp, pp, pAp, App, ApAp, LL, LAL, ALL, ALAL;
  // y space
  std::shared_ptr<TH2> Lp_y, LAp_y, ALp_y, ALAp_y, pp_y, pAp_y, App_y, ApAp_y, LL_y, LAL_y, ALL_y, ALAL_y;
  // rho1 eta
  std::shared_ptr<TH2> rho1_Proton, rho1_AntiProton, rho1_Lambda, rho1_AntiLambda;
  // rho1 y
  std::shared_ptr<TH2> rho1_Proton_y, rho1_AntiProton_y, rho1_Lambda_y, rho1_AntiLambda_y;
  // event counter
  std::shared_ptr<TH1> hEvents;
};

struct LambdaProtonBalanceFunction {
  HistogramRegistry registryLambda{"Lambda_invMass", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryLambdaExtended{"Lambda_invMassExtended", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryRho{"rho1ANDrho2_LP", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryCent{"rho1ANDrho2_LP_CentralityBased", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryOther{"Other_Hists", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryPid{"PID", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQaDetector{"QA_Detec", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryProtonPidCheck{"ProtonCounts_byPID_Check", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryCorrelationQA{"CorrelationQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // ── Dedicated registry for V0-daughter veto importance QA ─────────────
  // Folder: ProtonVetoQA/
  // Contains per-event distributions of N_before, N_after, N_removed, fraction.
  HistogramRegistry registryProtonVetoQA{"ProtonVetoQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // pT bins 0.6 → 3.6 in steps of 0.2 (used for per-bin rho2 / mass histograms)
  std::vector<float> ptEdges;

  // ════════════════════════════════════════════════════════════════════════
  // SECTION A — Invariant mass histograms
  // ════════════════════════════════════════════════════════════════════════
  std::vector<std::shared_ptr<TH1>> hMassLambdaPtBins;
  std::vector<std::shared_ptr<TH1>> hMassAntiLambdaPtBins;
  std::shared_ptr<TH1> hMassLambdaMerged;
  std::shared_ptr<TH1> hMassAntiLambdaMerged;

  // ════════════════════════════════════════════════════════════════════════
  // SECTION A2 — Extended invariant-mass histograms (y × pT grid)
  //   hMassLambdaExtended[iY][iPt]     → Lambda_invMassExtended/Lambda/
  //   hMassAntiLambdaExtended[iY][iPt] → Lambda_invMassExtended/AntiLambda/
  //   Rapidity bins (kExtNyBins = 2): [-0.8,-0.7), [-0.7,-0.6]
  //   pT bins (kExtNptBins = 25):     [0.0,0.1), [0.1,0.2), ..., [2.4,2.5]
  // ════════════════════════════════════════════════════════════════════════
  static constexpr int kExtNyBins = 2;
  static constexpr int kExtNptBins = 25;
  static constexpr float kExtYMin = -0.8f;
  static constexpr float kExtYMax = -0.6f;
  static constexpr float kExtYStep = 0.1f;
  static constexpr float kExtPtMin = 0.0f;
  static constexpr float kExtPtMax = 2.5f;
  static constexpr float kExtPtStep = 0.1f;

  std::vector<std::vector<std::shared_ptr<TH1>>> hMassLambdaExtended;
  std::vector<std::vector<std::shared_ptr<TH1>>> hMassAntiLambdaExtended;

  // ════════════════════════════════════════════════════════════════════════
  // SECTION B — rho2 eta-space histograms
  //   Naming: hRho2_<channel>         = per-pT-bin vector
  //           hRho2_<channel>_pT015toMaxDefined = merged (full pT range as
  //           set by the pT cuts: protons 0.5-3.6, Lambdas 0.6-3.6 GeV/c;
  //           "015" is a legacy name and does NOT reflect the actual lower cut)
  // ════════════════════════════════════════════════════════════════════════
  // Lambda–proton family
  std::shared_ptr<TH2> hRho2_Lp_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_LAp_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ALp_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ALAp_pT015toMaxDefined;
  // Proton–proton family
  std::shared_ptr<TH2> hRho2_pp_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_pAp_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_App_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ApAp_pT015toMaxDefined;
  // Lambda–Lambda family
  std::shared_ptr<TH2> hRho2_LL_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_LAL_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ALL_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ALAL_pT015toMaxDefined;

  // ════════════════════════════════════════════════════════════════════════
  // SECTION C — rho2 rapidity-space histograms
  //   Same channels as Section B but filled with unrolledIndexY(y, phi)
  //   Naming: hRho2_<channel>_y / hRho2_<channel>_y_pT015toMaxDefined
  // ════════════════════════════════════════════════════════════════════════
  // Lambda–proton family (y)
  std::shared_ptr<TH2> hRho2_Lp_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_LAp_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ALp_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ALAp_y_pT015toMaxDefined;
  // Proton–proton family (y)
  std::shared_ptr<TH2> hRho2_pp_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_pAp_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_App_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ApAp_y_pT015toMaxDefined;
  // Lambda–Lambda family (y)
  std::shared_ptr<TH2> hRho2_LL_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_LAL_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ALL_y_pT015toMaxDefined;
  std::shared_ptr<TH2> hRho2_ALAL_y_pT015toMaxDefined;

  // SECTION E — centrality-differential storage (struct CentRho2Set is defined above the task)
  std::vector<CentRho2Set> centSets;
  std::vector<float> centEdgesLocal;

  // ════════════════════════════════════════════════════════════════════════
  // SECTION D — pT-spectra and pair-count monitoring histograms
  // ════════════════════════════════════════════════════════════════════════
  std::vector<std::shared_ptr<TH1>> hPtPrimProton_Lambda;
  std::vector<std::shared_ptr<TH1>> hPtPrimAntiProton_Lambda;
  std::vector<std::shared_ptr<TH1>> hPtPrimProton_AntiLambda;
  std::vector<std::shared_ptr<TH1>> hPtPrimAntiProton_AntiLambda;

  std::shared_ptr<TH1> hPtSelectedPrimProton_pT015toMaxDefined;
  std::shared_ptr<TH1> hPtSelectedPrimAntiProton_pT015toMaxDefined;
  std::shared_ptr<TH1> hPtSelLambda_pT015toMaxDefined;
  std::shared_ptr<TH1> hPtSelAntiLambda_pT015toMaxDefined;

  // ── PID-check histograms (ProtonCounts_byPID_Check) ───────────────────
  std::shared_ptr<TH1> h_TPC;
  std::shared_ptr<TH1> h_TPCandTOF;
  std::shared_ptr<TH1> h_TPCorTOF;
  // ─────────────────────────────────────────────────────────────────────

  // ── Run-level accumulators for proton-veto importance QA ─────────────
  // These are summed over all processed events and printed in endOfStream().
  // They quantify how much Lambda feed-down was removed from the proton
  // sample before constructing future Balance Functions.
  int64_t runTotalProtonBeforeVeto = 0;  ///< sum of N_before_veto over all events
  int64_t runTotalProtonRemovedVeto = 0; ///< sum of N_removed_by_veto over all events
  int64_t runTotalAntiPBeforeVeto = 0;   ///< antiproton: sum of N_before_veto
  int64_t runTotalAntiPRemovedVeto = 0;  ///< antiproton: sum of N_removed_by_veto
  // ─────────────────────────────────────────────────────────────────────

  // ── Configurables ─────────────────────────────────────────────────────

  Configurable<float> cPPandV0ZVertexCut{"cPPandV0ZVertexCut", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> cFillEtaSpace{"cFillEtaSpace", false, "Fill eta-space rho1/rho2 histograms"};
  Configurable<bool> cUseQinvCut{"cUseQinvCut", true, "Apply q_inv > 0.01 GeV/c pair cut (false = no q_inv rejection)"};
  Configurable<bool> cFillCentHists{"cFillCentHists", true, "Fill centrality-differential rho1/rho2 histograms in separate rho1ANDrho2_LP_CentralityBased folder"};

  // Track quality cuts
  Configurable<int> pProtonTPCMinRows{"pProtonTPCMinRows", 70, "Minimum TPC crossed rows"};
  Configurable<float> pProtonTPCMinRowsOverFindable{"pProtonTPCMinRowsOverFindable", 0.8f, "Min TPC crossed rows / findable clusters"};
  Configurable<float> pProtonTPCMaxChi2PerCluster{"pProtonTPCMaxChi2PerCluster", 4.0f, "Max TPC chi2 / Ncls"};
  Configurable<int> pProtonITSMinClusters{"pProtonITSMinClusters", 5, "Minimum ITS clusters"};
  Configurable<float> pProtonITSMaxChi2PerCluster{"pProtonITSMaxChi2PerCluster", 36.0f, "Max ITS chi2 / Ncls"};
  Configurable<float> cPPMinTpcSignal{"cPPMinTpcSignal", 10.0f, "Min raw TPC signal (sanity cut)"};
  // DCA cuts
  Configurable<float> pProtonMaxDCAxy{"pProtonMaxDCAxy", 0.02f, "Max |DCAxy| (cm)"}; // tight cut for cleanest primary-proton sample
  Configurable<float> pProtonMaxDCAz{"pProtonMaxDCAz", 0.1f, "Max |DCAz| (cm)"};

  Configurable<float> pProtonTPCTOFSwitchP{"pProtonTPCTOFSwitchP", 0.700000238f, "p threshold (GeV/c) above which TOF is used if available"};
  // Configurable<float> cPPMinTofBeta{"cPPMinTofBeta", 0.4f, "Min TOF beta for proton selection (TOF branch only)"}; //BETACUTCommented
  // Configurable<float> cPPMaxTofBeta{"cPPMaxTofBeta", 0.96f, "Max TOF beta for proton selection (TOF branch only)"}; //BETACUTCommented

  // V0 topology cuts + V0 DCA Cuts
  Configurable<float> cV0MaxDcaDaughters{"cV0MaxDcaDaughters", 1.0f, "DCA V0 daughters < (cm)"};
  Configurable<float> lambdaV0DaughterProtonMinDCAToPV{"lambdaV0DaughterProtonMinDCAToPV", 0.01f, "DCA proton daughter to PV > (cm)"};
  Configurable<float> lambdaV0DaughterPionMinDCAToPV{"lambdaV0DaughterPionMinDCAToPV", 0.10f, "DCA pion daughter to PV > (cm)"};
  Configurable<double> lambdaV0MinCosPA{"lambdaV0MinCosPA", 0.995, "V0 CosPA >"};
  Configurable<float> lambdaV0MinDecayRadius{"lambdaV0MinDecayRadius", 0.5f, "V0 radius > (cm)"};
  Configurable<float> cV0MaxDcaToPV{"cV0MaxDcaToPV", 999.0f, "DCA V0 to PV < (cm)"};

  // PID cuts
  Configurable<float> lambdaV0DaughterPionTPCNsigma{"lambdaV0DaughterPionTPCNsigma", 2.0f, "lambdaV0DaughterPionTPCNsigma"};
  Configurable<float> pProtonTPCNsigma{"pProtonTPCNsigma", 2.0f, "pProtonTPCNsigma"};
  Configurable<float> lambdaV0DaughterProtonTPCNsigma{"lambdaV0DaughterProtonTPCNsigma", 2.0f, "lambdaV0DaughterProtonTPCNsigma"};
  Configurable<float> pProtonTOFNsigma{"pProtonTOFNsigma", 2.0f, "pProtonTOFNsigma"};

  // // V0 daughter beta cuts (particle-dependent) //BETACUTCommented
  // Configurable<float> cV0protonDauMinBeta{"cV0protonDauMinBeta", 0.4f, "Min TOF beta for (anti)proton V0 daughter"}; //BETACUTCommented
  // Configurable<float> cV0protonDauMaxBeta{"cV0protonDauMaxBeta", 0.98f, "Max TOF beta for (anti)proton V0 daughter"}; //BETACUTCommented
  // Configurable<float> cV0pionDauMinBeta{"cV0pionDauMinBeta", 0.91f, "Min TOF beta for (anti)pion V0 daughter"}; //BETACUTCommented
  // ─────────────────────────────────────────────────────────────────────

  // Proton Kinematic Cuts (from strange9)
  Configurable<float> pProtonMinP{"pProtonMinP", 0.5f, "Min p for primary protons (GeV/c)"};
  Configurable<float> pProtonMaxP{"pProtonMaxP", 3.6f, "Max p for primary protons (GeV/c)"};
  Configurable<float> pProtonMinPt{"pProtonMinPt", 0.5f, "Min pT for primary protons (GeV/c)"};
  Configurable<float> pProtonMaxPt{"pProtonMaxPt", 3.6f, "Max pT for primary protons (GeV/c)"};

  Configurable<float> pProtonMaxY{"pProtonMaxY", 0.5f, "Max |y| for primary protons"};

  // Lambda Kinematic Cuts (from strange9)
  Configurable<float> lambdaV0MinPt{"lambdaV0MinPt", 0.60f, "Min pT for Lambdas"};
  Configurable<float> lambdaV0MaxPt{"lambdaV0MaxPt", 3.60f, "Max pT for Lambdas"};
  Configurable<float> lambdaV0MaxY{"lambdaV0MaxY", 0.5f, "Max |y| for Lambdas"};

  // V0 Mass and Topology Cuts (from strange9)
  Configurable<float> lambdaV0MaxCTau{"lambdaV0MaxCTau", 30.0f, "Max V0 ctau (cm)"};
  Configurable<float> lambdaV0KShortRejectMassWindow{"lambdaV0KShortRejectMassWindow", 0.01f, "K0s rejection window"};
  Configurable<float> lambdaV0MassWindow{"lambdaV0MassWindow", 0.007f, "Lambda mass window (|m - mPDG| < this)"};

  // V0 Type and Daughter Quality Cuts (from strange9)
  Configurable<int> cV0TypeSelection{"cV0TypeSelection", 1, "V0 Type Selection"};
  Configurable<float> cV0DauMinPt{"cV0DauMinPt", 0.1f, "Daughter pT minimum"};
  Configurable<float> cV0DauMaxEta{"cV0DauMaxEta", 0.8f, "Daughter |eta| cut"};
  Configurable<int> cV0DauMinTpcCrossedRows{"cV0DauMinTpcCrossedRows", 70, "Daughter TPC min crossed rows"};
  // Rapidity-specific constants matched to pProtonMaxY's default value (0.5).
  // Note: Since these are compile-time constexpr, the axis will NOT auto-update
  // if pProtonMaxY is changed via configurable at runtime. The axis just needs to
  // stay >= the cut, not exactly equal.

  // Eta
  static constexpr int kRhoEtaBins = 40;
  static constexpr float kRhoMin = -0.5f;
  static constexpr float kRhoMax = 0.5f;

  // Rapidity
  static constexpr int kRhoYBins = 10;
  static constexpr float kRhoYMin = -0.5f;
  static constexpr float kRhoYMax = 0.5f;

  static constexpr int kRhoPhiBins = 72;
  static constexpr int kRhoUnrolledBins = kRhoEtaBins * kRhoPhiBins;
  static constexpr int kRhoUnrolledBinsY = kRhoYBins * kRhoPhiBins;

  void buildPtBins()
  {
    ptEdges.clear();
    static constexpr int kNPtEdges = 16; // 0.6, 0.8, ..., 3.6 GeV/c (15 bins)
    static constexpr float kPtEdgeMin = 0.6f;
    static constexpr float kPtEdgeStep = 0.2f;
    for (int i = 0; i < kNPtEdges; ++i) {
      ptEdges.push_back(kPtEdgeMin + kPtEdgeStep * static_cast<float>(i));
    }
  }

  int etaBinIndex(float eta) const
  {
    if (eta < kRhoMin || eta >= kRhoMax) {
      return -1;
    }
    const float binWidth = (kRhoMax - kRhoMin) / kRhoEtaBins;
    return static_cast<int>((eta - kRhoMin) / binWidth);
  }

  int phiBinIndex(float phi) const
  {
    constexpr float phiMin = 0.0f;
    constexpr float phiMax = o2::constants::math::TwoPI;
    if (phi < phiMin || phi >= phiMax) {
      return -1;
    }
    const float binWidth = (phiMax - phiMin) / kRhoPhiBins;
    return static_cast<int>((phi - phiMin) / binWidth);
  }

  int unrolledIndex(float eta, float phi) const
  {
    const int iEta = etaBinIndex(eta);
    const int iPhi = phiBinIndex(phi);
    if (iEta < 0 || iPhi < 0) {
      return -1;
    }
    return iEta * kRhoPhiBins + iPhi;
  }

  int yBinIndex(float y) const
  {
    if (y < kRhoYMin || y >= kRhoYMax) {
      return -1;
    }
    const float binWidth = (kRhoYMax - kRhoYMin) / kRhoYBins;
    return static_cast<int>((y - kRhoYMin) / binWidth);
  }

  int unrolledIndexY(float y, float phi) const
  {
    const int iY = yBinIndex(y);
    const int iPhi = phiBinIndex(phi);
    if (iY < 0 || iPhi < 0) {
      return -1;
    }
    return iY * kRhoPhiBins + iPhi;
  }

  // Convenience wrappers with particle-specific masses
  static constexpr float kMassProton = o2::constants::physics::MassProton;
  static constexpr float kMassLambda = o2::constants::physics::MassLambda;

  // Femtoscopic q_inv cut: pairs with q_inv <= kQinvCutLP are rejected from all rho2 histograms (pp, pAp, App, ApAp, LP-family, LL-family):
  static constexpr float kQinvCutLP = 0.01f;

  template <typename T>
  float protonRapidity(T const& track) const
  {
    return RecoDecay::y(std::array<float, 3>{track.px(), track.py(), track.pz()}, kMassProton);
  }
  // Fills both eta-space and rapidity-space rho2 histograms in one call.
  // Separated guards prevent filling when either index is out of acceptance.
  void fillRho2Pair(
    std::shared_ptr<TH2>& histEta, std::shared_ptr<TH2>& histY,
    int idxEta1, int idxEta2, int idxY1, int idxY2) const
  {
    if (idxEta1 >= 0 && idxEta2 >= 0) {
      histEta->Fill(idxEta1, idxEta2);
    }
    if (idxY1 >= 0 && idxY2 >= 0) {
      histY->Fill(idxY1, idxY2);
    }
  }

  float displayDeltaPhi(float dPhi) const
  {
    return RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI / 2.0f);
  }

  // k* = |p*| of one particle in the pair rest frame, from four-momenta.
  float computeKstar(float px1, float py1, float pz1, float m1,
                     float px2, float py2, float pz2, float m2) const
  {
    const float E1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + m1 * m1);
    const float E2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + m2 * m2);
    const float sx = px1 + px2;
    const float sy = py1 + py2;
    const float sz = pz1 + pz2;
    const float sE = E1 + E2;
    const float s = sE * sE - sx * sx - sy * sy - sz * sz;
    if (s <= 0.f) {
      return 0.f;
    }
    // k*^2 = [ (s - (m1+m2)^2)(s - (m1-m2)^2) ] / (4 s)
    const float a = s - (m1 + m2) * (m1 + m2);
    const float b = s - (m1 - m2) * (m1 - m2);
    return std::sqrt(std::max(a * b, 0.f) / (4.f * s));
  }

  // ── Femtoscopic q_inv computation ────────────────────────────────────────
  // Returns q_inv = sqrt( -(p1 - p2)^2 ) in natural units (GeV/c).
  // E1 and E2 must be computed by the caller using the appropriate mass
  // hypothesis for each track:
  //   - primary proton:  E = sqrt(p^2 + kMassProton^2)
  //   - V0 (Lambda leg): E = sqrt(p^2 + kMassLambda^2)
  float computeQinv(float px1, float py1, float pz1, float E1,
                    float px2, float py2, float pz2, float E2) const
  {
    const float dE = E1 - E2;
    const float dpx = px1 - px2;
    const float dpy = py1 - py2;
    const float dpz = pz1 - pz2;
    const float q2 = dE * dE - dpx * dpx - dpy * dpy - dpz * dpz;
    // q2 should be <= 0 for physical pairs; protect sqrt from numerical noise
    return std::sqrt(std::max(-q2, 0.0f));
  }

  // true => pair is rejected by the q_inv cut (never true when the cut is switched off)
  bool failsQinvCut(float qinv) const
  {
    return cUseQinvCut.value && (qinv <= kQinvCutLP);
  }

  // Returns index of centrality bin for cent, or -1 if outside all bins.
  int centBinIndex(float cent) const
  {
    for (size_t i = 0; i + 1 < centEdgesLocal.size(); ++i) {
      if (cent >= centEdgesLocal[i] && cent < centEdgesLocal[i + 1]) {
        return static_cast<int>(i);
      }
    }
    return -1;
  }

  // Fill helper: fills h (if valid) at (a,b) only when both indices are in range.
  static void fillIdx(std::shared_ptr<TH2> const& h, int a, int b)
  {
    if (h && a >= 0 && b >= 0) {
      h->Fill(a, b);
    }
  }
  // ─────────────────────────────────────────────────────────────────────

  // ── Two-region PID selection for (anti)protons ────────────────────────
  template <typename TTrack>
  bool passesPrimProtonPid(TTrack const& trk) const
  {
    const float pTPC = trk.tpcInnerParam(); // NOTE: pTPC, not trk.p() — matches reference exactly
    const float nsPr = trk.tpcNSigmaPr();

    if (std::abs(nsPr) >= pProtonTPCNsigma.value) {
      return false;
    }
    if (pTPC < pProtonTPCTOFSwitchP.value) {
      return true;
    }
    if (!trk.hasTOF()) {
      return false;
    }
    return std::abs(trk.tofNSigmaPr()) < pProtonTOFNsigma.value;
  }
  // ─────────────────────────────────────────────────────────────────────

  // ══════════════════════════════════════════════════════════════════════════
  // STAGED QA INFRASTRUCTURE — registryQaDetector ("QA_Detec")
  // ══════════════════════════════════════════════════════════════════════════

  /// Canonical analysis stages for (anti)proton QA fills.
  enum class ProtonQAStage : int {
    RawAfterTrackSel = 0, ///< after track-quality selection (TrackSelection table), before DCA cuts and PID
    TpcPid = 1,           ///< after DCA + TPC PID cut passes (all momenta)
    TpcTofPid = 2,        ///< after DCA + full PID: TPC everywhere, plus TOF required only above pProtonTPCTOFSwitchP
    FinalSelected = 3     ///< entered selectedPrimProtons / selectedPrimAntiProtons
  };

  // ── Staged fill helper: Kinematics ───────────────────────────────────────
  // IsProton=true → Proton subfolder ; IsProton=false → AntiProton subfolder
  template <bool IsProton, typename TTrack>
  void fillProtonKinematics(ProtonQAStage stage, TTrack const& trk)
  {
    const float p = trk.p();
    const float pt = trk.pt();
    const float eta = trk.eta();
    const float phi = trk.phi();
    const float y = protonRapidity(trk);
    if (stage == ProtonQAStage::RawAfterTrackSel) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/01_Kinematics/01_RAW_AfterTrackSel/h1f_p"), p);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/01_RAW_AfterTrackSel/h1f_pt"), pt);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/01_RAW_AfterTrackSel/h1f_eta"), eta);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/01_RAW_AfterTrackSel/h1f_phi"), phi);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/01_RAW_AfterTrackSel/h1f_rapidity"), y);
      } else {
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/01_RAW_AfterTrackSel/h1f_p"), p);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/01_RAW_AfterTrackSel/h1f_pt"), pt);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/01_RAW_AfterTrackSel/h1f_eta"), eta);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/01_RAW_AfterTrackSel/h1f_phi"), phi);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/01_RAW_AfterTrackSel/h1f_rapidity"), y);
      }
    } else if (stage == ProtonQAStage::TpcPid) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/01_Kinematics/02_TPC_PID/h1f_p"), p);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/02_TPC_PID/h1f_pt"), pt);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/02_TPC_PID/h1f_eta"), eta);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/02_TPC_PID/h1f_phi"), phi);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/02_TPC_PID/h1f_rapidity"), y);
      } else {
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/02_TPC_PID/h1f_p"), p);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/02_TPC_PID/h1f_pt"), pt);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/02_TPC_PID/h1f_eta"), eta);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/02_TPC_PID/h1f_phi"), phi);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/02_TPC_PID/h1f_rapidity"), y);
      }
    } else if (stage == ProtonQAStage::TpcTofPid) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/01_Kinematics/03_TPC_TOF_PID/h1f_p"), p);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/03_TPC_TOF_PID/h1f_pt"), pt);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/03_TPC_TOF_PID/h1f_eta"), eta);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/03_TPC_TOF_PID/h1f_phi"), phi);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/03_TPC_TOF_PID/h1f_rapidity"), y);
      } else {
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/03_TPC_TOF_PID/h1f_p"), p);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/03_TPC_TOF_PID/h1f_pt"), pt);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/03_TPC_TOF_PID/h1f_eta"), eta);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/03_TPC_TOF_PID/h1f_phi"), phi);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/03_TPC_TOF_PID/h1f_rapidity"), y);
      }
    } else { // FinalSelected
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/01_Kinematics/04_Final_SelectedProton/h1f_p"), p);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/04_Final_SelectedProton/h1f_pt"), pt);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/04_Final_SelectedProton/h1f_eta"), eta);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/04_Final_SelectedProton/h1f_phi"), phi);
        registryQaDetector.fill(HIST("Proton/01_Kinematics/04_Final_SelectedProton/h1f_rapidity"), y);
      } else {
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/04_Final_SelectedProton/h1f_p"), p);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/04_Final_SelectedProton/h1f_pt"), pt);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/04_Final_SelectedProton/h1f_eta"), eta);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/04_Final_SelectedProton/h1f_phi"), phi);
        registryQaDetector.fill(HIST("AntiProton/01_Kinematics/04_Final_SelectedProton/h1f_rapidity"), y);
      }
    }
  }

  // ── Staged fill helper: DetectorQuality ──────────────────────────────────
  template <bool IsProton, typename TTrack>
  void fillProtonDetectorQuality(ProtonQAStage stage, TTrack const& trk)
  {
    const float pt = trk.pt();
    const auto tpcCR = static_cast<float>(trk.tpcNClsCrossedRows());
    const float tpcCRoF = trk.tpcCrossedRowsOverFindableCls();
    const float tpcChi2 = trk.tpcChi2NCl();
    const auto itsN = static_cast<float>(trk.itsNCls());
    const float itsChi2 = trk.itsChi2NCl();
    const float dcaXY = trk.dcaXY();
    const float dcaZ = trk.dcaZ();
    if (stage == ProtonQAStage::RawAfterTrackSel) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_tpcCrossedRows"), tpcCR);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_tpcCrossedRowsOverFindable"), tpcCRoF);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_tpcChi2NCl"), tpcChi2);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_itsNCls"), itsN);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_itsChi2NCl"), itsChi2);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/01_RAW_AfterTrackSel/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/01_RAW_AfterTrackSel/h2f_dcaZ_vs_pt"), pt, dcaZ);
      } else {
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_tpcCrossedRows"), tpcCR);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_tpcCrossedRowsOverFindable"), tpcCRoF);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_tpcChi2NCl"), tpcChi2);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_itsNCls"), itsN);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/01_RAW_AfterTrackSel/h1f_itsChi2NCl"), itsChi2);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/01_RAW_AfterTrackSel/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/01_RAW_AfterTrackSel/h2f_dcaZ_vs_pt"), pt, dcaZ);
      }
    } else if (stage == ProtonQAStage::TpcPid) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/02_TPC_PID/h1f_tpcCrossedRows"), tpcCR);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/02_TPC_PID/h1f_tpcCrossedRowsOverFindable"), tpcCRoF);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/02_TPC_PID/h1f_tpcChi2NCl"), tpcChi2);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/02_TPC_PID/h1f_itsNCls"), itsN);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/02_TPC_PID/h1f_itsChi2NCl"), itsChi2);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/02_TPC_PID/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/02_TPC_PID/h2f_dcaZ_vs_pt"), pt, dcaZ);
      } else {
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/02_TPC_PID/h1f_tpcCrossedRows"), tpcCR);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/02_TPC_PID/h1f_tpcCrossedRowsOverFindable"), tpcCRoF);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/02_TPC_PID/h1f_tpcChi2NCl"), tpcChi2);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/02_TPC_PID/h1f_itsNCls"), itsN);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/02_TPC_PID/h1f_itsChi2NCl"), itsChi2);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/02_TPC_PID/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/02_TPC_PID/h2f_dcaZ_vs_pt"), pt, dcaZ);
      }
    } else if (stage == ProtonQAStage::TpcTofPid) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/03_TPC_TOF_PID/h1f_tpcCrossedRows"), tpcCR);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/03_TPC_TOF_PID/h1f_tpcCrossedRowsOverFindable"), tpcCRoF);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/03_TPC_TOF_PID/h1f_tpcChi2NCl"), tpcChi2);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/03_TPC_TOF_PID/h1f_itsNCls"), itsN);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/03_TPC_TOF_PID/h1f_itsChi2NCl"), itsChi2);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/03_TPC_TOF_PID/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/03_TPC_TOF_PID/h2f_dcaZ_vs_pt"), pt, dcaZ);
      } else {
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/03_TPC_TOF_PID/h1f_tpcCrossedRows"), tpcCR);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/03_TPC_TOF_PID/h1f_tpcCrossedRowsOverFindable"), tpcCRoF);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/03_TPC_TOF_PID/h1f_tpcChi2NCl"), tpcChi2);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/03_TPC_TOF_PID/h1f_itsNCls"), itsN);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/03_TPC_TOF_PID/h1f_itsChi2NCl"), itsChi2);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/03_TPC_TOF_PID/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/03_TPC_TOF_PID/h2f_dcaZ_vs_pt"), pt, dcaZ);
      }
    } else { // FinalSelected
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/04_Final_SelectedProton/h1f_tpcCrossedRows"), tpcCR);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/04_Final_SelectedProton/h1f_tpcCrossedRowsOverFindable"), tpcCRoF);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/04_Final_SelectedProton/h1f_tpcChi2NCl"), tpcChi2);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/04_Final_SelectedProton/h1f_itsNCls"), itsN);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/04_Final_SelectedProton/h1f_itsChi2NCl"), itsChi2);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/04_Final_SelectedProton/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("Proton/02_DetectorQuality/04_Final_SelectedProton/h2f_dcaZ_vs_pt"), pt, dcaZ);
      } else {
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/04_Final_SelectedProton/h1f_tpcCrossedRows"), tpcCR);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/04_Final_SelectedProton/h1f_tpcCrossedRowsOverFindable"), tpcCRoF);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/04_Final_SelectedProton/h1f_tpcChi2NCl"), tpcChi2);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/04_Final_SelectedProton/h1f_itsNCls"), itsN);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/04_Final_SelectedProton/h1f_itsChi2NCl"), itsChi2);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/04_Final_SelectedProton/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("AntiProton/02_DetectorQuality/04_Final_SelectedProton/h2f_dcaZ_vs_pt"), pt, dcaZ);
      }
    }
  }

  // ── Staged fill helper: TPC_QA ───────────────────────────────────────────
  template <bool IsProton, typename TTrack>
  void fillProtonTpcQA(ProtonQAStage stage, TTrack const& trk)
  {
    const float p = trk.p();                // total momentum (kept for vs_p plots)
    const float pTPC = trk.tpcInnerParam(); // TPC inner param = proton PID selection variable
    const float pt = trk.pt();
    const float dEdx = trk.tpcSignal();
    const float nsPr = trk.tpcNSigmaPr();
    const float nsPi = trk.tpcNSigmaPi();
    if (stage == ProtonQAStage::RawAfterTrackSel) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_dEdx_vs_p"), p, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_dEdx_vs_pTPC"), pTPC, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_dEdx_vs_pt"), pt, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_p"), p, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pTPC"), pTPC, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pt"), pt, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_nsigmaPi_vs_pt"), pt, nsPi);
      } else {
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_dEdx_vs_p"), p, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_dEdx_vs_pTPC"), pTPC, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_dEdx_vs_pt"), pt, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_p"), p, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pTPC"), pTPC, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pt"), pt, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/01_RAW_AfterTrackSel/h2f_nsigmaPi_vs_pt"), pt, nsPi);
      }
    } else if (stage == ProtonQAStage::TpcPid) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/02_TPC_PID/h2f_dEdx_vs_p"), p, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/02_TPC_PID/h2f_dEdx_vs_pTPC"), pTPC, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/02_TPC_PID/h2f_dEdx_vs_pt"), pt, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/02_TPC_PID/h2f_nsigmaPr_vs_p"), p, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/02_TPC_PID/h2f_nsigmaPr_vs_pTPC"), pTPC, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/02_TPC_PID/h2f_nsigmaPr_vs_pt"), pt, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/02_TPC_PID/h2f_nsigmaPi_vs_pt"), pt, nsPi);
      } else {
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/02_TPC_PID/h2f_dEdx_vs_p"), p, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/02_TPC_PID/h2f_dEdx_vs_pTPC"), pTPC, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/02_TPC_PID/h2f_dEdx_vs_pt"), pt, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/02_TPC_PID/h2f_nsigmaPr_vs_p"), p, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/02_TPC_PID/h2f_nsigmaPr_vs_pTPC"), pTPC, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/02_TPC_PID/h2f_nsigmaPr_vs_pt"), pt, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/02_TPC_PID/h2f_nsigmaPi_vs_pt"), pt, nsPi);
      }
    } else if (stage == ProtonQAStage::TpcTofPid) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/03_TPC_TOF_PID/h2f_dEdx_vs_p"), p, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/03_TPC_TOF_PID/h2f_dEdx_vs_pTPC"), pTPC, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/03_TPC_TOF_PID/h2f_dEdx_vs_pt"), pt, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_p"), p, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pTPC"), pTPC, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pt"), pt, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/03_TPC_TOF_PID/h2f_nsigmaPi_vs_pt"), pt, nsPi);
      } else {
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/03_TPC_TOF_PID/h2f_dEdx_vs_p"), p, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/03_TPC_TOF_PID/h2f_dEdx_vs_pTPC"), pTPC, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/03_TPC_TOF_PID/h2f_dEdx_vs_pt"), pt, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_p"), p, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pTPC"), pTPC, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pt"), pt, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/03_TPC_TOF_PID/h2f_nsigmaPi_vs_pt"), pt, nsPi);
      }
    } else { // FinalSelected
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/04_Final_SelectedProton/h2f_dEdx_vs_p"), p, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/04_Final_SelectedProton/h2f_dEdx_vs_pTPC"), pTPC, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/04_Final_SelectedProton/h2f_dEdx_vs_pt"), pt, dEdx);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_p"), p, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pTPC"), pTPC, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pt"), pt, nsPr);
        registryQaDetector.fill(HIST("Proton/03_TPC_QA/04_Final_SelectedProton/h2f_nsigmaPi_vs_pt"), pt, nsPi);
      } else {
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/04_Final_SelectedProton/h2f_dEdx_vs_p"), p, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/04_Final_SelectedProton/h2f_dEdx_vs_pTPC"), pTPC, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/04_Final_SelectedProton/h2f_dEdx_vs_pt"), pt, dEdx);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_p"), p, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pTPC"), pTPC, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pt"), pt, nsPr);
        registryQaDetector.fill(HIST("AntiProton/03_TPC_QA/04_Final_SelectedProton/h2f_nsigmaPi_vs_pt"), pt, nsPi);
      }
    }
  }

  // ── Staged fill helper: TOF_QA (call only when trk.hasTOF() == true) ─────
  template <bool IsProton, typename TTrack>
  void fillProtonTofQA(ProtonQAStage stage, TTrack const& trk)
  {
    const float p = trk.p();                // total momentum (kept for vs_p plots)
    const float pTPC = trk.tpcInnerParam(); // TPC inner param = proton PID selection variable
    const float pTOF = trk.tofExpMom();     // momentum at the TOF radius (TracksExtra::TOFExpMom)
    const float pt = trk.pt();
    const float beta = trk.beta();
    const float nsTOF = trk.tofNSigmaPr();
    if (stage == ProtonQAStage::RawAfterTrackSel) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_beta_vs_p"), p, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_beta_vs_pTPC"), pTPC, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_beta_vs_pt"), pt, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_beta_vs_pTOF"), pTOF, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_p"), p, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pTPC"), pTPC, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pt"), pt, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pTOF"), pTOF, nsTOF);
      } else {
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_beta_vs_p"), p, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_beta_vs_pTPC"), pTPC, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_beta_vs_pt"), pt, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_beta_vs_pTOF"), pTOF, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_p"), p, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pTPC"), pTPC, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pt"), pt, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/01_RAW_AfterTrackSel/h2f_nsigmaPr_vs_pTOF"), pTOF, nsTOF);
      }
    } else if (stage == ProtonQAStage::TpcPid) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/02_TPC_PID/h2f_beta_vs_p"), p, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/02_TPC_PID/h2f_beta_vs_pTPC"), pTPC, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/02_TPC_PID/h2f_beta_vs_pt"), pt, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/02_TPC_PID/h2f_beta_vs_pTOF"), pTOF, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/02_TPC_PID/h2f_nsigmaPr_vs_p"), p, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/02_TPC_PID/h2f_nsigmaPr_vs_pTPC"), pTPC, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/02_TPC_PID/h2f_nsigmaPr_vs_pt"), pt, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/02_TPC_PID/h2f_nsigmaPr_vs_pTOF"), pTOF, nsTOF);
      } else {
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/02_TPC_PID/h2f_beta_vs_p"), p, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/02_TPC_PID/h2f_beta_vs_pTPC"), pTPC, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/02_TPC_PID/h2f_beta_vs_pt"), pt, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/02_TPC_PID/h2f_beta_vs_pTOF"), pTOF, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/02_TPC_PID/h2f_nsigmaPr_vs_p"), p, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/02_TPC_PID/h2f_nsigmaPr_vs_pTPC"), pTPC, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/02_TPC_PID/h2f_nsigmaPr_vs_pt"), pt, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/02_TPC_PID/h2f_nsigmaPr_vs_pTOF"), pTOF, nsTOF);
      }
    } else if (stage == ProtonQAStage::TpcTofPid) {
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/03_TPC_TOF_PID/h2f_beta_vs_p"), p, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/03_TPC_TOF_PID/h2f_beta_vs_pTPC"), pTPC, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/03_TPC_TOF_PID/h2f_beta_vs_pt"), pt, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/03_TPC_TOF_PID/h2f_beta_vs_pTOF"), pTOF, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_p"), p, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pTPC"), pTPC, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pt"), pt, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pTOF"), pTOF, nsTOF);
      } else {
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/03_TPC_TOF_PID/h2f_beta_vs_p"), p, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/03_TPC_TOF_PID/h2f_beta_vs_pTPC"), pTPC, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/03_TPC_TOF_PID/h2f_beta_vs_pt"), pt, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/03_TPC_TOF_PID/h2f_beta_vs_pTOF"), pTOF, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_p"), p, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pTPC"), pTPC, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pt"), pt, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/03_TPC_TOF_PID/h2f_nsigmaPr_vs_pTOF"), pTOF, nsTOF);
      }
    } else { // FinalSelected
      if constexpr (IsProton) {
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/04_Final_SelectedProton/h2f_beta_vs_p"), p, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/04_Final_SelectedProton/h2f_beta_vs_pTPC"), pTPC, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/04_Final_SelectedProton/h2f_beta_vs_pt"), pt, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/04_Final_SelectedProton/h2f_beta_vs_pTOF"), pTOF, beta);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_p"), p, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pTPC"), pTPC, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pt"), pt, nsTOF);
        registryQaDetector.fill(HIST("Proton/04_TOF_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pTOF"), pTOF, nsTOF);
      } else {
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/04_Final_SelectedProton/h2f_beta_vs_p"), p, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/04_Final_SelectedProton/h2f_beta_vs_pTPC"), pTPC, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/04_Final_SelectedProton/h2f_beta_vs_pt"), pt, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/04_Final_SelectedProton/h2f_beta_vs_pTOF"), pTOF, beta);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_p"), p, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pTPC"), pTPC, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pt"), pt, nsTOF);
        registryQaDetector.fill(HIST("AntiProton/04_TOF_QA/04_Final_SelectedProton/h2f_nsigmaPr_vs_pTOF"), pTOF, nsTOF);
      }
    }
  }
  enum class LambdaQAStage { RawCandidate,
                             Final };
  enum class DaughterParent { Lambda,
                              AntiLambda };
  enum class DaughterType { Proton,
                            Pion,
                            AntiProton,
                            AntiPion };
  enum class DaughterQAStage { Raw,
                               Final };

  // ── Staged fill helper: Lambda QA ──────────────────────────────────────────
  template <bool IsLambda, typename TV0>
  void fillLambdaQA(LambdaQAStage stage, TV0 const& v0, float ctau, float bgctau)
  {
    const float pt = v0.pt();
    const float eta = v0.eta();
    const float phi = v0.phi();
    const float rap = IsLambda ? v0.rapidity(1) : v0.rapidity(2);
    const float cosPA = v0.v0cosPA();
    const float dcaToPV = v0.dcav0topv();
    const float dcaDau = v0.dcaV0daughters();

    // Need collision position for accurate decay coordinates, but we use the provided bgctau
    // and radius for topology, and we can pull alpha and qtarm.
    // VS Code: commented out unused variable.
    // constexpr float mLambdaPDG = 1.115683f;
    // const float ctau   = v0.distovertotmom(0, 0, 0) * mLambdaPDG; // simplified for this block, actual filled in process
    const float alpha = v0.alpha();
    const float qtarm = v0.qtarm();
    const float r = v0.v0radius();

    if constexpr (IsLambda) {
      if (stage == LambdaQAStage::RawCandidate) {
        registryQaDetector.fill(HIST("Lambda/01_Kinematics/01_RawCandidate/h1f_pt"), pt);
        registryQaDetector.fill(HIST("Lambda/01_Kinematics/01_RawCandidate/h1f_eta"), eta);
        registryQaDetector.fill(HIST("Lambda/01_Kinematics/01_RawCandidate/h1f_phi"), phi);

        registryQaDetector.fill(HIST("Lambda/02_Topology/01_RawCandidate/h1f_cospa"), cosPA);
        registryQaDetector.fill(HIST("Lambda/02_Topology/01_RawCandidate/h1f_cospa_zoom"), cosPA);
        registryQaDetector.fill(HIST("Lambda/02_Topology/01_RawCandidate/h1f_dcaV0toPV"), dcaToPV);
        registryQaDetector.fill(HIST("Lambda/02_Topology/01_RawCandidate/h2f_armpod"), alpha, qtarm);
        registryQaDetector.fill(HIST("Lambda/02_Topology/01_RawCandidate/h1f_radius"), r);
        registryQaDetector.fill(HIST("Lambda/02_Topology/01_RawCandidate/h1f_dcaDaughters"), dcaDau);
        registryQaDetector.fill(HIST("Lambda/02_Topology/01_RawCandidate/h1f_ctau"), ctau);
        registryQaDetector.fill(HIST("Lambda/02_Topology/01_RawCandidate/h1f_bgctau"), bgctau);
      } else {
        registryQaDetector.fill(HIST("Lambda/01_Kinematics/02_Final/h1f_pt"), pt);
        registryQaDetector.fill(HIST("Lambda/01_Kinematics/02_Final/h1f_eta"), eta);
        registryQaDetector.fill(HIST("Lambda/01_Kinematics/02_Final/h1f_phi"), phi);

        registryQaDetector.fill(HIST("Lambda/02_Topology/02_Final/h1f_cospa"), cosPA);
        registryQaDetector.fill(HIST("Lambda/02_Topology/02_Final/h1f_cospa_zoom"), cosPA);
        registryQaDetector.fill(HIST("Lambda/02_Topology/02_Final/h1f_dcaV0toPV"), dcaToPV);
        registryQaDetector.fill(HIST("Lambda/02_Topology/02_Final/h2f_armpod"), alpha, qtarm);
        registryQaDetector.fill(HIST("Lambda/02_Topology/02_Final/h1f_radius"), r);
        registryQaDetector.fill(HIST("Lambda/02_Topology/02_Final/h1f_dcaDaughters"), dcaDau);
        registryQaDetector.fill(HIST("Lambda/02_Topology/02_Final/h1f_ctau"), ctau);

        registryQaDetector.fill(HIST("Lambda/01_Kinematics/02_Final/h1f_rapidity"), rap);
        registryQaDetector.fill(HIST("Lambda/02_Topology/02_Final/h1f_bgctau"), bgctau);
      }
    } else {
      if (stage == LambdaQAStage::RawCandidate) {
        registryQaDetector.fill(HIST("AntiLambda/01_Kinematics/01_RawCandidate/h1f_pt"), pt);
        registryQaDetector.fill(HIST("AntiLambda/01_Kinematics/01_RawCandidate/h1f_eta"), eta);
        registryQaDetector.fill(HIST("AntiLambda/01_Kinematics/01_RawCandidate/h1f_phi"), phi);

        registryQaDetector.fill(HIST("AntiLambda/02_Topology/01_RawCandidate/h1f_cospa"), cosPA);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/01_RawCandidate/h1f_cospa_zoom"), cosPA);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/01_RawCandidate/h1f_dcaV0toPV"), dcaToPV);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/01_RawCandidate/h2f_armpod"), alpha, qtarm);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/01_RawCandidate/h1f_radius"), r);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/01_RawCandidate/h1f_dcaDaughters"), dcaDau);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/01_RawCandidate/h1f_ctau"), ctau);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/01_RawCandidate/h1f_bgctau"), bgctau);
      } else {
        registryQaDetector.fill(HIST("AntiLambda/01_Kinematics/02_Final/h1f_pt"), pt);
        registryQaDetector.fill(HIST("AntiLambda/01_Kinematics/02_Final/h1f_eta"), eta);
        registryQaDetector.fill(HIST("AntiLambda/01_Kinematics/02_Final/h1f_phi"), phi);

        registryQaDetector.fill(HIST("AntiLambda/02_Topology/02_Final/h1f_cospa"), cosPA);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/02_Final/h1f_cospa_zoom"), cosPA);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/02_Final/h1f_dcaV0toPV"), dcaToPV);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/02_Final/h2f_armpod"), alpha, qtarm);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/02_Final/h1f_radius"), r);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/02_Final/h1f_dcaDaughters"), dcaDau);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/02_Final/h1f_ctau"), ctau);

        registryQaDetector.fill(HIST("AntiLambda/01_Kinematics/02_Final/h1f_rapidity"), rap);
        registryQaDetector.fill(HIST("AntiLambda/02_Topology/02_Final/h1f_bgctau"), bgctau);
      }
    }
  }

  // ── Staged fill helper: Daughter QA ──────────────────────────────────────
  template <DaughterParent Parent, DaughterType DauType, typename TTrack>
  void fillDaughterQA(DaughterQAStage stage, TTrack const& trk)
  {
    const float p = trk.p();
    const float pt = trk.pt();
    const float dedx = trk.tpcSignal();
    const float dcaXY = trk.dcaXY();
    const float dcaZ = trk.dcaZ();

    // Determine nsigmaTPC based on compile-time DaughterType
    float nsigmaTPC = 0.0f;
    if constexpr (DauType == DaughterType::Proton || DauType == DaughterType::AntiProton) {
      nsigmaTPC = trk.tpcNSigmaPr();
    } else {
      nsigmaTPC = trk.tpcNSigmaPi();
    }

    if (stage == DaughterQAStage::Raw) {
      if constexpr (Parent == DaughterParent::Lambda && DauType == DaughterType::Proton) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_dEdx_vs_p"), p, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_dEdx_vs_pt"), pt, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_nsigmaTPC_vs_p"), p, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_nsigmaTPC_vs_pt"), pt, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_dcaZ_vs_pt"), pt, dcaZ);
        if (trk.hasTOF()) {
          const float beta = trk.beta();
          const float pTOF = trk.tofExpMom();
          const float nsigmaTOF = trk.tofNSigmaPr();
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_beta_vs_p"), p, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_beta_vs_pt"), pt, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_beta_vs_pTOF"), pTOF, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_nsigmaTOF_vs_p"), p, nsigmaTOF);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_nsigmaTOF_vs_pt"), pt, nsigmaTOF);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_nsigmaTOF_vs_pTOF"), pTOF, nsigmaTOF);
        }
        // ── Item 3: daughter tracking quality (Lambda/Proton/Raw) ──
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h1f_tpcCrossedRows"), static_cast<float>(trk.tpcNClsCrossedRows()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h1f_tpcCrossedRowsOverFindable"), trk.tpcCrossedRowsOverFindableCls());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h1f_tpcChi2NCl"), trk.tpcChi2NCl());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h1f_itsNCls"), static_cast<float>(trk.itsNCls()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h1f_itsChi2NCl"), trk.itsChi2NCl());
        // ──────────────────────────────────────────────────────────
      } else if constexpr (Parent == DaughterParent::Lambda && DauType == DaughterType::Pion) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_dEdx_vs_p"), p, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_dEdx_vs_pt"), pt, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_nsigmaTPC_vs_p"), p, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_nsigmaTPC_vs_pt"), pt, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_dcaZ_vs_pt"), pt, dcaZ);
        if (trk.hasTOF()) {
          const float beta = trk.beta();
          const float pTOF = trk.tofExpMom();
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_beta_vs_p"), p, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_beta_vs_pt"), pt, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_beta_vs_pTOF"), pTOF, beta);
        }
        // ── Item 3: daughter tracking quality (Lambda/Pion/Raw) ──
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h1f_tpcCrossedRows"), static_cast<float>(trk.tpcNClsCrossedRows()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h1f_tpcCrossedRowsOverFindable"), trk.tpcCrossedRowsOverFindableCls());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h1f_tpcChi2NCl"), trk.tpcChi2NCl());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h1f_itsNCls"), static_cast<float>(trk.itsNCls()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h1f_itsChi2NCl"), trk.itsChi2NCl());
        // ──────────────────────────────────────────────────────────
      } else if constexpr (Parent == DaughterParent::AntiLambda && DauType == DaughterType::AntiProton) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_dEdx_vs_p"), p, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_dEdx_vs_pt"), pt, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_nsigmaTPC_vs_p"), p, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_nsigmaTPC_vs_pt"), pt, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_dcaZ_vs_pt"), pt, dcaZ);
        if (trk.hasTOF()) {
          const float beta = trk.beta();
          const float pTOF = trk.tofExpMom();
          const float nsigmaTOF = trk.tofNSigmaPr();
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_beta_vs_p"), p, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_beta_vs_pt"), pt, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_beta_vs_pTOF"), pTOF, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_nsigmaTOF_vs_p"), p, nsigmaTOF);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_nsigmaTOF_vs_pt"), pt, nsigmaTOF);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_nsigmaTOF_vs_pTOF"), pTOF, nsigmaTOF);
        }
        // ── Item 3: daughter tracking quality (AntiLambda/AntiProton/Raw) ──
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h1f_tpcCrossedRows"), static_cast<float>(trk.tpcNClsCrossedRows()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h1f_tpcCrossedRowsOverFindable"), trk.tpcCrossedRowsOverFindableCls());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h1f_tpcChi2NCl"), trk.tpcChi2NCl());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h1f_itsNCls"), static_cast<float>(trk.itsNCls()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h1f_itsChi2NCl"), trk.itsChi2NCl());
        // ──────────────────────────────────────────────────────────────────
      } else if constexpr (Parent == DaughterParent::AntiLambda && DauType == DaughterType::AntiPion) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_dEdx_vs_p"), p, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_dEdx_vs_pt"), pt, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_nsigmaTPC_vs_p"), p, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_nsigmaTPC_vs_pt"), pt, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_dcaZ_vs_pt"), pt, dcaZ);
        if (trk.hasTOF()) {
          const float beta = trk.beta();
          const float pTOF = trk.tofExpMom();
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_beta_vs_p"), p, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_beta_vs_pt"), pt, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_beta_vs_pTOF"), pTOF, beta);
        }
        // ── Item 3: daughter tracking quality (AntiLambda/AntiPion/Raw) ──
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h1f_tpcCrossedRows"), static_cast<float>(trk.tpcNClsCrossedRows()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h1f_tpcCrossedRowsOverFindable"), trk.tpcCrossedRowsOverFindableCls());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h1f_tpcChi2NCl"), trk.tpcChi2NCl());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h1f_itsNCls"), static_cast<float>(trk.itsNCls()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h1f_itsChi2NCl"), trk.itsChi2NCl());
        // ──────────────────────────────────────────────────────────────
      }

      // ── Item 4: 2D TOF-matching efficiency map for V0 daughters (Raw) ──
      if constexpr (Parent == DaughterParent::Lambda && DauType == DaughterType::Proton) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/01_Raw/h2f_TOFMatchedFractionVsEtaPhi"),
                                trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      } else if constexpr (Parent == DaughterParent::Lambda && DauType == DaughterType::Pion) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/01_Raw/h2f_TOFMatchedFractionVsEtaPhi"),
                                trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      } else if constexpr (Parent == DaughterParent::AntiLambda && DauType == DaughterType::AntiProton) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/01_Raw/h2f_TOFMatchedFractionVsEtaPhi"),
                                trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      } else if constexpr (Parent == DaughterParent::AntiLambda && DauType == DaughterType::AntiPion) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/01_Raw/h2f_TOFMatchedFractionVsEtaPhi"),
                                trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      }
      // ───────────────────────────────────────────────────────────────────

    } else {
      if constexpr (Parent == DaughterParent::Lambda && DauType == DaughterType::Proton) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_dEdx_vs_p"), p, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_dEdx_vs_pt"), pt, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_nsigmaTPC_vs_p"), p, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_nsigmaTPC_vs_pt"), pt, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_dcaZ_vs_pt"), pt, dcaZ);
        if (trk.hasTOF()) {
          const float beta = trk.beta();
          const float pTOF = trk.tofExpMom();
          const float nsigmaTOF = trk.tofNSigmaPr();
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_beta_vs_p"), p, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_beta_vs_pt"), pt, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_beta_vs_pTOF"), pTOF, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_nsigmaTOF_vs_p"), p, nsigmaTOF);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_nsigmaTOF_vs_pt"), pt, nsigmaTOF);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_nsigmaTOF_vs_pTOF"), pTOF, nsigmaTOF);
        }
        // ── Item 3: daughter tracking quality (Lambda/Proton/Final) ──
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h1f_tpcCrossedRows"), static_cast<float>(trk.tpcNClsCrossedRows()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h1f_tpcCrossedRowsOverFindable"), trk.tpcCrossedRowsOverFindableCls());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h1f_tpcChi2NCl"), trk.tpcChi2NCl());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h1f_itsNCls"), static_cast<float>(trk.itsNCls()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h1f_itsChi2NCl"), trk.itsChi2NCl());
        // ──────────────────────────────────────────────────────────
      } else if constexpr (Parent == DaughterParent::Lambda && DauType == DaughterType::Pion) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_dEdx_vs_p"), p, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_dEdx_vs_pt"), pt, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_nsigmaTPC_vs_p"), p, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_nsigmaTPC_vs_pt"), pt, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_dcaZ_vs_pt"), pt, dcaZ);
        if (trk.hasTOF()) {
          const float beta = trk.beta();
          const float pTOF = trk.tofExpMom();
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_beta_vs_p"), p, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_beta_vs_pt"), pt, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_beta_vs_pTOF"), pTOF, beta);
        }
        // ── Item 3: daughter tracking quality (Lambda/Pion/Final) ──
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h1f_tpcCrossedRows"), static_cast<float>(trk.tpcNClsCrossedRows()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h1f_tpcCrossedRowsOverFindable"), trk.tpcCrossedRowsOverFindableCls());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h1f_tpcChi2NCl"), trk.tpcChi2NCl());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h1f_itsNCls"), static_cast<float>(trk.itsNCls()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h1f_itsChi2NCl"), trk.itsChi2NCl());
        // ──────────────────────────────────────────────────────────
      } else if constexpr (Parent == DaughterParent::AntiLambda && DauType == DaughterType::AntiProton) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_dEdx_vs_p"), p, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_dEdx_vs_pt"), pt, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_nsigmaTPC_vs_p"), p, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_nsigmaTPC_vs_pt"), pt, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_dcaZ_vs_pt"), pt, dcaZ);
        if (trk.hasTOF()) {
          const float beta = trk.beta();
          const float pTOF = trk.tofExpMom();
          const float nsigmaTOF = trk.tofNSigmaPr();
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_beta_vs_p"), p, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_beta_vs_pt"), pt, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_beta_vs_pTOF"), pTOF, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_nsigmaTOF_vs_p"), p, nsigmaTOF);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_nsigmaTOF_vs_pt"), pt, nsigmaTOF);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_nsigmaTOF_vs_pTOF"), pTOF, nsigmaTOF);
        }
        // ── Item 3: daughter tracking quality (AntiLambda/AntiProton/Final) ──
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h1f_tpcCrossedRows"), static_cast<float>(trk.tpcNClsCrossedRows()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h1f_tpcCrossedRowsOverFindable"), trk.tpcCrossedRowsOverFindableCls());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h1f_tpcChi2NCl"), trk.tpcChi2NCl());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h1f_itsNCls"), static_cast<float>(trk.itsNCls()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h1f_itsChi2NCl"), trk.itsChi2NCl());
        // ──────────────────────────────────────────────────────────────────
      } else if constexpr (Parent == DaughterParent::AntiLambda && DauType == DaughterType::AntiPion) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_dEdx_vs_p"), p, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_dEdx_vs_pt"), pt, dedx);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_nsigmaTPC_vs_p"), p, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_nsigmaTPC_vs_pt"), pt, nsigmaTPC);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_dcaXY_vs_pt"), pt, dcaXY);
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_dcaZ_vs_pt"), pt, dcaZ);
        if (trk.hasTOF()) {
          const float beta = trk.beta();
          const float pTOF = trk.tofExpMom();
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_beta_vs_p"), p, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_beta_vs_pt"), pt, beta);
          registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_beta_vs_pTOF"), pTOF, beta);
        }
        // ── Item 3: daughter tracking quality (AntiLambda/AntiPion/Final) ──
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h1f_tpcCrossedRows"), static_cast<float>(trk.tpcNClsCrossedRows()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h1f_tpcCrossedRowsOverFindable"), trk.tpcCrossedRowsOverFindableCls());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h1f_tpcChi2NCl"), trk.tpcChi2NCl());
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h1f_itsNCls"), static_cast<float>(trk.itsNCls()));
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h1f_itsChi2NCl"), trk.itsChi2NCl());
        // ──────────────────────────────────────────────────────────────
      }

      // ── Item 4: 2D TOF-matching efficiency map for V0 daughters (Final) ──
      if constexpr (Parent == DaughterParent::Lambda && DauType == DaughterType::Proton) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/01_Proton/02_Final/h2f_TOFMatchedFractionVsEtaPhi"),
                                trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      } else if constexpr (Parent == DaughterParent::Lambda && DauType == DaughterType::Pion) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/01_Lambda/02_Pion/02_Final/h2f_TOFMatchedFractionVsEtaPhi"),
                                trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      } else if constexpr (Parent == DaughterParent::AntiLambda && DauType == DaughterType::AntiProton) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/01_AntiProton/02_Final/h2f_TOFMatchedFractionVsEtaPhi"),
                                trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      } else if constexpr (Parent == DaughterParent::AntiLambda && DauType == DaughterType::AntiPion) {
        registryQaDetector.fill(HIST("LambdaDaughters_PID/02_AntiLambda/02_AntiPion/02_Final/h2f_TOFMatchedFractionVsEtaPhi"),
                                trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      }
      // ─────────────────────────────────────────────────────────────────────
    }
  }

  // Per-particle cache used by the pair loops (avoids recomputing per pair).
  struct PartInfo {
    int64_t gid;
    float px, py, pz, E;
    float eta, phi, y, pt;
    int idxEta;
    int idxY;
  };

  template <typename TTrack>
  PartInfo makeProtonInfo(TTrack const& t) const
  {
    PartInfo p;
    p.gid = t.globalIndex();
    p.px = t.px();
    p.py = t.py();
    p.pz = t.pz();
    p.E = std::sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz + kMassProton * kMassProton);
    p.eta = t.eta();
    p.phi = t.phi();
    p.pt = t.pt();
    p.y = protonRapidity(t);
    p.idxEta = unrolledIndex(p.eta, p.phi);
    p.idxY = unrolledIndexY(p.y, p.phi);
    return p;
  }

  void init(InitContext const&)
  {
    buildPtBins();

    // ── Common axis definitions ──────────────────────────────────────────
    AxisSpec lambdaMassAxis = {200, 1.08f, 1.15f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {100, -15.f, 15.f, "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pzAxis = {200, -10.0f, 10.0f, "p_{z} (GeV/#it{c})"};
    AxisSpec axisP = {200, 0.0f, 6.0f, "p (GeV/c)"};
    // axisPTPC: TPC inner momentum — the variable that drives the proton PID selection logic.
    // Used alongside axisP so both total-p and pTPC views are available for every diagnostic plot.
    AxisSpec axisPTPC = {200, 0.0f, 6.0f, "p_{TPC} (GeV/c)"};
    // p_TOF = trk.tofExpMom() (TracksExtra::TOFExpMom): momentum evaluated at the TOF radius,
    // i.e. the momentum actually used to build the TOF nsigma/beta response. Only valid when hasTOF()==true.
    AxisSpec axisPTOF = {200, 0.0f, 6.0f, "p_{TOF} (GeV/c)"};
    AxisSpec axisTPCdEdx = {200, 0.0f, 200.0f, "TPC dE/dx"};
    AxisSpec qaEtaAxis = {40, -0.8f, 0.8f, "#eta"}; // These should be -0.8 to 0.8 only prev it was -1 to 1
    AxisSpec qaPhiAxis = {72, 0.0f, o2::constants::math::TwoPI, "#varphi"};
    AxisSpec nsigmaAxis = {100, -10.0f, 10.0f, "n#sigma"};
    AxisSpec dcaAxis = {100, -1.0f, 1.0f, "DCAxy (cm)"};
    AxisSpec dcaZAxis = {100, -1.0f, 1.0f, "DCAz (cm)"};
    AxisSpec massAxis = {200, 1.08f, 1.15f, "M (GeV/#it{c}^{2})"};
    AxisSpec dcaAxisWide = {200, 0.0f, 5.0f, "DCA (cm)"};
    AxisSpec multAxis = {100, 0.0f, 100.0f, "FT0M Percentile (%)"};
    AxisSpec etaAxis = {kRhoEtaBins, kRhoMin, kRhoMax, "#eta"};
    AxisSpec phiAxis = {kRhoPhiBins, 0.0f, o2::constants::math::TwoPI, "#varphi"};
    AxisSpec unrolledAxis = {kRhoUnrolledBins, 0.0f, static_cast<float>(kRhoUnrolledBins), "index(#eta,#varphi)"};
    AxisSpec unrolledAxisY = {kRhoUnrolledBinsY, 0.0f, static_cast<float>(kRhoUnrolledBinsY), "index(y,#varphi)"};
    AxisSpec rapidityAxis = {kRhoYBins, kRhoYMin, kRhoYMax, "y"};
    AxisSpec yAxis = {kRhoYBins, kRhoYMin, kRhoYMax, "y"};

    // pT axis for ProtonCounts_byPID_Check: 0.15 → 6.0, 118 bins of 0.05 width
    AxisSpec pidCheckPtAxis = {118, 0.15f, 6.05f, "#it{p}_{T} (GeV/#it{c})"};
    // ─────────────────────────────────────────────────────────────────────

    // ── Other / event-level histograms ───────────────────────────────────
    registryOther.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    registryOther.add("hFT0MPercentile", "hFT0MPercentile", {HistType::kTH1F, {multAxis}});
    registryOther.add("hNSelectedPrimProtons", "Number of selected primary protons per event", {HistType::kTH1F, {{100, 0.0f, 100.0f}}});
    registryOther.add("hNSelectedPrimAntiProtons", "Number of selected primary antiprotons per event", {HistType::kTH1F, {{100, 0.0f, 100.0f}}});
    registryOther.add("hNSelLambda", "Number of selected Lambda candidates per event", {HistType::kTH1F, {{100, 0.0f, 100.0f}}});
    registryOther.add("hNSelAntiLambda", "Number of selected AntiLambda candidates per event", {HistType::kTH1F, {{100, 0.0f, 100.0f}}});
    // ── V0-daughter veto QA counters ─────────────────────────────────────
    // Proton veto counters
    registryOther.add("hNPrimProtonBeforeVeto", "Proton candidates before V0-daughter veto;N;Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPrimProtonRemovedByVeto", "Proton candidates removed by V0-daughter veto;N;Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPrimProtonAfterVeto", "Proton candidates after V0-daughter veto;N;Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    // AntiProton veto counters
    registryOther.add("hNPrimAntiProtonBeforeVeto", "AntiProton candidates before V0-daughter veto;N;Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPrimAntiProtonRemovedByVeto", "AntiProton candidates removed by V0-daughter veto;N;Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPrimAntiProtonAfterVeto", "AntiProton candidates after V0-daughter veto;N;Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    // ─────────────────────────────────────────────────────────────────────
    registryOther.add("hPtSelLambda", "Selected #Lambda p_{T} spectrum (raw counts);#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {ptAxis}});
    registryOther.add("hPtSelAntiLambda", "Selected #bar{#Lambda} p_{T} spectrum (raw counts);#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {ptAxis}});
    registryOther.add("hPtSelectedPrimProton", "Selected primary proton p_{T} spectrum (raw counts, after V0-daughter veto);#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {ptAxis}});
    registryOther.add("hPtSelectedPrimAntiProton", "Selected primary antiproton p_{T} spectrum (raw counts, after V0-daughter veto);#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {ptAxis}});
    registryOther.add("hNPairs_Lambda_PrimProton", "#Lambda-p pairs per event (after q_{inv} cut);N_{pairs};Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPairs_Lambda_PrimAntiProton", "#Lambda-#bar{p} pairs per event (after q_{inv} cut);N_{pairs};Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPairs_AntiLambda_PrimProton", "#bar{#Lambda}-p pairs per event (after q_{inv} cut);N_{pairs};Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPairs_AntiLambda_PrimAntiProton", "#bar{#Lambda}-#bar{p} pairs per event (after q_{inv} cut);N_{pairs};Events", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});

    // ── Event cutflow ─────────────────────────────────────────

    // ══════════════════════════════════════════════════════════════════════════
    // NEW STAGED QA_Detec REGISTRATION
    // ══════════════════════════════════════════════════════════════════════════

    // Axes used for the new hierarchy
    AxisSpec axisBeta = {110, 0.0f, 1.1f, "TOF #beta"};
    AxisSpec axisPtProton = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisCosPA_zoom = {200, 0.99f, 1.0f, "Cos(PA)"};
    AxisSpec axisDecayR = {100, 0.0f, 50.0f, "cm"};
    AxisSpec axisCtau = {100, 0.0f, 50.0f, "ctau (cm)"};
    AxisSpec axisCutFlow = {14, 0.5f, 14.5f, "Cut stage"};
    AxisSpec axisTOFMatchFrac = {50, 0.0f, 1.0f, "TOF matching fraction"};

    auto addProtonQAForStage = [&](const char* species, const char* stage) {
      // Kinematics
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_p", species, stage), "p;p (GeV/c);Counts", {HistType::kTH1F, {axisP}});
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_pt", species, stage), "p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {axisPtProton}});
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_eta", species, stage), "#eta;#eta;Counts", {HistType::kTH1F, {qaEtaAxis}});
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_phi", species, stage), "#phi;#phi;Counts", {HistType::kTH1F, {qaPhiAxis}});
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_rapidity", species, stage), "y;y;Counts", {HistType::kTH1F, {rapidityAxis}});

      // DetectorQuality
      registryQaDetector.add(Form("%s/02_DetectorQuality/%s/h1f_tpcCrossedRows", species, stage), "TPC Crossed Rows;N_{cr};Counts", {HistType::kTH1F, {{160, 0.0f, 160.0f, "TPC crossed rows"}}});
      registryQaDetector.add(Form("%s/02_DetectorQuality/%s/h1f_tpcCrossedRowsOverFindable", species, stage), "TPC Crossed Rows / Findable;Ratio;Counts", {HistType::kTH1F, {{100, 0.0f, 1.5f, "N_{cr}/N_{findable}"}}});
      registryQaDetector.add(Form("%s/02_DetectorQuality/%s/h1f_tpcChi2NCl", species, stage), "TPC #chi^{2}/NCl;#chi^{2}/NCl;Counts", {HistType::kTH1F, {{100, 0.0f, 10.0f, "TPC #chi^{2}/NCl"}}});
      registryQaDetector.add(Form("%s/02_DetectorQuality/%s/h1f_itsNCls", species, stage), "ITS NCls;N_{cls};Counts", {HistType::kTH1F, {{10, 0.0f, 10.0f, "ITS N_{cls}"}}});
      registryQaDetector.add(Form("%s/02_DetectorQuality/%s/h1f_itsChi2NCl", species, stage), "ITS #chi^{2}/NCl;#chi^{2}/NCl;Counts", {HistType::kTH1F, {{100, 0.0f, 40.0f, "ITS #chi^{2}/NCl"}}});
      registryQaDetector.add(Form("%s/02_DetectorQuality/%s/h2f_dcaXY_vs_pt", species, stage), "DCA_{xy} vs p_{T};p_{T};DCA_{xy}", {HistType::kTH2F, {axisPtProton, dcaAxis}});
      registryQaDetector.add(Form("%s/02_DetectorQuality/%s/h2f_dcaZ_vs_pt", species, stage), "DCA_{z} vs p_{T};p_{T};DCA_{z}", {HistType::kTH2F, {axisPtProton, dcaZAxis}});

      // TPC_QA — vs_p (total momentum) kept; vs_pTPC (TPC inner param = PID variable) added alongside
      registryQaDetector.add(Form("%s/03_TPC_QA/%s/h2f_dEdx_vs_p", species, stage), "TPC dE/dx vs p;p (GeV/c);dE/dx", {HistType::kTH2F, {axisP, axisTPCdEdx}});
      registryQaDetector.add(Form("%s/03_TPC_QA/%s/h2f_dEdx_vs_pTPC", species, stage), "TPC dE/dx vs p_{TPC};p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {axisPTPC, axisTPCdEdx}});
      registryQaDetector.add(Form("%s/03_TPC_QA/%s/h2f_dEdx_vs_pt", species, stage), "TPC dE/dx vs p_{T};p_{T} (GeV/c);dE/dx", {HistType::kTH2F, {axisPtProton, axisTPCdEdx}});
      registryQaDetector.add(Form("%s/03_TPC_QA/%s/h2f_nsigmaPr_vs_p", species, stage), "TPC n#sigma_{p} vs p;p (GeV/c);n#sigma", {HistType::kTH2F, {axisP, nsigmaAxis}});
      registryQaDetector.add(Form("%s/03_TPC_QA/%s/h2f_nsigmaPr_vs_pTPC", species, stage), "TPC n#sigma_{p} vs p_{TPC};p_{TPC} (GeV/c);n#sigma", {HistType::kTH2F, {axisPTPC, nsigmaAxis}});
      registryQaDetector.add(Form("%s/03_TPC_QA/%s/h2f_nsigmaPr_vs_pt", species, stage), "TPC n#sigma_{p} vs p_{T};p_{T} (GeV/c);n#sigma", {HistType::kTH2F, {axisPtProton, nsigmaAxis}});
      registryQaDetector.add(Form("%s/03_TPC_QA/%s/h2f_nsigmaPi_vs_pt", species, stage), "TPC n#sigma_{#pi} vs p_{T};p_{T} (GeV/c);n#sigma", {HistType::kTH2F, {axisPtProton, nsigmaAxis}});

      // TOF_QA — vs_p (total momentum) kept; vs_pTPC (TPC inner param = PID variable) added alongside
      registryQaDetector.add(Form("%s/04_TOF_QA/%s/h2f_beta_vs_p", species, stage), "TOF #beta vs p;p (GeV/c);#beta", {HistType::kTH2F, {axisP, axisBeta}});
      registryQaDetector.add(Form("%s/04_TOF_QA/%s/h2f_beta_vs_pTPC", species, stage), "TOF #beta vs p_{TPC};p_{TPC} (GeV/c);#beta", {HistType::kTH2F, {axisPTPC, axisBeta}});
      registryQaDetector.add(Form("%s/04_TOF_QA/%s/h2f_beta_vs_pt", species, stage), "TOF #beta vs p_{T};p_{T} (GeV/c);#beta", {HistType::kTH2F, {axisPtProton, axisBeta}});
      registryQaDetector.add(Form("%s/04_TOF_QA/%s/h2f_beta_vs_pTOF", species, stage), "TOF #beta vs p_{TOF};p_{TOF} (GeV/c);#beta", {HistType::kTH2F, {axisPTOF, axisBeta}});
      registryQaDetector.add(Form("%s/04_TOF_QA/%s/h2f_nsigmaPr_vs_p", species, stage), "TOF n#sigma_{p} vs p;p (GeV/c);n#sigma", {HistType::kTH2F, {axisP, nsigmaAxis}});
      registryQaDetector.add(Form("%s/04_TOF_QA/%s/h2f_nsigmaPr_vs_pTPC", species, stage), "TOF n#sigma_{p} vs p_{TPC};p_{TPC} (GeV/c);n#sigma", {HistType::kTH2F, {axisPTPC, nsigmaAxis}});
      registryQaDetector.add(Form("%s/04_TOF_QA/%s/h2f_nsigmaPr_vs_pt", species, stage), "TOF n#sigma_{p} vs p_{T};p_{T} (GeV/c);n#sigma", {HistType::kTH2F, {axisPtProton, nsigmaAxis}});
      registryQaDetector.add(Form("%s/04_TOF_QA/%s/h2f_nsigmaPr_vs_pTOF", species, stage), "TOF n#sigma_{p} vs p_{TOF};p_{TOF} (GeV/c);n#sigma", {HistType::kTH2F, {axisPTOF, nsigmaAxis}});
    };

    static constexpr std::array<const char*, 4> localProtonStages = {
      "01_RAW_AfterTrackSel", "02_TPC_PID", "03_TPC_TOF_PID", "04_Final_SelectedProton"};
    for (const auto& species : {"Proton", "AntiProton"}) {
      for (const auto& stage : localProtonStages) {
        addProtonQAForStage(species, stage);
      }
    }

    auto addLambdaQAForStage = [&](const char* species, const char* stage) {
      // Kinematics
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_pt", species, stage), "V0 p_{T} (RawCandidate = all V0s before any cut, both hypotheses);p_{T} (GeV/c);Counts", {HistType::kTH1F, {axisPtProton}});
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_eta", species, stage), "V0 #eta (RawCandidate = all V0s before any cut, both hypotheses);#eta;Counts", {HistType::kTH1F, {qaEtaAxis}});
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_phi", species, stage), "V0 #varphi (RawCandidate = all V0s before any cut, both hypotheses);#varphi;Counts", {HistType::kTH1F, {qaPhiAxis}});
      if (std::string(stage) == "02_Final") {
        registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_rapidity", species, stage), "y;y;Counts", {HistType::kTH1F, {rapidityAxis}});
      }

      // Topology
      registryQaDetector.add(Form("%s/02_Topology/%s/h1f_cospa", species, stage), "V0 cos(PA), full range;cos(PA);Counts", {HistType::kTH1F, {{100, 0.9f, 1.0f, "cos(PA)"}}});
      registryQaDetector.add(Form("%s/02_Topology/%s/h1f_cospa_zoom", species, stage), "V0 cos(PA), zoom 0.99-1.0;cos(PA);Counts", {HistType::kTH1F, {axisCosPA_zoom}});
      registryQaDetector.add(Form("%s/02_Topology/%s/h1f_dcaV0toPV", species, stage), "DCA V0 to PV;DCA (cm);Counts", {HistType::kTH1F, {dcaAxisWide}});
      registryQaDetector.add(Form("%s/02_Topology/%s/h1f_ctau", species, stage), "ctau;ctau (cm);Counts", {HistType::kTH1F, {axisCtau}});
      if (std::string(stage) == "01_RawCandidate" || std::string(stage) == "02_Final") {
        registryQaDetector.add(Form("%s/02_Topology/%s/h1f_bgctau", species, stage), "bgctau;bgctau (cm);Counts", {HistType::kTH1F, {{100, 0.0f, 200.0f, "bgctau (cm)"}}});
      }
      registryQaDetector.add(Form("%s/02_Topology/%s/h2f_armpod", species, stage), "Armenteros-Podolanski;#alpha;q_{T}", {HistType::kTH2F, {{80, -1.0f, 1.0f, "#alpha"}, {40, 0.0f, 0.4f, "q_{T}"}}});
      registryQaDetector.add(Form("%s/02_Topology/%s/h1f_radius", species, stage), "Decay Radius;R (cm);Counts", {HistType::kTH1F, {axisDecayR}});
      registryQaDetector.add(Form("%s/02_Topology/%s/h1f_dcaDaughters", species, stage), "DCA Daughters;DCA (cm);Counts", {HistType::kTH1F, {dcaAxisWide}});
    };

    for (const auto& species : {"Lambda", "AntiLambda"}) {
      registryQaDetector.add(Form("%s/InvariantMass/NoMassCut/h2f_mass_vs_pt", species), "Mass vs p_{T};p_{T};Mass", {HistType::kTH2F, {axisPtProton, lambdaMassAxis}});
      registryQaDetector.add(Form("%s/InvariantMass/MassCut/h2f_mass_vs_pt", species), "Mass vs p_{T};p_{T};Mass", {HistType::kTH2F, {axisPtProton, lambdaMassAxis}});
      // ── Item 2: K0S competing-mass hypothesis QA ───────────────────────
      registryQaDetector.add(Form("%s/InvariantMass/h1f_mK0Short_beforeCut", species),
                             "m_{K0S} before K0S rejection cut;m_{K0S} (GeV/c^{2});Counts",
                             {HistType::kTH1F, {{200, 0.40f, 0.60f, "m_{K0S} (GeV/c^{2})"}}});
      registryQaDetector.add(Form("%s/InvariantMass/h1f_mK0Short_afterLambdaSel", species),
                             "m_{K0S} for V0s passing Lambda/AntiLambda selection;m_{K0S} (GeV/c^{2});Counts",
                             {HistType::kTH1F, {{200, 0.40f, 0.60f, "m_{K0S} (GeV/c^{2})"}}});
      registryQaDetector.add(Form("%s/InvariantMass/h2f_mLambda_vs_mK0Short", species),
                             "m_{#Lambda} vs m_{K0S};m_{K0S} (GeV/c^{2});m_{#Lambda} (GeV/c^{2})",
                             {HistType::kTH2F, {{200, 0.40f, 0.60f, "m_{K0S}"}, lambdaMassAxis}});
      // ── Item 2: mass vs acceptance uniformity ──────────────────────────
      registryQaDetector.add(Form("%s/InvariantMass/h2f_mass_vs_eta", species),
                             "Mass vs #eta;#eta;m (GeV/c^{2})",
                             {HistType::kTH2F, {qaEtaAxis, lambdaMassAxis}});
      registryQaDetector.add(Form("%s/InvariantMass/h2f_mass_vs_phi", species),
                             "Mass vs #varphi;#varphi;m (GeV/c^{2})",
                             {HistType::kTH2F, {qaPhiAxis, lambdaMassAxis}});
      // ──────────────────────────────────────────────────────────────────
      for (const auto& stage : {"01_RawCandidate", "02_Final"}) {
        addLambdaQAForStage(species, stage);
      }
    }

    auto addDaughterQA = [&](const char* parent, const char* dauType, const char* stage) {
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_dEdx_vs_p", parent, dauType, stage), "TPC dE/dx vs p;p;dE/dx", {HistType::kTH2F, {axisP, axisTPCdEdx}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_dEdx_vs_pt", parent, dauType, stage), "TPC dE/dx vs p_{T};p_{T};dE/dx", {HistType::kTH2F, {axisPtProton, axisTPCdEdx}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_nsigmaTPC_vs_p", parent, dauType, stage), "TPC n#sigma vs p;p;n#sigma", {HistType::kTH2F, {axisP, nsigmaAxis}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_nsigmaTPC_vs_pt", parent, dauType, stage), "TPC n#sigma vs p_{T};p_{T};n#sigma", {HistType::kTH2F, {axisPtProton, nsigmaAxis}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_dcaXY_vs_pt", parent, dauType, stage), "DCA_{xy} vs p_{T};p_{T};DCA_{xy} (cm)", {HistType::kTH2F, {axisPtProton, dcaAxis}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_dcaZ_vs_pt", parent, dauType, stage), "DCA_{z} vs p_{T};p_{T};DCA_{z} (cm)", {HistType::kTH2F, {axisPtProton, dcaZAxis}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_beta_vs_p", parent, dauType, stage), "TOF #beta vs p;p;#beta", {HistType::kTH2F, {axisP, axisBeta}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_beta_vs_pt", parent, dauType, stage), "TOF #beta vs p_{T};p_{T};#beta", {HistType::kTH2F, {axisPtProton, axisBeta}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_beta_vs_pTOF", parent, dauType, stage), "TOF #beta vs p_{TOF};p_{TOF};#beta", {HistType::kTH2F, {axisPTOF, axisBeta}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_nsigmaTOF_vs_p", parent, dauType, stage), "TOF n#sigma vs p;p;n#sigma", {HistType::kTH2F, {axisP, nsigmaAxis}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_nsigmaTOF_vs_pt", parent, dauType, stage), "TOF n#sigma vs p_{T};p_{T};n#sigma", {HistType::kTH2F, {axisPtProton, nsigmaAxis}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_nsigmaTOF_vs_pTOF", parent, dauType, stage), "TOF n#sigma vs p_{TOF};p_{TOF};n#sigma", {HistType::kTH2F, {axisPTOF, nsigmaAxis}});
      // ── Item 3: Daughter tracking quality ────────────────────────────────
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h1f_tpcCrossedRows", parent, dauType, stage),
                             "TPC Crossed Rows;N_{cr};Counts", {HistType::kTH1F, {{160, 0.0f, 160.0f, "TPC crossed rows"}}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h1f_tpcCrossedRowsOverFindable", parent, dauType, stage),
                             "TPC Crossed Rows / Findable;Ratio;Counts", {HistType::kTH1F, {{100, 0.0f, 1.5f, "N_{cr}/N_{findable}"}}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h1f_tpcChi2NCl", parent, dauType, stage),
                             "TPC #chi^{2}/NCl;#chi^{2}/NCl;Counts", {HistType::kTH1F, {{100, 0.0f, 10.0f, "TPC #chi^{2}/NCl"}}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h1f_itsNCls", parent, dauType, stage),
                             "ITS NCls;N_{cls};Counts", {HistType::kTH1F, {{10, 0.0f, 10.0f, "ITS N_{cls}"}}});
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h1f_itsChi2NCl", parent, dauType, stage),
                             "ITS #chi^{2}/NCl;#chi^{2}/NCl;Counts", {HistType::kTH1F, {{100, 0.0f, 40.0f, "ITS #chi^{2}/NCl"}}});
      // ── Item 4: 2D TOF-matching efficiency maps (eta × phi) ──────────────
      registryQaDetector.add(Form("LambdaDaughters_PID/%s/%s/%s/h2f_TOFMatchedFractionVsEtaPhi", parent, dauType, stage),
                             "TOF matched fraction vs (#eta, #varphi);#eta;#varphi;Fraction", {HistType::kTProfile2D, {qaEtaAxis, qaPhiAxis}});
      // ─────────────────────────────────────────────────────────────────────
    };

    for (const auto& stage : {"01_Raw", "02_Final"}) {
      addDaughterQA("01_Lambda", "01_Proton", stage);
      addDaughterQA("01_Lambda", "02_Pion", stage);
      addDaughterQA("02_AntiLambda", "01_AntiProton", stage);
      addDaughterQA("02_AntiLambda", "02_AntiPion", stage);
    }

    registryQaDetector.add("Selection/TOF_Matching/hTOFMatchedFractionVsPt", "TOF-matched fraction of selected (anti)protons (before V0-daughter veto) vs p_{T};p_{T} (GeV/c);TOF-matched fraction", {HistType::kTProfile, {axisPtProton}});
    // ── Item 4: 2D TOF-matching efficiency maps (eta × phi) ──────────────────
    registryQaDetector.add("Selection/TOF_Matching/h2f_TOFMatchedFractionVsEtaPhi",
                           "TOF-matched fraction of selected (anti)protons (before V0-daughter veto) vs (#eta, #varphi);#eta;#varphi;TOF-matched fraction",
                           {HistType::kTProfile2D, {qaEtaAxis, qaPhiAxis}});
    // ─────────────────────────────────────────────────────────────────────────
    registryQaDetector.add("Selection/TOF_Matching/hCutFlow", "(Anti)proton track cutflow (both charges);Stage;Counts", {HistType::kTH1F, {axisCutFlow}});
    {
      auto hCF = registryQaDetector.get<TH1>(HIST("Selection/TOF_Matching/hCutFlow"));
      hCF->GetXaxis()->SetBinLabel(1, "All tracks");
      hCF->GetXaxis()->SetBinLabel(2, "hasTPC");
      hCF->GetXaxis()->SetBinLabel(3, "TPC crossed rows");
      hCF->GetXaxis()->SetBinLabel(4, "TPC rows/findable");
      hCF->GetXaxis()->SetBinLabel(5, "TPC chi2");
      hCF->GetXaxis()->SetBinLabel(6, "ITS Ncls");
      hCF->GetXaxis()->SetBinLabel(7, "ITS chi2");
      hCF->GetXaxis()->SetBinLabel(8, "TPC signal");
      hCF->GetXaxis()->SetBinLabel(9, "p_TPC & pT window");
      hCF->GetXaxis()->SetBinLabel(10, "|y| window");
      hCF->GetXaxis()->SetBinLabel(11, "DCA xy,z");
      hCF->GetXaxis()->SetBinLabel(12, "PID (TPC/TOF)");
      hCF->GetXaxis()->SetBinLabel(13, "V0-daughter veto");
      hCF->GetXaxis()->SetBinLabel(14, "Selected");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // EXISTING QA_Detec REGISTRATION (Preserved for backward compatibility)
    // ══════════════════════════════════════════════════════════════════════════

    // -- Event cutflow --
    registryOther.add("hEventCutflow", "Event Cutflow;Step;Counts",
                      {HistType::kTH1F, {{4, 0.5f, 4.5f}}});
    {
      auto hEC = registryOther.get<TH1>(HIST("hEventCutflow"));
      hEC->GetXaxis()->SetBinLabel(1, "sel8 + |z|<cut (Filter)");
      hEC->GetXaxis()->SetBinLabel(2, "numContrib >= 1");
      hEC->GetXaxis()->SetBinLabel(3, "centFT0M < 80");
      hEC->GetXaxis()->SetBinLabel(4, "INEL > 0");
    }

    // ── ProtonVetoQA histograms ───────────────────────────────────────────
    // Purpose: quantify the importance of the V0-daughter veto for the
    //          primary-proton sample used in rho1(p), rho2(pp), rho2(p pbar).
    // Definitions (protons only, sign > 0):
    //   h1f_nBeforeVeto    : N passing all quality+PID+DCA cuts, before veto.
    //   h1f_nRemovedByVeto : N rejected because globalIndex ∈ v0DaughterTrackIds.
    //   h1f_nAfterVeto     : N entering selectedPrimProtons (= before - removed).
    //   h1f_fractionRemoved: N_removed / N_before per event (0 when N_before = 0).
    //
    // All four histograms are filled ONCE PER EVENT.
    // Corresponding run-level totals are printed in endOfStream().
    registryProtonVetoQA.add(
      "h1f_nBeforeVeto",
      "Proton V0-daughter veto QA: N_{before veto} per event;N_{before veto};Events",
      {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryProtonVetoQA.add(
      "h1f_nRemovedByVeto",
      "Proton V0-daughter veto QA: N_{removed by veto} per event;N_{removed};Events",
      {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryProtonVetoQA.add(
      "h1f_nAfterVeto",
      "Proton V0-daughter veto QA: N_{after veto} per event;N_{after veto};Events",
      {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryProtonVetoQA.add(
      "h1f_fractionRemoved",
      "Proton V0-daughter veto QA: fraction_{removed} = N_{removed}/N_{before} per event;fraction_{removed};Events",
      {HistType::kTH1F, {{110, 0.0f, 1.1f}}});
    // ─────────────────────────────────────────────────────────────────────
    // ─────────────────────────────────────────────────────────────────────

    // ── Item 5: Event-level pileup rejection QA ───────────────────────────
    // hPileupFlags: each bin corresponds to a different event-selection bit.
    // Bin 1 = kNoSameBunchPileup, Bin 2 = kIsGoodZvtxFT0vsPV,
    // Bin 3 = kIsVertexITSTPC, Bin 4 = kIsVertexTOFmatched,
    // Bin 5 = kIsVertexTRDmatched, Bin 6 = kNoSameBunchPileup && kIsGoodZvtxFT0vsPV && kIsVertexITSTPC all set.
    registryOther.add("hPileupFlags",
                      "Event-level pile-up rejection flags;Flag;Events passing",
                      {HistType::kTH1F, {{6, 0.5f, 6.5f, "Pileup flag"}}});
    {
      auto hPU = registryOther.get<TH1>(HIST("hPileupFlags"));
      hPU->GetXaxis()->SetBinLabel(1, "NoSameBunchPileup");
      hPU->GetXaxis()->SetBinLabel(2, "GoodZvtxFT0vsPV");
      hPU->GetXaxis()->SetBinLabel(3, "VertexITSTPC");
      hPU->GetXaxis()->SetBinLabel(4, "VertexTOFmatched");
      hPU->GetXaxis()->SetBinLabel(5, "VertexTRDmatched");
      hPU->GetXaxis()->SetBinLabel(6, "NoSBPU&&GoodZvtx&&ITSTPC");
    }
    registryOther.add("h2f_NTracks_vs_Cent",
                      "N_{tracks}(|#eta|<0.8, all TracksIU, no quality cut) vs FT0M cent;FT0M (%);N_{tracks}(|#eta|<0.8)",
                      {HistType::kTH2F, {multAxis, {200, 0.0f, 200.0f, "N_{tracks}"}}});
    registryOther.add("h2f_NumContrib_vs_Cent",
                      "N_{contribs} vs FT0M percentile;FT0M (%);N_{contribs}",
                      {HistType::kTH2F, {multAxis, {200, 0.0f, 200.0f, "N_{contribs}"}}});
    // ─────────────────────────────────────────────────────────────────────

    // ── PID histograms ───────────────────────────────────────────────────
    // registryPid: vs total-p (kept unchanged — intentional total-momentum monitoring)
    registryPid.add("Proton/hTPCdEdxVsP_selected", "", {HistType::kTH2F, {axisP, axisTPCdEdx}});
    registryPid.add("Proton/hTPCdEdxVsP_TPConly", "", {HistType::kTH2F, {axisP, axisTPCdEdx}});
    registryPid.add("Proton/hTPCdEdxVsP_TPCTOF", "", {HistType::kTH2F, {axisP, axisTPCdEdx}});
    registryPid.add("AntiProton/hTPCdEdxVsP_selected", "", {HistType::kTH2F, {axisP, axisTPCdEdx}});
    registryPid.add("AntiProton/hTPCdEdxVsP_TPConly", "", {HistType::kTH2F, {axisP, axisTPCdEdx}});
    registryPid.add("AntiProton/hTPCdEdxVsP_TPCTOF", "", {HistType::kTH2F, {axisP, axisTPCdEdx}});
    // registryPid: vs pTPC (TPC inner param = PID selection variable) — added alongside for QA comparison
    registryPid.add("Proton/hTPCdEdxVsPTPC_selected", "", {HistType::kTH2F, {axisPTPC, axisTPCdEdx}});
    registryPid.add("Proton/hTPCdEdxVsPTPC_TPConly", "", {HistType::kTH2F, {axisPTPC, axisTPCdEdx}});
    registryPid.add("Proton/hTPCdEdxVsPTPC_TPCTOF", "", {HistType::kTH2F, {axisPTPC, axisTPCdEdx}});
    registryPid.add("AntiProton/hTPCdEdxVsPTPC_selected", "", {HistType::kTH2F, {axisPTPC, axisTPCdEdx}});
    registryPid.add("AntiProton/hTPCdEdxVsPTPC_TPConly", "", {HistType::kTH2F, {axisPTPC, axisTPCdEdx}});
    registryPid.add("AntiProton/hTPCdEdxVsPTPC_TPCTOF", "", {HistType::kTH2F, {axisPTPC, axisTPCdEdx}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Rho histograms ───────────────────────────────────────────────────
    registryRho.add("hEventCounter", "Events", {HistType::kTH1F, {{1, 0.0, 1.0}}});
    // eta-based rho1
    registryRho.add("h2_rho1_Proton", "#rho_{1}(p) in (#eta,#varphi) (after V0-daughter veto)", {HistType::kTH2F, {etaAxis, phiAxis}});
    registryRho.add("h2_rho1_AntiProton", "#rho_{1}(#bar{p}) in (#eta,#varphi) (after V0-daughter veto)", {HistType::kTH2F, {etaAxis, phiAxis}});
    registryRho.add("h2_rho1_Lambda", "#rho_{1}(#Lambda) in (#eta,#varphi)", {HistType::kTH2F, {etaAxis, phiAxis}});
    registryRho.add("h2_rho1_AntiLambda", "#rho_{1}(#bar{#Lambda}) in (#eta,#varphi)", {HistType::kTH2F, {etaAxis, phiAxis}});
    // rapidity-based rho1
    registryRho.add("h2_rho1_Proton_y", "#rho_{1}(p) in (y,#varphi) (after V0-daughter veto)", {HistType::kTH2F, {yAxis, phiAxis}});
    registryRho.add("h2_rho1_AntiProton_y", "#rho_{1}(#bar{p}) in (y,#varphi) (after V0-daughter veto)", {HistType::kTH2F, {yAxis, phiAxis}});
    registryRho.add("h2_rho1_Lambda_y", "#rho_{1}(#Lambda) in (y,#varphi)", {HistType::kTH2F, {yAxis, phiAxis}});
    registryRho.add("h2_rho1_AntiLambda_y", "#rho_{1}(#bar{#Lambda}) in (y,#varphi)", {HistType::kTH2F, {yAxis, phiAxis}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Autocorrelation / feed-down diagnostic histograms ────────────────
    // Purpose: check whether Lambda daughter protons enter the primary-proton
    //          sample, causing autocorrelation bias in the Lambda-p rho2.
    // Folder:  CorrelationQA/PrimaryProton_vs_LambdaDaughter/
    //
    // h2f_dEta_dPhi         : Δη vs Δφ for every (selected Lambda, selected
    //                         primary proton) pair. The primary-proton list is
    //                         ALREADY cleaned by the V0-daughter veto, so the same
    //                         track cannot appear; this checks for residual
    //                         proximity correlation (e.g. undetected feed-down).
    // h2f_dEta_dPhi_nonSameTrack : same, restricted to pairs where the primary
    //                         proton track is NOT the Lambda positive daughter.
    AxisSpec axisCorDEta = {64, -1.6f, 1.6f, "#Delta#eta"};
    AxisSpec axisCorDPhi = {72, -o2::constants::math::PI, o2::constants::math::PI, "#Delta#varphi"};
    AxisSpec axisCorDPhiDisplay = {72, -o2::constants::math::PI / 2.0f, 3.0f * o2::constants::math::PI / 2.0f, "#Delta#varphi"};
    registryCorrelationQA.add(
      "PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi",
      "Primary Proton (after V0-daughter veto) vs #Lambda Daughter Proton;#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_display",
      "Primary Proton (after V0-daughter veto) vs #Lambda Daughter Proton (display shifted #Delta#varphi);#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    registryCorrelationQA.add(
      "PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_nonSameTrack",
      "Primary Proton (after V0-daughter veto) vs #Lambda Daughter (non-same-track pairs only);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_nonSameTrack_display",
      "Primary Proton (after V0-daughter veto) vs #Lambda Daughter (non-same-track pairs only, display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Combo 2: Lambda daughter proton ↔ primary antiproton ─────────────
    // Purpose: check whether the Lambda pos-daughter proton (charge +1) could
    //          be confused with a primary antiproton (charge −1).
    // Same-track matching is impossible here (opposite charges).
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_LambdaDaughter/h2f_dEta_dPhi",
      "Primary AntiProton (after V0-daughter veto) vs #Lambda Daughter Proton;#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_LambdaDaughter/h2f_dEta_dPhi_display",
      "Primary AntiProton (after V0-daughter veto) vs #Lambda Daughter Proton (display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Combo 3: AntiLambda daughter antiproton ↔ primary proton ─────────
    // Purpose: check whether the AntiLambda neg-daughter antiproton (charge −1)
    //          could leak into the primary proton sample (charge +1).
    // Same-track matching is impossible here (opposite charges).
    registryCorrelationQA.add(
      "PrimaryProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi",
      "Primary Proton (after V0-daughter veto) vs #bar{#Lambda} Daughter AntiProton;#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_display",
      "Primary Proton (after V0-daughter veto) vs #bar{#Lambda} Daughter AntiProton (display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Combo 4: AntiLambda daughter antiproton ↔ primary antiproton ─────
    // Purpose: check whether the AntiLambda neg-daughter antiproton (charge −1)
    //          enters the primary antiproton sample, causing autocorrelation
    //          bias in the AntiLambda–pbar rho2. Analogous to Combo 1.
    // h2f_dEta_dPhi_nonSameTrack : pairs where the primary antiproton is NOT
    //          the AntiLambda neg-daughter (i.e. genuine different-track pairs).
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi",
      "Primary AntiProton (after V0-daughter veto) vs #bar{#Lambda} Daughter AntiProton;#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_display",
      "Primary AntiProton (after V0-daughter veto) vs #bar{#Lambda} Daughter AntiProton (display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_nonSameTrack",
      "Primary AntiProton (after V0-daughter veto) vs #bar{#Lambda} Daughter (non-same-track pairs only);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_nonSameTrack_display",
      "Primary AntiProton (after V0-daughter veto) vs #bar{#Lambda} Daughter (non-same-track pairs only, display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    // ─────────────────────────────────────────────────────────────────────

    // ── QA3: RawPairDensity — Δη vs Δφ for same-sign pairs ───────────────
    // Axes: Δη (400 bins, ±2.0) × Δφ (640 bins, [−π/2, 3π/2])
    AxisSpec axisQA3DEta = {400, -2.0f, 2.0f, "#Delta#eta"};
    AxisSpec axisQA3DPhi = {640, -o2::constants::math::PI / 2.0f, 3.0f * o2::constants::math::PI / 2.0f, "#Delta#varphi"};
    registryCorrelationQA.add(
      "QA3/RawPairDensity/h2f_dEta_dPhi_PP",
      "Raw pair density p-p;#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisQA3DEta, axisQA3DPhi}});
    registryCorrelationQA.add(
      "QA3/RawPairDensity/h2f_dEta_dPhi_LP",
      "Raw pair density #Lambda-p (all charge combos);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisQA3DEta, axisQA3DPhi}});
    registryCorrelationQA.add(
      "QA3/RawPairDensity/h2f_dEta_dPhi_LL",
      "Raw pair density #Lambda-#Lambda (all charge combos);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisQA3DEta, axisQA3DPhi}});
    // ─────────────────────────────────────────────────────────────────────

    // ── QA3: PairCleaning — proton-proton only ────────────────────────────
    registryCorrelationQA.add(
      "QA3/PairCleaning/h1f_pairCleaning",
      "p-p pair cleaning;Step;Counts",
      {HistType::kTH1F, {{4, 0.5f, 4.5f}}});
    // Bin labels are set after histogram registration
    {
      auto h = registryCorrelationQA.get<TH1>(HIST("QA3/PairCleaning/h1f_pairCleaning"));
      h->GetXaxis()->SetBinLabel(1, "All Candidate Pairs");
      h->GetXaxis()->SetBinLabel(2, "Rejected Same Track");
      h->GetXaxis()->SetBinLabel(3, "Rejected q_inv");
      h->GetXaxis()->SetBinLabel(4, "Accepted");
    }
    // ─────────────────────────────────────────────────────────────────────

    // ── QA3: SplitTrackQA — proton-proton close pairs ─────────────────────
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h2f_dEta_dPhi_ClosePairs",
      "Split-track QA: close p-p pairs;#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {{80, -0.02f, 0.02f, "#Delta#eta"}, {80, -0.05f, 0.05f, "#Delta#varphi"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_openingAngle_ClosePairs",
      "Split-track QA: opening angle of close p-p pairs;#theta (rad);Counts",
      {HistType::kTH1F, {{100, 0.0f, 0.2f}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_deltaPt_ClosePairs",
      "Split-track QA: |#Deltap_{T}| of close p-p pairs;|#Deltap_{T}| (GeV/#it{c});Counts",
      {HistType::kTH1F, {{100, 0.0f, 1.0f}}});
    // ─────────────────────────────────────────────────────────────────────

    // ── QA3: SplitTrackQA — q_inv and shared-cluster fraction ─────────────
    // "Before_qinvCut": filled for every accepted pair, prior to the q_inv
    // cut (kQinvCutLP), so they show the complete original q_inv distributions.
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_PP_Before_qinvCut",
      "Split-track QA: q_{inv} for p-p pairs, before q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_LP_Before_qinvCut",
      "Split-track QA: q_{inv} for #Lambda/#bar{#Lambda}-p/#bar{p} pairs (all four LP combos), before q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_LL_Before_qinvCut",
      "Split-track QA: q_{inv} for #Lambda/#bar{#Lambda}-#Lambda/#bar{#Lambda} pairs (all four LL combos), before q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    // "After_qinvCut": filled only for pairs surviving the q_inv cut
    // (q_inv > kQinvCutLP), applied uniformly to all pair types (PP, pAp, App, ApAp, LP-family, LL-family).
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_PP_After_qinvCut",
      "Split-track QA: q_{inv} for p-p pairs, after q_{inv} cut (q_{inv} > 0.01 GeV/c);q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_LP_After_qinvCut",
      "Split-track QA: q_{inv} for #Lambda/#bar{#Lambda}-p/#bar{p} pairs (all four LP combos), after q_{inv} cut (q_{inv} > 0.01 GeV/c);q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut",
      "Split-track QA: q_{inv} for #Lambda/#bar{#Lambda}-#Lambda/#bar{#Lambda} pairs (all four LL combos), after q_{inv} cut (q_{inv} > 0.01 GeV/c);q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    // q_inv Before/After QA for p-pbar, pbar-p and pbar-pbar pairs
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_pAp_Before_qinvCut",
      "Split-track QA: q_{inv} for p-#bar{p} pairs, before q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_pAp_After_qinvCut",
      "Split-track QA: q_{inv} for p-#bar{p} pairs, after q_{inv} cut (q_{inv} > 0.01 GeV/c);q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_App_Before_qinvCut",
      "Split-track QA: q_{inv} for #bar{p}-p pairs, before q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_App_After_qinvCut",
      "Split-track QA: q_{inv} for #bar{p}-p pairs, after q_{inv} cut (q_{inv} > 0.01 GeV/c);q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_ApAp_Before_qinvCut",
      "Split-track QA: q_{inv} for #bar{p}-#bar{p} pairs, before q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_ApAp_After_qinvCut",
      "Split-track QA: q_{inv} for #bar{p}-#bar{p} pairs, after q_{inv} cut (q_{inv} > 0.01 GeV/c);q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add("QA3/Kstar/h1f_kstar_Lp", "k* of #Lambda-p pairs (after q_{inv} cut);k* (GeV/c);Counts", {HistType::kTH1F, {{200, 0.0f, 2.0f}}});
    registryCorrelationQA.add("QA3/Kstar/h1f_kstar_LAp", "k* of #Lambda-#bar{p} pairs (after q_{inv} cut);k* (GeV/c);Counts", {HistType::kTH1F, {{200, 0.0f, 2.0f}}});
    registryCorrelationQA.add("QA3/Kstar/h1f_kstar_ALp", "k* of #bar{#Lambda}-p pairs (after q_{inv} cut);k* (GeV/c);Counts", {HistType::kTH1F, {{200, 0.0f, 2.0f}}});
    registryCorrelationQA.add("QA3/Kstar/h1f_kstar_ALAp", "k* of #bar{#Lambda}-#bar{p} pairs (after q_{inv} cut);k* (GeV/c);Counts", {HistType::kTH1F, {{200, 0.0f, 2.0f}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_sharedClsFraction_t1",
      "Split-track QA: TPC shared-cluster fraction, track 1 (p-p: first proton; #Lambda-p: V0 daughter) in close pairs;Fraction;Counts",
      {HistType::kTH1F, {{50, 0.0f, 1.0f}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_sharedClsFraction_t2",
      "Split-track QA: TPC shared-cluster fraction, track 2 (p-p: second proton; #Lambda-p: primary (anti)proton) in close pairs;Fraction;Counts",
      {HistType::kTH1F, {{50, 0.0f, 1.0f}}});
    // ─────────────────────────────────────────────────────────────────────

    // ── ProtonCounts_byPID_Check histograms (protons only, sign-guarded in process) ──
    h_TPC =
      registryProtonPidCheck.add<TH1>("h_TPC",
                                      "Proton candidates: TPC PID only (no TOF condition);#it{p}_{T} (GeV/#it{c});Counts",
                                      {HistType::kTH1F, {pidCheckPtAxis}});
    h_TPCandTOF =
      registryProtonPidCheck.add<TH1>("h_TPCandTOF",
                                      "Proton candidates: TPC AND TOF both required;#it{p}_{T} (GeV/#it{c});Counts",
                                      {HistType::kTH1F, {pidCheckPtAxis}});
    h_TPCorTOF =
      registryProtonPidCheck.add<TH1>("h_TPCorTOF",
                                      "Proton candidates: TPC required; if TOF present at any p, TOF n#sigma also required (differs from final selection, which requires TOF only above p_{TPC} switch);#it{p}_{T} (GeV/#it{c});Counts",
                                      {HistType::kTH1F, {pidCheckPtAxis}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Merged / pT015toMaxDefined histograms ───────────────────────────────────
    hRho2_Lp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_Lp_pT015toMaxDefined",
                           "#rho_{2}(#Lambda, p) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #Lambda (trigger);index(#eta,#varphi) of p (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_LAp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LAp_pT015toMaxDefined",
                           "#rho_{2}(#Lambda, #bar{p}) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #Lambda (trigger);index(#eta,#varphi) of #bar{p} (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ALp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALp_pT015toMaxDefined",
                           "#rho_{2}(#bar{#Lambda}, p) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #bar{#Lambda} (trigger);index(#eta,#varphi) of p (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ALAp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALAp_pT015toMaxDefined",
                           "#rho_{2}(#bar{#Lambda}, #bar{p}) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #bar{#Lambda} (trigger);index(#eta,#varphi) of #bar{p} (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_pp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_pp_pT015toMaxDefined",
                           "#rho_{2}(p, p) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of p (trigger);index(#eta,#varphi) of p (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_pAp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_pAp_pT015toMaxDefined",
                           "#rho_{2}(p, #bar{p}) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of p (trigger);index(#eta,#varphi) of #bar{p} (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_App_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_App_pT015toMaxDefined",
                           "#rho_{2}(#bar{p}, p) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #bar{p} (trigger);index(#eta,#varphi) of p (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ApAp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ApAp_pT015toMaxDefined",
                           "#rho_{2}(#bar{p}, #bar{p}) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #bar{p} (trigger);index(#eta,#varphi) of #bar{p} (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_LL_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LL_pT015toMaxDefined",
                           "#rho_{2}(#Lambda, #Lambda) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #Lambda (trigger);index(#eta,#varphi) of #Lambda (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_LAL_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LAL_pT015toMaxDefined",
                           "#rho_{2}(#Lambda, #bar{#Lambda}) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #Lambda (trigger);index(#eta,#varphi) of #bar{#Lambda} (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ALL_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALL_pT015toMaxDefined",
                           "#rho_{2}(#bar{#Lambda}, #Lambda) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #bar{#Lambda} (trigger);index(#eta,#varphi) of #Lambda (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ALAL_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALAL_pT015toMaxDefined",
                           "#rho_{2}(#bar{#Lambda}, #bar{#Lambda}) in (#eta,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(#eta,#varphi) of #bar{#Lambda} (trigger);index(#eta,#varphi) of #bar{#Lambda} (associate)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});

    // ── Rapidity merged rho2 histograms ─────────────────────────────────────
    // NOTE: these must use unrolledAxisY (kRhoYBins*kRhoPhiBins bins),
    // NOT the eta-space unrolledAxis (kRhoEtaBins*kRhoPhiBins bins).
    // GetUnrolledIndexY() below produces indices in [0, kRhoUnrolledBinsY),
    // so the booked axis must match that range or downstream R2/BF code that
    // expects y-space-sized histograms will see a dimension mismatch.
    hRho2_Lp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_Lp_y_pT015toMaxDefined", "#rho_{2}(#Lambda, p) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #Lambda (trigger);index(y,#varphi) of p (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_LAp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LAp_y_pT015toMaxDefined", "#rho_{2}(#Lambda, #bar{p}) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #Lambda (trigger);index(y,#varphi) of #bar{p} (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ALp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALp_y_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda}, p) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #bar{#Lambda} (trigger);index(y,#varphi) of p (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ALAp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALAp_y_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda}, #bar{p}) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #bar{#Lambda} (trigger);index(y,#varphi) of #bar{p} (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_pp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_pp_y_pT015toMaxDefined", "#rho_{2}(p, p) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of p (trigger);index(y,#varphi) of p (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_pAp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_pAp_y_pT015toMaxDefined", "#rho_{2}(p, #bar{p}) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of p (trigger);index(y,#varphi) of #bar{p} (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_App_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_App_y_pT015toMaxDefined", "#rho_{2}(#bar{p}, p) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #bar{p} (trigger);index(y,#varphi) of p (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ApAp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ApAp_y_pT015toMaxDefined", "#rho_{2}(#bar{p}, #bar{p}) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #bar{p} (trigger);index(y,#varphi) of #bar{p} (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_LL_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LL_y_pT015toMaxDefined", "#rho_{2}(#Lambda, #Lambda) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #Lambda (trigger);index(y,#varphi) of #Lambda (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_LAL_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LAL_y_pT015toMaxDefined", "#rho_{2}(#Lambda, #bar{#Lambda}) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #Lambda (trigger);index(y,#varphi) of #bar{#Lambda} (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ALL_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALL_y_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda}, #Lambda) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #bar{#Lambda} (trigger);index(y,#varphi) of #Lambda (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ALAL_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALAL_y_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda}, #bar{#Lambda}) in (y,#varphi) space, q_{inv} > 0.01 GeV/#it{c};index(y,#varphi) of #bar{#Lambda} (trigger);index(y,#varphi) of #bar{#Lambda} (associate)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    // ─────────────────────────────────────────────────────────────────────

    hPtSelectedPrimProton_pT015toMaxDefined =
      registryOther.add<TH1>("hPtSelectedPrimProton_pT015toMaxDefined",
                             "Selected primary proton p_{T} spectrum (integrated, 0.5-3.6 GeV/#it{c});#it{p}_{T} (GeV/#it{c});Counts",
                             {HistType::kTH1F, {ptAxis}});
    hPtSelectedPrimAntiProton_pT015toMaxDefined =
      registryOther.add<TH1>("hPtSelectedPrimAntiProton_pT015toMaxDefined",
                             "Selected primary antiproton p_{T} spectrum (integrated, 0.5-3.6 GeV/#it{c});#it{p}_{T} (GeV/#it{c});Counts",
                             {HistType::kTH1F, {ptAxis}});
    hPtSelLambda_pT015toMaxDefined =
      registryOther.add<TH1>("hPtSelLambda_pT015toMaxDefined",
                             "Selected #Lambda p_{T} spectrum (integrated, 0.6-3.6 GeV/#it{c});#it{p}_{T} (GeV/#it{c});Counts",
                             {HistType::kTH1F, {ptAxis}});
    hPtSelAntiLambda_pT015toMaxDefined =
      registryOther.add<TH1>("hPtSelAntiLambda_pT015toMaxDefined",
                             "Selected #bar{#Lambda} p_{T} spectrum (integrated, 0.6-3.6 GeV/#it{c});#it{p}_{T} (GeV/#it{c});Counts",
                             {HistType::kTH1F, {ptAxis}});

    hMassLambdaMerged =
      registryLambda.add<TH1>("hMassLambdaNoMassCut_pT015toMaxDefined",
                              "#Lambda mass (all selections except #Lambda mass window), integrated #it{p}_{T} 0.6-3.6 GeV/#it{c};#it{M}_{inv} [GeV/#it{c}^{2}];Counts",
                              {HistType::kTH1F, {lambdaMassAxis}});
    hMassAntiLambdaMerged =
      registryLambda.add<TH1>("hMassAntiLambdaNoMassCut_pT015toMaxDefined",
                              "#bar{#Lambda} mass (all selections except #Lambda mass window), integrated #it{p}_{T} 0.6-3.6 GeV/#it{c};#it{M}_{inv} [GeV/#it{c}^{2}];Counts",
                              {HistType::kTH1F, {lambdaMassAxis}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Per-pT-bin histograms (Lambda pT bins 0.6 → 6.0) ─────────────────
    hMassLambdaPtBins.clear();
    hMassAntiLambdaPtBins.clear();
    hPtPrimProton_Lambda.clear();
    hPtPrimAntiProton_Lambda.clear();
    hPtPrimProton_AntiLambda.clear();
    hPtPrimAntiProton_AntiLambda.clear();

    for (size_t i = 0; i < ptEdges.size() - 1; ++i) {
      const float ptLo = ptEdges[i];
      const float ptHi = ptEdges[i + 1];

      hPtPrimProton_Lambda.push_back(
        registryOther.add<TH1>(
          Form("hPtPrimProton_Lambda_ptBin%02zu", i),
          Form("Proton #it{p}_{T} in #Lambda-p pairs, #Lambda #it{p}_{T} in [%.1f, %.1f) GeV/#it{c};#it{p}_{T}^{p} (GeV/#it{c});Counts", ptLo, ptHi),
          {HistType::kTH1F, {ptAxis}}));
      hPtPrimAntiProton_Lambda.push_back(
        registryOther.add<TH1>(
          Form("hPtPrimAntiProton_Lambda_ptBin%02zu", i),
          Form("Antiproton #it{p}_{T} in #Lambda-#bar{p} pairs, #Lambda #it{p}_{T} in [%.1f, %.1f) GeV/#it{c};#it{p}_{T}^{#bar{p}} (GeV/#it{c});Counts", ptLo, ptHi),
          {HistType::kTH1F, {ptAxis}}));
      hPtPrimProton_AntiLambda.push_back(
        registryOther.add<TH1>(
          Form("hPtPrimProton_AntiLambda_ptBin%02zu", i),
          Form("Proton #it{p}_{T} in #bar{#Lambda}-p pairs, #bar{#Lambda} #it{p}_{T} in [%.1f, %.1f) GeV/#it{c};#it{p}_{T}^{p} (GeV/#it{c});Counts", ptLo, ptHi),
          {HistType::kTH1F, {ptAxis}}));
      hPtPrimAntiProton_AntiLambda.push_back(
        registryOther.add<TH1>(
          Form("hPtPrimAntiProton_AntiLambda_ptBin%02zu", i),
          Form("Antiproton #it{p}_{T} in #bar{#Lambda}-#bar{p} pairs, #bar{#Lambda} #it{p}_{T} in [%.1f, %.1f) GeV/#it{c};#it{p}_{T}^{#bar{p}} (GeV/#it{c});Counts", ptLo, ptHi),
          {HistType::kTH1F, {ptAxis}}));

      hMassLambdaPtBins.push_back(
        registryLambda.add<TH1>(
          Form("hMassLambdaNoMassCut_ptBin%02zu", i),
          Form("#Lambda mass (all selections except #Lambda mass window), %.1f #leq #it{p}_{T} < %.1f GeV/#it{c};#it{M}_{inv} [GeV/#it{c}^{2}];Counts", ptLo, ptHi),
          {HistType::kTH1F, {lambdaMassAxis}}));
      hMassAntiLambdaPtBins.push_back(
        registryLambda.add<TH1>(
          Form("hMassAntiLambdaNoMassCut_ptBin%02zu", i),
          Form("#bar{#Lambda} mass (all selections except #Lambda mass window), %.1f #leq #it{p}_{T} < %.1f GeV/#it{c};#it{M}_{inv} [GeV/#it{c}^{2}];Counts", ptLo, ptHi),
          {HistType::kTH1F, {lambdaMassAxis}}));
    }
    // ─────────────────────────────────────────────────────────────────────

    // ── Extended y × pT invariant-mass histograms ─────────────────────────
    // Rapidity bins: 2 bins from -0.8 to -0.6 in steps of 0.1
    // pT bins:       25 bins from 0.0 to 2.5 in steps of 0.1
    // Folder structure:
    //   Lambda_invMassExtended/Lambda/hMassLambda_yBinXX_ptBinYY
    //   Lambda_invMassExtended/AntiLambda/hMassAntiLambda_yBinXX_ptBinYY
    hMassLambdaExtended.clear();
    hMassAntiLambdaExtended.clear();
    hMassLambdaExtended.resize(kExtNyBins);
    hMassAntiLambdaExtended.resize(kExtNyBins);

    for (int iY = 0; iY < kExtNyBins; ++iY) {
      const float yLow = kExtYMin + iY * kExtYStep;
      const float yHigh = yLow + kExtYStep;
      // Last rapidity bin [-0.7,-0.6] has inclusive upper edge in filling logic
      const bool lastYBin = (iY == kExtNyBins - 1);
      const char* yUpBracket = lastYBin ? "<=" : "<";

      hMassLambdaExtended[iY].resize(kExtNptBins);
      hMassAntiLambdaExtended[iY].resize(kExtNptBins);

      for (int iPt = 0; iPt < kExtNptBins; ++iPt) {
        const float ptLow = kExtPtMin + iPt * kExtPtStep;
        const float ptHigh = ptLow + kExtPtStep;
        // Last pT bin [2.4,2.5] has inclusive upper edge in filling logic
        const bool lastPtBin = (iPt == kExtNptBins - 1);
        const char* ptUpBracket = lastPtBin ? "<=" : "<";

        // Lambda
        hMassLambdaExtended[iY][iPt] =
          registryLambdaExtended.add<TH1>(
            Form("Lambda/hMassLambda_yBin%02d_ptBin%02d", iY, iPt),
            Form("#Lambda mass (no mass-window cut) | %.1f #leq y %s %.1f | %.1f #leq #it{p}_{T} %s %.1f GeV/#it{c};#it{M}_{inv} [GeV/#it{c}^{2}];Counts",
                 yLow, yUpBracket, yHigh, ptLow, ptUpBracket, ptHigh),
            {HistType::kTH1F, {lambdaMassAxis}});

        // AntiLambda
        hMassAntiLambdaExtended[iY][iPt] =
          registryLambdaExtended.add<TH1>(
            Form("AntiLambda/hMassAntiLambda_yBin%02d_ptBin%02d", iY, iPt),
            Form("#bar{#Lambda} mass (no mass-window cut) | %.1f #leq y %s %.1f | %.1f #leq #it{p}_{T} %s %.1f GeV/#it{c};#it{M}_{inv} [GeV/#it{c}^{2}];Counts",
                 yLow, yUpBracket, yHigh, ptLow, ptUpBracket, ptHigh),
            {HistType::kTH1F, {lambdaMassAxis}});
      }
    }
    // ─────────────────────────────────────────────────────────────────────

    // ══════════════════════════════════════════════════════════════════════
    // Centrality-differential histograms → top-level folder "rho1ANDrho2_LP_CentralityBased"
    // ══════════════════════════════════════════════════════════════════════
    centSets.clear();
    centEdgesLocal = {0.f, 5.f, 10.f, 20.f, 40.f, 60.f, 80.f};
    if (cFillCentHists.value && centEdgesLocal.size() >= 2) {
      const size_t nCent = centEdgesLocal.size() - 1;
      centSets.resize(nCent);

      for (size_t ic = 0; ic < nCent; ++ic) {
        const float cLo = centEdgesLocal[ic];
        const float cHi = centEdgesLocal[ic + 1];
        const std::string dir = Form("Cent%02d_%02d", static_cast<int>(cLo), static_cast<int>(cHi));
        const std::string ctag = Form("FT0M %.0f-%.0f%%", cLo, cHi);
        auto& S = centSets[ic];

        S.hEvents = registryCent.add<TH1>((dir + "/hEventCounter").c_str(),
                                          ("Events, " + ctag).c_str(),
                                          {HistType::kTH1F, {{1, 0.0, 1.0}}});

        // ── rho1, y space (always) ──
        S.rho1_Proton_y = registryCent.add<TH2>((dir + "/h2_rho1_Proton_y").c_str(), ("#rho_{1}(p) in (y,#varphi), " + ctag).c_str(), {HistType::kTH2F, {yAxis, phiAxis}});
        S.rho1_AntiProton_y = registryCent.add<TH2>((dir + "/h2_rho1_AntiProton_y").c_str(), ("#rho_{1}(#bar{p}) in (y,#varphi), " + ctag).c_str(), {HistType::kTH2F, {yAxis, phiAxis}});
        S.rho1_Lambda_y = registryCent.add<TH2>((dir + "/h2_rho1_Lambda_y").c_str(), ("#rho_{1}(#Lambda) in (y,#varphi), " + ctag).c_str(), {HistType::kTH2F, {yAxis, phiAxis}});
        S.rho1_AntiLambda_y = registryCent.add<TH2>((dir + "/h2_rho1_AntiLambda_y").c_str(), ("#rho_{1}(#bar{#Lambda}) in (y,#varphi), " + ctag).c_str(), {HistType::kTH2F, {yAxis, phiAxis}});

        // ── rho2, y space (always) ──
        auto addY = [&](std::shared_ptr<TH2>& h, const char* name, const char* title) {
          h = registryCent.add<TH2>((dir + "/" + name).c_str(), (std::string(title) + ", " + ctag).c_str(), {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
        };
        addY(S.Lp_y, "h2_rho2_Lp_y_pT015toMaxDefined", "#rho_{2}(#Lambda,p) (y,#varphi)");
        addY(S.LAp_y, "h2_rho2_LAp_y_pT015toMaxDefined", "#rho_{2}(#Lambda,#bar{p}) (y,#varphi)");
        addY(S.ALp_y, "h2_rho2_ALp_y_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda},p) (y,#varphi)");
        addY(S.ALAp_y, "h2_rho2_ALAp_y_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda},#bar{p}) (y,#varphi)");
        addY(S.pp_y, "h2_rho2_pp_y_pT015toMaxDefined", "#rho_{2}(p,p) (y,#varphi)");
        addY(S.pAp_y, "h2_rho2_pAp_y_pT015toMaxDefined", "#rho_{2}(p,#bar{p}) (y,#varphi)");
        addY(S.App_y, "h2_rho2_App_y_pT015toMaxDefined", "#rho_{2}(#bar{p},p) (y,#varphi)");
        addY(S.ApAp_y, "h2_rho2_ApAp_y_pT015toMaxDefined", "#rho_{2}(#bar{p},#bar{p}) (y,#varphi)");
        addY(S.LL_y, "h2_rho2_LL_y_pT015toMaxDefined", "#rho_{2}(#Lambda,#Lambda) (y,#varphi)");
        addY(S.LAL_y, "h2_rho2_LAL_y_pT015toMaxDefined", "#rho_{2}(#Lambda,#bar{#Lambda}) (y,#varphi)");
        addY(S.ALL_y, "h2_rho2_ALL_y_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda},#Lambda) (y,#varphi)");
        addY(S.ALAL_y, "h2_rho2_ALAL_y_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda},#bar{#Lambda}) (y,#varphi)");

        // ── eta space (only if BOTH switches on: memory heavy) ──
        if (cFillEtaSpace.value) {
          S.rho1_Proton = registryCent.add<TH2>((dir + "/h2_rho1_Proton").c_str(), ("#rho_{1}(p) in (#eta,#varphi), " + ctag).c_str(), {HistType::kTH2F, {etaAxis, phiAxis}});
          S.rho1_AntiProton = registryCent.add<TH2>((dir + "/h2_rho1_AntiProton").c_str(), ("#rho_{1}(#bar{p}) in (#eta,#varphi), " + ctag).c_str(), {HistType::kTH2F, {etaAxis, phiAxis}});
          S.rho1_Lambda = registryCent.add<TH2>((dir + "/h2_rho1_Lambda").c_str(), ("#rho_{1}(#Lambda) in (#eta,#varphi), " + ctag).c_str(), {HistType::kTH2F, {etaAxis, phiAxis}});
          S.rho1_AntiLambda = registryCent.add<TH2>((dir + "/h2_rho1_AntiLambda").c_str(), ("#rho_{1}(#bar{#Lambda}) in (#eta,#varphi), " + ctag).c_str(), {HistType::kTH2F, {etaAxis, phiAxis}});

          auto addEta = [&](std::shared_ptr<TH2>& h, const char* name, const char* title) {
            h = registryCent.add<TH2>((dir + "/" + name).c_str(), (std::string(title) + ", " + ctag).c_str(), {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
          };
          addEta(S.Lp, "h2_rho2_Lp_pT015toMaxDefined", "#rho_{2}(#Lambda,p) (#eta,#varphi)");
          addEta(S.LAp, "h2_rho2_LAp_pT015toMaxDefined", "#rho_{2}(#Lambda,#bar{p}) (#eta,#varphi)");
          addEta(S.ALp, "h2_rho2_ALp_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda},p) (#eta,#varphi)");
          addEta(S.ALAp, "h2_rho2_ALAp_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda},#bar{p}) (#eta,#varphi)");
          addEta(S.pp, "h2_rho2_pp_pT015toMaxDefined", "#rho_{2}(p,p) (#eta,#varphi)");
          addEta(S.pAp, "h2_rho2_pAp_pT015toMaxDefined", "#rho_{2}(p,#bar{p}) (#eta,#varphi)");
          addEta(S.App, "h2_rho2_App_pT015toMaxDefined", "#rho_{2}(#bar{p},p) (#eta,#varphi)");
          addEta(S.ApAp, "h2_rho2_ApAp_pT015toMaxDefined", "#rho_{2}(#bar{p},#bar{p}) (#eta,#varphi)");
          addEta(S.LL, "h2_rho2_LL_pT015toMaxDefined", "#rho_{2}(#Lambda,#Lambda) (#eta,#varphi)");
          addEta(S.LAL, "h2_rho2_LAL_pT015toMaxDefined", "#rho_{2}(#Lambda,#bar{#Lambda}) (#eta,#varphi)");
          addEta(S.ALL, "h2_rho2_ALL_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda},#Lambda) (#eta,#varphi)");
          addEta(S.ALAL, "h2_rho2_ALAL_pT015toMaxDefined", "#rho_{2}(#bar{#Lambda},#bar{#Lambda}) (#eta,#varphi)");
        }
      }
    }
  }

  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cPPandV0ZVertexCut);

  // ── V0 topology cuts ──────────────────────────────────────────────────
  template <typename TV0>
  bool passTopoCuts(TV0 const& v0) const
  {
    if (v0.dcaV0daughters() > cV0MaxDcaDaughters.value) {
      return false;
    }
    if (v0.v0cosPA() < lambdaV0MinCosPA.value) {
      return false;
    }
    if (v0.v0radius() < lambdaV0MinDecayRadius.value) {
      return false;
    }
    if (std::abs(v0.dcav0topv()) > cV0MaxDcaToPV.value) {
      return false;
    }
    return true;
  }
  // ─────────────────────────────────────────────────────────────────────

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::MultsExtra, aod::CentFT0Ms>>::iterator const& collision,
               aod::V0Datas const& V0s,
               MyTracks const& tracks)
  {
    // ── Ordered-pair convention ────────────────────────────────────────────
    // pp   : trigger=proton1,   assoc=proton2   (self-pairs rejected)
    // pAp  : trigger=proton,    assoc=antiproton
    // App  : trigger=antiproton,assoc=proton
    // ApAp : trigger=antiproton1,assoc=antiproton2 (self-pairs rejected)
    // Lp   : trigger=Lambda,    assoc=proton    (daughter-track overlap rejected)
    // LAp  : trigger=Lambda,    assoc=antiproton
    // ALp  : trigger=AntiLambda,assoc=proton
    // ALAp : trigger=AntiLambda,assoc=antiproton
    // LL/LAL/ALL/ALAL: both V0s from respective selected lists (self-pairs rejected)
    //
    // Each channel has BOTH eta-space (unrolledIndex) and
    // rapidity-space (unrolledIndexY) versions filled in the SAME loop.
    // fillRho2Pair() handles both fills atomically.
    // ─────────────────────────────────────────────────────────────────────

    // ── Event cutflow bin 1 — enters process (after sel8 + posZ filters) ──
    registryOther.fill(HIST("hEventCutflow"), 1.0f);

    if (collision.numContrib() < 1) {
      return;
    }
    // ── Event cutflow bin 2 — after numContrib >= 1 ──────────────────────
    registryOther.fill(HIST("hEventCutflow"), 2.0f);

    static constexpr float kMaxCentFT0M = 80.0f;
    if (collision.centFT0M() > kMaxCentFT0M) {
      return;
    }
    // ── Event cutflow bin 3 — after centFT0M < 80 ────────────────────────
    registryOther.fill(HIST("hEventCutflow"), 3.0f);

    // INEL > 0: require at least one charged track in |eta| < 0.8
    int nTracksINEL = 0;
    static constexpr float kMaxInelEta = 0.8f;
    for (auto const& trk : tracks) {
      if (std::abs(trk.eta()) < kMaxInelEta) {
        ++nTracksINEL;
      }
    }
    if (nTracksINEL < 1) {
      return;
    }
    // ── Event cutflow bin 4 — after INEL > 0 ─────────────────────────────
    registryOther.fill(HIST("hEventCutflow"), 4.0f);

    registryRho.fill(HIST("hEventCounter"), 0.5);

    // ── Centrality bin for the centrality-differential histograms ─────────
    CentRho2Set* cs = nullptr;
    if (cFillCentHists.value) {
      const int iCent = centBinIndex(collision.centFT0M());
      if (iCent >= 0 && iCent < static_cast<int>(centSets.size())) {
        cs = &centSets[iCent];
        cs->hEvents->Fill(0.5);
      }
    }

    registryOther.fill(HIST("hVertexZRec"), collision.posZ());
    registryOther.fill(HIST("hFT0MPercentile"), collision.centFT0M());

    // ── Item 5: Event-level pileup rejection metrics ──────────────────────
    // Fill hPileupFlags: one bin per evsel bit, filled if the bit is set.
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      registryOther.fill(HIST("hPileupFlags"), 1.0f);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      registryOther.fill(HIST("hPileupFlags"), 2.0f);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      registryOther.fill(HIST("hPileupFlags"), 3.0f);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      registryOther.fill(HIST("hPileupFlags"), 4.0f);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      registryOther.fill(HIST("hPileupFlags"), 5.0f);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup) &&
        collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) &&
        collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      registryOther.fill(HIST("hPileupFlags"), 6.0f);
    }
    // N_{tracks} and N_{contribs} vs FT0M centrality — filled once per event
    {
      const float cent = collision.centFT0M();
      registryOther.fill(HIST("h2f_NTracks_vs_Cent"), cent, static_cast<float>(nTracksINEL));
      registryOther.fill(HIST("h2f_NumContrib_vs_Cent"), cent, static_cast<float>(collision.numContrib()));
    }
    // ─────────────────────────────────────────────────────────────────────

    // ====== Step 1: Initialise containers ======================================
    std::unordered_set<int64_t> v0DaughterTrackIds;
    v0DaughterTrackIds.reserve(64);
    std::vector<MyTracks::iterator> selectedPrimProtons;
    std::vector<MyTracks::iterator> selectedPrimAntiProtons;
    std::vector<aod::V0Datas::iterator> selectedLambdas;
    std::vector<aod::V0Datas::iterator> selectedAntiLambdas;

    int nSelLambda = 0;
    int nSelAntiLambda = 0;
    int nPairs_Lambda_PrimProton = 0;
    int nPairs_Lambda_PrimAntiProton = 0;
    int nPairs_AntiLambda_PrimProton = 0;
    int nPairs_AntiLambda_PrimAntiProton = 0;

    // ── Veto QA event counters (filled once per event) ───────────────────
    int nProtonBeforeVeto = 0;
    int nProtonRemovedVeto = 0;
    int nAntiPBeforeVeto = 0;
    int nAntiPRemovedVeto = 0;
    // ─────────────────────────────────────────────────────────────────────

    // ====== Step 2: V0 validation loop (BEFORE track loop) =================
    // Validates each V0, builds the veto set from validated daughters only,
    // fills single-particle histograms, and populates selectedLambdas/AntiLambdas.
    // Pairing with protons happens in Step 3 (after the track loop).
    constexpr float mLambdaPDG = o2::constants::physics::MassLambda;
    constexpr float mK0PDG = o2::constants::physics::MassK0Short;

    for (const auto& v0 : V0s) {
      const auto& v0PosDauCached = v0.posTrack_as<MyTracks>();
      const auto& v0NegDauCached = v0.negTrack_as<MyTracks>();

      // ── TRULY RAW QA Fills — before any V0 selection cut ─────────────
      // Placed here at the very top of the loop: before v0Type, daughter
      // kinematics, TPC PID, rapidity, pT, topology, DCA, CTau, K0S, or
      // mass-window cuts. Every V0 in the event is filled.
      // We double-fill under both the Lambda hypothesis (pos=proton, neg=pion)
      // and the AntiLambda hypothesis (pos=pion, neg=proton) since no PID
      // decision has been made yet. This gives the true combinatorial view.
      {
        const auto& rawPosDau = v0PosDauCached;
        const auto& rawNegDau = v0NegDauCached;
        const float rawCtau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * mLambdaPDG;
        const float rawBgctau = rawCtau * v0.p() / mLambdaPDG;
        // Lambda hypothesis (positive daughter assumed proton)
        fillLambdaQA<true>(LambdaQAStage::RawCandidate, v0, rawCtau, rawBgctau);
        fillDaughterQA<DaughterParent::Lambda, DaughterType::Proton>(DaughterQAStage::Raw, rawPosDau);
        fillDaughterQA<DaughterParent::Lambda, DaughterType::Pion>(DaughterQAStage::Raw, rawNegDau);
        // AntiLambda hypothesis (negative daughter assumed antiproton)
        fillLambdaQA<false>(LambdaQAStage::RawCandidate, v0, rawCtau, rawBgctau);
        fillDaughterQA<DaughterParent::AntiLambda, DaughterType::AntiProton>(DaughterQAStage::Raw, rawNegDau);
        fillDaughterQA<DaughterParent::AntiLambda, DaughterType::AntiPion>(DaughterQAStage::Raw, rawPosDau);
        // ── Item 2: K0S competing-mass hypothesis QA — TRULY RAW (before any cut) ──
        registryQaDetector.fill(HIST("Lambda/InvariantMass/h1f_mK0Short_beforeCut"), v0.mK0Short());
        registryQaDetector.fill(HIST("Lambda/InvariantMass/h2f_mLambda_vs_mK0Short"), v0.mK0Short(), v0.mLambda());
        registryQaDetector.fill(HIST("AntiLambda/InvariantMass/h1f_mK0Short_beforeCut"), v0.mK0Short());
        registryQaDetector.fill(HIST("AntiLambda/InvariantMass/h2f_mLambda_vs_mK0Short"), v0.mK0Short(), v0.mAntiLambda());
        // ─────────────────────────────────────────────────────────────────
      }
      // ─────────────────────────────────────────────────────────────────

      // ── V0 type selection (strange9 requires v0Type == 1) ─────────────
      if (v0.v0Type() != cV0TypeSelection.value) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      const auto& v0PosDau = v0PosDauCached;
      const auto& v0NegDau = v0NegDauCached;

      // ── Daughter kinematic quality cuts (from strange9) ───────────────
      if (v0PosDau.pt() <= cV0DauMinPt.value || v0NegDau.pt() <= cV0DauMinPt.value ||
          std::abs(v0PosDau.eta()) >= cV0DauMaxEta.value || std::abs(v0NegDau.eta()) >= cV0DauMaxEta.value ||
          v0PosDau.tpcNClsCrossedRows() <= cV0DauMinTpcCrossedRows.value ||
          v0NegDau.tpcNClsCrossedRows() <= cV0DauMinTpcCrossedRows.value) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // ── PID hypothesis ────────────────────────────────────────────────
      const bool v0DauPID_AsLambda =
        (std::abs(v0PosDau.tpcNSigmaPr()) < lambdaV0DaughterProtonTPCNsigma.value) &&
        (std::abs(v0NegDau.tpcNSigmaPi()) < lambdaV0DaughterPionTPCNsigma.value);
      const bool v0DauPID_AsAntiLambda =
        (std::abs(v0PosDau.tpcNSigmaPi()) < lambdaV0DaughterPionTPCNsigma.value) &&
        (std::abs(v0NegDau.tpcNSigmaPr()) < lambdaV0DaughterProtonTPCNsigma.value);

      // Reject ambiguous or unidentified V0s
      if (v0DauPID_AsLambda == v0DauPID_AsAntiLambda) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // ── TPC shared clusters cut on daughter tracks ────────────────────
      // if (v0DauPID_AsLambda) {
      //   if (v0PosDau.tpcFractionSharedCls() > 0.4f) { continue; }
      //   if (v0NegDau.tpcFractionSharedCls() > 0.4f) { continue; }
      // }
      // if (v0DauPID_AsAntiLambda) {
      //   if (v0PosDau.tpcFractionSharedCls() > 0.4f) { continue; }
      //   if (v0NegDau.tpcFractionSharedCls() > 0.4f) { continue; }
      // }
      // ─────────────────────────────────────────────────────────────────

      // ── Rapidity cuts ──────────────────────────────────────────────── // use lambdaV0MaxY configurable
      if (v0DauPID_AsLambda && std::abs(v0.rapidity(1)) > lambdaV0MaxY.value) {
        continue;
      }
      if (v0DauPID_AsAntiLambda && std::abs(v0.rapidity(2)) > lambdaV0MaxY.value) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // ── Eta cuts ──────────────────────────────────────────────────────
      // ETACUT-DISABLED (temporary): eta cut commented out, keeping rapidity cut only.
      // if (v0DauPID_AsLambda     && std::abs(v0.eta()) > lambdaV0MaxEta.value) { continue; }
      // ETACUT-DISABLED (temporary): eta cut commented out, keeping rapidity cut only.
      // if (v0DauPID_AsAntiLambda && std::abs(v0.eta()) > lambdaV0MaxEta.value) { continue; }
      // ─────────────────────────────────────────────────────────────────

      // ── pT acceptance ─────────────────────────────────────────────────
      if ((v0DauPID_AsLambda || v0DauPID_AsAntiLambda) &&
          (v0.pt() < lambdaV0MinPt.value || v0.pt() > lambdaV0MaxPt.value)) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // ── DCA of daughters to PV ────────────────────────────────────────
      const float dcaProton = v0DauPID_AsLambda ? std::abs(v0.dcapostopv())
                                                : std::abs(v0.dcanegtopv());
      const float dcaPion = v0DauPID_AsLambda ? std::abs(v0.dcanegtopv())
                                              : std::abs(v0.dcapostopv());

      // ── Lifetime (needed for Final QA fills and CTau cut below) ──────────
      const float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * mLambdaPDG;
      // FIX-3: compute lab decay length (gamma*c*tau) = ctau * |p| / m
      const float bgctau = ctau * v0.p() / mLambdaPDG;

      // ── (Raw QA fills moved to top of V0 loop — see TRULY RAW block above) ─

      // ── Apply V0 topology cuts ────────────────────────────────────────────
      if (!passTopoCuts(v0)) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // ── Apply daughter DCA-to-PV cuts ────────────────────────────────────
      if (dcaProton < lambdaV0DaughterProtonMinDCAToPV.value ||
          dcaPion < lambdaV0DaughterPionMinDCAToPV.value) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      if (ctau > lambdaV0MaxCTau.value) {
        continue;
      }

      if (std::abs(v0.mK0Short() - mK0PDG) <= lambdaV0KShortRejectMassWindow.value) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // ── Mass window selection: |m - mPDG| < lambdaV0MassWindow ─────────
      const bool selLambda = v0DauPID_AsLambda && (std::abs(v0.mLambda() - mLambdaPDG) < lambdaV0MassWindow.value);
      const bool selAntiLambda = v0DauPID_AsAntiLambda && (std::abs(v0.mAntiLambda() - mLambdaPDG) < lambdaV0MassWindow.value);
      // ── Item 2: K0S mass fill after K0S rejection + after Lambda/AntiLambda selection ──
      if (selLambda) {
        registryQaDetector.fill(HIST("Lambda/InvariantMass/h1f_mK0Short_afterLambdaSel"), v0.mK0Short());
      }
      if (selAntiLambda) {
        registryQaDetector.fill(HIST("AntiLambda/InvariantMass/h1f_mK0Short_afterLambdaSel"), v0.mK0Short());
      }
      // ─────────────────────────────────────────────────────────────────

      // ── pT bin index for per-bin histograms ──────────────────────────
      const float pt = v0.pt();
      // VS Code: commented out unused variable.
      // const float eta = v0.eta();
      // VS Code: commented out unused variable.
      // const float pz  = pt * std::sinh(eta);
      const float y = v0DauPID_AsLambda ? v0.rapidity(1) : v0.rapidity(2);
      int bin = -1;
      for (size_t i = 0; i < ptEdges.size() - 1; ++i) {
        if (pt >= ptEdges[i] && pt < ptEdges[i + 1]) {
          bin = static_cast<int>(i);
          break;
        }
      }
      // ─────────────────────────────────────────────────────────────────

      // ── NoMassCut QA fills (PID + topo only) ─────────────────────────
      // Note: fillLambdaQA/fillDaughterQA for 01_Raw were moved earlier in the
      // loop (before the daughter DCA-to-PV cut) to give a true pre-cut view.
      if (v0DauPID_AsLambda) {
        if (bin >= 0) {
          hMassLambdaPtBins[bin]->Fill(v0.mLambda());
        }
        hMassLambdaMerged->Fill(v0.mLambda());
        registryQaDetector.fill(HIST("Lambda/InvariantMass/NoMassCut/h2f_mass_vs_pt"), pt, v0.mLambda());
        // ── Item 2: mass vs eta/phi acceptance uniformity ──
        registryQaDetector.fill(HIST("Lambda/InvariantMass/h2f_mass_vs_eta"), v0.eta(), v0.mLambda());
        registryQaDetector.fill(HIST("Lambda/InvariantMass/h2f_mass_vs_phi"), v0.phi(), v0.mLambda());
        // ──────────────────────────────────────────────────
        if (selLambda) {
          registryQaDetector.fill(HIST("Lambda/InvariantMass/MassCut/h2f_mass_vs_pt"), pt, v0.mLambda());
        }

        // ── Extended y × pT Lambda invariant-mass fill ────────────────────
        // Use the same rapidity variable `y` computed above
        // and the same pT variable `pt` = v0.pt() already in scope.
        // Accept only -0.8 <= y <= -0.6 and 0.0 <= pT <= 2.5 GeV/c.
        if (y >= kExtYMin && y <= kExtYMax && pt >= kExtPtMin && pt <= kExtPtMax) {
          // Determine rapidity bin index
          int iYext = static_cast<int>((y - kExtYMin) / kExtYStep);
          if (iYext >= kExtNyBins) {
            iYext = kExtNyBins - 1;
          } // clamp upper edge of last bin
          // Determine pT bin index
          int iPtext = static_cast<int>((pt - kExtPtMin) / kExtPtStep);
          if (iPtext >= kExtNptBins) {
            iPtext = kExtNptBins - 1;
          } // clamp upper edge of last bin
          if (iYext >= 0 && iPtext >= 0) {
            hMassLambdaExtended[iYext][iPtext]->Fill(v0.mLambda());
          }
        }
        // ─────────────────────────────────────────────────────────────────
      }
      if (v0DauPID_AsAntiLambda) {
        if (bin >= 0) {
          hMassAntiLambdaPtBins[bin]->Fill(v0.mAntiLambda());
        }
        hMassAntiLambdaMerged->Fill(v0.mAntiLambda());
        registryQaDetector.fill(HIST("AntiLambda/InvariantMass/NoMassCut/h2f_mass_vs_pt"), pt, v0.mAntiLambda());
        // ── Item 2: mass vs eta/phi acceptance uniformity ──
        registryQaDetector.fill(HIST("AntiLambda/InvariantMass/h2f_mass_vs_eta"), v0.eta(), v0.mAntiLambda());
        registryQaDetector.fill(HIST("AntiLambda/InvariantMass/h2f_mass_vs_phi"), v0.phi(), v0.mAntiLambda());
        // ──────────────────────────────────────────────────
        if (selAntiLambda) {
          registryQaDetector.fill(HIST("AntiLambda/InvariantMass/MassCut/h2f_mass_vs_pt"), pt, v0.mAntiLambda());
        }

        // ── Extended y × pT AntiLambda invariant-mass fill ────────────────
        // Use the same rapidity variable `y` computed above
        // and the same pT variable `pt` = v0.pt() already in scope.
        // Accept only -0.8 <= y <= -0.6 and 0.0 <= pT <= 2.5 GeV/c.
        if (y >= kExtYMin && y <= kExtYMax && pt >= kExtPtMin && pt <= kExtPtMax) {
          int iYext = static_cast<int>((y - kExtYMin) / kExtYStep);
          if (iYext >= kExtNyBins) {
            iYext = kExtNyBins - 1;
          } // clamp upper edge of last bin
          int iPtext = static_cast<int>((pt - kExtPtMin) / kExtPtStep);
          if (iPtext >= kExtNptBins) {
            iPtext = kExtNptBins - 1;
          } // clamp upper edge of last bin
          if (iYext >= 0 && iPtext >= 0) {
            hMassAntiLambdaExtended[iYext][iPtext]->Fill(v0.mAntiLambda());
          }
        }
        // ─────────────────────────────────────────────────────────────────
      }
      // ─────────────────────────────────────────────────────────────────

      // // ── Daughter beta window veto (after NoMassCut fills, before Final) ─── //BETACUTCommented
      // // Raw and NoMassCut histograms above are filled for ALL topology-passing //BETACUTCommented
      // // candidates regardless of daughter beta. The veto here prevents V0s with //BETACUTCommented
      // // out-of-window daughter beta from reaching selLambda/selAntiLambda, the //BETACUTCommented
      // // Final QA fills, and the selectedLambdas / selectedAntiLambdas containers. //BETACUTCommented
      // // Consistent with how primary proton selection rejects incidental TOF hits //BETACUTCommented
      // // outside [cPPMinTofBeta, cPPMaxTofBeta], but using particle-dependent //BETACUTCommented
      // // beta windows for V0 daughters: [cV0protonDauMinBeta, cV0protonDauMaxBeta] //BETACUTCommented
      // // for (anti)proton daughters and >= cV0pionDauMinBeta (no upper bound) for //BETACUTCommented
      // // (anti)pion daughters. //BETACUTCommented
      // // Lambda:      pos = proton hypothesis, neg = pion hypothesis. //BETACUTCommented
      // // AntiLambda:  neg = antiproton hypothesis, pos = antipion hypothesis. //BETACUTCommented
      // if (v0DauPID_AsLambda) { //BETACUTCommented
      // // Proton daughter: reject if beta is outside [cV0protonDauMinBeta, cV0protonDauMaxBeta] //BETACUTCommented
      // if (v0PosDau.hasTOF() && (v0PosDau.beta() < cV0protonDauMinBeta.value || v0PosDau.beta() > cV0protonDauMaxBeta.value)) { continue; } //BETACUTCommented
      // // Pion daughter: only reject unphysically low beta; no upper bound (pions saturate near beta≈1) //BETACUTCommented
      // if (v0NegDau.hasTOF() && (v0NegDau.beta() < cV0pionDauMinBeta.value)) { continue; } //BETACUTCommented
      // } //BETACUTCommented
      // if (v0DauPID_AsAntiLambda) { //BETACUTCommented
      // // Antiproton daughter: reject if beta is outside [cV0protonDauMinBeta, cV0protonDauMaxBeta] //BETACUTCommented
      // if (v0NegDau.hasTOF() && (v0NegDau.beta() < cV0protonDauMinBeta.value || v0NegDau.beta() > cV0protonDauMaxBeta.value)) { continue; } //BETACUTCommented
      // // Antipion daughter: only reject unphysically low beta; no upper bound (pions saturate near beta≈1) //BETACUTCommented
      // if (v0PosDau.hasTOF() && (v0PosDau.beta() < cV0pionDauMinBeta.value)) { continue; } //BETACUTCommented
      // } //BETACUTCommented
      // ─────────────────────────────────────────────────────────────────

      // ── Step 3: Selected Lambda — populate veto set + rho1 + pT histos ─
      if (selLambda) {
        ++nSelLambda;
        selectedLambdas.push_back(v0);
        v0DaughterTrackIds.insert(v0.posTrackId());
        v0DaughterTrackIds.insert(v0.negTrackId());

        fillLambdaQA<true>(LambdaQAStage::Final, v0, ctau, bgctau);
        fillDaughterQA<DaughterParent::Lambda, DaughterType::Proton>(DaughterQAStage::Final, v0PosDau);
        fillDaughterQA<DaughterParent::Lambda, DaughterType::Pion>(DaughterQAStage::Final, v0NegDau);

        if (cFillEtaSpace.value) {
          registryRho.fill(HIST("h2_rho1_Lambda"), v0.eta(), v0.phi());
        }
        registryRho.fill(HIST("h2_rho1_Lambda_y"), v0.rapidity(1), v0.phi());
        if (cs) {
          cs->rho1_Lambda_y->Fill(v0.rapidity(1), v0.phi());
          if (cs->rho1_Lambda) {
            cs->rho1_Lambda->Fill(v0.eta(), v0.phi());
          }
        }
        registryOther.fill(HIST("hPtSelLambda"), pt);
        hPtSelLambda_pT015toMaxDefined->Fill(pt);
      }

      if (selAntiLambda) {
        ++nSelAntiLambda;
        selectedAntiLambdas.push_back(v0);
        v0DaughterTrackIds.insert(v0.posTrackId());
        v0DaughterTrackIds.insert(v0.negTrackId());

        fillLambdaQA<false>(LambdaQAStage::Final, v0, ctau, bgctau);
        fillDaughterQA<DaughterParent::AntiLambda, DaughterType::AntiProton>(DaughterQAStage::Final, v0NegDau);
        fillDaughterQA<DaughterParent::AntiLambda, DaughterType::AntiPion>(DaughterQAStage::Final, v0PosDau);

        if (cFillEtaSpace.value) {
          registryRho.fill(HIST("h2_rho1_AntiLambda"), v0.eta(), v0.phi());
        }
        registryRho.fill(HIST("h2_rho1_AntiLambda_y"), v0.rapidity(2), v0.phi());
        if (cs) {
          cs->rho1_AntiLambda_y->Fill(v0.rapidity(2), v0.phi());
          if (cs->rho1_AntiLambda) {
            cs->rho1_AntiLambda->Fill(v0.eta(), v0.phi());
          }
        }
        registryOther.fill(HIST("hPtSelAntiLambda"), pt);
        hPtSelAntiLambda_pT015toMaxDefined->Fill(pt);
      }
      // ─────────────────────────────────────────────────────────────────
    }
    // ====== End of Step 2+3 V0 validation loop ============================

    // ── Veto QA event counters (filled once per event) ───────────────────

    for (auto const& trk : tracks) {

      // ── RAW Stage QA Fills (truly before any selection — kinematics, detector quality,
      //    TPC, and TOF all filled here before any cut fires) ───────────────────────────
      if (trk.sign() > 0) {
        fillProtonKinematics<true>(ProtonQAStage::RawAfterTrackSel, trk);
        fillProtonDetectorQuality<true>(ProtonQAStage::RawAfterTrackSel, trk);
        if (trk.hasTPC()) {
          fillProtonTpcQA<true>(ProtonQAStage::RawAfterTrackSel, trk);
        }
        if (trk.hasTOF()) {
          fillProtonTofQA<true>(ProtonQAStage::RawAfterTrackSel, trk);
        }
      } else {
        fillProtonKinematics<false>(ProtonQAStage::RawAfterTrackSel, trk);
        fillProtonDetectorQuality<false>(ProtonQAStage::RawAfterTrackSel, trk);
        if (trk.hasTPC()) {
          fillProtonTpcQA<false>(ProtonQAStage::RawAfterTrackSel, trk);
        }
        if (trk.hasTOF()) {
          fillProtonTofQA<false>(ProtonQAStage::RawAfterTrackSel, trk);
        }
      }
      // ─────────────────────────────────────────────────────────────────

      // ── Track quality cuts ────────────────────────────────────────────
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 1.0f);
      if (!trk.hasTPC()) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 2.0f);
      if (trk.tpcNClsCrossedRows() < pProtonTPCMinRows.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 3.0f);
      if (trk.tpcCrossedRowsOverFindableCls() < pProtonTPCMinRowsOverFindable.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 4.0f);
      if (trk.tpcChi2NCl() > pProtonTPCMaxChi2PerCluster.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 5.0f);
      if (trk.itsNCls() < pProtonITSMinClusters.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 6.0f);
      if (trk.itsChi2NCl() > pProtonITSMaxChi2PerCluster.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 7.0f);
      if (trk.tpcSignal() < cPPMinTpcSignal.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 8.0f);
      // ─────────────────────────────────────────────────────────────────

      // ── p_TPC and pT acceptance cuts ─────────────────────────────────────────────
      const float pTPC = trk.tpcInnerParam();
      if (pTPC < pProtonMinP.value || pTPC > pProtonMaxP.value) {
        continue;
      }
      if (trk.pt() < pProtonMinPt.value || trk.pt() > pProtonMaxPt.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 9.0f);
      // ─────────────────────────────────────────────────────────────────

      // ── Proton kinematic acceptance ───────────────────────────────────
      // ETACUT-DISABLED (temporary): eta cut commented out, keeping rapidity cut only.
      // if (std::abs(trk.eta()) > pProtonMaxEta.value) { continue; }
      const float y = protonRapidity(trk);
      if (std::abs(y) > pProtonMaxY.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 10.0f);
      // ─────────────────────────────────────────────────────────────────
      // ── DCA cuts ─────────────────────────────────────────────────────
      const bool passDCA = (std::abs(trk.dcaXY()) < pProtonMaxDCAxy.value) &&
                           (std::abs(trk.dcaZ()) < pProtonMaxDCAz.value);
      if (!passDCA) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 11.0f);
      // ─────────────────────────────────────────────────────────────────

      // FIX-4: ProtonCounts_byPID_Check — protons only (sign > 0 guard added)
      if (trk.sign() > 0) {
        const float pt = trk.pt();
        const float nsTPC = std::abs(trk.tpcNSigmaPr());

        if (nsTPC < pProtonTPCNsigma.value) {
          h_TPC->Fill(pt);
        }
        if (trk.hasTOF() && nsTPC < pProtonTPCNsigma.value &&
            std::abs(trk.tofNSigmaPr()) < pProtonTOFNsigma.value) {
          h_TPCandTOF->Fill(pt);
        }
        if (nsTPC < pProtonTPCNsigma.value) {
          if (!trk.hasTOF() || std::abs(trk.tofNSigmaPr()) < pProtonTOFNsigma.value) {
            h_TPCorTOF->Fill(pt);
          }
        }
      }
      // ─────────────────────────────────────────────────────────────────

      // TPC PID pass QA fill (intermediate stage, before full selection)
      if (std::abs(trk.tpcNSigmaPr()) < pProtonTPCNsigma.value) {
        if (trk.sign() > 0) {
          fillProtonKinematics<true>(ProtonQAStage::TpcPid, trk);
          fillProtonDetectorQuality<true>(ProtonQAStage::TpcPid, trk);
          fillProtonTpcQA<true>(ProtonQAStage::TpcPid, trk);
          if (trk.hasTOF()) {
            fillProtonTofQA<true>(ProtonQAStage::TpcPid, trk);
          }
        } else {
          fillProtonKinematics<false>(ProtonQAStage::TpcPid, trk);
          fillProtonDetectorQuality<false>(ProtonQAStage::TpcPid, trk);
          fillProtonTpcQA<false>(ProtonQAStage::TpcPid, trk);
          if (trk.hasTOF()) {
            fillProtonTofQA<false>(ProtonQAStage::TpcPid, trk);
          }
        }
      }

      // TPC+TOF PID pass QA fill (intermediate stage, before full selection)
      if (std::abs(trk.tpcNSigmaPr()) < pProtonTPCNsigma.value && trk.hasTOF() && std::abs(trk.tofNSigmaPr()) < pProtonTOFNsigma.value) {
        if (trk.sign() > 0) {
          fillProtonKinematics<true>(ProtonQAStage::TpcTofPid, trk);
          fillProtonDetectorQuality<true>(ProtonQAStage::TpcTofPid, trk);
          fillProtonTpcQA<true>(ProtonQAStage::TpcTofPid, trk);
          fillProtonTofQA<true>(ProtonQAStage::TpcTofPid, trk);
        } else {
          fillProtonKinematics<false>(ProtonQAStage::TpcTofPid, trk);
          fillProtonDetectorQuality<false>(ProtonQAStage::TpcTofPid, trk);
          fillProtonTpcQA<false>(ProtonQAStage::TpcTofPid, trk);
          fillProtonTofQA<false>(ProtonQAStage::TpcTofPid, trk);
        }
      }

      // ── Main PID selection ────────────────────────────────────────────
      if (!passesPrimProtonPid(trk)) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 12.0f);

      // TOF Matching Fraction check on selected protons
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hTOFMatchedFractionVsPt"), trk.pt(), trk.hasTOF() ? 1.0f : 0.0f);
      // ── Item 4: 2D TOF-matching efficiency map for primary protons ──────
      registryQaDetector.fill(HIST("Selection/TOF_Matching/h2f_TOFMatchedFractionVsEtaPhi"),
                              trk.eta(), trk.phi(), trk.hasTOF() ? 1.0f : 0.0f);
      // ─────────────────────────────────────────────────────────────────

      // ── PID dEdx monitoring ───────────────────────────────────────────
      // vs total-p: unchanged — intentional total-momentum view
      // vs pTPC (tpcInnerParam): added alongside for proton-selection QA comparison
      if (trk.sign() > 0) {
        registryPid.fill(HIST("Proton/hTPCdEdxVsP_selected"), trk.p(), trk.tpcSignal());
        registryPid.fill(HIST("Proton/hTPCdEdxVsPTPC_selected"), trk.tpcInnerParam(), trk.tpcSignal());
        if (trk.hasTOF()) {
          registryPid.fill(HIST("Proton/hTPCdEdxVsP_TPCTOF"), trk.p(), trk.tpcSignal());
          registryPid.fill(HIST("Proton/hTPCdEdxVsPTPC_TPCTOF"), trk.tpcInnerParam(), trk.tpcSignal());
        } else {
          registryPid.fill(HIST("Proton/hTPCdEdxVsP_TPConly"), trk.p(), trk.tpcSignal());
          registryPid.fill(HIST("Proton/hTPCdEdxVsPTPC_TPConly"), trk.tpcInnerParam(), trk.tpcSignal());
        }
      } else {
        registryPid.fill(HIST("AntiProton/hTPCdEdxVsP_selected"), trk.p(), trk.tpcSignal());
        registryPid.fill(HIST("AntiProton/hTPCdEdxVsPTPC_selected"), trk.tpcInnerParam(), trk.tpcSignal());
        if (trk.hasTOF()) {
          registryPid.fill(HIST("AntiProton/hTPCdEdxVsP_TPCTOF"), trk.p(), trk.tpcSignal());
          registryPid.fill(HIST("AntiProton/hTPCdEdxVsPTPC_TPCTOF"), trk.tpcInnerParam(), trk.tpcSignal());
        } else {
          registryPid.fill(HIST("AntiProton/hTPCdEdxVsP_TPConly"), trk.p(), trk.tpcSignal());
          registryPid.fill(HIST("AntiProton/hTPCdEdxVsPTPC_TPConly"), trk.tpcInnerParam(), trk.tpcSignal());
        }
      }
      // ─────────────────────────────────────────────────────────────────

      // ── V0-daughter veto: reject tracks shared with any V0 daughter ───
      // Count before veto, check veto set, then classify or skip.
      // This MUST appear before push_back so that every pp histogram
      // (which iterates selectedPrimProtons / selectedPrimAntiProtons)
      // automatically uses the cleaned collection.
      const bool isV0Daughter = v0DaughterTrackIds.contains(trk.globalIndex());
      if (trk.sign() > 0) {
        ++nProtonBeforeVeto;
      } else {
        ++nAntiPBeforeVeto;
      }
      if (isV0Daughter) {
        // Track is a reconstructed V0 daughter — skip it.
        if (trk.sign() > 0) {
          ++nProtonRemovedVeto;
        } else {
          ++nAntiPRemovedVeto;
        }
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 13.0f);
      // ─────────────────────────────────────────────────────────────────

      // ── Classify into proton / antiproton lists ───────────────────────
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 14.0f);
      // VS Code: commented out unused variables.
      // const float pz = trk.pt() * std::sinh(trk.eta());
      // const float eta = trk.eta();
      if (trk.sign() > 0) {
        selectedPrimProtons.push_back(trk);
        fillProtonKinematics<true>(ProtonQAStage::FinalSelected, trk);
        fillProtonDetectorQuality<true>(ProtonQAStage::FinalSelected, trk);
        fillProtonTpcQA<true>(ProtonQAStage::FinalSelected, trk);
        // TOF QA only when track was accepted through the TPC+TOF branch:
        // p > PProtonTPC_TOFSwitchP guarantees hasTOF() and nSigmaTOF passed (beta window cut removed). //BETACUTCommented
        if (pTPC >= pProtonTPCTOFSwitchP.value) {
          fillProtonTofQA<true>(ProtonQAStage::FinalSelected, trk);
        }

        registryOther.fill(HIST("hPtSelectedPrimProton"), trk.pt());
        hPtSelectedPrimProton_pT015toMaxDefined->Fill(trk.pt());
      } else {
        selectedPrimAntiProtons.push_back(trk);
        fillProtonKinematics<false>(ProtonQAStage::FinalSelected, trk);
        fillProtonDetectorQuality<false>(ProtonQAStage::FinalSelected, trk);
        fillProtonTpcQA<false>(ProtonQAStage::FinalSelected, trk);
        // TOF QA only when track was accepted through the TPC+TOF branch:
        // p > pProtonTPCTOFSwitchP guarantees hasTOF() and nSigmaTOF passed (beta window cut removed). //BETACUTCommented
        if (pTPC >= pProtonTPCTOFSwitchP.value) {
          fillProtonTofQA<false>(ProtonQAStage::FinalSelected, trk);
        }

        registryOther.fill(HIST("hPtSelectedPrimAntiProton"), trk.pt());
        hPtSelectedPrimAntiProton_pT015toMaxDefined->Fill(trk.pt());
      }
      // ─────────────────────────────────────────────────────────────────
    }

    // ── Fill veto QA histograms (once per event, after the track loop) ──
    registryOther.fill(HIST("hNPrimProtonBeforeVeto"), static_cast<float>(nProtonBeforeVeto));
    registryOther.fill(HIST("hNPrimProtonRemovedByVeto"), static_cast<float>(nProtonRemovedVeto));
    registryOther.fill(HIST("hNPrimProtonAfterVeto"), static_cast<float>(nProtonBeforeVeto - nProtonRemovedVeto));
    registryOther.fill(HIST("hNPrimAntiProtonBeforeVeto"), static_cast<float>(nAntiPBeforeVeto));
    registryOther.fill(HIST("hNPrimAntiProtonRemovedByVeto"), static_cast<float>(nAntiPRemovedVeto));
    registryOther.fill(HIST("hNPrimAntiProtonAfterVeto"), static_cast<float>(nAntiPBeforeVeto - nAntiPRemovedVeto));
    // ─────────────────────────────────────────────────────────────────────

    // ── ProtonVetoQA: fill four dedicated per-event histograms ───────────
    // These are filled using PROTON (sign > 0) counters only.
    // N_before_veto  = nProtonBeforeVeto  (tracks passing quality+PID+DCA, before veto)
    // N_removed      = nProtonRemovedVeto (tracks rejected by the V0-daughter veto)
    // N_after_veto   = nProtonBeforeVeto - nProtonRemovedVeto
    //                = size of selectedPrimProtons, used by rho1(p), rho2(pp), rho2(p pbar)
    // fraction       = N_removed / N_before  (0 when N_before = 0, avoids divide-by-zero)
    {
      const int nBefore = nProtonBeforeVeto;
      const int nRemoved = nProtonRemovedVeto;
      const int nAfter = nBefore - nRemoved;
      const float fraction = (nBefore > 0)
                               ? static_cast<float>(nRemoved) / static_cast<float>(nBefore)
                               : 0.0f;
      registryProtonVetoQA.fill(HIST("h1f_nBeforeVeto"), static_cast<float>(nBefore));
      registryProtonVetoQA.fill(HIST("h1f_nRemovedByVeto"), static_cast<float>(nRemoved));
      registryProtonVetoQA.fill(HIST("h1f_nAfterVeto"), static_cast<float>(nAfter));
      registryProtonVetoQA.fill(HIST("h1f_fractionRemoved"), fraction);
      // Accumulate run-level totals for endOfStream() summary.
      runTotalProtonBeforeVeto += nBefore;
      runTotalProtonRemovedVeto += nRemoved;
      runTotalAntiPBeforeVeto += nAntiPBeforeVeto;
      runTotalAntiPRemovedVeto += nAntiPRemovedVeto;
    }
    // ─────────────────────────────────────────────────────────────────────

    // ── Single-particle density for protons (eta and y) ──────────────────
    for (auto const& primProton : selectedPrimProtons) {
      if (cFillEtaSpace.value) {
        registryRho.fill(HIST("h2_rho1_Proton"), primProton.eta(), primProton.phi());
      }
      const float yp = protonRapidity(primProton);
      registryRho.fill(HIST("h2_rho1_Proton_y"), yp, primProton.phi());
      if (cs) {
        cs->rho1_Proton_y->Fill(yp, primProton.phi());
        if (cs->rho1_Proton) {
          cs->rho1_Proton->Fill(primProton.eta(), primProton.phi());
        }
      }
    }
    for (auto const& primAntiProton : selectedPrimAntiProtons) {
      if (cFillEtaSpace.value) {
        registryRho.fill(HIST("h2_rho1_AntiProton"), primAntiProton.eta(), primAntiProton.phi());
      }
      const float yap = protonRapidity(primAntiProton);
      registryRho.fill(HIST("h2_rho1_AntiProton_y"), yap, primAntiProton.phi());
      if (cs) {
        cs->rho1_AntiProton_y->Fill(yap, primAntiProton.phi());
        if (cs->rho1_AntiProton) {
          cs->rho1_AntiProton->Fill(primAntiProton.eta(), primAntiProton.phi());
        }
      }
    }
    // ─────────────────────────────────────────────────────────────────────

    registryOther.fill(HIST("hNSelectedPrimProtons"), static_cast<float>(selectedPrimProtons.size()));
    registryOther.fill(HIST("hNSelectedPrimAntiProtons"), static_cast<float>(selectedPrimAntiProtons.size()));

    // ── Per-particle caches for the pair loops ───────────────────────────
    std::vector<PartInfo> protonInfo;
    std::vector<PartInfo> antiProtonInfo;
    protonInfo.reserve(selectedPrimProtons.size());
    antiProtonInfo.reserve(selectedPrimAntiProtons.size());
    for (auto const& t : selectedPrimProtons) {
      protonInfo.push_back(makeProtonInfo(t));
    }
    for (auto const& t : selectedPrimAntiProtons) {
      antiProtonInfo.push_back(makeProtonInfo(t));
    }

    // ── pp rho2 fills (eta and y) ───────────────────────────────────────────
    for (auto const& p1 : selectedPrimProtons) {
      const float y1_pp = protonRapidity(p1);
      for (auto const& p2 : selectedPrimProtons) {
        // QA3/PairCleaning: count all candidate pairs
        registryCorrelationQA.fill(HIST("QA3/PairCleaning/h1f_pairCleaning"), 1);
        if (p1.index() == p2.index()) {
          // QA3/PairCleaning: rejected as same track
          registryCorrelationQA.fill(HIST("QA3/PairCleaning/h1f_pairCleaning"), 2);
          continue;
        } // CHANGED from globalIndex()
        // QA3/PairCleaning: bins 3 (Rejected q_inv) and 4 (Accepted) are filled at the q_inv cut below

        // QA3/RawPairDensity: fill Δη vs Δφ for every accepted p-p pair
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          const float deta = p1.eta() - p2.eta();
          const float dphi = RecoDecay::constrainAngle((p1.phi() - p2.phi()), -PIHalf);
          registryCorrelationQA.fill(HIST("QA3/RawPairDensity/h2f_dEta_dPhi_PP"), deta, dphi);

          // QA3/SplitTrackQA: q_inv for every accepted p-p pair
          {
            const float p1mag = std::sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz());
            const float p2mag = std::sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz());
            const float E1 = std::sqrt(p1mag * p1mag + kMassProton * kMassProton);
            const float E2 = std::sqrt(p2mag * p2mag + kMassProton * kMassProton);
            const float qinvPP = computeQinv(p1.px(), p1.py(), p1.pz(), E1, p2.px(), p2.py(), p2.pz(), E2);
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_PP_Before_qinvCut"), qinvPP);
            if (!failsQinvCut(qinvPP)) {
              registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_PP_After_qinvCut"), qinvPP);
            }
          }

          // QA3/SplitTrackQA: flag close pairs for split-track diagnosis
          constexpr float cCloseDEta = 0.003f;
          constexpr float cCloseDPhi = 0.005f;
          if (std::abs(deta) < cCloseDEta && std::abs(dphi) < cCloseDPhi) {
            const float dotProd = p1.px() * p2.px() + p1.py() * p2.py() + p1.pz() * p2.pz();
            const float mag1 = std::sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz());
            const float mag2 = std::sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz());
            const float openingAngle = std::acos(std::clamp(dotProd / (mag1 * mag2), -1.f, 1.f));
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h2f_dEta_dPhi_ClosePairs"), deta, dphi);
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_openingAngle_ClosePairs"), openingAngle);
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_deltaPt_ClosePairs"), std::abs(p1.pt() - p2.pt()));
            // Shared-cluster fraction (only for close pairs)
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_sharedClsFraction_t1"), p1.tpcFractionSharedCls());
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_sharedClsFraction_t2"), p2.tpcFractionSharedCls());
          }
        }

        // q_inv cut (pp)
        {
          const float qE1 = std::sqrt(p1.px() * p1.px() + p1.py() * p1.py() + p1.pz() * p1.pz() + kMassProton * kMassProton);
          const float qE2 = std::sqrt(p2.px() * p2.px() + p2.py() * p2.py() + p2.pz() * p2.pz() + kMassProton * kMassProton);
          if (failsQinvCut(computeQinv(p1.px(), p1.py(), p1.pz(), qE1, p2.px(), p2.py(), p2.pz(), qE2))) {
            registryCorrelationQA.fill(HIST("QA3/PairCleaning/h1f_pairCleaning"), 3);
            continue;
          }
          registryCorrelationQA.fill(HIST("QA3/PairCleaning/h1f_pairCleaning"), 4);
        }
        const int idx1 = unrolledIndex(p1.eta(), p1.phi());
        const int idx2 = unrolledIndex(p2.eta(), p2.phi());
        const float y2_pp = protonRapidity(p2);
        const int idxY1 = unrolledIndexY(y1_pp, p1.phi());
        const int idxY2 = unrolledIndexY(y2_pp, p2.phi());
        // Physics pair correlations (rho2): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (cFillEtaSpace.value && idx1 >= 0 && idx2 >= 0) {
          hRho2_pp_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_pp_y_pT015toMaxDefined->Fill(idxY1, idxY2);
        }
        if (cs) {
          fillIdx(cs->pp, idx1, idx2);
          fillIdx(cs->pp_y, idxY1, idxY2);
        }
      }
    }
    for (auto const& a : protonInfo) {
      for (auto const& b : antiProtonInfo) {
        const float qinvPAp = computeQinv(a.px, a.py, a.pz, a.E, b.px, b.py, b.pz, b.E);
        registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_pAp_Before_qinvCut"), qinvPAp);
        if (failsQinvCut(qinvPAp)) {
          continue;
        }
        registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_pAp_After_qinvCut"), qinvPAp);
        if (cFillEtaSpace.value && a.idxEta >= 0 && b.idxEta >= 0) {
          hRho2_pAp_pT015toMaxDefined->Fill(a.idxEta, b.idxEta);
        }
        if (a.idxY >= 0 && b.idxY >= 0) {
          hRho2_pAp_y_pT015toMaxDefined->Fill(a.idxY, b.idxY);
        }
        if (cs) {
          fillIdx(cs->pAp, a.idxEta, b.idxEta);
          fillIdx(cs->pAp_y, a.idxY, b.idxY);
        }
      }
    }
    for (auto const& a : antiProtonInfo) {
      for (auto const& b : protonInfo) {
        const float qinvApp = computeQinv(a.px, a.py, a.pz, a.E, b.px, b.py, b.pz, b.E);
        registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_App_Before_qinvCut"), qinvApp);
        if (failsQinvCut(qinvApp)) {
          continue;
        }
        registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_App_After_qinvCut"), qinvApp);
        if (cFillEtaSpace.value && a.idxEta >= 0 && b.idxEta >= 0) {
          hRho2_App_pT015toMaxDefined->Fill(a.idxEta, b.idxEta);
        }
        if (a.idxY >= 0 && b.idxY >= 0) {
          hRho2_App_y_pT015toMaxDefined->Fill(a.idxY, b.idxY);
        }
        if (cs) {
          fillIdx(cs->App, a.idxEta, b.idxEta);
          fillIdx(cs->App_y, a.idxY, b.idxY);
        }
      }
    }
    for (auto const& a : antiProtonInfo) {
      for (auto const& b : antiProtonInfo) {
        if (a.gid == b.gid) {
          continue;
        }
        const float qinvApAp = computeQinv(a.px, a.py, a.pz, a.E, b.px, b.py, b.pz, b.E);
        registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_ApAp_Before_qinvCut"), qinvApAp);
        if (failsQinvCut(qinvApAp)) {
          continue;
        }
        registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_ApAp_After_qinvCut"), qinvApAp);
        if (cFillEtaSpace.value && a.idxEta >= 0 && b.idxEta >= 0) {
          hRho2_ApAp_pT015toMaxDefined->Fill(a.idxEta, b.idxEta);
        }
        if (a.idxY >= 0 && b.idxY >= 0) {
          hRho2_ApAp_y_pT015toMaxDefined->Fill(a.idxY, b.idxY);
        }
        if (cs) {
          fillIdx(cs->ApAp, a.idxEta, b.idxEta);
          fillIdx(cs->ApAp_y, a.idxY, b.idxY);
        }
      }
    }
    // ─────────────────────────────────────────────────────────────────────

    // ====== Step 5: Lambda-proton and AntiLambda-proton pairing ============
    // Both selectedLambdas/AntiLambdas and selectedPrimProtons/AntiProtons are
    // now fully populated. Veto set was built from validated daughters only.

    for (const auto& v0 : selectedLambdas) {
      const auto& v0PosDau = v0.posTrack_as<MyTracks>();
      const float pt = v0.pt();
      int bin = -1;
      for (size_t i = 0; i < ptEdges.size() - 1; ++i) {
        if (pt >= ptEdges[i] && pt < ptEdges[i + 1]) {
          bin = static_cast<int>(i);
          break;
        }
      }

      for (auto const& primProton : selectedPrimProtons) {
        if (primProton.globalIndex() == v0.posTrackId() || primProton.globalIndex() == v0.negTrackId()) {
          continue;
        }
        // QA3/SplitTrackQA: q_inv for every accepted Lambda-proton pair
        float qinvLP = 0.0f;
        {
          const float v0pmag = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
          const float pprmag = std::sqrt(primProton.px() * primProton.px() + primProton.py() * primProton.py() + primProton.pz() * primProton.pz());
          const float Ev0 = std::sqrt(v0pmag * v0pmag + kMassLambda * kMassLambda);
          const float Eppr = std::sqrt(pprmag * pprmag + kMassProton * kMassProton);
          qinvLP = computeQinv(v0.px(), v0.py(), v0.pz(), Ev0, primProton.px(), primProton.py(), primProton.pz(), Eppr);
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LP_Before_qinvCut"), qinvLP);
        }
        // QA3/SplitTrackQA: q_inv Before and After cut
        if (!failsQinvCut(qinvLP)) {
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LP_After_qinvCut"), qinvLP);
        }
        // QA3/SplitTrackQA: close-pair block for daughter-overlap diagnosed split tracks
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          const float deta_lp = v0.eta() - primProton.eta();
          const float dphi_lp = RecoDecay::constrainAngle((v0.phi() - primProton.phi()), -PIHalf);
          registryCorrelationQA.fill(
            HIST("QA3/RawPairDensity/h2f_dEta_dPhi_LP"),
            deta_lp, dphi_lp);
          constexpr float cCloseDEta = 0.003f;
          constexpr float cCloseDPhi = 0.005f;
          if (std::abs(deta_lp) < cCloseDEta && std::abs(dphi_lp) < cCloseDPhi) {
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_sharedClsFraction_t1"), v0PosDau.tpcFractionSharedCls());
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_sharedClsFraction_t2"), primProton.tpcFractionSharedCls());
          }
        }
        const int idxLambda = unrolledIndex(v0.eta(), v0.phi());
        const int idxPrimP = unrolledIndex(primProton.eta(), primProton.phi());
        const float yL = v0.rapidity(1);
        const float yP = protonRapidity(primProton);
        const int idxYL = unrolledIndexY(yL, v0.phi());
        const int idxYP = unrolledIndexY(yP, primProton.phi());
        // Physics pair correlations (rho2 / pT / pair-count): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (failsQinvCut(qinvLP)) {
          continue;
        }
        registryCorrelationQA.fill(HIST("QA3/Kstar/h1f_kstar_Lp"),
                                   computeKstar(v0.px(), v0.py(), v0.pz(), kMassLambda,
                                                primProton.px(), primProton.py(), primProton.pz(), kMassProton));
        if (idxLambda >= 0 && idxPrimP >= 0) {
          if (bin >= 0) {
            hPtPrimProton_Lambda[bin]->Fill(primProton.pt());
          }
          if (cFillEtaSpace.value) {
            hRho2_Lp_pT015toMaxDefined->Fill(idxLambda, idxPrimP);
          }
        }
        if (idxYL >= 0 && idxYP >= 0) {
          hRho2_Lp_y_pT015toMaxDefined->Fill(idxYL, idxYP);
        }
        if (cs) {
          fillIdx(cs->Lp, idxLambda, idxPrimP);
          fillIdx(cs->Lp_y, idxYL, idxYP);
        }
        ++nPairs_Lambda_PrimProton;
      }
      for (auto const& primAntiProton : selectedPrimAntiProtons) {
        if (primAntiProton.globalIndex() == v0.posTrackId() || primAntiProton.globalIndex() == v0.negTrackId()) {
          continue;
        }
        // QA3/SplitTrackQA: q_inv for every accepted Lambda-antiproton pair
        float qinvLAp = 0.0f;
        {
          const float v0pmag2 = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
          const float paprmag = std::sqrt(primAntiProton.px() * primAntiProton.px() + primAntiProton.py() * primAntiProton.py() + primAntiProton.pz() * primAntiProton.pz());
          const float Ev02 = std::sqrt(v0pmag2 * v0pmag2 + kMassLambda * kMassLambda);
          const float Epapr = std::sqrt(paprmag * paprmag + kMassProton * kMassProton);
          qinvLAp = computeQinv(v0.px(), v0.py(), v0.pz(), Ev02, primAntiProton.px(), primAntiProton.py(), primAntiProton.pz(), Epapr);
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LP_Before_qinvCut"), qinvLAp);
        }
        // QA3/SplitTrackQA: q_inv Before and After cut
        if (!failsQinvCut(qinvLAp)) {
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LP_After_qinvCut"), qinvLAp);
        }
        // QA3/RawPairDensity: Lambda-antiproton → same LP histogram
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          registryCorrelationQA.fill(
            HIST("QA3/RawPairDensity/h2f_dEta_dPhi_LP"),
            v0.eta() - primAntiProton.eta(),
            RecoDecay::constrainAngle((v0.phi() - primAntiProton.phi()), -PIHalf));
        }
        const int idxLambda = unrolledIndex(v0.eta(), v0.phi());
        const int idxPrimPbar = unrolledIndex(primAntiProton.eta(), primAntiProton.phi());
        const float yL2 = v0.rapidity(1);
        const float yAp = protonRapidity(primAntiProton);
        const int idxYL2 = unrolledIndexY(yL2, v0.phi());
        const int idxYAp = unrolledIndexY(yAp, primAntiProton.phi());
        // Physics pair correlations (rho2 / pT / pair-count): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (failsQinvCut(qinvLAp)) {
          continue;
        }
        registryCorrelationQA.fill(HIST("QA3/Kstar/h1f_kstar_LAp"),
                                   computeKstar(v0.px(), v0.py(), v0.pz(), kMassLambda,
                                                primAntiProton.px(), primAntiProton.py(), primAntiProton.pz(), kMassProton));
        if (idxLambda >= 0 && idxPrimPbar >= 0) {
          if (bin >= 0) {
            hPtPrimAntiProton_Lambda[bin]->Fill(primAntiProton.pt());
          }
          if (cFillEtaSpace.value) {
            hRho2_LAp_pT015toMaxDefined->Fill(idxLambda, idxPrimPbar);
          }
        }
        if (idxYL2 >= 0 && idxYAp >= 0) {
          hRho2_LAp_y_pT015toMaxDefined->Fill(idxYL2, idxYAp);
        }
        if (cs) {
          fillIdx(cs->LAp, idxLambda, idxPrimPbar);
          fillIdx(cs->LAp_y, idxYL2, idxYAp);
        }
        ++nPairs_Lambda_PrimAntiProton;
      }

      // ── Autocorrelation diagnostic: primary proton vs Lambda daughter proton──
      // Lambda proton daughter = v0PosDau (positive track, identified as proton
      // by |tpcNSigmaPr| < PProtonTPCNsigma, already in scope above).
      // Loop over ALL selected primary protons. NOTE: selectedPrimProtons is already
      // cleaned by the V0-daughter veto, so the same-track case cannot occur here;
      // this only checks for residual proximity correlations (e.g. undetected feed-down).
      // Fills inclusive h2f_dEta_dPhi for all pairs and
      // h2f_dEta_dPhi_nonSameTrack for every pair in this cleaned sample.
      {
        const float etaDau = v0PosDau.eta();
        const float phiDau = v0PosDau.phi();
        for (auto const& primProton : selectedPrimProtons) {
          const float dEta = primProton.eta() - etaDau;
          float dPhi = primProton.phi() - phiDau;
          // Wrap Δφ into [-π, +π]
          dPhi = RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI);

          // Always fill the inclusive histogram.
          registryCorrelationQA.fill(
            HIST("PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi"),
            dEta, dPhi);
          registryCorrelationQA.fill(
            HIST("PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_display"),
            dEta, displayDeltaPhi(dPhi));

          // Fill non-same-track histogram when primary proton ≠ Lambda pos-daughter.
          if (primProton.globalIndex() != v0.posTrackId()) {
            registryCorrelationQA.fill(
              HIST("PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_nonSameTrack"),
              dEta, dPhi);
            registryCorrelationQA.fill(
              HIST("PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_nonSameTrack_display"),
              dEta, displayDeltaPhi(dPhi));
          }
        }
      }
      // ─────────────────────────────────────────────────────────────────

      // ── Combo 2: Lambda daughter proton ↔ primary antiproton ─────────
      // The Lambda pos-daughter (charge +1) CANNOT appear in the primary
      // antiproton list (charge −1), so same-track matching is impossible.
      // We fill only the inclusive Δη vs Δφ histogram.
      {
        const float etaDauL = v0PosDau.eta();
        const float phiDauL = v0PosDau.phi();
        for (auto const& primAntiProton : selectedPrimAntiProtons) {
          const float dEta2 = primAntiProton.eta() - etaDauL;
          float dPhi2 = primAntiProton.phi() - phiDauL;
          dPhi2 = RecoDecay::constrainAngle(dPhi2, -o2::constants::math::PI);
          registryCorrelationQA.fill(
            HIST("PrimaryAntiProton_vs_LambdaDaughter/h2f_dEta_dPhi"),
            dEta2, dPhi2);
          registryCorrelationQA.fill(
            HIST("PrimaryAntiProton_vs_LambdaDaughter/h2f_dEta_dPhi_display"),
            dEta2, displayDeltaPhi(dPhi2));
        }
      }
      // ─────────────────────────────────────────────────────────────────
    }

    for (const auto& v0 : selectedAntiLambdas) {
      const auto& v0NegDau = v0.negTrack_as<MyTracks>();
      const float pt = v0.pt();
      int bin = -1;
      for (size_t i = 0; i < ptEdges.size() - 1; ++i) {
        if (pt >= ptEdges[i] && pt < ptEdges[i + 1]) {
          bin = static_cast<int>(i);
          break;
        }
      }

      for (auto const& primProton : selectedPrimProtons) {
        if (primProton.globalIndex() == v0.posTrackId() || primProton.globalIndex() == v0.negTrackId()) {
          continue;
        }
        // QA3/SplitTrackQA: q_inv for every accepted AntiLambda-proton pair
        float qinvALp = 0.0f;
        {
          const float alpmag = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
          const float pprmag3 = std::sqrt(primProton.px() * primProton.px() + primProton.py() * primProton.py() + primProton.pz() * primProton.pz());
          const float EalL = std::sqrt(alpmag * alpmag + kMassLambda * kMassLambda);
          const float Eppr3 = std::sqrt(pprmag3 * pprmag3 + kMassProton * kMassProton);
          qinvALp = computeQinv(v0.px(), v0.py(), v0.pz(), EalL, primProton.px(), primProton.py(), primProton.pz(), Eppr3);
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LP_Before_qinvCut"), qinvALp);
        }
        // QA3/SplitTrackQA: q_inv Before and After cut
        if (!failsQinvCut(qinvALp)) {
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LP_After_qinvCut"), qinvALp);
        }
        // QA3/RawPairDensity: AntiLambda-proton → same LP histogram
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          registryCorrelationQA.fill(
            HIST("QA3/RawPairDensity/h2f_dEta_dPhi_LP"),
            v0.eta() - primProton.eta(),
            RecoDecay::constrainAngle((v0.phi() - primProton.phi()), -PIHalf));
        }
        const int idxAL = unrolledIndex(v0.eta(), v0.phi());
        const int idxPrimP = unrolledIndex(primProton.eta(), primProton.phi());
        const float yAL = v0.rapidity(2);
        const float yP2 = protonRapidity(primProton);
        const int idxYAL = unrolledIndexY(yAL, v0.phi());
        const int idxYP2 = unrolledIndexY(yP2, primProton.phi());
        // Physics pair correlations (rho2 / pT / pair-count): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (failsQinvCut(qinvALp)) {
          continue;
        }
        registryCorrelationQA.fill(HIST("QA3/Kstar/h1f_kstar_ALp"),
                                   computeKstar(v0.px(), v0.py(), v0.pz(), kMassLambda,
                                                primProton.px(), primProton.py(), primProton.pz(), kMassProton));
        if (idxAL >= 0 && idxPrimP >= 0) {
          if (bin >= 0) {
            hPtPrimProton_AntiLambda[bin]->Fill(primProton.pt());
          }
          if (cFillEtaSpace.value) {
            hRho2_ALp_pT015toMaxDefined->Fill(idxAL, idxPrimP);
          }
        }
        if (idxYAL >= 0 && idxYP2 >= 0) {
          hRho2_ALp_y_pT015toMaxDefined->Fill(idxYAL, idxYP2);
        }
        if (cs) {
          fillIdx(cs->ALp, idxAL, idxPrimP);
          fillIdx(cs->ALp_y, idxYAL, idxYP2);
        }
        ++nPairs_AntiLambda_PrimProton;
      }
      for (auto const& primAntiProton : selectedPrimAntiProtons) {
        if (primAntiProton.globalIndex() == v0.posTrackId() || primAntiProton.globalIndex() == v0.negTrackId()) {
          continue;
        }
        // QA3/SplitTrackQA: q_inv for every accepted AntiLambda-antiproton pair
        float qinvALAp = 0.0f;
        {
          const float alpmag2 = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
          const float paprmag2 = std::sqrt(primAntiProton.px() * primAntiProton.px() + primAntiProton.py() * primAntiProton.py() + primAntiProton.pz() * primAntiProton.pz());
          const float EalL2 = std::sqrt(alpmag2 * alpmag2 + kMassLambda * kMassLambda);
          const float Epapr2 = std::sqrt(paprmag2 * paprmag2 + kMassProton * kMassProton);
          qinvALAp = computeQinv(v0.px(), v0.py(), v0.pz(), EalL2, primAntiProton.px(), primAntiProton.py(), primAntiProton.pz(), Epapr2);
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LP_Before_qinvCut"), qinvALAp);
        }
        // QA3/SplitTrackQA: q_inv Before and After cut
        if (!failsQinvCut(qinvALAp)) {
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LP_After_qinvCut"), qinvALAp);
        }
        // QA3/SplitTrackQA: close-pair block for daughter-overlap diagnosed split tracks
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          const float deta_alap = v0.eta() - primAntiProton.eta();
          const float dphi_alap = RecoDecay::constrainAngle((v0.phi() - primAntiProton.phi()), -PIHalf);
          registryCorrelationQA.fill(
            HIST("QA3/RawPairDensity/h2f_dEta_dPhi_LP"),
            deta_alap, dphi_alap);
          constexpr float cCloseDEta = 0.003f;
          constexpr float cCloseDPhi = 0.005f;
          if (std::abs(deta_alap) < cCloseDEta && std::abs(dphi_alap) < cCloseDPhi) {
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_sharedClsFraction_t1"), v0NegDau.tpcFractionSharedCls());
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_sharedClsFraction_t2"), primAntiProton.tpcFractionSharedCls());
          }
        }
        const int idxAL2 = unrolledIndex(v0.eta(), v0.phi());
        const int idxPrimPb = unrolledIndex(primAntiProton.eta(), primAntiProton.phi());
        const float yAL2 = v0.rapidity(2);
        const float yAp2 = protonRapidity(primAntiProton);
        const int idxYAL2 = unrolledIndexY(yAL2, v0.phi());
        const int idxYAp2 = unrolledIndexY(yAp2, primAntiProton.phi());
        // Physics pair correlations (rho2 / pT / pair-count): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (failsQinvCut(qinvALAp)) {
          continue;
        }
        registryCorrelationQA.fill(HIST("QA3/Kstar/h1f_kstar_ALAp"),
                                   computeKstar(v0.px(), v0.py(), v0.pz(), kMassLambda,
                                                primAntiProton.px(), primAntiProton.py(), primAntiProton.pz(), kMassProton));
        if (idxAL2 >= 0 && idxPrimPb >= 0) {
          if (bin >= 0) {
            hPtPrimAntiProton_AntiLambda[bin]->Fill(primAntiProton.pt());
          }
          if (cFillEtaSpace.value) {
            hRho2_ALAp_pT015toMaxDefined->Fill(idxAL2, idxPrimPb);
          }
        }
        if (idxYAL2 >= 0 && idxYAp2 >= 0) {
          hRho2_ALAp_y_pT015toMaxDefined->Fill(idxYAL2, idxYAp2);
        }
        if (cs) {
          fillIdx(cs->ALAp, idxAL2, idxPrimPb);
          fillIdx(cs->ALAp_y, idxYAL2, idxYAp2);
        }
        ++nPairs_AntiLambda_PrimAntiProton;
      }

      // ── Combo 3: AntiLambda daughter antiproton ↔ primary proton ─────
      // The AntiLambda neg-daughter (charge −1) CANNOT appear in the primary
      // proton list (charge +1), so same-track matching is impossible.
      // We fill only the inclusive Δη vs Δφ histogram.
      {
        const float etaDauAL = v0NegDau.eta();
        const float phiDauAL = v0NegDau.phi();
        for (auto const& primProton : selectedPrimProtons) {
          const float dEta3 = primProton.eta() - etaDauAL;
          float dPhi3 = primProton.phi() - phiDauAL;
          dPhi3 = RecoDecay::constrainAngle(dPhi3, -o2::constants::math::PI);
          registryCorrelationQA.fill(
            HIST("PrimaryProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi"),
            dEta3, dPhi3);
          registryCorrelationQA.fill(
            HIST("PrimaryProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_display"),
            dEta3, displayDeltaPhi(dPhi3));
        }
      }
      // ─────────────────────────────────────────────────────────────────

      // ── Combo 4: AntiLambda daughter antiproton ↔ primary antiproton ─
      // AntiLambda neg-daughter = v0NegDau (charge −1, identified as
      // antiproton). Primary antiproton = selectedPrimAntiProtons (charge −1).
      // Fills inclusive h2f_dEta_dPhi for all pairs and h2f_dEta_dPhi_nonSameTrack
      // for every pair in this cleaned sample (selectedPrimAntiProtons is already
      // veto-cleaned, so the primary antiproton can never be the AntiLambda
      // negative daughter; the nonSameTrack histogram will equal the inclusive one).
      {
        const float etaDauAL4 = v0NegDau.eta();
        const float phiDauAL4 = v0NegDau.phi();
        for (auto const& primAntiProton : selectedPrimAntiProtons) {
          const float dEta4 = primAntiProton.eta() - etaDauAL4;
          float dPhi4 = primAntiProton.phi() - phiDauAL4;
          dPhi4 = RecoDecay::constrainAngle(dPhi4, -o2::constants::math::PI);

          // Always fill the inclusive histogram.
          registryCorrelationQA.fill(
            HIST("PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi"),
            dEta4, dPhi4);
          registryCorrelationQA.fill(
            HIST("PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_display"),
            dEta4, displayDeltaPhi(dPhi4));

          // Fill non-same-track histogram when primary antiproton ≠ AntiLambda neg-daughter.
          if (primAntiProton.globalIndex() != v0.negTrackId()) {
            registryCorrelationQA.fill(
              HIST("PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_nonSameTrack"),
              dEta4, dPhi4);
            registryCorrelationQA.fill(
              HIST("PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_nonSameTrack_display"),
              dEta4, displayDeltaPhi(dPhi4));
          }
        }
      }
      // ─────────────────────────────────────────────────────────────────
    }
    // ====== End of Step 5 pairing loops ===================================

    // ── Lambda-Lambda rho2 fills ───────────────────────────────────────
    // LL
    for (auto const& v1 : selectedLambdas) {
      const float yv1L = v1.rapidity(1);
      for (auto const& v2 : selectedLambdas) {
        if (v1.index() == v2.index()) {
          continue;
        }
        if (v1.posTrackId() == v2.posTrackId()) {
          continue;
        }
        if (v1.posTrackId() == v2.negTrackId()) {
          continue;
        }
        if (v1.negTrackId() == v2.posTrackId()) {
          continue;
        }
        if (v1.negTrackId() == v2.negTrackId()) {
          continue;
        }
        // QA3/RawPairDensity: Lambda-Lambda → h2f_dEta_dPhi_LL
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          registryCorrelationQA.fill(
            HIST("QA3/RawPairDensity/h2f_dEta_dPhi_LL"),
            v1.eta() - v2.eta(),
            RecoDecay::constrainAngle((v1.phi() - v2.phi()), -PIHalf));
        }
        // QA3/SplitTrackQA: q_inv for every accepted Lambda-Lambda pair
        {
          const float v1pmag = std::sqrt(v1.px() * v1.px() + v1.py() * v1.py() + v1.pz() * v1.pz());
          const float v2pmag = std::sqrt(v2.px() * v2.px() + v2.py() * v2.py() + v2.pz() * v2.pz());
          const float Ev1 = std::sqrt(v1pmag * v1pmag + kMassLambda * kMassLambda);
          const float Ev2 = std::sqrt(v2pmag * v2pmag + kMassLambda * kMassLambda);
          const float qinvLL = computeQinv(v1.px(), v1.py(), v1.pz(), Ev1, v2.px(), v2.py(), v2.pz(), Ev2);
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_Before_qinvCut"), qinvLL);
          if (failsQinvCut(qinvLL)) {
            continue;
          }
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut"), qinvLL);
        }
        const int idx1 = unrolledIndex(v1.eta(), v1.phi());
        const int idx2 = unrolledIndex(v2.eta(), v2.phi());
        const int idxY1 = unrolledIndexY(yv1L, v1.phi());
        const int idxY2 = unrolledIndexY(v2.rapidity(1), v2.phi());
        // Physics pair correlations (rho2): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (cFillEtaSpace.value && idx1 >= 0 && idx2 >= 0) {
          hRho2_LL_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_LL_y_pT015toMaxDefined->Fill(idxY1, idxY2);
        }
        if (cs) {
          fillIdx(cs->LL, idx1, idx2);
          fillIdx(cs->LL_y, idxY1, idxY2);
        }
      }
    }
    // LAL
    for (auto const& v1 : selectedLambdas) {
      const float yv1LAL = v1.rapidity(1);
      for (auto const& v2 : selectedAntiLambdas) {
        if (v1.posTrackId() == v2.posTrackId()) {
          continue;
        }
        if (v1.posTrackId() == v2.negTrackId()) {
          continue;
        }
        if (v1.negTrackId() == v2.posTrackId()) {
          continue;
        }
        if (v1.negTrackId() == v2.negTrackId()) {
          continue;
        }
        // QA3/RawPairDensity: Lambda-AntiLambda → same h2f_dEta_dPhi_LL
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          registryCorrelationQA.fill(
            HIST("QA3/RawPairDensity/h2f_dEta_dPhi_LL"),
            v1.eta() - v2.eta(),
            RecoDecay::constrainAngle((v1.phi() - v2.phi()), -PIHalf));
        }
        // QA3/SplitTrackQA: q_inv for every accepted Lambda-AntiLambda pair
        {
          const float v1pmagL = std::sqrt(v1.px() * v1.px() + v1.py() * v1.py() + v1.pz() * v1.pz());
          const float v2pmagL = std::sqrt(v2.px() * v2.px() + v2.py() * v2.py() + v2.pz() * v2.pz());
          const float Ev1L = std::sqrt(v1pmagL * v1pmagL + kMassLambda * kMassLambda);
          const float Ev2L = std::sqrt(v2pmagL * v2pmagL + kMassLambda * kMassLambda);
          const float qinvLAL = computeQinv(v1.px(), v1.py(), v1.pz(), Ev1L, v2.px(), v2.py(), v2.pz(), Ev2L);
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_Before_qinvCut"), qinvLAL);
          if (failsQinvCut(qinvLAL)) {
            continue;
          }
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut"), qinvLAL);
        }
        const int idx1 = unrolledIndex(v1.eta(), v1.phi());
        const int idx2 = unrolledIndex(v2.eta(), v2.phi());
        const int idxY1 = unrolledIndexY(yv1LAL, v1.phi());
        const int idxY2 = unrolledIndexY(v2.rapidity(2), v2.phi());
        // Physics pair correlations (rho2): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (cFillEtaSpace.value && idx1 >= 0 && idx2 >= 0) {
          hRho2_LAL_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_LAL_y_pT015toMaxDefined->Fill(idxY1, idxY2);
        }
        if (cs) {
          fillIdx(cs->LAL, idx1, idx2);
          fillIdx(cs->LAL_y, idxY1, idxY2);
        }
      }
    }
    // ALL
    for (auto const& v1 : selectedAntiLambdas) {
      const float yv1ALL = v1.rapidity(2);
      for (auto const& v2 : selectedLambdas) {
        if (v1.posTrackId() == v2.posTrackId()) {
          continue;
        }
        if (v1.posTrackId() == v2.negTrackId()) {
          continue;
        }
        if (v1.negTrackId() == v2.posTrackId()) {
          continue;
        }
        if (v1.negTrackId() == v2.negTrackId()) {
          continue;
        }
        // QA3/RawPairDensity: AntiLambda-Lambda → same h2f_dEta_dPhi_LL
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          registryCorrelationQA.fill(
            HIST("QA3/RawPairDensity/h2f_dEta_dPhi_LL"),
            v1.eta() - v2.eta(),
            RecoDecay::constrainAngle((v1.phi() - v2.phi()), -PIHalf));
        }
        // QA3/SplitTrackQA: q_inv for every accepted AntiLambda-Lambda pair
        {
          const float v1pmagAL = std::sqrt(v1.px() * v1.px() + v1.py() * v1.py() + v1.pz() * v1.pz());
          const float v2pmagAL = std::sqrt(v2.px() * v2.px() + v2.py() * v2.py() + v2.pz() * v2.pz());
          const float Ev1AL = std::sqrt(v1pmagAL * v1pmagAL + kMassLambda * kMassLambda);
          const float Ev2AL = std::sqrt(v2pmagAL * v2pmagAL + kMassLambda * kMassLambda);
          const float qinvALL = computeQinv(v1.px(), v1.py(), v1.pz(), Ev1AL, v2.px(), v2.py(), v2.pz(), Ev2AL);
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_Before_qinvCut"), qinvALL);
          if (failsQinvCut(qinvALL)) {
            continue;
          }
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut"), qinvALL);
        }
        const int idx1 = unrolledIndex(v1.eta(), v1.phi());
        const int idx2 = unrolledIndex(v2.eta(), v2.phi());
        const int idxY1 = unrolledIndexY(yv1ALL, v1.phi());
        const int idxY2 = unrolledIndexY(v2.rapidity(1), v2.phi());
        // Physics pair correlations (rho2): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (cFillEtaSpace.value && idx1 >= 0 && idx2 >= 0) {
          hRho2_ALL_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_ALL_y_pT015toMaxDefined->Fill(idxY1, idxY2);
        }
        if (cs) {
          fillIdx(cs->ALL, idx1, idx2);
          fillIdx(cs->ALL_y, idxY1, idxY2);
        }
      }
    }
    // ALAL
    for (auto const& v1 : selectedAntiLambdas) {
      const float yv1ALAL = v1.rapidity(2);
      for (auto const& v2 : selectedAntiLambdas) {
        if (v1.index() == v2.index()) {
          continue;
        }
        if (v1.posTrackId() == v2.posTrackId()) {
          continue;
        }
        if (v1.posTrackId() == v2.negTrackId()) {
          continue;
        }
        if (v1.negTrackId() == v2.posTrackId()) {
          continue;
        }
        if (v1.negTrackId() == v2.negTrackId()) {
          continue;
        }
        // QA3/RawPairDensity: AntiLambda-AntiLambda → same h2f_dEta_dPhi_LL
        {
          constexpr float PIHalf = o2::constants::math::PI / 2.0f;
          registryCorrelationQA.fill(
            HIST("QA3/RawPairDensity/h2f_dEta_dPhi_LL"),
            v1.eta() - v2.eta(),
            RecoDecay::constrainAngle((v1.phi() - v2.phi()), -PIHalf));
        }
        // QA3/SplitTrackQA: q_inv for every accepted AntiLambda-AntiLambda pair
        {
          const float v1pmagAL2 = std::sqrt(v1.px() * v1.px() + v1.py() * v1.py() + v1.pz() * v1.pz());
          const float v2pmagAL2 = std::sqrt(v2.px() * v2.px() + v2.py() * v2.py() + v2.pz() * v2.pz());
          const float Ev1AL2 = std::sqrt(v1pmagAL2 * v1pmagAL2 + kMassLambda * kMassLambda);
          const float Ev2AL2 = std::sqrt(v2pmagAL2 * v2pmagAL2 + kMassLambda * kMassLambda);
          const float qinvALAL = computeQinv(v1.px(), v1.py(), v1.pz(), Ev1AL2, v2.px(), v2.py(), v2.pz(), Ev2AL2);
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_Before_qinvCut"), qinvALAL);
          if (failsQinvCut(qinvALAL)) {
            continue;
          }
          registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut"), qinvALAL);
        }
        const int idx1 = unrolledIndex(v1.eta(), v1.phi());
        const int idx2 = unrolledIndex(v2.eta(), v2.phi());
        const int idxY1 = unrolledIndexY(yv1ALAL, v1.phi());
        const int idxY2 = unrolledIndexY(v2.rapidity(2), v2.phi());
        // Physics pair correlations (rho2): q_inv cut applied (pairs with q_inv <= kQinvCutLP are rejected above)
        if (cFillEtaSpace.value && idx1 >= 0 && idx2 >= 0) {
          hRho2_ALAL_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_ALAL_y_pT015toMaxDefined->Fill(idxY1, idxY2);
        }
        if (cs) {
          fillIdx(cs->ALAL, idx1, idx2);
          fillIdx(cs->ALAL_y, idxY1, idxY2);
        }
      }
    }
    // ─────────────────────────────────────────────────────────────────────

    registryOther.fill(HIST("hNSelLambda"), static_cast<float>(nSelLambda));
    registryOther.fill(HIST("hNSelAntiLambda"), static_cast<float>(nSelAntiLambda));
    registryOther.fill(HIST("hNPairs_Lambda_PrimProton"), static_cast<float>(nPairs_Lambda_PrimProton));
    registryOther.fill(HIST("hNPairs_Lambda_PrimAntiProton"), static_cast<float>(nPairs_Lambda_PrimAntiProton));
    registryOther.fill(HIST("hNPairs_AntiLambda_PrimProton"), static_cast<float>(nPairs_AntiLambda_PrimProton));
    registryOther.fill(HIST("hNPairs_AntiLambda_PrimAntiProton"), static_cast<float>(nPairs_AntiLambda_PrimAntiProton));
  }
  // ── Run-level summary printed once at the end of the analysis job ─────
  // Called automatically by the O2 framework after all events are processed.
  // Prints totals accumulated in runTotalProtonBeforeVeto / runTotalProtonRemovedVeto.
  void endOfStream(EndOfStreamContext const&)
  {
    const int64_t totalBefore = runTotalProtonBeforeVeto;
    const int64_t totalRemoved = runTotalProtonRemovedVeto;
    const int64_t totalAfter = totalBefore - totalRemoved;
    const double fracRemoved = (totalBefore > 0)
                                 ? static_cast<double>(totalRemoved) / static_cast<double>(totalBefore)
                                 : 0.0;

    LOG(info) << "";
    LOG(info) << "================ Proton Veto QA ================";
    LOG(info) << "Total proton candidates before veto : " << totalBefore;
    LOG(info) << "Total proton candidates after veto  : " << totalAfter;
    LOG(info) << "Total removed by veto               : " << totalRemoved;
    LOG(info) << Form("Fraction removed                    : %.6f", fracRemoved);
    const double fracRemovedAp = (runTotalAntiPBeforeVeto > 0)
                                   ? static_cast<double>(runTotalAntiPRemovedVeto) / static_cast<double>(runTotalAntiPBeforeVeto)
                                   : 0.0;
    LOG(info) << "Total antiproton candidates before veto: " << runTotalAntiPBeforeVeto;
    LOG(info) << "Total antiproton removed by veto       : " << runTotalAntiPRemovedVeto;
    LOG(info) << Form("Antiproton fraction removed         : %.6f", fracRemovedAp);
    LOG(info) << "===============================================";
    LOG(info) << "";
  }
  // ─────────────────────────────────────────────────────────────────────
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LambdaProtonBalanceFunction>(cfgc)};
}
