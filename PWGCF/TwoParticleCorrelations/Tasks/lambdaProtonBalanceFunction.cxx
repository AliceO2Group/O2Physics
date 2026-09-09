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
/// \brief Single- and two-particle density correlations (rho1, rho2) and R2 correlations for Lambda-proton pairs to measure the balance function.
/// \author Anoop Poruthiyil <anoop.poruthiyil@cern.ch>

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <memory>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA,
                           aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPr,
                           aod::pidTOFbeta, aod::TrackSelection>;

struct LambdaProtonBalanceFunction {
  HistogramRegistry registryLambda{"Lambda_invMass", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryLambdaExtended{"Lambda_invMassExtended", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryRho{"rho1ANDrho2_LP", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryOther{"Other_Hists", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryPid{"PID", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQaDetector{"QA_Detec", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryProtonPidCheck{"ProtonCounts_byPID_Check", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryCorrelationQA{"CorrelationQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // ── Dedicated registry for V0-daughter veto importance QA ─────────────
  // Folder: ProtonVetoQA/
  // Contains per-event distributions of N_before, N_after, N_removed, fraction.
  HistogramRegistry registryProtonVetoQA{"ProtonVetoQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // pT bins 0.6 → 6.0 in steps of 0.2 (used for per-bin rho2 / mass histograms)
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
  //   Rapidity bins: [-1.0,-0.9), [-0.9,-0.8), [-0.8,-0.7), [-0.7,-0.6]
  //   pT bins:       [0.0,0.1), [0.1,0.2), ..., [2.4,2.5]
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
  //           hRho2_<channel>_pT015toMaxDefined = merged (full pT range)
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
  // ─────────────────────────────────────────────────────────────────────

  // ── Configurables ─────────────────────────────────────────────────────
  Configurable<int> cPPandV0NBins{"cPPandV0NBins", 100, "N bins in all histos"};
  Configurable<float> cPPandV0ZVertexCut{"cPPandV0ZVertexCut", 10.0f, "Accepted z-vertex range (cm)"};

  // Track quality cuts
  Configurable<int> pProtonTPCMinRows{"pProtonTPCMinRows", 70, "Minimum TPC crossed rows"};
  Configurable<float> pProtonTPCMinRowsOverFindable{"pProtonTPCMinRowsOverFindable", 0.8f, "Min TPC crossed rows / findable clusters"};
  Configurable<float> pProtonTPCMaxChi2PerCluster{"pProtonTPCMaxChi2PerCluster", 4.0f, "Max TPC chi2 / Ncls"};
  Configurable<int> pProtonITSMinClusters{"pProtonITSMinClusters", 5, "Minimum ITS clusters"};
  Configurable<float> pProtonITSMaxChi2PerCluster{"pProtonITSMaxChi2PerCluster", 36.0f, "Max ITS chi2 / Ncls"};
  Configurable<float> cPPMinTpcSignal{"cPPMinTpcSignal", 10.0f, "Min raw TPC signal (sanity cut)"};
  // DCA cuts
  Configurable<float> pProtonMaxDCAxy{"pProtonMaxDCAxy", 0.02f, "Max |DCAxy| (cm)"}; // changed 0.1-->0.05-->0.02 (tightened for cleanest primary-proton sample)
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
  Configurable<float> pProtonMinP{"pProtonMinP", 0.5f, "Min p_{TPC} for primary protons"};
  Configurable<float> pProtonMaxP{"pProtonMaxP", 3.5999999f, "Max p for primary protons"};
  Configurable<float> pProtonMaxEta{"pProtonMaxEta", 0.5f, "Max |eta| for primary protons"};
  Configurable<float> pProtonMaxY{"pProtonMaxY", 0.5f, "Max |y| for primary protons"};

  // Lambda Kinematic Cuts (from strange9)
  Configurable<float> lambdaV0MinPt{"lambdaV0MinPt", 0.60f, "Min pT for Lambdas"};
  Configurable<float> lambdaV0MaxPt{"lambdaV0MaxPt", 3.60f, "Max pT for Lambdas"};
  Configurable<float> lambdaV0MaxY{"lambdaV0MaxY", 0.5f, "Max |y| for Lambdas"};
  Configurable<float> lambdaV0MaxEta{"lambdaV0MaxEta", 0.5f, "Max |eta| for Lambdas"};

  // V0 Mass and Topology Cuts (from strange9)
  Configurable<float> lambdaV0MaxCTau{"lambdaV0MaxCTau", 30.0f, "Max V0 ctau (cm)"};
  Configurable<float> lambdaV0KShortRejectMassWindow{"lambdaV0KShortRejectMassWindow", 0.01f, "K0s rejection window"};
  Configurable<float> lambdaV0MassWindow{"lambdaV0MassWindow", 0.007f, "Lambda mass window (|m - mPDG| < this)"};

  // V0 Type and Daughter Quality Cuts (from strange9)
  Configurable<int> cV0TypeSelection{"cV0TypeSelection", 1, "V0 Type Selection"};
  Configurable<float> cV0DauMinPt{"cV0DauMinPt", 0.1f, "Daughter pT minimum"};
  Configurable<float> cV0DauMaxEta{"cV0DauMaxEta", 0.8f, "Daughter |eta| cut"};
  Configurable<int> cV0DauMinTpcCrossedRows{"cV0DauMinTpcCrossedRows", 70, "Daughter TPC min crossed rows"};
  // Rapidity-specific constants matched to PProtonMaxY's default value (0.8).
  // Note: Since these are compile-time constexpr, the axis will NOT auto-update
  // if PProtonMaxY is changed via configurable at runtime. The axis just needs to
  // stay >= the cut, not exactly equal.

  // When u change JSON, adjust this also (not needed is my new finding)

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

  void buildPtBins()
  {
    ptEdges.clear();
    static constexpr float kPtEdgeMax = 6.0001f;
    for (float pt = 0.6f; pt <= kPtEdgeMax; pt += 0.2f) {
      ptEdges.push_back(pt);
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

  // Femtoscopic q_inv cut (applied only to QA Before/After plots, not to physics rho2 histograms):
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
  // ─────────────────────────────────────────────────────────────────────

  // ── Two-region PID selection for (anti)protons ────────────────────────
  template <typename TTrack>
  bool passesPrimProtonPid(TTrack const& trk) const
  {
    const float pTPC = trk.tpcInnerParam(); // NOTE: pTPC, not trk.p() — matches reference exactly
    const float nsPr = trk.tpcNSigmaPr();

    if (pTPC < pProtonTPCTOFSwitchP.value) {
      if (std::abs(nsPr) >= pProtonTPCNsigma.value) {
        return false;
      }
      return true;
    } else {
      if (std::abs(nsPr) >= pProtonTPCNsigma.value) {
        return false;
      }
      if (!trk.hasTOF()) {
        return false;
      }
      if (std::abs(trk.tofNSigmaPr()) >= pProtonTOFNsigma.value) {
        return false;
      }
      return true;
    }
  }
  // ─────────────────────────────────────────────────────────────────────

  // ══════════════════════════════════════════════════════════════════════════
  // STAGED QA INFRASTRUCTURE — registryQaDetector ("QA_Detec")
  // ══════════════════════════════════════════════════════════════════════════

  /// Canonical analysis stages for (anti)proton QA fills.
  enum class ProtonQAStage : int {
    RawAfterTrackSel = 0, ///< truly RAW — before any track selection, DCA, or PID
    TpcPid = 1,           ///< after DCA + TPC PID cut passes
    TpcTofPid = 2,        ///< after DCA + TPC PID + TOF PID all pass
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
    const float tpcCR = static_cast<float>(trk.tpcNClsCrossedRows());
    const float tpcCRoF = trk.tpcCrossedRowsOverFindableCls();
    const float tpcChi2 = trk.tpcChi2NCl();
    const float itsN = static_cast<float>(trk.itsNCls());
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

  void init(InitContext const&)
  {
    buildPtBins();

    // ── Common axis definitions ──────────────────────────────────────────
    AxisSpec lambdaMassAxis = {200, 1.08f, 1.15f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {cPPandV0NBins, -15.f, 15.f, "vrtx_{Z} [cm]"};
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
    static constexpr int kRhoUnrolledBinsY = kRhoYBins * kRhoPhiBins;
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
    registryOther.add("hPtSelLambda", "#Lambda p_{T} yield", {HistType::kTH1F, {ptAxis}});
    registryOther.add("hPtSelAntiLambda", "#bar{#Lambda} p_{T} yield", {HistType::kTH1F, {ptAxis}});
    registryOther.add("hPtSelectedPrimProton", "Selected primary proton p_{T} yield", {HistType::kTH1F, {ptAxis}});
    registryOther.add("hPtSelectedPrimAntiProton", "Selected primary antiproton p_{T} yield", {HistType::kTH1F, {ptAxis}});
    registryOther.add("hNPairs_Lambda_PrimProton", "#Lambda-p pairs per event", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPairs_Lambda_PrimAntiProton", "#Lambda-#bar{p} pairs per event", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPairs_AntiLambda_PrimProton", "#bar{#Lambda}-p pairs per event", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});
    registryOther.add("hNPairs_AntiLambda_PrimAntiProton", "#bar{#Lambda}-#bar{p} pairs per event", {HistType::kTH1F, {{200, 0.0f, 200.0f}}});

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
    AxisSpec axisCutFlow = {5, 0.5f, 5.5f, "Cut stage"};
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
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_pt", species, stage), "p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {axisPtProton}});
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_eta", species, stage), "#eta;#eta;Counts", {HistType::kTH1F, {qaEtaAxis}});
      registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_phi", species, stage), "#phi;#phi;Counts", {HistType::kTH1F, {qaPhiAxis}});
      if (std::string(stage) == "02_Final") {
        registryQaDetector.add(Form("%s/01_Kinematics/%s/h1f_rapidity", species, stage), "y;y;Counts", {HistType::kTH1F, {rapidityAxis}});
      }

      // Topology
      registryQaDetector.add(Form("%s/02_Topology/%s/h1f_cospa", species, stage), "Cos(PA);Cos(PA);Counts", {HistType::kTH1F, {{100, 0.9f, 1.0f, "Cos(PA)"}}});
      registryQaDetector.add(Form("%s/02_Topology/%s/h1f_cospa_zoom", species, stage), "Cos(PA);Cos(PA);Counts", {HistType::kTH1F, {axisCosPA_zoom}});
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

    registryQaDetector.add("Selection/TOF_Matching/hTOFMatchedFractionVsPt", "TOF matched fraction vs p_{T};p_{T};Fraction", {HistType::kTProfile, {axisPtProton}});
    // ── Item 4: 2D TOF-matching efficiency maps (eta × phi) ──────────────────
    registryQaDetector.add("Selection/TOF_Matching/h2f_TOFMatchedFractionVsEtaPhi",
                           "Proton TOF matched fraction vs (#eta, #varphi);#eta;#varphi;Fraction",
                           {HistType::kTProfile2D, {qaEtaAxis, qaPhiAxis}});
    // ─────────────────────────────────────────────────────────────────────────
    registryQaDetector.add("Selection/TOF_Matching/hCutFlow", "Proton Track Cutflow;Stage;Counts", {HistType::kTH1F, {axisCutFlow}});

    // ══════════════════════════════════════════════════════════════════════════
    // EXISTING QA_Detec REGISTRATION (Preserved for backward compatibility)
    // ══════════════════════════════════════════════════════════════════════════

    // -- Proton (existing) --r INEL > 0
    registryOther.add("hEventCutflow", "Event Cutflow;Step;Counts",
                      {HistType::kTH1F, {{4, 0.5f, 4.5f}}});

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
    // Bin 5 = kIsVertexTRDmatched, Bin 6 = any pileup-rejection bit fired.
    registryOther.add("hPileupFlags",
                      "Event-level pile-up rejection flags (Bin 6: All 3 passed (AND));Flag;Events passing",
                      {HistType::kTH1F, {{6, 0.5f, 6.5f, "Pileup flag"}}});
    registryOther.add("h2f_NTracks_vs_Cent",
                      "N_{tracks}(|#eta|<0.8) vs FT0M cent;FT0M (%);N_{tracks}(|#eta|<0.8)",
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
    registryRho.add("h2_rho1_Proton", "rho1 proton", {HistType::kTH2F, {etaAxis, phiAxis}});
    registryRho.add("h2_rho1_AntiProton", "rho1 pbar", {HistType::kTH2F, {etaAxis, phiAxis}});
    registryRho.add("h2_rho1_Lambda", "rho1 lambda", {HistType::kTH2F, {etaAxis, phiAxis}});
    registryRho.add("h2_rho1_AntiLambda", "rho1 lambdabar", {HistType::kTH2F, {etaAxis, phiAxis}});
    // rapidity-based rho1
    registryRho.add("h2_rho1_Proton_y", "rho1 proton (y)", {HistType::kTH2F, {yAxis, phiAxis}});
    registryRho.add("h2_rho1_AntiProton_y", "rho1 pbar (y)", {HistType::kTH2F, {yAxis, phiAxis}});
    registryRho.add("h2_rho1_Lambda_y", "rho1 lambda (y)", {HistType::kTH2F, {yAxis, phiAxis}});
    registryRho.add("h2_rho1_AntiLambda_y", "rho1 lambdabar (y)", {HistType::kTH2F, {yAxis, phiAxis}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Autocorrelation / feed-down diagnostic histograms ────────────────
    // Purpose: check whether Lambda daughter protons enter the primary-proton
    //          sample, causing autocorrelation bias in the Lambda-p rho2.
    // Folder:  CorrelationQA/PrimaryProton_vs_LambdaDaughter/
    //
    // h2f_dEta_dPhi         : Δη vs Δφ for every (selected Lambda, selected
    //                         primary proton) pair.  No veto on shared tracks.
    //                         A peak at (0,0) indicates track sharing / feed-down.
    // h2f_dEta_dPhi_nonSameTrack : same, restricted to pairs where the primary
    //                         proton track is NOT the Lambda positive daughter.
    AxisSpec axisCorDEta = {64, -1.6f, 1.6f};
    AxisSpec axisCorDPhi = {72, -o2::constants::math::PI, o2::constants::math::PI, "#Delta#varphi"};
    AxisSpec axisCorDPhiDisplay = {72, -o2::constants::math::PI / 2.0f, 3.0f * o2::constants::math::PI / 2.0f, "#Delta#varphi"};
    registryCorrelationQA.add(
      "PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi",
      "Primary Proton vs #Lambda Daughter Proton;#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_display",
      "Primary Proton vs #Lambda Daughter Proton (display shifted #Delta#varphi);#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    registryCorrelationQA.add(
      "PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_nonSameTrack",
      "Primary Proton vs #Lambda Daughter (non-same-track pairs only);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryProton_vs_LambdaDaughter/h2f_dEta_dPhi_nonSameTrack_display",
      "Primary Proton vs #Lambda Daughter (non-same-track pairs only, display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Combo 2: Lambda daughter proton ↔ primary antiproton ─────────────
    // Purpose: check whether the Lambda pos-daughter proton (charge +1) could
    //          be confused with a primary antiproton (charge −1).
    // Same-track matching is impossible here (opposite charges).
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_LambdaDaughter/h2f_dEta_dPhi",
      "Primary AntiProton vs #Lambda Daughter Proton;#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_LambdaDaughter/h2f_dEta_dPhi_display",
      "Primary AntiProton vs #Lambda Daughter Proton (display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Combo 3: AntiLambda daughter antiproton ↔ primary proton ─────────
    // Purpose: check whether the AntiLambda neg-daughter antiproton (charge −1)
    //          could leak into the primary proton sample (charge +1).
    // Same-track matching is impossible here (opposite charges).
    registryCorrelationQA.add(
      "PrimaryProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi",
      "Primary Proton vs #bar{#Lambda} Daughter AntiProton;#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_display",
      "Primary Proton vs #bar{#Lambda} Daughter AntiProton (display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
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
      "Primary AntiProton vs #bar{#Lambda} Daughter AntiProton;#Delta#eta = #eta_{prim} - #eta_{dau};#Delta#varphi = #varphi_{prim} - #varphi_{dau}",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_display",
      "Primary AntiProton vs #bar{#Lambda} Daughter AntiProton (display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhiDisplay}});
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_nonSameTrack",
      "Primary AntiProton vs #bar{#Lambda} Daughter (non-same-track pairs only);#Delta#eta;#Delta#varphi",
      {HistType::kTH2F, {axisCorDEta, axisCorDPhi}});
    registryCorrelationQA.add(
      "PrimaryAntiProton_vs_AntiLambdaDaughter/h2f_dEta_dPhi_nonSameTrack_display",
      "Primary AntiProton vs #bar{#Lambda} Daughter (non-same-track pairs only, display shifted #Delta#varphi);#Delta#eta;#Delta#varphi",
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
      {HistType::kTH1F, {{3, 0.5f, 3.5f}}});
    // Bin labels are set after histogram registration
    {
      auto h = registryCorrelationQA.get<TH1>(HIST("QA3/PairCleaning/h1f_pairCleaning"));
      h->GetXaxis()->SetBinLabel(1, "All Candidate Pairs");
      h->GetXaxis()->SetBinLabel(2, "Rejected Same Track");
      h->GetXaxis()->SetBinLabel(3, "Accepted");
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
      "Split-track QA: q_{inv} for #Lambda-p pairs, before q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_LL_Before_qinvCut",
      "Split-track QA: q_{inv} for #Lambda-#Lambda pairs, before q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    // "After_qinvCut": filled only for pairs surviving the q_inv cut
    // (q_inv > kQinvCutLP), applied uniformly to PP, LP, and LL pairs.
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_PP_After_qinvCut",
      "Split-track QA: q_{inv} for p-p pairs, after q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_LP_After_qinvCut",
      "Split-track QA: q_{inv} for #Lambda-p pairs, after q_{inv} cut (q_{inv} > 0.01 GeV/c);q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut",
      "Split-track QA: q_{inv} for #Lambda-#Lambda pairs, after q_{inv} cut;q_{inv} (GeV/c);Counts",
      {HistType::kTH1F, {{600, 0.0f, 6.0f, "q_{inv} (GeV/c)"}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_sharedClsFraction_t1",
      "Split-track QA: TPC shared-cluster fraction (track 1 in close p-p pairs);Fraction;Counts",
      {HistType::kTH1F, {{50, 0.0f, 1.0f}}});
    registryCorrelationQA.add(
      "QA3/SplitTrackQA/h1f_sharedClsFraction_t2",
      "Split-track QA: TPC shared-cluster fraction (track 2 in close p-p pairs);Fraction;Counts",
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
                                      "Proton candidates: TPC mandatory, TOF optional if present;#it{p}_{T} (GeV/#it{c});Counts",
                                      {HistType::kTH1F, {pidCheckPtAxis}});
    // ─────────────────────────────────────────────────────────────────────

    // ── Merged / pT015toMaxDefined histograms ───────────────────────────────────
    hRho2_Lp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_Lp_pT015toMaxDefined",
                           "Lambda-p correlation (integrated pT range)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_LAp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LAp_pT015toMaxDefined",
                           "Lambda-pbar correlation (integrated pT range)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ALp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALp_pT015toMaxDefined",
                           "AntiLambda-p correlation (integrated pT range)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ALAp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALAp_pT015toMaxDefined",
                           "AntiLambda-pbar correlation (integrated pT range)",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_pp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_pp_pT015toMaxDefined",
                           "pp correlation",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_pAp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_pAp_pT015toMaxDefined",
                           "p-pbar correlation",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_App_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_App_pT015toMaxDefined",
                           "pbar-p correlation",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ApAp_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ApAp_pT015toMaxDefined",
                           "pbar-pbar correlation",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_LL_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LL_pT015toMaxDefined",
                           "Lambda-Lambda correlation",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_LAL_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LAL_pT015toMaxDefined",
                           "Lambda-AntiLambda correlation",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ALL_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALL_pT015toMaxDefined",
                           "AntiLambda-Lambda correlation",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});
    hRho2_ALAL_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALAL_pT015toMaxDefined",
                           "AntiLambda-AntiLambda correlation",
                           {HistType::kTH2F, {unrolledAxis, unrolledAxis}});

    // ── Rapidity merged rho2 histograms ─────────────────────────────────────
    // NOTE: these must use unrolledAxisY (kRhoYBins*kRhoPhiBins bins),
    // NOT the eta-space unrolledAxis (kRhoEtaBins*kRhoPhiBins bins).
    // GetUnrolledIndexY() below produces indices in [0, kRhoUnrolledBinsY),
    // so the booked axis must match that range or downstream R2/BF code that
    // expects y-space-sized histograms will see a dimension mismatch.
    hRho2_Lp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_Lp_y_pT015toMaxDefined", "Lambda-p (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_LAp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LAp_y_pT015toMaxDefined", "Lambda-pbar (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ALp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALp_y_pT015toMaxDefined", "AntiLambda-p (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ALAp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALAp_y_pT015toMaxDefined", "AntiLambda-pbar (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_pp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_pp_y_pT015toMaxDefined", "pp (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_pAp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_pAp_y_pT015toMaxDefined", "p-pbar (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_App_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_App_y_pT015toMaxDefined", "pbar-p (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ApAp_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ApAp_y_pT015toMaxDefined", "pbar-pbar (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_LL_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LL_y_pT015toMaxDefined", "Lambda-Lambda (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_LAL_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_LAL_y_pT015toMaxDefined", "Lambda-AntiLambda (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ALL_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALL_y_pT015toMaxDefined", "AntiLambda-Lambda (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    hRho2_ALAL_y_pT015toMaxDefined =
      registryRho.add<TH2>("h2_rho2_ALAL_y_pT015toMaxDefined", "AntiLambda-AntiLambda (y)", {HistType::kTH2F, {unrolledAxisY, unrolledAxisY}});
    // ─────────────────────────────────────────────────────────────────────

    hPtSelectedPrimProton_pT015toMaxDefined =
      registryOther.add<TH1>("hPtSelectedPrimProton_pT015toMaxDefined",
                             "Proton p_{T} yield (integrated pT range)",
                             {HistType::kTH1F, {ptAxis}});
    hPtSelectedPrimAntiProton_pT015toMaxDefined =
      registryOther.add<TH1>("hPtSelectedPrimAntiProton_pT015toMaxDefined",
                             "Antiproton p_{T} yield (integrated pT range)",
                             {HistType::kTH1F, {ptAxis}});
    hPtSelLambda_pT015toMaxDefined =
      registryOther.add<TH1>("hPtSelLambda_pT015toMaxDefined",
                             "#Lambda p_{T} yield (integrated pT range)",
                             {HistType::kTH1F, {ptAxis}});
    hPtSelAntiLambda_pT015toMaxDefined =
      registryOther.add<TH1>("hPtSelAntiLambda_pT015toMaxDefined",
                             "#bar{#Lambda} p_{T} yield (integrated pT range)",
                             {HistType::kTH1F, {ptAxis}});

    hMassLambdaMerged =
      registryLambda.add<TH1>("hMassLambdaNoMassCut_pT015toMaxDefined",
                              "Lambda Mass No Mass Cut (integrated pT range)",
                              {HistType::kTH1F, {lambdaMassAxis}});
    hMassAntiLambdaMerged =
      registryLambda.add<TH1>("hMassAntiLambdaNoMassCut_pT015toMaxDefined",
                              "AntiLambda Mass No Mass Cut (integrated pT range)",
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
          Form("Proton pT in #Lambda-p pairs %.1f < pT < %.1f", ptLo, ptHi),
          {HistType::kTH1F, {ptAxis}}));
      hPtPrimAntiProton_Lambda.push_back(
        registryOther.add<TH1>(
          Form("hPtPrimAntiProton_Lambda_ptBin%02zu", i),
          Form("Antiproton pT in #Lambda-#bar{p} pairs %.1f < pT < %.1f", ptLo, ptHi),
          {HistType::kTH1F, {ptAxis}}));
      hPtPrimProton_AntiLambda.push_back(
        registryOther.add<TH1>(
          Form("hPtPrimProton_AntiLambda_ptBin%02zu", i),
          Form("Proton pT in #bar{#Lambda}-p pairs %.1f < pT < %.1f", ptLo, ptHi),
          {HistType::kTH1F, {ptAxis}}));
      hPtPrimAntiProton_AntiLambda.push_back(
        registryOther.add<TH1>(
          Form("hPtPrimAntiProton_AntiLambda_ptBin%02zu", i),
          Form("Antiproton pT in #bar{#Lambda}-#bar{p} pairs %.1f < pT < %.1f", ptLo, ptHi),
          {HistType::kTH1F, {ptAxis}}));

      hMassLambdaPtBins.push_back(
        registryLambda.add<TH1>(
          Form("hMassLambdaNoMassCut_ptBin%02zu", i),
          Form("Lambda Mass No Mass Cut %.1f < pT < %.1f", ptLo, ptHi),
          {HistType::kTH1F, {lambdaMassAxis}}));
      hMassAntiLambdaPtBins.push_back(
        registryLambda.add<TH1>(
          Form("hMassAntiLambdaNoMassCut_ptBin%02zu", i),
          Form("AntiLambda Mass No Mass Cut %.1f < pT < %.1f", ptLo, ptHi),
          {HistType::kTH1F, {lambdaMassAxis}}));
    }
    // ─────────────────────────────────────────────────────────────────────

    // ── Extended y × pT invariant-mass histograms ─────────────────────────
    // Rapidity bins: 4 bins from -1.0 to -0.6 in steps of 0.1
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
            Form("Lambda Invariant Mass | %.1f <= y %s %.1f | %.1f <= pT %s %.1f GeV/c",
                 yLow, yUpBracket, yHigh, ptLow, ptUpBracket, ptHigh),
            {HistType::kTH1F, {lambdaMassAxis}});

        // AntiLambda
        hMassAntiLambdaExtended[iY][iPt] =
          registryLambdaExtended.add<TH1>(
            Form("AntiLambda/hMassAntiLambda_yBin%02d_ptBin%02d", iY, iPt),
            Form("AntiLambda Invariant Mass | %.1f <= y %s %.1f | %.1f <= pT %s %.1f GeV/c",
                 yLow, yUpBracket, yHigh, ptLow, ptUpBracket, ptHigh),
            {HistType::kTH1F, {lambdaMassAxis}});
      }
    }
    // ─────────────────────────────────────────────────────────────────────
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

    // INEL > 0: require at least one charged track in |eta| < 1.0
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

      // ── TRULY RAW QA Fills — before any V0 selection cut ─────────────
      // Placed here at the very top of the loop: before v0Type, daughter
      // kinematics, TPC PID, rapidity, pT, topology, DCA, CTau, K0S, or
      // mass-window cuts. Every V0 in the event is filled.
      // We double-fill under both the Lambda hypothesis (pos=proton, neg=pion)
      // and the AntiLambda hypothesis (pos=pion, neg=proton) since no PID
      // decision has been made yet. This gives the true combinatorial view.
      {
        const auto& rawPosDau = v0.posTrack_as<MyTracks>();
        const auto& rawNegDau = v0.negTrack_as<MyTracks>();
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

      const auto& v0PosDau = v0.posTrack_as<MyTracks>();
      const auto& v0NegDau = v0.negTrack_as<MyTracks>();

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
        // Accept only -1.0 <= y <= -0.6 and 0.0 <= pT <= 2.5 GeV/c.
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
        // Accept only -1.0 <= y <= -0.6 and 0.0 <= pT <= 2.5 GeV/c.
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

        registryRho.fill(HIST("h2_rho1_Lambda"), v0.eta(), v0.phi());
        registryRho.fill(HIST("h2_rho1_Lambda_y"), v0.rapidity(1), v0.phi());
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

        registryRho.fill(HIST("h2_rho1_AntiLambda"), v0.eta(), v0.phi());
        registryRho.fill(HIST("h2_rho1_AntiLambda_y"), v0.rapidity(2), v0.phi());
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
        if (trk.hasTPC())
          fillProtonTpcQA<true>(ProtonQAStage::RawAfterTrackSel, trk);
        if (trk.hasTOF())
          fillProtonTofQA<true>(ProtonQAStage::RawAfterTrackSel, trk);
      } else {
        fillProtonKinematics<false>(ProtonQAStage::RawAfterTrackSel, trk);
        fillProtonDetectorQuality<false>(ProtonQAStage::RawAfterTrackSel, trk);
        if (trk.hasTPC())
          fillProtonTpcQA<false>(ProtonQAStage::RawAfterTrackSel, trk);
        if (trk.hasTOF())
          fillProtonTofQA<false>(ProtonQAStage::RawAfterTrackSel, trk);
      }
      // ─────────────────────────────────────────────────────────────────

      // ── Track quality cuts ────────────────────────────────────────────
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 1.0f); // Quality cut starting point (after kinematics)
      if (!trk.hasTPC()) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 2.0f); // passed hasTPC
      if (trk.tpcNClsCrossedRows() < pProtonTPCMinRows.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 3.0f); // passed tpcNClsCrossedRows
      if (trk.tpcCrossedRowsOverFindableCls() < pProtonTPCMinRowsOverFindable.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 4.0f); // passed tpcCrossedRowsOverFindableCls
      if (trk.tpcChi2NCl() > pProtonTPCMaxChi2PerCluster.value) {
        continue;
      }
      registryQaDetector.fill(HIST("Selection/TOF_Matching/hCutFlow"), 5.0f); // passed tpcChi2NCl
      if (trk.itsNCls() < pProtonITSMinClusters.value) {
        continue;
      }
      if (trk.itsChi2NCl() > pProtonITSMaxChi2PerCluster.value) {
        continue;
      }
      if (trk.tpcSignal() < cPPMinTpcSignal.value) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // ── p_TPC and pT acceptance cuts ─────────────────────────────────────────────
      const float pTPC = trk.tpcInnerParam();
      if (pTPC < pProtonMinP.value || pTPC > pProtonMaxP.value) {
        continue;
      }
      if (trk.pt() < pProtonMinP.value || trk.pt() > pProtonMaxP.value) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // ── Proton kinematic acceptance ───────────────────────────────────
      // ETACUT-DISABLED (temporary): eta cut commented out, keeping rapidity cut only.
      // if (std::abs(trk.eta()) > pProtonMaxEta.value) { continue; }
      const float y = protonRapidity(trk);
      if (std::abs(y) > pProtonMaxY.value) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────
      // ── DCA cuts ─────────────────────────────────────────────────────
      const bool passDCA = (std::abs(trk.dcaXY()) < pProtonMaxDCAxy.value) &&
                           (std::abs(trk.dcaZ()) < pProtonMaxDCAz.value);
      if (!passDCA) {
        continue;
      }
      // ─────────────────────────────────────────────────────────────────

      // FIX-4: ProtonCounts_byPID_Check — protons only (sign > 0 guard added)
      if (passDCA && trk.sign() > 0) {
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
          if (trk.hasTOF())
            fillProtonTofQA<true>(ProtonQAStage::TpcPid, trk);
        } else {
          fillProtonKinematics<false>(ProtonQAStage::TpcPid, trk);
          fillProtonDetectorQuality<false>(ProtonQAStage::TpcPid, trk);
          fillProtonTpcQA<false>(ProtonQAStage::TpcPid, trk);
          if (trk.hasTOF())
            fillProtonTofQA<false>(ProtonQAStage::TpcPid, trk);
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
      const bool isV0Daughter = v0DaughterTrackIds.count(trk.globalIndex()) > 0;
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
      // ─────────────────────────────────────────────────────────────────

      // ── Classify into proton / antiproton lists ───────────────────────
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
        if (pTPC > pProtonTPCTOFSwitchP.value)
          fillProtonTofQA<true>(ProtonQAStage::FinalSelected, trk);

        registryOther.fill(HIST("hPtSelectedPrimProton"), trk.pt());
        hPtSelectedPrimProton_pT015toMaxDefined->Fill(trk.pt());
      } else {
        selectedPrimAntiProtons.push_back(trk);
        fillProtonKinematics<false>(ProtonQAStage::FinalSelected, trk);
        fillProtonDetectorQuality<false>(ProtonQAStage::FinalSelected, trk);
        fillProtonTpcQA<false>(ProtonQAStage::FinalSelected, trk);
        // TOF QA only when track was accepted through the TPC+TOF branch:
        // p > pProtonTPCTOFSwitchP guarantees hasTOF() and nSigmaTOF passed (beta window cut removed). //BETACUTCommented
        if (pTPC > pProtonTPCTOFSwitchP.value)
          fillProtonTofQA<false>(ProtonQAStage::FinalSelected, trk);

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
    }
    // ─────────────────────────────────────────────────────────────────────

    // ── Single-particle density for protons (eta and y) ──────────────────
    for (auto const& primProton : selectedPrimProtons) {
      registryRho.fill(HIST("h2_rho1_Proton"), primProton.eta(), primProton.phi());
      const float yp = protonRapidity(primProton);
      registryRho.fill(HIST("h2_rho1_Proton_y"), yp, primProton.phi());
    }
    for (auto const& primAntiProton : selectedPrimAntiProtons) {
      registryRho.fill(HIST("h2_rho1_AntiProton"), primAntiProton.eta(), primAntiProton.phi());
      const float yap = protonRapidity(primAntiProton);
      registryRho.fill(HIST("h2_rho1_AntiProton_y"), yap, primAntiProton.phi());
    }
    // ─────────────────────────────────────────────────────────────────────

    registryOther.fill(HIST("hNSelectedPrimProtons"), static_cast<float>(selectedPrimProtons.size()));
    registryOther.fill(HIST("hNSelectedPrimAntiProtons"), static_cast<float>(selectedPrimAntiProtons.size()));

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
        // QA3/PairCleaning: accepted pair
        registryCorrelationQA.fill(HIST("QA3/PairCleaning/h1f_pairCleaning"), 3);

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
            if (qinvPP > kQinvCutLP) {
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

        const int idx1 = unrolledIndex(p1.eta(), p1.phi());
        const int idx2 = unrolledIndex(p2.eta(), p2.phi());
        const float y2_pp = protonRapidity(p2);
        const int idxY1 = unrolledIndexY(y1_pp, p1.phi());
        const int idxY2 = unrolledIndexY(y2_pp, p2.phi());
        // Physics pair correlations (rho2): q_inv cut not applied to rho2
        if (idx1 >= 0 && idx2 >= 0) {
          hRho2_pp_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_pp_y_pT015toMaxDefined->Fill(idxY1, idxY2);
        }
      }
    }
    for (auto const& p1 : selectedPrimProtons) {
      const float y1_pAp = protonRapidity(p1);
      for (auto const& p2 : selectedPrimAntiProtons) {
        const int idx1 = unrolledIndex(p1.eta(), p1.phi());
        const int idx2 = unrolledIndex(p2.eta(), p2.phi());
        const float y2_pAp = protonRapidity(p2);
        const int idxY1 = unrolledIndexY(y1_pAp, p1.phi());
        const int idxY2 = unrolledIndexY(y2_pAp, p2.phi());
        if (idx1 >= 0 && idx2 >= 0) {
          hRho2_pAp_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_pAp_y_pT015toMaxDefined->Fill(idxY1, idxY2);
        }
      }
    }
    for (auto const& p1 : selectedPrimAntiProtons) {
      const float y1_App = protonRapidity(p1);
      for (auto const& p2 : selectedPrimProtons) {
        const int idx1 = unrolledIndex(p1.eta(), p1.phi());
        const int idx2 = unrolledIndex(p2.eta(), p2.phi());
        const float y2_App = protonRapidity(p2);
        const int idxY1 = unrolledIndexY(y1_App, p1.phi());
        const int idxY2 = unrolledIndexY(y2_App, p2.phi());
        if (idx1 >= 0 && idx2 >= 0) {
          hRho2_App_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_App_y_pT015toMaxDefined->Fill(idxY1, idxY2);
        }
      }
    }
    for (auto const& p1 : selectedPrimAntiProtons) {
      const float y1_ApAp = protonRapidity(p1);
      for (auto const& p2 : selectedPrimAntiProtons) {
        if (p1.index() == p2.index()) {
          continue;
        } // CHANGED from globalIndex()
        const int idx1 = unrolledIndex(p1.eta(), p1.phi());
        const int idx2 = unrolledIndex(p2.eta(), p2.phi());
        const float y2_ApAp = protonRapidity(p2);
        const int idxY1 = unrolledIndexY(y1_ApAp, p1.phi());
        const int idxY2 = unrolledIndexY(y2_ApAp, p2.phi());
        if (idx1 >= 0 && idx2 >= 0) {
          hRho2_ApAp_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_ApAp_y_pT015toMaxDefined->Fill(idxY1, idxY2);
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
        if (qinvLP > kQinvCutLP) {
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
        // Physics pair correlations (rho2 / pT / pair-count): q_inv cut not applied to rho2
        if (idxLambda >= 0 && idxPrimP >= 0) {
          if (bin >= 0) {
            hPtPrimProton_Lambda[bin]->Fill(primProton.pt());
          }
          hRho2_Lp_pT015toMaxDefined->Fill(idxLambda, idxPrimP);
        }
        if (idxYL >= 0 && idxYP >= 0) {
          hRho2_Lp_y_pT015toMaxDefined->Fill(idxYL, idxYP);
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
        if (qinvLAp > kQinvCutLP) {
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
        // Physics pair correlations (rho2 / pT / pair-count): q_inv cut not applied to rho2
        if (idxLambda >= 0 && idxPrimPbar >= 0) {
          if (bin >= 0) {
            hPtPrimAntiProton_Lambda[bin]->Fill(primAntiProton.pt());
          }
          hRho2_LAp_pT015toMaxDefined->Fill(idxLambda, idxPrimPbar);
        }
        if (idxYL2 >= 0 && idxYAp >= 0) {
          hRho2_LAp_y_pT015toMaxDefined->Fill(idxYL2, idxYAp);
        }
        ++nPairs_Lambda_PrimAntiProton;
      }

      // ── Autocorrelation diagnostic: primary proton vs Lambda daughter proton──
      // Lambda proton daughter = v0PosDau (positive track, identified as proton
      // by |tpcNSigmaPr| < PProtonTPCNsigma, already in scope above).
      // Loop over ALL selected primary protons — no daughter-track veto applied
      // intentionally (the purpose is to DETECT the autocorrelation).
      // Fills inclusive h2f_dEta_dPhi for all pairs and
      // h2f_dEta_dPhi_nonSameTrack for pairs where the primary proton is not
      // the Lambda positive daughter.
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
        if (qinvALp > kQinvCutLP) {
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
        // Physics pair correlations (rho2 / pT / pair-count): q_inv cut not applied to rho2
        if (idxAL >= 0 && idxPrimP >= 0) {
          if (bin >= 0) {
            hPtPrimProton_AntiLambda[bin]->Fill(primProton.pt());
          }
          hRho2_ALp_pT015toMaxDefined->Fill(idxAL, idxPrimP);
        }
        if (idxYAL >= 0 && idxYP2 >= 0) {
          hRho2_ALp_y_pT015toMaxDefined->Fill(idxYAL, idxYP2);
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
        if (qinvALAp > kQinvCutLP) {
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
        // Physics pair correlations (rho2 / pT / pair-count): q_inv cut not applied to rho2
        if (idxAL2 >= 0 && idxPrimPb >= 0) {
          if (bin >= 0) {
            hPtPrimAntiProton_AntiLambda[bin]->Fill(primAntiProton.pt());
          }
          hRho2_ALAp_pT015toMaxDefined->Fill(idxAL2, idxPrimPb);
        }
        if (idxYAL2 >= 0 && idxYAp2 >= 0) {
          hRho2_ALAp_y_pT015toMaxDefined->Fill(idxYAL2, idxYAp2);
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
      // Fills inclusive h2f_dEta_dPhi for all pairs and
      // h2f_dEta_dPhi_nonSameTrack for pairs where the primary antiproton is
      // not the AntiLambda negative daughter.
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
          if (qinvLL > kQinvCutLP) {
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut"), qinvLL);
          }
        }
        const int idx1 = unrolledIndex(v1.eta(), v1.phi());
        const int idx2 = unrolledIndex(v2.eta(), v2.phi());
        const int idxY1 = unrolledIndexY(yv1L, v1.phi());
        const int idxY2 = unrolledIndexY(v2.rapidity(1), v2.phi());
        // Physics pair correlations (rho2): q_inv cut not applied to rho2
        if (idx1 >= 0 && idx2 >= 0) {
          hRho2_LL_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_LL_y_pT015toMaxDefined->Fill(idxY1, idxY2);
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
          if (qinvLAL > kQinvCutLP) {
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut"), qinvLAL);
          }
        }
        const int idx1 = unrolledIndex(v1.eta(), v1.phi());
        const int idx2 = unrolledIndex(v2.eta(), v2.phi());
        const int idxY1 = unrolledIndexY(yv1LAL, v1.phi());
        const int idxY2 = unrolledIndexY(v2.rapidity(2), v2.phi());
        // Physics pair correlations (rho2): q_inv cut not applied to rho2
        if (idx1 >= 0 && idx2 >= 0) {
          hRho2_LAL_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_LAL_y_pT015toMaxDefined->Fill(idxY1, idxY2);
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
          if (qinvALL > kQinvCutLP) {
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut"), qinvALL);
          }
        }
        const int idx1 = unrolledIndex(v1.eta(), v1.phi());
        const int idx2 = unrolledIndex(v2.eta(), v2.phi());
        const int idxY1 = unrolledIndexY(yv1ALL, v1.phi());
        const int idxY2 = unrolledIndexY(v2.rapidity(1), v2.phi());
        // Physics pair correlations (rho2): q_inv cut not applied to rho2
        if (idx1 >= 0 && idx2 >= 0) {
          hRho2_ALL_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_ALL_y_pT015toMaxDefined->Fill(idxY1, idxY2);
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
          if (qinvALAL > kQinvCutLP) {
            registryCorrelationQA.fill(HIST("QA3/SplitTrackQA/h1f_qinv_LL_After_qinvCut"), qinvALAL);
          }
        }
        const int idx1 = unrolledIndex(v1.eta(), v1.phi());
        const int idx2 = unrolledIndex(v2.eta(), v2.phi());
        const int idxY1 = unrolledIndexY(yv1ALAL, v1.phi());
        const int idxY2 = unrolledIndexY(v2.rapidity(2), v2.phi());
        // Physics pair correlations (rho2): q_inv cut not applied to rho2
        if (idx1 >= 0 && idx2 >= 0) {
          hRho2_ALAL_pT015toMaxDefined->Fill(idx1, idx2);
        }
        if (idxY1 >= 0 && idxY2 >= 0) {
          hRho2_ALAL_y_pT015toMaxDefined->Fill(idxY1, idxY2);
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
    LOG(info) << "===============================================";
    LOG(info) << "";
  }
  // ─────────────────────────────────────────────────────────────────────
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LambdaProtonBalanceFunction>(cfgc)};
}

// What i did here is to change the eta acceptance from 0.8 to 1.0 for now for QA purposes, but later for final analysis we will change it back to 0.8. So please ignore the changes related to maxEta and qaEtaAxis for now. I have kept the maxEta as 1.0f so that it can be used for track acceptance in the code, but we will change it to 0.8f later for final analysis. The same goes for qaEtaAxis, I have set it to -1.0f to 1.0f for now, but we will change it to -0.8f to 0.8f later for final analysis.
//  Also we previously gave 0.8 rapidity cut for lambda and anti-lambda, but now we have changed it to 1.0 for now for QA purposes, but later for final analysis we will change it back to 0.8. So please ignore the changes related to rapidity cuts for now.
