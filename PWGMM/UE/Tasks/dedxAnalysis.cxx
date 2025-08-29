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
///
/// \author Paola Vargas Torres  (paola.vargas.torres@cern.ch)
/// \since January 8, 2025
/// \file dedxAnalysis.cxx
/// \brief  Analysis to do PID

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace constants::physics;

using PIDTracks = soa::Join<
  aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection,
  aod::pidTOFFullPi, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFbeta>;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

struct DedxAnalysis {

  // dE/dx for all charged particles
  HistogramRegistry registryDeDx{
    "registryDeDx",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  // Constant values
  static constexpr int kEtaIntervals = 8;
  static constexpr int kParticlesType = 4;
  float tpcCut = 0.6;
  float pionMin = 0.35;
  float pionMax = 0.45;
  float elTofCut = 0.1;
  float pionTofCut = 1.0;
  float invMassCut = 0.01;
  float invMassCutGamma = 0.0015;
  float magField = 1;
  float pTcut = 2.0;

  // Configurable Parameters
  // Tracks cuts
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 70.0f,
                                      "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of found TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f,
                                 "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f,
                                 "max chi2 per cluster ITS"};
  Configurable<float> maxZDistanceToIP{"maxZDistanceToIP", 10.0f,
                                       "max z distance to IP"};
  Configurable<float> etaMin{"etaMin", -0.8f, "etaMin"};
  Configurable<float> etaMax{"etaMax", +0.8f, "etaMax"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxDCAz{"maxDCAz", 2.f, "maxDCAz"};
  // v0 cuts
  Configurable<float> v0cospaMin{"v0cospaMin", 0.998f, "Minimum V0 CosPA"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f,
                                      "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100.0f,
                                      "Maximum V0 Radius"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f,
                                        "Maximum DCA Daughters"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", 3.0f, "Maximum nsigma TOF"};
  Configurable<float> minMassK0s{"minMassK0s", 0.4f, "Minimum Mass K0s"};
  Configurable<float> maxMassK0s{"maxMassK0s", 0.6f, "Maximum Mass K0s"};
  Configurable<float> minMassLambda{"minMassLambda", 1.1f,
                                    "Minimum Mass Lambda"};
  Configurable<float> maxMassLambda{"maxMassLambda", 1.2f,
                                    "Maximum Mass Lambda"};
  Configurable<float> minMassGamma{"minMassGamma", 0.000922f,
                                   "Minimum Mass Gamma"};
  Configurable<float> maxMassGamma{"maxMassGamma", 0.002022f,
                                   "Maximum Mass Gamma"};
  Configurable<float> nclCut{"nclCut", 135.0f,
                             "ncl Cut"};
  Configurable<bool> calibrationMode{"calibrationMode", false, "calibration mode"};
  Configurable<bool> additionalCuts{"additionalCuts", true, "additional cuts"};
  // Histograms names
  static constexpr std::string_view kDedxvsMomentumPos[kParticlesType] = {"dEdx_vs_Momentum_all_Pos", "dEdx_vs_Momentum_Pi_v0_Pos", "dEdx_vs_Momentum_Pr_v0_Pos", "dEdx_vs_Momentum_El_v0_Pos"};
  static constexpr std::string_view kDedxvsMomentumNeg[kParticlesType] = {"dEdx_vs_Momentum_all_Neg", "dEdx_vs_Momentum_Pi_v0_Neg", "dEdx_vs_Momentum_Pr_v0_Neg", "dEdx_vs_Momentum_El_v0_Neg"};
  static constexpr std::string_view kNclDedxMomentumNegBefore[kEtaIntervals] = {"Ncl_vs_dEdx_vs_Momentum_Neg_1_Before", "Ncl_vs_dEdx_vs_Momentum_Neg_2_Before", "Ncl_vs_dEdx_vs_Momentum_Neg_3_Before", "Ncl_vs_dEdx_vs_Momentum_Neg_4_Before", "Ncl_vs_dEdx_vs_Momentum_Neg_5_Before", "Ncl_vs_dEdx_vs_Momentum_Neg_6_Before", "Ncl_vs_dEdx_vs_Momentum_Neg_7_Before", "Ncl_vs_dEdx_vs_Momentum_Neg_8_Before"};
  static constexpr std::string_view kNclDedxMomentumPosBefore[kEtaIntervals] = {"Ncl_vs_dEdx_vs_Momentum_Pos_1_Before", "Ncl_vs_dEdx_vs_Momentum_Pos_2_Before", "Ncl_vs_dEdx_vs_Momentum_Pos_3_Before", "Ncl_vs_dEdx_vs_Momentum_Pos_4_Before", "Ncl_vs_dEdx_vs_Momentum_Pos_5_Before", "Ncl_vs_dEdx_vs_Momentum_Pos_6_Before", "Ncl_vs_dEdx_vs_Momentum_Pos_7_Before", "Ncl_vs_dEdx_vs_Momentum_Pos_8_Before"};
  static constexpr std::string_view kNclDedxMomentumNegAfter[kEtaIntervals] = {"Ncl_vs_dEdx_vs_Momentum_Neg_1_After", "Ncl_vs_dEdx_vs_Momentum_Neg_2_After", "Ncl_vs_dEdx_vs_Momentum_Neg_3_After", "Ncl_vs_dEdx_vs_Momentum_Neg_4_After", "Ncl_vs_dEdx_vs_Momentum_Neg_5_After", "Ncl_vs_dEdx_vs_Momentum_Neg_6_After", "Ncl_vs_dEdx_vs_Momentum_Neg_7_After", "Ncl_vs_dEdx_vs_Momentum_Neg_8_After"};
  static constexpr std::string_view kNclDedxMomentumPosAfter[kEtaIntervals] = {"Ncl_vs_dEdx_vs_Momentum_Pos_1_After", "Ncl_vs_dEdx_vs_Momentum_Pos_2_After", "Ncl_vs_dEdx_vs_Momentum_Pos_3_After", "Ncl_vs_dEdx_vs_Momentum_Pos_4_After", "Ncl_vs_dEdx_vs_Momentum_Pos_5_After", "Ncl_vs_dEdx_vs_Momentum_Pos_6_After", "Ncl_vs_dEdx_vs_Momentum_Pos_7_After", "Ncl_vs_dEdx_vs_Momentum_Pos_8_After"};
  static constexpr double EtaCut[kEtaIntervals + 1] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
  Configurable<std::vector<float>> calibrationFactorNeg{"calibrationFactorNeg", {50.4011, 50.4764, 50.186, 49.2955, 48.8222, 49.4273, 49.9292, 50.0556}, "negative calibration factors"};
  Configurable<std::vector<float>> calibrationFactorPos{"calibrationFactorPos", {50.5157, 50.6359, 50.3198, 49.3345, 48.9197, 49.4931, 50.0188, 50.1406}, "positive calibration factors"};
  ConfigurableAxis binP{"binP", {VARIABLE_WIDTH, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}, ""};

  // phi cut fits
  TF1* fphiCutHigh = nullptr;
  TF1* fphiCutLow = nullptr;

  TrackSelection myTrackSelection()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetPtRange(0.1f, 1e10f);
    selectedTracks.SetEtaRange(etaMin, etaMax);
    selectedTracks.SetRequireITSRefit(true);
    selectedTracks.SetRequireTPCRefit(true);
    selectedTracks.SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC);
    selectedTracks.SetMaxChi2PerClusterTPC(maxChi2TPC);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1});
    selectedTracks.SetMaxChi2PerClusterITS(maxChi2ITS);
    selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / std::pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(maxDCAz);
    selectedTracks.SetRequireGoldenChi2(true);

    return selectedTracks;
  }

  TrackSelection mySelectionPrim;

  void init(InitContext const&)
  {
    AxisSpec dedxAxis{100, 0.0, 100.0, "dE/dx (a. u.)"};
    AxisSpec ptAxis = {binP, "pT (GeV/c)"};
    AxisSpec etaAxis{8, -0.8, 0.8, "#eta"};
    AxisSpec pAxis = {binP, "#it{p}/Z (GeV/c)"};
    fphiCutLow = new TF1("StandardPhiCutLow", "0.119297/x/x+pi/18.0-0.000379693", 0, 50);
    fphiCutHigh = new TF1("StandardPhiCutHigh", "0.16685/x+pi/18.0+0.00981942", 0, 50);
    if (calibrationMode) {
      // MIP for pions
      registryDeDx.add(
        "hdEdx_vs_eta_Neg_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_Pos_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});
      // MIP for electrons
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      // Pions from TOF
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

    } else {
      // MIP for pions
      registryDeDx.add(
        "hdEdx_vs_eta_Neg_calibrated_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_Pos_calibrated_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});

      // MIP for electrons
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_calibrated_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_calibrated_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      // Pions from TOF
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_calibrated_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_calibrated_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      // pt vs p
      registryDeDx.add(
        "hp_vs_pt_all_Neg", "p_vs_pT", HistType::kTH2F,
        {{ptAxis}, {pAxis}});
      registryDeDx.add(
        "hp_vs_pt_all_Pos", "p_vs_pT", HistType::kTH2F,
        {{ptAxis}, {pAxis}});

      // De/Dx for ch and v0 particles
      for (int i = 0; i < kParticlesType; ++i) {
        registryDeDx.add(kDedxvsMomentumPos[i].data(), "dE/dx", HistType::kTH3F,
                         {{pAxis}, {dedxAxis}, {etaAxis}});
        registryDeDx.add(kDedxvsMomentumNeg[i].data(), "dE/dx", HistType::kTH3F,
                         {{pAxis}, {dedxAxis}, {etaAxis}});
      }
    }

    registryDeDx.add(
      "hdEdx_vs_phi", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 6.4, "#phi"}, {dedxAxis}});

    // phi cut
    registryDeDx.add(
      "hpt_vs_phi_Ncl_After", "phi cut", HistType::kTH3F,
      {{ptAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl}"}});

    registryDeDx.add(
      "hpt_vs_phi_Ncl_Before", "phi cut", HistType::kTH3F,
      {{ptAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl}"}});

    // Ncl vs de/dx

    for (int i = 0; i < kEtaIntervals; ++i) {
      registryDeDx.add(kNclDedxMomentumPosBefore[i].data(), "Ncl vs dE/dx vs Momentum Positive before", HistType::kTH3F,
                       {{100, 0, 160, "N_{cl}"}, {dedxAxis}, {pAxis}});
      registryDeDx.add(kNclDedxMomentumNegBefore[i].data(), "Ncl vs dE/dx vs Momentum Negative before", HistType::kTH3F,
                       {{100, 0, 160, "N_{cl}"}, {dedxAxis}, {pAxis}});

      registryDeDx.add(kNclDedxMomentumPosAfter[i].data(), "Ncl vs dE/dx vs Momentum Positive after", HistType::kTH3F,
                       {{100, 0, 160, "N_{cl}"}, {dedxAxis}, {pAxis}});
      registryDeDx.add(kNclDedxMomentumNegAfter[i].data(), "Ncl vs dE/dx vs Momentum Negative after", HistType::kTH3F,
                       {{100, 0, 160, "N_{cl}"}, {dedxAxis}, {pAxis}});
    }

    // beta plot
    registryDeDx.add(
      "hbeta_vs_p_Neg", "beta", HistType::kTH2F,
      {{pAxis}, {100, 0.0, 1.1, "#beta"}});

    registryDeDx.add(
      "hbeta_vs_p_Pos", "beta", HistType::kTH2F,
      {{pAxis}, {100, 0.0, 1.1, "#beta"}});

    // Event Counter
    registryDeDx.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{100, -20.0, +20.0, "z_{vtx} (cm)"}});

    mySelectionPrim = myTrackSelection();
  }

  // Single-Track Selection
  template <typename T1, typename C>
  bool passedSingleTrackSelection(const T1& track, const C& /*collision*/)
  {
    // Single-Track Selections
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;

    return true;
  }

  // General V0 Selections
  template <typename T1, typename C>
  bool passedV0Selection(const T1& v0, const C& /*collision*/)
  {
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    return true;
  }

  // K0s Selections
  template <typename T1, typename T2, typename C>
  bool passedK0Selection(const T1& v0, const T2& ntrack, const T2& ptrack,
                         const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mK0Short() < minMassK0s || v0.mK0Short() > maxMassK0s)
      return false;

    return true;
  }

  // Lambda Selections
  template <typename T1, typename T2, typename C>
  bool passedLambdaSelection(const T1& v0, const T2& ntrack, const T2& ptrack,
                             const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(ptrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mLambda() < minMassLambda || v0.mLambda() > maxMassLambda)
      return false;

    return true;
  }

  // AntiLambda Selections
  template <typename T1, typename T2, typename C>
  bool passedAntiLambdaSelection(const T1& v0, const T2& ntrack,
                                 const T2& ptrack, const C& collision)
  {

    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mAntiLambda() < minMassLambda || v0.mAntiLambda() > maxMassLambda)
      return false;

    return true;
  }

  // Gamma Selections
  template <typename T1, typename T2, typename C>
  bool passedGammaSelection(const T1& v0, const T2& ntrack, const T2& ptrack,
                            const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(ptrack.tofNSigmaEl()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaEl()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mGamma() < minMassGamma || v0.mGamma() > maxMassGamma)
      return false;

    return true;
  }

  // Phi cut
  template <typename T>
  bool passedPhiCut(const T& trk, float magField, const TF1& fphiCutLow, const TF1& fphiCutHigh)
  {
    float pt = trk.pt();
    float phi = trk.phi();
    int charge = trk.sign();
    float eta = trk.eta();
    auto nTPCCl = trk.tpcNClsFindable() - trk.tpcNClsFindableMinusFound();
    float sigP = trk.sign() * trk.tpcInnerParam();

    if (pt < pTcut)
      return true;

    if (magField < 0) // for negatve polarity field
      phi = o2::constants::math::TwoPI - phi;
    if (charge < 0) // for negatve charge
      phi = o2::constants::math::TwoPI - phi;

    // to center gap in the middle
    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);

    registryDeDx.fill(HIST("hpt_vs_phi_Ncl_Before"), pt, phi, nTPCCl);

    // cut phi
    if (phi < fphiCutHigh.Eval(pt) && phi > fphiCutLow.Eval(pt))
      return false; // reject track

    if (eta > EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegBefore[0]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosBefore[0]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegBefore[1]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosBefore[1]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegBefore[2]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosBefore[2]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegBefore[3]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosBefore[3]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegBefore[4]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosBefore[4]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegBefore[5]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosBefore[5]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegBefore[6]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosBefore[6]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegBefore[7]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosBefore[7]), nTPCCl, trk.tpcSignal(), sigP);
      }
    }

    // cut Ncl
    if (nTPCCl < nclCut)
      return false;

    registryDeDx.fill(HIST("hpt_vs_phi_Ncl_After"), pt, phi, nTPCCl);

    if (eta > EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegAfter[0]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosAfter[0]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegAfter[1]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosAfter[1]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegAfter[2]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosAfter[2]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegAfter[3]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosAfter[3]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegAfter[4]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosAfter[4]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegAfter[5]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosAfter[5]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegAfter[6]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosAfter[6]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclDedxMomentumNegAfter[7]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclDedxMomentumPosAfter[7]), nTPCCl, trk.tpcSignal(), sigP);
      }
    }

    return true;
  }

  // Phi cut Secondaries
  template <typename T>
  bool passedPhiCutSecondaries(const T& trk, float magField, const TF1& fphiCutLow, const TF1& fphiCutHigh)
  {
    float pt = trk.pt();
    float phi = trk.phi();
    int charge = trk.sign();
    auto nTPCCl = trk.tpcNClsFindable() - trk.tpcNClsFindableMinusFound();

    if (pt < pTcut)
      return true;

    if (magField < 0) // for negatve polarity field
      phi = o2::constants::math::TwoPI - phi;
    if (charge < 0) // for negatve charge
      phi = o2::constants::math::TwoPI - phi;

    // to center gap in the middle
    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);

    // cut phi
    if (phi < fphiCutHigh.Eval(pt) && phi > fphiCutLow.Eval(pt))
      return false; // reject track

    // cut Ncl
    if (nTPCCl < nclCut)
      return false;

    return true;
  }

  // Process Data
  void process(SelectedCollisions::iterator const& collision,
               aod::V0Datas const& fullV0s, PIDTracks const& tracks)
  {
    // Event Selection
    if (!collision.sel8())
      return;

    if (additionalCuts) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
        return;

      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
        return;

      if (std::abs(collision.posZ()) >= maxZDistanceToIP)
        return;

      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
        return;
    }

    // Event Counter
    registryDeDx.fill(HIST("histRecVtxZData"), collision.posZ());

    // Kaons
    for (const auto& trk : tracks) {

      // track Selection
      if (!passedSingleTrackSelection(trk, collision))
        continue;

      if (!mySelectionPrim.IsSelected(trk))
        continue;

      // phi and Ncl cut
      if (!passedPhiCut(trk, magField, *fphiCutLow, *fphiCutHigh))
        continue;

      float signedP = trk.sign() * trk.tpcInnerParam();

      // MIP calibration for pions
      if (trk.tpcInnerParam() >= pionMin && trk.tpcInnerParam() <= pionMax) {
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_Neg_Pi"), trk.eta(), trk.tpcSignal());
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_Pos_Pi"), trk.eta(), trk.tpcSignal());
          }

        } else {
          for (int i = 0; i < kEtaIntervals; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_Neg_calibrated_Pi"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_Pos_calibrated_Pi"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i));
              }
            }
          }
        }
      }
      // Beta from TOF
      if (signedP < 0) {
        registryDeDx.fill(HIST("hbeta_vs_p_Neg"), std::abs(signedP), trk.beta());
      } else {
        registryDeDx.fill(HIST("hbeta_vs_p_Pos"), signedP, trk.beta());
      }
      // Electrons from TOF
      if (std::abs(trk.beta() - 1) < elTofCut) { // beta cut
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_El"), trk.eta(), trk.tpcSignal(), std::abs(signedP));
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_El"), trk.eta(), trk.tpcSignal(), signedP);
          }
        } else {
          for (int i = 0; i < kEtaIntervals; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_calibrated_El"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), std::abs(signedP));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_calibrated_El"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i), signedP);
              }
            }
          }
        }
      }
      // pions from TOF
      if (trk.beta() > pionTofCut && trk.beta() < pionTofCut + 0.05) { // beta cut
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_TOF"), trk.eta(), trk.tpcSignal(), std::abs(signedP));
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_TOF"), trk.eta(), trk.tpcSignal(), signedP);
          }
        } else {
          for (int i = 0; i < kEtaIntervals; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_calibrated_TOF"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), std::abs(signedP));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_calibrated_TOF"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i), signedP);
              }
            }
          }
        }
      }

      registryDeDx.fill(HIST("hdEdx_vs_phi"), trk.phi(), trk.tpcSignal());

      if (!calibrationMode) {
        for (int i = 0; i < kEtaIntervals; ++i) {
          if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
            if (signedP > 0) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[0]), signedP, trk.tpcSignal() * 50 / calibrationFactorPos->at(i), trk.eta());
              registryDeDx.fill(HIST("hp_vs_pt_all_Pos"), trk.pt(), signedP);
            } else {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[0]), std::abs(signedP), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), trk.eta());
              registryDeDx.fill(HIST("hp_vs_pt_all_Neg"), trk.pt(), std::abs(signedP));
            }
          }
        }
      }
    }

    // Loop over Reconstructed V0s
    if (!calibrationMode) {
      for (const auto& v0 : fullV0s) {

        // Standard V0 Selections
        if (!passedV0Selection(v0, collision)) {
          continue;
        }

        if (v0.dcaV0daughters() > dcaV0DaughtersMax) {
          continue;
        }

        // Positive and Negative Tracks
        const auto& posTrack = v0.posTrack_as<PIDTracks>();
        const auto& negTrack = v0.negTrack_as<PIDTracks>();

        if (!posTrack.passedTPCRefit())
          continue;
        if (!negTrack.passedTPCRefit())
          continue;
        // phi and Ncl cut
        if (!passedPhiCutSecondaries(posTrack, magField, *fphiCutLow, *fphiCutHigh))
          continue;

        if (!passedPhiCutSecondaries(negTrack, magField, *fphiCutLow, *fphiCutHigh))
          continue;

        float signedPpos = posTrack.sign() * posTrack.tpcInnerParam();
        float signedPneg = negTrack.sign() * negTrack.tpcInnerParam();

        float pxPos = posTrack.px();
        float pyPos = posTrack.py();
        float pzPos = posTrack.pz();

        float pxNeg = negTrack.px();
        float pyNeg = negTrack.py();
        float pzNeg = negTrack.pz();

        const float gammaMass = 2 * MassElectron; // GeV/c^2

        // K0s Selection
        if (passedK0Selection(v0, negTrack, posTrack, collision)) {
          float ePosPi = posTrack.energy(MassPionCharged);
          float eNegPi = negTrack.energy(MassPionCharged);

          float invMass = std::sqrt((eNegPi + ePosPi) * (eNegPi + ePosPi) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassK0Short) > invMassCut) {
            continue;
          }

          for (int i = 0; i < kEtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[1]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[1]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // Lambda Selection
        if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {

          float ePosPr = posTrack.energy(MassProton);
          float eNegPi = negTrack.energy(MassPionCharged);

          float invMass = std::sqrt((eNegPi + ePosPr) * (eNegPi + ePosPr) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassLambda) > invMassCut) {
            continue;
          }

          for (int i = 0; i < kEtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[1]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[2]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // AntiLambda Selection
        if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {

          float ePosPi = posTrack.energy(MassPionCharged);
          float eNegPr = negTrack.energy(MassProton);

          float invMass = std::sqrt((eNegPr + ePosPi) * (eNegPr + ePosPi) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassLambda) > invMassCut) {
            continue;
          }

          for (int i = 0; i < kEtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[2]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[1]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // Gamma Selection
        if (passedGammaSelection(v0, negTrack, posTrack, collision)) {

          float ePosEl = posTrack.energy(MassElectron);
          float eNegEl = negTrack.energy(MassElectron);

          float invMass = std::sqrt((eNegEl + ePosEl) * (eNegEl + ePosEl) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - gammaMass) > invMassCutGamma) {
            continue;
          }

          for (int i = 0; i < kEtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[3]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[3]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DedxAnalysis>(cfgc)};
}
