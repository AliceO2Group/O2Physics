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
/// \brief this is a starting point for the Resonances tutorial
/// \author sourav kundu
/// \since 02/11/2023

#include "PWGLF/DataModel/ReducedDoublePhiTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct doublephimeson {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  struct : ConfigurableGroup {
    Configurable<int> strategyPID1{"strategyPID1", 0, "PID strategy 1"};
    Configurable<int> strategyPID2{"strategyPID2", 0, "PID strategy 2"};
    Configurable<float> daughterDeltaR{"daughterDeltaR", 0.0, "delta R of daughter"};
    Configurable<float> minPhiMass1{"minPhiMass1", 1.01, "Minimum phi mass1"};
    Configurable<float> maxPhiMass1{"maxPhiMass1", 1.03, "Maximum phi mass1"};
    Configurable<float> minPhiPt{"minPhiPt", 0, "Minimum phi Pt"};
    Configurable<float> maxPhiPt{"maxPhiPt", 100, "Maximum phi Pt"};
    Configurable<float> minPhiMass2{"minPhiMass2", 1.01, "Minimum phi mass2"};
    Configurable<float> maxPhiMass2{"maxPhiMass2", 1.03, "Maximum phi mass2"};
    Configurable<float> minExoticPt{"minExoticPt", 6.0, "Minimum Exotic Pt"};
    Configurable<float> minExoticMass{"minExoticMass", 2.0, "Minimum Exotic mass"};
    Configurable<float> maxExoticMass{"maxExoticMass", 3.6, "Maximum Exotic mass"};
    Configurable<bool> additionalEvsel{"additionalEvsel", true, "Additional event selection"};
    Configurable<bool> isDeep{"isDeep", true, "Store deep angle"};
    Configurable<float> cutMinNsigmaTPC{"cutMinNsigmaTPC", -2.5, "nsigma cut TPC"};
    Configurable<float> cutNsigmaTPC{"cutNsigmaTPC", 2.5, "nsigma cut TPC"};
    Configurable<float> cutNsigmaTOF{"cutNsigmaTOF", 2.5, "nsigma cut TOF"};
    Configurable<float> momTOFCut{"momTOFCut", 1.8, "minimum pT cut for madnatory TOF"};
    Configurable<float> maxKaonPt{"maxKaonPt", 100.0, "maximum kaon pt cut"};
    Configurable<float> cfgCrossPhiLow{"cfgCrossPhiLow", 1.01, "Lower edge of phi mass window for cross-pairing (ghost) veto"};
    Configurable<float> cfgCrossPhiHigh{"cfgCrossPhiHigh", 1.03, "Upper edge of phi mass window for cross-pairing (ghost) veto"};
    Configurable<bool> useParametrized{"useParametrized", false, "Use pT dependent mass peak and width"};
    Configurable<bool> useCrossPairRejection{"useCrossPairRejection", true, "Use cross pair phi signal compatibilaty"};
    Configurable<int> cRotations{"cRotations", 10, "Number of rotations for rotational background"};
    Configurable<bool> applyMomentumShift{"applyMomentumShift", false, "Apply momentum shift to kaons to check effect on phi mass peak"};
  } config;

  Configurable<bool> cfgFillDataDrivenPhiResolution{
    "cfgFillDataDrivenPhiResolution", true,
    "Fill the single-phi daughter-kinematics sparse used for data-driven X resolution and inter-dataset momentum-scale calibration"};

  Configurable<bool> cfgFillSelectedXResolutionSparse{
    "cfgFillSelectedXResolutionSparse", true,
    "Fill compact selected-X candidate sparse for data-driven X mass-resolution study"};
  Configurable<float> cfgSelectedXResolutionMassLow{
    "cfgSelectedXResolutionMassLow", 2.63f,
    "Lower M(phi-phi) edge for selected-X resolution sparse (GeV/c^2)"};
  Configurable<float> cfgSelectedXResolutionMassHigh{
    "cfgSelectedXResolutionMassHigh", 2.75f,
    "Upper M(phi-phi) edge for selected-X resolution sparse (GeV/c^2)"};
  Configurable<float> cfgSelectedXResolutionPtMin{
    "cfgSelectedXResolutionPtMin", 9.0f,
    "Minimum pT(phi-phi) for selected-X resolution sparse (GeV/c)"};

  // Optional 2025 -> 2026 kaon momentum-scale correction used only in
  // processopti5. The selected calibration model is a sigmoid plus a quadratic
  // tail:
  //
  // epsilon_corr(%) = p0 + p1/[1 + exp(-(pT_phi-p2)/p3)]
  //                   + p4*pT_phi + p5*pT_phi^2,
  //
  // with pT_phi, p2 and p3 in GeV/c. Both kaons belonging to a phi candidate
  // receive the same scale factor:
  // p_corr = (1 + epsilon_corr/100) * p_measured.
  // The numerical parameters should be supplied from the JSON file generated
  // by CalibratePhiMomentumScaleVsPt.C after the final calibration fit.
  Configurable<bool> cfgApplyKaonMomentumCorrection{
    "cfgApplyKaonMomentumCorrection", false,
    "Apply the pT-dependent kaon momentum-scale correction in processopti5"};
  Configurable<double> cfgMomCorrP0Percent{
    "cfgMomCorrP0Percent", 0.749950,
    "p0 of epsilon_corr(pT_phi) in percent"};
  Configurable<double> cfgMomCorrP1Percent{
    "cfgMomCorrP1Percent", 0.421824,
    "p1 sigmoid-step amplitude of epsilon_corr in percent"};
  Configurable<double> cfgMomCorrP2GeV{
    "cfgMomCorrP2GeV", 1.614360,
    "p2 sigmoid turn-on position in GeV/c"};
  Configurable<double> cfgMomCorrP3GeV{
    "cfgMomCorrP3GeV", 0.290147,
    "p3 sigmoid width in GeV/c; must be positive"};
  Configurable<double> cfgMomCorrP4PercentPerGeV{
    "cfgMomCorrP4PercentPerGeV", 0.0279158,
    "p4 linear-tail slope in percent per GeV/c"};
  Configurable<double> cfgMomCorrP5PercentPerGeV2{
    "cfgMomCorrP5PercentPerGeV2", 0.0,
    "p5 quadratic-tail coefficient in percent per (GeV/c)^2"};
  Configurable<double> cfgMomCorrPtMin{
    "cfgMomCorrPtMin", 0.8,
    "Minimum uncorrected phi pT for applying the momentum correction (GeV/c)"};
  Configurable<double> cfgMomCorrPtMax{
    "cfgMomCorrPtMax", 20.0,
    "Maximum uncorrected phi pT for applying the momentum correction (GeV/c)"};
  // ------------------------------------------------------------
  // pT-dependent phi mass peak and width from single-phi BW fits
  //
  // DeltaM = sqrt( ((m1 - mean(pt1))/width(pt1))^2
  //              + ((m2 - mean(pt2))/width(pt2))^2 )
  //
  // Units:
  //   pT    : GeV/c
  //   mean  : GeV/c^2
  //   width : GeV/c^2
  //
  // 35 pT-bin edges -> 34 mean/width values.
  // ------------------------------------------------------------

  Configurable<std::vector<float>> cfgPhiPtBins{
    "cfgPhiPtBins",
    std::vector<float>{
      0.4, 0.5, 0.6, 0.8, 0.9,
      1.0, 1.1, 1.2, 1.4, 1.6,
      1.8, 2.0, 2.2, 2.4, 2.6,
      2.8, 3.0, 3.5, 4.0, 4.5,
      5.0, 6.0, 7.0, 8.0, 9.0,
      10.0, 12.0, 14.0, 16.0, 18.0,
      20.0, 25.0, 30.0, 40.0, 100.0},
    "pT bin edges for phi mass peak and width calibration"};

  Configurable<std::vector<float>> cfgPhiMeanVsPt{
    "cfgPhiMeanVsPt",
    std::vector<float>{
      1.01779, 1.01826, 1.01897, 1.01924, 1.01934,
      1.01938, 1.01942, 1.01945, 1.01946, 1.01947,
      1.01953, 1.01960, 1.01965, 1.01967, 1.01971,
      1.01972, 1.01974, 1.01977, 1.01980, 1.01983,
      1.01987, 1.01988, 1.01989, 1.01990, 1.01992,
      1.01992, 1.01990, 1.01990, 1.01990, 1.01990,
      1.01990, 1.01990, 1.01990, 1.01990},
    "phi mass peak vs pT from single-phi Breit-Wigner fit"};

  Configurable<std::vector<float>> cfgPhiSigmaVsPt{
    "cfgPhiSigmaVsPt",
    std::vector<float>{
      0.00507421, 0.00554344, 0.00554447, 0.00562718, 0.00548766,
      0.00541166, 0.00539718, 0.00539168, 0.00532416, 0.00534485,
      0.00547586, 0.00567904, 0.00579755, 0.00583511, 0.00587245,
      0.00600147, 0.00600020, 0.00610160, 0.00628558, 0.00646913,
      0.00678283, 0.00720555, 0.00753636, 0.00784207, 0.00814550,
      0.00854121, 0.00890474, 0.00981538, 0.01011060, 0.01077820,
      0.01112900, 0.01203120, 0.01372720, 0.02443550},
    "phi Breit-Wigner width vs pT from single-phi fit"};

  Configurable<float> cfgDefaultPhiMean{
    "cfgDefaultPhiMean",
    1.019461f,
    "default phi mass peak if pT is outside calibration range"};

  Configurable<float> cfgDefaultPhiSigma{
    "cfgDefaultPhiSigma",
    0.0055f,
    "default phi width if pT is outside calibration range"};

  Configurable<float> cfgMinAllowedPhiSigma{
    "cfgMinAllowedPhiSigma",
    1.0e-6f,
    "protection against zero or negative phi width"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 1, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0, 20.0, 40.0, 60.0, 80.0, 500.0}, "Mixing bins - number of contributor"};

  // THnsparse bining
  ConfigurableAxis configThnAxisPtCorr{"configThnAxisPtCorr", {1000, 0.0, 100}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {400, 2.5, 2.9}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMassPhi{"configThnAxisInvMassPhi", {20, 1.01, 1.03}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMassDeltaPhi{"configThnAxisInvMassDeltaPhi", {80, 0.0, 0.08}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMassDeltaPhiSigma{"configThnAxisInvMassDeltaPhiSigma", {100, 0.0, 10}, "#it{M} (GeV/#it{c}^{2}) sigma"};
  ConfigurableAxis configThnAxisDaugherPt{"configThnAxisDaugherPt", {25, 0.0, 50.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {40, 0.0, 20.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisDaughterPt{"configThnAxisDaughterPt", {100, 0.0, 100.}, "daughter #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisKstar{"configThnAxisKstar", {200, 0.0, 2.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisDeltaR{"configThnAxisDeltaR", {VARIABLE_WIDTH, 0.0, 0.0001, 0.0003, 0.0005, 0.0007, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 10.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosTheta{"configThnAxisCosTheta", {10, 0.0, 1.0}, "cos #theta{*}"};
  ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {10, 0.0, 1.0}, "Rapidity"};
  ConfigurableAxis configThnAxisNumPhi{"configThnAxisNumPhi", {6, -0.5, 5.5}, "Number of Kaon daughter have TOF hit"};
  ConfigurableAxis configThnAxisDeltaPt{"configThnAxisDeltaPt", {100, 0.0, 1.0}, "delta pt"};
  Configurable<float> maxDeltaMPhi{"maxDeltaMPhi", 0.01f, "Delta-m cut on the two phi masses: sqrt((m1-mPDG)^2 + (m2-mPDG)^2) < maxDeltaMPhi (GeV/c^2)"};
  // --- NEW: steerable axes from JSON ---
  ConfigurableAxis configThnAxisDeltaRPhi{"configThnAxisDeltaRPhi", {120, 0.0, 6.0}, "ΔR(φ,φ)"};
  ConfigurableAxis configThnAxisZ{"configThnAxisZ", {100, 0.0, 1.0}, "z"};
  ConfigurableAxis configThnAxisA{"configThnAxisA", {10, 0.0, 1.0}, "A"};
  ConfigurableAxis configThnAxisPhiPtVertex{"configThnAxisPhiPtVertex", {100, 0.0, 100.0}, "phi pT (GeV/c)"};
  ConfigurableAxis configThnAxisDecayLength{"configThnAxisDecayLength", {200, 0.0, 1.0}, "3D decay length (cm)"};
  ConfigurableAxis configThnAxisFitChi2Ndf{"configThnAxisFitChi2Ndf", {200, 0.0, 100.0}, "four-kaon fit chi2/NDF"};
  ConfigurableAxis configThnAxisFitProbability{"configThnAxisFitProbability", {100, 0.0, 1.0}, "four-kaon fit probability"};
  ConfigurableAxis configThnAxisRmsDcaSig{"configThnAxisRmsDcaSig", {300, 0.0, 15.0}, "RMS DCA significance"};

  // Data-driven mass-resolution inputs.
  //
  // 1) PhiMassResolutionDataDriven:
  //    inclusive single-phi calibration, filled after daughter/PID selections
  //    but before the hard phi-mass window.
  //
  // 2) SelectedXResolutionDataDriven:
  //    compact candidate-level sparse filled only for the selected X window.
  //    It keeps the measured joint pT configuration of the two phis and all
  //    four kaons, so no pT(X) reweighting or regenerated daughter-pT topology
  //    is needed later.
  ConfigurableAxis cfgDDPhiMassAxis{
    "cfgDDPhiMassAxis",
    {100, 0.995, 1.045},
    "m(KK) for data-driven resolution (GeV/c^2)"};
  // The pT edges are deliberately fine around 1--3 GeV/c, where the momentum
  // response changes fastest, and coarser at high pT.  eta and DeltaPhi are
  // retained only as validation axes; the primary extraction integrates them.
  ConfigurableAxis cfgDDCalibKaonPtAxis{
    "cfgDDCalibKaonPtAxis",
    {VARIABLE_WIDTH, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
     1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0,
     3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0},
    "kaon pT in single-phi calibration (GeV/c)"};
  ConfigurableAxis cfgDDCalibKaonEtaAxis{
    "cfgDDCalibKaonEtaAxis",
    {10, -1.0, 1.0},
    "kaon eta in single-phi calibration"};
  ConfigurableAxis cfgDDCalibDeltaPhiAxis{
    "cfgDDCalibDeltaPhiAxis",
    {36, -3.141592653589793, 3.141592653589793},
    "Delta phi(K+,K-) in single-phi calibration"};
  ConfigurableAxis cfgDDPhiPtAxis{
    "cfgDDPhiPtAxis",
    {VARIABLE_WIDTH, 0.0, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2,
     1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
     3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
     12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 100.0},
    "phi pT for data-driven resolution (GeV/c)"};

  // Compact pT(X) axis.  The sparse itself is only filled for pT(X)>9 GeV/c,
  // so there is no need to allocate fine bins below the analysis threshold.
  ConfigurableAxis cfgDDSelectedXPtAxis{
    "cfgDDSelectedXPtAxis",
    {VARIABLE_WIDTH, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
     25.0, 30.0, 40.0, 50.0, 70.0, 100.0},
    "pT(phi-phi) for selected-X data-driven resolution (GeV/c)"};

  // The selected-X candidate sparse is intentionally restricted at filling
  // time to the narrow M(phi-phi) window and pT(X) threshold.  M(phi-phi) is
  // therefore NOT stored as an axis, keeping the 9D sparse compact.

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hnsigmaTPCKaonPlusBefore", "hnsigmaTPCKaonPlusBefore", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCKaonMinusBefore", "hnsigmaTPCKaonMinusBefore", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCTOFKaonBefore", "hnsigmaTPCTOFKaonBefore", kTH3F, {{500, -3.0, 3.0f}, {500, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCKaonPlus", "hnsigmaTPCKaonPlus", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCKaonMinus", "hnsigmaTPCKaonMinus", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCTOFKaon", "hnsigmaTPCTOFKaon", kTH3F, {{500, -3.0, 3.0f}, {500, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hPhiMassVsPt", "hPhiMassVsPt", kTH2F, {{40, 1.0, 1.04f}, {1000, 0.0f, 100.0f}});
    histos.add("hPhiMassVsPtShifted", "hPhiMassVsPtShifted", kTH2F, {{40, 1.0, 1.04f}, {1000, 0.0f, 100.0f}});
    histos.add("hPhiMass", "hPhiMass", kTH3F, {{40, 1.0, 1.04f}, {40, 1.0, 1.04f}, {250, 0.0f, 100.0f}});
    histos.add("hPhiMassShifted", "hPhiMassShifted", kTH3F, {{40, 1.0, 1.04f}, {40, 1.0, 1.04f}, {250, 0.0f, 100.0f}});
    histos.add("hPhiMassNormalized", "hPhiMassNormalized", kTH3F, {{100, -10.0, 10.0f}, {100, -10.0, 10.0f}, {250, 0.0f, 100.0f}});
    histos.add("hPhiMass2", "hPhiMass2", kTH2F, {{40, 1.0, 1.04f}, {40, 1.0f, 1.04f}});
    histos.add("hkPlusDeltaetaDeltaPhi", "hkPlusDeltaetaDeltaPhi", kTH2F, {{400, -2.0, 2.0}, {640, -2.0 * TMath::Pi(), 2.0 * TMath::Pi()}});
    histos.add("hkMinusDeltaetaDeltaPhi", "hkMinusDeltaetaDeltaPhi", kTH2F, {{400, -2.0, 2.0}, {640, -2.0 * TMath::Pi(), 2.0 * TMath::Pi()}});
    histos.add("hDeltaRkaonplus", "hDeltaRkaonplus", kTH1F, {{800, 0.0, 8.0}});
    histos.add("hDeltaRkaonminus", "hDeltaRkaonminus", kTH1F, {{800, 0.0, 8.0}});
    histos.add("hPtCorrelation", "hPtCorrelation", kTH2F, {{400, 0.0, 40.0}, {5000, 0.0, 100.0}});
    histos.add("hMassCent", "hMassCent", kTH3F, {{40, 1.0, 1.04f}, {40, 1.0, 1.04f}, {100, 0.0, 100.0}});
    const AxisSpec thnAxisdeltapt{configThnAxisDeltaPt, "Delta pt"};
    const AxisSpec thnAxisdaughterpt{configThnAxisDaughterPt, "Daughter pt"};
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisInvMassPhi{configThnAxisInvMassPhi, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassDeltaPhi{configThnAxisInvMassDeltaPhi, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassDeltaPhiSigma{configThnAxisInvMassDeltaPhiSigma, "#it{M} (GeV/#it{c}^{2}) sigma"};
    const AxisSpec thnAxisDeltaR{configThnAxisDeltaR, "#Delta R)"};
    const AxisSpec thnAxisCosTheta{configThnAxisCosTheta, "cos #theta"};
    const AxisSpec thnAxisRapidity{configThnAxisRapidity, "Rapidity"};
    const AxisSpec thnAxisNumPhi{configThnAxisNumPhi, "Number of phi meson"};
    const AxisSpec thnAxisPtCorr{configThnAxisPtCorr, "Pt Corr var"};
    const AxisSpec thnAxisDeltaRPhi{configThnAxisDeltaRPhi, "#Delta R(#phi,#phi)"};
    const AxisSpec thnAxisZ{configThnAxisZ, "z = p_{T1}/(p_{T1}+p_{T2})"};
    const AxisSpec thnAxisA{configThnAxisA, "A = |p_{T1}-p_{T2}|/(p_{T1}+p_{T2})"};
    AxisSpec axisDoublePhiPID{50, 0.0, 5.0, "max daughter n_{#sigma}^{comb}"};
    const AxisSpec thnAxisDecayLength{configThnAxisDecayLength, "#it{L}_{3D} (cm)"};
    const AxisSpec thnAxisFitChi2Ndf{configThnAxisFitChi2Ndf, "#chi^{2}/NDF"};
    const AxisSpec thnAxisFitProbability{configThnAxisFitProbability, "fit probability"};
    const AxisSpec thnAxisRmsDcaSig{configThnAxisRmsDcaSig, "RMS DCA significance"};

    const AxisSpec ddPhiMassAxis{cfgDDPhiMassAxis, "m_{K^{+}K^{-}} (GeV/c^{2})"};
    const AxisSpec ddCalibKaonPtAxis{cfgDDCalibKaonPtAxis, "p_{T,K} (GeV/c)"};
    const AxisSpec ddCalibKaonEtaAxis{cfgDDCalibKaonEtaAxis, "#eta_{K}"};
    const AxisSpec ddCalibDeltaPhiAxis{cfgDDCalibDeltaPhiAxis, "#Delta#varphi_{K^{+}K^{-}}"};
    const AxisSpec ddPhiPtAxis{cfgDDPhiPtAxis, "p_{T,#phi} (GeV/c)"};

    const AxisSpec ddSelectedXPtAxis{cfgDDSelectedXPtAxis, "p_{T,X} (GeV/c)"};

    histos.add("SEMassUnlike", "SEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisDeltaR, thnAxisPt, thnAxisDeltaR, thnAxisInvMassDeltaPhi, thnAxisPtCorr});
    // histos.add("SEMassLike", "SEMassLike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisDeltaR, thnAxisInvMassPhi, thnAxisInvMassPhi, thnAxisNumPhi});
    histos.add("MEMassUnlike", "MEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisDeltaR, thnAxisPt, thnAxisDeltaR, thnAxisInvMassDeltaPhi, thnAxisPtCorr});
    // --- NEW THnSparse storing (ΔRphi, z, A) WITHOUT applying cuts ---
    histos.add("SEMassUnlike_DeltaRZA", "SEMassUnlike_DeltaRZA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisDeltaRPhi, thnAxisZ, thnAxisA, thnAxisInvMassDeltaPhi});
    histos.add("MEMassUnlike_DeltaRZA", "MEMassUnlike_DeltaRZA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisDeltaRPhi, thnAxisZ, thnAxisA, thnAxisInvMassDeltaPhi});
    histos.add("SEMassUnlike_AllVars", "SEMassUnlike_AllVars", HistType::kTHnSparseF,
               {thnAxisInvMass,   // M(phi-phi)
                thnAxisPt,        // pT(phi-phi)
                thnAxisDeltaRPhi, // DeltaR(phi,phi)
                thnAxisDeltaR,    // min DeltaR(K,K)
                thnAxisInvMassPhi,
                thnAxisInvMassPhi,
                thnAxisInvMassDeltaPhi,
                thnAxisInvMassDeltaPhiSigma,
                thnAxisPtCorr,
                thnAxisNumPhi}); // pT correlation variable

    histos.add("MEMassUnlike_AllVars", "MEMassUnlike_AllVars", HistType::kTHnSparseF,
               {thnAxisInvMass,         // M(phi-phi)
                thnAxisPt,              // pT(phi-phi)
                thnAxisDeltaRPhi,       // DeltaR(phi,phi)
                thnAxisDeltaR,          // min DeltaR(K,K)
                thnAxisInvMassPhi,      // m(phi1)
                thnAxisInvMassPhi,      // m(phi2)
                thnAxisInvMassDeltaPhi, // DeltaM_phi
                thnAxisPtCorr,
                thnAxisNumPhi}); // pT correlation variable

    histos.add("SEMassDoublePhi", "SEMassDoublePhi", HistType::kTHnSparseF,
               {thnAxisInvMass, // M(phi-phi)
                thnAxisPt,      // pT(phi-phi)
                thnAxisA,
                thnAxisRapidity,
                thnAxisInvMassPhi,      // m(phi1)
                thnAxisInvMassPhi,      // m(phi2)
                thnAxisInvMassDeltaPhi, // DeltaM_phi
                thnAxisNumPhi,
                axisDoublePhiPID});

    histos.add("SEMassPhiPhi", "SEMassPhiPhi", HistType::kTHnSparseF,
               {
                 thnAxisInvMass,        // M(phi-phi)
                 thnAxisPt,             // pT(phi-phi)
                 thnAxisInvMassDeltaPhi // DeltaM_phi
               });

    histos.add("SEMassPhiPhiRefitted", "SEMassPhiPhiRefitted", HistType::kTHnSparseF,
               {
                 thnAxisInvMass,         // M(phi-phi)
                 thnAxisPt,              // pT(phi-phi)
                 thnAxisInvMassDeltaPhi, // DeltaM_phi
                 thnAxisFitChi2Ndf,      // chi2/NDF of the 4-kaon kinematic fit
                 thnAxisFitProbability,  // fit probability of the 4-kaon kinematic fit
                 thnAxisInvMassPhi,      // m(phi1)
                 thnAxisInvMassPhi       // m(phi2)
               });

    histos.add("SEMassPhiPhiShifted", "SEMassPhiPhiShifted", HistType::kTHnSparseF,
               {
                 thnAxisInvMass,         // M(phi-phi)
                 thnAxisPt,              // pT(phi-phi)
                 thnAxisInvMassDeltaPhi, // DeltaM_phi
                 thnAxisFitChi2Ndf,      // chi2/NDF of the 4-kaon kinematic fit
                 thnAxisFitProbability,  // fit probability of the 4-kaon kinematic fit
                 thnAxisInvMassPhi,      // m(phi1)
                 thnAxisInvMassPhi       // m(phi2)
               });

    histos.add("SEMassPhiPhiRotational", "SEMassPhiPhiRotational", HistType::kTHnSparseF,
               {
                 thnAxisInvMass, // M(phi-phi)
                 thnAxisPt       // pT(phi-phi)
               });

    histos.add("SEMassUnlike_VertexVars", "SEMassUnlike_VertexVars", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisInvMassDeltaPhi, thnAxisInvMassPhi, thnAxisInvMassPhi, thnAxisDecayLength, thnAxisFitChi2Ndf, thnAxisRmsDcaSig});
    histos.add("NPhiPerEvent", "NPhiPerEvent", HistType::kTH1F, {{20, 0, 20}});
    histos.add("NEvents", "NEvents", HistType::kTH1F, {{2, 0, 2}});

    histos.add("SEMassUnlike_VertexVars", "SEMassUnlike_VertexVars", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisInvMassDeltaPhi, thnAxisInvMassPhi, thnAxisInvMassPhi, thnAxisDecayLength, thnAxisFitChi2Ndf, thnAxisRmsDcaSig});

    // Single-phi calibration input. Axis order:
    //   0 m(KK), 1 pT(K+), 2 pT(K-), 3 eta(K+), 4 eta(K-),
    //   5 DeltaPhi(K+,K-), 6 pT(phi).
    // It is intentionally filled BEFORE minPhiMass/maxPhiMass are applied.
    histos.add("PhiMassResolutionDataDriven", "PhiMassResolutionDataDriven",
               HistType::kTHnSparseF,
               {ddPhiMassAxis, ddCalibKaonPtAxis, ddCalibKaonPtAxis,
                ddCalibKaonEtaAxis, ddCalibKaonEtaAxis,
                ddCalibDeltaPhiAxis, ddPhiPtAxis});

    // Selected-X candidate-level resolution sparse.  Filled only for
    //   cfgSelectedXResolutionMassLow < M(phi-phi) < cfgSelectedXResolutionMassHigh
    //   pT(phi-phi) > cfgSelectedXResolutionPtMin.
    //
    // Axis order:
    //   0 pT(X)
    //   1 m(phi1)
    //   2 m(phi2)
    //   3 pT(phi1)
    //   4 pT(phi2)
    //   5 pT(K+ from phi1)
    //   6 pT(K- from phi1)
    //   7 pT(K+ from phi2)
    //   8 pT(K- from phi2)
    //
    // No eta axes and no M(phi-phi) axis are stored.
    histos.add("SelectedXResolutionDataDriven", "SelectedXResolutionDataDriven",
               HistType::kTHnSparseF,
               {ddSelectedXPtAxis,
                ddPhiMassAxis, ddPhiMassAxis,
                ddPhiPtAxis, ddPhiPtAxis,
                ddCalibKaonPtAxis, ddCalibKaonPtAxis,
                ddCalibKaonPtAxis, ddCalibKaonPtAxis});
  }
  TRandom* rn = new TRandom();

  // get kstar
  TLorentzVector trackSum, PartOneCMS, PartTwoCMS, trackRelK;
  float getkstar(const TLorentzVector& part1,
                 const TLorentzVector& part2)
  {
    // const TLorentzVector trackSum = part1 + part2;
    trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    // TLorentzVector PartOneCMS(part1);
    // TLorentzVector PartTwoCMS(part2);
    PartOneCMS.SetXYZM(part1.Px(), part1.Py(), part1.Pz(), part1.M());
    PartTwoCMS.SetXYZM(part2.Px(), part2.Py(), part2.Pz(), part2.M());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    // const TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;
    trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  float deepangle2(const ROOT::Math::PtEtaPhiMVector& candidate1,
                   const ROOT::Math::PtEtaPhiMVector& candidate2)
  {
    const double pt1 = candidate1.Pt();
    const double pt2 = candidate2.Pt();
    const double pz1 = candidate1.Pz();
    const double pz2 = candidate2.Pz();
    const double p1 = candidate1.P();
    const double p2 = candidate2.P();
    const double angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    return angle;
  }

  float deepangle(const TLorentzVector& candidate1,
                  const TLorentzVector& candidate2)
  {
    const double pt1 = candidate1.Pt();
    const double pt2 = candidate2.Pt();
    const double pz1 = candidate1.Pz();
    const double pz2 = candidate2.Pz();
    const double p1 = candidate1.P();
    const double p2 = candidate2.P();
    const double angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    return angle;
  }

  // get cosTheta
  TLorentzVector daughterCMS;
  ROOT::Math::XYZVector threeVecDauCM, threeVecMother;
  float getCosTheta(const TLorentzVector& mother,
                    const TLorentzVector& daughter)
  {
    threeVecMother = mother.Vect();
    const float beta = mother.Beta();
    const float betax = beta * std::cos(mother.Phi()) * std::sin(mother.Theta());
    const float betay = beta * std::sin(mother.Phi()) * std::sin(mother.Theta());
    const float betaz = beta * std::cos(mother.Theta());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    daughterCMS = boostPRF(daughter);
    threeVecDauCM = daughterCMS.Vect();
    float cosThetaStar = TMath::Abs(threeVecDauCM.Dot(threeVecMother) / std::sqrt(threeVecMother.Mag2()) / std::sqrt(threeVecDauCM.Mag2()));
    return cosThetaStar;
  }

  bool selectionPID(float nsigmaTPC, float nsigmaTOF, int TOFHit, int PIDStrategy, float ptcand)
  {

    if (PIDStrategy == 2003) {
      constexpr float radius2Max = 2.5f * 2.5f;
      const bool hasTOF = (TOFHit == 1);

      if (!hasTOF) {
        // TPC-only branch
        if (ptcand < 0.5f) {
          return std::abs(nsigmaTPC) < config.cutNsigmaTPC;
        } else {
          return nsigmaTPC > -2.0f &&
                 nsigmaTPC < config.cutNsigmaTPC;
        }

      } else {
        // TPC+TOF branch
        const float radius2 =
          nsigmaTPC * nsigmaTPC +
          nsigmaTOF * nsigmaTOF;

        if (radius2 >= radius2Max) {
          return false;
        }

        // Circle only below 2 GeV/c
        if (ptcand < 2.0f) {
          return true;
        }

        float slope = 0.0f;
        float intercept = 0.0f;

        if (ptcand < 2.5f) {
          slope = 1.29f;
          intercept = 3.34f;
        } else if (ptcand < 3.0f) {
          slope = 0.94f;
          intercept = 2.45f;
        } else if (ptcand < 4.0f) {
          slope = 0.62f;
          intercept = 1.55f;
        } else if (ptcand < 5.0f) {
          slope = 0.32f;
          intercept = 1.45f;
        } else if (ptcand < 7.0f) {
          slope = 0.22f;
          intercept = 1.35f;
        } else {
          slope = 0.10f;
          intercept = 1.25f;
        }

        return nsigmaTPC <
               slope * nsigmaTOF + intercept;
      }
    }

    if (PIDStrategy == 2004) {

      const bool hasTOF = (TOFHit == 1);

      if (!hasTOF) {
        // No TOF hit
        return std::abs(nsigmaTPC) < 2.0f;

      } else {
        // TOF hit present

        if (ptcand < 2.0f) {
          // No TPC or TOF PID cut
          return true;
        }

        if (ptcand < 2.5f) {
          return !(nsigmaTPC > 1.5f &&
                   nsigmaTOF < -1.6f);
        }

        if (ptcand < 3.0f) {
          return !(nsigmaTPC > 1.0f &&
                   nsigmaTOF < -1.0f);
        }

        if (ptcand < 3.5f) {
          return !(nsigmaTPC > 1.0f &&
                   nsigmaTOF < -0.8f);
        }

        if (ptcand < 4.5f) {
          return !(nsigmaTPC > 1.0f &&
                   nsigmaTOF < -0.5f);
        }

        // pT >= 4.5 GeV/c
        return !(nsigmaTPC > 1.0f &&
                 nsigmaTOF < 1.0f);
      }
    }

    if (PIDStrategy == 1000) {
      if ((TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) ||
          (TOFHit != 1 && std::abs(nsigmaTPC) < 2.0)) {
        return true;
      }
    }

    if (PIDStrategy == 1001) {
      if ((TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) ||
          (TOFHit != 1 && ptcand < 2.5 && std::abs(nsigmaTPC) < 2.0) ||
          (TOFHit != 1 && ptcand >= 2.5 && nsigmaTPC > -2.0 && nsigmaTPC < 1.0)) {
        return true;
      }
    }

    if (PIDStrategy == 1002) {
      if ((TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) ||
          (TOFHit != 1 && std::abs(nsigmaTPC) < 2.0)) {
        return true;
      }
    }

    if (PIDStrategy == 1003) {
      if (ptcand < 0.5 && TOFHit != 1 && std::abs(nsigmaTPC) < config.cutNsigmaTPC) {
        return true;
      }
      if (ptcand < 0.5 && TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < config.cutNsigmaTPC) {
        return true;
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1 && nsigmaTPC > -2.0 && nsigmaTPC < config.cutNsigmaTPC) {
          return true;
        }
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < config.cutNsigmaTOF) {
          return true;
        }
      }
    }

    if (PIDStrategy == 1004) {
      if (ptcand < 0.5 && TOFHit != 1 && nsigmaTPC > -3.0 && nsigmaTPC < 3.0) {
        return true;
      }
      if (ptcand < 0.5 && TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
        return true;
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        }
      }
    }

    if (PIDStrategy == 1008) {
      if (ptcand < 0.5 && TOFHit != 1 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
        return true;
      }
      if (ptcand < 0.5 && TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
        return true;
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1 && ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.0) {
          return true;
        }
        if (TOFHit != 1 && ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (TOFHit != 1 && ptcand >= 0.7 && ptcand < 1.0 && nsigmaTPC > 0.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (TOFHit != 1 && ptcand >= 1.0 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        }
      }
    }

    if (PIDStrategy == 1005) {
      // low pT: TPC-only
      if (ptcand < 0.5 && nsigmaTPC > -2.0 && nsigmaTPC < 3.0) {
        return true;
      }

      // intermediate pT: require TOF + beta cut + combined TPC-TOF nsigma
      if (ptcand >= 0.5 && ptcand < 5.0) {
        if (TOFHit == 1 && std::sqrt(nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) < 2.5) {
          return true;
        }
      }
      // high pT: TPC-only
      if (ptcand >= 5.0 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
        return true;
      }
    }

    if (PIDStrategy == 1004) {
      if (ptcand < 1.2) {
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        } else if (TOFHit != 1) {
          if (ptcand < 0.5 && nsigmaTPC > -2.0 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > 0.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 1.0 && ptcand < 1.2 && nsigmaTPC > -2.5 && nsigmaTPC < 0.5) {
            return true;
          }
        }
      } else {
        if ((TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) ||
            (TOFHit != 1 && ptcand < 2.0 && nsigmaTPC > -2.5 && nsigmaTPC < 2.0) ||
            (TOFHit != 1 && ptcand > 2.0 && nsigmaTPC > -2.5 && nsigmaTPC < 1.0)) {
          return true;
        }
      }
    }

    if (PIDStrategy == 1007) {

      const bool hasTOF = (TOFHit == 1);
      const bool passTPC_low = (nsigmaTPC > -2.0 && nsigmaTPC < 3.0);
      const bool passTPC_midNoTOF = (nsigmaTPC > -1.5 && nsigmaTPC < 3.0);
      const bool passTPC_highNoTOF = (nsigmaTPC > -2.0 && nsigmaTPC < 2.0);

      const bool passTPC_TOF =
        hasTOF && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5;

      if (ptcand < 0.5) {
        if (passTPC_low)
          return true;
      } else if (ptcand < 0.7) {
        if (hasTOF && passTPC_TOF)
          return true;
        if (!hasTOF && passTPC_midNoTOF)
          return true;
      } else if (ptcand < 1.1) {
        if (passTPC_TOF)
          return true;
      } else {
        if (hasTOF && passTPC_TOF)
          return true;
        if (!hasTOF && passTPC_highNoTOF)
          return true;
      }
    }

    if (PIDStrategy == 100) {
      if (ptcand < 1.2) {
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        } else if (TOFHit != 1) {
          if (ptcand < 0.5 && nsigmaTPC > -2.0 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > 0.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 1.0 && ptcand < 1.2 && nsigmaTPC > -2.5 && nsigmaTPC < 0.5) {
            return true;
          }
        }
      } else {
        if ((TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) ||
            (TOFHit != 1 && nsigmaTPC > -2.5 && nsigmaTPC < 1.0)) {
          return true;
        }
      }
    }

    if (PIDStrategy == 101) {
      if (ptcand < 1.0) {
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        } else if (TOFHit != 1) {
          if (ptcand < 0.5 && nsigmaTPC > -2.0 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > 0.0 && nsigmaTPC < 2.0) {
            return true;
          }
        }
      } else if (ptcand >= 1.0 && ptcand < 2.0 && TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
        return true;
      } else if (ptcand > 2.0) {
        if ((TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) ||
            (TOFHit != 1 && nsigmaTPC > -2.5 && nsigmaTPC < 1.0)) {
          return true;
        }
      }
    }

    if (PIDStrategy == 102) {
      if (TOFHit != 1) {
        if (ptcand < 0.5 && nsigmaTPC > -2.0 && nsigmaTPC < 2.5) {
          return true;
        }
        if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.5) {
          return true;
        }
        if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 2.2 && nsigmaTPC > -2.5 && nsigmaTPC < 1.0) {
          return true;
        }
      }
      if (TOFHit == 1 && ptcand > 0.4 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
        return true;
      }
    }

    // optimized TPC TOF
    if (PIDStrategy == 0) {
      if (ptcand < 0.4) {
        if (nsigmaTPC > -3.0 && nsigmaTPC < 3.0) {
          return true;
        }
      } else if (ptcand >= 0.4 && ptcand < 0.5) {
        if (nsigmaTPC > -2.0 && nsigmaTPC < 3.0) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit == 1) {
        if (ptcand < 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        }
        if (ptcand >= 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit != 1) {
        if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > -0.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 1.0 && ptcand < 1.8 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 1.8 && ptcand < 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.5) {
          return true;
        }
        if (ptcand >= 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.0) {
          return true;
        }
      } else if (ptcand >= 5.0 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
        return true;
      }
    }
    // optimized TPC TOF combined
    if (PIDStrategy == 1) {
      if (ptcand < 0.4) {
        if (nsigmaTPC > config.cutMinNsigmaTPC && nsigmaTPC < config.cutNsigmaTPC) {
          return true;
        }
      } else if (ptcand >= 0.4 && ptcand < 0.5) {
        if (nsigmaTPC > -2.0 && nsigmaTPC < config.cutNsigmaTPC) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit == 1) {
        if (ptcand < 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        }
        if (ptcand >= 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
          return true;
        }
      } else if (ptcand >= 5.0 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
        return true;
      }
    }

    if (PIDStrategy == 2) {
      if (ptcand < 0.5) {
        if (nsigmaTPC > config.cutMinNsigmaTPC && nsigmaTPC < config.cutNsigmaTPC) {
          return true;
        }
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1 && ptcand < config.momTOFCut) {
          if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < config.cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < config.cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < config.cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > -0.0 && nsigmaTPC < config.cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 1.0 && ptcand < 1.8 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 1.8 && ptcand < 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.5) {
            return true;
          }
          if (ptcand >= 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.0) {
            return true;
          }
        }
        if (TOFHit == 1) {
          if (TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2.0) < config.cutNsigmaTOF) {
            return true;
          }
        }
      }
    }
    if (PIDStrategy == 3) {
      if (ptcand < 0.5) {
        if (nsigmaTPC > config.cutMinNsigmaTPC && nsigmaTPC < config.cutNsigmaTPC) {
          return true;
        }
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1) {
          if (nsigmaTPC > config.cutMinNsigmaTPC && nsigmaTPC < config.cutNsigmaTPC) {
            return true;
          }
        }
        if (TOFHit == 1) {
          if (TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2.0) < config.cutNsigmaTOF) {
            return true;
          }
        }
      }
    }
    if (PIDStrategy == 4) {
      // pT <= 0.5 GeV/c (TPC-only)
      if (ptcand <= 0.5f) {
        if (ptcand < 0.4f) {
          if (nsigmaTPC > -2.7f && nsigmaTPC < 2.3f)
            return true;
        } else {
          if (nsigmaTPC > -2.5f && nsigmaTPC < 2.5f)
            return true;
        }
      } else if (ptcand > 0.5f && ptcand < 1.0f) {
        // 0.5 < pT < 1.0 GeV/c
        if (TOFHit == 1) {
          if (TMath::Sqrt(nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) < 2.5f)
            return true;
        } else {
          if (ptcand >= 0.5f && ptcand < 0.6f) {
            if (nsigmaTPC > -1.5f && nsigmaTPC < 2.5f)
              return true;
          } else if (ptcand >= 0.6f && ptcand < 0.7f) {
            if (nsigmaTPC > -0.75f && nsigmaTPC < 2.5f)
              return true;
          } else if (ptcand >= 0.7f && ptcand < 0.8f) {
            if (nsigmaTPC > -0.3f && nsigmaTPC < 2.5f)
              return true;
          } else { // 0.8 - 1.0
            if (nsigmaTPC > 0.0f && nsigmaTPC < 2.5f)
              return true;
          }
        }
      } else if (ptcand >= 1.0f && ptcand < 1.8f) {
        // 1.0 <= pT < 1.8 GeV/c
        if (TOFHit == 1) {
          if (TMath::Sqrt(nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) < 2.5f)
            return true;
        } else {
          if (nsigmaTPC > -2.0f && nsigmaTPC < 0.0f)
            return true;
        }
      } else if (ptcand >= 1.8f && ptcand < 5.0f) {
        // 1.8 < pT < 5.0 GeV/c
        if (TOFHit == 1) {
          if (TMath::Sqrt(nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) < 2.0f)
            return true;
        } else {
          if (nsigmaTPC > -2.5f && nsigmaTPC < 1.0f)
            return true;
        }
      } else {
        // pT >= 5.0 GeV/c (TPC-only, low stats)
        if (nsigmaTPC > -2.5f && nsigmaTPC < 1.5f)
          return true;
      }
    }
    return false;
  }

  struct FourKFitResult {
    bool converged = false;

    double chi2 = -1.0;
    int ndf = 2;
    double probability = -1.0;

    double refittedMass = -1.0;

    TLorentzVector k11;
    TLorentzVector k12;
    TLorentzVector k21;
    TLorentzVector k22;

    double phi1Mass = -1.0;
    double phi2Mass = -1.0;
  };

  FourKFitResult fitFourKaons(
    const TLorentzVector& k11Init,
    const TLorentzVector& k12Init,
    const TLorentzVector& k21Init,
    const TLorentzVector& k22Init,
    double mK,
    double mPhi,
    double sigmaPFrac = 0.01,
    int maxIterations = 20,
    double tolerance = 1e-6)
  {
    FourKFitResult result;

    constexpr int N = 12;
    constexpr int NC = 2;

    // ============================================================
    // Initial state vector
    //
    // x = (px1,py1,pz1,
    //      px2,py2,pz2,
    //      px3,py3,pz3,
    //      px4,py4,pz4)
    // ============================================================

    TMatrixD x(N, 1);
    TMatrixD x0(N, 1);

    const TLorentzVector* pInit[4] = {&k11Init, &k12Init, &k21Init, &k22Init};
    for (int i = 0; i < 4; ++i) {

      x(3 * i + 0, 0) = pInit[i]->Px();
      x(3 * i + 1, 0) = pInit[i]->Py();
      x(3 * i + 2, 0) = pInit[i]->Pz();

      x0(3 * i + 0, 0) = pInit[i]->Px();
      x0(3 * i + 1, 0) = pInit[i]->Py();
      x0(3 * i + 2, 0) = pInit[i]->Pz();
    }

    // ============================================================
    // Approximate covariance matrix
    //
    // This is ONLY for Stage 1/2.
    // Later replace this with the actual O2 track covariance.
    // ============================================================

    TMatrixD V(N, N);
    V.Zero();

    for (int i = 0; i < 4; ++i) {

      const double p = pInit[i]->P();

      const double sigmaP = sigmaPFrac * std::max(p, 1e-6);

      V(3 * i + 0, 3 * i + 0) = sigmaP * sigmaP;
      V(3 * i + 1, 3 * i + 1) = sigmaP * sigmaP;
      V(3 * i + 2, 3 * i + 2) = sigmaP * sigmaP;
    }

    TMatrixD Vinv = V;
    Vinv.Invert();

    // ============================================================
    // Helper to construct a kaon 4-vector from x
    // ============================================================

    auto makeLV = [&](const TMatrixD& xx, int particle) {
      const double px = xx(3 * particle + 0, 0);
      const double py = xx(3 * particle + 1, 0);
      const double pz = xx(3 * particle + 2, 0);

      TLorentzVector p;
      p.SetXYZM(px, py, pz, mK);

      return p;
    };

    // ============================================================
    // Iterative constrained fit
    // ============================================================

    for (int iteration = 0; iteration < maxIterations; ++iteration) {

      TLorentzVector k11 = makeLV(x, 0);
      TLorentzVector k12 = makeLV(x, 1);
      TLorentzVector k21 = makeLV(x, 2);
      TLorentzVector k22 = makeLV(x, 3);
      TLorentzVector phi1 = k11 + k12;
      TLorentzVector phi2 = k21 + k22;
      // ----------------------------------------------------------
      // Constraints
      //
      // f1 = M(phi1)^2 - Mphi^2
      // f2 = M(phi2)^2 - Mphi^2
      // ----------------------------------------------------------

      TMatrixD f(NC, 1);

      f(0, 0) = phi1.M2() - mPhi * mPhi;

      f(1, 0) = phi2.M2() - mPhi * mPhi;

      // ----------------------------------------------------------
      // Check whether already sufficiently close
      // ----------------------------------------------------------

      const double scale = mPhi * mPhi;

      if (std::abs(f(0, 0)) < tolerance * scale &&
          std::abs(f(1, 0)) < tolerance * scale) {

        result.converged = true;
        break;
      }

      // ----------------------------------------------------------
      // Numerical Jacobian
      //
      // A(i,j) = df_i / dx_j
      // ----------------------------------------------------------

      TMatrixD A(NC, N);
      A.Zero();

      for (int j = 0; j < N; ++j) {

        TMatrixD xp = x;
        TMatrixD xm = x;

        const double step = 1e-5 * std::max(std::abs(x(j, 0)), 1.0);

        xp(j, 0) += step;
        xm(j, 0) -= step;

        TLorentzVector p1p = makeLV(xp, 0);
        TLorentzVector p2p = makeLV(xp, 1);
        TLorentzVector p3p = makeLV(xp, 2);
        TLorentzVector p4p = makeLV(xp, 3);
        TLorentzVector p1m = makeLV(xm, 0);
        TLorentzVector p2m = makeLV(xm, 1);
        TLorentzVector p3m = makeLV(xm, 2);
        TLorentzVector p4m = makeLV(xm, 3);

        const double f1p = (p1p + p2p).M2() - mPhi * mPhi;
        const double f1m = (p1m + p2m).M2() - mPhi * mPhi;
        const double f2p = (p3p + p4p).M2() - mPhi * mPhi;
        const double f2m = (p3m + p4m).M2() - mPhi * mPhi;
        A(0, j) = (f1p - f1m) / (2.0 * step);
        A(1, j) = (f2p - f2m) / (2.0 * step);
      }

      // ----------------------------------------------------------
      // C = A V A^T
      // ----------------------------------------------------------

      TMatrixD AT(TMatrixD::kTransposed, A);
      TMatrixD C = A * V * AT;
      TMatrixD Cinv = C;
      Cinv.Invert();

      // ----------------------------------------------------------
      // Delta x
      //
      // dx = -V A^T (A V A^T)^-1 f
      // ----------------------------------------------------------

      TMatrixD dx = -1.0 * V * AT * Cinv * f;

      x += dx;

      // ----------------------------------------------------------
      // Check convergence
      // ----------------------------------------------------------

      double maxCorrection = 0.0;

      for (int j = 0; j < N; ++j) {

        maxCorrection = std::max(maxCorrection, std::abs(dx(j, 0)));
      }

      if (maxCorrection < tolerance) {

        result.converged = true;
        break;
      }
    }

    // ============================================================
    // Construct final particles
    // ============================================================

    result.k11 = makeLV(x, 0);
    result.k12 = makeLV(x, 1);
    result.k21 = makeLV(x, 2);
    result.k22 = makeLV(x, 3);

    TLorentzVector phi1Fit = result.k11 + result.k12;
    TLorentzVector phi2Fit = result.k21 + result.k22;
    TLorentzVector pairFit = result.k11 + result.k12 + result.k21 + result.k22;

    result.phi1Mass = phi1Fit.M();
    result.phi2Mass = phi2Fit.M();
    result.refittedMass = pairFit.M();
    // ============================================================
    // chi2
    //
    // chi2 = (x-x0)^T V^-1 (x-x0)
    // ============================================================

    TMatrixD deltaX = x;
    deltaX -= x0;
    TMatrixD deltaXT(TMatrixD::kTransposed, deltaX);
    TMatrixD chi2Matrix = deltaXT * Vinv * deltaX;

    result.chi2 = chi2Matrix(0, 0);
    result.ndf = NC;

    if (result.converged) {

      result.probability = TMath::Prob(result.chi2, result.ndf);
    }

    return result;
  }

  TLorentzVector exotic, Phid1, Phid2;
  TLorentzVector Phi1kaonplus, Phi1kaonminus, Phi2kaonplus, Phi2kaonminus;
  // TLorentzVector exoticRot, Phid1Rot;

  int getPhiPtCalibBin(float pt)
  {
    const auto& bins = cfgPhiPtBins.value;

    if (bins.size() < 2) {
      return -1;
    }

    if (pt < bins.front()) {
      return -1;
    }

    for (size_t i = 0; i + 1 < bins.size(); ++i) {
      if (pt >= bins[i] && pt < bins[i + 1]) {
        return static_cast<int>(i);
      }
    }

    if (pt >= bins.back()) {
      return static_cast<int>(bins.size()) - 2;
    }

    return -1;
  }

  float getPhiMassPeakVsPt(float pt)
  {
    const int ibin = getPhiPtCalibBin(pt);
    const auto& means = cfgPhiMeanVsPt.value;

    if (ibin >= 0 && ibin < static_cast<int>(means.size())) {
      return means[ibin];
    }

    return cfgDefaultPhiMean.value;
  }

  float getPhiWidthVsPt(float pt)
  {
    const int ibin = getPhiPtCalibBin(pt);
    const auto& widths = cfgPhiSigmaVsPt.value;

    float width = cfgDefaultPhiSigma.value;

    if (ibin >= 0 && ibin < static_cast<int>(widths.size())) {
      width = widths[ibin];
    }

    if (width < cfgMinAllowedPhiSigma.value) {
      width = cfgDefaultPhiSigma.value;
    }

    return width;
  }

  float getNormalizedDeltaMPhi(float m1, float pt1, float m2, float pt2)
  {
    const float mean1 = getPhiMassPeakVsPt(pt1);
    const float mean2 = getPhiMassPeakVsPt(pt2);

    const float width1 = getPhiWidthVsPt(pt1);
    const float width2 = getPhiWidthVsPt(pt2);

    return TMath::Sqrt(
      TMath::Power((m1 - mean1) / width1, 2.0) +
      TMath::Power((m2 - mean2) / width2, 2.0));
  }

  float getDeltaMPhi(float m1, float pt1, float m2, float pt2)
  {
    const float mean1 = getPhiMassPeakVsPt(pt1);
    const float mean2 = getPhiMassPeakVsPt(pt2);

    return TMath::Sqrt(
      TMath::Power((m1 - mean1), 2.0) +
      TMath::Power((m2 - mean2), 2.0));
  }

  float getNormalizedMPhi(float m1, float pt1)
  {
    const float mean1 = getPhiMassPeakVsPt(pt1);
    const float width1 = getPhiWidthVsPt(pt1);
    return ((m1 - mean1) / width1);
  }

  void processSE(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (config.additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    int phimult = 0;
    for (auto const& phitrackd1 : phitracks) {
      if (phitrackd1.phiMass() < config.minPhiMass1 || phitrackd1.phiMass() > config.maxPhiMass1) {
        continue;
      }
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > config.maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > config.maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), config.strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), config.strategyPID2, kaonminusd1pt)) {
        continue;
      }
      phimult = phimult + 1;
    }
    for (auto const& phitrackd1 : phitracks) {
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > config.maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > config.maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), config.strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), config.strategyPID2, kaonminusd1pt)) {
        continue;
      }
      histos.fill(HIST("hnsigmaTPCTOFKaon"), phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), phitrackd1.phid1TPC(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), phitrackd1.phid2TPC(), kaonminusd1pt);
      histos.fill(HIST("hPhiMass2"), Phid1.M(), Phid1.Pt());
      auto phid1id = phitrackd1.index();
      Phid1.SetXYZM(phitrackd1.phiPx(), phitrackd1.phiPy(), phitrackd1.phiPz(), phitrackd1.phiMass());
      Phi1kaonplus.SetXYZM(phitrackd1.phid1Px(), phitrackd1.phid1Py(), phitrackd1.phid1Pz(), 0.493);
      Phi1kaonminus.SetXYZM(phitrackd1.phid2Px(), phitrackd1.phid2Py(), phitrackd1.phid2Pz(), 0.493);
      for (auto const& phitrackd2 : phitracks) {
        auto phid2id = phitrackd2.index();
        if (phid2id <= phid1id) {
          continue;
        }
        auto kaonplusd2pt = TMath::Sqrt(phitrackd2.phid1Px() * phitrackd2.phid1Px() + phitrackd2.phid1Py() * phitrackd2.phid1Py());
        auto kaonminusd2pt = TMath::Sqrt(phitrackd2.phid2Px() * phitrackd2.phid2Px() + phitrackd2.phid2Py() * phitrackd2.phid2Py());
        if (kaonplusd2pt > config.maxKaonPt) {
          continue;
        }
        if (kaonminusd2pt > config.maxKaonPt) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid1TPC(), phitrackd2.phid1TOF(), phitrackd2.phid1TOFHit(), config.strategyPID2, kaonplusd2pt)) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid2TPC(), phitrackd2.phid2TOF(), phitrackd2.phid2TOFHit(), config.strategyPID2, kaonminusd2pt)) {
          continue;
        }
        if (phitrackd1.phid1Index() == phitrackd2.phid1Index()) {
          continue;
        }
        if (phitrackd1.phid2Index() == phitrackd2.phid2Index()) {
          continue;
        }
        Phid2.SetXYZM(phitrackd2.phiPx(), phitrackd2.phiPy(), phitrackd2.phiPz(), phitrackd2.phiMass());
        Phi2kaonplus.SetXYZM(phitrackd2.phid1Px(), phitrackd2.phid1Py(), phitrackd2.phid1Pz(), 0.493);
        Phi2kaonminus.SetXYZM(phitrackd2.phid2Px(), phitrackd2.phid2Py(), phitrackd2.phid2Pz(), 0.493);

        /*
              // Like
              Phid1like = Phi1kaonplus + Phi2kaonplus;
              Phid2like = Phi1kaonminus + Phi2kaonminus;
              exoticlike = Phid1like + Phid2like;
              auto deltaRlike = TMath::Sqrt(TMath::Power(Phid1like.Phi() - Phid2like.Phi(), 2.0) + TMath::Power(Phid1like.Eta() - Phid2like.Eta(), 2.0));
              auto costhetalike = (Phid1like.Px() * Phid2like.Px() + Phid1like.Py() * Phid2like.Py() + Phid1like.Pz() * Phid2like.Pz()) / (Phid1like.P() * Phid2like.P());
              auto deltamlike = TMath::Sqrt(TMath::Power(Phid1like.M() - 1.0192, 2.0) + TMath::Power(Phid2like.M() - 1.0192, 2.0));
              if (!config.isDeep) {
                histos.fill(HIST("SEMassLike"), exoticlike.M(), exoticlike.Pt(), deltaRlike, costhetalike, deltamlike, phimult);
              }
              if (config.isDeep) {
                histos.fill(HIST("SEMassLike"), exoticlike.M(), exoticlike.Pt(), deltaRlike, deepangle(Phid1like, Phid2like), deltamlike, phimult);
              }
        */

        // Unlike
        // histos.fill(HIST("hPhiMass2"), Phid1.M(), Phid2.M());
        if (phitrackd2.phiMass() < config.minPhiMass2 || phitrackd2.phiMass() > config.maxPhiMass2) {
          continue;
        }
        if (phitrackd1.phiMass() < config.minPhiMass1 || phitrackd1.phiMass() > config.maxPhiMass1) {
          continue;
        }
        exotic = Phid1 + Phid2;
        if (exotic.M() < config.minExoticMass || exotic.M() > config.maxExoticMass) {
          continue;
        }
        histos.fill(HIST("hkPlusDeltaetaDeltaPhi"), Phi1kaonplus.Eta() - Phi2kaonplus.Eta(), Phi1kaonplus.Phi() - Phi2kaonplus.Phi());
        histos.fill(HIST("hkMinusDeltaetaDeltaPhi"), Phi1kaonminus.Eta() - Phi2kaonminus.Eta(), Phi1kaonminus.Phi() - Phi2kaonminus.Phi());
        // auto cosThetaStar = getCosTheta(exotic, Phid1);
        // auto kstar = getkstar(Phid1, Phid2);
        auto deltaR = TMath::Sqrt(TMath::Power(Phid1.Phi() - Phid2.Phi(), 2.0) + TMath::Power(Phid1.Eta() - Phid2.Eta(), 2.0));
        auto deltaRd1 = TMath::Sqrt(TMath::Power(Phi1kaonplus.Phi() - Phi2kaonplus.Phi(), 2.0) + TMath::Power(Phi1kaonplus.Eta() - Phi2kaonplus.Eta(), 2.0));
        auto deltaRd2 = TMath::Sqrt(TMath::Power(Phi1kaonminus.Phi() - Phi2kaonminus.Phi(), 2.0) + TMath::Power(Phi1kaonminus.Eta() - Phi2kaonminus.Eta(), 2.0));
        auto deltam = TMath::Sqrt(TMath::Power(Phid1.M() - 1.0192, 2.0) + TMath::Power(Phid2.M() - 1.0192, 2.0));
        if (deltaRd1 < config.daughterDeltaR) {
          continue;
        }
        if (deltaRd2 < config.daughterDeltaR) {
          continue;
        }
        if (!config.isDeep) {
          histos.fill(HIST("SEMassUnlike"), exotic.M(), std::abs(Phid1.Pt() - Phid2.Pt()) / exotic.Pt(), exotic.Pt(), deltaR, deltam, phimult);
        }
        if (config.isDeep) {
          histos.fill(HIST("SEMassUnlike"), exotic.M(), std::abs(Phid1.Pt() - Phid2.Pt()) / exotic.Pt(), exotic.Pt(), deltaR, deltam, phimult);
        }
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processSE, "Process Same Event", false);
  void processopti(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    std::vector<ROOT::Math::PtEtaPhiMVector> exoticresonance, phiresonanced1, phiresonanced2, kaonplus1, kaonplus2, kaonminus1, kaonminus2;
    std::vector<int> d1trackid = {};
    std::vector<int> d2trackid = {};
    std::vector<int> d3trackid = {};
    std::vector<int> d4trackid = {};
    if (config.additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    int phimult = 0;

    for (auto const& phitrackd1 : phitracks) {
      if (phitrackd1.phiMass() < config.minPhiMass1 || phitrackd1.phiMass() > config.maxPhiMass1) {
        continue;
      }
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > config.maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > config.maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), config.strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), config.strategyPID2, kaonminusd1pt)) {
        continue;
      }
      phimult = phimult + 1;
    }
    for (auto const& phitrackd1 : phitracks) {
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());

      histos.fill(HIST("hnsigmaTPCTOFKaonBefore"), phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlusBefore"), phitrackd1.phid1TPC(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinusBefore"), phitrackd1.phid2TPC(), kaonminusd1pt);

      if (kaonplusd1pt > config.maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > config.maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), config.strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), config.strategyPID2, kaonminusd1pt)) {
        continue;
      }
      histos.fill(HIST("hnsigmaTPCTOFKaon"), phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), phitrackd1.phid1TPC(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), phitrackd1.phid2TPC(), kaonminusd1pt);
      histos.fill(HIST("hPhiMass2"), Phid1.M(), Phid1.Pt());
      auto phid1id = phitrackd1.index();
      Phid1.SetXYZM(phitrackd1.phiPx(), phitrackd1.phiPy(), phitrackd1.phiPz(), phitrackd1.phiMass());
      Phi1kaonplus.SetXYZM(phitrackd1.phid1Px(), phitrackd1.phid1Py(), phitrackd1.phid1Pz(), 0.493);
      Phi1kaonminus.SetXYZM(phitrackd1.phid2Px(), phitrackd1.phid2Py(), phitrackd1.phid2Pz(), 0.493);
      for (auto const& phitrackd2 : phitracks) {
        auto phid2id = phitrackd2.index();
        if (phid2id <= phid1id) {
          continue;
        }
        auto kaonplusd2pt = TMath::Sqrt(phitrackd2.phid1Px() * phitrackd2.phid1Px() + phitrackd2.phid1Py() * phitrackd2.phid1Py());
        auto kaonminusd2pt = TMath::Sqrt(phitrackd2.phid2Px() * phitrackd2.phid2Px() + phitrackd2.phid2Py() * phitrackd2.phid2Py());
        if (kaonplusd2pt > config.maxKaonPt) {
          continue;
        }
        if (kaonminusd2pt > config.maxKaonPt) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid1TPC(), phitrackd2.phid1TOF(), phitrackd2.phid1TOFHit(), config.strategyPID2, kaonplusd2pt)) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid2TPC(), phitrackd2.phid2TOF(), phitrackd2.phid2TOFHit(), config.strategyPID2, kaonminusd2pt)) {
          continue;
        }
        if ((phitrackd1.phid1Index() == phitrackd2.phid1Index()) || (phitrackd1.phid2Index() == phitrackd2.phid2Index())) {
          continue;
        }
        Phid2.SetXYZM(phitrackd2.phiPx(), phitrackd2.phiPy(), phitrackd2.phiPz(), phitrackd2.phiMass());
        Phi2kaonplus.SetXYZM(phitrackd2.phid1Px(), phitrackd2.phid1Py(), phitrackd2.phid1Pz(), 0.493);
        Phi2kaonminus.SetXYZM(phitrackd2.phid2Px(), phitrackd2.phid2Py(), phitrackd2.phid2Pz(), 0.493);

        // unlike
        if (phitrackd1.phiMass() < config.minPhiMass1 || phitrackd1.phiMass() > config.maxPhiMass1) {
          continue;
        }
        if (phitrackd2.phiMass() < config.minPhiMass2 || phitrackd2.phiMass() > config.maxPhiMass2) {
          continue;
        }
        exotic = Phid1 + Phid2;
        if (exotic.M() < config.minExoticMass || exotic.M() > config.maxExoticMass) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector temp1(exotic.Pt(), exotic.Eta(), exotic.Phi(), exotic.M());
        ROOT::Math::PtEtaPhiMVector temp2(Phid1.Pt(), Phid1.Eta(), Phid1.Phi(), Phid1.M());
        ROOT::Math::PtEtaPhiMVector temp3(Phid2.Pt(), Phid2.Eta(), Phid2.Phi(), Phid2.M());
        exoticresonance.push_back(temp1);
        phiresonanced1.push_back(temp2);
        phiresonanced2.push_back(temp3);
        d1trackid.push_back(phitrackd1.phid1Index());
        d2trackid.push_back(phitrackd2.phid1Index());
        d3trackid.push_back(phitrackd1.phid2Index());
        d4trackid.push_back(phitrackd2.phid2Index());

        ROOT::Math::PtEtaPhiMVector temp4(Phi1kaonplus.Pt(), Phi1kaonplus.Eta(), Phi1kaonplus.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector temp5(Phi1kaonminus.Pt(), Phi1kaonminus.Eta(), Phi1kaonminus.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector temp6(Phi2kaonplus.Pt(), Phi2kaonplus.Eta(), Phi2kaonplus.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector temp7(Phi2kaonminus.Pt(), Phi2kaonminus.Eta(), Phi2kaonminus.Phi(), 0.493);
        kaonplus1.push_back(temp4);
        kaonplus2.push_back(temp6);
        kaonminus1.push_back(temp5);
        kaonminus2.push_back(temp7);
      }
    }
    if (exoticresonance.size() == 0) {
      return;
    }
    // LOGF(info, "Total number of exotic: %d", exoticresonance.size());
    if (exoticresonance.size() == 2) {
      for (auto if1 = exoticresonance.begin(); if1 != exoticresonance.end(); ++if1) {
        auto i5 = std::distance(exoticresonance.begin(), if1);

        const auto& exotic1phi1 = phiresonanced1.at(i5);
        const auto& exotic1phi2 = phiresonanced2.at(i5);
        const auto& exotic1 = exoticresonance.at(i5);

        const auto& exotic1kaonplus1 = kaonplus1.at(i5);
        const auto& exotic1kaonminus1 = kaonminus1.at(i5);
        const auto& exotic1kaonplus2 = kaonplus2.at(i5);
        const auto& exotic1kaonminus2 = kaonminus2.at(i5);
        auto deltaRkaonplus1 = TMath::Sqrt(TMath::Power(exotic1kaonplus1.Phi() - exotic1kaonplus2.Phi(), 2.0) + TMath::Power(exotic1kaonplus1.Eta() - exotic1kaonplus2.Eta(), 2.0));
        auto deltaRkaonminus1 = TMath::Sqrt(TMath::Power(exotic1kaonminus1.Phi() - exotic1kaonminus2.Phi(), 2.0) + TMath::Power(exotic1kaonminus1.Eta() - exotic1kaonminus2.Eta(), 2.0));
        histos.fill(HIST("hDeltaRkaonplus"), deltaRkaonplus1);
        histos.fill(HIST("hDeltaRkaonminus"), deltaRkaonminus1);

        auto deltam1 = TMath::Sqrt(TMath::Power(exotic1phi1.M() - 1.0192, 2.0) + TMath::Power(exotic1phi2.M() - 1.0192, 2.0));
        auto deltaR1 = TMath::Sqrt(TMath::Power(exotic1phi1.Phi() - exotic1phi2.Phi(), 2.0) + TMath::Power(exotic1phi1.Eta() - exotic1phi2.Eta(), 2.0));

        if (deltaRkaonplus1 < config.daughterDeltaR) {
          continue;
        }
        if (deltaRkaonminus1 < config.daughterDeltaR) {
          continue;
        }

        for (auto if2 = if1 + 1; if2 != exoticresonance.end(); ++if2) {
          auto i6 = std::distance(exoticresonance.begin(), if2);
          const auto& exotic2phi1 = phiresonanced1.at(i6);
          const auto& exotic2phi2 = phiresonanced2.at(i6);
          const auto& exotic2 = exoticresonance.at(i6);

          const auto& exotic2kaonplus1 = kaonplus1.at(i6);
          const auto& exotic2kaonminus1 = kaonminus1.at(i6);
          const auto& exotic2kaonplus2 = kaonplus2.at(i6);
          const auto& exotic2kaonminus2 = kaonminus2.at(i6);
          auto deltaRkaonplus2 = TMath::Sqrt(TMath::Power(exotic2kaonplus1.Phi() - exotic2kaonplus2.Phi(), 2.0) + TMath::Power(exotic2kaonplus1.Eta() - exotic2kaonplus2.Eta(), 2.0));
          auto deltaRkaonminus2 = TMath::Sqrt(TMath::Power(exotic2kaonminus1.Phi() - exotic2kaonminus2.Phi(), 2.0) + TMath::Power(exotic2kaonminus1.Eta() - exotic2kaonminus2.Eta(), 2.0));

          auto deltam2 = TMath::Sqrt(TMath::Power(exotic2phi1.M() - 1.0192, 2.0) + TMath::Power(exotic2phi2.M() - 1.0192, 2.0));
          auto deltaR2 = TMath::Sqrt(TMath::Power(exotic2phi1.Phi() - exotic2phi2.Phi(), 2.0) + TMath::Power(exotic2phi1.Eta() - exotic2phi2.Eta(), 2.0));

          if ((d1trackid.at(i5) == d1trackid.at(i6) || d1trackid.at(i5) == d2trackid.at(i6)) &&
              (d2trackid.at(i5) == d1trackid.at(i6) || d2trackid.at(i5) == d2trackid.at(i6)) &&
              (d3trackid.at(i5) == d3trackid.at(i6) || d3trackid.at(i5) == d4trackid.at(i6)) &&
              (d4trackid.at(i5) == d3trackid.at(i6) || d4trackid.at(i5) == d4trackid.at(i6))) {

            if (deltam2 < deltam1 && deltaRkaonplus2 > config.daughterDeltaR && deltaRkaonminus2 > config.daughterDeltaR) {
              histos.fill(HIST("SEMassUnlike"), exotic2.M(), std::abs(exotic2phi1.Pt() - exotic2phi2.Pt()) / exotic2.Pt(), exotic2.Pt(), deltaR2, deltam2, phimult);
              // LOGF(info, "Fill exotic Id %d which is pair of Id %d", i6, i5);
            } else {
              histos.fill(HIST("SEMassUnlike"), exotic1.M(), std::abs(exotic2phi1.Pt() - exotic2phi2.Pt()) / exotic1.Pt(), exotic1.Pt(), deltaR1, deltam1, phimult);
            }
          } else {
            histos.fill(HIST("SEMassUnlike"), exotic1.M(), std::abs(exotic2phi1.Pt() - exotic2phi2.Pt()) / exotic1.Pt(), exotic1.Pt(), deltaR1, deltam1, phimult);
          }
        }
      }
    } else {
      for (auto if1 = exoticresonance.begin(); if1 != exoticresonance.end(); ++if1) {
        auto i5 = std::distance(exoticresonance.begin(), if1);
        const auto& exotic1phi1 = phiresonanced1.at(i5);
        const auto& exotic1phi2 = phiresonanced2.at(i5);
        const auto& exotic1 = exoticresonance.at(i5);

        const auto& exotic1kaonplus1 = kaonplus1.at(i5);
        const auto& exotic1kaonminus1 = kaonminus1.at(i5);
        const auto& exotic1kaonplus2 = kaonplus2.at(i5);
        const auto& exotic1kaonminus2 = kaonminus2.at(i5);
        auto deltaRkaonplus1 = TMath::Sqrt(TMath::Power(exotic1kaonplus1.Phi() - exotic1kaonplus2.Phi(), 2.0) + TMath::Power(exotic1kaonplus1.Eta() - exotic1kaonplus2.Eta(), 2.0));
        auto deltaRkaonminus1 = TMath::Sqrt(TMath::Power(exotic1kaonminus1.Phi() - exotic1kaonminus2.Phi(), 2.0) + TMath::Power(exotic1kaonminus1.Eta() - exotic1kaonminus2.Eta(), 2.0));
        auto deltam1 = TMath::Sqrt(TMath::Power(exotic1phi1.M() - 1.0192, 2.0) + TMath::Power(exotic1phi2.M() - 1.0192, 2.0));
        auto deltaR1 = TMath::Sqrt(TMath::Power(exotic1phi1.Phi() - exotic1phi2.Phi(), 2.0) + TMath::Power(exotic1phi1.Eta() - exotic1phi2.Eta(), 2.0));

        histos.fill(HIST("hDeltaRkaonplus"), deltaRkaonplus1);
        histos.fill(HIST("hDeltaRkaonminus"), deltaRkaonminus1);

        if (deltaRkaonplus1 < config.daughterDeltaR) {
          continue;
        }
        if (deltaRkaonminus1 < config.daughterDeltaR) {
          continue;
        }

        histos.fill(HIST("SEMassUnlike"), exotic1.M(), std::abs(exotic1phi1.Pt() - exotic1phi2.Pt()) / exotic1.Pt(), exotic1.Pt(), deltaR1, deltam1, phimult);
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processopti, "Process Optimized same event", false);

  void processopti3(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (config.additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2))
      return;

    // --- φ multiplicity with PID ---
    int phimult = 0;
    for (auto const& t : phitracks) {
      const double kpluspt = std::hypot(t.phid1Px(), t.phid1Py());
      const double kminuspt = std::hypot(t.phid2Px(), t.phid2Py());
      // PID QA before
      histos.fill(HIST("hnsigmaTPCTOFKaonBefore"), t.phid1TPC(), t.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlusBefore"), t.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinusBefore"), t.phid2TPC(), kminuspt);
      if (t.phiMass() < config.minPhiMass1 || t.phiMass() > config.maxPhiMass1)
        continue;
      if (kpluspt > config.maxKaonPt || kminuspt > config.maxKaonPt)
        continue;
      if (!selectionPID(t.phid1TPC(), t.phid1TOF(), t.phid1TOFHit(), config.strategyPID1, kpluspt))
        continue;
      if (!selectionPID(t.phid2TPC(), t.phid2TOF(), t.phid2TOFHit(), config.strategyPID2, kminuspt))
        continue;
      // PID QA after
      histos.fill(HIST("hnsigmaTPCTOFKaon"), t.phid1TPC(), t.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), t.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), t.phid2TPC(), kminuspt);

      ++phimult;
    }
    if (phimult < 2)
      return;

    // --- helpers ---
    constexpr double mPhiPDG = 1.019461; // GeV/c^2

    const auto deltaMPhi = [=](double m1, double m2) {
      const double d1 = m1 - mPhiPDG, d2 = m2 - mPhiPDG;
      return std::sqrt(d1 * d1 + d2 * d2);
    };

    const auto deltaR = [](double phi1, double eta1, double phi2, double eta2) {
      const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
      const double deta = eta1 - eta2;
      return std::sqrt(dphi * dphi + deta * deta);
    };

    // minimum ΔR among all kaons in the candidate (4 kaons → 6 combinations)
    const auto minKaonDeltaR = [&](const ROOT::Math::PtEtaPhiMVector& kplusA,
                                   const ROOT::Math::PtEtaPhiMVector& kplusB,
                                   const ROOT::Math::PtEtaPhiMVector& kminusA,
                                   const ROOT::Math::PtEtaPhiMVector& kminusB) {
      // same-sign first (keep your QA histos)
      const double dRkplus = deltaR(kplusA.Phi(), kplusA.Eta(), kplusB.Phi(), kplusB.Eta());
      const double dRkminus = deltaR(kminusA.Phi(), kminusA.Eta(), kminusB.Phi(), kminusB.Eta());
      histos.fill(HIST("hDeltaRkaonplus"), dRkplus);
      histos.fill(HIST("hDeltaRkaonminus"), dRkminus);

      // all other combinations
      const double dR_k1p_k1m = deltaR(kplusA.Phi(), kplusA.Eta(), kminusA.Phi(), kminusA.Eta());
      const double dR_k1p_k2m = deltaR(kplusA.Phi(), kplusA.Eta(), kminusB.Phi(), kminusB.Eta());
      const double dR_k2p_k1m = deltaR(kplusB.Phi(), kplusB.Eta(), kminusA.Phi(), kminusA.Eta());
      const double dR_k2p_k2m = deltaR(kplusB.Phi(), kplusB.Eta(), kminusB.Phi(), kminusB.Eta());
      double minDR = dRkplus;
      minDR = std::min(minDR, dRkminus);
      minDR = std::min(minDR, dR_k1p_k1m);
      minDR = std::min(minDR, dR_k1p_k2m);
      minDR = std::min(minDR, dR_k2p_k1m);
      minDR = std::min(minDR, dR_k2p_k2m);
      return minDR;
    };

    // --- collect candidates once ---
    std::vector<ROOT::Math::PtEtaPhiMVector> pairV, phi1V, phi2V, kplus1V, kplus2V, kminus1V, kminus2V;
    std::vector<double> minDRV; // store minimum ΔR for each pair

    for (auto const& t1 : phitracks) {
      const double kplus1pt = std::hypot(t1.phid1Px(), t1.phid1Py());
      const double kminus1pt = std::hypot(t1.phid2Px(), t1.phid2Py());

      if (kplus1pt > config.maxKaonPt || kminus1pt > config.maxKaonPt)
        continue;
      if (!selectionPID(t1.phid1TPC(), t1.phid1TOF(), t1.phid1TOFHit(), config.strategyPID1, kplus1pt))
        continue;
      if (!selectionPID(t1.phid2TPC(), t1.phid2TOF(), t1.phid2TOFHit(), config.strategyPID2, kminus1pt))
        continue;

      TLorentzVector phi1, k1p, k1m;
      phi1.SetXYZM(t1.phiPx(), t1.phiPy(), t1.phiPz(), t1.phiMass());
      k1p.SetXYZM(t1.phid1Px(), t1.phid1Py(), t1.phid1Pz(), 0.493);
      k1m.SetXYZM(t1.phid2Px(), t1.phid2Py(), t1.phid2Pz(), 0.493);

      // φ mass windows
      if (t1.phiMass() < config.minPhiMass1 || t1.phiMass() > config.maxPhiMass1)
        continue;
      if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt)
        continue;

      const auto id1 = t1.index();

      for (auto const& t2 : phitracks) {
        const auto id2 = t2.index();
        if (id2 <= id1)
          continue;

        const double kplus2pt = std::hypot(t2.phid1Px(), t2.phid1Py());
        const double kminus2pt = std::hypot(t2.phid2Px(), t2.phid2Py());
        if (kplus2pt > config.maxKaonPt || kminus2pt > config.maxKaonPt)
          continue;
        if (!selectionPID(t2.phid1TPC(), t2.phid1TOF(), t2.phid1TOFHit(), config.strategyPID1, kplus2pt))
          continue;
        if (!selectionPID(t2.phid2TPC(), t2.phid2TOF(), t2.phid2TOFHit(), config.strategyPID2, kminus2pt))
          continue;

        // block shared same-sign daughters
        if ((t1.phid1Index() == t2.phid1Index()) || (t1.phid2Index() == t2.phid2Index()))
          continue;

        TLorentzVector phi2, k2p, k2m;
        phi2.SetXYZM(t2.phiPx(), t2.phiPy(), t2.phiPz(), t2.phiMass());
        k2p.SetXYZM(t2.phid1Px(), t2.phid1Py(), t2.phid1Pz(), 0.493);
        k2m.SetXYZM(t2.phid2Px(), t2.phid2Py(), t2.phid2Pz(), 0.493);
        if (t2.phiMass() < config.minPhiMass2 || t2.phiMass() > config.maxPhiMass2)
          continue;
        if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt)
          continue;
        // Δm cut (configurable)
        const double dM = deltaMPhi(phi1.M(), phi2.M());
        if (dM > maxDeltaMPhi)
          continue;

        TLorentzVector pair = phi1 + phi2;
        if (pair.M() < config.minExoticMass || pair.M() > config.maxExoticMass)
          continue;
        histos.fill(HIST("hPhiMass"), phi1.M(), phi2.M(), pair.Pt());
        // daughter ΔR QA and minΔR (NO CUT anymore)
        ROOT::Math::PtEtaPhiMVector k1pV(k1p.Pt(), k1p.Eta(), k1p.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector k1mV(k1m.Pt(), k1m.Eta(), k1m.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector k2pV(k2p.Pt(), k2p.Eta(), k2p.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector k2mV(k2m.Pt(), k2m.Eta(), k2m.Phi(), 0.493);
        const double minDR = minKaonDeltaR(k1pV, k2pV, k1mV, k2mV);

        // store for one-pass fill
        pairV.emplace_back(pair.Pt(), pair.Eta(), pair.Phi(), pair.M());
        phi1V.emplace_back(phi1.Pt(), phi1.Eta(), phi1.Phi(), phi1.M());
        phi2V.emplace_back(phi2.Pt(), phi2.Eta(), phi2.Phi(), phi2.M());
        kplus1V.emplace_back(k1p.Pt(), k1p.Eta(), k1p.Phi(), 0.493);
        kminus1V.emplace_back(k1m.Pt(), k1m.Eta(), k1m.Phi(), 0.493);
        kplus2V.emplace_back(k2p.Pt(), k2p.Eta(), k2p.Phi(), 0.493);
        kminus2V.emplace_back(k2m.Pt(), k2m.Eta(), k2m.Phi(), 0.493);
        minDRV.emplace_back(minDR); // per-candidate minimum ΔR of kaons
      }
    }

    if (pairV.empty())
      return;

    // --- fill the single THnSparse ---
    for (size_t i = 0; i < pairV.size(); ++i) {
      TLorentzVector p1, p2, pair;
      p1.SetPtEtaPhiM(phi1V[i].Pt(), phi1V[i].Eta(), phi1V[i].Phi(), phi1V[i].M());
      p2.SetPtEtaPhiM(phi2V[i].Pt(), phi2V[i].Eta(), phi2V[i].Phi(), phi2V[i].M());
      pair.SetPtEtaPhiM(pairV[i].Pt(), pairV[i].Eta(), pairV[i].Phi(), pairV[i].M());

      const double dM = deltaMPhi(p1.M(), p2.M());
      const double M = pair.M();
      const double dR = deltaR(p1.Phi(), p1.Eta(), p2.Phi(), p2.Eta());
      const double minDR = minDRV[i];
      double ptcorr = p1.Pt() / (pair.Pt() - p1.Pt());
      histos.fill(HIST("hPtCorrelation"), pair.Pt(), ptcorr);
      // NOTE: second axis is now minΔR(all kaons), ΔpT/pT has been removed
      histos.fill(HIST("SEMassUnlike"),
                  M,
                  minDR,
                  pair.Pt(),
                  dR,
                  dM,
                  ptcorr);
    }
  }
  PROCESS_SWITCH(doublephimeson, processopti3, "Process Optimized same event", false);

  void processopti4(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (config.additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2))
      return;

    // --- φ multiplicity with PID ---
    int phimult = 0;
    for (auto const& t : phitracks) {
      const double kpluspt = std::hypot(t.phid1Px(), t.phid1Py());
      const double kminuspt = std::hypot(t.phid2Px(), t.phid2Py());

      // PID QA before
      histos.fill(HIST("hnsigmaTPCTOFKaonBefore"), t.phid1TPC(), t.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlusBefore"), t.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinusBefore"), t.phid2TPC(), kminuspt);

      if (t.phiMass() < config.minPhiMass1 || t.phiMass() > config.maxPhiMass1)
        continue;
      if (kpluspt > config.maxKaonPt || kminuspt > config.maxKaonPt)
        continue;
      if (!selectionPID(t.phid1TPC(), t.phid1TOF(), t.phid1TOFHit(), config.strategyPID1, kpluspt))
        continue;
      if (!selectionPID(t.phid2TPC(), t.phid2TOF(), t.phid2TOFHit(), config.strategyPID2, kminuspt))
        continue;

      // PID QA after
      histos.fill(HIST("hnsigmaTPCTOFKaon"), t.phid1TPC(), t.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), t.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), t.phid2TPC(), kminuspt);

      ++phimult;
    }
    if (phimult < 2)
      return;

    // --- helpers ---
    constexpr double mPhiPDG = o2::constants::physics::MassPhi; // GeV/c^2
    constexpr double mKPDG = o2::constants::physics::MassKPlus; // GeV/c^2

    const auto deltaMPhi = [=](double m1, double m2) {
      const double d1 = m1 - mPhiPDG, d2 = m2 - mPhiPDG;
      return std::sqrt(d1 * d1 + d2 * d2);
    };

    const auto deltaR = [](double phi1, double eta1, double phi2, double eta2) {
      const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
      const double deta = eta1 - eta2;
      return std::sqrt(dphi * dphi + deta * deta);
    };

    // minimum ΔR among all kaons in the candidate (4 kaons → 6 combinations)
    const auto minKaonDeltaR = [&](const ROOT::Math::PtEtaPhiMVector& kplusA,
                                   const ROOT::Math::PtEtaPhiMVector& kplusB,
                                   const ROOT::Math::PtEtaPhiMVector& kminusA,
                                   const ROOT::Math::PtEtaPhiMVector& kminusB) {
      // same-sign first (keep your QA histos)
      const double dRkplus = deltaR(kplusA.Phi(), kplusA.Eta(), kplusB.Phi(), kplusB.Eta());
      const double dRkminus = deltaR(kminusA.Phi(), kminusA.Eta(), kminusB.Phi(), kminusB.Eta());
      histos.fill(HIST("hDeltaRkaonplus"), dRkplus);
      histos.fill(HIST("hDeltaRkaonminus"), dRkminus);

      // all other combinations
      const double dR_k1p_k1m = deltaR(kplusA.Phi(), kplusA.Eta(), kminusA.Phi(), kminusA.Eta());
      const double dR_k1p_k2m = deltaR(kplusA.Phi(), kplusA.Eta(), kminusB.Phi(), kminusB.Eta());
      const double dR_k2p_k1m = deltaR(kplusB.Phi(), kplusB.Eta(), kminusA.Phi(), kminusA.Eta());
      const double dR_k2p_k2m = deltaR(kplusB.Phi(), kplusB.Eta(), kminusB.Phi(), kminusB.Eta());

      double minDR = dRkplus;
      minDR = std::min(minDR, dRkminus);
      minDR = std::min(minDR, dR_k1p_k1m);
      minDR = std::min(minDR, dR_k1p_k2m);
      minDR = std::min(minDR, dR_k2p_k1m);
      minDR = std::min(minDR, dR_k2p_k2m);
      return minDR;
    };

    // --- collect candidates once ---
    std::vector<ROOT::Math::PtEtaPhiMVector> pairV, phi1V, phi2V;
    std::vector<double> minDRV; // store minimum ΔR for each pair

    // optional: swapped-cross-mass veto window (minimal new knobs)
    const double crossPhiLow = 1.01;  // or a dedicated config
    const double crossPhiHigh = 1.03; // or a dedicated config

    for (auto const& t1 : phitracks) {
      const double kplus1pt = std::hypot(t1.phid1Px(), t1.phid1Py());
      const double kminus1pt = std::hypot(t1.phid2Px(), t1.phid2Py());

      if (kplus1pt > config.maxKaonPt || kminus1pt > config.maxKaonPt)
        continue;
      if (!selectionPID(t1.phid1TPC(), t1.phid1TOF(), t1.phid1TOFHit(), config.strategyPID1, kplus1pt))
        continue;
      if (!selectionPID(t1.phid2TPC(), t1.phid2TOF(), t1.phid2TOFHit(), config.strategyPID2, kminus1pt))
        continue;

      TLorentzVector phi1, k1p, k1m;
      phi1.SetXYZM(t1.phiPx(), t1.phiPy(), t1.phiPz(), t1.phiMass());
      k1p.SetXYZM(t1.phid1Px(), t1.phid1Py(), t1.phid1Pz(), mKPDG);
      k1m.SetXYZM(t1.phid2Px(), t1.phid2Py(), t1.phid2Pz(), mKPDG);

      // φ1 mass window + φ1 pT
      if (t1.phiMass() < config.minPhiMass1 || t1.phiMass() > config.maxPhiMass1)
        continue;
      if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt)
        continue;

      const auto id1 = t1.index();

      for (auto const& t2 : phitracks) {
        const auto id2 = t2.index();
        if (id2 <= id1)
          continue;

        const double kplus2pt = std::hypot(t2.phid1Px(), t2.phid1Py());
        const double kminus2pt = std::hypot(t2.phid2Px(), t2.phid2Py());
        if (kplus2pt > config.maxKaonPt || kminus2pt > config.maxKaonPt)
          continue;
        if (!selectionPID(t2.phid1TPC(), t2.phid1TOF(), t2.phid1TOFHit(), config.strategyPID1, kplus2pt))
          continue;
        if (!selectionPID(t2.phid2TPC(), t2.phid2TOF(), t2.phid2TOFHit(), config.strategyPID2, kminus2pt))
          continue;

        // FIX + robust: block ANY shared daughter (4-way)
        if (t1.phid1Index() == t2.phid1Index() || t1.phid1Index() == t2.phid2Index() ||
            t1.phid2Index() == t2.phid1Index() || t1.phid2Index() == t2.phid2Index())
          continue;

        TLorentzVector phi2, k2p, k2m;
        phi2.SetXYZM(t2.phiPx(), t2.phiPy(), t2.phiPz(), t2.phiMass());
        k2p.SetXYZM(t2.phid1Px(), t2.phid1Py(), t2.phid1Pz(), mKPDG);
        k2m.SetXYZM(t2.phid2Px(), t2.phid2Py(), t2.phid2Pz(), mKPDG);

        // φ2 mass window + FIX: apply pT cut to phi2 (not phi1)
        if (t2.phiMass() < config.minPhiMass2 || t2.phiMass() > config.maxPhiMass2)
          continue;
        if (phi2.Pt() < config.minPhiPt || phi2.Pt() > config.maxPhiPt)
          continue;

        // NEW: cross (swapped) K+K- mass veto
        // veto if either cross-pair lands in φ window (tune to be looser/tighter)
        const double mCross12 = (k1p + k2m).M(); // K+1 + K-2
        const double mCross21 = (k2p + k1m).M(); // K+2 + K-1
        if ((mCross12 > crossPhiLow && mCross12 < crossPhiHigh) ||
            (mCross21 > crossPhiLow && mCross21 < crossPhiHigh)) {
          continue;
        }

        // Δm cut (configurable)
        const double dM = deltaMPhi(phi1.M(), phi2.M());
        if (dM > maxDeltaMPhi)
          continue;

        TLorentzVector pair = phi1 + phi2;
        if (pair.M() < config.minExoticMass || pair.M() > config.maxExoticMass)
          continue;

        histos.fill(HIST("hPhiMass"), phi1.M(), phi2.M(), pair.Pt());

        // daughter ΔR QA and minΔR (NO CUT)
        ROOT::Math::PtEtaPhiMVector k1pV(k1p.Pt(), k1p.Eta(), k1p.Phi(), mKPDG);
        ROOT::Math::PtEtaPhiMVector k1mV(k1m.Pt(), k1m.Eta(), k1m.Phi(), mKPDG);
        ROOT::Math::PtEtaPhiMVector k2pV(k2p.Pt(), k2p.Eta(), k2p.Phi(), mKPDG);
        ROOT::Math::PtEtaPhiMVector k2mV(k2m.Pt(), k2m.Eta(), k2m.Phi(), mKPDG);
        const double minDR = minKaonDeltaR(k1pV, k2pV, k1mV, k2mV);

        // store for one-pass fill
        pairV.emplace_back(pair.Pt(), pair.Eta(), pair.Phi(), pair.M());
        phi1V.emplace_back(phi1.Pt(), phi1.Eta(), phi1.Phi(), phi1.M());
        phi2V.emplace_back(phi2.Pt(), phi2.Eta(), phi2.Phi(), phi2.M());
        minDRV.emplace_back(minDR);
      }
    }

    if (pairV.empty())
      return;

    // --- fill the single THnSparse ---
    for (size_t i = 0; i < pairV.size(); ++i) {
      TLorentzVector p1, p2, pair;
      p1.SetPtEtaPhiM(phi1V[i].Pt(), phi1V[i].Eta(), phi1V[i].Phi(), phi1V[i].M());
      p2.SetPtEtaPhiM(phi2V[i].Pt(), phi2V[i].Eta(), phi2V[i].Phi(), phi2V[i].M());
      pair.SetPtEtaPhiM(pairV[i].Pt(), pairV[i].Eta(), pairV[i].Phi(), pairV[i].M());

      const double dM = deltaMPhi(p1.M(), p2.M());
      const double M = pair.M();
      const double dR = deltaR(p1.Phi(), p1.Eta(), p2.Phi(), p2.Eta());
      const double minDR = minDRV[i];

      // (optional but recommended) protect ptcorr from blow-ups
      const double denom = std::abs(pair.Pt() - p1.Pt());

      const double ptcorr = p1.Pt() / denom;

      histos.fill(HIST("hPtCorrelation"), pair.Pt(), ptcorr);

      // NOTE: second axis is minΔR(all kaons), ΔpT/pT has been removed
      histos.fill(HIST("SEMassUnlike"),
                  M,
                  minDR,
                  pair.Pt(),
                  dR,
                  dM,
                  ptcorr);

      // --- NEW: compute z and A from phi candidates (no cuts) ---
      const double pt1 = p1.Pt();
      const double pt2 = p2.Pt();
      const double ptsum = pt1 + pt2;
      if (ptsum <= 0.0)
        continue;
      const double z = pt1 / ptsum;
      const double A = std::abs(pt1 - pt2) / ptsum;
      // --- Fill NEW THnSparse (no cuts) ---
      histos.fill(HIST("SEMassUnlike_DeltaRZA"), M, pair.Pt(), pair.Pt() * dR, z, A, dM);
    }
  }
  PROCESS_SWITCH(doublephimeson, processopti4, "Process Optimized same event", true);
  double dMNominal = 100.0;

  void processopti5(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (config.additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }

    constexpr double mPhiPDG = o2::constants::physics::MassPhi; // GeV/c^2
    constexpr double mKPDG = o2::constants::physics::MassKPlus; // GeV/c^2

    // Build the phi and its daughter-kaon four-vectors. With the correction
    // disabled this reproduces the original opti5 construction exactly. With
    // it enabled, the scale is evaluated using the uncorrected phi pT, both
    // kaon three-momenta are scaled coherently, their energies are recomputed
    // using mKPDG, and the corrected phi is rebuilt from the two kaons.
    const auto buildPhiAndKaons = [&](const auto& t,
                                      TLorentzVector& phi,
                                      TLorentzVector& kplus,
                                      TLorentzVector& kminus) {
      kplus.SetXYZM(t.phid1Px(), t.phid1Py(), t.phid1Pz(), mKPDG);
      kminus.SetXYZM(t.phid2Px(), t.phid2Py(), t.phid2Pz(), mKPDG);
      phi.SetXYZM(t.phiPx(), t.phiPy(), t.phiPz(), t.phiMass());

      if (!cfgApplyKaonMomentumCorrection) {
        return;
      }

      const double uncorrectedPhiPt = std::hypot(t.phiPx(), t.phiPy());
      if (!std::isfinite(uncorrectedPhiPt) || uncorrectedPhiPt <= 0. ||
          uncorrectedPhiPt < cfgMomCorrPtMin ||
          uncorrectedPhiPt > cfgMomCorrPtMax) {
        return;
      }

      const double sigmoidWidth = cfgMomCorrP3GeV;
      if (!std::isfinite(sigmoidWidth) || sigmoidWidth <= 0.) {
        return;
      }

      const double sigmoid =
        1. / (1. + std::exp(-(uncorrectedPhiPt - cfgMomCorrP2GeV) /
                            sigmoidWidth));
      const double epsilonPercent =
        cfgMomCorrP0Percent +
        cfgMomCorrP1Percent * sigmoid +
        cfgMomCorrP4PercentPerGeV * uncorrectedPhiPt +
        cfgMomCorrP5PercentPerGeV2 * uncorrectedPhiPt * uncorrectedPhiPt;
      const double scale = 1. + 0.01 * epsilonPercent;

      if (!std::isfinite(scale) || scale <= 0.) {
        return;
      }

      kplus.SetXYZM(scale * t.phid1Px(),
                    scale * t.phid1Py(),
                    scale * t.phid1Pz(),
                    mKPDG);
      kminus.SetXYZM(scale * t.phid2Px(),
                     scale * t.phid2Py(),
                     scale * t.phid2Pz(),
                     mKPDG);
      phi = kplus + kminus;
    };

    // --- phi multiplicity with PID ---
    int phimult = 0;
    // int nBestPairingRejected = 0;
    for (auto const& t : phitracks) {
      const double kpluspt = std::hypot(t.phid1Px(), t.phid1Py());
      const double kminuspt = std::hypot(t.phid2Px(), t.phid2Py());

      histos.fill(HIST("hnsigmaTPCTOFKaonBefore"), t.phid1TPC(), t.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlusBefore"), t.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinusBefore"), t.phid2TPC(), kminuspt);

      TLorentzVector phi1;
      TLorentzVector kplus;
      TLorentzVector kminus;
      buildPhiAndKaons(t, phi1, kplus, kminus);
      const double phi1Mass = cfgApplyKaonMomentumCorrection ? phi1.M() : t.phiMass();

      // Apply the same daughter/PID requirements before filling the calibration
      // sparse.  The phi-pT requirement for that sparse is evaluated below from
      // the RAW daughters, so it remains independent of any optional correction.
      if (kpluspt > config.maxKaonPt || kminuspt > config.maxKaonPt) {
        continue;
      }
      if (!selectionPID(t.phid1TPC(), t.phid1TOF(), t.phid1TOFHit(), config.strategyPID1, kpluspt)) {
        continue;
      }
      if (!selectionPID(t.phid2TPC(), t.phid2TOF(), t.phid2TOFHit(), config.strategyPID2, kminuspt)) {
        continue;
      }

      if (cfgFillDataDrivenPhiResolution) {
        // IMPORTANT: always store the RAW reconstructed daughter momenta here,
        // independent of cfgApplyKaonMomentumCorrection.  This keeps the
        // resolution calibration data-driven and also allows two data sets to
        // be compared later to infer their relative momentum-scale shift
        // without circularly applying a pre-existing correction first.
        TLorentzVector kplusRaw;
        TLorentzVector kminusRaw;
        kplusRaw.SetXYZM(t.phid1Px(), t.phid1Py(), t.phid1Pz(), mKPDG);
        kminusRaw.SetXYZM(t.phid2Px(), t.phid2Py(), t.phid2Pz(), mKPDG);
        const TLorentzVector phiForResolution = kplusRaw + kminusRaw;
        const double dPhiKK = TVector2::Phi_mpi_pi(kplusRaw.Phi() - kminusRaw.Phi());
        if (phiForResolution.Pt() >= config.minPhiPt && phiForResolution.Pt() <= config.maxPhiPt) {
          histos.fill(HIST("PhiMassResolutionDataDriven"),
                      phiForResolution.M(),
                      kplusRaw.Pt(), kminusRaw.Pt(),
                      kplusRaw.Eta(), kminusRaw.Eta(),
                      dPhiKK, phiForResolution.Pt());
        }
      }

      // From here onward keep the original signal-phi definition unchanged.
      if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt) {
        continue;
      }
      if (phi1Mass < config.minPhiMass1 || phi1Mass > config.maxPhiMass1) {
        continue;
      }

      histos.fill(HIST("hnsigmaTPCTOFKaon"), t.phid1TPC(), t.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), t.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), t.phid2TPC(), kminuspt);
      histos.fill(HIST("hPhiMassVsPt"), phi1Mass, phi1.Pt());

      ++phimult;
    }

    if (phimult < 2) {
      return;
    }

    const auto deltaMPhiNominal = [=](double m1, double m2) {
      const double d1 = m1 - mPhiPDG;
      const double d2 = m2 - mPhiPDG;
      return std::sqrt(d1 * d1 + d2 * d2);
    };

    const auto deltaR = [](double phi1, double eta1, double phi2, double eta2) {
      const double dphi = std::abs(TVector2::Phi_mpi_pi(phi1 - phi2));
      const double deta = eta1 - eta2;
      return std::sqrt(dphi * dphi + deta * deta);
    };

    const auto minKaonDeltaR =
      [&](const ROOT::Math::PtEtaPhiMVector& kplus1,
          const ROOT::Math::PtEtaPhiMVector& kminus1,
          const ROOT::Math::PtEtaPhiMVector& kplus2,
          const ROOT::Math::PtEtaPhiMVector& kminus2) {
        const double dRpp =
          deltaR(kplus1.Phi(), kplus1.Eta(),
                 kplus2.Phi(), kplus2.Eta());

        const double dRmm =
          deltaR(kminus1.Phi(), kminus1.Eta(),
                 kminus2.Phi(), kminus2.Eta());

        const double dRpm =
          deltaR(kplus1.Phi(), kplus1.Eta(),
                 kminus2.Phi(), kminus2.Eta());

        const double dRmp =
          deltaR(kminus1.Phi(), kminus1.Eta(),
                 kplus2.Phi(), kplus2.Eta());

        histos.fill(HIST("hDeltaRkaonplus"), dRpp);
        histos.fill(HIST("hDeltaRkaonminus"), dRmm);
        return std::min({dRpp, dRmm, dRpm, dRmp});
      };

    const auto phiPtAsymmetry =
      [](const auto& phi1,
         const auto& phi2) -> double {
      const double sumPt = phi1.Pt() + phi2.Pt();

      if (sumPt <= 0.) {
        return -1.;
      }

      return std::abs(phi1.Pt() - phi2.Pt()) / sumPt;
    };

    const auto nKaonTOFHits =
      [](const auto& t1, const auto& t2) -> int {
      return static_cast<int>(t1.phid1TOFHit() == 1) +
             static_cast<int>(t1.phid2TOFHit() == 1) +
             static_cast<int>(t2.phid1TOFHit() == 1) +
             static_cast<int>(t2.phid2TOFHit() == 1);
    };

    const auto doublePhiPIDScore =
      [](const auto& t1, const auto& t2) -> double {
      const auto kaonPIDScore =
        [](double nSigmaTPC,
           double nSigmaTOF,
           int tofHit) -> double {
        // TOFHit is 1 when TOF is available and -1 otherwise.
        if (tofHit == 1) {
          return std::sqrt(
            (nSigmaTPC * nSigmaTPC +
             nSigmaTOF * nSigmaTOF));
        }

        // TPC-only track
        return std::abs(nSigmaTPC);
      };

      const double pid1 =
        kaonPIDScore(t1.phid1TPC(),
                     t1.phid1TOF(),
                     t1.phid1TOFHit());

      const double pid2 =
        kaonPIDScore(t1.phid2TPC(),
                     t1.phid2TOF(),
                     t1.phid2TOFHit());

      const double pid3 =
        kaonPIDScore(t2.phid1TPC(),
                     t2.phid1TOF(),
                     t2.phid1TOFHit());

      const double pid4 =
        kaonPIDScore(t2.phid2TPC(),
                     t2.phid2TOF(),
                     t2.phid2TOFHit());

      return std::max({pid1, pid2, pid3, pid4});
    };

    std::vector<ROOT::Math::PtEtaPhiMVector> pairV;
    std::vector<ROOT::Math::PtEtaPhiMVector> phi1V;
    std::vector<ROOT::Math::PtEtaPhiMVector> phi2V;
    std::vector<double> minDRV;
    std::vector<int> nTOFV;
    std::vector<double> pid4KV;

    for (auto const& t1 : phitracks) {
      const double kplus1pt = std::hypot(t1.phid1Px(), t1.phid1Py());
      const double kminus1pt = std::hypot(t1.phid2Px(), t1.phid2Py());

      if (kplus1pt > config.maxKaonPt || kminus1pt > config.maxKaonPt) {
        continue;
      }
      if (!selectionPID(t1.phid1TPC(), t1.phid1TOF(), t1.phid1TOFHit(), config.strategyPID1, kplus1pt)) {
        continue;
      }
      if (!selectionPID(t1.phid2TPC(), t1.phid2TOF(), t1.phid2TOFHit(), config.strategyPID2, kminus1pt)) {
        continue;
      }

      TLorentzVector phi1;
      TLorentzVector k1p;
      TLorentzVector k1m;

      buildPhiAndKaons(t1, phi1, k1p, k1m);
      const double phi1Mass = cfgApplyKaonMomentumCorrection ? phi1.M() : t1.phiMass();

      if (phi1Mass < config.minPhiMass1 || phi1Mass > config.maxPhiMass1) {
        continue;
      }
      if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt) {
        continue;
      }

      const auto id1 = t1.index();

      for (auto const& t2 : phitracks) {
        const auto id2 = t2.index();
        if (id2 <= id1) {
          // LOGF(info, "track reject %d %d %f %f", id1, id2, t1.phiMass(), t2.phiMass());
          continue;
        }
        const double kplus2pt = std::hypot(t2.phid1Px(), t2.phid1Py());
        const double kminus2pt = std::hypot(t2.phid2Px(), t2.phid2Py());

        if (kplus2pt > config.maxKaonPt || kminus2pt > config.maxKaonPt) {
          continue;
        }
        if (!selectionPID(t2.phid1TPC(), t2.phid1TOF(), t2.phid1TOFHit(), config.strategyPID1, kplus2pt)) {
          continue;
        }
        if (!selectionPID(t2.phid2TPC(), t2.phid2TOF(), t2.phid2TOFHit(), config.strategyPID2, kminus2pt)) {
          continue;
        }

        TLorentzVector phi2;
        TLorentzVector k2p;
        TLorentzVector k2m;

        buildPhiAndKaons(t2, phi2, k2p, k2m);
        const double phi2Mass = cfgApplyKaonMomentumCorrection ? phi2.M() : t2.phiMass();

        if (phi2Mass < config.minPhiMass1 || phi2Mass > config.maxPhiMass1) {
          continue;
        }
        if (phi2.Pt() < config.minPhiPt || phi2.Pt() > config.maxPhiPt) {
          continue;
        }

        TLorentzVector pair = phi1 + phi2;
        if (pair.Pt() < config.minExoticPt) {
          continue;
        }
        if (pair.M() < config.minExoticMass || pair.M() > config.maxExoticMass) {
          continue;
        }
        // reject any shared daughter between the two phi candidates
        if (t1.phid1Index() == t2.phid1Index() ||
            t1.phid1Index() == t2.phid2Index() ||
            t1.phid2Index() == t2.phid1Index() ||
            t1.phid2Index() == t2.phid2Index()) {
          LOGF(info, "track share %d %d %d %d", t1.phid1Index(), t1.phid2Index(), t2.phid1Index(), t2.phid2Index());
          continue;
        }

        auto cross12 = k1p + k2m;
        auto cross21 = k2p + k1m;
        bool alternativePairValid = cross12.M() > config.cfgCrossPhiLow && cross12.M() < config.cfgCrossPhiHigh && cross21.M() > config.cfgCrossPhiLow && cross21.M() < config.cfgCrossPhiHigh;
        if (alternativePairValid) {
          float scoreOriginal = deltaMPhiNominal(phi1.M(), phi2.M());
          float scoreCross = deltaMPhiNominal(cross12.M(), cross21.M());
          if (config.useCrossPairRejection && (scoreCross < scoreOriginal)) {
            LOGF(info, "Best-pairing rejected: original score = %3.4f, cross scoremPhi2 = %3.4f", scoreOriginal, scoreCross);
            continue; // another pairing of these four tracks is better
          }
        }

        // Compact candidate-level input for the data-driven X resolution.
        // Fill only in the requested X mass window and above the X pT threshold.
        // At this point the candidate has already passed:
        //   - daughter pT/PID selections,
        //   - phi mass and phi pT selections,
        //   - pair pT / broad pair-mass selections,
        //   - shared-track rejection,
        //   - cross-pairing rejection.
        //
        // Use the same four-vectors that define the selected candidate.  Thus,
        // if cfgApplyKaonMomentumCorrection is enabled, the sparse consistently
        // stores the corrected candidate kinematics; otherwise it stores raw
        // reconstructed kinematics.
        if (cfgFillSelectedXResolutionSparse &&
            pair.Pt() > cfgSelectedXResolutionPtMin &&
            pair.M() > cfgSelectedXResolutionMassLow &&
            pair.M() < cfgSelectedXResolutionMassHigh) {
          histos.fill(HIST("SelectedXResolutionDataDriven"),
                      pair.Pt(),
                      phi1.M(), phi2.M(),
                      phi1.Pt(), phi2.Pt(),
                      k1p.Pt(), k1m.Pt(),
                      k2p.Pt(), k2m.Pt());
        }

        histos.fill(HIST("hPhiMass"), phi1.M(), phi2.M(), pair.Pt());
        histos.fill(HIST("hPhiMassNormalized"), getNormalizedMPhi(phi1.M(), phi1.Pt()), getNormalizedMPhi(phi2.M(), phi2.Pt()), pair.Pt());

        ROOT::Math::PtEtaPhiMVector k1pV(k1p.Pt(), k1p.Eta(), k1p.Phi(), mKPDG);
        ROOT::Math::PtEtaPhiMVector k1mV(k1m.Pt(), k1m.Eta(), k1m.Phi(), mKPDG);
        ROOT::Math::PtEtaPhiMVector k2pV(k2p.Pt(), k2p.Eta(), k2p.Phi(), mKPDG);
        ROOT::Math::PtEtaPhiMVector k2mV(k2m.Pt(), k2m.Eta(), k2m.Phi(), mKPDG);

        const double minDR = minKaonDeltaR(k1pV, k1mV, k2pV, k2mV);
        const int nTOF = nKaonTOFHits(t1, t2);
        const double pid4K = doublePhiPIDScore(t1, t2);
        pairV.emplace_back(pair.Pt(), pair.Eta(), pair.Phi(), pair.M());
        phi1V.emplace_back(phi1.Pt(), phi1.Eta(), phi1.Phi(), phi1.M());
        phi2V.emplace_back(phi2.Pt(), phi2.Eta(), phi2.Phi(), phi2.M());
        minDRV.emplace_back(minDR);
        nTOFV.emplace_back(nTOF);
        pid4KV.emplace_back(pid4K);
      }
    }

    if (pairV.empty()) {
      return;
    }

    for (size_t i = 0; i < pairV.size(); ++i) {
      TLorentzVector p1;
      TLorentzVector p2;
      TLorentzVector pair;

      p1.SetPtEtaPhiM(phi1V[i].Pt(), phi1V[i].Eta(), phi1V[i].Phi(), phi1V[i].M());
      p2.SetPtEtaPhiM(phi2V[i].Pt(), phi2V[i].Eta(), phi2V[i].Phi(), phi2V[i].M());
      pair.SetPtEtaPhiM(pairV[i].Pt(), pairV[i].Eta(), pairV[i].Phi(), pairV[i].M());

      const double M = pair.M();
      const double pairPt = pair.Pt();
      const double dRphi = deltaR(p1.Phi(), p1.Eta(), p2.Phi(), p2.Eta());
      const double minDR = minDRV[i];
      const double combine4kpid = pid4KV[i];
      const double nkaonTOF = nTOFV[i];
      if (!config.useParametrized) {
        dMNominal = deltaMPhiNominal(p1.M(), p2.M());
      } else {
        dMNominal = getDeltaMPhi(p1.M(), p1.Pt(), p2.M(), p2.Pt());
      }
      const double dMNominalNsigma = getNormalizedDeltaMPhi(p1.M(), p1.Pt(), p2.M(), p2.Pt());
      const double denom = std::abs(pairPt - p1.Pt());

      const double ptcorr = p1.Pt() / denom;

      const double pt1 = p1.Pt();
      const double pt2 = p2.Pt();
      const double ptsum = pt1 + pt2;
      if (ptsum <= 0.0) {
        continue;
      }
      const double apt = phiPtAsymmetry(p1, p2);
      // const double absCosTheta = absCosThetaStar(p1, p2);

      if (pairPt > config.minExoticPt) {
        histos.fill(HIST("hPtCorrelation"), pairPt, ptcorr);
        // histos.fill(HIST("hMassCent"), p1.M(), p2.M(), collision.centrality());
        histos.fill(HIST("SEMassUnlike_AllVars"),
                    M,
                    pairPt,
                    dRphi,
                    minDR,
                    p1.M(),
                    p2.M(),
                    dMNominal,
                    dMNominalNsigma,
                    ptcorr,
                    pairV.size());

        histos.fill(HIST("SEMassDoublePhi"),
                    M,
                    pairPt,
                    apt,
                    std::abs(pair.Rapidity()),
                    p1.M(),
                    p2.M(),
                    dMNominal,
                    nkaonTOF,
                    combine4kpid);
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processopti5, "Process Optimized same event with all variables", true);

  void processPairOpti6(aod::PhiPhiPairs const& pairs)
  {
    constexpr double mPhiPDG = o2::constants::physics::MassPhi;
    constexpr double mKPDG = o2::constants::physics::MassKPlus;

    for (auto const& pair : pairs) {
      if (pair.k1Index() == pair.k2Index() || pair.k1Index() == pair.k3Index() || pair.k1Index() == pair.k4Index() || pair.k2Index() == pair.k3Index() || pair.k2Index() == pair.k4Index() || pair.k3Index() == pair.k4Index()) {
        continue;
      }
      if (std::abs(pair.k1DcaXY()) >= 0.02f || std::abs(pair.k1DcaZ()) >= 0.02f || std::abs(pair.k2DcaXY()) >= 0.02f || std::abs(pair.k2DcaZ()) >= 0.02f || std::abs(pair.k3DcaXY()) >= 0.02f || std::abs(pair.k3DcaZ()) >= 0.02f || std::abs(pair.k4DcaXY()) >= 0.02f || std::abs(pair.k4DcaZ()) >= 0.02f) {
        continue;
      }
      const double pairPt = std::hypot(pair.pairPx(), pair.pairPy());
      if (pairPt <= config.minExoticPt || pair.pairMass() < config.minExoticMass || pair.pairMass() > config.maxExoticMass) {
        continue;
      }

      const double pt1 = std::hypot(pair.phi1Px(), pair.phi1Py());
      const double pt2 = std::hypot(pair.phi2Px(), pair.phi2Py());
      if (pt1 < config.minPhiPt || pt1 > config.maxPhiPt || pt2 < config.minPhiPt || pt2 > config.maxPhiPt) {
        continue;
      }

      const ROOT::Math::PxPyPzMVector k1(pair.k1Px(), pair.k1Py(), pair.k1Pz(), mKPDG);
      const ROOT::Math::PxPyPzMVector k2(pair.k2Px(), pair.k2Py(), pair.k2Pz(), mKPDG);
      const ROOT::Math::PxPyPzMVector k3(pair.k3Px(), pair.k3Py(), pair.k3Pz(), mKPDG);
      const ROOT::Math::PxPyPzMVector k4(pair.k4Px(), pair.k4Py(), pair.k4Pz(), mKPDG);

      const double mCross14 = (k1 + k4).M();
      const double mCross32 = (k3 + k2).M();
      if ((mCross14 > config.cfgCrossPhiLow && mCross14 < config.cfgCrossPhiHigh) || (mCross32 > config.cfgCrossPhiLow && mCross32 < config.cfgCrossPhiHigh)) {
        continue;
      }

      double deltaM = std::hypot(pair.phi1Mass() - mPhiPDG, pair.phi2Mass() - mPhiPDG);
      if (config.useParametrized) {
        deltaM = getDeltaMPhi(pair.phi1Mass(), pt1, pair.phi2Mass(), pt2);
      }

      const double decayLength = pair.vertexL3DSig();
      const double fitChi2Ndf = pair.fitChi2Ndf();
      const double rmsDcaSig = pair.rmsDcaSig();

      if (!std::isfinite(pair.pairMass()) || !std::isfinite(deltaM) || !std::isfinite(pt1) || !std::isfinite(pt2) || !std::isfinite(decayLength) || !std::isfinite(fitChi2Ndf) || !std::isfinite(rmsDcaSig)) {
        continue;
      }

      histos.fill(HIST("SEMassUnlike_VertexVars"), pair.pairMass(), pairPt, deltaM, pair.phi1Mass(), pair.phi2Mass(), decayLength, fitChi2Ndf, rmsDcaSig);
    }
  }
  PROCESS_SWITCH(doublephimeson, processPairOpti6, "Process fitted phi-phi pairs with vertex variables", false);

  void processOpti7(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (config.additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    constexpr double mPhiPDG = o2::constants::physics::MassPhi;
    // constexpr double mKPDG = o2::constants::physics::MassKPlus;
    int phimult = 0;

    for (const auto& Phitrack : phitracks) {
      const double kpluspt = std::hypot(Phitrack.phid1Px(), Phitrack.phid1Py());
      const double kminuspt = std::hypot(Phitrack.phid2Px(), Phitrack.phid2Py());

      // pT cut on kaon tracks
      if (kpluspt > config.maxKaonPt || kminuspt > config.maxKaonPt) {
        continue;
      }

      // Mass window selection for phi candidates
      if (Phitrack.phiMass() < config.minPhiMass1 || Phitrack.phiMass() > config.maxPhiMass1) {
        continue;
      }

      TLorentzVector phi;
      phi.SetXYZM(Phitrack.phiPx(), Phitrack.phiPy(), Phitrack.phiPz(), Phitrack.phiMass());

      // pT cut on phi daughters
      if (phi.Pt() < config.minPhiPt || phi.Pt() > config.maxPhiPt) {
        continue;
      }

      // PID selection for kaon tracks
      if (!selectionPID(Phitrack.phid1TPC(), Phitrack.phid1TOF(), Phitrack.phid1TOFHit(), config.strategyPID1, kpluspt)) {
        continue;
      }
      if (!selectionPID(Phitrack.phid2TPC(), Phitrack.phid2TOF(), Phitrack.phid2TOFHit(), config.strategyPID2, kminuspt)) {
        continue;
      }
      phimult++;
    }

    if (phimult < 2) {
      return;
    }

    for (auto const& Phitrack1 : phitracks) {
      for (auto const& Phitrack2 : phitracks) {

        // Avoid double counting
        if (Phitrack2.index() <= Phitrack1.index()) {
          continue;
        }

        // pT cut kaon tracks
        const double kplus1pt = std::hypot(Phitrack1.phid1Px(), Phitrack1.phid1Py());
        const double kminus1pt = std::hypot(Phitrack1.phid2Px(), Phitrack1.phid2Py());
        const double kplus2pt = std::hypot(Phitrack2.phid1Px(), Phitrack2.phid1Py());
        const double kminus2pt = std::hypot(Phitrack2.phid2Px(), Phitrack2.phid2Py());

        if (kplus1pt > config.maxKaonPt || kminus1pt > config.maxKaonPt || kplus2pt > config.maxKaonPt || kminus2pt > config.maxKaonPt) {
          continue;
        }

        // Mass window selection for phi candidates
        if (Phitrack1.phiMass() < config.minPhiMass1 || Phitrack1.phiMass() > config.maxPhiMass1 || Phitrack2.phiMass() < config.minPhiMass1 || Phitrack2.phiMass() > config.maxPhiMass1) {
          continue;
        }

        // pT cut on phi daughters
        TLorentzVector phi1, phi2;
        phi1.SetXYZM(Phitrack1.phiPx(), Phitrack1.phiPy(), Phitrack1.phiPz(), Phitrack1.phiMass());
        phi2.SetXYZM(Phitrack2.phiPx(), Phitrack2.phiPy(), Phitrack2.phiPz(), Phitrack2.phiMass());
        if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt || phi2.Pt() < config.minPhiPt || phi2.Pt() > config.maxPhiPt) {
          continue;
        }

        // PID selection for kaon tracks
        if (!selectionPID(Phitrack1.phid1TPC(), Phitrack1.phid1TOF(), Phitrack1.phid1TOFHit(), config.strategyPID1, kplus1pt) ||
            !selectionPID(Phitrack1.phid2TPC(), Phitrack1.phid2TOF(), Phitrack1.phid2TOFHit(), config.strategyPID2, kminus1pt) ||
            !selectionPID(Phitrack2.phid1TPC(), Phitrack2.phid1TOF(), Phitrack2.phid1TOFHit(), config.strategyPID1, kplus2pt) ||
            !selectionPID(Phitrack2.phid2TPC(), Phitrack2.phid2TOF(), Phitrack2.phid2TOFHit(), config.strategyPID2, kminus2pt)) {
          continue;
        }

        // Check for shared daughters
        if (Phitrack1.phid1Index() == Phitrack2.phid1Index() ||
            Phitrack1.phid1Index() == Phitrack2.phid2Index() ||
            Phitrack1.phid2Index() == Phitrack2.phid1Index() ||
            Phitrack1.phid2Index() == Phitrack2.phid2Index()) {
          continue;
        }

        TLorentzVector pair = phi1 + phi2;
        // Mass window range for the phi-phi pair
        if (pair.Pt() < config.minExoticPt || pair.M() < config.minExoticMass || pair.M() > config.maxExoticMass) {
          continue;
        }

        double deltaM = std::hypot(Phitrack1.phiMass() - mPhiPDG, Phitrack2.phiMass() - mPhiPDG);

        histos.fill(HIST("SEMassPhiPhi"),
                    pair.M(),
                    pair.Pt(),
                    deltaM);
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processOpti7, "Process optimised save-event for cross-checks", false);

  void processOpti8(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    histos.fill(HIST("NEvents"), 0.5);
    if (config.additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    histos.fill(HIST("NEvents"), 1.5);
    constexpr double mPhiPDG = o2::constants::physics::MassPhi;
    constexpr double mKPDG = o2::constants::physics::MassKPlus;
    int phimult = 0;

    for (const auto& Phitrack : phitracks) {
      const double kpluspt = std::hypot(Phitrack.phid1Px(), Phitrack.phid1Py());
      const double kminuspt = std::hypot(Phitrack.phid2Px(), Phitrack.phid2Py());

      histos.fill(HIST("hnsigmaTPCTOFKaonBefore"), Phitrack.phid1TPC(), Phitrack.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlusBefore"), Phitrack.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinusBefore"), Phitrack.phid2TPC(), kminuspt);

      // pT cut on kaon tracks
      if (kpluspt > config.maxKaonPt || kminuspt > config.maxKaonPt) {
        continue;
      }

      // Mass window selection for phi candidates
      if (Phitrack.phiMass() < config.minPhiMass1 || Phitrack.phiMass() > config.maxPhiMass1) {
        continue;
      }

      TLorentzVector phi;
      phi.SetXYZM(Phitrack.phiPx(), Phitrack.phiPy(), Phitrack.phiPz(), Phitrack.phiMass());

      // pT cut on phi daughters
      if (phi.Pt() < config.minPhiPt || phi.Pt() > config.maxPhiPt) {
        continue;
      }

      // PID selection for kaon tracks
      if (!selectionPID(Phitrack.phid1TPC(), Phitrack.phid1TOF(), Phitrack.phid1TOFHit(), config.strategyPID1, kpluspt)) {
        continue;
      }
      if (!selectionPID(Phitrack.phid2TPC(), Phitrack.phid2TOF(), Phitrack.phid2TOFHit(), config.strategyPID2, kminuspt)) {
        continue;
      }

      histos.fill(HIST("hnsigmaTPCTOFKaon"), Phitrack.phid1TPC(), Phitrack.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), Phitrack.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), Phitrack.phid2TPC(), kminuspt);

      phimult++;
    }

    histos.fill(HIST("NPhiPerEvent"), phimult);

    if (phimult < 2) {
      return;
    }

    for (auto const& Phitrack1 : phitracks) {

      // pT cut kaon tracks
      const double kplus1pt = std::hypot(Phitrack1.phid1Px(), Phitrack1.phid1Py());
      const double kminus1pt = std::hypot(Phitrack1.phid2Px(), Phitrack1.phid2Py());
      if (kplus1pt > config.maxKaonPt || kminus1pt > config.maxKaonPt) {
        continue;
      }

      // PID selection for kaon tracks
      if (!selectionPID(Phitrack1.phid1TPC(), Phitrack1.phid1TOF(), Phitrack1.phid1TOFHit(), config.strategyPID1, kplus1pt) ||
          !selectionPID(Phitrack1.phid2TPC(), Phitrack1.phid2TOF(), Phitrack1.phid2TOFHit(), config.strategyPID2, kminus1pt)) {
        continue;
      }

      TLorentzVector phi1;
      phi1.SetXYZM(Phitrack1.phiPx(), Phitrack1.phiPy(), Phitrack1.phiPz(), Phitrack1.phiMass());
      if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt) {
        continue;
      }

      // Mass window selection for phi candidates
      if (Phitrack1.phiMass() < config.minPhiMass1 || Phitrack1.phiMass() > config.maxPhiMass1) {
        continue;
      }

      histos.fill(HIST("hPhiMassVsPt"), phi1.M(), phi1.Pt());

      for (auto const& Phitrack2 : phitracks) {

        // Avoid double counting
        if (Phitrack2.index() <= Phitrack1.index()) {
          continue;
        }

        // pT cut kaon tracks
        const double kplus2pt = std::hypot(Phitrack2.phid1Px(), Phitrack2.phid1Py());
        const double kminus2pt = std::hypot(Phitrack2.phid2Px(), Phitrack2.phid2Py());

        if (kplus2pt > config.maxKaonPt || kminus2pt > config.maxKaonPt) {
          continue;
        }

        // Mass window selection for phi candidates
        if (Phitrack2.phiMass() < config.minPhiMass1 || Phitrack2.phiMass() > config.maxPhiMass1) {
          continue;
        }

        // pT cut on phi
        TLorentzVector phi2;
        phi2.SetXYZM(Phitrack2.phiPx(), Phitrack2.phiPy(), Phitrack2.phiPz(), Phitrack2.phiMass());
        if (phi2.Pt() < config.minPhiPt || phi2.Pt() > config.maxPhiPt) {
          continue;
        }

        // PID selection for kaon tracks
        if (!selectionPID(Phitrack2.phid1TPC(), Phitrack2.phid1TOF(), Phitrack2.phid1TOFHit(), config.strategyPID1, kplus2pt) ||
            !selectionPID(Phitrack2.phid2TPC(), Phitrack2.phid2TOF(), Phitrack2.phid2TOFHit(), config.strategyPID2, kminus2pt)) {
          continue;
        }

        // Check for shared daughters
        if (Phitrack1.phid1Index() == Phitrack2.phid1Index() ||
            Phitrack1.phid1Index() == Phitrack2.phid2Index() ||
            Phitrack1.phid2Index() == Phitrack2.phid1Index() ||
            Phitrack1.phid2Index() == Phitrack2.phid2Index()) {
          continue;
        }

        for (int i = 0; i < config.cRotations; i++) {
          double thetaRot = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / 10, o2::constants::math::PI + o2::constants::math::PI / 10);

          TLorentzVector daughterRot;
          daughterRot.SetXYZM(phi1.Px() * std::cos(thetaRot) - phi1.Py() * std::sin(thetaRot), phi1.Px() * std::sin(thetaRot) + phi1.Py() * std::cos(thetaRot), phi1.Pz(), phi1.M());

          TLorentzVector pairRot = daughterRot + phi2;

          if (pairRot.Pt() < config.minExoticPt || pairRot.M() < config.minExoticMass || pairRot.M() > config.maxExoticMass) {
            continue;
          }

          histos.fill(HIST("SEMassPhiPhiRotational"),
                      pairRot.M(),
                      pairRot.Pt());
        }

        TLorentzVector pair = phi1 + phi2;
        // Mass window range for the phi-phi pair
        if (pair.Pt() < config.minExoticPt || pair.M() < config.minExoticMass || pair.M() > config.maxExoticMass) {
          continue;
        }
        histos.fill(HIST("hPhiMass"), phi1.M(), phi2.M(), pair.Pt());

        // =====================================================================
        //                  4-KAON KINEMATIC FIT
        //                  2 x M(KK) = M(phi) constraints
        // =====================================================================

        TLorentzVector k11, k12, k21, k22;
        k11.SetXYZM(Phitrack1.phid1Px(), Phitrack1.phid1Py(), Phitrack1.phid1Pz(), mKPDG);
        k12.SetXYZM(Phitrack1.phid2Px(), Phitrack1.phid2Py(), Phitrack1.phid2Pz(), mKPDG);
        k21.SetXYZM(Phitrack2.phid1Px(), Phitrack2.phid1Py(), Phitrack2.phid1Pz(), mKPDG);
        k22.SetXYZM(Phitrack2.phid2Px(), Phitrack2.phid2Py(), Phitrack2.phid2Pz(), mKPDG);

        FourKFitResult fitResult = fitFourKaons(k11, k12, k21, k22, mKPDG, mPhiPDG, 0.01, 20, 1e-6);

        // double refittedMass = pair.M();
        double fitChi2 = -1.0;
        double fitProb = -1.0;

        if (fitResult.converged) {
          // refittedMass = fitResult.refittedMass;
          fitChi2 = fitResult.chi2;
          fitProb = fitResult.probability;
        }

        double deltaM = std::hypot(Phitrack1.phiMass() - mPhiPDG, Phitrack2.phiMass() - mPhiPDG);

        if (fitResult.converged) {
          histos.fill(HIST("SEMassPhiPhiRefitted"),
                      pair.M(),
                      pair.Pt(),
                      deltaM,
                      fitChi2,
                      fitProb,
                      phi1.M(),
                      phi2.M());
        }
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processOpti8, "Process optimised save-event with phi-phi pairs with 4K fit", false);

  double getMomentumCorrection(double phiPt)
  {
    constexpr int kNBins = 14;
    constexpr double kPtEdges[kNBins + 1] = {0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    constexpr double kEpsilon[kNBins] = {0.0024146539379, 0.00762912304405, 0.00817386742107, 0.00982766792836, 0.0106788436324, 0.0112954252224, 0.0115583787246, 0.0117416086153, 0.012422150094, 0.011968584565, 0.0125561781092, 0.0139372556405, 0.014995530318, 0.0165687496159};

    if (phiPt < kPtEdges[0] || phiPt >= kPtEdges[kNBins])
      return 0.0;

    auto it = std::upper_bound(std::begin(kPtEdges), std::end(kPtEdges), phiPt);
    int binIndex = std::distance(std::begin(kPtEdges), it) - 1;

    return kEpsilon[binIndex];
  }

  TLorentzVector CorrectKaonMomentum(const TLorentzVector& kaon, double epsilon)
  {
    const double scale = 1.0 + epsilon;
    constexpr double mKPDG = o2::constants::physics::MassKPlus;
    TLorentzVector corrected;

    // Scale the 3-momentum
    TVector3 pCorr = scale * kaon.Vect();

    // Recalculate energy using fixed kaon mass
    corrected.SetVectM(pCorr, mKPDG);

    return corrected;
  }

  void processOpti9(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (config.additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    constexpr double mPhiPDG = o2::constants::physics::MassPhi;
    constexpr double mKPDG = o2::constants::physics::MassKPlus;
    int phimult = 0;

    for (const auto& Phitrack : phitracks) {
      const double kpluspt = std::hypot(Phitrack.phid1Px(), Phitrack.phid1Py());
      const double kminuspt = std::hypot(Phitrack.phid2Px(), Phitrack.phid2Py());

      histos.fill(HIST("hnsigmaTPCTOFKaonBefore"), Phitrack.phid1TPC(), Phitrack.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlusBefore"), Phitrack.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinusBefore"), Phitrack.phid2TPC(), kminuspt);

      // pT cut on kaon tracks
      if (kpluspt > config.maxKaonPt || kminuspt > config.maxKaonPt) {
        continue;
      }

      // Mass window selection for phi candidates
      if (Phitrack.phiMass() < config.minPhiMass1 || Phitrack.phiMass() > config.maxPhiMass1) {
        continue;
      }

      TLorentzVector phi;
      phi.SetXYZM(Phitrack.phiPx(), Phitrack.phiPy(), Phitrack.phiPz(), Phitrack.phiMass());

      if (phi.Pt() < config.minPhiPt || phi.Pt() > config.maxPhiPt) {
        continue;
      }

      // PID selection for kaon tracks
      if (!selectionPID(Phitrack.phid1TPC(), Phitrack.phid1TOF(), Phitrack.phid1TOFHit(), config.strategyPID1, kpluspt)) {
        continue;
      }
      if (!selectionPID(Phitrack.phid2TPC(), Phitrack.phid2TOF(), Phitrack.phid2TOFHit(), config.strategyPID2, kminuspt)) {
        continue;
      }

      histos.fill(HIST("hnsigmaTPCTOFKaon"), Phitrack.phid1TPC(), Phitrack.phid1TOF(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), Phitrack.phid1TPC(), kpluspt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), Phitrack.phid2TPC(), kminuspt);

      phimult++;
    }

    if (phimult < 2) {
      return;
    }

    for (auto const& Phitrack1 : phitracks) {
      // pT cut kaon tracks
      const double kplus1pt = std::hypot(Phitrack1.phid1Px(), Phitrack1.phid1Py());
      const double kminus1pt = std::hypot(Phitrack1.phid2Px(), Phitrack1.phid2Py());
      if (kplus1pt > config.maxKaonPt || kminus1pt > config.maxKaonPt) {
        continue;
      }

      // PID selection for kaon tracks
      if (!selectionPID(Phitrack1.phid1TPC(), Phitrack1.phid1TOF(), Phitrack1.phid1TOFHit(), config.strategyPID1, kplus1pt) ||
          !selectionPID(Phitrack1.phid2TPC(), Phitrack1.phid2TOF(), Phitrack1.phid2TOFHit(), config.strategyPID2, kminus1pt)) {
        continue;
      }

      // pT cut on Phi resonance
      TLorentzVector phi1;
      phi1.SetXYZM(Phitrack1.phiPx(), Phitrack1.phiPy(), Phitrack1.phiPz(), Phitrack1.phiMass());
      if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt) {
        continue;
      }

      // Mass window selection for phi candidates
      if (Phitrack1.phiMass() < config.minPhiMass1 || Phitrack1.phiMass() > config.maxPhiMass1) {
        continue;
      }
      histos.fill(HIST("hPhiMassVsPt"), phi1.M(), phi1.Pt());

      // Uncorrected phi pT for scale factor lookup
      const double phi1UncorrPt = std::hypot(Phitrack1.phiPx(), Phitrack1.phiPy());

      // Correction based on ORIGINAL phi pT
      const double epsilon1 = getMomentumCorrection(phi1UncorrPt);

      // Original kaons
      TLorentzVector kplus1, kminus1;
      kplus1.SetXYZM(Phitrack1.phid1Px(), Phitrack1.phid1Py(), Phitrack1.phid1Pz(), mKPDG);
      kminus1.SetXYZM(Phitrack1.phid2Px(), Phitrack1.phid2Py(), Phitrack1.phid2Pz(), mKPDG);

      // Corrected kaon momenta
      TLorentzVector kplus1Corr = CorrectKaonMomentum(kplus1, epsilon1);
      TLorentzVector kminus1Corr = CorrectKaonMomentum(kminus1, epsilon1);
      TLorentzVector phi1Corr = kplus1Corr + kminus1Corr;

      histos.fill(HIST("hPhiMassVsPtShifted"), phi1Corr.M(), phi1Corr.Pt());

      for (auto const& Phitrack2 : phitracks) {

        if (Phitrack2.index() <= Phitrack1.index()) {
          continue;
        }

        // pT cut kaon tracks
        const double kplus2pt = std::hypot(Phitrack2.phid1Px(), Phitrack2.phid1Py());
        const double kminus2pt = std::hypot(Phitrack2.phid2Px(), Phitrack2.phid2Py());

        if (kplus2pt > config.maxKaonPt || kminus2pt > config.maxKaonPt) {
          continue;
        }

        // Mass window selection for phi candidates
        if (Phitrack2.phiMass() < config.minPhiMass1 || Phitrack2.phiMass() > config.maxPhiMass1) {
          continue;
        }

        // Uncorrected phi pT for scale factor lookup
        const double phi2UncorrPt = std::hypot(Phitrack2.phiPx(), Phitrack2.phiPy());

        TLorentzVector phi2;
        phi2.SetXYZM(Phitrack2.phiPx(), Phitrack2.phiPy(), Phitrack2.phiPz(), Phitrack2.phiMass());
        if (phi2.Pt() < config.minPhiPt || phi2.Pt() > config.maxPhiPt) {
          continue;
        }

        // PID selection for kaon original tracks
        if (!selectionPID(Phitrack2.phid1TPC(), Phitrack2.phid1TOF(), Phitrack2.phid1TOFHit(), config.strategyPID1, kplus2pt) || !selectionPID(Phitrack2.phid2TPC(), Phitrack2.phid2TOF(), Phitrack2.phid2TOFHit(), config.strategyPID2, kminus2pt)) {
          continue;
        }

        // Check for shared daughters
        if (Phitrack1.phid1Index() == Phitrack2.phid1Index() ||
            Phitrack1.phid1Index() == Phitrack2.phid2Index() ||
            Phitrack1.phid2Index() == Phitrack2.phid1Index() ||
            Phitrack1.phid2Index() == Phitrack2.phid2Index()) {
          continue;
        }

        // Correction based on ORIGINAL phi pT
        const double epsilon2 = getMomentumCorrection(phi2UncorrPt);

        // Original kaons
        TLorentzVector kplus2, kminus2;

        kplus2.SetXYZM(Phitrack2.phid1Px(), Phitrack2.phid1Py(), Phitrack2.phid1Pz(), mKPDG);
        kminus2.SetXYZM(Phitrack2.phid2Px(), Phitrack2.phid2Py(), Phitrack2.phid2Pz(), mKPDG);

        // Correct kaon momenta
        TLorentzVector kplus2Corr = CorrectKaonMomentum(kplus2, epsilon2);
        TLorentzVector kminus2Corr = CorrectKaonMomentum(kminus2, epsilon2);

        // double kplus2CorrPt = kplus2Corr.Pt();
        // double kminus2CorrPt = kminus2Corr.Pt();

        // Corrected phi
        TLorentzVector phi2Corr = kplus2Corr + kminus2Corr;

        // Reconstruct double-phi pair from SHIFTED phi candidates
        TLorentzVector pair = phi1 + phi2;
        TLorentzVector pairShifted = phi1Corr + phi2Corr;

        if (pair.Pt() < config.minExoticPt || pair.M() < config.minExoticMass || pair.M() > config.maxExoticMass) {
          continue;
        }

        histos.fill(HIST("hPhiMass"), phi1.M(), phi2.M(), pair.Pt());
        histos.fill(HIST("hPhiMassShifted"), phi1Corr.M(), phi2Corr.M(), pairShifted.Pt());

        double deltaMShifted = std::hypot(phi1Corr.M() - mPhiPDG, phi2Corr.M() - mPhiPDG);
        // double deltaM = std::hypot(phi1.M() - mPhiPDG, phi2.M() - mPhiPDG);

        // 4-Kaon Kinematic Fit using SHIFTED kaon vectors
        FourKFitResult fitResult = fitFourKaons(kplus1Corr, kminus1Corr, kplus2Corr, kminus2Corr, mKPDG, mPhiPDG, 0.01, 20, 1e-6);

        // double refittedMass = pairShifted.M();
        double fitChi2 = -1.0;
        double fitProb = -1.0;

        if (fitResult.converged) {
          // refittedMass = fitResult.refittedMass;
          fitChi2 = fitResult.chi2;
          fitProb = fitResult.probability;
        }

        // Fill shifted double-phi THnSparse
        histos.fill(HIST("SEMassPhiPhiShifted"),
                    pairShifted.M(),
                    pairShifted.Pt(),
                    deltaMShifted,
                    fitChi2,
                    fitProb,
                    phi1Corr.M(),
                    phi2Corr.M());
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processOpti9, "Process optimised save-event for phi-phi pairs after kaon momentum shift", false);

  SliceCache cache;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;

  void processMixedEventopti5(aod::RedPhiEvents& collisions, aod::PhiTracks& phitracks)
  {
    auto tracksTuple = std::make_tuple(phitracks);

    BinningTypeVertexContributor binningOnPositions{{CfgVtxBins, CfgMultBins}, true};

    SameKindPair<aod::RedPhiEvents, aod::PhiTracks, BinningTypeVertexContributor> pair{
      binningOnPositions, nEvtMixing, -1, collisions, tracksTuple, &cache};

    constexpr double mKPDG = o2::constants::physics::MassKPlus; // GeV/c^2

    const auto deltaR = [](double phi1, double eta1, double phi2, double eta2) {
      const double dphi = std::abs(TVector2::Phi_mpi_pi(phi1 - phi2));
      const double deta = eta1 - eta2;
      return std::sqrt(dphi * dphi + deta * deta);
    };

    const auto minKaonDeltaR =
      [&](const ROOT::Math::PtEtaPhiMVector& kplus1,
          const ROOT::Math::PtEtaPhiMVector& kminus1,
          const ROOT::Math::PtEtaPhiMVector& kplus2,
          const ROOT::Math::PtEtaPhiMVector& kminus2) {
        const double dRkplus =
          deltaR(kplus1.Phi(), kplus1.Eta(), kplus2.Phi(), kplus2.Eta());

        const double dRkminus =
          deltaR(kminus1.Phi(), kminus1.Eta(), kminus2.Phi(), kminus2.Eta());

        histos.fill(HIST("hDeltaRkaonplus"), dRkplus);
        histos.fill(HIST("hDeltaRkaonminus"), dRkminus);

        double minDR = dRkplus;
        minDR = std::min(minDR, dRkminus);

        minDR = std::min(minDR,
                         deltaR(kplus1.Phi(), kplus1.Eta(),
                                kminus1.Phi(), kminus1.Eta()));

        minDR = std::min(minDR,
                         deltaR(kplus1.Phi(), kplus1.Eta(),
                                kminus2.Phi(), kminus2.Eta()));

        minDR = std::min(minDR,
                         deltaR(kplus2.Phi(), kplus2.Eta(),
                                kminus1.Phi(), kminus1.Eta()));

        minDR = std::min(minDR,
                         deltaR(kplus2.Phi(), kplus2.Eta(),
                                kminus2.Phi(), kminus2.Eta()));

        return minDR;
      };

    struct PhiCand {
      ROOT::Math::PtEtaPhiMVector phi;
      ROOT::Math::PtEtaPhiMVector kplus;
      ROOT::Math::PtEtaPhiMVector kminus;
    };

    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if (collision1.index() == collision2.index()) {
        continue;
      }

      if (config.additionalEvsel) {
        if (collision1.numPos() < 2 || collision1.numNeg() < 2) {
          continue;
        }
        if (collision2.numPos() < 2 || collision2.numNeg() < 2) {
          continue;
        }
      }

      std::vector<PhiCand> cands1;
      std::vector<PhiCand> cands2;

      // ======================================================
      // Build phi candidates from event 1
      // ======================================================
      for (auto const& t1 : tracks1) {
        const double kplus1pt = std::hypot(t1.phid1Px(), t1.phid1Py());
        const double kminus1pt = std::hypot(t1.phid2Px(), t1.phid2Py());

        if (kplus1pt > config.maxKaonPt || kminus1pt > config.maxKaonPt) {
          continue;
        }

        if (!selectionPID(t1.phid1TPC(), t1.phid1TOF(), t1.phid1TOFHit(),
                          config.strategyPID1, kplus1pt)) {
          continue;
        }

        if (!selectionPID(t1.phid2TPC(), t1.phid2TOF(), t1.phid2TOFHit(),
                          config.strategyPID2, kminus1pt)) {
          continue;
        }

        if (t1.phiMass() < config.minPhiMass1 || t1.phiMass() > config.maxPhiMass1) {
          continue;
        }

        TLorentzVector phi1;
        TLorentzVector k1p;
        TLorentzVector k1m;

        phi1.SetXYZM(t1.phiPx(), t1.phiPy(), t1.phiPz(), t1.phiMass());
        k1p.SetXYZM(t1.phid1Px(), t1.phid1Py(), t1.phid1Pz(), mKPDG);
        k1m.SetXYZM(t1.phid2Px(), t1.phid2Py(), t1.phid2Pz(), mKPDG);

        if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt) {
          continue;
        }

        PhiCand cand;
        cand.phi = ROOT::Math::PtEtaPhiMVector(phi1.Pt(), phi1.Eta(), phi1.Phi(), phi1.M());
        cand.kplus = ROOT::Math::PtEtaPhiMVector(k1p.Pt(), k1p.Eta(), k1p.Phi(), mKPDG);
        cand.kminus = ROOT::Math::PtEtaPhiMVector(k1m.Pt(), k1m.Eta(), k1m.Phi(), mKPDG);

        cands1.emplace_back(std::move(cand));
      }

      // ======================================================
      // Build phi candidates from event 2
      // ======================================================
      for (auto const& t2 : tracks2) {
        const double kplus2pt = std::hypot(t2.phid1Px(), t2.phid1Py());
        const double kminus2pt = std::hypot(t2.phid2Px(), t2.phid2Py());

        if (kplus2pt > config.maxKaonPt || kminus2pt > config.maxKaonPt) {
          continue;
        }

        if (!selectionPID(t2.phid1TPC(), t2.phid1TOF(), t2.phid1TOFHit(),
                          config.strategyPID1, kplus2pt)) {
          continue;
        }

        if (!selectionPID(t2.phid2TPC(), t2.phid2TOF(), t2.phid2TOFHit(),
                          config.strategyPID2, kminus2pt)) {
          continue;
        }

        if (t2.phiMass() < config.minPhiMass2 || t2.phiMass() > config.maxPhiMass2) {
          continue;
        }

        TLorentzVector phi2;
        TLorentzVector k2p;
        TLorentzVector k2m;

        phi2.SetXYZM(t2.phiPx(), t2.phiPy(), t2.phiPz(), t2.phiMass());
        k2p.SetXYZM(t2.phid1Px(), t2.phid1Py(), t2.phid1Pz(), mKPDG);
        k2m.SetXYZM(t2.phid2Px(), t2.phid2Py(), t2.phid2Pz(), mKPDG);

        if (phi2.Pt() < config.minPhiPt || phi2.Pt() > config.maxPhiPt) {
          continue;
        }

        PhiCand cand;
        cand.phi = ROOT::Math::PtEtaPhiMVector(phi2.Pt(), phi2.Eta(), phi2.Phi(), phi2.M());
        cand.kplus = ROOT::Math::PtEtaPhiMVector(k2p.Pt(), k2p.Eta(), k2p.Phi(), mKPDG);
        cand.kminus = ROOT::Math::PtEtaPhiMVector(k2m.Pt(), k2m.Eta(), k2m.Phi(), mKPDG);

        cands2.emplace_back(std::move(cand));
      }

      if (cands1.empty() || cands2.empty()) {
        continue;
      }

      // ======================================================
      // Build mixed-event phi-phi pairs
      // ======================================================
      for (auto const& c1 : cands1) {

        TLorentzVector phi1;
        TLorentzVector k1p;
        TLorentzVector k1m;

        phi1.SetPtEtaPhiM(c1.phi.Pt(), c1.phi.Eta(), c1.phi.Phi(), c1.phi.M());
        k1p.SetPtEtaPhiM(c1.kplus.Pt(), c1.kplus.Eta(), c1.kplus.Phi(), mKPDG);
        k1m.SetPtEtaPhiM(c1.kminus.Pt(), c1.kminus.Eta(), c1.kminus.Phi(), mKPDG);

        for (auto const& c2 : cands2) {

          TLorentzVector phi2;
          TLorentzVector k2p;
          TLorentzVector k2m;

          phi2.SetPtEtaPhiM(c2.phi.Pt(), c2.phi.Eta(), c2.phi.Phi(), c2.phi.M());
          k2p.SetPtEtaPhiM(c2.kplus.Pt(), c2.kplus.Eta(), c2.kplus.Phi(), mKPDG);
          k2m.SetPtEtaPhiM(c2.kminus.Pt(), c2.kminus.Eta(), c2.kminus.Phi(), mKPDG);

          const double dM = getDeltaMPhi(phi1.M(), phi1.Pt(),
                                         phi2.M(), phi2.Pt());

          TLorentzVector pairPhiPhi = phi1 + phi2;

          if (pairPhiPhi.Pt() < config.minExoticPt) {
            continue;
          }

          const double minDR =
            minKaonDeltaR(c1.kplus, c1.kminus, c2.kplus, c2.kminus);

          const double dRphi =
            deltaR(phi1.Phi(), phi1.Eta(), phi2.Phi(), phi2.Eta());

          const double denom = std::abs(pairPhiPhi.Pt() - phi1.Pt());
          if (denom < 1e-9) {
            continue;
          }

          const double ptcorr = phi1.Pt() / denom;

          const double pt1 = phi1.Pt();
          const double pt2 = phi2.Pt();
          const double ptsum = pt1 + pt2;

          if (ptsum <= 0.0) {
            continue;
          }

          const double z = pt1 / ptsum;
          const double A = std::abs(pt1 - pt2) / ptsum;

          histos.fill(HIST("MEMassUnlike"),
                      pairPhiPhi.M(),
                      minDR,
                      pairPhiPhi.Pt(),
                      dRphi,
                      dM,
                      ptcorr);

          histos.fill(HIST("MEMassUnlike_DeltaRZA"),
                      pairPhiPhi.M(),
                      pairPhiPhi.Pt(),
                      pairPhiPhi.Pt() * dRphi,
                      z,
                      A,
                      dM);

          histos.fill(HIST("MEMassUnlike_AllVars"),
                      pairPhiPhi.M(),
                      pairPhiPhi.Pt(),
                      dRphi,
                      minDR,
                      phi1.M(),
                      phi2.M(),
                      dM,
                      ptcorr);
        }
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processMixedEventopti5,
                 "Process EventMixing for combinatorial background", false);

  void processMixedEvent(aod::RedPhiEvents& collisions, aod::PhiTracks& phitracks)
  {
    auto tracksTuple = std::make_tuple(phitracks);
    BinningTypeVertexContributor binningOnPositions{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::RedPhiEvents, aod::PhiTracks, BinningTypeVertexContributor> pair{
      binningOnPositions, nEvtMixing, -1, collisions, tracksTuple, &cache};

    // --- helpers (same as in processopti3) ---
    constexpr double mPhiPDG = 1.019461; // GeV/c^2

    const auto deltaMPhi = [=](double m1, double m2) {
      const double d1 = m1 - mPhiPDG;
      const double d2 = m2 - mPhiPDG;
      return std::sqrt(d1 * d1 + d2 * d2);
    };

    const auto deltaR = [](double phi1, double eta1, double phi2, double eta2) {
      const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
      const double deta = eta1 - eta2;
      return std::sqrt(dphi * dphi + deta * deta);
    };

    const auto minKaonDeltaR =
      [&](const ROOT::Math::PtEtaPhiMVector& kplusA,
          const ROOT::Math::PtEtaPhiMVector& kplusB,
          const ROOT::Math::PtEtaPhiMVector& kminusA,
          const ROOT::Math::PtEtaPhiMVector& kminusB) {
        // same-sign first (keep QA as in SE)
        const double dRkplus =
          deltaR(kplusA.Phi(), kplusA.Eta(), kplusB.Phi(), kplusB.Eta());
        const double dRkminus =
          deltaR(kminusA.Phi(), kminusA.Eta(), kminusB.Phi(), kminusB.Eta());
        histos.fill(HIST("hDeltaRkaonplus"), dRkplus);
        histos.fill(HIST("hDeltaRkaonminus"), dRkminus);

        // all other combinations
        const double dR_k1p_k1m =
          deltaR(kplusA.Phi(), kplusA.Eta(), kminusA.Phi(), kminusA.Eta());
        const double dR_k1p_k2m =
          deltaR(kplusA.Phi(), kplusA.Eta(), kminusB.Phi(), kminusB.Eta());
        const double dR_k2p_k1m =
          deltaR(kplusB.Phi(), kplusB.Eta(), kminusA.Phi(), kminusA.Eta());
        const double dR_k2p_k2m =
          deltaR(kplusB.Phi(), kplusB.Eta(), kminusB.Phi(), kminusB.Eta());

        double minDR = dRkplus;
        minDR = std::min(minDR, dRkminus);
        minDR = std::min(minDR, dR_k1p_k1m);
        minDR = std::min(minDR, dR_k1p_k2m);
        minDR = std::min(minDR, dR_k2p_k1m);
        minDR = std::min(minDR, dR_k2p_k2m);
        return minDR;
      };

    struct PhiCand {
      ROOT::Math::PtEtaPhiMVector phi;
      ROOT::Math::PtEtaPhiMVector kplus;
      ROOT::Math::PtEtaPhiMVector kminus;
    };

    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {
      // safety: should never happen but keep it
      if (collision1.index() == collision2.index()) {
        continue;
      }

      // optional event-level selection (same idea as in SE)
      if (config.additionalEvsel) {
        if (collision1.numPos() < 2 || collision1.numNeg() < 2) {
          continue;
        }
        if (collision2.numPos() < 2 || collision2.numNeg() < 2) {
          continue;
        }
      }

      std::vector<PhiCand> cands1, cands2;

      // --- build φ candidates for event 1 (φ1) ---
      for (auto const& t1 : tracks1) {
        const double kplus1pt = std::hypot(t1.phid1Px(), t1.phid1Py());
        const double kminus1pt = std::hypot(t1.phid2Px(), t1.phid2Py());

        if (kplus1pt > config.maxKaonPt || kminus1pt > config.maxKaonPt)
          continue;
        if (!selectionPID(t1.phid1TPC(), t1.phid1TOF(), t1.phid1TOFHit(), config.strategyPID1, kplus1pt))
          continue;
        if (!selectionPID(t1.phid2TPC(), t1.phid2TOF(), t1.phid2TOFHit(), config.strategyPID2, kminus1pt))
          continue;

        TLorentzVector phi1, k1p, k1m;
        phi1.SetXYZM(t1.phiPx(), t1.phiPy(), t1.phiPz(), t1.phiMass());
        k1p.SetXYZM(t1.phid1Px(), t1.phid1Py(), t1.phid1Pz(), 0.493);
        k1m.SetXYZM(t1.phid2Px(), t1.phid2Py(), t1.phid2Pz(), 0.493);

        if (t1.phiMass() < config.minPhiMass1 || t1.phiMass() > config.maxPhiMass1)
          continue;
        if (phi1.Pt() < config.minPhiPt || phi1.Pt() > config.maxPhiPt)
          continue;

        PhiCand cand;
        cand.phi = ROOT::Math::PtEtaPhiMVector(phi1.Pt(), phi1.Eta(), phi1.Phi(), phi1.M());
        cand.kplus = ROOT::Math::PtEtaPhiMVector(k1p.Pt(), k1p.Eta(), k1p.Phi(), 0.493);
        cand.kminus = ROOT::Math::PtEtaPhiMVector(k1m.Pt(), k1m.Eta(), k1m.Phi(), 0.493);

        cands1.emplace_back(std::move(cand));
      }

      // --- build φ candidates for event 2 (φ2) ---
      for (auto const& t2 : tracks2) {
        const double kplus2pt = std::hypot(t2.phid1Px(), t2.phid1Py());
        const double kminus2pt = std::hypot(t2.phid2Px(), t2.phid2Py());

        if (kplus2pt > config.maxKaonPt || kminus2pt > config.maxKaonPt)
          continue;
        if (!selectionPID(t2.phid1TPC(), t2.phid1TOF(), t2.phid1TOFHit(), config.strategyPID1, kplus2pt))
          continue;
        if (!selectionPID(t2.phid2TPC(), t2.phid2TOF(), t2.phid2TOFHit(), config.strategyPID2, kminus2pt))
          continue;

        TLorentzVector phi2, k2p, k2m;
        phi2.SetXYZM(t2.phiPx(), t2.phiPy(), t2.phiPz(), t2.phiMass());
        k2p.SetXYZM(t2.phid1Px(), t2.phid1Py(), t2.phid1Pz(), 0.493);
        k2m.SetXYZM(t2.phid2Px(), t2.phid2Py(), t2.phid2Pz(), 0.493);

        if (t2.phiMass() < config.minPhiMass2 || t2.phiMass() > config.maxPhiMass2)
          continue;
        if (phi2.Pt() < config.minPhiPt || phi2.Pt() > config.maxPhiPt)
          continue;

        PhiCand cand;
        cand.phi = ROOT::Math::PtEtaPhiMVector(phi2.Pt(), phi2.Eta(), phi2.Phi(), phi2.M());
        cand.kplus = ROOT::Math::PtEtaPhiMVector(k2p.Pt(), k2p.Eta(), k2p.Phi(), 0.493);
        cand.kminus = ROOT::Math::PtEtaPhiMVector(k2m.Pt(), k2m.Eta(), k2m.Phi(), 0.493);

        cands2.emplace_back(std::move(cand));
      }

      if (cands1.empty() || cands2.empty()) {
        continue;
      }

      // --- build mixed-event pairs and fill MEMassUnlike ---
      for (auto const& c1 : cands1) {
        TLorentzVector phi1;
        phi1.SetPtEtaPhiM(c1.phi.Pt(), c1.phi.Eta(), c1.phi.Phi(), c1.phi.M());

        for (auto const& c2 : cands2) {
          TLorentzVector phi2;
          phi2.SetPtEtaPhiM(c2.phi.Pt(), c2.phi.Eta(), c2.phi.Phi(), c2.phi.M());

          const double dM = deltaMPhi(phi1.M(), phi2.M());
          if (dM > maxDeltaMPhi)
            continue;

          TLorentzVector pairPhiPhi = phi1 + phi2;
          if (pairPhiPhi.M() < config.minExoticMass || pairPhiPhi.M() > config.maxExoticMass)
            continue;

          const double minDR = minKaonDeltaR(c1.kplus, c2.kplus, c1.kminus, c2.kminus);
          const double dR = deltaR(phi1.Phi(), phi1.Eta(), phi2.Phi(), phi2.Eta());

          // same definition as SE
          const double ptcorr = (pairPhiPhi.Pt() - phi1.Pt() != 0.)
                                  ? phi1.Pt() / (pairPhiPhi.Pt() - phi1.Pt())
                                  : 0.;

          histos.fill(HIST("MEMassUnlike"),
                      pairPhiPhi.M(),  // M(phi-phi)
                      minDR,           // min ΔR among all kaon pairs
                      pairPhiPhi.Pt(), // pT(phi-phi)
                      dR,              // ΔR(phi1, phi2)
                      dM,              // Δm(phi)
                      ptcorr);         // pT correlation

          // --- NEW: compute z and A from phi candidates (no cuts) ---
          const double pt1 = phi1.Pt();
          const double pt2 = phi2.Pt();
          const double ptsum = pt1 + pt2;
          if (ptsum <= 0.0)
            continue;
          const double z = pt1 / ptsum;
          const double A = std::abs(pt1 - pt2) / ptsum;
          // --- Fill NEW THnSparse (no cuts) ---
          histos.fill(HIST("MEMassUnlike_DeltaRZA"), pairPhiPhi.M(), pairPhiPhi.Pt(), pairPhiPhi.Pt() * dR, z, A, dM);
        }
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processMixedEvent,
                 "Process EventMixing for combinatorial background", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<doublephimeson>(cfgc)}; }
