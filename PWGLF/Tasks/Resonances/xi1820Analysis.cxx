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

/// \file xi1820Analysis.cxx
/// \brief Invariant Mass Reconstruction of Xi(1820) Resonance
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <Math/Vector4D.h>
#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct Xi1820Analysis {
  SliceCache cache;
  Preslice<aod::ResoV0s> perResoCollisionV0 = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoTracks> perResoCollisionTrack = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoMicroTracks> perResoCollisionMicroTrack = aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Constants
  static constexpr float kSmallMomentumDenominator = 1e-10f; // Small value to avoid division by zero
  static constexpr float kMaxDcaToPv = 1.0f;                 // Maximum DCA to primary vertex
  static constexpr int kPdgXi1820 = 123314;                  // o2-linter: disable=pdg/explicit-code (Xi(1820) PDG code not available in PDG_t or o2::constants::physics::Pdg)
  static constexpr int kExpectedDaughters = 2;               // Expected number of daughters for two-body decay

  // Axes
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0}, "pT"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pT (QA)"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Centrality"};

  // Invariant mass range for Xi(1820) → Λ + K
  Configurable<float> cInvMassStart{"cInvMassStart", 1.6, "Invariant mass start (GeV/c^2)"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 2.2, "Invariant mass end (GeV/c^2)"};
  Configurable<int> cInvMassBins{"cInvMassBins", 600, "Invariant mass bins"};

  // Basic pre-selections
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Minimum pT for candidates"};
  Configurable<float> cMaxEtaCut{"cMaxEtaCut", 0.8, "Maximum |eta|"};

  // Kaon track selections
  Configurable<float> cKaonPtMin{"cKaonPtMin", 0.15, "Minimum kaon pT"};
  Configurable<float> cKaonEtaMax{"cKaonEtaMax", 0.8, "Maximum kaon |eta|"};
  Configurable<float> cKaonDCAxyMax{"cKaonDCAxyMax", 0.1, "Maximum kaon DCAxy to PV"};
  Configurable<float> cKaonDCAzMax{"cKaonDCAzMax", 0.2, "Maximum kaon DCAz to PV"};
  Configurable<bool> cKaonTPCRefit{"cKaonTPCRefit", true, "Require TPC refit for kaon"};
  Configurable<bool> cKaonITSRefit{"cKaonITSRefit", true, "Require ITS refit for kaon"};
  Configurable<int> cKaonTPCNClusMin{"cKaonTPCNClusMin", 70, "Minimum TPC clusters for kaon"};
  Configurable<int> cKaonITSNClusMin{"cKaonITSNClusMin", 2, "Minimum ITS clusters for kaon"};

  // Kaon PID selections
  Configurable<float> cKaonTPCNSigmaMax{"cKaonTPCNSigmaMax", 3.0, "Maximum TPC NSigma for kaon (if not using pT-dependent)"};
  Configurable<float> cKaonTOFNSigmaMax{"cKaonTOFNSigmaMax", 3.0, "Maximum TOF NSigma for kaon (if not using pT-dependent)"};
  Configurable<bool> cKaonUsePtDepPID{"cKaonUsePtDepPID", false, "Use pT-dependent PID cuts"};
  Configurable<std::vector<float>> cKaonPIDPtBins{"cKaonPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for PID cuts (N+1 values for N bins)"};
  Configurable<std::vector<float>> cKaonTPCNSigmaCuts{"cKaonTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (N values)"};
  Configurable<std::vector<float>> cKaonTOFNSigmaCuts{"cKaonTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (N values)"};
  Configurable<std::vector<int>> cKaonTOFRequired{"cKaonTOFRequired", {0, 0, 1, 1}, "Require TOF per pT bin (N values, 0=false, 1=true)"};

  // V0 (Lambda) selections
  Configurable<double> cV0MinCosPA{"cV0MinCosPA", 0.995, "V0 minimum pointing angle cosine"};
  Configurable<double> cV0MaxDaughDCA{"cV0MaxDaughDCA", 1.0, "V0 daughter DCA Maximum"};
  Configurable<double> cV0MassWindow{"cV0MassWindow", 0.005, "Mass window for Lambda selection"};
  Configurable<double> cMaxV0Etacut{"cMaxV0Etacut", 0.8, "V0 maximum eta cut"};
  Configurable<float> cV0RadiusMin{"cV0RadiusMin", 0.5, "V0 decay radius min"};
  Configurable<float> cV0RadiusMax{"cV0RadiusMax", 200.0, "V0 decay radius max"};
  Configurable<float> cV0DauPosDCAtoPVMin{"cV0DauPosDCAtoPVMin", 0.05, "V0 positive daughter DCA to PV min"};
  Configurable<float> cV0DauNegDCAtoPVMin{"cV0DauNegDCAtoPVMin", 0.05, "V0 negative daughter DCA to PV min"};
  Configurable<float> cV0ProperLifetimeMax{"cV0ProperLifetimeMax", 30.0, "Lambda proper lifetime max (cm/c)"};

  // K0s selections
  Configurable<double> cK0sMinCosPA{"cK0sMinCosPA", 0.97, "K0s minimum pointing angle cosine"};
  Configurable<double> cK0sMaxDaughDCA{"cK0sMaxDaughDCA", 1.0, "K0s daughter DCA Maximum"};
  Configurable<double> cK0sMassWindow{"cK0sMassWindow", 0.0043, "Mass window for K0s selection"};
  Configurable<float> cK0sProperLifetimeMax{"cK0sProperLifetimeMax", 20.0, "K0s proper lifetime max (cm/c)"};
  Configurable<float> cK0sArmenterosQtMin{"cK0sArmenterosQtMin", 0.0, "K0s Armenteros qt min"};
  Configurable<float> cK0sArmenterosAlphaMax{"cK0sArmenterosAlphaMax", 0.8, "K0s Armenteros alpha max"};
  Configurable<float> cK0sDauPosDCAtoPVMin{"cK0sDauPosDCAtoPVMin", 0.1, "K0s positive daughter DCA to PV min"};
  Configurable<float> cK0sDauNegDCAtoPVMin{"cK0sDauNegDCAtoPVMin", 0.1, "K0s negative daughter DCA to PV min"};
  Configurable<float> cK0sRadiusMin{"cK0sRadiusMin", 0.5, "K0s decay radius min"};
  Configurable<float> cK0sRadiusMax{"cK0sRadiusMax", 100.0, "K0s decay radius max"};
  Configurable<bool> cK0sCrossMassRejection{"cK0sCrossMassRejection", true, "Enable Lambda mass rejection for K0s"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - centrality"};

  // Track type selection
  Configurable<bool> cUseMicroTracks{"cUseMicroTracks", false, "Use ResoMicroTracks instead of ResoTracks"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  BinningTypeVertexContributor colBinning{{cfgVtxBins, cfgMultBins}, true};

  void init(InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^{2})"};
    AxisSpec lambdaMassAxis = {200, 1.08, 1.16, "#Lambda mass (GeV/#it{c}^{2})"};
    AxisSpec dcaAxis = {200, 0., 2.0, "DCA (cm)"};
    AxisSpec dcaxyAxis = {200, -1.0, 1.0, "DCA_{xy} (cm)"};
    AxisSpec dcazAxis = {200, -2.0, 2.0, "DCA_{z} (cm)"};
    AxisSpec cosPAAxis = {1000, 0.95, 1.0, "cos(PA)"};
    AxisSpec radiusAxis = {200, 0, 200, "Radius (cm)"};
    AxisSpec lifetimeAxis = {200, 0, 50, "Proper lifetime (cm/c)"};
    AxisSpec nsigmaAxis = {100, -5.0, 5.0, "N#sigma"};

    // Event QA histograms
    histos.add("Event/posZ", "Event vertex Z position", kTH1F, {{200, -20., 20., "V_{z} (cm)"}});
    histos.add("Event/centrality", "Event centrality distribution", kTH1F, {centAxis});
    histos.add("Event/posZvsCent", "Vertex Z vs Centrality", kTH2F, {{200, -20., 20., "V_{z} (cm)"}, centAxis});
    histos.add("Event/nV0s", "Number of V0s per event", kTH1F, {{200, 0., 200., "N_{V0s}"}});
    histos.add("Event/nKaons", "Number of kaons per event", kTH1F, {{200, 0., 200., "N_{kaons}"}});
    histos.add("Event/nV0sAfterCuts", "Number of V0s per event after cuts", kTH1F, {{100, 0., 100., "N_{V0s}"}});
    histos.add("Event/nKaonsAfterCuts", "Number of kaons per event after cuts", kTH1F, {{100, 0., 100., "N_{kaons}"}});

    // Lambda QA histograms
    histos.add("QAbefore/lambdaMass", "Lambda mass before cuts", kTH1F, {lambdaMassAxis});
    histos.add("QAbefore/lambdaPt", "Lambda pT before cuts", kTH1F, {ptAxisQA});
    histos.add("QAbefore/lambdaEta", "Lambda eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAbefore/lambdaCosPA", "Lambda CosPA before cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAbefore/lambdaRadius", "Lambda radius before cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAbefore/lambdaDauDCA", "Lambda daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAbefore/lambdaProperLifetime", "Lambda proper lifetime before cuts", kTH2F, {ptAxisQA, lifetimeAxis});
    histos.add("QAbefore/lambdaDauPosDCA", "Lambda positive daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAbefore/lambdaDauNegDCA", "Lambda negative daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});

    histos.add("QAafter/lambdaMass", "Lambda mass after cuts", kTH1F, {lambdaMassAxis});
    histos.add("QAafter/lambdaPt", "Lambda pT after cuts", kTH1F, {ptAxisQA});
    histos.add("QAafter/lambdaEta", "Lambda eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAafter/lambdaCosPA", "Lambda CosPA after cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAafter/lambdaRadius", "Lambda radius after cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAafter/lambdaDauDCA", "Lambda daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAafter/lambdaProperLifetime", "Lambda proper lifetime after cuts", kTH2F, {ptAxisQA, lifetimeAxis});
    histos.add("QAafter/lambdaDauPosDCA", "Lambda positive daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAafter/lambdaDauNegDCA", "Lambda negative daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});

    // Kaon QA histograms
    histos.add("QAbefore/kaonPt", "Kaon pT before cuts", kTH1F, {ptAxisQA});
    histos.add("QAbefore/kaonEta", "Kaon eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAbefore/kaonDCAxy", "Kaon DCAxy before cuts", kTH2F, {ptAxisQA, dcaxyAxis});
    histos.add("QAbefore/kaonDCAz", "Kaon DCAz before cuts", kTH2F, {ptAxisQA, dcazAxis});
    histos.add("QAbefore/kaonTPCNcls", "Kaon TPC clusters before cuts", kTH1F, {{160, 0, 160, "N_{TPC clusters}"}});
    histos.add("QAbefore/kaonITSNcls", "Kaon ITS clusters before cuts", kTH1F, {{10, 0, 10, "N_{ITS clusters}"}});
    histos.add("QAbefore/kaonTPCNSigma", "Kaon TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    histos.add("QAbefore/kaonTOFNSigma", "Kaon TOF NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});

    histos.add("QAafter/kaonPt", "Kaon pT after cuts", kTH1F, {ptAxisQA});
    histos.add("QAafter/kaonEta", "Kaon eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAafter/kaonDCAxy", "Kaon DCAxy after cuts", kTH2F, {ptAxisQA, dcaxyAxis});
    histos.add("QAafter/kaonDCAz", "Kaon DCAz after cuts", kTH2F, {ptAxisQA, dcazAxis});
    histos.add("QAafter/kaonTPCNcls", "Kaon TPC clusters after cuts", kTH1F, {{160, 0, 160, "N_{TPC clusters}"}});
    histos.add("QAafter/kaonITSNcls", "Kaon ITS clusters after cuts", kTH1F, {{10, 0, 10, "N_{ITS clusters}"}});
    histos.add("QAafter/kaonTPCNSigma", "Kaon TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    histos.add("QAafter/kaonTOFNSigma", "Kaon TOF NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});

    // Resonance histograms - 4 combinations
    // K+ Lambda
    histos.add("xi1820/kplus_lambda/hInvMassKplusLambda", "Invariant mass of Xi(1820) → K^{+} + #Lambda", kTH1F, {invMassAxis});
    histos.add("xi1820/kplus_lambda/hInvMassKplusLambda_Mix", "Mixed event Invariant mass of Xi(1820) → K^{+} + #Lambda", kTH1F, {invMassAxis});
    histos.add("xi1820/kplus_lambda/hMassPtCentKplusLambda", "Xi(1820) mass vs pT vs cent (K^{+}#Lambda)", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("xi1820/kplus_lambda/hMassPtCentKplusLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent (K^{+}#Lambda)", kTH3F, {invMassAxis, ptAxis, centAxis});

    // K+ Anti-Lambda
    histos.add("xi1820/kplus_antilambda/hInvMassKplusAntiLambda", "Invariant mass of Xi(1820) → K^{+} + #bar{#Lambda}", kTH1F, {invMassAxis});
    histos.add("xi1820/kplus_antilambda/hInvMassKplusAntiLambda_Mix", "Mixed event Invariant mass of Xi(1820) → K^{+} + #bar{#Lambda}", kTH1F, {invMassAxis});
    histos.add("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda", "Xi(1820) mass vs pT vs cent (K^{+}#bar{#Lambda})", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent (K^{+}#bar{#Lambda})", kTH3F, {invMassAxis, ptAxis, centAxis});

    // K- Lambda
    histos.add("xi1820/kminus_lambda/hInvMassKminusLambda", "Invariant mass of Xi(1820) → K^{-} + #Lambda", kTH1F, {invMassAxis});
    histos.add("xi1820/kminus_lambda/hInvMassKminusLambda_Mix", "Mixed event Invariant mass of Xi(1820) → K^{-} + #Lambda", kTH1F, {invMassAxis});
    histos.add("xi1820/kminus_lambda/hMassPtCentKminusLambda", "Xi(1820) mass vs pT vs cent (K^{-}#Lambda)", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("xi1820/kminus_lambda/hMassPtCentKminusLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent (K^{-}#Lambda)", kTH3F, {invMassAxis, ptAxis, centAxis});

    // K- Anti-Lambda
    histos.add("xi1820/kminus_antilambda/hInvMassKminusAntiLambda", "Invariant mass of Xi(1820) → K^{-} + #bar{#Lambda}", kTH1F, {invMassAxis});
    histos.add("xi1820/kminus_antilambda/hInvMassKminusAntiLambda_Mix", "Mixed event Invariant mass of Xi(1820) → K^{-} + #bar{#Lambda}", kTH1F, {invMassAxis});
    histos.add("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda", "Xi(1820) mass vs pT vs cent (K^{-}#bar{#Lambda})", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent (K^{-}#bar{#Lambda})", kTH3F, {invMassAxis, ptAxis, centAxis});

    // K0s + Lambda
    histos.add("xi1820/k0s_lambda/hInvMassK0sLambda", "Invariant mass of Xi(1820) → K^{0}_{S} + #Lambda", kTH1F, {invMassAxis});
    histos.add("xi1820/k0s_lambda/hInvMassK0sLambda_Mix", "Mixed event Invariant mass of Xi(1820) → K^{0}_{S} + #Lambda", kTH1F, {invMassAxis});
    histos.add("xi1820/k0s_lambda/hMassPtCentK0sLambda", "Xi(1820) mass vs pT vs cent (K^{0}_{S}#Lambda)", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("xi1820/k0s_lambda/hMassPtCentK0sLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent (K^{0}_{S}#Lambda)", kTH3F, {invMassAxis, ptAxis, centAxis});

    // K0s + Anti-Lambda
    histos.add("xi1820/k0s_antilambda/hInvMassK0sAntiLambda", "Invariant mass of Xi(1820) → K^{0}_{S} + #bar{#Lambda}", kTH1F, {invMassAxis});
    histos.add("xi1820/k0s_antilambda/hInvMassK0sAntiLambda_Mix", "Mixed event Invariant mass of Xi(1820) → K^{0}_{S} + #bar{#Lambda}", kTH1F, {invMassAxis});
    histos.add("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda", "Xi(1820) mass vs pT vs cent (K^{0}_{S}#bar{#Lambda})", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda_Mix", "Mixed event Xi(1820) mass vs pT vs cent (K^{0}_{S}#bar{#Lambda})", kTH3F, {invMassAxis, ptAxis, centAxis});

    // MC truth histograms
    AxisSpec etaAxis = {100, -2.0, 2.0, "#eta"};
    AxisSpec rapidityAxis = {100, -2.0, 2.0, "y"};

    histos.add("MC/hMCGenXi1820Pt", "MC Generated Xi(1820) pT", kTH1F, {ptAxis});
    histos.add("MC/hMCGenXi1820PtEta", "MC Generated Xi(1820) pT vs eta", kTH2F, {ptAxis, etaAxis});
    histos.add("MC/hMCGenXi1820Y", "MC Generated Xi(1820) rapidity", kTH1F, {rapidityAxis});
    histos.add("MC/hMCRecXi1820Pt", "MC Reconstructed Xi(1820) pT", kTH1F, {ptAxis});
    histos.add("MC/hMCRecXi1820PtEta", "MC Reconstructed Xi(1820) pT vs eta", kTH2F, {ptAxis, etaAxis});

    // MC truth invariant mass (from MC particles)
    histos.add("MC/hMCTruthInvMassKplusLambda", "MC Truth Inv Mass K^{+}#Lambda", kTH1F, {invMassAxis});
    histos.add("MC/hMCTruthInvMassKminusAntiLambda", "MC Truth Inv Mass K^{-}#bar{#Lambda}", kTH1F, {invMassAxis});
    histos.add("MC/hMCTruthInvMassK0sLambda", "MC Truth Inv Mass K^{0}_{S}#Lambda", kTH1F, {invMassAxis});
    histos.add("MC/hMCTruthInvMassK0sAntiLambda", "MC Truth Inv Mass K^{0}_{S}#bar{#Lambda}", kTH1F, {invMassAxis});

    // MC truth invariant mass vs pT (2D)
    histos.add("MC/hMCTruthMassPtKplusLambda", "MC Truth Mass vs pT K^{+}#Lambda", kTH2F, {invMassAxis, ptAxis});
    histos.add("MC/hMCTruthMassPtKminusAntiLambda", "MC Truth Mass vs pT K^{-}#bar{#Lambda}", kTH2F, {invMassAxis, ptAxis});
    histos.add("MC/hMCTruthMassPtK0sLambda", "MC Truth Mass vs pT K^{0}_{S}#Lambda", kTH2F, {invMassAxis, ptAxis});
    histos.add("MC/hMCTruthMassPtK0sAntiLambda", "MC Truth Mass vs pT K^{0}_{S}#bar{#Lambda}", kTH2F, {invMassAxis, ptAxis});

    // K0s QA histograms
    histos.add("QAbefore/k0sMass", "K0s mass before cuts", kTH1F, {{100, 0.4, 0.6, "K^{0}_{S} mass (GeV/#it{c}^{2})"}});
    histos.add("QAbefore/k0sPt", "K0s pT before cuts", kTH1F, {ptAxisQA});
    histos.add("QAbefore/k0sEta", "K0s eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAbefore/k0sCosPA", "K0s CosPA before cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAbefore/k0sRadius", "K0s radius before cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAbefore/k0sDauDCA", "K0s daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAbefore/k0sProperLifetime", "K0s proper lifetime before cuts", kTH2F, {ptAxisQA, lifetimeAxis});

    histos.add("QAafter/k0sMass", "K0s mass after cuts", kTH1F, {{100, 0.4, 0.6, "K^{0}_{S} mass (GeV/#it{c}^{2})"}});
    histos.add("QAafter/k0sPt", "K0s pT after cuts", kTH1F, {ptAxisQA});
    histos.add("QAafter/k0sEta", "K0s eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAafter/k0sCosPA", "K0s CosPA after cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAafter/k0sRadius", "K0s radius after cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAafter/k0sDauDCA", "K0s daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAafter/k0sProperLifetime", "K0s proper lifetime after cuts", kTH2F, {ptAxisQA, lifetimeAxis});
  }

  // Lambda/Anti-Lambda selection
  template <typename CollisionType, typename V0Type>
  bool v0Cut(const CollisionType& collision, const V0Type& v0, bool isLambda)
  {
    // Basic kinematic cuts
    if (std::abs(v0.eta()) > cMaxV0Etacut)
      return false;
    if (v0.pt() < cMinPtcut)
      return false;

    // Topological cuts
    if (v0.v0CosPA() < cV0MinCosPA)
      return false;
    if (v0.daughDCA() > cV0MaxDaughDCA)
      return false;

    // Daughter DCA to PV cuts
    if (std::abs(v0.dcapostopv()) < cV0DauPosDCAtoPVMin)
      return false;
    if (std::abs(v0.dcanegtopv()) < cV0DauNegDCAtoPVMin)
      return false;

    // Radius cuts
    auto radius = v0.transRadius();
    if (radius < cV0RadiusMin || radius > cV0RadiusMax)
      return false;

    // Proper lifetime cut
    float dx = v0.decayVtxX() - collision.posX();
    float dy = v0.decayVtxY() - collision.posY();
    float dz = v0.decayVtxZ() - collision.posZ();
    float l = std::sqrt(dx * dx + dy * dy + dz * dz);
    float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
    auto properLifetime = (l / (p + kSmallMomentumDenominator)) * MassLambda;
    if (properLifetime > cV0ProperLifetimeMax)
      return false;

    // Mass window
    if (isLambda) {
      if (std::abs(v0.mLambda() - MassLambda) > cV0MassWindow)
        return false;
    } else {
      if (std::abs(v0.mAntiLambda() - MassLambda) > cV0MassWindow)
        return false;
    }

    return true;
  }

  // K0s selection
  template <typename CollisionType, typename V0Type>
  bool k0sCut(const CollisionType& collision, const V0Type& v0)
  {
    // Basic kinematic cuts
    if (std::abs(v0.eta()) > cMaxV0Etacut)
      return false;
    if (v0.pt() < cMinPtcut)
      return false;

    // Topological cuts
    if (v0.v0CosPA() < cK0sMinCosPA)
      return false;
    if (v0.daughDCA() > cK0sMaxDaughDCA)
      return false;

    // Daughter DCA to PV cuts
    if (std::abs(v0.dcapostopv()) < cK0sDauPosDCAtoPVMin)
      return false;
    if (std::abs(v0.dcanegtopv()) < cK0sDauNegDCAtoPVMin)
      return false;

    // Radius cuts
    auto radius = v0.transRadius();
    if (radius < cK0sRadiusMin || radius > cK0sRadiusMax)
      return false;

    // DCA to PV
    if (std::abs(v0.dcav0topv()) > kMaxDcaToPv)
      return false;

    // Proper lifetime cut
    float dx = v0.decayVtxX() - collision.posX();
    float dy = v0.decayVtxY() - collision.posY();
    float dz = v0.decayVtxZ() - collision.posZ();
    float l = std::sqrt(dx * dx + dy * dy + dz * dz);
    float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
    auto properLifetime = (l / (p + kSmallMomentumDenominator)) * MassK0Short;
    if (properLifetime > cK0sProperLifetimeMax)
      return false;

    // Mass window
    if (std::abs(v0.mK0Short() - MassK0Short) > cK0sMassWindow)
      return false;

    // Competing V0 rejection: remove (Anti)Λ
    if (cK0sCrossMassRejection) {
      if (std::abs(v0.mLambda() - MassLambda) < cK0sMassWindow)
        return false;
      if (std::abs(v0.mAntiLambda() - MassLambda) < cK0sMassWindow)
        return false;
    }

    return true;
  }

  // Helper function to find pT bin index
  int getPtBinIndex(float pt)
  {
    auto ptBins = static_cast<std::vector<float>>(cKaonPIDPtBins);
    for (size_t i = 0; i < ptBins.size() - 1; i++) {
      if (pt >= ptBins[i] && pt < ptBins[i + 1]) {
        return i;
      }
    }
    return -1; // should not happen if bins are properly configured
  }

  // Kaon PID selection
  template <bool IsResoMicrotrack, typename TrackType>
  bool kaonPidCut(const TrackType& track)
  {
    float pt = track.pt();

    if constexpr (IsResoMicrotrack) {
      // For ResoMicroTracks - decode PID from flags
      float tpcNSigma = o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(track.pidNSigmaKaFlag());
      float tofNSigma = track.hasTOF() ? o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(track.pidNSigmaKaFlag()) : 999.f;

      if (cKaonUsePtDepPID) {
        // pT-dependent PID with binning
        int ptBin = getPtBinIndex(pt);
        if (ptBin < 0)
          return false; // safety check

        auto tpcCuts = static_cast<std::vector<float>>(cKaonTPCNSigmaCuts);
        auto tofCuts = static_cast<std::vector<float>>(cKaonTOFNSigmaCuts);
        auto tofRequired = static_cast<std::vector<int>>(cKaonTOFRequired);

        // Check array sizes
        if (ptBin >= static_cast<int>(tpcCuts.size()) ||
            ptBin >= static_cast<int>(tofCuts.size()) ||
            ptBin >= static_cast<int>(tofRequired.size())) {
          return false; // safety check
        }

        // Apply TPC cut
        if (std::abs(tpcNSigma) >= tpcCuts[ptBin])
          return false;

        // Apply TOF requirement and cut
        if (tofRequired[ptBin] != 0) {
          if (!track.hasTOF())
            return false;
          if (std::abs(tofNSigma) >= tofCuts[ptBin])
            return false;
        } else {
          // TOF optional but apply cut if present
          if (track.hasTOF() && std::abs(tofNSigma) >= tofCuts[ptBin])
            return false;
        }

        return true;
      } else {
        // Standard PID
        bool tpcPass = std::abs(tpcNSigma) < cKaonTPCNSigmaMax;
        bool tofPass = track.hasTOF() ? std::abs(tofNSigma) < cKaonTOFNSigmaMax : true;
        return tpcPass && tofPass;
      }
    } else {
      // For ResoTracks - direct access
      float tpcNSigma = track.tpcNSigmaKa();
      float tofNSigma = track.hasTOF() ? track.tofNSigmaKa() : 999.f;

      if (cKaonUsePtDepPID) {
        // pT-dependent PID with binning
        int ptBin = getPtBinIndex(pt);
        if (ptBin < 0)
          return false; // safety check

        auto tpcCuts = static_cast<std::vector<float>>(cKaonTPCNSigmaCuts);
        auto tofCuts = static_cast<std::vector<float>>(cKaonTOFNSigmaCuts);
        auto tofRequired = static_cast<std::vector<int>>(cKaonTOFRequired);

        // Check array sizes
        if (ptBin >= static_cast<int>(tpcCuts.size()) ||
            ptBin >= static_cast<int>(tofCuts.size()) ||
            ptBin >= static_cast<int>(tofRequired.size())) {
          return false; // safety check
        }

        // Apply TPC cut
        if (std::abs(tpcNSigma) >= tpcCuts[ptBin])
          return false;

        // Apply TOF requirement and cut
        if (tofRequired[ptBin] != 0) {
          if (!track.hasTOF())
            return false;
          if (std::abs(tofNSigma) >= tofCuts[ptBin])
            return false;
        } else {
          // TOF optional but apply cut if present
          if (track.hasTOF() && std::abs(tofNSigma) >= tofCuts[ptBin])
            return false;
        }

        return true;
      } else {
        // Standard PID
        bool tpcPass = std::abs(tpcNSigma) < cKaonTPCNSigmaMax;
        bool tofPass = track.hasTOF() ? std::abs(tofNSigma) < cKaonTOFNSigmaMax : true;
        return tpcPass && tofPass;
      }
    }
  }

  // Kaon track selection (for both ResoTracks and ResoMicroTracks)
  template <bool IsResoMicrotrack, typename TrackType>
  bool kaonCut(const TrackType& track)
  {
    // Basic kinematic cuts
    if (track.pt() < cKaonPtMin)
      return false;
    if (std::abs(track.eta()) > cKaonEtaMax)
      return false;

    // DCA cuts - different access for ResoMicroTracks
    if constexpr (IsResoMicrotrack) {
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags()) > cKaonDCAxyMax)
        return false;
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags()) > cKaonDCAzMax)
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cKaonDCAxyMax)
        return false;
      if (std::abs(track.dcaZ()) > cKaonDCAzMax)
        return false;
    }

    // Track quality cuts - check if fields are available (only for ResoTracks)
    if constexpr (!IsResoMicrotrack) {
      if constexpr (requires { track.tpcNClsFound(); }) {
        if (track.tpcNClsFound() < cKaonTPCNClusMin)
          return false;
      }
      if constexpr (requires { track.itsNCls(); }) {
        if (track.itsNCls() < cKaonITSNClusMin)
          return false;
      }
    }

    // PID selection
    if (!kaonPidCut<IsResoMicrotrack>(track))
      return false;

    return true;
  }

  template <bool IsMix, bool IsResoMicrotrack, typename CollisionT, typename V0sT, typename TracksT>
  void fill(const CollisionT& collision, const V0sT& v0s, const TracksT& tracks)
  {
    auto cent = collision.cent();

    // Fill event QA histograms (only for same-event)
    if constexpr (!IsMix) {
      histos.fill(HIST("Event/posZ"), collision.posZ());
      histos.fill(HIST("Event/centrality"), cent);
      histos.fill(HIST("Event/posZvsCent"), collision.posZ(), cent);
      histos.fill(HIST("Event/nV0s"), v0s.size());
      histos.fill(HIST("Event/nKaons"), tracks.size());
    }

    // Count candidates after cuts
    int nV0sAfterCuts = 0;
    int nKaonsAfterCuts = 0;

    // Loop over kaon candidates
    for (const auto& kaon : tracks) {
      // QA before cuts
      if constexpr (!IsMix) {
        histos.fill(HIST("QAbefore/kaonPt"), kaon.pt());
        histos.fill(HIST("QAbefore/kaonEta"), kaon.eta());
        if constexpr (IsResoMicrotrack) {
          histos.fill(HIST("QAbefore/kaonDCAxy"), kaon.pt(), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(kaon.trackSelectionFlags()));
          histos.fill(HIST("QAbefore/kaonDCAz"), kaon.pt(), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(kaon.trackSelectionFlags()));
          histos.fill(HIST("QAbefore/kaonTPCNSigma"), kaon.pt(), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(kaon.pidNSigmaKaFlag()));
          if (kaon.hasTOF()) {
            histos.fill(HIST("QAbefore/kaonTOFNSigma"), kaon.pt(), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(kaon.pidNSigmaKaFlag()));
          }
        } else {
          histos.fill(HIST("QAbefore/kaonDCAxy"), kaon.pt(), kaon.dcaXY());
          histos.fill(HIST("QAbefore/kaonDCAz"), kaon.pt(), kaon.dcaZ());
          histos.fill(HIST("QAbefore/kaonTPCNSigma"), kaon.pt(), kaon.tpcNSigmaKa());
          if (kaon.hasTOF()) {
            histos.fill(HIST("QAbefore/kaonTOFNSigma"), kaon.pt(), kaon.tofNSigmaKa());
          }
          if constexpr (requires { kaon.tpcNClsFound(); }) {
            histos.fill(HIST("QAbefore/kaonTPCNcls"), kaon.tpcNClsFound());
          }
          if constexpr (requires { kaon.itsNCls(); }) {
            histos.fill(HIST("QAbefore/kaonITSNcls"), kaon.itsNCls());
          }
        }
      }

      if (!kaonCut<IsResoMicrotrack>(kaon))
        continue;

      if constexpr (!IsMix) {
        nKaonsAfterCuts++;
        // QA after cuts
        histos.fill(HIST("QAafter/kaonPt"), kaon.pt());
        histos.fill(HIST("QAafter/kaonEta"), kaon.eta());
        if constexpr (IsResoMicrotrack) {
          histos.fill(HIST("QAafter/kaonDCAxy"), kaon.pt(), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(kaon.trackSelectionFlags()));
          histos.fill(HIST("QAafter/kaonDCAz"), kaon.pt(), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(kaon.trackSelectionFlags()));
          histos.fill(HIST("QAafter/kaonTPCNSigma"), kaon.pt(), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(kaon.pidNSigmaKaFlag()));
          if (kaon.hasTOF()) {
            histos.fill(HIST("QAafter/kaonTOFNSigma"), kaon.pt(), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(kaon.pidNSigmaKaFlag()));
          }
        } else {
          histos.fill(HIST("QAafter/kaonDCAxy"), kaon.pt(), kaon.dcaXY());
          histos.fill(HIST("QAafter/kaonDCAz"), kaon.pt(), kaon.dcaZ());
          histos.fill(HIST("QAafter/kaonTPCNSigma"), kaon.pt(), kaon.tpcNSigmaKa());
          if (kaon.hasTOF()) {
            histos.fill(HIST("QAafter/kaonTOFNSigma"), kaon.pt(), kaon.tofNSigmaKa());
          }
          if constexpr (requires { kaon.tpcNClsFound(); }) {
            histos.fill(HIST("QAafter/kaonTPCNcls"), kaon.tpcNClsFound());
          }
          if constexpr (requires { kaon.itsNCls(); }) {
            histos.fill(HIST("QAafter/kaonITSNcls"), kaon.itsNCls());
          }
        }
      }

      int kaonCharge = kaon.sign();

      // Loop over V0 candidates
      for (const auto& v0 : v0s) {
        // Lambda QA before cuts
        if constexpr (!IsMix) {
          histos.fill(HIST("QAbefore/lambdaMass"), v0.mLambda());
          histos.fill(HIST("QAbefore/lambdaPt"), v0.pt());
          histos.fill(HIST("QAbefore/lambdaEta"), v0.eta());
          histos.fill(HIST("QAbefore/lambdaCosPA"), v0.pt(), v0.v0CosPA());
          histos.fill(HIST("QAbefore/lambdaRadius"), v0.pt(), v0.transRadius());
          histos.fill(HIST("QAbefore/lambdaDauDCA"), v0.pt(), v0.daughDCA());
          histos.fill(HIST("QAbefore/lambdaDauPosDCA"), v0.pt(), std::abs(v0.dcapostopv()));
          histos.fill(HIST("QAbefore/lambdaDauNegDCA"), v0.pt(), std::abs(v0.dcanegtopv()));

          // Calculate proper lifetime manually
          float dx = v0.decayVtxX() - collision.posX();
          float dy = v0.decayVtxY() - collision.posY();
          float dz = v0.decayVtxZ() - collision.posZ();
          float l = std::sqrt(dx * dx + dy * dy + dz * dz);
          float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
          auto properLifetime = (l / (p + kSmallMomentumDenominator)) * MassLambda;
          histos.fill(HIST("QAbefore/lambdaProperLifetime"), v0.pt(), properLifetime);
        }

        // Try Lambda
        bool isLambda = v0Cut(collision, v0, true);
        // Try Anti-Lambda
        bool isAntiLambda = v0Cut(collision, v0, false);

        if (!isLambda && !isAntiLambda)
          continue;

        if constexpr (!IsMix) {
          nV0sAfterCuts++;
          // QA after cuts (fill for whichever passes)
          if (isLambda) {
            histos.fill(HIST("QAafter/lambdaMass"), v0.mLambda());
          }
          if (isAntiLambda) {
            histos.fill(HIST("QAafter/lambdaMass"), v0.mAntiLambda());
          }
          histos.fill(HIST("QAafter/lambdaPt"), v0.pt());
          histos.fill(HIST("QAafter/lambdaEta"), v0.eta());
          histos.fill(HIST("QAafter/lambdaCosPA"), v0.pt(), v0.v0CosPA());
          histos.fill(HIST("QAafter/lambdaRadius"), v0.pt(), v0.transRadius());
          histos.fill(HIST("QAafter/lambdaDauDCA"), v0.pt(), v0.daughDCA());
          histos.fill(HIST("QAafter/lambdaDauPosDCA"), v0.pt(), std::abs(v0.dcapostopv()));
          histos.fill(HIST("QAafter/lambdaDauNegDCA"), v0.pt(), std::abs(v0.dcanegtopv()));

          float dx = v0.decayVtxX() - collision.posX();
          float dy = v0.decayVtxY() - collision.posY();
          float dz = v0.decayVtxZ() - collision.posZ();
          float l = std::sqrt(dx * dx + dy * dy + dz * dz);
          float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
          auto properLifetime = (l / (p + kSmallMomentumDenominator)) * MassLambda;
          histos.fill(HIST("QAafter/lambdaProperLifetime"), v0.pt(), properLifetime);
        }

        // Build 4 combinations
        ROOT::Math::PxPyPzEVector pKaon, pLambda, pRes;
        pKaon = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(kaon.pt(), kaon.eta(), kaon.phi(), MassKaonCharged));

        // K+ Lambda
        if (kaonCharge > 0 && isLambda) {
          pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mLambda()));
          pRes = pKaon + pLambda;
          if constexpr (!IsMix) {
            histos.fill(HIST("xi1820/kplus_lambda/hInvMassKplusLambda"), pRes.M());
            histos.fill(HIST("xi1820/kplus_lambda/hMassPtCentKplusLambda"), pRes.M(), pRes.Pt(), cent);
          } else {
            histos.fill(HIST("xi1820/kplus_lambda/hInvMassKplusLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kplus_lambda/hMassPtCentKplusLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }
        }

        // K+ Anti-Lambda
        if (kaonCharge > 0 && isAntiLambda) {
          pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda()));
          pRes = pKaon + pLambda;
          if constexpr (!IsMix) {
            histos.fill(HIST("xi1820/kplus_antilambda/hInvMassKplusAntiLambda"), pRes.M());
            histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda"), pRes.M(), pRes.Pt(), cent);
          } else {
            histos.fill(HIST("xi1820/kplus_antilambda/hInvMassKplusAntiLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }
        }

        // K- Lambda
        if (kaonCharge < 0 && isLambda) {
          pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mLambda()));
          pRes = pKaon + pLambda;
          if constexpr (!IsMix) {
            histos.fill(HIST("xi1820/kminus_lambda/hInvMassKminusLambda"), pRes.M());
            histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentKminusLambda"), pRes.M(), pRes.Pt(), cent);
          } else {
            histos.fill(HIST("xi1820/kminus_lambda/hInvMassKminusLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentKminusLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }
        }

        // K- Anti-Lambda
        if (kaonCharge < 0 && isAntiLambda) {
          pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda()));
          pRes = pKaon + pLambda;
          if constexpr (!IsMix) {
            histos.fill(HIST("xi1820/kminus_antilambda/hInvMassKminusAntiLambda"), pRes.M());
            histos.fill(HIST("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda"), pRes.M(), pRes.Pt(), cent);
          } else {
            histos.fill(HIST("xi1820/kminus_antilambda/hInvMassKminusAntiLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }
        }
      }
    }

    // Fill event QA for after-cuts counters (only for same-event)
    if constexpr (!IsMix) {
      histos.fill(HIST("Event/nV0sAfterCuts"), nV0sAfterCuts);
      histos.fill(HIST("Event/nKaonsAfterCuts"), nKaonsAfterCuts);
    }
  }

  void processDummy(const aod::ResoCollision& /*collision*/)
  {
    // Dummy function to satisfy the compiler
  }
  PROCESS_SWITCH(Xi1820Analysis, processDummy, "Process Dummy", true);

  void processDataWithTracks(const aod::ResoCollision& collision,
                             aod::ResoV0s const& resov0s,
                             aod::ResoTracks const& resotracks)
  {
    fill<false, false>(collision, resov0s, resotracks);
  }
  PROCESS_SWITCH(Xi1820Analysis, processDataWithTracks, "Process Event with ResoTracks", false);

  void processDataWithMicroTracks(const aod::ResoCollision& collision,
                                  aod::ResoV0s const& resov0s,
                                  aod::ResoMicroTracks const& resomicrotracks)
  {
    fill<false, true>(collision, resov0s, resomicrotracks);
  }
  PROCESS_SWITCH(Xi1820Analysis, processDataWithMicroTracks, "Process Event with ResoMicroTracks", false);

  void processMixedEventWithTracks(const aod::ResoCollisions& collisions,
                                   aod::ResoV0s const& resov0s,
                                   aod::ResoTracks const& resotracks)
  {

    auto v0sTracksTuple = std::make_tuple(resov0s, resotracks);
    Pair<aod::ResoCollisions, aod::ResoV0s, aod::ResoTracks, BinningTypeVertexContributor> pairs{colBinning, nEvtMixing, -1, collisions, v0sTracksTuple, &cache};

    for (auto& [collision1, v0s1, collision2, tracks2] : pairs) { // o2-linter: disable=const-ref-in-for-loop (structured bindings from Pair iterator cannot be const)
      auto cent = collision1.cent();

      for (const auto& kaon : tracks2) {
        if (!kaonCut<false>(kaon))
          continue;
        int kaonCharge = kaon.sign();

        for (const auto& v0 : v0s1) {
          bool isLambda = v0Cut(collision1, v0, true);
          bool isAntiLambda = v0Cut(collision1, v0, false);

          if (!isLambda && !isAntiLambda)
            continue;

          ROOT::Math::PxPyPzEVector pKaon, pLambda, pRes;
          pKaon = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(kaon.pt(), kaon.eta(), kaon.phi(), MassKaonCharged));

          // K+ Lambda
          if (kaonCharge > 0 && isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mLambda()));
            pRes = pKaon + pLambda;
            histos.fill(HIST("xi1820/kplus_lambda/hInvMassKplusLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kplus_lambda/hMassPtCentKplusLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }

          // K+ Anti-Lambda
          if (kaonCharge > 0 && isAntiLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda()));
            pRes = pKaon + pLambda;
            histos.fill(HIST("xi1820/kplus_antilambda/hInvMassKplusAntiLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }

          // K- Lambda
          if (kaonCharge < 0 && isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mLambda()));
            pRes = pKaon + pLambda;
            histos.fill(HIST("xi1820/kminus_lambda/hInvMassKminusLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentKminusLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }

          // K- Anti-Lambda
          if (kaonCharge < 0 && isAntiLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda()));
            pRes = pKaon + pLambda;
            histos.fill(HIST("xi1820/kminus_antilambda/hInvMassKminusAntiLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processMixedEventWithTracks, "Process Mixed Event with ResoTracks", false);

  void processMixedEventWithMicroTracks(const aod::ResoCollisions& collisions,
                                        aod::ResoV0s const& resov0s,
                                        aod::ResoMicroTracks const& resomicrotracks)
  {

    auto v0sTracksTuple = std::make_tuple(resov0s, resomicrotracks);
    Pair<aod::ResoCollisions, aod::ResoV0s, aod::ResoMicroTracks, BinningTypeVertexContributor> pairs{colBinning, nEvtMixing, -1, collisions, v0sTracksTuple, &cache};

    for (auto& [collision1, v0s1, collision2, tracks2] : pairs) { // o2-linter: disable=const-ref-in-for-loop (structured bindings from Pair iterator cannot be const)
      auto cent = collision1.cent();

      for (const auto& kaon : tracks2) {
        if (!kaonCut<true>(kaon))
          continue;
        int kaonCharge = kaon.sign();

        for (const auto& v0 : v0s1) {
          bool isLambda = v0Cut(collision1, v0, true);
          bool isAntiLambda = v0Cut(collision1, v0, false);

          if (!isLambda && !isAntiLambda)
            continue;

          ROOT::Math::PxPyPzEVector pKaon, pLambda, pRes;
          pKaon = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(kaon.pt(), kaon.eta(), kaon.phi(), MassKaonCharged));

          // K+ Lambda
          if (kaonCharge > 0 && isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mLambda()));
            pRes = pKaon + pLambda;
            histos.fill(HIST("xi1820/kplus_lambda/hInvMassKplusLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kplus_lambda/hMassPtCentKplusLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }

          // K+ Anti-Lambda
          if (kaonCharge > 0 && isAntiLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda()));
            pRes = pKaon + pLambda;
            histos.fill(HIST("xi1820/kplus_antilambda/hInvMassKplusAntiLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kplus_antilambda/hMassPtCentKplusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }

          // K- Lambda
          if (kaonCharge < 0 && isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mLambda()));
            pRes = pKaon + pLambda;
            histos.fill(HIST("xi1820/kminus_lambda/hInvMassKminusLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kminus_lambda/hMassPtCentKminusLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }

          // K- Anti-Lambda
          if (kaonCharge < 0 && isAntiLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda()));
            pRes = pKaon + pLambda;
            histos.fill(HIST("xi1820/kminus_antilambda/hInvMassKminusAntiLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/kminus_antilambda/hMassPtCentKminusAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processMixedEventWithMicroTracks, "Process Mixed Event with ResoMicroTracks", false);

  // K0s + Lambda analysis
  void processK0sLambda(const aod::ResoCollision& collision,
                        aod::ResoV0s const& resov0s)
  {
    auto cent = collision.cent();

    // Fill event QA histograms
    histos.fill(HIST("Event/posZ"), collision.posZ());
    histos.fill(HIST("Event/centrality"), cent);
    histos.fill(HIST("Event/posZvsCent"), collision.posZ(), cent);
    histos.fill(HIST("Event/nV0s"), resov0s.size());

    // Loop over V0s for K0s
    for (const auto& k0s : resov0s) {
      // K0s QA before cuts
      histos.fill(HIST("QAbefore/k0sMass"), k0s.mK0Short());
      histos.fill(HIST("QAbefore/k0sPt"), k0s.pt());
      histos.fill(HIST("QAbefore/k0sEta"), k0s.eta());
      histos.fill(HIST("QAbefore/k0sCosPA"), k0s.pt(), k0s.v0CosPA());
      histos.fill(HIST("QAbefore/k0sRadius"), k0s.pt(), k0s.transRadius());
      histos.fill(HIST("QAbefore/k0sDauDCA"), k0s.pt(), k0s.daughDCA());

      float dx = k0s.decayVtxX() - collision.posX();
      float dy = k0s.decayVtxY() - collision.posY();
      float dz = k0s.decayVtxZ() - collision.posZ();
      float l = std::sqrt(dx * dx + dy * dy + dz * dz);
      float p = std::sqrt(k0s.px() * k0s.px() + k0s.py() * k0s.py() + k0s.pz() * k0s.pz());
      auto k0sProperLifetime = (l / (p + 1e-10)) * MassK0Short;
      histos.fill(HIST("QAbefore/k0sProperLifetime"), k0s.pt(), k0sProperLifetime);

      if (!k0sCut(collision, k0s))
        continue;

      // K0s QA after cuts
      histos.fill(HIST("QAafter/k0sMass"), k0s.mK0Short());
      histos.fill(HIST("QAafter/k0sPt"), k0s.pt());
      histos.fill(HIST("QAafter/k0sEta"), k0s.eta());
      histos.fill(HIST("QAafter/k0sCosPA"), k0s.pt(), k0s.v0CosPA());
      histos.fill(HIST("QAafter/k0sRadius"), k0s.pt(), k0s.transRadius());
      histos.fill(HIST("QAafter/k0sDauDCA"), k0s.pt(), k0s.daughDCA());
      histos.fill(HIST("QAafter/k0sProperLifetime"), k0s.pt(), k0sProperLifetime);

      // Loop over V0s for Lambda
      for (const auto& lambda : resov0s) {
        // Try Lambda
        bool isLambda = v0Cut(collision, lambda, true);
        // Try Anti-Lambda
        bool isAntiLambda = v0Cut(collision, lambda, false);

        if (!isLambda && !isAntiLambda)
          continue;

        // 4-vectors
        ROOT::Math::PxPyPzEVector pK0s, pLambda, pRes;
        pK0s = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(k0s.pt(), k0s.eta(), k0s.phi(), MassK0Short));

        // K0s + Lambda
        if (isLambda) {
          pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(lambda.pt(), lambda.eta(), lambda.phi(), lambda.mLambda()));
          pRes = pK0s + pLambda;
          histos.fill(HIST("xi1820/k0s_lambda/hInvMassK0sLambda"), pRes.M());
          histos.fill(HIST("xi1820/k0s_lambda/hMassPtCentK0sLambda"), pRes.M(), pRes.Pt(), cent);
        }

        // K0s + Anti-Lambda
        if (isAntiLambda) {
          pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(lambda.pt(), lambda.eta(), lambda.phi(), lambda.mAntiLambda()));
          pRes = pK0s + pLambda;
          histos.fill(HIST("xi1820/k0s_antilambda/hInvMassK0sAntiLambda"), pRes.M());
          histos.fill(HIST("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda"), pRes.M(), pRes.Pt(), cent);
        }
      }
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processK0sLambda, "Process K0s + Lambda", false);

  // K0s + Lambda mixed event analysis
  void processK0sLambdaMixedEvent(const aod::ResoCollisions& collisions,
                                  aod::ResoV0s const& resov0s)
  {

    auto v0sV0sTuple = std::make_tuple(resov0s, resov0s);
    Pair<aod::ResoCollisions, aod::ResoV0s, aod::ResoV0s, BinningTypeVertexContributor> pairs{colBinning, nEvtMixing, -1, collisions, v0sV0sTuple, &cache};

    for (auto& [collision1, k0s1, collision2, lambda2] : pairs) { // o2-linter: disable=const-ref-in-for-loop (structured bindings from Pair iterator cannot be const)
      auto cent = collision1.cent();

      for (const auto& k0s : k0s1) {
        if (!k0sCut(collision1, k0s))
          continue;

        for (const auto& lambda : lambda2) {
          bool isLambda = v0Cut(collision2, lambda, true);
          bool isAntiLambda = v0Cut(collision2, lambda, false);

          if (!isLambda && !isAntiLambda)
            continue;

          ROOT::Math::PxPyPzEVector pK0s, pLambda, pRes;
          pK0s = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(k0s.pt(), k0s.eta(), k0s.phi(), MassK0Short));

          // K0s + Lambda
          if (isLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(lambda.pt(), lambda.eta(), lambda.phi(), lambda.mLambda()));
            pRes = pK0s + pLambda;
            histos.fill(HIST("xi1820/k0s_lambda/hInvMassK0sLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/k0s_lambda/hMassPtCentK0sLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }

          // K0s + Anti-Lambda
          if (isAntiLambda) {
            pLambda = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(lambda.pt(), lambda.eta(), lambda.phi(), lambda.mAntiLambda()));
            pRes = pK0s + pLambda;
            histos.fill(HIST("xi1820/k0s_antilambda/hInvMassK0sAntiLambda_Mix"), pRes.M());
            histos.fill(HIST("xi1820/k0s_antilambda/hMassPtCentK0sAntiLambda_Mix"), pRes.M(), pRes.Pt(), cent);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processK0sLambdaMixedEvent, "Process K0s + Lambda Mixed Event", false);

  // MC processes - placeholder for future implementation
  void processMCWithTracks(const aod::ResoCollision& /*collision*/,
                           aod::ResoV0s const& /*resov0s*/,
                           aod::ResoTracks const& /*resotracks*/,
                           aod::McParticles const& /*mcParticles*/)
  {
    // TODO: Implement MC truth matching for K± + Lambda
    // - Match reconstructed kaons to MC kaons
    // - Match reconstructed Lambdas to MC Lambdas
    // - Fill MC truth histograms
    // - Fill reconstruction efficiency histograms
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCWithTracks, "Process MC with ResoTracks (placeholder)", false);

  void processMCWithMicroTracks(const aod::ResoCollision& /*collision*/,
                                aod::ResoV0s const& /*resov0s*/,
                                aod::ResoMicroTracks const& /*resomicrotracks*/,
                                aod::McParticles const& /*mcParticles*/)
  {
    // TODO: Implement MC truth matching for K± + Lambda with MicroTracks
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCWithMicroTracks, "Process MC with ResoMicroTracks (placeholder)", false);

  void processMCK0sLambda(const aod::ResoCollision& /*collision*/,
                          aod::ResoV0s const& /*resov0s*/,
                          aod::McParticles const& /*mcParticles*/)
  {
    // TODO: Implement MC truth matching for K0s + Lambda
    // - Match reconstructed K0s to MC K0s
    // - Match reconstructed Lambdas to MC Lambdas
    // - Fill MC truth histograms
    // - Fill reconstruction efficiency histograms
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCK0sLambda, "Process MC K0s + Lambda (placeholder)", false);

  void processMCGenerated(aod::McParticles const& mcParticles)
  {
    // Process MC generated particles (no reconstruction requirement)
    // Xi(1820)0 PDG code: 123314 (neutral, decays to K+ Lambda or K0s Lambda)
    // Note: PDG doesn't have separate codes for charge states in this case

    for (const auto& mcParticle : mcParticles) {
      // Look for Xi(1820) - PDG code can vary, check for resonance mass ~1820 MeV
      int pdg = mcParticle.pdgCode();

      // Xi(1820)0: PDG 123314
      // Check if it's Xi(1820) or similar resonance
      if (std::abs(pdg) != kPdgXi1820)
        continue;

      // Fill generated level histograms
      auto pt = mcParticle.pt();
      auto eta = mcParticle.eta();
      auto y = mcParticle.y();

      histos.fill(HIST("MC/hMCGenXi1820Pt"), pt);
      histos.fill(HIST("MC/hMCGenXi1820PtEta"), pt, eta);
      histos.fill(HIST("MC/hMCGenXi1820Y"), y);

      // Get daughters
      auto daughters = mcParticle.daughters_as<aod::McParticles>();
      if (daughters.size() != kExpectedDaughters)
        continue;

      int daughter1PDG = 0, daughter2PDG = 0;
      ROOT::Math::PxPyPzEVector p1, p2, pMother;

      int iDaughter = 0;
      for (const auto& daughter : daughters) {
        if (iDaughter == 0) {
          daughter1PDG = daughter.pdgCode();
          p1.SetPxPyPzE(daughter.px(), daughter.py(), daughter.pz(), daughter.e());
        } else {
          daughter2PDG = daughter.pdgCode();
          p2.SetPxPyPzE(daughter.px(), daughter.py(), daughter.pz(), daughter.e());
        }
        iDaughter++;
      }

      pMother = p1 + p2;

      // Check decay channels
      auto motherPt = pMother.Pt();
      auto motherM = pMother.M();

      // K+ + Lambda
      if ((daughter1PDG == PDG_t::kKPlus && daughter2PDG == PDG_t::kLambda0) ||
          (daughter1PDG == PDG_t::kLambda0 && daughter2PDG == PDG_t::kKPlus)) {
        histos.fill(HIST("MC/hMCTruthInvMassKplusLambda"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtKplusLambda"), motherM, motherPt);
      }

      // K- + Anti-Lambda
      if ((daughter1PDG == PDG_t::kKMinus && daughter2PDG == PDG_t::kLambda0Bar) ||
          (daughter1PDG == PDG_t::kLambda0Bar && daughter2PDG == PDG_t::kKMinus)) {
        histos.fill(HIST("MC/hMCTruthInvMassKminusAntiLambda"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtKminusAntiLambda"), motherM, motherPt);
      }

      // K0s + Lambda
      if ((daughter1PDG == PDG_t::kK0Short && daughter2PDG == PDG_t::kLambda0) ||
          (daughter1PDG == PDG_t::kLambda0 && daughter2PDG == PDG_t::kK0Short)) {
        histos.fill(HIST("MC/hMCTruthInvMassK0sLambda"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtK0sLambda"), motherM, motherPt);
      }

      // K0s + Anti-Lambda
      if ((daughter1PDG == PDG_t::kK0Short && daughter2PDG == PDG_t::kLambda0Bar) ||
          (daughter1PDG == PDG_t::kLambda0Bar && daughter2PDG == PDG_t::kK0Short)) {
        histos.fill(HIST("MC/hMCTruthInvMassK0sAntiLambda"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtK0sAntiLambda"), motherM, motherPt);
      }
    }
  }
  PROCESS_SWITCH(Xi1820Analysis, processMCGenerated, "Process MC generated particles", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Xi1820Analysis>(cfgc)};
}
