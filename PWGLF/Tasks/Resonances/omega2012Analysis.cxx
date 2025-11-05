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

/// \file omega2012Analysis.cxx
/// \brief Invariant Mass Reconstruction of Omega(2012) Resonance
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

#include <set>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct Omega2012Analysis {
  // Constants
  static constexpr float kSmallNumber = 1e-10f;       // Small number to avoid division by zero
  static constexpr float kMaxDCAV0ToPV = 1.0f;        // Maximum DCA of V0 to PV
  static constexpr int kNumExpectedDaughters = 2;     // Expected number of daughters for 2-body decay
  static constexpr int kPlaceholderPdgCode = 9999999; // o2-linter: disable=pdg/explicit-code (placeholder for generator-specific Omega(2012) PDG code)
  SliceCache cache;
  Preslice<aod::ResoCascades> perResoCollisionCasc = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoV0s> perResoCollisionV0 = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoTracks> perResoCollisionTrack = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoMicroTracks> perResoCollisionMicroTrack = aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Axes
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0}, "pT"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pT (QA)"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Centrality"};

  // Invariant mass range for Omega(2012) → Xi + K0s
  Configurable<float> cInvMassStart{"cInvMassStart", 1.6, "Invariant mass start (GeV/c^2)"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 2.2, "Invariant mass end (GeV/c^2)"};
  Configurable<int> cInvMassBins{"cInvMassBins", 600, "Invariant mass bins"};

  // Basic pre-selections (mirroring refs)
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Minimum pT for candidates"};
  Configurable<float> cMaxEtaCut{"cMaxEtaCut", 0.8, "Maximum |eta|"};
  // V0 selections (K0s) from k892pmanalysis
  Configurable<double> cV0MinCosPA{"cV0MinCosPA", 0.97, "V0 minimum pointing angle cosine"};
  Configurable<double> cV0MaxDaughDCA{"cV0MaxDaughDCA", 1.0, "V0 daughter DCA Maximum"};
  Configurable<double> cV0MassWindow{"cV0MassWindow", 0.0043, "Mass window for competing Lambda0 rejection"};
  Configurable<double> cMaxV0Etacut{"cMaxV0Etacut", 0.8, "V0 maximum eta cut"};

  // Xi (cascade) selections from xi1530Analysisqa.cxx
  Configurable<float> cDCAxyToPVByPtCascP0{"cDCAxyToPVByPtCascP0", 999., "Cascade DCAxy p0"};
  Configurable<float> cDCAxyToPVByPtCascExp{"cDCAxyToPVByPtCascExp", 1., "Cascade DCAxy exp"};
  Configurable<bool> cDCAxyToPVAsPtForCasc{"cDCAxyToPVAsPtForCasc", true, "Use pt-dep DCAxy cut (casc)"};

  Configurable<bool> cDCAzToPVAsPtForCasc{"cDCAzToPVAsPtForCasc", true, "Use pt-dep DCAz cut (casc)"};

  // V0 topology inside cascade (Λ)
  Configurable<float> cDCALambdaDaugtherscut{"cDCALambdaDaugtherscut", 0.7, "Λ daughters DCA cut"};
  Configurable<float> cDCALambdaToPVcut{"cDCALambdaToPVcut", 0.02, "Λ DCA to PV min"};
  Configurable<float> cDCAPionToPVcut{"cDCAPionToPVcut", 0.06, "π DCA to PV min"};
  Configurable<float> cDCAProtonToPVcut{"cDCAProtonToPVcut", 0.07, "p DCA to PV min"};
  Configurable<float> cV0CosPACutPtDepP0{"cV0CosPACutPtDepP0", 0.25, "V0 CosPA p0"};
  Configurable<float> cV0CosPACutPtDepP1{"cV0CosPACutPtDepP1", 0.022, "V0 CosPA p1"};
  Configurable<float> cMaxV0radiuscut{"cMaxV0radiuscut", 200., "V0 radius max"};
  Configurable<float> cMinV0radiuscut{"cMinV0radiuscut", 2.5, "V0 radius min"};
  Configurable<float> cMasswindowV0cut{"cMasswindowV0cut", 0.005, "Λ mass window for cascade V0"};

  // Cascade topology
  Configurable<float> cDCABachlorToPVcut{"cDCABachlorToPVcut", 0.06, "Bachelor DCA to PV min"};
  Configurable<float> cDCAXiDaugthersCutPtRangeLower{"cDCAXiDaugthersCutPtRangeLower", 1., "Xi pt low boundary"};
  Configurable<float> cDCAXiDaugthersCutPtRangeUpper{"cDCAXiDaugthersCutPtRangeUpper", 4., "Xi pt high boundary"};
  Configurable<float> cDCAXiDaugthersCutPtDepLower{"cDCAXiDaugthersCutPtDepLower", 0.8, "Xi daugh DCA (pt<low)"};
  Configurable<float> cDCAXiDaugthersCutPtDepMiddle{"cDCAXiDaugthersCutPtDepMiddle", 0.5, "Xi daugh DCA (low<=pt<high)"};
  Configurable<float> cDCAXiDaugthersCutPtDepUpper{"cDCAXiDaugthersCutPtDepUpper", 0.2, "Xi daugh DCA (pt>=high)"};
  Configurable<float> cCosPACascCutPtDepP0{"cCosPACascCutPtDepP0", 0.2, "Cascade CosPA p0"};
  Configurable<float> cCosPACascCutPtDepP1{"cCosPACascCutPtDepP1", 0.022, "Cascade CosPA p1"};
  Configurable<float> cMaxCascradiuscut{"cMaxCascradiuscut", 200., "Cascade radius max"};
  Configurable<float> cMinCascradiuscut{"cMinCascradiuscut", 1.1, "Cascade radius min"};
  Configurable<float> cMasswindowCasccut{"cMasswindowCasccut", 0.008, "Xi mass window"};
  Configurable<float> cMassXiminus{"cMassXiminus", 1.32171, "Xi mass (GeV/c^2)"}; // PDG

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - centrality"};

  // Enhanced K0s selections
  Configurable<float> cK0sProperLifetimeMax{"cK0sProperLifetimeMax", 20.0, "K0s proper lifetime max (cm/c)"};
  Configurable<float> cK0sArmenterosQtMin{"cK0sArmenterosQtMin", 0.0, "K0s Armenteros qt min"};
  Configurable<float> cK0sArmenterosAlphaMax{"cK0sArmenterosAlphaMax", 0.8, "K0s Armenteros alpha max"};
  Configurable<float> cK0sDauPosDCAtoPVMin{"cK0sDauPosDCAtoPVMin", 0.1, "K0s positive daughter DCA to PV min"};
  Configurable<float> cK0sDauNegDCAtoPVMin{"cK0sDauNegDCAtoPVMin", 0.1, "K0s negative daughter DCA to PV min"};
  Configurable<float> cK0sRadiusMin{"cK0sRadiusMin", 0.5, "K0s decay radius min"};
  Configurable<float> cK0sRadiusMax{"cK0sRadiusMax", 100.0, "K0s decay radius max"};
  Configurable<bool> cK0sCrossMassRejection{"cK0sCrossMassRejection", true, "Enable Lambda mass rejection for K0s"};

  // Pion track selections for 3-body decay
  Configurable<float> cPionPtMin{"cPionPtMin", 0.15, "Minimum pion pT"};
  Configurable<float> cPionEtaMax{"cPionEtaMax", 0.8, "Maximum pion |eta|"};
  Configurable<float> cPionDCAxyMax{"cPionDCAxyMax", 0.1, "Maximum pion DCAxy to PV"};
  Configurable<float> cPionDCAzMax{"cPionDCAzMax", 0.2, "Maximum pion DCAz to PV"};
  Configurable<int> cPionTPCNClusMin{"cPionTPCNClusMin", 70, "Minimum TPC clusters for pion"};

  // Pion PID selections
  Configurable<float> cPionTPCNSigmaMax{"cPionTPCNSigmaMax", 3.0, "Maximum TPC NSigma for pion"};
  Configurable<float> cPionTOFNSigmaMax{"cPionTOFNSigmaMax", 3.0, "Maximum TOF NSigma for pion"};
  Configurable<bool> cPionUsePtDepPID{"cPionUsePtDepPID", false, "Use pT-dependent PID cuts for pion"};
  Configurable<std::vector<float>> cPionPIDPtBins{"cPionPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for pion PID cuts"};
  Configurable<std::vector<float>> cPionTPCNSigmaCuts{"cPionTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (pion)"};
  Configurable<std::vector<float>> cPionTOFNSigmaCuts{"cPionTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (pion)"};
  Configurable<std::vector<int>> cPionTOFRequired{"cPionTOFRequired", {0, 0, 1, 1}, "Require TOF per pT bin (pion)"};

  // Xi1530 mass window cut
  Configurable<float> cXi1530Mass{"cXi1530Mass", 1.53, "Xi(1530) mass (GeV/c^2)"};
  Configurable<float> cXi1530MassWindow{"cXi1530MassWindow", 0.01, "Xi(1530) mass window (GeV/c^2)"};

  // PDG masses
  double massK0 = MassK0Short;

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  BinningTypeVertexContributor colBinning{{cfgVtxBins, cfgMultBins}, true};

  void init(InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^{2})"};
    AxisSpec xiMassAxis = {400, 1.25, 1.65, "#Xi mass (GeV/#it{c}^{2})"};
    AxisSpec k0sMassAxis = {100, 0.4, 0.6, "K^{0}_{S} mass (GeV/#it{c}^{2})"};
    AxisSpec dcaAxis = {200, 0., 2.0, "DCA (cm)"};
    AxisSpec dcaxyAxis = {200, -1.0, 1.0, "DCA_{xy} (cm)"};
    AxisSpec dcazAxis = {200, -2.0, 2.0, "DCA_{z} (cm)"};
    AxisSpec cosPAAxis = {1000, 0.95, 1.0, "cos(PA)"};
    AxisSpec radiusAxis = {200, 0, 200, "Radius (cm)"};
    AxisSpec lifetimeAxis = {200, 0, 50, "Proper lifetime (cm/c)"};
    AxisSpec armQtAxis = {100, 0, 0.3, "q_{T} (GeV/c)"};
    AxisSpec armAlphaAxis = {100, -1.0, 1.0, "#alpha"};

    // Event QA histograms
    histos.add("Event/posZ", "Event vertex Z position", kTH1F, {{200, -20., 20., "V_{z} (cm)"}});
    histos.add("Event/centrality", "Event centrality distribution", kTH1F, {centAxis});
    histos.add("Event/posZvsCent", "Vertex Z vs Centrality", kTH2F, {{200, -20., 20., "V_{z} (cm)"}, centAxis});
    histos.add("Event/nCascades", "Number of cascades per event", kTH1F, {{100, 0., 100., "N_{cascades}"}});
    histos.add("Event/nV0s", "Number of V0s per event", kTH1F, {{200, 0., 200., "N_{V0s}"}});
    histos.add("Event/nCascadesAfterCuts", "Number of cascades per event after cuts", kTH1F, {{50, 0., 50., "N_{cascades}"}});
    histos.add("Event/nV0sAfterCuts", "Number of V0s per event after cuts", kTH1F, {{100, 0., 100., "N_{V0s}"}});

    // Xi QA histograms
    histos.add("QAbefore/xiMass", "Xi mass before cuts", kTH1F, {xiMassAxis});
    histos.add("QAbefore/xiPt", "Xi pT before cuts", kTH1F, {ptAxisQA});
    histos.add("QAbefore/xiEta", "Xi eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAbefore/xiDCAxy", "Xi DCAxy before cuts", kTH2F, {ptAxisQA, dcaxyAxis});
    histos.add("QAbefore/xiDCAz", "Xi DCAz before cuts", kTH2F, {ptAxisQA, dcazAxis});
    histos.add("QAbefore/xiV0CosPA", "Xi V0 CosPA before cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAbefore/xiCascCosPA", "Xi Cascade CosPA before cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAbefore/xiV0Radius", "Xi V0 radius before cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAbefore/xiCascRadius", "Xi Cascade radius before cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAbefore/xiV0DauDCA", "Xi V0 daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAbefore/xiCascDauDCA", "Xi Cascade daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});

    histos.add("QAafter/xiMass", "Xi mass after cuts", kTH1F, {xiMassAxis});
    histos.add("QAafter/xiPt", "Xi pT after cuts", kTH1F, {ptAxisQA});
    histos.add("QAafter/xiEta", "Xi eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAafter/xiDCAxy", "Xi DCAxy after cuts", kTH2F, {ptAxisQA, dcaxyAxis});
    histos.add("QAafter/xiDCAz", "Xi DCAz after cuts", kTH2F, {ptAxisQA, dcazAxis});
    histos.add("QAafter/xiV0CosPA", "Xi V0 CosPA after cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAafter/xiCascCosPA", "Xi Cascade CosPA after cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAafter/xiV0Radius", "Xi V0 radius after cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAafter/xiCascRadius", "Xi Cascade radius after cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAafter/xiV0DauDCA", "Xi V0 daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAafter/xiCascDauDCA", "Xi Cascade daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});

    // K0s QA histograms
    histos.add("QAbefore/k0sMassPt", "K0s mass vs pT before cuts", kTH2F, {ptAxisQA, k0sMassAxis});
    histos.add("QAbefore/k0sPt", "K0s pT before cuts", kTH1F, {ptAxisQA});
    histos.add("QAbefore/k0sEta", "K0s eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAbefore/k0sCosPA", "K0s CosPA before cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAbefore/k0sRadius", "K0s radius before cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAbefore/k0sDauDCA", "K0s daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAbefore/k0sDCAtoPV", "K0s DCA to PV before cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAbefore/k0sProperLifetime", "K0s proper lifetime before cuts", kTH2F, {ptAxisQA, lifetimeAxis});
    histos.add("QAbefore/k0sArmenteros", "K0s Armenteros plot before cuts", kTH2F, {armAlphaAxis, armQtAxis});
    histos.add("QAbefore/k0sDauPosDCA", "K0s positive daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAbefore/k0sDauNegDCA", "K0s negative daughter DCA before cuts", kTH2F, {ptAxisQA, dcaAxis});

    histos.add("QAafter/k0sMassPt", "K0s mass vs pT after cuts", kTH2F, {ptAxisQA, k0sMassAxis});
    histos.add("QAafter/k0sPt", "K0s pT after cuts", kTH1F, {ptAxisQA});
    histos.add("QAafter/k0sEta", "K0s eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAafter/k0sCosPA", "K0s CosPA after cuts", kTH2F, {ptAxisQA, cosPAAxis});
    histos.add("QAafter/k0sRadius", "K0s radius after cuts", kTH2F, {ptAxisQA, radiusAxis});
    histos.add("QAafter/k0sDauDCA", "K0s daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAafter/k0sDCAtoPV", "K0s DCA to PV after cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAafter/k0sProperLifetime", "K0s proper lifetime after cuts", kTH2F, {ptAxisQA, lifetimeAxis});
    histos.add("QAafter/k0sArmenteros", "K0s Armenteros plot after cuts", kTH2F, {armAlphaAxis, armQtAxis});
    histos.add("QAafter/k0sDauPosDCA", "K0s positive daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});
    histos.add("QAafter/k0sDauNegDCA", "K0s negative daughter DCA after cuts", kTH2F, {ptAxisQA, dcaAxis});

    // Resonance (2-body decay: Xi + K0s)
    histos.add("omega2012/invmass", "Invariant mass of Omega(2012) → Xi + K0s", kTH1F, {invMassAxis});
    histos.add("omega2012/invmass_Mix", "Mixed event Invariant mass of Omega(2012) → Xi + K0s", kTH1F, {invMassAxis});
    histos.add("omega2012/massPtCent", "Omega(2012) mass vs pT vs cent", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("omega2012/massPtCent_Mix", "Mixed event Omega(2012) mass vs pT vs cent", kTH3F, {invMassAxis, ptAxis, centAxis});

    // 3-body decay: Xi + pi + K0s
    histos.add("omega2012_3body/invmass", "Invariant mass of Omega(2012) → Xi + #pi + K^{0}_{S}", kTH1F, {invMassAxis});
    histos.add("omega2012_3body/massPtCent", "Omega(2012) 3-body mass vs pT vs cent", kTH3F, {invMassAxis, ptAxis, centAxis});

    // Pion QA histograms for 3-body
    AxisSpec nsigmaAxis = {100, -5.0, 5.0, "N#sigma"};
    histos.add("QAbefore/pionPt", "Pion pT before cuts", kTH1F, {ptAxisQA});
    histos.add("QAbefore/pionEta", "Pion eta before cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAbefore/pionDCAxy", "Pion DCAxy before cuts", kTH2F, {ptAxisQA, dcaxyAxis});
    histos.add("QAbefore/pionDCAz", "Pion DCAz before cuts", kTH2F, {ptAxisQA, dcazAxis});
    histos.add("QAbefore/pionTPCNcls", "Pion TPC clusters before cuts", kTH1F, {{160, 0, 160, "N_{TPC clusters}"}});
    histos.add("QAbefore/pionTPCNSigma", "Pion TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    histos.add("QAbefore/pionTOFNSigma", "Pion TOF NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});

    histos.add("QAafter/pionPt", "Pion pT after cuts", kTH1F, {ptAxisQA});
    histos.add("QAafter/pionEta", "Pion eta after cuts", kTH1F, {{100, -2.0, 2.0, "#eta"}});
    histos.add("QAafter/pionDCAxy", "Pion DCAxy after cuts", kTH2F, {ptAxisQA, dcaxyAxis});
    histos.add("QAafter/pionDCAz", "Pion DCAz after cuts", kTH2F, {ptAxisQA, dcazAxis});
    histos.add("QAafter/pionTPCNcls", "Pion TPC clusters after cuts", kTH1F, {{160, 0, 160, "N_{TPC clusters}"}});
    histos.add("QAafter/pionTPCNSigma", "Pion TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    histos.add("QAafter/pionTOFNSigma", "Pion TOF NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});

    // MC truth histograms
    AxisSpec etaAxis = {100, -2.0, 2.0, "#eta"};
    AxisSpec rapidityAxis = {100, -2.0, 2.0, "y"};

    histos.add("MC/hMCGenOmega2012Pt", "MC Generated Omega(2012) pT", kTH1F, {ptAxis});
    histos.add("MC/hMCGenOmega2012PtEta", "MC Generated Omega(2012) pT vs eta", kTH2F, {ptAxis, etaAxis});
    histos.add("MC/hMCGenOmega2012Y", "MC Generated Omega(2012) rapidity", kTH1F, {rapidityAxis});
    histos.add("MC/hMCRecOmega2012Pt", "MC Reconstructed Omega(2012) pT", kTH1F, {ptAxis});
    histos.add("MC/hMCRecOmega2012PtEta", "MC Reconstructed Omega(2012) pT vs eta", kTH2F, {ptAxis, etaAxis});

    // MC truth invariant mass (from MC particles)
    histos.add("MC/hMCTruthInvMassXiK0s", "MC Truth Inv Mass Xi + K^{0}_{S}", kTH1F, {invMassAxis});
    histos.add("MC/hMCTruthMassPtXiK0s", "MC Truth Mass vs pT Xi + K^{0}_{S}", kTH2F, {invMassAxis, ptAxis});

    // MC reconstruction efficiency
    histos.add("MC/hMCRecXiPt", "MC Reconstructed Xi pT", kTH1F, {ptAxis});
    histos.add("MC/hMCRecK0sPt", "MC Reconstructed K0s pT", kTH1F, {ptAxis});
    histos.add("MC/hMCTrueXiPt", "MC True Xi pT", kTH1F, {ptAxis});
    histos.add("MC/hMCTrueK0sPt", "MC True K0s pT", kTH1F, {ptAxis});
  }

  // Enhanced V0 selection (K0s) with detailed criteria
  template <typename CollisionType, typename V0Type>
  bool v0CutEnhanced(const CollisionType& collision, const V0Type& v0)
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

    // Enhanced selections from chk892Flow
    // Daughter DCA to PV cuts
    if (std::abs(v0.dcapostopv()) < cK0sDauPosDCAtoPVMin)
      return false;
    if (std::abs(v0.dcanegtopv()) < cK0sDauNegDCAtoPVMin)
      return false;

    // Radius cuts - use transRadius instead of v0radius
    auto radius = v0.transRadius();
    if (radius < cK0sRadiusMin || radius > cK0sRadiusMax)
      return false;

    // DCA to PV
    if (std::abs(v0.dcav0topv()) > kMaxDCAV0ToPV)
      return false; // max DCA to PV

    // Proper lifetime cut - calculate manually
    float dx = v0.decayVtxX() - collision.posX();
    float dy = v0.decayVtxY() - collision.posY();
    float dz = v0.decayVtxZ() - collision.posZ();
    float l = std::sqrt(dx * dx + dy * dy + dz * dz);
    float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
    auto properLifetime = (l / (p + kSmallNumber)) * MassK0Short;
    if (properLifetime > cK0sProperLifetimeMax)
      return false;

    // Armenteros cut - skip for now as we don't have daughter momentum info in ResoV0s table
    // If needed, this would require accessing daughter tracks separately or using alternative cuts

    // Mass window
    if (std::abs(v0.mK0Short() - MassK0Short) > cV0MassWindow)
      return false;

    // Competing V0 rejection: remove (Anti)Λ
    if (cK0sCrossMassRejection) {
      if (std::abs(v0.mLambda() - MassLambda) < cV0MassWindow)
        return false;
      if (std::abs(v0.mAntiLambda() - MassLambda) < cV0MassWindow)
        return false;
    }

    return true;
  }

  // Original V0 selection for backward compatibility
  template <typename V0Type>
  bool v0Cut(const V0Type& v0)
  {
    if (std::abs(v0.eta()) > cMaxV0Etacut)
      return false;
    if (v0.v0CosPA() < cV0MinCosPA)
      return false;
    if (v0.daughDCA() > cV0MaxDaughDCA)
      return false;
    // Competing V0 rejection: remove (Anti)Λ
    if (std::abs(v0.mLambda() - MassLambda) < cV0MassWindow)
      return false;
    if (std::abs(v0.mAntiLambda() - MassLambda) < cV0MassWindow)
      return false;
    return true;
  }

  // Helper function to find pT bin index
  int getPtBinIndex(float pt)
  {
    auto ptBins = static_cast<std::vector<float>>(cPionPIDPtBins);
    for (size_t i = 0; i < ptBins.size() - 1; i++) {
      if (pt >= ptBins[i] && pt < ptBins[i + 1]) {
        return i;
      }
    }
    return -1;
  }

  // Pion PID selection
  template <bool IsResoMicrotrack, typename TrackType>
  bool pionPidCut(const TrackType& track)
  {
    float pt = track.pt();

    if constexpr (IsResoMicrotrack) {
      // For ResoMicroTracks - decode PID from flags
      float tpcNSigma = o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(track.pidNSigmaPiFlag());
      float tofNSigma = track.hasTOF() ? o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(track.pidNSigmaPiFlag()) : 999.f;

      if (cPionUsePtDepPID) {
        int ptBin = getPtBinIndex(pt);
        if (ptBin < 0)
          return false;

        auto tpcCuts = static_cast<std::vector<float>>(cPionTPCNSigmaCuts);
        auto tofCuts = static_cast<std::vector<float>>(cPionTOFNSigmaCuts);
        auto tofRequired = static_cast<std::vector<int>>(cPionTOFRequired);

        if (ptBin >= static_cast<int>(tpcCuts.size()) ||
            ptBin >= static_cast<int>(tofCuts.size()) ||
            ptBin >= static_cast<int>(tofRequired.size())) {
          return false;
        }

        if (std::abs(tpcNSigma) >= tpcCuts[ptBin])
          return false;

        if (tofRequired[ptBin] != 0) {
          if (!track.hasTOF())
            return false;
          if (std::abs(tofNSigma) >= tofCuts[ptBin])
            return false;
        } else {
          if (track.hasTOF() && std::abs(tofNSigma) >= tofCuts[ptBin])
            return false;
        }

        return true;
      } else {
        bool tpcPass = std::abs(tpcNSigma) < cPionTPCNSigmaMax;
        bool tofPass = track.hasTOF() ? std::abs(tofNSigma) < cPionTOFNSigmaMax : true;
        return tpcPass && tofPass;
      }
    } else {
      // For ResoTracks - direct access
      float tpcNSigma = track.tpcNSigmaPi();
      float tofNSigma = track.hasTOF() ? track.tofNSigmaPi() : 999.f;

      if (cPionUsePtDepPID) {
        int ptBin = getPtBinIndex(pt);
        if (ptBin < 0)
          return false;

        auto tpcCuts = static_cast<std::vector<float>>(cPionTPCNSigmaCuts);
        auto tofCuts = static_cast<std::vector<float>>(cPionTOFNSigmaCuts);
        auto tofRequired = static_cast<std::vector<int>>(cPionTOFRequired);

        if (ptBin >= static_cast<int>(tpcCuts.size()) ||
            ptBin >= static_cast<int>(tofCuts.size()) ||
            ptBin >= static_cast<int>(tofRequired.size())) {
          return false;
        }

        if (std::abs(tpcNSigma) >= tpcCuts[ptBin])
          return false;

        if (tofRequired[ptBin] != 0) {
          if (!track.hasTOF())
            return false;
          if (std::abs(tofNSigma) >= tofCuts[ptBin])
            return false;
        } else {
          if (track.hasTOF() && std::abs(tofNSigma) >= tofCuts[ptBin])
            return false;
        }

        return true;
      } else {
        bool tpcPass = std::abs(tpcNSigma) < cPionTPCNSigmaMax;
        bool tofPass = track.hasTOF() ? std::abs(tofNSigma) < cPionTOFNSigmaMax : true;
        return tpcPass && tofPass;
      }
    }
  }

  // Pion track selection (for both ResoTracks and ResoMicroTracks)
  template <bool IsResoMicrotrack, typename TrackType>
  bool pionCut(const TrackType& track)
  {
    // Basic kinematic cuts
    if (track.pt() < cPionPtMin)
      return false;
    if (std::abs(track.eta()) > cPionEtaMax)
      return false;

    // DCA cuts - different access for ResoMicroTracks
    if constexpr (IsResoMicrotrack) {
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags()) > cPionDCAxyMax)
        return false;
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags()) > cPionDCAzMax)
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cPionDCAxyMax)
        return false;
      if (std::abs(track.dcaZ()) > cPionDCAzMax)
        return false;
    }

    // Track quality cuts - only for ResoTracks
    if constexpr (!IsResoMicrotrack) {
      if constexpr (requires { track.tpcNClsFound(); }) {
        if (track.tpcNClsFound() < cPionTPCNClusMin)
          return false;
      }
    }

    // PID selection
    if (!pionPidCut<IsResoMicrotrack>(track))
      return false;

    return true;
  }

  // Xi1530 mass window cut
  template <typename XiType, typename PionType>
  bool xi1530MassCut(const XiType& xi, const PionType& pion)
  {
    // Calculate Xi + pion invariant mass
    ROOT::Math::PxPyPzEVector pXi, pPion, pXi1530;
    pXi = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(xi.pt(), xi.eta(), xi.phi(), xi.mXi()));
    pPion = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(pion.pt(), pion.eta(), pion.phi(), MassPionCharged));
    pXi1530 = pXi + pPion;

    // Check if mass is within Xi(1530) window
    float massDiff = std::abs(pXi1530.M() - cXi1530Mass);
    return massDiff < cXi1530MassWindow;
  }

  // Primary-level cascade kinematics
  template <typename CascT>
  bool cascprimaryTrackCut(const CascT& c)
  {
    if (std::abs(c.eta()) > cMaxEtaCut)
      return false;
    if (std::abs(c.pt()) < cMinPtcut)
      return false;
    if (cDCAxyToPVAsPtForCasc) {
      if (std::abs(c.dcaXYCascToPV()) > (cDCAxyToPVByPtCascP0 + cDCAxyToPVByPtCascExp * c.pt()))
        return false;
    }
    if (cDCAzToPVAsPtForCasc) {
      if (std::abs(c.dcaZCascToPV()) > (cDCAxyToPVByPtCascP0 + cDCAxyToPVByPtCascExp * std::pow(c.pt(), -1.1f)))
        return false;
    }
    return true;
  }

  // Cascade topological selections adapted from xi1530Analysisqa
  template <typename CascT>
  bool casctopCut(const CascT& c)
  {
    // V0 (Λ) topology inside cascade
    if (std::abs(c.daughDCA()) > cDCALambdaDaugtherscut)
      return false;
    if (std::abs(c.dcav0topv()) < cDCALambdaToPVcut)
      return false;

    if (c.sign() < 0) { // Xi-
      if (std::abs(c.dcanegtopv()) < cDCAPionToPVcut)
        return false;
      if (std::abs(c.dcapostopv()) < cDCAProtonToPVcut)
        return false;
    } else { // Anti-Xi
      if (std::abs(c.dcanegtopv()) < cDCAProtonToPVcut)
        return false;
      if (std::abs(c.dcapostopv()) < cDCAPionToPVcut)
        return false;
    }

    if (c.v0CosPA() < std::cos(cV0CosPACutPtDepP0 - cV0CosPACutPtDepP1 * c.pt()))
      return false;
    if (c.transRadius() > cMaxV0radiuscut || c.transRadius() < cMinV0radiuscut)
      return false;
    if (std::abs(c.mLambda() - MassLambda) > cMasswindowV0cut)
      return false;

    // Cascade-level topology
    if (std::abs(c.dcabachtopv()) < cDCABachlorToPVcut)
      return false;

    if (c.pt() < cDCAXiDaugthersCutPtRangeLower) {
      if (c.cascDaughDCA() > cDCAXiDaugthersCutPtDepLower)
        return false;
    } else if (c.pt() < cDCAXiDaugthersCutPtRangeUpper) {
      if (c.cascDaughDCA() > cDCAXiDaugthersCutPtDepMiddle)
        return false;
    } else {
      if (c.cascDaughDCA() > cDCAXiDaugthersCutPtDepUpper)
        return false;
    }

    if (c.cascCosPA() < std::cos(cCosPACascCutPtDepP0 - cCosPACascCutPtDepP1 * c.pt()))
      return false;
    if (c.cascTransRadius() > cMaxCascradiuscut || c.cascTransRadius() < cMinCascradiuscut)
      return false;
    if (std::abs(c.mXi() - cMassXiminus) > cMasswindowCasccut)
      return false;

    return true;
  }

  template <bool IsMix, typename CollisionT, typename CascadesT, typename V0sT>
  void fill(const CollisionT& collision, const CascadesT& cascades, const V0sT& v0s)
  {
    auto cent = collision.cent();

    // Fill event QA histograms (only for same-event)
    if constexpr (!IsMix) {
      histos.fill(HIST("Event/posZ"), collision.posZ());
      histos.fill(HIST("Event/centrality"), cent);
      histos.fill(HIST("Event/posZvsCent"), collision.posZ(), cent);
      histos.fill(HIST("Event/nCascades"), cascades.size());
      histos.fill(HIST("Event/nV0s"), v0s.size());
    }

    // Count candidates after cuts
    int nCascAfterCuts = 0;
    int nV0sAfterCuts = 0;

    for (const auto& xi : cascades) {
      // QA before Xi cuts - detailed histograms
      histos.fill(HIST("QAbefore/xiMass"), xi.mXi());
      histos.fill(HIST("QAbefore/xiPt"), xi.pt());
      histos.fill(HIST("QAbefore/xiEta"), xi.eta());
      histos.fill(HIST("QAbefore/xiDCAxy"), xi.pt(), xi.dcaXYCascToPV());
      histos.fill(HIST("QAbefore/xiDCAz"), xi.pt(), xi.dcaZCascToPV());
      histos.fill(HIST("QAbefore/xiV0CosPA"), xi.pt(), xi.v0CosPA());
      histos.fill(HIST("QAbefore/xiCascCosPA"), xi.pt(), xi.cascCosPA());
      histos.fill(HIST("QAbefore/xiV0Radius"), xi.pt(), xi.transRadius());
      histos.fill(HIST("QAbefore/xiCascRadius"), xi.pt(), xi.cascTransRadius());
      histos.fill(HIST("QAbefore/xiV0DauDCA"), xi.pt(), xi.daughDCA());
      histos.fill(HIST("QAbefore/xiCascDauDCA"), xi.pt(), xi.cascDaughDCA());

      if (!cascprimaryTrackCut(xi))
        continue;
      if (!casctopCut(xi))
        continue;

      // Count cascades passing cuts
      if constexpr (!IsMix) {
        nCascAfterCuts++;
      }

      // QA after Xi cuts - detailed histograms
      histos.fill(HIST("QAafter/xiMass"), xi.mXi());
      histos.fill(HIST("QAafter/xiPt"), xi.pt());
      histos.fill(HIST("QAafter/xiEta"), xi.eta());
      histos.fill(HIST("QAafter/xiDCAxy"), xi.pt(), xi.dcaXYCascToPV());
      histos.fill(HIST("QAafter/xiDCAz"), xi.pt(), xi.dcaZCascToPV());
      histos.fill(HIST("QAafter/xiV0CosPA"), xi.pt(), xi.v0CosPA());
      histos.fill(HIST("QAafter/xiCascCosPA"), xi.pt(), xi.cascCosPA());
      histos.fill(HIST("QAafter/xiV0Radius"), xi.pt(), xi.transRadius());
      histos.fill(HIST("QAafter/xiCascRadius"), xi.pt(), xi.cascTransRadius());
      histos.fill(HIST("QAafter/xiV0DauDCA"), xi.pt(), xi.daughDCA());
      histos.fill(HIST("QAafter/xiCascDauDCA"), xi.pt(), xi.cascDaughDCA());

      // Build Xi + K0s
      for (const auto& v0 : v0s) {
        // QA before K0s selection - detailed histograms
        histos.fill(HIST("QAbefore/k0sMassPt"), v0.pt(), v0.mK0Short());
        histos.fill(HIST("QAbefore/k0sPt"), v0.pt());
        histos.fill(HIST("QAbefore/k0sEta"), v0.eta());
        histos.fill(HIST("QAbefore/k0sCosPA"), v0.pt(), v0.v0CosPA());
        histos.fill(HIST("QAbefore/k0sRadius"), v0.pt(), v0.transRadius());
        histos.fill(HIST("QAbefore/k0sDauDCA"), v0.pt(), v0.daughDCA());
        histos.fill(HIST("QAbefore/k0sDCAtoPV"), v0.pt(), std::abs(v0.dcav0topv()));
        // Calculate proper lifetime manually
        float dx = v0.decayVtxX() - collision.posX();
        float dy = v0.decayVtxY() - collision.posY();
        float dz = v0.decayVtxZ() - collision.posZ();
        float l = std::sqrt(dx * dx + dy * dy + dz * dz);
        float p = std::sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
        auto properLifetime = (l / (p + 1e-10)) * MassK0Short;
        histos.fill(HIST("QAbefore/k0sProperLifetime"), v0.pt(), properLifetime);
        // Skip Armenteros plot for now - requires daughter momentum info
        // histos.fill(HIST("QAbefore/k0sArmenteros"), alpha, qt);
        histos.fill(HIST("QAbefore/k0sDauPosDCA"), v0.pt(), std::abs(v0.dcapostopv()));
        histos.fill(HIST("QAbefore/k0sDauNegDCA"), v0.pt(), std::abs(v0.dcanegtopv()));

        if (!v0CutEnhanced(collision, v0))
          continue;

        // Count V0s passing cuts
        if constexpr (!IsMix) {
          nV0sAfterCuts++;
        }

        // QA after K0s selection - detailed histograms
        histos.fill(HIST("QAafter/k0sMassPt"), v0.pt(), v0.mK0Short());
        histos.fill(HIST("QAafter/k0sPt"), v0.pt());
        histos.fill(HIST("QAafter/k0sEta"), v0.eta());
        histos.fill(HIST("QAafter/k0sCosPA"), v0.pt(), v0.v0CosPA());
        histos.fill(HIST("QAafter/k0sRadius"), v0.pt(), v0.transRadius());
        histos.fill(HIST("QAafter/k0sDauDCA"), v0.pt(), v0.daughDCA());
        histos.fill(HIST("QAafter/k0sDCAtoPV"), v0.pt(), std::abs(v0.dcav0topv()));
        histos.fill(HIST("QAafter/k0sProperLifetime"), v0.pt(), properLifetime);
        // Skip Armenteros plot for now - requires daughter momentum info
        // histos.fill(HIST("QAafter/k0sArmenteros"), alpha, qt);
        histos.fill(HIST("QAafter/k0sDauPosDCA"), v0.pt(), std::abs(v0.dcapostopv()));
        histos.fill(HIST("QAafter/k0sDauNegDCA"), v0.pt(), std::abs(v0.dcanegtopv()));

        // 4-vectors
        ROOT::Math::PxPyPzEVector pXi, pK0s, pRes;
        pXi = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(xi.pt(), xi.eta(), xi.phi(), xi.mXi()));
        pK0s = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), massK0));
        pRes = pXi + pK0s;

        if constexpr (!IsMix) {
          histos.fill(HIST("omega2012/invmass"), pRes.M());
          histos.fill(HIST("omega2012/massPtCent"), pRes.M(), pRes.Pt(), cent);
        } else {
          histos.fill(HIST("omega2012/invmass_Mix"), pRes.M());
          histos.fill(HIST("omega2012/massPtCent_Mix"), pRes.M(), pRes.Pt(), cent);
        }
      }
    }

    // Fill event QA for after-cuts counters (only for same-event)
    if constexpr (!IsMix) {
      histos.fill(HIST("Event/nCascadesAfterCuts"), nCascAfterCuts);
      histos.fill(HIST("Event/nV0sAfterCuts"), nV0sAfterCuts);
    }
  }

  void processDummy(aod::ResoCollision const& /*collision*/)
  {
    // Dummy function to satisfy the compiler
  }
  PROCESS_SWITCH(Omega2012Analysis, processDummy, "Process Dummy", true);

  void processData(const aod::ResoCollision& collision,
                   aod::ResoCascades const& resocasc,
                   aod::ResoV0s const& resov0s)
  {
    fill<false>(collision, resocasc, resov0s);
  }
  PROCESS_SWITCH(Omega2012Analysis, processData, "Process Event for data", false);

  void processMixedEvent(const aod::ResoCollisions& collisions,
                         aod::ResoCascades const& resocasc,
                         aod::ResoV0s const& resov0s)
  {

    auto cascV0sTuple = std::make_tuple(resocasc, resov0s);
    Pair<aod::ResoCollisions, aod::ResoCascades, aod::ResoV0s, BinningTypeVertexContributor> pairs{colBinning, nEvtMixing, -1, collisions, cascV0sTuple, &cache};

    for (auto& [collision1, casc1, collision2, v0s2] : pairs) { // o2-linter: disable=const-ref-in-for-loop (structured binding cannot be const in this context)
      auto cent = collision1.cent();

      for (const auto& xi : casc1) {
        if (!cascprimaryTrackCut(xi))
          continue;
        if (!casctopCut(xi))
          continue;

        for (const auto& v0 : v0s2) {
          if (!v0CutEnhanced(collision2, v0))
            continue;

          // 4-vectors for mixed event
          ROOT::Math::PxPyPzEVector pXi, pK0s, pRes;
          pXi = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(xi.pt(), xi.eta(), xi.phi(), xi.mXi()));
          pK0s = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), massK0));
          pRes = pXi + pK0s;

          histos.fill(HIST("omega2012/invmass_Mix"), pRes.M());
          histos.fill(HIST("omega2012/massPtCent_Mix"), pRes.M(), pRes.Pt(), cent);
        }
      }
    }
  }
  PROCESS_SWITCH(Omega2012Analysis, processMixedEvent, "Process Mixed Event", false);

  // MC processes - placeholder for future implementation
  void processMC(const aod::ResoCollision& /*collision*/,
                 aod::ResoCascades const& /*resocasc*/,
                 aod::ResoV0s const& /*resov0s*/,
                 aod::McParticles const& /*mcParticles*/)
  {
    // TODO: Implement MC truth matching for Xi + K0s
    // - Match reconstructed Xi to MC Xi
    // - Match reconstructed K0s to MC K0s
    // - Fill MC truth histograms
    // - Fill reconstruction efficiency histograms
    // - Check if the Xi and K0s come from same Omega(2012) mother
  }
  PROCESS_SWITCH(Omega2012Analysis, processMC, "Process MC with truth matching (placeholder)", false);

  void processMCGenerated(aod::McParticles const& mcParticles)
  {
    // Process MC generated particles (no reconstruction requirement)
    // Omega(2012) PDG code: Need to check - likely custom or using generator-specific codes
    // This resonance decays to Xi + K0s

    for (const auto& mcParticle : mcParticles) {
      // Look for Omega(2012) - PDG code may vary by generator
      int pdg = mcParticle.pdgCode();

      // TODO: Determine correct PDG code for Omega(2012)
      // Placeholder check - update with correct PDG code
      if (std::abs(pdg) != kPlaceholderPdgCode) // o2-linter: disable=pdg/explicit-code (placeholder for generator-specific PDG code)
        continue;

      // Fill generated level histograms
      auto pt = mcParticle.pt();
      auto eta = mcParticle.eta();
      auto y = mcParticle.y();

      histos.fill(HIST("MC/hMCGenOmega2012Pt"), pt);
      histos.fill(HIST("MC/hMCGenOmega2012PtEta"), pt, eta);
      histos.fill(HIST("MC/hMCGenOmega2012Y"), y);

      // Get daughters
      auto daughters = mcParticle.daughters_as<aod::McParticles>();
      if (daughters.size() != kNumExpectedDaughters)
        continue;

      int daughter1PDG = 0, daughter2PDG = 0;
      ROOT::Math::PxPyPzEVector p1, p2, pMother;

      int iDaughter = 0;
      for (const auto& daughter : daughters) {
        if (iDaughter == 0) {
          daughter1PDG = daughter.pdgCode();
          p1 = ROOT::Math::PxPyPzEVector(daughter.px(), daughter.py(), daughter.pz(), daughter.e());
        } else {
          daughter2PDG = daughter.pdgCode();
          p2 = ROOT::Math::PxPyPzEVector(daughter.px(), daughter.py(), daughter.pz(), daughter.e());
        }
        iDaughter++;
      }

      pMother = p1 + p2;

      // Check decay channels
      auto motherPt = pMother.Pt();
      auto motherM = pMother.M();

      // Xi- + K0s or Xi+ + K0s
      if ((std::abs(daughter1PDG) == kXiMinus && daughter2PDG == kK0Short) ||
          (std::abs(daughter2PDG) == kXiMinus && daughter1PDG == kK0Short)) {
        histos.fill(HIST("MC/hMCTruthInvMassXiK0s"), motherM);
        histos.fill(HIST("MC/hMCTruthMassPtXiK0s"), motherM, motherPt);
      }
    }
  }
  PROCESS_SWITCH(Omega2012Analysis, processMCGenerated, "Process MC generated particles (placeholder)", false);

  // Fill function for 3-body decay analysis
  template <bool IsResoMicrotrack, typename CollisionT, typename CascadesT, typename V0sT, typename TracksT>
  void fillThreeBody(const CollisionT& collision, const CascadesT& cascades, const V0sT& v0s, const TracksT& tracks)
  {
    auto cent = collision.cent();

    // Collect track IDs used in xi and v0 construction to exclude them from pion selection
    std::set<int> usedTrackIds;

    // Collect track IDs from xi cascades
    for (const auto& xi : cascades) {
      if (!cascprimaryTrackCut(xi))
        continue;
      if (!casctopCut(xi))
        continue;

      // Add xi daughter track IDs from cascadeIndices array (ordered: positive, negative, bachelor)
      auto cascIndices = xi.cascadeIndices();
      usedTrackIds.insert(cascIndices[0]); // positive track
      usedTrackIds.insert(cascIndices[1]); // negative track
      usedTrackIds.insert(cascIndices[2]); // bachelor track
    }

    // Collect track IDs from v0s
    for (const auto& v0 : v0s) {
      if (!v0CutEnhanced(collision, v0))
        continue;

      // Add v0 daughter track IDs from indices array
      auto v0Indices = v0.indices();
      usedTrackIds.insert(v0Indices[0]); // positive track
      usedTrackIds.insert(v0Indices[1]); // negative track
    }

    // First loop: xi + pion to check xi1530 mass window
    for (const auto& xi : cascades) {
      if (!cascprimaryTrackCut(xi))
        continue;
      if (!casctopCut(xi))
        continue;

      for (const auto& pion : tracks) {
        // Skip pion tracks that are already used in xi construction
        if (usedTrackIds.find(pion.globalIndex()) != usedTrackIds.end()) {
          continue;
        }

        // Pion QA before cuts
        histos.fill(HIST("QAbefore/pionPt"), pion.pt());
        histos.fill(HIST("QAbefore/pionEta"), pion.eta());

        if constexpr (IsResoMicrotrack) {
          histos.fill(HIST("QAbefore/pionDCAxy"), pion.pt(), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(pion.trackSelectionFlags()));
          histos.fill(HIST("QAbefore/pionDCAz"), pion.pt(), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(pion.trackSelectionFlags()));
          histos.fill(HIST("QAbefore/pionTPCNSigma"), pion.pt(), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(pion.pidNSigmaPiFlag()));
          if (pion.hasTOF()) {
            histos.fill(HIST("QAbefore/pionTOFNSigma"), pion.pt(), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(pion.pidNSigmaPiFlag()));
          }
        } else {
          histos.fill(HIST("QAbefore/pionDCAxy"), pion.pt(), pion.dcaXY());
          histos.fill(HIST("QAbefore/pionDCAz"), pion.pt(), pion.dcaZ());
          histos.fill(HIST("QAbefore/pionTPCNSigma"), pion.pt(), pion.tpcNSigmaPi());
          if (pion.hasTOF()) {
            histos.fill(HIST("QAbefore/pionTOFNSigma"), pion.pt(), pion.tofNSigmaPi());
          }
          if constexpr (requires { pion.tpcNClsFound(); }) {
            histos.fill(HIST("QAbefore/pionTPCNcls"), pion.tpcNClsFound());
          }
        }

        if (!pionCut<IsResoMicrotrack>(pion))
          continue;

        // Pion QA after cuts
        histos.fill(HIST("QAafter/pionPt"), pion.pt());
        histos.fill(HIST("QAafter/pionEta"), pion.eta());

        if constexpr (IsResoMicrotrack) {
          histos.fill(HIST("QAafter/pionDCAxy"), pion.pt(), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(pion.trackSelectionFlags()));
          histos.fill(HIST("QAafter/pionDCAz"), pion.pt(), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(pion.trackSelectionFlags()));
          histos.fill(HIST("QAafter/pionTPCNSigma"), pion.pt(), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(pion.pidNSigmaPiFlag()));
          if (pion.hasTOF()) {
            histos.fill(HIST("QAafter/pionTOFNSigma"), pion.pt(), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(pion.pidNSigmaPiFlag()));
          }
        } else {
          histos.fill(HIST("QAafter/pionDCAxy"), pion.pt(), pion.dcaXY());
          histos.fill(HIST("QAafter/pionDCAz"), pion.pt(), pion.dcaZ());
          histos.fill(HIST("QAafter/pionTPCNSigma"), pion.pt(), pion.tpcNSigmaPi());
          if (pion.hasTOF()) {
            histos.fill(HIST("QAafter/pionTOFNSigma"), pion.pt(), pion.tofNSigmaPi());
          }
          if constexpr (requires { pion.tpcNClsFound(); }) {
            histos.fill(HIST("QAafter/pionTPCNcls"), pion.tpcNClsFound());
          }
        }

        // Check xi1530 mass window cut
        if (!xi1530MassCut(xi, pion))
          continue;

        // Second loop: v0 for the selected xi-pion pair
        for (const auto& v0 : v0s) {
          if (!v0CutEnhanced(collision, v0))
            continue;

          // 4-vectors for 3-body decay: Xi + K0s + pion
          ROOT::Math::PxPyPzEVector pXi, pK0s, pPion, pRes;
          pXi = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(xi.pt(), xi.eta(), xi.phi(), xi.mXi()));
          pK0s = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(v0.pt(), v0.eta(), v0.phi(), massK0));
          pPion = ROOT::Math::PxPyPzEVector(ROOT::Math::PtEtaPhiMVector(pion.pt(), pion.eta(), pion.phi(), MassPionCharged));

          pRes = pXi + pK0s + pPion;

          histos.fill(HIST("omega2012_3body/invmass"), pRes.M());
          histos.fill(HIST("omega2012_3body/massPtCent"), pRes.M(), pRes.Pt(), cent);
        }
      }
    }
  }

  // 3-body decay analysis: Xi + pi + K0s with ResoTracks
  void processThreeBodyWithTracks(const aod::ResoCollision& collision,
                                  aod::ResoCascades const& resocasc,
                                  aod::ResoV0s const& resov0s,
                                  aod::ResoTracks const& resotracks)
  {
    fillThreeBody<false>(collision, resocasc, resov0s, resotracks);
  }
  PROCESS_SWITCH(Omega2012Analysis, processThreeBodyWithTracks, "Process 3-body decay with ResoTracks", false);

  // 3-body decay analysis: Xi + pi + K0s with ResoMicroTracks
  void processThreeBodyWithMicroTracks(const aod::ResoCollision& collision,
                                       aod::ResoCascades const& resocasc,
                                       aod::ResoV0s const& resov0s,
                                       aod::ResoMicroTracks const& resomicrotracks)
  {
    fillThreeBody<true>(collision, resocasc, resov0s, resomicrotracks);
  }
  PROCESS_SWITCH(Omega2012Analysis, processThreeBodyWithMicroTracks, "Process 3-body decay with ResoMicroTracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Omega2012Analysis>(cfgc)};
}
