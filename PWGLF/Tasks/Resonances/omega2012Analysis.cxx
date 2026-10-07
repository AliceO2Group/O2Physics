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
///
/// Two decay modes, analysed separately (never merged in any histogram or process function):
///  - mode A (XiK0s): Omega(2012)- -> Xi- K0S. processData, processMixedEvent, processMC, processMCGenerated.
///    Histograms: Event/, QAbefore/, QAafter/, omega2012/, MC/ (unchanged).
///  - mode B (Xi1530K): Omega(2012)- -> Xi(1530)0 K- -> Xi- pi+ K- with a charged kaon track.
///    processXi1530KMicro, processXi1530KTracks, processXi1530KMCMicro, processXi1530KMixedMicro.
///    Histograms: xi1530K/ (signal charge pattern) and xi1530K_wrongSign/ (charge-pattern controls).
/// Selection, candidate enumeration and truth classification: PWGLF/Core/Omega2012AnalysisCore.h.
///
/// Configuration keys: all mode-A keys are unchanged. The former three-body (Xi pi K0S) process functions
/// processThreeBodyWithTracks / processThreeBodyWithMicroTracks are removed (that final state cannot be an Omega(2012)-).
/// Mode-B pion/kaon track selection: the common resonance TrackCuts with the prefix "trk.":
///   cPionPtMin -> trk.cMinPtcut, cPionEtaMax -> trk.cMaxEtacut, cPionDCAxyMax -> trk.cMaxDCArToPVcut,
///   cPionDCAzMax -> trk.cMaxDCAzToPVcut (default 0.15 cm), cPionTPCNClusMin -> trk.cfgTPCcluster;
///   further trk.* keys as in PWGLF/Core/ResoAnalysisSelectionCore.h.
/// Pion PID keys cPion* are unchanged; the kaon PID keys cKaon* have the same shape.
/// New mode-B keys: cXi1530UseMassWindow, cXi1530KFillWrongSign, cByPassTOF.

#include "PWGLF/Core/Omega2012AnalysisCore.h"
#include "PWGLF/Core/Omega2012MlFeatures.h"
#include "PWGLF/Core/ResoAnalysisSelectionCore.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>

#include <cmath>
#include <cstdint>
#include <tuple>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::analysis::omega2012;

struct Omega2012Analysis {
  // Constants
  static constexpr int NumExpectedDaughters = 2; // Expected number of daughters for 2-body decay
  SliceCache cache;
  Preslice<aod::ResoCascades> perResoCollisionCasc = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoV0s> perResoCollisionV0 = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoTracks> perResoCollisionTrack = aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoMicroTracks_001> perResoCollisionMicroTrack = aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Axes
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0}, "pT"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pT (QA)"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Centrality"};

  // Invariant mass range for Omega(2012)
  Configurable<float> cInvMassStart{"cInvMassStart", 1.6, "Invariant mass start (GeV/c^2)"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 2.2, "Invariant mass end (GeV/c^2)"};
  Configurable<int> cInvMassBins{"cInvMassBins", 600, "Invariant mass bins"};

  // Selection shared with the Omega(2012) training-table task (plain JSON keys; the mode-B track cuts carry "trk.")
  o2::analysis::resonance::EventCuts eventCuts;
  XiCuts xiCuts;
  K0sCuts k0sCuts;
  Xi1530KCuts xi1530KCuts;
  o2::analysis::resonance::TrackCuts trackCuts = makeXi1530KTrackCuts();
  PionPidCuts pionPID;
  KaonPidCuts kaonPID;
  CandidateCuts candidateCuts;

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - centrality"};

  // Module-initializer collision and daughter tables (matches the resonance-module-initializer producer)
  using ResoCollisions = aod::ResoCollisions_001;
  using ResoMCCollisions = soa::Join<ResoCollisions, aod::ResoMCCollisions_001>;
  using ResoMicroTracks = aod::ResoMicroTracks_001;
  using ResoMCCascades = soa::Join<aod::ResoCascades, aod::ResoMCCascades>;
  using ResoMCV0s = soa::Join<aod::ResoV0s, aod::ResoMCV0s>;
  using ResoMCMicroTracks = soa::Join<ResoMicroTracks, aod::ResoMCMicroTracks_001>;

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  BinningTypeVertexContributor colBinning{{cfgVtxBins, cfgMultBins}};

  Omega2012AnalysisCore core;

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
    AxisSpec nsigmaAxis = {100, -5.0, 5.0, "N#sigma"};
    AxisSpec crossedRowsAxis = {160, 0, 160, "TPC crossed rows"};
    AxisSpec openingAngleAxis = {200, -0.15, 6.45, "#alpha_{oa}"};
    AxisSpec omegaKinPtAxis = {500, 0, 10, "#it{p}_{T} (GeV/#it{c})"};

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
    histos.add("QAbefore/k0sDauTPCNsigmaPosPi", "K0s positive daughter pion TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    histos.add("QAbefore/k0sDauTPCNsigmaNegPi", "K0s negative daughter pion TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    histos.add("QAbefore/k0sNCrossedRowsPos", "K0s positive daughter crossed rows before cuts", kTH2F, {ptAxisQA, crossedRowsAxis});
    histos.add("QAbefore/k0sNCrossedRowsNeg", "K0s negative daughter crossed rows before cuts", kTH2F, {ptAxisQA, crossedRowsAxis});

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
    histos.add("QAafter/k0sDauTPCNsigmaPosPi", "K0s positive daughter pion TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    histos.add("QAafter/k0sDauTPCNsigmaNegPi", "K0s negative daughter pion TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
    histos.add("QAafter/k0sNCrossedRowsPos", "K0s positive daughter crossed rows after cuts", kTH2F, {ptAxisQA, crossedRowsAxis});
    histos.add("QAafter/k0sNCrossedRowsNeg", "K0s negative daughter crossed rows after cuts", kTH2F, {ptAxisQA, crossedRowsAxis});

    // Resonance (2-body decay: Xi + K0s)
    histos.add("omega2012/invmass", "Invariant mass of Omega(2012) → Xi + K0s", kTH1F, {invMassAxis});
    histos.add("omega2012/invmass_Mix", "Mixed event Invariant mass of Omega(2012) → Xi + K0s", kTH1F, {invMassAxis});
    histos.add("omega2012/massPtCent", "Omega(2012) mass vs pT vs cent", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("omega2012/massPtCent_Mix", "Mixed event Omega(2012) mass vs pT vs cent", kTH3F, {invMassAxis, ptAxis, centAxis});
    histos.add("QAbefore/omegaAlphaVsPt", "#alpha_{oa} vs p_{T} before kinematic cuts", kTH2F, {omegaKinPtAxis, openingAngleAxis});
    histos.add("QAafter/omegaAlphaVsPt", "#alpha_{oa} vs p_{T} after kinematic cuts", kTH2F, {omegaKinPtAxis, openingAngleAxis});

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

    ProcessModes modes;
    modes.xiK0s = doprocessData || doprocessMixedEvent || doprocessMC;
    modes.xi1530K = doprocessXi1530KMicro || doprocessXi1530KTracks || doprocessXi1530KMCMicro || doprocessXi1530KMixedMicro;
    modes.microTracks = doprocessXi1530KMicro || doprocessXi1530KMCMicro || doprocessXi1530KMixedMicro;
    modes.mcReco = doprocessMC || doprocessXi1530KMCMicro;
    modes.mixing = doprocessMixedEvent || doprocessXi1530KMixedMicro;
    core.init(histos, eventCuts, xiCuts, k0sCuts, xi1530KCuts, trackCuts, pionPID, kaonPID, candidateCuts, modes);

    if (modes.xi1530K) {
      AxisSpec xiPiMassAxis = {300, 1.4, 1.7, "M_{#Xi#pi} (GeV/#it{c}^{2})"};
      AxisSpec patternAxis = {3, 0.5, 3.5, "charge pattern (1: wrong-sign #pi, 2: wrong-sign K, 3: both)"};
      AxisSpec etaAxisQA = {100, -2.0, 2.0, "#eta"};

      // Signal charge pattern: Xi- pi+ K- and Xi+ pi- K+
      histos.add("xi1530K/invmass", "Invariant mass of Omega(2012) → #Xi(1530)^{0} K → #Xi #pi K", kTH1F, {invMassAxis});
      histos.add("xi1530K/massPtCent", "Omega(2012) → #Xi(1530)^{0} K mass vs pT vs cent", kTH3F, {invMassAxis, ptAxis, centAxis});
      histos.add("xi1530K/invmass_Mix", "Mixed event invariant mass of Omega(2012) → #Xi(1530)^{0} K", kTH1F, {invMassAxis});
      histos.add("xi1530K/massPtCent_Mix", "Mixed event Omega(2012) → #Xi(1530)^{0} K mass vs pT vs cent", kTH3F, {invMassAxis, ptAxis, centAxis});
      histos.add("xi1530K/massXiPi", "#Xi #pi mass of selected candidates (before the rapidity cut)", kTH1F, {xiPiMassAxis});
      histos.add("xi1530K/massXiPiVsMass", "#Xi #pi mass vs #Xi #pi K mass", kTH2F, {invMassAxis, xiPiMassAxis});

      // QA of the mode-B inputs (once per object and collision)
      histos.add("xi1530K/QAbefore/xiMass", "Xi mass before cuts", kTH1F, {xiMassAxis});
      histos.add("xi1530K/QAbefore/xiPt", "Xi pT before cuts", kTH1F, {ptAxisQA});
      histos.add("xi1530K/QAbefore/xiEta", "Xi eta before cuts", kTH1F, {etaAxisQA});
      histos.add("xi1530K/QAafter/xiMass", "Xi mass after cuts", kTH1F, {xiMassAxis});
      histos.add("xi1530K/QAafter/xiPt", "Xi pT after cuts", kTH1F, {ptAxisQA});
      histos.add("xi1530K/QAafter/xiEta", "Xi eta after cuts", kTH1F, {etaAxisQA});
      histos.add("xi1530K/QAbefore/trackPt", "Track pT before cuts", kTH1F, {ptAxisQA});
      histos.add("xi1530K/QAbefore/trackEta", "Track eta before cuts", kTH1F, {etaAxisQA});
      histos.add("xi1530K/QAbefore/trackDCAxy", "Track DCAxy before cuts", kTH2F, {ptAxisQA, dcaxyAxis});
      histos.add("xi1530K/QAbefore/trackDCAz", "Track DCAz before cuts", kTH2F, {ptAxisQA, dcazAxis});
      histos.add("xi1530K/QAbefore/pionTPCNSigma", "Pion TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("xi1530K/QAbefore/pionTOFNSigma", "Pion TOF NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("xi1530K/QAbefore/kaonTPCNSigma", "Kaon TPC NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("xi1530K/QAbefore/kaonTOFNSigma", "Kaon TOF NSigma before cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("xi1530K/QAafter/pionPt", "Pion pT after cuts", kTH1F, {ptAxisQA});
      histos.add("xi1530K/QAafter/pionEta", "Pion eta after cuts", kTH1F, {etaAxisQA});
      histos.add("xi1530K/QAafter/pionDCAxy", "Pion DCAxy after cuts", kTH2F, {ptAxisQA, dcaxyAxis});
      histos.add("xi1530K/QAafter/pionDCAz", "Pion DCAz after cuts", kTH2F, {ptAxisQA, dcazAxis});
      histos.add("xi1530K/QAafter/pionTPCNSigma", "Pion TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("xi1530K/QAafter/pionTOFNSigma", "Pion TOF NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("xi1530K/QAafter/kaonPt", "Kaon pT after cuts", kTH1F, {ptAxisQA});
      histos.add("xi1530K/QAafter/kaonEta", "Kaon eta after cuts", kTH1F, {etaAxisQA});
      histos.add("xi1530K/QAafter/kaonDCAxy", "Kaon DCAxy after cuts", kTH2F, {ptAxisQA, dcaxyAxis});
      histos.add("xi1530K/QAafter/kaonDCAz", "Kaon DCAz after cuts", kTH2F, {ptAxisQA, dcazAxis});
      histos.add("xi1530K/QAafter/kaonTPCNSigma", "Kaon TPC NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});
      histos.add("xi1530K/QAafter/kaonTOFNSigma", "Kaon TOF NSigma after cuts", kTH2F, {ptAxisQA, nsigmaAxis});

      // Charge-pattern controls (never part of the signal)
      histos.add("xi1530K_wrongSign/invmassPattern", "Wrong-sign #Xi #pi K mass by charge pattern", kTH2F, {patternAxis, invMassAxis});
      histos.add("xi1530K_wrongSign/massPtPattern", "Wrong-sign #Xi #pi K mass vs pT by charge pattern", kTH3F, {invMassAxis, ptAxis, patternAxis});

      if (doprocessXi1530KMCMicro) {
        histos.add("xi1530K/MC/hMCRecOmega2012Pt", "MC reconstructed Omega(2012) → #Xi(1530)^{0} K pT", kTH1F, {ptAxis});
        histos.add("xi1530K/MC/hMCRecOmega2012PtEta", "MC reconstructed Omega(2012) → #Xi(1530)^{0} K pT vs eta", kTH2F, {ptAxis, etaAxis});
        histos.add("xi1530K/MC/hMCRecMass", "MC reconstructed Omega(2012) → #Xi(1530)^{0} K mass", kTH1F, {invMassAxis});
        histos.add("xi1530K/MC/hMCRecMassXiPi", "MC reconstructed #Xi #pi mass", kTH1F, {xiPiMassAxis});
        histos.add("xi1530K/MC/hMCRecXiPt", "MC reconstructed Xi pT", kTH1F, {ptAxis});
        histos.add("xi1530K/MC/hMCRecPionPt", "MC reconstructed pion pT", kTH1F, {ptAxis});
        histos.add("xi1530K/MC/hMCRecKaonPt", "MC reconstructed kaon pT", kTH1F, {ptAxis});
      }
    }

    LOG(info) << "Size of the histograms in Omega(2012) analysis task";
    histos.print();
  }

  // Mode A, same event (data): event QA, Xi and K0S QA, kinematic-cut QA and the Xi K0S mass
  template <typename CollisionT, typename CascadesT, typename V0sT>
  void fillXiK0s(const CollisionT& collision, const CascadesT& cascades, const V0sT& v0s)
  {
    auto cent = collision.cent();

    histos.fill(HIST("Event/posZ"), collision.posZ());
    histos.fill(HIST("Event/centrality"), cent);
    histos.fill(HIST("Event/posZvsCent"), collision.posZ(), cent);
    histos.fill(HIST("Event/nCascades"), cascades.size());
    histos.fill(HIST("Event/nV0s"), v0s.size());

    // Count candidates after cuts
    int nCascAfterCuts = 0;
    int nV0sAfterCuts = 0;

    auto onK0s = [&](auto const& v0, double properLifetime, bool selected) {
      histos.fill(HIST("QAbefore/k0sMassPt"), v0.pt(), v0.mK0Short());
      histos.fill(HIST("QAbefore/k0sPt"), v0.pt());
      histos.fill(HIST("QAbefore/k0sEta"), v0.eta());
      histos.fill(HIST("QAbefore/k0sCosPA"), v0.pt(), v0.v0CosPA());
      histos.fill(HIST("QAbefore/k0sRadius"), v0.pt(), v0.transRadius());
      histos.fill(HIST("QAbefore/k0sDauDCA"), v0.pt(), v0.daughDCA());
      histos.fill(HIST("QAbefore/k0sDCAtoPV"), v0.pt(), std::abs(v0.dcav0topv()));
      histos.fill(HIST("QAbefore/k0sProperLifetime"), v0.pt(), properLifetime);
      histos.fill(HIST("QAbefore/k0sArmenteros"), v0.alpha(), v0.qtarm());
      histos.fill(HIST("QAbefore/k0sDauPosDCA"), v0.pt(), std::abs(v0.dcapostopv()));
      histos.fill(HIST("QAbefore/k0sDauNegDCA"), v0.pt(), std::abs(v0.dcanegtopv()));
      histos.fill(HIST("QAbefore/k0sDauTPCNsigmaPosPi"), v0.pt(), v0.daughterTPCNSigmaPosPi());
      histos.fill(HIST("QAbefore/k0sDauTPCNsigmaNegPi"), v0.pt(), v0.daughterTPCNSigmaNegPi());
      histos.fill(HIST("QAbefore/k0sNCrossedRowsPos"), v0.pt(), v0.nCrossedRowsPos());
      histos.fill(HIST("QAbefore/k0sNCrossedRowsNeg"), v0.pt(), v0.nCrossedRowsNeg());
      if (!selected) {
        return;
      }
      nV0sAfterCuts++;
      histos.fill(HIST("QAafter/k0sMassPt"), v0.pt(), v0.mK0Short());
      histos.fill(HIST("QAafter/k0sPt"), v0.pt());
      histos.fill(HIST("QAafter/k0sEta"), v0.eta());
      histos.fill(HIST("QAafter/k0sCosPA"), v0.pt(), v0.v0CosPA());
      histos.fill(HIST("QAafter/k0sRadius"), v0.pt(), v0.transRadius());
      histos.fill(HIST("QAafter/k0sDauDCA"), v0.pt(), v0.daughDCA());
      histos.fill(HIST("QAafter/k0sDCAtoPV"), v0.pt(), std::abs(v0.dcav0topv()));
      histos.fill(HIST("QAafter/k0sProperLifetime"), v0.pt(), properLifetime);
      histos.fill(HIST("QAafter/k0sArmenteros"), v0.alpha(), v0.qtarm());
      histos.fill(HIST("QAafter/k0sDauPosDCA"), v0.pt(), std::abs(v0.dcapostopv()));
      histos.fill(HIST("QAafter/k0sDauNegDCA"), v0.pt(), std::abs(v0.dcanegtopv()));
      histos.fill(HIST("QAafter/k0sDauTPCNsigmaPosPi"), v0.pt(), v0.daughterTPCNSigmaPosPi());
      histos.fill(HIST("QAafter/k0sDauTPCNsigmaNegPi"), v0.pt(), v0.daughterTPCNSigmaNegPi());
      histos.fill(HIST("QAafter/k0sNCrossedRowsPos"), v0.pt(), v0.nCrossedRowsPos());
      histos.fill(HIST("QAafter/k0sNCrossedRowsNeg"), v0.pt(), v0.nCrossedRowsNeg());
    };

    auto onXi = [&](auto const& xi, bool selected) {
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
      if (!selected) {
        return;
      }
      nCascAfterCuts++;
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
    };

    const bool kinCutsOn = core.kinCutsEnabled();
    auto onCandidate = [&](auto const& /*xi*/, auto const& /*v0*/, XiK0sCandidateValues const& c) {
      if (kinCutsOn) {
        histos.fill(HIST("QAbefore/omegaAlphaVsPt"), c.omega.Pt(), c.alpha);
        if (!c.passesKinCut) {
          return;
        }
        histos.fill(HIST("QAafter/omegaAlphaVsPt"), c.omega.Pt(), c.alpha);
      }
      if (!c.inRapidity) {
        return;
      }
      histos.fill(HIST("omega2012/invmass"), c.omega.M());
      histos.fill(HIST("omega2012/massPtCent"), c.omega.M(), c.omega.Pt(), cent);
    };

    core.forEachXiK0sCandidate<false, false>(histos, collision, collision, cascades, v0s, false, onXi, onK0s, onCandidate);

    histos.fill(HIST("Event/nCascadesAfterCuts"), nCascAfterCuts);
    histos.fill(HIST("Event/nV0sAfterCuts"), nV0sAfterCuts);
  }

  void processDummy(aod::ResoCollision const& /*collision*/)
  {
    // Dummy function to satisfy the compiler
  }
  PROCESS_SWITCH(Omega2012Analysis, processDummy, "Process Dummy", true);

  void processData(ResoCollisions::iterator const& collision,
                   aod::ResoCascades const& resocasc,
                   aod::ResoV0s const& resov0s)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    fillXiK0s(collision, resocasc, resov0s);
  }
  PROCESS_SWITCH(Omega2012Analysis, processData, "Process Event for data", false);

  void processMixedEvent(ResoCollisions const& collisions,
                         aod::ResoCascades const& resocasc,
                         aod::ResoV0s const& resov0s)
  {
    auto cascV0sTuple = std::make_tuple(resocasc, resov0s);
    Pair<ResoCollisions, aod::ResoCascades, aod::ResoV0s, BinningTypeVertexContributor> pairs{colBinning, nEvtMixing, -1, collisions, cascV0sTuple, &cache};

    const bool kinCutsOn = core.kinCutsEnabled();
    for (const auto& [collision1, casc1, collision2, v0s2] : pairs) {
      if (!core.passesEventCuts(collision1) || !core.passesEventCuts(collision2)) {
        continue;
      }
      auto cent = collision1.cent();
      // Xi from collision 1, K0s from collision 2 (selected with the vertex of collision 2)
      auto onCandidate = [&](auto const& /*xi*/, auto const& /*v0*/, XiK0sCandidateValues const& c) {
        if (kinCutsOn && !c.passesKinCut) {
          return;
        }
        if (!c.inRapidity) {
          return;
        }
        histos.fill(HIST("omega2012/invmass_Mix"), c.omega.M());
        histos.fill(HIST("omega2012/massPtCent_Mix"), c.omega.M(), c.omega.Pt(), cent);
      };
      core.forEachXiK0sCandidate<false, true>(histos, collision1, collision2, casc1, v0s2, false, nullptr, nullptr, onCandidate);
    }
  }
  PROCESS_SWITCH(Omega2012Analysis, processMixedEvent, "Process Mixed Event", false);

  // MC reconstructed processing: match reconstructed Xi + K0s pairs to a common Omega(2012) mother
  void processMC(ResoMCCollisions::iterator const& collision,
                 ResoMCCascades const& resocasc,
                 ResoMCV0s const& resov0s)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    // No kinematic cut in the truth-matched spectra (as before the refactoring)
    auto onCandidate = [&](auto const& xi, auto const& v0, XiK0sCandidateValues const& c) {
      if (classifyXiK0sTruth(xi, v0) != XiK0sTruth::Matched) {
        return;
      }
      if (!c.inRapidity) {
        return;
      }
      histos.fill(HIST("MC/hMCRecOmega2012Pt"), c.omega.Pt());
      histos.fill(HIST("MC/hMCRecOmega2012PtEta"), c.omega.Pt(), c.omega.Eta());
      histos.fill(HIST("MC/hMCRecXiPt"), xi.pt());
      histos.fill(HIST("MC/hMCRecK0sPt"), v0.pt());
    };
    core.forEachXiK0sCandidate<true, false>(histos, collision, collision, resocasc, resov0s, false, nullptr, nullptr, onCandidate);
  }
  PROCESS_SWITCH(Omega2012Analysis, processMC, "Process MC with truth matching", false);

  void processMCGenerated(aod::McParticles const& mcParticles)
  {
    // Process MC generated particles (no reconstruction requirement)
    // This resonance decays to Xi + K0s

    for (const auto& mcParticle : mcParticles) {
      // Look for Omega(2012)
      int pdg = mcParticle.pdgCode();

      if (std::abs(pdg) != kOmega2012Minus)
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
      if (daughters.size() != NumExpectedDaughters)
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

        const bool isDaughter1Xi = std::abs(daughter1PDG) == kXiMinus;
        const auto& pXiTruth = isDaughter1Xi ? p1 : p2;
        const auto& pK0sTruth = isDaughter1Xi ? p2 : p1;
        histos.fill(HIST("MC/hMCTrueXiPt"), pXiTruth.Pt());
        histos.fill(HIST("MC/hMCTrueK0sPt"), pK0sTruth.Pt());
      }
    }
  }
  PROCESS_SWITCH(Omega2012Analysis, processMCGenerated, "Process MC generated particles", false);

  // Mode B: Xi(1530)0 K- -> Xi- pi+ K- (and charge conjugate) from one collision, or Xi and tracks from two mixed collisions
  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionT, typename CascadesT, typename TracksT, typename TrackIdsT>
  void fillXi1530K(const CollisionT& collision, float cent, const CascadesT& cascades, const TracksT& tracks, const TrackIdsT& trackIds)
  {
    // Input QA: same event only (each object once per collision)
    auto onXi = [&](auto const& xi, bool selected) {
      if constexpr (IsMix) {
        return;
      }
      histos.fill(HIST("xi1530K/QAbefore/xiMass"), xi.mXi());
      histos.fill(HIST("xi1530K/QAbefore/xiPt"), xi.pt());
      histos.fill(HIST("xi1530K/QAbefore/xiEta"), xi.eta());
      if (!selected) {
        return;
      }
      histos.fill(HIST("xi1530K/QAafter/xiMass"), xi.mXi());
      histos.fill(HIST("xi1530K/QAafter/xiPt"), xi.pt());
      histos.fill(HIST("xi1530K/QAafter/xiEta"), xi.eta());
    };

    auto onTrack = [&](auto const& track, int pionStage, int kaonStage) {
      if constexpr (IsMix) {
        return;
      }
      const bool hasTOF = track.hasTOF();
      histos.fill(HIST("xi1530K/QAbefore/trackPt"), track.pt());
      histos.fill(HIST("xi1530K/QAbefore/trackEta"), track.eta());
      histos.fill(HIST("xi1530K/QAbefore/trackDCAxy"), track.pt(), track.dcaXY());
      histos.fill(HIST("xi1530K/QAbefore/trackDCAz"), track.pt(), track.dcaZ());
      histos.fill(HIST("xi1530K/QAbefore/pionTPCNSigma"), track.pt(), track.tpcNSigmaPi());
      histos.fill(HIST("xi1530K/QAbefore/kaonTPCNSigma"), track.pt(), track.tpcNSigmaKa());
      if (hasTOF) {
        histos.fill(HIST("xi1530K/QAbefore/pionTOFNSigma"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("xi1530K/QAbefore/kaonTOFNSigma"), track.pt(), track.tofNSigmaKa());
      }
      if (pionStage == o2::analysis::resonance::kTrkPID) {
        histos.fill(HIST("xi1530K/QAafter/pionPt"), track.pt());
        histos.fill(HIST("xi1530K/QAafter/pionEta"), track.eta());
        histos.fill(HIST("xi1530K/QAafter/pionDCAxy"), track.pt(), track.dcaXY());
        histos.fill(HIST("xi1530K/QAafter/pionDCAz"), track.pt(), track.dcaZ());
        histos.fill(HIST("xi1530K/QAafter/pionTPCNSigma"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("xi1530K/QAafter/pionTOFNSigma"), track.pt(), track.tofNSigmaPi());
        }
      }
      if (kaonStage == o2::analysis::resonance::kTrkPID) {
        histos.fill(HIST("xi1530K/QAafter/kaonPt"), track.pt());
        histos.fill(HIST("xi1530K/QAafter/kaonEta"), track.eta());
        histos.fill(HIST("xi1530K/QAafter/kaonDCAxy"), track.pt(), track.dcaXY());
        histos.fill(HIST("xi1530K/QAafter/kaonDCAz"), track.pt(), track.dcaZ());
        histos.fill(HIST("xi1530K/QAafter/kaonTPCNSigma"), track.pt(), track.tpcNSigmaKa());
        if (hasTOF) {
          histos.fill(HIST("xi1530K/QAafter/kaonTOFNSigma"), track.pt(), track.tofNSigmaKa());
        }
      }
    };

    const bool fillWrongSign = core.fillWrongSign();
    auto onCandidate = [&](auto const& xi, auto const& pion, auto const& kaon, Xi1530KCandidateValues const& c) {
      if (c.chargePattern != o2::analysis::omega2012ml::kSignalPattern) {
        // Charge-pattern controls: same event only, never in the signal histograms
        if constexpr (!IsMix) {
          if (fillWrongSign && c.inRapidity) {
            histos.fill(HIST("xi1530K_wrongSign/invmassPattern"), c.chargePattern, c.omega.M());
            histos.fill(HIST("xi1530K_wrongSign/massPtPattern"), c.omega.M(), c.omega.Pt(), c.chargePattern);
          }
        }
        return;
      }
      if constexpr (IsMix) {
        if (c.inRapidity) {
          histos.fill(HIST("xi1530K/invmass_Mix"), c.omega.M());
          histos.fill(HIST("xi1530K/massPtCent_Mix"), c.omega.M(), c.omega.Pt(), cent);
        }
        return;
      }
      histos.fill(HIST("xi1530K/massXiPi"), c.massXiPi);
      if (!c.inRapidity) {
        return;
      }
      histos.fill(HIST("xi1530K/invmass"), c.omega.M());
      histos.fill(HIST("xi1530K/massPtCent"), c.omega.M(), c.omega.Pt(), cent);
      histos.fill(HIST("xi1530K/massXiPiVsMass"), c.omega.M(), c.massXiPi);
      if constexpr (IsMC) {
        if (classifyXi1530KTruth(xi, pion, kaon) != Xi1530KTruth::Matched) {
          return;
        }
        histos.fill(HIST("xi1530K/MC/hMCRecOmega2012Pt"), c.omega.Pt());
        histos.fill(HIST("xi1530K/MC/hMCRecOmega2012PtEta"), c.omega.Pt(), c.omega.Eta());
        histos.fill(HIST("xi1530K/MC/hMCRecMass"), c.omega.M());
        histos.fill(HIST("xi1530K/MC/hMCRecMassXiPi"), c.massXiPi);
        histos.fill(HIST("xi1530K/MC/hMCRecXiPt"), xi.pt());
        histos.fill(HIST("xi1530K/MC/hMCRecPionPt"), pion.pt());
        histos.fill(HIST("xi1530K/MC/hMCRecKaonPt"), kaon.pt());
      }
    };

    core.forEachXi1530KCandidate<IsMC, IsMix, IsResoMicrotrack>(histos, collision, cascades, tracks, trackIds, false, onXi, onTrack, onCandidate);
  }

  void processXi1530KMicro(ResoCollisions::iterator const& collision,
                           aod::ResoCascades const& resocasc,
                           ResoMicroTracks const& resomicrotracks)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    fillXi1530K<false, false, true>(collision, collision.cent(), resocasc, resomicrotracks, nullptr);
  }
  PROCESS_SWITCH(Omega2012Analysis, processXi1530KMicro, "Process Xi(1530)0 K mode with ResoMicroTracks_001", false);

  void processXi1530KTracks(ResoCollisions::iterator const& collision,
                            aod::ResoCascades const& resocasc,
                            aod::ResoTracks const& resotracks,
                            aod::ResoTrackTracks const& resotrackids)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    fillXi1530K<false, false, false>(collision, collision.cent(), resocasc, resotracks, resotrackids);
  }
  PROCESS_SWITCH(Omega2012Analysis, processXi1530KTracks, "Process Xi(1530)0 K mode with ResoTracks", false);

  void processXi1530KMCMicro(ResoMCCollisions::iterator const& collision,
                             ResoMCCascades const& resocasc,
                             ResoMCMicroTracks const& resomicrotracks)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    fillXi1530K<true, false, true>(collision, collision.cent(), resocasc, resomicrotracks, nullptr);
  }
  PROCESS_SWITCH(Omega2012Analysis, processXi1530KMCMicro, "Process Xi(1530)0 K mode with truth matching (ResoMicroTracks_001)", false);

  // Mixed events: Xi from collision 1, pion and kaon from collision 2
  void processXi1530KMixedMicro(ResoCollisions const& collisions,
                                aod::ResoCascades const& resocasc,
                                ResoMicroTracks const& resomicrotracks)
  {
    auto cascTracksTuple = std::make_tuple(resocasc, resomicrotracks);
    Pair<ResoCollisions, aod::ResoCascades, ResoMicroTracks, BinningTypeVertexContributor> pairs{colBinning, nEvtMixing, -1, collisions, cascTracksTuple, &cache};
    for (const auto& [collision1, casc1, collision2, tracks2] : pairs) {
      if (!core.passesEventCuts(collision1) || !core.passesEventCuts(collision2)) {
        continue;
      }
      fillXi1530K<false, true, true>(collision1, collision1.cent(), casc1, tracks2, nullptr);
    }
  }
  PROCESS_SWITCH(Omega2012Analysis, processXi1530KMixedMicro, "Process mixed events of the Xi(1530)0 K mode (ResoMicroTracks_001)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Omega2012Analysis>(cfgc)};
}
