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

/// \file jetFragmentation.cxx
/// \brief Task for jet fragmentation into V0s
///
/// \author Gijs van Weelden <g.van.weelden@cern.ch>

#include "JetDerivedDataUtilities.h"

#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using DataV0JetsWithConstituents = soa::Join<aod::V0ChargedJets, aod::V0ChargedJetConstituents>;

using CandidatesV0MCDWithLabels = soa::Join<aod::CandidatesV0MCD, aod::McV0Labels>;
using MatchedMCDV0Jets = soa::Join<aod::V0ChargedMCDetectorLevelJets, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;
using MatchedMCDV0JetsWithConstituents = soa::Join<aod::V0ChargedMCDetectorLevelJets, aod::V0ChargedMCDetectorLevelJetConstituents, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;

using MatchedMCPV0Jets = soa::Join<aod::V0ChargedMCParticleLevelJets, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;
using MatchedMCPV0JetsWithConstituents = soa::Join<aod::V0ChargedMCParticleLevelJets, aod::V0ChargedMCParticleLevelJetConstituents, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;

struct JetFragmentation {
  HistogramRegistry registry{"registry"}; // CallSumw2 = false?

  Configurable<std::string> evSel{"evSel", "sel8WithoutTimeFrameBorderCut", "choose event selection"};
  Configurable<std::string> trackSel{"trackSel", "globalTracks", "choose track selection"};
  Configurable<int> nV0Classes{"nV0Classes", 2, "Must be 2 or 4! Number of V0 signal/bkg classes"};
  Configurable<bool> doCorrectionWithTracks{"doCorrectionWithTracks", false, "add tracks during background subtraction"};
  Configurable<bool> fillHistsInclusiveV0s{"fillHistsInclusiveV0s", true, "Fill hists for inclusive V0s"};
  Configurable<bool> fillHistsJets{"fillHistsJets", true, "Fill hists for jets"};

  Configurable<std::vector<float>> ptBinsK0S{"ptBinsK0S", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0}, "K0S pt Vals"};
  Configurable<std::vector<float>> ptBinsLambda{"ptBinsLambda", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0}, "Lambda pt Vals"};
  Configurable<std::vector<float>> ptBinsAntiLambda{"ptBinsAntiLambda", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0}, "AntiLambda pt Vals"};

  // NB: these must be one shorter than ptbin vectors!
  Configurable<std::vector<float>> signalProbK0S{"signalProbK0S", {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}, "K0S signal probability per pt bin"};
  Configurable<std::vector<float>> signalProbLambda{"signalProbLambda", {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}, "Lambda signal probability per pt bin"};
  Configurable<std::vector<float>> signalProbAntiLambda{"signalProbAntiLambda", {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}, "AntiLambda signal probability per pt bin"};

  Configurable<float> vertexZCut{"vertexZCut", 10.f, "vertex z cut"};
  Configurable<float> v0EtaMin{"v0EtaMin", -0.75, "minimum data V0 eta"};
  Configurable<float> v0EtaMax{"v0EtaMax", 0.75, "maximum data V0 eta"};

  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {40, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {20, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, constants::math::TwoPI}, ""};
  ConfigurableAxis binZ{"binZ", {40, 0.0001f, 1.0001f}, ""};
  ConfigurableAxis binXi{"binXi", {50, 0.f, 10.f}, ""};
  ConfigurableAxis binTheta{"binTheta", {40, -0.05f, 0.395f}, ""};
  ConfigurableAxis binJetR{"binJetR", {6, 0.05f, 0.65f}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binVtxZ{"binVtxZ", {200, -20, 20}, ""};

  ConfigurableAxis binPtTrackDiff{"binPtTrackDiff", {121, -20.5f, 100.5f}, ""};
  ConfigurableAxis binPtDiff{"binPtDiff", {600, -299.5f, 300.5f}, ""};
  ConfigurableAxis binEtaDiff{"binEtaDiff", {40, -0.195f, 0.205f}, ""};
  ConfigurableAxis binPhiDiff{"binPhiDiff", {40, -0.195f, 0.205f}, ""};
  ConfigurableAxis binZDiff{"binZDiff", {80, -1.05f, 0.95f}, ""};
  ConfigurableAxis binXiDiff{"binXiDiff", {100, -9.5f, 10.5f}, ""};
  ConfigurableAxis binThetaDiff{"binThetaDiff", {80, -0.35f, 0.45f}, ""};
  ConfigurableAxis binPtRatio{"binPtRatio", {50, -0.5f, 9.5f}, ""};   // Ratio of pt, eta, phi
  ConfigurableAxis binMatchDist{"binMatchDist", {50, 0.f, 0.5f}, ""}; // Distance between matched jets

  ConfigurableAxis binPtRelDiff{"binPtRelDiff", {100, -9.5f, 10.5f}, ""};
  ConfigurableAxis binZRelDiff{"binZRelDiff", {100, -9.5f, 10.5f}, ""};

  ConfigurableAxis binCount{"binCount", {1, .5f, 1.5f}, ""};
  ConfigurableAxis jetCount{"jetCount", {20, -.5f, 19.5f}, ""};
  ConfigurableAxis trackCount{"trackCount", {1000, -.5f, 999.5f}, ""};
  ConfigurableAxis v0Count{"v0Count", {50, -.5f, 49.5f}, ""};
  ConfigurableAxis v0Weight{"v0Weight", {50, 0.f, 10.0f}, ""};

  ConfigurableAxis binV0Pt{"binV0Pt", {120, 0.0f, 60.0f}, ""};
  ConfigurableAxis binV0Eta{"binV0Eta", {20, -1.f, 1.f}, ""};
  ConfigurableAxis binV0Phi{"binV0Phi", {18 * 8, 0.f, constants::math::TwoPI}, ""};
  ConfigurableAxis binV0Ctau{"binV0Ctau", {200, 0.0f, 40.0f}, ""};
  ConfigurableAxis binV0Radius{"binV0Radius", {100, 0.0f, 100.0f}, ""};
  ConfigurableAxis binV0CosPA{"binV0CosPA", {100, 0.95f, 1.0f}, ""};
  ConfigurableAxis binV0DCA{"binV0DCA", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis binV0DCAp{"binV0DCAp", {100, -10.0f, 10.0f}, ""};
  ConfigurableAxis binV0DCAn{"binV0DCAn", {100, -10.0f, 10.0f}, ""};
  ConfigurableAxis binV0DCAd{"binV0DCAd", {100, 0.0f, 10.0f}, ""};

  ConfigurableAxis binK0SMass{"binK0SMass", {400, 0.400f, 0.600f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis binK0SMassWide{"binK0SMassWide", {400, 0.400f, 0.800f}, "Inv. Mass (GeV/c^{2})"}; // Wider version for high pt
  ConfigurableAxis binLambdaMass{"binLambdaMass", {200, 1.075f, 1.215f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis binLambdaMassDiff{"binLambdaMassDiff", {200, -0.199f, 0.201f}, "M(#Lambda) - M(#bar{#Lambda})"};
  ConfigurableAxis binLambdaMassRatio{"binLambdaMassRatio", {50, -0.05f, 4.95f}, "M(#bar{#Lambda}) / M(#Lambda)"};
  ConfigurableAxis binLambdaMassRelDiff{"binLambdaMassRelDiff", {200, -0.995f, 1.005f}, "(M(#Lambda) - M(#bar{#Lambda})) / M(#Lambda)"};

  // Binning for cut variation study
  ConfigurableAxis binV0RadiusCut{"binV0RadiusCut", {4, 1.0f, 1.4f}, "R"};
  ConfigurableAxis binV0CtauCut{"binV0CtauCut", {3, 15.0f, 30.0f}, "c#tau"};
  ConfigurableAxis binV0CosPACut{"binV0CosPACut", {4, 0.991f, 0.999f}, "cosPA"};
  ConfigurableAxis binV0DCApCut{"binV0DCApCut", {2, 0.05f, 0.15f}, "DCA pos"};
  ConfigurableAxis binV0DCAnCut{"binV0DCAnCut", {2, 0.05f, 0.15f}, "DCA neg"};
  ConfigurableAxis binV0DCAdCut{"binV0DCAdCut", {2, 0.5f, 1.5f}, "DCA daughters"};
  ConfigurableAxis binV0PtCut{"binV0PtCut", {60, 0.0f, 60.0f}, "p_{T, V0}"};
  ConfigurableAxis binK0SMassCut{"binK0SMassCut", {100, 0.4f, 0.6f}, "inv. mass, K0S hypothesis"};
  ConfigurableAxis binLambdaMassCut{"binLambdaMassCut", {100, 1.07f, 1.21f}, "inv. mass, Lambda hypothesis"};
  ConfigurableAxis binAntiLambdaMassCut{"binAntiLambdaMassCut", {100, 1.07f, 1.21f}, "inv. mass, AntiLambda hypothesis"};

  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(evSel));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSel));

    // Axes
    AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T}^{ jet}"}; // Data
    AxisSpec etaAxis = {binEta, "#eta"};
    AxisSpec phiAxis = {binPhi, "#phi"};
    AxisSpec zAxis = {binZ, "#it{z}"};
    AxisSpec xiAxis = {binXi, "#xi"};
    AxisSpec thetaAxis = {binTheta, "#theta"};

    AxisSpec detJetPtAxis = {binJetPt, "#it{p}_{T}^{ jet, det}"}; // MC detector level
    AxisSpec detEtaAxis = {binEta, "#eta^{ jet, det}"};
    AxisSpec detPhiAxis = {binPhi, "#phi^{ jet, det}"};
    AxisSpec detZAxis = {binZ, "#it{z}^{ det}"};
    AxisSpec detXiAxis = {binXi, "#xi^{ det}"};
    AxisSpec detThetaAxis = {binTheta, "#theta^{ det}"};

    AxisSpec partJetPtAxis = {binJetPt, "#it{p}_{T}^{ jet, part}"}; // MC particle level
    AxisSpec partEtaAxis = {binEta, "#eta^{ jet, part}"};
    AxisSpec partPhiAxis = {binPhi, "#phi^{ jet, part}"};
    AxisSpec partZAxis = {binZ, "#it{z}^{ part}"};
    AxisSpec partXiAxis = {binXi, "#xi^{ part}"};
    AxisSpec partThetaAxis = {binTheta, "#theta^{ part}"};

    AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{tr}"};
    AxisSpec ptTrackDiffAxis = {binPtTrackDiff, "#it{p}_{T}^{track} - #it{p}_{T}^{particle}"};
    AxisSpec ptDiffAxis = {binPtDiff, "#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part}"};
    AxisSpec etaDiffAxis = {binEtaDiff, "#eta^{jet, det} - #eta^{jet, part}"};
    AxisSpec phiDiffAxis = {binPhiDiff, "#varphi^{jet, det} - #varphi^{jet, part}"};
    AxisSpec zDiffAxis = {binZDiff, "#it{z}^{det} - #it{z}^{part}"};
    AxisSpec xiDiffAxis = {binXiDiff, "#xi^{det} - #xi^{part}"};
    AxisSpec thetaDiffAxis = {binThetaDiff, "#theta^{det} - #theta^{part}"};
    AxisSpec ptRatioAxis = {binPtRatio, ""};
    AxisSpec vtxZAxis = {binVtxZ, "Collision vertex z (cm)"};
    AxisSpec matchDistAxis = {binMatchDist, "#Delta"};

    AxisSpec ptJetRelDiffAxis = {binPtRelDiff, "(#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part})/#it{p}_{T, jet}^{part}"};
    AxisSpec ptTrackRelDiffAxis = {binPtRelDiff, "(#it{p}_{T}^{track, det} - #it{p}_{T}^{track, part})/#it{p}_{T, track}^{part}"};
    AxisSpec zRelDiffAxis = {binZRelDiff, "(#it{z}^{det} - #it{z}^{part})/#it{z}^{part}"};

    AxisSpec v0PtAxis = {binV0Pt, "#it{p}_{T}^{V0}"};
    AxisSpec v0PtRatioAxis = {binPtRatio, "#it{p}_{T}^{V0, det}/#it{p}_{T, V0}^{part}"};
    AxisSpec v0PtRelDiffAxis = {binPtRelDiff, "(#it{p}_{T}^{V0, det} - #it{p}_{T}^{V0, part})/#it{p}_{T, V0}^{part}"};
    AxisSpec v0EtaAxis = {binV0Eta, "#eta^{V0}"};
    AxisSpec v0PhiAxis = {binV0Phi, "#varphi^{V0}"};
    AxisSpec v0detPtAxis = {binV0Pt, "#it{p}_{T}^{V0, det}"};
    AxisSpec v0partPtAxis = {binV0Pt, "#it{p}_{T}^{V0, part}"};
    AxisSpec v0CtauAxis = {binV0Ctau, "c#tau (cm)"};
    AxisSpec v0RadiusAxis = {binV0Radius, "R (cm)"};
    AxisSpec v0CosPAAxis = {binV0CosPA, "cos(PA)"};
    AxisSpec v0DCApAxis = {binV0DCAp, "DCA pos (cm)"};
    AxisSpec v0DCAnAxis = {binV0DCAn, "DCA neg (cm)"};
    AxisSpec v0DCAdAxis = {binV0DCAd, "DCA daughters (cm^{2})"};

    AxisSpec k0SMassAxis = {binK0SMass, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec k0SWideAxis = {binK0SMassWide, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec lambdaMassAxis = {binLambdaMass, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec lambdaMassDiffAxis = {binLambdaMassDiff, "M(#Lambda) - M(#bar{#Lambda})"};
    AxisSpec lambdaMassRatioAxis = {binLambdaMassRatio, "M(#bar{#Lambda}) / M(#Lambda)"};
    AxisSpec lambdaMassRelDiffAxis = {binLambdaMassRelDiff, "(M(#Lambda) - M(#bar{#Lambda})) / M(#Lambda)"};

    // Cut variation study
    AxisSpec rCutAxis = {binV0RadiusCut, "R"};
    AxisSpec ctauCutAxis = {binV0CtauCut, "c#tau (K0S)"};
    AxisSpec cosPACutAxis = {binV0CosPACut, "cosPA"};
    AxisSpec dcapCutAxis = {binV0DCApCut, "DCA pos (cm)"};
    AxisSpec dcanCutAxis = {binV0DCAnCut, "DCA neg (cm)"};
    AxisSpec dcadCutAxis = {binV0DCAdCut, "DCA daughters (cm^{2})"};
    AxisSpec ptCutAxis = {binV0PtCut, "p_{T, V0}"};
    AxisSpec k0SMassCutAxis = {binK0SMassCut, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec lambdaMassCutAxis = {binLambdaMassCut, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec antiLambdaMassCutAxis = {binAntiLambdaMassCut, "Inv. mass (GeV/#it{c}^{2})"};

    if (doprocessDataV0) {
      registry.add("data/hEvents", "hEvents", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("data/V0/nV0sEvent", "nV0sEvent", HistType::kTH1D, {v0Count});
      registry.add("data/V0/nV0sEventAcc", "nV0s per event (accepted)", HistType::kTH1D, {v0Count});
      registry.add("data/V0/nV0sEventAccWeighted", "nV0s per event (accepted, weighted)", HistType::kTH1D, {v0Weight}, true);

      // Inclusive
      registry.add("data/V0/V0PtEtaPhi", "V0PtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/V0/V0PtCtau", "V0PtCtau", HistType::kTHnSparseD, {v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis});
      registry.add("data/V0/V0PtMass", "V0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/V0/V0PtMassWide", "V0PtMassWide", HistType::kTHnSparseD, {v0PtAxis, k0SWideAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/V0/V0PtLambdaMasses", "V0PtLambdaMasses", HistType::kTHnSparseD, {v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/V0/V0PtRadiusCosPA", "V0PtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/V0/V0PtDCAposneg", "V0PtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/V0/V0PtDCAd", "V0PtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis});

      // Inclusive Weighted
      registry.add("data/V0/V0PtEtaPhiWeighted", "V0PtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("data/V0/V0PtCtauWeighted", "V0PtCtau", HistType::kTHnSparseD, {v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("data/V0/V0PtMassWeighted", "V0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("data/V0/V0PtMassWideWeighted", "V0PtMassWide", HistType::kTHnSparseD, {v0PtAxis, k0SWideAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("data/V0/V0PtLambdaMassesWeighted", "V0PtLambdaMasses", HistType::kTHnSparseD, {v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("data/V0/V0PtRadiusCosPAWeighted", "V0PtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("data/V0/V0PtDCAposnegWeighted", "V0PtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/V0/V0PtDCAdWeighted", "V0PtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis}, true);

      // K0S
      registry.add("data/V0/K0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/V0/K0SPtRadiusCosPA", "K0SPtRadiusCosPA", HistType::kTH3D, {v0partPtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/V0/K0SPtDCAposneg", "K0SPtDCAposneg", HistType::kTH3D, {v0partPtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/V0/K0SPtDCAd", "K0SPtDCAd", HistType::kTH2D, {v0partPtAxis, v0DCAdAxis});
      registry.add("data/V0/K0SPtCtauMass", "K0SPtCtauMass", HistType::kTH3D, {v0partPtAxis, v0CtauAxis, k0SMassAxis});
      registry.add("data/V0/K0SPtRadiusMass", "K0SPtRadiusMass", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, k0SMassAxis});
      registry.add("data/V0/K0SPtCosPAMass", "K0SPtCosPAMass", HistType::kTH3D, {v0PtAxis, v0CosPAAxis, k0SMassAxis});
      registry.add("data/V0/K0SPtDCAposMass", "K0SPtDCAposMass", HistType::kTH3D, {v0PtAxis, v0DCApAxis, k0SMassAxis});
      registry.add("data/V0/K0SPtDCAnegMass", "K0SPtDCAnegMass", HistType::kTH3D, {v0PtAxis, v0DCAnAxis, k0SMassAxis});
      registry.add("data/V0/K0SPtDCAdMass", "K0SPtDCAdMass", HistType::kTH3D, {v0PtAxis, v0DCAdAxis, k0SMassAxis});

      // K0S Weighted
      registry.add("data/V0/K0SPtEtaPhiWeighted", "K0SPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("data/V0/K0SPtRadiusCosPAWeighted", "K0SPtRadiusCosPA", HistType::kTH3D, {v0partPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("data/V0/K0SPtDCAposnegWeighted", "K0SPtDCAposneg", HistType::kTH3D, {v0partPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/V0/K0SPtDCAdWeighted", "K0SPtDCAd", HistType::kTH2D, {v0partPtAxis, v0DCAdAxis}, true);
      registry.add("data/V0/K0SPtCtauMassWeighted", "K0SPtCtauMass", HistType::kTH3D, {v0partPtAxis, v0CtauAxis, k0SMassAxis}, true);
      registry.add("data/V0/K0SPtRadiusMassWeighted", "K0SPtRadiusMassWeighted", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, k0SMassAxis}, true);
      registry.add("data/V0/K0SPtCosPAMassWeighted", "K0SPtCosPAMassWeighted", HistType::kTH3D, {v0PtAxis, v0CosPAAxis, k0SMassAxis}, true);
      registry.add("data/V0/K0SPtDCAposMassWeighted", "K0SPtDCAposMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCApAxis, k0SMassAxis}, true);
      registry.add("data/V0/K0SPtDCAnegMassWeighted", "K0SPtDCAnegMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCAnAxis, k0SMassAxis}, true);
      registry.add("data/V0/K0SPtDCAdMassWeighted", "K0SPtDCAdMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCAdAxis, k0SMassAxis}, true);

      // Lambda
      registry.add("data/V0/LambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/V0/LambdaPtLambdaMasses", "LambdaPtLambdaMasses", HistType::kTHnSparseD, {v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/V0/LambdaPtRadiusCosPA", "LambdaPtRadiusCosPA", HistType::kTH3D, {v0partPtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/V0/LambdaPtDCAposneg", "LambdaPtDCAposneg", HistType::kTH3D, {v0partPtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/V0/LambdaPtDCAd", "LambdaPtDCAd", HistType::kTH2D, {v0partPtAxis, v0DCAdAxis});
      registry.add("data/V0/LambdaPtCtauMass", "LambdaPtCtauMass", HistType::kTH3D, {v0partPtAxis, v0CtauAxis, lambdaMassAxis});
      registry.add("data/V0/LambdaPtRadiusMass", "LambdaPtRadiusMass", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, lambdaMassAxis});
      registry.add("data/V0/LambdaPtCosPAMass", "LambdaPtCosPAMass", HistType::kTH3D, {v0PtAxis, v0CosPAAxis, lambdaMassAxis});
      registry.add("data/V0/LambdaPtDCAposMass", "LambdaPtDCAposMass", HistType::kTH3D, {v0PtAxis, v0DCApAxis, lambdaMassAxis});
      registry.add("data/V0/LambdaPtDCAnegMass", "LambdaPtDCAnegMass", HistType::kTH3D, {v0PtAxis, v0DCAnAxis, lambdaMassAxis});
      registry.add("data/V0/LambdaPtDCAdMass", "LambdaPtDCAdMass", HistType::kTH3D, {v0PtAxis, v0DCAdAxis, lambdaMassAxis});

      // Lambda Weighted
      registry.add("data/V0/LambdaPtEtaPhiWeighted", "LambdaPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("data/V0/LambdaPtLambdaMassesWeighted", "LambdaPtLambdaMasses", HistType::kTHnSparseD, {v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("data/V0/LambdaPtRadiusCosPAWeighted", "LambdaPtRadiusCosPA", HistType::kTH3D, {v0partPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("data/V0/LambdaPtDCAposnegWeighted", "LambdaPtDCAposneg", HistType::kTH3D, {v0partPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/V0/LambdaPtDCAdWeighted", "LambdaPtDCAd", HistType::kTH2D, {v0partPtAxis, v0DCAdAxis}, true);
      registry.add("data/V0/LambdaPtCtauMassWeighted", "LambdaPtCtauMass", HistType::kTH3D, {v0partPtAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("data/V0/LambdaPtRadiusMassWeighted", "LambdaPtRadiusMassWeighted", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("data/V0/LambdaPtCosPAMassWeighted", "LambdaPtCosPAMassWeighted", HistType::kTH3D, {v0PtAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("data/V0/LambdaPtDCAposMassWeighted", "LambdaPtDCAposMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("data/V0/LambdaPtDCAnegMassWeighted", "LambdaPtDCAnegMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("data/V0/LambdaPtDCAdMassWeighted", "LambdaPtDCAdMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCAdAxis, lambdaMassAxis}, true);

      // AntiLambda
      registry.add("data/V0/AntiLambdaPtEtaPhi", "AntiLambdaPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/V0/AntiLambdaPtLambdaMasses", "AntiLambdaPtLambdaMasses", HistType::kTHnSparseD, {v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/V0/AntiLambdaPtRadiusCosPA", "AntiLambdaPtRadiusCosPA", HistType::kTH3D, {v0partPtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/V0/AntiLambdaPtDCAposneg", "AntiLambdaPtDCAposneg", HistType::kTH3D, {v0partPtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/V0/AntiLambdaPtDCAd", "AntiLambdaPtDCAd", HistType::kTH2D, {v0partPtAxis, v0DCAdAxis});
      registry.add("data/V0/AntiLambdaPtCtauMass", "AntiLambdaPtCtauMass", HistType::kTH3D, {v0partPtAxis, v0CtauAxis, lambdaMassAxis});
      registry.add("data/V0/AntiLambdaPtRadiusMass", "AntiLambdaPtRadiusMass", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, lambdaMassAxis});
      registry.add("data/V0/AntiLambdaPtCosPAMass", "AntiLambdaPtCosPAMass", HistType::kTH3D, {v0PtAxis, v0CosPAAxis, lambdaMassAxis});
      registry.add("data/V0/AntiLambdaPtDCAposMass", "AntiLambdaPtDCAposMass", HistType::kTH3D, {v0PtAxis, v0DCApAxis, lambdaMassAxis});
      registry.add("data/V0/AntiLambdaPtDCAnegMass", "AntiLambdaPtDCAnegMass", HistType::kTH3D, {v0PtAxis, v0DCAnAxis, lambdaMassAxis});
      registry.add("data/V0/AntiLambdaPtDCAdMass", "AntiLambdaPtDCAdMass", HistType::kTH3D, {v0PtAxis, v0DCAdAxis, lambdaMassAxis});

      // AntiLambda Weighted
      registry.add("data/V0/AntiLambdaPtEtaPhiWeighted", "AntiLambdaPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("data/V0/AntiLambdaPtLambdaMassesWeighted", "AntiLambdaPtLambdaMasses", HistType::kTHnSparseD, {v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("data/V0/AntiLambdaPtRadiusCosPAWeighted", "AntiLambdaPtRadiusCosPA", HistType::kTH3D, {v0partPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("data/V0/AntiLambdaPtDCAposnegWeighted", "AntiLambdaPtDCAposneg", HistType::kTH3D, {v0partPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/V0/AntiLambdaPtDCAdWeighted", "AntiLambdaPtDCAd", HistType::kTH2D, {v0partPtAxis, v0DCAdAxis}, true);
      registry.add("data/V0/AntiLambdaPtCtauMassWeighted", "AntiLambdaPtCtauMass", HistType::kTH3D, {v0partPtAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("data/V0/AntiLambdaPtRadiusMassWeighted", "AntiLambdaPtRadiusMassWeighted", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("data/V0/AntiLambdaPtCosPAMassWeighted", "AntiLambdaPtCosPAMassWeighted", HistType::kTH3D, {v0PtAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("data/V0/AntiLambdaPtDCAposMassWeighted", "AntiLambdaPtDCAposMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("data/V0/AntiLambdaPtDCAnegMassWeighted", "AntiLambdaPtDCAnegMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("data/V0/AntiLambdaPtDCAdMassWeighted", "AntiLambdaPtDCAdMassWeighted", HistType::kTH3D, {v0PtAxis, v0DCAdAxis, lambdaMassAxis}, true);

      // Cut variation
      registry.add("data/V0/V0CutVariation", "V0CutVariation", HistType::kTHnSparseD, {ptCutAxis, k0SMassCutAxis, lambdaMassCutAxis, antiLambdaMassCutAxis, rCutAxis, ctauCutAxis, cosPACutAxis, dcapCutAxis, dcanCutAxis, dcadCutAxis});

      // Jets
      registry.add("data/jets/jetPtEtaPhi", "Jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {jetPtAxis, etaAxis, phiAxis});
      registry.add("data/jets/inclJetPtEtaPhi", "incljetPtEtaPhi", HistType::kTH3D, {jetPtAxis, etaAxis, phiAxis});
      registry.add("data/jets/V0/jetPtV0TrackProj", "jetPtV0TrackProj", HistType::kTH2D, {jetPtAxis, zAxis});
      registry.add("data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda", "jetPtnV0nK0SnLambdanAntiLambda", HistType::kTHnSparseD, {jetPtAxis, v0Count, v0Count, v0Count, v0Count});

      // Inclusive
      registry.add("data/jets/V0/jetPtV0PtEtaPhi", "jetPtV0PtEtaPhi", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/jets/V0/jetPtV0PtCtau", "jetPtV0PtCtau", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis});
      registry.add("data/jets/V0/jetPtV0PtMass", "jetPtV0PtMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtV0PtMassWide", "jetPtV0PtMassWide", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, k0SWideAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtV0PtLambdaMasses", "jetPtV0PtLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtV0PtRadiusCosPA", "jetPtV0PtRadiusCosPA", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/jets/V0/jetPtV0PtDCAposneg", "jetPtV0PtDCAposneg", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/jets/V0/jetPtV0PtDCAd", "jetPtV0PtDCAd", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0DCAdAxis});

      registry.add("data/jets/V0/jetPtV0TrackProjCtau", "jetPtV0TrackProjCtau", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjMass", "jetPtV0TrackProjMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjMassWide", "jetPtV0TrackProjMassWide", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SWideAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjLambdaMasses", "jetPtV0TrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjRadiusCosPA", "jetPtV0TrackProjRadiusCosPA", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjDCAposneg", "jetPtV0TrackProjDCAposneg", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjDCAd", "jetPtV0TrackProjDCAd", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis});

      // K0S
      registry.add("data/jets/V0/jetPtK0SPtCtau", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, c#tau", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0CtauAxis});
      registry.add("data/jets/V0/jetPtK0SPtMass", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, mass", HistType::kTH3D, {jetPtAxis, v0PtAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtAllMasses", "jetPtK0SPtAllMasses", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtRadius", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, radius", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0RadiusAxis});
      registry.add("data/jets/V0/jetPtK0SPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, cosPA", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0CosPAAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, DCA daughters", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0DCAdAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/jets/V0/jetPtK0SPtCtauMass", "jetPtK0SPtCtauMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0CtauAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtRadiusMass", "jetPtK0SPtRadiusMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0RadiusAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtCosPAMass", "jetPtK0SPtCosPAMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0CosPAAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAposMass", "jetPtK0SPtDCAposMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCApAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAnegMass", "jetPtK0SPtDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCAnAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAdMass", "jetPtK0SPtDCAdMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCAdAxis, k0SMassAxis});

      registry.add("data/jets/V0/jetPtK0STrackProjCtau", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, v0CtauAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjMass", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjAllMasses", "jetPtK0STrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjRadius", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, v0RadiusAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, v0CosPAAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjCtauMass", "jetPtK0STrackProjCtauMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjRadiusMass", "jetPtK0STrackProjRadiusMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjCosPAMass", "jetPtK0STrackProjCosPAMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAposMass", "jetPtK0STrackProjDCAposMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAnegMass", "jetPtK0STrackProjDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, k0SMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAdMass", "jetPtK0STrackProjDCAdMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, k0SMassAxis});

      // Lambda
      registry.add("data/jets/V0/jetPtLambdaPtCtau", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0CtauAxis});
      registry.add("data/jets/V0/jetPtLambdaPtMass", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, mass", HistType::kTH3D, {jetPtAxis, v0PtAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtAllMasses", "jetPtLambdaPtAllMasses", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtLambdaMasses", "jetPtLambdaPtLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtLambdaPtRadius", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, radius", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0RadiusAxis});
      registry.add("data/jets/V0/jetPtLambdaPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0CosPAAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0DCAdAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/jets/V0/jetPtLambdaPtCtauMass", "jetPtLambdaPtCtauMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0CtauAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtRadiusMass", "jetPtLambdaPtRadiusMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0RadiusAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtCosPAMass", "jetPtLambdaPtCosPAMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0CosPAAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAposMass", "jetPtLambdaPtDCAposMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCApAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAnegMass", "jetPtLambdaPtDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCAnAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAdMass", "jetPtLambdaPtDCAdMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCAdAxis, lambdaMassAxis});

      registry.add("data/jets/V0/jetPtLambdaTrackProjCtau", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, v0CtauAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjMass", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjAllMasses", "jetPtLambdaTrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjLambdaMasses", "jetPtLambdaTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjRadius", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, v0RadiusAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, v0CosPAAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjCtauMass", "jetPtLambdaTrackProjCtauMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjRadiusMass", "jetPtLambdaTrackProjRadiusMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjCosPAMass", "jetPtLambdaTrackProjCosPAMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAposMass", "jetPtLambdaTrackProjDCAposMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAnegMass", "jetPtLambdaTrackProjDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAdMass", "jetPtLambdaTrackProjDCAdMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, lambdaMassAxis});

      // AntiLambda
      registry.add("data/jets/V0/jetPtAntiLambdaPtCtau", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0CtauAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtMass", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, mass", HistType::kTH3D, {jetPtAxis, v0PtAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtAllMasses", "jetPtAntiLambdaPtAllMasses", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtLambdaMasses", "jetPtAntiLambdaPtLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtRadius", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, radius", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0RadiusAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0CosPAAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, v0PtAxis, v0DCAdAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtCtauMass", "jetPtAntiLambdaPtCtauMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0CtauAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtRadiusMass", "jetPtAntiLambdaPtRadiusMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0RadiusAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtCosPAMass", "jetPtAntiLambdaPtCosPAMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0CosPAAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAposMass", "jetPtAntiLambdaPtDCAposMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCApAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAnegMass", "jetPtAntiLambdaPtDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCAnAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAdMass", "jetPtAntiLambdaPtDCAdMass", HistType::kTHnSparseD, {jetPtAxis, v0PtAxis, v0DCAdAxis, lambdaMassAxis});

      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjCtau", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, v0CtauAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjMass", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjAllMasses", "jetPtAntiLambdaTrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjLambdaMasses", "jetPtAntiLambdaTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjRadius", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, v0RadiusAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, v0CosPAAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjCtauMass", "jetPtAntiLambdaTrackProjCtauMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjRadiusMass", "jetPtAntiLambdaTrackProjRadiusMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjCosPAMass", "jetPtAntiLambdaTrackProjCosPAMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAposMass", "jetPtAntiLambdaTrackProjDCAposMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAnegMass", "jetPtAntiLambdaTrackProjDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, lambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAdMass", "jetPtAntiLambdaTrackProjDCAdMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, lambdaMassAxis});

      // Weighted histograms
      registry.add("data/jets/weighted/jetPtEtaPhi", "Jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {jetPtAxis, etaAxis, phiAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtnV0nK0SnLambdanAntiLambda", "jetPtnV0nK0SnLambdanAntiLambda", HistType::kTHnSparseD, {jetPtAxis, v0Weight, v0Weight, v0Weight, v0Weight}, true);

      // Inclusive Weighted
      registry.add("data/jets/weighted/V0/jetPtV0TrackProjCtau", "jetPtV0TrackProjCtau", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtV0TrackProjMass", "jetPtV0TrackProjMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtV0TrackProjMassWide", "jetPtV0TrackProjMassWide", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SWideAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtV0TrackProjLambdaMasses", "jetPtV0TrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtV0TrackProjRadiusCosPA", "jetPtV0TrackProjRadiusCosPA", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtV0TrackProjDCAposneg", "jetPtV0TrackProjDCAposneg", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtV0TrackProjDCAd", "jetPtV0TrackProjDCAd", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis}, true);

      // K0S Weighted
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjMass", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjCtau", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, v0CtauAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjAllMasses", "jetPtK0STrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjRadius", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, v0RadiusAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, v0CosPAAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjCtauMass", "jetPtK0STrackProjCtauMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjRadiusMass", "jetPtK0STrackProjRadiusMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjCosPAMass", "jetPtK0STrackProjCosPAMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjDCAposMass", "jetPtK0STrackProjDCAposMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjDCAnegMass", "jetPtK0STrackProjDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtK0STrackProjDCAdMass", "jetPtK0STrackProjDCAdMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, k0SMassAxis}, true);

      // Lambda Weighted
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjMass", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjCtau", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, v0CtauAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjAllMasses", "jetPtLambdaTrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjLambdaMasses", "jetPtLambdaTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjRadius", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, v0RadiusAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, v0CosPAAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjCtauMass", "jetPtLambdaTrackProjCtauMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjRadiusMass", "jetPtLambdaTrackProjRadiusMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjCosPAMass", "jetPtLambdaTrackProjCosPAMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjDCAposMass", "jetPtLambdaTrackProjDCAposMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjDCAnegMass", "jetPtLambdaTrackProjDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtLambdaTrackProjDCAdMass", "jetPtLambdaTrackProjDCAdMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, lambdaMassAxis}, true);

      // AntiLambda Weighted
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjMass", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjCtau", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, v0CtauAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjAllMasses", "jetPtAntiLambdaTrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjLambdaMasses", "jetPtAntiLambdaTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjRadius", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, v0RadiusAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, v0CosPAAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjCtauMass", "jetPtAntiLambdaTrackProjCtauMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjRadiusMass", "jetPtAntiLambdaTrackProjRadiusMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjCosPAMass", "jetPtAntiLambdaTrackProjCosPAMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAposMass", "jetPtAntiLambdaTrackProjDCAposMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAnegMass", "jetPtAntiLambdaTrackProjDCAnegMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAdMass", "jetPtAntiLambdaTrackProjDCAdMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, lambdaMassAxis}, true);

      // Background Weighted
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjMass", "jetPtBkgTrackProjMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjCtau", "jetPtBkgTrackProjCtau", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjLambdaMasses", "jetPtBkgTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjRadiusCosPA", "jetPtBkgTrackProjRadiusCosPA", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAposneg", "jetPtBkgTrackProjDCAposneg", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAd", "jetPtBkgTrackProjDCAd", HistType::kTH3D, {jetPtAxis, zAxis, v0DCAdAxis}, true);

      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjCtauK0SMass", "jetPtBkgTrackProjCtauK0SMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjRadiusK0SMass", "jetPtBkgTrackProjRadiusK0SMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjCosPAK0SMass", "jetPtBkgTrackProjCosPAK0SMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAposK0SMass", "jetPtBkgTrackProjDCAposK0SMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAnegK0SMass", "jetPtBkgTrackProjDCAnegK0SMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, k0SMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAdK0SMass", "jetPtBkgTrackProjDCAdK0SMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, k0SMassAxis}, true);

      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjCtauLambdaMass", "jetPtBkgTrackProjCtauLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjRadiusLambdaMass", "jetPtBkgTrackProjRadiusLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjCosPALambdaMass", "jetPtBkgTrackProjCosPALambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAposLambdaMass", "jetPtBkgTrackProjDCAposLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAnegLambdaMass", "jetPtBkgTrackProjDCAnegLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAdLambdaMass", "jetPtBkgTrackProjDCAdLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, lambdaMassAxis}, true);

      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjCtauAntiLambdaMass", "jetPtBkgTrackProjCtauAntiLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjRadiusAntiLambdaMass", "jetPtBkgTrackProjRadiusAntiLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjCosPAAntiLambdaMass", "jetPtBkgTrackProjCosPAAntiLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAposAntiLambdaMass", "jetPtBkgTrackProjDCAposAntiLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAnegAntiLambdaMass", "jetPtBkgTrackProjDCAnegAntiLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("data/jets/weighted/V0/jetPtBkgTrackProjDCAdAntiLambdaMass", "jetPtBkgTrackProjDCAdAntiLambdaMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, v0DCAdAxis, lambdaMassAxis}, true);

      // Perp Cone
      registry.add("data/PC/JetPtEtaV0Pt", "JetPtEtaV0Pt", HistType::kTH3D, {jetPtAxis, etaAxis, v0PtAxis});
      registry.add("data/PC/V0PtEtaPhi", "V0 #it{p}_{T}, #eta, #phi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/PC/V0PtCtau", "V0PtCtau", HistType::kTHnSparseD, {v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis});
      registry.add("data/PC/V0PtMass", "V0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/PC/V0PtMassWide", "V0PtMassWide", HistType::kTHnSparseD, {v0PtAxis, k0SWideAxis, lambdaMassAxis, lambdaMassAxis});
      registry.add("data/PC/V0PtRadiusCosPA", "V0PtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/PC/V0PtDCAposneg", "V0PtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/PC/V0PtDCAd", "V0PtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis});

      // K0S
      registry.add("data/PC/JetPtK0SPtMass", "JetPtK0SPtMass", HistType::kTH3D, {jetPtAxis, v0PtAxis, k0SMassAxis});
      registry.add("data/PC/JetPtEtaK0SPt", "JetPtEtaK0SPt", HistType::kTH3D, {jetPtAxis, etaAxis, v0PtAxis});
      registry.add("data/PC/K0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/PC/K0SPtCtauMass", "K0SPtCtauMass", HistType::kTH3D, {v0PtAxis, v0CtauAxis, k0SMassAxis});
      registry.add("data/PC/K0SPtRadiusCosPA", "K0SPtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/PC/K0SPtDCAposneg", "K0SPtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/PC/K0SPtDCAd", "K0SPtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis});

      // Lambda
      registry.add("data/PC/JetPtLambdaPtMass", "JetPtLambdaPtMass", HistType::kTH3D, {jetPtAxis, v0PtAxis, lambdaMassAxis});
      registry.add("data/PC/JetPtEtaLambdaPt", "JetPtEtaLambdaPt", HistType::kTH3D, {jetPtAxis, etaAxis, v0PtAxis});
      registry.add("data/PC/LambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/PC/LambdaPtCtauMass", "LambdaPtCtauMass", HistType::kTH3D, {v0PtAxis, v0CtauAxis, lambdaMassAxis});
      registry.add("data/PC/LambdaPtRadiusCosPA", "LambdaPtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/PC/LambdaPtDCAposneg", "LambdaPtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/PC/LambdaPtDCAd", "LambdaPtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis});

      // AntiLambda
      registry.add("data/PC/JetPtAntiLambdaPtMass", "JetPtAntiLambdaPtMass", HistType::kTH3D, {jetPtAxis, v0PtAxis, lambdaMassAxis});
      registry.add("data/PC/JetPtEtaAntiLambdaPt", "JetPtEtaAntiLambdaPt", HistType::kTH3D, {jetPtAxis, etaAxis, v0PtAxis});
      registry.add("data/PC/AntiLambdaPtEtaPhi", "AntiLambdaPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis});
      registry.add("data/PC/AntiLambdaPtCtauMass", "AntiLambdaPtCtauMass", HistType::kTH3D, {v0PtAxis, v0CtauAxis, lambdaMassAxis});
      registry.add("data/PC/AntiLambdaPtRadiusCosPA", "AntiLambdaPtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis});
      registry.add("data/PC/AntiLambdaPtDCAposneg", "AntiLambdaPtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis});
      registry.add("data/PC/AntiLambdaPtDCAd", "AntiLambdaPtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis});

      registry.add("data/PC/nV0sConePtEta", "nV0sConePtEta", HistType::kTH3D, {v0Count, jetPtAxis, etaAxis});
      registry.add("data/PC/ConePtEtaPhi", "ConePtEtaPhi", HistType::kTH3D, {jetPtAxis, etaAxis, phiAxis});
      registry.add("data/PC/JetPtEtaConePt", "JetPtEtaConePt", HistType::kTH3D, {jetPtAxis, etaAxis, jetPtAxis});
    } // doprocessDataV0

    if (doprocessMcV0) {
      registry.add("matching/hEvents", "hEvents", {HistType::kTH1D, {{3, 0.0f, 3.0f}}});
      // MCP
      registry.add("mcp/V0/nV0sEventAcc", "nV0s per event (accepted)", HistType::kTH1D, {v0Count}, true);
      registry.add("mcp/V0/nV0sEventAccWeighted", "nV0s per event (accepted, weighted)", HistType::kTH1D, {v0Weight}, true);
      registry.add("mcp/V0/V0PtEtaPhi", "V0PtEtaPhi", HistType::kTH3D, {v0partPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcp/V0/K0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTH3D, {v0partPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcp/V0/LambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTH3D, {v0partPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcp/V0/AntiLambdaPtEtaPhi", "AntiLambdaPtEtaPhi", HistType::kTH3D, {v0partPtAxis, v0EtaAxis, v0PhiAxis}, true);

      // MCD Inclusive
      registry.add("mcd/V0/nV0sEventAcc", "nV0s per event (accepted)", HistType::kTH1D, {v0Count});
      registry.add("mcd/V0/nV0sEventAccWeighted", "nV0s per event (accepted, weighted)", HistType::kTH1D, {v0Weight}, true);
      registry.add("mcd/V0/V0PtEtaPhi", "V0PtEtaPhi", HistType::kTH3D, {v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcd/V0/V0PtCtau", "V0PtCtau", HistType::kTHnSparseD, {v0detPtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("mcd/V0/V0PtMass", "V0PtMass", HistType::kTHnSparseD, {v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/V0PtLambdaMasses", "V0PtLambdaMasses", HistType::kTHnSparseD, {v0detPtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("mcd/V0/V0PtRadiusCosPA", "V0PtRadiusCosPA", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("mcd/V0/V0PtDCAposneg", "V0PtDCAposneg", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("mcd/V0/V0PtDCAd", "V0PtDCAd", HistType::kTH2D, {v0detPtAxis, v0DCAdAxis}, true);

      // MCD K0S
      registry.add("mcd/V0/K0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTH3D, {v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcd/V0/K0SPtRadiusCosPA", "K0SPtRadiusCosPA", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("mcd/V0/K0SPtDCAposneg", "K0SPtDCAposneg", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("mcd/V0/K0SPtDCAd", "K0SPtDCAd", HistType::kTH2D, {v0detPtAxis, v0DCAdAxis}, true);
      registry.add("mcd/V0/K0SPtCtauMass", "K0SPtCtauMass", HistType::kTH3D, {v0detPtAxis, v0CtauAxis, k0SMassAxis}, true);
      registry.add("mcd/V0/K0SPtRadiusMass", "K0SPtRadiusMass", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, k0SMassAxis}, true);
      registry.add("mcd/V0/K0SPtCosPAMass", "K0SPtCosPAMass", HistType::kTH3D, {v0detPtAxis, v0CosPAAxis, k0SMassAxis}, true);
      registry.add("mcd/V0/K0SPtDCAposMass", "K0SPtDCAposMass", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, k0SMassAxis}, true);
      registry.add("mcd/V0/K0SPtDCAnegMass", "K0SPtDCAnegMass", HistType::kTH3D, {v0detPtAxis, v0DCAnAxis, k0SMassAxis}, true);
      registry.add("mcd/V0/K0SPtDCAdMass", "K0SPtDCAdMass", HistType::kTH3D, {v0detPtAxis, v0DCAdAxis, k0SMassAxis}, true);

      // MCD Lambda
      registry.add("mcd/V0/LambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTH3D, {v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcd/V0/LambdaPtLambdaMasses", "LambdaPtLambdaMasses", HistType::kTHnSparseD, {v0detPtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("mcd/V0/LambdaPtRadiusCosPA", "LambdaPtRadiusCosPA", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("mcd/V0/LambdaPtDCAposneg", "LambdaPtDCAposneg", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("mcd/V0/LambdaPtDCAd", "LambdaPtDCAd", HistType::kTH2D, {v0detPtAxis, v0DCAdAxis}, true);
      registry.add("mcd/V0/LambdaPtCtauMass", "LambdaPtCtauMass", HistType::kTH3D, {v0detPtAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/LambdaPtRadiusMass", "LambdaPtRadiusMass", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/LambdaPtCosPAMass", "LambdaPtCosPAMass", HistType::kTH3D, {v0detPtAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/LambdaPtDCAposMass", "LambdaPtDCAposMass", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/LambdaPtDCAnegMass", "LambdaPtDCAnegMass", HistType::kTH3D, {v0detPtAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/LambdaPtDCAdMass", "LambdaPtDCAdMass", HistType::kTH3D, {v0detPtAxis, v0DCAdAxis, lambdaMassAxis}, true);

      // MCD AntiLambda
      registry.add("mcd/V0/AntiLambdaPtEtaPhi", "AntiLambdaPtEtaPhi", HistType::kTH3D, {v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtLambdaMasses", "AntiLambdaPtLambdaMasses", HistType::kTHnSparseD, {v0detPtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtRadiusCosPA", "AntiLambdaPtRadiusCosPA", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtDCAposneg", "AntiLambdaPtDCAposneg", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtDCAd", "AntiLambdaPtDCAd", HistType::kTH2D, {v0detPtAxis, v0DCAdAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtCtauMass", "AntiLambdaPtCtauMass", HistType::kTH3D, {v0detPtAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtRadiusMass", "AntiLambdaPtRadiusMass", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtCosPAMass", "AntiLambdaPtCosPAMass", HistType::kTH3D, {v0detPtAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtDCAposMass", "AntiLambdaPtDCAposMass", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtDCAnegMass", "AntiLambdaPtDCAnegMass", HistType::kTH3D, {v0detPtAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("mcd/V0/AntiLambdaPtDCAdMass", "AntiLambdaPtDCAdMass", HistType::kTH3D, {v0detPtAxis, v0DCAdAxis, lambdaMassAxis}, true);

      // Matching - Misses
      registry.add("matching/V0/missV0PtEtaPhi", "missV0PtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/missK0SPtEtaPhi", "missK0SPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/missLambdaPtEtaPhi", "missLambdaPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/missAntiLambdaPtEtaPhi", "missAntiLambdaPtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);

      // Matching - Fakes Inclusive
      registry.add("matching/V0/fakeV0PtEtaPhi", "fakeV0PtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/fakeV0PtCtau", "fakeV0PtCtau", HistType::kTHnSparseD, {v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/V0/fakeV0PtMass", "fakeV0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeV0PtLambdaMasses", "fakeV0PtLambdaMasses", HistType::kTHnSparseD, {v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/V0/fakeV0PtRadiusCosPA", "fakeV0PtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/V0/fakeV0PtDCAposneg", "fakeV0PtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/V0/fakeV0PtDCAd", "fakeV0PtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis}, true);
      registry.add("matching/V0/fakeV0PosTrackPtEtaPhi", "fakeV0PosTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/V0/fakeV0NegTrackPtEtaPhi", "fakeV0NegTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis}, true);

      // Matching - Fakes K0S
      registry.add("matching/V0/fakeK0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTH3D, {v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/fakeK0SPtRadiusCosPA", "K0SPtRadiusCosPA", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/V0/fakeK0SPtDCAposneg", "K0SPtDCAposneg", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/V0/fakeK0SPtDCAd", "K0SPtDCAd", HistType::kTH2D, {v0detPtAxis, v0DCAdAxis}, true);
      registry.add("matching/V0/fakeK0SPtCtauMass", "K0SPtCtauMass", HistType::kTH3D, {v0detPtAxis, v0CtauAxis, k0SMassAxis}, true);
      registry.add("matching/V0/fakeK0SPtRadiusMass", "K0SPtRadiusMass", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, k0SMassAxis}, true);
      registry.add("matching/V0/fakeK0SPtCosPAMass", "K0SPtCosPAMass", HistType::kTH3D, {v0detPtAxis, v0CosPAAxis, k0SMassAxis}, true);
      registry.add("matching/V0/fakeK0SPtDCAposMass", "K0SPtDCAposMass", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, k0SMassAxis}, true);
      registry.add("matching/V0/fakeK0SPtDCAnegMass", "K0SPtDCAnegMass", HistType::kTH3D, {v0detPtAxis, v0DCAnAxis, k0SMassAxis}, true);
      registry.add("matching/V0/fakeK0SPtDCAdMass", "K0SPtDCAdMass", HistType::kTH3D, {v0detPtAxis, v0DCAdAxis, k0SMassAxis}, true);
      registry.add("matching/V0/fakeK0SPosTrackPtEtaPhi", "fakeK0SPosTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/V0/fakeK0SPosTrackPtMass", "fakeK0SPosTrackPtMass", HistType::kTH3D, {v0detPtAxis, trackPtAxis, k0SMassAxis}, true);
      registry.add("matching/V0/fakeK0SNegTrackPtEtaPhi", "fakeK0SNegTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/V0/fakeK0SNegTrackPtMass", "fakeK0SNegTrackPtMass", HistType::kTH3D, {v0detPtAxis, trackPtAxis, k0SMassAxis}, true);

      // Matching - Fakes Lambda
      registry.add("matching/V0/fakeLambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTH3D, {v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/fakeLambdaPtLambdaMasses", "LambdaPtLambdaMasses", HistType::kTHnSparseD, {v0detPtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/V0/fakeLambdaPtRadiusCosPA", "LambdaPtRadiusCosPA", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/V0/fakeLambdaPtDCAposneg", "LambdaPtDCAposneg", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/V0/fakeLambdaPtDCAd", "LambdaPtDCAd", HistType::kTH2D, {v0detPtAxis, v0DCAdAxis}, true);
      registry.add("matching/V0/fakeLambdaPtCtauMass", "LambdaPtCtauMass", HistType::kTH3D, {v0detPtAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeLambdaPtRadiusMass", "LambdaPtRadiusMass", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeLambdaPtCosPAMass", "LambdaPtCosPAMass", HistType::kTH3D, {v0detPtAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeLambdaPtDCAposMass", "LambdaPtDCAposMass", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeLambdaPtDCAnegMass", "LambdaPtDCAnegMass", HistType::kTH3D, {v0detPtAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeLambdaPtDCAdMass", "LambdaPtDCAdMass", HistType::kTH3D, {v0detPtAxis, v0DCAdAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeLambdaPosTrackPtEtaPhi", "fakeLambdaPosTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/V0/fakeLambdaPosTrackPtMass", "fakeLambdaPosTrackPtMass", HistType::kTH3D, {v0detPtAxis, trackPtAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeLambdaNegTrackPtEtaPhi", "fakeLambdaNegTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/V0/fakeLambdaNegTrackPtMass", "fakeLambdaNegTrackPtMass", HistType::kTH3D, {v0detPtAxis, trackPtAxis, lambdaMassAxis}, true);

      // Matching - Fakes AntiLambda
      registry.add("matching/V0/fakeAntiLambdaPtEtaPhi", "AntiLambdaPtEtaPhi", HistType::kTH3D, {v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtLambdaMasses", "AntiLambdaPtLambdaMasses", HistType::kTHnSparseD, {v0detPtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtRadiusCosPA", "AntiLambdaPtRadiusCosPA", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtDCAposneg", "AntiLambdaPtDCAposneg", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtDCAd", "AntiLambdaPtDCAd", HistType::kTH2D, {v0detPtAxis, v0DCAdAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtCtauMass", "AntiLambdaPtCtauMass", HistType::kTH3D, {v0detPtAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtRadiusMass", "AntiLambdaPtRadiusMass", HistType::kTH3D, {v0detPtAxis, v0RadiusAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtCosPAMass", "AntiLambdaPtCosPAMass", HistType::kTH3D, {v0detPtAxis, v0CosPAAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtDCAposMass", "AntiLambdaPtDCAposMass", HistType::kTH3D, {v0detPtAxis, v0DCApAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtDCAnegMass", "AntiLambdaPtDCAnegMass", HistType::kTH3D, {v0detPtAxis, v0DCAnAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPtDCAdMass", "AntiLambdaPtDCAdMass", HistType::kTH3D, {v0detPtAxis, v0DCAdAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPosTrackPtEtaPhi", "fakeAntiLambdaPosTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaPosTrackPtMass", "fakeAntiLambdaPosTrackPtMass", HistType::kTH3D, {v0detPtAxis, trackPtAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaNegTrackPtEtaPhi", "fakeAntiLambdaNegTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/V0/fakeAntiLambdaNegTrackPtMass", "fakeAntiLambdaNegTrackPtMass", HistType::kTH3D, {v0detPtAxis, trackPtAxis, lambdaMassAxis}, true);

      // Matching - Matched
      registry.add("matching/V0/nV0sEvent", "nV0sDet per event", HistType::kTH1D, {v0Count});
      registry.add("matching/V0/nV0sEventWeighted", "nV0sDet per event (weighted)", HistType::kTH1D, {v0Count}, true);
      registry.add("matching/V0/nV0sEventAcc", "nV0sDet per event (accepted, matched)", HistType::kTH1D, {v0Count}, true);
      registry.add("matching/V0/nV0sEventAccWeighted", "nV0sDet per event (accepted, matched, weighted)", HistType::kTH1D, {v0Weight}, true);

      // Matching - Matched Inclusive
      registry.add("matching/V0/V0PartPtDetPt", "V0PartPtDetPt", HistType::kTH2D, {v0partPtAxis, v0detPtAxis}, true);
      registry.add("matching/V0/V0PartPtRatioPtRelDiffPt", "V0PartPtRatioRelDiffPt", HistType::kTH3D, {v0partPtAxis, v0PtRatioAxis, v0PtRelDiffAxis}, true);

      // Matching - Matched K0S
      registry.add("matching/V0/K0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/K0SPtCtauMass", "K0SPtCtauMass", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0CtauAxis, k0SMassAxis}, true);
      registry.add("matching/V0/K0SPtRadiusCosPA", "K0SPtRadiusCosPA", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/V0/K0SPtDCAposneg", "K0SPtDCAposneg", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/V0/K0SPtDCAd", "K0SPtDCAd", HistType::kTH3D, {v0partPtAxis, v0detPtAxis, v0DCAdAxis}, true);
      registry.add("matching/V0/K0SPtMass", "K0SPtMass", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);

      // Matching - Matched Lambda
      registry.add("matching/V0/LambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/LambdaPtCtauMass", "LambdaPtCtauMass", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/LambdaPtRadiusCosPA", "LambdaPtRadiusCosPA", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/V0/LambdaPtDCAposneg", "LambdaPtDCAposneg", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/V0/LambdaPtDCAd", "LambdaPtDCAd", HistType::kTH3D, {v0partPtAxis, v0detPtAxis, v0DCAdAxis}, true);
      registry.add("matching/V0/LambdaPtMass", "LambdaPtMass", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/LambdaReflection", "pt, pt, mK, mL, maL, LambdaReflection", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis, lambdaMassAxis}, true);

      // Matching - Matched AntiLambda
      registry.add("matching/V0/AntiLambdaPtEtaPhi", "AntiLambdaPtEtaPhi", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/V0/AntiLambdaPtCtauMass", "AntiLambdaPtCtauMass", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0CtauAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/AntiLambdaPtRadiusCosPA", "AntiLambdaPtRadiusCosPA", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/V0/AntiLambdaPtDCAposneg", "AntiLambdaPtDCAposneg", HistType::kTHnSparseD, {v0partPtAxis, v0partPtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/V0/AntiLambdaPtDCAd", "AntiLambdaPtDCAd", HistType::kTH3D, {v0partPtAxis, v0detPtAxis, v0DCAdAxis}, true);
      registry.add("matching/V0/AntiLambdaPtMass", "AntiLambdaPtMass", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/AntiLambdaReflection", "pt, pt, mK, mL, maL, AntiLambdaReflection", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis, lambdaMassAxis}, true);

      // Matching - Matched Daughters
      registry.add("matching/V0/V0PosPartPtRatioPtRelDiffPt", "V0PosPartPtRatioRelDiffPt", HistType::kTH3D, {trackPtAxis, ptRatioAxis, ptTrackRelDiffAxis}, true);
      registry.add("matching/V0/V0NegPartPtRatioPtRelDiffPt", "V0NegPartPtRatioRelDiffPt", HistType::kTH3D, {trackPtAxis, ptRatioAxis, ptTrackRelDiffAxis}, true);
      registry.add("matching/V0/K0SPosNegPtMass", "K0SPosNegPtMass", HistType::kTHnSparseD, {v0partPtAxis, trackPtAxis, trackPtAxis, k0SMassAxis}, true);
      registry.add("matching/V0/LambdaPosNegPtMass", "LambdaPosNegPtMass", HistType::kTHnSparseD, {v0partPtAxis, trackPtAxis, trackPtAxis, lambdaMassAxis}, true);
      registry.add("matching/V0/AntiLambdaPosNegPtMass", "AntiLambdaPosNegPtMass", HistType::kTHnSparseD, {v0partPtAxis, trackPtAxis, trackPtAxis, lambdaMassAxis}, true);

      // Matching - Jets
      registry.add("mcp/jets/inclPartJetPtEtaPhi", "inclPartJetPtEtaPhi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis}, true);
      registry.add("mcp/jets/partJetPtEtaPhi", "Particle level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis}, true);
      registry.add("mcd/jets/inclDetJetPtEtaPhi", "inclDetJetPtEtaPhi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis}, true);
      registry.add("mcd/jets/detJetPtEtaPhi", "Detector level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis}, true);

      registry.add("matching/jets/fakeJetPtEtaPhi", "fakeJetPtEtaPhi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis}, true);
      registry.add("matching/jets/inclFakeJetPtEtaPhi", "inclFakeJetPtEtaPhi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis}, true);
      registry.add("matching/jets/inclMissJetPtEtaPhi", "inclMissJetPtEtaPhi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis}, true);
      registry.add("matching/jets/missJetPtEtaPhi", "missJetPtEtaPhi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis}, true);

      registry.add("matching/jets/matchDetJetPtEtaPhi", "Matched detector level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis}, true);
      registry.add("matching/jets/matchPartJetPtEtaPhi", "Matched particle level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis}, true);
      registry.add("matching/jets/matchPartJetPtEtaPhiMatchDist", "matchJetMatchDist", HistType::kTHnSparseD, {partJetPtAxis, partEtaAxis, partPhiAxis, matchDistAxis}, true);
      registry.add("matching/jets/matchPartJetPtEnergyScale", "jetEnergyScale", HistType::kTH2D, {partJetPtAxis, ptRatioAxis}, true);
      registry.add("matching/jets/matchDetJetPtPartJetPt", "matchDetJetPtPartJetPt", HistType::kTH2D, {detJetPtAxis, partJetPtAxis}, true);
      registry.add("matching/jets/matchPartJetPtDetJetEtaPartJetEta", "matchPartJetPtDetJetEtaPartJetEta", HistType::kTH3D, {partJetPtAxis, detEtaAxis, partEtaAxis}, true);
      registry.add("matching/jets/matchPartJetPtDetJetPhiPartJetPhi", "matchPartJetPtDetJetPhiPartJetPhi", HistType::kTH3D, {partJetPtAxis, detPhiAxis, partPhiAxis}, true);
      registry.add("matching/jets/matchPartJetPtResolutionPt", "#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part}", HistType::kTH2D, {partJetPtAxis, ptDiffAxis}, true);
      registry.add("matching/jets/matchPartJetPtResolutionEta", "#eta^{jet, det} - #eta^{jet, part}", HistType::kTH3D, {partJetPtAxis, partEtaAxis, etaDiffAxis}, true);
      registry.add("matching/jets/matchPartJetPtResolutionPhi", "#phi^{jet, det} - #phi^{jet, part}", HistType::kTH3D, {partJetPtAxis, partPhiAxis, phiDiffAxis}, true);
      registry.add("matching/jets/matchPartJetPtRelDiffPt", "#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part}", HistType::kTH2D, {partJetPtAxis, ptJetRelDiffAxis}, true);

      registry.add("matching/jets/inclMatchDetJetPtEtaPhi", "inclMatchDetJetPtEtaPhi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtEtaPhi", "inclMatchPartJetPtEtaPhi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtEtaPhiMatchDist", "matchJetMatchDist", HistType::kTHnSparseD, {partJetPtAxis, partEtaAxis, partPhiAxis, matchDistAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtEnergyScale", "jetEnergyScale", HistType::kTH2D, {partJetPtAxis, ptRatioAxis}, true);
      registry.add("matching/jets/inclMatchDetJetPtPartJetPt", "matchDetJetPtPartJetPt", HistType::kTH2D, {detJetPtAxis, partJetPtAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtDetJetEtaPartJetEta", "matchPartJetPtDetJetEtaPartJetEta", HistType::kTH3D, {partJetPtAxis, detEtaAxis, partEtaAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtDetJetPhiPartJetPhi", "matchPartJetPtDetJetPhiPartJetPhi", HistType::kTH3D, {partJetPtAxis, detPhiAxis, partPhiAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtResolutionPt", "inclMatchPartJetPtResolutionPt", HistType::kTH2D, {partJetPtAxis, ptDiffAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtResolutionEta", "inclMatchPartJetPtResolutionEta", HistType::kTH3D, {partJetPtAxis, partEtaAxis, etaDiffAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtResolutionPhi", "inclMatchPartJetPtResolutionPhi", HistType::kTH3D, {partJetPtAxis, partPhiAxis, phiDiffAxis}, true);
      registry.add("matching/jets/inclMatchPartJetPtRelDiffPt", "inclMatchPartJetPtRelDiffPt", HistType::kTH2D, {partJetPtAxis, ptJetRelDiffAxis}, true);

      registry.add("matching/jets/V0/jetPtnV0MatchednK0SnLambdanAntiLambda", "jet Pt, nV0 matched, nK0S nLambdan AntiLambda", HistType::kTHnSparseD, {detJetPtAxis, v0Count, v0Count, v0Count, v0Count}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetPt", "V0PartPtDetPt", HistType::kTH3D, {partJetPtAxis, v0partPtAxis, v0detPtAxis}, true);
      // registry.add("matching/jets/V0/partJetPtDetJetPtPartV0PtRatioPtRelDiffPt", "V0PartPtRatioRelDiffPt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0partPtAxis, v0PtRatioAxis, v0PtRelDiffAxis}, true);
      registry.add("matching/jets/V0/partJetPtDetJetPtPartV0PtRelDiffPt", "V0PartPtRelDiffPt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0partPtAxis, v0PtRelDiffAxis}, true);
      registry.add("matching/jets/V0/partJetPtDetJetPtPartV0ZRelDiffZ", "V0PartPtRelDiffZ", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, partZAxis, zRelDiffAxis}, true);

      // Matching - Matched V0s
      registry.add("matching/jets/V0/matchDetJetPtV0TrackProjPartJetPtV0TrackProj", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0Pt", "matched jet Pt, V0 Pt", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      // Matching - Matched V0s: pt
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauK0S", "matched jet Pt, V0 Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauLambda", "matched jet Pt, V0 Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauAntiLambda", "matched jet Pt, V0 Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassK0S", "matched jet Pt, V0 Pt, MassK0S", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassLambda", "matched jet Pt, V0 Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassAntiLambda", "matched jet Pt, V0 Pt, Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtRadius", "matched jet Pt, V0 Pt, Radius", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0RadiusAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCosPA", "matched jet Pt, V0 Pt, CosPA", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtDCAposneg", "matched jet Pt, V0 Pt, DCAposneg", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtDCAd", "matched jet Pt, V0 Pt, DCAd", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0DCAdAxis}, true);
      // Matching - Matched V0s: z
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauK0S", "matched jet Pt, V0 z, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauLambda", "matched jet Pt, V0 z, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauAntiLambda", "matched jet Pt, V0 z, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassK0S", "matched jet Pt, V0 z, MassK0S", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassLambda", "matched jet Pt, V0 z, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassAntiLambda", "matched jet Pt, V0 z, Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjRadius", "matched jet Pt, V0 z, Radius", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0RadiusAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCosPA", "matched jet Pt, V0 z, CosPA", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjDCAposneg", "matched jet Pt, V0 z, DCAposneg", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjDCAd", "matched jet Pt, V0 z, DCAd", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0DCAdAxis}, true);

      // Matching - Fake V0s: pt
      registry.add("matching/jets/V0/fakeJetPtV0PtEtaPhi", "fake jet Pt, V0 PtEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0TrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0PtCtau", "fake jet Pt, V0 PtCtau", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0PtMass", "fake jet Pt, V0 PtMass", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0PtLambdaMasses", "fake jet Pt, V0 PtLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0PtRadiusCosPA", "fake jet Pt, V0 PtRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0PtDCAposneg", "fake jet Pt, V0 PtDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0PtDCAd", "fake jet Pt, V0 PtDCAd", HistType::kTH3D, {detJetPtAxis, v0PtAxis, v0DCAdAxis}, true);
      // Matching - Fake V0s: z
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjCtau", "fake jet Pt, V0 zCtau", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjMass", "fake jet Pt, V0 zMass", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjLambdaMasses", "fake jet Pt, V0 zLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjRadiusCosPA", "fake jet Pt, V0 zRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjDCAposneg", "fake jet Pt, V0 zDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjDCAd", "fake jet Pt, V0 zDCAd", HistType::kTH3D, {detJetPtAxis, detZAxis, v0DCAdAxis}, true);
      // Matching - Missed V0s
      registry.add("matching/jets/V0/missJetPtV0PtEtaPhi", "miss jet Pt, V0 PtEtaPhi", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/missJetPtV0TrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis}, true);

      // Matching - Matched K0S
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPt", "Matched Pt and pt", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProj", "Matched pt and z", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtRightCollision", "Matched Pt, right collision", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjRightCollision", "Matched pt and z, right collision", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtWrongCollision", "Matched Pt, wrong collision", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjWrongCollision", "Matched pt and z, wrong collision", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      // Matching - Matched K0S: pt
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauLambda", "matched jet Pt, K^{0}_{S} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauAntiLambda", "matched jet Pt, K^{0}_{S} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassLambda", "matched jet Pt, K^{0}_{S} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassAntiLambda", "matched jet Pt, K^{0}_{S} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassK0S", "matched jet Pt, K^{0}_{S} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtAllMasses", "matched jet Pt, K^{0}_{S} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtRadius", "matched jet Pt, K^{0}_{S} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0RadiusAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCosPA", "matched jet Pt, K^{0}_{S} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtDCAposneg", "matched jet Pt, K^{0}_{S} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtDCAd", "matched jet Pt, K^{0}_{S} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0DCAdAxis}, true);
      // Matching - Matched K0S: z
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauLambda", "matched jet Pt, K^{0}_{S} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauAntiLambda", "matched jet Pt, K^{0}_{S} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassLambda", "matched jet Pt, K^{0}_{S} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassAntiLambda", "matched jet Pt, K^{0}_{S} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassK0S", "matched jet Pt, K^{0}_{S} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjAllMasses", "matched jet Pt, K^{0}_{S} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjRadius", "matched jet Pt, K^{0}_{S} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0RadiusAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCosPA", "matched jet Pt, K^{0}_{S} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjDCAposneg", "matched jet Pt, K^{0}_{S} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjDCAd", "matched jet Pt, K^{0}_{S} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0DCAdAxis}, true);
      // Matching - Fake K0S: pt
      registry.add("matching/jets/V0/fakeJetPtK0SPtEtaPhi", "fake jet Pt, K^{0}_{S} PtEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0STrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0SPtCtau", "fake jet Pt, K^{0}_{S} PtCtau", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0SPtMass", "fake jet Pt, K^{0}_{S} PtMass", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0SPtLambdaMasses", "fake jet Pt, K^{0}_{S} PtLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0SPtRadiusCosPA", "fake jet Pt, K^{0}_{S} PtRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0SPtDCAposneg", "fake jet Pt, K^{0}_{S} PtDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0SPtDCAd", "fake jet Pt, K^{0}_{S} PtDCAd", HistType::kTH3D, {detJetPtAxis, v0PtAxis, v0DCAdAxis}, true);
      // Matching - Fake K0S: z
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjCtau", "fake jet Pt, K^{0}_{S} zCtau", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjMass", "fake jet Pt, K^{0}_{S} zMass", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjLambdaMasses", "fake jet Pt, K^{0}_{S} zLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjRadiusCosPA", "fake jet Pt, K^{0}_{S} zRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjDCAposneg", "fake jet Pt, K^{0}_{S} zDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjDCAd", "fake jet Pt, K^{0}_{S} zDCAd", HistType::kTH3D, {detJetPtAxis, detZAxis, v0DCAdAxis}, true);
      // Matching - Missed K0S
      registry.add("matching/jets/V0/missJetPtK0SPtEtaPhi", "miss jet Pt, K^{0}_{S} PtEtaPhi", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/missJetPtK0STrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis}, true);

      // Matching - Matched Lambda
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPt", "matched Pt", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProj", "Matched pt and z", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtRightCollision", "matched Pt, right collision", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjRightCollision", "Matched pt and z, right collision", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtWrongCollision", "matched Pt, Wrong collision", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjWrongCollision", "Matched pt and z, Wrong collision", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      // Matching - Matched Lambda: pt
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtCtauLambda", "matched jet Pt, #Lambda^{0} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtCtauAntiLambda", "matched jet Pt, #Lambda^{0} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtMassK0S", "matched jet Pt, #Lambda^{0} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtMassLambda", "matched jet Pt, #Lambda^{0} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtMassAntiLambda", "matched jet Pt, #Lambda^{0} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtAllMasses", "matched jet Pt, #Lambda^{0} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtRadius", "matched jet Pt, #Lambda^{0} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0RadiusAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtCosPA", "matched jet Pt, #Lambda^{0} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtDCAposneg", "matched jet Pt, #Lambda^{0} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtDCAd", "matched jet Pt, #Lambda^{0} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0DCAdAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtLambdaReflection", "Lambda Reflection", HistType::kTHnSparseD, {partJetPtAxis, v0partPtAxis, detJetPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Matched Lambda: z
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjCtauLambda", "matched jet Pt, #Lambda^{0} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjCtauAntiLambda", "matched jet Pt, #Lambda^{0} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjMassLambda", "matched jet Pt, #Lambda^{0} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjMassAntiLambda", "matched jet Pt, #Lambda^{0} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjMassK0S", "matched jet Pt, #Lambda^{0} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjAllMasses", "matched jet Pt, #Lambda^{0} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjRadius", "matched jet Pt, #Lambda^{0} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0RadiusAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjCosPA", "matched jet Pt, #Lambda^{0} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjDCAposneg", "matched jet Pt, #Lambda^{0} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjDCAd", "matched jet Pt, #Lambda^{0} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0DCAdAxis}, true);
      registry.add("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjLambdaReflection", "Lambda Reflection", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Fake Lambda: pt
      registry.add("matching/jets/V0/fakeJetPtLambdaPtEtaPhi", "fake jet Pt, #Lambda^{0} PtEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaTrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaPtCtau", "fake jet Pt, #Lambda^{0} PtCtau", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaPtMass", "fake jet Pt, #Lambda^{0} PtMass", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaPtLambdaMasses", "fake jet Pt, #Lambda^{0} PtLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaPtRadiusCosPA", "fake jet Pt, #Lambda^{0} PtRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaPtDCAposneg", "fake jet Pt, #Lambda^{0} PtDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaPtDCAd", "fake jet Pt, #Lambda^{0} PtDCAd", HistType::kTH3D, {detJetPtAxis, v0PtAxis, v0DCAdAxis}, true);
      // Matching - Fake Lambda: z
      registry.add("matching/jets/V0/fakeJetPtLambdaTrackProjEtaPhi", "fake jet Pt, #Lambda^{0} zEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaTrackProjCtau", "fake jet Pt, #Lambda^{0} zCtau", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaTrackProjMass", "fake jet Pt, #Lambda^{0} zMass", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaTrackProjLambdaMasses", "fake jet Pt, #Lambda^{0} zLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaTrackProjRadiusCosPA", "fake jet Pt, #Lambda^{0} zRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaTrackProjDCAposneg", "fake jet Pt, #Lambda^{0} zDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtLambdaTrackProjDCAd", "fake jet Pt, #Lambda^{0} zDCAd", HistType::kTH3D, {detJetPtAxis, detZAxis, v0DCAdAxis}, true);
      // Matching - Missed Lambda
      registry.add("matching/jets/V0/missJetPtLambdaPtEtaPhi", "miss jet Pt, #Lambda^{0} PtEtaPhi", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/missJetPtLambdaTrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis}, true);

      // Matching - Matched AntiLambda
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPt", "matched Pt", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProj", "Matched pt and z", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtRightCollision", "matched Pt, right collision", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjRightCollision", "Matched pt and z, right collision", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtWrongCollision", "matched Pt, Wrong collision", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjWrongCollision", "Matched pt and z, Wrong collision", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis}, true);

      // Matching - Matched AntiLambda: pt
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtCtauLambda", "matched jet Pt, #bar{#Lambda}^{0} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtCtauAntiLambda", "matched jet Pt, #bar{#Lambda}^{0} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtMassK0S", "matched jet Pt, #bar{#Lambda}^{0} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtMassLambda", "matched jet Pt, #bar{#Lambda}^{0} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtMassAntiLambda", "matched jet Pt, #bar{#Lambda}^{0} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtAllMasses", "matched jet Pt, #bar{#Lambda}^{0} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtRadius", "matched jet Pt, #bar{#Lambda}^{0} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0RadiusAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtCosPA", "matched jet Pt, #bar{#Lambda}^{0} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtDCAposneg", "matched jet Pt, #bar{#Lambda}^{0} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtDCAd", "matched jet Pt, #bar{#Lambda}^{0} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, detJetPtAxis, v0PtAxis, v0DCAdAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtAntiLambdaReflection", "AntiLambda Reflection", HistType::kTHnSparseD, {partJetPtAxis, v0partPtAxis, detJetPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Matched AntiLambda: z
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjCtauLambda", "matched jet Pt, #bar{#Lambda}^{0} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjCtauAntiLambda", "matched jet Pt, #bar{#Lambda}^{0} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjMassLambda", "matched jet Pt, #bar{#Lambda}^{0} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjMassAntiLambda", "matched jet Pt, #bar{#Lambda}^{0} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjMassK0S", "matched jet Pt, #bar{#Lambda}^{0} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjAllMasses", "matched jet Pt, #bar{#Lambda}^{0} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjRadius", "matched jet Pt, #bar{#Lambda}^{0} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0RadiusAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjCosPA", "matched jet Pt, #bar{#Lambda}^{0} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjDCAposneg", "matched jet Pt, #bar{#Lambda}^{0} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjDCAd", "matched jet Pt, #bar{#Lambda}^{0} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, v0DCAdAxis}, true);
      registry.add("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjAntiLambdaReflection", "AntiLambda Reflection", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Fake AntiLambda: pt
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaPtEtaPhi", "fake jet Pt, #bar{#Lambda}^{0} PtEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaTrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaPtCtau", "fake jet Pt, #bar{#Lambda}^{0} PtCtau", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaPtMass", "fake jet Pt, #bar{#Lambda}^{0} PtMass", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaPtLambdaMasses", "fake jet Pt, #bar{#Lambda}^{0} PtLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaPtRadiusCosPA", "fake jet Pt, #bar{#Lambda}^{0} PtRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaPtDCAposneg", "fake jet Pt, #bar{#Lambda}^{0} PtDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaPtDCAd", "fake jet Pt, #bar{#Lambda}^{0} PtDCAd", HistType::kTH3D, {detJetPtAxis, v0PtAxis, v0DCAdAxis}, true);
      // Matching - Fake AntiLambda: z
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaTrackProjCtau", "fake jet Pt, #bar{#Lambda}^{0} zCtau", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaTrackProjMass", "fake jet Pt, #bar{#Lambda}^{0} zMass", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaTrackProjLambdaMasses", "fake jet Pt, #bar{#Lambda}^{0} zLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, lambdaMassDiffAxis, lambdaMassRatioAxis, lambdaMassRelDiffAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaTrackProjRadiusCosPA", "fake jet Pt, #bar{#Lambda}^{0} zRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaTrackProjDCAposneg", "fake jet Pt, #bar{#Lambda}^{0} zDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/jets/V0/fakeJetPtAntiLambdaTrackProjDCAd", "fake jet Pt, #bar{#Lambda}^{0} zDCAd", HistType::kTH3D, {detJetPtAxis, detZAxis, v0DCAdAxis}, true);
      // Matching - Missed AntiLambda
      registry.add("matching/jets/V0/missJetPtAntiLambdaPtEtaPhi", "miss jet Pt, #bar{#Lambda}^{0} PtEtaPhi", HistType::kTHnSparseD, {partJetPtAxis, v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/jets/V0/missJetPtAntiLambdaTrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis}, true);

      // Matching - Matched V0s: daughters
      registry.add("matching/jets/V0/partJetPtDetJetPtPartV0PtPosPtRatioPtRelDiffPt", "V0PtPosPartPtRatioRelDiffPt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0PtAxis, trackPtAxis, ptRatioAxis, ptTrackRelDiffAxis}, true);
      registry.add("matching/jets/V0/partJetPtDetJetPtPartV0PtNegPtRatioPtRelDiffPt", "V0PtNegPartPtRatioRelDiffPt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0PtAxis, trackPtAxis, ptRatioAxis, ptTrackRelDiffAxis}, true);
      // Matching - Fake V0s: daughters non-decayed (should be empty)
      registry.add("matching/V0/nonedecayedFakeV0PtMass", "nonedecayedFakeV0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/nonedecayedFakeV0PtMass", "nonedecayedFakeV0PtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/nonedecayedFakeV0TrackProjMass", "nonedecayedFakeV0TrackProjMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Fake V0s: daughters both-decayed (seems unlikely: e.g. K->pi pi->mu nu mu nu)
      registry.add("matching/V0/doubledecayedFakeV0PtMass", "doubledecayedFakeV0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/doubledecayedFakeV0PtMass", "doubledecayedFakeV0PtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/doubledecayedFakeV0TrackProjMass", "doubledecayedFakeV0TrackProjMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, zAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Fake V0s: daughters one-decayed (e.g. K->pi pi->pi mu nu)
      // Matching - Fake K0S: daughters one-decayed
      registry.add("matching/V0/decayedK0SV0PtMass", "decayedK0SV0PtMass", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/decayedK0SV0PtMass", "decayedK0SV0PtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/decayedK0SV0TrackProjMass", "decayedK0SV0TrackProjMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, partZAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Fake Lambda: daughters one-decayed
      registry.add("matching/V0/decayedLambdaV0PtMass", "decayedLambdaV0PtMass", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/decayedLambdaV0PtMass", "decayedLambdaV0PtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/decayedLambdaV0TrackProjMass", "decayedLambdaV0TrackProjMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, partZAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Fake AntiLambda: daughters one-decayed
      registry.add("matching/V0/decayedAntiLambdaV0PtMass", "decayedAntiLambdaV0PtMass", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/decayedAntiLambdaV0PtMass", "decayedAntiLambdaV0PtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/decayedAntiLambdaV0TrackProjMass", "decayedAntiLambdaV0TrackProjMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, partZAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      // Matching - Fake Other: daughters one-decayed
      registry.add("matching/V0/decayedOtherPtV0PtMass", "decayedOtherPtV0PtMass", HistType::kTHnSparseD, {v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/decayedOtherPtV0PtMass", "decayedOtherPtV0PtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0partPtAxis, v0detPtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/jets/V0/decayedOtherPtV0TrackProjMass", "decayedOtherPtV0TrackProjMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, partZAxis, detZAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);

      registry.add("mcd/V0/nV0sEvent", "NV0s in event", HistType::kTH1D, {v0Count});
      registry.add("mcd/V0/nV0sEventWeighted", "NV0s in event weighted", HistType::kTH1D, {v0Count}, true);

      // PerpCone - Fake V0s in MCD jets
      registry.add("mcd/PC/jetPtEtaFakeV0Pt", "JetPtEtaFakeV0Pt", HistType::kTH3D, {detJetPtAxis, etaAxis, v0PtAxis}, true);
      registry.add("mcd/PC/fakeV0PtEtaPhi", "fakeV0PtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcd/PC/fakeV0PtCtau", "fakeV0PtCtau", HistType::kTHnSparseD, {v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("mcd/PC/fakeV0PtMass", "fakeV0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("mcd/PC/fakeV0PtRadiusCosPA", "fakeV0PtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("mcd/PC/fakeV0PtDCAposneg", "fakeV0PtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("mcd/PC/fakeV0PtDCAd", "fakeV0PtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis}, true);
      // PerpCone - Fake V0s in cone from MCD jets
      registry.add("mcd/PC/fakenV0sConePtEta", "fakenV0sConePtEta", HistType::kTH3D, {v0Count, detJetPtAxis, etaAxis}, true);
      registry.add("mcd/PC/fakeConePtEtaPhi", "fakeConePtEtaPhi", HistType::kTH3D, {detJetPtAxis, etaAxis, phiAxis}, true);
      registry.add("mcd/PC/fakeJetPtEtaConePt", "fakeJetPtEtaConePt", HistType::kTH3D, {detJetPtAxis, etaAxis, detJetPtAxis}, true);
      // PerpCone - Matched V0s in MCD jets
      registry.add("mcd/PC/jetPtEtaMatchedV0Pt", "JetPtEtaMatchedV0Pt", HistType::kTH3D, {detJetPtAxis, etaAxis, v0PtAxis}, true);
      registry.add("mcd/PC/matchedV0PtEtaPhi", "matchedV0PtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("mcd/PC/matchedV0PtCtau", "matchedV0PtCtau", HistType::kTHnSparseD, {v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("mcd/PC/matchedV0PtMass", "matchedV0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("mcd/PC/matchedV0PtRadiusCosPA", "matchedV0PtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("mcd/PC/matchedV0PtDCAposneg", "matchedV0PtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("mcd/PC/matchedV0PtDCAd", "matchedV0PtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis}, true);
      // PerpCone - Matched V0s in cone from MCD jets
      registry.add("mcd/PC/matchednV0sConePtEta", "matchednV0sConePtEta", HistType::kTH3D, {v0Count, detJetPtAxis, etaAxis}, true);
      registry.add("mcd/PC/matchedConePtEtaPhi", "matchedConePtEtaPhi", HistType::kTH3D, {detJetPtAxis, etaAxis, phiAxis}, true);
      registry.add("mcd/PC/matchedJetPtEtaConePt", "matchedJetPtEtaConePt", HistType::kTH3D, {detJetPtAxis, etaAxis, detJetPtAxis}, true);
      // PerpCone - Matched K0S in MCD jets
      registry.add("mcd/PC/matchedJetPtK0SPtMass", "matchedJetPtK0SPtMass", HistType::kTH3D, {detJetPtAxis, v0PtAxis, k0SMassAxis}, true);
      // PerpCone - Matched Lambda in MCD jets
      registry.add("mcd/PC/matchedJetPtLambdaPtMass", "matchedJetPtLambdaPtMass", HistType::kTH3D, {detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      // PerpCone - Matched AntiLambda in MCD jets
      registry.add("mcd/PC/matchedJetPtAntiLambdaPtMass", "matchedJetPtAntiLambdaPtMass", HistType::kTH3D, {detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      // PerpCone - Fake V0s in Matched jets
      registry.add("matching/PC/jetPtEtaFakeV0Pt", "JetPtEtaFakeV0Pt", HistType::kTH3D, {detJetPtAxis, etaAxis, v0PtAxis}, true);
      registry.add("matching/PC/jetsPtFakeV0Pt", "jetsPtFakeV0Pt", HistType::kTH3D, {partJetPtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/PC/fakeV0PtEtaPhi", "fakeV0PtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/PC/fakeV0PtCtau", "fakeV0PtCtau", HistType::kTHnSparseD, {v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/PC/fakeV0PtMass", "fakeV0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/PC/fakeV0PtRadiusCosPA", "fakeV0PtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/PC/fakeV0PtDCAposneg", "fakeV0PtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/PC/fakeV0PtDCAd", "fakeV0PtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis}, true);
      // PerpCone - Fake V0s in cone from Matched jets
      registry.add("matching/PC/fakenV0sConePtEta", "fakenV0sConePtEta", HistType::kTH3D, {v0Count, detJetPtAxis, etaAxis}, true);
      registry.add("matching/PC/fakeConePtEtaPhi", "fakeConePtEtaPhi", HistType::kTH3D, {detJetPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/PC/fakeJetPtEtaConePt", "fakeJetPtEtaConePt", HistType::kTH3D, {detJetPtAxis, etaAxis, detJetPtAxis}, true);
      registry.add("matching/PC/fakeJetsPtEtaConePt", "fakeJetsPtEtaConePt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, etaAxis, detJetPtAxis}, true);
      // PerpCone - Matched V0s in Matched jets
      registry.add("matching/PC/jetPtEtaMatchedV0Pt", "jetPtEtaMatchedV0Pt", HistType::kTH3D, {detJetPtAxis, etaAxis, v0PtAxis}, true);
      registry.add("matching/PC/jetsPtMatchedV0Pt", "jetsPtMatchedV0Pt", HistType::kTH3D, {partJetPtAxis, detJetPtAxis, v0PtAxis}, true);
      registry.add("matching/PC/matchedV0PtEtaPhi", "matchedV0PtEtaPhi", HistType::kTH3D, {v0PtAxis, v0EtaAxis, v0PhiAxis}, true);
      registry.add("matching/PC/matchedV0PtCtau", "matchedV0PtCtau", HistType::kTHnSparseD, {v0PtAxis, v0CtauAxis, v0CtauAxis, v0CtauAxis}, true);
      registry.add("matching/PC/matchedV0PtMass", "matchedV0PtMass", HistType::kTHnSparseD, {v0PtAxis, k0SMassAxis, lambdaMassAxis, lambdaMassAxis}, true);
      registry.add("matching/PC/matchedV0PtRadiusCosPA", "matchedV0PtRadiusCosPA", HistType::kTH3D, {v0PtAxis, v0RadiusAxis, v0CosPAAxis}, true);
      registry.add("matching/PC/matchedV0PtDCAposneg", "matchedV0PtDCAposneg", HistType::kTH3D, {v0PtAxis, v0DCApAxis, v0DCAnAxis}, true);
      registry.add("matching/PC/matchedV0PtDCAd", "matchedV0PtDCAd", HistType::kTH2D, {v0PtAxis, v0DCAdAxis}, true);
      // PerpCone - Matched V0s in cone from Matched jets
      registry.add("matching/PC/matchednV0sConePtEta", "matchednV0sConePtEta", HistType::kTH3D, {v0Count, detJetPtAxis, etaAxis}, true);
      registry.add("matching/PC/matchedConePtEtaPhi", "matchedConePtEtaPhi", HistType::kTH3D, {detJetPtAxis, etaAxis, phiAxis}, true);
      registry.add("matching/PC/matchedJetPtEtaConePt", "matchedJetPtEtaConePt", HistType::kTH3D, {detJetPtAxis, etaAxis, detJetPtAxis}, true);
      registry.add("matching/PC/matchedJetsPtEtaConePt", "matchedJetsPtEtaConePt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, etaAxis, detJetPtAxis}, true);
      // PerpCone - Matched K0S in Matched jets
      registry.add("matching/PC/matchedJetPtK0SPtMass", "matchedJetPtK0SPtMass", HistType::kTH3D, {detJetPtAxis, v0PtAxis, k0SMassAxis}, true);
      registry.add("matching/PC/matchedJetsPtK0SPtMass", "matchedJetsPtK0SPtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0PtAxis, k0SMassAxis}, true);
      // PerpCone - Matched Lambda in Matched jets
      registry.add("matching/PC/matchedJetPtLambdaPtMass", "matchedJetPtLambdaPtMass", HistType::kTH3D, {detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/PC/matchedJetsPtLambdaPtMass", "matchedJetsPtLambdaPtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      // PerpCone - Matched AntiLambda in Matched jets
      registry.add("matching/PC/matchedJetPtAntiLambdaPtMass", "matchedJetPtAntiLambdaPtMass", HistType::kTH3D, {detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
      registry.add("matching/PC/matchedJetsPtAntiLambdaPtMass", "matchedJetsPtAntiLambdaPtMass", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, v0PtAxis, lambdaMassAxis}, true);
    } // doprocessMcV0
  } // init

  // ---------------------------------------------------
  // Implementation of background subtraction at runtime
  // ---------------------------------------------------

  int getPtBin(float pt, const std::vector<float>& ptBins)
  {
    if (pt < ptBins.at(0))
      return -1;
    if (pt > ptBins.at(ptBins.size() - 1))
      return -2;

    for (unsigned int i = 0; i < ptBins.size() - 1; i++) {
      if (pt >= ptBins.at(i) && pt < ptBins.at(i + 1)) {
        return i;
      }
    }
    return -3;
  }

  // Return probability for a V0 to be signal
  // Assumes V0 can only be of one type!
  // Assumes V0 is not rejected!
  template <typename V>
  float getV0SignalProb(V const& v0)
  {
    double purity = 0.;
    if (v0.isK0SCandidate()) {
      int ptBin = getPtBin(v0.pt(), ptBinsK0S);
      if (ptBin >= 0) {
        purity = signalProbK0S->at(ptBin);
      }
    } else if (v0.isLambdaCandidate()) {
      int ptBin = getPtBin(v0.pt(), ptBinsLambda);
      if (ptBin >= 0) {
        purity = signalProbLambda->at(ptBin);
      }
    } else if (v0.isAntiLambdaCandidate()) {
      int ptBin = getPtBin(v0.pt(), ptBinsAntiLambda);
      if (ptBin >= 0) {
        purity = signalProbAntiLambda->at(ptBin);
      }
    }
    return purity;
  }
  // Return a 2-length std::vector of probabilities for a particle to correspond to signal or background
  template <typename V>
  std::vector<double> getV0SignalProbVector2Classes(V const& v0)
  {
    // 0: bkg, 1: signal
    if (v0.isRejectedCandidate())
      return {1., 0.};

    double purity = getV0SignalProb(v0);
    return {1. - purity, purity};
  }
  // Return a 4-length std::vector of probabilities for a particle to correspond to signal types
  template <typename V>
  std::vector<double> getV0SignalProbVector4Classes(V const& v0)
  {
    // 0: bkg, 1: K0S, 2: Lambda, 3: AntiLambda
    if (v0.isRejectedCandidate())
      return {1., 0., 0., 0.};

    double purity = getV0SignalProb(v0);
    if (v0.isK0SCandidate())
      return {1. - purity, purity, 0., 0.};
    if (v0.isLambdaCandidate())
      return {1. - purity, 0., purity, 0.};
    if (v0.isAntiLambdaCandidate())
      return {1. - purity, 0., 0., purity};

    return {1., 0., 0., 0.};
  }
  // Convert state from uint32_t to std::vector<int> containing the particle classes for that weight
  std::vector<int> convertState(uint32_t state, int nParticles, int nClasses = 4)
  {
    std::vector<int> v(nParticles, nClasses);
    int nStates = std::round(std::pow(static_cast<double>(nClasses), static_cast<double>(nParticles)));
    int nBitsPerParticle = std::round(std::log2(nClasses));
    int nBitsPerInt = sizeof(uint32_t) * 8;

    if (nClasses <= 0) {
      LOGF(warning, "Number of classes (%d) must be greater than 0", nClasses);
      return v;
    }
    if ((nClasses & (nClasses - 1)) != 0) {
      // It's likely possible to make this work for non-power of 2 classes, but it's not needed and therefore not implemented
      LOGF(warning, "Number of classes (%d) must be a power of 2", nClasses);
      return v;
    }
    // Check for overflow in number of states
    if (nStates <= 0) {
      LOGF(warning, "Illegal number of states (%d)! %s", nStates, (nStates == 0) ? "" : "Max = 2^31");
      return v;
    }
    // Check if we have enough bits to store the state
    if (nParticles * nBitsPerParticle > nBitsPerInt) {
      LOGF(warning, "Number of bits required to parse the state (%d * %d = %d) is too large for %d bits per int!", nParticles, nBitsPerParticle, nParticles * nBitsPerParticle, nBitsPerInt);
      return v;
    }
    // Check if the state is valid. This should never be triggered
    if (state >= static_cast<uint32_t>(nStates)) {
      LOGF(warning, "Illegal state! State %d >= %d", state, nStates);
      return v;
    }

    for (int ip = 0; ip < nParticles; ip++) {
      int value = 0;
      int startBit = ip * nBitsPerParticle;
      for (int ib = 0; ib < nBitsPerParticle; ib++) {
        uint32_t bit = startBit + ib;
        int bitVal = ((state & (1 << bit)) > 0);
        value += bitVal * std::round(std::pow(2., static_cast<double>(ib)));
      } // loop over bits for particle ip
      v[ip] = value;
    }
    return v;
  }
  // Return the probability associated with a given outcome
  double stateWeight(const std::vector<int>& state, const std::vector<std::vector<double>>& weights)
  {
    double w = 1.;
    for (uint32_t ip = 0; ip < state.size(); ip++) {
      w *= weights[ip][state[ip]];
    }
    return w;
  }
  // Return the corrected values for z and ptjet for a given state
  // Scale values by the fraction of jet momentum carried by removed V0s
  std::vector<double> correctedValues(const std::vector<int>& state, const std::vector<double>& values)
  {
    // Assumes values = (z1, z2, ..., zn, ptjet)
    std::vector<double> v(values);
    double r = 0;
    int nParticles = state.size();

    if (values.size() != static_cast<uint32_t>(nParticles + 1)) {
      LOGF(warning, "Number of values (%d) must be equal to the number of particles (%d) + 1!", values.size(), nParticles);
      return v;
    }
    for (int ip = 0; ip < nParticles; ip++) {
      if (state[ip] == 0) {
        r += values[ip];
      }
    }
    for (int ip = 0; ip < nParticles; ip++) {
      if (state[ip] == 0) {
        v[ip] = values[ip] / (1 - r);
      }
    }
    v[nParticles] = values[nParticles] * (1 - r);
    return v;
  }
  // Return the corrected values for z and ptjet for a given state
  // Take into account tracks that would have been included in the jet regardless of the V0s
  template <typename T, typename U, typename V>
  std::vector<double> correctedValuesPlusTracks(const std::vector<int>& state, V const& jet)
  {
    int ip = 0;
    double ptToSubtract = 0., ptToAdd = 0.;

    for (const auto& v0 : jet.template candidates_as<T>()) {
      if (v0.isRejectedCandidate())
        continue;

      // Background
      if (state[ip] == 0) {
        ptToSubtract += v0.pt();

        // TODO: This is okay in pt-scheme jets, but should add 4-momentum for E-scheme
        auto negTrack = v0.template negTrack_as<U>();
        if (jetderiveddatautilities::selectTrack(negTrack, trackSelection) && jetutilities::deltaR(jet, negTrack) < jet.r() * 1e-2)
          ptToAdd += negTrack.pt();

        auto posTrack = v0.template posTrack_as<U>();
        if (jetderiveddatautilities::selectTrack(posTrack, trackSelection) && jetutilities::deltaR(jet, negTrack) < jet.r() * 1e-2)
          ptToAdd += posTrack.pt();
      }
      ip++;
    }

    double ptjet = jet.pt() - ptToSubtract + ptToAdd;
    std::vector<double> values;

    ip = 0;
    for (const auto& v0 : jet.template candidates_as<T>()) {
      if (v0.isRejectedCandidate())
        continue;

      if (state[ip] == 0) // Background
        values.push_back(v0.pt() / jet.pt());
      else
        values.push_back(v0.pt() / ptjet);
      ip++;
    }
    values.push_back(ptjet);
    return values;
  }

  // ---------------------------------------------------
  // Helper functions
  // ---------------------------------------------------
  template <typename JetType>
  bool jetContainsV0s(JetType const& jet)
  {
    return (jet.candidatesIds().size() > 0);
  }
  template <typename T, typename U, typename V>
  // bool v0sAreMatched(T const& v0, U const& particle, V const& /*tracks*/)
  bool v0sAreMatched(U const& v0, V const& particle)
  {
    auto negId = v0.template negTrack_as<T>().mcParticleId();
    auto posId = v0.template posTrack_as<T>().mcParticleId();
    auto daughters = particle.daughtersIds();
    return ((negId == daughters[0] && posId == daughters[1]) || (posId == daughters[0] && negId == daughters[1]));
  }
  template <typename V0Type>
  double getReflectedMass(V0Type const& v0, bool isLambda)
  {
    // If V0 is Lambda, posTrack = proton, negTrack = pion
    // In that case, we assign pion mass to posTrack and proton mass to negTrack to calculate the reflection
    // Vice versa for AntiLambda
    float negM = (isLambda ? constants::physics::MassProton : constants::physics::MassPionCharged);
    float posM = (isLambda ? constants::physics::MassPionCharged : constants::physics::MassProton);
    std::array<std::array<float, 3>, 2> momenta = {std::array{v0.pxpos(), v0.pypos(), v0.pzpos()}, std::array{v0.pxneg(), v0.pyneg(), v0.pzneg()}};
    std::array<float, 2> masses = {posM, negM};
    return RecoDecay::m(momenta, masses);
  }
  template <typename Jet, typename Constituent>
  double getMomFrac(Jet const& jet, Constituent const& constituent)
  {
    double divByZeroProtect = 1e-5;
    if (jet.pt() < divByZeroProtect)
      return -1.;
    else
      return constituent.pt() / jet.pt();
  }
  template <typename Jet, typename Constituent>
  double getMomProj(Jet const& jet, Constituent const& constituent)
  {
    double divByZeroProtect = 1e-5;
    if (jet.p() < divByZeroProtect)
      return -1.;

    double trackProj = constituent.px() * jet.px() + constituent.py() * jet.py() + constituent.pz() * jet.pz();
    trackProj /= (jet.p() * jet.p());
    return trackProj;
  }

  // ---------------------------------------------------
  // Histogram filling functions
  // ---------------------------------------------------
  // Data - Counts
  template <typename CollisionType, typename V0Type>
  void fillDataV0sInclusive(CollisionType const& coll, V0Type const& V0s)
  {
    // Fill histograms unweighted. Hists will be filled with V0 counts
    float nV0s = 0;
    for (const auto& v0 : V0s) {
      if (v0.isRejectedCandidate())
        continue;

      nV0s += 1;
      double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
      double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
      double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

      double massDiff = v0.mLambda() - v0.mAntiLambda();
      double massRatio = v0.mAntiLambda() / v0.mLambda();
      double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

      registry.fill(HIST("data/V0/V0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
      registry.fill(HIST("data/V0/V0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda);
      registry.fill(HIST("data/V0/V0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/V0/V0PtMassWide"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/V0/V0PtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/V0/V0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
      registry.fill(HIST("data/V0/V0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/V0/V0PtDCAd"), v0.pt(), v0.dcaV0daughters());

      registry.fill(HIST("data/V0/V0CutVariation"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), v0.v0radius(), ctauK0s, v0.v0cosPA(), std::abs(v0.dcapostopv()), std::abs(v0.dcanegtopv()), v0.dcaV0daughters());

      if (v0.isK0SCandidate()) {
        registry.fill(HIST("data/V0/K0SPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/V0/K0SPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/V0/K0SPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/V0/K0SPtDCAd"), v0.pt(), v0.dcaV0daughters());

        registry.fill(HIST("data/V0/K0SPtCtauMass"), v0.pt(), ctauK0s, v0.mK0Short());
        registry.fill(HIST("data/V0/K0SPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mK0Short());
        registry.fill(HIST("data/V0/K0SPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mK0Short());
        registry.fill(HIST("data/V0/K0SPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mK0Short());
        registry.fill(HIST("data/V0/K0SPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mK0Short());
        registry.fill(HIST("data/V0/K0SPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mK0Short());
      }
      if (v0.isLambdaCandidate()) {
        registry.fill(HIST("data/V0/LambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/V0/LambdaPtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff);
        registry.fill(HIST("data/V0/LambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/V0/LambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/V0/LambdaPtDCAd"), v0.pt(), v0.dcaV0daughters());

        registry.fill(HIST("data/V0/LambdaPtCtauMass"), v0.pt(), ctauLambda, v0.mLambda());
        registry.fill(HIST("data/V0/LambdaPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mLambda());
        registry.fill(HIST("data/V0/LambdaPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mLambda());
        registry.fill(HIST("data/V0/LambdaPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mLambda());
        registry.fill(HIST("data/V0/LambdaPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mLambda());
        registry.fill(HIST("data/V0/LambdaPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mLambda());
      }
      if (v0.isAntiLambdaCandidate()) {
        registry.fill(HIST("data/V0/AntiLambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/V0/AntiLambdaPtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff);
        registry.fill(HIST("data/V0/AntiLambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/V0/AntiLambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/V0/AntiLambdaPtDCAd"), v0.pt(), v0.dcaV0daughters());

        registry.fill(HIST("data/V0/AntiLambdaPtCtauMass"), v0.pt(), ctauAntiLambda, v0.mAntiLambda());
        registry.fill(HIST("data/V0/AntiLambdaPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mAntiLambda());
        registry.fill(HIST("data/V0/AntiLambdaPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mAntiLambda());
        registry.fill(HIST("data/V0/AntiLambdaPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mAntiLambda());
        registry.fill(HIST("data/V0/AntiLambdaPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mAntiLambda());
        registry.fill(HIST("data/V0/AntiLambdaPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mAntiLambda());
      }
    } // for v0
    registry.fill(HIST("data/V0/nV0sEventAcc"), nV0s);
  }
  template <typename T>
  void fillDataJet(T const& jet)
  {
    registry.fill(HIST("data/jets/inclJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
    if (jetContainsV0s(jet))
      registry.fill(HIST("data/jets/jetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
  }
  void fillDataJetWeighted(double jetpt, double jeteta, double jetphi, double weight = 1.)
  {
    registry.fill(HIST("data/jets/weighted/jetPtEtaPhi"), jetpt, jeteta, jetphi, weight);
  }
  template <typename CollisionType, typename JetType, typename V0Type>
  void fillDataV0sInJet(CollisionType const& coll, JetType const& jet, V0Type const& v0)
  {
    double trackProj = getMomFrac(jet, v0);
    double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

    double massDiff = v0.mLambda() - v0.mAntiLambda();
    double massRatio = v0.mAntiLambda() / v0.mLambda();
    double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

    registry.fill(HIST("data/jets/V0/jetPtV0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi());
    registry.fill(HIST("data/jets/V0/jetPtV0PtCtau"), jet.pt(), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda);
    registry.fill(HIST("data/jets/V0/jetPtV0PtMass"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
    registry.fill(HIST("data/jets/V0/jetPtV0PtMassWide"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
    registry.fill(HIST("data/jets/V0/jetPtV0PtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff);
    registry.fill(HIST("data/jets/V0/jetPtV0PtRadiusCosPA"), jet.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA());
    registry.fill(HIST("data/jets/V0/jetPtV0PtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
    registry.fill(HIST("data/jets/V0/jetPtV0PtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());

    registry.fill(HIST("data/jets/V0/jetPtV0TrackProj"), jet.pt(), trackProj);
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjCtau"), jet.pt(), trackProj, ctauK0s, ctauLambda, ctauAntiLambda);
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjMass"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjMassWide"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff);
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjRadiusCosPA"), jet.pt(), trackProj, v0.v0radius(), v0.v0cosPA());
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());

    if (v0.isK0SCandidate()) {
      registry.fill(HIST("data/jets/V0/jetPtK0SPtCtau"), jet.pt(), v0.pt(), ctauK0s);
      registry.fill(HIST("data/jets/V0/jetPtK0SPtMass"), jet.pt(), v0.pt(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtAllMasses"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtCtauMass"), jet.pt(), v0.pt(), ctauK0s, v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtRadiusMass"), jet.pt(), v0.pt(), v0.v0radius(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtCosPAMass"), jet.pt(), v0.pt(), v0.v0cosPA(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAposMass"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAnegMass"), jet.pt(), v0.pt(), v0.dcanegtopv(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAdMass"), jet.pt(), v0.pt(), v0.dcaV0daughters(), v0.mK0Short());

      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjCtau"), jet.pt(), trackProj, ctauK0s);
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjMass"), jet.pt(), trackProj, v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjAllMasses"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjRadius"), jet.pt(), trackProj, v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjCosPA"), jet.pt(), trackProj, v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjCtauMass"), jet.pt(), trackProj, ctauK0s, v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjRadiusMass"), jet.pt(), trackProj, v0.v0radius(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjCosPAMass"), jet.pt(), trackProj, v0.v0cosPA(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAposMass"), jet.pt(), trackProj, v0.dcapostopv(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAnegMass"), jet.pt(), trackProj, v0.dcanegtopv(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAdMass"), jet.pt(), trackProj, v0.dcaV0daughters(), v0.mK0Short());
    }
    if (v0.isLambdaCandidate()) {
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtCtau"), jet.pt(), v0.pt(), ctauLambda);
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtMass"), jet.pt(), v0.pt(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtAllMasses"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtCtauMass"), jet.pt(), v0.pt(), ctauLambda, v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtRadiusMass"), jet.pt(), v0.pt(), v0.v0radius(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtCosPAMass"), jet.pt(), v0.pt(), v0.v0cosPA(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAposMass"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAnegMass"), jet.pt(), v0.pt(), v0.dcanegtopv(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAdMass"), jet.pt(), v0.pt(), v0.dcaV0daughters(), v0.mLambda());

      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjCtau"), jet.pt(), trackProj, ctauLambda);
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjMass"), jet.pt(), trackProj, v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjAllMasses"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjRadius"), jet.pt(), trackProj, v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjCosPA"), jet.pt(), trackProj, v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjCtauMass"), jet.pt(), trackProj, ctauLambda, v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjRadiusMass"), jet.pt(), trackProj, v0.v0radius(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjCosPAMass"), jet.pt(), trackProj, v0.v0cosPA(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAposMass"), jet.pt(), trackProj, v0.dcapostopv(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAnegMass"), jet.pt(), trackProj, v0.dcanegtopv(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAdMass"), jet.pt(), trackProj, v0.dcaV0daughters(), v0.mLambda());
    }
    if (v0.isAntiLambdaCandidate()) {
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtCtau"), jet.pt(), v0.pt(), ctauAntiLambda);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtMass"), jet.pt(), v0.pt(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtAllMasses"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtCtauMass"), jet.pt(), v0.pt(), ctauAntiLambda, v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtRadiusMass"), jet.pt(), v0.pt(), v0.v0radius(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtCosPAMass"), jet.pt(), v0.pt(), v0.v0cosPA(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAposMass"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAnegMass"), jet.pt(), v0.pt(), v0.dcanegtopv(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAdMass"), jet.pt(), v0.pt(), v0.dcaV0daughters(), v0.mAntiLambda());

      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjCtau"), jet.pt(), trackProj, ctauAntiLambda);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjMass"), jet.pt(), trackProj, v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjAllMasses"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjRadius"), jet.pt(), trackProj, v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjCosPA"), jet.pt(), trackProj, v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjCtauMass"), jet.pt(), trackProj, ctauAntiLambda, v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjRadiusMass"), jet.pt(), trackProj, v0.v0radius(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjCosPAMass"), jet.pt(), trackProj, v0.v0cosPA(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjDCAposMass"), jet.pt(), trackProj, v0.dcapostopv(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjDCAnegMass"), jet.pt(), trackProj, v0.dcanegtopv(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjDCAdMass"), jet.pt(), trackProj, v0.dcaV0daughters(), v0.mAntiLambda());
    }
  }
  template <typename T, typename U, typename V>
  void fillDataV0sInPerpCone(T const& coll, U const& jet, V const& v0s)
  {
    const int nCones = 2;
    double perpConeR = jet.r() * 1e-2;
    double conePhi[nCones] = {RecoDecay::constrainAngle(jet.phi() - constants::math::PIHalf, -constants::math::PI),
                              RecoDecay::constrainAngle(jet.phi() + constants::math::PIHalf, -constants::math::PI)};
    double conePt[nCones] = {0., 0.};
    int nV0sinCone[nCones] = {0, 0};
    for (const auto& v0 : v0s) {
      if (v0.isRejectedCandidate())
        continue;

      bool v0InCones = false;
      double dEta = v0.eta() - jet.eta();
      double dPhi[nCones] = {RecoDecay::constrainAngle(v0.phi() - conePhi[0], -constants::math::PI),
                             RecoDecay::constrainAngle(v0.phi() - conePhi[1], -constants::math::PI)};
      for (int i = 0; i < nCones; i++) {
        if (std::sqrt(dEta * dEta + dPhi[i] * dPhi[i]) < perpConeR) {
          conePt[i] += v0.pt();
          nV0sinCone[i]++;
          v0InCones = true;
        }
      }
      if (!v0InCones) {
        continue;
      }

      double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
      double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
      double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

      registry.fill(HIST("data/PC/JetPtEtaV0Pt"), jet.pt(), jet.eta(), v0.pt());
      registry.fill(HIST("data/PC/V0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
      registry.fill(HIST("data/PC/V0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda);
      registry.fill(HIST("data/PC/V0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/PC/V0PtMassWide"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/PC/V0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
      registry.fill(HIST("data/PC/V0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/PC/V0PtDCAd"), v0.pt(), v0.dcaV0daughters());

      if (v0.isLambdaCandidate()) {
        registry.fill(HIST("data/PC/JetPtLambdaPtMass"), jet.pt(), v0.pt(), v0.mLambda());

        registry.fill(HIST("data/PC/JetPtEtaLambdaPt"), jet.pt(), jet.eta(), v0.pt());
        registry.fill(HIST("data/PC/LambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/PC/LambdaPtCtauMass"), v0.pt(), ctauLambda, v0.mLambda());
        registry.fill(HIST("data/PC/LambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/PC/LambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/PC/LambdaPtDCAd"), v0.pt(), v0.dcaV0daughters());
      }
      if (v0.isAntiLambdaCandidate()) {
        registry.fill(HIST("data/PC/JetPtAntiLambdaPtMass"), jet.pt(), v0.pt(), v0.mAntiLambda());

        registry.fill(HIST("data/PC/JetPtEtaAntiLambdaPt"), jet.pt(), jet.eta(), v0.pt());
        registry.fill(HIST("data/PC/AntiLambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/PC/AntiLambdaPtCtauMass"), v0.pt(), ctauAntiLambda, v0.mAntiLambda());
        registry.fill(HIST("data/PC/AntiLambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/PC/AntiLambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/PC/AntiLambdaPtDCAd"), v0.pt(), v0.dcaV0daughters());
      }
      if (v0.isK0SCandidate()) {
        registry.fill(HIST("data/PC/JetPtK0SPtMass"), jet.pt(), v0.pt(), v0.mK0Short());

        registry.fill(HIST("data/PC/JetPtEtaK0SPt"), jet.pt(), jet.eta(), v0.pt());
        registry.fill(HIST("data/PC/K0SPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/PC/K0SPtCtauMass"), v0.pt(), ctauK0s, v0.mK0Short());
        registry.fill(HIST("data/PC/K0SPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/PC/K0SPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/PC/K0SPtDCAd"), v0.pt(), v0.dcaV0daughters());
      }
    }
    // Fill hist for nCones: nv0s, conePt, coneEta, conePhi
    for (int i = 0; i < nCones; i++) {
      registry.fill(HIST("data/PC/nV0sConePtEta"), nV0sinCone[i], conePt[i], jet.eta());
      registry.fill(HIST("data/PC/ConePtEtaPhi"), conePt[i], jet.eta(), conePhi[i]);
      registry.fill(HIST("data/PC/JetPtEtaConePt"), jet.pt(), jet.eta(), conePt[i]);
    }
  }
  // Data - Weighted
  template <typename CollisionType, typename V0Type>
  void fillDataV0sInclusiveWeighted(CollisionType const& coll, V0Type const& V0s)
  {
    // Fill histograms with V0 signal weights
    float nV0s = 0;
    for (const auto& v0 : V0s) {
      if (v0.isRejectedCandidate())
        continue;

      float weight = getV0SignalProb(v0);
      nV0s += weight; // Sum weights (purity) of V0s

      double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
      double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
      double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

      double massDiff = v0.mLambda() - v0.mAntiLambda();
      double massRatio = v0.mAntiLambda() / v0.mLambda();
      double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

      registry.fill(HIST("data/V0/V0PtEtaPhiWeighted"), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("data/V0/V0PtCtauWeighted"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("data/V0/V0PtMassWeighted"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("data/V0/V0PtMassWideWeighted"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("data/V0/V0PtLambdaMassesWeighted"), v0.pt(), massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("data/V0/V0PtRadiusCosPAWeighted"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("data/V0/V0PtDCAposnegWeighted"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("data/V0/V0PtDCAdWeighted"), v0.pt(), v0.dcaV0daughters(), weight);

      if (v0.isK0SCandidate()) {
        registry.fill(HIST("data/V0/K0SPtEtaPhiWeighted"), v0.pt(), v0.eta(), v0.phi(), weight);
        registry.fill(HIST("data/V0/K0SPtRadiusCosPAWeighted"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
        registry.fill(HIST("data/V0/K0SPtDCAposnegWeighted"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
        registry.fill(HIST("data/V0/K0SPtDCAdWeighted"), v0.pt(), v0.dcaV0daughters(), weight);

        registry.fill(HIST("data/V0/K0SPtCtauMassWeighted"), v0.pt(), ctauK0s, v0.mK0Short(), weight);
        registry.fill(HIST("data/V0/K0SPtRadiusMassWeighted"), v0.pt(), v0.v0radius(), v0.mK0Short(), weight);
        registry.fill(HIST("data/V0/K0SPtCosPAMassWeighted"), v0.pt(), v0.v0cosPA(), v0.mK0Short(), weight);
        registry.fill(HIST("data/V0/K0SPtDCAposMassWeighted"), v0.pt(), v0.dcapostopv(), v0.mK0Short(), weight);
        registry.fill(HIST("data/V0/K0SPtDCAnegMassWeighted"), v0.pt(), v0.dcanegtopv(), v0.mK0Short(), weight);
        registry.fill(HIST("data/V0/K0SPtDCAdMassWeighted"), v0.pt(), v0.dcaV0daughters(), v0.mK0Short(), weight);
      }
      if (v0.isLambdaCandidate()) {
        registry.fill(HIST("data/V0/LambdaPtEtaPhiWeighted"), v0.pt(), v0.eta(), v0.phi(), weight);
        registry.fill(HIST("data/V0/LambdaPtLambdaMassesWeighted"), v0.pt(), massDiff, massRatio, massRelDiff, weight);
        registry.fill(HIST("data/V0/LambdaPtRadiusCosPAWeighted"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
        registry.fill(HIST("data/V0/LambdaPtDCAposnegWeighted"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
        registry.fill(HIST("data/V0/LambdaPtDCAdWeighted"), v0.pt(), v0.dcaV0daughters(), weight);

        registry.fill(HIST("data/V0/LambdaPtCtauMassWeighted"), v0.pt(), ctauLambda, v0.mLambda(), weight);
        registry.fill(HIST("data/V0/LambdaPtRadiusMassWeighted"), v0.pt(), v0.v0radius(), v0.mLambda(), weight);
        registry.fill(HIST("data/V0/LambdaPtCosPAMassWeighted"), v0.pt(), v0.v0cosPA(), v0.mLambda(), weight);
        registry.fill(HIST("data/V0/LambdaPtDCAposMassWeighted"), v0.pt(), v0.dcapostopv(), v0.mLambda(), weight);
        registry.fill(HIST("data/V0/LambdaPtDCAnegMassWeighted"), v0.pt(), v0.dcanegtopv(), v0.mLambda(), weight);
        registry.fill(HIST("data/V0/LambdaPtDCAdMassWeighted"), v0.pt(), v0.dcaV0daughters(), v0.mLambda(), weight);
      }
      if (v0.isAntiLambdaCandidate()) {
        registry.fill(HIST("data/V0/AntiLambdaPtEtaPhiWeighted"), v0.pt(), v0.eta(), v0.phi(), weight);
        registry.fill(HIST("data/V0/AntiLambdaPtLambdaMassesWeighted"), v0.pt(), massDiff, massRatio, massRelDiff, weight);
        registry.fill(HIST("data/V0/AntiLambdaPtRadiusCosPAWeighted"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
        registry.fill(HIST("data/V0/AntiLambdaPtDCAposnegWeighted"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
        registry.fill(HIST("data/V0/AntiLambdaPtDCAdWeighted"), v0.pt(), v0.dcaV0daughters(), weight);

        registry.fill(HIST("data/V0/AntiLambdaPtCtauMassWeighted"), v0.pt(), ctauAntiLambda, v0.mAntiLambda(), weight);
        registry.fill(HIST("data/V0/AntiLambdaPtRadiusMassWeighted"), v0.pt(), v0.v0radius(), v0.mAntiLambda(), weight);
        registry.fill(HIST("data/V0/AntiLambdaPtCosPAMassWeighted"), v0.pt(), v0.v0cosPA(), v0.mAntiLambda(), weight);
        registry.fill(HIST("data/V0/AntiLambdaPtDCAposMassWeighted"), v0.pt(), v0.dcapostopv(), v0.mAntiLambda(), weight);
        registry.fill(HIST("data/V0/AntiLambdaPtDCAnegMassWeighted"), v0.pt(), v0.dcanegtopv(), v0.mAntiLambda(), weight);
        registry.fill(HIST("data/V0/AntiLambdaPtDCAdMassWeighted"), v0.pt(), v0.dcaV0daughters(), v0.mAntiLambda(), weight);
      }
    } // for v0
    registry.fill(HIST("data/V0/nV0sEventAccWeighted"), nV0s);
  }
  template <typename C, typename J>
  void fillDataV0sInJetWeighted(C const& coll, J const& jet, const std::vector<int>& state, const std::vector<double>& values, double weight)
  {
    double jetpt = values[values.size() - 1];
    int ip = 0;
    for (const auto& v0 : jet.template candidates_as<aod::CandidatesV0Data>()) {
      if (v0.isRejectedCandidate())
        continue;

      double z = values[ip];
      ip++;

      double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;
      double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
      double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;

      double massDiff = v0.mLambda() - v0.mAntiLambda();
      double massRatio = v0.mAntiLambda() / v0.mLambda();
      double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

      switch (state[ip]) {
        case 0: // Background
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjCtau"), jetpt, z, ctauK0s, ctauLambda, ctauAntiLambda, weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjMass"), jetpt, z, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjRadiusCosPA"), jetpt, z, v0.v0radius(), v0.v0cosPA(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAposneg"), jetpt, z, v0.dcapostopv(), v0.dcanegtopv(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAd"), jetpt, z, v0.dcaV0daughters(), weight);

          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjCtauK0SMass"), jetpt, z, ctauK0s, v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjRadiusK0SMass"), jetpt, z, v0.v0radius(), v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjCosPAK0SMass"), jetpt, z, v0.v0cosPA(), v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAposK0SMass"), jetpt, z, v0.dcapostopv(), v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAnegK0SMass"), jetpt, z, v0.dcanegtopv(), v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAdK0SMass"), jetpt, z, v0.dcaV0daughters(), v0.mK0Short(), weight);

          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjCtauLambdaMass"), jetpt, z, ctauLambda, v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjRadiusLambdaMass"), jetpt, z, v0.v0radius(), v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjCosPALambdaMass"), jetpt, z, v0.v0cosPA(), v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAposLambdaMass"), jetpt, z, v0.dcapostopv(), v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAnegLambdaMass"), jetpt, z, v0.dcanegtopv(), v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAdLambdaMass"), jetpt, z, v0.dcaV0daughters(), v0.mLambda(), weight);

          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjCtauAntiLambdaMass"), jetpt, z, ctauAntiLambda, v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjRadiusAntiLambdaMass"), jetpt, z, v0.v0radius(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjCosPAAntiLambdaMass"), jetpt, z, v0.v0cosPA(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAposAntiLambdaMass"), jetpt, z, v0.dcapostopv(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAnegAntiLambdaMass"), jetpt, z, v0.dcanegtopv(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtBkgTrackProjDCAdAntiLambdaMass"), jetpt, z, v0.dcaV0daughters(), v0.mAntiLambda(), weight);
          break;
        case 1: // K0S
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjMass"), jetpt, z, v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjCtau"), jetpt, z, ctauK0s, weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjAllMasses"), jetpt, z, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjRadius"), jetpt, z, v0.v0radius(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjCosPA"), jetpt, z, v0.v0cosPA(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjDCAd"), jetpt, z, v0.dcaV0daughters(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjDCAposneg"), jetpt, z, v0.dcapostopv(), v0.dcanegtopv(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjCtauMass"), jetpt, z, ctauK0s, v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjRadiusMass"), jetpt, z, v0.v0radius(), v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjCosPAMass"), jetpt, z, v0.v0cosPA(), v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjDCAposMass"), jetpt, z, v0.dcapostopv(), v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjDCAnegMass"), jetpt, z, v0.dcanegtopv(), v0.mK0Short(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtK0STrackProjDCAdMass"), jetpt, z, v0.dcaV0daughters(), v0.mK0Short(), weight);
          break;
        case 2: // Lambda
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjMass"), jetpt, z, v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjCtau"), jetpt, z, ctauLambda, weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjAllMasses"), jetpt, z, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjLambdaMasses"), jetpt, z, massDiff, massRatio, massRelDiff, weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjRadius"), jetpt, z, v0.v0radius(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjCosPA"), jetpt, z, v0.v0cosPA(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjDCAd"), jetpt, z, v0.dcaV0daughters(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjDCAposneg"), jetpt, z, v0.dcapostopv(), v0.dcanegtopv(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjCtauMass"), jetpt, z, ctauLambda, v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjRadiusMass"), jetpt, z, v0.v0radius(), v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjCosPAMass"), jetpt, z, v0.v0cosPA(), v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjDCAposMass"), jetpt, z, v0.dcapostopv(), v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjDCAnegMass"), jetpt, z, v0.dcanegtopv(), v0.mLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtLambdaTrackProjDCAdMass"), jetpt, z, v0.dcaV0daughters(), v0.mLambda(), weight);
          break;
        case 3: // AntiLambda
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjMass"), jetpt, z, v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjCtau"), jetpt, z, ctauAntiLambda, weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjAllMasses"), jetpt, z, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjLambdaMasses"), jetpt, z, massDiff, massRatio, massRelDiff, weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjRadius"), jetpt, z, v0.v0radius(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjCosPA"), jetpt, z, v0.v0cosPA(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAd"), jetpt, z, v0.dcaV0daughters(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAposneg"), jetpt, z, v0.dcapostopv(), v0.dcanegtopv(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjCtauMass"), jetpt, z, ctauAntiLambda, v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjRadiusMass"), jetpt, z, v0.v0radius(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjCosPAMass"), jetpt, z, v0.v0cosPA(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAposMass"), jetpt, z, v0.dcapostopv(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAnegMass"), jetpt, z, v0.dcanegtopv(), v0.mAntiLambda(), weight);
          registry.fill(HIST("data/jets/weighted/V0/jetPtAntiLambdaTrackProjDCAdMass"), jetpt, z, v0.dcaV0daughters(), v0.mAntiLambda(), weight);
          break;
      }
      registry.fill(HIST("data/jets/weighted/V0/jetPtV0TrackProjMass"), jetpt, z, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("data/jets/weighted/V0/jetPtV0TrackProjCtau"), jetpt, z, ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("data/jets/weighted/V0/jetPtV0TrackProjMassWide"), jetpt, z, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("data/jets/weighted/V0/jetPtV0TrackProjLambdaMasses"), jetpt, z, massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("data/jets/weighted/V0/jetPtV0TrackProjRadiusCosPA"), jetpt, z, v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("data/jets/weighted/V0/jetPtV0TrackProjDCAposneg"), jetpt, z, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("data/jets/weighted/V0/jetPtV0TrackProjDCAd"), jetpt, z, v0.dcaV0daughters(), weight);
    } // v0 loop
  }

  // MC - Counts (or event weights)
  template <typename T>
  void fillMcpV0sInclusive(T const& pV0s, double weight = 1.)
  {
    float nV0s = 0;
    for (const auto& pv0 : pV0s) {
      int pdg = pv0.pdgCode();
      nV0s += 1;
      registry.fill(HIST("mcp/V0/V0PtEtaPhi"), pv0.pt(), pv0.eta(), pv0.phi(), weight);
      if (std::abs(pdg) == PDG_t::kK0Short)
        registry.fill(HIST("mcp/V0/K0SPtEtaPhi"), pv0.pt(), pv0.eta(), pv0.phi(), weight);

      if (pdg == PDG_t::kLambda0)
        registry.fill(HIST("mcp/V0/LambdaPtEtaPhi"), pv0.pt(), pv0.eta(), pv0.phi(), weight);

      if (pdg == PDG_t::kLambda0Bar)
        registry.fill(HIST("mcp/V0/AntiLambdaPtEtaPhi"), pv0.pt(), pv0.eta(), pv0.phi(), weight);
    }
    registry.fill(HIST("mcp/V0/nV0sEventAcc"), nV0s);
    registry.fill(HIST("mcp/V0/nV0sEventAccWeighted"), nV0s, weight);
  }
  template <typename T>
  void fillMcpJet(T const& jet, double weight = 1.)
  {
    registry.fill(HIST("mcp/jets/inclPartJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    if (jetContainsV0s(jet)) {
      registry.fill(HIST("mcp/jets/partJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    }
  }
  template <typename T, typename U>
  void fillMcdV0sInclusive(T const& coll, U const& V0s, double weight = 1.)
  {
    float nV0s = 0;
    for (const auto& v0 : V0s) {
      if (v0.isRejectedCandidate())
        continue;

      nV0s += 1;
      double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
      double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
      double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

      double massDiff = v0.mLambda() - v0.mAntiLambda();
      double massRatio = v0.mAntiLambda() / v0.mLambda();
      double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

      registry.fill(HIST("mcd/V0/V0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("mcd/V0/V0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("mcd/V0/V0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("mcd/V0/V0PtLambdaMasses"), v0.pt(), v0.mLambda() - v0.mAntiLambda(), v0.mAntiLambda() / v0.mLambda(), (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda(), weight);
      registry.fill(HIST("mcd/V0/V0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("mcd/V0/V0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("mcd/V0/V0PtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

      if (v0.isK0SCandidate()) {
        registry.fill(HIST("mcd/V0/K0SPtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
        registry.fill(HIST("mcd/V0/K0SPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
        registry.fill(HIST("mcd/V0/K0SPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
        registry.fill(HIST("mcd/V0/K0SPtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

        registry.fill(HIST("mcd/V0/K0SPtCtauMass"), v0.pt(), ctauK0s, v0.mK0Short(), weight);
        registry.fill(HIST("mcd/V0/K0SPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mK0Short(), weight);
        registry.fill(HIST("mcd/V0/K0SPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mK0Short(), weight);
        registry.fill(HIST("mcd/V0/K0SPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mK0Short(), weight);
        registry.fill(HIST("mcd/V0/K0SPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mK0Short(), weight);
        registry.fill(HIST("mcd/V0/K0SPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mK0Short(), weight);
      }
      if (v0.isLambdaCandidate()) {
        registry.fill(HIST("mcd/V0/LambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
        registry.fill(HIST("mcd/V0/LambdaPtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff, weight);
        registry.fill(HIST("mcd/V0/LambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
        registry.fill(HIST("mcd/V0/LambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
        registry.fill(HIST("mcd/V0/LambdaPtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

        registry.fill(HIST("mcd/V0/LambdaPtCtauMass"), v0.pt(), ctauLambda, v0.mLambda(), weight);
        registry.fill(HIST("mcd/V0/LambdaPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mLambda(), weight);
        registry.fill(HIST("mcd/V0/LambdaPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mLambda(), weight);
        registry.fill(HIST("mcd/V0/LambdaPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mLambda(), weight);
        registry.fill(HIST("mcd/V0/LambdaPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mLambda(), weight);
        registry.fill(HIST("mcd/V0/LambdaPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mLambda(), weight);
      }
      if (v0.isAntiLambdaCandidate()) {
        registry.fill(HIST("mcd/V0/AntiLambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff, weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

        registry.fill(HIST("mcd/V0/AntiLambdaPtCtauMass"), v0.pt(), ctauAntiLambda, v0.mAntiLambda(), weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mAntiLambda(), weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mAntiLambda(), weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mAntiLambda(), weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mAntiLambda(), weight);
        registry.fill(HIST("mcd/V0/AntiLambdaPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mAntiLambda(), weight);
      }
    } // for v0
    registry.fill(HIST("mcd/V0/nV0sEventAcc"), nV0s);
    registry.fill(HIST("mcd/V0/nV0sEventAccWeighted"), nV0s, weight);
  }
  template <typename T, typename U, typename V, typename W, typename X>
  void fillMatchingV0sInclusive(V const& coll, W const& V0s, X const& pV0s, double weight = 1.)
  {
    float nV0s = 0;
    for (const auto& v0 : V0s) {
      if (v0.isRejectedCandidate())
        continue;

      if (!v0.has_mcParticle()) {
        fillFakeV0Inclusive(coll, v0, weight);
        fillFakeV0DaughtersInclusive<T>(v0, weight);
        fillFakeV0DecayedInclusive<T, U>(v0, weight);
        continue;
      }
      for (const auto& pv0 : pV0s) {
        if (v0sAreMatched<T>(v0, pv0)) {
          nV0s += 1;
          fillMatchedV0Inclusive(coll, v0, pv0, weight);
          fillMatchedV0DaughtersInclusive<T, U>(v0, pv0, weight);
        }
      }
    } // Reconstructed V0s
    for (const auto& pv0 : pV0s) {
      for (const auto& v0 : V0s) {
        if (v0sAreMatched<T>(v0, pv0))
          continue;

        fillMissV0Inclusive(pv0);
      }
    }
    registry.fill(HIST("matching/V0/nV0sEventAcc"), nV0s);
    registry.fill(HIST("matching/V0/nV0sEventAccWeighted"), nV0s, weight);
  }
  template <typename T>
  void fillMcdJet(T const& jet, double weight = 1.)
  {
    registry.fill(HIST("mcd/jets/inclDetJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    if (jetContainsV0s(jet)) {
      registry.fill(HIST("mcd/jets/detJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    }
  }
  // Reconstructed V0s in the cone of MCD jets
  template <typename T, typename U, typename V, typename W>
  void fillMcV0sInPerpCone(T const& coll, U const& mcdjet, V const& v0s, W const& /* V0 particles */, double weight = 1.)
  {
    const int nCones = 2;
    double perpConeR = mcdjet.r() * 1e-2;
    double conePhi[nCones] = {RecoDecay::constrainAngle(mcdjet.phi() - constants::math::PIHalf, -constants::math::PI),
                              RecoDecay::constrainAngle(mcdjet.phi() + constants::math::PIHalf, -constants::math::PI)};
    double coneMatchedPt[nCones] = {0., 0.};
    double coneFakePt[nCones] = {0., 0.};
    int nMatchedV0sinCone[nCones] = {0, 0};
    int nFakeV0sinCone[nCones] = {0, 0};

    for (const auto& v0 : v0s) {
      double dEta = v0.eta() - mcdjet.eta();
      double dPhi[nCones] = {RecoDecay::constrainAngle(v0.phi() - conePhi[0], -constants::math::PI),
                             RecoDecay::constrainAngle(v0.phi() - conePhi[1], -constants::math::PI)};
      for (int i = 0; i < nCones; i++) {
        if (std::sqrt(dEta * dEta + dPhi[i] * dPhi[i]) > perpConeR) {
          continue;
        }

        double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
        double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
        double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

        if (!v0.has_mcParticle()) { // The V0 is combinatorial background
          coneFakePt[i] += v0.pt();
          nFakeV0sinCone[i]++;
          registry.fill(HIST("mcd/PC/jetPtEtaFakeV0Pt"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);

          registry.fill(HIST("mcd/PC/fakeV0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
          registry.fill(HIST("mcd/PC/fakeV0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
          registry.fill(HIST("mcd/PC/fakeV0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          registry.fill(HIST("mcd/PC/fakeV0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
          registry.fill(HIST("mcd/PC/fakeV0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
          registry.fill(HIST("mcd/PC/fakeV0PtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);
        } else {
          coneMatchedPt[i] += v0.pt();
          nMatchedV0sinCone[i]++;
          registry.fill(HIST("mcd/PC/jetPtEtaMatchedV0Pt"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);

          registry.fill(HIST("mcd/PC/matchedV0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
          registry.fill(HIST("mcd/PC/matchedV0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
          registry.fill(HIST("mcd/PC/matchedV0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          registry.fill(HIST("mcd/PC/matchedV0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
          registry.fill(HIST("mcd/PC/matchedV0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
          registry.fill(HIST("mcd/PC/matchedV0PtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

          auto particle = v0.template mcParticle_as<W>();
          if (std::abs(particle.pdgCode()) == PDG_t::kK0Short) { // K0S
            registry.fill(HIST("mcd/PC/matchedJetPtK0SPtMass"), mcdjet.pt(), v0.pt(), v0.mK0Short(), weight);
          } else if (particle.pdgCode() == PDG_t::kLambda0) { // Lambda
            registry.fill(HIST("mcd/PC/matchedJetPtLambdaPtMass"), mcdjet.pt(), v0.pt(), v0.mLambda(), weight);
          } else if (particle.pdgCode() == PDG_t::kLambda0Bar) {
            registry.fill(HIST("mcd/PC/matchedJetPtAntiLambdaPtMass"), mcdjet.pt(), v0.pt(), v0.mAntiLambda(), weight);
          }
        } // if v0 has mcParticle
      } // for cone
    } // for v0s
    for (int i = 0; i < nCones; i++) {
      registry.fill(HIST("mcd/PC/matchednV0sConePtEta"), nMatchedV0sinCone[i], coneMatchedPt[i], mcdjet.eta(), weight);
      registry.fill(HIST("mcd/PC/matchedConePtEtaPhi"), coneMatchedPt[i], mcdjet.eta(), conePhi[i], weight);
      registry.fill(HIST("mcd/PC/matchedJetPtEtaConePt"), mcdjet.pt(), mcdjet.eta(), coneMatchedPt[i], weight);

      registry.fill(HIST("mcd/PC/fakenV0sConePtEta"), nFakeV0sinCone[i], coneFakePt[i], mcdjet.eta(), weight);
      registry.fill(HIST("mcd/PC/fakeConePtEtaPhi"), coneFakePt[i], mcdjet.eta(), conePhi[i], weight);
      registry.fill(HIST("mcd/PC/fakeJetPtEtaConePt"), mcdjet.pt(), mcdjet.eta(), coneFakePt[i], weight);
    }
  }
  // Reconstructed V0s in the cone of matched jets
  template <typename T, typename U, typename V, typename W, typename X>
  void fillMcV0sInMatchedPerpCone(T const& coll, U const& mcdjet, V const& mcpjet, W const& v0s, X const& /* V0 particles */, double weight = 1.)
  {
    const int nCones = 2;
    double perpConeR = mcdjet.r() * 1e-2;
    double conePhi[nCones] = {RecoDecay::constrainAngle(mcdjet.phi() - constants::math::PIHalf, -constants::math::PI),
                              RecoDecay::constrainAngle(mcdjet.phi() + constants::math::PIHalf, -constants::math::PI)};
    double coneMatchedPt[nCones] = {0., 0.};
    double coneFakePt[nCones] = {0., 0.};
    int nMatchedV0sinCone[nCones] = {0, 0};
    int nFakeV0sinCone[nCones] = {0, 0};

    for (const auto& v0 : v0s) {
      double dEta = v0.eta() - mcdjet.eta();
      double dPhi[nCones] = {RecoDecay::constrainAngle(v0.phi() - conePhi[0], -constants::math::PI),
                             RecoDecay::constrainAngle(v0.phi() - conePhi[1], -constants::math::PI)};
      for (int i = 0; i < nCones; i++) {
        if (std::sqrt(dEta * dEta + dPhi[i] * dPhi[i]) > perpConeR) {
          continue;
        }

        double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
        double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
        double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

        if (!v0.has_mcParticle()) { // The V0 is combinatorial background
          coneFakePt[i] += v0.pt();
          nFakeV0sinCone[i]++;
          registry.fill(HIST("matching/PC/jetPtEtaFakeV0Pt"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
          registry.fill(HIST("matching/PC/jetsPtFakeV0Pt"), mcpjet.pt(), mcdjet.pt(), v0.pt(), weight);

          registry.fill(HIST("matching/PC/fakeV0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
          registry.fill(HIST("matching/PC/fakeV0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
          registry.fill(HIST("matching/PC/fakeV0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          registry.fill(HIST("matching/PC/fakeV0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
          registry.fill(HIST("matching/PC/fakeV0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
          registry.fill(HIST("matching/PC/fakeV0PtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);
        } else {
          coneMatchedPt[i] += v0.pt();
          nMatchedV0sinCone[i]++;
          registry.fill(HIST("matching/PC/jetPtEtaMatchedV0Pt"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
          registry.fill(HIST("matching/PC/jetsPtMatchedV0Pt"), mcpjet.pt(), mcdjet.pt(), v0.pt(), weight);

          registry.fill(HIST("matching/PC/matchedV0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
          registry.fill(HIST("matching/PC/matchedV0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
          registry.fill(HIST("matching/PC/matchedV0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          registry.fill(HIST("matching/PC/matchedV0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
          registry.fill(HIST("matching/PC/matchedV0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
          registry.fill(HIST("matching/PC/matchedV0PtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

          auto particle = v0.template mcParticle_as<X>();
          if (std::abs(particle.pdgCode()) == PDG_t::kK0Short) { // K0S
            registry.fill(HIST("matching/PC/matchedJetPtK0SPtMass"), mcdjet.pt(), v0.pt(), v0.mK0Short(), weight);
            registry.fill(HIST("matching/PC/matchedJetsPtK0SPtMass"), mcpjet.pt(), mcdjet.pt(), v0.pt(), v0.mK0Short(), weight);
          } else if (particle.pdgCode() == PDG_t::kLambda0) { // Lambda
            registry.fill(HIST("matching/PC/matchedJetPtLambdaPtMass"), mcdjet.pt(), v0.pt(), v0.mLambda(), weight);
            registry.fill(HIST("matching/PC/matchedJetsPtLambdaPtMass"), mcpjet.pt(), mcdjet.pt(), v0.pt(), v0.mLambda(), weight);
          } else if (particle.pdgCode() == PDG_t::kLambda0Bar) {
            registry.fill(HIST("matching/PC/matchedJetPtAntiLambdaPtMass"), mcdjet.pt(), v0.pt(), v0.mAntiLambda(), weight);
            registry.fill(HIST("matching/PC/matchedJetsPtAntiLambdaPtMass"), mcpjet.pt(), mcdjet.pt(), v0.pt(), v0.mAntiLambda(), weight);
          }
        } // if v0 has mcParticle
      } // for cone
    } // for v0s
    for (int i = 0; i < nCones; i++) {
      registry.fill(HIST("matching/PC/matchednV0sConePtEta"), nMatchedV0sinCone[i], coneMatchedPt[i], mcdjet.eta(), weight);
      registry.fill(HIST("matching/PC/matchedConePtEtaPhi"), coneMatchedPt[i], mcdjet.eta(), conePhi[i], weight);
      registry.fill(HIST("matching/PC/matchedJetPtEtaConePt"), mcdjet.pt(), mcdjet.eta(), coneMatchedPt[i], weight);
      registry.fill(HIST("matching/PC/matchedJetsPtEtaConePt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), coneMatchedPt[i], weight);

      registry.fill(HIST("matching/PC/fakenV0sConePtEta"), nFakeV0sinCone[i], coneFakePt[i], mcdjet.eta(), weight);
      registry.fill(HIST("matching/PC/fakeConePtEtaPhi"), coneFakePt[i], mcdjet.eta(), conePhi[i], weight);
      registry.fill(HIST("matching/PC/fakeJetPtEtaConePt"), mcdjet.pt(), mcdjet.eta(), coneFakePt[i], weight);
      registry.fill(HIST("matching/PC/fakeJetsPtEtaConePt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), coneFakePt[i], weight);
    }
  }
  // Matched - Counts (or event weights)
  template <typename CollisionType, typename V0Type, typename particleType> // Reconstructed signal for inclusive V0s
  void fillMatchedV0Inclusive(CollisionType const& coll, V0Type const& v0, particleType const& particle, double weight = 1.)
  {
    double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

    registry.fill(HIST("matching/V0/V0PartPtDetPt"), particle.pt(), v0.pt(), weight);
    registry.fill(HIST("matching/V0/V0PartPtRatioPtRelDiffPt"), particle.pt(), v0.pt() / particle.pt(), (v0.pt() - particle.pt()) / particle.pt(), weight);

    if (std::abs(particle.pdgCode()) == PDG_t::kK0Short) { // K0S
      registry.fill(HIST("matching/V0/K0SPtEtaPhi"), particle.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/V0/K0SPtCtauMass"), particle.pt(), v0.pt(), ctauK0s, v0.mK0Short(), weight);
      registry.fill(HIST("matching/V0/K0SPtRadiusCosPA"), particle.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/V0/K0SPtDCAposneg"), particle.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/V0/K0SPtDCAd"), particle.pt(), v0.pt(), v0.dcaV0daughters(), weight);
      registry.fill(HIST("matching/V0/K0SPtMass"), particle.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
    } else if (particle.pdgCode() == PDG_t::kLambda0) { // Lambda
      registry.fill(HIST("matching/V0/LambdaPtEtaPhi"), particle.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/V0/LambdaPtCtauMass"), particle.pt(), v0.pt(), ctauLambda, v0.mLambda(), weight);
      registry.fill(HIST("matching/V0/LambdaPtRadiusCosPA"), particle.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/V0/LambdaPtDCAposneg"), particle.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/V0/LambdaPtDCAd"), particle.pt(), v0.pt(), v0.dcaV0daughters(), weight);
      registry.fill(HIST("matching/V0/LambdaPtMass"), particle.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);

      // Reflection
      double reflectedMass = getReflectedMass(v0, true);
      registry.fill(HIST("matching/V0/LambdaReflection"), particle.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
    } else if (particle.pdgCode() == PDG_t::kLambda0Bar) { // AntiLambda
      registry.fill(HIST("matching/V0/AntiLambdaPtEtaPhi"), particle.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/V0/AntiLambdaPtCtauMass"), particle.pt(), v0.pt(), ctauAntiLambda, v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/V0/AntiLambdaPtRadiusCosPA"), particle.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/V0/AntiLambdaPtDCAposneg"), particle.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/V0/AntiLambdaPtDCAd"), particle.pt(), v0.pt(), v0.dcaV0daughters(), weight);
      registry.fill(HIST("matching/V0/AntiLambdaPtMass"), particle.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);

      // Reflection
      double reflectedMass = getReflectedMass(v0, false);
      registry.fill(HIST("matching/V0/AntiLambdaReflection"), particle.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
    }
  }
  template <typename V0DaughterType, typename ParticleDaughterType, typename V0Type, typename ParticleType> // Reconstructed signal for inclusive V0s: daughters
  void fillMatchedV0DaughtersInclusive(V0Type const& v0, ParticleType const& pv0, double weight = 1.)
  {
    auto negTrack = v0.template negTrack_as<V0DaughterType>();
    auto posTrack = v0.template posTrack_as<V0DaughterType>();
    auto negPart = negTrack.template mcParticle_as<ParticleDaughterType>();
    auto posPart = posTrack.template mcParticle_as<ParticleDaughterType>();
    registry.fill(HIST("matching/V0/V0PosPartPtRatioPtRelDiffPt"), posPart.pt(), posTrack.pt() / posPart.pt(), (posTrack.pt() - posPart.pt()) / posPart.pt(), weight);
    registry.fill(HIST("matching/V0/V0NegPartPtRatioPtRelDiffPt"), negPart.pt(), negTrack.pt() / negPart.pt(), (negTrack.pt() - negPart.pt()) / negPart.pt(), weight);

    if (std::abs(v0.pdgCode()) == PDG_t::kK0Short) { // K0S
      registry.fill(HIST("matching/V0/K0SPosNegPtMass"), pv0.pt(), posPart.pt(), negPart.pt(), v0.mK0Short(), weight);
    } else if (v0.pdgCode() == PDG_t::kLambda0) { // Lambda
      registry.fill(HIST("matching/V0/LambdaPosNegPtMass"), pv0.pt(), posPart.pt(), negPart.pt(), v0.mLambda(), weight);
    } else if (v0.pdgCode() == PDG_t::kLambda0Bar) { // AntiLambda
      registry.fill(HIST("matching/V0/AntiLambdaPosNegPtMass"), pv0.pt(), posPart.pt(), negPart.pt(), v0.mAntiLambda(), weight);
    }
  }
  template <typename DetJet, typename PartJet> // Reconstructed jets
  void fillMatchedJet(DetJet const& detJet, PartJet const& partJet, double weight = 1.)
  {
    double deltaEta = detJet.eta() - partJet.eta();
    double deltaPhi = RecoDecay::constrainAngle(detJet.phi() - partJet.phi(), -constants::math::PI);
    double dR = jetutilities::deltaR(detJet, partJet);

    registry.fill(HIST("matching/jets/inclMatchDetJetPtEtaPhi"), detJet.pt(), detJet.eta(), detJet.phi(), weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtEtaPhi"), partJet.pt(), partJet.eta(), partJet.phi(), weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtEtaPhiMatchDist"), partJet.pt(), partJet.eta(), partJet.phi(), dR, weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtEnergyScale"), partJet.pt(), detJet.pt() / partJet.pt(), weight);
    registry.fill(HIST("matching/jets/inclMatchDetJetPtPartJetPt"), detJet.pt(), partJet.pt(), weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtDetJetEtaPartJetEta"), partJet.pt(), detJet.eta(), partJet.eta(), weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtDetJetPhiPartJetPhi"), partJet.pt(), detJet.phi(), partJet.phi(), weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtResolutionPt"), partJet.pt(), (detJet.pt() - partJet.pt()), weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtResolutionEta"), partJet.pt(), partJet.eta(), deltaEta, weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtResolutionPhi"), partJet.pt(), partJet.phi(), deltaPhi, weight);
    registry.fill(HIST("matching/jets/inclMatchPartJetPtRelDiffPt"), partJet.pt(), (detJet.pt() - partJet.pt()) / partJet.pt(), weight);

    if (!jetContainsV0s(detJet))
      return;

    registry.fill(HIST("matching/jets/matchDetJetPtEtaPhi"), detJet.pt(), detJet.eta(), detJet.phi(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtEtaPhi"), partJet.pt(), partJet.eta(), partJet.phi(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtEtaPhiMatchDist"), partJet.pt(), partJet.eta(), partJet.phi(), dR, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtEnergyScale"), partJet.pt(), detJet.pt() / partJet.pt(), weight);
    registry.fill(HIST("matching/jets/matchDetJetPtPartJetPt"), detJet.pt(), partJet.pt(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtDetJetEtaPartJetEta"), partJet.pt(), detJet.eta(), partJet.eta(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtDetJetPhiPartJetPhi"), partJet.pt(), detJet.phi(), partJet.phi(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionPt"), partJet.pt(), (detJet.pt() - partJet.pt()), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionEta"), partJet.pt(), partJet.eta(), deltaEta, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionPhi"), partJet.pt(), partJet.phi(), deltaPhi, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtRelDiffPt"), partJet.pt(), (detJet.pt() - partJet.pt()) / partJet.pt(), weight);
  }
  template <typename CollisionType, typename DetJetType, typename PartJetType, typename V0Type, typename ParticleType> // Reconstructed signal for in-jet V0s
  void fillMatchedV0InJet(CollisionType const& coll, DetJetType const& detJet, PartJetType const& partJet, V0Type const& v0, ParticleType const& particle, double weight = 1.)
  {
    bool correctCollision = (coll.mcCollisionId() == particle.mcCollisionId());
    double detTrackProj = getMomFrac(detJet, v0);
    double partTrackProj = getMomFrac(partJet, particle);

    double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

    registry.fill(HIST("matching/jets/V0/matchDetJetPtV0TrackProjPartJetPtV0TrackProj"), detJet.pt(), detTrackProj, partJet.pt(), partTrackProj, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0Pt"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetPt"), partJet.pt(), particle.pt(), detJet.pt(), weight);
    // registry.fill(HIST("matching/jets/V0/partJetPtDetJetPtPartV0PtRatioPtRelDiffPt"), partJet.pt(), detJet.pt(), particle.pt(), v0.pt() / particle.pt(), (v0.pt() - particle.pt()) / particle.pt(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtDetJetPtPartV0PtRelDiffPt"), partJet.pt(), detJet.pt(), particle.pt(), (v0.pt() - particle.pt()) / particle.pt(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtDetJetPtPartV0ZRelDiffZ"), partJet.pt(), detJet.pt(), partTrackProj, (detTrackProj - partTrackProj) / partTrackProj, weight);

    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauLambda, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauAntiLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauAntiLambda, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauK0s, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mLambda(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassAntiLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mAntiLambda(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtRadius"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0radius(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCosPA"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0cosPA(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtDCAposneg"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtDCAd"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauLambda, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauAntiLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauAntiLambda, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauK0s, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mLambda(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassAntiLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mAntiLambda(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjRadius"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0radius(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCosPA"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0cosPA(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjDCAposneg"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjDCAd"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcaV0daughters(), weight);

    if (std::abs(particle.pdgCode()) == PDG_t::kK0Short) { // K0S
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPt"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProj"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      if (correctCollision) {
        registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtRightCollision"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
        registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjRightCollision"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      } else {
        registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtWrongCollision"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
        registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjWrongCollision"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      }

      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauAntiLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassAntiLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtAllMasses"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtRadius"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCosPA"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtDCAposneg"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtDCAd"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauAntiLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassAntiLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjAllMasses"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjRadius"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCosPA"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjDCAposneg"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjDCAd"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcaV0daughters(), weight);
    } else if (particle.pdgCode() == PDG_t::kLambda0) { // Lambda
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPt"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProj"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      if (correctCollision) {
        registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtRightCollision"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
        registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjRightCollision"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      } else {
        registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtWrongCollision"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
        registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjWrongCollision"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      }

      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtCtauLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtCtauAntiLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtCtauK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtMassLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtMassAntiLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtMassK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtAllMasses"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtRadius"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtCosPA"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtDCAposneg"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtDCAd"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjCtauLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjCtauAntiLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjCtauK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjMassLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjMassAntiLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjMassK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjAllMasses"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjRadius"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjCosPA"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjDCAposneg"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjDCAd"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcaV0daughters(), weight);

      // Reflection
      double reflectedMass = getReflectedMass(v0, true);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaPtDetJetPtLambdaPtLambdaReflection"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambdaTrackProjDetJetPtLambdaTrackProjLambdaReflection"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
    } else if (particle.pdgCode() == PDG_t::kLambda0Bar) { // AntiLambda
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPt"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProj"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      if (correctCollision) {
        registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtRightCollision"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
        registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjRightCollision"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      } else {
        registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtWrongCollision"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
        registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjWrongCollision"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, weight);
      }

      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtCtauLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtCtauAntiLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtCtauK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtMassLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtMassAntiLambda"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtMassK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtAllMasses"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtRadius"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtCosPA"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtDCAposneg"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtDCAd"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjCtauLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjCtauAntiLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjCtauK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjMassLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjMassAntiLambda"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjMassK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjAllMasses"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjRadius"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjCosPA"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjDCAposneg"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjDCAd"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcaV0daughters(), weight);

      // Reflection
      double reflectedMass = getReflectedMass(v0, false);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaPtDetJetPtAntiLambdaPtAntiLambdaReflection"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambdaTrackProjDetJetPtAntiLambdaTrackProjAntiLambdaReflection"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
    } // AntiLambda
  }
  template <typename V0DaughterType, typename ParticleDaughterType, typename DetJetType, typename PartJetType, typename V0Type, typename ParticleType> // Reconstructed signal for in-jet V0s: daughters
  void fillMatchedV0DaughtersInJet(DetJetType const& detJet, PartJetType const& partJet, V0Type const& v0, ParticleType const& particle, double weight = 1.)
  {
    auto negTrack = v0.template negTrack_as<V0DaughterType>();
    auto posTrack = v0.template posTrack_as<V0DaughterType>();
    auto negPart = negTrack.template mcParticle_as<ParticleDaughterType>();
    auto posPart = posTrack.template mcParticle_as<ParticleDaughterType>();
    registry.fill(HIST("matching/jets/V0/partJetPtDetJetPtPartV0PtPosPtRatioPtRelDiffPt"), partJet.pt(), detJet.pt(), particle.pt(), posPart.pt(), posTrack.pt() / posPart.pt(), (posTrack.pt() - posPart.pt()) / posPart.pt(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtDetJetPtPartV0PtNegPtRatioPtRelDiffPt"), partJet.pt(), detJet.pt(), particle.pt(), negPart.pt(), negTrack.pt() / negPart.pt(), (negTrack.pt() - negPart.pt()) / negPart.pt(), weight);
  }
  // Misses - Counts (or event weights)
  template <typename T>
  void fillMissV0Inclusive(T const& pv0, double weight = 1.)
  {
    int pdg = pv0.pdgCode();
    registry.fill(HIST("matching/V0/missV0PtEtaPhi"), pv0.pt(), pv0.eta(), pv0.phi(), weight);
    if (std::abs(pdg) == PDG_t::kK0Short) { // K0S
      registry.fill(HIST("matching/V0/missK0SPtEtaPhi"), pv0.pt(), pv0.eta(), pv0.phi(), weight);
    } else if (pdg == PDG_t::kLambda0) { // Lambda
      registry.fill(HIST("matching/V0/missLambdaPtEtaPhi"), pv0.pt(), pv0.eta(), pv0.phi(), weight);
    } else if (pdg == PDG_t::kLambda0Bar) { // AntiLambda
      registry.fill(HIST("matching/V0/missAntiLambdaPtEtaPhi"), pv0.pt(), pv0.eta(), pv0.phi(), weight);
    }
  }
  template <typename T>
  void fillMissJet(T const& jet, double weight = 1.)
  {
    registry.fill(HIST("matching/jets/inclMissJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    if (!jetContainsV0s(jet))
      return;

    registry.fill(HIST("matching/jets/missJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
  }
  template <typename JetType, typename V0Type>
  void fillMissV0InJet(JetType const& jet, V0Type const& v0, double weight = 1.)
  {
    double trackProj = getMomFrac(jet, v0);

    registry.fill(HIST("matching/jets/V0/missJetPtV0TrackProj"), jet.pt(), trackProj, weight);
    registry.fill(HIST("matching/jets/V0/missJetPtV0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
    if (std::abs(v0.pdgCode()) == PDG_t::kK0Short) { // K0S
      registry.fill(HIST("matching/jets/V0/missJetPtK0SPtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/missJetPtK0STrackProj"), jet.pt(), trackProj, weight);
    } else if (v0.pdgCode() == PDG_t::kLambda0) { // Lambda
      registry.fill(HIST("matching/jets/V0/missJetPtLambdaPtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/missJetPtLambdaTrackProj"), jet.pt(), trackProj, weight);
    } else if (v0.pdgCode() == PDG_t::kLambda0Bar) { // AntiLambda
      registry.fill(HIST("matching/jets/V0/missJetPtAntiLambdaPtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/missJetPtAntiLambdaTrackProj"), jet.pt(), trackProj, weight);
    }
  }
  // Fakes - Counts (or event weights)
  template <typename T, typename U>
  void fillFakeV0Inclusive(T const& coll, U const& v0, double weight = 1.)
  {
    double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;

    double massDiff = v0.mLambda() - v0.mAntiLambda();
    double massRatio = v0.mAntiLambda() / v0.mLambda();
    double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

    registry.fill(HIST("matching/V0/fakeV0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
    registry.fill(HIST("matching/V0/fakeV0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
    registry.fill(HIST("matching/V0/fakeV0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
    registry.fill(HIST("matching/V0/fakeV0PtLambdaMasses"), v0.pt(), v0.mLambda() - v0.mAntiLambda(), v0.mAntiLambda() / v0.mLambda(), (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda(), weight);
    registry.fill(HIST("matching/V0/fakeV0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
    registry.fill(HIST("matching/V0/fakeV0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
    registry.fill(HIST("matching/V0/fakeV0PtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

    if (v0.isK0SCandidate()) {
      registry.fill(HIST("matching/V0/fakeK0SPtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/V0/fakeK0SPtCtauMass"), v0.pt(), ctauK0s, v0.mK0Short(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mK0Short(), weight);
    }
    if (v0.isLambdaCandidate()) {
      registry.fill(HIST("matching/V0/fakeLambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/V0/fakeLambdaPtCtauMass"), v0.pt(), ctauLambda, v0.mLambda(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mLambda(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mLambda(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mLambda(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mLambda(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mLambda(), weight);
    }
    if (v0.isAntiLambdaCandidate()) {
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtDCAd"), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/V0/fakeAntiLambdaPtCtauMass"), v0.pt(), ctauAntiLambda, v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtRadiusMass"), v0.pt(), v0.v0radius(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtCosPAMass"), v0.pt(), v0.v0cosPA(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtDCAposMass"), v0.pt(), v0.dcapostopv(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtDCAnegMass"), v0.pt(), v0.dcanegtopv(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPtDCAdMass"), v0.pt(), v0.dcaV0daughters(), v0.mAntiLambda(), weight);
    }
  }
  template <typename T, typename U>
  void fillFakeV0DaughtersInclusive(U const& v0, double weight = 1.)
  {
    auto negTrack = v0.template negTrack_as<T>();
    auto posTrack = v0.template posTrack_as<T>();
    registry.fill(HIST("matching/V0/fakeV0PosTrackPtEtaPhi"), posTrack.pt(), posTrack.eta(), posTrack.phi(), weight);
    registry.fill(HIST("matching/V0/fakeV0NegTrackPtEtaPhi"), negTrack.pt(), negTrack.eta(), negTrack.phi(), weight);

    if (v0.isK0SCandidate()) {
      registry.fill(HIST("matching/V0/fakeK0SPosTrackPtEtaPhi"), posTrack.pt(), posTrack.eta(), posTrack.phi(), weight);
      registry.fill(HIST("matching/V0/fakeK0SPosTrackPtMass"), v0.pt(), posTrack.pt(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/V0/fakeK0SNegTrackPtEtaPhi"), negTrack.pt(), negTrack.eta(), negTrack.phi(), weight);
      registry.fill(HIST("matching/V0/fakeK0SNegTrackPtMass"), v0.pt(), negTrack.pt(), v0.mK0Short(), weight);
    }
    if (v0.isLambdaCandidate()) {
      registry.fill(HIST("matching/V0/fakeLambdaPosTrackPtEtaPhi"), posTrack.pt(), posTrack.eta(), posTrack.phi(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaPosTrackPtMass"), v0.pt(), posTrack.pt(), v0.mLambda(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaNegTrackPtEtaPhi"), negTrack.pt(), negTrack.eta(), negTrack.phi(), weight);
      registry.fill(HIST("matching/V0/fakeLambdaNegTrackPtMass"), v0.pt(), negTrack.pt(), v0.mLambda(), weight);
    }
    if (v0.isAntiLambdaCandidate()) {
      registry.fill(HIST("matching/V0/fakeAntiLambdaPosTrackPtEtaPhi"), posTrack.pt(), posTrack.eta(), posTrack.phi(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaPosTrackPtMass"), v0.pt(), posTrack.pt(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaNegTrackPtEtaPhi"), negTrack.pt(), negTrack.eta(), negTrack.phi(), weight);
      registry.fill(HIST("matching/V0/fakeAntiLambdaNegTrackPtMass"), v0.pt(), negTrack.pt(), v0.mAntiLambda(), weight);
    }
  }
  template <typename T, typename U, typename V> // Check if inclusive V0 was missed because daughter decayed
  void fillFakeV0DecayedInclusive(V const& v0, double weight = 1.)
  {
    // Check if decayed daughter
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();

    auto posPart = posTrack.template mcParticle_as<U>();
    auto negPart = negTrack.template mcParticle_as<U>();

    auto posMom = posPart.template mothers_first_as<U>();
    auto negMom = negPart.template mothers_first_as<U>();

    bool posDecayed = false;
    bool negDecayed = false;

    // This should not happen. They should have been matched
    if (posMom == negMom) {
      registry.fill(HIST("matching/V0/nonedecayedFakeV0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      return;
    }

    if (posMom.has_mothers()) {
      auto posGrandMom = posMom.template mothers_first_as<U>();
      if (posGrandMom == negMom) {
        posDecayed = true;
      }
    }
    if (negMom.has_mothers()) {
      auto negGrandMom = negMom.template mothers_first_as<U>();
      if (negGrandMom == posMom) {
        negDecayed = true;
      }
    }

    // This shouldn't happen
    if (posDecayed && negDecayed) {
      registry.fill(HIST("matching/V0/doubledecayedFakeV0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      return;
    }
    if (posDecayed || negDecayed) {
      double pt = posDecayed ? negMom.pt() : posMom.pt();
      int pdg = posDecayed ? negMom.pdgCode() : posMom.pdgCode();

      if (std::abs(pdg) == PDG_t::kK0Short) {
        registry.fill(HIST("matching/V0/decayedK0SV0PtMass"), pt, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      } else if (pdg == PDG_t::kLambda0) {
        registry.fill(HIST("matching/V0/decayedLambdaV0PtMass"), pt, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      } else if (pdg == PDG_t::kLambda0Bar) {
        registry.fill(HIST("matching/V0/decayedAntiLambdaV0PtMass"), pt, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      } else {
        registry.fill(HIST("matching/V0/decayedOtherPtV0PtMass"), pt, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      }
    }
  }
  template <typename T>
  void fillFakeJet(T const& jet, double weight = 1.)
  {
    registry.fill(HIST("matching/jets/fakeJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    if (!jetContainsV0s(jet))
      return;

    registry.fill(HIST("matching/jets/inclFakeJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
  }
  template <typename CollisionType, typename JetType, typename V0Type>
  void fillFakeV0InJet(CollisionType const& coll, JetType const& jet, V0Type const& v0, double weight = 1.)
  {
    double trackProj = getMomFrac(jet, v0);
    double ctauLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * constants::physics::MassK0Short;
    double massDiff = v0.mLambda() - v0.mAntiLambda();
    double massRatio = v0.mAntiLambda() / v0.mLambda();
    double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

    registry.fill(HIST("matching/jets/V0/fakeJetPtV0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0TrackProj"), jet.pt(), trackProj, weight);

    registry.fill(HIST("matching/jets/V0/fakeJetPtV0PtCtau"), jet.pt(), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0PtMass"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0PtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff, weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0PtRadiusCosPA"), jet.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0PtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0PtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

    registry.fill(HIST("matching/jets/V0/fakeJetPtV0TrackProjCtau"), jet.pt(), trackProj, ctauK0s, ctauLambda, ctauAntiLambda, weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0TrackProjMass"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0TrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff, weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0TrackProjRadiusCosPA"), jet.pt(), trackProj, v0.v0radius(), v0.v0cosPA(), weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0TrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
    registry.fill(HIST("matching/jets/V0/fakeJetPtV0TrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters(), weight);

    if (v0.isK0SCandidate()) {
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0SPtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0STrackProj"), jet.pt(), trackProj, weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtK0SPtCtau"), jet.pt(), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0SPtMass"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0SPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0SPtRadiusCosPA"), jet.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0SPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0SPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtK0STrackProjCtau"), jet.pt(), trackProj, ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0STrackProjMass"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0STrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0STrackProjRadiusCosPA"), jet.pt(), trackProj, v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0STrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtK0STrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters(), weight);
    }
    if (v0.isLambdaCandidate()) {
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaPtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaTrackProj"), jet.pt(), trackProj, weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaPtCtau"), jet.pt(), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaPtMass"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaPtRadiusCosPA"), jet.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaTrackProjEtaPhi"), jet.pt(), trackProj, v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaTrackProjCtau"), jet.pt(), trackProj, ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaTrackProjMass"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaTrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaTrackProjRadiusCosPA"), jet.pt(), trackProj, v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaTrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambdaTrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters(), weight);
    }
    if (v0.isAntiLambdaCandidate()) {
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaPtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaTrackProj"), jet.pt(), trackProj, weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaPtCtau"), jet.pt(), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaPtMass"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaPtRadiusCosPA"), jet.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaTrackProjCtau"), jet.pt(), trackProj, ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaTrackProjMass"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaTrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaTrackProjRadiusCosPA"), jet.pt(), trackProj, v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaTrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambdaTrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters(), weight);
    }
  }
  template <typename T, typename U, typename V, typename W, typename X> // Check if V0 in jet was missed because daughter decayed
  void fillFakeV0DecayedInJet(V const& partJet, W const& detJet, X const& v0, double weight = 1.)
  {
    // Check if decayed daughter
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();

    auto posPart = posTrack.template mcParticle_as<U>();
    auto negPart = negTrack.template mcParticle_as<U>();

    auto posMom = posPart.template mothers_first_as<U>();
    auto negMom = negPart.template mothers_first_as<U>();

    bool posDecayed = false;
    bool negDecayed = false;

    double zv0 = getMomFrac(detJet, v0);

    // This should not happen. They should have been matched
    if (posMom == negMom) {
      registry.fill(HIST("matching/jets/V0/nonedecayedFakeV0PtMass"), partJet.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/nonedecayedFakeV0TrackProjMass"), partJet.pt(), detJet.pt(), zv0, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      return;
    }

    if (posMom.has_mothers()) {
      auto posGrandMom = posMom.template mothers_first_as<U>();
      if (posGrandMom == negMom) {
        posDecayed = true;
      }
    }
    if (negMom.has_mothers()) {
      auto negGrandMom = negMom.template mothers_first_as<U>();
      if (negGrandMom == posMom) {
        negDecayed = true;
      }
    }

    // This shouldn't happen
    if (posDecayed && negDecayed) {
      registry.fill(HIST("matching/jets/V0/doubledecayedFakeV0PtMass"), partJet.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/doubledecayedFakeV0TrackProjMass"), partJet.pt(), detJet.pt(), zv0, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      return;
    }
    if (posDecayed || negDecayed) {
      double pt = posDecayed ? negMom.pt() : posMom.pt();
      int pdg = posDecayed ? negMom.pdgCode() : posMom.pdgCode();

      double z = 0.;
      bool partIsInJet = false;
      for (auto const& part : partJet.template tracks_as<U>()) {
        if (posDecayed && (part == negMom)) {
          partIsInJet = true;
          z = getMomFrac(partJet, part);
          break;
        }
        if (negDecayed && (part == posMom)) {
          partIsInJet = true;
          z = getMomFrac(partJet, part);
          break;
        }
      }

      if (std::abs(pdg) == PDG_t::kK0Short) {
        registry.fill(HIST("matching/jets/V0/decayedK0SV0PtMass"), partJet.pt(), detJet.pt(), pt, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
        if (partIsInJet) {
          registry.fill(HIST("matching/jets/V0/decayedK0SV0TrackProjMass"), partJet.pt(), detJet.pt(), z, zv0, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
        }
      } else if (pdg == PDG_t::kLambda0) {
        registry.fill(HIST("matching/jets/V0/decayedLambdaV0PtMass"), partJet.pt(), detJet.pt(), pt, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
        if (partIsInJet) {
          registry.fill(HIST("matching/jets/V0/decayedLambdaV0TrackProjMass"), partJet.pt(), detJet.pt(), z, zv0, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
        }
      } else if (pdg == PDG_t::kLambda0Bar) {
        registry.fill(HIST("matching/jets/V0/decayedAntiLambdaV0PtMass"), partJet.pt(), detJet.pt(), pt, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
        if (partIsInJet) {
          registry.fill(HIST("matching/jets/V0/decayedAntiLambdaV0TrackProjMass"), partJet.pt(), detJet.pt(), z, zv0, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
        }
      } else {
        registry.fill(HIST("matching/jets/V0/decayedOtherPtV0PtMass"), partJet.pt(), detJet.pt(), pt, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
        if (partIsInJet) {
          registry.fill(HIST("matching/jets/V0/decayedOtherPtV0TrackProjMass"), partJet.pt(), detJet.pt(), z, zv0, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
        }
      }
    }
  }

  // ---------------------------------------------------
  // Processes
  // ---------------------------------------------------
  void processDummy(aod::JetTracks const&) {}
  PROCESS_SWITCH(JetFragmentation, processDummy, "Dummy process function turned on by default", true);

  void processDataV0(soa::Filtered<aod::JetCollisions>::iterator const& coll, DataV0JetsWithConstituents const& jets, aod::CandidatesV0Data const& V0s, aod::JetTracks const&)
  {
    registry.fill(HIST("data/hEvents"), 0.5);
    if (!jetderiveddatautilities::selectCollision(coll, eventSelectionBits))
      return;

    registry.fill(HIST("data/hEvents"), 1.5);
    registry.fill(HIST("data/V0/nV0sEvent"), V0s.size());

    if (fillHistsInclusiveV0s) {
      fillDataV0sInclusive(coll, V0s);
      fillDataV0sInclusiveWeighted(coll, V0s);
    }

    if (!fillHistsJets)
      return;

    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, -99., -99., v0EtaMin, v0EtaMax))
        continue;

      fillDataJet(jet);
      fillDataV0sInPerpCone(coll, jet, V0s);

      int nV0inJet = 0, nLambdainJet = 0, nAntiLambdainJet = 0, nK0SinJet = 0;   // Counts
      float wV0inJet = 0, wLambdainJet = 0, wAntiLambdainJet = 0, wK0SinJet = 0; // Weights

      std::vector<double> values;
      std::vector<std::vector<double>> weights;
      for (const auto& v0 : jet.candidates_as<aod::CandidatesV0Data>()) {
        if (v0.isRejectedCandidate())
          continue;

        float signalProb = getV0SignalProb(v0);
        nV0inJet++;
        wV0inJet += signalProb;
        if (v0.isK0SCandidate()) {
          nK0SinJet++;
          wK0SinJet += signalProb;
        }
        if (v0.isLambdaCandidate()) {
          nLambdainJet++;
          wLambdainJet += signalProb;
        }
        if (v0.isAntiLambdaCandidate()) {
          nAntiLambdainJet++;
          wAntiLambdainJet += signalProb;
        }

        fillDataV0sInJet(coll, jet, v0);
        double z = getMomFrac(jet, v0);
        std::vector<double> w;

        if (nV0Classes == 2)
          w = getV0SignalProbVector2Classes(v0);
        else if (nV0Classes == 4)
          w = getV0SignalProbVector4Classes(v0);
        else
          return;

        values.push_back(z);
        weights.push_back(w);
      }
      values.push_back(jet.pt());
      registry.fill(HIST("data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda"), jet.pt(), nV0inJet, nK0SinJet, nLambdainJet, nAntiLambdainJet);
      registry.fill(HIST("data/jets/weighted/V0/jetPtnV0nK0SnLambdanAntiLambda"), jet.pt(), wV0inJet, wK0SinJet, wLambdainJet, wAntiLambdainJet);

      if (nV0inJet == 0)
        continue;

      int nStates = std::round(std::pow(static_cast<double>(nV0Classes), static_cast<double>(nV0inJet)));
      for (int M = 0; M < nStates; M++) {
        std::vector<int> state = convertState(M, nV0inJet, nV0Classes);
        std::vector<double> corrected;
        if (doCorrectionWithTracks)
          corrected = correctedValuesPlusTracks<aod::CandidatesV0Data, aod::JetTracks>(state, jet);
        else
          corrected = correctedValues(state, values);

        double ws = stateWeight(state, weights);
        double jetpt = corrected[nV0inJet];
        fillDataJetWeighted(jetpt, jet.eta(), jet.phi(), ws);
        fillDataV0sInJetWeighted(coll, jet, state, corrected, ws);
      }
    }
  }
  PROCESS_SWITCH(JetFragmentation, processDataV0, "Data V0", false);

  void processMcV0(soa::Filtered<aod::JetCollisionsMCD>::iterator const& coll, aod::JetMcCollisions const&, MatchedMCDV0JetsWithConstituents const& v0jetsMCD, MatchedMCPV0JetsWithConstituents const& v0jetsMCP, CandidatesV0MCDWithLabels const& V0s, aod::CandidatesV0MCP const& pV0s, aod::JetTracksMCD const&, aod::JetParticles const&, aod::McParticles const& particles)
  {
    registry.fill(HIST("matching/hEvents"), 0.5);
    if (!coll.has_mcCollision())
      return;
    if (!jetderiveddatautilities::selectCollision(coll, eventSelectionBits))
      return;

    registry.fill(HIST("matching/hEvents"), 1.5);
    double weight = coll.mcCollision().weight();
    registry.fill(HIST("matching/hEvents"), 2.5, weight);
    registry.fill(HIST("matching/V0/nV0sEvent"), V0s.size());
    registry.fill(HIST("matching/V0/nV0sEventWeighted"), V0s.size(), weight);

    if (fillHistsInclusiveV0s) {
      fillMcdV0sInclusive(coll, V0s, weight);
      fillMcpV0sInclusive(pV0s, weight);
      fillMatchingV0sInclusive<aod::JetTracksMCD, aod::JetParticles>(coll, V0s, pV0s, weight);
    }

    if (!fillHistsJets)
      return;

    for (const auto& detJet : v0jetsMCD) {
      if (!jetfindingutilities::isInEtaAcceptance(detJet, -99., -99., v0EtaMin, v0EtaMax))
        continue;

      fillMcdJet(detJet, weight);
      fillMcV0sInPerpCone(coll, detJet, V0s, particles, weight);

      int nV0inJet = 0, nLambdainJet = 0, nAntiLambdainJet = 0, nK0SinJet = 0;
      if (!detJet.has_matchedJetGeo()) {
        fillFakeJet(detJet, weight);
        for (const auto& v0 : detJet.candidates_as<CandidatesV0MCDWithLabels>()) {
          fillFakeV0InJet(coll, detJet, v0, weight);
        }
        continue;
      } // if jet not matched

      for (const auto& partJet : detJet.template matchedJetGeo_as<MatchedMCPV0JetsWithConstituents>()) {
        fillMatchedJet(detJet, partJet, weight);
        fillMcV0sInMatchedPerpCone(coll, detJet, partJet, V0s, particles, weight);

        for (const auto& detV0 : detJet.candidates_as<CandidatesV0MCDWithLabels>()) {
          if (!detV0.has_mcParticle()) {
            fillFakeV0InJet(coll, detJet, detV0, weight);
            fillFakeV0DecayedInJet<aod::JetTracksMCD, aod::JetParticles>(partJet, detJet, detV0, weight);
            continue;
          }

          bool isV0Matched = false;
          for (const auto& partV0 : partJet.template candidates_as<aod::CandidatesV0MCP>()) {
            if (!v0sAreMatched<aod::JetTracksMCD>(detV0, partV0))
              continue;

            isV0Matched = true;
            nV0inJet++;
            fillMatchedV0InJet(coll, detJet, partJet, detV0, partV0, weight);
            fillMatchedV0DaughtersInJet<aod::JetTracksMCD, aod::JetParticles>(detJet, partJet, detV0, partV0, weight);

            if (std::abs(partV0.pdgCode()) == PDG_t::kK0Short) {
              nK0SinJet++;
            } else if (partV0.pdgCode() == PDG_t::kLambda0) {
              nLambdainJet++;
            } else if (partV0.pdgCode() == PDG_t::kLambda0Bar) {
              nAntiLambdainJet++;
            }
            break;
          } // partV0 loop

          if (!isV0Matched) {
            fillFakeV0InJet(coll, detJet, detV0, weight);
          }
        } // detV0 loop
        registry.fill(HIST("matching/jets/V0/jetPtnV0MatchednK0SnLambdanAntiLambda"), partJet.pt(), nV0inJet, nK0SinJet, nLambdainJet, nAntiLambdainJet, weight);
      } // Matched partJet loop
    } // detJet loop

    for (const auto& partJet : v0jetsMCP) {
      fillMcpJet(partJet, weight);

      if (!partJet.has_matchedJetGeo()) {
        fillMissJet(partJet, weight);
        for (const auto& partV0 : partJet.candidates_as<aod::CandidatesV0MCP>()) {
          fillMissV0InJet(partJet, partV0, weight);
        }
        continue;
      } // if jet not matched

      bool isJetMatched = false;
      for (const auto& detJet : partJet.template matchedJetGeo_as<MatchedMCDV0JetsWithConstituents>()) {
        if (!jetfindingutilities::isInEtaAcceptance(detJet, -99., -99., v0EtaMin, v0EtaMax))
          continue;

        isJetMatched = true;
        for (const auto& partV0 : partJet.candidates_as<aod::CandidatesV0MCP>()) {
          bool isV0Matched = false;
          for (const auto& detV0 : detJet.candidates_as<CandidatesV0MCDWithLabels>()) {
            if (v0sAreMatched<aod::JetTracksMCD>(detV0, partV0)) {
              isV0Matched = true;
              break;
            }
          } // detV0 loop

          // If V0 is matched, it has already been filled in the mcdjet loop
          if (!isV0Matched)
            fillMissV0InJet(partJet, partV0, weight);
        } // partV0 loop
      } // detJet loop

      // To account for matched jets where the detector level jet is outside of the eta range (cut applied within this task)
      if (!isJetMatched) {
        for (const auto& partV0 : partJet.candidates_as<aod::CandidatesV0MCP>()) {
          fillMissV0InJet(partJet, partV0, weight);
        }
      }
    } // partJet loop
  }
  PROCESS_SWITCH(JetFragmentation, processMcV0, "MC V0", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetFragmentation>(cfgc)};
}
