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

// jet trigger QA task
//
/// \author Gijs van Weelden <g.van.weelden@cern.ch>
//

#include "TH1F.h"
#include "TTree.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using McDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
using McPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

struct JetFragmentation {
  HistogramRegistry registry{"registry"};

  std::vector<int> pdgVector = {211, 321, 2212, 111, 130, 310, 311, 3122};
  std::vector<std::string> hadronVector = {"#pi^{#pm}", "#it{K}^{#pm}", "#it{p}^{#pm}", "#pi^{0}", "#it{K}^{0}_{L}", "#it{K}^{0}_{S}", "#it{K}^{0}", "#Lambda^{0}"};

  Configurable<double> matchedDetJetEtaMin{"matchedDetJetEtaMin", -0.5, "minimum matchedDetJet eta"};
  Configurable<double> matchedDetJetEtaMax{"matchedDetJetEtaMax", 0.5, "maximum matchedDetJet eta"};
  Configurable<double> dataJetEtaMin{"dataJetEtaMin", -0.5, "minimum data jet eta"};
  Configurable<double> dataJetEtaMax{"dataJetEtaMax", 0.5, "maximum data jet eta"};
  Configurable<double> v0EtaMin{"v0EtaMin", -0.75, "minimum data V0 eta"};
  Configurable<double> v0EtaMax{"v0EtaMax", 0.75, "maximum data V0 eta"};

  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<double> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<double> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<double> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<double> v0radius{"v0radius", 5.0, "V0 Radius"};
  Configurable<double> lifetimecutK0S{"lifetimecutK0S", 20., "lifetimecutK0S"};
  Configurable<double> lifetimecutLambda{"lifetimecutLambda", 25., "lifetimecutLambda"};
  Configurable<double> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5., "TpcPidNsigmaCut"};

  Configurable<double> k0sMassAccWindow{"k0sMassAccWindow", 0.03, "k0sMassAccWindow"};
  Configurable<double> lambdaMassAccWindow{"lambdaMassAccWindow", 0.01, "lambdaMassAccWindow"};
  Configurable<double> antilambdaMassAccWindow{"antilambdaMassAccWindow", 0.01, "antilambdaMassAccWindow"};
  Configurable<double> k0sMassRejWindow{"k0sMassRejWindow", 0.01, "k0sMassRejWindow"};
  Configurable<double> lambdaMassRejWindow{"lambdaMassRejWindow", 0.005, "lambdaMassRejWindow"};
  Configurable<double> antilambdaMassRejWindow{"antilambdaMassRejWindow", 0.005, "antilambdaMassRejWindow"};

  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {40, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {20, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binZ{"binZ", {40, 0.0001f, 1.0001f}, ""};
  ConfigurableAxis binXi{"binXi", {50, 0.f, 10.f}, ""};
  ConfigurableAxis binTheta{"binTheta", {40, -0.05f, 0.395f}, ""};
  ConfigurableAxis binJetR{"binJetR", {6, 0.05f, 0.65f}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binPDG{"binPDG", {static_cast<double>(pdgVector.size()), -0.5f, static_cast<double>(pdgVector.size()) - 0.5f}, ""};
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

  ConfigurableAxis binV0Pt{"binV0Pt", {200, 0.0f, 10.0f}, ""};
  ConfigurableAxis binV0Eta{"binV0Eta", {20, -1.f, 1.f}, ""};
  ConfigurableAxis binV0Phi{"binV0Phi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binV0Ctau{"binV0Ctau", {200, 0.0f, 40.0f}, ""};
  ConfigurableAxis binV0Radius{"binV0Radius", {100, 0.0f, 100.0f}, ""};
  ConfigurableAxis binV0CosPA{"binV0CosPA", {100, 0.95f, 1.0f}, ""};
  ConfigurableAxis binV0DCA{"binV0DCA", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis binV0DCAp{"binV0DCAp", {100, -10.0f, 10.0f}, ""};
  ConfigurableAxis binV0DCAn{"binV0DCAn", {100, -10.0f, 10.0f}, ""};
  ConfigurableAxis binV0DCAd{"binV0DCAd", {100, 0.0f, 10.0f}, ""};

  ConfigurableAxis binK0SMass{"binK0SMass", {400, 0.400f, 0.600f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis binLambdaMass{"binLambdaMass", {200, 1.015f, 1.215f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis binLambdaMassDiff{"binLambdaMassDiff", {200, -0.199f, 0.201f}, "M(#Lambda) - M(#bar{#Lambda})"};
  ConfigurableAxis binLambdaMassRatio{"binLambdaMassRatio", {50, -0.05f, 4.95f}, "M(#bar{#Lambda}) / M(#Lambda)"};
  ConfigurableAxis binLambdaMassRelDiff{"binLambdaMassRelDiff", {200, -0.995f, 1.005f}, "(M(#Lambda) - M(#bar{#Lambda})) / M(#Lambda)"};

  Preslice<MyTracks> TracksPerCollision = aod::track::collisionId;
  Preslice<aod::V0Datas> V0sPerCollision = aod::v0data::collisionId;

  void init(InitContext& initContext)
  {
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

    AxisSpec pdgAxis = {binPDG, ""};
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
    AxisSpec ptTrackRelDiffAxis = {binPtRelDiff, "(#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part})/#it{p}_{T, jet}^{part}"};
    AxisSpec zRelDiffAxis = {binZRelDiff, "(#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part})/#it{p}_{T, jet}^{part}"};

    AxisSpec V0PtAxis = {binV0Pt, "#it{p}_{T}^{V0}"};
    AxisSpec V0EtaAxis = {binV0Eta, "#eta^{V0}"};
    AxisSpec V0PhiAxis = {binV0Phi, "#varphi^{V0}"};
    AxisSpec V0detPtAxis = {binV0Pt, "#it{p}_{T}^{V0, det}"};
    AxisSpec V0partPtAxis = {binV0Pt, "#it{p}_{T}^{V0, part}"};
    AxisSpec V0CtauAxis = {binV0Ctau, "c#tau (cm)"};
    AxisSpec V0RadiusAxis = {binV0Radius, "R (cm)"};
    AxisSpec V0CosPAAxis = {binV0CosPA, "cos(PA)"};
    AxisSpec V0DCApAxis = {binV0DCAp, "DCA pos (cm)"};
    AxisSpec V0DCAnAxis = {binV0DCAn, "DCA neg (cm)"};
    AxisSpec V0DCAdAxis = {binV0DCAd, "DCA daughters (cm^{2})"};

    AxisSpec K0SMassAxis = {binK0SMass, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec LambdaMassAxis = {binLambdaMass, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec LambdaMassDiffAxis = {binLambdaMassDiff, "M(#Lambda) - M(#bar{#Lambda})"};
    AxisSpec LambdaMassRatioAxis = {binLambdaMassRatio, "M(#bar{#Lambda}) / M(#Lambda)"};
    AxisSpec LambdaMassRelDiffAxis = {binLambdaMassRelDiff, "(M(#Lambda) - M(#bar{#Lambda})) / M(#Lambda)"};

    if (doprocessDataRun3 || doprocessDataV0Frag) {
      registry.add("data/nJetsnTracks", "nJetsnTracks; nJets; nTracks", HistType::kTH2D, {jetCount, trackCount});
      registry.add("data/collision/collisionVtxZ", "Collision vertex z (cm)", HistType::kTH1D, {binVtxZ});
      registry.add("data/tracks/trackPtEtaPhi", "trackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis});

      registry.add("data/jets/jetPtEtaPhi", "Jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {jetPtAxis, etaAxis, phiAxis});
      registry.add("data/jets/jetPtTrackPt", "Jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2D, {jetPtAxis, trackPtAxis});
      registry.add("data/jets/jetTrackPtEtaPhi", "Tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis});
      registry.add("data/jets/jetPtFrag", "Jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2D, {jetPtAxis, zAxis});
      registry.add("data/jets/jetPtTrackProj", "Jet #it{p}_{T}, #it{z}", HistType::kTH2D, {jetPtAxis, zAxis});
      registry.add("data/jets/jetPtXi", "Jet #it{p}_{T}, #xi", HistType::kTH2D, {jetPtAxis, xiAxis});
      registry.add("data/jets/jetPtTheta", "Jet #it{p}_{T}, #theta", HistType::kTH2D, {jetPtAxis, thetaAxis});
      registry.add("data/jets/jetPtXiTheta", "Jet #it{p}_{T}, #xi, #theta", HistType::kTH3D, {jetPtAxis, xiAxis, thetaAxis});
      registry.add("data/jets/jetPtZTheta", "Jet #it{p}_{T}, z, #theta", HistType::kTH3D, {jetPtAxis, zAxis, thetaAxis});
    }

    if (doprocessDataV0 || doprocessDataV0Frag) {
      registry.add("data/V0/nV0sEvent", "nV0sEvent", HistType::kTH1D, {v0Count});

      // Unidentified
      registry.add("data/V0/V0PtEtaPhi", "V0PtEtaPhi", HistType::kTH3D, {V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("data/V0/V0PtCtau", "V0PtCtau", HistType::kTHnSparseD, {V0PtAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("data/V0/V0PtMass", "V0PtMass", HistType::kTHnSparseD, {V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/V0/V0PtLambdaMasses", "V0PtLambdaMasses", HistType::kTHnSparseD, {V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/V0/V0PtRadiusCosPA", "V0PtRadiusCosPA", HistType::kTH3D, {V0PtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("data/V0/V0PtDCAposneg", "V0PtDCAposneg", HistType::kTH3D, {V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("data/V0/V0PtDCAd", "V0PtDCAd", HistType::kTH2D, {V0PtAxis, V0DCAdAxis});

      // Identified
      registry.add("data/V0/K0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTH3D, {V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("data/V0/K0SPtCtauMass", "K0SPtCtauMass", HistType::kTH3D, {V0partPtAxis, V0CtauAxis, K0SMassAxis});
      registry.add("data/V0/K0SPtRadiusCosPA", "K0SPtRadiusCosPA", HistType::kTH3D, {V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("data/V0/K0SPtDCAposneg", "K0SPtDCAposneg", HistType::kTH3D, {V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("data/V0/K0SPtDCAd", "K0SPtDCAd", HistType::kTH2D, {V0partPtAxis, V0DCAdAxis});

      registry.add("data/V0/LambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTH3D, {V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("data/V0/LambdaPtCtauMass", "LambdaPtCtauMass", HistType::kTH3D, {V0partPtAxis, V0CtauAxis, LambdaMassAxis});
      registry.add("data/V0/LambdaPtLambdaMasses", "LambdaPtLambdaMasses", HistType::kTHnSparseD, {V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/V0/LambdaPtRadiusCosPA", "LambdaPtRadiusCosPA", HistType::kTH3D, {V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("data/V0/LambdaPtDCAposneg", "LambdaPtDCAposneg", HistType::kTH3D, {V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("data/V0/LambdaPtDCAd", "LambdaPtDCAd", HistType::kTH2D, {V0partPtAxis, V0DCAdAxis});

      registry.add("data/V0/antiLambdaPtEtaPhi", "antiLambdaPtEtaPhi", HistType::kTH3D, {V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("data/V0/antiLambdaPtCtauMass", "antiLambdaPtCtauMass", HistType::kTH3D, {V0partPtAxis, V0CtauAxis, LambdaMassAxis});
      registry.add("data/V0/antiLambdaPtLambdaMasses", "antiLambdaPtLambdaMasses", HistType::kTHnSparseD, {V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/V0/antiLambdaPtRadiusCosPA", "antiLambdaPtRadiusCosPA", HistType::kTH3D, {V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("data/V0/antiLambdaPtDCAposneg", "antiLambdaPtDCAposneg", HistType::kTH3D, {V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("data/V0/antiLambdaPtDCAd", "antiLambdaPtDCAd", HistType::kTH2D, {V0partPtAxis, V0DCAdAxis});
    }

    if (doprocessDataV0Frag) {
      registry.add("data/jets/V0/jetCorrectedPtEtaPhi", "Jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {jetPtAxis, etaAxis, phiAxis});
      registry.add("data/jets/V0/jetPtnV0", "jetPtnV0", HistType::kTH2D, {jetPtAxis, v0Count});
      registry.add("data/jets/V0/jetCorrectedPtV0TrackProj", "jetCorrectedPtV0TrackProj", HistType::kTH2D, {jetPtAxis, zAxis});
      registry.add("data/jets/V0/jetPtV0TrackProj", "jetPtV0TrackProj", HistType::kTH2D, {jetPtAxis, zAxis});
      registry.add("data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda", "jetPtnV0nK0SnLambdanAntiLambda", HistType::kTHnSparseD, {jetPtAxis, v0Count, v0Count, v0Count, v0Count});

      registry.add("data/jets/V0/jetPtV0PtEtaPhi", "jetPtV0PtEtaPhi", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("data/jets/V0/jetPtV0PtCtau", "jetPtV0PtCtau", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtV0PtMass", "jetPtV0PtMass", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtV0PtLambdaMasses", "jetPtV0PtLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtV0PtRadiusCosPA", "jetPtV0PtRadiusCosPA", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtV0PtDCAposneg", "jetPtV0PtDCAposneg", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("data/jets/V0/jetPtV0PtDCAd", "jetPtV0PtDCAd", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0DCAdAxis});

      registry.add("data/jets/V0/jetPtV0TrackProjCtau", "jetPtV0TrackProjCtau", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjMass", "jetPtV0TrackProjMass", HistType::kTHnSparseD, {jetPtAxis, zAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjLambdaMasses", "jetPtV0TrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjRadiusCosPA", "jetPtV0TrackProjRadiusCosPA", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjDCAposneg", "jetPtV0TrackProjDCAposneg", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("data/jets/V0/jetPtV0TrackProjDCAd", "jetPtV0TrackProjDCAd", HistType::kTH3D, {jetPtAxis, zAxis, V0DCAdAxis});

      // Identified
      registry.add("data/jets/V0/jetPtnLambda", "jetPtnLambda", HistType::kTH2D, {jetPtAxis, trackCount});
      registry.add("data/jets/V0/jetPtLambdaPtCtau", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtLambdaPtMass", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, mass", HistType::kTH3D, {jetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtLambdaMasses", "jetPtLambdaPtLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtLambdaPtRadius", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, radius", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtLambdaPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtLambdaTrackProjCtau", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjMass", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjLambdaMasses", "jetPtLambdaTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjRadius", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtnAntiLambda", "jetPtnAntiLambda", HistType::kTH2D, {jetPtAxis, trackCount});
      registry.add("data/jets/V0/jetPtAntiLambdaPtCtau", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtMass", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, mass", HistType::kTH3D, {jetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtLambdaMasses", "jetPtAntiLambdaPtLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtRadius", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, radius", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjCtau", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjMass", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjLambdaMasses", "jetPtAntiLambdaTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjRadius", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtnK0S", "jetPtnK0S", HistType::kTH2D, {jetPtAxis, trackCount});
      registry.add("data/jets/V0/jetPtK0SPtCtau", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, c#tau", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtK0SPtMass", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, mass", HistType::kTH3D, {jetPtAxis, V0PtAxis, K0SMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtRadius", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, radius", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtK0SPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, cosPA", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, DCA daughters", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtK0STrackProjCtau", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjMass", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, K0SMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjRadius", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0DCApAxis, V0DCAnAxis});
    }

    if (doprocessMcP) {
      registry.add("particle-level/nJetsnTracks", "nJetsnTracks; nJets; nTracks", HistType::kTH2D, {jetCount, trackCount});
      registry.add("particle-level/collision/partCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1D, {binVtxZ});
      registry.add("particle-level/tracks/partTrackPtEtaPhi", "partTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis});
      registry.add("particle-level/jets/partJetPtEtaPhi", "Particle level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis});
      registry.add("particle-level/jets/partJetPtTrackPt", "Particle level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2D, {partJetPtAxis, trackPtAxis});
      registry.add("particle-level/jets/partJetTrackPtEtaPhi", "Particle level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3D, {trackPtAxis, partEtaAxis, partPhiAxis});
      registry.add("particle-level/jets/partJetPtFrag", "Particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("particle-level/jets/partJetPtTrackProj", "Particle level jet #it{p}_{T}, #it{z}", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("particle-level/jets/partJetPtXi", "Particle level jet #it{p}_{T}, #xi", HistType::kTH2D, {partJetPtAxis, partXiAxis});
      registry.add("particle-level/jets/partJetPtTheta", "Particle level jet #it{p}_{T}, #theta", HistType::kTH2D, {partJetPtAxis, partThetaAxis});
      registry.add("particle-level/jets/partJetPtXiTheta", "Particle level jet #it{p}_{T}, #xi, #theta", HistType::kTH3D, {partJetPtAxis, partXiAxis, partThetaAxis});
      registry.add("particle-level/jets/partJetPtZTheta", "Particle level jet #it{p}_{T}, z, #theta", HistType::kTH3D, {partJetPtAxis, partZAxis, partThetaAxis});
    }

    if (doprocessMcD) {
      registry.add("detector-level/nJetsnTracks", "nJetsnTracks; nJets; nTracks", HistType::kTH2D, {jetCount, trackCount});
      registry.add("detector-level/collision/detCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1D, {binVtxZ});
      registry.add("detector-level/tracks/detTrackPtEtaPhi", "detTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis});
      registry.add("detector-level/jets/detJetPtEtaPhi", "Detector level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis});
      registry.add("detector-level/jets/detJetPtTrackPt", "Detector level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2D, {detJetPtAxis, trackPtAxis});
      registry.add("detector-level/jets/detJetTrackPtEtaPhi", "Detector level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3D, {trackPtAxis, detEtaAxis, detPhiAxis});
      registry.add("detector-level/jets/detJetPtFrag", "Detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("detector-level/jets/detJetPtTrackProj", "Detector level jet #it{p}_{T}, #it{z}", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("detector-level/jets/detJetPtXi", "Detector level jet #it{p}_{T}, #xi", HistType::kTH2D, {detJetPtAxis, detXiAxis});
      registry.add("detector-level/jets/detJetPtTheta", "Detector level jet #it{p}_{T}, #theta", HistType::kTH2D, {detJetPtAxis, detThetaAxis});
      registry.add("detector-level/jets/detJetPtXiTheta", "Detector level jet #it{p}_{T}, #xi, #theta", HistType::kTH3D, {detJetPtAxis, detXiAxis, detThetaAxis});
      registry.add("detector-level/jets/detJetPtZTheta", "Detector level jet #it{p}_{T}, z, #theta", HistType::kTH3D, {detJetPtAxis, detZAxis, detThetaAxis});
    }

    if (doprocessMcMatched) {
      registry.add("matching/collision/matchCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1D, {binVtxZ});
      registry.add("matching/tracks/matchDetTrackPtEtaPhi", "matchDetTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis});
      registry.add("matching/tracks/matchPartTrackPtEtaPhi", "matchPartTrackPtEtaPhi", HistType::kTH3D, {trackPtAxis, etaAxis, phiAxis});
      registry.add("matching/tracks/matchDetTrackPtPartTrackPt", "matchDetTrackPtPartTrackPt", HistType::kTH2D, {trackPtAxis, trackPtAxis});
      registry.add("matching/tracks/matchDetTrackEtaPartTrackEta", "matchDetTrackEtaPartTrackEta", HistType::kTH2D, {etaAxis, etaAxis});
      registry.add("matching/tracks/matchDetTrackPhiPartTrackPhi", "matchDetTrackPhiPartTrackPhi", HistType::kTH2D, {phiAxis, phiAxis});
      registry.add("matching/tracks/trackResolutionPt", "trackResolutionPt", HistType::kTH2D, {trackPtAxis, ptDiffAxis});
      registry.add("matching/tracks/trackResolutionEta", "trackResolutionEta", HistType::kTH2D, {etaAxis, etaDiffAxis});
      registry.add("matching/tracks/trackResolutionPhi", "trackResolutionPhi", HistType::kTH2D, {phiAxis, phiDiffAxis});
      // Detector level jets with a match
      registry.add("matching/jets/matchDetJetPtEtaPhi", "Matched detector level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis});
      registry.add("matching/jets/matchDetJetPtTrackPt", "Matched detector level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2D, {detJetPtAxis, trackPtAxis});
      registry.add("matching/jets/matchDetJetTrackPtEtaPhi", "Matched detector level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3D, {trackPtAxis, detEtaAxis, detPhiAxis});
      registry.add("matching/jets/matchDetJetPtFrag", "Matched detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/matchDetJetPtTrackProj", "Matched detector level jet #it{p}_{T}, #it{z}", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/matchDetJetPtXi", "Matched detector level jet #it{p}_{T}, #xi", HistType::kTH2D, {detJetPtAxis, detXiAxis});
      registry.add("matching/jets/matchDetJetPtTheta", "Matched detector level jet #it{p}_{T}, #theta", HistType::kTH2D, {detJetPtAxis, detThetaAxis});
      registry.add("matching/jets/matchDetJetPtXiTheta", "Matched detector level jet #it{p}_{T}, #xi, #theta", HistType::kTH3D, {detJetPtAxis, detXiAxis, detThetaAxis});
      registry.add("matching/jets/matchDetJetPtZTheta", "Matched detector level jet #it{p}_{T}, z, #theta", HistType::kTH3D, {detJetPtAxis, detZAxis, detThetaAxis});
      // Particle level jets with a match
      registry.add("matching/jets/matchPartJetPtEtaPhi", "Matched particle level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis});
      registry.add("matching/jets/matchPartJetPtTrackPt", "Matched particle level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2D, {partJetPtAxis, trackPtAxis});
      registry.add("matching/jets/matchPartJetTrackPtEtaPhi", "Matched particle level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3D, {trackPtAxis, partEtaAxis, partPhiAxis});
      registry.add("matching/jets/matchPartJetPtFrag", "Matched particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/matchPartJetPtTrackProj", "Matched particle level jet #it{p}_{T}, #it{z}", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/matchPartJetPtXi", "Matched particle level jet #it{p}_{T}, #xi", HistType::kTH2D, {partJetPtAxis, partXiAxis});
      registry.add("matching/jets/matchPartJetPtTheta", "Matched particle level jet #it{p}_{T}, #theta", HistType::kTH2D, {partJetPtAxis, partThetaAxis});
      registry.add("matching/jets/matchPartJetPtXiTheta", "Matched particle level jet #it{p}_{T}, #xi, #theta", HistType::kTH3D, {partJetPtAxis, partXiAxis, partThetaAxis});
      registry.add("matching/jets/matchPartJetPtZTheta", "Matched particle level jet #it{p}_{T}, z, #theta", HistType::kTH3D, {partJetPtAxis, partZAxis, partThetaAxis});
      // Combined information of matched jets
      registry.add("matching/jets/matchDetJetPtPartJetPt", "matchDetJetPtPartJetPt", HistType::kTH2D, {detJetPtAxis, partJetPtAxis});
      registry.add("matching/jets/matchPartJetPtDetJetEtaPartJetEta", "matchPartJetPtDetJetEtaPartJetEta", HistType::kTH3D, {partJetPtAxis, detEtaAxis, partEtaAxis});
      registry.add("matching/jets/matchPartJetPtDetJetPhiPartJetPhi", "matchPartJetPtDetJetPhiPartJetPhi", HistType::kTH3D, {partJetPtAxis, detPhiAxis, partPhiAxis});
      registry.add("matching/jets/matchPartJetPtResolutionPt", "#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part}", HistType::kTH2D, {partJetPtAxis, ptDiffAxis});
      registry.add("matching/jets/matchPartJetPtRelDiffPt", "#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part}", HistType::kTH2D, {partJetPtAxis, ptJetRelDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionEta", "#eta^{jet, det} - #eta^{jet, part}", HistType::kTH3D, {partJetPtAxis, partEtaAxis, etaDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionPhi", "#phi^{jet, det} - #phi^{jet, part}", HistType::kTH3D, {partJetPtAxis, partPhiAxis, phiDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionChargeFrag", "Resolution #it{p}_{T}^{tr} / #it{p}_{T}^{jet}", HistType::kTH3D, {partJetPtAxis, partZAxis, zDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionTrackPt", "Resolution #it{p}_{T}^{track}", HistType::kTH3D, {partJetPtAxis, trackPtAxis, ptTrackDiffAxis});
      registry.add("matching/jets/matching/jets/matchPartJetPtRelDiffTrackPt", "Rel. diff #it{p}_{T}^{track}", HistType::kTHnSparseD, {partJetPtAxis, ptRatioAxis, trackPtAxis, ptTrackRelDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionTrackProj", "Resolution #it{p}^{proj} / #it{p}^{jet}", HistType::kTH3D, {partJetPtAxis, partZAxis, zDiffAxis});
      registry.add("matching/jets/matchPartJetPtRelDiffTrackProj", "Rel. diff #it{p}^{proj} / #it{p}^{jet}", HistType::kTHnSparseD, {partJetPtAxis, ptRatioAxis, partZAxis, zRelDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionXi", "Resolution ln(1/#it{z})", HistType::kTH3D, {partJetPtAxis, partXiAxis, xiDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionTheta", "Resolution #theta", HistType::kTH3D, {partJetPtAxis, partThetaAxis, thetaDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionXiResolutionTheta", "Resolution #xi, #theta", HistType::kTHnSparseD, {partJetPtAxis, partXiAxis, xiDiffAxis, partThetaAxis, thetaDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionZResolutionTheta", "Resolution #it{z}, #theta", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, zDiffAxis, partThetaAxis, thetaDiffAxis});
      registry.add("matching/jets/matchPartJetPtEtaPhiMatchDist", "matchJetMatchDist", HistType::kTHnSparseD, {partJetPtAxis, partEtaAxis, partPhiAxis, matchDistAxis});
      registry.add("matching/jets/matchPartJetPtEnergyScale", "jetEnergyScale", HistType::kTH2D, {partJetPtAxis, ptRatioAxis});
      // QA histograms for fakes, misses
      registry.add("matching/jets/fakeDetJetPtEtaPhi", "Fakes", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis});
      registry.add("matching/jets/missPartJetPtEtaPhi", "Misses", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis});
      // Response matrix, fakes, misses
      registry.add("matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
      registry.add("matching/jets/fakeDetJetPtTrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/missPartJetPtTrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis});

      registry.add("matching/jets/matchDetJetPtXiPartJetPtXi", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detXiAxis, partJetPtAxis, partXiAxis});
      registry.add("matching/jets/fakeDetJetPtXi", "Fakes", HistType::kTH2D, {detJetPtAxis, detXiAxis});
      registry.add("matching/jets/missPartJetPtXi", "Misses", HistType::kTH2D, {partJetPtAxis, partXiAxis});

      registry.add("matching/jets/matchDetJetPtFragPartJetPtFrag", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
      registry.add("matching/jets/fakeDetJetPtFrag", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/missPartJetPtFrag", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis});

      registry.add("matching/jets/matchDetJetPtThetaPartJetPtTheta", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detThetaAxis, partJetPtAxis, partThetaAxis});
      registry.add("matching/jets/fakeDetJetPtTheta", "Fakes", HistType::kTH2D, {detJetPtAxis, detThetaAxis});
      registry.add("matching/jets/missPartJetPtTheta", "Misses", HistType::kTH2D, {partJetPtAxis, partThetaAxis});

      // TODO: Maybe need different axes here to have less granularity. Histogram may become too large.
      registry.add("matching/jets/matchDetJetPtXiThetaPartJetPtXiTheta", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detXiAxis, detThetaAxis, partJetPtAxis, partXiAxis, partThetaAxis});
      registry.add("matching/jets/fakeDetJetPtXiTheta", "Fakes", HistType::kTH3D, {detJetPtAxis, detXiAxis, detThetaAxis});
      registry.add("matching/jets/missPartJetPtXiTheta", "Misses", HistType::kTH3D, {partJetPtAxis, partXiAxis, partThetaAxis});

      registry.add("matching/jets/matchDetJetPtZThetaPartJetPtZTheta", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, detThetaAxis, partJetPtAxis, partZAxis, partThetaAxis});
      registry.add("matching/jets/fakeDetJetPtZTheta", "Fakes", HistType::kTH3D, {detJetPtAxis, detZAxis, detThetaAxis});
      registry.add("matching/jets/missPartJetPtZTheta", "Misses", HistType::kTH3D, {partJetPtAxis, partZAxis, partThetaAxis});
    }

    if (doprocessMcV0) {
      registry.add("matching/V0/V0PartPtDetPt", "V0PartPtDetPt", HistType::kTH2D, {V0partPtAxis, V0detPtAxis});

      registry.add("matching/V0/K0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTH3D, {V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/V0/K0SPtCtauMass", "K0SPtCtauMass", HistType::kTH3D, {V0partPtAxis, V0CtauAxis, K0SMassAxis});
      registry.add("matching/V0/K0SPtRadiusCosPA", "K0SPtRadiusCosPA", HistType::kTH3D, {V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/V0/K0SPtDCAposneg", "K0SPtDCAposneg", HistType::kTH3D, {V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/V0/K0SPtDCAd", "K0SPtDCAd", HistType::kTH2D, {V0partPtAxis, V0DCAdAxis});

      registry.add("matching/V0/LambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTH3D, {V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/V0/LambdaPtCtauMass", "LambdaPtCtauMass", HistType::kTH3D, {V0partPtAxis, V0CtauAxis, LambdaMassAxis});
      registry.add("matching/V0/LambdaPtRadiusCosPA", "LambdaPtRadiusCosPA", HistType::kTH3D, {V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/V0/LambdaPtDCAposneg", "LambdaPtDCAposneg", HistType::kTH3D, {V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/V0/LambdaPtDCAd", "LambdaPtDCAd", HistType::kTH2D, {V0partPtAxis, V0DCAdAxis});

      registry.add("matching/V0/antiLambdaPtEtaPhi", "antiLambdaPtEtaPhi", HistType::kTH3D, {V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/V0/antiLambdaPtCtauMass", "antiLambdaPtCtauMass", HistType::kTH3D, {V0partPtAxis, V0CtauAxis, LambdaMassAxis});
      registry.add("matching/V0/antiLambdaPtRadiusCosPA", "antiLambdaPtRadiusCosPA", HistType::kTH3D, {V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/V0/antiLambdaPtDCAposneg", "antiLambdaPtDCAposneg", HistType::kTH3D, {V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/V0/antiLambdaPtDCAd", "antiLambdaPtDCAd", HistType::kTH2D, {V0partPtAxis, V0DCAdAxis});
    }
  } // init

  // TODO: Can we move most/all of this stuff into a filter?
  template <typename V0Type>
  bool IsV0Candidate(V0Type v0)
  {
    if (v0.eta() < v0EtaMin || v0.eta() > v0EtaMax) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < dcanegtopv) {
      return false;
    }
    if (TMath::Abs(v0.dcapostopv()) < dcapostopv) {
      return false;
    }
    if (v0.dcaV0daughters() > dcav0dau) {
      return false;
    }
    if (v0.v0radius() < v0radius) {
      return false;
    }
    if (v0.v0cosPA() < v0cospa) {
      return false;
    }
    return true;
  }
  template <typename CollisionType, typename V0Type>
  bool IsK0SCandidate(CollisionType const& collision, V0Type v0)
  {
    if (!IsV0Candidate(v0)) {
      return false;
    }
    double ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    if (ctauK0s > lifetimecutK0S) {
      return false;
    }
    bool k0sMassCondition = (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < k0sMassAccWindow);
    bool lambdaMassCondition = (TMath::Abs(v0.mLambda() - o2::constants::physics::MassLambda0) < lambdaMassRejWindow);
    bool antilambdaMassCondition = (TMath::Abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0Bar) < antilambdaMassRejWindow);
    if (!k0sMassCondition) {
      return false;
    }
    if (lambdaMassCondition || antilambdaMassCondition) {
      return false;
    }
    return true;
  }
  template <typename CollisionType, typename V0Type>
  bool IsLambdaCandidate(CollisionType const& collision, V0Type v0)
  {
    if (!IsV0Candidate(v0)) {
      return false;
    }
    double ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (ctauLambda > lifetimecutLambda) {
      return false;
    }
    bool k0sMassCondition = (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < k0sMassRejWindow);
    bool lambdaMassCondition = (TMath::Abs(v0.mLambda() - o2::constants::physics::MassLambda0) < lambdaMassAccWindow);
    bool antilambdaMassCondition = (TMath::Abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0Bar) < antilambdaMassRejWindow);
    if (!lambdaMassCondition) {
      return false;
    }
    if (k0sMassCondition || antilambdaMassCondition) {
      return false;
    }
    return true;
  }
  template <typename CollisionType, typename V0Type>
  bool IsAntiLambdaCandidate(CollisionType const& collision, V0Type v0)
  {
    if (!IsV0Candidate(v0)) {
      return false;
    }
    double ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    if (ctauAntiLambda > lifetimecutLambda) {
      return false;
    }
    bool k0sMassCondition = (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < k0sMassRejWindow);
    bool lambdaMassCondition = (TMath::Abs(v0.mLambda() - o2::constants::physics::MassLambda0) < lambdaMassRejWindow);
    bool antilambdaMassCondition = (TMath::Abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0Bar) < antilambdaMassAccWindow);
    if (!antilambdaMassCondition) {
      return false;
    }
    if (k0sMassCondition || lambdaMassCondition) {
      return false;
    }
    return true;
  }

  template <typename Jet, typename Constituent>
  double ChargeFrag(Jet const& jet, Constituent const& constituent)
  {
    double chargeFrag = -1.;
    chargeFrag = constituent.pt() / jet.pt();
    return chargeFrag;
  }

  template <typename Jet, typename Constituent>
  double Theta(Jet const& jet, Constituent const& constituent)
  {
    double theta = -1.;
    theta = jetutilities::deltaR(jet, constituent);
    return theta;
  }

  template <typename Jet, typename Constituent>
  double TrackProj(Jet const& jet, Constituent const& constituent)
  {
    double trackProj = -1.;
    trackProj = constituent.px() * jet.px() + constituent.py() * jet.py() + constituent.pz() * jet.pz();
    trackProj /= (jet.p() * jet.p());
    return trackProj;
  }

  template <typename Jet, typename Constituent>
  double Xi(Jet const& jet, Constituent const& constituent)
  {
    double xi = -1., trackProj = -1.;
    trackProj = TrackProj(jet, constituent);
    if (trackProj > 0) {
      xi = TMath::Log(1. / trackProj);
    }
    return xi;
  }

  template <typename collisionType, typename v0Type, typename trackType, typename particleType>
  void fillMcV0Histograms(collisionType const& collision, v0Type const& v0, trackType const& tracks, particleType const& particles, double weight = 1.)
  {
    auto negTrack = v0.template negTrack_as<trackType>();
    auto posTrack = v0.template posTrack_as<trackType>();
    if (!negTrack.has_mcParticle() || !posTrack.has_mcParticle()) {
      return;
    }
    auto mcNegTrack = negTrack.template mcParticle_as<particleType>();
    auto mcPosTrack = posTrack.template mcParticle_as<particleType>();
    if (!mcNegTrack.has_mothers() || !mcPosTrack.has_mothers()) {
      return;
    }
    double ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    // Can tracks have multiple mothers?
    for (auto& particleMotherOfNeg : mcNegTrack.template mothers_as<aod::McParticles>()) {
      for (auto& particleMotherOfPos : mcPosTrack.template mothers_as<aod::McParticles>()) {
        if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos) {
          double ptPartV0 = particleMotherOfNeg.pt();
          int pdg = particleMotherOfNeg.pdgCode();
          registry.fill(HIST("matching/V0/V0PartPtDetPt"), ptPartV0, v0.pt());

          if (pdg == 310) { // K0S
            registry.fill(HIST("matching/V0/K0SPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
            registry.fill(HIST("matching/V0/K0SPtCtauMass"), ptPartV0, ctauK0s, v0.mK0Short(), weight);
            registry.fill(HIST("matching/V0/K0SPtRadiusCosPA"), ptPartV0, v0.v0radius(), v0.v0cosPA(), weight);
            registry.fill(HIST("matching/V0/K0SPtDCAposneg"), ptPartV0, v0.dcapostopv(), v0.dcanegtopv(), weight);
            registry.fill(HIST("matching/V0/K0SPtDCAd"), ptPartV0, v0.dcaV0daughters(), weight);
          } else if (pdg == 3122) { // Lambda
            registry.fill(HIST("matching/V0/LambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
            registry.fill(HIST("matching/V0/LambdaPtCtauMass"), ptPartV0, ctauLambda, v0.mLambda(), weight);
            registry.fill(HIST("matching/V0/LambdaPtRadiusCosPA"), ptPartV0, v0.v0radius(), v0.v0cosPA(), weight);
            registry.fill(HIST("matching/V0/LambdaPtDCAposneg"), ptPartV0, v0.dcapostopv(), v0.dcanegtopv(), weight);
            registry.fill(HIST("matching/V0/LambdaPtDCAd"), ptPartV0, v0.dcaV0daughters(), weight);
          } else if (pdg == -3122) { // AntiLambda
            registry.fill(HIST("matching/V0/antiLambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
            registry.fill(HIST("matching/V0/antiLambdaPtCtauMass"), ptPartV0, ctauAntiLambda, v0.mAntiLambda(), weight);
            registry.fill(HIST("matching/V0/antiLambdaPtRadiusCosPA"), ptPartV0, v0.v0radius(), v0.v0cosPA(), weight);
            registry.fill(HIST("matching/V0/antiLambdaPtDCAposneg"), ptPartV0, v0.dcapostopv(), v0.dcanegtopv(), weight);
            registry.fill(HIST("matching/V0/antiLambdaPtDCAd"), ptPartV0, v0.dcaV0daughters(), weight);
          }
        } // if mothers match
      }   // for mothers of pos
    }     // for mothers of neg
  }

  template <typename T>
  void fillDataRun3Histograms(T const& jet)
  {
    registry.fill(HIST("data/jets/jetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
    for (const auto& track : jet.template tracks_as<JetTracks>()) {
      double chargeFrag = -1., trackProj = -1., xi = -1.;
      double theta = -1.;
      chargeFrag = ChargeFrag(jet, track);
      trackProj = TrackProj(jet, track);
      theta = Theta(jet, track);
      xi = Xi(jet, track);

      registry.fill(HIST("data/jets/jetPtTrackPt"), jet.pt(), track.pt());
      registry.fill(HIST("data/jets/jetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      registry.fill(HIST("data/jets/jetPtFrag"), jet.pt(), chargeFrag);
      registry.fill(HIST("data/jets/jetPtTrackProj"), jet.pt(), trackProj);
      registry.fill(HIST("data/jets/jetPtXi"), jet.pt(), xi);
      registry.fill(HIST("data/jets/jetPtTheta"), jet.pt(), theta);
      registry.fill(HIST("data/jets/jetPtXiTheta"), jet.pt(), xi, theta);
      registry.fill(HIST("data/jets/jetPtZTheta"), jet.pt(), trackProj, theta);
    }
  }

  template <typename CollisionType, typename V0Type, typename TrackType>
  void fillDataV0Histograms(CollisionType collision, V0Type V0s, TrackType tracks)
  {
    for (const auto& v0 : V0s) {
      double ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      double ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
      double ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

      double massDiff = v0.mLambda() - v0.mAntiLambda();
      double massRatio = v0.mAntiLambda() / v0.mLambda();
      double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

      registry.fill(HIST("data/V0/V0PtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
      registry.fill(HIST("data/V0/V0PtCtau"), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda);
      registry.fill(HIST("data/V0/V0PtMass"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/V0/V0PtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/V0/V0PtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
      registry.fill(HIST("data/V0/V0PtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
      registry.fill(HIST("data/V0/V0PtDCAd"), v0.pt(), v0.dcaV0daughters());

      if (IsLambdaCandidate(collision, v0)) {
        registry.fill(HIST("data/V0/LambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/V0/LambdaPtCtauMass"), v0.pt(), ctauLambda, v0.mLambda());
        registry.fill(HIST("data/V0/LambdaPtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff);
        registry.fill(HIST("data/V0/LambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/V0/LambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/V0/LambdaPtDCAd"), v0.pt(), v0.dcaV0daughters());
      }
      if (IsAntiLambdaCandidate(collision, v0)) {
        registry.fill(HIST("data/V0/antiLambdaPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/V0/antiLambdaPtCtauMass"), v0.pt(), ctauAntiLambda, v0.mAntiLambda());
        registry.fill(HIST("data/V0/antiLambdaPtLambdaMasses"), v0.pt(), massDiff, massRatio, massRelDiff);
        registry.fill(HIST("data/V0/antiLambdaPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/V0/antiLambdaPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/V0/antiLambdaPtDCAd"), v0.pt(), v0.dcaV0daughters());
      }
      if (IsK0SCandidate(collision, v0)) {
        registry.fill(HIST("data/V0/K0SPtEtaPhi"), v0.pt(), v0.eta(), v0.phi());
        registry.fill(HIST("data/V0/K0SPtCtauMass"), v0.pt(), ctauK0s, v0.mK0Short());
        registry.fill(HIST("data/V0/K0SPtRadiusCosPA"), v0.pt(), v0.v0radius(), v0.v0cosPA());
        registry.fill(HIST("data/V0/K0SPtDCAposneg"), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
        registry.fill(HIST("data/V0/K0SPtDCAd"), v0.pt(), v0.dcaV0daughters());
      }
    } // for v0
  }

  template <typename CollisionType, typename JetType, typename V0Type>
  void fillDataV0FragHistograms(CollisionType collision, JetType jet, V0Type v0)
  {
    double trackProj = TrackProj(jet, v0);
    double ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

    double massDiff = v0.mLambda() - v0.mAntiLambda();
    double massRatio = v0.mAntiLambda() / v0.mLambda();
    double massRelDiff = (v0.mLambda() - v0.mAntiLambda()) / v0.mLambda();

    registry.fill(HIST("data/jets/V0/jetPtV0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi());
    registry.fill(HIST("data/jets/V0/jetPtV0PtCtau"), jet.pt(), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda);
    registry.fill(HIST("data/jets/V0/jetPtV0PtMass"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
    registry.fill(HIST("data/jets/V0/jetPtV0PtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff);
    registry.fill(HIST("data/jets/V0/jetPtV0PtRadiusCosPA"), jet.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA());
    registry.fill(HIST("data/jets/V0/jetPtV0PtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());
    registry.fill(HIST("data/jets/V0/jetPtV0PtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());

    registry.fill(HIST("data/jets/V0/jetPtV0TrackProj"), jet.pt(), trackProj);
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjCtau"), jet.pt(), trackProj, ctauK0s, ctauLambda, ctauAntiLambda);
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjMass"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff);
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjRadiusCosPA"), jet.pt(), trackProj, v0.v0radius(), v0.v0cosPA());
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
    registry.fill(HIST("data/jets/V0/jetPtV0TrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());

    if (IsK0SCandidate(collision, v0)) {
      registry.fill(HIST("data/jets/V0/jetPtK0SPtCtau"), jet.pt(), v0.pt(), ctauK0s);
      registry.fill(HIST("data/jets/V0/jetPtK0SPtMass"), jet.pt(), v0.pt(), v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());

      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjCtau"), jet.pt(), trackProj, ctauK0s);
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjMass"), jet.pt(), trackProj, v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjRadius"), jet.pt(), trackProj, v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjCosPA"), jet.pt(), trackProj, v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
    }
    if (IsLambdaCandidate(collision, v0)) {
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtCtau"), jet.pt(), v0.pt(), ctauLambda);
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtMass"), jet.pt(), v0.pt(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());

      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjCtau"), jet.pt(), trackProj, ctauLambda);
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjMass"), jet.pt(), trackProj, v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjRadius"), jet.pt(), trackProj, v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjCosPA"), jet.pt(), trackProj, v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
    }
    if (IsAntiLambdaCandidate(collision, v0)) {
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtCtau"), jet.pt(), v0.pt(), ctauAntiLambda);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtMass"), jet.pt(), v0.pt(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());

      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjCtau"), jet.pt(), trackProj, ctauAntiLambda);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjMass"), jet.pt(), trackProj, v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjRadius"), jet.pt(), trackProj, v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjCosPA"), jet.pt(), trackProj, v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
    }
  }

  template <typename DetJet, typename PartJet>
  void fillMatchingHistogramsJet(DetJet const& detJet, PartJet const& partJet, double weight = 1.)
  {
    double deltaEta = detJet.eta() - partJet.eta();
    double deltaPhi = RecoDecay::constrainAngle(detJet.phi() - partJet.phi(), -M_PI);
    double dR = jetutilities::deltaR(detJet, partJet);

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

  template <typename DetJet, typename PartJet, typename Track, typename Particle>
  void fillMatchingHistogramsConstituent(DetJet const& detJet, PartJet const& partJet, Track const& track, Particle const& particle, double weight = 1.)
  {
    double detChargeFrag = -1., detTrackProj = -1., detTheta = -1., detXi = -1.;
    double partChargeFrag = -1., partTrackProj = -1., partTheta = -1., partXi = -1.;

    detChargeFrag = ChargeFrag(detJet, track);
    detTrackProj = TrackProj(detJet, track);
    detTheta = Theta(detJet, track);
    detXi = Xi(detJet, track);

    partChargeFrag = ChargeFrag(partJet, particle);
    partTrackProj = TrackProj(partJet, particle);
    partTheta = Theta(partJet, particle);
    partXi = Xi(partJet, particle);

    // Detector level
    registry.fill(HIST("matching/jets/matchDetJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi(), weight);
    registry.fill(HIST("matching/jets/matchDetJetPtTrackPt"), detJet.pt(), track.pt(), weight);
    registry.fill(HIST("matching/jets/matchDetJetPtFrag"), detJet.pt(), detChargeFrag, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtTrackProj"), detJet.pt(), detTrackProj, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtXi"), detJet.pt(), detXi, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtTheta"), detJet.pt(), detTheta, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtXiTheta"), detJet.pt(), detXi, detTheta, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtZTheta"), detJet.pt(), detTrackProj, detTheta, weight);

    // Particle level
    registry.fill(HIST("matching/jets/matchPartJetTrackPtEtaPhi"), particle.pt(), particle.eta(), particle.phi(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtTrackPt"), partJet.pt(), particle.pt(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtFrag"), partJet.pt(), partChargeFrag, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtTrackProj"), partJet.pt(), partTrackProj, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtXi"), partJet.pt(), partXi, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtTheta"), partJet.pt(), partTheta, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtXiTheta"), partJet.pt(), partXi, partTheta, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtZTheta"), partJet.pt(), partTrackProj, partTheta, weight);

    // Resolution
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionTrackPt"), partJet.pt(), particle.pt(), (particle.pt() - track.pt()), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionChargeFrag"), partJet.pt(), partChargeFrag, (detChargeFrag - partChargeFrag), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionTrackProj"), partJet.pt(), partTrackProj, (detTrackProj - partTrackProj), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionXi"), partJet.pt(), partXi, (detXi - partXi), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionTheta"), partJet.pt(), partTheta, (detTheta - partTheta), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionXiResolutionTheta"), partJet.pt(), partXi, (detXi - partXi), partTheta, (detTheta - partTheta), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionZResolutionTheta"), partJet.pt(), partTrackProj, (detTrackProj - partTrackProj), partTheta, (detTheta - partTheta), weight);

    // Relative difference
    registry.fill(HIST("matching/jets/matching/jets/matchPartJetPtRelDiffTrackPt"), partJet.pt(), detJet.pt() / partJet.pt(), particle.pt(), (track.pt() - particle.pt()) / particle.pt(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtRelDiffTrackProj"), partJet.pt(), detJet.pt() / partJet.pt(), partTrackProj, (detTrackProj - partTrackProj) / partTrackProj, weight);

    // Response
    registry.fill(HIST("matching/jets/matchDetJetPtFragPartJetPtFrag"), detJet.pt(), detChargeFrag, partJet.pt(), partChargeFrag, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj"), detJet.pt(), detTrackProj, partJet.pt(), partTrackProj, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtXiPartJetPtXi"), detJet.pt(), detXi, partJet.pt(), partXi, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtThetaPartJetPtTheta"), detJet.pt(), detTheta, partJet.pt(), partTheta, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtXiThetaPartJetPtXiTheta"), detJet.pt(), detXi, detTheta, partJet.pt(), partXi, partTheta, weight);
    registry.fill(HIST("matching/jets/matchDetJetPtZThetaPartJetPtZTheta"), detJet.pt(), detTrackProj, detTheta, partJet.pt(), partTrackProj, partTheta, weight);
  }

  template <typename Jet, typename Constituent>
  void fillMatchingFakeOrMiss(Jet const& jet, Constituent const& constituent, bool isFake, double weight = 1.)
  {
    double chargeFrag = -1., trackProj = -1., theta = -1., xi = -1.;
    chargeFrag = ChargeFrag(jet, constituent);
    trackProj = TrackProj(jet, constituent);
    theta = Theta(jet, constituent);
    xi = Xi(jet, constituent);

    if (isFake) {
      registry.fill(HIST("matching/jets/fakeDetJetPtFrag"), jet.pt(), chargeFrag, weight);
      registry.fill(HIST("matching/jets/fakeDetJetPtTrackProj"), jet.pt(), trackProj, weight);
      registry.fill(HIST("matching/jets/fakeDetJetPtXi"), jet.pt(), xi, weight);
      registry.fill(HIST("matching/jets/fakeDetJetPtTheta"), jet.pt(), theta, weight);
      registry.fill(HIST("matching/jets/fakeDetJetPtXiTheta"), jet.pt(), xi, theta, weight);
      registry.fill(HIST("matching/jets/fakeDetJetPtZTheta"), jet.pt(), trackProj, theta, weight);
    } else {
      registry.fill(HIST("matching/jets/missPartJetPtFrag"), jet.pt(), chargeFrag, weight);
      registry.fill(HIST("matching/jets/missPartJetPtTrackProj"), jet.pt(), trackProj, weight);
      registry.fill(HIST("matching/jets/missPartJetPtXi"), jet.pt(), xi, weight);
      registry.fill(HIST("matching/jets/missPartJetPtTheta"), jet.pt(), theta, weight);
      registry.fill(HIST("matching/jets/missPartJetPtXiTheta"), jet.pt(), xi, theta, weight);
      registry.fill(HIST("matching/jets/missPartJetPtZTheta"), jet.pt(), trackProj, theta, weight);
    }
  }

  template <typename Jet>
  void fillMCDHistograms(Jet const& jet, double weight = 1.)
  {
    registry.fill(HIST("detector-level/jets/detJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    for (const auto& track : jet.template tracks_as<JetTracks>()) {
      double chargeFrag = -1., trackProj = -1., theta = -1., xi = -1.;
      chargeFrag = ChargeFrag(jet, track);
      trackProj = TrackProj(jet, track);
      theta = Theta(jet, track);
      xi = Xi(jet, track);

      registry.fill(HIST("detector-level/jets/detJetPtTrackPt"), jet.pt(), track.pt(), weight);
      registry.fill(HIST("detector-level/jets/detJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi(), weight);
      registry.fill(HIST("detector-level/jets/detJetPtFrag"), jet.pt(), chargeFrag, weight);
      registry.fill(HIST("detector-level/jets/detJetPtTrackProj"), jet.pt(), trackProj, weight);
      registry.fill(HIST("detector-level/jets/detJetPtXi"), jet.pt(), xi, weight);
      registry.fill(HIST("detector-level/jets/detJetPtTheta"), jet.pt(), theta, weight);
      registry.fill(HIST("detector-level/jets/detJetPtXiTheta"), jet.pt(), xi, theta, weight);
      registry.fill(HIST("detector-level/jets/detJetPtZTheta"), jet.pt(), trackProj, theta, weight);
    }
  }

  template <typename Jet>
  void fillMCPHistograms(Jet const& jet, double weight = 1.)
  {
    registry.fill(HIST("particle-level/jets/partJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    for (const auto& track : jet.template tracks_as<JetParticles>()) {
      double chargeFrag = -1., trackProj = -1., theta = -1., xi = -1.;
      chargeFrag = ChargeFrag(jet, track);
      trackProj = TrackProj(jet, track);
      theta = Theta(jet, track);
      xi = Xi(jet, track);

      registry.fill(HIST("particle-level/jets/partJetPtTrackPt"), jet.pt(), track.pt(), weight);
      registry.fill(HIST("particle-level/jets/partJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi(), weight);
      registry.fill(HIST("particle-level/jets/partJetPtFrag"), jet.pt(), chargeFrag, weight);
      registry.fill(HIST("particle-level/jets/partJetPtTrackProj"), jet.pt(), trackProj, weight);
      registry.fill(HIST("particle-level/jets/partJetPtXi"), jet.pt(), xi, weight);
      registry.fill(HIST("particle-level/jets/partJetPtTheta"), jet.pt(), theta, weight);
      registry.fill(HIST("particle-level/jets/partJetPtXiTheta"), jet.pt(), xi, theta, weight);
      registry.fill(HIST("particle-level/jets/partJetPtZTheta"), jet.pt(), trackProj, theta, weight);
    }
  }

  void processDummy(JetTracks const& tracks) {}
  PROCESS_SWITCH(JetFragmentation, processDummy, "Dummy process function turned on by default", true);

  void processMcD(JetCollisionsMCD::iterator const& collision,
                  JetMcCollisions const& mcCollisions,
                  McDJets const& jets,
                  JetTracks const& tracks)
  {
    double nJets = 0, nTracks = 0;
    double weight = collision.mcCollision().weight();
    for (const auto& track : tracks) {
      if (track.pt() > 0.1) {
        nTracks++;
        registry.fill(HIST("detector-level/tracks/detTrackPtEtaPhi"), track.pt(), track.eta(), track.phi(), weight);
      }
    }
    for (const auto& jet : jets) {
      nJets++;
      fillMCDHistograms(jet, weight);
    }
    registry.fill(HIST("detector-level/nJetsnTracks"), nJets, nTracks, weight);
  }
  PROCESS_SWITCH(JetFragmentation, processMcD, "Monte Carlo detector level", false);

  void processMcP(JetMcCollision const& mcCollision, // Add some form of event selection?
                  McPJets const& jets,
                  JetParticles const& particles)
  {
    double nJets = 0, nTracks = 0;
    double weight = mcCollision.weight();
    for (const auto& particle : particles) {
      if (particle.pt() > 0.1) {
        nTracks++;
        registry.fill(HIST("particle-level/tracks/partTrackPtEtaPhi"), particle.pt(), particle.eta(), particle.phi(), weight);
      }
    }
    for (const auto& jet : jets) {
      nJets++;
      fillMCPHistograms(jet, weight);
    }
    registry.fill(HIST("particle-level/nJetsnTracks"), nJets, nTracks, weight);
  }

  PROCESS_SWITCH(JetFragmentation, processMcP, "Monte Carlo particle level", false);

  void processDataRun3(JetCollision const& collision,
                       soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                       JetTracks const& tracks)
  {
    double nJets = 0, nTracks = 0;
    for (const auto& track : tracks) {
      if (track.pt() > 0.1) {
        nTracks++;
        registry.fill(HIST("data/tracks/trackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      }
    }
    for (const auto& jet : jets) {
      nJets++;
      fillDataRun3Histograms(jet);
    }
    registry.fill(HIST("data/nJetsnTracks"), nJets, nTracks);
  }
  PROCESS_SWITCH(JetFragmentation, processDataRun3, "Run 3 Data", false);

  void processMcMatched(JetCollisionsMCD::iterator const& collision,
                        soa::Join<McDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& mcDetJets,
                        JetTracksMCD const& tracks,
                        JetMcCollisions const& mcCollisions,
                        McPJets const& mcPartJets,
                        JetParticles const& mcParticles)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    double weight = collision.mcCollision().weight();
    bool isFake = false;
    for (const auto& detJet : mcDetJets) {
      if (detJet.eta() < matchedDetJetEtaMin || detJet.eta() > matchedDetJetEtaMax) {
        continue; // TODO: should be done in filter
      }
      for (auto& partJet : detJet.template matchedJetGeo_as<McPJets>()) {
        fillMatchingHistogramsJet(detJet, partJet, weight);

        for (const auto& track : detJet.tracks_as<JetTracksMCD>()) {
          bool isTrackMatched = false;
          for (const auto& particle : partJet.tracks_as<JetParticles>()) {
            if (track.has_mcParticle() && particle.globalIndex() == track.template mcParticle_as<JetParticles>().globalIndex()) {
              isTrackMatched = true;
              fillMatchingHistogramsConstituent(detJet, partJet, track, particle, weight);
              break; // No need to inspect other particles
            }        // if track has mcParticle and particle is in matched jet
          }          // for particle in matched partJet
          if (!isTrackMatched) {
            isFake = true;
            fillMatchingFakeOrMiss(detJet, track, isFake, weight);
          } // if track is not matched
        }   // for detJet tracks
      }
      if (!detJet.has_matchedJetGeo()) {
        isFake = true;
        registry.fill(HIST("matching/jets/fakeDetJetPtEtaPhi"), detJet.pt(), detJet.eta(), detJet.phi(), weight);
        for (const auto& track : detJet.tracks_as<JetTracksMCD>()) {
          fillMatchingFakeOrMiss(detJet, track, isFake, weight);
        }
      } // if detJet does not have a match
    }   // for det jet
    // for (const auto& partJet : mcPartJets) {
    //   // We've already treated the matched jets in the previous loop
    //   if (!partJet.has_matchedJetGeo()) {
    //     isFake = false;
    //     registry.fill(HIST("matching/jets/missPartJetPtEtaPhi"), partJet.pt(), partJet.eta(), partJet.phi(), weight);
    //     for (const auto& particle : partJet.tracks_as<aod::McParticles>()) {
    //       fillMatchingFakeOrMiss(partJet, particle, isFake, weight);
    //     }
    //   }
    // }
  }
  PROCESS_SWITCH(JetFragmentation, processMcMatched, "Monte Carlo particle and detector level", false);

  void processMcV0(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
                   aod::McCollisions const& mcCollisions,
                   aod::V0Datas const& V0s,
                   soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels> const& tracks,
                   aod::McParticles const& mcParticles)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    double weight = collision.mcCollision().weight();
    for (const auto& v0 : V0s) {
      fillMcV0Histograms(collision, v0, tracks, mcParticles, weight);
    }
  }
  PROCESS_SWITCH(JetFragmentation, processMcV0, "Monte Carlo V0", false);

  void processDataV0(aod::Collision const& collision,
                     aod::V0Datas const& V0s,
                     MyTracks const& tracks)
  {
    registry.fill(HIST("data/V0/nV0sEvent"), V0s.size());
    fillDataV0Histograms(collision, V0s, tracks);
  }
  PROCESS_SWITCH(JetFragmentation, processDataV0, "Data V0", false);

  void processDataV0Frag(soa::Join<JetCollisions, aod::JCollisionPIs>::iterator const& jcoll,
                         soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                         JetTracks const& jtracks,
                         aod::Collisions const& collisions,
                         aod::V0Datas const& allV0s,
                         MyTracks const& allTracks)
  {
    // This is necessary, because jets are linked to JetCollisions, but V0s are linked to Collisions
    const auto& collision = jcoll.collision_as<aod::Collisions>();
    const auto& tracks = allTracks.sliceBy(TracksPerCollision, collision.globalIndex()); // Will use in future
    const auto& v0s = allV0s.sliceBy(V0sPerCollision, collision.globalIndex());

    int kNV0s = v0s.size();
    bool isV0Used[kNV0s];
    for (int i = 0; i < kNV0s; i++) {
      isV0Used[i] = false;
    }
    registry.fill(HIST("data/V0/nV0sEvent"), kNV0s);

    fillDataV0Histograms(collision, v0s, tracks);
    for (const auto& jet : jets) {
      if ((jet.eta() < v0EtaMin + jet.r() * 1e-2) || (jet.eta() > v0EtaMax - jet.r() * 1e-2)) {
        continue; // TODO: make filter?
      }
      fillDataRun3Histograms(jet);
      // fastjet::PseudoJet newjet(jet.px(), jet.py(), jet.pz(), jet.e()); // Jet with corrections from V0
      int iv0 = -1;
      int nV0inJet = 0, nLambdainJet = 0, nAntiLambdainJet = 0, nK0SinJet = 0;

      // Jets are pt-sorted, so we prioritise matching V0s with high pt jets
      // Correct jet momentum (currently only corrects for v0 in jet, not v0 outside jet, is this an issue?)
      for (const auto& v0 : v0s) {
        iv0++;
        if (isV0Used[iv0]) {
          continue;
        }
        double dR = jetutilities::deltaR(jet, v0);
        if (dR < jet.r() * 1e-2) {
          // fastjet::PseudoJet pjv0(v0.px(), v0.py(), v0.pz(), v0.e());
          // newjet += pjv0;
        }
      }
      // Loop over V0s and fill histograms
      iv0 = -1;
      for (const auto& v0 : v0s) {
        iv0++;
        if (isV0Used[iv0]) {
          continue;
        }
        double dR = jetutilities::deltaR(jet, v0);
        if (dR < jet.r() * 1e-2) {
          isV0Used[iv0] = true;
          nV0inJet++;
          fillDataV0FragHistograms(collision, jet, v0);
          if (IsK0SCandidate(collision, v0)) {
            nK0SinJet++;
          }
          if (IsLambdaCandidate(collision, v0)) {
            nLambdainJet++;
          }
          if (IsAntiLambdaCandidate(collision, v0)) {
            nAntiLambdainJet++;
          }
          // double newTrackProj = TrackProj(newjet, v0); // TODO: Does this work?
          // registry.fill(HIST("data/jets/V0/jetCorrectedPtV0TrackProj"), newjet.pt(), newTrackProj);
        }
      } // v0 loop
      registry.fill(HIST("data/jets/V0/jetPtnV0"), jet.pt(), nV0inJet);
      registry.fill(HIST("data/jets/V0/jetPtnLambda"), jet.pt(), nLambdainJet);
      registry.fill(HIST("data/jets/V0/jetPtnAntiLambda"), jet.pt(), nAntiLambdainJet);
      registry.fill(HIST("data/jets/V0/jetPtnK0S"), jet.pt(), nK0SinJet);
      registry.fill(HIST("data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda"), jet.pt(), nV0inJet, nK0SinJet, nLambdainJet, nAntiLambdainJet);

      // registry.fill(HIST("data/jets/V0/jetCorrectedPtEtaPhi"), newjet.pt(), newjet.eta(), newjet.phi());
    }
  }
  PROCESS_SWITCH(JetFragmentation, processDataV0Frag, "Data V0 fragmentation", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetFragmentation>(cfgc, TaskName{"jet-fragmentation"})};
}
