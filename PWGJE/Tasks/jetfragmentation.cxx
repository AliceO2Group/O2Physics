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

// jet V0 fragmentation
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
#include "PWGJE/Core/JetFindingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Charged jets
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using ChargedJetsWithConstituents = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;

using MCDJets = aod::ChargedMCDetectorLevelJets;
using MCDJetsWithConstituents = soa::Join<MCDJets, aod::ChargedMCDetectorLevelJetConstituents>;
using MatchedMCDJets = soa::Join<MCDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;
using MatchedMCDJetsWithConstituents = soa::Join<MCDJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;

using MCPJets = aod::ChargedMCParticleLevelJets;
using MatchedMCPJets = soa::Join<MCPJets, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;
using MCPJetsWithConstituents = soa::Join<MCPJets, aod::ChargedMCParticleLevelJetConstituents>;
using MatchedMCPJetsWithConstituents = soa::Join<MCPJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;

// V0 jets
using MCDV0Jets = aod::V0ChargedMCDetectorLevelJets;
using MCDV0JetsWithConstituents = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetConstituents>;
using MatchedMCDV0Jets = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;
using MatchedMCDV0JetsWithConstituents = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetConstituents, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;

using MCPV0Jets = aod::V0ChargedMCParticleLevelJets;
using MCPV0JetsWithConstituents = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetConstituents>;
using MatchedMCPV0Jets = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;
using MatchedMCPV0JetsWithConstituents = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetConstituents, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;

struct JetFragmentation {
  HistogramRegistry registry{"registry"};

  Configurable<std::string> evSel{"evSel", "sel8WithoutTimeFrameBorderCut", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.f, "vertex z cut"};

  Configurable<float> matchedDetJetEtaMin{"matchedDetJetEtaMin", -0.5, "minimum matchedDetJet eta"};
  Configurable<float> matchedDetJetEtaMax{"matchedDetJetEtaMax", 0.5, "maximum matchedDetJet eta"};
  Configurable<float> dataJetEtaMin{"dataJetEtaMin", -0.5, "minimum data jet eta"};
  Configurable<float> dataJetEtaMax{"dataJetEtaMax", 0.5, "maximum data jet eta"};
  Configurable<float> v0EtaMin{"v0EtaMin", -0.75, "minimum data V0 eta"};
  Configurable<float> v0EtaMax{"v0EtaMax", 0.75, "maximum data V0 eta"};

  Configurable<double> v0cospaMin{"v0cospaMin", 0.995, "V0 CosPA"};
  Configurable<double> dcav0dauMax{"dcav0dauMax", 1.0, "DCA V0 Daughters"};
  Configurable<double> dcaprMin{"dcaprMin", 0.1, "DCA proton To PV"};
  Configurable<double> dcapiMin{"dcapiMin", 0.1, "DCA pion To PV"};
  Configurable<double> v0radiusMin{"v0radiusMin", 1.2, "V0 Radius"};
  Configurable<double> lifetimeK0SMax{"lifetimeK0SMax", 20., "lifetimeK0SMax"};
  Configurable<double> lifetimeLambdaMax{"lifetimeLambdaMax", 25., "lifetimeLambdaMax"};

  Configurable<double> k0sMassAccWindow{"k0sMassAccWindow", 0.03, "k0sMassAccWindow"};
  Configurable<double> lambdaMassAccWindow{"lambdaMassAccWindow", 0.01, "lambdaMassAccWindow"};
  Configurable<double> antilambdaMassAccWindow{"antilambdaMassAccWindow", 0.01, "antilambdaMassAccWindow"};

  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {40, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {20, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
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

  ConfigurableAxis binV0Pt{"binV0Pt", {120, 0.0f, 60.0f}, ""};
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

  // Binning for cut variation study
  ConfigurableAxis binV0RadiusCut{"binV0RadiusCut", {4, 1.0f, 1.4f}, "R"};
  ConfigurableAxis binV0CtauCut{"binV0CtauCut", {3, 15.0f, 30.0f}, "c#tau"};
  ConfigurableAxis binV0CosPACut{"binV0CosPACut", {4, 0.991f, 0.999f}, "cosPA"};
  ConfigurableAxis binV0DCApCut{"binV0DCApCut", {2, 0.05f, 0.15f}, "DCA pos"};
  ConfigurableAxis binV0DCAnCut{"binV0DCAnCut", {2, 0.05f, 0.15f}, "DCA neg"};
  ConfigurableAxis binV0DCAdCut{"binV0DCAdCut", {2, 0.5f, 1.5f}, "DCA daughters"};
  ConfigurableAxis binV0PtCut{"binV0PtCut", {60, 0.0f, 60.0f}, "p_{T, V0}"};
  ConfigurableAxis binK0SMassCut{"binK0SMassCut", {100, 0.4f, 0.6f}, "inv. mass, K0S hypothesis"};
  ConfigurableAxis binLambda0MassCut{"binLambda0MassCut", {100, 1.07f, 1.21f}, "inv. mass, Lambda0 hypothesis"};
  ConfigurableAxis binAntiLambda0MassCut{"binAntiLambda0MassCut", {100, 1.07f, 1.21f}, "inv. mass, AntiLambda0 hypothesis"};

  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;

  Partition<MatchedMCDJetsWithConstituents> detJetEtaPartition = (aod::jet::eta > matchedDetJetEtaMin) && (aod::jet::eta < matchedDetJetEtaMax);
  Partition<MatchedMCDJetsWithConstituents> detJetEtaV0Partition = (aod::jet::eta > v0EtaMin + aod::jet::r * 0.01f) && (aod::jet::eta < v0EtaMax - aod::jet::r * 0.01f);

  Preslice<MyTracks> TracksPerCollision = aod::track::collisionId;
  Preslice<aod::V0Datas> V0sPerCollision = aod::v0data::collisionId;
  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> McV0sPerCollision = aod::v0data::collisionId;
  Preslice<MCPJetsWithConstituents> PartJetsPerCollision = aod::jet::mcCollisionId;
  Preslice<JetParticles> JetParticlesPerCollision = aod::jmcparticle::mcCollisionId;
  Preslice<aod::McParticles> ParticlesPerCollision = aod::mcparticle::mcCollisionId;

  int eventSelection = -1;

  void init(InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(evSel));

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
    AxisSpec zRelDiffAxis = {binZRelDiff, "(#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part})/#it{p}_{T, jet}^{part}"};

    AxisSpec V0PtAxis = {binV0Pt, "#it{p}_{T}^{V0}"};
    AxisSpec V0PtRatioAxis = {binPtRatio, "#it{p}_{T}^{V0, det}/#it{p}_{T, V0}^{part}"};
    AxisSpec V0PtRelDiffAxis = {binPtRelDiff, "(#it{p}_{T}^{V0, det} - #it{p}_{T}^{V0, part})/#it{p}_{T, V0}^{part}"};
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

    // Cut variation study
    AxisSpec RcutAxis = {binV0RadiusCut, "R"};
    AxisSpec ctauCutAxis = {binV0CtauCut, "c#tau (K0S)"};
    AxisSpec cosPACutAxis = {binV0CosPACut, "cosPA"};
    AxisSpec DCApCutAxis = {binV0DCApCut, "DCA pos (cm)"};
    AxisSpec DCAnCutAxis = {binV0DCAnCut, "DCA neg (cm)"};
    AxisSpec DCAdCutAxis = {binV0DCAdCut, "DCA daughters (cm^{2})"};
    AxisSpec PtCutAxis = {binV0PtCut, "p_{T, V0}"};
    AxisSpec K0SMassCutAxis = {binK0SMassCut, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec LambdaMassCutAxis = {binLambda0MassCut, "Inv. mass (GeV/#it{c}^{2})"};
    AxisSpec AntiLambdaMassCutAxis = {binAntiLambda0MassCut, "Inv. mass (GeV/#it{c}^{2})"};

    if (doprocessDataRun3 || doprocessDataV0Frag || doprocessDataV0JetsFrag) {
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
    } // doprocessDataRun3 || doprocessDataV0Frag

    if (doprocessDataV0 || doprocessDataV0Frag || doprocessDataV0JetsFrag) {
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

      registry.add("data/V0/V0CutVariation", "V0CutVariation", HistType::kTHnSparseD, {PtCutAxis, K0SMassCutAxis, LambdaMassCutAxis, AntiLambdaMassCutAxis, RcutAxis, ctauCutAxis, cosPACutAxis, DCApCutAxis, DCAnCutAxis, DCAdCutAxis});
    } // doprocessDataV0 || doprocessDataV0Frag

    if (doprocessDataV0Frag || doprocessDataV0JetsFrag) {
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
      registry.add("data/jets/V0/jetPtLambdaPtAllMasses", "jetPtLambdaPtAllMasses", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaPtLambdaMasses", "jetPtLambdaPtLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtLambdaPtRadius", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, radius", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtLambdaPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtLambdaPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, #Lambda^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtLambdaTrackProjCtau", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjMass", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjAllMasses", "jetPtLambdaTrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjLambdaMasses", "jetPtLambdaTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjRadius", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtLambdaTrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{#Lambda^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtnAntiLambda", "jetPtnAntiLambda", HistType::kTH2D, {jetPtAxis, trackCount});
      registry.add("data/jets/V0/jetPtAntiLambdaPtCtau", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtMass", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, mass", HistType::kTH3D, {jetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtAllMasses", "jetPtAntiLambdaPtAllMasses", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtLambdaMasses", "jetPtAntiLambdaPtLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtRadius", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, radius", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, #bar{#Lambda}^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjCtau", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjMass", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjAllMasses", "jetPtAntiLambdaTrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjLambdaMasses", "jetPtAntiLambdaTrackProjLambdaMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjRadius", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtAntiLambdaTrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{#bar{#Lambda}^{0}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtnK0S", "jetPtnK0S", HistType::kTH2D, {jetPtAxis, trackCount});
      registry.add("data/jets/V0/jetPtK0SPtCtau", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, c#tau", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtK0SPtMass", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, mass", HistType::kTH3D, {jetPtAxis, V0PtAxis, K0SMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtAllMasses", "jetPtK0SPtAllMasses", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtK0SPtRadius", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, radius", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtK0SPtCosPA", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, cosPA", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAd", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, DCA daughters", HistType::kTH3D, {jetPtAxis, V0PtAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtK0SPtDCAposneg", "Jet #it{p}_{T}, #it{p}_{T, K^{0}_{S}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});

      registry.add("data/jets/V0/jetPtK0STrackProjCtau", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, c#tau", HistType::kTH3D, {jetPtAxis, zAxis, V0CtauAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjMass", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, mass", HistType::kTH3D, {jetPtAxis, zAxis, K0SMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjAllMasses", "jetPtK0STrackProjAllMasses", HistType::kTHnSparseD, {jetPtAxis, zAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjRadius", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, radius", HistType::kTH3D, {jetPtAxis, zAxis, V0RadiusAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjCosPA", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, cosPA", HistType::kTH3D, {jetPtAxis, zAxis, V0CosPAAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAd", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, DCA daughters", HistType::kTH3D, {jetPtAxis, zAxis, V0DCAdAxis});
      registry.add("data/jets/V0/jetPtK0STrackProjDCAposneg", "Jet #it{p}_{T}, #it{z}_{K^{0}_{S}}, DCA#pm", HistType::kTHnSparseD, {jetPtAxis, zAxis, V0DCApAxis, V0DCAnAxis});
    } // doprocessDataV0Frag

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
    } // doprocessMcP

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
    } // doprocessMcD

    if (doprocessMcMatched || doprocessMcMatchedV0Frag || doprocessMcMatchedV0JetsFrag) {
      registry.add("matching/jets/matchDetJetPtEtaPhi", "Matched detector level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {detJetPtAxis, detEtaAxis, detPhiAxis});
      registry.add("matching/jets/matchPartJetPtEtaPhi", "Matched particle level jet #it{p}_{T}, #eta, #phi", HistType::kTH3D, {partJetPtAxis, partEtaAxis, partPhiAxis});
      registry.add("matching/jets/matchDetJetPtPartJetPt", "matchDetJetPtPartJetPt", HistType::kTH2D, {detJetPtAxis, partJetPtAxis});
      registry.add("matching/jets/matchPartJetPtDetJetEtaPartJetEta", "matchPartJetPtDetJetEtaPartJetEta", HistType::kTH3D, {partJetPtAxis, detEtaAxis, partEtaAxis});
      registry.add("matching/jets/matchPartJetPtDetJetPhiPartJetPhi", "matchPartJetPtDetJetPhiPartJetPhi", HistType::kTH3D, {partJetPtAxis, detPhiAxis, partPhiAxis});
      registry.add("matching/jets/matchPartJetPtResolutionPt", "#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part}", HistType::kTH2D, {partJetPtAxis, ptDiffAxis});
      registry.add("matching/jets/matchPartJetPtRelDiffPt", "#it{p}_{T}^{jet, det} - #it{p}_{T}^{jet, part}", HistType::kTH2D, {partJetPtAxis, ptJetRelDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionEta", "#eta^{jet, det} - #eta^{jet, part}", HistType::kTH3D, {partJetPtAxis, partEtaAxis, etaDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionPhi", "#phi^{jet, det} - #phi^{jet, part}", HistType::kTH3D, {partJetPtAxis, partPhiAxis, phiDiffAxis});
      registry.add("matching/jets/matchPartJetPtEtaPhiMatchDist", "matchJetMatchDist", HistType::kTHnSparseD, {partJetPtAxis, partEtaAxis, partPhiAxis, matchDistAxis});
      registry.add("matching/jets/matchPartJetPtEnergyScale", "jetEnergyScale", HistType::kTH2D, {partJetPtAxis, ptRatioAxis});
    } // doprocessMcMatched || doprocessMcMatchedV0Frag

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
      registry.add("matching/jets/matchDetJetPtTrackPt", "Matched detector level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2D, {detJetPtAxis, trackPtAxis});
      registry.add("matching/jets/matchDetJetTrackPtEtaPhi", "Matched detector level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3D, {trackPtAxis, detEtaAxis, detPhiAxis});
      registry.add("matching/jets/matchDetJetPtFrag", "Matched detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/matchDetJetPtTrackProj", "Matched detector level jet #it{p}_{T}, #it{z}", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/matchDetJetPtXi", "Matched detector level jet #it{p}_{T}, #xi", HistType::kTH2D, {detJetPtAxis, detXiAxis});
      registry.add("matching/jets/matchDetJetPtTheta", "Matched detector level jet #it{p}_{T}, #theta", HistType::kTH2D, {detJetPtAxis, detThetaAxis});
      registry.add("matching/jets/matchDetJetPtXiTheta", "Matched detector level jet #it{p}_{T}, #xi, #theta", HistType::kTH3D, {detJetPtAxis, detXiAxis, detThetaAxis});
      registry.add("matching/jets/matchDetJetPtZTheta", "Matched detector level jet #it{p}_{T}, z, #theta", HistType::kTH3D, {detJetPtAxis, detZAxis, detThetaAxis});
      // Particle level jets with a match
      registry.add("matching/jets/matchPartJetPtTrackPt", "Matched particle level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2D, {partJetPtAxis, trackPtAxis});
      registry.add("matching/jets/matchPartJetTrackPtEtaPhi", "Matched particle level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3D, {trackPtAxis, partEtaAxis, partPhiAxis});
      registry.add("matching/jets/matchPartJetPtFrag", "Matched particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/matchPartJetPtTrackProj", "Matched particle level jet #it{p}_{T}, #it{z}", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/matchPartJetPtXi", "Matched particle level jet #it{p}_{T}, #xi", HistType::kTH2D, {partJetPtAxis, partXiAxis});
      registry.add("matching/jets/matchPartJetPtTheta", "Matched particle level jet #it{p}_{T}, #theta", HistType::kTH2D, {partJetPtAxis, partThetaAxis});
      registry.add("matching/jets/matchPartJetPtXiTheta", "Matched particle level jet #it{p}_{T}, #xi, #theta", HistType::kTH3D, {partJetPtAxis, partXiAxis, partThetaAxis});
      registry.add("matching/jets/matchPartJetPtZTheta", "Matched particle level jet #it{p}_{T}, z, #theta", HistType::kTH3D, {partJetPtAxis, partZAxis, partThetaAxis});
      // Combined information of matched jets
      registry.add("matching/jets/matchPartJetPtResolutionChargeFrag", "Resolution #it{p}_{T}^{tr} / #it{p}_{T}^{jet}", HistType::kTH3D, {partJetPtAxis, partZAxis, zDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionTrackPt", "Resolution #it{p}_{T}^{track}", HistType::kTH3D, {partJetPtAxis, trackPtAxis, ptTrackDiffAxis});
      registry.add("matching/jets/matching/jets/matchPartJetPtRelDiffTrackPt", "Rel. diff #it{p}_{T}^{track}", HistType::kTHnSparseD, {partJetPtAxis, ptRatioAxis, trackPtAxis, ptTrackRelDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionTrackProj", "Resolution #it{p}^{proj} / #it{p}^{jet}", HistType::kTH3D, {partJetPtAxis, partZAxis, zDiffAxis});
      registry.add("matching/jets/matchPartJetPtRelDiffTrackProj", "Rel. diff #it{p}^{proj} / #it{p}^{jet}", HistType::kTHnSparseD, {partJetPtAxis, ptRatioAxis, partZAxis, zRelDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionXi", "Resolution ln(1/#it{z})", HistType::kTH3D, {partJetPtAxis, partXiAxis, xiDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionTheta", "Resolution #theta", HistType::kTH3D, {partJetPtAxis, partThetaAxis, thetaDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionXiResolutionTheta", "Resolution #xi, #theta", HistType::kTHnSparseD, {partJetPtAxis, partXiAxis, xiDiffAxis, partThetaAxis, thetaDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionZResolutionTheta", "Resolution #it{z}, #theta", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, zDiffAxis, partThetaAxis, thetaDiffAxis});
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

      registry.add("matching/jets/matchDetJetPtXiThetaPartJetPtXiTheta", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detXiAxis, detThetaAxis, partJetPtAxis, partXiAxis, partThetaAxis});
      registry.add("matching/jets/fakeDetJetPtXiTheta", "Fakes", HistType::kTH3D, {detJetPtAxis, detXiAxis, detThetaAxis});
      registry.add("matching/jets/missPartJetPtXiTheta", "Misses", HistType::kTH3D, {partJetPtAxis, partXiAxis, partThetaAxis});

      registry.add("matching/jets/matchDetJetPtZThetaPartJetPtZTheta", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, detThetaAxis, partJetPtAxis, partZAxis, partThetaAxis});
      registry.add("matching/jets/fakeDetJetPtZTheta", "Fakes", HistType::kTH3D, {detJetPtAxis, detZAxis, detThetaAxis});
      registry.add("matching/jets/missPartJetPtZTheta", "Misses", HistType::kTH3D, {partJetPtAxis, partZAxis, partThetaAxis});
    } // doprocessMcMatched

    if (doprocessMcMatchedV0 || doprocessMcMatchedV0Frag || doprocessMcMatchedV0JetsFrag) {
      registry.add("matching/V0/nV0sEvent", "nV0sDet per event", HistType::kTH1D, {v0Count});
    } // doprocessMcMatchedV0 || doprocessMcMatchedV0Frag

    if (doprocessMcMatchedV0) {
      registry.add("matching/V0/V0PartPtDetPt", "V0PartPtDetPt", HistType::kTH2D, {V0partPtAxis, V0detPtAxis});
      registry.add("matching/V0/V0PartPtRatioPtRelDiffPt", "V0PartPtRatioRelDiffPt", HistType::kTH3D, {V0partPtAxis, V0PtRatioAxis, V0PtRelDiffAxis});
      registry.add("matching/V0/V0PosPartPtRatioPtRelDiffPt", "V0PosPartPtRatioRelDiffPt", HistType::kTH3D, {trackPtAxis, ptRatioAxis, ptTrackRelDiffAxis});
      registry.add("matching/V0/V0NegPartPtRatioPtRelDiffPt", "V0NegPartPtRatioRelDiffPt", HistType::kTH3D, {trackPtAxis, ptRatioAxis, ptTrackRelDiffAxis});

      registry.add("matching/V0/K0SPtEtaPhi", "K0SPtEtaPhi", HistType::kTHnSparseD, {V0partPtAxis, V0detPtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/V0/K0SPtCtauMass", "K0SPtCtauMass", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0CtauAxis, K0SMassAxis});
      registry.add("matching/V0/K0SPtRadiusCosPA", "K0SPtRadiusCosPA", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/V0/K0SPtDCAposneg", "K0SPtDCAposneg", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/V0/K0SPtDCAd", "K0SPtDCAd", HistType::kTH3D, {V0partPtAxis, V0detPtAxis, V0DCAdAxis});
      registry.add("matching/V0/K0SPtMass", "K0SPtMass", HistType::kTHnSparseD, {V0partPtAxis, V0detPtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});

      registry.add("matching/V0/LambdaPtEtaPhi", "LambdaPtEtaPhi", HistType::kTHnSparseD, {V0partPtAxis, V0detPtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/V0/LambdaPtCtauMass", "LambdaPtCtauMass", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0CtauAxis, LambdaMassAxis});
      registry.add("matching/V0/LambdaPtRadiusCosPA", "LambdaPtRadiusCosPA", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/V0/LambdaPtDCAposneg", "LambdaPtDCAposneg", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/V0/LambdaPtDCAd", "LambdaPtDCAd", HistType::kTH3D, {V0partPtAxis, V0detPtAxis, V0DCAdAxis});
      registry.add("matching/V0/LambdaPtMass", "LambdaPtMass", HistType::kTHnSparseD, {V0partPtAxis, V0detPtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});

      registry.add("matching/V0/antiLambdaPtEtaPhi", "antiLambdaPtEtaPhi", HistType::kTHnSparseD, {V0partPtAxis, V0detPtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/V0/antiLambdaPtCtauMass", "antiLambdaPtCtauMass", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0CtauAxis, LambdaMassAxis});
      registry.add("matching/V0/antiLambdaPtRadiusCosPA", "antiLambdaPtRadiusCosPA", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/V0/antiLambdaPtDCAposneg", "antiLambdaPtDCAposneg", HistType::kTHnSparseD, {V0partPtAxis, V0partPtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/V0/antiLambdaPtDCAd", "antiLambdaPtDCAd", HistType::kTH3D, {V0partPtAxis, V0detPtAxis, V0DCAdAxis});
      registry.add("matching/V0/antiLambdaPtMass", "antiLambdaPtMass", HistType::kTHnSparseD, {V0partPtAxis, V0detPtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});

      // Reflections
      registry.add("matching/V0/Lambda0Reflection", "pt, pt, mK, mL, maL, Lambda0Reflection", HistType::kTHnSparseD, {V0partPtAxis, V0detPtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/V0/antiLambda0Reflection", "pt, pt, mK, mL, maL, antiLambda0Reflection", HistType::kTHnSparseD, {V0partPtAxis, V0detPtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis, LambdaMassAxis});
    } // doprocessMcMatchedV0

    if (doprocessMcMatchedV0Frag || doprocessMcMatchedV0JetsFrag) {
      registry.add("matching/jets/V0/jetPtnV0Matched", "jet pt, nV0 matched", HistType::kTH2D, {detJetPtAxis, v0Count});
      registry.add("matching/jets/V0/jetPtnV0MatchednK0SnLambdanAntiLambda", "jet Pt, nV0 matched, nK0S nLambdan AntiLambda", HistType::kTHnSparseD, {detJetPtAxis, v0Count, v0Count, v0Count, v0Count});
      registry.add("matching/jets/V0/partJetPtV0PtDetPt", "V0PartPtDetPt", HistType::kTH3D, {partJetPtAxis, V0partPtAxis, V0detPtAxis});
      registry.add("matching/jets/V0/partJetPtDetJetPtPartV0PtRatioPtRelDiffPt", "V0PartPtRatioRelDiffPt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, V0partPtAxis, V0PtRatioAxis, V0PtRelDiffAxis});
      registry.add("matching/jets/V0/partJetPtDetJetPtPartV0PtPosPtRatioPtRelDiffPt", "V0PtPosPartPtRatioRelDiffPt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, V0PtAxis, trackPtAxis, ptRatioAxis, ptTrackRelDiffAxis});
      registry.add("matching/jets/V0/partJetPtDetJetPtPartV0PtNegPtRatioPtRelDiffPt", "V0PtNegPartPtRatioRelDiffPt", HistType::kTHnSparseD, {partJetPtAxis, detJetPtAxis, V0PtAxis, trackPtAxis, ptRatioAxis, ptTrackRelDiffAxis});

      // -----------------------------
      // Unidentified V0s
      // -----------------------------
      registry.add("matching/jets/V0/matchDetJetPtV0TrackProjPartJetPtV0TrackProj", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0Pt", "matched jet Pt, V0 Pt", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis});
      // Matched V0: pt
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauLambda0", "matched jet Pt, V0 Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauAntiLambda0", "matched jet Pt, V0 Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauK0S", "matched jet Pt, V0 Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassLambda0", "matched jet Pt, V0 Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassAntiLambda0", "matched jet Pt, V0 Pt, Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassK0S", "matched jet Pt, V0 Pt, MassK0S", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, K0SMassAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtRadius", "matched jet Pt, V0 Pt, Radius", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCosPA", "matched jet Pt, V0 Pt, CosPA", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtDCAposneg", "matched jet Pt, V0 Pt, DCAposneg", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/partJetPtV0PtDetJetPtV0PtDCAd", "matched jet Pt, V0 Pt, DCAd", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0DCAdAxis});
      // Matched Lambda0: z
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauLambda0", "matched jet Pt, V0 z, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauAntiLambda0", "matched jet Pt, V0 z, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauK0S", "matched jet Pt, V0 z, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassLambda0", "matched jet Pt, V0 z, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassAntiLambda0", "matched jet Pt, V0 z, Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassK0S", "matched jet Pt, V0 z, MassK0S", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjRadius", "matched jet Pt, V0 z, Radius", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0RadiusAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCosPA", "matched jet Pt, V0 z, CosPA", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjDCAposneg", "matched jet Pt, V0 z, DCAposneg", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjDCAd", "matched jet Pt, V0 z, DCAd", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0DCAdAxis});
      // Fakes
      registry.add("matching/jets/V0/fakeJetPtV0TrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/V0/fakeJetPtV0PtEtaPhi", "fake jet Pt, V0 PtEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/jets/V0/fakeJetPtV0PtCtau", "fake jet Pt, V0 PtCtau", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("matching/jets/V0/fakeJetPtV0PtMass", "fake jet Pt, V0 PtMass", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/fakeJetPtV0PtLambdaMasses", "fake jet Pt, V0 PtLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("matching/jets/V0/fakeJetPtV0PtRadiusCosPA", "fake jet Pt, V0 PtRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/fakeJetPtV0PtDCAposneg", "fake jet Pt, V0 PtDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/fakeJetPtV0PtDCAd", "fake jet Pt, V0 PtDCAd", HistType::kTH3D, {detJetPtAxis, V0PtAxis, V0DCAdAxis});
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjCtau", "fake jet Pt, V0 zCtau", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjMass", "fake jet Pt, V0 zMass", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjLambdaMasses", "fake jet Pt, V0 zLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjRadiusCosPA", "fake jet Pt, V0 zRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjDCAposneg", "fake jet Pt, V0 zDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/fakeJetPtV0TrackProjDCAd", "fake jet Pt, V0 zDCAd", HistType::kTH3D, {detJetPtAxis, detZAxis, V0DCAdAxis});
      // Misses
      registry.add("matching/jets/V0/missJetPtV0PtEtaPhi", "miss jet Pt, V0 PtEtaPhi", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/jets/V0/missJetPtV0TrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis});

      // -----------------------------
      // Lambda0
      // -----------------------------
      registry.add("matching/jets/V0/matchDetJetPtLambda0TrackProjPartJetPtLambda0TrackProj", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0Pt", "matched jet Pt, #Lambda^{0} Pt", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis});
      // Matched Lambda0: pt
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtCtauLambda0", "matched jet Pt, #Lambda^{0} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtCtauAntiLambda0", "matched jet Pt, #Lambda^{0} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtMassLambda0", "matched jet Pt, #Lambda^{0} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtMassAntiLambda0", "matched jet Pt, #Lambda^{0} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtMassK0S", "matched jet Pt, #Lambda^{0} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, K0SMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtAllMasses", "matched jet Pt, #Lambda^{0} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtRadius", "matched jet Pt, #Lambda^{0} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtCosPA", "matched jet Pt, #Lambda^{0} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtDCAposneg", "matched jet Pt, #Lambda^{0} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtDCAd", "matched jet Pt, #Lambda^{0} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0DCAdAxis});
      // Matched Lambda0: z
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjCtauLambda0", "matched jet Pt, #Lambda^{0} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjCtauAntiLambda0", "matched jet Pt, #Lambda^{0} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjMassLambda0", "matched jet Pt, #Lambda^{0} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjMassAntiLambda0", "matched jet Pt, #Lambda^{0} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjMassK0S", "matched jet Pt, #Lambda^{0} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjAllMasses", "matched jet Pt, #Lambda^{0} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjRadius", "matched jet Pt, #Lambda^{0} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0RadiusAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjCosPA", "matched jet Pt, #Lambda^{0} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjDCAposneg", "matched jet Pt, #Lambda^{0} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjDCAd", "matched jet Pt, #Lambda^{0} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0DCAdAxis});
      // Fake Lambda0
      registry.add("matching/jets/V0/fakeJetPtLambda0TrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0PtEtaPhi", "fake jet Pt, #Lambda^{0} PtEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0PtCtau", "fake jet Pt, #Lambda^{0} PtCtau", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0PtMass", "fake jet Pt, #Lambda^{0} PtMass", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0PtLambdaMasses", "fake jet Pt, #Lambda^{0} PtLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0PtRadiusCosPA", "fake jet Pt, #Lambda^{0} PtRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0PtDCAposneg", "fake jet Pt, #Lambda^{0} PtDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0PtDCAd", "fake jet Pt, #Lambda^{0} PtDCAd", HistType::kTH3D, {detJetPtAxis, V0PtAxis, V0DCAdAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0TrackProjEtaPhi", "fake jet Pt, #Lambda^{0} zEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0TrackProjCtau", "fake jet Pt, #Lambda^{0} zCtau", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0TrackProjMass", "fake jet Pt, #Lambda^{0} zMass", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0TrackProjLambdaMasses", "fake jet Pt, #Lambda^{0} zLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0TrackProjRadiusCosPA", "fake jet Pt, #Lambda^{0} zRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0TrackProjDCAposneg", "fake jet Pt, #Lambda^{0} zDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/fakeJetPtLambda0TrackProjDCAd", "fake jet Pt, #Lambda^{0} zDCAd", HistType::kTH3D, {detJetPtAxis, detZAxis, V0DCAdAxis});
      // Missed Lambda0
      registry.add("matching/jets/V0/missJetPtLambda0TrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/V0/missJetPtLambda0PtEtaPhi", "miss jet Pt, #Lambda^{0} PtEtaPhi", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});

      // -----------------------------
      // AntiLambda0
      // -----------------------------
      registry.add("matching/jets/V0/matchDetJetPtAntiLambda0TrackProjPartJetPtAntiLambda0TrackProj", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0Pt", "matched jet Pt, #bar{#Lambda}^{0} Pt", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis});
      // Matched AntiLambda0: pt
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtCtauLambda0", "matched jet Pt, #bar{#Lambda}^{0} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtCtauAntiLambda0", "matched jet Pt, #bar{#Lambda}^{0} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtMassLambda0", "matched jet Pt, #bar{#Lambda}^{0} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtMassAntiLambda0", "matched jet Pt, #bar{#Lambda}^{0} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtMassK0S", "matched jet Pt, #bar{#Lambda}^{0} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, K0SMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtAllMasses", "matched jet Pt, #bar{#Lambda}^{0} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtRadius", "matched jet Pt, #bar{#Lambda}^{0} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtCosPA", "matched jet Pt, #bar{#Lambda}^{0} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtDCAposneg", "matched jet Pt, #bar{#Lambda}^{0} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtDCAd", "matched jet Pt, #bar{#Lambda}^{0} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0DCAdAxis});
      // Matched AntiLambda0: z
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjCtauLambda0", "matched jet Pt, #bar{#Lambda}^{0} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjCtauAntiLambda0", "matched jet Pt, #bar{#Lambda}^{0} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjMassLambda0", "matched jet Pt, #bar{#Lambda}^{0} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjMassAntiLambda0", "matched jet Pt, #bar{#Lambda}^{0} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjMassK0S", "matched jet Pt, #bar{#Lambda}^{0} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjAllMasses", "matched jet Pt, #bar{#Lambda}^{0} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjRadius", "matched jet Pt, #bar{#Lambda}^{0} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0RadiusAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjCosPA", "matched jet Pt, #bar{#Lambda}^{0} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjDCAposneg", "matched jet Pt, #bar{#Lambda}^{0} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjDCAd", "matched jet Pt, #bar{#Lambda}^{0} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0DCAdAxis});
      // Fake AntiLambda0
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0TrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0PtEtaPhi", "fake jet Pt, #bar{#Lambda}^{0} PtEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0PtCtau", "fake jet Pt, #bar{#Lambda}^{0} PtCtau", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0PtMass", "fake jet Pt, #bar{#Lambda}^{0} PtMass", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0PtLambdaMasses", "fake jet Pt, #bar{#Lambda}^{0} PtLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0PtRadiusCosPA", "fake jet Pt, #bar{#Lambda}^{0} PtRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0PtDCAposneg", "fake jet Pt, #bar{#Lambda}^{0} PtDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0PtDCAd", "fake jet Pt, #bar{#Lambda}^{0} PtDCAd", HistType::kTH3D, {detJetPtAxis, V0PtAxis, V0DCAdAxis});

      registry.add("matching/jets/V0/fakeJetPtAntiLambda0TrackProjCtau", "fake jet Pt, #bar{#Lambda}^{0} zCtau", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0TrackProjMass", "fake jet Pt, #bar{#Lambda}^{0} zMass", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0TrackProjLambdaMasses", "fake jet Pt, #bar{#Lambda}^{0} zLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0TrackProjRadiusCosPA", "fake jet Pt, #bar{#Lambda}^{0} zRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0TrackProjDCAposneg", "fake jet Pt, #bar{#Lambda}^{0} zDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/fakeJetPtAntiLambda0TrackProjDCAd", "fake jet Pt, #bar{#Lambda}^{0} zDCAd", HistType::kTH3D, {detJetPtAxis, detZAxis, V0DCAdAxis});
      // Missed AntiLambda0
      registry.add("matching/jets/V0/missJetPtAntiLambda0TrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/V0/missJetPtAntiLambda0PtEtaPhi", "miss jet Pt, #bar{#Lambda}^{0} PtEtaPhi", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});

      // -----------------------------
      // K0S
      // -----------------------------
      registry.add("matching/jets/V0/matchDetJetPtK0STrackProjPartJetPtK0STrackProj", "Matched", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPt", "matched jet Pt, K_{S}^{0} Pt", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis});
      // Matched K0S: pt
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauLambda0", "matched jet Pt, K^{0}_{S} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauAntiLambda0", "matched jet Pt, K^{0}_{S} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassLambda0", "matched jet Pt, K^{0}_{S} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassAntiLambda0", "matched jet Pt, K^{0}_{S} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassK0S", "matched jet Pt, K^{0}_{S} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, K0SMassAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtAllMasses", "matched jet Pt, K^{0}_{S} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtRadius", "matched jet Pt, K^{0}_{S} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0RadiusAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCosPA", "matched jet Pt, K^{0}_{S} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtDCAposneg", "matched jet Pt, K^{0}_{S} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtDCAd", "matched jet Pt, K^{0}_{S} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, detJetPtAxis, V0PtAxis, V0DCAdAxis});
      // Matched K0S: z
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauLambda0", "matched jet Pt, K^{0}_{S} Pt, Ctau #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauAntiLambda0", "matched jet Pt, K^{0}_{S} Pt, Ctau #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauK0S", "matched jet Pt, #{K}^{0}_{S} Pt, Ctau #{K}^{0}_{S}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CtauAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassLambda0", "matched jet Pt, K^{0}_{S} Pt, Mass #Lambda^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassAntiLambda0", "matched jet Pt, K^{0}_{S} Pt Mass #bar{#Lambda}^{0}", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassK0S", "matched jet Pt, K^{0}_{S} PtMassK0S", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjAllMasses", "matched jet Pt, K^{0}_{S} Pt, Masses", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjRadius", "matched jet Pt, K^{0}_{S} Pt Radius", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0RadiusAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCosPA", "matched jet Pt, K^{0}_{S} Pt CosPA", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjDCAposneg", "matched jet Pt, K^{0}_{S} PtDCAposneg", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjDCAd", "matched jet Pt, K^{0}_{S} PtDCAd", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, V0DCAdAxis});
      // Fake K0S
      registry.add("matching/jets/V0/fakeJetPtK0STrackProj", "Fakes", HistType::kTH2D, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/V0/fakeJetPtK0SPtEtaPhi", "fake jet Pt, K^{0}_{S} PtEtaPhi", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});
      registry.add("matching/jets/V0/fakeJetPtK0SPtCtau", "fake jet Pt, K^{0}_{S} PtCtau", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("matching/jets/V0/fakeJetPtK0SPtMass", "fake jet Pt, K^{0}_{S} PtMass", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/fakeJetPtK0SPtLambdaMasses", "fake jet Pt, K^{0}_{S} PtLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("matching/jets/V0/fakeJetPtK0SPtRadiusCosPA", "fake jet Pt, K^{0}_{S} PtRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/fakeJetPtK0SPtDCAposneg", "fake jet Pt, K^{0}_{S} PtDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, V0PtAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/fakeJetPtK0SPtDCAd", "fake jet Pt, K^{0}_{S} PtDCAd", HistType::kTH3D, {detJetPtAxis, V0PtAxis, V0DCAdAxis});

      registry.add("matching/jets/V0/fakeJetPtK0STrackProjCtau", "fake jet Pt, K^{0}_{S} zCtau", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0CtauAxis, V0CtauAxis, V0CtauAxis});
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjMass", "fake jet Pt, K^{0}_{S} zMass", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjLambdaMasses", "fake jet Pt, K^{0}_{S} zLambdaMasses", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, LambdaMassDiffAxis, LambdaMassRatioAxis, LambdaMassRelDiffAxis});
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjRadiusCosPA", "fake jet Pt, K^{0}_{S} zRadiusCosPA", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0RadiusAxis, V0CosPAAxis});
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjDCAposneg", "fake jet Pt, K^{0}_{S} zDCAposneg", HistType::kTHnSparseD, {detJetPtAxis, detZAxis, V0DCApAxis, V0DCAnAxis});
      registry.add("matching/jets/V0/fakeJetPtK0STrackProjDCAd", "fake jet Pt, K^{0}_{S} zDCAd", HistType::kTH3D, {detJetPtAxis, detZAxis, V0DCAdAxis});
      // Missed K0S
      registry.add("matching/jets/V0/missJetPtK0STrackProj", "Misses", HistType::kTH2D, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/V0/missJetPtK0SPtEtaPhi", "miss jet Pt, K^{0}_{S} PtEtaPhi", HistType::kTHnSparseD, {partJetPtAxis, V0PtAxis, V0EtaAxis, V0PhiAxis});

      // Reflections
      registry.add("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtLambda0Reflection", "Lambda0 Reflection", HistType::kTHnSparseD, {partJetPtAxis, V0partPtAxis, detJetPtAxis, V0detPtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjLambda0Reflection", "Lambda0 Reflection", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtAntiLambda0Reflection", "antiLambda0 Reflection", HistType::kTHnSparseD, {partJetPtAxis, V0partPtAxis, detJetPtAxis, V0detPtAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis, LambdaMassAxis});
      registry.add("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjAntiLambda0Reflection", "antiLambda0 Reflection", HistType::kTHnSparseD, {partJetPtAxis, partZAxis, detJetPtAxis, detZAxis, K0SMassAxis, LambdaMassAxis, LambdaMassAxis, LambdaMassAxis});
    } // doprocessMcMatchedV0Frag
  }   // init

  template <typename JetType>
  bool JetContainsV0s(JetType const& jet)
  {
    return (jet.hfcandidatesIds().size() > 0);
  }
  template <typename T, typename U, typename V>
  bool V0sAreMatched(T const& v0, U const& particle, V const& /*tracks*/)
  {
    auto negId = v0.template negTrack_as<V>().template mcParticleId();
    auto posId = v0.template posTrack_as<V>().template mcParticleId();
    auto daughters = particle.daughtersIds();

    return ((negId == daughters[0] && posId == daughters[1]) || (posId == daughters[0] && negId == daughters[1]));
  }
  template <typename V0Type>
  bool IsV0Candidate(V0Type const& v0)
  {
    if (v0.eta() < v0EtaMin || v0.eta() > v0EtaMax) { // TODO: Should be rapidity, mass matters!
      return false;
    }
    if (v0.dcaV0daughters() > dcav0dauMax) {
      return false;
    }
    if (v0.v0radius() < v0radiusMin) {
      return false;
    }
    if (v0.v0cosPA() < v0cospaMin) {
      return false;
    }
    return true;
  }
  template <typename CollisionType, typename V0Type>
  bool IsK0SCandidate(CollisionType const& collision, V0Type const& v0)
  {
    if (!IsV0Candidate(v0)) {
      return false;
    }
    if (v0.dcanegtopv() < dcapiMin || v0.dcapostopv() < dcapiMin) {
      return false;
    }
    double ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    if (ctauK0s > lifetimeK0SMax) {
      return false;
    }
    bool k0sMassCondition = (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < k0sMassAccWindow);
    if (!k0sMassCondition) {
      return false;
    }
    return true;
  }
  template <typename CollisionType, typename V0Type>
  bool IsLambdaCandidate(CollisionType const& collision, V0Type const& v0)
  {
    if (!IsV0Candidate(v0)) {
      return false;
    }
    if (v0.dcanegtopv() < dcapiMin || v0.dcapostopv() < dcaprMin) {
      return false;
    }
    double ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (ctauLambda > lifetimeLambdaMax) {
      return false;
    }
    bool lambdaMassCondition = (TMath::Abs(v0.mLambda() - o2::constants::physics::MassLambda0) < lambdaMassAccWindow);
    if (!lambdaMassCondition) {
      return false;
    }
    return true;
  }
  template <typename CollisionType, typename V0Type>
  bool IsAntiLambdaCandidate(CollisionType const& collision, V0Type const& v0)
  {
    if (!IsV0Candidate(v0)) {
      return false;
    }
    if (v0.dcanegtopv() < dcaprMin || v0.dcapostopv() < dcapiMin) {
      return false;
    }
    double ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    if (ctauAntiLambda > lifetimeLambdaMax) {
      return false;
    }
    bool antilambdaMassCondition = (TMath::Abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0Bar) < antilambdaMassAccWindow);
    if (!antilambdaMassCondition) {
      return false;
    }
    return true;
  }

  template <typename V0Type>
  double ReflectedMass(V0Type const& v0, bool isLambda)
  {
    // If V0 is Lambda, posTrack = proton, negTrack = pion
    // In that case, we assign pion mass to posTrack and proton mass to negTrack to calculate the reflection
    // Vice versa for AntiLambda
    double negM = (isLambda ? o2::constants::physics::MassProton : o2::constants::physics::MassPionCharged);
    double posM = (isLambda ? o2::constants::physics::MassPionCharged : o2::constants::physics::MassProton);
    double negPsq = v0.pxneg() * v0.pxneg() + v0.pyneg() * v0.pyneg() + v0.pzneg() * v0.pzneg();
    double posPsq = v0.pxpos() * v0.pxpos() + v0.pypos() * v0.pypos() + v0.pzpos() * v0.pzpos();
    double negE = TMath::Sqrt(negM * negM + negPsq);
    double posE = TMath::Sqrt(posM * posM + posPsq);
    double Esquared = (negE + posE) * (negE + posE);
    double psquared = v0.p() * v0.p();
    return TMath::Sqrt(Esquared - psquared);
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

  // TODO: Can probably be made simpler/shorter by using V0MCLabels
  template <typename CollisionType, typename V0Type, typename trackType, typename particleType>
  void fillMcMatchedV0Histograms(CollisionType const& collision, V0Type const& v0, trackType const&, particleType const&, double weight = 1.)
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
          registry.fill(HIST("matching/V0/V0PartPtRatioPtRelDiffPt"), ptPartV0, v0.pt() / ptPartV0, (v0.pt() - ptPartV0) / ptPartV0);

          if (pdg == 310) { // K0S
            registry.fill(HIST("matching/V0/K0SPtEtaPhi"), ptPartV0, v0.pt(), v0.eta(), v0.phi());
            registry.fill(HIST("matching/V0/K0SPtCtauMass"), ptPartV0, v0.pt(), ctauK0s, v0.mK0Short(), weight);
            registry.fill(HIST("matching/V0/K0SPtRadiusCosPA"), ptPartV0, v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
            registry.fill(HIST("matching/V0/K0SPtDCAposneg"), ptPartV0, v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
            registry.fill(HIST("matching/V0/K0SPtDCAd"), ptPartV0, v0.pt(), v0.dcaV0daughters(), weight);
            registry.fill(HIST("matching/V0/K0SPtMass"), ptPartV0, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
          } else if (pdg == 3122) { // Lambda
            registry.fill(HIST("matching/V0/LambdaPtEtaPhi"), ptPartV0, v0.pt(), v0.eta(), v0.phi());
            registry.fill(HIST("matching/V0/LambdaPtCtauMass"), ptPartV0, v0.pt(), ctauLambda, v0.mLambda(), weight);
            registry.fill(HIST("matching/V0/LambdaPtRadiusCosPA"), ptPartV0, v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
            registry.fill(HIST("matching/V0/LambdaPtDCAposneg"), ptPartV0, v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
            registry.fill(HIST("matching/V0/LambdaPtDCAd"), ptPartV0, v0.pt(), v0.dcaV0daughters(), weight);
            registry.fill(HIST("matching/V0/LambdaPtMass"), ptPartV0, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);

            // Reflection
            double reflectedMass = ReflectedMass(v0, true);
            registry.fill(HIST("matching/V0/Lambda0Reflection"), ptPartV0, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
          } else if (pdg == -3122) { // AntiLambda
            registry.fill(HIST("matching/V0/antiLambdaPtEtaPhi"), ptPartV0, v0.pt(), v0.eta(), v0.phi());
            registry.fill(HIST("matching/V0/antiLambdaPtCtauMass"), ptPartV0, v0.pt(), ctauAntiLambda, v0.mAntiLambda(), weight);
            registry.fill(HIST("matching/V0/antiLambdaPtRadiusCosPA"), ptPartV0, v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
            registry.fill(HIST("matching/V0/antiLambdaPtDCAposneg"), ptPartV0, v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
            registry.fill(HIST("matching/V0/antiLambdaPtDCAd"), ptPartV0, v0.pt(), v0.dcaV0daughters(), weight);
            registry.fill(HIST("matching/V0/antiLambdaPtMass"), ptPartV0, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);

            // Reflection
            double reflectedMass = ReflectedMass(v0, false);
            registry.fill(HIST("matching/V0/antiLambda0Reflection"), ptPartV0, v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
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
      double chargeFrag = -1., trackProj = -1., xi = -1., theta = -1.;
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

  template <typename CollisionType, typename V0Type>
  void fillDataV0Histograms(CollisionType const& collision, V0Type const& V0s)
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

      registry.fill(HIST("data/V0/V0CutVariation"), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), v0.v0radius(), ctauK0s, v0.v0cosPA(), TMath::Abs(v0.dcapostopv()), TMath::Abs(v0.dcanegtopv()), v0.dcaV0daughters());

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
  void fillDataV0FragHistograms(CollisionType const& collision, JetType const& jet, V0Type const& v0)
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
      registry.fill(HIST("data/jets/V0/jetPtK0SPtAllMasses"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtK0SPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());

      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjCtau"), jet.pt(), trackProj, ctauK0s);
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjMass"), jet.pt(), trackProj, v0.mK0Short());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjAllMasses"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjRadius"), jet.pt(), trackProj, v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjCosPA"), jet.pt(), trackProj, v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtK0STrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
    }
    if (IsLambdaCandidate(collision, v0)) {
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtCtau"), jet.pt(), v0.pt(), ctauLambda);
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtMass"), jet.pt(), v0.pt(), v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtAllMasses"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtLambdaPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());

      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjCtau"), jet.pt(), trackProj, ctauLambda);
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjMass"), jet.pt(), trackProj, v0.mLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjAllMasses"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjRadius"), jet.pt(), trackProj, v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjCosPA"), jet.pt(), trackProj, v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtLambdaTrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv());
    }
    if (IsAntiLambdaCandidate(collision, v0)) {
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtCtau"), jet.pt(), v0.pt(), ctauAntiLambda);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtMass"), jet.pt(), v0.pt(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtAllMasses"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtRadius"), jet.pt(), v0.pt(), v0.v0radius());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtCosPA"), jet.pt(), v0.pt(), v0.v0cosPA());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaPtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv());

      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjCtau"), jet.pt(), trackProj, ctauAntiLambda);
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjMass"), jet.pt(), trackProj, v0.mAntiLambda());
      registry.fill(HIST("data/jets/V0/jetPtAntiLambdaTrackProjAllMasses"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda());
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

  template <typename JetType, typename V0Type>
  void fillMatchingV0Miss(JetType const& jet, V0Type const& v0, double weight = 1.)
  {
    double trackProj = TrackProj(jet, v0);

    registry.fill(HIST("matching/jets/V0/missJetPtV0TrackProj"), jet.pt(), trackProj, weight);
    registry.fill(HIST("matching/jets/V0/missJetPtV0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
    if (v0.pdgCode() == 310) { // K0S
      registry.fill(HIST("matching/jets/V0/missJetPtK0SPtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/missJetPtK0STrackProj"), jet.pt(), trackProj, weight);
    } else if (v0.pdgCode() == 3122) { // Lambda
      registry.fill(HIST("matching/jets/V0/missJetPtLambda0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/missJetPtLambda0TrackProj"), jet.pt(), trackProj, weight);
    } else if (v0.pdgCode() == -3122) { // AntiLambda
      registry.fill(HIST("matching/jets/V0/missJetPtAntiLambda0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/missJetPtAntiLambda0TrackProj"), jet.pt(), trackProj, weight);
    }
  }

  template <typename CollisionType, typename JetType, typename V0Type>
  void fillMatchingV0Fake(CollisionType const& collision, JetType const& jet, V0Type const& v0, double weight = 1.)
  {
    double trackProj = TrackProj(jet, v0);
    double ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
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

    if (IsLambdaCandidate(collision, v0)) {
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0TrackProj"), jet.pt(), trackProj, weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0PtCtau"), jet.pt(), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0PtMass"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0PtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0PtRadiusCosPA"), jet.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0PtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0PtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0TrackProjCtau"), jet.pt(), trackProj, ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0TrackProjMass"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0TrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0TrackProjRadiusCosPA"), jet.pt(), trackProj, v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0TrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtLambda0TrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters(), weight);
    }
    if (IsAntiLambdaCandidate(collision, v0)) {
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0PtEtaPhi"), jet.pt(), v0.pt(), v0.eta(), v0.phi(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0TrackProj"), jet.pt(), trackProj, weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0PtCtau"), jet.pt(), v0.pt(), ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0PtMass"), jet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0PtLambdaMasses"), jet.pt(), v0.pt(), massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0PtRadiusCosPA"), jet.pt(), v0.pt(), v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0PtDCAposneg"), jet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0PtDCAd"), jet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0TrackProjCtau"), jet.pt(), trackProj, ctauK0s, ctauLambda, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0TrackProjMass"), jet.pt(), trackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0TrackProjLambdaMasses"), jet.pt(), trackProj, massDiff, massRatio, massRelDiff, weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0TrackProjRadiusCosPA"), jet.pt(), trackProj, v0.v0radius(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0TrackProjDCAposneg"), jet.pt(), trackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/fakeJetPtAntiLambda0TrackProjDCAd"), jet.pt(), trackProj, v0.dcaV0daughters(), weight);
    }
    if (IsK0SCandidate(collision, v0)) {
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
  }

  template <typename CollisionType, typename DetJetType, typename PartJetType, typename V0Type, typename ParticleType>
  // template <typename CollisionType, typename DetJetType, typename PartJetType, typename V0Type, typename ParticleType, typename V0DaughterType, typename ParticleDaughterType>
  void fillMatchingHistogramsV0(CollisionType const& collision, DetJetType const& detJet, PartJetType const& partJet, V0Type const& v0, ParticleType const& particle, double weight = 1.)
  {
    double detTrackProj = TrackProj(detJet, v0);
    double partTrackProj = TrackProj(partJet, particle);

    double ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    double ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    double ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

    registry.fill(HIST("matching/jets/V0/matchDetJetPtV0TrackProjPartJetPtV0TrackProj"), detJet.pt(), detTrackProj, partJet.pt(), partTrackProj, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0Pt"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetPt"), partJet.pt(), particle.pt(), detJet.pt(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtDetJetPtPartV0PtRatioPtRelDiffPt"), partJet.pt(), detJet.pt(), particle.pt(), v0.pt() / particle.pt(), (v0.pt() - particle.pt()) / particle.pt(), weight);

    // V0CandidatesMCP contains daughter list. How do we know which is pos?
    // auto negTrack = v0.template negTrack_as<V0DaughterType>();
    // auto posTrack = v0.template posTrack_as<V0DaughterType>();
    // auto negPart = negTrack.template particle_as<ParticleDaughterType>();
    // auto posPart = posTrack.template particle_as<ParticleDaughterType>();
    // registry.fill(HIST("matching/jets/V0/partJetPtDetJetPtPartV0PtPosPtRatioPtRelDiffPt"), partJet.pt(), detJet.pt(), particle.pt(), posTrack.pt(), posTrack.pt() / posPart.pt(), (posTrack.pt() - posPart.pt()) / posPart.pt(), weight);
    // registry.fill(HIST("matching/jets/V0/partJetPtDetJetPtPartV0PtNegPtRatioPtRelDiffPt"), partJet.pt(), detJet.pt(), particle.pt(), negTrack.pt(), negTrack.pt() / negPart.pt(), (negTrack.pt() - negPart.pt()) / negPart.pt(), weight);

    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauLambda, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauAntiLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauAntiLambda, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCtauK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauK0s, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mLambda(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassAntiLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mAntiLambda(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtMassK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtRadius"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0radius(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtCosPA"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0cosPA(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtDCAposneg"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0PtDetJetPtV0PtDCAd"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauLambda, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauAntiLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauAntiLambda, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCtauK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauK0s, weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mLambda(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassAntiLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mAntiLambda(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjMassK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjRadius"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0radius(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjCosPA"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0cosPA(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjDCAposneg"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
    registry.fill(HIST("matching/jets/V0/partJetPtV0TrackProjDetJetPtV0TrackProjDCAd"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcaV0daughters(), weight);

    if (particle.pdgCode() == 310) { // K0S
      registry.fill(HIST("matching/jets/V0/matchDetJetPtK0STrackProjPartJetPtK0STrackProj"), detJet.pt(), detTrackProj, partJet.pt(), partTrackProj, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPt"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauAntiLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCtauK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassAntiLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtMassK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtAllMasses"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtRadius"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtCosPA"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtDCAposneg"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0SPtDetJetPtK0SPtDCAd"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauAntiLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCtauK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassAntiLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjMassK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjAllMasses"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjRadius"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjCosPA"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjDCAposneg"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtK0STrackProjDetJetPtK0STrackProjDCAd"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcaV0daughters(), weight);
    } else if (particle.pdgCode() == 3122) { // Lambda
      registry.fill(HIST("matching/jets/V0/matchDetJetPtLambda0TrackProjPartJetPtLambda0TrackProj"), detJet.pt(), detTrackProj, partJet.pt(), partTrackProj, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0Pt"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtCtauLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtCtauAntiLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtCtauK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtMassLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtMassAntiLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtMassK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtAllMasses"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtRadius"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtCosPA"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtDCAposneg"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtDCAd"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjCtauLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjCtauAntiLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjCtauK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjMassLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjMassAntiLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjMassK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjAllMasses"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjRadius"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjCosPA"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjDCAposneg"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjDCAd"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcaV0daughters(), weight);

      // Reflection
      double reflectedMass = ReflectedMass(v0, true);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0PtLambda0Reflection"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtLambda0TrackProjDetJetPtLambda0TrackProjLambda0Reflection"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
    } else if (particle.pdgCode() == -3122) { // AntiLambda
      registry.fill(HIST("matching/jets/V0/matchDetJetPtAntiLambda0TrackProjPartJetPtAntiLambda0TrackProj"), detJet.pt(), detTrackProj, partJet.pt(), partTrackProj, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0Pt"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtCtauLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtCtauAntiLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtCtauK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtMassLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtMassAntiLambda0"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtMassK0S"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtAllMasses"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtRadius"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtCosPA"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtDCAposneg"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtDCAd"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.dcaV0daughters(), weight);

      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjCtauLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjCtauAntiLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauAntiLambda, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjCtauK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, ctauK0s, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjMassLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjMassAntiLambda0"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjMassK0S"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjAllMasses"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjRadius"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0radius(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjCosPA"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.v0cosPA(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjDCAposneg"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcapostopv(), v0.dcanegtopv(), weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjDCAd"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.dcaV0daughters(), weight);

      // Reflection
      double reflectedMass = ReflectedMass(v0, false);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0PtAntiLambda0Reflection"), partJet.pt(), particle.pt(), detJet.pt(), v0.pt(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
      registry.fill(HIST("matching/jets/V0/partJetPtAntiLambda0TrackProjDetJetPtAntiLambda0TrackProjAntiLambda0Reflection"), partJet.pt(), partTrackProj, detJet.pt(), detTrackProj, v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(), reflectedMass, weight);
    } // AntiLambda
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

  void processDummy(JetTracks const&) {}
  PROCESS_SWITCH(JetFragmentation, processDummy, "Dummy process function turned on by default", true);

  void processMcD(soa::Filtered<JetCollisionsMCD>::iterator const& collision,
                  JetMcCollisions const&,
                  MCDJetsWithConstituents const&,
                  JetTracks const& tracks)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    double nJets = 0, nTracks = 0;
    double weight = collision.mcCollision().weight();
    for (const auto& track : tracks) {
      if (track.pt() > 0.1) {
        nTracks++;
        registry.fill(HIST("detector-level/tracks/detTrackPtEtaPhi"), track.pt(), track.eta(), track.phi(), weight);
      }
    }
    for (const auto& jet : detJetEtaPartition) {
      nJets++;
      fillMCDHistograms(jet, weight);
    }
    registry.fill(HIST("detector-level/nJetsnTracks"), nJets, nTracks, weight);
  }
  PROCESS_SWITCH(JetFragmentation, processMcD, "Monte Carlo detector level", false);

  void processMcP(JetMcCollision const& mcCollision,
                  MCPJetsWithConstituents const& jets,
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

  void processDataRun3(soa::Filtered<JetCollisions>::iterator const& collision,
                       ChargedJetsWithConstituents const& jets,
                       JetTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    double nJets = 0, nTracks = 0;
    for (const auto& track : tracks) {
      if (track.pt() > 0.1) {
        nTracks++;
        registry.fill(HIST("data/tracks/trackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      }
    }
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, dataJetEtaMin, dataJetEtaMax)) {
        continue;
      }
      nJets++;
      fillDataRun3Histograms(jet);
    }
    registry.fill(HIST("data/nJetsnTracks"), nJets, nTracks);
  }
  PROCESS_SWITCH(JetFragmentation, processDataRun3, "Run 3 Data", false);

  void processMcMatched(soa::Filtered<JetCollisionsMCD>::iterator const& collision,
                        MatchedMCDJetsWithConstituents const&,
                        JetTracksMCD const&,
                        JetMcCollisions const&,
                        MatchedMCPJetsWithConstituents const& allMcPartJets,
                        JetParticles const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    double weight = collision.mcCollision().weight();
    const auto& mcPartJets = allMcPartJets.sliceBy(PartJetsPerCollision, collision.mcCollision().globalIndex()); // Only jets from the same collision
    bool isFake = false;
    for (const auto& detJet : detJetEtaPartition) {
      for (auto& partJet : detJet.template matchedJetGeo_as<MatchedMCPJetsWithConstituents>()) {
        fillMatchingHistogramsJet(detJet, partJet, weight);

        for (const auto& track : detJet.tracks_as<JetTracksMCD>()) {
          bool isTrackMatched = false;
          if (!track.has_mcParticle()) {
            isFake = true;
            fillMatchingFakeOrMiss(detJet, track, isFake, weight);
            continue;
          }
          for (const auto& particle : partJet.tracks_as<JetParticles>()) {
            if (particle.globalIndex() == track.template mcParticle_as<JetParticles>().globalIndex()) {
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
    for (const auto& partJet : mcPartJets) {
      for (const auto& detJet : partJet.template matchedJetGeo_as<MatchedMCDJetsWithConstituents>()) {
        // Check if the matched detector level jet is outside the allowed eta range
        if ((detJet.eta() <= matchedDetJetEtaMin) || (detJet.eta() >= matchedDetJetEtaMax)) {
          for (const auto& particle : partJet.tracks_as<JetParticles>()) {
            isFake = false;
            fillMatchingFakeOrMiss(partJet, particle, isFake, weight);
          }
          continue;
        }
        // If the jets are properly matched, we can check the particles
        for (const auto& particle : partJet.tracks_as<JetParticles>()) {
          bool isParticleMatched = false;
          for (const auto& track : detJet.tracks_as<JetTracksMCD>()) {
            if (!track.has_mcParticle()) {
              continue;
            }
            if (particle.globalIndex() == track.template mcParticle_as<JetParticles>().globalIndex()) {
              isParticleMatched = true;
            }
          }
          // Ignore matched particles. They have been handled in the previous loop
          if (!isParticleMatched) {
            isFake = false;
            fillMatchingFakeOrMiss(partJet, particle, isFake, weight);
          }
        } // for particle
      }   // for matched det jet
      if (!partJet.has_matchedJetGeo()) {
        isFake = false;
        registry.fill(HIST("matching/jets/missPartJetPtEtaPhi"), partJet.pt(), partJet.eta(), partJet.phi(), weight);
        for (const auto& particle : partJet.tracks_as<JetParticles>()) {
          fillMatchingFakeOrMiss(partJet, particle, isFake, weight);
        }
      } // if no matched jet
    }   // for part jet
  }
  PROCESS_SWITCH(JetFragmentation, processMcMatched, "Monte Carlo particle and detector level", false);

  // Should take in JCollisions?
  void processMcMatchedV0(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels>>::iterator const& collision,
                          aod::McCollisions const&,
                          soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                          soa::Join<MyTracks, aod::McTrackLabels> const& tracks,
                          aod::McParticles const& mcParticles)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    // if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
    //   return;
    // }
    double weight = collision.mcCollision().weight();
    for (const auto& v0 : V0s) {
      if (!v0.has_mcParticle()) {
        continue;
      }
      fillMcMatchedV0Histograms(collision, v0, tracks, mcParticles, weight);
    }
  }
  PROCESS_SWITCH(JetFragmentation, processMcMatchedV0, "Monte Carlo V0", false);

  void processMcMatchedV0Frag(soa::Filtered<soa::Join<JetCollisionsMCD, aod::JCollisionPIs>>::iterator const& jcoll,
                              MatchedMCDJetsWithConstituents const&,
                              JetTracksMCD const&,
                              soa::Join<aod::V0Datas, aod::McV0Labels> const& allV0s,
                              JetMcCollisions const&,
                              MatchedMCPJetsWithConstituents const& allMcPartJets,
                              JetParticles const&,
                              aod::McCollisions const&,
                              aod::McParticles const& allMcParticles,
                              aod::Collisions const&)
  {
    if (!jcoll.has_mcCollision()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelection)) {
      return;
    }
    double weight = jcoll.mcCollision().weight();
    // This is necessary, because jets are linked to JetCollisions, but V0s are linked to Collisions
    const auto& collision = jcoll.collision_as<aod::Collisions>();
    const auto& v0s = allV0s.sliceBy(V0sPerCollision, collision.globalIndex());
    const auto& mcPartJets = allMcPartJets.sliceBy(PartJetsPerCollision, jcoll.mcCollision().globalIndex());
    const auto& mcParticles = allMcParticles.sliceBy(ParticlesPerCollision, jcoll.mcCollision().globalIndex());

    int kNV0s = v0s.size();
    bool isV0Used[kNV0s];
    for (int i = 0; i < kNV0s; i++) {
      isV0Used[i] = false;
    }
    registry.fill(HIST("matching/V0/nV0sEvent"), kNV0s);

    int kNParticles = mcParticles.size();
    bool isParticleUsed[kNParticles];
    for (int i = 0; i < kNParticles; i++) {
      isParticleUsed[i] = false;
    }

    for (const auto& detJet : detJetEtaV0Partition) {
      int iv0 = -1;
      int nV0inJet = 0, nLambdainJet = 0, nAntiLambdainJet = 0, nK0SinJet = 0;

      for (auto& partJet : detJet.template matchedJetGeo_as<MatchedMCPJetsWithConstituents>()) {
        fillMatchingHistogramsJet(detJet, partJet, weight);
        // Jets are pt-sorted, so we prioritise matching V0s with high pt jets
        for (const auto& v0 : v0s) {
          iv0++;
          if (isV0Used[iv0]) {
            continue;
          }
          double dR = jetutilities::deltaR(detJet, v0);
          if (dR >= detJet.r() * 1e-2) {
            continue;
          }
          isV0Used[iv0] = true;
          if (!v0.has_mcParticle()) {
            fillMatchingV0Fake(collision, detJet, v0, weight);
            continue;
          }
          const auto& particle = v0.template mcParticle_as<aod::McParticles>();
          if (!((TMath::Abs(particle.pdgCode()) == 310) || (TMath::Abs(particle.pdgCode()) == 3122))) {
            fillMatchingV0Fake(collision, detJet, v0, weight);
            continue;
          }
          // Found a matched V0 in the jet
          nV0inJet++;
          fillMatchingHistogramsV0(collision, detJet, partJet, v0, particle, weight);
          if (TMath::Abs(particle.pdgCode()) == 310) {
            nK0SinJet++;
          } else if (particle.pdgCode() == 3122) {
            nLambdainJet++;
          } else if (particle.pdgCode() == -3122) {
            nAntiLambdainJet++;
          }
        } // v0 loop
        registry.fill(HIST("matching/jets/V0/jetPtnV0Matched"), partJet.pt(), nV0inJet, weight);
        registry.fill(HIST("matching/jets/V0/jetPtnV0MatchednK0SnLambdanAntiLambda"), partJet.pt(), nV0inJet, nK0SinJet, nLambdainJet, nAntiLambdainJet, weight);
      } // for partJet in matched detJet
      iv0 = -1;
      if (!detJet.has_matchedJetGeo()) {
        for (const auto& v0 : v0s) {
          iv0++;
          if (isV0Used[iv0]) {
            continue;
          }
          double dR = jetutilities::deltaR(detJet, v0);
          if (dR >= detJet.r() * 1e-2) {
            continue;
          }
          isV0Used[iv0] = true;
          fillMatchingV0Fake(collision, detJet, v0, weight);
        } // v0 loop
      }   // if no matched jet
    }     // det jet loop
    for (const auto& partJet : mcPartJets) {
      int iparticle = -1;
      for (const auto& particle : mcParticles) {
        iparticle++;
        if (isParticleUsed[iparticle]) {
          continue;
        }
        // Check if particle is primary and is a particle of interest that has not been used yet
        // If it doesn't pass these selections, set isParticleUsed to true to skip it in the future
        if (!particle.isPhysicalPrimary()) {
          isParticleUsed[iparticle] = true;
          continue;
        }
        if (!((TMath::Abs(particle.pdgCode()) == 310) || TMath::Abs((particle.pdgCode()) == 3122))) {
          isParticleUsed[iparticle] = true;
          continue;
        }
        // If the particle has been used or it is not a particle of interest, skip it
        if (isParticleUsed[iparticle]) {
          continue;
        }
        if (jetutilities::deltaR(partJet, particle) >= partJet.r() * 1e-2) {
          continue;
        }
        // Particle may be a miss, but we need to check if it is matched with a V0 in a detector level jet
        // If it is, it has been treated in the loop over detector level jets above
        if (!partJet.has_matchedJetGeo()) {
          isParticleUsed[iparticle] = true;
          fillMatchingV0Miss(partJet, particle, weight);
          continue;
        }
        for (const auto& detJet : partJet.template matchedJetGeo_as<MatchedMCDJetsWithConstituents>()) {
          if ((detJet.eta() <= v0EtaMin + detJet.r() * 1e-2) || (detJet.eta() >= v0EtaMax - detJet.r() * 1e-2)) {
            continue;
          }
          for (const auto& v0 : v0s) {
            if (!v0.has_mcParticle()) {
              continue;
            }
            if (v0.template mcParticle_as<aod::McParticles>().globalIndex() == particle.globalIndex()) {
              if (jetutilities::deltaR(detJet, v0) < detJet.r() * 1e-2) {
                // The particle is matched with a V0 and we ignore it
                isParticleUsed[iparticle] = true;
              }
            }
          } // v0 loop
        }   // detJet loop
        if (!isParticleUsed[iparticle]) {
          isParticleUsed[iparticle] = true;
          fillMatchingV0Miss(partJet, particle, weight);
        }
      } // particle loop
    }   // part jet loop
  }
  PROCESS_SWITCH(JetFragmentation, processMcMatchedV0Frag, "Monte Carlo V0 fragmentation", false);

  void processDataV0(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                     aod::V0Datas const& V0s)
  {
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("data/V0/nV0sEvent"), V0s.size());
    fillDataV0Histograms(collision, V0s);
  }
  PROCESS_SWITCH(JetFragmentation, processDataV0, "Data V0", false);

  void processDataV0Frag(soa::Filtered<soa::Join<JetCollisions, aod::JCollisionPIs>>::iterator const& jcoll,
                         ChargedJetsWithConstituents const& jets,
                         JetTracks const&,
                         aod::Collisions const&,
                         aod::V0Datas const& allV0s)
  {
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelection)) {
      return;
    }
    // This is necessary, because jets are linked to JetCollisions, but V0s are linked to Collisions
    const auto& collision = jcoll.collision_as<aod::Collisions>();
    const auto& v0s = allV0s.sliceBy(V0sPerCollision, collision.globalIndex());

    int kNV0s = v0s.size();
    bool isV0Used[kNV0s];
    for (int i = 0; i < kNV0s; i++) {
      isV0Used[i] = false;
    }
    registry.fill(HIST("data/V0/nV0sEvent"), kNV0s);

    fillDataV0Histograms(collision, v0s);
    for (const auto& jet : jets) {
      if ((jet.eta() < v0EtaMin + jet.r() * 1e-2) || (jet.eta() > v0EtaMax - jet.r() * 1e-2)) {
        continue;
      }
      fillDataRun3Histograms(jet);
      // fastjet::PseudoJet newjet(jet.px(), jet.py(), jet.pz(), jet.e()); // Jet with corrections from V0
      int iv0 = -1;
      int nV0inJet = 0, nLambdainJet = 0, nAntiLambdainJet = 0, nK0SinJet = 0;

      // Jets are pt-sorted, so we prioritise matching V0s with high pt jets
      // Correct jet momentum (currently only corrects for v0 in jet, not v0 outside jet, is this an issue?)
      // for (const auto& v0 : v0s) {
      //   iv0++;
      //   if (isV0Used[iv0]) {
      //     continue;
      //   }
      //   double dR = jetutilities::deltaR(jet, v0);
      //   if (dR < jet.r() * 1e-2) {
      //     // fastjet::PseudoJet pjv0(v0.px(), v0.py(), v0.pz(), v0.e());
      //     // newjet += pjv0;
      //   }
      // }
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

  // ---------------- V0 jets ----------------
  void processDataV0JetsFrag(soa::Filtered<JetCollisions>::iterator const& jcoll, soa::Join<aod::V0ChargedJets, aod::V0ChargedJetConstituents> const& v0jets, CandidatesV0Data const& v0s)
  {
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelection)) {
      return;
    }
    registry.fill(HIST("data/V0/nV0sEvent"), v0s.size());
    fillDataV0Histograms(jcoll, v0s);

    for (const auto& jet : v0jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, -99., -99., v0EtaMin, v0EtaMax)) {
        continue;
      }
      // Double check if the jet contains V0s
      if (!JetContainsV0s(jet)) {
        continue;
      }

      int nV0inJet = 0, nLambdainJet = 0, nAntiLambdainJet = 0, nK0SinJet = 0;
      for (const auto& v0 : jet.hfcandidates_as<CandidatesV0Data>()) {
        nV0inJet++;
        fillDataV0FragHistograms(jcoll, jet, v0);
        if (IsK0SCandidate(jcoll, v0)) {
          nK0SinJet++;
        }
        if (IsLambdaCandidate(jcoll, v0)) {
          nLambdainJet++;
        }
        if (IsAntiLambdaCandidate(jcoll, v0)) {
          nAntiLambdainJet++;
        }
      }
      registry.fill(HIST("data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda"), jet.pt(), nV0inJet, nK0SinJet, nLambdainJet, nAntiLambdainJet);
    } // Jet loop
  }
  PROCESS_SWITCH(JetFragmentation, processDataV0JetsFrag, "Data V0 jets fragmentation", false);

  void processMcMatchedV0JetsFrag(soa::Filtered<JetCollisionsMCD>::iterator const& jcoll, JetMcCollisions const&, MatchedMCDV0JetsWithConstituents const& v0jetsMCD, MatchedMCPV0JetsWithConstituents const& v0jetsMCP, soa::Join<CandidatesV0MCD, aod::McV0Labels> const& v0s, CandidatesV0MCP const&, JetTracksMCD const& jTracks)
  {
    if (!jcoll.has_mcCollision()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelection)) {
      return;
    }
    double weight = jcoll.mcCollision().weight();
    registry.fill(HIST("matching/V0/nV0sEvent"), v0s.size());

    for (const auto& detJet : v0jetsMCD) {
      if (!jetfindingutilities::isInEtaAcceptance(detJet, -99., -99., v0EtaMin, v0EtaMax)) {
        continue;
      }
      // Double check if the jet contains V0s
      if (!JetContainsV0s(detJet)) {
        continue;
      }
      int nV0inJet = 0, nLambdainJet = 0, nAntiLambdainJet = 0, nK0SinJet = 0;
      if (!detJet.has_matchedJetGeo()) {
        for (const auto& detV0 : detJet.hfcandidates_as<soa::Join<CandidatesV0MCD, aod::McV0Labels>>()) {
          fillMatchingV0Fake(jcoll, detJet, detV0, weight);
        }
        continue;
      } // if jet not matched

      for (const auto& partJet : detJet.template matchedJetGeo_as<MatchedMCPV0JetsWithConstituents>()) {
        fillMatchingHistogramsJet(detJet, partJet, weight);
        for (const auto& detV0 : detJet.hfcandidates_as<soa::Join<CandidatesV0MCD, aod::McV0Labels>>()) {
          if (!detV0.has_mcParticle()) {
            fillMatchingV0Fake(jcoll, detJet, detV0, weight);
            continue;
          }
          bool isV0Matched = false;
          for (const auto& partV0 : partJet.template hfcandidates_as<CandidatesV0MCP>()) {
            if (V0sAreMatched(detV0, partV0, jTracks)) {
              isV0Matched = true;
              nV0inJet++;
              fillMatchingHistogramsV0(jcoll, detJet, partJet, detV0, partV0, weight);
              if (TMath::Abs(partV0.pdgCode()) == 310) {
                nK0SinJet++;
              } else if (partV0.pdgCode() == 3122) {
                nLambdainJet++;
              } else if (partV0.pdgCode() == -3122) {
                nAntiLambdainJet++;
              }
              break;
            } // if matched
          }   // partV0 loop

          if (!isV0Matched) {
            fillMatchingV0Fake(jcoll, detJet, detV0, weight);
          }
        } // detV0 loop
        registry.fill(HIST("matching/jets/V0/jetPtnV0MatchednK0SnLambdanAntiLambda"), partJet.pt(), nV0inJet, nK0SinJet, nLambdainJet, nAntiLambdainJet, weight);
      } // Matched partJet loop
    }   // detJet loop

    for (const auto& partJet : v0jetsMCP) {
      if (!JetContainsV0s(partJet)) {
        continue;
      }
      if (!partJet.has_matchedJetGeo()) {
        for (const auto& partV0 : partJet.hfcandidates_as<CandidatesV0MCP>()) {
          fillMatchingV0Miss(partJet, partV0, weight);
        }
        continue;
      } // if jet not matched

      bool isJetMatched = false;
      for (const auto& detJet : partJet.template matchedJetGeo_as<MatchedMCDV0JetsWithConstituents>()) {
        if (!jetfindingutilities::isInEtaAcceptance(detJet, -99., -99., v0EtaMin, v0EtaMax)) {
          continue;
        }
        isJetMatched = true;
        for (const auto& partV0 : partJet.hfcandidates_as<CandidatesV0MCP>()) {
          bool isV0Matched = false;
          for (const auto& detV0 : detJet.hfcandidates_as<soa::Join<CandidatesV0MCD, aod::McV0Labels>>()) {
            if (V0sAreMatched(detV0, partV0, jTracks)) {
              isV0Matched = true;
              break;
            }
          } // detV0 loop
          if (!isV0Matched) {
            fillMatchingV0Miss(partJet, partV0, weight);
          }
        } // partV0 loop
      }   // detJet loop

      // To account for matched jets where the detector level jet is outside of the eta range (cut applied within this task)
      if (!isJetMatched) {
        for (const auto& partV0 : partJet.hfcandidates_as<CandidatesV0MCP>()) {
          fillMatchingV0Miss(partJet, partV0, weight);
        }
      }
    } // partJet loop
  }
  PROCESS_SWITCH(JetFragmentation, processMcMatchedV0JetsFrag, "Matched V0 jets fragmentation", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetFragmentation>(cfgc, TaskName{"jet-fragmentation"})};
}
