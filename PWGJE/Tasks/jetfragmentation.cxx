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

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using McTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>>;
using McTracks = soa::Join<aod::JTracks, aod::JMcTrackLbs>;
using McDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
// using MatchedMcDJets = soa::Filtered<soa::Join<McDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
using McPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

struct JetFragmentation {
  HistogramRegistry registry{"registry"};

  std::vector<int> pdgVector = {211, 321, 2212, 111, 130, 310, 311, 3122};
  std::vector<std::string> hadronVector = {"#pi^{#pm}", "#it{K}^{#pm}", "#it{p}^{#pm}", "#pi^{0}", "#it{K}^{0}_{L}", "#it{K}^{0}_{S}", "#it{K}^{0}", "#Lambda^{0}"};

  Configurable<float> matchedDetJetEtaMin{"matchedDetJetEtaMin", -0.5, "minimum matchedDetJet eta"};
  Configurable<float> matchedDetJetEtaMax{"matchedDetJetEtaMax", 0.5, "maximum matchedDetJet eta"};

  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {40, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {20, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binZ{"binZ", {40, 0.f, 1.f}, ""};
  ConfigurableAxis binXi{"binXi", {50, 0.f, 10.f}, ""};
  ConfigurableAxis binTheta{"binTheta", {41, -0.05f, 0.405f}, ""};
  ConfigurableAxis binJetR{"binJetR", {6, 0.05f, 0.65f}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binPDG{"binPDG", {static_cast<double>(pdgVector.size()), -0.5f, static_cast<double>(pdgVector.size()) - 0.5f}, ""};
  ConfigurableAxis binVtxZ{"binVtxZ", {200, -20, 20}, ""};

  ConfigurableAxis binPtDiff{"binPtDiff", {121, -20.5f, 100.5f}, ""};
  ConfigurableAxis binEtaDiff{"binEtaDiff", {41, -0.205f, 0.205f}, ""};
  ConfigurableAxis binPhiDiff{"binPhiDiff", {41, -0.205f, 0.205f}, ""};
  ConfigurableAxis binZDiff{"binZDiff", {21, -1.05f, 1.05f}, ""};
  ConfigurableAxis binXiDiff{"binXiDiff", {201, -10.5f, 10.5f}, ""};
  ConfigurableAxis binThetaDiff{"binThetaDiff", {201, -10.5f, 10.5f}, ""};
  ConfigurableAxis binPtRatio{"binPtRatio", {100, -0.5f, 9.5f}, ""};  // Ratio of pt, eta, phi
  ConfigurableAxis binMatchDist{"binMatchDist", {50, 0.f, 0.5f}, ""}; // Distance between matched jets

  ConfigurableAxis binCount{"binCount", {1, .5f, 1.5f}, ""};
  ConfigurableAxis jetCount{"jetCount", {20, -.5f, 19.5f}, ""};
  ConfigurableAxis trackCount{"trackCount", {1000, -.5f, 999.5f}, ""};

  Preslice<McPJets> perMcPJet = aod::jet::mcCollisionId;
  // Filter matchedDetJetFilter = (aod::chargedmcdetectorleveljets::eta >= matchedDetJetEtaMin && aod::chargedmcdetectorleveljets::eta <= matchedDetJetEtaMax);

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
    AxisSpec ptDiffAxis = {binPtDiff, "#it{p}_{T}^{jet, part} - #it{p}_{T}^{jet, det}"};
    AxisSpec etaDiffAxis = {binEtaDiff, "#eta^{jet, part} - #eta^{jet, det}"};
    AxisSpec phiDiffAxis = {binPhiDiff, "#varphi^{jet, part} - #varphi^{jet, det}"};
    AxisSpec zDiffAxis = {binZDiff, "#it{z}^{part} - #it{z}^{det}"};
    AxisSpec xiDiffAxis = {binXiDiff, "#xi^{part} - #xi^{det}"};
    AxisSpec thetaDiffAxis = {binThetaDiff, "#theta^{part} - #theta^{det}"};
    AxisSpec ptRatioAxis = {binPtRatio, ""};
    AxisSpec vtxZAxis = {binVtxZ, "Collision vertex z (cm)"};
    AxisSpec matchDistAxis = {binMatchDist, "#Delta"};

    if (doprocessDataRun3) {
      registry.add("data/nJetsnTracks", "nJetsnTracks; nJets; nTracks", HistType::kTH2F, {jetCount, trackCount});
      registry.add("data/collision/collisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});
      registry.add("data/tracks/trackPtEtaPhi", "trackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
      registry.add("data/jets/jetPtEtaPhi", "Jet #it{p}_{T}, #eta, #phi", HistType::kTH3F, {jetPtAxis, etaAxis, phiAxis});
      registry.add("data/jets/jetPtTrackPt", "Jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2F, {jetPtAxis, trackPtAxis});
      registry.add("data/jets/jetTrackPtEtaPhi", "Tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
      registry.add("data/jets/jetPtFrag", "Jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {jetPtAxis, zAxis});
      registry.add("data/jets/jetPtTrackProj", "Jet #it{p}_{T}, #it{z}", HistType::kTH2F, {jetPtAxis, zAxis});
      registry.add("data/jets/jetPtXi", "Jet #it{p}_{T}, #xi", HistType::kTH2F, {jetPtAxis, xiAxis});
      registry.add("data/jets/jetPtTheta", "Jet #it{p}_{T}, #theta", HistType::kTH2F, {jetPtAxis, thetaAxis});
      registry.add("data/jets/jetPtXiTheta", "Jet #it{p}_{T}, #xi, #theta", HistType::kTH3F, {jetPtAxis, xiAxis, thetaAxis});
      registry.add("data/jets/jetPtZTheta", "Jet #it{p}_{T}, z, #theta", HistType::kTH3F, {jetPtAxis, zAxis, thetaAxis});
    }

    if (doprocessMcP) {
      registry.add("particle-level/nJetsnTracks", "nJetsnTracks; nJets; nTracks", HistType::kTH2F, {jetCount, trackCount});
      registry.add("particle-level/collision/partCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});
      registry.add("particle-level/tracks/partTrackPtEtaPhi", "partTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
      registry.add("particle-level/jets/partJetPtEtaPhi", "Particle level jet #it{p}_{T}, #eta, #phi", HistType::kTH3F, {partJetPtAxis, partEtaAxis, partPhiAxis});
      registry.add("particle-level/jets/partJetPtTrackPt", "Particle level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2F, {partJetPtAxis, trackPtAxis});
      registry.add("particle-level/jets/partJetTrackPtEtaPhi", "Particle level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3F, {trackPtAxis, partEtaAxis, partPhiAxis});
      registry.add("particle-level/jets/partJetPtFrag", "Particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {partJetPtAxis, partZAxis});
      registry.add("particle-level/jets/partJetPtTrackProj", "Particle level jet #it{p}_{T}, #it{z}", HistType::kTH2F, {partJetPtAxis, partZAxis});
      registry.add("particle-level/jets/partJetPtXi", "Particle level jet #it{p}_{T}, #xi", HistType::kTH2F, {partJetPtAxis, partXiAxis});
      registry.add("particle-level/jets/partJetPtTheta", "Particle level jet #it{p}_{T}, #theta", HistType::kTH2F, {partJetPtAxis, partThetaAxis});
      registry.add("particle-level/jets/partJetPtXiTheta", "Particle level jet #it{p}_{T}, #xi, #theta", HistType::kTH3F, {partJetPtAxis, partXiAxis, partThetaAxis});
      registry.add("particle-level/jets/partJetPtZTheta", "Particle level jet #it{p}_{T}, z, #theta", HistType::kTH3F, {partJetPtAxis, partZAxis, partThetaAxis});
    }

    if (doprocessMcD) {
      registry.add("detector-level/nJetsnTracks", "nJetsnTracks; nJets; nTracks", HistType::kTH2F, {jetCount, trackCount});
      registry.add("detector-level/collision/detCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});
      registry.add("detector-level/tracks/detTrackPtEtaPhi", "detTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
      registry.add("detector-level/jets/detJetPtEtaPhi", "Detector level jet #it{p}_{T}, #eta, #phi", HistType::kTH3F, {detJetPtAxis, detEtaAxis, detPhiAxis});
      registry.add("detector-level/jets/detJetPtTrackPt", "Detector level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2F, {detJetPtAxis, trackPtAxis});
      registry.add("detector-level/jets/detJetTrackPtEtaPhi", "Detector level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3F, {trackPtAxis, detEtaAxis, detPhiAxis});
      registry.add("detector-level/jets/detJetPtFrag", "Detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {detJetPtAxis, detZAxis});
      registry.add("detector-level/jets/detJetPtTrackProj", "Detector level jet #it{p}_{T}, #it{z}", HistType::kTH2F, {detJetPtAxis, detZAxis});
      registry.add("detector-level/jets/detJetPtXi", "Detector level jet #it{p}_{T}, #xi", HistType::kTH2F, {detJetPtAxis, detXiAxis});
      registry.add("detector-level/jets/detJetPtTheta", "Detector level jet #it{p}_{T}, #theta", HistType::kTH2F, {detJetPtAxis, detThetaAxis});
      registry.add("detector-level/jets/detJetPtXiTheta", "Detector level jet #it{p}_{T}, #xi, #theta", HistType::kTH3F, {detJetPtAxis, detXiAxis, detThetaAxis});
      registry.add("detector-level/jets/detJetPtZTheta", "Detector level jet #it{p}_{T}, z, #theta", HistType::kTH3F, {detJetPtAxis, detZAxis, detThetaAxis});
    }

    if (doprocessMcMatched) {
      registry.add("matching/collision/matchCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});
      registry.add("matching/tracks/matchDetTrackPtEtaPhi", "matchDetTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
      registry.add("matching/tracks/matchPartTrackPtEtaPhi", "matchPartTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
      registry.add("matching/tracks/matchDetTrackPtPartTrackPt", "matchDetTrackPtPartTrackPt", HistType::kTH2F, {trackPtAxis, trackPtAxis});
      registry.add("matching/tracks/matchDetTrackEtaPartTrackEta", "matchDetTrackEtaPartTrackEta", HistType::kTH2F, {etaAxis, etaAxis});
      registry.add("matching/tracks/matchDetTrackPhiPartTrackPhi", "matchDetTrackPhiPartTrackPhi", HistType::kTH2F, {phiAxis, phiAxis});
      registry.add("matching/tracks/trackResolutionPt", "trackResolutionPt", HistType::kTH2F, {trackPtAxis, ptDiffAxis});
      registry.add("matching/tracks/trackResolutionEta", "trackResolutionEta", HistType::kTH2F, {etaAxis, etaDiffAxis});
      registry.add("matching/tracks/trackResolutionPhi", "trackResolutionPhi", HistType::kTH2F, {phiAxis, phiDiffAxis});
      // Detector level jets with a match
      registry.add("matching/jets/matchDetJetPtEtaPhi", "Matched detector level jet #it{p}_{T}, #eta, #phi", HistType::kTH3F, {detJetPtAxis, detEtaAxis, detPhiAxis});
      registry.add("matching/jets/matchDetJetPtTrackPt", "Matched detector level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2F, {detJetPtAxis, trackPtAxis});
      registry.add("matching/jets/matchDetJetTrackPtEtaPhi", "Matched detector level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3F, {trackPtAxis, detEtaAxis, detPhiAxis});
      registry.add("matching/jets/matchDetJetPtFrag", "Matched detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/matchDetJetPtTrackProj", "Matched detector level jet #it{p}_{T}, #it{z}", HistType::kTH2F, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/matchDetJetPtXi", "Matched detector level jet #it{p}_{T}, #xi", HistType::kTH2F, {detJetPtAxis, detXiAxis});
      registry.add("matching/jets/matchDetJetPtTheta", "Matched detector level jet #it{p}_{T}, #theta", HistType::kTH2F, {detJetPtAxis, detThetaAxis});
      registry.add("matching/jets/matchDetJetPtXiTheta", "Matched detector level jet #it{p}_{T}, #xi, #theta", HistType::kTH3F, {detJetPtAxis, detXiAxis, detThetaAxis});
      registry.add("matching/jets/matchDetJetPtZTheta", "Matched detector level jet #it{p}_{T}, z, #theta", HistType::kTH3F, {detJetPtAxis, detZAxis, detThetaAxis});
      // Particle level jets with a match
      registry.add("matching/jets/matchPartJetPtEtaPhi", "Matched particle level jet #it{p}_{T}, #eta, #phi", HistType::kTH3F, {partJetPtAxis, partEtaAxis, partPhiAxis});
      registry.add("matching/jets/matchPartJetPtTrackPt", "Matched particle level jet #it{p}_{T}, track #it{p}_{T}", HistType::kTH2F, {partJetPtAxis, trackPtAxis});
      registry.add("matching/jets/matchPartJetTrackPtEtaPhi", "Matched particle level tracks in jets #it{p}_{T}, #eta, #phi", HistType::kTH3F, {trackPtAxis, partEtaAxis, partPhiAxis});
      registry.add("matching/jets/matchPartJetPtFrag", "Matched particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/matchPartJetPtTrackProj", "Matched particle level jet #it{p}_{T}, #it{z}", HistType::kTH2F, {partJetPtAxis, partZAxis});
      registry.add("matching/jets/matchPartJetPtXi", "Matched particle level jet #it{p}_{T}, #xi", HistType::kTH2F, {partJetPtAxis, partXiAxis});
      registry.add("matching/jets/matchPartJetPtTheta", "Matched particle level jet #it{p}_{T}, #theta", HistType::kTH2F, {partJetPtAxis, partThetaAxis});
      registry.add("matching/jets/matchPartJetPtXiTheta", "Matched particle level jet #it{p}_{T}, #xi, #theta", HistType::kTH3F, {partJetPtAxis, partXiAxis, partThetaAxis});
      registry.add("matching/jets/matchPartJetPtZTheta", "Matched particle level jet #it{p}_{T}, z, #theta", HistType::kTH3F, {partJetPtAxis, partZAxis, partThetaAxis});
      // Combined information of matched jets
      registry.add("matching/jets/matchDetJetPtPartJetPt", "matchDetJetPtPartJetPt", HistType::kTH2F, {detJetPtAxis, partJetPtAxis});
      registry.add("matching/jets/matchPartJetPtDetJetEtaPartJetEta", "matchPartJetPtDetJetEtaPartJetEta", HistType::kTH3F, {partJetPtAxis, detEtaAxis, partEtaAxis});
      registry.add("matching/jets/matchPartJetPtDetJetPhiPartJetPhi", "matchPartJetPtDetJetPhiPartJetPhi", HistType::kTH3F, {partJetPtAxis, detPhiAxis, partPhiAxis});
      registry.add("matching/jets/matchPartJetPtResolutionPt", "#it{p}_{T}^{jet, part} - #it{p}_{T}^{jet, det}", HistType::kTH2F, {partJetPtAxis, ptDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionEta", "#eta^{jet, part} - #eta^{jet, det}", HistType::kTH3F, {partJetPtAxis, partEtaAxis, etaDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionPhi", "#phi^{jet, part} - #phi^{jet, det}", HistType::kTH3F, {partJetPtAxis, partPhiAxis, phiDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionChargeFrag", "Resolution #it{p}_{T}^{tr} / #it{p}_{T}^{jet}", HistType::kTH3F, {partJetPtAxis, partZAxis, zDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionTrackProj", "Resolution #it{p}^{proj} / #it{p}^{jet}", HistType::kTH3F, {partJetPtAxis, partZAxis, zDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionXi", "Resolution ln(1/#it{z})", HistType::kTH3F, {partJetPtAxis, partXiAxis, xiDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionTheta", "Resolution #theta", HistType::kTH3F, {partJetPtAxis, partThetaAxis, thetaDiffAxis});
      // TODO: does it make sense to add 2D resolution plots?
      registry.add("matching/jets/matchPartJetPtResolutionXiResolutionTheta", "Resolution #xi, #theta", HistType::kTHnSparseF, {partJetPtAxis, partXiAxis, xiDiffAxis, partThetaAxis, thetaDiffAxis});
      registry.add("matching/jets/matchPartJetPtResolutionZResolutionTheta", "Resolution #it{z}, #theta", HistType::kTHnSparseF, {partJetPtAxis, partZAxis, zDiffAxis, partThetaAxis, thetaDiffAxis});
      registry.add("matching/jets/matchPartJetPtEtaPhiMatchDist", "matchJetMatchDist", HistType::kTHnSparseF, {partJetPtAxis, partEtaAxis, partPhiAxis, matchDistAxis});
      registry.add("matching/jets/matchPartJetPtEnergyScale", "jetEnergyScale", HistType::kTH2F, {partJetPtAxis, ptRatioAxis});
      // QA histograms for fakes, misses
      registry.add("matching/jets/fakeDetJetPtEtaPhi", "Fakes", HistType::kTH3F, {detJetPtAxis, detEtaAxis, detPhiAxis});
      registry.add("matching/jets/missPartJetPtEtaPhi", "Misses", HistType::kTH3F, {partJetPtAxis, partEtaAxis, partPhiAxis});
      // Response matrix, fakes, misses
      registry.add("matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj", "Matched", HistType::kTHnSparseF, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
      registry.add("matching/jets/fakeDetJetPtTrackProj", "Fakes", HistType::kTH2F, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/missPartJetPtTrackProj", "Misses", HistType::kTH2F, {partJetPtAxis, partZAxis});

      registry.add("matching/jets/matchDetJetPtXiPartJetPtXi", "Matched", HistType::kTHnSparseF, {detJetPtAxis, detXiAxis, partJetPtAxis, partXiAxis});
      registry.add("matching/jets/fakeDetJetPtXi", "Fakes", HistType::kTH2F, {detJetPtAxis, detXiAxis});
      registry.add("matching/jets/missPartJetPtXi", "Misses", HistType::kTH2F, {partJetPtAxis, partXiAxis});

      registry.add("matching/jets/matchDetJetPtFragPartJetPtFrag", "Matched", HistType::kTHnSparseF, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
      registry.add("matching/jets/fakeDetJetPtFrag", "Fakes", HistType::kTH2F, {detJetPtAxis, detZAxis});
      registry.add("matching/jets/missPartJetPtFrag", "Misses", HistType::kTH2F, {partJetPtAxis, partZAxis});

      registry.add("matching/jets/matchDetJetPtThetaPartJetPtTheta", "Matched", HistType::kTHnSparseF, {detJetPtAxis, detThetaAxis, partJetPtAxis, partThetaAxis});
      registry.add("matching/jets/fakeDetJetPtTheta", "Fakes", HistType::kTH2F, {detJetPtAxis, detThetaAxis});
      registry.add("matching/jets/missPartJetPtTheta", "Misses", HistType::kTH2F, {partJetPtAxis, partThetaAxis});

      // TODO: Maybe need different axes here to have less granularity. Histogram may become too large.
      registry.add("matching/jets/matchDetJetPtXiThetaPartJetPtXiTheta", "Matched", HistType::kTHnSparseF, {detJetPtAxis, detXiAxis, detThetaAxis, partJetPtAxis, partXiAxis, partThetaAxis});
      registry.add("matching/jets/fakeDetJetPtXiTheta", "Fakes", HistType::kTH3F, {detJetPtAxis, detXiAxis, detThetaAxis});
      registry.add("matching/jets/missPartJetPtXiTheta", "Misses", HistType::kTH3F, {partJetPtAxis, partXiAxis, partThetaAxis});

      registry.add("matching/jets/matchDetJetPtZThetaPartJetPtZTheta", "Matched", HistType::kTHnSparseF, {detJetPtAxis, detZAxis, detThetaAxis, partJetPtAxis, partZAxis, partThetaAxis});
      registry.add("matching/jets/fakeDetJetPtZTheta", "Fakes", HistType::kTH3F, {detJetPtAxis, detZAxis, detThetaAxis});
      registry.add("matching/jets/missPartJetPtZTheta", "Misses", HistType::kTH3F, {partJetPtAxis, partZAxis, partThetaAxis});
    }
  } // init

  double CheckDphi(double dphi)
  {
    if (dphi > TMath::Pi()) {
      return (dphi - 2 * TMath::Pi());
    } else if (dphi < -1 * TMath::Pi()) {
      return (dphi + 2 * TMath::Pi());
    }
    return dphi;
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
    double dEta = -1., dPhi = -1., theta = -1.;
    dEta = constituent.eta() - jet.eta();
    dPhi = constituent.phi() - jet.phi();
    dPhi = CheckDphi(dPhi);
    theta = TMath::Sqrt(dEta * dEta + dPhi * dPhi);
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

  template <typename T>
  void fillDataRun3Histograms(T const& jet)
  {
    registry.fill(HIST("data/jets/jetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
    for (const auto& track : jet.template tracks_as<aod::JTracks>()) {
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

  template <typename DetJet, typename PartJet>
  void fillMatchingHistogramsJet(DetJet const& detJet, PartJet const& partJet, double weight = 1.)
  {
    double deltaEta = partJet.eta() - detJet.eta();
    double deltaPhi = partJet.phi() - detJet.phi();
    deltaPhi = CheckDphi(deltaPhi);
    double deltaR = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

    registry.fill(HIST("matching/jets/matchDetJetPtEtaPhi"), detJet.pt(), detJet.eta(), detJet.phi(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtEtaPhi"), partJet.pt(), partJet.eta(), partJet.phi(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtEtaPhiMatchDist"), partJet.pt(), partJet.eta(), partJet.phi(), deltaR, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtEnergyScale"), partJet.pt(), detJet.pt() / partJet.pt(), weight);
    registry.fill(HIST("matching/jets/matchDetJetPtPartJetPt"), detJet.pt(), partJet.pt(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtDetJetEtaPartJetEta"), partJet.pt(), detJet.eta(), partJet.eta(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtDetJetPhiPartJetPhi"), partJet.pt(), detJet.phi(), partJet.phi(), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionPt"), partJet.pt(), (partJet.pt() - detJet.pt()), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionEta"), partJet.pt(), partJet.eta(), deltaEta, weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionPhi"), partJet.pt(), partJet.phi(), deltaPhi, weight);
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

    partChargeFrag = ChargeFrag(partJet, track);
    partTrackProj = TrackProj(partJet, track);
    partTheta = Theta(partJet, track);
    partXi = Xi(partJet, track);

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
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionChargeFrag"), partJet.pt(), partChargeFrag, (partChargeFrag - detChargeFrag), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionTrackProj"), partJet.pt(), partTrackProj, (partTrackProj - detTrackProj), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionXi"), partJet.pt(), partXi, (partXi - detXi), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionTheta"), partJet.pt(), partTheta, (partTheta - detTheta), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionXiResolutionTheta"), partJet.pt(), partXi, (partXi - detXi), partTheta, (partTheta - detTheta), weight);
    registry.fill(HIST("matching/jets/matchPartJetPtResolutionZResolutionTheta"), partJet.pt(), partTrackProj, (partTrackProj - detTrackProj), partTheta, (partTheta - detTheta), weight);

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
    for (const auto& track : jet.template tracks_as<aod::JTracks>()) {
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
    for (const auto& track : jet.template tracks_as<aod::JMcParticles>()) {
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

  void processDummy(aod::JTracks const& track) {}
  PROCESS_SWITCH(JetFragmentation, processDummy, "Dummy process function turned on by default", true);

  void processMcD(soa::Join<aod::JCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                  aod::JMcCollisions const& mcCollisions,
                  McDJets const& jets,
                  aod::JTracks const& tracks)
  {
    float nJets = 0, nTracks = 0;
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

  void processMcP(aod::JMcCollision const& mcCollision, // Add some form of event selection?
                  McPJets const& jets,
                  aod::JMcParticles const& particles)
  {
    float nJets = 0, nTracks = 0;
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

  void processDataRun3(aod::JCollision const& collision,
                       soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                       aod::JTracks const& tracks)
  {
    float nJets = 0, nTracks = 0;
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

  void processMcMatched(soa::Join<aod::JCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                        soa::Join<McDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& mcDetJets,
                        McTracks const& tracks,
                        aod::JMcCollisions const& mcCollisions,
                        McPJets const& mcPartJets,
                        aod::JMcParticles const& mcParticles)
  {
    double weight = collision.mcCollision().weight();
    bool isFake = false;
    for (const auto& detJet : mcDetJets) {
      if (detJet.eta() < matchedDetJetEtaMin || detJet.eta() > matchedDetJetEtaMax) {
        continue; // TODO: should be done in filter
      }
      for (auto& partJet : detJet.template matchedJetGeo_as<McPJets>()) {
        fillMatchingHistogramsJet(detJet, partJet, weight);

        for (const auto& track : detJet.tracks_as<McTracks>()) {
          bool isTrackMatched = false;
          for (const auto& particle : partJet.tracks_as<aod::JMcParticles>()) {
            if (track.has_mcParticle() && particle.globalIndex() == track.template mcParticle_as<aod::JMcParticles>().globalIndex()) {
              isTrackMatched = true;
              fillMatchingHistogramsConstituent(detJet, partJet, track, particle);
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
        for (const auto& track : detJet.tracks_as<McTracks>()) {
          fillMatchingFakeOrMiss(detJet, track, isFake, weight);
        }
      } // if detJet does not have a match
    }   // for det jet
    // TODO: how to deal with misses?
    //       -> ParticleToDetector table is currently bugged: size does not correspond to ParticleLevelJets (04.06.2023)
  }
  PROCESS_SWITCH(JetFragmentation, processMcMatched, "Monte Carlo particle and detector level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetFragmentation>(cfgc, TaskName{"jet-fragmentation"})};
}
