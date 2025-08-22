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

/// \file v0QA.cxx
/// \brief QA task for V0s in the jets framework, based on the LF v0cascadesqa task
///
/// \author Gijs van Weelden <g.van.weelden@cern.ch>

#include "JetDerivedDataUtilities.h"

#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/V0SelectorTables.h"

#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// V0 jets
using MCDV0Jets = aod::V0ChargedMCDetectorLevelJets;
using MCDV0JetsWithConstituents = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetConstituents>;
using MatchedMCDV0Jets = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;
using MatchedMCDV0JetsWithConstituents = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetConstituents, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;

using CandidatesV0MCDWithFlags = soa::Join<aod::CandidatesV0MCD, aod::McV0Labels>;

using MCPV0Jets = aod::V0ChargedMCParticleLevelJets;
using MCPV0JetsWithConstituents = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetConstituents>;
using MatchedMCPV0Jets = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;
using MatchedMCPV0JetsWithConstituents = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetConstituents, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;

using JetMcCollisionsWithPIs = soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>;

struct V0QA {
  HistogramRegistry registry{"registry"};

  Configurable<std::string> evSel{"evSel", "sel8WithoutTimeFrameBorderCut", "choose event selection"};
  Configurable<float> yPartMax{"yPartMax", 0.5, "Maximum rapidity of particles"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "Vertex Z cut"};
  Configurable<float> v0Fraction{"v0Fraction", 1.0, "Fraction of V0s to be kept inside jets"};

  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  ConfigurableAxis binPtJet{"binPtJet", {100., 0.0f, 50.0f}, ""};
  ConfigurableAxis binPtV0{"binPtV0", {100., 0.0f, 50.0f}, ""};
  ConfigurableAxis binZV0{"binZV0", {100., 1e-3f, 1 + 1e-3f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.0f, 1.0f}, ""};
  ConfigurableAxis binPhi{"binPhi", {constants::math::PI * 10 / 2, 0.0f, constants::math::TwoPI}, ""};

  ConfigurableAxis binInvMassK0S{"binInvMassK0S", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis binInvMassLambda{"binInvMassLambda", {200, 1.07f, 1.17f}, ""};
  ConfigurableAxis binV0Radius{"binV0Radius", {100., 0.0f, 50.0f}, ""};
  ConfigurableAxis binV0CosPA{"binV0CosPA", {50., 0.95f, 1.0f}, ""};

  ConfigurableAxis binsDcaXY{"binsDcaXY", {100, -0.5f, 0.5f}, ""};
  ConfigurableAxis binsDcaZ{"binsDcaZ", {100, -5.f, 5.f}, ""};
  ConfigurableAxis binPtDiff{"binPtDiff", {200., -49.5f, 50.5f}, ""};
  ConfigurableAxis binPtRelDiff{"binPtRelDiff", {100., -1.0f, 1.0f}, ""};
  ConfigurableAxis binITSNCl{"binITSNCl", {8, -0.5, 7.5}, ""};
  ConfigurableAxis binITSChi2NCl{"binITSChi2NCl", {100, 0, 40}, ""};

  ConfigurableAxis binTPCNCl{"binTPCNCl", {165, -0.5, 164.5}, ""};
  ConfigurableAxis binTPCChi2NCl{"binTPCChi2NCl", {100, 0, 10}, ""};
  ConfigurableAxis binTPCNClSharedFraction{"binTPCNClSharedFraction", {100, 0., 1.}, ""};
  ConfigurableAxis binTPCCrossedRowsOverFindableCl{"binTPCCrossedRowsOverFindableCl", {120, 0.0, 1.2}, ""};

  std::vector<int> eventSelectionBits;

  void init(InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(evSel));

    const AxisSpec axisJetPt{binPtJet, "Jet Pt (GeV/c)"};
    const AxisSpec axisV0Pt{binPtV0, "V0 Pt (GeV/c)"};
    const AxisSpec axisV0Z{binZV0, "z_{V0} = #it{p}_{T, V0} / #it{p}_{T, jet}"};
    const AxisSpec axisEta{binEta, "Eta"};
    const AxisSpec axisPhi{binPhi, "Phi"};
    const AxisSpec axisV0Radius{binV0Radius, "V0 Radius (cm)"};
    const AxisSpec axisV0CosPA{binV0CosPA, "V0 CosPA"};
    const AxisSpec axisK0SM{binInvMassK0S, "M(#pi^{+} #pi^{-}) (GeV/c^{2})"};
    const AxisSpec axisLambdaM{binInvMassLambda, "M(p #pi^{-}) (GeV/c^{2})"};
    const AxisSpec axisAntiLambdaM{binInvMassLambda, "M(#bar{p} #pi^{+}) (GeV/c^{2})"};

    const AxisSpec axisPtDiff{binPtDiff, "Pt difference (GeV/c)"};
    const AxisSpec axisPtRelDiff{binPtRelDiff, "Pt relative difference"};
    const AxisSpec axisDcaXY{binsDcaXY, "DCA_{xy} (cm)"};
    const AxisSpec axisDcaZ{binsDcaZ, "DCA_{z} (cm)"};
    const AxisSpec axisITSNCl{binITSNCl, "# clusters ITS"};
    const AxisSpec axisITSChi2NCl{binITSChi2NCl, "Chi2 / cluster ITS"};

    const AxisSpec axisNClFindable{binTPCNCl, "# findable clusters TPC"};
    const AxisSpec axisNClFound{binTPCNCl, "# found clusters TPC"};
    const AxisSpec axisNClShared{binTPCNCl, "# shared clusters TPC"};
    const AxisSpec axisNClCrossedRows{binTPCNCl, "# crossed rows TPC"};
    const AxisSpec axisTPCChi2NCl{binTPCChi2NCl, "Chi2 / cluster TPC"};
    const AxisSpec axisSharedFraction{binTPCNClSharedFraction, "Fraction shared clusters TPC"};
    const AxisSpec axisCrossedRowsOverFindable{binTPCCrossedRowsOverFindableCl, "Crossed rows / findable clusters TPC"};

    const bool doSumw2 = true;

    if (doprocessFlags) {
      registry.add("inclusive/V0Flags", "V0Flags", HistType::kTH2D, {{5, -0.5, 4.5}, {5, -0.5, 4.5}});
    }
    if (doprocessMcD) {
      registry.add("inclusive/hEvents", "Events", {HistType::kTH1D, {{3, 0.0f, 3.0f}}}, doSumw2);

      registry.add("inclusive/K0SPtEtaMass", "K0S Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisK0SM}, doSumw2);
      registry.add("inclusive/InvMassK0STrue", "Invariant mass of K0S", HistType::kTH3D, {axisV0Pt, axisV0Radius, axisK0SM}, doSumw2);
      registry.add("inclusive/K0SPtEtaMassWrongCollision", "K0S Pt, Eta, Mass (Wrong Collision)", HistType::kTH3D, {axisV0Pt, axisEta, axisK0SM}, doSumw2);
      registry.add("inclusive/LambdaPtEtaMass", "Lambda Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisLambdaM}, doSumw2);
      registry.add("inclusive/InvMassLambdaTrue", "Invariant mass of Lambda", HistType::kTH3D, {axisV0Pt, axisV0Radius, axisLambdaM}, doSumw2);
      registry.add("inclusive/LambdaPtEtaMassWrongCollision", "Lambda Pt, Eta, Mass (Wrong Collision)", HistType::kTH3D, {axisV0Pt, axisEta, axisLambdaM}, doSumw2);
      registry.add("inclusive/AntiLambdaPtEtaMass", "AntiLambda Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisAntiLambdaM}, doSumw2);
      registry.add("inclusive/InvMassAntiLambdaTrue", "Invariant mass of AntiLambda", HistType::kTH3D, {axisV0Pt, axisV0Radius, axisAntiLambdaM}, doSumw2);
      registry.add("inclusive/AntiLambdaPtEtaMassWrongCollision", "AntiLambda Pt, Eta, Mass (Wrong Collision)", HistType::kTH3D, {axisV0Pt, axisEta, axisAntiLambdaM}, doSumw2);

      registry.add("jets/JetPtEtaPhi", "Jet Pt, Eta, Phi", HistType::kTH3D, {axisJetPt, axisEta, axisPhi}, doSumw2);
      registry.add("jets/JetPtEtaK0SPt", "Jet Pt, Eta, K0S Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetPtEtaK0SZ", "Jet Pt, Eta, K0S Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetPtEtaK0SPtWrongCollision", "Jet Pt, Eta, K0S Pt (Wrong Collision)", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetPtEtaK0SZWrongCollision", "Jet Pt, Eta, K0S Z (Wrong Collision)", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetPtEtaLambdaPt", "Jet Pt, Eta, Lambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetPtEtaLambdaZ", "Jet Pt, Eta, Lambda Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetPtEtaLambdaPtWrongCollision", "Jet Pt, Eta, Lambda Pt (Wrong Collision)", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetPtEtaLambdaZWrongCollision", "Jet Pt, Eta, Lambda Z (Wrong Collision)", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetPtEtaAntiLambdaPt", "Jet Pt, Eta, AntiLambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetPtEtaAntiLambdaZ", "Jet Pt, Eta, AntiLambda Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetPtEtaAntiLambdaPtWrongCollision", "Jet Pt, Eta, AntiLambda Pt (Wrong Collision)", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetPtEtaAntiLambdaZWrongCollision", "Jet Pt, Eta, AntiLambda Z (Wrong Collision)", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);

      registry.add("jets/JetsPtEta", "Matched Jet Pt, Eta", HistType::kTH3D, {axisJetPt, axisJetPt, axisEta}, doSumw2);
      registry.add("jets/JetsPtEtaK0SPt", "Matched Jet Pt, Eta, K0S Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetsPtEtaK0SZ", "Matched Jet Pt, Eta, K0S Z", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetsPtEtaK0SPtWrongCollision", "Matched Jet Pt, Eta, K0S Pt (Wrong Collision)", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetsPtEtaK0SZWrongCollision", "Matched Jet Pt, Eta, K0S Z (Wrong Collision)", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetsPtEtaLambdaPt", "Matched Jet Pt, Eta, Lambda Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetsPtEtaLambdaZ", "Matched Jet Pt, Eta, Lambda Z", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetsPtEtaLambdaPtWrongCollision", "Matched Jet Pt, Eta, Lambda Pt (Wrong Collision)", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetsPtEtaLambdaZWrongCollision", "Matched Jet Pt, Eta, Lambda Z (Wrong Collision)", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetsPtEtaAntiLambdaPt", "Matched Jet Pt, Eta, AntiLambda Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetsPtEtaAntiLambdaZ", "Matched Jet Pt, Eta, AntiLambda Z", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/JetsPtEtaAntiLambdaPtWrongCollision", "Matched Jet Pt, Eta, AntiLambda Pt (Wrong Collision)", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/JetsPtEtaAntiLambdaZWrongCollision", "Matched Jet Pt, Eta, AntiLambda Z (Wrong Collision)", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z}, doSumw2);
    }
    if (doprocessMcP) {
      registry.add("inclusive/hMcEvents", "MC Events", {HistType::kTH1D, {{3, 0.0f, 3.0f}}}, doSumw2);
      registry.add("inclusive/GeneratedK0S", "Generated K0S", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Radius}, doSumw2);
      registry.add("inclusive/GeneratedLambda", "Generated Lambda", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Radius}, doSumw2);
      registry.add("inclusive/GeneratedAntiLambda", "Generated AntiLambda", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Radius}, doSumw2);

      registry.add("jets/GeneratedJets", "Generated Jets", HistType::kTH3D, {axisJetPt, axisEta, axisPhi}, doSumw2);
      registry.add("jets/GeneratedJetK0S", "Generated Jet K0S", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/GeneratedJetK0SFrag", "Generated Jet K0S", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/GeneratedJetLambda", "Generated Jet Lambda", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/GeneratedJetLambdaFrag", "Generated Jet Lambda", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);
      registry.add("jets/GeneratedJetAntiLambda", "Generated Jet AntiLambda", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt}, doSumw2);
      registry.add("jets/GeneratedJetAntiLambdaFrag", "Generated Jet AntiLambda", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z}, doSumw2);
    }
    if (doprocessCollisionAssociation) {
      registry.add("collisions/V0PtEta", "V0 Pt, Eta", HistType::kTH2D, {axisV0Pt, axisEta});
      registry.add("collisions/V0PtEtaWrongColl", "V0 Pt, Eta, (Wrong Collision)", HistType::kTH2D, {axisV0Pt, axisEta});
      registry.add("collisions/K0SPtEtaMass", "K0S Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisK0SM});
      registry.add("collisions/K0SPtEtaMassWrongColl", "K0S Pt, Eta, Mass, (Wrong Collision)", HistType::kTH3D, {axisV0Pt, axisEta, axisK0SM});
      registry.add("collisions/LambdaPtEtaMass", "Lambda Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisLambdaM});
      registry.add("collisions/LambdaPtEtaMassWrongColl", "Lambda Pt, Eta, Mass, (Wrong Collision)", HistType::kTH3D, {axisV0Pt, axisEta, axisLambdaM});
      registry.add("collisions/AntiLambdaPtEtaMass", "AntiLambda Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisAntiLambdaM});
      registry.add("collisions/AntiLambdaPtEtaMassWrongColl", "AntiLambda Pt, Eta, Mass, (Wrong Collision)", HistType::kTH3D, {axisV0Pt, axisEta, axisAntiLambdaM});

      registry.add("collisions/XiMinusPtYLambdaPt", "#Xi^{-} Pt, Y, #Lambda Pt", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Pt});
      registry.add("collisions/XiMinusPtYLambdaPtWrongColl", "#Xi^{-} Pt, Y, #Lambda Pt, (Wrong Collision)", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Pt});
      registry.add("collisions/XiPlusPtYAntiLambdaPt", "#Xi^{+} Pt, Y, #bar{#Lambda} Pt", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Pt});
      registry.add("collisions/XiPlusPtYAntiLambdaPtWrongColl", "#Xi^{+} Pt, Y, #bar{#Lambda} Pt, (Wrong Collision)", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Pt});
    }
    if (doprocessCollisionAssociationJets) {
      registry.add("collisions/JetPtEtaV0Pt", "Jet Pt, Eta, V0 Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("collisions/JetPtEtaV0PtWrongColl", "Jet Pt, Eta, V0 Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("collisions/JetPtEtaK0SPtMass", "Jet Pt, Eta, K0S Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisK0SM});
      registry.add("collisions/JetPtEtaK0SFragMass", "Jet Pt, Eta, K0S Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Z, axisK0SM});
      registry.add("collisions/JetPtEtaK0SPtMassWrongColl", "Jet Pt, Eta, K0S Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisK0SM});
      registry.add("collisions/JetPtEtaK0SFragMassWrongColl", "Jet Pt, Eta, K0S Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Z, axisK0SM});
      registry.add("collisions/JetPtEtaLambdaPtMass", "Jet Pt, Eta, Lambda Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisLambdaM});
      registry.add("collisions/JetPtEtaLambdaFragMass", "Jet Pt, Eta, Lambda Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Z, axisLambdaM});
      registry.add("collisions/JetPtEtaLambdaPtMassWrongColl", "Jet Pt, Eta, Lambda Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisLambdaM});
      registry.add("collisions/JetPtEtaLambdaFragMassWrongColl", "Jet Pt, Eta, Lambda Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Z, axisLambdaM});
      registry.add("collisions/JetPtEtaAntiLambdaPtMass", "Jet Pt, Eta, AntiLambda Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisAntiLambdaM});
      registry.add("collisions/JetPtEtaAntiLambdaFragMass", "Jet Pt, Eta, AntiLambda Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Z, axisAntiLambdaM});
      registry.add("collisions/JetPtEtaAntiLambdaPtMassWrongColl", "Jet Pt, Eta, AntiLambda Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisAntiLambdaM});
      registry.add("collisions/JetPtEtaAntiLambdaFragMassWrongColl", "Jet Pt, Eta, AntiLambda Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Z, axisAntiLambdaM});

      registry.add("collisions/JetPtEtaXiMinusPtLambdaPt", "Jet Pt, #Xi^{-} Pt, #Lambda Pt", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisV0Pt});
      registry.add("collisions/JetPtEtaXiMinusPtLambdaPtWrongColl", "Jet Pt, #Xi^{-} Pt, #Lambda Pt", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisV0Pt});
      registry.add("collisions/JetPtEtaXiPlusPtAntiLambdaPt", "Jet Pt, #Xi^{+} Pt, #bar{#Lambda} Pt", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisV0Pt});
      registry.add("collisions/JetPtEtaXiPlusPtAntiLambdaPtWrongColl", "Jet Pt, #Xi^{+} Pt, #bar{#Lambda} Pt", HistType::kTHnSparseD, {axisJetPt, axisEta, axisV0Pt, axisV0Pt});
    }
    if (doprocessCollisionAssociationMatchedJets) {
      registry.add("collisions/JetsPtEtaV0Pt", "Jets Pt, Eta, V0 Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt});
      registry.add("collisions/JetsPtEtaV0PtWrongColl", "Jets Pt, Eta, V0 Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt});
      registry.add("collisions/JetsPtEtaK0SPtMass", "Jets Pt, Eta, K0S Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisK0SM});
      registry.add("collisions/JetsPtEtaK0SFragMass", "Jets Pt, Eta, K0S Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z, axisK0SM});
      registry.add("collisions/JetsPtEtaK0SFragMassWrongColl", "Jets Pt, Eta, K0S Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z, axisK0SM});
      registry.add("collisions/JetsPtEtaLambdaPtMass", "Jets Pt, Eta, Lambda Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisLambdaM});
      registry.add("collisions/JetsPtEtaLambdaFragMass", "Jets Pt, Eta, Lambda Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z, axisLambdaM});
      registry.add("collisions/JetsPtEtaLambdaPtMassWrongColl", "Jets Pt, Eta, Lambda Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisLambdaM});
      registry.add("collisions/JetsPtEtaLambdaFragMassWrongColl", "Jets Pt, Eta, Lambda Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z, axisLambdaM});
      registry.add("collisions/JetsPtEtaAntiLambdaPtMass", "Jets Pt, Eta, AntiLambda Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisAntiLambdaM});
      registry.add("collisions/JetsPtEtaAntiLambdaFragMass", "Jets Pt, Eta, AntiLambda Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z, axisAntiLambdaM});
      registry.add("collisions/JetsPtEtaAntiLambdaPtMassWrongColl", "Jets Pt, Eta, AntiLambda Pt Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisAntiLambdaM});
      registry.add("collisions/JetsPtEtaAntiLambdaFragMassWrongColl", "Jets Pt, Eta, AntiLambda Frag Mass", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Z, axisAntiLambdaM});

      registry.add("collisions/JetsPtEtaXiMinusPtLambdaPt", "Jets Pt, Eta, #Xi^{-} Pt, #Lambda Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisV0Pt});
      registry.add("collisions/JetsPtEtaXiMinusPtLambdaPtWrongColl", "Jets Pt, Eta, #Xi^{-} Pt, #Lambda Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisV0Pt});
      registry.add("collisions/JetsPtEtaXiPlusPtAntiLambdaPt", "Jets Pt, Eta, #Xi^{+} Pt, #bar{#Lambda} Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisV0Pt});
      registry.add("collisions/JetsPtEtaXiPlusPtAntiLambdaPtWrongColl", "Jets Pt, Eta, #Xi^{+} Pt, #bar{#Lambda} Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt, axisV0Pt});
    }
    if (doprocessFeeddown) {
      registry.add("feeddown/XiMinusPtYLambdaPt", "#Xi^{-} Pt, Y, #Lambda Pt", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Pt});
      registry.add("feeddown/XiPlusPtYAntiLambdaPt", "#Xi^{-} Pt, Y, #Lambda Pt", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Pt});
    }
    if (doprocessFeeddownJets) {
      registry.add("feeddown/JetPtXiMinusPtLambdaPt", "Jets Pt, #Xi^{-} Pt, #Lambda Pt", HistType::kTH3D, {axisJetPt, axisJetPt, axisV0Pt});
      registry.add("feeddown/JetPtXiPlusPtAntiLambdaPt", "Jets Pt, #Xi^{+} Pt, #Lambda Pt", HistType::kTH3D, {axisJetPt, axisJetPt, axisV0Pt});
    }
    if (doprocessFeeddownMatchedJets) {
      registry.add("feeddown/JetsPtXiMinusPtLambdaPt", "Jets Pt, #Xi^{-} Pt, #Lambda Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisV0Pt, axisV0Pt});
      registry.add("feeddown/JetsPtXiPlusPtAntiLambdaPt", "Jets Pt, #Xi^{+} Pt, #bar{#Lambda} Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisV0Pt, axisV0Pt});
    }
    if (doprocessTestWeightedJetFinder) {
      registry.add("tests/weighted/hEvents", "Events", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("tests/weighted/JetPtEtaPhi", "Jet Pt, Eta, Phi", HistType::kTH3D, {axisJetPt, axisEta, axisPhi});
      registry.add("tests/weighted/JetPtEtaV0Pt", "Jet Pt, Eta, V0 Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/weighted/JetPtEtaV0Z", "Jet Pt, Eta, V0 Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/weighted/JetPtEtaK0SPt", "Jet Pt, Eta, K0S Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/weighted/JetPtEtaK0SZ", "Jet Pt, Eta, K0S Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/weighted/JetPtEtaLambdaPt", "Jet Pt, Eta, Lambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/weighted/JetPtEtaLambdaZ", "Jet Pt, Eta, Lambda Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/weighted/JetPtEtaAntiLambdaPt", "Jet Pt, Eta, AntiLambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/weighted/JetPtEtaAntiLambdaZ", "Jet Pt, Eta, AntiLambda Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
    }
    if (doprocessTestSubtractedJetFinder) {
      registry.add("tests/hEvents", "Events", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("tests/nosub/JetPtEtaPhi", "Jet Pt, Eta, Phi", HistType::kTH3D, {axisJetPt, axisEta, axisPhi});
      registry.add("tests/nosub/JetPtEtaV0Pt", "Jet Pt, Eta, V0 Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/nosub/JetPtEtaV0Z", "Jet Pt, Eta, V0 Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/nosub/JetPtEtaK0SPt", "Jet Pt, Eta, K0S Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/nosub/JetPtEtaK0SZ", "Jet Pt, Eta, K0S Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/nosub/JetPtEtaLambdaPt", "Jet Pt, Eta, Lambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/nosub/JetPtEtaLambdaZ", "Jet Pt, Eta, Lambda Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/nosub/JetPtEtaAntiLambdaPt", "Jet Pt, Eta, AntiLambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/nosub/JetPtEtaAntiLambdaZ", "Jet Pt, Eta, AntiLambda Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});

      registry.add("tests/sub/JetPtEtaPhi", "Jet Pt, Eta, Phi", HistType::kTH3D, {axisJetPt, axisEta, axisPhi});
      registry.add("tests/sub/JetPtEtaV0Pt", "Jet Pt, Eta, V0 Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/sub/JetPtEtaV0Z", "Jet Pt, Eta, V0 Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/sub/JetPtEtaK0SPt", "Jet Pt, Eta, K0S Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/sub/JetPtEtaK0SZ", "Jet Pt, Eta, K0S Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/sub/JetPtEtaLambdaPt", "Jet Pt, Eta, Lambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/sub/JetPtEtaLambdaZ", "Jet Pt, Eta, Lambda Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
      registry.add("tests/sub/JetPtEtaAntiLambdaPt", "Jet Pt, Eta, AntiLambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("tests/sub/JetPtEtaAntiLambdaZ", "Jet Pt, Eta, AntiLambda Z", HistType::kTH3D, {axisJetPt, axisEta, axisV0Z});
    }
    if (doprocessV0TrackQA) {
      registry.add("tracks/Pos", "pos", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisEta, axisPhi});
      registry.add("tracks/Neg", "neg", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisEta, axisPhi});
      registry.add("tracks/Pt", "pt", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff, axisPtRelDiff});
      registry.add("tracks/PtMass", "pt mass", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM, axisLambdaM, axisAntiLambdaM});
      registry.add("tracks/PtDiffMass", "ptdiff mass", HistType::kTHnSparseD, {axisV0Pt, axisPtDiff, axisK0SM, axisLambdaM, axisAntiLambdaM});

      registry.add("tracks/DCAxy", "dcaxy", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisDcaXY, axisDcaXY, axisPtDiff});
      registry.add("tracks/DCAxyMassK0S", "dcaxy mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisDcaXY, axisDcaXY, axisK0SM});
      registry.add("tracks/DCAxyMassLambda0", "dcaxy mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisDcaXY, axisDcaXY, axisLambdaM});
      registry.add("tracks/DCAxyMassAntiLambda0", "dcaxy mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisDcaXY, axisDcaXY, axisAntiLambdaM});

      registry.add("tracks/DCAz", "dcaz", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisDcaZ, axisDcaZ, axisPtDiff});
      registry.add("tracks/DCAzMassK0S", "dcaz mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisDcaZ, axisDcaZ, axisK0SM});
      registry.add("tracks/DCAzMassLambda0", "dcaz mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisDcaZ, axisDcaZ, axisLambdaM});
      registry.add("tracks/DCAzMassAntiLambda0", "dcaz mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisDcaZ, axisDcaZ, axisAntiLambdaM});

      registry.add("tracks/V0Radius", "v0 radius", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisV0Radius, axisPtDiff});
      registry.add("tracks/V0RadiusMassK0S", "v0 radius mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisV0Radius, axisK0SM});
      registry.add("tracks/V0RadiusMassLambda0", "v0 radius mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisV0Radius, axisLambdaM});
      registry.add("tracks/V0RadiusMassAntiLambda0", "v0 radius mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisV0Radius, axisAntiLambdaM});

      registry.add("tracks/V0CosPa", "v0 cos pa", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisV0CosPA, axisPtDiff});
      registry.add("tracks/V0CosPaMassK0S", "v0 cos pa mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisV0CosPA, axisK0SM});
      registry.add("tracks/V0CosPaMassLambda0", "v0 cos pa mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisV0CosPA, axisLambdaM});
      registry.add("tracks/V0CosPaMassAntiLambda0", "v0 cos pa mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisV0CosPA, axisAntiLambdaM});

      // TRD
      registry.add("tracks/posTRDPt", "pos trd pt", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/posTRDPtMass", "pos trd pt mass", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM, axisLambdaM, axisAntiLambdaM});
      registry.add("tracks/posNoTRDPt", "pos no trd pt", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/posNoTRDPtMass", "pos no trd pt mass", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM, axisLambdaM, axisAntiLambdaM});

      registry.add("tracks/negTRDPt", "neg trd pt", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/negTRDPtMass", "neg trd pt mass", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM, axisLambdaM, axisAntiLambdaM});
      registry.add("tracks/negNoTRDPt", "neg no trd pt", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/negNoTRDPtMass", "neg no trd pt mass", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM, axisLambdaM, axisAntiLambdaM});

      // ITS: positive track
      registry.add("tracks/ITS/posLayer1", "pos layer 1", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer1MassK0S", "pos layer 1 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer1MassLambda0", "pos layer 1 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer1MassAntiLambda0", "pos layer 1 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer2", "pos layer 2", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer2MassK0S", "pos layer 2 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer2MassLambda0", "pos layer 2 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer2MassAntiLambda0", "pos layer 2 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer3", "pos layer 3", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer3MassK0S", "pos layer 3 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer3MassLambda0", "pos layer 3 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer3MassAntiLambda0", "pos layer 3 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer4", "pos layer 4", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer4MassK0S", "pos layer 4 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer4MassLambda0", "pos layer 4 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer4MassAntiLambda0", "pos layer 4 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer5", "pos layer 5", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer5MassK0S", "pos layer 5 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer5MassLambda0", "pos layer 5 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer5MassAntiLambda0", "pos layer 5 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer6", "pos layer 6", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer6MassK0S", "pos layer 6 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer6MassLambda0", "pos layer 6 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer6MassAntiLambda0", "pos layer 6 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer7", "pos layer 7", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer7MassK0S", "pos layer 7 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer7MassLambda0", "pos layer 7 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer7MassAntiLambda0", "pos layer 7 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer56", "pos layer 56", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer56MassK0S", "pos layer 56 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer56MassLambda0", "pos layer 56 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer56MassAntiLambda0", "pos layer 56 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer67", "pos layer 67", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer67MassK0S", "pos layer 67 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer67MassLambda0", "pos layer 67 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer67MassAntiLambda0", "pos layer 67 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer57", "pos layer 57", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer57MassK0S", "pos layer 57 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer57MassLambda0", "pos layer 57 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer57MassAntiLambda0", "pos layer 57 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posLayer567", "pos layer 567", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/posLayer567MassK0S", "pos layer 567 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/posLayer567MassLambda0", "pos layer 567 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/posLayer567MassAntiLambda0", "pos layer 567 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/posNCl", "pos ncl", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisITSNCl});
      registry.add("tracks/ITS/posChi2NCl", "pos chi2ncl", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisITSChi2NCl});

      // ITS: Negative track
      registry.add("tracks/ITS/negLayer1", "neg layer 1", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer1MassK0S", "neg layer 1 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer1MassLambda0", "neg layer 1 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer1MassAntiLambda0", "neg layer 1 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer2", "neg layer 2", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer2MassK0S", "neg layer 2 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer2MassLambda0", "neg layer 2 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer2MassAntiLambda0", "neg layer 2 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer3", "neg layer 3", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer3MassK0S", "neg layer 3 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer3MassLambda0", "neg layer 3 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer3MassAntiLambda0", "neg layer 3 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer4", "neg layer 4", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer4MassK0S", "neg layer 4 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer4MassLambda0", "neg layer 4 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer4MassAntiLambda0", "neg layer 4 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer5", "neg layer 5", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer5MassK0S", "neg layer 5 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer5MassLambda0", "neg layer 5 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer5MassAntiLambda0", "neg layer 5 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer6", "neg layer 6", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer6MassK0S", "neg layer 6 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer6MassLambda0", "neg layer 6 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer6MassAntiLambda0", "neg layer 6 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer7", "neg layer 7", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer7MassK0S", "neg layer 7 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer7MassLambda0", "neg layer 7 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer7MassAntiLambda0", "neg layer 7 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer56", "neg layer 56", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer56MassK0S", "neg layer 56 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer56MassLambda0", "neg layer 56 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer56MassAntiLambda0", "neg layer 56 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer67", "neg layer 67", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer67MassK0S", "neg layer 67 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer67MassLambda0", "neg layer 67 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer67MassAntiLambda0", "neg layer 67 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer57", "neg layer 57", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer57MassK0S", "neg layer 57 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer57MassLambda0", "neg layer 57 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer57MassAntiLambda0", "neg layer 57 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negLayer567", "neg layer 567", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisPtDiff});
      registry.add("tracks/ITS/negLayer567MassK0S", "neg layer 567 mass K0S", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisK0SM});
      registry.add("tracks/ITS/negLayer567MassLambda0", "neg layer 567 mass Lambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisLambdaM});
      registry.add("tracks/ITS/negLayer567MassAntiLambda0", "neg layer 567 mass AntiLambda0", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisAntiLambdaM});

      registry.add("tracks/ITS/negNCl", "neg ncl", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisITSNCl});
      registry.add("tracks/ITS/negChi2NCl", "neg chi2ncl", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisITSChi2NCl});

      // TPC information
      registry.add("tracks/TPC/posNClFindable", "pos ncl findable", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisNClFindable});
      registry.add("tracks/TPC/posNClsFound", "pos ncl found", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisNClFound});
      registry.add("tracks/TPC/posNClsShared", "pos ncl shared", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisNClShared});
      registry.add("tracks/TPC/posNClsCrossedRows", "pos ncl crossed rows", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisNClCrossedRows});
      registry.add("tracks/TPC/posNClsCrossedRowsOverFindableCls", "pos ncl crossed rows over findable cls", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisCrossedRowsOverFindable});
      registry.add("tracks/TPC/posFractionSharedCls", "pos fraction shared cls", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisSharedFraction});
      registry.add("tracks/TPC/posChi2NCl", "pos chi2ncl", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisTPCChi2NCl});

      registry.add("tracks/TPC/negNClFindable", "neg ncl findable", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisNClFindable});
      registry.add("tracks/TPC/negNClsFound", "neg ncl found", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisNClFound});
      registry.add("tracks/TPC/negNClsShared", "neg ncl shared", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisNClShared});
      registry.add("tracks/TPC/negNClsCrossedRows", "neg ncl crossed rows", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisNClCrossedRows});
      registry.add("tracks/TPC/negNClsCrossedRowsOverFindableCls", "neg ncl crossed rows over findable cls", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisCrossedRowsOverFindable});
      registry.add("tracks/TPC/negFractionSharedCls", "neg fraction shared cls", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisSharedFraction});
      registry.add("tracks/TPC/negChi2NCl", "neg chi2ncl", HistType::kTHnSparseD, {axisV0Pt, axisV0Pt, axisV0Pt, axisTPCChi2NCl});
    } // doprocessV0TrackQA
  } // init

  template <typename T, typename U>
  bool isCollisionReconstructed(T const& collision, U const& eventSelectionBits)
  {
    if (!collision.has_mcCollision()) {
      return false;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return false;
    }
    return true;
  }
  template <typename T, typename U, typename V>
  bool v0sAreMatched(T const& v0, U const& particle, V const& /*tracks*/)
  {
    // This is necessary, because the V0Labels table points to aod::McParticles, not to aod::CandidatesV0MCP
    auto negId = v0.template negTrack_as<V>().mcParticleId();
    auto posId = v0.template posTrack_as<V>().mcParticleId();
    auto daughters = particle.daughtersIds();
    return ((negId == daughters[0] && posId == daughters[1]) || (posId == daughters[0] && negId == daughters[1]));
  }

  template <typename T>
  bool hasITSHit(T const& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  template <typename T>
  void fillMcDV0(T const& v0, bool correctCollision, double weight)
  {
    int pdg = v0.mcParticle().pdgCode();
    if (std::abs(pdg) == PDG_t::kK0Short) {
      registry.fill(HIST("inclusive/K0SPtEtaMass"), v0.pt(), v0.eta(), v0.mK0Short(), weight);
      registry.fill(HIST("inclusive/InvMassK0STrue"), v0.pt(), v0.v0radius(), v0.mK0Short(), weight);
      if (!correctCollision)
        registry.fill(HIST("inclusive/K0SPtEtaMassWrongCollision"), v0.pt(), v0.eta(), v0.mK0Short(), weight);
    } else if (pdg == PDG_t::kLambda0) {
      registry.fill(HIST("inclusive/LambdaPtEtaMass"), v0.pt(), v0.eta(), v0.mLambda(), weight);
      registry.fill(HIST("inclusive/InvMassLambdaTrue"), v0.pt(), v0.v0radius(), v0.mLambda(), weight);
      if (!correctCollision)
        registry.fill(HIST("inclusive/LambdaPtEtaMassWrongCollision"), v0.pt(), v0.eta(), v0.mLambda(), weight);
    } else if (pdg == PDG_t::kLambda0Bar) {
      registry.fill(HIST("inclusive/AntiLambdaPtEtaMass"), v0.pt(), v0.eta(), v0.mAntiLambda(), weight);
      registry.fill(HIST("inclusive/InvMassAntiLambdaTrue"), v0.pt(), v0.v0radius(), v0.mAntiLambda(), weight);
      if (!correctCollision)
        registry.fill(HIST("inclusive/AntiLambdaPtEtaMassWrongCollision"), v0.pt(), v0.eta(), v0.mAntiLambda(), weight);
    }
  }

  template <typename T, typename U>
  void fillMcDJets(U const& mcdjet, double weight)
  {
    registry.fill(HIST("jets/JetPtEtaPhi"), mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), weight);
    for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<T>()) {
      registry.fill(HIST("jets/JetsPtEta"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), weight);
    }
  }
  template <typename T, typename U>
  void fillMcDV0InJets(T const& mcdjet, U const& v0, bool correctCollision, double weight)
  {
    int pdg = v0.mcParticle().pdgCode();
    double z = v0.pt() / mcdjet.pt();
    if (std::abs(pdg) == PDG_t::kK0Short) {
      registry.fill(HIST("jets/JetPtEtaK0SPt"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
      registry.fill(HIST("jets/JetPtEtaK0SZ"), mcdjet.pt(), mcdjet.eta(), z, weight);
      if (!correctCollision) {
        registry.fill(HIST("jets/JetPtEtaK0SPtWrongCollision"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
        registry.fill(HIST("jets/JetPtEtaK0SZWrongCollision"), mcdjet.pt(), mcdjet.eta(), z, weight);
      }
    } else if (pdg == PDG_t::kLambda0) {
      registry.fill(HIST("jets/JetPtEtaLambdaPt"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
      registry.fill(HIST("jets/JetPtEtaLambdaZ"), mcdjet.pt(), mcdjet.eta(), z, weight);
      if (!correctCollision) {
        registry.fill(HIST("jets/JetPtEtaLambdaPtWrongCollision"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
        registry.fill(HIST("jets/JetPtEtaLambdaZWrongCollision"), mcdjet.pt(), mcdjet.eta(), z, weight);
      }
    } else if (pdg == PDG_t::kLambda0Bar) {
      registry.fill(HIST("jets/JetPtEtaAntiLambdaPt"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
      registry.fill(HIST("jets/JetPtEtaAntiLambdaZ"), mcdjet.pt(), mcdjet.eta(), z, weight);
      if (!correctCollision) {
        registry.fill(HIST("jets/JetPtEtaAntiLambdaPtWrongCollision"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
        registry.fill(HIST("jets/JetPtEtaAntiLambdaZWrongCollision"), mcdjet.pt(), mcdjet.eta(), z, weight);
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillMcDV0InMatchedJets(T const& mcpjet, U const& mcdjet, V const& v0, bool correctCollision, double weight)
  {
    int pdg = v0.mcParticle().pdgCode();
    double z = v0.pt() / mcdjet.pt();
    if (std::abs(pdg) == PDG_t::kK0Short) {
      registry.fill(HIST("jets/JetsPtEtaK0SPt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
      registry.fill(HIST("jets/JetsPtEtaK0SZ"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), z, weight);
      if (!correctCollision) {
        registry.fill(HIST("jets/JetsPtEtaK0SPtWrongCollision"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
        registry.fill(HIST("jets/JetsPtEtaK0SZWrongCollision"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), z, weight);
      }
    } else if (pdg == PDG_t::kLambda0) {
      registry.fill(HIST("jets/JetsPtEtaLambdaPt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
      registry.fill(HIST("jets/JetsPtEtaLambdaZ"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), z, weight);
      if (!correctCollision) {
        registry.fill(HIST("jets/JetsPtEtaLambdaPtWrongCollision"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
        registry.fill(HIST("jets/JetsPtEtaLambdaZWrongCollision"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), z, weight);
      }
    } else if (pdg == PDG_t::kLambda0Bar) {
      registry.fill(HIST("jets/JetsPtEtaAntiLambdaPt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
      registry.fill(HIST("jets/JetsPtEtaAntiLambdaZ"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), z, weight);
      if (!correctCollision) {
        registry.fill(HIST("jets/JetsPtEtaAntiLambdaPtWrongCollision"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
        registry.fill(HIST("jets/JetsPtEtaAntiLambdaZWrongCollision"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), z, weight);
      }
    }
  }

  template <typename T>
  void fillMcPV0(T const& pv0, double weight)
  {
    // Can calculate this from aod::CandidatesV0MCD (contains decay vertex)
    double rDecay = 1.0;

    if (pv0.pdgCode() == PDG_t::kK0Short) {
      registry.fill(HIST("inclusive/GeneratedK0S"), pv0.pt(), pv0.eta(), rDecay, weight);
    }
    if (pv0.pdgCode() == PDG_t::kLambda0) {
      registry.fill(HIST("inclusive/GeneratedLambda"), pv0.pt(), pv0.eta(), rDecay, weight);
    }
    if (pv0.pdgCode() == PDG_t::kLambda0Bar) {
      registry.fill(HIST("inclusive/GeneratedAntiLambda"), pv0.pt(), pv0.eta(), rDecay, weight);
    }
  }

  template <typename T>
  void fillMcPJets(T const& jet, double weight)
  {
    registry.fill(HIST("jets/GeneratedJets"), jet.pt(), jet.eta(), jet.phi(), weight);
  }

  template <typename T, typename U>
  void fillMcPV0InJets(T const& jet, U const& pv0, double weight)
  {
    double z = pv0.pt() / jet.pt();
    if (pv0.pdgCode() == PDG_t::kK0Short) {
      registry.fill(HIST("jets/GeneratedJetK0S"), jet.pt(), jet.eta(), pv0.pt(), weight);
      registry.fill(HIST("jets/GeneratedJetK0SFrag"), jet.pt(), jet.eta(), z, weight);
    }
    if (pv0.pdgCode() == PDG_t::kLambda0) {
      registry.fill(HIST("jets/GeneratedJetLambda"), jet.pt(), jet.eta(), pv0.pt(), weight);
      registry.fill(HIST("jets/GeneratedJetLambdaFrag"), jet.pt(), jet.eta(), z, weight);
    }
    if (pv0.pdgCode() == PDG_t::kLambda0Bar) {
      registry.fill(HIST("jets/GeneratedJetAntiLambda"), jet.pt(), jet.eta(), pv0.pt(), weight);
      registry.fill(HIST("jets/GeneratedJetAntiLambdaFrag"), jet.pt(), jet.eta(), z, weight);
    }
  }

  template <typename T, typename U, typename V>
  void fillTrackQa(V const& v0)
  {
    auto posTrack = v0.template posTrack_as<T>().template track_as<U>();
    auto negTrack = v0.template negTrack_as<T>().template track_as<U>();

    double mK = v0.mK0Short();
    double mL = v0.mLambda();
    double mAL = v0.mAntiLambda();

    double vPt = v0.pt();
    double pPt = posTrack.pt();
    double nPt = negTrack.pt();
    double dPt = posTrack.pt() - negTrack.pt();

    registry.fill(HIST("tracks/Pos"), vPt, pPt, posTrack.eta(), posTrack.phi());
    registry.fill(HIST("tracks/Neg"), vPt, nPt, negTrack.eta(), negTrack.phi());
    registry.fill(HIST("tracks/Pt"), vPt, pPt, nPt, dPt, dPt / vPt);
    registry.fill(HIST("tracks/PtMass"), vPt, pPt, nPt, mK, mL, mAL);
    registry.fill(HIST("tracks/PtDiffMass"), vPt, dPt, mK, mL, mAL);

    registry.fill(HIST("tracks/DCAxy"), vPt, pPt, nPt, posTrack.dcaXY(), negTrack.dcaXY(), dPt);
    registry.fill(HIST("tracks/DCAxyMassK0S"), vPt, pPt, nPt, posTrack.dcaXY(), negTrack.dcaXY(), mK);
    registry.fill(HIST("tracks/DCAxyMassLambda0"), vPt, pPt, nPt, posTrack.dcaXY(), negTrack.dcaXY(), mL);
    registry.fill(HIST("tracks/DCAxyMassAntiLambda0"), vPt, pPt, nPt, posTrack.dcaXY(), negTrack.dcaXY(), mAL);

    registry.fill(HIST("tracks/DCAz"), vPt, pPt, nPt, posTrack.dcaZ(), negTrack.dcaZ(), dPt);
    registry.fill(HIST("tracks/DCAzMassK0S"), vPt, pPt, nPt, posTrack.dcaZ(), negTrack.dcaZ(), mK);
    registry.fill(HIST("tracks/DCAzMassLambda0"), vPt, pPt, nPt, posTrack.dcaZ(), negTrack.dcaZ(), mL);
    registry.fill(HIST("tracks/DCAzMassAntiLambda0"), vPt, pPt, nPt, posTrack.dcaZ(), negTrack.dcaZ(), mAL);

    registry.fill(HIST("tracks/V0Radius"), vPt, pPt, nPt, v0.v0radius(), dPt);
    registry.fill(HIST("tracks/V0RadiusMassK0S"), vPt, pPt, nPt, v0.v0radius(), mK);
    registry.fill(HIST("tracks/V0RadiusMassLambda0"), vPt, pPt, nPt, v0.v0radius(), mL);
    registry.fill(HIST("tracks/V0RadiusMassAntiLambda0"), vPt, pPt, nPt, v0.v0radius(), mAL);

    registry.fill(HIST("tracks/V0CosPa"), vPt, pPt, nPt, v0.v0cosPA(), dPt);
    registry.fill(HIST("tracks/V0CosPaMassK0S"), vPt, pPt, nPt, v0.v0cosPA(), mK);
    registry.fill(HIST("tracks/V0CosPaMassLambda0"), vPt, pPt, nPt, v0.v0cosPA(), mL);
    registry.fill(HIST("tracks/V0CosPaMassAntiLambda0"), vPt, pPt, nPt, v0.v0cosPA(), mAL);

    // Has TRD or not
    if (posTrack.hasTRD()) {
      registry.fill(HIST("tracks/posTRDPt"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/posTRDPtMass"), vPt, pPt, nPt, mK, mL, mAL);
    } else {
      registry.fill(HIST("tracks/posNoTRDPt"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/posNoTRDPtMass"), vPt, pPt, nPt, mK, mL, mAL);
    }
    if (negTrack.hasTRD()) {
      registry.fill(HIST("tracks/negTRDPt"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/negTRDPtMass"), vPt, pPt, nPt, mK, mL, mAL);
    } else {
      registry.fill(HIST("tracks/negNoTRDPt"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/negNoTRDPtMass"), vPt, pPt, nPt, mK, mL, mAL);
    }

    // ITS information
    if (hasITSHit(posTrack, 1)) {
      registry.fill(HIST("tracks/ITS/posLayer1"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer1MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer1MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer1MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 2)) {
      registry.fill(HIST("tracks/ITS/posLayer2"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer2MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer2MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer2MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 3)) {
      registry.fill(HIST("tracks/ITS/posLayer3"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer3MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer3MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer3MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 4)) {
      registry.fill(HIST("tracks/ITS/posLayer4"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer4MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer4MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer4MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 5)) {
      registry.fill(HIST("tracks/ITS/posLayer5"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer5MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer5MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer5MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 6)) {
      registry.fill(HIST("tracks/ITS/posLayer6"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer6MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer6MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer6MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 7)) {
      registry.fill(HIST("tracks/ITS/posLayer7"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer7MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer7MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer7MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 5) && hasITSHit(posTrack, 6)) {
      registry.fill(HIST("tracks/ITS/posLayer56"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer56MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer56MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer56MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 6) && hasITSHit(posTrack, 7)) {
      registry.fill(HIST("tracks/ITS/posLayer67"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer67MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer67MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer67MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 5) && hasITSHit(posTrack, 7)) {
      registry.fill(HIST("tracks/ITS/posLayer57"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer57MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer57MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer57MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(posTrack, 5) && hasITSHit(posTrack, 6) && hasITSHit(posTrack, 7)) {
      registry.fill(HIST("tracks/ITS/posLayer567"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/posLayer567MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/posLayer567MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/posLayer567MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    registry.fill(HIST("tracks/ITS/posNCl"), vPt, pPt, nPt, posTrack.itsNCls());
    registry.fill(HIST("tracks/ITS/posChi2NCl"), vPt, pPt, nPt, posTrack.itsChi2NCl());

    if (hasITSHit(negTrack, 1)) {
      registry.fill(HIST("tracks/ITS/negLayer1"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer1MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer1MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer1MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 2)) {
      registry.fill(HIST("tracks/ITS/negLayer2"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer2MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer2MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer2MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 3)) {
      registry.fill(HIST("tracks/ITS/negLayer3"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer3MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer3MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer3MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 4)) {
      registry.fill(HIST("tracks/ITS/negLayer4"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer4MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer4MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer4MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 5)) {
      registry.fill(HIST("tracks/ITS/negLayer5"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer5MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer5MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer5MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 6)) {
      registry.fill(HIST("tracks/ITS/negLayer6"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer6MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer6MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer6MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 7)) {
      registry.fill(HIST("tracks/ITS/negLayer7"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer7MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer7MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer7MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 5) && hasITSHit(negTrack, 6)) {
      registry.fill(HIST("tracks/ITS/negLayer56"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer56MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer56MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer56MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 6) && hasITSHit(negTrack, 7)) {
      registry.fill(HIST("tracks/ITS/negLayer67"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer67MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer67MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer67MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 5) && hasITSHit(negTrack, 7)) {
      registry.fill(HIST("tracks/ITS/negLayer57"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer57MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer57MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer57MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    if (hasITSHit(negTrack, 5) && hasITSHit(negTrack, 6) && hasITSHit(negTrack, 7)) {
      registry.fill(HIST("tracks/ITS/negLayer567"), vPt, pPt, nPt, dPt);
      registry.fill(HIST("tracks/ITS/negLayer567MassK0S"), vPt, pPt, nPt, mK);
      registry.fill(HIST("tracks/ITS/negLayer567MassLambda0"), vPt, pPt, nPt, mL);
      registry.fill(HIST("tracks/ITS/negLayer567MassAntiLambda0"), vPt, pPt, nPt, mAL);
    }
    registry.fill(HIST("tracks/ITS/negNCl"), vPt, pPt, nPt, negTrack.itsNCls());
    registry.fill(HIST("tracks/ITS/negChi2NCl"), vPt, pPt, nPt, negTrack.itsChi2NCl());

    // TPC information
    registry.fill(HIST("tracks/TPC/posNClFindable"), vPt, pPt, nPt, posTrack.tpcNClsFindable());
    registry.fill(HIST("tracks/TPC/posNClsFound"), vPt, pPt, nPt, posTrack.tpcNClsFound());
    registry.fill(HIST("tracks/TPC/posChi2NCl"), vPt, pPt, nPt, posTrack.tpcChi2NCl());
    registry.fill(HIST("tracks/TPC/posNClsShared"), vPt, pPt, nPt, posTrack.tpcNClsShared());
    registry.fill(HIST("tracks/TPC/posFractionSharedCls"), vPt, pPt, nPt, posTrack.tpcFractionSharedCls());
    registry.fill(HIST("tracks/TPC/posNClsCrossedRows"), vPt, pPt, nPt, posTrack.tpcNClsCrossedRows());
    registry.fill(HIST("tracks/TPC/posNClsCrossedRowsOverFindableCls"), vPt, pPt, nPt, posTrack.tpcCrossedRowsOverFindableCls());

    registry.fill(HIST("tracks/TPC/negNClFindable"), vPt, pPt, nPt, negTrack.tpcNClsFindable());
    registry.fill(HIST("tracks/TPC/negNClsFound"), vPt, pPt, nPt, negTrack.tpcNClsFound());
    registry.fill(HIST("tracks/TPC/negChi2NCl"), vPt, pPt, nPt, negTrack.tpcChi2NCl());
    registry.fill(HIST("tracks/TPC/negNClsShared"), vPt, pPt, nPt, negTrack.tpcNClsShared());
    registry.fill(HIST("tracks/TPC/negFractionSharedCls"), vPt, pPt, nPt, negTrack.tpcFractionSharedCls());
    registry.fill(HIST("tracks/TPC/negNClsCrossedRows"), vPt, pPt, nPt, negTrack.tpcNClsCrossedRows());
    registry.fill(HIST("tracks/TPC/negNClsCrossedRowsOverFindableCls"), vPt, pPt, nPt, negTrack.tpcCrossedRowsOverFindableCls());
  }

  void processDummy(aod::CandidatesV0MCD const&) {}
  PROCESS_SWITCH(V0QA, processDummy, "Dummy process function turned on by default", true);

  void processFlags(soa::Join<aod::V0Datas, aod::V0SignalFlags>::iterator const& v0)
  {
    int isK0S = static_cast<int>(v0.isK0SCandidate());
    int isLambda = static_cast<int>((v0.isLambdaCandidate()));
    int isAntiLambda = static_cast<int>(v0.isAntiLambdaCandidate());
    int isRejected = static_cast<int>(v0.isRejectedCandidate());

    registry.fill(HIST("inclusive/V0Flags"), 0, 0, isRejected);
    registry.fill(HIST("inclusive/V0Flags"), 1, 1, isK0S);
    registry.fill(HIST("inclusive/V0Flags"), 2, 2, isLambda);
    registry.fill(HIST("inclusive/V0Flags"), 3, 3, isAntiLambda);

    registry.fill(HIST("inclusive/V0Flags"), 0, 1, isRejected * isK0S);
    registry.fill(HIST("inclusive/V0Flags"), 1, 0, isRejected * isK0S);
    registry.fill(HIST("inclusive/V0Flags"), 0, 2, isRejected * isLambda);
    registry.fill(HIST("inclusive/V0Flags"), 2, 0, isRejected * isLambda);
    registry.fill(HIST("inclusive/V0Flags"), 0, 3, isRejected * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 3, 0, isRejected * isAntiLambda);

    registry.fill(HIST("inclusive/V0Flags"), 1, 2, isK0S * isLambda);
    registry.fill(HIST("inclusive/V0Flags"), 2, 1, isK0S * isLambda);
    registry.fill(HIST("inclusive/V0Flags"), 1, 3, isK0S * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 3, 1, isK0S * isAntiLambda);

    registry.fill(HIST("inclusive/V0Flags"), 2, 3, isLambda * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 3, 2, isLambda * isAntiLambda);

    // V0 satisfies 3+ classes
    registry.fill(HIST("inclusive/V0Flags"), 0, 4, isRejected * isK0S * isLambda);
    registry.fill(HIST("inclusive/V0Flags"), 4, 0, isRejected * isK0S * isLambda);
    registry.fill(HIST("inclusive/V0Flags"), 1, 4, isRejected * isK0S * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 4, 1, isRejected * isK0S * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 2, 4, isRejected * isLambda * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 4, 2, isRejected * isLambda * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 3, 4, isRejected * isK0S * isLambda * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 4, 3, isRejected * isK0S * isLambda * isAntiLambda);
    registry.fill(HIST("inclusive/V0Flags"), 4, 4, isK0S * isLambda * isAntiLambda);
  }
  PROCESS_SWITCH(V0QA, processFlags, "V0 flags", false);

  void processMcD(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, CandidatesV0MCDWithFlags const& v0s, aod::McParticles const&, MatchedMCDV0JetsWithConstituents const& mcdjets, MatchedMCPV0JetsWithConstituents const&, aod::CandidatesV0MCP const&, aod::JetTracksMCD const& jTracks, JetMcCollisionsWithPIs const&, aod::McCollisions const&)
  {
    registry.fill(HIST("inclusive/hEvents"), 0.5);
    if (!isCollisionReconstructed(jcoll, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("inclusive/hEvents"), 1.5);
    auto mcColl = jcoll.template mcCollision_as<JetMcCollisionsWithPIs>();
    double weight = mcColl.weight();
    registry.fill(HIST("inclusive/hEvents"), 2.5, weight);

    for (const auto& v0 : v0s) {
      if (!v0.has_mcParticle() || v0.isRejectedCandidate())
        continue;

      bool correctCollision = (mcColl.mcCollisionId() == v0.mcParticle().mcCollisionId());
      fillMcDV0(v0, correctCollision, weight);
    } // v0 loop

    for (const auto& mcdjet : mcdjets) {
      fillMcDJets<MatchedMCPV0JetsWithConstituents>(mcdjet, weight);
      for (const auto& v0 : mcdjet.template candidates_as<CandidatesV0MCDWithFlags>()) {
        if (!v0.has_mcParticle() || v0.isRejectedCandidate())
          continue;

        bool correctCollision = (mcColl.mcCollisionId() == v0.mcParticle().mcCollisionId());
        fillMcDV0InJets(mcdjet, v0, correctCollision, weight);

        for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<MatchedMCPV0JetsWithConstituents>()) {
          for (const auto& pv0 : mcpjet.template candidates_as<aod::CandidatesV0MCP>()) {
            if (!v0sAreMatched(v0, pv0, jTracks))
              continue;

            fillMcDV0InMatchedJets(mcpjet, mcdjet, v0, correctCollision, weight);
          } // v0 particle loop
        } // mcpjet loop
      } // v0 loop
    } // mcd jet loop
  }
  PROCESS_SWITCH(V0QA, processMcD, "Reconstructed true V0s", false);

  void processMcP(aod::JetMcCollision const& mccoll, soa::SmallGroups<aod::JetCollisionsMCD> const& collisions, MCPV0JetsWithConstituents const& jets, aod::CandidatesV0MCP const& pv0s)
  {
    registry.fill(HIST("inclusive/hMcEvents"), 0.5);
    bool isReconstructed = false;

    for (const auto& collision : collisions) {
      if (!isCollisionReconstructed(collision, eventSelectionBits))
        continue;

      if (collision.mcCollision().globalIndex() != mccoll.globalIndex())
        continue;

      isReconstructed = true;
      break;
    }
    if (!isReconstructed)
      return;

    registry.fill(HIST("inclusive/hMcEvents"), 1.5);
    double weight = mccoll.weight();
    registry.fill(HIST("inclusive/hMcEvents"), 2.5, weight);

    for (const auto& pv0 : pv0s) {
      if (!pv0.has_daughters() || !pv0.isPhysicalPrimary())
        continue;
      if (std::abs(pv0.y()) > yPartMax)
        continue;

      fillMcPV0(pv0, weight);
    }

    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, -99., -99., -1. * yPartMax, yPartMax))
        continue;

      fillMcPJets(jet, weight);

      for (const auto& pv0 : jet.template candidates_as<aod::CandidatesV0MCP>()) {
        if (!pv0.has_daughters() || !pv0.isPhysicalPrimary())
          continue;

        fillMcPV0InJets(jet, pv0, weight);
      }
    }
  }
  PROCESS_SWITCH(V0QA, processMcP, "Particle level V0s", false);

  void processCollisionAssociation(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, CandidatesV0MCDWithFlags const& v0s, JetMcCollisionsWithPIs const&, aod::McCollisions const&, aod::McParticles const&)
  {
    // Based on PWGLF/Tasks/Strangeness/derivedlambdakzeroanalysis.cxx
    if (!jcoll.has_mcCollision())
      return;

    auto mcColl = jcoll.template mcCollision_as<JetMcCollisionsWithPIs>();
    double weight = mcColl.weight();

    for (const auto& v0 : v0s) {
      if (!v0.has_mcParticle())
        continue;

      auto pv0 = v0.mcParticle();
      bool correctCollision = (mcColl.mcCollisionId() == v0.mcParticle().mcCollisionId());
      int pdg = v0.mcParticle().pdgCode();

      // Check V0 decay kinematics
      if (v0.isRejectedCandidate())
        continue;

      registry.fill(HIST("collisions/V0PtEta"), pv0.pt(), pv0.eta(), weight);
      if (!correctCollision) {
        registry.fill(HIST("collisions/V0PtEtaWrongColl"), pv0.pt(), pv0.eta(), weight);
      }
      if (std::abs(pdg) == PDG_t::kK0Short) {
        registry.fill(HIST("collisions/K0SPtEtaMass"), pv0.pt(), pv0.eta(), v0.mK0Short(), weight);
        if (!correctCollision) {
          registry.fill(HIST("collisions/K0SPtEtaMassWrongColl"), pv0.pt(), pv0.eta(), v0.mK0Short(), weight);
        }
      }
      if (pdg == PDG_t::kLambda0) {
        registry.fill(HIST("collisions/LambdaPtEtaMass"), pv0.pt(), pv0.eta(), v0.mLambda(), weight);
        if (!correctCollision) {
          registry.fill(HIST("collisions/LambdaPtEtaMassWrongColl"), pv0.pt(), pv0.eta(), v0.mLambda(), weight);
        }
      }
      if (pdg == PDG_t::kLambda0Bar) {
        registry.fill(HIST("collisions/AntiLambdaPtEtaMass"), pv0.pt(), pv0.eta(), v0.mAntiLambda(), weight);
        if (!correctCollision) {
          registry.fill(HIST("collisions/AntiLambdaPtEtaMassWrongColl"), pv0.pt(), pv0.eta(), v0.mAntiLambda(), weight);
        }
      }
      // Feed-down from Xi
      if (!v0.has_mcMotherParticle())
        continue;

      auto mother = v0.mcMotherParticle();
      pdg = mother.pdgCode();
      correctCollision = (mcColl.mcCollisionId() == mother.mcCollisionId());

      if (pdg == PDG_t::kXiMinus) {
        registry.fill(HIST("collisions/XiMinusPtYLambdaPt"), mother.pt(), mother.y(), pv0.pt(), weight);
        if (!correctCollision) {
          registry.fill(HIST("collisions/XiMinusPtYLambdaPtWrongColl"), mother.pt(), mother.y(), pv0.pt(), weight);
        }
      }
      if (pdg == PDG_t::kXiPlusBar) {
        registry.fill(HIST("collisions/XiPlusPtYAntiLambdaPt"), mother.pt(), mother.y(), pv0.pt(), weight);
        if (!correctCollision) {
          registry.fill(HIST("collisions/XiPlusPtYAntiLambdaPtWrongColl"), mother.pt(), mother.y(), pv0.pt(), weight);
        }
      }
    }
  }
  PROCESS_SWITCH(V0QA, processCollisionAssociation, "V0 collision association", false);

  void processCollisionAssociationJets(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, MCDV0JetsWithConstituents const& mcdjets, CandidatesV0MCDWithFlags const&, JetMcCollisionsWithPIs const&, aod::McCollisions const&, aod::McParticles const&)
  {
    if (!isCollisionReconstructed(jcoll, eventSelectionBits))
      return;

    auto mcColl = jcoll.template mcCollision_as<JetMcCollisionsWithPIs>();
    double weight = mcColl.weight();

    for (const auto& mcdjet : mcdjets) {
      // Eta cut?
      for (const auto& v0 : mcdjet.template candidates_as<CandidatesV0MCDWithFlags>()) {
        if (!v0.has_mcParticle())
          continue;

        auto pv0 = v0.mcParticle();
        bool correctCollision = (mcColl.mcCollisionId() == pv0.mcCollisionId());
        int pdg = pv0.pdgCode();
        double z = v0.pt() / mcdjet.pt();

        // Check V0 decay kinematics
        if (v0.isRejectedCandidate())
          continue;

        registry.fill(HIST("collisions/JetPtEtaV0Pt"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
        if (!correctCollision) {
          registry.fill(HIST("collisions/JetPtEtaV0PtWrongColl"), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
        }
        if (std::abs(pdg) == PDG_t::kK0Short) {
          registry.fill(HIST("collisions/JetPtEtaK0SPtMass"), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mK0Short(), weight);
          registry.fill(HIST("collisions/JetPtEtaK0SFragMass"), mcdjet.pt(), mcdjet.eta(), z, v0.mK0Short(), weight);
          if (!correctCollision) {
            registry.fill(HIST("collisions/JetPtEtaK0SPtMassWrongColl"), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mK0Short(), weight);
            registry.fill(HIST("collisions/JetPtEtaK0SFragMassWrongColl"), mcdjet.pt(), mcdjet.eta(), z, v0.mK0Short(), weight);
          }
        }
        if (pdg == PDG_t::kLambda0) {
          registry.fill(HIST("collisions/JetPtEtaLambdaPtMass"), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mLambda(), weight);
          registry.fill(HIST("collisions/JetPtEtaLambdaFragMass"), mcdjet.pt(), mcdjet.eta(), z, v0.mLambda(), weight);
          if (!correctCollision) {
            registry.fill(HIST("collisions/JetPtEtaLambdaPtMassWrongColl"), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mLambda(), weight);
            registry.fill(HIST("collisions/JetPtEtaLambdaFragMassWrongColl"), mcdjet.pt(), mcdjet.eta(), z, v0.mLambda(), weight);
          }
        }
        if (pdg == PDG_t::kLambda0Bar) {
          registry.fill(HIST("collisions/JetPtEtaAntiLambdaPtMass"), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mAntiLambda(), weight);
          registry.fill(HIST("collisions/JetPtEtaAntiLambdaFragMass"), mcdjet.pt(), mcdjet.eta(), z, v0.mAntiLambda(), weight);
          if (!correctCollision) {
            registry.fill(HIST("collisions/JetPtEtaAntiLambdaPtMassWrongColl"), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mAntiLambda(), weight);
            registry.fill(HIST("collisions/JetPtEtaAntiLambdaFragMassWrongColl"), mcdjet.pt(), mcdjet.eta(), z, v0.mAntiLambda(), weight);
          }
        }

        if (!v0.has_mcMotherParticle())
          continue;

        auto mother = v0.mcMotherParticle();
        pdg = mother.pdgCode();
        correctCollision = (mcColl.mcCollisionId() == mother.mcCollisionId());
        if (pdg == PDG_t::kXiMinus) {
          registry.fill(HIST("collisions/JetPtEtaXiMinusPtLambdaPt"), mcdjet.pt(), mcdjet.eta(), mother.pt(), v0.pt(), weight);
          if (!correctCollision) {
            registry.fill(HIST("collisions/JetPtEtaXiMinusPtLambdaPtWrongColl"), mcdjet.pt(), mcdjet.eta(), mother.pt(), v0.pt(), weight);
          }
        }
        if (pdg == PDG_t::kXiPlusBar) {
          registry.fill(HIST("collisions/JetPtEtaXiPlusPtAntiLambdaPt"), mcdjet.pt(), mcdjet.eta(), mother.pt(), v0.pt(), weight);
          if (!correctCollision) {
            registry.fill(HIST("collisions/JetPtEtaXiPlusPtAntiLambdaPtWrongColl"), mcdjet.pt(), mcdjet.eta(), mother.pt(), v0.pt(), weight);
          }
        }
      } // for v0s
    } // for mcdjets
  }
  PROCESS_SWITCH(V0QA, processCollisionAssociationJets, "V0 in jets collision association", false);

  void processCollisionAssociationMatchedJets(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, MatchedMCDV0JetsWithConstituents const& mcdjets, MatchedMCPV0JetsWithConstituents const&, CandidatesV0MCDWithFlags const&, aod::CandidatesV0MCP const&, JetMcCollisionsWithPIs const&, aod::McCollisions const&, aod::McParticles const&, aod::JetTracksMCD const& jTracks)
  {
    if (!jcoll.has_mcCollision())
      return;

    auto mcColl = jcoll.template mcCollision_as<JetMcCollisionsWithPIs>();
    double weight = mcColl.weight();

    for (const auto& mcdjet : mcdjets) {
      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<MatchedMCPV0JetsWithConstituents>()) {
        for (const auto& v0 : mcdjet.template candidates_as<CandidatesV0MCDWithFlags>()) {
          if (!v0.has_mcParticle())
            continue;

          for (const auto& pv0 : mcpjet.template candidates_as<aod::CandidatesV0MCP>()) {
            if (!v0sAreMatched(v0, pv0, jTracks))
              continue;

            int pdg = pv0.pdgCode();
            bool correctCollision = (mcColl.mcCollisionId() == pv0.mcCollisionId());
            double z = v0.pt() / mcdjet.pt();

            // Check V0 decay kinematics
            if (v0.isRejectedCandidate())
              continue;

            registry.fill(HIST("collisions/JetsPtEtaV0Pt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
            if (!correctCollision) {
              registry.fill(HIST("collisions/JetsPtEtaV0PtWrongColl"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
            }
            if (std::abs(pdg) == PDG_t::kK0Short) {
              registry.fill(HIST("collisions/JetsPtEtaK0SPtMass"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mK0Short(), weight);
              registry.fill(HIST("collisions/JetsPtEtaK0SFragMass"), mcpjet.pt(), mcdjet.eta(), z, v0.mK0Short(), weight);
              if (!correctCollision) {
                registry.fill(HIST("collisions/JetsPtEtaK0SPtMassWrongColl"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mK0Short(), weight);
                registry.fill(HIST("collisions/JetsPtEtaK0SFragMassWrongColl"), mcpjet.pt(), mcdjet.eta(), z, v0.mK0Short(), weight);
              }
            }
            if (pdg == PDG_t::kLambda0) {
              registry.fill(HIST("collisions/JetsPtEtaLambdaPtMass"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mLambda(), weight);
              registry.fill(HIST("collisions/JetsPtEtaLambdaFragMass"), mcpjet.pt(), mcdjet.eta(), z, v0.mLambda(), weight);
              if (!correctCollision) {
                registry.fill(HIST("collisions/JetsPtEtaLambdaPtMassWrongColl"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mLambda(), weight);
                registry.fill(HIST("collisions/JetsPtEtaLambdaFragMassWrongColl"), mcpjet.pt(), mcdjet.eta(), z, v0.mLambda(), weight);
              }
            }
            if (pdg == PDG_t::kLambda0Bar) {
              registry.fill(HIST("collisions/JetsPtEtaAntiLambdaPtMass"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mAntiLambda(), weight);
              registry.fill(HIST("collisions/JetsPtEtaAntiLambdaFragMass"), mcpjet.pt(), mcdjet.eta(), z, v0.mAntiLambda(), weight);
              if (!correctCollision) {
                registry.fill(HIST("collisions/JetsPtEtaAntiLambdaPtMassWrongColl"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), v0.mAntiLambda(), weight);
                registry.fill(HIST("collisions/JetsPtEtaAntiLambdaFragMassWrongColl"), mcpjet.pt(), mcdjet.eta(), z, v0.mAntiLambda(), weight);
              }
            }

            if (!v0.has_mcMotherParticle())
              continue;

            auto mother = v0.mcMotherParticle();
            pdg = mother.pdgCode();
            correctCollision = (mcColl.mcCollisionId() == mother.mcCollisionId());
            if (pdg == PDG_t::kXiMinus) {
              registry.fill(HIST("collisions/JetsPtEtaXiMinusPtLambdaPt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), mother.pt(), v0.pt(), weight);
              if (!correctCollision) {
                registry.fill(HIST("collisions/JetsPtEtaXiMinusPtLambdaPtWrongColl"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), mother.pt(), v0.pt(), weight);
              }
            }
            if (pdg == PDG_t::kXiPlusBar) {
              registry.fill(HIST("collisions/JetsPtEtaXiPlusPtAntiLambdaPt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), mother.pt(), v0.pt(), weight);
              if (!correctCollision) {
                registry.fill(HIST("collisions/JetsPtEtaXiPlusPtAntiLambdaPtWrongColl"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), mother.pt(), v0.pt(), weight);
              }
            }
          } // for pv0
        } // for v0
      } // for mcpjet
    } // for mcdjet
  }
  PROCESS_SWITCH(V0QA, processCollisionAssociationMatchedJets, "V0 in matched jets collision association", false);

  void processFeeddown(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, CandidatesV0MCDWithFlags const& v0s, aod::CandidatesV0MCP const&, JetMcCollisionsWithPIs const&, aod::McCollisions const&, aod::McParticles const&)
  {
    // Based on PWGLF/Tasks/Strangeness/derivedlambdakzeroanalysis.cxx
    if (!jcoll.has_mcCollision())
      return;

    auto mcColl = jcoll.template mcCollision_as<JetMcCollisionsWithPIs>();
    double weight = mcColl.weight();

    for (const auto& v0 : v0s) {
      if (!v0.has_mcParticle())
        continue;

      int pdg = v0.mcParticle().pdgCode();

      // Check V0 decay kinematics
      if (v0.isRejectedCandidate())
        continue;
      // Feed-down from Xi
      if (!v0.has_mcMotherParticle())
        continue;

      auto pv0 = v0.mcParticle();
      auto mother = v0.mcMotherParticle();
      pdg = mother.pdgCode();

      if (pdg == PDG_t::kXiMinus) {
        registry.fill(HIST("feeddown/XiMinusPtYLambdaPt"), mother.pt(), mother.y(), pv0.pt(), weight);
      }
      if (pdg == PDG_t::kXiPlusBar) {
        registry.fill(HIST("feeddown/XiPlusPtYAntiLambdaPt"), mother.pt(), mother.y(), pv0.pt(), weight);
      }
    }
  }
  PROCESS_SWITCH(V0QA, processFeeddown, "Inclusive feeddown", false);

  void processFeeddownJets(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, MCDV0JetsWithConstituents const& mcdjets, CandidatesV0MCDWithFlags const&, aod::CandidatesV0MCP const&, JetMcCollisionsWithPIs const&, aod::McCollisions const&, aod::McParticles const&)
  {
    // Based on PWGLF/Tasks/Strangeness/derivedlambdakzeroanalysis.cxx
    if (!jcoll.has_mcCollision())
      return;

    auto mcColl = jcoll.template mcCollision_as<JetMcCollisionsWithPIs>();
    double weight = mcColl.weight();

    for (const auto& mcdjet : mcdjets) {
      for (const auto& v0 : mcdjet.template candidates_as<CandidatesV0MCDWithFlags>()) {
        if (!v0.has_mcParticle())
          continue;

        int pdg = v0.mcParticle().pdgCode();

        // Check V0 decay kinematics
        if (v0.isRejectedCandidate())
          continue;
        // Feed-down from Xi
        if (!v0.has_mcMotherParticle())
          continue;

        auto pv0 = v0.mcParticle();
        auto mother = v0.mcMotherParticle();
        pdg = mother.pdgCode();

        if (pdg == PDG_t::kXiMinus) {
          registry.fill(HIST("feeddown/JetPtXiMinusPtLambdaPt"), mcdjet.pt(), mother.pt(), pv0.pt(), weight);
        }
        if (pdg == PDG_t::kXiPlusBar) {
          registry.fill(HIST("feeddown/JetPtXiPlusPtAntiLambdaPt"), mcdjet.pt(), mother.pt(), pv0.pt(), weight);
        }
      }
    }
  }
  PROCESS_SWITCH(V0QA, processFeeddownJets, "Jets feeddown", false);

  void processFeeddownMatchedJets(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, MatchedMCDV0JetsWithConstituents const& mcdjets, aod::JetTracksMCD const& jTracks, MatchedMCPV0JetsWithConstituents const&, CandidatesV0MCDWithFlags const&, aod::CandidatesV0MCP const&, JetMcCollisionsWithPIs const&, aod::McCollisions const&, aod::McParticles const&)
  {
    // Based on PWGLF/Tasks/Strangeness/derivedlambdakzeroanalysis.cxx
    if (!jcoll.has_mcCollision())
      return;

    auto mcColl = jcoll.template mcCollision_as<JetMcCollisionsWithPIs>();
    double weight = mcColl.weight();

    for (const auto& mcdjet : mcdjets) {
      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<MatchedMCPV0JetsWithConstituents>()) {
        for (const auto& v0 : mcdjet.template candidates_as<CandidatesV0MCDWithFlags>()) {
          if (!v0.has_mcParticle())
            continue;
          if (!v0.has_mcMotherParticle())
            continue;

          for (const auto& pv0 : mcpjet.template candidates_as<aod::CandidatesV0MCP>()) {
            if (!v0sAreMatched(v0, pv0, jTracks))
              continue;

            int pdg = v0.mcParticle().pdgCode();

            // Check V0 decay kinematics
            if (v0.isRejectedCandidate())
              continue;

            auto mother = v0.mcMotherParticle();
            pdg = mother.pdgCode();
            if (pdg == PDG_t::kXiMinus) {
              registry.fill(HIST("feeddown/JetsPtXiMinusPtLambdaPt"), mcpjet.pt(), mcdjet.pt(), mother.pt(), pv0.pt(), weight);
            }
            if (pdg == PDG_t::kXiPlusBar) {
              registry.fill(HIST("feeddown/JetsPtXiPlusPtAntiLambdaPt"), mcpjet.pt(), mcdjet.pt(), mother.pt(), pv0.pt(), weight);
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(V0QA, processFeeddownMatchedJets, "Jets feeddown", false);

  // Test the difference between excluding V0s from jet finding and subtracting V0s from jets afterwards
  void processTestWeightedJetFinder(soa::Filtered<aod::JetCollisions>::iterator const& jcoll, soa::Join<aod::V0ChargedJets, aod::V0ChargedJetConstituents> const& jets, aod::CandidatesV0Data const&)
  {
    registry.fill(HIST("tests/weighted/hEvents"), 0.5);
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelectionBits))
      return;

    registry.fill(HIST("tests/weighted/hEvents"), 1.5);

    for (const auto& jet : jets) {
      registry.fill(HIST("tests/weighted/JetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());

      for (const auto& v0 : jet.template candidates_as<aod::CandidatesV0Data>()) {
        if (v0.isRejectedCandidate())
          continue;

        double z = v0.pt() / jet.pt();

        registry.fill(HIST("tests/weighted/JetPtEtaV0Pt"), jet.pt(), jet.eta(), v0.pt());
        registry.fill(HIST("tests/weighted/JetPtEtaV0Z"), jet.pt(), jet.eta(), z);

        if (v0.isK0SCandidate()) {
          registry.fill(HIST("tests/weighted/JetPtEtaK0SPt"), jet.pt(), jet.eta(), v0.pt());
          registry.fill(HIST("tests/weighted/JetPtEtaK0SZ"), jet.pt(), jet.eta(), z);
        }
        if (v0.isLambdaCandidate()) {
          registry.fill(HIST("tests/weighted/JetPtEtaLambdaPt"), jet.pt(), jet.eta(), v0.pt());
          registry.fill(HIST("tests/weighted/JetPtEtaLambdaZ"), jet.pt(), jet.eta(), z);
        }
        if (v0.isAntiLambdaCandidate()) {
          registry.fill(HIST("tests/weighted/JetPtEtaAntiLambdaPt"), jet.pt(), jet.eta(), v0.pt());
          registry.fill(HIST("tests/weighted/JetPtEtaAntiLambdaZ"), jet.pt(), jet.eta(), z);
        }
      }
    }
  }
  PROCESS_SWITCH(V0QA, processTestWeightedJetFinder, "Test weighted jet finder", false);

  void processTestSubtractedJetFinder(soa::Filtered<aod::JetCollisions>::iterator const& jcoll, soa::Join<aod::V0ChargedJets, aod::V0ChargedJetConstituents> const& jets, aod::CandidatesV0Data const&)
  {
    registry.fill(HIST("tests/hEvents"), 0.5);
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelectionBits))
      return;

    registry.fill(HIST("tests/hEvents"), 1.5);

    for (const auto& jet : jets) {
      registry.fill(HIST("tests/nosub/JetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());

      std::vector<double> v0Pt;
      std::vector<int> v0Type;
      double ptjetsub = jet.pt();

      for (const auto& v0 : jet.template candidates_as<aod::CandidatesV0Data>()) {
        double z = v0.pt() / jet.pt();

        registry.fill(HIST("tests/nosub/JetPtEtaV0Pt"), jet.pt(), jet.eta(), v0.pt());
        registry.fill(HIST("tests/nosub/JetPtEtaV0Z"), jet.pt(), jet.eta(), z);

        if (v0.isK0SCandidate()) {
          registry.fill(HIST("tests/nosub/JetPtEtaK0SPt"), jet.pt(), jet.eta(), v0.pt());
          registry.fill(HIST("tests/nosub/JetPtEtaK0SZ"), jet.pt(), jet.eta(), z);
        }
        if (v0.isLambdaCandidate()) {
          registry.fill(HIST("tests/nosub/JetPtEtaLambdaPt"), jet.pt(), jet.eta(), v0.pt());
          registry.fill(HIST("tests/nosub/JetPtEtaLambdaZ"), jet.pt(), jet.eta(), z);
        }
        if (v0.isAntiLambdaCandidate()) {
          registry.fill(HIST("tests/nosub/JetPtEtaAntiLambdaPt"), jet.pt(), jet.eta(), v0.pt());
          registry.fill(HIST("tests/nosub/JetPtEtaAntiLambdaZ"), jet.pt(), jet.eta(), z);
        }

        if (gRandom->Uniform() > v0Fraction) { // Rejected V0
          ptjetsub -= v0.pt();
        } else { // Accepted V0
          v0Pt.push_back(v0.pt());
          if (v0.isK0SCandidate()) {
            v0Type.push_back(PDG_t::kK0Short);
          } else if (v0.isLambdaCandidate()) {
            v0Type.push_back(PDG_t::kLambda0);
          } else if (v0.isAntiLambdaCandidate()) {
            v0Type.push_back(PDG_t::kLambda0Bar);
          }
        }
      } // V0s in jet loop

      registry.fill(HIST("tests/sub/JetPtEtaPhi"), ptjetsub, jet.eta(), jet.phi());
      for (unsigned int i = 0; i < v0Pt.size(); ++i) {
        int type = v0Type[i];
        double pt = v0Pt[i];
        double z = pt / ptjetsub;

        registry.fill(HIST("tests/sub/JetPtEtaV0Pt"), ptjetsub, jet.eta(), pt);
        registry.fill(HIST("tests/sub/JetPtEtaV0Z"), ptjetsub, jet.eta(), z);

        if (type == PDG_t::kK0Short) {
          registry.fill(HIST("tests/sub/JetPtEtaK0SPt"), ptjetsub, jet.eta(), pt);
          registry.fill(HIST("tests/sub/JetPtEtaK0SZ"), ptjetsub, jet.eta(), z);
        } else if (type == PDG_t::kLambda0) {
          registry.fill(HIST("tests/sub/JetPtEtaLambdaPt"), ptjetsub, jet.eta(), pt);
          registry.fill(HIST("tests/sub/JetPtEtaLambdaZ"), ptjetsub, jet.eta(), z);
        } else if (type == PDG_t::kLambda0Bar) {
          registry.fill(HIST("tests/sub/JetPtEtaAntiLambdaPt"), ptjetsub, jet.eta(), pt);
          registry.fill(HIST("tests/sub/JetPtEtaAntiLambdaZ"), ptjetsub, jet.eta(), z);
        }
      } // Accepted V0s in jet loop
    } // Jets loop
  }
  PROCESS_SWITCH(V0QA, processTestSubtractedJetFinder, "Test subtracted jet finder", false);

  using DaughterJTracks = soa::Join<aod::JetTracks, aod::JTrackPIs>;
  using DaughterTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::TracksCov>;
  void processV0TrackQA(aod::JetCollision const& /*jcoll*/, aod::CandidatesV0Data const& v0s, DaughterJTracks const&, DaughterTracks const&)
  {
    //   if (!jetderiveddatautilities::selectCollision(jcoll, eventSelectionBits)) {
    //     return;
    //   }
    for (const auto& v0 : v0s) {
      if (v0.isRejectedCandidate())
        continue;

      fillTrackQa<DaughterJTracks, DaughterTracks>(v0);
    }
  }
  PROCESS_SWITCH(V0QA, processV0TrackQA, "V0 track QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0QA>(cfgc, TaskName{"jet-v0qa"})};
}
