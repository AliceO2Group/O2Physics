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

/// \brief QA task for V0s in the jets framework, based on the LF v0cascadesqa task
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

// V0 jets
using MCDV0Jets = aod::V0ChargedMCDetectorLevelJets;
using MCDV0JetsWithConstituents = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetConstituents>;
using MatchedMCDV0Jets = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;
using MatchedMCDV0JetsWithConstituents = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetConstituents, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;

using MCPV0Jets = aod::V0ChargedMCParticleLevelJets;
using MCPV0JetsWithConstituents = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetConstituents>;
using MatchedMCPV0Jets = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;
using MatchedMCPV0JetsWithConstituents = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetConstituents, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;

struct V0QA {
  HistogramRegistry registry{"registry"}; // CallSumw2 = false?

  Configurable<std::string> evSel{"evSel", "sel8WithoutTimeFrameBorderCut", "choose event selection"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.995, "Minimum V0 cosine of pointing angle"};
  Configurable<float> v0radiusMin{"v0radiusMin", 0.5, "Minimum V0 radius (cm)"};
  Configurable<float> dcav0dauMax{"dcav0dauMax", 1.0, "Maximum DCA between V0 daughters (cm)"};
  Configurable<float> dcapiMin{"dcapiMin", 0.1, "Minimum DCA of pion daughter to PV (cm)"};
  Configurable<float> dcaprMin{"dcaprMin", 0.1, "Minimum DCA of proton daughter to PV (cm)"};
  Configurable<float> yK0SMax{"yK0SMax", 0.5, "Maximum rapidity of K0S"};
  Configurable<float> yLambdaMax{"yLambdaMax", 0.5, "Maximum rapidity of Lambda(bar)"};
  Configurable<float> lifetimeK0SMax{"lifetimeK0SMax", 20.0, "Maximum lifetime of K0S (cm)"};
  Configurable<float> lifetimeLambdaMax{"lifetimeLambdaMax", 30.0, "Maximum lifetime of Lambda (cm)"};
  Configurable<float> yPartMax{"yPartMax", 0.5, "Maximum rapidity of particles"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "Vertex Z cut"};

  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  // PresliceUnsorted<JetCollisionsMCD> CollisionsPerMcCollision = aod::jv0indices::mcCollisionId;

  ConfigurableAxis binPtJet{"ptJet", {100., 0.0f, 50.0f}, ""};
  ConfigurableAxis binPtV0{"ptV0", {100., 0.0f, 50.0f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.0f, 1.0f}, ""};
  ConfigurableAxis binPhi{"binPhi", {static_cast<int>(TMath::Pi()) * 10 / 2, 0.0f, 2. * static_cast<int>(TMath::Pi())}, ""};

  ConfigurableAxis binInvMassK0S{"binInvMassK0S", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis binInvMassLambda{"binInvMassLambda", {200, 1.07f, 1.17f}, ""};
  ConfigurableAxis binV0Radius{"R", {100., 0.0f, 50.0f}, ""};

  int eventSelection = -1;

  void init(InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(evSel));

    const AxisSpec axisJetPt{binPtJet, "Jet Pt (GeV/c)"};
    const AxisSpec axisV0Pt{binPtV0, "V0 Pt (GeV/c)"};
    const AxisSpec axisEta{binEta, "Eta"};
    const AxisSpec axisPhi{binPhi, "Phi"};
    const AxisSpec axisV0Radius{binV0Radius, "V0 Radius (cm)"};
    const AxisSpec axisK0SM{binInvMassK0S, "M(#pi^{+} #pi^{-}) (GeV/c^{2})"};
    const AxisSpec axisLambdaM{binInvMassLambda, "M(p #pi^{-}) (GeV/c^{2})"};
    const AxisSpec axisAntiLambdaM{binInvMassLambda, "M(#bar{p} #pi^{+}) (GeV/c^{2})"};

    if (doprocessMcD) {
      registry.add("hEvents", "Events", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("K0SPtEtaMass", "K0S Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisK0SM});
      registry.add("InvMassK0STrue", "Invariant mass of K0S", HistType::kTH3D, {axisV0Pt, axisV0Radius, axisK0SM});
      registry.add("InvMassLambdaTrue", "Invariant mass of Lambda", HistType::kTH3D, {axisV0Pt, axisV0Radius, axisLambdaM});
      registry.add("LambdaPtEtaMass", "Lambda Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisLambdaM});
      registry.add("InvMassAntiLambdaTrue", "Invariant mass of AntiLambda", HistType::kTH3D, {axisV0Pt, axisV0Radius, axisAntiLambdaM});
      registry.add("AntiLambdaPtEtaMass", "AntiLambda Pt, Eta, Mass", HistType::kTH3D, {axisV0Pt, axisEta, axisAntiLambdaM});
    }
    if (doprocessMcP) {
      registry.add("hMcEvents", "MC Events", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("GeneratedK0S", "Generated K0S", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Radius});
      registry.add("GeneratedLambda", "Generated Lambda", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Radius});
      registry.add("GeneratedAntiLambda", "Generated AntiLambda", HistType::kTH3D, {axisV0Pt, axisEta, axisV0Radius});
    }
    if (doprocessMcDJets) {
      registry.add("hJetEvents", "Jet Events", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("JetPtEtaK0SPt", "Jet Pt, Eta, K0S Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("InvMassJetK0STrue", "Invariant mass of K0S in jets", HistType::kTH3D, {axisJetPt, axisV0Pt, axisK0SM});
      registry.add("JetPtEtaLambdaPt", "Jet Pt, Eta, Lambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("InvMassJetLambdaTrue", "Invariant mass of Lambda in jets", HistType::kTH3D, {axisJetPt, axisV0Pt, axisLambdaM});
      registry.add("JetPtEtaAntiLambdaPt", "Jet Pt, Eta, AntiLambda Pt", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("InvMassJetAntiLambdaTrue", "Invariant mass of AntiLambda in jets", HistType::kTH3D, {axisJetPt, axisV0Pt, axisAntiLambdaM});
    }
    if (doprocessMcDMatchedJets) {
      registry.add("hMatchedJetEvents", "Matched Jet Events", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("JetsPtEtaK0SPt", "Matched Jet Pt, Eta, K0S Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt});
      registry.add("InvMassJetsK0STrue", "Invariant mass of K0S in matched jets", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisV0Pt, axisK0SM});
      registry.add("JetsPtEtaLambdaPt", "Matched Jet Pt, Eta, Lambda Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt});
      registry.add("InvMassJetsLambdaTrue", "Invariant mass of Lambda in matched jets", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisV0Pt, axisLambdaM});
      registry.add("JetsPtEtaAntiLambdaPt", "Matched Jet Pt, Eta, AntiLambda Pt", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisEta, axisV0Pt});
      registry.add("InvMassJetsAntiLambdaTrue", "Invariant mass of AntiLambda in matched jets", HistType::kTHnSparseD, {axisJetPt, axisJetPt, axisV0Pt, axisAntiLambdaM});
    }
    if (doprocessMcPJets) {
      registry.add("hMcJetEvents", "MC Jet Events", {HistType::kTH1D, {{2, 0.0f, 2.0f}}});
      registry.add("GeneratedJetK0S", "Generated Jet K0S", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("GeneratedJetLambda", "Generated Jet Lambda", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
      registry.add("GeneratedJetAntiLambda", "Generated Jet AntiLambda", HistType::kTH3D, {axisJetPt, axisEta, axisV0Pt});
    }
  } // init

  template <typename T, typename U>
  bool isCollisionReconstructed(T const& collision, U const& eventSelection)
  {
    if (!collision.has_mcCollision()) {
      return false;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return false;
    }
    return true;
  }
  template <typename T, typename U, typename V>
  bool V0sAreMatched(T const& v0, U const& particle, V const& /*tracks*/)
  {
    auto negId = v0.template negTrack_as<V>().mcParticleId();
    auto posId = v0.template posTrack_as<V>().mcParticleId();
    auto daughters = particle.daughtersIds();
    return ((negId == daughters[0] && posId == daughters[1]) || (posId == daughters[0] && negId == daughters[1]));
  }

  void processDummy(CandidatesV0MCD const&) {}
  PROCESS_SWITCH(V0QA, processDummy, "Dummy process function turned on by default", true);

  void processMcD(soa::Filtered<JetCollisionsMCD>::iterator const& jcoll, JetMcCollisions const&, soa::Join<CandidatesV0MCD, aod::McV0Labels> const& v0s, aod::McParticles const&)
  {
    registry.fill(HIST("hEvents"), 0.5);
    if (!isCollisionReconstructed(jcoll, eventSelection)) {
      return;
    }
    registry.fill(HIST("hEvents"), 1.5);
    double weight = jcoll.mcCollision().weight();

    for (const auto& v0 : v0s) {
      if (!v0.has_mcParticle()) {
        continue;
      }
      int pdg = v0.mcParticle().pdgCode();

      if (v0.v0cosPA() < v0cospaMin)
        continue;
      if (v0.v0radius() < v0radiusMin)
        continue;
      if (v0.dcaV0daughters() > dcav0dauMax)
        continue;

      // K0S
      if (TMath::Abs(pdg) == 310) {
        if (TMath::Abs(v0.dcapostopv()) < dcapiMin)
          continue;
        if (TMath::Abs(v0.dcanegtopv()) < dcapiMin)
          continue;
        if (TMath::Abs(v0.yK0Short()) > yK0SMax)
          continue;
        float ctauK0S = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassK0Short;
        if (ctauK0S > lifetimeK0SMax)
          continue;

        registry.fill(HIST("K0SPtEtaMass"), v0.pt(), v0.eta(), v0.mK0Short(), weight);
        registry.fill(HIST("InvMassK0STrue"), v0.pt(), v0.v0radius(), v0.mK0Short(), weight);
      }
      // Lambda
      if (pdg == 3122) {
        if (TMath::Abs(v0.dcapostopv()) < dcaprMin)
          continue;
        if (TMath::Abs(v0.dcanegtopv()) < dcapiMin)
          continue;
        if (TMath::Abs(v0.yLambda()) > yLambdaMax)
          continue;
        float ctauLambda = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassLambda0;
        if (ctauLambda > lifetimeLambdaMax)
          continue;

        registry.fill(HIST("LambdaPtEtaMass"), v0.pt(), v0.eta(), v0.mLambda(), weight);
        registry.fill(HIST("InvMassLambdaTrue"), v0.pt(), v0.v0radius(), v0.mLambda(), weight);
      }
      if (pdg == -3122) {
        if (TMath::Abs(v0.dcapostopv()) < dcapiMin)
          continue;
        if (TMath::Abs(v0.dcanegtopv()) < dcaprMin)
          continue;
        if (TMath::Abs(v0.yLambda()) > yLambdaMax)
          continue;
        float ctauAntiLambda = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassLambda0Bar;
        if (ctauAntiLambda > lifetimeLambdaMax)
          continue;

        registry.fill(HIST("AntiLambdaPtEtaMass"), v0.pt(), v0.eta(), v0.mAntiLambda(), weight);
        registry.fill(HIST("InvMassAntiLambdaTrue"), v0.pt(), v0.v0radius(), v0.mAntiLambda(), weight);
      }
    }
  }
  PROCESS_SWITCH(V0QA, processMcD, "Reconstructed true V0s", false);

  void processMcP(JetMcCollision const& mccoll, CandidatesV0MCP const& pv0s, soa::SmallGroups<JetCollisionsMCD> const& collisions)
  {
    registry.fill(HIST("hMcEvents"), 0.5);
    bool isReconstructed = false;

    for (auto collision : collisions) {
      if (!isCollisionReconstructed(collision, eventSelection)) {
        continue;
      }
      if (collision.mcCollision().globalIndex() != mccoll.globalIndex()) {
        continue;
      }
      isReconstructed = true;
      break;
    }
    if (!isReconstructed) {
      return;
    }

    registry.fill(HIST("hMcEvents"), 1.5);
    double weight = mccoll.weight();

    for (auto& pv0 : pv0s) {
      if (!pv0.has_daughters())
        continue;
      if (!pv0.isPhysicalPrimary())
        continue;
      if (TMath::Abs(pv0.y() > yPartMax))
        continue;

      // Can calculate this from CandidatesV0MCD (contains decay vertex)
      double R_Decay = 1.0;

      if (pv0.pdgCode() == 310) {
        registry.fill(HIST("GeneratedK0S"), pv0.pt(), pv0.eta(), R_Decay, weight);
      }
      if (pv0.pdgCode() == 3122) {
        registry.fill(HIST("GeneratedLambda"), pv0.pt(), pv0.eta(), R_Decay, weight);
      }
      if (pv0.pdgCode() == -3122) {
        registry.fill(HIST("GeneratedAntiLambda"), pv0.pt(), pv0.eta(), R_Decay, weight);
      }
    }
  }
  PROCESS_SWITCH(V0QA, processMcP, "Particle level V0s", false);

  void processMcDJets(soa::Filtered<JetCollisionsMCD>::iterator const& jcoll, JetMcCollisions const&, MCDV0JetsWithConstituents const& jets, soa::Join<CandidatesV0MCD, aod::McV0Labels> const&, aod::McParticles const&)
  {
    registry.fill(HIST("hJetEvents"), 0.5);
    if (!isCollisionReconstructed(jcoll, eventSelection)) {
      return;
    }
    registry.fill(HIST("hJetEvents"), 1.5);
    double weight = jcoll.mcCollision().weight();

    for (const auto& jet : jets) {
      // if (!jetfindingutilities::isInEtaAcceptance(jet, -99., -99., v0EtaMin, v0EtaMax))
      for (const auto& v0 : jet.template candidates_as<soa::Join<CandidatesV0MCD, aod::McV0Labels>>()) {
        if (!v0.has_mcParticle()) {
          continue;
        }
        int pdg = v0.mcParticle().pdgCode();

        if (v0.v0cosPA() < v0cospaMin)
          continue;
        if (v0.v0radius() < v0radiusMin)
          continue;
        if (v0.dcaV0daughters() > dcav0dauMax)
          continue;

        // K0S
        if (TMath::Abs(pdg) == 310) {
          if (TMath::Abs(v0.dcapostopv()) < dcapiMin)
            continue;
          if (TMath::Abs(v0.dcanegtopv()) < dcapiMin)
            continue;
          if (TMath::Abs(v0.yK0Short()) > yK0SMax)
            continue;
          float ctauK0S = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassK0Short;
          if (ctauK0S > lifetimeK0SMax)
            continue;

          registry.fill(HIST("JetPtEtaK0SPt"), jet.pt(), jet.eta(), v0.pt(), weight);
          registry.fill(HIST("InvMassJetK0STrue"), jet.pt(), v0.pt(), v0.mK0Short(), weight);
        }
        // Lambda
        if (pdg == 3122) {
          if (TMath::Abs(v0.dcapostopv()) < dcaprMin)
            continue;
          if (TMath::Abs(v0.dcanegtopv()) < dcapiMin)
            continue;
          if (TMath::Abs(v0.yLambda()) > yLambdaMax)
            continue;
          float ctauLambda = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassLambda0;
          if (ctauLambda > lifetimeLambdaMax)
            continue;

          registry.fill(HIST("JetPtEtaLambdaPt"), jet.pt(), jet.eta(), v0.pt(), weight);
          registry.fill(HIST("InvMassJetLambdaTrue"), jet.pt(), v0.pt(), v0.mLambda(), weight);
        }
        if (pdg == -3122) {
          if (TMath::Abs(v0.dcapostopv()) < dcapiMin)
            continue;
          if (TMath::Abs(v0.dcanegtopv()) < dcaprMin)
            continue;
          if (TMath::Abs(v0.yLambda()) > yLambdaMax)
            continue;
          float ctauAntiLambda = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassLambda0Bar;
          if (ctauAntiLambda > lifetimeLambdaMax)
            continue;

          registry.fill(HIST("JetPtEtaAntiLambdaPt"), jet.pt(), jet.eta(), v0.pt(), weight);
          registry.fill(HIST("InvMassJetAntiLambdaTrue"), jet.pt(), v0.pt(), v0.mAntiLambda(), weight);
        }
      }
    }
  }
  PROCESS_SWITCH(V0QA, processMcDJets, "Reconstructed true V0s in jets", false);

  void processMcDMatchedJets(soa::Filtered<JetCollisionsMCD>::iterator const& jcoll, JetMcCollisions const&, MatchedMCDV0JetsWithConstituents const& mcdjets, MatchedMCPV0JetsWithConstituents const& mcpjets, soa::Join<CandidatesV0MCD, aod::McV0Labels> const&, CandidatesV0MCP const&, JetTracksMCD const& jTracks, aod::McParticles const&)
  {
    registry.fill(HIST("hMatchedJetEvents"), 0.5);
    if (!isCollisionReconstructed(jcoll, eventSelection)) {
      return;
    }
    registry.fill(HIST("hMatchedJetEvents"), 1.5);
    double weight = jcoll.mcCollision().weight();

    for (const auto& mcdjet : mcdjets) {
      // if (!jetfindingutilities::isInEtaAcceptance(mcdjet, -99., -99., v0EtaMin, v0EtaMax))
      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<MatchedMCPV0JetsWithConstituents>()) {
        for (const auto& v0 : mcdjet.template candidates_as<soa::Join<CandidatesV0MCD, aod::McV0Labels>>()) {
          if (!v0.has_mcParticle())
            continue;
          if (v0.v0cosPA() < v0cospaMin)
            continue;
          if (v0.v0radius() < v0radiusMin)
            continue;
          if (v0.dcaV0daughters() > dcav0dauMax)
            continue;

          for (const auto& pv0 : mcpjet.template candidates_as<CandidatesV0MCP>()) {
            if (!V0sAreMatched(v0, pv0, jTracks))
              continue;
            int pdg = pv0.pdgCode();

            // K0S
            if (TMath::Abs(pdg) == 310) {
              if (TMath::Abs(v0.dcapostopv()) < dcapiMin)
                continue;
              if (TMath::Abs(v0.dcanegtopv()) < dcapiMin)
                continue;
              if (TMath::Abs(v0.yK0Short()) > yK0SMax)
                continue;
              float ctauK0S = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassK0Short;
              if (ctauK0S > lifetimeK0SMax)
                continue;

              registry.fill(HIST("JetsPtEtaK0SPt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
              registry.fill(HIST("InvMassJetsK0STrue"), mcpjet.pt(), mcdjet.pt(), v0.pt(), v0.mK0Short(), weight);
            }
            // Lambda
            if (pdg == 3122) {
              if (TMath::Abs(v0.dcapostopv()) < dcaprMin)
                continue;
              if (TMath::Abs(v0.dcanegtopv()) < dcapiMin)
                continue;
              if (TMath::Abs(v0.yLambda()) > yLambdaMax)
                continue;
              float ctauLambda = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassLambda0;
              if (ctauLambda > lifetimeLambdaMax)
                continue;

              registry.fill(HIST("JetsPtEtaLambdaPt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
              registry.fill(HIST("InvMassJetsLambdaTrue"), mcpjet.pt(), mcdjet.pt(), v0.pt(), v0.mLambda(), weight);
            }
            if (pdg == -3122) {
              if (TMath::Abs(v0.dcapostopv()) < dcapiMin)
                continue;
              if (TMath::Abs(v0.dcanegtopv()) < dcaprMin)
                continue;
              if (TMath::Abs(v0.yLambda()) > yLambdaMax)
                continue;
              float ctauAntiLambda = v0.distovertotmom(jcoll.posX(), jcoll.posY(), jcoll.posZ()) * o2::constants::physics::MassLambda0Bar;
              if (ctauAntiLambda > lifetimeLambdaMax)
                continue;

              registry.fill(HIST("JetsPtEtaAntiLambdaPt"), mcpjet.pt(), mcdjet.pt(), mcdjet.eta(), v0.pt(), weight);
              registry.fill(HIST("InvMassJetsAntiLambdaTrue"), mcpjet.pt(), mcdjet.pt(), v0.pt(), v0.mAntiLambda(), weight);
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(V0QA, processMcDMatchedJets, "Reconstructed true V0s in jets", false);

  void processMcPJets(JetMcCollision const& mccoll, soa::SmallGroups<JetCollisionsMCD> const& collisions, MCPV0JetsWithConstituents const& jets, CandidatesV0MCP const&)
  {
    registry.fill(HIST("hMcJetEvents"), 0.5);
    bool isReconstructed = false;

    for (auto collision : collisions) {
      if (!isCollisionReconstructed(collision, eventSelection)) {
        continue;
      }
      if (collision.mcCollision().globalIndex() != mccoll.globalIndex()) {
        continue;
      }
      isReconstructed = true;
      break;
    }
    if (!isReconstructed) {
      return;
    }

    registry.fill(HIST("hMcJetEvents"), 1.5);
    double weight = mccoll.weight();

    for (auto& jet : jets) {
      // if (!jetfindingutilities::isInEtaAcceptance(jet, -99., -99., v0EtaMin, v0EtaMax))
      for (const auto& pv0 : jet.template candidates_as<CandidatesV0MCP>()) {
        if (!pv0.has_daughters())
          continue;
        if (!pv0.isPhysicalPrimary())
          continue;
        if (TMath::Abs(pv0.y() > yPartMax))
          continue; // TODO: Should actually check the jets

        if (pv0.pdgCode() == 310) {
          registry.fill(HIST("GeneratedJetK0S"), jet.pt(), jet.eta(), pv0.pt(), weight);
        }
        if (pv0.pdgCode() == 3122) {
          registry.fill(HIST("GeneratedJetLambda"), jet.pt(), jet.eta(), pv0.pt(), weight);
        }
        if (pv0.pdgCode() == -3122) {
          registry.fill(HIST("GeneratedJetAntiLambda"), jet.pt(), jet.eta(), pv0.pt(), weight);
        }
      }
    }
  }
  PROCESS_SWITCH(V0QA, processMcPJets, "Particle level V0s in jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0QA>(cfgc, TaskName{"jet-v0qa"})};
}
