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
//
// This is a task that employs the standard V0 tables and attempts to combine
// two V0s into a Sigma0 -> Lambda + gamma candidate.

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFSigmaTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CCDB/BasicCCDBManager.h"
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using std::cout;
using std::endl;
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct sigmaanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  // Configurable<bool> fSaveTTree{"fSaveTTree", false, "Save TTree with Sigma0 candidates info for ML-based analysis"};
  // Configurable<bool> analyseSigma{"analyseSigma", false, "process Sigma-like candidates"};
  // Configurable<bool> analyseAntiSigma{"analyseAntiSigma", false, "process AntiSigma-like candidates"};
  // Base selection criteria

  // Selection criteria: acceptance
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
  Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

  // Photon standard criteria:
  Configurable<float> minPhotonRadius{"minPhotonRadius", 5.0, "minimum photon conversion radius (cm)"};
  Configurable<float> maxPhotonRadius{"maxPhotonRadius", 180, "maximum photon conversion radius (cm)"};
  Configurable<float> minPhotonPt{"minPhotonPt", 0.02, "minimum photon pT (GeV/c)"};
  Configurable<float> minPhotonEta{"minPhotonEta", 0.9, "minimum photon eta"};
  Configurable<float> maxPhotonMass{"maxPhotonMass", 0.1, "Maximum photon mass (GeV/c^{2})"};
  Configurable<float> maxPhotonqt{"maxPhotonqt", 0.06, "Maximum photon qt value (AP plot) (GeV/c)"};
  Configurable<float> Photonalpha{"Photonalpha", 0.95, "Max photon alpha absolute value (AP plot)"};

  // PID (TPC)
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 4, "TpcPidNsigmaCut"};
  Configurable<bool> allowTPConly{"allowTPConly", false, "Accept V0s that are TPC-only"};

  // Track quality
  Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};

  // Lambda:
  Configurable<float> lambdaWindow{"lambdaWindow", 0.01, "Accept +/- this wrt Lambda mass (GeV/c^{2})"};
  Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
  Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};

  // Axis
  // base properties
  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisRadius{"axisRadius", {200, 0.0f, 100.0f}, "V0 radius (cm)"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {200, 1.16f, 1.23f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.08f, 1.16f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, -0.1f, 0.1f}, "M_{#Gamma}"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // ML
  ConfigurableAxis MLProb{"MLOutput", {100, 0.0f, 1.0f}, ""};

  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});
    histos.add("hEventVertexZMC", "hEventVertexZMC", kTH1F, {vertexZ});

    // Number of reconstructed sigma per collision as QA
    histos.add("hNSigmaCandidates", "hNSigmaCandidates", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("hNSigmaCandidatesMC", "hNSigmaCandidatesMC", kTH1F, {{100, -0.5f, 99.5f}});

    // V0 Radius
    histos.add("h2dLambdaRadiusVsPt", "hLambdaRadiusVsPt", {HistType::kTH2F, {axisPt, axisRadius}});
    histos.add("h2dPhotonRadiusVsPt", "hPhotonRadiusVsPt", {HistType::kTH2F, {axisPt, axisRadius}});
    histos.add("h2dLambdaRadiusVsPtMC", "hLambdaRadiusVsPtMC", {HistType::kTH2F, {axisPt, axisRadius}});
    histos.add("h2dPhotonRadiusVsPtMC", "hPhotonRadiusVsPtMC", {HistType::kTH2F, {axisPt, axisRadius}});

    // Invariant Mass
    histos.add("h2dSigmaMassVsPt", "hSigmaMassVsPt", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("h2dLambdaMassVsPt", "hLambdaMassVsPt", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("h2dPhotonMassVsPt", "hPhotonMassVsPt", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("h2dSigmaMassVsPtMC", "hSigmaMassVsPtMC", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("h2dLambdaMassVsPtMC", "hLambdaMassVsPtMC", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("h2dPhotonMassVsPtMC", "hPhotonMassVsPtMC", {HistType::kTH2F, {axisPt, axisPhotonMass}});

    // Exploratory Analysis with MC:
    /// Armenteros-Polanski plot:
    histos.add("h2dMCArmenterosPolanski", "h2dMCArmenterosPolanski", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});

    /// Lambda
    histos.add("hMCLambdaCosPA", "hMCLambdaCosPA", kTH1F, {{100, 0.9f, 1.0f}});
    histos.add("hMCLambdaDCA_V0Dau", "hMCLambdaDCA_V0Dau", kTH1F, {{100, 0.0f, 3.0f}});
    histos.add("hMCLambdaDCA_V0Pos", "hMCLambdaDCA_V0Pos", kTH1F, {{100, 0.0f, 2.0f}});
    histos.add("hMCLambdaDCA_V0Neg", "hMCLambdaDCA_V0Neg", kTH1F, {{100, 0.0f, 2.0f}});
    histos.add("hMCLambda_V0Radius", "hMCLambda_V0Radius", kTH1F, {{100, 0.0f, 40.0f}});

    /// Photon:
    histos.add("hMCPhoton_ConversionRadius", "hMCPhoton_ConversionRadius", kTH1F, {{100, 0.0f, 150.0f}});
    histos.add("hMCPhotonCosPA", "hMCPhotonCosPA", kTH1F, {{100, 0.9f, 1.0f}});
    histos.add("hMCPhotonDCA_V0Dau", "hMCPhotonDCA_V0Dau", kTH1F, {{100, 0.0f, 5.0}});
    histos.add("hMCPhotonDCA_V0Pos", "hMCPhotonDCA_V0Pos", kTH1F, {{100, 0.0f, 5.0f}});
    histos.add("hMCPhotonDCA_V0Neg", "hMCPhotonDCA_V0Neg", kTH1F, {{100, 0.0f, 5.0f}});

    // ML Analysis
    histos.add("hMLOutputLambda", "hMLOutputLambda", kTH1F, {MLProb});
    histos.add("hMLOutputGamma", "hMLOutputGamma", kTH1F, {MLProb});

    if (doprocessCounterQA) {
      histos.add("hGammaIndices", "hGammaIndices", {HistType::kTH1F, {{4000, 0.0f, 400000.0f}}});
      histos.add("hCollIndices", "hCollIndices", {HistType::kTH1F, {{4000, 0.0f, 4000.0f}}});
      histos.add("h2dIndices", "h2dIndices", {HistType::kTH2F, {{4000, 0.0f, 40000.0f}, {4000, 0.0f, 400000.0f}}});
    }
  }

  // Helper struct to pass v0 information
  struct {
    float mass;
    float pT;
  } sigmaCandidate;

  // Process sigma candidate and store properties in object
  template <typename TV0Object, typename TCollision>
  bool processSigmaCandidate(TCollision const&, TV0Object const& lambda, TV0Object const& gamma)
  {
    // FIXME: this should be at the single particle level, preferably partitions
    if (!allowTPConly && (lambda.v0Type() > 1 || gamma.v0Type() > 1))
      return false;

    // Gamma selection criteria:
    if (gamma.mGamma() > maxPhotonMass)
      return false;
    if (gamma.v0radius() < minPhotonRadius)
      return false;
    if (gamma.v0radius() > maxPhotonRadius)
      return false;
    if (gamma.pt() < minPhotonPt)
      return false;
    if (gamma.qtarm() > maxPhotonqt)
      return false;
    if (TMath::Abs(gamma.alpha()) > Photonalpha)
      return false;

    // Lambda selection criteria:
    if (TMath::Abs(lambda.mLambda() - 1.115683) > lambdaWindow)
      return false;
    if (lambda.v0radius() < v0radius)
      return false;
    if (lambda.v0cosPA() < v0cospa)
      return false;
    if (lambda.dcapostopv() < dcapostopv)
      return false;
    if (lambda.dcanegtopv() < dcanegtopv)
      return false;
    if (lambda.dcaV0daughters() > dcav0dau)
      return false;

    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
    auto arrMom = std::array{pVecPhotons, pVecLambda};
    sigmaCandidate.mass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    sigmaCandidate.pT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
    return true;
  }

  // This process function cross-checks index correctness
  void processCounterQA(soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas> const& v0s)
  {
    for (auto& gamma : v0s) {
      histos.fill(HIST("hGammaIndices"), gamma.globalIndex());
      histos.fill(HIST("hCollIndices"), gamma.straCollisionId());
      histos.fill(HIST("h2dIndices"), gamma.straCollisionId(), gamma.globalIndex());
    }
  }

  void processMonteCarlo(aod::StraCollision const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas> const& v0s)
  {
    int SigmaCounter = 0;

    histos.fill(HIST("hEventVertexZMC"), coll.posZ());

    for (auto& gamma : v0s) { // selecting photons from Sigma0

      if ((gamma.pdgCode() != 22 || gamma.pdgCodeMother() != 3212) && (gamma.pdgCode() != -22 || gamma.pdgCodeMother() != -3212))
        continue;

      for (auto& lambda : v0s) { // selecting lambdas from Sigma0
        if ((lambda.pdgCode() != 3122 || lambda.pdgCodeMother() != 3212) && (lambda.pdgCode() != -3122 || lambda.pdgCodeMother() != -3212))
          continue;
        if (gamma.motherMCPartId() != lambda.motherMCPartId())
          continue; // selecting pair from exactly the same mother

        // Exploratory Analysis histograms:

        // Lambda
        histos.fill(HIST("hMCLambdaCosPA"), lambda.v0cosPA());
        histos.fill(HIST("hMCLambdaDCA_V0Dau"), TMath::Abs(lambda.dcaV0daughters()));
        histos.fill(HIST("hMCLambdaDCA_V0Pos"), TMath::Abs(lambda.dcapostopv()));
        histos.fill(HIST("hMCLambdaDCA_V0Neg"), TMath::Abs(lambda.dcanegtopv()));
        histos.fill(HIST("hMCLambda_V0Radius"), lambda.v0radius());

        // Photon
        histos.fill(HIST("hMCPhoton_ConversionRadius"), gamma.v0radius());
        histos.fill(HIST("hMCPhotonCosPA"), gamma.v0cosPA());
        histos.fill(HIST("hMCPhotonDCA_V0Dau"), TMath::Abs(gamma.dcaV0daughters()));
        histos.fill(HIST("hMCPhotonDCA_V0Pos"), TMath::Abs(gamma.dcapostopv()));
        histos.fill(HIST("hMCPhotonDCA_V0Neg"), TMath::Abs(gamma.dcanegtopv()));

        // Armenteros-Polanski plot:
        histos.fill(HIST("h2dMCArmenterosPolanski"), gamma.alpha(), gamma.qtarm());
        histos.fill(HIST("h2dMCArmenterosPolanski"), lambda.alpha(), lambda.qtarm());

        if (!processSigmaCandidate(coll, lambda, gamma))
          continue;
        SigmaCounter++;

        histos.fill(HIST("h2dLambdaRadiusVsPtMC"), sigmaCandidate.pT, lambda.v0radius());
        histos.fill(HIST("h2dPhotonRadiusVsPtMC"), sigmaCandidate.pT, gamma.v0radius());

        // Inv Mass
        histos.fill(HIST("h2dLambdaMassVsPtMC"), sigmaCandidate.pT, lambda.mLambda());
        histos.fill(HIST("h2dPhotonMassVsPtMC"), sigmaCandidate.pT, gamma.mGamma());
        histos.fill(HIST("h2dSigmaMassVsPtMC"), sigmaCandidate.pT, sigmaCandidate.mass);
      }
    }
    histos.fill(HIST("hNSigmaCandidatesMC"), SigmaCounter);
  }

  void processRealData(aod::StraCollision const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0LambdaMLScores, aod::V0GammaMLScores> const& v0s)
  {
    int SigmaCounter = 0;

    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    for (auto& gamma : v0s) {    // selecting photons from Sigma0
      for (auto& lambda : v0s) { // selecting lambdas from Sigma0
        if (!processSigmaCandidate(coll, lambda, gamma))
          continue;
        SigmaCounter++;

        histos.fill(HIST("h2dLambdaRadiusVsPt"), sigmaCandidate.pT, lambda.v0radius());
        histos.fill(HIST("h2dPhotonRadiusVsPt"), sigmaCandidate.pT, gamma.v0radius());

        // Inv Mass
        histos.fill(HIST("h2dLambdaMassVsPt"), sigmaCandidate.pT, lambda.mLambda());
        histos.fill(HIST("h2dPhotonMassVsPt"), sigmaCandidate.pT, gamma.mGamma());
        histos.fill(HIST("h2dSigmaMassVsPt"), sigmaCandidate.pT, sigmaCandidate.mass);

        histos.fill(HIST("hMLOutputLambda"), lambda.lambdaBDTScore());
        histos.fill(HIST("hMLOutputGamma"), gamma.gammaBDTScore());
      }
    }
    histos.fill(HIST("hNSigmaCandidates"), SigmaCounter);
  }

  PROCESS_SWITCH(sigmaanalysis, processCounterQA, "Check standard counter correctness", true);
  PROCESS_SWITCH(sigmaanalysis, processMonteCarlo, "Do Monte-Carlo-based analysis", true);
  PROCESS_SWITCH(sigmaanalysis, processRealData, "Do real data analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigmaanalysis>(cfgc)};
}
