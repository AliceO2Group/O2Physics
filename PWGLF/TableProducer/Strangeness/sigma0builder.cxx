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
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using std::cout;
using std::endl;

struct sigma0builder{
  Produces<aod::V0MLSigmaCandidates> v0MLSigmas; // save sigma0 candidates for analysis
  Produces<aod::V0MLSigmaMCCandidates> v0MLMCSigmas; // save sigma0 candidates for analysis

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas of sigma0 candidates"};
  Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas of sigma0 candidates"};
  // Axis	
  // base properties
  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};
  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {200, 1.16f, 1.23f}, "M_{#Sigma^{0}} (GeV/c^{2})"};

  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});
    histos.add("hEventVertexZMC", "hEventVertexZMC", kTH1F, {vertexZ});
    
    // Number of reconstructed sigma per collision as QA
    histos.add("hNSigmaCandidates", "hNSigmaCandidates", kTH1F, {{100,-0.5f, 99.5f}});
    histos.add("hNSigmaCandidatesMC", "hNSigmaCandidatesMC", kTH1F, {{100,-0.5f, 99.5f}});
    histos.add("hNGammaCandidatesMC", "hNGammaCandidatesMC", kTH1F, {{100,-0.5f, 99.5f}});
    histos.add("hNLambdaCandidatesMC", "hNLambdaCandidatesMC", kTH1F, {{100,-0.5f, 99.5f}});

    // Invariant Mass
    histos.add("InvMassSigmaMC", "InvMassSigmaMC", kTH1F, {axisSigmaMass});

    if(doprocessCounterQA){
      histos.add("hGammaIndices", "hGammaIndices", {HistType::kTH1F, {{4000,0.0f, 400000.0f}}});
      histos.add("hCollIndices", "hCollIndices", {HistType::kTH1F, {{4000,0.0f, 4000.0f}}});
      histos.add("h2dIndices", "h2dIndices", {HistType::kTH2F, {{4000,0.0f, 40000.0f},{4000,0.0f, 400000.0f}}});
    }
  }

  // Helper struct to pass v0 information
  struct {
    float mass;
    float pT;
  } sigmaCandidate;
  
  // Process sigma candidate and store properties in object
  template <typename TV0Object, typename TCollision>
  bool processSigmaCandidate(TCollision const& collision, TV0Object const& lambda, TV0Object const& gamma)
  {    
    // Gamma selection:
    if( gamma.gammaBDTScore() <= Gamma_MLThreshold) 
      return false;

    // Lambda selection:
    if( lambda.lambdaBDTScore() <= Lambda_MLThreshold) 
      return false;
      
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(),lambda.py(), lambda.pz()};
    auto arrMom = std::array{pVecPhotons, pVecLambda};
    sigmaCandidate.mass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    sigmaCandidate.pT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
    return true;
  }

  // This process function cross-checks index correctness 
  void processCounterQA(soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas> const& v0s)
  {
   for (auto& gamma: v0s) {
     histos.fill(HIST("hGammaIndices"), gamma.globalIndex());
     histos.fill(HIST("hCollIndices"), gamma.straCollisionId());
     histos.fill(HIST("h2dIndices"), gamma.straCollisionId(), gamma.globalIndex());
   }
 }

  void processMonteCarlo(aod::StraCollision const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas, aod::V0LambdaMLScores, aod::V0GammaMLScores> const& v0s)
  {
    int SigmaCounter = 0;
    int GammaCounter = 0;
    int LambdaCounter = 0;

    histos.fill(HIST("hEventVertexZMC"), coll.posZ());

    for (auto& gamma: v0s) { // selecting photons from Sigma0

      if ((gamma.pdgCode()!=22) || (gamma.pdgCodeMother()!=3212)) continue;
      GammaCounter++;
      for (auto& lambda: v0s){ // selecting lambdas from Sigma0
        if ((lambda.pdgCode()!=3122) || (lambda.pdgCodeMother()!=3212)) continue;
        LambdaCounter++;
        
        //if (gamma.motherMCPartId()!=lambda.motherMCPartId()) continue; 
        if(!processSigmaCandidate(coll, lambda,gamma)) continue;  
	      SigmaCounter ++;

        bool fIsSigma = (gamma.motherMCPartId()==lambda.motherMCPartId());

        // Sigma related
        float fSigmapT = sigmaCandidate.pT;
        float fSigmaMass = sigmaCandidate.mass; 
      
        // Daughters related
        /// Photon
        float fPhotonPt = gamma.pt();
        float fPhotonMass = gamma.mGamma();
        float fPhotonQt = gamma.qtarm();
        float fPhotonAlpha = gamma.alpha();
        float fPhotonRadius = gamma.v0radius(); 
        float fPhotonCosPA = gamma.v0cosPA();
        float fPhotonDCADau = gamma.dcaV0daughters();
        float fPhotonDCANegPV = gamma.dcanegtopv();
        float fPhotonDCAPosPV = gamma.dcapostopv();
        float fPhotonZconv = gamma.z(); 
      
        // Lambda
        float fLambdaPt = lambda.pt();
        float fLambdaMass = lambda.mLambda();
        float fLambdaQt = lambda.qtarm();
        float fLambdaAlpha = lambda.alpha();
        float fLambdaRadius = lambda.v0radius();
        float fLambdaCosPA = lambda.v0cosPA();
        float fLambdaDCADau = lambda.dcaV0daughters();
        float fLambdaDCANegPV = lambda.dcanegtopv();
        float fLambdaDCAPosPV = lambda.dcapostopv();

        // Filling TTree for ML analysis
        v0MLMCSigmas(fSigmapT, fSigmaMass, fPhotonPt, fPhotonMass, fPhotonQt, fPhotonAlpha, fPhotonRadius, fPhotonCosPA, fPhotonDCADau, fPhotonDCANegPV, fPhotonDCAPosPV, fPhotonZconv, fLambdaPt, fLambdaMass, fLambdaQt, fLambdaAlpha, fLambdaRadius, fLambdaCosPA, fLambdaDCADau, fLambdaDCANegPV, fLambdaDCAPosPV, fIsSigma, gamma.gammaBDTScore(), lambda.lambdaBDTScore());
        histos.fill(HIST("InvMassSigmaMC"), fSigmaMass);
      }
    }
    histos.fill(HIST("hNSigmaCandidatesMC"), SigmaCounter);
    histos.fill(HIST("hNGammaCandidatesMC"), GammaCounter);
    histos.fill(HIST("hNLambdaCandidatesMC"), LambdaCounter);
  }

    void processRealData(aod::StraCollision const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0LambdaMLScores, aod::V0GammaMLScores> const& v0s)
  {
    int SigmaCounter = 0;

    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    for (auto& gamma: v0s) { // selecting photons from Sigma0
      for (auto& lambda: v0s){ // selecting lambdas from Sigma0
        if(!processSigmaCandidate(coll, lambda,gamma))
          continue; 
        SigmaCounter ++; 
        
        //bool fIsSigma = false;

        // Sigma related
        float fSigmapT = sigmaCandidate.pT;
        float fSigmaMass = sigmaCandidate.mass; 
      
        // Daughters related
        /// Photon
        float fPhotonPt = gamma.pt();
        float fPhotonMass = gamma.mGamma();
        float fPhotonQt = gamma.qtarm();
        float fPhotonAlpha = gamma.alpha();
        float fPhotonRadius = gamma.v0radius(); 
        float fPhotonCosPA = gamma.v0cosPA();
        float fPhotonDCADau = gamma.dcaV0daughters();
        float fPhotonDCANegPV = gamma.dcanegtopv();
        float fPhotonDCAPosPV = gamma.dcapostopv();
        float fPhotonZconv = gamma.z(); 
      
        // Lambda
        float fLambdaPt = lambda.pt();
        float fLambdaMass = lambda.mLambda();
        float fLambdaQt = lambda.qtarm();
        float fLambdaAlpha = lambda.alpha();
        float fLambdaRadius = lambda.v0radius();
        float fLambdaCosPA = lambda.v0cosPA();
        float fLambdaDCADau = lambda.dcaV0daughters();
        float fLambdaDCANegPV = lambda.dcanegtopv();
        float fLambdaDCAPosPV = lambda.dcapostopv();

        // Filling TTree for ML analysis
        v0MLSigmas(fSigmapT, fSigmaMass, fPhotonPt, fPhotonMass, fPhotonQt, fPhotonAlpha, fPhotonRadius, fPhotonCosPA, fPhotonDCADau, fPhotonDCANegPV, fPhotonDCAPosPV, fPhotonZconv, fLambdaPt, fLambdaMass, fLambdaQt, fLambdaAlpha, fLambdaRadius, fLambdaCosPA, fLambdaDCADau, fLambdaDCANegPV, fLambdaDCAPosPV, gamma.gammaBDTScore(), lambda.lambdaBDTScore());

      }
    }
    histos.fill(HIST("hNSigmaCandidates"), SigmaCounter);
  }

  PROCESS_SWITCH(sigma0builder, processCounterQA, "Check standard counter correctness", true);
  PROCESS_SWITCH(sigma0builder, processMonteCarlo, "Do Monte-Carlo-based analysis", true);
  PROCESS_SWITCH(sigma0builder, processRealData, "Do real data analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigma0builder>(cfgc)};
  
}