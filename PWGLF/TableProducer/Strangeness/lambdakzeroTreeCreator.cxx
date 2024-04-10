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
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//  V0 ML Tree Creator task
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//
//  This task loops over a set of V0 indices and
//  creates a TTree for ML training and testing.
//
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//


#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
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
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

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

struct lambdakzeroTreeCreator{
    Produces<aod::V0MLCandidates> v0MLCandidates;
    HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

    // selection switches: please, use just one at a time 
    Configurable<bool> fGetAllCandidates{"fGetAllCandidates", true, "If True, create a Tree containing all available candidates"};
    Configurable<bool> fGetLambdaOnly{"fGetLambdaOnly", false, "If True, apply cuts to select only lambda signal and bkg candidates"};
    Configurable<bool> fGetAntiLambdaOnly{"fGetAntiLambdaOnly", false, "If True, apply cuts to select only lambda signal and bkg candidates"};
    Configurable<bool> fGetGammaOnly{"fGetGammaOnly", false, "If True, apply cuts to select only Gamma signal and bkg candidates"};
    Configurable<bool> fGetK0ShortOnly{"fGetK0ShortOnly", false, "If True, apply cuts to select only K0Short signal and bkg candidates"};
    
    // Base selection criteria

    /// Selection criteria: acceptance
    Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
    Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

    // PID (TPC)
    Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 4, "TpcPidNsigmaCut"};
    Configurable<bool> allowTPConly{"allowTPConly", false, "Accept V0s that are TPC-only"};
    
    // Track quality
    Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};

    // Lambda standard criteria::
    //Configurable<float> lambdaWindow{"lambdaWindow", 0.01, "Accept +/- this wrt Lambda mass (GeV/c^{2})"};
    Configurable<float> Lambdav0cospa{"Lambdav0cospa", 0.90, "min V0 CosPA"};
    Configurable<float> Lambdadcav0dau{"Lambdadcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> Lambdadcanegtopv{"Lambdadcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> Lambdadcapostopv{"Lambdadcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> Lambdav0radius{"Lambdav0radius", 1.5, "minimum V0 radius (cm)"};

    // Anti-Lambda standard criteria::
    Configurable<float> AntiLambdav0cospa{"AntiLambdav0cospa", 0.90, "min V0 CosPA"};
    Configurable<float> AntiLambdadcav0dau{"AntiLambdadcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> AntiLambdadcanegtopv{"AntiLambdadcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> AntiLambdadcapostopv{"AntiLambdadcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> AntiLambdav0radius{"AntiLambdav0radius", 1.5, "minimum V0 radius (cm)"};

    // Photon standard criteria:
    Configurable<float> PhotonMinRadius{"PhotonMinRadius", 1.0, "minimum photon conversion radius (cm)"};
    Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 200, "maximum photon conversion radius (cm)"};
    Configurable<float> PhotonMinPt{"PhotonMinPt", 0.001, "minimum photon pT (GeV/c)"};
    Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.2, "Maximum photon mass (GeV/c^{2})"};
    Configurable<float> PhotonMaxqt{"PhotonMaxqt", 0.1, "Maximum photon qt value (AP plot) (GeV/c)"};
    Configurable<float> PhotonMaxalpha{"PhotonMaxalpha", 1.0, "Max photon alpha absolute value (AP plot)"};

    // TODO: Include here K0Short standard criteria

    // Axis:
    ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};

    void init(InitContext const&)
    {
      histos.add("hEventVertexZMC", "hEventVertexZMC", kTH1F, {vertexZ});
    }

    // Helper struct to pass v0 information
  struct {
    float LambdaMass;
    float AntiLambdaMass;
    float GammaMass;
    float KZeroShortMass;
    float pT;
    float qt;
    float alpha;
    float v0radius;
    float v0cosPA;
    float dcapostopv;
    float dcanegtopv;
    float dcaV0daughters;
    bool isLambda;
    bool isAntiLambda;
    bool isGamma;
    bool isKZeroShort;
  } Candidate;
  
  // Process candidate and store properties in object
  template <typename TV0Object>
  bool processCandidate(TV0Object const& cand)
  {    
    Candidate.LambdaMass = cand.mLambda();
    Candidate.AntiLambdaMass = cand.mAntiLambda();
    Candidate.GammaMass = cand.mGamma();
    Candidate.KZeroShortMass = cand.mK0Short();
    Candidate.pT = cand.pt();
    Candidate.qt = cand.qtarm();
    Candidate.alpha = cand.alpha();
    Candidate.v0radius = cand.v0radius();
    Candidate.v0cosPA = cand.v0cosPA();
    Candidate.dcapostopv = cand.dcapostopv();
    Candidate.dcanegtopv = cand.dcanegtopv();
    Candidate.dcaV0daughters = cand.dcaV0daughters();
    Candidate.isLambda = (cand.pdgCode()==3122);
    Candidate.isAntiLambda = (cand.pdgCode()==-3122);
    Candidate.isGamma = (cand.pdgCode()==22);
    Candidate.isKZeroShort = (cand.pdgCode()==310);

    return true;
  }

  // Process lambda candidate 
  template <typename TV0Object>
  bool processLambdaCandidate(TV0Object const& lambda)
  {    
    // FIXME: there are smarter ways to perform this selection
    // Lambda base selection criteria:

    if( lambda.v0radius() < Lambdav0radius) 
      return false;
    if( lambda.v0cosPA() < Lambdav0cospa) 
      return false;
    if(TMath::Abs(lambda.dcapostopv()) < Lambdadcapostopv) 
      return false;
    if(TMath::Abs(lambda.dcanegtopv()) < Lambdadcanegtopv) 
      return false;
    if( lambda.dcaV0daughters() > Lambdadcav0dau) 
      return false;
  
    return processCandidate(lambda);
  }

  // Process antilambda candidate 
  template <typename TV0Object>
  bool processAntiLambdaCandidate(TV0Object const& antilambda)
  {    
    // FIXME: there are smarter ways to perform this selection
    // AntiLambda base selection criteria:

    if( antilambda.v0radius() < AntiLambdav0radius) 
      return false;
    if( antilambda.v0cosPA() < AntiLambdav0cospa) 
      return false;
    if(TMath::Abs(antilambda.dcapostopv()) < AntiLambdadcapostopv) 
      return false;
    if(TMath::Abs(antilambda.dcanegtopv()) < AntiLambdadcanegtopv) 
      return false;
    if( antilambda.dcaV0daughters() > AntiLambdadcav0dau) 
      return false;
  
    return processCandidate(antilambda);
  }

  // Process gamma candidate
  template <typename TV0Object>
  bool processGammaCandidate(TV0Object const& gamma)
  {    
    // FIXME: there are smarter ways to perform this selection
    // Gamma selection criteria:
    if( gamma.mGamma() > PhotonMaxMass ) 
      return false;
    if( gamma.v0radius() < PhotonMinRadius ) 
      return false;
    if( gamma.v0radius() > PhotonMaxRadius ) 
      return false;
    if( gamma.pt() < PhotonMinPt) 
      return false;
    if( gamma.qtarm() > PhotonMaxqt) 
      return false;
    if(TMath::Abs(gamma.alpha()) > PhotonMaxalpha) 
      return false;
    
    return processCandidate(gamma);
  }

  // Process k0short candidate
  template <typename TV0Object>
  bool processK0ShortCandidate(TV0Object const& kzero)
  {    
    // FIXME: there are smarter ways to perform this selection
    // TODO: Include KZeroShort selection criteria
    
    return processCandidate(kzero);
  }

  void process(aod::StraCollision const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas> const& v0s)
  {
    histos.fill(HIST("hEventVertexZMC"), coll.posZ());
    for (auto& cand: v0s){ // looping over lambdas 

      if(fGetLambdaOnly){ 
        if (!processLambdaCandidate(cand))
           continue;
      }
      if(fGetAntiLambdaOnly){ 
        if (!processAntiLambdaCandidate(cand))
           continue;
      }
      if(fGetGammaOnly){
        if (!processGammaCandidate(cand))
           continue;
      }
      if(fGetK0ShortOnly){
        if(!processK0ShortCandidate(cand))
           continue;
      }
      if(fGetAllCandidates) { 
         if(!processCandidate(cand))
            continue;
      }

      // Filling TTree for ML analysis
      v0MLCandidates(Candidate.LambdaMass, Candidate.AntiLambdaMass, Candidate.GammaMass, Candidate.KZeroShortMass, 
      Candidate.pT, Candidate.qt, Candidate.alpha, Candidate.v0radius, Candidate.v0cosPA,
      Candidate.dcapostopv, Candidate.dcanegtopv, Candidate.dcaV0daughters, 
      Candidate.isLambda, Candidate.isAntiLambda, Candidate.isGamma, Candidate.isKZeroShort);
  }

}

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdakzeroTreeCreator>(cfgc)};
  
}