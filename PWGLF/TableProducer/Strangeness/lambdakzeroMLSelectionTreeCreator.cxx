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

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct lambdakzeroMLSelectionTreeCreator{
    Produces<aod::V0MLCandidates> v0MLCandidates;
    HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

    // selection switches: please, use just one at a time! 
    Configurable<bool> fGetAllCandidates{"fGetAllCandidates", true, "If True, create a Tree containing all available candidates"};
    Configurable<bool> fGetLambdaOnly{"fGetLambdaOnly", false, "If True, apply cuts to select only lambda signal and bkg candidates"};
    Configurable<bool> fGetAntiLambdaOnly{"fGetAntiLambdaOnly", false, "If True, apply cuts to select only antilambda signal and bkg candidates"};
    Configurable<bool> fGetGammaOnly{"fGetGammaOnly", false, "If True, apply cuts to select only Gamma signal and bkg candidates"};
    Configurable<bool> fGetK0ShortOnly{"fGetK0ShortOnly", false, "If True, apply cuts to select only K0Short signal and bkg candidates"};
    
    // Base selection criteria

    // Lambda standard criteria::
    Configurable<float> Lambdav0cospa{"Lambdav0cospa", 0.90, "min V0 CosPA"};
    Configurable<float> Lambdadcav0dau{"Lambdadcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> Lambdadcanegtopv{"Lambdadcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> Lambdadcapostopv{"Lambdadcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> Lambdav0radius{"Lambdav0radius", 1.5, "minimum V0 radius (cm)"};
    Configurable<float> LambdaWindow{"LambdaWindow", 0.01, "Mass window around expected (in GeV/c2)"};

    // Anti-Lambda standard criteria::
    Configurable<float> AntiLambdav0cospa{"AntiLambdav0cospa", 0.90, "min V0 CosPA"};
    Configurable<float> AntiLambdadcav0dau{"AntiLambdadcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> AntiLambdadcanegtopv{"AntiLambdadcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> AntiLambdadcapostopv{"AntiLambdadcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> AntiLambdav0radius{"AntiLambdav0radius", 1.5, "minimum V0 radius (cm)"};
    Configurable<float> AntiLambdaWindow{"AntiLambdaWindow", 0.01, "Mass window around expected (in GeV/c2)"};

    // Photon standard criteria:
    Configurable<float> PhotonMinRadius{"PhotonMinRadius", 1.0, "minimum photon conversion radius (cm)"};
    Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 200, "maximum photon conversion radius (cm)"};
    Configurable<float> PhotonMinPt{"PhotonMinPt", 0.001, "minimum photon pT (GeV/c)"};
    Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.2, "Maximum photon mass (GeV/c^{2})"};
    Configurable<float> PhotonMaxqt{"PhotonMaxqt", 0.1, "Maximum photon qt value (AP plot) (GeV/c)"};
    Configurable<float> PhotonMaxalpha{"PhotonMaxalpha", 1.0, "Max photon alpha absolute value (AP plot)"};
    Configurable<float> PhotonWindow{"PhotonWindow", 0.01, "Mass window around expected (in GeV/c2)"};

    // KZeroShort standard criteria:
    // TODO: update criteria
    Configurable<float> K0Shortv0cospa{"K0Shortv0cospa", 0.90, "min V0 CosPA"};
    Configurable<float> K0Shortdcav0dau{"K0Shortdcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> K0Shortdcanegtopv{"K0Shortdcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> K0Shortdcapostopv{"K0Shortdcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> K0Shortv0radius{"K0Shortv0radius", 1.5, "minimum V0 radius (cm)"};
    Configurable<float> K0ShortWindow{"K0ShortWindow", 0.01, "Mass window around expected (in GeV/c2)"};

    // Axis:
    ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};

    void init(InitContext const&)
    {
      histos.add("hEventVertexZMC", "hEventVertexZMC", kTH1F, {vertexZ});
    }

    // Helper struct to pass v0 information
  struct {
  int posITSCls;
  int negITSCls;
  int posITSClSize;
  int negITSClSize;
  float posTPCRows;
  float negTPCRows;
  float posTPCSigmaPi;
  float negTPCSigmaPi;
  float posTPCSigmaPr;
  float negTPCSigmaPr;
  float posTPCSigmaEl;
  float negTPCSigmaEl;
  float TOFSigmaLaPr;
  float TOFSigmaLaPi;
  float TOFSigmaALaPi;
  float TOFSigmaALaPr;
  float TOFSigmaK0PiPlus;
  float TOFSigmaK0PiMinus;
  float LambdaMass;
  float AntiLambdaMass;
  float GammaMass;
  float KZeroShortMass;
  float pT;
  float qt;
  float alpha;
  float posEta;
  float negEta;
  float v0Eta;
  float Z;
  float v0radius;
  float PA;
  float dcapostopv;
  float dcanegtopv;
  float dcaV0daughters;
  float dcav0topv;
  float PsiPair;
  uint8_t v0type;
  bool isLambda;
  bool isAntiLambda;
  bool isGamma;
  bool isKZeroShort; 
  } Candidate;
  
  // Process candidate and store properties in object
  template <typename TV0Object>
  bool processCandidate(TV0Object const& cand)
  {    
    auto posTrackExtra = cand.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = cand.template negTrackExtra_as<dauTracks>();

    // Track quality
    Candidate.posITSCls = posTrackExtra.itsNCls();
    Candidate.negITSCls = negTrackExtra.itsNCls();
    Candidate.posITSClSize = posTrackExtra.itsClusterSizes();
    Candidate.negITSClSize = negTrackExtra.itsClusterSizes();
    Candidate.posTPCRows = posTrackExtra.tpcCrossedRows();
    Candidate.negTPCRows = negTrackExtra.tpcCrossedRows();

    // TPC PID
    Candidate.posTPCSigmaPi = posTrackExtra.tpcNSigmaPi();
    Candidate.negTPCSigmaPi = negTrackExtra.tpcNSigmaPi();
    Candidate.posTPCSigmaPr = posTrackExtra.tpcNSigmaPr();
    Candidate.negTPCSigmaPr = negTrackExtra.tpcNSigmaPr();
    Candidate.posTPCSigmaEl = posTrackExtra.tpcNSigmaEl();
    Candidate.negTPCSigmaEl = negTrackExtra.tpcNSigmaEl();

    // TOF PID: 
    Candidate.TOFSigmaLaPr = cand.tofNSigmaLaPr();
    Candidate.TOFSigmaLaPi = cand.tofNSigmaLaPi();
    Candidate.TOFSigmaALaPi = cand.tofNSigmaALaPi();
    Candidate.TOFSigmaALaPr = cand.tofNSigmaALaPr();
    Candidate.TOFSigmaK0PiPlus = cand.tofNSigmaK0PiPlus();
    Candidate.TOFSigmaK0PiMinus = cand.tofNSigmaK0PiMinus();

    // General  
    Candidate.LambdaMass = cand.mLambda();
    Candidate.AntiLambdaMass = cand.mAntiLambda();
    Candidate.GammaMass = cand.mGamma();
    Candidate.KZeroShortMass = cand.mK0Short();
    Candidate.pT = cand.pt();
    Candidate.qt = cand.qtarm();
    Candidate.alpha = cand.alpha();
    Candidate.posEta = cand.positiveeta();
    Candidate.negEta = cand.negativeeta();
    Candidate.v0Eta = cand.eta();
    Candidate.Z = cand.z();

    // Topological 
    Candidate.v0radius = cand.v0radius();
    Candidate.PA = TMath::ACos(cand.v0cosPA());
    Candidate.dcapostopv = cand.dcapostopv();
    Candidate.dcanegtopv = cand.dcanegtopv();
    Candidate.dcaV0daughters = cand.dcaV0daughters();
    Candidate.dcav0topv = cand.dcav0topv();
    Candidate.PsiPair = cand.psipair();

    // Debug  
    Candidate.v0type = cand.v0Type();

    // MC flags
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
    if ((std::abs(lambda.mLambda() - 1.115683) > LambdaWindow) || (lambda.v0radius() < Lambdav0radius) || ( lambda.v0cosPA() < Lambdav0cospa) || (TMath::Abs(lambda.dcapostopv()) < Lambdadcapostopv) || (TMath::Abs(lambda.dcanegtopv()) < Lambdadcanegtopv) || ( lambda.dcaV0daughters() > Lambdadcav0dau))
      return false;
    return processCandidate(lambda);
  }

  // Process antilambda candidate 
  template <typename TV0Object>
  bool processAntiLambdaCandidate(TV0Object const& antilambda)
  {    
    // FIXME: there are smarter ways to perform this selection
    // AntiLambda base selection criteria:
    if ((std::abs(antilambda.mAntiLambda() - 1.115683) > AntiLambdaWindow) || (antilambda.v0radius() < AntiLambdav0radius) || (antilambda.v0cosPA() < AntiLambdav0cospa) || (TMath::Abs(antilambda.dcapostopv()) < AntiLambdadcapostopv) || (TMath::Abs(antilambda.dcanegtopv()) < AntiLambdadcanegtopv) || (antilambda.dcaV0daughters() > AntiLambdadcav0dau))
      return false;
    return processCandidate(antilambda);
  }

  // Process gamma candidate
  template <typename TV0Object>
  bool processGammaCandidate(TV0Object const& gamma)
  {    
    // FIXME: there are smarter ways to perform this selection
    // Gamma selection criteria:
    if ((std::abs(gamma.mGamma()) > PhotonWindow) || (gamma.mGamma() > PhotonMaxMass) || (gamma.v0radius() < PhotonMinRadius) || (gamma.v0radius() > PhotonMaxRadius) || (gamma.pt() < PhotonMinPt) || (gamma.qtarm() > PhotonMaxqt) || (TMath::Abs(gamma.alpha()) > PhotonMaxalpha))
      return false;    
    return processCandidate(gamma);
  }

  // Process k0short candidate
  template <typename TV0Object>
  bool processK0ShortCandidate(TV0Object const& kzero)
  {    
    // FIXME: there are smarter ways to perform this selection
    // TODO: Update KZeroShort selection criteria
    if ((std::abs(kzero.mK0Short() - 0.497) > K0ShortWindow) || (kzero.v0radius() < K0Shortv0radius) || (kzero.v0cosPA() < K0Shortv0cospa) || (TMath::Abs(kzero.dcapostopv()) < K0Shortdcapostopv) || (TMath::Abs(kzero.dcanegtopv()) < K0Shortdcanegtopv) || (kzero.dcaV0daughters() > K0Shortdcav0dau))
      return false;    
    return processCandidate(kzero);
  }

  void process(aod::StraCollision const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas, aod::V0TOFNSigmas> const& v0s)
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
      v0MLCandidates(Candidate.posITSCls, Candidate.negITSCls, Candidate.posITSClSize, Candidate.negITSClSize, Candidate.posTPCRows, Candidate.negTPCRows, 
      Candidate.posTPCSigmaPi, Candidate.negTPCSigmaPi, Candidate.posTPCSigmaPr, Candidate.negTPCSigmaPr, 
      Candidate.posTPCSigmaEl, Candidate.negTPCSigmaEl, Candidate.TOFSigmaLaPr, Candidate.TOFSigmaLaPi, 
      Candidate.TOFSigmaALaPi, Candidate.TOFSigmaALaPr, Candidate.TOFSigmaK0PiPlus, Candidate.TOFSigmaK0PiMinus, 
      Candidate.LambdaMass, Candidate.AntiLambdaMass, Candidate.GammaMass, Candidate.KZeroShortMass, Candidate.pT, 
      Candidate.qt, Candidate.alpha, Candidate.posEta, Candidate.negEta, Candidate.v0Eta, Candidate.Z, 
      Candidate.v0radius, Candidate.PA, Candidate.dcapostopv, Candidate.dcanegtopv, Candidate.dcaV0daughters, Candidate.dcav0topv, Candidate.PsiPair, 
      Candidate.v0type, Candidate.isLambda, Candidate.isAntiLambda, Candidate.isGamma, Candidate.isKZeroShort);
  }

}

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdakzeroMLSelectionTreeCreator>(cfgc)};
  
}