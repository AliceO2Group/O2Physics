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
// Task to add a table of track parameters propagated to the primary vertex
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "TRandom3.h"

// This task goes straight from a combination of track table and mcParticles 
// and a custom TOF configuration to a table of TOF NSigmas for the particles
// being analysed. It currently contemplates 5 particle types: 
// electrons, pions, kaons, protons and muons
//
// More particles could be added but would have to be added to the LUT 
// being used in the onTheFly tracker task. 

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace upgradetof
{
DECLARE_SOA_COLUMN(NSigmaElectronInner, nSigmaElectronInner, float); //! NSigma electron iTOF
DECLARE_SOA_COLUMN(NSigmaMuonInner, nSigmaMuonInner, float); //! NSigma muon iTOF
DECLARE_SOA_COLUMN(NSigmaPionInner, nSigmaPionInner, float); //! NSigma pion iTOF
DECLARE_SOA_COLUMN(NSigmaKaonInner, nSigmaKaonInner, float); //! NSigma kaon iTOF
DECLARE_SOA_COLUMN(NSigmaProtonInner, nSigmaProtonInner, float); //! NSigma proton iTOF
DECLARE_SOA_COLUMN(NSigmaElectronOuter, nSigmaElectronOuter, float); //! NSigma electron oTOF
DECLARE_SOA_COLUMN(NSigmaMuonOuter, nSigmaMuonOuter, float); //! NSigma muon oTOF
DECLARE_SOA_COLUMN(NSigmaPionOuter, nSigmaPionOuter, float); //! NSigma pion oTOF
DECLARE_SOA_COLUMN(NSigmaKaonOuter, nSigmaKaonOuter, float); //! NSigma kaon oTOF
DECLARE_SOA_COLUMN(NSigmaProtonOuter, nSigmaProtonOuter, float); //! NSigma proton oTOF
}
DECLARE_SOA_TABLE(UpgradeTof, "AOD", "UPGRADETOF",
                  upgradetof::NSigmaElectronInner,
                  upgradetof::NSigmaMuonInner,
                  upgradetof::NSigmaPionInner,
                  upgradetof::NSigmaKaonInner,
                  upgradetof::NSigmaProtonInner,
                  upgradetof::NSigmaElectronOuter,
                  upgradetof::NSigmaMuonOuter,
                  upgradetof::NSigmaPionOuter,
                  upgradetof::NSigmaKaonOuter,
                  upgradetof::NSigmaProtonOuter
                  );
} // namespace o2::aod

struct onTheFlyTOFPID {
  Produces<aod::UpgradeTof> upgradetof;

  // necessary for particle charges 
  Service<O2DatabasePDG> pdg;

  // these are the settings governing the TOF layers to be used
  // note that there are two layers foreseen for now: inner and outer TOF 
  // more could be added (especially a disk TOF at a certain z?)
  // in the evolution of this effort
  Configurable<float> dBz{"dBz", 20, "magnetic field (kilogauss)"};
  Configurable<float> iTOFRadius{"iTOFRadius", 20, "barrel iTOF radius (cm)"};
  Configurable<float> oTOFRadius{"oTOFRadius", 80, "barrel oTOF radius (cm)"};
  Configurable<float> iTOFTimeReso{"iTOFTimeReso", 20, "barrel iTOF time error (ps)"};
  Configurable<float> oTOFTimeReso{"oTOFTimeReso", 20, "barrel oTOF time error (ps)"};
  Configurable<int> lNStepsLIntegrator{"lNStepsLIntegrator", 200, "number of steps in length integrator"};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // needed: random number generator for smearing
  TRandom3 lPRNG;
  
  void init(o2::framework::InitContext& initContext) {
    lPRNG.SetSeed(0); // fully randomize
  }

  template <typename mcParticleType>
  void convertMCParticleToO2Track(mcParticleType& particle, o2::track::TrackParCov& o2track)
  {
    // FIXME: this is a fundamentally important piece of code. 
    // It could be placed in a utility file instead of here. 
    std::array<float, 3> xyz = {static_cast<float>(particle.vx()), static_cast<float>(particle.vy()), static_cast<float>(particle.vz())};
    std::array<float, 3> ptetaphi = {static_cast<float>(particle.pt()), static_cast<float>(particle.eta()), static_cast<float>(particle.phi())};
    auto pdgInfo = pdg->GetParticle(particle.pdgCode());
    int charge = 0; 
    if (pdgInfo != nullptr)
        charge = pdgInfo->Charge();
    std::array<float, 5> params;
    std::array<float, 15> covm = {
          0.,
          0., 0.,
          0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0., 0.};
    float s, c, x;
    o2::math_utils::sincos(ptetaphi[2], s, c);
    o2::math_utils::rotateZInv(xyz[0], xyz[1], x, params[0], s, c);
    params[1] = xyz[2];
    params[2] = 0.; // since alpha = phi
    auto theta = 2. * atan(exp(-ptetaphi[1]));  
    params[3] = 1. / tan(theta);
    params[4] = charge / ptetaphi[0];
    
    // Initialize TrackParCov in-place
    new (&o2track) (o2::track::TrackParCov)(x, ptetaphi[2], params, covm);
  }

  Double_t TrackLength( o2::track::TrackParCov track, Double_t lX0, Double_t lX1 , Double_t lMagneticField ){
    std::array<float, 3> lPointN;
    std::array<float, 3> lPointNplus;
    Double_t lLength = 0.0;
    track.propagateTo(lX0, lMagneticField);
    for(Int_t iStep=1; iStep<lNStepsLIntegrator; iStep++){
      track.getXYZGlo(lPointN);
      Float_t lPosition = lX0 + (lX1-lX0)*((Float_t)(iStep))/((Float_t)(lNStepsLIntegrator-1));
      track.propagateTo(lPosition, lMagneticField);
      track.getXYZGlo(lPointNplus);
      lLength += std::hypot( lPointNplus[0]-lPointN[0], lPointNplus[1]-lPointN[1], lPointNplus[2]-lPointN[2] );
    }
    return lLength;
  }

  Double_t Velocity(Double_t lMomentum, Double_t lMass){
    //Momentum p and mass m -> returns speed in centimeters per picosecond
    //Useful for TOF calculations
    Double_t lA = TMath::Power(lMomentum / lMass, 2);
    return 0.0299792458*TMath::Sqrt(lA/(1+lA));
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks, aod::McParticles const&)
  {
    o2::math_utils::Point3D<float> pvpos{collision.posX(),collision.posY(),collision.posZ()};
    std::array<float, 6> pvcov;
    pvcov[0] = collision.covXX();
    pvcov[1] = collision.covXY();
    pvcov[2] = collision.covYY();
    pvcov[3] = collision.covXZ();
    pvcov[4] = collision.covXY();
    pvcov[5] = collision.covXZ();
    o2::dataformats::VertexBase pvvtx(pvpos, pvcov);

    auto mcCollision = collision.mcCollision();
    o2::math_utils::Point3D<float> mcpvpos{mcCollision.posX(),mcCollision.posY(),mcCollision.posZ()};
    std::array<float, 6> mcpvcov; // dummy! 
    for(int ii=0; ii<6; ii++) mcpvcov[ii] = 1e-6;
    o2::dataformats::VertexBase mcpvvtx(mcpvpos, mcpvcov);

    for (auto& track : tracks) {
      // first step: find precise arrival time (if any)
      // --- convert track into perfect track 
      if ( !track.has_mcParticle() ) // should always be OK but check please
        LOG(error) << "Oh no! No mcParticle label for this track! This shouldn't happen!";

      o2::track::TrackParCov o2track; 
      auto mcParticle = track.mcParticle(); 
      convertMCParticleToO2Track( mcParticle, o2track );
  
      float lX_PV=-100, lX_iTOF=-100, lX_oTOF=-100, lThisTrackLength_iTOF=-1, lThisTrackLength_oTOF=-1;
      if (o2track.propagateToDCA(mcpvvtx, dBz))
        lX_PV = o2track.getX();
      if (!o2track.getXatLabR(iTOFRadius,lX_iTOF,dBz,o2::track::DirOutward)) 
        lX_iTOF = -100;
      if (!o2track.getXatLabR(oTOFRadius,lX_oTOF,dBz,o2::track::DirOutward)) 
        lX_oTOF = -100;
      if(lX_PV>-99.&&lX_iTOF>-99.) lThisTrackLength_iTOF = TrackLength(o2track, lX_PV, lX_iTOF, dBz);
      if(lX_PV>-99.&&lX_oTOF>-99.) lThisTrackLength_oTOF = TrackLength(o2track, lX_PV, lX_oTOF, dBz);

      // get mass to calculate velocity 
      auto pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      float lExpectedTime_iTOF = lThisTrackLength_iTOF / Velocity(o2track.getP(), pdgInfo->Mass() );
      float lExpectedTime_oTOF = lThisTrackLength_oTOF / Velocity(o2track.getP(), pdgInfo->Mass() );

      //Smear with expected resolutions
      float lMeasuredTime_iTOF = lPRNG.Gaus(lExpectedTime_iTOF, iTOFTimeReso);
      float lMeasuredTime_oTOF = lPRNG.Gaus(lExpectedTime_oTOF, iTOFTimeReso);

      // Now we calculate the expected arrival time following certain mass hypotheses
      // and the (imperfect!) reconstructed track parametrizations
      float lThisTrackLengthReco_iTOF=-1, lThisTrackLengthReco_oTOF=-1;
      auto recoTrack = getTrackParCov(track);
      if (recoTrack.propagateToDCA(pvvtx, dBz))
        lX_PV = recoTrack.getX();
      if (!recoTrack.getXatLabR(iTOFRadius,lX_iTOF,dBz,o2::track::DirOutward)) 
        lX_iTOF = -100;
      if (!recoTrack.getXatLabR(oTOFRadius,lX_oTOF,dBz,o2::track::DirOutward)) 
        lX_oTOF = -100;
      if(lX_PV>-99.&&lX_iTOF>-99.) lThisTrackLengthReco_iTOF = TrackLength(recoTrack, lX_PV, lX_iTOF, dBz);
      if(lX_PV>-99.&&lX_oTOF>-99.) lThisTrackLengthReco_oTOF = TrackLength(recoTrack, lX_PV, lX_oTOF, dBz);

      // Straight to Nsigma
      float lDeltaTime_iTOF[5], lNSigma_iTOF[5]; 
      float lDeltaTime_oTOF[5], lNSigma_oTOF[5];
      int lpdg_array[5] = {11, 13, 211, 321, 2212}; 
      float lMasses[5]; 

      for(int ii=0; ii<5; ii++){ 
        lNSigma_iTOF[ii] = -100;
        lNSigma_oTOF[ii] = -100;

        auto pdgInfoThis = pdg->GetParticle(lpdg_array[ii]);
        lMasses[ii] = pdgInfoThis->Mass();
        lDeltaTime_iTOF[ii] = lThisTrackLengthReco_iTOF / Velocity(recoTrack.getP(), lMasses[ii]) - lMeasuredTime_iTOF;
        lDeltaTime_oTOF[ii] = lThisTrackLengthReco_oTOF / Velocity(recoTrack.getP(), lMasses[ii]) - lMeasuredTime_oTOF;

        //Fixme: assumes dominant resolution effect is the TOF resolution 
        //and not the tracking itself. It's *probably* a fair assumption 
        //but it should be tested further! 
        if( lThisTrackLength_iTOF > 0 && lThisTrackLengthReco_iTOF)
          lNSigma_iTOF[ii] = lDeltaTime_iTOF[ii]/iTOFTimeReso;
        if( lThisTrackLength_oTOF > 0 && lThisTrackLengthReco_oTOF)
          lNSigma_oTOF[ii] = lDeltaTime_oTOF[ii]/oTOFTimeReso;
      }

      // Sigmas have been fully calculated. Please populate the NSigma helper table (once per track)
      upgradetof(lNSigma_iTOF[0], lNSigma_iTOF[1], lNSigma_iTOF[2], lNSigma_iTOF[3], lNSigma_iTOF[4], 
                 lNSigma_oTOF[0], lNSigma_oTOF[1], lNSigma_oTOF[2], lNSigma_oTOF[3], lNSigma_oTOF[4]);

    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<onTheFlyTOFPID>(cfgc)};
  return workflow;
}
