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

#ifndef PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_
#define PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_

#include <cstdlib>
#include <cmath>
#include <array>
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"

namespace o2
{
namespace pwglf
{
//__________________________________________
// V0 information storage
struct v0candidate{
  // indexing
  int collisionId = -1;
  int negativeTrack = -1;
  int positiveTrack = -1;

  // daughter properties
  std::array<float, 3> positiveMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativeMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> positivePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativePosition = {0.0f, 0.0f, 0.0f};
  float positiveTrackX = 0.0f; 
  float negativeTrackX = 0.0f;
  float positiveDCAxy = 0.0f;
  float negativeDCAxy = 0.0f;
  
  // V0 properties
  std::array<float, 3> position = {0.0f, 0.0f, 0.0f};
  float daughterDCA = 1000.0f;
  float pointingAngle = 0.0f;
  float dcaXY = 0.0f;
  
  // calculated masses for convenience
  float massGamma;
  float massK0Short;
  float massLambda;
  float massAntiLambda;

  // stored for decay chains
  float positionCovariance[6];
  float momentumCovariance[6];
};

//__________________________________________
// Cascade information storage
struct cascadeCandidate {
  // indexing
  int collisionId = -1;
  int v0Id = -1;
  int negativeTrack = -1;
  int positiveTrack = -1;
  int bachelorTrack = -1;

  // daughter properties
  std::array<float, 3> positiveMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativeMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> bachelorMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> positivePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> bachelorPosition = {0.0f, 0.0f, 0.0f};
  float positiveDCAxy = 0.0f;
  float negativeDCAxy = 0.0f;
  float bachelorDCAxy = 0.0f;
  float positiveTrackX = 0.0f; 
  float negativeTrackX = 0.0f;
  float bachelorTrackX = 0.0f;

  // cascade properties
  std::array<float, 3> position = {0.0f, 0.0f, 0.0f};
  int charge = -1; // default: []Minus
  float cascadeDaughterDCA = 1000.0f;
  float v0DaughterDCA = 1000.0f;
  
  float pointingAngle = 0.0f;
  float cascadeDCAxy = 0.0f;
  float cascadeDCAz = 0.0f;
  std::array<float, 3> v0Position = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> v0Momentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> cascadePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> cascadeMomentum = {0.0f, 0.0f, 0.0f};

  float bachBaryonCosPA = 0.0f;
  float bachBaryonDCAxyToPV = 1e+3;
  float massXi = 0.0f; 
  float massOmega = 0.0f; 

  // stored for decay chains
  float covariance[21];
};

//__________________________________________
// builder helper class
class strangenessBuilderHelper
{
  public:
    strangenessBuilderHelper(){
      // standards hardcoded in builder ...
      // ...but can be changed easily since fitter is public 
      fitter.setPropagateToPCA(true);
      fitter.setMaxR(200.);
      fitter.setMinParamChange(1e-3);
      fitter.setMinRelChi2Change(0.9);
      fitter.setMaxDZIni(1e9);
      fitter.setMaxDXYIni(4.0f);
      fitter.setMaxChi2(1e9);
      fitter.setUseAbsDCA(true);
      fitter.setWeightedFinalPCA(false);

      // LUT has to be loaded later
      lut = nullptr;
      fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrLUT);

      // mag field has to be set later
      fitter.setBz(-999.9f); // will NOT make sense if not changed
    };
    
    template <typename TTrack>
    bool buildV0Candidate(o2::aod::Collision const& collision,
                          TTrack const& positiveTrack, 
                          TTrack const& negativeTrack, 
                          bool useCollinearFit = false, 
                          bool calculateCovariance = false){
      // Calculate DCA with respect to the collision associated to the V0, not individual tracks
      gpu::gpustd::array<float, 2> dcaInfo;

      auto posTrackPar = getTrackPar(positiveTrack);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, posTrackPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
      v0.positiveDCAxy = dcaInfo[0];

      auto negTrackPar = getTrackPar(negativeTrack);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, negTrackPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
      v0.negativeDCAxy = dcaInfo[0];

      o2::track::TrackParCov positiveTrackParam = getTrackParCov(positiveTrack);
      o2::track::TrackParCov negativeTrackParam = getTrackParCov(negativeTrack);

      // Perform DCA fit
      int nCand = 0;
      fitter.setCollinear(useCollinearFit);
      try {
        nCand = fitter.process(positiveTrackParam, negativeTrackParam);
      } catch (...) {
        return false;
      }
      if (nCand == 0) {
        return false;
      }

      v0.positiveTrackX = fitter.getTrack(0).getX();
      v0.negativeTrackX = fitter.getTrack(1).getX();
      positiveTrackParam = fitter.getTrack(0);
      negativeTrackParam = fitter.getTrack(1);
      positiveTrackParam.getPxPyPzGlo(v0.positiveMomentum);
      negativeTrackParam.getPxPyPzGlo(v0.negativeMomentum);
      positiveTrackParam.getXYZGlo(v0.positivePosition);
      negativeTrackParam.getXYZGlo(v0.negativePosition);

      // get decay vertex coordinates
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        v0.position[i] = vtx[i];
      }

      v0.daughterDCA = TMath::Sqrt(fitter.getChi2AtPCACandidate());
      v0.pointingAngle = TMath::ACos(RecoDecay::cpa(
        std::array{collision.posX(), collision.posY(), collision.posZ()}, 
        std::array{v0.position[0], v0.position[1], v0.position[2]}, 
        std::array{v0.positiveMomentum[0] + v0.negativeMomentum[0], v0.positiveMomentum[1] + v0.negativeMomentum[1], v0.positiveMomentum[2] + v0.negativeMomentum[2]}
      ));

      v0.dcaXY = CalculateDCAStraightToPV(
        v0.position[0], v0.position[1], v0.position[2],
        v0.positiveMomentum[0] + v0.negativeMomentum[0],
        v0.positiveMomentum[1] + v0.negativeMomentum[1],
        v0.positiveMomentum[2] + v0.negativeMomentum[2],
        collision.posX(), collision.posY(), collision.posZ());

      // Calculate masses
      v0.massGamma = RecoDecay::m(std::array{
        std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]}, 
        std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}}, 
        std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
      v0.massK0Short = RecoDecay::m(std::array{
        std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]}, 
        std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}}, 
        std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
      v0.massLambda = RecoDecay::m(std::array{
        std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]}, 
        std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}}, 
        std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
      v0.massAntiLambda = RecoDecay::m(std::array{
        std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]}, 
        std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}}, 
        std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});

      // calculate covariance if requested
      if (calculateCovariance) {
        // Calculate position covariance matrix
        auto covVtxV = fitter.calcPCACovMatrix(0);
        // std::array<float, 6> positionCovariance;
        v0.positionCovariance[0] = covVtxV(0, 0);
        v0.positionCovariance[1] = covVtxV(1, 0);
        v0.positionCovariance[2] = covVtxV(1, 1);
        v0.positionCovariance[3] = covVtxV(2, 0);
        v0.positionCovariance[4] = covVtxV(2, 1);
        v0.positionCovariance[5] = covVtxV(2, 2);
        std::array<float, 21> covTpositive = {0.};
        std::array<float, 21> covTnegative = {0.};
        // std::array<float, 6> momentumCovariance;
        float momentumCovariance[6];
        positiveTrackParam.getCovXYZPxPyPzGlo(covTpositive);
        negativeTrackParam.getCovXYZPxPyPzGlo(covTnegative);
        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 6; i++) {
          v0.momentumCovariance[i] = covTpositive[MomInd[i]] + covTnegative[MomInd[i]];
        }
      }

      // information validated, V0 built successfully. Signal OK
      return true;   
    }

    // cascade builder creating a cascade from plain tracks
    template <typename TTrack>
    bool buildCascadeCandidate(o2::aod::Collision const& collision,
                              TTrack const& positiveTrack, 
                              TTrack const& negativeTrack, 
                              TTrack const& bachelorTrack, 
                              bool calculateBachelorBaryonVariables = false, 
                              bool useCascadeMomentumAtPV = false, 
                              bool processCovariances = false)
    {
      if(!buildV0Candidate(collision, positiveTrack, negativeTrack, false, processCovariances)){ 
        return false;
      }
      if(!buildCascadeCandidate(collision, v0, positiveTrack, negativeTrack, bachelorTrack, calculateBachelorBaryonVariables, useCascadeMomentumAtPV, processCovariances)){
        return false;
      }
      return true;
    }

    // cascade builder using pre-fabricated information, thus not calling
    // the DCAfitter again for the V0 contained in the cascade
    // if generating from scratch, prefer the other variant
    template <typename TTrack>
    bool buildCascadeCandidate(o2::aod::Collision const& collision,
                              v0candidate const& v0input, 
                              TTrack const& positiveTrack, 
                              TTrack const& negativeTrack, 
                              TTrack const& bachelorTrack, 
                              bool calculateBachelorBaryonVariables = false, 
                              bool useCascadeMomentumAtPV = false,
                              bool processCovariances = false)
    {
      if (calculateBachelorBaryonVariables) {
        // Calculates properties of the V0 comprised of bachelor and baryon in the cascade
        // baryon: distinguished via bachelor charge
        if (bachelorTrack.sign() < 0) {
          processBachBaryonVariables(collision, bachelorTrack, positiveTrack);
        } else {
          processBachBaryonVariables(collision, bachelorTrack, negativeTrack);
        }
      }

      // Overall cascade charge
      cascade.charge = bachelorTrack.signed1Pt() > 0 ? +1 : -1;

      // bachelor DCA track to PV
      // Calculate DCA with respect to the collision associated to the V0, not individual tracks
      gpu::gpustd::array<float, 2> dcaInfo;

      auto bachTrackPar = getTrackPar(bachelorTrack);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, bachTrackPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
      cascade.bachelorDCAxy = dcaInfo[0];

      // Do actual minimization
      auto lBachelorTrack = getTrackParCov(bachelorTrack);

      // Set up covariance matrices (should in fact be optional)
      std::array<float, 21> covV = {0.};
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covV[MomInd[i]] = v0input.momentumCovariance[i];
        covV[i] = v0input.positionCovariance[i];
      }
      auto lV0Track = o2::track::TrackParCov(
        {v0input.position[0], v0input.position[1], v0input.position[2]},
        {v0input.positiveMomentum[0] + v0input.negativeMomentum[0], v0input.positiveMomentum[1] + v0input.negativeMomentum[1], v0input.positiveMomentum[2] + v0input.negativeMomentum[2]},
        covV, 0, true);
      lV0Track.setAbsCharge(0);
      lV0Track.setPID(o2::track::PID::Lambda);

      //---/---/---/
      // Move close to minima
      int nCand = 0;
      try {
        nCand = fitter.process(lV0Track, lBachelorTrack);
      } catch (...) {
        return false;
      }
      if (nCand == 0)
        return false;

      lV0Track = fitter.getTrack(0);
      lBachelorTrack = fitter.getTrack(1);

      // DCA between cascade daughters
      cascade.cascadeDaughterDCA = TMath::Sqrt(fitter.getChi2AtPCACandidate());

      lBachelorTrack.getPxPyPzGlo(cascade.bachelorMomentum);
      // get decay vertex coordinates
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        cascade.position[i] = vtx[i];
      }

      cascade.pointingAngle = TMath::ACos(RecoDecay::cpa(
        std::array{collision.posX(), collision.posY(), collision.posZ()},
        std::array{cascade.position[0], cascade.position[1], cascade.position[2]},
        std::array{v0input.positiveMomentum[0] + v0input.negativeMomentum[0] + cascade.bachelorMomentum[0], v0input.positiveMomentum[0] + v0input.negativeMomentum[1] + cascade.bachelorMomentum[1], v0input.positiveMomentum[2] + v0input.negativeMomentum[2] + cascade.bachelorMomentum[2]}));

      // Calculate DCAxy of the cascade (with bending)
      auto lCascadeTrack = fitter.createParentTrackParCov();
      lCascadeTrack.setAbsCharge(cascade.charge); // to be sure
      lCascadeTrack.setPID(o2::track::PID::XiMinus);       // FIXME: not OK for omegas
      dcaInfo[0] = 999;
      dcaInfo[1] = 999;

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lCascadeTrack, 2.f, fitter.getMatCorrType(), &dcaInfo);
      cascade.cascadeDCAxy = dcaInfo[0];
      cascade.cascadeDCAz = dcaInfo[1];

      // Calculate masses a priori
      cascade.massXi = RecoDecay::m(std::array{std::array{cascade.bachelorMomentum[0], cascade.bachelorMomentum[1], cascade.bachelorMomentum[2]}, 
                                              std::array{v0input.positiveMomentum[0] + v0input.negativeMomentum[0], v0input.positiveMomentum[1] + v0input.negativeMomentum[1], v0input.positiveMomentum[2] + v0input.negativeMomentum[2]}}, 
                                              std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
      cascade.massOmega = RecoDecay::m(std::array{std::array{cascade.bachelorMomentum[0], cascade.bachelorMomentum[1], cascade.bachelorMomentum[2]}, 
                                                  std::array{v0input.positiveMomentum[0] + v0input.negativeMomentum[0], v0input.positiveMomentum[1] + v0input.negativeMomentum[1], v0input.positiveMomentum[2] + v0input.negativeMomentum[2]}}, 
                                                  std::array{o2::constants::physics::MassKaonCharged, o2::constants::physics::MassLambda});

      // Populate information
      // cascadecandidate.v0Id = v0index.globalIndex();
      cascade.positiveTrack = positiveTrack.globalIndex();
      cascade.negativeTrack = negativeTrack.globalIndex();
      cascade.bachelorTrack = bachelorTrack.globalIndex();
      cascade.positiveTrackX = v0input.positiveTrackX;
      cascade.negativeTrackX = v0input.negativeTrackX;
      cascade.bachelorTrackX = lBachelorTrack.getX(); // from this minimization
      cascade.v0Position[0] = v0input.position[0];
      cascade.v0Position[1] = v0input.position[1];
      cascade.v0Position[2] = v0input.position[2];
      cascade.v0Momentum[0] = v0input.positiveMomentum[0] + v0input.negativeMomentum[0];
      cascade.v0Momentum[1] = v0input.positiveMomentum[1] + v0input.negativeMomentum[1];
      cascade.v0Momentum[2] = v0input.positiveMomentum[2] + v0input.negativeMomentum[2];
      cascade.positiveMomentum[0] = v0input.positiveMomentum[0];
      cascade.positiveMomentum[1] = v0input.positiveMomentum[1];
      cascade.positiveMomentum[2] = v0input.positiveMomentum[2];
      cascade.negativeMomentum[0] = v0input.negativeMomentum[0];
      cascade.negativeMomentum[1] = v0input.negativeMomentum[1];
      cascade.negativeMomentum[2] = v0input.negativeMomentum[2];
      cascade.v0DaughterDCA = v0input.daughterDCA;
      cascade.positiveDCAxy = v0input.positiveDCAxy;
      cascade.negativeDCAxy = v0input.negativeDCAxy;

      if (useCascadeMomentumAtPV) {
        lCascadeTrack.getPxPyPzGlo(cascade.cascadeMomentum);
      } else {
        cascade.cascadeMomentum[0] = cascade.bachelorMomentum[0] + v0input.positiveMomentum[0] + v0input.negativeMomentum[0];
        cascade.cascadeMomentum[1] = cascade.bachelorMomentum[1] + v0input.positiveMomentum[1] + v0input.negativeMomentum[1];
        cascade.cascadeMomentum[2] = cascade.bachelorMomentum[2] + v0input.positiveMomentum[2] + v0input.negativeMomentum[2];
      }

      if(processCovariances){ 
        // Calculate position covariance matrix
        auto covVtxV = fitter.calcPCACovMatrix(0);
        // std::array<float, 6> positionCovariance;
        float positionCovariance[6];
        positionCovariance[0] = covVtxV(0, 0);
        positionCovariance[1] = covVtxV(1, 0);
        positionCovariance[2] = covVtxV(1, 1);
        positionCovariance[3] = covVtxV(2, 0);
        positionCovariance[4] = covVtxV(2, 1);
        positionCovariance[5] = covVtxV(2, 2);
        // store momentum covariance matrix
        std::array<float, 21> covTv0 = {0.};
        std::array<float, 21> covTbachelor = {0.};
        float covCascade[21];
        // std::array<float, 6> momentumCovariance;
        lV0Track.getCovXYZPxPyPzGlo(covTv0);
        lBachelorTrack.getCovXYZPxPyPzGlo(covTbachelor);
        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 21; i++) {
          cascade.covariance[i] = 0.0f;
        }
        for (int i = 0; i < 6; i++) {
          cascade.covariance[i] = positionCovariance[i];
          cascade.covariance[MomInd[i]] = covTv0[MomInd[i]] + covTbachelor[MomInd[i]];
        }
      }

      // Final outcome is YES if I got here!
      return true;
    }

    o2::base::MatLayerCylSet* lut; // material LUT for DCA fitter
    o2::vertexing::DCAFitterN<2> fitter; // 2-prong o2 dca fitter

    v0candidate v0;  // storage for V0 candidate properties
    cascadeCandidate cascade; // storage for cascade candidate properties

  private: 
    // internal helper to calculate DCAxy of a straight line to a given PV analytically
    float CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ){
      return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
    }

    template <typename TCollision, typename TTrack>
    void processBachBaryonVariables(TCollision const& collision, TTrack const& track1, TTrack const& track2)
    {
      cascade.bachBaryonCosPA = 0;        // would ordinarily accept all
      cascade.bachBaryonDCAxyToPV = 1e+3; // would ordinarily accept all

      // create tracks from table rows
      o2::track::TrackParCov tr1 = getTrackParCov(track1);
      o2::track::TrackParCov tr2 = getTrackParCov(track2);

      //---/---/---/
      // Move close to minima
      int nCand = 0;
      try {
        nCand = fitter.process(tr1, tr2);
      } catch (...) {
        return;
      }
      if (nCand == 0)
        return; // variables are such that candidate is accepted (not obvious...)

      // Calculate DCAxy of the cascade (with bending)
      o2::track::TrackPar wrongV0 = fitter.createParentTrackPar();
      wrongV0.setAbsCharge(0); // charge zero
      gpu::gpustd::array<float, 2> dcaInfo;
      dcaInfo[0] = 999;
      dcaInfo[1] = 999;

      // bachelor-baryon DCAxy to PV
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, wrongV0, 2.f, fitter.getMatCorrType(), &dcaInfo);
      cascade.bachBaryonDCAxyToPV = dcaInfo[0];

      const auto& vtx = fitter.getPCACandidate();
      if (!fitter.isPropagateTracksToVertexDone())
        return;

      std::array<float, 3> tr1p;
      std::array<float, 3> tr2p;

      fitter.getTrack(1).getPxPyPzGlo(tr1p);
      fitter.getTrack(2).getPxPyPzGlo(tr2p);

      // bachelor-baryon CosPA
      cascade.bachBaryonCosPA = RecoDecay::cpa(
        std::array{collision.posX(), collision.posY(), collision.posZ()},
        std::array{vtx[0], vtx[1], vtx[2]},
        std::array{tr1p[0] + tr2p[0], tr1p[1] + tr2p[1], tr1p[2] + tr2p[2]});

      // Potentially also to be considered: bachelor-baryon DCA (between the two tracks)
      // to be added here as complementary information in the future
    }
};

} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_
