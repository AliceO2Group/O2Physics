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

#ifndef PWGLF_UTILS_DECAY3BODYBUILDERHELPER_H_
#define PWGLF_UTILS_DECAY3BODYBUILDERHELPER_H_

#include <cstdlib>
#include <cmath>
#include <array>
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "Tools/KFparticle/KFUtilities.h"

#ifndef HomogeneousField
#define HomogeneousField
#endif

/// includes KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

namespace o2
{
namespace pwglf
{

//_______________________________________________________________________
// Deca3body information storage
struct decay3bodyCandidate {
  // indexing
  int collisionID = -1;
  int decay3bodyID = -1;
  int trackPosID = -1;
  int trackNegID = -1;
  int trackBachID = -1;

  // daughter properties
  std::array<float, 3> trackPosMom = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> trackNegMom = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> trackBachMom = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> trackPosPos = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> trackNegPos = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> trackBachPos = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> trackDCAxyToPV = {0.0f, 0.0f, 0.0f};   // 0 - pos, 1 - neg, 2 - bach
  std::array<float, 3> trackDCAzToPV = {0.0f, 0.0f, 0.0f};    // 0 - pos, 1 - neg, 2 - bach
  std::array<float, 4> tpcNsigma = {0.0f, 0.0f, 0.0f, 0.0f};  // 0 - proton, 1 - pion, 2 - deuteron, 3 - bach with pion hyp
  float tofNsigmaDeuteron = 0.0f;
  float averageITSClusSizeDeuteron = 0.0f;
  int pidForTrackingDeuteron = 0;

  // vertex properties
  int charge = -1;
  float momentum[3];
  float position[3];
  // std::array<float, 3> momentum = {0.0f, 0.0f, 0.0f};
  // std::array<float, 3> position = {0.0f, 0.0f, 0.0f};
  // float dcaToPV = 0.0f;
  // float dcaxyToPV = 0.0f;
  float chi2 = 0.0f;
  float trackedClSize = 0.0f;
  float daughterDCAatSV = 0.0f;                               // quadratic sum of DCA between daughters at SV
  std::array<float, 3> daughterDCAtoSV = {0.0f, 0.0f, 0.0f};  // 0 - pos, 1 - neg, 2 - bach

  // covariance matrix
  float trackPosCov[21] = {0.0f};
  float trackNegCov[21] = {0.0f};
  float trackBachCov[21] = {0.0f};
  float covariance[21] = {0.0f};
};

//_______________________________________________________________________
// builder helper class
class Decay3BodyBuilderHelper
{
 public:
  strangenessBuilderHelper()
    {
    // standards hardcoded in builder ...
    // ...but can be changed easily since fitter is public
    fitter3body.setPropagateToPCA(true);
    fitter3body.setMaxR(200.);
    fitter3body.setMinParamChange(1e-3);
    fitter3body.setMinRelChi2Change(0.9);
    fitter3body.setMaxDZIni(1e9);
    fitter3body.setMaxDXYIni(4.0f);
    fitter3body.setMaxChi2(1e9);
    fitter3body.setUseAbsDCA(true);
    fitter3body.setWeightedFinalPCA(false);

    // LUT has to be loaded later
    lut = nullptr;
    fitter3body.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrLUT);

    // mag field has to be set later
    fitter3body.setBz(-999.9f); // will NOT make sense if not changed
  };

  o2::base::MatLayerCylSet* lut;            // material LUT for DCA fitter
  o2::vertexing::DCAFitterN<2> fitterV0;    // 2-prong o2 dca fitter
  o2::vertexing::DCAFitterN<3> fitter3body; // 3-prong o2 dca fitter

  decay3bodyCandidate decay3body; // storage for Decay3body candidate properties

  // decay3body candidate criteria
  struct {
    // daughter tracks
    float maxEtaDaughters;
    int minTPCNClProton;
    int minTPCNClPion;
    int minTPCNClBach;
    float minDCAProtonToPV;
    float minDCAPionToPV;
    float minDCABachToPV;
    float minPtProton;
    float minPtPion;
    float minPtBach;
    float maxPtProton;
    float maxPtPion;
    float maxPtBach;
    float maxTPCnSigma;
    float minTOFnSigmaDeuteron;
    float maxTOFnSigmaDeuteron;
    float minPBachUseTOF;
    float maxDCADauAtSV;
    // candidate
    float maxRapidity;
    float minPt;
    float maxPt;
    float minMass;
    float maxMass;
    float minCtau;
    float maxCtau;
    float minCosPA;
    float maxChi2;
  } decay3bodyselections;

  //_______________________________________________________________________
  // build Decay3body from three tracks, including V0 building.
  template <typename TCollision, typename TTrack>
  bool buildDecay3BodyCandidate(TCollision const& collision,
                                TTrack const& trackPos, 
                                TTrack const& trackNeg,
                                TTrack const& trackBach,
                                int decay3bodyIndex, 
                                float todNsigmaDeuteron, 
                                float trackedClSize,
                                int bachelorcharge = 1,
                                bool useKFParticle = false,
                                bool kfSetTopologicalConstraint = false,
                                bool useSelections = true,
                                bool useTPCforPion = false,
                                bool acceptTPCOnly = false,
                                bool calculateCovariance = true)
  {
    int collisionIndex = collision.globalIndex();
    float pvX = collision.posX();
    float pvY = collision.posY();
    float pvZ = collision.posZ();
    bool isMatter = trackBach.sign() > 0 ? true : false;

    auto trackParCovPos = getTrackParCov(trackPos);
    auto trackParCovNeg = getTrackParCov(trackNeg);
    auto trackParCovBach = getTrackParCov(trackBach);

    decay3body.collisionID = collisionIndex;
    decay3body.decay3bodyID = decay3bodyIndex;
    decay3body.trackPosID = trackPos.globalIndex();
    decay3body.trackNegID = trackNeg.globalIndex();
    decay3body.trackBachID = trackBach.globalIndex();

    //_______________________________________________________________________
    // track selections
    if constexpr (useSelections) {
      // proton track quality
      if (isMatter && trackPos.tpcNCls() < decay3bodyselections.minTPCNClProton) {
        decay3body = {};
        return false;
      } else if (!isMatter && trackNeg.tpcNCls() < decay3bodyselections.minTPCNClProton) {
        decay3body = {};
        return false;
      }
      // pion track quality
      if (useTPCforPion) {
        if (isMatter && trackNeg.tpcNCls() < decay3bodyselections.minTPCNClPion) {
          decay3body = {};
          return false;
        } else if (!isMatter && trackPos.tpcNCls() < decay3bodyselections.minTPCNClPion) {
          decay3body = {};
          return false;
        }
      }
      // bachelor track quality
      if (trackBach.tpcNCls() < decay3bodyselections.minTPCNClBach) {
        decay3body = {};
        return false;
      }

      // track signs
      if (trackPos.sign() != +1 || trackNeg.sign() != -1) {
        decay3body = {};
        return false;
      }

      // track eta
      if (std::fabs(trackPos.eta()) > decay3bodyselections.maxEtaDaughters) {
        decay3body = {};
        return false;
      }
      if (std::fabs(trackNeg.eta()) > decay3bodyselections.maxEtaDaughters) {
        decay3body = {};
        return false;
      }
      if (std::fabs(trackBach.eta()) > decay3bodyselections.maxEtaDaughters) {
        decay3body = {};
        return false;
      }

      // TPC only
      if (!acceptTPCOnly && !trackPos.hasITS()) {
        decay3body = {};
        return false;
      }
      if (!acceptTPCOnly && !trackNeg.hasITS()) {
        decay3body = {};
        return false;
      }
      if (!acceptTPCOnly && !trackBach.hasITS()) {
        decay3body = {};
        return false;
      }

      // daughter TPC PID
      if (isMatter) {
        if (std::fabs(trackPos.tpcNsigmaPr()) > decay3bodyselections.maxTPCnSigma) {
          decay3body = {};
          return false;
        } else if(useTPCforPion && std::fabs(trackNeg.tpcNsigmaPi()) > decay3bodyselections.maxTPCnSigma) {
          decay3body = {};
          return false;
        }
      } else if (!isMatter) {
        if (std::fabs(trackNeg.tpcNsigmaPr()) > decay3bodyselections.maxTPCnSigma) {
          decay3body = {};
          return false;
        } else if (useTPCforPion && std::fabs(trackPos.tpcNsigmaPi()) > decay3bodyselections.maxTPCnSigma) {
          decay3body = {};
          return false;
        }
      } 
      if (std::fabs(trackBach.tpcNsigmaDeuteron()) > decay3bodyselections.maxTPCnSigma) {
        decay3body = {};
        return false;
      }

      // deuteron TOF PID
      if ((tofNsigmaDeuteron < decay3bodyselections.minTOFnSigmaDeuteron || tofNsigmaDeuteron > decay3bodyselections.maxTOFnSigmaDeuteron) && trackBach.p() > decay3bodyselections.minPBachUseTOF) {
        decay3body = {};
        return false;
      }
    } // end of selections

    //_______________________________________________________________________
    // daughter track DCA to PV associated with decay3body
    o2::dataformats::VertexBase mPV;
    o2::dataformats::DCA mDcaInfoCov;
    auto trackParCovPosCopy = trackParCovPos;
    auto trackParCovNegCopy = trackParCovNeg;
    auto trackParCovBachCopy = trackParCovBach;
    mPV.setPos({pvX, pvY, pvZ});
    mPV.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

    // positive track
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovPosCopy, 2.f, matCorr, &mDcaInfoCov);
    decay3body.trackDCAxyToPV[0] = mDcaInfoCov.getY();
    decay3body.trackDCAzToPV[0] = mDcaInfoCov.getZ();
    auto trackPosDCAToPV = std::sqrt(decay3body.trackDCAxyToPV[0] * decay3body.trackDCAxyToPV[0] + decay3body.trackDCAzToPV[0] * decay3body.trackDCAzToPV[0]);
    if constexpr (useSelections) {
      if (isMatter && trackPosDCAToPV < decay3bodyselections.minDCAProtonToPV) {
        decay3body = {};
        return false;
      } else if (!isMatter && trackPosDCAToPV < decay3bodyselections.minDCAPionToPV) {
        decay3body = {};
        return false;
      }
    }
    // negative track
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParNegPosCopy, 2.f, matCorr, &mDcaInfoCov);
    decay3body.trackDCAxyToPV[1] = mDcaInfoCov.getY();
    decay3body.trackDCAzToPV[1] = mDcaInfoCov.getZ();
    auto trackNegDCAToPV = std::sqrt(decay3body.trackDCAxyToPV[1] * decay3body.trackDCAxyToPV[1] + decay3body.trackDCAzToPV[1] * decay3body.trackDCAzToPV[1]);
    if constexpr (useSelections) {
      if (isMatter && trackNegDCAToPV < decay3bodyselections.minDCAPionToPV) {
        decay3body = {};
        return false;
      } else if (!isMatter && trackNegDCAToPV < decay3bodyselections.minDCAProtonToPV) {
        decay3body = {};
        return false;
      }
    }
    // bachelor track
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovBachCopy, 2.f, matCorr, &mDcaInfoCov);
    decay3body.trackDCAxyToPV[2] = mDcaInfoCov.getY();
    decay3body.trackDCAzToPV[2] = mDcaInfoCov.getZ();
    auto trackBachDCAToPV = std::sqrt(decay3body.trackDCAxyToPV[2] * decay3body.trackDCAxyToPV[2] + decay3body.trackDCAzToPV[2] * decay3body.trackDCAzToPV[2]);
    if constexpr (useSelections) {
      if (trackBachDCAToPV < decay3bodyselections.minDCABachToPV) {
        decay3body = {};
        return false;
      }
    }

    //_______________________________________________________________________
    // fit 3body vertex
    if (!useKFParticle) {
      fitVertexWithDCAFitter(collision, trackParCovPos, trackParCovNeg, trackParCovBach, bachelorcharge, calculateCovariance);
    } else {
      fitVertexWithKF(collision, trackParCovPos, trackParCovNeg, trackParCovBach, bachelorcharge, kfSetTopologicalConstraint, calculateCovariance);
    }

    //_______________________________________________________________________
    // get vertex information
    // daughter pT
    auto trackPosPt = std::sqrt(decay3body.trackPosMom[0] * decay3body.trackPosMom[0] + decay3body.trackPosMom[1] * decay3body.trackPosMom[1]);
    auto trackNegPt = std::sqrt(decay3body.trackNegMom[0] * decay3body.trackNegMom[0] + decay3body.trackNegMom[1] * decay3body.trackNegMom[1]);
    auto trackBachPt = std::sqrt(decay3body.trackBachMom[0] * decay3body.trackBachMom[0] + decay3body.trackBachMom[1] * decay3body.trackBachMom[1]);

    // DCA between daughters at SV
    decay3body.daughterDCAatSV = std::hypot(
                                 std::hypot(decay3body.trackPosPos[0] - decay3body.trackNegPos[0],
                                            decay3body.trackPosPos[1] - decay3body.trackNegPos[1],
                                            decay3body.trackPosPos[2] - decay3body.trackNegPos[2]),
                                 std::hypot(decay3body.trackPosPos[0] - decay3body.trackBachPos[0],
                                            decay3body.trackPosPos[1] - decay3body.trackBachPos[1],
                                            decay3body.trackPosPos[2] - decay3body.trackBachPos[2]),
                                 std::hypot(decay3body.trackNegPos[0] - decay3body.trackBachPos[0],
                                            decay3body.trackNegPos[1] - decay3body.trackBachPos[1],
                                            decay3body.trackNegPos[2] - decay3body.trackBachPos[2]));
    
    // daughter DCA to SV
    // positive daughter
    decay3body.daughterDCAtoSV[0] = std::hypot(
      decay3body.trackPosPos[0] - decay3body.position[0],
      decay3body.trackPosPos[1] - decay3body.position[1],
      decay3body.trackPosPos[2] - decay3body.position[2]);
    // negative daughter
    decay3body.daughterDCAtoSV[1] = std::hypot(
      decay3body.trackNegPos[0] - decay3body.position[0],
      decay3body.trackNegPos[1] - decay3body.position[1],
      decay3body.trackNegPos[2] - decay3body.position[2]);
    // bachelor daughter
    decay3body.daughterDCAtoSV[2] = std::hypot(
      decay3body.trackBachPos[0] - decay3body.position[0],
      decay3body.trackBachPos[1] - decay3body.position[1],
      decay3body.trackBachPos[2] - decay3body.position[2]);
    
    //_____________________________________________________
    // selections after vertex fit
    if constexpr (useSelections) {
      // daughter pT
      // proton
      if (isMatter && (trackPosPt < decay3bodyselections.minPtProton || trackPosPt > decay3bodyselections.maxPtProton)) {
        decay3body = {};
        return false;
      } else if (!isMatter && (trackNegPt < decay3bodyselections.minPtProton || trackNegPt > decay3bodyselections.maxPtProton)) {
        decay3body = {};
        return false;
      }
      // pion
      if (isMatter && (trackNegPt < decay3bodyselections.minPtPion || trackNegPt > decay3bodyselections.maxPtPion)) {
        decay3body = {};
        return false;
      } else if (!isMatter && (trackPosPt < decay3bodyselections.minPtPion || trackPosPt > decay3bodyselections.maxPtPion)) {
        decay3body = {};
        return false;
      }
      // bachelor
      if (trackBachPt < decay3bodyselections.minPtBach || trackBachPt > decay3bodyselections.maxPtBach) {
        decay3body = {};
        return false;
      }

      // daughter DCAs at SV
      if (decay3body.daughterDCAatSV > decay3bodyselections.maxDCADauAtSV) {
        decay3body = {};
        return false;
      }

      // rapidity
      float rapidity = RecoDecay::y(std::array{decay3body.momentum[0], decay3body.momentum[1], decay3body.momentum[2]}, o2::constants::physics::MassHyperTriton);
      if (std::fabs(rapidity) > decay3bodyselections.maxRapidity) {
        decay3body = {};
        return false;
      }

      // pT
      float pT = RecoDecay::pt(std::array{decay3body.momentum[0], decay3body.momentum[1], decay3body.momentum[2]});
      if (pT < decay3bodyselections.minPt || pT > decay3bodyselections.maxPt) {
        decay3body = {};
        return false;
      }

      // mass window
      if (decay3body.mass < decay3bodyselections.minMass || decay3body.mass > decay3bodyselections.maxMass) {
        decay3body = {};
        return false;
      }

      // pointing angle
      float cpa = RecoDecay::cpa(std::array{pvX, pvY, pvZ}, std::array{decay3body.position[0], decay3body.position[1], decay3body.position[2]}, std::array{decay3body.momentum[0], decay3body.momentum[1], decay3body.momentum[2]});
      if (cpa < decay3bodyselections.minCosPA) {
        decay3body = {};
        return false;
      }

      // vertex chi2
      if (decay3body.chi2 > decay3bodyselections.maxChi2) {
        decay3body = {};
        return false;
      }

      // ctau
      float P = RecoDecay::sqrtSumOfSquares(decay3body.momentum[0], decay3body.momentum[1], decay3body.momentum[2]);
      float ctau = std::sqrt(std::pow(decay3body.position[0] - pvX, 2) + std::pow(decay3body.position[1] - pvY, 2) + std::pow(decay3body.position[2] - pvZ, 2)) / (P + 1E-10) * o2::constants::physics::MassHyperTriton;
      if (ctau < decay3bodyselections.minCtau || ctau > decay3bodyselections.maxCtau) {
        decay3body = {};
        return false;
      }
    }

    //_______________________________________________________________________
    // fill remaining candidate information
    // daughter PID
    if (isMatter) {
      decay3body.tpcNsigma[0] = trackPos.tpcNsigmaPr();
      decay3body.tpcNsigma[1] = trackNeg.tpcNsigmaPi();
    } else if (!isMatter) {
      decay3body.tpcNsigma[0] = trackNeg.tpcNsigmaPr();
      decay3body.tpcNsigma[1] = trackPos.tpcNsigmaPi();
    }
    decay3body.tpcNsigma[2] = trackBach.tpcNsigmaDeuteron();
    decay3body.tpcNsigma[3] = trackBach.tpcNsigmaPi();
    // recalculated bachelor TOF PID
    decay3body.tofNsigmaDeuteron = todNsigmaDeuteron;

    // average ITS cluster size of deuteron track
    double averageClusterSizeDeuteron(0);
    int nCls(0);
    for (int i = 0; i < 7; i++) {
      int clusterSize = trackBach.itsClsSizeInLayer(i);
      averageClusterSizeDeuteron += static_cast<double>(clusterSize);
      if (clusterSize > 0)
        nCls++;
    }
    averageClusterSizeDeuteron = averageClusterSizeDeuteron / static_cast<double>(nCls);
    decay3body.averageITSClusSizeDeuteron = averageITSClusSizeDeuteron;

    // PID for tracking of deuteron track
    decay3body.pidForTrackingDeuteron = trackBach.pidForTracking();

    // tracked cluster size
    decay3body.trackedClSize = trackedClSize;

    return true;
  }

  //_______________________________________________________________________
  // build Decay3body from V0 and bachelor
  // template <>
  // bool buildDecay3BodyCandidate()
  // {
  //   return true;
  // }

  //___________________________________________________________________________________
  // functionality to fit 3body vertex with KFParticle and fill vertex information
  template <typename TCollision, typename TTrackParCov>
  void fitVertexWithKF(TCollision const& collision,
                       TTrackParCov const& trackParCovPos,
                       TTrackParCov const& trackParCovNeg,
                       TTrackParCov const& trackParCovBach,
                       int bachelorcharge = 1,
                       bool kfSetTopologicalConstraint = true,
                       bool calculateCovariance = true)
  {
    bool isMatter = trackParCovBach.sign() > 0 ? true : false;

    // initialise KF primary vertex
    KPFVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle kfpv(kfpVertex);

    // create KFParticle objects
    KFParticle kfpProton, kfpPion, kfpDeuteron;
    if (isMatter) {
      kfpProton = createKFParticleFromTrackParCov(trackParCovPos, trackParCovPos.sign(), constants::physics::MassProton);
      kfpPion = createKFParticleFromTrackParCov(trackParCovNeg, trackParCovNeg.sign(), constants::physics::MassPionCharged);
    } else if (!isMatter) {
      kfpProton = createKFParticleFromTrackParCov(trackParCovNeg, trackParCovNeg.sign(), constants::physics::MassProton);
      kfpPion = createKFParticleFromTrackParCov(trackParCovPos, trackParCovPos.sign(), constants::physics::MassPionCharged);
    }
    kfpDeuteron = createKFParticleFromTrackParCov(trackParCovBach, trackParCovBach.sign(), constants::physics::MassDeuteron);

    // construct vertex
    KFParticle KFH3L;
    int nDaughters3body = 3;
    const KFParticle* Daughters3body[3] = {&kfpProton, &kfpPion, &kfpDeuteron};
    KFH3L.SetConstructMethod(2);
    try {
      KFH3L.Construct(Daughters3body, nDaughters3body);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create Hyper triton 3-body vertex." << e.what();
      return;
    }

    // topological constraint
    if (kfSetTopologicalConstraint) {
      KFH3L.SetProductionVertex(kfpv);
      KFH3L.TransportToDecayVertex();
    }

    // get vertex position and momentum
    decay3body.position[0] = KFH3L.GetX();
    decay3body.position[1] = KFH3L.GetY();
    decay3body.position[2] = KFH3L.GetZ();
    decay3body.momentum[0] = KFH3L.GetPx();
    decay3body.momentum[1] = KFH3L.GetPy();
    decay3body.momentum[2] = KFH3L.GetPz();

    // get charge
    decay3body.charge = KFH3L.GetQ();

    // transport all daughter tracks to hypertriton vertex
    kfpProton.TransportToPoint(decay3body.position);
    kfpPion.TransportToPoint(decay3body.position);
    kfpDeuteron.TransportToPoint(decay3body.position);

    // daughter momenta and positions
    if (isMatter) {
      decay3body.trackPosMom[0] = kfpProton.GetPx();
      decay3body.trackPosMom[1] = kfpProton.GetPy();
      decay3body.trackPosMom[2] = kfpProton.GetPz();
      decay3body.trackNegMom[0] = kfpPion.GetPx();
      decay3body.trackNegMom[1] = kfpPion.GetPy();
      decay3body.trackNegMom[2] = kfpPion.GetPz();
      decay3body.trackPosPos[0] = kfpProton.GetX();
      decay3body.trackPosPos[1] = kfpProton.GetY();
      decay3body.trackPosPos[2] = kfpProton.GetZ();
      decay3body.trackNegPos[0] = kfpPion.GetX();
      decay3body.trackNegPos[1] = kfpPion.GetY();
      decay3body.trackNegPos[2] = kfpPion.GetZ();
    } else {
      decay3body.trackPosMom[0] = kfpPion.GetPx();
      decay3body.trackPosMom[1] = kfpPion.GetPy();
      decay3body.trackPosMom[2] = kfpPion.GetPz();
      decay3body.trackNegMom[0] = kfpProton.GetPx();
      decay3body.trackNegMom[1] = kfpProton.GetPy();
      decay3body.trackNegMom[2] = kfpProton.GetPz();
      decay3body.trackPosPos[0] = kfpPion.GetX();
      decay3body.trackPosPos[1] = kfpPion.GetY();
      decay3body.trackPosPos[2] = kfpPion.GetZ();
      decay3body.trackNegPos[0] = kfpProton.GetX();
      decay3body.trackNegPos[1] = kfpProton.GetY();
      decay3body.trackNegPos[2] = kfpProton.GetZ();
    }
    decay3body.trackBachMom[0] = kfpDeuteron.GetPx();
    decay3body.trackBachMom[1] = kfpDeuteron.GetPy();
    decay3body.trackBachMom[2] = kfpDeuteron.GetPz();
    decay3body.trackBachPos[0] = kfpDeuteron.GetX();
    decay3body.trackBachPos[1] = kfpDeuteron.GetY();
    decay3body.trackBachPos[2] = kfpDeuteron.GetZ();
    for (int i = 0; i < 3; i++) {
      trackBachMom[i] *= bachelorcharge;
    }

    // candidate mass
    float mass, massErr;
    KFH3L.GetMass(mass, massErr);
    decay3body.mass = mass;
    
    // vertex chi2
    decay3body.chi2 = KFH3L.GetChi2() / KFH3L.GetNDF();

    // caluclate covariance matrices
    if (calculateCovariance) {
      // candidate covariance matrix
      o2::gpu::gpustd::array<float, 21> covKF;
      for (int i = 0; i < 21; i++) { // get covariance matrix elements (lower triangle)
        covKF[i] = KFH3L.GetCovariance(i);
        decay3body.covariance[i] = covKF[i];
      }
      // daughter track covariance matrices
      for (int i = 0; i < 21; i++) { // get covariance matrix elements (lower triangle)
        decay3body.trackPosCov[i] = kfpProton.GetCovariance(i);
        decay3body.trackNegCov[i] = kfpPion.GetCovariance(i);
        decay3body.trackBachCov[i] = kfpDeuteron.GetCovariance(i);
      }
    }

    return;
  }

  //_______________________________________________________________________
  // functionality to fit 3body vertex with DCAFitter
  template <typename TCollision, typename TTrackParCov>
  void fitVertexWithDCAFitter(TCollision const& collision,
                              TTrackParCov const& trackParCovPos,
                              TTrackParCov const& trackParCovNeg,
                              TTrackParCov const& trackParCovBach
                              int bachelorcharge = 1,
                              bool calculateCovariance = true)
  {
    bool isMatter = trackParCovBach.sign() > 0 ? true : false;

    // fit the vertex
    int n3bodyVtx = fitter3body.process(trackParCovPos, trackParCovNeg, trackParCovBach);
    if (n3bodyVtx == 0) { // discard this pair
      return;
    }

    // get vertex position
    const auto& vtxXYZ = fitter3body.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      decay3body.position[i] = vtxXYZ[i];
    }

    // get daughter momenta
    const auto& propagatedTrackPos = fitter3body.getTrack(0);
    const auto& propagatedTrackNeg = fitter3body.getTrack(1);
    const auto& propagatedTrackBach = fitter3body.getTrack(2);
    propagatedTrackPos.getPxPyPzGlo(trackPosMom);
    propagatedTrackNeg.getPxPyPzGlo(trackNegMom);
    propagatedTrackBach.getPxPyPzGlo(trackBachMom);
    for (int i = 0; i < 3; i++) {
      trackBachMom[i] *= bachelorcharge;
    }

    // calculate candidate momentum
    decay3body.momentum[0] = trackPosMom[0] + trackNegMom[0] + trackBachMom[0];
    decay3body.momentum[1] = trackPosMom[1] + trackNegMom[1] + trackBachMom[1];
    decay3body.momentum[2] = trackPosMom[2] + trackNegMom[2] + trackBachMom[2];

    // candidate mass
    if (isMatter) {
      decay3body.mass = RecoDecay::m(
        std::array{std::array{trackPosMom[0], trackPosMom[1], trackPosMom[2]}, 
                    std::array{trackNegMom[0], trackNegMom[1], trackNegMom[2]}, 
                    std::array{trackBachMom[0], trackBachMom[1], trackBachMom[2]}}, 
        std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
    } else {
      decay3body.mass = RecoDecay::m(
        std::array{std::array{trackPosMom[0], trackPosMom[1], trackPosMom[2]}, 
                   std::array{trackNegMom[0], trackNegMom[1], trackNegMom[2]}, 
                   std::array{trackBachMom[0], trackBachMom[1], trackBachMom[2]}}, 
        std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron});
    }
    
    // vertex chi2 at PCA
    decay3body.chi2 = fitter3body.getChhi2AtPCACandidate():

    // caluclate covariance matrices
    if (calculateCovariance) {
      // candidate position covariance matrix
      auto covVtxV = fitter3body.calcPCACovMatrix(0);
      decay3body.covariance[0] = covVtxV(0, 0);
      decay3body.covariance[1] = covVtxV(1, 0);
      decay3body.covariance[2] = covVtxV(1, 1);
      decay3body.covariance[3] = covVtxV(2, 0);
      decay3body.covariance[4] = covVtxV(2, 1);
      decay3body.covariance[5] = covVtxV(2, 2);
      // daughter covariance matrices
      std::array<float, 21> covTpos = {0.};
      std::array<float, 21> covTneg = {0.};
      std::array<float, 21> covTbach = {0.};
      propagatedTrackPos.getCovXYZPxPyPzGlo(covTpos);
      propagatedTrackNeg.getCovXYZPxPyPzGlo(covTneg);
      propagatedTrackBach.getCovXYZPxPyPzGlo(covTbach);
      for (int i = 0; i < 21; i++) {
        decay3body.trackPosCov[i] = covTpos[i];
        decay3body.trackNegCov[i] = covTneg[i];
        decay3body.trackBachCov[i] = covTbach[i];
      }
      // candidate momentum covairance matrix
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        decay3body.covariance[MomInd[i]] = covTpos[MomInd[i]] + covTneg[MomInd[i]] + covTbach[MomInd[i]];
      }
      /// WARNING: position-momentum covariances are not calculated in the DCAFitter - remain zero
    }

    return;
  }

 private:
  // internal helper to calculate DCA (3D) of a straight line to a given PV analytically
  float CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }

};

} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_DECAY3BODYBUILDERHELPER_H_