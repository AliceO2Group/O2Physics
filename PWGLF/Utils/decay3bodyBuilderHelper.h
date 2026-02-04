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

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsVertexing/SVertexHypothesis.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/Track.h"

#include <array>
#include <cmath>
#include <cstdlib>

#ifndef HomogeneousField
#define HomogeneousField
#endif

/// includes KFParticle
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
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
  int protonID = -1;
  int pionID = -1;
  int deuteronID = -1;

  // daughter properties
  std::array<float, 3> momProton = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> momPion = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> momDeuteron = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> posProton = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> posPion = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> posDeuteron = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> trackDCAxyToPV = {0.0f, 0.0f, 0.0f};     // 0 - proton, 1 - pion, 2 - deuteron
  std::array<float, 3> trackDCAToPV = {0.0f, 0.0f, 0.0f};       // 0 - proton, 1 - pion, 2 - deuteron
  std::array<float, 3> trackDCAxyToPVprop = {0.0f, 0.0f, 0.0f}; // 0 - proton, 1 - pion, 2 - deuteron
  std::array<float, 3> trackDCAToPVprop = {0.0f, 0.0f, 0.0f};   // 0 - proton, 1 - pion, 2 - deuteron
  std::array<float, 4> tpcNsigma = {0.0f, 0.0f, 0.0f, 0.0f};    // 0 - proton, 1 - pion, 2 - deuteron, 3 - bach with pion hyp
  double tofNsigmaDeuteron = 0.0f;
  std::array<float, 3> averageITSClSize = {0.0f, 0.0f, 0.0f}; // 0 - proton, 1 - pion, 2 - deuteron
  std::array<float, 3> tpcNCl = {0.0f, 0.0f, 0.0f};           // 0 - proton, 1 - pion, 2 - deuteron
  int pidForTrackingDeuteron = 0;

  // vertex properties
  float mass;
  float massV0;
  int sign;
  float momentum[3];
  float position[3];
  float chi2 = 0.0f;
  float trackedClSize = 0.0f;
  float cosPA = 0.0f;                                        // cosine of pointing angle
  float ctau = 0.0f;                                         // ctau of the candidate
  float daughterDCAtoSVaverage = 0.0f;                       // average of quadratic sum of daughter DCAs to SV
  std::array<float, 3> daughterDCAtoSV = {0.0f, 0.0f, 0.0f}; // 0 - pos, 1 - neg, 2 - bach

  // covariance matrix
  float covProton[21] = {0.0f};
  float covPion[21] = {0.0f};
  float covDeuteron[21] = {0.0f};
  float covariance[21] = {0.0f};
};

//_______________________________________________________________________
// builder helper class
class decay3bodyBuilderHelper
{
 public:
  decay3bodyBuilderHelper()
  {
    fitter3body.setPropagateToPCA(true);
    fitter3body.setMaxR(200.); //->maxRIni3body
    fitter3body.setMinParamChange(1e-3);
    fitter3body.setMinRelChi2Change(0.9);
    fitter3body.setMaxDZIni(1e9);
    fitter3body.setMaxDXYIni(4.0f);
    fitter3body.setMaxChi2(1e9);
    fitter3body.setUseAbsDCA(true);

    fitterV0.setPropagateToPCA(true);
    fitterV0.setMaxR(200.);
    fitterV0.setMinParamChange(1e-3);
    fitterV0.setMinRelChi2Change(0.9);
    fitterV0.setMaxDZIni(1e9);
    fitterV0.setMaxChi2(1e9);
    fitterV0.setUseAbsDCA(true);

    // mag field has to be set later
    fitter3body.setBz(-999.9f); // will NOT make sense if not changed
  };

  o2::vertexing::DCAFitterN<2> fitterV0;    // 2-prong o2 dca fitter
  o2::vertexing::DCAFitterN<3> fitter3body; // 3-prong o2 dca fitter

  decay3bodyCandidate decay3body; // storage for Decay3body candidate properties

  o2::dataformats::VertexBase mMeanVertex{{0., 0., 0.}, {0.1 * 0.1, 0., 0.1 * 0.1, 0., 0., 6. * 6.}};
  o2::vertexing::SVertexHypothesis mV0Hyps; // 0 - Lambda, 1 - AntiLambda

  // decay3body candidate criteria
  struct {
    // daughter tracks
    float maxEtaDaughters;
    int minTPCNClProton;
    int minTPCNClPion;
    int minTPCNClDeuteron;
    float minDCAProtonToPV;
    float minDCAPionToPV;
    float minDCADeuteronToPV;
    float minDCAProtonToPVprop;
    float minDCAPionToPVprop;
    float minDCADeuteronToPVprop;
    float minPtProton;
    float minPtPion;
    float minPtDeuteron;
    float maxPtProton;
    float maxPtPion;
    float maxPtDeuteron;
    float maxTPCnSigma;
    double minTOFnSigmaDeuteron;
    double maxTOFnSigmaDeuteron;
    float minPDeuteronUseTOF;
    float maxDCADauToSVaverage;
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

  // SVertexer selection criteria
  struct {
    float minPt2V0;
    float maxTgl2V0;
    float maxDCAXY2ToMeanVertex3bodyV0;
    float minCosPAXYMeanVertex3bodyV0;
    float minCosPA3bodyV0;
    float maxRDiffV03body;
    float minPt3Body;
    float maxTgl3Body;
    float maxDCAXY3Body;
    float maxDCAZ3Body;
  } svertexerselections;

  //_______________________________________________________________________
  // build Decay3body from three tracks, including V0 building.
  template <typename TCollision, typename TTrack>
  bool buildDecay3BodyCandidate(TCollision const& collision,
                                TTrack const& trackProton,
                                TTrack const& trackPion,
                                TTrack const& trackDeuteron,
                                int decay3bodyIndex,
                                double tofNsigmaDeuteron,
                                float trackedClSize,
                                bool useKFParticle = false,
                                bool kfSetTopologicalConstraint = false,
                                bool useSelections = true,
                                bool useChi2Selection = true,
                                bool useTPCforPion = false,
                                bool acceptTPCOnly = false,
                                bool askOnlyITSMatch = true,
                                bool calculateCovariance = true,
                                bool isEventMixing = false,
                                bool doApplySVertexerCuts = false)
  {
    int collisionIndex = collision.globalIndex();
    float pvX = collision.posX();
    float pvY = collision.posY();
    float pvZ = collision.posZ();

    auto trackParCovProton = getTrackParCov(trackProton);
    auto trackParCovPion = getTrackParCov(trackPion);
    auto trackParCovDeuteron = getTrackParCov(trackDeuteron);

    decay3body.collisionID = collisionIndex;
    decay3body.decay3bodyID = decay3bodyIndex;
    decay3body.protonID = trackProton.globalIndex();
    decay3body.pionID = trackPion.globalIndex();
    decay3body.deuteronID = trackDeuteron.globalIndex();

    //_______________________________________________________________________
    // track selections
    if (useSelections) {
      // proton track quality
      if (trackProton.tpcNClsFound() < decay3bodyselections.minTPCNClProton) {
        decay3body = {};
        return false;
      }
      // pion track quality
      if (useTPCforPion) {
        if (trackPion.tpcNClsFound() < decay3bodyselections.minTPCNClPion) {
          decay3body = {};
          return false;
        }
      }
      // deuteron track quality
      if (trackDeuteron.tpcNClsFound() < decay3bodyselections.minTPCNClDeuteron) {
        decay3body = {};
        return false;
      }

      // track eta
      if (std::fabs(trackProton.eta()) > decay3bodyselections.maxEtaDaughters) {
        decay3body = {};
        return false;
      }
      if (std::fabs(trackPion.eta()) > decay3bodyselections.maxEtaDaughters) {
        decay3body = {};
        return false;
      }
      if (std::fabs(trackDeuteron.eta()) > decay3bodyselections.maxEtaDaughters) {
        decay3body = {};
        return false;
      }

      // TPC only
      if (!acceptTPCOnly) {
        if (askOnlyITSMatch) {
          if (!trackProton.hasITS() || !trackPion.hasITS() || !trackDeuteron.hasITS()) {
            decay3body = {};
            return false;
          }
        } else {
          bool isProtonTPCOnly = !trackProton.hasITS() && !trackProton.hasTOF() && !trackProton.hasTRD();
          bool isPionTPCOnly = !trackPion.hasITS() && !trackPion.hasTOF() && !trackPion.hasTRD();
          bool isDeuteronTPCOnly = !trackDeuteron.hasITS() && !trackDeuteron.hasTOF() && !trackDeuteron.hasTRD();
          if (isProtonTPCOnly || isPionTPCOnly || isDeuteronTPCOnly) {
            decay3body = {};
            return false;
          }
        }
      }

      // daughter TPC PID
      if (std::fabs(trackProton.tpcNSigmaPr()) > decay3bodyselections.maxTPCnSigma) {
        decay3body = {};
        return false;
      }
      if (useTPCforPion && std::fabs(trackPion.tpcNSigmaPi()) > decay3bodyselections.maxTPCnSigma) {
        decay3body = {};
        return false;
      }
      if (std::fabs(trackDeuteron.tpcNSigmaDe()) > decay3bodyselections.maxTPCnSigma) {
        decay3body = {};
        return false;
      }

      // deuteron TOF PID
      if ((tofNsigmaDeuteron < decay3bodyselections.minTOFnSigmaDeuteron || tofNsigmaDeuteron > decay3bodyselections.maxTOFnSigmaDeuteron) && trackDeuteron.p() > decay3bodyselections.minPDeuteronUseTOF) {
        decay3body = {};
        return false;
      }
    } // end of selections

    //_______________________________________________________________________
    // daughter track DCA to PV associated with decay3body --> computed with KFParticle
    float pvXY[2] = {pvX, pvY};
    float pv[3] = {pvX, pvY, pvZ};
    auto trackParCovProtonCopy = trackParCovProton;
    auto trackParCovPionCopy = trackParCovPion;
    auto trackParCovDeuteronCopy = trackParCovDeuteron;
    KFParticle kfproton = createKFParticleFromTrackParCov(trackParCovProtonCopy, trackProton.sign(), constants::physics::MassProton);
    KFParticle kfpion = createKFParticleFromTrackParCov(trackParCovPionCopy, trackPion.sign(), constants::physics::MassPionCharged);
    KFParticle kfdeuteron = createKFParticleFromTrackParCov(trackParCovDeuteronCopy, trackDeuteron.sign(), constants::physics::MassDeuteron);

    // proton DCA to PV
    decay3body.trackDCAxyToPV[0] = kfproton.GetDistanceFromVertexXY(pvXY);
    decay3body.trackDCAToPV[0] = kfproton.GetDistanceFromVertex(pv);
    // pion DCA to PV
    decay3body.trackDCAxyToPV[1] = kfpion.GetDistanceFromVertexXY(pvXY);
    decay3body.trackDCAToPV[1] = kfpion.GetDistanceFromVertex(pv);
    // deuteron DCA to PV
    decay3body.trackDCAxyToPV[2] = kfdeuteron.GetDistanceFromVertexXY(pvXY);
    decay3body.trackDCAToPV[2] = kfdeuteron.GetDistanceFromVertex(pv);
    // selection
    if (useSelections) {
      if (decay3body.trackDCAToPV[0] < decay3bodyselections.minDCAProtonToPV) {
        decay3body = {};
        return false;
      }
      if (decay3body.trackDCAToPV[1] < decay3bodyselections.minDCAPionToPV) {
        decay3body = {};
        return false;
      }
      if (decay3body.trackDCAToPV[2] < decay3bodyselections.minDCADeuteronToPV) {
        decay3body = {};
        return false;
      }
    }

    //_______________________________________________________________________
    // daughter track DCA to PV associated with decay3body --> with O2 Propagator
    o2::dataformats::VertexBase mPV;
    o2::dataformats::DCA mDcaInfoCov;
    auto trackParCovProtonCopyProp = trackParCovProton;
    auto trackParCovPionCopyProp = trackParCovPion;
    auto trackParCovDeuteronCopyProp = trackParCovDeuteron;
    mPV.setPos({pvX, pvY, pvZ});
    mPV.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

    // proton track
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovProtonCopyProp, 2.f, fitter3body.getMatCorrType(), &mDcaInfoCov);
    decay3body.trackDCAxyToPVprop[0] = mDcaInfoCov.getY();
    auto trackProtonDCAzToPVprop = mDcaInfoCov.getZ();
    decay3body.trackDCAToPVprop[0] = std::sqrt(decay3body.trackDCAxyToPVprop[0] * decay3body.trackDCAxyToPVprop[0] + trackProtonDCAzToPVprop * trackProtonDCAzToPVprop);
    if (useSelections) {
      if (decay3body.trackDCAToPVprop[0] < decay3bodyselections.minDCAProtonToPVprop) {
        decay3body = {};
        return false;
      }
    }
    // pion track
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovPionCopyProp, 2.f, fitter3body.getMatCorrType(), &mDcaInfoCov);
    decay3body.trackDCAxyToPVprop[1] = mDcaInfoCov.getY();
    auto trackPionDCAzToPVprop = mDcaInfoCov.getZ();
    decay3body.trackDCAToPVprop[1] = std::sqrt(decay3body.trackDCAxyToPVprop[1] * decay3body.trackDCAxyToPVprop[1] + trackPionDCAzToPVprop * trackPionDCAzToPVprop);
    if (useSelections) {
      if (decay3body.trackDCAToPVprop[1] < decay3bodyselections.minDCAPionToPVprop) {
        decay3body = {};
        return false;
      }
    }
    // deuteron track
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, trackParCovDeuteronCopyProp, 2.f, fitter3body.getMatCorrType(), &mDcaInfoCov);
    decay3body.trackDCAxyToPVprop[2] = mDcaInfoCov.getY();
    auto trackDeuteronDCAzToPVprop = mDcaInfoCov.getZ();
    decay3body.trackDCAToPVprop[2] = std::sqrt(decay3body.trackDCAxyToPVprop[2] * decay3body.trackDCAxyToPVprop[2] + trackDeuteronDCAzToPVprop * trackDeuteronDCAzToPVprop);
    if (useSelections) {
      if (decay3body.trackDCAToPVprop[2] < decay3bodyselections.minDCADeuteronToPVprop) {
        decay3body = {};
        return false;
      }
    }

    //_______________________________________________________________________
    // fit 3body vertex
    if (!useKFParticle) {
      fitVertexWithDCAFitter(trackProton, trackPion, trackDeuteron, calculateCovariance);
    } else {
      fitVertexWithKF(collision, trackProton, trackPion, trackDeuteron, kfSetTopologicalConstraint, calculateCovariance);
    }

    //_______________________________________________________________________
    // get vertex information
    // daughter pT
    auto trackProtonPt = std::sqrt(decay3body.momProton[0] * decay3body.momProton[0] + decay3body.momProton[1] * decay3body.momProton[1]);
    auto trackPionPt = std::sqrt(decay3body.momPion[0] * decay3body.momPion[0] + decay3body.momPion[1] * decay3body.momPion[1]);
    auto trackDeuteronPt = std::sqrt(decay3body.momDeuteron[0] * decay3body.momDeuteron[0] + decay3body.momDeuteron[1] * decay3body.momDeuteron[1]);

    // daughter DCA to SV
    // proton daughter
    decay3body.daughterDCAtoSV[0] = std::hypot(
      decay3body.posProton[0] - decay3body.position[0],
      decay3body.posProton[1] - decay3body.position[1],
      decay3body.posProton[2] - decay3body.position[2]);
    // pion daughter
    decay3body.daughterDCAtoSV[1] = std::hypot(
      decay3body.posPion[0] - decay3body.position[0],
      decay3body.posPion[1] - decay3body.position[1],
      decay3body.posPion[2] - decay3body.position[2]);
    // deuteron daughter
    decay3body.daughterDCAtoSV[2] = std::hypot(
      decay3body.posDeuteron[0] - decay3body.position[0],
      decay3body.posDeuteron[1] - decay3body.position[1],
      decay3body.posDeuteron[2] - decay3body.position[2]);

    // DCA daughters to SV average of quadratic sum
    decay3body.daughterDCAtoSVaverage = (decay3body.daughterDCAtoSV[0] * decay3body.daughterDCAtoSV[0] +
                                         decay3body.daughterDCAtoSV[1] * decay3body.daughterDCAtoSV[1] +
                                         decay3body.daughterDCAtoSV[2] * decay3body.daughterDCAtoSV[2]) /
                                        3;

    //_____________________________________________________
    // selections after vertex fit
    if (useSelections) {
      // daughter pT
      // proton
      if (trackProtonPt < decay3bodyselections.minPtProton || trackProtonPt > decay3bodyselections.maxPtProton) {
        decay3body = {};
        return false;
      }
      // pion
      if (trackPionPt < decay3bodyselections.minPtPion || trackPionPt > decay3bodyselections.maxPtPion) {
        decay3body = {};
        return false;
      }
      // deuteron
      if (trackDeuteronPt < decay3bodyselections.minPtDeuteron || trackDeuteronPt > decay3bodyselections.maxPtDeuteron) {
        decay3body = {};
        return false;
      }

      // daughter DCAs at SV
      if (decay3body.daughterDCAtoSVaverage > decay3bodyselections.maxDCADauToSVaverage) {
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

      // vertex chi2
      if (useChi2Selection && decay3body.chi2 > decay3bodyselections.maxChi2) {
        decay3body = {};
        return false;
      }
    }

    // pointing angle
    float cpa = RecoDecay::cpa(std::array{pvX, pvY, pvZ}, std::array{decay3body.position[0], decay3body.position[1], decay3body.position[2]}, std::array{decay3body.momentum[0], decay3body.momentum[1], decay3body.momentum[2]});
    if (useSelections) {
      if (cpa < decay3bodyselections.minCosPA) {
        decay3body = {};
        return false;
      }
    }
    decay3body.cosPA = cpa;

    // ctau
    float P = RecoDecay::sqrtSumOfSquares(decay3body.momentum[0], decay3body.momentum[1], decay3body.momentum[2]);
    float ctau = std::sqrt(std::pow(decay3body.position[0] - pvX, 2) + std::pow(decay3body.position[1] - pvY, 2) + std::pow(decay3body.position[2] - pvZ, 2)) / (P + 1E-10) * o2::constants::physics::MassHyperTriton;
    if (useSelections) {
      if (ctau < decay3bodyselections.minCtau || ctau > decay3bodyselections.maxCtau) {
        decay3body = {};
        return false;
      }
    }
    decay3body.ctau = ctau;

    //_______________________________________________________________________
    // SVertexer selections in case of event mixing
    if (isEventMixing && doApplySVertexerCuts) {
      applySVertexerCuts(collision, trackProton, trackPion, trackDeuteron, /*applyV0Cut = */ true);
    }

    //_______________________________________________________________________
    // fill remaining candidate information
    // daughter PID
    decay3body.tpcNsigma[0] = trackProton.tpcNSigmaPr();
    decay3body.tpcNsigma[1] = trackPion.tpcNSigmaPi();
    decay3body.tpcNsigma[2] = trackDeuteron.tpcNSigmaDe();
    decay3body.tpcNsigma[3] = trackDeuteron.tpcNSigmaPi();
    // recalculated bachelor TOF PID
    decay3body.tofNsigmaDeuteron = tofNsigmaDeuteron;

    // average ITS cluster size of daughter tracks
    double averageClusterSizeProton(0), averageClusterSizePion(0), averageClusterSizeDeuteron(0);
    int nClsProton(0), nClsPion(0), nClsDeuteron(0);
    for (int i = 0; i < 7; i++) {
      int clusterSizePr = trackProton.itsClsSizeInLayer(i);
      int clusterSizePi = trackPion.itsClsSizeInLayer(i);
      int clusterSizeDe = trackDeuteron.itsClsSizeInLayer(i);
      averageClusterSizeProton += static_cast<double>(clusterSizePr);
      averageClusterSizePion += static_cast<double>(clusterSizePi);
      averageClusterSizeDeuteron += static_cast<double>(clusterSizeDe);
      if (clusterSizePr > 0)
        nClsProton++;
      if (clusterSizePi > 0)
        nClsPion++;
      if (clusterSizeDe > 0)
        nClsDeuteron++;
    }
    averageClusterSizeProton = averageClusterSizeProton / static_cast<double>(nClsProton);
    averageClusterSizePion = averageClusterSizePion / static_cast<double>(nClsPion);
    averageClusterSizeDeuteron = averageClusterSizeDeuteron / static_cast<double>(nClsDeuteron);
    decay3body.averageITSClSize[0] = averageClusterSizeProton;
    decay3body.averageITSClSize[1] = averageClusterSizePion;
    decay3body.averageITSClSize[2] = averageClusterSizeDeuteron;

    // number of TPC clusters
    decay3body.tpcNCl[0] = trackProton.tpcNClsFound();
    decay3body.tpcNCl[1] = trackPion.tpcNClsFound();
    decay3body.tpcNCl[2] = trackDeuteron.tpcNClsFound();

    // PID for tracking of deuteron track
    decay3body.pidForTrackingDeuteron = trackDeuteron.pidForTracking();

    // tracked cluster size
    decay3body.trackedClSize = trackedClSize;

    return true;
  }

  //___________________________________________________________________________________
  // functionality to fit 3body vertex with KFParticle and fill vertex information
  template <typename TCollision, typename TTrack>
  void fitVertexWithKF(TCollision const& collision,
                       TTrack const& trackProton,
                       TTrack const& trackPion,
                       TTrack const& trackDeuteron,
                       bool kfSetTopologicalConstraint = true,
                       bool calculateCovariance = true)
  {
    // get TrackParCov daughters
    auto trackParCovProton = getTrackParCov(trackProton);
    auto trackParCovPion = getTrackParCov(trackPion);
    auto trackParCovDeuteron = getTrackParCov(trackDeuteron);

    // initialise KF primary vertex
    KFVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle kfpv(kfpVertex);

    // create KFParticle objects
    KFParticle kfpProton, kfpPion, kfpDeuteron;
    kfpProton = createKFParticleFromTrackParCov(trackParCovProton, trackProton.sign(), constants::physics::MassProton);
    kfpPion = createKFParticleFromTrackParCov(trackParCovPion, trackPion.sign(), constants::physics::MassPionCharged);
    kfpDeuteron = createKFParticleFromTrackParCov(trackParCovDeuteron, trackDeuteron.sign(), constants::physics::MassDeuteron);

    // construct V0 vertex
    KFParticle KFV0;
    int nDaughtersV0 = 2;
    const KFParticle* DaughtersV0[2] = {&kfpProton, &kfpPion};
    KFV0.SetConstructMethod(2);
    try {
      KFV0.Construct(DaughtersV0, nDaughtersV0);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create V0 vertex." << e.what();
      return;
    }

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

    // get sign
    decay3body.sign = KFH3L.GetQ() / std::abs(KFH3L.GetQ());

    // transport all daughter tracks to hypertriton vertex
    kfpProton.TransportToPoint(decay3body.position);
    kfpPion.TransportToPoint(decay3body.position);
    kfpDeuteron.TransportToPoint(decay3body.position);

    // daughter positions
    decay3body.posProton[0] = kfpProton.GetX();
    decay3body.posProton[1] = kfpProton.GetY();
    decay3body.posProton[2] = kfpProton.GetZ();
    decay3body.posPion[0] = kfpPion.GetX();
    decay3body.posPion[1] = kfpPion.GetY();
    decay3body.posPion[2] = kfpPion.GetZ();
    decay3body.posDeuteron[0] = kfpDeuteron.GetX();
    decay3body.posDeuteron[1] = kfpDeuteron.GetY();
    decay3body.posDeuteron[2] = kfpDeuteron.GetZ();

    // daughter momenta
    decay3body.momProton[0] = kfpProton.GetPx();
    decay3body.momProton[1] = kfpProton.GetPy();
    decay3body.momProton[2] = kfpProton.GetPz();
    decay3body.momPion[0] = kfpPion.GetPx();
    decay3body.momPion[1] = kfpPion.GetPy();
    decay3body.momPion[2] = kfpPion.GetPz();
    decay3body.momDeuteron[0] = kfpDeuteron.GetPx();
    decay3body.momDeuteron[1] = kfpDeuteron.GetPy();
    decay3body.momDeuteron[2] = kfpDeuteron.GetPz();

    // candidate mass
    float mass, massErr;
    KFH3L.GetMass(mass, massErr);
    decay3body.mass = mass;

    // V0 mass
    float massV0, massV0Err;
    KFV0.GetMass(massV0, massV0Err);
    decay3body.massV0 = massV0;

    // vertex chi2
    decay3body.chi2 = KFH3L.GetChi2() / KFH3L.GetNDF();

    // caluclate covariance matrices
    if (calculateCovariance) {
      // candidate covariance matrix
      std::array<float, 21> covKF;
      for (int i = 0; i < 21; i++) { // get covariance matrix elements (lower triangle)
        covKF[i] = KFH3L.GetCovariance(i);
        decay3body.covariance[i] = covKF[i];
      }
      // daughter track covariance matrices
      for (int i = 0; i < 21; i++) { // get covariance matrix elements (lower triangle)
        decay3body.covProton[i] = kfpProton.GetCovariance(i);
        decay3body.covPion[i] = kfpPion.GetCovariance(i);
        decay3body.covDeuteron[i] = kfpDeuteron.GetCovariance(i);
      }
    }

    return;
  }

  //_______________________________________________________________________
  // functionality to fit 3body vertex with DCAFitter
  template <typename TTrack>
  void fitVertexWithDCAFitter(TTrack const& trackProton,
                              TTrack const& trackPion,
                              TTrack const& trackDeuteron,
                              bool calculateCovariance = true)
  {
    // get TrackParCov daughters
    auto trackParCovProton = getTrackParCov(trackProton);
    auto trackParCovPion = getTrackParCov(trackPion);
    auto trackParCovDeuteron = getTrackParCov(trackDeuteron);

    // fit the vertex
    int n3bodyVtx = fitter3body.process(trackParCovProton, trackParCovPion, trackParCovDeuteron);
    if (n3bodyVtx == 0) { // discard this pair
      return;
    }

    // get vertex position
    const auto& vtxXYZ = fitter3body.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      decay3body.position[i] = vtxXYZ[i];
    }

    // get daughter momenta
    const auto& propagatedTrackProton = fitter3body.getTrack(0);
    const auto& propagatedTrackPion = fitter3body.getTrack(1);
    const auto& propagatedTrackDeuteron = fitter3body.getTrack(2);
    propagatedTrackProton.getPxPyPzGlo(decay3body.momProton);
    propagatedTrackPion.getPxPyPzGlo(decay3body.momPion);
    propagatedTrackDeuteron.getPxPyPzGlo(decay3body.momDeuteron);
    propagatedTrackProton.getXYZGlo(decay3body.posProton);
    propagatedTrackPion.getXYZGlo(decay3body.posPion);
    propagatedTrackDeuteron.getXYZGlo(decay3body.posDeuteron);

    // get daughter positions at vertex

    // calculate candidate momentum
    decay3body.momentum[0] = decay3body.momProton[0] + decay3body.momPion[0] + decay3body.momDeuteron[0];
    decay3body.momentum[1] = decay3body.momProton[1] + decay3body.momPion[1] + decay3body.momDeuteron[1];
    decay3body.momentum[2] = decay3body.momProton[2] + decay3body.momPion[2] + decay3body.momDeuteron[2];

    // candidate and V0 mass
    decay3body.mass = RecoDecay::m(
      std::array{std::array{decay3body.momProton[0], decay3body.momProton[1], decay3body.momProton[2]},
                 std::array{decay3body.momPion[0], decay3body.momPion[1], decay3body.momPion[2]},
                 std::array{decay3body.momDeuteron[0], decay3body.momDeuteron[1], decay3body.momDeuteron[2]}},
      std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
    decay3body.massV0 = RecoDecay::m(
      std::array{std::array{decay3body.momProton[0], decay3body.momProton[1], decay3body.momProton[2]},
                 std::array{decay3body.momPion[0], decay3body.momPion[1], decay3body.momPion[2]}},
      std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});

    // vertex chi2 at PCA
    decay3body.chi2 = fitter3body.getChi2AtPCACandidate();

    // candidate sign
    decay3body.sign = trackDeuteron.sign();

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
      std::array<float, 21> covTproton = {0.};
      std::array<float, 21> covTpion = {0.};
      std::array<float, 21> covTdeuteron = {0.};
      propagatedTrackProton.getCovXYZPxPyPzGlo(covTproton);
      propagatedTrackPion.getCovXYZPxPyPzGlo(covTpion);
      propagatedTrackDeuteron.getCovXYZPxPyPzGlo(covTdeuteron);
      for (int i = 0; i < 21; i++) {
        decay3body.covProton[i] = covTproton[i];
        decay3body.covPion[i] = covTpion[i];
        decay3body.covDeuteron[i] = covTdeuteron[i];
      }
      // candidate momentum covairance matrix
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        decay3body.covariance[MomInd[i]] = covTproton[MomInd[i]] + covTpion[MomInd[i]] + covTdeuteron[MomInd[i]];
      }
      /// WARNING: position-momentum covariances are not calculated in the DCAFitter - remain zero
    }

    return;
  }

  //_______________________________________________________________________
  // functionality to apply SVertexer cuts in case of event mixing
  template <typename TCollision, typename TTrack>
  void applySVertexerCuts(TCollision const& collision,
                          TTrack const& trackProton,
                          TTrack const& trackPion,
                          TTrack const& trackDeuteron,
                          bool applyV0Cut = true)
  {
    // get TrackParCov daughters
    auto trackParCovProton = getTrackParCov(trackProton);
    auto trackParCovPion = getTrackParCov(trackPion);
    auto trackParCovDeuteron = getTrackParCov(trackDeuteron);

    const float pidCutsLambda[o2::vertexing::SVertexHypothesis::NPIDParams] = {0., 20, 0., 5.0, 0.0, 1.09004e-03, 2.62291e-04, 8.93179e-03, 2.83121}; // Lambda
    mV0Hyps.set(o2::track::PID::Lambda, o2::track::PID::Proton, o2::track::PID::Pion, pidCutsLambda, fitter3body.getBz());

    int nV0 = fitterV0.process(trackParCovProton, trackParCovPion);
    if (nV0 == 0) {
      return;
    }

    std::array<float, 3> v0pos = {0.};
    const auto& v0vtxXYZ = fitterV0.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      v0pos[i] = v0vtxXYZ[i];
    }
    const int cand = 0;
    if (!fitterV0.isPropagateTracksToVertexDone(cand) && !fitterV0.propagateTracksToVertex(cand)) {
      return;
    }

    const auto& trProtonProp = fitterV0.getTrack(0, cand);
    const auto& trPionProp = fitterV0.getTrack(1, cand);
    std::array<float, 3> pProtonV0{}, pPionV0{};
    trProtonProp.getPxPyPzGlo(pProtonV0);
    trPionProp.getPxPyPzGlo(pPionV0);
    std::array<float, 3> pV0 = {pProtonV0[0] + pPionV0[0], pProtonV0[1] + pPionV0[1], pProtonV0[2] + pPionV0[2]};
    // Cut for Virtual V0
    float dxv0 = v0pos[0] - mMeanVertex.getX(), dyv0 = v0pos[1] - mMeanVertex.getY(), r2v0 = dxv0 * dxv0 + dyv0 * dyv0;
    float rv0 = std::sqrt(r2v0);
    float pt2V0 = pV0[0] * pV0[0] + pV0[1] * pV0[1], prodXYv0 = dxv0 * pV0[0] + dyv0 * pV0[1], tDCAXY = prodXYv0 / pt2V0;
    if (applyV0Cut && pt2V0 <= svertexerselections.minPt2V0) {
      return;
    }
    if (applyV0Cut && pV0[2] * pV0[2] / pt2V0 > svertexerselections.maxTgl2V0) { // tgLambda cut
      return;
    }

    float p2V0 = pt2V0 + pV0[2] * pV0[2], ptV0 = std::sqrt(pt2V0);
    // apply mass selections
    float p2Proton = pProtonV0[0] * pProtonV0[0] + pProtonV0[1] * pProtonV0[1] + pProtonV0[2] * pProtonV0[2], p2Pion = pPionV0[0] * pPionV0[0] + pPionV0[1] * pPionV0[1] + pPionV0[2] * pPionV0[2];
    bool good3bodyV0Hyp = false;
    float massForLambdaHyp = mV0Hyps.calcMass(p2Proton, p2Pion, p2V0);
    if (massForLambdaHyp - mV0Hyps.getMassV0Hyp() < mV0Hyps.getMargin(ptV0)) {
      good3bodyV0Hyp = true;
    }
    if (applyV0Cut && !good3bodyV0Hyp) {
      return;
    }

    float dcaX = dxv0 - pV0[0] * tDCAXY, dcaY = dyv0 - pV0[1] * tDCAXY, dca2 = dcaX * dcaX + dcaY * dcaY;
    float cosPAXY = prodXYv0 / rv0 * ptV0;
    if (applyV0Cut && dca2 > svertexerselections.maxDCAXY2ToMeanVertex3bodyV0) {
      return;
    }
    // FIXME: V0 cosPA cut to be investigated
    if (applyV0Cut && cosPAXY < svertexerselections.minCosPAXYMeanVertex3bodyV0) {
      return;
    }
    // Check: CosPA Cut of Virtual V0 may not be used since the V0 may be based on another PV
    float dx = v0pos[0] - collision.posX();
    float dy = v0pos[1] - collision.posY();
    float dz = v0pos[2] - collision.posZ();
    float prodXYZv0 = dx * pV0[0] + dy * pV0[1] + dz * pV0[2];
    float v0CosPA = prodXYZv0 / std::sqrt((dx * dx + dy * dy + dz * dz) * p2V0);
    if (applyV0Cut && v0CosPA < svertexerselections.minCosPA3bodyV0) {
      return;
    }

    // 3body vertex
    int n3bodyVtx = fitter3body.process(trackParCovProton, trackParCovPion, trackParCovDeuteron);
    if (n3bodyVtx == 0) { // discard this pair
      return;
    }
    const auto& vertexXYZ = fitter3body.getPCACandidatePos();
    std::array<float, 3> pos = {0.};
    for (int i = 0; i < 3; i++) {
      pos[i] = vertexXYZ[i];
    }

    std::array<float, 3> pProton = {0.}, pPion = {0.}, pDeuteron{0.};
    const auto& propagatedTrackProton = fitter3body.getTrack(0);
    const auto& propagatedTrackPion = fitter3body.getTrack(1);
    const auto& propagatedTrackDeuteron = fitter3body.getTrack(2);
    propagatedTrackProton.getPxPyPzGlo(pProton);
    propagatedTrackPion.getPxPyPzGlo(pPion);
    propagatedTrackDeuteron.getPxPyPzGlo(pDeuteron);
    std::array<float, 3> p3B = {pProton[0] + pPion[0] + pDeuteron[0], pProton[1] + pPion[1] + pDeuteron[1], pProton[2] + pPion[2] + pDeuteron[2]};

    float r3body = std::hypot(pos[0], pos[1]);
    if (r3body < 0.5) {
      return;
    }

    // Cut for the compatibility of V0 and 3body vertex
    float deltaR = std::abs(rv0 - r3body);
    if (deltaR > svertexerselections.maxRDiffV03body) {
      return;
    }

    float pt3B = std::hypot(p3B[0], p3B[1]);
    if (pt3B < svertexerselections.minPt3Body) { // pt cut
      return;
    }
    if (p3B[2] / pt3B > svertexerselections.maxTgl3Body) { // tgLambda cut
      return;
    }

    // H3L DCA Check
    auto track3B = o2::track::TrackParCov(vertexXYZ, p3B, trackDeuteron.sign());
    o2::dataformats::DCA dca;
    if (!track3B.propagateToDCA({{collision.posX(), collision.posY(), collision.posZ()}, {collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()}}, fitter3body.getBz(), &dca, 5.) ||
        std::abs(dca.getY()) > svertexerselections.maxDCAXY3Body || std::abs(dca.getZ()) > svertexerselections.maxDCAZ3Body) {
      return;
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
