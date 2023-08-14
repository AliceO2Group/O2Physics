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

/// \file KFUtilities.h
/// \brief Utilities needed for the KFParticle package
///
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt

#ifndef TOOLS_KFPARTICLE_KFUTILITIES_H_
#define TOOLS_KFPARTICLE_KFUTILITIES_H_

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

#include "Common/Core/RecoDecay.h"

/// @brief Function to create a KFPVertex from the collision table in the AO2Ds.
/// The Multiplicity table is required to set the number of real PV Contributors
/// This function works only for Run 3 data.
/// In Run 2 converted data the number of real PV contruíbutors is not available. Switch to numContrib.
/// @tparam T
/// @param collision Collision from aod::Collisions, aod::Mults
/// @return
template <typename T>
KFPVertex createKFPVertexFromCollision(const T& collision)
{
  KFPVertex kfpVertex;
  kfpVertex.SetXYZ(collision.posX(), collision.posY(), collision.posZ());
  kfpVertex.SetCovarianceMatrix(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
  kfpVertex.SetChi2(collision.chi2());
  kfpVertex.SetNDF(2 * collision.numContrib() - 3);
  kfpVertex.SetNContributors(collision.numContrib());
  return kfpVertex;
}

/// @brief Function to create a KFPTrack from o2::track::TrackParametrizationWithError tracks. The Covariance matrix is needed.
/// @param track Track from o2::track::TrackParametrizationWithError
/// @return KFPTrack
KFPTrack createKFPTrack(const o2::track::TrackParametrizationWithError<float>& trackparCov,
                        int16_t trackSign,
                        int16_t tpcNClsFound,
                        float tpcChi2NCl)
{
  std::array<float, 3> trkpos_par;
  std::array<float, 3> trkmom_par;
  std::array<float, 21> trk_cov;
  trackparCov.getXYZGlo(trkpos_par);
  trackparCov.getPxPyPzGlo(trkmom_par);
  trackparCov.getCovXYZPxPyPzGlo(trk_cov);
  float trkpar_KF[6] = {trkpos_par[0], trkpos_par[1], trkpos_par[2],
                        trkmom_par[0], trkmom_par[1], trkmom_par[2]};
  float trkcov_KF[21];
  for (int i = 0; i < 21; i++) {
    trkcov_KF[i] = trk_cov[i];
  }
  KFPTrack kfpTrack;
  kfpTrack.SetParameters(trkpar_KF);
  kfpTrack.SetCovarianceMatrix(trkcov_KF);
  kfpTrack.SetCharge(trackSign);
  kfpTrack.SetNDF(tpcNClsFound - 5);
  kfpTrack.SetChi2(tpcNClsFound * tpcChi2NCl);
  return kfpTrack;
}

/// @brief Function to create a KFPTrack from AO2D tracks. The Covariance matrix is needed.
/// @tparam T
/// @param track Track from aod::Tracks, aod::TracksExtra, aod::TracksCov
/// @return KFPTrack
template <typename T>
KFPTrack createKFPTrackFromTrack(const T& track)
{
  o2::track::TrackParametrizationWithError trackparCov;
  trackparCov = getTrackParCov(track);

  KFPTrack kfpTrack = createKFPTrack(trackparCov, track.sign(), track.tpcNClsFound(), track.tpcChi2NCl());
  return kfpTrack;
}

/// @brief Function to create a KFPTrack from o2::track::TrackParametrizationWithError tracks. The Covariance matrix is needed.
/// @param track Track from o2::track::TrackParametrizationWithError
/// @return KFPTrack
KFPTrack createKFPTrackFromTrackParCov(const o2::track::TrackParametrizationWithError<float>& trackparCov,
                                       int16_t trackSign,
                                       int16_t tpcNClsFound,
                                       float tpcChi2NCl)
{
  KFPTrack kfpTrack = createKFPTrack(trackparCov, trackSign, tpcNClsFound, tpcChi2NCl);
  return kfpTrack;
}

/// @brief Cosine of pointing angle from KFParticles
/// @param kfp KFParticle
/// @param PV KFParticle primary vertex
/// @return cpa
float cpaFromKF(KFParticle kfp, KFParticle PV)
{
  float xVtxP, yVtxP, zVtxP, xVtxS, yVtxS, zVtxS, px, py, pz = 0.;

  xVtxP = PV.GetX();
  yVtxP = PV.GetY();
  zVtxP = PV.GetZ();

  xVtxS = kfp.GetX();
  yVtxS = kfp.GetY();
  zVtxS = kfp.GetZ();

  px = kfp.GetPx();
  py = kfp.GetPy();
  pz = kfp.GetPz();

  float cpa = RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz});
  return cpa;
}

/// @brief Cosine of pointing angle in xy plane from KFParticles
/// @param kfp KFParticle
/// @param PV KFParticle primary vertex
/// @return cpa in xy
float cpaXYFromKF(KFParticle kfp, KFParticle PV)
{
  float xVtxP, yVtxP, xVtxS, yVtxS, px, py = 0.;

  xVtxP = PV.GetX();
  yVtxP = PV.GetY();

  xVtxS = kfp.GetX();
  yVtxS = kfp.GetY();

  px = kfp.GetPx();
  py = kfp.GetPy();

  float cpaXY = RecoDecay::cpaXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{px, py});
  return cpaXY;
}

/// @brief Cosine theta star KFParticles
/// @param ip Daughter
/// @param pdgvtx pdg code of the mother
/// @param pdgprong0 pdg code prong 0
/// @param pdgprong1 pdg code prong 1
/// @param kfpvtx KFPartice mother
/// @param kfpprong0 KFParticle Prong 0
/// @param kfpprong1 KFParticele Prong 1
/// @return cos theta star
float cosThetaStarFromKF(int ip, int pdgvtx, int pdgprong0, int pdgprong1, KFParticle kfpprong0, KFParticle kfpprong1)
{
  float px0, py0, pz0, px1, py1, pz1 = 0.;

  px0 = kfpprong0.GetPx();
  py0 = kfpprong0.GetPy();
  pz0 = kfpprong0.GetPz();

  px1 = kfpprong1.GetPx();
  py1 = kfpprong1.GetPy();
  pz1 = kfpprong1.GetPz();
  std::array<double, 2> m = {0., 0.};
  m[0] = TDatabasePDG::Instance()->GetParticle(pdgprong0)->Mass();
  m[1] = TDatabasePDG::Instance()->GetParticle(pdgprong1)->Mass();
  double mTot = TDatabasePDG::Instance()->GetParticle(pdgvtx)->Mass();
  int iProng = ip;

  float cosThetastar = RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m, mTot, iProng);
  return cosThetastar;
}

/// Calculates impact parameter in the bending plane of the particle w.r.t. a point

/// @brief Impact parameter in the bending plane of the particle w.r.t. a point
/// @param kfpParticle KFParticle
/// @param Vertex KFParticle vertex
/// @return impact parameter
float impParXYFromKF(KFParticle kfpParticle, KFParticle Vertex)
{
  float xVtxP, yVtxP, zVtxP, xVtxS, yVtxS, zVtxS, px, py, pz = 0.;

  xVtxP = Vertex.GetX();
  yVtxP = Vertex.GetY();
  zVtxP = Vertex.GetZ();

  xVtxS = kfpParticle.GetX();
  yVtxS = kfpParticle.GetY();
  zVtxS = kfpParticle.GetZ();

  px = kfpParticle.GetPx();
  py = kfpParticle.GetPy();
  pz = kfpParticle.GetPz();

  float impParXY = RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz});
  return impParXY;
}

/// @brief distance between production vertex and decay vertex normalised by the uncertainty
/// @param kfpParticle KFParticle
/// @param PV KFParticle primary vertex
/// @return l/delta l
float ldlFromKF(KFParticle kfpParticle, KFParticle PV)
{
  float dx_particle = PV.GetX() - kfpParticle.GetX();
  float dy_particle = PV.GetY() - kfpParticle.GetY();
  float dz_particle = PV.GetZ() - kfpParticle.GetZ();
  float l_particle = sqrt(dx_particle * dx_particle + dy_particle * dy_particle + dz_particle * dz_particle);
  float dl_particle = (PV.GetCovariance(0) + kfpParticle.GetCovariance(0)) * dx_particle * dx_particle + (PV.GetCovariance(2) + kfpParticle.GetCovariance(2)) * dy_particle * dy_particle + (PV.GetCovariance(5) + kfpParticle.GetCovariance(5)) * dz_particle * dz_particle + 2 * ((PV.GetCovariance(1) + kfpParticle.GetCovariance(1)) * dx_particle * dy_particle + (PV.GetCovariance(3) + kfpParticle.GetCovariance(3)) * dx_particle * dz_particle + (PV.GetCovariance(4) + kfpParticle.GetCovariance(4)) * dy_particle * dz_particle);
  if (fabs(l_particle) < 1.e-8f)
    l_particle = 1.e-8f;
  dl_particle = dl_particle < 0. ? 1.e8f : sqrt(dl_particle) / l_particle;
  if (dl_particle == 0.)
    return 9999.;
  return l_particle / dl_particle;
}

/// @brief distance between production vertex and decay vertex normalised by the uncertainty in xy plane
/// @param kfpParticle KFParticle
/// @param PV KFParticle primary vertex
/// @return l/delta l in xy plane
float ldlXYFromKF(KFParticle kfpParticle, KFParticle PV)
{
  float dx_particle = PV.GetX() - kfpParticle.GetX();
  float dy_particle = PV.GetY() - kfpParticle.GetY();
  float l_particle = sqrt(dx_particle * dx_particle + dy_particle * dy_particle);
  float dl_particle = (PV.GetCovariance(0) + kfpParticle.GetCovariance(0)) * dx_particle * dx_particle + (PV.GetCovariance(2) + kfpParticle.GetCovariance(2)) * dy_particle * dy_particle + 2 * ((PV.GetCovariance(1) + kfpParticle.GetCovariance(1)) * dx_particle * dy_particle);
  if (fabs(l_particle) < 1.e-8f)
    l_particle = 1.e-8f;
  dl_particle = dl_particle < 0. ? 1.e8f : sqrt(dl_particle) / l_particle;
  if (dl_particle == 0.)
    return 9999.;
  return l_particle / dl_particle;
}

#endif // TOOLS_KFPARTICLE_KFUTILITIES_H_
