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
/// \author Carolina Reetz <c.reetz@cern.ch>, Heidelberg University

#ifndef TOOLS_KFPARTICLE_KFUTILITIES_H_
#define TOOLS_KFPARTICLE_KFUTILITIES_H_

#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (required by KFParticle)
#endif

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <Framework/Logger.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <TDatabasePDG.h> // FIXME

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>
#include <KFParticleBase.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>

constexpr float ArbitrarySmallNumber{1e-8f};
constexpr float ArbitraryHugeNumber{1e8f};
constexpr float UndefLdLValue{9999.f};
constexpr int NumberOfMomentumComponents{3};
constexpr int NumberOfCovMatrixComponents{21};

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
  std::array<float, 3> trkPositionPar{};
  std::array<float, 3> trkMomentumPar{};
  std::array<float, NumberOfCovMatrixComponents> trkCovMatrix{};
  trackparCov.getXYZGlo(trkPositionPar);
  trackparCov.getPxPyPzGlo(trkMomentumPar);
  trackparCov.getCovXYZPxPyPzGlo(trkCovMatrix);
  const float trkParKf[6] = {trkPositionPar[0], trkPositionPar[1], trkPositionPar[2],
                             trkMomentumPar[0], trkMomentumPar[1], trkMomentumPar[2]};
  KFPTrack kfpTrack;
  kfpTrack.SetParameters(trkParKf);
  kfpTrack.SetCovarianceMatrix(trkCovMatrix.data());
  kfpTrack.SetCharge(trackSign);
  kfpTrack.SetNDF(tpcNClsFound - 5);
  kfpTrack.SetChi2(static_cast<float>(tpcNClsFound) * tpcChi2NCl);
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

/// @brief Function to create a KFParticle from a o2::track::TrackParametrizationWithError track
/// @tparam T
/// @param trackparCov TrackParCov
/// @param charge charg of track
/// @param mass mass hypothesis
/// @return KFParticle
template <typename T>
KFParticle createKFParticleFromTrackParCov(const o2::track::TrackParametrizationWithError<T>& trackparCov, int charge, float mass)
{
  std::array<T, 3> xyz, pxpypz;
  float xyzpxpypz[6];
  trackparCov.getPxPyPzGlo(pxpypz);
  trackparCov.getXYZGlo(xyz);
  for (int i{0}; i < NumberOfMomentumComponents; ++i) {
    xyzpxpypz[i] = xyz[i];
    xyzpxpypz[i + 3] = pxpypz[i];
  }

  std::array<float, NumberOfCovMatrixComponents> cv{};
  try {
    trackparCov.getCovXYZPxPyPzGlo(cv);
  } catch (std::runtime_error& e) {
    LOG(debug) << "Failed to get cov matrix from TrackParCov" << e.what();
  }

  KFParticle kfPart;
  float mInit{}, sigmaMInit{}, m{}, sigmaM{};
  kfPart.GetMass(mInit, sigmaMInit);
  LOG(debug) << "Daughter KFParticle mass before creation: " << mInit << " +- " << sigmaMInit;

  try {
    kfPart.Create(xyzpxpypz, cv.data(), charge, mass);
  } catch (std::runtime_error& e) {
    LOG(debug) << "Failed to create KFParticle from daughter TrackParCov" << e.what();
  }

  kfPart.GetMass(m, sigmaM);
  LOG(debug) << "Daughter KFParticle mass after creation: " << m << " +- " << sigmaM;
  return kfPart;
}

/// @brief Function to create a o2::track::TrackParametrizationWithError track from a KFParticle
/// @param kfParticle KFParticle to transform
/// @param pid PID hypothesis
/// @param sign sign of the particle
/// @return o2::track::TrackParametrizationWithError track
o2::track::TrackParCov getTrackParCovFromKFP(const KFParticle& kfParticle, const o2::track::PID pid, const int sign)
{
  std::array<float, 3> xyz{}, pxpypz{};
  std::array<float, NumberOfCovMatrixComponents> cv{};

  // get parameters from kfParticle
  xyz[0] = kfParticle.GetX();
  xyz[1] = kfParticle.GetY();
  xyz[2] = kfParticle.GetZ();
  pxpypz[0] = kfParticle.GetPx();
  pxpypz[1] = kfParticle.GetPy();
  pxpypz[2] = kfParticle.GetPz();

  // set covariance matrix elements (lower triangle)
  for (int i = 0; i < NumberOfCovMatrixComponents; i++) {
    cv[i] = kfParticle.GetCovariance(i);
  }

  // create TrackParCov track
  o2::track::TrackParCov track = o2::track::TrackParCov(xyz, pxpypz, cv, sign, true, pid);
  return track;
}

/// @brief Cosine of pointing angle from KFParticles
/// @param kfp KFParticle
/// @param PV KFParticle primary vertex
/// @return cpa
float cpaFromKF(KFParticle kfp, KFParticle PV)
{
  float xVtxP{}, yVtxP{}, zVtxP{}, xVtxS{}, yVtxS{}, zVtxS{}, px{}, py{}, pz{};

  xVtxP = PV.GetX();
  yVtxP = PV.GetY();
  zVtxP = PV.GetZ();

  xVtxS = kfp.GetX();
  yVtxS = kfp.GetY();
  zVtxS = kfp.GetZ();

  px = kfp.GetPx();
  py = kfp.GetPy();
  pz = kfp.GetPz();

  const float cpa = static_cast<float>(RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}));
  return cpa;
}

/// @brief Cosine of pointing angle in xy plane from KFParticles
/// @param kfp KFParticle
/// @param PV KFParticle primary vertex
/// @return cpa in xy
float cpaXYFromKF(KFParticle kfp, KFParticle PV)
{
  float xVtxP{}, yVtxP{}, xVtxS{}, yVtxS{}, px{}, py{};

  xVtxP = PV.GetX();
  yVtxP = PV.GetY();

  xVtxS = kfp.GetX();
  yVtxS = kfp.GetY();

  px = kfp.GetPx();
  py = kfp.GetPy();

  const float cpaXY = static_cast<float>(RecoDecay::cpaXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{px, py}));
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
float cosThetaStarFromKF(int iProng, int pdgvtx, int pdgprong0, int pdgprong1, KFParticle kfpprong0, KFParticle kfpprong1)
{
  float px0{}, py0{}, pz0{}, px1{}, py1{}, pz1{};

  px0 = kfpprong0.GetPx();
  py0 = kfpprong0.GetPy();
  pz0 = kfpprong0.GetPz();

  px1 = kfpprong1.GetPx();
  py1 = kfpprong1.GetPy();
  pz1 = kfpprong1.GetPz();
  std::array<double, 2> m = {0., 0.};
  m[0] = TDatabasePDG::Instance()->GetParticle(pdgprong0)->Mass();           // FIXME: Get from the PDG service of the common header
  m[1] = TDatabasePDG::Instance()->GetParticle(pdgprong1)->Mass();           // FIXME: Get from the PDG service of the common header
  const double mTot = TDatabasePDG::Instance()->GetParticle(pdgvtx)->Mass(); // FIXME: Get from the PDG service of the common header

  const float cosThetastar = static_cast<float>(RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m, mTot, iProng));
  return cosThetastar;
}

/// Calculates impact parameter in the bending plane of the particle w.r.t. a point

/// @brief Impact parameter in the bending plane of the particle w.r.t. a point
/// @param kfpParticle KFParticle
/// @param Vertex KFParticle vertex
/// @return impact parameter
float impParXYFromKF(KFParticle kfpParticle, KFParticle Vertex)
{
  float xVtxP{}, yVtxP{}, zVtxP{}, xVtxS{}, yVtxS{}, zVtxS{}, px{}, py{}, pz{};

  xVtxP = Vertex.GetX();
  yVtxP = Vertex.GetY();
  zVtxP = Vertex.GetZ();

  xVtxS = kfpParticle.GetX();
  yVtxS = kfpParticle.GetY();
  zVtxS = kfpParticle.GetZ();

  px = kfpParticle.GetPx();
  py = kfpParticle.GetPy();
  pz = kfpParticle.GetPz();

  const float impParXY = static_cast<float>(RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}));
  return impParXY;
}

/// @brief distance between production vertex and decay vertex normalised by the uncertainty
/// @param kfpParticle KFParticle
/// @param PV KFParticle primary vertex
/// @return l/delta l
float ldlFromKF(KFParticle kfpParticle, KFParticle PV)
{
  const float dxParticle = PV.GetX() - kfpParticle.GetX();
  const float dyParticle = PV.GetY() - kfpParticle.GetY();
  const float dzParticle = PV.GetZ() - kfpParticle.GetZ();
  float lParticle = std::sqrt(dxParticle * dxParticle + dyParticle * dyParticle + dzParticle * dzParticle);
  float dlParticle = (PV.GetCovariance(0) + kfpParticle.GetCovariance(0)) * dxParticle * dxParticle + (PV.GetCovariance(2) + kfpParticle.GetCovariance(2)) * dyParticle * dyParticle + (PV.GetCovariance(5) + kfpParticle.GetCovariance(5)) * dzParticle * dzParticle + 2 * ((PV.GetCovariance(1) + kfpParticle.GetCovariance(1)) * dxParticle * dyParticle + (PV.GetCovariance(3) + kfpParticle.GetCovariance(3)) * dxParticle * dzParticle + (PV.GetCovariance(4) + kfpParticle.GetCovariance(4)) * dyParticle * dzParticle);
  if (std::fabs(lParticle) < ArbitrarySmallNumber)
    lParticle = ArbitrarySmallNumber;
  dlParticle = dlParticle < 0.f ? ArbitraryHugeNumber : std::sqrt(dlParticle) / lParticle;
  if (dlParticle == 0.)
    return UndefLdLValue;
  return lParticle / dlParticle;
}

/// @brief distance between production vertex and decay vertex normalised by the uncertainty in xy plane
/// @param kfpParticle KFParticle
/// @param PV KFParticle primary vertex
/// @return l/delta l in xy plane
float ldlXYFromKF(KFParticle kfpParticle, KFParticle PV)
{
  const float dxParticle = PV.GetX() - kfpParticle.GetX();
  const float dyParticle = PV.GetY() - kfpParticle.GetY();
  float lParticle = std::sqrt(dxParticle * dxParticle + dyParticle * dyParticle);
  float dlParticle = (PV.GetCovariance(0) + kfpParticle.GetCovariance(0)) * dxParticle * dxParticle + (PV.GetCovariance(2) + kfpParticle.GetCovariance(2)) * dyParticle * dyParticle + 2 * ((PV.GetCovariance(1) + kfpParticle.GetCovariance(1)) * dxParticle * dyParticle);
  if (std::fabs(lParticle) < ArbitrarySmallNumber)
    lParticle = ArbitrarySmallNumber;
  dlParticle = dlParticle < 0.f ? ArbitraryHugeNumber : std::sqrt(dlParticle) / lParticle;
  if (dlParticle == 0.)
    return UndefLdLValue;
  return lParticle / dlParticle;
}

/// @brief squared distance between track and primary vertex normalised by its uncertainty evaluated in matrix form
/// @param track KFParticle track (must be passed as a copy)
/// @param vtx KFParticle primary vertex
/// @return chi2 to primary vertex
float kfCalculateChi2ToPrimaryVertex(KFParticle track, const KFParticle& vtx)
{
  const float pvPoint[3] = {vtx.X(), vtx.Y(), vtx.Z()};

  track.TransportToPoint(pvPoint);
  return track.GetDeviationFromVertex(vtx);
}

/// @brief prong's momentum in the secondary (decay) vertex
/// @param track KFParticle track (must be passed as a copy)
/// @param vtx KFParticle secondary vertex
/// @return array with components of prong's momentum in the secondary (decay) vertex
std::array<float, 3> kfCalculateProngMomentumInSecondaryVertex(KFParticle track, const KFParticle& vtx)
{
  const float svPoint[3] = {vtx.X(), vtx.Y(), vtx.Z()};

  track.TransportToPoint(svPoint);
  return {track.GetPx(), track.GetPy(), track.GetPz()};
}

/// @brief distance of closest approach between two tracks, cm
/// @param track1 KFParticle first track (must be passed as a copy)
/// @param track2 KFParticle second track (must be passed as a copy)
/// @return DCA [cm] in the PCA
float kfCalculateDistanceBetweenParticles(KFParticle track1, KFParticle track2)
{
  float dS[2];
  float dsdr[4][6];
  float params1[8], params2[8];
  float covs1[36], covs2[36];
  track1.GetDStoParticle(track2, dS, dsdr);
  track1.Transport(dS[0], dsdr[0], params1, covs1);
  track2.Transport(dS[1], dsdr[3], params2, covs2);
  const float dx = params1[0] - params2[0];
  const float dy = params1[1] - params2[1];
  const float dz = params1[2] - params2[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/// @brief squared distance between two tracks normalised by its uncertainty evaluated in matrix form
/// @param track1 KFParticle first track (must be passed as a copy)
/// @param track2 KFParticle second track (must be passed as a copy)
/// @return chi2 in PCA
float kfCalculateChi2geoBetweenParticles(KFParticle track1, KFParticle track2)
{
  KFParticle kfPair;
  const KFParticle* kfDaughters[2] = {&track1, &track2};
  kfPair.SetConstructMethod(2);
  kfPair.Construct(kfDaughters, 2);

  return kfPair.Chi2() / static_cast<float>(kfPair.NDF());
}

/// @brief signed distance between primary and secondary vertex and its uncertainty, cm
/// @param candidate KFParticle decay candidate (must be passed as a copy)
/// @param vtx KFParticle primary vertex
/// @return pair of l and delta l
std::pair<float, float> kfCalculateLdL(KFParticle candidate, const KFParticle& vtx)
{
  float l{}, dl{};
  candidate.SetProductionVertex(vtx);
  candidate.KFParticleBase::GetDecayLength(l, dl);

  return std::make_pair(l, dl);
}

/// @brief Z projection of the impact parameter from the track to the primary vertex, cm
/// @param candidate KFParticle prong
/// @param vtx KFParticle primary vertex
/// @return pair of impact parameter and its error
std::pair<float, float> kfCalculateImpactParameterZ(const KFParticle& candidate, const KFParticle& vtx)
{
  float distanceToVertexXY{}, errDistanceToVertexXY{};
  candidate.GetDistanceFromVertexXY(vtx, distanceToVertexXY, errDistanceToVertexXY);
  const float distanceToVertex = candidate.GetDistanceFromVertex(vtx);
  const float chi2ToVertex = candidate.GetDeviationFromVertex(vtx);
  const float distanceToVertexZ2 = distanceToVertex * distanceToVertex - distanceToVertexXY * distanceToVertexXY;
  const float distanceToVertexZ = distanceToVertexZ2 > 0 ? std::sqrt(distanceToVertexZ2) : -std::sqrt(-distanceToVertexZ2);
  const float errDistanceToVertexZ2 = (distanceToVertex * distanceToVertex * distanceToVertex * distanceToVertex / chi2ToVertex - distanceToVertexXY * distanceToVertexXY * errDistanceToVertexXY * errDistanceToVertexXY) / distanceToVertexZ2;
  const float errDistanceToVertexZ = errDistanceToVertexZ2 > 0 ? std::sqrt(errDistanceToVertexZ2) : -std::sqrt(-errDistanceToVertexZ2);
  return std::make_pair(distanceToVertexZ, errDistanceToVertexZ);
}

#endif // TOOLS_KFPARTICLE_KFUTILITIES_H_
