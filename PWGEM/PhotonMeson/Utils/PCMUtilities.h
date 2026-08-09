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

/// \file PCMUtilities.h
/// \brief helper functions commonly used for PCM analyses.
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CommonConstants/MathConstants.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/Concepts.h>
#include <ReconstructionDataFormats/HelixHelper.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <Math/GenVector/DisplacementVector2D.h> // IWYU pragma: keep (for rotate)
#include <Math/Vector2D.h>                       // IWYU pragma: keep (do not replace with Math/Vector2Dfwd.h)
#include <Math/Vector2Dfwd.h>

#include <GPUROOTCartesianFwd.h>

#include <array>
#include <cmath>

//_______________________________________________________________________
inline bool checkAP(const float alpha, const float qt, const float alpha_max = 0.95, const float qt_max = 0.05)
{
  float ellipse = std::pow(alpha / alpha_max, 2) + std::pow(qt / qt_max, 2);
  return (ellipse < 1.0);
}
//_______________________________________________________________________
inline float v0_alpha(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
{
  float momTot = RecoDecay::p(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  float lQlNeg = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
  float lQlPos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
  return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // longitudinal momentum asymmetry of v0
}
//_______________________________________________________________________
inline float v0_qt(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
{
  float momTot = RecoDecay::p2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  float dp = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
  return std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qt of v0
}
//_______________________________________________________________________
template <typename TrackPrecision = float>
inline void Vtx_recalculationParCov(o2::base::Propagator* prop, const o2::track::TrackParametrizationWithError<TrackPrecision>& trackPosInformation, const o2::track::TrackParametrizationWithError<TrackPrecision>& trackNegInformation, std::array<float, 3>& xyz, o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE)
{
  float bz = prop->getNominalBz();

  o2::track::TrackAuxPar helixPos(trackPosInformation, bz); // This object is a descendant of a CircleXY and stores cirlce information with respect to the magnetic field. This object uses functions and information of the o2::track::TrackParametrizationWithError<TrackPrecision> object (positive)
  o2::track::TrackAuxPar helixNeg(trackNegInformation, bz); // This object is a descendant of a CircleXY and stores cirlce information with respect to the magnetic field. This object uses functions and information of the o2::track::TrackParametrizationWithError<TrackPrecision> object (negative)

  xyz[0] = (helixPos.xC * helixNeg.rC + helixNeg.xC * helixPos.rC) / (helixPos.rC + helixNeg.rC); // This calculates the coordinates of the conversion point as an weighted average of the two helix centers. xC and yC should be the global coordinates for the helix center as far as I understand. But you can double check the code of trackPosInformation.getCircleParamsLoc
  xyz[1] = (helixPos.yC * helixNeg.rC + helixNeg.yC * helixPos.rC) / (helixPos.rC + helixNeg.rC); // If this calculation doesn't work check if the rotateZ function, because the "documentation" says I get global coordinates but maybe i don't.

  // I am unsure about the Z calculation but this is how it is done in AliPhysics as far as I understand
  auto trackPosInformationCopy = trackPosInformation;
  auto trackNegInformationCopy = trackNegInformation;

  // I think this calculation gets the closest point on the track to the conversion point
  // This alpha is a different alpha than the usual alpha and I think it is the angle between X axis and conversion point
  float alphaPos = o2::constants::math::PI + std::atan2(-(xyz[1] - helixPos.yC), -(xyz[0] - helixPos.xC));
  float alphaNeg = o2::constants::math::PI + std::atan2(-(xyz[1] - helixNeg.yC), -(xyz[0] - helixNeg.xC));

  float vertexXPos = helixPos.xC + helixPos.rC * std::cos(alphaPos);
  float vertexYPos = helixPos.yC + helixPos.rC * std::sin(alphaPos);
  float vertexXNeg = helixNeg.xC + helixNeg.rC * std::cos(alphaNeg);
  float vertexYNeg = helixNeg.yC + helixNeg.rC * std::sin(alphaNeg);

  ROOT::Math::XYVector vertexPos(vertexXPos, vertexYPos);
  ROOT::Math::XYVector vertexNeg(vertexXNeg, vertexYNeg);

  // Convert to local coordinate system
  vertexPos.Rotate(-trackPosInformationCopy.getAlpha());
  vertexNeg.Rotate(-trackNegInformationCopy.getAlpha());

  prop->propagateToX(trackPosInformationCopy,
                     vertexPos.X(),
                     bz,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_SIN_PHI,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_STEP,
                     matCorr);
  prop->propagateToX(trackNegInformationCopy,
                     vertexNeg.X(),
                     bz,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_SIN_PHI,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_STEP,
                     matCorr);

  xyz[2] = (trackPosInformationCopy.getZ() * helixNeg.rC + trackNegInformationCopy.getZ() * helixPos.rC) / (helixPos.rC + helixNeg.rC);
}

//_______________________________________________________________________
/// \brief Calculate DCA for tracks using the track helix and the inclination angl tan(lambda)
/// \param trk track parameterization to obtain helix parameters
/// \param vtx primary vertex position
/// \param magField magnetic field strenght of L3
/// \return DCAxy, DCAz
template <typename TrackPrecision = float>
std::array<float, 2> CalculateDCAFast(const o2::track::TrackParametrizationWithError<TrackPrecision>& trk, const o2::math_utils::Point3D<float>& vtx, const float magField)
{

  std::array<float, 2> dca{};

  // obtain circle from track in x-y plane
  const o2::track::TrackAuxPar helixPos(trk, magField);

  // obtain position in global coordinates and tan(lambda)
  const auto posTrack = trk.getXYZGlo();
  const float trX = posTrack.X();
  const float trY = posTrack.Y();
  const float trZ = posTrack.Z();
  const float tangentLambda = trk.getTgl();

  // Calculate DCAxy
  // Use distance in x and y between circle center and vtx, afterwards subtract radius of circle
  const float distX = helixPos.xC - vtx.X();
  const float distY = helixPos.yC - vtx.Y();
  const float trackCircCenter = std::sqrt(distX * distX + distY * distY);
  dca[0] = helixPos.rC - trackCircCenter;

  // Calculate DCAz
  // First step: Calculate arc lenght of circle between current position and primary vertex in x-y
  const float theta0 = std::atan2(trY - helixPos.yC, trX - helixPos.xC);
  const float thetav = std::atan2(vtx.Y() - helixPos.yC, vtx.X() - helixPos.xC);

  // Make sure angle is between -pi and pi
  const auto dtheta = RecoDecay::constrainAngle<float>(thetav - theta0, -o2::constants::math::PI);

  // arc-lenght along helix
  const float arcLenght = std::fabs(helixPos.rC * dtheta);

  // get global z-position at DCA
  const float ZPosGlo = trZ - arcLenght * tangentLambda;

  // DCA calculated from
  dca[1] = ZPosGlo - vtx.Z();

  return dca;
}

//_______________________________________________________________________
/// \brief Calculate the track momentum at a different place of the track Helix.
/// \param s arc-lenght where track should be propagated to
/// \param track particle track
/// \param trHelix track helix param. for circle approximation
/// \param bz magnetic field strenght of L3
/// \param addPhi optional additional rotation of the track in phi-direction
/// \return track momentum vector at new position
template <o2::soa::is_iterator TTrack>
inline std::array<float, 3> getPropMomentumFromTrackHelix(const float s, const TTrack& track, const o2::track::TrackAuxPar& trHelix, const float bz, float addPhi = 0.f)
{

  // Calculate the change in phi considering the track radius and the track arc lenght s
  const float dphi = -track.sign() * 0.3f * bz * s / 100. / track.pt(); // s in cm
  const float dotProd = std::cos(track.phi() + addPhi) * track.pt() * trHelix.xC + std::sin(track.phi() + addPhi) * track.pt() * trHelix.yC;
  if (dotProd < 0) {
    addPhi -= o2::constants::math::PI;
  }

  // Calculate the phi at the secondary vertex
  const auto phi = RecoDecay::constrainAngle<float>(track.phi() + dphi + addPhi);

  // Calculate px,y,z at the new propagated vertex
  std::array<float, 3> trackP{};
  trackP[0] = std::cos(phi) * track.pt();
  trackP[1] = std::sin(phi) * track.pt();
  trackP[2] = track.tgl() * track.pt();

  return trackP;
}

//_______________________________________________________________________
template <typename TrackPrecision = float, o2::soa::is_iterator T1, o2::soa::is_iterator T2>
inline void Vtx_recalculation(o2::base::Propagator* prop, T1 lTrackPos, T2 lTrackNeg, std::array<float, 3>& xyz, o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE)
{

  // o2::track::TrackParametrizationWithError<TrackPrecision> = TrackParCov, I use the full version to have control over the data type
  o2::track::TrackParametrizationWithError<TrackPrecision> trackPosInformation = getTrackParCov(lTrackPos); // first get an object that stores Track information (positive)
  o2::track::TrackParametrizationWithError<TrackPrecision> trackNegInformation = getTrackParCov(lTrackNeg); // first get an object that stores Track information (negative)

  Vtx_recalculationParCov<TrackPrecision>(prop, trackPosInformation, trackNegInformation, xyz, matCorr);
}

//_______________________________________________________________________
/// \brief Function to calculate a score based on cosPA and PCA. The smaller the score the better the V0 Candidate
/// \param cosPA cosine of pointing angle
/// \param pca point of closest approach in cm
/// \param weight how much is cosPA weighted (for pca its 1-weight). Weight has to be between 0 and 1
/// \return final score
inline float getScoreV0(float cosPA, float pca, float weight)
{
  float cosScore = 60 * std::acos(cosPA); // pointing angle in degrees, the smaller the better
  float pcaScore = pca / 3.f;             // assume pca is between 0 and 3
  float wCos = weight;
  float wPca = 1.f - wCos; // random values for now
  float score = wCos * cosScore + wPca * pcaScore;
  return score;
}

#endif // PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_
