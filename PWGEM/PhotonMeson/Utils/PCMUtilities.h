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

/// \commonly used for PCM analyses.
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_

#include <TVector2.h>
#include "DCAFitter/HelixHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"

//_______________________________________________________________________
inline bool checkAP(const float alpha, const float qt, const float alpha_max = 0.95, const float qt_max = 0.05)
{
  float ellipse = pow(alpha / alpha_max, 2) + pow(qt / qt_max, 2);
  if (ellipse < 1.0) {
    return true;
  } else {
    return false;
  }
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
inline void Vtx_recalculationParCov(o2::base::Propagator* prop, const o2::track::TrackParametrizationWithError<TrackPrecision>& trackPosInformation, const o2::track::TrackParametrizationWithError<TrackPrecision>& trackNegInformation, float xyz[3], o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE)
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
  float alphaPos = M_PI + std::atan2(-(xyz[1] - helixPos.yC), -(xyz[0] - helixPos.xC));
  float alphaNeg = M_PI + std::atan2(-(xyz[1] - helixNeg.yC), -(xyz[0] - helixNeg.xC));

  float vertexXPos = helixPos.xC + helixPos.rC * std::cos(alphaPos);
  float vertexYPos = helixPos.yC + helixPos.rC * std::sin(alphaPos);
  float vertexXNeg = helixNeg.xC + helixNeg.rC * std::cos(alphaNeg);
  float vertexYNeg = helixNeg.yC + helixNeg.rC * std::sin(alphaNeg);

  TVector2 vertexPos(vertexXPos, vertexYPos);
  TVector2 vertexNeg(vertexXNeg, vertexYNeg);

  // Convert to local coordinate system
  TVector2 vertexPosRot = vertexPos.Rotate(-trackPosInformationCopy.getAlpha());
  TVector2 vertexNegRot = vertexNeg.Rotate(-trackNegInformationCopy.getAlpha());

  prop->propagateToX(trackPosInformationCopy,
                     vertexPosRot.X(),
                     bz,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_SIN_PHI,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_STEP,
                     matCorr);
  prop->propagateToX(trackNegInformationCopy,
                     vertexNegRot.X(),
                     bz,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_SIN_PHI,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_STEP,
                     matCorr);

  xyz[2] = (trackPosInformationCopy.getZ() * helixNeg.rC + trackNegInformationCopy.getZ() * helixPos.rC) / (helixPos.rC + helixNeg.rC);
}
//_______________________________________________________________________
template <typename TrackPrecision = float, typename T1, typename T2>
inline void Vtx_recalculation(o2::base::Propagator* prop, T1 lTrackPos, T2 lTrackNeg, float xyz[3], o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE)
{
  // o2::track::TrackParametrizationWithError<TrackPrecision> = TrackParCov, I use the full version to have control over the data type
  o2::track::TrackParametrizationWithError<TrackPrecision> trackPosInformation = getTrackParCov(lTrackPos); // first get an object that stores Track information (positive)
  o2::track::TrackParametrizationWithError<TrackPrecision> trackNegInformation = getTrackParCov(lTrackNeg); // first get an object that stores Track information (negative)

  Vtx_recalculationParCov<TrackPrecision>(prop, trackPosInformation, trackNegInformation, xyz, matCorr);
}
#endif // PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_
