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
#include "Framework/AnalysisTask.h"

//_______________________________________________________________________
bool checkAP(const float alpha, const float qt, const float alpha_max = 0.95, const float qt_max = 0.05)
{
  float ellipse = pow(alpha / alpha_max, 2) + pow(qt / qt_max, 2);
  if (ellipse < 1.0) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
template <typename TrackPrecision = float, typename T>
void Vtx_recalculation(o2::base::Propagator* prop, T lTrackPos, T lTrackNeg, float xyz[3])
{
  float bz = prop->getNominalBz();

  //*******************************************************

  // o2::track::TrackParametrizationWithError<TrackPrecision> = TrackParCov, I use the full version to have control over the data type
  o2::track::TrackParametrizationWithError<TrackPrecision> trackPosInformation = getTrackParCov(lTrackPos); // first get an object that stores Track information (positive)
  o2::track::TrackParametrizationWithError<TrackPrecision> trackNegInformation = getTrackParCov(lTrackNeg); // first get an object that stores Track information (negative)

  o2::track::TrackAuxPar helixPos(trackPosInformation, bz); // This object is a descendant of a CircleXY and stores cirlce information with respect to the magnetic field. This object uses functions and information of the o2::track::TrackParametrizationWithError<TrackPrecision> object (positive)
  o2::track::TrackAuxPar helixNeg(trackNegInformation, bz); // This object is a descendant of a CircleXY and stores cirlce information with respect to the magnetic field. This object uses functions and information of the o2::track::TrackParametrizationWithError<TrackPrecision> object (negative)

  xyz[0] = (helixPos.xC * helixNeg.rC + helixNeg.xC * helixPos.rC) / (helixPos.rC + helixNeg.rC); // This calculates the coordinates of the conversion point as an weighted average of the two helix centers. xC and yC should be the global coordinates for the helix center as far as I understand. But you can double check the code of trackPosInformation.getCircleParamsLoc
  xyz[1] = (helixPos.yC * helixNeg.rC + helixNeg.yC * helixPos.rC) / (helixPos.rC + helixNeg.rC); // If this calculation doesn't work check if the rotateZ function, because the "documentation" says I get global coordinates but maybe i don't.

  // I am unsure about the Z calculation but this is how it is done in AliPhysics as far as I understand
  o2::track::TrackParametrizationWithError<TrackPrecision> trackPosInformationCopy = o2::track::TrackParametrizationWithError<TrackPrecision>(trackPosInformation);
  o2::track::TrackParametrizationWithError<TrackPrecision> trackNegInformationCopy = o2::track::TrackParametrizationWithError<TrackPrecision>(trackNegInformation);

  // I think this calculation gets the closest point on the track to the conversion point
  // This alpha is a different alpha than the usual alpha and I think it is the angle between X axis and conversion point
  Double_t alphaPos = TMath::Pi() + TMath::ATan2(-(xyz[1] - helixPos.yC), -(xyz[0] - helixPos.xC));
  Double_t alphaNeg = TMath::Pi() + TMath::ATan2(-(xyz[1] - helixNeg.yC), -(xyz[0] - helixNeg.xC));

  Double_t vertexXPos = helixPos.xC + helixPos.rC * TMath::Cos(alphaPos);
  Double_t vertexYPos = helixPos.yC + helixPos.rC * TMath::Sin(alphaPos);
  Double_t vertexXNeg = helixNeg.xC + helixNeg.rC * TMath::Cos(alphaNeg);
  Double_t vertexYNeg = helixNeg.yC + helixNeg.rC * TMath::Sin(alphaNeg);

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
                     o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrNONE);
  // o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrLUT);
  prop->propagateToX(trackNegInformationCopy,
                     vertexNegRot.X(),
                     bz,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_SIN_PHI,
                     o2::base::PropagatorImpl<TrackPrecision>::MAX_STEP,
                     o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrNONE);
  // o2::base::PropagatorImpl<TrackPrecision>::MatCorrType::USEMatCorrLUT);

  // TODO: This is still off and needs to be checked...
  xyz[2] = (trackPosInformationCopy.getZ() * helixNeg.rC + trackNegInformationCopy.getZ() * helixPos.rC) / (helixPos.rC + helixNeg.rC);
}
//_______________________________________________________________________
float getPhivPair(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg, int cpos, int cneg, float bz)
{
  // cos(phiv) = w*a /|w||a|
  // with w = u x v
  // and  a = u x z / |u x z|   , unit vector perpendicular to v12 and z-direction (magnetic field)
  // u = v12 / |v12|            , the unit vector of v12
  // v = v1 x v2 / |v1 x v2|    , unit vector perpendicular to v1 and v2

  // momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
  // vector product of pep X pem
  // std::array<float, 3> arr_pos{t1.GetPx(), t1.GetPy(), t1.GetPz()};
  // std::array<float, 3> arr_ele{t2.GetPx(), t2.GetPy(), t2.GetPz()};
  std::array<float, 3> arr_pos{pxpos, pypos, pzpos};
  std::array<float, 3> arr_ele{pxneg, pyneg, pzneg};
  std::array<double, 3> pos_x_ele{0, 0, 0};
  // LOGF(info, "Q1 = %d , Q2 = %d", cpos, cneg);

  if (cpos * cneg > 0) { // Like Sign
    if (bz < 0) {
      if (cpos > 0) {
        pos_x_ele = RecoDecay::crossProd(arr_pos, arr_ele);
      } else {
        pos_x_ele = RecoDecay::crossProd(arr_ele, arr_pos);
      }
    } else {
      if (cpos > 0) {
        pos_x_ele = RecoDecay::crossProd(arr_ele, arr_pos);
      } else {
        pos_x_ele = RecoDecay::crossProd(arr_pos, arr_ele);
      }
    }
  } else { // Unlike Sign
    if (bz > 0) {
      if (cpos > 0) {
        pos_x_ele = RecoDecay::crossProd(arr_pos, arr_ele);
      } else {
        pos_x_ele = RecoDecay::crossProd(arr_ele, arr_pos);
      }
    } else {
      if (cpos > 0) {
        pos_x_ele = RecoDecay::crossProd(arr_ele, arr_pos);
      } else {
        pos_x_ele = RecoDecay::crossProd(arr_pos, arr_ele);
      }
    }
  }

  // unit vector of pep X pem
  float vx = pos_x_ele[0] / RecoDecay::sqrtSumOfSquares(pos_x_ele[0], pos_x_ele[1], pos_x_ele[2]);
  float vy = pos_x_ele[1] / RecoDecay::sqrtSumOfSquares(pos_x_ele[0], pos_x_ele[1], pos_x_ele[2]);
  float vz = pos_x_ele[2] / RecoDecay::sqrtSumOfSquares(pos_x_ele[0], pos_x_ele[1], pos_x_ele[2]);

  // unit vector of (pep+pem)
  float ux = (pxpos + pxneg) / RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  float uy = (pypos + pyneg) / RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  float uz = (pzpos + pzneg) / RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);

  float ax = uy / TMath::Sqrt(ux * ux + uy * uy);
  float ay = -ux / TMath::Sqrt(ux * ux + uy * uy);

  // The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
  float wx = uy * vz - uz * vy;
  float wy = uz * vx - ux * vz;
  // by construction, (wx,wy,wz) must be a unit vector. Measure angle between (wx,wy,wz) and (ax,ay,0).
  // The angle between them should be small if the pair is conversion. This function then returns values close to pi!
  return TMath::ACos(wx * ax + wy * ay); // phiv in [0,pi] //cosPhiV = wx * ax + wy * ay;
}
//_______________________________________________________________________
float getPsiPair(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
{
  // float pxpos = t1.GetPx();
  // float pypos = t1.GetPy();
  // float pzpos = t1.GetPz();
  // float pxneg = t2.GetPx();
  // float pyneg = t2.GetPy();
  // float pzneg = t2.GetPz();

  auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };
  float ptot2 = RecoDecay::p2(pxpos, pypos, pzpos) * RecoDecay::p2(pxneg, pyneg, pzneg);
  float argcos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}) / std::sqrt(ptot2);
  float thetaPos = std::atan2(RecoDecay::sqrtSumOfSquares(pxpos, pypos), pzpos);
  float thetaNeg = std::atan2(RecoDecay::sqrtSumOfSquares(pxneg, pyneg), pzneg);
  float argsin = (thetaNeg - thetaPos) / std::acos(clipToPM1(argcos));
  return std::asin(clipToPM1(argsin));
}
//_______________________________________________________________________
#endif // PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_
