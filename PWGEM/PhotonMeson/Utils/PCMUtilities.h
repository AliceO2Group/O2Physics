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
#include <Framework/AnalysisDataModel.h>
#include <Framework/Concepts.h>
#include <Framework/DataTypes.h>
#include <ReconstructionDataFormats/HelixHelper.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <Math/GenVector/DisplacementVector2D.h> // IWYU pragma: keep (for rotate)
#include <Math/MatrixFunctions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <Math/Vector2D.h> // IWYU pragma: keep (do not replace with Math/Vector2Dfwd.h)
#include <Math/Vector2Dfwd.h>
#include <TPDGCode.h>

#include <GPUROOTCartesianFwd.h>

#include <array>
#include <cmath>
#include <cstdint>

//_______________________________________________________________________
/// \brief Armenteros-Podolanski cut of photon conversions: (alpha, qT) inside the ellipse with half axes alpha_max and qt_max
/// \param alpha Armenteros alpha
/// \param qt Armenteros qT (GeV/c)
/// \param alpha_max half axis in alpha
/// \param qt_max half axis in qT (GeV/c)
/// \return true if inside the ellipse
inline bool checkAP(const float alpha, const float qt, const float alpha_max = 0.95, const float qt_max = 0.05)
{
  float ellipse = std::pow(alpha / alpha_max, 2) + std::pow(qt / qt_max, 2);
  return (ellipse < 1.0);
}
//_______________________________________________________________________
/// \brief Armenteros alpha of a V0: (pL+ - pL-) / (pL+ + pL-), pL = momentum of a leg along the V0 momentum
/// \param pxpos px of the positive leg (GeV/c), same for pypos, pzpos
/// \param pxneg px of the negative leg (GeV/c), same for pyneg, pzneg
/// \return alpha
inline float v0Alpha(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
{
  float momTot = RecoDecay::p(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  float lQlNeg = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
  float lQlPos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
  return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // longitudinal momentum asymmetry of v0
}
//_______________________________________________________________________
/// \brief Armenteros qT of a V0: momentum of the negative leg transverse to the V0 momentum
/// \param pxpos px of the positive leg (GeV/c), same for pypos, pzpos
/// \param pxneg px of the negative leg (GeV/c), same for pyneg, pzneg
/// \return qT (GeV/c)
inline float v0Qt(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
{
  float momTot = RecoDecay::p2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  float dp = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
  return std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qt of v0
}
//_______________________________________________________________________
/// \brief conversion point of two legs (Run 2 recipe): xy = touching point of the two helix circles (centres weighted with the
/// radius of the other circle), z = radius-weighted mean of the leg z after propagating each leg to its point closest to it
/// \param prop propagator (also gives the nominal field)
/// \param trackPosInformation positive leg
/// \param trackNegInformation negative leg
/// \param xyz output: conversion point (cm)
/// \param matCorr material correction of the propagation
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
/// \brief conversion point of two legs from their track rows (Vtx_recalculationParCov)
/// \param prop propagator (also gives the nominal field)
/// \param lTrackPos positive leg
/// \param lTrackNeg negative leg
/// \param xyz output: conversion point (cm)
/// \param matCorr material correction of the propagation
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
//_______________________________________________________________________
/// \brief truth class of a V0 candidate, from the mothers of its two legs
enum V0TruthClass {
  kV0True = 0,         // both legs come from the SAME true photon
  kV0CrossLegFake = 1, // legs come from TWO different true photons
  kV0OtherFake = 2     // at least one leg is not a photon conversion leg
};
//_______________________________________________________________________
/// \brief global index of the photon a track descends from
/// \param track track with an MC label
/// \param mcparticles the full McParticles table
/// \return index of the mother if it is a photon, -1 otherwise
template <o2::soa::is_iterator TTrack, o2::soa::is_table TMCParticles>
inline int64_t getPhotonMotherId(TTrack const& track, TMCParticles const& mcparticles)
{
  if (!track.has_mcParticle()) {
    return -1;
  }
  const auto mcp = mcparticles.iteratorAt(track.mcParticleId());
  if (!mcp.has_mothers()) {
    return -1;
  }
  const auto& mothers = mcp.mothersIds();
  if (mothers.empty() || mothers[0] < 0) {
    return -1;
  }
  return (mcparticles.iteratorAt(mothers[0]).pdgCode() == PDG_t::kGamma) ? static_cast<int64_t>(mothers[0]) : -1;
}

//_______________________________________________________________________
/// \brief classify a V0 candidate by the photon mothers of its two legs
/// \param motherPos photon mother of the positive leg, -1 if none
/// \param motherEle photon mother of the negative leg, -1 if none
/// \return kV0True, kV0CrossLegFake or kV0OtherFake
template <o2::soa::is_iterator TTrack, o2::soa::is_table TMCParticles>
inline V0TruthClass classifyV0Truth(TTrack const& pos, TTrack const& ele, TMCParticles const& mcparticles,
                                    int64_t& motherPos, int64_t& motherEle)
{
  motherPos = getPhotonMotherId(pos, mcparticles);
  motherEle = getPhotonMotherId(ele, mcparticles);
  if (motherPos >= 0 && motherPos == motherEle) {
    return kV0True;
  }
  if (motherPos >= 0 && motherEle >= 0) {
    return kV0CrossLegFake;
  }
  return kV0OtherFake;
}

//_______________________________________________________________________
/// \brief electron-PID compatibility of the two legs: mean |n sigma_e| + 2 * (n sigma_e+ + n sigma_e-) in the TPC
/// \param pos positive leg
/// \param ele negative leg
/// \return kappa
template <o2::soa::is_iterator TLeg>
float getV0Kappa(TLeg const& pos, TLeg const& ele)
{
  float kappa = (std::fabs(pos.tpcNSigmaEl()) + std::fabs(ele.tpcNSigmaEl())) / 2.f + (2.f * (pos.tpcNSigmaEl() + ele.tpcNSigmaEl()));
  return kappa;
}

//_______________________________________________________________________
/// \brief gap between the helix circles of two legs in xy
/// \param cPos circle of the positive leg
/// \param cEle circle of the negative leg
/// \return distance of the centres minus the sum of the radii (cm); < 0: the circles overlap
inline float circleGap(o2::track::TrackAuxPar const& cPos, o2::track::TrackAuxPar const& cEle)
{
  return std::hypot(cEle.xC - cPos.xC, cEle.yC - cPos.yC) - (cPos.rC + cEle.rC);
}

/// \brief return value of tpcDriftSide if the side cannot be determined
constexpr int TPCSideUndetermined = 2;

//_______________________________________________________________________
/// \brief TPC side in which a track drifted, from the TPC side flags, else from z, tgl and the time asymmetry
/// \param track track row (flags, z, tgl, tpcDeltaTFwd, tpcDeltaTBwd)
/// \return 1 A side, -1 C side, 0 both sides (crosses the central electrode), TPCSideUndetermined
template <o2::soa::is_iterator TTrack>
inline int tpcDriftSide(TTrack const& track)
{
  constexpr float MaxCECrossingDeltaTNS = 1000.f;
  const auto sides = track.flags() & (o2::aod::track::TPCSideA | o2::aod::track::TPCSideC);
  if (sides == o2::aod::track::TPCSideA) {
    return 1;
  }
  if (sides == o2::aod::track::TPCSideC) {
    return -1;
  }
  if (sides != 0) {
    return 0;
  }
  const float dFwd = track.tpcDeltaTFwd();
  const float dBwd = track.tpcDeltaTBwd();
  if (dFwd == dBwd && dFwd < MaxCECrossingDeltaTNS) {
    return 0;
  }
  const bool zPositive = track.z() > 0.f;
  const bool tglPositive = track.tgl() > 0.f;
  if (zPositive == tglPositive) {
    return tglPositive ? 1 : -1;
  }
  if (dBwd == 0.f && dFwd > 0.f) {
    return zPositive ? 1 : -1;
  }
  if (dBwd > 0.f) {
    return tglPositive ? 1 : -1;
  }
  return TPCSideUndetermined;
}

/// \brief z and momentum of a track helix at a point of its circle (helixPointAt)
struct HelixPoint {
  bool ok = false;
  bool inward = false;
  float z = 0.f;
  std::array<float, 3> p{};
};

//_______________________________________________________________________
/// \brief z and momentum of a track helix at an xy point of its circle
/// \param par track parameters
/// \param xP x of the point on the circle (cm)
/// \param yP y of the point on the circle (cm)
/// \param xB x of the beam position (cm)
/// \param yB y of the beam position (cm)
/// \param bz magnetic field (kG)
/// \return ok = false for a degenerate circle; inward = the track moves towards the beam at the point (more than 10 cm from
/// the beam), then the momentum is reversed; z and momentum at the point
inline HelixPoint helixPointAt(o2::track::TrackPar const& par, const float xP, const float yP, const float xB, const float yB, const float bz)
{
  HelixPoint h;
  const o2::track::TrackAuxPar c(par, bz);
  std::array<float, 3> r{}, p{};
  par.getXYZGlo(r);
  par.getPxPyPzGlo(p);
  const float pt = std::hypot(p[0], p[1]);
  if (!(c.rC > 0.f) || !(pt > 0.f)) {
    return h;
  }
  const float sense = (r[0] - c.xC) * p[1] - (r[1] - c.yC) * p[0] > 0.f ? 1.f : -1.f;
  float back = RecoDecay::constrainAngle(sense * (std::atan2(r[1] - c.yC, r[0] - c.xC) - std::atan2(yP - c.yC, xP - c.xC)), 0.f);
  const float tx = -sense * (yP - c.yC) / c.rC;
  const float ty = sense * (xP - c.xC) / c.rC;
  constexpr float MinRadiusOrientation = 10.f;
  h.inward = std::hypot(xP - xB, yP - yB) > MinRadiusOrientation && tx * (xP - xB) + ty * (yP - yB) < 0.f;
  if (!h.inward) {
    if (back > o2::constants::math::PI && (r[0] - xB) * p[0] + (r[1] - yB) * p[1] > 0.f) {
      back -= o2::constants::math::TwoPI; // o2-linter: disable=two-pi-add-subtract (signed arc length to a point behind the track, not an angle range)
    }
    h.z = r[2] - par.getTgl() * c.rC * back;
    h.p = {pt * tx, pt * ty, pt * par.getTgl()};
  } else {
    h.z = r[2] + par.getTgl() * c.rC * (back > 0.f ? o2::constants::math::TwoPI - back : 0.f);
    h.p = {-pt * tx, -pt * ty, -pt * par.getTgl()};
  }
  h.ok = true;
  return h;
}

//_______________________________________________________________________
/// \brief compatibility of two track pieces with one track
/// \param trA first track
/// \param trB second track
/// \param bz magnetic field (kG)
/// \return chi2 of (y, snp, tgl, q/pt) at a common x with the summed covariances (4 dof); -1 if the propagation or the inversion fails
template <o2::soa::is_iterator TTrack>
inline float splitPieceChi2(TTrack const& trA, TTrack const& trB, const float bz)
{
  auto a = getTrackParCov(trA);
  auto b = getTrackParCov(trB);
  if (!a.rotate(b.getAlpha()) || !a.propagateTo(b.getX(), bz)) {
    a = getTrackParCov(trA);
    b = getTrackParCov(trB);
    if (!b.rotate(a.getAlpha()) || !b.propagateTo(a.getX(), bz)) {
      return -1.f;
    }
  }
  ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepSym<double, 4>> c;
  c(0, 0) = a.getSigmaY2() + b.getSigmaY2();
  c(1, 0) = a.getSigmaSnpY() + b.getSigmaSnpY();
  c(1, 1) = a.getSigmaSnp2() + b.getSigmaSnp2();
  c(2, 0) = a.getSigmaTglY() + b.getSigmaTglY();
  c(2, 1) = a.getSigmaTglSnp() + b.getSigmaTglSnp();
  c(2, 2) = a.getSigmaTgl2() + b.getSigmaTgl2();
  c(3, 0) = a.getSigma1PtY() + b.getSigma1PtY();
  c(3, 1) = a.getSigma1PtSnp() + b.getSigma1PtSnp();
  c(3, 2) = a.getSigma1PtTgl() + b.getSigma1PtTgl();
  c(3, 3) = a.getSigma1Pt2() + b.getSigma1Pt2();
  if (!c.Invert()) {
    return -1.f;
  }
  const ROOT::Math::SVector<double, 4> d(a.getY() - b.getY(), a.getSnp() - b.getSnp(), a.getTgl() - b.getTgl(), a.getQ2Pt() - b.getQ2Pt());
  return static_cast<float>(ROOT::Math::Similarity(c, d));
}

#endif // PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_
