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
///
/// \file   TrackUtilities.cxx
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
/// \brief  Set of utilities for the ALICE3 track handling
/// \since  May 21, 2025
///

#include "TrackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <MathUtils/Primitive2D.h>
#include <MathUtils/Utils.h>
#include <ReconstructionDataFormats/Track.h>

#include <TLorentzVector.h>

#include <array>
#include <cmath>
#include <vector>

void o2::upgrade::convertTLorentzVectorToO2Track(const int charge,
                                                 const TLorentzVector particle,
                                                 const std::vector<double> productionVertex,
                                                 o2::track::TrackParCov& o2track)
{
  std::array<float, 5> params;
  std::array<float, 15> covm = {0.};
  float s, c, x;
  o2::math_utils::sincos(static_cast<float>(particle.Phi()), s, c);
  o2::math_utils::rotateZInv(static_cast<float>(productionVertex[0]), static_cast<float>(productionVertex[1]), x, params[0], s, c);
  params[1] = static_cast<float>(productionVertex[2]);
  params[2] = 0.; // since alpha = phi
  const auto theta = 2. * std::atan(std::exp(-particle.PseudoRapidity()));
  params[3] = 1. / std::tan(theta);
  params[4] = charge / particle.Pt();

  // Initialize TrackParCov in-place
  new (&o2track)(o2::track::TrackParCov)(x, particle.Phi(), params, covm);
}

float o2::upgrade::computeParticleVelocity(float momentum, float mass)
{
  const float a = momentum / mass;
  // uses light speed in cm/ps so output is in those units
  return o2::constants::physics::LightSpeedCm2PS * a / std::sqrt((1.f + a * a));
}

float o2::upgrade::computeTrackLength(o2::track::TrackParCov track, float radius, float magneticField)
{
  // don't make use of the track parametrization
  float length = -100;

  o2::math_utils::CircleXYf_t trcCircle;
  float sna, csa;
  track.getCircleParams(magneticField, trcCircle, sna, csa);

  // distance between circle centers (one circle is at origin -> easy)
  const float centerDistance = std::hypot(trcCircle.xC, trcCircle.yC);

  // condition of circles touching - if not satisfied returned length will be -100
  if (centerDistance < trcCircle.rC + radius && centerDistance > std::fabs(trcCircle.rC - radius)) {
    length = 0.0f;

    // base radical direction
    const float ux = trcCircle.xC / centerDistance;
    const float uy = trcCircle.yC / centerDistance;
    // calculate perpendicular vector (normalized) for +/- displacement
    const float vx = -uy;
    const float vy = +ux;
    // calculate coordinate for radical line
    const float radical = (centerDistance * centerDistance - trcCircle.rC * trcCircle.rC + radius * radius) / (2.0f * centerDistance);
    // calculate absolute displacement from center-to-center axis
    const float displace = (0.5f / centerDistance) * std::sqrt(
                                                       (-centerDistance + trcCircle.rC - radius) *
                                                       (-centerDistance - trcCircle.rC + radius) *
                                                       (-centerDistance + trcCircle.rC + radius) *
                                                       (centerDistance + trcCircle.rC + radius));

    // possible intercept points of track and TOF layer in 2D plane
    const float point1[2] = {radical * ux + displace * vx, radical * uy + displace * vy};
    const float point2[2] = {radical * ux - displace * vx, radical * uy - displace * vy};

    // decide on correct intercept point
    std::array<float, 3> mom;
    track.getPxPyPzGlo(mom);
    const float scalarProduct1 = point1[0] * mom[0] + point1[1] * mom[1];
    const float scalarProduct2 = point2[0] * mom[0] + point2[1] * mom[1];

    // get start point
    std::array<float, 3> startPoint;
    track.getXYZGlo(startPoint);

    float cosAngle = -1000, modulus = -1000;

    if (scalarProduct1 > scalarProduct2) {
      modulus = std::hypot(point1[0] - trcCircle.xC, point1[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
      cosAngle = (point1[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point1[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
    } else {
      modulus = std::hypot(point2[0] - trcCircle.xC, point2[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
      cosAngle = (point2[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point2[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
    }
    cosAngle /= modulus;
    length = trcCircle.rC * std::acos(cosAngle);
    length *= std::sqrt(1.0f + track.getTgl() * track.getTgl());
  }
  return length;
}
