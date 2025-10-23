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

/// \file trackUtilities.h
/// \brief Utilities for manipulating parameters of tracks and vertices
///
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#ifndef COMMON_CORE_TRACKUTILITIES_H_
#define COMMON_CORE_TRACKUTILITIES_H_

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <ReconstructionDataFormats/TrackParametrization.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <GPUROOTCartesianFwd.h>

#include <array>
#include <utility> // std::move

/// Extracts track parameters from a track.
template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrization<TrackPrecision> getTrackPar(const T& track)
{
  std::array<TrackPrecision, o2::track::kNParams> arraypar = {track.y(), track.z(), track.snp(),
                                                              track.tgl(), track.signed1Pt()};
  return o2::track::TrackParametrization<TrackPrecision>(track.x(), track.alpha(), std::move(arraypar));
}

/// Extracts track parameters from a track and sets a TrackParametrization object.
template <typename TrackPrecision = float, typename T>
void setTrackPar(const T& track, o2::track::TrackParametrization<TrackPrecision>& trackPar)
{
  std::array<TrackPrecision, o2::track::kNParams> arraypar = {track.y(), track.z(), track.snp(),
                                                              track.tgl(), track.signed1Pt()};
  trackPar.set(track.x(), track.alpha(), std::move(arraypar));
}

/// Extracts track parameters and covariance matrix from a track.
template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrizationWithError<TrackPrecision> getTrackParCov(const T& track)
{
  std::array<TrackPrecision, o2::track::kNParams> arraypar = {track.y(), track.z(), track.snp(),
                                                              track.tgl(), track.signed1Pt()};
  std::array<TrackPrecision, o2::track::kCovMatSize> covpar = {track.cYY(), track.cZY(), track.cZZ(),
                                                               track.cSnpY(), track.cSnpZ(),
                                                               track.cSnpSnp(), track.cTglY(), track.cTglZ(),
                                                               track.cTglSnp(), track.cTglTgl(),
                                                               track.c1PtY(), track.c1PtZ(), track.c1PtSnp(),
                                                               track.c1PtTgl(), track.c1Pt21Pt2()};
  return o2::track::TrackParametrizationWithError<TrackPrecision>(track.x(), track.alpha(), std::move(arraypar), std::move(covpar));
}

/// Extracts track parameters and covariance matrix from a track and sets a TrackParametrizationWithError object.
template <typename TrackPrecision = float, typename T>
void setTrackParCov(const T& track, o2::track::TrackParametrizationWithError<TrackPrecision>& trackParCov)
{
  std::array<TrackPrecision, o2::track::kNParams> arraypar = {track.y(), track.z(), track.snp(),
                                                              track.tgl(), track.signed1Pt()};
  std::array<TrackPrecision, o2::track::kCovMatSize> covpar = {track.cYY(), track.cZY(), track.cZZ(),
                                                               track.cSnpY(), track.cSnpZ(),
                                                               track.cSnpSnp(), track.cTglY(), track.cTglZ(),
                                                               track.cTglSnp(), track.cTglTgl(),
                                                               track.c1PtY(), track.c1PtZ(), track.c1PtSnp(),
                                                               track.c1PtTgl(), track.c1Pt21Pt2()};
  trackParCov.set(track.x(), track.alpha(), std::move(arraypar), std::move(covpar));
}

/// Extracts primary vertex position and covariance matrix from a collision.
template <typename T>
o2::dataformats::VertexBase getPrimaryVertex(const T& collision)
{
  o2::math_utils::Point3D<float> vtxXYZ(collision.posX(), collision.posY(), collision.posZ());
  std::array<float, 6> vtxCov{collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};
  return o2::dataformats::VertexBase{std::move(vtxXYZ), std::move(vtxCov)};
}

/// Calculates direction of one point w.r.t another point.
/// \param point1,point2  points with {x, y, z} coordinates accessible by index
/// \param phi  azimuth angle in the {x, y} plane (taken w.r.t. the x-axis towards the y-axis)
/// \param theta  angle w.r.t the {x, y} plane towards the z-axis
/// \return phi,theta of point2 w.r.t. point1
template <typename T, typename U, typename V, typename W>
void getPointDirection(const T& point1, const U& point2, V& phi, W& theta)
{
  phi = std::atan2(point2[1] - point1[1], point2[0] - point1[0]);
  // auto x1 = point1[0] * std::cos(phi) + point1[1] * std::sin(phi);
  // auto x2 = point2[0] * std::cos(phi) + point2[1] * std::sin(phi);
  // theta = std::atan2(point2[2] - point1[2], x2 - x1);
  theta = std::atan2(point2[2] - point1[2], RecoDecay::distanceXY(point1, point2));
}

/// Calculates the XX element of a XYZ covariance matrix after rotation of the coordinate system
/// by phi around the z-axis and by minus theta around the new y-axis.
/// \param matrix  matrix
/// \param phi  azimuth angle in the {x, y} plane (taken w.r.t. the x-axis towards the y-axis)
/// \param theta  angle w.r.t the {x, y} plane towards the z-axis
/// \return XX element of the rotated covariance matrix
template <typename T, typename U, typename V>
auto getRotatedCovMatrixXX(const T& matrix, U phi, V theta)
{
  auto cp = std::cos(phi);
  auto sp = std::sin(phi);
  auto ct = std::cos(theta);
  auto st = std::sin(theta);
  return matrix[0] * cp * cp * ct * ct        // covXX
         + matrix[1] * 2. * cp * sp * ct * ct // covXY
         + matrix[2] * sp * sp * ct * ct      // covYY
         + matrix[3] * 2. * cp * ct * st      // covXZ
         + matrix[4] * 2. * sp * ct * st      // covYZ
         + matrix[5] * st * st;               // covZZ
}

/// Method to compute px, py, pz from the track parameters
/// \param trackPars are the track parameters
/// \param pVec is the momentum array to be filled
template <typename T, typename U>
void getPxPyPz(T const& trackPars, U& pVec)
{
  auto ptTrack = 1.f / std::abs(trackPars.getQ2Pt());
  float cs = cosf(trackPars.getAlpha()), sn = sinf(trackPars.getAlpha());
  auto r = std::sqrt((1.f - trackPars.getSnp()) * (1.f + trackPars.getSnp()));
  pVec[0] = ptTrack * (r * cs - trackPars.getSnp() * sn);
  pVec[1] = ptTrack * (trackPars.getSnp() * cs + r * sn);
  pVec[2] = ptTrack * trackPars.getTgl();
}

/// Calculates DCA XYZ of a track w.r.t. the primary vertex and its uncertainty if required.
/// \param track  track from a table containing `o2::aod::TracksCov, o2::aod::TracksDCA`.
/// \param sigmaDcaXYZ2  pointer to the sigma^2 of the impact parameter in XYZ to be calculated
/// \return impact parameter in XYZ of the track w.r.t the primary vertex
template <typename T>
float getDcaXYZ(T const& track, float* sigmaDcaXYZ2 = nullptr)
{
  float dcaXY = track.dcaXY();
  float dcaZ = track.dcaZ();
  float dcaXYZ = std::sqrt(dcaXY * dcaXY + dcaZ * dcaZ);
  if (sigmaDcaXYZ2) {
    if (dcaXYZ < o2::constants::math::Almost0) {
      *sigmaDcaXYZ2 = o2::constants::math::VeryBig; // Protection against division by zero
    } else {
      float dFdxy = 2.f * dcaXY / dcaXYZ;
      float dFdz = 2.f * dcaZ / dcaXYZ;
      *sigmaDcaXYZ2 = track.cYY() * dFdxy * dFdxy + track.cZZ() * dFdz * dFdz + 2.f * track.cZY() * dFdxy * dFdz;
    }
  }
  return dcaXYZ;
}

#endif // COMMON_CORE_TRACKUTILITIES_H_
