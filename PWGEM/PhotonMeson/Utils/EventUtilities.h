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

/// \file EventUtilities.h
/// \brief Utilities for event level variables (e.g. Q vectors) used in PWGEM/PhotonMeson
/// \author M. Hemmer, marvin.hemmer@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_EVENTUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_EVENTUTILITIES_H_

#include <Framework/Concepts.h>
#include <Framework/Logger.h>

#include <cmath>
#include <span>
#include <utility>
#include <vector>

namespace o2::aod::pwgem::pm::utils::event
{

/// Detectors that provide a Q vector in the EMEventsQvec tables
enum QvecEstimator {
  FT0M = 0,
  FT0A = 1,
  FT0C = 2,
  TPCPos = 3,
  TPCNeg = 4,
  TPCTot = 5,
  FV0A = 6
};

/// Supported flow harmonics
enum Harmonics {
  kNone = 0,
  kDirect = 1,
  kElliptic = 2,
  kTriangluar = 3,
  kQuadrangular = 4,
  kPentagonal = 5,
  kHexagonal = 6,
  kHeptagonal = 7,
  kOctagonal = 8
};

/// Value returned for Q vectors that are not available
constexpr float InvalidQvec = -999.f;

// Compile-time checks for the availability of Q vector columns.
// EMEventsQvec_000 has q2 and q3 without FV0A, EMEventsQvec_001 has q2 and q3 with FV0A,
// EMEventsQvec2 has only q2 (with FV0A) and EMEventsQvec3 has only q3 (with FV0A).
template <typename T>
concept HasQvec2 = requires(T const& c) { c.q2xft0m(); };

template <typename T>
concept HasQvec3 = requires(T const& c) { c.q3xft0m(); };

template <typename T>
concept HasQvec2FV0A = requires(T const& c) { c.q2xfv0a(); };

template <typename T>
concept HasQvec3FV0A = requires(T const& c) { c.q3xfv0a(); };

template <typename T>
concept HasQvec = HasQvec2<T> || HasQvec3<T>;

// Compile-time checks for the availability of the occupancy columns (evsel::NumTracksInTimeRange and evsel::SumAmpFT0CInTimeRange)
template <typename T>
concept HasTrackOccupancy = requires(T const& c) { c.trackOccupancyInTimeRange(); };

template <typename T>
concept HasFT0COccupancy = requires(T const& c) { c.ft0cOccupancyInTimeRange(); };

/// \brief Get the Q vector of a given detector and harmonic
/// \param collision collision with the Q vector information
/// \param detector detector used for the Q vector (see QvecEstimator)
/// \param harmonic harmonic of the Q vector (see Harmonics)
/// \return pair of (Qx, Qy). {InvalidQvec, InvalidQvec} if not available in the table or for an unsupported harmonic
template <o2::soa::is_iterator TCollision>
std::pair<float, float> getQvec(TCollision const& collision, int detector, int harmonic)
{
  if (harmonic == kElliptic) {
    if constexpr (HasQvec2<TCollision>) {
      switch (detector) {
        case QvecEstimator::FT0M:
          return {collision.q2xft0m(), collision.q2yft0m()};
        case QvecEstimator::FT0A:
          return {collision.q2xft0a(), collision.q2yft0a()};
        case QvecEstimator::FT0C:
          return {collision.q2xft0c(), collision.q2yft0c()};
        case QvecEstimator::TPCPos:
          return {collision.q2xbpos(), collision.q2ybpos()};
        case QvecEstimator::TPCNeg:
          return {collision.q2xbneg(), collision.q2ybneg()};
        case QvecEstimator::TPCTot:
          return {collision.q2xbtot(), collision.q2ybtot()};
        case QvecEstimator::FV0A:
          if constexpr (HasQvec2FV0A<TCollision>) {
            return {collision.q2xfv0a(), collision.q2yfv0a()};
          } else {
            return {InvalidQvec, InvalidQvec};
          }
        default:
          LOG(warning) << "Q vector estimator not valid. Falling back to FT0M";
          return {collision.q2xft0m(), collision.q2yft0m()};
      }
    }
  } else if (harmonic == kTriangluar) {
    if constexpr (HasQvec3<TCollision>) {
      switch (detector) {
        case QvecEstimator::FT0M:
          return {collision.q3xft0m(), collision.q3yft0m()};
        case QvecEstimator::FT0A:
          return {collision.q3xft0a(), collision.q3yft0a()};
        case QvecEstimator::FT0C:
          return {collision.q3xft0c(), collision.q3yft0c()};
        case QvecEstimator::TPCPos:
          return {collision.q3xbpos(), collision.q3ybpos()};
        case QvecEstimator::TPCNeg:
          return {collision.q3xbneg(), collision.q3ybneg()};
        case QvecEstimator::TPCTot:
          return {collision.q3xbtot(), collision.q3ybtot()};
        case QvecEstimator::FV0A:
          if constexpr (HasQvec3FV0A<TCollision>) {
            return {collision.q3xfv0a(), collision.q3yfv0a()};
          } else {
            return {InvalidQvec, InvalidQvec};
          }
        default:
          LOG(warning) << "Q vector estimator not valid. Falling back to FT0M";
          return {collision.q3xft0m(), collision.q3yft0m()};
      }
    }
  }
  return {InvalidQvec, InvalidQvec};
}

/// \brief Check if all given Q vector components are within a reasonable range
/// \param qvecs Q vector components (e.g. {Qx_1, Qy_1, Qx_2, Qy_2, ...})
/// \param maxQvec maximum allowed absolute value of each component
inline bool isQvecGood(std::span<const float> qvecs, float maxQvec)
{
  for (const auto& qvec : qvecs) {
    if (std::fabs(qvec) > maxQvec) {
      return false;
    }
  }
  return true;
}

/// \brief Check if the Q vectors of all given detectors are within a reasonable range
/// \param collision collision with the Q vector information
/// \param detectors detectors to check (see QvecEstimator)
/// \param harmonic harmonic of the Q vectors (see Harmonics)
/// \param maxQvec maximum allowed absolute value of Qx and Qy
template <o2::soa::is_iterator TCollision>
bool isQvecGood(TCollision const& collision, std::vector<int> const& detectors, int harmonic, float maxQvec)
{
  for (const auto& detector : detectors) {
    auto [xQvec, yQvec] = getQvec(collision, detector, harmonic);
    if (std::fabs(xQvec) > maxQvec || std::fabs(yQvec) > maxQvec) {
      return false;
    }
  }
  return true;
}

} // namespace o2::aod::pwgem::pm::utils::event

#endif // PWGEM_PHOTONMESON_UTILS_EVENTUTILITIES_H_
