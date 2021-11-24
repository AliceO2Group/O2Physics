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
/// \file   PIDTOF.h
/// \author Nicolo' Jacazio
/// \since  02/07/2020
/// \brief  Implementation of the TOF detector response for PID
///

#ifndef O2_FRAMEWORK_PIDTOF_H_
#define O2_FRAMEWORK_PIDTOF_H_

// ROOT includes
#include "Rtypes.h"
#include "TMath.h"

// O2 includes
#include "Framework/Logger.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/PID/DetectorResponse.h"

namespace o2::pid::tof
{

// Utility values
static constexpr float kCSPEED = TMath::C() * 1.0e2f * 1.0e-12f; /// Speed of light in TOF units (cm/ps)

/// \brief Class to handle the the TOF detector response for the TOF beta measurement
template <typename TrackType, o2::track::PID::ID id>
class Beta
{
 public:
  Beta() = default;
  ~Beta() = default;

  /// Computes the beta of a track given a length, a time measurement and an event time (in ps)
  /// \param length Length in cm of the track
  /// \param tofSignal TOF signal in ps for the track
  /// \param collisionTime collision time in ps for the event of the track
  static float GetBeta(const float length, const float tofSignal, const float collisionTime);

  /// Gets the beta for the track of interest
  float GetBeta(const TrackType& trk) const { return GetBeta(trk.length(), trk.tofSignal(), trk.collision().collisionTime() * 1000.f); }

  /// Computes the expected uncertainty on the beta measurement
  static float GetExpectedSigma(const float& length, const float& tofSignal, const float& collisionTime, const float& time_reso);

  /// Gets the expected uncertainty on the beta measurement of the track of interest
  float GetExpectedSigma(const TrackType& trk) const { return GetExpectedSigma(trk.length(), trk.tofSignal(), trk.collision().collisionTime() * 1000.f, mExpectedResolution); }

  /// Gets the expected beta for a given mass hypothesis (no energy loss taken into account)
  static float GetExpectedSignal(const float& mom, const float& mass);

  /// Gets the expected beta given the particle index (no energy loss taken into account)
  float GetExpectedSignal(const TrackType& trk) const { return GetExpectedSignal(trk.p(), o2::track::PID::getMass2Z(id)); }

  /// Gets the number of sigmas with respect the approximate beta (no energy loss taken into account)
  float GetSeparation(const TrackType& trk) const { return (GetBeta(trk) - GetExpectedSignal(trk)) / GetExpectedSigma(trk); }

  float mExpectedResolution = 80; /// Expected time resolution
};

//_________________________________________________________________________
template <typename TrackType, o2::track::PID::ID id>
float Beta<TrackType, id>::GetBeta(const float length, const float tofSignal, const float collisionTime)
{
  if (tofSignal <= 0) {
    return -999.f;
  }
  return length / (tofSignal - collisionTime) / kCSPEED;
}

//_________________________________________________________________________
template <typename TrackType, o2::track::PID::ID id>
float Beta<TrackType, id>::GetExpectedSigma(const float& length, const float& tofSignal, const float& collisionTime, const float& time_reso)
{
  if (tofSignal <= 0) {
    return -999.f;
  }
  return GetBeta(length, tofSignal, collisionTime) / (tofSignal - collisionTime) * time_reso;
}

//_________________________________________________________________________
template <typename TrackType, o2::track::PID::ID id>
float Beta<TrackType, id>::GetExpectedSignal(const float& mom, const float& mass)
{
  if (mom > 0) {
    return mom / TMath::Sqrt(mom * mom + mass * mass);
  }
  return 0;
}

/// \brief Class to handle the the TOF detector response for the expected time
template <typename TrackType, o2::track::PID::ID id>
class ExpTimes
{
 public:
  ExpTimes() = default;
  ~ExpTimes() = default;

  /// Computes the expected time of a track, given it TOF expected momentum
  static float ComputeExpectedTime(const float& tofExpMom, const float& length, const float& massZ)
  {
    if (tofExpMom <= 0.f) {
      return 0.f;
    }
    const float energy = sqrt((massZ * massZ) + (tofExpMom * tofExpMom));
    const float exp = length * energy / (kCSPEED * tofExpMom);
    return exp >= 0.f ? exp : 0.f;
  }

  /// Gets the expected signal of the track of interest under the PID assumption
  static float GetExpectedSignal(const TrackType& trk) { return ComputeExpectedTime(trk.tofExpMom() / kCSPEED, trk.length(), o2::track::PID::getMass2Z(id)); }

  /// Gets the expected resolution of the measurement
  static float GetExpectedSigma(const DetectorResponse& response, const TrackType& trk)
  {
    if (!trk.hasTOF()) {
      return -999.f;
    }
    if (!trk.has_collision()) {
      return -999.f;
    }
    const float x[7] = {trk.p(), trk.tofSignal(), trk.collision().collisionTimeRes() * 1000.f, o2::track::PID::getMass2Z(id), trk.length(), trk.sigma1Pt(), trk.pt()};
    const float reso = response(response.kSigma, x);
    return reso >= 0.f ? reso : 0.f;
  }

  /// Gets the number of sigmas with respect the expected time
  static float GetSeparation(const DetectorResponse& response, const TrackType& trk) { return trk.has_collision() ? (trk.hasTOF() ? (trk.tofSignal() - trk.collision().collisionTime() * 1000.f - GetExpectedSignal(trk)) / GetExpectedSigma(response, trk) : -999.f) : -999.f; }

  /// Gets the expected resolution of the measurement from the track time
  float GetExpectedSigmaFromTrackTime(const DetectorResponse& response, const TrackType& trk) const;

  /// Gets the number of sigmas with respect the expected time from the track time
  float GetSeparationFromTrackTime(const DetectorResponse& response, const TrackType& trk) const;
};

/// \brief Class to convert the trackTime to the tofSignal used for PID
template <typename TrackType>
class TOFSignal
{
 public:
  TOFSignal() = default;
  ~TOFSignal() = default;

  template <o2::track::PID::ID pid>
  using ResponseImplementation = tof::ExpTimes<TrackType, pid>;
  static constexpr auto responseEl = ResponseImplementation<o2::track::PID::Electron>();
  static constexpr auto responseMu = ResponseImplementation<o2::track::PID::Muon>();
  static constexpr auto responsePi = ResponseImplementation<o2::track::PID::Pion>();
  static constexpr auto responseKa = ResponseImplementation<o2::track::PID::Kaon>();
  static constexpr auto responsePr = ResponseImplementation<o2::track::PID::Proton>();
  static constexpr auto responseDe = ResponseImplementation<o2::track::PID::Deuteron>();
  static constexpr auto responseTr = ResponseImplementation<o2::track::PID::Triton>();
  static constexpr auto responseHe = ResponseImplementation<o2::track::PID::Helium3>();
  static constexpr auto responseAl = ResponseImplementation<o2::track::PID::Alpha>();

  /// Computes the expected time of a track, given its TOF expected time
  /// \param trackTime trackTime in ns
  /// \param response response to use to compute the expected time for the propagation from the vertex to TOF
  static float ComputeTOFSignal(const float& trackTime, const float& expTime) { return trackTime * 1e+3 + expTime; }

  /// Returns the expected time of a track, given its TOF response to use to compute the expected time
  /// \param trk input track
  template <o2::track::PID::ID pid>
  static float GetTOFSignal(const TrackType& trk, const ResponseImplementation<pid>& response)
  {
    return ComputeTOFSignal(trk.trackTime(), response.GetExpectedSignal(trk));
  }

  /// Returns the expected time of a track considering the PID hypothesis used in tracking
  /// \param trk input track
  static float GetTOFSignal(const TrackType& trk)
  {
    switch (trk.pidForTracking()) {
      case 0:
        return GetTOFSignal(trk, responseEl);
      case 1:
        return GetTOFSignal(trk, responseMu);
      case 2:
        return GetTOFSignal(trk, responsePi);
      case 3:
        return GetTOFSignal(trk, responseKa);
      case 4:
        return GetTOFSignal(trk, responsePr);
      case 5:
        return GetTOFSignal(trk, responseDe);
      case 6:
        return GetTOFSignal(trk, responseTr);
      case 7:
        return GetTOFSignal(trk, responseHe);
      case 8:
        return GetTOFSignal(trk, responseAl);
      default:
        return 0.f;
        break;
    }
  }
};

//_________________________________________________________________________
template <typename TrackType, o2::track::PID::ID id>
float ExpTimes<TrackType, id>::GetExpectedSigmaFromTrackTime(const DetectorResponse& response, const TrackType& trk) const
{
  if (!trk.hasTOF()) {
    return -999.f;
  }
  if (!trk.has_collision()) {
    return -999.f;
  }
  const float x[7] = {trk.p(), o2::pid::tof::TOFSignal<TrackType>::GetTOFSignal(trk), trk.collision().collisionTimeRes() * 1000.f, o2::track::PID::getMass2Z(id), trk.length(), trk.sigma1Pt(), trk.pt()};
  const float reso = response(response.kSigma, x);
  return reso >= 0.f ? reso : 0.f;
}

//_________________________________________________________________________
template <typename TrackType, o2::track::PID::ID id>
float ExpTimes<TrackType, id>::GetSeparationFromTrackTime(const DetectorResponse& response, const TrackType& trk) const
{
  return trk.has_collision() ? (trk.hasTOF() ? (o2::pid::tof::TOFSignal<TrackType>::GetTOFSignal(trk) - trk.collision().collisionTime() * 1000.f - GetExpectedSignal(trk)) / GetExpectedSigmaFromTrackTime(response, trk) : -999.f) : -999.f;
}

} // namespace o2::pid::tof

#endif // O2_FRAMEWORK_PIDTOF_H_
