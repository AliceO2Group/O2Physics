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
#include "Framework/DataTypes.h"
#include "Common/Core/PID/DetectorResponse.h"

namespace o2::pid::tof
{

// Utility values
static constexpr float kCSPEED = TMath::C() * 1.0e2f * 1.0e-12f; /// Speed of light in TOF units (cm/ps)
static constexpr float defaultReturnValue = -999.f;              /// Default return value in case TOF measurement is not available

/// \brief Class to handle the the TOF detector response for the TOF beta measurement
template <typename TrackType>
class Beta
{
 public:
  Beta() = default;
  ~Beta() = default;

  /// Computes the beta of a track given a length, a time measurement and an event time (in ps)
  /// \param length Length in cm of the track
  /// \param tofSignal TOF signal in ps for the track
  /// \param collisionTime collision time in ps for the event of the track
  static float GetBeta(const float& length, const float& tofSignal, const float& collisionTime) { return length / (tofSignal - collisionTime) / kCSPEED; }

  /// Gets the beta for the track of interest
  /// \param track Track of interest
  /// \param collisionTime Collision time
  static float GetBeta(const TrackType& track, const float& collisionTime) { return track.hasTOF() ? GetBeta(track.length(), track.tofSignal(), collisionTime) : defaultReturnValue; }

  /// Gets the beta for the track of interest
  /// \param track Track of interest
  static float GetBeta(const TrackType& track) { return track.isEvTimeDefined() ? GetBeta(track, track.tofEvTime()) : defaultReturnValue; }

  /// Computes the expected uncertainty on the beta measurement
  /// \param length Length in cm of the track
  /// \param tofSignal TOF signal in ps for the track
  /// \param collisionTime collision time in ps for the event of the track
  /// \param time_reso expected time resolution
  static float GetExpectedSigma(const float& length, const float& tofSignal, const float& collisionTime, const float& expectedResolution) { return GetBeta(length, tofSignal, collisionTime) / (tofSignal - collisionTime) * expectedResolution; }

  /// Gets the expected uncertainty on the beta measurement of the track of interest
  /// \param track Track of interest
  float GetExpectedSigma(const TrackType& track) const { return track.isEvTimeDefined() ? GetExpectedSigma(track.length(), track.tofSignal(), track.tofEvTime(), mExpectedResolution) : defaultReturnValue; }

  /// Gets the expected beta for a given mass hypothesis (no energy loss taken into account)
  /// \param momentum momentum in GeV/c of the track
  /// \param mass mass in GeV/c2 of the particle of interest
  static float GetExpectedSignal(const float& momentum, const float& mass) { return momentum > 0 ? momentum / std::sqrt(momentum * momentum + mass * mass) : 0.f; }

  /// Gets the expected beta given the particle index (no energy loss taken into account) of the track of interest
  /// \param track Track of interest
  template <o2::track::PID::ID id>
  float GetExpectedSignal(const TrackType& track) const
  {
    return GetExpectedSignal(track.p(), o2::track::PID::getMass2Z(id));
  }

  /// Gets the number of sigmas with respect the approximate beta (no energy loss taken into account) of the track of interest
  /// \param track Track of interest
  template <o2::track::PID::ID id>
  float GetSeparation(const TrackType& track) const
  {
    return (GetBeta(track) - GetExpectedSignal<id>(track)) / GetExpectedSigma(track);
  }

  float mExpectedResolution = 80; /// Expected time resolution
};

/// \brief Class to handle the the TOF detector response for the TOF mass measurement
template <typename TrackType>
class TOFMass
{
 public:
  TOFMass() = default;
  ~TOFMass() = default;

  /// Computes the TOF mass of a track given a momentum, a beta measurement
  /// \param momentum momentum of the track
  /// \param beta TOF beta measurement
  static float GetTOFMass(const float& momentum, const float& beta) { return (momentum / beta) * std::sqrt(std::abs(1.f - beta * beta)); }

  /// Gets the TOF mass for the track of interest
  /// \param track Track of interest
  static float GetTOFMass(const TrackType& track, const float beta) { return track.hasTOF() ? GetTOFMass(track.p(), beta) : defaultReturnValue; }

  /// Gets the TOF mass for the track of interest
  /// \param track Track of interest
  static float GetTOFMass(const TrackType& track) { return track.hasTOF() ? GetTOFMass(track.p(), Beta<TrackType>::GetBeta(track)) : defaultReturnValue; }
};

/// \brief Class to handle the the TOF detector response for the expected time
template <typename TrackType, o2::track::PID::ID id>
class ExpTimes
{
 public:
  ExpTimes() = default;
  ~ExpTimes() = default;

  /// Computes the expected time of a track, given it TOF expected momentum
  static float ComputeExpectedTime(const float& tofExpMom, const float& length, const float& massZ) { return length * sqrt((massZ * massZ) + (tofExpMom * tofExpMom)) / (kCSPEED * tofExpMom); }

  /// Gets the expected signal of the track of interest under the PID assumption
  /// \param track Track of interest
  static float GetExpectedSignal(const TrackType& track)
  {
    if (!track.hasTOF()) {
      return defaultReturnValue;
    }
    if (track.trackType() == o2::aod::track::Run2Track) {
      return ComputeExpectedTime(track.tofExpMom() / kCSPEED, track.length(), o2::track::PID::getMass2Z(id));
    }
    return ComputeExpectedTime(track.tofExpMom(), track.length(), o2::track::PID::getMass2Z(id));
  }

  /// Gets the expected resolution of the t-texp-t0
  /// Given a TOF signal and collision time resolutions
  /// \param response Detector response with parameters
  /// \param track Track of interest
  /// \param tofSignal TOF signal of the track of interest
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetExpectedSigma(const DetectorResponse& response, const TrackType& track, const float& tofSignal, const float& collisionTimeRes)
  {
    if (!track.hasTOF()) {
      return defaultReturnValue;
    }
    // const float x[7] = {track.p(), tofSignal, collisionTimeRes, o2::track::PID::getMass2Z(id), track.length(), track.sigma1Pt(), track.pt()};
    const float x[4] = {track.p(), tofSignal, collisionTimeRes, o2::track::PID::getMass2Z(id)};
    const float reso = response(response.kSigma, x);
    return reso >= 0.f ? reso : 0.f;
  }

  /// Gets the expected resolution of the t-texp-t0
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetExpectedSigma(const DetectorResponse& response, const TrackType& track) { return track.isEvTimeDefined() ? GetExpectedSigma(response, track, track.tofSignal(), track.tofEvTime()) : defaultReturnValue; }

  /// Gets the expected resolution of the time measurement, uses the expected time and no event time resolution
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetExpectedSigmaTracking(const DetectorResponse& response, const TrackType& track) { return GetExpectedSigma(response, track, GetExpectedSignal(track), 0.f); }

  /// Gets the separation between the measured signal and the expected one
  /// \param track Track of interest
  /// \param collisionTime Collision time
  static float GetDelta(const TrackType& track, const float& collisionTime) { return track.hasTOF() ? (track.tofSignal() - collisionTime - GetExpectedSignal(track)) : defaultReturnValue; }

  /// Gets the number of sigmas with respect the expected time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  /// \param collisionTime Collision time
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetSeparation(const DetectorResponse& response, const TrackType& track, const float& collisionTime, const float& collisionTimeRes) { return track.hasTOF() ? GetDelta(track, collisionTime) / GetExpectedSigma(response, track, track.tofSignal(), collisionTimeRes) : defaultReturnValue; }

  /// Gets the number of sigmas with respect the expected time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetSeparation(const DetectorResponse& response, const TrackType& track) { return track.isEvTimeDefined() ? GetSeparation(response, track, track.tofEvTime(), track.tofEvTimeErr()) : defaultReturnValue; }

  /// Gets the expected resolution of the measurement from the track time and from the collision time, explicitly passed as argument
  /// \param response Detector response with parameters
  /// \param track Track of interest
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetExpectedSigmaFromTrackTime(const DetectorResponse& response, const TrackType& track, const float& collisionTimeRes);

  /// Gets the expected resolution of the measurement from the track time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetExpectedSigmaFromTrackTime(const DetectorResponse& response, const TrackType& track) { return track.isEvTimeDefined() ? GetExpectedSigmaFromTrackTime(response, track, track.tofEvTimeErr()) : defaultReturnValue; }

  /// Gets the number of sigmas with respect the expected time from the track time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  /// \param collisionTime Collision time
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetSeparationFromTrackTime(const DetectorResponse& response, const TrackType& track, const float& collisionTime, const float& collisionTimeRes);

  /// Gets the number of sigmas with respect the expected time from the track time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetSeparationFromTrackTime(const DetectorResponse& response, const TrackType& track) { return track.isEvTimeDefined() ? GetSeparationFromTrackTime(response, track, track.tofEvTime(), track.tofEvTimeErr()) : defaultReturnValue; }
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
  static float ComputeTOFSignal(const float& trackTime, const float& expTime) { return trackTime * 1000.f + expTime; }

  /// Returns the expected time of a track, given its TOF response to use to compute the expected time
  /// \param track input track
  template <o2::track::PID::ID pid>
  static float GetTOFSignal(const TrackType& track, const ResponseImplementation<pid>& response)
  {
    return ComputeTOFSignal(track.trackTime(), response.GetExpectedSignal(track));
  }

  /// Returns the expected time of a track considering the PID hypothesis used in tracking
  /// \param track input track
  static float GetTOFSignal(const TrackType& track)
  {
    switch (track.pidForTracking()) {
      case 0:
        return GetTOFSignal(track, responseEl);
      case 1:
        return GetTOFSignal(track, responseMu);
      case 2:
        return GetTOFSignal(track, responsePi);
      case 3:
        return GetTOFSignal(track, responseKa);
      case 4:
        return GetTOFSignal(track, responsePr);
      case 5:
        return GetTOFSignal(track, responseDe);
      case 6:
        return GetTOFSignal(track, responseTr);
      case 7:
        return GetTOFSignal(track, responseHe);
      case 8:
        return GetTOFSignal(track, responseAl);
      default:
        return 0.f;
        break;
    }
  }
};

//_________________________________________________________________________
template <typename TrackType, o2::track::PID::ID id>
float ExpTimes<TrackType, id>::GetExpectedSigmaFromTrackTime(const DetectorResponse& response, const TrackType& track, const float& collisionTimeRes)
{
  return GetExpectedSigma(response, track, o2::pid::tof::TOFSignal<TrackType>::GetTOFSignal(track), collisionTimeRes);
}

//_________________________________________________________________________
template <typename TrackType, o2::track::PID::ID id>
float ExpTimes<TrackType, id>::GetSeparationFromTrackTime(const DetectorResponse& response, const TrackType& track, const float& collisionTime, const float& collisionTimeRes)
{
  return track.hasTOF() ? (o2::pid::tof::TOFSignal<TrackType>::GetTOFSignal(track) - collisionTime - GetExpectedSignal(track)) / GetExpectedSigmaFromTrackTime(response, track, collisionTimeRes) : defaultReturnValue;
}

} // namespace o2::pid::tof

#endif // O2_FRAMEWORK_PIDTOF_H_
