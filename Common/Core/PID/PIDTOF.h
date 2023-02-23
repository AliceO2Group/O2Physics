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
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  02/07/2020
/// \brief  Implementation of the TOF detector response for PID
///

#ifndef COMMON_CORE_PID_PIDTOF_H_
#define COMMON_CORE_PID_PIDTOF_H_

#include <string>

// ROOT includes
#include "Rtypes.h"
#include "TMath.h"

// O2 includes
#include "DataFormatsTOF/ParameterContainers.h"
#include "Framework/Logger.h"
#include "ReconstructionDataFormats/PID.h"
#include "Framework/DataTypes.h"

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
  static float GetBeta(const float length, const float tofSignal, const float collisionTime) { return length / (tofSignal - collisionTime) / kCSPEED; }

  /// Gets the beta for the track of interest
  /// \param track Track of interest
  /// \param collisionTime Collision time
  static float GetBeta(const TrackType& track, const float collisionTime) { return track.hasTOF() ? GetBeta(track.length(), track.tofSignal(), collisionTime) : defaultReturnValue; }

  /// Gets the beta for the track of interest
  /// \param track Track of interest
  static float GetBeta(const TrackType& track) { return GetBeta(track, track.tofEvTime()); }

  /// Computes the expected uncertainty on the beta measurement
  /// \param length Length in cm of the track
  /// \param tofSignal TOF signal in ps for the track
  /// \param collisionTime collision time in ps for the event of the track
  /// \param time_reso expected time resolution
  static float GetExpectedSigma(const float length, const float tofSignal, const float collisionTime, const float expectedResolution) { return GetBeta(length, tofSignal, collisionTime) / (tofSignal - collisionTime) * expectedResolution; }

  /// Gets the expected uncertainty on the beta measurement of the track of interest
  /// \param track Track of interest
  float GetExpectedSigma(const TrackType& track) const { return GetExpectedSigma(track.length(), track.tofSignal(), track.tofEvTime(), mExpectedResolution); }

  /// Gets the expected beta for a given mass hypothesis (no energy loss taken into account)
  /// \param momentum momentum in GeV/c of the track
  /// \param mass mass in GeV/c2 of the particle of interest
  static float GetExpectedSignal(const float momentum, const float mass) { return momentum > 0 ? momentum / std::sqrt(momentum * momentum + mass * mass) : 0.f; }

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
  static float GetTOFMass(const float momentum, const float beta) { return (momentum / beta) * std::sqrt(std::abs(1.f - beta * beta)); }

  /// Gets the TOF mass for the track of interest
  /// \param track Track of interest
  static float GetTOFMass(const TrackType& track, const float beta) { return track.hasTOF() ? GetTOFMass(track.p(), beta) : defaultReturnValue; }

  /// Gets the TOF mass for the track of interest
  /// \param track Track of interest
  static float GetTOFMass(const TrackType& track) { return track.hasTOF() ? GetTOFMass(track.p(), Beta<TrackType>::GetBeta(track)) : defaultReturnValue; }
};

/// \brief Implementation class to store TOF response parameters for exp. times
class TOFResoParams : public o2::pid::PidParameters<5>
{
 public:
  TOFResoParams() : PidParameters("TOFResoParams") { SetParameters(std::array<float, 5>{0.008, 0.008, 0.002, 40.0, 60.0}); } // Default constructor with default parameters
  ~TOFResoParams() = default;
  ClassDef(TOFResoParams, 1);
};

/// \brief Next implementation class to store TOF response parameters for exp. times
class TOFResoParamsV2 : public o2::tof::Parameters<5>
{
 public:
  TOFResoParamsV2() : Parameters(std::array<std::string, 5>{"p0", "p1", "p2", "p3", "time_resolution"},
                                 "ParamExample") { setParameters(std::array<float, 5>{0.008, 0.008, 0.002, 40.0, 60.0}); } // Default constructor with default parameters

  ~TOFResoParamsV2() = default;
};

/// \brief Class to handle the the TOF detector response for the expected time
template <typename TrackType, o2::track::PID::ID id>
class ExpTimes
{
 public:
  ExpTimes() = default;
  ~ExpTimes() = default;
  static constexpr float mMassZ = o2::track::pid_constants::sMasses2Z[id]; /// Mass hypothesis divided by the charge (in units of e): M/z. Equivalent to o2::track::PID::getMass2Z(id)
  static constexpr float mMassZSqared = mMassZ * mMassZ;                   /// (M/z)^2

  /// Computes the expected time of a track, given it TOF expected momentum
  static float ComputeExpectedTime(const float tofExpMom, const float length) { return length * sqrt((mMassZSqared) + (tofExpMom * tofExpMom)) / (kCSPEED * tofExpMom); }

  /// Gets the expected signal of the track of interest under the PID assumption
  /// \param track Track of interest
  static float GetExpectedSignal(const TrackType& track)
  {
    if (!track.hasTOF()) {
      return defaultReturnValue;
    }
    if (track.trackType() == o2::aod::track::Run2Track) {
      return ComputeExpectedTime(track.tofExpMom() / kCSPEED, track.length());
    }
    return ComputeExpectedTime(track.tofExpMom(), track.length());
  }

  /// Gets the expected resolution of the t-texp-t0
  /// Given a TOF signal and collision time resolutions
  /// \param response Detector response with parameters
  /// \param track Track of interest
  /// \param tofSignal TOF signal of the track of interest
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetExpectedSigma(const DetectorResponse& response, const TrackType& track, const float tofSignal, const float collisionTimeRes)
  {
    if (!track.hasTOF()) {
      return defaultReturnValue;
    }
    // const float x[7] = {track.p(), tofSignal, collisionTimeRes, o2::track::PID::getMass2Z(id), track.length(), track.sigma1Pt(), track.pt()};
    const float x[4] = {track.p(), tofSignal, collisionTimeRes, mMassZ};
    const float reso = response(response.kSigma, x);
    return reso >= 0.f ? reso : 0.f;
  }

  /// Gets the expected resolution of the t-texp-t0
  /// Given a TOF signal and collision time resolutions
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  /// \param tofSignal TOF signal of the track of interest
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetExpectedSigma(const TOFResoParams& parameters, const TrackType& track, const float tofSignal, const float collisionTimeRes)
  {
    const float mom = track.p();
    if (mom <= 0) {
      return -999.f;
    }
    const float dpp = parameters[0] + parameters[1] * mom + parameters[2] * mMassZ / mom; // mean relative pt resolution;
    const float sigma = dpp * tofSignal / (1. + mom * mom / (mMassZSqared));
    return std::sqrt(sigma * sigma + parameters[3] * parameters[3] / mom / mom + parameters[4] * parameters[4] + collisionTimeRes * collisionTimeRes);
  }

  /// Gets the expected resolution of the t-texp-t0
  /// Given a TOF signal and collision time resolutions
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  /// \param tofSignal TOF signal of the track of interest
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetExpectedSigma(const TOFResoParamsV2& parameters, const TrackType& track, const float tofSignal, const float collisionTimeRes)
  {
    const float mom = track.p();
    if (mom <= 0) {
      return -999.f;
    }
    const float dpp = parameters[0] + parameters[1] * mom + parameters[2] * mMassZ / mom; // mean relative pt resolution;
    const float sigma = dpp * tofSignal / (1. + mom * mom / (mMassZSqared));
    return std::sqrt(sigma * sigma + parameters[3] * parameters[3] / mom / mom + parameters[4] * parameters[4] + collisionTimeRes * collisionTimeRes);
  }

  /// Gets the expected resolution of the t-texp-t0
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetExpectedSigma(const DetectorResponse& response, const TrackType& track) { return GetExpectedSigma(response, track, track.tofSignal(), track.tofEvTime()); }

  /// Gets the expected resolution of the t-texp-t0
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  static float GetExpectedSigma(const TOFResoParams& parameters, const TrackType& track) { return GetExpectedSigma(parameters, track, track.tofSignal(), track.tofEvTime()); }

  /// Gets the expected resolution of the t-texp-t0
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  static float GetExpectedSigma(const TOFResoParamsV2& parameters, const TrackType& track) { return GetExpectedSigma(parameters, track, track.tofSignal(), track.tofEvTime()); }

  /// Gets the expected resolution of the time measurement, uses the expected time and no event time resolution
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetExpectedSigmaTracking(const DetectorResponse& response, const TrackType& track) { return GetExpectedSigma(response, track, GetExpectedSignal(track), 0.f); }

  /// Gets the expected resolution of the time measurement, uses the expected time and no event time resolution
  /// \param parameters Parameters to use to compute the expected resolution
  /// \param track Track of interest
  static float GetExpectedSigmaTracking(const TOFResoParamsV2& parameters, const TrackType& track) { return GetExpectedSigma(parameters, track, GetExpectedSignal(track), 0.f); }

  /// Gets the separation between the measured signal and the expected one
  /// \param track Track of interest
  /// \param collisionTime Collision time
  static float GetDelta(const TrackType& track, const float collisionTime) { return track.hasTOF() ? (track.tofSignal() - collisionTime - GetExpectedSignal(track)) : defaultReturnValue; }

  /// Gets the number of sigmas with respect the expected time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  /// \param collisionTime Collision time
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetSeparation(const DetectorResponse& response, const TrackType& track, const float collisionTime, const float collisionTimeRes) { return track.hasTOF() ? GetDelta(track, collisionTime) / GetExpectedSigma(response, track, track.tofSignal(), collisionTimeRes) : defaultReturnValue; }

  /// Gets the number of sigmas with respect the expected time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetSeparation(const DetectorResponse& response, const TrackType& track) { return GetSeparation(response, track, track.tofEvTime(), track.tofEvTimeErr()); }

  /// Gets the number of sigmas with respect the expected time
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  /// \param collisionTime Collision time
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetSeparation(const TOFResoParams& parameters, const TrackType& track, const float collisionTime, const float collisionTimeRes) { return track.hasTOF() ? GetDelta(track, collisionTime) / GetExpectedSigma(parameters, track, track.tofSignal(), collisionTimeRes) : defaultReturnValue; }

  /// Gets the number of sigmas with respect the expected time
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  static float GetSeparation(const TOFResoParams& parameters, const TrackType& track) { return GetSeparation(parameters, track, track.tofEvTime(), track.tofEvTimeErr()); }

  /// Gets the number of sigmas with respect the expected time
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  /// \param collisionTime Collision time
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetSeparation(const TOFResoParamsV2& parameters, const TrackType& track, const float collisionTime, const float collisionTimeRes) { return track.hasTOF() ? GetDelta(track, collisionTime) / GetExpectedSigma(parameters, track, track.tofSignal(), collisionTimeRes) : defaultReturnValue; }

  /// Gets the number of sigmas with respect the expected time
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  static float GetSeparation(const TOFResoParamsV2& parameters, const TrackType& track) { return GetSeparation(parameters, track, track.tofEvTime(), track.tofEvTimeErr()); }

  /// Gets the expected resolution of the measurement from the track time and from the collision time, explicitly passed as argument
  /// \param response Detector response with parameters
  /// \param track Track of interest
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetExpectedSigmaFromTrackTime(const DetectorResponse& response, const TrackType& track, const float collisionTimeRes);

  /// Gets the expected resolution of the measurement from the track time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetExpectedSigmaFromTrackTime(const DetectorResponse& response, const TrackType& track) { return GetExpectedSigmaFromTrackTime(response, track, track.tofEvTimeErr()); }

  /// Gets the number of sigmas with respect the expected time from the track time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  /// \param collisionTime Collision time
  /// \param collisionTimeRes Collision time resolution of the track of interest
  static float GetSeparationFromTrackTime(const DetectorResponse& response, const TrackType& track, const float collisionTime, const float collisionTimeRes);

  /// Gets the number of sigmas with respect the expected time from the track time
  /// \param response Detector response with parameters
  /// \param track Track of interest
  static float GetSeparationFromTrackTime(const DetectorResponse& response, const TrackType& track) { return GetSeparationFromTrackTime(response, track, track.tofEvTime(), track.tofEvTimeErr()); }
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
  static float ComputeTOFSignal(const float trackTime, const float expTime) { return trackTime * 1000.f + expTime; }

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
float ExpTimes<TrackType, id>::GetExpectedSigmaFromTrackTime(const DetectorResponse& response, const TrackType& track, const float collisionTimeRes)
{
  return GetExpectedSigma(response, track, o2::pid::tof::TOFSignal<TrackType>::GetTOFSignal(track), collisionTimeRes);
}

//_________________________________________________________________________
template <typename TrackType, o2::track::PID::ID id>
float ExpTimes<TrackType, id>::GetSeparationFromTrackTime(const DetectorResponse& response, const TrackType& track, const float collisionTime, const float collisionTimeRes)
{
  return track.hasTOF() ? (o2::pid::tof::TOFSignal<TrackType>::GetTOFSignal(track) - collisionTime - GetExpectedSignal(track)) / GetExpectedSigmaFromTrackTime(response, track, collisionTimeRes) : defaultReturnValue;
}

} // namespace o2::pid::tof

#endif // COMMON_CORE_PID_PIDTOF_H_
