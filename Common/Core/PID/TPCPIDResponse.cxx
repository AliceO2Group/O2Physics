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
/// \file   TPCPIDResponse.h
/// \author Annalena Kalteyer, Christian Sonnabend
/// \brief

#include "Common/Core/PID/TPCPIDResponse.h"

using namespace o2::pid::tpc;

/// @param bg Beta*Gamma of the incident particle
/// @param kp* Parameters for the ALICE TPC
/// @return Bethe-Bloch value in MIP units

template <typename TrackType, o2::track::PID::ID id>
template <typename T>
inline T TPCPIDResponse<TrackType, id>::BetheBlochAleph(T bg, T kp1, T kp2, T kp3, T kp4, T kp5)
{
  T beta = bg / std::sqrt(static_cast<T>(1.) + bg * bg);
  T aa = std::pow(beta, kp4);
  T bb = std::pow(static_cast<T>(1.) / bg, kp5);
  bb = std::log(kp3 + bb);
  return (kp2 - aa - bb) * kp1 / aa;
}

/// Get expected Signal of the measurement
template <typename TrackType, o2::track::PID::ID id>
inline float TPCPIDResponse<TrackType, id>::GetExpectedSignal(const TrackType& track) const
{
  const float bethe = mChargeFactor * BetheBlochAleph(track.tpcInnerParam() / o2::track::pid_constants::sMasses[id], mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * pow((float)o2::track::pid_constants::sCharges[id], mMIP);
  return bethe >= 0.f ? bethe : 0.f;
}

/// Gets the expected resolution of the measurement
template <typename TrackType, o2::track::PID::ID id>
inline float TPCPIDResponse<TrackType, id>::GetExpectedSigma(const TrackType& track) const
{
  const float reso = track.tpcSignal() * mResolutionParams[0] * ((float)track.tpcNClsFound() > 0 ? sqrt(1. + mResolutionParams[1] / (float)track.tpcNClsFound()) : 1.f);
  return reso >= 0.f ? reso : 0.f;
}

/// Gets the separation between the actual signal and the expected signal
template <typename TrackType, o2::track::PID::ID id>
inline float TPCPIDResponse<TrackType, id>::GetSeparation(const TrackType& trk) const
{
  return ((trk.tpcSignal() - GetExpectedSignal(trk)) / GetExpectedSigma(trk));
}; // namespace o2::pid::tpc