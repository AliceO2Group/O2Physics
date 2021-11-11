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

#ifndef O2_PID_TPC_RESPONSE_H_
#define O2_PID_TPC_RESPONSE_H_

#include <array>
#include <vector>
#include "Framework/Logger.h"
#include "ReconstructionDataFormats/PID.h"
#include "TMath.h"
#include "TPCSimulation/Detector.h"

namespace o2::pid::tpc
{

/// \brief Class to handle the TPC PID response

class Response
{

 public:
  Response() = default;
  ~Response() = default;

  /// Setter and Getter for the private parameters
  void SetBetheBlochParams(const std::array<float, 5>& betheBlochParams) { mBetheBlochParams = betheBlochParams; };
  void SetResolutionParams(const std::array<float, 2>& resolutionParams) { mResolutionParams = resolutionParams; };
  void SetMIP(const float mip) { mMIP = mip; };
  void SetChargeFactor(const float chargeFactor) { mChargeFactor = chargeFactor; };

  const std::array<float, 5> GetBetheBlochParams() const { return mBetheBlochParams; };
  const std::array<float, 2> GetResolutionParams() const { return mResolutionParams; };
  const float GetMIP() const { return mMIP; };
  const float GetChargeFactor() const { return mChargeFactor; };

  /// Gets the expected signal of the track
  template <typename TrackType>
  float GetExpectedSignal(const TrackType& track, const o2::track::PID::ID id) const;
  /// Gets the expected resolution of the track
  template <typename TrackType>
  float GetExpectedSigma(const TrackType& trk, const o2::track::PID::ID id) const;
  /// Gets the number of sigmas with respect the expected value
  template <typename TrackType>
  float GetNumberOfSigma(const TrackType& trk, const o2::track::PID::ID id) const;
  /// Gets the deviation to the expected signal
  template <typename TrackType>
  float GetSignalDelta(const TrackType& trk, const o2::track::PID::ID id) const;
  void PrintAll() const;

 private:
  std::array<float, 5> mBetheBlochParams = {0.0320981, 19.9768, 2.52666e-16, 2.72123, 6.08092};
  std::array<float, 2> mResolutionParams = {0.07, 0.0};
  float mMIP = 50.f;
  float mChargeFactor = 2.3f;

  ClassDefNV(Response, 1);

}; // class Response

/// Get expected Signal of the measurement
template <typename TrackType>
inline float Response::GetExpectedSignal(const TrackType& track, const o2::track::PID::ID id) const
{
  const float bethe = mMIP * o2::tpc::Detector::BetheBlochAleph(track.tpcInnerParam() / o2::track::pid_constants::sMasses[id], mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow((float)o2::track::pid_constants::sCharges[id], mChargeFactor);
  return bethe >= 0.f ? bethe : 0.f;
}

/// Gets the expected resolution of the measurement
template <typename TrackType>
inline float Response::GetExpectedSigma(const TrackType& track, const o2::track::PID::ID id) const
{
  const float reso = track.tpcSignal() * mResolutionParams[0] * ((float)track.tpcNClsFound() > 0 ? std::sqrt(1. + mResolutionParams[1] / (float)track.tpcNClsFound()) : 1.f);
  return reso >= 0.f ? reso : 0.f;
}

/// Gets the number of sigma between the actual signal and the expected signal
template <typename TrackType>
inline float Response::GetNumberOfSigma(const TrackType& trk, const o2::track::PID::ID id) const
{
  return ((trk.tpcSignal() - GetExpectedSignal(trk, id)) / GetExpectedSigma(trk, id));
}

/// Gets the deviation between the actual signal and the expected signal
template <typename TrackType>
inline float Response::GetSignalDelta(const TrackType& trk, const o2::track::PID::ID id) const
{
  return (trk.tpcSignal() - GetExpectedSignal(trk, id));
}

inline void Response::PrintAll() const
{
  LOGP(info, "==== TPC PID response parameters: ====");
  for (int i = 0; i < int(mBetheBlochParams.size()); i++)
    LOGP(info, "BB param [{}] = {}", i, mBetheBlochParams[i]);
  for (int i = 0; i < int(mResolutionParams.size()); i++)
    LOGP(info, "Resolution param [{}] = {}", i, mResolutionParams[i]);
  LOGP(info, "mMIP = {}", mMIP);
  LOGP(info, "mChargeFactor = {}", mChargeFactor);
}

} // namespace o2::pid::tpc

#endif // O2_FRAMEWORK_PIDRESPONSE_H_