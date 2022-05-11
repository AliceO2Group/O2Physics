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
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend
/// \brief

#ifndef O2_PID_TPC_RESPONSE_H_
#define O2_PID_TPC_RESPONSE_H_

#include <array>
#include <vector>
#include <cmath>
#include "Framework/Logger.h"
// O2 includes
#include "ReconstructionDataFormats/PID.h"
#include "DataFormatsTPC/BetheBlochAleph.h"

namespace o2::pid::tpc
{

/// \brief Class to handle the TPC PID response

class Response
{

 public:
  Response() = default;
  ~Response() = default;

  /// Setter and Getter for the private parameters
  void SetBetheBlochParams(const std::array<float, 5>& betheBlochParams) { mBetheBlochParams = betheBlochParams; }
  void SetResolutionParamsDefault(const std::array<float, 2>& resolutionParamsDefault) { mResolutionParamsDefault = resolutionParamsDefault; }
  void SetResolutionParams(const std::vector<double>& resolutionParams) { mResolutionParams = resolutionParams; }
  void SetMIP(const float mip) { mMIP = mip; }
  void SetChargeFactor(const float chargeFactor) { mChargeFactor = chargeFactor; }
  void SetMultiplicityNormalization(const float multNormalization) { mMultNormalization = multNormalization; }
  void SetUseDefaultResolutionParam(const bool useDefault) { mUseDefaultResolutionParam = useDefault; }
  void SetParameters(const Response* response)
  {
    mBetheBlochParams = response->GetBetheBlochParams();
    mResolutionParamsDefault = response->GetResolutionParamsDefault();
    mResolutionParams = response->GetResolutionParams();
    mMIP = response->GetMIP();
    mChargeFactor = response->GetChargeFactor();
    mMultNormalization = response->GetMultiplicityNormalization();
    mUseDefaultResolutionParam = response->GetUseDefaultResolutionParam();
  }

  const std::array<float, 5> GetBetheBlochParams() const { return mBetheBlochParams; }
  const std::array<float, 2> GetResolutionParamsDefault() const { return mResolutionParamsDefault; }
  const std::vector<double> GetResolutionParams() const { return mResolutionParams; }
  const float GetMIP() const { return mMIP; }
  const float GetChargeFactor() const { return mChargeFactor; }
  const float GetMultiplicityNormalization() const { return mMultNormalization; }
  const bool GetUseDefaultResolutionParam() const { return mUseDefaultResolutionParam; }

  /// Gets the expected signal of the track
  template <typename TrackType>
  float GetExpectedSignal(const TrackType& track, const o2::track::PID::ID id) const;
  /// Gets the expected resolution of the track
  template <typename CollisionType, typename TrackType>
  float GetExpectedSigma(const CollisionType& collision, const TrackType& trk, const o2::track::PID::ID id) const;
  /// Gets the number of sigmas with respect the expected value
  template <typename CollisionType, typename TrackType>
  float GetNumberOfSigma(const CollisionType& collision, const TrackType& trk, const o2::track::PID::ID id) const;
  /// Gets the deviation to the expected signal
  template <typename TrackType>
  float GetSignalDelta(const TrackType& trk, const o2::track::PID::ID id) const;
  /// Gets relative dEdx resolution contribution due to relative pt resolution
  float GetRelativeResolutiondEdx(const float p, const float mass, const float charge, const float resol) const;

  void PrintAll() const;

 private:
  std::array<float, 5> mBetheBlochParams = {0.03209809958934784, 19.9768009185791, 2.5266601063857674e-16, 2.7212300300598145, 6.080920219421387};
  std::array<float, 2> mResolutionParamsDefault = {0.07, 0.0};
  std::vector<double> mResolutionParams = {5.43799e-7, 0.053044, 0.667584, 0.0142667, 0.00235175, 1.22482, 2.3501e-7, 0.031585};
  float mMIP = 50.f;
  float mChargeFactor = 2.299999952316284f;
  float mMultNormalization = 11000.;
  bool mUseDefaultResolutionParam = true;

  ClassDefNV(Response, 2);

}; // class Response

/// Get expected Signal of the measurement
template <typename TrackType>
inline float Response::GetExpectedSignal(const TrackType& track, const o2::track::PID::ID id) const
{
  const float bethe = mMIP * o2::tpc::BetheBlochAleph(track.tpcInnerParam() / o2::track::pid_constants::sMasses[id], mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow((float)o2::track::pid_constants::sCharges[id], mChargeFactor);
  return bethe >= 0.f ? bethe : 0.f;
}

/// Gets the expected resolution of the measurement
template <typename CollisionType, typename TrackType>
inline float Response::GetExpectedSigma(const CollisionType& collision, const TrackType& track, const o2::track::PID::ID id) const
{
  float resolution = 0.;
  if (mUseDefaultResolutionParam) {
    const float reso = track.tpcSignal() * mResolutionParamsDefault[0] * ((float)track.tpcNClsFound() > 0 ? std::sqrt(1. + mResolutionParamsDefault[1] / (float)track.tpcNClsFound()) : 1.f);
    reso >= 0.f ? resolution = reso : resolution = 0.f;
  } else {

    const double ncl = 159. / track.tpcNClsFound(); //
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const double dEdx = o2::tpc::BetheBlochAleph((float)bg, mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow((float)o2::track::pid_constants::sCharges[id], mChargeFactor);
    const double relReso = GetRelativeResolutiondEdx(p, mass, o2::track::pid_constants::sCharges[id], mResolutionParams[3]);

    const std::vector<double> values{1.f / dEdx, track.tgl(), std::sqrt(ncl), relReso, track.signed1Pt(), collision.multTPC() / mMultNormalization};

    const float reso = sqrt(pow(mResolutionParams[0], 2) * values[0] + pow(mResolutionParams[1], 2) * (values[2] * mResolutionParams[5]) * pow(values[0] / sqrt(1 + pow(values[1], 2)), mResolutionParams[2]) + values[2] * pow(values[3], 2) + pow(mResolutionParams[4] * values[4], 2) + pow(values[5] * mResolutionParams[6], 2) + pow(values[5] * (values[0] / sqrt(1 + pow(values[1], 2))) * mResolutionParams[7], 2)) * dEdx * mMIP;
    reso >= 0.f ? resolution = reso : resolution = 0.f;
  }
  return resolution;
}

/// Gets the number of sigma between the actual signal and the expected signal
template <typename CollisionType, typename TrackType>
inline float Response::GetNumberOfSigma(const CollisionType& collision, const TrackType& trk, const o2::track::PID::ID id) const
{
  return ((trk.tpcSignal() - GetExpectedSignal(trk, id)) / GetExpectedSigma(collision, trk, id));
}

/// Gets the deviation between the actual signal and the expected signal
template <typename TrackType>
inline float Response::GetSignalDelta(const TrackType& trk, const o2::track::PID::ID id) const
{
  return (trk.tpcSignal() - GetExpectedSignal(trk, id));
}

//// Gets relative dEdx resolution contribution due relative pt resolution
inline float Response::GetRelativeResolutiondEdx(const float p, const float mass, const float charge, const float resol) const
{
  const float bg = p / mass;
  const float dEdx = o2::tpc::BetheBlochAleph(bg, mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow(charge, mChargeFactor);
  const float deltaP = resol * std::sqrt(dEdx);
  const float bgDelta = p * (1 + deltaP) / mass;
  const float dEdx2 = o2::tpc::BetheBlochAleph(bgDelta, mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow(charge, mChargeFactor);
  const float deltaRel = std::abs(dEdx2 - dEdx) / dEdx;
  return deltaRel;
}

inline void Response::PrintAll() const
{
  LOGP(info, "==== TPC PID response parameters: ====");
  for (int i = 0; i < int(mBetheBlochParams.size()); i++) {
    LOGP(info, "BB param [{}] = {}", i, mBetheBlochParams[i]);
  }
  LOGP(info, "use default resolution parametrization = {}", mUseDefaultResolutionParam);
  if (mUseDefaultResolutionParam) {
    LOGP(info, "Default Resolution parametrization: ");
    for (int i = 0; i < int(mResolutionParamsDefault.size()); i++) {
      LOGP(info, "Resolution param [{}] = {}", i, mResolutionParamsDefault[i]);
    }
  } else {
    for (int i = 0; i < int(mResolutionParams.size()); i++) {
      LOGP(info, "Resolution param [{}] = {}", i, mResolutionParams[i]);
    }
  }
  LOGP(info, "mMIP = {}", mMIP);
  LOGP(info, "mChargeFactor = {}", mChargeFactor);
  LOGP(info, "mMultNormalization = {}", mMultNormalization);
}

} // namespace o2::pid::tpc

#endif // O2_FRAMEWORK_PIDRESPONSE_H_