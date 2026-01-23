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

#ifndef COMMON_CORE_PID_TPCPIDRESPONSE_H_
#define COMMON_CORE_PID_TPCPIDRESPONSE_H_

#include <MathUtils/BetheBlochAleph.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/PID.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <vector>

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
  void SetNClNormalization(const float nclnorm) { nClNorm = nclnorm; }
  void SetUseDefaultResolutionParam(const bool useDefault) { mUseDefaultResolutionParam = useDefault; }
  void SetParameters(const Response* response)
  {
    mBetheBlochParams = response->GetBetheBlochParams();
    mResolutionParamsDefault = response->GetResolutionParamsDefault();
    mResolutionParams = response->GetResolutionParams();
    mMIP = response->GetMIP();
    nClNorm = response->GetNClNormalization();
    mChargeFactor = response->GetChargeFactor();
    mMultNormalization = response->GetMultiplicityNormalization();
    mUseDefaultResolutionParam = response->GetUseDefaultResolutionParam();
  }

  const std::array<float, 5> GetBetheBlochParams() const { return mBetheBlochParams; }
  const std::array<float, 2> GetResolutionParamsDefault() const { return mResolutionParamsDefault; }
  const std::vector<double> GetResolutionParams() const { return mResolutionParams; }
  float GetMIP() const { return mMIP; }
  float GetNClNormalization() const { return nClNorm; }
  float GetChargeFactor() const { return mChargeFactor; }
  float GetMultiplicityNormalization() const { return mMultNormalization; }
  bool GetUseDefaultResolutionParam() const { return mUseDefaultResolutionParam; }

  /// Gets the expected signal of the track
  template <typename TrackType>
  float GetExpectedSignal(const TrackType& track, const o2::track::PID::ID id) const;
  /// Gets the expected resolution of the track
  template <typename CollisionType, typename TrackType>
  float GetExpectedSigma(const CollisionType& collision, const TrackType& trk, const o2::track::PID::ID id) const;
  /// Gets the expected resolution of the track with multTPC explicitly provided
  template <typename TrackType>
  float GetExpectedSigmaAtMultiplicity(const long multTPC, const TrackType& trk, const o2::track::PID::ID id) const;
  /// Gets the number of sigmas with respect the expected value
  template <typename CollisionType, typename TrackType>
  float GetNumberOfSigma(const CollisionType& collision, const TrackType& trk, const o2::track::PID::ID id) const;
  // Number of sigmas with respect to expected for MC, defining a tune-on-data signal value
  template <typename CollisionType, typename TrackType>
  float GetNumberOfSigmaMCTuned(const CollisionType& collision, const TrackType& trk, const o2::track::PID::ID id, float mcTunedTPCSignal) const;
  // Number of sigmas with respect to expected for MC, defining a tune-on-data signal value, explicit multTPC
  template <typename TrackType>
  float GetNumberOfSigmaMCTunedAtMultiplicity(const long multTPC, const TrackType& trk, const o2::track::PID::ID id, float mcTunedTPCSignal) const;
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
  float nClNorm = 152.f;

  ClassDefNV(Response, 3);

}; // class Response

/// Get expected Signal of the measurement
template <typename TrackType>
inline float Response::GetExpectedSignal(const TrackType& track, const o2::track::PID::ID id) const
{
  if (!track.hasTPC()) {
    return -999.f;
  }
  const float bethe = mMIP * o2::common::BetheBlochAleph(track.tpcInnerParam() / o2::track::pid_constants::sMasses[id], mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow(static_cast<float>(o2::track::pid_constants::sCharges[id]), mChargeFactor);
  return bethe >= 0.f ? bethe : -999.f;
}

/// Gets the expected resolution of the measurement
template <typename CollisionType, typename TrackType>
inline float Response::GetExpectedSigma(const CollisionType& collision, const TrackType& track, const o2::track::PID::ID id) const
{
  // use multTPC (legacy behaviour) if multTPC not provided
  return Response::GetExpectedSigmaAtMultiplicity(collision.multTPC(), track, id);
}

/// Gets the expected resolution of the measurement
template <typename TrackType>
inline float Response::GetExpectedSigmaAtMultiplicity(const long multTPC, const TrackType& track, const o2::track::PID::ID id) const
{
  if (!track.hasTPC()) {
    return -999.f;
  }
  float resolution = 0.;
  if (mUseDefaultResolutionParam) {
    const float reso = GetExpectedSignal(track, id) * mResolutionParamsDefault[0] * (static_cast<float>(track.tpcNClsFound()) > 0 ? std::sqrt(1. + mResolutionParamsDefault[1] / static_cast<float>(track.tpcNClsFound())) : 1.f);
    reso >= 0.f ? resolution = reso : resolution = -999.f;
  } else {

    const double ncl = nClNorm / track.tpcNClsFound(); //
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const double dEdx = o2::common::BetheBlochAleph(static_cast<float>(bg), mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow(static_cast<float>(o2::track::pid_constants::sCharges[id]), mChargeFactor);
    const double relReso = GetRelativeResolutiondEdx(p, mass, o2::track::pid_constants::sCharges[id], mResolutionParams[3]);

    const std::vector<double> values{1.f / dEdx, track.tgl(), std::sqrt(ncl), relReso, track.signed1Pt(), multTPC / mMultNormalization};

    const float reso = sqrt(pow(mResolutionParams[0], 2) * values[0] + pow(mResolutionParams[1], 2) * (values[2] * mResolutionParams[5]) * pow(values[0] / sqrt(1 + pow(values[1], 2)), mResolutionParams[2]) + values[2] * pow(values[3], 2) + pow(mResolutionParams[4] * values[4], 2) + pow(values[5] * mResolutionParams[6], 2) + pow(values[5] * (values[0] / sqrt(1 + pow(values[1], 2))) * mResolutionParams[7], 2)) * dEdx * mMIP;
    reso >= 0.f ? resolution = reso : resolution = -999.f;
  }
  return resolution;
}

/// Gets the number of sigma between the actual signal and the expected signal
template <typename CollisionType, typename TrackType>
inline float Response::GetNumberOfSigma(const CollisionType& collision, const TrackType& trk, const o2::track::PID::ID id) const
{
  if (GetExpectedSigma(collision, trk, id) < 0.) {
    return -999.f;
  }
  if (GetExpectedSignal(trk, id) < 0.) {
    return -999.f;
  }
  if (!trk.hasTPC()) {
    return -999.f;
  }
  return ((trk.tpcSignal() - GetExpectedSignal(trk, id)) / GetExpectedSigma(collision, trk, id));
}

template <typename CollisionType, typename TrackType>
inline float Response::GetNumberOfSigmaMCTuned(const CollisionType& collision, const TrackType& trk, const o2::track::PID::ID id, float mcTunedTPCSignal) const
{
  return Response::GetNumberOfSigmaMCTunedAtMultiplicity(collision.multTPC(), trk, id, mcTunedTPCSignal);
}

template <typename TrackType>
inline float Response::GetNumberOfSigmaMCTunedAtMultiplicity(const long multTPC, const TrackType& trk, const o2::track::PID::ID id, float mcTunedTPCSignal) const
{
  if (GetExpectedSigmaAtMultiplicity(multTPC, trk, id) < 0.) {
    return -999.f;
  }
  if (GetExpectedSignal(trk, id) < 0.) {
    return -999.f;
  }
  if (!trk.hasTPC()) {
    return -999.f;
  }
  return ((mcTunedTPCSignal - GetExpectedSignal(trk, id)) / GetExpectedSigmaAtMultiplicity(multTPC, trk, id));
}

/// Gets the deviation between the actual signal and the expected signal
template <typename TrackType>
inline float Response::GetSignalDelta(const TrackType& trk, const o2::track::PID::ID id) const
{
  if (GetExpectedSignal(trk, id) < 0.) {
    return -999.f;
  }
  if (!trk.hasTPC()) {
    return -999.f;
  }
  return (trk.tpcSignal() - GetExpectedSignal(trk, id));
}

//// Gets relative dEdx resolution contribution due relative pt resolution
inline float Response::GetRelativeResolutiondEdx(const float p, const float mass, const float charge, const float resol) const
{
  const float bg = p / mass;
  const float dEdx = o2::common::BetheBlochAleph(bg, mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow(charge, mChargeFactor);
  const float deltaP = resol * std::sqrt(dEdx);
  const float bgDelta = p * (1 + deltaP) / mass;
  const float dEdx2 = o2::common::BetheBlochAleph(bgDelta, mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow(charge, mChargeFactor);
  const float deltaRel = std::abs(dEdx2 - dEdx) / dEdx;
  return deltaRel;
}

inline void Response::PrintAll() const
{
  LOGP(info, "==== TPC PID response parameters: ====");
  for (int i = 0; i < static_cast<int>(mBetheBlochParams.size()); i++) {
    LOGP(info, "BB param [{}] = {}", i, mBetheBlochParams[i]);
  }
  LOGP(info, "use default resolution parametrization = {}", mUseDefaultResolutionParam);
  if (mUseDefaultResolutionParam) {
    LOGP(info, "Default Resolution parametrization: ");
    for (int i = 0; i < static_cast<int>(mResolutionParamsDefault.size()); i++) {
      LOGP(info, "Resolution param [{}] = {}", i, mResolutionParamsDefault[i]);
    }
  } else {
    for (int i = 0; i < static_cast<int>(mResolutionParams.size()); i++) {
      LOGP(info, "Resolution param [{}] = {}", i, mResolutionParams[i]);
    }
  }
  LOGP(info, "mMIP = {}", mMIP);
  LOGP(info, "mChargeFactor = {}", mChargeFactor);
  LOGP(info, "mMultNormalization = {}", mMultNormalization);
  LOGP(info, "nClNorm = {}", nClNorm);
}

} // namespace o2::pid::tpc

#endif // COMMON_CORE_PID_TPCPIDRESPONSE_H_
