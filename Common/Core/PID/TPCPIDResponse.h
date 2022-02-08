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
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "TFormula.h"
#include <iostream>

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
  void SetResolutionParams(const std::vector<double>& resolutionParams, const int readoutchamber = 0) { mResolutionParams[readoutchamber] = resolutionParams; };
  void SetMIP(const float mip) { mMIP = mip; };
  void SetChargeFactor(const float chargeFactor) { mChargeFactor = chargeFactor; };
  void SetMultiplicityNormalization(const float multNormalization) { mMultNormalization = multNormalization; };
  void SetMaxNumberClusters(const float maxNumClusters) { mMaxClusters = maxNumClusters; };
  void SetResolutionParametrization(TFormula* sigmaParametrization) { fSigmaParametrization = sigmaParametrization; };
  void SetUseDefaultResolutionParam(const bool useDefault) { useDefaultResolutionParam = useDefault; };
  void SetReadOutChamber(const int readoutchamber) { mReadOutChamber = readoutchamber; };

  const std::array<float, 5> GetBetheBlochParams() const { return mBetheBlochParams; };
  const std::vector<double> GetResolutionParams() const { return mResolutionParams[mReadOutChamber]; };
  const float GetMIP() const { return mMIP; };
  const float GetChargeFactor() const { return mChargeFactor; };
  const float GetMultiplicityNormalization() const { return mMultNormalization; };
  const float GetMaxNumberClusters() const { return mMaxClusters; };
  TFormula* GetResolutionParametrization() { return fSigmaParametrization; };
  const int GetReadOutChamber() const
  {
    if (mReadOutChamber == 0) {
      LOGP(info, "Readout chamber Global");
    } else if (mReadOutChamber == 1) {
      LOGP(info, "Readout chamber IROC");
    } else if (mReadOutChamber == 2) {
      LOGP(info, "Readout chamber OROC1");
    } else if (mReadOutChamber == 3) {
      LOGP(info, "Readout chamber OROC2");
    } else if (mReadOutChamber == 4) {
      LOGP(info, "Readout chamber OROC3");
    }
    return mReadOutChamber;
  };

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
  std::array<float, 5> mBetheBlochParams = {0.0320981, 19.9768, 2.52666e-16, 2.72123, 6.08092};
  std::array<float, 2> mResolutionParamsDefault = {0.07, 0.0};
  std::array<std::vector<double>, 5> mResolutionParams = {{{5.43799e-7, 0.053044, 0.667584, 0.0142667, 0.00235175, 1.22482, 2.3501e-7, 0.031585}, {5.43799e-7, 0.053044, 0.667584, 0.0142667, 0.00235175, 1.22482, 2.3501e-7, 0.031585}, {5.43799e-7, 0.053044, 0.667584, 0.0142667, 0.00235175, 1.22482, 2.3501e-7, 0.031585}, {5.43799e-7, 0.053044, 0.667584, 0.0142667, 0.00235175, 1.22482, 2.3501e-7, 0.031585}, {5.43799e-7, 0.053044, 0.667584, 0.0142667, 0.00235175, 1.22482, 2.3501e-7, 0.031585}}};
  /// Choose between global parameters and different readout chambers:
  /// 0 = Global, 1 = IROC, 2 = OROC1, 3 = OROC2, 4 = OROC3
  int mReadOutChamber = 0;
  float mMIP = 50.f;
  float mChargeFactor = 2.3f;
  float mMultNormalization = 11000.;
  float mMaxClusters = 63.;
  bool useDefaultResolutionParam = true;
  TFormula* fSigmaParametrization = new TFormula("fSigmaParametrization", "sqrt(([0]**2)*x[0]+(([1]**2)*(x[2]*[5])*(x[0]/sqrt(1+x[1]**2))**[2])+x[2]*x[3]**2+([4]*x[4])**2 +((x[5]*[6])**2)+(x[5]*(x[0]/sqrt(1+x[1]**2))*[7])**2)");

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
  // if (useDefaultResolutionParam == true){
  // std::cout << "Use Default " << std::endl;
  const float reso = track.tpcSignal() * mResolutionParamsDefault[0] * ((float)track.tpcNClsFound() > 0 ? std::sqrt(1. + mResolutionParamsDefault[1] / (float)track.tpcNClsFound()) : 1.f);
  return reso >= 0.f ? reso : 0.f;
  //}
  // else{

  //   std::vector<double> values;
  //   const double ncl = track.tpcNClsFound();
  //   const double p = track.tpcInnerParam();
  //   const double mass = o2::track::pid_constants::sMasses[id];
  //   const double bg =  p/mass;
  //   const double dEdx = mMIP * o2::tpc::BetheBlochAleph((float)bg, mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4]) * std::pow((float)o2::track::pid_constants::sCharges[id], mChargeFactor);
  //   const double relReso = o2::pid::tpc::Response::GetRelativeResolutiondEdx(p,mass,o2::track::pid_constants::sCharges[id],mResolutionParams[mReadOutChamber][3]);

  //   values.push_back((1./dEdx));
  //   values.push_back((track.tgl()));
  //   values.push_back((std::sqrt(mMaxClusters/ncl)));
  //   values.push_back(relReso);
  //   values.push_back((track.signed1Pt()));
  //   values.push_back((collision.multTPC()) / mMultNormalization);

  //   fSigmaParametrization->SetParameters(mResolutionParams[mReadOutChamber].data());

  //   return fSigmaParametrization->EvalPar(values.data())*(dEdx);
  // }
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
  const float deltaRel = TMath::Abs(dEdx2 - dEdx) / dEdx;
  return deltaRel;
}

inline void Response::PrintAll() const
{
  LOGP(info, "==== TPC PID response parameters: ====");
  for (int i = 0; i < int(mBetheBlochParams.size()); i++)
    LOGP(info, "BB param [{}] = {}", i, mBetheBlochParams[i]);
  if (useDefaultResolutionParam) {
    LOGP(info, "Default Resolution parametrization: ");
    for (int i = 0; i < int(mResolutionParamsDefault.size()); i++)
      LOGP(info, "Resolution param [{}] = {}", i, mResolutionParamsDefault[i]);
  } else {
    for (int i = 0; i < int(mResolutionParams[mReadOutChamber].size()); i++)
      LOGP(info, "Resolution param [{}] = {}", i, mResolutionParams[mReadOutChamber][i]);
  }
  switch (mReadOutChamber) {
    case 0:
      LOGP(info, "No specific readout chamber selected. Use Global resolution parameters.");
      break;
    case 1:
      LOGP(info, "Use resolution parameters for IROC.");
      break;
    case 2:
      LOGP(info, "Use resolution parameters for OROC1.");
      break;
    case 3:
      LOGP(info, "Use resolution parameters for OROC2.");
      break;
    case 4:
      LOGP(info, "Use resolution parameters for OROC3.");
      break;
    default:
      break;
  }
  LOGP(info, "mMIP = {}", mMIP);
  LOGP(info, "mChargeFactor = {}", mChargeFactor);
  LOGP(info, "mMultNormalization = {}", mMultNormalization);
  LOGP(info, "mMaxClusters = {}", mMaxClusters);
}

} // namespace o2::pid::tpc

#endif // O2_FRAMEWORK_PIDRESPONSE_H_