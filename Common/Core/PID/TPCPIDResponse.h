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

// ROOT includes
#include "Rtypes.h"
#include "TMath.h"
#include "TNamed.h"

// O2 includes
#include "ReconstructionDataFormats/PID.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::pid::tpc
{

/// \brief Class to handle the TPC PID response
template <typename TrackType, o2::track::PID::ID id>
class TPCPIDResponse
{

 public:
  TPCPIDResponse() = default;
  ~TPCPIDResponse() = default;

  /// Setter and Getter for the private parameters
  void SetBetheBlochParams(const std::array<float, 5>& betheBlochParams) { mBetheBlochParams = betheBlochParams; };
  void SetResolutionParams(const std::array<float, 2>& resolutionParams) { mResolutionParams = resolutionParams; };
  void SetMIP(const float mip) { mMIP = mip; };
  void SetChargeFactor(const float chargeFactor) { mChargeFactor = chargeFactor; };
  /// A: I would either rename or remove this function. From the name it is not clear what it does. The ones before already set everything
  void Set(const std::array<float, 5>& betheBlochParameters,
           const std::array<float, 2>& resolutionParameters,
           const float chargeFactor = 2.3f,
           const float mip = 50.f)
  {
    mBetheBlochParams = betheBlochParameters;
    mResolutionParams = resolutionParameters;
    mMIP = mip;
    mChargeFactor = chargeFactor;
  };

  const std::array<float, 5> GetBetheBlochParams() const { return mBetheBlochParams; };
  const std::array<float, 2> GetResolutionParams() const { return mResolutionParams; };
  const float GetMIP() const { return mMIP; };
  const float GetChargeFactor() const { return mChargeFactor; };

  /////////// Actual calculation on the tracks //////////////

  /// Empirical ALEPH parameterization of the Bethe-Bloch formula, normalized to 1 at the minimum.
  template <typename T>
  static T BetheBlochAleph(T bg, T kp1, T kp2, T kp3, T kp4, T kp5);

  /// Gets the expected signal of the track
  float GetExpectedSignal(const TrackType& track) const;
  /// Gets the expected resolution of the track
  float GetExpectedSigma(const TrackType& trk) const;
  /// Gets the number of sigmas with respect the expected value
  float GetSeparation(const TrackType& trk) const;

  /// Save to ccdb
  /// Get object from ccdb

 private:
  std::array<float, 5> mBetheBlochParams = {0.0320981, 19.9768, 2.52666e-16, 2.72123, 6.08092};
  std::array<float, 2> mResolutionParams = {0.07, 0.0};
  float mMIP = 50.f;
  float mChargeFactor = 2.3f;

  ClassDef(TPCPIDResponse, 1);

}; // class TPCPIDResponse

} // namespace o2::pid::tpc

#endif // O2_FRAMEWORK_PIDRESPONSE_H_