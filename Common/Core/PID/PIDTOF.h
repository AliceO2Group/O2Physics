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

#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsTOF/ParameterContainers.h>
#include <Framework/DataTypes.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/PID.h>

#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::pid::tof
{

// Utility values
static constexpr float defaultReturnValue = -999.f; /// Default return value in case TOF measurement is not available

/// \brief Next implementation class to store TOF response parameters for exp. times
class TOFResoParamsV2 : public o2::tof::Parameters<13>
{
 public:
  TOFResoParamsV2() : Parameters(std::array<std::string, 13>{"TrkRes.Pi.P0", "TrkRes.Pi.P1", "TrkRes.Pi.P2", "TrkRes.Pi.P3", "time_resolution",
                                                             "TrkRes.Ka.P0", "TrkRes.Ka.P1", "TrkRes.Ka.P2", "TrkRes.Ka.P3",
                                                             "TrkRes.Pr.P0", "TrkRes.Pr.P1", "TrkRes.Pr.P2", "TrkRes.Pr.P3"},
                                 "TOFResoParamsV2")
  {
    setParameters(std::array<float, 13>{0.008, 0.008, 0.002, 40.0, 60.0,
                                        0.008, 0.008, 0.002, 40.0,
                                        0.008, 0.008, 0.002, 40.0});
  } // Default constructor with default parameters

  ~TOFResoParamsV2() = default;

  template <o2::track::PID::ID pid>
  float getResolution(const float, const float) const
  {
    return -1.f;
  }
  // Momentum shift for charge calibration
  void setMomentumChargeShiftParameters(std::unordered_map<std::string, float> const& pars)
  {
    if (pars.count("Shift.etaN") == 0) { // If the map does not contain the number of eta bins, we assume that no correction has to be applied
      mEtaN = 0;
      return;
    }
    mEtaN = static_cast<int>(pars.at("Shift.etaN"));
    if (mEtaN <= 0) {
      LOG(fatal) << "TOFResoParamsV2 shift: etaN must be positive";
    }
    mEtaStart = pars.at("Shift.etaStart");
    mEtaStop = pars.at("Shift.etaStop");
    if (mEtaStart >= mEtaStop) {
      LOG(fatal) << "TOFResoParamsV2 shift: etaStart must be smaller than etaStop";
    }
    mInvEtaWidth = 1.f / ((mEtaStop - mEtaStart) / mEtaN);
    mContent.clear();
    mContent.resize(mEtaN);
    for (int i = 0; i < mEtaN; ++i) {
      mContent[i] = pars.at(Form("Shift.etaC%i", i));
    }
  }
  // To remove
  void setShiftParameters(std::unordered_map<std::string, float> const& pars)
  {
    setMomentumChargeShiftParameters(pars);
  }

  float getMomentumChargeShift(float eta) const
  {
    if (mEtaN == 0) { // No correction
      // LOG(info) << "TOFResoParamsV2 shift: no correction mEtaN is " << mEtaN;
      return 0.f;
    }
    const int& etaIndex = (eta <= mEtaStart) ? 0 : (eta >= mEtaStop ? (mEtaN - 1) : (eta - mEtaStart) * mInvEtaWidth);
    // LOG(info) << "TOFResoParamsV2 shift: correction for eta " << eta << " is for index " << etaIndex << " = " << shift;
    return mContent.at(etaIndex);
  }
  // To remove
  float getShift(float eta) const
  {
    return getMomentumChargeShift(eta);
  }

  void printMomentumChargeShiftParameters() const
  {
    LOG(info) << "TOF momentum shift parameters";
    LOG(info) << "etaN: " << mEtaN;
    LOG(info) << "etaStart: " << mEtaStart;
    LOG(info) << "etaStop: " << mEtaStop;
    LOG(info) << "content size " << mContent.size();
    for (int i = 0; i < mEtaN; ++i) {
      LOG(info) << "etaC" << i << ": " << mContent[i];
    }
  }
  // To remove
  void printShiftParameters() const
  {
    printMomentumChargeShiftParameters();
  }

  void setTimeShiftParameters(std::string const& filename, std::string const& objname, bool positive)
  {
    TFile f(filename.c_str(), "READ");
    if (f.IsOpen()) {
      if (positive) {
        f.GetObject(objname.c_str(), gPosEtaTimeCorr);
      } else {
        f.GetObject(objname.c_str(), gNegEtaTimeCorr);
      }
      f.Close();
    }
    LOG(info) << "Set the Time Shift parameters from file " << filename << " and object " << objname << " for " << (positive ? "positive" : "negative") << " example of shift at eta 0: " << getTimeShift(0, positive);
  }
  void setTimeShiftParameters(TGraph* g, bool positive)
  {
    if (positive) {
      gPosEtaTimeCorr = g;
    } else {
      gNegEtaTimeCorr = g;
    }
    LOG(info) << "Set the Time Shift parameters from object " << g->GetName() << " " << g->GetTitle() << " for " << (positive ? "positive" : "negative");
  }
  float getTimeShift(float eta, int16_t sign) const
  {
    if (sign > 0) {
      if (!gPosEtaTimeCorr) {
        return 0.f;
      }
      return gPosEtaTimeCorr->Eval(eta);
    }
    if (!gNegEtaTimeCorr) {
      return 0.f;
    }
    return gNegEtaTimeCorr->Eval(eta);
  }

 private:
  // Charge calibration
  int mEtaN = 0; // Number of eta bins, 0 means no correction
  float mEtaStart = 0.f;
  float mEtaStop = 0.f;
  float mInvEtaWidth = 9999.f;
  std::vector<float> mContent;

  // Time shift for post calibration
  TGraph* gPosEtaTimeCorr = nullptr; /// Time shift correction for positive tracks
  TGraph* gNegEtaTimeCorr = nullptr; /// Time shift correction for negative tracks
};

/// \brief Next implementation class to store TOF response parameters for exp. times
class TOFResoParamsV3 : public o2::tof::Parameters<13>
{
 public:
  TOFResoParamsV3() : Parameters(std::array<std::string, 13>{"time_resolution", "time_resolution", "time_resolution", "time_resolution", "time_resolution",
                                                             "time_resolution", "time_resolution", "time_resolution", "time_resolution",
                                                             "time_resolution", "time_resolution", "time_resolution", "time_resolution"},
                                 "TOFResoParamsV3")
  {
    setParameters(std::array<float, 13>{60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0});
  } // Default constructor with default parameters

  ~TOFResoParamsV3() = default;

  // Momentum shift for charge calibration
  void setMomentumChargeShiftParameters(std::unordered_map<std::string, float> const& pars)
  {
    if (pars.count("Shift.etaN") == 0) { // If the map does not contain the number of eta bins, we assume that no correction has to be applied
      mEtaN = 0;
      return;
    }
    mEtaN = static_cast<int>(pars.at("Shift.etaN"));
    if (mEtaN <= 0) {
      LOG(fatal) << "TOFResoParamsV3 shift: etaN must be positive";
    }
    mEtaStart = pars.at("Shift.etaStart");
    mEtaStop = pars.at("Shift.etaStop");
    if (mEtaStart >= mEtaStop) {
      LOG(fatal) << "TOFResoParamsV3 shift: etaStart must be smaller than etaStop";
    }
    mInvEtaWidth = 1.f / ((mEtaStop - mEtaStart) / mEtaN);
    mContent.clear();
    mContent.resize(mEtaN);
    for (int i = 0; i < mEtaN; ++i) {
      mContent[i] = pars.at(Form("Shift.etaC%i", i));
    }
  }

  float getMomentumChargeShift(float eta) const
  {
    if (mEtaN == 0) { // No correction
      // LOG(info) << "TOFResoParamsV3 shift: no correction mEtaN is " << mEtaN;
      return 0.f;
    }
    const int& etaIndex = (eta <= mEtaStart) ? 0 : (eta >= mEtaStop ? (mEtaN - 1) : (eta - mEtaStart) * mInvEtaWidth);
    // LOG(info) << "TOFResoParamsV3 shift: correction for eta " << eta << " is for index " << etaIndex << " = " << shift;
    return mContent.at(etaIndex);
  }

  void printMomentumChargeShiftParameters() const
  {
    LOG(info) << "TOF momentum shift parameters";
    LOG(info) << "etaN: " << mEtaN;
    LOG(info) << "etaStart: " << mEtaStart;
    LOG(info) << "etaStop: " << mEtaStop;
    LOG(info) << "content size " << mContent.size();
    for (int i = 0; i < mEtaN; ++i) {
      LOG(info) << "etaC" << i << ": " << mContent[i];
    }
  }

  // Time shift for post calibration to realign as a function of eta
  void setTimeShiftParameters(std::unordered_map<std::string, float> const& pars, const bool positive);
  void setTimeShiftParameters(std::string const& filename, std::string const& objname, const bool positive);
  void setTimeShiftParameters(TGraph* g, const bool positive);
  float getTimeShift(float eta, int16_t sign) const;

  void printTimeShiftParameters() const
  {
    if (gPosEtaTimeCorr) {
      LOG(info) << "Using a time shift for Pos " << gPosEtaTimeCorr->GetName() << " " << gPosEtaTimeCorr->GetTitle() << " value at 0: " << gPosEtaTimeCorr->Eval(0) << " vs correction " << getTimeShift(0, 1);
    } else {
      LOG(info) << "Using no time shift for Pos vs correction " << getTimeShift(0, 1);
    }
    if (gNegEtaTimeCorr) {
      LOG(info) << "Using a time shift for Neg " << gNegEtaTimeCorr->GetName() << " " << gNegEtaTimeCorr->GetTitle() << " value at 0: " << gNegEtaTimeCorr->Eval(0) << " vs correction " << getTimeShift(0, -1);
    } else {
      LOG(info) << "Using no time shift for Neg vs correction " << getTimeShift(0, -1);
    }
  }

  void setResolutionParametrization(std::unordered_map<std::string, float> const& pars)
  {
    for (int i = 0; i < 9; i++) {
      const std::string baseOpt = Form("tofResTrack.%s_", particleNames[i]);
      // Check if a key begins with a string
      for (const auto& [key, value] : pars) {
        if (key.find(baseOpt) == 0) {
          // Remove from the key the baseOpt
          const std::string fun = key.substr(baseOpt.size());
          if (mResolution[i]) {
            delete mResolution[i];
          }
          mResolution[i] = new TF2(baseOpt.c_str(), fun.c_str(), 0., 20, -1, 1.);
          LOG(info) << "Set the resolution function for " << particleNames[i] << " with formula " << mResolution[i]->GetFormula()->GetExpFormula();
          break;
        }
      }
    }
    // Print a summary
    for (int i = 0; i < 9; i++) {
      if (!mResolution[i]) {
        LOG(info) << "Resolution function for " << particleNames[i] << " not provided, using default " << mDefaultResoParams[i];
        mResolution[i] = new TF2(Form("tofResTrack.%s_Default", particleNames[i]), mDefaultResoParams[i], 0., 20, -1, 1.);
      }
      LOG(info) << "Resolution function for " << particleNames[i] << " is " << mResolution[i]->GetName() << " with formula " << mResolution[i]->GetFormula()->GetExpFormula();
    }
  }

  void setResolutionParametrizationRun2(std::unordered_map<std::string, float> const& pars);

  template <o2::track::PID::ID pid>
  float getResolution(const float p, const float eta) const
  {
    return mResolution[pid]->Eval(p, eta);
  }

  void printResolution() const
  {
    // Print a summary
    for (int i = 0; i < 9; i++) {
      if (!mResolution[i]) {
        LOG(info) << "Resolution function for " << particleNames[i] << " is not defined yet";
        continue;
      }
      LOG(info) << "Resolution function for " << particleNames[i] << " is " << mResolution[i]->GetName() << " with formula " << mResolution[i]->GetFormula()->GetExpFormula();
    }
  }
  void printFullConfig() const
  {
    print();
    printMomentumChargeShiftParameters();
    printTimeShiftParameters();
    printResolution();
  }

 private:
  // Charge calibration
  int mEtaN = 0; // Number of eta bins, 0 means no correction
  float mEtaStart = 0.f;
  float mEtaStop = 0.f;
  float mInvEtaWidth = 9999.f;
  std::vector<float> mContent;
  std::array<TF2*, 9> mResolution{nullptr};
  static constexpr std::array<const char*, 9> mDefaultResoParams{"14.3*TMath::Power((TMath::Max(x-0.319,0.1))*(1-0.4235*y*y),-0.8467)",
                                                                 "14.3*TMath::Power((TMath::Max(x-0.319,0.1))*(1-0.4235*y*y),-0.8467)",
                                                                 "14.3*TMath::Power((TMath::Max(x-0.319,0.1))*(1-0.4235*y*y),-0.8467)",
                                                                 "42.66*TMath::Power((TMath::Max(x-0.417,0.1))*(1-0.4235*y*y),-0.7145)",
                                                                 "99.46*TMath::Power((TMath::Max(x-0.447,0.1))*(1-0.4235*y*y),-0.8094)",
                                                                 "216*TMath::Power((TMath::Max(x-0.647,0.1))*(1-0.4235*y*y),-0.76)",
                                                                 "315*TMath::Power((TMath::Max(x-0.811,0.1))*(1-0.4235*y*y),-0.783)",
                                                                 "157*TMath::Power((TMath::Max(x-0.556,0.1))*(1-0.4235*y*y),-0.783)",
                                                                 "216*TMath::Power((TMath::Max(x-0.647,0.1))*(1-0.4235*y*y),-0.76)"};
  static constexpr std::array<const char*, 9> particleNames = {"El", "Mu", "Pi", "Ka", "Pr", "De", "Tr", "He", "Al"};

  // Time shift for post calibration
  TGraph* gPosEtaTimeCorr = nullptr; /// Time shift correction for positive tracks
  TGraph* gNegEtaTimeCorr = nullptr; /// Time shift correction for negative tracks
};

/// \brief Class to handle the the TOF detector response for the TOF beta measurement
class Beta
{
 public:
  Beta() = default;
  ~Beta() = default;

  /// Computes the beta of a track given a length, a time measurement and an event time (in ps)
  /// \param length Length in cm of the track
  /// \param tofSignal TOF signal in ps for the track
  /// \param collisionTime collision time in ps for the event of the track
  static float GetBeta(const float length, const float tofSignal, const float collisionTime) { return length / (tofSignal - collisionTime) * o2::constants::physics::invLightSpeedCm2PS; }

  /// Gets the beta for the track of interest
  /// \param track Track of interest
  /// \param collisionTime Collision time
  template <typename TrackType>
  static float GetBeta(const TrackType& track, const float collisionTime)
  {
    return track.hasTOF() ? GetBeta(track.length(), track.tofSignal(), collisionTime) : defaultReturnValue;
  }

  /// Gets the beta for the track of interest
  /// \param track Track of interest
  template <typename TrackType>
  static float GetBeta(const TrackType& track)
  {
    return GetBeta(track, track.tofEvTime());
  }

  /// Computes the expected uncertainty on the beta measurement
  /// \param length Length in cm of the track
  /// \param tofSignal TOF signal in ps for the track
  /// \param collisionTime collision time in ps for the event of the track
  /// \param time_reso expected time resolution
  static float GetExpectedSigma(const float length, const float tofSignal, const float collisionTime, const float expectedResolution) { return GetBeta(length, tofSignal, collisionTime) / (tofSignal - collisionTime) * expectedResolution; }

  /// Gets the expected uncertainty on the beta measurement of the track of interest
  /// \param track Track of interest
  template <typename TrackType>
  float GetExpectedSigma(const TrackType& track) const
  {
    return GetExpectedSigma(track.length(), track.tofSignal(), track.tofEvTime(), mExpectedResolution);
  }

  /// Gets the expected beta for a given mass hypothesis (no energy loss taken into account)
  /// \param momentum momentum in GeV/c of the track
  /// \param mass mass in GeV/c2 of the particle of interest
  static float GetExpectedBeta(const float momentum, const float mass) { return momentum > 0 ? momentum / std::sqrt(momentum * momentum + mass * mass) : 0.f; }

  /// Gets the expected beta given the particle index (no energy loss taken into account) of the track of interest
  /// \param track Track of interest
  template <o2::track::PID::ID id, typename TrackType>
  float GetExpectedBeta(const TrackType& track) const
  {
    return GetExpectedBeta(track.p(), o2::track::PID::getMass2Z(id));
  }

  /// Gets the number of sigmas with respect the approximate beta (no energy loss taken into account) of the track of interest
  /// \param track Track of interest
  template <o2::track::PID::ID id, typename TrackType>
  float GetSeparation(const TrackType& track) const
  {
    return (GetBeta(track) - GetExpectedBeta<id>(track)) / GetExpectedSigma(track);
  }

  float mExpectedResolution = 80; /// Expected time resolution
};

/// \brief Class to handle the the TOF detector response for the TOF mass measurement
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
  template <typename TrackType>
  static float GetTOFMass(const TrackType& track, const float beta)
  {
    return track.hasTOF() ? GetTOFMass(track.p(), beta) : defaultReturnValue;
  }

  /// Gets the TOF mass for the track of interest
  /// \param track Track of interest
  template <typename TrackType>
  static float GetTOFMass(const TrackType& track)
  {
    return track.hasTOF() ? GetTOFMass(track.p(), Beta::GetBeta<TrackType>(track)) : defaultReturnValue;
  }
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
  static float ComputeExpectedTime(const float tofExpMom, const float length) { return length * std::sqrt((mMassZSqared) + (tofExpMom * tofExpMom)) / (o2::constants::physics::LightSpeedCm2PS * tofExpMom); }

  /// Gets the expected signal of the track of interest under the PID assumption
  /// \param track Track of interest
  static float GetExpectedSignal(const TrackType& track)
  {
    if (!track.hasTOF()) {
      return defaultReturnValue;
    }
    if (track.trackType() == o2::aod::track::Run2Track) {
      return ComputeExpectedTime(track.tofExpMom() * o2::constants::physics::invLightSpeedCm2PS, track.length());
    }
    return ComputeExpectedTime(track.tofExpMom(), track.length());
  }

  /// Gets the expected signal of the track of interest under the PID assumption corrected for shifts in expected momentum
  /// \param parameters Parameters to correct for the momentum shift
  /// \param track Track of interest
  template <typename ParamType>
  static float GetCorrectedExpectedSignal(const ParamType& parameters, const TrackType& track)
  {
    if (!track.hasTOF()) {
      return defaultReturnValue;
    }
    if (track.trackType() == o2::aod::track::Run2Track) {
      return ComputeExpectedTime(track.tofExpMom() * o2::constants::physics::invLightSpeedCm2PS / (1.f + track.sign() * parameters.getMomentumChargeShift(track.eta())), track.length());
    }
    LOG(debug) << "TOF exp. mom. " << track.tofExpMom() << " shifted = " << track.tofExpMom() / (1.f + track.sign() * parameters.getMomentumChargeShift(track.eta()));
    return ComputeExpectedTime(track.tofExpMom() / (1.f + track.sign() * parameters.getMomentumChargeShift(track.eta())), track.length()) + parameters.getTimeShift(track.eta(), track.sign());
  }

  /// Gets the expected resolution of the t-texp-t0
  /// Given a TOF signal and collision time resolutions
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  /// \param tofSignal TOF signal of the track of interest
  /// \param collisionTimeRes Collision time resolution of the track of interest
  template <typename ParamType>
  static float GetExpectedSigma(const ParamType& parameters, const TrackType& track, const float tofSignal, const float collisionTimeRes)
  {
    const float& mom = track.p();
    const float& etaTrack = track.eta();
    if (mom <= 0) {
      return -999.f;
    }
    const float reso = parameters.template getResolution<id>(mom, etaTrack);
    if (reso > 0) {
      return std::sqrt(reso * reso + parameters[4] * parameters[4] + collisionTimeRes * collisionTimeRes);
    }
    if constexpr (id <= o2::track::PID::Pion) {
      LOG(debug) << "Using parameters for the pion hypothesis and ID " << id;
      const float dpp = parameters[0] + parameters[1] * mom + parameters[2] * mMassZ / mom; // mean relative pt resolution;
      const float sigma = dpp * tofSignal / (1. + mom * mom / (mMassZSqared));
      return std::sqrt(sigma * sigma + parameters[3] * parameters[3] / mom / mom + parameters[4] * parameters[4] + collisionTimeRes * collisionTimeRes);
    } else if constexpr (id == o2::track::PID::Kaon) {
      LOG(debug) << "Using parameters for the kaon hypothesis and ID " << id;
      const float dpp = parameters[5] + parameters[6] * mom + parameters[7] * mMassZ / mom; // mean relative pt resolution;
      const float sigma = dpp * tofSignal / (1. + mom * mom / (mMassZSqared));
      return std::sqrt(sigma * sigma + parameters[8] * parameters[8] / mom / mom + parameters[4] * parameters[4] + collisionTimeRes * collisionTimeRes);
    }
    LOG(debug) << "Using parameters for the proton hypothesis and ID " << id;
    const float dpp = parameters[9] + parameters[10] * mom + parameters[11] * mMassZ / mom; // mean relative pt resolution;
    const float sigma = dpp * tofSignal / (1. + mom * mom / (mMassZSqared));
    return std::sqrt(sigma * sigma + parameters[12] * parameters[12] / mom / mom + parameters[4] * parameters[4] + collisionTimeRes * collisionTimeRes);
  }

  /// Gets the expected resolution of the t-texp-t0
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  template <typename ParamType>
  static float GetExpectedSigma(const ParamType& parameters, const TrackType& track)
  {
    return GetExpectedSigma(parameters, track, track.tofSignal(), track.tofEvTimeErr());
  }

  /// Gets the expected resolution of the time measurement, uses the expected time and no event time resolution
  /// \param parameters Parameters to use to compute the expected resolution
  /// \param track Track of interest
  template <typename ParamType>
  static float GetExpectedSigmaTracking(const ParamType& parameters, const TrackType& track)
  {
    return GetExpectedSigma(parameters, track, GetCorrectedExpectedSignal(parameters, track), 0.f);
  }

  /// Gets the number of sigmas with respect the expected time
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  /// \param collisionTime Collision time
  /// \param collisionTimeRes Collision time resolution of the track of interest
  template <typename ParamType>
  static float GetSeparation(const ParamType& parameters, const TrackType& track, const float collisionTime, const float resolution)
  {
    return track.hasTOF() ? (track.tofSignal() - collisionTime - GetCorrectedExpectedSignal(parameters, track)) / resolution : defaultReturnValue;
  }

  /// Gets the number of sigmas with respect the expected time
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  template <typename ParamType>
  static float GetSeparation(const ParamType& parameters, const TrackType& track, const float resolution)
  {
    return GetSeparation(parameters, track, track.tofEvTime(), resolution);
  }

  /// Gets the number of sigmas with respect the expected time
  /// \param parameters Detector response parameters
  /// \param track Track of interest
  template <typename ParamType>
  static float GetSeparation(const ParamType& parameters, const TrackType& track)
  {
    return GetSeparation(parameters, track, track.tofEvTime(), GetExpectedSigma(parameters, track));
  }
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
  static float GetTOFSignalForParticleHypothesis(const TrackType& track, const ResponseImplementation<pid>& response)
  {
    return ComputeTOFSignal(track.trackTime(), response.GetExpectedSignal(track));
  }

  /// Returns the expected time of a track considering the PID hypothesis used in tracking
  /// \param track input track
  static float GetTOFSignal(const TrackType& track)
  {
    switch (track.pidForTracking()) {
      case 0:
        return GetTOFSignalForParticleHypothesis(track, responseEl);
      case 1:
        return GetTOFSignalForParticleHypothesis(track, responseMu);
      case 2:
        return GetTOFSignalForParticleHypothesis(track, responsePi);
      case 3:
        return GetTOFSignalForParticleHypothesis(track, responseKa);
      case 4:
        return GetTOFSignalForParticleHypothesis(track, responsePr);
      case 5:
        return GetTOFSignalForParticleHypothesis(track, responseDe);
      case 6:
        return GetTOFSignalForParticleHypothesis(track, responseTr);
      case 7:
        return GetTOFSignalForParticleHypothesis(track, responseHe);
      case 8:
        return GetTOFSignalForParticleHypothesis(track, responseAl);
      default:
        return 0.f;
        break;
    }
  }
};

} // namespace o2::pid::tof

#endif // COMMON_CORE_PID_PIDTOF_H_
