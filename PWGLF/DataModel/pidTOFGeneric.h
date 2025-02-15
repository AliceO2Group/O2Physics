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

#ifndef PWGLF_DATAMODEL_PIDTOFGENERIC_H_
#define PWGLF_DATAMODEL_PIDTOFGENERIC_H_
#include "CommonDataFormat/InteractionRecord.h"
#include "Common/Core/PID/PIDTOF.h"

namespace o2::aod
{
namespace evtime
{

DECLARE_SOA_COLUMN(EvTime, evTime, float);             //! Event time. Can be obtained via a combination of detectors e.g. TOF, FT0A, FT0C
DECLARE_SOA_COLUMN(EvTimeErr, evTimeErr, float);       //! Error of event time. Can be obtained via a combination of detectors e.g. TOF, FT0A, FT0C
DECLARE_SOA_COLUMN(EvTimeTOF, evTimeTOF, float);       //! Event time computed with the TOF detector
DECLARE_SOA_COLUMN(EvTimeTOFErr, evTimeTOFErr, float); //! Error of the event time computed with the TOF detector
DECLARE_SOA_COLUMN(EvTimeFT0, evTimeFT0, float);       //! Event time computed with the FT0 detector
DECLARE_SOA_COLUMN(EvTimeFT0Err, evTimeFT0Err, float); //! Error of the event time computed with the FT0 detector
} // namespace evtime

DECLARE_SOA_TABLE(EvTimeTOFFT0, "AOD", "EvTimeTOFFT0", //! Table of the event time. One entry per collision.
                  evtime::EvTime,
                  evtime::EvTimeErr,
                  evtime::EvTimeTOF,
                  evtime::EvTimeTOFErr,
                  evtime::EvTimeFT0,
                  evtime::EvTimeFT0Err);

namespace tracktime
{

DECLARE_SOA_COLUMN(EvTimeForTrack, evTimeForTrack, float);       //! Event time. Removed the bias for the specific track
DECLARE_SOA_COLUMN(EvTimeErrForTrack, evTimeErrForTrack, float); //! Error of event time. Removed the bias for the specific track
} // namespace tracktime

DECLARE_SOA_TABLE(EvTimeTOFFT0ForTrack, "AOD", "EvTimeForTrack", //! Table of the event time. One entry per track.
                  tracktime::EvTimeForTrack,
                  tracktime::EvTimeErrForTrack);

namespace pidtofgeneric
{

static constexpr float kCSPEED = TMath::C() * 1.0e2f * 1.0e-12f; // c in cm/ps

template <typename TTrack>
class TofPidNewCollision
{
 public:
  TofPidNewCollision() = default;
  ~TofPidNewCollision() = default;

  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;
  o2::track::PID::ID pidType;

  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<TTrack, pid>;
  static constexpr auto responseEl = ResponseImplementation<o2::track::PID::Electron>();
  static constexpr auto responseMu = ResponseImplementation<o2::track::PID::Muon>();
  static constexpr auto responsePi = ResponseImplementation<o2::track::PID::Pion>();
  static constexpr auto responseKa = ResponseImplementation<o2::track::PID::Kaon>();
  static constexpr auto responsePr = ResponseImplementation<o2::track::PID::Proton>();
  static constexpr auto responseDe = ResponseImplementation<o2::track::PID::Deuteron>();
  static constexpr auto responseTr = ResponseImplementation<o2::track::PID::Triton>();
  static constexpr auto responseHe = ResponseImplementation<o2::track::PID::Helium3>();
  static constexpr auto responseAl = ResponseImplementation<o2::track::PID::Alpha>();

  void SetParams(o2::pid::tof::TOFResoParamsV2 const& para)
  {
    mRespParamsV2.setParameters(para);
  }

  void SetPidType(o2::track::PID::ID pidId)
  {
    pidType = pidId;
  }

  template <typename TCollision>
  float GetTOFNSigma(o2::track::PID::ID pidId, TTrack const& track, TCollision const& originalcol, TCollision const& correctedcol, bool EnableBCAO2D = true);

  template <typename TCollision>
  float GetTOFNSigma(TTrack const& track, TCollision const& originalcol, TCollision const& correctedcol, bool EnableBCAO2D = true);

  float GetTOFNSigma(TTrack const& track);
  float GetTOFNSigma(o2::track::PID::ID pidId, TTrack const& track);

  float CalculateTOFNSigma(o2::track::PID::ID pidId, TTrack const& track, double tofsignal, double evTime, double evTimeErr)
  {

    float expSigma, tofNsigma = -999;

    switch (pidId) {
      case 0:
        expSigma = responseEl.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseEl.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      case 1:
        expSigma = responseMu.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseMu.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      case 2:
        expSigma = responsePi.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responsePi.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      case 3:
        expSigma = responseKa.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseKa.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      case 4:
        expSigma = responsePr.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responsePr.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      case 5:
        expSigma = responseDe.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseDe.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      case 6:
        expSigma = responseTr.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseTr.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      case 7:
        expSigma = responseHe.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseHe.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      case 8:
        expSigma = responseAl.GetExpectedSigma(mRespParamsV2, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseAl.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
        break;
      default:
        LOG(fatal) << "Wrong particle ID in TofPidSecondary class";
        return -999;
    }

    return tofNsigma;
  }
};

template <typename TTrack>
template <typename TCollision>
float TofPidNewCollision<TTrack>::GetTOFNSigma(o2::track::PID::ID pidId, TTrack const& track, TCollision const& originalcol, TCollision const& correctedcol, bool EnableBCAO2D)
{

  if (!track.has_collision() || !track.hasTOF()) {
    return -999;
  }

  float mMassHyp = o2::track::pid_constants::sMasses2Z[track.pidForTracking()];
  float expTime = track.length() * sqrt((mMassHyp * mMassHyp) + (track.tofExpMom() * track.tofExpMom())) / (kCSPEED * track.tofExpMom()); // L*E/(p*c) = L/v

  float evTime = correctedcol.evTime();
  float evTimeErr = correctedcol.evTimeErr();
  float tofsignal = track.trackTime() * 1000 + expTime; // in ps

  if (originalcol.globalIndex() == correctedcol.globalIndex()) {
    evTime = track.evTimeForTrack();
    evTimeErr = track.evTimeErrForTrack();
  } else {
    if (EnableBCAO2D) {
      auto originalbc = originalcol.template bc_as<o2::aod::BCsWithTimestamps>();
      auto correctedbc = correctedcol.template bc_as<o2::aod::BCsWithTimestamps>();
      o2::InteractionRecord originalIR = o2::InteractionRecord::long2IR(originalbc.globalBC());
      o2::InteractionRecord correctedIR = o2::InteractionRecord::long2IR(correctedbc.globalBC());
      tofsignal += originalIR.differenceInBCNS(correctedIR) * 1000;
    } else {
      auto originalbc = originalcol.template foundBC_as<o2::aod::BCsWithTimestamps>();
      auto correctedbc = correctedcol.template foundBC_as<o2::aod::BCsWithTimestamps>();
      o2::InteractionRecord originalIR = o2::InteractionRecord::long2IR(originalbc.globalBC());
      o2::InteractionRecord correctedIR = o2::InteractionRecord::long2IR(correctedbc.globalBC());
      tofsignal += originalIR.differenceInBCNS(correctedIR) * 1000;
    }
  }

  float tofNsigma = CalculateTOFNSigma(pidId, track, tofsignal, evTime, evTimeErr);
  return tofNsigma;
}

template <typename TTrack>
template <typename TCollision>
float TofPidNewCollision<TTrack>::GetTOFNSigma(TTrack const& track, TCollision const& originalcol, TCollision const& correctedcol, bool EnableBCAO2D)
{
  return GetTOFNSigma(pidType, track, originalcol, correctedcol, EnableBCAO2D);
}

template <typename TTrack>
float TofPidNewCollision<TTrack>::GetTOFNSigma(o2::track::PID::ID pidId, TTrack const& track)
{

  if (!track.has_collision() || !track.hasTOF()) {
    return -999;
  }

  float mMassHyp = o2::track::pid_constants::sMasses2Z[track.pidForTracking()];
  float expTime = track.length() * sqrt((mMassHyp * mMassHyp) + (track.tofExpMom() * track.tofExpMom())) / (kCSPEED * track.tofExpMom()); // L*E/(p*c) = L/v

  float evTime = track.evTimeForTrack();
  float evTimeErr = track.evTimeErrForTrack();
  float tofsignal = track.trackTime() * 1000 + expTime; // in ps

  float tofNsigma = CalculateTOFNSigma(pidId, track, tofsignal, evTime, evTimeErr);
  return tofNsigma;
}

template <typename TTrack>
float TofPidNewCollision<TTrack>::GetTOFNSigma(TTrack const& track)
{
  return GetTOFNSigma(pidType, track);
}

} // namespace pidtofgeneric
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_PIDTOFGENERIC_H_
