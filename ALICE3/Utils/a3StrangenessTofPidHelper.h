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
/// \file a3StrangenessTofPidHelper.h
/// \brief Helper to calculate expected arrival time for decay products of weakly decaying particles
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>
///

#ifndef ALICE3_UTILS_UTILSSTRANGENESSTOFPID_H_
#define ALICE3_UTILS_UTILSSTRANGENESSTOFPID_H_

#include "ALICE3/Core/TrackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <TRandom3.h>

#include <array>
#include <cmath>

namespace o2::upgrade::stratofpid
{

enum class Topology { V0,
                      Cascade };

enum class TrackType { Positive,
                       Negative,
                       V0,
                       Bachelor,
                       Cascade };

enum class CandidateType { K0S,
                           Lambda,
                           AntiLambda,
                           Xi,
                           AntiXi,
                           Omega,
                           AntiOmega };

struct TrackTofResult {
  bool hasInnerTof, hasOuterTof;
  float nSigmaInner, nSigmaOuter;
  float expectedTimeInner, expectedTimeOuter;
  float measuredTimeInner, measuredTimeOuter;
};

template <Topology T>
struct StrangenessTofResults;

template <>
struct StrangenessTofResults<Topology::V0> {
  TrackTofResult pos;
  TrackTofResult neg;
};

template <>
struct StrangenessTofResults<Topology::Cascade> {
  TrackTofResult pos;
  TrackTofResult neg;
  TrackTofResult bach;
};

template <CandidateType TType>
struct CandidateInfo;

template <>
struct CandidateInfo<CandidateType::K0S> {
  static constexpr float V0Mass = o2::constants::physics::MassKaonNeutral;
  static constexpr float PosMass = o2::constants::physics::MassPionCharged;
  static constexpr float NegMass = o2::constants::physics::MassPionCharged;
};
template <>
struct CandidateInfo<CandidateType::Lambda> {
  static constexpr float V0Mass = o2::constants::physics::MassLambda0;
  static constexpr float PosMass = o2::constants::physics::MassProton;
  static constexpr float NegMass = o2::constants::physics::MassPionCharged;
};
template <>
struct CandidateInfo<CandidateType::AntiLambda> {
  static constexpr float V0Mass = o2::constants::physics::MassLambda0;
  static constexpr float PosMass = o2::constants::physics::MassPionCharged;
  static constexpr float NegMass = o2::constants::physics::MassProton;
};
template <>
struct CandidateInfo<CandidateType::Xi> {
  static constexpr float V0Mass = o2::constants::physics::MassLambda0;
  static constexpr float PosMass = o2::constants::physics::MassProton;
  static constexpr float NegMass = o2::constants::physics::MassPionCharged;
  static constexpr float BachMass = o2::constants::physics::MassPionCharged;
  static constexpr float CascMass = o2::constants::physics::MassXiMinus;
};
template <>
struct CandidateInfo<CandidateType::AntiXi> {
  static constexpr float V0Mass = o2::constants::physics::MassLambda0;
  static constexpr float PosMass = o2::constants::physics::MassPionCharged;
  static constexpr float NegMass = o2::constants::physics::MassProton;
  static constexpr float BachMass = o2::constants::physics::MassPionCharged;
  static constexpr float CascMass = o2::constants::physics::MassXiMinus;
};
template <>
struct CandidateInfo<CandidateType::Omega> {
  static constexpr float V0Mass = o2::constants::physics::MassLambda0;
  static constexpr float PosMass = o2::constants::physics::MassProton;
  static constexpr float NegMass = o2::constants::physics::MassPionCharged;
  static constexpr float BachMass = o2::constants::physics::MassKaonCharged;
  static constexpr float CascMass = o2::constants::physics::MassOmegaMinus;
};
template <>
struct CandidateInfo<CandidateType::AntiOmega> {
  static constexpr float V0Mass = o2::constants::physics::MassLambda0;
  static constexpr float PosMass = o2::constants::physics::MassPionCharged;
  static constexpr float NegMass = o2::constants::physics::MassProton;
  static constexpr float BachMass = o2::constants::physics::MassKaonCharged;
  static constexpr float CascMass = o2::constants::physics::MassOmegaMinus;
};

class StrangenessTofPidBase
{
 public:
  StrangenessTofPidBase() : mRand(0) {}
  void setRandomSeed(const int seed) { mRand.SetSeed(seed); }
  void setMagneticField(const float b) { mMagneticField = b; }
  void setInnerResolution(const float res) { mInnerResolution = res; }
  void setOuterResolution(const float res) { mOuterResolution = res; }
  void setResolution(const float ires, const float ores)
  {
    setInnerResolution(ires);
    setOuterResolution(ores);
  }
  void setResolution(const float res)
  {
    setInnerResolution(res);
    setOuterResolution(res);
  }
  void setInnerRadius(const float radius) { mInnerRadius = radius; }
  void setOuterRadius(const float radius) { mOuterRadius = radius; }
  void setRadius(const float iradius, const float oradius)
  {
    setInnerRadius(iradius);
    setOuterRadius(oradius);
  }
  void setTrack(const TrackType type, const o2::track::TrackParCov& track, const float vt = -1.f, const float mass = -1.f)
  {
    switch (type) {
      case TrackType::Positive:
        mPosTrack = track;
        mPosTrackTime = vt * NanoToPico;
        mTruePosCandMass = mass;
        break;
      case TrackType::Negative:
        mNegTrack = track;
        mNegTrackTime = vt * NanoToPico;
        mTrueNegCandMass = mass;
        break;
      case TrackType::V0:
        mV0Track = track;
        break;

      default:
        break;
    }
  }

  bool propagateTrackToOriginalRadius(o2::track::TrackParCov& track, const std::array<float, 3> originalVtx) const
  {
    static constexpr float Tolerance = 1e-2;
    std::array<float, 3> trackXYZ{};
    track.getXYZGlo(trackXYZ);

    // Check if track is already in the right place
    if (std::abs(std::hypot(trackXYZ[0], trackXYZ[1]) - std::hypot(originalVtx[0], originalVtx[1])) < Tolerance) {
      return true;
    }

    float targetX = 1e+3;
    if (!track.getXatLabR(std::hypot(originalVtx[0], originalVtx[1]), targetX, mMagneticField)) {
      return false;
    }

    return track.propagateTo(targetX, mMagneticField);
  }

 [[nodiscard]] float calculateArrivalTime(const o2::track::TrackParCov& track, const std::array<float, 3> vtx, const float radius, const float mass) const
  {
    o2::track::TrackParCov propagatedTrack(track);
    if (!propagateTrackToOriginalRadius(propagatedTrack, vtx)) {
      return FailedToPropagate;
    }

    const float trackLength = o2::upgrade::computeTrackLength(propagatedTrack, radius, mMagneticField);
    if (trackLength < 0) {
      return FailedToPropagate;
    }

    const float velocity = o2::upgrade::computeParticleVelocity(propagatedTrack.getP(), mass);
    return trackLength / velocity;
  }

  // Special case for V0s
  [[nodiscard]] float calculateArrivalTimeV0(const o2::track::TrackParCov& track, const std::array<float, 3> pv, const std::array<float, 3> sv, const float mass) const
  {
    o2::track::TrackParCov propagatedTrack(track);
    if (!propagateTrackToOriginalRadius(propagatedTrack, pv)) {
      return FailedToPropagate;
    }

    const float dx = pv[0] - sv[0];
    const float dy = pv[1] - sv[1];
    const float dz = pv[2] - sv[2];
    const float trackLength = std::hypot(dx, dy, dz);
    const float velocity = o2::upgrade::computeParticleVelocity(propagatedTrack.getP(), mass);
    return trackLength / velocity;
  }

  float smear(const float time, const float sigma)
  {
    return mRand.Gaus(time, sigma);
  }

  template <CandidateType TCand>
  StrangenessTofResults<Topology::V0> findNSigmas(const std::array<float, 3> vtx, const std::array<float, 3> vtxV0, const float timeOffset = 0)
  {
    const float v0ExpectedTime = timeOffset + calculateArrivalTimeV0(mV0Track, vtx, vtxV0, CandidateInfo<TCand>::V0Mass);
    const float posExpectedTimeInner = v0ExpectedTime + calculateArrivalTime(mPosTrack, vtxV0, mInnerRadius, CandidateInfo<TCand>::PosMass);
    const float posExpectedTimeOuter = v0ExpectedTime + calculateArrivalTime(mPosTrack, vtxV0, mOuterRadius, CandidateInfo<TCand>::PosMass);
    const float negExpectedTimeInner = v0ExpectedTime + calculateArrivalTime(mNegTrack, vtxV0, mInnerRadius, CandidateInfo<TCand>::NegMass);
    const float negExpectedTimeOuter = v0ExpectedTime + calculateArrivalTime(mNegTrack, vtxV0, mOuterRadius, CandidateInfo<TCand>::NegMass);

    if (!mMeasuredTimeSmeared) {
      mPosMeasuredTimeInner = smear(mPosTrackTime + calculateArrivalTime(mPosTrack, vtxV0, mInnerRadius, mTruePosCandMass), mInnerResolution);
      mPosMeasuredTimeOuter = smear(mPosTrackTime + calculateArrivalTime(mPosTrack, vtxV0, mOuterRadius, mTruePosCandMass), mOuterResolution);
      mNegMeasuredTimeInner = smear(mNegTrackTime + calculateArrivalTime(mNegTrack, vtxV0, mInnerRadius, mTrueNegCandMass), mInnerResolution);
      mNegMeasuredTimeOuter = smear(mNegTrackTime + calculateArrivalTime(mNegTrack, vtxV0, mOuterRadius, mTrueNegCandMass), mOuterResolution);
      mMeasuredTimeSmeared = true;
    }

    return StrangenessTofResults<Topology::V0>{
      .pos = {
        .hasInnerTof = (mPosMeasuredTimeInner > 0 && posExpectedTimeInner > 0),
        .hasOuterTof = (mPosMeasuredTimeOuter > 0 && posExpectedTimeOuter > 0),
        .nSigmaInner = (mPosMeasuredTimeInner - posExpectedTimeInner) / mInnerResolution,
        .nSigmaOuter = (mPosMeasuredTimeOuter - posExpectedTimeOuter) / mOuterResolution,
        .expectedTimeInner = posExpectedTimeInner * PicoToNano,
        .expectedTimeOuter = posExpectedTimeOuter * PicoToNano,
        .measuredTimeInner = mPosMeasuredTimeInner * PicoToNano,
        .measuredTimeOuter = mPosMeasuredTimeOuter * PicoToNano,
      },
      .neg = {
        .hasInnerTof = (mNegMeasuredTimeInner > 0 && negExpectedTimeInner > 0),
        .hasOuterTof = (mNegMeasuredTimeOuter > 0 && negExpectedTimeOuter > 0),
        .nSigmaInner = (mNegMeasuredTimeInner - negExpectedTimeInner) / mInnerResolution,
        .nSigmaOuter = (mNegMeasuredTimeOuter - negExpectedTimeOuter) / mOuterResolution,
        .expectedTimeInner = negExpectedTimeInner * PicoToNano,
        .expectedTimeOuter = negExpectedTimeOuter * PicoToNano,
        .measuredTimeInner = mNegMeasuredTimeInner * PicoToNano,
        .measuredTimeOuter = mNegMeasuredTimeOuter * PicoToNano,
      }};
  }

  void reset()
  {
    mMeasuredTimeSmeared = false;
    mPosTrackTime = 0.f;
    mNegTrackTime = 0.f;
    mPosMeasuredTimeInner = 0.f;
    mPosMeasuredTimeOuter = 0.f;
    mNegMeasuredTimeInner = 0.f;
    mNegMeasuredTimeOuter = 0.f;
    mTruePosCandMass = 0.f;
    mTrueNegCandMass = 0.f;
    mPosTrack = o2::track::TrackParCov();
    mNegTrack = o2::track::TrackParCov();
    mV0Track = o2::track::TrackParCov();
  }

 protected:
  static constexpr float NanoToPico = 1e+3;
  static constexpr float PicoToNano = 1e-3;
  static constexpr float FailedToPropagate = -1e+8;

  // Detector config
  float mInnerResolution = 20.f; // ps
  float mOuterResolution = 20.f; // ps
  float mInnerRadius = 21.f;     // cm
  float mOuterRadius = 92.f;     // cm
  float mMagneticField = 20.f;   // kG

  // Track properties
  bool mMeasuredTimeSmeared = false;
  float mPosTrackTime{}, mNegTrackTime{};
  float mPosMeasuredTimeInner{}, mPosMeasuredTimeOuter{};
  float mNegMeasuredTimeInner{}, mNegMeasuredTimeOuter{};
  float mTruePosCandMass{}, mTrueNegCandMass{};
  o2::track::TrackParCov mPosTrack, mNegTrack, mV0Track;

  // For smearing
  TRandom3 mRand;
};

template <Topology T>
class StrangenessTofPid;

template <>
class StrangenessTofPid<Topology::V0> : public StrangenessTofPidBase
{
};

template <>
class StrangenessTofPid<Topology::Cascade> : public StrangenessTofPidBase
{
 public:
  void setTrack(const TrackType type, const o2::track::TrackParCov& track, const float vt = -1.f, const float mass = -1.f)
  {
    switch (type) {
      case TrackType::Bachelor:
        mBachTrack = track;
        mBachTrackTime = vt * NanoToPico;
        mTrueBachCandMass = mass;
        break;
      case TrackType::Cascade:
        mCascTrack = track;
        break;
      default:
        StrangenessTofPidBase::setTrack(type, track, vt, mass);
        break;
    }
  }

  template <CandidateType TCand>
  StrangenessTofResults<Topology::Cascade> findNSigmas(const std::array<float, 3> vtx, const std::array<float, 3> vtxCascade, const std::array<float, 3> vtxV0)
  {
    const float cascExpectedTime = calculateArrivalTime(mCascTrack, vtx, std::hypot(vtxCascade[0], vtxCascade[1]), CandidateInfo<TCand>::CascMass);
    const float bachExpectedTimeInner = cascExpectedTime + calculateArrivalTime(mBachTrack, vtxCascade, mInnerRadius, CandidateInfo<TCand>::BachMass);
    const float bachExpectedTimeOuter = cascExpectedTime + calculateArrivalTime(mBachTrack, vtxCascade, mOuterRadius, CandidateInfo<TCand>::BachMass);
  
    if (!mMeasuredTimeSmeared) {
      mBachMeasuredTimeInner = smear(mBachTrackTime + calculateArrivalTime(mBachTrack, vtxCascade, mInnerRadius, mTrueBachCandMass), mInnerResolution);
      mBachMeasuredTimeOuter = smear(mBachTrackTime + calculateArrivalTime(mBachTrack, vtxCascade, mOuterRadius, mTrueBachCandMass), mOuterResolution);
    }

    const auto straTofResultsV0 = StrangenessTofPidBase::findNSigmas<TCand>(vtxCascade, vtxV0, cascExpectedTime);
    return StrangenessTofResults<Topology::Cascade>{
      .pos = straTofResultsV0.pos,
      .neg = straTofResultsV0.neg,
      .bach = {
        .hasInnerTof = (mBachMeasuredTimeInner > 0 && bachExpectedTimeInner > 0),
        .hasOuterTof = (mBachMeasuredTimeOuter > 0 && bachExpectedTimeOuter > 0),
        .nSigmaInner = (mBachMeasuredTimeInner - bachExpectedTimeInner) / mInnerResolution,
        .nSigmaOuter = (mBachMeasuredTimeOuter - bachExpectedTimeOuter) / mOuterResolution,
        .expectedTimeInner = bachExpectedTimeInner * PicoToNano,
        .expectedTimeOuter = bachExpectedTimeOuter * PicoToNano,
        .measuredTimeInner = mBachMeasuredTimeInner * PicoToNano,
        .measuredTimeOuter = mBachMeasuredTimeOuter * PicoToNano,
      }};
  }

  void reset()
  {
    mBachTrackTime = 0.f;
    mTrueBachCandMass = 0.f;
    mBachMeasuredTimeInner = 0.f;
    mBachMeasuredTimeOuter = 0.f;
    mBachTrack = o2::track::TrackParCov();
    mCascTrack = o2::track::TrackParCov();
    StrangenessTofPidBase::reset();
  }

 protected:
  // Track properties
  float mBachTrackTime{}, mTrueBachCandMass{};
  float mBachMeasuredTimeInner{}, mBachMeasuredTimeOuter{};
  o2::track::TrackParCov mBachTrack, mCascTrack;
};

} // namespace o2::upgrade::stratofpid
#endif // ALICE3_UTILS_UTILSSTRANGENESSTOFPID_H_
