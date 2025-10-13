// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file closePairRejection.h
/// \brief Definition of ClosePairRejection class
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_CLOSEPAIRREJECTION_H_
#define PWGCF_FEMTO_CORE_CLOSEPAIRREJECTION_H_

#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/histManager.h"

#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace o2::analysis::femto
{
namespace closepairrejection
{
// enum for track histograms
enum CprHist {
  // kinemtics
  kAverage,
  kRadius0,
  kRadius1,
  kRadius2,
  kRadius3,
  kRadius4,
  kRadius5,
  kRadius6,
  kRadius7,
  kRadius8,
  kCprHistogramLast
};

struct ConfCpr : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("ClosePairRejection");
  o2::framework::Configurable<bool> on{"on", true, "Turn on CPR"};
  o2::framework::Configurable<float> detaMax{"detaMax", 0.01f, "Maximium deta"};
  o2::framework::Configurable<float> dphistarMax{"dphistarMax", 0.01f, "Maximum dphistar"};
  o2::framework::ConfigurableAxis binningDeta{"binningDeta", {{500, -0.5, 0.5}}, "deta"};
  o2::framework::ConfigurableAxis binningDphistar{"binningDphistar", {{500, -0.5, 0.5}}, "dphi"};
};

// tpc radii for computing phistar
constexpr int kNradii = 9;
constexpr std::array<float, kNradii> kTpcRadius = {85., 105., 125., 145., 165., 185., 205., 225., 245.}; // in cm

// directory names
constexpr char PrefixTrackTrackSe[] = "CPR_TrackTrack/SE/";
constexpr char PrefixTrackTrackMe[] = "CPR_TrackTrack/ME/";
constexpr char PrefixTrackV0Se[] = "CPR_TrackV0Daughter/SE/";
constexpr char PrefixTrackV0Me[] = "CPR_TrackV0Daughter/ME/";
constexpr char PrefixV0V0PosSe[] = "CPR_V0V0_PosDau/SE/";
constexpr char PrefixV0V0NegSe[] = "CPR_V0V0_NegDau/SE/";
constexpr char PrefixV0V0PosMe[] = "CPR_V0V0_PosDau/ME/";
constexpr char PrefixV0V0NegMe[] = "CPR_V0V0_NegDau/ME/";
constexpr char PrefixTrackTwoTrackResonanceSe[] = "CPR_TrackResonanceDaughter/SE/";
constexpr char PrefixTrackTwoTrackResonnaceMe[] = "CPR_TrackResonanceDaughter/ME/";
constexpr char PrefixTrackCascadeSe[] = "CPR_TrackCascadeBachelor/SE/";
constexpr char PrefixTrackCascadeMe[] = "CPR_TrackCascadeBachelor/ME/";
constexpr char PrefixTrackKinkSe[] = "CPR_TrackKink/SE/";
constexpr char PrefixTrackKinkMe[] = "CPR_TrackKink/ME/";

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<CprHist>, kCprHistogramLast> HistTable = {
  {{kAverage, o2::framework::kTH2F, "hAverage", "#Delta #eta vs #Delta #phi* (averaged over all radii); #Delta #eta; #Delta #phi*"},
   {kRadius0, o2::framework::kTH2F, "hRadius0", "Radius 0: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius1, o2::framework::kTH2F, "hRadius1", "Radius 1: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius2, o2::framework::kTH2F, "hRadius2", "Radius 2: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius3, o2::framework::kTH2F, "hRadius3", "Radius 3: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius4, o2::framework::kTH2F, "hRadius4", "Radius 4: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius5, o2::framework::kTH2F, "hRadius5", "Radius 5: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius6, o2::framework::kTH2F, "hRadius6", "Radius 6: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius7, o2::framework::kTH2F, "hRadius7", "Radius 7: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
   {kRadius8, o2::framework::kTH2F, "hRadius8", "Radius 8: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"}}};

template <typename T>
auto makeCprHistSpecMap(const T& confCpr)
{
  return std::map<CprHist, std::vector<framework::AxisSpec>>{
    {kAverage, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius0, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius1, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius2, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius3, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius4, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius5, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius6, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius7, {confCpr.binningDeta, confCpr.binningDphistar}},
    {kRadius8, {confCpr.binningDeta, confCpr.binningDphistar}}};
};

template <const char* prefix>
class CloseTrackRejection
{
 public:
  CloseTrackRejection() = default;
  virtual ~CloseTrackRejection() = default;

  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, int chargeAbsTrack1, int chargeAbsTrack2)
  {
    mDetaMax = detaMax;
    mDphistarMax = dphistarMax;

    mChargeAbsTrack1 = chargeAbsTrack1;
    mChargeAbsTrack2 = chargeAbsTrack2;

    mHistogramRegistry = registry;

    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kAverage, HistTable), getHistDesc(kAverage, HistTable), getHistType(kAverage, HistTable), {specs.at(kAverage)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius0, HistTable), getHistDesc(kRadius0, HistTable), getHistType(kRadius0, HistTable), {specs.at(kRadius0)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius1, HistTable), getHistDesc(kRadius1, HistTable), getHistType(kRadius1, HistTable), {specs.at(kRadius1)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius2, HistTable), getHistDesc(kRadius2, HistTable), getHistType(kRadius2, HistTable), {specs.at(kRadius2)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius3, HistTable), getHistDesc(kRadius3, HistTable), getHistType(kRadius3, HistTable), {specs.at(kRadius3)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius4, HistTable), getHistDesc(kRadius4, HistTable), getHistType(kRadius4, HistTable), {specs.at(kRadius4)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius5, HistTable), getHistDesc(kRadius5, HistTable), getHistType(kRadius5, HistTable), {specs.at(kRadius5)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius6, HistTable), getHistDesc(kRadius6, HistTable), getHistType(kRadius6, HistTable), {specs.at(kRadius6)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius7, HistTable), getHistDesc(kRadius7, HistTable), getHistType(kRadius7, HistTable), {specs.at(kRadius7)});
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kRadius8, HistTable), getHistDesc(kRadius8, HistTable), getHistType(kRadius8, HistTable), {specs.at(kRadius8)});
  }

  void setMagField(float magField) { mMagField = magField; }

  template <typename T1, typename T2>
  void compute(const T1& track1, const T2& track2)
  {
    // reset values
    mAverageDphistar = 0.f;
    mDeta = 0.f;
    mDphistar.fill(0.f);

    mDeta = track1.eta() - track2.eta();
    for (size_t i = 0; i < kTpcRadius.size(); i++) {
      auto phistar1 = utils::dphistar(mMagField, kTpcRadius[i], mChargeAbsTrack1 * track1.signedPt(), track1.phi());
      auto phistar2 = utils::dphistar(mMagField, kTpcRadius[i], mChargeAbsTrack2 * track2.signedPt(), track2.phi());
      if (phistar1 && phistar2) {
        // if the calculation for one phistar fails, keep the default value, which is 0
        // this makes it more likelier for the pair to be rejected sind the averave will be biased towards lower values
        mDphistar.at(i) = phistar1.value() - phistar2.value();
      }
    }
    mAverageDphistar = std::accumulate(mDphistar.begin(), mDphistar.end(), 0.f) / mDphistar.size();
  }

  void fill()
  {
    // fill average hist
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kAverage, HistTable)), mDeta, mAverageDphistar);

    // fill radii hists
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius0, HistTable)), mDeta, mDphistar.at(0));
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius1, HistTable)), mDeta, mDphistar.at(1));
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius2, HistTable)), mDeta, mDphistar.at(2));
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius3, HistTable)), mDeta, mDphistar.at(3));
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius4, HistTable)), mDeta, mDphistar.at(4));
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius5, HistTable)), mDeta, mDphistar.at(5));
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius6, HistTable)), mDeta, mDphistar.at(6));
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius7, HistTable)), mDeta, mDphistar.at(7));
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius8, HistTable)), mDeta, mDphistar.at(8));
  }

  bool isClosePair() const
  {
    return std::hypot(mAverageDphistar / mDphistarMax, mDeta / mDetaMax) < 1.f;
  }

 private:
  int mChargeAbsTrack1 = 0;
  int mChargeAbsTrack2 = 0;
  float mMagField = 0.f;
  float mAverageDphistar = 0.f;
  float mDeta = 0.f;
  float mDetaMax = 0.f;
  float mDphistarMax = 0.f;
  std::array<float, kNradii> mDphistar = {0.f};

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
};

template <const char* prefix>
class ClosePairRejectionTrackTrack
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, int absChargeTrack1, int absChargeTrack2, bool isActivated)
  {
    mIsActivated = isActivated;
    if (mIsActivated) {
      mCtr.init(registry, specs, detaMax, dphistarMax, absChargeTrack1, absChargeTrack2);
    }
  }

  void setMagField(float magField) { mCtr.setMagField(magField); }
  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& track1, T2 const& track2, T3 const& /*tracks*/)
  {
    mCtr.compute(track1, track2);
  }
  bool isClosePair() const { return mCtr.isClosePair(); }
  void fill() { mCtr.fill(); }
  bool isActivated() const { return mIsActivated; }

 private:
  CloseTrackRejection<prefix> mCtr;
  bool mIsActivated = true;
};

template <const char* prefixPosDaus, const char* prefixNegDaus>
class ClosePairRejectionV0V0
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, bool isActivated)
  {
    mIsActivated = isActivated;
    if (mIsActivated) {
      mCtrPos.init(registry, specs, detaMax, dphistarMax, 1, 1);
      mCtrNeg.init(registry, specs, detaMax, dphistarMax, 1, 1);
    }
  }

  void setMagField(float magField)
  {
    mCtrPos.setMagField(magField);
    mCtrNeg.setMagField(magField);
  }
  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& v01, T2 const& v02, T3 const& /*tracks*/)
  {
    auto posDau1 = v01.template posDau_as<T3>();
    auto negDau1 = v01.template posDau_as<T3>();

    auto posDau2 = v02.template posDau_as<T3>();
    auto negDau2 = v02.template posDau_as<T3>();

    mCtrPos.compute(posDau1, posDau2);
    mCtrNeg.compute(negDau1, negDau2);
  }
  bool isClosePair() const { return mCtrPos.isClosePair() && mCtrNeg.isClosePair(); }
  void fill()
  {
    mCtrPos.fill();
    mCtrNeg.fill();
  }
  bool isActivated() const { return mIsActivated; }

 private:
  CloseTrackRejection<prefixPosDaus> mCtrPos;
  CloseTrackRejection<prefixNegDaus> mCtrNeg;
  bool mIsActivated = true;
};

template <const char* prefix>
class ClosePairRejectionTrackV0 // can also be used for any particle type that has pos/neg daughters, like resonances
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, int absChargeTrack, bool isActivated)
  {
    mIsActivated = isActivated;
    // initialize CPR with charge of the track and the same charge for the daughter particle
    // absolute charge of the daughter track will be 1, so we just pass 1
    if (mIsActivated) {
      mCtr.init(registry, specs, detaMax, dphistarMax, absChargeTrack, 1);
    }
  }

  void setMagField(float magField)
  {
    mCtr.setMagField(magField);
  }
  template <typename T1, typename T2, typename T3>
  void setPair(const T1& track, const T2& v0, const T3 /*trackTable*/)
  {
    if (track.signedPt() > 0) {
      auto daughter = v0.template posDau_as<T3>();
      mCtr.compute(track, daughter);
    } else {
      auto daughter = v0.template negDau_as<T3>();
      mCtr.compute(track, daughter);
    }
  }

  bool isClosePair() const { return mCtr.isClosePair(); }
  void fill()
  {
    mCtr.fill();
  }
  bool isActivated() const { return mIsActivated; }

 private:
  CloseTrackRejection<prefix> mCtr;
  bool mIsActivated = true;
};

template <const char* prefix>
class ClosePairRejectionTrackCascade
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, int absChargeTrack, bool isActivated)
  {
    mIsActivated = isActivated;
    if (mIsActivated) {
      // charge of cascade is always 1
      mCtr.init(registry, specs, detaMax, dphistarMax, absChargeTrack, 1);
    }
  }

  void setMagField(float magField)
  {
    mCtr.setMagField(magField);
  }
  template <typename T1, typename T2, typename T3>
  void setPair(const T1& track, const T2& cascade, const T3 /*trackTable*/)
  {
    auto bachelor = cascade.template posDau_as<T3>();
    mCtr.compute(track, bachelor);
  }

  bool isClosePair() const { return mCtr.isClosePair(); }
  void fill()
  {
    mCtr.fill();
  }
  bool isActivated() const { return mIsActivated; }

 private:
  CloseTrackRejection<prefix> mCtr;
  bool mIsActivated = true;
};

template <const char* prefix>
class ClosePairRejectionTrackKink
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, int absChargeTrack, bool isActivated)
  {
    mIsActivated = isActivated;
    // The charged daughter has absolute charge of 1, so we can pass 1 directly
    if (mIsActivated) {
      mCtr.init(registry, specs, detaMax, dphistarMax, absChargeTrack, 1);
    }
  }

  void setMagField(float magField)
  {
    mCtr.setMagField(magField);
  }

  template <typename T1, typename T2, typename T3>
  void setPair(const T1& track, const T2& kink, const T3 /*trackTable*/)
  {
    auto daughter = kink.template chaDau_as<T3>();
    mCtr.compute(track, daughter);
  }

  bool isClosePair() const { return mCtr.isClosePair(); }
  void fill()
  {
    mCtr.fill();
  }
  bool isActivated() const { return mIsActivated; }

 private:
  CloseTrackRejection<prefix> mCtr;
  bool mIsActivated = true;
};

}; // namespace closepairrejection
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_CLOSEPAIRREJECTION_H_
