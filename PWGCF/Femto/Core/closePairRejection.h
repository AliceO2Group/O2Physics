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

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/HistogramRegistry.h"

#include <array>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_map>
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
  o2::framework::ConfigurableAxis binningDeta{"binningDeta", {{200, -0.2, 0.2}}, "deta"};
  o2::framework::ConfigurableAxis binningDphistar{"binningDphistar", {{200, -0.2, 0.2}}, "dphi"};
};

// tpc radii for computing phistar
constexpr int kNradii = 9;
constexpr std::array<float, kNradii> kTpcRadius = {85., 105., 125., 145., 165., 185., 205., 225., 245.};

// directory names
constexpr char PrefixTrackTrackSe[] = "CPR_TrackTrack/SE/";
constexpr char PrefixTrackTrackMe[] = "CPR_TrackTrack/ME/";
constexpr char PrefixTrackPosDauSe[] = "CPR_TrackPosDau/SE/";
constexpr char PrefixTrackNegDauSe[] = "CPR_TrackNegDau/SE/";
constexpr char PrefixTrackPosDauMe[] = "CPR_TrackPosDau/ME/";
constexpr char PrefixTrackNegDauMe[] = "CPR_TrackNegDau/ME/";

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
  std::map<CprHist, std::vector<framework::AxisSpec>> specs;
  for (int i = 0; i < kCprHistogramLast; ++i) {
    specs[static_cast<CprHist>(i)] = {confCpr.binningDeta, confCpr.binningDphistar};
  }
  return specs;
};

template <const char* prefix>
class CloseTrackRejection
{
 public:
  CloseTrackRejection() {}

  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax)
  {
    mHistogramRegistry = registry;
    mDetaMax = detaMax;
    mDphistarMax = dphistarMax;

    for (int i = 0; i < kCprHistogramLast; ++i) {
      auto hist = static_cast<CprHist>(i);
      mHistogramRegistry->add(std::string(prefix) + GetHistNamev2(hist, HistTable), GetHistDesc(hist, HistTable), GetHistType(hist, HistTable), {specs.at(hist)});
    }
  }

  void setMagField(float magField) { mMagField = magField; }

  void reset()
  {
    mSameCharge = false;
    mAverageDphistar = 0.f;
    mDeta = 0.f;
    mDphistar.fill(0.f);
  }

  template <typename T1, typename T2>
  void compute(const T1& track1, const T2& track2)
  {
    reset();
    if (track1.sign() != track2.sign()) {
      return;
    }
    mSameCharge = true;
    mDeta = track1.eta() - track2.eta();
    for (size_t i = 0; i < kTpcRadius.size(); ++i) {
      auto phistar1 = utils::dphistar(mMagField, kTpcRadius[i], track1.sign(), track1.pt(), track1.phi());
      auto phistar2 = utils::dphistar(mMagField, kTpcRadius[i], track2.sign(), track2.pt(), track2.phi());

      if (phistar1 && phistar2) {
        // if the calculation for one phistar fails, keep the default value, which is 0
        // this makes it more likelier for the pair to be rejected sind the averave will be biased towards lower values
        mDphistar[i] = phistar1.value() - phistar2.value();
      }
    }
    mAverageDphistar = std::accumulate(mDphistar.begin(), mDphistar.end(), 0.f) / mDphistar.size();
  }

  void fill()
  {
    // fill average hist
    mHistogramRegistry->fill(HIST(prefix) + HIST(GetHistName(kAverage, HistTable)), mDeta, mAverageDphistar);

    // fill radii hists
    for (int i = 0; i < kNradii; ++i) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(GetHistName(static_cast<CprHist>(kRadius0 + 1), HistTable)), mDeta, mDphistar.at(i));
    }
  }

  bool isClosePair() const
  {
    if (!mSameCharge) {
      return false;
    } else {
      return ((mAverageDphistar * mAverageDphistar) / (mDphistarMax * mDphistarMax) + (mDeta * mDeta) / (mDetaMax * mDetaMax)) < 1.f;
    }
  }

 private:
  float mMagField = 0.f;
  float mAverageDphistar = 0.f;
  bool mSameCharge = false;
  float mDeta = 0.f;
  float mDetaMax = 0.f;
  float mDphistarMax = 0.f;
  std::array<float, kNradii> mDphistar = {};

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
};

template <const char* prefix>
class ClosePairRejectionTrackTrack
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, bool isActivated)
  {
    mIsActivated = isActivated;
    mCtr.init(registry, specs, detaMax, dphistarMax);
  }

  void setMagField(float magField) { mCtr.setMagField(magField); }
  template <typename T1, typename T2>
  void setPair(const T1& track1, const T2& track2)
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

template <const char* prefixPosDau, const char* prefixNegDau>
class ClosePairRejectionTrackV0 // can also be used for any particle type that has pos/neg daughters, like resonances
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, bool isActivated)
  {
    mIsActivated = isActivated;
    mCtrPosDau.init(registry, specs, detaMax, dphistarMax);
    mCtrNegDau.init(registry, specs, detaMax, dphistarMax);
  }

  void setMagField(float magField)
  {
    mCtrPosDau.setMagField(magField);
    mCtrNegDau.setMagField(magField);
  }
  template <typename T1, typename T2, typename T3>
  void setPair(const T1& track, const T2& v0, const T3 /*trackTable*/)
  {
    auto posDaughter = v0.template posDau_as<T3>();
    auto negDaughter = v0.template negDau_as<T3>();
    mCtrPosDau.compute(track, posDaughter);
    mCtrNegDau.compute(track, negDaughter);
  }
  bool isClosePair() const { return mCtrPosDau.isClosePair() || mCtrNegDau.isClosePair(); }
  void fill()
  {
    mCtrPosDau.fill();
    mCtrNegDau.fill();
  }
  bool isActivated() const { return mIsActivated; }

 private:
  CloseTrackRejection<prefixPosDau> mCtrPosDau;
  CloseTrackRejection<prefixNegDau> mCtrNegDau;
  bool mIsActivated = true;
};

}; // namespace closepairrejection
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_CLOSEPAIRREJECTION_H_
