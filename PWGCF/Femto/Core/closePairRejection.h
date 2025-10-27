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

// default config, cpr between two charged tracks (or track vs particle which decays into one charged track, like sigma)
struct ConfCpr : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("ClosePairRejection");
  o2::framework::Configurable<bool> on{"on", true, "Turn on CPR"};
  o2::framework::Configurable<bool> plotAllRadii{"plotAllRadii", false, "Plot deta-dphi distribution at all radii"};
  o2::framework::Configurable<bool> plotAverage{"plotAverage", true, "Plot average deta dphi distribution"};
  o2::framework::Configurable<float> detaMax{"detaMax", 0.01f, "Maximium deta"};
  o2::framework::Configurable<float> dphistarMax{"dphistarMax", 0.01f, "Maximum dphistar"};
  o2::framework::ConfigurableAxis binningDeta{"binningDeta", {{500, -0.5, 0.5}}, "deta"};
  o2::framework::ConfigurableAxis binningDphistar{"binningDphistar", {{500, -0.5, 0.5}}, "dphi"};
};

struct ConfCprTrackV0 : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("ClosePairRejectionTrackV0");
  o2::framework::Configurable<bool> onSameCharge{"onSameCharge", true, "Turn on CPR for track and same charge daughter"};
  o2::framework::Configurable<bool> onOppositeCharge{"onOppositeCharge", false, "Turn on CPR for track and opposite charge daughter"};
  o2::framework::Configurable<bool> plotAllRadii{"plotAllRadii", false, "Plot deta-dphi distribution at all radii"};
  o2::framework::Configurable<bool> plotAverage{"plotAverage", true, "Plot average deta dphi distribution"};
  o2::framework::Configurable<float> detaMaxSameCharge{"detaMaxSameCharge", 0.01f, "Maximium deta between track and same charge daughter"};
  o2::framework::Configurable<float> dphistarMaxSameCharge{"dphistarMaxSameCharge", 0.01f, "Maximium dphistar between track and same charge daughter"};
  o2::framework::Configurable<float> detaMaxOppositeCharge{"detaMaxOppositeCharge", 0.01f, "Maximum deta between track and opposite charge daughter"};
  o2::framework::Configurable<float> dphistarMaxOppositeCharge{"dphistarMaxOppositeCharge", 0.01f, "Maximum dphistar between track and opposite charge daughter"};
  o2::framework::ConfigurableAxis binningDeta{"binningDeta", {{500, -0.5, 0.5}}, "deta"};
  o2::framework::ConfigurableAxis binningDphistar{"binningDphistar", {{500, -0.5, 0.5}}, "dphi"};
};

struct ConfCprV0V0 : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("ClosePairRejectionV0V0");
  o2::framework::Configurable<bool> on{"on", true, "Turn on CPR"};
  o2::framework::Configurable<bool> plotAllRadii{"plotAllRadii", false, "Plot deta-dphi distribution at all radii"};
  o2::framework::Configurable<bool> plotAverage{"plotAverage", true, "Plot average deta dphi distribution"};
  o2::framework::Configurable<float> detaMaxPosDau{"detaMaxPosDau", 0.01f, "Maximium deta between positive daughters"};
  o2::framework::Configurable<float> dphistarMaxPosDau{"dphistarMaxPosDau", 0.01f, "Maximium dphistar between positive daughters"};
  o2::framework::Configurable<float> detaMaxNegDau{"detaMaxNegDau", 0.01f, "Maximum deta between negative daughters"};
  o2::framework::Configurable<float> dphistarMaxNegDau{"dphistarMaxNegDau", 0.01f, "Maximum dphistar between negative daughters"};
  o2::framework::ConfigurableAxis binningDeta{"binningDeta", {{500, -0.5, 0.5}}, "deta"};
  o2::framework::ConfigurableAxis binningDphistar{"binningDphistar", {{500, -0.5, 0.5}}, "dphi"};
};

struct ConfCprTrrackCascade : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("ClosePairRejectionTrackCascade");
  o2::framework::Configurable<bool> onBachelor{"onBachelor", true, "Turn on CPR for track and bachelor"};
  o2::framework::Configurable<bool> onSameCharge{"onSameCharge", false, "Turn on CPR for track and same charge V0 daughter"};
  o2::framework::Configurable<bool> onOppositeCharge{"onOppositeCharge", false, "Turn on CPR for track and opposite charge V0 daughter"};
  o2::framework::Configurable<bool> plotAllRadii{"plotAllRadii", false, "Plot deta-dphi distribution at all radii"};
  o2::framework::Configurable<bool> plotAverage{"plotAverage", true, "Plot average deta dphi distribution"};
  o2::framework::Configurable<float> detaMaxBachelor{"detaMaxBachelor", 0.01f, "Maximium deta between track and bachelor"};
  o2::framework::Configurable<float> dphistarMaxBachelor{"dphistarMaxBachelor", 0.01f, "Maximium dphistar between track and bachelor"};
  o2::framework::Configurable<float> detaMaxSameCharge{"detaMaxSameCharge", 0.01f, "Maximium deta between track and same charge daughter"};
  o2::framework::Configurable<float> dphistarMaxSameCharge{"dphistarMaxSameCharge", 0.01f, "Maximium dphistar between track and same charge daughter"};
  o2::framework::Configurable<float> detaMaxOppositeCharge{"detaMaxOppositeCharge", 0.01f, "Maximum deta between track and opposite charge daughter"};
  o2::framework::Configurable<float> dphistarMaxOppositeCharge{"dphistarMaxOppositeCharge", 0.01f, "Maximum dphistar between track and opposite charge daughter"};
  o2::framework::ConfigurableAxis binningDeta{"binningDeta", {{500, -0.5, 0.5}}, "deta"};
  o2::framework::ConfigurableAxis binningDphistar{"binningDphistar", {{500, -0.5, 0.5}}, "dphi"};
  o2::framework::Configurable<bool> on{"on", true, "Turn on CPR"};
};

// tpc radii for computing phistar
constexpr int Nradii = 9;
constexpr std::array<float, Nradii> TpcRadii = {85., 105., 125., 145., 165., 185., 205., 225., 245.}; // in cm

// directory names
constexpr char PrefixTrackTrackSe[] = "CPR_TrackTrack/SE/";
constexpr char PrefixTrackTrackMe[] = "CPR_TrackTrack/ME/";
constexpr char PrefixTrackV0SameChargeSe[] = "CPR_TrackV0DauSameCharge/SE/";
constexpr char PrefixTrackV0SameChargeMe[] = "CPR_TrackV0DauSameCharge/ME/";
constexpr char PrefixTrackV0OppositeChargeSe[] = "CPR_TrackV0DauOppositeCharge/SE/";
constexpr char PrefixTrackV0OppositeChargeMe[] = "CPR_TrackV0DauOppositeCharge/ME/";
constexpr char PrefixV0V0PosSe[] = "CPR_V0V0_PosDau/SE/";
constexpr char PrefixV0V0NegSe[] = "CPR_V0V0_NegDau/SE/";
constexpr char PrefixV0V0PosMe[] = "CPR_V0V0_PosDau/ME/";
constexpr char PrefixV0V0NegMe[] = "CPR_V0V0_NegDau/ME/";
constexpr char PrefixTrackTwoTrackResonanceSe[] = "CPR_TrackResonanceDaughter/SE/";
constexpr char PrefixTrackTwoTrackResonnaceMe[] = "CPR_TrackResonanceDaughter/ME/";
constexpr char PrefixTrackCascadeBachelorSe[] = "CPR_TrackCascadeBachelor/SE/";
constexpr char PrefixTrackCascadeBachelorMe[] = "CPR_TrackCascadeBachelor/ME/";
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
  ~CloseTrackRejection() = default;

  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            bool plotAverage,
            bool plotAllRadii,
            float detaMax,
            float dphistarMax,
            int chargeAbsTrack1,
            int chargeAbsTrack2)
  {
    mDetaMax = detaMax;
    mDphistarMax = dphistarMax;

    // check the limits
    if (mDetaMax <= 0 || mDphistarMax <= 0) {
      LOG(warn) << "Close Pair Rejection configured with 0 or negative limits. Histograms will be filled, but no CPR cut will be applied!";
      mPlotOnly = true;
    } else {
      mPlotOnly = false;
    }

    mChargeAbsTrack1 = chargeAbsTrack1;
    mChargeAbsTrack2 = chargeAbsTrack2;

    mHistogramRegistry = registry;

    mPlotAverage = plotAverage;
    mPlotAllRadii = plotAllRadii;

    if (mPlotAverage) {
      mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kAverage, HistTable), getHistDesc(kAverage, HistTable), getHistType(kAverage, HistTable), {specs.at(kAverage)});
    }
    if (mPlotAllRadii) {
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
  }

  void setMagField(float magField) { mMagField = magField; }

  template <typename T1, typename T2>
  void compute(T1 const& track1, T2 const& track2)
  {
    // reset values
    mAverageDphistar = 0.f;
    mDeta = 0.f;
    mDphistar.fill(0.f);

    mDeta = track1.eta() - track2.eta();
    for (size_t i = 0; i < TpcRadii.size(); i++) {
      auto phistar1 = utils::dphistar(mMagField, TpcRadii[i], mChargeAbsTrack1 * track1.signedPt(), track1.phi());
      auto phistar2 = utils::dphistar(mMagField, TpcRadii[i], mChargeAbsTrack2 * track2.signedPt(), track2.phi());
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
    if (mPlotAverage) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kAverage, HistTable)), mDeta, mAverageDphistar);
    }

    // fill radii hists
    if (mPlotAllRadii) {
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
  }

  bool isClosePair() const
  {
    return !mPlotOnly && std::hypot(mAverageDphistar / mDphistarMax, mDeta / mDetaMax) < 1.f;
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mPlotAllRadii = false;
  bool mPlotAverage = true;

  int mChargeAbsTrack1 = 0;
  int mChargeAbsTrack2 = 0;
  float mMagField = 0.f;
  float mDetaMax = 0.f;
  float mDphistarMax = 0.f;

  float mAverageDphistar = 0.f;
  float mDeta = 0.f;
  std::array<float, Nradii> mDphistar = {0.f};

  bool mPlotOnly = true;
};

template <const char* prefix>
class ClosePairRejectionTrackTrack
{
 public:
  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr,
            int absChargeTrack1,
            int absChargeTrack2)
  {
    mIsActivated = confCpr.on.value;
    if (mIsActivated) {
      mCtr.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMax.value, confCpr.dphistarMax.value, absChargeTrack1, absChargeTrack2);
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
  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr)
  {
    mIsActivated = confCpr.on.value;
    if (mIsActivated) {
      mCtrPos.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMaxPosDau.value, confCpr.dphistarMaxPosDau.value, 1, 1);
      mCtrNeg.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMaxNegDau.value, confCpr.dphistarMaxNegDau.value, 1, 1);
    }
  }

  void setMagField(float magField)
  {
    mCtrPos.setMagField(magField);
    mCtrNeg.setMagField(magField);
  }
  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& v01, T2 const& v02, T3 const& tracks)
  {
    auto posDau1 = tracks.rawIteratorAt(v01.posDauId() - tracks.offset());
    auto negDau1 = tracks.rawIteratorAt(v01.negDauId() - tracks.offset());

    auto posDau2 = tracks.rawIteratorAt(v02.posDauId() - tracks.offset());
    auto negDau2 = tracks.rawIteratorAt(v02.negDauId() - tracks.offset());

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

template <const char* prefixSameCharge, const char* prefixOppositeCharge>
class ClosePairRejectionTrackV0 // can also be used for any particle type that has pos/neg daughters, like resonances
{
 public:
  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr,
            int absChargeTrack)
  {
    mIsActivatedSameCharge = confCpr.onSameCharge.value;
    if (mIsActivatedSameCharge) {
      mCtrSameCharge.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMaxSameCharge.value, confCpr.dphistarMaxSameCharge.value, absChargeTrack, 1);
    }

    mIsActivatedOppositeCharge = confCpr.onOppositeCharge.value;
    if (mIsActivatedOppositeCharge) {
      mCtrOppositeCharge.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMaxOppositeCharge.value, confCpr.dphistarMaxOppositeCharge.value, absChargeTrack, 1);
    }

    mIsActivated = mIsActivatedSameCharge || mIsActivatedOppositeCharge;
  }

  void setMagField(float magField)
  {
    if (mIsActivatedSameCharge) {
      mCtrSameCharge.setMagField(magField);
    }
    if (mIsActivatedOppositeCharge) {
      mCtrOppositeCharge.setMagField(magField);
    }
  }

  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& track, T2 const& v0, T3 const& trackTable)
  {
    auto posDau = trackTable.rawIteratorAt(v0.posDauId() - trackTable.offset());
    auto negDau = trackTable.rawIteratorAt(v0.negDauId() - trackTable.offset());
    if (track.sign() > 0) {
      if (mIsActivatedSameCharge) {
        mCtrSameCharge.compute(track, posDau);
      }
      if (mIsActivatedOppositeCharge) {
        mCtrOppositeCharge.compute(track, negDau);
      }
    } else {
      if (mIsActivatedSameCharge) {
        mCtrSameCharge.compute(track, negDau);
      }
      if (mIsActivatedOppositeCharge) {
        mCtrOppositeCharge.compute(track, posDau);
      }
    }
  }

  bool isClosePair() const
  {
    bool cprSameCharge = mIsActivatedSameCharge && mCtrSameCharge.isClosePair();
    bool cprOppositeCharrge = mIsActivatedOppositeCharge && mCtrOppositeCharge.isClosePair();
    return cprSameCharge || cprOppositeCharrge;
  }

  void fill()
  {
    if (mIsActivatedSameCharge) {
      mCtrSameCharge.fill();
    }
    if (mIsActivatedOppositeCharge) {
      mCtrOppositeCharge.fill();
    }
  }
  bool isActivated() const { return mIsActivated; }

 private:
  CloseTrackRejection<prefixSameCharge> mCtrSameCharge;
  CloseTrackRejection<prefixOppositeCharge> mCtrOppositeCharge;
  bool mIsActivated = true;
  bool mIsActivatedSameCharge = true;
  bool mIsActivatedOppositeCharge = false;
};

template <const char* prefixBachelor, const char* prefixSameCharge, const char* prefixOppositeCharge>
class ClosePairRejectionTrackCascade
{
 public:
  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr,
            int absChargeTrack)
  {
    mIsActivatedBachelor = confCpr.onBachelor.value;
    if (mIsActivatedBachelor) {
      mCtrBachelor.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMaxBachelor.value, confCpr.dphistarMaxBachelor.value, absChargeTrack, 1);
    }

    mIsActivatedSameChargeV0Daughter = confCpr.onSameCharge.value;
    if (mIsActivatedSameChargeV0Daughter) {
      mCtrSameChargeV0Daughter.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMaxSameCharge.value, confCpr.dphistarMaxSameCharge.value, absChargeTrack, 1);
    }

    mIsActivatedOppositeChargeV0Daughter = confCpr.onOppositeCharge.value;
    if (mIsActivatedOppositeChargeV0Daughter) {
      mCtrOppositeChargeV0Daughter.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMaxOppositeCharge.value, confCpr.dphistarMaxOppositeCharge.value, absChargeTrack, 1);
    }

    mIsActivated = mIsActivatedBachelor || mIsActivatedSameChargeV0Daughter || mIsActivatedOppositeChargeV0Daughter;
  }

  void setMagField(float magField)
  {
    if (mIsActivatedBachelor) {
      mCtrBachelor.setMagField(magField);
    }
    if (mIsActivatedSameChargeV0Daughter) {
      mCtrSameChargeV0Daughter.setMagField(magField);
    }
    if (mIsActivatedOppositeChargeV0Daughter) {
      mCtrOppositeChargeV0Daughter.setMagField(magField);
    }
  }
  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& track, T2 const& cascade, T3 const& trackTable)
  {
    auto bachelor = trackTable.rawIteratorAt(cascade.bachelorId() - trackTable.offset());
    auto posDau = trackTable.rawIteratorAt(cascade.posDauId() - trackTable.offset());
    auto negDau = trackTable.rawIteratorAt(cascade.negDauId() - trackTable.offset());

    if (mIsActivatedBachelor) {
      mCtrBachelor.compute(track, bachelor);
    }

    if (track.sign() > 0) {
      if (mIsActivatedSameChargeV0Daughter) {
        mCtrSameChargeV0Daughter.compute(track, posDau);
      }
      if (mIsActivatedOppositeChargeV0Daughter) {
        mCtrOppositeChargeV0Daughter.compute(track, negDau);
      }
    } else {
      if (mIsActivatedSameChargeV0Daughter) {
        mCtrSameChargeV0Daughter.compute(track, negDau);
      }
      if (mIsActivatedOppositeChargeV0Daughter) {
        mCtrOppositeChargeV0Daughter.compute(track, posDau);
      }
    }
  }

  bool isClosePair() const
  {
    bool cprBachelor = mIsActivatedBachelor && mCtrBachelor.isClosePair();
    bool cprSameCharge = mIsActivatedSameChargeV0Daughter && mCtrSameChargeV0Daughter.isClosePair();
    bool cprOppositeCharrge = mIsActivatedOppositeChargeV0Daughter && mCtrOppositeChargeV0Daughter.isClosePair();
    return cprBachelor || cprSameCharge || cprOppositeCharrge;
  }

  void fill()
  {
    if (mIsActivatedBachelor) {
      mCtrBachelor.fill();
    }
    if (mIsActivatedSameChargeV0Daughter) {
      mCtrSameChargeV0Daughter.fill();
    }
    if (mIsActivatedOppositeChargeV0Daughter) {
      mCtrOppositeChargeV0Daughter.fill();
    }
  }

  bool isActivated() const { return mIsActivated; }

 private:
  CloseTrackRejection<prefixBachelor> mCtrBachelor;
  CloseTrackRejection<prefixSameCharge> mCtrSameChargeV0Daughter;
  CloseTrackRejection<prefixOppositeCharge> mCtrOppositeChargeV0Daughter;
  bool mIsActivated = true;
  bool mIsActivatedBachelor = false;
  bool mIsActivatedSameChargeV0Daughter = false;
  bool mIsActivatedOppositeChargeV0Daughter = false;
};

template <const char* prefix>
class ClosePairRejectionTrackKink
{
 public:
  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr,
            int absChargeTrack)
  {
    mIsActivated = confCpr.on.value;
    if (mIsActivated) {
      mCtr.init(registry, specs, confCpr.plotAverage.value, confCpr.plotAllRadii.value, confCpr.detaMax.value, confCpr.dphistarMax.value, absChargeTrack, 1);
    }
  }

  void setMagField(float magField)
  {
    mCtr.setMagField(magField);
  }

  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& track, T2 const& kink, T3 const& trackTable)
  {
    auto daughter = trackTable.rawIteratorAt(kink.chaDauId() - trackTable.offset());
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
