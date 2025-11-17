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

#include "RecoDecay.h"

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

// template configurable group for Cpr
template <const char* Prefix>
struct ConfCpr : o2::framework::ConfigurableGroup {
  std::string prefix = std::string(Prefix);
  o2::framework::Configurable<bool> cutAverage{"cutAverage", true, "Apply CPR if the average deta-dphistar is below the configured values"};
  o2::framework::Configurable<bool> cutAnyRadius{"cutAnyRadius", false, "Apply CPR if the deta-dphistar is below the configured values at any radius"};
  o2::framework::Configurable<bool> plotAllRadii{"plotAllRadii", true, "Plot deta-dphi distribution at all radii"};
  o2::framework::Configurable<bool> plotAverage{"plotAverage", true, "Plot average deta dphi distribution"};
  o2::framework::Configurable<float> detaMax{"detaMax", 0.01f, "Maximium deta"};
  o2::framework::Configurable<float> dphistarMax{"dphistarMax", 0.01f, "Maximum dphistar"};
  o2::framework::Configurable<float> detaCenter{"detaCenter", 0.f, "Center of deta cut"};
  o2::framework::Configurable<float> dphistarCenter{"dphistarCenter", 0.f, "Center of dphistar cut"};
  o2::framework::Configurable<float> kstarMin{"kstarMin", -1.f, "Minimum kstar of pair for plotting (Set to negative value to turn off the cut)"};
  o2::framework::Configurable<float> kstarMax{"kstarMax", -1.f, "Maximum kstar of pair for plotting (Set to negative value to turn off the cut)"};
  o2::framework::ConfigurableAxis binningDeta{"binningDeta", {{250, -0.5, 0.5}}, "deta"};
  o2::framework::ConfigurableAxis binningDphistar{"binningDphistar", {{250, -0.5, 0.5}}, "dphi"};
};

constexpr const char PrefixCprTrackTrack[] = "CprTrackTrack";
constexpr const char PrefixCprTrackV0Daughter[] = "CprTrackV0Daughter";
constexpr const char PrefixCprTrackResonanceDaughter[] = "CprTrackResonanceDaughter";
constexpr const char PrefixCprTrackKinkDaughter[] = "CprTrackKinkDaughter";
constexpr const char PrefixCprV0DaughterV0DaughterPos[] = "CprV0DaughterV0DaughterPos";
constexpr const char PrefixCprV0DaughterV0DaughterNeg[] = "CprV0DaughterV0DaughterNeg";
constexpr const char PrefixCprTrackCascadeBachelor[] = "CprTrackCascadeBachelor";

using ConfCprTrackTrack = ConfCpr<PrefixCprTrackTrack>;
using ConfCprTrackV0Daughter = ConfCpr<PrefixCprTrackV0Daughter>;
using ConfCprTrackResonanceDaughter = ConfCpr<PrefixCprTrackResonanceDaughter>;
using ConfCprTrackKinkDaughter = ConfCpr<PrefixCprTrackKinkDaughter>;
using ConfCprV0DaugherV0DaughterPos = ConfCpr<PrefixCprV0DaughterV0DaughterPos>;
using ConfCprV0DaugherV0DaughterNeg = ConfCpr<PrefixCprV0DaughterV0DaughterNeg>;
using ConfCprTrackCascadeBachelor = ConfCpr<PrefixCprTrackCascadeBachelor>;

// tpc radii for computing phistar
constexpr int Nradii = 9;
constexpr std::array<float, Nradii> TpcRadii = {85., 105., 125., 145., 165., 185., 205., 225., 245.}; // in cm

// directory names
constexpr char PrefixTrackTrackSe[] = "CPR_TrackTrack/SE/";
constexpr char PrefixTrackTrackMe[] = "CPR_TrackTrack/ME/";
constexpr char PrefixTrackV0DaughterSe[] = "CPR_TrackV0Dau/SE/";
constexpr char PrefixTrackV0DaughterMe[] = "CPR_TrackV0Dau/ME/";
constexpr char PrefixV0V0PosSe[] = "CPR_V0V0_PosDau/SE/";
constexpr char PrefixV0V0NegSe[] = "CPR_V0V0_NegDau/SE/";
constexpr char PrefixV0V0PosMe[] = "CPR_V0V0_PosDau/ME/";
constexpr char PrefixV0V0NegMe[] = "CPR_V0V0_NegDau/ME/";
constexpr char PrefixTrackTwoTrackResonanceSe[] = "CPR_TrackResonanceDau/SE/";
constexpr char PrefixTrackTwoTrackResonanceMe[] = "CPR_TrackResonanceDau/ME/";
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

  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr,
            int chargeAbsTrack1,
            int chargeAbsTrack2)
  {
    mDetaMax = confCpr.detaMax.value;
    mDphistarMax = confCpr.dphistarMax.value;

    // check the limits
    if (mDetaMax <= 0 || mDphistarMax <= 0) {
      LOG(fatal) << "Limits for Close Pair Rejection are invalid (0 or negative). Breaking...";
    }

    mDetaCenter = confCpr.detaCenter.value;
    mDphistarCenter = confCpr.dphistarCenter.value;

    mChargeAbsTrack1 = std::abs(chargeAbsTrack1);
    mChargeAbsTrack2 = std::abs(chargeAbsTrack2);

    mCutAverage = confCpr.cutAverage.value;
    mCutAnyRadius = confCpr.cutAnyRadius.value;

    mKstarMin = confCpr.kstarMin.value;
    mKstarMax = confCpr.kstarMax.value;

    mPlotAverage = confCpr.plotAverage.value;
    mPlotAllRadii = confCpr.plotAllRadii.value;

    // check if we need to apply any cut a plot is requested
    mIsActivated = mCutAverage || mCutAnyRadius || mPlotAverage || mPlotAllRadii;

    mHistogramRegistry = registry;

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
    if (!mIsActivated) {
      return;
    }
    // reset values
    mAverageDphistar = 0.f;
    int count = 0;
    mDeta = 0.f;
    mDphistar.fill(0.f);
    mDphistarMask.fill(false);

    mDeta = track1.eta() - track2.eta();

    for (size_t i = 0; i < TpcRadii.size(); i++) {
      auto phistar1 = utils::dphistar(mMagField, TpcRadii[i], mChargeAbsTrack1 * track1.signedPt(), track1.phi());
      auto phistar2 = utils::dphistar(mMagField, TpcRadii[i], mChargeAbsTrack2 * track2.signedPt(), track2.phi());
      if (phistar1 && phistar2) {
        mDphistar.at(i) = RecoDecay::constrainAngle(phistar1.value() - phistar2.value(), -o2::constants::math::PI); // constrain angular difference between -pi and pi
        mDphistarMask.at(i) = true;
        count++;
      }
    }
    // for small momemeta the calculation of phistar might fail, if the particle did not reach a certain radius
    mAverageDphistar = std::accumulate(mDphistar.begin(), mDphistar.end(), 0.f) / count; // only average values if phistar could be computed
  }

  void fill(float kstar)
  {
    if (!mIsActivated) {
      return;
    }

    if (mKstarMin > 0.f && kstar < mKstarMin) {
      return;
    }

    if (mKstarMax > 0.f && kstar > mKstarMax) {
      return;
    }

    // fill average hist
    if (mPlotAverage) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kAverage, HistTable)), mDeta, mAverageDphistar);
    }

    // fill radii hists
    if (mPlotAllRadii) {
      if (mDphistarMask.at(0)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius0, HistTable)), mDeta, mDphistar.at(0));
      }
      if (mDphistarMask.at(1)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius1, HistTable)), mDeta, mDphistar.at(1));
      }
      if (mDphistarMask.at(2)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius2, HistTable)), mDeta, mDphistar.at(2));
      }
      if (mDphistarMask.at(3)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius3, HistTable)), mDeta, mDphistar.at(3));
      }
      if (mDphistarMask.at(4)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius4, HistTable)), mDeta, mDphistar.at(4));
      }
      if (mDphistarMask.at(5)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius5, HistTable)), mDeta, mDphistar.at(5));
      }
      if (mDphistarMask.at(6)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius6, HistTable)), mDeta, mDphistar.at(6));
      }
      if (mDphistarMask.at(7)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius7, HistTable)), mDeta, mDphistar.at(7));
      }
      if (mDphistarMask.at(8)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kRadius8, HistTable)), mDeta, mDphistar.at(8));
      }
    }
  }

  bool isClosePair() const
  {
    if (!mIsActivated) {
      return false;
    }
    bool isCloseAverage = false;
    bool isCloseAnyRadius = false;

    if (mCutAverage) {
      isCloseAverage = std::hypot((mAverageDphistar - mDphistarCenter) / mDphistarMax, (mDeta - mDetaCenter) / mDetaMax) < 1.f;
    }

    if (mCutAnyRadius) {
      for (size_t i = 0; i < TpcRadii.size(); i++) {
        if (isCloseAnyRadius) {
          break;
        }
        if (mDphistarMask.at(i)) {
          isCloseAnyRadius = std::hypot((mDphistar.at(i) - mDphistarCenter) / mDphistarMax, (mDeta - mDetaCenter) / mDetaMax) < 1.f;
        }
      }
    }
    return isCloseAverage || isCloseAnyRadius;
  }

  bool isActivated() const { return mIsActivated; }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mPlotAllRadii = false;
  bool mPlotAverage = false;

  float mKstarMin = -1.f;
  float mKstarMax = -1.f;

  bool mCutAverage = false;
  bool mCutAnyRadius = false;

  bool mIsActivated = false;

  int mChargeAbsTrack1 = 0;
  int mChargeAbsTrack2 = 0;
  float mMagField = 0.f;
  float mDetaMax = 0.f;
  float mDphistarMax = 0.f;
  float mDetaCenter = 0.f;
  float mDphistarCenter = 0.f;

  float mAverageDphistar = 0.f;
  float mDeta = 0.f;
  std::array<float, Nradii> mDphistar = {0.f};
  std::array<bool, Nradii> mDphistarMask = {false};
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
    mCtr.init(registry, specs, confCpr, absChargeTrack1, absChargeTrack2);
  }

  void setMagField(float magField) { mCtr.setMagField(magField); }
  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& track1, T2 const& track2, T3 const& /*tracks*/)
  {
    mCtr.compute(track1, track2);
  }
  bool isClosePair() const { return mCtr.isClosePair(); }
  void fill(float kstar) { mCtr.fill(kstar); }

 private:
  CloseTrackRejection<prefix> mCtr;
};

template <const char* prefixPosDaus, const char* prefixNegDaus>
class ClosePairRejectionV0V0
{
 public:
  template <typename T1, typename T2>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specsPos,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specsNeg,
            T1 const& confCprPos,
            T2 const& confCprNeg)
  {
    mCtrPos.init(registry, specsPos, confCprPos, 1, 1);
    mCtrNeg.init(registry, specsNeg, confCprNeg, 1, 1);
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
    auto posDau2 = tracks.rawIteratorAt(v02.posDauId() - tracks.offset());
    mCtrPos.compute(posDau1, posDau2);

    auto negDau1 = tracks.rawIteratorAt(v01.negDauId() - tracks.offset());
    auto negDau2 = tracks.rawIteratorAt(v02.negDauId() - tracks.offset());
    mCtrNeg.compute(negDau1, negDau2);
  }

  bool isClosePair() const { return mCtrPos.isClosePair() || mCtrNeg.isClosePair(); }

  void fill(float kstar)
  {
    mCtrPos.fill(kstar);
    mCtrNeg.fill(kstar);
  }

 private:
  CloseTrackRejection<prefixPosDaus> mCtrPos;
  CloseTrackRejection<prefixNegDaus> mCtrNeg;
};

template <const char* prefixTrackV0>
class ClosePairRejectionTrackV0 // can also be used for any particle type that has pos/neg daughters, like resonances
{
 public:
  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr,
            int absChargeTrack)
  {
    mCtr.init(registry, specs, confCpr, absChargeTrack, 1);
  }

  void setMagField(float magField) { mCtr.setMagField(magField); }

  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& track, T2 const& v0, T3 const& trackTable)
  {
    if (track.sign() > 0) {
      auto posDau = trackTable.rawIteratorAt(v0.posDauId() - trackTable.offset());
      mCtr.compute(track, posDau);
    } else {
      auto negDau = trackTable.rawIteratorAt(v0.negDauId() - trackTable.offset());
      mCtr.compute(track, negDau);
    }
  }

  bool isClosePair() const { return mCtr.isClosePair(); }

  void fill(float kstar) { mCtr.fill(kstar); }

 private:
  CloseTrackRejection<prefixTrackV0> mCtr;
};

template <const char* prefixBachelor, const char* prefixV0Daughter>
class ClosePairRejectionTrackCascade
{
 public:
  template <typename T1, typename T2>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specsBachelor,
            std::map<CprHist, std::vector<o2::framework::AxisSpec>> const& specsV0Daughter,
            T1 const& confCprBachelor,
            T2 const& confCprV0Daughter,
            int absChargeTrack)
  {
    mCtrBachelor.init(registry, specsBachelor, confCprBachelor, absChargeTrack, 1);
    mCtrV0Daughter.init(registry, specsV0Daughter, confCprV0Daughter, absChargeTrack, 1);
  }

  void setMagField(float magField)
  {
    mCtrBachelor.setMagField(magField);
    mCtrV0Daughter.setMagField(magField);
  }

  template <typename T1, typename T2, typename T3>
  void setPair(T1 const& track, T2 const& cascade, T3 const& trackTable)
  {
    auto bachelor = trackTable.rawIteratorAt(cascade.bachelorId() - trackTable.offset());
    mCtrBachelor.compute(track, bachelor);

    if (track.sign() > 0) {
      auto posDau = trackTable.rawIteratorAt(cascade.posDauId() - trackTable.offset());
      mCtrV0Daughter.compute(track, posDau);
    } else {
      auto negDau = trackTable.rawIteratorAt(cascade.negDauId() - trackTable.offset());
      mCtrV0Daughter.compute(track, negDau);
    }
  }

  bool
    isClosePair() const
  {
    return mCtrBachelor.isClosePair() || mCtrBachelor.isClosePair();
  }

  void fill(float kstar)
  {
    mCtrBachelor.fill(kstar);
    mCtrV0Daughter.fill(kstar);
  }

 private:
  CloseTrackRejection<prefixBachelor> mCtrBachelor;
  CloseTrackRejection<prefixV0Daughter> mCtrV0Daughter;
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
    mCtr.init(registry, specs, confCpr, absChargeTrack, 1);
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
  void fill(float kstar) { mCtr.fill(kstar); }

 private:
  CloseTrackRejection<prefix> mCtr;
};

}; // namespace closepairrejection
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_CLOSEPAIRREJECTION_H_
