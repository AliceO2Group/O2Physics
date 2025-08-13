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

#ifndef PWGCF_FEMTOUNITED_CORE_CLOSEPAIRREJECTION_H_
#define PWGCF_FEMTOUNITED_CORE_CLOSEPAIRREJECTION_H_

#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/Core/femtoUtils.h"
#include "PWGCF/FemtoUnited/Core/histManager.h"
#include "PWGCF/FemtoUnited/Core/modes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "Framework/HistogramRegistry.h"

#include <array>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femtounited
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
constexpr char PrefixTrackVzeroSe[] = "CPR_TrackVzero/SE/";
constexpr char PrefixTrackVzeroMe[] = "CPR_TrackVzero/ME/";

constexpr char PrefixTrackTrack[] = "TrackTrack/";
constexpr char PrefixTrackPosDau[] = "TrackPosDau/";
constexpr char PrefixTrackNegDau[] = "TrackNegDau/";
constexpr char PrefixTrackBachelor[] = "TrackBachelorDau/";

constexpr std::string_view AnalysisDir = "Analysis/";

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

class CloseTrackRejection
{
 public:
  CloseTrackRejection(float detaMax, float dphistarMax) : mDetaMax(detaMax), mDphistarMax(dphistarMax) {}

  void setMagField(float magField) { mMagField = magField; }

  void reset()
  {
    mSameCharge = true;
    mAverageDphistar = 0.f;
    mDeta = 0.f;
    mDphistar.fill(0.f);
  }

  template <typename T1, typename T2>
  void compute(const T1& track1, const T2& track2)
  {
    reset();
    if (track1.sign() != track2.sign()) {
      mSameCharge = false;
      return;
    }
    mDeta = track1.eta() - track2.eta();
    for (size_t i = 0; i < kTpcRadius.size(); ++i) {
      auto dphi1 = utils::dphistar(mMagField, kTpcRadius[i], track1.sign(), track1.pt(), track1.phi());
      auto dphi2 = utils::dphistar(mMagField, kTpcRadius[i], track2.sign(), track2.pt(), track2.phi());

      if (dphi1 && dphi2) {
        mDphistar[i] = dphi1.value() - dphi2.value();
      }
    }
    mAverageDphistar = std::accumulate(mDphistar.begin(), mDphistar.end(), 0.f) / mDphistar.size();
  }

  bool isClosePair() const { return mSameCharge && (std::hypot(mAverageDphistar / mDphistarMax, mDeta / mDetaMax) < 1.f); }
  float getAverageDphistar() const { return mAverageDphistar; }
  float getDeta() const { return mDeta; }
  const std::array<float, kNradii>& getDphistarArray() const { return mDphistar; }

 private:
  float mMagField = 0.f;
  float mAverageDphistar = 0.f;
  float mDeta = 0.f;
  float mDetaMax = 0.f;
  float mDphistarMax = 0.f;
  std::array<float, kNradii> mDphistar = {};
  bool mSameCharge = true;
};

template <const char* prefix, modes::Mode mode, modes::Pairs pair>
class ClosePairRejection
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>>& specs, float detaMax, float dphistarMax, bool isActivated)
  {
    mHistogramRegistry = registry;
    mIsActivated = isActivated;

    if constexpr (modes::isEqual(pair, modes::Pairs::kTrackTrack)) {
      mRejections[modes::TrackPairs::kTrackTrack] = std::make_unique<CloseTrackRejection>(detaMax, dphistarMax);
      std::string dir = std::string(prefix) + std::string(PrefixTrackTrack);
      mHistogramRegistry->add(dir + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius0, HistTable), GetHistDesc(kRadius0, HistTable), GetHistType(kRadius0, HistTable), {specs[kRadius0]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius1, HistTable), GetHistDesc(kRadius1, HistTable), GetHistType(kRadius1, HistTable), {specs[kRadius1]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius2, HistTable), GetHistDesc(kRadius2, HistTable), GetHistType(kRadius2, HistTable), {specs[kRadius2]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius3, HistTable), GetHistDesc(kRadius3, HistTable), GetHistType(kRadius3, HistTable), {specs[kRadius3]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius4, HistTable), GetHistDesc(kRadius4, HistTable), GetHistType(kRadius4, HistTable), {specs[kRadius4]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius5, HistTable), GetHistDesc(kRadius5, HistTable), GetHistType(kRadius5, HistTable), {specs[kRadius5]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius6, HistTable), GetHistDesc(kRadius6, HistTable), GetHistType(kRadius6, HistTable), {specs[kRadius6]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius7, HistTable), GetHistDesc(kRadius7, HistTable), GetHistType(kRadius7, HistTable), {specs[kRadius7]});
      mHistogramRegistry->add(dir + GetHistNamev2(kRadius8, HistTable), GetHistDesc(kRadius8, HistTable), GetHistType(kRadius8, HistTable), {specs[kRadius8]});
    }
    if constexpr (modes::isEqual(pair, modes::Pairs::kTrackV0) ||
                  modes::isEqual(pair, modes::Pairs::kTrackResonance)) {
      mRejections[modes::TrackPairs::kTrackPosDaughter] = std::make_unique<CloseTrackRejection>(detaMax, dphistarMax);
      mRejections[modes::TrackPairs::kTrackNegDaughter] = std::make_unique<CloseTrackRejection>(detaMax, dphistarMax);

      std::string dir1 = std::string(prefix) + std::string(PrefixTrackPosDau);
      mHistogramRegistry->add(dir1 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius0, HistTable), GetHistDesc(kRadius0, HistTable), GetHistType(kRadius0, HistTable), {specs[kRadius0]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius1, HistTable), GetHistDesc(kRadius1, HistTable), GetHistType(kRadius1, HistTable), {specs[kRadius1]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius2, HistTable), GetHistDesc(kRadius2, HistTable), GetHistType(kRadius2, HistTable), {specs[kRadius2]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius3, HistTable), GetHistDesc(kRadius3, HistTable), GetHistType(kRadius3, HistTable), {specs[kRadius3]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius4, HistTable), GetHistDesc(kRadius4, HistTable), GetHistType(kRadius4, HistTable), {specs[kRadius4]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius5, HistTable), GetHistDesc(kRadius5, HistTable), GetHistType(kRadius5, HistTable), {specs[kRadius5]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius6, HistTable), GetHistDesc(kRadius6, HistTable), GetHistType(kRadius6, HistTable), {specs[kRadius6]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius7, HistTable), GetHistDesc(kRadius7, HistTable), GetHistType(kRadius7, HistTable), {specs[kRadius7]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius8, HistTable), GetHistDesc(kRadius8, HistTable), GetHistType(kRadius8, HistTable), {specs[kRadius8]});
      std::string dir2 = std::string(prefix) + std::string(PrefixTrackNegDau);
      mHistogramRegistry->add(dir2 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius0, HistTable), GetHistDesc(kRadius0, HistTable), GetHistType(kRadius0, HistTable), {specs[kRadius0]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius1, HistTable), GetHistDesc(kRadius1, HistTable), GetHistType(kRadius1, HistTable), {specs[kRadius1]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius2, HistTable), GetHistDesc(kRadius2, HistTable), GetHistType(kRadius2, HistTable), {specs[kRadius2]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius3, HistTable), GetHistDesc(kRadius3, HistTable), GetHistType(kRadius3, HistTable), {specs[kRadius3]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius4, HistTable), GetHistDesc(kRadius4, HistTable), GetHistType(kRadius4, HistTable), {specs[kRadius4]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius5, HistTable), GetHistDesc(kRadius5, HistTable), GetHistType(kRadius5, HistTable), {specs[kRadius5]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius6, HistTable), GetHistDesc(kRadius6, HistTable), GetHistType(kRadius6, HistTable), {specs[kRadius6]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius7, HistTable), GetHistDesc(kRadius7, HistTable), GetHistType(kRadius7, HistTable), {specs[kRadius7]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius8, HistTable), GetHistDesc(kRadius8, HistTable), GetHistType(kRadius8, HistTable), {specs[kRadius8]});
    }
    if constexpr (modes::isEqual(pair, modes::Pairs::kTrackCascade)) {
      mRejections[modes::TrackPairs::kTrackPosDaughter] = std::make_unique<CloseTrackRejection>(detaMax, dphistarMax);
      mRejections[modes::TrackPairs::kTrackNegDaughter] = std::make_unique<CloseTrackRejection>(detaMax, dphistarMax);
      mRejections[modes::TrackPairs::kTrackBachelor] = std::make_unique<CloseTrackRejection>(detaMax, dphistarMax);

      std::string dir1 = std::string(prefix) + std::string(PrefixTrackPosDau);
      mHistogramRegistry->add(dir1 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius0, HistTable), GetHistDesc(kRadius0, HistTable), GetHistType(kRadius0, HistTable), {specs[kRadius0]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius1, HistTable), GetHistDesc(kRadius1, HistTable), GetHistType(kRadius1, HistTable), {specs[kRadius1]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius2, HistTable), GetHistDesc(kRadius2, HistTable), GetHistType(kRadius2, HistTable), {specs[kRadius2]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius3, HistTable), GetHistDesc(kRadius3, HistTable), GetHistType(kRadius3, HistTable), {specs[kRadius3]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius4, HistTable), GetHistDesc(kRadius4, HistTable), GetHistType(kRadius4, HistTable), {specs[kRadius4]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius5, HistTable), GetHistDesc(kRadius5, HistTable), GetHistType(kRadius5, HistTable), {specs[kRadius5]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius6, HistTable), GetHistDesc(kRadius6, HistTable), GetHistType(kRadius6, HistTable), {specs[kRadius6]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius7, HistTable), GetHistDesc(kRadius7, HistTable), GetHistType(kRadius7, HistTable), {specs[kRadius7]});
      mHistogramRegistry->add(dir1 + GetHistNamev2(kRadius8, HistTable), GetHistDesc(kRadius8, HistTable), GetHistType(kRadius8, HistTable), {specs[kRadius8]});
      std::string dir2 = std::string(prefix) + std::string(PrefixTrackNegDau);
      mHistogramRegistry->add(dir2 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius0, HistTable), GetHistDesc(kRadius0, HistTable), GetHistType(kRadius0, HistTable), {specs[kRadius0]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius1, HistTable), GetHistDesc(kRadius1, HistTable), GetHistType(kRadius1, HistTable), {specs[kRadius1]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius2, HistTable), GetHistDesc(kRadius2, HistTable), GetHistType(kRadius2, HistTable), {specs[kRadius2]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius3, HistTable), GetHistDesc(kRadius3, HistTable), GetHistType(kRadius3, HistTable), {specs[kRadius3]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius4, HistTable), GetHistDesc(kRadius4, HistTable), GetHistType(kRadius4, HistTable), {specs[kRadius4]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius5, HistTable), GetHistDesc(kRadius5, HistTable), GetHistType(kRadius5, HistTable), {specs[kRadius5]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius6, HistTable), GetHistDesc(kRadius6, HistTable), GetHistType(kRadius6, HistTable), {specs[kRadius6]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius7, HistTable), GetHistDesc(kRadius7, HistTable), GetHistType(kRadius7, HistTable), {specs[kRadius7]});
      mHistogramRegistry->add(dir2 + GetHistNamev2(kRadius8, HistTable), GetHistDesc(kRadius8, HistTable), GetHistType(kRadius8, HistTable), {specs[kRadius8]});
      std::string dir3 = std::string(prefix) + std::string(PrefixTrackBachelor);
      mHistogramRegistry->add(dir3 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {specs[kAverage]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius0, HistTable), GetHistDesc(kRadius0, HistTable), GetHistType(kRadius0, HistTable), {specs[kRadius0]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius1, HistTable), GetHistDesc(kRadius1, HistTable), GetHistType(kRadius1, HistTable), {specs[kRadius1]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius2, HistTable), GetHistDesc(kRadius2, HistTable), GetHistType(kRadius2, HistTable), {specs[kRadius2]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius3, HistTable), GetHistDesc(kRadius3, HistTable), GetHistType(kRadius3, HistTable), {specs[kRadius3]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius4, HistTable), GetHistDesc(kRadius4, HistTable), GetHistType(kRadius4, HistTable), {specs[kRadius4]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius5, HistTable), GetHistDesc(kRadius5, HistTable), GetHistType(kRadius5, HistTable), {specs[kRadius5]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius6, HistTable), GetHistDesc(kRadius6, HistTable), GetHistType(kRadius6, HistTable), {specs[kRadius6]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius7, HistTable), GetHistDesc(kRadius7, HistTable), GetHistType(kRadius7, HistTable), {specs[kRadius7]});
      mHistogramRegistry->add(dir3 + GetHistNamev2(kRadius8, HistTable), GetHistDesc(kRadius8, HistTable), GetHistType(kRadius8, HistTable), {specs[kRadius8]});
    }
  }

  void setMagField(float magField)
  {
    for (auto const& [_, rejection] : mRejections) {
      rejection->setMagField(magField);
    }
  }

  template <typename T1, typename T2>
  void setPair(const T1& track1, const T2& track2)
  {
    if constexpr (modes::isEqual(pair, modes::Pairs::kTrackTrack)) {
      mRejections.at(modes::TrackPairs::kTrackTrack)->compute(track1, track2);
    }
  }

  template <typename T1, typename T2, typename T3>
  void setPair(const T1& particle, const T2& track, const T3& tracks)
  {
    if constexpr (modes::isEqual(pair, modes::Pairs::kTrackV0) ||
                  modes::isEqual(pair, modes::Pairs::kTrackResonance)) {
      auto pos = particle.template posDau_as<T3>();
      auto neg = particle.template negDau_as<T3>();
      mRejections.at(modes::TrackPairs::kTrackPosDaughter)->compute(track, pos);
      mRejections.at(modes::TrackPairs::kTrackNegDaughter)->compute(track, neg);
    }

    if constexpr (modes::isEqual(pair, modes::Pairs::kTrackCascade)) {
      auto pos = particle.template posDau_as<T3>();
      auto neg = particle.template negDau_as<T3>();
      auto bach = particle.template bach_as<T3>();
      mRejections.at(modes::TrackPairs::kTrackPosDaughter)->compute(track, pos);
      mRejections.at(modes::TrackPairs::kTrackNegDaughter)->compute(track, neg);
      mRejections.at(modes::TrackPairs::kTrackBachelor)->compute(track, bach);
    }
  }

  bool isClosePair() const
  {
    for (const auto& [_, rejection] : mRejections) {
      if (rejection && rejection->isClosePair()) {
        return true;
      }
    }
    return false;
  }

  void fill()
  {
    if constexpr (modes::isEqual(pair, modes::Pairs::kTrackTrack)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kAverage, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getAverageDphistar());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius0, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(0));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius1, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(1));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius2, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(2));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius3, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(3));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius4, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(4));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius5, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(5));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius6, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(6));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius7, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(7));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackTrack) + HIST(GetHistName(kRadius8, HistTable)),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackTrack]->getDphistarArray().at(8));
    }
    if constexpr (modes::isEqual(pair, modes::Pairs::kTrackV0) ||
                  modes::isEqual(pair, modes::Pairs::kTrackResonance)) {
      // for pos daughter-track
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kAverage, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getAverageDphistar());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius0, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(0));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius1, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(1));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius2, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(2));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius3, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(3));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius4, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(4));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius5, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(5));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius6, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(6));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius7, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(7));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackPosDau) + HIST(GetHistName(kRadius8, HistTable)),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackPosDaughter]->getDphistarArray().at(8));

      // For TrackNegDaughter
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kAverage, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getAverageDphistar());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius0, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(0));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius1, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(1));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius2, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(2));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius3, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(3));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius4, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(4));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius5, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(5));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius6, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(6));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius7, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(7));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PrefixTrackNegDau) + HIST(GetHistName(kRadius8, HistTable)),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDeta(),
                               mRejections[modes::TrackPairs::kTrackNegDaughter]->getDphistarArray().at(8));
    }
  }

  bool isActivated() const { return mIsActivated; }

 private:
  std::unordered_map<modes::TrackPairs, std::unique_ptr<CloseTrackRejection>> mRejections;
  std::unordered_map<modes::TrackPairs, std::string> mPrefixMap;
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mIsActivated = true;
};

}; // namespace closepairrejection
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_CLOSEPAIRREJECTION_H_
