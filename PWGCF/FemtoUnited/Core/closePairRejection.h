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
#include <numeric>
#include <string>
#include <type_traits>
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
  o2::framework::Configurable<bool> on{"on", true, "Trun on CPR"};
  o2::framework::Configurable<float> detaMax{"detaMax", 0.01f, "Maximium deta"};
  o2::framework::Configurable<float> dphistarMax{"dphistarMax", 0.01f, "Maximum dphistar"};
  o2::framework::ConfigurableAxis binningDeta{"binningDeta", {{200, -0.2, 0.2}}, "deta"};
  o2::framework::ConfigurableAxis binningDphistar{"binningDphistar", {{200, -0.2, 0.2}}, "dphi"};
};

// tpc radii for computing phistar
constexpr int kNradii = 9;
constexpr std::array<float, kNradii> kTpcRadius = {85., 105., 125., 145., 165., 185., 205., 225., 245.};

// directory names
constexpr char PrefixTrackTrackSe[] = "TrackTrack/CPR/SE/";
constexpr char PrefixTrackTrackMe[] = "TrackTrack/CPR/ME/";
constexpr char PrefixTrackVzeroSe[] = "TrackVzero/CPR/SE/";
constexpr char PrefixTrackVzeroMe[] = "TrackVzero/CPR/ME/";
constexpr std::string_view AnalysisDir = "Analysis/";

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<CprHist>, kCprHistogramLast> HistTable = {
  {{kAverage, o2::framework::kTH2F, "hAverage", "Average: #Delta #eta vs #Delta #phi*; #Delta #eta; #Delta #phi*"},
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

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix>
class ClosePairRejection
{
 public:
  /// Destructor
  virtual ~ClosePairRejection() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<CprHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kAverage, HistTable), GetHistDesc(kAverage, HistTable), GetHistType(kAverage, HistTable), {Specs[kAverage]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius0, HistTable), GetHistDesc(kRadius0, HistTable), GetHistType(kRadius0, HistTable), {Specs[kRadius0]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius1, HistTable), GetHistDesc(kRadius1, HistTable), GetHistType(kRadius1, HistTable), {Specs[kRadius1]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius2, HistTable), GetHistDesc(kRadius2, HistTable), GetHistType(kRadius2, HistTable), {Specs[kRadius2]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius3, HistTable), GetHistDesc(kRadius3, HistTable), GetHistType(kRadius3, HistTable), {Specs[kRadius3]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius4, HistTable), GetHistDesc(kRadius4, HistTable), GetHistType(kRadius4, HistTable), {Specs[kRadius4]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius5, HistTable), GetHistDesc(kRadius5, HistTable), GetHistType(kRadius5, HistTable), {Specs[kRadius5]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius6, HistTable), GetHistDesc(kRadius6, HistTable), GetHistType(kRadius6, HistTable), {Specs[kRadius6]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius7, HistTable), GetHistDesc(kRadius7, HistTable), GetHistType(kRadius7, HistTable), {Specs[kRadius7]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kRadius8, HistTable), GetHistDesc(kRadius8, HistTable), GetHistType(kRadius8, HistTable), {Specs[kRadius8]});
    }

    // if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
    // std::string qaDir = std::string(prefix) + std::string(QaDir);
    // }
  }

  void activate(bool activate) { mIsActivated = activate; }
  bool isActivated() { return mIsActivated; }
  void setMagField(float magField) { mMagField = magField; }
  void setLimits(float detaMax, float dphistarMax)
  {
    mDetaMax = detaMax;
    mDphistarMax = dphistarMax;
  };

  template <typename T1, typename T2>
  void setPair(T1 const& particle1, T2 const& particle2)
  {
    if constexpr (std::is_same_v<T1, o2::aod::FUTracks::iterator> && std::is_same_v<T2, o2::aod::FUTracks::iterator>) {
      // set deta
      mDeta = particle1.eta() - particle2.eta();
      // set dphi at primary vertex
      mDphi = particle1.phi() - particle2.phi();
      // compute dphistar at different TPC radii
      for (int i = 0; i < kNradii; i++) {
        mDphistar.at(i) = (utils::dphistar(mMagField, kTpcRadius.at(i), particle1.sign(), particle1.pt(), particle1.phi()) -
                           utils::dphistar(mMagField, kTpcRadius.at(i), particle2.sign(), particle2.pt(), particle2.phi()));
        ;
      }
      // get average
      mAverageDphistar = std::accumulate(mDphistar.begin(), mDphistar.end(), 0.f) / mDphistar.size();
    }
  }

  bool isClosePair()
  {
    return std::pow(mAverageDphistar / mDphistarMax, 2) + std::pow(mDeta / mDetaMax, 2) < 1.;
  }

  template <modes::Mode mode>
  void fill()
  {
    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kAverage, HistTable)), mDeta, mAverageDphistar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius0, HistTable)), mDeta, mDphistar.at(0));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius1, HistTable)), mDeta, mDphistar.at(1));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius2, HistTable)), mDeta, mDphistar.at(2));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius3, HistTable)), mDeta, mDphistar.at(3));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius4, HistTable)), mDeta, mDphistar.at(4));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius5, HistTable)), mDeta, mDphistar.at(5));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius6, HistTable)), mDeta, mDphistar.at(6));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius7, HistTable)), mDeta, mDphistar.at(7));
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kRadius8, HistTable)), mDeta, mDphistar.at(8));
    }
    // to be implemented
    // if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
    // }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
  float mMagField = 0.f;

  float mDphi = 0;
  float mDphistarMax = 0.f;
  float mAverageDphistar = 0;
  std::array<float, kNradii> mDphistar{0};

  float mDetaMax = 0.f;
  float mDeta = 0.f;

  bool mIsActivated = false;
};
}; // namespace closepairrejection
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_CLOSEPAIRREJECTION_H_
