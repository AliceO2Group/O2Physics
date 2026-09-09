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

/// \file closeTripletRejection.h
/// \brief Definition of CloseTripletRejection class
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_CLOSETRIPLETREJECTION_H_
#define PWGCF_FEMTO_CORE_CLOSETRIPLETREJECTION_H_

#include "PWGCF/Femto/Core/closePairRejection.h"

#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <map>
#include <vector>

namespace o2::analysis::femto::closetripletrejection
{
constexpr const char PrefixCtrTrackTrackTrack[] = "CtrTrackTrackTrack";
constexpr const char PrefixCtrTrackTrackV0[] = "CtrTrackTrackV0";
constexpr const char PrefixCtrTrackTrackCascade[] = "CtrTrackTrackCascade";

using ConfCtrTrackTrackTrack = closepairrejection::ConfCpr<PrefixCtrTrackTrackTrack>;
using ConfCtrTrackTrackV0 = closepairrejection::ConfCpr<PrefixCtrTrackTrackV0>;
using ConfCtrTrackTrackCascade = closepairrejection::ConfCpr<PrefixCtrTrackTrackCascade>;

// directory names
constexpr char PrefixTrack1Track2Se[] = "CPR_Track1Track2/SE/";
constexpr char PrefixTrack2Track3Se[] = "CPR_Track2Track3/SE/";
constexpr char PrefixTrack1Track3Se[] = "CPR_Track1Track3/SE/";
constexpr char PrefixTrack1Track2Me[] = "CPR_Track1Track2/ME/";
constexpr char PrefixTrack2Track3Me[] = "CPR_Track2Track3/ME/";
constexpr char PrefixTrack1Track3Me[] = "CPR_Track1Track3/ME/";

constexpr char PrefixTrack1V0Se[] = "CPR_Track1V0/SE/";
constexpr char PrefixTrack2V0Se[] = "CPR_Track2V0/SE/";
constexpr char PrefixTrack1V0Me[] = "CPR_Track1V0/ME/";
constexpr char PrefixTrack2V0Me[] = "CPR_Track2V0/ME/";

constexpr char PrefixTrack1CascadeSe[] = "CPR_Track1Cascade/SE/";
constexpr char PrefixTrack2CascadeSe[] = "CPR_Track2Cascade/SE/";
constexpr char PrefixTrack1CascadeMe[] = "CPR_Track1Cascade/ME/";
constexpr char PrefixTrack2CascadeMe[] = "CPR_Track2Cascade/ME/";

constexpr char PrefixTrack1V0DaughterSe[] = "CPR_Track1V0Dau/SE/";
constexpr char PrefixTrack2V0DaughterSe[] = "CPR_Track2V0Dau/SE/";
constexpr char PrefixTrack1V0DaughterMe[] = "CPR_Track1V0Dau/ME/";
constexpr char PrefixTrack2V0DaughterMe[] = "CPR_Track2V0Dau/ME/";

constexpr char PrefixTrack1CascadeBachelorSe[] = "CPR_Track1CascadeBachelor/SE/";
constexpr char PrefixTrack2CascadeBachelorSe[] = "CPR_Track2CascadeBachelor/SE/";
constexpr char PrefixTrack1CascadeBachelorMe[] = "CPR_Track1CascadeBachelor/ME/";
constexpr char PrefixTrack2CascadeBachelorMe[] = "CPR_Track2CascadeBachelor/ME/";

template <auto& prefixTrack1Track2,
          auto& prefixTrack2Track3,
          auto& prefixTrack1Track3>
class CloseTripletRejectionTrackTrackTrack
{
 public:
  CloseTripletRejectionTrackTrackTrack() = default;
  ~CloseTripletRejectionTrackTrackTrack() = default;

  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr,
            int absChargeTrack1,
            int absChargeTrack2,
            int absChargeTrack3)
  {
    mCtrTrack12.init(registry, specs, confCpr, absChargeTrack1, absChargeTrack2);
    mCtrTrack23.init(registry, specs, confCpr, absChargeTrack2, absChargeTrack3);
    mCtrTrack13.init(registry, specs, confCpr, absChargeTrack1, absChargeTrack3);
  }

  void setMagField(float magField)
  {
    mCtrTrack12.setMagField(magField);
    mCtrTrack23.setMagField(magField);
    mCtrTrack13.setMagField(magField);
  }

  // checks all three constituent pairs of the triplet; fills the deta-dphi/kinematic
  // histograms of each pair internally and returns whether the triplet is rejected.
  // tripletHistManager must expose getKinematic() (same interface CloseTrackRejection expects).
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  [[nodiscard]] bool isCloseTriplet(T1 const& track1, T2 const& track2, T3 const& track3, T4 const& trackTable, T5 const& tripletHistManager)
  {
    bool isClose12 = mCtrTrack12.isClosePair(track1, track2, trackTable, tripletHistManager);
    bool isClose23 = mCtrTrack23.isClosePair(track2, track3, trackTable, tripletHistManager);
    bool isClose13 = mCtrTrack13.isClosePair(track1, track3, trackTable, tripletHistManager);
    return isClose12 || isClose23 || isClose13;
  }

 private:
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack1Track2> mCtrTrack12;
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack2Track3> mCtrTrack23;
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack1Track3> mCtrTrack13;
};

template <auto& prefixTrack1Track2,
          auto& prefixTrack1V0,
          auto& prefixTrack2V0>
class CloseTripletRejectionTrackTrackV0
{
 public:
  CloseTripletRejectionTrackTrackV0() = default;
  ~CloseTripletRejectionTrackTrackV0() = default;

  template <typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            T const& confCpr,
            int absChargeTrack1,
            int absChargeTrack2)
  {
    mCtrTrack12.init(registry, specs, confCpr, absChargeTrack1, absChargeTrack2);
    mCtrTrack1V0.init(registry, specs, confCpr, absChargeTrack1);
    mCtrTrack2V0.init(registry, specs, confCpr, absChargeTrack2);
  }

  void setMagField(float magField)
  {
    mCtrTrack12.setMagField(magField);
    mCtrTrack1V0.setMagField(magField);
    mCtrTrack2V0.setMagField(magField);
  }

  // checks track1-track2, track1-v0 and track2-v0; fills the deta-dphi/kinematic
  // histograms of each pair internally and returns whether the triplet is rejected.
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  [[nodiscard]] bool isCloseTriplet(T1 const& track1, T2 const& track2, T3 const& v0, T4 const& trackTable, T5 const& tripletHistManager)
  {
    bool isClose12 = mCtrTrack12.isClosePair(track1, track2, trackTable, tripletHistManager);
    bool isClose1V0 = mCtrTrack1V0.isClosePair(track1, v0, trackTable, tripletHistManager);
    bool isClose2V0 = mCtrTrack2V0.isClosePair(track2, v0, trackTable, tripletHistManager);
    return isClose12 || isClose1V0 || isClose2V0;
  }

 private:
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack1Track2> mCtrTrack12;
  closepairrejection::ClosePairRejectionTrackV0<prefixTrack1V0> mCtrTrack1V0;
  closepairrejection::ClosePairRejectionTrackV0<prefixTrack2V0> mCtrTrack2V0;
};

template <auto& prefixTrack1Track2,
          auto& prefixTrack1Bachelor,
          auto& prefixTrack1V0Daughter,
          auto& prefixTrack2Bachelor,
          auto& prefixTrack2V0Daughter>
class CloseTripletRejectionTrackTrackCascade
{
 public:
  CloseTripletRejectionTrackTrackCascade() = default;
  ~CloseTripletRejectionTrackTrackCascade() = default;

  template <typename T1, typename T2, typename T3>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> const& specs,
            std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> const& specsBachelor,
            std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> const& specsV0Daughter,
            T1 const& confCpr,
            T2 const& confCprBachelor,
            T3 const& confCprV0Daughter,
            int absChargeTrack1,
            int absChargeTrack2)
  {
    mCtrTrack12.init(registry, specs, confCpr, absChargeTrack1, absChargeTrack2);
    mCtrTrack1Cascade.init(registry, specsBachelor, specsV0Daughter, confCprBachelor, confCprV0Daughter, absChargeTrack1);
    mCtrTrack2Cascade.init(registry, specsBachelor, specsV0Daughter, confCprBachelor, confCprV0Daughter, absChargeTrack2);
  }

  void setMagField(float magField)
  {
    mCtrTrack12.setMagField(magField);
    mCtrTrack1Cascade.setMagField(magField);
    mCtrTrack2Cascade.setMagField(magField);
  }

  // checks track1-track2, track1-cascade and track2-cascade; fills the deta-dphi/kinematic
  // histograms of each pair internally and returns whether the triplet is rejected.
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  [[nodiscard]] bool isCloseTriplet(T1 const& track1, T2 const& track2, T3 const& cascade, T4 const& trackTable, T5 const& tripletHistManager)
  {
    bool isClose12 = mCtrTrack12.isClosePair(track1, track2, trackTable, tripletHistManager);
    bool isClose1Cascade = mCtrTrack1Cascade.isClosePair(track1, cascade, trackTable, tripletHistManager);
    bool isClose2Cascade = mCtrTrack2Cascade.isClosePair(track2, cascade, trackTable, tripletHistManager);
    return isClose12 || isClose1Cascade || isClose2Cascade;
  }

 private:
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack1Track2> mCtrTrack12;
  closepairrejection::ClosePairRejectionTrackCascade<prefixTrack1Bachelor, prefixTrack1V0Daughter> mCtrTrack1Cascade;
  closepairrejection::ClosePairRejectionTrackCascade<prefixTrack2Bachelor, prefixTrack2V0Daughter> mCtrTrack2Cascade;
};

} // namespace o2::analysis::femto::closetripletrejection
#endif // PWGCF_FEMTO_CORE_CLOSETRIPLETREJECTION_H_
