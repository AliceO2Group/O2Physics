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

#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include <cmath>
#include <map>
#include <vector>

namespace o2::analysis::femto
{
namespace closetripletrejection
{

constexpr const char PrefixCtrTrackTrackTrack[] = "CtrTrackTrackTrack";
using ConfCtrTrackTrackTrack = closepairrejection::ConfCpr<PrefixCtrTrackTrackTrack>;

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

template <const char* prefixTrack1Track2,
          const char* prefixTrack2Track3,
          const char* prefixTrack1Track3>
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
  template <typename T1, typename T2, typename T3, typename T4>
  void setTriplet(T1 const& track1, T2 const& track2, T3 const& track3, T4 const& trackTable)
  {
    mCtrTrack12.setPair(track1, track2, trackTable);
    mCtrTrack23.setPair(track2, track3, trackTable);
    mCtrTrack13.setPair(track1, track3, trackTable);
  }
  bool isCloseTriplet() const
  {
    return mCtrTrack12.isClosePair() || mCtrTrack23.isClosePair() || mCtrTrack13.isClosePair();
  }

  void fill(float q3)
  {
    mCtrTrack12.fill(q3);
    mCtrTrack23.fill(q3);
    mCtrTrack13.fill(q3);
  }

 private:
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack1Track2> mCtrTrack12;
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack2Track3> mCtrTrack23;
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack1Track3> mCtrTrack13;
};

template <const char* prefixTrack1Track2,
          const char* prefixTrack1V0,
          const char* prefixTrack2V0>
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
  template <typename T1, typename T2, typename T3, typename T4>
  void setTriplet(T1 const& track1, T2 const& track2, T3 const& v0, T4 const& trackTable)
  {
    mCtrTrack12.setPair(track1, track2, trackTable);
    mCtrTrack1V0.setPair(track1, v0, trackTable);
    mCtrTrack2V0.setPair(track2, v0, trackTable);
  }
  bool isCloseTriplet() const
  {
    return mCtrTrack12.isClosePair() || mCtrTrack1V0.isClosePair() || mCtrTrack2V0.isClosePair();
  }

  void fill(float q3)
  {
    mCtrTrack12.fill(q3);
    mCtrTrack1V0.fill(q3);
    mCtrTrack2V0.fill(q3);
  }

 private:
  closepairrejection::ClosePairRejectionTrackTrack<prefixTrack1Track2> mCtrTrack12;
  closepairrejection::ClosePairRejectionTrackV0<prefixTrack1V0> mCtrTrack1V0;
  closepairrejection::ClosePairRejectionTrackV0<prefixTrack2V0> mCtrTrack2V0;
};

}; // namespace closetripletrejection
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_CLOSETRIPLETREJECTION_H_
