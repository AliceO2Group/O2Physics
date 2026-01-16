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
    mCtr1.init(registry, specs, confCpr, absChargeTrack1, absChargeTrack2);
    mCtr2.init(registry, specs, confCpr, absChargeTrack2, absChargeTrack3);
    mCtr3.init(registry, specs, confCpr, absChargeTrack1, absChargeTrack3);
  }

  void setMagField(float magField)
  {
    mCtr1.setMagField(magField);
    mCtr2.setMagField(magField);
    mCtr3.setMagField(magField);
  }
  template <typename T1, typename T2, typename T3, typename T4>
  void setTriplet(T1 const& track1, T2 const& track2, T3 const& track3, T4 const& /*tracks*/)
  {
    mCtr1.compute(track1, track2);
    mCtr2.compute(track2, track3);
    mCtr3.compute(track1, track3);
  }
  bool isCloseTriplet() const
  {
    return mCtr1.isClosePair() || mCtr2.isClosePair() || mCtr3.isClosePair();
  }

  void fill(float q3)
  {
    mCtr1.fill(q3);
    mCtr2.fill(q3);
    mCtr3.fill(q3);
  }

 private:
  closepairrejection::CloseTrackRejection<prefixTrack1Track2> mCtr1;
  closepairrejection::CloseTrackRejection<prefixTrack2Track3> mCtr2;
  closepairrejection::CloseTrackRejection<prefixTrack1Track3> mCtr3;
};

}; // namespace closetripletrejection
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_CLOSETRIPLETREJECTION_H_
