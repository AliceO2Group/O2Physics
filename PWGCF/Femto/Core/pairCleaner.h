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

/// \file pairCleaner.h
/// \brief pair cleaner class
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PAIRCLEANER_H_
#define PWGCF_FEMTO_CORE_PAIRCLEANER_H_

namespace o2::analysis::femto
{
namespace paircleaner
{

class BasePairCleaner
{
 public:
  BasePairCleaner() = default;
  virtual ~BasePairCleaner() = default;

 protected:
  template <typename T1, typename T2>
  bool isCleanTrackPair(T1 const& track1, T2 const& track2) const
  {
    return track1.globalIndex() != track2.globalIndex();
  };
};

class TrackTrackPairCleaner : public BasePairCleaner
{
 public:
  TrackTrackPairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(T1 const& track1, T2 const& track2, T3 const& /*trackTable*/) const
  {
    return this->isCleanTrackPair(track1, track2);
  }
};

class V0V0PairCleaner : public BasePairCleaner
{
 public:
  V0V0PairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& v01, const T2& v02, const T3& /*tracks*/) const
  {
    auto posDaughter1 = v01.template posDau_as<T3>();
    auto negDaughter1 = v01.template negDau_as<T3>();
    auto posDaughter2 = v02.template posDau_as<T3>();
    auto negDaughter2 = v02.template negDau_as<T3>();
    return this->isCleanTrackPair(posDaughter1, posDaughter2) && this->isCleanTrackPair(negDaughter1, negDaughter2);
  }
};

class TrackV0PairCleaner : public BasePairCleaner // also works for particles decaying into a positive and negative daughter, like resonances
{
 public:
  TrackV0PairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& track, const T2& v0, const T3& /*trackTable*/) const
  {
    auto posDaughter = v0.template posDau_as<T3>();
    auto negDaughter = v0.template negDau_as<T3>();
    return (this->isCleanTrackPair(posDaughter, track) && this->isCleanTrackPair(negDaughter, track));
  }
};

class TrackKinkPairCleaner : public BasePairCleaner
{
 public:
  TrackKinkPairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& track, const T2& kink, const T3& /*trackTable*/) const
  {
    auto chaDaughter = kink.template chaDau_as<T3>();
    return this->isCleanTrackPair(chaDaughter, track);
  }
};

class TrackCascadePairCleaner : public BasePairCleaner
{
 public:
  TrackCascadePairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& track, const T2& cascade, const T3& /*trackTable*/) const
  {
    auto bachelor = cascade.template bachelor_as<T3>();
    auto posDaughter = cascade.template posDau_as<T3>();
    auto negDaughter = cascade.template negDau_as<T3>();
    return (this->isCleanTrackPair(bachelor, track) && this->isCleanTrackPair(posDaughter, track) && this->isCleanTrackPair(negDaughter, track));
  }
};
} // namespace paircleaner
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_PAIRCLEANER_H_
