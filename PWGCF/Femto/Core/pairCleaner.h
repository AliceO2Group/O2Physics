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

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

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
  bool isCleanTrackPair(const T1& track1, const T2& track2) const
  {
    return track1.globalIndex() != track2.globalIndex();
  };
};

class TrackTrackPairCleaner : public BasePairCleaner
{
 public:
  TrackTrackPairCleaner() = default;
  template <typename T>
  bool isCleanPair(const T& track1, const T& track2) const
  {
    return this->isCleanTrackPair(track1, track2);
  }
};

class TrackV0PairCleaner : public BasePairCleaner
{
 public:
  TrackV0PairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& track, const T2& v0, const T3& /*trackTable */) const
  {
    auto posDaughter = v0.template posDau_as<T3>();
    auto negDaughter = v0.template negDau_as<T3>();
    return this->isCleanTrackPair(posDaughter, track) && this->isCleanTrackPair(negDaughter, track);
  }
};

// template <femtomodes::Mode mode>
// template <modes::Pairs pair>
// class PairCleaner
// {
//  public:
//   /// Destructor
//   virtual ~PairCleaner() = default;
//
//   template <typename T1, typename T2>
//   bool isCleanPair(const T1& particle1, const T2& particle2)
//   {
//     if constexpr (modes::isEqual(pair, modes::Pairs::kTrackTrack)) {
//       return particle1.globalIndex() != particle2.globalIndex();
//     }
//     return true;
//   };
//
//   template <typename T1, typename T2, typename T3>
//   bool isCleanPair(const T1& particle1, const T2& particle2, const T3& /*trackTable*/)
//   {
//     if constexpr (modes::isEqual(pair, modes::Pairs::kTrackV0) || modes::isEqual(pair, modes::Pairs::kTrackResonance)) {
//       auto posDaughter = particle2.template posDau_as<T3>();
//       auto negDaughter = particle2.template negDau_as<T3>();
//       return (particle1.globalIndex() != posDaughter.globalIndex() && particle1.globalIndex() != negDaughter.globalIndex());
//     }
//     if constexpr (modes::isEqual(pair, modes::Pairs::kTrackCascade)) {
//       auto posDaughter = particle2.template posDau_as<T3>();
//       auto negDaughter = particle2.template negDau_as<T3>();
//       auto bachelor = particle2.template bachelor_as<T3>();
//       return (particle1.globalIndex() != posDaughter.globalIndex() &&
//               particle1.globalIndex() != negDaughter.globalIndex() &&
//               particle1.globalIndex() != bachelor.globalIndex());
//     }
//     return true;
//   };
//
//  private:
// };
}; // namespace paircleaner
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTOUNITED_CORE_PAIRCLEANER_H_
