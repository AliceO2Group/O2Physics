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

#ifndef PWGCF_FEMTOUNITED_CORE_PAIRCLEANER_H_
#define PWGCF_FEMTOUNITED_CORE_PAIRCLEANER_H_

#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/Core/histManager.h"
#include "PWGCF/FemtoUnited/Core/modes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

namespace o2::analysis::femtounited
{
namespace paircleaner
{
/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
class PairCleaner
{
 public:
  /// Destructor
  virtual ~PairCleaner() = default;

  template <typename T1, typename T2>
  bool isCleanPair(T1 particle1, T2 particle2)
  {
    if constexpr (std::is_same_v<T1, o2::aod::FUTracks::iterator> && std::is_same_v<T2, o2::aod::FUTracks::iterator>) {
      return particle1.globalIndex() != particle2.globalIndex();
    }
    return true;
  };

 private:
};
}; // namespace paircleaner
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_PAIRCLEANER_H_
