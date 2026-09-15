// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GFWPowerArray.h
/// \brief Class to compute necessary powers of Q-vectors based on input correlations
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#ifndef PWGCF_GENERICFRAMEWORK_CORE_GFWPOWERARRAY_H_
#define PWGCF_GENERICFRAMEWORK_CORE_GFWPOWERARRAY_H_

#include <cstddef>
#include <vector>

typedef std::vector<int> HarSet;
class GFWPowerArray
{
 public:
  static HarSet GetPowerArray(std::vector<HarSet> inHarmonics); // o2-linter: disable=name/function-variable (preserve existing public API)
  static void PowerArrayTest();                                 // o2-linter: disable=name/function-variable (preserve existing public API)

 private:
  static int getHighestHarmonic(const HarSet& inhar);
  static HarSet trimVec(const HarSet& hars, int ind);
  static HarSet addConstant(const HarSet& hars, int offset);
  static void flushVectorToMaster(HarSet& masterVector, const HarSet& comVec, int maxPower);
  static void recursiveFunction(HarSet& masterVector, const HarSet& hars, int offset, int maxPower, std::size_t startIndex);
  static void printVector(const HarSet& singleSet);
};
#endif // PWGCF_GENERICFRAMEWORK_CORE_GFWPOWERARRAY_H_
