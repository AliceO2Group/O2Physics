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

/// \file GFWPowerArray.cxx
/// \brief Compute required powers of Q-vectors for harmonic sets
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#include "GFWPowerArray.h"

#include <Framework/Logger.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

int GFWPowerArray::getHighestHarmonic(const HarSet& inhar)
{
  // Highest possible harmonic: sum of same-sign harmonics
  std::int64_t maxPos = 0, maxNeg = 0;
  for (const int& val : inhar) {
    if (val > 0)
      maxPos += val;
    else
      maxNeg += std::abs(static_cast<std::int64_t>(val));
    if (maxPos >= std::numeric_limits<int>::max() || maxNeg >= std::numeric_limits<int>::max())
      throw std::overflow_error("Harmonic sum exceeds the supported range");
  }
  return static_cast<int>(std::max(maxPos, maxNeg));
};
HarSet GFWPowerArray::trimVec(const HarSet& hars, int ind)
{
  HarSet retVec = hars;
  retVec.erase(retVec.begin() + ind);
  return retVec;
};
HarSet GFWPowerArray::addConstant(const HarSet& hars, int offset)
{
  HarSet retVec = hars;
  for (int& val : retVec) // o2-linter: disable=const-ref-in-for-loop (updates each harmonic)
    val += offset;
  return retVec;
};
void GFWPowerArray::flushVectorToMaster(HarSet& masterVector, const HarSet& comVec, int maxPower)
{
  int nPartLoc = maxPower - static_cast<int>(comVec.size()) + 1;
  for (const int& val : comVec) {
    int absVal = std::abs(val);
    if (masterVector.at(absVal) < nPartLoc) {
      masterVector.at(absVal) = nPartLoc;
    }
  }
};
void GFWPowerArray::recursiveFunction(HarSet& masterVector, const HarSet& hars, int offset, int maxPower, std::size_t startIndex)
{
  HarSet compVec = addConstant(hars, offset);
  flushVectorToMaster(masterVector, compVec, maxPower);
  for (std::size_t i = startIndex; i < hars.size(); i++)
    recursiveFunction(masterVector, trimVec(hars, static_cast<int>(i)), offset + hars.at(i), maxPower, i);
};
void GFWPowerArray::printVector(const HarSet& singleSet)
{
  if (singleSet.empty()) {
    LOGF(info, "Vector is empty!");
    return;
  }
  std::string output = "{" + std::to_string(singleSet[0]);
  for (std::size_t i = 1; i < singleSet.size(); i++)
    output += ", " + std::to_string(singleSet[i]);
  output += "}";
  LOGF(info, "%s", output.c_str());
}
HarSet GFWPowerArray::GetPowerArray(const std::vector<HarSet>& inHarmonics) // o2-linter: disable=name/function-variable (preserve existing public API)
{
  // First, find maximum number of particle correlations ( = max power) and maximum (sum of) harmonics
  int maxHar = 0;
  int maxPart = 0;
  for (const HarSet& singleSet : inHarmonics) {
    int harSum = getHighestHarmonic(singleSet);
    maxHar = harSum > maxHar ? harSum : maxHar;
  }
  // Make a vector with maxHar+1 entries (entry 0 for sum=0)
  HarSet retVec = HarSet(maxHar + 1);
  // Then loop over all combinations and calculate max powers
  for (const HarSet& singleSet : inHarmonics) {
    if (singleSet.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::overflow_error("Harmonic set exceeds the supported particle count");
    int lNPart = static_cast<int>(singleSet.size()); // Total number of particles correlated
    recursiveFunction(retVec, singleSet, 0, lNPart, 0);
    // Harmonic sum = 0 is a special case. In principle all 0 cases with non-zero harmonics are captured by the function above, but to calculate normalization, we set all harmonics to 0. This means that sum=0 power is the max number of harmonics/particles being correlated
    maxPart = (lNPart > maxPart) ? lNPart : maxPart;
  }
  // Override the sum=0 power with the number of correlated particles
  if (retVec[0] < maxPart)
    retVec[0] = maxPart;
  // Need an extra power ( = 0) for all non-zero powers
  for (int& val : retVec) // o2-linter: disable=const-ref-in-for-loop (increments the required powers)
    if (val != 0)
      val++;
  return retVec;
};
void GFWPowerArray::PowerArrayTest() // o2-linter: disable=name/function-variable (preserve existing public API)
{
  std::vector<HarSet> allHars = {
    HarSet{2},
    HarSet{3},
    HarSet{2, 2},
    HarSet{3, 3}};
  LOGF(info, "Input harmonics are:");
  for (const HarSet& inSet : allHars)
    printVector(inSet);
  LOGF(info, "The configuration of powers must then be:");
  auto vc = GetPowerArray(allHars);
  printVector(vc);
};
