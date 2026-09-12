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

#include "PWGDQ/Core/MixingHandler.h"

#include "PWGDQ/Core/VarManager.h"

#include <TNamed.h>

#include <Rtypes.h>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <vector>
using namespace std;

//_________________________________________________________________________
MixingHandler::MixingHandler() : TNamed(),
                                 fIsInitialized(false),
                                 fVariableLimits(),
                                 fVariables(),
                                 fPoolDepth(0),
                                 fPools()
{
  //
  // default constructor
  //
}

//_________________________________________________________________________
MixingHandler::MixingHandler(const char* name, const char* title) : TNamed(name, title),
                                                                    fIsInitialized(false),
                                                                    fVariableLimits(),
                                                                    fVariables(),
                                                                    fPoolDepth(0),
                                                                    fPools()
{
  //
  // Named constructor
  //
}

//_________________________________________________________________________
MixingHandler::~MixingHandler()
{
  //
  // destructor
  //
}

//_________________________________________________________________________
void MixingHandler::AddMixingVariable(int var, const std::vector<float>& binLims)
{
  fVariables[var] = fVariableLimits.size();
  fVariableLimits.push_back(binLims);
  // FillEvent() only fills variables marked as used
  VarManager::SetUseVariable(var);
}

/*
//_________________________________________________________________________
int MixingHandler::GetMixingVariable(VarManager::Variables var)
{
  int i = 0;
  for (auto v = fVariables.begin(); v != fVariables.end(); v++, i++) {
    if (*v == var) {
      return i;
    }
  }
  return -1;
}
*/

/*
//_________________________________________________________________________
std::vector<float> MixingHandler::GetMixingVariableLimits(VarManager::Variables var)
{
  std::vector<float> binLimits;
  int i = 0;
  for (auto v = fVariables.begin(); v != fVariables.end(); v++, i++) {
    if (*v == var) {
      for (int iBin = 0; iBin < fVariableLimits[i].GetSize(); ++iBin) {
        binLimits.push_back(fVariableLimits[i].At(iBin));
      }
      break;
    }
  }
  return binLimits;
}*/

//_________________________________________________________________________
void MixingHandler::Init()
{
  // loop over all variables and create a mixing pool for each category defined by the binning of the variables
  int nCategories = 1;
  for (auto& var : fVariables) {
    nCategories *= (fVariableLimits[var.second].size() - 1);
  }
  // add elements in the map for each category (the key is the category and the value is an empty pool)
  for (int i = 0; i < nCategories; i++) {
    fPools[i] = MixingPool();
  }
  fIsInitialized = true;
}

//_________________________________________________________________________
int MixingHandler::FindEventCategory(float* values)
{
  //
  // Find the event category corresponding to the added mixing variables
  //
  if (fVariables.size() == 0) {
    return -1;
  }
  if (!fIsInitialized) {
    Init();
  }

  // loop over the variables and find out in which bin the value of the variable for the event is located
  std::vector<int> bin;
  // number of bins per variable in the iteration order of fVariables (fVariableLimits is in insertion order)
  std::vector<int> nBins;
  for (auto [var, pos] : fVariables) {
    // check that the value is within limits, if not return -1 to exclude the event from mixing
    size_t binValue = std::distance(fVariableLimits[pos].begin(), std::upper_bound(fVariableLimits[pos].begin(), fVariableLimits[pos].end(), values[var]));
    if (binValue == 0 || binValue == fVariableLimits[pos].size()) {
      return -1; // all variables must be inside limits
    }
    bin.push_back(binValue - 1);
    nBins.push_back(fVariableLimits[pos].size() - 1);
  }

  // Hash the bin values to define a unique category
  // The hashing is done such that the original bin values can be retrieved from the category
  // For example, for 3 variables with n1, n2, n3 bins respectively, the category for bin values (b1, b2, b3) would be:
  // category = b1*(n2*n3) + b2*(n3) + b3
  int category = 0;
  int tempCategory = 1;
  int iv1 = 0;
  int iv2 = 0;
  for (auto v1 = fVariables.begin(); v1 != fVariables.end(); v1++, iv1++) {
    tempCategory = 1;
    iv2 = iv1;
    for (auto v2 = v1; v2 != fVariables.end(); v2++, iv2++) {
      if (iv2 == iv1) {
        tempCategory *= bin[iv2];
      } else {
        tempCategory *= nBins[iv2];
      }
    }
    category += tempCategory;
  }
  return category;
}

//_________________________________________________________________________
int MixingHandler::GetBinFromCategory(VarManager::Variables var, int category) const
{
  //
  // find the bin in variable var for the n-dimensional "category"
  //
  if (fVariables.size() == 0) {
    return -1;
  }

  // number of bins and position of var in the iteration order of fVariables, as used by FindEventCategory()
  std::vector<int> nBins;
  int ivar = -1;
  for (auto const& [v, pos] : fVariables) {
    if (v == var) {
      ivar = static_cast<int>(nBins.size());
    }
    nBins.push_back(fVariableLimits[pos].size() - 1);
  }
  if (ivar < 0) {
    return -1;
  }

  // extract the bin position in variable "var" from the category
  int norm = 1;
  for (size_t i = nBins.size() - 1; i > static_cast<size_t>(ivar); --i) {
    norm *= nBins[i];
  }
  int truncatedCategory = category - (category % norm);
  truncatedCategory /= norm;
  return truncatedCategory % nBins[ivar];
}

//_________________________________________________________________________
void MixingHandler::SetCategoryBinCenters(int category, float* values) const
{
  //
  // set the mixing variables to the bin centers of this category (used for the leftover mixing)
  //
  for (auto const& [var, pos] : fVariables) {
    int bin = GetBinFromCategory(static_cast<VarManager::Variables>(var), category);
    if (bin < 0) {
      continue;
    }
    values[var] = 0.5 * (fVariableLimits[pos][bin] + fVariableLimits[pos][bin + 1]);
  }
}
