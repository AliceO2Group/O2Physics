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

#include <TArrayF.h>
#include <TMathBase.h>
#include <TNamed.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <vector>
using namespace std;

ClassImp(MixingHandler);

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
void MixingHandler::AddMixingVariable(int var, std::vector<float> binLims)
{
  fVariables[var] = fVariableLimits.size();
  fVariableLimits.push_back(binLims);
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
  for (auto [var, pos] : fVariables) {
    // check that the value is within limits, if not return -1 to exclude the event from mixing
    size_t binValue = std::distance(fVariableLimits[pos].begin(), std::upper_bound(fVariableLimits[pos].begin(), fVariableLimits[pos].end(), values[var]));
    if (binValue == 0 || binValue == fVariableLimits[pos].size()) {
      return -1; // all variables must be inside limits
    }
    bin.push_back(binValue - 1);
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
        tempCategory *= (fVariableLimits[iv2].size() - 1);
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

  // Search for the position of the variable "var" in the internal variable list of the handler
  int ivar = fVariables.at(var);

  // extract the bin position in variable "var" from the category
  int norm = 1;
  for (int i = fVariables.size() - 1; i > ivar; --i) {
    norm *= (fVariableLimits[i].size() - 1);
  }
  int truncatedCategory = category - (category % norm);
  truncatedCategory /= norm;
  return truncatedCategory % (fVariableLimits[ivar].size() - 1);
}
