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

#include <fstream>
#include <iostream>
using namespace std;

#include <TMath.h>
#include <TRandom.h>
#include <TTimeStamp.h>

ClassImp(MixingHandler);

//_________________________________________________________________________
MixingHandler::MixingHandler() : TNamed(),
                                 fIsInitialized(kFALSE),
                                 fVariableLimits(),
                                 fVariables()
{
  //
  // default constructor
  //
}

//_________________________________________________________________________
MixingHandler::MixingHandler(const char* name, const char* title) : TNamed(name, title),
                                                                    fIsInitialized(kFALSE),
                                                                    fVariableLimits(),
                                                                    fVariables()
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
void MixingHandler::AddMixingVariable(int var, int nBins, float* binLims)
{
  //
  // add a mixing variable
  //
  fVariables.push_back(var);
  TArrayF varBins;
  varBins.Set(nBins, binLims);
  fVariableLimits.push_back(varBins);
  VarManager::SetUseVariable(var);
}

//_________________________________________________________________________
void MixingHandler::AddMixingVariable(int var, int nBins, std::vector<float> binLims)
{

  float* bins = new float[nBins];
  for (int i = 0; i < nBins; ++i) {
    bins[i] = binLims[i];
  }
  AddMixingVariable(var, nBins, bins);
}

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
}

//_________________________________________________________________________
void MixingHandler::Init()
{
  //
  // Initialization of pools
  //       The correct event category will be retrieved using the function FindEventCategory()
  //
  int size = 1;
  for (auto v : fVariableLimits) {
    size *= (v.GetSize() - 1);
  }
  (void)size;
  fIsInitialized = kTRUE;
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

  std::vector<int> bin;
  int iVar = 0;
  for (auto v = fVariableLimits.begin(); v != fVariableLimits.end(); v++, iVar++) {
    int binValue = TMath::BinarySearch((*v).GetSize(), (*v).GetArray(), values[fVariables[iVar]]);
    bin.push_back(binValue);
    if (bin[iVar] == -1 || bin[iVar] == (*v).GetSize() - 1) {
      return -1; // all variables must be inside limits
    }
  }

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
        tempCategory *= (fVariableLimits[iv2].GetSize() - 1);
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
  int tempVar = 0;
  for (auto v = fVariables.begin(); v != fVariables.end(); v++, tempVar++) {
    if (*v == var) {
      break;
    }
  }

  // extract the bin position in variable "var" from the category
  int norm = 1;
  for (int i = fVariables.size() - 1; i > tempVar; --i) {
    norm *= (fVariableLimits[i].GetSize() - 1);
  }
  int truncatedCategory = category - (category % norm);
  truncatedCategory /= norm;
  return truncatedCategory % (fVariableLimits[tempVar].GetSize() - 1);
}
