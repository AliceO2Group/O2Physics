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

#include <iostream>
#include <fstream>
using namespace std;

#include <TMath.h>
#include <TTimeStamp.h>
#include <TRandom.h>

ClassImp(MixingHandler);

//_________________________________________________________________________
MixingHandler::MixingHandler() : TNamed(),
                                 fIsInitialized(kFALSE),
                                 fVariableLimits(),
                                 fVariables(),
                                 fVariableNames(),
                                 fNMixingVariables(0)
{
  //
  // default constructor
  //
  float dummyRange[2] = {-99999., +99999.};
  for (int iVar = 0; iVar < kNMaxVariables; ++iVar) {
    fVariableLimits[iVar].Set(2, dummyRange);
    fVariables[iVar] = VarManager::kNothing;
    fVariableNames[iVar] = "";
  }
}

//_________________________________________________________________________
MixingHandler::MixingHandler(const char* name, const char* title) : TNamed(name, title),
                                                                    fIsInitialized(kFALSE),
                                                                    fVariableLimits(),
                                                                    fVariables(),
                                                                    fNMixingVariables(0)
{
  //
  // Named constructor
  //
  float dummyRange[2] = {-99999., +99999.};
  for (int iVar = 0; iVar < kNMaxVariables; ++iVar) {
    fVariableLimits[iVar].Set(2, dummyRange);
    fVariables[iVar] = VarManager::kNothing;
    fVariableNames[iVar] = "";
  }
}

//_________________________________________________________________________
MixingHandler::~MixingHandler()
{
  //
  // destructor
  //
}

//_________________________________________________________________________
void MixingHandler::AddMixingVariable(int var, int nBins, float* binLims, TString varName)
{
  //
  // add a mixing variable
  //
  if (fNMixingVariables >= kNMaxVariables) {
    cout << "MixingHandler::AddMixingVariable(): ERROR Too many variables for the mixing!" << endl;
    cout << "                  Maximum number of variables: " << kNMaxVariables << endl;
    return;
  }
  fVariables[fNMixingVariables] = var;
  fVariableLimits[fNMixingVariables].Set(nBins, binLims);
  fVariableNames[fNMixingVariables] += varName;
  fNMixingVariables++;
  VarManager::SetUseVariable(var);
}

//_________________________________________________________________________
void MixingHandler::AddMixingVariable(int var, int nBins, std::vector<float> binLims, TString varName)
{

  float* bins = new float[nBins];
  for (int i = 0; i < nBins; ++i) {
    bins[i] = binLims[i];
  }
  AddMixingVariable(var, nBins, bins, varName);
}

//_________________________________________________________________________
int MixingHandler::GetMixingVariable(TString vars)
{
  int varNum = -1;
  for (int iVar = 0; iVar < fNMixingVariables; ++iVar) {
    if (!fVariableNames[iVar].CompareTo(vars)) {
      varNum = iVar;
    }
  }
  return varNum;
}

//_________________________________________________________________________
std::vector<float> MixingHandler::GetMixingVariableLimits(TString vars)
{
  std::vector<float> binLimits;
  for (int iVar = 0; iVar < fNMixingVariables; ++iVar) {
    if (!fVariableNames[iVar].CompareTo(vars)) {
      for (int iBin = 0; iBin < fVariableLimits[iVar].GetSize(); ++iBin) {
        binLimits.push_back(fVariableLimits[iVar].At(iBin));
      }
    }
  }
  return binLimits;
}

//_________________________________________________________________________
void MixingHandler::Init()
{
  //
  // Initialization of pools
  // NOTE: The master array holding tracks is a 1D array but it will be represented as an n-dim array
  //       The size of the array will be N_1 x N_2 x ... x N_n
  //       The correct event category will be retrieved using the function FindEventCategory()
  //
  int size = 1;
  for (int iVar = 0; iVar < fNMixingVariables; ++iVar) {
    size *= (fVariableLimits[iVar].GetSize() - 1);
  }

  fIsInitialized = kTRUE;
}

//_________________________________________________________________________
int MixingHandler::FindEventCategory(float* values)
{
  //
  // Find the event category corresponding to the centrality, vtxz and ep values
  //
  if (fNMixingVariables == 0) {
    return -1;
  }
  if (!fIsInitialized) {
    Init();
  }

  int bin[kNMaxVariables];
  for (int i = 0; i < fNMixingVariables; ++i) {
    bin[i] = TMath::BinarySearch(fVariableLimits[i].GetSize(), fVariableLimits[i].GetArray(), values[fVariables[i]]);
    if (bin[i] == -1 || bin[i] == fVariableLimits[i].GetSize() - 1) {
      return -1; // all variables must be inside limits
    }
  }

  int category = 0;
  for (int iVar = 0; iVar < fNMixingVariables; ++iVar) {
    int tempCategory = 1;
    for (int iVar2 = iVar; iVar2 < fNMixingVariables; ++iVar2) {
      if (iVar2 == iVar) {
        tempCategory *= bin[iVar2];
      } else {
        tempCategory *= (fVariableLimits[iVar2].GetSize() - 1);
      }
    }
    category += tempCategory;
  }
  return category;
}

//_________________________________________________________________________
int MixingHandler::GetBinFromCategory(int iVar, int category) const
{
  //
  // find the bin in variable var for the n-dimensional "category"
  //
  if (fNMixingVariables == 0) {
    return -1;
  }
  int norm = 1;
  for (int i = fNMixingVariables - 1; i > iVar; --i) {
    norm *= (fVariableLimits[i].GetSize() - 1);
  }
  int truncatedCategory = category - (category % norm);
  truncatedCategory /= norm;
  return truncatedCategory % (fVariableLimits[iVar].GetSize() - 1);
}
