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

/// \file   GFWWeightsList.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Dec/25/2024
/// \brief  one object to hold a list of GFWWeights objects,

#include <utility>
#include <cstdio>
#include "GFWWeightsList.h"

GFWWeightsList::GFWWeightsList() : TNamed("", ""), list(0)
{
  runNumerMap.clear();
}

GFWWeightsList::GFWWeightsList(const char* name) : TNamed(name, name), list(0)
{
  runNumerMap.clear();
}

GFWWeightsList::~GFWWeightsList()
{
  delete list;
  runNumerMap.clear();
}

void GFWWeightsList::init(const char* listName)
{
  list = new TObjArray();
  list->SetName(listName);
  list->SetOwner(kTRUE);
}

void GFWWeightsList::addGFWWeightsByName(const char* weightName, int nPtBins, double* ptBins, bool addData, bool addMC)
{
  if (!list) {
    init("weightList");
  }
  if (reinterpret_cast<GFWWeights*>(list->FindObject(weightName))) {
    return;
  }
  GFWWeights* weight = new GFWWeights(weightName);
  weight->SetPtBins(nPtBins, ptBins);
  weight->Init(addData, addMC);
  list->Add(weight);
}

GFWWeights* GFWWeightsList::getGFWWeightsByName(const char* weightName)
{
  if (!list) {
    printf("Error: weight list is not initialized\n");
    return nullptr;
  }
  return reinterpret_cast<GFWWeights*>(list->FindObject(weightName));
}

void GFWWeightsList::addGFWWeightsByRun(int runNumber, int nPtBins, double* ptBins, bool addData, bool addMC)
{
  if (!list) {
    init("weightList");
  }
  if (runNumerMap.contains(runNumber)) {
    return;
  }
  GFWWeights* weight = new GFWWeights(Form("weight_%d", runNumber));
  weight->SetPtBins(nPtBins, ptBins);
  weight->Init(addData, addMC);
  list->Add(weight);
  runNumerMap.insert(std::make_pair(runNumber, weight));
}

GFWWeights* GFWWeightsList::getGFWWeightsByRun(int runNumber)
{
  if (!list) {
    printf("Error: weight list is not initialized\n");
    return nullptr;
  }
  if (!runNumerMap.contains(runNumber)) {
    printf("Error: weight for run %d is not found\n", runNumber);
    return nullptr;
  }
  return runNumerMap.at(runNumber);
}
