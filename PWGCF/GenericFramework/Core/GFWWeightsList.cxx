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
#include "GFWWeightsList.h"

GFWWeightsList::GFWWeightsList() : TNamed("", ""), list(0)
{
  runNumberMap.clear();
  runNumberPIDMap.clear();
}

GFWWeightsList::GFWWeightsList(const char* name) : TNamed(name, name), list(0)
{
  runNumberMap.clear();
  runNumberPIDMap.clear();
}

GFWWeightsList::~GFWWeightsList()
{
  delete list;
  runNumberMap.clear();
  runNumberPIDMap.clear();
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
  weight->setPtBins(nPtBins, ptBins);
  weight->init(addData, addMC);
  list->Add(weight);
}

GFWWeights* GFWWeightsList::getGFWWeightsByName(const char* weightName)
{
  if (!list) {
    LOGF(error, "weight list is not initialized\n");
    return nullptr;
  }
  return reinterpret_cast<GFWWeights*>(list->FindObject(weightName));
}

void GFWWeightsList::addGFWWeightsByRun(int runNumber, int nPtBins, double* ptBins, bool addData, bool addMC)
{
  if (!list) {
    init("weightList");
  }
  if (runNumberMap.contains(runNumber)) {
    return;
  }
  GFWWeights* weight = new GFWWeights(Form("weight_%d", runNumber));
  weight->setPtBins(nPtBins, ptBins);
  weight->init(addData, addMC);
  list->Add(weight);
  runNumberMap.insert(std::make_pair(runNumber, weight));
}

GFWWeights* GFWWeightsList::getGFWWeightsByRun(int runNumber)
{
  if (!list) {
    LOGF(error, "weight list is not initialized\n");
    return nullptr;
  }
  if (!runNumberMap.contains(runNumber)) {
    LOGF(error, "weight for run %d is not found\n", runNumber);
    return nullptr;
  }
  return runNumberMap.at(runNumber);
}

void GFWWeightsList::addPIDGFWWeightsByName(const char* weightName, int nPtBins, double* ptBins, double ptrefup, bool addData, bool addMC)
{
  if (!list) {
    init("weightList");
  }

  std::vector<double> ptbins(ptBins, ptBins + nPtBins + 1);
  auto it = std::find(ptbins.begin(), ptbins.end(), ptrefup);
  std::vector<double> refpt(ptbins.begin(), it + 1);

  for (auto& type : species) {
    if (reinterpret_cast<GFWWeights*>(list->FindObject((static_cast<std::string>(weightName) + type).c_str()))) {
      continue;
    }
    GFWWeights* weight = new GFWWeights(Form("%s", (static_cast<std::string>(weightName) + type).c_str()));
    if (!type.compare("_ref"))
      weight->setPtBins(refpt.size() - 1, &(refpt[0]));
    else
      weight->setPtBins(nPtBins, ptBins);
    weight->init(addData, addMC);
    list->Add(weight);
  }
}
GFWWeights* GFWWeightsList::getPIDGFWWeightsByName(const char* weightName, int pidIndex)
{
  if (static_cast<size_t>(pidIndex) >= species.size())
    return nullptr;
  if (!list) {
    LOGF(error, "weight list is not initialized\n");
    return nullptr;
  }
  return reinterpret_cast<GFWWeights*>(list->FindObject((static_cast<std::string>(weightName) + species[pidIndex]).c_str()));
}
void GFWWeightsList::addPIDGFWWeightsByRun(int runNumber, int nPtBins, double* ptBins, double ptrefup, bool addData, bool addMC)
{
  if (!list) {
    init("weightList");
  }

  if (runNumberPIDMap.contains(runNumber))
    return;
  std::vector<double> ptbins(ptBins, ptBins + nPtBins + 1);
  auto it = std::find(ptbins.begin(), ptbins.end(), ptrefup);
  std::vector<double> refpt(ptbins.begin(), it + 1);

  std::vector<GFWWeights*> weights;
  for (auto& type : species) {
    GFWWeights* weight = new GFWWeights(Form("weight_%d%s", runNumber, type.c_str()));
    if (!type.compare("_ref"))
      weight->setPtBins(refpt.size() - 1, &(refpt[0]));
    else
      weight->setPtBins(nPtBins, ptBins);
    weight->init(addData, addMC);
    list->Add(weight);
    weights.push_back(weight);
  }
  LOGF(info, "Adding weights for run %d\n", runNumber);
  runNumberPIDMap.insert(std::make_pair(runNumber, weights));
  return;
}

GFWWeights* GFWWeightsList::getPIDGFWWeightsByRun(int runNumber, int pidIndex)
{
  if (!list) {
    LOGF(error, "weight list is not initialized\n");
    return nullptr;
  }
  if (!runNumberPIDMap.contains(runNumber)) {
    LOGF(error, "PID weights for run %d is not found\n", runNumber);
    return nullptr;
  }
  return runNumberPIDMap.at(runNumber)[pidIndex];
}
Long64_t GFWWeightsList::Merge(TCollection* collist)
{
  Long64_t nmerged = 0;
  if (!list) {
    list = new TObjArray();
    list->SetName("weightList");
    list->SetOwner(kTRUE);
  }
  TIter allWeights(collist);
  GFWWeightsList* lWeight = 0;
  while ((lWeight = (reinterpret_cast<GFWWeightsList*>(allWeights())))) {
    addArray(list, lWeight->getList());
    nmerged++;
  }
  return nmerged;
}
void GFWWeightsList::addArray(TObjArray* target, TObjArray* source)
{
  if (!source) {
    return;
  }
  for (int i = 0; i < source->GetEntries(); i++) {
    GFWWeights* sourw = reinterpret_cast<GFWWeights*>(source->At(i));
    GFWWeights* targw = reinterpret_cast<GFWWeights*>(target->FindObject(sourw->GetName()));
    if (!targw) {
      targw = reinterpret_cast<GFWWeights*>(sourw->Clone(sourw->GetName()));
      target->Add(targw);
    } else {
      targw->mergeWeights(sourw);
    }
  }
};
