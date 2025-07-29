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

/// \file   GFWWeightsList.h
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Dec/25/2024
/// \brief  one object to hold a list of GFWWeights objects,

#ifndef PWGCF_GENERICFRAMEWORK_CORE_GFWWEIGHTSLIST_H_
#define PWGCF_GENERICFRAMEWORK_CORE_GFWWEIGHTSLIST_H_
#include <map>
#include <cstdio>
#include <string>
#include <vector>

#include "Framework/Logger.h"

#include "TObjArray.h"
#include "GFWWeights.h"

class GFWWeightsList : public TNamed
{
 public:
  GFWWeightsList();
  explicit GFWWeightsList(const char* name);
  ~GFWWeightsList();
  void init(const char* listName);
  void addGFWWeightsByName(const char* weightName, int nPtBins, double* ptBins, bool addData = kTRUE, bool addMC = kTRUE);
  GFWWeights* getGFWWeightsByName(const char* weightName);
  void addGFWWeightsByRun(int runNumber, int nPtBins, double* ptBins, bool addData = kTRUE, bool addMC = kTRUE);
  GFWWeights* getGFWWeightsByRun(int runNumber);
  void addPIDGFWWeightsByName(const char* weightName, int nPtBins, double* ptBins, double ptrefup, bool addData = kTRUE, bool addMC = kTRUE);
  GFWWeights* getPIDGFWWeightsByName(const char* weightName, int pidIndex);
  void addPIDGFWWeightsByRun(int runNumber, int nPtBins, double* ptBins, double ptrefup, bool addData = kTRUE, bool addMC = kTRUE);
  GFWWeights* getPIDGFWWeightsByRun(int runNumber, int pidIndex);
  void printRuns()
  {
    for (auto& el : runNumberPIDMap)
      printf("%i\n", el.first);
  }

  TObjArray* getList() const { return list; }
  Long64_t Merge(TCollection* collist);

 private:
  TObjArray* list;
  std::vector<std::string> species = {"_ref", "_ch", "_pi", "_ka", "_pr"}; //!
  std::map<int, GFWWeights*> runNumberMap;
  std::map<int, std::vector<GFWWeights*>> runNumberPIDMap;
  void addArray(TObjArray* target, TObjArray* source);

  ClassDef(GFWWeightsList, 1);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_GFWWEIGHTSLIST_H_
