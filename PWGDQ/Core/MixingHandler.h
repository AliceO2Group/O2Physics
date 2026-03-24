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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
// Class to define and fill histograms
//

#ifndef PWGDQ_CORE_MIXINGHANDLER_H_
#define PWGDQ_CORE_MIXINGHANDLER_H_

#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/VarManager.h"

#include <TArrayF.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TNamed.h>
#include <TString.h>

#include <vector>

class MixingHandler : public TNamed
{

 public:
  MixingHandler();
  MixingHandler(const char* name, const char* title);
  virtual ~MixingHandler();

  // setters
  void AddMixingVariable(int var, int nBins, float* binLims);
  void AddMixingVariable(int var, int nBins, std::vector<float> binLims);

  // getters
  int GetNMixingVariables() const { return fVariables.size(); }
  int GetMixingVariable(VarManager::Variables var); // returns the position in the internal varible list of the handler. Useful for checks, mostly
  std::vector<float> GetMixingVariableLimits(VarManager::Variables var);

  void Init();
  int FindEventCategory(float* values);
  int GetBinFromCategory(VarManager::Variables var, int category) const;

 private:
  MixingHandler(const MixingHandler& handler);
  MixingHandler& operator=(const MixingHandler& handler);

  // User options
  bool fIsInitialized; // check if the mixing handler is initialized

  std::vector<TArrayF> fVariableLimits;
  std::vector<int> fVariables;

  ClassDef(MixingHandler, 1);
};

#endif // PWGDQ_CORE_MIXINGHANDLER_H_
