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

#ifndef MixingHandler_H
#define MixingHandler_H

#include <TNamed.h>
#include <TArrayF.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TString.h>

#include "PWGDQ/Core/HistogramManager.h"

class MixingHandler : public TNamed
{

 public:
  enum Constants { // event mixing for correlations
    kNMaxVariables = 10
  };

 public:
  MixingHandler();
  MixingHandler(const char* name, const char* title);
  virtual ~MixingHandler();

  // setters
  void AddMixingVariable(int var, int nBins, float* binLims, TString varName);
  void AddMixingVariable(int var, int nBins, std::vector<float> binLims, TString varName);

  // getters
  int GetNMixingVariables() const { return fNMixingVariables; }
  int GetMixingVariable(TString vars);
  std::vector<float> GetMixingVariableLimits(TString vars);

  void Init();
  int FindEventCategory(float* values);
  int GetBinFromCategory(int iVar, int category) const;

 private:
  MixingHandler(const MixingHandler& handler);
  MixingHandler& operator=(const MixingHandler& handler);

  // User options
  bool fIsInitialized; // check if the mixing handler is initialized

  TArrayF fVariableLimits[kNMaxVariables];
  int fVariables[kNMaxVariables];
  TString fVariableNames[kNMaxVariables]; //! variable names
  int fNMixingVariables;

  ClassDef(MixingHandler, 4);
};

#endif
