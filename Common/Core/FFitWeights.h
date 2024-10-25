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

/// \file FFitWeights.h
/// \brief Class for handling fit weights for ESE framework. Hold methods for loading and calculating all ESE splines. Supports FT0C, in the future it will support FT0A, FV0A and TPC.
///
/// \author Joachim C. K. B. Hansen

#ifndef COMMON_CORE_FFITWEIGHTS_H_
#define COMMON_CORE_FFITWEIGHTS_H_

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <complex>
#include <memory>

#include "TNamed.h"
#include "TObjArray.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCollection.h"
#include "TString.h"
#include "TMath.h"

class FFitWeights : public TNamed
{
 public:
  FFitWeights();
  explicit FFitWeights(const char* name);
  ~FFitWeights();

  void Init();
  void Fill(float centrality, float qn, int nh, const char* pf = "");
  TObjArray* GetDataArray() { return fW_data; }

  void SetCentBin(int bin) { CentBin = bin; }
  void SetBinAxis(int bin, float min, float max)
  {
    qAxis = new TAxis(bin, min, max);
  }
  TAxis* GetqVecAx() { return qAxis; }

  // Long64_t Merge(TCollection* collist);
  void qSelectionSpline(std::vector<int> nhv, std::vector<std::string> stv);
  float EvalSplines(float centr, const float& dqn, const int nh, const char* pf = "");
  void SetResolution(int res) { nResolution = res; }
  void SetQnType(std::vector<std::pair<int, std::string>> qninp) { qnTYPE = qninp; }

 private:
  TObjArray* fW_data;

  int CentBin;
  TAxis* qAxis; //!
  int nResolution;

  std::vector<std::pair<int, std::string>> qnTYPE; 

  const char* GetQName(const int nh, const char* pf = "")
  {
    return Form("q%i%s", nh, pf);
  };
  const char* GetAxisName(const int nh, const char* pf = "")
  {
    return Form(";Centrality;q_{%i}^{%s}", nh, pf);
  };
  void AddArray(TObjArray* targ, TObjArray* sour);

  ClassDef(FFitWeights, 1); // calibration class
};
#endif // COMMON_CORE_FFITWEIGHTS_H_
