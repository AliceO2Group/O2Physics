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

#include <TAxis.h>
#include <TCollection.h>
#include <TH2.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <string>
#include <utility>
#include <vector>

class FFitWeights : public TNamed
{
 public:
  FFitWeights();
  explicit FFitWeights(const char* name);
  ~FFitWeights();

  void init();
  void fillWeights(float centrality, float qn, int nh, const char* pf = "");
  void fillPt(float centrality, float pt, float weight, bool first);
  float getPtMult(float centrality);
  TObjArray* getDataArray() { return fW_data; }

  void setCentBin(int bin) { centBin = bin; }
  void setBinAxis(int bin, float min, float max)
  {
    qAxis = new TAxis(bin, min, max);
  }
  TAxis* getqVecAx() { return qAxis; }

  void setPtBin(int bin) { ptBin = bin; }
  void setPtAxis(int bin, float min, float max)
  {
    ptAxis = new TAxis(bin, min, max);
  }
  TAxis* getPtAx() { return ptAxis; }

  Long64_t Merge(TCollection* collist);
  void qSelection(const std::vector<int>& nhv, const std::vector<std::string>& stv);
  float eval(float centr, const float& dqn, const int nh, const char* pf = "");
  float evalPt(float centr, const float& mpt);
  void setResolution(int res) { nResolution = res; }
  int getResolution() const { return nResolution; }
  void setQnType(const std::vector<std::pair<int, std::string>>& qninp) { qnTYPE = qninp; }

  void mptSel();

 private:
  TObjArray* fW_data;
  TProfile* ptProfCent; //!
  TH2D* h2ptCent;       //!

  int centBin;
  TAxis* qAxis; //!
  int nResolution;
  int ptBin;
  TAxis* ptAxis; //!

  std::vector<std::pair<int, std::string>> qnTYPE;

  const char* getQName(const int nh, const char* pf = "")
  {
    return Form("q%i%s", nh, pf);
  };
  const char* getAxisName(const int nh, const char* pf = "")
  {
    return Form(";Centrality;q_{%i}^{%s}", nh, pf);
  };
  void addArray(TObjArray* targ, TObjArray* sour);

  static constexpr int NumberSp = 90;
  static constexpr float MaxTol = 100.05;

  float internalEval(float centr, const float& val, const char* name);

  ClassDef(FFitWeights, 1); // calibration class
};
#endif // COMMON_CORE_FFITWEIGHTS_H_
