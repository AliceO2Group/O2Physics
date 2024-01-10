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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_FLOWPTCONTAINER_H_
#define PWGCF_GENERICFRAMEWORK_CORE_FLOWPTCONTAINER_H_

#include <algorithm>
#include <vector>
#include "BootstrapProfile.h"
#include "TNamed.h"
#include "TList.h"
#include "TCollection.h"
#include "Framework/HistogramSpec.h"
#include "GFW.h"
#include "GFWConfig.h"

namespace o2::analysis::genericframework::eventweight
{
enum kEventWeight {
  kUnity,
  kTuples
};
};

using namespace o2::analysis::genericframework;
using namespace o2::analysis::genericframework::eventweight;

class FlowPtContainer : public TNamed
{
 public:
  FlowPtContainer();
  explicit FlowPtContainer(const char* name);
  ~FlowPtContainer();
  FlowPtContainer(const char* name, const char* title, int nbinsx, double* xbins, const int& m, const GFWCorrConfigs& configs);
  FlowPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, const int& m, const GFWCorrConfigs& configs);
  void Initialise(const o2::framework::AxisSpec axis, const int& m, const GFWCorrConfigs& configs, const int& nsub = 10);
  void Initialise(int nbinsx, double* xbins, const int& m, const GFWCorrConfigs& configs, const int& nsub = 10);
  void Initialise(int nbinsx, double xlow, double xhigh, const int& m, const GFWCorrConfigs& configs, const int& nsub = 10);
  void Fill(const double& w, const double& pt);
  int GetVectorIndex(const int i, const int j) { return j * (mpar + 1) + i; }
  void CalculateCorrelations();
  void CalculateCMTerms();
  void FillPtProfiles(const Double_t& lMult, const Double_t& rn);
  void FillVnPtProfiles(const double& lMult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask);
  void FillCMProfiles(const double& lMult, const double& rn);
  TList* GetCorrList() { return fCorrList; }
  TList* GetCMTermList() { return fCMTermList; }
  void SetEventWeight(const unsigned int& lWeight) { fEventWeight = lWeight; }
  void RebinMulti(Int_t nbins);
  void RebinMulti(Int_t nbins, double* binedges);
  TH1* getCentralMomentHist(int ind, int m);
  TH1* getCumulantHist(int ind, int m);
  TH1* getCorrHist(int ind, int m);
  Int_t getMpar() { return mpar; }
  Long64_t Merge(TCollection* collist);
  Double_t OrderedAddition(std::vector<double> vec);
  void CreateCentralMomentList();
  void CalculateCentralMomentHists(std::vector<TH1*> inh, int ind, int m, TH1* hMpt);
  void CreateCumulantList();
  void CalculateCumulantHists(std::vector<TH1*> inh, Int_t ind);
  void ClearVector()
  {
    sumP.clear();
    sumP.resize((mpar + 1) * (mpar + 1));
    fillCounter = 0;
  };

 private:
  TList* fCMTermList;
  TList* fCorrList;
  TList* fCovList;
  TList* fCumulantList;
  TList* fCentralMomentList;
  int mpar;
  int fillCounter;
  unsigned int fEventWeight;
  void MergeBSLists(TList* source, TList* target);
  TH1* raiseHistToPower(TH1* inh, double p);
  std::vector<double> sumP;    //!
  std::vector<double> corrNum; //!
  std::vector<double> corrDen; //!

  static constexpr float fFactorial[9] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320.};
  static constexpr int fSign[9] = {1, -1, 1, -1, 1, -1, 1, -1, 1};
  ClassDef(FlowPtContainer, 1);
};
#endif // PWGCF_GENERICFRAMEWORK_CORE_FLOWPTCONTAINER_H_
