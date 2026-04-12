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

/// \file GFW.h/.cxx
/// \brief Container to store correlations and compute common cumulants
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#ifndef PWGCF_GENERICFRAMEWORK_CORE_FLOWCONTAINER_H_
#define PWGCF_GENERICFRAMEWORK_CORE_FLOWCONTAINER_H_
#include <vector>
#include "TH3F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNamed.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "TString.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TString.h"
#include "TCollection.h"
#include "TAxis.h"
#include "ProfileSubset.h"
#include "Framework/HistogramSpec.h"

class FlowContainer : public TNamed
{
 public:
  FlowContainer();
  explicit FlowContainer(const char* name);
  ~FlowContainer();
  enum StatisticsType { kSingleSample,
                        kJackKnife,
                        kBootstrap };
  void Initialize(TObjArray* inputList, const o2::framework::AxisSpec axis, int nRandomized = 0);
  void Initialize(TObjArray* inputList, int nMultiBins, double MultiMin, double MultiMax, int nRandomized = 0);
  bool CreateBinsFromAxis(TAxis* inax);
  void SetXAxis(TAxis* inax);
  void SetXAxis();
  void RebinMulti(int rN)
  {
    if (fProf)
      fProf->RebinX(rN);
  };
  int GetNMultiBins() { return fProf->GetNbinsX(); }
  double GetMultiAtBin(int bin) { return fProf->GetXaxis()->GetBinCenter(bin); }
  int FillProfile(const char* hname, double multi, double y, double w, double rn);
  TProfile2D* GetProfile() { return fProf; }
  void OverrideProfileErrors(TProfile2D* inpf);
  void ReadAndMerge(const char* infile);
  void PickAndMerge(TFile* tfi);
  bool OverrideBinsWithZero(int xb1, int yb1, int xb2, int yb2);
  bool OverrideMainWithSub(int subind, bool ExcludeChosen);
  bool RandomizeProfile(int nSubsets = 0);
  bool CreateStatisticsProfile(StatisticsType StatType, int arg);
  TObjArray* GetSubProfiles() { return fProfRand; }
  Long64_t Merge(TCollection* collist);
  void SetIDName(TString newname); //! do not store
  void SetPtRebin(int newval) { fPtRebin = newval; }
  void SetPtRebin(int nbins, double* binedges);
  void SetMultiRebin(int nbins, double* binedges);
  double* GetMultiRebin(int& nBins);
  void SetPropagateErrors(bool newval) { fPropagateErrors = newval; }
  TProfile* GetCorrXXVsMulti(const char* order, int l_pti = 0);                             // pti = 0 for pt-integrated
  TH1D* GetCorrXXVsPt(const char* order, double lminmulti = -1, double lmaxmulti = -1);     // 0 for multi. integrated
  TH1D* GetHistCorrXXVsMulti(const char* order, int l_pti = 0);                             // pti = 0 for pt-integrated
  TH1D* GetHistCorrXXVsPt(const char* order, double lminmulti = -1, double lmaxmulti = -1); // 0 for multi. integrated

  TH1D* GetVN2VsMulti(int n = 2, int l_pta = 0) { return GetVN2VsX(n, kFALSE, l_pta); }
  TH1D* GetVN2VsPt(int n = 2, double min = -1, double max = -1) { return GetVN2VsX(n, kTRUE, min, max); }
  TH1D* GetCN4VsMulti(int n = 2, int pti = 0) { return GetCN4VsX(n, kFALSE, pti); }
  TH1D* GetCN4VsPt(int n = 2, double min = -1, double max = -1) { return GetCN4VsX(n, kTRUE, min, max); }

  TH1D* GetVN4VsMulti(int n = 2, int pti = 0) { return GetVN4VsX(n, kFALSE, pti); }
  TH1D* GetVN4VsPt(int n = 2, double min = -1, double max = -1) { return GetVN4VsX(n, kTRUE, min, max); }

  TH1D* GetVN6VsMulti(int n = 2, int pti = 0) { return GetVN6VsX(n, kFALSE, pti); }
  TH1D* GetVN6VsPt(int n = 2, double min = -1, double max = -1) { return GetVN6VsX(n, kTRUE, min, max); }

  TH1D* GetVN8VsMulti(int n = 2, int pti = 0) { return GetVN8VsX(n, kFALSE, pti); }
  TH1D* GetVN8VsPt(int n = 2, double min = -1, double max = -1) { return GetVN8VsX(n, kTRUE, min, max); }

  TH1D* GetCNN(int n = 2, int c = 2, bool onPt = kTRUE, double arg1 = -1, double arg2 = -1);
  TH1D* GetVNN(int n = 2, int c = 2, bool onPt = kTRUE, double arg1 = -1, double arg2 = -1);

  // private:

  double CN2Value(double cor2);  // This is redundant, but adding for completeness
  double CN2Error(double cor2e); // Also redundant
  double VN2Value(double cor2);
  double VN2Error(double cor2, double cor2e);
  double VDN2Value(double cor2d, double cor2);
  double VDN2Error(double cor2d, double cor2de, double cor2, double cor2e);

  double CN4Value(double cor4, double cor2);
  double CN4Error(double cor4e, double cor2, double cor2e);
  double DN4Value(double cor4d, double cor2d, double cor2);
  double DN4Error(double cor4de, double cor2d, double cor2de, double cor2, double cor2e);
  double VN4Value(double c4);
  double VN4Error(double c4, double c4e);
  double VDN4Value(double d4, double c4);
  double VDN4Error(double d4, double d4e, double c4, double c4e);

  double CN6Value(double cor6, double cor4, double cor2);
  double CN6Error(double cor6e, double cor4, double cor4e, double cor2, double cor2e);

  double DN6Value(double cor6d, double cor4d, double cor2d, double cor4, double cor2);
  double DN6Error(double d6e, double d4, double d4e, double d2,
                  double d2e, double c4, double c4e, double c2,
                  double c2e);
  double VN6Value(double c6);
  double VN6Error(double c6, double c6e);
  double VDN6Value(double d6, double c6);
  double VDN6Error(double d6, double d6e, double c6, double c6e);

  double CN8Value(double cor8, double cor6, double cor4, double cor2);
  double CN8Error(double cor8e, double cor6, double cor6e,
                  double cor4, double cor4e, double cor2, double cor2e);
  double DN8Value(double cor8d, double cor6d, double cor4d, double cor2d, double cor6, double cor4, double cor2);
  double DN8Error(double d8e, double d6, double d6e, double d4,
                  double d4e, double d2, double d2e, double c6,
                  double c6e, double c4, double c4e, double c2,
                  double c2e);
  double VN8Value(double c8);
  double VN8Error(double c8, double c8e);
  double VDN8Value(double d8, double c8);
  double VDN8Error(double d8, double d8e, double c8, double c8e);

  TH1D* GetCN2VsX(int n = 2, bool onPt = kTRUE, double larg1 = -1, double larg2 = -1); // This one is redundant
  TH1D* GetVN2VsX(int n = 2, bool onPt = kTRUE, double larg1 = -1, double larg2 = -1);
  TH1D* GetCN4VsX(int n = 2, bool onPt = kTRUE, double larg1 = -1, double larg2 = -1);
  TH1D* GetVN4VsX(int n = 2, bool onPt = kTRUE, double larg1 = -1, double larg2 = -1);
  TH1D* GetCN6VsX(int n = 2, bool onPt = kTRUE, double larg1 = -1, double larg2 = -1);
  TH1D* GetVN6VsX(int n = 2, bool onPt = kTRUE, double larg1 = -1, double larg2 = -1);
  TH1D* GetCN8VsX(int n = 2, bool onPt = kTRUE, double larg1 = -1, double larg2 = -1);
  TH1D* GetVN8VsX(int n = 2, bool onPt = kTRUE, double larg1 = -1, double larg2 = -1);

  TH1D* GetVN2(TH1D* cn2);
  TH1D* GetVN4(TH1D* inh);
  TH1D* GetVN6(TH1D* inh);
  TH1D* GetVN8(TH1D* inh);
  TH1D* GetCN2(TH1D* corrN2);
  TH1D* GetCN4(TH1D* corrN4, TH1D* corrN2);
  TH1D* GetCN6(TH1D* corrN6, TH1D* corrN4, TH1D* corrN2);
  TH1D* GetCN8(TH1D* corrN8, TH1D* corrN6, TH1D* corrN4, TH1D* corrN2);
  TH1D* ProfToHist(TProfile* inpf);
  TProfile2D* fProf;
  TObjArray* fProfRand;
  int fNRandom;
  TString fIDName;
  int fPtRebin;             //! do not store
  double* fPtRebinEdges;    //! do not store
  int fMultiRebin;          //! do not store
  double* fMultiRebinEdges; //! do not store
  TAxis* fXAxis;
  int fNbinsPt;          //! Do not store; stored in the fXAxis
  double* fbinsPt;       //! Do not store; stored in fXAxis
  bool fPropagateErrors; //! do not store
  TProfile* GetRefFlowProfile(const char* order, double m1 = -1, double m2 = -1);
  ClassDef(FlowContainer, 2);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_FLOWCONTAINER_H_
