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

#ifndef PWGCF_GENERICFRAMEWORK_GFWWEIGHTS_H_
#define PWGCF_GENERICFRAMEWORK_GFWWEIGHTS_H_
#include "TObjArray.h"
#include "TNamed.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCollection.h"
#include "TString.h"

class GFWWeights : public TNamed
{
 public:
  GFWWeights();
  explicit GFWWeights(const char* name);
  ~GFWWeights();
  void Init(bool AddData = kTRUE, bool AddM = kTRUE);
  void Fill(double phi, double eta, double vz, double pt, double cent, int htype, double weight = 1); // htype: 0 for data, 1 for mc rec, 2 for mc gen
  double GetWeight(double phi, double eta, double vz, double pt, double cent, int htype);             // htype: 0 for data, 1 for mc rec, 2 for mc gen
  double GetNUA(double phi, double eta, double vz);                                                   // This just fetches correction from integrated NUA, should speed up
  double GetNUE(double pt, double eta, double vz);                                                    // fetches weight from fEffInt
  bool IsDataFilled() { return fDataFilled; }
  bool IsMCFilled() { return fMCFilled; }
  double FindMax(TH3D* inh, int& ix, int& iy, int& iz);
  void MCToEfficiency();
  TObjArray* GetRecArray() { return fW_mcrec; }
  TObjArray* GetGenArray() { return fW_mcgen; }
  TObjArray* GetDataArray() { return fW_data; }
  void CreateNUA(bool IntegrateOverCentAndPt = kTRUE);
  void CreateNUE(bool IntegrateOverCentrality = kTRUE);
  TH1D* GetIntegratedEfficiencyHist();
  bool CalculateIntegratedEff();
  double GetIntegratedEfficiency(double pt);
  void SetDataFilled(bool newval) { fDataFilled = newval; }
  void SetMCFilled(bool newval) { fMCFilled = newval; }
  void ReadAndMerge(TString filelinks, TString listName = "OutputList", bool addData = kTRUE, bool addRec = kTRUE, bool addGen = kTRUE);
  void SetPtBins(int Nbins, double* bins);
  Long64_t Merge(TCollection* collist);
  void RebinNUA(int nX = 1, int nY = 2, int nZ = 5);
  void OverwriteNUA();
  TH1D* GetdNdPhi();
  TH1D* GetEfficiency(double etamin, double etamax, double vzmin, double vzmax);

 private:
  bool fDataFilled;
  bool fMCFilled;
  TObjArray* fW_data;
  TObjArray* fW_mcrec;
  TObjArray* fW_mcgen;
  TH3D* fEffInt;   //!
  TH1D* fIntEff;   //!
  TH3D* fAccInt;   //!
  int fNbinsPt;    //! do not store
  double* fbinsPt; //! do not store
  void AddArray(TObjArray* targ, TObjArray* sour);
  const char* GetBinName(double ptv, double v0mv, const char* pf = "")
  {
    int ptind = 0;  // GetPtBin(ptv);
    int v0mind = 0; // GetV0MBin(v0mv);
    return Form("Bin%s_weights%i_%i", pf, ptind, v0mind);
  };

  ClassDef(GFWWeights, 1);
};

#endif // PWGCF_GENERICFRAMEWORK_GFWWEIGHTS_H_
