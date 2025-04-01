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

/// \file GFWWeights.h
/// \brief Class to store corrections for the Generic Framework
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#ifndef PWGCF_GENERICFRAMEWORK_CORE_GFWWEIGHTS_H_
#define PWGCF_GENERICFRAMEWORK_CORE_GFWWEIGHTS_H_

#include "Framework/Logger.h"

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
  void init(bool AddData = kTRUE, bool AddM = kTRUE);
  void fill(double phi, double eta, double vz, double pt, double cent, int htype, double weight = 1); // htype: 0 for data, 1 for mc rec, 2 for mc gen
  double getWeight(double phi, double eta, double vz, double pt, double cent, int htype);             // htype: 0 for data, 1 for mc rec, 2 for mc gen
  double getNUA(double phi, double eta, double vz);                                                   // This just fetches correction from integrated NUA, should speed up
  double getNUE(double pt, double eta, double vz);                                                    // fetches weight from fEffInt
  bool isDataFilled() { return fDataFilled; }
  bool isMCFilled() { return fMCFilled; }
  double findMax(TH3D* inh, int& ix, int& iy, int& iz);
  void mcToEfficiency();
  TObjArray* getRecArray() { return fW_mcrec; }
  TObjArray* getGenArray() { return fW_mcgen; }
  TObjArray* getDataArray() { return fW_data; }
  void createNUA(bool IntegrateOverCentAndPt = kTRUE);
  void createNUE(bool IntegrateOverCentrality = kTRUE);
  TH1D* getIntegratedEfficiencyHist();
  bool calculateIntegratedEff();
  double getIntegratedEfficiency(double pt);
  void setDataFilled(bool newval) { fDataFilled = newval; }
  void setMCFilled(bool newval) { fMCFilled = newval; }
  void readAndMerge(TString filelinks, TString listName = "OutputList", bool addData = kTRUE, bool addRec = kTRUE, bool addGen = kTRUE);
  void setPtBins(int Nbins, double* bins);
  Long64_t Merge(TCollection* collist);
  void rebinNUA(int nX = 1, int nY = 2, int nZ = 5);
  void overwriteNUA();
  TH1D* getdNdPhi();
  TH1D* getEfficiency(double etamin, double etamax, double vzmin, double vzmax);
  void mergeWeights(GFWWeights* other);
  void setTH3D(TH3D* th3d);

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
  void addArray(TObjArray* targ, TObjArray* sour);
  const char* getBinName(double /*ptv*/, double /*v0mv*/, const char* pf = "")
  {
    int ptind = 0;  // GetPtBin(ptv);
    int v0mind = 0; // GetV0MBin(v0mv);
    return Form("Bin%s_weights%i_%i", pf, ptind, v0mind);
  };

  ClassDef(GFWWeights, 1);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_GFWWEIGHTS_H_
