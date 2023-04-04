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

#ifndef GFWCUMULANT__H
#define GFWCUMULANT__H
#include "TComplex.h"
#include "TNamed.h"
#include "TMath.h"
#include "TAxis.h"
using std::vector;
class GFWCumulant
{
 public:
  GFWCumulant();
  ~GFWCumulant();
  void ResetQs();
  void FillArray(double eta, int ptin, double phi, double weight = 1, double SecondWeight = -1);
  enum UsedFlags_t { kBlank = 0,
                     kFull = 1,
                     kPt = 2 };
  void SetType(unsigned int infl)
  {
    DestroyComplexVectorArray();
    fUsed = infl;
  };
  void Inc() { fNEntries++; };
  int GetN() { return fNEntries; };
  // protected:
  TComplex*** fQvector;
  unsigned int fUsed;
  int fNEntries;
  // Q-vectors. Could be done recursively, but maybe defining each one of them explicitly is easier to read
  TComplex Vec(int, int, int ptbin = 0); // envelope class to summarize pt-dif. Q-vec getter
  int fN;                                //! Harmonics
  int fPow;                              //! Power
  vector<int> fPowVec;                   //! Powers array
  int fPt;                               //! fPt bins
  bool* fFilledPts;
  bool fInitialized; // Arrays are initialized
  void CreateComplexVectorArray(int N = 1, int P = 1, int Pt = 1);
  void CreateComplexVectorArrayVarPower(int N = 1, vector<int> Pvec = {1}, int Pt = 1);
  int PW(int ind) { return fPowVec.at(ind); }; // No checks to speed up, be careful!!!
  void DestroyComplexVectorArray();
  bool IsPtBinFilled(int ptb)
  {
    if (!fFilledPts)
      return kFALSE;
    return fFilledPts[ptb];
  };
};

#endif