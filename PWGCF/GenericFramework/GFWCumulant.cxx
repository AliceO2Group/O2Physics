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

#include "GFWCumulant.h"

GFWCumulant::GFWCumulant() : fQvector(0),
                             fUsed(kBlank),
                             fNEntries(-1),
                             fN(1),
                             fPow(1),
                             fPt(1),
                             fFilledPts(0),
                             fInitialized(kFALSE){};

GFWCumulant::~GFWCumulant(){
  // printf("Destructor (?) for some reason called?\n");
  // DestroyComplexVectorArray();
};
void GFWCumulant::FillArray(double eta, int ptin, double phi, double weight, double SecondWeight)
{
  if (!fInitialized)
    CreateComplexVectorArray(1, 1, 1);
  if (fPt == 1)
    ptin = 0; // If one bin, then just fill it straight; otherwise, if ptin is out-of-range, do not fill
  else if (ptin < 0 || ptin >= fPt)
    return;
  fFilledPts[ptin] = kTRUE;
  for (int lN = 0; lN < fN; lN++) {
    double lSin = TMath::Sin(lN * phi); // No need to recalculate for each power
    double lCos = TMath::Cos(lN * phi); // No need to recalculate for each power
    for (int lPow = 0; lPow < PW(lN); lPow++) {
      double lPrefactor = 0;
      // Dont calculate it twice; multiplication is cheaper that power
      // Also, if second weight is specified, then keep the first weight with power no more than 1, and us the other weight otherwise
      // this is important when POIs are a subset of REFs and have different weights than REFs
      if (SecondWeight > 0 && lPow > 1)
        lPrefactor = TMath::Power(SecondWeight, lPow - 1) * weight;
      else
        lPrefactor = TMath::Power(weight, lPow);
      double qsin = lPrefactor * lSin;
      double qcos = lPrefactor * lCos;
      fQvector[ptin][lN][lPow](fQvector[ptin][lN][lPow].Re() + qcos, fQvector[ptin][lN][lPow].Im() + qsin); //+=TComplex(qcos,qsin);
    };
  };
  Inc();
};
void GFWCumulant::ResetQs()
{
  if (!fNEntries)
    return; // If 0 entries, then no need to reset. Otherwise, if -1, then just initialized and need to set to 0.
  for (int i = 0; i < fPt; i++) {
    fFilledPts[i] = kFALSE;
    for (int lN = 0; lN < fN; lN++) {
      for (int lPow = 0; lPow < PW(lN); lPow++) {
        fQvector[i][lN][lPow](0., 0.);
      };
    };
  };
  fNEntries = 0;
};
void GFWCumulant::DestroyComplexVectorArray()
{
  if (!fInitialized)
    return;
  for (int l_n = 0; l_n < fN; l_n++) {
    for (int i = 0; i < fPt; i++) {
      delete[] fQvector[i][l_n];
    };
  };
  for (int i = 0; i < fPt; i++) {
    delete[] fQvector[i];
  };
  delete[] fQvector;
  delete[] fFilledPts;
  fInitialized = kFALSE;
  fNEntries = -1;
};

void GFWCumulant::CreateComplexVectorArray(int N, int Pow, int Pt)
{
  DestroyComplexVectorArray();
  vector<int> pwv;
  for (int i = 0; i < N; i++)
    pwv.push_back(Pow);
  CreateComplexVectorArrayVarPower(N, pwv, Pt);
};
void GFWCumulant::CreateComplexVectorArrayVarPower(int N, vector<int> PowVec, int Pt)
{
  DestroyComplexVectorArray();
  fN = N;
  fPow = 0;
  fPt = Pt;
  fFilledPts = new bool[Pt];
  fPowVec = PowVec;
  fQvector = new TComplex**[fPt];
  for (int i = 0; i < fPt; i++) {
    fQvector[i] = new TComplex*[fN];
  };
  for (int l_n = 0; l_n < fN; l_n++) {
    for (int i = 0; i < fPt; i++) {
      fQvector[i][l_n] = new TComplex[PW(l_n)];
    };
  };
  ResetQs();
  fInitialized = kTRUE;
};
TComplex GFWCumulant::Vec(int n, int p, int ptbin)
{
  if (!fInitialized)
    return 0;
  if (ptbin >= fPt || ptbin < 0)
    ptbin = 0;
  if (n >= 0)
    return fQvector[ptbin][n][p];
  return TComplex::Conjugate(fQvector[ptbin][-n][p]);
};