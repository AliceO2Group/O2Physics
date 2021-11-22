/*
Author: Vytautas Vislavicius
Ported to O2: Emil Gorm Nielsen
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
*/
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
  void SetType(uint infl)
  {
    DestroyComplexVectorArray();
    fUsed = infl;
  };
  void Inc() { fNEntries++; };
  int GetN() { return fNEntries; };
  // protected:
  TComplex*** fQvector;
  uint fUsed;
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
  int PW(int ind) { return fPowVec.at(ind); }; // No checks to speed up, be carefull!!!
  void DestroyComplexVectorArray();
  bool IsPtBinFilled(int ptb)
  {
    if (!fFilledPts)
      return kFALSE;
    return fFilledPts[ptb];
  };
};

#endif