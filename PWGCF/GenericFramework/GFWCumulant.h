/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.)
A part of <GFW.cxx/h>
A container to store Q vectors for one subevent with an extra layer to recursively calculate particle correlations.
If used, modified, or distributed, please aknowledge the author of this code.
*/
#ifndef ALIGFWCUMULANT__H
#define ALIGFWCUMULANT__H
#include <cmath>
#include <complex>
#include <vector>
using std::complex;
using std::vector;
class GFWCumulant
{
 public:
  GFWCumulant();
  ~GFWCumulant();
  void ResetQs();
  void FillArray(int ptin, double phi, double weight = 1, double SecondWeight = -1);
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
  bool IsPtBinFilled(int ptb);
  void CreateComplexVectorArray(int N = 1, int P = 1, int Pt = 1);
  void CreateComplexVectorArrayVarPower(int N = 1, vector<int> Pvec = {1}, int Pt = 1);
  int PW(int ind) { return fPowVec.at(ind); }; // No checks to speed up, be carefull!!!
  void DestroyComplexVectorArray();
  complex<double> Vec(int, int, int ptbin = 0); // envelope class to summarize pt-dif. Q-vec getter
 protected:
  complex<double>*** fQvector;
  uint fUsed;
  int fNEntries;
  // Q-vectors. Could be done recursively, but maybe defining each one of them explicitly is easier to read
  int fN;              //! Harmonics
  int fPow;            //! Power
  vector<int> fPowVec; //! Powers array
  int fPt;             //! fPt bins
  bool* fFilledPts;
  bool fInitialized; // Arrays are initialized
  complex<double> fNullQ = 0;
};

#endif
