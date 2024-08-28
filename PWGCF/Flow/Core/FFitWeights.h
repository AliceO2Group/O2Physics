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

/*
Author: Joachim Hansen

*/
#ifndef PWGCF_FLOW_CORE_FFITWEIGHTS_H_
#define PWGCF_FLOW_CORE_FFITWEIGHTS_H_

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <complex>

#include "TNamed.h"
#include "TObjArray.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCollection.h"
#include "TString.h"
#include "TMath.h"

struct BinRange {
  BinRange(int b, int r) : bins{b}, range{r} {};
  int bins{0};
  int range{0};
};

class FFitWeights : public TNamed
{
 public:
  FFitWeights();
  explicit FFitWeights(const char* name);
  ~FFitWeights();

  void Init();
  void FillFT0(std::size_t iCh, float amplitude, float GainCst);
  void FillQ(float mult, float vec, int nHarm, const char* coord, const char* qType);
  TObjArray* GetDataArray() { return fW_data; }
  // double GetGain(double phi, double eta, double vz);
  void CreateGain();
  std::vector<float> GetGain();
  void SetChIdBin(int bin) { ChIDBin = bin; }
  void SetCentBin(int bin) { CentBin = bin; }
  void SetAmplBin(int bin, int Range)
  {
    sAmpl.bins = bin;
    sAmpl.range = Range;
  }
  void SetqVBin(int bin, int Range)
  {
    sqVec.bins = bin;
    sqVec.range = Range;
  }
  void SetqCorVBin(int bin, int Range)
  {
    sqCorVec.bins = bin;
    sqCorVec.range = Range;
  }

  void CreateRecenter(const char* xy);
  float GetRecVal(int cent, const char* xy, const int nHarm);
  void CreateRMS(const char* xy);
  float GetRMSVal(int cent, const char* xy, const int nHarm);

  template <typename CollType>
  static float EventPlane(const CollType& coll, int nHarm)
  {
    auto x = coll.qFT0CRe();
    auto y = coll.qFT0CIm();
    return 1 / nHarm * TMath::ATan2(y[nHarm - 2], x[nHarm - 2]);
  }
  // static float EventPlane(const float& x, const float& y, const float& nHarm);

 private:
  TObjArray* fW_data;
  std::vector<float> vGain;

  int CentBin;
  int ChIDBin;

  BinRange sAmpl;
  BinRange sqVec;
  BinRange sqCorVec;

  // TH2F *FT0Ampl;

  ClassDef(FFitWeights, 1); // calibration class
};
#endif // PWGCF_FLOW_CORE_FFITWEIGHTS_H_
