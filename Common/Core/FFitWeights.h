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
/// \brief Class for handling fit weights. Right now holds FT0, will hold methods for loading and calculating all ESE splines in the future.
///
/// \author Joachim C. K. B. Hansen, Lund University

#ifndef COMMON_CORE_FFITWEIGHTS_H_
#define COMMON_CORE_FFITWEIGHTS_H_

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <complex>
#include <memory>

#include "TNamed.h"
#include "TObjArray.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCollection.h"
#include "TString.h"
#include "TMath.h"

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

  void SetBinAxis(int bin, float min, float max, int axisLevel)
  {
    if (axisLevel == 0)
      sAmpl = std::shared_ptr<TAxis>(new TAxis(bin, min, max));
    else if (axisLevel == 1)
      sqVec = std::shared_ptr<TAxis>(new TAxis(bin, min, max));
    else if (axisLevel == 2)
      sqCorVec = std::shared_ptr<TAxis>(new TAxis(bin, min, max));
    else
      printf("something went wrong assigning axes");
  }
  std::shared_ptr<TAxis> GetAmplAx() { return sAmpl; }
  std::shared_ptr<TAxis> GetqVecAx() { return sqVec; }
  std::shared_ptr<TAxis> GetqCorVecAx() { return sqCorVec; }

  void CreateRecenter(const char* xy);
  float GetRecVal(int cent, const char* xy, const int nHarm);
  TH1F* GetRecHist(const char* xy, const int nHarm);
  void CreateRMS(const char* xy);
  float GetRMSVal(int cent, const char* xy, const int nHarm);
  TH1F* GetRmsHist(const char* xy, const int nHarm);

  template <typename CollType>
  static float EventPlane(const CollType& coll, float nHarm)
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

  std::shared_ptr<TAxis> sAmpl;    //!
  std::shared_ptr<TAxis> sqVec;    //!
  std::shared_ptr<TAxis> sqCorVec; //!

  // TH2F *FT0Ampl;

  ClassDef(FFitWeights, 2); // calibration class
};
#endif // COMMON_CORE_FFITWEIGHTS_H_
