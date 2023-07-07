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
//
//
// Class to produce smeared pt,eta,phi

#ifndef PWGEM_DILEPTON_UTILS_MOMENTUMSMEARER_H_
#define PWGEM_DILEPTON_UTILS_MOMENTUMSMEARER_H_

#include <TH1D.h>
#include <TString.h>
#include <TGrid.h>
#include <TObjArray.h>
#include <TFile.h>
#include "Framework/Logger.h"

class MomentumSmearer
{
 public:
  /// Default constructor
  MomentumSmearer() = default;

  /// Constructor with resolution histograms
  MomentumSmearer(TString resFileName, TString resPtHistName, TString resEtaHistName, TString resPhiPosHistName, TString resPhiNegHistName)
  {
    setResFileName(resFileName);
    setResPtHistName(resPtHistName);
    setResEtaHistName(resEtaHistName);
    setResPhiPosHistName(resPhiPosHistName);
    setResPhiNegHistName(resPhiNegHistName);
    init();
  }

  /// Default destructor
  ~MomentumSmearer() = default;

  void init()
  {
    if (fInitialized)
      return;

    if (fResFileName.BeginsWith("alien://")) {
      TGrid::Connect("alien://");
    }

    // get resolution histo
    LOGP(info, "Set Resolution histo");
    // Get Resolution map
    TFile* fFile = TFile::Open(fResFileName);
    if (!fFile) {
      LOGP(error, "Could not open Resolution file {}", fResFileName);
      return;
    }
    TObjArray* ArrResoPt = nullptr;
    if (fFile->GetListOfKeys()->Contains(fResPtHistName)) {
      ArrResoPt = reinterpret_cast<TObjArray*>(fFile->Get(fResPtHistName));
    } else {
      LOGP(error, "Could not open {} from file {}", fResPtHistName, fResFileName);
    }

    TObjArray* ArrResoEta = nullptr;
    if (fFile->GetListOfKeys()->Contains(fResEtaHistName)) {
      ArrResoEta = reinterpret_cast<TObjArray*>(fFile->Get(fResEtaHistName));
    } else {
      LOGP(error, "Could not open {} from file {}", fResEtaHistName, fResFileName);
    }

    TObjArray* ArrResoPhi_Pos = nullptr;
    if (fFile->GetListOfKeys()->Contains(TString(fResPhiPosHistName))) {
      ArrResoPhi_Pos = reinterpret_cast<TObjArray*>(fFile->Get(fResPhiPosHistName));
    } else {
      LOGP(error, "Could not open {} from file {}", fResPhiPosHistName, fResFileName);
    }

    TObjArray* ArrResoPhi_Neg = nullptr;
    if (fFile->GetListOfKeys()->Contains(TString(fResPhiNegHistName))) {
      ArrResoPhi_Neg = reinterpret_cast<TObjArray*>(fFile->Get(fResPhiNegHistName));
    } else {
      LOGP(error, "Could not open {} from file {}", fResPhiNegHistName, fResFileName);
    }

    fArrResoPt = ArrResoPt;
    fArrResoEta = ArrResoEta;
    fArrResoPhi_Pos = ArrResoPhi_Pos;
    fArrResoPhi_Neg = ArrResoPhi_Neg;
    fFile->Close();

    fInitialized = true;
  }

  void applySmearing(const int ch, const float ptgen, const float etagen, const float phigen, float& ptsmeared, float& etasmeared, float& phismeared)
  {
    // smear pt
    int ptbin = reinterpret_cast<TH2D*>(fArrResoPt->At(0))->GetXaxis()->FindBin(ptgen);
    if (ptbin < 1) {
      ptbin = 1;
    }
    if (ptbin > fArrResoPt->GetLast()) {
      ptbin = fArrResoPt->GetLast();
    }
    float smearing = 0.;
    TH1D* thisHist = reinterpret_cast<TH1D*>(fArrResoPt->At(ptbin));
    if (thisHist->GetEntries() > 0) {
      smearing = thisHist->GetRandom() * ptgen;
    }
    ptsmeared = ptgen - smearing;

    // smear eta
    ptbin = reinterpret_cast<TH2D*>(fArrResoEta->At(0))->GetXaxis()->FindBin(ptgen);
    if (ptbin < 1) {
      ptbin = 1;
    }
    if (ptbin > fArrResoEta->GetLast()) {
      ptbin = fArrResoEta->GetLast();
    }
    smearing = 0.;
    thisHist = reinterpret_cast<TH1D*>(fArrResoEta->At(ptbin));
    if (thisHist->GetEntries() > 0) {
      smearing = thisHist->GetRandom();
    }
    etasmeared = etagen - smearing;

    // smear phi
    ptbin = reinterpret_cast<TH2D*>(fArrResoPhi_Pos->At(0))->GetXaxis()->FindBin(ptgen);
    if (ptbin < 1) {
      ptbin = 1;
    }
    if (ptbin > fArrResoPhi_Pos->GetLast()) {
      ptbin = fArrResoPhi_Pos->GetLast();
    }
    smearing = 0.;
    if (ch < 0) {
      thisHist = reinterpret_cast<TH1D*>(fArrResoPhi_Neg->At(ptbin));
    } else {
      thisHist = reinterpret_cast<TH1D*>(fArrResoPhi_Pos->At(ptbin));
    }
    if (thisHist->GetEntries() > 0) {
      smearing = thisHist->GetRandom();
    }
    phismeared = phigen - smearing;
  }

  // setters
  void setResFileName(TString resFileName) { fResFileName = resFileName; }
  void setResPtHistName(TString resPtHistName) { fResPtHistName = resPtHistName; }
  void setResEtaHistName(TString resEtaHistName) { fResEtaHistName = resEtaHistName; }
  void setResPhiPosHistName(TString resPhiPosHistName) { fResPhiPosHistName = resPhiPosHistName; }
  void setResPhiNegHistName(TString resPhiNegHistName) { fResPhiNegHistName = resPhiNegHistName; }

  // getters
  TString getResFileName() { return fResFileName; }
  TString getResPtHistName() { return fResPtHistName; }
  TString getResEtaHistName() { return fResEtaHistName; }
  TString getResPhiPosHistName() { return fResPhiPosHistName; }
  TString getResPhiNegHistName() { return fResPhiNegHistName; }
  TObjArray* getArrResoPt() { return fArrResoPt; }
  TObjArray* getArrResoEta() { return fArrResoEta; }
  TObjArray* getArrResoPhiPos() { return fArrResoPhi_Pos; }
  TObjArray* getArrResoPhiNeg() { return fArrResoPhi_Neg; }

 private:
  bool fInitialized = false;
  TString fResFileName;
  TString fResPtHistName;
  TString fResEtaHistName;
  TString fResPhiPosHistName;
  TString fResPhiNegHistName;
  TObjArray* fArrResoPt;
  TObjArray* fArrResoEta;
  TObjArray* fArrResoPhi_Pos;
  TObjArray* fArrResoPhi_Neg;
};

#endif // PWGEM_DILEPTON_UTILS_MOMENTUMSMEARER_H_
