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
    setResPhiNegHistName("");
    setResPhiNegHistName("");
    init();
  }

  /// Constructor with resolution histograms and efficiency
  MomentumSmearer(TString resFileName, TString resPtHistName, TString resEtaHistName, TString resPhiPosHistName, TString resPhiNegHistName, TString effFileName, TString effHistName)
  {
    setResFileName(resFileName);
    setResPtHistName(resPtHistName);
    setResEtaHistName(resEtaHistName);
    setResPhiPosHistName(resPhiPosHistName);
    setResPhiNegHistName(resPhiNegHistName);
    setResPhiNegHistName(resPhiNegHistName);
    setResPhiNegHistName(resPhiNegHistName);
    init();
  }

  /// Default destructor
  ~MomentumSmearer() = default;

  void init()
  {
    if (fInitialized)
      return;

    // ToDo: make it possible to access .root files from CCDB

    if (fResFileName.BeginsWith("alien://") || fEffFileName.BeginsWith("alien://")) {
      TGrid::Connect("alien://");
    }

    // get resolution histo
    if (fResFileName.CompareTo("") != 0) {
      LOGP(info, "Set Resolution histo");
      // Get Resolution map
      TFile* fFile = TFile::Open(fResFileName);
      if (!fFile) {
        LOGP(fatal, "Could not open Resolution file {}", fResFileName.Data());
        return;
      }
      TObjArray* ArrResoPt = nullptr;
      if (fFile->GetListOfKeys()->Contains(fResPtHistName)) {
        ArrResoPt = reinterpret_cast<TObjArray*>(fFile->Get(fResPtHistName));
      } else {
        LOGP(fatal, "Could not open {} from file {}", fResPtHistName.Data(), fResFileName.Data());
      }

      TObjArray* ArrResoEta = nullptr;
      if (fFile->GetListOfKeys()->Contains(fResEtaHistName)) {
        ArrResoEta = reinterpret_cast<TObjArray*>(fFile->Get(fResEtaHistName));
      } else {
        LOGP(fatal, "Could not open {} from file {}", fResEtaHistName.Data(), fResFileName.Data());
      }

      TObjArray* ArrResoPhi_Pos = nullptr;
      if (fFile->GetListOfKeys()->Contains(TString(fResPhiPosHistName))) {
        ArrResoPhi_Pos = reinterpret_cast<TObjArray*>(fFile->Get(fResPhiPosHistName));
      } else {
        LOGP(fatal, "Could not open {} from file {}", fResPhiPosHistName.Data(), fResFileName.Data());
      }

      TObjArray* ArrResoPhi_Neg = nullptr;
      if (fFile->GetListOfKeys()->Contains(TString(fResPhiNegHistName))) {
        ArrResoPhi_Neg = reinterpret_cast<TObjArray*>(fFile->Get(fResPhiNegHistName));
      } else {
        LOGP(fatal, "Could not open {} from file {}", fResPhiNegHistName.Data(), fResFileName.Data());
      }

      fArrResoPt = ArrResoPt;
      fArrResoEta = ArrResoEta;
      fArrResoPhi_Pos = ArrResoPhi_Pos;
      fArrResoPhi_Neg = ArrResoPhi_Neg;
      fFile->Close();
    }

    // get efficiency histo
    fEffType = 0;
    if (fEffFileName.CompareTo("") != 0) {
      LOGP(info, "Set Efficiency histo");
      TFile* fEffFile = TFile::Open(fEffFileName);
      if (!fEffFile) {
        LOGP(fatal, "Could not open efficiency file {}", fEffFileName.Data());
        return;
      }
      if (fEffFile->GetListOfKeys()->Contains(fEffHistName.Data())) {
        fArrEff = reinterpret_cast<TObject*>(fEffFile->Get(fEffHistName.Data()));
        // check which type is used
        if (dynamic_cast<TH3*>(fArrEff)) {
          fEffType = 3;
          LOGP(info, "Use 3d efficiency histo (pt, eta, phi)");
        } else if (dynamic_cast<TH2*>(fArrEff)) {
          fEffType = 2;
          LOGP(info, "Use 2d efficiency histo (pt, eta)");
        } else if (dynamic_cast<TH1*>(fArrEff)) {
          fEffType = 1;
          LOGP(info, "Use 1d efficiency histo (pt)");
        } else {
          LOGP(fatal, "Could not identify type of histogram {}", fEffHistName.Data());
        }
      } else {
        LOGP(fatal, "Could not find histogram {} in file {}", fEffHistName.Data(), fEffFileName.Data());
      }
    }

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

  float getEfficiency(float pt, float eta, float phi)
  {

    if (fEffType == 0) {
      return 1.;
    }

    if (fEffType == 1) {
      TH1F* hist = reinterpret_cast<TH1F*>(fArrEff);
      int ptbin = hist->GetXaxis()->FindBin(pt);
      int ptbin_max = hist->GetXaxis()->GetNbins();
      // make sure that no underflow or overflow bins are used
      if (ptbin < 1)
        ptbin = 1;
      else if (ptbin > ptbin_max)
        ptbin = ptbin_max;
      return hist->GetBinContent(ptbin);
    }

    if (fEffType == 2) {
      TH2F* hist = reinterpret_cast<TH2F*>(fArrEff);
      int ptbin = hist->GetXaxis()->FindBin(pt);
      int ptbin_max = hist->GetXaxis()->GetNbins();
      int etabin = hist->GetYaxis()->FindBin(eta);
      int etabin_max = hist->GetYaxis()->GetNbins();
      // make sure that no underflow or overflow bins are used
      if (ptbin < 1)
        ptbin = 1;
      else if (ptbin > ptbin_max)
        ptbin = ptbin_max;
      if (etabin < 1)
        etabin = 1;
      else if (etabin > etabin_max)
        etabin = etabin_max;
      return hist->GetBinContent(ptbin, etabin);
    }

    if (fEffType == 3) {
      TH3F* hist = reinterpret_cast<TH3F*>(fArrEff);
      int ptbin = hist->GetXaxis()->FindBin(pt);
      int ptbin_max = hist->GetXaxis()->GetNbins();
      int etabin = hist->GetYaxis()->FindBin(eta);
      int etabin_max = hist->GetYaxis()->GetNbins();
      int phibin = hist->GetZaxis()->FindBin(phi);
      int phibin_max = hist->GetZaxis()->GetNbins();
      // make sure that no underflow or overflow bins are used
      if (ptbin < 1)
        ptbin = 1;
      else if (ptbin > ptbin_max)
        ptbin = ptbin_max;
      if (etabin < 1)
        etabin = 1;
      else if (etabin > etabin_max)
        etabin = etabin_max;
      if (phibin < 1)
        phibin = 1;
      else if (phibin > phibin_max)
        phibin = phibin_max;
      return hist->GetBinContent(ptbin, etabin, phibin);
    }

    return 1.;
  }

  // setters
  void setResFileName(TString resFileName) { fResFileName = resFileName; }
  void setResPtHistName(TString resPtHistName) { fResPtHistName = resPtHistName; }
  void setResEtaHistName(TString resEtaHistName) { fResEtaHistName = resEtaHistName; }
  void setResPhiPosHistName(TString resPhiPosHistName) { fResPhiPosHistName = resPhiPosHistName; }
  void setResPhiNegHistName(TString resPhiNegHistName) { fResPhiNegHistName = resPhiNegHistName; }
  void setEffFileName(TString effFileName) { fEffFileName = effFileName; }
  void setEffHistName(TString effHistName) { fEffHistName = effHistName; }

  // getters
  TString getResFileName() { return fResFileName; }
  TString getResPtHistName() { return fResPtHistName; }
  TString getResEtaHistName() { return fResEtaHistName; }
  TString getResPhiPosHistName() { return fResPhiPosHistName; }
  TString getResPhiNegHistName() { return fResPhiNegHistName; }
  TString getEffFileName() { return fEffFileName; }
  TString getEffHistName() { return fEffHistName; }
  TObjArray* getArrResoPt() { return fArrResoPt; }
  TObjArray* getArrResoEta() { return fArrResoEta; }
  TObjArray* getArrResoPhiPos() { return fArrResoPhi_Pos; }
  TObjArray* getArrResoPhiNeg() { return fArrResoPhi_Neg; }
  TObject* getArrEff() { return fArrEff; }

 private:
  bool fInitialized = false;
  TString fResFileName;
  TString fResPtHistName;
  TString fResEtaHistName;
  TString fResPhiPosHistName;
  TString fResPhiNegHistName;
  TString fEffFileName;
  TString fEffHistName;
  int fEffType = 0;
  TObjArray* fArrResoPt;
  TObjArray* fArrResoEta;
  TObjArray* fArrResoPhi_Pos;
  TObjArray* fArrResoPhi_Neg;
  TObject* fArrEff;
};

#endif // PWGEM_DILEPTON_UTILS_MOMENTUMSMEARER_H_
