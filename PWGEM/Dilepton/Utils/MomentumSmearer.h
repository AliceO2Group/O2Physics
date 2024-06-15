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
#include <TKey.h>
#include <CCDB/BasicCCDBManager.h>
#include "Framework/Logger.h"

using namespace o2::framework;
using namespace o2;

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

  // template <typename T>
  void getResolutionHistos(TList* list)
  {

    fArrResoPt = reinterpret_cast<TObjArray*>(list->FindObject(fResPtHistName));
    if (!fArrResoPt) {
      LOGP(fatal, "Could not open {} from file {}", fResPtHistName.Data(), fResFileName.Data());
    }

    fArrResoEta = reinterpret_cast<TObjArray*>(list->FindObject(fResEtaHistName));
    if (!fArrResoEta) {
      LOGP(fatal, "Could not open {} from file {}", fResEtaHistName.Data(), fResFileName.Data());
    }

    fArrResoPhi_Pos = reinterpret_cast<TObjArray*>(list->FindObject(fResPhiPosHistName));
    if (!fArrResoPhi_Pos) {
      LOGP(fatal, "Could not open {} from file {}", fResPhiPosHistName.Data(), fResFileName.Data());
    }

    fArrResoPhi_Neg = reinterpret_cast<TObjArray*>(list->FindObject(fResPhiNegHistName));
    if (!fArrResoPhi_Neg) {
      LOGP(fatal, "Could not open {} from file {}", fResPhiNegHistName.Data(), fResFileName.Data());
    }
  }

  /*template <typename T>
  void getEfficiencyHistos(T dir){

  }*/

  void init()
  {
    if (fInitialized)
      return;

    if ((fResFileName.BeginsWith("alien://") || fEffFileName.BeginsWith("alien://")) && (!fFromCcdb)) {
      TGrid::Connect("alien://");
    }

    LOGP(info, "Set resolution histos");
    TList* listRes = new TList();
    if (fFromCcdb) {
      if (fCcdbPathRes.CompareTo("") != 0) {
        fResType = 1;
        listRes = fCcdb->getForTimeStamp<TList>(fCcdbPathRes.Data(), fTimestamp);
        if (!listRes) {
          LOGP(fatal, "Could not get resolution file from CCDB");
          return;
        }
      }
    } else {
      if (fResFileName.CompareTo("") != 0) {
        fResType = 1;
        TFile* fFile = TFile::Open(fResFileName);
        if (!fFile) {
          LOGP(fatal, "Could not open resolution file {}", fResFileName.Data());
          return;
        }
        if (fFile->GetListOfKeys()->Contains("ccdb_object")) {
          listRes = reinterpret_cast<TList*>(fFile->Get("ccdb_object"));
        } else {
          for (TObject* keyAsObj : *(fFile->GetListOfKeys())) {
            auto key = dynamic_cast<TKey*>(keyAsObj);
            TObject* arr = nullptr;
            fFile->GetObject(key->GetName(), arr);
            listRes->Add(arr);
          }
        }
        fFile->Close();
      }
    }
    if (fResType != 0) {
      fArrResoPt = reinterpret_cast<TObjArray*>(listRes->FindObject(fResPtHistName));
      if (!fArrResoPt) {
        LOGP(fatal, "Could not open {} from file {}", fResPtHistName.Data(), fResFileName.Data());
      }

      fArrResoEta = reinterpret_cast<TObjArray*>(listRes->FindObject(fResEtaHistName));
      if (!fArrResoEta) {
        LOGP(fatal, "Could not open {} from file {}", fResEtaHistName.Data(), fResFileName.Data());
      }

      fArrResoPhi_Pos = reinterpret_cast<TObjArray*>(listRes->FindObject(fResPhiPosHistName));
      if (!fArrResoPhi_Pos) {
        LOGP(fatal, "Could not open {} from file {}", fResPhiPosHistName.Data(), fResFileName.Data());
      }

      fArrResoPhi_Neg = reinterpret_cast<TObjArray*>(listRes->FindObject(fResPhiNegHistName));
      if (!fArrResoPhi_Neg) {
        LOGP(fatal, "Could not open {} from file {}", fResPhiNegHistName.Data(), fResFileName.Data());
      }
    }
    delete listRes;

    LOGP(info, "Set efficiency histos");
    TList* listEff = new TList();
    if (fFromCcdb) {
      if (fCcdbPathEff.CompareTo("") != 0) {
        fEffType = 1;
        listEff = fCcdb->getForTimeStamp<TList>(fCcdbPathEff.Data(), fTimestamp);
        if (!listEff) {
          LOGP(fatal, "Could not get efficiency file from CCDB");
          return;
        }
      }
    } else {
      if (fEffFileName.CompareTo("") != 0) {
        fEffType = 1;
        TFile* fFile = TFile::Open(fEffFileName);
        if (!fFile) {
          LOGP(fatal, "Could not open efficiency file {}", fEffFileName.Data());
          return;
        }
        if (fFile->GetListOfKeys()->Contains("ccdb_object")) {
          listEff = reinterpret_cast<TList*>(fFile->Get("ccdb_object"));
        } else {
          for (TObject* keyAsObj : *(fFile->GetListOfKeys())) {
            auto key = dynamic_cast<TKey*>(keyAsObj);
            TObject* hist = nullptr;
            fFile->GetObject(key->GetName(), hist);
            listEff->Add(hist);
          }
        }
        fFile->Close();
      }
    }
    if (fEffType != 0) {
      fArrEff = reinterpret_cast<TObject*>(listEff->FindObject(fEffHistName));
      if (!fArrEff) {
        LOGP(fatal, "Could not open {} from file {}", fEffHistName.Data(), fEffFileName.Data());
      }
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
    }
    delete listEff;

    fInitialized = true;
  }

  void applySmearing(const int ch, const float ptgen, const float etagen, const float phigen, float& ptsmeared, float& etasmeared, float& phismeared)
  {
    if (fResType == 0) {
      ptsmeared = ptgen;
      etasmeared = etagen;
      phismeared = phigen;
      return;
    }
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
  void setCcdbPathRes(TString ccdbPathRes) { fCcdbPathRes = ccdbPathRes; }
  void setCcdbPathEff(TString ccdbPathEff) { fCcdbPathEff = ccdbPathEff; }
  void setCcdb(Service<ccdb::BasicCCDBManager> ccdb)
  {
    fCcdb = ccdb;
    fFromCcdb = true;
  }
  void setTimestamp(int64_t timestamp) { fTimestamp = timestamp; }

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
  TString getCcdbPathRes() { return fCcdbPathRes; }
  TString getCcdbPathEff() { return fCcdbPathEff; }

 private:
  bool fInitialized = false;
  TString fResFileName;
  TString fResPtHistName;
  TString fResEtaHistName;
  TString fResPhiPosHistName;
  TString fResPhiNegHistName;
  TString fEffFileName;
  TString fEffHistName;
  TString fCcdbPathRes;
  TString fCcdbPathEff;
  int fEffType = 0;
  int fResType = 0;
  TObjArray* fArrResoPt;
  TObjArray* fArrResoEta;
  TObjArray* fArrResoPhi_Pos;
  TObjArray* fArrResoPhi_Neg;
  TObject* fArrEff;
  int64_t fTimestamp;
  bool fFromCcdb = false;
  Service<ccdb::BasicCCDBManager> fCcdb;
};

#endif // PWGEM_DILEPTON_UTILS_MOMENTUMSMEARER_H_
