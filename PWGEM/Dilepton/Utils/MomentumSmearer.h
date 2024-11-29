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
#include <TFile.h>
#include <TKey.h>
#include <CCDB/BasicCCDBManager.h>
#include <vector>
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
    setEffFileName("");
    setEffHistName("");
    setDCAFileName("");
    setDCAHistName("");
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
    setEffFileName(effFileName);
    setEffHistName(effHistName);
    setDCAFileName("");
    setDCAHistName("");
    init();
  }

  /// Constructor with resolution histograms and efficiency and dca
  MomentumSmearer(TString resFileName, TString resPtHistName, TString resEtaHistName, TString resPhiPosHistName, TString resPhiNegHistName, TString effFileName, TString effHistName, TString dcaFileName, TString dcaHistName)
  {
    setResFileName(resFileName);
    setResPtHistName(resPtHistName);
    setResEtaHistName(resEtaHistName);
    setResPhiPosHistName(resPhiPosHistName);
    setResPhiNegHistName(resPhiNegHistName);
    setEffFileName(effFileName);
    setEffHistName(effHistName);
    setDCAFileName(dcaFileName);
    setDCAHistName(dcaHistName);
    init();
  }

  /// Default destructor
  ~MomentumSmearer() = default;

  // template <typename T>
  void getResolutionHistos(TList* list)
  {

    fResoPt = reinterpret_cast<TH2F*>(list->FindObject(fResPtHistName));
    if (!fResoPt) {
      LOGP(fatal, "Could not open {} from file {}", fResPtHistName.Data(), fResFileName.Data());
    }

    fResoEta = reinterpret_cast<TH2F*>(list->FindObject(fResEtaHistName));
    if (!fResoEta) {
      LOGP(fatal, "Could not open {} from file {}", fResEtaHistName.Data(), fResFileName.Data());
    }

    fResoPhi_Pos = reinterpret_cast<TH2F*>(list->FindObject(fResPhiPosHistName));
    if (!fResoPhi_Pos) {
      LOGP(fatal, "Could not open {} from file {}", fResPhiPosHistName.Data(), fResFileName.Data());
    }

    fResoPhi_Neg = reinterpret_cast<TH2F*>(list->FindObject(fResPhiNegHistName));
    if (!fResoPhi_Neg) {
      LOGP(fatal, "Could not open {} from file {}", fResPhiNegHistName.Data(), fResFileName.Data());
    }
  }

  void fillVecReso(TH2F* fReso, std::vector<TH1F*>& fVecReso)
  {
    TAxis* axisPt = fReso->GetXaxis();
    int nBinsPt = axisPt->GetNbins();
    for (int i = 1; i <= nBinsPt; i++) {
      fVecReso.push_back(reinterpret_cast<TH1F*>(fReso->ProjectionY("", i, i)));
    }
  }

  void init()
  {
    if (fInitialized)
      return;

    if ((fResFileName.BeginsWith("alien://") || fEffFileName.BeginsWith("alien://") || fDCAFileName.BeginsWith("alien://")) && (!fFromCcdb)) {
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
      fResoPt = reinterpret_cast<TH2F*>(listRes->FindObject(fResPtHistName));
      if (!fResoPt) {
        LOGP(fatal, "Could not open {} from file {}", fResPtHistName.Data(), fResFileName.Data());
      }

      fResoEta = reinterpret_cast<TH2F*>(listRes->FindObject(fResEtaHistName));
      if (!fResoEta) {
        LOGP(fatal, "Could not open {} from file {}", fResEtaHistName.Data(), fResFileName.Data());
      }

      fResoPhi_Pos = reinterpret_cast<TH2F*>(listRes->FindObject(fResPhiPosHistName));
      if (!fResoPhi_Pos) {
        LOGP(fatal, "Could not open {} from file {}", fResPhiPosHistName.Data(), fResFileName.Data());
      }

      fResoPhi_Neg = reinterpret_cast<TH2F*>(listRes->FindObject(fResPhiNegHistName));
      if (!fResoPhi_Neg) {
        LOGP(fatal, "Could not open {} from file {}", fResPhiNegHistName.Data(), fResFileName.Data());
      }
      fillVecReso(fResoPt, fVecResoPt);
      fillVecReso(fResoEta, fVecResoEta);
      fillVecReso(fResoPhi_Pos, fVecResoPhi_Pos);
      fillVecReso(fResoPhi_Neg, fVecResoPhi_Neg);
    }

    if (!fFromCcdb)
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
      fEff = reinterpret_cast<TObject*>(listEff->FindObject(fEffHistName));
      if (!fEff) {
        LOGP(fatal, "Could not open {} from file {}", fEffHistName.Data(), fEffFileName.Data());
      }
      // check which type is used
      if (dynamic_cast<TH3*>(fEff)) {
        fEffType = 3;
        LOGP(info, "Use 3d efficiency histo (pt, eta, phi)");
      } else if (dynamic_cast<TH2*>(fEff)) {
        fEffType = 2;
        LOGP(info, "Use 2d efficiency histo (pt, eta)");
      } else if (dynamic_cast<TH1*>(fEff)) {
        fEffType = 1;
        LOGP(info, "Use 1d efficiency histo (pt)");
      } else {
        LOGP(fatal, "Could not identify type of histogram {}", fEffHistName.Data());
      }
    }

    if (!fFromCcdb)
      delete listEff;

    LOGP(info, "Set DCA histos");
    TList* listDCA = new TList();
    if (fFromCcdb) {
      if (fCcdbPathDCA.CompareTo("") != 0) {
        fDCAType = 1;
        listDCA = fCcdb->getForTimeStamp<TList>(fCcdbPathDCA.Data(), fTimestamp);
        if (!listDCA) {
          LOGP(fatal, "Could not get DCA file from CCDB");
          return;
        }
      }
    } else {
      if (fDCAFileName.CompareTo("") != 0) {
        fDCAType = 1;
        TFile* fFile = TFile::Open(fDCAFileName);
        if (!fFile) {
          LOGP(fatal, "Could not open DCA file {}", fDCAFileName.Data());
          return;
        }
        if (fFile->GetListOfKeys()->Contains("ccdb_object")) {
          listDCA = reinterpret_cast<TList*>(fFile->Get("ccdb_object"));
        } else {
          for (TObject* keyAsObj : *(fFile->GetListOfKeys())) {
            auto key = dynamic_cast<TKey*>(keyAsObj);
            TObject* arr = nullptr;
            fFile->GetObject(key->GetName(), arr);
            listDCA->Add(arr);
          }
        }
        fFile->Close();
      }
    }
    if (fDCAType != 0) {
      fDCA = reinterpret_cast<TH2F*>(listDCA->FindObject(fDCAHistName));
      if (!fDCA) {
        LOGP(fatal, "Could not open {} from file {}", fDCAHistName.Data(), fDCAFileName.Data());
      }
      fillVecReso(fDCA, fVecDCA);
    }

    if (!fFromCcdb)
      delete listDCA;

    fInitialized = true;
  }

  void applySmearing(float ptgen, float vargen, float multiply, float& varsmeared, TH2F* fReso, std::vector<TH1F*>& fVecReso)
  {
    TAxis* axisPt = fReso->GetXaxis();
    int nBinsPt = axisPt->GetNbins();
    int ptbin = axisPt->FindBin(ptgen);
    if (ptbin < 1) {
      ptbin = 1;
    }
    if (ptbin > nBinsPt) {
      ptbin = nBinsPt;
    }
    float smearing = 0.;
    if (fVecReso[ptbin - 1]->GetEntries() > 0) {
      smearing = fVecReso[ptbin - 1]->GetRandom() * multiply;
    }
    varsmeared = vargen - smearing;
  }

  void applySmearing(const int ch, const float ptgen, const float etagen, const float phigen, float& ptsmeared, float& etasmeared, float& phismeared)
  {
    if (fResType == 0) {
      ptsmeared = ptgen;
      etasmeared = etagen;
      phismeared = phigen;
      return;
    }
    applySmearing(ptgen, ptgen, ptgen, ptsmeared, fResoPt, fVecResoPt);
    applySmearing(ptgen, etagen, 1., etasmeared, fResoEta, fVecResoEta);
    if (ch > 0) {
      applySmearing(ptgen, phigen, 1., phismeared, fResoPhi_Pos, fVecResoPhi_Pos);
    } else {
      applySmearing(ptgen, phigen, 1., phismeared, fResoPhi_Neg, fVecResoPhi_Neg);
    }
  }

  float getEfficiency(float pt, float eta, float phi)
  {

    if (fEffType == 0) {
      return 1.;
    }

    if (fEffType == 1) {
      TH1F* hist = reinterpret_cast<TH1F*>(fEff);
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
      TH2F* hist = reinterpret_cast<TH2F*>(fEff);
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
      TH3F* hist = reinterpret_cast<TH3F*>(fEff);
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

  float getDCA(float ptsmeared)
  {
    if (fDCAType == 0) {
      return 0.;
    }

    TAxis* axisPt = fDCA->GetXaxis();
    int nBinsPt = axisPt->GetNbins();
    int ptbin = axisPt->FindBin(ptsmeared);
    if (ptbin < 1) {
      ptbin = 1;
    }
    if (ptbin > nBinsPt) {
      ptbin = nBinsPt;
    }
    float dca = 0.;
    if (fVecDCA[ptbin - 1]->GetEntries() > 0) {
      dca = fVecDCA[ptbin - 1]->GetRandom();
    }
    return dca;
  }

  // setters
  void setResFileName(TString resFileName) { fResFileName = resFileName; }
  void setResPtHistName(TString resPtHistName) { fResPtHistName = resPtHistName; }
  void setResEtaHistName(TString resEtaHistName) { fResEtaHistName = resEtaHistName; }
  void setResPhiPosHistName(TString resPhiPosHistName) { fResPhiPosHistName = resPhiPosHistName; }
  void setResPhiNegHistName(TString resPhiNegHistName) { fResPhiNegHistName = resPhiNegHistName; }
  void setEffFileName(TString effFileName) { fEffFileName = effFileName; }
  void setEffHistName(TString effHistName) { fEffHistName = effHistName; }
  void setDCAFileName(TString dcaFileName) { fDCAFileName = dcaFileName; }
  void setDCAHistName(TString dcaHistName) { fDCAHistName = dcaHistName; }
  void setCcdbPathRes(TString ccdbPathRes) { fCcdbPathRes = ccdbPathRes; }
  void setCcdbPathEff(TString ccdbPathEff) { fCcdbPathEff = ccdbPathEff; }
  void setCcdbPathDCA(TString ccdbPathDCA) { fCcdbPathDCA = ccdbPathDCA; }
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
  TString getDCAFileName() { return fDCAFileName; }
  TString getDCaHistName() { return fDCAHistName; }
  TH2F* getHistResoPt() { return fResoPt; }
  TH2F* getHistResoEta() { return fResoEta; }
  TH2F* getHistResoPhiPos() { return fResoPhi_Pos; }
  TH2F* getHistResoPhiNeg() { return fResoPhi_Neg; }
  TObject* getHistEff() { return fEff; }
  TString getCcdbPathRes() { return fCcdbPathRes; }
  TString getCcdbPathEff() { return fCcdbPathEff; }
  TString getCcdbPathDCA() { return fCcdbPathDCA; }

 private:
  bool fInitialized = false;
  TString fResFileName;
  TString fResPtHistName;
  TString fResEtaHistName;
  TString fResPhiPosHistName;
  TString fResPhiNegHistName;
  TString fEffFileName;
  TString fEffHistName;
  TString fDCAFileName;
  TString fDCAHistName;
  TString fCcdbPathRes;
  TString fCcdbPathEff;
  TString fCcdbPathDCA;
  int fEffType = 0;
  int fResType = 0;
  int fDCAType = 0;
  TH2F* fResoPt;
  TH2F* fResoEta;
  TH2F* fResoPhi_Pos;
  TH2F* fResoPhi_Neg;
  std::vector<TH1F*> fVecResoPt;
  std::vector<TH1F*> fVecResoEta;
  std::vector<TH1F*> fVecResoPhi_Pos;
  std::vector<TH1F*> fVecResoPhi_Neg;
  TObject* fEff;
  TH2F* fDCA;
  std::vector<TH1F*> fVecDCA;
  int64_t fTimestamp;
  bool fFromCcdb = false;
  Service<ccdb::BasicCCDBManager> fCcdb;
};

#endif // PWGEM_DILEPTON_UTILS_MOMENTUMSMEARER_H_
