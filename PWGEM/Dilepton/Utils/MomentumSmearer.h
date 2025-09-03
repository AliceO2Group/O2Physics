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

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"

#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TString.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

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

  /// Constructor with resolution ND sparse histogram
  MomentumSmearer(TString resFileName, TString resNDHistName)
  {
    setResFileName(resFileName);
    setResNDHistName(resNDHistName);
    setResPtHistName("");
    setResEtaHistName("");
    setResPhiPosHistName("");
    setResPhiNegHistName("");
    setEffFileName("");
    setEffHistName("");
    setDCAFileName("");
    setDCAHistName("");
    fDoNDSmearing = true;
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

  void fillVecReso(TH2F* fReso, std::vector<TH1F*>& fVecReso, const char* suffix)
  {
    TAxis* axisPt = fReso->GetXaxis(); // be careful! This works only for variable bin width.
    int nBinsPt = axisPt->GetNbins();
    fVecReso.resize(nBinsPt);
    for (int i = 0; i < nBinsPt; i++) {
      auto h1 = reinterpret_cast<TH1F*>(fReso->ProjectionY(Form("h1reso%s_pt%d", suffix, i), i + 1, i + 1));
      h1->Scale(1.f, "width"); // convert ntrack to probability density
      fVecReso[i] = h1;
    }
  }

  void fillVecResoND(THnSparseF* hs_reso)
  {
    LOGP(info, "prepare TH3D");
    fNCenBins = hs_reso->GetAxis(0)->GetNbins();
    fNPtBins = hs_reso->GetAxis(1)->GetNbins();
    fNEtaBins = hs_reso->GetAxis(2)->GetNbins();
    fNPhiBins = hs_reso->GetAxis(3)->GetNbins();
    fNChBins = hs_reso->GetAxis(4)->GetNbins();
    LOGF(info, "ncen = %d, npt = %d, neta = %d, nphi = %d, nch = %d without under- and overflow bins", fNCenBins, fNPtBins, fNEtaBins, fNPhiBins, fNChBins);
    // fVecResoND.reserve(npt * neta * nphi * nch);

    fVecResoND.resize(fNCenBins, std::vector<std::vector<std::vector<std::vector<TH3D*>>>>(fNPtBins, std::vector<std::vector<std::vector<TH3D*>>>(fNEtaBins, std::vector<std::vector<TH3D*>>(fNPhiBins, std::vector<TH3D*>(fNChBins)))));
    // fVecResoND.resize(fNPtBins, std::vector<std::vector<std::vector<TH3D*>>>(fNEtaBins, std::vector<std::vector<TH3D*>>(fNPhiBins, std::vector<TH3D*>(fNChBins))));
    //  auto h3 = reinterpret_cast<TH3D*>(hs_reso->Projection(4, 5, 6));
    //  h3->SetName(Form("h3reso_pt%d_eta%d_phi%d_ch%d", 0, 0, 0, 0));
    //  fVecResoND[0][0][0][0] = h3;

    for (int icen = 0; icen < fNCenBins; icen++) {
      hs_reso->GetAxis(0)->SetRange(icen + 1, icen + 1);
      for (int ipt = 0; ipt < fNPtBins; ipt++) {
        hs_reso->GetAxis(1)->SetRange(ipt + 1, ipt + 1);
        for (int ieta = 0; ieta < fNEtaBins; ieta++) {
          hs_reso->GetAxis(2)->SetRange(ieta + 1, ieta + 1);
          for (int iphi = 0; iphi < fNPhiBins; iphi++) {
            hs_reso->GetAxis(3)->SetRange(iphi + 1, iphi + 1);
            for (int ich = 0; ich < fNChBins; ich++) {
              if (-0.5 < hs_reso->GetAxis(4)->GetBinCenter(ich + 1) && hs_reso->GetAxis(4)->GetBinCenter(ich + 1) < 0.5) {
                continue;
              }
              hs_reso->GetAxis(4)->SetRange(ich + 1, ich + 1);
              auto h3 = reinterpret_cast<TH3D*>(hs_reso->Projection(5, 6, 7));
              h3->SetName(Form("h3reso_cen%d_pt%d_eta%d_phi%d_ch%d", icen, ipt, ieta, iphi, ich));
              h3->Scale(1.f, "width"); // convert ntrack to probability density
              fVecResoND[icen][ipt][ieta][iphi][ich] = h3;
            } // end of charge loop
          } // end of phi loop
        } // end of eta loop
      } // end of pt loop
    } // end of centrality loop
  }

  void init()
  {
    if (fInitialized) {
      return;
    }

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

    LOGF(info, "Apply ND-correlated smearing: fDoNDSmearing = %d , fResNDHistName = %s", fDoNDSmearing, fResNDHistName.Data());

    if (fDoNDSmearing) {
      if (fResType != 0) {
        fResoND = reinterpret_cast<THnSparseF*>(listRes->FindObject(fResNDHistName));
        if (!fResoND) {
          LOGP(fatal, "Could not open {} from file {}", fResNDHistName.Data(), fResFileName.Data());
        }
        fillVecResoND(fResoND);
      }
    } else {
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
        fillVecReso(fResoPt, fVecResoPt, "_reldpt");
        fillVecReso(fResoEta, fVecResoEta, "_deta");
        fillVecReso(fResoPhi_Pos, fVecResoPhi_Pos, "_dphi_pos");
        fillVecReso(fResoPhi_Neg, fVecResoPhi_Neg, "_dphi_neg");
      }
    }

    if (!fFromCcdb) {
      delete listRes;
    }

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

    if (!fFromCcdb) {
      delete listEff;
    }

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
      fillVecReso(fDCA, fVecDCA, "_dca");
    }

    if (!fFromCcdb) {
      delete listDCA;
    }

    fInitialized = true;
  }

  void applySmearing(const float ptgen, const float vargen, const float multiply, float& varsmeared, TH2F* fReso, std::vector<TH1F*>& fVecReso)
  {
    float ptgen_tmp = ptgen > fMinPtGen ? ptgen : fMinPtGen;
    TAxis* axisPt = fReso->GetXaxis();
    int nBinsPt = axisPt->GetNbins();
    int ptbin = axisPt->FindBin(ptgen_tmp);
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

  void applySmearing(const float centrality, const int ch, const float ptgen, const float etagen, const float phigen, float& ptsmeared, float& etasmeared, float& phismeared)
  {
    if (fResType == 0) {
      ptsmeared = ptgen;
      etasmeared = etagen;
      phismeared = phigen;
      return;
    }

    if (fDoNDSmearing) {
      if (centrality < 0) {
        ptsmeared = ptgen;
        etasmeared = etagen;
        phismeared = phigen;
        return;
      }
      applySmearingND(centrality, ch, ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
    } else {
      applySmearing(ptgen, ptgen, ptgen, ptsmeared, fResoPt, fVecResoPt);
      applySmearing(ptgen, etagen, 1., etasmeared, fResoEta, fVecResoEta);
      if (ch > 0) {
        applySmearing(ptgen, phigen, 1., phismeared, fResoPhi_Pos, fVecResoPhi_Pos);
      } else {
        applySmearing(ptgen, phigen, 1., phismeared, fResoPhi_Neg, fVecResoPhi_Neg);
      }
    }
  }

  void applySmearingND(const float centrality, const int ch, const float ptgen, const float etagen, const float phigen, float& ptsmeared, float& etasmeared, float& phismeared)
  {
    float ptgen_tmp = ptgen > fMinPtGen ? ptgen : fMinPtGen;
    int cenbin = fResoND->GetAxis(0)->FindBin(centrality);
    int ptbin = fResoND->GetAxis(1)->FindBin(ptgen_tmp);
    int etabin = fResoND->GetAxis(2)->FindBin(etagen);
    int phibin = fResoND->GetAxis(3)->FindBin(phigen);
    int chbin = fResoND->GetAxis(4)->FindBin(ch);

    // protection
    if (cenbin < 1) {
      cenbin = 1;
    } else if (cenbin > fNCenBins) {
      cenbin = fNCenBins;
    }

    // protection
    if (ptbin < 1) {
      ptbin = 1;
    } else if (ptbin > fNPtBins) {
      ptbin = fNPtBins;
    }

    // protection
    if (etabin < 1) {
      etabin = 1;
    } else if (etabin > fNEtaBins) {
      etabin = fNEtaBins;
    }

    // protection
    if (phibin < 1) {
      phibin = 1;
    } else if (phibin > fNPhiBins) {
      phibin = fNPhiBins;
    }

    // protection
    if (chbin < 1) {
      chbin = 1;
    } else if (chbin > fNChBins) {
      chbin = fNChBins;
    }

    double dpt_rel = 0, deta = 0, dphi = 0;
    if (fVecResoND[cenbin - 1][ptbin - 1][etabin - 1][phibin - 1][chbin - 1]->GetEntries() > 0) {
      fVecResoND[cenbin - 1][ptbin - 1][etabin - 1][phibin - 1][chbin - 1]->GetRandom3(dpt_rel, deta, dphi);
    }
    ptsmeared = ptgen - dpt_rel * ptgen;
    etasmeared = etagen - deta;
    phismeared = phigen - dphi;
    // LOGF(info, "ptgen = %f (GeV/c), etagen = %f, phigen = %f (rad.), ptsmeared = %f (GeV/c), etasmeared = %f, phismeared = %f (rad.)", ptgen, etagen, phigen, ptsmeared, etasmeared, phismeared);
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
  void setNDSmearing(bool flag) { fDoNDSmearing = flag; }
  void setResFileName(TString resFileName) { fResFileName = resFileName; }
  void setResNDHistName(TString resNDHistName) { fResNDHistName = resNDHistName; }
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
  void setMinPt(float minpt) { fMinPtGen = minpt; }

  // getters
  bool getNDSmearing() { return fDoNDSmearing; }
  TString getResFileName() { return fResFileName; }
  TString getResNDHistName() { return fResNDHistName; }
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
  float getMinPt() { return fMinPtGen; }

 private:
  bool fInitialized = false;
  bool fDoNDSmearing = false;
  TString fResFileName;
  TString fResNDHistName;
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
  THnSparseF* fResoND;
  TH2F* fResoPt;
  TH2F* fResoEta;
  TH2F* fResoPhi_Pos;
  TH2F* fResoPhi_Neg;
  std::vector<std::vector<std::vector<std::vector<std::vector<TH3D*>>>>> fVecResoND;
  int fNCenBins = 1;
  int fNPtBins = 1;
  int fNEtaBins = 1;
  int fNPhiBins = 1;
  int fNChBins = 1;
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
  float fMinPtGen = -1.f;
};

#endif // PWGEM_DILEPTON_UTILS_MOMENTUMSMEARER_H_
