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

#include "GFWWeights.h"
#include "TMath.h"
#include <cstdio>

GFWWeights::GFWWeights() : TNamed("", ""),
                           fDataFilled(kFALSE),
                           fMCFilled(kFALSE),
                           fW_data(0),
                           fW_mcrec(0),
                           fW_mcgen(0),
                           fEffInt(0),
                           fIntEff(0),
                           fAccInt(0),
                           fNbinsPt(0),
                           fbinsPt(0) {}
GFWWeights::GFWWeights(const char* name) : TNamed(name, name),
                                           fDataFilled(kFALSE),
                                           fMCFilled(kFALSE),
                                           fW_data(0),
                                           fW_mcrec(0),
                                           fW_mcgen(0),
                                           fEffInt(0),
                                           fIntEff(0),
                                           fAccInt(0),
                                           fNbinsPt(0),
                                           fbinsPt(0) {}
GFWWeights::~GFWWeights()
{
  delete fW_data;
  delete fW_mcrec;
  delete fW_mcgen;
  delete fEffInt;
  delete fIntEff;
  delete fAccInt;
  if (fbinsPt)
    delete[] fbinsPt;
};
void GFWWeights::setPtBins(int Nbins, double* bins)
{
  if (fbinsPt)
    delete[] fbinsPt;
  fNbinsPt = Nbins;
  fbinsPt = new double[fNbinsPt + 1];
  for (int i = 0; i <= fNbinsPt; ++i)
    fbinsPt[i] = bins[i];
};
void GFWWeights::init(bool AddData, bool AddMC)
{
  if (!fbinsPt) { // If pT bins not initialized, set to default (-1 to 1e6) to accept everything
    fNbinsPt = 1;
    fbinsPt = new double[2];
    fbinsPt[0] = -1;
    fbinsPt[1] = 1e6;
  }
  if (AddData) {
    fW_data = new TObjArray();
    fW_data->SetName("GFWWeights_Data");
    fW_data->SetOwner(kTRUE);
    const char* tnd = getBinName(0, 0, Form("data_%s", this->GetName()));
    fW_data->Add(new TH3D(tnd, ";#varphi;#eta;v_{z}", 60, 0, TMath::TwoPi(), 64, -1.6, 1.6, 40, -10, 10));
    fDataFilled = kTRUE;
  }
  if (AddMC) {
    fW_mcrec = new TObjArray();
    fW_mcrec->SetName("GFWWeights_MCRec");
    fW_mcgen = new TObjArray();
    fW_mcgen->SetName("GFWWeights_MCGen");
    fW_mcrec->SetOwner(kTRUE);
    fW_mcgen->SetOwner(kTRUE);
    const char* tnr = getBinName(0, 0, "mcrec"); // all integrated over cent. anyway
    const char* tng = getBinName(0, 0, "mcgen"); // all integrated over cent. anyway
    fW_mcrec->Add(new TH3D(tnr, ";#it{p}_{T};#eta;v_{z}", fNbinsPt, 0, 20, 64, -1.6, 1.6, 40, -10, 10));
    fW_mcgen->Add(new TH3D(tng, ";#it{p}_{T};#eta;v_{z}", fNbinsPt, 0, 20, 64, -1.6, 1.6, 40, -10, 10));
    reinterpret_cast<TH3D*>(fW_mcrec->At(fW_mcrec->GetEntries() - 1))->GetXaxis()->Set(fNbinsPt, fbinsPt);
    reinterpret_cast<TH3D*>(fW_mcgen->At(fW_mcgen->GetEntries() - 1))->GetXaxis()->Set(fNbinsPt, fbinsPt);
    fMCFilled = kTRUE;
  }
};

void GFWWeights::fill(double phi, double eta, double vz, double pt, double /*cent*/, int htype, double weight)
{
  TObjArray* tar = 0;
  const char* pf = "";
  if (htype == 0) {
    tar = fW_data;
    pf = Form("data_%s", this->GetName());
  }
  if (htype == 1) {
    tar = fW_mcrec;
    pf = "mcrec";
  }
  if (htype == 2) {
    tar = fW_mcgen;
    pf = "mcgen";
  }
  if (!tar)
    return;
  TH3D* th3 = reinterpret_cast<TH3D*>(tar->FindObject(getBinName(0, 0, pf))); // pT bin 0, V0M bin 0, since all integrated
  if (!th3) {
    if (!htype)
      tar->Add(new TH3D(getBinName(0, 0, pf), ";#varphi;#eta;v_{z}", 60, 0, TMath::TwoPi(), 64, -1.6, 1.6, 40, -10, 10)); // 0,0 since all integrated
    th3 = reinterpret_cast<TH3D*>(tar->At(tar->GetEntries() - 1));
  }
  th3->Fill(htype ? pt : phi, eta, vz, weight);
};
double GFWWeights::getWeight(double phi, double eta, double vz, double pt, double /*cent*/, int htype)
{
  TObjArray* tar = 0;
  const char* pf = "";
  if (htype == 0) {
    tar = fW_data;
    pf = "data";
  }
  if (htype == 1) {
    tar = fW_mcrec;
    pf = "mcrec";
  }
  if (htype == 2) {
    tar = fW_mcgen;
    pf = "mcgen";
  }
  if (!tar)
    return 1;
  TH3D* th3 = reinterpret_cast<TH3D*>(tar->FindObject(getBinName(0, 0, pf)));
  if (!th3)
    return 1; //-1;
  int xind = th3->GetXaxis()->FindBin(htype ? pt : phi);
  int etaind = th3->GetYaxis()->FindBin(eta);
  int vzind = th3->GetZaxis()->FindBin(vz);
  double weight = th3->GetBinContent(xind, etaind, vzind);
  if (weight != 0)
    return 1. / weight;
  return 1;
};
double GFWWeights::getNUA(double phi, double eta, double vz)
{
  if (!fAccInt)
    createNUA();
  int xind = fAccInt->GetXaxis()->FindBin(phi);
  int etaind = fAccInt->GetYaxis()->FindBin(eta);
  int vzind = fAccInt->GetZaxis()->FindBin(vz);
  double weight = fAccInt->GetBinContent(xind, etaind, vzind);
  if (weight != 0)
    return 1. / weight;
  return 1;
}
double GFWWeights::getNUE(double pt, double eta, double vz)
{
  if (!fEffInt)
    createNUE();
  int xind = fEffInt->GetXaxis()->FindBin(pt);
  int etaind = fEffInt->GetYaxis()->FindBin(eta);
  int vzind = fEffInt->GetZaxis()->FindBin(vz);
  double weight = fEffInt->GetBinContent(xind, etaind, vzind);
  if (weight != 0)
    return 1. / weight;
  return 1;
}
double GFWWeights::findMax(TH3D* inh, int& ix, int& iy, int& iz)
{
  double maxv = inh->GetBinContent(1, 1, 1);
  for (int i = 1; i <= inh->GetNbinsX(); i++)
    for (int j = 1; j <= inh->GetNbinsY(); j++)
      for (int k = 1; k <= inh->GetNbinsZ(); k++)
        if (inh->GetBinContent(i, j, k) > maxv) {
          ix = i;
          iy = j;
          iz = k;
          maxv = inh->GetBinContent(i, j, k);
        }
  return maxv;
};
void GFWWeights::mcToEfficiency()
{
  if (fW_mcgen->GetEntries() < 1) {
    LOGF(info, "MC gen. array empty. This is probably because effs. have been calculated and the generated particle histograms have been cleared out!\n");
    return;
  }
  for (int i = 0; i < fW_mcrec->GetEntries(); i++) {
    TH3D* hr = reinterpret_cast<TH3D*>(fW_mcrec->At(i));
    TH3D* hg = reinterpret_cast<TH3D*>(fW_mcgen->At(i));
    hr->Sumw2();
    hg->Sumw2();
    hr->Divide(hg);
  }
  fW_mcgen->Clear();
};
void GFWWeights::rebinNUA(int nX, int nY, int nZ)
{
  if (fW_data->GetEntries() < 1)
    return;
  for (int i = 0; i < fW_data->GetEntries(); i++) {
    reinterpret_cast<TH3D*>(fW_data->At(i))->RebinX(nX);
    reinterpret_cast<TH3D*>(fW_data->At(i))->RebinY(nY);
    reinterpret_cast<TH3D*>(fW_data->At(i))->RebinZ(nZ);
  }
};
void GFWWeights::createNUA(bool IntegrateOverCentAndPt)
{
  if (!IntegrateOverCentAndPt) {
    LOGF(info, "Method is outdated! NUA is integrated over centrality and pT. Quit now, or the behaviour will be bad\n");
    return;
  }
  TH1D* h1;
  if (fW_data->GetEntries() < 1)
    return;
  if (IntegrateOverCentAndPt) {
    if (fAccInt)
      delete fAccInt;
    fAccInt = reinterpret_cast<TH3D*>(fW_data->At(0)->Clone("IntegratedAcceptance"));
    fAccInt->Sumw2();
    for (int etai = 1; etai <= fAccInt->GetNbinsY(); etai++) {
      fAccInt->GetYaxis()->SetRange(etai, etai);
      if (fAccInt->Integral() < 1)
        continue;
      for (int vzi = 1; vzi <= fAccInt->GetNbinsZ(); vzi++) {
        fAccInt->GetZaxis()->SetRange(vzi, vzi);
        if (fAccInt->Integral() < 1)
          continue;
        h1 = reinterpret_cast<TH1D*>(fAccInt->Project3D("x"));
        double maxv = h1->GetMaximum();
        for (int phii = 1; phii <= h1->GetNbinsX(); phii++) {
          fAccInt->SetBinContent(phii, etai, vzi, fAccInt->GetBinContent(phii, etai, vzi) / maxv);
          fAccInt->SetBinError(phii, etai, vzi, fAccInt->GetBinError(phii, etai, vzi) / maxv);
        }
        delete h1;
      }
      fAccInt->GetZaxis()->SetRange(1, fAccInt->GetNbinsZ());
    }
    fAccInt->GetYaxis()->SetRange(1, fAccInt->GetNbinsY());
    return;
  }
};
TH1D* GFWWeights::getdNdPhi()
{
  TH3D* temph = reinterpret_cast<TH3D*>(fW_data->At(0)->Clone("tempH3"));
  TH1D* reth = reinterpret_cast<TH1D*>(temph->Project3D("x"));
  reth->SetName("RetHist");
  delete temph;
  double max = reth->GetMaximum();
  if (max == 0)
    return 0;
  for (int phi = 1; phi <= reth->GetNbinsX(); phi++) {
    if (reth->GetBinContent(phi) == 0)
      continue;
    reth->SetBinContent(phi, reth->GetBinContent(phi) / max);
    reth->SetBinError(phi, reth->GetBinError(phi) / max);
  }
  return reth;
}
void GFWWeights::createNUE(bool IntegrateOverCentrality)
{
  if (!IntegrateOverCentrality) {
    LOGF(info, "Method is outdated! NUE is integrated over centrality. Quit now, or the behaviour will be bad\n");
    return;
  }
  TH3D* num = 0;
  TH3D* den = 0;
  if (fW_mcrec->GetEntries() < 1 || fW_mcgen->GetEntries() < 1)
    return;
  if (IntegrateOverCentrality) {
    num = reinterpret_cast<TH3D*>(fW_mcrec->At(0));
    den = reinterpret_cast<TH3D*>(fW_mcgen->At(0));
    num->Sumw2();
    den->Sumw2();
    num->RebinY(2);
    den->RebinY(2);
    num->RebinZ(5);
    den->RebinZ(5);
    fEffInt = reinterpret_cast<TH3D*>(num->Clone("Efficiency_Integrated"));
    fEffInt->Divide(den);
    return;
  }
};
void GFWWeights::readAndMerge(TString filelinks, TString listName, bool addData, bool addRec, bool addGen)
{
  FILE* flist = fopen(filelinks.Data(), "r");
  char str[150];
  int nFiles = 0;
  while (fscanf(flist, "%s\n", str) == 1)
    nFiles++;
  rewind(flist);
  if (nFiles == 0) {
    LOGF(info, "No files to read!\n");
    return;
  }
  if (!fW_data && addData) {
    fW_data = new TObjArray();
    fW_data->SetName("Weights_Data");
    fW_data->SetOwner(kTRUE);
  }
  if (!fW_mcrec && addRec) {
    fW_mcrec = new TObjArray();
    fW_mcrec->SetName("Weights_MCRec");
    fW_mcrec->SetOwner(kTRUE);
  }
  if (!fW_mcgen && addGen) {
    fW_mcgen = new TObjArray();
    fW_mcgen->SetName("Weights_MCGen");
    fW_mcgen->SetOwner(kTRUE);
  }
  TFile* tf = 0;
  for (int i = 0; i < nFiles; i++) {
    auto retVal = fscanf(flist, "%s\n", str);
    (void)retVal;
    tf = new TFile(str, "READ");
    if (tf->IsZombie()) {
      LOGF(warning, "Could not open file %s!\n", str);
      tf->Close();
      continue;
    }
    TList* tl = reinterpret_cast<TList*>(tf->Get(listName.Data()));
    GFWWeights* tw = reinterpret_cast<GFWWeights*>(tl->FindObject(this->GetName()));
    if (!tw) {
      LOGF(warning, "Could not fetch weights object from %s\n", str);
      tf->Close();
      continue;
    }
    if (addData)
      addArray(fW_data, tw->getDataArray());
    if (addRec)
      addArray(fW_mcrec, tw->getRecArray());
    if (addGen)
      addArray(fW_mcgen, tw->getGenArray());
    tf->Close();
    delete tw;
  }
};
void GFWWeights::addArray(TObjArray* targ, TObjArray* sour)
{
  if (!sour) {
    LOGF(info, "Source array does not exist!\n");
    return;
  }
  for (int i = 0; i < sour->GetEntries(); i++) {
    TH3D* sourh = reinterpret_cast<TH3D*>(sour->At(i));
    TH3D* targh = reinterpret_cast<TH3D*>(targ->FindObject(sourh->GetName()));
    if (!targh) {
      targh = reinterpret_cast<TH3D*>(sourh->Clone(sourh->GetName()));
      targh->SetDirectory(0);
      targ->Add(targh);
    } else {
      targh->Add(sourh);
    }
  }
};
void GFWWeights::overwriteNUA()
{
  if (!fAccInt)
    createNUA();
  TString ts(fW_data->At(0)->GetName());
  TH3D* trash = reinterpret_cast<TH3D*>(fW_data->RemoveAt(0));
  delete trash;
  fW_data->Add(reinterpret_cast<TH3D*>(fAccInt->Clone(ts.Data())));
  delete fAccInt;
}
Long64_t GFWWeights::Merge(TCollection* collist)
{
  Long64_t nmerged = 0;
  if (!fW_data) {
    fW_data = new TObjArray();
    fW_data->SetName("Weights_Data");
    fW_data->SetOwner(kTRUE);
  }
  if (!fW_mcrec) {
    fW_mcrec = new TObjArray();
    fW_mcrec->SetName("Weights_MCRec");
    fW_mcrec->SetOwner(kTRUE);
  }
  if (!fW_mcgen) {
    fW_mcgen = new TObjArray();
    fW_mcgen->SetName("Weights_MCGen");
    fW_mcgen->SetOwner(kTRUE);
  }
  GFWWeights* l_w = 0;
  TIter all_w(collist);
  while ((l_w = (reinterpret_cast<GFWWeights*>(all_w())))) {
    addArray(fW_data, l_w->getDataArray());
    addArray(fW_mcrec, l_w->getRecArray());
    addArray(fW_mcgen, l_w->getGenArray());
    nmerged++;
  }
  return nmerged;
};
TH1D* GFWWeights::getIntegratedEfficiencyHist()
{
  if (!fW_mcgen) {
    LOGF(warning, "MCGen array does not exist!\n");
    return 0;
  }
  if (!fW_mcrec) {
    LOGF(warning, "MCRec array does not exist!\n");
    return 0;
  }
  if (!fW_mcgen->GetEntries()) {
    LOGF(warning, "MCGen array is empty!\n");
    return 0;
  }
  if (!fW_mcrec->GetEntries()) {
    LOGF(warning, "MCRec array is empty!\n");
    return 0;
  }
  TH3D* num = reinterpret_cast<TH3D*>(fW_mcrec->At(0)->Clone("Numerator"));
  for (int i = 1; i < fW_mcrec->GetEntries(); i++)
    num->Add(reinterpret_cast<TH3D*>(fW_mcrec->At(i)));
  TH3D* den = reinterpret_cast<TH3D*>(fW_mcgen->At(0)->Clone("Denominator"));
  for (int i = 1; i < fW_mcgen->GetEntries(); i++)
    den->Add(reinterpret_cast<TH3D*>(fW_mcgen->At(i)));
  TH1D* num1d = reinterpret_cast<TH1D*>(num->Project3D("x"));
  num1d->SetName("retHist");
  num1d->Sumw2();
  TH1D* den1d = reinterpret_cast<TH1D*>(den->Project3D("x"));
  den1d->Sumw2();
  num1d->Divide(den1d);
  delete num;
  delete den;
  delete den1d;
  return num1d;
}
bool GFWWeights::calculateIntegratedEff()
{
  if (fIntEff)
    delete fIntEff;
  fIntEff = getIntegratedEfficiencyHist();
  if (!fIntEff) {
    return kFALSE;
  }
  fIntEff->SetName("IntegratedEfficiency");
  return kTRUE;
}
double GFWWeights::getIntegratedEfficiency(double pt)
{
  if (!fIntEff)
    if (!calculateIntegratedEff())
      return 0;
  return fIntEff->GetBinContent(fIntEff->FindBin(pt));
}
TH1D* GFWWeights::getEfficiency(double etamin, double etamax, double vzmin, double vzmax)
{
  TH3D* num = reinterpret_cast<TH3D*>(fW_mcrec->At(0)->Clone("Numerator"));
  for (int i = 1; i < fW_mcrec->GetEntries(); i++)
    num->Add(reinterpret_cast<TH3D*>(fW_mcrec->At(i)));
  TH3D* den = reinterpret_cast<TH3D*>(fW_mcgen->At(0)->Clone("Denominator"));
  for (int i = 1; i < fW_mcgen->GetEntries(); i++)
    den->Add(reinterpret_cast<TH3D*>(fW_mcgen->At(i)));
  int eb1 = num->GetYaxis()->FindBin(etamin + 1e-6);
  int eb2 = num->GetYaxis()->FindBin(etamax - 1e-6);
  int vz1 = num->GetZaxis()->FindBin(vzmin + 1e-6);
  int vz2 = num->GetZaxis()->FindBin(vzmax - 1e-6);
  num->GetYaxis()->SetRange(eb1, eb2);
  num->GetZaxis()->SetRange(vz1, vz2);
  den->GetYaxis()->SetRange(eb1, eb2);
  den->GetZaxis()->SetRange(vz1, vz2);
  TH1D* num1d = reinterpret_cast<TH1D*>(num->Project3D("x"));
  TH1D* den1d = reinterpret_cast<TH1D*>(den->Project3D("x"));
  delete num;
  delete den;
  num1d->Sumw2();
  den1d->Sumw2();
  num1d->Divide(den1d);
  delete den1d;
  return num1d;
}
void GFWWeights::mergeWeights(GFWWeights* other)
{
  if (!fW_data) {
    fW_data = new TObjArray();
    fW_data->SetName("Weights_Data");
    fW_data->SetOwner(kTRUE);
  }
  addArray(fW_data, other->getDataArray());
  return;
}
void GFWWeights::setTH3D(TH3D* th3d)
{
  if (!fW_data) {
    fW_data = new TObjArray();
    fW_data->SetName("GFWWeights_Data");
    fW_data->SetOwner(kTRUE);
    fW_data->Add(th3d);
    return;
  }
  TString ts(fW_data->At(0)->GetName());
  TH3D* trash = reinterpret_cast<TH3D*>(fW_data->RemoveAt(0));
  delete trash;
  fW_data->Add(reinterpret_cast<TH3D*>(th3d->Clone(ts.Data())));
}
