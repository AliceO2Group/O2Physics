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

#include "FlowContainer.h"

#include <cstdio>
#include <vector>

ClassImp(FlowContainer);

FlowContainer::FlowContainer() : TNamed("", ""),
                                 fProf(0),
                                 fProfRand(0),
                                 fNRandom(0),
                                 fIDName("MidV"),
                                 fPtRebin(1),
                                 fPtRebinEdges(0),
                                 fMultiRebin(0),
                                 fMultiRebinEdges(0),
                                 fXAxis(0),
                                 fNbinsPt(0),
                                 fbinsPt(0),
                                 fPropagateErrors(kFALSE) {}
FlowContainer::FlowContainer(const char* name) : TNamed(name, name),
                                                 fProf(0),
                                                 fProfRand(0),
                                                 fNRandom(0),
                                                 fIDName("MidV"),
                                                 fPtRebin(1),
                                                 fPtRebinEdges(0),
                                                 fMultiRebin(0),
                                                 fMultiRebinEdges(0),
                                                 fXAxis(0),
                                                 fNbinsPt(0),
                                                 fbinsPt(0),
                                                 fPropagateErrors(kFALSE) {}
FlowContainer::~FlowContainer()
{
  delete fProf;
  delete fProfRand;
};
void FlowContainer::Initialize(TObjArray* inputList, const o2::framework::AxisSpec axis, int nRandom)
{
  std::vector<double> multiBins = axis.binEdges;
  int nMultiBins = axis.nBins.value_or(0);
  if (nMultiBins <= 0)
    nMultiBins = multiBins.size() - 1;
  if (nMultiBins <= 0) {
    printf("Multiplicity axis does not exist");
    return;
  }
  if (!inputList) {
    printf("Input list not specified\n");
    return;
  }
  if (inputList->GetEntries() < 1) {
    printf("Input list empty!\n");
    return;
  }
  fProf = new TProfile2D(Form("%s_CorrProfile", this->GetName()), "CorrProfile", nMultiBins, &multiBins[0], inputList->GetEntries(), 0.5, inputList->GetEntries() + 0.5);
  for (int i = 0; i < inputList->GetEntries(); i++)
    fProf->GetYaxis()->SetBinLabel(i + 1, inputList->At(i)->GetName());
  fProf->Sumw2();
  if (nRandom) {
    fNRandom = nRandom;
    fProfRand = new TObjArray();
    fProfRand->SetOwner(kTRUE);
    for (int i = 0; i < nRandom; i++) {
      fProfRand->Add(dynamic_cast<TProfile2D*>(fProf->Clone(Form("%s_Rand_%i", fProf->GetName(), i))));
      dynamic_cast<TProfile2D*>(fProfRand->At(i))->Sumw2();
    }
  }
};
void FlowContainer::Initialize(TObjArray* inputList, int nMultiBins, double MultiMin, double MultiMax, int nRandom)
{
  if (!inputList) {
    printf("Input list not specified\n");
    return;
  }
  if (inputList->GetEntries() < 1) {
    printf("Input list empty!\n");
    return;
  }
  fProf = new TProfile2D(Form("%s_CorrProfile", this->GetName()), "CorrProfile", nMultiBins, MultiMin, MultiMax, inputList->GetEntries(), 0.5, inputList->GetEntries() + 0.5);
  fProf->SetDirectory(0);
  fProf->Sumw2();
  for (int i = 0; i < inputList->GetEntries(); i++)
    fProf->GetYaxis()->SetBinLabel(i + 1, inputList->At(i)->GetName());
  if (nRandom) {
    fNRandom = nRandom;
    fProfRand = new TObjArray();
    fProfRand->SetOwner(kTRUE);
    for (int i = 0; i < nRandom; i++) {
      fProfRand->Add(dynamic_cast<TProfile2D*>(fProf->Clone(Form("%s_Rand_%i", fProf->GetName(), i))));
      dynamic_cast<TProfile2D*>(fProfRand->At(i))->Sumw2();
    }
  }
};
bool FlowContainer::CreateBinsFromAxis(TAxis* inax)
{
  if (!inax)
    return kFALSE;
  fNbinsPt = inax->GetNbins();
  fbinsPt = new double[fNbinsPt + 1];
  inax->GetLowEdge(fbinsPt);
  fbinsPt[fNbinsPt] = inax->GetBinUpEdge(fNbinsPt);
  return kTRUE;
}
void FlowContainer::SetXAxis(TAxis* inax)
{
  fXAxis = dynamic_cast<TAxis*>(inax->Clone("pTAxis"));
  bool success = CreateBinsFromAxis(fXAxis);
  if (!success)
    printf("Something went wrong setting the x axis!\n");
}
void FlowContainer::SetXAxis()
{
  if (!CreateBinsFromAxis(fXAxis)) { // Legacy; if fXAxis not defined, then setup default one
    const int NbinsPtForV2 = 24;
    double binsPtForV2[NbinsPtForV2 + 1] = {
      0.2, 0.4, 0.6, 0.8, 1.0,
      1.2, 1.4, 1.6, 1.8, 2.0,
      2.2, 2.4, 2.6, 3.0, 3.4,
      3.8, 4.2, 4.6, 5.2, 5.8,
      6.6, 8.0, 12.0, 16.0, 20.0};
    TAxis* tempax = new TAxis(NbinsPtForV2, binsPtForV2);
    SetXAxis(tempax);
    delete tempax;
  }
}
int FlowContainer::FillProfile(const char* hname, double multi, double corr, double w, double rn)
{
  if (!fProf)
    return -1;
  int yin = fProf->GetYaxis()->FindBin(hname);
  if (!yin) {
    printf("Could not find bin %s\n", hname);
    return -1;
  }
  fProf->Fill(multi, yin, corr, w);
  if (fNRandom) {
    double rnind = rn * fNRandom;
    dynamic_cast<TProfile2D*>(fProfRand->At(static_cast<int>(rnind)))->Fill(multi, yin, corr, w);
  }
  return 0;
};
void FlowContainer::OverrideProfileErrors(TProfile2D* inpf)
{
  int nBinsX = fProf->GetNbinsX();
  int nBinsY = fProf->GetNbinsY();
  if ((inpf->GetNbinsX() != nBinsX) || (inpf->GetNbinsY() != nBinsY)) {
    printf("Number of bins in two profiles do not match, not doing anything\n");
    return;
  }
  if (!inpf->GetBinSumw2()->fArray) {
    printf("Input profile has no BinSumw2()! Returning\n");
    return;
  }
  if (!fProf->GetBinSumw2()->fArray)
    fProf->Sumw2();
  double* sumw2Prof = fProf->GetSumw2()->fArray;
  double* sumw2Targ = inpf->GetSumw2()->fArray;
  double* binsw2Prof = fProf->GetBinSumw2()->fArray;
  double* binsw2Targ = inpf->GetBinSumw2()->fArray;
  double* farrProf = fProf->fArray;
  for (int ix = 1; ix <= nBinsX; ix++) {
    double xval = fProf->GetXaxis()->GetBinCenter(ix);
    printf("Processing x-bin %i\n", ix);
    for (int iy = 1; iy <= nBinsY; iy++) {
      double yval = fProf->GetYaxis()->GetBinCenter(iy);
      int binno = fProf->FindBin(xval, yval);
      double h = fProf->GetBinContent(binno);
      double lEnt = inpf->GetBinEntries(binno);
      fProf->SetBinEntries(binno, lEnt);
      sumw2Prof[binno] = sumw2Targ[binno];
      binsw2Prof[binno] = binsw2Targ[binno];
      farrProf[binno] = h * lEnt;
    }
  }
}

Long64_t FlowContainer::Merge(TCollection* collist)
{
  Long64_t nmerged = 0;
  FlowContainer* l_FC = 0;
  TIter all_FC(collist);
  while ((l_FC = dynamic_cast<FlowContainer*>(all_FC()))) {
    if (!fProf)
      continue;
    TProfile2D* tpro = GetProfile();
    TProfile2D* spro = l_FC->GetProfile();
    if (!tpro) {
      fProf = dynamic_cast<TProfile2D*>(spro->Clone(spro->GetName()));
      fProf->SetDirectory(0);
    } else {
      tpro->Add(spro);
    }
    nmerged++;
    TObjArray* tarr = l_FC->GetSubProfiles();
    if (!tarr)
      continue;
    if (!fProfRand) {
      fProfRand = new TObjArray();
      fProfRand->SetOwner(kTRUE);
    }
    for (int i = 0; i < tarr->GetEntries(); i++) {
      if (!(fProfRand->FindObject(tarr->At(i)->GetName()))) {
        fProfRand->Add(dynamic_cast<TProfile2D*>(tarr->At(i)->Clone(tarr->At(i)->GetName())));
        dynamic_cast<TProfile2D*>(fProfRand->At(fProfRand->GetEntries() - 1))->SetDirectory(0);
      } else {
        dynamic_cast<TProfile2D*>(fProfRand->FindObject(tarr->At(i)->GetName()))->Add(dynamic_cast<TProfile2D*>(tarr->At(i)));
      }
    }
  }
  return nmerged;
}

void FlowContainer::ReadAndMerge(const char* filelist)
{
  FILE* flist = fopen(filelist, "r");
  char str[150];
  int nFiles = 0;
  while (fscanf(flist, "%s\n", str) == 1)
    nFiles++;
  rewind(flist);
  if (nFiles == 0) {
    printf("No files to read!\n");
    return;
  }
  for (int i = 0; i < nFiles; i++) {
    auto retVal = fscanf(flist, "%s\n", str);
    (void)retVal;
    TFile* tf = new TFile(str, "READ");
    if (tf->IsZombie()) {
      printf("Could not open file %s!\n", str);
      tf->Close();
      continue;
    }
    PickAndMerge(tf);
    tf->Close();
  }
}
void FlowContainer::PickAndMerge(TFile* tfi)
{
  FlowContainer* lfc = dynamic_cast<FlowContainer*>(tfi->Get(this->GetName()));
  if (!lfc) {
    printf("Could not pick up the %s from %s\n", this->GetName(), tfi->GetName());
    return;
  }
  TProfile2D* spro = lfc->GetProfile();
  TProfile2D* tpro = GetProfile();
  if (!tpro) {
    fProf = dynamic_cast<TProfile2D*>(spro->Clone(spro->GetName()));
    fProf->SetDirectory(0);
  } else {
    tpro->Add(spro);
  }
  TObjArray* tarr = lfc->GetSubProfiles();
  if (!tarr) {
    return;
  }
  if (!fProfRand) {
    fProfRand = new TObjArray();
    fProfRand->SetOwner(kTRUE);
  }
  for (int i = 0; i < tarr->GetEntries(); i++) {
    if (!(fProfRand->FindObject(tarr->At(i)->GetName()))) {
      fProfRand->Add(dynamic_cast<TProfile2D*>(tarr->At(i)->Clone(tarr->At(i)->GetName())));
      dynamic_cast<TProfile2D*>(fProfRand->At(fProfRand->GetEntries() - 1))->SetDirectory(0);
    } else {
      dynamic_cast<TProfile2D*>(fProfRand->FindObject(tarr->At(i)->GetName()))->Add(dynamic_cast<TProfile2D*>(tarr->At(i)));
    }
  }
}
bool FlowContainer::OverrideBinsWithZero(int xb1, int yb1, int xb2, int yb2)
{
  ProfileSubset* t_apf = new ProfileSubset(*fProf);
  if (!t_apf->OverrideBinsWithZero(xb1, yb1, xb2, yb2)) {
    delete t_apf;
    return kFALSE;
  }
  delete fProf;
  fProf = dynamic_cast<TProfile2D*>(t_apf);
  return kTRUE;
}
bool FlowContainer::OverrideMainWithSub(int ind, bool ExcludeChosen)
{
  if (!fProfRand) {
    printf("Cannot override main profile with a randomized one. Random profile array does not exist.\n");
    return kFALSE;
  }
  if (!ExcludeChosen) {
    TProfile2D* tarprof = dynamic_cast<TProfile2D*>(fProfRand->At(ind));
    if (!tarprof) {
      printf("Target random histogram does not exist.\n");
      return kFALSE;
    }
    TString ts(fProf->GetName());
    delete fProf;
    fProf = dynamic_cast<TProfile2D*>(tarprof->Clone(ts.Data()));
    return kTRUE;
  } else {
    TString ts(fProf->GetName());
    delete fProf;
    fProf = 0;
    for (int i = 0; i < fProfRand->GetEntries(); i++) {
      if (i == ind)
        continue;
      TProfile2D* tarprof = dynamic_cast<TProfile2D*>(fProfRand->At(i));
      if (!fProf)
        fProf = dynamic_cast<TProfile2D*>(tarprof->Clone(ts.Data()));
      else
        fProf->Add(tarprof);
    }
    return kTRUE;
  }
}
bool FlowContainer::RandomizeProfile(int nSubsets)
{
  if (!fProfRand) {
    printf("Cannot randomize profile, random array does not exist.\n");
    return kFALSE;
  }
  int l_Subsets = nSubsets ? nSubsets : fProfRand->GetEntries();
  TRandom* rndm = new TRandom(0);
  for (int i = 0; i < l_Subsets; i++) {
    int rInd = TMath::FloorNint(rndm->Rndm() * fProfRand->GetEntries());
    if (!i) {
      TString ts(fProf->GetName());
      delete fProf;
      fProf = dynamic_cast<TProfile2D*>(fProfRand->At(rInd)->Clone(ts.Data()));
    } else {
      fProf->Add(dynamic_cast<TProfile2D*>(fProfRand->At(rInd)));
    }
  }
  return kTRUE;
}
bool FlowContainer::CreateStatisticsProfile(StatisticsType StatType, int arg)
{
  switch (StatType) {
    case kSingleSample:
      return OverrideMainWithSub(arg, kFALSE);
      break;
    case kJackKnife:
      return OverrideMainWithSub(arg, kTRUE);
      break;
    case kBootstrap:
      return RandomizeProfile(arg);
      break;
    default:
      return kFALSE;
      break;
  }
}
void FlowContainer::SetIDName(TString newname)
{
  fIDName = newname;
}
TProfile* FlowContainer::GetCorrXXVsMulti(const char* order, int l_pti)
{
  TProfile* retSubset = 0;
  TString l_name("");
  Ssiz_t l_pos = 0;
  while (fIDName.Tokenize(l_name, l_pos)) {
    const char* ptpf = l_pti > 0 ? Form("_pt_%i", l_pti) : "";
    const char* ybinlab = Form("%s%s%s", l_name.Data(), order, ptpf);
    int ybinno = fProf->GetYaxis()->FindBin(ybinlab);
    if (ybinno < 0) {
      printf("Could not find %s!\n", ybinlab);
      return 0;
    }
    TProfile* rethist = dynamic_cast<TProfile*>(fProf->ProfileX("temp_prof", ybinno, ybinno));
    rethist->SetTitle(Form(";multi.;#LT#LT%s#GT#GT", order));
    if (!retSubset) {
      retSubset = dynamic_cast<TProfile*>(rethist->Clone(Form("corr_%s", order)));
    } else {
      retSubset->Add(rethist);
    }
    delete rethist;
  }
  if (fMultiRebin > 0) {
    TString temp_name(retSubset->GetName());
    TProfile* tempprof = dynamic_cast<TProfile*>(retSubset->Clone("tempProfile"));
    delete retSubset;
    retSubset = dynamic_cast<TProfile*>(tempprof->Rebin(fMultiRebin, temp_name.Data(), fMultiRebinEdges));
    delete tempprof;
  }
  return retSubset;
};
TH1D* FlowContainer::GetCorrXXVsPt(const char* order, double lminmulti, double lmaxmulti)
{
  int minm = 1;
  int maxm = fProf->GetXaxis()->GetNbins();
  if (!fbinsPt)
    SetXAxis();
  if (lminmulti > 0) {
    minm = fProf->GetXaxis()->FindBin(lminmulti + 0.001);
    maxm = minm;
  }
  if (lmaxmulti > lminmulti)
    maxm = fProf->GetXaxis()->FindBin(lmaxmulti - 0.001);
  ProfileSubset* rhProfSub = new ProfileSubset(*fProf);
  TString l_name("");
  Ssiz_t l_pos = 0;
  while (fIDName.Tokenize(l_name, l_pos)) {
    TString ybl1(Form("%s%s_pt_1", l_name.Data(), order));
    TString ybl2(Form("%s%s_pt_%i", l_name.Data(), order, fNbinsPt));
    int ybn1 = fProf->GetYaxis()->FindBin(ybl1.Data());
    int ybn2 = fProf->GetYaxis()->FindBin(ybl2.Data());
    if (fNbinsPt != (ybn2 - ybn1 + 1)) {
      printf("fNbinsPt is not matching the num of found histograms");
      return nullptr;
    }
    TProfile* profY = rhProfSub->ProfileY("profY", minm, maxm);
    TH1D* histY = ProfToHist(profY);
    TH1D* hist = new TH1D("temphist", "temphist", fNbinsPt, fbinsPt);
    for (int ibin = 1; ibin < hist->GetNbinsX(); ibin++) {
      TString bLabel = rhProfSub->GetYaxis()->GetBinLabel(ibin + ybn1 - 1);
      hist->GetXaxis()->SetBinLabel(ibin, bLabel.Data());
      hist->SetBinContent(ibin, histY->GetBinContent(ibin + ybn1 - 1));
      hist->SetBinError(ibin, histY->GetBinError(ibin + ybn1 - 1));
    }
    delete histY;
    delete rhProfSub;
    return hist;
  }
  return nullptr;
};
TH1D* FlowContainer::ProfToHist(TProfile* inpf)
{
  int nbins = inpf->GetNbinsX();
  double* xbs = new double[nbins + 1];
  inpf->GetLowEdge(xbs);
  xbs[nbins] = xbs[nbins - 1] + inpf->GetBinWidth(nbins);
  TH1D* rethist = new TH1D(Form("%s_hist", inpf->GetName()), inpf->GetTitle(), nbins, xbs);
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    if (inpf->GetBinContent(i) != 0) {
      rethist->SetBinContent(i, inpf->GetBinContent(i));
      rethist->SetBinError(i, inpf->GetBinError(i));
    }
  }
  return rethist;
}
TH1D* FlowContainer::GetHistCorrXXVsMulti(const char* order, int l_pti)
{
  TProfile* tpf = GetCorrXXVsMulti(order, l_pti);
  TH1D* rethist = ProfToHist(tpf);
  delete tpf;
  return rethist;
}
TH1D* FlowContainer::GetHistCorrXXVsPt(const char* order, double lminmulti, double lmaxmulti)
{
  TH1D* rethist = GetCorrXXVsPt(order, lminmulti, lmaxmulti);
  if (!rethist) {
    printf("GetCorrXXVsPt return nullptr!");
    return nullptr;
  }
  TProfile* refflow = GetRefFlowProfile(order, lminmulti, lmaxmulti);
  if (refflow) {
    refflow->RebinX(refflow->GetNbinsX());
    rethist->SetBinContent(0, refflow->GetBinContent(1));
    rethist->SetBinError(0, refflow->GetBinError(1));
  }
  delete refflow;
  return rethist;
}
TH1D* FlowContainer::GetVN2(TH1D* cn2)
{
  TH1D* rethist = dynamic_cast<TH1D*>(cn2->Clone(Form("vn2_%s", cn2->GetName())));
  rethist->Reset();
  double rf2 = cn2->GetBinContent(0);
  double rf2e = cn2->GetBinError(0);
  bool OnPt = ((!rf2) == 0);
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    double d2 = cn2->GetBinContent(i);
    double d2e = cn2->GetBinError(i);
    if (d2 > 0) {
      rethist->SetBinContent(i, OnPt ? VDN2Value(d2, rf2) : VN2Value(d2));
      rethist->SetBinError(i, OnPt ? VDN2Error(d2, d2e, rf2, rf2e) : VN2Error(d2, d2e));
    }
  }
  return rethist;
}
TH1D* FlowContainer::GetCN2VsX(int n, bool onPt, double arg1, double arg2)
{
  TH1D* corrN2;
  if (onPt)
    corrN2 = GetHistCorrXXVsPt(Form("%i2", n), arg1, arg2);
  else
    corrN2 = GetHistCorrXXVsMulti(Form("%i2", n), static_cast<int>(arg1));
  corrN2->SetName(Form("Corr_%s", corrN2->GetName()));
  TH1D* rethist = GetCN2(corrN2);
  TString* nam = new TString(corrN2->GetName());
  delete corrN2;
  rethist->SetName(nam->Data());
  if (onPt) {
    int bins = fProf->GetXaxis()->FindBin(arg1);
    int bins2 = fProf->GetXaxis()->FindBin(arg2);
    double bv1 = fProf->GetXaxis()->GetBinLowEdge(bins);
    double bv2 = fProf->GetXaxis()->GetBinUpEdge(bins2);
    rethist->SetTitle(Form("%2.0f - %4.0f;#it{p}_{T} (GeV/#it{c}); c_{%i}{2}", bv1, bv2, n));
  } else {
    rethist->SetTitle(Form(";#it{N}_{tr};c_{%i}{2}", n));
  }
  return rethist;
}

TH1D* FlowContainer::GetVN2VsX(int n, bool onPt, double arg1, double arg2)
{
  TH1D* corrh = GetCN2VsX(n, onPt, arg1, arg2);
  TString* nam = new TString(corrh->GetName());
  TH1D* rethist = GetVN2(corrh);
  delete corrh;
  rethist->SetName(nam->Data());
  if (onPt) {
    int bins = fProf->GetXaxis()->FindBin(arg1);
    int bins2 = fProf->GetXaxis()->FindBin(arg2);
    double bv1 = fProf->GetXaxis()->GetBinCenter(bins);
    double bv2 = fProf->GetXaxis()->GetBinCenter(bins2);
    rethist->SetTitle(Form("%2.0f - %4.0f;#it{p}_{T} (GeV/#it{c}); v_{%i}{2}", bv1, bv2, n));
  } else {
    rethist->SetTitle(Form(";#it{N}_{tr};v_{%i}{2}", n));
  }
  delete nam;
  return rethist;
}

TH1D* FlowContainer::GetCN2(TH1D* corrN2)
{
  double rf2 = corrN2->GetBinContent(0);
  double rf2e = corrN2->GetBinError(0);
  bool OnPt = (rf2 != 0);
  TH1D* rethist = dynamic_cast<TH1D*>(corrN2->Clone(Form("cN2_%s", corrN2->GetName())));
  rethist->Reset();
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    double cor2v = corrN2->GetBinContent(i);
    double cor2e = corrN2->GetBinError(i);
    rethist->SetBinContent(i, cor2v);
    rethist->SetBinError(i, cor2e);
  }
  if (OnPt) {
    rethist->SetBinContent(0, rf2);
    rethist->SetBinError(0, rf2e);
  } else {
    rethist->SetBinContent(0, 0);
    rethist->SetBinError(0, 0);
  }
  return rethist;
}

TH1D* FlowContainer::GetCN4(TH1D* corrN4, TH1D* corrN2)
{
  double rf2 = corrN2->GetBinContent(0);
  double rf2e = corrN2->GetBinError(0);
  double rf4 = corrN4->GetBinContent(0);
  double rf4e = corrN4->GetBinError(0);
  bool OnPt = (rf2 != 0);
  TH1D* rethist = dynamic_cast<TH1D*>(corrN4->Clone(Form("cN4_%s", corrN4->GetName())));
  rethist->Reset();
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    double cor4v = corrN4->GetBinContent(i);
    double cor4e = corrN4->GetBinError(i);
    double cor2v = corrN2->GetBinContent(i);
    double cor2e = corrN2->GetBinError(i);
    rethist->SetBinContent(i, OnPt ? DN4Value(cor4v, cor2v, rf2) : CN4Value(cor4v, cor2v));
    rethist->SetBinError(i, OnPt ? DN4Error(cor4e, cor2v, cor2e, rf2, rf2e) : CN4Error(cor4e, cor2v, cor2e));
  }
  if (OnPt) {
    rethist->SetBinContent(0, CN4Value(rf4, rf2));
    rethist->SetBinError(0, CN4Error(rf4e, rf2, rf2e));
  } else {
    rethist->SetBinContent(0, 0);
    rethist->SetBinError(0, 0);
  }
  return rethist;
}

TH1D* FlowContainer::GetCN6(TH1D* corrN6, TH1D* corrN4, TH1D* corrN2)
{
  TH1D* tn2 = dynamic_cast<TH1D*>(corrN2->Clone(Form("tn2_%s", corrN2->GetName())));
  TH1D* tn4 = dynamic_cast<TH1D*>(corrN4->Clone(Form("tn4_%s", corrN4->GetName())));
  TH1D* tn6 = dynamic_cast<TH1D*>(corrN6->Clone(Form("tn6_%s", corrN6->GetName())));

  double rf2 = corrN2->GetBinContent(0);
  double rf2e = corrN2->GetBinError(0);
  double rf4 = corrN4->GetBinContent(0);
  double rf4e = corrN4->GetBinError(0);
  double rf6 = corrN6->GetBinContent(0);
  double rf6e = corrN6->GetBinError(0);
  bool OnPt = (rf2 != 0);
  TH1D* rethist = dynamic_cast<TH1D*>(corrN6->Clone(Form("cN6_%s", corrN6->GetName())));
  rethist->Reset();
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    double cor6v = corrN6->GetBinContent(i);
    double cor6e = corrN6->GetBinError(i);
    double cor4v = corrN4->GetBinContent(i);
    double cor4e = corrN4->GetBinError(i);
    double cor2v = corrN2->GetBinContent(i);
    double cor2e = corrN2->GetBinError(i);
    rethist->SetBinContent(i, OnPt ? DN6Value(cor6v, cor4v, cor2v, rf4, rf2) : CN6Value(cor6v, cor4v, cor2v));
    rethist->SetBinError(i, OnPt ? DN6Error(cor6e, cor4v, cor4e, cor2v, cor2e, rf4, rf4e, rf2, rf2e) : CN6Error(cor6e, cor4v, cor4e, cor2v, cor2e));
  }
  if (OnPt) {
    rethist->SetBinContent(0, CN6Value(rf6, rf4, rf2));
    rethist->SetBinError(0, CN6Error(rf6e, rf4, rf4e, rf2, rf2e));
  } else {
    rethist->SetBinContent(0, 0);
    rethist->SetBinError(0, 0);
  }
  delete tn2;
  delete tn4;
  delete tn6;
  return rethist;
}
TH1D* FlowContainer::GetCN8(TH1D* corrN8, TH1D* corrN6, TH1D* corrN4, TH1D* corrN2)
{
  TH1D* tn2 = dynamic_cast<TH1D*>(corrN2->Clone(Form("tn2_%s", corrN2->GetName())));
  TH1D* tn4 = dynamic_cast<TH1D*>(corrN4->Clone(Form("tn4_%s", corrN4->GetName())));
  TH1D* tn6 = dynamic_cast<TH1D*>(corrN6->Clone(Form("tn6_%s", corrN6->GetName())));
  TH1D* tn8 = dynamic_cast<TH1D*>(corrN8->Clone(Form("tn8_%s", corrN6->GetName())));

  double rf2 = corrN2->GetBinContent(0);
  double rf2e = corrN2->GetBinError(0);
  double rf4 = corrN4->GetBinContent(0);
  double rf4e = corrN4->GetBinError(0);
  double rf6 = corrN6->GetBinContent(0);
  double rf6e = corrN6->GetBinError(0);
  double rf8 = corrN8->GetBinContent(0);
  double rf8e = corrN8->GetBinError(0);
  bool OnPt = (rf2 != 0);
  TH1D* rethist = dynamic_cast<TH1D*>(corrN8->Clone(Form("cN8_%s", corrN8->GetName())));
  rethist->Reset();
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    double cor8v = corrN8->GetBinContent(i);
    double cor8e = corrN8->GetBinError(i);
    double cor6v = corrN6->GetBinContent(i);
    double cor6e = corrN6->GetBinError(i);
    double cor4v = corrN4->GetBinContent(i);
    double cor4e = corrN4->GetBinError(i);
    double cor2v = corrN2->GetBinContent(i);
    double cor2e = corrN2->GetBinError(i);
    rethist->SetBinContent(i, OnPt ? DN8Value(cor8v, cor6v, cor4v, cor2v, rf6, rf4, rf2) : CN8Value(cor8v, cor6v, cor4v, cor2v));
    rethist->SetBinError(i, OnPt ? DN8Error(cor8e, cor6v, cor6e, cor4v, cor4e, cor2v, cor2e, rf6, rf6e, rf4, rf4e, rf2, rf2e) : CN8Error(cor8e, cor6v, cor6e, cor4v, cor4e, cor2v, cor2e));
  }
  if (OnPt) {
    rethist->SetBinContent(0, CN8Value(rf8, rf6, rf4, rf2));
    rethist->SetBinError(0, CN8Error(rf8e, rf6, rf6e, rf4, rf4e, rf2, rf2e));
  } else {
    rethist->SetBinContent(0, 0);
    rethist->SetBinError(0, 0);
  }
  delete tn2;
  delete tn4;
  delete tn6;
  delete tn8;
  return rethist;
}

TH1D* FlowContainer::GetCN4VsX(int n, bool onPt, double arg1, double arg2)
{
  TH1D *corrN2, *corrN4;
  if (onPt) {
    corrN2 = GetHistCorrXXVsPt(Form("%i2", n), arg1, arg2);
    corrN2->SetName(Form("Corr_%s", corrN2->GetName()));
    corrN4 = GetHistCorrXXVsPt(Form("%i4", n), arg1, arg2);
    corrN4->SetName(Form("Corr_%s", corrN4->GetName()));
  } else {
    corrN2 = GetHistCorrXXVsMulti(Form("%i2", n), static_cast<int>(arg1));
    corrN2->SetName(Form("Corr_%s", corrN2->GetName()));
    corrN4 = GetHistCorrXXVsMulti(Form("%i4", n), static_cast<int>(arg1));
    corrN4->SetName(Form("Corr_%s", corrN4->GetName()));
  }
  TH1D* rethist = GetCN4(corrN4, corrN2);
  TString* nam = new TString(corrN4->GetName());
  delete corrN2;
  delete corrN4;
  rethist->SetName(nam->Data());
  if (onPt) {
    int bins = fProf->GetXaxis()->FindBin(arg1);
    int bins2 = fProf->GetXaxis()->FindBin(arg2);
    double bv1 = fProf->GetXaxis()->GetBinLowEdge(bins);
    double bv2 = fProf->GetXaxis()->GetBinUpEdge(bins2);
    rethist->SetTitle(Form("%2.0f - %4.0f;#it{p}_{T} (GeV/#it{c}); c_{%i}{4}", bv1, bv2, n));
  } else {
    rethist->SetTitle(Form(";#it{N}_{tr};c_{%i}{4}", n));
  }
  return rethist;
}
TH1D* FlowContainer::GetVN4(TH1D* inh)
{
  TH1D* rethist = dynamic_cast<TH1D*>(inh->Clone(Form("v24_%s", inh->GetName())));
  double c4 = inh->GetBinContent(0);
  double c4e = inh->GetBinError(0);
  bool OnPt = ((!c4) == 0);
  rethist->Reset();
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    double d4 = inh->GetBinContent(i);
    double d4e = inh->GetBinError(i);
    if (OnPt && c4 >= 0)
      continue;
    // if(d4>=0) continue;
    rethist->SetBinContent(i, OnPt ? VDN4Value(d4, c4) : VN4Value(d4));
    rethist->SetBinError(i, OnPt ? VDN4Error(d4, d4e, c4, c4e) : VN4Error(d4, d4e));
  }
  return rethist;
}
TH1D* FlowContainer::GetVN6(TH1D* inh)
{
  TH1D* rethist = dynamic_cast<TH1D*>(inh->Clone(Form("v26_%s", inh->GetName())));
  double c6 = inh->GetBinContent(0);
  double c6e = inh->GetBinError(0);
  bool OnPt = ((!c6) == 0);
  rethist->Reset();
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    double d6 = inh->GetBinContent(i);
    double d6e = inh->GetBinError(i);
    if (OnPt && c6 <= 0)
      continue;
    // if(d6<=0) continue;
    rethist->SetBinContent(i, OnPt ? VDN6Value(d6, c6) : VN6Value(d6));
    rethist->SetBinError(i, OnPt ? VDN6Error(d6, d6e, c6, c6e) : VN6Error(d6, d6e));
  }
  return rethist;
}
TH1D* FlowContainer::GetVN8(TH1D* inh)
{
  TH1D* rethist = dynamic_cast<TH1D*>(inh->Clone(Form("v28_%s", inh->GetName())));
  double c8 = inh->GetBinContent(0);
  double c8e = inh->GetBinError(0);
  bool OnPt = ((!c8) == 0);
  rethist->Reset();
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    double d8 = inh->GetBinContent(i);
    double d8e = inh->GetBinError(i);
    if (OnPt && c8 > 0)
      continue;
    // if(d8>0) continue;
    rethist->SetBinContent(i, OnPt ? VDN8Value(d8, c8) : VN8Value(d8));
    rethist->SetBinError(i, OnPt ? VDN8Error(d8, d8e, c8, c8e) : VN8Error(d8, d8e));
  }
  return rethist;
}

TH1D* FlowContainer::GetVN4VsX(int n, bool onPt, double arg1, double arg2)
{
  TH1D* temph = GetCN4VsX(n, onPt, arg1, arg2);
  TH1D* rethist = GetVN4(temph);
  TString* nam = new TString(temph->GetName());
  delete temph;
  rethist->SetName(nam->Data());
  if (onPt) {
    int bins = fProf->GetXaxis()->FindBin(arg1);
    int bins2 = fProf->GetXaxis()->FindBin(arg2);
    double bv1 = fProf->GetXaxis()->GetBinLowEdge(bins);
    double bv2 = fProf->GetXaxis()->GetBinUpEdge(bins2);
    rethist->SetTitle(Form("%2.0f - %4.0f;#it{p}_{T} (GeV/#it{c}); v_{%i}{4}", bv1, bv2, n));
  } else {
    rethist->SetTitle(Form(";#it{N}_{tr};v_{%i}{4}", n));
  }
  return rethist;
}
TH1D* FlowContainer::GetCN6VsX(int n, bool onPt, double arg1, double arg2)
{
  TH1D *corrN2, *corrN4, *corrN6;
  if (onPt) {
    corrN2 = GetHistCorrXXVsPt(Form("%i2", n), arg1, arg2);
    corrN2->SetName(Form("Corr_%s", corrN2->GetName()));
    corrN4 = GetHistCorrXXVsPt(Form("%i4", n), arg1, arg2);
    corrN4->SetName(Form("Corr_%s", corrN4->GetName()));
    corrN6 = GetHistCorrXXVsPt(Form("%i6", n), arg1, arg2);
    corrN6->SetName(Form("Corr_%s", corrN6->GetName()));
  } else {
    corrN2 = GetHistCorrXXVsMulti(Form("%i2", n), static_cast<int>(arg1));
    corrN4 = GetHistCorrXXVsMulti(Form("%i4", n), static_cast<int>(arg1));
    corrN6 = GetHistCorrXXVsMulti(Form("%i6", n), static_cast<int>(arg1));
  }
  TH1D* rethist = GetCN6(corrN6, corrN4, corrN2);
  delete corrN2;
  delete corrN4;
  TString* nam = new TString(corrN6->GetName());
  delete corrN6;
  rethist->SetName(nam->Data());
  if (onPt) {
    int bins = fProf->GetXaxis()->FindBin(arg1);
    int bins2 = fProf->GetXaxis()->FindBin(arg2);
    double bv1 = fProf->GetXaxis()->GetBinLowEdge(bins);
    double bv2 = fProf->GetXaxis()->GetBinUpEdge(bins2);
    rethist->SetTitle(Form("%2.0f - %4.0f;#it{p}_{T} (GeV/#it{c}); c_{%i}{6}", bv1, bv2, n));
  } else {
    rethist->SetTitle(Form(";#it{N}_{tr};c_{%i}{6}", n));
  }
  return rethist;
}
TH1D* FlowContainer::GetVN6VsX(int n, bool onPt, double arg1, double arg2)
{
  TH1D* temph = GetCN6VsX(n, onPt, arg1, arg2);
  TH1D* rethist = GetVN6(temph);
  TString* nam = new TString(temph->GetName());
  delete temph;
  rethist->SetName(nam->Data());
  if (onPt) {
    int bins = fProf->GetXaxis()->FindBin(arg1);
    int bins2 = fProf->GetXaxis()->FindBin(arg2);
    double bv1 = fProf->GetXaxis()->GetBinLowEdge(bins);
    double bv2 = fProf->GetXaxis()->GetBinUpEdge(bins2);
    rethist->SetTitle(Form("%2.0f - %4.0f;#it{p}_{T} (GeV/#it{c}); v_{%i}{6}", bv1, bv2, n));
  } else {
    rethist->SetTitle(Form(";#it{N}_{tr};v_{%i}{6}", n));
  }
  return rethist;
}
TH1D* FlowContainer::GetCN8VsX(int n, bool onPt, double arg1, double arg2)
{
  TH1D *corrN2, *corrN4, *corrN6, *corrN8;
  if (onPt) {
    corrN2 = GetHistCorrXXVsPt(Form("%i2", n), arg1, arg2);
    corrN2->SetName(Form("Corr_%s", corrN2->GetName()));
    corrN4 = GetHistCorrXXVsPt(Form("%i4", n), arg1, arg2);
    corrN4->SetName(Form("Corr_%s", corrN4->GetName()));
    corrN6 = GetHistCorrXXVsPt(Form("%i6", n), arg1, arg2);
    corrN6->SetName(Form("Corr_%s", corrN6->GetName()));
    corrN8 = GetHistCorrXXVsPt(Form("%i8", n), arg1, arg2);
    corrN8->SetName(Form("Corr_%s", corrN8->GetName()));
  } else {
    corrN2 = GetHistCorrXXVsMulti(Form("%i2", n), static_cast<int>(arg1));
    corrN4 = GetHistCorrXXVsMulti(Form("%i4", n), static_cast<int>(arg1));
    corrN6 = GetHistCorrXXVsMulti(Form("%i6", n), static_cast<int>(arg1));
    corrN8 = GetHistCorrXXVsMulti(Form("%i8", n), static_cast<int>(arg1));
  }
  TH1D* rethist = GetCN8(corrN8, corrN6, corrN4, corrN2);
  delete corrN2;
  delete corrN4;
  delete corrN6;
  TString* nam = new TString(corrN8->GetName());
  delete corrN8;
  rethist->SetName(nam->Data());
  if (onPt) {
    int bins = fProf->GetXaxis()->FindBin(arg1);
    int bins2 = fProf->GetXaxis()->FindBin(arg2);
    double bv1 = fProf->GetXaxis()->GetBinLowEdge(bins);
    double bv2 = fProf->GetXaxis()->GetBinUpEdge(bins2);
    rethist->SetTitle(Form("%2.0f - %4.0f;#it{p}_{T} (GeV/#it{c}); c_{%i}{8}", bv1, bv2, n));
  } else {
    rethist->SetTitle(Form(";#it{N}_{tr};c_{%i}{8}", n));
  }
  return rethist;
}
TH1D* FlowContainer::GetVN8VsX(int n, bool onPt, double arg1, double arg2)
{
  TH1D* temph = GetCN8VsX(n, onPt, arg1, arg2);
  TH1D* rethist = GetVN8(temph);
  TString* nam = new TString(temph->GetName());
  delete temph;
  rethist->SetName(nam->Data());
  if (onPt) {
    int bins = fProf->GetXaxis()->FindBin(arg1);
    int bins2 = fProf->GetXaxis()->FindBin(arg2);
    double bv1 = fProf->GetXaxis()->GetBinLowEdge(bins);
    double bv2 = fProf->GetXaxis()->GetBinUpEdge(bins2);
    rethist->SetTitle(Form("%2.0f - %4.0f;#it{p}_{T} (GeV/#it{c}); v_{%i}{8}", bv1, bv2, n));
  } else {
    rethist->SetTitle(Form(";#it{N}_{tr};v_{%i}{8}", n));
  }
  return rethist;
}
TH1D* FlowContainer::GetCNN(int n, int c, bool onPt, double arg1, double arg2)
{
  if (c == 8)
    return GetCN8VsX(n, onPt, arg1, arg2);
  if (c == 6)
    return GetCN6VsX(n, onPt, arg1, arg2);
  if (c == 4)
    return GetCN4VsX(n, onPt, arg1, arg2);
  return GetCN2VsX(n, onPt, arg1, arg2);
};
TH1D* FlowContainer::GetVNN(int n, int c, bool onPt, double arg1, double arg2)
{
  if (c == 8)
    return GetVN8VsX(n, onPt, arg1, arg2);
  if (c == 6)
    return GetVN6VsX(n, onPt, arg1, arg2);
  if (c == 4)
    return GetVN4VsX(n, onPt, arg1, arg2);
  return GetVN2VsX(n, onPt, arg1, arg2);
};
TProfile* FlowContainer::GetRefFlowProfile(const char* order, double m1, double m2)
{
  int nStartBin = fProf->GetXaxis()->FindBin(m1 + 0.001);
  int nStopBin = fProf->GetXaxis()->FindBin(m2 - 0.001);
  if (nStartBin == 0)
    nStartBin = 1;
  if (nStopBin < nStartBin)
    nStopBin = fProf->GetXaxis()->GetNbins();
  int nBins = nStopBin - nStartBin + 1;
  double* l_bins = new double[nBins + 1];
  for (int i = 0; i <= nBins; i++)
    l_bins[i] = i;
  TProfile* retpf = 0;
  TString l_name("");
  Ssiz_t l_pos = 0;
  ProfileSubset* rhSubset = new ProfileSubset(*fProf);
  rhSubset->GetXaxis()->SetRange(nStartBin, nStopBin);
  while (fIDName.Tokenize(l_name, l_pos)) {
    l_name.Append(order);
    int ybin = fProf->GetYaxis()->FindBin(l_name.Data());
    TProfile* tempprof = rhSubset->GetSubset(kTRUE, "tempprof", ybin, ybin, nBins, l_bins);
    if (!retpf)
      retpf = dynamic_cast<TProfile*>(tempprof->Clone("RefFlowProf"));
    else
      retpf->Add(tempprof);
    delete tempprof;
  }
  delete rhSubset;
  if (!retpf) {
    LOGF(error, "Reference flow profile is null");
    return nullptr;
  } else {
    retpf->RebinX(nBins);
  }
  return retpf;
};

//{2} particle correlations
double FlowContainer::CN2Value(double cor2)
{
  if (!fPropagateErrors)
    return 0;
  return cor2;
};
double FlowContainer::CN2Error(double cor2e)
{
  if (!fPropagateErrors)
    return 0;
  return cor2e;
};
double FlowContainer::VN2Value(double cor2)
{
  if (cor2 < 0)
    return -2;
  return TMath::Sqrt(cor2);
};
double FlowContainer::VN2Error(double cor2, double cor2e)
{
  if (cor2 < 0)
    return 0;
  if (!fPropagateErrors)
    return 0;
  return 0.5 * cor2e / TMath::Sqrt(cor2);
};
double FlowContainer::VDN2Value(double cor2d, double cor2)
{
  if (cor2 < 0)
    return -2;
  return cor2d / TMath::Sqrt(cor2);
};
double FlowContainer::VDN2Error(double cor2d, double cor2de, double cor2, double cor2e)
{
  if (!fPropagateErrors)
    return 0;
  double sqrtv = cor2de * cor2de / cor2 + 0.25 * cor2d * cor2d * cor2e * cor2e / (cor2 * cor2 * cor2);
  if (sqrtv < 0)
    return 0;
  return TMath::Sqrt(sqrtv);
};

// C{4} and V{4} calculations
double FlowContainer::CN4Value(double cor4, double cor2)
{
  return cor4 - 2 * TMath::Power(cor2, 2);
};
double FlowContainer::CN4Error(double cor4e, double cor2, double cor2e)
{
  if (!fPropagateErrors)
    return 0;
  return TMath::Sqrt(cor4e * cor4e + 16 * cor2 * cor2 * cor2e * cor2e);
};
double FlowContainer::DN4Value(double cor4d, double cor2d, double cor2)
{
  return cor4d - 2 * cor2d * cor2;
};
double FlowContainer::DN4Error(double cor4de, double cor2d, double cor2de, double cor2, double cor2e)
{
  if (!fPropagateErrors)
    return 0;
  return TMath::Sqrt(cor4de * cor4de + 4 * cor2 * cor2 * cor2de * cor2de + 4 * cor2d * cor2d * cor2e * cor2e);
};
double FlowContainer::VN4Value(double c4)
{
  if (c4 > 0)
    return -2;
  return TMath::Power(-c4, 1. / 4);
};
double FlowContainer::VN4Error(double c4, double c4e)
{
  if (c4 > 0)
    return 0;
  if (!fPropagateErrors)
    return 0;
  return TMath::Power(-c4, (-3. / 4)) * c4e / 4;
};
double FlowContainer::VDN4Value(double d4, double c4)
{
  if (c4 > 0)
    return -2;
  return -d4 * TMath::Power(-c4, (-3. / 4));
};
double FlowContainer::VDN4Error(double d4, double d4e, double c4, double c4e)
{
  if (!fPropagateErrors)
    return 0;
  if (c4 > 0)
    return 0;
  return TMath::Sqrt(TMath::Power(-c4, -6. / 4) * d4e * d4e +
                     TMath::Power(-c4, -14. / 4) * d4 * d4 * c4e * c4e * 9. / 16);
};

//{6} particle correlations

double FlowContainer::CN6Value(double cor6, double cor4, double cor2)
{
  return cor6 - 9 * cor2 * cor4 + 12 * cor2 * cor2 * cor2;
};
double FlowContainer::CN6Error(double cor6e, double cor4, double cor4e, double cor2, double cor2e)
{
  if (!fPropagateErrors)
    return 0;
  double inters[3];
  inters[0] = cor6e;
  inters[1] = -9 * cor2 * cor4e;
  inters[2] = (-9 * cor4 + 36 * cor2 * cor2) * cor2e;
  double sum = 0;
  for (int i = 0; i < 3; i++)
    sum += (inters[i] * inters[i]);
  return TMath::Sqrt(sum);
};
double FlowContainer::DN6Value(double cor6d, double cor4d, double cor2d, double cor4, double cor2)
{
  return cor6d - 6 * cor4d * cor2 - 3 * cor2d * cor4 + 12 * cor2d * cor2 * cor2;
};
double FlowContainer::DN6Error(double d6e, double d4, double d4e, double d2,
                               double d2e, double c4, double c4e, double c2,
                               double c2e)
{
  if (!fPropagateErrors)
    return 0;
  double inters[5];
  inters[0] = d6e;
  inters[1] = -6 * c2 * d4e;
  inters[2] = (-3 * c4 + 12 * c2 * c2) * d2e;
  inters[3] = -3 * d2 * c4e;
  inters[4] = (-6 * d4 + 24 * d2 * c2) * c2e;
  double sum = 0;
  for (int i = 0; i < 5; i++)
    sum += (inters[i] * inters[i]);
  return TMath::Sqrt(sum);
};
double FlowContainer::VN6Value(double c6)
{
  if (c6 < 0)
    return -2;
  return TMath::Power(c6 / 4, 1. / 6);
};
double FlowContainer::VN6Error(double c6, double c6e)
{
  if (c6 < 0)
    return 0;
  if (!fPropagateErrors)
    return 0;
  return c6e / 6 * TMath::Power(4, -1. / 6) * TMath::Power(c6, -5. / 6);
};
double FlowContainer::VDN6Value(double d6, double c6)
{
  if (c6 < 0)
    return -2;
  return d6 * TMath::Power(4, -1. / 6) * TMath::Power(c6, -5. / 6);
};
double FlowContainer::VDN6Error(double d6, double d6e, double c6, double c6e)
{
  if (!fPropagateErrors)
    return 0;
  if (c6 < 0)
    return 0;
  if (d6 == 0)
    return 0;
  double vdn6 = VDN6Value(d6, c6);
  double dp = d6e / d6;
  double cp = 5 * c6e / 6;
  return vdn6 * TMath::Sqrt(dp * dp + cp * cp);
};

// {8} particle correlations

double FlowContainer::CN8Value(double cor8, double cor6, double cor4, double cor2)
{
  return cor8 - 16 * cor6 * cor2 - 18 * cor4 * cor4 + 144 * cor4 * cor2 * cor2 - 144 * cor2 * cor2 * cor2 * cor2;
};
double FlowContainer::CN8Error(double cor8e, double cor6, double cor6e,
                               double cor4, double cor4e, double cor2, double cor2e)
{
  if (!fPropagateErrors)
    return 0;
  double parts[4];
  parts[0] = cor8e;
  parts[1] = -16 * cor2 * cor6e;
  parts[2] = (-36 * cor4 + 144 * cor2 * cor2) * cor4e;
  parts[3] = (-16 * cor6 + 288 * cor4 * cor2 + 576 * cor2 * cor2 * cor2) * cor2e;
  double retval = 0;
  for (int i = 0; i < 4; i++)
    retval += TMath::Power(parts[i], 2);
  return TMath::Sqrt(retval);
};
double FlowContainer::DN8Value(double cor8d, double cor6d, double cor4d, double cor2d, double cor6, double cor4, double cor2)
{
  return cor8d - 12 * cor6d * cor2 - 4 * cor2d * cor6 - 18 * cor4d * cor4 + 72 * cor4d * cor2 * cor2 + 72 * cor4 * cor2 * cor2d - 144 * cor2d * cor2 * cor2 * cor2;
};
double FlowContainer::DN8Error(double d8e, double d6, double d6e, double d4,
                               double d4e, double d2, double d2e, double c6,
                               double c6e, double c4, double c4e, double c2,
                               double c2e)
{
  if (!fPropagateErrors)
    return 0;
  double parts[7];
  parts[0] = d8e;                             // d/d8'
  parts[1] = -12 * c2 * d6e;                  // d/d6'
  parts[2] = -4 * d2 * c6e;                   // d/d6
  parts[3] = (-16 * c4 + 72 * d2) * d4e;      // d/d4'
  parts[4] = (-16 * d4 + 72 * c2 * d2) * c4e; // d/d4
  parts[5] = (-4 * c6 + 72 * c4 * c2 - 144 * c2 * c2 * c2) * d2e;
  parts[6] = (-12 * d6 + 144 * d4 * c2 + 72 * c4 * d2 - 432 * d2 * c2 * c2) * c2e;
  double retval = 0;
  for (int i = 0; i < 7; i++)
    retval += TMath::Power(parts[i], 2);
  return TMath::Sqrt(retval);
};
double FlowContainer::VN8Value(double c8)
{
  if (c8 > 0)
    return -2; // Return -2 if not ok
  return TMath::Power(-c8 / 33, 1. / 8);
};
double FlowContainer::VN8Error(double c8, double c8e)
{
  if (c8 > 0)
    return 0;
  if (!fPropagateErrors)
    return 0;
  return c8e * 1. / (8 * c8) * VN8Value(c8);
};
double FlowContainer::VDN8Value(double d8, double c8)
{
  if (c8 > 0)
    return -2; // Return -2 if not OK
  return d8 / c8 * VN8Value(c8);
};
double FlowContainer::VDN8Error(double d8, double d8e, double c8, double c8e)
{
  if (c8 > 0)
    return 1;
  if (d8 == 0)
    return 1;
  if (!fPropagateErrors)
    return 0;
  double vdn8v = VDN8Value(d8, c8);
  double dd = d8e / d8;
  double dc = -7 * c8e / (8 * c8);
  return vdn8v * TMath::Sqrt(dd * dd + dc * dc);
};
void FlowContainer::SetPtRebin(int nbins, double* binedges)
{
  fPtRebin = nbins;
  fPtRebinEdges = binedges;
  return;
  int fPtRebin = 0;
  // double *lPtRebinEdges=binedges;
  if (!fbinsPt)
    SetXAxis();
  for (int i = 0; i < nbins; i++)
    if (binedges[i] < fbinsPt[0] || binedges[i] > fbinsPt[fNbinsPt - 1])
      continue;
    else
      fPtRebin++;
  if (fPtRebinEdges)
    delete[] fPtRebinEdges;
  fPtRebinEdges = new double[fPtRebin];
  fPtRebin = 0;
  for (int i = 0; i < nbins; i++)
    if (binedges[i] < fbinsPt[0] || binedges[i] > fbinsPt[fNbinsPt])
      continue;
    else
      fPtRebinEdges[fPtRebin++] = binedges[i];
  // fPtRebin--;
}
void FlowContainer::SetMultiRebin(int nbins, double* binedges)
{
  if (fMultiRebinEdges) {
    delete[] fMultiRebinEdges;
    fMultiRebinEdges = 0;
  }
  if (nbins <= 0) {
    fMultiRebin = 0;
    return;
  }
  fMultiRebin = nbins;
  fMultiRebinEdges = new double[nbins + 1];
  for (int i = 0; i <= fMultiRebin; i++)
    fMultiRebinEdges[i] = binedges[i];
}
double* FlowContainer::GetMultiRebin(int& nbins)
{
  if (fMultiRebin <= 0) {
    nbins = 0;
    return 0;
  }
  nbins = fMultiRebin;
  double* retBins = new double[fMultiRebin + 1];
  for (int i = 0; i <= nbins; i++)
    retBins[i] = fMultiRebinEdges[i];
  return fMultiRebinEdges;
}
