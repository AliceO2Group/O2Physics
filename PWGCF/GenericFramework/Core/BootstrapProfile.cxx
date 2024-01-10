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

#include "BootstrapProfile.h"
BootstrapProfile::BootstrapProfile() : TProfile(),
                                       fListOfEntries(0),
                                       fProfInitialized(kFALSE),
                                       fNSubs(0),
                                       fMultiRebin(0),
                                       fMultiRebinEdges(0),
                                       fPresetWeights(0) {}
BootstrapProfile::~BootstrapProfile()
{
  delete fListOfEntries;
};
BootstrapProfile::BootstrapProfile(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins) : TProfile(name, title, nbinsx, xbins),
                                                                                                               fListOfEntries(0),
                                                                                                               fProfInitialized(kTRUE),
                                                                                                               fNSubs(0),
                                                                                                               fMultiRebin(0),
                                                                                                               fMultiRebinEdges(0),
                                                                                                               fPresetWeights(0) {}
BootstrapProfile::BootstrapProfile(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup) : TProfile(name, title, nbinsx, xlow, xup),
                                                                                                                     fListOfEntries(0),
                                                                                                                     fProfInitialized(kFALSE),
                                                                                                                     fNSubs(0),
                                                                                                                     fMultiRebin(0),
                                                                                                                     fMultiRebinEdges(0),
                                                                                                                     fPresetWeights(0) {}
void BootstrapProfile::InitializeSubsamples(Int_t nSub)
{
  if (nSub < 1) {
    printf("Number of subprofiles has to be > 0!\n");
    return;
  }
  if (fListOfEntries)
    delete fListOfEntries;
  fListOfEntries = new TList();
  fListOfEntries->SetOwner(kTRUE);
  TProfile* dummyPF = reinterpret_cast<TProfile*>(this);
  for (Int_t i = 0; i < nSub; i++) {
    fListOfEntries->Add(reinterpret_cast<TProfile*>(dummyPF->Clone(Form("%s_Subpf%i", dummyPF->GetName(), i))));
    reinterpret_cast<TProfile*>(fListOfEntries->At(i))->Reset();
  }
  fNSubs = nSub;
}
void BootstrapProfile::FillProfile(const Double_t& xv, const Double_t& yv, const Double_t& w, const Double_t& rn)
{
  TProfile::Fill(xv, yv, w);
  if (!fNSubs)
    return;
  Int_t targetInd = rn * fNSubs;
  if (targetInd >= fNSubs)
    targetInd = 0;
  reinterpret_cast<TProfile*>(fListOfEntries->At(targetInd))->Fill(xv, yv, w);
}
void BootstrapProfile::FillProfile(const Double_t& xv, const Double_t& yv, const Double_t& w)
{
  TProfile::Fill(xv, yv, w);
}
void BootstrapProfile::RebinMulti(Int_t nbins)
{
  this->RebinX(nbins);
  if (!fListOfEntries)
    return;
  for (Int_t i = 0; i < fListOfEntries->GetEntries(); i++)
    reinterpret_cast<TProfile*>(fListOfEntries->At(i))->RebinX(nbins);
}
TH1* BootstrapProfile::getHist(Int_t ind)
{
  if (fPresetWeights && fMultiRebin > 0)
    return getWeightBasedRebin(ind);
  if (ind < 0) {
    if (reinterpret_cast<TProfile*>(this)) {
      return getHistRebinned(reinterpret_cast<TProfile*>(this)); //((TProfile*)this)->ProjectionX(Form("%s_hist",this->GetName()));
    } else {
      printf("Empty BootstrapProfile addressed, cannot get a histogram\n");
      return 0;
    }
  } else {
    if (!fListOfEntries) {
      printf("No subprofiles exist!\n");
      return 0;
    }
    if (ind < fNSubs) {
      return getHistRebinned(reinterpret_cast<TProfile*>(fListOfEntries->At(ind))); ////((TProfile*)fListOfEntries->At(ind))->ProjectionX(Form("%s_sub%i",((TProfile*)fListOfEntries->At(ind))->GetName(),ind));
    } else {
      printf("Trying to fetch subprofile no %i out of %i, not possible\n", ind, fNSubs);
      return 0;
    }
  }
  return 0;
}
TProfile* BootstrapProfile::getProfile(Int_t ind)
{
  if (ind < 0) {
    if (reinterpret_cast<TProfile*>(this)) {
      return reinterpret_cast<TProfile*>(this);
    } else {
      printf("Empty BootstrapProfile addressed, cannot get a histogram\n");
      return 0;
    }
  } else {
    if (!fListOfEntries) {
      printf("No subprofiles exist!\n");
      return 0;
    }
    if (ind < fNSubs) {
      return reinterpret_cast<TProfile*>(fListOfEntries->At(ind));
    } else {
      printf("Trying to fetch subprofile no %i out of %i, not possible\n", ind, fNSubs);
      return 0;
    }
  }
}
Long64_t BootstrapProfile::Merge(TCollection* collist)
{
  Long64_t nmergedpf = TProfile::Merge(collist);
  Long64_t nmerged = 0;
  BootstrapProfile* l_PBS = 0;
  TIter all_PBS(collist);
  while ((l_PBS = reinterpret_cast<BootstrapProfile*>(all_PBS()))) {
    reinterpret_cast<TProfile*>(this)->Add(reinterpret_cast<TProfile*>(l_PBS));
    TList* tarL = l_PBS->fListOfEntries;
    if (!tarL)
      continue;
    if (!fListOfEntries) {
      fListOfEntries = reinterpret_cast<TList*>(tarL->Clone());
      for (Int_t i = 0; i < fListOfEntries->GetEntries(); i++)
        reinterpret_cast<TProfile*>(fListOfEntries->At(i))->Reset();
    }
    for (Int_t i = 0; i < fListOfEntries->GetEntries(); i++)
      reinterpret_cast<TProfile*>(fListOfEntries->At(i))->Add(reinterpret_cast<TProfile*>(tarL->At(i)));
    nmerged++;
  }
  return nmergedpf;
};
void BootstrapProfile::RebinMulti(Int_t nbins, Double_t* binedges)
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
  fMultiRebinEdges = new Double_t[nbins + 1];
  for (Int_t i = 0; i <= fMultiRebin; i++)
    fMultiRebinEdges[i] = binedges[i];
}
TH1* BootstrapProfile::getHistRebinned(TProfile* inpf)
{
  if (!inpf)
    return 0;
  if (fMultiRebin <= 0)
    return reinterpret_cast<TProfile*>(inpf)->ProjectionX(Form("%s_hist", inpf->GetName()));
  TProfile* temppf = reinterpret_cast<TProfile*>(inpf->Rebin(fMultiRebin, "tempProfile", fMultiRebinEdges));
  TH1* reth = reinterpret_cast<TH1*>(temppf->ProjectionX(Form("%s_hist", inpf->GetName())));
  delete temppf;
  return reth;
}
TH1* BootstrapProfile::getWeightBasedRebin(Int_t ind)
{
  if (!fPresetWeights) {
    printf("Weights are not preset!\n");
    return 0;
  }
  TProfile* lProf = getProfile(ind);
  TH1* reth = getHistRebinned(lProf);
  reth->Reset();
  TProfile* lW = fPresetWeights->getProfile(ind);
  if (!lW) {
    printf("Weight profile could not be found!\n");
    return 0;
  }
  for (Int_t i = 1; i <= lW->GetNbinsX(); i++) {
    Int_t i_n = reth->FindBin(lW->GetBinCenter(i));
    Double_t bc2 = lProf->GetBinContent(i);
    Double_t be2 = lW->GetBinEntries(i);
    Double_t bc1 = reth->GetBinContent(i_n);
    Double_t be1 = reth->GetBinError(i_n);
    if (be2 == 0)
      continue;
    reth->SetBinContent(i_n, bc1 + bc2 * be2);
    reth->SetBinError(i_n, be1 + be2);
  }
  for (Int_t i = 1; i <= reth->GetNbinsX(); i++) {
    Double_t bc1 = reth->GetBinContent(i);
    Double_t be1 = reth->GetBinError(i);
    if (be1 == 0)
      continue;
    reth->SetBinContent(i, bc1 / be1);
    reth->SetBinError(i, 1. / TMath::Sqrt(be1));
  }
  return reth;
}
void BootstrapProfile::MergeBS(BootstrapProfile* target)
{
  this->Add(target);
  TList* tarL = target->fListOfEntries;
  if (!fListOfEntries) {
    if (!target->fListOfEntries)
      return;
    fListOfEntries = reinterpret_cast<TList*>(tarL->Clone());
    for (Int_t i = 0; i < fListOfEntries->GetEntries(); i++)
      reinterpret_cast<TProfile*>(fListOfEntries->At(i))->Reset();
  }
  for (Int_t i = 0; i < fListOfEntries->GetEntries(); i++)
    reinterpret_cast<TProfile*>(fListOfEntries->At(i))->Add(reinterpret_cast<TProfile*>(tarL->At(i)));
}
TProfile* BootstrapProfile::getSummedProfiles()
{
  if (!fListOfEntries || !fListOfEntries->GetEntries()) {
    printf("No subprofiles initialized for the BootstrapProfile.\n");
    return 0;
  }
  TProfile* retpf = reinterpret_cast<TProfile*>(fListOfEntries->At(0)->Clone("SummedProfile"));
  for (Int_t i = 1; i < fListOfEntries->GetEntries(); i++)
    retpf->Add(reinterpret_cast<TProfile*>(fListOfEntries->At(i)));
  return retpf;
}
void BootstrapProfile::OverrideMainWithSub()
{
  TProfile* sum = getSummedProfiles();
  if (!sum)
    return;
  reinterpret_cast<TProfile*>(this)->Reset();
  reinterpret_cast<TProfile*>(this)->Add(sum);
  delete sum;
}
