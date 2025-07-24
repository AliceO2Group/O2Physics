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

/// \file BootstrapProfile.h/.cxx
/// \brief Derived class from TProfile that stores extra TProfiles for bootstrap samples
/// \author Emil Gorm Nielsen (ack. V. Vislavicius), NBI, emil.gorm.nielsen@cern.ch

#ifndef PWGCF_GENERICFRAMEWORK_CORE_BOOTSTRAPPROFILE_H_
#define PWGCF_GENERICFRAMEWORK_CORE_BOOTSTRAPPROFILE_H_

#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TCollection.h"
#include "TMath.h"

class BootstrapProfile : public TProfile
{
 public:
  BootstrapProfile();
  ~BootstrapProfile();
  BootstrapProfile(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins);
  BootstrapProfile(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup);
  TList* fListOfEntries;
  void MergeBS(BootstrapProfile* target);
  void InitializeSubsamples(Int_t nSub);
  void FillProfile(const Double_t& xv, const Double_t& yv, const Double_t& w, const Double_t& rn);
  void FillProfile(const Double_t& xv, const Double_t& yv, const Double_t& w);
  Long64_t Merge(TCollection* collist);
  void RebinMulti(Int_t nbins);
  void RebinMulti(Int_t nbins, Double_t* binedges);
  TH1* getHist(Int_t ind = -1);
  TProfile* getProfile(Int_t ind = -1);
  TProfile* getSummedProfiles();
  void OverrideMainWithSub();
  Int_t getNSubs() { return fListOfEntries->GetEntries(); }
  void PresetWeights(BootstrapProfile* targetBS) { fPresetWeights = targetBS; }
  void ResetBin(Int_t nbin)
  {
    ResetBin(reinterpret_cast<TProfile*>(this), nbin);
    for (Int_t i = 0; i < fListOfEntries->GetEntries(); i++)
      ResetBin(reinterpret_cast<TProfile*>(fListOfEntries->At(i)), nbin);
  };
  ClassDef(BootstrapProfile, 2);

 protected:
  TH1* getHistRebinned(TProfile* inpf); // Performs rebinning, if required, and returns a projection of profile
  TH1* getWeightBasedRebin(Int_t ind = -1);
  Bool_t fProfInitialized;
  Int_t fNSubs;
  Int_t fMultiRebin;                //! externaly set runtime, no need to store
  Double_t* fMultiRebinEdges;       //! externaly set runtime, no need to store
  BootstrapProfile* fPresetWeights; //! BootstrapProfile whose weights we should copy
  void ResetBin(TProfile* tpf, Int_t nbin)
  {
    tpf->SetBinEntries(nbin, 0);
    tpf->SetBinContent(nbin, 0);
    tpf->SetBinError(nbin, 0);
  };
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_BOOTSTRAPPROFILE_H_
