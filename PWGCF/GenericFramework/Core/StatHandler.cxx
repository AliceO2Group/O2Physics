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

// NOLINTBEGIN(readability/casting, readability/braces)
#include "StatHandler.h"

#include "TRandom.h"

#include <cstdio>
#include <memory>
#include <vector>

StatHandler::StatHandler() : nSamples(10)
{
  methodMap[StatMethod::kBootstrap] = &StatHandler::bootstrap;
  methodMap[StatMethod::kVarianceOfSamples] = &StatHandler::varianceOfSamples;
  methodMap[StatMethod::kJackknife] = &StatHandler::jackknife;
}
StatHandler::StatHandler(int n) : nSamples(n)
{
  methodMap[StatMethod::kBootstrap] = &StatHandler::bootstrap;
  methodMap[StatMethod::kVarianceOfSamples] = &StatHandler::varianceOfSamples;
  methodMap[StatMethod::kJackknife] = &StatHandler::jackknife;
}

StatHandler::~StatHandler() {};

TH1D* StatHandler::getNominal(const std::vector<TH1D*>& samplehists, StatMethod method)
{
  auto stat_profile = (this->*methodMap[method])(samplehists);
  auto reth = samplehists[0];
  applyStatErrors(reth, stat_profile);
  delete stat_profile;
  return reth;
}
void StatHandler::applyStatErrors(TH1D* nominal, TProfile*& stat_profile)
{
  if (nominal->GetNbinsX() != stat_profile->GetNbinsX()) {
    printf("Bootstrap warning: number of bins in histogram and profile is not the same! Not applying...\n");
    return;
  };
  for (Int_t i = 1; i <= nominal->GetNbinsX(); i++) {
    if (nominal->GetBinContent(i) == 0)
      continue;
    nominal->SetBinError(i, stat_profile->GetBinError(i));
  };
  return;
}
TProfile* StatHandler::varianceOfSamples(const std::vector<TH1D*>& samplehists)
{
  setCurrentAxis(samplehists[0]);
  auto retpf = new TProfile("Statistics_Profile", "StatProf", nxbins, &(xbins[0]));
  for (Int_t i = 1; i <= nxbins; i++)
    for (Int_t j = 1; j <= nSamples; j++)
      retpf->Fill(samplehists[0]->GetBinCenter(i), samplehists[j]->GetBinContent(i));
  return retpf;
}
TProfile* StatHandler::bootstrap(const std::vector<TH1D*>& samplehists)
{
  setCurrentAxis(samplehists[0]);
  auto rnd = new TRandom();
  auto retpf = new TProfile("Statistics_Profile", "StatProf", nxbins, &(xbins[0]));
  for (int i(0); i < 1000; ++i) {
    std::vector<TH1D*> random_sample_hists;
    for (Int_t j = 0; j < nSamples; ++j) {
      int rnd_int = (int)rnd->Uniform(nSamples) + 1;
      auto random_sample = samplehists[rnd_int];
      random_sample->SetName(Form("%s_sample%i", random_sample->GetName(), j));
      random_sample_hists.push_back(random_sample);
    }
    auto tmppf = std::make_unique<TProfile>("tmp_Profile", "tmpProf", nxbins, &(xbins[0]));
    for (Int_t i = 1; i <= nxbins; ++i)
      for (Int_t j = 0; j < nSamples; j++)
        tmppf->Fill(samplehists[0]->GetBinCenter(i), random_sample_hists[j]->GetBinContent(i));
    for (Int_t i = 1; i <= nxbins; ++i) {
      retpf->Fill(samplehists[0]->GetBinCenter(i), tmppf->GetBinContent(i));
    }
  }

  delete rnd;
  return retpf;
}
TProfile* StatHandler::jackknife(const std::vector<TH1D*>& samplehists)
{
  setCurrentAxis(samplehists[0]);
  auto retpf = new TProfile("Statistics_Profile", "StatProf", nxbins, &(xbins[0]));
  // calculate jackknife means
  for (int i = 0; i < nSamples; ++i) {
    std::vector<TH1D*> jackknife_hists = samplehists;
    // remove the nominal hist
    jackknife_hists.erase(jackknife_hists.begin());
    // remove ith subsample hist
    jackknife_hists.erase(jackknife_hists.begin() + i);
    auto tmppf = std::make_unique<TProfile>("tmp_Profile", "tmpProf", nxbins, &(xbins[0]));
    for (int bin = 1; bin <= nxbins; ++bin) {
      for (size_t j = 0; j < jackknife_hists.size(); ++j) {
        tmppf->Fill(samplehists[0]->GetBinCenter(bin), jackknife_hists[j]->GetBinContent(bin));
      }
    }
    // jackknife estimator of mean as mean of jackknife subsample means
    for (int bin = 1; bin <= nxbins; ++bin) {
      retpf->Fill(samplehists[0]->GetBinCenter(bin), tmppf->GetBinContent(bin));
    }
  }
  return retpf;
}
void StatHandler::setCurrentAxis(TH1D* tmph)
{
  nxbins = tmph->GetNbinsX();
  double* tmp_xbins = new double[nxbins + 1];
  tmph->GetXaxis()->GetLowEdge(tmp_xbins);
  tmp_xbins[nxbins] = tmph->GetXaxis()->GetBinUpEdge(nxbins);
  xbins.clear();
  xbins.insert(xbins.end(), &tmp_xbins[0], &tmp_xbins[nxbins + 1]);
  delete[] tmp_xbins;
  return;
}
// NOLINTEND(readability/casting, readability/braces)
