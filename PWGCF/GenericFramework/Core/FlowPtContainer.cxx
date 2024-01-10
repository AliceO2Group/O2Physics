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

#include "FlowPtContainer.h"

FlowPtContainer::FlowPtContainer() : TNamed("name", "name"),
                                     fCMTermList(0),
                                     fCorrList(0),
                                     fCovList(0),
                                     fCumulantList(0),
                                     fCentralMomentList(0),
                                     mpar(0),
                                     fillCounter(0),
                                     fEventWeight(kEventWeight::kUnity),
                                     corrNum(),
                                     corrDen() {}
FlowPtContainer::~FlowPtContainer()
{
  delete fCMTermList;
  delete fCorrList;
};
FlowPtContainer::FlowPtContainer(const char* name) : TNamed(name, name),
                                                     fCMTermList(0),
                                                     fCorrList(0),
                                                     fCovList(0),
                                                     fCumulantList(0),
                                                     fCentralMomentList(0),
                                                     mpar(0),
                                                     fillCounter(0),
                                                     fEventWeight(kEventWeight::kUnity),
                                                     corrNum(),
                                                     corrDen() {}
FlowPtContainer::FlowPtContainer(const char* name, const char* title, int nbinsx, double* xbins, const int& m, const GFWCorrConfigs& configs) : TNamed(name, title),
                                                                                                                                                fCMTermList(0),
                                                                                                                                                fCorrList(0),
                                                                                                                                                fCovList(0),
                                                                                                                                                fCumulantList(0),
                                                                                                                                                fCentralMomentList(0),
                                                                                                                                                mpar(m),
                                                                                                                                                fillCounter(0),
                                                                                                                                                fEventWeight(kEventWeight::kUnity),
                                                                                                                                                corrNum(),
                                                                                                                                                corrDen()
{
  Initialise(nbinsx, xbins, m, configs);
};
FlowPtContainer::FlowPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, const int& m, const GFWCorrConfigs& configs) : TNamed(name, title),
                                                                                                                                                            fCMTermList(0),
                                                                                                                                                            fCorrList(0),
                                                                                                                                                            fCovList(0),
                                                                                                                                                            fCumulantList(0),
                                                                                                                                                            fCentralMomentList(0),
                                                                                                                                                            mpar(m),
                                                                                                                                                            fillCounter(0),
                                                                                                                                                            fEventWeight(kEventWeight::kUnity),
                                                                                                                                                            corrNum(),
                                                                                                                                                            corrDen()
{
  Initialise(nbinsx, xlow, xhigh, m, configs);
};
void FlowPtContainer::Initialise(const o2::framework::AxisSpec axis, const int& m, const GFWCorrConfigs& configs, const int& nsub)
{
  if (!mpar)
    mpar = m;
  std::vector<double> multiBins = axis.binEdges;
  int nMultiBins = axis.nBins.value_or(0);
  if (nMultiBins <= 0)
    nMultiBins = multiBins.size() - 1;
  if (nMultiBins <= 0) {
    printf("Multiplicity axis does not exist");
    return;
  }
  if (fCMTermList)
    delete fCMTermList;
  fCMTermList = new TList();
  fCMTermList->SetOwner(kTRUE);
  if (fCorrList)
    delete fCorrList;
  fCorrList = new TList();
  fCorrList->SetOwner(kTRUE);
  if (fCovList)
    delete fCovList;
  fCovList = new TList();
  fCovList->SetOwner(kTRUE);
  for (int m = 0; m < mpar; ++m)
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("corr_%ipar", m + 1), nMultiBins, &multiBins[0]));
  for (int m = 0; m < 4; ++m) {
    for (int i = 0; i <= m; ++i)
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nMultiBins, &multiBins[0]));
  }
  for (int i = 0; i < configs.GetSize(); ++i) {
    for (auto m(1); m <= mpar; ++m) {
      if (!(configs.GetpTCorrMasks()[i] & (1 << (m - 1))))
        continue;
      fCovList->Add(new BootstrapProfile(Form("%s_mpt%i", configs.GetHeads()[i].c_str(), m), Form("%s_mpt%i", configs.GetHeads()[i].c_str(), m), nMultiBins, &multiBins[0]));
    }
  }
  if (nsub) {
    for (int i = 0; i < fCorrList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->InitializeSubsamples(nsub);
    for (int i = 0; i < fCMTermList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->InitializeSubsamples(nsub);
    for (int i = 0; i < fCovList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->InitializeSubsamples(nsub);
  }
  printf("Container %s initialized with m = %i\n and %i subsamples", this->GetName(), mpar, nsub);
  return;
};
void FlowPtContainer::Initialise(int nbinsx, double* xbins, const int& m, const GFWCorrConfigs& configs, const int& nsub)
{
  if (!mpar)
    mpar = m;
  if (fCMTermList)
    delete fCMTermList;
  fCMTermList = new TList();
  fCMTermList->SetOwner(kTRUE);
  if (fCorrList)
    delete fCorrList;
  fCorrList = new TList();
  fCorrList->SetOwner(kTRUE);
  for (int m = 0; m < mpar; ++m)
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("mpt%i", m + 1), nbinsx, xbins));
  for (int m = 0; m < 4; ++m) {
    for (int i = 0; i <= m; ++i)
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nbinsx, xbins));
  }
  for (int i = 0; i < configs.GetSize(); ++i) {
    for (auto m(1); m <= mpar; ++m) {
      if (!(configs.GetpTCorrMasks()[i] & (1 << (m - 1))))
        continue;
      fCovList->Add(new BootstrapProfile(Form("%s_mpt%i", configs.GetHeads()[i].c_str(), m + 1), Form("%s_mpt%i", configs.GetHeads()[i].c_str(), m + 1), nbinsx, xbins));
    }
  }
  if (nsub) {
    for (int i = 0; i < fCorrList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->InitializeSubsamples(nsub);
    for (int i = 0; i < fCMTermList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->InitializeSubsamples(nsub);
    for (int i = 0; i < fCovList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->InitializeSubsamples(nsub);
  }
  printf("Container %s initialized with m = %i\n", this->GetName(), mpar);
};
void FlowPtContainer::Initialise(int nbinsx, double xlow, double xhigh, const int& m, const GFWCorrConfigs& configs, const int& nsub)
{
  if (!mpar)
    mpar = m;
  if (fCMTermList)
    delete fCMTermList;
  fCMTermList = new TList();
  fCMTermList->SetOwner(kTRUE);
  if (fCorrList)
    delete fCorrList;
  fCorrList = new TList();
  fCorrList->SetOwner(kTRUE);
  for (int m = 0; m < mpar; ++m)
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("mpt%i", m + 1), nbinsx, xlow, xhigh));
  for (int m = 0; m < 4; ++m) {
    for (int i = 0; i <= m; ++i)
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nbinsx, xlow, xhigh));
  }
  for (int i = 0; i < configs.GetSize(); ++i) {
    for (auto m(1); m <= mpar; ++m) {
      if (!(configs.GetpTCorrMasks()[i] & (1 << (m - 1))))
        continue;
      fCovList->Add(new BootstrapProfile(Form("%s_mpt%i", configs.GetHeads()[i].c_str(), m + 1), Form("%s_mpt%i", configs.GetHeads()[i].c_str(), m + 1), nbinsx, xlow, xhigh));
    }
  }
  if (nsub) {
    for (int i = 0; i < fCorrList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->InitializeSubsamples(nsub);
    for (int i = 0; i < fCMTermList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->InitializeSubsamples(nsub);
    for (int i = 0; i < fCovList->GetEntries(); ++i)
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->InitializeSubsamples(nsub);
  }
  printf("Container %s initialized with m = %i\n", this->GetName(), mpar);
};
void FlowPtContainer::Fill(const double& w, const double& pt)
{
  for (auto i = 0; i < sumP.size(); ++i) {
    sumP[i] += pow(w, i % (mpar + 1)) * pow(pt, i / (mpar + 1));
  }
  return;
}
void FlowPtContainer::CalculateCorrelations()
{
  corrNum.clear();
  corrNum.resize(mpar + 1, 0);
  corrNum[0] = 1.0;
  corrDen.clear();
  corrDen.resize(mpar + 1, 0);
  corrDen[0] = 1.0;
  double sumNum = 0;
  double sumDenum = 0;
  std::vector<double> valNum;
  std::vector<double> valDenum;
  for (int m(1); m <= mpar; ++m) {
    for (int k(1); k <= m; ++k) {
      valNum.push_back(fSign[k - 1] * corrNum[m - k] * (fFactorial[m - 1] / fFactorial[m - k]) * sumP[GetVectorIndex(k, k)]);
      valDenum.push_back(fSign[k - 1] * corrDen[m - k] * (fFactorial[m - 1] / fFactorial[m - k]) * sumP[GetVectorIndex(k, 0)]);
    }
    sumNum = OrderedAddition(valNum);
    sumDenum = OrderedAddition(valDenum);
    valNum.clear();
    valDenum.clear();

    corrNum[m] = sumNum;
    corrDen[m] = sumDenum;
  }
  return;
}
void FlowPtContainer::FillPtProfiles(const double& centmult, const double& rn)
{
  for (int m = 1; m <= mpar; ++m) {
    if (corrDen[m] != 0)
      dynamic_cast<BootstrapProfile*>(fCorrList->At(m - 1))->FillProfile(centmult, corrNum[m] / corrDen[m], (fEventWeight == kEventWeight::kUnity) ? 1.0 : corrDen[m], rn);
  }
  return;
}
void FlowPtContainer::FillVnPtProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
{
  if (!mask)
    return;
  for (auto m(1); m <= mpar; ++m) {
    if (!(mask & (1 << (m - 1))))
      continue;
    if (corrDen[m] != 0)
      dynamic_cast<BootstrapProfile*>(fCovList->At(fillCounter))->FillProfile(centmult, flowval * corrNum[m] / corrDen[m], (fEventWeight == kUnity) ? 1.0 : flowtuples * corrDen[m], rn);
    ++fillCounter;
  }
  return;
}
void FlowPtContainer::FillCMProfiles(const double& centmult, const double& rn)
{
  if (sumP[GetVectorIndex(0, 0)] == 0)
    return;
  double tau1 = sumP[GetVectorIndex(2, 0)] / pow(sumP[GetVectorIndex(1, 0)], 2);
  double tau2 = sumP[GetVectorIndex(3, 0)] / pow(sumP[GetVectorIndex(1, 0)], 3);
  double tau3 = sumP[GetVectorIndex(4, 0)] / pow(sumP[GetVectorIndex(1, 0)], 4);
  // double tau4 = sumP[GetVectorIndex(5,0)]/pow(sumP[GetVectorIndex(1,0)],5);
  double weight1 = 1 - tau1;
  double weight2 = 1 - 3 * tau1 + 2 * tau2;
  double weight3 = 1 - 6 * tau1 + 3 * tau1 * tau1 + 8 * tau2 - 6 * tau3;
  // double weight4 = 1 - 10*tau1 + 15*tau1*tau1 + 20*tau2 - 20*tau1*tau2 - 30*tau3 + 24*tau4;
  if (mpar < 1 || sumP[GetVectorIndex(1, 0)] == 0)
    return;
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(0))->FillProfile(centmult, sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)], (fEventWeight == kEventWeight::kUnity) ? 1.0 : sumP[GetVectorIndex(1, 0)], rn);
  if (mpar < 2 || sumP[GetVectorIndex(2, 0)] == 0 || weight1 == 0)
    return;
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(1))->FillProfile(centmult, 1 / weight1 * (sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight1, rn);
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(2))->FillProfile(centmult, 1 / weight1 * (-2 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 2 * tau1 * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight1, rn);
  if (mpar < 3 || sumP[GetVectorIndex(3, 0)] == 0 || weight2 == 0)
    return;
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(3))->FillProfile(centmult, 1 / weight2 * (sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 3 * tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 2 * tau2 * sumP[GetVectorIndex(3, 3)] / sumP[GetVectorIndex(3, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight2, rn);
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(4))->FillProfile(centmult, 1 / weight2 * (-3 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 3 * tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] + 6 * tau1 * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 6 * tau2 * sumP[GetVectorIndex(3, 2)] / sumP[GetVectorIndex(3, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight2, rn);
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(5))->FillProfile(centmult, 1 / weight2 * (3 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 6 * tau1 * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] - 3 * tau1 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 6 * tau2 * sumP[GetVectorIndex(3, 1)] / sumP[GetVectorIndex(3, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight2, rn);
  if (mpar < 4 || sumP[GetVectorIndex(4, 0)] == 0 || weight3 == 0)
    return;
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(6))->FillProfile(centmult, 1 / weight3 * (sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 6 * tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 3 * tau1 * tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] + 8 * tau2 * sumP[GetVectorIndex(3, 3)] / sumP[GetVectorIndex(3, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 6 * tau3 * sumP[GetVectorIndex(4, 4)] / sumP[GetVectorIndex(4, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight3, rn);
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(7))->FillProfile(centmult, 1 / weight3 * (-4 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 12 * tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 12 * tau1 * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 12 * tau1 * tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] - 8 * tau2 * sumP[GetVectorIndex(3, 3)] / sumP[GetVectorIndex(3, 0)] - 24 * tau2 * sumP[GetVectorIndex(3, 2)] / sumP[GetVectorIndex(3, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 24 * tau3 * sumP[GetVectorIndex(4, 3)] / sumP[GetVectorIndex(4, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight3, rn);
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(8))->FillProfile(centmult, 1 / weight3 * (6 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 6 * tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] - 24 * tau1 * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 6 * tau1 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 6 * tau1 * tau1 * sumP[GetVectorIndex(2, 2)] / sumP[GetVectorIndex(2, 0)] + 12 * tau1 * tau1 * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] + 24 * tau2 * sumP[GetVectorIndex(3, 2)] / sumP[GetVectorIndex(3, 0)] + 24 * tau2 * sumP[GetVectorIndex(3, 1)] / sumP[GetVectorIndex(3, 0)] * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 36 * tau3 * sumP[GetVectorIndex(4, 2)] / sumP[GetVectorIndex(4, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight3, rn);
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(9))->FillProfile(centmult, 1 / weight3 * (-4 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 12 * tau1 * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] + 12 * tau1 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] - 12 * tau1 * tau1 * sumP[GetVectorIndex(2, 1)] / sumP[GetVectorIndex(2, 0)] - 24 * tau2 * sumP[GetVectorIndex(3, 1)] / sumP[GetVectorIndex(3, 0)] - 8 * tau2 * sumP[GetVectorIndex(1, 1)] / sumP[GetVectorIndex(1, 0)] + 24 * tau3 * sumP[GetVectorIndex(4, 1)] / sumP[GetVectorIndex(4, 0)]), (fEventWeight == kEventWeight::kUnity) ? 1.0 : weight3, rn);
  return;
}
double FlowPtContainer::OrderedAddition(std::vector<double> vec)
{
  double sum = 0;
  std::sort(vec.begin(), vec.end());
  for (int i = 0; i < vec.size(); i++) {
    sum += vec[i];
  }
  return sum;
}
void FlowPtContainer::RebinMulti(Int_t nbins)
{
  if (fCMTermList) {
    for (Int_t i = 0; i < fCMTermList->GetEntries(); i++)
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->RebinMulti(nbins);
  }
  if (fCorrList) {
    for (Int_t i = 0; i < fCorrList->GetEntries(); i++)
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->RebinMulti(nbins);
  }
  if (fCovList) {
    for (Int_t i = 0; i < fCovList->GetEntries(); i++)
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->RebinMulti(nbins);
  }
  return;
}
void FlowPtContainer::RebinMulti(Int_t nbins, Double_t* binedges)
{
  if (fCMTermList) {
    for (Int_t i = 0; i < fCMTermList->GetEntries(); i++)
      dynamic_cast<BootstrapProfile*>(fCMTermList->At(i))->RebinMulti(nbins, binedges);
  }
  if (fCorrList) {
    for (Int_t i = 0; i < fCorrList->GetEntries(); i++)
      dynamic_cast<BootstrapProfile*>(fCorrList->At(i))->RebinMulti(nbins, binedges);
  }
  if (fCovList) {
    for (Int_t i = 0; i < fCovList->GetEntries(); i++)
      dynamic_cast<BootstrapProfile*>(fCovList->At(i))->RebinMulti(nbins, binedges);
  }
  return;
}
TH1* FlowPtContainer::getCorrHist(int ind, int m)
{
  return dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(Form("corr_%ipar", m)))->getHist(ind);
}
TH1* FlowPtContainer::getCentralMomentHist(int ind, int m)
{
  if (!fCentralMomentList)
    CreateCentralMomentList();
  if (!fCentralMomentList)
    return 0;
  if (ind + 1 < fCentralMomentList->GetEntries())
    return dynamic_cast<TH1*>(fCentralMomentList->FindObject(Form("cm%i_%i", m + 1, ind)));
  return 0;
}
void FlowPtContainer::CreateCentralMomentList()
{
  if (fCentralMomentList)
    delete fCentralMomentList;
  fCentralMomentList = new TList();
  fCentralMomentList->SetOwner();
  for (auto m(1); m <= 4; ++m) {
    for (int i = -1; i < reinterpret_cast<BootstrapProfile*>(fCMTermList->At(0))->getNSubs(); ++i) {
      TH1* hMpt = reinterpret_cast<BootstrapProfile*>(fCMTermList->At(0))->getHist(i);
      std::vector<TH1*> hTs;
      for (int j = 0; j < m; ++j) {
        dynamic_cast<BootstrapProfile*>(fCMTermList->FindObject(Form("cm%i_Mpt%i", m, j)))->SetErrorOption("g");
        hTs.push_back(reinterpret_cast<BootstrapProfile*>(fCMTermList->FindObject(Form("cm%i_Mpt%i", m, j)))->getHist(i));
      }
      CalculateCentralMomentHists(hTs, i, m, hMpt);
    }
  }
  return;
}
void FlowPtContainer::CalculateCentralMomentHists(std::vector<TH1*> inh, int ind, int m, TH1* hMpt)
{
  TH1* reth = reinterpret_cast<TH1*>(inh[0]->Clone(Form("cm%i_%i", m, ind)));
  for (auto i(1); i < m; ++i) {
    TH1* mptPow = raiseHistToPower(hMpt, i);
    inh[i]->Multiply(mptPow);
    reth->Add(inh[i]);
  }
  TH1* mptLast = raiseHistToPower(hMpt, m);
  reth->Add(mptLast, (m % 2) ? (-1) : 1);
  fCentralMomentList->Add(reth);
  return;
}
TH1* FlowPtContainer::getCumulantHist(int ind, int m)
{
  if (!fCumulantList)
    CreateCumulantList();
  if (!fCumulantList)
    return 0;
  if (ind + 1 < fCumulantList->GetEntries())
    return reinterpret_cast<TH1*>(fCumulantList->At((ind + 1) * mpar + m - 1));
}
void FlowPtContainer::CreateCumulantList()
{
  if (fCumulantList)
    delete fCumulantList;
  fCumulantList = new TList();
  fCumulantList->SetOwner();
  //((BootstrapProfile*)fCorrList->At(0))->PresetWeights((BootstrapProfile*)fCorrList->At(mpar-1));
  for (int i = -1; i < reinterpret_cast<BootstrapProfile*>(fCorrList->At(0))->getNSubs(); ++i) {
    std::vector<TH1*> hTs;
    for (int j = 0; j < mpar; ++j) {
      dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(Form("corr_%ipar", j + 1)))->SetErrorOption("g");
      hTs.push_back(reinterpret_cast<BootstrapProfile*>(fCorrList->FindObject(Form("corr_%ipar", j + 1)))->getHist(i));
    }
    CalculateCumulantHists(hTs, i);
  }
  //((BootstrapProfile*)fCorrList->At(0))->PresetWeights(0);
  return;
}
void FlowPtContainer::CalculateCumulantHists(std::vector<TH1*> inh, int ind)
{
  auto binomial = [&](const int n, const int m) { assert(n >= m); return fFactorial[n]/(fFactorial[m]*fFactorial[n-m]); };
  for (int m = 1; m <= mpar; ++m) {
    TH1* reth = dynamic_cast<TH1*>(inh[m - 1]->Clone(Form("reth%i_%i", m, ind)));
    // TH1* hWeights = (TH1*)inh[m-1]->Clone(Form("hWeights%i_%i",m,ind));
    for (int k = 1; k < m; ++k) {
      TH1* corrh = dynamic_cast<TH1*>(inh[m - k - 1]->Clone(Form("hcorr%i%i_%i", m, k, ind)));
      corrh->Multiply(dynamic_cast<TH1*>(fCumulantList->At((ind + 1) * mpar + k - 1)));
      corrh->Scale(binomial(m - 1, k - 1));
      reth->Add(corrh, -1);
      delete corrh;
    }
    // for(int i=1;i<=hWeights->GetNbinsX();++i) reth->SetBinError(i,hWeights->GetBinError(i));
    // delete hWeights;
    fCumulantList->Add(dynamic_cast<TH1*>(reth->Clone(Form("kappa%i_%i", m, ind))));
  }
  return;
}
Long64_t FlowPtContainer::Merge(TCollection* collist)
{
  if (!fCorrList || !fCMTermList)
    return 0;
  Long64_t nmerged = 0;
  TIter all_PTC(collist);
  FlowPtContainer* l_PTC = 0;
  while ((l_PTC = dynamic_cast<FlowPtContainer*>(all_PTC()))) {
    TList* t_CMTerm = l_PTC->fCMTermList;
    TList* t_Corr = l_PTC->fCorrList;
    TList* t_Cum = l_PTC->fCumulantList;
    TList* t_CM = l_PTC->fCentralMomentList;
    if (t_CMTerm) {
      if (!fCMTermList)
        fCMTermList = dynamic_cast<TList*>(t_CMTerm->Clone());
      else
        MergeBSLists(fCMTermList, t_CMTerm);
      nmerged++;
    }
    if (t_Corr) {
      if (!fCorrList)
        fCorrList = dynamic_cast<TList*>(t_Corr->Clone());
      else
        MergeBSLists(fCorrList, t_Corr);
    }
    if (t_Cum) {
      if (!fCumulantList)
        fCumulantList = dynamic_cast<TList*>(t_Cum->Clone());
      else
        MergeBSLists(fCumulantList, t_Cum);
    }
    if (t_CM) {
      if (!fCentralMomentList)
        fCentralMomentList = dynamic_cast<TList*>(t_CM->Clone());
      else
        MergeBSLists(fCentralMomentList, t_CM);
    }
  }
  return nmerged;
}
void FlowPtContainer::MergeBSLists(TList* source, TList* target)
{
  if (source->GetEntries() != target->GetEntries()) {
    printf("Number in lists to be merged are not the same, skipping...\n");
    return;
  }
  for (Int_t i = 0; i < source->GetEntries(); i++) {
    BootstrapProfile* l_obj = dynamic_cast<BootstrapProfile*>(source->At(i));
    BootstrapProfile* t_obj = dynamic_cast<BootstrapProfile*>(target->At(i));
    l_obj->MergeBS(t_obj);
  }
}
TH1* FlowPtContainer::raiseHistToPower(TH1* inh, double p)
{
  TH1D* reth = dynamic_cast<TH1D*>(inh->Clone("reth"));
  reth->SetName(Form("power%.2f_%s", p, inh->GetName()));
  for (int i = 1; i <= inh->GetNbinsX(); i++) {
    if (inh->GetBinContent(i) >= 0 || std::floor(p) == p) {
      reth->SetBinContent(i, pow(inh->GetBinContent(i), p));
      reth->SetBinError(i, p * pow(reth->GetBinContent(i), p - 1) * inh->GetBinError(i));
    } else {
      reth->SetBinContent(i, -999);
      reth->SetBinError(i, 0.000000001);
    }
  }
  return reth;
}
