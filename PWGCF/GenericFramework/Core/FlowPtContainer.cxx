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
                                     fUseCentralMoments(true),
                                     fUseGap(false),
                                     sumP(),
                                     corrNum(),
                                     corrDen(),
                                     cmVal(),
                                     cmDen() {}
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
                                                     fUseCentralMoments(true),
                                                     fUseGap(false),
                                                     sumP(),
                                                     corrNum(),
                                                     corrDen(),
                                                     cmVal(),
                                                     cmDen() {}
FlowPtContainer::FlowPtContainer(const char* name, const char* title) : TNamed(name, title),
                                                                        fCMTermList(0),
                                                                        fCorrList(0),
                                                                        fCovList(0),
                                                                        fCumulantList(0),
                                                                        fCentralMomentList(0),
                                                                        mpar(0),
                                                                        fillCounter(0),
                                                                        fEventWeight(kEventWeight::kUnity),
                                                                        fUseCentralMoments(true),
                                                                        fUseGap(false),
                                                                        sumP(),
                                                                        corrNum(),
                                                                        corrDen(),
                                                                        cmVal(),
                                                                        cmDen() {}
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
  for (int m = 0; m < mpar; ++m) {
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("mpt%i", m + 1), nMultiBins, &multiBins[0]));
  }
  for (int m = 0; m < 4; ++m) {
    for (int i = 0; i <= m; ++i) {
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nMultiBins, &multiBins[0]));
    }
  }
  if (fUseGap) {
    for (int i = 0; i < configs.GetSize(); ++i) {
      for (auto m(1); m <= mpar; ++m) {
        if (!(configs.GetpTCorrMasks()[i] & (1 << (m - 1))))
          continue;
        if (fUseCentralMoments) {
          for (auto j = 0; j <= m; ++j) {
            fCovList->Add(new BootstrapProfile(Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), nMultiBins, &multiBins[0]));
          }
        } else {
          fCovList->Add(new BootstrapProfile(Form("%spt%i", configs.GetHeads()[i].c_str(), m), Form("%spt%i", configs.GetHeads()[i].c_str(), m), nMultiBins, &multiBins[0]));
        }
      }
    }
  } else {
    if (fUseCentralMoments) {
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt0", "ChFull24pt2_Mpt0", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt1", "ChFull24pt2_Mpt1", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt2", "ChFull24pt2_Mpt2", nMultiBins, &multiBins[0]));

      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt0", "ChFull24pt1_Mpt0", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt1", "ChFull24pt1_Mpt1", nMultiBins, &multiBins[0]));

      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt0", "ChFull22pt2_Mpt0", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt1", "ChFull22pt2_Mpt1", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt2", "ChFull22pt2_Mpt2", nMultiBins, &multiBins[0]));

      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt0", "ChFull22pt1_Mpt0", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt1", "ChFull22pt1_Mpt1", nMultiBins, &multiBins[0]));
    } else {
      fCovList->Add(new BootstrapProfile("ChFull24pt2", "ChFull24pt2", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull24pt1", "ChFull24pt1", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull22pt2", "ChFull22pt2", nMultiBins, &multiBins[0]));
      fCovList->Add(new BootstrapProfile("ChFull22pt1", "ChFull22pt1", nMultiBins, &multiBins[0]));
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
  for (int m = 0; m < mpar; ++m) {
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("mpt%i", m + 1), nbinsx, xbins));
  }
  for (int m = 0; m < 4; ++m) {
    for (int i = 0; i <= m; ++i) {
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nbinsx, xbins));
    }
  }
  if (fUseGap) {
    for (int i = 0; i < configs.GetSize(); ++i) {
      for (auto m(1); m <= mpar; ++m) {
        if (!(configs.GetpTCorrMasks()[i] & (1 << (m - 1))))
          continue;
        if (fUseCentralMoments) {
          for (auto j = 0; j <= m; ++j) {
            fCovList->Add(new BootstrapProfile(Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), nbinsx, xbins));
          }
        } else {
          fCovList->Add(new BootstrapProfile(Form("%spt%i", configs.GetHeads()[i].c_str(), m), Form("%spt%i", configs.GetHeads()[i].c_str(), m), nbinsx, xbins));
        }
      }
    }
  } else {
    if (fUseCentralMoments) {
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt0", "ChFull24pt2_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt1", "ChFull24pt2_Mpt1", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt2", "ChFull24pt2_Mpt2", nbinsx, xbins));

      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt0", "ChFull24pt1_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt1", "ChFull24pt1_Mpt1", nbinsx, xbins));

      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt0", "ChFull22pt2_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt1", "ChFull22pt2_Mpt1", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt2", "ChFull22pt2_Mpt2", nbinsx, xbins));

      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt0", "ChFull22pt1_Mpt0", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt1", "ChFull22pt1_Mpt1", nbinsx, xbins));
    } else {
      fCovList->Add(new BootstrapProfile("ChFull24pt2", "ChFull24pt2", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull24pt1", "ChFull24pt1", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt2", "ChFull22pt2", nbinsx, xbins));
      fCovList->Add(new BootstrapProfile("ChFull22pt1", "ChFull22pt1", nbinsx, xbins));
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
  for (int m = 0; m < mpar; ++m) {
    fCorrList->Add(new BootstrapProfile(Form("mpt%i", m + 1), Form("mpt%i", m + 1), nbinsx, xlow, xhigh));
  }
  for (int m = 0; m < 4; ++m) {
    for (int i = 0; i <= m; ++i) {
      fCMTermList->Add(new BootstrapProfile(Form("cm%i_Mpt%i", m + 1, i), Form("cm%i_Mpt%i", m + 1, i), nbinsx, xlow, xhigh));
    }
  }
  if (fUseGap) {
    for (int i = 0; i < configs.GetSize(); ++i) {
      for (auto m(1); m <= mpar; ++m) {
        if (!(configs.GetpTCorrMasks()[i] & (1 << (m - 1))))
          continue;
        if (fUseCentralMoments) {
          for (auto j = 0; j <= m; ++j) {
            fCovList->Add(new BootstrapProfile(Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), Form("%spt%i_Mpt%i", configs.GetHeads()[i].c_str(), m, j), nbinsx, xlow, xhigh));
          }
        } else {
          fCovList->Add(new BootstrapProfile(Form("%spt%i", configs.GetHeads()[i].c_str(), m), Form("%spt%i", configs.GetHeads()[i].c_str(), m), nbinsx, xlow, xhigh));
        }
      }
    }
  } else {
    if (fUseCentralMoments) {
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt0", "ChFull24pt2_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt1", "ChFull24pt2_Mpt1", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull24pt2_Mpt2", "ChFull24pt2_Mpt2", nbinsx, xlow, xhigh));

      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt0", "ChFull24pt1_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull24pt1_Mpt1", "ChFull24pt1_Mpt1", nbinsx, xlow, xhigh));

      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt0", "ChFull22pt2_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt1", "ChFull22pt2_Mpt1", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt2_Mpt2", "ChFull22pt2_Mpt2", nbinsx, xlow, xhigh));

      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt0", "ChFull22pt1_Mpt0", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt1_Mpt1", "ChFull22pt1_Mpt1", nbinsx, xlow, xhigh));
    } else {
      fCovList->Add(new BootstrapProfile("ChFull24pt2", "ChFull24pt2", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull24pt1", "ChFull24pt1", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt2", "ChFull22pt2", nbinsx, xlow, xhigh));
      fCovList->Add(new BootstrapProfile("ChFull22pt1", "ChFull22pt1", nbinsx, xlow, xhigh));
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
  for (size_t i = 0; i < sumP.size(); ++i) {
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
    if (corrDen[m] != 0) {
      dynamic_cast<BootstrapProfile*>(fCorrList->At(m - 1))->FillProfile(centmult, corrNum[m] / corrDen[m], (fEventWeight == kEventWeight::kUnity) ? 1.0 : corrDen[m], rn);
    }
  }
  return;
}
void FlowPtContainer::FillVnPtCorrProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
{
  if (!mask) {
    return;
  }
  for (auto m(1); m <= mpar; ++m) {
    if (!(mask & (1 << (m - 1)))) {
      continue;
    }
    if (corrDen[m] != 0) {
      dynamic_cast<BootstrapProfile*>(fCovList->At(fillCounter))->FillProfile(centmult, flowval * corrNum[m] / corrDen[m], (fEventWeight == kUnity) ? 1.0 : flowtuples * corrDen[m], rn);
    }
    ++fillCounter;
  }
  return;
}
void FlowPtContainer::FillVnDeltaPtProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
{
  if (!mask) {
    return;
  }
  for (auto m(1); m <= mpar; ++m) {
    if (!(mask & (1 << (m - 1))))
      continue;
    for (auto i = 0; i <= m; ++i) {
      if (cmDen[m] != 0) {
        dynamic_cast<BootstrapProfile*>(fCovList->At(fillCounter))->FillProfile(centmult, flowval * ((i == m) ? cmVal[0] : cmVal[m * (m - 1) / 2 + (m - i)]), (fEventWeight == kUnity) ? 1.0 : flowtuples * cmDen[m], rn);
      }
      ++fillCounter;
    }
  }
  return;
}
void FlowPtContainer::FillVnPtCorrStdProfiles(const double& centmult, const double& rn)
{
  double wAABBCC = getStdAABBCC(warr);
  if (wAABBCC != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(0))->FillProfile(centmult, getStdAABBCC(arr) / wAABBCC, (fEventWeight == kUnity) ? 1.0 : wAABBCC, rn);
  double wAABBC = getStdAABBC(warr);
  if (wAABBC != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(1))->FillProfile(centmult, getStdAABBCC(arr) / wAABBC, (fEventWeight == kUnity) ? 1.0 : wAABBC, rn);
  double wABCC = getStdAABBC(warr);
  if (wABCC != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(2))->FillProfile(centmult, getStdABCC(arr) / wABCC, (fEventWeight == kUnity) ? 1.0 : wABCC, rn);
  double wABC = getStdABC(warr);
  if (wABC != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(3))->FillProfile(centmult, getStdABC(arr) / wABC, (fEventWeight == kUnity) ? 1.0 : wABC, rn);
  return;
}
void FlowPtContainer::FillVnDeltaPtStdProfiles(const double& centmult, const double& rn)
{
  double wAABBCC = getStdAABBCC(warr);
  if (wAABBCC != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(0))->FillProfile(centmult, getStdAABBCC(arr) / wAABBCC, (fEventWeight == kUnity) ? 1.0 : wAABBCC, rn);
  double wAABBCD = getStdAABBCD(warr);
  if (wAABBCD != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(1))->FillProfile(centmult, getStdAABBCD(arr) / wAABBCD, (fEventWeight == kUnity) ? 1.0 : wAABBCD, rn);
  double wAABBDD = getStdAABBDD(warr);
  if (wAABBDD != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(2))->FillProfile(centmult, getStdAABBDD(arr) / wAABBDD, (fEventWeight == kUnity) ? 1.0 : wAABBDD, rn);

  double wAABBC = getStdAABBC(warr);
  if (wAABBC != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(3))->FillProfile(centmult, getStdAABBC(arr) / wAABBC, (fEventWeight == kUnity) ? 1.0 : wAABBC, rn);
  double wAABBD = getStdAABBD(warr);
  if (wAABBD != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(4))->FillProfile(centmult, getStdAABBD(arr) / wAABBD, (fEventWeight == kUnity) ? 1.0 : wAABBD, rn);

  double wABCC = getStdABCC(warr);
  if (wABCC != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(5))->FillProfile(centmult, getStdABCC(arr) / wABCC, (fEventWeight == kUnity) ? 1.0 : wABCC, rn);
  double wABCD = getStdABCD(warr);
  if (wABCD != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(6))->FillProfile(centmult, getStdABCD(arr) / wABCD, (fEventWeight == kUnity) ? 1.0 : wABCD, rn);
  double wABDD = getStdABDD(warr);
  if (wABDD != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(7))->FillProfile(centmult, getStdABDD(arr) / wABDD, (fEventWeight == kUnity) ? 1.0 : wABDD, rn);

  double wABC = getStdABC(warr);
  if (wABC != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(8))->FillProfile(centmult, getStdABC(arr) / wABC, (fEventWeight == kUnity) ? 1.0 : wABC, rn);
  double wABD = getStdABD(warr);
  if (wABD != 0)
    dynamic_cast<BootstrapProfile*>(fCovList->At(9))->FillProfile(centmult, getStdABD(arr) / wABD, (fEventWeight == kUnity) ? 1.0 : wABD, rn);
  return;
}
void FlowPtContainer::FillCMProfiles(const double& centmult, const double& rn)
{
  if (sumP[GetVectorIndex(0, 0)] == 0)
    return;
  // 0th order correlation
  cmDen.push_back(1.);
  cmVal.push_back(1.);

  cmDen.push_back(sumP[GetVectorIndex(1, 0)]);
  cmDen.push_back(sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] - sumP[GetVectorIndex(2, 0)]);
  cmDen.push_back(sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] - 3 * sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 0)] + 2 * sumP[GetVectorIndex(3, 0)]);
  cmDen.push_back(sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] - 6 * sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] + 8 * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(3, 0)] + 3 * sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(2, 0)] - 6 * sumP[GetVectorIndex(4, 0)]);
  if (mpar < 1 || cmDen[1] == 0)
    return;
  cmVal.push_back(sumP[GetVectorIndex(1, 1)] / cmDen[1]);
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(0))->FillProfile(centmult, cmVal[1], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[0], rn);
  if (mpar < 2 || sumP[GetVectorIndex(2, 0)] == 0 || cmDen[2] == 0)
    return;
  cmVal.push_back(1 / cmDen[2] * (sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] - sumP[GetVectorIndex(2, 2)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(1))->FillProfile(centmult, cmVal[2], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[1], rn);
  cmVal.push_back(-2 * 1 / cmDen[2] * (sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 1)] - sumP[GetVectorIndex(2, 1)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(2))->FillProfile(centmult, cmVal[3], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[1], rn);
  if (mpar < 3 || sumP[GetVectorIndex(3, 0)] == 0 || cmDen[3] == 0)
    return;
  cmVal.push_back(1 / cmDen[3] * (sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] - 3 * sumP[GetVectorIndex(2, 2)] * sumP[GetVectorIndex(1, 1)] + 2 * sumP[GetVectorIndex(3, 3)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(3))->FillProfile(centmult, cmVal[4], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[2], rn);
  cmVal.push_back(-3 * 1 / cmDen[3] * (sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 0)] - 2 * sumP[GetVectorIndex(2, 1)] * sumP[GetVectorIndex(1, 1)] + 2 * sumP[GetVectorIndex(3, 2)] - sumP[GetVectorIndex(2, 2)] * sumP[GetVectorIndex(1, 0)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(4))->FillProfile(centmult, cmVal[5], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[2], rn);
  cmVal.push_back(3 * 1 / cmDen[3] * (sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] - 2 * sumP[GetVectorIndex(2, 1)] * sumP[GetVectorIndex(1, 0)] + 2 * sumP[GetVectorIndex(3, 1)] - sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(2, 0)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(5))->FillProfile(centmult, cmVal[6], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[2], rn);
  if (mpar < 4 || sumP[GetVectorIndex(4, 0)] == 0 || cmDen[4] == 0)
    return;
  cmVal.push_back(1 / cmDen[4] * (sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] - 6 * sumP[GetVectorIndex(2, 2)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] + 3 * sumP[GetVectorIndex(2, 2)] * sumP[GetVectorIndex(2, 2)] + 8 * sumP[GetVectorIndex(3, 3)] * sumP[GetVectorIndex(1, 1)] - 6 * sumP[GetVectorIndex(4, 4)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(6))->FillProfile(centmult, cmVal[7], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[3], rn);
  cmVal.push_back(-4 * 1 / cmDen[4] * (sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 0)] - 3 * sumP[GetVectorIndex(2, 2)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 0)] - 3 * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(2, 1)] + 3 * sumP[GetVectorIndex(2, 2)] * sumP[GetVectorIndex(2, 1)] + 6 * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(3, 2)] - 6 * sumP[GetVectorIndex(4, 3)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(7))->FillProfile(centmult, cmVal[8], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[3], rn);
  cmVal.push_back(6 * 1 / cmDen[4] * (sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] - sumP[GetVectorIndex(2, 2)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] - sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 1)] + sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(2, 2)] - 4 * sumP[GetVectorIndex(2, 1)] * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 0)] + 4 * sumP[GetVectorIndex(3, 2)] * sumP[GetVectorIndex(1, 0)] + 4 * sumP[GetVectorIndex(3, 1)] * sumP[GetVectorIndex(1, 1)] + 2 * sumP[GetVectorIndex(2, 1)] * sumP[GetVectorIndex(2, 1)] - 6 * sumP[GetVectorIndex(4, 2)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(8))->FillProfile(centmult, cmVal[9], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[3], rn);
  cmVal.push_back(-4 * 1 / cmDen[4] * (sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] - 3 * sumP[GetVectorIndex(2, 1)] * sumP[GetVectorIndex(1, 0)] * sumP[GetVectorIndex(1, 0)] - 3 * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(2, 0)] * sumP[GetVectorIndex(1, 0)] + 3 * sumP[GetVectorIndex(2, 1)] * sumP[GetVectorIndex(2, 0)] + 2 * sumP[GetVectorIndex(1, 1)] * sumP[GetVectorIndex(3, 0)] + 6 * sumP[GetVectorIndex(3, 1)] * sumP[GetVectorIndex(1, 0)] - 6 * sumP[GetVectorIndex(4, 1)]));
  dynamic_cast<BootstrapProfile*>(fCMTermList->At(9))->FillProfile(centmult, cmVal[10], (fEventWeight == kEventWeight::kUnity) ? 1.0 : cmDen[3], rn);
  return;
}
void FlowPtContainer::FillArray(FillType a, FillType b, double c, double d)
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          if (std::holds_alternative<std::complex<double>>(a) && std::holds_alternative<std::complex<double>>(b)) {
            arr[i][j][k][l] += pow(std::get<0>(a), i) * pow(std::get<0>(b), j) * pow(c, k) * pow(d, l);
          } else if (std::holds_alternative<double>(a) && std::holds_alternative<double>(b)) {
            warr[i][j][k][l] += pow(std::get<1>(a), i) * pow(std::get<1>(b), j) * pow(c, k) * pow(d, l);
          } else {
            LOGF(error, "FillType variant should hold same type for a and b during single function c");
          }
        }
      }
    }
  }
  return;
}
void FlowPtContainer::ClearArray()
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          arr[i][j][k][l] = {0.0, 0.0};
          warr[i][j][k][l] = 0.0;
        }
      }
    }
  }

  return;
}
template <typename T>
double FlowPtContainer::getStdAABBCC(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> c = inarr[0][0][1][0];
  std::complex<double> aa = inarr[2][0][0][0];
  std::complex<double> bb = inarr[0][2][0][0];
  std::complex<double> cc = inarr[0][0][2][0];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ac = inarr[1][0][1][0];
  std::complex<double> bc = inarr[0][1][1][0];
  std::complex<double> aab = inarr[2][1][0][0];
  std::complex<double> aac = inarr[2][0][1][0];
  std::complex<double> abb = inarr[1][2][0][0];
  std::complex<double> acc = inarr[1][0][2][0];
  std::complex<double> abc = inarr[1][1][1][0];
  std::complex<double> bbc = inarr[0][2][1][0];
  std::complex<double> bcc = inarr[0][1][2][0];
  std::complex<double> aabb = inarr[2][2][0][0];
  std::complex<double> aacc = inarr[2][0][2][0];
  std::complex<double> aabc = inarr[2][1][1][0];
  std::complex<double> abbc = inarr[1][2][1][0];
  std::complex<double> abcc = inarr[1][1][2][0];
  std::complex<double> bbcc = inarr[0][2][2][0];
  std::complex<double> aabbc = inarr[2][2][1][0];
  std::complex<double> aabcc = inarr[2][1][2][0];
  std::complex<double> abbcc = inarr[0][0][0][0];
  std::complex<double> aabbcc = inarr[2][2][2][0];
  return (a * a * b * b * c * c - aa * b * b * c * c - a * a * bb * c * c - a * a * b * b * cc - 4. * a * ab * b * c * c -
          4. * a * ac * b * b * c - 4. * a * a * b * bc * c + 4. * aab * b * c * c + 4. * aac * b * b * c +
          4. * a * abb * c * c + 4. * a * acc * b * b + 4. * a * a * bbc * c + 4. * a * a * b * bcc +
          16. * a * abc * b * c + aa * bb * c * c + aa * b * b * cc + a * a * bb * cc + 2. * ab * ab * c * c +
          2. * ac * ac * b * b + 2. * a * a * bc * bc + 4. * aa * b * bc * c + 4. * a * ac * bb * c +
          4. * a * ab * b * cc + 8. * ab * ac * b * c + 8. * a * ab * bc * c + 8. * a * ac * b * bc - 6. * aabb * c * c -
          24. * aabc * b * c - 6. * aacc * b * b - 24. * abbc * a * c - 24. * abcc * a * b - 6. * bbcc * a * a -
          8. * aab * bc * c - 8. * aac * b * bc - 4. * aac * bb * c - 4. * aab * b * cc - 8. * abb * ac * c -
          4. * abb * a * cc - 8. * acc * ab * b - 4. * acc * a * bb - 8. * bbc * a * ac - 4. * bbc * aa * c -
          8. * bcc * a * ab - 4. * bcc * aa * b - 16. * abc * ab * c - 16. * abc * ac * b - 16. * abc * a * bc -
          aa * bb * cc - 2. * ab * ab * cc - 2. * ac * ac * bb - 2. * bc * bc * aa - 8. * ab * ac * bc +
          48. * aabbc * c + 48. * aabcc * b + 48. * abbcc * a + 6. * aabb * cc + 6. * aacc * bb +
          6. * bbcc * aa + 24. * aabc * bc + 24. * abbc * ac + 24. * abcc * ab + 8. * aab * bcc +
          8. * aac * bbc + 8. * abb * acc + 16. * abc * abc - 120. * aabbcc)
    .real();
}
template <typename T>
double FlowPtContainer::getStdAABBCD(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> c = inarr[0][0][1][0];
  std::complex<double> d = inarr[0][0][0][1];
  std::complex<double> aa = inarr[2][0][0][0];
  std::complex<double> bb = inarr[0][2][0][0];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ac = inarr[1][0][1][0];
  std::complex<double> ad = inarr[1][0][0][1];
  std::complex<double> bc = inarr[0][1][1][0];
  std::complex<double> bd = inarr[0][1][0][1];
  std::complex<double> cd = inarr[0][0][1][1];
  std::complex<double> aab = inarr[2][1][0][0];
  std::complex<double> aac = inarr[2][0][1][0];
  std::complex<double> aad = inarr[2][0][0][1];
  std::complex<double> abb = inarr[1][2][0][0];
  std::complex<double> abc = inarr[1][1][1][0];
  std::complex<double> abd = inarr[1][1][0][1];
  std::complex<double> acd = inarr[1][0][1][1];
  std::complex<double> bbc = inarr[0][2][1][0];
  std::complex<double> bbd = inarr[0][2][0][1];
  std::complex<double> bcd = inarr[0][1][1][1];
  std::complex<double> aabb = inarr[2][2][0][0];
  std::complex<double> aabc = inarr[2][1][1][0];
  std::complex<double> aabd = inarr[2][1][0][1];
  std::complex<double> aacd = inarr[2][0][1][1];
  std::complex<double> abbc = inarr[1][2][1][0];
  std::complex<double> abbd = inarr[1][2][0][1];
  std::complex<double> abcd = inarr[0][1][1][1];
  std::complex<double> bbcd = inarr[0][2][1][1];
  std::complex<double> aabbc = inarr[2][2][1][0];
  std::complex<double> aabbd = inarr[2][2][0][1];
  std::complex<double> aabcd = inarr[2][1][1][1];
  std::complex<double> abbcd = inarr[1][2][1][1];
  std::complex<double> aabbcd = inarr[2][2][1][1];
  return (-120. * aabbcd + 48. * a * abbcd + 24. * ab * abcd + 16. * abc * abd + 12. * abbd * ac +
          8. * abb * acd + 12. * abbc * ad + 48. * aabcd * b - 24. * a * abcd * b - 8. * abd * ac * b -
          8. * ab * acd * b - 8. * abc * ad * b - 6. * aacd * b * b + 4. * a * acd * b * b + 2. * ac * ad * b * b +
          6. * aacd * bb - 4. * a * acd * bb - 2. * ac * ad * bb + 4. * aad * bbc - 4. * a * ad * bbc -
          6. * a * a * bbcd + 6. * aa * bbcd + 4. * aac * bbd - 4. * a * ac * bbd + 12. * aabd * bc -
          8. * a * abd * bc - 4. * ab * ad * bc - 4. * aad * b * bc + 4. * a * ad * b * bc + 8. * aab * bcd -
          8. * a * ab * bcd + 4. * a * a * b * bcd - 4. * aa * b * bcd + 12. * aabc * bd - 8. * a * abc * bd -
          4. * ab * ac * bd - 4. * aac * b * bd + 4. * a * ac * b * bd + 2. * a * a * bc * bd - 2. * aa * bc * bd +
          24. * aabbd * c - 12. * a * abbd * c - 8. * ab * abd * c - 4. * abb * ad * c - 12. * aabd * b * c +
          8. * a * abd * b * c + 4. * ab * ad * b * c + 2. * aad * b * b * c - 2. * a * ad * b * b * c -
          2. * aad * bb * c + 2. * a * ad * bb * c + 2. * a * a * bbd * c - 2. * aa * bbd * c - 4. * aab * bd * c +
          4. * a * ab * bd * c - 2. * a * a * b * bd * c + 2. * aa * b * bd * c + 6. * aabb * cd - 2. * ab * ab * cd -
          4. * a * abb * cd - 4. * aab * b * cd + 4. * a * ab * b * cd - a * a * b * b * cd + aa * b * b * cd +
          a * a * bb * cd - aa * bb * cd + 24. * aabbc * d - 12. * a * abbc * d - 8. * ab * abc * d -
          4. * abb * ac * d - 12. * aabc * b * d + 8. * a * abc * b * d + 4. * ab * ac * b * d + 2. * aac * b * b * d -
          2. * a * ac * b * b * d - 2. * aac * bb * d + 2. * a * ac * bb * d + 2. * a * a * bbc * d - 2. * aa * bbc * d -
          4. * aab * bc * d + 4. * a * ab * bc * d - 2. * a * a * b * bc * d + 2. * aa * b * bc * d - 6. * aabb * c * d +
          2. * ab * ab * c * d + 4. * a * abb * c * d + 4. * aab * b * c * d - 4. * a * ab * b * c * d +
          a * a * b * b * c * d - aa * b * b * c * d - a * a * bb * c * d + aa * bb * c * d)
    .real();
}
template <typename T>
double FlowPtContainer::getStdAABBDD(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> d = inarr[0][0][1][1];
  std::complex<double> aa = inarr[2][0][0][0];
  std::complex<double> bb = inarr[0][2][0][0];
  std::complex<double> dd = inarr[0][0][0][2];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ad = inarr[1][0][0][1];
  std::complex<double> bd = inarr[0][1][0][1];
  std::complex<double> aab = inarr[2][1][0][0];
  std::complex<double> aad = inarr[2][0][0][1];
  std::complex<double> abb = inarr[1][2][0][0];
  std::complex<double> add = inarr[1][0][0][2];
  std::complex<double> abd = inarr[1][1][0][1];
  std::complex<double> bbd = inarr[0][2][0][1];
  std::complex<double> bdd = inarr[0][1][0][2];
  std::complex<double> aabb = inarr[2][2][0][0];
  std::complex<double> aadd = inarr[2][0][0][2];
  std::complex<double> aabd = inarr[2][1][0][1];
  std::complex<double> abbd = inarr[1][2][0][1];
  std::complex<double> abdd = inarr[1][1][0][2];
  std::complex<double> bbdd = inarr[0][2][0][2];
  std::complex<double> aabbd = inarr[2][2][0][1];
  std::complex<double> aabdd = inarr[2][1][0][2];
  std::complex<double> abbdd = inarr[0][0][0][2];
  std::complex<double> aabbdd = inarr[2][2][0][2];
  return (-120. * aabbdd + 48. * a * abbdd + 16. * abd * abd + 24. * ab * abdd + 24. * abbd * ad +
          8. * abb * add + 48. * aabdd * b - 24. * a * abdd * b - 16. * abd * ad * b - 8. * ab * add * b -
          6. * aadd * b * b + 2. * ad * ad * b * b + 4. * a * add * b * b + 6. * aadd * bb - 2. * ad * ad * bb -
          4. * a * add * bb + 8. * aad * bbd - 8. * a * ad * bbd - 6. * a * a * bbdd + 6. * aa * bbdd +
          24. * aabd * bd - 16. * a * abd * bd - 8. * ab * ad * bd - 8. * aad * b * bd + 8. * a * ad * b * bd +
          2. * a * a * bd * bd - 2. * aa * bd * bd + 8. * aab * bdd - 8. * a * ab * bdd + 4. * a * a * b * bdd -
          4. * aa * b * bdd + 48. * aabbd * d - 24. * a * abbd * d - 16. * ab * abd * d - 8. * abb * ad * d -
          24. * aabd * b * d + 16. * a * abd * b * d + 8. * ab * ad * b * d + 4. * aad * b * b * d -
          4. * a * ad * b * b * d - 4. * aad * bb * d + 4. * a * ad * bb * d + 4. * a * a * bbd * d - 4. * aa * bbd * d -
          8. * aab * bd * d + 8. * a * ab * bd * d - 4. * a * a * b * bd * d + 4. * aa * b * bd * d - 6. * aabb * d * d +
          2. * ab * ab * d * d + 4. * a * abb * d * d + 4. * aab * b * d * d - 4. * a * ab * b * d * d +
          a * a * b * b * d * d - aa * b * b * d * d - a * a * bb * d * d + aa * bb * d * d + 6. * aabb * dd -
          2. * ab * ab * dd - 4. * a * abb * dd - 4. * aab * b * dd + 4. * a * ab * b * dd - a * a * b * b * dd +
          aa * b * b * dd + a * a * bb * dd - aa * bb * dd)
    .real();
}
template <typename T>
double FlowPtContainer::getStdAABBC(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> c = inarr[0][0][1][0];
  std::complex<double> aa = inarr[2][0][0][0];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ac = inarr[1][0][1][0];
  std::complex<double> bb = inarr[0][2][0][0];
  std::complex<double> bc = inarr[0][1][1][0];
  std::complex<double> aab = inarr[2][1][0][0];
  std::complex<double> aac = arr[2][0][1][0];
  std::complex<double> abb = inarr[1][2][0][0];
  std::complex<double> abc = inarr[1][1][1][0];
  std::complex<double> bbc = inarr[0][2][1][0];
  std::complex<double> aabb = inarr[2][2][0][0];
  std::complex<double> aabc = inarr[2][1][1][0];
  std::complex<double> abbc = inarr[1][2][1][0];
  std::complex<double> aabbc = inarr[2][2][1][0];
  return (a * a * b * b * c - aa * b * b * c - a * a * bb * c - 4. * ab * a * b * c - 2. * a * ac * b * b - 2. * a * a * bc * b + 2. * ab * ab * c + 4. * ab * ac * b + 4. * ab * bc * a + 8. * abc * a * b + 4. * aab * b * c + 2. * aac * b * b + 4. * abb * a * c + 2. * bbc * a * a + aa * bb * c + 2. * aa * b * bc + 2. * bb * a * ac - 12. * aabc * b - 12. * abbc * a - 6. * aabb * c - 8. * abc * ab - 2. * bbc * aa - 2. * aac * bb - 4. * aab * bc - 4. * abb * ac + 24. * aabbc).real();
}
template <typename T>
double FlowPtContainer::getStdAABBD(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> d = inarr[0][0][1][0];
  std::complex<double> aa = inarr[2][0][0][0];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ad = inarr[1][0][1][0];
  std::complex<double> bb = inarr[0][2][0][0];
  std::complex<double> bd = inarr[0][1][1][0];
  std::complex<double> aab = inarr[2][1][0][0];
  std::complex<double> aad = arr[2][0][1][0];
  std::complex<double> abb = inarr[1][2][0][0];
  std::complex<double> abd = inarr[1][1][1][0];
  std::complex<double> bbd = inarr[0][2][1][0];
  std::complex<double> aabb = inarr[2][2][0][0];
  std::complex<double> aabd = inarr[2][1][1][0];
  std::complex<double> abbd = inarr[1][2][1][0];
  std::complex<double> aabbd = inarr[2][2][1][0];
  return (a * a * b * b * d - aa * b * b * d - a * a * bb * d - 4. * ab * a * b * d - 2. * a * ad * b * b - 2. * a * a * bd * b + 2. * ab * ab * d + 4. * ab * ad * b + 4. * ab * bd * a + 8. * abd * a * b + 4. * aab * b * d + 2. * aad * b * b + 4. * abb * a * d + 2. * bbd * a * a + aa * bb * d + 2. * aa * b * bd + 2. * bb * a * ad - 12. * aabd * b - 12. * abbd * a - 6. * aabb * d - 8. * abd * ab - 2. * bbd * aa - 2. * aad * bb - 4. * aab * bd - 4. * abb * ad + 24. * aabbd).real();
}
template <typename T>
double FlowPtContainer::getStdABCC(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> c = inarr[0][0][1][0];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ac = inarr[1][0][1][0];
  std::complex<double> bc = inarr[0][1][1][0];
  std::complex<double> cc = inarr[0][0][2][0];
  std::complex<double> abc = inarr[1][1][1][0];
  std::complex<double> acc = inarr[1][0][2][0];
  std::complex<double> bcc = inarr[0][1][2][0];
  std::complex<double> abcc = inarr[1][1][2][0];
  return (a * b * c * c - a * b * cc - 2. * a * bc * c - 2. * ac * b * c - ab * c * c + 2. * acc * b + 2. * a * bcc + 4. * abc * c + ab * cc + 2. * ac * bc - 6. * abcc).real();
}
template <typename T>
double FlowPtContainer::getStdABCD(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> c = inarr[0][0][1][0];
  std::complex<double> d = inarr[0][0][0][1];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ac = inarr[1][0][1][0];
  std::complex<double> ad = inarr[1][0][0][1];
  std::complex<double> bc = inarr[0][1][1][0];
  std::complex<double> bd = inarr[0][1][0][1];
  std::complex<double> cd = inarr[0][0][1][1];
  std::complex<double> abc = inarr[1][1][1][0];
  std::complex<double> abd = inarr[1][1][0][1];
  std::complex<double> acd = inarr[1][0][1][1];
  std::complex<double> bcd = inarr[0][1][1][1];
  std::complex<double> abcd = inarr[1][1][0][1];
  return (-6. * abcd + 2. * acd * b + ad * bc + 2. * a * bcd + ac * bd + 2. * abd * c - ad * b * c -
          a * bd * c + ab * cd - a * b * cd + 2. * abc * d - ac * b * d - a * bc * d - ab * c * d +
          a * b * c * d)
    .real();
}
template <typename T>
double FlowPtContainer::getStdABDD(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> d = inarr[0][0][0][1];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ad = inarr[1][0][0][1];
  std::complex<double> bd = inarr[0][1][0][1];
  std::complex<double> dd = inarr[0][0][0][2];
  std::complex<double> abd = inarr[1][1][0][1];
  std::complex<double> add = inarr[1][0][0][2];
  std::complex<double> bdd = inarr[0][1][0][2];
  std::complex<double> abdd = inarr[1][1][0][2];
  return (a * b * d * d - a * b * dd - 2. * a * bd * d - 2. * ad * b * d - ab * d * d + 2. * add * b + 2. * a * bdd + 4. * abd * d + ab * dd + 2. * ad * bd - 6. * abdd).real();
}
template <typename T>
double FlowPtContainer::getStdABC(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> c = inarr[0][0][1][0];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ac = inarr[1][0][1][0];
  std::complex<double> bc = inarr[0][1][1][0];
  std::complex<double> abc = inarr[1][1][1][0];
  return (a * b * c - ab * c - ac * b - a * bc + 2. * abc).real();
}
template <typename T>
double FlowPtContainer::getStdABD(T& inarr)
{
  std::complex<double> a = inarr[1][0][0][0];
  std::complex<double> b = inarr[0][1][0][0];
  std::complex<double> d = inarr[0][0][0][1];
  std::complex<double> ab = inarr[1][1][0][0];
  std::complex<double> ad = inarr[1][0][0][1];
  std::complex<double> bd = inarr[0][1][0][1];
  std::complex<double> abd = inarr[1][1][0][1];
  return (a * b * d - ab * d - ad * b - a * bd + 2. * abd).real();
}
double FlowPtContainer::OrderedAddition(std::vector<double> vec)
{
  double sum = 0;
  std::sort(vec.begin(), vec.end());
  for (size_t i = 0; i < vec.size(); i++) {
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
  return dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(Form("mpt%i", m + 1)))->getHist(ind);
}
TH1* FlowPtContainer::getCentralMomentHist(int ind, int m)
{
  if (!fCentralMomentList)
    CreateCentralMomentList();
  if (!fCentralMomentList)
    return 0;
  if (ind + 1 < fCentralMomentList->GetEntries())
    return dynamic_cast<TH1*>(fCentralMomentList->FindObject(Form("cm%i_%i", m, ind)));
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
  return 0;
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
      dynamic_cast<BootstrapProfile*>(fCorrList->FindObject(Form("mpt%i", j + 1)))->SetErrorOption("g");
      hTs.push_back(reinterpret_cast<BootstrapProfile*>(fCorrList->FindObject(Form("mpt%i", j + 1)))->getHist(i));
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
