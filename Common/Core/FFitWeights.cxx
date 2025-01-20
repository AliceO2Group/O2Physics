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

/// \file FFitWeights.cxx
/// \brief Implementation file for FFitWeights.h, see the header for more information
///
/// \author Joachim C. K. B. Hansen

#include "FFitWeights.h"

#include <string>
#include <cstdio>
#include <vector>

#include <TSpline.h>

ClassImp(FFitWeights)

  FFitWeights::FFitWeights() : TNamed("", ""),
                               fW_data{nullptr},
                               CentBin{100},
                               qAxis{nullptr},
                               nResolution{3000},
                               qnTYPE{0}
{
}

FFitWeights::FFitWeights(const char* name) : TNamed(name, name),
                                             fW_data{nullptr},
                                             CentBin{100},
                                             qAxis{nullptr},
                                             nResolution{3000},
                                             qnTYPE{0} {}

FFitWeights::~FFitWeights()
{
  delete fW_data;
  if (qAxis)
    delete qAxis;
};

void FFitWeights::init()
{
  fW_data = new TObjArray();
  fW_data->SetName("FFitWeights_Data");
  fW_data->SetOwner(kTRUE);

  if (!qAxis)
    this->setBinAxis(500, 0, 25);
  for (const auto& qn : qnTYPE) {
    fW_data->Add(new TH2D(this->getQName(qn.first, qn.second.c_str()), this->getAxisName(qn.first, qn.second.c_str()), CentBin, 0, CentBin, qAxis->GetNbins(), qAxis->GetXmin(), qAxis->GetXmax()));
  }
};

void FFitWeights::fillWeights(float centrality, float qn, int nh, const char* pf)
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return;

  TH2D* th2 = reinterpret_cast<TH2D*>(tar->FindObject(this->getQName(nh, pf)));
  if (!th2) {
    tar->Add(new TH2D(this->getQName(nh, pf), this->getAxisName(nh, pf), CentBin, 0, CentBin, qAxis->GetNbins(), qAxis->GetXmin(), qAxis->GetXmax()));
    th2 = reinterpret_cast<TH2D*>(tar->At(tar->GetEntries() - 1));
  }
  th2->Fill(centrality, qn);
};

Long64_t FFitWeights::Merge(TCollection* collist)
{
  Long64_t nmerged = 0;
  if (!fW_data) {
    fW_data = new TObjArray();
    fW_data->SetName("FFitWeights_Data");
    fW_data->SetOwner(kTRUE);
  }
  FFitWeights* l_w = 0;
  TIter all_w(collist);
  while ((l_w = (reinterpret_cast<FFitWeights*>(all_w())))) {
    addArray(fW_data, l_w->getDataArray());
    nmerged++;
  }
  return nmerged;
};
void FFitWeights::addArray(TObjArray* targ, TObjArray* sour)
{
  if (!sour) {
    printf("Source array does not exist!\n");
    return;
  }
  for (int i = 0; i < sour->GetEntries(); i++) {
    TH2D* sourh = reinterpret_cast<TH2D*>(sour->At(i));
    TH2D* targh = reinterpret_cast<TH2D*>(targ->FindObject(sourh->GetName()));
    if (!targh) {
      targh = reinterpret_cast<TH2D*>(sourh->Clone(sourh->GetName()));
      targh->SetDirectory(0);
      targ->Add(targh);
    } else {
      targh->Add(sourh);
    }
  }
};

void FFitWeights::qSelection(std::vector<int> nhv, std::vector<std::string> stv) /* only execute OFFLINE */
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return;

  for (const auto& pf : stv) {
    for (const auto& nh : nhv) {
      TH2D* th2{reinterpret_cast<TH2D*>(tar->FindObject(this->getQName(nh, pf.c_str())))};
      if (!th2) {
        printf("qh not found!\n");
        return;
      }

      TH1D* tmp{nullptr};
      TGraph* tmpgr{nullptr};
      // TSpline3* spline = nullptr;
      for (int iSP{0}; iSP < 90; iSP++) {
        tmp = th2->ProjectionY(Form("q%i_%i_%i", nh, iSP, iSP + 1), iSP + 1, iSP + 1);
        std::vector<double> xq(nResolution);
        std::vector<double> yq(nResolution);
        for (int i{0}; i < nResolution; i++)
          xq[i] = static_cast<double>(i + 1) / static_cast<double>(nResolution);
        tmp->GetQuantiles(nResolution, yq.data(), xq.data());
        tmpgr = new TGraph(nResolution, yq.data(), xq.data());
        tmpgr->SetName(Form("sp_q%i%s_%i", nh, pf.c_str(), iSP));
        // spline = new TSpline3(Form("sp_q%i%s_%i", nh, pf.c_str(), iSP), tmpgr);
        // spline->SetName(Form("sp_q%i%s_%i", nh, pf.c_str(), iSP));
        fW_data->Add(tmpgr);
      }
    }
  }
};

float FFitWeights::eval(float centr, const float& dqn, const int nh, const char* pf)
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return -1;

  int isp{static_cast<int>(centr)};
  if (isp < 0 || isp > 90) {
    return -1;
  }

  TGraph* spline{nullptr};
  spline = reinterpret_cast<TGraph*>(tar->FindObject(Form("sp_q%i%s_%i", nh, pf, isp)));
  if (!spline) {
    return -1;
  }

  float qn_val{static_cast<float>(100. * spline->Eval(dqn))};
  if (qn_val < 0 || qn_val > 100.05) {
    return -1;
  }

  return qn_val;
};
