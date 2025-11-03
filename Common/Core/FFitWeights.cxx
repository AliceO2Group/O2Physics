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

#include <TCollection.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TSpline.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdio>
#include <string>
#include <vector>

ClassImp(FFitWeights)

  FFitWeights::FFitWeights() : TNamed("", ""),
                               fW_data{nullptr},
                               centBin{100},
                               qAxis{nullptr},
                               nResolution{3000},
                               ptBin{100},
                               ptAxis{nullptr},
                               qnTYPE{0}
{
}

FFitWeights::FFitWeights(const char* name) : TNamed(name, name),
                                             fW_data{nullptr},
                                             centBin{100},
                                             qAxis{nullptr},
                                             nResolution{3000},
                                             ptBin{100},
                                             ptAxis{nullptr},
                                             qnTYPE{0} {}

FFitWeights::~FFitWeights()
{
  delete fW_data;
  if (qAxis)
    delete qAxis;
  if (ptAxis)
    delete ptAxis;
};

void FFitWeights::init()
{
  fW_data = new TObjArray();
  fW_data->SetName("FFitWeights_Data");
  fW_data->SetOwner(kTRUE);

  if (!qAxis)
    this->setBinAxis(500, 0, 25);
  for (const auto& qn : qnTYPE) {
    fW_data->Add(new TH2D(this->getQName(qn.first, qn.second.c_str()), this->getAxisName(qn.first, qn.second.c_str()), centBin, 0, centBin, qAxis->GetNbins(), qAxis->GetXmin(), qAxis->GetXmax()));
  }

  if (!ptAxis)
    this->setPtAxis(3000, -3, 3);
  // fW_data->Add(new TH2D("hPtWeight", "", centBin, 0, centBin, ptBin, ptAxis->GetXmin(), ptAxis->GetXmax()));
  fW_data->Add(new TProfile("pMeanPt", "", centBin, 0, centBin));
  fW_data->Add(new TH2D("hPtWeight", "", centBin, 0, centBin, ptBin, ptAxis->GetXmin(), ptAxis->GetXmax()));
};

void FFitWeights::fillWeights(float centrality, float qn, int nh, const char* pf)
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return;

  TH2D* th2 = reinterpret_cast<TH2D*>(tar->FindObject(this->getQName(nh, pf)));
  if (!th2) {
    tar->Add(new TH2D(this->getQName(nh, pf), this->getAxisName(nh, pf), centBin, 0, centBin, qAxis->GetNbins(), qAxis->GetXmin(), qAxis->GetXmax()));
    th2 = reinterpret_cast<TH2D*>(tar->At(tar->GetEntries() - 1));
  }
  th2->Fill(centrality, qn);
};
void FFitWeights::fillPt(float centrality, float pt, bool first)
{
  TObjArray* tar{nullptr};
  tar = fW_data;
  if (!tar)
    return;
  if (first) {
    auto tp = reinterpret_cast<TProfile*>(tar->FindObject("pMeanPt"));
    if (!tp) {
      tar->Add(new TProfile("pMeanPt", "", centBin, 0, centBin));
      tp = reinterpret_cast<TProfile*>(tar->At(tar->GetEntries() - 1));
    }
    tp->Fill(centrality, pt);
  } else {
    auto th2 = reinterpret_cast<TH2D*>(tar->FindObject("hPtWeight"));
    if (!th2) {
      tar->Add(new TH2D("hPtWeight", "", centBin, 0, centBin, ptBin, ptAxis->GetXmin(), ptAxis->GetXmax()));
      th2 = reinterpret_cast<TH2D*>(tar->At(tar->GetEntries() - 1));
    }
    th2->Fill(centrality, pt);
  }
};
float FFitWeights::getPtMult(float centrality)
{
  TObjArray* tar{nullptr};
  tar = fW_data;
  if (!tar)
    return -1;

  auto tp = reinterpret_cast<TProfile*>(tar->FindObject("pMeanPt"));
  if (!tp) {
    return -1;
  }
  return tp->GetBinContent(tp->FindBin(centrality));
};

Long64_t FFitWeights::Merge(TCollection* collist)
{
  Long64_t nmerged = 0;
  if (!fW_data) {
    fW_data = new TObjArray();
    fW_data->SetName("FFitWeights_Data");
    fW_data->SetOwner(kTRUE);
  }
  FFitWeights* lW = 0;
  TIter allW(collist);
  while ((lW = (reinterpret_cast<FFitWeights*>(allW())))) {
    addArray(fW_data, lW->getDataArray());
    nmerged++;
  }
  return nmerged;
};
void FFitWeights::addArray(TObjArray* targ, TObjArray* sour)
{
  if (!sour) {
    // printf("Source array does not exist!\n");
    // LOGF(info, "FFitWeights source array does not exist!");
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

void FFitWeights::mptSel()
{
  TObjArray* tar{nullptr};
  tar = fW_data;
  if (!tar)
    return;

  TH2D* th2 = reinterpret_cast<TH2D*>(tar->FindObject("hPtWeight"));
  if (!th2) {
    return;
  }

  TH1D* tmp{nullptr};
  TGraph* tmpgr{nullptr};
  for (int iSP{0}; iSP < NumberSp; iSP++) {
    tmp = th2->ProjectionY(Form("mpt_%i_%i", iSP, iSP + 1), iSP + 1, iSP + 1);
    std::vector<double> xq(nResolution);
    std::vector<double> yq(nResolution);
    for (int i{0}; i < nResolution; i++)
      xq[i] = static_cast<double>(i + 1) / static_cast<double>(nResolution);
    tmp->GetQuantiles(nResolution, yq.data(), xq.data());
    tmpgr = new TGraph(nResolution, yq.data(), xq.data());
    tmpgr->SetName(Form("sp_mpt_%i", iSP));
    fW_data->Add(tmpgr);
  }
}

void FFitWeights::qSelection(const std::vector<int>& nhv, const std::vector<std::string>& stv) /* only execute OFFLINE */
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return;

  for (const auto& pf : stv) {
    for (const auto& nh : nhv) {
      TH2D* th2{reinterpret_cast<TH2D*>(tar->FindObject(this->getQName(nh, pf.c_str())))};
      if (!th2) {
        // printf("qh not found!\n");
        // LOGF(info, "FFitWeights qh not found!");
        return;
      }

      TH1D* tmp{nullptr};
      TGraph* tmpgr{nullptr};
      // TSpline3* spline = nullptr;
      for (int iSP{0}; iSP < NumberSp; iSP++) {
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
  if (!tar) {
    return -1;
  }

  int isp{static_cast<int>(centr)};
  if (isp < 0 || isp > NumberSp) {
    return -1;
  }

  TGraph* spline{nullptr};
  spline = reinterpret_cast<TGraph*>(tar->FindObject(Form("sp_q%i%s_%i", nh, pf, isp)));
  if (!spline) {
    return -1;
  }

  float qnVal{static_cast<float>(100. * spline->Eval(dqn))};
  if (qnVal < 0 || qnVal > MaxTol) {
    return -1;
  }

  return qnVal;
};

float FFitWeights::evalPt(float centr, const float& mpt)
{
  TObjArray* tar{nullptr};
  tar = fW_data;
  if (!tar) {
    return -1;
  }

  int isp{static_cast<int>(centr)};
  if (isp < 0 || isp > NumberSp) {
    return -1;
  }

  TGraph* spline{nullptr};
  spline = reinterpret_cast<TGraph*>(tar->FindObject(Form("sp_mpt_%i", isp)));
  if (!spline) {
    return -1;
  }
  float ptVal{static_cast<float>(100. * spline->Eval(mpt))};
  if (ptVal < 0 || ptVal > MaxTol) {
    return -1;
  }
  return ptVal;
};
