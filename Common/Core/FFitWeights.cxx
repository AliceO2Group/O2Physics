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

#include "Framework/Logger.h"

#include <TCollection.h>
#include <TH1.h>
#include <TNamed.h>
#include <TObjArray.h>
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
                               ptProfCent{nullptr},
                               h2ptCent{nullptr},
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
                                             ptProfCent{nullptr},
                                             h2ptCent{nullptr},
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
  fW_data->Add(new TProfile("pMeanPt", "", centBin, 0, centBin));
  fW_data->Add(new TH2D("hPtWeight", "", centBin, 0, centBin, ptBin, ptAxis->GetXmin(), ptAxis->GetXmax()));

  ptProfCent = reinterpret_cast<TProfile*>(fW_data->FindObject("pMeanPt"));
  h2ptCent = reinterpret_cast<TH2D*>(fW_data->FindObject("hPtWeight"));
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
  if (first) {
    ptProfCent->Fill(centrality, pt);
  } else {
    h2ptCent->Fill(centrality, pt);
  }
};
float FFitWeights::getPtMult(float centrality)
{
  if (!ptProfCent) {
    ptProfCent = reinterpret_cast<TProfile*>(fW_data->FindObject("pMeanPt"));
  }
  return ptProfCent->GetBinContent(ptProfCent->FindBin(centrality));
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
  if (!sour)
    return;

  for (int i = 0; i < sour->GetEntries(); i++) {
    auto* obj = sour->At(i);
    if (!obj)
      continue;

    auto* tObj = targ->FindObject(obj->GetName());
    if (!tObj) {
      auto* clone = static_cast<TObject*>(obj->Clone(obj->GetName()));
      if (auto* h = dynamic_cast<TH1*>(clone))
        h->SetDirectory(0);
      targ->Add(clone);
    } else if (auto* h1 = dynamic_cast<TH1*>(tObj)) {
      if (auto* h2 = dynamic_cast<TH1*>(obj))
        h1->Add(h2);
    }
  }
}

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
        LOGF(info, "FFitWeights qh not found!");
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
        fW_data->Add(tmpgr);
      }
    }
  }
};
float FFitWeights::internalEval(float centr, const float& val, const char* name)
{
  if (!fW_data) {
    return -1;
  }
  int isp = static_cast<int>(centr);
  if (isp < 0 || isp > NumberSp) {
    return -1;
  }

  auto* spline = dynamic_cast<TGraph*>(fW_data->FindObject(Form(name, isp)));
  if (!spline) {
    return -1;
  }

  float perc = 100.f * spline->Eval(val);
  return (perc < 0 || perc > MaxTol) ? -1 : perc;
};

float FFitWeights::eval(float centr, const float& dqn, int nh, const char* pf)
{
  return internalEval(centr, dqn, Form("sp_q%i%s_%%i", nh, pf));
};

float FFitWeights::evalPt(float centr, const float& mpt)
{
  return internalEval(centr, mpt, "sp_mpt_%i");
};
