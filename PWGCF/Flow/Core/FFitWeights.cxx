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

/*
Author: Joachim Hansen
*/

#include "FFitWeights.h"

// using std::complex;
// using std::pair;
// using std::string;
// using std::vector;

ClassImp(FFitWeights)

  FFitWeights::FFitWeights() : TNamed("", ""),
                               fW_data{nullptr},
                               vGain{0},
                               CentBin{100},
                               ChIDBin{220},
                               sAmpl{1000, 5000},
                               sqVec{250, 3500},
                               sqCorVec{250, 250}
{
}

FFitWeights::FFitWeights(const char* name) : TNamed(name, name),
                                             fW_data{nullptr},
                                             vGain{0},
                                             CentBin{100},
                                             ChIDBin{220},
                                             sAmpl{1000, 5000},
                                             sqVec{250, 3500},
                                             sqCorVec{250, 250} {}

FFitWeights::~FFitWeights()
{
  delete fW_data;
};

void FFitWeights::Init()
{
  fW_data = new TObjArray();
  fW_data->SetName("FFitWeights_Data");
  fW_data->SetOwner(kTRUE);

  const char* tnd = "FT0Ampl";
  fW_data->Add(new TH2F(tnd, ";channel;amplitude", ChIDBin, 0, ChIDBin, sAmpl.bins, 0, sAmpl.range));
  fW_data->Add(new TH2F(Form("%sCorr", tnd), ";channel;amplitude", ChIDBin, 0, ChIDBin, sAmpl.bins, 0, sAmpl.range));

  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C", "x", 2), ";Centrality;Qx_{2}", CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C", "y", 2), ";Centrality;Qy_{2}", CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C", "x", 3), ";Centrality;Qx_{3}", CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C", "y", 3), ";Centrality;Qy_{3}", CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));

  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C_Rec", "x", 2), ";Centrality;Qx_{2}", CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C_Rec", "y", 2), ";Centrality;Qy_{2}", CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C_Rec", "x", 3), ";Centrality;Qx_{3}", CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C_Rec", "y", 3), ";Centrality;Qy_{3}", CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));

  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C_RecTot", "x", 2), ";Centrality;Qx_{2}", CentBin, 0, CentBin, sqCorVec.bins, -sqCorVec.range, sqCorVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C_RecTot", "y", 2), ";Centrality;Qy_{2}", CentBin, 0, CentBin, sqCorVec.bins, -sqCorVec.range, sqCorVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C_RecTot", "x", 3), ";Centrality;Qx_{3}", CentBin, 0, CentBin, sqCorVec.bins, -sqCorVec.range, sqCorVec.range));
  fW_data->Add(new TH2F(Form("hQ%s%i_FT0C_RecTot", "y", 3), ";Centrality;Qy_{3}", CentBin, 0, CentBin, sqCorVec.bins, -sqCorVec.range, sqCorVec.range));
};

void FFitWeights::FillFT0(std::size_t iCh, float amplitude, float GainCst)
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return;

  TH2F* th2 = reinterpret_cast<TH2F*>(tar->FindObject("FT0Ampl"));
  if (!th2) {
    tar->Add(new TH2F("FT0Ampl", ";channel;amplitude", ChIDBin, 0, ChIDBin, sAmpl.bins, 0, sAmpl.range));
    th2 = reinterpret_cast<TH2F*>(tar->At(tar->GetEntries() - 1));
  }
  th2->Fill(iCh, amplitude);

  TH2F* th2Cor = reinterpret_cast<TH2F*>(tar->FindObject("FT0AmplCorr"));
  if (!th2Cor) {
    tar->Add(new TH2F("FT0AmplCorr", ";channel;amplitude", ChIDBin, 0, ChIDBin, sAmpl.bins, 0, sAmpl.range));
    th2Cor = reinterpret_cast<TH2F*>(tar->At(tar->GetEntries() - 1));
  }
  th2Cor->Fill(iCh, amplitude / GainCst);
};

void FFitWeights::FillQ(float mult, float vec, int nHarm, const char* coord, const char* qType)
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return;

  TH2F* th2 = reinterpret_cast<TH2F*>(tar->FindObject(Form("hQ%s%i_FT0C%s", coord, nHarm, qType)));
  if (!th2) {
    tar->Add(new TH2F(Form("hQ%s%i_FT0C%s", coord, nHarm, qType), Form(";Centrality;Q%s_{%i}", coord, nHarm), CentBin, 0, CentBin, sqVec.bins, -sqVec.range, sqVec.range));
    th2 = reinterpret_cast<TH2F*>(tar->At(tar->GetEntries() - 1));
  }
  th2->Fill(mult, vec);
};

void FFitWeights::CreateGain()
{
  vGain.clear();

  TH1D* h1;
  if (fW_data->GetEntries() < 1)
    return;

  TH2F* hGain = reinterpret_cast<TH2F*>(fW_data->At(0)->Clone("FT0Ampl"));
  double vMean = hGain->GetMean(2);

  for (int iCh{0}; iCh < hGain->GetNbinsX(); iCh++) {
    h1 = static_cast<TH1D*>(hGain->ProjectionY(Form("proj%i", iCh), iCh, iCh));
    double mean = h1->GetMean();
    double meanErr = h1->GetMeanError();

    double fWeight = mean / vMean;

    if (fWeight > 0) {
      vGain.push_back(fWeight);
    } else {
      vGain.push_back(1.0);
    }

    TObjArray* tar = fW_data;
    if (!tar)
      return;

    fW_data->Add(new TH1F("FT0MultCorr", ";channel", ChIDBin, 0, ChIDBin));
    TH1F* htmp = reinterpret_cast<TH1F*>(tar->FindObject("FT0MultCorr"));
    if (!htmp)
      return;

    htmp->SetBinContent(iCh, mean);
    htmp->SetBinError(iCh, meanErr);
  }

  // note to self if FT0A is implemented this has to be done differently
};

std::vector<float> FFitWeights::GetGain()
{
  return vGain;
};

void FFitWeights::CreateRecenter(const char* xy)
{
  if (fW_data->GetEntries() < 1)
    return;

  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return;

  for (int nHarm{2}; nHarm <= 3; nHarm++) {
    fW_data->Add(new TH1F(Form("havgQ%s%i_FT0C", xy, nHarm), "", CentBin, 0, CentBin));
    TH1F* hRec = reinterpret_cast<TH1F*>(tar->FindObject(Form("havgQ%s%i_FT0C", xy, nHarm)));
    if (!hRec)
      return;

    TH2F* hQ = reinterpret_cast<TH2F*>(fW_data->At(0)->Clone(Form("hQ%s%i_FT0C", xy, nHarm)));
    TH1D* h1;

    for (int i{1}; i < hQ->GetXaxis()->GetNbins() + 1; i++) {
      h1 = static_cast<TH1D*>(hQ->ProjectionX(Form("proj%i_Q%s%is", i, xy, nHarm), i, i));

      double mean = h1->GetMean();
      double meanErr = h1->GetMeanError();
      int binc = hRec->GetXaxis()->GetBinCenter(i);

      hRec->SetBinContent(binc, mean);
      hRec->SetBinError(binc, meanErr);
    }
  }
};

void FFitWeights::CreateRMS(const char* xy)
{
  if (fW_data->GetEntries() < 1)
    return;

  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return;

  for (int nHarm{2}; nHarm <= 3; nHarm++) {
    fW_data->Add(new TH1F(Form("hrmsQ%s%i_FT0C", xy, nHarm), "", CentBin, 0, CentBin));
    TH1F* hRec = reinterpret_cast<TH1F*>(tar->FindObject(Form("hrmsQ%s%i_FT0C", xy, nHarm)));
    if (!hRec)
      return;

    TH2F* hQ = reinterpret_cast<TH2F*>(fW_data->At(0)->Clone(Form("hQ%s%i_FT0C", xy, nHarm)));
    TH1D* h1;

    for (int i{1}; i < hQ->GetXaxis()->GetNbins() + 1; i++) {
      h1 = static_cast<TH1D*>(hQ->ProjectionX(Form("proj%i_Q%s%is", i, xy, nHarm), i, i));

      double mean = h1->GetRMS();
      double meanErr = h1->GetRMSError();
      int binc = hRec->GetXaxis()->GetBinCenter(i);

      hRec->SetBinContent(binc, mean);
      hRec->SetBinError(binc, meanErr);
    }
  }
};

float FFitWeights::GetRecVal(int cent, const char* xy, const int nHarm)
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return -999;

  TH1F* htmp = reinterpret_cast<TH1F*>(tar->FindObject(Form("havgQ%s%i_FT0C", xy, nHarm)));

  return htmp->GetBinContent(cent);
};

float FFitWeights::GetRMSVal(int cent, const char* xy, const int nHarm)
{
  TObjArray* tar{nullptr};

  tar = fW_data;
  if (!tar)
    return -999;

  TH1F* htmp = reinterpret_cast<TH1F*>(tar->FindObject(Form("hrmsQ%s%i_FT0C", xy, nHarm)));

  return htmp->GetBinContent(cent);
};

// float FFitWeights::EventPlane(const float& x, const float& y, const float& nHarm) {
//   return 1/nHarm * TMath::ATan2(y, x);
// };
