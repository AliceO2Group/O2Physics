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
/// \author Junlee Kim (jikim1290@gmail.com)
/// \since November 2021
// code for fitting the luminous region for each step

TFile* fin = new TFile("mergedOutput.root", "read");
TTree* tin;

ULong64_t fTimeStamp;
double fVertexX;
double fVertexY;
double fVertexZ;
double fVertexXY;
double fVertexXX;
double fVertexYY;
double fVertexChi2;
int fNContrib;

double tmin;
double tmax;

int Countig;

void chi2(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t flag)
{
  TMatrixD sigma(2, 2);
  double results = 0.;

  tin = (TTree*)fin->Get("EventInfo_merged");

  tin->SetBranchAddress("fTimeStamp", &fTimeStamp);

  tin->SetBranchAddress("fVertexX", &fVertexX);
  tin->SetBranchAddress("fVertexY", &fVertexY);
  tin->SetBranchAddress("fVertexZ", &fVertexZ);

  tin->SetBranchAddress("fVertexXY", &fVertexXY);
  tin->SetBranchAddress("fVertexXX", &fVertexXX);
  tin->SetBranchAddress("fVertexYY", &fVertexYY);

  tin->SetBranchAddress("fVertexChi2", &fVertexChi2);
  tin->SetBranchAddress("fNContrib", &fNContrib);

  Countig = 0;
  for (int i = 0; i < tin->GetEntries(); i++) {
    tin->GetEntry(i);

    if (fVertexChi2 / fNContrib > 4 || fVertexChi2 < 0)
      continue;
    if (fNContrib < 20 || fNContrib > 2000)
      continue;
    if (fTimeStamp < tmin || fTimeStamp > tmax)
      continue;

    Countig++;
    sigma(0, 0) = par[3] * par[3];
    sigma(1, 1) = par[4] * par[4];
    sigma(0, 1) = par[3] * par[4] * par[6];
    sigma(1, 0) = par[3] * par[4] * par[6];

    sigma(0, 0) += par[9] * par[9] * fVertexXX; // vtx_cxx[i];
    sigma(1, 1) += par[9] * par[9] * fVertexYY; // vtx_cyy[i];
    sigma(0, 1) += par[9] * par[9] * fVertexXY; // vtx_cxy[i];
    sigma(1, 0) += par[9] * par[9] * fVertexXY; // vtx_cxy[i];

    Double_t det(0);
    sigma.InvertFast(&det);
    Double_t sum = TMath::Log(TMath::Power(TMath::TwoPi(), 1.5));
    sum += TMath::Log(det);
    sum += TMath::Log(TMath::Abs(par[5]));
    const Double_t x = fVertexX - par[0] - par[7] * (fVertexZ - par[2]);
    const Double_t y = fVertexY - par[1] - par[8] * (fVertexZ - par[2]);
    const Double_t z = fVertexZ - par[2];

    sum += 0.5 * (x * x * sigma(0, 0) + y * y * sigma(1, 1) +
                  x * y * sigma(0, 1) + x * y * sigma(1, 0));
    sum += 0.5 * z * z / (par[5] * par[5]);

    results += sum;
  }
  f = results;
}

void non_fac()
{
  const int nstep = 6;
  double TMIN[nstep] = {6.15e5, 6.9e5, 7.65e5, 8.65e5, 9.45e5, 10.4e5};
  double TMAX[nstep] = {6.25e5, 7.0e5, 7.75e5, 8.75e5, 9.55e5, 10.5e5};

  double FitRes_cntl[nstep][10];
  double FitRes_stat[nstep][10];
  TMinuit* FunMinuit[nstep];
  int ierflg = 0;
  Double_t arglist[10];

  TH1D* hResults[nstep];

  TH1D* hPreDistPos[nstep][3];
  TH1D* hPreDistSig[nstep][2];
  for (int i = 0; i < nstep; i++) {
    for (int j = 0; j < 2; j++) {
      hPreDistPos[i][j] =
        new TH1D(Form("hPreDistPos_%d_%d", i, j), "", 1000, -0.5, 0.5);
      hPreDistSig[i][j] =
        new TH1D(Form("hPreDistSig_%d_%d", i, j), "", 1000, -1e-4, 9.9e-3);
    }
    hPreDistPos[i][2] =
      new TH1D(Form("hPreDistPos_%d_%d", i, 2), "", 1000, -100, 100);
  }

  tin = (TTree*)fin->Get("EventInfo_merged");

  tin->SetBranchAddress("fTimeStamp", &fTimeStamp);

  tin->SetBranchAddress("fVertexX", &fVertexX);
  tin->SetBranchAddress("fVertexY", &fVertexY);
  tin->SetBranchAddress("fVertexZ", &fVertexZ);

  tin->SetBranchAddress("fVertexXY", &fVertexXY);
  tin->SetBranchAddress("fVertexXX", &fVertexXX);
  tin->SetBranchAddress("fVertexYY", &fVertexYY);

  tin->SetBranchAddress("fVertexChi2", &fVertexChi2);
  tin->SetBranchAddress("fNContrib", &fNContrib);

  for (int i = 0; i < tin->GetEntries(); i++) {
    tin->GetEntry(i);

    if (fVertexChi2 / fNContrib > 4 || fVertexChi2 < 0)
      continue;
    if (fNContrib < 20 || fNContrib > 2000)
      continue;

    for (int s = 0; s < nstep; s++) {
      if (fTimeStamp < TMIN[s] || fTimeStamp > TMAX[s])
        continue;
      hPreDistPos[s][0]->Fill(fVertexX);
      hPreDistPos[s][1]->Fill(fVertexY);
      hPreDistPos[s][2]->Fill(fVertexZ);

      hPreDistSig[s][0]->Fill(sqrt(fVertexXY));
      hPreDistSig[s][1]->Fill(sqrt(fVertexXX));
    }
  }

  for (int i = 0; i < nstep; i++) {
    tmin = TMIN[i];
    tmax = TMAX[i];

    FunMinuit[i] = new TMinuit(10);
    FunMinuit[i]->SetFCN(chi2);
    FunMinuit[i]->SetPrintLevel(1);
    FunMinuit[i]->mnexcm("SET ERR", arglist, 1, ierflg);

    FunMinuit[i]->mnparm(
      0, "", hPreDistPos[i][0]->GetMean(),
      fabs(hPreDistPos[i][0]->GetMean()) * 0.01,
      hPreDistPos[i][0]->GetMean() - hPreDistPos[i][0]->GetRMS() * 10.0,
      hPreDistPos[i][0]->GetMean() + hPreDistPos[i][0]->GetRMS() * 10.0,
      ierflg);
    FunMinuit[i]->mnparm(
      1, "", hPreDistPos[i][1]->GetMean(),
      fabs(hPreDistPos[i][1]->GetMean()) * 0.01,
      hPreDistPos[i][1]->GetMean() - hPreDistPos[i][1]->GetRMS() * 10.0,
      hPreDistPos[i][1]->GetMean() + hPreDistPos[i][1]->GetRMS() * 10.0,
      ierflg);
    FunMinuit[i]->mnparm(2, "", hPreDistPos[i][2]->GetMean(),
                         fabs(hPreDistPos[i][2]->GetMean()) * 0.01, -10, 10,
                         ierflg);

    FunMinuit[i]->mnparm(3, "", hPreDistSig[i][0]->GetMean(),
                         hPreDistSig[i][0]->GetMean() * 0.01, 0,
                         hPreDistSig[i][0]->GetMean() * 10.0, ierflg);
    FunMinuit[i]->mnparm(4, "", hPreDistSig[i][1]->GetMean(),
                         hPreDistSig[i][1]->GetMean() * 0.01, 0,
                         hPreDistSig[i][1]->GetMean() * 10.0, ierflg);
    FunMinuit[i]->mnparm(5, "", 1, 1e-2, 0.0, 200, ierflg);

    FunMinuit[i]->mnparm(6, "", 0.1, 1e-2, -1, 1, ierflg);
    FunMinuit[i]->mnparm(7, "", 0.01, 1e-4, -1, 1, ierflg);
    FunMinuit[i]->mnparm(8, "", 0.01, 1e-4, -1, 1, ierflg);
    FunMinuit[i]->mnparm(9, "", 1.0, 1e-2, 0.5, 5, ierflg);

    FunMinuit[i]->mnexcm("MIGRAD", arglist, 2, ierflg);
    for (int p = 0; p < 10; p++) {
      FunMinuit[i]->GetParameter(p, FitRes_cntl[i][p], FitRes_stat[i][p]);
      cout << FitRes_cntl[i][p] << ", " << FitRes_stat[i][p] << endl;
    }
    hResults[i] = new TH1D(Form("hResults%d", i), "", 12, 0, 12);
    for (int p = 0; p < 10; p++) {
      hResults[i]->SetBinContent(p + 1, FitRes_cntl[i][p]);
      hResults[i]->SetBinError(p + 1, FitRes_stat[i][p]);
    }
    hResults[i]->SetBinContent(11, Countig);
    hResults[i]->SetBinContent(12, FunMinuit[i]->fAmin);

    FunMinuit[i]->DeleteArrays();
  }

  TH1D* hParam[10];
  for (int i = 0; i < 10; i++) {
    hParam[i] = new TH1D(Form("hParam%d", i), "", nstep, 0, nstep);
    for (int j = 0; j < nstep; j++) {
      hParam[i]->SetBinContent(j + 1, FitRes_cntl[j][i]);
      hParam[i]->SetBinError(j + 1, FitRes_stat[j][i]);
    }
  }
  TFile* fout = new TFile("lumiRegionOut.root", "recreate");
  for (int i = 0; i < nstep; i++) {
    hResults[i]->Write();
    for (int j = 0; j < 2; j++) {
      hPreDistPos[i][j]->Write();
      hPreDistSig[i][j]->Write();
    }
    hPreDistPos[i][2]->Write();
  }
  for (int i = 0; i < 10; i++) {
    hParam[i]->Write();
  }
}
