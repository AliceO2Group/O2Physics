// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/// \author Junlee Kim (jikim1290@gmail.com)
/// \since November 2021
// usage: o2-analysis-timestamp -b --aod-file AO2D.root --configuration
// json://./config.json | o2-analysis-trackextension -b |
// o2-analysis-trackselection -b --isRun3 <0, 1> | o2-analysis-mm-lumi -b
// --configuration json://./config.json


void DrawVtx_PbPb(){
  const char fname[1000] = {"upd_full_trk"};
  // org, upd

  TFile *ftot = new TFile(
      Form("../data/PbPb_Run2conv/AnalysisResults_%s.root", fname), "read"); //change me
  TH2F *hTime_VtxX = (TH2F *)ftot->Get("lumi/vertexx_Refitted_timestamp");
  TH2F *hTime_VtxY = (TH2F *)ftot->Get("lumi/vertexy_Refitted_timestamp");

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.03);
  gPad->SetTicks(1);
  gStyle->SetOptStat(0);

  hTime_VtxX->GetXaxis()->SetRangeUser(0, 3200e3);
  hTime_VtxX->GetYaxis()->SetRangeUser(0.08, 0.125);
  hTime_VtxX->GetYaxis()->SetTitle("x (cm)");
  hTime_VtxX->GetXaxis()->SetTitle("t (ms)");
  hTime_VtxX->Draw("colz");
  gPad->Update();

  TPaletteAxis *paletteX =
      (TPaletteAxis *)hTime_VtxX->GetListOfFunctions()->FindObject("palette");
  paletteX->SetY1NDC(0.2);
  c->Update();
  hTime_VtxX->Draw("colz");
  c->SaveAs(Form("../figures_%s/vertexX_297325_full.pdf", fname));

  const int nstage = 3;
  double StepLength[nstage] = {0.0021079, 0.0042159, 0.003243};

  const int nsteps[nstage] = {6, 3, 4};

  double trange_lower_x[3][6] = {{580, 680, 760, 840, 920, 1010},
                                 {1120, 1220, 1300, -1, -1, -1},
                                 {1400, 1500, 1600, 1680, -1, -1}};
  double trange_upper_x[3][6] = {{620, 720, 800, 880, 960, 1050},
                                 {1160, 1260, 1340, -1, -1, -1},
                                 {1440, 1540, 1640, 1720, -1, -1}};

  TLine *line_xl_1st[3][6];
  TLine *line_xu_1st[3][6];

  for (int s = 0; s < nstage; s++) {
    if (s == 0) {
      hTime_VtxX->GetXaxis()->SetRangeUser(500e3, 1200e3);
      hTime_VtxX->Draw("colz");
    } else if (s == 1) {
      hTime_VtxX->GetXaxis()->SetRangeUser(1050e3, 1400e3);
      hTime_VtxX->Draw("colz");
    } else if (s == 2) {
      hTime_VtxX->GetXaxis()->SetRangeUser(1350e3, 1800e3);
      hTime_VtxX->Draw("colz");
    }
    for (int i = 0; i < nsteps[s]; i++) {
      trange_lower_x[s][i] += 20;
      //		trange_upper_x[s][i] -= 5;
      line_xl_1st[s][i] = new TLine(1e3 * trange_lower_x[s][i], 0.08,
                                    1e3 * trange_lower_x[s][i], 0.11);
      line_xu_1st[s][i] = new TLine(1e3 * trange_upper_x[s][i], 0.08,
                                    1e3 * trange_upper_x[s][i], 0.11);

      line_xl_1st[s][i]->SetLineColor(kRed);
      line_xu_1st[s][i]->SetLineColor(kRed);

      line_xl_1st[s][i]->SetLineWidth(3);
      line_xu_1st[s][i]->SetLineWidth(3);

      line_xl_1st[s][i]->Draw("same");
      line_xu_1st[s][i]->Draw("same");
    }
    c->SaveAs(
        Form("../figures_%s/vertexX_297325_full_zoom%d.pdf", fname, s + 1));
  }

  //////////

  hTime_VtxY->GetXaxis()->SetRangeUser(0, 3200e3);
  hTime_VtxY->GetYaxis()->SetRangeUser(0.355, 0.400);
  hTime_VtxY->GetYaxis()->SetTitle("y (cm)");
  hTime_VtxY->GetXaxis()->SetTitle("t (ms)");

  hTime_VtxY->Draw("colz");
  gPad->Update();

  TPaletteAxis *paletteY =
      (TPaletteAxis *)hTime_VtxY->GetListOfFunctions()->FindObject("palette");
  paletteY->SetY1NDC(0.2);
  c->Update();
  hTime_VtxY->Draw("colz");
  c->SaveAs(Form("../figures_%s/vertexY_297325_full.pdf", fname));

  double trange_lower_y[3][6] = {{1820, 1920, 2000, 2100, 2180, 2260},
                                 {2360, 2460, 2540, -1, -1, -1},
                                 {2640, 2740, 2840, 2920, -1, -1}};
  double trange_upper_y[3][6] = {{1860, 1960, 2040, 2140, 2220, 2300},
                                 {2400, 2500, 2580, -1, -1, -1},
                                 {2680, 2780, 2880, 2960, -1, -1}};

  TLine *line_yl_1st[3][6];
  TLine *line_yu_1st[3][6];

  for (int s = 0; s < nstage; s++) {
    if (s == 0) {
      hTime_VtxY->GetXaxis()->SetRangeUser(1700e3, 2350e3);
      hTime_VtxY->Draw("colz");
    } else if (s == 1) {
      hTime_VtxY->GetXaxis()->SetRangeUser(2300e3, 2650e3);
      hTime_VtxY->Draw("colz");
    } else if (s == 2) {
      hTime_VtxY->GetXaxis()->SetRangeUser(2600e3, 3000e3);
      hTime_VtxY->Draw("colz");
    }
    for (int i = 0; i < nsteps[s]; i++) {
      trange_lower_y[s][i] += 20;

      line_yl_1st[s][i] = new TLine(1e3 * trange_lower_y[s][i], 0.355,
                                    1e3 * trange_lower_y[s][i], 0.385);
      line_yu_1st[s][i] = new TLine(1e3 * trange_upper_y[s][i], 0.355,
                                    1e3 * trange_upper_y[s][i], 0.385);

      line_yl_1st[s][i]->SetLineColor(kRed);
      line_yu_1st[s][i]->SetLineColor(kRed);

      line_yl_1st[s][i]->SetLineWidth(3);
      line_yu_1st[s][i]->SetLineWidth(3);

      line_yl_1st[s][i]->Draw("same");
      line_yu_1st[s][i]->Draw("same");
    }
    c->SaveAs(
        Form("../figures_%s/vertexY_297325_full_zoom%d.pdf", fname, s + 1));
  }

  //////////////////

  TH1D *hVtxX_projected[3][6];
  TH1D *hVtxY_projected[3][6];

  double MeanXCntl[3][6];
  double MeanXStat[3][6];
  double MeanYCntl[3][6];
  double MeanYStat[3][6];

  TF1 *f1 = new TF1("f1", "gaus", -10, 10);

  for (int s = 0; s < nstage; s++) {
    for (int i = 0; i < nsteps[s]; i++) {
      hVtxX_projected[s][i] = (TH1D *)hTime_VtxX->ProjectionY(
          Form("hVtxX_projected_%d_%d", s, i),
          hTime_VtxX->GetXaxis()->FindBin(trange_lower_x[s][i] * 1e3),
          hTime_VtxX->GetXaxis()->FindBin(trange_upper_x[s][i] * 1e3), "s");
      hVtxX_projected[s][i]->Rebin(4);

      hVtxX_projected[s][i]->Fit(
          f1, "", "",
          hVtxX_projected[s][i]->GetMean() - hVtxX_projected[s][i]->GetRMS(),
          hVtxX_projected[s][i]->GetMean() + hVtxX_projected[s][i]->GetRMS());

      MeanXCntl[s][i] = f1->GetParameter(1);
      MeanXStat[s][i] = f1->GetParError(1);

      hVtxY_projected[s][i] = (TH1D *)hTime_VtxY->ProjectionY(
          Form("hVtxY_projected_%d_%d", s, i),
          hTime_VtxY->GetXaxis()->FindBin(trange_lower_y[s][i] * 1e3),
          hTime_VtxY->GetXaxis()->FindBin(trange_upper_y[s][i] * 1e3), "s");
      hVtxY_projected[s][i]->Rebin(4);

      hVtxY_projected[s][i]->Fit(
          f1, "", "",
          hVtxY_projected[s][i]->GetMean() - hVtxY_projected[s][i]->GetRMS(),
          hVtxY_projected[s][i]->GetMean() + hVtxY_projected[s][i]->GetRMS());

      MeanYCntl[s][i] = f1->GetParameter(1);
      MeanYStat[s][i] = f1->GetParError(1);
    }
  }

  TLegend *legMoving = new TLegend(0.7, 0.5, 0.93, 0.93);
  legMoving->SetLineWidth(0.0);
  legMoving->SetFillColorAlpha(0, 0);
  gPad->SetLogy();
  for (int s = 0; s < nstage; s++) {
    for (int i = 0; i < nsteps[s]; i++) {
      hVtxX_projected[s][i]->SetMarkerStyle(20);
      hVtxX_projected[s][i]->SetMarkerColor(40 + i);
      hVtxX_projected[s][i]->SetLineColor(40 + i);
      hVtxX_projected[s][i]->GetXaxis()->SetTitle("x (cm)");
      if (i == 0)
        hVtxX_projected[s][i]->Draw("p");
      hVtxX_projected[s][i]->Draw("same,p");
      legMoving->AddEntry(hVtxX_projected[s][i], Form("%d step", i), "lp");
    }
    legMoving->Draw();
    c->SaveAs(Form("../figures_%s/vertexX_1D_297325_%dstg.pdf", fname, s + 1));
    legMoving->Clear();
  }

  for (int s = 0; s < nstage; s++) {
    for (int i = 0; i < nsteps[s]; i++) {
      hVtxY_projected[s][i]->SetMarkerStyle(20);
      hVtxY_projected[s][i]->SetMarkerColor(40 + i);
      hVtxY_projected[s][i]->SetLineColor(40 + i);
      hVtxY_projected[s][i]->GetXaxis()->SetTitle("y (cm)");
      if (i == 0)
        hVtxY_projected[s][i]->Draw("p");
      hVtxY_projected[s][i]->Draw("same,p");
      legMoving->AddEntry(hVtxX_projected[s][i], Form("%d step", i), "lp");
    }
    legMoving->Draw();
    c->SaveAs(Form("../figures_%s/vertexY_1D_297325_%dstg.pdf", fname, s + 1));
    legMoving->Clear();
  }

  TF1 *flinear = new TF1("f1", "pol1", -1, 1);

  TGraphErrors *gX[nstage];
  TGraphErrors *gY[nstage];

  gPad->SetLogy(0);

  TLegend *legx = new TLegend(0.150, 0.549, 0.526, 0.888);
  legx->SetFillColorAlpha(0, 0);
  legx->SetLineWidth(0.0);

  TLegend *legy = new TLegend(0.50, 0.549, 0.876, 0.888);
  legy->SetFillColorAlpha(0, 0);
  legy->SetLineWidth(0.0);

  for (int s = 0; s < nstage; s++) {
    gX[s] = new TGraphErrors();
    gY[s] = new TGraphErrors();

    gX[s]->SetMarkerStyle(20);
    gX[s]->SetMarkerSize(1.5);

    gY[s]->SetMarkerStyle(20);
    gY[s]->SetMarkerSize(1.5);

    gX[s]->GetYaxis()->SetTitle("x (cm)");
    gX[s]->GetXaxis()->SetTitle("LHC Step");

    gY[s]->GetYaxis()->SetTitle("y (cm)");
    gY[s]->GetXaxis()->SetTitle("LHC Step");
    for (int i = 0; i < nsteps[s]; i++) {
      gX[s]->SetPoint(i, StepLength[s] * (double)i, MeanXCntl[s][i]);
      gX[s]->SetPointError(i, 0.0, MeanXStat[s][i]);

      gY[s]->SetPoint(i, StepLength[s] * (double)i, MeanYCntl[s][i]);
      gY[s]->SetPointError(i, 0.0, MeanYStat[s][i]);
    }
    gX[s]->Fit(flinear);
    legx->AddEntry(gX[s], "Data", "pl");
    legx->AddEntry(
        flinear,
        Form("p_{0} + p_{1}*x, #chi^{2} = %.1lf", flinear->GetChisquare()),
        "l");
    legx->AddEntry((TObject *)0,
                   Form("p_{0} = %.3lf #pm %.3lf", flinear->GetParameter(0),
                        flinear->GetParError(0)),
                   "");
    legx->AddEntry((TObject *)0,
                   Form("p_{1} = %.3lf #pm %.3lf", flinear->GetParameter(1),
                        flinear->GetParError(1)),
                   "");

    gX[s]->Draw("AP");
    legx->Draw();
    c->SaveAs(Form("../figures_%s/linearX_%dstg.pdf", fname, s + 1));
    legx->Clear();

    gY[s]->Fit(flinear);
    legy->AddEntry(gY[s], "Data", "pl");
    legy->AddEntry(
        flinear,
        Form("p_{0} + p_{1}*x, #chi^{2} = %.1lf", flinear->GetChisquare()),
        "l");
    legy->AddEntry((TObject *)0,
                   Form("p_{0} = %.3lf #pm %.3lf", flinear->GetParameter(0),
                        flinear->GetParError(0)),
                   "");
    legy->AddEntry((TObject *)0,
                   Form("p_{1} = %.3lf #pm %.3lf", flinear->GetParameter(1),
                        flinear->GetParError(1)),
                   "");

    gY[s]->Draw("AP");
    legy->Draw();
    c->SaveAs(Form("../figures_%s/linearY_%dstg.pdf", fname, s + 1));
    legy->Clear();
  }
  TFile *fout = new TFile(Form("figOut_%s.root", fname), "recreate");
  for (int s = 0; s < nstage; s++) {
    gX[s]->SetName(Form("gX_%s_%d", fname, s));
    gY[s]->SetName(Form("gY_%s_%d", fname, s));

    gX[s]->Write();
    gY[s]->Write();
    for (int i = 0; i < nsteps[s]; i++) {
      hVtxX_projected[s][i]->SetName(
          Form("hVtxX_projected_%s_%d_%d", fname, s, i));
      hVtxY_projected[s][i]->SetName(
          Form("hVtxY_projected_%s_%d_%d", fname, s, i));

      hVtxX_projected[s][i]->Write();
      hVtxY_projected[s][i]->Write();
    }
  }
}
