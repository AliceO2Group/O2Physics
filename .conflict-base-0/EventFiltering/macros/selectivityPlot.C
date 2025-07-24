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

void selectivityPlot(int runNumber = 550781, TString inputfile = "AnalysisResults_550781.root", TString outputfolder = "")
{
  TCanvas* Canvas_1 = new TCanvas("Canvas_1", "Canvas_1", 928, 592);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  Canvas_1->Range(-6.079755, -8.734506, 67.09509, -2.808365);
  Canvas_1->SetFillColor(0);
  Canvas_1->SetBorderMode(0);
  Canvas_1->SetBorderSize(2);
  Canvas_1->SetLogy();
  Canvas_1->SetGridx();
  Canvas_1->SetGridy();
  Canvas_1->SetTickx(1);
  Canvas_1->SetTicky(1);
  Canvas_1->SetLeftMargin(0.08991826);
  Canvas_1->SetRightMargin(0.02179836);
  Canvas_1->SetTopMargin(0.05);
  Canvas_1->SetBottomMargin(0.2711864);
  Canvas_1->SetFrameBorderMode(0);
  Canvas_1->SetFrameBorderMode(0);

  TFile* file = new TFile(inputfile.Data());
  TH1D* mFiltered = (TH1D*)file->Get("central-event-filter-task/scalers/mFiltered");
  double nEvents = mFiltered->GetBinContent(1);
  mFiltered->Scale(1. / nEvents);
  mFiltered->SetLineColor(2);
  mFiltered->SetLineWidth(2);
  mFiltered->SetMarkerColor(2);
  mFiltered->SetMarkerStyle(4);
  mFiltered->SetMarkerSize(0.6);
  mFiltered->GetXaxis()->SetRange(2, mFiltered->GetXaxis()->GetNbins());
  for (int i = 1; i <= mFiltered->GetNbinsX(); i++) {
    std::string label = mFiltered->GetXaxis()->GetBinLabel(i);
    if (label[0] == 'f') {
      mFiltered->GetXaxis()->SetBinLabel(i, label.substr(1).c_str());
    }
    mFiltered->SetBinContent(i, mFiltered->GetBinContent(i) / mFiltered->GetBinContent(1));
    mFiltered->SetBinError(i, mFiltered->GetBinError(i) / mFiltered->GetBinContent(1));
  }
  mFiltered->GetXaxis()->SetBinLabel(mFiltered->GetNbinsX(), "Total");
  mFiltered->GetXaxis()->SetLabelFont(42);
  mFiltered->GetXaxis()->SetTitleOffset(1);
  mFiltered->GetXaxis()->SetTitleFont(42);
  mFiltered->GetYaxis()->SetTitle("Async. trigger selectivity");
  mFiltered->GetYaxis()->SetLabelFont(42);
  mFiltered->GetYaxis()->SetTitleFont(42);
  mFiltered->GetZaxis()->SetLabelFont(42);
  mFiltered->GetZaxis()->SetTitleOffset(1);
  mFiltered->GetZaxis()->SetTitleFont(42);
  mFiltered->DrawClone("");
  TLine* line = new TLine(mFiltered->GetXaxis()->GetBinLowEdge(2), 5.e-05, mFiltered->GetXaxis()->GetBinLowEdge(mFiltered->GetNbinsX()), 5.e-05);
  line->SetLineStyle(2);
  line->Draw();
  TArrow* arrow = new TArrow(mFiltered->GetXaxis()->GetXmax(), 0.0005, mFiltered->GetXaxis()->GetXmax() + 1, 0.0005, 0.01, "<|");

  Int_t ci;      // for color index setting
  TColor* color; // for color definition with alpha
  ci = TColor::GetColor("#0000ff");
  arrow->SetFillColor(ci);
  arrow->SetFillStyle(1001);
  TLatex* tex = new TLatex();
  tex->SetTextSize(0.035);
  tex->SetTextFont(42);
  tex->DrawLatexNDC(0.09, 0.96, Form("ALICE Internal, pp #sqrt{s} = 13.6 TeV, Run %d, #it{L}_{int} #approx %.1f nb^{-1}", runNumber, nEvents / 78.e6));
  ci = TColor::GetColor("#0000ff");
  arrow->SetLineColor(ci);
  arrow->Draw();
  Canvas_1->Modified();
  Canvas_1->SetSelected(Canvas_1);

  if (outputfolder.IsNull()) {
    return;
  }
  gROOT->SetBatch(true);
  TCanvas Canvas_2("Canvas_2", "Canvas_2", Canvas_1->GetWw() * 2, Canvas_1->GetWh() * 2);
  Canvas_1->DrawClonePad();
  Canvas_2.SaveAs(Form("%s/asyncTriggerSelectivity_%d.png", outputfolder.Data(), runNumber));
  gROOT->SetBatch(false);
}
