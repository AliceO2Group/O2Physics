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

const int nHists = 7;
TH1F* mee[nHists];
TH1F* mee_orig[nHists];
TString histNames[nHists] = {"", "_Pi0", "_Eta", "_EtaP",
                             "_Rho", "_Omega", "_Phi"};
TString histLegends[nHists] = {
  "Cocktail sum",
  "#pi^{0}#rightarrow#gammae^{+}e^{-}",
  "#eta#rightarrow#gammae^{+}e^{-}",
  "#eta'#rightarrow#gammae^{+}e^{-}, #eta'#rightarrow#omegae^{+}e^{-}",
  "#rho#rightarrowe^{+}e^{-}",
  "#omega#rightarrow#pi^{0}e^{+}e^{-}, #omega#rightarrowe^{+}e^{-}",
  "#phi#rightarrow#etae^{+}e^{-}, #phi#rightarrow#pi^{0}e^{+}e^{-}, "
  "#phi#rightarrowe^{+}e^{-}"};

void loadHistos(TFile* file, TH1F* hists[], TString name_extra, int rebin,
                int nEvents)
{
  for (int i = 0; i < nHists; i++) {
    hists[i] = (TH1F*)file->GetDirectory("em-lmee-lf-cocktail")
                 ->Get(TString("mee") + name_extra + histNames[i]);
    hists[i]->Rebin(rebin);
    hists[i]->Scale(1. / nEvents);
    hists[i]->SetTitle(histLegends[i]);
  }
}

void plotLFCocktail(TString filename = "AnalysisResults.root", int rebin = 1)
{

  TFile* file = TFile::Open(filename.Data());

  TH1I* hNEvents =
    (TH1I*)file->GetDirectory("em-lmee-lf-cocktail")->Get("NEvents");
  int nEvents = hNEvents->GetBinContent(1);

  loadHistos(file, mee, "", rebin, nEvents);
  loadHistos(file, mee_orig, "_orig", rebin, nEvents);

  auto canvas = new TCanvas("LF Cocktail", "LF Cocktail", 1500, 500);
  canvas->Divide(2, 1);

  mee[0]->SetStats(0);
  mee[0]->GetXaxis()->SetRangeUser(0.0, 1.1);
  mee_orig[0]->SetStats(0);
  mee_orig[0]->GetXaxis()->SetRangeUser(0.0, 1.1);

  canvas->cd(1);
  gPad->SetLogy();
  for (int i = 0; i < nHists; i++) {
    mee_orig[i]->SetLineColor(i + 1);
    mee_orig[i]->Draw("hist same");
  }
  gPad->BuildLegend(0.62, 0.9, 0.9, 0.6);
  mee_orig[0]->SetTitle("before resolution and acceptance");
  canvas->cd(2);
  gPad->SetLogy();
  for (int i = 0; i < nHists; i++) {
    mee[i]->SetLineColor(i + 1);
    mee[i]->Draw("hist same");
  }
  gPad->BuildLegend(0.62, 0.9, 0.9, 0.6);
  mee[0]->SetTitle("after resolution and acceptance");
}