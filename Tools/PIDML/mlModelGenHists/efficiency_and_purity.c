#include "TH1F.h"
#include <TCanvas.h>
#include <TFile.h>
#include <iostream>

void eap_draw_hist(TH1F* hist, int markerStyle, int color, const char* title, double maxX, double maxY);

void efficiency_and_purity() {
  const TString myfile[] = {"AnalysisResults.root"};
  const TString label[] = {"ML"};
  const char *histGenName = "ml-model-gen-hists/hPtMCPositive";
  const char *histGenTrackedName = "ml-model-gen-hists/hPtMCTracked";
  const char *histIdenName = "ml-model-gen-hists/hPtMLPositive";
  const char *histTPosIdenName = "ml-model-gen-hists/hPtMLTruePositive";
  const int filesCount = sizeof(myfile)/sizeof(myfile[0]);

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->Divide(2, filesCount);
  TLegend *pleg = new TLegend(0.4, 0.1); 

  for (int i = 0; i < filesCount; i++) {
    TFile *f = new TFile(myfile[i].Data());
    if (f->IsZombie()) {
      printf("Failed to open file %s\n", myfile[i].Data());
      return;
    }
    
    TH1F *pHistGen = (TH1F *)f->Get(histGenName);
    TH1F *pHistTracked = (TH1F *)f->Get(histGenTrackedName);
    TH1F *pHistIden = (TH1F *)f->Get(histIdenName);
    TH1F *pHistTPosIden = (TH1F *)f->Get(histTPosIdenName);

    TH1F *phFullTrackEff = (TH1F *)pHistTPosIden->Clone();
    phFullTrackEff->Sumw2();
    phFullTrackEff->Divide(pHistGen);

    TH1F *phPidTrackEff = (TH1F *)pHistTPosIden->Clone();
    phPidTrackEff->Sumw2();
    phPidTrackEff->Divide(pHistTracked);

    TH1F *phReconstructionTrackEff = (TH1F *)pHistTracked->Clone();
    phReconstructionTrackEff->Sumw2();
    phReconstructionTrackEff->Divide(pHistGen);

    TH1F *hPurity = (TH1F *)pHistTPosIden->Clone();
    pHistIden->Sumw2();
    hPurity->Divide(pHistIden);
    hPurity->Sumw2();
    pleg->AddEntry(hPurity, "Purity");

    if (i == 0) {
      pleg->AddEntry(phReconstructionTrackEff, "Reco eff.");
      pleg->AddEntry(phFullTrackEff, "Full eff.");
      pleg->AddEntry(phPidTrackEff, "PID eff.");
    }

    canvas->cd(i*2 + 1);

    TString effTitle = Form("Efficiency %s", label[i].Data());

    eap_draw_hist(phReconstructionTrackEff, 20, 2, effTitle, 3.1, 2.5);
    eap_draw_hist(phFullTrackEff, 21, 3, effTitle, 3.1, 2.5);
    eap_draw_hist(phPidTrackEff, 22, 4, effTitle, 3.1, 2.5);
   
    canvas->cd(i*2 + 2);

    eap_draw_hist(hPurity, 23, 3, Form("Purity %s", label[i].Data()), 3.0, 1.0);
}

  canvas->cd(1);
  pleg->Draw();
}

void eap_draw_hist(TH1F* hist, int markerStyle, int color, const char* title,
    double maxX, double maxY) {
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(1.8);
    hist->SetTitle(title);
    hist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hist->GetYaxis()->SetTitle(title);
    hist->GetXaxis()->SetRangeUser(0.0, maxX);
    hist->GetYaxis()->SetRangeUser(0.0, maxY);
    hist->GetXaxis()->SetTitleSize(0.03);
    hist->GetYaxis()->SetTitleSize(0.03);
    hist->GetXaxis()->SetLabelSize(0.03);
    hist->GetYaxis()->SetLabelSize(0.03);
    hist->SetStats(0);
    hist->Draw("he,same");
}

