#include "TH1F.h"
#include "TH2F.h"
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <iostream>

void eap_calculate_efficencies_and_purity(TH1F *pHistIden, TH1F *pHistTPosIden, TH1F *pHistGen,
    TH1F *pHistTracked, TString label, TCanvas *canvas, TLegend *pleg, int style, int color);
void eap_draw_hist(TH1* hist, int markerStyle, int color, const char* title, const char* xLabel,
    const char* yLabel, double minX, double maxX, double minY, double maxY, const char* drawOption);

void efficiency_and_purity(TString filename) {
  const TString analysisFilePath = "AnalysisResults.root";
  const TString mlLabel = "ML";
  const TString nSigmaLabel = "nSigma";

  // Monte Carlo variables
  const char *histGenName = "ml-model-gen-hists/hPtMCPositive";
  const char *histGenTrackedName = "ml-model-gen-hists/hPtMCTracked";

  // ML variables
  const char *histMlIdenName = "ml-model-gen-hists/hPtMLPositive";
  const char *histMlTPosIdenName = "ml-model-gen-hists/hPtMLTruePositive";

  // nSigma variables
  const char *histNSigmaIdenName = "ml-model-gen-hists/hPtNSigmaPositive";
  const char *histNSigmaTPosIdenName = "ml-model-gen-hists/hPtNSigmaTruePositive";

  // Context detectors's data
  const char *histTOFNSigmaName = "ml-model-gen-hists/hPtTOFNSigma";
  const char *histTOFBetaName = "ml-model-gen-hists/hPtTOFBeta";
  const char *histTPCNSigmaName = "ml-model-gen-hists/hPtTPCNSigma";
  const char *histTPCSignalName = "ml-model-gen-hists/hPtTPCSignal";

  TCanvas *canvas = new TCanvas("canvas", "canvas", 0, 0, 1920, 2160);
  canvas->Divide(2, 4);
  TLegend *pleg = new TLegend(0.15, 0.1); 

  gStyle->SetPalette(kRainBow);

  TFile *f = new TFile(analysisFilePath.Data());
  if (f->IsZombie()) {
    printf("Failed to open file %s\n", analysisFilePath.Data());
    return;
  }

  // Monte Carlo variables
  TH1F *pHistGen = (TH1F *)f->Get(histGenName);
  TH1F *pHistTracked = (TH1F *)f->Get(histGenTrackedName);

  // ML variables
  TH1F *pHistMlIden = (TH1F *)f->Get(histMlIdenName);
  TH1F *pHistMlTPosIden = (TH1F *)f->Get(histMlTPosIdenName);

  // nSigma variables
  TH1F *pHistNSigmaIden = (TH1F *)f->Get(histNSigmaIdenName);
  TH1F *pHistNSigmaTPosIden = (TH1F *)f->Get(histNSigmaTPosIdenName);

  eap_calculate_efficencies_and_purity(pHistMlIden, pHistMlTPosIden, pHistGen, pHistTracked, mlLabel, canvas, pleg, 20, 38);
  eap_calculate_efficencies_and_purity(pHistNSigmaIden, pHistNSigmaTPosIden, pHistGen, pHistTracked, nSigmaLabel, canvas, pleg, 23, 46);

  int contextStyle = 24, contextColor = 2;

  TH2F *pHistTOFNSigma = (TH2F *)f->Get(histTOFNSigmaName);
  canvas->cd(5);
  eap_draw_hist(pHistTOFNSigma, contextStyle, contextColor, "TOFNSigma", "#it{p}_{T} (GeV/#it{c})", "TOF NSigma", 0.0, 3.1, -5.0, 5.0, "colz");

  TH2F *pHistTPCNSigma = (TH2F *)f->Get(histTPCNSigmaName);
  canvas->cd(6);
  eap_draw_hist(pHistTPCNSigma, contextStyle, contextColor, "TPCNSigma", "#it{p}_{T} (GeV/#it{c})", "TPC NSigma", 0.0, 3.1, -5.0, 5.0, "colz");
  gPad->SetLogz();

  TH2F *pHistTOFBeta = (TH2F *)f->Get(histTOFBetaName);
  canvas->cd(7);
  eap_draw_hist(pHistTOFBeta, contextStyle, contextColor, "TOFBeta", "#it{p}_{T} (GeV/#it{c})", "TOF Beta", 0.0, 3.1, 0.0, 1.2, "colz");
  gPad->SetLogz();

  TH2F *pHistTPCSignal = (TH2F *)f->Get(histTPCSignalName);
  canvas->cd(8);
  eap_draw_hist(pHistTPCSignal, contextStyle, contextColor, "TPCSignal", "#it{p}_{T} (GeV/#it{c})", "TPC Signal", 0.0, 3.1, 20.0, 120.0, "colz");
  gPad->SetLogz();

  canvas->cd(2);
  pleg->Draw();
  canvas->Print(Form("./graphs/%s.png", filename.Data()));
}

void eap_calculate_efficencies_and_purity(TH1F *pHistIden, TH1F *pHistTPosIden,
    TH1F *pHistGen, TH1F *pHistTracked, TString label, TCanvas *canvas, TLegend *pleg, int style, int color) {
  const char* xTitle = "#it{p}_{T} (GeV/#it{c})";

  // Full Efficiency
  canvas->cd(1);

  TH1F *phFullTrackEff = (TH1F *)pHistTPosIden->Clone();
  phFullTrackEff->Sumw2();
  phFullTrackEff->Divide(pHistGen);

  eap_draw_hist(phFullTrackEff, style, color, "Full efficiency", xTitle, "Full efficiency", 0.0, 3.1, 0.0, 1.0, "he,same");

  // PID Efficiency
  canvas->cd(2);

  TH1F *phPidTrackEff = (TH1F *)pHistTPosIden->Clone();
  phPidTrackEff->Sumw2();
  phPidTrackEff->Divide(pHistTracked);

  eap_draw_hist(phPidTrackEff, style, color, "PID efficiency", xTitle, "PID efficiency", 0.0, 3.1, 0.0, 1.0, "he,same");

  // Reconstruction Efficiency
  canvas->cd(3);

  TH1F *phReconstructionTrackEff = (TH1F *)pHistTracked->Clone();
  phReconstructionTrackEff->Sumw2();
  phReconstructionTrackEff->Divide(pHistGen);

  eap_draw_hist(phReconstructionTrackEff, style, color, "Reconstruction Efficiency", xTitle, "Reconstruction Efficiency", 0.0, 3.1, 0.0, 1.0, "he,same");

  // Purity
  canvas->cd(4);

  TH1F *hPurity = (TH1F *)pHistTPosIden->Clone();
  pHistIden->Sumw2();
  hPurity->Divide(pHistIden);
  hPurity->Sumw2();

  eap_draw_hist(hPurity, style, color, "Purity", xTitle, "Purity", 0.0, 3.1, 0.0, 1.0, "he,same");

  // Add to legend
  pleg->AddEntry(hPurity, label);
}

void eap_draw_hist(TH1* hist, int markerStyle, int color, const char* title, const char* xLabel,
    const char* yLabel, double minX, double maxX, double minY, double maxY, const char* drawOption) {
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(1.4);
    hist->SetTitle(title);
    hist->GetXaxis()->SetTitle(xLabel);
    hist->GetYaxis()->SetTitle(yLabel);
    hist->GetXaxis()->SetRangeUser(minX, maxX);
    hist->GetYaxis()->SetRangeUser(minY, maxY);
    hist->GetXaxis()->SetTitleSize(0.03);
    hist->GetYaxis()->SetTitleSize(0.03);
    hist->GetXaxis()->SetLabelSize(0.03);
    hist->GetYaxis()->SetLabelSize(0.03);
    hist->SetStats(0);
    hist->SetContour(1000);
    hist->Draw(drawOption);
}

