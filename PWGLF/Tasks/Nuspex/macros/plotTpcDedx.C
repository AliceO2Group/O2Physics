// Publication-style drawing of the dE/dx histograms written by tpcDedxAntiLi6Search.
//
// Usage:
//   root -l -b -q 'plotTpcDedx.C("AnalysisResults.root", "tpc-dedx-antili6-search", "dEdxVsPLog", "dedx_log.pdf")'
//
// Optional overlay: pass the 5 Bethe-Bloch parameters (+ MIP, charge factor, mass, Z) of a species to draw its curve.

#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>

#include <cmath>
#include <string>

namespace
{
double betheBloch(double bg, const double* k)
{
  const double beta = bg / std::sqrt(1. + bg * bg);
  const double aa = std::pow(beta, k[3]);
  const double bb = std::log(k[2] + std::pow(1. / bg, k[4]));
  return (k[1] - aa - bb) * k[0] / aa;
}
} // namespace

void plotTpcDedx(const char* fileName = "AnalysisResults.root",
                 const char* taskDir = "tpc-dedx-antili6-search",
                 const char* histName = "dEdxVsPLinear",
                 const char* outName = "dedx.pdf",
                 const char* label = "ALICE Performance",
                 const char* system = "Run 3 pp, #sqrt{#it{s}} = 13.6 TeV",
                 bool drawLi6Curve = false,
                 double bb0 = 0.0320981, double bb1 = 19.9768, double bb2 = 2.52666e-16, double bb3 = 2.72123, double bb4 = 6.08092,
                 double mip = 50., double chargeFactor = 2.3, double mass = 5.60152, double charge = 3.)
{
  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    printf("Cannot open %s\n", fileName);
    return;
  }
  auto* h = dynamic_cast<TH2*>(f->Get(Form("%s/%s", taskDir, histName)));
  if (!h) {
    printf("Histogram %s/%s not found\n", taskDir, histName);
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(255);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetTitleFont(42, "xyz");

  auto* c = new TCanvas("cDedx", "cDedx", 950, 650);
  c->SetLogz();
  c->SetTicks(1, 1);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.13);
  c->SetTopMargin(0.03);

  const bool isLog = std::string(histName).find("Log") != std::string::npos;
  if (isLog) {
    // signed logarithmic bins cannot be shown on a log axis: draw as is, bin edges are logarithmic per sign
    h->GetXaxis()->SetTitle("#it{p} / #it{Z} (GeV/#it{c})");
  }
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->SetMinimum(1);
  h->Draw("COL");

  if (drawLi6Curve) {
    const double k[5] = {bb0, bb1, bb2, bb3, bb4};
    const double xmax = h->GetXaxis()->GetXmax();
    for (int sign : {-1, 1}) {
      auto* fn = new TF1(Form("fLi6_%d", sign), [=](double* x, double*) {
        const double p = std::abs(x[0]);
        return mip * betheBloch(charge * p / mass, k) * std::pow(charge, chargeFactor); }, sign < 0 ? -xmax : 0.05, sign < 0 ? -0.05 : xmax, 0);
      fn->SetNpx(2000);
      fn->SetLineColor(kBlack);
      fn->SetLineWidth(2);
      fn->Draw("same");
    }
  }

  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(0.04);
  t.DrawLatex(0.62, 0.90, label);
  t.SetTextSize(0.034);
  t.DrawLatex(0.62, 0.85, system);

  c->SaveAs(outName);
}
