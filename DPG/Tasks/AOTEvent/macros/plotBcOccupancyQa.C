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

/// \file plotBcOccupancyQa.C
/// \brief Plots and summary of the BC occupancy QA (output of o2-analysis-bc-occupancy-qa)
///
/// Usage:
///   root -l -b -q 'plotBcOccupancyQa.C("AnalysisResults.root")'
///   root -l -b -q 'plotBcOccupancyQa.C("AnalysisResults.root", "bcOccupancyQa.pdf", "bc-occupancy-qa-task", "bcOccupancyQa.root")'
/// Output: multi-page PDF, and a ROOT file with the canvases and the histograms (by default: name of the PDF with .root)
///
/// \author Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace bcoccqa
{
TDirectory* gDir = nullptr;

// histogram from the task directory (or from one of its subdirectories), detached from the file
template <typename T>
T* getHist(const char* name)
{
  TObject* obj = gDir->Get(name);
  if (!obj) {
    for (TObject* key : *gDir->GetListOfKeys()) {
      if (auto* sub = dynamic_cast<TDirectory*>(static_cast<TKey*>(key)->ReadObj())) {
        obj = sub->Get(name);
        if (obj) {
          break;
        }
      }
    }
  }
  auto* h = dynamic_cast<T*>(obj);
  if (!h) {
    std::printf("histogram %s not found\n", name);
    return nullptr;
  }
  h->SetDirectory(nullptr);
  h->SetStats(false); // also when the saved canvases are drawn later
  return h;
}

// fraction of BCs per bin: numerator / denominator (binomial errors)
TH1D* makeRatio(const TH1* num, const TH1* den, const char* name, const char* title)
{
  auto* r = static_cast<TH1D*>(num->Clone(name));
  r->SetTitle(title);
  r->Divide(num, den, 1., 1., "B");
  r->SetStats(false);
  return r;
}

// upper edge of the last filled bin (with 10% margin), to zoom the occupancy axes (e.g. for pp)
double upperEdgeOfFilledRange(const TH1* h)
{
  int last = 1;
  for (int i = h->GetNbinsX(); i >= 1; i--) {
    if (h->GetBinContent(i) > 0) {
      last = i;
      break;
    }
  }
  return std::min(h->GetXaxis()->GetXmax(), 1.1 * h->GetXaxis()->GetBinUpEdge(last));
}
} // namespace bcoccqa

void plotBcOccupancyQa(const char* fileName = "AnalysisResults.root", const char* outPdf = "bcOccupancyQa.pdf", const char* dirName = "bc-occupancy-qa-task", const char* outRoot = "")
{
  using namespace bcoccqa;
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  std::unique_ptr<TFile> f(TFile::Open(fileName, "READ"));
  if (!f || f->IsZombie()) {
    std::printf("cannot open %s\n", fileName);
    return;
  }
  gDir = f->GetDirectory(dirName);
  if (!gDir) {
    std::printf("directory %s not found in %s\n", dirName, fileName);
    return;
  }

  auto* hBcCounter = getHist<TH1>("hBcCounter");
  auto* hNumCollsInBC = getHist<TH1>("hNumCollsInBC");
  auto* hBcInTF = getHist<TH1>("hBcInTF");
  auto* hBcInTFDefined = getHist<TH1>("hBcInTFOccupancyDefined");
  auto* hOccTracks = getHist<TH1>("hOccupancyByTracks");
  auto* hOccFT0C = getHist<TH1>("hOccupancyByFT0C");
  auto* hOccTracksVsFT0C = getHist<TH2>("hOccupancyByTracksVsFT0C");
  auto* hMedianTime = getHist<TH1>("hMedianTime");
  auto* pOccVsBcInTF = getHist<TProfile>("pOccupancyByTracksVsBcInTF");
  auto* pOccVsBcInOrbit = getHist<TProfile>("pOccupancyByTracksVsBcInOrbit");
  auto* hFlags = getHist<TH1>("hOccupancyFlags");
  auto* hFlagsVsOcc = getHist<TH2>("hOccupancyFlagsVsOccupancyByTracks");
  auto* hBcIndexCheck = getHist<TH1>("hBcIndexCheck");

  // ### summary
  std::printf("\n=== BC occupancy QA: %s ===\n", fileName);
  if (hBcCounter) {
    double nAll = hBcCounter->GetBinContent(1), nDef = hBcCounter->GetBinContent(2);
    std::printf("BCs: %.0f, with calculated occupancy: %.0f (%.2f%%)\n", nAll, nDef, nAll > 0 ? 100. * nDef / nAll : 0.);
  }
  if (hOccTracks && hOccFT0C) {
    std::printf("mean occupancy: by ITS tracks %.1f, by FT0C %.1f\n", hOccTracks->GetMean(), hOccFT0C->GetMean());
  }
  if (hNumCollsInBC) {
    double nAll = hNumCollsInBC->GetEntries();
    std::printf("BCs with >= 1 collision: %.4f%%, with >= 2 collisions: %.4f%%\n",
                nAll > 0 ? 100. * hNumCollsInBC->Integral(hNumCollsInBC->FindBin(1), hNumCollsInBC->GetNbinsX() + 1) / nAll : 0.,
                nAll > 0 ? 100. * hNumCollsInBC->Integral(hNumCollsInBC->FindBin(2), hNumCollsInBC->GetNbinsX() + 1) / nAll : 0.);
  }
  if (hFlags && hFlags->GetBinContent(1) > 0) {
    std::printf("fraction of BCs (with calculated occupancy) passing each flag:\n");
    for (int i = 2; i <= hFlags->GetNbinsX(); i++) {
      std::printf("  %-28s %.4f\n", hFlags->GetXaxis()->GetBinLabel(i), hFlags->GetBinContent(i) / hFlags->GetBinContent(1));
    }
  }

  if (hBcIndexCheck) {
    // all rows and all found BCs should be found at the right place
    const bool ok = hBcIndexCheck->GetBinContent(2) == hBcIndexCheck->GetBinContent(1) && hBcIndexCheck->GetBinContent(5) == hBcIndexCheck->GetBinContent(4);
    std::printf("lookup by global BC (BcOccAllIndex): rows found %.0f / %.0f, BCs of aod::BCs found %.0f / %.0f (not found %.0f), with the same global BC %.0f -> %s\n",
                hBcIndexCheck->GetBinContent(2), hBcIndexCheck->GetBinContent(1), hBcIndexCheck->GetBinContent(4), hBcIndexCheck->GetBinContent(3),
                hBcIndexCheck->GetBinContent(6), hBcIndexCheck->GetBinContent(5), ok ? "OK" : "MISMATCH");
  }

  // ### plots
  const double maxOccTracks = hOccTracks ? upperEdgeOfFilledRange(hOccTracks) : 0.;
  const double maxOccFT0C = hOccFT0C ? upperEdgeOfFilledRange(hOccFT0C) : 0.;
  // one canvas per page, with the size adapted to the number of pads
  std::vector<TCanvas*> canvases;
  std::vector<TH1*> derivedHists; // histograms made here (fractions), also saved to the ROOT file
  auto newPage = [&canvases](const char* name, const char* title, int width, int height, int nx, int ny) {
    auto* canvas = new TCanvas(name, title, width, height);
    canvas->Divide(nx, ny);
    canvases.push_back(canvas);
    return canvas;
  };
  TCanvas* c = nullptr;

  // page 1: BC counters and the cross-check of the lookup by global BC
  c = newPage("cBcCounters", "BC counters", 1500, 550, 3, 1);
  if (hBcCounter) {
    c->cd(1)->SetLeftMargin(0.15);
    hBcCounter->SetTitle("number of BCs");
    hBcCounter->Draw("hist text");
  }
  if (hNumCollsInBC) {
    c->cd(2)->SetLogy();
    hNumCollsInBC->SetTitle("all BCs;number of collisions found in the BC;N BCs");
    hNumCollsInBC->Draw("hist");
  }
  if (hBcIndexCheck) {
    c->cd(3)->SetBottomMargin(0.3);
    c->cd(3)->SetLogy(); // rows of BcOccsAll vs (much fewer) BCs of aod::BCs
    hBcIndexCheck->SetTitle("lookup in BcOccsAll by global BC;;counts");
    hBcIndexCheck->GetXaxis()->LabelsOption("v");
    hBcIndexCheck->SetMinimum(0.5);
    hBcIndexCheck->Draw("hist text");
  }

  // page 2: occupancy distributions
  c = newPage("cOccupancy", "occupancy distributions", 1200, 800, 2, 2);
  if (hOccTracks) {
    c->cd(1)->SetLogy();
    hOccTracks->SetTitle("occupancy by ITS tracks;occupancy by ITS tracks;N BCs");
    hOccTracks->GetXaxis()->SetRangeUser(0., maxOccTracks);
    hOccTracks->Draw("hist");
  }
  if (hOccFT0C) {
    c->cd(2)->SetLogy();
    hOccFT0C->SetTitle("occupancy by FT0C;occupancy by FT0C amplitudes;N BCs");
    hOccFT0C->GetXaxis()->SetRangeUser(0., maxOccFT0C);
    hOccFT0C->Draw("hist");
  }
  if (hOccTracksVsFT0C) {
    c->cd(3)->SetLogz();
    hOccTracksVsFT0C->SetTitle("occupancy: ITS tracks vs FT0C;occupancy by FT0C amplitudes;occupancy by ITS tracks");
    hOccTracksVsFT0C->GetXaxis()->SetRangeUser(0., maxOccFT0C);
    hOccTracksVsFT0C->GetYaxis()->SetRangeUser(0., maxOccTracks);
    hOccTracksVsFT0C->Draw("colz");
  }
  if (hMedianTime) {
    c->cd(4)->SetLogy();
    hMedianTime->SetTitle("median time of the occupancy;#Deltat (#mus);N BCs");
    hMedianTime->Draw("hist");
  }

  // page 3: occupancy flags
  c = newPage("cOccupancyFlags", "occupancy flags", 1200, 550, 2, 1);
  if (hFlags && hFlags->GetBinContent(1) > 0) {
    c->cd(1)->SetBottomMargin(0.3);
    auto* hFlagsFrac = static_cast<TH1*>(hFlags->Clone("hFlagsFrac"));
    hFlagsFrac->Scale(1. / hFlags->GetBinContent(1));
    hFlagsFrac->SetTitle("fraction of BCs passing each flag;;fraction of BCs");
    hFlagsFrac->GetXaxis()->LabelsOption("v");
    hFlagsFrac->SetMinimum(0.);
    hFlagsFrac->SetMaximum(1.1);
    hFlagsFrac->Draw("hist text");
    derivedHists.push_back(hFlagsFrac);
  }
  if (hFlagsVsOcc) {
    // fraction of BCs passing each flag vs occupancy
    c->cd(2);
    TLegend* leg = new TLegend(0.45, 0.6, 0.89, 0.89);
    TH1D* hAll = hFlagsVsOcc->ProjectionY("hFlagsAll", 1, 1);
    const int colors[] = {kRed + 1, kOrange + 7, kGreen + 2, kAzure + 2, kBlue + 2, kMagenta + 1, kBlack};
    bool first = true;
    for (int i = 2; i <= hFlagsVsOcc->GetNbinsX(); i++) {
      TH1D* hPass = hFlagsVsOcc->ProjectionY(Form("hFlagPass%d", i), i, i);
      TH1D* hFrac = makeRatio(hPass, hAll, Form("hFracVsOccupancy_%s", hFlagsVsOcc->GetXaxis()->GetBinLabel(i)), "fraction of BCs passing the flags vs occupancy;occupancy by ITS tracks;fraction of BCs");
      hFrac->SetLineColor(colors[(i - 2) % 7]);
      hFrac->SetMarkerColor(colors[(i - 2) % 7]);
      hFrac->SetMarkerStyle(20);
      hFrac->SetMarkerSize(0.4);
      hFrac->SetMinimum(0.);
      hFrac->SetMaximum(1.1);
      hFrac->GetXaxis()->SetRangeUser(0., maxOccTracks);
      hFrac->Draw(first ? "p" : "p same");
      derivedHists.push_back(hFrac);
      leg->AddEntry(hFrac, hFlagsVsOcc->GetXaxis()->GetBinLabel(i), "p");
      first = false;
    }
    leg->Draw();
  }

  // page 4: BC position in the TF and in the orbit
  c = newPage("cBcPosition", "BC position in the TF and in the orbit", 1200, 800, 1, 3);
  if (hBcInTF && hBcInTFDefined) {
    c->cd(1);
    TH1D* hFracDefined = makeRatio(hBcInTFDefined, hBcInTF, "hFracDefined", "fraction of BCs with calculated occupancy;BC in TF;fraction");
    hFracDefined->SetMinimum(0.);
    hFracDefined->SetMaximum(1.1);
    hFracDefined->Draw("hist");
    derivedHists.push_back(hFracDefined);
  }
  if (pOccVsBcInTF) {
    c->cd(2);
    pOccVsBcInTF->SetTitle("mean occupancy vs BC in TF;BC in TF;mean occupancy by ITS tracks");
    pOccVsBcInTF->Draw();
  }
  if (pOccVsBcInOrbit) {
    c->cd(3);
    pOccVsBcInOrbit->SetTitle("mean occupancy vs BC in orbit;BC in orbit;mean occupancy by ITS tracks");
    pOccVsBcInOrbit->Draw("hist"); // 3564 bins: error bars would hide the pattern
  }

  // PDF: each page with the aspect ratio of its canvas
  for (size_t i = 0; i < canvases.size(); i++) {
    gStyle->SetPaperSize(30., 30. * canvases[i]->GetWh() / canvases[i]->GetWw());
    if (i == 0) {
      canvases[i]->Print(Form("%s[", outPdf));
    }
    canvases[i]->Print(outPdf);
    if (i == canvases.size() - 1) {
      canvases[i]->Print(Form("%s]", outPdf));
    }
  }
  std::printf("plots saved in %s\n", outPdf);

  // ROOT file: canvases and histograms (from the QA task and the fractions made here)
  TString rootName = outRoot;
  if (rootName.IsNull()) {
    rootName = outPdf;
    if (rootName.EndsWith(".pdf")) {
      rootName.Resize(rootName.Length() - 4);
    }
    rootName += ".root";
  }
  std::unique_ptr<TFile> fOut(TFile::Open(rootName, "RECREATE"));
  if (!fOut || fOut->IsZombie()) {
    std::printf("cannot create %s\n", rootName.Data());
    return;
  }
  fOut->mkdir("canvases")->cd();
  for (auto* canvas : canvases) {
    canvas->Write();
  }
  fOut->mkdir("histograms")->cd();
  for (TH1* h : {static_cast<TH1*>(hBcCounter), hNumCollsInBC, hBcInTF, hBcInTFDefined, hOccTracks, hOccFT0C, static_cast<TH1*>(hOccTracksVsFT0C), hMedianTime,
                 static_cast<TH1*>(pOccVsBcInTF), static_cast<TH1*>(pOccVsBcInOrbit), hFlags, static_cast<TH1*>(hFlagsVsOcc), hBcIndexCheck}) {
    if (h) {
      h->Write();
    }
  }
  for (TH1* h : derivedHists) {
    h->Write();
  }
  std::printf("canvases and histograms saved in %s\n", rootName.Data());
}
