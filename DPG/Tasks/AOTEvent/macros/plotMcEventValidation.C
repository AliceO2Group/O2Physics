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

/// \file plotMcEventValidation.C
/// \brief Plots and summary of the MC event validation (output of o2-analysis-mc-event-validation):
///        Nch_gen vs FT0C, BC association, vertex reconstruction efficiency and purity
///
/// Usage:
///   root -l -b -q 'plotMcEventValidation.C("AnalysisResults.root")'
///   root -l -b -q 'plotMcEventValidation.C("AnalysisResults.root", "mcEventValidation.pdf", "mc-event-validation-task", 4)'
///   (the last argument: rebinning of the distance axis for the efficiency and purity vs distance to the nearest other MC collision)
/// Output: multi-page PDF, and a ROOT file with the canvases and the 2D projections (name of the PDF with .root)
///
/// \author Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include <TAttMarker.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TError.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TObject.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TStyle.h>
#include <TVirtualPad.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <memory>
#include <utility>
#include <vector>

namespace mcevval
{
// histogram from the task directory (e.g. "Reco/hNchVsFT0C"), detached from the file
template <typename T>
T* getHist(TDirectory* dir, const char* name)
{
  auto* h = dynamic_cast<T*>(dir->Get(name));
  if (!h) {
    std::printf("histogram %s not found\n", name);
    return nullptr;
  }
  h->SetDirectory(nullptr);
  h->SetStats(false);
  return h;
}

// 2D projection (x: FT0C, y: multiplicity) of a TH3 for the categories [catFirst, catLast] (bins of the z axis, from 1)
TH2* project(TH3* h3, int catFirst, int catLast, const char* name, const char* title)
{
  h3->GetZaxis()->SetRange(catFirst, catLast);
  auto* h2 = dynamic_cast<TH2*>(h3->Project3D("yx"));
  h3->GetZaxis()->SetRange(0, 0);
  if (!h2) {
    return nullptr;
  }
  h2->SetName(name);
  h2->SetTitle(title);
  h2->SetStats(false);
  return h2;
}

// draw a 2D histogram with log z in pad i
void draw2D(TCanvas* c, int i, TH2* h)
{
  if (!h) {
    return;
  }
  c->cd(i)->SetLogz();
  c->cd(i)->SetRightMargin(0.13);
  h->Draw("colz");
}

// efficiency reco / gen (Bayesian 68% intervals; bins without generated collisions are skipped)
TGraphAsymmErrors* efficiency(TH1* reco, TH1* gen, const char* name, int color)
{
  auto* eff = new TGraphAsymmErrors();
  const int errorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError; // no warning about skipped empty bins
  eff->Divide(reco, gen, "cl=0.683 b(1,1) mode");
  gErrorIgnoreLevel = errorLevel;
  eff->SetName(name);
  eff->SetLineColor(color);
  eff->SetMarkerColor(color);
  eff->SetMarkerStyle(kFullCircle);
  eff->SetMarkerSize(0.6);
  return eff;
}

// mean of a profile over all bins (weighted with the bin entries)
double profileMean(const TProfile* p)
{
  double sum = 0, n = 0;
  for (int i = 1; i <= p->GetNbinsX(); i++) {
    sum += p->GetBinContent(i) * p->GetBinEntries(i);
    n += p->GetBinEntries(i);
  }
  return n > 0 ? sum / n : -1.;
}
// profile as a graph: bins without entries are skipped
TGraphErrors* profileGraph(const TProfile* p, const char* name, int color, int marker)
{
  auto* g = new TGraphErrors();
  g->SetName(name);
  for (int i = 1; i <= p->GetNbinsX(); i++) {
    if (p->GetBinEntries(i) > 0) {
      const int n = g->GetN();
      g->SetPoint(n, p->GetBinCenter(i), p->GetBinContent(i));
      g->SetPointError(n, p->GetBinWidth(i) / 2., p->GetBinError(i));
    }
  }
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(0.6);
  return g;
}
} // namespace mcevval

void plotMcEventValidation(const char* fileName = "AnalysisResults.root", const char* outPdf = "mcEventValidation.pdf", const char* dirName = "mc-event-validation-task", int rebinDist = 4)
{
  using namespace mcevval;
  TH1::AddDirectory(false);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetPadLeftMargin(0.14);

  std::unique_ptr<TFile> f(TFile::Open(fileName, "READ"));
  if (!f || f->IsZombie()) {
    std::printf("cannot open %s\n", fileName);
    return;
  }
  TDirectory* dir = f->GetDirectory(dirName);
  if (!dir) {
    std::printf("directory %s not found in %s\n", dirName, fileName);
    return;
  }

  // MC truth
  auto* hNchVsB = getHist<TH2>(dir, "MC/hNch");
  auto* hNchFromVertexVsB = getHist<TH2>(dir, "MC/hNchFromVertex");
  auto* pFracDisplaced = getHist<TProfile>(dir, "MC/pFractionDisplacedVsB");
  auto* hNchVsNchFT0C = getHist<TH2>(dir, "MC/hNchVsNchFT0C");
  auto* hNchVsNchFT0CFromVertex = getHist<TH2>(dir, "MC/hNchVsNchFT0CFromVertex");
  auto* hDispPt = getHist<TH1>(dir, "MC/Displaced/hPt");
  auto* hDispDz = getHist<TH1>(dir, "MC/Displaced/hDeltaZ");
  auto* hDispR = getHist<TH1>(dir, "MC/Displaced/hR");
  auto* hDispSpecies = getHist<TH1>(dir, "MC/Displaced/hSpecies");
  auto* hDispGenStatus = getHist<TH1>(dir, "MC/Displaced/hGenStatus");
  auto* hVtxGenStatus = getHist<TH1>(dir, "MC/FromVertex/hGenStatus");
  auto* hNumMcCollsInBC = getHist<TH1>(dir, "MC/hNumMcCollsInBC");
  auto* hNumRecoPerMcColl = getHist<TH1>(dir, "MC/hNumRecoPerMcColl");
  // reconstructed collisions
  auto* hCounter = getHist<TH1>(dir, "Reco/hCounter");
  auto* hCategory = getHist<TH1>(dir, "Reco/hCategory");
  auto* hDeltaBC = getHist<TH1>(dir, "Reco/hDeltaBC");
  auto* hPurity = getHist<TH2>(dir, "Reco/hPurity");
  auto* h3NchVsFT0C = getHist<TH3>(dir, "Reco/hNchVsFT0C");
  auto* h3NchFromVertexVsFT0C = getHist<TH3>(dir, "Reco/hNchFromVertexVsFT0C");
  auto* h3NGlobalVsFT0C = getHist<TH3>(dir, "Reco/hNGlobalVsFT0C");
  auto* h3NGlobalVsNch = getHist<TH3>(dir, "Reco/hNGlobalVsNchFromVertex");
  // vertex reconstruction efficiency and purity
  auto* hNchGen = getHist<TH1>(dir, "Efficiency/hNchGen");
  auto* hNchGenReco = getHist<TH1>(dir, "Efficiency/hNchGenReco");
  auto* hNchGenWide = getHist<TH1>(dir, "Efficiency/hNchGenWide");
  auto* hNchGenRecoWide = getHist<TH1>(dir, "Efficiency/hNchGenRecoWide");
  auto* hDistGen = getHist<TH2>(dir, "Efficiency/hDistGen");
  auto* hDistGenReco = getHist<TH2>(dir, "Efficiency/hDistGenReco");
  auto* pPurityPVVsNch = getHist<TProfile>(dir, "Purity/pPurityPVVsNch");
  auto* pPurityAssocVsNch = getHist<TProfile>(dir, "Purity/pPurityAssocVsNch");
  auto* hPurityPVVsNch = getHist<TH2>(dir, "Purity/hPurityPVVsNch");
  auto* hPurityAssocVsNch = getHist<TH2>(dir, "Purity/hPurityAssocVsNch");
  auto* pPurityPVVsDist = getHist<TProfile2D>(dir, "Purity/pPurityPVVsDist");
  auto* pPurityAssocVsDist = getHist<TProfile2D>(dir, "Purity/pPurityAssocVsDist");

  // ### summary
  std::printf("\n=== MC event validation: %s ===\n", fileName);
  if (hNchVsB && hNchFromVertexVsB && hNchVsB->GetEntries() > 0) {
    const double meanAll = hNchVsB->GetMean(2), meanVtx = hNchFromVertexVsB->GetMean(2);
    std::printf("MC collisions: %.0f; mean Nch(|eta|<0.5): all charged physical primaries %.1f, produced at the vertex %.1f -> %.1f%% displaced\n",
                hNchVsB->GetEntries(), meanAll, meanVtx, meanAll > 0 ? 100. * (1. - meanVtx / meanAll) : 0.);
  }
  if (pFracDisplaced) {
    for (const auto& [bMin, bMax] : {std::array<double, 2>{0., 5.}, std::array<double, 2>{5., 10.}, std::array<double, 2>{10., 20.}}) {
      double sum = 0, n = 0;
      for (int i = pFracDisplaced->FindBin(bMin + 1e-6); i <= pFracDisplaced->FindBin(bMax - 1e-6); i++) {
        sum += pFracDisplaced->GetBinContent(i) * pFracDisplaced->GetBinEntries(i);
        n += pFracDisplaced->GetBinEntries(i);
      }
      std::printf("  fraction of Nch(|eta|<0.5) not produced at the vertex, %4.1f < b < %4.1f fm: %.3f\n", bMin, bMax, n > 0 ? sum / n : 0.);
    }
  }
  if (hNchVsNchFT0C && hNchVsNchFT0CFromVertex) {
    std::printf("correlation factor Nch(|eta|<0.5) vs Nch(FT0C acceptance): all primaries %.4f, produced at the vertex %.4f\n",
                hNchVsNchFT0C->GetCorrelationFactor(), hNchVsNchFT0CFromVertex->GetCorrelationFactor());
  }
  if (hDispSpecies && hDispSpecies->GetEntries() > 0) {
    std::printf("displaced charged physical primaries: e- %.1f%%, e+ %.1f%%, other %.1f%%; mean pT %.4f GeV/c\n", 100. * hDispSpecies->GetBinContent(1) / hDispSpecies->GetEntries(),
                100. * hDispSpecies->GetBinContent(2) / hDispSpecies->GetEntries(), 100. * hDispSpecies->GetBinContent(3) / hDispSpecies->GetEntries(), hDispPt ? hDispPt->GetMean() : 0.);
  }
  if (hNumMcCollsInBC && hNumMcCollsInBC->GetEntries() > 0) {
    std::printf("MC collisions sharing their BC with other MC collisions (same-bunch pile-up in MC): %.2f%%\n", 100. * (1. - hNumMcCollsInBC->GetBinContent(1) / hNumMcCollsInBC->GetEntries()));
  }
  if (hNumRecoPerMcColl && hNumRecoPerMcColl->GetEntries() > 0) {
    const double n = hNumRecoPerMcColl->GetEntries();
    std::printf("MC collisions with 0 / 1 / >=2 reconstructed collisions: %.2f%% / %.2f%% / %.2f%%\n", 100. * hNumRecoPerMcColl->GetBinContent(1) / n, 100. * hNumRecoPerMcColl->GetBinContent(2) / n,
                100. * hNumRecoPerMcColl->Integral(3, hNumRecoPerMcColl->GetNbinsX() + 1) / n);
  }
  if (hCounter) {
    std::printf("reconstructed collisions:");
    for (int i = 1; i <= hCounter->GetNbinsX(); i++) {
      std::printf(" [%s] %.0f", hCounter->GetXaxis()->GetBinLabel(i), hCounter->GetBinContent(i));
    }
    std::printf("\n");
  }
  if (hCategory && hCategory->GetEntries() > 0) {
    std::printf("selected reconstructed collisions by found BC:");
    for (int i = 1; i <= hCategory->GetNbinsX(); i++) {
      std::printf(" [%s] %.2f%%", hCategory->GetXaxis()->GetBinLabel(i), 100. * hCategory->GetBinContent(i) / hCategory->GetEntries());
    }
    std::printf("\n");
  }
  if (hPurity) {
    for (int i = 1; i <= hPurity->GetNbinsX(); i++) {
      std::unique_ptr<TH1D> p(hPurity->ProjectionY(Form("purity_%d", i), i, i));
      if (p->GetEntries() > 0) {
        std::printf("  [%s] mean fraction of PV contributors from the labelled MC collision: %.3f\n", hPurity->GetXaxis()->GetBinLabel(i), p->GetMean());
      }
    }
  }
  if (hNchGenWide && hNchGenRecoWide && hNchGenWide->GetEntries() > 0) {
    std::printf("vertex reconstruction efficiency (matched reconstructed collision): all MC collisions %.4f (%.0f / %.0f)\n", hNchGenRecoWide->GetEntries() / hNchGenWide->GetEntries(),
                hNchGenRecoWide->GetEntries(), hNchGenWide->GetEntries());
  }
  if (hNchGen && hNchGenReco && hNchGenWide && hNchGenRecoWide) {
    // classes by Nch (|eta|<0.8), as in the task: [0, 10), [10, 20), [20, 30), >= 30
    const std::array<int, 5> classEdges = {0, 10, 20, 30, -1}; // o2-linter: disable=magic-number (Nch classes, as in the task)
    for (size_t i = 0; i + 1 < classEdges.size(); i++) {
      double nGen = 0, nReco = 0;
      if (classEdges[i + 1] > 0) {
        nGen = hNchGen->Integral(hNchGen->FindBin(classEdges[i]), hNchGen->FindBin(classEdges[i + 1] - 1));
        nReco = hNchGenReco->Integral(hNchGenReco->FindBin(classEdges[i]), hNchGenReco->FindBin(classEdges[i + 1] - 1));
        std::printf("  Nch(|eta|<0.8) in [%d, %d): ", classEdges[i], classEdges[i + 1]);
      } else { // last class: all MC collisions minus those below the edge
        nGen = hNchGenWide->GetEntries() - hNchGen->Integral(1, hNchGen->FindBin(classEdges[i] - 1));
        nReco = hNchGenRecoWide->GetEntries() - hNchGenReco->Integral(1, hNchGenReco->FindBin(classEdges[i] - 1));
        std::printf("  Nch(|eta|<0.8) >= %d: ", classEdges[i]);
      }
      std::printf("efficiency %.4f (%.0f / %.0f)\n", nGen > 0 ? nReco / nGen : 0., nReco, nGen);
    }
  }
  if (hPurityPVVsNch && hPurityAssocVsNch) {
    std::printf("mean purity of reconstructed collisions wrt the labelled MC collision: PV contributors %.4f, associated tracks %.4f\n", hPurityPVVsNch->GetMean(2), hPurityAssocVsNch->GetMean(2));
  }
  if (pPurityPVVsDist && pPurityAssocVsDist) {
    for (int i = 1; i <= pPurityPVVsDist->GetNbinsY(); i++) {
      std::unique_ptr<TProfile> pPV(pPurityPVVsDist->ProfileX(Form("pPV_%d", i), i, i));
      std::unique_ptr<TProfile> pAssoc(pPurityAssocVsDist->ProfileX(Form("pAssoc_%d", i), i, i));
      std::printf("  [%s] mean purity: PV contributors %.4f, associated tracks %.4f (%.0f collisions)\n", pPurityPVVsDist->GetYaxis()->GetBinLabel(i), profileMean(pPV.get()),
                  profileMean(pAssoc.get()), pPV->GetEntries());
    }
  }

  // ### plots
  std::vector<TCanvas*> canvases;
  auto newPage = [&canvases](const char* name, const char* title, int width, int height, int nx, int ny) {
    auto* canvas = new TCanvas(name, title, width, height);
    canvas->Divide(nx, ny);
    canvases.push_back(canvas);
    return canvas;
  };
  const int nCat = h3NchVsFT0C ? h3NchVsFT0C->GetNbinsZ() : 0;
  constexpr int NCategories = 3; // categories of the found BC: correct, shared with other MC collisions, wrong

  // page 1: MC truth, Nch(|eta|<0.5) vs Nch in the FT0C acceptance and vs impact parameter
  TCanvas* c = newPage("cMcTruth", "MC truth", 1500, 1000, 2, 2);
  if (hNchVsNchFT0C) {
    hNchVsNchFT0C->SetTitle("MC: all charged physical primaries (as in MultModule)");
  }
  if (hNchVsNchFT0CFromVertex) {
    hNchVsNchFT0CFromVertex->SetTitle("MC: charged physical primaries produced at the vertex");
  }
  draw2D(c, 1, hNchVsNchFT0C);
  draw2D(c, 2, hNchVsNchFT0CFromVertex);
  draw2D(c, 3, hNchVsB);
  if (pFracDisplaced) {
    c->cd(4);
    pFracDisplaced->SetMinimum(0.);
    pFracDisplaced->SetMaximum(1.);
    pFracDisplaced->Draw();
  }

  // page 2: displaced charged physical primaries
  c = newPage("cDisplaced", "displaced charged physical primaries", 1500, 1000, 3, 2);
  int pad = 1;
  for (const auto& h : {hDispPt, hDispDz, hDispR}) {
    if (h) {
      c->cd(pad)->SetLogy();
      h->Draw("hist");
    }
    pad++;
  }
  if (hDispSpecies) {
    c->cd(4);
    hDispSpecies->SetMinimum(0.);
    hDispSpecies->Draw("hist text");
  }
  if (hDispGenStatus && hVtxGenStatus) {
    c->cd(5)->SetLogy();
    hVtxGenStatus->SetTitle("generator status of charged physical primaries;generator status code;particles");
    hVtxGenStatus->SetLineColor(kBlue + 1);
    hDispGenStatus->SetLineColor(kRed + 1);
    hVtxGenStatus->SetMinimum(0.5);
    hVtxGenStatus->SetMaximum(10. * std::max(hVtxGenStatus->GetMaximum(), hDispGenStatus->GetMaximum()));
    hVtxGenStatus->Draw("hist");
    hDispGenStatus->Draw("hist same");
    auto* leg = new TLegend(0.45, 0.75, 0.89, 0.89);
    leg->AddEntry(hVtxGenStatus, "produced at the vertex", "l");
    leg->AddEntry(hDispGenStatus, "displaced", "l");
    leg->Draw();
  }

  // page 3: reconstructed collisions, Nch_gen vs FT0C
  c = newPage("cNchVsFT0C", "Nch_gen vs FT0C", 1500, 1000, 2, 2);
  std::vector<TObject*> derived; // projections, profiles and efficiencies saved in the output file
  if (h3NchVsFT0C && h3NchFromVertexVsFT0C && nCat >= NCategories) {
    std::array<TH2*, 4> hs = {project(h3NchVsFT0C, 1, nCat, "hNchVsFT0C_all", "all primaries (as in MultModule), all collisions;FT0C amplitude (found BC);N_{ch} (|#eta|<0.5)"),
                              project(h3NchFromVertexVsFT0C, 1, nCat, "hNchFromVertexVsFT0C_all", "primaries from the vertex, all collisions;FT0C amplitude (found BC);N_{ch} (|#eta|<0.5)"),
                              project(h3NchFromVertexVsFT0C, 1, 1, "hNchFromVertexVsFT0C_correctBC", "primaries from the vertex, correct BC;FT0C amplitude (found BC);N_{ch} (|#eta|<0.5)"),
                              project(h3NchFromVertexVsFT0C, 2, 3, "hNchFromVertexVsFT0C_sharedOrWrongBC", "primaries from the vertex, shared or wrong BC;FT0C amplitude (found BC);N_{ch} (|#eta|<0.5)")};
    for (int i = 0; i < 4; i++) { // o2-linter: disable=magic-number (number of pads)
      draw2D(c, i + 1, hs[i]);
      derived.push_back(hs[i]);
    }
  }

  // page 4: reconstructed collisions, global-like tracks vs FT0C and vs Nch_gen
  c = newPage("cTracksVsFT0C", "tracks vs FT0C", 1500, 1000, 3, 2);
  if (h3NGlobalVsFT0C && nCat >= NCategories) {
    std::array<TH2*, 4> hs = {project(h3NGlobalVsFT0C, 1, nCat, "hNGlobalVsFT0C_all", "all collisions;FT0C amplitude (found BC);global-like tracks"),
                              project(h3NGlobalVsFT0C, 1, 1, "hNGlobalVsFT0C_correctBC", "correct BC;FT0C amplitude (found BC);global-like tracks"),
                              project(h3NGlobalVsFT0C, 2, 2, "hNGlobalVsFT0C_sharedBC", "BC shared with other MC collisions;FT0C amplitude (found BC);global-like tracks"),
                              project(h3NGlobalVsFT0C, 3, 3, "hNGlobalVsFT0C_wrongBC", "wrong BC;FT0C amplitude (found BC);global-like tracks")};
    for (int i = 0; i < 4; i++) { // o2-linter: disable=magic-number (number of pads)
      draw2D(c, i + 1, hs[i]);
      derived.push_back(hs[i]);
    }
  }
  if (h3NGlobalVsNch && nCat >= NCategories) {
    auto* h = project(h3NGlobalVsNch, 1, nCat, "hNGlobalVsNchFromVertex_all", "all collisions;N_{ch} (|#eta|<0.5) from the vertex;global-like tracks");
    draw2D(c, 5, h);
    derived.push_back(h);
  }

  // page 5: BC association, same-bunch pile-up in MC, vertex purity and splitting
  c = newPage("cAssociation", "BC association and pile-up", 1500, 1000, 3, 2);
  if (hCategory) {
    c->cd(1)->SetBottomMargin(0.25);
    hCategory->SetMinimum(0.);
    hCategory->Draw("hist text");
  }
  if (hDeltaBC) {
    c->cd(2)->SetLogy();
    hDeltaBC->Draw("hist");
  }
  if (hNumMcCollsInBC) {
    c->cd(3)->SetLogy();
    hNumMcCollsInBC->Draw("hist text");
  }
  if (hPurity) {
    c->cd(4)->SetLogy();
    auto* leg = new TLegend(0.15, 0.65, 0.6, 0.89);
    const std::array<int, 3> colors = {kBlue + 1, kOrange + 7, kRed + 1};
    for (int i = 1; i <= hPurity->GetNbinsX(); i++) {
      auto* p = hPurity->ProjectionY(Form("hPurity_%d", i), i, i);
      p->SetTitle("fraction of PV contributors from the labelled MC collision;purity;collisions");
      p->SetLineColor(colors[(i - 1) % 3]);
      p->SetStats(false);
      p->Draw(i == 1 ? "hist" : "hist same");
      leg->AddEntry(p, hPurity->GetXaxis()->GetBinLabel(i), "l");
      derived.push_back(p);
    }
    leg->Draw();
  }
  if (hNumRecoPerMcColl) {
    c->cd(5)->SetLogy();
    hNumRecoPerMcColl->Draw("hist text");
  }

  // page 6: vertex reconstruction efficiency
  c = newPage("cEfficiency", "vertex reconstruction efficiency", 1500, 1000, 2, 2);
  const std::array<int, 4> classColors = {kRed + 1, kOrange + 7, kGreen + 2, kBlue + 1};
  pad = 1;
  for (const auto& [reco, gen] : {std::pair<TH1*, TH1*>{hNchGenReco, hNchGen}, std::pair<TH1*, TH1*>{hNchGenRecoWide, hNchGenWide}}) {
    if (reco && gen) {
      c->cd(pad)->SetGridy();
      gPad->DrawFrame(gen->GetXaxis()->GetXmin(), 0., gen->GetXaxis()->GetXmax(), 1.1,
                      pad == 1 ? "vertex reconstruction efficiency;N_{ch} (|#eta|<0.8) from the vertex;efficiency" : "vertex reconstruction efficiency (full range);N_{ch} (|#eta|<0.8) from the vertex;efficiency");
      auto* eff = efficiency(reco, gen, pad == 1 ? "gEffVsNch" : "gEffVsNchWide", kBlack);
      eff->Draw("p");
      derived.push_back(eff);
    }
    pad++;
  }
  if (hDistGen && hDistGenReco) {
    auto* leg = new TLegend(0.55, 0.15, 0.89, 0.4);
    auto* legGen = new TLegend(0.55, 0.65, 0.89, 0.89);
    TH1* genFirst = nullptr; // first drawn class: its maximum sets the y range for all classes
    double genMax = 0;
    for (int i = 1; i <= hDistGen->GetNbinsY(); i++) {
      auto* gen = hDistGen->ProjectionX(Form("hDistGen_%d", i), i, i);
      auto* reco = hDistGenReco->ProjectionX(Form("hDistGenReco_%d", i), i, i);
      gen->Rebin(rebinDist);
      reco->Rebin(rebinDist);
      const int color = classColors[(i - 1) % classColors.size()];
      c->cd(3)->SetGridy();
      if (i == 1) {
        gPad->DrawFrame(gen->GetXaxis()->GetXmin(), 0., gen->GetXaxis()->GetXmax(), 1.1, "vertex reconstruction efficiency;distance to the nearest other MC collision (BCs);efficiency");
      }
      auto* eff = efficiency(reco, gen, Form("gEffVsDist_%d", i), color);
      eff->Draw("p");
      leg->AddEntry(eff, hDistGen->GetYaxis()->GetBinLabel(i), "lp");
      c->cd(4)->SetLogy();
      gen->SetTitle("MC collisions;distance to the nearest other MC collision (BCs);MC collisions");
      gen->SetStats(false);
      gen->SetLineColor(color);
      gen->SetMinimum(0.5);
      gen->Draw(i == 1 ? "hist" : "hist same");
      genFirst = genFirst ? genFirst : gen;
      genMax = std::max(genMax, gen->GetMaximum());
      legGen->AddEntry(gen, hDistGen->GetYaxis()->GetBinLabel(i), "l");
      derived.push_back(eff);
      derived.push_back(gen);
      derived.push_back(reco);
    }
    c->cd(3);
    leg->Draw();
    c->cd(4);
    if (genFirst) {
      genFirst->SetMaximum(5. * genMax);
    }
    legGen->Draw();
  }

  // page 7: purity of the reconstructed collisions wrt the labelled MC collision
  c = newPage("cPurity", "purity of reconstructed collisions", 1500, 1000, 3, 2);
  if (pPurityPVVsNch && pPurityAssocVsNch) {
    c->cd(1)->SetGridy();
    gPad->DrawFrame(pPurityPVVsNch->GetXaxis()->GetXmin(), 0., pPurityPVVsNch->GetXaxis()->GetXmax(), 1.1, "purity of reconstructed collisions;N_{ch} (|#eta|<0.8) of the labelled MC collision;purity");
    auto* gPV = profileGraph(pPurityPVVsNch, "gPurityPVVsNch", kBlue + 1, kFullCircle);
    auto* gAssoc = profileGraph(pPurityAssocVsNch, "gPurityAssocVsNch", kRed + 1, kOpenCircle);
    gPV->Draw("p");
    gAssoc->Draw("p");
    auto* leg = new TLegend(0.45, 0.15, 0.89, 0.3);
    leg->AddEntry(gPV, "PV contributors", "lp");
    leg->AddEntry(gAssoc, "all associated tracks", "lp");
    leg->Draw();
    derived.push_back(gPV);
    derived.push_back(gAssoc);
  }
  if (hPurityPVVsNch) {
    hPurityPVVsNch->SetTitle("purity: PV contributors;N_{ch} (|#eta|<0.8) of the labelled MC collision;purity");
  }
  if (hPurityAssocVsNch) {
    hPurityAssocVsNch->SetTitle("purity: all associated tracks;N_{ch} (|#eta|<0.8) of the labelled MC collision;purity");
  }
  draw2D(c, 2, hPurityPVVsNch);
  draw2D(c, 3, hPurityAssocVsNch);
  pad = 4;
  for (const auto& [p2, what] : {std::pair<TProfile2D*, const char*>{pPurityPVVsDist, "PV contributors"}, std::pair<TProfile2D*, const char*>{pPurityAssocVsDist, "all associated tracks"}}) {
    if (p2) {
      c->cd(pad)->SetGridy();
      gPad->DrawFrame(p2->GetXaxis()->GetXmin(), 0., p2->GetXaxis()->GetXmax(), 1.1, Form("purity: %s;distance to the nearest other MC collision (BCs);purity", what));
      auto* leg = new TLegend(0.42, 0.15, 0.89, 0.4);
      for (int i = 1; i <= p2->GetNbinsY(); i++) {
        std::unique_ptr<TProfile> p(p2->ProfileX(Form("%s_%d", p2->GetName(), i), i, i));
        p->Rebin(rebinDist);
        auto* g = profileGraph(p.get(), Form("g%s_%d", p2->GetName(), i), classColors[(i - 1) % classColors.size()], kFullCircle);
        g->Draw("p");
        leg->AddEntry(g, p2->GetYaxis()->GetBinLabel(i), "lp");
        derived.push_back(g);
      }
      leg->Draw();
    }
    pad++;
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

  // ROOT file with the canvases
  TString rootName = outPdf;
  if (rootName.EndsWith(".pdf")) {
    rootName.Resize(rootName.Length() - 4);
  }
  rootName += ".root";
  std::unique_ptr<TFile> fOut(TFile::Open(rootName, "RECREATE"));
  if (fOut && !fOut->IsZombie()) {
    fOut->mkdir("canvases")->cd();
    for (const auto& canvas : canvases) {
      canvas->Write();
    }
    fOut->mkdir("projections")->cd();
    for (const auto& h : derived) {
      h->Write();
    }
    std::printf("canvases and projections saved in %s\n", rootName.Data());
  }
}
