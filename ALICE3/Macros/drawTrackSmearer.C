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

/// \file drawTrackSmearer.C
/// \brief Draw the ALICE3 track-smearing performance curves for key particle species.

#include "FlatTrackSmearer.h"
#include "TrackUtilities.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TPDGCode.h>

namespace
{
constexpr int kElectronPdg = static_cast<int>(PDG_t::kElectron);
constexpr int kMuonPdg = static_cast<int>(PDG_t::kMuonMinus);
constexpr int kPionPdg = static_cast<int>(PDG_t::kPiPlus);
constexpr int kProtonPdg = static_cast<int>(PDG_t::kProton);
constexpr int kKaonPdg = static_cast<int>(PDG_t::kKPlus);

const std::vector<std::pair<int, std::string>> kDefaultLutFiles = {{kElectronPdg, "/tmp/lut/lutCov.acts.11.2T.dNdEta5.dat"},
                                                                   {kMuonPdg, "/tmp/lut/lutCov.acts.13.2T.dNdEta5.dat"},
                                                                   {kPionPdg, "/tmp/lut/lutCov.acts.211.2T.dNdEta5.dat"},
                                                                   {kProtonPdg, "/tmp/lut/lutCov.acts.2212.2T.dNdEta5.dat"},
                                                                   {kKaonPdg, "/tmp/lut/lutCov.acts.321.2T.dNdEta5.dat"}};
} // namespace

void drawTrackSmearer(const std::vector<std::pair<int, std::string>>& filenames = kDefaultLutFiles)
{

  o2::delphes::TrackSmearer trackSmearer;

  TCanvas* cPtReso = new TCanvas("cPtReso", "cPtReso", 800, 600);
  TCanvas* cPtEff = new TCanvas("cPtEff", "cPtEff", 800, 600);

  for (const auto& [pdg, filename] : filenames) {
    trackSmearer.loadTable(pdg, filename.c_str(), true);
    TGraph* gPt = new TGraph();
    gPt->SetName(Form("gPt_%d", pdg));
    gPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    gPt->GetYaxis()->SetTitle("pt resolution");

    TGraph* gPtEff = new TGraph();
    gPtEff->SetName(Form("gPtEff_%d", pdg));
    gPtEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    gPtEff->GetYaxis()->SetTitle("efficiency");

    for (int i = 0; i < trackSmearer.getLUTHeader(pdg)->ptmap.nbins; i++) {
      const float pt = trackSmearer.getLUTHeader(pdg)->ptmap.eval(i);
      const float res = trackSmearer.getPtRes(pdg, 0, 0., pt);
      gPt->AddPoint(pt, res / pt);
      const float eff = trackSmearer.getEfficiency(pdg, 0, 0., pt);
      gPtEff->AddPoint(pt, eff);
    }

    int color = 0;
    if (pdg == kElectronPdg) {
      color = TColor::GetColor("#e41a1c");
    } else if (pdg == kMuonPdg) {
      color = TColor::GetColor("#377eb8");
    } else if (pdg == kPionPdg) {
      color = TColor::GetColor("#4daf4a");
    } else if (pdg == kProtonPdg) {
      color = TColor::GetColor("#984ea3");
    } else if (pdg == kKaonPdg) {
      color = TColor::GetColor("#ff7f00");
    }
    gPt->SetLineColor(color);
    gPtEff->SetLineColor(color);

    cPtReso->cd();
    if (cPtReso->GetListOfPrimitives()->GetEntries() == 0) {
      gPt->Draw("ALP");
    } else {
      gPt->Draw("LP SAME");
    }
    cPtEff->cd();
    if (cPtEff->GetListOfPrimitives()->GetEntries() == 0) {
      gPtEff->Draw("ALP");
    } else {
      gPtEff->Draw("LP SAME");
    }
    gPt->SaveAs("/tmp/gPt.root");
  }
}
