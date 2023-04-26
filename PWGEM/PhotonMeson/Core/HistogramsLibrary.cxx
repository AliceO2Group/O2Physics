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
//
// Contact: daiki.sekihata@cern.ch
//

#include <iostream>
#include <memory>
#include <fstream>
using namespace std;

#include <TObject.h>
#include <TObjArray.h>
#include <THashList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <THn.h>
#include <THnSparse.h>
#include <TIterator.h>
#include <TClass.h>
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"

void o2::aod::emphotonhistograms::DefineHistograms(THashList* list, const char* histClass, const char* subGroup)
{
  if (TString(histClass) == "Event") {
    list->Add(new TH1F("hCollisionCounter", "hCollisionCounter", 5, 0.5f, 5.5f));
    list->Add(new TH1F("hZvtx_before", "vertex z; Zvtx (cm)", 100, -50, +50));
    list->Add(new TH1F("hZvtx_after", "vertex z; Zvtx (cm)", 100, -50, +50));
  }
  if (TString(histClass) == "Track") {
    list->Add(new TH1F("hPt", "pT", 1000, 0.0f, 10));
    list->Add(new TH1F("hQoverPt", "Q/pT;Q/p_{T} (GeV/c)^{-1}", 1000, -50, 50));
    list->Add(new TH2F("hEtaPhi", "#eta vs. #varphi", 180, 0, TMath::TwoPi(), 40, -2.0f, 2.0f));
    list->Add(new TH2F("hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", 100, -5.0f, 5.0f, 100, -5.0f, 5.0f));
    list->Add(new TH1F("hNclsTPC", "number of TPC clusters", 161, -0.5, 160.5));
    list->Add(new TH1F("hNcrTPC", "number of TPC crossed rows", 161, -0.5, 160.5));
    list->Add(new TH1F("hChi2TPC", "chi2/number of TPC clusters", 100, 0, 10));
    list->Add(new TH2F("hTPCdEdx", "TPC dE/dx", 1000, 0, 10, 200, 0, 200));
    list->Add(new TH2F("hTPCNsigmaEl", "TPC n sigma el", 1000, 0, 10, 100, -5, +5));
    list->Add(new TH2F("hTPCNsigmaPi", "TPC n sigma pi", 1000, 0, 10, 100, -5, +5));
    list->Add(new TH1F("hTPCNcr2Nf", "TPC Ncr/Nfindable", 200, 0, 2));
    list->Add(new TH1F("hNclsITS", "number of ITS clusters", 8, -0.5, 7.5));
    list->Add(new TH1F("hChi2ITS", "chi2/number of ITS clusters", 36, 0, 36));
  }
  if (TString(histClass) == "V0") {
    list->Add(new TH1F("hPt", "pT;p_{T} (GeV/c)", 1000, 0.0f, 10));
    list->Add(new TH2F("hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, TMath::TwoPi(), 40, -2.0f, 2.0f));
    list->Add(new TH2F("hRadius", "V0Radius; radius in Z (cm);radius in XY (cm)", 500, -250, 250, 500, 0.0f, 250.0f));
    list->Add(new TH2F("hRadius_recalc", "V0Radius; radius in Z (cm);radius in XY (cm)", 500, -250, 250, 500, 0.0f, 250.0f));
    list->Add(new TH1F("hCosPA", "V0CosPA;cosine pointing angle", 100, 0.9f, 1.0f));
    list->Add(new TH1F("hPCA", "distance between 2 legs; PCA (cm)", 100, 0.0f, 10.0f));
    list->Add(new TH2F("hAPplot", "AP plot;#alpha;q_{T} (GeV/c)", 200, -1.0f, +1.0f, 250, 0.0f, 0.25f));
    list->Add(new TH2F("hGammaPsiPair", "#psi_{pair} for photon conversion;#psi_{pair} (rad.);m_{ee} (GeV/c^{2})", 150, 0, TMath::PiOver2(), 100, 0.0f, 0.1f));
    list->Add(new TH2F("hMassGamma", "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", 200, 0.0f, 100.0f, 100, 0.0f, 0.1f));
    list->Add(new TH2F("hMassGamma_recalc", "recalc. hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", 200, 0.0f, 100.0f, 100, 0.0f, 0.1f));
    list->Add(new TH2F("hGammaRxy", "conversion point in XY;V_{x} (cm);V_{y} (cm)", 400, -100.0f, 100.0f, 400, -100.0f, 100.0f));
    list->Add(new TH2F("hGammaRxy_recalc", "recalc. conversion point in XY;V_{x} (cm);V_{y} (cm)", 400, -100.0f, 100.0f, 400, -100.0f, 100.0f));
    list->Add(new TH2F("hKFChi2vsR_recalc", "KF chi2 vs. recalc. conversion point in XY;R_{xy} (cm);KF chi2/NDF", 250, 0.0f, 250.0f, 500, 0.f, 5000.0f));
    list->Add(new TH2F("hKFChi2vsZ_recalc", "KF chi2 vs. recalc. conversion point in Z;Z (cm);KF chi2/NDF", 500, -250.0f, 250.0f, 500, 0.f, 5000.0f));
    list->Add(new TH1F("hNgamma", "Number of #gamma candidates per collision", 101, -0.5f, 100.5f));
    list->Add(new TH1F("hPt_Primary", "pT;p_{T} (GeV/c)", 1000, 0.0f, 10));                                                  // for MC efficiency
    list->Add(new TH2F("hEtaPhi_Primary", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, TMath::TwoPi(), 40, -2.0f, 2.0f)); // for MC efficiency
    list->Add(new TH1F("hPt_FromWD", "pT;p_{T} (GeV/c)", 1000, 0.0f, 10));                                                   // for MC feed down correction
    list->Add(new TH2F("hEtaPhi_FromWD", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, TMath::TwoPi(), 40, -2.0f, 2.0f));  // for MC feed down correction

    list->Add(new TH2F("hConvPoint_diffX", "conversion point diff X MC;X_{MC} (cm);X_{rec} - X_{MC} (cm)", 500, -250, +250, 100, -50.0f, 50.0f));
    list->Add(new TH2F("hConvPoint_diffY", "conversion point diff Y MC;Y_{MC} (cm);Y_{rec} - Y_{MC} (cm)", 500, -250, +250, 100, -50.0f, 50.0f));
    list->Add(new TH2F("hConvPoint_diffZ", "conversion point diff Z MC;Z_{MC} (cm);Z_{rec} - Z_{MC} (cm)", 500, -250, +250, 100, -50.0f, 50.0f));

    list->Add(new TH2F("hConvPoint_diffX_recalc", "conversion point diff X MC;X_{MC} (cm);X_{rec}^{recalc} - X_{MC} (cm)", 500, -250, +250, 100, -50.0f, 50.0f));
    list->Add(new TH2F("hConvPoint_diffY_recalc", "conversion point diff Y MC;Y_{MC} (cm);Y_{rec}^{recalc} - Y_{MC} (cm)", 500, -250, +250, 100, -50.0f, 50.0f));
    list->Add(new TH2F("hConvPoint_diffZ_recalc", "conversion point diff Z MC;Z_{MC} (cm);Z_{rec}^{recalc} - Z_{MC} (cm)", 500, -250, +250, 100, -50.0f, 50.0f));
  }

  if (TString(histClass) == "gammagamma_mass_pt") {
    // LOGF(info, "Add 2 photon histograms");
    list->Add(new TH2F("hMggPt_Same", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", 400, 0, 0.8, 400, 0.0f, 40));
    list->Add(new TH2F("hMggPt_Mixed", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", 400, 0, 0.8, 400, 0.0f, 40));
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Same"))->Sumw2();
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Mixed"))->Sumw2();
    // registry.add("EMCEMC/h2MggPt_Rotated", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/#it{c}^{2});p_{T,#gamma#gamma} (GeV/#it{c})", HistType::kTH2F, {{400, 0, 0.8}, {400, 0.0f, 40}}, true);
  }

  if (TString(histClass) == "Generated") {
    list->Add(new TH1F("hCollisionCounter", "hCollisionCounter", 5, 0.5f, 5.5f));
    list->Add(new TH1F("hZvtx_before", "vertex z; Zvtx (cm)", 100, -50, +50));
    list->Add(new TH1F("hZvtx_after", "vertex z; Zvtx (cm)", 100, -50, +50));

    if (TString(subGroup) == "ConversionStudy") {
      list->Add(new TH2F("hPhotonRxy", "conversion point in XY MC;V_{x} (cm);V_{y} (cm)", 2000, -100.0f, 100.0f, 2000, -100.0f, 100.0f));
      list->Add(new TH2F("hPhotonRZ", "conversion point in RZ MC;V_{z} (cm);R_{xy} (cm)", 5000, -250.0f, 250.0f, 1000, 0.f, 100.0f));
      list->Add(new TH2F("hPhotonPhivsRxy", "conversion point of #varphi vs. R_{xy} MC;#varphi (rad.);R_{xy} (cm);N_{e}", 360, 0.0f, TMath::TwoPi(), 900, 0, 90));
    }

    // Generated, particles
    if (TString(subGroup) == "Photon") {
      list->Add(new TH1F("hPt_Photon", "photon pT;p_{T} (GeV/c)", 1000, 0.0f, 10));
      list->Add(new TH1F("hY_Photon", "photon y;rapidity y", 40, -2.0f, 2.0f));
      list->Add(new TH1F("hPhi_Photon", "photon #varphi;#varphi (rad.)", 180, 0, TMath::TwoPi()));
    }

    ////Generated, particles
    // if (TString(subGroup) == "Pi0Eta") {
    // }
  }

  if (TString(histClass) == "tagged_photon") {
    list->Add(new TH2F("hMggPt_Same", "m_{ee#gamma} vs. p_{T,ee};m_{ee#gamma} (GeV/c^{2});p_{T,ee} (GeV/c)", 200, 0, 0.4, 100, 0.0f, 10));
    list->Add(new TH2F("hMggPt_Mixed", "m_{ee#gamma} vs. p_{T,ee};m_{ee#gamma} (GeV/c^{2});p_{T,ee} (GeV/c)", 200, 0, 0.4, 100, 0.0f, 10));
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Same"))->Sumw2();
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Mixed"))->Sumw2();
  }

  if (TString(histClass) == "photon_hbt") {
    const int ndim = 5; // qinv, qlong, qout, qside, kt
    const int nbins[ndim] = {40, 80, 80, 80, 20};
    const double xmin[ndim] = {0.0, -0.4, -0.4, -0.4, 0.0};
    const double xmax[ndim] = {0.4, +0.4, +0.4, +0.4, 1.0};

    THnSparseF* hs_q_same = new THnSparseF("hs_q_same", "hs_q_same;q_{inv} (GeV/c);q_{long} (GeV/c);q_{out} (GeV/c);q_{side} (GeV/c);k_{T} (GeV/c)", ndim, nbins, xmin, xmax);
    THnSparseF* hs_q_mix = new THnSparseF("hs_q_mix", "hs_q_mix;q_{inv} (GeV/c);q_{long} (GeV/c);q_{out} (GeV/c);q_{side} (GeV/c);k_{T} (GeV/c)", ndim, nbins, xmin, xmax);
    hs_q_same->Sumw2();
    hs_q_mix->Sumw2();
    list->Add(hs_q_same);
    list->Add(hs_q_mix);
  }
}
void o2::aod::emphotonhistograms::AddHistClass(THashList* list, const char* histClass)
{
  if (list->FindObject(histClass)) {
    std::cout << "Warning in HistogramsLibrary::AddHistClass(): Cannot add histogram class " << histClass << " because it already exists." << std::endl;
    return;
  }

  THashList* sublist = new THashList();
  sublist->SetOwner(true);
  sublist->SetName(histClass);
  list->Add(sublist);
}
