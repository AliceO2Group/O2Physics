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
#include "Framework/Logger.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"

void o2::aod::emphotonhistograms::DefineHistograms(THashList* list, const char* histClass, const char* subGroup)
{
  if (TString(histClass) == "Event") {
    list->Add(new TH1F("hCollisionCounter", "hCollisionCounter", 5, 0.5f, 5.5f));
    list->Add(new TH1F("hZvtx_before", "vertex z; Zvtx (cm)", 100, -50, +50));
    list->Add(new TH1F("hZvtx_after", "vertex z; Zvtx (cm)", 100, -50, +50));
    list->Add(new TH1F("hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", 6001, -0.5, 6000.5));
    list->Add(new TH1F("hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", 6001, -0.5, 6000.5));
    list->Add(new TH2F("hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", 300, 0, 6000, 300, 0, 6000));
    list->Add(new TH1F("hCentFT0A", "hCentFT0A;centrality FT0A (%)", 110, 0, 110));
    list->Add(new TH1F("hCentFT0C", "hCentFT0C;centrality FT0C (%)", 110, 0, 110));
    list->Add(new TH1F("hCentFT0M", "hCentFT0M;centrality FT0M (%)", 110, 0, 110));
    list->Add(new TH2F("hCentFT0MvsMultNTracksPV", "hCentFT0MvsMultNTracksPV;centrality FT0M (%);N_{track} to PV", 110, 0, 110, 600, 0, 6000));
    list->Add(new TH2F("hMultFT0MvsMultNTracksPV", "hMultFT0MvsMultNTracksPV;mult. FT0M;N_{track} to PV", 600, 0, 6000, 600, 0, 6000));
  }
  if (TString(histClass) == "V0Leg") {
    list->Add(new TH1F("hPt", "pT;p_{T} (GeV/c)", 1000, 0.0f, 10));
    list->Add(new TH1F("hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", 400, -20, 20));
    list->Add(new TH2F("hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, 2 * M_PI, 40, -2.0f, 2.0f));
    list->Add(new TH2F("hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", 200, -50.0f, 50.0f, 200, -50.0f, 50.0f));
    list->Add(new TH1F("hNclsTPC", "number of TPC clusters", 161, -0.5, 160.5));
    list->Add(new TH1F("hNcrTPC", "number of TPC crossed rows", 161, -0.5, 160.5));
    list->Add(new TH1F("hChi2TPC", "chi2/number of TPC clusters", 100, 0, 10));
    list->Add(new TH2F("hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", 1000, 0, 10, 200, 0, 200));
    list->Add(new TH2F("hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", 1000, 0, 10, 100, -5, +5));
    list->Add(new TH2F("hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", 1000, 0, 10, 100, -5, +5));
    list->Add(new TH1F("hTPCNcr2Nf", "TPC Ncr/Nfindable", 200, 0, 2));
    list->Add(new TH1F("hTPCNcls2Nf", "TPC Ncls/Nfindable", 200, 0, 2));
    list->Add(new TH1F("hNclsITS", "number of ITS clusters", 8, -0.5, 7.5));
    list->Add(new TH1F("hChi2ITS", "chi2/number of ITS clusters", 100, 0, 10));
    list->Add(new TH1F("hITSClusterMap", "ITS cluster map", 128, -0.5, 127.5));
    list->Add(new TH1F("hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", 160, 0, 16));
    list->Add(new TH2F("hXY", "X vs. Y;X (cm);Y (cm)", 100, 0, 100, 100, -50, 50));
    list->Add(new TH2F("hZX", "Z vs. X;Z (cm);X (cm)", 200, -100, 100, 100, 0, 100));
    list->Add(new TH2F("hZY", "Z vs. Y;Z (cm);Y (cm)", 200, -100, 100, 100, -50, 50));
    if (TString(subGroup) == "mc") {
      list->Add(new TH2F("hPtGen_DeltaPtOverPtGen", "electron p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", 1000, 0, 10, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hPtGen_DeltaEta", "electron #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", 1000, 0, 10, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hPtGen_DeltaPhi", "electron #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", 1000, 0, 10, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hEtaRec_DeltaPtOverPtGen", "electron p_{T} resolution;#eta^{rec} of conversion point;(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", 400, -2, +2, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hEtaRec_DeltaEta", "electron #eta resolution;#eta^{rec} of conversion point;#eta^{rec} - #eta^{gen}", 400, -2, +2, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hEtaRec_DeltaPhi", "electron #varphi resolution;#eta^{rec} of conversion point;#varphi^{rec} - #varphi^{gen} (rad.)", 400, -2, +2, 400, -1.0f, 1.0f));
    }
  }
  if (TString(histClass) == "V0") {
    list->Add(new TH1F("hPt", "pT;p_{T} (GeV/c)", 2000, 0.0f, 20));
    list->Add(new TH2F("hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, 2 * M_PI, 40, -2.0f, 2.0f));
    list->Add(new TH2F("hRadius", "V0Radius; radius in Z (cm);radius in XY (cm)", 200, -100, 100, 200, 0.0f, 100.0f));
    list->Add(new TH1F("hCosPA", "V0CosPA;cosine pointing angle", 100, 0.9f, 1.0f));
    list->Add(new TH2F("hCosPA_Rxy", "cos PA vs. R_{xy};R_{xy} (cm);cosine pointing angle", 200, 0.f, 100.f, 100, 0.9f, 1.0f));
    list->Add(new TH1F("hPCA", "distance between 2 legs;PCA (cm)", 500, 0.0f, 5.0f));
    list->Add(new TH2F("hPCA_Rxy", "distance between 2 legs vs. R_{xy};R_{xy} (cm);PCA (cm)", 200, 0.f, 100.f, 500, 0.0f, 5.0f));
    list->Add(new TH2F("hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", 200, -5.f, +5.f, 200, -5.f, +5.f));
    list->Add(new TH2F("hAPplot", "AP plot;#alpha;q_{T} (GeV/c)", 200, -1.0f, +1.0f, 250, 0.0f, 0.25f));
    list->Add(new TH2F("hMassGamma", "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", 200, 0.0f, 100.0f, 100, 0.0f, 0.1f));
    list->Add(new TH2F("hGammaRxy", "conversion point in XY;V_{x} (cm);V_{y} (cm)", 400, -100.0f, 100.0f, 400, -100.0f, 100.0f));
    list->Add(new TH2F("hKFChi2vsM", "KF chi2 vs. m_{ee};m_{ee} (GeV/c^{2});KF chi2/NDF", 100, 0.0f, 0.1f, 100, 0.f, 100.0f));
    list->Add(new TH2F("hKFChi2vsR", "KF chi2 vs. conversion point in XY;R_{xy} (cm);KF chi2/NDF", 200, 0.0f, 100.0f, 100, 0.f, 100.0f));
    list->Add(new TH2F("hKFChi2vsX", "KF chi2 vs. conversion point in X;X (cm);KF chi2/NDF", 200, -100.0f, 100.0f, 100, 0.f, 100.0f));
    list->Add(new TH2F("hKFChi2vsY", "KF chi2 vs. conversion point in Y;Y (cm);KF chi2/NDF", 200, -100.0f, 100.0f, 100, 0.f, 100.0f));
    list->Add(new TH2F("hKFChi2vsZ", "KF chi2 vs. conversion point in Z;Z (cm);KF chi2/NDF", 200, -100.0f, 100.0f, 100, 0.f, 100.0f));
    list->Add(new TH1F("hNgamma", "Number of #gamma candidates per collision", 101, -0.5f, 100.5f));

    if (TString(subGroup) == "mc") {
      list->Add(new TH1F("hPt_Photon_Candidate", "pT of photon candidate;p_{T} (GeV/c)", 2000, 0.0f, 20));                                             // for denominator of purity
      list->Add(new TH2F("hEtaPhi_Photon_Candidate", "#eta vs. #varphi of photon candidate ;#varphi (rad.);#eta", 180, 0, 2 * M_PI, 40, -2.0f, 2.0f)); // for denominator of purity

      list->Add(new TH1F("hPt_Photon_Primary", "pT;p_{T} (GeV/c)", 2000, 0.0f, 20));                                            // for MC efficiency
      list->Add(new TH2F("hEtaPhi_Photon_Primary", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, 2 * M_PI, 40, -2.0f, 2.0f)); // for MC efficiency
      list->Add(new TH2F("hXY_Photon_Primary", "X vs. Y of photon rec.;X (cm);Y (cm)", 400, -100.0f, +100, 400, -100, +100));
      list->Add(new TH2F("hXY_Photon_Primary_MC", "X vs. Y of photon gen.;X (cm);Y (cm)", 400, -100.0f, +100, 400, -100, +100));
      list->Add(new TH2F("hRZ_Photon_Primary", "R vs. Z of photon rec.;Z (cm);R_{xy} (cm)", 200, -100.0f, +100, 100, 0, 100));
      list->Add(new TH2F("hRZ_Photon_Primary_MC", "R vs. Z of photon gen;Z (cm);R_{xy} (cm)", 200, -100.0f, +100, 100, 0, 100));

      list->Add(new TH1F("hPt_Photon_FromWD", "pT;p_{T} (GeV/c)", 2000, 0.0f, 20));                                            // for MC feed down correction
      list->Add(new TH2F("hEtaPhi_Photon_FromWD", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, 2 * M_PI, 40, -2.0f, 2.0f)); // for MC feed down correction

      list->Add(new TH2F("hRZ_Photon_hs", "R vs. Z of photon from hadronic shower in materials;Z (cm);R_{xy} (cm)", 200, -100.0f, +100, 100, 0, 100));
      list->Add(new TH1F("hPt_Photon_hs", "pT of photon from hadronic shower in materials;p_{T} (GeV/c)", 1000, 0.0f, 10));
      list->Add(new TH1F("hEta_Photon_hs", "rapidity of photon from hadronic shower in materials;y", 40, -2.0f, 2.0f));
      list->Add(new TH1F("hPhi_Photon_hs", "#varphi of photon from hadronic shower in materials;#varphi (rad.);", 72, 0, 2 * M_PI));

      list->Add(new TH2F("hConvPoint_diffX", "conversion point diff X MC;X_{MC} (cm);X_{rec} - X_{MC} (cm)", 200, -100, +100, 100, -50.0f, 50.0f));
      list->Add(new TH2F("hConvPoint_diffY", "conversion point diff Y MC;Y_{MC} (cm);Y_{rec} - Y_{MC} (cm)", 200, -100, +100, 100, -50.0f, 50.0f));
      list->Add(new TH2F("hConvPoint_diffZ", "conversion point diff Z MC;Z_{MC} (cm);Z_{rec} - Z_{MC} (cm)", 200, -100, +100, 100, -50.0f, 50.0f));

      list->Add(new TH2F("hPtGen_DeltaPtOverPtGen", "photon p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", 1000, 0, 10, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hPtGen_DeltaEta", "photon #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", 1000, 0, 10, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hPtGen_DeltaPhi", "photon #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", 1000, 0, 10, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hEtaRec_DeltaPtOverPtGen", "photon p_{T} resolution;#eta^{rec} of conversion point;(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", 400, -2, +2, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hEtaRec_DeltaEta", "photon #eta resolution;#eta^{rec} of conversion point;#eta^{rec} - #eta^{gen}", 400, -2, +2, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hEtaRec_DeltaPhi", "photon #varphi resolution;#eta^{rec} of conversion point;#varphi^{rec} - #varphi^{gen} (rad.)", 400, -2, +2, 400, -1.0f, 1.0f));
    } // end of mc
  }   // end of V0

  if (TString(histClass).Contains("Dalitz")) {
    THnSparseF* hs_dilepton_uls_same = nullptr;
    THnSparseF* hs_dilepton_lspp_same = nullptr;
    THnSparseF* hs_dilepton_lsmm_same = nullptr;
    THnSparseF* hs_dilepton_uls_dca_same = nullptr;
    THnSparseF* hs_dilepton_lspp_dca_same = nullptr;
    THnSparseF* hs_dilepton_lsmm_dca_same = nullptr;

    const int ndca = 27;
    double dca[ndca] = {0.f};
    for (int i = 0; i < 20; i++) {
      dca[i] = 0.1 * i;
    }
    for (int i = 20; i < ndca; i++) {
      dca[i] = 0.5 * (i - 20) + 2.0;
    }

    if (TString(histClass).Contains("EE")) {
      const int nm = 167;
      double mee[nm] = {0.f};
      for (int i = 0; i < 110; i++) {
        mee[i] = 0.01 * (i - 0) + 0.0; // every 0.01 GeV/c2 up to 1.1 GeV/c2
      }
      for (int i = 110; i < 128; i++) {
        mee[i] = 0.1 * (i - 110) + 1.1; // every 0.1 GeV/c2 from 1.1 to 2.9 GeV/c2
      }
      for (int i = 128; i < 158; i++) {
        mee[i] = 0.01 * (i - 128) + 2.9; // every 0.01 GeV/c2 from 2.9 to 3.2 GeV/c2
      }
      for (int i = 158; i < nm; i++) {
        mee[i] = 0.1 * (i - 158) + 3.2; // every 0.01 GeV/c2 from 3.2 to 3.5 GeV/c2
      }

      const int npt = 61;
      double pt[npt] = {0.f};
      for (int i = 0; i < 50; i++) {
        pt[i] = 0.1 * i;
      }
      for (int i = 50; i < npt; i++) {
        pt[i] = 0.5 * (i - 50) + 5.0;
      }

      const int ndim = 4; // m, pt, dca, phiv
      const int nbins[ndim] = {nm - 1, npt - 1, ndca - 1, 18};
      const double xmin[ndim] = {0.0, 0.0, 0.0, 0.0};
      const double xmax[ndim] = {4.0, 10.0, 5.0, M_PI};

      hs_dilepton_uls_same = new THnSparseF("hs_dilepton_uls_same", "hs_dilepton_uls;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);DCA_{ee}^{3D} (#sigma);#varphi_{V} (rad.);", ndim, nbins, xmin, xmax);
      hs_dilepton_uls_same->SetBinEdges(0, mee);
      hs_dilepton_uls_same->SetBinEdges(1, pt);
      hs_dilepton_uls_same->SetBinEdges(2, dca);
      hs_dilepton_uls_same->Sumw2();
      list->Add(hs_dilepton_uls_same);

      hs_dilepton_lspp_same = new THnSparseF("hs_dilepton_lspp_same", "hs_dilepton_lspp;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);DCA_{ee}^{3D} (#sigma);#varphi_{V} (rad.);", ndim, nbins, xmin, xmax);
      hs_dilepton_lspp_same->SetBinEdges(0, mee);
      hs_dilepton_lspp_same->SetBinEdges(1, pt);
      hs_dilepton_lspp_same->SetBinEdges(2, dca);
      hs_dilepton_lspp_same->Sumw2();
      list->Add(hs_dilepton_lspp_same);

      hs_dilepton_lsmm_same = new THnSparseF("hs_dilepton_lsmm_same", "hs_dilepton_lsmm;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);DCA_{ee}^{3D} (#sigma);#varphi_{V} (rad.);", ndim, nbins, xmin, xmax);
      hs_dilepton_lsmm_same->SetBinEdges(0, mee);
      hs_dilepton_lsmm_same->SetBinEdges(1, pt);
      hs_dilepton_lsmm_same->SetBinEdges(2, dca);
      hs_dilepton_lsmm_same->Sumw2();
      list->Add(hs_dilepton_lsmm_same);

      if (TString(subGroup).Contains("mix")) {
        THnSparseF* hs_dilepton_uls_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_uls_same->Clone("hs_dilepton_uls_mix"));
        THnSparseF* hs_dilepton_lspp_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_lspp_same->Clone("hs_dilepton_lspp_mix"));
        THnSparseF* hs_dilepton_lsmm_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_lsmm_same->Clone("hs_dilepton_lsmm_mix"));
        list->Add(hs_dilepton_uls_mix);
        list->Add(hs_dilepton_lspp_mix);
        list->Add(hs_dilepton_lsmm_mix);
      }

      if (TString(subGroup).Contains("dca")) {
        const int ndim_dca = 4; // mee, dca1, dca2, dca12
        const int nbins_dca[ndim_dca] = {nm - 1, ndca - 1, ndca - 1, ndca - 1};
        const double xmin_dca[ndim_dca] = {0.0, 0.0, 0.0, 0.0};
        const double xmax_dca[ndim_dca] = {4.0, 5.0, 5.0, 5.0};

        hs_dilepton_uls_dca_same = new THnSparseF("hs_dilepton_uls_dca_same", "hs_dilepton_uls_dca;m_{ee} (GeV/c^{2});DCA_{e^{+}}^{3D} (#sigma);DCA_{e^{-}}^{3D} (#sigma);DCA_{ee}^{3D} (#sigma);", ndim_dca, nbins_dca, xmin_dca, xmax_dca);
        hs_dilepton_uls_dca_same->SetBinEdges(0, mee);
        hs_dilepton_uls_dca_same->SetBinEdges(1, dca);
        hs_dilepton_uls_dca_same->SetBinEdges(2, dca);
        hs_dilepton_uls_dca_same->SetBinEdges(3, dca);
        hs_dilepton_uls_dca_same->Sumw2();
        list->Add(hs_dilepton_uls_dca_same);

        hs_dilepton_lspp_dca_same = reinterpret_cast<THnSparseF*>(hs_dilepton_uls_dca_same->Clone("hs_dilepton_lspp_dca_same"));
        hs_dilepton_lsmm_dca_same = reinterpret_cast<THnSparseF*>(hs_dilepton_uls_dca_same->Clone("hs_dilepton_lsmm_dca_same"));
        list->Add(hs_dilepton_lspp_dca_same);
        list->Add(hs_dilepton_lsmm_dca_same);

        THnSparseF* hs_dilepton_uls_dca_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_uls_dca_same->Clone("hs_dilepton_uls_dca_mix"));
        THnSparseF* hs_dilepton_lspp_dca_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_lspp_dca_same->Clone("hs_dilepton_lspp_dca_mix"));
        THnSparseF* hs_dilepton_lsmm_dca_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_lsmm_dca_same->Clone("hs_dilepton_lsmm_dca_mix"));
        list->Add(hs_dilepton_uls_dca_mix);
        list->Add(hs_dilepton_lspp_dca_mix);
        list->Add(hs_dilepton_lsmm_dca_mix);
      }

      if (TString(subGroup).Contains("mc")) {
        // create phiv template
        list->Add(new TH2F("hMvsPhiV_Pi0", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", 72, 0, M_PI, 100, 0.0f, 0.1f));    // ee from pi0 dalitz decay
        list->Add(new TH2F("hMvsPhiV_Eta", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", 72, 0, M_PI, 100, 0.0f, 0.1f));    // ee from eta dalitz decay
        list->Add(new TH2F("hMvsPhiV_Photon", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", 72, 0, M_PI, 100, 0.0f, 0.1f)); // ee from photon conversion

        list->Add(new TH2F("hMvsOPA_Pi0", "m_{ee} vs. opening angle;opening angle (rad.);m_{ee} (GeV/c^{2})", 100, 0, 0.1, 100, 0.0f, 0.1f));    // ee from pi0 dalitz decay
        list->Add(new TH2F("hMvsOPA_Eta", "m_{ee} vs. opening angle;opening angle (rad.);m_{ee} (GeV/c^{2})", 100, 0, 0.1, 100, 0.0f, 0.1f));    // ee from eta dalitz decay
        list->Add(new TH2F("hMvsOPA_Photon", "m_{ee} vs. opening angle;opening angle (rad.);m_{ee} (GeV/c^{2})", 100, 0, 0.1, 100, 0.0f, 0.1f)); // ee from photon conversion
      }                                                                                                                                          // end of mc
    } else if (TString(histClass).Contains("MuMu")) {
      const int ndim = 4; // m, pt, dca, phiv
      const int nbins[ndim] = {90, 20, ndca - 1, 1};
      const double xmin[ndim] = {0.2, 0.0, 0.0, 0.0};
      const double xmax[ndim] = {1.1, 2.0, 5.0, 3.2};

      hs_dilepton_uls_same = new THnSparseF("hs_dilepton_uls_same", "hs_dilepton_uls;m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c);DCA_{#mu#mu}^{3D} (#sigma);#varphi_{V} (rad.);", ndim, nbins, xmin, xmax);
      hs_dilepton_uls_same->Sumw2();
      hs_dilepton_uls_same->SetBinEdges(2, dca);
      list->Add(hs_dilepton_uls_same);

      hs_dilepton_lspp_same = new THnSparseF("hs_dilepton_lspp_same", "hs_dilepton_lspp;m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c);DCA_{#mu#mu}^{3D} (#sigma);#varphi_{V} (rad.);", ndim, nbins, xmin, xmax);
      hs_dilepton_lspp_same->Sumw2();
      hs_dilepton_lspp_same->SetBinEdges(2, dca);
      list->Add(hs_dilepton_lspp_same);

      hs_dilepton_lsmm_same = new THnSparseF("hs_dilepton_lsmm_same", "hs_dilepton_lsmm;m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c);DCA_{#mu#mu}^{3D} (#sigma);#varphi_{V} (rad.);", ndim, nbins, xmin, xmax);
      hs_dilepton_lsmm_same->Sumw2();
      hs_dilepton_lsmm_same->SetBinEdges(2, dca);
      list->Add(hs_dilepton_lsmm_same);

      if (TString(subGroup).Contains("mix")) {
        THnSparseF* hs_dilepton_uls_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_uls_same->Clone("hs_dilepton_uls_mix"));
        THnSparseF* hs_dilepton_lspp_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_lspp_same->Clone("hs_dilepton_lspp_mix"));
        THnSparseF* hs_dilepton_lsmm_mix = reinterpret_cast<THnSparseF*>(hs_dilepton_lsmm_same->Clone("hs_dilepton_lsmm_mix"));
        list->Add(hs_dilepton_uls_mix);
        list->Add(hs_dilepton_lspp_mix);
        list->Add(hs_dilepton_lsmm_mix);
      }
    } else {
      LOGF(info, "EE or MuMu are supported.");
    }

    list->Add(new TH1F("hNpair_uls", "Number of ULS pairs per collision", 101, -0.5f, 100.5f));
    list->Add(new TH1F("hNpair_lspp", "Number of LS++ pairs per collision", 101, -0.5f, 100.5f));
    list->Add(new TH1F("hNpair_lsmm", "Number of LS-- pairs per collision", 101, -0.5f, 100.5f));

  } // end of Dalitz
  if (TString(histClass) == "Track") {
    float maxP = 10.f;
    if (TString(subGroup).Contains("Mu")) {
      maxP = 1.f;
    }

    list->Add(new TH1F("hPt", "pT;p_{T} (GeV/c)", 1000, 0.0f, maxP));
    list->Add(new TH1F("hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", 400, -20, 20));
    list->Add(new TH2F("hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, 2 * M_PI, 40, -2.0f, 2.0f));
    list->Add(new TH2F("hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", 200, -1.0f, 1.0f, 200, -1.0f, 1.0f));
    list->Add(new TH2F("hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", 200, -10.0f, 10.0f, 200, -10.0f, 10.0f));
    list->Add(new TH2F("hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", 1000, 0, maxP, 100, 0., 1000));
    list->Add(new TH2F("hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", 1000, 0, maxP, 100, 0., 1000));
    list->Add(new TH1F("hNclsTPC", "number of TPC clusters", 161, -0.5, 160.5));
    list->Add(new TH1F("hNcrTPC", "number of TPC crossed rows", 161, -0.5, 160.5));
    list->Add(new TH1F("hChi2TPC", "chi2/number of TPC clusters", 100, 0, 10));
    list->Add(new TH2F("hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", 1000, 0, maxP, 200, 0, 200));
    list->Add(new TH2F("hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTOFbeta", "TOF beta;p_{in} (GeV/c);TOF #beta", 1000, 0, maxP, 600, 0, 1.2));
    list->Add(new TH2F("hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTOFNsigmaMu", "TOF n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH2F("hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", 1000, 0, maxP, 100, -5, +5));
    list->Add(new TH1F("hTPCNcr2Nf", "TPC Ncr/Nfindable", 200, 0, 2));
    list->Add(new TH1F("hTPCNcls2Nf", "TPC Ncls/Nfindable", 200, 0, 2));
    list->Add(new TH1F("hNclsITS", "number of ITS clusters", 8, -0.5, 7.5));
    list->Add(new TH1F("hChi2ITS", "chi2/number of ITS clusters", 100, 0, 10));
    list->Add(new TH1F("hITSClusterMap", "ITS cluster map", 128, -0.5, 127.5));
    list->Add(new TH1F("hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", 160, 0, 16));
    if (TString(subGroup).Contains("mc")) {
      list->Add(new TH2F("hPtGen_DeltaPtOverPtGen", "electron p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", 1000, 0, maxP, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hPtGen_DeltaEta", "electron #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", 1000, 0, maxP, 400, -1.0f, 1.0f));
      list->Add(new TH2F("hPtGen_DeltaPhi", "electron #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", 1000, 0, maxP, 400, -1.0f, 1.0f));
    }
  }

  if (TString(histClass) == "Cluster") {
    list->Add(new TH1F("hPt", "pT;p_{T} (GeV/c)", 1000, 0.0f, 10));
    list->Add(new TH2F("hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", 180, 0, 2 * M_PI, 400, -2.0f, 2.0f));
    list->Add(new TH1F("hNgamma", "Number of #gamma candidates per collision", 101, -0.5f, 100.5f));

    if (TString(subGroup) == "PHOS") {
      list->Add(new TH2F("hEvsNcell", "E_{cluster} vs. N_{cell};E_{cluster} (GeV);N_{cell}", 200, 0, 20, 51, -0.5, 50.5f));
      list->Add(new TH2F("hEvsM02", "E_{cluster} vs. M02;E_{cluster} (GeV);M02 (cm)", 200, 0, 20, 100, 0, 10));
      list->Add(new TH2F("hEvsM20", "E_{cluster} vs. M20;E_{cluster} (GeV);M20 (cm)", 200, 0, 20, 100, 0, 10));
      list->Add(new TH1F("hDistToBC", "distance to bad channel", 100, 0, 10));

      const int nmod = 4;
      for (int i = 1; i <= nmod; i++) {
        list->Add(new TH2F(Form("hClusterXZM%d", i), Form("cluster (X,Z) on M%d;X;Z", i), 64, 0.5, 64.5, 56, 0.5, 56.5));
      } // end of module loop
    }
  }

  if (TString(histClass) == "singlephoton") {
    list->Add(new TH1F("hPt", "pT of photon candidates;p_{T} (GeV/c)", 2000, 0.0f, 20));
    list->Add(new TH1F("hY", "rapidity of photon candidates;y", 40, -2.0f, 2.0f));
    list->Add(new TH1F("hPhi", "azimuthal angle of photon candidates;#varphi (rad.)", 180, 0, 2 * M_PI));
    if (TString(subGroup) == "mc") {
      list->Add(new TH1F("hPt_Photon_Primary", "pT;p_{T} (GeV/c)", 2000, 0.0f, 20));           // for MC efficiency
      list->Add(new TH1F("hY_Photon_Primary", "rapidity;y", 40, -2.0f, 2.0f));                 // for MC efficiency
      list->Add(new TH1F("hPhi_Photon_Primary", "#varphi;#varphi (rad.);", 180, 0, 2 * M_PI)); // for MC efficiency
      list->Add(new TH1F("hPt_Photon_FromWD", "pT;p_{T} (GeV/c)", 2000, 0.0f, 20));            // for Feed down correction
      list->Add(new TH1F("hY_Photon_FromWD", "rapidity;y", 40, -2.0f, 2.0f));                  // for Feed down correction
      list->Add(new TH1F("hPhi_Photon_FromWD", "#varphi;#varphi (rad.);", 180, 0, 2 * M_PI));  // for Feed down correction
      list->Add(new TH1F("hPt_Photon_hs", "pT of photon from hadronic shower in materials;p_{T} (GeV/c)", 2000, 0.0f, 20));
      list->Add(new TH1F("hY_Photon_hs", "rapidity from hadronic shower in materials;y", 40, -2.0f, 2.0f));
      list->Add(new TH1F("hPhi_Photon_hs", "#varphi from hadronic shower in materials;#varphi (rad.);", 180, 0, 2 * M_PI));
    }
  }

  const int nmgg = 401;
  float mgg[nmgg] = {};
  for (int i = 0; i < nmgg; i++) {
    mgg[i] = 0.002 * i;
  }
  const int npTgg = 71;
  float pTgg[npTgg] = {};
  for (int i = 0; i < 50; i++) {
    pTgg[i] = 0.1 * (i - 0) + 0.0; // from 0 to 5 GeV/c, every 0.1 GeV/c
  }
  for (int i = 50; i < 60; i++) {
    pTgg[i] = 0.5 * (i - 50) + 5.0; // from 5 to 10 GeV/c, evety 0.5 GeV/c
  }
  for (int i = 60; i < npTgg; i++) {
    pTgg[i] = 1.0 * (i - 60) + 10.0; // from 10 to 20 GeV/c, evety 1 GeV/c
  }
  if (TString(histClass) == "gammagamma_mass_pt") {
    list->Add(new TH2F("hMggPt_Same", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", nmgg - 1, mgg, npTgg - 1, pTgg));
    list->Add(new TH2F("hMggPt_Mixed", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", nmgg - 1, mgg, npTgg - 1, pTgg));
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Same"))->Sumw2();
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Mixed"))->Sumw2();

    if (TString(subGroup) == "PCMPHOS" || TString(subGroup) == "PCMEMC") {
      list->Add(new TH2F("hdEtadPhi", "#Delta#eta vs. #Delta#varphi;#Delta#varphi (rad.);#Delta#eta", 200, -1, +1, 200, -1, +1));
      list->Add(new TH2F("hdEtaPt", "#Delta#eta vs. p_{T}^{leg};p_{T}^{leg} (GeV/c);#Delta#eta", 100, 0, 10, 200, -1, +1));
      list->Add(new TH2F("hdPhiPt", "#Delta#varphi vs. p_{T}^{leg};p_{T}^{leg} (GeV/c);#Delta#varphi", 100, 0, 10, 200, -1, +1));
      list->Add(new TH2F("hEp_E", "E/p vs. matched E_{cluster};E_{cluster} (GeV);E/p", 100, 0, 10, 200, 0, 2));
    }

    if (TString(subGroup) == "EMCEMC") {
      list->Add(new TH2F("hMggPt_Same_RotatedBkg", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", nmgg - 1, mgg, npTgg - 1, pTgg));
      reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Same_RotatedBkg"))->Sumw2();
    }
  }
  if (TString(histClass) == "gammagamma_mass_pt_mc") {
    list->Add(new TH2F("hMggPt_Pi0_Primary", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", nmgg - 1, mgg, npTgg - 1, pTgg));
    list->Add(new TH2F("hMggPt_Pi0_FromWD", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", nmgg - 1, mgg, npTgg - 1, pTgg));
    list->Add(new TH2F("hMggPt_Eta_Primary", "m_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", nmgg - 1, mgg, npTgg - 1, pTgg));
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Pi0_Primary"))->Sumw2();
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Pi0_FromWD"))->Sumw2();
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Eta_Primary"))->Sumw2();
  }

  if (TString(histClass) == "material_budget_study") {
    const int npt = 70;
    double pt[npt] = {0.f};
    for (int i = 0; i < 10; i++) {
      pt[i] = 0.01 * (i - 0) + 0.0;
    }
    for (int i = 10; i < 59; i++) {
      pt[i] = 0.1 * (i - 10) + 0.1;
    }
    for (int i = 59; i < npt; i++) {
      pt[i] = 0.5 * (i - 59) + 5.0;
    }
    if (TString(subGroup) == "V0") {
      const int ndim = 4; // pt, r, phi, eta
      const int nbins[ndim] = {npt - 1, 90, 72, 40};
      const double xmin[ndim] = {0.0, 0.0, 0.0, -2.0};
      const double xmax[ndim] = {10.0, 90.0, 2 * M_PI, +2.0};
      THnSparseF* hs_conv_point = new THnSparseF("hs_conv_point", "hs_conv_point;p_{T,#gamma} (GeV/c);R_{xy} (cm);#varphi (rad.);#eta;", ndim, nbins, xmin, xmax);
      hs_conv_point->SetBinEdges(0, pt);
      hs_conv_point->Sumw2();
      list->Add(hs_conv_point);
    } else if (TString(subGroup) == "Pair") {
      const int ndim = 6; // mgg, pT1, pT2, rxy2, eta2, phi2
      const int nbins[ndim] = {200, npt - 1, npt - 1, 90, 72, 40};
      const double xmin[ndim] = {0.0, 0.0, 0.0, 0, 0, -2};
      const double xmax[ndim] = {0.4, 10.0, 10.0, 90.0, 2 * M_PI, +2};

      THnSparseF* hs_conv_point_same = new THnSparseF("hs_conv_point_same", "hs_conv_point;m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma}^{tag} (GeV/c);p_{T,#gamma}^{probe} (GeV/c);R_{xy}^{probe} (cm);#varphi^{probe} (rad.);#eta^{probe};", ndim, nbins, xmin, xmax);
      hs_conv_point_same->SetBinEdges(1, pt);
      hs_conv_point_same->SetBinEdges(2, pt);
      hs_conv_point_same->Sumw2();
      list->Add(hs_conv_point_same);

      THnSparseF* hs_conv_point_mix = new THnSparseF("hs_conv_point_mix", "hs_conv_point;m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma}^{tag} (GeV/c);p_{T,#gamma}^{probe} (GeV/c);R_{xy}^{probe} (cm);#varphi^{probe} (rad.);#eta^{probe};", ndim, nbins, xmin, xmax);
      hs_conv_point_mix->SetBinEdges(1, pt);
      hs_conv_point_mix->SetBinEdges(2, pt);
      hs_conv_point_mix->Sumw2();
      list->Add(hs_conv_point_mix);
    } // end of pair
  }   // end of material budget study

  if (TString(histClass) == "Generated") {
    list->Add(new TH1F("hCollisionCounter", "hCollisionCounter", 5, 0.5f, 5.5f));
    list->Add(new TH1F("hZvtx_before", "vertex z; Zvtx (cm)", 100, -50, +50));
    list->Add(new TH1F("hZvtx_after", "vertex z; Zvtx (cm)", 100, -50, +50));
    list->Add(new TH1F("hNrecPerMCCollision", "Nrec per mc collisions;N_{rec} collisions per MC collisions", 101, -0.5f, 100.5f));

    if (TString(subGroup) == "Photon") {
      list->Add(new TH1F("hPt_Photon", "photon pT;p_{T} (GeV/c)", 2000, 0.0f, 20));
      list->Add(new TH1F("hY_Photon", "photon y;rapidity y", 40, -2.0f, 2.0f));
      list->Add(new TH1F("hPhi_Photon", "photon #varphi;#varphi (rad.)", 180, 0, 2 * M_PI));
      list->Add(new TH1F("hPt_ConvertedPhoton", "converted photon pT;p_{T} (GeV/c)", 2000, 0.0f, 20));
      list->Add(new TH1F("hY_ConvertedPhoton", "converted photon y;rapidity y", 40, -2.0f, 2.0f));
      list->Add(new TH1F("hPhi_ConvertedPhoton", "converted photon #varphi;#varphi (rad.)", 180, 0, 2 * M_PI));
      list->Add(new TH2F("hPhotonRxy", "conversion point in XY MC;V_{x} (cm);V_{y} (cm)", 2000, -100.0f, 100.0f, 2000, -100.0f, 100.0f));
      list->Add(new TH2F("hPhotonRZ", "conversion point in RZ MC;V_{z} (cm);R_{xy} (cm)", 2000, -100.0f, 100.0f, 1000, 0.f, 100.0f));
      list->Add(new TH2F("hPhotonPhivsRxy", "conversion point of #varphi vs. R_{xy} MC;#varphi (rad.);R_{xy} (cm);N_{e}", 360, 0.0f, 2 * M_PI, 100, 0, 100));
    }

    if (TString(subGroup) == "Pi0Eta") {
      static constexpr std::string_view parnames[2] = {"Pi0", "Eta"};
      for (int i = 0; i < 2; i++) {
        list->Add(new TH1F(Form("hPt_%s", parnames[i].data()), Form("%s pT;p_{T} (GeV/c)", parnames[i].data()), 2000, 0.0f, 20));
        list->Add(new TH1F(Form("hY_%s", parnames[i].data()), Form("%s y;rapidity y", parnames[i].data()), 40, -2.0f, 2.0f));
        list->Add(new TH1F(Form("hPhi_%s", parnames[i].data()), Form("%s #varphi;#varphi (rad.)", parnames[i].data()), 180, 0, 2 * M_PI));
      }
    }
    if (TString(subGroup) == "dielectron") {
      TH2F* hMvsPt = new TH2F("hMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", 110, 0, 1.1f, 1000, 0, 10.f);
      hMvsPt->Sumw2();
      list->Add(hMvsPt);
    } else if (TString(subGroup) == "dimuon") {
      TH2F* hMvsPt = new TH2F("hMvsPt", "m_{#mu#mu} vs. p_{T,#mu#mu};m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c)", 90, 0.2, 1.1f, 200, 0, 2.f);
      hMvsPt->Sumw2();
      list->Add(hMvsPt);
    }
  }

  const int nmgg04 = 201;
  float mgg04[nmgg04] = {};
  for (int i = 0; i < nmgg04; i++) {
    mgg04[i] = 0.002 * i;
  }
  const int npTgg10 = 61;
  float pTgg10[npTgg10] = {};
  for (int i = 0; i < 50; i++) {
    pTgg10[i] = 0.1 * (i - 0) + 0.0; // from 0 to 5 GeV/c, every 0.1 GeV/c
  }
  for (int i = 50; i < npTgg10; i++) {
    pTgg10[i] = 0.5 * (i - 50) + 5.0; // from 5 to 10 GeV/c, evety 0.5 GeV/c
  }
  if (TString(histClass) == "tagging_pi0") {
    list->Add(new TH2F("hMggPt_Same", "m_{ee#gamma} vs. p_{T,#gamma};m_{ee#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c)", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10));
    list->Add(new TH2F("hMggPt_Mixed", "m_{ee#gamma} vs. p_{T,#gamma};m_{ee#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c)", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10));
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Same"))->Sumw2();
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Mixed"))->Sumw2();
  }
  if (TString(histClass) == "tagging_pi0_mc") {
    if (TString(subGroup) == "pcm") {
      list->Add(new TH1F("hPt_v0photon_Pi0_Primary", "reconstcuted v0 photon from primary #pi^{0};p_{T,#gamma} (GeV/c);N_{#gamma}^{#pi^{0}}", npTgg10 - 1, pTgg10)); // denominator for conditional probability
      reinterpret_cast<TH1F*>(list->FindObject("hPt_v0photon_Pi0_Primary"))->Sumw2();
      list->Add(new TH1F("hPt_v0photon_Pi0_FromWD", "reconstcuted v0 photon from #pi^{0} from WD;p_{T,#gamma} (GeV/c);N_{#gamma}^{#pi^{0}}", npTgg10 - 1, pTgg10)); // denominator for conditional probability
      reinterpret_cast<TH1F*>(list->FindObject("hPt_v0photon_Pi0_FromWD"))->Sumw2();
      list->Add(new TH1F("hPt_v0photon_Pi0_hs", "reconstcuted v0 photon from #pi^{0} from hadronic shower in materials;p_{T,#gamma} (GeV/c);N_{#gamma}^{#pi^{0}}", npTgg10 - 1, pTgg10)); // denominator for conditional probability
      reinterpret_cast<TH1F*>(list->FindObject("hPt_v0photon_Pi0_hs"))->Sumw2();
    } else if (TString(subGroup) == "pair") {
      list->Add(new TH2F("hMggPt_Pi0_Primary", "reconstructed m_{ee#gamma} vs. p_{T,#gamma} from primary #pi^{0};m_{ee#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c);N_{#gamma}^{tagged #pi^{0}}", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10)); // numerator for conditional probability
      reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Pi0_Primary"))->Sumw2();
      list->Add(new TH2F("hMggPt_Pi0_FromWD", "reconstructed m_{ee#gamma} vs. p_{T,#gamma} from #pi^{0} from WD;m_{ee#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c);N_{#gamma}^{tagged #pi^{0}}", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10)); // numerator for conditional probability
      reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Pi0_FromWD"))->Sumw2();
      list->Add(new TH2F("hMggPt_Pi0_hs", "reconstructed m_{ee#gamma} vs. p_{T,#gamma} from #pi^{0} from hadronic shower in material;m_{ee#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c);N_{#gamma}^{tagged #pi^{0}}", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10)); // numerator for conditional probability
      reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Pi0_hs"))->Sumw2();
    }
  }

  if (TString(histClass) == "tag_and_probe") {
    list->Add(new TH2F("hMggPt_Probe_Same", "m_{#gamma#gamma} vs. p_{T,#gamma};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c)", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10));
    list->Add(new TH2F("hMggPt_Probe_Mixed", "m_{#gamma#gamma} vs. p_{T,#gamma};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c)", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10));
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Probe_Same"))->Sumw2();
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_Probe_Mixed"))->Sumw2();
    list->Add(new TH2F("hMggPt_PassingProbe_Same", "m_{#gamma#gamma} vs. p_{T,#gamma};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c)", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10));
    list->Add(new TH2F("hMggPt_PassingProbe_Mixed", "m_{#gamma#gamma} vs. p_{T,#gamma};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c)", nmgg04 - 1, mgg04, npTgg10 - 1, pTgg10));
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_PassingProbe_Same"))->Sumw2();
    reinterpret_cast<TH2F*>(list->FindObject("hMggPt_PassingProbe_Mixed"))->Sumw2();
  }

  if (TString(histClass) == "photon_hbt") {
    const int nm_hbt = 6;
    double m_hbt[nm_hbt] = {0.0, 0.14, 0.5, 1.1, 2.0, 2.5};
    THnSparseF* hs_q_same = nullptr;
    THnSparseF* hs_q_mix = nullptr;

    if (TString(subGroup) == "1d") {
      const int ndim_1d = 4; // m1, m2, kt, qinv
      const int nbins_1d[ndim_1d] = {nm_hbt - 1, nm_hbt - 1, 10, 40};
      const double xmin_1d[ndim_1d] = {0.0, 0.0, 0.0, 0.0};
      const double xmax_1d[ndim_1d] = {2.5, 2.5, 1.0, 0.4};

      hs_q_same = new THnSparseF("hs_q_same", "hs_q_same;m_{1} (GeV/c^{2});m_{2} (GeV/c^{2});k_{T} (GeV/c);q_{inv} (GeV/c);q_{long}^{CMS} (GeV/c);q_{out}^{CMS} (GeV/c);q_{side}^{CMS} (GeV/c);q_{long}^{LCMS} (GeV/c);", ndim_1d, nbins_1d, xmin_1d, xmax_1d);
      hs_q_same->Sumw2();
      hs_q_same->SetBinEdges(0, m_hbt);
      hs_q_same->SetBinEdges(1, m_hbt);
      hs_q_mix = reinterpret_cast<THnSparseF*>(hs_q_same->Clone("hs_q_mix"));
      list->Add(hs_q_same);
      list->Add(hs_q_mix);
    } else if (TString(subGroup) == "3d") {
      const int ndim_3d = 8; // m1, m2, kt, qinv, qlong_cms, qout_cms, qside_cms, qlong_lcms
      const int nbins_3d[ndim_3d] = {nm_hbt - 1, nm_hbt - 1, 10, 40, 80, 80, 80, 80};
      const double xmin_3d[ndim_3d] = {0.0, 0.0, 0.0, 0.0, -0.4, -0.4, -0.4, -0.4};
      const double xmax_3d[ndim_3d] = {2.5, 2.5, 1.0, 0.4, +0.4, +0.4, +0.4, +0.4};

      hs_q_same = new THnSparseF("hs_q_same", "hs_q_same;m_{1} (GeV/c^{2});m_{2} (GeV/c^{2});k_{T} (GeV/c);q_{inv} (GeV/c);q_{long}^{CMS} (GeV/c);q_{out}^{CMS} (GeV/c);q_{side}^{CMS} (GeV/c);q_{long}^{LCMS} (GeV/c);", ndim_3d, nbins_3d, xmin_3d, xmax_3d);
      hs_q_same->Sumw2();
      hs_q_same->SetBinEdges(0, m_hbt);
      hs_q_same->SetBinEdges(1, m_hbt);
      hs_q_mix = reinterpret_cast<THnSparseF*>(hs_q_same->Clone("hs_q_mix"));
      list->Add(hs_q_same);
      list->Add(hs_q_mix);
    }
  }
}
THashList* o2::aod::emphotonhistograms::AddHistClass(THashList* list, const char* histClass)
{
  if (list->FindObject(histClass)) {
    LOGF(info, "HistogramsLibrary::AddHistClass(): Cannot add histogram class %s because it already exists.", histClass);
    return static_cast<THashList*>(list->FindObject(histClass));
  }

  auto* sublist = new THashList();
  sublist->SetOwner(true);
  sublist->SetName(histClass);
  list->Add(sublist);
  return sublist;
}
