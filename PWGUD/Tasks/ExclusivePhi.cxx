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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
#include <TString.h>
#include "TLorentzVector.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGUD/Core/SGSelector.h"
using std::array;
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
#define mpion 0.1396   // mass of pion
#define mkaon 0.493696 // mass of kaon
// #define mmuon 0.1057 // mass of muon

/// \brief Exclusive phi without PID
/// \author Simone Ragoni, Creighton
/// \author Anisa Khatun, Kansas University
/// \date 14/2/2024

struct ExclusivePhi {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //_____________________________________________________________________________
  Double_t CosThetaHelicityFrame(TLorentzVector pionPositive,
                                 TLorentzVector pionNegative,
                                 TLorentzVector possibleRhoZero)
  {

    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    TVector3 beta = (-1. / possibleRhoZero.E()) * possibleRhoZero.Vect();
    TLorentzVector pPi1Dipion = pionPositive;
    TLorentzVector pPi2Dipion = pionNegative;
    TLorentzVector pProjDipion = pProjCM;
    TLorentzVector pTargDipion = pTargCM;

    pPi1Dipion.Boost(beta);
    pPi2Dipion.Boost(beta);
    pProjDipion.Boost(beta);
    pTargDipion.Boost(beta);

    TVector3 zaxis = (possibleRhoZero.Vect()).Unit();

    Double_t CosThetaHE = zaxis.Dot((pPi1Dipion.Vect()).Unit());
    return CosThetaHE;
  }
  //------------------------------------------------------------------------------------------------------
  Double_t PhiHelicityFrame(TLorentzVector muonPositive, TLorentzVector muonNegative, TLorentzVector possibleJPsi)
  {

    // Half of the energy per pair of the colliding nucleons.
    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    // Translate the dimuon parameters in the dimuon rest frame
    TVector3 beta = (-1. / possibleJPsi.E()) * possibleJPsi.Vect();
    TLorentzVector pMu1Dimu = muonPositive;
    TLorentzVector pMu2Dimu = muonNegative;
    TLorentzVector pProjDimu = pProjCM;
    TLorentzVector pTargDimu = pTargCM;
    pMu1Dimu.Boost(beta);
    pMu2Dimu.Boost(beta);
    pProjDimu.Boost(beta);
    pTargDimu.Boost(beta);

    // Axes
    TVector3 zaxis = (possibleJPsi.Vect()).Unit();
    TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
    TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
    //
    // --- Calculation of the azimuthal angle (Helicity)
    //
    Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis), (pMu1Dimu.Vect()).Dot(xaxis));
    return phi;
  }
  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {
    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});

    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hITSCluster", "N_{cluster}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hChi2ITSTrkSegment", "N_{cluster}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTPCCluster", "N_{cluster}", kTH1F, {{200, -0.5, 199.5}});
    registry.add("hTracksKaons", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon2", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon3", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hNsigEvsKa1", "NSigmaKa(t1) vs NSigmaKa (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigEvsKa2", "NSigmaKa(t1) vs NSigmaKa (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hMomentum", "p_{#ka};#it{p_{trk}}, GeV/c;", kTH1F, {{100, 0., 3.}});
    registry.add("hEta1", "#eta_{#ka};#it{#eta_{trk}}, GeV/c;", kTH1F, {{100, -2., 2.}});
    registry.add("hEta2", "#eta_{#ka};#it{#eta_{trk}}, GeV/c;", kTH1F, {{100, -2., 2.}});
    registry.add("hPtPhi", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("hMassPhi", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("hMassPtPhi", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{10, 0., 10.}});

    TString SelectionCuts[9] = {"NoSelection", "GAPcondition", "PVtracks", "TPCcrossrow", "|nsigmaka|<2", "doublegap", "two tracks", "Phi-peak", "pt<0.2 GeV/c"};
    // now we can set BinLabel in histogram Registry

    for (int i = 0; i < 9; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    registry.add("hMassPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("hMassPtPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("hMassPhiWrongMomentumWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});

    // Phi peak region
    registry.add("PHI/hPtPhiWithoutPID", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hMassVsPt", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PHI/hRapidityPhiWithoutPID", "Rapidity;#it{y_{KK}};", kTH1F, {{100, -2., 2.}});
    registry.add("PHI/hCosThetaPhiWithoutPID", "CosTheta;cos(#theta);", kTH1F, {{100, -2., 2.}});
    registry.add("PHI/hPhiPhiWithoutPID", "Phi;#varphi;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PHI/hPtKaonVsKaon", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("PHI/hCostheta_Phi", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {100, -2, 2}});
    registry.add("PHI/hMassPhiWithoutPIDPionHypothesis", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("PHI/hPtPhiWithoutPIDPionHypothesis", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});

    registry.add("PHI/hMassLike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hMassUnlike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hUnlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hCoherentPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hInCoherentPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hCoherentMassLike", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hInCoherentMassLike", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});

    // High Mass region
    registry.add("PHIHIGH/hPtPhiWithoutPID1", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHIHIGH/hMassVsPt1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PHIHIGH/hRapidityPhiWithoutPID1", "Rapidity;#it{y_{KK}};", kTH1F, {{100, -2., 2.}});
    registry.add("PHIHIGH/hCosThetaPhiWithoutPID1", "CosTheta;cos(#theta);", kTH1F, {{100, -2., 2.}});
    registry.add("PHIHIGH/hPhiPhiWithoutPID1", "Phi;#varphi;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PHIHIGH/hPtKaonVsKaon1", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("PHIHIGH/hCostheta_Phi1", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {100, -2, 2}});
    registry.add("PHIHIGH/hMassPhiWithoutPIDPionHypothesis1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("PHIHIGH/hPtPhiWithoutPIDPionHypothesis1", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});

    registry.add("PHIHIGH/hMassLike1", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hMassUnlike1", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hlikePt1", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hUnlikePt1", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hCoherentPhiWithoutPID1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hInCoherentPhiWithoutPID1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hCoherentMassLike1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHIHIGH/hInCoherentMassLike1", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});

    // Low Mass region

    registry.add("PHILOW/hPtPhiWithoutPID2", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHILOW/hMassVsPt2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PHILOW/hRapidityPhiWithoutPID2", "Rapidity;#it{y_{KK}};", kTH1F, {{100, -2., 2.}});
    registry.add("PHILOW/hCosThetaPhiWithoutPID2", "CosTheta;cos(#theta);", kTH1F, {{100, -2., 2.}});
    registry.add("PHILOW/hPhiPhiWithoutPID2", "Phi;#varphi;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PHILOW/hPtKaonVsKaon2", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("PHILOW/hCostheta_Phi2", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {100, -2, 2}});
    registry.add("PHILOW/hMassPhiWithoutPIDPionHypothesis2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("PHILOW/hPtPhiWithoutPIDPionHypothesis2", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});

    registry.add("PHILOW/hMassLike2", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hMassUnlike2", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hlikePt2", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hUnlikePt2", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hCoherentPhiWithoutPID2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hInCoherentPhiWithoutPID2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hCoherentMassLike2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHILOW/hInCoherentMassLike2", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  // Main process
  void process(UDCollisions::iterator const& collision, udtracksfull const& tracks)
  //  void process(UDCollisionsFull::iterator const& collision, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);

    /*int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;*/
    // if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)){
    // if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    registry.fill(HIST("hSelectionCounter"), 1);

    /* int truegapSide = sgSelector.trueGap(collision, FV0_cut, ZDC_cut);
     registry.fill(HIST("GapSide"), gapSide);
     registry.fill(HIST("TrueGapSide"), truegapSide);
     gapSide = truegapSide;*/

    TLorentzVector phi, phiWithoutPID, phiWithKaonPID, phiWrongMomentaWithoutPID, phiWithoutPIDPionHypothesis; // lorentz vectors of tracks and the mother

    // ===================================
    // Task for phi WITH PID FROM TPC
    // ===================================
    // Here the tracks are allowed to have
    // at most 30 TPC clusters
    // TPC PID is then possible, although
    // not truly reliable
    // NB: from MC simulations coherent
    // phi SHOULD NOT leave clusters
    // in the TPC at all
    std::vector<TLorentzVector> onlyKaonTracks;
    std::vector<float> onlyKaonSigma;
    std::vector<decltype(tracks.begin())> rawKaonTracks;

    for (auto trk : tracks) {
      if (!trk.isPVContributor()) {
        continue;
      }
      registry.fill(HIST("hSelectionCounter"), 2);

      int NFindable = trk.tpcNClsFindable();
      int NMinusFound = trk.tpcNClsFindableMinusFound();
      int NCluster = NFindable - NMinusFound;
      registry.fill(HIST("hTPCCluster"), NCluster);
      registry.fill(HIST("hITSCluster"), trk.itsNCls());
      registry.fill(HIST("hChi2ITSTrkSegment"), trk.itsChi2NCl());

      registry.fill(HIST("hSelectionCounter"), 3);

      double momentum = TMath::Sqrt(trk.px() * trk.px() + trk.py() * trk.py() + trk.pz() * trk.pz());
      double dEdx = trk.tpcSignal();
      registry.fill(HIST("hdEdx"), momentum, dEdx);

      /* if(trk.pt() > 0.180){
               continue;
           }*/

      TLorentzVector kaon;
      kaon.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassKaonCharged);
      auto nSigmaKa = trk.tpcNSigmaKa();

      if (fabs(nSigmaKa) < 2.) {
        onlyKaonTracks.push_back(kaon);
        onlyKaonSigma.push_back(nSigmaKa);
        rawKaonTracks.push_back(trk);
        registry.fill(HIST("hdEdxKaon"), momentum, dEdx);
        registry.fill(HIST("hSelectionCounter"), 4);
      }

    } // trk loop

    // registry.fill(HIST("hTracksKaons"), rawKaonTracks.size());

    // Creating phis using Kaon PID
    //  if (gapSide == 2) {
    registry.fill(HIST("hSelectionCounter"), 5);
    if (onlyKaonTracks.size() == 2) {
      registry.fill(HIST("hSelectionCounter"), 6);

      for (auto kaon : onlyKaonTracks) {
        phi += kaon;
      }

      registry.fill(HIST("hNsigEvsKa1"), rawKaonTracks[0].tpcNSigmaKa(), rawKaonTracks[1].tpcNSigmaKa());
      registry.fill(HIST("hMassPtPhi"), phi.M(), phi.Pt());
      // if(rawKaonTracks[0].sign() != rawKaonTracks[1].sign()){
      if ((phi.M() > 0.98) && (phi.M() < 1.06)) {
        registry.fill(HIST("hSelectionCounter"), 7);
        registry.fill(HIST("hPtPhi"), phi.Pt());
      }
      if (phi.Pt() < 0.2) {
        registry.fill(HIST("hMassPhi"), phi.M());
        registry.fill(HIST("hSelectionCounter"), 8);
      }

      // }
    }
    //  } // double gap

    // ===================================
    // Task for phi WITHOUT PID
    // ===================================

    // ====================================
    // Selections for events to be stored
    // ------------------------------------
    // - PV: only PV contributors
    // - only two PV contributors
    // - both PV have t.hasITS() ON
    // - only ITS tracks
    //_____________________________________
    // Create kaons WITHOUT PID
    std::vector<decltype(tracks.begin())> onlyTwoTracks;
    std::vector<TLorentzVector> allTracksAreKaons;
    std::vector<TLorentzVector> allTracksArePions;
    std::vector<TLorentzVector> allTracksAreKaonsWrongMomentum;

    int counter = 0;
    for (auto t : tracks) {
      if (!t.isPVContributor()) {
        continue;
      }

      registry.fill(HIST("hTracks"), t.size());

      /* if(t.itsNCls() < 4) {
           continue;
       }*/

      double momentum = TMath::Sqrt(t.px() * t.px() + t.py() * t.py() + t.pz() * t.pz());
      double dEdx = t.tpcSignal();

      registry.fill(HIST("hMomentum"), momentum);

      /* if ((!t.hasITS()) && (t.hasTPC()) && (t.hasTOF())) {
           continue;
       } // not working*/

      int NFindable = t.tpcNClsFindable();
      int NMinusFound = t.tpcNClsFindableMinusFound();
      int NCluster = NFindable - NMinusFound;
      // registry.fill(HIST("hTPCCluster"), NCluster);

      if (NCluster > 50) {
        continue;
      }

      registry.fill(HIST("hdEdxKaon2"), momentum, dEdx);
      if (t.pt() > 0.180) {
        continue;
      }

      registry.fill(HIST("hdEdxKaon3"), momentum, dEdx);

      onlyTwoTracks.push_back(t);

      TLorentzVector a;
      a.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassKaonCharged);
      TLorentzVector b;
      b.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassPionCharged);
      TLorentzVector a2;
      a2.SetXYZM(-1. * t.px(), -1. * t.py(), -1. * t.pz(), o2::constants::physics::MassKaonCharged);

      if (-0.6 > a.Eta() || a.Eta() > 0.6)
        continue;

      allTracksAreKaons.push_back(a);
      allTracksArePions.push_back(b);

      if (counter < 1) {
        allTracksAreKaonsWrongMomentum.push_back(a);
      } else {
        allTracksAreKaonsWrongMomentum.push_back(a2);
      }
      counter += 1;
    } // track loop

    //_____________________________________
    // Creating phis and saving all the information
    // in the case that there are ONLY 2 PV
    /* int hasITS[2] = {-1, -1};
     int hasTPC[2] = {-1, -1};
     int hasTOF[2] = {-1, -1};
     bool goodFirstTrack = false;
     bool goodSecondTrack = false;
     bool goodEvent = false;
     if (onlyTwoTracks.size() == 2) {
         hasITS[0] = onlyTwoTracks[0].hasITS();
         hasTPC[0] = onlyTwoTracks[0].hasTPC();
         hasTOF[0] = onlyTwoTracks[0].hasTOF();
         hasITS[1] = onlyTwoTracks[1].hasITS();
         hasTPC[1] = onlyTwoTracks[1].hasTPC();
         hasTOF[1] = onlyTwoTracks[1].hasTOF();
         if ((hasITS[0]) && (hasTPC[0]) && (!hasTOF[0])) {
             goodFirstTrack = true;
         }

         if ((hasITS[1]) && (hasTPC[1]) && (!hasTOF[1])) {
             goodSecondTrack = true;
         }

         if (goodFirstTrack || goodSecondTrack) {
            goodEvent = true;
        }

     }*/

    //  if (gapSide == 2) {
    if (allTracksAreKaons.size() == 2) {

      for (auto kaon : allTracksAreKaons) {
        phiWithoutPID += kaon;
      }
      registry.fill(HIST("hTracksKaons"), allTracksAreKaons.size());
      // kaon mass hypothesis with wrong momentum for one track
      for (auto kaon : allTracksAreKaonsWrongMomentum) {
        phiWrongMomentaWithoutPID += kaon;
      }
      // pion mass hypothesis
      for (auto pion : allTracksArePions) {
        phiWithoutPIDPionHypothesis += pion;
      }

      // TLorentzVector kaons[2], kaonsLikeSign[2], phiLikeSignWithoutPID;
      registry.fill(HIST("hNsigEvsKa2"), onlyTwoTracks[0].tpcNSigmaKa(), onlyTwoTracks[1].tpcNSigmaKa());
      registry.fill(HIST("hEta1"), allTracksAreKaons[0].Eta());
      registry.fill(HIST("hEta2"), allTracksAreKaons[1].Eta());

      /* if (-0.6 > allTracksAreKaons[0].Eta() || allTracksAreKaons[0].Eta() > 0.6)
          return;
        if (-0.6 > allTracksAreKaons[1].Eta() || allTracksAreKaons[1].Eta() > 0.6)
          return;*/

      // kaons[0].SetXYZM(onlyTwoTracks[0].px(), onlyTwoTracks[0].py(), onlyTwoTracks[0].pz(), o2::constants::physics::MassKaonCharged);
      // kaons[1].SetXYZM(onlyTwoTracks[1].px(), onlyTwoTracks[1].py(), onlyTwoTracks[1].pz(), o2::constants::physics::MassKaonCharged);

      // phiWithoutPID += kaons[0];
      // phiWithoutPID += kaons[1];
      /*phiLikeSignWithoutPID += kaonsLikeSign[0];
      phiLikeSignWithoutPID += kaonsLikeSign[1];*/

      auto costhetaPhi = CosThetaHelicityFrame(allTracksAreKaons[0], allTracksAreKaons[1], phiWithoutPID);
      auto phiPhi = 1. * TMath::Pi() + PhiHelicityFrame(allTracksAreKaons[0], allTracksAreKaons[1], phiWithoutPID);

      // All invariant mass region
      registry.fill(HIST("hMassPhiWithoutPID"), phiWithoutPID.M());
      registry.fill(HIST("hMassPtPhiWithoutPID"), phiWithoutPID.M(), phiWithoutPID.Pt());
      registry.fill(HIST("hMassPhiWrongMomentumWithoutPID"), phiWrongMomentaWithoutPID.M());

      // Phi peak region
      // if ((phiWithoutPID.Rapidity()>-0.6)&&(phiWithoutPID.Rapidity()<0.6)){
      // if ((phiWithoutPID.M() > 0.98) && (phiWithoutPID.M() < 1.06)) {
      registry.fill(HIST("PHI/hPtPhiWithoutPID"), phiWithoutPID.Pt());
      registry.fill(HIST("PHI/hMassVsPt"), phiWithoutPID.M(), phiWithoutPID.Pt());
      registry.fill(HIST("PHI/hRapidityPhiWithoutPID"), phiWithoutPID.Rapidity());
      registry.fill(HIST("PHI/hPtKaonVsKaon"), allTracksAreKaons[0].Pt(), allTracksAreKaons[1].Pt());

      registry.fill(HIST("PHI/hPtPhiWithoutPIDPionHypothesis"), phiWithoutPIDPionHypothesis.Pt());
      registry.fill(HIST("PHI/hMassPhiWithoutPIDPionHypothesis"), phiWithoutPIDPionHypothesis.M());

      // unlike-sign
      if (onlyTwoTracks[0].sign() != onlyTwoTracks[1].sign()) {
        registry.fill(HIST("PHI/hCosThetaPhiWithoutPID"), costhetaPhi);
        registry.fill(HIST("PHI/hPhiPhiWithoutPID"), phiPhi);
        registry.fill(HIST("PHI/hCostheta_Phi"), phiPhi, costhetaPhi);
        registry.fill(HIST("PHI/hUnlikePt"), phiWithoutPID.Pt());
        registry.fill(HIST("PHI/hMassUnlike"), phiWithoutPID.M());
        if (phiWithoutPID.Pt() < 0.2) {
          registry.fill(HIST("PHI/hCoherentPhiWithoutPID"), phiWithoutPID.M());
        }
        if (phiWithoutPID.Pt() > 0.2) {
          registry.fill(HIST("PHI/hInCoherentPhiWithoutPID"), phiWithoutPID.M());
        }
      }
      //}//Rapidity
      // Likesign quantities
      if (onlyTwoTracks[0].sign() == onlyTwoTracks[1].sign()) {
        registry.fill(HIST("PHI/hMassLike"), phiWithoutPID.M());
        registry.fill(HIST("PHI/hlikePt"), phiWithoutPID.Pt());

        if (phiWithoutPID.Pt() < 0.2) {
          registry.fill(HIST("PHI/hCoherentMassLike"), phiWithoutPID.M());
        }
        if (phiWithoutPID.Pt() > 0.2) {
          registry.fill(HIST("PHI/hInCoherentMassLike"), phiWithoutPID.M());
        }
      }

      // Side band above phi mass region
      if (phiWithoutPID.M() > 1.06) {

        registry.fill(HIST("PHIHIGH/hPtPhiWithoutPID1"), phiWithoutPID.Pt());
        registry.fill(HIST("PHIHIGH/hMassVsPt1"), phiWithoutPID.M(), phiWithoutPID.Pt());
        registry.fill(HIST("PHIHIGH/hRapidityPhiWithoutPID1"), phiWithoutPID.Rapidity());

        registry.fill(HIST("PHIHIGH/hPtKaonVsKaon1"), allTracksAreKaons[0].Pt(), allTracksAreKaons[1].Pt());
        registry.fill(HIST("PHIHIGH/hPtPhiWithoutPIDPionHypothesis1"), phiWithoutPIDPionHypothesis.Pt());
        registry.fill(HIST("PHIHIGH/hMassPhiWithoutPIDPionHypothesis1"), phiWithoutPIDPionHypothesis.M());

        // unlike-sign
        if (onlyTwoTracks[0].sign() != onlyTwoTracks[1].sign()) {
          registry.fill(HIST("PHIHIGH/hCosThetaPhiWithoutPID1"), costhetaPhi);
          registry.fill(HIST("PHIHIGH/hPhiPhiWithoutPID1"), phiPhi);
          registry.fill(HIST("PHIHIGH/hCostheta_Phi1"), phiPhi, costhetaPhi);

          registry.fill(HIST("PHIHIGH/hUnlikePt1"), phiWithoutPID.Pt());
          registry.fill(HIST("PHIHIGH/hMassUnlike1"), phiWithoutPID.M());
          if (phiWithoutPID.Pt() < 0.2) {
            registry.fill(HIST("PHIHIGH/hCoherentPhiWithoutPID1"), phiWithoutPID.M());
          }
          if (phiWithoutPID.Pt() > 0.2) {
            registry.fill(HIST("PHIHIGH/hInCoherentPhiWithoutPID1"), phiWithoutPID.M());
          }
        }

        // Like sign
        if (onlyTwoTracks[0].sign() == onlyTwoTracks[1].sign()) {
          // Likesign quantities
          registry.fill(HIST("PHIHIGH/hMassLike1"), phiWithoutPID.M());
          registry.fill(HIST("PHIHIGH/hlikePt1"), phiWithoutPID.Pt());

          if (phiWithoutPID.Pt() < 0.2) {
            registry.fill(HIST("PHIHIGH/hCoherentMassLike1"), phiWithoutPID.M());
          }
          if (phiWithoutPID.Pt() > 0.2) {
            registry.fill(HIST("PHIHIGH/hInCoherentMassLike1"), phiWithoutPID.M());
          }
        }
      }

      // Side band below phi mass region
      if (phiWithoutPID.M() < 0.98) {

        registry.fill(HIST("PHILOW/hPtPhiWithoutPID2"), phiWithoutPID.Pt());
        registry.fill(HIST("PHILOW/hMassVsPt2"), phiWithoutPID.M(), phiWithoutPID.Pt());
        registry.fill(HIST("PHILOW/hRapidityPhiWithoutPID2"), phiWithoutPID.Rapidity());

        registry.fill(HIST("PHILOW/hPtKaonVsKaon2"), allTracksAreKaons[0].Pt(), allTracksAreKaons[1].Pt());

        registry.fill(HIST("PHILOW/hPtPhiWithoutPIDPionHypothesis2"), phiWithoutPIDPionHypothesis.Pt());
        registry.fill(HIST("PHILOW/hMassPhiWithoutPIDPionHypothesis2"), phiWithoutPIDPionHypothesis.M());

        // unlike-sign
        registry.fill(HIST("PHILOW/hCosThetaPhiWithoutPID2"), costhetaPhi);
        registry.fill(HIST("PHILOW/hPhiPhiWithoutPID2"), phiPhi);
        registry.fill(HIST("PHILOW/hCostheta_Phi2"), phiPhi, costhetaPhi);
        registry.fill(HIST("PHILOW/hUnlikePt2"), phiWithoutPID.Pt());
        registry.fill(HIST("PHILOW/hMassUnlike2"), phiWithoutPID.M());
        if (phiWithoutPID.Pt() < 0.2) {
          registry.fill(HIST("PHILOW/hCoherentPhiWithoutPID2"), phiWithoutPID.M());
        }
        if (phiWithoutPID.Pt() > 0.2) {
          registry.fill(HIST("PHILOW/hInCoherentPhiWithoutPID2"), phiWithoutPID.M());
        }
        // like-sign
        if (onlyTwoTracks[0].sign() == onlyTwoTracks[1].sign()) {
          // Likesign quantities
          registry.fill(HIST("PHI/hMassLike2"), phiWithoutPID.M());
          registry.fill(HIST("PHI/hlikePt2"), phiWithoutPID.Pt());

          if (phiWithoutPID.Pt() < 0.2) {
            registry.fill(HIST("PHI/hCoherentMassLike2"), phiWithoutPID.M());
          }

          if (phiWithoutPID.Pt() > 0.2) {
            registry.fill(HIST("PHI/hInCoherentMassLike2"), phiWithoutPID.M());
          }
        }
      }
    } // end of two tracks only loop

    // } // double gap
    // }
  } // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusivePhi>(cfgc)};
}
