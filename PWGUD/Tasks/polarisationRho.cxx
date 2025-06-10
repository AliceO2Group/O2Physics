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
#include <iostream>
#include "PWGUD/DataModel/UDTables.h"
#include <TString.h>
#include "TLorentzVector.h"
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
#define mpion 0.1396   // mass of pion
#define mkaon 0.493696 // mass of kaon
// #define mmuon 0.1057 // mass of muon

/// \brief Rho polarisation task, exclusive four pions, exclusive two particles without PID
/// \author Simone Ragoni, Creighton
/// \author Anisa Khatun, Kansas University
/// \date 14/2/2024

struct PolarisationRho {

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

    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hMassLike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 20.}});
    registry.add("hMassUnlike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 20.}});
    registry.add("hMassUnlikePt", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 20.}});
    registry.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});
    registry.add("hdEdxPion", "p_{#pi} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});
    registry.add("hdEdxKaon", "p_{#pi} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{10, 0., 10.}});

    TString SelectionCuts[7] = {"NoSelection", "2pioncandidates", "tpcNClsCrossedRows", "opposite_charge", "|nsigmapi|<3", "nsigmael", "track_momenta>0.3GeV/c"};
    // now we can set BinLabel in histogram Registry

    for (int i = 0; i < 7; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    registry.add("hNsigEvsPi1", "NSigmaPi(t1) vs NSigmapi (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigEvsKa1", "NSigmaKa(t1) vs NSigmaKa (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigEvsPi2", "NSigmaEl(t1) vs NSigmaEl (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigPivsPt1", "Pt vs NSigmaPi (t1);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("hNsigPivsPt2", "Pt vs NSigmaPi (t2);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});

    registry.add("hNsigElvsPt1", "Pt vs NSigmaEl (t1);#it{p_{t}}, GeV/c;n#sigma_{#e}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("hNsigElvsPt2", "Pt vs NSigmaEl (t2);#it{p_{t}}, GeV/c;n#sigma_{#e}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});

    registry.add("hNsigMuvsPt1", "Pt vs NSigmaMu (t1);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("hNsigMuvsPt2", "Pt vs NSigmaMu (t2);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});

    registry.add("hCharge", "Charge;#it{charge};", kTH1F, {{500, -10., 10.}});

    registry.add("hPtsingle_track1", "Pt t1;#it{p_{t}}, GeV/c;", kTH1F, {{600, 0., 3.}});
    registry.add("hPtsingle_track2", "Pt t2;#it{p_{t}}, GeV/c;", kTH1F, {{600, 0., 3.}});

    registry.add("hP1", "P vs TPC signal;#it{P_{track}}, GeV/c; signal_{TPC} t1", kTH2F, {{100, 0., 2.}, {300, 0, 150}});
    registry.add("hTPCsig", "TPC signal;signal_{TPC} t2; signal_{TPC} t2", kTH2F, {{300, 0., 150.}, {300, 0, 150}});
    registry.add("hMPt1", "Inv.M vs track Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 2.}, {100, 0., 2.}});

    registry.add("hPhi1", "Phi (t1);#it{#phi};", kTH1F, {{120, 0., 6.28}});
    registry.add("hPhi2", "Phi (t1);#it{#phi};", kTH1F, {{120, 0., 6.28}});

    registry.add("hEta_t1", "Eta of t1;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("hEta_t2", "Eta of t2;#it{#eta};", kTH1F, {{100, -5., 5.}});

    registry.add("hCostheta_Phi", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {120, -1, 1}});

    registry.add("hMass", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassFourPions", "Raw Inv.M;#it{m_{4#pi}}, GeV/c^{2};", kTH1D, {{600, 0., 6.}});
    registry.add("hMassSixPions", "Raw Inv.M;#it{m_{6#pi}}, GeV/c^{2};", kTH1D, {{600, 0., 6.}});
    registry.add("hMassEightPions", "Raw Inv.M;#it{m_{8#pi}}, GeV/c^{2};", kTH1D, {{600, 0., 6.}});
    registry.add("hMassFourPionsRightSign", "Raw Inv.M;#it{m_{4#pi}}, GeV/c^{2};", kTH1D, {{600, 0., 6.}});
    registry.add("hMassSixPionsRightSign", "Raw Inv.M;#it{m_{6#pi}}, GeV/c^{2};", kTH1D, {{600, 0., 6.}});
    registry.add("hPhiEtaSixPionsRightSign", "Phi vs Eta;#it{#phi};#it{#eta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {120, -2, 2}});
    registry.add("hMassPhi", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("hMassPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("hMassPhiWrongMomentumWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("hRapidityPhiWithoutPID", "Rapidity;#it{y_{KK}};", kTH1D, {{100, -2., 2.}});
    registry.add("hCosThetaPhiWithoutPID", "CosTheta;cos(#theta);", kTH1D, {{100, -2., 2.}});
    registry.add("hPhiPhiWithoutPID", "Phi;#varphi;", kTH1D, {{100, -2. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("hPtKaonVsKaon", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("hMassPhiWrongMomentaWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("hMassPhiWithoutPIDPionHypothesis", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1D, {{400, 0., 4.}});
    registry.add("hMass1", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMass2", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});

    registry.add("hMassCosTheta_0", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_1", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_2", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_3", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_4", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_5", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_6", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_7", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_8", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_9", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_10", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_11", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_12", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_13", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_14", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_15", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassCosTheta_16", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});

    registry.add("hMassPhi_0", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_1", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_2", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_3", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_4", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_5", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_6", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_7", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_8", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_9", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_10", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_11", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hMassPhi_12", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});

    registry.add("hMasseta2", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
    registry.add("hPt", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
    registry.add("hPtPhi", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
    registry.add("hPtPhiWithoutPID", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
    registry.add("hPtPhiWithoutPIDPionHypothesis", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
    registry.add("hPt1", "Pt1;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
    registry.add("hPt2", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
    registry.add("hEta", "Eta;#it{#eta};", kTH1F, {{500, -10., 10.}});
    registry.add("hRap", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
    registry.add("hPhi", "Phi;#it{#Phi};", kTH1F, {{120, 0, 2. * TMath::Pi()}});
    registry.add("hCostheta", "Costheta;#it{Cos#Theta};", kTH1F, {{100, -1., 1.}});
    registry.add("hMPt", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 2.}, {100, 0., 2.}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;

  //__________________________________________________________________________
  // Main process
  void process(UDCollisions::iterator const& /*collision*/, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);

    TLorentzVector p, phi, phiWithoutPID, phiWrongMomentaWithoutPID, phiWithoutPIDPionHypothesis; // lorentz vectors of tracks and the mother

    //_____________________________________
    // Create kaons and apply TPC Ka PID
    std::vector<TLorentzVector> onlyKaonTracks;
    std::vector<float> onlyKaonSigma;
    std::vector<decltype(tracks.begin())> rawKaonTracks;
    for (auto t : tracks) {
      if (!t.isPVContributor()) {
        continue;
      }
      double momentum = TMath::Sqrt(t.px() * t.px() + t.py() * t.py() + t.pz() * t.pz());
      double dEdx = t.tpcSignal();
      TLorentzVector a;
      a.SetXYZM(t.px(), t.py(), t.pz(), mkaon);
      float nSigmaKa = t.tpcNSigmaKa();
      if (fabs(nSigmaKa) < 2.) {
        onlyKaonTracks.push_back(a);
        onlyKaonSigma.push_back(nSigmaKa);
        rawKaonTracks.push_back(t);
        registry.fill(HIST("hdEdxKaon"), momentum, dEdx);
      }
    }
    //_____________________________________
    // Creating phis
    if (onlyKaonTracks.size() == 2) {
      for (auto kaon : onlyKaonTracks) {
        phi += kaon;
      }
      registry.fill(HIST("hNsigEvsKa1"), rawKaonTracks[0].tpcNSigmaKa(), rawKaonTracks[1].tpcNSigmaKa());
      if (phi.M() < 1.1) {
        registry.fill(HIST("hPtPhi"), phi.Pt());
      }
      registry.fill(HIST("hMassPhi"), phi.M());
    }

    //_____________________________________
    // Create kaons WITHOUT PID
    std::vector<TLorentzVector> allTracksAreKaons;
    std::vector<TLorentzVector> allTracksAreKaonsWrongMomentum;
    std::vector<TLorentzVector> allTracksArePions;
    int counter = 0;
    for (auto t : tracks) {
      if (!t.isPVContributor()) {
        continue;
      }
      int NFindable = t.tpcNClsFindable();
      int NMinusFound = t.tpcNClsFindableMinusFound();
      int NCluster = NFindable - NMinusFound;
      if (NCluster > 30) {
        continue;
      }
      TLorentzVector a;
      a.SetXYZM(t.px(), t.py(), t.pz(), mkaon);
      TLorentzVector b;
      b.SetXYZM(t.px(), t.py(), t.pz(), mpion);
      TLorentzVector a2;
      a2.SetXYZM(-1. * t.px(), -1. * t.py(), -1. * t.pz(), mkaon);
      allTracksAreKaons.push_back(a);
      allTracksArePions.push_back(b);
      if (counter < 1) {
        allTracksAreKaonsWrongMomentum.push_back(a);
      } else {
        allTracksAreKaonsWrongMomentum.push_back(a2);
      }
      counter += 1;
    }
    //_____________________________________
    // Creating phis
    if (allTracksAreKaons.size() == 2) {
      for (auto kaon : allTracksAreKaons) {
        phiWithoutPID += kaon;
      }
      for (auto kaon : allTracksAreKaonsWrongMomentum) {
        phiWrongMomentaWithoutPID += kaon;
      }
      for (auto pion : allTracksArePions) {
        phiWithoutPIDPionHypothesis += pion;
      }
      if (phiWithoutPID.M() < 1.05) {
        auto costhetaPhi = CosThetaHelicityFrame(allTracksAreKaons[0], allTracksAreKaons[1], phiWithoutPID);
        auto phiPhi = 1. * TMath::Pi() + PhiHelicityFrame(allTracksAreKaons[0], allTracksAreKaons[1], phiWithoutPID);
        registry.fill(HIST("hPtPhiWithoutPID"), phiWithoutPID.Pt());
        registry.fill(HIST("hCosThetaPhiWithoutPID"), costhetaPhi);
        registry.fill(HIST("hPhiPhiWithoutPID"), phiPhi);
        registry.fill(HIST("hPtPhiWithoutPIDPionHypothesis"), phiWithoutPIDPionHypothesis.Pt());
        registry.fill(HIST("hMassPhiWithoutPIDPionHypothesis"), phiWithoutPIDPionHypothesis.M());
        registry.fill(HIST("hRapidityPhiWithoutPID"), phiWithoutPID.Rapidity());
        registry.fill(HIST("hPtKaonVsKaon"), allTracksAreKaons[0].Pt(), allTracksAreKaons[1].Pt());
      }
      registry.fill(HIST("hMassPhiWithoutPID"), phiWithoutPID.M());
      registry.fill(HIST("hMassPhiWrongMomentumWithoutPID"), phiWrongMomentaWithoutPID.M());
    }

    //_____________________________________
    // Create pions and apply TPC Pi PID
    std::vector<TLorentzVector> allTracks;
    std::vector<TLorentzVector> onlyPionTracks;
    std::vector<float> onlyPionSigma;
    std::vector<decltype(tracks.begin())> rawPionTracks;
    registry.fill(HIST("hTracks"), tracks.size());
    float sign = 1.;
    for (auto t : tracks) {
      if (!t.isPVContributor()) {
        continue;
      }
      double momentum = TMath::Sqrt(t.px() * t.px() + t.py() * t.py() + t.pz() * t.pz());
      double dEdx = t.tpcSignal();
      if (momentum < 0.1) {
        continue;
      }
      registry.fill(HIST("hdEdx"), momentum, dEdx);
      TLorentzVector a;
      a.SetXYZM(t.px(), t.py(), t.pz(), mpion);
      allTracks.push_back(a);
      float nSigmaPi = t.tpcNSigmaPi();
      if (fabs(nSigmaPi) < 5.) {
        onlyPionTracks.push_back(a);
        onlyPionSigma.push_back(nSigmaPi);
        rawPionTracks.push_back(t);
        sign *= t.sign();
        registry.fill(HIST("hdEdxPion"), momentum, dEdx);
      }
    }
    registry.fill(HIST("hTracksPions"), onlyPionTracks.size());
    //_____________________________________
    // Creating rhos
    for (auto pion : onlyPionTracks) {
      p += pion;
    }
    //_____________________________________
    // Four pions analysis
    if (onlyPionTracks.size() == 4) {
      registry.fill(HIST("hMassFourPions"), p.M());
      int signSum = 0;
      for (auto rawPion : rawPionTracks) {
        if (rawPion.sign() > 0) {
          signSum += 1;
        } else if (rawPion.sign() < 0) {
          signSum -= 1;
        }
      }
      if (signSum == 0) {
        registry.fill(HIST("hMassFourPionsRightSign"), p.M());
      }
    }
    //_____________________________________
    // Six pions analysis
    if (onlyPionTracks.size() == 6) {
      registry.fill(HIST("hMassSixPions"), p.M());
      int signSum = 0;
      for (auto rawPion : rawPionTracks) {
        if (rawPion.sign() > 0) {
          signSum += 1;
        } else if (rawPion.sign() < 0) {
          signSum -= 1;
        }
      }
      if (signSum == 0) {
        registry.fill(HIST("hMassSixPionsRightSign"), p.M());
        for (auto pion : onlyPionTracks) {
          registry.fill(HIST("hPhiEtaSixPionsRightSign"), pion.Phi(), pion.Eta());
        }
      }
    }
    //_____________________________________
    // Eight pions analysis
    if (onlyPionTracks.size() == 8) {
      registry.fill(HIST("hMassEightPions"), p.M());
    }
    //_____________________________________
    // CUTS
    // registry.fill(HIST("hMassUnlike"), p.M());
    if (onlyPionTracks.size() != 2) {
      return;
    }
    if ((onlyPionSigma[0] * onlyPionSigma[0] + onlyPionSigma[1] * onlyPionSigma[1]) > 25) {
      return;
    }
    // Removing electrons
    if (fabs(rawPionTracks[0].tpcNSigmaEl()) > 7.0) {
      return;
    }
    if (onlyPionTracks[0].P() < 0.3 || onlyPionTracks[1].P() < 0.3) {
      return;
    }
    if (p.Rapidity() < -0.8 || p.Rapidity() > 0.8) {
      return;
    }
    if (rawPionTracks[0].sign() == rawPionTracks[1].sign()) {
      registry.fill(HIST("hMassLike"), p.M());
      return;
    }
    registry.fill(HIST("hMassUnlike"), p.M());
    TLorentzVector tplus, tminus;
    if (rawPionTracks[0].sign() > 0) {
      tplus = onlyPionTracks[0];
    } else {
      tminus = onlyPionTracks[0];
    }
    if (rawPionTracks[1].sign() < 0) {
      tminus = onlyPionTracks[1];
    } else {
      tplus = onlyPionTracks[1];
    }

    //______________________________
    // Quality Control
    registry.fill(HIST("hNsigEvsPi1"), rawPionTracks[0].tpcNSigmaPi(), rawPionTracks[1].tpcNSigmaPi());
    registry.fill(HIST("hNsigEvsPi2"), rawPionTracks[0].tpcNSigmaEl(), rawPionTracks[1].tpcNSigmaEl());
    registry.fill(HIST("hNsigPivsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaPi());
    registry.fill(HIST("hNsigPivsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaPi());
    registry.fill(HIST("hNsigElvsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaEl());
    registry.fill(HIST("hNsigElvsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaEl());
    registry.fill(HIST("hNsigMuvsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaMu());
    registry.fill(HIST("hNsigMuvsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaMu());

    registry.fill(HIST("hCharge"), sign);
    registry.fill(HIST("hPtsingle_track1"), onlyPionTracks[0].Pt());
    registry.fill(HIST("hPtsingle_track2"), onlyPionTracks[1].Pt());
    registry.fill(HIST("hP1"), onlyPionTracks[0].P(), rawPionTracks[0].tpcSignal());
    registry.fill(HIST("hTPCsig"), rawPionTracks[0].tpcSignal(), rawPionTracks[1].tpcSignal());
    registry.fill(HIST("hMPt1"), p.M(), onlyPionTracks[0].Pt());

    registry.fill(HIST("hPt"), p.Pt());
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hEta"), p.Eta());
    registry.fill(HIST("hRap"), p.Rapidity());
    // registry.fill(HIST("hPhi"),  p.Phi());
    registry.fill(HIST("hMPt"), p.M(), p.Pt());
    if (p.Pt() < 0.2) {
      registry.fill(HIST("hMassUnlikePt"), p.M());
    }

    //____________________________
    // Polarisation
    if (p.Pt() < 0.2) {

      auto costheta = CosThetaHelicityFrame(tplus, tminus, p);
      registry.fill(HIST("hCostheta"), costheta);
      auto phihel = 1. * TMath::Pi() + PhiHelicityFrame(tplus, tminus, p);
      registry.fill(HIST("hPhi"), phihel);
      registry.fill(HIST("hCostheta_Phi"), phihel, costheta);

      if (costheta > -1. && costheta < -0.7) {
        registry.fill(HIST("hMassCosTheta_0"), p.M());
      }
      if (costheta > -0.7 && costheta < -0.55) {
        registry.fill(HIST("hMassCosTheta_1"), p.M());
      }
      if (costheta > -0.55 && costheta < -0.45) {
        registry.fill(HIST("hMassCosTheta_2"), p.M());
      }
      if (costheta > -0.45 && costheta < -0.35) {
        registry.fill(HIST("hMassCosTheta_3"), p.M());
      }
      if (costheta > -0.35 && costheta < -0.25) {
        registry.fill(HIST("hMassCosTheta_4"), p.M());
      }
      if (costheta > -0.25 && costheta < -0.15) {
        registry.fill(HIST("hMassCosTheta_5"), p.M());
      }
      if (costheta > -0.15 && costheta < -0.05) {
        registry.fill(HIST("hMassCosTheta_6"), p.M());
      }
      if (costheta > -0.05 && costheta < 0.) {
        registry.fill(HIST("hMassCosTheta_7"), p.M());
      }
      if (costheta > 0. && costheta < 0.05) {
        registry.fill(HIST("hMassCosTheta_8"), p.M());
      }
      if (costheta > 0.05 && costheta < 0.15) {
        registry.fill(HIST("hMassCosTheta_9"), p.M());
      }
      if (costheta > 0.15 && costheta < 0.25) {
        registry.fill(HIST("hMassCosTheta_10"), p.M());
      }
      if (costheta > 0.25 && costheta < 0.35) {
        registry.fill(HIST("hMassCosTheta_11"), p.M());
      }
      if (costheta > 0.35 && costheta < 0.45) {
        registry.fill(HIST("hMassCosTheta_12"), p.M());
      }
      if (costheta > 0.45 && costheta < 0.55) {
        registry.fill(HIST("hMassCosTheta_13"), p.M());
      }
      if (costheta > 0.55 && costheta < 0.7) {
        registry.fill(HIST("hMassCosTheta_14"), p.M());
      }
      if (costheta > 0.7 && costheta < 1.) {
        registry.fill(HIST("hMassCosTheta_15"), p.M());
      }

      if (phihel > 2. * TMath::Pi() / 12. * 0. && phihel < 2. * TMath::Pi() / 12. * 1.) {
        registry.fill(HIST("hMassPhi_0"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 1. && phihel < 2. * TMath::Pi() / 12. * 2.) {
        registry.fill(HIST("hMassPhi_1"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 2. && phihel < 2. * TMath::Pi() / 12. * 3.) {
        registry.fill(HIST("hMassPhi_2"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 3. && phihel < 2. * TMath::Pi() / 12. * 4.) {
        registry.fill(HIST("hMassPhi_3"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 4. && phihel < 2. * TMath::Pi() / 12. * 5.) {
        registry.fill(HIST("hMassPhi_4"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 5. && phihel < 2. * TMath::Pi() / 12. * 6.) {
        registry.fill(HIST("hMassPhi_5"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 6. && phihel < 2. * TMath::Pi() / 12. * 7.) {
        registry.fill(HIST("hMassPhi_6"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 7. && phihel < 2. * TMath::Pi() / 12. * 8.) {
        registry.fill(HIST("hMassPhi_7"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 8. && phihel < 2. * TMath::Pi() / 12. * 9.) {
        registry.fill(HIST("hMassPhi_8"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 9. && phihel < 2. * TMath::Pi() / 12. * 10.) {
        registry.fill(HIST("hMassPhi_9"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 10. && phihel < 2. * TMath::Pi() / 12. * 11.) {
        registry.fill(HIST("hMassPhi_10"), p.M());
      }
      if (phihel > 2. * TMath::Pi() / 12. * 11. && phihel < 2. * TMath::Pi() / 12. * 12.) {
        registry.fill(HIST("hMassPhi_11"), p.M());
      }
    }

  } // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PolarisationRho>(cfgc)};
}
