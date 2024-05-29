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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
#include <TString.h>
#include "TLorentzVector.h"
#include "PWGUD/Core/SGSelector.h"
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// \brief  UPC pion study
/// \author Anisa Khatun
/// \date 11.11.2022

struct UPCPionAnalysis {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 200., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 100., "FT0C threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};
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
  Double_t PhiHelicityFrame(TLorentzVector piPositive, TLorentzVector piNegative, TLorentzVector possibleRho)
  {

    // Half of the energy per pair of the colliding nucleons.
    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    // Translate the dipion parameters in the dipion rest frame
    TVector3 beta = (-1. / possibleRho.E()) * possibleRho.Vect();
    TLorentzVector pPi1Dipi = piPositive;
    TLorentzVector pPi2Dipi = piNegative;
    TLorentzVector pProjDipi = pProjCM;
    TLorentzVector pTargDipi = pTargCM;
    pPi1Dipi.Boost(beta);
    pPi2Dipi.Boost(beta);
    pProjDipi.Boost(beta);
    pTargDipi.Boost(beta);

    // Axes
    TVector3 zaxis = (possibleRho.Vect()).Unit();
    TVector3 yaxis = ((pProjDipi.Vect()).Cross(pTargDipi.Vect())).Unit();
    TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
    //
    // --- Calculation of the azimuthal angle (Helicity)
    //
    Double_t phi = TMath::ATan2((pPi1Dipi.Vect()).Dot(yaxis), (pPi1Dipi.Vect()).Dot(xaxis));
    return phi;
  }
  //____________________________________________________________________________________________
  float momentum(float px, float py, float pz)
  // Just a simple function to return momentum
  {
    return std::sqrt(px * px + py * py + pz * pz);
  }

  //______________________________________________________________________________________________________
  float eta(float px, float py, float pz)
  // Just a simple function to return pseudorapidity
  {
    float arg = -2.; // outside valid range for std::atanh
    float mom = momentum(px, py, pz);
    if (mom != 0)
      arg = pz / mom;
    if (-1. < arg && arg < 1.)
      return std::atanh(arg); // definition of eta
    return -999.;
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {

    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{10, 0., 10.}});

    TString SelectionCuts[7] = {"NoSelection", "2pioncandidates", "tpcNClsCrossedRows", "opposite_charge", "|nsigmapi|<3", "nsigmael", "track_momenta>0.3GeV/c"};
    // now we can set BinLabel in histogram Registry

    for (int i = 0; i < 7; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("h2TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("h4TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("h6TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("h8TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});

    registry.add("hMassLike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 20.}});
    registry.add("hMassUnlike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 20.}});
    registry.add("hMassUnlikeCoherent", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 20.}});
    registry.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});
    registry.add("hdEdxPion", "p_{#pi} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});

    registry.add("hNsigPi1vsPi2", "NSigmaPi(t1) vs NSigmapi (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigEl1vsEl2", "NSigmaEl(t1) vs NSigmaEl (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigPivsPt1", "Pt vs NSigmaPi (t1);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("hNsigPivsPt2", "Pt vs NSigmaPi (t2);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("hNsigElvsPt1", "Pt vs NSigmaEl (t1);#it{p_{t}}, GeV/c;n#sigma_{#e}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("hNsigElvsPt2", "Pt vs NSigmaEl (t2);#it{p_{t}}, GeV/c;n#sigma_{#e}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("hNsigMuvsPt1", "Pt vs NSigmaMu (t1);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("hNsigMuvsPt2", "Pt vs NSigmaMu (t2);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});

    registry.add("hEta_t1", "Eta of t1;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("hEta_t2", "Eta of t2;#it{#eta};", kTH1F, {{100, -5., 5.}});

    registry.add("hCharge", "Charge;#it{charge};", kTH1F, {{500, -10., 10.}});
    registry.add("hPtsingle_track1", "Pt t1;#it{p_{t}}, GeV/c;", kTH1F, {{600, 0., 3.}});
    registry.add("hPtsingle_track2", "Pt t2;#it{p_{t}}, GeV/c;", kTH1F, {{600, 0., 3.}});
    registry.add("hP1", "P vs TPC signal;#it{P_{track}}, GeV/c; signal_{TPC} t1", kTH2F, {{100, 0., 2.}, {300, 0, 150}});
    registry.add("hTPCsig", "TPC signal;signal_{TPC} t2; signal_{TPC} t2", kTH2F, {{300, 0., 150.}, {300, 0, 150}});
    registry.add("hP2", "P vs TPC signal;#it{P_{track}}, GeV/c; signal_{TPC} t1", kTH2F, {{100, 0., 2.}, {300, 0, 150}});
    registry.add("hTPCsig1", "TPC signal;signal_{TPC} t2; signal_{TPC} t2", kTH2F, {{300, 0., 150.}, {300, 0, 150}});

    registry.add("hMassFourPionsCoherent", "Raw Inv.M;#it{M_{#pi#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});
    registry.add("hMassSixPionsCoherent", "Raw Inv.M;#it{M_{6#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});
    registry.add("hMassEightPionsCoherent", "Raw Inv.M;#it{M_{6#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});
    registry.add("hMassFourPions", "Raw Inv.M;#it{M_{4#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});
    registry.add("hMassSixPions", "Raw Inv.M;#it{M_{6#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});
    registry.add("hMassEightPions", "Raw Inv.M;#it{M_{8#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});
    registry.add("hMassFourPionsRightSign", "Raw Inv.M;#it{M_{4#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});
    registry.add("hMassSixPionsRightSign", "Raw Inv.M;#it{M_{6#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});
    registry.add("hMass8PionsRightSign", "Raw Inv.M;#it{M_{8#pi}} (GeV/c^{2});", kTH1D, {{1000, 0., 10.}});

    registry.add("hPhiEtaFourPionsRightSign", "Phi vs Eta;#it{#phi};#it{#eta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {120, -2, 2}});
    registry.add("hPhiEtaSixPionsRightSign", "Phi vs Eta;#it{#phi};#it{#eta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {120, -2, 2}});
    registry.add("hPhiEta8PionsRightSign", "Phi vs Eta;#it{#phi};#it{#eta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {120, -2, 2}});

    registry.add("hPt", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt4Pion", "Pt1;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt6Pion", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt8Pion", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt4PionRightSign", "Pt1;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt6PionRightSign", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt8PionRightSign", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});

    registry.add("hEta", "Eta;#it{#eta};", kTH1F, {{500, -10., 10.}});
    registry.add("hEta4Pion", "Eta of four pion;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("hEta6Pion", "Eta of six pion;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("hEta8Pion", "Eta of six pion;#it{#eta};", kTH1F, {{100, -5., 5.}});

    registry.add("hRap", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
    registry.add("hRap4Pion", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
    registry.add("hRap6Pion", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
    registry.add("hRap8pion", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});

    registry.add("hPhi", "Phi;#it{#Phi};", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("hPhi4Pion", "Phi;#it{#Phi};", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("hPhi6Pion", "Phi;#it{#Phi};", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("hPhi8Pion", "Phi;#it{#Phi};", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});

    registry.add("hCostheta", "Costheta;#it{Cos#Theta};", kTH1F, {{100, -1., 1.}});
    registry.add("hCostheta4Pion", "Costheta;#it{Cos#Theta};", kTH1F, {{100, -1., 1.}});
    registry.add("hCostheta6Pion", "Costheta;#it{Cos#Theta};", kTH1F, {{100, -1., 1.}});
    registry.add("hCostheta8Pion", "Costheta;#it{Cos#Theta};", kTH1F, {{100, -1., 1.}});

    registry.add("hMPt", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 10.}, {100, 0., 10.}});
    registry.add("hMPt1", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 10.}, {100, 0., 10.}});
    registry.add("hMPt2", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 10.}, {100, 0., 10.}});
    registry.add("hMPt3", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 20.}, {100, 0., 10.}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  // Main process
  void process(UDCollisionsFull::iterator const& collision, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);

    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;

    // int truegapSide = sgSelector.trueGap(collision, FV0_cut, ZDC_cut);
    int truegapSide = sgSelector.trueGap(collision, FV0_cut, FT0A_cut, FT0C_cut, ZDC_cut);
    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;

    if (gapSide == gap_Side) {
      //_____________________________________
      // Create pions and apply TPC Pion PID
      std::vector<TLorentzVector> allTracks;
      std::vector<TLorentzVector> onlyPionTracks;
      std::vector<float> onlyPionSigma;
      std::vector<decltype(tracks.begin())> rawPionTracks;
      TLorentzVector p;
      registry.fill(HIST("hTracks"), tracks.size());

      // if(collision.numContrib() > 10) return;

      for (auto t : tracks) {
        if (!t.isPVContributor()) {
          continue;
        }

        double dEdx = t.tpcSignal();

        if (TMath::Abs(eta(t.px(), t.py(), t.pz()) > 0.9)) {
          continue;
        }

        if (t.pt() < 0.15) {
          continue;
        }

        registry.fill(HIST("hdEdx"), t.tpcInnerParam() / t.sign(), dEdx);
        TLorentzVector a;
        a.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassPionCharged);
        allTracks.push_back(a);
        auto nSigmaPi = t.tpcNSigmaPi();

        if (fabs(nSigmaPi) < 3.) {
          onlyPionTracks.push_back(a);
          onlyPionSigma.push_back(nSigmaPi);
          rawPionTracks.push_back(t);
          registry.fill(HIST("hdEdxPion"), t.tpcInnerParam() / t.sign(), dEdx);
        }
      }
      registry.fill(HIST("hTracksPions"), onlyPionTracks.size());
      //_____________________________________
      // Creating rhos
      for (auto pion : onlyPionTracks) {
        p += pion;
      }
      //_____________________________________
      if (collision.numContrib() == 4) {
        // Four pions analysis
        if ((rawPionTracks.size() == 4) && (onlyPionTracks.size() == 4)) {

          registry.fill(HIST("hEta_t1"), onlyPionTracks[0].Eta());
          registry.fill(HIST("hEta_t2"), onlyPionTracks[1].Eta());

          /*if ((onlyPionSigma[0] * onlyPionSigma[0] + onlyPionSigma[1] * onlyPionSigma[1]) > 9) {
           return;
           }
           if ((onlyPionSigma[2] * onlyPionSigma[2] + onlyPionSigma[3] * onlyPionSigma[3]) > 9) {
           return;
           }*/

          // Quality Control
          registry.fill(HIST("hNsigPi1vsPi2"), rawPionTracks[0].tpcNSigmaPi(), rawPionTracks[1].tpcNSigmaPi());
          registry.fill(HIST("hNsigEl1vsEl2"), rawPionTracks[2].tpcNSigmaPi(), rawPionTracks[3].tpcNSigmaPi());
          registry.fill(HIST("hNsigPivsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaPi());
          registry.fill(HIST("hNsigPivsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaPi());
          registry.fill(HIST("hNsigElvsPt1"), onlyPionTracks[2].Pt(), rawPionTracks[2].tpcNSigmaPi());
          registry.fill(HIST("hNsigElvsPt2"), onlyPionTracks[3].Pt(), rawPionTracks[3].tpcNSigmaPi());

          registry.fill(HIST("hNsigMuvsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaMu());
          registry.fill(HIST("hNsigMuvsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaMu());

          registry.fill(HIST("hPtsingle_track1"), onlyPionTracks[0].Pt());
          registry.fill(HIST("hPtsingle_track2"), onlyPionTracks[1].Pt());
          registry.fill(HIST("hP1"), onlyPionTracks[0].P(), rawPionTracks[0].tpcSignal());
          registry.fill(HIST("hTPCsig"), rawPionTracks[0].tpcSignal(), rawPionTracks[1].tpcSignal());
          registry.fill(HIST("hP2"), onlyPionTracks[2].P(), rawPionTracks[2].tpcSignal());
          registry.fill(HIST("hTPCsig1"), rawPionTracks[2].tpcSignal(), rawPionTracks[3].tpcSignal());

          int sign = 0;
          TLorentzVector piplus, piminus;
          for (auto rawPion : rawPionTracks) {
            sign += rawPion.sign();
            if (rawPion.sign() > 0) {
              piplus = onlyPionTracks[0];
              piplus = onlyPionTracks[1];
              piplus = onlyPionTracks[2];
              piplus = onlyPionTracks[3];
            } else if (rawPion.sign() < 0) {
              piminus = onlyPionTracks[0];
              piminus = onlyPionTracks[1];
              piminus = onlyPionTracks[2];
              piminus = onlyPionTracks[3];
            }
          }

          auto costheta = CosThetaHelicityFrame(piplus, piminus, p);
          registry.fill(HIST("hCostheta4Pion"), costheta);
          auto phihel = 1. * TMath::Pi() + PhiHelicityFrame(piplus, piminus, p);
          registry.fill(HIST("hPhi4Pion"), phihel);

          if (p.Rapidity() > -0.9 && p.Rapidity() < 0.9) {

            if (sign != 0) {
              registry.fill(HIST("hMassFourPions"), p.M());
              registry.fill(HIST("hPt4Pion"), p.Pt());
            }

            if (sign == 0) {
              registry.fill(HIST("hCharge"), sign);
              registry.fill(HIST("h4TracksPions"), onlyPionTracks.size());
              if (p.Pt() < 0.15) {
                registry.fill(HIST("hMassFourPionsCoherent"), p.M());
              }
              // if ((p.M() > 0.8) && (p.M() < 2.5)) {
              registry.fill(HIST("hPt4PionRightSign"), p.Pt());
              // }
              registry.fill(HIST("hMassFourPionsRightSign"), p.M());
              registry.fill(HIST("hMPt1"), p.M(), p.Pt());
              registry.fill(HIST("hRap4Pion"), p.Rapidity());
              registry.fill(HIST("hEta4Pion"), p.Eta());

              for (auto pion : onlyPionTracks) {
                registry.fill(HIST("hPhiEtaFourPionsRightSign"), pion.Phi(), pion.Eta());
              }
            }
          }
        }
      }
      //_____________________________________
      // Six pions analysis
      if (collision.numContrib() == 6) {
        if ((rawPionTracks.size() == 6) && (onlyPionTracks.size() == 6)) {

          int sign = 0;
          TLorentzVector piplus, piminus;
          for (auto rawPion : rawPionTracks) {
            sign += rawPion.sign();
            if (rawPion.sign() > 0) {
              piplus = onlyPionTracks[0];
              piplus = onlyPionTracks[1];
              piplus = onlyPionTracks[2];
              piplus = onlyPionTracks[3];
              piplus = onlyPionTracks[4];
              piplus = onlyPionTracks[5];

            } else if (rawPion.sign() < 0) {
              piminus = onlyPionTracks[0];
              piminus = onlyPionTracks[1];
              piminus = onlyPionTracks[2];
              piminus = onlyPionTracks[3];
              piminus = onlyPionTracks[4];
              piminus = onlyPionTracks[5];
            }
          }

          auto costheta = CosThetaHelicityFrame(piplus, piminus, p);
          registry.fill(HIST("hCostheta6Pion"), costheta);
          auto phihel = 1. * TMath::Pi() + PhiHelicityFrame(piplus, piminus, p);
          registry.fill(HIST("hPhi6Pion"), phihel);

          if (p.Rapidity() > -0.9 && p.Rapidity() < 0.9) {
            if (sign != 0) {
              registry.fill(HIST("hMassSixPions"), p.M());
              registry.fill(HIST("hPt6Pion"), p.Pt());
            }

            if (sign == 0) {
              registry.fill(HIST("h6TracksPions"), onlyPionTracks.size());
              if (p.Pt() < 0.15) {
                registry.fill(HIST("hMassSixPionsCoherent"), p.M());
              }
              registry.fill(HIST("hMassSixPionsRightSign"), p.M());
              registry.fill(HIST("hMPt2"), p.M(), p.Pt());
              registry.fill(HIST("hRap6Pion"), p.Rapidity());
              registry.fill(HIST("hEta6Pion"), p.Eta());
              registry.fill(HIST("hPt6PionRightSign"), p.Pt());
              for (auto pion : onlyPionTracks) {
                registry.fill(HIST("hPhiEtaSixPionsRightSign"), pion.Phi(), pion.Eta());
              }
            }
          }
        }
      }
      //_____________________________________
      // Eight pions analysis
      if (collision.numContrib() == 8) {
        if ((rawPionTracks.size() == 8) && (onlyPionTracks.size() == 8)) {
          TLorentzVector piplus, piminus;
          int sign = 0;
          for (auto rawPion : rawPionTracks) {
            sign += rawPion.sign();

            if (rawPion.sign() > 0) {
              piplus = onlyPionTracks[0];
              piplus = onlyPionTracks[1];
              piplus = onlyPionTracks[2];
              piplus = onlyPionTracks[3];
              piplus = onlyPionTracks[4];
              piplus = onlyPionTracks[5];
              piplus = onlyPionTracks[6];
              piplus = onlyPionTracks[7];

            } else if (rawPion.sign() < 0) {
              piminus = onlyPionTracks[0];
              piminus = onlyPionTracks[1];
              piminus = onlyPionTracks[2];
              piminus = onlyPionTracks[3];
              piminus = onlyPionTracks[4];
              piminus = onlyPionTracks[5];
              piminus = onlyPionTracks[6];
              piminus = onlyPionTracks[7];
            }
          }

          auto costheta = CosThetaHelicityFrame(piplus, piminus, p);
          registry.fill(HIST("hCostheta8Pion"), costheta);
          auto phihel = 1. * TMath::Pi() + PhiHelicityFrame(piplus, piminus, p);
          registry.fill(HIST("hPhi8Pion"), phihel);

          if (p.Rapidity() > -0.9 && p.Rapidity() < 0.9) {
            if (sign != 0) {
              registry.fill(HIST("hMassEightPions"), p.M());
              registry.fill(HIST("hPt8Pion"), p.Pt());
            }
            if (sign == 0) {
              registry.fill(HIST("h8TracksPions"), onlyPionTracks.size());
              if (p.Pt() < 0.15) {
                registry.fill(HIST("hMassEightPionsCoherent"), p.M());
              }
              registry.fill(HIST("hMass8PionsRightSign"), p.M());
              registry.fill(HIST("hPt8PionRightSign"), p.Pt());
              registry.fill(HIST("hMPt3"), p.M(), p.Pt());
              registry.fill(HIST("hRap8pion"), p.Rapidity());
              registry.fill(HIST("hEta8Pion"), p.Eta());
              for (auto pion : onlyPionTracks) {
                registry.fill(HIST("hPhiEta8PionsRightSign"), pion.Phi(), pion.Eta());
              }
            }
          }
        }
      }
      //_____________________________________
      // CUTS
      if (collision.numContrib() == 2) {
        if (onlyPionTracks.size() != 2) {
          return;
        }
        registry.fill(HIST("h2TracksPions"), onlyPionTracks.size());
        if ((onlyPionSigma[0] * onlyPionSigma[0] + onlyPionSigma[1] * onlyPionSigma[1]) > 9) {
          return;
        }

        if (onlyPionTracks[0].P() < 0.3 || onlyPionTracks[1].P() < 0.3) {
          return;
        }
        if (p.Rapidity() < -0.9 || p.Rapidity() > 0.9) {
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

        auto costheta = CosThetaHelicityFrame(tplus, tminus, p);
        registry.fill(HIST("hCostheta"), costheta);
        auto phihel = 1. * TMath::Pi() + PhiHelicityFrame(tplus, tminus, p);
        registry.fill(HIST("hPhi"), phihel);

        registry.fill(HIST("hPt"), p.Pt());
        registry.fill(HIST("hEta"), p.Eta());
        registry.fill(HIST("hRap"), p.Rapidity());
        // registry.fill(HIST("hPhi"),  p.Phi());
        registry.fill(HIST("hMPt"), p.M(), p.Pt());

        if (p.Pt() < 0.2) {
          registry.fill(HIST("hMassUnlikeCoherent"), p.M());
        }
      }
    } // double gap
  }   // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCPionAnalysis>(cfgc)};
}
