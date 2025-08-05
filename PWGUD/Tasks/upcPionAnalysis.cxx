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
#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include <iostream>
#include "PWGUD/DataModel/UDTables.h"
#include <TString.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "Common/Core/RecoDecay.h"
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
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};

  // Collision selection
  Configurable<float> collcontrib_cut{"collcontrib_cut", 10, "no. of PV contributor per collsion"};
  Configurable<float> Zvtx_cut{"Zvtx_cut", 15, "z-vertex selection"};

  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  Configurable<float> TPC_cluster{"TPC_cluster", 50, "No.of TPC cluster"};

  // Kinmatic cuts
  Configurable<float> PID_cut{"PID_cut", 5, "TPC PID"};
  Configurable<float> Rap_cut{"Rap_cut", 0.9, "Track rapidity"};
  Configurable<float> Mass_Max{"Mass_Max", 10, "Invariant Mass range high"};
  Configurable<float> Mass_Min{"Mass_Min", 0, "Invariant Mass range low"};
  Configurable<float> Pt_coherent{"Pt_coherent", 0.15, "Coherent selection"};

  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //_____________________________________________________________________________________________
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
  //-----------------------------------------------------------------------------------------------
  float* AngCorrelation(TLorentzVector* lv1, TLorentzVector* lv2, TLorentzVector* lv)
  {
    TLorentzVector pa(1., 0., 0, 1);  // projectile
    TLorentzVector pb(1., 0., 0, -1); // target

    float* q = new float[3]; // to save values

    // Accoplanarity angle
    Float_t deltaPhi;
    deltaPhi = lv1->Phi() - lv2->Phi();
    float accoplCut = 1. - std::abs(deltaPhi) / o2::constants::math::PI;
    // z
    TLorentzVector z;
    Float_t part1, part2;

    // Dot product: v1*v2 = t1*t2-x1*x2-y1*y2-z1*z2

    part1 = lv->Dot(pb);
    part2 = lv->Dot(pa);

    Float_t part3x = pa.X() * part1;
    Float_t part3y = pa.Y() * part1;
    Float_t part3z = pa.Z() * part1;
    Float_t part3e = pa.T() * part1;

    Float_t part4x = pb.X() * part2;
    Float_t part4y = pb.Y() * part2;
    Float_t part4z = pb.Z() * part2;
    Float_t part4e = pb.T() * part2;

    TLorentzVector part3(TVector3(part3x, part3y, part3z), part3e);
    TLorentzVector part4(TVector3(part4x, part4y, part4z), part4e);

    // Un-normalized Z
    z = part3 - part4;

    // Normalized z
    Float_t normz = std::sqrt(-z * z);
    Float_t znx = z.X() / normz;
    Float_t zny = z.Y() / normz;
    Float_t znz = z.Z() / normz;
    Float_t zne = z.E() / normz;

    // Normalized z
    TLorentzVector zhat(TVector3(znx, zny, znz), zne);

    // calculate x
    TLorentzVector x;

    Float_t constant1 = (lv->Dot(*lv)) / (2 * (lv->Dot(pa)));
    Float_t constant2 = (lv->Dot(*lv)) / (2 * (lv->Dot(pb)));

    Float_t comp1x = pa.X() * constant1;
    Float_t comp1y = pa.Y() * constant1;
    Float_t comp1z = pa.Z() * constant1;
    Float_t comp1e = pa.T() * constant1;

    TLorentzVector comp1(TVector3(comp1x, comp1y, comp1z), comp1e);

    Float_t comp2x = pb.X() * constant2;
    Float_t comp2y = pb.Y() * constant2;
    Float_t comp2z = pb.Z() * constant2;
    Float_t comp2e = pb.T() * constant2;

    TLorentzVector comp2(TVector3(comp2x, comp2y, comp2z), comp2e);

    // Un-normalized x
    x = *lv - comp1 - comp2;
    // normalize x
    Float_t normx = std::sqrt(-x * x);
    Float_t xnx = x.X() / normx;
    Float_t xny = x.Y() / normx;
    Float_t xnz = x.Z() / normx;
    Float_t xne = x.E() / normx;

    // Normalized x
    TLorentzVector xhat(TVector3(xnx, xny, xnz), xne);

    // calculate y
    // TLorentzVector y;
    Float_t yone = pa.Y() * pb.Z() * lv->E() - pa.Z() * pb.Y() * lv->E() + pa.Z() * pb.E() * lv->Y() + pa.E() * pb.Y() * lv->Z() - pa.Y() * pb.E() * lv->Z() - pa.E() * pb.Z() * lv->Y();
    Float_t ytwo = -pa.Z() * pb.E() * lv->X() + pa.Z() * pb.X() * lv->E() - pa.X() * pb.Z() * lv->E() + pa.X() * pb.E() * lv->Z() - pa.E() * pb.X() * lv->Z() + pa.E() * pb.Z() * lv->X();
    Float_t ythree = pa.X() * pb.Y() * lv->E() - pa.Y() * pb.X() * lv->E() + pa.Y() * pb.E() * lv->X() - pa.X() * pb.E() * lv->Y() + pa.E() * pb.X() * lv->Y() - pa.E() * pb.Y() * lv->X();
    Float_t yfour = -pa.X() * pb.Y() * lv->Z() + pa.X() * pb.Z() * lv->Y() - pa.Z() * pb.X() * lv->Y() + pa.Z() * pb.Y() * lv->X() - pa.Y() * pb.Z() * lv->X() + pa.Y() * pb.X() * lv->Z();

    // Un-normalized y
    TLorentzVector y(TVector3(yone, ytwo, ythree), yfour);

    // normalize y
    Float_t normy = std::sqrt(-y * y);
    Float_t ynx = y.X() / normy;
    Float_t yny = y.Y() / normy;
    Float_t ynz = y.Z() / normy;
    Float_t yne = y.E() / normy;

    // normalized y
    TLorentzVector yhat(TVector3(ynx, yny, ynz), yne);

    // Lepton momentum difference
    TLorentzVector diff;
    diff = (*lv1 - *lv2);
    Float_t diff2x = diff.X() / 2.;
    Float_t diff2y = diff.Y() / 2.;
    Float_t diff2z = diff.Z() / 2.;
    Float_t diff2e = diff.E() / 2.;
    TLorentzVector diff2(TVector3(diff2x, diff2y, diff2z), diff2e);

    // Normalize diff2
    Float_t norm2 = std::sqrt(-diff2 * diff2);
    Float_t diff3x = diff2.X() / norm2;
    Float_t diff3y = diff2.Y() / norm2;
    Float_t diff3z = diff2.Z() / norm2;
    Float_t diff3e = diff2.E() / norm2;

    TLorentzVector diff3(TVector3(diff3x, diff3y, diff3z), diff3e);

    // computing the angles
    float cosThetaCS = zhat * diff3;
    Double_t SinThetaCosPhiCS = xhat * diff3;
    Double_t SinThetaSinPhiCS = yhat * diff3;
    //**************************************

    float phi = atan2(SinThetaSinPhiCS, SinThetaCosPhiCS);
    // if (phi>=0) phi = phi;
    if (phi < 0)
      phi = phi + o2::constants::math::TwoPI;

    q[0] = accoplCut;
    q[1] = phi;
    q[2] = cosThetaCS;

    return q;
  }

  double DeltaPhi(TLorentzVector lv1, TLorentzVector lv2)
  {
    TLorentzVector lv_sum = lv1 + lv2;
    TLorentzVector lv_diff = lv1 - lv2;

    double dp = lv_sum.DeltaPhi(lv_diff);

    return dp;
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

    registry.add("posx", "Vertex position in x", kTH1F, {{100, -0.5, 0.5}});
    registry.add("posy", "Vertex position in y", kTH1F, {{100, -0.5, 0.5}});
    registry.add("posz", "Vertex position in z", kTH1F, {{1000, -100., 100.}});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{20, 0., 20.}});

    TString SelectionCuts[18] = {"NoSelection", "gapside", "goodtracks", "truegap", "ncollcontrib ", "zvtx", "2collcontrib", "2goodtrk", "TPCPID", "Rap_cut", "unlikesign", "mass_cut", "coherent", "incoherent", "likesign", "mass_cut", "coherent", "incoherent"};
    // now we can set BinLabel in histogram Registry

    for (int i = 0; i < 18; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }
    // tracks
    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("TwoPion/h2TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("h4TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("h6TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("h8TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});

    registry.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});
    registry.add("hdEdxPion", "p_{#pi} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});

    registry.add("TwoPion/hNsigPi1vsPi2", "NSigmaPi(t1) vs NSigmapi (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("TwoPion/hNsigEl1vsEl2", "NSigmaEl(t1) vs NSigmaEl (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("TwoPion/hNsigPivsPt1", "Pt vs NSigmaPi (t1);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("TwoPion/hNsigPivsPt2", "Pt vs NSigmaPi (t2);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("TwoPion/hNsigElvsPt1", "Pt vs NSigmaEl (t1);#it{p_{t}}, GeV/c;n#sigma_{#e}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("TwoPion/hNsigElvsPt2", "Pt vs NSigmaEl (t2);#it{p_{t}}, GeV/c;n#sigma_{#e}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("TwoPion/hNsigMuvsPt1", "Pt vs NSigmaMu (t1);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("TwoPion/hNsigMuvsPt2", "Pt vs NSigmaMu (t2);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5}, {100, -15., 15}});
    registry.add("trackpt", "Track pT per track; Track pT; Track bit; Tracks", {HistType::kTH2F, {{4, 0.5, 4.5}, {100, 0.0, 10.}}});
    registry.add("tracketa", "Track eta per track; Track eta; Track bit; Tracks", {HistType::kTH2F, {{4, 0.5, 4.5}, {100, -5, 5.}}});
    registry.add("trackphi", "Track phi per track; Track phi; Track bit; Tracks", {HistType::kTH2F, {{4, 0.5, 4.5}, {100, 0. * TMath::Pi(), 2. * TMath::Pi()}}});
    registry.add("tracksign", "Track sign per track; Track sign; Track bit; Tracks", {HistType::kTH2F, {{4, 0.5, 4.5}, {200, -10., 10.}}});

    registry.add("hCharge", "Charge;#it{charge};", kTH1F, {{500, -10., 10.}});
    registry.add("TwoPion/hPtsingle_track1", "Pt t1;#it{p_{t}}, GeV/c;", kTH1F, {{600, 0., 3.}});
    registry.add("TwoPion/hPtsingle_track2", "Pt t2;#it{p_{t}}, GeV/c;", kTH1F, {{600, 0., 3.}});
    registry.add("TwoPion/hEta_t1", "Eta of t1;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("TwoPion/hEta_t2", "Eta of t2;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("TwoPion/hP1", "P vs TPC signal;#it{P_{track}}, GeV/c; signal_{TPC} t1", kTH2F, {{100, 0., 2.}, {300, 0, 150}});
    registry.add("TwoPion/hTPCsig", "TPC signal;signal_{TPC} t2; signal_{TPC} t2", kTH2F, {{300, 0., 150.}, {300, 0, 150}});
    registry.add("TwoPion/hP2", "P vs TPC signal;#it{P_{track}}, GeV/c; signal_{TPC} t1", kTH2F, {{100, 0., 2.}, {300, 0, 150}});
    registry.add("TwoPion/hTPCsig1", "TPC signal;signal_{TPC} t2; signal_{TPC} t2", kTH2F, {{300, 0., 150.}, {300, 0, 150}});

    registry.add("TwoPion/hMassLike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{20000, 0., 20.}});
    registry.add("TwoPion/hMassUnlike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{20000, 0., 20.}});
    registry.add("TwoPion/Coherent/hMassUnlikeCoherent", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{20000, 0., 20.}});
    registry.add("TwoPion/Coherent/hMassLikeCoherent", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{20000, 0., 20.}});
    registry.add("TwoPion/Incoherent/hMassUnlikeInCoherent", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{20000, 0., 20.}});
    registry.add("TwoPion/Incoherent/hMassLikeInCoherent", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{20000, 0., 20.}});
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

    registry.add("TwoPion/hPt", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("TwoPion/hPtLike", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt4Pion", "Pt1;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt6Pion", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt8Pion", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt4PionRightSign", "Pt1;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt6PionRightSign", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("hPt8PionRightSign", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});

    registry.add("TwoPion/hEta", "Eta;#it{#eta};", kTH1F, {{500, -10., 10.}});
    registry.add("hEta4Pion", "Eta of four pion;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("hEta6Pion", "Eta of six pion;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("hEta8Pion", "Eta of six pion;#it{#eta};", kTH1F, {{100, -5., 5.}});

    registry.add("TwoPion/hRap", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
    registry.add("hRap4Pion", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
    registry.add("hRap6Pion", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
    registry.add("hRap8pion", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});

    registry.add("TwoPion/hPhiSystem", "Phi;#it{#Phi};", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});

    registry.add("TwoPion/hMPt", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 10.}, {100, 0., 10.}});
    registry.add("hMPt1", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 10.}, {100, 0., 10.}});
    registry.add("hMPt2", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 10.}, {100, 0., 10.}});
    registry.add("hMPt3", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 20.}, {100, 0., 10.}});

    // Uisng Polarisaion method
    registry.add("TwoPion/Coherent/hPhi", "Phi;#it{#Phi};", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("TwoPion/Incoherent/hPhiInCoherent", "Phi;#it{#Phi};", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("hPhi4Pion", "Phi;#it{#Phi};", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("hPhi6Pion", "Phi;#it{#Phi};", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("hPhi8Pion", "Phi;#it{#Phi};", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});

    registry.add("TwoPion/Coherent/hCostheta", "Costheta;#it{Cos#Theta};", kTH1F, {{300, -1.5, 1.5}});
    registry.add("TwoPion/Incoherent/hCosthetaInCoherent", "Costheta;#it{Cos#Theta};", kTH1F, {{300, -1.5, 1.5}});
    registry.add("hCostheta4Pion", "Costheta;#it{Cos#Theta};", kTH1F, {{300, -1.5, 1.5}});
    registry.add("hCostheta6Pion", "Costheta;#it{Cos#Theta};", kTH1F, {{300, -1.5, 1.5}});
    registry.add("hCostheta8Pion", "Costheta;#it{Cos#Theta};", kTH1F, {{300, -1.5, 1.5}});

    // Using Angular Correlation method

    registry.add("TwoPion/Coherent/AccoplAngle", "AccoplAngle", kTH1F, {{250, -0.2, 0.2}});
    registry.add("TwoPion/Coherent/CosTheta", "CosTheta", kTH1F, {{300, -1.5, 1.5}});
    registry.add("TwoPion/Coherent/Phi", "Phi", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("TwoPion/Coherent/Phi1", "Phi1", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("TwoPion/Coherent/Phi2", "Phi2", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("TwoPion/Coherent/CosThetaPhi", "CosThetaPhi", kTH2F, {{300, -1.5, 1.5}, {250, 0. * TMath::Pi(), 2. * TMath::Pi()}});

    registry.add("TwoPion/Incoherent/AccoplAngle", "AccoplAngle", kTH1F, {{250, -0.2, 0.2}});
    registry.add("TwoPion/Incoherent/CosTheta", "CosTheta", kTH1F, {{300, -1.5, 1.5}});
    registry.add("TwoPion/Incoherent/Phi", "Phi", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("TwoPion/Incoherent/Phi1", "Phi1", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("TwoPion/Incoherent/Phi2", "Phi2", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("TwoPion/Incoherent/CosThetaPhi", "CosThetaPhi", kTH2F, {{300, -1.5, 1.5}, {250, 0., 2. * TMath::Pi()}});

    // Asymmetry histograms
    registry.add("TwoPion/Coherent/DeltaPhi", "DeltaPhi", kTH1F, {{360, -TMath::Pi(), TMath::Pi()}});
    registry.add("TwoPion/Incoherent/DeltaPhi", "DeltaPhi", kTH1F, {{360, -TMath::Pi(), TMath::Pi()}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  // Main process
  void process(UDCollisionsFull::iterator const& collision, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);
    // LOGF(info, " BC ID %d",collision.gapSide());
    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;

    registry.fill(HIST("hSelectionCounter"), 1);

    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    // int truegapSide = sgSelector.trueGap(collision, FV0_cut, ZDC_cut);
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    registry.fill(HIST("hSelectionCounter"), 2);

    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;

    if (gapSide == gap_Side) {

      registry.fill(HIST("hSelectionCounter"), 3);
      //_____________________________________
      // Create pions and apply TPC Pion PID
      std::vector<TLorentzVector> allTracks;
      std::vector<TLorentzVector> onlyPionTracks;
      std::vector<float> onlyPionSigma;
      std::vector<decltype(tracks.begin())> rawPionTracks;
      std::vector<float> trackpt;
      std::vector<float> tracketa;
      std::vector<float> trackphi;
      std::vector<float> tracksign;
      TLorentzVector p;
      registry.fill(HIST("hTracks"), tracks.size());

      if (collision.numContrib() > collcontrib_cut)
        return;

      registry.fill(HIST("hSelectionCounter"), 4);
      if ((collision.posZ() < -(Zvtx_cut)) || (collision.posZ() > Zvtx_cut))
        return;
      registry.fill(HIST("hSelectionCounter"), 5);

      for (auto t : tracks) {

        /*if (!t.isPVContributor()) {
          continue;
        }*/

        if (!trackselector(t, parameters))
          continue;

        // int NFindable = t.tpcNClsFindable();
        // int NMinusFound = t.tpcNClsFindableMinusFound();
        // int NCluster = NFindable - NMinusFound;

        /*if (NCluster < TPC_cluster) {
           continue;
         }*/

        double dEdx = t.tpcSignal();

        registry.fill(HIST("hdEdx"), t.tpcInnerParam() / t.sign(), dEdx);
        TLorentzVector a;
        a.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassPionCharged);
        allTracks.push_back(a);
        auto nSigmaPi = t.tpcNSigmaPi();

        if (fabs(nSigmaPi) < PID_cut) {
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

          registry.get<TH2>(HIST("trackpt"))->Fill(1., onlyPionTracks[0].Pt());
          registry.get<TH2>(HIST("trackpt"))->Fill(2., onlyPionTracks[1].Pt());
          registry.get<TH2>(HIST("trackpt"))->Fill(3., onlyPionTracks[2].Pt());
          registry.get<TH2>(HIST("trackpt"))->Fill(4., onlyPionTracks[3].Pt());

          registry.get<TH2>(HIST("tracketa"))->Fill(1., onlyPionTracks[0].Eta());
          registry.get<TH2>(HIST("tracketa"))->Fill(2., onlyPionTracks[1].Eta());
          registry.get<TH2>(HIST("tracketa"))->Fill(3., onlyPionTracks[2].Eta());
          registry.get<TH2>(HIST("tracketa"))->Fill(4., onlyPionTracks[3].Eta());

          registry.get<TH2>(HIST("trackphi"))->Fill(1., onlyPionTracks[0].Phi());
          registry.get<TH2>(HIST("trackphi"))->Fill(2., onlyPionTracks[1].Phi());
          registry.get<TH2>(HIST("trackphi"))->Fill(3., onlyPionTracks[2].Phi());
          registry.get<TH2>(HIST("trackphi"))->Fill(4., onlyPionTracks[3].Phi());

          registry.get<TH2>(HIST("tracksign"))->Fill(1., rawPionTracks[0].sign());
          registry.get<TH2>(HIST("tracksign"))->Fill(2., rawPionTracks[1].sign());
          registry.get<TH2>(HIST("tracksign"))->Fill(3., rawPionTracks[2].sign());
          registry.get<TH2>(HIST("tracksign"))->Fill(4., rawPionTracks[3].sign());

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

          if ((p.Rapidity() > -Rap_cut) && (p.Rapidity() < Rap_cut)) {
            if ((p.M() > 0.8) && (p.M() < 2.5)) {
              if (sign != 0) {
                registry.fill(HIST("hMassFourPions"), p.M());
                registry.fill(HIST("hPt4Pion"), p.Pt());
              }

              if (sign == 0) {
                registry.fill(HIST("hCharge"), sign);
                registry.fill(HIST("h4TracksPions"), onlyPionTracks.size());

                if (p.Pt() < Pt_coherent) {
                  registry.fill(HIST("hMassFourPionsCoherent"), p.M());
                }

                registry.fill(HIST("hPt4PionRightSign"), p.Pt());

                registry.fill(HIST("hMassFourPionsRightSign"), p.M());
                registry.fill(HIST("hMPt1"), p.M(), p.Pt());
                registry.fill(HIST("hRap4Pion"), p.Rapidity());
                registry.fill(HIST("hEta4Pion"), p.Eta());
              }
              for (auto pion : onlyPionTracks) {
                registry.fill(HIST("hPhiEtaFourPionsRightSign"), pion.Phi(), pion.Eta());
              }
            }
          }
        }
      }
      //_____________________________________________________________________________________________________
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

          if ((p.Rapidity() > -Rap_cut) && (p.Rapidity() < Rap_cut)) {
            if (sign != 0) {
              registry.fill(HIST("hMassSixPions"), p.M());
              registry.fill(HIST("hPt6Pion"), p.Pt());
            }

            if (sign == 0) {
              registry.fill(HIST("h6TracksPions"), onlyPionTracks.size());
              if (p.Pt() < Pt_coherent) {
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

          if ((p.Rapidity() > -Rap_cut) && (p.Rapidity() < Rap_cut)) {
            if (sign != 0) {
              registry.fill(HIST("hMassEightPions"), p.M());
              registry.fill(HIST("hPt8Pion"), p.Pt());
            }
            if (sign == 0) {
              registry.fill(HIST("h8TracksPions"), onlyPionTracks.size());
              if (p.Pt() < Pt_coherent) {
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
        registry.fill(HIST("hSelectionCounter"), 6);

        if ((rawPionTracks.size() == 2) && (onlyPionTracks.size() == 2)) {
          registry.fill(HIST("hSelectionCounter"), 7);

          registry.fill(HIST("TwoPion/h2TracksPions"), onlyPionTracks.size());
          registry.fill(HIST("TwoPion/hNsigPivsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaPi());
          registry.fill(HIST("TwoPion/hNsigPivsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaPi());
          registry.fill(HIST("TwoPion/hTPCsig"), rawPionTracks[0].tpcSignal(), rawPionTracks[1].tpcSignal());
          registry.fill(HIST("TwoPion/hNsigPi1vsPi2"), rawPionTracks[0].tpcNSigmaPi(), rawPionTracks[1].tpcNSigmaPi());

          if ((onlyPionSigma[0] * onlyPionSigma[0] + onlyPionSigma[1] * onlyPionSigma[1]) > (PID_cut * PID_cut)) {
            return;
          }
          registry.fill(HIST("hSelectionCounter"), 8);
          if ((p.Rapidity() < -Rap_cut) || (p.Rapidity() > Rap_cut)) {
            return;
          }
          registry.fill(HIST("hSelectionCounter"), 9);

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
          auto phihel = 1. * TMath::Pi() + PhiHelicityFrame(tplus, tminus, p);

          // Unlike sign
          if (rawPionTracks[0].sign() != rawPionTracks[1].sign()) {

            registry.fill(HIST("hSelectionCounter"), 10);
            registry.fill(HIST("TwoPion/hMassUnlike"), p.M());

            if ((p.M() > Mass_Min) && (p.M() < Mass_Max)) {

              registry.fill(HIST("hSelectionCounter"), 11);

              registry.fill(HIST("TwoPion/hPt"), p.Pt());
              registry.fill(HIST("TwoPion/hEta"), p.Eta());
              registry.fill(HIST("TwoPion/hRap"), p.Rapidity());
              registry.fill(HIST("TwoPion/hPhiSystem"), p.Phi());
              registry.fill(HIST("TwoPion/hMPt"), p.M(), p.Pt());

              if (p.Pt() < Pt_coherent) {

                registry.fill(HIST("hSelectionCounter"), 12);

                // Quality Control
                registry.fill(HIST("TwoPion/hEta_t1"), onlyPionTracks[0].Eta());
                registry.fill(HIST("TwoPion/hEta_t2"), onlyPionTracks[1].Eta());
                registry.fill(HIST("TwoPion/hPtsingle_track1"), onlyPionTracks[0].Pt());
                registry.fill(HIST("TwoPion/hPtsingle_track2"), onlyPionTracks[1].Pt());

                registry.fill(HIST("TwoPion/hNsigMuvsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaPi());
                registry.fill(HIST("TwoPion/hNsigMuvsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaPi());
                registry.fill(HIST("TwoPion/hNsigElvsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaEl());
                registry.fill(HIST("TwoPion/hNsigElvsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaEl());
                registry.fill(HIST("TwoPion/hNsigEl1vsEl2"), rawPionTracks[0].tpcNSigmaPi(), rawPionTracks[1].tpcNSigmaPi());

                registry.fill(HIST("TwoPion/hP1"), onlyPionTracks[0].P(), rawPionTracks[0].tpcSignal());
                registry.fill(HIST("TwoPion/hP2"), onlyPionTracks[1].P(), rawPionTracks[1].tpcSignal());
                registry.fill(HIST("TwoPion/hTPCsig1"), rawPionTracks[0].tpcSignal(), rawPionTracks[1].tpcSignal());

                registry.fill(HIST("posx"), collision.posX());
                registry.fill(HIST("posy"), collision.posY());
                registry.fill(HIST("posz"), collision.posZ());

                registry.fill(HIST("TwoPion/Coherent/hCostheta"), costheta);
                registry.fill(HIST("TwoPion/Coherent/hPhi"), phihel);
                registry.fill(HIST("TwoPion/Coherent/hMassUnlikeCoherent"), p.M());

                // angular correlations
                float* q = AngCorrelation(&onlyPionTracks[0], &onlyPionTracks[1], &p);
                registry.fill(HIST("TwoPion/Coherent/Phi"), q[1], 1.);
                registry.fill(HIST("TwoPion/Coherent/Phi1"), RecoDecay::phi(onlyPionTracks[0]), 1.);
                registry.fill(HIST("TwoPion/Coherent/Phi2"), RecoDecay::phi(onlyPionTracks[1]), 1.);
                registry.fill(HIST("TwoPion/Coherent/CosTheta"), q[2], 1.);
                registry.fill(HIST("TwoPion/Coherent/CosThetaPhi"), q[2], q[1]);
                registry.fill(HIST("TwoPion/Coherent/AccoplAngle"), q[0], 1.);

                double delphi = DeltaPhi(onlyPionTracks[0], onlyPionTracks[1]);
                registry.fill(HIST("TwoPion/Coherent/DeltaPhi"), delphi);
              }

              if (p.Pt() > Pt_coherent) {

                registry.fill(HIST("hSelectionCounter"), 13);

                registry.fill(HIST("TwoPion/Incoherent/hCosthetaInCoherent"), costheta);
                registry.fill(HIST("TwoPion/Incoherent/hPhiInCoherent"), phihel);
                registry.fill(HIST("TwoPion/Incoherent/hMassUnlikeInCoherent"), p.M());

                // angular correlations
                float* q = AngCorrelation(&onlyPionTracks[0], &onlyPionTracks[1], &p);
                registry.fill(HIST("TwoPion/Incoherent/Phi"), q[1], 1.);
                registry.fill(HIST("TwoPion/Incoherent/Phi1"), RecoDecay::phi(onlyPionTracks[0]), 1.);
                registry.fill(HIST("TwoPion/Incoherent/Phi2"), RecoDecay::phi(onlyPionTracks[1]), 1.);
                registry.fill(HIST("TwoPion/Incoherent/CosTheta"), q[2], 1.);
                registry.fill(HIST("TwoPion/Incoherent/CosThetaPhi"), q[2], q[1]);
                registry.fill(HIST("TwoPion/Incoherent/AccoplAngle"), q[0], 1.);

                double delphi = DeltaPhi(onlyPionTracks[0], onlyPionTracks[1]);
                registry.fill(HIST("TwoPion/Incoherent/DeltaPhi"), delphi);
              }
            }
          }

          // Like sign
          if (rawPionTracks[0].sign() == rawPionTracks[1].sign()) {

            registry.fill(HIST("hSelectionCounter"), 14);
            registry.fill(HIST("TwoPion/hMassLike"), p.M());
            if ((p.M() > Mass_Min) && (p.M() < Mass_Max)) {
              registry.fill(HIST("hSelectionCounter"), 15);
              registry.fill(HIST("TwoPion/hPtLike"), p.Pt());
              if (p.Pt() < Pt_coherent) {
                registry.fill(HIST("hSelectionCounter"), 16);
                registry.fill(HIST("TwoPion/Coherent/hMassLikeCoherent"), p.M());
              }
              if (p.Pt() > Pt_coherent) {
                registry.fill(HIST("hSelectionCounter"), 17);
                registry.fill(HIST("TwoPion/Incoherent/hMassLikeInCoherent"), p.M());
              }
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCPionAnalysis>(cfgc)};
}
