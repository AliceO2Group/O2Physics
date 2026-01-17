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
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TLorentzVector.h"
#include <TString.h>

#include <iostream>
using std::array;
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// \brief Exclusive proton+proton task
/// \author Simone Ragoni, Creighton
/// \date 29/4/2024

struct ExclusiveTwoProtonsSG {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 200., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 100., "FT0C threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};
  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  //_____________________________________________________________________________
  Double_t CosThetaHelicityFrame(TLorentzVector posDaughter, TLorentzVector negDaughter, TLorentzVector mother)
  {

    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    TVector3 beta = (-1. / mother.E()) * mother.Vect();
    TLorentzVector pPi1Dipion = posDaughter;
    TLorentzVector pPi2Dipion = negDaughter;
    TLorentzVector pProjDipion = pProjCM;
    TLorentzVector pTargDipion = pTargCM;

    pPi1Dipion.Boost(beta);
    pPi2Dipion.Boost(beta);
    pProjDipion.Boost(beta);
    pTargDipion.Boost(beta);

    TVector3 zaxis = (mother.Vect()).Unit();

    Double_t CosThetaHE = zaxis.Dot((pPi1Dipion.Vect()).Unit());
    return CosThetaHE;
  }
  //------------------------------------------------------------------------------------------------------
  Double_t PhiHelicityFrame(TLorentzVector posDaughter, TLorentzVector negDaughter, TLorentzVector mother)
  {

    // Half of the energy per pair of the colliding nucleons.
    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    // Translate the dimuon parameters in the dimuon rest frame
    TVector3 beta = (-1. / mother.E()) * mother.Vect();
    TLorentzVector pMu1Dimu = posDaughter;
    TLorentzVector pMu2Dimu = negDaughter;
    TLorentzVector pProjDimu = pProjCM;
    TLorentzVector pTargDimu = pTargCM;
    pMu1Dimu.Boost(beta);
    pMu2Dimu.Boost(beta);
    pProjDimu.Boost(beta);
    pTargDimu.Boost(beta);

    // Axes
    TVector3 zaxis = (mother.Vect()).Unit();
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
    registry.add("posx", "Vertex position in x", kTH1F, {{100, -0.5, 0.5}});
    registry.add("posy", "Vertex position in y", kTH1F, {{100, -0.5, 0.5}});
    registry.add("posz", "Vertex position in z", kTH1F, {{1000, -100., 100.}});
    registry.add("hITSCluster", "N_{cluster}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hChi2ITSTrkSegment", "N_{cluster}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTPCCluster", "N_{cluster}", kTH1F, {{200, -0.5, 199.5}});
    registry.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdx2", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdSigmaProton", "p vs dE/dx sigma proton", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hdSigmaProton2", "p vs dE/dx sigma proton", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hdSigmaProton3", "p vs dE/dx sigma proton", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hNsigEvsKa1", "NSigma(t1) vs NSigma (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigEvsKa2", "NSigma(t1) vs NSigma (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigEvsKa3", "NSigma(t1) vs NSigma (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hMomentum", "p_{#ka};#it{p_{trk}}, GeV/c;", kTH1F, {{100, 0., 3.}});
    registry.add("hClusterSizeAllTracks", "ClusterSizeAllTracks;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeProtonsTPC", "ClusterSizeProtonsTPC;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeProtonsTOF", "ClusterSizeProtonsTOF;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hEta1", "#eta_{#ka};#it{#eta_{trk}}, GeV/c;", kTH1F, {{100, -2., 2.}});
    registry.add("hPtLikeSignProton", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("hMassLikeSignProton", "Raw Inv.M;#it{m_{pp}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("hMassPtLikeSignProton", "Raw Inv.M;#it{m_{PP}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{1000, 0., 10.}, {400, 0., 4.}});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{10, 0., 10.}});

    TString SelectionCuts[9] = {"NoSelection", "GAPcondition", "PVtracks", "Good TPC-ITS track", "TPC/TOF PID track", "End trk loop", "Exactly 2p", "Like-sign ev", "Unlike-sign ev"};
    // now we can set BinLabel in histogram Registry

    for (int i = 0; i < 9; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    // Unlike sign pp
    registry.add("PP/hRapidity", "Rapidity;#it{y_{pp}};", kTH1F, {{100, -2., 2.}});
    registry.add("PP/hPtProtonVsProton", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("PP/hMassPtUnlikeSignProton", "Raw Inv.M;#it{m_{pp}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PP/hMassUnlike", "#it{m_{pp}} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hMassHighPt0150", "#it{m_{pp}} [GeV/#it{c}^{2}]", kTH1F, {{2000, 0., 20.}});
    registry.add("PP/hMassHighPt0200", "#it{m_{pp}} [GeV/#it{c}^{2}]", kTH1F, {{2000, 0., 20.}});
    registry.add("PP/hMassHighPt0300", "#it{m_{pp}} [GeV/#it{c}^{2}]", kTH1F, {{2000, 0., 20.}});
    registry.add("PP/hMassHighPt0500", "#it{m_{pp}} [GeV/#it{c}^{2}]", kTH1F, {{2000, 0., 20.}});
    registry.add("PP/hMassHighPt1000", "#it{m_{pp}} [GeV/#it{c}^{2}]", kTH1F, {{2000, 0., 20.}});
    registry.add("PP/hUnlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hCoherentMass", "Raw Inv.M;#it{m_{pp}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hCoherentMassWithoutDCAcuts", "Raw Inv.M;#it{m_{pp}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hCoherentMassWithoutDCAxycut", "Raw Inv.M;#it{m_{pp}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hCoherentMassWithoutDCAzcut", "Raw Inv.M;#it{m_{pp}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hUnlikePtLowBand", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePtLowBand24to275", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePtLowBand275to300", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePtHighBand", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePtJpsi", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});

    registry.add("PP/hAngularDstribLab", "Angular distrib in the lab system;#it{#varphi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {100, -2, 2}});
    registry.add("PP/hCosThetaLab", "Angular distrib in the lab system;#it{Cos#Theta};", kTH1F, {{100, -2, 2}});
    registry.add("PP/hAngularDistribHelicity", "Angular distrib helicity frame;#it{#varphi};#it{Cos#Theta};", kTH2F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}, {100, -2, 2}});
    registry.add("PP/hCosThetaHelicity", "Angular distrib helicity frame;#it{Cos#Theta};", kTH1F, {{100, -2, 2}});
    registry.add("PP/hPhiHelicity", "Angular distrib helicity frame;#it{#varphi};", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  // Main process
  void process(UDCollisionsFull::iterator const& collision, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);
    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;

    registry.fill(HIST("posx"), collision.posX());
    registry.fill(HIST("posy"), collision.posY());
    registry.fill(HIST("posz"), collision.posZ());
    int truegapSide = sgSelector.trueGap(collision, FV0_cut, FT0A_cut, FT0C_cut, ZDC_cut);
    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    if (gapSide != gap_Side) {
      return;
    }
    TLorentzVector resonance; // lorentz vectors of tracks and the mother

    // ===================================
    // Task for p+pbar pairs with PID
    // ===================================
    std::vector<TLorentzVector> onlyProtonTracks;
    std::vector<TLorentzVector> onlyProtonTracksTOF;
    std::vector<float> onlyProtonSigma;
    std::vector<float> onlyProtonSigmaTOF;
    std::vector<decltype(tracks.begin())> rawProtonTracks;
    std::vector<decltype(tracks.begin())> rawProtonTracksTOF;

    for (auto trk : tracks) {
      if (!trk.isPVContributor()) {
        continue;
      }
      registry.fill(HIST("hSelectionCounter"), 2);

      int NFindable = trk.tpcNClsFindable();
      int NMinusFound = trk.tpcNClsFindableMinusFound();
      int NCluster = NFindable - NMinusFound;
      registry.fill(HIST("hTPCCluster"), NCluster);
      registry.fill(HIST("hChi2ITSTrkSegment"), trk.itsChi2NCl());
      if (NCluster < 70) {
        continue;
      }
      if (trk.itsNCls() < 7) {
        continue;
      }
      if (trk.pt() < 0.7) {
        continue;
      }
      if (!(std::abs(trk.dcaZ()) < 2.)) {
        continue;
      }
      double dcaLimit = 0.0105 + 0.035 / pow(trk.pt(), 1.1);
      if (!(std::abs(trk.dcaXY()) < dcaLimit)) {
        continue;
      }

      registry.fill(HIST("hSelectionCounter"), 3);
      registry.fill(HIST("hITSCluster"), trk.itsNCls());

      int clusterSize[7];
      double averageClusterSize = 0.;
      for (int i = 0; i < 7; i++) { // info stored in 4 bits
        clusterSize[i] = (((1 << 4) - 1) & (trk.itsClusterSizes() >> 4 * i));
        averageClusterSize += static_cast<double>(clusterSize[i]);
      }
      averageClusterSize /= 7.;
      registry.fill(HIST("hClusterSizeAllTracks"), averageClusterSize);

      double momentum = TMath::Sqrt(trk.px() * trk.px() + trk.py() * trk.py() + trk.pz() * trk.pz());
      double dEdx = trk.tpcSignal();
      registry.fill(HIST("hdEdx"), momentum, dEdx);

      TLorentzVector proton;
      proton.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassProton);
      TLorentzVector protonbar;
      protonbar.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassProtonBar);
      auto nSigmaPr = trk.tpcNSigmaPr();
      auto nSigmaPrTOF = trk.tofNSigmaPr();

      if (trk.hasTOF()) {
        registry.fill(HIST("hdSigmaProton"), momentum, nSigmaPrTOF);
      }
      if (fabs(nSigmaPr) < 4.) {
        registry.fill(HIST("hdEdx2"), momentum, dEdx);
        registry.fill(HIST("hdSigmaProton2"), momentum, nSigmaPr);
        registry.fill(HIST("hMomentum"), momentum);
        if (trk.hasTOF() && fabs(nSigmaPrTOF) < 4.) {
          registry.fill(HIST("hdSigmaProton3"), momentum, nSigmaPrTOF);
          registry.fill(HIST("hClusterSizeProtonsTOF"), averageClusterSize);
          onlyProtonTracksTOF.push_back(proton);
          onlyProtonSigmaTOF.push_back(nSigmaPrTOF);
          rawProtonTracksTOF.push_back(trk);
        } else if (!trk.hasTOF()) {
          onlyProtonTracks.push_back(proton);
          onlyProtonSigma.push_back(nSigmaPr);
          rawProtonTracks.push_back(trk);
        }
        registry.fill(HIST("hSelectionCounter"), 4);
        registry.fill(HIST("hClusterSizeProtonsTPC"), averageClusterSize);
      }

    } // trk loop

    registry.fill(HIST("hSelectionCounter"), 5);
    if ((onlyProtonTracksTOF.size() > 0) && (onlyProtonTracks.size() + onlyProtonTracksTOF.size()) == 2) {
      registry.fill(HIST("hSelectionCounter"), 6);

      int signSum = -999.;
      double sigmaTotal = -999.;
      TLorentzVector a, b;
      int positiveFlag = -1;
      // double dcaZ[2] = {-99., -99.};
      // double dcaXY[2] = {-99., -99.};
      // int dcaZbool = -1;
      // int dcaXYbool = -1;
      int dcaZbool = 1;
      int dcaXYbool = 1;
      if (onlyProtonTracksTOF.size() == 1) {

        if (!(fabs(onlyProtonTracks[0].Eta()) < 0.8 && fabs(onlyProtonTracksTOF[0].Eta()) < 0.8)) {
          return;
        }
        registry.fill(HIST("hEta1"), onlyProtonTracks[0].Eta());
        registry.fill(HIST("hEta1"), onlyProtonTracksTOF[0].Eta());
        resonance += onlyProtonTracks[0];
        resonance += onlyProtonTracksTOF[0];
        a = onlyProtonTracks[0];
        b = onlyProtonTracksTOF[0];
        sigmaTotal = 0;
        sigmaTotal = onlyProtonSigma[0] * onlyProtonSigma[0] + onlyProtonSigmaTOF[0] * onlyProtonSigmaTOF[0];
        ;
        registry.fill(HIST("hNsigEvsKa1"), onlyProtonSigma[0], onlyProtonSigmaTOF[0]);
        signSum = rawProtonTracks[0].sign() + rawProtonTracksTOF[0].sign();
        if (signSum == 0) {
          registry.fill(HIST("PP/hPtProtonVsProton"), onlyProtonTracks[0].Pt(), onlyProtonTracksTOF[0].Pt());
          if (rawProtonTracks[0].sign() == 1) {
            positiveFlag = 1;
          } else {
            positiveFlag = 2;
          }
        }

        // // DCA checks
        // dcaZ[0] = rawProtonTracks[0].dcaZ();
        // dcaZ[1] = rawProtonTracksTOF[0].dcaZ();
        // dcaXY[0] = rawProtonTracks[0].dcaXY();
        // dcaXY[1] = rawProtonTracksTOF[0].dcaXY();
        // if (std::abs(dcaZ[0]) < 2. && std::abs(dcaZ[1]) < 2.) {
        //   dcaZbool = 1;
        // } else {
        //   dcaZbool = 0;
        // }
        // double dcaLimit[2] = {-99., -99.};
        // dcaLimit[0] = 0.0105 + 0.035 / pow(a.Pt(), 1.1);
        // dcaLimit[1] = 0.0105 + 0.035 / pow(b.Pt(), 1.1);
        // if (std::abs(dcaXY[0]) < dcaLimit[0] && std::abs(dcaXY[1]) < dcaLimit[1]) {
        //   dcaXYbool = 1;
        // } else {
        //   dcaXYbool = 0;
        // }

      } else if (onlyProtonTracksTOF.size() == 2) {

        if (!(fabs(onlyProtonTracksTOF[0].Eta()) < 0.8 && fabs(onlyProtonTracksTOF[1].Eta()) < 0.8)) {
          return;
        }
        registry.fill(HIST("hEta1"), onlyProtonTracksTOF[0].Eta());
        registry.fill(HIST("hEta1"), onlyProtonTracksTOF[1].Eta());
        resonance += onlyProtonTracksTOF[0];
        resonance += onlyProtonTracksTOF[1];
        a = onlyProtonTracksTOF[0];
        b = onlyProtonTracksTOF[1];
        sigmaTotal = 0;
        sigmaTotal = onlyProtonSigmaTOF[0] * onlyProtonSigmaTOF[0] + onlyProtonSigmaTOF[1] * onlyProtonSigmaTOF[1];
        ;
        registry.fill(HIST("hNsigEvsKa1"), onlyProtonSigmaTOF[0], onlyProtonSigmaTOF[1]);
        signSum = rawProtonTracksTOF[0].sign() + rawProtonTracksTOF[1].sign();
        if (signSum == 0) {
          registry.fill(HIST("PP/hPtProtonVsProton"), onlyProtonTracksTOF[0].Pt(), onlyProtonTracksTOF[1].Pt());
          if (rawProtonTracksTOF[0].sign() == 1) {
            positiveFlag = 1;
          } else {
            positiveFlag = 2;
          }
        }

        // // DCA checks
        // dcaZ[0] = rawProtonTracksTOF[0].dcaZ();
        // dcaZ[1] = rawProtonTracksTOF[1].dcaZ();
        // dcaXY[0] = rawProtonTracksTOF[0].dcaXY();
        // dcaXY[1] = rawProtonTracksTOF[1].dcaXY();
        // if (std::abs(dcaZ[0]) < 2. && std::abs(dcaZ[1]) < 2.) {
        //   dcaZbool = 1;
        // } else {
        //   dcaZbool = 0;
        // }
        // double dcaLimit[2] = {-99., -99.};
        // dcaLimit[0] = 0.0105 + 0.035 / pow(a.Pt(), 1.1);
        // dcaLimit[1] = 0.0105 + 0.035 / pow(b.Pt(), 1.1);
        // if (std::abs(dcaXY[0]) < dcaLimit[0] && std::abs(dcaXY[1]) < dcaLimit[1]) {
        //   dcaXYbool = 1;
        // } else {
        //   dcaXYbool = 0;
        // }
      }

      if (sigmaTotal > 16.) {
        return;
      }
      if (onlyProtonTracksTOF.size() == 1) {
        registry.fill(HIST("hNsigEvsKa2"), onlyProtonSigma[0], onlyProtonSigmaTOF[0]);
      } else if (onlyProtonTracksTOF.size() == 2) {
        registry.fill(HIST("hNsigEvsKa2"), onlyProtonSigmaTOF[0], onlyProtonSigmaTOF[1]);
      }

      if (signSum != 0) {
        registry.fill(HIST("hMassPtLikeSignProton"), resonance.M(), resonance.Pt());
        registry.fill(HIST("hSelectionCounter"), 7);
        registry.fill(HIST("hPtLikeSignProton"), resonance.Pt());
        registry.fill(HIST("hMassLikeSignProton"), resonance.M());
      } else {
        registry.fill(HIST("PP/hMassPtUnlikeSignProton"), resonance.M(), resonance.Pt());
        registry.fill(HIST("hSelectionCounter"), 8);
        registry.fill(HIST("PP/hUnlikePt"), resonance.Pt());
        registry.fill(HIST("PP/hMassUnlike"), resonance.M());
        registry.fill(HIST("PP/hRapidity"), resonance.Rapidity());
        if (resonance.Pt() < 0.15) {
          registry.fill(HIST("PP/hCoherentMassWithoutDCAcuts"), resonance.M());
          if (dcaZbool == 1) {
            registry.fill(HIST("PP/hCoherentMassWithoutDCAxycut"), resonance.M());
          }
          if (dcaXYbool == 1) {
            registry.fill(HIST("PP/hCoherentMassWithoutDCAzcut"), resonance.M());
          }
        }
        if (resonance.Pt() < 0.15 && dcaZbool == 1 && dcaXYbool == 1) {
          registry.fill(HIST("PP/hCoherentMass"), resonance.M());
          if (resonance.M() < 3.0) {
            registry.fill(HIST("PP/hAngularDstribLab"), a.Phi() + TMath::Pi(), a.CosTheta());
            registry.fill(HIST("PP/hAngularDstribLab"), b.Phi() + TMath::Pi(), b.CosTheta());
            registry.fill(HIST("PP/hCosThetaLab"), a.CosTheta());
            registry.fill(HIST("PP/hCosThetaLab"), b.CosTheta());
          }
          double costhetaHel = -999.;
          double phiHel = -999.;
          if (positiveFlag == 1 && (resonance.M() > 2.4 && resonance.M() < 3.)) {
            costhetaHel = CosThetaHelicityFrame(a, b, resonance);
            phiHel = 1. * TMath::Pi() + PhiHelicityFrame(a, b, resonance);
            registry.fill(HIST("PP/hAngularDistribHelicity"), phiHel, costhetaHel);
            registry.fill(HIST("PP/hCosThetaHelicity"), costhetaHel);
            registry.fill(HIST("PP/hPhiHelicity"), phiHel);
          } else if (positiveFlag == 2 && (resonance.M() > 2.4 && resonance.M() < 3.)) {
            costhetaHel = CosThetaHelicityFrame(b, a, resonance);
            phiHel = 1. * TMath::Pi() + PhiHelicityFrame(b, a, resonance);
            registry.fill(HIST("PP/hAngularDistribHelicity"), phiHel, costhetaHel);
            registry.fill(HIST("PP/hCosThetaHelicity"), costhetaHel);
            registry.fill(HIST("PP/hPhiHelicity"), phiHel);
          }
        }
        // outside the hard pT cut, but with opposite charges
        if (dcaZbool == 1 && dcaXYbool == 1) {
          if (resonance.M() > 2.4 && resonance.M() < 2.75) {
            registry.fill(HIST("PP/hUnlikePtLowBand24to275"), resonance.Pt());
            registry.fill(HIST("PP/hUnlikePtLowBand"), resonance.Pt());
          } else if (resonance.M() > 2.75 && resonance.M() < 3.) {
            registry.fill(HIST("PP/hUnlikePtLowBand275to300"), resonance.Pt());
            registry.fill(HIST("PP/hUnlikePtLowBand"), resonance.Pt());
          } else if (resonance.M() > 3.0 && resonance.M() < 3.15) {
            registry.fill(HIST("PP/hUnlikePtJpsi"), resonance.Pt());
          } else if (resonance.M() > 3.15 && resonance.M() < 3.4) {
            registry.fill(HIST("PP/hUnlikePtHighBand"), resonance.Pt());
          }

          // high mass states
          if (resonance.M() > 4.) {
            if (onlyProtonTracksTOF.size() == 1) {
              registry.fill(HIST("hNsigEvsKa3"), onlyProtonSigma[0], onlyProtonSigmaTOF[0]);
            } else if (onlyProtonTracksTOF.size() == 2) {
              registry.fill(HIST("hNsigEvsKa3"), onlyProtonSigmaTOF[0], onlyProtonSigmaTOF[1]);
            }
            if (resonance.Pt() < 0.15) {
              registry.fill(HIST("PP/hMassHighPt0150"), resonance.M());
            }
            if (resonance.Pt() < 0.20) {
              registry.fill(HIST("PP/hMassHighPt0200"), resonance.M());
            }
            if (resonance.Pt() < 0.30) {
              registry.fill(HIST("PP/hMassHighPt0300"), resonance.M());
            }
            if (resonance.Pt() < 0.50) {
              registry.fill(HIST("PP/hMassHighPt0500"), resonance.M());
            }
            if (resonance.Pt() < 1.00) {
              registry.fill(HIST("PP/hMassHighPt1000"), resonance.M());
            }
          }
        }
      }
    }

  } // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusiveTwoProtonsSG>(cfgc)};
}
