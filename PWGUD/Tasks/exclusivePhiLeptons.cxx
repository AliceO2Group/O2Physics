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
#include <vector>
using std::array;
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// \brief Exclusive phi->ee task
/// \author Simone Ragoni, Creighton
/// \date 5/5/2024

struct ExclusivePhiLeptons {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 200., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 100., "FT0C threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};
  Configurable<float> pid2d_cut{"PID2D", 2., "PID cut in 2D"};
  Configurable<float> pid_cut{"PID", 2., "PID cut in 1D"};
  Configurable<std::size_t> electronsInTOF{"eTOF", 2, "electrons in TOF"};
  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

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
    registry.add("hdSigmaElectron", "p vs dE/dx sigma electron", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hdSigmaElectron2", "p vs dE/dx sigma electron", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hdSigmaElectron3", "p vs dE/dx sigma electron", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hNsigEvsKa1", "NSigma(t1) vs NSigma (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigEvsKa2", "NSigma(t1) vs NSigma (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hMomentum", "p_{#ka};#it{p_{trk}}, GeV/c;", kTH1F, {{100, 0., 3.}});
    registry.add("hClusterSizeAllTracks", "ClusterSizeAllTracks;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeElectronsTPC", "ClusterSizeElectronsTPC;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeElectronsTOF", "ClusterSizeElectronsTOF;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hEta1", "#eta_{#ka};#it{#eta_{trk}}, GeV/c;", kTH1F, {{100, -2., 2.}});
    registry.add("hPtLikeSignElectron", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("hMassLikeSignElectron", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("hMassPtLikeSignElectron", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{1000, 0., 10.}, {400, 0., 4.}});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{10, 0., 10.}});

    TString SelectionCuts[9] = {"NoSelection", "GAPcondition", "PVtracks", "Good TPC-ITS track", "TPC/TOF PID track", "End trk loop", "Exactly 2p", "Like-sign ev", "Unlike-sign ev"};
    // now we can set BinLabel in histogram Registry

    for (int i = 0; i < 9; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    // Unlike sign pp
    registry.add("PP/hRapidity", "Rapidity;#it{y_{ee}};", kTH1F, {{100, -2., 2.}});
    registry.add("PP/hPtElectronVsElectron", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("PP/hMassPtUnlikeSignElectron", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PP/hMassUnlike", "m_{ee} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hUnlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePt2", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePt3", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePt4", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePt5", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hUnlikePt6", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hCoherentMass", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hCoherentMass2", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hCoherentMass3", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hCoherentMass4", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hCoherentMass5", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hCoherentMass6", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hIncoherentMass", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hIncoherentMass2", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hIncoherentMass3", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hIncoherentMass4", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hIncoherentMass5", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hIncoherentMass6", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hAngle", "Angular distrib helicity frame;#theta;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PP/hAngle2", "Angular distrib helicity frame;#theta;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PP/hAngle3", "Angular distrib helicity frame;#theta;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PP/hAngle4", "Angular distrib helicity frame;#theta;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PP/hAngle5", "Angular distrib helicity frame;#theta;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("PP/hAngle6", "Angular distrib helicity frame;#theta;", kTH1F, {{100, 0. * TMath::Pi(), 2. * TMath::Pi()}});
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
    std::vector<TLorentzVector> onlyElectronTracks;
    std::vector<TLorentzVector> onlyElectronTracksTOF;
    std::vector<float> onlyElectronSigma;
    std::vector<float> onlyElectronSigmaTOF;
    std::vector<decltype(tracks.begin())> rawElectronTracks;
    std::vector<decltype(tracks.begin())> rawElectronTracksTOF;

    int counterPV = 0;
    for (auto trk : tracks) {
      if (!trk.isPVContributor()) {
        continue;
      }
      counterPV += 1;
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
      if (trk.pt() < 0.300) {
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

      TLorentzVector electron;
      electron.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassElectron);
      auto nSigmaEl = trk.tpcNSigmaEl();
      auto nSigmaElTOF = trk.tofNSigmaEl();

      if (trk.hasTOF()) {
        registry.fill(HIST("hdSigmaElectron"), momentum, nSigmaElTOF);
      }
      if (fabs(nSigmaEl) < pid_cut) {
        registry.fill(HIST("hdEdx2"), momentum, dEdx);
        registry.fill(HIST("hdSigmaElectron2"), momentum, nSigmaEl);
        registry.fill(HIST("hMomentum"), momentum);
        if (trk.hasTOF() && fabs(nSigmaElTOF) < pid_cut) {
          registry.fill(HIST("hdSigmaElectron3"), momentum, nSigmaElTOF);
          registry.fill(HIST("hClusterSizeElectronsTOF"), averageClusterSize);
          onlyElectronTracksTOF.push_back(electron);
          onlyElectronSigmaTOF.push_back(nSigmaElTOF);
          rawElectronTracksTOF.push_back(trk);
        } else if (!trk.hasTOF()) {
          onlyElectronTracks.push_back(electron);
          onlyElectronSigma.push_back(nSigmaEl);
          rawElectronTracks.push_back(trk);
        }
        registry.fill(HIST("hSelectionCounter"), 4);
        registry.fill(HIST("hClusterSizeElectronsTPC"), averageClusterSize);
      }

    } // trk loop

    if (counterPV != 2) {
      return;
    }
    registry.fill(HIST("hSelectionCounter"), 5);
    if ((onlyElectronTracksTOF.size() >= electronsInTOF) && (onlyElectronTracks.size() + onlyElectronTracksTOF.size()) == 2) {
      registry.fill(HIST("hSelectionCounter"), 6);

      int signSum = -999.;
      double sigmaTotal = -999.;
      TVector3 momentum[2];
      int restrictedMomenta = -1;
      // two electrons in the TPC
      if (onlyElectronTracksTOF.size() == 0) {

        if (!(fabs(onlyElectronTracks[0].Eta()) < 0.8 && fabs(onlyElectronTracks[1].Eta()) < 0.8)) {
          return;
        }
        if (!(onlyElectronTracks[0].P() < 0.8 && onlyElectronTracks[0].P() > 0.3 && onlyElectronTracks[1].P() < 0.8 && onlyElectronTracks[1].P() > 0.3)) {
          return;
        }
        if (!(onlyElectronTracks[0].P() < 0.6 && onlyElectronTracks[0].P() > 0.4 && onlyElectronTracks[1].P() < 0.6 && onlyElectronTracks[1].P() > 0.4)) {
          restrictedMomenta = 1;
        }
        momentum[0] = onlyElectronTracks[0].Vect();
        momentum[1] = onlyElectronTracks[1].Vect();
        registry.fill(HIST("hEta1"), onlyElectronTracks[0].Eta());
        registry.fill(HIST("hEta1"), onlyElectronTracks[1].Eta());
        resonance += onlyElectronTracks[0];
        resonance += onlyElectronTracks[1];
        sigmaTotal = 0;
        sigmaTotal = onlyElectronSigma[0] * onlyElectronSigma[0] + onlyElectronSigma[1] * onlyElectronSigma[1];
        ;
        registry.fill(HIST("hNsigEvsKa1"), onlyElectronSigma[0], onlyElectronSigma[1]);
        signSum = rawElectronTracks[0].sign() + rawElectronTracks[1].sign();
        if (signSum == 0) {
          registry.fill(HIST("PP/hPtElectronVsElectron"), onlyElectronTracks[0].Pt(), onlyElectronTracks[1].Pt());
        }

      } else if (onlyElectronTracksTOF.size() == 1) {

        if (!(fabs(onlyElectronTracks[0].Eta()) < 0.8 && fabs(onlyElectronTracksTOF[0].Eta()) < 0.8)) {
          return;
        }
        if (!(onlyElectronTracks[0].P() < 0.8 && onlyElectronTracks[0].P() > 0.3 && onlyElectronTracksTOF[0].P() < 0.8 && onlyElectronTracksTOF[0].P() > 0.3)) {
          return;
        }
        if (!(onlyElectronTracks[0].P() < 0.6 && onlyElectronTracks[0].P() > 0.4 && onlyElectronTracksTOF[0].P() < 0.6 && onlyElectronTracksTOF[0].P() > 0.4)) {
          restrictedMomenta = 1;
        }
        momentum[0] = onlyElectronTracks[0].Vect();
        momentum[1] = onlyElectronTracksTOF[0].Vect();
        registry.fill(HIST("hEta1"), onlyElectronTracks[0].Eta());
        registry.fill(HIST("hEta1"), onlyElectronTracksTOF[0].Eta());
        resonance += onlyElectronTracks[0];
        resonance += onlyElectronTracksTOF[0];
        sigmaTotal = 0;
        sigmaTotal = onlyElectronSigma[0] * onlyElectronSigma[0] + onlyElectronSigmaTOF[0] * onlyElectronSigmaTOF[0];
        ;
        registry.fill(HIST("hNsigEvsKa1"), onlyElectronSigma[0], onlyElectronSigmaTOF[0]);
        signSum = rawElectronTracks[0].sign() + rawElectronTracksTOF[0].sign();
        if (signSum == 0) {
          registry.fill(HIST("PP/hPtElectronVsElectron"), onlyElectronTracks[0].Pt(), onlyElectronTracksTOF[0].Pt());
        }

      } else if (onlyElectronTracksTOF.size() == 2) {

        if (!(fabs(onlyElectronTracksTOF[0].Eta()) < 0.8 && fabs(onlyElectronTracksTOF[1].Eta()) < 0.8)) {
          return;
        }
        if (!(onlyElectronTracksTOF[0].P() < 0.8 && onlyElectronTracksTOF[0].P() > 0.3 && onlyElectronTracksTOF[1].P() < 0.8 && onlyElectronTracksTOF[1].P() > 0.3)) {
          return;
        }
        if (!(onlyElectronTracksTOF[0].P() < 0.6 && onlyElectronTracksTOF[0].P() > 0.4 && onlyElectronTracksTOF[1].P() < 0.6 && onlyElectronTracksTOF[1].P() > 0.4)) {
          restrictedMomenta = 1;
        }
        momentum[0] = onlyElectronTracksTOF[0].Vect();
        momentum[1] = onlyElectronTracksTOF[1].Vect();
        registry.fill(HIST("hEta1"), onlyElectronTracksTOF[0].Eta());
        registry.fill(HIST("hEta1"), onlyElectronTracksTOF[1].Eta());
        resonance += onlyElectronTracksTOF[0];
        resonance += onlyElectronTracksTOF[1];
        sigmaTotal = 0;
        sigmaTotal = onlyElectronSigmaTOF[0] * onlyElectronSigmaTOF[0] + onlyElectronSigmaTOF[1] * onlyElectronSigmaTOF[1];
        ;
        registry.fill(HIST("hNsigEvsKa1"), onlyElectronSigmaTOF[0], onlyElectronSigmaTOF[1]);
        signSum = rawElectronTracksTOF[0].sign() + rawElectronTracksTOF[1].sign();
        if (signSum == 0) {
          registry.fill(HIST("PP/hPtElectronVsElectron"), onlyElectronTracksTOF[0].Pt(), onlyElectronTracksTOF[1].Pt());
        }
      }

      if (sigmaTotal > pid2d_cut * pid2d_cut) {
        return;
      }
      if (onlyElectronTracksTOF.size() == 1) {
        registry.fill(HIST("hNsigEvsKa2"), onlyElectronSigma[0], onlyElectronSigmaTOF[0]);
      } else if (onlyElectronTracksTOF.size() == 2) {
        registry.fill(HIST("hNsigEvsKa2"), onlyElectronSigmaTOF[0], onlyElectronSigmaTOF[1]);
      }

      auto angleBetweenMomenta = momentum[0].Angle(momentum[1]);
      if (signSum == 0) {
        registry.fill(HIST("PP/hAngle"), angleBetweenMomenta);
      }

      // if (fabs(angleBetweenMomenta - TMath::Pi()) > TMath::Pi() / 18.) {
      //   return;
      // }

      if (signSum != 0) {
        registry.fill(HIST("hMassPtLikeSignElectron"), resonance.M(), resonance.Pt());
        registry.fill(HIST("hSelectionCounter"), 7);
        registry.fill(HIST("hPtLikeSignElectron"), resonance.Pt());
        registry.fill(HIST("hMassLikeSignElectron"), resonance.M());
      } else {
        registry.fill(HIST("PP/hMassPtUnlikeSignElectron"), resonance.M(), resonance.Pt());
        registry.fill(HIST("hSelectionCounter"), 8);
        if (fabs(angleBetweenMomenta - TMath::Pi()) < TMath::Pi() / 18.) {
          registry.fill(HIST("PP/hAngle2"), angleBetweenMomenta);
          if (resonance.M() > 1.02 && resonance.M() < 1.06) {
            registry.fill(HIST("PP/hUnlikePt"), resonance.Pt());
          }
          registry.fill(HIST("PP/hMassUnlike"), resonance.M());
          registry.fill(HIST("PP/hRapidity"), resonance.Rapidity());
          if (resonance.Pt() > 0.1) {
            registry.fill(HIST("PP/hIncoherentMass"), resonance.M());
          } else {
            registry.fill(HIST("PP/hCoherentMass"), resonance.M());
          }
        }
        if (fabs(angleBetweenMomenta - TMath::Pi()) < TMath::Pi() / 90.) {
          // two degs
          registry.fill(HIST("PP/hAngle3"), angleBetweenMomenta);
          if (resonance.M() > 1.02 && resonance.M() < 1.06) {
            registry.fill(HIST("PP/hUnlikePt3"), resonance.Pt());
          }
          if (resonance.Pt() > 0.1) {
            registry.fill(HIST("PP/hIncoherentMass3"), resonance.M());
          } else {
            registry.fill(HIST("PP/hCoherentMass3"), resonance.M());
          }
        }
        if (fabs(angleBetweenMomenta - TMath::Pi()) < TMath::Pi() / 180.) {
          // one deg
          registry.fill(HIST("PP/hAngle4"), angleBetweenMomenta);
          if (resonance.M() > 1.02 && resonance.M() < 1.06) {
            registry.fill(HIST("PP/hUnlikePt4"), resonance.Pt());
          }
          if (resonance.Pt() > 0.1) {
            registry.fill(HIST("PP/hIncoherentMass4"), resonance.M());
          } else {
            registry.fill(HIST("PP/hCoherentMass4"), resonance.M());
          }
        }
        if (fabs(angleBetweenMomenta - TMath::Pi()) < TMath::Pi() / 180.) {
          // one deg
          if (resonance.M() > 1.02 && resonance.M() < 1.06) {
            registry.fill(HIST("PP/hUnlikePt4"), resonance.Pt());
          }
          if (resonance.Pt() > 0.1) {
            registry.fill(HIST("PP/hIncoherentMass4"), resonance.M());
          } else {
            registry.fill(HIST("PP/hCoherentMass4"), resonance.M());
          }
        }
        if (fabs(angleBetweenMomenta - TMath::Pi()) < TMath::Pi() / 18. && restrictedMomenta == 1) {
          // ten degs
          registry.fill(HIST("PP/hAngle5"), angleBetweenMomenta);
          if (resonance.M() > 1.02 && resonance.M() < 1.06) {
            registry.fill(HIST("PP/hUnlikePt5"), resonance.Pt());
          }
          if (resonance.Pt() > 0.1) {
            registry.fill(HIST("PP/hIncoherentMass5"), resonance.M());
          } else {
            registry.fill(HIST("PP/hCoherentMass5"), resonance.M());
          }
        }
        if (fabs(angleBetweenMomenta - TMath::Pi()) < TMath::Pi() / 90. && restrictedMomenta == 1) {
          // two degs
          registry.fill(HIST("PP/hAngle6"), angleBetweenMomenta);
          if (resonance.M() > 1.02 && resonance.M() < 1.06) {
            registry.fill(HIST("PP/hUnlikePt6"), resonance.Pt());
          }
          if (resonance.Pt() > 0.1) {
            registry.fill(HIST("PP/hIncoherentMass6"), resonance.M());
          } else {
            registry.fill(HIST("PP/hCoherentMass6"), resonance.M());
          }
        }
        if (resonance.M() > 1.02 && resonance.M() < 1.06) {
          registry.fill(HIST("PP/hUnlikePt2"), resonance.Pt());
        }
      }
    }

  } // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusivePhiLeptons>(cfgc)};
}
