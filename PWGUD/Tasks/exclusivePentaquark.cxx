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

/// \brief Exclusive proton+Jpsi task
/// \author Simone Ragoni, Creighton
/// \date 29/4/2024

struct ExclusivePentaquark {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 200., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 100., "FT0C threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};

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
    registry.add("hdEdx3", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdSigmaElectron", "p vs dE/dx sigma electron", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hdSigmaProton", "p vs dE/dx sigma proton", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hdSigmaProton2", "p vs dE/dx sigma proton", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hdSigmaProton3", "p vs dE/dx sigma proton", kTH2F, {{100, 0.0, 3.0}, {1000, -500.0, 500.0}});
    registry.add("hNsigEvsKa1", "NSigma(t1) vs NSigma (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
    registry.add("hNsigEvsKa2", "NSigma(t1) vs NSigma (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.}, {100, -15., 15}});
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
    registry.add("PP/hMassUnlike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 10.}});
    registry.add("PP/hUnlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PP/hCoherentMass", "Raw Inv.M;#it{m_{pp}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
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
    // if (gapSide != gap_Side) {
    //   return;
    // }

    TLorentzVector resonance; // lorentz vectors of tracks and the mother

    // ===================================
    // Task for p+Jpsi pairs with PID
    // ===================================
    std::vector<TLorentzVector> onlyLeptonTracks;
    std::vector<TLorentzVector> onlyProtonTracks;
    std::vector<TLorentzVector> onlyProtonTracksTOF;
    std::vector<float> onlyProtonSigma;
    std::vector<float> onlyProtonSigmaTOF;
    std::vector<decltype(tracks.begin())> rawLeptonTracks;
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
      if (trk.pt() < 0.1) {
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
      TLorentzVector muon;
      muon.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassMuon);
      TLorentzVector proton;
      proton.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassProton);
      TLorentzVector protonbar;
      protonbar.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassProtonBar);
      auto nSigmaPr = trk.tpcNSigmaPr();
      auto nSigmaPrTOF = trk.tofNSigmaPr();
      auto nSigmaElTOF = trk.tofNSigmaEl();
      auto nSigmaMuTOF = trk.tofNSigmaMu();

      if (trk.hasTOF()) {
        registry.fill(HIST("hdSigmaProton"), momentum, nSigmaPrTOF);
      }
      if (fabs(nSigmaPr) < 2.) {
        registry.fill(HIST("hdEdx2"), momentum, dEdx);
        registry.fill(HIST("hdSigmaProton2"), momentum, nSigmaPr);
        registry.fill(HIST("hMomentum"), momentum);
        if (trk.hasTOF() && fabs(nSigmaPrTOF) < 2.) {
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

      if (trk.hasTOF() && (fabs(nSigmaElTOF) < 2. || fabs(nSigmaMuTOF) < 2.) && trk.pt() > 0.5) {
        if (fabs(trk.tpcNSigmaEl()) < fabs(trk.tpcNSigmaMu())) {
          onlyLeptonTracks.push_back(electron);
        } else {
          onlyLeptonTracks.push_back(muon);
        }
        rawLeptonTracks.push_back(trk);
        registry.fill(HIST("hdEdx3"), momentum, dEdx);
        registry.fill(HIST("hdSigmaElectron"), momentum, nSigmaElTOF);
      }

    } // trk loop

    registry.fill(HIST("hSelectionCounter"), 5);

    if (onlyLeptonTracks.size() == 2) {
      int signSum = rawLeptonTracks[0].sign() + rawLeptonTracks[1].sign();
      TLorentzVector resonance;
      resonance += onlyLeptonTracks[0];
      resonance += onlyLeptonTracks[1];
      if (signSum == 0) {
        registry.fill(HIST("hMassPtLikeSignProton"), resonance.M(), resonance.Pt());
        registry.fill(HIST("hSelectionCounter"), 7);
        registry.fill(HIST("hPtLikeSignProton"), resonance.Pt());
        registry.fill(HIST("hMassLikeSignProton"), resonance.M());
      }

      TLorentzVector pentaquark;
      pentaquark += resonance;
      if (onlyProtonTracksTOF.size() == 1) {
        pentaquark += onlyProtonTracksTOF[0];
        if (resonance.M() > 3. && resonance.M() < 3.13) {
          registry.fill(HIST("PP/hMassPtUnlikeSignProton"), pentaquark.M(), pentaquark.Pt());
          registry.fill(HIST("hSelectionCounter"), 8);
          registry.fill(HIST("PP/hUnlikePt"), pentaquark.Pt());
          registry.fill(HIST("PP/hMassUnlike"), pentaquark.M());
          registry.fill(HIST("PP/hRapidity"), pentaquark.Rapidity());
          if (pentaquark.Pt() < 0.3) {
            registry.fill(HIST("PP/hCoherentMass"), pentaquark.M());
          }
        }
      }
    }

    // if ((onlyProtonTracksTOF.size() > 0) && (onlyProtonTracks.size() + onlyProtonTracksTOF.size()) == 2) {
    //   registry.fill(HIST("hSelectionCounter"), 6);

    //   int signSum = -999.;
    //   double sigmaTotal = -999.;
    //   if (onlyProtonTracksTOF.size() == 1) {

    //     if (!(fabs(onlyProtonTracks[0].Eta()) < 0.8 && fabs(onlyProtonTracksTOF[0].Eta()) < 0.8)) {
    //       return;
    //     }
    //     registry.fill(HIST("hEta1"), onlyProtonTracks[0].Eta());
    //     registry.fill(HIST("hEta1"), onlyProtonTracksTOF[0].Eta());
    //     resonance += onlyProtonTracks[0];
    //     resonance += onlyProtonTracksTOF[0];
    //     sigmaTotal = 0;
    //     sigmaTotal = onlyProtonSigma[0] * onlyProtonSigma[0] + onlyProtonSigmaTOF[0] * onlyProtonSigmaTOF[0];
    //     ;
    //     registry.fill(HIST("hNsigEvsKa1"), onlyProtonSigma[0], onlyProtonSigmaTOF[0]);
    //     signSum = rawProtonTracks[0].sign() + rawProtonTracksTOF[0].sign();
    //     if (signSum == 0) {
    //       registry.fill(HIST("PP/hPtProtonVsProton"), onlyProtonTracks[0].Pt(), onlyProtonTracksTOF[0].Pt());
    //     }

    //   } else if (onlyProtonTracksTOF.size() == 2) {

    //     if (!(fabs(onlyProtonTracksTOF[0].Eta()) < 0.8 && fabs(onlyProtonTracksTOF[1].Eta()) < 0.8)) {
    //       return;
    //     }
    //     registry.fill(HIST("hEta1"), onlyProtonTracksTOF[0].Eta());
    //     registry.fill(HIST("hEta1"), onlyProtonTracksTOF[1].Eta());
    //     resonance += onlyProtonTracksTOF[0];
    //     resonance += onlyProtonTracksTOF[1];
    //     sigmaTotal = 0;
    //     sigmaTotal = onlyProtonSigmaTOF[0] * onlyProtonSigmaTOF[0] + onlyProtonSigmaTOF[1] * onlyProtonSigmaTOF[1];
    //     ;
    //     registry.fill(HIST("hNsigEvsKa1"), onlyProtonSigmaTOF[0], onlyProtonSigmaTOF[1]);
    //     signSum = rawProtonTracksTOF[0].sign() + rawProtonTracksTOF[1].sign();
    //     if (signSum == 0) {
    //       registry.fill(HIST("PP/hPtProtonVsProton"), onlyProtonTracksTOF[0].Pt(), onlyProtonTracksTOF[1].Pt());
    //     }
    //   }

    //   if (sigmaTotal > 16.) {
    //     return;
    //   }
    //   if (onlyProtonTracksTOF.size() == 1) {
    //     registry.fill(HIST("hNsigEvsKa2"), onlyProtonSigma[0], onlyProtonSigmaTOF[0]);
    //   } else if (onlyProtonTracksTOF.size() == 2) {
    //     registry.fill(HIST("hNsigEvsKa2"), onlyProtonSigmaTOF[0], onlyProtonSigmaTOF[1]);
    //   }

    //   if (signSum != 0) {
    //     registry.fill(HIST("hMassPtLikeSignProton"), resonance.M(), resonance.Pt());
    //     registry.fill(HIST("hSelectionCounter"), 7);
    //     registry.fill(HIST("hPtLikeSignProton"), resonance.Pt());
    //     registry.fill(HIST("hMassLikeSignProton"), resonance.M());
    //   } else {
    //     registry.fill(HIST("PP/hMassPtUnlikeSignProton"), resonance.M(), resonance.Pt());
    //     registry.fill(HIST("hSelectionCounter"), 8);
    //     registry.fill(HIST("PP/hUnlikePt"), resonance.Pt());
    //     registry.fill(HIST("PP/hMassUnlike"), resonance.M());
    //     registry.fill(HIST("PP/hRapidity"), resonance.Rapidity());
    //     if (resonance.Pt() < 0.15) {
    //       registry.fill(HIST("PP/hCoherentMass"), resonance.M());
    //     }
    //   }
    // }

  } // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusivePentaquark>(cfgc)};
}
