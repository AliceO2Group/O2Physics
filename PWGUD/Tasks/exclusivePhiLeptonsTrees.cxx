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

/// \brief Exclusive phi->ee tree producer, for ML applications, DG-based
/// \author Simone Ragoni, Creighton
/// \date 11/11/2024

namespace o2::aod
{
namespace tree
{
// track tables
DECLARE_SOA_COLUMN(PX1, px1, float);
DECLARE_SOA_COLUMN(PY1, py1, float);
DECLARE_SOA_COLUMN(PZ1, pz1, float);
DECLARE_SOA_COLUMN(PE1, pE1, float);
DECLARE_SOA_COLUMN(PX2, px2, float);
DECLARE_SOA_COLUMN(PY2, py2, float);
DECLARE_SOA_COLUMN(PZ2, pz2, float);
DECLARE_SOA_COLUMN(PE2, pE2, float);
DECLARE_SOA_COLUMN(NCOUNTERPV, nCounterPV, int);
DECLARE_SOA_COLUMN(NELECTRONSTOF, nElectronsTOF, int);
} // namespace tree

DECLARE_SOA_TABLE(TREE, "AOD", "Tree", tree::PX1, tree::PY1, tree::PZ1, tree::PE1, tree::PX2, tree::PY2, tree::PZ2, tree::PE2, tree::NCOUNTERPV, tree::NELECTRONSTOF);
} // namespace o2::aod

struct ExclusivePhiLeptonsTrees {
  Produces<aod::TREE> tree;
  Configurable<float> gap_Side{"gap", 2, "gap selection"};
  Configurable<float> pid2d_cut{"PID2D", 2., "PID cut in 2D"};
  Configurable<float> pid_cut{"PID", 2., "PID cut in 1D"};
  Configurable<std::size_t> electronsInTOF{"eTOF", 2, "electrons in TOF"};
  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {
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
    registry.add("hEta1", "#eta_{#ka};#it{#eta_{trk}}, GeV/c;", kTH1F, {{100, -2., 2.}});
    registry.add("hPtLikeSignElectron", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("hMassLikeSignElectron", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("hMassPtLikeSignElectron", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{1000, 0., 10.}, {400, 0., 4.}});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{10, 0., 10.}});

    TString SelectionCuts[9] = {"NoSelection", "GAPcondition", "PVtracks", "Good TPC-ITS track", "TPC/TOF PID track", "End trk loop", "Exactly 2e", "Like-sign ev", "Unlike-sign ev"};
    // now we can set BinLabel in histogram Registry

    for (int i = 0; i < 9; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    // Unlike sign pp
    registry.add("ee/hRapidity", "Rapidity;#it{y_{ee}};", kTH1F, {{100, -2., 2.}});
    registry.add("ee/hPtElectronVsElectron", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});
    registry.add("ee/hMassPtUnlikeSignElectron", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("ee/hMassUnlike", "m_{ee} [GeV/#it{c}^{2}]", kTH1F, {{1000, 0., 10.}});
    registry.add("ee/hUnlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("ee/hCoherentMass", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
    registry.add("ee/hIncoherentMass", "Raw Inv.M;#it{m_{ee}}, GeV/c^{2};", kTH1F, {{1000, 0., 10.}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  // using UDCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  // Main process
  void process(UDCollisions::iterator const& collision, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);
    registry.fill(HIST("posx"), collision.posX());
    registry.fill(HIST("posy"), collision.posY());
    registry.fill(HIST("posz"), collision.posZ());
    TLorentzVector resonance; // lorentz vectors of tracks and the mother
    // ===================================
    // Task for ee pairs with PID
    // Topology:
    // - 2 TOF ee
    // - 1 TOF e + 1 TPC e
    // ===================================
    std::vector<TLorentzVector> onlyElectronTracks;
    std::vector<TLorentzVector> onlyElectronTracksTOF;
    std::vector<float> onlyElectronSigma;
    std::vector<float> onlyElectronSigmaTOF;
    std::vector<decltype(tracks.begin())> rawElectronTracks;
    std::vector<decltype(tracks.begin())> rawElectronTracksTOF;

    // -------------------------------------------
    // TO BE SAVED:
    // - counterPV
    // - electronsTOF (0,1,2) = 2 - electronsTPC
    // - (px,py,pz,E)1
    // - (px,py,pz,E)2
    int counterPV = 0;
    for (auto trk : tracks) {
      // ----------------------------------------
      // SELECTIONS:
      // - PV track
      // - at least 70 TPC clusters
      // - at least 6 ITS clusters
      // - 0.3 < pT < 0.65 GeV from STARlight
      // - DCAxy, DCAz
      // - Nsigma^2 < 2^2
      // - |track eta| < 0.8
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
      if (trk.itsNCls() < 6) {
        continue;
      }
      if (trk.pt() < 0.300) {
        continue;
      }
      if (trk.pt() > 0.650) {
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

      double momentum = TMath::Sqrt(trk.px() * trk.px() + trk.py() * trk.py() + trk.pz() * trk.pz());
      double dEdx = trk.tpcSignal();
      registry.fill(HIST("hdEdx"), momentum, dEdx);

      TLorentzVector electron;
      electron.SetXYZM(trk.px(), trk.py(), trk.pz(), o2::constants::physics::MassElectron);
      if (fabs(electron.Eta()) > 0.8) {
        return;
      }
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
          onlyElectronTracksTOF.push_back(electron);
          onlyElectronSigmaTOF.push_back(nSigmaElTOF);
          rawElectronTracksTOF.push_back(trk);
        } else if (!trk.hasTOF()) {
          onlyElectronTracks.push_back(electron);
          onlyElectronSigma.push_back(nSigmaEl);
          rawElectronTracks.push_back(trk);
        }
        registry.fill(HIST("hSelectionCounter"), 4);
      }

    } // trk loop

    registry.fill(HIST("hSelectionCounter"), 5);
    if ((onlyElectronTracksTOF.size() >= electronsInTOF) && (onlyElectronTracks.size() + onlyElectronTracksTOF.size()) == 2) {
      registry.fill(HIST("hSelectionCounter"), 6);

      int signSum = -999.;
      double sigmaTotal = -999.;
      TLorentzVector a, b;
      // two electrons in the TPC
      if (onlyElectronTracksTOF.size() == 0) {

        registry.fill(HIST("hEta1"), onlyElectronTracks[0].Eta());
        registry.fill(HIST("hEta1"), onlyElectronTracks[1].Eta());
        resonance += onlyElectronTracks[0];
        resonance += onlyElectronTracks[1];
        a += onlyElectronTracks[0];
        b += onlyElectronTracks[1];
        sigmaTotal = 0;
        sigmaTotal = onlyElectronSigma[0] * onlyElectronSigma[0] + onlyElectronSigma[1] * onlyElectronSigma[1];
        ;
        registry.fill(HIST("hNsigEvsKa1"), onlyElectronSigma[0], onlyElectronSigma[1]);
        signSum = rawElectronTracks[0].sign() + rawElectronTracks[1].sign();
        if (signSum == 0) {
          registry.fill(HIST("ee/hPtElectronVsElectron"), onlyElectronTracks[0].Pt(), onlyElectronTracks[1].Pt());
        }

      } else if (onlyElectronTracksTOF.size() == 1) {

        registry.fill(HIST("hEta1"), onlyElectronTracks[0].Eta());
        registry.fill(HIST("hEta1"), onlyElectronTracksTOF[0].Eta());
        resonance += onlyElectronTracks[0];
        resonance += onlyElectronTracksTOF[0];
        a += onlyElectronTracks[0];
        b += onlyElectronTracksTOF[0];
        sigmaTotal = 0;
        sigmaTotal = onlyElectronSigma[0] * onlyElectronSigma[0] + onlyElectronSigmaTOF[0] * onlyElectronSigmaTOF[0];
        ;
        registry.fill(HIST("hNsigEvsKa1"), onlyElectronSigma[0], onlyElectronSigmaTOF[0]);
        signSum = rawElectronTracks[0].sign() + rawElectronTracksTOF[0].sign();
        if (signSum == 0) {
          registry.fill(HIST("ee/hPtElectronVsElectron"), onlyElectronTracks[0].Pt(), onlyElectronTracksTOF[0].Pt());
        }

      } else if (onlyElectronTracksTOF.size() == 2) {

        registry.fill(HIST("hEta1"), onlyElectronTracksTOF[0].Eta());
        registry.fill(HIST("hEta1"), onlyElectronTracksTOF[1].Eta());
        resonance += onlyElectronTracksTOF[0];
        resonance += onlyElectronTracksTOF[1];
        a += onlyElectronTracksTOF[0];
        b += onlyElectronTracksTOF[1];
        sigmaTotal = 0;
        sigmaTotal = onlyElectronSigmaTOF[0] * onlyElectronSigmaTOF[0] + onlyElectronSigmaTOF[1] * onlyElectronSigmaTOF[1];
        ;
        registry.fill(HIST("hNsigEvsKa1"), onlyElectronSigmaTOF[0], onlyElectronSigmaTOF[1]);
        signSum = rawElectronTracksTOF[0].sign() + rawElectronTracksTOF[1].sign();
        if (signSum == 0) {
          registry.fill(HIST("ee/hPtElectronVsElectron"), onlyElectronTracksTOF[0].Pt(), onlyElectronTracksTOF[1].Pt());
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

      if (signSum != 0) {
        registry.fill(HIST("hMassPtLikeSignElectron"), resonance.M(), resonance.Pt());
        registry.fill(HIST("hSelectionCounter"), 7);
        registry.fill(HIST("hPtLikeSignElectron"), resonance.Pt());
        registry.fill(HIST("hMassLikeSignElectron"), resonance.M());
      } else {
        registry.fill(HIST("ee/hMassPtUnlikeSignElectron"), resonance.M(), resonance.Pt());
        registry.fill(HIST("hSelectionCounter"), 8);
        registry.fill(HIST("ee/hMassUnlike"), resonance.M());
        registry.fill(HIST("ee/hRapidity"), resonance.Rapidity());
        if (resonance.Pt() > 0.1) {
          registry.fill(HIST("ee/hIncoherentMass"), resonance.M());
        } else {
          registry.fill(HIST("ee/hCoherentMass"), resonance.M());
        }
        if (resonance.M() > 1.01 && resonance.M() < 1.03) {
          registry.fill(HIST("ee/hUnlikePt"), resonance.Pt());
        }
      }
      // Filling tree, make to be consistent with the declared tables
      tree(a.Px(), a.Py(), a.Pz(), a.E(), b.Px(), b.Py(), b.Pz(), b.E(), counterPV, onlyElectronTracksTOF.size());
    }
  } // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusivePhiLeptonsTrees>(cfgc)};
}
