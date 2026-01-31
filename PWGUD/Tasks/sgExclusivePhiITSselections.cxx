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

/// \brief Exclusive phi without PID
/// \author Simone Ragoni, Creighton
/// \date 8/8/2024

struct sgExclusivePhiITSselections {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 200., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 100., "FT0C threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};
  Configurable<bool> DGactive{"DGactive", false, "DG active"};
  Configurable<bool> SGactive{"SGactive", true, "SG active"};
  // Configurable<int> DGorSG{"DGorSG", 1, "SG = 1, DG = 2"};
  Configurable<int> NofITShits{"NofITShits", 7, "ITS layers hit"};
  Configurable<float> pt_threshold{"pt_threshold", 0.180, "pT threshold"};
  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {
    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{12, 0., 12.}});
    TString SelectionCuts[12] = {"NoSelection", "Trackloop", "PVtracks", "|nsigmaka|<3", "|nsigmapi|>3", "|nsigmael|>3", "|nsigmamu|>3", "two tracks", "Phi-peak", "pt<0.2 GeV/c", "pt>0.2 GeV/c"};

    for (int i = 0; i < 12; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});

    registry.add("posx", "Vertex position in x", kTH1F, {{100, -0.5, 0.5}});
    registry.add("posy", "Vertex position in y", kTH1F, {{100, -0.5, 0.5}});
    registry.add("posz", "Vertex position in z", kTH1F, {{1000, -100., 100.}});

    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksITSonly", "N_{tracks ITS only}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hITSCluster", "N_{cluster}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hChi2ITSTrkSegment", "N_{cluster}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTPCCluster", "N_{cluster}", kTH1F, {{200, -0.5, 199.5}});
    registry.add("hTracksKaons", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon1", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon2", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon3", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon4", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon5", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon6", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon7", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon8", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    registry.add("hdEdxKaon9", "p_{#ka} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});

    registry.add("hNsigEvsKa2", "NSigmaKa(t1) vs NSigmaKa (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, 0., 1000.}, {100, 0., 1000}});
    registry.add("hMomentum", "p_{#ka};#it{p_{trk}}, GeV/c;", kTH1F, {{100, 0., 3.}});
    registry.add("hClusterSizeAllTracks", "ClusterSizeAllTracks;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeMomentumCut", "ClusterSizeMomentumCut;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeOnlyIdentifiedKaons", "ClusterSizeOnlyIdentifiedKaons;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeOnlyITS", "ClusterSizeOnlyITS;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeOnlyITS2", "ClusterSizeOnlyITS;Average cls size in the ITS layers;", kTH1F, {{1000, 0., 100.}});
    registry.add("hClusterSizeAllTracksVsP", "ClusterSizeAllTracks vs p; p; Average cls size in the ITS layers", kTH2F, {{100, 0.0, 10.0}, {1000, 0., 100.}});
    registry.add("hClusterSizeMomentumCutVsP", "hClusterSizeMomentumCut vs p; p; Average cls size in the ITS layers", kTH2F, {{100, 0.0, 10.0}, {1000, 0., 100.}});
    registry.add("hClusterSizeOnlyITSVsP", "hClusterSizeOnlyITS vs p; p; Average cls size in the ITS layers", kTH2F, {{100, 0.0, 10.0}, {1000, 0., 100.}});
    registry.add("hEta1", "#eta_{#ka};#it{#eta_{trk}}, GeV/c;", kTH1F, {{100, -2., 2.}});

    auto hSelectionCounter2 = registry.add<TH1>("hSelectionCounter2", "hSelectionCounter;;NEvents", HistType::kTH1I, {{12, 0., 12.}});
    TString SelectionCuts2[12] = {"NoSelection", "Trackloop", "PVtracks", "|nTPCCluster|<50", " track pt<0.180 GeV/c", "Kaon Band", "ITSCluster<6", "two tracks", "Phi-peak", "pt<0.2 GeV/c", "pt>0.2 GeV/c"};

    for (int i = 0; i < 12; i++) {
      hSelectionCounter2->GetXaxis()->SetBinLabel(i + 1, SelectionCuts2[i].Data());
    }

    registry.add("hMassPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("hMassPtPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});

    // Phi peak region
    registry.add("PHI/hPtPhiWithoutPID", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hMassVsPt", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("PHI/hRapidityPhiWithoutPID", "Rapidity;#it{y_{KK}};", kTH1F, {{100, -2., 2.}});
    registry.add("PHI/hPtKaonVsKaon", "Pt1 vs Pt2;p_{T};p_{T};", kTH2F, {{100, 0., 3.}, {100, 0., 3.}});

    // DIfferent phi topologies, 2 identified kaons, 1 identified kaon + 1 ITS track with correct selections and N ITS clusters
    registry.add("KaonBandPHI/hMassPtPhiIdentifiedKaons", "Raw Inv.M;#it{m_{KK}} (GeV/c^{2});Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("KaonBandPHI/hMassPhiIdentifiedKaons", "Raw Inv.M;#it{M_{KK}} (GeV/c^{2});", kTH1F, {{400, 0., 4.}});
    registry.add("KaonBandPHI/hPtPhiIdentifiedKaons", "Pt;#it{p_{t}} (GeV/c);", kTH1F, {{400, 0., 4.}});
    registry.add("KaonBandPHI/hMassPtPhiIdentifiedKaonAndITSkaon", "Raw Inv.M;#it{m_{KK}} (GeV/c^{2});Pt;#it{p_{t}}, GeV/c;", kTH2F, {{400, 0., 4.}, {400, 0., 4.}});
    registry.add("KaonBandPHI/hMassPhiIdentifiedKaonAndITSkaon", "Raw Inv.M;#it{m_{KK}} (GeV/c^{2});", kTH1F, {{400, 0., 4.}});
    registry.add("KaonBandPHI/hPtPhiIdentifiedKaonAndITSkaon", "Pt;#it{p_{t}} (GeV/c);", kTH1F, {{400, 0., 4.}});

    registry.add("PHI/hMassLike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hMassUnlike", "m_{#pi#pi} [GeV/#it{c}^{2}]", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hUnlikePt", "Pt;#it{p_{t}}, GeV/c;", kTH1F, {{500, 0., 5.}});
    registry.add("PHI/hCoherentPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hInCoherentPhiWithoutPID", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hCoherentMassLike", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
    registry.add("PHI/hInCoherentMassLike", "Raw Inv.M;#it{m_{KK}}, GeV/c^{2};", kTH1F, {{400, 0., 4.}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  // Main process
  void eventprocessing(int gapSide, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);

    if (gapSide == gap_Side) {
      TLorentzVector phi, phiWithoutPID, phiWithKaonPID; // lorentz vectors of tracks and the mother

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
      std::vector<decltype(tracks.begin())> onlyKaonBandPID;
      std::vector<decltype(tracks.begin())> onlyITS;
      std::vector<TLorentzVector> allTracksAreKaons;
      std::vector<TLorentzVector> allTracksAreKaonsBandPID;
      std::vector<TLorentzVector> allTracksAreITSonlyAndFourITSclusters;
      for (auto t : tracks) {
        registry.fill(HIST("hSelectionCounter2"), 0);
        if (!t.isPVContributor()) {
          continue;
        }
        registry.fill(HIST("hSelectionCounter2"), 1);
        registry.fill(HIST("hTracks"), t.size());

        double momentum = TMath::Sqrt(t.px() * t.px() + t.py() * t.py() + t.pz() * t.pz());
        double dEdx = t.tpcSignal();

        int clusterSize[7];
        double averageClusterSize = 0.;
        double activeClusters = 0.;
        for (int i = 0; i < 7; i++) { // info stored in 4 bits
          clusterSize[i] = (((1 << 4) - 1) & (t.itsClusterSizes() >> 4 * i));
          auto clusterSizeValue = static_cast<double>(clusterSize[i]);
          if (clusterSizeValue != 0) {
            averageClusterSize += clusterSizeValue;
            activeClusters += 1;
          }
          averageClusterSize += static_cast<double>(clusterSize[i]);
        }
        if (activeClusters != 0) {
          averageClusterSize /= activeClusters;
        } else {
          averageClusterSize = -999.;
        }
        registry.fill(HIST("hClusterSizeAllTracks"), averageClusterSize);
        registry.fill(HIST("hClusterSizeAllTracksVsP"), momentum, averageClusterSize);

        int NFindable = t.tpcNClsFindable();
        int NMinusFound = t.tpcNClsFindableMinusFound();
        int NCluster = NFindable - NMinusFound;

        if (NCluster > 50) {
          continue;
        }
        registry.fill(HIST("hSelectionCounter2"), 3);
        registry.fill(HIST("hdEdxKaon5"), t.tpcInnerParam() / t.sign(), dEdx);

        if (t.pt() > pt_threshold) {
          continue;
        }
        if (!(std::abs(t.dcaZ()) < 2.)) {
          continue;
        }
        double dcaLimit = 0.0105 + 0.035 / pow(t.pt(), 1.1);
        if (!(std::abs(t.dcaXY()) < dcaLimit)) {
          continue;
        }

        registry.fill(HIST("hSelectionCounter2"), 4);
        registry.fill(HIST("hMomentum"), momentum);
        registry.fill(HIST("hdEdxKaon6"), t.tpcInnerParam() / t.sign(), dEdx);
        registry.fill(HIST("hdEdxKaon2"), momentum, dEdx);
        registry.fill(HIST("hClusterSizeMomentumCut"), averageClusterSize);
        registry.fill(HIST("hClusterSizeMomentumCutVsP"), momentum, averageClusterSize);

        onlyTwoTracks.push_back(t);

        TLorentzVector a;
        a.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassKaonCharged);
        allTracksAreKaons.push_back(a);

        bool kaonBand = false;
        if ((momentum < 0.150) && (dEdx > 300)) {
          kaonBand = true;
        } else if ((momentum > 0.150) && (momentum < 0.220) && (dEdx > 250)) {
          kaonBand = true;
        } else if ((momentum > 0.220) && (momentum < 0.300) && (dEdx > 180)) {
          kaonBand = true;
        } else if ((momentum > 0.300) && (momentum < 0.500) && (dEdx > 110)) {
          kaonBand = true;
        }

        if (kaonBand == true) {
          registry.fill(HIST("hSelectionCounter2"), 5);
          allTracksAreKaonsBandPID.push_back(a);
          onlyKaonBandPID.push_back(t);
          registry.fill(HIST("hdEdxKaon7"), t.tpcInnerParam() / t.sign(), dEdx);
          registry.fill(HIST("hClusterSizeOnlyIdentifiedKaons"), averageClusterSize);
        }

        if (NFindable < 1 && t.itsNCls() < NofITShits) {
          allTracksAreITSonlyAndFourITSclusters.push_back(a);
          onlyITS.push_back(t);
          registry.fill(HIST("hdEdxKaon8"), t.tpcInnerParam() / t.sign(), dEdx);
          registry.fill(HIST("hSelectionCounter2"), 6);
          registry.fill(HIST("hClusterSizeOnlyITS2"), averageClusterSize);
        }

      } // track loop

      //_____________________________________
      // Creating phis and saving all the information
      // in the case that there are ONLY 2 PV
      if (allTracksAreKaons.size() == 2) {
        registry.fill(HIST("hSelectionCounter2"), 7);
        for (auto kaon : allTracksAreKaons) {
          phiWithoutPID += kaon;
        }
        registry.fill(HIST("hTracksKaons"), allTracksAreKaons.size());
        registry.fill(HIST("hNsigEvsKa2"), onlyTwoTracks[0].tpcSignal(), onlyTwoTracks[1].tpcSignal());
        registry.fill(HIST("hEta1"), allTracksAreKaons[0].Eta());
        registry.fill(HIST("hEta1"), allTracksAreKaons[1].Eta());

        // All invariant mass region
        registry.fill(HIST("hMassPhiWithoutPID"), phiWithoutPID.M());
        registry.fill(HIST("hMassPtPhiWithoutPID"), phiWithoutPID.M(), phiWithoutPID.Pt());

        // Phi peak region
        registry.fill(HIST("PHI/hMassVsPt"), phiWithoutPID.M(), phiWithoutPID.Pt());
        if ((phiWithoutPID.M() > 0.98) && (phiWithoutPID.M() < 1.05)) {
          registry.fill(HIST("PHI/hPtPhiWithoutPID"), phiWithoutPID.Pt());
          registry.fill(HIST("PHI/hRapidityPhiWithoutPID"), phiWithoutPID.Rapidity());
          registry.fill(HIST("PHI/hPtKaonVsKaon"), allTracksAreKaons[0].Pt(), allTracksAreKaons[1].Pt());

          // unlike-sign
          if (onlyTwoTracks[0].sign() != onlyTwoTracks[1].sign()) {
            registry.fill(HIST("hSelectionCounter2"), 8);
            registry.fill(HIST("PHI/hUnlikePt"), phiWithoutPID.Pt());
            registry.fill(HIST("PHI/hMassUnlike"), phiWithoutPID.M());
            if (phiWithoutPID.Pt() < 0.1) {
              registry.fill(HIST("hSelectionCounter2"), 9);
              registry.fill(HIST("PHI/hCoherentPhiWithoutPID"), phiWithoutPID.M());
            } else {
              registry.fill(HIST("hSelectionCounter2"), 10);
              registry.fill(HIST("PHI/hInCoherentPhiWithoutPID"), phiWithoutPID.M());
            }
          } else {
            // Likesign quantities
            registry.fill(HIST("PHI/hMassLike"), phiWithoutPID.M());
            registry.fill(HIST("PHI/hlikePt"), phiWithoutPID.Pt());
            if (phiWithoutPID.Pt() < 0.2) {
              registry.fill(HIST("PHI/hCoherentMassLike"), phiWithoutPID.M());
            } else {
              registry.fill(HIST("PHI/hInCoherentMassLike"), phiWithoutPID.M());
            }
          }
        } // Mass cut
      } // end of two tracks only loop

      if (allTracksAreKaonsBandPID.size() == 2) {
        registry.fill(HIST("hTracksKaons"), allTracksAreKaonsBandPID.size() + 10);
        TLorentzVector reallyPhi;
        for (auto kaon : allTracksAreKaonsBandPID) {
          reallyPhi += kaon;
        }
        registry.fill(HIST("KaonBandPHI/hPtPhiIdentifiedKaons"), reallyPhi.Pt());
        registry.fill(HIST("KaonBandPHI/hMassPtPhiIdentifiedKaons"), reallyPhi.M(), reallyPhi.Pt());
        if (reallyPhi.Pt() < 0.2) {
          registry.fill(HIST("KanonBandPHI/hMassPhiIdentifiedKaons"), reallyPhi.M());
        }
      }

      if (allTracksAreKaonsBandPID.size() == 1 && allTracksAreITSonlyAndFourITSclusters.size() > 0) {

        registry.fill(HIST("hTracksKaons"), allTracksAreKaonsBandPID.size() + 20);
        registry.fill(HIST("hTracksKaons"), allTracksAreITSonlyAndFourITSclusters.size() + 40);

        double momentum = TMath::Sqrt(onlyKaonBandPID[0].px() * onlyKaonBandPID[0].px() + onlyKaonBandPID[0].py() * onlyKaonBandPID[0].py() + onlyKaonBandPID[0].pz() * onlyKaonBandPID[0].pz());
        double dEdx = onlyKaonBandPID[0].tpcSignal();
        registry.fill(HIST("hdEdxKaon9"), momentum, dEdx);
        registry.fill(HIST("hTracksITSonly"), allTracksAreITSonlyAndFourITSclusters.size());

        for (std::size_t kaon = 0; kaon < allTracksAreITSonlyAndFourITSclusters.size(); kaon++) {

          int clusterSize[7];
          double averageClusterSize = 0.;
          double activeClusters = 0.;
          for (int i = 0; i < 7; i++) { // info stored in 4 bits
            clusterSize[i] = (((1 << 4) - 1) & (onlyITS[kaon].itsClusterSizes() >> 4 * i));
            auto clusterSizeValue = static_cast<double>(clusterSize[i]);
            if (clusterSizeValue != 0) {
              averageClusterSize += clusterSizeValue;
              activeClusters += 1;
            }
            averageClusterSize += static_cast<double>(clusterSize[i]);
          }
          if (activeClusters != 0) {
            averageClusterSize /= activeClusters;
          } else {
            averageClusterSize = -999.;
          }
          registry.fill(HIST("hClusterSizeOnlyITS"), averageClusterSize);
          registry.fill(HIST("hClusterSizeOnlyITSVsP"), momentum, averageClusterSize);

          TLorentzVector reallyPhi;
          reallyPhi += allTracksAreKaonsBandPID[0];
          reallyPhi += allTracksAreITSonlyAndFourITSclusters[kaon];
          registry.fill(HIST("KaonBandPHI/hMassPtPhiIdentifiedKaonAndITSkaon"), reallyPhi.M(), reallyPhi.Pt());
          registry.fill(HIST("KaonBandPHI/hPtPhiIdentifiedKaonAndITSkaon"), reallyPhi.Pt());
          if (reallyPhi.Pt() < 0.2) {
            registry.fill(HIST("KaonBandPHI/hMassPhiIdentifiedKaonAndITSkaon"), reallyPhi.M());
          }
        }
      } // Kaon Band

    } // double gap
  } // end of process

  void processSG(UDCollisionsFull::iterator const& collision, udtracksfull const& tracks)
  // process function subscribing to SG data
  {
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

    eventprocessing(gapSide, tracks);
  }
  PROCESS_SWITCH(sgExclusivePhiITSselections, processSG, "Process SG data", SGactive);

  void processDG(aod::UDCollisions::iterator const& collision, udtracksfull const& tracks)
  // process function subscribing to DG data
  {
    int gapSide = 2;
    registry.fill(HIST("posx"), collision.posX());
    registry.fill(HIST("posy"), collision.posY());
    registry.fill(HIST("posz"), collision.posZ());
    registry.fill(HIST("GapSide"), gapSide);
    eventprocessing(gapSide, tracks);
  }
  PROCESS_SWITCH(sgExclusivePhiITSselections, processDG, "Process DG data", DGactive);

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<sgExclusivePhiITSselections>(cfgc)};
}
