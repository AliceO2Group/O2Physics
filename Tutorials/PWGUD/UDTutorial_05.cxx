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

#include "iostream"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"

// using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// \brief  Example task to study two pions candidates using SG derive data
/// \author Anisa Khatun
/// \date 10.10.2024

// Struct to define the analysis task
struct UDTutorial05 {
  // SGSelector object to manage track and collision selections
  SGSelector sgSelector;

  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};

  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.01, "Track Pt"};
  Configurable<float> TPC_cluster{"TPC_cluster", 50, "No.of TPC cluster"};

  // Kinmatic cuts
  Configurable<float> PID_cut{"PID_cut", 5, "TPC PID"};
  Configurable<float> Rap_cut{"Rap_cut", 0.9, "Track rapidity"};
  Configurable<float> Mass_Max{"Mass_Max", 10, "Invariant Mass range high"};
  Configurable<float> Mass_Min{"Mass_Min", 0, "Invariant Mass range low"};
  Configurable<float> Pt_coherent{"Pt_coherent", 0.15, "Coherent selection"};

  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {

    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});

    // Fill counter to see effect of each selection criteria
    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{20, 0., 20.}});

    TString SelectionCuts[16] = {"NoSelection", "gapside", "goodtracks", "truegap", "2collcontrib", "2goodtrk", "TPCPID", "Rap_cut", "unlikesign", "mass_cut", "coherent", "incoherent", "likesign", "mass_cut", "coherent", "incoherent"};

    for (int i = 0; i < 16; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }
    // tracks
    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("TwoPion/h2TracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});

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

    registry.add("TwoPion/hPt", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("TwoPion/hPtLike", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{1000, 0., 10.}});
    registry.add("TwoPion/hEta", "Eta;#it{#eta};", kTH1F, {{500, -10., 10.}});
    registry.add("TwoPion/hRap", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
    registry.add("TwoPion/hPhiSystem", "Phi;#it{#Phi};", kTH1F, {{250, 0. * TMath::Pi(), 2. * TMath::Pi()}});
    registry.add("TwoPion/hMPt", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 10.}, {100, 0., 10.}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;

  //__________________________________________________________________________
  // Main process
  void process(UDCollisionsFull::iterator const& collision, udtracksfull const& tracks)
  {
    // No selection criteria
    registry.fill(HIST("hSelectionCounter"), 0);

    // Accessing gap sides
    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;

    registry.fill(HIST("hSelectionCounter"), 1);

    // Accessing FIT information for further exclusivity and/or inclusivity
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);

    // Intiating track parameters to select good tracks, values to be optimized in the configurables, parameters will be taken from SGtrackselector.h task included in the header
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};

    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);

    // Gap side to be selected in the configuables
    gapSide = truegapSide;

    if (gapSide == gap_Side) {

      registry.fill(HIST("hSelectionCounter"), 2);

      //____________________________________________________________________________________

      // Create LorentzVector to store all tracks, Pion tracks and TPC Pion PID
      std::vector<TLorentzVector> allTracks;
      std::vector<TLorentzVector> onlyPionTracks;
      std::vector<float> onlyPionSigma;
      std::vector<decltype(tracks.begin())> rawPionTracks;
      TLorentzVector p;

      registry.fill(HIST("hTracks"), tracks.size());

      for (auto t : tracks) {
        // Apply good track selection criteria
        if (!trackselector(t, parameters))
          continue;

        double dEdx = t.tpcSignal();

        registry.fill(HIST("hdEdx"), t.tpcInnerParam() / t.sign(), dEdx);

        // Creating Lorenz vector to store raw tracks and piontracks
        TLorentzVector a;
        a.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassPionCharged);
        allTracks.push_back(a);

        // Apply TPC pion sigma
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
      // Adding all onlypiontracks
      for (auto pion : onlyPionTracks) {
        p += pion;
      }

      //_____________________________________
      // Selecting collisions with Two PV contributors
      if (collision.numContrib() == 2) {

        registry.fill(HIST("hSelectionCounter"), 3);

        // Selecting only Two good tracks
        if ((rawPionTracks.size() == 2) && (onlyPionTracks.size() == 2)) {

          registry.fill(HIST("hSelectionCounter"), 4);

          registry.fill(HIST("TwoPion/h2TracksPions"), onlyPionTracks.size());

          registry.fill(HIST("TwoPion/hNsigPivsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaPi());
          registry.fill(HIST("TwoPion/hNsigPivsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaPi());
          registry.fill(HIST("TwoPion/hTPCsig"), rawPionTracks[0].tpcSignal(), rawPionTracks[1].tpcSignal());
          registry.fill(HIST("TwoPion/hNsigPi1vsPi2"), rawPionTracks[0].tpcNSigmaPi(), rawPionTracks[1].tpcNSigmaPi());

          // Make sure two good tracks are with TPC pion sigma limit
          if ((onlyPionSigma[0] * onlyPionSigma[0] + onlyPionSigma[1] * onlyPionSigma[1]) > (PID_cut * PID_cut)) {
            return;
          }

          registry.fill(HIST("hSelectionCounter"), 5);

          // Rapidity of midrapidity acceptance
          if ((p.Rapidity() < -Rap_cut) || (p.Rapidity() > Rap_cut)) {
            return;
          }

          registry.fill(HIST("hSelectionCounter"), 6);

          // Two oppsotite sign tracks
          if (rawPionTracks[0].sign() != rawPionTracks[1].sign()) {

            registry.fill(HIST("hSelectionCounter"), 7);
            registry.fill(HIST("TwoPion/hMassUnlike"), p.M());

            // Flexible mass limits, can be selected in the configurable
            if ((p.M() > Mass_Min) && (p.M() < Mass_Max)) {

              registry.fill(HIST("hSelectionCounter"), 8);

              registry.fill(HIST("TwoPion/hPt"), p.Pt());
              registry.fill(HIST("TwoPion/hEta"), p.Eta());
              registry.fill(HIST("TwoPion/hRap"), p.Rapidity());
              registry.fill(HIST("TwoPion/hPhiSystem"), p.Phi());
              registry.fill(HIST("TwoPion/hMPt"), p.M(), p.Pt());

              // flexible pt limit for selecting coherent Rho(0)
              if (p.Pt() < Pt_coherent) {

                registry.fill(HIST("hSelectionCounter"), 9);

                // Quality Control plots after coherent Rho(0) selection
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

                registry.fill(HIST("TwoPion/Coherent/hMassUnlikeCoherent"), p.M());
              }
              // Incoherent Rho(0) selection
              if (p.Pt() > Pt_coherent) {
                registry.fill(HIST("hSelectionCounter"), 10);
                registry.fill(HIST("TwoPion/Incoherent/hMassUnlikeInCoherent"), p.M());
              }
            }
          }

          // Same charge particles
          if (rawPionTracks[0].sign() == rawPionTracks[1].sign()) {

            registry.fill(HIST("hSelectionCounter"), 11);
            registry.fill(HIST("TwoPion/hMassLike"), p.M());

            // Mass limit
            if ((p.M() > Mass_Min) && (p.M() < Mass_Max)) {

              registry.fill(HIST("hSelectionCounter"), 12);
              registry.fill(HIST("TwoPion/hPtLike"), p.Pt());

              // Coherent Rho(0) selection
              if (p.Pt() < Pt_coherent) {

                registry.fill(HIST("hSelectionCounter"), 13);
                registry.fill(HIST("TwoPion/Coherent/hMassLikeCoherent"), p.M());
              }
              // Incoherent Rho(0) selection
              if (p.Pt() > Pt_coherent) {

                registry.fill(HIST("hSelectionCounter"), 14);
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
    adaptAnalysisTask<UDTutorial05>(cfgc, TaskName{"udtutorial05"})};
}
