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
/// \brief  Task to study two pions candidates using SG derive data
/// \author Levi Van Ryder (based on Anisa Khatun's UD Tutorial 5 example)
/// \file sgExclusiveJpsiMidrapidity.cxx

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TMath.h"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using LorentzVector = ROOT::Math::PxPyPzMVector;

// Struct to define the analysis task
struct SgExclusiveJpsiMidrapidity {
  // SGSelector object to manage track and collision selections
  SGSelector sgSelector;

  // Number of Selection cuts
  static constexpr int numSelectionCuts = 16;

  // Gapside rejection
  static constexpr int gapSideLow = 0;
  static constexpr int gapSideHigh = 2;

  // Numbers for selections
  static constexpr int two = 2;

  Configurable<float> fv0Cut{"fv0Cut", 100., "fv0aThreshold"};
  Configurable<float> ft0aCut{"ft0aCut", 100., "ft0aThreshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "ft0cThreshold"};
  Configurable<float> fddaCut{"fddaCut", 10000., "fddaThreshold"};
  Configurable<float> fddcCut{"fddcCut", 10000., "fddcThreshold"};
  Configurable<float> zdcCut{"zdcCut", 10., "zdcThreshold"};
  Configurable<float> selectedGapSide{"selectedGapSide", 2, "gapSelection"};

  // Track Selections
  Configurable<float> pvCut{"pvCut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcazCut{"dcazCut", 2.0, "dcaZ cut"};
  Configurable<float> dcaxyCut{"dcaxyCut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2Cut{"itsChi2Cut", 36, "Max itsChi2NCl"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> ptCut{"ptCut", 0.01, "Track Pt"};
  Configurable<float> tpcCluster{"tpcCluster", 50, "No.of TPC cluster"};

  // Kinematic cuts
  Configurable<float> pidCut{"pidCut", 5, "TPC PID"};
  Configurable<float> rapCut{"rapCut", 0.9, "Track rapidity"};
  Configurable<float> massMax{"massMax", 10, "Invariant Mass range high"};
  Configurable<float> massMin{"massMin", 0, "Invariant Mass range low"};
  Configurable<float> ptCoherent{"ptCoherent", 0.15, "Coherent selection"};

  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {

    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});

    // Fill counter to see effect of each selection criteria
    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{20, 0., 20.}});

    TString selectionCuts[numSelectionCuts] = {"NoSelection", "gapside", "goodtracks", "truegap", "2collcontrib", "2goodtrk", "TPCPID", "rapCut", "unlikesign", "mass_cut", "coherent", "incoherent", "likesign", "mass_cut", "coherent", "incoherent"};

    for (int i = 0; i < numSelectionCuts; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, selectionCuts[i].Data());
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
    registry.add("TwoPion/hTPCsig1", "TPC signal;signal_{TPC} t2; signal_{TPC} t2", kTH2F, {{300, 0, 150.}, {300, 0, 150}});

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
    registry.add("TwoPion/hPhiSystem", "Phi;#it{#Phi};", kTH1F, {{250, 0., o2::constants::math::TwoPI}});
    registry.add("TwoPion/hMPt", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 10.}, {100, 0., 10.}});
  }

  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;

  //__________________________________________________________________________
  // Main process
  void process(UDCollisionsFull::iterator const& collision, UDTracksFull const& tracks)
  {
    // No selection criteria
    registry.fill(HIST("hSelectionCounter"), 0);

    // Accessing gap sides
    int gapSide = collision.gapSide();
    if (gapSide < gapSideLow || gapSide > gapSideHigh)
      return;

    registry.fill(HIST("hSelectionCounter"), 1);

    // Accessing FIT information for further exclusivity and/or inclusivity
    int truegapSide = sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut);

    // Initiating track parameters to select good tracks, values to be optimized in the configurables, parameters will be taken from SGtrackselector.h task included in the header
    std::vector<float> parameters = {pvCut, dcazCut, dcaxyCut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, ptCut};

    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);

    // Gap side to be selected in the configurables
    gapSide = truegapSide;

    if (gapSide == selectedGapSide) {

      registry.fill(HIST("hSelectionCounter"), 2);

      if (collision.flags() != 1)
        return; // UPC setting vs std setting
      //____________________________________________________________________________________

      // Create LorentzVector to store all tracks, Pion tracks and TPC Pion PID
      std::vector<LorentzVector> onlyPionTracks;
      std::vector<float> onlyPionSigma;
      std::vector<decltype(tracks.begin())> rawPionTracks;

      // initialize pair 4-vector to zero before accumulation
      LorentzVector p(0.0, 0.0, 0.0, 0.0);

      registry.fill(HIST("hTracks"), tracks.size());

      for (const auto& t : tracks) {
        // Apply good track selection criteria
        if (!trackselector(t, parameters))
          continue;

        double dEdx = t.tpcSignal();

        registry.fill(HIST("hdEdx"), t.tpcInnerParam() / t.sign(), dEdx);

        // Create Lorentz vector for this track (use constructor, portable)
        LorentzVector a(t.px(), t.py(), t.pz(), o2::constants::physics::MassPionCharged);

        // Apply TPC pion sigma
        auto nSigmaPi = t.tpcNSigmaPi();
        if (std::fabs(nSigmaPi) < pidCut) {
          onlyPionTracks.push_back(a);
          onlyPionSigma.push_back(nSigmaPi);
          rawPionTracks.push_back(t);
          registry.fill(HIST("hdEdxPion"), t.tpcInnerParam() / t.sign(), dEdx);
        }
      }

      registry.fill(HIST("hTracksPions"), onlyPionTracks.size());

      //_____________________________________
      // Add all onlyPionTracks into p
      for (const auto& pion : onlyPionTracks) {
        p += pion;
      }

      //_____________________________________
      // Selecting collisions with Two PV contributors
      if (collision.numContrib() == two) {

        registry.fill(HIST("hSelectionCounter"), 3);

        // Selecting only Two good tracks
        if ((static_cast<int>(rawPionTracks.size()) == two) && (static_cast<int>(onlyPionTracks.size()) == two)) {

          registry.fill(HIST("hSelectionCounter"), 4);

          registry.fill(HIST("TwoPion/h2TracksPions"), onlyPionTracks.size());

          registry.fill(HIST("TwoPion/hNsigPivsPt1"), onlyPionTracks[0].Pt(), rawPionTracks[0].tpcNSigmaPi());
          registry.fill(HIST("TwoPion/hNsigPivsPt2"), onlyPionTracks[1].Pt(), rawPionTracks[1].tpcNSigmaPi());
          registry.fill(HIST("TwoPion/hTPCsig"), rawPionTracks[0].tpcSignal(), rawPionTracks[1].tpcSignal());
          registry.fill(HIST("TwoPion/hNsigPi1vsPi2"), rawPionTracks[0].tpcNSigmaPi(), rawPionTracks[1].tpcNSigmaPi());

          // Make sure two good tracks are within TPC pion sigma limit
          if ((onlyPionSigma[0] * onlyPionSigma[0] + onlyPionSigma[1] * onlyPionSigma[1]) > (pidCut * pidCut)) {
            return;
          }

          registry.fill(HIST("hSelectionCounter"), 5);

          // Rapidity of midrapidity acceptance
          if ((p.Rapidity() < -rapCut) || (p.Rapidity() > rapCut)) {
            return;
          }

          registry.fill(HIST("hSelectionCounter"), 6);

          //  opposite sign tracks
          if (rawPionTracks[0].sign() != rawPionTracks[1].sign()) {

            registry.fill(HIST("hSelectionCounter"), 7);
            registry.fill(HIST("TwoPion/hMassUnlike"), p.M());

            // Flexible mass limits, can be selected in the configurable
            if ((p.M() > massMin) && (p.M() < massMax)) {

              registry.fill(HIST("hSelectionCounter"), 8);

              registry.fill(HIST("TwoPion/hPt"), p.Pt());
              registry.fill(HIST("TwoPion/hEta"), p.Eta());
              registry.fill(HIST("TwoPion/hRap"), p.Rapidity());
              registry.fill(HIST("TwoPion/hPhiSystem"), p.Phi());
              registry.fill(HIST("TwoPion/hMPt"), p.M(), p.Pt());

              // flexible pt limit for selecting coherent Rho(0)
              if (p.Pt() < ptCoherent) {

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
              if (p.Pt() > ptCoherent) {
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
            if ((p.M() > massMin) && (p.M() < massMax)) {

              registry.fill(HIST("hSelectionCounter"), 12);
              registry.fill(HIST("TwoPion/hPtLike"), p.Pt());

              // Coherent Rho(0) selection
              if (p.Pt() < ptCoherent) {

                registry.fill(HIST("hSelectionCounter"), 13);
                registry.fill(HIST("TwoPion/Coherent/hMassLikeCoherent"), p.M());
              }
              // Incoherent Rho(0) selection
              if (p.Pt() > ptCoherent) {

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
    adaptAnalysisTask<SgExclusiveJpsiMidrapidity>(cfgc)};
}
