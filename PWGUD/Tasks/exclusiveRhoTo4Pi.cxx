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

/// \file exclusiveRhoTo4Pi.cxx
/// \brief Task for analyzing exclusive rho decays to 4 pions
/// \author Anantha Padmanabhan M Nair

#include <cstdlib>
#include <vector>
#include <cmath>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "Common/DataModel/PIDResponse.h"
#include <TString.h>
#include "TLorentzVector.h"
#include <TMath.h>
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/Boost.h"

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct exclusiveRhoTo4Pi { // o2-linter: disable=name/workflow-file,name/struct
  SGSelector sgSelector;

  HistogramRegistry histosData{"histosData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosMCgen{"histosMCgen", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosMCreco{"histosMCreco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Configurable<float> fv0Cut{"fv0Cut", 50., "FV0A threshold"};
  Configurable<float> ft0aCut{"ft0aCut", 150., "FT0A threshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "FT0C threshold"};
  Configurable<float> fddaCut{"fddaCut", 10000., "FDDA threshold"};
  Configurable<float> fddcCut{"fddcCut", 10000., "FDDC threshold"};
  Configurable<float> zdcCut{"zdcCut", 10., "ZDC threshold"};

  Configurable<float> pvCut{"pvCut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZcut{"dcaZcut", 3.2, "dcaZ cut"};
  Configurable<float> dcaXYcut{"dcaXYcut", 2.4, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2Cut{"itsChi2Cut", 36, "Max itsChi2NCl"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pTcut{"pTcut", 0.15, "Track Pt"};

  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 3, "TPC cut"};
  Configurable<float> nSigmaTOFcut{"nSigmaTOFcut", 3, "TOF cut"};
  Configurable<bool> strictEventSelection{"strictEventSelection", true, "Event Selection"};
  Configurable<bool> ifDataAnalysis{"ifDataAnalysis", false, "Data Analysis"};
  Configurable<bool> ifMCAnalysis{"ifMCAnalysis", true, "MC Analysis"};

  Configurable<int> nBinsPt{"nBinsPt", 1000, "Number of bins for pT"};
  Configurable<int> nBinsInvariantMass{"nBinsInvariantMass", 1000, "Number of bins for Invariant Mass"};
  Configurable<int> nBinsRapidity{"nBinsRapidity", 1000, "Number of bins for Rapidity"};
  Configurable<int> nBinsPhi{"nBinsPhi", 360, "Number of bins for Phi"};
  Configurable<int> nBinsCosTheta{"nBinsCosTheta", 360, "Number of bins for cos Theta"};

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Begin of Init Function-----------------------------------------------------------------------------------------------------------------------------------------------------
  void init(InitContext const&)
  {

    if (ifDataAnalysis) {
      histosData.add("GapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
      histosData.add("TrueGapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
      histosData.add("EventCounts", "Total Events; Events", kTH1F, {{10, 0, 10}});

      // TPC nSigma
      histosData.add("tpcNSigmaPi_WOTS", "TPC nSigma Pion without track selection; Events", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tpcNSigmaPi_WTS", "TPC nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tpcNSigmaPi_WTS_PID_Pi", "TPC nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});

      // TPC nSigma of other particles with selected pion tracks
      histosData.add("tpcNSigmaKa_WTS_PID_Pi", "TPC nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tpcNSigmaPr_WTS_PID_Pi", "TPC nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tpcNSigmaEl_WTS_PID_Pi", "TPC nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tpcNSigmaMu_WTS_PID_Pi", "TPC nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});

      // TOF nSigma
      histosData.add("tofNSigmaPi_WTS", "TOF nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tofNSigmaPi_WOTS", "TOF nSigma Pion without track selection; Events", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tofNSigmaPi_WTS_PID_Pi", "TOF nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});

      // TOF nSigma of other particles with selected pion tracks
      histosData.add("tofNSigmaKa_WTS_PID_Pi", "TOF nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tofNSigmaPr_WTS_PID_Pi", "TOF nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tofNSigmaEl_WTS_PID_Pi", "TOF nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosData.add("tofNSigmaMu_WTS_PID_Pi", "TOF nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});

      // Track Transverse Momentum
      histosData.add("pT_track_WOTS", "pT without track selection; pT [GeV/c]; Counts", kTH1F, {{nBinsPt, 0, 2}});
      histosData.add("pT_track_WTS", "pT with track selection; pT [GeV/c]; Counts", kTH1F, {{nBinsPt, 0, 2}});
      histosData.add("pT_track_WTS_PID_Pi", "pT with track selection and PID selection of Pi; pT [GeV/c]; Events", kTH1F, {{nBinsPt, 0, 2}});

      // Zero charge Event Transverse Momentum
      histosData.add("pT_event_0charge_WTS_PID_Pi", "Event pT in 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {{nBinsPt, 0, 2}});

      // Non Zero charge Event Transverse Momentum
      histosData.add("pT_event_non0charge_WTS_PID_Pi", "Event pT in Non 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {{nBinsPt, 0, 2}});

      // Rapidity of 0 charge Events
      histosData.add("rapidity_event_0charge_WTS_PID_Pi_domainA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});
      histosData.add("rapidity_event_0charge_WTS_PID_Pi_domainB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});
      histosData.add("rapidity_event_0charge_WTS_PID_Pi_domainC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});

      // Rapidity of non 0 charge Events
      histosData.add("rapidity_event_non0charge_WTS_PID_Pi_domainA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Events", kTH1F, {{nBinsRapidity, -2.5, 2.5}});
      histosData.add("rapidity_event_non0charge_WTS_PID_Pi_domainB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c$; y; Events", kTH1F, {{nBinsRapidity, -2.5, 2.5}});
      histosData.add("rapidity_event_non0charge_WTS_PID_Pi_domainC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Events", kTH1F, {{nBinsRapidity, -2.5, 2.5}});

      // Invariant Mass of 0 charge events
      histosData.add("invMass_event_0charge_WTS_PID_Pi_domainA", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});       // pT < 0.15GeV
      histosData.add("invMass_event_0charge_WTS_PID_Pi_domainB", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}}); // 0.15GeV < pT < 0.8GeV
      histosData.add("invMass_event_0charge_WTS_PID_Pi_domainC", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});       // 0.8GeV < pT

      // Invariant mass of non 0 charge events
      histosData.add("invMass_event_non0charge_WTS_PID_Pi_domainA", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});       // pT < 0.15GeV
      histosData.add("invMass_event_non0charge_WTS_PID_Pi_domainB", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}}); // 0.15GeV < pT < 0.8GeV
      histosData.add("invMass_event_non0charge_WTS_PID_Pi_domainC", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});       // 0.8GeV < pT

      // tpc signal
      histosData.add("tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
      histosData.add("tpcSignal_Pi", "TPC dEdx vs p for pions; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});

      // tof beta
      histosData.add("tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
      histosData.add("tofBeta_Pi", "TOF beta vs p for pions; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});

      // Other signals
      histosData.add("FT0A", "T0A amplitude", kTH1F, {{200, 0.0, 500.0}});
      histosData.add("FT0C", "T0C amplitude", kTH1F, {{200, 0.0, 500.0}});
      histosData.add("ZDC_A", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
      histosData.add("ZDC_C", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
      histosData.add("V0A", "V0A amplitude", kTH1F, {{1000, 0.0, 100}});

      // Collin Soper Theta and Phi
      histosData.add("CS_phi_pair_1", "#phi Distribution; #phi; Events", kTH1F, {{nBinsPhi, -3.2, 3.2}});
      histosData.add("CS_phi_pair_2", "#phi Distribution; #phi; Events", kTH1F, {{nBinsPhi, -3.2, 3.2}});
      histosData.add("CS_costheta_pair_1", "#theta Distribution;cos(#theta); Counts", kTH1F, {{nBinsCosTheta, -1, 1}});
      histosData.add("CS_costheta_pair_2", "#theta Distribution;cos(#theta); Counts", kTH1F, {{nBinsCosTheta, -1, 1}});
    }

    // MC Gen Stuff
    if (ifMCAnalysis) {
      // counts
      histosMCgen.add("rhoPrimeCounts", "Total Rho prime Events; Events", kTH1F, {{10, 0, 10}});

      // Track Stuff
      histosMCgen.add("MCgen_particle_pT", "Generated pT; pT [GeV/c]; Events", kTH1F, {{nBinsPt, 0, 10}});
      histosMCgen.add("MCgen_particle_rapidity", "Generated Rapidity; y; Events", kTH1F, {{nBinsRapidity, -2.5, 2.5}});

      // Generated Transverse Momentum, Rapidty and Invariant Mass
      histosMCgen.add("MCgen_4pion_pT", "Generated pT; pT [GeV/c]; Events", kTH1F, {{nBinsPt, 0, 2}});
      histosMCgen.add("MCgen_4pion_rapidity", "Generated Rapidity; y; Events", kTH1F, {{nBinsRapidity, -2.5, 2.5}});
      histosMCgen.add("MCgen_4pion_invmass", "Invariant Mass of 4-Pions; m(4-pion); Events", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});

      // Collin Soper Theta and Phi
      histosMCgen.add("MCgen_CS_phi_pair_1", "#phi Distribution; #phi; Events", kTH1F, {{nBinsPhi, -3.2, 3.2}});
      histosMCgen.add("MCgen_CS_phi_pair_2", "#phi Distribution; #phi; Events", kTH1F, {{nBinsPhi, -3.2, 3.2}});
      histosMCgen.add("MCgen_CS_costheta_pair_1", "#theta Distribution;cos(#theta); Events", kTH1F, {{nBinsCosTheta, -1, 1}});
      histosMCgen.add("MCgen_CS_costheta_pair_2", "#theta Distribution;cos(#theta); Events", kTH1F, {{nBinsCosTheta, -1, 1}});

      // MC Reco Stuff

      histosMCreco.add("GapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
      histosMCreco.add("TrueGapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
      histosMCreco.add("EventCounts", "Total Events; Events", kTH1F, {{10, 0, 10}});

      // TPC nSigma
      histosMCreco.add("tpcNSigmaPi_WOTS", "TPC nSigma Pion without track selection; Events", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tpcNSigmaPi_WTS", "TPC nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tpcNSigmaPi_WTS_PID_Pi", "TPC nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});

      // TPC nSigma of other particles with selected pion tracks
      histosMCreco.add("tpcNSigmaKa_WTS_PID_Pi", "TPC nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tpcNSigmaPr_WTS_PID_Pi", "TPC nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tpcNSigmaEl_WTS_PID_Pi", "TPC nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tpcNSigmaMu_WTS_PID_Pi", "TPC nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});

      // TOF nSigma
      histosMCreco.add("tofNSigmaPi_WTS", "TOF nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tofNSigmaPi_WOTS", "TOF nSigma Pion without track selection; Events", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tofNSigmaPi_WTS_PID_Pi", "TOF nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});

      // TOF nSigma of other particles with selected pion tracks
      histosMCreco.add("tofNSigmaKa_WTS_PID_Pi", "TOF nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tofNSigmaPr_WTS_PID_Pi", "TOF nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tofNSigmaEl_WTS_PID_Pi", "TOF nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});
      histosMCreco.add("tofNSigmaMu_WTS_PID_Pi", "TOF nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {nBinsPt, 0, 10}});

      // Track Transverse Momentum
      histosMCreco.add("pT_track_WOTS", "pT without track selection; pT [GeV/c]; Counts", kTH1F, {{nBinsPt, 0, 2}});
      histosMCreco.add("pT_track_WTS", "pT with track selection; pT [GeV/c]; Counts", kTH1F, {{nBinsPt, 0, 2}});
      histosMCreco.add("pT_track_WTS_PID_Pi", "pT with track selection and PID selection of Pi; pT [GeV/c]; Events", kTH1F, {{nBinsPt, 0, 2}});

      // Zero charge Event Transverse Momentum
      histosMCreco.add("pT_event_0charge_WTS_PID_Pi", "Event pT in 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {{nBinsPt, 0, 2}});

      // Non Zero charge Event Transverse Momentum
      histosMCreco.add("pT_event_non0charge_WTS_PID_Pi", "Event pT in Non 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {{nBinsPt, 0, 2}});

      // Rapidity of 0 charge Events
      histosMCreco.add("rapidity_event_0charge_WTS_PID_Pi_domainA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});
      histosMCreco.add("rapidity_event_0charge_WTS_PID_Pi_domainB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});
      histosMCreco.add("rapidity_event_0charge_WTS_PID_Pi_domainC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});

      // Rapidity of non 0 charge Events
      histosMCreco.add("rapidity_event_non0charge_WTS_PID_Pi_domainA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Events", kTH1F, {{nBinsRapidity, -2.5, 2.5}});
      histosMCreco.add("rapidity_event_non0charge_WTS_PID_Pi_domainB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c$; y; Events", kTH1F, {{nBinsRapidity, -2.5, 2.5}});
      histosMCreco.add("rapidity_event_non0charge_WTS_PID_Pi_domainC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Events", kTH1F, {{nBinsRapidity, -2.5, 2.5}});

      // Invariant Mass of 0 charge events
      histosMCreco.add("invMass_event_0charge_WTS_PID_Pi_domainA", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});       // pT < 0.15GeV
      histosMCreco.add("invMass_event_0charge_WTS_PID_Pi_domainB", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}}); // 0.15GeV < pT < 0.8GeV
      histosMCreco.add("invMass_event_0charge_WTS_PID_Pi_domainC", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});       // 0.8GeV < pT

      // Invariant mass of non 0 charge events
      histosMCreco.add("invMass_event_non0charge_WTS_PID_Pi_domainA", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});       // pT < 0.15GeV
      histosMCreco.add("invMass_event_non0charge_WTS_PID_Pi_domainB", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}}); // 0.15GeV < pT < 0.8GeV
      histosMCreco.add("invMass_event_non0charge_WTS_PID_Pi_domainC", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBinsInvariantMass, 0.8, 2.5}});       // 0.8GeV < pT

      // tpc signal
      histosMCreco.add("tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
      histosMCreco.add("tpcSignal_Pi", "TPC dEdx vs p for pions; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});

      // tof beta
      histosMCreco.add("tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
      histosMCreco.add("tofBeta_Pi", "TOF beta vs p for pions; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});

      // Other signals
      histosMCreco.add("FT0A", "T0A amplitude", kTH1F, {{200, 0.0, 500.0}});
      histosMCreco.add("FT0C", "T0C amplitude", kTH1F, {{200, 0.0, 500.0}});
      histosMCreco.add("ZDC_A", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
      histosMCreco.add("ZDC_C", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
      histosMCreco.add("V0A", "V0A amplitude", kTH1F, {{1000, 0.0, 100}});

      // Collin Soper Theta and Phi
      histosMCreco.add("CS_phi_pair_1", "#phi Distribution; #phi; Events", kTH1F, {{nBinsPhi, -3.2, 3.2}});
      histosMCreco.add("CS_phi_pair_2", "#phi Distribution; #phi; Events", kTH1F, {{nBinsPhi, -3.2, 3.2}});
      histosMCreco.add("CS_costheta_pair_1", "#theta Distribution;cos(#theta); Counts", kTH1F, {{nBinsCosTheta, -1, 1}});
      histosMCreco.add("CS_costheta_pair_2", "#theta Distribution;cos(#theta); Counts", kTH1F, {{nBinsCosTheta, -1, 1}});
    }

  } // End of init function
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Calculate the Collins-Soper Frame----------------------------------------------------------------------------------------------------------------------------
  double cosThetaCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1, ROOT::Math::PtEtaPhiMVector pair2, ROOT::Math::PtEtaPhiMVector fourpion)
  {
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);

    TLorentzVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target

    //  TVector3 beta = (-1. / fourpion.E()) * fourpion.Vect();
    ROOT::Math::PtEtaPhiMVector v1 = pair1;
    ROOT::Math::PtEtaPhiMVector v2 = pair2;
    ROOT::Math::PtEtaPhiMVector v12 = fourpion;

    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam1CM{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam2CM{(boostv12(pTargCM).Vect()).Unit()};

    // Axes
    ROOT::Math::XYZVectorF zaxisCS{((beam1CM.Unit() - beam2CM.Unit()).Unit())};

    double cosThetaCS = zaxisCS.Dot((v1CM));
    return cosThetaCS;
  } // End of cosThetaCollinsSoperFrame function------------------------------------------------------------------------------------------------------------------------

  // Calculate Phi in Collins-Soper Frame------------------------------------------------------------------------------------------------------------------------
  double phiCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1, ROOT::Math::PtEtaPhiMVector pair2, ROOT::Math::PtEtaPhiMVector fourpion)
  {
    // Half of the energy per pair of the colliding nucleons.
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);

    TLorentzVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target
    ROOT::Math::PtEtaPhiMVector v1 = pair1;
    ROOT::Math::PtEtaPhiMVector v2 = pair2;
    ROOT::Math::PtEtaPhiMVector v12 = fourpion;
    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam1CM{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam2CM{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxisCS{((beam1CM.Unit() - beam2CM.Unit()).Unit())};
    ROOT::Math::XYZVectorF yaxisCS{(beam1CM.Cross(beam2CM)).Unit()};
    ROOT::Math::XYZVectorF xaxisCS{(yaxisCS.Cross(zaxisCS)).Unit()};

    double phi = std::atan2(yaxisCS.Dot(v1CM), xaxisCS.Dot(v1CM));
    return phi;
  } // End of phiCollinsSoperFrame function------------------------------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  using UDtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Begin of Process function--------------------------------------------------------------------------------------------------------------------------------------------------
  void processData(UDCollisionFull const& collision, UDtracksfull const& tracks)
  {

    if (std::abs(collision.posZ()) > 10) {
      return;
    }

    int gapSide = collision.gapSide();
    float fitCuts[5] = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    std::vector<float> parameters = {pvCut, dcaZcut, dcaXYcut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, pTcut};
    int truegapSide = sgSelector.trueGap(collision, fitCuts[0], fitCuts[1], fitCuts[2], zdcCut);
    histosData.fill(HIST("GapSide"), gapSide);
    histosData.fill(HIST("TrueGapSide"), truegapSide);
    histosData.fill(HIST("EventCounts"), 1);
    gapSide = truegapSide;

    if ((gapSide != 2)) {
      return;
    }

    histosData.fill(HIST("V0A"), collision.totalFV0AmplitudeA());
    histosData.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
    histosData.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
    histosData.fill(HIST("ZDC_A"), collision.energyCommonZNA());
    histosData.fill(HIST("ZDC_C"), collision.energyCommonZNC());

    if (strictEventSelection) {
      if (collision.numContrib() != 4) {
        return;
      }
    } else {
      if (collision.numContrib() >= 10) {
        return;
      }
    }

    std::vector<decltype(tracks.begin())> WOTS_tracks;
    std::vector<decltype(tracks.begin())> WTS_tracks;
    std::vector<decltype(tracks.begin())> WTS_PID_Pi_tracks;
    std::vector<decltype(tracks.begin())> Pi_plus_tracks;
    std::vector<decltype(tracks.begin())> Pi_minus_tracks;

    for (const auto& t0 : tracks) {

      WOTS_tracks.push_back(t0);

      if (trackselector(t0, parameters)) {
        WTS_tracks.push_back(t0);

        if (selectionPIDPion(t0, true, nSigmaTPCcut, nSigmaTOFcut)) {
          WTS_PID_Pi_tracks.push_back(t0);
          if (t0.sign() == 1) {
            Pi_plus_tracks.push_back(t0);
          }
          if (t0.sign() == -1) {
            Pi_minus_tracks.push_back(t0);
          }
        } // End of Selection PID Pion

      } // End of track selections

    } // End of loop over tracks

    int numTracksWOTS = static_cast<int>(WOTS_tracks.size());
    int numTracksWTS = static_cast<int>(WTS_tracks.size());
    int numTracksWTSandPIDpi = static_cast<int>(WTS_PID_Pi_tracks.size());
    int numPiPlusTracks = static_cast<int>(Pi_plus_tracks.size());
    int numPionMinusTRacks = static_cast<int>(Pi_minus_tracks.size());

    for (int i = 0; i < numTracksWOTS; i++) {
      histosData.fill(HIST("tpcNSigmaPi_WOTS"), WOTS_tracks[i].tpcNSigmaPi(), std::sqrt(WOTS_tracks[i].px() * WOTS_tracks[i].px() + WOTS_tracks[i].py() * WOTS_tracks[i].py()));
      histosData.fill(HIST("tofNSigmaPi_WOTS"), WOTS_tracks[i].tofNSigmaPi(), std::sqrt(WOTS_tracks[i].px() * WOTS_tracks[i].px() + WOTS_tracks[i].py() * WOTS_tracks[i].py()));
      histosData.fill(HIST("pT_track_WOTS"), std::sqrt(WOTS_tracks[i].px() * WOTS_tracks[i].px() + WOTS_tracks[i].py() * WOTS_tracks[i].py()));
    } // End of loop over tracks without selection

    for (int i = 0; i < numTracksWTS; i++) {
      histosData.fill(HIST("tpcSignal"), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py() + WTS_tracks[i].pz() * WTS_tracks[i].pz()), WTS_tracks[i].tpcSignal());
      histosData.fill(HIST("tofBeta"), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py() + WTS_tracks[i].pz() * WTS_tracks[i].pz()), WTS_tracks[i].beta());
      histosData.fill(HIST("tpcNSigmaPi_WTS"), WTS_tracks[i].tpcNSigmaPi(), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py()));
      histosData.fill(HIST("tofNSigmaPi_WTS"), WTS_tracks[i].tofNSigmaPi(), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py()));
      histosData.fill(HIST("pT_track_WTS"), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py()));
    } // End of loop over tracks with selection only

    for (int i = 0; i < numTracksWTSandPIDpi; i++) {

      histosData.fill(HIST("tpcSignal_Pi"), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py() + WTS_PID_Pi_tracks[i].pz() * WTS_PID_Pi_tracks[i].pz()), WTS_PID_Pi_tracks[i].tpcSignal());
      histosData.fill(HIST("tofBeta_Pi"), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py() + WTS_PID_Pi_tracks[i].pz() * WTS_PID_Pi_tracks[i].pz()), WTS_PID_Pi_tracks[i].beta());

      histosData.fill(HIST("tpcNSigmaPi_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaPi(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosData.fill(HIST("tpcNSigmaKa_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaKa(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosData.fill(HIST("tpcNSigmaPr_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaPr(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosData.fill(HIST("tpcNSigmaEl_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaEl(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosData.fill(HIST("tpcNSigmaMu_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaMu(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));

      histosData.fill(HIST("tofNSigmaPi_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaPi(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosData.fill(HIST("tofNSigmaKa_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaKa(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosData.fill(HIST("tofNSigmaPr_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaPr(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosData.fill(HIST("tofNSigmaEl_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaEl(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosData.fill(HIST("tofNSigmaMu_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaMu(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));

      histosData.fill(HIST("pT_track_WTS_PID_Pi"), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
    } // End of loop over tracks with selection and PID selection of Pions

    if (numTracksWTSandPIDpi != 4) {
      return;
    }

    // Selecting Events with net charge = 0
    if (numPionMinusTRacks == 2 && numPiPlusTracks == 2) {

      TLorentzVector p1, p2, p3, p4, p1234;
      ROOT::Math::PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

      p1.SetXYZM(Pi_plus_tracks[0].px(), Pi_plus_tracks[0].py(), Pi_plus_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p2.SetXYZM(Pi_plus_tracks[1].px(), Pi_plus_tracks[1].py(), Pi_plus_tracks[1].pz(), o2::constants::physics::MassPionCharged);
      p3.SetXYZM(Pi_minus_tracks[0].px(), Pi_minus_tracks[0].py(), Pi_minus_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p4.SetXYZM(Pi_minus_tracks[1].px(), Pi_minus_tracks[1].py(), Pi_minus_tracks[1].pz(), o2::constants::physics::MassPionCharged);

      k1.SetCoordinates(p1.Pt(), p1.Eta(), p1.Phi(), o2::constants::physics::MassPionCharged);
      k2.SetCoordinates(p2.Pt(), p2.Eta(), p2.Phi(), o2::constants::physics::MassPionCharged);
      k3.SetCoordinates(p3.Pt(), p3.Eta(), p3.Phi(), o2::constants::physics::MassPionCharged);
      k4.SetCoordinates(p4.Pt(), p4.Eta(), p4.Phi(), o2::constants::physics::MassPionCharged);

      p1234 = p1 + p2 + p3 + p4;
      k1234 = k1 + k2 + k3 + k4;

      k13 = k1 + k3;
      k14 = k1 + k4;
      k23 = k2 + k3;
      k24 = k2 + k4;

      if (std::fabs(p1234.Rapidity()) < 0.5) {
        histosData.fill(HIST("pT_event_0charge_WTS_PID_Pi"), p1234.Pt());
        if (p1234.Pt() < 0.15) {
          histosData.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainA"), p1234.Rapidity());
          histosData.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainA"), p1234.M());

          auto phiPair1 = phiCollinsSoperFrame(k13, k24, k1234);
          auto phiPair2 = phiCollinsSoperFrame(k14, k23, k1234);
          auto cosThetaPair1 = cosThetaCollinsSoperFrame(k13, k24, k1234);
          auto cosThetaPair2 = cosThetaCollinsSoperFrame(k14, k23, k1234);

          histosData.fill(HIST("CS_phi_pair_1"), phiPair1);
          histosData.fill(HIST("CS_phi_pair_2"), phiPair2);
          histosData.fill(HIST("CS_costheta_pair_1"), cosThetaPair1);
          histosData.fill(HIST("CS_costheta_pair_2"), cosThetaPair2);
        }
        if (p1234.Pt() > 0.15 && p1234.Pt() < 0.80) {
          histosData.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainB"), p1234.Rapidity());
          histosData.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainB"), p1234.M());
        }
        if (p1234.Pt() > 0.80) {
          histosData.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainC"), p1234.Rapidity());
          histosData.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for 0 charge events

    // Selecting Events with net charge != 0 for estimation of background
    if (numPionMinusTRacks != 2 && numPiPlusTracks != 2) {

      TLorentzVector p1, p2, p3, p4, p1234;
      p1.SetXYZM(WTS_PID_Pi_tracks[0].px(), WTS_PID_Pi_tracks[0].py(), WTS_PID_Pi_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p2.SetXYZM(WTS_PID_Pi_tracks[1].px(), WTS_PID_Pi_tracks[1].py(), WTS_PID_Pi_tracks[1].pz(), o2::constants::physics::MassPionCharged);
      p3.SetXYZM(WTS_PID_Pi_tracks[2].px(), WTS_PID_Pi_tracks[2].py(), WTS_PID_Pi_tracks[2].pz(), o2::constants::physics::MassPionCharged);
      p4.SetXYZM(WTS_PID_Pi_tracks[3].px(), WTS_PID_Pi_tracks[3].py(), WTS_PID_Pi_tracks[3].pz(), o2::constants::physics::MassPionCharged);

      p1234 = p1 + p2 + p3 + p4;

      if (std::fabs(p1234.Rapidity()) < 0.5) {
        histosData.fill(HIST("pT_event_non0charge_WTS_PID_Pi"), p1234.Pt());

        if (p1234.Pt() < 0.15) {
          histosData.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainA"), p1234.Rapidity());
          histosData.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainA"), p1234.M());
        }
        if (p1234.Pt() > 0.15 && p1234.Pt() < 0.80) {
          histosData.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainB"), p1234.Rapidity());
          histosData.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainB"), p1234.M());
        }
        if (p1234.Pt() > 0.80) {
          histosData.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainC"), p1234.Rapidity());
          histosData.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for non 0 charge events

  } // End of 4 Pion Analysis Process function for Data
  PROCESS_SWITCH(exclusiveRhoTo4Pi, processData, "The Process for 4 Pion Analysis from data", ifDataAnalysis);
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Begin of MC Generation function-----------------------------------------------------------------------------------------------------------------------------------------------
  void processMCgen(aod::UDMcCollisions::iterator const&, aod::UDMcParticles const& mcParts)
  {
    std::vector<TLorentzVector> piPlusvectors;
    std::vector<TLorentzVector> piMinusvectors;
    TLorentzVector tempVector, p1, p2, p3, p4;
    TLorentzVector p1234;

    bool flag = false;

    for (const auto& particle : mcParts) {
      tempVector.SetXYZM(particle.px(), particle.py(), particle.pz(), o2::constants::physics::MassPionCharged);

      if (!particle.has_mothers()) {
        continue;
      }

      for (const auto& mother : particle.mothers_as<aod::UDMcParticles>()) {
        if (mother.pdgCode() == 30113) {

          if (flag == false) {
            histosMCgen.fill(HIST("rhoPrimeCounts"), 5);
          }
          flag = true;

          if (particle.pdgCode() == 211) {
            piPlusvectors.push_back(tempVector);
            histosMCgen.fill(HIST("MCgen_particle_pT"), tempVector.Pt());
            histosMCgen.fill(HIST("MCgen_particle_rapidity"), tempVector.Rapidity());
          }
          if (particle.pdgCode() == -211) {
            piMinusvectors.push_back(tempVector);
            histosMCgen.fill(HIST("MCgen_particle_pT"), tempVector.Pt());
            histosMCgen.fill(HIST("MCgen_particle_rapidity"), tempVector.Rapidity());
          }
        } // End of Mother ID 30113 rho prime
      } // End of loop over mothers
    } // End of loop over MC particles

    if (piPlusvectors.size() != 2 || piMinusvectors.size() != 2) {
      return;
    }

    p1234 = piPlusvectors[0] + piPlusvectors[1] + piMinusvectors[0] + piMinusvectors[1];
    histosMCgen.fill(HIST("MCgen_4pion_pT"), p1234.Pt());
    histosMCgen.fill(HIST("MCgen_4pion_rapidity"), p1234.Rapidity());
    histosMCgen.fill(HIST("MCgen_4pion_invmass"), p1234.M());

    ROOT::Math::PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

    k1.SetCoordinates(piPlusvectors[0].Pt(), piPlusvectors[0].Eta(), piPlusvectors[0].Phi(), o2::constants::physics::MassPionCharged);
    k2.SetCoordinates(piPlusvectors[1].Pt(), piPlusvectors[1].Eta(), piPlusvectors[1].Phi(), o2::constants::physics::MassPionCharged);
    k3.SetCoordinates(piMinusvectors[0].Pt(), piMinusvectors[0].Eta(), piMinusvectors[0].Phi(), o2::constants::physics::MassPionCharged);
    k4.SetCoordinates(piMinusvectors[1].Pt(), piMinusvectors[1].Eta(), piMinusvectors[1].Phi(), o2::constants::physics::MassPionCharged);

    k1234 = k1 + k2 + k3 + k4;

    k13 = k1 + k3;
    k14 = k1 + k4;
    k23 = k2 + k3;
    k24 = k2 + k4;

    auto phiPair1 = phiCollinsSoperFrame(k13, k24, k1234);
    auto phiPair2 = phiCollinsSoperFrame(k14, k23, k1234);
    auto cosThetaPair1 = cosThetaCollinsSoperFrame(k13, k24, k1234);
    auto cosThetaPair2 = cosThetaCollinsSoperFrame(k14, k23, k1234);

    histosMCgen.fill(HIST("MCgen_CS_phi_pair_1"), phiPair1);
    histosMCgen.fill(HIST("MCgen_CS_phi_pair_2"), phiPair2);
    histosMCgen.fill(HIST("MCgen_CS_costheta_pair_1"), cosThetaPair1);
    histosMCgen.fill(HIST("MCgen_CS_costheta_pair_2"), cosThetaPair2);

  } // End of 4 Pion MC Generation Process function
  PROCESS_SWITCH(exclusiveRhoTo4Pi, processMCgen, "The Process for 4 Pion Analysis from MC Generation", ifMCAnalysis);

  using CollisionStuff = soa::Join<aod::UDCollisions_001, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDMcCollsLabels>; //
  using CollisionTotal = CollisionStuff::iterator;
  using TrackStuff = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels>;

  void processMCrec(CollisionTotal const& collision, TrackStuff const& tracks)
  {

    if (std::abs(collision.posZ()) > 10) {
      return;
    }

    if (!collision.has_udMcCollision()) {
      return;
    }

    int gapSide = collision.gapSide();
    float fitCuts[5] = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    std::vector<float> parameters = {pvCut, dcaZcut, dcaXYcut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, pTcut};
    int truegapSide = sgSelector.trueGap(collision, fitCuts[0], fitCuts[1], fitCuts[2], zdcCut);
    histosMCreco.fill(HIST("GapSide"), gapSide);
    histosMCreco.fill(HIST("TrueGapSide"), truegapSide);
    histosMCreco.fill(HIST("EventCounts"), 1);
    gapSide = truegapSide;

    if ((gapSide != 2)) {
      return;
    }

    histosMCreco.fill(HIST("V0A"), collision.totalFV0AmplitudeA());
    histosMCreco.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
    histosMCreco.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
    histosMCreco.fill(HIST("ZDC_A"), collision.energyCommonZNA());
    histosMCreco.fill(HIST("ZDC_C"), collision.energyCommonZNC());

    if (strictEventSelection) {
      if (collision.numContrib() != 4) {
        return;
      }
    } else {
      if (collision.numContrib() >= 10) {
        return;
      }
    }

    std::vector<decltype(tracks.begin())> WOTS_tracks;
    std::vector<decltype(tracks.begin())> WTS_tracks;
    std::vector<decltype(tracks.begin())> WTS_PID_Pi_tracks;
    std::vector<decltype(tracks.begin())> Pi_plus_tracks;
    std::vector<decltype(tracks.begin())> Pi_minus_tracks;

    for (const auto& t0 : tracks) {

      WOTS_tracks.push_back(t0);

      if (trackselector(t0, parameters) && t0.has_udMcParticle()) {
        WTS_tracks.push_back(t0);

        if (selectionPIDPion(t0, true, nSigmaTPCcut, nSigmaTOFcut)) {
          WTS_PID_Pi_tracks.push_back(t0);
          if (t0.sign() == 1) {
            Pi_plus_tracks.push_back(t0);
          }
          if (t0.sign() == -1) {
            Pi_minus_tracks.push_back(t0);
          }
        } // End of Selection PID Pion

      } // End of track selections

    } // End of loop over tracks

    int numTracksWOTS = static_cast<int>(WOTS_tracks.size());
    int numTracksWTS = static_cast<int>(WTS_tracks.size());
    int numTracksWTSandPIDpi = static_cast<int>(WTS_PID_Pi_tracks.size());
    int numPiPlusTracks = static_cast<int>(Pi_plus_tracks.size());
    int numPionMinusTRacks = static_cast<int>(Pi_minus_tracks.size());

    for (int i = 0; i < numTracksWOTS; i++) {
      histosMCreco.fill(HIST("tpcNSigmaPi_WOTS"), WOTS_tracks[i].tpcNSigmaPi(), std::sqrt(WOTS_tracks[i].px() * WOTS_tracks[i].px() + WOTS_tracks[i].py() * WOTS_tracks[i].py()));
      histosMCreco.fill(HIST("tofNSigmaPi_WOTS"), WOTS_tracks[i].tofNSigmaPi(), std::sqrt(WOTS_tracks[i].px() * WOTS_tracks[i].px() + WOTS_tracks[i].py() * WOTS_tracks[i].py()));
      histosMCreco.fill(HIST("pT_track_WOTS"), std::sqrt(WOTS_tracks[i].px() * WOTS_tracks[i].px() + WOTS_tracks[i].py() * WOTS_tracks[i].py()));
    } // End of loop over tracks without selection

    for (int i = 0; i < numTracksWTS; i++) {
      histosMCreco.fill(HIST("tpcSignal"), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py() + WTS_tracks[i].pz() * WTS_tracks[i].pz()), WTS_tracks[i].tpcSignal());
      histosMCreco.fill(HIST("tofBeta"), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py() + WTS_tracks[i].pz() * WTS_tracks[i].pz()), WTS_tracks[i].beta());
      histosMCreco.fill(HIST("tpcNSigmaPi_WTS"), WTS_tracks[i].tpcNSigmaPi(), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py()));
      histosMCreco.fill(HIST("tofNSigmaPi_WTS"), WTS_tracks[i].tofNSigmaPi(), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py()));
      histosMCreco.fill(HIST("pT_track_WTS"), std::sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py()));
    } // End of loop over tracks with selection only

    for (int i = 0; i < numTracksWTSandPIDpi; i++) {

      histosMCreco.fill(HIST("tpcSignal_Pi"), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py() + WTS_PID_Pi_tracks[i].pz() * WTS_PID_Pi_tracks[i].pz()), WTS_PID_Pi_tracks[i].tpcSignal());
      histosMCreco.fill(HIST("tofBeta_Pi"), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py() + WTS_PID_Pi_tracks[i].pz() * WTS_PID_Pi_tracks[i].pz()), WTS_PID_Pi_tracks[i].beta());

      histosMCreco.fill(HIST("tpcNSigmaPi_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaPi(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosMCreco.fill(HIST("tpcNSigmaKa_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaKa(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosMCreco.fill(HIST("tpcNSigmaPr_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaPr(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosMCreco.fill(HIST("tpcNSigmaEl_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaEl(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosMCreco.fill(HIST("tpcNSigmaMu_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaMu(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));

      histosMCreco.fill(HIST("tofNSigmaPi_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaPi(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosMCreco.fill(HIST("tofNSigmaKa_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaKa(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosMCreco.fill(HIST("tofNSigmaPr_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaPr(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosMCreco.fill(HIST("tofNSigmaEl_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaEl(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
      histosMCreco.fill(HIST("tofNSigmaMu_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaMu(), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));

      histosMCreco.fill(HIST("pT_track_WTS_PID_Pi"), std::sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
    } // End of loop over tracks with selection and PID selection of Pions

    if (numTracksWTSandPIDpi != 4) {
      return;
    }

    // Selecting Events with net charge = 0
    if (numPionMinusTRacks == 2 && numPiPlusTracks == 2) {

      TLorentzVector p1, p2, p3, p4, p1234;
      ROOT::Math::PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

      p1.SetXYZM(Pi_plus_tracks[0].px(), Pi_plus_tracks[0].py(), Pi_plus_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p2.SetXYZM(Pi_plus_tracks[1].px(), Pi_plus_tracks[1].py(), Pi_plus_tracks[1].pz(), o2::constants::physics::MassPionCharged);
      p3.SetXYZM(Pi_minus_tracks[0].px(), Pi_minus_tracks[0].py(), Pi_minus_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p4.SetXYZM(Pi_minus_tracks[1].px(), Pi_minus_tracks[1].py(), Pi_minus_tracks[1].pz(), o2::constants::physics::MassPionCharged);

      k1.SetCoordinates(p1.Pt(), p1.Eta(), p1.Phi(), o2::constants::physics::MassPionCharged);
      k2.SetCoordinates(p2.Pt(), p2.Eta(), p2.Phi(), o2::constants::physics::MassPionCharged);
      k3.SetCoordinates(p3.Pt(), p3.Eta(), p3.Phi(), o2::constants::physics::MassPionCharged);
      k4.SetCoordinates(p4.Pt(), p4.Eta(), p4.Phi(), o2::constants::physics::MassPionCharged);

      p1234 = p1 + p2 + p3 + p4;
      k1234 = k1 + k2 + k3 + k4;

      k13 = k1 + k3;
      k14 = k1 + k4;
      k23 = k2 + k3;
      k24 = k2 + k4;

      if (std::fabs(p1234.Rapidity()) < 0.5) {
        histosMCreco.fill(HIST("pT_event_0charge_WTS_PID_Pi"), p1234.Pt());
        if (p1234.Pt() < 0.15) {
          histosMCreco.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainA"), p1234.Rapidity());
          histosMCreco.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainA"), p1234.M());

          auto phiPair1 = phiCollinsSoperFrame(k13, k24, k1234);
          auto phiPair2 = phiCollinsSoperFrame(k14, k23, k1234);
          auto cosThetaPair1 = cosThetaCollinsSoperFrame(k13, k24, k1234);
          auto cosThetaPair2 = cosThetaCollinsSoperFrame(k14, k23, k1234);

          histosMCreco.fill(HIST("CS_phi_pair_1"), phiPair1);
          histosMCreco.fill(HIST("CS_phi_pair_2"), phiPair2);
          histosMCreco.fill(HIST("CS_costheta_pair_1"), cosThetaPair1);
          histosMCreco.fill(HIST("CS_costheta_pair_2"), cosThetaPair2);
        }
        if (p1234.Pt() > 0.15 && p1234.Pt() < 0.80) {
          histosMCreco.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainB"), p1234.Rapidity());
          histosMCreco.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainB"), p1234.M());
        }
        if (p1234.Pt() > 0.80) {
          histosMCreco.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainC"), p1234.Rapidity());
          histosMCreco.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for 0 charge events

    // Selecting Events with net charge != 0 for estimation of background
    if (numPionMinusTRacks != 2 && numPiPlusTracks != 2) {

      TLorentzVector p1, p2, p3, p4, p1234;
      p1.SetXYZM(WTS_PID_Pi_tracks[0].px(), WTS_PID_Pi_tracks[0].py(), WTS_PID_Pi_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p2.SetXYZM(WTS_PID_Pi_tracks[1].px(), WTS_PID_Pi_tracks[1].py(), WTS_PID_Pi_tracks[1].pz(), o2::constants::physics::MassPionCharged);
      p3.SetXYZM(WTS_PID_Pi_tracks[2].px(), WTS_PID_Pi_tracks[2].py(), WTS_PID_Pi_tracks[2].pz(), o2::constants::physics::MassPionCharged);
      p4.SetXYZM(WTS_PID_Pi_tracks[3].px(), WTS_PID_Pi_tracks[3].py(), WTS_PID_Pi_tracks[3].pz(), o2::constants::physics::MassPionCharged);

      p1234 = p1 + p2 + p3 + p4;

      if (std::fabs(p1234.Rapidity()) < 0.5) {
        histosMCreco.fill(HIST("pT_event_non0charge_WTS_PID_Pi"), p1234.Pt());

        if (p1234.Pt() < 0.15) {
          histosMCreco.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainA"), p1234.Rapidity());
          histosMCreco.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainA"), p1234.M());
        }
        if (p1234.Pt() > 0.15 && p1234.Pt() < 0.80) {
          histosMCreco.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainB"), p1234.Rapidity());
          histosMCreco.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainB"), p1234.M());
        }
        if (p1234.Pt() > 0.80) {
          histosMCreco.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainC"), p1234.Rapidity());
          histosMCreco.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for non 0 charge events

  } // End of 4 Pion Analysis Process function for MC Reconstruction
  PROCESS_SWITCH(exclusiveRhoTo4Pi, processMCrec, "The Process for 4 Pion Analysis from MC Reconstruction", ifMCAnalysis);

}; // End of Struct exclusiveRhoTo4Pi
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<exclusiveRhoTo4Pi>(cfgc)};
}
