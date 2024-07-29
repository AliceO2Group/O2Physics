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
///
/// \brief  task for an analysis of rho photoproduction in UPCs, intended usage is with UD tables
///         includes event tagging based on ZN information, track selection, reconstruction,
///         and also some basic stuff for phi anisotropy studies
/// \author Jakub Juracka, jakub.juracka@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h" // has some useful funtions for stuff not available from the tables

// ROOT headers
#include <Math/Vector4D.h> // these should apparently be used instead of TLorentzVector
#include <Math/Vector2D.h>

#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <chrono>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>::iterator;
using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;

struct upcRhoAnalysis {
  // configurables
  double PcEtaCut = 0.9; // cut on track eta as per physics coordination recommendation
  Configurable<bool> tracksRequireTOF{"tracksRequireTOF", false, "requireTOF"};
  Configurable<double> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "tpcNSigmaPiCut"};
  Configurable<double> tracksTpcNSigmaElCut{"tracksTpcNSigmaElCut", 3.0, "tpcNSigmaElCut"};
  Configurable<double> tracksTofNSigmaPiCut{"tracksTofNSigmaPiCut", 3.0, "tofNSigmaPiCut"};
  Configurable<double> tracksPtMaxCut{"tracksPtMaxCut", 1.0, "ptMaxCut"};
  Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "dcaMaxCut"};
  Configurable<double> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "posZMaxCut"};
  Configurable<double> ZNcommonEnergyCut{"ZNcommonEnergyCut", 0.0, "ZNcommonEnergyCut"};
  Configurable<double> ZNtimeCut{"ZNtimeCut", 2.0, "ZNtimeCut"};
  Configurable<double> systemMassMinCut{"2systemMassMinCut", 0.5, "2systemMassMinCut"};
  Configurable<double> systemMassMaxCut{"2systemMassMaxCut", 1.1, "2systemMassMaxCut"};
  Configurable<double> systemPtCut{"2systemPtMaxCut", 0.12, "2systemPtMaxCut"};
  Configurable<double> systemYCut{"2systemYCut", 0.9, "2systemYCut"};

  ConfigurableAxis mAxis{"mAxis", {150, 0.0, 1.5}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {500, 0.0, 5.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {2000, 0.0, 0.2}, "#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis etaAxis{"etaAxis", {180, -0.9, 0.9}, "#eta"};
  ConfigurableAxis yAxis{"yAxis", {180, -0.9, 0.9}, "y"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, 2.0*o2::constants::math::PI}, "#phi"};
  ConfigurableAxis phiAssymAxis{"phiAssymAxis", {180, 0, o2::constants::math::PI}, "#phi"};
  ConfigurableAxis nTracksAxis{"nTracksAxis", {101, -0.5, 100.5}, "N_{tracks}"};
  ConfigurableAxis tpcNSigmaPiAxis{"tpcNSigmaPiAxis", {400, -10.0, 30.0}, "TPC n#sigma_{#pi}"};
  ConfigurableAxis tofNSigmaPiAxis{"tofNSigmaPiAxis", {400, -20.0, 20.0}, "TOF n#sigma_{#pi}"};
  ConfigurableAxis dcaAxis{"dcaXYAxis", {1000, -5.0, 5.0}, "DCA (cm)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&) {
    // QA //
    // collisions
    registry.add("QC/collisions/hNetCharge", ";net charge;counts", kTH1D, {{11, -5.5, 5.5}});
    registry.add("QC/collisions/hNumContributors", ";N_{contributors};counts", kTH1D, {{11, -0.5, 10.5}});
    registry.add("QC/collisions/hPosXY", ";x (cm);y (cm);counts", kTH2D, {{2000, -0.1, 0.1},{2000, -0.1, 0.1}});
    registry.add("QC/collisions/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    registry.add("QC/collisions/hZDCcommonEnergy", ";ZNA common energy;ZNC common energy;counts", kTH2D, {{250, -5.0, 20.0},{250, -5.0, 20.0}});
    registry.add("QC/collisions/hZDCtime", ";ZNA time;ZNC time;counts", kTH2D, {{200, -10.0, 10.0},{200, -10.0, 10.0}});
    // all tracks
    registry.add("QC/allTracks/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {tpcNSigmaPiAxis});
    registry.add("QC/allTracks/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {tofNSigmaPiAxis});
    registry.add("QC/allTracks/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {dcaAxis, dcaAxis});
    registry.add("QC/allTracks/hTpcNsigmaPi2D", ";TPC n#sigma_{#pi^{+}};TPC n#sigma_{#pi^{-}};counts", kTH2D, {tpcNSigmaPiAxis, tpcNSigmaPiAxis});
    // tracks passing selections
    registry.add("QC/cutTracks/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {tpcNSigmaPiAxis});
    registry.add("QC/cutTracks/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {tofNSigmaPiAxis});
    registry.add("QC/cutTracks/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {dcaAxis, dcaAxis});
    registry.add("QC/cutTracks/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    registry.add("QC/cutTracks/hTpcNsigmaPi2D", ";TPC n#sigma_{#pi^{+}};TPC n#sigma_{#pi^{-}};counts", kTH2D, {tpcNSigmaPiAxis, tpcNSigmaPiAxis});
    
    // RECO HISTOS //
    // PIONS
    // no selection
    registry.add("reco/pions/no-selection/unlike-sign/hPt", ";p_{T}(#pi^{+}) (GeV/#it{c});p_{T}(#pi^{-}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/no-selection/unlike-sign/hEta", ";#eta(#pi^{+});#eta(#pi^{-});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/no-selection/unlike-sign/hPhi", ";#phi(#pi^{+});#phi(#pi^{-});counts", kTH2D, {phiAxis, phiAxis});
    registry.add("reco/pions/no-selection/like-sign/hPt", ";p_{T}(#pi_{1}) (GeV/#it{c});p_{T}(#pi_{2}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/no-selection/like-sign/hEta", ";#eta(#pi_{1});#eta(#pi_{2});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/no-selection/like-sign/hPhi", ";#phi(#pi_{1});#phi(#pi_{2});counts", kTH2D, {phiAxis, phiAxis});
    // 0n0n
    registry.add("reco/pions/0n0n/unlike-sign/hPt", ";p_{T}(#pi^{+}) (GeV/#it{c});p_{T}(#pi^{-}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/0n0n/unlike-sign/hEta", ";#eta(#pi^{+});#eta(#pi^{-});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/0n0n/unlike-sign/hPhi", ";#phi(#pi^{+});#phi(#pi^{-});counts", kTH2D, {phiAxis, phiAxis});
    registry.add("reco/pions/0n0n/like-sign/hPt", ";p_{T}(#pi_{1}) (GeV/#it{c});p_{T}(#pi_{2}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/0n0n/like-sign/hEta", ";#eta(#pi_{1});#eta(#pi_{2});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/0n0n/like-sign/hPhi", ";#phi(#pi_{1});#phi(#pi_{2});counts", kTH2D, {phiAxis, phiAxis});
    // Xn0n
    registry.add("reco/pions/Xn0n/unlike-sign/hPt", ";p_{T}(#pi^{+}) (GeV/#it{c});p_{T}(#pi^{-}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/Xn0n/unlike-sign/hEta", ";#eta(#pi^{+});#eta(#pi^{-});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/Xn0n/unlike-sign/hPhi", ";#phi(#pi^{+});#phi(#pi^{-});counts", kTH2D, {phiAxis, phiAxis});
    registry.add("reco/pions/Xn0n/like-sign/hPt", ";p_{T}(#pi_{1}) (GeV/#it{c});p_{T}(#pi_{2}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/Xn0n/like-sign/hEta", ";#eta(#pi_{1});#eta(#pi_{2});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/Xn0n/like-sign/hPhi", ";#phi(#pi_{1});#phi(#pi_{2});counts", kTH2D, {phiAxis, phiAxis});
    // XnXn
    registry.add("reco/pions/XnXn/unlike-sign/hPt", ";p_{T}(#pi^{+}) (GeV/#it{c});p_{T}(#pi^{-}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/XnXn/unlike-sign/hEta", ";#eta(#pi^{+});#eta(#pi^{-});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/XnXn/unlike-sign/hPhi", ";#phi(#pi^{+});#phi(#pi^{-});counts", kTH2D, {phiAxis, phiAxis});
    registry.add("reco/pions/XnXn/like-sign/hPt", ";p_{T}(#pi_{1}) (GeV/#it{c});p_{T}(#pi_{2}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/XnXn/like-sign/hEta", ";#eta(#pi_{1});#eta(#pi_{2});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/XnXn/like-sign/hPhi", ";#phi(#pi_{1});#phi(#pi_{2});counts", kTH2D, {phiAxis, phiAxis});
    
    // RHOS
    // no selection
    registry.add("reco/system/2pi/no-selection/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/no-selection/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/no-selection/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/no-selection/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/no-selection/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/no-selection/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/no-selection/unlike-sign/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/no-selection/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/positive/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/no-selection/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/no-selection/like-sign/negative/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    // 0n0n
    registry.add("reco/system/2pi/0n0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/0n0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/0n0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/0n0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/0n0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/0n0n/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/0n0n/unlike-sign/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/0n0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/positive/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/0n0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/0n0n/like-sign/negative/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    // Xn0n
    registry.add("reco/system/2pi/Xn0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/Xn0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/Xn0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/Xn0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/Xn0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/Xn0n/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/Xn0n/unlike-sign/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/Xn0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/positive/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/Xn0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/Xn0n/like-sign/negative/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    // XnXn
    registry.add("reco/system/2pi/XnXn/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/XnXn/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/XnXn/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/XnXn/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/XnXn/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/XnXn/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/XnXn/unlike-sign/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/XnXn/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/positive/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/XnXn/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/XnXn/like-sign/negative/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});

    // SELECTED RHOS
    // no selection
    registry.add("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    // 0n0n
    registry.add("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    // Xn0n
    registry.add("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    // XnXn
    registry.add("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPhiAssym", ";#phi;counts", kTH1D, {phiAssymAxis});
    registry.add("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPtPhiProjection", ";p_{T}cos#phi;p_{T}sin#phi;counts", kTH2D, {{500, -0.25, 0.25},{250, 0.0, 0.25}});
    
    // 4PI AND 6PI SYSTEM
    registry.add("reco/system/4pi/net-zero/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {{800, 0.0, 8.0}});
    registry.add("reco/system/4pi/net-zero/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/4pi/net-zero/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/4pi/net-zero/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/4pi/net-zero/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/4pi/non-net-zero/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {{800, 0.0, 8.0}});
    registry.add("reco/system/4pi/non-net-zero/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/4pi/non-net-zero/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/4pi/non-net-zero/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/4pi/non-net-zero/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/6pi/net-zero/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {{800, 0.0, 8.0}});
    registry.add("reco/system/6pi/net-zero/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/6pi/net-zero/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/6pi/net-zero/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/6pi/net-zero/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/6pi/non-net-zero/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {{800, 0.0, 8.0}});
    registry.add("reco/system/6pi/non-net-zero/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/6pi/non-net-zero/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/6pi/non-net-zero/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/6pi/non-net-zero/hPhi", ";#phi;counts", kTH1D, {phiAxis});
  }

  template <typename T>
  bool trackPassesCuts(T const& track) { // check if track passes preliminary cuts (PID done separately)
    if (!track.isPVContributor()) return false;
    if (!track.hasITS()) return false;
    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > tracksDcaMaxCut) return false;
    if (std::abs(eta(track.px(), track.py(), track.pz())) > PcEtaCut) return false;
    if (std::abs(track.pt()) > tracksPtMaxCut) return false;
    return true;
  }

  template <typename T>
  bool tracksPassPiPID(const T& cutTracks) { // check if pre-cut tracks pass PID cut (fall within n-dimensional hypersphere)
    double radius = 0.0;
    for (const auto& track : cutTracks) radius += std::pow(track.tpcNSigmaPi(), 2);
    return radius < std::pow(tracksTpcNSigmaPiCut, 2);
  }

  template <typename T>
  double tracksTotalCharge(const T& cutTracks) { // calculate total charge of selected tracks
    double charge = 0.0;
    for (const auto& track : cutTracks) charge += track.sign();
    return charge;
  }

  template <typename T>
  bool systemPassCuts(const T& system) { // check if system passes system cuts
    if (system.M() < systemMassMinCut || system.M() > systemMassMaxCut) return false;
    if (system.Pt() > systemPtCut) return false;
    if (std::abs(system.Rapidity()) > systemYCut) return false;
    return true;
  }

  ROOT::Math::PxPyPzMVector reconstructSystem(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracks4Vecs) { // reconstruct system from 4-vectors
    ROOT::Math::PxPyPzMVector system;
    for (const auto& track4Vec : cutTracks4Vecs) system += track4Vec;
    return system;
  }

  template <typename T>
  double getPhiAssym(const T& cutTracks) { // phi anisotropy based on STAR + ALICE papers
    // assign the tracks randomly
    std::vector<int> indices = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // get time-based seed
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    ROOT::Math::XYVector pOne(cutTracks[indices[0]].px(), cutTracks[indices[0]].py());
    ROOT::Math::XYVector pTwo(cutTracks[indices[1]].px(), cutTracks[indices[1]].py());
    auto pPlus = pOne + pTwo;
    auto pMinus = pOne - pTwo;
    // no method for direct calculation of angle -> use dot product formula
    double cosPhi = (pPlus.Dot(pMinus))/(std::sqrt(pPlus.Mag2()) * std::sqrt(pMinus.Mag2()));
    return std::acos(cosPhi);
  }

  void processReco(FullUDCollision const& collision, FullUDTracks const& tracks, aod::UDZdcsReduced const& zdcs) {
    // QC histograms
    registry.fill(HIST("QC/collisions/hNetCharge"), collision.netCharge());
    registry.fill(HIST("QC/collisions/hNumContributors"), collision.numContrib());
    registry.fill(HIST("QC/collisions/hPosXY"), collision.posX(), collision.posY());
    registry.fill(HIST("QC/collisions/hPosZ"), collision.posZ());
    registry.fill(HIST("QC/collisions/hZDCcommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    registry.fill(HIST("QC/collisions/hZDCtime"), collision.timeZNA(), collision.timeZNC());

    // apply some cuts on collisions
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut) return;
    // event tagging
    bool XnXn = false, OnOn = false, XnOn = false; // note: OnOn == 0n0n
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() < ZNcommonEnergyCut) OnOn = true;
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut &&
        collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) XnXn = true;
    if ((collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) || 
        (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() < ZNcommonEnergyCut)) XnOn = true;

    // vectors for storing selected tracks and their 4-vectors
    std::vector<decltype(tracks.begin())> cutTracks;
    std::vector<ROOT::Math::PxPyPzMVector> cutTracks4Vecs;

    for (const auto& track : tracks) {
      registry.fill(HIST("QC/allTracks/hTpcNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("QC/allTracks/hTofNSigmaPi"), track.tofNSigmaPi());
      registry.fill(HIST("QC/allTracks/hDcaXYZ"), track.dcaZ(), track.dcaXY());

      if (!trackPassesCuts(track)) continue; // apply cuts
      cutTracks.push_back(track);
      cutTracks4Vecs.push_back(ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged)); // apriori assume pion mass
      
      registry.fill(HIST("QC/cutTracks/hTpcNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("QC/cutTracks/hTofNSigmaPi"), track.tofNSigmaPi());
      registry.fill(HIST("QC/cutTracks/hDcaXYZ"), track.dcaZ(), track.dcaXY());
      registry.fill(HIST("QC/cutTracks/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
    }

    if (!tracksPassPiPID(cutTracks)) return;
    auto system = reconstructSystem(cutTracks4Vecs); // reconstruct system from 4-vectors

    if (cutTracks.size() == 2) {
      ROOT::Math::PxPyPzMVector piPos, piNeg;
      decltype(tracks.begin()) posTrack, negTrack;

      // unlike-sign tracks
      if (tracksTotalCharge(cutTracks) == 0) {
        if (cutTracks[0].sign() > 0) {
          piPos = cutTracks4Vecs[0];
          posTrack = cutTracks[0];
          piNeg = cutTracks4Vecs[1];
          negTrack = cutTracks[1];
        } else {
          piNeg = cutTracks4Vecs[0];
          negTrack = cutTracks[0];
          piPos = cutTracks4Vecs[1];
          posTrack = cutTracks[1];
        }

        // no selection on ZDC energy
        registry.fill(HIST("reco/pions/no-selection/unlike-sign/hPt"), piPos.Pt(), piNeg.Pt());
        registry.fill(HIST("reco/pions/no-selection/unlike-sign/hEta"), piPos.Eta(), piNeg.Eta());
        registry.fill(HIST("reco/pions/no-selection/unlike-sign/hPhi"), piPos.Phi() + o2::constants::math::PI, piNeg.Phi() + o2::constants::math::PI);
        registry.fill(HIST("reco/system/2pi/no-selection/unlike-sign/hM"), system.M());
        registry.fill(HIST("reco/system/2pi/no-selection/unlike-sign/hPt"), system.Pt());
        registry.fill(HIST("reco/system/2pi/no-selection/unlike-sign/hPt2"), system.Pt()*system.Pt());
        registry.fill(HIST("reco/system/2pi/no-selection/unlike-sign/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/2pi/no-selection/unlike-sign/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/2pi/no-selection/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
        registry.fill(HIST("reco/system/2pi/no-selection/unlike-sign/hPhiAssym"), getPhiAssym(cutTracks));
        if (systemPassCuts(system)) { // ful cuts
          registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPhiAssym"), getPhiAssym(cutTracks));
          registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
        }
        if (system.M() > systemMassMinCut && system.M() < systemMassMaxCut && std::abs(system.Rapidity()) < systemYCut) // only m and y cut
          registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/unlike-sign/hPt2"), system.Pt()*system.Pt());

        // selection on ZDC energy
        if (OnOn) {
          registry.fill(HIST("reco/pions/0n0n/unlike-sign/hPt"), piPos.Pt(), piNeg.Pt());
          registry.fill(HIST("reco/pions/0n0n/unlike-sign/hEta"), piPos.Eta(), piNeg.Eta());
          registry.fill(HIST("reco/pions/0n0n/unlike-sign/hPhi"), piPos.Phi() + o2::constants::math::PI, piNeg.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/0n0n/unlike-sign/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/0n0n/unlike-sign/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/0n0n/unlike-sign/hPt2"), system.Pt()*system.Pt());
          registry.fill(HIST("reco/system/2pi/0n0n/unlike-sign/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/0n0n/unlike-sign/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/0n0n/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/0n0n/unlike-sign/hPhiAssym"), getPhiAssym(cutTracks));
          if (systemPassCuts(system)) {
            registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPhiAssym"), getPhiAssym(cutTracks));
            registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/unlike-sign/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
          }
        } else if (XnOn) {
          registry.fill(HIST("reco/pions/Xn0n/unlike-sign/hPt"), piPos.Pt(), piNeg.Pt());
          registry.fill(HIST("reco/pions/Xn0n/unlike-sign/hEta"), piPos.Eta(), piNeg.Eta());
          registry.fill(HIST("reco/pions/Xn0n/unlike-sign/hPhi"), piPos.Phi() + o2::constants::math::PI, piNeg.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/Xn0n/unlike-sign/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/Xn0n/unlike-sign/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/Xn0n/unlike-sign/hPt2"), system.Pt()*system.Pt());
          registry.fill(HIST("reco/system/2pi/Xn0n/unlike-sign/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/Xn0n/unlike-sign/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/Xn0n/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/Xn0n/unlike-sign/hPhiAssym"), getPhiAssym(cutTracks));
          if (systemPassCuts(system)) {
            registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPhiAssym"), getPhiAssym(cutTracks));
            registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/unlike-sign/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
          }
        } else if (XnXn) {
          registry.fill(HIST("reco/pions/XnXn/unlike-sign/hPt"), piPos.Pt(), piNeg.Pt());
          registry.fill(HIST("reco/pions/XnXn/unlike-sign/hEta"), piPos.Eta(), piNeg.Eta());
          registry.fill(HIST("reco/pions/XnXn/unlike-sign/hPhi"), piPos.Phi() + o2::constants::math::PI, piNeg.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/XnXn/unlike-sign/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/XnXn/unlike-sign/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/XnXn/unlike-sign/hPt2"), system.Pt()*system.Pt());
          registry.fill(HIST("reco/system/2pi/XnXn/unlike-sign/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/XnXn/unlike-sign/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/XnXn/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/XnXn/unlike-sign/hPhiAssym"), getPhiAssym(cutTracks));
          if (systemPassCuts(system)) {
            registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPhiAssym"), getPhiAssym(cutTracks));
            registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/unlike-sign/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
          }
        }
      }

      // like-sign tracks PIDed as pions
      if (tracksTotalCharge(cutTracks) != 0) {
        // no selection on ZDC energy
        registry.fill(HIST("reco/pions/no-selection/like-sign/hPt"), cutTracks4Vecs[0].Pt(), cutTracks4Vecs[1].Pt());
        registry.fill(HIST("reco/pions/no-selection/like-sign/hEta"), cutTracks4Vecs[0].Eta(), cutTracks4Vecs[1].Eta());
        registry.fill(HIST("reco/pions/no-selection/like-sign/hPhi"), cutTracks4Vecs[0].Phi() + o2::constants::math::PI, cutTracks4Vecs[1].Phi() + o2::constants::math::PI);
        if (tracksTotalCharge(cutTracks) > 0) {
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/positive/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/positive/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/positive/hPt2"), system.Pt()*system.Pt());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/positive/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/positive/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/positive/hPhi"), system.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/positive/hPhiAssym"), getPhiAssym(cutTracks));
          if (systemPassCuts(system)) {
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPhiAssym"), getPhiAssym(cutTracks));
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/positive/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
          }
        } else {
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/negative/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/negative/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/negative/hPt2"), system.Pt()*system.Pt());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/negative/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/negative/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/negative/hPhi"), system.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2pi/no-selection/like-sign/negative/hPhiAssym"), getPhiAssym(cutTracks));
          if (systemPassCuts(system)) {
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPhiAssym"), getPhiAssym(cutTracks));
            registry.fill(HIST("reco/system/2pi/pass-cuts/no-selection/like-sign/negative/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
          }
        }

        // selection on ZDC energy
        if (OnOn) {
          registry.fill(HIST("reco/pions/0n0n/like-sign/hPt"), cutTracks4Vecs[0].Pt(), cutTracks4Vecs[1].Pt());
          registry.fill(HIST("reco/pions/0n0n/like-sign/hEta"), cutTracks4Vecs[0].Eta(), cutTracks4Vecs[1].Eta());
          registry.fill(HIST("reco/pions/0n0n/like-sign/hPhi"), cutTracks4Vecs[0].Phi() + o2::constants::math::PI, cutTracks4Vecs[1].Phi() + o2::constants::math::PI);
          if (tracksTotalCharge(cutTracks) > 0) {
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/positive/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/positive/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/positive/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/positive/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/positive/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/positive/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/positive/hPhiAssym"), getPhiAssym(cutTracks));
            if (systemPassCuts(system)) {
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hM"), system.M());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPt"), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPt2"), system.Pt()*system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPtVsM"), system.M(), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hY"), system.Rapidity());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPhi"), system.Phi() + o2::constants::math::PI);
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPhiAssym"), getPhiAssym(cutTracks));
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/positive/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
            }
          } else {
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/negative/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/negative/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/negative/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/negative/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/negative/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/negative/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/0n0n/like-sign/negative/hPhiAssym"), getPhiAssym(cutTracks));
            if (systemPassCuts(system)) {
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hM"), system.M());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPt"), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPt2"), system.Pt()*system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPtVsM"), system.M(), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hY"), system.Rapidity());
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPhi"), system.Phi() + o2::constants::math::PI);
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPhiAssym"), getPhiAssym(cutTracks));
              registry.fill(HIST("reco/system/2pi/pass-cuts/0n0n/like-sign/negative/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
            }
          }
        } else if (XnOn) {
          registry.fill(HIST("reco/pions/Xn0n/like-sign/hPt"), cutTracks4Vecs[0].Pt(), cutTracks4Vecs[1].Pt());
          registry.fill(HIST("reco/pions/Xn0n/like-sign/hEta"), cutTracks4Vecs[0].Eta(), cutTracks4Vecs[1].Eta());
          registry.fill(HIST("reco/pions/Xn0n/like-sign/hPhi"), cutTracks4Vecs[0].Phi() + o2::constants::math::PI, cutTracks4Vecs[1].Phi() + o2::constants::math::PI);
          if (tracksTotalCharge(cutTracks) > 0) {
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/positive/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/positive/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/positive/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/positive/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/positive/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/positive/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/positive/hPhiAssym"), getPhiAssym(cutTracks));
            if (systemPassCuts(system)) {
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hM"), system.M());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPt"), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPt2"), system.Pt()*system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPtVsM"), system.M(), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hY"), system.Rapidity());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPhi"), system.Phi() + o2::constants::math::PI);
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPhiAssym"), getPhiAssym(cutTracks));
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/positive/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
            }
          } else {
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/negative/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/negative/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/negative/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/negative/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/negative/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/negative/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/Xn0n/like-sign/negative/hPhiAssym"), getPhiAssym(cutTracks));
            if (systemPassCuts(system)) {
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hM"), system.M());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPt"), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPt2"), system.Pt()*system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPtVsM"), system.M(), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hY"), system.Rapidity());
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPhi"), system.Phi() + o2::constants::math::PI);
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPhiAssym"), getPhiAssym(cutTracks));
              registry.fill(HIST("reco/system/2pi/pass-cuts/Xn0n/like-sign/negative/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
            }
          }
        } else if (XnXn) {
          registry.fill(HIST("reco/pions/XnXn/like-sign/hPt"), cutTracks4Vecs[0].Pt(), cutTracks4Vecs[1].Pt());
          registry.fill(HIST("reco/pions/XnXn/like-sign/hEta"), cutTracks4Vecs[0].Eta(), cutTracks4Vecs[1].Eta());
          registry.fill(HIST("reco/pions/XnXn/like-sign/hPhi"), cutTracks4Vecs[0].Phi() + o2::constants::math::PI, cutTracks4Vecs[1].Phi() + o2::constants::math::PI);
          if (tracksTotalCharge(cutTracks) > 0) {
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/positive/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/positive/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/positive/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/positive/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/positive/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/positive/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/positive/hPhiAssym"), getPhiAssym(cutTracks));
            if (systemPassCuts(system)) {
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hM"), system.M());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPt"), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPt2"), system.Pt()*system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPtVsM"), system.M(), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hY"), system.Rapidity());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPhi"), system.Phi() + o2::constants::math::PI);
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPhiAssym"), getPhiAssym(cutTracks));
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/positive/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
            }
          } else {
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/negative/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/negative/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/negative/hPt2"), system.Pt()*system.Pt());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/negative/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/negative/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/negative/hPhi"), system.Phi() + o2::constants::math::PI);
            registry.fill(HIST("reco/system/2pi/XnXn/like-sign/negative/hPhiAssym"), getPhiAssym(cutTracks));
            if (systemPassCuts(system)) {
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hM"), system.M());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPt"), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPt2"), system.Pt()*system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPtVsM"), system.M(), system.Pt());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hY"), system.Rapidity());
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPhi"), system.Phi() + o2::constants::math::PI);
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPhiAssym"), getPhiAssym(cutTracks));
              registry.fill(HIST("reco/system/2pi/pass-cuts/XnXn/like-sign/negative/hPtPhiProjection"), system.Pt()*std::cos(getPhiAssym(cutTracks)),system.Pt()*std::sin(getPhiAssym(cutTracks)));
            }
          }
        }
      }
    }

    if (cutTracks.size() == 4) {
      if (tracksTotalCharge(cutTracks) == 0) {
        registry.fill(HIST("reco/system/4pi/net-zero/hM"), system.M());
        registry.fill(HIST("reco/system/4pi/net-zero/hPt"), system.Pt());
        registry.fill(HIST("reco/system/4pi/net-zero/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/4pi/net-zero/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/4pi/net-zero/hPhi"), system.Phi() + o2::constants::math::PI);
      } else {
        registry.fill(HIST("reco/system/4pi/non-net-zero/hM"), system.M());
        registry.fill(HIST("reco/system/4pi/non-net-zero/hPt"), system.Pt());
        registry.fill(HIST("reco/system/4pi/non-net-zero/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/4pi/non-net-zero/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/4pi/non-net-zero/hPhi"), system.Phi() + o2::constants::math::PI);
      }
    }

    if (cutTracks.size() == 6) {
      if (tracksTotalCharge(cutTracks) == 0) {
        registry.fill(HIST("reco/system/6pi/net-zero/hM"), system.M());
        registry.fill(HIST("reco/system/6pi/net-zero/hPt"), system.Pt());
        registry.fill(HIST("reco/system/6pi/net-zero/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/6pi/net-zero/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/6pi/net-zero/hPhi"), system.Phi() + o2::constants::math::PI);
      } else {
        registry.fill(HIST("reco/system/6pi/non-net-zero/hM"), system.M());
        registry.fill(HIST("reco/system/6pi/non-net-zero/hPt"), system.Pt());
        registry.fill(HIST("reco/system/6pi/non-net-zero/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/6pi/non-net-zero/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/6pi/non-net-zero/hPhi"), system.Phi() + o2::constants::math::PI);
      }
    }
  } PROCESS_SWITCH(upcRhoAnalysis, processReco, "analyse reco tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<upcRhoAnalysis>(cfgc)
  };
}