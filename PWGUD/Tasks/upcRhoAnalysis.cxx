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
/// \brief  task for analysis of rho in UPCs using UD tables (from SG producer)
///         includes event tagging based on ZN information, track selection, reconstruction,
///         and also some basic stuff for decay phi anisotropy studies
/// \author Jakub Juracka, jakub.juracka@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h" // used instead of TLorentzVector
#include "Math/Vector2D.h"
#include "random"

#include "Common/DataModel/PIDResponse.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUDSgCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::SGCollisions>::iterator;
using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;

struct upcRhoAnalysis {
  double PcEtaCut = 0.9; // physics coordination recommendation
  Configurable<bool> specifyGapSide{"specifyGapSide", true, "specify gap side for SG/DG produced data"};
  Configurable<double> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "TPC nSigma pion cut"};
  Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "max DCA cut on tracks"};
  Configurable<double> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "max Z position cut on collisions"};
  Configurable<double> ZNcommonEnergyCut{"ZNcommonEnergyCut", 0.0, "ZN common energy cut"};
  Configurable<double> ZNtimeCut{"ZNtimeCut", 2.0, "ZN time cut"};
  Configurable<double> systemMassMinCut{"systemMassMinCut", 0.5, "min M cut for reco system"};
  Configurable<double> systemMassMaxCut{"systemMassMaxCut", 1.2, "max M cut for reco system"};
  Configurable<double> systemPtCut{"systemPtMaxCut", 0.3, "max pT cut for reco system"};
  Configurable<double> systemYCut{"systemYCut", 0.9, "rapiditiy cut for reco system"};

  ConfigurableAxis mAxis{"mAxis", {1000, 0.0, 10.0}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis mCutAxis{"mCutAxis", {70, 0.5, 1.2}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {1000, 0.0, 10.0}, "p_{T} (GeV/#it{c})"};
  ConfigurableAxis ptCutAxis{"ptCutAxis", {30, 0.0, 0.3}, "p_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {90, 0.0, 0.09}, "p_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis etaAxis{"etaAxis", {180, -0.9, 0.9}, "#eta"};
  ConfigurableAxis yAxis{"yAxis", {180, -0.9, 0.9}, "y"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, 2.0 * o2::constants::math::PI}, "#phi"};
  ConfigurableAxis phiAsymmAxis{"phiAsymmAxis", {364, 0, o2::constants::math::PI}, "#phi"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // QA //
    // collisions
    registry.add("QC/collisions/hPosXY", ";x (cm);y (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    registry.add("QC/collisions/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    registry.add("QC/collisions/hNumContrib", ";number of contributors;counts", kTH1D, {{36, -0.5, 35.5}});
    registry.add("QC/collisions/hZdcCommonEnergy", ";ZNA common energy;ZNC common energy;counts", kTH2D, {{250, -5.0, 20.0}, {250, -5.0, 20.0}});
    registry.add("QC/collisions/hZdcTime", ";ZNA time (ns);ZNC time (ns);counts", kTH2D, {{200, -10.0, 10.0}, {200, -10.0, 10.0}});
    registry.add("QC/collisions/hZnaTimeVsCommonEnergy", ";ZNA common energy;ZNA time (ns);counts", kTH2D, {{250, -5.0, 20.0}, {200, -10.0, 10.0}});
    registry.add("QC/collisions/hZncTimeVsCommonEnergy", ";ZNC common energy;ZNC time (ns);counts", kTH2D, {{250, -5.0, 20.0}, {200, -10.0, 10.0}});
    // all tracks
    registry.add("QC/allTracks/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {{400, -10.0, 30.0}});
    registry.add("QC/allTracks/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {{400, -20.0, 20.0}});
    registry.add("QC/allTracks/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {1000, -5.0, 5.0}});
    // tracks passing selections
    registry.add("QC/cutTracks/hTpcNsigmaPi2D", ";TPC n#sigma(#pi_{1});TPC n#sigma(#pi_{2});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    registry.add("QC/cutTracks/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    registry.add("QC/cutTracks/hRemainingTracks", ";remaining tracks;counts", kTH1D, {{21, -0.5, 20.5}});

    // RECO HISTOS //
    // PIONS
    // no selection
    registry.add("reco/pions/unlike-sign/hPt", ";p_{T}(#pi_{1}) (GeV/#it{c});p_{T}(#pi_{2}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/unlike-sign/hEta", ";#eta(#pi_{1});#eta(#pi_{2});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/unlike-sign/hPhi", ";#phi(#pi_{1});#phi(#pi_{2});counts", kTH2D, {phiAxis, phiAxis});
    registry.add("reco/pions/like-sign/hPt", ";p_{T}(#pi_{1}) (GeV/#it{c});p_{T}(#pi_{2}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/like-sign/hEta", ";#eta(#pi_{1});#eta(#pi_{2});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/like-sign/hPhi", ";#phi(#pi_{1});#phi(#pi_{2});counts", kTH2D, {phiAxis, phiAxis});

    // RAW RHOS
    registry.add("reco/system/2pi/raw/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/raw/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/raw/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/raw/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/raw/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/raw/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/raw/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/raw/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/raw/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/raw/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/raw/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/raw/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});

    // SELECTED RHOS
    // no selection
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/no-selection/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    // 0n0n
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0n0n/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    // Xn0n
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/Xn0n/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    // 0nXn
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/0nXn/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    // XnXn
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("reco/system/2pi/cut/XnXn/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});

    // 4PI AND 6PI SYSTEM
    registry.add("reco/system/4pi/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/4pi/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/4pi/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/4pi/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/6pi/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/6pi/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/6pi/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/6pi/hY", ";y;counts", kTH1D, {yAxis});
  }

  template <typename T>
  bool collisionPassesCuts(T const& collision) // collision cuts
  {
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut)
      return false;
    if (specifyGapSide && collision.gapSide() != 2)
      return false;
    return true;
  }

  template <typename T>
  bool trackPassesCuts(T const& track) // track cuts (PID done separately)
  {
    if (!track.isPVContributor())
      return false;
    if (!track.hasITS())
      return false;
    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > tracksDcaMaxCut)
      return false;
    if (std::abs(eta(track.px(), track.py(), track.pz())) > PcEtaCut)
      return false;
    return true;
  }

  template <typename T>
  bool tracksPassPiPID(const T& cutTracks) // n-dimensional PID cut
  {
    double radius = 0.0;
    for (const auto& track : cutTracks)
      radius += std::pow(track.tpcNSigmaPi(), 2);
    return radius < std::pow(tracksTpcNSigmaPiCut, 2);
  }

  template <typename T>
  double tracksTotalCharge(const T& cutTracks) // total charge of selected tracks
  {
    double charge = 0.0;
    for (const auto& track : cutTracks)
      charge += track.sign();
    return charge;
  }

  template <typename T>
  bool systemPassCuts(const T& system) // system cuts
  {
    if (system.M() < systemMassMinCut || system.M() > systemMassMaxCut)
      return false;
    if (system.Pt() > systemPtCut)
      return false;
    if (std::abs(system.Rapidity()) > systemYCut)
      return false;
    return true;
  }

  ROOT::Math::PxPyPzMVector reconstructSystem(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracks4Vecs) // reconstruct system from 4-vectors
  {
    ROOT::Math::PxPyPzMVector system;
    for (const auto& track4Vec : cutTracks4Vecs)
      system += track4Vec;
    return system;
  }

  template <typename T>
  double getPhiRandom(const T& cutTracks) // decay phi anisotropy
  {                                       // two possible definitions of phi: randomize the tracks
    std::vector<int> indices = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // get time-based seed
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    ROOT::Math::XYVector pOne(cutTracks[indices[0]].px(), cutTracks[indices[0]].py());
    ROOT::Math::XYVector pTwo(cutTracks[indices[1]].px(), cutTracks[indices[1]].py());
    auto pPlus = pOne + pTwo;
    auto pMinus = pOne - pTwo;
    // no method for direct calculation of angle -> use dot product formula
    double cosPhi = (pPlus.Dot(pMinus)) / (std::sqrt(pPlus.Mag2()) * std::sqrt(pMinus.Mag2()));
    return std::acos(cosPhi);
  }

  template <typename T>
  double getPhiCharge(const T& cutTracks)
  { // two possible definitions of phi: charge-based assignment
    ROOT::Math::XYVector pOne, pTwo;
    if (cutTracks[0].sign() > 0) {
      pOne.SetXY(cutTracks[0].px(), cutTracks[0].py());
      pTwo.SetXY(cutTracks[1].px(), cutTracks[1].py());
    } else {
      pOne.SetXY(cutTracks[1].px(), cutTracks[1].py());
      pTwo.SetXY(cutTracks[0].px(), cutTracks[0].py());
    }
    auto pPlus = pOne + pTwo;
    auto pMinus = pOne - pTwo;
    double cosPhi = (pPlus.Dot(pMinus)) / (std::sqrt(pPlus.Mag2()) * std::sqrt(pMinus.Mag2()));
    return std::acos(cosPhi);
  }

  void processReco(FullUDSgCollision const& collision, FullUDTracks const& tracks)
  {
    // QC histograms
    registry.fill(HIST("QC/collisions/hPosXY"), collision.posX(), collision.posY());
    registry.fill(HIST("QC/collisions/hPosZ"), collision.posZ());
    registry.fill(HIST("QC/collisions/hZdcCommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    registry.fill(HIST("QC/collisions/hZdcTime"), collision.timeZNA(), collision.timeZNC());
    registry.fill(HIST("QC/collisions/hZnaTimeVsCommonEnergy"), collision.energyCommonZNA(), collision.timeZNA());
    registry.fill(HIST("QC/collisions/hZncTimeVsCommonEnergy"), collision.energyCommonZNC(), collision.timeZNC());
    registry.fill(HIST("QC/collisions/hNumContrib"), collision.numContrib());

    if (!collisionPassesCuts(collision))
      return;

    // event tagging
    bool XnXn = false, OnOn = false, XnOn = false, OnXn = false; // note: On == 0n...
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() < ZNcommonEnergyCut)
      OnOn = true;
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut &&
        collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut)
      XnXn = true;
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() < ZNcommonEnergyCut)
      XnOn = true;
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut)
      OnXn = true;
    // vectors for storing selected tracks and their 4-vectors
    std::vector<decltype(tracks.begin())> cutTracks;
    std::vector<ROOT::Math::PxPyPzMVector> cutTracks4Vecs;

    int trackCounter = 0;
    for (const auto& track : tracks) {
      registry.fill(HIST("QC/allTracks/hTpcNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("QC/allTracks/hTofNSigmaPi"), track.tofNSigmaPi());
      registry.fill(HIST("QC/allTracks/hDcaXYZ"), track.dcaZ(), track.dcaXY());

      if (!trackPassesCuts(track))
        continue;
      trackCounter++;
      cutTracks.push_back(track);
      cutTracks4Vecs.push_back(ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged)); // apriori assume pion mass
      registry.fill(HIST("QC/cutTracks/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
    }
    registry.fill(HIST("QC/cutTracks/hRemainingTracks"), trackCounter);

    if (cutTracks.size() == 2)
      registry.fill(HIST("QC/cutTracks/hTpcNsigmaPi2D"), cutTracks[0].tpcNSigmaPi(), cutTracks[1].tpcNSigmaPi());

    if (!tracksPassPiPID(cutTracks))
      return;
    // reonstruct system and calculate total charge
    auto system = reconstructSystem(cutTracks4Vecs);
    int totalCharge = tracksTotalCharge(cutTracks);
    int nTracks = cutTracks.size();

    if (nTracks == 2) {
      // fill raw histograms according to the total charge
      if (totalCharge == 0) {
        registry.fill(HIST("reco/pions/unlike-sign/hPt"), cutTracks4Vecs[0].Pt(), cutTracks4Vecs[1].Pt());
        registry.fill(HIST("reco/pions/unlike-sign/hEta"), cutTracks4Vecs[0].Eta(), cutTracks4Vecs[1].Eta());
        registry.fill(HIST("reco/pions/unlike-sign/hPhi"), cutTracks4Vecs[0].Phi() + o2::constants::math::PI, cutTracks4Vecs[1].Phi() + o2::constants::math::PI);
        registry.fill(HIST("reco/system/2pi/raw/unlike-sign/hM"), system.M());
        registry.fill(HIST("reco/system/2pi/raw/unlike-sign/hPt"), system.Pt());
        registry.fill(HIST("reco/system/2pi/raw/unlike-sign/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/2pi/raw/unlike-sign/hY"), system.Rapidity());
      } else {
        registry.fill(HIST("reco/pions/like-sign/hPt"), cutTracks4Vecs[0].Pt(), cutTracks4Vecs[1].Pt());
        registry.fill(HIST("reco/pions/like-sign/hEta"), cutTracks4Vecs[0].Eta(), cutTracks4Vecs[1].Eta());
        registry.fill(HIST("reco/pions/like-sign/hPhi"), cutTracks4Vecs[0].Phi() + o2::constants::math::PI, cutTracks4Vecs[1].Phi() + o2::constants::math::PI);
        if (totalCharge == 2) {
          registry.fill(HIST("reco/system/2pi/raw/like-sign/positive/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/raw/like-sign/positive/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/raw/like-sign/positive/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/raw/like-sign/positive/hY"), system.Rapidity());
        } else if (totalCharge == -2) {
          registry.fill(HIST("reco/system/2pi/raw/like-sign/negative/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/raw/like-sign/negative/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/raw/like-sign/negative/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/raw/like-sign/negative/hY"), system.Rapidity());
        }
      }
      // apply cuts to system
      if (!systemPassCuts(system))
        return;
      // fill histograms for system passing cuts
      switch (totalCharge) {
        case 0:
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hPt2"), system.Pt() * system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hPhiRandom"), getPhiRandom(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hPhiCharge"), getPhiCharge(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/unlike-sign/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          if (OnOn) {
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/unlike-sign/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (XnOn) {
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/unlike-sign/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (OnXn) {
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/unlike-sign/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (XnXn) {
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/unlike-sign/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          }
          break;

        case 2:
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hPt2"), system.Pt() * system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hPhiRandom"), getPhiRandom(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hPhiCharge"), getPhiCharge(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/positive/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          if (OnOn) {
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/positive/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (XnOn) {
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/positive/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (OnXn) {
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/positive/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (XnXn) {
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/positive/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          }
          break;

        case -2:
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hPt2"), system.Pt() * system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hPhiRandom"), getPhiRandom(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hPhiCharge"), getPhiCharge(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
          registry.fill(HIST("reco/system/2pi/cut/no-selection/like-sign/negative/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          if (OnOn) {
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0n0n/like-sign/negative/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (XnOn) {
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/Xn0n/like-sign/negative/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (OnXn) {
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/0nXn/like-sign/negative/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          } else if (XnXn) {
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hM"), system.M());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hPt"), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hPt2"), system.Pt() * system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hPtVsM"), system.M(), system.Pt());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hY"), system.Rapidity());
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hPhiRandom"), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hPhiCharge"), getPhiCharge(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hPhiRandomVsM"), system.M(), getPhiRandom(cutTracks));
            registry.fill(HIST("reco/system/2pi/cut/XnXn/like-sign/negative/hPhiChargeVsM"), system.M(), getPhiCharge(cutTracks));
          }
          break;

        default:
          break;
      }
    } else if (nTracks == 4 && tracksTotalCharge(cutTracks) == 0) { // 4pi system
      registry.fill(HIST("reco/system/4pi/hM"), system.M());
      registry.fill(HIST("reco/system/4pi/hPt"), system.Pt());
      registry.fill(HIST("reco/system/4pi/hPtVsM"), system.M(), system.Pt());
      registry.fill(HIST("reco/system/4pi/hY"), system.Rapidity());
    } else if (nTracks == 6 && tracksTotalCharge(cutTracks) == 0) { // 6pi system
      registry.fill(HIST("reco/system/6pi/hM"), system.M());
      registry.fill(HIST("reco/system/6pi/hPt"), system.Pt());
      registry.fill(HIST("reco/system/6pi/hPtVsM"), system.M(), system.Pt());
      registry.fill(HIST("reco/system/6pi/hY"), system.Rapidity());
    }
  }
  PROCESS_SWITCH(upcRhoAnalysis, processReco, "analyse reco tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<upcRhoAnalysis>(cfgc)};
}
