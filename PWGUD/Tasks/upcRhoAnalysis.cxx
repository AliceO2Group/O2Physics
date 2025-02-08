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
/// \file   upcRhoAnalysis.cxx

#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "random"
#include "TLorentzVector.h"

#include "Common/DataModel/PIDResponse.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUdSgCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::SGCollisions>::iterator;
using FullUdDgCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>::iterator;
using FullUdTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;

namespace o2::aod
{
namespace reco_tree
{
// misc event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(LocalBC, localBC, int);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int);
// event vertex
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
// FIT info
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFDDAmplitudeA, totalFDDAmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFDDAmplitudeC, totalFDDAmplitudeC, float);
DECLARE_SOA_COLUMN(TimeFT0A, timeFT0A, float);
DECLARE_SOA_COLUMN(TimeFT0C, timeFT0C, float);
DECLARE_SOA_COLUMN(TimeFV0A, timeFV0A, float);
DECLARE_SOA_COLUMN(TimeFDDA, timeFDDA, float);
DECLARE_SOA_COLUMN(TimeFDDC, timeFDDC, float);
// ZDC info
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
// DECLARE_SOA_COLUMN(NeutronClass, neutronClass, int);
// Rhos
DECLARE_SOA_COLUMN(TotalCharge, totalCharge, int);
// store things for reconstruction of lorentz vectors
DECLARE_SOA_COLUMN(RhoPt, rhoPt, float);
DECLARE_SOA_COLUMN(RhoY, rhoY, float);
DECLARE_SOA_COLUMN(RhoPhi, rhoPhi, float);
DECLARE_SOA_COLUMN(RhoM, rhoM, float);
// other stuff
DECLARE_SOA_COLUMN(RhoPhiRandom, rhoPhiRandom, float);
DECLARE_SOA_COLUMN(RhoPhiCharge, rhoPhiCharge, float);
// Pion tracks
DECLARE_SOA_COLUMN(TrackSign, trackSign, std::vector<int>);
// for lorentz vector reconstruction
DECLARE_SOA_COLUMN(TrackPt, trackPt, std::vector<float>);
DECLARE_SOA_COLUMN(TrackEta, trackEta, std::vector<float>);
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, std::vector<float>);
DECLARE_SOA_COLUMN(TrackM, trackM, std::vector<float>);
// other stuff
DECLARE_SOA_COLUMN(TrackPiPID, trackPiPID, std::vector<float>);
DECLARE_SOA_COLUMN(TrackElPID, trackElPID, std::vector<float>);
DECLARE_SOA_COLUMN(TrackDcaXY, trackDcaXY, std::vector<float>);
DECLARE_SOA_COLUMN(TrackDcaZ, trackDcaZ, std::vector<float>);
DECLARE_SOA_COLUMN(TrackTpcSignal, trackTpcSignal, std::vector<float>);
} // namespace recoTree
DECLARE_SOA_TABLE(RecoTree, "AOD", "RECOTREE",
                  reco_tree::RunNumber, reco_tree::LocalBC, reco_tree::NumContrib,
                  reco_tree::PosX, reco_tree::PosY, reco_tree::PosZ, reco_tree::TotalFT0AmplitudeA, reco_tree::TotalFT0AmplitudeC, reco_tree::TotalFV0AmplitudeA, reco_tree::TotalFDDAmplitudeA, reco_tree::TotalFDDAmplitudeC,
                  reco_tree::TimeFT0A, reco_tree::TimeFT0C, reco_tree::TimeFV0A, reco_tree::TimeFDDA, reco_tree::TimeFDDC,
                  reco_tree::EnergyCommonZNA, reco_tree::EnergyCommonZNC, reco_tree::TimeZNA, reco_tree::TimeZNC, /* reco_tree::NeutronClass, */
                  reco_tree::TotalCharge, reco_tree::RhoPt, reco_tree::RhoY, reco_tree::RhoPhi, reco_tree::RhoM, reco_tree::RhoPhiRandom, reco_tree::RhoPhiCharge,
                  reco_tree::TrackSign, reco_tree::TrackPt, reco_tree::TrackEta, reco_tree::TrackPhi, reco_tree::TrackM, reco_tree::TrackPiPID, reco_tree::TrackElPID, reco_tree::TrackDcaXY, reco_tree::TrackDcaZ, reco_tree::TrackTpcSignal);

namespace mc_tree
{
// misc event info
DECLARE_SOA_COLUMN(LocalBc, localBc, int);
// event vertex
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
// rho lorentz vector reconstruction
DECLARE_SOA_COLUMN(RhoPt, rhoPt, float);
DECLARE_SOA_COLUMN(RhoY, rhoY, float);
DECLARE_SOA_COLUMN(RhoPhi, rhoPhi, float);
DECLARE_SOA_COLUMN(RhoM, rhoM, float);
// other stuff
DECLARE_SOA_COLUMN(RhoPhiRandom, rhoPhiRandom, float);
DECLARE_SOA_COLUMN(RhoPhiCharge, rhoPhiCharge, float);
// Pion tracks
DECLARE_SOA_COLUMN(TrackSign, trackSign, std::vector<int>);
// for lorentz vector reconstruction
DECLARE_SOA_COLUMN(TrackPt, trackPt, std::vector<float>);
DECLARE_SOA_COLUMN(TrackEta, trackEta, std::vector<float>);
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, std::vector<float>);
DECLARE_SOA_COLUMN(TrackM, trackM, std::vector<float>);
} // namespace mcTree
DECLARE_SOA_TABLE(McTree, "AOD", "MCTREE",
                  mc_tree::LocalBc,
                  mc_tree::PosX, mc_tree::PosY, mc_tree::PosZ, mc_tree::RhoPt, mc_tree::RhoY, mc_tree::RhoPhi, mc_tree::RhoM, mc_tree::RhoPhiRandom, mc_tree::RhoPhiCharge,
                  mc_tree::TrackSign, mc_tree::TrackPt, mc_tree::TrackEta, mc_tree::TrackPhi, mc_tree::TrackM);
} // namespace o2::aod

struct upcRhoAnalysis {
  Produces<o2::aod::RecoTree> recoTree;
  Produces<o2::aod::McTree> mcTree;

  float pcEtaCut = 0.9; // physics coordination recommendation
  Configurable<bool> requireTof{"requireTof", false, "require TOF signal"};

  Configurable<bool> doPions{"doPions", true, "do analysis with PIDed pions"};
  Configurable<bool> doElectrons{"doElectrons", true, "do analysis with PIDed electrons"};
  Configurable<bool> doKaons{"doKaons", true, "do analysis with PIDed kaons"};

  Configurable<float> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "max Z position cut on collisions"};
  Configurable<int> collisionsNumContribsMaxCut{"collisionsNumContribsMaxCut", 7, "max number of contributors cut on collisions"};
  Configurable<float> znCommonEnergyCut{"znCommonEnergyCut", 0.0, "ZN common energy cut"};
  Configurable<float> znTimeCut{"znTimeCut", 2.0, "ZN time cut"};

  Configurable<float> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "TPC nSigma pion cut"};
  Configurable<float> tracksTpcNSigmaElCut{"tracksTpcNSigmaElCut", 3.0, "TPC nSigma electron cut"};
  Configurable<float> tracksTpcNSigmaKaCut{"tracksTpcNSigmaKaCut", 3.0, "TPC nSigma kaon cut"};
  Configurable<float> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "max DCA cut on tracks"};
  Configurable<int> tracksMinItsNClsCut{"tracksMinItsNClsCut", 6, "min ITS clusters cut"};
  Configurable<float> tracksMaxItsChi2NClCut{"tracksMaxItsChi2NClCut", 3.0, "max ITS chi2/Ncls cut"};
  Configurable<int> tracksMinTpcNClsCut{"tracksMinTpcNClsCut", 120, "min TPC clusters cut"};
  Configurable<int> tracksMinTpcNClsCrossedRowsCut{"tracksMinTpcNClsCrossedRowsCut", 140, "min TPC crossed rows cut"};
  Configurable<float> tracksMinTpcChi2NClCut{"tracksMinTpcChi2NClCut", 1.0, "min TPC chi2/Ncls cut"};
  Configurable<float> tracksMaxTpcChi2NClCut{"tracksMaxTpcChi2NClCut", 1.8, "max TPC chi2/Ncls cut"};
  Configurable<float> tracksMinTpcNClsCrossedOverFindableCut{"tracksMinTpcNClsCrossedOverFindableCut", 1.05, "min TPC crossed rows / findable clusters cut"};
  Configurable<float> tracksMinPtCut{"tracksMinPtCut", 0.2, "min pT cut on tracks"};

  Configurable<float> systemMassMinCut{"systemMassMinCut", 0.4, "min M cut for reco system"};
  Configurable<float> systemMassMaxCut{"systemMassMaxCut", 1.2, "max M cut for reco system"};
  Configurable<float> systemPtCut{"systemPtMaxCut", 0.1, "max pT cut for reco system"};
  Configurable<float> systemYCut{"systemYCut", 0.9, "rapiditiy cut for reco system"};

  ConfigurableAxis mAxis{"mAxis", {1000, 0.0, 10.0}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis mCutAxis{"mCutAxis", {160, 0.4, 1.2}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {1000, 0.0, 10.0}, "p_{T} (GeV/#it{c})"};
  ConfigurableAxis ptCutAxis{"ptCutAxis", {100, 0.0, 0.1}, "p_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {100, 0.0, 0.01}, "p_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis etaAxis{"etaAxis", {800, -4.0, 4.0}, "#eta"};
  ConfigurableAxis etaCutAxis{"etaCutAxis", {180, -0.9, 0.9}, "#eta"};
  ConfigurableAxis yAxis{"yAxis", {400, -4.0, 4.0}, "y"};
  ConfigurableAxis yCutAxis{"yCutAxis", {180, -0.9, 0.9}, "y"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, o2::constants::math::TwoPI}, "#phi"};
  ConfigurableAxis phiAsymmAxis{"phiAsymmAxis", {182, -o2::constants::math::PI, o2::constants::math::PI}, "#phi"};
  ConfigurableAxis momentumFromPhiAxis{"momentumFromPhiAxis", {400, -0.1, 0.1}, "p (GeV/#it{c})"};
  ConfigurableAxis znCommonEnergyAxis{"znCommonEnergyAxis", {250, -5.0, 20.0}, "ZN common energy (TeV)"};
  ConfigurableAxis znTimeAxis{"znTimeAxis", {200, -10.0, 10.0}, "ZN time (ns)"};
  
  HistogramRegistry rQC{"rQC", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rTracks{"rTracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rSystem{"rSystem", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rMC{"rMC", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // QA //
    // collisions
    rQC.add("QC/collisions/all/hPosXY", ";x (cm);y (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    rQC.add("QC/collisions/all/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    rQC.add("QC/collisions/all/hNumContrib", ";number of contributors;counts", kTH1D, {{36, -0.5, 35.5}});
    rQC.add("QC/collisions/all/hZdcCommonEnergy", ";ZNA common energy (TeV);ZNC common energy (TeV);counts", kTH2D, {znCommonEnergyAxis, znCommonEnergyAxis});
    rQC.add("QC/collisions/all/hZdcTime", ";ZNA time (ns);ZNC time (ns);counts", kTH2D, {znTimeAxis, znTimeAxis});
    rQC.add("QC/collisions/all/hTotalFT0AmplitudeA", ";FT0A amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    rQC.add("QC/collisions/all/hTotalFT0AmplitudeC", ";FT0C amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    rQC.add("QC/collisions/all/hTotalFV0AmplitudeA", ";FV0A amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    rQC.add("QC/collisions/all/hTotalFDDAmplitudeA", ";FDDA amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    rQC.add("QC/collisions/all/hTotalFDDAmplitudeC", ";FDDC amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    rQC.add("QC/collisions/all/hTimeFT0A", ";FT0A time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    rQC.add("QC/collisions/all/hTimeFT0C", ";FT0C time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    rQC.add("QC/collisions/all/hTimeFV0A", ";FV0A time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    rQC.add("QC/collisions/all/hTimeFDDA", ";FDDA time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    rQC.add("QC/collisions/all/hTimeFDDC", ";FDDC time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    // events with selected rho candidates
    rQC.addClone("QC/collisions/all/", "QC/collisions/selected/"); // clone "all" histograms as "selected"
    rQC.add("QC/collisions/selected/hNumContribVsSystemPt", ";p_{T} (GeV/#it{c});number of contributors;counts", kTH2D, {ptAxis,{36, -0.5, 35.5}});
    rQC.add("QC/collisions/selected/hNumContribVsSystemY", ";y;number of contributors;counts", kTH2D, {yAxis,{36, -0.5, 35.5}});
    rQC.add("QC/collisions/selected/hNumContribVsSystemM", ";m (GeV/#it{c}^{2});number of contributors;counts", kTH2D, {mAxis,{36, -0.5, 35.5}});

    // tracks
    rQC.add("QC/tracks/raw/hTpcNSigmaPi", ";TPC n#sigma(#pi);counts", kTH1D, {{400, -10.0, 30.0}});
    rQC.add("QC/tracks/raw/hTpcNSigmaEl", ";TPC n#sigma(e);counts", kTH1D, {{400, -10.0, 30.0}});
    rQC.add("QC/tracks/raw/hTpcNSigmaKa", ";TPC n#sigma(K);counts", kTH1D, {{400, -10.0, 30.0}});
    rQC.add("QC/tracks/raw/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {1000, -5.0, 5.0}});
    rQC.add("QC/tracks/raw/hItsNCls", ";ITS N_{cls};counts", kTH1D, {{11, -0.5, 10.5}});
    rQC.add("QC/tracks/raw/hItsChi2NCl", ";ITS #chi^{2}/N_{cls};counts", kTH1D, {{1000, 0.0, 100.0}});
    rQC.add("QC/tracks/raw/hTpcChi2NCl", ";TPC #chi^{2}/N_{cls};counts", kTH1D, {{1000, 0.0, 100.0}});
    rQC.add("QC/tracks/raw/hTpcNCls", ";TPC N_{cls} found;counts", kTH1D, {{200, 0.0, 200.0}});
    rQC.add("QC/tracks/raw/hTpcNClsCrossedRows", ";TPC crossed rows;counts", kTH1D, {{200, 0.0, 200.0}});
    rQC.add("QC/tracks/raw/hTpcNClsCrossedRowsOverNClsFindable", ";TPC crossed rows/findable N_{cls};counts", kTH1D, {{100, 0.0, 10.0}});
    rQC.add("QC/tracks/raw/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rQC.add("QC/tracks/raw/hEta", ";y;counts", kTH1D, {etaAxis});
    rQC.add("QC/tracks/raw/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    rQC.add("QC/tracks/raw/hTpcSignalVsP", ";p (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    rQC.add("QC/tracks/raw/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    // tracks passing selections
    rQC.addClone("QC/tracks/raw/", "QC/tracks/cut/"); // clone "raw" histograms as "cut"
    rQC.add("QC/tracks/cut/hRemainingTracks", ";remaining tracks;counts", kTH1D, {{21, -0.5, 20.5}});
    rQC.add("QC/tracks/cut/hTpcNSigmaPi2D", ";TPC n#sigma(#pi_{leading});TPC n#sigma(#pi_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    rQC.add("QC/tracks/cut/hTpcNSigmaEl2D", ";TPC n#sigma(e_{leading});TPC n#sigma(e_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    rQC.add("QC/tracks/cut/hTpcNSigmaKa2D", ";TPC n#sigma(K_{leading});TPC n#sigma(K_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    // selection counter
    std::vector<std::string> selectionCounterLabels = {"all tracks", "PV contributor", "ITS hit", "ITS N_{clusters}", "ITS #chi^{2}/N_{clusters}", "TPC hit", "TPC N_{clusters} found", "TPC #chi^{2}/N_{clusters}", "TPC crossed rows",
                                                       "TPC crossed rows/N_{clusters}",
                                                       "TOF requirement",
                                                       "p_{T}", "DCA", "#eta", "exactly 2 tracks", "PID"};
    rQC.add("QC/tracks/hSelectionCounter", ";;tracks passing selections", kTH1D, {{static_cast<int>(selectionCounterLabels.size()), -0.5, static_cast<float>(selectionCounterLabels.size()) - 0.5}});
    for (int i = 0; i < static_cast<int>(selectionCounterLabels.size()); ++i)
      rQC.get<TH1>(HIST("QC/tracks/hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, selectionCounterLabels[i].c_str());
    rQC.add("QC/tracks/hTofHitCheck", ";leading track TOF hit;subleading track TOF hit;counts", kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});

    // TRACKS (2D)
    // pions
    if (doPions) {
      rTracks.add("tracks/pions/no-selection/unlike-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
      rTracks.add("tracks/pions/no-selection/unlike-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaCutAxis, etaCutAxis});
      rTracks.add("tracks/pions/no-selection/unlike-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
      rTracks.addClone("tracks/pions/no-selection/unlike-sign/", "tracks/pions/no-selection/like-sign/");
      rTracks.addClone("tracks/pions/no-selection/", "tracks/pions/selected/");
    }
    // electrons
    if (doElectrons) {
      rTracks.add("tracks/electrons/no-selection/unlike-sign/hPt", ";p_{T}(e_{leading}) (GeV/#it{c});p_{T}(e_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
      rTracks.add("tracks/electrons/no-selection/unlike-sign/hEta", ";#eta(e_{leading});#eta(e_{subleading});counts", kTH2D, {etaCutAxis, etaCutAxis});
      rTracks.add("tracks/electrons/no-selection/unlike-sign/hPhi", ";#phi(e_{leading});#phi(e_{subleading});counts", kTH2D, {phiAxis, phiAxis});
      rTracks.addClone("tracks/electrons/no-selection/unlike-sign/", "tracks/electrons/no-selection/like-sign/");
      rTracks.addClone("tracks/electrons/no-selection/", "tracks/electrons/selected/");
    }
    // kaons
    if (doKaons) {
      rTracks.add("tracks/kaons/no-selection/unlike-sign/hPt", ";p_{T}(K_{leading}) (GeV/#it{c});p_{T}(K_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
      rTracks.add("tracks/kaons/no-selection/unlike-sign/hEta", ";#eta(K_{leading});#eta(K_{subleading});counts", kTH2D, {etaCutAxis, etaCutAxis});
      rTracks.add("tracks/kaons/no-selection/unlike-sign/hPhi", ";#phi(K_{leading});#phi(K_{subleading});counts", kTH2D, {phiAxis, phiAxis});
      rTracks.addClone("tracks/kaons/no-selection/unlike-sign/", "tracks/kaons/no-selection/like-sign/");
      rTracks.addClone("tracks/kaons/no-selection/", "tracks/kaons/selected/");
    }

    // SYSTEM
    rSystem.add("system/raw/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    rSystem.add("system/raw/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rSystem.add("system/raw/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    rSystem.add("system/raw/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    rSystem.add("system/raw/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    // clones for like-sign
    rSystem.addClone("system/raw/unlike-sign/", "system/raw/like-sign/positive/");
    rSystem.addClone("system/raw/unlike-sign/", "system/raw/like-sign/negative/");
    // selected rhos
    rSystem.add("system/cut/no-selection/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    rSystem.add("system/cut/no-selection/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    // clones for like-sign
    rSystem.addClone("system/cut/no-selection/unlike-sign/", "system/cut/no-selection/like-sign/positive/");
    rSystem.addClone("system/cut/no-selection/unlike-sign/", "system/cut/no-selection/like-sign/negative/");
    // clones for neutron classes
    rSystem.addClone("system/cut/no-selection/", "system/cut/0n0n/");
    rSystem.addClone("system/cut/no-selection/", "system/cut/Xn0n/");
    rSystem.addClone("system/cut/no-selection/", "system/cut/0nXn/");
    rSystem.addClone("system/cut/no-selection/", "system/cut/XnXn/");
    
    // MC
    // collisions
    rMC.add("MC/collisions/hPosXY", ";x (cm);y (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    rMC.add("MC/collisions/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    rMC.add("MC/collisions/hNPions", ";number of pions;counts", kTH1D, {{11, -0.5, 10.5}});
    rMC.add("MC/collisions/hNumOfCollisionRecos", ";number of collision reconstructions;counts", kTH1D, {{11, -0.5, 10.5}});
    // tracks
    rMC.add("MC/tracks/all/hPdgCode", ";pdg code;counts", kTH1D, {{2001, -1000.5, 1000.5}});
    rMC.add("MC/tracks/all/hProducedByGenerator", ";produced by generator;counts", kTH1D, {{2, -0.5, 1.5}});
    rMC.add("MC/tracks/all/hIsPhysicalPrimary", ";is physical primary;counts", kTH1D, {{2, -0.5, 1.5}});
    rMC.add("MC/tracks/all/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rMC.add("MC/tracks/all/hEta", ";#eta;counts", kTH1D, {etaAxis});
    rMC.add("MC/tracks/all/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    rMC.add("MC/tracks/pions/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    rMC.add("MC/tracks/pions/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaAxis, etaAxis});
    rMC.add("MC/tracks/pions/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    // system
    rMC.add("MC/system/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    rMC.add("MC/system/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rMC.add("MC/system/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    rMC.add("MC/system/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    rMC.add("MC/system/hY", ";y;counts", kTH1D, {yAxis});
    rMC.add("MC/system/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    rMC.add("MC/system/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    rMC.add("MC/system/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
  }

  template <int cuts, typename C>
  void fillCollisionQcHistos(const C& collision) // fills collision QC histograms before/after cuts
  {
    static constexpr std::string_view Mode[2] = {"all", "selected"};

    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hPosXY"), collision.posX(), collision.posY());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hPosZ"), collision.posZ());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hZdcCommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hZdcTime"), collision.timeZNA(), collision.timeZNC());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hNumContrib"), collision.numContrib());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTotalFT0AmplitudeA"), collision.totalFT0AmplitudeA());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTotalFT0AmplitudeC"), collision.totalFT0AmplitudeC());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTotalFV0AmplitudeA"), collision.totalFV0AmplitudeA());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTotalFDDAmplitudeA"), collision.totalFDDAmplitudeA());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTotalFDDAmplitudeC"), collision.totalFDDAmplitudeC());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTimeFT0A"), collision.timeFT0A());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTimeFT0C"), collision.timeFT0C());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTimeFV0A"), collision.timeFV0A());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTimeFDDA"), collision.timeFDDA());
    rQC.fill(HIST("QC/collisions/") + HIST(Mode[cuts]) + HIST("/hTimeFDDC"), collision.timeFDDC());
  }

  template <typename C>
  bool collisionPassesCuts(const C& collision) // collision cuts
  {
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut)
      return false;
    // if (collision.numContrib() > collisionsNumContribsMaxCut)
    //   return false;
    return true;
  }

  template <typename T>
  bool trackPassesCuts(const T& track) // track cuts (PID done separately)
  {
    if (!track.isPVContributor())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 1);

    if (!track.hasITS())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 2);

    if (track.itsNCls() < tracksMinItsNClsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 3);

    if (track.itsChi2NCl() > tracksMaxItsChi2NClCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 4);

    if (!track.hasTPC())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 5);

    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < tracksMinTpcNClsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 6);

    if (track.tpcChi2NCl() > tracksMaxTpcChi2NClCut || track.tpcChi2NCl() < tracksMinTpcChi2NClCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 7);

    if (track.tpcNClsCrossedRows() < tracksMinTpcNClsCrossedRowsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 8);

    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < tracksMinTpcNClsCrossedOverFindableCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 9);

    if (requireTof && !track.hasTOF())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 10);

    if (track.pt() < tracksMinPtCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 11);

    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > (0.0105 + 0.0350 / std::pow(track.pt(), 1.01)))
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 12);

    if (std::abs(eta(track.px(), track.py(), track.pz())) > pcEtaCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 13);
    // if all selections passed
    return true;
  }

  template <typename T>
  int mostProbablePID(const T& cutTracks) {
    // calculate NSigma "radii"
    float squareRadiusPi = 0.0, squareRadiusEl = 0.0, squareRadiusKa = 0.0;  
    for (const auto& track : cutTracks) {
      if (doPions)
        squareRadiusPi += std::pow(track.tpcNSigmaPi(), 2);
      if (doElectrons)
        squareRadiusEl += std::pow(track.tpcNSigmaEl(), 2);
      if (doKaons)
        squareRadiusKa += std::pow(track.tpcNSigmaKa(), 2);
    }
    // select most probable PID based on lowest NSigma radius
    if (doPions && (!doElectrons || squareRadiusPi < squareRadiusEl) && (!doKaons || squareRadiusPi < squareRadiusKa) && squareRadiusPi < std::pow(tracksTpcNSigmaPiCut, 2))
      return 0;
    else if (doElectrons && (!doPions || squareRadiusEl < squareRadiusPi) && (!doKaons || squareRadiusEl < squareRadiusKa) && squareRadiusEl < std::pow(tracksTpcNSigmaElCut, 2))
      return 1;
    else if (doKaons && (!doPions || squareRadiusKa < squareRadiusPi) && (!doElectrons || squareRadiusKa < squareRadiusEl) && squareRadiusKa < std::pow(tracksTpcNSigmaKaCut, 2))
      return 2;  
    else
      return -1;
  }

  template <typename T>
  bool tracksPassPiPID(const T& cutTracks) // n-dimensional pion PID cut
  {
    float radius = 0.0;
    for (const auto& track : cutTracks)
      radius += std::pow(track.tpcNSigmaPi(), 2);
    return radius < std::pow(tracksTpcNSigmaPiCut, 2);
  }

  template <typename T>
  bool tracksPassElPID(const T& cutTracks) // n-dimensional electron PID cut
  {
    float radius = 0.0;
    for (const auto& track : cutTracks)
      radius += std::pow(track.tpcNSigmaEl(), 2);
    return radius < std::pow(tracksTpcNSigmaElCut, 2);
  }

  template <typename T>
  int tracksTotalCharge(const T& cutTracks) // total charge of selected tracks
  {
    int charge = 0;
    for (const auto& track : cutTracks)
      charge += track.sign();
    return charge;
  }

  template <typename T>
  int tracksTotalChargeMC(const T& cutTracks) // total charge of selected tracks
  {
    int charge = 0;
    for (const auto& track : cutTracks)
      charge += track.pdgCode();
    return charge;
  }

  bool systemPassesCuts(const TLorentzVector& system) // system cuts
  {
    if (system.M() < systemMassMinCut || system.M() > systemMassMaxCut)
      return false;
    if (system.Pt() > systemPtCut)
      return false;
    if (std::abs(system.Rapidity()) > systemYCut)
      return false;
    return true;
  }

  TLorentzVector reconstructSystem(const std::vector<TLorentzVector>& cutTracks4Vecs) // reconstruct system from 4-vectors
  {
    TLorentzVector system;
    for (const auto& track4Vec : cutTracks4Vecs)
      system += track4Vec;
    return system;
  }

  float getPhiRandom(const std::vector<TLorentzVector>& cutTracks4Vecs) // decay phi anisotropy
  {                                                                      // two possible definitions of phi: randomize the tracks
    std::vector<int> indices = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // get time-based seed
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    TLorentzVector pOne = cutTracks4Vecs[indices[0]];
    TLorentzVector pTwo = cutTracks4Vecs[indices[1]];
    TLorentzVector pPlus = pOne + pTwo;
    TLorentzVector pMinus = pOne - pTwo;
    return pPlus.DeltaPhi(pMinus);
  }

  template <typename T>
  float getPhiCharge(const T& cutTracks, const std::vector<TLorentzVector>& cutTracks4Vecs)
  { // two possible definitions of phi: charge-based assignment
    TLorentzVector pOne, pTwo;
    pOne = (cutTracks[0].sign() > 0) ? cutTracks4Vecs[0] : cutTracks4Vecs[1];
    pTwo = (cutTracks[0].sign() > 0) ? cutTracks4Vecs[1] : cutTracks4Vecs[0];
    TLorentzVector pPlus = pOne + pTwo;
    TLorentzVector pMinus = pOne - pTwo;
    return pPlus.DeltaPhi(pMinus);
  }

  template <typename T>
  float getPhiChargeMC(const T& cutTracks, const std::vector<TLorentzVector>& cutTracks4Vecs)
  { // the same as for data but using pdg code instead of charge
    TLorentzVector pOne, pTwo;
    pOne = (cutTracks[0].pdgCode() > 0) ? cutTracks4Vecs[0] : cutTracks4Vecs[1];
    pTwo = (cutTracks[0].pdgCode() > 0) ? cutTracks4Vecs[1] : cutTracks4Vecs[0];
    TLorentzVector pPlus = pOne + pTwo;
    TLorentzVector pMinus = pOne - pTwo;
    return pPlus.DeltaPhi(pMinus);
  }

  template <typename C, typename T>
  void processReco(C const& collision, T const& tracks)
  {
    fillCollisionQcHistos<0>(collision); // fill QC histograms before cuts
    if (!collisionPassesCuts(collision))
      return;

    // event tagging
    bool xnxn = false, onon = false, xnon = false, onxn = false; // note: On == 0n...
    int neutronClass = -1;
    if (collision.energyCommonZNA() < znCommonEnergyCut && collision.energyCommonZNC() < znCommonEnergyCut) {
      onon = true;
      neutronClass = 0;
    }
    if (collision.energyCommonZNA() > znCommonEnergyCut && std::abs(collision.timeZNA()) < znTimeCut && collision.energyCommonZNC() < znCommonEnergyCut) {
      xnon = true;
      neutronClass = 1;
    }
    if (collision.energyCommonZNA() < znCommonEnergyCut && collision.energyCommonZNC() > znCommonEnergyCut && std::abs(collision.timeZNC()) < znTimeCut) {
      onxn = true;
      neutronClass = 2;
    }
    if (collision.energyCommonZNA() > znCommonEnergyCut && std::abs(collision.timeZNA()) < znTimeCut &&
        collision.energyCommonZNC() > znCommonEnergyCut && std::abs(collision.timeZNC()) < znTimeCut) {
      xnxn = true;
      neutronClass = 3;
    }
    
    // vectors for storing selected tracks and their 4-vectors
    std::vector<decltype(tracks.begin())> cutTracks;
    std::vector<TLorentzVector> cutTracks4Vecs;

    for (const auto& track : tracks) {
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 0);
      // float p = momentum(track.px(), track.py(), track.pz());
      rQC.fill(HIST("QC/tracks/raw/hPt"), track.pt());
      rQC.fill(HIST("QC/tracks/raw/hEta"), eta(track.px(), track.py(), track.pz()));
      rQC.fill(HIST("QC/tracks/raw/hPhi"), phi(track.px(), track.py()));
      rQC.fill(HIST("QC/tracks/raw/hTpcNSigmaPi"), track.tpcNSigmaPi());
      rQC.fill(HIST("QC/tracks/raw/hTpcNSigmaEl"), track.tpcNSigmaEl());
      rQC.fill(HIST("QC/tracks/raw/hTpcNSigmaKa"), track.tpcNSigmaKa());
      rQC.fill(HIST("QC/tracks/raw/hDcaXYZ"), track.dcaZ(), track.dcaXY());
      rQC.fill(HIST("QC/tracks/raw/hItsNCls"), track.itsNCls());
      rQC.fill(HIST("QC/tracks/raw/hItsChi2NCl"), track.itsChi2NCl());
      rQC.fill(HIST("QC/tracks/raw/hTpcChi2NCl"), track.tpcChi2NCl());
      rQC.fill(HIST("QC/tracks/raw/hTpcNCls"), (track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()));
      rQC.fill(HIST("QC/tracks/raw/hTpcNClsCrossedRows"), track.tpcNClsCrossedRows());
      rQC.fill(HIST("QC/tracks/raw/hTpcNClsCrossedRowsOverNClsFindable"), (static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())));

      if (!trackPassesCuts(track))
        continue;

      cutTracks.push_back(track);
      TLorentzVector track4Vec;
      track4Vec.SetXYZM(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged); // apriori assume pion mass
      cutTracks4Vecs.push_back(track4Vec);
      rQC.fill(HIST("QC/tracks/cut/hTpcNSigmaPi"), track.tpcNSigmaPi());
      rQC.fill(HIST("QC/tracks/cut/hTpcNSigmaEl"), track.tpcNSigmaEl());
      rQC.fill(HIST("QC/tracks/cut/hTpcNSigmaKa"), track.tpcNSigmaKa());
      rQC.fill(HIST("QC/tracks/cut/hTpcSignalVsP"), momentum(track.px(), track.py(), track.pz()), track.tpcSignal());
      rQC.fill(HIST("QC/tracks/cut/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
      rQC.fill(HIST("QC/tracks/cut/hDcaXYZ"), track.dcaZ(), track.dcaXY());
    }
    rQC.fill(HIST("QC/tracks/cut/hRemainingTracks"), cutTracks.size());
    if (cutTracks.size() != cutTracks4Vecs.size()) { // sanity check
      LOG(error);
      return;
    }

    // further consider only two pion systems
    if (cutTracks.size() != 2)
      return;
    for (int i = 0; i < static_cast<int>(cutTracks.size()); i++)
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 14);
    rQC.fill(HIST("QC/tracks/cut/hTpcNSigmaPi2D"), cutTracks[0].tpcNSigmaPi(), cutTracks[1].tpcNSigmaPi());
    rQC.fill(HIST("QC/tracks/cut/hTpcNSigmaEl2D"), cutTracks[0].tpcNSigmaEl(), cutTracks[1].tpcNSigmaEl());
    rQC.fill(HIST("QC/tracks/cut/hTpcNSigmaKa2D"), cutTracks[0].tpcNSigmaKa(), cutTracks[1].tpcNSigmaKa());

    if (!tracksPassPiPID(cutTracks))
      return;
    for (int i = 0; i < static_cast<int>(cutTracks.size()); i++)
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 15);

    // reonstruct system and calculate total charge, save commonly used values into variables
    TLorentzVector system = reconstructSystem(cutTracks4Vecs);
    int totalCharge = tracksTotalCharge(cutTracks);
    float mass = system.M();
    float pT = system.Pt();
    float pTsquare = pT * pT;
    float rapidity = system.Rapidity();
    float systemPhi = system.Phi() + o2::constants::math::PI;
    float phiRandom = getPhiRandom(cutTracks4Vecs);
    float phiCharge = getPhiCharge(cutTracks, cutTracks4Vecs);

    // differentiate leading- and subleading-momentum tracks
    auto leadingMomentumTrack = momentum(cutTracks[0].px(), cutTracks[0].py(), cutTracks[0].pz()) > momentum(cutTracks[1].px(), cutTracks[1].py(), cutTracks[1].pz()) ? cutTracks[0] : cutTracks[1];
    auto subleadingMomentumTrack = (leadingMomentumTrack == cutTracks[0]) ? cutTracks[1] : cutTracks[0];
    
    float leadingPt = leadingMomentumTrack.pt();
    float subleadingPt = subleadingMomentumTrack.pt();
    float leadingEta = eta(leadingMomentumTrack.px(), leadingMomentumTrack.py(), leadingMomentumTrack.pz());
    float subleadingEta = eta(subleadingMomentumTrack.px(), subleadingMomentumTrack.py(), subleadingMomentumTrack.pz());
    float leadingPhi = phi(leadingMomentumTrack.px(), leadingMomentumTrack.py());
    float subleadingPhi = phi(subleadingMomentumTrack.px(), subleadingMomentumTrack.py());

    // fill TOF hit checker
    rQC.fill(HIST("QC/tracks/hTofHitCheck"), leadingMomentumTrack.hasTOF(), subleadingMomentumTrack.hasTOF());

    // fill recoTree
    int localBc = collision.globalBC() % o2::constants::lhc::LHCMaxBunches;
    std::vector<int> trackSigns = {leadingMomentumTrack.sign(), subleadingMomentumTrack.sign()};
    std::vector<float> trackPts = {leadingPt, subleadingPt};
    std::vector<float> trackEtas = {leadingEta, subleadingEta};
    std::vector<float> trackPhis = {leadingPhi, subleadingPhi};
    std::vector<float> trackMs = {o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged};
    std::vector<float> trackPiPIDs = {leadingMomentumTrack.tpcNSigmaPi(), subleadingMomentumTrack.tpcNSigmaPi()};
    std::vector<float> trackElPIDs = {leadingMomentumTrack.tpcNSigmaEl(), subleadingMomentumTrack.tpcNSigmaEl()};
    std::vector<float> trackDcaXYs = {leadingMomentumTrack.dcaXY(), subleadingMomentumTrack.dcaXY()};
    std::vector<float> trackDcaZs = {leadingMomentumTrack.dcaZ(), subleadingMomentumTrack.dcaZ()};
    std::vector<float> trackTpcSignals = {leadingMomentumTrack.tpcSignal(), subleadingMomentumTrack.tpcSignal()};
    recoTree(collision.runNumber(), localBc, collision.numContrib(),
         collision.posX(), collision.posY(), collision.posZ(),
         collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
         collision.timeFT0A(), collision.timeFT0C(), collision.timeFV0A(), collision.timeFDDA(), collision.timeFDDC(),
         collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC(), /* neutronClass, */
         totalCharge, pT, rapidity, system.Phi(), mass, phiRandom, phiCharge,
         trackSigns, trackPts, trackEtas, trackPhis, trackMs, trackPiPIDs, trackElPIDs, trackDcaXYs, trackDcaZs, trackTpcSignals);

    // fill raw histograms according to the total charge
    switch (totalCharge) {
      case 0:
        rTracks.fill(HIST("tracks/pions/no-selection/unlike-sign/hPt"), leadingPt, subleadingPt);
        rTracks.fill(HIST("tracks/pions/no-selection/unlike-sign/hEta"), leadingEta, subleadingEta);
        rTracks.fill(HIST("tracks/pions/no-selection/unlike-sign/hPhi"), leadingPhi, subleadingPhi);
        rSystem.fill(HIST("system/raw/unlike-sign/hM"), mass);
        rSystem.fill(HIST("system/raw/unlike-sign/hPt"), pT);
        rSystem.fill(HIST("system/raw/unlike-sign/hPtVsM"), mass, pT);
        rSystem.fill(HIST("system/raw/unlike-sign/hY"), rapidity);
        rSystem.fill(HIST("system/raw/unlike-sign/hPhi"), systemPhi);
        rQC.fill(HIST("QC/collisions/selected/hNumContribVsSystemPt"), pT, collision.numContrib());
        rQC.fill(HIST("QC/collisions/selected/hNumContribVsSystemY"), rapidity, collision.numContrib());
        rQC.fill(HIST("QC/collisions/selected/hNumContribVsSystemM"), mass, collision.numContrib());
        break;

      case 2:
        rTracks.fill(HIST("tracks/pions/no-selection/like-sign/hPt"), leadingPt, subleadingPt);
        rTracks.fill(HIST("tracks/pions/no-selection/like-sign/hEta"), leadingEta, subleadingEta);
        rTracks.fill(HIST("tracks/pions/no-selection/like-sign/hPhi"), leadingPhi, subleadingPhi);
        rSystem.fill(HIST("system/raw/like-sign/positive/hM"), mass);
        rSystem.fill(HIST("system/raw/like-sign/positive/hPt"), pT);
        rSystem.fill(HIST("system/raw/like-sign/positive/hPtVsM"), mass, pT);
        rSystem.fill(HIST("system/raw/like-sign/positive/hY"), rapidity);
        rSystem.fill(HIST("system/raw/like-sign/positive/hPhi"), systemPhi);
        break;

      case -2:
        rTracks.fill(HIST("tracks/pions/no-selection/like-sign/hPt"), leadingPt, subleadingPt);
        rTracks.fill(HIST("tracks/pions/no-selection/like-sign/hEta"), leadingEta, subleadingEta);
        rTracks.fill(HIST("tracks/pions/no-selection/like-sign/hPhi"), leadingPhi, subleadingPhi);
        rSystem.fill(HIST("system/raw/like-sign/negative/hM"), mass);
        rSystem.fill(HIST("system/raw/like-sign/negative/hPt"), pT);
        rSystem.fill(HIST("system/raw/like-sign/negative/hPtVsM"), mass, pT);
        rSystem.fill(HIST("system/raw/like-sign/negative/hY"), rapidity);
        rSystem.fill(HIST("system/raw/like-sign/negative/hPhi"), systemPhi);
        break;

      default:
        break;
    }

    // apply cuts to system
    if (!systemPassesCuts(system))
      return;
    fillCollisionQcHistos<1>(collision); // fill QC histograms for collisions with selected system

    // fill histograms for system passing cuts
    switch (totalCharge) {
      case 0:
        rTracks.fill(HIST("tracks/pions/selected/unlike-sign/hPt"), leadingPt, subleadingPt);
        rTracks.fill(HIST("tracks/pions/selected/unlike-sign/hEta"), leadingEta, subleadingEta);
        rTracks.fill(HIST("tracks/pions/selected/unlike-sign/hPhi"), leadingPhi, subleadingPhi);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hM"), mass);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPt"), pT);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPt2"), pTsquare);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPtVsM"), mass, pT);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hY"), rapidity);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPhi"), systemPhi);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPhiRandom"), phiRandom);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPhiCharge"), phiCharge);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        rSystem.fill(HIST("system/cut/no-selection/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        if (onon) {
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hM"), mass);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPt"), pT);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hY"), rapidity);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/0n0n/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (xnon) {
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hM"), mass);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPt"), pT);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hY"), rapidity);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/Xn0n/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (onxn) {
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hM"), mass);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPt"), pT);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hY"), rapidity);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/0nXn/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (xnxn) {
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hM"), mass);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPt"), pT);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hY"), rapidity);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/XnXn/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      case 2:
        rTracks.fill(HIST("tracks/pions/selected/like-sign/hPt"), leadingPt, subleadingPt);
        rTracks.fill(HIST("tracks/pions/selected/like-sign/hEta"), leadingEta, subleadingEta);
        rTracks.fill(HIST("tracks/pions/selected/like-sign/hPhi"), leadingPhi, subleadingPhi);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hM"), mass);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPt"), pT);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPt2"), pTsquare);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPtVsM"), mass, pT);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hY"), rapidity);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPhi"), systemPhi);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiRandom"), phiRandom);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiCharge"), phiCharge);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        rSystem.fill(HIST("system/cut/no-selection/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        if (onon) {
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hM"), mass);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPt"), pT);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hY"), rapidity);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/0n0n/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (xnon) {
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hM"), mass);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPt"), pT);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hY"), rapidity);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (onxn) {
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hM"), mass);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPt"), pT);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hY"), rapidity);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/0nXn/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (xnxn) {
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hM"), mass);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPt"), pT);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hY"), rapidity);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/XnXn/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      case -2:
        rTracks.fill(HIST("tracks/pions/selected/like-sign/hPt"), leadingPt, subleadingPt);
        rTracks.fill(HIST("tracks/pions/selected/like-sign/hEta"), leadingEta, subleadingEta);
        rTracks.fill(HIST("tracks/pions/selected/like-sign/hPhi"), leadingPhi, subleadingPhi);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hM"), mass);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPt"), pT);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPt2"), pTsquare);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPtVsM"), mass, pT);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hY"), rapidity);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPhi"), systemPhi);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiRandom"), phiRandom);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiCharge"), phiCharge);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        rSystem.fill(HIST("system/cut/no-selection/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        if (onon) {
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hM"), mass);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPt"), pT);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hY"), rapidity);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/0n0n/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (xnon) {
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hM"), mass);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPt"), pT);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hY"), rapidity);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/Xn0n/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (onxn) {
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hM"), mass);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPt"), pT);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hY"), rapidity);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/0nXn/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (xnxn) {
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hM"), mass);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPt"), pT);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPt2"), pTsquare);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPtVsM"), mass, pT);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hY"), rapidity);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPhi"), systemPhi);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiRandom"), phiRandom);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiCharge"), phiCharge);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          rSystem.fill(HIST("system/cut/XnXn/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      default:
        break;
    }
  }

  template <typename C, typename T>
  void processMC(C const& mcCollision, T const& mcParticles)
  {
    rMC.fill(HIST("MC/collisions/hPosXY"), mcCollision.posX(), mcCollision.posY());
    rMC.fill(HIST("MC/collisions/hPosZ"), mcCollision.posZ());

    std::vector<decltype(mcParticles.begin())> cutMcParticles;
    std::vector<TLorentzVector> mcParticles4Vecs;

    for (auto const& mcParticle : mcParticles) {
      rMC.fill(HIST("MC/tracks/all/hPdgCode"), mcParticle.pdgCode());
      rMC.fill(HIST("MC/tracks/all/hProducedByGenerator"), mcParticle.producedByGenerator());
      rMC.fill(HIST("MC/tracks/all/hIsPhysicalPrimary"), mcParticle.isPhysicalPrimary());
      rMC.fill(HIST("MC/tracks/all/hPt"), pt(mcParticle.px(), mcParticle.py()));
      rMC.fill(HIST("MC/tracks/all/hEta"), eta(mcParticle.px(), mcParticle.py(), mcParticle.pz()));
      rMC.fill(HIST("MC/tracks/all/hPhi"), phi(mcParticle.px(), mcParticle.py()));
      if (!mcParticle.isPhysicalPrimary() || std::abs(mcParticle.pdgCode()) != 211)
        continue;
      cutMcParticles.push_back(mcParticle);
      TLorentzVector pion4Vec;
      pion4Vec.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());
      mcParticles4Vecs.push_back(pion4Vec);
    }
    rMC.fill(HIST("MC/collisions/hNPions"), cutMcParticles.size());
    
    if (cutMcParticles.size() != 2)
      return;
    if (mcParticles4Vecs.size() != cutMcParticles.size())
      return;
    if (tracksTotalChargeMC(cutMcParticles) != 0) // shouldn't happen in theory
      return;

    TLorentzVector system = reconstructSystem(mcParticles4Vecs);
    float mass = system.M();
    float pT = system.Pt();
    float pTsquare = pT * pT;
    float rapidity = system.Rapidity();
    float systemPhi = system.Phi() + o2::constants::math::PI;
    float phiRandom = getPhiRandom(mcParticles4Vecs);
    float phiCharge = getPhiChargeMC(cutMcParticles, mcParticles4Vecs);

    auto leadingMomentumPion = momentum(cutMcParticles[0].px(), cutMcParticles[0].py(), cutMcParticles[0].pz()) > momentum(cutMcParticles[1].px(), cutMcParticles[1].py(), cutMcParticles[1].pz()) ? cutMcParticles[0] : cutMcParticles[1];
    auto subleadingMomentumPion = (leadingMomentumPion == cutMcParticles[0]) ? cutMcParticles[1] : cutMcParticles[0];
    rMC.fill(HIST("MC/tracks/pions/hPt"), pt(leadingMomentumPion.px(), leadingMomentumPion.py()), pt(subleadingMomentumPion.px(), subleadingMomentumPion.py()));
    rMC.fill(HIST("MC/tracks/pions/hEta"), eta(leadingMomentumPion.px(), leadingMomentumPion.py(), leadingMomentumPion.pz()), eta(subleadingMomentumPion.px(), subleadingMomentumPion.py(), subleadingMomentumPion.pz()));
    rMC.fill(HIST("MC/tracks/pions/hPhi"), phi(leadingMomentumPion.px(), leadingMomentumPion.py()), phi(subleadingMomentumPion.px(), subleadingMomentumPion.py()));

    rMC.fill(HIST("MC/system/hM"), mass);
    rMC.fill(HIST("MC/system/hPt"), pT);
    rMC.fill(HIST("MC/system/hPtVsM"), mass, pT);
    rMC.fill(HIST("MC/system/hPt2"), pTsquare);
    rMC.fill(HIST("MC/system/hY"), rapidity);
    rMC.fill(HIST("MC/system/hPhi"), systemPhi);
    rMC.fill(HIST("MC/system/hPhiRandom"), phiRandom);
    rMC.fill(HIST("MC/system/hPhiCharge"), phiCharge);

    // fill mcTree
    int localBc = mcCollision.globalBC() % o2::constants::lhc::LHCMaxBunches;
    std::vector<int> trackSigns = {leadingMomentumPion.pdgCode() / std::abs(leadingMomentumPion.pdgCode()), subleadingMomentumPion.pdgCode() / std::abs(subleadingMomentumPion.pdgCode())};
    std::vector<float> trackPts = {pt(leadingMomentumPion.px(), leadingMomentumPion.py()), pt(subleadingMomentumPion.px(), subleadingMomentumPion.py())};
    std::vector<float> trackEtas = {eta(leadingMomentumPion.px(), leadingMomentumPion.py(), leadingMomentumPion.pz()), eta(subleadingMomentumPion.px(), subleadingMomentumPion.py(), subleadingMomentumPion.pz())};
    std::vector<float> trackPhis = {phi(leadingMomentumPion.px(), leadingMomentumPion.py()), phi(subleadingMomentumPion.px(), subleadingMomentumPion.py())};
    std::vector<float> trackMs = {o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged};
    mcTree(localBc,
           mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), pT, rapidity, systemPhi, mass, phiRandom, phiCharge,
           trackSigns, trackPts, trackEtas, trackPhis, trackMs);
  }

  template <typename C>
  void checkNumberOfCollisionReconstructions(C const& collisions)
  {
    rMC.fill(HIST("MC/collisions/hNumOfCollisionRecos"), collisions.size());
  }

  void processSGdata(FullUdSgCollision const& collision, FullUdTracks const& tracks)
  {
    if (collision.gapSide() != 2)
      return;
    processReco(collision, tracks);
  }
  PROCESS_SWITCH(upcRhoAnalysis, processSGdata, "analyse SG data", true);

  void processDGdata(FullUdDgCollision const& collision, FullUdTracks const& tracks)
  {
    processReco(collision, tracks);
  }
  PROCESS_SWITCH(upcRhoAnalysis, processDGdata, "analyse DG data", false);

  void processMCdata(aod::UDMcCollision const& mcCollision, aod::UDMcParticles const& mcParticles)
  {
    processMC(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(upcRhoAnalysis, processMCdata, "analyse MC data", false);

  void processCollisionRecoCheck(aod::McCollision const& /* mcCollision */, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions)
  {
    checkNumberOfCollisionReconstructions(collisions);
  }
  PROCESS_SWITCH(upcRhoAnalysis, processCollisionRecoCheck, "check number of collision reconstructions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<upcRhoAnalysis>(cfgc)};
}
