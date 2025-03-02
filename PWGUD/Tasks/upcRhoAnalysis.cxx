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
#include <string_view>
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
// event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(LocalBC, localBC, int);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int);
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
// pion tracks
DECLARE_SOA_COLUMN(PhiRandom, phiRandom, float);
DECLARE_SOA_COLUMN(PhiCharge, phiCharge, float);
DECLARE_SOA_COLUMN(TrackSign, trackSign, int[2]);
DECLARE_SOA_COLUMN(TrackPt, trackPt, float[2]);
DECLARE_SOA_COLUMN(TrackEta, trackEta, float[2]);
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, float[2]);
DECLARE_SOA_COLUMN(TrackPiPID, trackPiPID, float[2]);
DECLARE_SOA_COLUMN(TrackElPID, trackElPID, float[2]);
DECLARE_SOA_COLUMN(TrackKaPID, trackKaPID, float[2]);
DECLARE_SOA_COLUMN(TrackDcaXY, trackDcaXY, float[2]);
DECLARE_SOA_COLUMN(TrackDcaZ, trackDcaZ, float[2]);
DECLARE_SOA_COLUMN(TrackTpcSignal, trackTpcSignal, float[2]);
} // namespace reco_tree
DECLARE_SOA_TABLE(RecoTree, "AOD", "RECOTREE",
                  reco_tree::RunNumber, reco_tree::LocalBC, reco_tree::NumContrib, reco_tree::PosX, reco_tree::PosY, reco_tree::PosZ,
                  reco_tree::TotalFT0AmplitudeA, reco_tree::TotalFT0AmplitudeC, reco_tree::TotalFV0AmplitudeA, reco_tree::TotalFDDAmplitudeA, reco_tree::TotalFDDAmplitudeC,
                  reco_tree::TimeFT0A, reco_tree::TimeFT0C, reco_tree::TimeFV0A, reco_tree::TimeFDDA, reco_tree::TimeFDDC,
                  reco_tree::EnergyCommonZNA, reco_tree::EnergyCommonZNC, reco_tree::TimeZNA, reco_tree::TimeZNC,
                  reco_tree::PhiRandom, reco_tree::PhiCharge, reco_tree::TrackSign, reco_tree::TrackPt, reco_tree::TrackEta, reco_tree::TrackPhi, reco_tree::TrackPiPID, reco_tree::TrackElPID, reco_tree::TrackKaPID, reco_tree::TrackDcaXY, reco_tree::TrackDcaZ, reco_tree::TrackTpcSignal);

namespace mc_tree
{
// misc event info
DECLARE_SOA_COLUMN(LocalBc, localBc, int);
// event vertex
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
// pion tracks
DECLARE_SOA_COLUMN(PhiRandom, phiRandom, float);
DECLARE_SOA_COLUMN(PhiCharge, phiCharge, float);
DECLARE_SOA_COLUMN(TrackSign, trackSign, int[2]);
DECLARE_SOA_COLUMN(TrackPt, trackPt, float[2]);
DECLARE_SOA_COLUMN(TrackEta, trackEta, float[2]);
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, float[2]);
} // namespace mc_tree
DECLARE_SOA_TABLE(McTree, "AOD", "MCTREE",
                  mc_tree::LocalBc,
                  mc_tree::PosX, mc_tree::PosY, mc_tree::PosZ,
                  mc_tree::PhiRandom, mc_tree::PhiCharge, mc_tree::TrackSign, mc_tree::TrackPt, mc_tree::TrackEta, mc_tree::TrackPhi);
} // namespace o2::aod

struct UpcRhoAnalysis {
  Produces<o2::aod::RecoTree> recoTree;
  Produces<o2::aod::McTree> mcTree;

  Configurable<bool> savePions{"savePions", true, "save pion tracks into derived tables"};
  Configurable<bool> saveElectrons{"saveElectrons", false, "save electron tracks into derived tables"};
  Configurable<bool> saveKaons{"saveKaons", false, "save kaon tracks into derived tables"};

  float pcEtaCut = 0.9; // physics coordination recommendation
  Configurable<bool> requireTof{"requireTof", false, "require TOF signal"};

  Configurable<float> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "max Z position cut on collisions"};
  Configurable<int> collisionsNumContribsMaxCut{"collisionsNumContribsMaxCut", 4, "max number of contributors cut on collisions"};
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
  Configurable<float> systemPtCut{"systemPtCut", 0.1, "max pT cut for reco system"};
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

    // tracks
    rQC.add("QC/tracks/all/hTpcNSigmaPi", ";TPC n#sigma(#pi);counts", kTH1D, {{400, -10.0, 30.0}});
    rQC.add("QC/tracks/all/hTpcNSigmaEl", ";TPC n#sigma(e);counts", kTH1D, {{400, -10.0, 30.0}});
    rQC.add("QC/tracks/all/hTpcNSigmaKa", ";TPC n#sigma(K);counts", kTH1D, {{400, -10.0, 30.0}});
    rQC.add("QC/tracks/all/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {1000, -5.0, 5.0}});
    rQC.add("QC/tracks/all/hItsNCls", ";ITS N_{cls};counts", kTH1D, {{11, -0.5, 10.5}});
    rQC.add("QC/tracks/all/hItsChi2NCl", ";ITS #chi^{2}/N_{cls};counts", kTH1D, {{1000, 0.0, 100.0}});
    rQC.add("QC/tracks/all/hTpcChi2NCl", ";TPC #chi^{2}/N_{cls};counts", kTH1D, {{1000, 0.0, 100.0}});
    rQC.add("QC/tracks/all/hTpcNCls", ";TPC N_{cls} found;counts", kTH1D, {{200, 0.0, 200.0}});
    rQC.add("QC/tracks/all/hTpcNClsCrossedRows", ";TPC crossed rows;counts", kTH1D, {{200, 0.0, 200.0}});
    rQC.add("QC/tracks/all/hTpcNClsCrossedRowsOverNClsFindable", ";TPC crossed rows/findable N_{cls};counts", kTH1D, {{100, 0.0, 10.0}});
    rQC.add("QC/tracks/all/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rQC.add("QC/tracks/all/hEta", ";y;counts", kTH1D, {etaAxis});
    rQC.add("QC/tracks/all/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    rQC.add("QC/tracks/all/hTpcSignalVsP", ";p (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    rQC.add("QC/tracks/all/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    // tracks passing selections
    rQC.addClone("QC/tracks/all/", "QC/tracks/selected/"); // clone "raw" histograms as "cut"
    rQC.add("QC/tracks/selected/hRemainingTracks", ";remaining tracks;counts", kTH1D, {{21, -0.5, 20.5}});
    rQC.add("QC/tracks/selected/hTpcNSigmaPi2D", ";TPC n#sigma(#pi_{leading});TPC n#sigma(#pi_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    rQC.add("QC/tracks/selected/hTpcNSigmaEl2D", ";TPC n#sigma(e_{leading});TPC n#sigma(e_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    rQC.add("QC/tracks/selected/hTpcNSigmaKa2D", ";TPC n#sigma(K_{leading});TPC n#sigma(K_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    // selection counter
    std::vector<std::string> selectionCounterLabels = {"all tracks", "PV contributor", "ITS hit", "ITS N_{clusters}", "ITS #chi^{2}/N_{clusters}", "TPC hit", "TPC N_{clusters} found", "TPC #chi^{2}/N_{clusters}", "TPC crossed rows",
                                                       "TPC crossed rows/N_{clusters}",
                                                       "TOF requirement",
                                                       "p_{T}", "DCA", "#eta", "exactly 2 tracks"};
    rQC.add("QC/tracks/hSelectionCounter", ";;tracks passing selections", kTH1D, {{static_cast<int>(selectionCounterLabels.size()), -0.5, static_cast<float>(selectionCounterLabels.size()) - 0.5}});
    for (int i = 0; i < static_cast<int>(selectionCounterLabels.size()); ++i)
      rQC.get<TH1>(HIST("QC/tracks/hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, selectionCounterLabels[i].c_str());
    rQC.add("QC/tracks/hTofHitCheck", ";leading track TOF hit;subleading track TOF hit;counts", kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});

    // TRACKS (2D)
    rTracks.add("tracks/all/unlike-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    rTracks.add("tracks/all/unlike-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaCutAxis, etaCutAxis});
    rTracks.add("tracks/all/unlike-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    rTracks.addClone("tracks/all/unlike-sign/", "tracks/all/like-sign/positive/");
    rTracks.addClone("tracks/all/unlike-sign/", "tracks/all/like-sign/negative/");
    rTracks.addClone("tracks/all/", "tracks/selected/");

    // SYSTEM
    rSystem.add("system/all/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    rSystem.add("system/all/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rSystem.add("system/all/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    rSystem.add("system/all/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    rSystem.add("system/all/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    rSystem.add("system/all/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    rSystem.add("system/all/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    rSystem.add("system/all/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    // clones for like-sign
    rSystem.addClone("system/all/unlike-sign/", "system/all/like-sign/positive/");
    rSystem.addClone("system/all/unlike-sign/", "system/all/like-sign/negative/");
    // selected rhos
    rSystem.addClone("system/all/", "system/selected/no-selection/");
    // clones for neutron classes
    rSystem.addClone("system/selected/no-selection/", "system/selected/0n0n/");
    rSystem.addClone("system/selected/no-selection/", "system/selected/Xn0n/");
    rSystem.addClone("system/selected/no-selection/", "system/selected/0nXn/");
    rSystem.addClone("system/selected/no-selection/", "system/selected/XnXn/");

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
    rMC.add("MC/tracks/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    rMC.add("MC/tracks/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaAxis, etaAxis});
    rMC.add("MC/tracks/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    // system
    rMC.add("MC/system/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    rMC.add("MC/system/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rMC.add("MC/system/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    rMC.add("MC/system/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    rMC.add("MC/system/hY", ";y;counts", kTH1D, {yAxis});
    rMC.add("MC/system/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    rMC.add("MC/system/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    rMC.add("MC/system/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    rMC.addClone("MC/system/", "MC/system/selected/");
  }

  static constexpr std::string_view AppliedSelections[2] = {"all/", "selected/"};
  static constexpr std::string_view ChargeLabel[3] = {"unlike-sign/", "like-sign/positive/", "like-sign/negative/"};
  static constexpr std::string_view NeutronClass[5] = {"no-selection/", "0n0n/", "Xn0n/", "0nXn/", "XnXn/"};

  template <int cuts, typename C>
  void fillCollisionQcHistos(const C& collision) // fills collision QC histograms before/after cuts
  {
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hPosXY"), collision.posX(), collision.posY());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hPosZ"), collision.posZ());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hZdcCommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hZdcTime"), collision.timeZNA(), collision.timeZNC());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hNumContrib"), collision.numContrib());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTotalFT0AmplitudeA"), collision.totalFT0AmplitudeA());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTotalFT0AmplitudeC"), collision.totalFT0AmplitudeC());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTotalFV0AmplitudeA"), collision.totalFV0AmplitudeA());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTotalFDDAmplitudeA"), collision.totalFDDAmplitudeA());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTotalFDDAmplitudeC"), collision.totalFDDAmplitudeC());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTimeFT0A"), collision.timeFT0A());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTimeFT0C"), collision.timeFT0C());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTimeFV0A"), collision.timeFV0A());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTimeFDDA"), collision.timeFDDA());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hTimeFDDC"), collision.timeFDDC());
  }

  template <int cuts, typename T>
  void fillTrackQcHistos(const T& track)
  {
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hPt"), track.pt());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hEta"), eta(track.px(), track.py(), track.pz()));
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hPhi"), phi(track.px(), track.py()));
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcNSigmaPi"), track.tpcNSigmaPi());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcNSigmaEl"), track.tpcNSigmaEl());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcNSigmaKa"), track.tpcNSigmaKa());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hDcaXYZ"), track.dcaZ(), track.dcaXY());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hItsNCls"), track.itsNCls());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hItsChi2NCl"), track.itsChi2NCl());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcChi2NCl"), track.tpcChi2NCl());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcNCls"), (track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()));
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcNClsCrossedRows"), track.tpcNClsCrossedRows());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcNClsCrossedRowsOverNClsFindable"), (static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())));
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcSignalVsP"), momentum(track.px(), track.py(), track.pz()), track.tpcSignal());
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcSignalVsPt"), track.pt(), track.tpcSignal());
  }

  template <int cuts, int charge>
  void fillTrack2dHistos(float leadingPt, float subleadingPt, float leadingEta, float subleadingEta, float leadingPhi, float subleadingPhi)
  {
    rTracks.fill(HIST("tracks/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPt"), leadingPt, subleadingPt);
    rTracks.fill(HIST("tracks/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hEta"), leadingEta, subleadingEta);
    rTracks.fill(HIST("tracks/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPhi"), leadingPhi, subleadingPhi);
  }

  template <int cuts, int neutronClass, int charge>
  void fillSystemHistos(float mass, float pt, float rapidity, float phi, float phiRandom, float phiCharge)
  {
    if (cuts == 0) {
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hM"), mass);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPt"), pt);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPt2"), pt * pt);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPtVsM"), mass, pt);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hY"), rapidity);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPhi"), phi);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPhiRandom"), phiRandom);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPhiCharge"), phiCharge);
    } else {
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hM"), mass);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPt"), pt);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPt2"), pt * pt);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPtVsM"), mass, pt);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hY"), rapidity);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPhi"), phi);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPhiRandom"), phiRandom);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPhiCharge"), phiCharge);
    }
  }

  template <typename C>
  bool collisionPassesCuts(const C& collision) // collision cuts
  {
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut)
      return false;
    if (collision.numContrib() > collisionsNumContribsMaxCut)
      return false;
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
  bool tracksPassKaPID(const T& cutTracks) // n-dimensional kaon PID cut
  {
    float radius = 0.0;
    for (const auto& track : cutTracks)
      radius += std::pow(track.tpcNSigmaKa(), 2);
    return radius < std::pow(tracksTpcNSigmaKaCut, 2);
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
  int tracksTotalChargeMC(const T& cutTracks) // total charge of selected MC tracks
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

  TLorentzVector reconstructSystem(const std::vector<TLorentzVector>& cutTracksLVs) // reconstruct system from 4-vectors
  {
    TLorentzVector system;
    for (const auto& trackLV : cutTracksLVs)
      system += trackLV;
    return system;
  }

  float getPhiRandom(const std::vector<TLorentzVector>& cutTracksLVs)    // decay phi anisotropy
  {                                                                      // two possible definitions of phi: randomize the tracks
    int indices[2] = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // get time-based seed
    std::shuffle(std::begin(indices), std::end(indices), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    TLorentzVector pOne = cutTracksLVs[indices[0]];
    TLorentzVector pTwo = cutTracksLVs[indices[1]];
    TLorentzVector pPlus = pOne + pTwo;
    TLorentzVector pMinus = pOne - pTwo;
    return pPlus.DeltaPhi(pMinus);
  }

  template <typename T>
  float getPhiCharge(const T& cutTracks, const std::vector<TLorentzVector>& cutTracksLVs)
  { // two possible definitions of phi: charge-based assignment
    TLorentzVector pOne, pTwo;
    pOne = (cutTracks[0].sign() > 0) ? cutTracksLVs[0] : cutTracksLVs[1];
    pTwo = (cutTracks[0].sign() > 0) ? cutTracksLVs[1] : cutTracksLVs[0];
    TLorentzVector pPlus = pOne + pTwo;
    TLorentzVector pMinus = pOne - pTwo;
    return pPlus.DeltaPhi(pMinus);
  }

  template <typename T>
  float getPhiChargeMC(const T& cutTracks, const std::vector<TLorentzVector>& cutTracksLVs)
  { // the same as for data but using pdg code instead of charge
    TLorentzVector pOne, pTwo;
    pOne = (cutTracks[0].pdgCode() > 0) ? cutTracksLVs[0] : cutTracksLVs[1];
    pTwo = (cutTracks[0].pdgCode() > 0) ? cutTracksLVs[1] : cutTracksLVs[0];
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

    bool xnxn = false, onon = false, xnon = false, onxn = false; // note: On == 0n...
    if (collision.energyCommonZNA() < znCommonEnergyCut && collision.energyCommonZNC() < znCommonEnergyCut)
      onon = true;
    if (collision.energyCommonZNA() > znCommonEnergyCut && std::abs(collision.timeZNA()) < znTimeCut && collision.energyCommonZNC() < znCommonEnergyCut)
      xnon = true;
    if (collision.energyCommonZNA() < znCommonEnergyCut && collision.energyCommonZNC() > znCommonEnergyCut && std::abs(collision.timeZNC()) < znTimeCut)
      onxn = true;
    if (collision.energyCommonZNA() > znCommonEnergyCut && std::abs(collision.timeZNA()) < znTimeCut &&
        collision.energyCommonZNC() > znCommonEnergyCut && std::abs(collision.timeZNC()) < znTimeCut)
      xnxn = true;

    std::vector<decltype(tracks.begin())> cutTracks; // store selected tracks
    for (const auto& track : tracks) {
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 0);
      fillTrackQcHistos<0>(track); // fill QC histograms before cuts

      if (!trackPassesCuts(track)) // apply track cuts
        continue;

      fillTrackQcHistos<1>(track); // fill QC histograms after cuts
      cutTracks.push_back(track);
    }
    rQC.fill(HIST("QC/tracks/selected/hRemainingTracks"), cutTracks.size());

    if (cutTracks.size() != 2) // further consider only two pion systems
      return;
    for (int i = 0; i < static_cast<int>(cutTracks.size()); i++)
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 14);
    rQC.fill(HIST("QC/tracks/selected/hTpcNSigmaPi2D"), cutTracks[0].tpcNSigmaPi(), cutTracks[1].tpcNSigmaPi());
    rQC.fill(HIST("QC/tracks/selected/hTpcNSigmaEl2D"), cutTracks[0].tpcNSigmaEl(), cutTracks[1].tpcNSigmaEl());
    rQC.fill(HIST("QC/tracks/selected/hTpcNSigmaKa2D"), cutTracks[0].tpcNSigmaKa(), cutTracks[1].tpcNSigmaKa());

    // create a vector of 4-vectors for selected tracks
    std::vector<TLorentzVector> cutTracksLVs;
    for (const auto& track : cutTracks) {
      TLorentzVector trackLV;
      trackLV.SetXYZM(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged); // apriori assume pion mass
      cutTracksLVs.push_back(trackLV);
    }

    // differentiate leading- and subleading-momentum tracks
    auto leadingMomentumTrack = momentum(cutTracks[0].px(), cutTracks[0].py(), cutTracks[0].pz()) > momentum(cutTracks[1].px(), cutTracks[1].py(), cutTracks[1].pz()) ? cutTracks[0] : cutTracks[1];
    auto subleadingMomentumTrack = (leadingMomentumTrack == cutTracks[0]) ? cutTracks[1] : cutTracks[0];
    rQC.fill(HIST("QC/tracks/hTofHitCheck"), leadingMomentumTrack.hasTOF(), subleadingMomentumTrack.hasTOF());

    float leadingPt = leadingMomentumTrack.pt();
    float subleadingPt = subleadingMomentumTrack.pt();
    float leadingEta = eta(leadingMomentumTrack.px(), leadingMomentumTrack.py(), leadingMomentumTrack.pz());
    float subleadingEta = eta(subleadingMomentumTrack.px(), subleadingMomentumTrack.py(), subleadingMomentumTrack.pz());
    float leadingPhi = phi(leadingMomentumTrack.px(), leadingMomentumTrack.py());
    float subleadingPhi = phi(subleadingMomentumTrack.px(), subleadingMomentumTrack.py());
    float phiRandom = getPhiRandom(cutTracksLVs);
    float phiCharge = getPhiCharge(cutTracks, cutTracksLVs);

    // fill recoTree
    int localBc = collision.globalBC() % o2::constants::lhc::LHCMaxBunches;
    int trackSigns[2] = {leadingMomentumTrack.sign(), subleadingMomentumTrack.sign()};
    float trackPts[2] = {leadingPt, subleadingPt};
    float trackEtas[2] = {leadingEta, subleadingEta};
    float trackPhis[2] = {leadingPhi, subleadingPhi};
    float trackPiPIDs[2] = {leadingMomentumTrack.tpcNSigmaPi(), subleadingMomentumTrack.tpcNSigmaPi()};
    float trackElPIDs[2] = {leadingMomentumTrack.tpcNSigmaEl(), subleadingMomentumTrack.tpcNSigmaEl()};
    float trackKaPIDs[2] = {leadingMomentumTrack.tpcNSigmaKa(), subleadingMomentumTrack.tpcNSigmaKa()};
    float trackDcaXYs[2] = {leadingMomentumTrack.dcaXY(), subleadingMomentumTrack.dcaXY()};
    float trackDcaZs[2] = {leadingMomentumTrack.dcaZ(), subleadingMomentumTrack.dcaZ()};
    float trackTpcSignals[2] = {leadingMomentumTrack.tpcSignal(), subleadingMomentumTrack.tpcSignal()};
    if ((savePions && tracksPassPiPID(cutTracks)) || (saveElectrons && tracksPassElPID(cutTracks)) || (saveKaons && tracksPassKaPID(cutTracks)))
      recoTree(collision.runNumber(), localBc, collision.numContrib(), collision.posX(), collision.posY(), collision.posZ(),
               collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
               collision.timeFT0A(), collision.timeFT0C(), collision.timeFV0A(), collision.timeFDDA(), collision.timeFDDC(),
               collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC(),
               phiRandom, phiCharge, trackSigns, trackPts, trackEtas, trackPhis, trackPiPIDs, trackElPIDs, trackKaPIDs, trackDcaXYs, trackDcaZs, trackTpcSignals);

    if (!tracksPassPiPID(cutTracks)) // apply PID cut
      return;
    TLorentzVector system = reconstructSystem(cutTracksLVs);
    int totalCharge = tracksTotalCharge(cutTracks);
    float mass = system.M();
    float pT = system.Pt();
    float rapidity = system.Rapidity();
    float systemPhi = system.Phi() + o2::constants::math::PI;

    // fill raw histograms according to total charge
    switch (totalCharge) {
      case 0:
        fillTrack2dHistos<0, 0>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<0, 0, 0>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        break;

      case 2:
        fillTrack2dHistos<0, 1>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<0, 0, 1>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        break;

      case -2:
        fillTrack2dHistos<0, 2>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<0, 0, 2>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
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
        fillTrack2dHistos<1, 0>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<1, 0, 0>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (onon)
          fillSystemHistos<1, 1, 0>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (xnon)
          fillSystemHistos<1, 2, 0>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (onxn)
          fillSystemHistos<1, 3, 0>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (xnxn)
          fillSystemHistos<1, 4, 0>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        break;

      case 2:
        fillTrack2dHistos<1, 1>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<1, 0, 1>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (onon)
          fillSystemHistos<1, 1, 1>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (xnon)
          fillSystemHistos<1, 2, 1>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (onxn)
          fillSystemHistos<1, 3, 1>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (xnxn)
          fillSystemHistos<1, 4, 1>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        break;

      case -2:
        fillTrack2dHistos<1, 2>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<1, 0, 2>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (onon)
          fillSystemHistos<1, 1, 2>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (xnon)
          fillSystemHistos<1, 2, 2>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (onxn)
          fillSystemHistos<1, 3, 2>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        if (xnxn)
          fillSystemHistos<1, 4, 2>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
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
    std::vector<TLorentzVector> mcParticlesLVs;

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
      TLorentzVector pionLV;
      pionLV.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());
      mcParticlesLVs.push_back(pionLV);
    }
    rMC.fill(HIST("MC/collisions/hNPions"), cutMcParticles.size());

    if (cutMcParticles.size() != 2)
      return;
    if (mcParticlesLVs.size() != cutMcParticles.size())
      return;
    if (tracksTotalChargeMC(cutMcParticles) != 0) // shouldn't happen in theory
      return;

    TLorentzVector system = reconstructSystem(mcParticlesLVs);
    float mass = system.M();
    float pT = system.Pt();
    float rapidity = system.Rapidity();
    float systemPhi = system.Phi() + o2::constants::math::PI;
    float phiRandom = getPhiRandom(mcParticlesLVs);
    float phiCharge = getPhiChargeMC(cutMcParticles, mcParticlesLVs);

    auto leadingMomentumPion = momentum(cutMcParticles[0].px(), cutMcParticles[0].py(), cutMcParticles[0].pz()) > momentum(cutMcParticles[1].px(), cutMcParticles[1].py(), cutMcParticles[1].pz()) ? cutMcParticles[0] : cutMcParticles[1];
    auto subleadingMomentumPion = (leadingMomentumPion == cutMcParticles[0]) ? cutMcParticles[1] : cutMcParticles[0];
    rMC.fill(HIST("MC/tracks/hPt"), pt(leadingMomentumPion.px(), leadingMomentumPion.py()), pt(subleadingMomentumPion.px(), subleadingMomentumPion.py()));
    rMC.fill(HIST("MC/tracks/hEta"), eta(leadingMomentumPion.px(), leadingMomentumPion.py(), leadingMomentumPion.pz()), eta(subleadingMomentumPion.px(), subleadingMomentumPion.py(), subleadingMomentumPion.pz()));
    rMC.fill(HIST("MC/tracks/hPhi"), phi(leadingMomentumPion.px(), leadingMomentumPion.py()), phi(subleadingMomentumPion.px(), subleadingMomentumPion.py()));

    rMC.fill(HIST("MC/system/hM"), mass);
    rMC.fill(HIST("MC/system/hPt"), pT);
    rMC.fill(HIST("MC/system/hPtVsM"), mass, pT);
    rMC.fill(HIST("MC/system/hPt2"), pT * pT);
    rMC.fill(HIST("MC/system/hY"), rapidity);
    rMC.fill(HIST("MC/system/hPhi"), systemPhi);
    rMC.fill(HIST("MC/system/hPhiRandom"), phiRandom);
    rMC.fill(HIST("MC/system/hPhiCharge"), phiCharge);

    if (systemPassesCuts(system)) {
      rMC.fill(HIST("MC/system/selected/hM"), mass);
      rMC.fill(HIST("MC/system/selected/hPt"), pT);
      rMC.fill(HIST("MC/system/selected/hPtVsM"), mass, pT);
      rMC.fill(HIST("MC/system/selected/hPt2"), pT * pT);
      rMC.fill(HIST("MC/system/selected/hY"), rapidity);
      rMC.fill(HIST("MC/system/selected/hPhi"), systemPhi);
      rMC.fill(HIST("MC/system/selected/hPhiRandom"), phiRandom);
      rMC.fill(HIST("MC/system/selected/hPhiCharge"), phiCharge);
    }

    // fill mcTree
    int localBc = mcCollision.globalBC() % o2::constants::lhc::LHCMaxBunches;
    int trackSigns[2] = {leadingMomentumPion.pdgCode() / std::abs(leadingMomentumPion.pdgCode()), subleadingMomentumPion.pdgCode() / std::abs(subleadingMomentumPion.pdgCode())};
    float trackPts[2] = {pt(leadingMomentumPion.px(), leadingMomentumPion.py()), pt(subleadingMomentumPion.px(), subleadingMomentumPion.py())};
    float trackEtas[2] = {eta(leadingMomentumPion.px(), leadingMomentumPion.py(), leadingMomentumPion.pz()), eta(subleadingMomentumPion.px(), subleadingMomentumPion.py(), subleadingMomentumPion.pz())};
    float trackPhis[2] = {phi(leadingMomentumPion.px(), leadingMomentumPion.py()), phi(subleadingMomentumPion.px(), subleadingMomentumPion.py())};
    mcTree(localBc,
           mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
           phiRandom, phiCharge, trackSigns, trackPts, trackEtas, trackPhis);
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
  PROCESS_SWITCH(UpcRhoAnalysis, processSGdata, "analyse SG data", true);

  void processDGdata(FullUdDgCollision const& collision, FullUdTracks const& tracks)
  {
    processReco(collision, tracks);
  }
  PROCESS_SWITCH(UpcRhoAnalysis, processDGdata, "analyse DG data", false);

  void processMCdata(aod::UDMcCollision const& mcCollision, aod::UDMcParticles const& mcParticles)
  {
    processMC(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(UpcRhoAnalysis, processMCdata, "analyse MC data", false);

  void processCollisionRecoCheck(aod::McCollision const& /* mcCollision */, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions)
  {
    checkNumberOfCollisionReconstructions(collisions);
  }
  PROCESS_SWITCH(UpcRhoAnalysis, processCollisionRecoCheck, "check number of collision reconstructions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<UpcRhoAnalysis>(cfgc)};
}
