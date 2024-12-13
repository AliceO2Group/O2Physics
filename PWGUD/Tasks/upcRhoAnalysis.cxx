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
namespace tree
{
// misc event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int);
// event vertex
DECLARE_SOA_COLUMN(PosX, posX, double);
DECLARE_SOA_COLUMN(PosY, posY, double);
DECLARE_SOA_COLUMN(PosZ, posZ, double);
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
DECLARE_SOA_COLUMN(NeutronClass, neutronClass, int);
// Rhos
DECLARE_SOA_COLUMN(TotalCharge, totalCharge, int);
// store things for reconstruction of lorentz vectors
DECLARE_SOA_COLUMN(RhoPt, rhoPt, double);
DECLARE_SOA_COLUMN(RhoEta, rhoEta, double);
DECLARE_SOA_COLUMN(RhoPhi, rhoPhi, double);
DECLARE_SOA_COLUMN(RhoM, rhoM, double);
// other stuff
DECLARE_SOA_COLUMN(RhoPhiRandom, rhoPhiRandom, double);
DECLARE_SOA_COLUMN(RhoPhiCharge, rhoPhiCharge, double);
// Pion tracks
DECLARE_SOA_COLUMN(TrackSign, trackSign, std::vector<int>);
// for lorentz vector reconstruction
DECLARE_SOA_COLUMN(TrackPt, trackPt, std::vector<double>);
DECLARE_SOA_COLUMN(TrackEta, trackEta, std::vector<double>);
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, std::vector<double>);
DECLARE_SOA_COLUMN(TrackM, trackM, std::vector<double>);
// other stuff
DECLARE_SOA_COLUMN(TrackPiPID, trackPiPID, std::vector<double>);
DECLARE_SOA_COLUMN(TrackElPID, trackElPID, std::vector<double>);
DECLARE_SOA_COLUMN(TrackDcaXY, trackDcaXY, std::vector<double>);
DECLARE_SOA_COLUMN(TrackDcaZ, trackDcaZ, std::vector<double>);
DECLARE_SOA_COLUMN(TrackTpcSignal, trackTpcSignal, std::vector<double>);
} // namespace tree
DECLARE_SOA_TABLE(Tree, "AOD", "TREE",
                  tree::RunNumber, tree::GlobalBC, tree::NumContrib,
                  tree::PosX, tree::PosY, tree::PosZ, tree::TotalFT0AmplitudeA, tree::TotalFT0AmplitudeC, tree::TotalFV0AmplitudeA, tree::TotalFDDAmplitudeA, tree::TotalFDDAmplitudeC,
                  tree::TimeFT0A, tree::TimeFT0C, tree::TimeFV0A, tree::TimeFDDA, tree::TimeFDDC,
                  tree::EnergyCommonZNA, tree::EnergyCommonZNC, tree::TimeZNA, tree::TimeZNC, tree::NeutronClass,
                  tree::TotalCharge, tree::RhoPt, tree::RhoEta, tree::RhoPhi, tree::RhoM, tree::RhoPhiRandom, tree::RhoPhiCharge,
                  tree::TrackSign, tree::TrackPt, tree::TrackEta, tree::TrackPhi, tree::TrackM, tree::TrackPiPID, tree::TrackElPID, tree::TrackDcaXY, tree::TrackDcaZ, tree::TrackTpcSignal);
} // namespace o2::aod

struct upcRhoAnalysis {
  Produces<o2::aod::Tree> Tree;

  double PcEtaCut = 0.9; // physics coordination recommendation
  Configurable<bool> requireTof{"requireTof", false, "require TOF signal"};
  Configurable<bool> do4pi{"do4pi", true, "do 4pi analysis"};

  Configurable<double> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "max Z position cut on collisions"};
  Configurable<double> ZNcommonEnergyCut{"ZNcommonEnergyCut", 0.0, "ZN common energy cut"};
  Configurable<double> ZNtimeCut{"ZNtimeCut", 2.0, "ZN time cut"};

  Configurable<double> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "TPC nSigma pion cut"};
  Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "max DCA cut on tracks"};
  Configurable<int> tracksMinItsNClsCut{"tracksMinItsNClsCut", 6, "min ITS clusters cut"};
  Configurable<double> tracksMaxItsChi2NClCut{"tracksMaxItsChi2NClCut", 3.0, "max ITS chi2/Ncls cut"};
  Configurable<int> tracksMinTpcNClsCut{"tracksMinTpcNClsCut", 120, "min TPC clusters cut"};
  Configurable<int> tracksMinTpcNClsCrossedRowsCut{"tracksMinTpcNClsCrossedRowsCut", 140, "min TPC crossed rows cut"};
  Configurable<double> tracksMinTpcChi2NClCut{"tracksMinTpcChi2NClCut", 1.0, "min TPC chi2/Ncls cut"};
  Configurable<double> tracksMaxTpcChi2NClCut{"tracksMaxTpcChi2NClCut", 1.8, "max TPC chi2/Ncls cut"};
  Configurable<double> tracksMinTpcNClsCrossedOverFindableCut{"tracksMinTpcNClsCrossedOverFindableCut", 1.05, "min TPC crossed rows / findable clusters cut"};
  Configurable<double> tracksMinPtCut{"tracksMinPtCut", 0.2, "min pT cut on tracks"};

  Configurable<double> systemMassMinCut{"systemMassMinCut", 0.4, "min M cut for reco system"};
  Configurable<double> systemMassMaxCut{"systemMassMaxCut", 1.2, "max M cut for reco system"};
  Configurable<double> systemPtCut{"systemPtMaxCut", 0.1, "max pT cut for reco system"};
  Configurable<double> systemYCut{"systemYCut", 0.9, "rapiditiy cut for reco system"};

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
  // ConfigurableAxis ptQuantileAxis{"ptQuantileAxis", {0, 0.0181689, 0.0263408, 0.0330488, 0.0390369, 0.045058, 0.0512604, 0.0582598, 0.066986, 0.0788085, 0.1}, "p_{T} (GeV/#it{c})"};

  HistogramRegistry QC{"QC", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry Pions{"Pions", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry System{"System", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MC{"MC", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry FourPiQA{"4piQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // QA //
    // collisions
    QC.add("QC/collisions/all/hPosXY", ";x (cm);y (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    QC.add("QC/collisions/all/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    QC.add("QC/collisions/all/hNumContrib", ";number of contributors;counts", kTH1D, {{36, -0.5, 35.5}});
    QC.add("QC/collisions/all/hZdcCommonEnergy", ";ZNA common energy;ZNC common energy;counts", kTH2D, {{250, -5.0, 20.0}, {250, -5.0, 20.0}});
    QC.add("QC/collisions/all/hZdcTime", ";ZNA time (ns);ZNC time (ns);counts", kTH2D, {{200, -10.0, 10.0}, {200, -10.0, 10.0}});
    QC.add("QC/collisions/all/hTotalFT0AmplitudeA", ";FT0A amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/all/hTotalFT0AmplitudeC", ";FT0C amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/all/hTotalFV0AmplitudeA", ";FV0A amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/all/hTotalFDDAmplitudeA", ";FDDA amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/all/hTotalFDDAmplitudeC", ";FDDC amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/all/hTimeFT0A", ";FT0A time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    QC.add("QC/collisions/all/hTimeFT0C", ";FT0C time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    QC.add("QC/collisions/all/hTimeFV0A", ";FV0A time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    QC.add("QC/collisions/all/hTimeFDDA", ";FDDA time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    QC.add("QC/collisions/all/hTimeFDDC", ";FDDC time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    // events with selected rho candidates
    QC.add("QC/collisions/selected/hPosXY", ";x (cm);y (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    QC.add("QC/collisions/selected/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    QC.add("QC/collisions/selected/hNumContrib", ";number of contributors;counts", kTH1D, {{36, -0.5, 35.5}});
    QC.add("QC/collisions/selected/hZdcCommonEnergy", ";ZNA common energy;ZNC common energy;counts", kTH2D, {{250, -5.0, 20.0}, {250, -5.0, 20.0}});
    QC.add("QC/collisions/selected/hZdcTime", ";ZNA time (ns);ZNC time (ns);counts", kTH2D, {{200, -10.0, 10.0}, {200, -10.0, 10.0}});
    QC.add("QC/collisions/selected/hTotalFT0AmplitudeA", ";FT0A amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/selected/hTotalFT0AmplitudeC", ";FT0C amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/selected/hTotalFV0AmplitudeA", ";FV0A amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/selected/hTotalFDDAmplitudeA", ";FDDA amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/selected/hTotalFDDAmplitudeC", ";FDDC amplitude;counts", kTH1D, {{1000, 0.0, 1000.0}});
    QC.add("QC/collisions/selected/hTimeFT0A", ";FT0A time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    QC.add("QC/collisions/selected/hTimeFT0C", ";FT0C time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    QC.add("QC/collisions/selected/hTimeFV0A", ";FV0A time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    QC.add("QC/collisions/selected/hTimeFDDA", ";FDDA time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    QC.add("QC/collisions/selected/hTimeFDDC", ";FDDC time (ns);counts", kTH1D, {{200, -100.0, 100.0}});
    // all tracks
    QC.add("QC/tracks/raw/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {{400, -10.0, 30.0}});
    QC.add("QC/tracks/raw/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {{400, -20.0, 20.0}});
    QC.add("QC/tracks/raw/hTpcNSigmaEl", ";TPC n#sigma_{e};counts", kTH1D, {{400, -10.0, 30.0}});
    QC.add("QC/tracks/raw/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {1000, -5.0, 5.0}});
    QC.add("QC/tracks/raw/hItsNCls", ";ITS N_{cls};counts", kTH1D, {{11, -0.5, 10.5}});
    QC.add("QC/tracks/raw/hItsChi2NCl", ";ITS #chi^{2}/N_{cls};counts", kTH1D, {{1000, 0.0, 100.0}});
    QC.add("QC/tracks/raw/hTpcChi2NCl", ";TPC #chi^{2}/N_{cls};counts", kTH1D, {{1000, 0.0, 100.0}});
    QC.add("QC/tracks/raw/hTpcNCls", ";TPC N_{cls} found;counts", kTH1D, {{200, 0.0, 200.0}});
    QC.add("QC/tracks/raw/hTpcNClsCrossedRows", ";TPC crossed rows;counts", kTH1D, {{200, 0.0, 200.0}});
    QC.add("QC/tracks/raw/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    QC.add("QC/tracks/raw/hEta", ";y;counts", kTH1D, {etaAxis});
    QC.add("QC/tracks/raw/hPhi", ";#phi;counts", kTH1D, {phiAxis}); // tracks passing selections
    QC.add("QC/tracks/cut/hTpcNSigmaPi2D", ";TPC n#sigma(#pi_{leading});TPC n#sigma(#pi_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    QC.add("QC/tracks/cut/hTpcNSigmaEl2D", ";TPC n#sigma(e_{leading});TPC n#sigma(e_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    QC.add("QC/tracks/cut/hTpcSignalVsP", ";p (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    QC.add("QC/tracks/cut/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    QC.add("QC/tracks/cut/hRemainingTracks", ";remaining tracks;counts", kTH1D, {{21, -0.5, 20.5}});
    QC.add("QC/tracks/cut/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {1000, -5.0, 5.0}});
    // selection counter
    std::vector<std::string> selectionCounterLabels = {"all tracks", "PV contributor", "ITS hit", "ITS N_{clusters}", "ITS #chi^{2}/N_{clusters}", "TPC hit", "TPC N_{clusters} found", "TPC #chi^{2}/N_{clusters}", "TPC crossed rows",
                                                       "TPC crossed rows/N_{clusters}",
                                                       "TOF requirement",
                                                       "p_{T}", "DCA", "#eta", "exactly 2 tracks", "PID"};
    auto hSelectionCounter = QC.add<TH1>("QC/tracks/hSelectionCounter", ";;counts", kTH1D, {{static_cast<int>(selectionCounterLabels.size()), -0.5, static_cast<double>(selectionCounterLabels.size()) - 0.5}});
    for (int i = 0; i < static_cast<int>(selectionCounterLabels.size()); ++i)
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, selectionCounterLabels[i].c_str());
    // TOF hit check
    QC.add("QC/tracks/hTofHitCheck", ";leading track TOF hit;subleading track TOF hit;counts", kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
    // RECO HISTOS //
    // PIONS
    // no selection
    Pions.add("pions/no-selection/unlike-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    Pions.add("pions/no-selection/unlike-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaCutAxis, etaCutAxis});
    Pions.add("pions/no-selection/unlike-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    Pions.add("pions/no-selection/like-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    Pions.add("pions/no-selection/like-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaCutAxis, etaCutAxis});
    Pions.add("pions/no-selection/like-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    // selected
    Pions.add("pions/selected/unlike-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    Pions.add("pions/selected/unlike-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaCutAxis, etaCutAxis});
    Pions.add("pions/selected/unlike-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    Pions.add("pions/selected/like-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    Pions.add("pions/selected/like-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaCutAxis, etaCutAxis});
    Pions.add("pions/selected/like-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});

    // RAW RHOS
    System.add("system/raw/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    System.add("system/raw/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    System.add("system/raw/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    System.add("system/raw/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/raw/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    System.add("system/raw/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    System.add("system/raw/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    System.add("system/raw/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    System.add("system/raw/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/raw/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    System.add("system/raw/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    System.add("system/raw/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    System.add("system/raw/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    System.add("system/raw/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/raw/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});

    // SELECTED RHOS
    // no selection
    System.add("system/cut/no-selection/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/no-selection/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/no-selection/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/no-selection/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/no-selection/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/no-selection/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    System.add("system/cut/no-selection/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/no-selection/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/no-selection/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/no-selection/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/no-selection/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/no-selection/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/no-selection/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/no-selection/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/no-selection/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/no-selection/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/no-selection/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/no-selection/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/no-selection/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/no-selection/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    // 0n0n
    System.add("system/cut/0n0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/0n0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/0n0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/0n0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/0n0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/0n0n/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    System.add("system/cut/0n0n/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0n0n/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0n0n/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0n0n/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0n0n/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/0n0n/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/0n0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/0n0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/0n0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/0n0n/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/0n0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/0n0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/0n0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/0n0n/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    // Xn0n
    System.add("system/cut/Xn0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/Xn0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/Xn0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/Xn0n/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/Xn0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/Xn0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/Xn0n/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/Xn0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/Xn0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/Xn0n/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    // 0nXn
    System.add("system/cut/0nXn/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/0nXn/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/0nXn/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/0nXn/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/0nXn/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/0nXn/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0nXn/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0nXn/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/0nXn/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/0nXn/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/0nXn/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/0nXn/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/0nXn/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/0nXn/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/0nXn/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/0nXn/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/0nXn/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    // XnXn
    System.add("system/cut/XnXn/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/XnXn/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/XnXn/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/XnXn/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/XnXn/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/XnXn/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/XnXn/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/XnXn/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/XnXn/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/XnXn/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/XnXn/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/XnXn/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/XnXn/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    System.add("system/cut/XnXn/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    System.add("system/cut/XnXn/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    System.add("system/cut/XnXn/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPhi", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    System.add("system/cut/XnXn/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    // MC
    MC.add("MC/collisions/hPosXY", ";x (cm);y (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    MC.add("MC/collisions/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    MC.add("MC/collisions/hNPions", ";number of pions;counts", kTH1D, {{11, -0.5, 10.5}});
    MC.add("MC/collisions/hNumOfCollisionRecos", ";number of collision reconstructions;counts", kTH1D, {{11, -0.5, 10.5}});

    MC.add("MC/tracks/all/hPdgCode", ";pdg code;counts", kTH1D, {{2001, -1000.5, 1000.5}});
    MC.add("MC/tracks/all/hProducedByGenerator", ";produced by generator;counts", kTH1D, {{2, -0.5, 1.5}});
    MC.add("MC/tracks/all/hIsPhysicalPrimary", ";is physical primary;counts", kTH1D, {{2, -0.5, 1.5}});
    MC.add("MC/tracks/all/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    MC.add("MC/tracks/all/hEta", ";#eta;counts", kTH1D, {etaAxis});
    MC.add("MC/tracks/all/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    MC.add("MC/tracks/pions/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    MC.add("MC/tracks/pions/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaAxis, etaAxis});
    MC.add("MC/tracks/pions/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});

    MC.add("MC/system/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    MC.add("MC/system/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    MC.add("MC/system/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    MC.add("MC/system/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    MC.add("MC/system/hY", ";y;counts", kTH1D, {yAxis});
    MC.add("MC/system/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    MC.add("MC/system/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    MC.add("MC/system/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});

    // 4 pi QA
    if (do4pi) {
      FourPiQA.add("FourPiQA/reco/tracks/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
      FourPiQA.add("FourPiQA/reco/tracks/hEta", ";#eta;counts", kTH1D, {etaAxis});
      FourPiQA.add("FourPiQA/reco/tracks/hPhi", ";#phi;counts", kTH1D, {phiAxis});
      FourPiQA.add("FourPiQA/reco/system/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
      FourPiQA.add("FourPiQA/reco/system/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
      FourPiQA.add("FourPiQA/reco/system/hY", ";y;counts", kTH1D, {yAxis});
      FourPiQA.add("FourPiQA/reco/system/hPhi", ";#phi;counts", kTH1D, {phiAxis});
      FourPiQA.add("FourPiQA/MC/tracks/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
      FourPiQA.add("FourPiQA/MC/tracks/hEta", ";#eta;counts", kTH1D, {etaAxis});
      FourPiQA.add("FourPiQA/MC/tracks/hPhi", ";#phi;counts", kTH1D, {phiAxis});
      FourPiQA.add("FourPiQA/MC/system/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
      FourPiQA.add("FourPiQA/MC/system/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
      FourPiQA.add("FourPiQA/MC/system/hY", ";y;counts", kTH1D, {yAxis});
      FourPiQA.add("FourPiQA/MC/system/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    }
  }

  template <typename T>
  bool trackPassesCuts(const T& track) // track cuts (PID done separately)
  {
    if (!track.isPVContributor())
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 1);

    if (!track.hasITS())
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 2);

    if (track.itsNCls() < tracksMinItsNClsCut)
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 3);

    if (track.itsChi2NCl() > tracksMaxItsChi2NClCut)
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 4);

    if (!track.hasTPC())
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 5);

    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < tracksMinTpcNClsCut)
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 6);

    if (track.tpcChi2NCl() > tracksMaxTpcChi2NClCut || track.tpcChi2NCl() < tracksMinTpcChi2NClCut)
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 7);

    if (track.tpcNClsCrossedRows() < tracksMinTpcNClsCrossedRowsCut)
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 8);

    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < tracksMinTpcNClsCrossedOverFindableCut)
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 9);

    if (requireTof && !track.hasTOF())
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 10);

    if (track.pt() < tracksMinPtCut)
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 11);

    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > (0.0105 + 0.0350 / std::pow(track.pt(), 1.01)))
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 12);

    if (std::abs(eta(track.px(), track.py(), track.pz())) > PcEtaCut)
      return false;
    QC.fill(HIST("QC/tracks/hSelectionCounter"), 13);
    // if all selections passed
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

  bool systemPassCuts(const TLorentzVector& system) // system cuts
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

  double getPhiRandom(const std::vector<TLorentzVector>& cutTracks4Vecs) // decay phi anisotropy
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
  double getPhiCharge(const T& cutTracks, const std::vector<TLorentzVector>& cutTracks4Vecs)
  { // two possible definitions of phi: charge-based assignment
    TLorentzVector pOne, pTwo;
    pOne = (cutTracks[0].sign() > 0) ? cutTracks4Vecs[0] : cutTracks4Vecs[1];
    pTwo = (cutTracks[0].sign() > 0) ? cutTracks4Vecs[1] : cutTracks4Vecs[0];
    TLorentzVector pPlus = pOne + pTwo;
    TLorentzVector pMinus = pOne - pTwo;
    return pPlus.DeltaPhi(pMinus);
  }

  template <typename T>
  double getPhiChargeMC(const T& cutTracks, const std::vector<TLorentzVector>& cutTracks4Vecs)
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
    // QC histograms
    QC.fill(HIST("QC/collisions/all/hPosXY"), collision.posX(), collision.posY());
    QC.fill(HIST("QC/collisions/all/hPosZ"), collision.posZ());
    QC.fill(HIST("QC/collisions/all/hZdcCommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    QC.fill(HIST("QC/collisions/all/hZdcTime"), collision.timeZNA(), collision.timeZNC());
    QC.fill(HIST("QC/collisions/all/hNumContrib"), collision.numContrib());
    QC.fill(HIST("QC/collisions/all/hTotalFT0AmplitudeA"), collision.totalFT0AmplitudeA());
    QC.fill(HIST("QC/collisions/all/hTotalFT0AmplitudeC"), collision.totalFT0AmplitudeC());
    QC.fill(HIST("QC/collisions/all/hTotalFV0AmplitudeA"), collision.totalFV0AmplitudeA());
    QC.fill(HIST("QC/collisions/all/hTotalFDDAmplitudeA"), collision.totalFDDAmplitudeA());
    QC.fill(HIST("QC/collisions/all/hTotalFDDAmplitudeC"), collision.totalFDDAmplitudeC());
    QC.fill(HIST("QC/collisions/all/hTimeFT0A"), collision.timeFT0A());
    QC.fill(HIST("QC/collisions/all/hTimeFT0C"), collision.timeFT0C());
    QC.fill(HIST("QC/collisions/all/hTimeFV0A"), collision.timeFV0A());
    QC.fill(HIST("QC/collisions/all/hTimeFDDA"), collision.timeFDDA());
    QC.fill(HIST("QC/collisions/all/hTimeFDDC"), collision.timeFDDC());

    // vertex z-position cut
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut)
      return;

    // event tagging
    bool XnXn = false, OnOn = false, XnOn = false, OnXn = false; // note: On == 0n...
    int neutronClass = -1;
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() < ZNcommonEnergyCut) {
      OnOn = true;
      neutronClass = 0;
    }
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() < ZNcommonEnergyCut) {
      XnOn = true;
      neutronClass = 1;
    }
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) {
      OnXn = true;
      neutronClass = 2;
    }
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut &&
        collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) {
      XnXn = true;
      neutronClass = 3;
    }
    // vectors for storing selected tracks and their 4-vectors
    std::vector<decltype(tracks.begin())> cutTracks;
    std::vector<TLorentzVector> cutTracks4Vecs;

    for (const auto& track : tracks) {
      // double p = momentum(track.px(), track.py(), track.pz());
      QC.fill(HIST("QC/tracks/raw/hPt"), track.pt());
      QC.fill(HIST("QC/tracks/raw/hEta"), eta(track.px(), track.py(), track.pz()));
      QC.fill(HIST("QC/tracks/raw/hPhi"), phi(track.px(), track.py()));
      QC.fill(HIST("QC/tracks/raw/hTpcNSigmaPi"), track.tpcNSigmaPi());
      QC.fill(HIST("QC/tracks/raw/hTofNSigmaPi"), track.tofNSigmaPi());
      QC.fill(HIST("QC/tracks/raw/hTpcNSigmaEl"), track.tpcNSigmaEl());
      QC.fill(HIST("QC/tracks/raw/hDcaXYZ"), track.dcaZ(), track.dcaXY());
      QC.fill(HIST("QC/tracks/raw/hItsNCls"), track.itsNCls());
      QC.fill(HIST("QC/tracks/raw/hItsChi2NCl"), track.itsChi2NCl());
      QC.fill(HIST("QC/tracks/raw/hTpcChi2NCl"), track.tpcChi2NCl());
      QC.fill(HIST("QC/tracks/raw/hTpcNCls"), (track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()));
      QC.fill(HIST("QC/tracks/raw/hTpcNClsCrossedRows"), track.tpcNClsCrossedRows());
      QC.fill(HIST("QC/tracks/hSelectionCounter"), 0);

      if (!trackPassesCuts(track))
        continue;

      cutTracks.push_back(track);
      TLorentzVector track4Vec;
      track4Vec.SetXYZM(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged); // apriori assume pion mass
      cutTracks4Vecs.push_back(track4Vec);
      QC.fill(HIST("QC/tracks/cut/hTpcSignalVsP"), momentum(track.px(), track.py(), track.pz()), track.tpcSignal());
      QC.fill(HIST("QC/tracks/cut/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
      QC.fill(HIST("QC/tracks/cut/hDcaXYZ"), track.dcaZ(), track.dcaXY());
    }
    QC.fill(HIST("QC/tracks/cut/hRemainingTracks"), cutTracks.size());
    if (cutTracks.size() != cutTracks4Vecs.size()) { // sanity check
      LOG(error);
      return;
    }
    // reonstruct system and calculate total charge, save commonly used values into variables
    TLorentzVector system = reconstructSystem(cutTracks4Vecs);
    int totalCharge = tracksTotalCharge(cutTracks);
    double mass = system.M();
    double pT = system.Pt();
    double pTsquare = pT * pT;
    double rapidity = system.Rapidity();
    double systemPhi = system.Phi() + o2::constants::math::PI;

    if (do4pi && cutTracks.size() == 4 && totalCharge == 0) {
      // fill out some 4pi QC histograms
      for (int i = 0; i < static_cast<int>(cutTracks.size()); i++) {
        FourPiQA.fill(HIST("FourPiQA/reco/tracks/hPt"), cutTracks[i].pt());
        FourPiQA.fill(HIST("FourPiQA/reco/tracks/hEta"), eta(cutTracks[i].px(), cutTracks[i].py(), cutTracks[i].pz()));
        FourPiQA.fill(HIST("FourPiQA/reco/tracks/hPhi"), phi(cutTracks[i].px(), cutTracks[i].py()));
      }
      FourPiQA.fill(HIST("FourPiQA/reco/system/hM"), mass);
      FourPiQA.fill(HIST("FourPiQA/reco/system/hPt"), pT);
      FourPiQA.fill(HIST("FourPiQA/reco/system/hY"), rapidity);
      FourPiQA.fill(HIST("FourPiQA/reco/system/hPhi"), systemPhi);
    }

    // further consider only two pion systems
    if (cutTracks.size() != 2)
      return;
    for (int i = 0; i < static_cast<int>(cutTracks.size()); i++)
      QC.fill(HIST("QC/tracks/hSelectionCounter"), 14);

    QC.fill(HIST("QC/tracks/cut/hTpcNSigmaPi2D"), cutTracks[0].tpcNSigmaPi(), cutTracks[1].tpcNSigmaPi());
    QC.fill(HIST("QC/tracks/cut/hTpcNSigmaEl2D"), cutTracks[0].tpcNSigmaEl(), cutTracks[1].tpcNSigmaEl());

    if (!tracksPassPiPID(cutTracks))
      return;
    for (int i = 0; i < static_cast<int>(cutTracks.size()); i++)
      QC.fill(HIST("QC/tracks/hSelectionCounter"), 15);

    double phiRandom = getPhiRandom(cutTracks4Vecs);
    double phiCharge = getPhiCharge(cutTracks, cutTracks4Vecs);

    // differentiate leading- and subleading-momentum tracks
    auto leadingMomentumTrack = momentum(cutTracks[0].px(), cutTracks[0].py(), cutTracks[0].pz()) > momentum(cutTracks[1].px(), cutTracks[1].py(), cutTracks[1].pz()) ? cutTracks[0] : cutTracks[1];
    auto subleadingMomentumTrack = (leadingMomentumTrack == cutTracks[0]) ? cutTracks[1] : cutTracks[0];
    double leadingPt = leadingMomentumTrack.pt();
    double subleadingPt = subleadingMomentumTrack.pt();
    double leadingEta = eta(leadingMomentumTrack.px(), leadingMomentumTrack.py(), leadingMomentumTrack.pz());
    double subleadingEta = eta(subleadingMomentumTrack.px(), subleadingMomentumTrack.py(), subleadingMomentumTrack.pz());
    double leadingPhi = phi(leadingMomentumTrack.px(), leadingMomentumTrack.py());
    double subleadingPhi = phi(subleadingMomentumTrack.px(), subleadingMomentumTrack.py());
    // fill TOF hit checker
    QC.fill(HIST("QC/tracks/hTofHitCheck"), leadingMomentumTrack.hasTOF(), subleadingMomentumTrack.hasTOF());

    // fill tree
    std::vector<int> trackSigns = {leadingMomentumTrack.sign(), subleadingMomentumTrack.sign()};
    std::vector<double> trackPts = {leadingPt, subleadingPt};
    std::vector<double> trackEtas = {leadingEta, subleadingEta};
    std::vector<double> trackPhis = {leadingPhi, subleadingPhi};
    std::vector<double> trackMs = {o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged};
    std::vector<double> trackPiPIDs = {leadingMomentumTrack.tpcNSigmaPi(), subleadingMomentumTrack.tpcNSigmaPi()};
    std::vector<double> trackElPIDs = {leadingMomentumTrack.tpcNSigmaEl(), subleadingMomentumTrack.tpcNSigmaEl()};
    std::vector<double> trackDcaXYs = {leadingMomentumTrack.dcaXY(), subleadingMomentumTrack.dcaXY()};
    std::vector<double> trackDcaZs = {leadingMomentumTrack.dcaZ(), subleadingMomentumTrack.dcaZ()};
    std::vector<double> trackTpcSignals = {leadingMomentumTrack.tpcSignal(), subleadingMomentumTrack.tpcSignal()};
    Tree(collision.runNumber(), collision.globalBC(), collision.numContrib(),
         collision.posX(), collision.posY(), collision.posZ(),
         collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
         collision.timeFT0A(), collision.timeFT0C(), collision.timeFV0A(), collision.timeFDDA(), collision.timeFDDC(),
         collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC(), neutronClass,
         totalCharge, pT, system.Eta(), system.Phi(), mass, phiRandom, phiCharge,
         trackSigns, trackPts, trackEtas, trackPhis, trackMs, trackPiPIDs, trackElPIDs, trackDcaXYs, trackDcaZs, trackTpcSignals);
    // fill raw histograms according to the total charge
    switch (totalCharge) {
      case 0:
        Pions.fill(HIST("pions/no-selection/unlike-sign/hPt"), leadingPt, subleadingPt);
        Pions.fill(HIST("pions/no-selection/unlike-sign/hEta"), leadingEta, subleadingEta);
        Pions.fill(HIST("pions/no-selection/unlike-sign/hPhi"), leadingPhi, subleadingPhi);
        System.fill(HIST("system/raw/unlike-sign/hM"), mass);
        System.fill(HIST("system/raw/unlike-sign/hPt"), pT);
        System.fill(HIST("system/raw/unlike-sign/hPtVsM"), mass, pT);
        System.fill(HIST("system/raw/unlike-sign/hY"), rapidity);
        System.fill(HIST("system/raw/unlike-sign/hPhi"), systemPhi);
        break;

      case 2:
        Pions.fill(HIST("pions/no-selection/like-sign/hPt"), leadingPt, subleadingPt);
        Pions.fill(HIST("pions/no-selection/like-sign/hEta"), leadingEta, subleadingEta);
        Pions.fill(HIST("pions/no-selection/like-sign/hPhi"), leadingPhi, subleadingPhi);
        System.fill(HIST("system/raw/like-sign/positive/hM"), mass);
        System.fill(HIST("system/raw/like-sign/positive/hPt"), pT);
        System.fill(HIST("system/raw/like-sign/positive/hPtVsM"), mass, pT);
        System.fill(HIST("system/raw/like-sign/positive/hY"), rapidity);
        System.fill(HIST("system/raw/like-sign/positive/hPhi"), systemPhi);
        break;

      case -2:
        Pions.fill(HIST("pions/no-selection/like-sign/hPt"), leadingPt, subleadingPt);
        Pions.fill(HIST("pions/no-selection/like-sign/hEta"), leadingEta, subleadingEta);
        Pions.fill(HIST("pions/no-selection/like-sign/hPhi"), leadingPhi, subleadingPhi);
        System.fill(HIST("system/raw/like-sign/negative/hM"), mass);
        System.fill(HIST("system/raw/like-sign/negative/hPt"), pT);
        System.fill(HIST("system/raw/like-sign/negative/hPtVsM"), mass, pT);
        System.fill(HIST("system/raw/like-sign/negative/hY"), rapidity);
        System.fill(HIST("system/raw/like-sign/negative/hPhi"), systemPhi);
        break;

      default:
        break;
    }

    // apply cuts to system
    if (!systemPassCuts(system))
      return;

    QC.fill(HIST("QC/collisions/selected/hPosXY"), collision.posX(), collision.posY());
    QC.fill(HIST("QC/collisions/selected/hPosZ"), collision.posZ());
    QC.fill(HIST("QC/collisions/selected/hZdcCommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    QC.fill(HIST("QC/collisions/selected/hZdcTime"), collision.timeZNA(), collision.timeZNC());
    QC.fill(HIST("QC/collisions/selected/hNumContrib"), collision.numContrib());
    QC.fill(HIST("QC/collisions/selected/hTotalFT0AmplitudeA"), collision.totalFT0AmplitudeA());
    QC.fill(HIST("QC/collisions/selected/hTotalFT0AmplitudeC"), collision.totalFT0AmplitudeC());
    QC.fill(HIST("QC/collisions/selected/hTotalFV0AmplitudeA"), collision.totalFV0AmplitudeA());
    QC.fill(HIST("QC/collisions/selected/hTotalFDDAmplitudeA"), collision.totalFDDAmplitudeA());
    QC.fill(HIST("QC/collisions/selected/hTotalFDDAmplitudeC"), collision.totalFDDAmplitudeC());
    QC.fill(HIST("QC/collisions/selected/hTimeFT0A"), collision.timeFT0A());
    QC.fill(HIST("QC/collisions/selected/hTimeFT0C"), collision.timeFT0C());
    QC.fill(HIST("QC/collisions/selected/hTimeFV0A"), collision.timeFV0A());
    QC.fill(HIST("QC/collisions/selected/hTimeFDDA"), collision.timeFDDA());
    QC.fill(HIST("QC/collisions/selected/hTimeFDDC"), collision.timeFDDC());

    // fill histograms for system passing cuts
    switch (totalCharge) {
      case 0:
        Pions.fill(HIST("pions/selected/unlike-sign/hPt"), leadingPt, subleadingPt);
        Pions.fill(HIST("pions/selected/unlike-sign/hEta"), leadingEta, subleadingEta);
        Pions.fill(HIST("pions/selected/unlike-sign/hPhi"), leadingPhi, subleadingPhi);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hM"), mass);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPt"), pT);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPt2"), pTsquare);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPtVsM"), mass, pT);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hY"), rapidity);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPhi"), systemPhi);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPhiRandom"), phiRandom);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPhiCharge"), phiCharge);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        System.fill(HIST("system/cut/no-selection/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        if (OnOn) {
          System.fill(HIST("system/cut/0n0n/unlike-sign/hM"), mass);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPt"), pT);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPt2"), pTsquare);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hY"), rapidity);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPhi"), systemPhi);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/0n0n/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnOn) {
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hM"), mass);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPt"), pT);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPt2"), pTsquare);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hY"), rapidity);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPhi"), systemPhi);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/Xn0n/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (OnXn) {
          System.fill(HIST("system/cut/0nXn/unlike-sign/hM"), mass);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPt"), pT);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPt2"), pTsquare);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hY"), rapidity);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPhi"), systemPhi);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/0nXn/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnXn) {
          System.fill(HIST("system/cut/XnXn/unlike-sign/hM"), mass);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPt"), pT);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPt2"), pTsquare);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hY"), rapidity);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPhi"), systemPhi);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/XnXn/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      case 2:
        Pions.fill(HIST("pions/selected/like-sign/hPt"), leadingPt, subleadingPt);
        Pions.fill(HIST("pions/selected/like-sign/hEta"), leadingEta, subleadingEta);
        Pions.fill(HIST("pions/selected/like-sign/hPhi"), leadingPhi, subleadingPhi);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hM"), mass);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPt"), pT);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPt2"), pTsquare);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPtVsM"), mass, pT);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hY"), rapidity);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPhi"), systemPhi);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiRandom"), phiRandom);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiCharge"), phiCharge);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        System.fill(HIST("system/cut/no-selection/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        if (OnOn) {
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hM"), mass);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPt"), pT);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPt2"), pTsquare);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hY"), rapidity);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPhi"), systemPhi);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/0n0n/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnOn) {
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hM"), mass);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPt"), pT);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPt2"), pTsquare);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hY"), rapidity);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhi"), systemPhi);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/Xn0n/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (OnXn) {
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hM"), mass);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPt"), pT);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPt2"), pTsquare);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hY"), rapidity);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPhi"), systemPhi);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/0nXn/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnXn) {
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hM"), mass);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPt"), pT);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPt2"), pTsquare);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hY"), rapidity);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPhi"), systemPhi);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/XnXn/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      case -2:
        Pions.fill(HIST("pions/selected/like-sign/hPt"), leadingPt, subleadingPt);
        Pions.fill(HIST("pions/selected/like-sign/hEta"), leadingEta, subleadingEta);
        Pions.fill(HIST("pions/selected/like-sign/hPhi"), leadingPhi, subleadingPhi);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hM"), mass);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPt"), pT);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPt2"), pTsquare);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPtVsM"), mass, pT);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hY"), rapidity);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPhi"), systemPhi);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiRandom"), phiRandom);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiCharge"), phiCharge);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        System.fill(HIST("system/cut/no-selection/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        if (OnOn) {
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hM"), mass);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPt"), pT);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPt2"), pTsquare);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hY"), rapidity);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPhi"), systemPhi);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/0n0n/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnOn) {
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hM"), mass);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPt"), pT);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPt2"), pTsquare);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hY"), rapidity);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhi"), systemPhi);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/Xn0n/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (OnXn) {
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hM"), mass);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPt"), pT);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPt2"), pTsquare);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hY"), rapidity);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPhi"), systemPhi);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/0nXn/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnXn) {
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hM"), mass);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPt"), pT);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPt2"), pTsquare);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPtVsM"), mass, pT);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hY"), rapidity);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPhi"), systemPhi);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiRandom"), phiRandom);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiCharge"), phiCharge);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          System.fill(HIST("system/cut/XnXn/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      default:
        break;
    }
  }

  template <typename C, typename T>
  void processMC(C const& mcCollision, T const& mcParticles)
  {
    MC.fill(HIST("MC/collisions/hPosXY"), mcCollision.posX(), mcCollision.posY());
    MC.fill(HIST("MC/collisions/hPosZ"), mcCollision.posZ());

    std::vector<decltype(mcParticles.begin())> cutMcParticles;
    std::vector<TLorentzVector> mcParticles4Vecs;

    for (auto const& mcParticle : mcParticles) {
      MC.fill(HIST("MC/tracks/all/hPdgCode"), mcParticle.pdgCode());
      MC.fill(HIST("MC/tracks/all/hProducedByGenerator"), mcParticle.producedByGenerator());
      MC.fill(HIST("MC/tracks/all/hIsPhysicalPrimary"), mcParticle.isPhysicalPrimary());
      MC.fill(HIST("MC/tracks/all/hPt"), pt(mcParticle.px(), mcParticle.py()));
      MC.fill(HIST("MC/tracks/all/hEta"), eta(mcParticle.px(), mcParticle.py(), mcParticle.pz()));
      MC.fill(HIST("MC/tracks/all/hPhi"), phi(mcParticle.px(), mcParticle.py()));
      if (!mcParticle.isPhysicalPrimary() || std::abs(mcParticle.pdgCode()) != 211)
        continue;
      cutMcParticles.push_back(mcParticle);
      TLorentzVector pion4Vec;
      pion4Vec.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());
      mcParticles4Vecs.push_back(pion4Vec);
    }
    MC.fill(HIST("MC/collisions/hNPions"), cutMcParticles.size());

    if (mcParticles4Vecs.size() != cutMcParticles.size())
      return;
    if (tracksTotalChargeMC(cutMcParticles) != 0) // shouldn't happen in theory
      return;

    TLorentzVector system = reconstructSystem(mcParticles4Vecs);
    double mass = system.M();
    double pT = system.Pt();
    double pTsquare = pT * pT;
    double rapidity = system.Rapidity();
    double systemPhi = system.Phi() + o2::constants::math::PI;

    if (do4pi && cutMcParticles.size() == 4) {
      for (int i = 0; i < static_cast<int>(cutMcParticles.size()); i++) {
        FourPiQA.fill(HIST("FourPiQA/MC/tracks/hPt"), pt(cutMcParticles[i].px(), cutMcParticles[i].py()));
        FourPiQA.fill(HIST("FourPiQA/MC/tracks/hEta"), eta(cutMcParticles[i].px(), cutMcParticles[i].py(), cutMcParticles[i].pz()));
        FourPiQA.fill(HIST("FourPiQA/MC/tracks/hPhi"), phi(cutMcParticles[i].px(), cutMcParticles[i].py()));
      }
      FourPiQA.fill(HIST("FourPiQA/MC/system/hM"), mass);
      FourPiQA.fill(HIST("FourPiQA/MC/system/hPt"), pT);
      FourPiQA.fill(HIST("FourPiQA/MC/system/hY"), rapidity);
      FourPiQA.fill(HIST("FourPiQA/MC/system/hPhi"), systemPhi);
    }

    if (cutMcParticles.size() != 2)
      return;

    auto leadingMomentumPion = momentum(cutMcParticles[0].px(), cutMcParticles[0].py(), cutMcParticles[0].pz()) > momentum(cutMcParticles[1].px(), cutMcParticles[1].py(), cutMcParticles[1].pz()) ? cutMcParticles[0] : cutMcParticles[1];
    auto subleadingMomentumPion = (leadingMomentumPion == cutMcParticles[0]) ? cutMcParticles[1] : cutMcParticles[0];
    MC.fill(HIST("MC/tracks/pions/hPt"), pt(leadingMomentumPion.px(), leadingMomentumPion.py()), pt(subleadingMomentumPion.px(), subleadingMomentumPion.py()));
    MC.fill(HIST("MC/tracks/pions/hEta"), eta(leadingMomentumPion.px(), leadingMomentumPion.py(), leadingMomentumPion.pz()), eta(subleadingMomentumPion.px(), subleadingMomentumPion.py(), subleadingMomentumPion.pz()));
    MC.fill(HIST("MC/tracks/pions/hPhi"), phi(leadingMomentumPion.px(), leadingMomentumPion.py()), phi(subleadingMomentumPion.px(), subleadingMomentumPion.py()));

    double phiRandom = getPhiRandom(mcParticles4Vecs);
    double phiCharge = getPhiChargeMC(cutMcParticles, mcParticles4Vecs);

    MC.fill(HIST("MC/system/hM"), mass);
    MC.fill(HIST("MC/system/hPt"), pT);
    MC.fill(HIST("MC/system/hPtVsM"), mass, pT);
    MC.fill(HIST("MC/system/hPt2"), pTsquare);
    MC.fill(HIST("MC/system/hY"), rapidity);
    MC.fill(HIST("MC/system/hPhi"), systemPhi);
    MC.fill(HIST("MC/system/hPhiRandom"), phiRandom);
    MC.fill(HIST("MC/system/hPhiCharge"), phiCharge);
  }

  template <typename C>
  void checkNumberOfCollisionReconstructions(C const& collisions)
  {
    MC.fill(HIST("MC/collisions/hNumOfCollisionRecos"), collisions.size());
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
