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
/// \brief  Task for analysis of rho in UPCs using UD tables (from SG producer).
///         Includes event tagging based on ZN information, track selection, reconstruction,
///         and also some basic stuff for decay phi anisotropy studies
/// \author Jakub Juracka, jakub.juracka@cern.ch
/// \file   upcRhoAnalysis.cxx

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TPDGCode.h"

#include <random>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUdSgCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::UDZdcsReduced, aod::SGCollisions>::iterator;
using FullUdDgCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::UDZdcsReduced>::iterator;
using FullUdTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;
using FullMcUdCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDMcCollsLabels>::iterator;

namespace o2::aod
{
namespace reco_tree
{
// event info
DECLARE_SOA_COLUMN(RecoSetting, recoSetting, uint16_t);
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
DECLARE_SOA_COLUMN(NeutronClass, neutronClass, int);
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
                  reco_tree::RecoSetting, reco_tree::RunNumber, reco_tree::LocalBC, reco_tree::NumContrib, reco_tree::PosX, reco_tree::PosY, reco_tree::PosZ,
                  reco_tree::TotalFT0AmplitudeA, reco_tree::TotalFT0AmplitudeC, reco_tree::TotalFV0AmplitudeA, reco_tree::TotalFDDAmplitudeA, reco_tree::TotalFDDAmplitudeC,
                  reco_tree::TimeFT0A, reco_tree::TimeFT0C, reco_tree::TimeFV0A, reco_tree::TimeFDDA, reco_tree::TimeFDDC,
                  reco_tree::EnergyCommonZNA, reco_tree::EnergyCommonZNC, reco_tree::TimeZNA, reco_tree::TimeZNC, reco_tree::NeutronClass,
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

  SGSelector sgSelector;

  float pcEtaCut = 0.9; // physics coordination recommendation

  Configurable<int> numPions{"numPions", 2, "required number of pions in the event"};

  Configurable<bool> cutGapSide{"cutGapSide", true, "apply gap side cut"};
  Configurable<int> gapSide{"gapSide", 2, "required gap side"};
  Configurable<bool> useTrueGap{"useTrueGap", false, "use true gap"};
  Configurable<float> cutTrueGapSideFV0{"cutTrueGapSideFV0", 180000, "FV0A threshold for SG selector"};
  Configurable<float> cutTrueGapSideFT0A{"cutTrueGapSideFT0A", 150., "FT0A threshold for SG selector"};
  Configurable<float> cutTrueGapSideFT0C{"cutTrueGapSideFT0C", 50., "FT0C threshold for SG selector"};
  Configurable<float> cutTrueGapSideZDC{"cutTrueGapSideZDC", 10000., "ZDC threshold for SG selector. 0 is <1n, 4.2 is <2n, 6.7 is <3n, 9.5 is <4n, 12.5 is <5n"};

  Configurable<bool> requireTof{"requireTof", false, "require TOF signal"};
  Configurable<bool> onlyGoldenRuns{"onlyGoldenRuns", false, "process only golden runs"};
  Configurable<bool> useRecoFlag{"useRecoFlag", false, "use reco flag for event selection"};
  Configurable<int> cutRecoFlag{"cutRecoFlag", 1, "0 = std mode, 1 = upc mode"};
  Configurable<bool> useRctFlag{"useRctFlag", false, "use RCT flags for event selection"};
  Configurable<int> cutRctFlag{"cutRctFlag", 0, "0 = off, 1 = CBT, 2 = CBT+ZDC, 3 = CBThadron, 4 = CBThadron+ZDC"};

  Configurable<float> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "max Z position cut on collisions"};
  Configurable<bool> cutNumContribs{"cutNumContribs", true, "cut on number of contributors"};
  Configurable<int> collisionsNumContribsMaxCut{"collisionsNumContribsMaxCut", 2, "max number of contributors cut on collisions"};
  Configurable<float> znCommonEnergyCut{"znCommonEnergyCut", 0.0, "ZN common energy cut"};
  Configurable<float> znTimeCut{"znTimeCut", 2.0, "ZN time cut"};

  Configurable<float> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "TPC nSigma pion cut"};
  Configurable<float> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "max DCA cut on tracks"};
  Configurable<int> tracksMinItsNClsCut{"tracksMinItsNClsCut", 4, "min ITS clusters cut"};
  Configurable<float> tracksMaxItsChi2NClCut{"tracksMaxItsChi2NClCut", 3.0, "max ITS chi2/Ncls cut"};
  Configurable<int> tracksMinTpcNClsCut{"tracksMinTpcNClsCut", 120, "min TPC clusters cut"};
  Configurable<int> tracksMinTpcNClsCrossedRowsCut{"tracksMinTpcNClsCrossedRowsCut", 130, "min TPC crossed rows cut"};
  Configurable<float> tracksMinTpcChi2NClCut{"tracksMinTpcChi2NClCut", 1.0, "min TPC chi2/Ncls cut"};
  Configurable<float> tracksMaxTpcChi2NClCut{"tracksMaxTpcChi2NClCut", 3.0, "max TPC chi2/Ncls cut"};
  Configurable<float> tracksMinTpcNClsCrossedOverFindableCut{"tracksMinTpcNClsCrossedOverFindableCut", 1.0, "min TPC crossed rows / findable clusters cut"};
  Configurable<float> tracksMinPtCut{"tracksMinPtCut", 0.1, "min pT cut on tracks"};

  Configurable<float> systemMassMinCut{"systemMassMinCut", 0.5, "min M cut for reco system"};
  Configurable<float> systemMassMaxCut{"systemMassMaxCut", 1.0, "max M cut for reco system"};
  Configurable<float> systemPtCut{"systemPtCut", 0.1, "max pT cut for reco system"};
  Configurable<float> systemYCut{"systemYCut", 0.9, "rapiditiy cut for reco system"};

  ConfigurableAxis mAxis{"mAxis", {1000, 0.0, 10.0}, "#it{m} (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {1000, 0.0, 10.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {1000, 0.0, 0.1}, "#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis etaAxis{"etaAxis", {300, -1.5, 1.5}, "#it{#eta}"};
  ConfigurableAxis yAxis{"yAxis", {400, -4.0, 4.0}, "#it{y}"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, o2::constants::math::TwoPI}, "#it{#phi} (rad)"};
  ConfigurableAxis deltaPhiAxis{"deltaPhiAxis", {182, -o2::constants::math::PI, o2::constants::math::PI}, "#Delta#it{#phi} (rad)"};
  ConfigurableAxis znCommonEnergyAxis{"znCommonEnergyAxis", {250, -5.0, 20.0}, "ZN common energy (TeV)"};
  ConfigurableAxis znTimeAxis{"znTimeAxis", {200, -10.0, 10.0}, "ZN time (ns)"};
  ConfigurableAxis runNumberAxis{"runNumberAxis", {1355, 544012.5, 545367.5}, "run number"};
  ConfigurableAxis nSigmaAxis{"nSigmaAxis", {600, -30.0, 30.0}, "TPC #it{n#sigma}"};

  HistogramRegistry rQC{"rQC", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rTracks{"rTracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rSystem{"rSystem", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rMC{"rMC", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // QA //
    // collisions
    rQC.add("QC/collisions/all/hPosXY", ";vertex #it{x} (cm);vertex #it{y} (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    rQC.add("QC/collisions/all/hPosZ", ";vertex #it{z} (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    rQC.add("QC/collisions/all/hNumContrib", ";number of PV contributors;counts", kTH1D, {{36, -0.5, 35.5}});
    rQC.add("QC/collisions/all/hZdcCommonEnergy", ";ZNA common energy (TeV);ZNC common energy (TeV);counts", kTH2D, {znCommonEnergyAxis, znCommonEnergyAxis});
    rQC.add("QC/collisions/all/hZdcTime", ";ZNA time (ns);ZNC time (ns);counts", kTH2D, {znTimeAxis, znTimeAxis});
    rQC.add("QC/collisions/all/hTotalFT0AmplitudeA", ";FT0A amplitude;counts", kTH1D, {{160, 0.0, 160.0}});
    rQC.add("QC/collisions/all/hTotalFT0AmplitudeC", ";FT0C amplitude;counts", kTH1D, {{160, 0.0, 160.0}});
    rQC.add("QC/collisions/all/hTotalFV0AmplitudeA", ";FV0A amplitude;counts", kTH1D, {{300, 0.0, 300.0}});
    rQC.add("QC/collisions/all/hTotalFDDAmplitudeA", ";FDDA amplitude;counts", kTH1D, {{160, 0.0, 160.0}});
    rQC.add("QC/collisions/all/hTotalFDDAmplitudeC", ";FDDC amplitude;counts", kTH1D, {{50, 0.0, 50.0}});
    rQC.add("QC/collisions/all/hTimeFT0A", ";FT0A time (ns);counts", kTH1D, {{500, -10.0, 40.0}});
    rQC.add("QC/collisions/all/hTimeFT0C", ";FT0C time (ns);counts", kTH1D, {{500, -10.0, 40.0}});
    rQC.add("QC/collisions/all/hTimeFV0A", ";FV0A time (ns);counts", kTH1D, {{500, -10.0, 40.0}});
    rQC.add("QC/collisions/all/hTimeFDDA", ";FDDA time (ns);counts", kTH1D, {{500, -10.0, 40.0}});
    rQC.add("QC/collisions/all/hTimeFDDC", ";FDDC time (ns);counts", kTH1D, {{500, -10.0, 40.0}});
    // events with selected rho candidates
    rQC.addClone("QC/collisions/all/", "QC/collisions/trackSelections/"); // clone "all" histograms as "selected"
    rQC.addClone("QC/collisions/all/", "QC/collisions/systemSelections/");

    // tracks
    rQC.add("QC/tracks/all/hTpcNSigmaPi", ";TPC #it{n#sigma}(#pi);counts", kTH1D, {nSigmaAxis});
    rQC.add("QC/tracks/all/hTpcNSigmaEl", ";TPC #it{n#sigma}(e);counts", kTH1D, {nSigmaAxis});
    rQC.add("QC/tracks/all/hTpcNSigmaKa", ";TPC #it{n#sigma}(K);counts", kTH1D, {nSigmaAxis});
    rQC.add("QC/tracks/all/hDcaXYZ", ";track #it{DCA}_{z} (cm);track #it{DCA}_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {400, -2.0, 2.0}});
    rQC.add("QC/tracks/all/hItsNCls", ";ITS #it{N}_{cls};counts", kTH1D, {{11, -0.5, 10.5}});
    rQC.add("QC/tracks/all/hItsChi2NCl", ";ITS #it{#chi}^{2}/#it{N}_{cls};counts", kTH1D, {{200, 0.0, 20.0}});
    rQC.add("QC/tracks/all/hTpcChi2NCl", ";TPC #it{#chi}^{2}/#it{N}_{cls};counts", kTH1D, {{200, 0.0, 20.0}});
    rQC.add("QC/tracks/all/hTpcNCls", ";found TPC #it{N}_{cls};counts", kTH1D, {{160, 0.0, 160.0}}); // tpcNClsFindable() - track.tpcNClsFindableMinusFound
    rQC.add("QC/tracks/all/hTpcNClsCrossedRows", ";TPC crossed rows;counts", kTH1D, {{160, 0.0, 160.0}});
    rQC.add("QC/tracks/all/hTpcNClsCrossedRowsOverNClsFindable", ";TPC crossed rows/findable #it{N}_{cls};counts", kTH1D, {{300, 0.5, 2.5}});
    rQC.add("QC/tracks/all/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rQC.add("QC/tracks/all/hEta", ";#it{#eta};counts", kTH1D, {etaAxis});
    rQC.add("QC/tracks/all/hPhi", ";#it{#phi};counts", kTH1D, {phiAxis});
    rQC.add("QC/tracks/all/hTpcSignalVsP", ";|#it{p}| (GeV/#it{c});TPC d#it{E}/d#it{x} signal (arb. units);counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    rQC.add("QC/tracks/all/hTpcSignalVsPt", ";#it{p}_{T} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (arb. units);counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    // tracks passing selections
    rQC.addClone("QC/tracks/all/", "QC/tracks/trackSelections/"); // clone "raw" histograms as "cut"
    rQC.addClone("QC/tracks/all/", "QC/tracks/systemSelections/");
    rQC.add("QC/tracks/trackSelections/hRemainingTracks", ";remaining tracks;counts", kTH1D, {{21, -0.5, 20.5}});
    rQC.add("QC/tracks/trackSelections/hTpcNSigmaPi2D", ";TPC #it{n#sigma}(#pi)_{leading};TPC #it{n#sigma}(#pi)_{subleading};counts", kTH2D, {nSigmaAxis, nSigmaAxis});
    rQC.add("QC/tracks/trackSelections/hTpcNSigmaEl2D", ";TPC #it{n#sigma}(e)_{leading};TPC #it{n#sigma}(e)_{subleading};counts", kTH2D, {nSigmaAxis, nSigmaAxis});
    rQC.add("QC/tracks/trackSelections/hTpcNSigmaKa2D", ";TPC #it{n#sigma}(K)_{leading};TPC #it{n#sigma}(K)_{subleading};counts", kTH2D, {nSigmaAxis, nSigmaAxis});
    // selection counter
    std::vector<std::string> selectionCounterLabels = {"all tracks", "PV contributor", "ITS hit", "ITS #it{N}_{cls}", "itsClusterMap check", "ITS #it{#chi}^{2}/#it{N}_{cls}", "TPC hit", "found TPC #it{N}_{cls}", "TPC #it{#chi}^{2}/#it{N}_{cls}", "TPC crossed rows",
                                                       "TPC crossed rows/#it{N}_{cls}",
                                                       "TOF requirement",
                                                       "#it{p}_{T}", "#it{DCA}", "#it{#eta}", "exactly 2 tracks", "PID"};
    rQC.add("QC/tracks/hSelectionCounter", ";;tracks passing selections", kTH1D, {{static_cast<int>(selectionCounterLabels.size()), -0.5, static_cast<float>(selectionCounterLabels.size()) - 0.5}});
    rQC.add("QC/tracks/hSelectionCounterPerRun", ";;run number;tracks passing selections", kTH2D, {{static_cast<int>(selectionCounterLabels.size()), -0.5, static_cast<float>(selectionCounterLabels.size()) - 0.5}, runNumberAxis});
    for (int i = 0; i < static_cast<int>(selectionCounterLabels.size()); ++i) {
      rQC.get<TH1>(HIST("QC/tracks/hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, selectionCounterLabels[i].c_str());
      rQC.get<TH2>(HIST("QC/tracks/hSelectionCounterPerRun"))->GetXaxis()->SetBinLabel(i + 1, selectionCounterLabels[i].c_str());
    }
    rQC.add("QC/tracks/hTofHitCheck", ";leading track TOF hit;subleading track TOF hit;counts", kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
    rQC.get<TH2>(HIST("QC/tracks/hTofHitCheck"))->GetXaxis()->SetBinLabel(1, "no hit");
    rQC.get<TH2>(HIST("QC/tracks/hTofHitCheck"))->GetXaxis()->SetBinLabel(2, "hit");
    rQC.get<TH2>(HIST("QC/tracks/hTofHitCheck"))->GetYaxis()->SetBinLabel(1, "no hit");
    rQC.get<TH2>(HIST("QC/tracks/hTofHitCheck"))->GetYaxis()->SetBinLabel(2, "hit");

    // TRACKS (2D)
    rTracks.add("tracks/trackSelections/unlike-sign/hPt", ";#it{p}_{T leading} (GeV/#it{c});#it{p}_{T subleading} (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    rTracks.add("tracks/trackSelections/unlike-sign/hEta", ";#it{#eta}_{leading};#it{#eta}_{subleading};counts", kTH2D, {etaAxis, etaAxis});
    rTracks.add("tracks/trackSelections/unlike-sign/hPhi", ";#it{#phi}_{leading};#it{#phi}_{subleading};counts", kTH2D, {phiAxis, phiAxis});
    rTracks.addClone("tracks/trackSelections/unlike-sign/", "tracks/trackSelections/like-sign/positive/");
    rTracks.addClone("tracks/trackSelections/unlike-sign/", "tracks/trackSelections/like-sign/negative/");
    rTracks.addClone("tracks/trackSelections/", "tracks/systemSelections/");

    // SYSTEM
    rSystem.add("system/all/unlike-sign/hM", ";#it{m} (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    rSystem.add("system/all/unlike-sign/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rSystem.add("system/all/unlike-sign/hPt2", ";#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    rSystem.add("system/all/unlike-sign/hPtVsM", ";#it{m} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    rSystem.add("system/all/unlike-sign/hY", ";#it{y};counts", kTH1D, {yAxis});
    rSystem.add("system/all/unlike-sign/hPhi", ";#it{#phi};counts", kTH1D, {phiAxis});
    rSystem.add("system/all/unlike-sign/hPhiRandom", ";#Delta#it{#phi}_{random};counts", kTH1D, {deltaPhiAxis});
    rSystem.add("system/all/unlike-sign/hPhiCharge", ";#Delta#it{#phi}_{charge};counts", kTH1D, {deltaPhiAxis});
    rSystem.add("system/all/unlike-sign/hPhiRandomVsM", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi}_{random};counts", kTH2D, {mAxis, deltaPhiAxis});
    rSystem.add("system/all/unlike-sign/hPhiChargeVsM", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi}_{charge};counts", kTH2D, {mAxis, deltaPhiAxis});
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
    rMC.add("MC/collisions/hPosXY", ";vertex #it{x} (cm);vertex #it{y} (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    rMC.add("MC/collisions/hPosZ", ";vertex #it{z} (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    rMC.add("MC/collisions/hNPions", ";number of pions;counts", kTH1D, {{11, -0.5, 10.5}});
    rMC.add("MC/collisions/hNumOfCollisionRecos", ";number of collision reconstructions;counts", kTH1D, {{6, -0.5, 5.5}});
    rMC.add("MC/collisions/hRunNumberVsNumOfCollisionRecos", ";number of collision reconstructions;run number;counts", kTH2D, {{6, -0.5, 5.5}, runNumberAxis});
    // tracks
    rMC.add("MC/tracks/all/hPdgCode", ";pdg code;counts", kTH1D, {{2001, -1000.5, 1000.5}});
    rMC.add("MC/tracks/all/hProducedByGenerator", ";produced by generator;counts", kTH1D, {{2, -0.5, 1.5}});
    rMC.add("MC/tracks/all/hIsPhysicalPrimary", ";is physical primary;counts", kTH1D, {{2, -0.5, 1.5}});
    rMC.add("MC/tracks/all/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rMC.add("MC/tracks/all/hEta", ";#it{#eta};counts", kTH1D, {etaAxis});
    rMC.add("MC/tracks/all/hPhi", ";#it{#phi};counts", kTH1D, {phiAxis});
    rMC.add("MC/tracks/hPt", ";#it{p}_{T leading} (GeV/#it{c});#it{p}_{T subleading} (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    rMC.add("MC/tracks/hEta", ";#it{#eta}_{leading};#it{#eta}_{subleading};counts", kTH2D, {etaAxis, etaAxis});
    rMC.add("MC/tracks/hPhi", ";#it{#phi}_{leading};#it{#phi}_{subleading};counts", kTH2D, {phiAxis, phiAxis});
    // system
    rMC.add("MC/system/hM", ";#it{m} (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    rMC.add("MC/system/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    rMC.add("MC/system/hPt2", ";#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    rMC.add("MC/system/hPtVsM", ";#it{m} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    rMC.add("MC/system/hY", ";#it{y};counts", kTH1D, {yAxis});
    rMC.add("MC/system/hPhi", ";#it{#phi};counts", kTH1D, {phiAxis});
    rMC.add("MC/system/hPhiRandom", ";#Delta#it{#phi}_{random};counts", kTH1D, {deltaPhiAxis});
    rMC.add("MC/system/hPhiCharge", ";#Delta#it{#phi}_{charge};counts", kTH1D, {deltaPhiAxis});
    rMC.add("MC/system/hPhiRandomVsM", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi};counts", kTH2D, {mAxis, deltaPhiAxis});
    rMC.add("MC/system/hPhiChargeVsM", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi};counts", kTH2D, {mAxis, deltaPhiAxis});
    rMC.addClone("MC/system/", "MC/system/selected/");
  }

  std::unordered_set<int> goldenRuns = {544491, 544474, 544123, 544098, 544121, 544389, 544032, 544454, 544122,
                                        544510, 544476, 544091, 544095, 544490, 544124, 544508, 544391, 544013,
                                        544390, 544184, 544451, 544116, 544185, 544492, 544475, 544392, 544477, 544028};

  static constexpr std::string_view AppliedSelections[3] = {"all/", "trackSelections/", "systemSelections/"};
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
    rQC.fill(HIST("QC/tracks/") + HIST(AppliedSelections[cuts]) + HIST("hTpcSignalVsP"), std::abs(momentum(track.px(), track.py(), track.pz())), track.tpcSignal());
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
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPhiRandomVsM"), mass, phiRandom);
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPhiChargeVsM"), mass, phiCharge);
    } else {
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hM"), mass);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPt"), pt);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPt2"), pt * pt);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPtVsM"), mass, pt);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hY"), rapidity);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPhi"), phi);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPhiRandom"), phiRandom);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPhiCharge"), phiCharge);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPhiRandomVsM"), mass, phiRandom);
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPhiChargeVsM"), mass, phiCharge);
    }
  }

  bool cutItsLayers(uint8_t itsClusterMap) const
  {
    std::vector<std::pair<int8_t, std::array<uint8_t, 3>>> requiredITSHits{};
    requiredITSHits.push_back(std::make_pair(1, std::array<uint8_t, 3>{0, 1, 2})); // at least one hit in the innermost layer
    constexpr uint8_t kBit = 1;
    for (const auto& itsRequirement : requiredITSHits) {
      auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return itsClusterMap & (kBit << requiredLayer); });

      if ((itsRequirement.first == -1) && (hits > 0)) {
        return false; // no hits were required in specified layers
      } else if (hits < itsRequirement.first) {
        return false; // not enough hits found in specified layers
      }
    }
    return true;
  }

  template <typename C>
  bool isGoodRctFlag(const C& collision)
  {
    switch (cutRctFlag) {
      case 1:
        return sgSelector.isCBTOk(collision);
      case 2:
        return sgSelector.isCBTZdcOk(collision);
      case 3:
        return sgSelector.isCBTHadronOk(collision);
      case 4:
        return sgSelector.isCBTHadronZdcOk(collision);
      default:
        return true;
    }
  }

  template <typename C>
  bool collisionPassesCuts(const C& collision) // collision cuts
  {
    if (!collision.vtxITSTPC() || !collision.sbp() || !collision.itsROFb() || !collision.tfb()) // not applied automatically in 2023 Pb-Pb pass5
      return false;

    if (std::abs(collision.posZ()) > collisionsPosZMaxCut)
      return false;
    if (cutNumContribs && (collision.numContrib() > collisionsNumContribsMaxCut))
      return false;

    if (useRctFlag && !isGoodRctFlag(collision)) // check RCT flags
      return false;
    if (useRecoFlag && (collision.flags() != cutRecoFlag)) // check reconstruction mode
      return false;

    return true;
  }

  template <typename T, typename C>
  bool trackPassesCuts(const T& track, const C& collision) // track cuts (PID done separately)
  {
    if (!track.isPVContributor())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 1);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 1, collision.runNumber());

    if (!track.hasITS())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 2);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 2, collision.runNumber());

    if (track.itsNCls() < tracksMinItsNClsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 3);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 3, collision.runNumber());

    if (!cutItsLayers(track.itsClusterMap()))
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 4);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 4, collision.runNumber());

    if (track.itsChi2NCl() > tracksMaxItsChi2NClCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 5);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 5, collision.runNumber());

    if (!track.hasTPC())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 6);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 6, collision.runNumber());

    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < tracksMinTpcNClsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 7);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 7, collision.runNumber());

    if (track.tpcChi2NCl() > tracksMaxTpcChi2NClCut || track.tpcChi2NCl() < tracksMinTpcChi2NClCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 8);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 8, collision.runNumber());

    if (track.tpcNClsCrossedRows() < tracksMinTpcNClsCrossedRowsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 9);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 9, collision.runNumber());

    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < tracksMinTpcNClsCrossedOverFindableCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 10);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 10, collision.runNumber());

    if (requireTof && !track.hasTOF())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 11);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 11, collision.runNumber());

    if (track.pt() < tracksMinPtCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 12);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 12, collision.runNumber());

    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > (0.0105 + 0.0350 / std::pow(track.pt(), 1.01)))
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 13);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 13, collision.runNumber());

    if (std::abs(eta(track.px(), track.py(), track.pz())) > pcEtaCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 14);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 14, collision.runNumber());
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

  bool systemPassesCuts(const ROOT::Math::PxPyPzMVector& system) // system cuts
  {
    if (system.M() < systemMassMinCut || system.M() > systemMassMaxCut)
      return false;
    if (system.Pt() > systemPtCut)
      return false;
    if (std::abs(system.Rapidity()) > systemYCut)
      return false;
    return true;
  }

  ROOT::Math::PxPyPzMVector reconstructSystem(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracksLVs) // reconstruct system from 4-vectors
  {
    ROOT::Math::PxPyPzMVector system;
    for (const auto& trackLV : cutTracksLVs)
      system += trackLV;
    return system;
  }

  double deltaPhi(const ROOT::Math::PxPyPzMVector& p1, const ROOT::Math::PxPyPzMVector& p2)
  {
    double dPhi = p1.Phi() - p2.Phi();
    dPhi = std::fmod(dPhi + o2::constants::math::PI, o2::constants::math::TwoPI) - o2::constants::math::PI; // normalize to (-pi, pi)
    return dPhi;
  }

  float getPhiRandom(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracksLVs) // decay phi anisotropy
  {                                                                              // two possible definitions of phi: randomize the tracks
    int indices[2] = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();            // get time-based seed
    std::shuffle(std::begin(indices), std::end(indices), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    ROOT::Math::PxPyPzMVector pOne = cutTracksLVs[indices[0]];
    ROOT::Math::PxPyPzMVector pTwo = cutTracksLVs[indices[1]];
    ROOT::Math::PxPyPzMVector pPlus = pOne + pTwo;
    ROOT::Math::PxPyPzMVector pMinus = pOne - pTwo;
    return deltaPhi(pPlus, pMinus);
  }

  template <typename T>
  float getPhiCharge(const T& cutTracks, const std::vector<ROOT::Math::PxPyPzMVector>& cutTracksLVs)
  { // two possible definitions of phi: charge-based assignment
    ROOT::Math::PxPyPzMVector pOne, pTwo;
    pOne = (cutTracks[0].sign() > 0) ? cutTracksLVs[0] : cutTracksLVs[1];
    pTwo = (cutTracks[0].sign() > 0) ? cutTracksLVs[1] : cutTracksLVs[0];
    ROOT::Math::PxPyPzMVector pPlus = pOne + pTwo;
    ROOT::Math::PxPyPzMVector pMinus = pOne - pTwo;
    return deltaPhi(pPlus, pMinus);
  }

  template <typename T>
  float getPhiChargeMC(const T& cutTracks, const std::vector<ROOT::Math::PxPyPzMVector>& cutTracksLVs)
  { // the same as for data but using pdg code instead of charge
    ROOT::Math::PxPyPzMVector pOne, pTwo;
    pOne = (cutTracks[0].pdgCode() > 0) ? cutTracksLVs[0] : cutTracksLVs[1];
    pTwo = (cutTracks[0].pdgCode() > 0) ? cutTracksLVs[1] : cutTracksLVs[0];
    ROOT::Math::PxPyPzMVector pPlus = pOne + pTwo;
    ROOT::Math::PxPyPzMVector pMinus = pOne - pTwo;
    return deltaPhi(pPlus, pMinus);
  }

  template <typename C, typename T>
  void processReco(C const& collision, T const& tracks)
  {
    // check if the collision run number is contained within the goldenRuns set
    if (onlyGoldenRuns && !goldenRuns.contains(collision.runNumber()))
      return;

    fillCollisionQcHistos<0>(collision); // fill QC histograms before cuts
    if (!collisionPassesCuts(collision))
      return;

    int neutronClass = -1;
    bool xnxn = false, onon = false, xnon = false, onxn = false;
    float energyCommonZNA = collision.energyCommonZNA(), energyCommonZNC = collision.energyCommonZNC();
    float timeZNA = collision.timeZNA(), timeZNC = collision.timeZNC();
    if (std::isinf(energyCommonZNA))
      energyCommonZNA = -999;
    if (std::isinf(energyCommonZNC))
      energyCommonZNC = -999;
    if (std::isinf(timeZNA))
      timeZNA = -999;
    if (std::isinf(timeZNC))
      timeZNC = -999;

    if (energyCommonZNA <= znCommonEnergyCut && energyCommonZNC <= znCommonEnergyCut) {
      onon = true;
      neutronClass = 0;
    }
    if (energyCommonZNA > znCommonEnergyCut && std::abs(timeZNA) <= znTimeCut && energyCommonZNC <= znCommonEnergyCut) {
      xnon = true;
      neutronClass = 1;
    }
    if (energyCommonZNA <= znCommonEnergyCut && energyCommonZNC > znCommonEnergyCut && std::abs(timeZNC) <= znTimeCut) {
      onxn = true;
      neutronClass = 2;
    }
    if (energyCommonZNA > znCommonEnergyCut && std::abs(timeZNA) <= znTimeCut &&
        energyCommonZNC > znCommonEnergyCut && std::abs(timeZNC) <= znTimeCut) {
      xnxn = true;
      neutronClass = 3;
    }

    std::vector<decltype(tracks.begin())> cutTracks; // store selected tracks
    for (const auto& track : tracks) {
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 0);
      rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 0, collision.runNumber());
      fillTrackQcHistos<0>(track); // fill QC histograms before cuts

      if (!trackPassesCuts(track, collision)) // apply track cuts
        continue;
      cutTracks.push_back(track);
    }
    rQC.fill(HIST("QC/tracks/trackSelections/hRemainingTracks"), cutTracks.size());

    if (static_cast<int>(cutTracks.size()) != numPions) // further consider only two pion systems
      return;
    for (int i = 0; i < numPions; i++) {
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 15);
      rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 15, collision.runNumber());
    }
    rQC.fill(HIST("QC/tracks/trackSelections/hTpcNSigmaPi2D"), cutTracks[0].tpcNSigmaPi(), cutTracks[1].tpcNSigmaPi());
    rQC.fill(HIST("QC/tracks/trackSelections/hTpcNSigmaEl2D"), cutTracks[0].tpcNSigmaEl(), cutTracks[1].tpcNSigmaEl());
    rQC.fill(HIST("QC/tracks/trackSelections/hTpcNSigmaKa2D"), cutTracks[0].tpcNSigmaKa(), cutTracks[1].tpcNSigmaKa());

    // create a vector of 4-vectors for selected tracks
    std::vector<ROOT::Math::PxPyPzMVector> cutTracksLVs;
    for (const auto& cutTrack : cutTracks) {
      cutTracksLVs.push_back(ROOT::Math::PxPyPzMVector(cutTrack.px(), cutTrack.py(), cutTrack.pz(), o2::constants::physics::MassPionCharged)); // apriori assume pion mass
    }

    // differentiate leading- and subleading-momentum tracks
    auto leadingMomentumTrack = momentum(cutTracks[0].px(), cutTracks[0].py(), cutTracks[0].pz()) > momentum(cutTracks[1].px(), cutTracks[1].py(), cutTracks[1].pz()) ? cutTracks[0] : cutTracks[1];
    auto subleadingMomentumTrack = (leadingMomentumTrack == cutTracks[0]) ? cutTracks[1] : cutTracks[0];

    auto positiveTrack = cutTracks[0].sign() > 0 ? cutTracks[0] : cutTracks[1];
    auto negativeTrack = cutTracks[0].sign() > 0 ? cutTracks[1] : cutTracks[0];

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
    int trackSigns[2] = {positiveTrack.sign(), negativeTrack.sign()};
    float trackPts[2] = {positiveTrack.pt(), negativeTrack.pt()};
    float trackEtas[2] = {eta(positiveTrack.px(), positiveTrack.py(), positiveTrack.pz()), eta(negativeTrack.px(), negativeTrack.py(), negativeTrack.pz())};
    float trackPhis[2] = {phi(positiveTrack.px(), positiveTrack.py()), phi(negativeTrack.px(), negativeTrack.py())};
    float trackPiPIDs[2] = {positiveTrack.tpcNSigmaPi(), negativeTrack.tpcNSigmaPi()};
    float trackElPIDs[2] = {positiveTrack.tpcNSigmaEl(), negativeTrack.tpcNSigmaEl()};
    float trackKaPIDs[2] = {positiveTrack.tpcNSigmaKa(), negativeTrack.tpcNSigmaKa()};
    float trackDcaXYs[2] = {positiveTrack.dcaXY(), negativeTrack.dcaXY()};
    float trackDcaZs[2] = {positiveTrack.dcaZ(), negativeTrack.dcaZ()};
    float trackTpcSignals[2] = {positiveTrack.tpcSignal(), negativeTrack.tpcSignal()};
    recoTree(collision.flags(), collision.runNumber(), localBc, collision.numContrib(), collision.posX(), collision.posY(), collision.posZ(),
             collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
             collision.timeFT0A(), collision.timeFT0C(), collision.timeFV0A(), collision.timeFDDA(), collision.timeFDDC(),
             energyCommonZNA, energyCommonZNC, timeZNA, timeZNC, neutronClass,
             phiRandom, phiCharge, trackSigns, trackPts, trackEtas, trackPhis, trackPiPIDs, trackElPIDs, trackKaPIDs, trackDcaXYs, trackDcaZs, trackTpcSignals);

    if (!tracksPassPiPID(cutTracks)) // apply PID cut
      return;

    for (const auto& cutTrack : cutTracks) {
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 16);
      rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 16, collision.runNumber());
      fillTrackQcHistos<1>(cutTrack); // fill QC histograms after cuts
    }
    rQC.fill(HIST("QC/tracks/hTofHitCheck"), leadingMomentumTrack.hasTOF(), subleadingMomentumTrack.hasTOF());
    fillCollisionQcHistos<1>(collision); // fill QC histograms after track selections

    ROOT::Math::PxPyPzMVector system = reconstructSystem(cutTracksLVs);
    int totalCharge = tracksTotalCharge(cutTracks);
    float mass = system.M();
    float pT = system.Pt();
    float rapidity = system.Rapidity();
    float systemPhi = system.Phi() + o2::constants::math::PI;

    // fill raw histograms according to total charge
    switch (totalCharge) {
      case 0:
        fillTrack2dHistos<1, 0>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<0, 0, 0>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        break;

      case 2:
        fillTrack2dHistos<1, 1>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<0, 0, 1>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        break;

      case -2:
        fillTrack2dHistos<1, 2>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
        fillSystemHistos<0, 0, 2>(mass, pT, rapidity, systemPhi, phiRandom, phiCharge);
        break;

      default:
        break;
    }

    // apply cuts to system
    if (!systemPassesCuts(system))
      return;

    // fill histograms for system passing cuts
    switch (totalCharge) {
      case 0:
        fillCollisionQcHistos<2>(collision);
        for (const auto& cutTrack : cutTracks)
          fillTrackQcHistos<2>(cutTrack);
        fillTrack2dHistos<2, 0>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
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
        fillTrack2dHistos<2, 1>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
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
        fillTrack2dHistos<2, 2>(leadingPt, subleadingPt, leadingEta, subleadingEta, leadingPhi, subleadingPhi);
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
    std::vector<ROOT::Math::PxPyPzMVector> mcParticlesLVs;

    for (auto const& mcParticle : mcParticles) {
      rMC.fill(HIST("MC/tracks/all/hPdgCode"), mcParticle.pdgCode());
      rMC.fill(HIST("MC/tracks/all/hProducedByGenerator"), mcParticle.producedByGenerator());
      rMC.fill(HIST("MC/tracks/all/hIsPhysicalPrimary"), mcParticle.isPhysicalPrimary());
      rMC.fill(HIST("MC/tracks/all/hPt"), pt(mcParticle.px(), mcParticle.py()));
      rMC.fill(HIST("MC/tracks/all/hEta"), eta(mcParticle.px(), mcParticle.py(), mcParticle.pz()));
      rMC.fill(HIST("MC/tracks/all/hPhi"), phi(mcParticle.px(), mcParticle.py()));
      if (!mcParticle.isPhysicalPrimary() || std::abs(mcParticle.pdgCode()) != kPiPlus)
        continue;
      cutMcParticles.push_back(mcParticle);
      ROOT::Math::PxPyPzMVector pionLV;
      pionLV.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());
      mcParticlesLVs.push_back(pionLV);
    }
    rMC.fill(HIST("MC/collisions/hNPions"), cutMcParticles.size());

    if (static_cast<int>(cutMcParticles.size()) != numPions)
      return;
    if (mcParticlesLVs.size() != cutMcParticles.size())
      return;
    if (tracksTotalChargeMC(cutMcParticles) != 0) // shouldn't happen in theory
      return;

    ROOT::Math::PxPyPzMVector system = reconstructSystem(mcParticlesLVs);
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
    rMC.fill(HIST("MC/system/hPhiRandomVsM"), mass, phiRandom);
    rMC.fill(HIST("MC/system/hPhiChargeVsM"), mass, phiCharge);

    if (systemPassesCuts(system)) {
      rMC.fill(HIST("MC/system/selected/hM"), mass);
      rMC.fill(HIST("MC/system/selected/hPt"), pT);
      rMC.fill(HIST("MC/system/selected/hPtVsM"), mass, pT);
      rMC.fill(HIST("MC/system/selected/hPt2"), pT * pT);
      rMC.fill(HIST("MC/system/selected/hY"), rapidity);
      rMC.fill(HIST("MC/system/selected/hPhi"), systemPhi);
      rMC.fill(HIST("MC/system/selected/hPhiRandom"), phiRandom);
      rMC.fill(HIST("MC/system/selected/hPhiCharge"), phiCharge);
      rMC.fill(HIST("MC/system/selected/hPhiRandomVsM"), mass, phiRandom);
      rMC.fill(HIST("MC/system/selected/hPhiChargeVsM"), mass, phiCharge);
    }

    // fill mcTree
    auto positivePion = cutMcParticles[0].pdgCode() > 0 ? cutMcParticles[0] : cutMcParticles[1];
    auto negativePion = cutMcParticles[0].pdgCode() > 0 ? cutMcParticles[1] : cutMcParticles[0];
    int localBc = mcCollision.globalBC() % o2::constants::lhc::LHCMaxBunches;
    int trackSigns[2] = {positivePion.pdgCode() / std::abs(positivePion.pdgCode()), negativePion.pdgCode() / std::abs(negativePion.pdgCode())};
    float trackPts[2] = {pt(positivePion.px(), positivePion.py()), pt(negativePion.px(), negativePion.py())};
    float trackEtas[2] = {eta(positivePion.px(), positivePion.py(), positivePion.pz()), eta(negativePion.px(), negativePion.py(), negativePion.pz())};
    float trackPhis[2] = {phi(positivePion.px(), positivePion.py()), phi(negativePion.px(), negativePion.py())};
    mcTree(localBc,
           mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
           phiRandom, phiCharge, trackSigns, trackPts, trackEtas, trackPhis);
  }

  template <typename C>
  void checkNumberOfCollisionReconstructions(C const& collisions)
  {
    rMC.fill(HIST("MC/collisions/hNumOfCollisionRecos"), collisions.size());
    rMC.fill(HIST("MC/collisions/hRunNumberVsNumOfCollisionRecos"), collisions.size(), collisions.begin().runNumber());
  }

  void processSGdata(FullUdSgCollision const& collision, FullUdTracks const& tracks)
  {
    if (cutGapSide && collision.gapSide() != gapSide)
      return;
    if (useTrueGap && (collision.gapSide() != sgSelector.trueGap(collision, cutTrueGapSideFV0, cutTrueGapSideFT0A, cutTrueGapSideFT0C, cutTrueGapSideZDC))) // check true gap side
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

  void processCollisionRecoCheck(aod::UDMcCollision const& /* mcCollision */, soa::SmallGroups<soa::Join<aod::UDMcCollsLabels, aod::UDCollisions>> const& collisions)
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
