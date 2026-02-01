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
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(OccupancyInTime, occupancyInTime, float);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, float);
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
// tracks
DECLARE_SOA_COLUMN(LeadingTrackSign, leadingTrackSign, int);
DECLARE_SOA_COLUMN(SubleadingTrackSign, subleadingTrackSign, int);
DECLARE_SOA_COLUMN(LeadingTrackPt, leadingTrackPt, float);
DECLARE_SOA_COLUMN(SubleadingTrackPt, subleadingTrackPt, float);
DECLARE_SOA_COLUMN(LeadingTrackEta, leadingTrackEta, float);
DECLARE_SOA_COLUMN(SubleadingTrackEta, subleadingTrackEta, float);
DECLARE_SOA_COLUMN(LeadingTrackPhi, leadingTrackPhi, float);
DECLARE_SOA_COLUMN(SubleadingTrackPhi, subleadingTrackPhi, float);
DECLARE_SOA_COLUMN(LeadingTrackPiPID, leadingTrackPiPID, float);
DECLARE_SOA_COLUMN(SubleadingTrackPiPID, subleadingTrackPiPID, float);
DECLARE_SOA_COLUMN(LeadingTrackElPID, leadingTrackElPID, float);
DECLARE_SOA_COLUMN(SubleadingTrackElPID, subleadingTrackElPID, float);
DECLARE_SOA_COLUMN(LeadingTrackKaPID, leadingTrackKaPID, float);
DECLARE_SOA_COLUMN(SubleadingTrackKaPID, subleadingTrackKaPID, float);
DECLARE_SOA_COLUMN(LeadingTrackPrPID, leadingTrackPrPID, float);
DECLARE_SOA_COLUMN(SubleadingTrackPrPID, subleadingTrackPrPID, float);
} // namespace reco_tree
DECLARE_SOA_TABLE(RecoTree, "AOD", "RECOTREE",
                  reco_tree::RecoSetting, reco_tree::RunNumber, reco_tree::PosZ, reco_tree::OccupancyInTime, reco_tree::HadronicRate,
                  reco_tree::TotalFT0AmplitudeA, reco_tree::TotalFT0AmplitudeC, reco_tree::TotalFV0AmplitudeA, reco_tree::TotalFDDAmplitudeA, reco_tree::TotalFDDAmplitudeC,
                  reco_tree::TimeFT0A, reco_tree::TimeFT0C, reco_tree::TimeFV0A, reco_tree::TimeFDDA, reco_tree::TimeFDDC,
                  reco_tree::EnergyCommonZNA, reco_tree::EnergyCommonZNC, reco_tree::TimeZNA, reco_tree::TimeZNC, reco_tree::NeutronClass,
                  reco_tree::LeadingTrackSign, reco_tree::SubleadingTrackSign,
                  reco_tree::LeadingTrackPt, reco_tree::SubleadingTrackPt,
                  reco_tree::LeadingTrackEta, reco_tree::SubleadingTrackEta,
                  reco_tree::LeadingTrackPhi, reco_tree::SubleadingTrackPhi,
                  reco_tree::LeadingTrackPiPID, reco_tree::SubleadingTrackPiPID,
                  reco_tree::LeadingTrackElPID, reco_tree::SubleadingTrackElPID,
                  reco_tree::LeadingTrackKaPID, reco_tree::SubleadingTrackKaPID,
                  reco_tree::LeadingTrackPrPID, reco_tree::SubleadingTrackPrPID);

namespace mc_tree
{
// misc event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
// tracks
DECLARE_SOA_COLUMN(LeadingTrackSign, leadingTrackSign, int);
DECLARE_SOA_COLUMN(SubleadingTrackSign, subleadingTrackSign, int);
DECLARE_SOA_COLUMN(LeadingTrackPt, leadingTrackPt, float);
DECLARE_SOA_COLUMN(SubleadingTrackPt, subleadingTrackPt, float);
DECLARE_SOA_COLUMN(LeadingTrackEta, leadingTrackEta, float);
DECLARE_SOA_COLUMN(SubleadingTrackEta, subleadingTrackEta, float);
DECLARE_SOA_COLUMN(LeadingTrackPhi, leadingTrackPhi, float);
DECLARE_SOA_COLUMN(SubleadingTrackPhi, subleadingTrackPhi, float);
} // namespace mc_tree
DECLARE_SOA_TABLE(McTree, "AOD", "MCTREE",
                  mc_tree::RunNumber, mc_tree::PosZ,
                  mc_tree::LeadingTrackSign, mc_tree::SubleadingTrackSign,
                  mc_tree::LeadingTrackPt, mc_tree::SubleadingTrackPt,
                  mc_tree::LeadingTrackEta, mc_tree::SubleadingTrackEta,
                  mc_tree::LeadingTrackPhi, mc_tree::SubleadingTrackPhi);
} // namespace o2::aod

struct UpcRhoAnalysis {
  Produces<o2::aod::RecoTree> recoTree;
  Produces<o2::aod::McTree> mcTree;

  SGSelector sgSelector;

  const float pcEtaCut = 0.9; // physics coordination recommendation
  const std::vector<int> runNumbers = {544013, 544028, 544032, 544091, 544095, 544098, 544116, 544121, 544122, 544123, 544124, 544184, 544185, 544389, 544390, 544391, 544392, 544451, 544454, 544474, 544475, 544476, 544477, 544490, 544491, 544492, 544508, 544510, 544511, 544512, 544514, 544515, 544518, 544548, 544549, 544550, 544551, 544564, 544565, 544567, 544568, 544580, 544582, 544583, 544585, 544614, 544640, 544652, 544653, 544672, 544674, 544692, 544693, 544694, 544696, 544739, 544742, 544754, 544767, 544794, 544795, 544797, 544813, 544868, 544886, 544887, 544896, 544911, 544913, 544914, 544917, 544931, 544947, 544961, 544963, 544964, 544968, 544991, 544992, 545004, 545008, 545009, 545041, 545042, 545044, 545047, 545060, 545062, 545063, 545064, 545066, 545086, 545103, 545117, 545171, 545184, 545185, 545210, 545222, 545223, 545246, 545249, 545262, 545289, 545291, 545294, 545295, 545296, 545311, 545312, 545332, 545345, 545367};
  AxisSpec runNumberAxis = {static_cast<int>(runNumbers.size()), 0.5, static_cast<double>(runNumbers.size()) + 0.5, "run number"};

  Configurable<bool> isPO{"isPO", false, "process proton-oxygen data"};

  Configurable<bool> cutGapSide{"cutGapSide", true, "apply gap side cut"};
  Configurable<int> gapSide{"gapSide", 2, "required gap side"};
  Configurable<bool> useTrueGap{"useTrueGap", false, "use true gap"};
  Configurable<float> cutTrueGapSideFV0{"cutTrueGapSideFV0", 180000, "FV0A threshold for SG selector"};
  Configurable<float> cutTrueGapSideFT0A{"cutTrueGapSideFT0A", 150., "FT0A threshold for SG selector"};
  Configurable<float> cutTrueGapSideFT0C{"cutTrueGapSideFT0C", 50., "FT0C threshold for SG selector"};
  Configurable<float> cutTrueGapSideZDC{"cutTrueGapSideZDC", 10000., "ZDC threshold for SG selector. 0 is <1n, 4.2 is <2n, 6.7 is <3n, 9.5 is <4n, 12.5 is <5n"};

  Configurable<bool> requireTof{"requireTof", false, "require TOF signal"};
  Configurable<bool> useRecoFlag{"useRecoFlag", false, "use reco flag for event selection"};
  Configurable<int> cutRecoFlag{"cutRecoFlag", 1, "0 = std mode, 1 = upc mode"};
  Configurable<bool> useRctFlag{"useRctFlag", false, "use RCT flags for event selection"};
  Configurable<int> cutRctFlag{"cutRctFlag", 0, "0 = off, 1 = CBT, 2 = CBT+ZDC, 3 = CBThadron, 4 = CBThadron+ZDC"};

  Configurable<bool> selectRuns{"selectRuns", false, "select runs from the list"};
  Configurable<std::vector<int>> selectedRuns{"selectedRuns", {544013, 544028, 544032, 544091, 544095, 544098, 544116, 544121, 544122, 544123, 544124, 544184, 544185, 544389, 544390, 544391, 544392, 544451, 544454, 544474, 544475, 544476, 544477, 544490, 544491, 544492, 544508, 544510, 544511, 544512, 544514, 544515, 544518, 544548, 544549, 544550, 544551, 544564, 544565, 544567, 544568, 544580, 544582, 544583, 544585, 544614, 544640, 544652, 544653, 544672, 544674, 544692, 544693, 544694, 544696, 544739, 544742, 544754, 544767, 544794, 544795, 544797, 544813, 544868, 544886, 544887, 544896, 544913, 544914, 544917, 544931, 544947, 544961, 544963, 544964, 544968, 544992, 545009, 545044, 545047, 545063, 545064, 545066, 545185, 545210, 545223, 545249, 545291, 545294, 545295, 545296, 545312}, "list of selected runs"};

  Configurable<float> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "max Z position cut on collisions"};
  Configurable<bool> cutNumContribs{"cutNumContribs", true, "cut on number of contributors"};
  Configurable<int> collisionsNumContribsMaxCut{"collisionsNumContribsMaxCut", 2, "max number of contributors cut on collisions"};
  Configurable<float> znCommonEnergyCut{"znCommonEnergyCut", 0.0, "ZN common energy cut"};
  Configurable<float> znTimeCut{"znTimeCut", 2.0, "ZN time cut"};

  Configurable<float> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "TPC nSigma pion cut"};
  Configurable<bool> rejectLowerProbPairs{"rejectLowerProbPairs", true, "reject track pairs with lower El or Ka PID radii"};
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

  ConfigurableAxis mAxis{"mAxis", {400, 0.0, 4.0}, "#it{m} (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {400, 0.0, 4.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {1000, 0.0, 1.0}, "#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis etaAxis{"etaAxis", {300, -1.5, 1.5}, "#it{#eta}"};
  ConfigurableAxis yAxis{"yAxis", {300, -1.5, 1.5}, "#it{y}"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, o2::constants::math::TwoPI}, "#it{#phi} (rad)"};
  ConfigurableAxis deltaPhiAxis{"deltaPhiAxis", {182, -o2::constants::math::PI, o2::constants::math::PI}, "#Delta#it{#phi} (rad)"};
  ConfigurableAxis znCommonEnergyAxis{"znCommonEnergyAxis", {250, -5.0, 20.0}, "ZN common energy (TeV)"};
  ConfigurableAxis znTimeAxis{"znTimeAxis", {200, -10.0, 10.0}, "ZN time (ns)"};
  ConfigurableAxis nSigmaAxis{"nSigmaAxis", {600, -30.0, 30.0}, "TPC #it{n#sigma}"};

  HistogramRegistry rQC{"rQC", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rTracks{"rTracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rSystem{"rSystem", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rMC{"rMC", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry rResolution{"rResolution", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processSGdata") || context.mOptions.get<bool>("processDGdata")) {
      // QA //
      // collisions
      rQC.add("QC/collisions/all/hPosXY", ";vertex #it{x} (cm);vertex #it{y} (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
      rQC.add("QC/collisions/all/hPosZ", ";vertex #it{z} (cm);counts", kTH1D, {{400, -20.0, 20.0}});
      rQC.add("QC/collisions/all/hNumContrib", ";number of PV contributors;counts", kTH1D, {{36, -0.5, 35.5}});
      rQC.add("QC/collisions/all/hZdcCommonEnergy", ";ZNA common energy (TeV);ZNC common energy (TeV);counts", kTH2D, {znCommonEnergyAxis, znCommonEnergyAxis});
      rQC.add("QC/collisions/all/hZdcTime", ";ZNA time (ns);ZNC time (ns);counts", kTH2D, {znTimeAxis, znTimeAxis});
      rQC.add("QC/collisions/all/hZNTimeVsZNCommonEnergy", ";ZNA/C common energy (TeV);ZNA/C time (ns);counts", kTH2D, {znCommonEnergyAxis, znTimeAxis});
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
      rQC.add("QC/collisions/all/hOccupancyInTime", ";occupancy in time;counts", kTH1D, {{1100, -100.0, 1000.0}});
      // events with selected rho candidates
      rQC.addClone("QC/collisions/all/", "QC/collisions/trackSelections/");
      rQC.addClone("QC/collisions/all/", "QC/collisions/systemSelections/");

      std::vector<std::string> collisionSelectionCounterLabels = {"all collisions", "rapidity gap", "ITS-TPC vertex", "same bunch pile-up", "ITS ROF border", "TF border", "#it{z} position", "number of contributors", "RCT selections", "reco flag selection"};
      rQC.add("QC/collisions/hSelectionCounter", ";;collisions passing selections", kTH1D, {{static_cast<int>(collisionSelectionCounterLabels.size()), -0.5, static_cast<float>(collisionSelectionCounterLabels.size()) - 0.5}});
      rQC.add("QC/collisions/hSelectionCounterPerRun", ";;run number;collisions passing selections", kTH2D, {{static_cast<int>(collisionSelectionCounterLabels.size()), -0.5, static_cast<float>(collisionSelectionCounterLabels.size()) - 0.5}, runNumberAxis});
      for (int i = 0; i < static_cast<int>(collisionSelectionCounterLabels.size()); ++i) {
        rQC.get<TH1>(HIST("QC/collisions/hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, collisionSelectionCounterLabels[i].c_str());
        rQC.get<TH2>(HIST("QC/collisions/hSelectionCounterPerRun"))->GetXaxis()->SetBinLabel(i + 1, collisionSelectionCounterLabels[i].c_str());
      }
      // tracks
      rQC.add("QC/tracks/all/hTpcNSigmaPi", ";TPC #it{n#sigma}(#pi);counts", kTH1D, {nSigmaAxis});
      rQC.add("QC/tracks/all/hTpcNSigmaEl", ";TPC #it{n#sigma}(e);counts", kTH1D, {nSigmaAxis});
      rQC.add("QC/tracks/all/hTpcNSigmaKa", ";TPC #it{n#sigma}(K);counts", kTH1D, {nSigmaAxis});
      rQC.add("QC/tracks/all/hTpcNSigmaPr", ";TPC #it{n#sigma}(p);counts", kTH1D, {nSigmaAxis});
      rQC.add("QC/tracks/all/hDcaXYZ", ";track #it{DCA}_{z} (cm);track #it{DCA}_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {400, -2.0, 2.0}});
      rQC.add("QC/tracks/all/hItsNCls", ";ITS #it{N}_{cls};counts", kTH1D, {{11, -0.5, 10.5}});
      rQC.add("QC/tracks/all/hItsChi2NCl", ";ITS #it{#chi}^{2}/#it{N}_{cls};counts", kTH1D, {{200, 0.0, 20.0}});
      rQC.add("QC/tracks/all/hTpcChi2NCl", ";TPC #it{#chi}^{2}/#it{N}_{cls};counts", kTH1D, {{200, 0.0, 20.0}});
      rQC.add("QC/tracks/all/hTpcNCls", ";found TPC #it{N}_{cls};counts", kTH1D, {{160, 0.0, 160.0}}); // tpcNClsFindable() - track.tpcNClsFindableMinusFound
      rQC.add("QC/tracks/all/hTpcNClsCrossedRows", ";TPC crossed rows;counts", kTH1D, {{160, 0.0, 160.0}});
      rQC.add("QC/tracks/all/hTpcNClsCrossedRowsOverNClsFindable", ";TPC crossed rows/findable #it{N}_{cls};counts", kTH1D, {{300, 0.5, 2.5}});
      rQC.add("QC/tracks/all/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
      rQC.add("QC/tracks/all/hEta", ";#it{#eta};counts", kTH1D, {etaAxis});
      rQC.add("QC/tracks/all/hPhi", ";#it{#phi} (rad);counts", kTH1D, {phiAxis});
      rQC.add("QC/tracks/all/hTpcSignalVsP", ";|#it{p}| (GeV/#it{c});TPC d#it{E}/d#it{x} signal (arb. units);counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
      rQC.add("QC/tracks/all/hTpcSignalVsPt", ";#it{p}_{T} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (arb. units);counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
      // tracks passing selections
      rQC.addClone("QC/tracks/all/", "QC/tracks/trackSelections/");
      rQC.addClone("QC/tracks/all/", "QC/tracks/systemSelections/");
      rQC.add("QC/tracks/trackSelections/hRemainingTracks", ";remaining tracks;counts", kTH1D, {{21, -0.5, 20.5}});
      rQC.add("QC/tracks/trackSelections/hTpcNSigmaPi2D", ";TPC #it{n#sigma}(#pi)_{leading};TPC #it{n#sigma}(#pi)_{subleading};counts", kTH2D, {nSigmaAxis, nSigmaAxis});
      rQC.add("QC/tracks/trackSelections/hTpcNSigmaEl2D", ";TPC #it{n#sigma}(e)_{leading};TPC #it{n#sigma}(e)_{subleading};counts", kTH2D, {nSigmaAxis, nSigmaAxis});
      rQC.add("QC/tracks/trackSelections/hTpcNSigmaKa2D", ";TPC #it{n#sigma}(K)_{leading};TPC #it{n#sigma}(K)_{subleading};counts", kTH2D, {nSigmaAxis, nSigmaAxis});
      rQC.add("QC/tracks/trackSelections/hTpcNSigmaPr2D", ";TPC #it{n#sigma}(p)_{leading};TPC #it{n#sigma}(p)_{subleading};counts", kTH2D, {nSigmaAxis, nSigmaAxis});
      // selection counter
      std::vector<std::string> trackSelectionCounterLabels = {"all tracks", "PV contributor", "ITS hit", "ITS #it{N}_{cls}", "itsClusterMap check", "ITS #it{#chi}^{2}/#it{N}_{cls}", "TPC hit", "found TPC #it{N}_{cls}", "TPC #it{#chi}^{2}/#it{N}_{cls}", "TPC crossed rows",
                                                              "TPC crossed rows/#it{N}_{cls}",
                                                              "TOF requirement",
                                                              "#it{p}_{T}", "#it{DCA}", "#it{#eta}", "exactly 2 tracks", "PID"};
      rQC.add("QC/tracks/hSelectionCounter", ";;tracks passing selections", kTH1D, {{static_cast<int>(trackSelectionCounterLabels.size()), -0.5, static_cast<float>(trackSelectionCounterLabels.size()) - 0.5}});
      rQC.add("QC/tracks/hSelectionCounterPerRun", ";;run number;tracks passing selections", kTH2D, {{static_cast<int>(trackSelectionCounterLabels.size()), -0.5, static_cast<float>(trackSelectionCounterLabels.size()) - 0.5}, runNumberAxis});
      for (int i = 0; i < static_cast<int>(trackSelectionCounterLabels.size()); ++i) {
        rQC.get<TH1>(HIST("QC/tracks/hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, trackSelectionCounterLabels[i].c_str());
        rQC.get<TH2>(HIST("QC/tracks/hSelectionCounterPerRun"))->GetXaxis()->SetBinLabel(i + 1, trackSelectionCounterLabels[i].c_str());
      }
      for (int i = 0; i < static_cast<int>(runNumbers.size()); ++i)
        rQC.get<TH2>(HIST("QC/tracks/hSelectionCounterPerRun"))->GetYaxis()->SetBinLabel(i + 1, std::to_string(runNumbers[i]).c_str());
      
      rQC.add("QC/tracks/hTofHitCheck", ";leading track TOF hit;subleading track TOF hit;counts", kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
      rQC.get<TH2>(HIST("QC/tracks/hTofHitCheck"))->GetXaxis()->SetBinLabel(1, "no hit");
      rQC.get<TH2>(HIST("QC/tracks/hTofHitCheck"))->GetXaxis()->SetBinLabel(2, "hit");
      rQC.get<TH2>(HIST("QC/tracks/hTofHitCheck"))->GetYaxis()->SetBinLabel(1, "no hit");
      rQC.get<TH2>(HIST("QC/tracks/hTofHitCheck"))->GetYaxis()->SetBinLabel(2, "hit");
      // PID "radii" plots
      rQC.add("QC/tracks/hPiPIDRadius", ";#it{n#sigma}(#pi) radius;counts", kTH1D, {{1000, 0.0, 10.0}});
      rQC.add("QC/tracks/hElPIDRadius", ";#it{n#sigma}(e) radius;counts", kTH1D, {{1000, 0.0, 10.0}});
      rQC.add("QC/tracks/hKaPIDRadius", ";#it{n#sigma}(K) radius;counts", kTH1D, {{1000, 0.0, 10.0}});
      rQC.add("QC/tracks/hPrPIDRadius", ";#it{n#sigma}(p) radius;counts", kTH1D, {{1000, 0.0, 10.0}});

      // TRACKS (2D)
      rTracks.add("tracks/trackSelections/unlike-sign/hPt", ";#it{p}_{T leading} (GeV/#it{c});#it{p}_{T subleading} (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
      rTracks.add("tracks/trackSelections/unlike-sign/hEta", ";#it{#eta}_{leading};#it{#eta}_{subleading};counts", kTH2D, {etaAxis, etaAxis});
      rTracks.add("tracks/trackSelections/unlike-sign/hPhi", ";#it{#phi}_{leading} (rad);#it{#phi}_{subleading} (rad);counts", kTH2D, {phiAxis, phiAxis});
      rTracks.addClone("tracks/trackSelections/unlike-sign/", "tracks/trackSelections/like-sign/positive/");
      rTracks.addClone("tracks/trackSelections/unlike-sign/", "tracks/trackSelections/like-sign/negative/");
      rTracks.addClone("tracks/trackSelections/", "tracks/systemSelections/");

      // SYSTEM
      rSystem.add("system/all/unlike-sign/hM", ";#it{m} (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
      rSystem.add("system/all/unlike-sign/hRecoSettingVsM", ";#it{m} (GeV/#it{c}^{2});reco setting;counts", kTH2D, {mAxis, {2, -0.5, 1.5}});
      rSystem.add("system/all/unlike-sign/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
      rSystem.add("system/all/unlike-sign/hPt2", ";#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
      rSystem.add("system/all/unlike-sign/hPtVsM", ";#it{m} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
      rSystem.add("system/all/unlike-sign/hPt2VsM", ";#it{m} (GeV/#it{c}^{2});#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH2D, {mAxis, pt2Axis});
      rSystem.add("system/all/unlike-sign/hY", ";#it{y};counts", kTH1D, {yAxis});
      rSystem.add("system/all/unlike-sign/hPhi", ";#it{#phi} (rad);counts", kTH1D, {phiAxis});
      rSystem.add("system/all/unlike-sign/hPhiRandom", ";#Delta#it{#phi}_{random} (rad);counts", kTH1D, {deltaPhiAxis});
      rSystem.add("system/all/unlike-sign/hPhiCharge", ";#Delta#it{#phi}_{charge} (rad);counts", kTH1D, {deltaPhiAxis});
      rSystem.add("system/all/unlike-sign/hPhiRandomVsM", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi}_{random} (rad);counts", kTH2D, {mAxis, deltaPhiAxis});
      rSystem.add("system/all/unlike-sign/hPhiChargeVsM", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi}_{charge} (rad);counts", kTH2D, {mAxis, deltaPhiAxis});
      // clones for like-sign
      rSystem.addClone("system/all/unlike-sign/", "system/all/like-sign/positive/");
      rSystem.addClone("system/all/unlike-sign/", "system/all/like-sign/negative/");
      // selected rhos
      rSystem.addClone("system/all/", "system/selected/AnAn/");
      // clones for neutron classes
      rSystem.addClone("system/selected/AnAn/", "system/selected/0n0n/");
      rSystem.addClone("system/selected/AnAn/", "system/selected/Xn0n/");
      rSystem.addClone("system/selected/AnAn/", "system/selected/0nXn/");
      rSystem.addClone("system/selected/AnAn/", "system/selected/XnXn/");
    }

    if (context.mOptions.get<bool>("processMCdata") || context.mOptions.get<bool>("processMCdataWithBCs")) {
      // MC
      // collisions
      rMC.add("MC/collisions/hPosXY", ";vertex #it{x} (cm);vertex #it{y} (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
      rMC.add("MC/collisions/hPosZ", ";vertex #it{z} (cm);counts", kTH1D, {{400, -20.0, 20.0}});
      rMC.add("MC/collisions/hNPions", ";number of pions;counts", kTH1D, {{11, -0.5, 10.5}});
      // tracks
      rMC.add("MC/tracks/all/hPdgCode", ";pdg code;counts", kTH1D, {{2001, -1000.5, 1000.5}});
      rMC.add("MC/tracks/all/hMotherPdgCode", ";mother pdg code;counts", kTH1D, {{2001, -1000.5, 1000.5}});
      rMC.add("MC/tracks/all/hProducedByGenerator", ";produced by generator;counts", kTH1D, {{2, -0.5, 1.5}});
      rMC.add("MC/tracks/all/hIsPhysicalPrimary", ";is physical primary;counts", kTH1D, {{2, -0.5, 1.5}});
      rMC.add("MC/tracks/all/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
      rMC.add("MC/tracks/all/hEta", ";#it{#eta};counts", kTH1D, {etaAxis});
      rMC.add("MC/tracks/all/hPhi", ";#it{#phi} (rad);counts", kTH1D, {phiAxis});
      rMC.addClone("MC/tracks/all/", "MC/tracks/primaries/");
      rMC.addClone("MC/tracks/all/", "MC/tracks/prodByGen/");
      rMC.add("MC/tracks/hPt", ";#it{p}_{T leading} (GeV/#it{c});#it{p}_{T subleading} (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
      rMC.add("MC/tracks/hEta", ";#it{#eta}_{leading};#it{#eta}_{subleading};counts", kTH2D, {etaAxis, etaAxis});
      rMC.add("MC/tracks/hPhi", ";#it{#phi}_{leading} (rad);#it{#phi}_{subleading} (rad);counts", kTH2D, {phiAxis, phiAxis});
      // resolution
      rMC.add("MC/resolution/tracks/hPt", ";#it{p}_{T, reco} - #it{p}_{T, true} (GeV/#it{c});counts", kTH1D, {{200, -1.0, 1.0}});
      rMC.add("MC/resolution/tracks/hEta", ";#it{#eta}_{reco} - #it{#eta}_{true};counts", kTH1D, {{200, -0.2, 0.2}});
      rMC.add("MC/resolution/tracks/hPhi", ";#it{#phi}_{reco} - #it{#phi}_{true} (rad);counts", kTH1D, {{200, -0.2, 0.2}});
      // system
      rMC.add("MC/system/hM", ";#it{m} (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
      rMC.add("MC/system/hPt", ";#it{p}_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
      rMC.add("MC/system/hPt2", ";#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
      rMC.add("MC/system/hPtVsM", ";#it{m} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
      rMC.add("MC/system/hPt2VsM", ";#it{m} (GeV/#it{c}^{2});#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH2D, {mAxis, pt2Axis});
      rMC.add("MC/system/hY", ";#it{y};counts", kTH1D, {yAxis});
      rMC.add("MC/system/hPhi", ";#it{#phi} (rad);counts", kTH1D, {phiAxis});
      rMC.add("MC/system/hPhiRandom", ";#Delta#it{#phi}_{random} (rad);counts", kTH1D, {deltaPhiAxis});
      rMC.add("MC/system/hPhiCharge", ";#Delta#it{#phi}_{charge} (rad);counts", kTH1D, {deltaPhiAxis});
      rMC.add("MC/system/hPhiRandomVsM", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi} (rad);counts", kTH2D, {mAxis, deltaPhiAxis});
      rMC.add("MC/system/hPhiChargeVsM", ";#it{m} (GeV/#it{c}^{2});#Delta#it{#phi} (rad);counts", kTH2D, {mAxis, deltaPhiAxis});
      rMC.addClone("MC/system/", "MC/system/selected/");
    }

    if (context.mOptions.get<bool>("processCollisionRecoCheck"))
      rMC.add("MC/collisions/hNumOfCollisionRecos", ";number of collision reconstructions;counts", kTH1D, {{6, -0.5, 5.5}});

    if (context.mOptions.get<bool>("processResolution")) {
      // collision matching
      rResolution.add("MC/resolution/collisions/hMatch", ";matched;counts", kTH1D, {{2, -0.5, 1.5}});
      // track matching and resolutions
      rResolution.add("MC/resolution/tracks/hMatch", ";matched;counts", kTH1D, {{2, -0.5, 1.5}});
      rResolution.add("MC/resolution/tracks/hPt", ";#it{p}_{T, reco} - #it{p}_{T, true} (GeV/#it{c});counts", kTH1D, {{200, -1.0, 1.0}});
      rResolution.add("MC/resolution/tracks/hEta", ";#it{#eta}_{reco} - #it{#eta}_{true};counts", kTH1D, {{200, -0.2, 0.2}});
      rResolution.add("MC/resolution/tracks/hPhi", ";#it{#phi}_{reco} - #it{#phi}_{true} (rad);counts", kTH1D, {{200, -0.2, 0.2}});
      // dipion system resolutions (1D and 2D)
      rResolution.add("MC/resolution/system/1D/hM", ";#it{m}_{reco} - #it{m}_{true} (GeV/#it{c}^{2});counts", kTH1D, {{200, -1.0, 1.0}});
      rResolution.add("MC/resolution/system/2D/hMVsM", ";#it{m}_{true} (GeV/#it{c}^{2});#it{m}_{reco} (GeV/#it{c}^{2});counts", kTH2D, {mAxis, mAxis});
      rResolution.add("MC/resolution/system/1D/hPt", ";#it{p}_{T, reco} - #it{p}_{T, true} (GeV/#it{c});counts", kTH1D, {{200, -1.0, 1.0}});
      rResolution.add("MC/resolution/system/2D/hPtVsPt", ";#it{p}_{T, true} (GeV/#it{c});#it{p}_{T, reco} (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
      rResolution.add("MC/resolution/system/1D/hY", ";#it{y}_{reco} - #it{y}_{true};counts", kTH1D, {{200, -0.2, 0.2}});
      rResolution.add("MC/resolution/system/2D/hYVsY", ";#it{y}_{true};#it{y}_{reco};counts", kTH2D, {yAxis, yAxis});
      rResolution.add("MC/resolution/system/1D/hDeltaPhi", ";#Delta#it{#phi}_{reco} - #Delta#it{#phi}_{true} (rad);counts", kTH1D, {{2000, -1.0, 1.0}});
      rResolution.add("MC/resolution/system/2D/hDeltaPhiVsDeltaPhi", ";#Delta#it{#phi}_{true} (rad);#Delta#it{#phi}_{reco} (rad);counts", kTH2D, {deltaPhiAxis, deltaPhiAxis});
    }
  }

  static constexpr std::string_view AppliedSelections[3] = {"all/", "trackSelections/", "systemSelections/"};
  static constexpr std::string_view ChargeLabel[3] = {"unlike-sign/", "like-sign/positive/", "like-sign/negative/"};
  static constexpr std::string_view NeutronClass[5] = {"AnAn/", "0n0n/", "Xn0n/", "0nXn/", "XnXn/"};

  template <int cuts, typename C>
  void fillCollisionQcHistos(const C& collision) // fills collision QC histograms before/after cuts
  {
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hPosXY"), collision.posX(), collision.posY());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hPosZ"), collision.posZ());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hZdcCommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hZdcTime"), collision.timeZNA(), collision.timeZNC());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hZNTimeVsZNCommonEnergy"), collision.energyCommonZNA(), collision.timeZNA());
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hZNTimeVsZNCommonEnergy"), collision.energyCommonZNC(), collision.timeZNC());
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
    rQC.fill(HIST("QC/collisions/") + HIST(AppliedSelections[cuts]) + HIST("hOccupancyInTime"), collision.occupancyInTime());
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
      rSystem.fill(HIST("system/") + HIST(AppliedSelections[cuts]) + HIST(ChargeLabel[charge]) + HIST("hPt2VsM"), mass, pt * pt);
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
      rSystem.fill(HIST("system/") + HIST("selected/") + HIST(NeutronClass[neutronClass]) + HIST(ChargeLabel[charge]) + HIST("hPt2VsM"), mass, pt * pt);
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
  bool collisionPassesCuts(const C& collision, int runIndex) // collision cuts
  {
    if (!isPO) {
      if (!collision.vtxITSTPC())
        return false;
      rQC.fill(HIST("QC/collisions/hSelectionCounter"), 2);
      rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 2, runIndex);

      if (!collision.sbp())
        return false;
      rQC.fill(HIST("QC/collisions/hSelectionCounter"), 3);
      rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 3, runIndex);
    }

    if (!collision.itsROFb())
      return false;
    rQC.fill(HIST("QC/collisions/hSelectionCounter"), 4);
    rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 4, runIndex);

    if (!collision.tfb())
      return false;
    rQC.fill(HIST("QC/collisions/hSelectionCounter"), 5);
    rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 5, runIndex);

    if (std::abs(collision.posZ()) > collisionsPosZMaxCut)
      return false;
    rQC.fill(HIST("QC/collisions/hSelectionCounter"), 6);
    rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 6, runIndex);

    if (cutNumContribs) {
      if (collision.numContrib() > collisionsNumContribsMaxCut)
        return false;
      rQC.fill(HIST("QC/collisions/hSelectionCounter"), 7);
      rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 7, runIndex);
    }

    if (useRctFlag) {
      if (!isGoodRctFlag(collision)) // check RCT flags
        return false;
      rQC.fill(HIST("QC/collisions/hSelectionCounter"), 8);
      rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 8, runIndex);
    }

    if (useRecoFlag) {
      if (collision.flags() != cutRecoFlag) // check reconstruction mode
        return false;
      rQC.fill(HIST("QC/collisions/hSelectionCounter"), 9);
      rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 9, runIndex);
    }

    // if all selections passed
    return true;
  }

  template <typename T>
  bool trackPassesCuts(const T& track, int runIndex) // track cuts (PID done separately)
  {
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 0);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 0, runIndex);

    if (!track.isPVContributor())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 1);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 1, runIndex);

    if (!track.hasITS())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 2);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 2, runIndex);

    if (track.itsNCls() < tracksMinItsNClsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 3);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 3, runIndex);

    if (!cutItsLayers(track.itsClusterMap()))
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 4);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 4, runIndex);

    if (track.itsChi2NCl() > tracksMaxItsChi2NClCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 5);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 5, runIndex);

    if (!track.hasTPC())
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 6);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 6, runIndex);

    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < tracksMinTpcNClsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 7);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 7, runIndex);

    if (track.tpcChi2NCl() > tracksMaxTpcChi2NClCut || track.tpcChi2NCl() < tracksMinTpcChi2NClCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 8);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 8, runIndex);

    if (track.tpcNClsCrossedRows() < tracksMinTpcNClsCrossedRowsCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 9);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 9, runIndex);

    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < tracksMinTpcNClsCrossedOverFindableCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 10);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 10, runIndex);

    if (requireTof) {
      if (!track.hasTOF())
        return false;
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 11);
      rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 11, runIndex);
    }

    if (track.pt() < tracksMinPtCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 12);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 12, runIndex);

    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > (0.0105 + 0.0350 / std::pow(track.pt(), 1.01)))
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 13);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 13, runIndex);

    if (std::abs(eta(track.px(), track.py(), track.pz())) > pcEtaCut)
      return false;
    rQC.fill(HIST("QC/tracks/hSelectionCounter"), 14);
    rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 14, runIndex);

    // if all selections passed
    return true;
  }

  template <typename T>
  bool tracksPassPID(const T& cutTracks) // n-dimensional pion PID cut
  {
    float radiusPi = 0.0, radiusEl = 0.0, radiusKa = 0.0, radiusPr = 0.0;
    for (const auto& track : cutTracks) {
      radiusEl += std::pow(track.tpcNSigmaEl(), 2);
      radiusKa += std::pow(track.tpcNSigmaKa(), 2);
      radiusPi += std::pow(track.tpcNSigmaPi(), 2);
      radiusPr += std::pow(track.tpcNSigmaPr(), 2);
    }
    rQC.fill(HIST("QC/tracks/hPiPIDRadius"), std::sqrt(radiusPi));
    rQC.fill(HIST("QC/tracks/hElPIDRadius"), std::sqrt(radiusEl));
    rQC.fill(HIST("QC/tracks/hKaPIDRadius"), std::sqrt(radiusKa));
    rQC.fill(HIST("QC/tracks/hPrPIDRadius"), std::sqrt(radiusPr));
    if (rejectLowerProbPairs)
      return ((radiusPi < std::pow(tracksTpcNSigmaPiCut, 2)) && (radiusPi < radiusEl) && (radiusPi < radiusKa) && (radiusPi < radiusPr));
    else
      return radiusPi < std::pow(tracksTpcNSigmaPiCut, 2);
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
      charge += track.pdgCode() / std::abs(track.pdgCode());
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
    while (dPhi >= o2::constants::math::PI)
      dPhi -= o2::constants::math::TwoPI;
    while (dPhi < -o2::constants::math::PI)
      dPhi += o2::constants::math::TwoPI;
    return dPhi;
  }

  float getPhiRandom(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracksLVs) // decay phi anisotropy
  {                                                                              // two possible definitions of phi: randomize the tracks
    int indices[2] = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();            // get time-based seed
    std::shuffle(std::begin(indices), std::end(indices), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    ROOT::Math::PxPyPzMVector p1 = cutTracksLVs[indices[0]], p2 = cutTracksLVs[indices[1]];
    ROOT::Math::PxPyPzMVector pPlus = p1 + p2, pMinus = p1 - p2;
    return deltaPhi(pPlus, pMinus);
  }

  template <typename T>
  float getPhiCharge(const T& cutTracks, const std::vector<ROOT::Math::PxPyPzMVector>& cutTracksLVs)
  { // two possible definitions of phi: charge-based assignment
    ROOT::Math::PxPyPzMVector p1, p2;
    p1 = (cutTracks[0].sign() > 0) ? cutTracksLVs[0] : cutTracksLVs[1];
    p2 = (cutTracks[0].sign() > 0) ? cutTracksLVs[1] : cutTracksLVs[0];
    ROOT::Math::PxPyPzMVector pPlus = p1 + p2, pMinus = p1 - p2;
    return deltaPhi(pPlus, pMinus);
  }

  template <typename T>
  float getPhiChargeMC(const T& cutTracks, const std::vector<ROOT::Math::PxPyPzMVector>& cutTracksLVs)
  { // the same as for data but using pdg code instead of charge
    ROOT::Math::PxPyPzMVector p1, p2;
    p1 = (cutTracks[0].pdgCode() > 0) ? cutTracksLVs[0] : cutTracksLVs[1];
    p2 = (cutTracks[0].pdgCode() > 0) ? cutTracksLVs[1] : cutTracksLVs[0];
    ROOT::Math::PxPyPzMVector pPlus = p1 + p2, pMinus = p1 - p2;
    return deltaPhi(pPlus, pMinus);
  }

  // function to obtain index of run from the run number vector
  // search for passed run number in the vector and return its index +1 to use in the filling of a histogram
  int getRunIndex(int runNumber, const std::vector<int>& runNumbers)
  {
    auto it = std::find(runNumbers.begin(), runNumbers.end(), runNumber);
    if (it != runNumbers.end()) {
      return std::distance(runNumbers.begin(), it) + 1; // +1 to avoid 0 bin in histogram
    } else {
      return 0; // return 0 if run number not found
    }
  }

  template <typename C, typename T>
  void processReco(C const& collision, T const& tracks, const int runIndex)
  {
    // check if the collision run number is contained within the selectedRuns vector
    if (selectRuns && getRunIndex(collision.runNumber(), selectedRuns) == 0)
      return;

    fillCollisionQcHistos<0>(collision); // fill QC histograms before cuts
    if (!collisionPassesCuts(collision, runIndex)) // apply collision cuts
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
      fillTrackQcHistos<0>(track); // fill QC histograms before cuts

      if (!trackPassesCuts(track, runIndex)) // apply track cuts
        continue;
      cutTracks.push_back(track);
    }
    rQC.fill(HIST("QC/tracks/trackSelections/hRemainingTracks"), cutTracks.size());

    if (static_cast<int>(cutTracks.size()) != 2) // further consider only two pion systems
      return;
    for (int i = 0; i < 2; i++) {
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 15);
      rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 15, runIndex);
    }
    rQC.fill(HIST("QC/tracks/trackSelections/hTpcNSigmaPi2D"), cutTracks[0].tpcNSigmaPi(), cutTracks[1].tpcNSigmaPi());
    rQC.fill(HIST("QC/tracks/trackSelections/hTpcNSigmaEl2D"), cutTracks[0].tpcNSigmaEl(), cutTracks[1].tpcNSigmaEl());
    rQC.fill(HIST("QC/tracks/trackSelections/hTpcNSigmaKa2D"), cutTracks[0].tpcNSigmaKa(), cutTracks[1].tpcNSigmaKa());
    rQC.fill(HIST("QC/tracks/trackSelections/hTpcNSigmaPr2D"), cutTracks[0].tpcNSigmaPr(), cutTracks[1].tpcNSigmaPr());

    // create a vector of 4-vectors for selected tracks
    std::vector<ROOT::Math::PxPyPzMVector> cutTracksLVs;
    for (const auto& cutTrack : cutTracks) {
      cutTracksLVs.push_back(ROOT::Math::PxPyPzMVector(cutTrack.px(), cutTrack.py(), cutTrack.pz(), o2::constants::physics::MassPionCharged)); // apriori assume pion mass
    }

    // differentiate leading- and subleading-momentum tracks
    auto leadingTrack = momentum(cutTracks[0].px(), cutTracks[0].py(), cutTracks[0].pz()) > momentum(cutTracks[1].px(), cutTracks[1].py(), cutTracks[1].pz()) ? cutTracks[0] : cutTracks[1];
    auto subleadingTrack = (leadingTrack == cutTracks[0]) ? cutTracks[1] : cutTracks[0];

    float leadingPt = leadingTrack.pt();
    float subleadingPt = subleadingTrack.pt();
    float leadingEta = eta(leadingTrack.px(), leadingTrack.py(), leadingTrack.pz());
    float subleadingEta = eta(subleadingTrack.px(), subleadingTrack.py(), subleadingTrack.pz());
    float leadingPhi = phi(leadingTrack.px(), leadingTrack.py());
    float subleadingPhi = phi(subleadingTrack.px(), subleadingTrack.py());
    float phiRandom = getPhiRandom(cutTracksLVs);
    float phiCharge = getPhiCharge(cutTracks, cutTracksLVs);
    
    // fill recoTree
    recoTree(collision.flags(), collision.runNumber(), collision.posZ(), collision.occupancyInTime(), collision.hadronicRate(),
             collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
             collision.timeFT0A(), collision.timeFT0C(), collision.timeFV0A(), collision.timeFDDA(), collision.timeFDDC(),
             energyCommonZNA, energyCommonZNC, timeZNA, timeZNC, neutronClass,
             leadingTrack.sign(), subleadingTrack.sign(),
             leadingPt, subleadingPt,
             leadingEta, subleadingEta,
             leadingPhi, subleadingPhi,
             leadingTrack.tpcNSigmaPi(), subleadingTrack.tpcNSigmaPi(),
             leadingTrack.tpcNSigmaEl(), subleadingTrack.tpcNSigmaEl(),
             leadingTrack.tpcNSigmaKa(), subleadingTrack.tpcNSigmaKa(),
             leadingTrack.tpcNSigmaPr(), subleadingTrack.tpcNSigmaPr());
    
    if (!tracksPassPID(cutTracks)) // apply PID cut
      return;
    
    for (const auto& cutTrack : cutTracks) {
      rQC.fill(HIST("QC/tracks/hSelectionCounter"), 16);
      rQC.fill(HIST("QC/tracks/hSelectionCounterPerRun"), 16, runIndex);
      fillTrackQcHistos<1>(cutTrack); // fill QC histograms after cuts
    }
    rQC.fill(HIST("QC/tracks/hTofHitCheck"), leadingTrack.hasTOF(), subleadingTrack.hasTOF());
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
        rSystem.fill(HIST("system/all/unlike-sign/hRecoSettingVsM"), mass, collision.flags());
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
  void processMC(C const& mcCollision, T const& mcParticles, const int runNumber)
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
      if (mcParticle.producedByGenerator()) {
        rMC.fill(HIST("MC/tracks/prodByGen/hPdgCode"), mcParticle.pdgCode());
        rMC.fill(HIST("MC/tracks/prodByGen/hProducedByGenerator"), mcParticle.producedByGenerator());
        rMC.fill(HIST("MC/tracks/prodByGen/hIsPhysicalPrimary"), mcParticle.isPhysicalPrimary());
        rMC.fill(HIST("MC/tracks/prodByGen/hPt"), pt(mcParticle.px(), mcParticle.py()));
        rMC.fill(HIST("MC/tracks/prodByGen/hEta"), eta(mcParticle.px(), mcParticle.py(), mcParticle.pz()));
        rMC.fill(HIST("MC/tracks/prodByGen/hPhi"), phi(mcParticle.px(), mcParticle.py()));
      }
      if (mcParticle.isPhysicalPrimary()) {
        rMC.fill(HIST("MC/tracks/primaries/hPdgCode"), mcParticle.pdgCode());
        rMC.fill(HIST("MC/tracks/primaries/hProducedByGenerator"), mcParticle.producedByGenerator());
        rMC.fill(HIST("MC/tracks/primaries/hIsPhysicalPrimary"), mcParticle.isPhysicalPrimary());
        rMC.fill(HIST("MC/tracks/primaries/hPt"), pt(mcParticle.px(), mcParticle.py()));
        rMC.fill(HIST("MC/tracks/primaries/hEta"), eta(mcParticle.px(), mcParticle.py(), mcParticle.pz()));
        rMC.fill(HIST("MC/tracks/primaries/hPhi"), phi(mcParticle.px(), mcParticle.py()));
      }
      if (mcParticle.has_daughters()) {
        rMC.fill(HIST("MC/tracks/all/hMotherPdgCode"), mcParticle.pdgCode());
        if (mcParticle.pdgCode() != kRho770_0)
          continue; // consider only rho0s
        for (const auto& daughter : mcParticle.template daughters_as<T>()) {
          if (!daughter.isPhysicalPrimary() || std::abs(daughter.pdgCode()) != kPiPlus)
            continue;
          cutMcParticles.push_back(daughter);
          ROOT::Math::PxPyPzMVector pionLV;
          pionLV.SetPxPyPzE(daughter.px(), daughter.py(), daughter.pz(), daughter.e());
          mcParticlesLVs.push_back(pionLV);
        }
      }
    }
    rMC.fill(HIST("MC/collisions/hNPions"), cutMcParticles.size());

    if (static_cast<int>(cutMcParticles.size()) != 2)
      return;
    if (mcParticlesLVs.size() != cutMcParticles.size()) // sanity check
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

    auto leadingPion = momentum(cutMcParticles[0].px(), cutMcParticles[0].py(), cutMcParticles[0].pz()) > momentum(cutMcParticles[1].px(), cutMcParticles[1].py(), cutMcParticles[1].pz()) ? cutMcParticles[0] : cutMcParticles[1];
    auto subleadingPion = (leadingPion == cutMcParticles[0]) ? cutMcParticles[1] : cutMcParticles[0];
    rMC.fill(HIST("MC/tracks/hPt"), pt(leadingPion.px(), leadingPion.py()), pt(subleadingPion.px(), subleadingPion.py()));
    rMC.fill(HIST("MC/tracks/hEta"), eta(leadingPion.px(), leadingPion.py(), leadingPion.pz()), eta(subleadingPion.px(), subleadingPion.py(), subleadingPion.pz()));
    rMC.fill(HIST("MC/tracks/hPhi"), phi(leadingPion.px(), leadingPion.py()), phi(subleadingPion.px(), subleadingPion.py()));

    rMC.fill(HIST("MC/system/hM"), mass);
    rMC.fill(HIST("MC/system/hPt"), pT);
    rMC.fill(HIST("MC/system/hPtVsM"), mass, pT);
    rMC.fill(HIST("MC/system/hPt2VsM"), mass, pT * pT);
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
    mcTree(runNumber, mcCollision.posZ(),
           leadingPion.pdgCode() / std::abs(leadingPion.pdgCode()), subleadingPion.pdgCode() / std::abs(subleadingPion.pdgCode()),
           pt(leadingPion.px(), leadingPion.py()), pt(subleadingPion.px(), subleadingPion.py()),
           eta(leadingPion.px(), leadingPion.py(), leadingPion.pz()), eta(subleadingPion.px(), subleadingPion.py(), subleadingPion.pz()),
           phi(leadingPion.px(), leadingPion.py()), phi(subleadingPion.px(), subleadingPion.py()));
  }

  template <typename C>
  void checkNumberOfCollisionReconstructions(C const& collisions)
  {
    rMC.fill(HIST("MC/collisions/hNumOfCollisionRecos"), collisions.size());
  }

  void processSGdata(FullUdSgCollision const& collision, FullUdTracks const& tracks)
  {
    int runIndex = getRunIndex(collision.runNumber(), runNumbers);
    rQC.fill(HIST("QC/collisions/hSelectionCounter"), 0); // all collisions
    rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 0, runIndex);

    if (cutGapSide && collision.gapSide() != gapSide)
      return;
    if (useTrueGap && (collision.gapSide() != sgSelector.trueGap(collision, cutTrueGapSideFV0, cutTrueGapSideFT0A, cutTrueGapSideFT0C, cutTrueGapSideZDC))) // check true gap side
      return;
    rQC.fill(HIST("QC/collisions/hSelectionCounter"), 1); // only double-gap collisions
    rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 1, runIndex);

    processReco(collision, tracks, runIndex);
  }
  PROCESS_SWITCH(UpcRhoAnalysis, processSGdata, "analyse SG data", true);

  void processDGdata(FullUdDgCollision const& collision, FullUdTracks const& tracks)
  {
    int runIndex = getRunIndex(collision.runNumber(), runNumbers);
    rQC.fill(HIST("QC/collisions/hSelectionCounter"), 1); // no single-gap collisions in dataset
    rQC.fill(HIST("QC/collisions/hSelectionCounterPerRun"), 1, runIndex);

    processReco(collision, tracks, runIndex);
  }
  PROCESS_SWITCH(UpcRhoAnalysis, processDGdata, "analyse DG data", false);

  void processMCdata(aod::UDMcCollision const& mcCollision, aod::UDMcParticles const& mcParticles)
  {
    processMC(mcCollision, mcParticles, -1);
  }
  PROCESS_SWITCH(UpcRhoAnalysis, processMCdata, "analyse MC data", false);

  void processMCdataWithBCs(aod::UDMcCollision const& mcCollision, aod::UDMcParticles const& mcParticles, aod::BCs const& bcs)
  {
    int runNumber = -1;
    if (bcs.size() != 0) {
      auto bc = bcs.begin();
      runNumber = bc.runNumber();
    }
    processMC(mcCollision, mcParticles, runNumber);
  }
  PROCESS_SWITCH(UpcRhoAnalysis, processMCdataWithBCs, "analyse MC data with BCs (only with on-the-fly skimming)", false);

  void processResolution(soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDMcCollsLabels>::iterator const& collision, soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags, aod::UDMcTrackLabels> const& tracks, aod::UDMcCollisions const&, aod::UDMcParticles const&)
  {
    rResolution.fill(HIST("MC/resolution/collisions/hMatch"), 0);
    if (!collision.has_udMcCollision())
      return;
    rResolution.fill(HIST("MC/resolution/collisions/hMatch"), 1);

    std::vector<decltype(tracks.begin().udMcParticle())> trueTracks;
    std::vector<decltype(tracks.begin())> recoTracks;
    std::vector<ROOT::Math::PxPyPzMVector> truePionLVs, recoPionLVs;

    for (const auto& track : tracks) {
      rResolution.fill(HIST("MC/resolution/tracks/hMatch"), 0);
      if (!track.has_udMcParticle())
        continue;
      rResolution.fill(HIST("MC/resolution/tracks/hMatch"), 1);
      auto mcParticle = track.udMcParticle();
      rResolution.fill(HIST("MC/resolution/tracks/hPt"), pt(track.px(), track.py()) - pt(mcParticle.px(), mcParticle.py()));
      rResolution.fill(HIST("MC/resolution/tracks/hEta"), eta(track.px(), track.py(), track.pz()) - eta(mcParticle.px(), mcParticle.py(), mcParticle.pz()));
      rResolution.fill(HIST("MC/resolution/tracks/hPhi"), phi(track.px(), track.py()) - phi(mcParticle.px(), mcParticle.py()));
      if (std::abs(mcParticle.pdgCode()) != kPiPlus && !mcParticle.isPhysicalPrimary())
        continue;
      truePionLVs.push_back(ROOT::Math::PxPyPzMVector(mcParticle.px(), mcParticle.py(), mcParticle.pz(), o2::constants::physics::MassPionCharged));
      trueTracks.push_back(mcParticle);
      recoPionLVs.push_back(ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged));
      recoTracks.push_back(track);
    }

    if (truePionLVs.size() != 2 || recoPionLVs.size() != 2)
      return;
    
    ROOT::Math::PxPyPzMVector trueSystem = reconstructSystem(truePionLVs);
    const float trueDeltaPhi = getPhiChargeMC(trueTracks, truePionLVs);
    ROOT::Math::PxPyPzMVector recoSystem = reconstructSystem(recoPionLVs);
    const float recoDeltaPhi = getPhiCharge(recoTracks, recoPionLVs);
    
    rResolution.fill(HIST("MC/resolution/system/1D/hM"), recoSystem.M() - trueSystem.M());
    rResolution.fill(HIST("MC/resolution/system/2D/hMVsM"), trueSystem.M(), recoSystem.M());
    rResolution.fill(HIST("MC/resolution/system/1D/hPt"), recoSystem.Pt() - trueSystem.Pt());
    rResolution.fill(HIST("MC/resolution/system/2D/hPtVsPt"), trueSystem.Pt(), recoSystem.Pt());
    rResolution.fill(HIST("MC/resolution/system/1D/hY"), recoSystem.Rapidity() - trueSystem.Rapidity());
    rResolution.fill(HIST("MC/resolution/system/2D/hYVsY"), trueSystem.Rapidity(), recoSystem.Rapidity());
    rResolution.fill(HIST("MC/resolution/system/1D/hDeltaPhi"), recoDeltaPhi - trueDeltaPhi);
    rResolution.fill(HIST("MC/resolution/system/2D/hDeltaPhiVsDeltaPhi"), trueDeltaPhi, recoDeltaPhi);
  }
  PROCESS_SWITCH(UpcRhoAnalysis, processResolution, "check resolution of kinematic variables", false);

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
