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
/// \brief
/// \author Sara Haidlova, sara.haidlova@cern.ch
/// \since March 2024

#include <vector>
#include <utility>

// O2 headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"

// O2Physics headers
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/UPCJpsiCentralBarrelCorrHelper.h"
#include "PWGUD/Core/SGSelector.h"

// ROOT headers
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;

SGSelector sgSelector;

namespace o2::aod
{
namespace tree
{
// misc event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
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

// J/Psi
DECLARE_SOA_COLUMN(JPsiPt, jpsiPt, double);
DECLARE_SOA_COLUMN(JPsiEta, jpsiEta, double);
DECLARE_SOA_COLUMN(JPsiPhi, jpsiPhi, double);
DECLARE_SOA_COLUMN(JPsiM, jpsiM, double);
DECLARE_SOA_COLUMN(JPsiRap, jpsiRap, double);

// asymmetry and correlation
DECLARE_SOA_COLUMN(JPsiPhiRandom, jpsiPhiRandom, double);
DECLARE_SOA_COLUMN(JPsiPhiCharge, jpsiPhiCharge, double);
DECLARE_SOA_COLUMN(JPsiPhiCS, jpsiPhiCS, double);
DECLARE_SOA_COLUMN(JPsiCosThetaCS, jpsiCosThetaCS, double);

// muon tracks
DECLARE_SOA_COLUMN(TrackSign1, trackSign1, double);
DECLARE_SOA_COLUMN(TrackPt1, trackPt1, double);
DECLARE_SOA_COLUMN(TrackEta1, trackEta1, double);
DECLARE_SOA_COLUMN(TrackPhi1, trackPhi1, double);
DECLARE_SOA_COLUMN(TrackSign2, trackSign2, double);
DECLARE_SOA_COLUMN(TrackPt2, trackPt2, double);
DECLARE_SOA_COLUMN(TrackEta2, trackEta2, double);
DECLARE_SOA_COLUMN(TrackPhi2, trackPhi2, double);
} // namespace tree
DECLARE_SOA_TABLE(Tree, "AOD", "TREE",
                  tree::RunNumber, tree::GlobalBC,
                  tree::PosX, tree::PosY, tree::PosZ, tree::TotalFT0AmplitudeA, tree::TotalFT0AmplitudeC, tree::TotalFV0AmplitudeA, tree::TotalFDDAmplitudeA, tree::TotalFDDAmplitudeC,
                  tree::TimeFT0A, tree::TimeFT0C, tree::TimeFV0A, tree::TimeFDDA, tree::TimeFDDC,
                  tree::EnergyCommonZNA, tree::EnergyCommonZNC, tree::TimeZNA, tree::TimeZNC, tree::NeutronClass,
                  tree::JPsiPt, tree::JPsiEta, tree::JPsiPhi, tree::JPsiRap, tree::JPsiM, tree::JPsiPhiRandom, tree::JPsiPhiCharge, tree::JPsiPhiCS, tree::JPsiCosThetaCS,
                  tree::TrackSign1, tree::TrackPt1, tree::TrackEta1, tree::TrackPhi1,
                  tree::TrackSign2, tree::TrackPt2, tree::TrackEta2, tree::TrackPhi2);
namespace tree_mc
{
// misc event info
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
// event vertex
DECLARE_SOA_COLUMN(PosX, posX, double);
DECLARE_SOA_COLUMN(PosY, posY, double);
DECLARE_SOA_COLUMN(PosZ, posZ, double);

// J/Psi
DECLARE_SOA_COLUMN(JPsiPt, jpsiPt, double);
DECLARE_SOA_COLUMN(JPsiEta, jpsiEta, double);
DECLARE_SOA_COLUMN(JPsiPhi, jpsiPhi, double);
DECLARE_SOA_COLUMN(JPsiM, jpsiM, double);
DECLARE_SOA_COLUMN(JPsiRap, jpsiRap, double);

// asymmetry and correlation
DECLARE_SOA_COLUMN(JPsiPhiRandom, jpsiPhiRandom, double);
DECLARE_SOA_COLUMN(JPsiPhiCharge, jpsiPhiCharge, double);
DECLARE_SOA_COLUMN(JPsiPhiCS, jpsiPhiCS, double);
DECLARE_SOA_COLUMN(JPsiCosThetaCS, jpsiCosThetaCS, double);

// muon tracks
DECLARE_SOA_COLUMN(TrackSign1, trackSign1, double);
DECLARE_SOA_COLUMN(TrackPt1, trackPt1, double);
DECLARE_SOA_COLUMN(TrackEta1, trackEta1, double);
DECLARE_SOA_COLUMN(TrackPhi1, trackPhi1, double);
DECLARE_SOA_COLUMN(TrackSign2, trackSign2, double);
DECLARE_SOA_COLUMN(TrackPt2, trackPt2, double);
DECLARE_SOA_COLUMN(TrackEta2, trackEta2, double);
DECLARE_SOA_COLUMN(TrackPhi2, trackPhi2, double);
} // namespace tree_mc
DECLARE_SOA_TABLE(TreeMC, "AOD", "TREEMC",
                  tree_mc::GlobalBC,
                  tree_mc::JPsiPt, tree_mc::JPsiEta, tree_mc::JPsiPhi, tree_mc::JPsiRap, tree_mc::JPsiM, tree_mc::JPsiPhiRandom, tree_mc::JPsiPhiCharge, tree_mc::JPsiPhiCS, tree_mc::JPsiCosThetaCS,
                  tree_mc::TrackSign1, tree_mc::TrackPt1, tree_mc::TrackEta1, tree_mc::TrackPhi1,
                  tree_mc::TrackSign2, tree_mc::TrackPt2, tree_mc::TrackEta2, tree_mc::TrackPhi2);
//} // namespace tree_mc
} // namespace o2::aod

struct UpcJpsiCentralBarrel {
  Produces<o2::aod::Tree> Tree;
  Produces<o2::aod::TreeMC> TreeMC;

  // configurable axes
  ConfigurableAxis axisIVM{"axisIVM", {500.0f, 2.0f, 4.5f}, "M_#it{inv} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisIVMWide{"axisIVMWide", {350.0f, 0.0f, 4.5f}, "M_#it{inv} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisPt{"axisPt", {250.0f, 0.1f, 3.0f}, "#it{p}_T (GeV/#it{c})"};
  ConfigurableAxis axisPt2{"axisPt2", {250.0f, 0.0f, 0.01f}, "#it{p}_T (GeV/#it{c})"};
  ConfigurableAxis axisPt2wide{"axisPt2wide", {250.0f, 0.0f, 0.04f}, "#it{p}_T (GeV/#it{c})"};
  ConfigurableAxis axisP{"axisP", {250.0f, 0.1f, 3.0f}, "#it{p} (GeV/#it{c})"};
  ConfigurableAxis axisEta{"axisEta", {250.0f, -1.5f, 1.5f}, "#eta (-)"};
  ConfigurableAxis axisCounter{"axisCounter", {20.0f, 0.0f, 20.0f}, "Number of events (-)"};
  ConfigurableAxis axisPhi{"axisPhi", {250.0f, 0, TwoPI}, "#phi (rad)"};
  ConfigurableAxis axisAccAngle{"axisAccAngle", {250.0f, -0.2f, 0.2f}, "accAngle"};
  ConfigurableAxis axisAngTheta{"axisAngTheta", {250.0f, -1.5f, 1.5f}, "cos #theta (-)"};
  ConfigurableAxis axisTPC{"axisTPC", {1000.0f, 0, 200.0f}, "TPC d#it{E}/d#it{x}"};
  ConfigurableAxis axisBetaTOF{"axisBetaTOF", {100.0f, 0, 1.5}, "TOF #beta"};
  ConfigurableAxis axisSigma{"axisSigma", {50, -25, 25}, "#sigma"};
  ConfigurableAxis axisZDCEnergy{"axisZDCEnergy", {250, -5.0, 20.0}, "ZDC energy"};
  ConfigurableAxis axisZDCTime{"axisZDCTime", {200, -10.0, 10.0}, "ZDC time"};
  ConfigurableAxis axisDCA{"axisDCA", {1000, -20.0, 20.0}, "DCA"};
  ConfigurableAxis axisChi2{"axisChi2", {200, -10.0, 40}, "Chi2"};
  ConfigurableAxis axisIVMSel{"axisIVMSel", {1000, 0.0, 10.0}, "IVM"};
  ConfigurableAxis axisCounterSel{"axisCounterSel", {1000, 0.0, 200.0}, "Selection"};
  ConfigurableAxis axisITSNCls{"axisITSNCls", {10, 0.0, 10.0}, "ITSNCls"};
  ConfigurableAxis axisTPCNCls{"axisTPCNCls", {170, -1, 160.0}, "TPCNCls"};
  ConfigurableAxis axisTPCCrossed{"axisTPCCrossed", {300, 20, 170.0}, "TPCCrossedRows"};

  // configurable cuts (modify in json)
  // track quality cuts
  Configurable<int> TPCNClsCrossedRows{"TPCNClsCrossedRows", 70, "number of crossed rows in TPC"};
  Configurable<bool> TOFAtLeastOneProton{"TOFAtLeastOneProton", false, "at least one candidate track has TOF hits"};
  Configurable<bool> TOFBothProtons{"TOFBothProtons", false, "both candidate protons have TOF hits"};
  Configurable<bool> TOFOneProton{"TOFOneProton", false, "one candidate proton has TOF hits"};
  Configurable<bool> DCAcut{"DCAcut", false, "DCA cut from run2."};
  Configurable<bool> newCutTPC{"newCutTPC", false, "New cuts for TPC quality tracks."};
  Configurable<float> EtaCut{"EtaCut", 0.9f, "acceptance cut per track"};
  Configurable<float> cutPtTrack{"cutPtTrack", 0.1f, "pT cut per track"};
  Configurable<float> cutVertexZ{"cutVertexZ", 10.0f, "cut on vertex position in Z"};
  Configurable<float> RapCut{"RapCut", 0.9f, "choose event in midrapidity"};
  Configurable<float> dcaZCut{"dcaZCut", 2, "cut on the impact parameter in z of the track to the PV"};
  Configurable<float> dcaXYCut{"dcaXYCut", 1e10, "cut on the impact parameter in xy of the track to the PV"};
  Configurable<int> ITSNClsCut{"ITSNClsCut", 4, "minimal number of ITS clusters"};
  Configurable<int> ITSChi2NClsCut{"ITSChi2NClsCut", 36, "minimal Chi2/cluster for the ITS track"};
  Configurable<int> TPCNClsCrossedRowsCut{"TPCNClsCrossedRowsCut", 70, "minimal number of crossed TPC rows"};
  Configurable<int> TPCChi2NCls{"TPCChi2NCls", 4, "minimal Chi2/cluster for the TPC track"};
  Configurable<float> TPCMinNCls{"TPCMinNCls", 3, "minimum number of TPC clusters"};
  Configurable<float> TPCCrossedOverFindable{"TPCCrossedOverFindable", 3, "number of TPC crossed rows over findable clusters"};

  // ZDC classes cuts
  Configurable<double> ZNenergyCut{"ZNenergyCut", 0.0, "ZN common energy cut"};
  Configurable<double> ZNtimeCut{"ZNtimeCut", 2.0, "ZN time cut"};

  // Analysis cuts
  Configurable<float> maxJpsiMass{"maxJpsiMass", 3.18, "Maximum of the jpsi peak for peak cut"};
  Configurable<float> minJpsiMass{"minJpsiMass", 3.0, "Minimum of the jpsi peak for peak cut"};

  // SG cuts
  Configurable<int> whichGapSide{"whichGapSide", 2, {"0 for side A, 1 for side C, 2 for both sides"}};
  Configurable<bool> useTrueGap{"useTrueGap", true, {"Calculate gapSide for a given FV0/FT0/ZDC thresholds"}};
  Configurable<float> cutMyGapSideFV0{"FV0", 100, "FV0A threshold for SG selector"};
  Configurable<float> cutMyGapSideFT0A{"FT0A", 200., "FT0A threshold for SG selector"};
  Configurable<float> cutMyGapSideFT0C{"FT0C", 100., "FT0C threshold for SG selector"};
  Configurable<float> cutMyGapSideZDC{"ZDC", 1000., "ZDC threshold for SG selector"};

  // process cuts
  Configurable<bool> doMuons{"doMuons", true, "Provide muon plots."};
  Configurable<bool> doElectrons{"doElectrons", true, "Provide electron plots."};
  Configurable<bool> doProtons{"doProtons", true, "Provide proton plots."};
  Configurable<bool> chargeOrdered{"chargeOrdered", false, "Order tracks based on charge."};
  Configurable<bool> doOnlyUPC{"doOnlyUPC", false, "Analyse only collisions with UPC settings."};
  Configurable<bool> doOnlyStandard{"doOnlyStandard", false, "Analyse only collisions with standard settings."};

  // initialize histogram registry
  HistogramRegistry Statistics{
    "Statistics",
    {}};

  HistogramRegistry RawData{
    "RawData",
    {}};

  HistogramRegistry Selections{
    "Selections",
    {}};

  HistogramRegistry PVContributors{
    "PVContributors",
    {}};

  HistogramRegistry TG{
    "TG",
    {}};

  HistogramRegistry TGdaug{
    "TGdaug",
    {}};

  HistogramRegistry TGdaugCand{
    "TGdaugCand",
    {}};

  HistogramRegistry JPsiToDaug{
    "JPsiToDaug",
    {}};

  HistogramRegistry Correlation{
    "Correlation",
    {}};

  HistogramRegistry Asymmetry{
    "Asymmetry",
    {}};

  HistogramRegistry MC{
    "MC",
    {}};

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDTrackFull = UDTracksFull::iterator;
  using SGUDCollisionFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>::iterator;
  using MCUDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags, aod::UDMcTrackLabels>;
  using MCUDCollisionFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDMcCollsLabels>::iterator;
  using MCSGUDCollisionFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDMcCollsLabels>::iterator;

  void init(InitContext&)
  {

    // statistics histograms for counters
    Statistics.add("Statistics/hNumberOfCollisions", "hNumberOfCollisions", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberOfTracks", "hNumberOfTracks", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGT", "hNumberGT", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hCutCounterCollisions", "hCutCounterCollisions", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hCutCounterTracks", "hCutCounterTracks", {HistType::kTH1F, {axisCounter}});

    // raw data histograms
    RawData.add("RawData/hTrackPt", "hTrackPt", {HistType::kTH1F, {axisPt}});
    RawData.add("RawData/hTrackEta", "hTrackEta", {HistType::kTH1F, {axisEta}});
    RawData.add("RawData/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});
    RawData.add("RawData/hTrackDCAXYZ", "hTrackDCAXYZ", {HistType::kTH2F, {axisDCA, axisDCA}});
    RawData.add("RawData/hTPCNClsFindable", "hTPCNClsFindable", {HistType::kTH1F, {axisTPC}});
    RawData.add("RawData/hTPCNClsFindableMinusFound", "hTPCNClsFindableMinusFound", {HistType::kTH1F, {axisTPC}});
    RawData.add("RawData/hITSNCls", "hITSNCls", {HistType::kTH1F, {axisCounter}});
    RawData.add("RawData/hTPCNCls", "hTPCNCls", {HistType::kTH1F, {axisCounter}});
    RawData.add("RawData/hITSChi2NCls", "hITSChi2NCls", {HistType::kTH1F, {axisChi2}});
    RawData.add("RawData/hTPCChi2NCls", "hTPCChi2NCls", {HistType::kTH1F, {axisChi2}});
    RawData.add("RawData/hPositionZ", "hPositionZ", {HistType::kTH1F, {axisSigma}});
    RawData.add("RawData/hPositionX", "hPositionX", {HistType::kTH1F, {axisSigma}});
    RawData.add("RawData/hPositionY", "hPositionY", {HistType::kTH1F, {axisSigma}});
    RawData.add("RawData/hPositionXY", "hPositionXY", {HistType::kTH2F, {axisSigma, axisSigma}});
    RawData.add("RawData/QA/ZDC/hZNACommonEnergy", "hZNACommonEnergy", {HistType::kTH1F, {axisZDCEnergy}});
    RawData.add("RawData/QA/ZDC/hZNCCommonEnergy", "hZNCCommonEnergy", {HistType::kTH1F, {axisZDCEnergy}});
    RawData.add("RawData/QA/ZDC/hZNAvsZNCCommonEnergy", "hZNAvsZNCCommonEnergy", {HistType::kTH2F, {axisZDCEnergy, axisZDCEnergy}});
    RawData.add("RawData/QA/ZDC/hZNATime", "hZNATime", {HistType::kTH1F, {axisZDCTime}});
    RawData.add("RawData/QA/ZDC/hZNCTime", "hZNCTime", {HistType::kTH1F, {axisZDCTime}});
    RawData.add("RawData/QA/ZDC/hZNAvsZNCTime", "hZNAvsZNCTime", {HistType::kTH2F, {axisZDCTime, axisZDCTime}});
    RawData.add("RawData/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    RawData.add("RawData/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    RawData.add("RawData/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    RawData.add("RawData/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    RawData.add("RawData/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    RawData.add("RawData/QA/FIT/hTotFT0AmplitudeA", "hTotFT0AmplitudeA", {HistType::kTH1D, {{1000, 0.0, 1000}}});
    RawData.add("RawData/QA/FIT/hTotFT0AmplitudeC", "hTotFT0AmplitudeC", {HistType::kTH1D, {{1000, 0.0, 1000}}});
    RawData.add("RawData/QA/FIT/hTotFV0AmplitudeA", "hTotFV0AmplitudeA", {HistType::kTH1D, {{1000, 0.0, 1000}}});
    RawData.add("RawData/QA/FIT/hTotFDDAmplitudeA", "hTotFDDAmplitudeA", {HistType::kTH1D, {{1000, 0.0, 1000}}});
    RawData.add("RawData/QA/FIT/hTotFDDAmplitudeC", "hTotFDDAmplitudeC", {HistType::kTH1D, {{1000, 0.0, 1000}}});
    RawData.add("RawData/QA/FIT/hTimeFT0A", "hTimeFT0A", {HistType::kTH1D, {{200, -100, 100}}});
    RawData.add("RawData/QA/FIT/hTimeFT0C", "hTimeFT0C", {HistType::kTH1D, {{200, -100, 100}}});
    RawData.add("RawData/QA/FIT/hTimeFV0A", "hTimeFV0A", {HistType::kTH1D, {{200, -100, 100}}});
    RawData.add("RawData/QA/FIT/hTimeFDDA", "hTimeFDDA", {HistType::kTH1D, {{200, -100, 100}}});
    RawData.add("RawData/QA/FIT/hTimeFDDC", "hTimeFDDC", {HistType::kTH1D, {{200, -100, 100}}});

    // Selection checks
    Selections.add("Selections/Electron/Mass/Leading/hITSNCls", "hITSNCls", {HistType::kTH2F, {axisIVMSel, axisITSNCls}});
    Selections.add("Selections/Electron/Mass/Leading/hITSChi2NCls", "hITSChi2NCls", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCNCls", "hTPCNCls", {HistType::kTH2F, {axisIVMSel, axisTPCNCls}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCChi2NCls", "hTPCChi2NCls", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCNClsCrossedRows", "hTPCNClsCrossedRows", {HistType::kTH2F, {axisIVMSel, axisTPCCrossed}});
    Selections.addClone("Selections/Electron/Mass/Leading/", "Selections/Electron/Mass/Subleading/");
    Selections.addClone("Selections/Electron/Mass/", "Selections/Electron/Rapidity/");
    Selections.addClone("Selections/Electron/Mass/", "Selections/Electron/Pt/");
    Selections.addClone("Selections/Electron/", "Selections/Muon/");

    // PVContributors histograms
    PVContributors.add("PVContributors/hTrackPt", "hTrackPt", {HistType::kTH1F, {axisPt}});
    PVContributors.add("PVContributors/hTrackEta", "hTrackEta", {HistType::kTH1F, {axisEta}});
    PVContributors.add("PVContributors/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});
    PVContributors.add("PVContributors/hTPCNClsFindable", "hTPCNClsFindable", {HistType::kTH1F, {axisTPC}});
    PVContributors.add("PVContributors/hTPCNClsFindableMinusFound", "hTPCNClsFindableMinusFound", {HistType::kTH1F, {axisTPC}});
    PVContributors.add("PVContributors/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    PVContributors.add("PVContributors/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    PVContributors.add("PVContributors/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    PVContributors.add("PVContributors/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    PVContributors.add("PVContributors/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // TG histograms
    TG.add("TG/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TG.add("TG/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TG.add("TG/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TG.add("TG/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TG.add("TG/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TG.add("TG/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TG.add("TG/hTPCNClsFindable", "hTPCNClsFindable", {HistType::kTH1F, {axisTPC}});
    TG.add("TG/hTPCNClsFindableMinusFound", "hTPCNClsFindableMinusFound", {HistType::kTH1F, {axisTPC}});
    TG.add("TG/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TG.add("TG/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TG.add("TG/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TG.add("TG/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TG.add("TG/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TG.add("TG/TPC/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});

    // TGdaug histograms
    TGdaug.add("TGdaug/Muon/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGdaug.add("TGdaug/Muon/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGdaug.add("TGdaug/Muon/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGdaug.add("TGdaug/Muon/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGdaug.add("TGdaug/Muon/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGdaug.add("TGdaug/Muon/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGdaug.add("TGdaug/Muon/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGdaug.add("TGdaug/Muon/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGdaug.add("TGdaug/Muon/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGdaug.add("TGdaug/Muon/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGdaug.add("TGdaug/Muon/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGdaug.add("TGdaug/Muon/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TGdaug.addClone("TGdaug/Muon/", "TGdaug/Electron/");
    TGdaug.addClone("TGdaug/Muon/", "TGdaug/Proton/");

    // TGmuCand histograms
    TGdaugCand.add("TGdaugCand/Muon/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGdaugCand.add("TGdaugCand/Muon/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGdaugCand.add("TGdaugCand/Muon/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGdaugCand.add("TGdaugCand/Muon/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGdaugCand.add("TGdaugCand/Muon/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGdaugCand.add("TGdaugCand/Muon/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGdaugCand.add("TGdaugCand/Muon/hPairPt", "hPairPt", {HistType::kTH1F, {axisPt}});
    TGdaugCand.add("TGdaugCand/Muon/hPairIVM", "hPairIVM", {HistType::kTH1F, {axisIVMWide}});
    TGdaugCand.add("TGdaugCand/Muon/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axisPt}});
    TGdaugCand.add("TGdaugCand/Muon/hJpsiRap", "hJpsiRap", {HistType::kTH1F, {axisEta}});
    TGdaugCand.add("TGdaugCand/Muon/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGdaugCand.add("TGdaugCand/Muon/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGdaugCand.add("TGdaugCand/Muon/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGdaugCand.add("TGdaugCand/Muon/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGdaugCand.add("TGdaugCand/Muon/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGdaugCand.add("TGdaugCand/Muon/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TGdaugCand.addClone("TGdaugCand/Muon/", "TGdaugCand/Electron/");
    TGdaugCand.addClone("TGdaugCand/Muon/", "TGdaugCand/Proton/");

    // JPsiToEl histograms
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/XnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToDaug.add("JPsiToDaug/Electron/Coherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    JPsiToDaug.add("JPsiToDaug/Electron/NotCoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToDaug.add("JPsiToDaug/Electron/NotCoherent/XnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToDaug.add("JPsiToDaug/Electron/NotCoherent/OnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToDaug.add("JPsiToDaug/Electron/NotCoherent/XnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToDaug.add("JPsiToDaug/Electron/NotCoherent/OnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToDaug.addClone("JPsiToDaug/Electron/Coherent/XnXn/", "JPsiToDaug/Electron/Coherent/OnXn/");
    JPsiToDaug.addClone("JPsiToDaug/Electron/Coherent/XnXn/", "JPsiToDaug/Electron/Coherent/OnOn/");
    JPsiToDaug.addClone("JPsiToDaug/Electron/Coherent/XnXn/", "JPsiToDaug/Electron/Coherent/XnOn/");
    JPsiToDaug.addClone("JPsiToDaug/Electron/Coherent/", "JPsiToDaug/Electron/Incoherent/");
    JPsiToDaug.addClone("JPsiToDaug/Electron/", "JPsiToDaug/Muon/");
    JPsiToDaug.addClone("JPsiToDaug/Electron/", "JPsiToDaug/Proton/");

    // Correlation histograms
    Correlation.add("Correlation/Muon/Coherent/AccoplAngle", "AccoplAngle", {HistType::kTH1F, {axisAccAngle}});
    Correlation.add("Correlation/Muon/Coherent/CosTheta", "CosTheta", {HistType::kTH1F, {axisAngTheta}});
    Correlation.add("Correlation/Muon/Coherent/Phi", "Phi", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Coherent/Phi1", "Phi1", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Coherent/Phi2", "Phi2", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Coherent/CosThetaPhi", "CosThetaPhi", {HistType::kTH2F, {{axisAngTheta}, {axisPhi}}});
    Correlation.addClone("Correlation/Muon/Coherent/", "Correlation/Muon/Incoherent/");
    Correlation.addClone("Correlation/Muon/", "Correlation/Electron/");
    Correlation.addClone("Correlation/Muon/", "Correlation/Proton/");

    // Asymmetry histograms
    Asymmetry.add("Asymmetry/Muon/Coherent/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/XnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/XnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.addClone("Asymmetry/Muon/Coherent/XnXn/", "Asymmetry/Muon/Coherent/OnXn/");
    Asymmetry.addClone("Asymmetry/Muon/Coherent/XnXn/", "Asymmetry/Muon/Coherent/OnOn/");
    Asymmetry.addClone("Asymmetry/Muon/Coherent/XnXn/", "Asymmetry/Muon/Coherent/XnOn/");
    Asymmetry.addClone("Asymmetry/Muon/Coherent/", "Asymmetry/Muon/Incoherent/");
    Asymmetry.addClone("Asymmetry/Muon/", "Asymmetry/Electron/");

    // MC histograms
    MC.add("MC/hNumberOfMCCollisions", "hNumberOfCollisions", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/hNumberOfRecoCollisions", "hNumberOfRecoCollisions", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/hNumberOfTrueCollisions", "hNumberOfTrueCollisions", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/Muon/hNumberOfMCTracks", "hNumberOfMCTracks", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/hNumberOfMatchedMCTracks", "hNumberOfMatchedMCTracks", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/hPosZ", "hPosZ", {HistType::kTH1F, {{60, -15, 15}}});
    MC.add("MC/hMothersPdg", "hMothersPdg", {HistType::kTH1F, {{1000, -500, 500}}});
    MC.add("MC/hPdg", "hPdg", {HistType::kTH1F, {{1000, -500, 500}}});
    MC.add("MC/Muon/hEta1", "hEta1", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Muon/hEta2", "hEta2", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Muon/hEta", "hEta", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Muon/hPhi1", "hPhi1", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Muon/hPhi2", "hPhi2", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Muon/hPhi", "hPhi", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Muon/hIVM", "hIVM", {HistType::kTH1F, {axisIVM}});
    MC.add("MC/Muon/hRapidity", "hRapidity", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Muon/hPt1", "hPt1", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Muon/hPt2", "hPt2", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Muon/hPt", "hPt", {HistType::kTH1F, {axisPt}});
    MC.add("MC/hResolution", "hResolution", {HistType::kTH1F, {{100, -10, 10}}});
    MC.add("MC/hResolutionPhi", "hResolutionPhi", {HistType::kTH1F, {{100, -0.01, 0.01}}});
  }

  bool cutITSLayers(uint8_t itsClusterMap) const
  {
    std::vector<std::pair<int8_t, std::array<uint8_t, 3>>> requiredITSHits{};
    requiredITSHits.push_back(std::make_pair(1, std::array<uint8_t, 3>{0, 1, 2})); // at least one hit in the innermost layer
    constexpr uint8_t bit = 1;
    for (auto& itsRequirement : requiredITSHits) {
      auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return itsClusterMap & (bit << requiredLayer); });

      if ((itsRequirement.first == -1) && (hits > 0)) {
        return false; // no hits were required in specified layers
      } else if (hits < itsRequirement.first) {
        return false; // not enough hits found in specified layers
      }
    }
    return true;
  }

  template <typename T>
  bool GoodTrackCuts(T const& track)
  {
    // choose only PV contributors
    if (!track.isPVContributor()) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(1);
      return false;
    }
    // pT cut to choose only tracks contributing to J/psi
    if (track.pt() < cutPtTrack) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(2);
      return false;
    }
    // acceptance cut (TPC)
    if (std::abs(RecoDecay::eta(std::array{track.px(), track.py(), track.pz()})) > EtaCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(3);
      return false;
    }
    // DCA
    if (std::abs(track.dcaZ()) > dcaZCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(4);
      return false;
    }
    if (DCAcut) {
      float dcaXYPtCut = 0.0105f + 0.0350f / pow(track.pt(), 1.1f);
      if (std::abs(track.dcaXY()) > dcaXYPtCut) {
        Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(5);
        return false;
      }
    } else {
      if (std::abs(track.dcaXY()) > dcaXYCut) {
        Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(5);
        return false;
      }
    }
    // ITS
    if (!track.hasITS()) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(6);
      return false;
    }
    if (track.itsNCls() < ITSNClsCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(7);
      return false;
    }
    if (!cutITSLayers(track.itsClusterMap())) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(8);
      return false;
    }
    if (track.itsChi2NCl() > ITSChi2NClsCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(9);
      return false;
    }
    //  TPC
    if (!track.hasTPC()) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(10);
      return false;
    }
    if (track.tpcNClsCrossedRows() < TPCNClsCrossedRowsCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(11);
      return false;
    }
    if (track.tpcChi2NCl() > TPCChi2NCls) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(12);
      return false; // TPC chi2
    }
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < TPCMinNCls) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(13);
      return false;
    }
    if (newCutTPC) {
      if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < TPCCrossedOverFindable) {
        Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(14);
        return false;
      }
    }

    return true;
  }

  bool CandidateCuts(float massJpsi, float rapJpsi)
  {
    if (std::abs(rapJpsi) > RapCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(15);
      return false;
    }

    if (massJpsi < 2.0f) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(16);
      return false;
    }

    return true;
  }

  static constexpr std::string_view DaugType[3] = {"Electron/", "Muon/", "Proton/"};
  static constexpr std::string_view ProdType[2] = {"Coherent/", "Incoherent/"};
  static constexpr std::string_view NeutronClass[5] = {"", "XnXn/", "OnXn/", "OnOn/", "XnOn/"};

  template <int daug, typename T>
  void fillSelections(const T& leadingP, const T& subleadingP, float massJpsi, float rapJpsi, std::array<double, 3> mother)
  {
    // selections
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Leading/hITSNCls"))->Fill(massJpsi, leadingP.itsNCls());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Leading/hITSChi2NCls"))->Fill(massJpsi, leadingP.itsChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Leading/hTPCNCls"))->Fill(massJpsi, leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Leading/hTPCChi2NCls"))->Fill(massJpsi, leadingP.tpcChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Leading/hTPCNClsCrossedRows"))->Fill(massJpsi, subleadingP.tpcNClsCrossedRows());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Subleading/hITSNCls"))->Fill(massJpsi, subleadingP.itsNCls());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Subleading/hITSChi2NCls"))->Fill(massJpsi, subleadingP.itsChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Subleading/hTPCNCls"))->Fill(massJpsi, subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Subleading/hTPCChi2NCls"))->Fill(massJpsi, subleadingP.tpcChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Mass/Subleading/hTPCNClsCrossedRows"))->Fill(massJpsi, subleadingP.tpcNClsCrossedRows());

    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Leading/hITSNCls"))->Fill(rapJpsi, leadingP.itsNCls());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Leading/hITSChi2NCls"))->Fill(rapJpsi, leadingP.itsChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Leading/hTPCNCls"))->Fill(rapJpsi, leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Leading/hTPCChi2NCls"))->Fill(rapJpsi, leadingP.tpcChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Leading/hTPCNClsCrossedRows"))->Fill(rapJpsi, subleadingP.tpcNClsCrossedRows());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Subleading/hITSNCls"))->Fill(rapJpsi, subleadingP.itsNCls());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Subleading/hITSChi2NCls"))->Fill(rapJpsi, subleadingP.itsChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Subleading/hTPCNCls"))->Fill(rapJpsi, subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Subleading/hTPCChi2NCls"))->Fill(rapJpsi, subleadingP.tpcChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Rapidity/Subleading/hTPCNClsCrossedRows"))->Fill(rapJpsi, subleadingP.tpcNClsCrossedRows());

    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Leading/hITSNCls"))->Fill(RecoDecay::pt(mother), leadingP.itsNCls());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Leading/hITSChi2NCls"))->Fill(RecoDecay::pt(mother), leadingP.itsChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Leading/hTPCNCls"))->Fill(RecoDecay::pt(mother), leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Leading/hTPCChi2NCls"))->Fill(RecoDecay::pt(mother), leadingP.tpcChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Leading/hTPCNClsCrossedRows"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsCrossedRows());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Subleading/hITSNCls"))->Fill(RecoDecay::pt(mother), subleadingP.itsNCls());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Subleading/hITSChi2NCls"))->Fill(RecoDecay::pt(mother), subleadingP.itsChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Subleading/hTPCNCls"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Subleading/hTPCChi2NCls"))->Fill(RecoDecay::pt(mother), subleadingP.tpcChi2NCl());
    Selections.get<TH2>(HIST("Selections/") + HIST(DaugType[daug]) + HIST("Pt/Subleading/hTPCNClsCrossedRows"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsCrossedRows());
  }

  template <int daug, typename T>
  void fillTGdaug(const T& trkDaughter1, const T& trkDaughter2, std::array<double, 3> daughter1, std::array<double, 3> daughter2)
  {
    TGdaug.get<TH1>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("hTrackPt1"))->Fill(trkDaughter1.pt());
    TGdaug.get<TH1>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("hTrackPt2"))->Fill(trkDaughter2.pt());
    TGdaug.get<TH1>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
    TGdaug.get<TH1>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
    TGdaug.get<TH1>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
    TGdaug.get<TH1>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));

    if (trkDaughter1.hasTPC()) {
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
      if (trkDaughter1.sign() < 0) {
        TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
      } else {
        TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
      }
    }
    if (trkDaughter2.hasTPC()) {
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
    }
    if (trkDaughter1.hasTOF()) {
      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
    }
    if (trkDaughter2.hasTOF()) {

      TGdaug.get<TH2>(HIST("TGdaug/") + HIST(DaugType[daug]) + HIST("PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
    }
  }
  template <int daug, typename T>
  void fillCand(const T& trkDaughter1, const T& trkDaughter2, std::array<double, 3> daughter1, std::array<double, 3> daughter2, std::array<double, 3> mother, float massJpsi, bool XnXn, bool OnOn, bool OnXn, bool XnOn)
  {
    TGdaugCand.get<TH1>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("hTrackPt1"))->Fill(trkDaughter1.pt());
    TGdaugCand.get<TH1>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("hTrackPt2"))->Fill(trkDaughter2.pt());
    TGdaugCand.get<TH1>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
    TGdaugCand.get<TH1>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
    TGdaugCand.get<TH1>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
    TGdaugCand.get<TH1>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));
    TGdaugCand.get<TH1>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("hPairPt"))->Fill(RecoDecay::pt(mother));
    TGdaugCand.get<TH1>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("hPairIVM"))->Fill(massJpsi);

    if (trkDaughter1.hasTPC()) {
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
    }
    if (trkDaughter2.hasTPC()) {
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
    }
    if (trkDaughter1.hasTOF()) {
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
    }
    if (trkDaughter2.hasTOF()) {
      TGdaugCand.get<TH2>(HIST("TGdaugCand/") + HIST(DaugType[daug]) + HIST("PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
    }

    if (RecoDecay::pt(mother) < 0.2f) {
      JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Coherent/hIVM"))->Fill(massJpsi);
      if (XnXn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Coherent/XnXn/hIVM"))->Fill(massJpsi);
      } else if (OnXn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Coherent/OnXn/hIVM"))->Fill(massJpsi);
      } else if (XnOn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Coherent/XnOn/hIVM"))->Fill(massJpsi);
      } else if (OnOn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Coherent/OnOn/hIVM"))->Fill(massJpsi);
      }
    } else if (RecoDecay::pt(mother) > 0.4f) {
      JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("NotCoherent/hIVM"))->Fill(massJpsi);
      if (XnXn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("NotCoherent/XnXn/hIVM"))->Fill(massJpsi);
      } else if (OnXn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("NotCoherent/OnXn/hIVM"))->Fill(massJpsi);
      } else if (XnOn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("NotCoherent/XnOn/hIVM"))->Fill(massJpsi);
      } else if (OnOn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("NotCoherent/OnOn/hIVM"))->Fill(massJpsi);
      }
    }
    if (RecoDecay::pt(mother) > 0.2f) {
      JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Incoherent/hIVM"))->Fill(massJpsi);
      if (XnXn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Incoherent/XnXn/hIVM"))->Fill(massJpsi);
      } else if (OnXn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Incoherent/OnXn/hIVM"))->Fill(massJpsi);
      } else if (XnOn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Incoherent/XnOn/hIVM"))->Fill(massJpsi);
      } else if (OnOn) {
        JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST("Incoherent/OnOn/hIVM"))->Fill(massJpsi);
      }
    }
  }

  template <int daug, int prod, int neutron, typename T>
  void fillPeak(const T& trkDaughter1, const T& trkDaughter2, std::array<double, 3> daughter1, std::array<double, 3> daughter2, std::array<double, 3> mother, float rapJpsi)
  {
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hPt1"))->Fill(trkDaughter1.pt());
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hPt2"))->Fill(trkDaughter2.pt());
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hEta1"))->Fill(RecoDecay::eta(daughter1));
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hEta2"))->Fill(RecoDecay::eta(daughter2));
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hPhi1"))->Fill(RecoDecay::phi(daughter1));
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hPhi2"))->Fill(RecoDecay::phi(daughter2));
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hPt"))->Fill(RecoDecay::pt(mother));
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hEta"))->Fill(RecoDecay::eta(mother));
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hRap"))->Fill(rapJpsi);
    JPsiToDaug.get<TH1>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("hPhi"))->Fill(RecoDecay::phi(mother));

    if (trkDaughter1.hasTPC()) {
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
    }
    if (trkDaughter2.hasTPC()) {
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
    }
    if (trkDaughter1.hasTOF()) {
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
    }
    if (trkDaughter2.hasTOF()) {
      JPsiToDaug.get<TH2>(HIST("JPsiToDaug/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
    }
  }

  template <int daug, int prod, int neutron>
  void fillCorrAsy(std::array<double, 3> daughter1, std::array<double, 3> daughter2, double dp, double dpRandom, float* q)
  {
    if (neutron == 0) {
      Correlation.get<TH1>(HIST("Correlation/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("Phi1"))->Fill(RecoDecay::phi(daughter1), 1.);
      Correlation.get<TH1>(HIST("Correlation/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("Phi2"))->Fill(RecoDecay::phi(daughter2), 1.);
      Correlation.get<TH1>(HIST("Correlation/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("Phi"))->Fill(q[1], 1.);
      Correlation.get<TH1>(HIST("Correlation/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("CosTheta"))->Fill(q[2], 1.);
      Correlation.get<TH1>(HIST("Correlation/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("AccoplAngle"))->Fill(q[0], 1.);
      Correlation.get<TH2>(HIST("Correlation/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST("CosThetaPhi"))->Fill(q[2], q[1]);
    }

    Asymmetry.get<TH1>(HIST("Asymmetry/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("DeltaPhi"))->Fill(dp);
    Asymmetry.get<TH1>(HIST("Asymmetry/") + HIST(DaugType[daug]) + HIST(ProdType[prod]) + HIST(NeutronClass[neutron]) + HIST("DeltaPhiRandom"))->Fill(dpRandom);
  }

  template <typename C, typename Ts>
  void fillHistograms(C collision, Ts tracks)
  {
    Statistics.get<TH1>(HIST("Statistics/hCutCounterCollisions"))->Fill(0); // number of collisions without any cuts

    // check UPC vs standard
    if (doOnlyUPC) {
      if (collision.flags() == 0) {
        return;
      }
    }
    if (doOnlyStandard) {
      if (collision.flags() == 1) {
        return;
      }
    }

    if (std::abs(collision.posZ()) > cutVertexZ) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterCollisions"))->Fill(1);
      return;
    }

    // distinguish ZDC classes
    bool XnXn = false, OnOn = false, XnOn = false, OnXn = false;
    int neutronClass = -1;
    if (collision.energyCommonZNA() < ZNenergyCut && collision.energyCommonZNC() < ZNenergyCut) {
      OnOn = true;
      neutronClass = 0;
    }
    if (collision.energyCommonZNA() > ZNenergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() > ZNenergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) {
      XnXn = true;
      neutronClass = 1;
    }
    if (collision.energyCommonZNA() > ZNenergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() < ZNenergyCut) {
      XnOn = true;
      neutronClass = 2;
    }
    if (collision.energyCommonZNA() < ZNenergyCut && collision.energyCommonZNC() > ZNenergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) {
      OnXn = true;
      neutronClass = 3;
    }

    // loop over tracks without selections
    for (auto& track : tracks) {
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();

      Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(0);
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(0);
      if (track.isPVContributor() == 1) {
        Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(1);
        PVContributors.get<TH1>(HIST("PVContributors/hTrackPt"))->Fill(track.pt());
        PVContributors.get<TH1>(HIST("PVContributors/hTrackEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}));
        PVContributors.get<TH1>(HIST("PVContributors/hTrackPhi"))->Fill(RecoDecay::phi(trkPx, trkPy));

        if (track.hasTPC()) {
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.tpcSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTPCVsPt"))->Fill(track.pt(), track.tpcSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTPCVsEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}), track.tpcSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(trkPx, trkPy), track.tpcSignal());
          PVContributors.get<TH1>(HIST("PVContributors/hTPCNClsFindable"))->Fill(track.tpcNClsFindable());
          PVContributors.get<TH1>(HIST("PVContributors/hTPCNClsFindableMinusFound"))->Fill(track.tpcNClsFindableMinusFound());
        }

        if (track.hasTOF()) {
          PVContributors.get<TH2>(HIST("PVContributors/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.beta());
        }
      }

      RawData.get<TH1>(HIST("RawData/hTrackPt"))->Fill(track.pt());
      RawData.get<TH1>(HIST("RawData/hTrackEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}));
      RawData.get<TH1>(HIST("RawData/hTrackPhi"))->Fill(RecoDecay::phi(trkPx, trkPy));
      RawData.get<TH2>(HIST("RawData/hTrackDCAXYZ"))->Fill(track.dcaXY(), track.dcaZ());
      RawData.get<TH1>(HIST("RawData/hTPCNClsFindable"))->Fill(track.tpcNClsFindable());
      RawData.get<TH1>(HIST("RawData/hTPCNClsFindableMinusFound"))->Fill(track.tpcNClsFindableMinusFound());
      RawData.get<TH1>(HIST("RawData/hITSNCls"))->Fill(track.itsNCls());
      RawData.get<TH1>(HIST("RawData/hTPCNCls"))->Fill(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
      RawData.get<TH1>(HIST("RawData/hITSChi2NCls"))->Fill(track.itsChi2NCl());
      RawData.get<TH1>(HIST("RawData/hTPCChi2NCls"))->Fill(track.tpcChi2NCl());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTotFT0AmplitudeA"))->Fill(collision.totalFT0AmplitudeA());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTotFT0AmplitudeC"))->Fill(collision.totalFT0AmplitudeC());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTotFV0AmplitudeA"))->Fill(collision.totalFV0AmplitudeA());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTotFDDAmplitudeA"))->Fill(collision.totalFDDAmplitudeA());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTotFDDAmplitudeC"))->Fill(collision.totalFDDAmplitudeC());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTimeFT0A"))->Fill(collision.timeFT0A());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTimeFT0C"))->Fill(collision.timeFT0C());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTimeFV0A"))->Fill(collision.timeFV0A());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTimeFDDA"))->Fill(collision.timeFDDA());
      RawData.get<TH1>(HIST("RawData/QA/FIT/hTimeFDDC"))->Fill(collision.timeFDDC());

      if (track.hasTPC()) {
        RawData.get<TH2>(HIST("RawData/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.tpcSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTPCVsPt"))->Fill(track.pt(), track.tpcSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTPCVsEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}), track.tpcSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(trkPx, trkPy), track.tpcSignal());
      }

      if (track.hasTOF()) {
        RawData.get<TH2>(HIST("RawData/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.beta());
      }
    }

    int countGT = 0;
    std::vector<int> trkIdx;
    // loop over tracks with selections

    RawData.get<TH1>(HIST("RawData/hPositionX"))->Fill(collision.posX());
    RawData.get<TH1>(HIST("RawData/hPositionY"))->Fill(collision.posY());
    RawData.get<TH1>(HIST("RawData/hPositionZ"))->Fill(collision.posZ());
    RawData.get<TH2>(HIST("RawData/hPositionXY"))->Fill(collision.posX(), collision.posY());
    RawData.get<TH1>(HIST("RawData/QA/ZDC/hZNACommonEnergy"))->Fill(collision.energyCommonZNA());
    RawData.get<TH1>(HIST("RawData/QA/ZDC/hZNCCommonEnergy"))->Fill(collision.energyCommonZNC());
    RawData.get<TH2>(HIST("RawData/QA/ZDC/hZNAvsZNCCommonEnergy"))->Fill(collision.energyCommonZNA(), collision.energyCommonZNC());
    RawData.get<TH1>(HIST("RawData/QA/ZDC/hZNCTime"))->Fill(collision.timeZNC());
    RawData.get<TH1>(HIST("RawData/QA/ZDC/hZNATime"))->Fill(collision.timeZNA());
    RawData.get<TH2>(HIST("RawData/QA/ZDC/hZNAvsZNCTime"))->Fill(collision.timeZNA(), collision.timeZNC());

    for (auto& track : tracks) {

      Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(2.);

      // select good tracks
      if (GoodTrackCuts(track) != 1) {
        continue;
      }

      countGT++;
      trkIdx.push_back(track.index());
    }

    Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(3., countGT);
    Statistics.get<TH1>(HIST("Statistics/hNumberGT"))->Fill(countGT);

    float massEl = o2::constants::physics::MassElectron;
    float massMu = o2::constants::physics::MassMuonMinus;
    float massPr = o2::constants::physics::MassProton;

    if (countGT == 2) {
      TLorentzVector mom, daughter[2];
      auto trkDaughter1 = tracks.iteratorAt(trkIdx[0]);
      auto trkDaughter2 = tracks.iteratorAt(trkIdx[1]);

      if ((trkDaughter1.sign() * trkDaughter2.sign()) != -1) {
        return;
      }

      if (chargeOrdered) {
        if (tracks.iteratorAt(trkIdx[0]).sign() < 0) {
          trkDaughter1 = tracks.iteratorAt(trkIdx[0]);
          trkDaughter2 = tracks.iteratorAt(trkIdx[1]);
        } else {
          trkDaughter1 = tracks.iteratorAt(trkIdx[1]);
          trkDaughter2 = tracks.iteratorAt(trkIdx[0]);
        }
      }

      std::array<double, 3> daughter1 = {trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()};
      std::array<double, 3> daughter2 = {trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()};

      TG.get<TH1>(HIST("TG/hTrackPt1"))->Fill(trkDaughter1.pt());
      TG.get<TH1>(HIST("TG/hTrackPt2"))->Fill(trkDaughter2.pt());
      TG.get<TH1>(HIST("TG/hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
      TG.get<TH1>(HIST("TG/hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
      TG.get<TH1>(HIST("TG/hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
      TG.get<TH1>(HIST("TG/hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));

      if (trkDaughter1.hasTPC()) {
        TG.get<TH2>(HIST("TG/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
        TG.get<TH2>(HIST("TG/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
        TG.get<TH2>(HIST("TG/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
        TG.get<TH2>(HIST("TG/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
        TG.get<TH1>(HIST("TG/hTPCNClsFindable"))->Fill(trkDaughter1.tpcNClsFindable());
        TG.get<TH1>(HIST("TG/hTPCNClsFindableMinusFound"))->Fill(trkDaughter1.tpcNClsFindableMinusFound());
        if (trkDaughter1.sign() < 0) {
          TG.get<TH2>(HIST("TG/TPC/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
        } else {
          TG.get<TH2>(HIST("TG/TPC/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
        }
      }
      if (trkDaughter2.hasTPC()) {
        TG.get<TH2>(HIST("TG/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
        TG.get<TH2>(HIST("TG/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
        TG.get<TH2>(HIST("TG/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
        TG.get<TH2>(HIST("TG/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
        TG.get<TH1>(HIST("TG/hTPCNClsFindable"))->Fill(trkDaughter2.tpcNClsFindable());
        TG.get<TH1>(HIST("TG/hTPCNClsFindableMinusFound"))->Fill(trkDaughter2.tpcNClsFindableMinusFound());
      }
      if (trkDaughter1.hasTOF()) {
        TG.get<TH2>(HIST("TG/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
      }
      if (trkDaughter2.hasTOF()) {
        TG.get<TH2>(HIST("TG/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
      }

      if (doElectrons) {
        if (RecoDecay::sumOfSquares(trkDaughter1.tpcNSigmaMu(), trkDaughter2.tpcNSigmaMu()) > RecoDecay::sumOfSquares(trkDaughter1.tpcNSigmaEl(), trkDaughter2.tpcNSigmaEl())) {

          auto ene1 = RecoDecay::e(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), massEl);
          auto ene2 = RecoDecay::e(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), massEl);
          daughter[0].SetPxPyPzE(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), ene1);
          daughter[1].SetPxPyPzE(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), ene2);
          mom = daughter[0] + daughter[1];

          std::array<double, 3> mother = {trkDaughter1.px() + trkDaughter2.px(), trkDaughter1.py() + trkDaughter2.py(), trkDaughter1.pz() + trkDaughter2.pz()};

          auto arrMom = std::array{daughter1, daughter2};
          float massJpsi = RecoDecay::m(arrMom, std::array{massEl, massEl});
          float rapJpsi = RecoDecay::y(mother, massJpsi);

          auto leadingP = RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()) > RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()) ? trkDaughter1 : trkDaughter2;
          auto subleadingP = (leadingP == trkDaughter1) ? trkDaughter1 : trkDaughter2;

          // TGdaug
          fillTGdaug<0>(trkDaughter1, trkDaughter2, daughter1, daughter2);

          // selections
          fillSelections<0>(leadingP, subleadingP, massJpsi, rapJpsi, mother);

          // TGdaugCand
          fillCand<0>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, massJpsi, XnXn, OnOn, OnXn, XnOn);

          double dp = DeltaPhi(daughter[0], daughter[1]);
          float* q = correlation(&daughter[0], &daughter[1], &mom);
          double dpRandom = DeltaPhiRandom(daughter[0], daughter[1]);

          if ((massJpsi < maxJpsiMass) && (massJpsi > minJpsiMass)) {
            TGdaugCand.get<TH1>(HIST("TGdaugCand/Electron/hJpsiPt"))->Fill(RecoDecay::pt(mother));
            TGdaugCand.get<TH1>(HIST("TGdaugCand/Electron/hJpsiRap"))->Fill(rapJpsi);

            if (trkDaughter1.sign() < 0) {
              TGdaugCand.get<TH2>(HIST("TGdaugCand/Electron/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
            } else {
              TGdaugCand.get<TH2>(HIST("TGdaugCand/Electron/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
            }

            if (RecoDecay::pt(mother) < 0.2f) {
              fillPeak<0, 0, 0>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);

              if (XnXn) {
                fillPeak<0, 0, 1>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<0, 0, 1>(daughter1, daughter2, dp, dpRandom, q);
              } else if (OnOn) {
                fillPeak<0, 0, 3>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<0, 0, 3>(daughter1, daughter2, dp, dpRandom, q);
              } else if (XnOn) {
                fillPeak<0, 0, 4>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<0, 0, 4>(daughter1, daughter2, dp, dpRandom, q);
              } else if (OnXn) {
                fillPeak<0, 0, 2>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<0, 0, 2>(daughter1, daughter2, dp, dpRandom, q);
              }
              fillCorrAsy<0, 0, 0>(daughter1, daughter2, dp, dpRandom, q);

            } // end coherent electrons
            if (RecoDecay::pt(mother) > 0.2f) {
              fillPeak<0, 1, 0>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
              fillCorrAsy<0, 1, 0>(daughter1, daughter2, dp, dpRandom, q);

              if (XnXn) {
                fillPeak<0, 1, 1>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
              } else if (OnOn) {
                fillPeak<0, 1, 3>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<0, 1, 3>(daughter1, daughter2, dp, dpRandom, q);
              } else if (XnOn) {
                fillPeak<0, 1, 4>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<0, 1, 4>(daughter1, daughter2, dp, dpRandom, q);
              } else if (OnXn) {
                fillPeak<0, 1, 2>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<0, 1, 2>(daughter1, daughter2, dp, dpRandom, q);
              }

            } // end incoherent electrons
          } // end mass cut
          delete[] q;
        }
      } // end electrons
      if (doMuons) {
        if (RecoDecay::sumOfSquares(trkDaughter1.tpcNSigmaEl(), trkDaughter2.tpcNSigmaEl()) > RecoDecay::sumOfSquares(trkDaughter1.tpcNSigmaMu(), trkDaughter2.tpcNSigmaMu())) {

          auto ene1 = RecoDecay::e(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), massMu);
          auto ene2 = RecoDecay::e(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), massMu);
          daughter[0].SetPxPyPzE(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), ene1);
          daughter[1].SetPxPyPzE(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), ene2);
          mom = daughter[0] + daughter[1];

          std::array<double, 3> mother = {trkDaughter1.px() + trkDaughter2.px(), trkDaughter1.py() + trkDaughter2.py(), trkDaughter1.pz() + trkDaughter2.pz()};

          auto arrMom = std::array{daughter1, daughter2};
          float massJpsi = RecoDecay::m(arrMom, std::array{massMu, massMu});
          float rapJpsi = RecoDecay::y(mother, massJpsi);

          auto leadingP = RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()) > RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()) ? trkDaughter1 : trkDaughter2;
          auto subleadingP = (leadingP == trkDaughter1) ? trkDaughter1 : trkDaughter2;

          // TGdaug
          fillTGdaug<1>(trkDaughter1, trkDaughter2, daughter1, daughter2);
          // selections
          fillSelections<1>(leadingP, subleadingP, massJpsi, rapJpsi, mother);

          if (CandidateCuts(massJpsi, rapJpsi) != 1) {
            return;
          }

          // TGdaugCand
          fillCand<1>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, massJpsi, XnXn, OnOn, OnXn, XnOn);

          double dp = DeltaPhi(daughter[0], daughter[1]);
          float* q = correlation(&daughter[0], &daughter[1], &mom);
          double dpRandom = DeltaPhiRandom(daughter[0], daughter[1]);
          // fill tree with muon J/Psi candidates
          Tree(collision.runNumber(), collision.globalBC(),
               collision.posX(), collision.posY(), collision.posZ(),
               collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
               collision.timeFT0A(), collision.timeFT0C(), collision.timeFV0A(), collision.timeFDDA(), collision.timeFDDC(),
               collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC(), neutronClass,
               RecoDecay::pt(mother), RecoDecay::eta(mother), RecoDecay::phi(mother), rapJpsi, massJpsi, dpRandom, dp, q[1], q[2],
               trkDaughter1.sign(), trkDaughter1.pt(), RecoDecay::eta(daughter1), RecoDecay::phi(daughter1),
               trkDaughter2.sign(), trkDaughter2.pt(), RecoDecay::eta(daughter2), RecoDecay::phi(daughter2));

          if ((massJpsi < maxJpsiMass) && (massJpsi > minJpsiMass)) {
            TGdaugCand.get<TH1>(HIST("TGdaugCand/Muon/hJpsiPt"))->Fill(RecoDecay::pt(mother));
            TGdaugCand.get<TH1>(HIST("TGdaugCand/Muon/hJpsiRap"))->Fill(rapJpsi);
            if (trkDaughter1.sign() < 0) {
              TGdaugCand.get<TH2>(HIST("TGdaugCand/Muon/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
            } else {
              TGdaugCand.get<TH2>(HIST("TGdaugCand/Muon/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
            }

            if (RecoDecay::pt(mother) < 0.2f) {
              fillPeak<1, 0, 0>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
              fillCorrAsy<1, 0, 0>(daughter1, daughter2, dp, dpRandom, q);

              if (XnXn) {
                fillPeak<1, 0, 1>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<1, 0, 1>(daughter1, daughter2, dp, dpRandom, q);
              } else if (OnOn) {
                fillPeak<1, 0, 3>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<1, 0, 3>(daughter1, daughter2, dp, dpRandom, q);
              } else if (XnOn) {
                fillPeak<1, 0, 4>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<1, 0, 4>(daughter1, daughter2, dp, dpRandom, q);
              } else if (OnXn) {
                fillPeak<1, 0, 2>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<1, 0, 2>(daughter1, daughter2, dp, dpRandom, q);
              }
            }
            if (RecoDecay::pt(mother) > 0.2f) {
              fillPeak<1, 1, 0>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
              fillCorrAsy<1, 1, 0>(daughter1, daughter2, dp, dpRandom, q);

              if (XnXn) {
                fillPeak<1, 1, 1>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<1, 1, 1>(daughter1, daughter2, dp, dpRandom, q);
              } else if (OnOn) {
                fillPeak<1, 1, 3>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<1, 1, 3>(daughter1, daughter2, dp, dpRandom, q);
              } else if (XnOn) {
                fillPeak<1, 1, 4>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<1, 1, 4>(daughter1, daughter2, dp, dpRandom, q);
              } else if (OnXn) {
                fillPeak<1, 1, 2>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
                fillCorrAsy<1, 1, 2>(daughter1, daughter2, dp, dpRandom, q);
              }
            } // end incoherent
          } // end mass cut
          delete[] q;
        }
      } // end muons
      if (doProtons) {
        if (RecoDecay::sumOfSquares((trkDaughter1.hasTOF() ? trkDaughter1.tofNSigmaPr() : trkDaughter1.tpcNSigmaPr()), (trkDaughter2.hasTOF() ? trkDaughter2.tofNSigmaPr() : trkDaughter2.tpcNSigmaPr()) < 4)) {

          auto ene1 = RecoDecay::e(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), massPr);
          auto ene2 = RecoDecay::e(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), massPr);
          daughter[0].SetPxPyPzE(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), ene1);
          daughter[1].SetPxPyPzE(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), ene2);
          mom = daughter[0] + daughter[1];

          std::array<double, 3> mother = {trkDaughter1.px() + trkDaughter2.px(), trkDaughter1.py() + trkDaughter2.py(), trkDaughter1.pz() + trkDaughter2.pz()};

          if (TOFBothProtons) {
            if (!trkDaughter1.hasTOF() || !trkDaughter2.hasTOF())
              return;
          }
          if (TOFOneProton) {
            if ((trkDaughter1.hasTOF() && trkDaughter2.hasTOF()) || (!trkDaughter1.hasTOF() && !trkDaughter2.hasTOF()))
              return;
          }
          if (TOFAtLeastOneProton) {
            if (!trkDaughter1.hasTOF() && !trkDaughter2.hasTOF())
              return;
          }

          auto arrMom = std::array{daughter1, daughter2};
          float massJpsi = RecoDecay::m(arrMom, std::array{massPr, massPr});
          float rapJpsi = RecoDecay::y(mother, massJpsi);

          // TGdaug
          fillTGdaug<2>(trkDaughter1, trkDaughter2, daughter1, daughter2);

          if (CandidateCuts(massJpsi, rapJpsi) != 1) {
            return;
          }

          // TGdaugCand
          fillCand<2>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, massJpsi, XnXn, OnOn, OnXn, XnOn);

          if (RecoDecay::pt(mother) < 0.2f) {
            JPsiToDaug.get<TH1>(HIST("JPsiToDaug/Proton/Coherent/hIVM"))->Fill(massJpsi);
          } else {
            JPsiToDaug.get<TH1>(HIST("JPsiToDaug/Proton/Incoherent/hIVM"))->Fill(massJpsi);
          }

          if (massJpsi < maxJpsiMass && massJpsi > minJpsiMass) {
            TGdaugCand.get<TH1>(HIST("TGdaugCand/Proton/hJpsiPt"))->Fill(RecoDecay::pt(mother));
            TGdaugCand.get<TH1>(HIST("TGdaugCand/Proton/hJpsiRap"))->Fill(rapJpsi);
            if (RecoDecay::pt(mother) < 0.2f) {

              fillPeak<2, 0, 0>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);

            } // end coherent
            if (RecoDecay::pt(mother) > 0.2f) {
              // fill track histos
              fillPeak<1, 1, 0>(trkDaughter1, trkDaughter2, daughter1, daughter2, mother, rapJpsi);
            } // end incoherent
          } // end mass cut
        }
      } // end protons
    } // end two tracks
  } // end reco process

  template <typename C, typename T>
  void processMC(C const& mcCollision, T const& mcParticles)
  {

    MC.get<TH1>(HIST("MC/hNumberOfMCCollisions"))->Fill(1.);
    MC.get<TH1>(HIST("MC/hPosZ"))->Fill(mcCollision.posZ());

    std::array<float, 3> daughPart1Mu;
    std::array<float, 3> daughPart2Mu;
    std::array<float, 3> motherPart;
    float energyMother;
    float daughPart1pdg;
    float daughPart2pdg;

    // fill number of particles
    for (auto const& mcParticle : mcParticles) {
      MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(1.);
      MC.get<TH1>(HIST("MC/hPdg"))->Fill(mcParticle.pdgCode());
      if (mcParticle.has_daughters()) {
        MC.get<TH1>(HIST("MC/hMothersPdg"))->Fill(mcParticle.pdgCode());
        MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(3.);
        if (mcParticle.pdgCode() == 443) {
          MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(4.);
          int count = 0;
          for (auto& daughter : mcParticle.template daughters_as<T>()) {
            if ((daughter.pdgCode() == -13) && daughter.isPhysicalPrimary()) {
              daughPart1Mu = {daughter.px(), daughter.py(), daughter.pz()};
              daughPart1pdg = daughter.pdgCode();

              MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(5.);
              count++;
            }
            if ((daughter.pdgCode() == 13) && daughter.isPhysicalPrimary()) {
              daughPart2Mu = {daughter.px(), daughter.py(), daughter.pz()};
              daughPart2pdg = daughter.pdgCode();

              MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(6.);
              count++;
            }
          }
          if (count == 2) {
            MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(7.);
            motherPart = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};
            energyMother = mcParticle.e();
          }
        }
      }
    }

    // calculate needed distributions

    MC.get<TH1>(HIST("MC/Muon/hEta1"))->Fill(RecoDecay::eta(daughPart1Mu));
    MC.get<TH1>(HIST("MC/Muon/hEta2"))->Fill(RecoDecay::eta(daughPart2Mu));
    MC.get<TH1>(HIST("MC/Muon/hEta"))->Fill(RecoDecay::eta(motherPart));
    MC.get<TH1>(HIST("MC/Muon/hPhi1"))->Fill(RecoDecay::phi(daughPart1Mu));
    MC.get<TH1>(HIST("MC/Muon/hPhi2"))->Fill(RecoDecay::phi(daughPart2Mu));
    MC.get<TH1>(HIST("MC/Muon/hPhi"))->Fill(RecoDecay::phi(motherPart));
    MC.get<TH1>(HIST("MC/Muon/hPt1"))->Fill(RecoDecay::pt(daughPart1Mu));
    MC.get<TH1>(HIST("MC/Muon/hPt2"))->Fill(RecoDecay::pt(daughPart2Mu));
    MC.get<TH1>(HIST("MC/Muon/hPt"))->Fill(RecoDecay::pt(motherPart));
    MC.get<TH1>(HIST("MC/Muon/hIVM"))->Fill(RecoDecay::m(motherPart, energyMother));
    MC.get<TH1>(HIST("MC/Muon/hRapidity"))->Fill(RecoDecay::y(motherPart, RecoDecay::m(motherPart, energyMother)));

    float massMu = o2::constants::physics::MassMuonMinus;
    auto ene1 = RecoDecay::e(daughPart1Mu, massMu);
    auto ene2 = RecoDecay::e(daughPart2Mu, massMu);

    TLorentzVector mom, daughter[2];
    daughter[0].SetPxPyPzE(daughPart1Mu[0], daughPart1Mu[1], daughPart1Mu[2], ene1);
    daughter[1].SetPxPyPzE(daughPart2Mu[0], daughPart2Mu[1], daughPart2Mu[2], ene2);
    mom = daughter[0] + daughter[1];

    double dp = DeltaPhi(daughter[0], daughter[1]);
    float* q = correlation(&daughter[0], &daughter[1], &mom);
    double dpRandom = DeltaPhiRandom(daughter[0], daughter[1]);

    // fill tree with muon J/Psi candidates
    TreeMC(mcCollision.globalBC(),
           RecoDecay::pt(motherPart), RecoDecay::eta(motherPart), RecoDecay::phi(motherPart), RecoDecay::y(motherPart, RecoDecay::m(motherPart, energyMother)), RecoDecay::m(motherPart, energyMother), dpRandom, dp, q[1], q[2],
           daughPart1pdg, RecoDecay::pt(daughPart1Mu), RecoDecay::eta(daughPart1Mu), RecoDecay::phi(daughPart1Mu),
           daughPart2pdg, RecoDecay::pt(daughPart2Mu), RecoDecay::eta(daughPart2Mu), RecoDecay::phi(daughPart2Mu));

  } // end MC skimmed process

  void processMCFull(MCUDCollisionFull const& collision, MCUDTracksFull const& tracks, aod::UDMcCollisions const& mcCollision, aod::UDMcParticles const&)
  {
    MC.get<TH1>(HIST("MC/hNumberOfRecoCollisions"))->Fill(collision.size());
    MC.get<TH1>(HIST("MC/hNumberOfTrueCollisions"))->Fill(mcCollision.size());
    MC.get<TH1>(HIST("MC/hNumberOfMatchedMCTracks"))->Fill(1.);
    if (collision.has_udMcCollision()) {
      std::array<float, 3> recoTrack;
      std::array<float, 3> truePart;
      for (auto& track : tracks) {
        MC.get<TH1>(HIST("MC/hNumberOfMatchedMCTracks"))->Fill(1.);
        if (track.has_udMcParticle()) {
          MC.get<TH1>(HIST("MC/hNumberOfMatchedMCTracks"))->Fill(2.);
          auto mcParticle = track.udMcParticle();
          MC.fill(HIST("MC/hResolution"), mcParticle.px() - track.px());
          if (std::abs(mcParticle.pdgCode()) == 13 && mcParticle.isPhysicalPrimary()) {
            MC.get<TH1>(HIST("MC/hNumberOfMatchedMCTracks"))->Fill(3.);
            if (mcParticle.has_mothers()) {
              MC.get<TH1>(HIST("MC/hNumberOfMatchedMCTracks"))->Fill(4.);
              auto const& mother = mcParticle.mothers_first_as<aod::UDMcParticles>();
              if (mother.pdgCode() == 443) {
                MC.get<TH1>(HIST("MC/hNumberOfMatchedMCTracks"))->Fill(5.);
                recoTrack = {track.px(), track.py(), track.pz()};
                truePart = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};
              }
            }
          }
          MC.fill(HIST("MC/hResolutionPhi"), RecoDecay::phi(recoTrack) - RecoDecay::phi(truePart));
        }
      }
    }
  } // end MC Full process

  void processDGrecoLevel(UDCollisionFull const& collision, UDTracksFull const& tracks)
  {
    fillHistograms(collision, tracks);
  } // end DG process

  void processSGrecoLevel(SGUDCollisionFull const& collision, UDTracksFull const& tracks)
  {
    int gapSide = collision.gapSide();
    int trueGapSide = sgSelector.trueGap(collision, cutMyGapSideFV0, cutMyGapSideFT0A, cutMyGapSideFT0C, cutMyGapSideZDC);
    if (useTrueGap) {
      gapSide = trueGapSide;
    }
    if (gapSide != whichGapSide) {
      return;
    }
    fillHistograms(collision, tracks);

  } // end SG process

  void processMCtruth(aod::UDMcCollision const& mcCollision, aod::UDMcParticles const& mcParticles)
  {
    processMC(mcCollision, mcParticles);
  }

  PROCESS_SWITCH(UpcJpsiCentralBarrel, processDGrecoLevel, "Iterate over DG skimmed data.", false);
  PROCESS_SWITCH(UpcJpsiCentralBarrel, processSGrecoLevel, "Iterate over SG skimmed data.", true);
  PROCESS_SWITCH(UpcJpsiCentralBarrel, processMCtruth, "Iterate of MC true data.", true);
  PROCESS_SWITCH(UpcJpsiCentralBarrel, processMCFull, "Iterate over both true and reco.", true);
}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcJpsiCentralBarrel>(cfgc, TaskName{"upc-jpsi-corr"}),
  };
}
