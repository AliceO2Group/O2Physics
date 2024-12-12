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

struct UpcJpsiCentralBarrel {
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
  Configurable<float> cutPtTrack{"cutPtTrack", 0.7f, "pT cut per track"};
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

  HistogramRegistry TGmu{
    "TGmu",
    {}};

  HistogramRegistry TGmuCand{
    "TGmuCand",
    {}};

  HistogramRegistry TGel{
    "TGel",
    {}};

  HistogramRegistry TGelCand{
    "TGelCand",
    {}};

  HistogramRegistry TGp{
    "TGp",
    {}};

  HistogramRegistry TGpCand{
    "TGpCand",
    {}};

  HistogramRegistry JPsiToEl{
    "JPsiToEl",
    {}};

  HistogramRegistry JPsiToMu{
    "JPsiToMu",
    {}};

  HistogramRegistry JPsiToP{
    "JPsiToP",
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
    RawData.add("RawData/hZNACommonEnergy", "hZNACommonEnergy", {HistType::kTH1F, {axisZDCEnergy}});
    RawData.add("RawData/hZNCCommonEnergy", "hZNCCommonEnergy", {HistType::kTH1F, {axisZDCEnergy}});
    RawData.add("RawData/hZNAvsZNCCommonEnergy", "hZNAvsZNCCommonEnergy", {HistType::kTH2F, {axisZDCEnergy, axisZDCEnergy}});
    RawData.add("RawData/hZNATime", "hZNATime", {HistType::kTH1F, {axisZDCTime}});
    RawData.add("RawData/hZNCTime", "hZNCTime", {HistType::kTH1F, {axisZDCTime}});
    RawData.add("RawData/hZNAvsZNCTime", "hZNAvsZNCTime", {HistType::kTH2F, {axisZDCTime, axisZDCTime}});
    RawData.add("RawData/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    RawData.add("RawData/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    RawData.add("RawData/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    RawData.add("RawData/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    RawData.add("RawData/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // Selection checks
    Selections.add("Selections/Electron/Mass/Leading/hITSNClsVsM", "hITSNClsVsM", {HistType::kTH2F, {axisIVMSel, axisITSNCls}});
    Selections.add("Selections/Electron/Mass/Leading/hITSChi2NClsVsM", "hITSChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCNClsVsM", "hTPCNClsVsM", {HistType::kTH2F, {axisIVMSel, axisTPCNCls}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCChi2NClsVsM", "hTPCChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCNClsCrossedRowsVsM", "hTPCNClsCrossedRowsVsM", {HistType::kTH2F, {axisIVMSel, axisTPCCrossed}});
    Selections.add("Selections/Electron/Mass/Subleading/hITSNClsVsM", "hITSNClsVsM", {HistType::kTH2F, {axisIVMSel, axisITSNCls}});
    Selections.add("Selections/Electron/Mass/Subleading/hITSChi2NClsVsM", "hITSChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Electron/Mass/Subleading/hTPCNClsVsM", "hTPCNClsVsM", {HistType::kTH2F, {axisIVMSel, axisTPCNCls}});
    Selections.add("Selections/Electron/Mass/Subleading/hTPCChi2NClsVsM", "hTPCChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Electron/Mass/Subleading/hTPCNClsCrossedRowsVsM", "hTPCNClsCrossedRowsVsM", {HistType::kTH2F, {axisIVMSel, axisTPCCrossed}});

    Selections.add("Selections/Electron/Rapidity/Leading/hITSNClsVsY", "hITSNClsVsY", {HistType::kTH2F, {axisEta, axisITSNCls}});
    Selections.add("Selections/Electron/Rapidity/Leading/hITSChi2NClsVsY", "hITSChi2NClsVsY", {HistType::kTH2F, {axisEta, axisChi2}});
    Selections.add("Selections/Electron/Rapidity/Leading/hTPCNClsVsY", "hTPCNClsVsY", {HistType::kTH2F, {axisEta, axisTPCNCls}});
    Selections.add("Selections/Electron/Rapidity/Leading/hTPCChi2NClsVsY", "hTPCChi2NClsVsY", {HistType::kTH2F, {axisEta, axisChi2}});
    Selections.add("Selections/Electron/Rapidity/Leading/hTPCNClsCrossedRowsVsY", "hTPCNClsCrossedRowsVsY", {HistType::kTH2F, {axisEta, axisTPCCrossed}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hITSNClsVsY", "hITSNClsVsY", {HistType::kTH2F, {axisEta, axisITSNCls}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hITSChi2NClsVsY", "hITSChi2NClsVsY", {HistType::kTH2F, {axisEta, axisChi2}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hTPCNClsVsY", "hTPCNClsVsY", {HistType::kTH2F, {axisEta, axisTPCNCls}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hTPCChi2NClsVsY", "hTPCChi2NClsVsY", {HistType::kTH2F, {axisEta, axisChi2}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hTPCNClsCrossedRowsVsY", "hTPCNClsCrossedRowsVsY", {HistType::kTH2F, {axisEta, axisTPCCrossed}});

    Selections.add("Selections/Electron/Pt/Leading/hITSNClsVsPt", "hITSNClsVsPt", {HistType::kTH2F, {axisPt, axisITSNCls}});
    Selections.add("Selections/Electron/Pt/Leading/hITSChi2NClsVsPt", "hITSChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisChi2}});
    Selections.add("Selections/Electron/Pt/Leading/hTPCNClsVsPt", "hTPCNClsVsPt", {HistType::kTH2F, {axisPt, axisTPCNCls}});
    Selections.add("Selections/Electron/Pt/Leading/hTPCChi2NClsVsPt", "hTPCChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisChi2}});
    Selections.add("Selections/Electron/Pt/Leading/hTPCNClsCrossedRowsVsPt", "hTPCNClsCrossedRowsVsPt", {HistType::kTH2F, {axisPt, axisTPCCrossed}});
    Selections.add("Selections/Electron/Pt/Subleading/hITSNClsVsPt", "hITSNClsVsPt", {HistType::kTH2F, {axisPt, axisITSNCls}});
    Selections.add("Selections/Electron/Pt/Subleading/hITSChi2NClsVsPt", "hITSChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisChi2}});
    Selections.add("Selections/Electron/Pt/Subleading/hTPCNClsVsPt", "hTPCNClsVsPt", {HistType::kTH2F, {axisPt, axisTPCNCls}});
    Selections.add("Selections/Electron/Pt/Subleading/hTPCChi2NClsVsPt", "hTPCChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisChi2}});
    Selections.add("Selections/Electron/Pt/Subleading/hTPCNClsCrossedRowsVsPt", "hTPCNClsCrossedRowsVsPt", {HistType::kTH2F, {axisPt, axisTPCCrossed}});

    Selections.add("Selections/Muon/Mass/Leading/hITSNClsVsM", "hITSNClsVsM", {HistType::kTH2F, {axisIVMSel, axisITSNCls}});
    Selections.add("Selections/Muon/Mass/Leading/hITSChi2NClsVsM", "hITSChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Muon/Mass/Leading/hTPCNClsVsM", "hTPCNClsVsM", {HistType::kTH2F, {axisIVMSel, axisTPCNCls}});
    Selections.add("Selections/Muon/Mass/Leading/hTPCChi2NClsVsM", "hTPCChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Muon/Mass/Leading/hTPCNClsCrossedRowsVsM", "hTPCNClsCrossedRowsVsM", {HistType::kTH2F, {axisIVMSel, axisTPCCrossed}});
    Selections.add("Selections/Muon/Mass/Subleading/hITSNClsVsM", "hITSNClsVsM", {HistType::kTH2F, {axisIVMSel, axisITSNCls}});
    Selections.add("Selections/Muon/Mass/Subleading/hITSChi2NClsVsM", "hITSChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Muon/Mass/Subleading/hTPCNClsVsM", "hTPCNClsVsM", {HistType::kTH2F, {axisIVMSel, axisTPCNCls}});
    Selections.add("Selections/Muon/Mass/Subleading/hTPCChi2NClsVsM", "hTPCChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisChi2}});
    Selections.add("Selections/Muon/Mass/Subleading/hTPCNClsCrossedRowsVsM", "hTPCNClsCrossedRowsVsM", {HistType::kTH2F, {axisIVMSel, axisTPCCrossed}});

    Selections.add("Selections/Muon/Rapidity/Leading/hITSNClsVsY", "hITSNClsVsY", {HistType::kTH2F, {axisEta, axisITSNCls}});
    Selections.add("Selections/Muon/Rapidity/Leading/hITSChi2NClsVsY", "hITSChi2NClsVsY", {HistType::kTH2F, {axisEta, axisChi2}});
    Selections.add("Selections/Muon/Rapidity/Leading/hTPCNClsVsY", "hTPCNClsVsY", {HistType::kTH2F, {axisEta, axisTPCNCls}});
    Selections.add("Selections/Muon/Rapidity/Leading/hTPCChi2NClsVsY", "hTPCChi2NClsVsY", {HistType::kTH2F, {axisEta, axisChi2}});
    Selections.add("Selections/Muon/Rapidity/Leading/hTPCNClsCrossedRowsVsY", "hTPCNClsCrossedRowsVsY", {HistType::kTH2F, {axisEta, axisTPCCrossed}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hITSNClsVsY", "hITSNClsVsY", {HistType::kTH2F, {axisEta, axisITSNCls}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hITSChi2NClsVsY", "hITSChi2NClsVsY", {HistType::kTH2F, {axisEta, axisChi2}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hTPCNClsVsY", "hTPCNClsVsY", {HistType::kTH2F, {axisEta, axisTPCNCls}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hTPCChi2NClsVsY", "hTPCChi2NClsVsY", {HistType::kTH2F, {axisEta, axisChi2}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hTPCNClsCrossedRowsVsY", "hTPCNClsCrossedRowsVsY", {HistType::kTH2F, {axisEta, axisTPCCrossed}});

    Selections.add("Selections/Muon/Pt/Leading/hITSNClsVsPt", "hITSNClsVsPt", {HistType::kTH2F, {axisPt, axisITSNCls}});
    Selections.add("Selections/Muon/Pt/Leading/hITSChi2NClsVsPt", "hITSChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisChi2}});
    Selections.add("Selections/Muon/Pt/Leading/hTPCNClsVsPt", "hTPCNClsVsPt", {HistType::kTH2F, {axisPt, axisTPCNCls}});
    Selections.add("Selections/Muon/Pt/Leading/hTPCChi2NClsVsPt", "hTPCChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisChi2}});
    Selections.add("Selections/Muon/Pt/Leading/hTPCNClsCrossedRowsVsPt", "hTPCNClsCrossedRowsVsPt", {HistType::kTH2F, {axisPt, axisTPCCrossed}});
    Selections.add("Selections/Muon/Pt/Subleading/hITSNClsVsPt", "hITSNClsVsPt", {HistType::kTH2F, {axisPt, axisITSNCls}});
    Selections.add("Selections/Muon/Pt/Subleading/hITSChi2NClsVsPt", "hITSChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisChi2}});
    Selections.add("Selections/Muon/Pt/Subleading/hTPCNClsVsPt", "hTPCNClsVsPt", {HistType::kTH2F, {axisPt, axisTPCNCls}});
    Selections.add("Selections/Muon/Pt/Subleading/hTPCChi2NClsVsPt", "hTPCChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisChi2}});
    Selections.add("Selections/Muon/Pt/Subleading/hTPCNClsCrossedRowsVsPt", "hTPCNClsCrossedRowsVsPt", {HistType::kTH2F, {axisPt, axisTPCCrossed}});

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

    // TGmu histograms
    TGmu.add("TGmu/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGmu.add("TGmu/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGmu.add("TGmu/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGmu.add("TGmu/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGmu.add("TGmu/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGmu.add("TGmu/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGmu.add("TGmu/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGmu.add("TGmu/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGmu.add("TGmu/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGmu.add("TGmu/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGmu.add("TGmu/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGmu.add("TGmu/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // TGmuCand histograms
    TGmuCand.add("TGmuCand/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGmuCand.add("TGmuCand/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGmuCand.add("TGmuCand/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGmuCand.add("TGmuCand/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGmuCand.add("TGmuCand/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGmuCand.add("TGmuCand/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGmuCand.add("TGmuCand/hPairPt", "hPairPt", {HistType::kTH1F, {axisPt}});
    TGmuCand.add("TGmuCand/hPairIVM", "hPairIVM", {HistType::kTH1F, {axisIVMWide}});
    TGmuCand.add("TGmuCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axisPt}});
    TGmuCand.add("TGmuCand/hJpsiPt2", "hJpsiPt2", {HistType::kTH1F, {axisPt2}});
    TGmuCand.add("TGmuCand/XnXn/hJpsiPt2", "hJpsiPt2", {HistType::kTH1F, {axisPt2}});
    TGmuCand.add("TGmuCand/OnOn/hJpsiPt2", "hJpsiPt2", {HistType::kTH1F, {axisPt2}});
    TGmuCand.add("TGmuCand/OnXn/hJpsiPt2", "hJpsiPt2", {HistType::kTH1F, {axisPt2}});
    TGmuCand.add("TGmuCand/XnOn/hJpsiPt2", "hJpsiPt2", {HistType::kTH1F, {axisPt2}});
    TGmuCand.add("TGmuCand/hJpsiPt2wide", "hJpsiPt2wide", {HistType::kTH1F, {axisPt2}});
    TGmuCand.add("TGmuCand/XnXn/hJpsiPt2wide", "hJpsiPt2wide", {HistType::kTH1F, {axisPt2wide}});
    TGmuCand.add("TGmuCand/OnOn/hJpsiPt2wide", "hJpsiPt2wide", {HistType::kTH1F, {axisPt2wide}});
    TGmuCand.add("TGmuCand/OnXn/hJpsiPt2wide", "hJpsiPt2wide", {HistType::kTH1F, {axisPt2wide}});
    TGmuCand.add("TGmuCand/XnOn/hJpsiPt2wide", "hJpsiPt2wide", {HistType::kTH1F, {axisPt2wide}});
    TGmuCand.add("TGmuCand/hJpsiRap", "hJpsiRap", {HistType::kTH1F, {axisEta}});
    TGmuCand.add("TGmuCand/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGmuCand.add("TGmuCand/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGmuCand.add("TGmuCand/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGmuCand.add("TGmuCand/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGmuCand.add("TGmuCand/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGmuCand.add("TGmuCand/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // TGel histograms
    TGel.add("TGel/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGel.add("TGel/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGel.add("TGel/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGel.add("TGel/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGel.add("TGel/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGel.add("TGel/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGel.add("TGel/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGel.add("TGel/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGel.add("TGel/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGel.add("TGel/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGel.add("TGel/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGel.add("TGel/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // TGelCand histograms
    TGelCand.add("TGelCand/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGelCand.add("TGelCand/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGelCand.add("TGelCand/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGelCand.add("TGelCand/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGelCand.add("TGelCand/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGelCand.add("TGelCand/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGelCand.add("TGelCand/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGelCand.add("TGelCand/hPairPt", "hPairPt", {HistType::kTH1F, {axisPt}});
    TGelCand.add("TGelCand/hPairIVM", "hPairIVM", {HistType::kTH1F, {axisIVMWide}});
    TGelCand.add("TGelCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axisPt}});
    TGelCand.add("TGelCand/hJpsiRap", "hJpsiRap", {HistType::kTH1F, {axisEta}});
    TGelCand.add("TGelCand/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGelCand.add("TGelCand/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGelCand.add("TGelCand/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGelCand.add("TGelCand/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGelCand.add("TGelCand/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // TGp histograms
    TGp.add("TGp/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGp.add("TGp/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGp.add("TGp/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGp.add("TGp/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGp.add("TGp/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGp.add("TGp/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGp.add("TGp/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGp.add("TGp/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGp.add("TGp/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGp.add("TGp/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGp.add("TGp/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // TGpCand histograms
    TGpCand.add("TGpCand/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGpCand.add("TGpCand/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGpCand.add("TGpCand/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGpCand.add("TGpCand/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGpCand.add("TGpCand/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGpCand.add("TGpCand/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGpCand.add("TGpCand/hPairPt", "hPairPt", {HistType::kTH1F, {axisPt}});
    TGpCand.add("TGpCand/hPairIVM", "hPairIVM", {HistType::kTH1F, {axisIVMWide}});
    TGpCand.add("TGpCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axisPt}});
    TGpCand.add("TGpCand/hJpsiRap", "hJpsiRap", {HistType::kTH1F, {axisEta}});
    TGpCand.add("TGpCand/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGpCand.add("TGpCand/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGpCand.add("TGpCand/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGpCand.add("TGpCand/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGpCand.add("TGpCand/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // JPsiToEl histograms
    JPsiToEl.add("JPsiToEl/Coherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Coherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/XnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/OnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/OnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Coherent/XnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    JPsiToEl.add("JPsiToEl/Incoherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/NotCoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Incoherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/NotCoherent/XnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/NotCoherent/OnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/NotCoherent/OnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/OnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/NotCoherent/XnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // JPsiToMu histograms
    JPsiToMu.add("JPsiToMu/Coherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Coherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/XnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/OnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/OnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Coherent/XnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    JPsiToMu.add("JPsiToMu/Incoherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/NotCoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Incoherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/NotCoherent/XnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/NotCoherent/OnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/NotCoherent/OnXn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/OnXn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/NotCoherent/XnOn/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // JPsiToP histograms
    JPsiToP.add("JPsiToP/Coherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToP.add("JPsiToP/Coherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToP.add("JPsiToP/Coherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToP.add("JPsiToP/Coherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToP.add("JPsiToP/Coherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToP.add("JPsiToP/Coherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Coherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Coherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToP.add("JPsiToP/Coherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToP.add("JPsiToP/Coherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToP.add("JPsiToP/Coherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToP.add("JPsiToP/Coherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    JPsiToP.add("JPsiToP/Incoherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToP.add("JPsiToP/Incoherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToP.add("JPsiToP/Incoherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToP.add("JPsiToP/Incoherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToP.add("JPsiToP/Incoherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToP.add("JPsiToP/Incoherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Incoherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Incoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToP.add("JPsiToP/Incoherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToP.add("JPsiToP/Incoherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToP.add("JPsiToP/Incoherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});

    // Correlation histograms
    Correlation.add("Correlation/Muon/Coherent/AccoplAngle", "AccoplAngle", {HistType::kTH1F, {axisAccAngle}});
    Correlation.add("Correlation/Muon/Coherent/CosTheta", "CosTheta", {HistType::kTH1F, {axisAngTheta}});
    Correlation.add("Correlation/Muon/Coherent/Phi", "Phi", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Coherent/Phi1", "Phi1", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Coherent/Phi2", "Phi2", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Coherent/CosThetaPhi", "CosThetaPhi", {HistType::kTH2F, {{axisAngTheta}, {axisPhi}}});

    Correlation.add("Correlation/Muon/Incoherent/AccoplAngle", "AccoplAngle", {HistType::kTH1F, {axisAccAngle}});
    Correlation.add("Correlation/Muon/Incoherent/CosTheta", "CosTheta", {HistType::kTH1F, {axisAngTheta}});
    Correlation.add("Correlation/Muon/Incoherent/Phi", "Phi", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Incoherent/Phi1", "Phi1", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Incoherent/Phi2", "Phi2", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Muon/Incoherent/CosThetaPhi", "CosThetaPhi", {HistType::kTH2F, {{axisAngTheta}, {axisPhi}}});

    Correlation.add("Correlation/Electron/Coherent/AccoplAngle", "AccoplAngle", {HistType::kTH1F, {axisAccAngle}});
    Correlation.add("Correlation/Electron/Coherent/CosTheta", "CosTheta", {HistType::kTH1F, {axisAngTheta}});
    Correlation.add("Correlation/Electron/Coherent/Phi", "Phi", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Electron/Coherent/Phi1", "Phi1", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Electron/Coherent/Phi2", "Phi2", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Electron/Coherent/CosThetaPhi", "CosThetaPhi", {HistType::kTH2F, {{axisAngTheta}, {axisPhi}}});

    Correlation.add("Correlation/Electron/Incoherent/AccoplAngle", "AccoplAngle", {HistType::kTH1F, {axisAccAngle}});
    Correlation.add("Correlation/Electron/Incoherent/CosTheta", "CosTheta", {HistType::kTH1F, {axisAngTheta}});
    Correlation.add("Correlation/Electron/Incoherent/Phi", "Phi", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Electron/Incoherent/Phi1", "Phi1", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Electron/Incoherent/Phi2", "Phi2", {HistType::kTH1F, {axisPhi}});
    Correlation.add("Correlation/Electron/Incoherent/CosThetaPhi", "CosThetaPhi", {HistType::kTH2F, {{axisAngTheta}, {axisPhi}}});

    // Asymmetry histograms
    Asymmetry.add("Asymmetry/Muon/Coherent/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/XnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/OnOn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/XnOn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/OnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/XnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/OnOn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/XnOn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/OnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/XnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/OnOn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/XnOn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/OnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/XnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/OnOn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/XnOn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/OnXn/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/XnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/OnOn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/XnOn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Coherent/OnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/XnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/OnOn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/XnOn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Muon/Incoherent/OnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/XnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/OnOn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/XnOn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/OnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/XnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/OnOn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/XnOn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/OnXn/DeltaPhiRandom", "DeltaPhiRandom", {HistType::kTH1F, {{180, -PI, PI}}});

    // MC histograms
    MC.add("MC/hNumberOfMCCollisions", "hNumberOfCollisions", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/hNumberOfRecoCollisions", "hNumberOfRecoCollisions", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/Muon/hNumberOfMCTracks", "hNumberOfMCTracks", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/Electron/hNumberOfMCTracks", "hNumberOfMCTracks", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/Proton/hNumberOfMCTracks", "hNumberOfMCTracks", {HistType::kTH1F, {{10, 0, 10}}});
    MC.add("MC/hPosZ", "hPosZ", {HistType::kTH1F, {{60, -15, 15}}});
    MC.add("MC/hPDG", "hPDG", {HistType::kTH1F, {{200000, -100000, 100000}}});
    MC.add("MC/Electron/hEta1", "hEta1", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Electron/hEta2", "hEta2", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Electron/hPhi1", "hPhi1", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Electron/hPhi2", "hPhi2", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Electron/hIVM", "hIVM", {HistType::kTH1F, {axisIVM}});
    MC.add("MC/Electron/hRapidity", "hRapidity", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Electron/hPt1", "hPt1", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Electron/hPt2", "hPt2", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Electron/hPt", "hPt", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Muon/hEta1", "hEta1", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Muon/hEta2", "hEta2", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Muon/hPhi1", "hPhi1", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Muon/hPhi2", "hPhi2", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Muon/hIVM", "hIVM", {HistType::kTH1F, {axisIVM}});
    MC.add("MC/Muon/hRapidity", "hRapidity", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Muon/hPt1", "hPt1", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Muon/hPt2", "hPt2", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Muon/hPt", "hPt", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Proton/hEta1", "hEta1", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Proton/hEta2", "hEta2", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Proton/hPhi1", "hPhi1", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Proton/hPhi2", "hPhi2", {HistType::kTH1F, {axisPhi}});
    MC.add("MC/Proton/hIVM", "hIVM", {HistType::kTH1F, {axisIVM}});
    MC.add("MC/Proton/hRapidity", "hRapidity", {HistType::kTH1F, {axisEta}});
    MC.add("MC/Proton/hPt1", "hPt1", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Proton/hPt2", "hPt2", {HistType::kTH1F, {axisPt}});
    MC.add("MC/Proton/hPt", "hPt", {HistType::kTH1F, {axisPt}});
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

  // template <typename C>
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

    // distinguish ZDC classes
    bool XnXn = false, OnOn = false, XnOn = false, OnXn = false;
    if (collision.energyCommonZNA() < ZNenergyCut && collision.energyCommonZNC() < ZNenergyCut) {
      OnOn = true;
    }
    if (collision.energyCommonZNA() > ZNenergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() > ZNenergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) {
      XnXn = true;
    }
    if (collision.energyCommonZNA() > ZNenergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() < ZNenergyCut) {
      XnOn = true;
    }
    if (collision.energyCommonZNA() < ZNenergyCut && collision.energyCommonZNC() > ZNenergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) {
      OnXn = true;
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
    if (std::abs(collision.posZ()) > cutVertexZ) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterCollisions"))->Fill(1);
      return;
    }

    RawData.get<TH1>(HIST("RawData/hPositionX"))->Fill(collision.posX());
    RawData.get<TH1>(HIST("RawData/hPositionY"))->Fill(collision.posY());
    RawData.get<TH1>(HIST("RawData/hPositionZ"))->Fill(collision.posZ());
    RawData.get<TH2>(HIST("RawData/hPositionXY"))->Fill(collision.posX(), collision.posY());
    RawData.get<TH1>(HIST("RawData/hZNACommonEnergy"))->Fill(collision.energyCommonZNA());
    RawData.get<TH1>(HIST("RawData/hZNCCommonEnergy"))->Fill(collision.energyCommonZNC());
    RawData.get<TH2>(HIST("RawData/hZNAvsZNCCommonEnergy"))->Fill(collision.energyCommonZNA(), collision.energyCommonZNC());
    RawData.get<TH1>(HIST("RawData/hZNCTime"))->Fill(collision.timeZNC());
    RawData.get<TH1>(HIST("RawData/hZNATime"))->Fill(collision.timeZNA());
    RawData.get<TH2>(HIST("RawData/hZNAvsZNCTime"))->Fill(collision.timeZNA(), collision.timeZNC());

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

          TGel.get<TH1>(HIST("TGel/hTrackPt1"))->Fill(trkDaughter1.pt());
          TGel.get<TH1>(HIST("TGel/hTrackPt2"))->Fill(trkDaughter2.pt());
          TGel.get<TH1>(HIST("TGel/hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
          TGel.get<TH1>(HIST("TGel/hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
          TGel.get<TH1>(HIST("TGel/hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
          TGel.get<TH1>(HIST("TGel/hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));

          // selections
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Leading/hITSNClsVsM"))->Fill(massJpsi, leadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Leading/hITSChi2NClsVsM"))->Fill(massJpsi, leadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Leading/hTPCNClsVsM"))->Fill(massJpsi, leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Leading/hTPCChi2NClsVsM"))->Fill(massJpsi, leadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Leading/hTPCNClsCrossedRowsVsM"))->Fill(massJpsi, subleadingP.tpcNClsCrossedRows());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Subleading/hITSNClsVsM"))->Fill(massJpsi, subleadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Subleading/hITSChi2NClsVsM"))->Fill(massJpsi, subleadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Subleading/hTPCNClsVsM"))->Fill(massJpsi, subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Subleading/hTPCChi2NClsVsM"))->Fill(massJpsi, subleadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Mass/Subleading/hTPCNClsCrossedRowsVsM"))->Fill(massJpsi, subleadingP.tpcNClsCrossedRows());

          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Leading/hITSNClsVsY"))->Fill(rapJpsi, leadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Leading/hITSChi2NClsVsY"))->Fill(rapJpsi, leadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Leading/hTPCNClsVsY"))->Fill(rapJpsi, leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Leading/hTPCChi2NClsVsY"))->Fill(rapJpsi, leadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Leading/hTPCNClsCrossedRowsVsY"))->Fill(rapJpsi, subleadingP.tpcNClsCrossedRows());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Subleading/hITSNClsVsY"))->Fill(rapJpsi, subleadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Subleading/hITSChi2NClsVsY"))->Fill(rapJpsi, subleadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Subleading/hTPCNClsVsY"))->Fill(rapJpsi, subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Subleading/hTPCChi2NClsVsY"))->Fill(rapJpsi, subleadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Rapidity/Subleading/hTPCNClsCrossedRowsVsY"))->Fill(rapJpsi, subleadingP.tpcNClsCrossedRows());

          Selections.get<TH2>(HIST("Selections/Electron/Pt/Leading/hITSNClsVsPt"))->Fill(RecoDecay::pt(mother), leadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Leading/hITSChi2NClsVsPt"))->Fill(RecoDecay::pt(mother), leadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Leading/hTPCNClsVsPt"))->Fill(RecoDecay::pt(mother), leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Leading/hTPCChi2NClsVsPt"))->Fill(RecoDecay::pt(mother), leadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Leading/hTPCNClsCrossedRowsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsCrossedRows());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Subleading/hITSNClsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Subleading/hITSChi2NClsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Subleading/hTPCNClsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Subleading/hTPCChi2NClsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Electron/Pt/Subleading/hTPCNClsCrossedRowsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsCrossedRows());

          if (trkDaughter1.hasTPC()) {
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          if (trkDaughter1.hasTOF()) {
            TGel.get<TH2>(HIST("TGel/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
          }
          if (trkDaughter2.hasTOF()) {
            TGel.get<TH2>(HIST("TGel/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
          }

          if (CandidateCuts(massJpsi, rapJpsi) != 1) {
            return;
          }

          TGelCand.get<TH1>(HIST("TGelCand/hTrackPt1"))->Fill(trkDaughter1.pt());
          TGelCand.get<TH1>(HIST("TGelCand/hTrackPt2"))->Fill(trkDaughter2.pt());
          TGelCand.get<TH1>(HIST("TGelCand/hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
          TGelCand.get<TH1>(HIST("TGelCand/hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
          TGelCand.get<TH1>(HIST("TGelCand/hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
          TGelCand.get<TH1>(HIST("TGelCand/hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));
          TGelCand.get<TH1>(HIST("TGelCand/hPairPt"))->Fill(RecoDecay::pt(mother));
          TGelCand.get<TH1>(HIST("TGelCand/hPairIVM"))->Fill(massJpsi);

          if (trkDaughter1.hasTPC()) {
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
            if (trkDaughter1.sign() < 0) {
              TGelCand.get<TH2>(HIST("TGelCand/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
            } else {
              TGelCand.get<TH2>(HIST("TGelCand/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
            }
          }
          if (trkDaughter2.hasTPC()) {
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          if (trkDaughter1.hasTOF()) {
            TGelCand.get<TH2>(HIST("TGelCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
          }
          if (trkDaughter2.hasTOF()) {
            TGelCand.get<TH2>(HIST("TGelCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
          }

          if (RecoDecay::pt(mother) < 0.2f) {
            JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hIVM"))->Fill(massJpsi);
            if (XnXn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hIVM"))->Fill(massJpsi);
            } else if (OnXn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hIVM"))->Fill(massJpsi);
            } else if (XnOn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hIVM"))->Fill(massJpsi);
            } else if (OnOn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hIVM"))->Fill(massJpsi);
            }
          } else if (RecoDecay::pt(mother) > 0.4f) {
            JPsiToEl.get<TH1>(HIST("JPsiToEl/NotCoherent/hIVM"))->Fill(massJpsi);
            if (XnXn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/NotCoherent/XnXn/hIVM"))->Fill(massJpsi);
            } else if (OnXn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/NotCoherent/OnXn/hIVM"))->Fill(massJpsi);
            } else if (XnOn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/NotCoherent/XnOn/hIVM"))->Fill(massJpsi);
            } else if (OnOn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/NotCoherent/OnOn/hIVM"))->Fill(massJpsi);
            }
          }
          if (RecoDecay::pt(mother) > 0.2f) {
            JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hIVM"))->Fill(massJpsi);
            if (XnXn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hIVM"))->Fill(massJpsi);
            } else if (OnXn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hIVM"))->Fill(massJpsi);
            } else if (XnOn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hIVM"))->Fill(massJpsi);
            } else if (OnOn) {
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hIVM"))->Fill(massJpsi);
            }
          }

          if ((massJpsi < maxJpsiMass) && (massJpsi > minJpsiMass)) {
            TGelCand.get<TH1>(HIST("TGelCand/hJpsiPt"))->Fill(RecoDecay::pt(mother));
            TGelCand.get<TH1>(HIST("TGelCand/hJpsiRap"))->Fill(rapJpsi);

            if (trkDaughter1.sign() < 0) {
              TGelCand.get<TH2>(HIST("TGelCand/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
            } else {
              TGelCand.get<TH2>(HIST("TGelCand/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
            }

            if (RecoDecay::pt(mother) < 0.2f) {
              // fill track histos
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPt1"))->Fill(trkDaughter1.pt());
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPt2"))->Fill(trkDaughter2.pt());
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));

              if (XnXn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (OnOn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (XnOn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (OnXn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              }

              if (trkDaughter1.hasTPC()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
              }
              if (trkDaughter2.hasTPC()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
              }

              if (trkDaughter1.hasTOF()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
              }
              // fill J/psi histos
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPt"))->Fill(RecoDecay::pt(mother));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hEta"))->Fill(RecoDecay::eta(mother));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPhi"))->Fill(RecoDecay::phi(mother));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hRap"))->Fill(rapJpsi);

              if (XnXn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnXn/hRap"))->Fill(rapJpsi);
              } else if (OnOn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnOn/hRap"))->Fill(rapJpsi);
              } else if (OnXn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/OnXn/hRap"))->Fill(rapJpsi);
              } else if (XnOn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/XnOn/hRap"))->Fill(rapJpsi);
              }

              float* q = correlation(&daughter[0], &daughter[1], &mom);
              Correlation.get<TH1>(HIST("Correlation/Electron/Coherent/Phi1"))->Fill(RecoDecay::phi(daughter1), 1.);
              Correlation.get<TH1>(HIST("Correlation/Electron/Coherent/Phi2"))->Fill(RecoDecay::phi(daughter2), 1.);
              Correlation.get<TH1>(HIST("Correlation/Electron/Coherent/Phi"))->Fill(q[1], 1.);
              Correlation.get<TH1>(HIST("Correlation/Electron/Coherent/CosTheta"))->Fill(q[2], 1.);
              Correlation.get<TH1>(HIST("Correlation/Electron/Coherent/AccoplAngle"))->Fill(q[0], 1.);
              Correlation.get<TH2>(HIST("Correlation/Electron/Coherent/CosThetaPhi"))->Fill(q[2], q[1]);

              double dp = DeltaPhi(daughter[0], daughter[1]);
              Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/DeltaPhi"))->Fill(dp);
              if (XnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/XnXn/DeltaPhi"))->Fill(dp);
              } else if (OnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/OnOn/DeltaPhi"))->Fill(dp);
              } else if (XnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/XnOn/DeltaPhi"))->Fill(dp);
              } else if (OnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/OnXn/DeltaPhi"))->Fill(dp);
              }

              double dpRandom = DeltaPhiRandom(daughter[0], daughter[1]);
              Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/DeltaPhiRandom"))->Fill(dpRandom);
              if (XnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/XnXn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (OnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/OnOn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (XnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/XnOn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (OnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Coherent/OnXn/DeltaPhiRandom"))->Fill(dpRandom);
              }

              delete[] q;
            } // end coherent electrons
            if (RecoDecay::pt(mother) > 0.2f) {
              // fill track histos
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPt1"))->Fill(trkDaughter1.pt());
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPt2"))->Fill(trkDaughter2.pt());
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));

              if (XnXn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (OnOn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (XnOn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (OnXn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              }

              if (trkDaughter1.hasTPC()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
              }
              if (trkDaughter2.hasTPC()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
              }

              if (trkDaughter1.hasTOF()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
              }
              // fill J/psi histos
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPt"))->Fill(RecoDecay::pt(mother));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hEta"))->Fill(RecoDecay::eta(mother));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPhi"))->Fill(RecoDecay::phi(mother));
              JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hRap"))->Fill(rapJpsi);

              if (XnXn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnXn/hRap"))->Fill(rapJpsi);
              } else if (OnOn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnOn/hRap"))->Fill(rapJpsi);
              } else if (OnXn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/OnXn/hRap"))->Fill(rapJpsi);
              } else if (XnOn) {
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/XnOn/hRap"))->Fill(rapJpsi);
              }

              float* q = correlation(&daughter[0], &daughter[1], &mom);
              Correlation.get<TH1>(HIST("Correlation/Electron/Incoherent/Phi1"))->Fill(RecoDecay::phi(daughter1), 1.);
              Correlation.get<TH1>(HIST("Correlation/Electron/Incoherent/Phi2"))->Fill(RecoDecay::phi(daughter2), 1.);
              Correlation.get<TH1>(HIST("Correlation/Electron/Incoherent/Phi"))->Fill(q[1], 1.);
              Correlation.get<TH1>(HIST("Correlation/Electron/Incoherent/CosTheta"))->Fill(q[2], 1.);
              Correlation.get<TH1>(HIST("Correlation/Electron/Incoherent/AccoplAngle"))->Fill(q[0], 1.);
              Correlation.get<TH2>(HIST("Correlation/Electron/Incoherent/CosThetaPhi"))->Fill(q[2], q[1]);

              double dp = DeltaPhi(daughter[0], daughter[1]);
              Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/DeltaPhi"))->Fill(dp);
              if (XnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/XnXn/DeltaPhi"))->Fill(dp);
              } else if (OnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/OnOn/DeltaPhi"))->Fill(dp);
              } else if (XnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/XnOn/DeltaPhi"))->Fill(dp);
              } else if (OnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/OnXn/DeltaPhi"))->Fill(dp);
              }

              double dpRandom = DeltaPhiRandom(daughter[0], daughter[1]);
              Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/DeltaPhiRandom"))->Fill(dpRandom);
              if (XnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/XnXn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (OnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/OnOn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (XnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/XnOn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (OnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Electron/Incoherent/OnXn/DeltaPhiRandom"))->Fill(dpRandom);
              }

              delete[] q;
            } // end incoherent electrons
          } // end mass cut
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

          TGmu.get<TH1>(HIST("TGmu/hTrackPt1"))->Fill(trkDaughter1.pt());
          TGmu.get<TH1>(HIST("TGmu/hTrackPt2"))->Fill(trkDaughter2.pt());
          TGmu.get<TH1>(HIST("TGmu/hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
          TGmu.get<TH1>(HIST("TGmu/hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
          TGmu.get<TH1>(HIST("TGmu/hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
          TGmu.get<TH1>(HIST("TGmu/hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));

          // selections
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Leading/hITSNClsVsM"))->Fill(massJpsi, leadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Leading/hITSChi2NClsVsM"))->Fill(massJpsi, leadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Leading/hTPCNClsVsM"))->Fill(massJpsi, leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Leading/hTPCChi2NClsVsM"))->Fill(massJpsi, leadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Leading/hTPCNClsCrossedRowsVsM"))->Fill(massJpsi, subleadingP.tpcNClsCrossedRows());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Subleading/hITSNClsVsM"))->Fill(massJpsi, subleadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Subleading/hITSChi2NClsVsM"))->Fill(massJpsi, subleadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Subleading/hTPCNClsVsM"))->Fill(massJpsi, subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Subleading/hTPCChi2NClsVsM"))->Fill(massJpsi, subleadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Mass/Subleading/hTPCNClsCrossedRowsVsM"))->Fill(massJpsi, subleadingP.tpcNClsCrossedRows());

          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Leading/hITSNClsVsY"))->Fill(rapJpsi, leadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Leading/hITSChi2NClsVsY"))->Fill(rapJpsi, leadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Leading/hTPCNClsVsY"))->Fill(rapJpsi, leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Leading/hTPCChi2NClsVsY"))->Fill(rapJpsi, leadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Leading/hTPCNClsCrossedRowsVsY"))->Fill(rapJpsi, subleadingP.tpcNClsCrossedRows());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Subleading/hITSNClsVsY"))->Fill(rapJpsi, subleadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Subleading/hITSChi2NClsVsY"))->Fill(rapJpsi, subleadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Subleading/hTPCNClsVsY"))->Fill(rapJpsi, subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Subleading/hTPCChi2NClsVsY"))->Fill(rapJpsi, subleadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Rapidity/Subleading/hTPCNClsCrossedRowsVsY"))->Fill(rapJpsi, subleadingP.tpcNClsCrossedRows());

          Selections.get<TH2>(HIST("Selections/Muon/Pt/Leading/hITSNClsVsPt"))->Fill(RecoDecay::pt(mother), leadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Leading/hITSChi2NClsVsPt"))->Fill(RecoDecay::pt(mother), leadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Leading/hTPCNClsVsPt"))->Fill(RecoDecay::pt(mother), leadingP.tpcNClsFindable() - leadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Leading/hTPCChi2NClsVsPt"))->Fill(RecoDecay::pt(mother), leadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Leading/hTPCNClsCrossedRowsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsCrossedRows());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Subleading/hITSNClsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.itsNCls());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Subleading/hITSChi2NClsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.itsChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Subleading/hTPCNClsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsFindable() - subleadingP.tpcNClsFindableMinusFound());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Subleading/hTPCChi2NClsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.tpcChi2NCl());
          Selections.get<TH2>(HIST("Selections/Muon/Pt/Subleading/hTPCNClsCrossedRowsVsPt"))->Fill(RecoDecay::pt(mother), subleadingP.tpcNClsCrossedRows());

          if (trkDaughter1.hasTPC()) {
            TGmu.get<TH2>(HIST("TGmu/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
            if (trkDaughter1.sign() < 0) {
              TGmu.get<TH2>(HIST("TGmu/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
            } else {
              TGmu.get<TH2>(HIST("TGmu/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
            }
          }
          if (trkDaughter2.hasTPC()) {
            TGmu.get<TH2>(HIST("TGmu/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          if (trkDaughter1.hasTOF()) {
            TGmu.get<TH2>(HIST("TGmu/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
          }
          if (trkDaughter2.hasTOF()) {

            TGmu.get<TH2>(HIST("TGmu/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
          }

          if (CandidateCuts(massJpsi, rapJpsi) != 1) {
            return;
          }

          TGmuCand.get<TH1>(HIST("TGmuCand/hTrackPt1"))->Fill(trkDaughter1.pt());
          TGmuCand.get<TH1>(HIST("TGmuCand/hTrackPt2"))->Fill(trkDaughter2.pt());
          TGmuCand.get<TH1>(HIST("TGmuCand/hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
          TGmuCand.get<TH1>(HIST("TGmuCand/hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
          TGmuCand.get<TH1>(HIST("TGmuCand/hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
          TGmuCand.get<TH1>(HIST("TGmuCand/hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));
          TGmuCand.get<TH1>(HIST("TGmuCand/hPairPt"))->Fill(RecoDecay::pt(mother));
          TGmuCand.get<TH1>(HIST("TGmuCand/hPairIVM"))->Fill(massJpsi);

          if (trkDaughter1.hasTPC()) {
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          if (trkDaughter1.hasTOF()) {
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
          }
          if (trkDaughter2.hasTOF()) {
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
          }

          if (RecoDecay::pt(mother) < 0.2f) {
            JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hIVM"))->Fill(massJpsi);
            if (XnXn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hIVM"))->Fill(massJpsi);
            } else if (OnXn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hIVM"))->Fill(massJpsi);
            } else if (XnOn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hIVM"))->Fill(massJpsi);
            } else if (OnOn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hIVM"))->Fill(massJpsi);
            }
          } else if (RecoDecay::pt(mother) > 0.4f) {
            JPsiToMu.get<TH1>(HIST("JPsiToMu/NotCoherent/hIVM"))->Fill(massJpsi);
            if (XnXn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/NotCoherent/XnXn/hIVM"))->Fill(massJpsi);
            } else if (OnXn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/NotCoherent/OnXn/hIVM"))->Fill(massJpsi);
            } else if (XnOn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/NotCoherent/XnOn/hIVM"))->Fill(massJpsi);
            } else if (OnOn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/NotCoherent/OnOn/hIVM"))->Fill(massJpsi);
            }
          }
          if (RecoDecay::pt(mother) > 0.2f) {
            JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hIVM"))->Fill(massJpsi);
            if (XnXn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hIVM"))->Fill(massJpsi);
            } else if (OnXn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hIVM"))->Fill(massJpsi);
            } else if (XnOn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hIVM"))->Fill(massJpsi);
            } else if (OnOn) {
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hIVM"))->Fill(massJpsi);
            }
          }

          if ((massJpsi < maxJpsiMass) && (massJpsi > minJpsiMass)) {
            TGmuCand.get<TH1>(HIST("TGmuCand/hJpsiPt"))->Fill(RecoDecay::pt(mother));
            if (RecoDecay::pt(mother) < 0.1f) {
              TGmuCand.get<TH1>(HIST("TGmuCand/hJpsiPt2"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              if (XnXn) {
                TGmuCand.get<TH1>(HIST("TGmuCand/XnXn/hJpsiPt2"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              } else if (XnOn) {
                TGmuCand.get<TH1>(HIST("TGmuCand/XnOn/hJpsiPt2"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              } else if (OnXn) {
                TGmuCand.get<TH1>(HIST("TGmuCand/OnXn/hJpsiPt2"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              } else if (OnOn) {
                TGmuCand.get<TH1>(HIST("TGmuCand/OnOn/hJpsiPt2"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              }
            }
            TGmuCand.get<TH1>(HIST("TGmuCand/hJpsiRap"))->Fill(rapJpsi);
            if (trkDaughter1.sign() < 0) {
              TGmuCand.get<TH2>(HIST("TGmuCand/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
            } else {
              TGmuCand.get<TH2>(HIST("TGmuCand/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
            }

            if (RecoDecay::pt(mother) < 0.2f) {
              // fill track histos
              TGmuCand.get<TH1>(HIST("TGmuCand/hJpsiPt2wide"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPt1"))->Fill(trkDaughter1.pt());
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPt2"))->Fill(trkDaughter2.pt());
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));

              if (XnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
                TGmuCand.get<TH1>(HIST("TGmuCand/XnXn/hJpsiPt2wide"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              } else if (OnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
                TGmuCand.get<TH1>(HIST("TGmuCand/OnOn/hJpsiPt2wide"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              } else if (XnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
                TGmuCand.get<TH1>(HIST("TGmuCand/XnOn/hJpsiPt2wide"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              } else if (OnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
                TGmuCand.get<TH1>(HIST("TGmuCand/OnXn/hJpsiPt2wide"))->Fill(RecoDecay::pt(mother) * RecoDecay::pt(mother));
              }

              if (trkDaughter1.hasTPC()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
              }
              if (trkDaughter2.hasTPC()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
              }
              if (trkDaughter1.hasTOF()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
              }
              // fill J/psi histos
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPt"))->Fill(RecoDecay::pt(mother));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hEta"))->Fill(RecoDecay::eta(mother));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPhi"))->Fill(RecoDecay::phi(mother));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hRap"))->Fill(rapJpsi);

              if (XnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnXn/hRap"))->Fill(rapJpsi);
              } else if (OnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hRap"))->Fill(rapJpsi);
              } else if (OnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hRap"))->Fill(rapJpsi);
              } else if (XnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hRap"))->Fill(rapJpsi);
              }

              float* q = correlation(&daughter[0], &daughter[1], &mom);
              Correlation.get<TH1>(HIST("Correlation/Muon/Coherent/Phi1"))->Fill(RecoDecay::phi(daughter1), 1.);
              Correlation.get<TH1>(HIST("Correlation/Muon/Coherent/Phi2"))->Fill(RecoDecay::phi(daughter2), 1.);
              Correlation.get<TH1>(HIST("Correlation/Muon/Coherent/Phi"))->Fill(q[1], 1.);
              Correlation.get<TH1>(HIST("Correlation/Muon/Coherent/CosTheta"))->Fill(q[2], 1.);
              Correlation.get<TH1>(HIST("Correlation/Muon/Coherent/AccoplAngle"))->Fill(q[0], 1.);
              Correlation.get<TH2>(HIST("Correlation/Muon/Coherent/CosThetaPhi"))->Fill(q[2], q[1]);

              double dp = DeltaPhi(daughter[0], daughter[1]);
              Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/DeltaPhi"))->Fill(dp);
              if (XnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/XnXn/DeltaPhi"))->Fill(dp);
              } else if (OnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/OnOn/DeltaPhi"))->Fill(dp);
              } else if (XnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/XnOn/DeltaPhi"))->Fill(dp);
              } else if (OnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/OnXn/DeltaPhi"))->Fill(dp);
              }
              double dpRandom = DeltaPhiRandom(daughter[0], daughter[1]);
              Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/DeltaPhiRandom"))->Fill(dpRandom);
              if (XnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/XnXn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (OnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/OnOn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (XnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/XnOn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (OnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Coherent/OnXn/DeltaPhiRandom"))->Fill(dpRandom);
              }

              delete[] q;
            }
            if (RecoDecay::pt(mother) > 0.2f) {
              // fill track histos
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hPt1"))->Fill(trkDaughter1.pt());
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hPt2"))->Fill(trkDaughter2.pt());
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));

              if (XnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (OnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (XnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (OnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              }

              if (trkDaughter1.hasTPC()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
              }
              if (trkDaughter2.hasTPC()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
              }
              if (trkDaughter1.hasTOF()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
              }

              // fill J/psi histos
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hPt"))->Fill(RecoDecay::pt(mother));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hEta"))->Fill(RecoDecay::eta(mother));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hPhi"))->Fill(RecoDecay::phi(mother));
              JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hRap"))->Fill(rapJpsi);

              if (XnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnXn/hRap"))->Fill(rapJpsi);
              } else if (OnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnOn/hRap"))->Fill(rapJpsi);
              } else if (OnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/OnXn/hRap"))->Fill(rapJpsi);
              } else if (XnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hPt"))->Fill(RecoDecay::pt(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hEta"))->Fill(RecoDecay::eta(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hPhi"))->Fill(RecoDecay::phi(mother));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/XnOn/hRap"))->Fill(rapJpsi);
              }

              float* q = correlation(&daughter[0], &daughter[1], &mom);
              Correlation.get<TH1>(HIST("Correlation/Muon/Incoherent/Phi1"))->Fill(RecoDecay::phi(daughter1), 1.);
              Correlation.get<TH1>(HIST("Correlation/Muon/Incoherent/Phi2"))->Fill(RecoDecay::phi(daughter2), 1.);
              Correlation.get<TH1>(HIST("Correlation/Muon/Incoherent/Phi"))->Fill(q[1], 1.);
              Correlation.get<TH1>(HIST("Correlation/Muon/Incoherent/CosTheta"))->Fill(q[2], 1.);
              Correlation.get<TH1>(HIST("Correlation/Muon/Incoherent/AccoplAngle"))->Fill(q[0], 1.);
              Correlation.get<TH2>(HIST("Correlation/Muon/Incoherent/CosThetaPhi"))->Fill(q[2], q[1]);

              double dp = DeltaPhi(daughter[0], daughter[1]);
              Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/DeltaPhi"))->Fill(dp);
              if (XnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/XnXn/DeltaPhi"))->Fill(dp);
              } else if (OnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/OnOn/DeltaPhi"))->Fill(dp);
              } else if (XnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/XnOn/DeltaPhi"))->Fill(dp);
              } else if (OnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/OnXn/DeltaPhi"))->Fill(dp);
              }

              double dpRandom = DeltaPhiRandom(daughter[0], daughter[1]);
              Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/DeltaPhiRandom"))->Fill(dpRandom);
              if (XnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/XnXn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (OnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/OnOn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (XnOn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/XnOn/DeltaPhiRandom"))->Fill(dpRandom);
              } else if (OnXn) {
                Asymmetry.get<TH1>(HIST("Asymmetry/Muon/Incoherent/OnXn/DeltaPhiRandom"))->Fill(dpRandom);
              }

              delete[] q;
            }
          }
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

          TGp.get<TH1>(HIST("TGp/hTrackPt1"))->Fill(trkDaughter1.pt());
          TGp.get<TH1>(HIST("TGp/hTrackPt2"))->Fill(trkDaughter2.pt());
          TGp.get<TH1>(HIST("TGp/hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
          TGp.get<TH1>(HIST("TGp/hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
          TGp.get<TH1>(HIST("TGp/hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
          TGp.get<TH1>(HIST("TGp/hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));

          if (trkDaughter1.hasTPC()) {
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          if (trkDaughter1.hasTOF()) {
            TGp.get<TH2>(HIST("TGp/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
          }
          if (trkDaughter2.hasTOF()) {
            TGp.get<TH2>(HIST("TGp/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
          }

          if (CandidateCuts(massJpsi, rapJpsi) != 1) {
            return;
          }

          TGpCand.get<TH1>(HIST("TGpCand/hTrackPt1"))->Fill(trkDaughter1.pt());
          TGpCand.get<TH1>(HIST("TGpCand/hTrackPt2"))->Fill(trkDaughter2.pt());
          TGpCand.get<TH1>(HIST("TGpCand/hTrackEta1"))->Fill(RecoDecay::eta(daughter1));
          TGpCand.get<TH1>(HIST("TGpCand/hTrackEta2"))->Fill(RecoDecay::eta(daughter2));
          TGpCand.get<TH1>(HIST("TGpCand/hTrackPhi1"))->Fill(RecoDecay::phi(daughter1));
          TGpCand.get<TH1>(HIST("TGpCand/hTrackPhi2"))->Fill(RecoDecay::phi(daughter2));
          TGpCand.get<TH1>(HIST("TGpCand/hPairPt"))->Fill(RecoDecay::pt(mother));
          TGpCand.get<TH1>(HIST("TGpCand/hPairIVM"))->Fill(massJpsi);

          if (trkDaughter1.hasTPC()) {
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          if (trkDaughter1.hasTOF()) {
            TGpCand.get<TH2>(HIST("TGpCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
          }
          if (trkDaughter2.hasTOF()) {
            TGpCand.get<TH2>(HIST("TGpCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
          }

          if (RecoDecay::pt(mother) < 0.2f) {
            JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hIVM"))->Fill(massJpsi);
          } else {
            JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hIVM"))->Fill(massJpsi);
          }

          if (massJpsi < maxJpsiMass && massJpsi > minJpsiMass) {
            TGpCand.get<TH1>(HIST("TGpCand/hJpsiPt"))->Fill(RecoDecay::pt(mother));
            TGpCand.get<TH1>(HIST("TGpCand/hJpsiRap"))->Fill(rapJpsi);
            if (RecoDecay::pt(mother) < 0.2f) {
              // fill track histos
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPt1"))->Fill(trkDaughter1.pt());
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPt2"))->Fill(trkDaughter2.pt());
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              if (trkDaughter1.hasTPC()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
              }
              if (trkDaughter2.hasTPC()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
              }
              if (trkDaughter1.hasTOF()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
              }
              // fill J/psi histos
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPt"))->Fill(RecoDecay::pt(mother));
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hEta"))->Fill(RecoDecay::eta(mother));
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPhi"))->Fill(RecoDecay::phi(mother));
              JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hRap"))->Fill(rapJpsi);
            } // end coherent
            if (RecoDecay::pt(mother) > 0.2f) {
              // fill track histos
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPt1"))->Fill(trkDaughter1.pt());
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPt2"))->Fill(trkDaughter2.pt());
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              if (trkDaughter1.hasTPC()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
              }
              if (trkDaughter2.hasTPC()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
              }
              if (trkDaughter1.hasTOF()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
              }
              // fill J/psi histos
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPt"))->Fill(RecoDecay::pt(mother));
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hEta"))->Fill(RecoDecay::eta(mother));
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPhi"))->Fill(RecoDecay::phi(mother));
              JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hRap"))->Fill(rapJpsi);
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
    std::array<float, 3> daughPart1El;
    std::array<float, 3> daughPart2El;
    std::array<float, 3> daughPart1Pr;
    std::array<float, 3> daughPart2Pr;

    // fill number of particles
    for (auto const& mcParticle : mcParticles) {
      MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(1.);
      MC.get<TH1>(HIST("MC/Electron/hNumberOfMCTracks"))->Fill(1.);
      MC.get<TH1>(HIST("MC/Proton/hNumberOfMCTracks"))->Fill(1.);
      MC.get<TH1>(HIST("MC/hPDG"))->Fill(mcParticle.pdgCode());

      // check if only muons and physical primaries are present
      if ((mcParticle.pdgCode() == 13) && mcParticle.isPhysicalPrimary()) {
        daughPart1Mu = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};

        MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(2.);
      }
      if ((mcParticle.pdgCode() == -13) && mcParticle.isPhysicalPrimary()) {
        daughPart2Mu = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};

        MC.get<TH1>(HIST("MC/Muon/hNumberOfMCTracks"))->Fill(3.);
      }
      if ((mcParticle.pdgCode() == 11) && mcParticle.isPhysicalPrimary()) {
        daughPart1El = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};

        MC.get<TH1>(HIST("MC/Electron/hNumberOfMCTracks"))->Fill(2.);
      }
      if ((mcParticle.pdgCode() == -11) && mcParticle.isPhysicalPrimary()) {
        daughPart2El = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};

        MC.get<TH1>(HIST("MC/Electron/hNumberOfMCTracks"))->Fill(3.);
      }
      if ((mcParticle.pdgCode() == 2212) && mcParticle.isPhysicalPrimary()) {
        daughPart1Pr = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};

        MC.get<TH1>(HIST("MC/Proton/hNumberOfMCTracks"))->Fill(2.);
      }
      if ((mcParticle.pdgCode() == -2212) && mcParticle.isPhysicalPrimary()) {
        daughPart2Pr = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};

        MC.get<TH1>(HIST("MC/Proton/hNumberOfMCTracks"))->Fill(3.);
      }
    }

    // calculate J/psi system and needed distributions
    std::array<double, 3> motherPartMu = {daughPart1Mu[0] + daughPart2Mu[0], daughPart1Mu[1] + daughPart2Mu[1], daughPart1Mu[2] + daughPart2Mu[2]};
    std::array<double, 3> motherPartEl = {daughPart1El[0] + daughPart2El[0], daughPart1El[1] + daughPart2El[1], daughPart1El[2] + daughPart2El[2]};
    std::array<double, 3> motherPartPr = {daughPart1Pr[0] + daughPart2Pr[0], daughPart1Pr[1] + daughPart2Pr[1], daughPart1Pr[2] + daughPart2Pr[2]};

    float massEl = o2::constants::physics::MassElectron;
    float massMu = o2::constants::physics::MassMuonMinus;
    float massPr = o2::constants::physics::MassProton;

    auto arrMomMu = std::array{daughPart1Mu, daughPart2Mu};
    float massJpsiMu = RecoDecay::m(arrMomMu, std::array{massMu, massMu});
    auto arrMomEl = std::array{daughPart1El, daughPart2El};
    float massJpsiEl = RecoDecay::m(arrMomEl, std::array{massEl, massEl});
    auto arrMomPr = std::array{daughPart1Pr, daughPart2Pr};
    float massJpsiPr = RecoDecay::m(arrMomPr, std::array{massPr, massPr});

    MC.get<TH1>(HIST("MC/Muon/hEta1"))->Fill(RecoDecay::eta(daughPart1Mu));
    MC.get<TH1>(HIST("MC/Muon/hEta2"))->Fill(RecoDecay::eta(daughPart2Mu));
    MC.get<TH1>(HIST("MC/Muon/hPhi1"))->Fill(RecoDecay::phi(daughPart1Mu));
    MC.get<TH1>(HIST("MC/Muon/hPhi2"))->Fill(RecoDecay::phi(daughPart2Mu));
    MC.get<TH1>(HIST("MC/Muon/hPt1"))->Fill(RecoDecay::pt(daughPart1Mu));
    MC.get<TH1>(HIST("MC/Muon/hPt2"))->Fill(RecoDecay::pt(daughPart2Mu));
    MC.get<TH1>(HIST("MC/Muon/hPt"))->Fill(RecoDecay::pt(motherPartMu));
    MC.get<TH1>(HIST("MC/Muon/hIVM"))->Fill(massJpsiMu);
    MC.get<TH1>(HIST("MC/Muon/hRapidity"))->Fill(RecoDecay::y(motherPartMu, massJpsiMu));

    MC.get<TH1>(HIST("MC/Electron/hEta1"))->Fill(RecoDecay::eta(daughPart1El));
    MC.get<TH1>(HIST("MC/Electron/hEta2"))->Fill(RecoDecay::eta(daughPart2El));
    MC.get<TH1>(HIST("MC/Electron/hPhi1"))->Fill(RecoDecay::phi(daughPart1El));
    MC.get<TH1>(HIST("MC/Electron/hPhi2"))->Fill(RecoDecay::phi(daughPart2El));
    MC.get<TH1>(HIST("MC/Electron/hPt1"))->Fill(RecoDecay::pt(daughPart1El));
    MC.get<TH1>(HIST("MC/Electron/hPt2"))->Fill(RecoDecay::pt(daughPart2El));
    MC.get<TH1>(HIST("MC/Electron/hPt"))->Fill(RecoDecay::pt(motherPartEl));
    MC.get<TH1>(HIST("MC/Electron/hIVM"))->Fill(massJpsiEl);
    MC.get<TH1>(HIST("MC/Electron/hRapidity"))->Fill(RecoDecay::y(motherPartEl, massJpsiEl));

    MC.get<TH1>(HIST("MC/Proton/hEta1"))->Fill(RecoDecay::eta(daughPart1Pr));
    MC.get<TH1>(HIST("MC/Proton/hEta2"))->Fill(RecoDecay::eta(daughPart2Pr));
    MC.get<TH1>(HIST("MC/Proton/hPhi1"))->Fill(RecoDecay::phi(daughPart1Pr));
    MC.get<TH1>(HIST("MC/Proton/hPhi2"))->Fill(RecoDecay::phi(daughPart2Pr));
    MC.get<TH1>(HIST("MC/Proton/hPt1"))->Fill(RecoDecay::pt(daughPart1Pr));
    MC.get<TH1>(HIST("MC/Proton/hPt2"))->Fill(RecoDecay::pt(daughPart2Pr));
    MC.get<TH1>(HIST("MC/Proton/hPt"))->Fill(RecoDecay::pt(motherPartPr));
    MC.get<TH1>(HIST("MC/Proton/hIVM"))->Fill(massJpsiPr);
    MC.get<TH1>(HIST("MC/Proton/hRapidity"))->Fill(RecoDecay::y(motherPartPr, massJpsiPr));

  } // end MC skimmed process

  template <typename C>
  void processMCU(C const& collisions)
  {
    MC.fill(HIST("MC/hNumberOfRecoCollisions"), collisions.size());
  } // end MC unskimmed process

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

  void processMCUnskimmed(aod::McCollision const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions /*, aod::McParticles const& mcParticles, aod::TracksIU const& tracks*/)
  {
    processMCU(collisions);
  }

  PROCESS_SWITCH(UpcJpsiCentralBarrel, processDGrecoLevel, "Iterate over DG skimmed data.", false);
  PROCESS_SWITCH(UpcJpsiCentralBarrel, processSGrecoLevel, "Iterate over SG skimmed data.", false);
  PROCESS_SWITCH(UpcJpsiCentralBarrel, processMCtruth, "Iterate of MC true data.", true);
  PROCESS_SWITCH(UpcJpsiCentralBarrel, processMCUnskimmed, "Iterate over unskimmed data.", true);
}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcJpsiCentralBarrel>(cfgc, TaskName{"upc-jpsi-corr"}),
  };
}
