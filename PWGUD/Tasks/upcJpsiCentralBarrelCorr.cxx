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
  ConfigurableAxis axisP{"axisP", {250.0f, 0.1f, 3.0f}, "#it{p} (GeV/#it{c})"};
  ConfigurableAxis axisEta{"axisEta", {250.0f, -1.5f, 1.5f}, "#eta (-)"};
  ConfigurableAxis axisCounter{"axisCounter", {20.0f, 0.0f, 20.0f}, "Number of events (-)"};
  ConfigurableAxis axisPhi{"axisPhi", {250.0f, 0, TwoPI}, "#phi (rad)"};
  ConfigurableAxis axisAccAngle{"axisAccAngle", {250.0f, -0.2f, 0.2f}, "accAngle"};
  ConfigurableAxis axisAngTheta{"axisAngTheta", {250.0f, -1.5f, 1.5f}, "cos #theta (-)"};
  ConfigurableAxis axisTPC{"axisTPC", {1000.0f, 0, 200.0f}, "TPC d#it{E}/d#it{x}"};
  ConfigurableAxis axisTOF{"axisTOF", {1000.0f, 0, 200.0f}, "TOF d#it{E}/d#it{x}"};
  ConfigurableAxis axisBetaTOF{"axisBetaTOF", {100.0f, 0, 1.5}, "TOF #beta"};
  ConfigurableAxis axisSigma{"axisSigma", {50, -25, 25}, "#sigma"};
  ConfigurableAxis axisZDCEnergy{"axisZDCEnergy", {250, -5.0, 20.0}, "ZDC energy"};
  ConfigurableAxis axisZDCTime{"axisZDCTime", {200, -10.0, 10.0}, "ZDC time"};
  ConfigurableAxis axisDCA{"axisDCA", {1000, -20.0, 20.0}, "DCA"};
  ConfigurableAxis axisChi2{"axisChi2", {1000, 0.0, 100.0}, "Chi2"};
  ConfigurableAxis axisIVMSel{"axisIVMSel", {1000, 0.0, 10.0}, "IVM"};
  ConfigurableAxis axisCounterSel{"axisCounterSel", {1000, 0.0, 200.0}, "Selection"};

  // configurable cuts (modify in json)
  // track quality cuts
  Configurable<int> TPCNClsCrossedRows{"TPCNClsCrossedRows", 70, "number of crossed rows in TPC"};
  Configurable<bool> TOFAtLeastOneProton{"TOFAtLeastOneProton", false, "at least one candidate track has TOF hits"};
  Configurable<bool> TOFBothProtons{"TOFBothProtons", false, "both candidate protons have TOF hits"};
  Configurable<bool> TOFOneProton{"TOFOneProton", false, "one candidate proton has TOF hits"};
  Configurable<bool> DCAcut{"DCAcut", false, "DCA cut from run2."};
  Configurable<bool> newCutTPC{"newCutTPC", false, "New cuts for TPC quality tracks."};
  Configurable<float> TPCNSigmaMu{"TPCNSigmaMu", 3, "PID for TPC Mu track"};
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
    Statistics.add("Statistics/hNumberGTselected", "hNumberGTselected", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTel", "hNumberGTel", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTmu", "hNumberGTmu", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTp", "hNumberGTp", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTelSigma", "hNumberGTelSigma", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTmuSigma", "hNumberGTmuSigma", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTpSigma", "hNumberGTpSigma", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTpSigmaTOF", "hNumberGTpSigmaTOF", {HistType::kTH1F, {axisCounter}});
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
    RawData.add("RawData/hTPCNCls", "hITSNCls", {HistType::kTH1F, {axisCounter}});
    RawData.add("RawData/hITSChi2NCls", "hITSChi2NCls", {HistType::kTH1F, {axisChi2}});
    RawData.add("RawData/hTPCChi2NCls", "hTPCChi2NCls", {HistType::kTH1F, {axisChi2}});
    RawData.add("RawData/hPositionZ", "hPositionZ", {HistType::kTH1F, {axisSigma}});
    RawData.add("RawData/hPositionX", "hPositionX", {HistType::kTH1F, {axisSigma}});
    RawData.add("RawData/hPositionY", "hPositionY", {HistType::kTH1F, {axisSigma}});
    RawData.add("RawData/hPositionXY", "hPositionXY", {HistType::kTH2F, {axisSigma, axisSigma}});
    RawData.add("RawData/hZNACommonEnergy", "hZNACommonEnergy", {HistType::kTH1F, {axisZDCEnergy}});
    RawData.add("RawData/hZNCCommonEnergy", "hZNCCommonEnergy", {HistType::kTH1F, {axisZDCEnergy}});
    RawData.add("RawData/hZNATime", "hZNATime", {HistType::kTH1F, {axisZDCTime}});
    RawData.add("RawData/hZNCTime", "hZNCTime", {HistType::kTH1F, {axisZDCTime}});
    RawData.add("RawData/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    RawData.add("RawData/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    RawData.add("RawData/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    RawData.add("RawData/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    RawData.add("RawData/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    RawData.add("RawData/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    RawData.add("RawData/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    RawData.add("RawData/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    RawData.add("RawData/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

    // Selection checks
    Selections.add("Selections/Electron/Mass/Leading/hITSNClsVsM", "hITSNClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Leading/hITSChi2NClsVsM", "hITSChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCNClsVsM", "hTPCNClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCChi2NClsVsM", "hTPCChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Leading/hTPCNClsCrossedRowsVsM", "hTPCNClsCrossedRowsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Subleading/hITSNClsVsM", "hITSNClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Subleading/hITSChi2NClsVsM", "hITSChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Subleading/hTPCNClsVsM", "hTPCNClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Subleading/hTPCChi2NClsVsM", "hTPCChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Electron/Mass/Subleading/hTPCNClsCrossedRowsVsM", "hTPCNClsCrossedRowsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});

    Selections.add("Selections/Electron/Rapidity/Leading/hITSNClsVsY", "hITSNClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Leading/hITSChi2NClsVsY", "hITSChi2NClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Leading/hTPCNClsVsY", "hTPCNClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Leading/hTPCChi2NClsVsY", "hTPCChi2NClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Leading/hTPCNClsCrossedRowsVsY", "hTPCNClsCrossedRowsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hITSNClsVsY", "hITSNClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hITSChi2NClsVsY", "hITSChi2NClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hTPCNClsVsY", "hTPCNClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hTPCChi2NClsVsY", "hTPCChi2NClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Electron/Rapidity/Subleading/hTPCNClsCrossedRowsVsY", "hTPCNClsCrossedRowsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});

    Selections.add("Selections/Electron/Pt/Leading/hITSNClsVsPt", "hITSNClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Leading/hITSChi2NClsVsPt", "hITSChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Leading/hTPCNClsVsPt", "hTPCNClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Leading/hTPCChi2NClsVsPt", "hTPCChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Leading/hTPCNClsCrossedRowsVsPt", "hTPCNClsCrossedRowsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Subleading/hITSNClsVsPt", "hITSNClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Subleading/hITSChi2NClsVsPt", "hITSChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Subleading/hTPCNClsVsPt", "hTPCNClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Subleading/hTPCChi2NClsVsPt", "hTPCChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Electron/Pt/Subleading/hTPCNClsCrossedRowsVsPt", "hTPCNClsCrossedRowsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});

    Selections.add("Selections/Muon/Mass/Leading/hITSNClsVsM", "hITSNClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Leading/hITSChi2NClsVsM", "hITSChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Leading/hTPCNClsVsM", "hTPCNClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Leading/hTPCChi2NClsVsM", "hTPCChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Leading/hTPCNClsCrossedRowsVsM", "hTPCNClsCrossedRowsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Subleading/hITSNClsVsM", "hITSNClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Subleading/hITSChi2NClsVsM", "hITSChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Subleading/hTPCNClsVsM", "hTPCNClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Subleading/hTPCChi2NClsVsM", "hTPCChi2NClsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});
    Selections.add("Selections/Muon/Mass/Subleading/hTPCNClsCrossedRowsVsM", "hTPCNClsCrossedRowsVsM", {HistType::kTH2F, {axisIVMSel, axisCounterSel}});

    Selections.add("Selections/Muon/Rapidity/Leading/hITSNClsVsY", "hITSNClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Leading/hITSChi2NClsVsY", "hITSChi2NClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Leading/hTPCNClsVsY", "hTPCNClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Leading/hTPCChi2NClsVsY", "hTPCChi2NClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Leading/hTPCNClsCrossedRowsVsY", "hTPCNClsCrossedRowsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hITSNClsVsY", "hITSNClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hITSChi2NClsVsY", "hITSChi2NClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hTPCNClsVsY", "hTPCNClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hTPCChi2NClsVsY", "hTPCChi2NClsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});
    Selections.add("Selections/Muon/Rapidity/Subleading/hTPCNClsCrossedRowsVsY", "hTPCNClsCrossedRowsVsY", {HistType::kTH2F, {axisEta, axisCounterSel}});

    Selections.add("Selections/Muon/Pt/Leading/hITSNClsVsPt", "hITSNClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Leading/hITSChi2NClsVsPt", "hITSChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Leading/hTPCNClsVsPt", "hTPCNClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Leading/hTPCChi2NClsVsPt", "hTPCChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Leading/hTPCNClsCrossedRowsVsPt", "hTPCNClsCrossedRowsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Subleading/hITSNClsVsPt", "hITSNClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Subleading/hITSChi2NClsVsPt", "hITSChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Subleading/hTPCNClsVsPt", "hTPCNClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Subleading/hTPCChi2NClsVsPt", "hTPCChi2NClsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});
    Selections.add("Selections/Muon/Pt/Subleading/hTPCNClsCrossedRowsVsPt", "hTPCNClsCrossedRowsVsPt", {HistType::kTH2F, {axisPt, axisCounterSel}});

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
    PVContributors.add("PVContributors/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    PVContributors.add("PVContributors/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    PVContributors.add("PVContributors/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    PVContributors.add("PVContributors/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    PVContributors.add("PVContributors/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

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
    TG.add("TG/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    TG.add("TG/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TG.add("TG/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    TG.add("TG/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    TG.add("TG/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});
    TG.add("TG/TPC/hNsigmaMu", "hNsigmaMu", HistType::kTH1F, {axisSigma});
    TG.add("TG/TPC/hNsigmaEl", "hNsigmaEl", HistType::kTH1F, {axisSigma});
    TG.add("TG/TPC/hNsigmaPr", "hNsigmaPr", HistType::kTH1F, {axisSigma});
    TG.add("TG/TPC/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TG.add("TG/TOF/hNsigmaMu", "hNsigmaMu", HistType::kTH1F, {axisSigma});
    TG.add("TG/TOF/hNsigmaEl", "hNsigmaEl", HistType::kTH1F, {axisSigma});
    TG.add("TG/TOF/hNsigmaPr", "hNsigmaPr", HistType::kTH1F, {axisSigma});

    // TGmu histograms
    TGmu.add("TGmu/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGmu.add("TGmu/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGmu.add("TGmu/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGmu.add("TGmu/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGmu.add("TGmu/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGmu.add("TGmu/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGmu.add("TGmu/hNsigmaMu", "hNsigmaMu", HistType::kTH1F, {axisSigma});
    TGmu.add("TGmu/hNsigmaMuTOF", "hNsigmaMuTOF", HistType::kTH1F, {axisSigma});
    TGmu.add("TGmu/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGmu.add("TGmu/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGmu.add("TGmu/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGmu.add("TGmu/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGmu.add("TGmu/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGmu.add("TGmu/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    TGmu.add("TGmu/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TGmu.add("TGmu/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    TGmu.add("TGmu/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    TGmu.add("TGmu/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

    // TGmuCand histograms
    TGmuCand.add("TGmuCand/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGmuCand.add("TGmuCand/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGmuCand.add("TGmuCand/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGmuCand.add("TGmuCand/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGmuCand.add("TGmuCand/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGmuCand.add("TGmuCand/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGmuCand.add("TGmuCand/hTrackITSNcls1", "hTrackITSNcls1", {HistType::kTH1F, {axisCounter}});
    TGmuCand.add("TGmuCand/hTrackITSNcls2", "hTrackITSNcls2", {HistType::kTH1F, {axisCounter}});
    TGmuCand.add("TGmuCand/hPairPt", "hPairPt", {HistType::kTH1F, {axisPt}});
    TGmuCand.add("TGmuCand/hPairIVM", "hPairIVM", {HistType::kTH1F, {axisIVMWide}});
    TGmuCand.add("TGmuCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axisPt}});
    TGmuCand.add("TGmuCand/hJpsiRap", "hJpsiRap", {HistType::kTH1F, {axisEta}});
    TGmuCand.add("TGmuCand/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGmuCand.add("TGmuCand/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGmuCand.add("TGmuCand/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGmuCand.add("TGmuCand/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGmuCand.add("TGmuCand/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGmuCand.add("TGmuCand/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    TGmuCand.add("TGmuCand/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TGmuCand.add("TGmuCand/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    TGmuCand.add("TGmuCand/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    TGmuCand.add("TGmuCand/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

    // TGel histograms
    TGel.add("TGel/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGel.add("TGel/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGel.add("TGel/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGel.add("TGel/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGel.add("TGel/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGel.add("TGel/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGel.add("TGel/hNsigmaEl", "hNsigmaEl", HistType::kTH1F, {axisSigma});
    TGel.add("TGel/hNsigmaElTOF", "hNsigmaElTOF", HistType::kTH1F, {axisSigma});
    TGel.add("TGel/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGel.add("TGel/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGel.add("TGel/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGel.add("TGel/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGel.add("TGel/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGel.add("TGel/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    TGel.add("TGel/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TGel.add("TGel/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    TGel.add("TGel/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    TGel.add("TGel/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

    // TGelCand histograms
    TGelCand.add("TGelCand/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGelCand.add("TGelCand/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGelCand.add("TGelCand/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGelCand.add("TGelCand/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGelCand.add("TGelCand/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGelCand.add("TGelCand/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGelCand.add("TGelCand/TPCNegVsPosSignal", "TPCNegVsPosSignal", HistType::kTH2F, {axisTPC, axisTPC});
    TGelCand.add("TGelCand/hTrackITSNcls1", "hTrackITSNcls1", {HistType::kTH1F, {axisCounter}});
    TGelCand.add("TGelCand/hTrackITSNcls2", "hTrackITSNcls2", {HistType::kTH1F, {axisCounter}});
    TGelCand.add("TGelCand/hPairPt", "hPairPt", {HistType::kTH1F, {axisPt}});
    TGelCand.add("TGelCand/hPairIVM", "hPairIVM", {HistType::kTH1F, {axisIVMWide}});
    TGelCand.add("TGelCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axisPt}});
    TGelCand.add("TGelCand/hJpsiRap", "hJpsiRap", {HistType::kTH1F, {axisEta}});
    TGelCand.add("TGelCand/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGelCand.add("TGelCand/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGelCand.add("TGelCand/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGelCand.add("TGelCand/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGelCand.add("TGelCand/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    TGelCand.add("TGelCand/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TGelCand.add("TGelCand/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    TGelCand.add("TGelCand/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    TGelCand.add("TGelCand/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

    // TGp histograms
    TGp.add("TGp/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGp.add("TGp/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGp.add("TGp/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGp.add("TGp/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGp.add("TGp/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGp.add("TGp/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGp.add("TGp/hNsigmaMu", "hNsigmaMu", HistType::kTH1F, {axisSigma});
    TGp.add("TGp/hNsigmaMuTOF", "hNsigmaMuTOF", HistType::kTH1F, {axisSigma});
    TGp.add("TGp/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGp.add("TGp/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGp.add("TGp/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGp.add("TGp/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGp.add("TGp/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    TGp.add("TGp/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TGp.add("TGp/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    TGp.add("TGp/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    TGp.add("TGp/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

    // TGpCand histograms
    TGpCand.add("TGpCand/hTrackPt1", "hTrackPt1", {HistType::kTH1F, {axisPt}});
    TGpCand.add("TGpCand/hTrackEta1", "hTrackEta1", {HistType::kTH1F, {axisEta}});
    TGpCand.add("TGpCand/hTrackPhi1", "hTrackPhi1", {HistType::kTH1F, {axisPhi}});
    TGpCand.add("TGpCand/hTrackPt2", "hTrackPt2", {HistType::kTH1F, {axisPt}});
    TGpCand.add("TGpCand/hTrackEta2", "hTrackEta2", {HistType::kTH1F, {axisEta}});
    TGpCand.add("TGpCand/hTrackPhi2", "hTrackPhi2", {HistType::kTH1F, {axisPhi}});
    TGpCand.add("TGpCand/hTrackITSNcls1", "hTrackITSNcls1", {HistType::kTH1F, {axisCounter}});
    TGpCand.add("TGpCand/hTrackITSNcls2", "hTrackITSNcls2", {HistType::kTH1F, {axisCounter}});
    TGpCand.add("TGpCand/hPairPt", "hPairPt", {HistType::kTH1F, {axisPt}});
    TGpCand.add("TGpCand/hPairIVM", "hPairIVM", {HistType::kTH1F, {axisIVMWide}});
    TGpCand.add("TGpCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axisPt}});
    TGpCand.add("TGpCand/hJpsiRap", "hJpsiRap", {HistType::kTH1F, {axisEta}});
    TGpCand.add("TGpCand/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    TGpCand.add("TGpCand/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    TGpCand.add("TGpCand/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    TGpCand.add("TGpCand/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    TGpCand.add("TGpCand/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    TGpCand.add("TGpCand/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    TGpCand.add("TGpCand/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    TGpCand.add("TGpCand/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    TGpCand.add("TGpCand/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

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
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

    JPsiToEl.add("JPsiToEl/Incoherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
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
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToEl.add("JPsiToEl/Incoherent/XnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

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
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

    JPsiToMu.add("JPsiToMu/Incoherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axisPt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
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
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axisEta}});
    JPsiToMu.add("JPsiToMu/Incoherent/XnOn/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisP, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axisPt, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axisEta, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

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
    JPsiToP.add("JPsiToP/Coherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    JPsiToP.add("JPsiToP/Coherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

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
    JPsiToP.add("JPsiToP/Incoherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisP, axisTOF}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hBetaTOFVsP", "hBetaTOFVsP", {HistType::kTH2F, {axisP, axisBetaTOF}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axisPt, axisTOF}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axisEta, axisTOF}});

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
    Asymmetry.add("Asymmetry/Muon/Incoherent/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Coherent/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
    Asymmetry.add("Asymmetry/Electron/Incoherent/DeltaPhi", "DeltaPhi", {HistType::kTH1F, {{180, -PI, PI}}});
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
    // acceptance
    if (std::abs(RecoDecay::eta(std::array{track.px(), track.py(), track.pz()})) > EtaCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(3);
      return false;
    }
    // DCA
    if (std::abs(track.dcaZ()) > dcaZCut) {
      return false;
    }
    if (DCAcut) {
      float dcaXYPtCut = 0.0105f + 0.0350f / pow(track.pt(), 1.1f);
      if (std::abs(track.dcaXY()) > dcaXYPtCut) {
        Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(4);
        return false;
      }
    } else {
      if (std::abs(track.dcaXY()) > dcaXYCut) {
        Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(4);
        return false;
      }
    }
    // ITS
    if (!track.hasITS()) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(5);
      return false;
    }
    if (track.itsNCls() < ITSNClsCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(6);
      return false;
    }
    if (track.itsChi2NCl() > ITSChi2NClsCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(7);
      return false;
    }
    //  TPC
    if (!track.hasTPC()) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(8);
      return false;
    }
    if (track.tpcNClsCrossedRows() < TPCNClsCrossedRowsCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(9);
      return false;
    }
    if (track.tpcChi2NCl() > TPCChi2NCls) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(10);
      return false; // TPC chi2
    }
    if (newCutTPC) {
      if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < TPCMinNCls) {
        Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(11);
        return false;
      }
      if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < TPCCrossedOverFindable) {
        Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(12);
        return false;
      }
    }

    return true;
  }

  // template <typename C>
  bool CandidateCuts(float massJpsi, float rapJpsi)
  {
    if (std::abs(rapJpsi) > RapCut) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(13);
      return false;
    }

    if (massJpsi < 2.0f) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterTracks"))->Fill(14);
      return false;
    }

    return true;
  }

  template <typename C, typename Ts>
  void fillHistograms(C collision, Ts tracks)
  {
    Statistics.get<TH1>(HIST("Statistics/hCutCounterCollisions"))->Fill(0); // number of collisions without any cuts
    RawData.get<TH1>(HIST("RawData/hPositionX"))->Fill(collision.posX());
    RawData.get<TH1>(HIST("RawData/hPositionY"))->Fill(collision.posY());
    RawData.get<TH1>(HIST("RawData/hPositionZ"))->Fill(collision.posZ());
    RawData.get<TH2>(HIST("RawData/hPositionXY"))->Fill(collision.posX(), collision.posY());
    RawData.get<TH1>(HIST("RawData/hZNACommonEnergy"))->Fill(collision.energyCommonZNA());
    RawData.get<TH1>(HIST("RawData/hZNCCommonEnergy"))->Fill(collision.energyCommonZNC());
    RawData.get<TH1>(HIST("RawData/hZNCTime"))->Fill(collision.timeZNC());
    RawData.get<TH1>(HIST("RawData/hZNATime"))->Fill(collision.timeZNA());

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
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.tofSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.beta());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTOFVsPt"))->Fill(track.pt(), track.tofSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTOFVsEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}), track.tofSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(trkPx, trkPy), track.tofSignal());
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
        RawData.get<TH2>(HIST("RawData/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.tofSignal());
        RawData.get<TH2>(HIST("RawData/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.beta());
        RawData.get<TH2>(HIST("RawData/PID/hTOFVsPt"))->Fill(track.pt(), track.tofSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTOFVsEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}), track.tofSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(trkPx, trkPy), track.tofSignal());
      }
    }

    int countGT = 0;
    int countGTMuSigma = 0;
    int countGTElSigma = 0;
    int countGTPSigma = 0;
    int countGTPSigmaTOF = 0;
    std::vector<int> trkIdx;
    // loop over tracks with selections
    if (std::abs(collision.posZ()) > cutVertexZ) {
      Statistics.get<TH1>(HIST("Statistics/hCutCounterCollisions"))->Fill(1);
      return;
    }

    for (auto& track : tracks) {

      // select good tracks
      if (GoodTrackCuts(track) != 1) {
        continue;
      }

      countGT++;
      trkIdx.push_back(track.index());

      if (std::abs(track.tpcNSigmaMu()) <= 3) {
        countGTMuSigma++;
      }
      if (std::abs(track.tpcNSigmaEl()) <= 3) {
        countGTElSigma++;
      }
      if (std::abs(track.tpcNSigmaPr()) <= 3) {
        countGTPSigma++;
      }
      if (std::abs(track.tofNSigmaPr()) <= 3) {
        countGTPSigmaTOF++;
      }
    }

    Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(2., countGT);
    Statistics.get<TH1>(HIST("Statistics/hNumberGT"))->Fill(countGT);
    Statistics.get<TH1>(HIST("Statistics/hNumberGTelSigma"))->Fill(countGTElSigma);
    Statistics.get<TH1>(HIST("Statistics/hNumberGTmuSigma"))->Fill(countGTMuSigma);
    Statistics.get<TH1>(HIST("Statistics/hNumberGTpSigma"))->Fill(countGTPSigma);
    Statistics.get<TH1>(HIST("Statistics/hNumberGTpSigmaTOF"))->Fill(countGTPSigmaTOF);

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
        TG.get<TH1>(HIST("TG/TPC/hNsigmaMu"))->Fill(trkDaughter1.tpcNSigmaMu());
        TG.get<TH1>(HIST("TG/TPC/hNsigmaEl"))->Fill(trkDaughter1.tpcNSigmaEl());
        TG.get<TH1>(HIST("TG/TPC/hNsigmaPr"))->Fill(trkDaughter1.tpcNSigmaPr());
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
        TG.get<TH1>(HIST("TG/TPC/hNsigmaMu"))->Fill(trkDaughter2.tpcNSigmaMu());
        TG.get<TH1>(HIST("TG/TPC/hNsigmaEl"))->Fill(trkDaughter2.tpcNSigmaEl());
        TG.get<TH1>(HIST("TG/TPC/hNsigmaPr"))->Fill(trkDaughter2.tpcNSigmaPr());
        TG.get<TH1>(HIST("TG/hTPCNClsFindable"))->Fill(trkDaughter2.tpcNClsFindable());
        TG.get<TH1>(HIST("TG/hTPCNClsFindableMinusFound"))->Fill(trkDaughter2.tpcNClsFindableMinusFound());
      }
      if (trkDaughter1.hasTOF()) {
        TG.get<TH2>(HIST("TG/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
        TG.get<TH2>(HIST("TG/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
        TG.get<TH2>(HIST("TG/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
        TG.get<TH2>(HIST("TG/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
        TG.get<TH1>(HIST("TG/TOF/hNsigmaMu"))->Fill(trkDaughter1.tofNSigmaMu());
        TG.get<TH1>(HIST("TG/TOF/hNsigmaEl"))->Fill(trkDaughter1.tofNSigmaEl());
        TG.get<TH1>(HIST("TG/TOF/hNsigmaPr"))->Fill(trkDaughter1.tofNSigmaPr());
      }
      if (trkDaughter2.hasTOF()) {
        TG.get<TH2>(HIST("TG/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
        TG.get<TH2>(HIST("TG/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
        TG.get<TH2>(HIST("TG/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
        TG.get<TH2>(HIST("TG/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
        TG.get<TH2>(HIST("TG/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
        TG.get<TH1>(HIST("TG/TOF/hNsigmaMu"))->Fill(trkDaughter2.tofNSigmaMu());
        TG.get<TH1>(HIST("TG/TOF/hNsigmaEl"))->Fill(trkDaughter2.tofNSigmaEl());
        TG.get<TH1>(HIST("TG/TOF/hNsigmaPr"))->Fill(trkDaughter2.tofNSigmaPr());
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
            TGel.get<TH1>(HIST("TGel/hNsigmaEl"))->Fill(trkDaughter1.tpcNSigmaEl());
          }
          if (trkDaughter2.hasTPC()) {
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
            TGel.get<TH1>(HIST("TGel/hNsigmaEl"))->Fill(trkDaughter2.tpcNSigmaEl());
            if (trkDaughter1.sign() < 0) {
              TGel.get<TH2>(HIST("TGel/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
            } else {
              TGel.get<TH2>(HIST("TGel/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
            }
          }
          if (trkDaughter1.hasTOF()) {
            TGel.get<TH2>(HIST("TGel/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
            TGel.get<TH1>(HIST("TGel/hNsigmaElTOF"))->Fill(trkDaughter1.tofNSigmaEl());
          }
          if (trkDaughter2.hasTOF()) {
            TGel.get<TH2>(HIST("TGel/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
            TGel.get<TH2>(HIST("TGel/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
            TGel.get<TH2>(HIST("TGel/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
            TGel.get<TH2>(HIST("TGel/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
            TGel.get<TH1>(HIST("TGel/hNsigmaElTOF"))->Fill(trkDaughter2.tofNSigmaEl());
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
          TGelCand.get<TH1>(HIST("TGelCand/hTrackITSNcls1"))->Fill(trkDaughter1.itsNCls());
          TGelCand.get<TH1>(HIST("TGelCand/hTrackITSNcls2"))->Fill(trkDaughter2.itsNCls());
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
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
          }
          if (trkDaughter2.hasTOF()) {
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
            TGelCand.get<TH2>(HIST("TGelCand/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
          }

          if (RecoDecay::pt(mother) < 0.2f) {
            JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hIVM"))->Fill(massJpsi);
          } else {
            JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hIVM"))->Fill(massJpsi);
          }

          if ((massJpsi < maxJpsiMass) && (massJpsi > minJpsiMass)) {
            TGelCand.get<TH1>(HIST("TGelCand/hJpsiPt"))->Fill(RecoDecay::pt(mother));
            TGelCand.get<TH1>(HIST("TGelCand/hJpsiRap"))->Fill(rapJpsi);
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
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
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
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
                JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
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
            TGmu.get<TH1>(HIST("TGmu/hNsigmaMu"))->Fill(trkDaughter1.tpcNSigmaMu());
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
            TGmu.get<TH1>(HIST("TGmu/hNsigmaMu"))->Fill(trkDaughter2.tpcNSigmaMu());
          }
          if (trkDaughter1.hasTOF()) {
            TGmu.get<TH2>(HIST("TGmu/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
            TGmu.get<TH1>(HIST("TGmu/hNsigmaMuTOF"))->Fill(trkDaughter1.tofNSigmaMu());
          }
          if (trkDaughter2.hasTOF()) {
            TGmu.get<TH2>(HIST("TGmu/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
            TGmu.get<TH2>(HIST("TGmu/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
            TGmu.get<TH2>(HIST("TGmu/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
            TGmu.get<TH1>(HIST("TGmu/hNsigmaMuTOF"))->Fill(trkDaughter1.tofNSigmaMu());
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
          TGmuCand.get<TH1>(HIST("TGmuCand/hTrackITSNcls1"))->Fill(trkDaughter1.itsNCls());
          TGmuCand.get<TH1>(HIST("TGmuCand/hTrackITSNcls2"))->Fill(trkDaughter2.itsNCls());
          TGmuCand.get<TH1>(HIST("TGmuCand/hPairPt"))->Fill(RecoDecay::pt(mother));
          TGmuCand.get<TH1>(HIST("TGmuCand/hPairIVM"))->Fill(massJpsi);

          if (trkDaughter1.hasTPC()) {
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
            if (trkDaughter1.sign() < 0) {
              TGmuCand.get<TH2>(HIST("TGmuCand/TPCNegVsPosSignal"))->Fill(trkDaughter1.tpcSignal(), trkDaughter2.tpcSignal());
            } else {
              TGmuCand.get<TH2>(HIST("TGmuCand/TPCNegVsPosSignal"))->Fill(trkDaughter2.tpcSignal(), trkDaughter1.tpcSignal());
            }
          }
          if (trkDaughter2.hasTPC()) {
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          if (trkDaughter1.hasTOF()) {
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
          }
          if (trkDaughter2.hasTOF()) {
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
            TGmuCand.get<TH2>(HIST("TGmuCand/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
          }

          if (RecoDecay::pt(mother) < 0.2f) {
            JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hIVM"))->Fill(massJpsi);
          } else {
            JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hIVM"))->Fill(massJpsi);
          }

          if ((massJpsi < maxJpsiMass) && (massJpsi > minJpsiMass)) {
            TGmuCand.get<TH1>(HIST("TGmuCand/hJpsiPt"))->Fill(RecoDecay::pt(mother));
            TGmuCand.get<TH1>(HIST("TGmuCand/hJpsiRap"))->Fill(rapJpsi);

            if (RecoDecay::pt(mother) < 0.2f) {
              // fill track histos
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
              } else if (OnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (XnOn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/XnOn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
              } else if (OnXn) {
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPt1"))->Fill(trkDaughter1.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPt2"))->Fill(trkDaughter2.pt());
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hEta1"))->Fill(RecoDecay::eta(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hEta2"))->Fill(RecoDecay::eta(daughter2));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPhi1"))->Fill(RecoDecay::phi(daughter1));
                JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/OnXn/hPhi2"))->Fill(RecoDecay::phi(daughter2));
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
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
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
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
                JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
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
            TGp.get<TH1>(HIST("TGp/hNsigmaMu"))->Fill(trkDaughter1.tpcNSigmaPr());
          }
          if (trkDaughter2.hasTPC()) {
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
            TGp.get<TH1>(HIST("TGp/hNsigmaMu"))->Fill(trkDaughter2.tpcNSigmaPr());
          }
          if (trkDaughter1.hasTOF()) {
            TGp.get<TH2>(HIST("TGp/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
            TGp.get<TH1>(HIST("TGp/hNsigmaMuTOF"))->Fill(trkDaughter1.tofNSigmaPr());
          }
          if (trkDaughter2.hasTOF()) {
            TGp.get<TH2>(HIST("TGp/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
            TGp.get<TH2>(HIST("TGp/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
            TGp.get<TH2>(HIST("TGp/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
            TGp.get<TH2>(HIST("TGp/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
            TGp.get<TH1>(HIST("TGp/hNsigmaMuTOF"))->Fill(trkDaughter2.tofNSigmaPr());
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
          TGpCand.get<TH1>(HIST("TGpCand/hTrackITSNcls1"))->Fill(trkDaughter1.itsNCls());
          TGpCand.get<TH1>(HIST("TGpCand/hTrackITSNcls2"))->Fill(trkDaughter2.itsNCls());
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
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
          }
          if (trkDaughter2.hasTOF()) {
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
            TGpCand.get<TH2>(HIST("TGpCand/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
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
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
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
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()), trkDaughter1.beta());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
              }
              if (trkDaughter2.hasTOF()) {
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hBetaTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()), trkDaughter2.beta());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
                JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
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
  } // end process

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

  PROCESS_SWITCH(UpcJpsiCentralBarrel, processDGrecoLevel, "Iterate over DG skimmed data.", true);
  PROCESS_SWITCH(UpcJpsiCentralBarrel, processSGrecoLevel, "Iterate over SG skimmed data.", false);
}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcJpsiCentralBarrel>(cfgc, TaskName{"upc-jpsi-corr"}),
  };
}
