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

// ROOT headers
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;

struct upcJpsiCentralBarrel {
  // configurable axes
  ConfigurableAxis IVMAxis{"IVMAxis", {350.0f, 0.0f, 4.5f}, "M_#it{inv} (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {250.0f, 0.1f, 3.0f}, "#it{p}_T (GeV/#it{c})"};
  ConfigurableAxis pAxis{"pAxis", {250.0f, 0.1f, 3.0f}, "#it{p} (GeV/#it{c})"};
  ConfigurableAxis etaAxis{"etaAxis", {250.0f, -1.5f, 1.5f}, "#eta (-)"};
  ConfigurableAxis countAxis{"countAxis", {10.0f, 0.0f, 10.0f}, "Number of events (-)"};
  ConfigurableAxis phiAxis{"phiAxis", {250.0f, 0, TwoPI}, "#phi (rad)"};
  ConfigurableAxis accoplAxis{"accoplAxis", {250.0f, -0.2f, 0.2f}, "accAngle"};
  ConfigurableAxis thetaAxis{"thetaAxis", {250.0f, -1.5f, 1.5f}, "cos #theta (-)"};
  ConfigurableAxis sigTPCAxis{"sigTPCAxis", {100.0f, 0, 200.0f}, "TPC d#it{E}/d#it{x}"};
  ConfigurableAxis sigTOFAxis{"sigTOFAxis", {100.0f, 0, 200.0f}, "TOF d#it{E}/d#it{x}"};

  // configurable cuts (modify in json)
  Configurable<int> TPCNClsCrossedRows{"TPCNClsCrossedRows", 70, "number of crossed rows in TPC"};
  Configurable<float> TPCNSigmaMu{"TPCNSigmaMu", 5, "PID for TPC Mu track"};
  Configurable<float> EtaCut{"EtaCut", 0.9f, "acceptance cut per track"};
  Configurable<float> RapCut{"RapCut", 0.8f, "choose event in midrapidity"};
  Configurable<float> dcaZCut{"dcaZCut", 2, "cut on the impact parameter in z of the track to the PV"};
  Configurable<float> dcaXYCut{"dcaXYCut", 1e10, "cut on the impact parameter in xy of the track to the PV"};
  Configurable<int> ITSNClsCut{"ITSNClsCut", 1, "minimal number of ITS clusters"};
  Configurable<int> ITSChi2NClsCut{"ITSChi2NClsCut", 36, "minimal Chi2/cluster for the ITS track"};
  Configurable<int> TPCNClsCrossedRowsCut{"TPCNClsCrossedRowsCut", 70, "minimal number of crossed TPC rows"};
  Configurable<int> TPCChi2NCls{"TPCChi2NCls", 4, "minimal Chi2/cluster for the TPC track"};

  // initialize histogram registry
  HistogramRegistry Statistics{
    "Statistics",
    {}};

  HistogramRegistry RawData{
    "RawData",
    {}};

  HistogramRegistry PVContributors{
    "PVContributors",
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

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDTrackFull = UDTracksFull::iterator;

  void init(InitContext&)
  {

    const AxisSpec axisIVM{IVMAxis, "IVM axis"};
    const AxisSpec axispt{ptAxis, "pt axis"};
    const AxisSpec axiseta{etaAxis, "eta axis"};
    const AxisSpec axisCounter{countAxis, "counter axis"};
    const AxisSpec axisAccAngle{accoplAxis, "accAngle"};
    const AxisSpec axisAngTheta{thetaAxis, "cosTheta"};
    const AxisSpec axisPhi{phiAxis, "phi"};
    const AxisSpec axisp{pAxis, "p axis"};
    const AxisSpec axisTPC{sigTPCAxis, ""};
    const AxisSpec axisTOF{sigTOFAxis, ""};

    // statistics histograms for counters
    Statistics.add("Statistics/hNumberOfCollisions", "hNumberOfCollisions", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberOfTracks", "hNumberOfTracks", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGT", "hNumberGT", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTselected", "hNumberGTselected", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTel", "hNumberGTel", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTmu", "hNumberGTmu", {HistType::kTH1F, {axisCounter}});
    Statistics.add("Statistics/hNumberGTp", "hNumberGTp", {HistType::kTH1F, {axisCounter}});

    // raw data histograms
    RawData.add("RawData/hTrackPt", "hTrackPt", {HistType::kTH1F, {axispt}});
    RawData.add("RawData/hTrackEta", "hTrackEta", {HistType::kTH1F, {axiseta}});
    RawData.add("RawData/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});
    RawData.add("RawData/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisp, axisTPC}});
    RawData.add("RawData/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axispt, axisTPC}});
    RawData.add("RawData/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    RawData.add("RawData/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axiseta, axisTPC}});
    RawData.add("RawData/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisp, axisTOF}});
    RawData.add("RawData/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axispt, axisTOF}});
    RawData.add("RawData/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    RawData.add("RawData/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axiseta, axisTOF}});

    // PVContributors histograms
    PVContributors.add("PVContributors/hTrackPt", "hTrackPt", {HistType::kTH1F, {axispt}});
    PVContributors.add("PVContributors/hTrackEta", "hTrackEta", {HistType::kTH1F, {axiseta}});
    PVContributors.add("PVContributors/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});
    PVContributors.add("PVContributors/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisp, axisTPC}});
    PVContributors.add("PVContributors/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axispt, axisTPC}});
    PVContributors.add("PVContributors/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    PVContributors.add("PVContributors/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axiseta, axisTPC}});
    PVContributors.add("PVContributors/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisp, axisTOF}});
    PVContributors.add("PVContributors/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axispt, axisTOF}});
    PVContributors.add("PVContributors/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    PVContributors.add("PVContributors/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axiseta, axisTOF}});

    // TGmu histograms
    TGmu.add("TGmu/hTrackPt", "hTrackPt", {HistType::kTH1F, {axispt}});
    TGmu.add("TGmu/hTrackEta", "hTrackEta", {HistType::kTH1F, {axiseta}});
    TGmu.add("TGmu/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});

    // TGmuCand histograms
    TGmuCand.add("TGmuCand/hTrackPt", "hTrackPt", {HistType::kTH1F, {axispt}});
    TGmuCand.add("TGmuCand/hTrackEta", "hTrackEta", {HistType::kTH1F, {axiseta}});
    TGmuCand.add("TGmuCand/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});
    TGmuCand.add("TGmuCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axispt}});

    // TGel histograms
    TGel.add("TGel/hTrackPt", "hTrackPt", {HistType::kTH1F, {axispt}});
    TGel.add("TGel/hTrackEta", "hTrackEta", {HistType::kTH1F, {axiseta}});
    TGel.add("TGel/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});

    // TGelCand histograms
    TGelCand.add("TGelCand/hTrackPt", "hTrackPt", {HistType::kTH1F, {axispt}});
    TGelCand.add("TGelCand/hTrackEta", "hTrackEta", {HistType::kTH1F, {axiseta}});
    TGelCand.add("TGelCand/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});
    TGelCand.add("TGelCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axispt}});

    // TGp histograms
    TGp.add("TGp/hTrackPt", "hTrackPt", {HistType::kTH1F, {axispt}});
    TGp.add("TGp/hTrackEta", "hTrackEta", {HistType::kTH1F, {axiseta}});
    TGp.add("TGp/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});

    // TGpCand histograms
    TGpCand.add("TGpCand/hTrackPt", "hTrackPt", {HistType::kTH1F, {axispt}});
    TGpCand.add("TGpCand/hTrackEta", "hTrackEta", {HistType::kTH1F, {axiseta}});
    TGpCand.add("TGpCand/hTrackPhi", "hTrackPhi", {HistType::kTH1F, {axisPhi}});
    TGpCand.add("TGpCand/hJpsiPt", "hJpsiPt", {HistType::kTH1F, {axispt}});

    // JPsiToEl histograms
    JPsiToEl.add("JPsiToEl/Coherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToEl.add("JPsiToEl/Coherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToEl.add("JPsiToEl/Coherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToEl.add("JPsiToEl/Coherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToEl.add("JPsiToEl/Coherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToEl.add("JPsiToEl/Coherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Coherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axiseta}});
    JPsiToEl.add("JPsiToEl/Coherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToEl.add("JPsiToEl/Coherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisp, axisTPC}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axispt, axisTPC}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axiseta, axisTPC}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisp, axisTOF}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axispt, axisTOF}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToEl.add("JPsiToEl/Coherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axiseta, axisTOF}});

    JPsiToEl.add("JPsiToEl/Incoherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToEl.add("JPsiToEl/Incoherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToEl.add("JPsiToEl/Incoherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axiseta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToEl.add("JPsiToEl/Incoherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisp, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axispt, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axiseta, axisTPC}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisp, axisTOF}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axispt, axisTOF}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToEl.add("JPsiToEl/Incoherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axiseta, axisTOF}});

    // JPsiToMu histograms
    JPsiToMu.add("JPsiToMu/Coherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToMu.add("JPsiToMu/Coherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToMu.add("JPsiToMu/Coherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToMu.add("JPsiToMu/Coherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToMu.add("JPsiToMu/Coherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToMu.add("JPsiToMu/Coherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Coherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axiseta}});
    JPsiToMu.add("JPsiToMu/Coherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToMu.add("JPsiToMu/Coherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisp, axisTPC}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axispt, axisTPC}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axiseta, axisTPC}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisp, axisTOF}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axispt, axisTOF}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToMu.add("JPsiToMu/Coherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axiseta, axisTOF}});

    JPsiToMu.add("JPsiToMu/Incoherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToMu.add("JPsiToMu/Incoherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToMu.add("JPsiToMu/Incoherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axiseta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToMu.add("JPsiToMu/Incoherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisp, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axispt, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axiseta, axisTPC}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisp, axisTOF}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axispt, axisTOF}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToMu.add("JPsiToMu/Incoherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axiseta, axisTOF}});

    // JPsiToP histograms
    JPsiToP.add("JPsiToP/Coherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToP.add("JPsiToP/Coherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToP.add("JPsiToP/Coherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToP.add("JPsiToP/Coherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToP.add("JPsiToP/Coherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToP.add("JPsiToP/Coherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Coherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Coherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToP.add("JPsiToP/Coherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axiseta}});
    JPsiToP.add("JPsiToP/Coherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToP.add("JPsiToP/Coherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisp, axisTPC}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axispt, axisTPC}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axiseta, axisTPC}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisp, axisTOF}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axispt, axisTOF}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToP.add("JPsiToP/Coherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axiseta, axisTOF}});

    JPsiToP.add("JPsiToP/Incoherent/hPt", "Pt of J/Psi ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToP.add("JPsiToP/Incoherent/hPt1", "pT of track 1 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToP.add("JPsiToP/Incoherent/hPt2", "pT of track 2 ; p_{T} {GeV/c]", {HistType::kTH1F, {axispt}});
    JPsiToP.add("JPsiToP/Incoherent/hEta1", "eta of track 1 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToP.add("JPsiToP/Incoherent/hEta2", "eta of track 2 ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToP.add("JPsiToP/Incoherent/hPhi1", "phi of track 1 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Incoherent/hPhi2", "phi of track 2 ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Incoherent/hIVM", "J/Psi Invariant Mass ; m {GeV]", {HistType::kTH1F, {axisIVM}});
    JPsiToP.add("JPsiToP/Incoherent/hRap", "Rap of J/Psi ; y {-]", {HistType::kTH1F, {axiseta}});
    JPsiToP.add("JPsiToP/Incoherent/hEta", "Eta of J/Psi ; #eta {-]", {HistType::kTH1F, {axiseta}});
    JPsiToP.add("JPsiToP/Incoherent/hPhi", "Phi of J/Psi ; #phi {-]", {HistType::kTH1F, {axisPhi}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTPCVsP", "hTPCVsP", {HistType::kTH2F, {axisp, axisTPC}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTPCVsPt", "hTPCVsPt", {HistType::kTH2F, {axispt, axisTPC}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTPCVsPhi", "hTPCVsPhi", {HistType::kTH2F, {axisPhi, axisTPC}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTPCVsEta", "hTPCVsEta", {HistType::kTH2F, {axiseta, axisTPC}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTOFVsP", "hTOFVsP", {HistType::kTH2F, {axisp, axisTOF}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTOFVsPt", "hTOFVsPt", {HistType::kTH2F, {axispt, axisTOF}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTOFVsPhi", "hTOFVsPhi", {HistType::kTH2F, {axisPhi, axisTOF}});
    JPsiToP.add("JPsiToP/Incoherent/PID/hTOFVsEta", "hTOFVsEta", {HistType::kTH2F, {axiseta, axisTOF}});

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
    // kinematics
    if (std::abs(RecoDecay::eta(std::array{track.px(), track.py(), track.pz()})) > EtaCut) {
      return false;
    }
    // DCA
    if (track.dcaZ() > dcaZCut) {
      return false;
    }
    if (track.dcaXY() > dcaXYCut) {
      return false;
    }
    // ITS
    if (!track.hasITS()) {
      return false;
    }
    if (track.itsNCls() < ITSNClsCut) {
      return false;
    }
    if (track.itsChi2NCl() > ITSChi2NClsCut) {
      return false;
    }
    //  TPC
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < TPCNClsCrossedRowsCut) {
      return false;
    }
    if (track.tpcChi2NCl() > TPCChi2NCls) {
      return false; // TPC chi2
    }

    return true;
  }

  // template <typename C>
  bool CandidateCuts(float massJpsi, float rapJpsi)
  {
    if (rapJpsi > RapCut) {
      return false;
    }

    if (massJpsi < 2.5f) {
      return false;
    }

    return true;
  }

  void process(UDCollisionFull const&, UDTracksFull const& tracks)
  {
    Statistics.get<TH1>(HIST("Statistics/hNumberOfCollisions"))->Fill(0); // number of collisions without any cuts

    // loop over tracks without selections
    for (auto& track : tracks) {
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();

      Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(0);
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
        }

        if (track.hasTOF()) {
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.tofSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTOFVsPt"))->Fill(track.pt(), track.tofSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTOFVsEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}), track.tofSignal());
          PVContributors.get<TH2>(HIST("PVContributors/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(trkPx, trkPy), track.tofSignal());
        }
      }

      RawData.get<TH1>(HIST("RawData/hTrackPt"))->Fill(track.pt());
      RawData.get<TH1>(HIST("RawData/hTrackEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}));
      RawData.get<TH1>(HIST("RawData/hTrackPhi"))->Fill(RecoDecay::phi(trkPx, trkPy));

      if (track.hasTPC()) {
        RawData.get<TH2>(HIST("RawData/PID/hTPCVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.tpcSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTPCVsPt"))->Fill(track.pt(), track.tpcSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTPCVsEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}), track.tpcSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(trkPx, trkPy), track.tpcSignal());
      }

      if (track.hasTOF()) {
        RawData.get<TH2>(HIST("RawData/PID/hTOFVsP"))->Fill(RecoDecay::sqrtSumOfSquares(trkPx, trkPy, trkPz), track.tofSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTOFVsPt"))->Fill(track.pt(), track.tofSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTOFVsEta"))->Fill(RecoDecay::eta(std::array{trkPx, trkPy, trkPz}), track.tofSignal());
        RawData.get<TH2>(HIST("RawData/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(trkPx, trkPy), track.tofSignal());
      }
    }

    int countGT = 0;
    int countGTselected = 0;
    int countGTel = 0;
    int countGTmu = 0;
    int countGTp = 0;
    std::vector<int> trkIdx;
    // loop over tracks with selections
    for (auto& track : tracks) {
      // select primary vertex contributors
      if (track.isPVContributor() != 1) {
        return;
      }
      // select good tracks
      if (GoodTrackCuts(track) != 1) {
        return;
      }
      countGT++;
      int hypoID = testPIDhypo(track);
      if (hypoID == P_ELECTRON || hypoID == P_MUON) {
        countGTselected++;
        trkIdx.push_back(track.index());
        if (hypoID == P_ELECTRON) {
          countGTel++;
        }
        if (hypoID == P_MUON) {
          countGTmu++;
        }
        if (hypoID == P_PROTON) {
          countGTp++;
        }
      }
    }

    Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(2., countGT);
    Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(3., countGTselected);
    Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(4., countGTel);
    Statistics.get<TH1>(HIST("Statistics/hNumberOfTracks"))->Fill(5., countGTmu);
    Statistics.get<TH1>(HIST("Statistics/hNumberGT"))->Fill(countGT);
    Statistics.get<TH1>(HIST("Statistics/hNumberGTselected"))->Fill(countGTselected);
    Statistics.get<TH1>(HIST("Statistics/hNumberGTel"))->Fill(countGTel);
    Statistics.get<TH1>(HIST("Statistics/hNumberGTmu"))->Fill(countGTmu);
    Statistics.get<TH1>(HIST("Statistics/hNumberGTp"))->Fill(countGTp);

    float massEl = o2::constants::physics::MassElectron;
    float massMu = o2::constants::physics::MassMuonMinus;
    float massPr = o2::constants::physics::MassProton;

    if (countGT == 2) {
      if (countGTel == 2) {
        TLorentzVector mom, daughter[2];
        auto trkDaughter1 = tracks.iteratorAt(trkIdx[0]);
        auto trkDaughter2 = tracks.iteratorAt(trkIdx[1]);
        if ((trkDaughter1.sign() * trkDaughter2.sign()) > 0) {
          return;
        }
        if (!(RecoDecay::sumOfSquares(trkDaughter1.tpcNSigmaMu(), trkDaughter2.tpcNSigmaMu()) < RecoDecay::sumOfSquares(trkDaughter1.tpcNSigmaEl(), trkDaughter2.tpcNSigmaEl()))) {
          return;
        }
        auto ene1 = RecoDecay::e(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), massEl);
        auto ene2 = RecoDecay::e(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), massEl);
        daughter[0].SetPxPyPzE(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), ene1);
        daughter[1].SetPxPyPzE(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), ene2);
        mom = daughter[0] + daughter[1];

        std::array<double, 3> daughter1 = {trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()};
        std::array<double, 3> daughter2 = {trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()};

        std::array<double, 3> mother = {trkDaughter1.px() + trkDaughter2.px(), trkDaughter1.py() + trkDaughter2.py(), trkDaughter1.pz() + trkDaughter2.pz()};

        auto arrMom = std::array{daughter1, daughter2};
        float massJpsi = RecoDecay::m(arrMom, std::array{massEl, massEl});
        float rapJpsi = RecoDecay::y(mother, massJpsi);

        TGel.get<TH1>(HIST("TGel/hTrackPt"))->Fill(trkDaughter1.pt());
        TGel.get<TH1>(HIST("TGel/hTrackPt"))->Fill(trkDaughter2.pt());
        TGel.get<TH1>(HIST("TGel/hTrackEta"))->Fill(RecoDecay::eta(daughter1));
        TGel.get<TH1>(HIST("TGel/hTrackEta"))->Fill(RecoDecay::eta(daughter2));
        TGel.get<TH1>(HIST("TGel/hTrackPhi"))->Fill(RecoDecay::phi(daughter1));
        TGel.get<TH1>(HIST("TGel/hTrackPhi"))->Fill(RecoDecay::phi(daughter2));

        if (CandidateCuts(massJpsi, rapJpsi) != 1) {
          return;
        }

        TGelCand.get<TH1>(HIST("TGelCand/hTrackPt"))->Fill(trkDaughter1.pt());
        TGelCand.get<TH1>(HIST("TGelCand/hTrackPt"))->Fill(trkDaughter2.pt());
        TGelCand.get<TH1>(HIST("TGelCand/hTrackEta"))->Fill(RecoDecay::eta(daughter1));
        TGelCand.get<TH1>(HIST("TGelCand/hTrackEta"))->Fill(RecoDecay::eta(daughter2));
        TGelCand.get<TH1>(HIST("TGelCand/hTrackPhi"))->Fill(RecoDecay::phi(daughter1));
        TGelCand.get<TH1>(HIST("TGelCand/hTrackPhi"))->Fill(RecoDecay::phi(daughter2));
        TGelCand.get<TH1>(HIST("TGelCand/hJpsiPt"))->Fill(RecoDecay::pt(mother));
        if (RecoDecay::pt(mother) < 0.2f) {
          // fill track histos
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPt1"))->Fill(trkDaughter1.pt());
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPt2"))->Fill(trkDaughter2.pt());
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));
          if (trkDaughter1.hasTPC()) {
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }

          if (trkDaughter1.hasTOF()) {
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
          }
          if (trkDaughter2.hasTOF()) {
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Coherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
          }
          // fill J/psi histos
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPt"))->Fill(RecoDecay::pt(mother));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hEta"))->Fill(RecoDecay::eta(mother));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hPhi"))->Fill(RecoDecay::phi(mother));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hRap"))->Fill(rapJpsi);
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Coherent/hIVM"))->Fill(massJpsi);

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
        }
        if (RecoDecay::pt(mother) > 0.2f) {
          // fill track histos
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPt1"))->Fill(trkDaughter1.pt());
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPt2"))->Fill(trkDaughter2.pt());
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));
          if (trkDaughter1.hasTPC()) {
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }

          if (trkDaughter1.hasTOF()) {
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tofSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tofSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tofSignal());
          }
          if (trkDaughter2.hasTOF()) {
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tofSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tofSignal());
            JPsiToEl.get<TH2>(HIST("JPsiToEl/Incoherent/PID/hTOFVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tofSignal());
          }
          // fill J/psi histos
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPt"))->Fill(RecoDecay::pt(mother));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hEta"))->Fill(RecoDecay::eta(mother));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hPhi"))->Fill(RecoDecay::phi(mother));
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hRap"))->Fill(rapJpsi);
          JPsiToEl.get<TH1>(HIST("JPsiToEl/Incoherent/hIVM"))->Fill(massJpsi);

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
        }
      } // end electrons
      if (countGTmu == 2) {
        TLorentzVector mom, daughter[2];
        auto trkDaughter1 = tracks.iteratorAt(trkIdx[0]);
        auto trkDaughter2 = tracks.iteratorAt(trkIdx[1]);
        if ((trkDaughter1.sign() * trkDaughter2.sign()) > 0) {
          return;
        }
        if (!(RecoDecay::sumOfSquares(trkDaughter1.tpcNSigmaEl(), trkDaughter2.tpcNSigmaEl() < RecoDecay::sumOfSquares(trkDaughter1.tpcNSigmaMu(), trkDaughter2.tpcNSigmaMu())))) {
          return;
        }
        auto ene1 = RecoDecay::e(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), massMu);
        auto ene2 = RecoDecay::e(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), massMu);
        daughter[0].SetPxPyPzE(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), ene1);
        daughter[1].SetPxPyPzE(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), ene2);
        mom = daughter[0] + daughter[1];

        std::array<double, 3> daughter1 = {trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()};
        std::array<double, 3> daughter2 = {trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()};

        std::array<double, 3> mother = {trkDaughter1.px() + trkDaughter2.px(), trkDaughter1.py() + trkDaughter2.py(), trkDaughter1.pz() + trkDaughter2.pz()};

        auto arrMom = std::array{daughter1, daughter2};
        float massJpsi = RecoDecay::m(arrMom, std::array{massMu, massMu});
        float rapJpsi = RecoDecay::y(mother, massJpsi);

        TGmu.get<TH1>(HIST("TGmu/hTrackPt"))->Fill(trkDaughter1.pt());
        TGmu.get<TH1>(HIST("TGmu/hTrackPt"))->Fill(trkDaughter2.pt());
        TGmu.get<TH1>(HIST("TGmu/hTrackEta"))->Fill(RecoDecay::eta(daughter1));
        TGmu.get<TH1>(HIST("TGmu/hTrackEta"))->Fill(RecoDecay::eta(daughter2));
        TGmu.get<TH1>(HIST("TGmu/hTrackPhi"))->Fill(RecoDecay::phi(daughter1));
        TGmu.get<TH1>(HIST("TGmu/hTrackPhi"))->Fill(RecoDecay::phi(daughter2));

        if (CandidateCuts(massJpsi, rapJpsi) != 1) {
          return;
        }

        TGmuCand.get<TH1>(HIST("TGmuCand/hTrackPt"))->Fill(trkDaughter1.pt());
        TGmuCand.get<TH1>(HIST("TGmuCand/hTrackPt"))->Fill(trkDaughter2.pt());
        TGmuCand.get<TH1>(HIST("TGmuCand/hTrackEta"))->Fill(RecoDecay::eta(daughter1));
        TGmuCand.get<TH1>(HIST("TGmuCand/hTrackEta"))->Fill(RecoDecay::eta(daughter2));
        TGmuCand.get<TH1>(HIST("TGmuCand/hTrackPhi"))->Fill(RecoDecay::phi(daughter1));
        TGmuCand.get<TH1>(HIST("TGmuCand/hTrackPhi"))->Fill(RecoDecay::phi(daughter2));
        TGmuCand.get<TH1>(HIST("TGmuCand/hJpsiPt"))->Fill(RecoDecay::pt(mother));

        if (RecoDecay::pt(mother) < 0.2f) {
          // fill track histos
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPt1"))->Fill(trkDaughter1.pt());
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPt2"))->Fill(trkDaughter2.pt());
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));
          if (trkDaughter1.hasTPC()) {
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          // fill J/psi histos
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPt"))->Fill(RecoDecay::pt(mother));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hEta"))->Fill(RecoDecay::eta(mother));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hPhi"))->Fill(RecoDecay::phi(mother));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hRap"))->Fill(rapJpsi);
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Coherent/hIVM"))->Fill(massJpsi);

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
          if (trkDaughter1.hasTPC()) {
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            JPsiToMu.get<TH2>(HIST("JPsiToMu/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          // fill J/psi histos
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hPt"))->Fill(RecoDecay::pt(mother));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hEta"))->Fill(RecoDecay::eta(mother));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hPhi"))->Fill(RecoDecay::phi(mother));
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hRap"))->Fill(rapJpsi);
          JPsiToMu.get<TH1>(HIST("JPsiToMu/Incoherent/hIVM"))->Fill(massJpsi);

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
      } // end muons
      if (countGTp == 2) {
        TLorentzVector mom, daughter[2];
        auto trkDaughter1 = tracks.iteratorAt(trkIdx[0]);
        auto trkDaughter2 = tracks.iteratorAt(trkIdx[1]);
        if ((trkDaughter1.sign() * trkDaughter2.sign()) > 0) {
          return;
        }

        auto ene1 = RecoDecay::e(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), massPr);
        auto ene2 = RecoDecay::e(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), massPr);
        daughter[0].SetPxPyPzE(trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz(), ene1);
        daughter[1].SetPxPyPzE(trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz(), ene2);
        mom = daughter[0] + daughter[1];

        std::array<double, 3> daughter1 = {trkDaughter1.px(), trkDaughter1.py(), trkDaughter1.pz()};
        std::array<double, 3> daughter2 = {trkDaughter2.px(), trkDaughter2.py(), trkDaughter2.pz()};

        std::array<double, 3> mother = {trkDaughter1.px() + trkDaughter2.px(), trkDaughter1.py() + trkDaughter2.py(), trkDaughter1.pz() + trkDaughter2.pz()};

        auto arrMom = std::array{daughter1, daughter2};
        float massJpsi = RecoDecay::m(arrMom, std::array{massPr, massPr});
        float rapJpsi = RecoDecay::y(mother, massJpsi);

        TGp.get<TH1>(HIST("TGp/hTrackPt"))->Fill(trkDaughter1.pt());
        TGp.get<TH1>(HIST("TGp/hTrackPt"))->Fill(trkDaughter2.pt());
        TGp.get<TH1>(HIST("TGp/hTrackEta"))->Fill(RecoDecay::eta(daughter1));
        TGp.get<TH1>(HIST("TGp/hTrackEta"))->Fill(RecoDecay::eta(daughter2));
        TGp.get<TH1>(HIST("TGp/hTrackPhi"))->Fill(RecoDecay::phi(daughter1));
        TGp.get<TH1>(HIST("TGp/hTrackPhi"))->Fill(RecoDecay::phi(daughter2));

        if (CandidateCuts(massJpsi, rapJpsi) != 1) {
          return;
        }

        TGpCand.get<TH1>(HIST("TGpCand/hTrackPt"))->Fill(trkDaughter1.pt());
        TGpCand.get<TH1>(HIST("TGpCand/hTrackPt"))->Fill(trkDaughter2.pt());
        TGpCand.get<TH1>(HIST("TGpCand/hTrackEta"))->Fill(RecoDecay::eta(daughter1));
        TGpCand.get<TH1>(HIST("TGpCand/hTrackEta"))->Fill(RecoDecay::eta(daughter2));
        TGpCand.get<TH1>(HIST("TGpCand/hTrackPhi"))->Fill(RecoDecay::phi(daughter1));
        TGpCand.get<TH1>(HIST("TGpCand/hTrackPhi"))->Fill(RecoDecay::phi(daughter2));
        TGpCand.get<TH1>(HIST("TGpCand/hJpsiPt"))->Fill(RecoDecay::pt(mother));

        if (RecoDecay::pt(mother) < 0.2f) {
          // fill track histos
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPt1"))->Fill(trkDaughter1.pt());
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPt2"))->Fill(trkDaughter2.pt());
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));
          if (trkDaughter1.hasTPC()) {
            JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            JPsiToP.get<TH2>(HIST("JPsiToP/Coherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          // fill J/psi histos
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPt"))->Fill(RecoDecay::pt(mother));
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hEta"))->Fill(RecoDecay::eta(mother));
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hPhi"))->Fill(RecoDecay::phi(mother));
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hRap"))->Fill(rapJpsi);
          JPsiToP.get<TH1>(HIST("JPsiToP/Coherent/hIVM"))->Fill(massJpsi);
        }
        if (RecoDecay::pt(mother) > 0.2f) {
          // fill track histos
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPt1"))->Fill(trkDaughter1.pt());
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPt2"))->Fill(trkDaughter2.pt());
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hEta1"))->Fill(RecoDecay::eta(daughter1));
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hEta2"))->Fill(RecoDecay::eta(daughter2));
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPhi1"))->Fill(RecoDecay::phi(daughter1));
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPhi2"))->Fill(RecoDecay::phi(daughter2));
          if (trkDaughter1.hasTPC()) {
            JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter1.pt(), trkDaughter1.tpcSignal());
            JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter1), trkDaughter1.tpcSignal());
            JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter1), trkDaughter1.tpcSignal());
          }
          if (trkDaughter2.hasTPC()) {
            JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsPt"))->Fill(trkDaughter2.pt(), trkDaughter2.tpcSignal());
            JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsEta"))->Fill(RecoDecay::eta(daughter2), trkDaughter2.tpcSignal());
            JPsiToP.get<TH2>(HIST("JPsiToP/Incoherent/PID/hTPCVsPhi"))->Fill(RecoDecay::phi(daughter2), trkDaughter2.tpcSignal());
          }
          // fill J/psi histos
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPt"))->Fill(RecoDecay::pt(mother));
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hEta"))->Fill(RecoDecay::eta(mother));
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hPhi"))->Fill(RecoDecay::phi(mother));
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hRap"))->Fill(rapJpsi);
          JPsiToP.get<TH1>(HIST("JPsiToP/Incoherent/hIVM"))->Fill(massJpsi);
        }
      } // end protons
    }   // end two tracks
  }     // end process
};      // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<upcJpsiCentralBarrel>(cfgc, TaskName{"upc-jpsi-corr"}),
  };
}
