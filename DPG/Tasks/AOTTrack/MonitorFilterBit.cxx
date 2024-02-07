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

/// \file   MonitorFilterBit.cxx
/// \author Andrea Rossi <andrea.rossi@cern.ch>
///
/// \brief Task performing basic checks on filter-bit selections.
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::track;
using namespace o2::aod::mctracklabel;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace trackPairForEff
{
DECLARE_SOA_COLUMN(PtTPCtr, ptTPCtr, float);
DECLARE_SOA_COLUMN(EtaTPCtr, etaTPCtr, float);
DECLARE_SOA_COLUMN(PhiTPCtr, phiTPCtr, float);
DECLARE_SOA_COLUMN(PtITStr, ptITStr, float);
DECLARE_SOA_COLUMN(EtaITStr, etaITStr, float);
DECLARE_SOA_COLUMN(PhiITStr, phiITStr, float);
DECLARE_SOA_COLUMN(NClustITS, nClustITS, uint8_t);
DECLARE_SOA_COLUMN(NClustTPC, nClustTPC, int16_t);
DECLARE_SOA_COLUMN(McPtIfisSamePart, mcPtIfisSamePart, float);
DECLARE_SOA_COLUMN(PairType, pairType, uint8_t);
} // namespace trackPairForEff
DECLARE_SOA_TABLE(TrackPairForEffPP, "AOD", "TRACKPAIREFFPP",
                  trackPairForEff::PtTPCtr, trackPairForEff::EtaTPCtr, trackPairForEff::PhiTPCtr,
                  trackPairForEff::PtITStr, trackPairForEff::EtaITStr, trackPairForEff::PhiITStr,
                  trackPairForEff::NClustITS, trackPairForEff::NClustTPC, trackPairForEff::McPtIfisSamePart, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(TrackPairForEffNN, "AOD", "TRACKPAIREFFNN",
                  trackPairForEff::PtTPCtr, trackPairForEff::EtaTPCtr, trackPairForEff::PhiTPCtr,
                  trackPairForEff::PtITStr, trackPairForEff::EtaITStr, trackPairForEff::PhiITStr,
                  trackPairForEff::NClustITS, trackPairForEff::NClustTPC, trackPairForEff::McPtIfisSamePart, o2::soa::Marker<2>);
DECLARE_SOA_TABLE(TrackPairForEffPN, "AOD", "TRACKPAIREFFPN",
                  trackPairForEff::PtTPCtr, trackPairForEff::EtaTPCtr, trackPairForEff::PhiTPCtr,
                  trackPairForEff::PtITStr, trackPairForEff::EtaITStr, trackPairForEff::PhiITStr,
                  trackPairForEff::NClustITS, trackPairForEff::NClustTPC, trackPairForEff::McPtIfisSamePart, o2::soa::Marker<3>);
DECLARE_SOA_TABLE(TrackPairForEffNP, "AOD", "TRACKPAIREFFNP",
                  trackPairForEff::PtTPCtr, trackPairForEff::EtaTPCtr, trackPairForEff::PhiTPCtr,
                  trackPairForEff::PtITStr, trackPairForEff::EtaITStr, trackPairForEff::PhiITStr,
                  trackPairForEff::NClustITS, trackPairForEff::NClustTPC, trackPairForEff::McPtIfisSamePart, o2::soa::Marker<4>);
} // namespace o2::aod

struct CheckFilterBit {

  Produces<aod::TrackPairForEffPP> trackPairForEffTablePP;
  Produces<aod::TrackPairForEffNN> trackPairForEffTableNN;
  Produces<aod::TrackPairForEffNP> trackPairForEffTableNP;
  Produces<aod::TrackPairForEffPN> trackPairForEffTablePN;

  // Binning
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};
  ConfigurableAxis binsEta{"binsEta", {30, -1.5, 1.5}, ""};
  ConfigurableAxis binsLength{"binsLength", {VARIABLE_WIDTH, 0.0, 10., 20., 50., 80., 120., 160., 200., 240., 280., 320., 380, 420, 460., 500., 550}, ""};
  Configurable<float> zVtxCut{"zVtxCut", 10., "Primary Vtx z cut"};
  ConfigurableAxis binsPhi{"binsPhi", {180, 0., 2 * M_PI}, "Phi binning"};
  ConfigurableAxis binsTPCITSmatching{"binsTPCITSmatching", {2, 0.5, 2.5}, "ITSTPCmatching"};
  ConfigurableAxis binsNclustTPC{"binsNclustTPC", {VARIABLE_WIDTH, -0.5, 0.5, 10, 30, 50, 60, 70, 80, 90, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160}, ""};
  ConfigurableAxis binsRxy{"binsRxy", {VARIABLE_WIDTH, 0.0, 0.5, 1., 1.5, 2., 2.5, 3., 4., 10., 10., 20., 50., 80., 120., 160., 200., 240., 280., 320., 380, 420, 460., 500., 550}, ""};
  ConfigurableAxis binsIsPropagaed{"binsIsPropaged", {2, -0.5, 1.5}, "isTrackWithinBeamPipe"};

  HistogramRegistry histos;
  Int_t ncollisionCounter = 0;
  float fzero = 0.;
  using Tracksextension = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA>;
  using TracksextensionMC = soa::Join<Tracksextension, aod::McTrackLabels>;
  SliceCache cache;
  Partition<Tracksextension> positiveTPConlyTracks = o2::aod::track::signed1Pt > fzero&& o2::aod::track::tpcNClsFindable > (uint8_t)0 && o2::aod::track::itsChi2NCl < (float_t)0 && o2::aod::track::trackType == (uint8_t)o2::aod::track::TrackTypeEnum::Track;  //&& (o2::aod::track::detectorMap & o2::aod::track::TPC) ==o2::aod::track::TPC && (o2::aod::track::detectorMap & o2::aod::track::ITS) ==0;
  Partition<Tracksextension> negativeTPConlyTracks = o2::aod::track::signed1Pt < fzero && o2::aod::track::tpcNClsFindable > (uint8_t)0 && o2::aod::track::itsChi2NCl < (float_t)0 && o2::aod::track::trackType == (uint8_t)o2::aod::track::TrackTypeEnum::Track; //&& (o2::aod::track::detectorMap & o2::aod::track::TPC) ==o2::aod::track::TPC && (o2::aod::track::detectorMap & o2::aod::track::ITS) ==0;
  Partition<Tracksextension> positiveITSonlyTracks = o2::aod::track::signed1Pt > fzero&& o2::aod::track::tpcChi2NCl<(float_t)0 && o2::aod::track::itsChi2NCl>(float_t) 0 && o2::aod::track::trackType == (uint8_t)o2::aod::track::TrackTypeEnum::Track;          // && (o2::aod::track::detectorMap & o2::aod::track::TPC) ==0 && (o2::aod::track::detectorMap & o2::aod::track::ITS) ==o2::aod::track::ITS && o2::aod::track::passedITSNCls==true;// && o2::aod::track::itsNCls==7;
  Partition<Tracksextension> negativeITSonlyTracks = o2::aod::track::signed1Pt < fzero && o2::aod::track::tpcChi2NCl < (float_t)0 && o2::aod::track::itsChi2NCl > (float_t)0 && o2::aod::track::trackType == (uint8_t)o2::aod::track::TrackTypeEnum::Track;      // && (o2::aod::track::detectorMap & o2::aod::track::TPC)  ==0 &&(o2::aod::track::detectorMap & o2::aod::track::ITS) ==o2::aod::track::ITS && o2::aod::track::passedITSNCls==true;// && o2::aod::track::itsNCls==7;

  Partition<TracksextensionMC> positiveTPConlyTracksMC = o2::aod::track::signed1Pt > fzero&& o2::aod::track::tpcNClsFindable > (uint8_t)0 && o2::aod::track::itsChi2NCl < (float_t)0 && o2::aod::track::trackType == (uint8_t)o2::aod::track::TrackTypeEnum::Track;  //&& (o2::aod::track::detectorMap & o2::aod::track::TPC) ==o2::aod::track::TPC && (o2::aod::track::detectorMap & o2::aod::track::ITS) ==0;
  Partition<TracksextensionMC> negativeTPConlyTracksMC = o2::aod::track::signed1Pt < fzero && o2::aod::track::tpcNClsFindable > (uint8_t)0 && o2::aod::track::itsChi2NCl < (float_t)0 && o2::aod::track::trackType == (uint8_t)o2::aod::track::TrackTypeEnum::Track; //&& (o2::aod::track::detectorMap & o2::aod::track::TPC) ==o2::aod::track::TPC && (o2::aod::track::detectorMap & o2::aod::track::ITS) ==0;
  Partition<TracksextensionMC> positiveITSonlyTracksMC = o2::aod::track::signed1Pt > fzero&& o2::aod::track::tpcChi2NCl<(float_t)0 && o2::aod::track::itsChi2NCl>(float_t) 0 && o2::aod::track::trackType == (uint8_t)o2::aod::track::TrackTypeEnum::Track;          // && (o2::aod::track::detectorMap & o2::aod::track::TPC) ==0 && (o2::aod::track::detectorMap & o2::aod::track::ITS) ==o2::aod::track::ITS && o2::aod::track::passedITSNCls==true;// && o2::aod::track::itsNCls==7;
  Partition<TracksextensionMC> negativeITSonlyTracksMC = o2::aod::track::signed1Pt < fzero && o2::aod::track::tpcChi2NCl < (float_t)0 && o2::aod::track::itsChi2NCl > (float_t)0 && o2::aod::track::trackType == (uint8_t)o2::aod::track::TrackTypeEnum::Track;      // && (o2::aod::track::detectorMap & o2::aod::track::TPC)  ==0 &&(o2::aod::track::detectorMap & o2::aod::track::ITS) ==o2::aod::track::ITS && o2::aod::track::passedITSNCls==true;// && o2::aod::track::itsNCls==7;

  void init(InitContext const&)
  {

    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisEta{binsEta, "#it{#eta}"};
    const AxisSpec axisPhi{binsPhi, "#it{#varphi}"};
    const AxisSpec axisNclustTPC{binsNclustTPC, "NclustTPC"};
    const AxisSpec axisTPCITSmatching{binsNclustTPC, "NclustTPC"};
    const AxisSpec axisLength{binsLength, "track length (cm)"};
    const AxisSpec axisRxy{binsRxy, "track Rxy (cm)"};
    const AxisSpec axisIsPropagated{binsIsPropagaed, "is track within BP"};
    histos.add("EventProp/histMCcollZ", "MC coll Z (cm); #it{z_{MCcoll}} (cm)", kTH1D, {{100, -20., 20.}});
    histos.add("EventProp/histDatacollZ", "MC coll Z (cm); #it{z_{MCcoll}} (cm)", kTH1D, {{100, -20., 20.}});
    histos.add("EventProp/histPtTrackNegCollID", "pt", kTH1D, {axisPt});

    histos.add("Tracks/Reco/histptAll", "pt", kTH1D, {axisPt});
    histos.add("Tracks/Reco/histpt3DAll", "All tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB0", "FB0 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB1", "FB1 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB2", "FB2 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB3", "FB3 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB4", "FB4 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB5", "FB5 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DITSonly", "ITSonly tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DTPConly", "TPConly tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DTPConlyPtLengthRadiusNclust", "TPConly PtLengthRadiusNclust;#it{p}_{T} (GeV/#it{c});length (cm);Rx,y (cm);nclust;isPropagated", HistType::kTHnF, {axisPt, axisLength, axisRxy, axisNclustTPC, axisIsPropagated});
    histos.add("Tracks/Reco/histpt3DFB0PtLengthRadiusNclust", "TPConly PtLengthRadiusNclust;#it{p}_{T} (GeV/#it{c});length (cm);Rx,y (cm);nclust;isPropagated", HistType::kTHnF, {axisPt, axisLength, axisRxy, axisNclustTPC, axisIsPropagated});

    histos.add("Tracks/MCgen/histMCgenpt", "pt", kTH1D, {axisPt});
    histos.add("Tracks/MCgen/histMCgen3dPhysPrimary", "MC Phys. Prim.;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/MCgen/histMCgen3dChargedProdRad1to15cm", "MC Prod Rad_xy 1 to 15 cm;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/MCgen/histMCgen3dChargedProdRad1mumto5mm", "MC Prod Rad_xy 1#mum to 5 mm ;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/MCgen/histMCgen3dChargedfromHFdecay", "MC Phys. Prim from HF decay ;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/RecoMCPhysPrimCollMatch/histpt", "pt;#it{p}_{T}^{MC} (GeV/#it{c})", kTH1D, {axisPt});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB0", "FB0;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB1", "FB1;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB2", "FB2;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB3", "FB3;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB4", "FB4;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB5", "FB5;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptITSonly", "ITSonly;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptTPConly", "TPConly;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCFB0", "FB0;#it{p}_{T}^{MC} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCFB1", "FB1;#it{p}_{T}^{MC} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCFB2", "FB2;#it{p}_{T}^{MC} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCFB3", "FB3;#it{p}_{T}^{MC} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCFB4", "FB4;#it{p}_{T}^{MC} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCFB5", "FB5;#it{p}_{T}^{MC} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCITSonly", "ITSonly;#it{p}_{T}^{MC} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCTPConly", "TPConly;#it{p}_{T}^{MC} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptMCTPConlyWithClusters", "TPConlyWithClusters;#it{p}_{T}^{gen} (GeV/#it{c});#it{#eta};#it{#varphi};NclustTPC", HistType::kTHnF, {axisPt, axisEta, axisPhi, axisNclustTPC});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptTPConlyWithClusters", "TPConlyWithClusters;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi};NclustTPC", HistType::kTHnF, {axisPt, axisEta, axisPhi, axisNclustTPC});

    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB0", "FB0;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB1", "FB1;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB2", "FB2;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB3", "FB3;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB4", "FB4;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB5", "FB5;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptITSonly", "ITSonly;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptTPConly", "TPConly;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB0", "FB0;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB1", "FB1;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB2", "FB2;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB3", "FB3;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB4", "FB4;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB5", "FB5;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptITSonly", "ITSonly;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptTPConly", "TPConly;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB0", "FB0;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB1", "FB1;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB2", "FB2;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB3", "FB3;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB4", "FB4;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB5", "FB5;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptITSonly", "ITSonly;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptTPConly", "TPConly;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/Reco/histNclustTPC", "N clusters TPC;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi};NclustTPCl;TPCITSmatching", HistType::kTHnF, {axisPt, axisEta, axisPhi, axisNclustTPC, binsTPCITSmatching});
    histos.add("Tracks/RecoMCVariablesPrimary/histNclustTPC", "N clusters TPC;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi};NclustTPC;TPCITSmatching", HistType::kTHnF, {axisPt, axisEta, axisPhi, axisNclustTPC, binsTPCITSmatching});
  }

  void FillDataHistoStd(const Tracksextension::iterator& track)
  {
    if (std::abs(track.eta()) < 0.9) {
      histos.fill(HIST("Tracks/Reco/histptAll"), track.pt());
    }
    histos.fill(HIST("Tracks/Reco/histpt3DAll"), track.pt(), track.eta(), track.phi());
    Int_t hasITS = 0;
    if (track.itsNCls() > 0)
      hasITS++;
    if (track.itsNClsInnerBarrel() > 0)
      hasITS++;
    histos.fill(HIST("Tracks/Reco/histNclustTPC"), track.pt(), track.eta(), track.phi(), track.tpcNClsFound(), hasITS);
    if (track.isGlobalTrack()) {
      histos.fill(HIST("Tracks/Reco/histpt3DFB0"), track.pt(), track.eta(), track.phi());
      histos.fill(HIST("Tracks/Reco/histpt3DFB0PtLengthRadiusNclust"), track.pt(), track.length(), std::sqrt(track.x() * track.x() + track.y() * track.y()), track.tpcNClsFound(), track.isWithinBeamPipe() ? 1 : 0);
    }
    if (track.trackCutFlagFb1())
      histos.fill(HIST("Tracks/Reco/histpt3DFB1"), track.pt(), track.eta(), track.phi());
    if (track.trackCutFlagFb2())
      histos.fill(HIST("Tracks/Reco/histpt3DFB2"), track.pt(), track.eta(), track.phi());
    if (track.trackCutFlagFb3())
      histos.fill(HIST("Tracks/Reco/histpt3DFB3"), track.pt(), track.eta(), track.phi());
    if (track.trackCutFlagFb4())
      histos.fill(HIST("Tracks/Reco/histpt3DFB4"), track.pt(), track.eta(), track.phi());
    if (track.trackCutFlagFb5())
      histos.fill(HIST("Tracks/Reco/histpt3DFB5"), track.pt(), track.eta(), track.phi());
    if (track.itsChi2NCl() > 0. && track.tpcChi2NCl() < 0.)
      histos.fill(HIST("Tracks/Reco/histpt3DITSonly"), track.pt(), track.eta(), track.phi());
    if (track.itsChi2NCl() < 0. && track.tpcChi2NCl() > 0.)
      histos.fill(HIST("Tracks/Reco/histpt3DTPConly"), track.pt(), track.eta(), track.phi());
    if (std::abs(track.eta()) < 0.8) {
      histos.fill(HIST("Tracks/Reco/histpt3DTPConlyPtLengthRadiusNclust"), track.pt(), track.length(), std::sqrt(track.x() * track.x() + track.y() * track.y()), track.tpcNClsFound(), track.isWithinBeamPipe() ? 1 : 0);
      // Printf("track length %f, radius %f, isPropagated %d, tracktype %d",track.length(),std::sqrt(track.x()*track.x()+track.y()*track.y()),track.isWithinBeamPipe()?1:0,track.trackType());
    }
  }

  void processData(Tracksextension const& tracks)
  {
    for (auto& track : tracks) {
      FillDataHistoStd(track);
    }
  }
  PROCESS_SWITCH(CheckFilterBit, processData, "process data", true);

  void processDataCombineTracks(o2::aod::Collision const& collision, Tracksextension const& tracks)
  {
    auto positiveITSonlyTracksThisColl = positiveITSonlyTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negativeITSonlyTracksThisColl = negativeITSonlyTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto positiveTPConlyTracksThisColl = positiveTPConlyTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negativeTPConlyTracksThisColl = negativeTPConlyTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    processPair<false>(collision, positiveTPConlyTracksThisColl, positiveITSonlyTracksThisColl, negativeTPConlyTracksThisColl, negativeITSonlyTracksThisColl);
  }
  PROCESS_SWITCH(CheckFilterBit, processDataCombineTracks, "process data combined tracks", false);

  void processMCCombineTracks(soa::Join<aod::Collisions, o2::aod::McCollisionLabels>::iterator const& collision, TracksextensionMC const& tracks, aod::McParticles const& mcParticles)
  {

    if (std::abs(collision.posZ()) > 10.) {
      return;
    }

    auto positiveITSonlyTracksThisColl = positiveITSonlyTracksMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negativeITSonlyTracksThisColl = negativeITSonlyTracksMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto positiveTPConlyTracksThisColl = positiveTPConlyTracksMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negativeTPConlyTracksThisColl = negativeTPConlyTracksMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    processPair<true>(collision, positiveTPConlyTracksThisColl, positiveITSonlyTracksThisColl, negativeTPConlyTracksThisColl, negativeITSonlyTracksThisColl);
  }
  PROCESS_SWITCH(CheckFilterBit, processMCCombineTracks, "process data combined tracks", true);

  template <bool IS_MC, typename Tcollisions, typename TpTPC, typename TpITS, typename TnTPC, typename TnITS>
  void processPair(Tcollisions const& collision, const TpTPC& positiveTPConlyTracksThisColl, const TpITS& positiveITSonlyTracksThisColl, const TnTPC& negativeTPConlyTracksThisColl, const TnITS& negativeITSonlyTracksThisColl)
  {
    float issamemcpt = -9999.;
    // run first over global tracks
    /*  for (auto& track : tracks) {
      if (track.isGlobalTrack()) {
  track.pt();
      }
      else{
  track.pt();
      }
    }
    */
    for (auto& [track0, track1] : combinations(soa::CombinationsFullIndexPolicy(positiveITSonlyTracksThisColl, positiveTPConlyTracksThisColl))) {
      issamemcpt = -9999.;
      if constexpr (IS_MC) {

        if (track0.mcParticleId() == track1.mcParticleId()) {
          if (track0.has_mcParticle()) {
            /// the track is not fake
            auto mcparticle = track0.mcParticle();
            auto mcCollID_recoColl = collision.mcCollisionId();
            auto mcCollID_particle = mcparticle.mcCollisionId();
            bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);

            int partpdg = std::abs(mcparticle.pdgCode());
            if (indexMatchOK && mcparticle.isPhysicalPrimary() && (partpdg == 211 || partpdg == 321 || partpdg == 2212 || partpdg == 11 || partpdg == 13)) {
              issamemcpt = mcparticle.pt();
            } else {
              issamemcpt = -mcparticle.pt();
            }
          }
        }
      }

      trackPairForEffTablePP(track1.pt(), track1.eta(), track1.phi(), track0.pt(), track0.eta(), track0.phi(), track0.itsNCls(), track1.tpcNClsFound(), issamemcpt);
    }
    for (auto& [track0, track1] : combinations(soa::CombinationsFullIndexPolicy(negativeITSonlyTracksThisColl, negativeTPConlyTracksThisColl))) {

      issamemcpt = -9999.;
      if constexpr (IS_MC) {
        if (track0.mcParticleId() == track1.mcParticleId()) {
          if (track0.has_mcParticle()) {
            /// the track is not fake
            auto mcparticle = track0.mcParticle();
            auto mcCollID_recoColl = collision.mcCollisionId();
            auto mcCollID_particle = mcparticle.mcCollisionId();
            bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);

            int partpdg = std::abs(mcparticle.pdgCode());
            if (indexMatchOK && mcparticle.isPhysicalPrimary() && (partpdg == 211 || partpdg == 321 || partpdg == 2212 || partpdg == 11 || partpdg == 13)) {
              issamemcpt = mcparticle.pt();
            } else {
              issamemcpt = -mcparticle.pt();
            }
          }
        }
      }

      trackPairForEffTableNN(track1.pt(), track1.eta(), track1.phi(), track0.pt(), track0.eta(), track0.phi(), track0.itsNCls(), track1.tpcNClsFound(), issamemcpt);
    }
    for (auto& [track0, track1] : combinations(soa::CombinationsFullIndexPolicy(negativeITSonlyTracksThisColl, positiveTPConlyTracksThisColl))) {
      issamemcpt = -9999.;
      if constexpr (IS_MC) {
        if (track0.mcParticleId() == track1.mcParticleId()) {
          if (track0.has_mcParticle()) {
            /// the track is not fake
            auto mcparticle = track0.mcParticle();
            auto mcCollID_recoColl = collision.mcCollisionId();
            auto mcCollID_particle = mcparticle.mcCollisionId();
            bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);

            int partpdg = std::abs(mcparticle.pdgCode());
            if (indexMatchOK && mcparticle.isPhysicalPrimary() && (partpdg == 211 || partpdg == 321 || partpdg == 2212 || partpdg == 11 || partpdg == 13)) {
              issamemcpt = mcparticle.pt();
            } else {
              issamemcpt = -mcparticle.pt();
            }
          }
        }
      }

      trackPairForEffTableNP(track1.pt(), track1.eta(), track1.phi(), track0.pt(), track0.eta(), track0.phi(), track0.itsNCls(), track1.tpcNClsFound(), issamemcpt);
    }
    for (auto& [track0, track1] : combinations(soa::CombinationsFullIndexPolicy(positiveITSonlyTracksThisColl, negativeTPConlyTracksThisColl))) {
      issamemcpt = -9999.;
      if constexpr (IS_MC) {
        if (track0.mcParticleId() == track1.mcParticleId()) {
          if (track0.has_mcParticle()) {
            /// the track is not fake
            auto mcparticle = track0.mcParticle();
            auto mcCollID_recoColl = collision.mcCollisionId();
            auto mcCollID_particle = mcparticle.mcCollisionId();
            bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);

            int partpdg = std::abs(mcparticle.pdgCode());
            if (indexMatchOK && mcparticle.isPhysicalPrimary() && (partpdg == 211 || partpdg == 321 || partpdg == 2212 || partpdg == 11 || partpdg == 13)) {
              issamemcpt = mcparticle.pt();
            } else {
              issamemcpt = -mcparticle.pt();
            }
          }
        }
      }
      trackPairForEffTablePN(track1.pt(), track1.eta(), track1.phi(), track0.pt(), track0.eta(), track0.phi(), track0.itsNCls(), track1.tpcNClsFound(), issamemcpt);
    }
  }

  template <typename T>
  int isFromStrangeDecay(const T& particlesMC, const typename T::iterator& particle)
  {

    std::vector<std::vector<int64_t>> arrayIds{};
    std::vector<int64_t> initVec{particle.globalIndex()};
    arrayIds.push_back(initVec);

    int stage = 0;
    int strangeness = 0;
    while (stage < 4 && arrayIds[stage].size() > 0 && strangeness == 0) {
      std::vector<int64_t> arrayIdsStage{};
      for (auto& iPart : arrayIds[stage]) {
        auto particleMother = particlesMC.rawIteratorAt(iPart - particlesMC.offset());
        if (particleMother.has_mothers()) {
          for (auto iMother = particleMother.mothersIds().front(); iMother <= particleMother.mothersIds().back(); ++iMother) {
            if (std::find(arrayIdsStage.begin(), arrayIdsStage.end(), iMother) != arrayIdsStage.end()) {
              continue;
            }
            auto mother = particlesMC.rawIteratorAt(iMother - particlesMC.offset());
            auto motherPDG = std::abs(mother.pdgCode());
            if (motherPDG / 1000 == 3) { // strange baryon
              strangeness = 1;
              break;
            } else if (motherPDG / 1000 == 0 && motherPDG / 100 == 3) {
              strangeness = 1;
              break;
            }
            arrayIdsStage.push_back(iMother);
          }
        }
        arrayIds.push_back(arrayIdsStage);
      }
      stage++;
    }
    return strangeness;
  }

  void processRecoMC(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, TracksextensionMC const& tracks, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  { // this will loop over data (PV) collisions

    histos.fill(HIST("EventProp/histDatacollZ"), collision.posZ());
    if (std::abs(collision.posZ()) > 10.)
      return;
    ncollisionCounter++;
    for (auto& track : tracks) {

      if (track.collisionId() < 0)
        histos.fill(HIST("EventProp/histPtTrackNegCollID"), track.pt());
      bool has_MCparticle = track.has_mcParticle();
      if (has_MCparticle) {
        /// the track is not fake
        auto mcparticle = track.mcParticle();
        // auto collReco = track.collision_as<CollisionTableMC>();
        auto collMC = mcparticle.mcCollision();
        auto mcCollID_recoColl = collision.mcCollisionId();
        auto mcCollID_particle = mcparticle.mcCollisionId();
        bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);
        // double pvZdiff = collision.posZ() - collMC.posZ();
        if (indexMatchOK) {
          double prodRadius2 = (mcparticle.vx() - collMC.posX()) * (mcparticle.vx() - collMC.posX()) + (mcparticle.vy() - collMC.posY()) * (mcparticle.vy() - collMC.posY());
          int partpdg = std::abs(mcparticle.pdgCode());
          if (partpdg == 211 || partpdg == 321 || partpdg == 2212 || partpdg == 11 || partpdg == 13) {
            int isFromStrange = isFromStrangeDecay(mcParticles, mcparticle);
            if (mcparticle.isPhysicalPrimary()) {
              int isHF = RecoDecay::getCharmHadronOrigin(mcParticles, mcparticle, false);
              if (std::abs(track.eta()) < 0.9) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histpt"), mcparticle.pt()); // note: one needs to avoid double counting of tracks reco both in TPC and ITS but not matched
              }
              Int_t hasITS = 0;
              if (track.itsNCls() > 0.)
                hasITS++;
              if (track.itsNClsInnerBarrel() > 0.)
                hasITS++;
              histos.fill(HIST("Tracks/RecoMCVariablesPrimary/histNclustTPC"), track.pt(), track.eta(), track.phi(), track.tpcNClsFound(), hasITS);
              if (track.isGlobalTrack()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB0"), track.pt(), track.eta(), track.phi());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCFB0"), mcparticle.pt(), track.eta(), track.phi());
              }
              if (track.itsChi2NCl() > 0. && track.tpcChi2NCl() < 0.) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptITSonly"), track.pt(), track.eta(), track.phi());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCITSonly"), mcparticle.pt(), track.eta(), track.phi());
              } else if (track.itsChi2NCl() < 0. && track.tpcChi2NCl() > 0.) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptTPConly"), track.pt(), track.eta(), track.phi());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCTPConly"), mcparticle.pt(), track.eta(), track.phi());

                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCTPConlyWithClusters"), mcparticle.pt(), track.eta(), track.phi(), track.tpcNClsFound());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptTPConlyWithClusters"), track.pt(), track.eta(), track.phi(), track.tpcNClsFound());
              }
              if (track.trackCutFlagFb1()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB1"), track.pt(), track.eta(), track.phi());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCFB1"), mcparticle.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb2()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB2"), track.pt(), track.eta(), track.phi());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCFB2"), mcparticle.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb3()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB3"), track.pt(), track.eta(), track.phi());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCFB3"), mcparticle.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb4()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB4"), track.pt(), track.eta(), track.phi());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCFB4"), mcparticle.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb5()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB5"), track.pt(), track.eta(), track.phi());
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptMCFB5"), mcparticle.pt(), track.eta(), track.phi());
              }
              if (isHF == RecoDecay::OriginType::Prompt || isHF == RecoDecay::OriginType::NonPrompt) {
                if (track.isGlobalTrack()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB0"), track.pt(), track.eta(), track.phi());
                }
                if (track.itsChi2NCl() > 0. && track.tpcChi2NCl() < 0.) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptITSonly"), track.pt(), track.eta(), track.phi());
                } else if (track.itsChi2NCl() < 0. && track.tpcChi2NCl() > 0.) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptTPConly"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb1()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB1"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb2()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB2"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb3()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB3"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb4()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB4"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb5()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB5"), track.pt(), track.eta(), track.phi());
                }
              }
            } else if (prodRadius2 > 1. && prodRadius2 < 225. && isFromStrange) {
              if (track.isGlobalTrack())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB0"), track.pt(), track.eta(), track.phi());
              if (track.itsChi2NCl() > 0. && track.tpcChi2NCl() < 0.) {
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptITSonly"), track.pt(), track.eta(), track.phi());
              } else if (track.itsChi2NCl() < 0. && track.tpcChi2NCl() > 0.) {
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptTPConly"), track.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb1())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB1"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb2())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB2"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb3())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB3"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb4())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB4"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb5())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB5"), track.pt(), track.eta(), track.phi());
            }
            if (prodRadius2 > 1.e-8 && prodRadius2 < 0.25) {
              if (track.isGlobalTrack())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB0"), track.pt(), track.eta(), track.phi());
              if (track.itsChi2NCl() > 0. && track.tpcChi2NCl() < 0.) {
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptITSonly"), track.pt(), track.eta(), track.phi());
              } else if (track.itsChi2NCl() < 0. && track.tpcChi2NCl() > 0.) {
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptTPConly"), track.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb1())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB1"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb2())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB2"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb3())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB3"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb4())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB4"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb5())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB5"), track.pt(), track.eta(), track.phi());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(CheckFilterBit, processRecoMC, "process mc", true);

  void processMC(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {

    histos.fill(HIST("EventProp/histMCcollZ"), mcCollision.posZ());
    if (std::abs(mcCollision.posZ()) > 10.)
      return;
    ncollisionCounter++;
    for (auto& mcpart : mcParticles) {
      // if(!mcpart.producedByGenerator())continue;
      double prodRadius2 = (mcpart.vx() - mcCollision.posX()) * (mcpart.vx() - mcCollision.posX()) + (mcpart.vy() - mcCollision.posY()) * (mcpart.vy() - mcCollision.posY());
      if (std::abs(mcpart.pdgCode()) == 211 || std::abs(mcpart.pdgCode()) == 321 || std::abs(mcpart.pdgCode()) == 2212 || std::abs(mcpart.pdgCode()) == 11 || std::abs(mcpart.pdgCode()) == 13) {
        int isHF = RecoDecay::getCharmHadronOrigin(mcParticles, mcpart, false);
        int isFromStrange = isFromStrangeDecay(mcParticles, mcpart);
        if (mcpart.isPhysicalPrimary()) {
          if (std::abs(mcpart.eta()) < 0.9) {
            histos.fill(HIST("Tracks/MCgen/histMCgenpt"), mcpart.pt());
          }
          histos.fill(HIST("Tracks/MCgen/histMCgen3dPhysPrimary"), mcpart.pt(), mcpart.eta(), mcpart.phi());
          if (isHF == RecoDecay::OriginType::Prompt || isHF == RecoDecay::OriginType::NonPrompt)
            histos.fill(HIST("Tracks/MCgen/histMCgen3dChargedfromHFdecay"), mcpart.pt(), mcpart.eta(), mcpart.phi());
        } else if (prodRadius2 > 1. && prodRadius2 < 225. && isFromStrange) {
          histos.fill(HIST("Tracks/MCgen/histMCgen3dChargedProdRad1to15cm"), mcpart.pt(), mcpart.eta(), mcpart.phi());
        }
        if (prodRadius2 > 1.e-8 && prodRadius2 < 0.25) {
          histos.fill(HIST("Tracks/MCgen/histMCgen3dChargedProdRad1mumto5mm"), mcpart.pt(), mcpart.eta(), mcpart.phi());
        }
      }
    }
  }
  PROCESS_SWITCH(CheckFilterBit, processMC, "process mc", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CheckFilterBit>(cfgc)

  };
}
