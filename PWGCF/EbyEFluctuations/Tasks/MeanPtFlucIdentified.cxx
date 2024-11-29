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

/// \file MeanPtFlucIdentified.cxx
/// \brief Calculate EbyE <pt> fluctuations with cumulant method.
///        For charged particles and identified particles.
///        For RUN-3
///
/// \author Tanu Gahlaut <tanu.gahlaut@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/RecoDecay.h"

#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace std;

struct meanPtFlucId {
  Configurable<int> nPtBins{"nPtBins", 300, ""};
  Configurable<int> nPartBins{"nPartBins", 250, ""};
  Configurable<int> nCentBins{"nCentBins", 101, ""};
  Configurable<int> nEtaBins{"nEtaBins", 100, ""};
  Configurable<int> nPhiBins{"nPhiBins", 100, ""};
  Configurable<float> cfgCutPtMax{"cfgCutPtMax", 3.0, "maximum pT"};
  Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.15, "minimum pT"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut"};
  Configurable<float> cfgCutRap{"cfgCutRap", 0.5, "Rapidity Cut"};
  Configurable<float> cfgCutDcaXY{"cfgCutDcaXY", 0.12, "DCAxy cut"};
  Configurable<float> cfgCutDcaZ{"cfgCutDcaZ", 0.3, "DCAz cut"};
  Configurable<float> cfgCutPosZ{"cfgCutPosZ", 10.0, "cut for vertex Z"};
  Configurable<float> cfgGammaCut{"cfgGammaCut", 0.003, "Gamma inv Mass Cut for electron-positron rejection"};
  Configurable<float> cfgCutNSig2{"cfgCutNSig2", 2.0, "nSigma cut (2)"};
  Configurable<float> cfgCutNSig3{"cfgCutNSig3", 3.0, "nSigma cut (3)"};
  Configurable<float> cfgCutPiPtMin{"cfgCutPiPtMin", 0.2, "Minimum pion p_{T} cut"};
  Configurable<float> cfgCutKaPtMin{"cfgCutKaPtMin", 0.3, "Minimum kaon p_{T} cut"};
  Configurable<float> cfgCutPrPtMin{"cfgCutPrPtMin", 0.5, "Minimum proton p_{T} cut"};
  Configurable<float> cfgCutPiThrsldP{"cfgCutPiThrsldP", 0.6, "Threshold p cut pion"};
  Configurable<float> cfgCutKaThrsldP{"cfgCutKaThrsldP", 0.6, "Threshold p cut kaon"};
  Configurable<float> cfgCutPrThrsldP{"cfgCutPrThrsldP", 1.0, "Threshold p cut proton "};
  Configurable<float> cfgCutPiP1{"cfgCutPiP1", 0.5, "pion p cut-1"};
  Configurable<float> cfgCutPiP2{"cfgCutPiP2", 0.6, "pion p cut-2"};
  Configurable<float> cfgCutKaP1{"cfgCutKaP1", 0.4, "kaon p cut-1"};
  Configurable<float> cfgCutKaP2{"cfgCutKaP2", 0.6, "kaon p cut-2"};
  Configurable<float> cfgCutKaP3{"cfgCutKaP3", 1.2, "kaon p cut-3"};
  Configurable<float> cfgCutPrP1{"cfgCutPrP1", 0.9, "proton p cut-1"};
  Configurable<float> cfgCutPrP2{"cfgCutPrP2", 1.0, "proton p cut-2"};
  Configurable<float> cfgMcTpcShiftEl{"cfgMcTpcShiftEl", 0., "Electron Shift in TPC (MC data) "};
  Configurable<float> cfgMcTpcShiftPi{"cfgMcTpcShiftPi", 0., "Pion Shift in TPC (MC data) "};
  Configurable<float> cfgMcTpcShiftKa{"cfgMcTpcShiftKa", 0., "Kaon Shift in TPC (MC data) "};
  Configurable<float> cfgMcTpcShiftPr{"cfgMcTpcShiftPr", 0., "Proton Shift in TPC (MC data) "};
  Configurable<float> cfgMcTofShiftPi{"cfgMcTofShiftPi", 0., "Pion Shift in TOF (MC data) "};
  Configurable<float> cfgMcTofShiftKa{"cfgMcTofShiftKa", 0., "Kaon Shift in TOF (MC data) "};
  Configurable<float> cfgMcTofShiftPr{"cfgMcTofShiftPr", 0., "Proton Shift in TOF (MC data) "};
  Configurable<bool> cfgMCReco{"cfgMCReco", false, ""};
  Configurable<bool> cfgMCTruth{"cfgMCTruth", false, ""};
  Configurable<bool> cfgPosZ{"cfgPosZ", true, "Position Z"};
  Configurable<bool> cfgSel8{"cfgSel8", true, "Sel8 trigger"};
  Configurable<bool> cfgEvSel1{"cfgEvSel1", true, "kNoSameBunchPileup"};
  Configurable<bool> cfgEvSel2{"cfgEvSel2", true, "kIsGoodZvtxFT0vsPV"};
  Configurable<bool> cfgEvSel3{"cfgEvSel3", true, "kIsVertexITSTPC"};
  Configurable<bool> cfgInvMass{"cfgInvMass", true, "electron Inv Mass cut selection"};
  Configurable<bool> cfgSelOR{"cfgSelOR", true, "Low OR High momentum "};
  Configurable<bool> cfgSelAND{"cfgSelAND", false, "Low AND High momentum"};
  Configurable<bool> cfgSelLow{"cfgSelLow", true, "PID selection cut for Low momentum"};
  Configurable<bool> cfgSelHigh{"cfgSelHigh", true, "PID selection cut for High momentum"};
  ConfigurableAxis multTPCBins{"multTPCBins", {150, 0, 150}, "TPC Multiplicity bins"};
  ConfigurableAxis multFT0MBins{"multFT0MBins", {400, 0, 4000}, "Forward Multiplicity bins"};
  ConfigurableAxis multFT0MMCBins{"multFT0MMCBins", {250, 0, 250}, "Forward Multiplicity bins"};
  ConfigurableAxis dcaXYBins{"dcaXYBins", {100, -0.15, 0.15}, "dcaXY bins"};
  ConfigurableAxis dcaZBins{"dcaZBins", {100, -1.2, 1.2}, "dcaZ bins"};
  ConfigurableAxis QnBins{"QnBins", {1000, 0., 100.}, "nth moments bins"};
  ConfigurableAxis TpNBins{"TpNBins", {300, 0., 3000.}, ""};
  ConfigurableAxis TpDBins{"TpDBins", {100, 0., 2000.}, ""};

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                                aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                                aod::pidTOFbeta, aod::pidTOFmass>;
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Cs>;
  using MyMCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Cs, aod::McCollisionLabels>;
  using MyMCTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                               aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                               aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                               aod::pidTOFbeta, aod::pidTOFmass, aod::McTrackLabels>;

  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    const AxisSpec axisEvents{6, 0, 6, "Counts"};
    const AxisSpec axisEta{nEtaBins, -1., +1., "#eta"};
    const AxisSpec axisPhi{nPhiBins, 0., +7., "#phi (rad)"};
    const AxisSpec axisY{nEtaBins, -1., +1., "y"};
    const AxisSpec axisPt{nPtBins, 0., 3., "p_{T} (GeV/c)"};
    const AxisSpec axisP{nPtBins, 0., 3., "p (GeV/c)"};
    const AxisSpec axisInnerParam{nPtBins, 0., 3., "p_{InnerParam } (GeV/c)"};
    const AxisSpec axisPart{nPartBins, 0., 18., " "};
    const AxisSpec axisQn{QnBins, ""};
    const AxisSpec axisTpN{TpNBins, "(Q_{1}^{2} - Q_{2})"};
    const AxisSpec axisTpD{TpDBins, "N_{pairs}"};
    const AxisSpec axisDeno{100, 1., 2.0, "#frac{1}{#sqrt{1 - #frac{1}{N}}}"};
    const AxisSpec axisMeanPt{100, 0., 3., "M(p_{T}) (GeV/c)"};
    const AxisSpec axisMult{100, 0, 100, "N_{ch}"};
    const AxisSpec axisMultTPC{multTPCBins, "N_{TPC} "};
    const AxisSpec axisMultFT0M{multFT0MBins, "N_{FT0M}"};
    const AxisSpec axisMultFT0MMC{multFT0MMCBins, "N_{FT0M}"};
    const AxisSpec axisCentFT0C{nCentBins, 0, 101, "FT0C (%)"};
    const AxisSpec axisVtxZ{80, -20., 20., "V_{Z} (cm)"};
    const AxisSpec axisDCAz{dcaZBins, "DCA_{Z} (cm)"};
    const AxisSpec axisDCAxy{dcaXYBins, "DCA_{XY} (cm)"};
    const AxisSpec axisTPCNsigma{500, -5., 5., "n #sigma_{TPC}"};
    const AxisSpec axisTOFNsigma{500, -5., 5., "n #sigma_{TOF}"};
    const AxisSpec axisTPCSignal{100, 20., 500., "#frac{dE}{dx}"};
    const AxisSpec axisTOFSignal{200, 0.2, 1.2, "TOF #beta"};
    const AxisSpec axisChi2{40, 0., 40., "Chi2"};
    const AxisSpec axisCrossedTPC{300, 0, 300, "Crossed TPC"};
    const AxisSpec axisM2{100, 0., 1.4, "#it{m}^{2} (GeV/#it{c}^{2})^{2}"};
    const AxisSpec axisMass{1000, 0., 0.1, "M_{inv} (GeV/#it{c}^2)"};

    HistogramConfigSpec QnHist({HistType::kTHnSparseD, {axisMultTPC, axisQn, axisMultFT0M}});
    HistogramConfigSpec PartHist({HistType::kTHnSparseD, {axisMultTPC, axisPart, axisMultFT0M}});
    HistogramConfigSpec DenoHist({HistType::kTHnSparseD, {axisMultTPC, axisDeno, axisMultFT0M}});
    HistogramConfigSpec QnMCHist({HistType::kTHnSparseD, {axisMultTPC, axisQn, axisMultFT0MMC}});
    HistogramConfigSpec PartMCHist({HistType::kTHnSparseD, {axisMultTPC, axisPart, axisMultFT0MMC}});
    HistogramConfigSpec DenoMCHist({HistType::kTHnSparseD, {axisMultTPC, axisDeno, axisMultFT0MMC}});
    HistogramConfigSpec TOFnSigmaHist({HistType::kTH2D, {axisP, axisTOFNsigma}});
    HistogramConfigSpec TOFSignalHist({HistType::kTH2D, {axisP, axisTOFSignal}});
    HistogramConfigSpec TPCnSigmaHist({HistType::kTH2D, {axisP, axisTPCNsigma}});
    HistogramConfigSpec TPCSignalHist({HistType::kTH2D, {axisP, axisTPCSignal}});
    HistogramConfigSpec TPCTOFHist({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});
    HistogramConfigSpec PvsM2Hist({HistType::kTH2D, {axisM2, axisP}});

    HistogramConfigSpec TOFnSigmaHist1({HistType::kTH2D, {axisInnerParam, axisTOFNsigma}});
    HistogramConfigSpec TOFSignalHist1({HistType::kTH2D, {axisInnerParam, axisTOFSignal}});
    HistogramConfigSpec TPCnSigmaHist1({HistType::kTH2D, {axisInnerParam, axisTPCNsigma}});
    HistogramConfigSpec TPCSignalHist1({HistType::kTH2D, {axisInnerParam, axisTPCSignal}});
    HistogramConfigSpec TPCTOFHist1({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});
    HistogramConfigSpec PvsM2Hist1({HistType::kTH2D, {axisM2, axisInnerParam}});

    // QA Plots:
    hist.add("QA/before/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("QA/before/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});
    hist.add("QA/before/h_TPCChi2perCluster", "TPC #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/before/h_ITSChi2perCluster", "ITS #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/before/h_crossedTPC", "Crossed TPC", kTH1D, {axisCrossedTPC});
    hist.add("QA/before/h_Pt", "p_{T}", kTH1D, {axisPt});
    hist.add("QA/before/h2_PvsPinner", "p_{InnerParam} vs p", kTH2D, {{axisP}, {axisInnerParam}});
    hist.add("QA/before/h2_PtofvsPinner", "p_{InnerParam} vs p_{TOF}", kTH2D, {{axisP}, {axisInnerParam}});
    hist.add("QA/before/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/before/h_Phi", "#phi ", kTH1D, {axisPhi});
    hist.add("QA/before/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/before/h_DcaZ", "DCA_{Z}", kTH1D, {axisDCAz});
    hist.add("QA/before/h_DcaXY", "DCA_{XY}", kTH1D, {axisDCAxy});
    hist.add("QA/before/h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/before/h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/before/h_NTPC", "N_{TPC}", kTH1D, {axisMultTPC});
    hist.add("QA/before/h_NFT0M", "FT0M Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/before/h_NFT0C", "FT0C Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/before/h_Cent", "FT0C (%)", kTH1D, {axisCentFT0C});
    hist.add("QA/before/h2_NTPC_Cent", "N_{TPC} vs FT0C(%)", kTH2D, {{axisCentFT0C}, {axisMultTPC}});
    hist.add("QA/before/h2_NTPC_NFT0M", "N_{TPC} vs N_{FT0M}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});
    hist.add("QA/before/h2_NTPC_NFT0C", "N_{TPC} vs N_{FT0C}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});

    hist.add("QA/before/h2_TPCSignal", "TPC Signal", TPCSignalHist);
    hist.add("QA/before/h2_TOFSignal", "TOF Signal", TOFSignalHist);
    hist.add("QA/before/h2_pvsm2", "p vs m^{2}", PvsM2Hist);

    hist.add("QA/before/innerParam/h2_TPCSignal", "TPC Signal", TPCSignalHist1);
    hist.add("QA/before/innerParam/h2_TOFSignal", "TOF Signal", TOFSignalHist1);
    hist.add("QA/before/innerParam/h2_pvsm2", "p vs m^{2}", PvsM2Hist1);

    hist.addClone("QA/before/", "QA/after/");

    hist.add("QA/after/p_NTPC_NFT0M", "N_{TPC} vs N_{FT0M} (Profile)", kTProfile, {axisMultFT0M});
    hist.add("QA/after/p_NTPC_NFT0C", "N_{TPC} vs N_{FT0C} (Profile)", kTProfile, {axisMultFT0M});
    hist.add("QA/after/p_NTPC_Cent", "N_{TPC} vs FT0C(%) (Profile)", kTProfile, {axisCentFT0C});
    hist.add("QA/after/h2_NTPC_Nch", "N_{ch} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMult}});
    hist.add("QA/after/h_invMass_gamma", "Inv Mass of #gamma", kTH1D, {axisMass});
    hist.add("QA/after/counts_evSelCuts", "Event selection cuts", kTH1D, {axisEvents});
    hist.add("QA/after/h_vtxZ_sim", "Simulated Vertex Z", kTH1D, {axisVtxZ});
    hist.add("QA/after/h_NSim", "Truth Multiplicity TPC", kTH1D, {axisMultTPC});
    hist.add("QA/after/h2_NTPC_NSim", "Reco vs Truth Multiplicty TPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("QA/after/h2_NchSim_NSim", "Truth Multiplicty Nch vs NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("QA/after/h2_NFT0C_NFT0CSim", "Reco vs Truth Multplicity FT0C", kTH2D, {{axisMultFT0MMC}, {axisMultFT0M}});

    hist.add("QA/Pion/h_Pt", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtPos", "p_{T} (positive) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtNeg", "p_{T} (negative) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtPosTruth", "p_{T} (positive) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtNegTruth", "p_{T} (negative) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_rap", "y ", kTH1D, {axisY});
    hist.add("QA/Pion/h_Eta", "Pseudorapidity ", kTH1D, {axisEta});
    hist.add("QA/Pion/h_Phi", "Azimuthal Distribution ", kTH1D, {axisPhi});
    hist.add("QA/Pion/h2_Pt_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h_DcaZ", "DCA_{z}", kTH1D, {axisDCAz});
    hist.add("QA/Pion/h_DcaXY", "DCA_{xy}", kTH1D, {axisDCAxy});
    hist.add("QA/Pion/h2_DcaZ", "DCA_{z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/Pion/h2_DcaXY", "DCA_{xy}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/Pion/h2_P_Pinner", "p_{TPCinner} vs p", kTH2D, {{axisP}, {axisInnerParam}});
    hist.add("QA/Pion/h2_Pt_Pinner", "p_{TPCinner} vs p_{T}", kTH2D, {{axisPt}, {axisInnerParam}});

    hist.add("QA/Pion/before/h2_TPCNsigma", "n #sigma_{TPC}", TPCnSigmaHist);
    hist.add("QA/Pion/before/h2_TOFNsigma", "n #sigma_{TOF}", TOFnSigmaHist);
    hist.add("QA/Pion/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", TPCTOFHist);
    hist.add("QA/Pion/h2_TPCNsigma", "n #sigma_{TPC}", TPCnSigmaHist);
    hist.add("QA/Pion/h2_TPCNsigma_El", "n #sigma_{TPC, El}", TPCnSigmaHist);
    hist.add("QA/Pion/h2_TOFNsigma_El", "n #sigma_{TOF, El}", TOFnSigmaHist);
    hist.add("QA/Pion/h2_TOFNsigma", "n #sigma_{TOF}", TOFnSigmaHist);
    hist.add("QA/Pion/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", TPCTOFHist);
    hist.add("QA/Pion/h2_TPCSignal", "TPC Signal ", TPCSignalHist);
    hist.add("QA/Pion/h2_TOFSignal", "TOF Signal", TOFSignalHist);
    hist.add("QA/Pion/h2_pvsm2", "p vs m^{2}", PvsM2Hist);

    hist.add("QA/Pion/innerParam/before/h2_TPCNsigma", "n #sigma_{TPC}", TPCnSigmaHist1);
    hist.add("QA/Pion/innerParam/before/h2_TOFNsigma", "n #sigma_{TOF}", TOFnSigmaHist1);
    hist.add("QA/Pion/innerParam/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", TPCTOFHist1);
    hist.add("QA/Pion/innerParam/h2_TPCNsigma", "n #sigma_{TPC}", TPCnSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TPCNsigma_El", "n #sigma_{TPC, El}", TPCnSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TOFNsigma_El", "n #sigma_{TOF, El}", TOFnSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TOFNsigma", "n #sigma_{TOF}", TOFnSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", TPCTOFHist1);
    hist.add("QA/Pion/innerParam/h2_TPCSignal", "TPC Signal ", TPCSignalHist1);
    hist.add("QA/Pion/innerParam/h2_TOFSignal", "TOF Signal", TOFSignalHist1);
    hist.add("QA/Pion/innerParam/h2_pvsm2", "p vs m^{2}", PvsM2Hist1);

    hist.addClone("QA/Pion/", "QA/Kaon/");
    hist.addClone("QA/Pion/", "QA/Proton/");

    // Analysis Plots:
    hist.add("Analysis/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Analysis/Charged/h_Q1", "Q1", QnHist);
    hist.add("Analysis/Charged/h_Q2", "Q2", QnHist);
    hist.add("Analysis/Charged/h_Q3", "Q3", QnHist);
    hist.add("Analysis/Charged/h_Q4", "Q4", QnHist);
    hist.add("Analysis/Charged/h_mean_pT", " <p_{T}> ", kTH1D, {axisMeanPt});
    hist.add("Analysis/Charged/p_mean_pT_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/p_CheckNCH", " 1/denominator vs N_{TPC} ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/h_CheckNCH", " 1/denominator vs N_{TPC} ", DenoHist);
    hist.add("Analysis/Charged/h_Q1_var", "Q1 vs N_{TPC}", QnHist);
    hist.add("Analysis/Charged/h_N_var", "N vs N_{TPC}", kTHnSparseD, {axisMultTPC, axisMult, axisMultFT0M});
    hist.add("Analysis/Charged/h_twopart_nume_Mult_var", "twopart numerator", kTHnSparseD, {axisMultTPC, axisTpN, axisMultFT0M});
    hist.add("Analysis/Charged/h_twopart_deno_Mult_var", "twopart denominator", kTHnSparseD, {axisMultTPC, axisTpD, axisMultFT0M});
    hist.add("Analysis/Charged/h_mean_pT_Mult_var", " <p_{T}> vs N_{TPC} ", PartHist);
    hist.add("Analysis/Charged/h_mean_pT_Mult_skew", " <p_{T}> vs N_{TPC} ", PartHist);
    hist.add("Analysis/Charged/h_mean_pT_Mult_kurto", " <p_{T}> vs N_{TPC} ", PartHist);
    hist.add("Analysis/Charged/h_twopart_Mult_var", "Twopart vs N_{TPC} ", PartHist);
    hist.add("Analysis/Charged/h_twopart_Mult_skew", "Twopart vs N_{TPC} ", PartHist);
    hist.add("Analysis/Charged/h_twopart_Mult_kurto", "Twopart vs N_{TPC} ", PartHist);
    hist.add("Analysis/Charged/h_threepart_Mult_skew", "Threepart vs N_{TPC} ", PartHist);
    hist.add("Analysis/Charged/h_threepart_Mult_kurto", "Threepart vs N_{TPC} ", PartHist);
    hist.add("Analysis/Charged/h_fourpart_Mult_kurto", "Fourpart vs N_{TPC} ", PartHist);

    hist.addClone("Analysis/Charged/", "Analysis/Pion/");
    hist.addClone("Analysis/Charged/", "Analysis/Kaon/");
    hist.addClone("Analysis/Charged/", "Analysis/Proton/");

    // MC Generated
    hist.add("Gen/Counts", "Counts", kTH1D, {axisEvents});
    hist.add("Gen/vtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/NTPC", "Mid rapidity Multiplicity", kTH1D, {axisMultTPC});
    hist.add("Gen/NFT0C", "Forward Multiplicity", kTH1D, {axisMultFT0MMC});
    hist.add("Gen/h2_NTPC_NFT0C", "N_{TPC} vs N_{FT0C}", kTH2D, {{axisMultFT0MMC}, {axisMultTPC}});
    hist.add("Gen/Charged/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("Gen/Charged/h_PtPosTruth", "p_{T} (Positive)", kTH1D, {axisPt});
    hist.add("Gen/Charged/h_PtNegTruth", "p_{T} (negative)", kTH1D, {axisPt});

    hist.add("Gen/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Gen/Charged/h_mean_pT", " <p_{T}> ", kTH1D, {axisMeanPt});

    hist.add("Gen/Charged/h_Q1", "Q1", QnMCHist);
    hist.add("Gen/Charged/h_Q2", "Q2", QnMCHist);
    hist.add("Gen/Charged/h_Q3", "Q3", QnMCHist);
    hist.add("Gen/Charged/h_Q4", "Q4", QnMCHist);
    hist.add("Gen/Charged/h_Q1_var", "Q1 vs N_{TPC}", QnMCHist);
    hist.add("Gen/Charged/h_N_var", "N vs N_{TPC}", kTHnSparseD, {axisMultTPC, axisMult, axisMultFT0MMC});
    hist.add("Gen/Charged/h_twopart_nume_Mult_var", "twopart numerator", kTHnSparseD, {axisMultTPC, axisTpN, axisMultFT0MMC});
    hist.add("Gen/Charged/h_twopart_deno_Mult_var", "twopart denominator", kTHnSparseD, {axisMultTPC, axisTpD, axisMultFT0MMC});

    hist.add("Gen/Charged/p_mean_pT_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
    hist.add("Gen/Charged/p_CheckNCH", " 1/denominator vs N_{TPC} ", kTProfile, {axisMultTPC});
    hist.add("Gen/Charged/h_CheckNCH", " 1/denominator vs N_{TPC} ", DenoMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_var", " <p_{T}> vs N_{TPC} ", PartMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_skew", " <p_{T}> vs N_{TPC} ", PartMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_kurto", " <p_{T}> vs N_{TPC} ", PartMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_var", "Twopart vs N_{TPC} ", PartMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_skew", "Twopart vs N_{TPC} ", PartMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_kurto", "Twopart vs N_{TPC} ", PartMCHist);
    hist.add("Gen/Charged/h_threepart_Mult_skew", "Threepart vs N_{TPC} ", PartMCHist);
    hist.add("Gen/Charged/h_threepart_Mult_kurto", "Threepart vs N_{TPC} ", PartMCHist);
    hist.add("Gen/Charged/h_fourpart_Mult_kurto", "Fourpart vs N_{TPC} ", PartMCHist);

    hist.addClone("Gen/Charged/", "Gen/Pion/");
    hist.addClone("Gen/Charged/", "Gen/Kaon/");
    hist.addClone("Gen/Charged/", "Gen/Proton/");
  }

  enum mode {
    QA_Pion = 0,
    QA_Kaon,
    QA_Proton,
    Analysis_Charged,
    Analysis_Pion,
    Analysis_Kaon,
    Analysis_Proton,
    Gen_Charged,
    Gen_Pion,
    Gen_Kaon,
    Gen_Proton
  };

  static constexpr std::string_view dire[] = {
    "QA/Pion/",
    "QA/Kaon/",
    "QA/Proton/",
    "Analysis/Charged/",
    "Analysis/Pion/",
    "Analysis/Kaon/",
    "Analysis/Proton/",
    "Gen/Charged/",
    "Gen/Pion/",
    "Gen/Kaon/",
    "Gen/Proton/"};

  // Event selection cuts:
  template <typename T>
  bool selRun3Col(T const& col)
  {
    hist.fill(HIST("QA/after/counts_evSelCuts"), 0);

    if (cfgPosZ && std::abs(col.posZ()) > cfgCutPosZ)
      return false;
    hist.fill(HIST("QA/after/counts_evSelCuts"), 1);

    if (cfgSel8 && !col.sel8())
      return false;
    hist.fill(HIST("QA/after/counts_evSelCuts"), 2);

    if (cfgEvSel1 && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return false;
    hist.fill(HIST("QA/after/counts_evSelCuts"), 3);

    if (cfgEvSel2 && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    hist.fill(HIST("QA/after/counts_evSelCuts"), 4);

    if (cfgEvSel3 && !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return false;
    hist.fill(HIST("QA/after/counts_evSelCuts"), 5);

    return true;
  }

  // Track selection cuts:
  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrackWoPtEta())
      return false;

    if (track.pt() < cfgCutPtMin)
      return false;

    if (track.pt() > cfgCutPtMax)
      return false;

    if (track.sign() == 0)
      return false;

    if (std::fabs(track.dcaZ()) > cfgCutDcaZ)
      return false;

    if (std::fabs(track.dcaXY()) > cfgCutDcaXY)
      return false;

    return true;
  }

  // Cuts to reject the tracks
  template <typename T>
  bool rejectTracks(T const& track)
  {
    if (((track.tpcNSigmaEl() - cfgMcTpcShiftEl) > -3. &&
         (track.tpcNSigmaEl() - cfgMcTpcShiftEl) < 5.) &&
        (std::fabs(track.tpcNSigmaPi() - cfgMcTpcShiftPi) > 3 &&
         std::fabs(track.tpcNSigmaKa() - cfgMcTpcShiftKa) > 3 &&
         std::fabs(track.tpcNSigmaPr() - cfgMcTpcShiftPr) > 3)) {
      return true;
    }

    return false;
  }

  template <typename T>
  bool selElectrons(T const& track)
  {
    if (std::fabs(track.tpcNSigmaEl() - cfgMcTpcShiftEl) < cfgCutNSig3) {
      return true;
    }

    return false;
  }

  // PID selction cuts for Low momentum Pions
  template <typename T>
  bool selLowPi(T const& track, double p)
  {
    if (track.pt() >= cfgCutPiPtMin &&
        track.p() <= cfgCutPiThrsldP &&
        ((!track.hasTOF() &&
          ((std::fabs(track.tpcNSigmaPi() - cfgMcTpcShiftPi) < cfgCutNSig3 && p <= cfgCutPiP1) ||
           (std::fabs(track.tpcNSigmaPi() - cfgMcTpcShiftPi) < cfgCutNSig2 && p > cfgCutPiP1 && p <= cfgCutPiP2))) ||
         (track.hasTOF() &&
          std::fabs(track.tpcNSigmaPi() - cfgMcTpcShiftPi) < cfgCutNSig3 &&
          std::fabs(track.tofNSigmaPi() - cfgMcTofShiftPi) < cfgCutNSig3))) {

      if (std::abs(track.rapidity(MassPiPlus)) < cfgCutRap) {
        return true;
      }
    }
    return false;
  }

  // PID selction cuts for Low momentum Kaons
  template <typename T>
  bool selLowKa(T const& track, double p)
  {
    if (track.pt() >= cfgCutKaPtMin &&
        track.p() <= cfgCutKaThrsldP &&
        ((!track.hasTOF() &&
          ((std::fabs(track.tpcNSigmaKa() - cfgMcTpcShiftKa) < cfgCutNSig3 && p <= cfgCutKaP1) ||
           (std::fabs(track.tpcNSigmaKa() - cfgMcTpcShiftKa) < cfgCutNSig2 && p > cfgCutKaP1 && p <= cfgCutKaP2))) ||
         (track.hasTOF() &&
          std::fabs(track.tpcNSigmaKa() - cfgMcTpcShiftKa) < cfgCutNSig3 &&
          std::fabs(track.tofNSigmaKa() - cfgMcTofShiftKa) < cfgCutNSig3))) {

      if (std::abs(track.rapidity(MassKPlus)) < cfgCutRap) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for Low momentum Protons
  template <typename T>
  bool selLowPr(T const& track, double p)
  {
    if (track.pt() >= cfgCutPrPtMin &&
        track.p() <= cfgCutPrThrsldP &&
        ((!track.hasTOF() &&
          ((std::fabs(track.tpcNSigmaPr() - cfgMcTpcShiftPr) < cfgCutNSig3 && p <= cfgCutPrP1) ||
           (std::fabs(track.tpcNSigmaPr() - cfgMcTpcShiftPr) < cfgCutNSig2 && p > cfgCutPrP1 && p <= cfgCutPrP2))) ||
         (track.hasTOF() &&
          std::fabs(track.tpcNSigmaPr() - cfgMcTpcShiftPr) < cfgCutNSig3 &&
          std::fabs(track.tofNSigmaPr() - cfgMcTofShiftPr) < cfgCutNSig3))) {

      if (std::abs(track.rapidity(MassProton)) < cfgCutRap) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for High momentum Protons
  template <typename T>
  bool selHighPi(T const& track)
  {
    if (track.hasTOF() &&
        track.p() > cfgCutPiThrsldP &&
        std::fabs(track.tpcNSigmaPi() - cfgMcTpcShiftPi) < cfgCutNSig3 &&
        std::fabs(track.tofNSigmaPi() - cfgMcTofShiftPi) < cfgCutNSig3) {

      if (std::abs(track.rapidity(MassPiPlus)) < cfgCutRap) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for High momentum Kaons
  template <typename T>
  bool selHighKa(T const& track)
  {
    if (track.hasTOF() &&
        track.p() > cfgCutKaThrsldP &&
        std::fabs(track.tpcNSigmaKa() - cfgMcTpcShiftKa) < cfgCutNSig3 &&
        ((std::fabs(track.tofNSigmaKa() - cfgMcTofShiftKa) < cfgCutNSig3 && track.p() <= cfgCutKaP3) ||
         (std::fabs(track.tofNSigmaKa() - cfgMcTofShiftKa) < cfgCutNSig2 && track.p() > cfgCutKaP3))) {

      if (std::abs(track.rapidity(MassKPlus)) < cfgCutRap) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for High momentum Protons
  template <typename T>
  bool selHighPr(T const& track)
  {
    if (track.hasTOF() &&
        track.p() > cfgCutPrThrsldP &&
        std::fabs(track.tpcNSigmaPr() - cfgMcTpcShiftPr) < cfgCutNSig3 &&
        std::fabs(track.tofNSigmaPr() - cfgMcTofShiftPr) < cfgCutNSig3) {

      if (std::abs(track.rapidity(MassProton)) < cfgCutRap) {
        return true;
      }
    }

    return false;
  }

  // Fill hist before selection cuts:
  template <typename T, typename U>
  void FillBeforeQAHistos(T const& col, U const& tracks)
  {
    for (auto& track : tracks) {
      hist.fill(HIST("QA/before/h_Eta"), track.eta());
      hist.fill(HIST("QA/before/h_Phi"), track.phi());
      hist.fill(HIST("QA/before/h_Pt"), track.pt());
      hist.fill(HIST("QA/before/h2_PvsPinner"), track.p(), track.tpcInnerParam());
      hist.fill(HIST("QA/before/h2_Pt_Eta"), track.eta(), track.pt());
      hist.fill(HIST("QA/before/h_TPCChi2perCluster"), track.tpcChi2NCl());
      hist.fill(HIST("QA/before/h_ITSChi2perCluster"), track.itsChi2NCl());
      hist.fill(HIST("QA/before/h_crossedTPC"), track.tpcNClsCrossedRows());
      hist.fill(HIST("QA/before/h_DcaXY"), track.dcaXY());
      hist.fill(HIST("QA/before/h_DcaZ"), track.dcaZ());
      hist.fill(HIST("QA/before/h2_DcaXY"), track.pt(), track.dcaXY());
      hist.fill(HIST("QA/before/h2_DcaZ"), track.pt(), track.dcaZ());

      if (track.hasTOF())
        hist.fill(HIST("QA/before/h2_PtofvsPinner"), track.p(), track.tpcInnerParam());
    }
    hist.fill(HIST("QA/before/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/before/h_Counts"), 2);

    int NTPC = col.multNTracksHasTPC();
    int N_FT0M = col.multFT0M();
    int N_FT0C = col.multFT0C();
    double cent_FT0C = col.centFT0C();

    if (NTPC != 0 && N_FT0M != 0) {
      hist.fill(HIST("QA/before/h_NTPC"), NTPC);
      hist.fill(HIST("QA/before/h_Cent"), cent_FT0C);
      hist.fill(HIST("QA/before/h_NFT0M"), N_FT0M);
      hist.fill(HIST("QA/before/h_NFT0C"), N_FT0M);
      hist.fill(HIST("QA/before/h2_NTPC_NFT0M"), N_FT0M, NTPC);
      hist.fill(HIST("QA/before/h2_NTPC_NFT0C"), N_FT0C, NTPC);
      hist.fill(HIST("QA/before/h2_NTPC_Cent"), cent_FT0C, NTPC);
    }
  }

  // Fill hist after selection cuts:
  template <typename T>
  void FillAfterQAHistos(T const& col)
  {
    int NTPC = col.multNTracksHasTPC();
    int N_FT0M = col.multFT0M();
    int N_FT0C = col.multFT0C();
    double cent_FT0C = col.centFT0C();

    hist.fill(HIST("QA/after/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/after/h_Counts"), 2);
    if (NTPC != 0 && N_FT0M != 0) {
      hist.fill(HIST("QA/after/h_NTPC"), NTPC);
      hist.fill(HIST("QA/after/h_Cent"), cent_FT0C);
      hist.fill(HIST("QA/after/h_NFT0M"), N_FT0M);
      hist.fill(HIST("QA/after/h_NFT0C"), N_FT0C);
      hist.fill(HIST("QA/after/h2_NTPC_NFT0M"), N_FT0M, NTPC);
      hist.fill(HIST("QA/after/h2_NTPC_NFT0C"), N_FT0C, NTPC);
      hist.fill(HIST("QA/after/h2_NTPC_Cent"), cent_FT0C, NTPC);
      hist.fill(HIST("QA/after/p_NTPC_Cent"), cent_FT0C, NTPC);
      hist.fill(HIST("QA/after/p_NTPC_NFT0M"), N_FT0M, NTPC);
      hist.fill(HIST("QA/after/p_NTPC_NFT0C"), N_FT0C, NTPC);
    }
  }

  // Fill Charged particles QA:
  template <typename T>
  void FillChargedQAHistos(T const& track)
  {
    hist.fill(HIST("QA/after/h_Eta"), track.eta());
    hist.fill(HIST("QA/after/h_Phi"), track.phi());
    hist.fill(HIST("QA/after/h_Pt"), track.pt());
    hist.fill(HIST("QA/after/h2_PvsPinner"), track.p(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/h2_Pt_Eta"), track.eta(), track.pt());
    hist.fill(HIST("QA/after/h_DcaZ"), track.dcaZ());
    hist.fill(HIST("QA/after/h_DcaXY"), track.dcaXY());
    hist.fill(HIST("QA/after/h2_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST("QA/after/h2_DcaZ"), track.pt(), track.dcaZ());

    hist.fill(HIST("QA/after/h_TPCChi2perCluster"), track.tpcChi2NCl());
    hist.fill(HIST("QA/after/h_ITSChi2perCluster"), track.itsChi2NCl());
    hist.fill(HIST("QA/after/h_crossedTPC"), track.tpcNClsCrossedRows());

    if (track.hasTOF())
      hist.fill(HIST("QA/after/h2_PtofvsPinner"), track.p(), track.tpcInnerParam());
  }

  // Fill before PID cut QA hist:
  template <typename T>
  void FillBeforePIDQAHistos(T const& track)
  {
    hist.fill(HIST("QA/before/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/before/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/before/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST("QA/Pion/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPi());
    hist.fill(HIST("QA/Pion/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/before/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Proton/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPr());
    hist.fill(HIST("QA/Proton/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/before/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Kaon/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaKa());
    hist.fill(HIST("QA/Kaon/before/h2_TOFNsigma"), track.p(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/before/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());

    hist.fill(HIST("QA/before/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/before/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/before/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());

    hist.fill(HIST("QA/Pion/innerParam/before/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaPi());
    hist.fill(HIST("QA/Pion/innerParam/before/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/innerParam/before/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Proton/innerParam/before/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaPr());
    hist.fill(HIST("QA/Proton/innerParam/before/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/innerParam/before/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Kaon/innerParam/before/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaKa());
    hist.fill(HIST("QA/Kaon/innerParam/before/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/innerParam/before/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());
  }

  // Moments Calculation:
  void moments(double pt, double& Q1, double& Q2, double& Q3, double& Q4)
  {
    Q1 += pt;
    Q2 += pt * pt;
    Q3 += pt * pt * pt;
    Q4 += pt * pt * pt * pt;
  }

  // Fill after PID cut QA hist:
  template <int mode, typename T>
  void FillIdParticleQAHistos(T const& track, double rap, double nSigmaTPC, double nSigmaTOF, int& N, double& Q1, double& Q2, double& Q3, double& Q4)
  {
    N++;
    moments(track.pt(), Q1, Q2, Q3, Q4);

    hist.fill(HIST(dire[mode]) + HIST("h_Pt"), track.pt());
    if (track.sign() > 0)
      hist.fill(HIST(dire[mode]) + HIST("h_PtPos"), track.pt());

    if (track.sign() < 0)
      hist.fill(HIST(dire[mode]) + HIST("h_PtNeg"), track.pt());

    hist.fill(HIST(dire[mode]) + HIST("h_Eta"), track.eta());
    hist.fill(HIST(dire[mode]) + HIST("h_Phi"), track.phi());
    hist.fill(HIST(dire[mode]) + HIST("h_rap"), rap);
    hist.fill(HIST(dire[mode]) + HIST("h2_Pt_rap"), rap, track.pt());
    hist.fill(HIST(dire[mode]) + HIST("h_DcaZ"), track.dcaZ());
    hist.fill(HIST(dire[mode]) + HIST("h_DcaXY"), track.dcaXY());
    hist.fill(HIST(dire[mode]) + HIST("h2_DcaZ"), track.pt(), track.dcaZ());
    hist.fill(HIST(dire[mode]) + HIST("h2_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST(dire[mode]) + HIST("h2_Pt_Pinner"), track.tpcInnerParam(), track.pt());
    hist.fill(HIST(dire[mode]) + HIST("h2_P_Pinner"), track.tpcInnerParam(), track.p());

    hist.fill(HIST(dire[mode]) + HIST("h2_TPCNsigma_El"), track.p(), track.tpcNSigmaEl());
    hist.fill(HIST(dire[mode]) + HIST("h2_TOFNsigma_El"), track.p(), track.tofNSigmaEl());
    hist.fill(HIST(dire[mode]) + HIST("h2_TPCNsigma"), track.p(), nSigmaTPC);
    hist.fill(HIST(dire[mode]) + HIST("h2_TOFNsigma"), track.p(), nSigmaTOF);
    hist.fill(HIST(dire[mode]) + HIST("h2_TpcTofNsigma"), nSigmaTPC, nSigmaTOF);
    hist.fill(HIST(dire[mode]) + HIST("h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST(dire[mode]) + HIST("h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST(dire[mode]) + HIST("h2_pvsm2"), track.mass() * track.mass(), track.p());
    hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/after/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST(dire[mode]) + HIST("innerParam/h2_TPCNsigma_El"), track.tpcInnerParam(), track.tpcNSigmaEl());
    hist.fill(HIST(dire[mode]) + HIST("innerParam/h2_TOFNsigma_El"), track.tpcInnerParam(), track.tofNSigmaEl());
    hist.fill(HIST(dire[mode]) + HIST("innerParam/h2_TPCNsigma"), track.tpcInnerParam(), nSigmaTPC);
    hist.fill(HIST(dire[mode]) + HIST("innerParam/h2_TOFNsigma"), track.tpcInnerParam(), nSigmaTOF);
    hist.fill(HIST(dire[mode]) + HIST("innerParam/h2_TpcTofNsigma"), nSigmaTPC, nSigmaTOF);
    hist.fill(HIST(dire[mode]) + HIST("innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST(dire[mode]) + HIST("innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST(dire[mode]) + HIST("innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/after/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/after/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
  }

  template <int mode>
  void FillPtMCHist(double pt, int PID, int pdgCodePos, int pdgCodeNeg)
  {
    hist.fill(HIST(dire[mode]) + HIST("h_PtTruth"), pt);
    if (PID == pdgCodePos) {
      hist.fill(HIST(dire[mode]) + HIST("h_PtPosTruth"), pt);
    }
    if (PID == pdgCodeNeg) {
      hist.fill(HIST(dire[mode]) + HIST("h_PtNegTruth"), pt);
    }
  }

  template <int mode>
  void FillAnalysisHistos(int NTPC, int N_FT0M, int N, double Q1, double Q2, double Q3, double Q4)
  {
    if (N == 0 /* || Q1 == 0 || Q2 == 0 || Q3 == 0 || Q4 == 0 */) {
      return;
    }
    double twopart1 = ((Q1 * Q1) - Q2);
    double threepart1 = ((Q1 * Q1 * Q1) - (3 * Q2 * Q1) + 2 * Q3);
    double fourpart1 = ((Q1 * Q1 * Q1 * Q1) - (6 * Q2 * Q1 * Q1) + (3 * Q2 * Q2) + (8 * Q3 * Q1) - 6 * Q4);

    hist.fill(HIST(dire[mode]) + HIST("h_Mult"), N);
    hist.fill(HIST(dire[mode]) + HIST("h_Q1"), NTPC, Q1, N_FT0M);
    hist.fill(HIST(dire[mode]) + HIST("h_Q2"), NTPC, Q2, N_FT0M);
    hist.fill(HIST(dire[mode]) + HIST("h_Q3"), NTPC, Q3, N_FT0M);
    hist.fill(HIST(dire[mode]) + HIST("h_Q4"), NTPC, Q4, N_FT0M);

    if (N > 1) {
      double mean_pT = Q1 / static_cast<double>(N);
      double N_pair = (static_cast<double>(N) * (static_cast<double>(N) - 1));
      double twopart = twopart1 / N_pair;
      double checkN_deno_var = (1 / std::sqrt(1 - (1 / static_cast<double>(N))));
      hist.fill(HIST(dire[mode]) + HIST("h_mean_pT"), mean_pT);
      hist.fill(HIST(dire[mode]) + HIST("p_mean_pT_Mult_var"), NTPC, mean_pT);

      hist.fill(HIST(dire[mode]) + HIST("h_Q1_var"), NTPC, Q1, N_FT0M);
      hist.fill(HIST(dire[mode]) + HIST("h_N_var"), NTPC, N, N_FT0M);
      hist.fill(HIST(dire[mode]) + HIST("h_twopart_nume_Mult_var"), NTPC, twopart1, N_FT0M);
      hist.fill(HIST(dire[mode]) + HIST("h_twopart_deno_Mult_var"), NTPC, N_pair, N_FT0M);
      hist.fill(HIST(dire[mode]) + HIST("h_mean_pT_Mult_var"), NTPC, mean_pT, N_FT0M);
      hist.fill(HIST(dire[mode]) + HIST("h_twopart_Mult_var"), NTPC, twopart, N_FT0M);
      hist.fill(HIST(dire[mode]) + HIST("p_CheckNCH"), NTPC, checkN_deno_var);
      hist.fill(HIST(dire[mode]) + HIST("h_CheckNCH"), NTPC, checkN_deno_var, N_FT0M);

      if (N > 2) {
        double N_triplet = (static_cast<double>(N) * (static_cast<double>(N) - 1) * (static_cast<double>(N) - 2));
        double threepart = threepart1 / N_triplet;
        hist.fill(HIST(dire[mode]) + HIST("h_mean_pT_Mult_skew"), NTPC, mean_pT, N_FT0M);
        hist.fill(HIST(dire[mode]) + HIST("h_twopart_Mult_skew"), NTPC, twopart, N_FT0M);
        hist.fill(HIST(dire[mode]) + HIST("h_threepart_Mult_skew"), NTPC, threepart, N_FT0M);

        if (N > 3) {
          double N_quad = (static_cast<double>(N) * (static_cast<double>(N) - 1) * (static_cast<double>(N) - 2) * (static_cast<double>(N) - 3));
          double fourpart = fourpart1 / N_quad;
          hist.fill(HIST(dire[mode]) + HIST("h_mean_pT_Mult_kurto"), NTPC, mean_pT, N_FT0M);
          hist.fill(HIST(dire[mode]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart, N_FT0M);
          hist.fill(HIST(dire[mode]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart, N_FT0M);
          hist.fill(HIST(dire[mode]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart, N_FT0M);
        }
      }
    }
  }

  template <bool DataFlag, bool RecoFlag, typename T, typename U>
  void FillHistos(T const& col, U const& tracks)
  {
    int Nch = 0, N_TPC = 0, N_FT0M = 0;
    int N_Pi = 0, N_Ka = 0, N_Pr = 0;
    double pt_ch = 0, Q1_ch = 0, Q2_ch = 0, Q3_ch = 0, Q4_ch = 0;
    double pt_Pi = 0, Q1_Pi = 0, Q2_Pi = 0, Q3_Pi = 0, Q4_Pi = 0;
    double pt_Pr = 0, Q1_Pr = 0, Q2_Pr = 0, Q3_Pr = 0, Q4_Pr = 0;
    double pt_Ka = 0, Q1_Ka = 0, Q2_Ka = 0, Q3_Ka = 0, Q4_Ka = 0;

    int Nch_sim = 0, N_sim = 0, NFT0C_sim = 0;
    int N_Pi_sim = 0, N_Ka_sim = 0, N_Pr_sim = 0;
    double pt_ch_sim = 0, Q1_ch_sim = 0, Q2_ch_sim = 0, Q3_ch_sim = 0, Q4_ch_sim = 0;
    double pt_Pi_sim = 0, Q1_Pi_sim = 0, Q2_Pi_sim = 0, Q3_Pi_sim = 0, Q4_Pi_sim = 0;
    double pt_Pr_sim = 0, Q1_Pr_sim = 0, Q2_Pr_sim = 0, Q3_Pr_sim = 0, Q4_Pr_sim = 0;
    double pt_Ka_sim = 0, Q1_Ka_sim = 0, Q2_Ka_sim = 0, Q3_Ka_sim = 0, Q4_Ka_sim = 0;

    array<float, 3> p1, p2;
    double invMassGamma = 0.0;

    for (auto const& [trkEl, trkPos] : soa::combinations(soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (trkEl.index() == trkPos.index())
        continue;

      if (!selTrack(trkEl) || !selTrack(trkPos))
        continue;

      if (!selElectrons(trkEl) || !selElectrons(trkPos))
        continue;

      p1 = array{trkEl.px(), trkEl.py(), trkEl.pz()};
      p2 = array{trkPos.px(), trkPos.py(), trkPos.pz()};

      invMassGamma = RecoDecay::m(array{p1, p2}, array{MassElectron, MassElectron});
      hist.fill(HIST("QA/after/h_invMass_gamma"), invMassGamma);
    }

    if constexpr (DataFlag) {
      for (auto& track : tracks) {
        if (!selTrack(track)) {
          continue;
        }

        double nSigmaTPCPi = track.tpcNSigmaPi();
        double nSigmaTPCKa = track.tpcNSigmaKa();
        double nSigmaTPCPr = track.tpcNSigmaPr();
        double nSigmaTOFPi = track.tofNSigmaPi();
        double nSigmaTOFKa = track.tofNSigmaKa();
        double nSigmaTOFPr = track.tofNSigmaPr();
        double rapPi = track.rapidity(MassPiPlus);
        double rapKa = track.rapidity(MassKPlus);
        double rapPr = track.rapidity(MassProton);
        double innerParam = track.tpcInnerParam();

        if (std::fabs(track.eta()) < 0.8) {
          Nch++;
          pt_ch = track.pt();
          moments(pt_ch, Q1_ch, Q2_ch, Q3_ch, Q4_ch);
          FillChargedQAHistos(track);
        }

        FillBeforePIDQAHistos(track);

        if (rejectTracks(track)) {
          return;
        }

        if (cfgInvMass == true && invMassGamma < cfgGammaCut) {
          continue;
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPi(track, innerParam) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
            FillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, N_Pi, Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPi(track, innerParam) == cfgSelLow && selHighPi(track) == cfgSelHigh) {
            FillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, N_Pi, Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi);
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowKa(track, innerParam) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
            FillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, N_Ka, Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowKa(track, innerParam) == cfgSelLow && selHighKa(track) == cfgSelHigh) {
            FillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, N_Ka, Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka);
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPr(track, innerParam) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
            FillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, N_Pr, Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPr(track, innerParam) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
            FillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, N_Pr, Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr);
          }
        }
      }
    } else if constexpr (RecoFlag) {
      if (!col.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }

      for (auto& track : tracks) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC Particle for this track, skip...");
          continue;
        }
        auto mcPart = track.mcParticle();
        int PID = mcPart.pdgCode();
        if (!mcPart.isPhysicalPrimary()) {
          continue;
        }

        //___________________________________Truth Level____________________________________________________//
        auto charge = 0.;
        auto* pd = pdg->GetParticle(PID);
        if (pd != nullptr) {
          charge = pd->Charge();
        }
        if (std::fabs(charge) < 1e-3) {
          continue;
        }
        if (std::abs(mcPart.eta()) < 0.8) {
          N_sim++;
        }
        if (-3.3 < mcPart.eta() && mcPart.eta() < -2.1) {
          NFT0C_sim++;
        }

        if (mcPart.pt() > cfgCutPtMin && mcPart.pt() < cfgCutPtMax) {
          if (std::abs(mcPart.y()) > cfgCutRap) {
            continue;
          }

          Nch_sim++;
          pt_ch_sim = mcPart.pt();
          moments(pt_ch_sim, Q1_ch_sim, Q2_ch_sim, Q3_ch_sim, Q4_ch_sim);
          hist.fill(HIST("Gen/Charged/h_PtTruth"), mcPart.pt());

          if (std::abs(PID) == kPiPlus && mcPart.pt() >= cfgCutPiPtMin) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if (mcPart.p() <= cfgCutPiThrsldP || mcPart.p() > cfgCutPiThrsldP) {
                N_Pi_sim++;
                pt_Pi_sim = mcPart.pt();
                moments(pt_Pi_sim, Q1_Pi_sim, Q2_Pi_sim, Q3_Pi_sim, Q4_Pi_sim);
                FillPtMCHist<Gen_Pion>(pt_Pi_sim, PID, kPiPlus, kPiMinus);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPiThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutPiThrsldP)) {
                N_Pi_sim++;
                pt_Pi_sim = mcPart.pt();
                moments(pt_Pi_sim, Q1_Pi_sim, Q2_Pi_sim, Q3_Pi_sim, Q4_Pi_sim);
                FillPtMCHist<Gen_Pion>(pt_Pi_sim, PID, kPiPlus, kPiMinus);
              }
            }
          }

          if (std::abs(PID) == kKPlus && mcPart.pt() >= cfgCutKaPtMin) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPiThrsldP) || (cfgSelHigh == true && mcPart.p() > cfgCutPiThrsldP)) {
                N_Ka_sim++;
                pt_Ka_sim = mcPart.pt();
                moments(pt_Ka_sim, Q1_Ka_sim, Q2_Ka_sim, Q3_Ka_sim, Q4_Ka_sim);
                FillPtMCHist<Gen_Kaon>(pt_Ka_sim, PID, kKPlus, kKMinus);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutKaThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutKaThrsldP)) {
                N_Ka_sim++;
                pt_Ka_sim = mcPart.pt();
                moments(pt_Ka_sim, Q1_Ka_sim, Q2_Ka_sim, Q3_Ka_sim, Q4_Ka_sim);
                FillPtMCHist<Gen_Kaon>(pt_Ka_sim, PID, kKPlus, kKMinus);
              }
            }
          }

          if (std::abs(PID) == kProton && mcPart.pt() >= cfgCutPrPtMin) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPrThrsldP) || (cfgSelHigh == true && mcPart.p() > cfgCutPrThrsldP)) {
                N_Pr_sim++;
                pt_Pr_sim = mcPart.pt();
                moments(pt_Pr_sim, Q1_Pr_sim, Q2_Pr_sim, Q3_Pr_sim, Q4_Pr_sim);
                FillPtMCHist<Gen_Proton>(pt_Pr_sim, PID, kProton, kProtonBar);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPrThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutPrThrsldP)) {
                N_Pr_sim++;
                pt_Pr_sim = mcPart.pt();
                moments(pt_Pr_sim, Q1_Pr_sim, Q2_Pr_sim, Q3_Pr_sim, Q4_Pr_sim);
                FillPtMCHist<Gen_Proton>(pt_Pr_sim, PID, kProton, kProtonBar);
              }
            }
          }
        }

        //______________________________Reconstructed Level____________________________________________________//

        if (!selTrack(track)) {
          continue;
        }

        double nSigmaTPCPi = track.tpcNSigmaPi();
        double nSigmaTPCKa = track.tpcNSigmaKa();
        double nSigmaTPCPr = track.tpcNSigmaPr();
        double nSigmaTOFPi = track.tofNSigmaPi();
        double nSigmaTOFKa = track.tofNSigmaKa();
        double nSigmaTOFPr = track.tofNSigmaPr();
        double rapPi = track.rapidity(MassPiPlus);
        double rapKa = track.rapidity(MassKPlus);
        double rapPr = track.rapidity(MassProton);
        double innerParam = track.tpcInnerParam();

        if (std::fabs(track.eta()) < 0.8) {
          Nch++;
          pt_ch = track.pt();
          moments(pt_ch, Q1_ch, Q2_ch, Q3_ch, Q4_ch);
          FillChargedQAHistos(track);
        }
        FillBeforePIDQAHistos(track);

        if (rejectTracks(track)) {
          return;
        }

        if (cfgInvMass == true && invMassGamma < cfgGammaCut) {
          continue;
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPi(track, innerParam) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
            pt_Pi = track.pt();
            FillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, N_Pi, Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi);

            if (std::abs(PID) == kPiPlus) {
              FillPtMCHist<QA_Pion>(pt_Pi, PID, kPiPlus, kPiMinus);
            }
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPi(track, innerParam) == cfgSelLow && selHighPi(track) == cfgSelHigh) {
            pt_Pi = track.pt();
            FillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, N_Pi, Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi);

            if (std::abs(PID) == kPiPlus) {
              FillPtMCHist<QA_Pion>(pt_Pi, PID, kPiPlus, kPiMinus);
            }
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowKa(track, innerParam) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
            pt_Ka = track.pt();
            FillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, N_Ka, Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka);
            if (std::abs(PID) == kKPlus) {
              FillPtMCHist<QA_Kaon>(pt_Ka, PID, kKPlus, kKMinus);
            }
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowKa(track, innerParam) == cfgSelLow && selHighKa(track) == cfgSelHigh) {
            pt_Ka = track.pt();
            FillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, N_Ka, Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka);

            if (std::abs(PID) == kKPlus) {
              FillPtMCHist<QA_Kaon>(pt_Ka, PID, kKPlus, kKMinus);
            }
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPr(track, innerParam) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
            pt_Pr = track.pt();
            FillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, N_Pr, Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr);

            if (std::abs(PID) == kProton) {
              FillPtMCHist<QA_Proton>(pt_Pr, PID, kProton, kProtonBar);
            }
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPr(track, innerParam) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
            pt_Pr = track.pt();
            FillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, N_Pr, Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr);

            if (std::abs(PID) == kProton) {
              FillPtMCHist<QA_Proton>(pt_Pr, PID, kProton, kProtonBar);
            }
          }
        }
      }
      hist.fill(HIST("QA/after/h_vtxZ_sim"), col.mcCollision().posZ());
    }

    N_TPC = col.multNTracksHasTPC();
    N_FT0M = col.multFT0M();

    if (cfgMCTruth) {
      if (N_sim != 0)
        hist.fill(HIST("QA/after/h_NSim"), N_sim);

      if (N_sim != 0 && Nch != 0)
        hist.fill(HIST("QA/after/h2_NchSim_NSim"), N_sim, Nch_sim);

      if (N_sim != 0 && N_TPC != 0)
        hist.fill(HIST("QA/after/h2_NTPC_NSim"), N_sim, N_TPC);

      int NFT0C = col.multFT0C();
      if (NFT0C != 0 && NFT0C_sim != 0)
        hist.fill(HIST("QA/after/h2_NFT0C_NFT0CSim"), NFT0C_sim, NFT0C);

      N_TPC = N_sim;

      hist.fill(HIST("Gen/NTPC"), N_TPC);
      hist.fill(HIST("Gen/NFT0C"), NFT0C_sim);
      hist.fill(HIST("Gen/h2_NTPC_NFT0C"), NFT0C_sim, N_TPC);

      FillAnalysisHistos<Gen_Charged>(N_TPC, N_FT0M, Nch_sim, Q1_ch_sim, Q2_ch_sim, Q3_ch_sim, Q4_ch_sim);
      FillAnalysisHistos<Gen_Pion>(N_TPC, N_FT0M, N_Pi_sim, Q1_Pi_sim, Q2_Pi_sim, Q3_Pi_sim, Q4_Pi_sim);
      FillAnalysisHistos<Gen_Kaon>(N_TPC, N_FT0M, N_Ka_sim, Q1_Ka_sim, Q2_Ka_sim, Q3_Ka_sim, Q4_Ka_sim);
      FillAnalysisHistos<Gen_Proton>(N_TPC, N_FT0M, N_Pr_sim, Q1_Pr_sim, Q2_Pr_sim, Q3_Pr_sim, Q4_Pr_sim);
    }

    FillAfterQAHistos(col);
    if (N_TPC != 0 && Nch != 0)
      hist.fill(HIST("QA/after/h2_NTPC_Nch"), N_TPC, Nch);

    FillAnalysisHistos<Analysis_Charged>(N_TPC, N_FT0M, Nch, Q1_ch, Q2_ch, Q3_ch, Q4_ch);
    FillAnalysisHistos<Analysis_Pion>(N_TPC, N_FT0M, N_Pi, Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi);
    FillAnalysisHistos<Analysis_Kaon>(N_TPC, N_FT0M, N_Ka, Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka);
    FillAnalysisHistos<Analysis_Proton>(N_TPC, N_FT0M, N_Pr, Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr);
  }

  void process_Run3(MyCollisions::iterator const& col, MyAllTracks const& tracks)
  {
    // Before Collision and Track Cuts:
    FillBeforeQAHistos(col, tracks);

    // After Collision and Track Cuts:
    if (selRun3Col(col)) {
      FillHistos<true, false>(col, tracks);
    }
  }
  PROCESS_SWITCH(meanPtFlucId, process_Run3, "Process for Run3", false);

  void process_MCRecoSimRun3(MyMCCollisions::iterator const& col, aod::McCollisions const&, MyMCTracks const& tracks, aod::McParticles const&)
  {
    // Before Collision and Track Cuts:
    FillBeforeQAHistos(col, tracks);

    // After Collision and Track Cuts:
    if (selRun3Col(col)) {
      FillHistos<false, true>(col, tracks);
    }
  }
  PROCESS_SWITCH(meanPtFlucId, process_MCRecoSimRun3, "process MC Reconstructed & Truth Run-3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<meanPtFlucId>(cfgc)};
}
