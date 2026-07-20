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

/// \file meanPtFlucId.cxx
/// \brief Calculate EbyE <pt> fluctuations with cumulant method.
///        For charged particles and identified particles.
///        For RUN-3
///
/// \author Tanu Gahlaut <tanu.gahlaut@cern.ch>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace std;

struct MeanPtFlucId {
  Configurable<int> nPBins{"nPBins", 300, ""};
  Configurable<int> nPartBins{"nPartBins", 100, ""};
  Configurable<int> nPhiBins{"nPhiBins", 100, ""};
  Configurable<float> cfgCutPtMax{"cfgCutPtMax", 2.0, "maximum pT"};
  Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.2, "minimum pT"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut"};
  Configurable<float> cfgCutDcaZ{"cfgCutDcaZ", 0.15, "DCAz cut"};
  Configurable<float> cfgCutPosZ{"cfgCutPosZ", 7.0, "cut for vertex Z"};
  Configurable<bool> cfgPosZ{"cfgPosZ", true, "Position Z"};
  Configurable<float> cfgCutNSig2{"cfgCutNSig2", 2.0, "nSigma cut: 2"};
  Configurable<float> cfgCutNSig3{"cfgCutNSig3", 3.0, "nSigma cut: 3"};
  Configurable<float> cfgCutNSig5{"cfgCutNSig5", 5.0, "nSigma cut: 5"};
  Configurable<float> cfgSelCutNSigPi{"cfgSelCutNSigPi", 2.0, "nSigma cut for pion selection"};
  Configurable<float> cfgSelCutNSigKa{"cfgSelCutNSigKa", 3.0, "nSigma cut for kaon selection"};
  Configurable<float> cfgSelCutNSigPr{"cfgSelCutNSigPr", 2.0, "nSigma cut for proton selection"};
  Configurable<float> cfgCutPiPtMin{"cfgCutPiPtMin", 0.2, "Minimum pion p_{T} cut"};
  Configurable<float> cfgCutKaPtMin{"cfgCutKaPtMin", 0.3, "Minimum kaon p_{T} cut"};
  Configurable<float> cfgCutPrPtMin{"cfgCutPrPtMin", 0.5, "Minimum proton p_{T} cut"};
  Configurable<float> cfgCutPiThrsldP{"cfgCutPiThrsldP", 0.6, "Threshold p cut pion"};
  Configurable<float> cfgCutKaThrsldP{"cfgCutKaThrsldP", 0.6, "Threshold p cut kaon"};
  Configurable<float> cfgCutPrThrsldP{"cfgCutPrThrsldP", 1.0, "Threshold p cut proton "};
  Configurable<float> cfgP0U{"cfgP0U", 6.36269, "p_{0} for upper cut for N_{TPC} vs N{sim}"};
  Configurable<float> cfgP1U{"cfgP1U", 2.22108, "p_{1} for upper cut for N_{TPC} vs N{sim}"};
  Configurable<float> cfgP2U{"cfgP2U", -0.0118451, "p_{2} for upper cut for N_{TPC} vs N{sim}"};
  Configurable<float> cfgP0L{"cfgP0L", 23.8213, "p_{0} for lower cut for N_{TPC} vs N{sim}"};
  Configurable<float> cfgP1L{"cfgP1L", -0.402328, "p_{1} for lower cut for N_{TPC} vs N{sim}"};
  Configurable<float> cfgP2L{"cfgP2L", 0.0128423, "p_{2} for lower cut for N_{TPC} vs N{sim}"};
  Configurable<float> cfgP0L1{"cfgP0L1", -10.5149, "p_{0} for lower cut for N_{TPC} vs N{sim}"};
  Configurable<float> cfgP1L1{"cfgP1L1", 0.936596, "p_{1} for lower cut for N_{TPC} vs N{sim}"};
  Configurable<float> cfgP0corr{"cfgP0corr", 0.0001, "p_{0} for corrected N_{TPC} "};
  Configurable<float> cfgP1corr{"cfgP1corr", 0.0001, "p_{1} for corrected N_{TPC} "};
  Configurable<float> cfgP2corr{"cfgP2corr", 0.0001, "p_{2} for corrected N_{TPC} "};
  Configurable<bool> cfgSel8{"cfgSel8", true, "Sel8 trigger"};
  Configurable<bool> cfgNoSameBunchPileup{"cfgNoSameBunchPileup", true, "kNoSameBunchPileup"};
  Configurable<bool> cfgIsVertexITSTPC{"cfgIsVertexITSTPC", true, "kIsVertexITSTPC"};
  Configurable<bool> cfgLightIonCuts{"cfgLightIonCuts", false, "Light Ion Cuts"};
  Configurable<bool> cfgIsGoodITSLayersAll{"cfgIsGoodITSLayersAll", true, "kIsGoodITSLayersAll"};
  Configurable<bool> cfgIsGoodZvtxFT0vsPV{"cfgIsGoodZvtxFT0vsPV", true, "kIsGoodZvtxFT0vsPV"};
  Configurable<bool> cfgRejTrk{"cfgRejTrk", true, "Rejected Tracks"};
  ConfigurableAxis multTPCBins{"multTPCBins", {150, 0, 150}, "TPC Multiplicity bins"};
  ConfigurableAxis multMCBins{"multMCBins", {300, 0, 300}, "MC Multiplicity bins"};
  ConfigurableAxis multCorrBins{"multCorrBins", {100, 0., 150.}, "Corrected TPC Multiplicity bins"};
  ConfigurableAxis multFT0MBins{"multFT0MBins", {1000, 0, 5000}, "Forward Multiplicity bins"};
  ConfigurableAxis dcaXYBins{"dcaXYBins", {100, -0.15, 0.15}, "dcaXY bins"};
  ConfigurableAxis dcaZBins{"dcaZBins", {500, -1.2, 1.2}, "dcaZ bins"};
  ConfigurableAxis centBins{"centBins", {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100}, "Centrality bins"};

  Configurable<std::vector<double>> ptBins{"ptBins", {0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00}, "p_{T} bins"};
  Configurable<std::vector<double>> etaBins{"etaBins", {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8}, "#eta bins"};

  Configurable<std::vector<double>> ptBinsPi{"ptBinsPi", {0.20, 0.225, 0.25, 0.275, 0.30, 0.325, 0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50, 0.525, 0.55, 0.575, 0.60, 0.70, 0.90, 1.1, 1.3, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00}, "p_{T} bins for pions"};
  Configurable<std::vector<double>> ptBinsKa{"ptBinsKa", {0.30, 0.325, 0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.8, 0.85, 0.90, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.70, 1.80, 1.90, 2.00}, "p_{T} bins for kaons"};
  Configurable<std::vector<double>> ptBinsPr{"ptBinsPr", {0.40, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.60, 0.625, 0.65, 0.675, 0.70, 0.725, 0.75, 0.775, 0.80, 0.825, 0.85, 0.875, 0.90, 0.925, 0.95, 1.0, 1.025, 1.05, 1.075, 1.1, 1.15, 1.2, 1.4, 1.6, 1.8, 2.0}, "p_{T} bins for protons"};

  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};

  enum CollisionLabels {
    kTotCol = 1,
    kPassSelCol
  };

  enum TrackLabels {
    kTracksBeforeHasMcParticle = 1,
    kAllTracks,
    kAllSelPassed
  };

  enum SelCollisionLabels {
    kBeforeSelCol = 1,
    kSelColPosZ,
    kSelColSel8,
    kSelColNoSameBunchPileup,
    kSelColIsVertexITSTPC,
    kSelColIsGoodITSLayersAll,
    kSelColIsGoodZvtxFT0vsPV
  };

  enum class PIDType {
    kNone,
    kPions,
    kKaons,
    kProtons
  };

  enum Mode {
    QA_Charged = 0,
    QA_Pion,
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

  static constexpr std::string_view Dire[] = {
    "QA/after/",
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

  struct TrackInfo {
    double pt;
    double eta;
  };

  std::vector<TrackInfo> chargedTracks;
  std::vector<TrackInfo> pionTracks;
  std::vector<TrackInfo> kaonTracks;
  std::vector<TrackInfo> protonTracks;

  std::vector<TrackInfo> chargedTracksGen;
  std::vector<TrackInfo> pionTracksGen;
  std::vector<TrackInfo> kaonTracksGen;
  std::vector<TrackInfo> protonTracksGen;

  void init(InitContext const&)
  {

    const int nEta = etaBins->size() - 1;
    const int nPtCh = ptBins->size() - 1;
    const int nPtPi = ptBinsPi->size() - 1;
    const int nPtKa = ptBinsKa->size() - 1;
    const int nPtPr = ptBinsPr->size() - 1;

    const AxisSpec axisKCh{nEta * nPtCh, 0, static_cast<double>(nEta * nPtCh), "k_{ch}"};
    const AxisSpec axisKPi{nEta * nPtPi, 0, static_cast<double>(nEta * nPtPi), "k_{#pi}"};
    const AxisSpec axisKKa{nEta * nPtKa, 0, static_cast<double>(nEta * nPtKa), "k_{K}"};
    const AxisSpec axisKPr{nEta * nPtPr, 0, static_cast<double>(nEta * nPtPr), "k_{p}"};

    const AxisSpec axisCol{3, 1, 4, ""};
    const AxisSpec axisTrack{5, 1, 6, ""};
    const AxisSpec axisEvents{10, 1, 11, "Counts"};
    const AxisSpec axisEta{etaBins, "#eta"};
    const AxisSpec axisPhi{nPhiBins, 0., +7., "#phi (rad)"};
    const AxisSpec axisY{100, -0.6, 0.6, "y"};
    const AxisSpec axisPt{ptBins, "p_{T} (GeV/c)"};
    const AxisSpec axisPtPi{ptBinsPi, "p_{T} (GeV/c)"};
    const AxisSpec axisPtKa{ptBinsKa, "p_{T}(GeV/c)"};
    const AxisSpec axisPtPr{ptBinsPr, "p_{T} (GeV/c)"};
    const AxisSpec axisP{nPBins, 0., 3., "p (GeV/c)"};
    const AxisSpec axisInnerParam{nPBins, 0., 3., "p_{InnerParam} (GeV/c)"};
    const AxisSpec axisPart{nPartBins, 0., 18., " "};
    const AxisSpec axisMult{100, 0, 100, "N_{ch}"};
    const AxisSpec axisMultTPC{multTPCBins, "N_{TPC} "};
    const AxisSpec axisMultCorr{multCorrBins, "N_{Corr} "};
    const AxisSpec axisMultMC{multMCBins, "N_{Gen} "};
    const AxisSpec axisMultFT0M{multFT0MBins, "N_{FT0M}"};
    const AxisSpec axisCentFT0M{101, 0, 101, "FT0M (%)"};
    const AxisSpec axisCent{centBins, "FT0M (%)"};
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
    const AxisSpec axisPid{300, 0, 3000, "PID"};

    HistogramConfigSpec partHistCent({HistType::kTHnSparseD, {axisCentFT0M, axisPart}});
    HistogramConfigSpec partMCHist({HistType::kTHnSparseD, {axisCentFT0M, axisPart}});
    HistogramConfigSpec tofNSigmaHist({HistType::kTH2D, {axisP, axisTOFNsigma}});
    HistogramConfigSpec tofSignalHist({HistType::kTH2D, {axisP, axisTOFSignal}});
    HistogramConfigSpec tpcNSigmaHist({HistType::kTH2D, {axisP, axisTPCNsigma}});
    HistogramConfigSpec tpcSignalHist({HistType::kTH2D, {axisP, axisTPCSignal}});
    HistogramConfigSpec tpcTofHist({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});
    HistogramConfigSpec pvsM2Hist({HistType::kTH2D, {axisM2, axisP}});
    HistogramConfigSpec tpcSignalHist1({HistType::kTH2D, {axisInnerParam, axisTPCSignal}});
    HistogramConfigSpec pvsM2Hist1({HistType::kTH2D, {axisM2, axisInnerParam}});

    // QA Plots
    hist.add("QA/before/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("QA/before/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});
    hist.add("QA/before/h_NTPC", "N_{TPC}", kTH1D, {axisMultTPC});
    hist.add("QA/before/h_NFT0M", "FT0M Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/before/h_CentM", "FT0M (%)", kTH1D, {axisCentFT0M});

    hist.add("QA/before/h2_TPCSignal", "TPC Signal", tpcSignalHist);
    hist.add("QA/before/h2_TOFSignal", "TOF Signal", tofSignalHist);
    hist.add("QA/before/h2_pvsm2", "p vs m^{2}", pvsM2Hist);

    hist.addClone("QA/before/", "QA/after/");

    hist.add("QA/before/h_Pt", "p_{T}", kTH1D, {axisPt});
    hist.add("QA/before/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/before/h_Phi", "#phi ", kTH1D, {axisPhi});
    hist.add("QA/before/h_DcaZ", "DCA_{Z}", kTH1D, {axisDCAz});
    hist.add("QA/before/h_DcaXY", "DCA_{XY}", kTH1D, {axisDCAxy});
    hist.add("QA/before/h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/before/h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});

    hist.add("QA/after/h_Ncorr", "N_{corr}", kTH1D, {axisMultCorr});
    hist.add("QA/after/h_counts_evSelCuts", "Event selection cuts", kTH1D, {axisEvents});
    hist.add("QA/after/h_TPCChi2perCluster", "TPC #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/after/h_ITSChi2perCluster", "ITS #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/after/h_crossedTPC", "Crossed TPC", kTH1D, {axisCrossedTPC});
    hist.add("QA/after/h2_NTPC_CentM", "N_{TPC} vs FT0M(%)", kTH2D, {{axisCentFT0M}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NFT0M", "N_{TPC} vs N_{FT0M}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});
    hist.add("QA/after/p_NTPC_NFT0M", "N_{TPC} vs N_{FT0M} (Profile)", kTProfile, {axisMultFT0M});
    hist.add("QA/after/p_NTPC_CentM", "N_{TPC} vs FT0M(%) (Profile)", kTProfile, {axisCentFT0M});
    hist.add("QA/after/h_DCAxy_primary", "DCA_{XY} (Primary)", kTH1D, {axisDCAxy});
    hist.add("QA/after/h_DCAz_primary", "DCA_{Z} (Primary)", kTH1D, {axisDCAz});
    hist.add("QA/after/h_DCAxy_secondary", "DCA_{XY} (Secondary)", kTH1D, {axisDCAxy});
    hist.add("QA/after/h_DCAz_secondary", "DCA_{Z} (Secondary)", kTH1D, {axisDCAz});
    hist.add("QA/after/innerParam/h2_TPCSignal", "TPC Signal", tpcSignalHist1);
    hist.add("QA/after/h2_NSim_NTPC", "Reco vs Truth Multiplicty TPC", kTH2D, {{axisMultMC}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NSim", "Truth vs Reco Multiplicty TPC", kTH2D, {{axisMultTPC}, {axisMultMC}});
    hist.add("QA/after/p_NTPC_NSim", "Truth vs Reco Multiplicty TPC", kTProfile, {axisMultTPC});
    hist.add("QA/after/h2_Ncorr_NSim", "Truth vs Corrected Multiplicty TPC", kTH2D, {{axisMultCorr}, {axisMultMC}});
    hist.add("QA/after/h2_Ncorr_CentM", "N_{corr} vs FT0M(%)", kTH2D, {{axisCentFT0M}, {axisMultCorr}});

    hist.add("QA/Charged/h_Pt", "p_{T}", kTH1D, {axisPt});
    hist.add("QA/Charged/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/Charged/h_Phi", "#phi ", kTH1D, {axisPhi});
    hist.add("QA/Charged/h_DcaZ", "DCA_{Z}", kTH1D, {axisDCAz});
    hist.add("QA/Charged/h_DcaXY", "DCA_{XY}", kTH1D, {axisDCAxy});
    hist.add("QA/Charged/h2_Pt_PtMC", "p_{T} vs p_{T} (MC)", kTH2D, {{axisPt}, {axisPt}});
    hist.add("QA/Charged/h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/Charged/h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/Charged/h2_Pt_centFT0M", "p_{T} in centrality Classes ", kTH2D, {{axisCentFT0M}, {axisPt}});
    hist.add("QA/Charged/h3_Pt_Eta_Phi", "p_{T}, #eta, #phi ", kTHnSparseD, {{axisPt}, {axisEta}, {axisPhi}});

    hist.addClone("QA/Charged/", "QA/Pion/");

    hist.add("QA/Pion/before/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/before/h2_TPCNsigma_nottof", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/before/h2_TPCNsigma_tof", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/before/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist);
    hist.add("QA/Pion/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist);
    hist.add("QA/Pion/h_Rap", "y ", kTH1D, {axisY});
    hist.add("QA/Pion/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist);
    hist.add("QA/Pion/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist);
    hist.add("QA/Pion/h2_TPCSignal", "TPC Signal ", tpcSignalHist);
    hist.add("QA/Pion/h2_TOFSignal", "TOF Signal", tofSignalHist);
    hist.add("QA/Pion/innerParam/h2_TPCSignal", "TPC Signal", tpcSignalHist1);
    hist.add("QA/Pion/h2_pvsm2", "p vs m^{2}", pvsM2Hist);
    hist.addClone("QA/Pion/", "QA/Kaon/");
    hist.addClone("QA/Pion/", "QA/Proton/");

    hist.add("QA/Charged/h_PtMC", "p_{T} (MC)", kTH1D, {axisPt});
    hist.add("QA/Charged/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/Charged/h3_Pt_Eta_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPt}, {axisCent}});
    hist.add("QA/Charged/h3_Pt_EtaMC_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPt}, {axisCent}});
    hist.add("QA/Charged/h2_Pt_EtaMC", "p_{T} vs #eta (MC)", kTH2D, {{axisEta}, {axisPt}});

    hist.add("QA/Pion/h_Pt_e", "p_{T}", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h_PtMC", "p_{T} (MC)", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPtPi}});
    hist.add("QA/Pion/h3_Pt_Eta_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtPi}, {axisCent}});
    hist.add("QA/Pion/h3_Pt_EtaMC_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtPi}, {axisCent}});
    hist.add("QA/Pion/h2_Pt_EtaMC", "p_{T} vs #eta (MC)", kTH2D, {{axisEta}, {axisPtPi}});
    hist.add("QA/Pion/h_PtPos", "p_{T} (positive) ", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h_PtNeg", "p_{T} (negative) ", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h_PtTruth", "p_{T} (Truth)", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h_PtPosTruth", "p_{T} (positive) (Truth)", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h_PtNegTruth", "p_{T} (negative) (Truth) ", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h_PtTruth_primary", "p_{T} (Truth Primary)", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h_PtTruth_secondary", "p_{T} (Truth Secondary)", kTH1D, {axisPtPi});
    hist.add("QA/Pion/h2_Pt_EtaTruth", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPtPi}});

    hist.add("QA/Kaon/h_Pt_e", "p_{T}", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h_PtMC", "p_{T} (MC)", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPtKa}});
    hist.add("QA/Kaon/h3_Pt_Eta_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtKa}, {axisCent}});
    hist.add("QA/Kaon/h3_Pt_EtaMC_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtKa}, {axisCent}});
    hist.add("QA/Kaon/h2_Pt_EtaMC", "p_{T} vs #eta (MC)", kTH2D, {{axisEta}, {axisPtKa}});
    hist.add("QA/Kaon/h_PtPos", "p_{T} (positive) ", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h_PtNeg", "p_{T} (negative) ", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h_PtTruth", "p_{T} (Truth)", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h_PtPosTruth", "p_{T} (positive) (Truth)", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h_PtNegTruth", "p_{T} (negative) (Truth) ", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h_PtTruth_primary", "p_{T} (Truth Primary)", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h_PtTruth_secondary", "p_{T} (Truth Secondary)", kTH1D, {axisPtKa});
    hist.add("QA/Kaon/h2_Pt_EtaTruth", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPtKa}});

    hist.add("QA/Proton/h_Pt_e", "p_{T}", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h_PtMC", "p_{T} (MC)", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPtPr}});
    hist.add("QA/Proton/h3_Pt_Eta_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtPr}, {axisCent}});
    hist.add("QA/Proton/h3_Pt_EtaMC_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtPr}, {axisCent}});
    hist.add("QA/Proton/h2_Pt_EtaMC", "p_{T} vs #eta (MC)", kTH2D, {{axisEta}, {axisPtPr}});
    hist.add("QA/Proton/h_PtPos", "p_{T} (positive) ", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h_PtNeg", "p_{T} (negative) ", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h_PtTruth", "p_{T} (Truth)", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h_PtPosTruth", "p_{T} (positive) (Truth)", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h_PtNegTruth", "p_{T} (negative) (Truth) ", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h_PtTruth_primary", "p_{T} (Truth Primary)", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h_PtTruth_secondary", "p_{T} (Truth Secondary)", kTH1D, {axisPtPr});
    hist.add("QA/Proton/h2_Pt_EtaTruth", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPtPr}});

    // AnalysisPlots
    hist.add("Analysis/Charged/h_mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Analysis/Charged/p_Q1_Mult", " Q1 vs Mult ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/p_Q2_Mult", " Q2 vs Mult ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/p_mean_pT_Mult", " <p_{T}> vs Mult ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/p_mean_pT_Cent", " <p_{T}> vs CentFT0M ", kTProfile, {axisCentFT0M});
    hist.add("Analysis/Charged/p_mean_pT_Cent_var", " <p_{T}> vs CentFT0M ", kTProfile, {axisCentFT0M});
    hist.add("Analysis/Charged/p_twopart_Cent_var", "Twopart vs CentFT0M ", kTProfile, {axisCentFT0M});
    hist.add("Analysis/Charged/h_mean_pT_Cent", " <p_{T}> vs CentFT0M ", kTH1D, {axisCentFT0M});
    hist.add("Analysis/Charged/h_mean_pT_Cent_var", " <p_{T}> vs CentFT0M ", kTH1D, {axisCentFT0M});
    hist.add("Analysis/Charged/h_twopart_Cent_var", "Twopart vs CentFT0M ", kTH1D, {axisCentFT0M});
    hist.add("Analysis/Charged/p_mean_pT_Num_Cent_var", " <p_{T}> Numerator vs CentFT0M ", kTProfile, {axisCentFT0M});
    hist.add("Analysis/Charged/p_mean_pT_Den_Cent_var", "<p_{T}> Denominator vs CentFT0M ", kTProfile, {axisCentFT0M});
    hist.add("Analysis/Charged/p_twopart_Num_Cent_var", "Twopart Numerator vs CentFT0M ", kTProfile, {axisCentFT0M});
    hist.add("Analysis/Charged/p_twopart_Den_Cent_var", "Twopart Denominator vs CentFT0M ", kTProfile, {axisCentFT0M});
    hist.add("Analysis/Charged/h_mean_pT_Num_Cent_var", " <p_{T}> Numerator vs CentFT0M ", kTH1D, {axisCentFT0M});
    hist.add("Analysis/Charged/h_mean_pT_Den_Cent_var", "<p_{T}> Denominator vs CentFT0M ", kTH1D, {axisCentFT0M});
    hist.add("Analysis/Charged/h_twopart_Num_Cent_var", "Twopart Numerator vs CentFT0M ", kTH1D, {axisCentFT0M});
    hist.add("Analysis/Charged/h_twopart_Den_Cent_var", "Twopart Denominator vs CentFT0M ", kTH1D, {axisCentFT0M});

    hist.addClone("Analysis/Charged/", "Analysis/Pion/");
    hist.addClone("Analysis/Charged/", "Analysis/Kaon/");
    hist.addClone("Analysis/Charged/", "Analysis/Proton/");

    hist.add("Analysis/Charged/hN1Matrix", "", kTH2D, {{axisCent}, {axisKCh}});
    hist.add("Analysis/Charged/hPt1Matrix", "", kTH2D, {{axisCent}, {axisKCh}});
    hist.add("Analysis/Charged/hN2Matrix", "", kTH3D, {{axisCent}, {axisKCh}, {axisKCh}});
    hist.add("Analysis/Charged/hPtPtMatrix", "", kTH3D, {{axisCent}, {axisKCh}, {axisKCh}});

    hist.add("Analysis/Pion/hN1Matrix", "", kTH2D, {{axisCent}, {axisKPi}});
    hist.add("Analysis/Pion/hPt1Matrix", "", kTH2D, {{axisCent}, {axisKPi}});
    hist.add("Analysis/Pion/hN2Matrix", "", kTH3D, {{axisCent}, {axisKPi}, {axisKPi}});
    hist.add("Analysis/Pion/hPtPtMatrix", "", kTH3D, {{axisCent}, {axisKPi}, {axisKPi}});

    hist.add("Analysis/Kaon/hN1Matrix", "", kTH2D, {{axisCent}, {axisKKa}});
    hist.add("Analysis/Kaon/hPt1Matrix", "", kTH2D, {{axisCent}, {axisKKa}});
    hist.add("Analysis/Kaon/hN2Matrix", "", kTH3D, {{axisCent}, {axisKKa}, {axisKKa}});
    hist.add("Analysis/Kaon/hPtPtMatrix", "", kTH3D, {{axisCent}, {axisKKa}, {axisKKa}});

    hist.add("Analysis/Proton/hN1Matrix", "", kTH2D, {{axisCent}, {axisKPr}});
    hist.add("Analysis/Proton/hPt1Matrix", "", kTH2D, {{axisCent}, {axisKPr}});
    hist.add("Analysis/Proton/hN2Matrix", "", kTH3D, {{axisCent}, {axisKPr}, {axisKPr}});
    hist.add("Analysis/Proton/hPtPtMatrix", "", kTH3D, {{axisCent}, {axisKPr}, {axisKPr}});

    // MC Generated
    hist.add("Gen/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("Gen/h_VtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/h_VtxZ_b", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/h_NSim", "Truth Multiplicity TPC", kTH1D, {axisMultMC});
    hist.add("Gen/h2_NSim_NTPC", "Reco vs Truth Multiplicty TPC", kTH2D, {{axisMultMC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NSim", "Truth vs Reco Multiplicty TPC", kTH2D, {{axisMultTPC}, {axisMultMC}});

    hist.add("Gen/Charged/h_EtaTruth", "#eta ", kTH1D, {axisEta});
    hist.add("Gen/Charged/h_PhiTruth", "#phi ", kTH1D, {axisPhi});
    hist.add("Gen/Charged/h2_PtTruth_centFT0M", "p_{T} in centrality Classes ", kTH2D, {{axisCentFT0M}, {axisPt}});
    hist.add("Gen/Charged/h3_Pt_Eta_PhiTruth", "p_{T}, #eta, #phi ", kTHnSparseD, {{axisPt}, {axisEta}, {axisPhi}});

    hist.addClone("Gen/Charged/", "Gen/Pion/");

    hist.add("Gen/Pion/h_RapTruth", "y", kTH1D, {axisY});

    hist.addClone("Gen/Pion/", "Gen/Kaon/");
    hist.addClone("Gen/Pion/", "Gen/Proton/");

    hist.addClone("Analysis/Charged/", "Gen/Charged/");
    hist.addClone("Analysis/Pion/", "Gen/Pion/");
    hist.addClone("Analysis/Kaon/", "Gen/Kaon/");
    hist.addClone("Analysis/Proton/", "Gen/Proton/");

    hist.add("Gen/Charged/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("Gen/Charged/h2_Pt_EtaTruth", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPt}});
    hist.add("Gen/Charged/h3_Pt_EtaTruth_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPt}, {axisCent}});

    hist.add("Gen/Pion/h_PtTruth", "p_{T} ", kTH1D, {axisPtPi});
    hist.add("Gen/Pion/h_PtPosTruth", "p_{T} (positive) ", kTH1D, {axisPtPi});
    hist.add("Gen/Pion/h_PtNegTruth", "p_{T} (negative) ", kTH1D, {axisPtPi});
    hist.add("Gen/Pion/h2_Pt_EtaTruth", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPtPi}});
    hist.add("Gen/Pion/h3_Pt_EtaTruth_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtPi}, {axisCent}});

    hist.add("Gen/Kaon/h_PtTruth", "p_{T} ", kTH1D, {axisPtKa});
    hist.add("Gen/Kaon/h_PtPosTruth", "p_{T} (positive) ", kTH1D, {axisPtKa});
    hist.add("Gen/Kaon/h_PtNegTruth", "p_{T} (negative) ", kTH1D, {axisPtKa});
    hist.add("Gen/Kaon/h2_Pt_EtaTruth", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPtKa}});
    hist.add("Gen/Kaon/h3_Pt_EtaTruth_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtKa}, {axisCent}});

    hist.add("Gen/Proton/h_PtTruth", "p_{T} ", kTH1D, {axisPtPr});
    hist.add("Gen/Proton/h_PtPosTruth", "p_{T} (positive) ", kTH1D, {axisPtPr});
    hist.add("Gen/Proton/h_PtNegTruth", "p_{T} (negative) ", kTH1D, {axisPtPr});
    hist.add("Gen/Proton/h2_Pt_EtaTruth", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPtPr}});
    hist.add("Gen/Proton/h3_Pt_EtaTruth_centFT0M", "p_{T} vs #eta in centrality classes", kTH3D, {{axisEta}, {axisPtPr}, {axisCent}});

    hist.add("QA/h_collisions_info", "Collisions info", kTH1D, {axisCol});
    hist.add("Gen/h_collisions_info", "Collisions info", kTH1D, {axisCol});
    hist.add("Gen/h_collision_recgen", "Number of Collisions ", kTH1D, {axisCol});
    hist.add("Gen/h2_collision_posZ", "Reco vs truth posZ ", kTH2D, {{axisVtxZ}, {axisVtxZ}});
    hist.add("Tracks/h_tracks_info", "Track info", kTH1D, {axisTrack});
    hist.add("Tracks/h2_tracks_pid_before_sel", "Track pid info before selection", kTH2D, {{axisPid}, {axisPt}});

    hist.get<TH1>(HIST("QA/h_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotCol, "kTotCol");
    hist.get<TH1>(HIST("QA/h_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kPassSelCol, "kPassSelCol");
    hist.get<TH1>(HIST("Gen/h_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotCol, "kTotCol");
    hist.get<TH1>(HIST("Gen/h_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kPassSelCol, "kPassSelCol");
    hist.get<TH1>(HIST("Tracks/h_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kTracksBeforeHasMcParticle, "kTracksBeforeHasMcParticle");
    hist.get<TH1>(HIST("Tracks/h_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kAllTracks, "kAllTracks");
    hist.get<TH1>(HIST("Tracks/h_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kAllSelPassed, "kAllSelPassed");
    hist.get<TH1>(HIST("QA/after/h_counts_evSelCuts"))->GetXaxis()->SetBinLabel(SelCollisionLabels::kBeforeSelCol, "kBeforeSelCol");
    hist.get<TH1>(HIST("QA/after/h_counts_evSelCuts"))->GetXaxis()->SetBinLabel(SelCollisionLabels::kSelColPosZ, "kSelColPosZ");
    hist.get<TH1>(HIST("QA/after/h_counts_evSelCuts"))->GetXaxis()->SetBinLabel(SelCollisionLabels::kSelColSel8, "kSelColSel8");
    hist.get<TH1>(HIST("QA/after/h_counts_evSelCuts"))->GetXaxis()->SetBinLabel(SelCollisionLabels::kSelColNoSameBunchPileup, "kSelColNoSameBunchPileup");
    hist.get<TH1>(HIST("QA/after/h_counts_evSelCuts"))->GetXaxis()->SetBinLabel(SelCollisionLabels::kSelColIsVertexITSTPC, "kSelColIsVertexITSTPC");
    if (cfgLightIonCuts) {
      hist.get<TH1>(HIST("QA/after/h_counts_evSelCuts"))->GetXaxis()->SetBinLabel(SelCollisionLabels::kSelColIsGoodITSLayersAll, "kSelColIsGoodITSLayersAll");
      hist.get<TH1>(HIST("QA/after/h_counts_evSelCuts"))->GetXaxis()->SetBinLabel(SelCollisionLabels::kSelColIsGoodZvtxFT0vsPV, "kSelColIsGoodZvtxFT0vsPV");
    }
  }

  float centFT0M = 0.;
  int nTPC = 0, nFT0M = 0;
  int nSim = 0;

  // Event selection cuts:
  template <typename T>
  bool selRun3Col(T const& col)
  {

    centFT0M = col.centFT0M();
    nTPC = col.multNTracksHasTPC();
    nFT0M = col.multFT0M();
    hist.fill(HIST("QA/after/h_counts_evSelCuts"), kBeforeSelCol);

    if (cfgPosZ) {
      if (std::abs(col.posZ()) > cfgCutPosZ) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), kSelColPosZ);
    }

    if (cfgSel8) {
      if (!col.sel8()) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), kSelColSel8);
    }
    if (cfgNoSameBunchPileup) {
      if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), kSelColNoSameBunchPileup);
    }

    if (cfgIsVertexITSTPC) {
      if (!col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), kSelColIsVertexITSTPC);
    }
    if (cfgLightIonCuts) {
      if (cfgIsGoodITSLayersAll) {
        if (!col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
          return false;
        }
        hist.fill(HIST("QA/after/h_counts_evSelCuts"), kSelColIsGoodITSLayersAll);
      }

      if (cfgIsGoodZvtxFT0vsPV) {
        if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
          return false;
        }
        hist.fill(HIST("QA/after/h_counts_evSelCuts"), kSelColIsGoodZvtxFT0vsPV);
      }
    }

    return true;
  }

  // Track selection cuts:
  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrack())
      return false;

    if (track.pt() < cfgCutPtMin)
      return false;

    if (track.pt() >= cfgCutPtMax)
      return false;

    if (track.sign() == 0)
      return false;

    if (std::fabs(track.dcaZ()) > cfgCutDcaZ)
      return false;

    if (std::fabs(track.dcaZ()) > (0.0105 + 0.035 / std::pow(track.p(), 1.1)))
      return false;

    if (std::abs(track.eta()) >= cfgCutEta)
      return false;

    return true;
  }

  // Rejection cuts to reject the tracks
  template <typename T>
  bool rejectTracks(T const& track)
  {
    if (((track.tpcNSigmaEl()) > -cfgCutNSig3 &&
         (track.tpcNSigmaEl()) < cfgCutNSig5) &&
        (std::fabs(track.tpcNSigmaPi()) > cfgCutNSig3 &&
         std::fabs(track.tpcNSigmaKa()) > cfgCutNSig3 &&
         std::fabs(track.tpcNSigmaPr()) > cfgCutNSig3)) {
      return true;
    }

    return false;
  }

  // PID cuts for identified particles
  template <typename T>
  bool identifyParticle(T const& track, PIDType species, PIDType reject1, PIDType reject2, float p, float momThreshold)
  {
    std::vector<float> vTpcNSigma = {-999.f, track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};

    std::vector<float> vTofNSigma = {-999.f, track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};

    const int sp = static_cast<int>(species);
    const int sq = static_cast<int>(reject1);
    const int sr = static_cast<int>(reject2);

    bool isTofPidFlag = false;
    bool isTpcPidFlag = false;

    float nSigmaSelCut = 0.f;

    switch (species) {
      case PIDType::kPions:
        nSigmaSelCut = cfgSelCutNSigPi;
        break;

      case PIDType::kKaons:
        nSigmaSelCut = cfgSelCutNSigKa;
        break;

      case PIDType::kProtons:
        nSigmaSelCut = cfgSelCutNSigPr;
        break;

      default:
        return false;
    }

    if (track.hasTOF()) {

      if (std::abs(vTofNSigma[sp]) < nSigmaSelCut &&
          std::abs(vTofNSigma[sq]) > std::abs(vTofNSigma[sp]) &&
          std::abs(vTofNSigma[sr]) > std::abs(vTofNSigma[sp])) {
        isTofPidFlag = true;
      }

      if (std::abs(vTpcNSigma[sp]) < cfgCutNSig2) {
        isTpcPidFlag = true;
      }

    } else {

      if (p >= momThreshold) {
        return false;
      }

      if (std::abs(vTpcNSigma[sp]) < nSigmaSelCut &&
          std::abs(vTpcNSigma[sq]) > std::abs(vTpcNSigma[sp]) &&
          std::abs(vTpcNSigma[sr]) > std::abs(vTpcNSigma[sp])) {
        isTofPidFlag = true;
        isTpcPidFlag = true;
      }
    }

    return (isTofPidFlag && isTpcPidFlag);
  }

  // Fill QA histograms before selection cuts:
  template <typename T, typename U>
  void fillBeforeQAHistos(T const& col, U const& tracks)
  {
    for (const auto& track : tracks) {
      hist.fill(HIST("QA/before/h_Eta"), track.eta());
      hist.fill(HIST("QA/before/h_Phi"), track.phi());
      hist.fill(HIST("QA/before/h_Pt"), track.pt());
      hist.fill(HIST("QA/before/h_DcaXY"), track.dcaXY());
      hist.fill(HIST("QA/before/h_DcaZ"), track.dcaZ());
      hist.fill(HIST("QA/before/h2_DcaXY"), track.pt(), track.dcaXY());
      hist.fill(HIST("QA/before/h2_DcaZ"), track.pt(), track.dcaZ());
    }
    hist.fill(HIST("QA/before/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/before/h_Counts"), 2);
    hist.fill(HIST("QA/before/h_NTPC"), col.multNTracksHasTPC());
    hist.fill(HIST("QA/before/h_CentM"), col.centFT0M());
    hist.fill(HIST("QA/before/h_NFT0M"), col.multFT0M());
  }

  // Fill QA histograms after selection cuts:
  template <typename T>
  void fillAfterQAHistos(T const& col)
  {
    hist.fill(HIST("QA/after/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/after/h_Counts"), 2);
    hist.fill(HIST("QA/after/h_NTPC"), nTPC);
    hist.fill(HIST("QA/after/h_CentM"), centFT0M);
    hist.fill(HIST("QA/after/h_NFT0M"), nFT0M);
    hist.fill(HIST("QA/after/h2_NTPC_NFT0M"), nFT0M, nTPC);
    hist.fill(HIST("QA/after/h2_NTPC_CentM"), centFT0M, nTPC);
    hist.fill(HIST("QA/after/p_NTPC_CentM"), centFT0M, nTPC);
    hist.fill(HIST("QA/after/p_NTPC_NFT0M"), nFT0M, nTPC);
    hist.fill(HIST("QA/after/h2_NSim_NTPC"), nSim, nTPC);
    hist.fill(HIST("QA/after/h2_NTPC_NSim"), nTPC, nSim);
    hist.fill(HIST("QA/after/p_NTPC_NSim"), nTPC, nSim);
  }

  // Fill Charged particles QA histograms after selection cuts:
  template <typename T>
  void fillChargedQAHistos(T const& track, float centFT0M)
  {
    hist.fill(HIST("QA/Charged/h_Eta"), track.eta());
    hist.fill(HIST("QA/Charged/h_Phi"), track.phi());
    hist.fill(HIST("QA/Charged/h2_Pt_centFT0M"), centFT0M, track.pt());
    hist.fill(HIST("QA/Charged/h3_Pt_Eta_Phi"), track.pt(), track.eta(), track.phi());
    hist.fill(HIST("QA/Charged/h3_Pt_Eta_centFT0M"), track.eta(), track.pt(), centFT0M);
    hist.fill(HIST("QA/Charged/h2_Pt_Eta"), track.eta(), track.pt());
    hist.fill(HIST("QA/Charged/h_DcaZ"), track.dcaZ());
    hist.fill(HIST("QA/Charged/h_DcaXY"), track.dcaXY());
    hist.fill(HIST("QA/Charged/h2_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST("QA/Charged/h2_DcaZ"), track.pt(), track.dcaZ());

    hist.fill(HIST("QA/after/h_TPCChi2perCluster"), track.tpcChi2NCl());
    hist.fill(HIST("QA/after/h_ITSChi2perCluster"), track.itsChi2NCl());
    hist.fill(HIST("QA/after/h_crossedTPC"), track.tpcNClsCrossedRows());
  }

  // Fill identified QA histograms before PID cuts:
  template <typename T>
  void fillBeforePIDQAHistos(T const& track)
  {
    hist.fill(HIST("QA/before/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/before/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/before/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST("QA/Pion/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPi());
    if (!track.hasTOF())
      hist.fill(HIST("QA/Pion/before/h2_TPCNsigma_nottof"), track.p(), track.tpcNSigmaPi());
    if (track.hasTOF())
      hist.fill(HIST("QA/Pion/before/h2_TPCNsigma_tof"), track.p(), track.tpcNSigmaPi());
    hist.fill(HIST("QA/Pion/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/before/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());

    hist.fill(HIST("QA/Kaon/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaKa());
    if (!track.hasTOF())
      hist.fill(HIST("QA/Kaon/before/h2_TPCNsigma_nottof"), track.p(), track.tpcNSigmaKa());
    if (track.hasTOF())
      hist.fill(HIST("QA/Kaon/before/h2_TPCNsigma_tof"), track.p(), track.tpcNSigmaKa());
    hist.fill(HIST("QA/Kaon/before/h2_TOFNsigma"), track.p(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/before/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());

    hist.fill(HIST("QA/Proton/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPr());
    if (!track.hasTOF())
      hist.fill(HIST("QA/Proton/before/h2_TPCNsigma_nottof"), track.p(), track.tpcNSigmaPr());
    if (track.hasTOF())
      hist.fill(HIST("QA/Proton/before/h2_TPCNsigma_tof"), track.p(), track.tpcNSigmaPr());
    hist.fill(HIST("QA/Proton/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/before/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
  }

  // Fill identified QA histograms after PID cuts:
  template <int Mode, typename T>
  void fillIdParticleQAHistos(T const& track, float ptMC, float etaMC, float rap, float nSigmaTPC, float nSigmaTOF, std::vector<TrackInfo>& Tracks)
  {

    float pt = track.pt();
    float eta = track.eta();
    float phi = track.phi();

    Tracks.push_back({pt, eta});

    if (track.sign() > 0) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtPos"), pt);
    }
    if (track.sign() < 0) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtNeg"), pt);
    }

    hist.fill(HIST(Dire[Mode]) + HIST("h_Pt"), pt);
    hist.fill(HIST(Dire[Mode]) + HIST("h_PtMC"), ptMC);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_PtMC"), ptMC, pt);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_centFT0M"), centFT0M, pt);
    hist.fill(HIST(Dire[Mode]) + HIST("h3_Pt_Eta_Phi"), pt, eta, phi);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Eta"), eta, pt);
    hist.fill(HIST(Dire[Mode]) + HIST("h3_Pt_Eta_centFT0M"), eta, pt, centFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_EtaMC"), etaMC, ptMC);
    hist.fill(HIST(Dire[Mode]) + HIST("h3_Pt_EtaMC_centFT0M"), etaMC, ptMC, centFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Eta"), eta);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Phi"), phi);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Rap"), rap);
    hist.fill(HIST(Dire[Mode]) + HIST("h_DcaZ"), track.dcaZ());
    hist.fill(HIST(Dire[Mode]) + HIST("h_DcaXY"), track.dcaXY());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_DcaZ"), pt, track.dcaZ());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_DcaXY"), pt, track.dcaXY());

    hist.fill(HIST(Dire[Mode]) + HIST("h2_TPCNsigma"), track.p(), nSigmaTPC);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TOFNsigma"), track.p(), nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TpcTofNsigma"), nSigmaTPC, nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/after/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/after/h2_pvsm2"), track.mass() * track.mass(), track.p());
  }

  // Fill Truth identified particles histograms
  template <int Mode>
  void fillPtMCHist(bool cfgGen, float pt, float eta, float rap, float phi, int pid, int pdgCodePos, int pdgCodeNeg)
  {
    if (cfgGen) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtTruth"), pt);
      hist.fill(HIST(Dire[Mode]) + HIST("h_EtaTruth"), eta);
      hist.fill(HIST(Dire[Mode]) + HIST("h_RapTruth"), rap);
      hist.fill(HIST(Dire[Mode]) + HIST("h_PhiTruth"), phi);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtTruth_centFT0M"), centFT0M, pt);
      hist.fill(HIST(Dire[Mode]) + HIST("h3_Pt_EtaTruth_centFT0M"), eta, pt, centFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h3_Pt_Eta_PhiTruth"), pt, eta, phi);
    }

    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_EtaTruth"), eta, pt);

    if (pid == pdgCodePos) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtPosTruth"), pt);
    }
    if (pid == pdgCodeNeg) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtNegTruth"), pt);
    }
  }

  int getPtBin(double pt, const std::vector<double>& ptBins)
  {
    for (size_t i = 0; i < ptBins.size() - 1; ++i) {
      if (pt >= ptBins[i] && pt < ptBins[i + 1]) {
        return i;
      }
    }
    return -1;
  }

  int getEtaBin(double eta)
  {
    for (size_t i = 0; i < etaBins->size() - 1; ++i) {
      if (eta >= etaBins->at(i) && eta < etaBins->at(i + 1)) {
        return i;
      }
    }
    return -1;
  }

  // fill differential matrices for two-particle analysis
  template <int Mode>
  void fillDifferentialMatrices(const std::vector<TrackInfo>& tracks, const std::vector<double>& ptB)
  {
    const int nPtBins = ptB.size() - 1;

    auto getK = [&](int ipt, int ieta) {
      return ieta * nPtBins + ipt;
    };

    for (const auto& trk : tracks) {
      int ipt = getPtBin(trk.pt, ptB);
      int ieta = getEtaBin(trk.eta);

      if (ipt < 0 || ieta < 0)
        continue;

      int k = getK(ipt, ieta);

      hist.fill(HIST(Dire[Mode]) + HIST("hN1Matrix"), centFT0M, k, 1.0);
      hist.fill(HIST(Dire[Mode]) + HIST("hPt1Matrix"), centFT0M, k, trk.pt);
    }

    for (size_t i = 0; i < tracks.size(); ++i) {
      int ipt1 = getPtBin(tracks[i].pt, ptB);
      int ieta1 = getEtaBin(tracks[i].eta);

      if (ipt1 < 0 || ieta1 < 0)
        continue;

      int k1 = getK(ipt1, ieta1);

      for (size_t j = i + 1; j < tracks.size(); ++j) {
        int ipt2 = getPtBin(tracks[j].pt, ptB);
        int ieta2 = getEtaBin(tracks[j].eta);

        if (ipt2 < 0 || ieta2 < 0)
          continue;

        int k2 = getK(ipt2, ieta2);

        hist.fill(HIST(Dire[Mode]) + HIST("hN2Matrix"), centFT0M, k1, k2, 1.0);
        hist.fill(HIST(Dire[Mode]) + HIST("hN2Matrix"), centFT0M, k2, k1, 1.0);
        hist.fill(HIST(Dire[Mode]) + HIST("hPtPtMatrix"), centFT0M, k1, k2, tracks[i].pt * tracks[j].pt);
        hist.fill(HIST(Dire[Mode]) + HIST("hPtPtMatrix"), centFT0M, k2, k1, tracks[j].pt * tracks[i].pt);
      }
    }
  }

  // Calculate moments for two-particle correlations
  void moments(double pt, int& n, double& q1, double& q2)
  {
    n++;
    q1 += pt;
    q2 += pt * pt;
  }

  // Fill terms calculted from moments method for two-particle correlations
  template <int Mode>
  void fillAnalysis(int n, double q1, double q2)
  {

    double npair = n * (n - 1);

    if (n > 0) {
      double meanPt = q1 / n;

      hist.fill(HIST(Dire[Mode]) + HIST("h_mult"), n);

      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Cent"), centFT0M, meanPt);
      hist.fill(HIST(Dire[Mode]) + HIST("p_mean_pT_Cent"), centFT0M, meanPt);

      hist.fill(HIST(Dire[Mode]) + HIST("p_mean_pT_Num_Cent_var"), centFT0M, q1);
      hist.fill(HIST(Dire[Mode]) + HIST("p_mean_pT_Den_Cent_var"), centFT0M, n);

      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Num_Cent_var"), centFT0M, q1);
      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Den_Cent_var"), centFT0M, n);

      hist.fill(HIST(Dire[Mode]) + HIST("p_Q1_Mult"), n, q1);
      hist.fill(HIST(Dire[Mode]) + HIST("p_Q2_Mult"), n, q2);
      hist.fill(HIST(Dire[Mode]) + HIST("p_mean_pT_Mult"), n, meanPt);
    }

    if (npair != 0) {
      double meanPt = q1 / n;
      double twopart = (q1 * q1 - q2);
      double varPt = twopart / npair;
      hist.fill(HIST(Dire[Mode]) + HIST("p_twopart_Num_Cent_var"), centFT0M, twopart);
      hist.fill(HIST(Dire[Mode]) + HIST("p_twopart_Den_Cent_var"), centFT0M, npair);

      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Num_Cent_var"), centFT0M, twopart);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Den_Cent_var"), centFT0M, npair);

      hist.fill(HIST(Dire[Mode]) + HIST("p_mean_pT_Cent_var"), centFT0M, meanPt);
      hist.fill(HIST(Dire[Mode]) + HIST("p_twopart_Cent_var"), centFT0M, varPt);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Cent_var"), centFT0M, varPt);
      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Cent_var"), centFT0M, meanPt);
    }
  }

  // Fill reconstructed histograms for charged and identified particles
  template <bool DataFlag, bool RecoFlag, typename C, typename T>
  void fillRecoHistos(C const& col, T const& tracks)
  {
    chargedTracks.clear();
    pionTracks.clear();
    kaonTracks.clear();
    protonTracks.clear();

    int nCh = 0, nPi = 0, nKa = 0, nPr = 0;
    float pt = 0., eta = 0., phi = 0., p = 0.;
    double q1Ch = 0., q2Ch = 0., q1Pi = 0., q2Pi = 0., q1Ka = 0., q2Ka = 0., q1Pr = 0., q2Pr = 0.;

    hist.fill(HIST("QA/h_collisions_info"), kTotCol);

    if constexpr (DataFlag) {
      if (!selRun3Col(col)) {
        return;
      }
    }

    if constexpr (RecoFlag) {
      // float lower1 =  cfgP1L1 * nSim + cfgP0L1;
      float lower2 = cfgP2L * nSim * nSim + cfgP1L * nSim - cfgP0L;
      float upper = cfgP2U * nSim * nSim + cfgP1U * nSim + cfgP0U;
      if (nTPC < lower2 || nTPC > upper) {
        return; // Reject event
      }
    }

    hist.fill(HIST("QA/h_collisions_info"), kPassSelCol);

    fillAfterQAHistos(col);

    double nCorr = cfgP2corr * nTPC * nTPC + cfgP1corr * nTPC + cfgP0corr;
    hist.fill(HIST("QA/after/h_Ncorr"), nCorr);
    hist.fill(HIST("QA/after/h2_Ncorr_NSim"), nCorr, nSim);
    hist.fill(HIST("QA/after/h2_Ncorr_CentM"), centFT0M, nCorr);

    for (auto const& track : tracks) {
      float nSigmaTPCPi = track.tpcNSigmaPi();
      float nSigmaTPCKa = track.tpcNSigmaKa();
      float nSigmaTPCPr = track.tpcNSigmaPr();
      float nSigmaTOFPi = track.tofNSigmaPi();
      float nSigmaTOFKa = track.tofNSigmaKa();
      float nSigmaTOFPr = track.tofNSigmaPr();
      float rapPi = track.rapidity(MassPiPlus);
      float rapKa = track.rapidity(MassKPlus);
      float rapPr = track.rapidity(MassProton);

      if constexpr (RecoFlag) {
        if (!track.has_mcParticle()) {
          hist.fill(HIST("Tracks/h_tracks_info"), kTracksBeforeHasMcParticle);
          continue;
        }
      }
      hist.fill(HIST("Tracks/h_tracks_info"), kAllTracks);

      if (!selTrack(track)) {
        continue;
      }
      hist.fill(HIST("Tracks/h_tracks_info"), kAllSelPassed);

      pt = track.pt();
      eta = track.eta();
      phi = track.phi();
      p = track.p();
      float ptMC = pt;
      float etaMC = eta;
      if constexpr (RecoFlag) {
        hist.fill(HIST("Tracks/h2_tracks_pid_before_sel"), track.mcParticle().pdgCode(), track.pt());

        auto mc = track.template mcParticle_as<aod::McParticles>();
        ptMC = mc.pt();
        etaMC = mc.eta();
        p = mc.p();

        if (mc.isPhysicalPrimary()) {
          hist.fill(HIST("QA/after/h_DCAxy_primary"), track.dcaXY());
          hist.fill(HIST("QA/after/h_DCAz_primary"), track.dcaZ());
        } else {
          hist.fill(HIST("QA/after/h_DCAxy_secondary"), track.dcaXY());
          hist.fill(HIST("QA/after/h_DCAz_secondary"), track.dcaZ());
        }
      }

      nCh++;
      chargedTracks.push_back({pt, eta});
      moments(pt, nCh, q1Ch, q2Ch);

      hist.fill(HIST("QA/Charged/h_Pt"), pt);
      hist.fill(HIST("QA/Charged/h_PtMC"), ptMC);
      hist.fill(HIST("QA/Charged/h2_Pt_PtMC"), ptMC, pt);
      hist.fill(HIST("QA/Charged/h2_Pt_EtaMC"), etaMC, ptMC);
      hist.fill(HIST("QA/Charged/h3_Pt_EtaMC_centFT0M"), etaMC, ptMC, centFT0M);

      fillChargedQAHistos(track, centFT0M);

      fillBeforePIDQAHistos(track);

      if (cfgRejTrk && rejectTracks(track)) {
        continue;
      }
      auto selIDPion = identifyParticle(track, PIDType::kPions, PIDType::kKaons, PIDType::kProtons, p, cfgCutPiThrsldP);
      auto selIDKaon = identifyParticle(track, PIDType::kKaons, PIDType::kPions, PIDType::kProtons, p, cfgCutKaThrsldP);
      auto selIDProton = identifyParticle(track, PIDType::kProtons, PIDType::kPions, PIDType::kKaons, p, cfgCutPrThrsldP);

      if (selIDPion && track.pt() >= cfgCutPiPtMin) {
        moments(pt, nPi, q1Pi, q2Pi);
        fillIdParticleQAHistos<QA_Pion>(track, ptMC, etaMC, rapPi, nSigmaTPCPi, nSigmaTOFPi, pionTracks);
      }
      if (selIDKaon && track.pt() >= cfgCutKaPtMin) {
        moments(pt, nKa, q1Ka, q2Ka);
        fillIdParticleQAHistos<QA_Kaon>(track, ptMC, etaMC, rapKa, nSigmaTPCKa, nSigmaTOFKa, kaonTracks);
      }
      if (selIDProton && track.pt() >= cfgCutPrPtMin) {
        moments(pt, nPr, q1Pr, q2Pr);
        fillIdParticleQAHistos<QA_Proton>(track, ptMC, etaMC, rapPr, nSigmaTPCPr, nSigmaTOFPr, protonTracks);
      }

      if constexpr (RecoFlag) {
        auto mc = track.template mcParticle_as<aod::McParticles>();
        int pid = mc.pdgCode();
        if (selIDPion && track.pt() >= cfgCutPiPtMin) {
          if (std::abs(pid) == kPiPlus) {
            hist.fill(HIST("QA/Pion/h_PtTruth"), track.pt());
            fillPtMCHist<QA_Pion>(false, pt, eta, rapPi, phi, pid, kPiPlus, kPiMinus);
            if (mc.isPhysicalPrimary()) {
              hist.fill(HIST("QA/Pion/h_PtTruth_primary"), track.pt());
            } else {
              hist.fill(HIST("QA/Pion/h_PtTruth_secondary"), track.pt());
            }
          }
        }
        if (selIDKaon && track.pt() >= cfgCutKaPtMin) {
          if (std::abs(pid) == kKPlus) {
            hist.fill(HIST("QA/Kaon/h_PtTruth"), track.pt());
            fillPtMCHist<QA_Kaon>(false, pt, eta, rapKa, phi, pid, kKPlus, kKMinus);
            if (mc.isPhysicalPrimary()) {
              hist.fill(HIST("QA/Kaon/h_PtTruth_primary"), track.pt());
            } else {
              hist.fill(HIST("QA/Kaon/h_PtTruth_secondary"), track.pt());
            }
          }
        }
        if (selIDProton && track.pt() >= cfgCutPrPtMin) {
          if (std::abs(pid) == kProton) {
            hist.fill(HIST("QA/Proton/h_PtTruth"), track.pt());
            fillPtMCHist<QA_Proton>(false, pt, eta, rapPr, phi, pid, kProton, kProtonBar);
            if (mc.isPhysicalPrimary()) {
              hist.fill(HIST("QA/Proton/h_PtTruth_primary"), track.pt());
            } else {
              hist.fill(HIST("QA/Proton/h_PtTruth_secondary"), track.pt());
            }
          }
        }
      }
    }
    fillDifferentialMatrices<Analysis_Charged>(chargedTracks, ptBins.value);
    fillDifferentialMatrices<Analysis_Pion>(pionTracks, ptBinsPi.value);
    fillDifferentialMatrices<Analysis_Kaon>(kaonTracks, ptBinsKa.value);
    fillDifferentialMatrices<Analysis_Proton>(protonTracks, ptBinsPr.value);
    fillAnalysis<Analysis_Charged>(nCh, q1Ch, q2Ch);
    fillAnalysis<Analysis_Pion>(nPi, q1Pi, q2Pi);
    fillAnalysis<Analysis_Kaon>(nKa, q1Ka, q2Ka);
    fillAnalysis<Analysis_Proton>(nPr, q1Pr, q2Pr);
  }

  template <typename C, typename M>
  void fillGenHistos(C const& mcCol, M const& mcParticles)
  {
    chargedTracksGen.clear();
    pionTracksGen.clear();
    kaonTracksGen.clear();
    protonTracksGen.clear();
    nSim = 0;
    int nChSim = 0, nPiSim = 0, nKaSim = 0, nPrSim = 0;
    float pt = 0., eta = 0, phi = 0., rap = 0.;
    double q1Ch = 0., q2Ch = 0., q1Pi = 0., q2Pi = 0., q1Ka = 0., q2Ka = 0., q1Pr = 0., q2Pr = 0.;
    for (auto const& mcPart : mcParticles) {
      if (!mcPart.isPhysicalPrimary()) {
        continue;
      }

      auto* particle = pdg->GetParticle(mcPart.pdgCode());

      if (!particle)
        continue;

      if (particle->Charge() == 0)
        continue;

      pt = mcPart.pt();
      eta = mcPart.eta();
      phi = mcPart.phi();

      if (std::abs(eta) < cfgCutEta) {
        nSim++;
      }

      auto pid = mcPart.pdgCode();
      if (std::abs(pid) != kElectron && std::abs(pid) != kMuonMinus && std::abs(pid) != kPiPlus && std::abs(pid) != kKPlus && std::abs(pid) != kProton) {
        continue;
      }

      if (pt >= cfgCutPtMin && pt < cfgCutPtMax && std::abs(eta) < cfgCutEta) {
        moments(pt, nChSim, q1Ch, q2Ch);
        chargedTracksGen.push_back({pt, eta});
        hist.fill(HIST("Gen/Charged/h_PtTruth"), pt);
        hist.fill(HIST("Gen/Charged/h_EtaTruth"), eta);
        hist.fill(HIST("Gen/Charged/h_PhiTruth"), phi);
        hist.fill(HIST("Gen/Charged/h2_Pt_EtaTruth"), eta, pt);
        hist.fill(HIST("Gen/Charged/h2_PtTruth_centFT0M"), centFT0M, pt);
        hist.fill(HIST("Gen/Charged/h3_Pt_EtaTruth_centFT0M"), eta, pt, centFT0M);
        hist.fill(HIST("Gen/Charged/h3_Pt_Eta_PhiTruth"), pt, eta, phi);

        rap = mcPart.y();
        if (std::abs(pid) == kPiPlus && pt >= cfgCutPiPtMin) {
          moments(pt, nPiSim, q1Pi, q2Pi);
          pionTracksGen.push_back({pt, eta});
          fillPtMCHist<Gen_Pion>(true, pt, eta, rap, phi, pid, kPiMinus, kPiMinus);
        }
        if (std::abs(pid) == kKPlus && pt >= cfgCutKaPtMin) {
          moments(pt, nKaSim, q1Ka, q2Ka);
          kaonTracksGen.push_back({pt, eta});
          fillPtMCHist<Gen_Kaon>(true, pt, eta, rap, phi, pid, kKMinus, kKMinus);
        }
        if (std::abs(pid) == kProton && pt >= cfgCutPrPtMin) {
          moments(pt, nPrSim, q1Pr, q2Pr);
          protonTracksGen.push_back({pt, eta});
          fillPtMCHist<Gen_Proton>(true, pt, eta, rap, phi, pid, kProtonBar, kProtonBar);
        }
      }
    }

    fillDifferentialMatrices<Gen_Charged>(chargedTracksGen, ptBins.value);
    fillDifferentialMatrices<Gen_Pion>(pionTracksGen, ptBinsPi.value);
    fillDifferentialMatrices<Gen_Kaon>(kaonTracksGen, ptBinsKa.value);
    fillDifferentialMatrices<Gen_Proton>(protonTracksGen, ptBinsPr.value);
    fillAnalysis<Gen_Charged>(nChSim, q1Ch, q2Ch);
    fillAnalysis<Gen_Pion>(nPiSim, q1Pi, q2Pi);
    fillAnalysis<Gen_Kaon>(nKaSim, q1Ka, q2Ka);
    fillAnalysis<Gen_Proton>(nPrSim, q1Pr, q2Pr);

    hist.fill(HIST("Gen/h_Counts"), 2);
    hist.fill(HIST("Gen/h_VtxZ"), mcCol.posZ());
    hist.fill(HIST("Gen/h_NSim"), nSim);
    hist.fill(HIST("Gen/h2_NSim_NTPC"), nSim, nTPC);
    hist.fill(HIST("Gen/h2_NTPC_NSim"), nTPC, nSim);
  }

  template <typename M, typename C, typename T, typename P>
  void analyzeMC(M const& mcCol, C const& cols, T const& tracks, P const& mcParts)
  {
    hist.fill(HIST("Gen/h_VtxZ_b"), mcCol.posZ());
    int nRecCols = cols.size();
    if (nRecCols == 0) {
      hist.fill(HIST("Gen/h_collision_recgen"), nRecCols);
    }
    // Do not analyze if more than one reco collision is accociated to one mc gen collision
    if (nRecCols != 1) {
      return;
    }
    hist.fill(HIST("Gen/h_collisions_info"), kTotCol);

    // Check the reco collision
    if (!cols.begin().has_mcCollision() || !selRun3Col(cols.begin()) || cols.begin().mcCollisionId() != mcCol.globalIndex()) {
      return;
    }

    hist.fill(HIST("Gen/h_collisions_info"), kPassSelCol);
    hist.fill(HIST("Gen/h2_collision_posZ"), mcCol.posZ(), cols.begin().posZ());

    auto sTracks = tracks.sliceBy(perCollision, cols.begin().globalIndex());
    fillGenHistos(mcCol, mcParts);
    fillRecoHistos<false, true>(cols.begin(), sTracks);
  }

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                                aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                                aod::pidTOFbeta, aod::pidTOFmass>;
  using MyRun3Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Ms, aod::CentFT0Ms>;
  using MyRun3MCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Ms, aod::CentFT0Ms, aod::McCollisionLabels>;
  using MyMCTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                               aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                               aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                               aod::pidTOFbeta, aod::pidTOFmass, aod::McTrackLabels>;
  SliceCache cache;
  Preslice<MyMCTracks> perCollision = aod::track::collisionId;

  void processRun3(MyRun3Collisions::iterator const& col, MyAllTracks const& tracks)
  {
    fillBeforeQAHistos(col, tracks);
    fillRecoHistos<true, false>(col, tracks);
  }
  PROCESS_SWITCH(MeanPtFlucId, processRun3, "Process for Run-3", false);

  void processMCRecoSimRun3(aod::McCollisions::iterator const& mcCol, soa::SmallGroups<MyRun3MCCollisions> const& cols, MyMCTracks const& tracks, aod::McParticles const& mcParts)
  {
    analyzeMC(mcCol, cols, tracks, mcParts);
  }
  PROCESS_SWITCH(MeanPtFlucId, processMCRecoSimRun3, "process MC Reconstructed & Truth Run-3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MeanPtFlucId>(cfgc)};
}
