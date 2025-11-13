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

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <string>
#include <utility>
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
  Configurable<float> cfgCutRap{"cfgCutRap", 0.5, "Rapidity Cut"};
  Configurable<float> cfgCutDcaZ{"cfgCutDcaZ", 0.15, "DCAz cut"};
  Configurable<float> cfgCutPosZ{"cfgCutPosZ", 7.0, "cut for vertex Z"};
  Configurable<bool> cfgPosZ{"cfgPosZ", true, "Position Z"};
  Configurable<float> cfgCutNSig2{"cfgCutNSig2", 2.0, "nSigma cut: 2"};
  Configurable<float> cfgCutNSig3{"cfgCutNSig3", 3.0, "nSigma cut: 3"};
  Configurable<float> cfgCutNSig5{"cfgCutNSig5", 5.0, "nSigma cut: 5"};
  Configurable<float> cfgSelCutNSigPi{"cfgSelCutNSigPi", 2.0, "nSigma cut for pion selection"};
  Configurable<float> cfgSelCutNSigKa{"cfgSelCutNSigKa", 3.0, "nSigma cut for kaon selection"};
  Configurable<float> cfgSelCutNSigPr{"cfgSelCutNSigPr", 2.0, "nSigma cut for proton selection"};
  Configurable<float> cfgRejCutNSigPi{"cfgRejCutNSigPi", 3.0, "nSigma cut for rejection of other particles while selecting pion"};
  Configurable<float> cfgCutPiPtMin{"cfgCutPiPtMin", 0.2, "Minimum pion p_{T} cut"};
  Configurable<float> cfgCutKaPtMin{"cfgCutKaPtMin", 0.3, "Minimum kaon p_{T} cut"};
  Configurable<float> cfgCutPrPtMin{"cfgCutPrPtMin", 0.5, "Minimum proton p_{T} cut"};
  Configurable<float> cfgCutPiThrsldP{"cfgCutPiThrsldP", 0.6, "Threshold p cut pion"};
  Configurable<float> cfgCutKaThrsldP{"cfgCutKaThrsldP", 0.6, "Threshold p cut kaon"};
  Configurable<float> cfgCutPrThrsldP{"cfgCutPrThrsldP", 1.0, "Threshold p cut proton "};
  Configurable<bool> cfgSel8{"cfgSel8", true, "Sel8 trigger"};
  Configurable<float> cfgMinWeight{"cfgMinWeight", 1e-6, "Minimum weight for efficiency correction"};
  Configurable<bool> cfgNoSameBunchPileup{"cfgNoSameBunchPileup", true, "kNoSameBunchPileup"};
  Configurable<bool> cfgIsVertexITSTPC{"cfgIsVertexITSTPC", true, "kIsVertexITSTPC"};
  Configurable<bool> cfgRejTrk{"cfgRejTrk", true, "Rejected Tracks"};
  Configurable<bool> cfgLoadEff{"cfgLoadEff", true, "Load efficiency"};
  Configurable<bool> cfgCorrection{"cfgCorrection", true, "Correction"};
  Configurable<bool> cfgWeightPtCh{"cfgWeightPtCh", true, "Efficiency correction (pT) for charged particles"};
  Configurable<bool> cfgWeightPtId{"cfgWeightPtId", false, "Efficiency correction (pT) "};
  Configurable<bool> cfgWeightPtEtaId{"cfgWeightPtEtaId", false, "Efficiency correction (pT) "};
  Configurable<bool> cfgPurityId{"cfgPurityId", false, "Purity correction"};
  ConfigurableAxis multTPCBins{"multTPCBins", {150, 0, 150}, "TPC Multiplicity bins"};
  ConfigurableAxis multFT0MBins{"multFT0MBins", {1000, 0, 5000}, "Forward Multiplicity bins"};
  ConfigurableAxis qNBins{"qNBins", {1000, 0., 100.}, "nth moments bins"};
  ConfigurableAxis nPairBins{"nPairBins", {2000, 0, 10000}, "nPair bins"};
  ConfigurableAxis dcaXYBins{"dcaXYBins", {100, -0.15, 0.15}, "dcaXY bins"};
  ConfigurableAxis dcaZBins{"dcaZBins", {500, -1.2, 1.2}, "dcaZ bins"};

  Configurable<std::vector<double>> ptBins{"ptBins", {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00}, "p_{T} bins"};
  Configurable<std::vector<double>> etaBins{"etaBins", {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8}, "#eta bins"};

  Configurable<std::string> cfgUrlCCDB{"cfgUrlCCDB", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cfgPathCCDB{"cfgPathCCDB", "Users/t/tgahlaut/weightCorr/", "Path for ccdb-object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};

  TH1D* hWeightPt = nullptr;
  TH1D* hPurePt = nullptr;
  TH1D* hWeightPtPi = nullptr;
  TH1D* hWeightPtKa = nullptr;
  TH1D* hWeightPtPr = nullptr;
  TH1D* hPurePtPi = nullptr;
  TH1D* hPurePtKa = nullptr;
  TH1D* hPurePtPr = nullptr;

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

  void init(InitContext const&)
  {
    if (cfgLoadEff) {
      // Set CCDB url
      ccdb->setURL(cfgUrlCCDB.value);
      ccdb->setCaching(true);

      TList* lst = ccdb->getForTimeStamp<TList>(cfgPathCCDB.value, -1);
      hWeightPt = reinterpret_cast<TH1D*>(lst->FindObject("hWeightPt"));
      hWeightPtPi = reinterpret_cast<TH1D*>(lst->FindObject("hWeightPtPi"));
      hWeightPtKa = reinterpret_cast<TH1D*>(lst->FindObject("hWeightPtKa"));
      hWeightPtPr = reinterpret_cast<TH1D*>(lst->FindObject("hWeightPtPr"));
      hPurePtPi = reinterpret_cast<TH1D*>(lst->FindObject("hPurePtPi"));
      hPurePtKa = reinterpret_cast<TH1D*>(lst->FindObject("hPurePtKa"));
      hPurePtPr = reinterpret_cast<TH1D*>(lst->FindObject("hPurePtPr"));

      if (!hWeightPt || !hWeightPtPi || !hWeightPtKa || !hWeightPtPr || !hPurePtPi || !hPurePtKa || !hPurePtPr) {
        LOGF(info, "FATAL!! Could not find required histograms in CCDB");
      }
    }

    const AxisSpec axisCol{3, 1, 4, ""};
    const AxisSpec axisTrack{5, 1, 6, ""};
    const AxisSpec axisEvents{10, 1, 11, "Counts"};
    const AxisSpec axisEta{etaBins, "#eta"};
    const AxisSpec axisPhi{nPhiBins, 0., +7., "#phi (rad)"};
    const AxisSpec axisY{100, -0.6, 0.6, "y"};
    const AxisSpec axisPt{ptBins, "p_{T} (GeV/c)"};
    const AxisSpec axisP{nPBins, 0., 3., "p (GeV/c)"};
    const AxisSpec axisInnerParam{nPBins, 0., 3., "p_{InnerParam } (GeV/c)"};
    const AxisSpec axisPart{nPartBins, 0., 18., " "};
    const AxisSpec axisQn{qNBins, ""};
    const AxisSpec axisNpair{nPairBins, "N_{pairs}"};
    const AxisSpec axisMeanPt{100, 0., 3., "M(p_{T}) (GeV/c)"};
    const AxisSpec axisMult{100, 0, 100, "N_{ch}"};
    const AxisSpec axisMultTPC{multTPCBins, "N_{TPC} "};
    const AxisSpec axisMultFT0M{multFT0MBins, "N_{FT0M}"};
    const AxisSpec axisCentFT0M{101, 0, 101, "FT0M (%)"};
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

    HistogramConfigSpec qNHist({HistType::kTHnSparseD, {axisCentFT0M, axisQn}});
    HistogramConfigSpec partHist({HistType::kTHnSparseD, {axisCentFT0M, axisPart}});
    HistogramConfigSpec qNMCHist({HistType::kTHnSparseD, {axisCentFT0M, axisQn}});
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

    hist.add("QA/Charged/h_Pt", "p_{T}", kTH1D, {axisPt});
    hist.add("QA/Charged/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/Charged/h_Phi", "#phi ", kTH1D, {axisPhi});
    hist.add("QA/Charged/h_DcaZ", "DCA_{Z}", kTH1D, {axisDCAz});
    hist.add("QA/Charged/h_DcaXY", "DCA_{XY}", kTH1D, {axisDCAxy});
    hist.add("QA/Charged/h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/Charged/h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/Charged/h_Pt_weighted", "weighted pT distribution", kTH1D, {axisPt});
    hist.add("QA/Charged/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/Charged/h2_Pt_centFT0M", "p_{T} in centrality Classes ", kTH2D, {{axisCentFT0M}, {axisPt}});
    hist.add("QA/Charged/h3_PtEtaPhi", "p_{T}, #eta, #phi ", kTHnSparseD, {{axisPt}, {axisEta}, {axisPhi}});

    hist.addClone("QA/Charged/", "QA/Pion/");

    hist.add("QA/Pion/before/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/before/h2_TPCNsigma_tof", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/before/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist);
    hist.add("QA/Pion/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist);

    hist.add("QA/Pion/h_Rap", "y ", kTH1D, {axisY});
    hist.add("QA/Pion/h_PtPos", "p_{T} (positive) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtNeg", "p_{T} (negative) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtTruth", "p_{T} (Truth)", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtPosTruth", "p_{T} (positive) (Truth)", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtNegTruth", "p_{T} (negative) (Truth) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist);
    hist.add("QA/Pion/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist);
    hist.add("QA/Pion/h2_TPCSignal", "TPC Signal ", tpcSignalHist);
    hist.add("QA/Pion/h2_TOFSignal", "TOF Signal", tofSignalHist);
    hist.add("QA/Pion/innerParam/h2_TPCSignal", "TPC Signal", tpcSignalHist1);
    hist.add("QA/Pion/h2_pvsm2", "p vs m^{2}", pvsM2Hist);
    hist.add("QA/Pion/h_PtTruth_primary", "p_{T} (Truth Primary)", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtTruth_secondary", "p_{T} (Truth Secondary)", kTH1D, {axisPt});
    hist.addClone("QA/Pion/", "QA/Kaon/");
    hist.addClone("QA/Pion/", "QA/Proton/");

    // AnalysisPlots
    hist.add("Analysis/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Analysis/Charged/h_N_CentFT0M", "Multiplicity vs CentFT0M", kTHnSparseD, {axisCentFT0M, axisMult});
    hist.add("Analysis/Charged/h_Npair_CentFT0M", "Npair vs CentFT0M", kTHnSparseD, {axisCentFT0M, axisNpair});
    hist.add("Analysis/Charged/h_Q1_CentFT0M", "Q1 vs CentFT0M", qNHist);
    hist.add("Analysis/Charged/h_Q2_CentFT0M", "Q2 vs CentFT0M", qNHist);
    hist.add("Analysis/Charged/h_twopart1_CentFT0M", "twopart (neum)", qNMCHist);
    hist.add("Analysis/Charged/h_mean_pT_Mult_var", " <p_{T}> vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_twopart_Mult_var", "Twopart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/p_twopart_CentFT0M", "Twopart vs cent_{FT0M} ", kTProfile, {axisCentFT0M});
    hist.add("Analysis/Charged/p_mean_pT_CentFT0M", " <p_{T}> vs cent_{FT0M} ", kTProfile, {axisCentFT0M});

    hist.addClone("Analysis/Charged/", "Analysis/Pion/");
    hist.addClone("Analysis/Charged/", "Analysis/Kaon/");
    hist.addClone("Analysis/Charged/", "Analysis/Proton/");

    // MC Generated
    hist.add("Gen/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("Gen/h_VtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/h_VtxZ_b", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/h_NSim", "Truth Multiplicity TPC", kTH1D, {axisMultTPC});
    hist.add("Gen/h2_NTPC_NSim", "Reco vs Truth Multiplicty TPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});

    hist.add("Gen/Charged/h_EtaTruth", "#eta ", kTH1D, {axisEta});
    hist.add("Gen/Charged/h_PhiTruth", "#phi ", kTH1D, {axisPhi});
    hist.add("Gen/Charged/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("Gen/Charged/h2_Pt_EtaTruth", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPt}});
    hist.add("Gen/Charged/h2_PtTruth_centFT0M", "p_{T} in centrality Classes ", kTH2D, {{axisCentFT0M}, {axisPt}});
    hist.add("Gen/Charged/h_PtEtaPhiTruth", "p_{T}, #eta, #phi ", kTHnSparseD, {{axisPt}, {axisEta}, {axisPhi}});

    hist.add("Gen/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Gen/Charged/h_N_CentFT0M", "Multiplicity vs CentFT0M", kTHnSparseD, {axisCentFT0M, axisMult});
    hist.add("Gen/Charged/h_Npair_CentFT0M", "Npair vs CentFT0M", kTHnSparseD, {axisCentFT0M, axisNpair});
    hist.add("Gen/Charged/h_Q1_CentFT0M", "Q1", qNMCHist);
    hist.add("Gen/Charged/h_Q2_CentFT0M", "Q2", qNMCHist);
    hist.add("Gen/Charged/h_twopart1_CentFT0M", "twopart (neum)", qNMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_var", " <p_{T}> vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_var", "Twopart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/p_twopart_CentFT0M", "Twopart vs CentFT0M ", kTProfile, {axisCentFT0M});
    hist.add("Gen/Charged/p_mean_pT_CentFT0M", " <p_{T}> vs CentFT0M ", kTProfile, {axisCentFT0M});

    hist.addClone("Gen/Charged/", "Gen/Pion/");

    hist.add("Gen/Pion/h_RapTruth", "y", kTH1D, {axisY});
    hist.add("Gen/Pion/h_PtPosTruth", "p_{T} (positive) ", kTH1D, {axisPt});
    hist.add("Gen/Pion/h_PtNegTruth", "p_{T} (negative) ", kTH1D, {axisPt});

    hist.addClone("Gen/Pion/", "Gen/Kaon/");
    hist.addClone("Gen/Pion/", "Gen/Proton/");

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
  }

  float centFT0M = 0.;
  int nTPC = 0, nFT0M = 0;

  // Event selection cuts:
  template <typename T>
  bool selRun3Col(T const& col)
  {
    hist.fill(HIST("QA/after/h_counts_evSelCuts"), kBeforeSelCol);

    if (cfgPosZ) {
      if (std::abs(col.posZ()) > cfgCutPosZ) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), kSelColPosZ);
    }

    centFT0M = col.centFT0M();
    nTPC = col.multNTracksHasTPC();
    nFT0M = col.multFT0M();

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

  // Cuts to reject the tracks
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

  template <PIDType s1, PIDType s2, PIDType s3, typename T>
  bool identifyParticle(T const& track, float momThreshold)
  {

    const int sp = static_cast<int>(s1);
    const int sq = static_cast<int>(s2);
    const int sr = static_cast<int>(s3);
    std::vector<float> vTpcNSigma = {-999., track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::vector<float> vTofNSigma = {-999., track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    bool isTofPidFlag = false, isTpcPidFlag = false;

    float nSigmaSelCut = 0.;
    float nSigmaTofRejCut = 0.;
    float nSigmaTpcRejCut = 0.;
    if constexpr (s1 == PIDType::kProtons) {
      nSigmaSelCut = cfgSelCutNSigPr;
      nSigmaTofRejCut = std::fabs(vTofNSigma[sp]);
      nSigmaTpcRejCut = std::fabs(vTpcNSigma[sp]);
    } else if constexpr (s1 == PIDType::kKaons) {
      nSigmaSelCut = cfgSelCutNSigKa;
      nSigmaTofRejCut = std::fabs(vTofNSigma[sp]);
      nSigmaTpcRejCut = std::fabs(vTpcNSigma[sp]);
    } else if constexpr (s1 == PIDType::kPions) {
      nSigmaSelCut = cfgSelCutNSigPi;
      nSigmaTofRejCut = cfgRejCutNSigPi;
      nSigmaTpcRejCut = cfgRejCutNSigPi;
    }

    if (track.hasTOF()) {
      if (std::fabs(vTofNSigma[sp]) < nSigmaSelCut &&
          std::fabs(vTofNSigma[sq]) > nSigmaTofRejCut &&
          std::fabs(vTofNSigma[sr]) > nSigmaTofRejCut) {
        isTofPidFlag = true;
      }
      if (std::fabs(vTpcNSigma[sp]) < cfgCutNSig2) {
        isTpcPidFlag = true;
      }
    } else { // select from TPC Only
      if (track.p() >= momThreshold) {
        return false;
      }
      if (std::fabs(vTpcNSigma[sp]) < nSigmaSelCut &&
          std::fabs(vTpcNSigma[sq]) > nSigmaTpcRejCut &&
          std::fabs(vTpcNSigma[sr]) > nSigmaTpcRejCut) {
        isTofPidFlag = true;
        isTpcPidFlag = true;
      }
    }

    if (isTofPidFlag && isTpcPidFlag) {
      return true; // Track is identified as one of the particles
    }

    return false; // Track is not identified as any of the particles
  }

  // Get corrected weight for the track:
  template <typename T1>
  float getCorrectedWeight(T1 hWeightPt, T1 hPurePt, float pt, bool cfgWeightPt, bool cfgPurity)
  {
    float weight = 1.0;
    float purity = 1.0;

    if (cfgPurity) {
      purity = hPurePt->GetBinContent(hPurePt->FindBin(pt));
    }

    if (cfgWeightPt) {
      float weightPt = hWeightPt->GetBinContent(hWeightPt->FindBin(pt));
      weight = purity * weightPt;
    }

    return weight;
  }

  // Fill hist before selection cuts:
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

  // Fill hist after selection cuts:
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
  }

  // Fill Charged particles QA:
  template <typename T>
  void fillChargedQAHistos(T const& track, float centFT0M)
  {
    hist.fill(HIST("QA/Charged/h_Eta"), track.eta());
    hist.fill(HIST("QA/Charged/h_Phi"), track.phi());
    hist.fill(HIST("QA/Charged/h2_Pt_centFT0M"), centFT0M, track.pt());
    hist.fill(HIST("QA/Charged/h3_PtEtaPhi"), track.pt(), track.eta(), track.phi());
    hist.fill(HIST("QA/Charged/h2_Pt_Eta"), track.eta(), track.pt());
    hist.fill(HIST("QA/Charged/h_DcaZ"), track.dcaZ());
    hist.fill(HIST("QA/Charged/h_DcaXY"), track.dcaXY());
    hist.fill(HIST("QA/Charged/h2_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST("QA/Charged/h2_DcaZ"), track.pt(), track.dcaZ());

    hist.fill(HIST("QA/after/h_TPCChi2perCluster"), track.tpcChi2NCl());
    hist.fill(HIST("QA/after/h_ITSChi2perCluster"), track.itsChi2NCl());
    hist.fill(HIST("QA/after/h_crossedTPC"), track.tpcNClsCrossedRows());
  }

  // Fill before PID cut QA hist:
  template <typename T>
  void fillBeforePIDQAHistos(T const& track)
  {
    hist.fill(HIST("QA/before/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/before/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/before/h2_pvsm2"), track.mass() * track.mass(), track.p());

    if (!track.hasTOF())
      hist.fill(HIST("QA/Pion/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPi());
    if (track.hasTOF())
      hist.fill(HIST("QA/Pion/before/h2_TPCNsigma_tof"), track.p(), track.tpcNSigmaPi());
    hist.fill(HIST("QA/Pion/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/before/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
    if (!track.hasTOF())
      hist.fill(HIST("QA/Proton/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPr());
    if (track.hasTOF())
      hist.fill(HIST("QA/Proton/before/h2_TPCNsigma_tof"), track.p(), track.tpcNSigmaPr());
    hist.fill(HIST("QA/Proton/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/before/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
    if (!track.hasTOF())
      hist.fill(HIST("QA/Kaon/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaKa());
    if (track.hasTOF())
      hist.fill(HIST("QA/Kaon/before/h2_TPCNsigma_tof"), track.p(), track.tpcNSigmaKa());
    hist.fill(HIST("QA/Kaon/before/h2_TOFNsigma"), track.p(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/before/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());
  }

  // Moments Calculation:
  void moments(float pt, float weight, double& Q1, double& Q2)
  {
    Q1 += pt * weight;
    Q2 += pt * pt * weight * weight;
  }

  template <int Mode, typename T, typename T1>
  void fillIdParticleQAHistos(T const& track, float rap, float nSigmaTPC, float nSigmaTOF, float centFT0M, T1 hWeightPt, T1 hPurePt, bool cfgWeightPtId, bool cfgPurityId, double& NW, double& NW2, double& Q1, double& Q2)
  {
    float pt = track.pt();
    float eta = track.eta();
    float phi = track.phi();
    float weight = getCorrectedWeight(hWeightPt, hPurePt, pt, cfgWeightPtId, cfgPurityId);

    if (std::abs(weight) < cfgMinWeight)
      return;

    NW += weight;
    NW2 += weight * weight;
    moments(pt, weight, Q1, Q2);

    if (track.sign() > 0) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtPos"), pt);
    }
    if (track.sign() < 0) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtNeg"), pt);
    }

    if (cfgWeightPtEtaId)
      hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Eta_weighted"), eta, pt, weight);

    hist.fill(HIST(Dire[Mode]) + HIST("h_Pt_weighted"), pt, weight);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_centFT0M"), centFT0M, pt);
    hist.fill(HIST(Dire[Mode]) + HIST("h3_PtEtaPhi"), pt, eta, phi);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Eta"), eta, pt);
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

  template <int Mode>
  void fillPtMCHist(bool cfgGen, float pt, float eta, float rap, float phi, float centFT0M, int pid, int pdgCodePos, int pdgCodeNeg)
  {
    if (cfgGen) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_EtaTruth"), eta);
      hist.fill(HIST(Dire[Mode]) + HIST("h_RapTruth"), rap);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtTruth_centFT0M"), centFT0M, pt);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_EtaTruth"), eta, pt);
      hist.fill(HIST(Dire[Mode]) + HIST("h_PhiTruth"), phi);
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtEtaPhiTruth"), pt, eta, phi);
    }

    if (pid == pdgCodePos) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtPosTruth"), pt);
    }
    if (pid == pdgCodeNeg) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtNegTruth"), pt);
    }
  }

  template <int Mode>
  void fillAnalysisHistos(bool cfgCorrection, float centFT0M, double nW, double nW2, double Q1, double Q2)
  {
    if (nW == 0) {
      return;
    }

    hist.fill(HIST(Dire[Mode]) + HIST("h_Mult"), nW);

    double nPair = 0., meanPt = 0., twopart1 = 0., twopart = 0.;

    if (cfgCorrection) {
      nPair = nW * nW - nW2;
    } else {
      if (nW > 1) {
        nPair = nW * (nW - 1);
      }
    }

    if (nPair > 0) {
      meanPt = Q1 / nW;
      twopart1 = (Q1 * Q1 - Q2);
      twopart = twopart1 / nPair;

      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Mult_var"), centFT0M, meanPt);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Mult_var"), centFT0M, twopart);
      hist.fill(HIST(Dire[Mode]) + HIST("p_mean_pT_CentFT0M"), centFT0M, meanPt);
      hist.fill(HIST(Dire[Mode]) + HIST("p_twopart_CentFT0M"), centFT0M, twopart);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart1_CentFT0M"), centFT0M, twopart1);
    }

    hist.fill(HIST(Dire[Mode]) + HIST("h_Mult"), nW);
    hist.fill(HIST(Dire[Mode]) + HIST("h_N_CentFT0M"), centFT0M, nW);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Npair_CentFT0M"), centFT0M, nPair);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q1_CentFT0M"), centFT0M, Q1);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q2_CentFT0M"), centFT0M, Q2);
  }

  template <bool DataFlag, bool RecoFlag, typename C, typename T>
  void fillRecoHistos(C const& col, T const& tracks)
  {
    double nChW = 0., nChW2 = 0., q1Ch = 0., q2Ch = 0.;
    float wghtCh = 1.0;
    double nPiW = 0., nPiW2 = 0., q1Pi = 0., q2Pi = 0.;
    double nKaW = 0., nKaW2 = 0., q1Ka = 0., q2Ka = 0.;
    double nPrW = 0., nPrW2 = 0., q1Pr = 0., q2Pr = 0.;
    float pt = 0., eta = 0., phi = 0.;

    hist.fill(HIST("QA/h_collisions_info"), kTotCol);

    if constexpr (DataFlag) {
      if (!selRun3Col(col)) {
        return;
      }
    }

    hist.fill(HIST("QA/h_collisions_info"), kPassSelCol);

    fillAfterQAHistos(col);

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
      if constexpr (RecoFlag) {
        hist.fill(HIST("Tracks/h2_tracks_pid_before_sel"), track.mcParticle().pdgCode(), track.pt());

        auto mc = track.template mcParticle_as<aod::McParticles>();
        pt = mc.pt();
        eta = mc.eta();
        phi = mc.phi();

        if (mc.isPhysicalPrimary()) {
          hist.fill(HIST("QA/after/h_DCAxy_primary"), track.dcaXY());
          hist.fill(HIST("QA/after/h_DCAz_primary"), track.dcaZ());
        } else {
          hist.fill(HIST("QA/after/h_DCAxy_secondary"), track.dcaXY());
          hist.fill(HIST("QA/after/h_DCAz_secondary"), track.dcaZ());
        }
      }

      // Charged particles:
      wghtCh = getCorrectedWeight(hWeightPt, hPurePt, pt, cfgWeightPtCh, false);

      if (std::abs(wghtCh) < cfgMinWeight)
        return;

      nChW += wghtCh;
      nChW2 += wghtCh * wghtCh;
      moments(pt, wghtCh, q1Ch, q2Ch);
      hist.fill(HIST("QA/Charged/h_Pt"), track.pt());
      fillChargedQAHistos(track, centFT0M);

      fillBeforePIDQAHistos(track);

      // identified particles:
      if (cfgRejTrk && rejectTracks(track)) {
        continue;
      }
      auto selIDPion = identifyParticle<PIDType::kPions, PIDType::kKaons, PIDType::kProtons>(track, cfgCutPiThrsldP);
      auto selIDKaon = identifyParticle<PIDType::kKaons, PIDType::kPions, PIDType::kProtons>(track, cfgCutKaThrsldP);
      auto selIDProton = identifyParticle<PIDType::kProtons, PIDType::kPions, PIDType::kKaons>(track, cfgCutPrThrsldP);
      if (selIDPion && pt >= cfgCutPiPtMin) {
        hist.fill(HIST("QA/Pion/h_Pt"), track.pt());
        fillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, centFT0M, hWeightPtPi, hPurePtPi, cfgWeightPtId, cfgPurityId, nPiW, nPiW2, q1Pi, q2Pi);
      }
      if (selIDKaon && pt >= cfgCutKaPtMin) {
        hist.fill(HIST("QA/Kaon/h_Pt"), track.pt());
        fillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, centFT0M, hWeightPtKa, hPurePtKa, cfgWeightPtId, cfgPurityId, nKaW, nKaW2, q1Ka, q2Ka);
      }
      if (selIDProton && pt >= cfgCutPrPtMin) {
        hist.fill(HIST("QA/Proton/h_Pt"), track.pt());
        fillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, centFT0M, hWeightPtPr, hPurePtPr, cfgWeightPtId, cfgPurityId, nPrW, nPrW2, q1Pr, q2Pr);
      }

      if constexpr (RecoFlag) {
        auto mc = track.template mcParticle_as<aod::McParticles>();
        int pid = mc.pdgCode();
        if (selIDPion && pt >= cfgCutPiPtMin) {
          if (std::abs(pid) == kPiPlus) {
            hist.fill(HIST("QA/Pion/h_PtTruth"), pt);
            fillPtMCHist<QA_Pion>(false, pt, eta, rapPi, phi, centFT0M, pid, kPiPlus, kPiMinus);
            if (mc.isPhysicalPrimary()) {
              hist.fill(HIST("QA/Pion/h_PtTruth_primary"), pt);
            } else {
              hist.fill(HIST("QA/Pion/h_PtTruth_secondary"), pt);
            }
          }
        }
        if (selIDKaon && pt >= cfgCutKaPtMin) {
          if (std::abs(pid) == kKPlus) {
            hist.fill(HIST("QA/Kaon/h_PtTruth"), pt);
            fillPtMCHist<QA_Kaon>(false, pt, eta, rapKa, phi, centFT0M, pid, kKPlus, kKMinus);
            if (mc.isPhysicalPrimary()) {
              hist.fill(HIST("QA/Kaon/h_PtTruth_primary"), pt);
            } else {
              hist.fill(HIST("QA/Kaon/h_PtTruth_secondary"), pt);
            }
          }
        }
        if (selIDProton && pt >= cfgCutPrPtMin) {
          if (std::abs(pid) == kProton) {
            hist.fill(HIST("QA/Proton/h_PtTruth"), pt);
            fillPtMCHist<QA_Proton>(false, pt, eta, rapPr, phi, centFT0M, pid, kProton, kProtonBar);
            if (mc.isPhysicalPrimary()) {
              hist.fill(HIST("QA/Proton/h_PtTruth_primary"), pt);
            } else {
              hist.fill(HIST("QA/Proton/h_PtTruth_secondary"), pt);
            }
          }
        }
      }
    }
    fillAnalysisHistos<Analysis_Charged>(cfgCorrection, centFT0M, nChW, nChW2, q1Ch, q2Ch);
    fillAnalysisHistos<Analysis_Pion>(cfgCorrection, centFT0M, nPiW, nPiW2, q1Pi, q2Pi);
    fillAnalysisHistos<Analysis_Kaon>(cfgCorrection, centFT0M, nKaW, nKaW2, q1Ka, q2Ka);
    fillAnalysisHistos<Analysis_Proton>(cfgCorrection, centFT0M, nPrW, nPrW2, q1Pr, q2Pr);
  }

  template <typename C, typename M>
  void fillGenHistos(C const& mcCol, M const& mcParticles)
  {
    int nSim = 0;
    double nChSim = 0., q1ChSim = 0., q2ChSim = 0.;
    double nPiSim = 0., q1PiSim = 0., q2PiSim = 0.;
    double nKaSim = 0., q1KaSim = 0., q2KaSim = 0.;
    double nPrSim = 0., q1PrSim = 0., q2PrSim = 0.;
    float pt = 0, eta = 0, phi = 0, rap = 0;
    for (auto const& mcPart : mcParticles) {
      if (!mcPart.isPhysicalPrimary()) {
        continue;
      }

      auto pid = mcPart.pdgCode();
      if (std::abs(pid) != kElectron && std::abs(pid) != kMuonMinus && std::abs(pid) != kPiPlus && std::abs(pid) != kKPlus && std::abs(pid) != kProton) {
        continue;
      }

      pt = mcPart.pt();
      eta = mcPart.eta();
      phi = mcPart.phi();
      if (std::abs(eta) < cfgCutEta) {
        nSim++;
      }

      if (pt >= cfgCutPtMin && pt < cfgCutPtMax && std::abs(eta) < cfgCutEta) {
        nChSim++;
        moments(pt, 1.0, q1ChSim, q2ChSim);
        hist.fill(HIST("Gen/Charged/h_PtTruth"), pt);
        hist.fill(HIST("Gen/Charged/h_EtaTruth"), eta);
        hist.fill(HIST("Gen/Charged/h_PhiTruth"), phi);
        hist.fill(HIST("Gen/Charged/h2_Pt_EtaTruth"), eta, pt);
        hist.fill(HIST("Gen/Charged/h2_PtTruth_centFT0M"), centFT0M, pt);
        hist.fill(HIST("Gen/Charged/h_PtEtaPhiTruth"), pt, eta, phi);

        rap = mcPart.y();
        if (std::abs(pid) == kPiPlus && mcPart.pt() > cfgCutPiPtMin) {
          nPiSim++;
          moments(pt, 1.0, q1PiSim, q2PiSim);
          hist.fill(HIST("Gen/Pion/h_PtTruth"), pt);
          fillPtMCHist<Gen_Pion>(true, pt, eta, rap, phi, centFT0M, pid, kPiMinus, kPiMinus);
        }
        if (std::abs(pid) == kKPlus && mcPart.pt() > cfgCutKaPtMin) {
          nKaSim++;
          moments(pt, 1.0, q1KaSim, q2KaSim);
          hist.fill(HIST("Gen/Kaon/h_PtTruth"), pt);
          fillPtMCHist<Gen_Kaon>(true, pt, eta, rap, phi, centFT0M, pid, kKMinus, kKMinus);
        }
        if (std::abs(pid) == kProton && mcPart.pt() > cfgCutPrPtMin) {
          nPrSim++;
          moments(pt, 1.0, q1PrSim, q2PrSim);
          hist.fill(HIST("Gen/Proton/h_PtTruth"), pt);
          fillPtMCHist<Gen_Proton>(true, pt, eta, rap, phi, centFT0M, pid, kProtonBar, kProtonBar);
        }
      }
    }
    fillAnalysisHistos<Gen_Charged>(false, centFT0M, nChSim, nChSim, q1ChSim, q2ChSim);
    fillAnalysisHistos<Gen_Pion>(false, centFT0M, nPiSim, nPiSim, q1PiSim, q2PiSim);
    fillAnalysisHistos<Gen_Kaon>(false, centFT0M, nKaSim, nKaSim, q1KaSim, q2KaSim);
    fillAnalysisHistos<Gen_Proton>(false, centFT0M, nPrSim, nPrSim, q1PrSim, q2PrSim);
    hist.fill(HIST("Gen/h_Counts"), 2);
    hist.fill(HIST("Gen/h_VtxZ"), mcCol.posZ());
    hist.fill(HIST("Gen/h_NSim"), nSim);
    hist.fill(HIST("Gen/h2_NTPC_NSim"), nSim, nTPC);
  }

  template <typename M, typename C, typename T, typename P>
  void analyzeMC(M const& mcCol, C const& cols, T const& tracks, P const& mcParts)
  {
    // fillBeforeQAHistos(cols.begin(), tracks);
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
    fillRecoHistos<false, true>(cols.begin(), sTracks);
    fillGenHistos(mcCol, mcParts);
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
