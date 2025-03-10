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

#include <utility>
#include <vector>
#include <string>
#include <TPDGCode.h>

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
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace std;

struct MeanPtFlucId {
  Configurable<int> nPBins{"nPBins", 300, ""};
  Configurable<int> nPartBins{"nPartBins", 250, ""};
  Configurable<int> nCentBins{"nCentBins", 101, ""};
  Configurable<int> nRapBins{"nRapBins", 100, ""};
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
  Configurable<bool> cfgLoadEff{"cfgLoadEff", true, "Load efficiency"};
  Configurable<bool> cfgEffPtCh{"cfgEffPtCh", true, "Efficiency correction (pT) for charged particles"};
  Configurable<bool> cfgEffPtPi{"cfgEffPtPi", true, "Efficiency correction (pT) for pions"};
  Configurable<bool> cfgEffPtKa{"cfgEffPtKa", true, "Efficiency correction (pT) for kaons"};
  Configurable<bool> cfgEffPtPr{"cfgEffPtPr", true, "Efficiency correction (pT) for protons"};
  Configurable<bool> cfgEffPtYPi{"cfgEffPtYPi", true, "Efficiency correction (pT, y) for pions"};
  Configurable<bool> cfgEffPtYKa{"cfgEffPtYKa", true, "Efficiency correction (pT, y) for kaons"};
  Configurable<bool> cfgEffPtYPr{"cfgEffPtYPr", true, "Efficiency correction (pT, y) for protons"};
  Configurable<bool> cfgPidCut{"cfgPidCut", false, ""};
  Configurable<bool> cfgPDGCodeOnly{"cfgPDGCodeOnly", true, ""};
  Configurable<bool> cfgMCReco{"cfgMCReco", false, ""};
  Configurable<bool> cfgMCTruth{"cfgMCTruth", false, ""};
  Configurable<bool> cfgPosZ{"cfgPosZ", true, "Position Z"};
  Configurable<bool> cfgSel8{"cfgSel8", true, "Sel8 trigger"};
  Configurable<bool> cfgNoSameBunchPileup{"cfgNoSameBunchPileup", true, "kNoSameBunchPileup"};
  Configurable<bool> cfgIsVertexITSTPC{"cfgIsVertexITSTPC", true, "kIsVertexITSTPC"};
  Configurable<bool> cfgRejTrk{"cfgRejTrk", true, "Rejected Tracks"};
  Configurable<bool> cfgInvMass{"cfgInvMass", true, "electron Inv Mass cut selection"};
  Configurable<bool> cfgSelOR{"cfgSelOR", true, "Low OR High momentum "};
  Configurable<bool> cfgSelAND{"cfgSelAND", false, "Low AND High momentum"};
  Configurable<bool> cfgSelLow{"cfgSelLow", true, "PID selection cut for Low momentum"};
  Configurable<bool> cfgSelHigh{"cfgSelHigh", true, "PID selection cut for High momentum"};
  ConfigurableAxis multTPCBins{"multTPCBins", {150, 0, 150}, "TPC Multiplicity bins"};
  ConfigurableAxis multFT0MBins{"multFT0MBins", {1000, 0, 5000}, "Forward Multiplicity bins"};
  ConfigurableAxis multFT0MMCBins{"multFT0MMCBins", {250, 0, 250}, "Forward Multiplicity bins"};
  ConfigurableAxis dcaXYBins{"dcaXYBins", {100, -0.15, 0.15}, "dcaXY bins"};
  ConfigurableAxis dcaZBins{"dcaZBins", {100, -1.2, 1.2}, "dcaZ bins"};
  ConfigurableAxis qNBins{"qNBins", {1000, 0., 100.}, "nth moments bins"};
  ConfigurableAxis tpNBins{"tpNBins", {300, 0., 3000.}, ""};
  ConfigurableAxis tpDBins{"tpDBins", {100, 0., 2000.}, ""};
  Configurable<std::vector<double>> ptBins{"ptBins", {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00}, "p_{T} bins"};
  Configurable<std::vector<double>> etaBins{"etaBins", {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8}, "#eta bins"};
  Configurable<std::vector<double>> rapBins{"rapBins", {-0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6}, "#rap bins"};

  Configurable<std::string> cfgUrlCCDB{"cfgUrlCCDB", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cfgPathCCDB{"cfgPathCCDB", "Users/t/tgahlaut/", "Path for ccdb-object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  std::vector<std::string> vEff = {"heffPt", "heffPtPi", "heffPtKa", "heffPtPr", "heffPtRapPi", "heffPtRapKa", "heffPtRapPr"};

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                                aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                                aod::pidTOFbeta, aod::pidTOFmass>;
  using MyRun3Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Cs, aod::CentFT0Ms>;
  using MyRun3MCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Cs, aod::CentFT0Ms, aod::McCollisionLabels>;
  using MyMCTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                               aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                               aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                               aod::pidTOFbeta, aod::pidTOFmass, aod::McTrackLabels>;

  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};

  TH1D* heffPt = nullptr;
  TH1D* heffPtPi = nullptr;
  TH1D* heffPtKa = nullptr;
  TH1D* heffPtPr = nullptr;
  TH2D* heffPtRapPi = nullptr;
  TH2D* heffPtRapKa = nullptr;
  TH2D* heffPtRapPr = nullptr;

  void init(InitContext const&)
  {
    if (cfgLoadEff) {
      // Set CCDB url
      ccdb->setURL(cfgUrlCCDB.value);
      ccdb->setCaching(true);

      TList* lst = ccdb->getForTimeStamp<TList>(cfgPathCCDB.value, -1);
      heffPt = reinterpret_cast<TH1D*>(lst->FindObject("heffPt"));
      heffPtPi = reinterpret_cast<TH1D*>(lst->FindObject("heffPtPi"));
      heffPtKa = reinterpret_cast<TH1D*>(lst->FindObject("heffPtKa"));
      heffPtPr = reinterpret_cast<TH1D*>(lst->FindObject("heffPtPr"));
      heffPtRapPi = reinterpret_cast<TH2D*>(lst->FindObject("heffPtRapPi"));
      heffPtRapKa = reinterpret_cast<TH2D*>(lst->FindObject("heffPtRapKa"));
      heffPtRapPr = reinterpret_cast<TH2D*>(lst->FindObject("heffPtRapPr"));

      if (!heffPt || !heffPtPi || !heffPtKa || !heffPtPr || !heffPtRapPi || !heffPtRapKa || !heffPtRapPr) {
        LOGF(info, "FATAL!! Could not find required histograms in CCDB");
      }
    }

    const AxisSpec axisEvents{10, 0, 10, "Counts"};
    const AxisSpec axisEta{etaBins, "#eta"};
    const AxisSpec axisPhi{nPhiBins, 0., +7., "#phi (rad)"};
    const AxisSpec axisY{rapBins, "y"};
    const AxisSpec axisPt{ptBins, "p_{T} (GeV/c)"};
    const AxisSpec axisP{nPBins, 0., 3., "p (GeV/c)"};
    const AxisSpec axisInnerParam{nPBins, 0., 3., "p_{InnerParam } (GeV/c)"};
    const AxisSpec axisPart{nPartBins, 0., 18., " "};
    const AxisSpec axisQn{qNBins, ""};
    const AxisSpec axisTpN{tpNBins, "(Q_{1}^{2} - Q_{2})"};
    const AxisSpec axisTpD{tpDBins, "N_{pairs}"};
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

    HistogramConfigSpec qNHist({HistType::kTHnSparseD, {axisMultTPC, axisQn, axisMultFT0M}});
    HistogramConfigSpec partHist({HistType::kTHnSparseD, {axisMultTPC, axisPart, axisMultFT0M}});
    HistogramConfigSpec denoHist({HistType::kTHnSparseD, {axisMultTPC, axisDeno, axisMultFT0M}});
    HistogramConfigSpec qNMCHist({HistType::kTHnSparseD, {axisMultTPC, axisQn, axisMultFT0M}});
    HistogramConfigSpec partMCHist({HistType::kTHnSparseD, {axisMultTPC, axisPart, axisMultFT0M}});
    HistogramConfigSpec denoMCHist({HistType::kTHnSparseD, {axisMultTPC, axisDeno, axisMultFT0M}});
    HistogramConfigSpec tofNSigmaHist({HistType::kTH2D, {axisP, axisTOFNsigma}});
    HistogramConfigSpec tofSignalHist({HistType::kTH2D, {axisP, axisTOFSignal}});
    HistogramConfigSpec tpcNSigmaHist({HistType::kTH2D, {axisP, axisTPCNsigma}});
    HistogramConfigSpec tpcSignalHist({HistType::kTH2D, {axisP, axisTPCSignal}});
    HistogramConfigSpec tpcTofHist({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});
    HistogramConfigSpec pvsM2Hist({HistType::kTH2D, {axisM2, axisP}});

    HistogramConfigSpec tofNSigmaHist1({HistType::kTH2D, {axisInnerParam, axisTOFNsigma}});
    HistogramConfigSpec tofSignalHist1({HistType::kTH2D, {axisInnerParam, axisTOFSignal}});
    HistogramConfigSpec tpcNSigmaHist1({HistType::kTH2D, {axisInnerParam, axisTPCNsigma}});
    HistogramConfigSpec tpcSignalHist1({HistType::kTH2D, {axisInnerParam, axisTPCSignal}});
    HistogramConfigSpec tpcTofHist1({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});
    HistogramConfigSpec pvsM2Hist1({HistType::kTH2D, {axisM2, axisInnerParam}});

    // QA Plots:
    hist.add("QA/before/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("QA/before/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});
    hist.add("QA/before/h_Pt", "p_{T}", kTH1D, {axisPt});
    hist.add("QA/before/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/before/h_Phi", "#phi ", kTH1D, {axisPhi});
    hist.add("QA/before/h_DcaZ", "DCA_{Z}", kTH1D, {axisDCAz});
    hist.add("QA/before/h_DcaXY", "DCA_{XY}", kTH1D, {axisDCAxy});
    hist.add("QA/before/h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/before/h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/before/h_NTPC", "N_{TPC}", kTH1D, {axisMultTPC});
    hist.add("QA/before/h_NFT0M", "FT0M Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/before/h_Cent", "FT0C (%)", kTH1D, {axisCentFT0C});
    hist.add("QA/before/h_CentM", "FT0M (%)", kTH1D, {axisCentFT0C});

    hist.add("QA/before/h2_TPCSignal", "TPC Signal", tpcSignalHist);
    hist.add("QA/before/h2_TOFSignal", "TOF Signal", tofSignalHist);
    hist.add("QA/before/h2_pvsm2", "p vs m^{2}", pvsM2Hist);

    hist.add("QA/before/innerParam/h2_TPCSignal", "TPC Signal", tpcSignalHist1);
    hist.add("QA/before/innerParam/h2_TOFSignal", "TOF Signal", tofSignalHist1);
    hist.add("QA/before/innerParam/h2_pvsm2", "p vs m^{2}", pvsM2Hist1);

    hist.addClone("QA/before/", "QA/after/");

    hist.add("QA/after/h_TPCChi2perCluster", "TPC #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/after/h_ITSChi2perCluster", "ITS #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/after/h_crossedTPC", "Crossed TPC", kTH1D, {axisCrossedTPC});
    hist.add("QA/after/h_NFT0C", "FT0C Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/after/h_invMass_gamma", "Inv Mass of #gamma", kTH1D, {axisMass});
    hist.add("QA/after/h_counts_evSelCuts", "Event selection cuts", kTH1D, {axisEvents});
    hist.add("QA/after/h_VtxZReco", "Simulated Vertex Z", kTH1D, {axisVtxZ});

    hist.add("QA/after/h2_PvsPinner", "p_{InnerParam} vs p", kTH2D, {{axisP}, {axisInnerParam}});
    hist.add("QA/after/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/after/h2_NTPC_Cent", "N_{TPC} vs FT0C(%)", kTH2D, {{axisCentFT0C}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NFT0M", "N_{TPC} vs N_{FT0M}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NFT0C", "N_{TPC} vs N_{FT0C}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NCh", "N_{ch} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NPi", "N_{Pi} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NKa", "N_{Ka} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("QA/after/h2_NTPC_NPr", "N_{Pr} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMultTPC}});

    hist.add("QA/after/p_NTPC_NFT0M", "N_{TPC} vs N_{FT0M} (Profile)", kTProfile, {axisMultFT0M});
    hist.add("QA/after/p_NTPC_NFT0C", "N_{TPC} vs N_{FT0C} (Profile)", kTProfile, {axisMultFT0M});
    hist.add("QA/after/p_NTPC_Cent", "N_{TPC} vs FT0C(%) (Profile)", kTProfile, {axisCentFT0C});

    hist.add("QA/after/h_Pt_weighted", "weighted pT distribution", kTH1D, {axisPt});
    hist.add("QA/after/h2_Pt_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/after/h2_pt_nch", "Truth", kTH2D, {{axisMult}, {axisPt}});
    hist.add("QA/after/h3_nft0m_pt_nch", "Reco", kTHnSparseD, {{axisMult}, {axisPt}, {axisMultFT0M}});
    hist.add("QA/after/h2_pt_nch_prof", "Truth", kTProfile, {axisMult});

    hist.add("QA/Pion/before/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/before/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist);
    hist.add("QA/Pion/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist);

    hist.add("QA/Pion/h_Rap", "y ", kTH1D, {axisY});
    hist.add("QA/Pion/h_RapTruth", "y ", kTH1D, {axisY});
    hist.add("QA/Pion/h_Eta", "Pseudorapidity ", kTH1D, {axisEta});
    hist.add("QA/Pion/h_EtaTruth", "Pseudorapidity (Reco Truth) ", kTH1D, {axisEta});
    hist.add("QA/Pion/h_Phi", "Azimuthal Distribution ", kTH1D, {axisPhi});
    hist.add("QA/Pion/h_DcaZ", "DCA_{z}", kTH1D, {axisDCAz});
    hist.add("QA/Pion/h_DcaXY", "DCA_{xy}", kTH1D, {axisDCAxy});
    hist.add("QA/Pion/h_Pt", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtPos", "p_{T} (positive) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtNeg", "p_{T} (negative) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtPosTruth", "p_{T} (positive) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_PtNegTruth", "p_{T} (negative) ", kTH1D, {axisPt});
    hist.add("QA/Pion/h_Pt_weighted", "weighted pT distribution", kTH1D, {axisPt});

    hist.add("QA/Pion/h2_Pt_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtPos_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtNeg_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtPosTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_PtNegTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_Pt_Eta", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/Pion/h2_DcaZ", "DCA_{z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/Pion/h2_DcaXY", "DCA_{xy}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/Pion/h2_Pt_Rap_weighted", "p_{T} vs #eta weighted", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h2_Pt_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtPos_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtNeg_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtPosTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_PtNegTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_pt_nch", "Reco", kTH2D, {{axisMult}, {axisPt}});
    hist.add("QA/Pion/h3_nft0m_pt_nch", "Reco", kTHnSparseD, {{axisMult}, {axisPt}, {axisMultFT0M}});
    hist.add("QA/Pion/h2_pt_nch_prof", "Reco", kTProfile, {axisMult});

    hist.add("QA/Pion/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist);
    hist.add("QA/Pion/h2_TPCNsigma_El", "n #sigma_{TPC, El}", tpcNSigmaHist);
    hist.add("QA/Pion/h2_TOFNsigma_El", "n #sigma_{TOF, El}", tofNSigmaHist);
    hist.add("QA/Pion/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist);
    hist.add("QA/Pion/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist);
    hist.add("QA/Pion/h2_TPCSignal", "TPC Signal ", tpcSignalHist);
    hist.add("QA/Pion/h2_TOFSignal", "TOF Signal", tofSignalHist);
    hist.add("QA/Pion/h2_pvsm2", "p vs m^{2}", pvsM2Hist);

    hist.add("QA/Pion/innerParam/before/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist1);
    hist.add("QA/Pion/innerParam/before/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist1);
    hist.add("QA/Pion/innerParam/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist1);
    hist.add("QA/Pion/innerParam/h2_TPCNsigma", "n #sigma_{TPC}", tpcNSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TPCNsigma_El", "n #sigma_{TPC, El}", tpcNSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TOFNsigma_El", "n #sigma_{TOF, El}", tofNSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TOFNsigma", "n #sigma_{TOF}", tofNSigmaHist1);
    hist.add("QA/Pion/innerParam/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", tpcTofHist1);
    hist.add("QA/Pion/innerParam/h2_TPCSignal", "TPC Signal ", tpcSignalHist1);
    hist.add("QA/Pion/innerParam/h2_TOFSignal", "TOF Signal", tofSignalHist1);
    hist.add("QA/Pion/innerParam/h2_pvsm2", "p vs m^{2}", pvsM2Hist1);

    hist.addClone("QA/Pion/", "QA/Kaon/");
    hist.addClone("QA/Pion/", "QA/Proton/");

    // Analysis Plots:
    hist.add("Analysis/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Analysis/Charged/h_Mult_weighted", "Multiplicity", kTH1D, {axisMult});
    hist.add("Analysis/Charged/h_Q1", "Q1", qNHist);
    hist.add("Analysis/Charged/h_Q2", "Q2", qNHist);
    hist.add("Analysis/Charged/h_Q3", "Q3", qNHist);
    hist.add("Analysis/Charged/h_Q4", "Q4", qNHist);
    hist.add("Analysis/Charged/h_mean_pT", " <p_{T}> ", kTH1D, {axisMeanPt});
    hist.add("Analysis/Charged/p_mean_pT_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/p_CheckNCh", " 1/denominator vs N_{TPC} ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/h_CheckNCh", " 1/denominator vs N_{TPC} ", denoHist);
    hist.add("Analysis/Charged/h_Q1_var", "Q1 vs N_{TPC}", qNHist);
    hist.add("Analysis/Charged/h_N_var", "N vs N_{TPC}", kTHnSparseD, {axisMultTPC, axisMult, axisMultFT0M});
    hist.add("Analysis/Charged/h_twopart_nume_Mult_var", "twopart numerator", kTHnSparseD, {axisMultTPC, axisTpN, axisMultFT0M});
    hist.add("Analysis/Charged/h_twopart_deno_Mult_var", "twopart denominator", kTHnSparseD, {axisMultTPC, axisTpD, axisMultFT0M});
    hist.add("Analysis/Charged/h_mean_pT_Mult_var", " <p_{T}> vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_mean_pT_Mult_skew", " <p_{T}> vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_mean_pT_Mult_kurto", " <p_{T}> vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_twopart_Mult_var", "Twopart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_twopart_Mult_skew", "Twopart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_twopart_Mult_kurto", "Twopart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_threepart_Mult_skew", "Threepart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_threepart_Mult_kurto", "Threepart vs N_{TPC} ", partHist);
    hist.add("Analysis/Charged/h_fourpart_Mult_kurto", "Fourpart vs N_{TPC} ", partHist);

    hist.addClone("Analysis/Charged/", "Analysis/Pion/");
    hist.addClone("Analysis/Charged/", "Analysis/Kaon/");
    hist.addClone("Analysis/Charged/", "Analysis/Proton/");

    // MC Generated
    hist.add("Gen/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("Gen/h_VtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/h_VtxZ_b", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/h_NTPC", "Mid rapidity Multiplicity", kTH1D, {axisMultTPC});
    hist.add("Gen/h_NFT0C", "Forward Multiplicity", kTH1D, {axisMultFT0MMC});
    hist.add("Gen/h2_NTPC_NFT0C", "N_{TPC} vs N_{FT0C}", kTH2D, {{axisMultFT0MMC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NFT0M", "N_{TPC} vs N_{FT0M} Reco", kTH2D, {{axisMultFT0M}, {axisMultTPC}});

    hist.add("Gen/h_NSim", "Truth Multiplicity TPC", kTH1D, {axisMultTPC});
    hist.add("Gen/h2_NTPC_NSim", "Reco vs Truth Multiplicty TPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NChSim_NSim", "Truth Multiplicty NCh vs NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NChSim", "Truth Multiplicty NCh vs Reco NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NPiSim", "Truth Multiplicty NPi vs Reco NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NKaSim", "Truth Multiplicty NKa vs Reco NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NTPC_NPrSim", "Truth Multiplicty NPr vs Reco NTPC", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/h2_NFT0C_NFT0CSim", "Reco vs Truth Multplicity FT0C", kTH2D, {{axisMultFT0MMC}, {axisMultFT0M}});

    hist.add("Gen/Charged/h2_Nid_NidSim", "reco vs truth multiplicity", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/Charged/h_EtaTruth", "#eta ", kTH1D, {axisEta});
    hist.add("Gen/Charged/h_PhiTruth", "#phi ", kTH1D, {axisPhi});
    hist.add("Gen/Charged/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("Gen/Charged/h2_PtTruth_Eta", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPt}});
    hist.add("Gen/Charged/h2_PtTruth_NFT0M", "p_{T} in Multiplicity Classes", kTH2D, {{axisPt}, {axisMultFT0M}});

    hist.add("Gen/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Gen/Charged/h_Mult_weighted", "Multiplicity", kTH1D, {axisMult});
    hist.add("Gen/Charged/h2_pt_nch", "Truth", kTH2D, {{axisMult}, {axisPt}});
    hist.add("Gen/Charged/h3_nft0m_pt_nch", "Truth", kTHnSparseD, {{axisMult}, {axisPt}, {axisMultFT0M}});
    hist.add("Gen/Charged/h2_pt_nch_prof", "Truth", kTProfile, {axisMult});
    hist.add("Gen/Charged/h_mean_pT", " <p_{T}> ", kTH1D, {axisMeanPt});

    hist.add("Gen/Charged/h_Q1", "Q1", qNMCHist);
    hist.add("Gen/Charged/h_Q2", "Q2", qNMCHist);
    hist.add("Gen/Charged/h_Q3", "Q3", qNMCHist);
    hist.add("Gen/Charged/h_Q4", "Q4", qNMCHist);
    hist.add("Gen/Charged/h_Q1_var", "Q1 vs N_{TPC}", qNMCHist);
    hist.add("Gen/Charged/h_N_var", "N vs N_{TPC}", kTHnSparseD, {axisMultTPC, axisMult, axisMultFT0M});
    hist.add("Gen/Charged/h_twopart_nume_Mult_var", "twopart numerator", kTHnSparseD, {axisMultTPC, axisTpN, axisMultFT0M});
    hist.add("Gen/Charged/h_twopart_deno_Mult_var", "twopart denominator", kTHnSparseD, {axisMultTPC, axisTpD, axisMultFT0M});

    hist.add("Gen/Charged/p_mean_pT_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
    hist.add("Gen/Charged/p_CheckNCh", " 1/denominator vs N_{TPC} ", kTProfile, {axisMultTPC});
    hist.add("Gen/Charged/h_CheckNCh", " 1/denominator vs N_{TPC} ", denoMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_var", " <p_{T}> vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_skew", " <p_{T}> vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_mean_pT_Mult_kurto", " <p_{T}> vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_var", "Twopart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_skew", "Twopart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_twopart_Mult_kurto", "Twopart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_threepart_Mult_skew", "Threepart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_threepart_Mult_kurto", "Threepart vs N_{TPC} ", partMCHist);
    hist.add("Gen/Charged/h_fourpart_Mult_kurto", "Fourpart vs N_{TPC} ", partMCHist);

    hist.addClone("Gen/Charged/", "Gen/Pion/");

    hist.add("Gen/Pion/h_RapTruth", "y", kTH1D, {axisY});
    hist.add("Gen/Pion/h2_PtTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("Gen/Pion/h2_PtPosTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("Gen/Pion/h2_PtNegTruth_Rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("Gen/Pion/h_PtPosTruth", "p_{T} (Positive)", kTH1D, {axisPt});
    hist.add("Gen/Pion/h_PtNegTruth", "p_{T} (negative)", kTH1D, {axisPt});
    hist.add("Gen/Pion/h2_PtPosTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});
    hist.add("Gen/Pion/h2_PtNegTruth_NFT0M", "p_{T} in Multiplicity Classes ", kTH2D, {{axisPt}, {axisMultFT0M}});

    hist.addClone("Gen/Pion/", "Gen/Kaon/");
    hist.addClone("Gen/Pion/", "Gen/Proton/");
  }

  enum Mode {
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

  static constexpr std::string_view Dire[] = {
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
    hist.fill(HIST("QA/after/h_counts_evSelCuts"), 0);

    if (cfgPosZ) {
      if (std::abs(col.posZ()) > cfgCutPosZ) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), 1);
    }

    if (cfgSel8) {
      if (!col.sel8()) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), 2);
    }
    if (cfgNoSameBunchPileup) {
      if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), 4);
    }

    if (cfgIsVertexITSTPC) {
      if (!col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        return false;
      }
      hist.fill(HIST("QA/after/h_counts_evSelCuts"), 5);
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

    if (track.pt() > cfgCutPtMax)
      return false;

    if (track.sign() == 0)
      return false;

    if (std::fabs(track.dcaZ()) > cfgCutDcaZ)
      return false;

    return true;
  }

  // Cuts to reject the tracks
  template <typename T>
  bool rejectTracks(T const& track)
  {
    if (((track.tpcNSigmaEl()) > -3. &&
         (track.tpcNSigmaEl()) < 5.) &&
        (std::fabs(track.tpcNSigmaPi()) > 3 &&
         std::fabs(track.tpcNSigmaKa()) > 3 &&
         std::fabs(track.tpcNSigmaPr()) > 3)) {
      return true;
    }

    return false;
  }

  template <typename T>
  bool selElectrons(T const& track)
  {
    if (std::fabs(track.tpcNSigmaEl()) < cfgCutNSig3) {
      return true;
    }

    return false;
  }

  // PID selction cuts for Low momentum Pions
  template <typename T>
  bool selLowPi(T const& track)
  {
    if (track.pt() >= cfgCutPiPtMin &&
        track.p() <= cfgCutPiThrsldP &&
        std::abs(track.rapidity(MassPiPlus)) < cfgCutRap) {
      if (!track.hasTOF() &&
          std::fabs(track.tpcNSigmaPi()) < cfgCutNSig2) {
        return true;
      }

      if (track.hasTOF() &&
          std::fabs(track.tpcNSigmaPi()) < cfgCutNSig2 &&
          std::fabs(track.tofNSigmaPi()) < cfgCutNSig3) {
        return true;
      }
    }
    return false;
  }

  // PID selction cuts for Low momentum Kaons
  template <typename T>
  bool selLowKa(T const& track)
  {
    if (track.pt() >= cfgCutKaPtMin &&
        track.p() <= cfgCutKaThrsldP &&
        std::abs(track.rapidity(MassKPlus)) < cfgCutRap) {
      if (!track.hasTOF() &&
          std::fabs(track.tpcNSigmaKa()) < cfgCutNSig2) {
        return true;
      }

      if (track.hasTOF() &&
          std::fabs(track.tpcNSigmaKa()) < cfgCutNSig2 &&
          std::fabs(track.tofNSigmaKa()) < cfgCutNSig3) {
        return true;
      }
    }

    return false;
  }

  // PID selction cuts for Low momentum Protons
  template <typename T>
  bool selLowPr(T const& track)
  {
    if (track.pt() >= cfgCutPrPtMin &&
        track.p() <= cfgCutPrThrsldP &&
        std::abs(track.rapidity(MassProton)) < cfgCutRap) {
      if (!track.hasTOF() &&
          std::fabs(track.tpcNSigmaPr()) < cfgCutNSig2) {
        return true;
      }

      if (track.hasTOF() &&
          std::fabs(track.tpcNSigmaPr()) < cfgCutNSig2 &&
          std::fabs(track.tofNSigmaPr()) < cfgCutNSig3) {
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
        std::fabs(track.tpcNSigmaPi()) < cfgCutNSig3 &&
        std::fabs(track.tofNSigmaPi()) < cfgCutNSig3) {

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
        std::fabs(track.tpcNSigmaKa()) < cfgCutNSig3 &&
        ((std::fabs(track.tofNSigmaKa()) < cfgCutNSig3 && track.p() <= cfgCutKaP3) ||
         (std::fabs(track.tofNSigmaKa()) < cfgCutNSig2 && track.p() > cfgCutKaP3))) {

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
        std::fabs(track.tpcNSigmaPr()) < cfgCutNSig3 &&
        std::fabs(track.tofNSigmaPr()) < cfgCutNSig3) {

      if (std::abs(track.rapidity(MassProton)) < cfgCutRap) {
        return true;
      }
    }

    return false;
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
    hist.fill(HIST("QA/before/h_Cent"), col.centFT0C());
    hist.fill(HIST("QA/after/h_CentM"), col.centFT0M());
    hist.fill(HIST("QA/before/h_NFT0M"), col.multFT0M());
  }

  // Fill hist after selection cuts:
  template <typename T>
  void fillAfterQAHistos(T const& col)
  {
    hist.fill(HIST("QA/after/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/after/h_Counts"), 2);
    hist.fill(HIST("QA/after/h_NTPC"), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/h_Cent"), col.centFT0C());
    hist.fill(HIST("QA/after/h_CentM"), col.centFT0M());
    hist.fill(HIST("QA/after/h_NFT0M"), col.multFT0M());
    hist.fill(HIST("QA/after/h_NFT0C"), col.multFT0C());
    hist.fill(HIST("QA/after/h2_NTPC_NFT0M"), col.multFT0M(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/h2_NTPC_NFT0C"), col.multFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/h2_NTPC_Cent"), col.centFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/p_NTPC_Cent"), col.centFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/p_NTPC_NFT0M"), col.multFT0M(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/p_NTPC_NFT0C"), col.multFT0C(), col.multNTracksHasTPC());
  }

  // Fill Charged particles QA:
  template <typename T>
  void fillChargedQAHistos(T const& track, int nFT0M)
  {
    hist.fill(HIST("QA/after/h_Eta"), track.eta());
    hist.fill(HIST("QA/after/h_Phi"), track.phi());
    hist.fill(HIST("QA/after/h_Pt"), track.pt());
    hist.fill(HIST("QA/after/h2_Pt_NFT0M"), track.pt(), nFT0M);
    hist.fill(HIST("QA/after/h2_PvsPinner"), track.p(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/h2_Pt_Eta"), track.eta(), track.pt());
    hist.fill(HIST("QA/after/h_DcaZ"), track.dcaZ());
    hist.fill(HIST("QA/after/h_DcaXY"), track.dcaXY());
    hist.fill(HIST("QA/after/h2_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST("QA/after/h2_DcaZ"), track.pt(), track.dcaZ());

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
  void moments(double pt, double weight, double& Q1, double& Q2, double& Q3, double& Q4)
  {
    Q1 += pt * weight;
    Q2 += pt * pt * weight;
    Q3 += pt * pt * pt * weight;
    Q4 += pt * pt * pt * pt * weight;
  }

  float getCorrectedWeight(double pt, double rap, bool cfgEffPtCh, bool cfgEffPtPi, bool cfgEffPtKa, bool cfgEffPtPr, bool cfgEffPtYPi, bool cfgEffPtYKa, bool cfgEffPtYPr)
  {
    float weight = 1.0;
    if (cfgEffPtCh) {
      float effPt = heffPt->GetBinContent(heffPt->FindBin(pt));
      weight = 1.0 / effPt;
    } else if (cfgEffPtPi) {
      float effPtPi = heffPtPi->GetBinContent(heffPtPi->FindBin(pt));
      weight = 1.0 / effPtPi;
    } else if (cfgEffPtKa) {
      float effPtKa = heffPtKa->GetBinContent(heffPtKa->FindBin(pt));
      weight = 1.0 / effPtKa;
    } else if (cfgEffPtPr) {
      float effPtPr = heffPtPr->GetBinContent(heffPtPr->FindBin(pt));
      weight = 1.0 / effPtPr;
    } else if (cfgEffPtYPi) {
      float effPtYPi = heffPtRapPi->GetBinContent(heffPtRapPi->FindBin(rap, pt));
      weight = 1.0 / effPtYPi;
    } else if (cfgEffPtYKa) {
      float effPtYKa = heffPtRapKa->GetBinContent(heffPtRapKa->FindBin(rap, pt));
      weight = 1.0 / effPtYKa;
    } else if (cfgEffPtYPr) {
      float effPtYPr = heffPtRapPr->GetBinContent(heffPtRapPr->FindBin(rap, pt));
      weight = 1.0 / effPtYPr;
    } else {
      weight = 1.0;
    }

    return weight;
  }

  // Fill after PID cut QA hist:
  template <int Mode, typename T>
  void fillIdParticleQAHistos(T const& track, double rap, double nSigmaTPC, double nSigmaTOF, int nFT0M, bool cfgEffPtPi, bool cfgEffPtKa, bool cfgEffPtPr, bool cfgEffPtYPi, bool cfgEffPtYKa, bool cfgEffPtYPr, int& N, int& NW, double& Q1, double& Q2, double& Q3, double& Q4)
  {
    double pt = track.pt();
    double weight = getCorrectedWeight(pt, rap, false, cfgEffPtPi, cfgEffPtKa, cfgEffPtPr, cfgEffPtYPi, cfgEffPtYKa, cfgEffPtYPr);

    NW += weight;
    N++;
    moments(pt, weight, Q1, Q2, Q3, Q4);

    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Rap_weighted"), rap, pt, weight);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Pt_weighted"), pt, weight);

    hist.fill(HIST(Dire[Mode]) + HIST("h_Pt"), track.pt());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_NFT0M"), track.pt(), nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Eta"), track.eta(), track.pt());
    if (track.sign() > 0) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtPos"), track.pt());
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtPos_NFT0M"), track.pt(), nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtPos_rap"), rap, track.pt());
    }
    if (track.sign() < 0) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtNeg"), track.pt());
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtNeg_NFT0M"), track.pt(), nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtNeg_rap"), rap, track.pt());
    }

    hist.fill(HIST(Dire[Mode]) + HIST("h_Eta"), track.eta());
    hist.fill(HIST(Dire[Mode]) + HIST("h_Phi"), track.phi());
    hist.fill(HIST(Dire[Mode]) + HIST("h_Rap"), rap);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_Pt_Rap"), rap, track.pt());
    hist.fill(HIST(Dire[Mode]) + HIST("h_DcaZ"), track.dcaZ());
    hist.fill(HIST(Dire[Mode]) + HIST("h_DcaXY"), track.dcaXY());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_DcaZ"), track.pt(), track.dcaZ());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_DcaXY"), track.pt(), track.dcaXY());

    hist.fill(HIST(Dire[Mode]) + HIST("h2_TPCNsigma_El"), track.p(), track.tpcNSigmaEl());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TOFNsigma_El"), track.p(), track.tofNSigmaEl());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TPCNsigma"), track.p(), nSigmaTPC);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TOFNsigma"), track.p(), nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TpcTofNsigma"), nSigmaTPC, nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST(Dire[Mode]) + HIST("h2_pvsm2"), track.mass() * track.mass(), track.p());
    hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/after/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TPCNsigma_El"), track.tpcInnerParam(), track.tpcNSigmaEl());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TOFNsigma_El"), track.tpcInnerParam(), track.tofNSigmaEl());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TPCNsigma"), track.tpcInnerParam(), nSigmaTPC);
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TOFNsigma"), track.tpcInnerParam(), nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TpcTofNsigma"), nSigmaTPC, nSigmaTOF);
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST(Dire[Mode]) + HIST("innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/after/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/after/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
  }

  template <int Mode>
  void fillPtMCHist(double pt, double eta, double rap, int nFT0M, int pid, int pdgCodePos, int pdgCodeNeg)
  {
    hist.fill(HIST(Dire[Mode]) + HIST("h_PtTruth"), pt);
    hist.fill(HIST(Dire[Mode]) + HIST("h_EtaTruth"), eta);
    hist.fill(HIST(Dire[Mode]) + HIST("h_RapTruth"), rap);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_PtTruth_NFT0M"), pt, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h2_PtTruth_Rap"), rap, pt);

    if (pid == pdgCodePos) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtPosTruth"), pt);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtPosTruth_NFT0M"), pt, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtPosTruth_Rap"), rap, pt);
    }
    if (pid == pdgCodeNeg) {
      hist.fill(HIST(Dire[Mode]) + HIST("h_PtNegTruth"), pt);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtNegTruth_NFT0M"), pt, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h2_PtNegTruth_Rap"), rap, pt);
    }
  }

  template <int Mode>
  void fillAnalysisHistos(int nTPC, int nFT0M, int N, int NW, double Q1, double Q2, double Q3, double Q4)
  {
    if (N == 0) {
      return;
    }
    double twopart1 = ((Q1 * Q1) - Q2);
    double threepart1 = ((Q1 * Q1 * Q1) - (3 * Q2 * Q1) + 2 * Q3);
    double fourpart1 = ((Q1 * Q1 * Q1 * Q1) - (6 * Q2 * Q1 * Q1) + (3 * Q2 * Q2) + (8 * Q3 * Q1) - 6 * Q4);

    hist.fill(HIST(Dire[Mode]) + HIST("h_Mult"), N);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Mult_weighted"), NW);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q1"), nTPC, Q1, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q2"), nTPC, Q2, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q3"), nTPC, Q3, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q4"), nTPC, Q4, nFT0M);

    if (N > 1) {
      double meanPt = Q1 / static_cast<double>(NW);
      double nPair = (static_cast<double>(NW) * (static_cast<double>(NW) - 1));
      double twopart = twopart1 / nPair;
      double checkNDenoVar = (1 / std::sqrt(1 - (1 / static_cast<double>(NW))));
      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT"), meanPt);
      hist.fill(HIST(Dire[Mode]) + HIST("p_mean_pT_Mult_var"), nTPC, meanPt);

      hist.fill(HIST(Dire[Mode]) + HIST("h_Q1_var"), nTPC, Q1, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_N_var"), nTPC, N, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_nume_Mult_var"), nTPC, twopart1, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_deno_Mult_var"), nTPC, nPair, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Mult_var"), nTPC, meanPt, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Mult_var"), nTPC, twopart, nFT0M);
      hist.fill(HIST(Dire[Mode]) + HIST("p_CheckNCh"), nTPC, checkNDenoVar);
      hist.fill(HIST(Dire[Mode]) + HIST("h_CheckNCh"), nTPC, checkNDenoVar, nFT0M);

      if (N > 2) {
        double nTriplet = (static_cast<double>(N) * (static_cast<double>(N) - 1) * (static_cast<double>(N) - 2));
        double threepart = threepart1 / nTriplet;
        hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Mult_skew"), nTPC, meanPt, nFT0M);
        hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Mult_skew"), nTPC, twopart, nFT0M);
        hist.fill(HIST(Dire[Mode]) + HIST("h_threepart_Mult_skew"), nTPC, threepart, nFT0M);

        if (N > 3) {
          double nQuad = (static_cast<double>(N) * (static_cast<double>(N) - 1) * (static_cast<double>(N) - 2) * (static_cast<double>(N) - 3));
          double fourpart = fourpart1 / nQuad;
          hist.fill(HIST(Dire[Mode]) + HIST("h_mean_pT_Mult_kurto"), nTPC, meanPt, nFT0M);
          hist.fill(HIST(Dire[Mode]) + HIST("h_twopart_Mult_kurto"), nTPC, twopart, nFT0M);
          hist.fill(HIST(Dire[Mode]) + HIST("h_threepart_Mult_kurto"), nTPC, threepart, nFT0M);
          hist.fill(HIST(Dire[Mode]) + HIST("h_fourpart_Mult_kurto"), nTPC, fourpart, nFT0M);
        }
      }
    }
  }

  template <bool DataFlag, bool RecoFlag, typename T, typename U>
  void fillHistos(T const& col, U const& tracks)
  {
    int nCh = 0, nTPC = 0, nFT0M = 0, nFT0C = 0;
    int nChW = 0;

    int nPi = 0, nKa = 0, nPr = 0;
    int nPiW = 0, nKaW = 0, nPrW = 0;
    double ptCh = 0, q1Ch = 0, q2Ch = 0, q3Ch = 0, q4Ch = 0;
    double ptPi = 0, q1Pi = 0, q2Pi = 0, q3Pi = 0, q4Pi = 0;
    double ptPr = 0, q1Pr = 0, q2Pr = 0, q3Pr = 0, q4Pr = 0;
    double ptKa = 0, q1Ka = 0, q2Ka = 0, q3Ka = 0, q4Ka = 0;

    int nChSim = 0, nSim = 0, nFT0CSim = 0;
    int nPiSim = 0, nKaSim = 0, nPrSim = 0;
    double eta = 0, etaSim = -999, rapSim = -999;
    double ptChSim = 0, q1ChSim = 0, q2ChSim = 0, q3ChSim = 0, q4ChSim = 0;
    double ptPiSim = 0, q1PiSim = 0, q2PiSim = 0, q3PiSim = 0, q4PiSim = 0;
    double ptPrSim = 0, q1PrSim = 0, q2PrSim = 0, q3PrSim = 0, q4PrSim = 0;
    double ptKaSim = 0, q1KaSim = 0, q2KaSim = 0, q3KaSim = 0, q4KaSim = 0;

    double wghtCh = 1.0, wghtPi = 1.0, wghtKa = 1.0, wghtPr = 1.0;

    array<float, 3> p1, p2;
    double invMassGamma = 0.0;

    for (const auto& [trkEl, trkPos] : soa::combinations(soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (trkEl.index() == trkPos.index())
        continue;

      if (!selTrack(trkEl) || !selTrack(trkPos))
        continue;

      if (!selElectrons(trkEl) || !selElectrons(trkPos))
        continue;

      p1 = std::array{trkEl.px(), trkEl.py(), trkEl.pz()};
      p2 = std::array{trkPos.px(), trkPos.py(), trkPos.pz()};

      invMassGamma = RecoDecay::m(std::array{p1, p2}, std::array{MassElectron, MassElectron});
      hist.fill(HIST("QA/after/h_invMass_gamma"), invMassGamma);
    }

    fillAfterQAHistos(col);

    if constexpr (DataFlag) {
      nTPC = col.multNTracksHasTPC();
      nFT0M = col.multFT0M();
      nFT0C = col.multFT0C();

      for (const auto& track : tracks) {
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

        if (std::fabs(track.eta()) < 0.8) {
          ptCh = track.pt();
          double weight = getCorrectedWeight(ptCh, 0.0, cfgEffPtCh, false, false, false, false, false, false);
          nChW += weight;
          nCh++;
          hist.fill(HIST("QA/after/h_Pt_weighted"), ptCh, weight);
          moments(ptCh, weight, q1Ch, q2Ch, q3Ch, q4Ch);

          fillChargedQAHistos(track, nFT0M);
        }

        fillBeforePIDQAHistos(track);

        if (rejectTracks(track)) {
          return;
        }

        if (cfgInvMass == true && invMassGamma < cfgGammaCut) {
          continue;
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPi(track) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, cfgEffPtPi, false, false, cfgEffPtYPi, false, false, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPi(track) == cfgSelLow && selHighPi(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, cfgEffPtPi, false, false, cfgEffPtYPi, false, false, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, false, cfgEffPtKa, false, false, cfgEffPtYKa, false, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowKa(track) == cfgSelLow && selHighKa(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, false, cfgEffPtKa, false, false, cfgEffPtYKa, false, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, false, false, cfgEffPtPr, false, false, cfgEffPtYPr, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPr(track) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, false, false, cfgEffPtPr, false, false, cfgEffPtYPr, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
          }
        }
      }
    } else if constexpr (RecoFlag) {
      if (!col.has_mcCollision()) {
        LOGF(warning, "No MC collision for this collision, skip...");
        return;
      }
      nTPC = col.multNTracksHasTPC();
      nFT0M = col.multFT0M();
      nFT0C = col.multFT0C();

      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC Particle for this track, skip...");
          continue;
        }
        auto mcPart = track.mcParticle();
        int pid = mcPart.pdgCode();
        if (!mcPart.isPhysicalPrimary()) {
          continue;
        }

        if (std::abs(track.eta()) < 0.8) {
          nTPC++;
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

        //______________________________Reconstructed Level____________________________________________________//

        if (selTrack(track)) {

          if (std::fabs(track.eta()) < 0.8) {
            ptCh = track.pt();
            double weight = getCorrectedWeight(ptCh, 0.0, cfgEffPtCh, false, false, false, false, false, false);
            hist.fill(HIST("QA/after/h_Pt_weighted"), ptCh, weight);
            nChW += weight;
            nCh++;
            moments(ptCh, weight, q1Ch, q2Ch, q3Ch, q4Ch);
            fillChargedQAHistos(track, nFT0M);
          }
          fillBeforePIDQAHistos(track);

          if (cfgRejTrk == true && rejectTracks(track)) {
            return;
          }

          if (cfgInvMass == true && invMassGamma < cfgGammaCut) {
            continue;
          }

          eta = track.eta();
          if (cfgPDGCodeOnly == true) {
            if (std::abs(pid) == kPiPlus && std::abs(rapPi) < cfgCutRap && track.pt() >= cfgCutPiPtMin) {
              ptPi = track.pt();
              fillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, cfgEffPtPi, false, false, cfgEffPtYPi, false, false, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
              fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
            }

            if (std::abs(pid) == kKPlus && std::abs(rapKa) < cfgCutRap && track.pt() >= cfgCutKaPtMin) {
              ptKa = track.pt();
              fillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, false, cfgEffPtKa, false, false, cfgEffPtYKa, false, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
              fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
            }

            if (std::abs(pid) == kProton && std::abs(rapPr) < cfgCutRap && track.pt() >= cfgCutPrPtMin) {
              ptPr = track.pt();
              fillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, false, false, cfgEffPtPr, false, false, cfgEffPtYPr, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
              fillPtMCHist<QA_Proton>(ptPr, eta, rapPr, nFT0M, pid, kProton, kProtonBar);
            }
          }

          if (cfgPidCut == true) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowPi(track) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
                ptPi = track.pt();
                fillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, cfgEffPtPi, false, false, cfgEffPtYPi, false, false, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
                if (std::abs(pid) == kPiPlus) {
                  fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowPi(track) == cfgSelLow && selHighPi(track) == cfgSelHigh) {
                ptPi = track.pt();
                fillIdParticleQAHistos<QA_Pion>(track, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, cfgEffPtPi, false, false, cfgEffPtYPi, false, false, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
                if (std::abs(pid) == kPiPlus) {
                  fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
                }
              }
            }

            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
                ptKa = track.pt();
                fillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, false, cfgEffPtKa, false, false, cfgEffPtYKa, false, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
                if (std::abs(pid) == kKPlus) {
                  fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowKa(track) == cfgSelLow && selHighKa(track) == cfgSelHigh) {
                ptKa = track.pt();
                fillIdParticleQAHistos<QA_Kaon>(track, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, false, cfgEffPtKa, false, false, cfgEffPtYKa, false, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
                if (std::abs(pid) == kKPlus) {
                  fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
                }
              }
            }

            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
                ptPr = track.pt();
                fillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, false, false, cfgEffPtPr, false, false, cfgEffPtYPr, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
                if (std::abs(pid) == kProton) {
                  fillPtMCHist<QA_Proton>(ptPr, eta, rapPr, nFT0M, pid, kProton, kProtonBar);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowPr(track) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
                ptPr = track.pt();
                fillIdParticleQAHistos<QA_Proton>(track, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, false, false, cfgEffPtPr, false, false, cfgEffPtYPr, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
                if (std::abs(pid) == kProton) {
                  fillPtMCHist<QA_Proton>(ptPr, eta, rapPr, nFT0M, pid, kProton, kProtonBar);
                }
              }
            }
          }
        }

        //___________________________________Truth Level____________________________________________________//
        auto charge = 0.;
        auto* pd = pdg->GetParticle(pid);
        if (pd != nullptr) {
          charge = pd->Charge();
        }
        if (std::fabs(charge) < 1e-3) {
          continue;
        }
        if (std::abs(pid) != kElectron && std::abs(pid) != kMuonMinus && std::abs(pid) != kPiPlus && std::abs(pid) != kKPlus && std::abs(pid) != kProton) {
          continue;
        }

        if (std::fabs(mcPart.eta()) < 0.8) {
          nSim++;
        }

        if (mcPart.eta() > -3.3 || mcPart.eta() < -2.1) {
          nFT0CSim++;
        }

        if (mcPart.pt() > cfgCutPtMin && mcPart.pt() < cfgCutPtMax) {

          if (std::abs(mcPart.eta()) < 0.8) {
            nChSim++;
            ptChSim = mcPart.pt();
            moments(ptChSim, 1.0, q1ChSim, q2ChSim, q3ChSim, q4ChSim);
            hist.fill(HIST("Gen/Charged/h_PtTruth"), mcPart.pt());
            hist.fill(HIST("Gen/Charged/h2_PtTruth_NFT0M"), mcPart.pt(), nFT0M);
            hist.fill(HIST("Gen/Charged/h2_PtTruth_Eta"), mcPart.eta(), mcPart.pt());
            hist.fill(HIST("Gen/Charged/h_EtaTruth"), mcPart.eta());
            hist.fill(HIST("Gen/Charged/h_PhiTruth"), mcPart.phi());
          }

          if (std::abs(mcPart.y()) > cfgCutRap) {
            continue;
          }

          if (std::abs(pid) == kPiPlus && mcPart.pt() >= cfgCutPiPtMin) {
            etaSim = mcPart.eta();
            rapSim = mcPart.y();

            if (cfgSelOR == true && cfgSelAND == false) {
              if (mcPart.p() <= cfgCutPiThrsldP || mcPart.p() > cfgCutPiThrsldP) {
                nPiSim++;
                ptPiSim = mcPart.pt();
                moments(ptPiSim, 1.0, q1PiSim, q2PiSim, q3PiSim, q4PiSim);
                fillPtMCHist<Gen_Pion>(ptPiSim, etaSim, rapSim, nFT0M, pid, kPiPlus, kPiMinus);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPiThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutPiThrsldP)) {
                nPiSim++;
                ptPiSim = mcPart.pt();
                moments(ptPiSim, 1.0, q1PiSim, q2PiSim, q3PiSim, q4PiSim);
                fillPtMCHist<Gen_Pion>(ptPiSim, etaSim, rapSim, nFT0M, pid, kPiPlus, kPiMinus);
              }
            }
            hist.fill(HIST("Gen/Pion/h_PhiTruth"), mcPart.phi());
          }

          if (std::abs(pid) == kKPlus && mcPart.pt() >= cfgCutKaPtMin) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPiThrsldP) || (cfgSelHigh == true && mcPart.p() > cfgCutPiThrsldP)) {
                nKaSim++;
                ptKaSim = mcPart.pt();
                moments(ptKaSim, 1.0, q1KaSim, q2KaSim, q3KaSim, q4KaSim);
                fillPtMCHist<Gen_Kaon>(ptKaSim, etaSim, rapSim, nFT0M, pid, kKPlus, kKMinus);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutKaThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutKaThrsldP)) {
                nKaSim++;
                ptKaSim = mcPart.pt();
                moments(ptKaSim, 1.0, q1KaSim, q2KaSim, q3KaSim, q4KaSim);
                fillPtMCHist<Gen_Kaon>(ptKaSim, etaSim, rapSim, nFT0M, pid, kKPlus, kKMinus);
              }
            }
            hist.fill(HIST("Gen/Kaon/h_PhiTruth"), mcPart.phi());
          }

          if (std::abs(pid) == kProton && mcPart.pt() >= cfgCutPrPtMin) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPrThrsldP) || (cfgSelHigh == true && mcPart.p() > cfgCutPrThrsldP)) {
                nPrSim++;
                ptPrSim = mcPart.pt();
                moments(ptPrSim, 1.0, q1PrSim, q2PrSim, q3PrSim, q4PrSim);
                fillPtMCHist<Gen_Proton>(ptPrSim, etaSim, rapSim, nFT0M, pid, kProton, kProtonBar);
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if ((cfgSelLow == true && mcPart.p() <= cfgCutPrThrsldP) && (cfgSelHigh == true && mcPart.p() > cfgCutPrThrsldP)) {
                nPrSim++;
                ptPrSim = mcPart.pt();
                moments(ptPrSim, 1.0, q1PrSim, q2PrSim, q3PrSim, q4PrSim);
                fillPtMCHist<Gen_Proton>(ptPrSim, etaSim, rapSim, nFT0M, pid, kProton, kProtonBar);
              }
            }
            hist.fill(HIST("Gen/Proton/h_PhiTruth"), mcPart.phi());
          }
        }
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////

      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC Particle for this track, skip...");
          continue;
        }
        auto mcPart = track.mcParticle();
        int pid = mcPart.pdgCode();
        if (!mcPart.isPhysicalPrimary()) {
          continue;
        }

        double rapPi = track.rapidity(MassPiPlus);
        double rapKa = track.rapidity(MassKPlus);
        double rapPr = track.rapidity(MassProton);

        if (selTrack(track)) {
          if (std::abs(track.eta()) < 0.8) {
            double pt = track.pt();
            wghtCh = getCorrectedWeight(pt, 0.0, cfgEffPtCh, false, false, false, false, false, false);
            hist.fill(HIST("QA/after/h2_pt_nch"), nCh, pt, wghtCh);
            hist.fill(HIST("QA/after/h3_nft0m_pt_nch"), nCh, pt, nFT0M, wghtCh);
            hist.fill(HIST("QA/after/h2_pt_nch_prof"), nCh, pt, wghtCh);
          }

          if (selLowPi(track) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
            ptPi = track.pt();
            wghtPi = getCorrectedWeight(ptPi, rapPi, false, cfgEffPtPi, false, false, cfgEffPtYPi, false, false);
            hist.fill(HIST("QA/Pion/h2_pt_nch"), nPi, ptPi, wghtPi);
            hist.fill(HIST("QA/Pion/h3_nft0m_pt_nch"), nPi, ptPi, nFT0M, wghtPi);
            hist.fill(HIST("QA/Pion/h2_pt_nch_prof"), nPi, ptPi, wghtPi);
          }

          if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
            ptKa = track.pt();
            wghtKa = getCorrectedWeight(ptKa, rapKa, false, false, cfgEffPtKa, false, false, cfgEffPtYKa, false);
            hist.fill(HIST("QA/Kaon/h2_pt_nch"), nKa, ptKa, wghtKa);
            hist.fill(HIST("QA/Kaon/h3_nft0m_pt_nch"), nKa, ptKa, nFT0M, wghtKa);
            hist.fill(HIST("QA/Kaon/h2_pt_nch_prof"), nKa, ptKa, wghtKa);
          }

          if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
            ptPr = track.pt();
            wghtPr = getCorrectedWeight(ptPr, rapPr, false, false, false, cfgEffPtPr, false, false, cfgEffPtYPr);
            hist.fill(HIST("QA/Proton/h2_pt_nch"), nPr, ptPr, wghtPr);
            hist.fill(HIST("QA/Proton/h3_nft0m_pt_nch"), nPr, ptPr, nFT0M, wghtPr);
            hist.fill(HIST("QA/Proton/h2_pt_nch_prof"), nPr, ptPr, wghtPr);
          }
        }

        auto charge = 0.;
        auto* pd = pdg->GetParticle(pid);
        if (pd != nullptr) {
          charge = pd->Charge();
        }
        if (std::fabs(charge) < 1e-3) {
          continue;
        }
        if (std::abs(pid) != kElectron && std::abs(pid) != kMuonMinus && std::abs(pid) != kPiPlus && std::abs(pid) != kKPlus && std::abs(pid) != kProton) {
          continue;
        }
        if (mcPart.pt() > cfgCutPtMin && mcPart.pt() < cfgCutPtMax) {
          if (std::abs(mcPart.eta()) < 0.8) {
            double pt = mcPart.pt();
            hist.fill(HIST("Gen/Charged/h2_pt_nch"), nChSim, pt);
            hist.fill(HIST("Gen/Charged/h3_nft0m_pt_nch"), nChSim, pt, nFT0M);
            hist.fill(HIST("Gen/Charged/h2_pt_nch_prof"), nChSim, pt);
          }

          if (std::abs(mcPart.y()) >= 0.5)
            continue;

          if (std::abs(pid) == kPiPlus && mcPart.pt() >= cfgCutPiPtMin) {
            hist.fill(HIST("Gen/Pion/h2_pt_nch"), nPiSim, mcPart.pt());
            hist.fill(HIST("Gen/Pion/h3_nft0m_pt_nch"), nPiSim, mcPart.pt(), nFT0M);
            hist.fill(HIST("Gen/Pion/h2_pt_nch_prof"), nPiSim, mcPart.pt());
          }
          if (std::abs(pid) == kKPlus && mcPart.pt() >= cfgCutKaPtMin) {
            hist.fill(HIST("Gen/Kaon/h2_pt_nch"), nKaSim, mcPart.pt());
            hist.fill(HIST("Gen/Kaon/h3_nft0m_pt_nch"), nKaSim, mcPart.pt(), nFT0M);
            hist.fill(HIST("Gen/Kaon/h2_pt_nch_prof"), nKaSim, mcPart.pt());
          }
          if (std::abs(pid) == kProton && mcPart.pt() >= cfgCutPrPtMin) {
            hist.fill(HIST("Gen/Proton/h2_pt_nch"), nPrSim, mcPart.pt());
            hist.fill(HIST("Gen/Proton/h3_nft0m_pt_nch"), nPrSim, mcPart.pt(), nFT0M);
            hist.fill(HIST("Gen/Proton/h2_pt_nch_prof"), nPrSim, mcPart.pt());
          }
        }
      }

      hist.fill(HIST("Gen/h_Counts"), 2);
      hist.fill(HIST("QA/after/h_VtxZReco"), col.posZ());
      hist.fill(HIST("Gen/h_VtxZ"), col.mcCollision().posZ());

      if (nSim > 0)
        hist.fill(HIST("Gen/h_NSim"), nSim);

      if (nSim > 0 && nChSim > 0)
        hist.fill(HIST("Gen/h2_NChSim_NSim"), nSim, nChSim);

      if (nSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NSim"), nSim, nTPC);

      if (nChSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NChSim"), nTPC, nChSim);

      if (nPiSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NPiSim"), nTPC, nPiSim);

      if (nKaSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NKaSim"), nTPC, nKaSim);

      if (nPrSim > 0 && nTPC > 0)
        hist.fill(HIST("Gen/h2_NTPC_NPrSim"), nTPC, nPrSim);

      if (nChSim > 0 && nCh > 0)
        hist.fill(HIST("Gen/Charged/h2_Nid_NidSim"), nChSim, nCh, wghtCh);

      if (nPi > 0 && nPiSim > 0)
        hist.fill(HIST("Gen/Pion/h2_Nid_NidSim"), nPiSim, nPi, wghtPi);

      if (nKa > 0 && nKaSim > 0)
        hist.fill(HIST("Gen/Kaon/h2_Nid_NidSim"), nKaSim, nKa, wghtKa);

      if (nPr > 0 && nPrSim > 0)
        hist.fill(HIST("Gen/Proton/h2_Nid_NidSim"), nPrSim, nPr, wghtPr);

      hist.fill(HIST("Gen/h_NTPC"), nTPC);
      hist.fill(HIST("Gen/h_NFT0C"), nFT0CSim);
      hist.fill(HIST("Gen/h2_NTPC_NFT0C"), nFT0CSim, nTPC);
      hist.fill(HIST("Gen/h2_NTPC_NFT0M"), nFT0M, nTPC);

      if (nFT0C != 0 && nFT0CSim != 0)
        hist.fill(HIST("Gen/h2_NFT0C_NFT0CSim"), nFT0CSim, nFT0C);

      fillAnalysisHistos<Gen_Charged>(nTPC, nFT0M, nChSim, nChSim, q1ChSim, q2ChSim, q3ChSim, q4ChSim);
      fillAnalysisHistos<Gen_Pion>(nTPC, nFT0M, nPiSim, nPiSim, q1PiSim, q2PiSim, q3PiSim, q4PiSim);
      fillAnalysisHistos<Gen_Kaon>(nTPC, nFT0M, nKaSim, nKaSim, q1KaSim, q2KaSim, q3KaSim, q4KaSim);
      fillAnalysisHistos<Gen_Proton>(nTPC, nFT0M, nPrSim, nPrSim, q1PrSim, q2PrSim, q3PrSim, q4PrSim);
    }

    if (nTPC > 0 && nCh > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NCh"), nTPC, nCh, wghtCh);

    if (nPi > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NPi"), nTPC, nPi, wghtPi);

    if (nKa > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NKa"), nTPC, nKa, wghtKa);

    if (nPr > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NPr"), nTPC, nPr, wghtPr);

    fillAnalysisHistos<Analysis_Charged>(nTPC, nFT0M, nCh, nChW, q1Ch, q2Ch, q3Ch, q4Ch);
    fillAnalysisHistos<Analysis_Pion>(nTPC, nFT0M, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
    fillAnalysisHistos<Analysis_Kaon>(nTPC, nFT0M, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
    fillAnalysisHistos<Analysis_Proton>(nTPC, nFT0M, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
  }

  void processRun3(MyRun3Collisions::iterator const& col, MyAllTracks const& tracks)
  {
    // Before Collision and Track Cuts:
    fillBeforeQAHistos(col, tracks);

    // After Collision and Track Cuts:
    if (selRun3Col(col)) {
      fillHistos<true, false>(col, tracks);
    }
  }
  PROCESS_SWITCH(MeanPtFlucId, processRun3, "Process for Run-3", false);

  void processMCRecoSimRun3(MyRun3MCCollisions::iterator const& col, aod::McCollisions const&, MyMCTracks const& tracks, aod::McParticles const&)
  {
    // Before Collision and Track Cuts:
    fillBeforeQAHistos(col, tracks);

    hist.fill(HIST("Gen/h_VtxZ_b"), col.mcCollision().posZ());

    // After Collision and Track Cuts:
    if (selRun3Col(col)) {
      fillHistos<false, true>(col, tracks);
    }
  }
  PROCESS_SWITCH(MeanPtFlucId, processMCRecoSimRun3, "process MC Reconstructed & Truth Run-3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MeanPtFlucId>(cfgc)};
}
