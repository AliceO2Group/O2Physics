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
// #include "CCDB/BasicCCDBManager.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
double massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
double massPr = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

struct meanPtFlucId {
  Configurable<int> nPtBins{"nPtBins", 300, ""};
  Configurable<int> nPartBins{"nPartBins", 250, ""};
  Configurable<int> nCentBins{"nCentBins", 101, ""};
  Configurable<int> nEtaBins{"nEtaBins", 100, ""};
  Configurable<float> cfgCutPtMax{"cfgCutPtMax", 2.0, "maximum pT"};
  Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.15, "minimum pT"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut"};
  Configurable<float> cfgCutRap{"cfgCutRap", 0.5, "Rapidity Cut"};
  Configurable<float> cfgCutDcaXY{"cfgCutDcaXY", 0.12, "DCAxy cut"};
  Configurable<float> cfgCutDcaZ{"cfgCutDcaZ", 1.0, "DCAz cut"};
  Configurable<float> cfgCutPosZ{"cfgCutPosZ", 10.0, "cut for vertex Z"};
  Configurable<float> cfgCutNSigTpcEl{"cfgCutNSigTpcEl", 1.5, "TPC nSigma Electron veto cut"};
  Configurable<float> cfgCutNSigTofEl{"cfgCutNSigTofEl", 1.5, "TOF nSigma Electron veto cut"};
  Configurable<float> cfgCutNSig2{"cfgCutNSig2", 2.0, "nSigma cut (2)"};
  Configurable<float> cfgCutNSig3{"cfgCutNSig3", 3.0, "nSigma cut (3)"};
  Configurable<float> cfgCutPiPtMin{"cfgCutPiPtMin", 0.2, "Minimum pion p_{T} cut"};
  Configurable<float> cfgCutKaPtMin{"cfgCutKaPtMin", 0.3, "Minimum kaon p_{T} cut"};
  Configurable<float> cfgCutPrPtMin{"cfgCutPrPtMin", 0.5, "Minimum proton p_{T} cut"};
  Configurable<float> cfgCutPiP1{"cfgCutPiP1", 0.65, "pion p cut-1"};
  Configurable<float> cfgCutPiP2{"cfgCutPiP2", 0.75, "pion p cut-2"};
  // Configurable<float> cfgCutPiP3{"cfgCutPiP3", 1.20, "pion p cut-3"};
  Configurable<float> cfgCutKaP1{"cfgCutKaP1", 0.50, "kaon p cut-1"};
  Configurable<float> cfgCutKaP2{"cfgCutKaP2", 0.60, "kaon p cut-2"};
  Configurable<float> cfgCutKaP3{"cfgCutKaP3", 1.60, "kaon p cut-3"};
  Configurable<float> cfgCutPrP1{"cfgCutPrP1", 0.90, "proton p cut-1"};
  Configurable<float> cfgCutPrP2{"cfgCutPrP2", 1.00, "proton p cut-2"};
  Configurable<float> cfgMCPi{"cfgMCPi", 1.0, "Adjust Pi nSigmaTPC for MC data"};
  Configurable<float> cfgMCKa{"cfgMCKa", 1.4, "Adjust Ka nSigmaTPC for MC data"};
  Configurable<float> cfgMCPr{"cfgMCPr", 1.7, "Adjust Pr nSigmaTPC for MC data"};
  Configurable<bool> cfgSelPi{"cfgSelPi", true, "PID selection cut for Pions"};
  Configurable<bool> cfgSelKa{"cfgSelKa", true, "PID selection cut for Kaons"};
  Configurable<bool> cfgSelPr{"cfgSelPr", true, "PID selection cut for Protons"};
  Configurable<bool> cfgSelPiInnerParam{"cfgSelPiInnerParam", false, "PID selection cut for Pions by using Momentum at inner wall of the TPC"};
  Configurable<bool> cfgSelKaInnerParam{"cfgSelKaInnerParam", false, "PID selection cut for Kaons by using Momentum at inner wall of the TPC"};
  Configurable<bool> cfgSelPrInnerParam{"cfgSelPrInnerParam", false, "PID selection cut for Protons by using Momentum at inner wall of the TPC"};
  ConfigurableAxis multTPCBins{"multTPCBins", {150, 0, 150}, "TPC Multiplicity bins"};
  ConfigurableAxis multFT0CBins{"multFT0CBins", {150, 0, 1500}, "Forward Multiplicity bins"};
  ConfigurableAxis dcaXYBins{"dcaXYBins", {100, -0.15, 0.15}, "dcaXY bins"};
  ConfigurableAxis dcaZBins{"dcaZBins", {100, -1.2, 1.2}, "dcaZ bins"};

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

  // Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    const AxisSpec axisEvents{5, 0, 5, "Counts"};
    const AxisSpec axisEta{nEtaBins, -1., +1., "#eta"};
    const AxisSpec axisY{nEtaBins, -1., +1., "Rapidity"};
    const AxisSpec axisPt{nPtBins, 0., 3., "p_{T} (GeV/c)"};
    const AxisSpec axisP{nPtBins, 0., 3., "p (GeV/c)"};
    const AxisSpec axisInnerParam{nPtBins, 0., 3., "p_{InnerParam } (GeV/c)"};
    const AxisSpec axisPart{nPartBins, 0., 20., " "};
    const AxisSpec axisDeno{100, 1., 2.0, "#frac{1}{#sqrt{1 - #frac{1}{N}}}"};
    const AxisSpec axisMeanPt{100, 0., 3., "M(p_{T}) (GeV/c)"};
    const AxisSpec axisMult{100, 0, 100, "N_{ch}"};
    const AxisSpec axisMultTPC{multTPCBins, "N_{TPC} "};
    const AxisSpec axisMultFT0C{multFT0CBins, "N_{FT0C}"};
    const AxisSpec axisCentFT0C{nCentBins, 0, 101, "FT0C (%)"};
    const AxisSpec axisVtxZ{80, -20., 20., "V_{Z} (cm)"};
    const AxisSpec axisDCAz{dcaZBins, "DCA_{Z} (cm)"};
    const AxisSpec axisDCAxy{dcaXYBins, "DCA_{XY} (cm)"};
    const AxisSpec axisTPCNsigma{500, -5., 5., "n #sigma_{TPC}"};
    const AxisSpec axisTOFNsigma{500, -5., 5., "n #sigma_{TOF}"};
    const AxisSpec axisTPCSignal{180, 20., 200., "#frac{dE}{dx}"};
    const AxisSpec axisTOFSignal{200, 0.2, 1.2, "TOF #beta"};
    const AxisSpec axisChi2{40, 0., 40., "Chi2"};
    const AxisSpec axisCrossedTPC{300, 0, 300, "Crossed TPC"};
    const AxisSpec axisM2{100, 0., 1.4, "#it{m}^{2} (GeV/#it{c}^{2})^{2}"};

    HistogramConfigSpec QnHist({HistType::kTHnSparseD, {axisMultTPC, axisPart, axisMultFT0C}});
    HistogramConfigSpec DenoHist({HistType::kTHnSparseD, {axisMultTPC, axisDeno, axisMultFT0C}});
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
    hist.add("QA/before/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/before/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/before/h_DcaZ", "DCA_{Z}", kTH1D, {axisDCAz});
    hist.add("QA/before/h_DcaXY", "DCA_{XY}", kTH1D, {axisDCAxy});
    hist.add("QA/before/h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/before/h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/before/h_NTPC", "N_{TPC}", kTH1D, {axisMultTPC});
    hist.add("QA/before/h_NFT0C", "FT0C Multiplicity", kTH1D, {axisMultFT0C});
    hist.add("QA/before/h_Cent", "FT0C (%)", kTH1D, {axisCentFT0C});
    hist.add("QA/before/h2_NTPC_Cent", "N_{TPC} vs FT0C(%)", kTH2D, {{axisCentFT0C}, {axisMultTPC}});
    hist.add("QA/before/h2_NTPC_NFT0C", "N_{TPC} vs N_{FT0C}", kTH2D, {{axisMultFT0C}, {axisMultTPC}});

    hist.add("QA/before/h2_TPCSignal", "TPC Signal", TPCSignalHist);
    hist.add("QA/before/h2_TOFSignal", "TOF Signal", TOFSignalHist);
    hist.add("QA/before/h2_pvsm2", "p vs m^{2}", PvsM2Hist);

    hist.add("QA/before/innerParam/h2_TPCSignal", "TPC Signal", TPCSignalHist1);
    hist.add("QA/before/innerParam/h2_TOFSignal", "TOF Signal", TOFSignalHist1);
    hist.add("QA/before/innerParam/h2_pvsm2", "p vs m^{2}", PvsM2Hist1);

    hist.addClone("QA/before/", "QA/after/");

    hist.add("QA/after/p_NTPC_NFT0C", "N_{TPC} vs N_{FT0C} (Profile)", kTProfile, {{axisMultFT0C}});
    hist.add("QA/after/p_NTPC_Cent", "N_{TPC} vs FT0C(%) (Profile)", kTProfile, {{axisCentFT0C}});
    hist.add("QA/after/h2_NTPC_Nch", "N_{ch} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMult}});

    hist.add("QA/Pion/h_Pt", "p_{T} (TPC & TPC+TOF)", kTH1D, {axisPt});
    hist.add("QA/Pion/h_rap", "y (TPC & TPC+TOF)", kTH1D, {axisY});
    hist.add("QA/Pion/h_Eta", "Pseudorapidity (TPC & TPC+TOF)", kTH1D, {axisY});
    hist.add("QA/Pion/h2_Pt_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/Pion/h_DcaZ", "DCA_{z}", kTH1D, {axisDCAz});
    hist.add("QA/Pion/h_DcaXY", "DCA_{xy}", kTH1D, {axisDCAxy});
    hist.add("QA/Pion/h2_DcaZ", "DCA_{z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/Pion/h2_DcaXY", "DCA_{xy}", kTH2D, {{axisPt}, {axisDCAxy}});

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
    hist.add("QA/Pion/h2_ExpTPCSignal", "Expected TPC Signal", TPCSignalHist);
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
    hist.add("QA/Pion/innerParam/h2_ExpTPCSignal", "Expected TPC Signal", TPCSignalHist1);
    hist.add("QA/Pion/innerParam/h2_pvsm2", "p vs m^{2}", PvsM2Hist1);

    hist.addClone("QA/Pion/", "QA/Kaon/");
    hist.addClone("QA/Pion/", "QA/Proton/");

    // Analysis Plots:
    hist.add("Analysis/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Analysis/Charged/h_mean_Q1", " <p_{T}> ", kTH1D, {axisMeanPt});
    hist.add("Analysis/Charged/p_mean_Q1_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/p_CheckNCH", " 1/denominator vs N_{ch} ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/h_CheckNCH", " 1/denominator vs N_{ch} ", DenoHist);
    hist.add("Analysis/Charged/h_mean_Q1_Mult_var", " <p_{T}> vs N_{ch} ", QnHist);
    hist.add("Analysis/Charged/h_mean_Q1_Mult_skew", " <p_{T}> vs N_{ch} ", QnHist);
    hist.add("Analysis/Charged/h_mean_Q1_Mult_kurto", " <p_{T}> vs N_{ch} ", QnHist);
    hist.add("Analysis/Charged/p_twopart_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
    hist.add("Analysis/Charged/h_twopart_Mult_var", "Twopart vs N_{ch} ", QnHist);
    hist.add("Analysis/Charged/h_twopart_Mult_skew", "Twopart vs N_{ch} ", QnHist);
    hist.add("Analysis/Charged/h_twopart_Mult_kurto", "Twopart vs N_{ch} ", QnHist);
    hist.add("Analysis/Charged/h_threepart_Mult_skew", "Threepart vs N_{ch} ", QnHist);
    hist.add("Analysis/Charged/h_threepart_Mult_kurto", "Threepart vs N_{ch} ", QnHist);
    hist.add("Analysis/Charged/h_fourpart_Mult_kurto", "Fourpart vs N_{ch} ", QnHist);

    hist.addClone("Analysis/Charged/", "Analysis/Pion/");
    hist.addClone("Analysis/Charged/", "Analysis/Kaon/");
    hist.addClone("Analysis/Charged/", "Analysis/Proton/");

    hist.add("QA/Reco/Pion/h_allPt", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Reco/Kaon/h_allPt", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Reco/Proton/h_allPt", "p_{T} ", kTH1D, {axisPt});

    // MC Generated
    hist.add("Gen/Counts", "Counts", kTH1D, {axisEvents});
    hist.add("Gen/vtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/NTPC", "Mid rapidity Multiplicity", kTH1D, {axisMultTPC});
    hist.add("Gen/NFT0C", "Forward Multiplicity", kTH1D, {axisMultFT0C});
    hist.add("Gen/Charged/h_Pt", "p_{T} ", kTH1D, {axisPt});

    hist.addClone("Analysis/Charged/", "Gen/Charged/");

    hist.addClone("Gen/Charged/", "Gen/Pion/");
    hist.addClone("Gen/Charged/", "Gen/Kaon/");
    hist.addClone("Gen/Charged/", "Gen/Proton/");
  }

  // Event selection cuts:
  template <typename T>
  bool selRun3Col(T const& col)
  {
    if (std::abs(col.posZ()) > cfgCutPosZ)
      return false;

    if (!col.sel8())
      return false;

    if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return false;

    // if (!col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
    //   return false;

    if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;

    if (!col.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return false;

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

    if (std::abs(track.dcaZ()) > cfgCutDcaZ)
      return false;

    if (std::abs(track.dcaXY()) > cfgCutDcaXY)
      return false;

    return true;
  }

  // Cuts to rejects the tracks
  template <typename T>
  bool rejectTracks(T const& track)
  {
    if (track.tpcNSigmaEl() > -3. && track.tpcNSigmaEl() < 5. && (track.tpcNSigmaPi() - cfgMCPi) < -3. && (track.tpcNSigmaKa() - cfgMCKa) < -3 && (track.tpcNSigmaPr() - cfgMCPr) < -3 && (track.tpcNSigmaPi() - cfgMCPi) > 3 && (track.tpcNSigmaKa() - cfgMCKa) > 3 && (track.tpcNSigmaPr() - cfgMCPr) > 3) {
      return true;
    }
    return false;
  }

  // PID selction cuts for Pions
  template <typename T>
  bool selPions(T const& track)
  {
    if ((!track.hasTOF() && std::abs(track.tpcNSigmaEl()) > cfgCutNSigTpcEl &&
         ((std::abs(track.tpcNSigmaPi() - cfgMCPi) < cfgCutNSig3 && track.pt() >= cfgCutPiPtMin && track.p() <= cfgCutPiP1) || (std::abs(track.tpcNSigmaPi() - cfgMCPi) < cfgCutNSig2 && track.p() > cfgCutPiP1 && track.p() <= cfgCutPiP2))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaPi() - cfgMCPi) < cfgCutNSig3 && std::abs(track.tofNSigmaEl()) > cfgCutNSigTofEl && (std::abs(track.tofNSigmaPi()) < cfgCutNSig3 && track.pt() >= cfgCutPiPtMin))) {
      if (abs(track.rapidity(massPi)) < cfgCutRap)
        return true;
    }
    return false;
  }

  // PID selction cuts for Kaons
  template <typename T>
  bool selKaons(T const& track)
  {
    if (((!track.hasTOF()) && std::abs(track.tpcNSigmaEl()) > cfgCutNSigTpcEl &&
         ((std::abs(track.tpcNSigmaKa() - cfgMCKa) < cfgCutNSig3 && track.pt() >= cfgCutKaPtMin && track.p() <= cfgCutKaP1) || (std::abs(track.tpcNSigmaKa() - cfgMCKa) < cfgCutNSig2 && track.p() > cfgCutKaP1 && track.p() <= cfgCutKaP2))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaKa() - cfgMCKa) < cfgCutNSig3 && std::abs(track.tofNSigmaEl()) > cfgCutNSigTofEl &&
         ((std::abs(track.tofNSigmaKa()) < cfgCutNSig3 && track.pt() >= cfgCutKaPtMin && track.p() <= cfgCutKaP3) || (std::abs(track.tofNSigmaKa()) < cfgCutNSig2 && track.p() > cfgCutKaP3)))) {
      if (abs(track.rapidity(massKa)) < cfgCutRap)
        return true;
    }

    return false;
  }

  // PID selction cuts for Protons
  template <typename T>
  bool selProtons(T const& track)
  {
    if (((!track.hasTOF()) && std::abs(track.tpcNSigmaEl()) > cfgCutNSigTpcEl &&
         ((std::abs(track.tpcNSigmaPr() - cfgMCPr) < cfgCutNSig3 && track.pt() >= cfgCutPrPtMin && track.p() <= cfgCutPrP1) || (std::abs(track.tpcNSigmaPr() - cfgMCPr) < cfgCutNSig2 && track.p() > cfgCutPrP1 && track.p() <= cfgCutPrP2))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaPr() - cfgMCPr) < cfgCutNSig3 && std::abs(track.tofNSigmaEl()) > cfgCutNSigTofEl && std::abs(track.tofNSigmaPr()) < cfgCutNSig3 && track.pt() > cfgCutPrPtMin)) {
      if (abs(track.rapidity(massPr)) < cfgCutRap)
        return true;
    }

    return false;
  }

  // PID selction cuts for Pions by using momentum at inner wall of TPC
  template <typename T>
  bool selPionsInnerParam(T const& track)
  {
    if ((!track.hasTOF() && std::abs(track.tpcNSigmaEl()) > cfgCutNSigTpcEl &&
         ((std::abs(track.tpcNSigmaPi() - cfgMCPi) < cfgCutNSig3 && track.pt() >= cfgCutPiPtMin && track.tpcInnerParam() <= cfgCutPiP1) || (std::abs(track.tpcNSigmaPi() - cfgMCPi) < cfgCutNSig2 && track.tpcInnerParam() > cfgCutPiP1 && track.tpcInnerParam() <= cfgCutPiP2))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaPi() - cfgMCPi) < cfgCutNSig3 && std::abs(track.tofNSigmaEl()) > cfgCutNSigTofEl && (std::abs(track.tofNSigmaPi()) < cfgCutNSig3 && track.pt() >= cfgCutPiPtMin))) {
      if (abs(track.rapidity(massPi)) < cfgCutRap)
        return true;
    }
    return false;
  }

  // PID selction cuts for Kaons by using momentum at inner wall of TPC
  template <typename T>
  bool selKaonsInnerParam(T const& track)
  {
    if (((!track.hasTOF()) && std::abs(track.tpcNSigmaEl()) > cfgCutNSigTpcEl &&
         ((std::abs(track.tpcNSigmaKa() - cfgMCKa) < cfgCutNSig3 && track.pt() >= cfgCutKaPtMin && track.tpcInnerParam() <= cfgCutKaP1) || (std::abs(track.tpcNSigmaKa() - cfgMCKa) < cfgCutNSig2 && track.tpcInnerParam() > cfgCutKaP1 && track.tpcInnerParam() <= cfgCutKaP2))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaKa() - cfgMCKa) < cfgCutNSig3 && std::abs(track.tofNSigmaEl()) > cfgCutNSigTofEl &&
         ((std::abs(track.tofNSigmaKa()) < cfgCutNSig3 && track.pt() >= cfgCutKaPtMin && track.tpcInnerParam() <= cfgCutKaP3) || (std::abs(track.tofNSigmaKa()) < cfgCutNSig2 && track.tpcInnerParam() > cfgCutKaP3)))) {
      if (abs(track.rapidity(massKa)) < cfgCutRap)
        return true;
    }

    return false;
  }

  // PID selction cuts for Protons by using momentum at inner wall of TPC
  template <typename T>
  bool selProtonsInnerParam(T const& track)
  {
    if (((!track.hasTOF()) && std::abs(track.tpcNSigmaEl()) > cfgCutNSigTpcEl &&
         ((std::abs(track.tpcNSigmaPr() - cfgMCPr) < cfgCutNSig3 && track.pt() >= cfgCutPrPtMin && track.tpcInnerParam() <= cfgCutPrP1) || (std::abs(track.tpcNSigmaPr() - cfgMCPr) < cfgCutNSig2 && track.tpcInnerParam() > cfgCutPrP1 && track.tpcInnerParam() <= cfgCutPrP2))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaPr() - cfgMCPr) < cfgCutNSig3 && std::abs(track.tofNSigmaEl()) > cfgCutNSigTofEl && std::abs(track.tofNSigmaPr()) < cfgCutNSig3 && track.pt() > cfgCutPrPtMin)) {
      if (abs(track.rapidity(massPr)) < cfgCutRap)
        return true;
    }

    return false;
  }

  // Moments Calculation:
  void moments(double pt, double* Q1, double* Q2, double* Q3, double* Q4)
  {
    *Q1 += pt;
    *Q2 += pt * pt;
    *Q3 += pt * pt * pt;
    *Q4 += pt * pt * pt * pt;
  }

  // Cumulant parts Calculation:
  void parts(double Q1, double Q2, double Q3, double Q4, int N, double* mean_Q1, double* twopart, double* threepart, double* fourpart)
  {
    if (N > 1) {
      *mean_Q1 = Q1 / static_cast<double>(N);
      *twopart = ((Q1 * Q1) - Q2) / (static_cast<double>(N) * (static_cast<double>(N) - 1));
    }
    if (N > 2) {
      *threepart = ((Q1 * Q1 * Q1) - (3 * Q2 * Q1) + 2 * Q3) / (static_cast<double>(N) * (static_cast<double>(N) - 1) * (static_cast<double>(N) - 2));
    }
    if (N > 3) {
      *fourpart = ((Q1 * Q1 * Q1 * Q1) - (6 * Q2 * Q1 * Q1) + (3 * Q2 * Q2) + (8 * Q3 * Q1) - 6 * Q4) / (static_cast<double>(N) * (static_cast<double>(N) - 1) * (static_cast<double>(N) - 2) * (static_cast<double>(N) - 3));
    }
  }

  // Fill hist before selection cuts:
  template <typename T, typename U>
  void FillBeforeQAHistos(T const& col, U const& tracks)
  {
    for (auto& myTrack : tracks) {
      hist.fill(HIST("QA/before/h_Eta"), myTrack.eta());
      hist.fill(HIST("QA/before/h_Pt"), myTrack.pt());
      hist.fill(HIST("QA/before/h2_PvsPinner"), myTrack.p(), myTrack.tpcInnerParam());
      hist.fill(HIST("QA/before/h2_Pt_Eta"), myTrack.eta(), myTrack.pt());
      hist.fill(HIST("QA/before/h_TPCChi2perCluster"), myTrack.tpcChi2NCl());
      hist.fill(HIST("QA/before/h_ITSChi2perCluster"), myTrack.itsChi2NCl());
      hist.fill(HIST("QA/before/h_crossedTPC"), myTrack.tpcNClsCrossedRows());
      hist.fill(HIST("QA/before/h_DcaXY"), myTrack.dcaXY());
      hist.fill(HIST("QA/before/h_DcaZ"), myTrack.dcaZ());
      hist.fill(HIST("QA/before/h2_DcaXY"), myTrack.pt(), myTrack.dcaXY());
      hist.fill(HIST("QA/before/h2_DcaZ"), myTrack.pt(), myTrack.dcaZ());
    }
    hist.fill(HIST("QA/before/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/before/h_Counts"), 2);

    hist.fill(HIST("QA/before/h_NTPC"), col.multTPC());
    hist.fill(HIST("QA/before/h_Cent"), col.centFT0C());
    hist.fill(HIST("QA/before/h_NFT0C"), col.multFT0C());
    hist.fill(HIST("QA/before/h2_NTPC_NFT0C"), col.multFT0C(), col.multTPC());
    hist.fill(HIST("QA/before/h2_NTPC_Cent"), col.centFT0C(), col.multTPC());
  }

  // Fill hist after selection cuts:
  template <typename T>
  void FillAfterQAHistos(T const& col)
  {
    hist.fill(HIST("QA/after/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/after/h_Counts"), 2);
    hist.fill(HIST("QA/after/h_NTPC"), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/h_Cent"), col.centFT0C());
    hist.fill(HIST("QA/after/h_NFT0C"), col.multFT0C());
    hist.fill(HIST("QA/after/h2_NTPC_NFT0C"), col.multFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/h2_NTPC_Cent"), col.centFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/p_NTPC_Cent"), col.centFT0C(), col.multNTracksHasTPC());
    hist.fill(HIST("QA/after/p_NTPC_NFT0C"), col.multFT0C(), col.multNTracksHasTPC());
  }

  // Fill Charged particles QA:
  template <typename T>
  void FillChargedQAHistos(T const& track)
  {
    hist.fill(HIST("QA/after/h_Eta"), track.eta());
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

  // Fill Pion QA hist:
  template <typename T>
  void FillPionQAHistos(T const& track)
  {
    hist.fill(HIST("QA/Pion/h_Pt"), track.pt());
    hist.fill(HIST("QA/Pion/h_Eta"), track.eta());
    hist.fill(HIST("QA/Pion/h_rap"), track.rapidity(massPi));
    hist.fill(HIST("QA/Pion/h2_Pt_rap"), track.rapidity(massPi), track.pt());
    hist.fill(HIST("QA/Pion/h_DcaZ"), track.dcaZ());
    hist.fill(HIST("QA/Pion/h_DcaXY"), track.dcaXY());
    hist.fill(HIST("QA/Pion/h2_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST("QA/Pion/h2_DcaZ"), track.pt(), track.dcaZ());

    hist.fill(HIST("QA/Pion/h2_TPCNsigma_El"), track.p(), track.tpcNSigmaEl());
    hist.fill(HIST("QA/Pion/h2_TOFNsigma_El"), track.p(), track.tofNSigmaEl());
    hist.fill(HIST("QA/Pion/h2_TPCNsigma"), track.p(), track.tpcNSigmaPi());
    hist.fill(HIST("QA/Pion/h2_TOFNsigma"), track.p(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/Pion/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/Pion/h2_ExpTPCSignal"), track.p(), track.tpcExpSignalPi(track.tpcSignal()));
    hist.fill(HIST("QA/Pion/h2_pvsm2"), track.mass() * track.mass(), track.p());
    hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/after/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST("QA/Pion/innerParam/h2_TPCNsigma_El"), track.tpcInnerParam(), track.tpcNSigmaEl());
    hist.fill(HIST("QA/Pion/innerParam/h2_TOFNsigma_El"), track.tpcInnerParam(), track.tofNSigmaEl());
    hist.fill(HIST("QA/Pion/innerParam/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaPi());
    hist.fill(HIST("QA/Pion/innerParam/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/innerParam/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
    hist.fill(HIST("QA/Pion/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/Pion/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/Pion/innerParam/h2_ExpTPCSignal"), track.tpcInnerParam(), track.tpcExpSignalPi(track.tpcSignal()));
    hist.fill(HIST("QA/Pion/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/after/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/after/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
  }

  // Fill Kaon QA hist:
  template <typename T>
  void FillKaonQAHistos(T const& track)
  {
    hist.fill(HIST("QA/Kaon/h_Pt"), track.pt());
    hist.fill(HIST("QA/Kaon/h_Eta"), track.eta());
    hist.fill(HIST("QA/Kaon/h_rap"), track.rapidity(massKa));
    hist.fill(HIST("QA/Kaon/h2_Pt_rap"), track.rapidity(massKa), track.pt());
    hist.fill(HIST("QA/Kaon/h_DcaZ"), track.dcaZ());
    hist.fill(HIST("QA/Kaon/h_DcaXY"), track.dcaXY());
    hist.fill(HIST("QA/Kaon/h2_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST("QA/Kaon/h2_DcaZ"), track.pt(), track.dcaZ());

    hist.fill(HIST("QA/Kaon/h2_TPCNsigma_El"), track.p(), track.tpcNSigmaEl());
    hist.fill(HIST("QA/Kaon/h2_TOFNsigma_El"), track.p(), track.tofNSigmaEl());
    hist.fill(HIST("QA/Kaon/h2_TPCNsigma"), track.p(), track.tpcNSigmaKa());
    hist.fill(HIST("QA/Kaon/h2_TOFNsigma"), track.p(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/Kaon/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/Kaon/h2_ExpTPCSignal"), track.p(), track.tpcExpSignalKa(track.tpcSignal()));
    hist.fill(HIST("QA/Kaon/h2_pvsm2"), track.mass() * track.mass(), track.p());
    hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/after/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST("QA/Kaon/innerParam/h2_TPCNsigma_El"), track.tpcInnerParam(), track.tpcNSigmaEl());
    hist.fill(HIST("QA/Kaon/innerParam/h2_TOFNsigma_El"), track.tpcInnerParam(), track.tofNSigmaEl());
    hist.fill(HIST("QA/Kaon/innerParam/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaKa());
    hist.fill(HIST("QA/Kaon/innerParam/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/innerParam/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    hist.fill(HIST("QA/Kaon/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/Kaon/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/Kaon/innerParam/h2_ExpTPCSignal"), track.tpcInnerParam(), track.tpcExpSignalKa(track.tpcSignal()));
    hist.fill(HIST("QA/Kaon/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/after/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/after/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
  }

  // Fill Proton QA hist:
  template <typename T>
  void FillProtonQAHistos(T const& track)
  {
    hist.fill(HIST("QA/Proton/h_Pt"), track.pt());
    hist.fill(HIST("QA/Proton/h_Eta"), track.eta());
    hist.fill(HIST("QA/Proton/h_rap"), track.rapidity(massPr));
    hist.fill(HIST("QA/Proton/h2_Pt_rap"), track.rapidity(massPr), track.pt());
    hist.fill(HIST("QA/Proton/h_DcaZ"), track.dcaZ());
    hist.fill(HIST("QA/Proton/h_DcaXY"), track.dcaXY());
    hist.fill(HIST("QA/Proton/h2_DcaZ"), track.pt(), track.dcaZ());
    hist.fill(HIST("QA/Proton/h2_DcaXY"), track.pt(), track.dcaXY());

    hist.fill(HIST("QA/Proton/h2_TPCNsigma_El"), track.p(), track.tpcNSigmaEl());
    hist.fill(HIST("QA/Proton/h2_TOFNsigma_El"), track.p(), track.tofNSigmaEl());
    hist.fill(HIST("QA/Proton/h2_TPCNsigma"), track.p(), track.tpcNSigmaPr());
    hist.fill(HIST("QA/Proton/h2_TOFNsigma"), track.p(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/Proton/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/Proton/h2_ExpTPCSignal"), track.p(), track.tpcExpSignalPr(track.tpcSignal()));
    hist.fill(HIST("QA/Proton/h2_pvsm2"), track.mass() * track.mass(), track.p());
    hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
    hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
    hist.fill(HIST("QA/after/h2_pvsm2"), track.mass() * track.mass(), track.p());

    hist.fill(HIST("QA/Proton/innerParam/h2_TPCNsigma_El"), track.tpcInnerParam(), track.tpcNSigmaEl());
    hist.fill(HIST("QA/Proton/innerParam/h2_TOFNsigma_El"), track.tpcInnerParam(), track.tofNSigmaEl());
    hist.fill(HIST("QA/Proton/innerParam/h2_TPCNsigma"), track.tpcInnerParam(), track.tpcNSigmaPr());
    hist.fill(HIST("QA/Proton/innerParam/h2_TOFNsigma"), track.tpcInnerParam(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/innerParam/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
    hist.fill(HIST("QA/Proton/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/Proton/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/Proton/innerParam/h2_ExpTPCSignal"), track.tpcInnerParam(), track.tpcExpSignalPr(track.tpcSignal()));
    hist.fill(HIST("QA/Proton/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
    hist.fill(HIST("QA/after/innerParam/h2_TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST("QA/after/innerParam/h2_TOFSignal"), track.tpcInnerParam(), track.beta());
    hist.fill(HIST("QA/after/innerParam/h2_pvsm2"), track.mass() * track.mass(), track.tpcInnerParam());
  }

  template <bool DataFlag, bool RecoFlag, typename T, typename U>
  void FillHistos(T const& col, U const& tracks)
  {
    int N_Pi = 0, N_Ka = 0, N_Pr = 0;
    int Nch = 0, NTPC = 0, N_FT0C = 0;
    // double Cent_FT0C = 0;
    double pt_ch = 0, Q1_ch = 0, Q2_ch = 0, Q3_ch = 0, Q4_ch = 0;
    double pt_Pi = 0, Q1_Pi = 0, Q2_Pi = 0, Q3_Pi = 0, Q4_Pi = 0;
    double pt_Pr = 0, Q1_Pr = 0, Q2_Pr = 0, Q3_Pr = 0, Q4_Pr = 0;
    double pt_Ka = 0, Q1_Ka = 0, Q2_Ka = 0, Q3_Ka = 0, Q4_Ka = 0;
    double mean_Q1_Ch = 0, mean_Q1_Pi = 0, mean_Q1_Ka = 0, mean_Q1_Pr = 0;
    double twopart_Ch = 0, twopart_Pi = 0, twopart_Ka = 0, twopart_Pr = 0;
    double threepart_Ch = 0, threepart_Pi = 0, threepart_Ka = 0, threepart_Pr = 0;
    double fourpart_Ch = 0, fourpart_Pi = 0, fourpart_Ka = 0, fourpart_Pr = 0;

    for (auto& track : tracks) {
      if (!selTrack(track)) {
        continue;
      }

      if constexpr (DataFlag) {
        if (std::abs(track.eta()) < 0.8) {
          Nch++;
          pt_ch = track.pt();
          moments(pt_ch, &Q1_ch, &Q2_ch, &Q3_ch, &Q4_ch);
          FillChargedQAHistos(track);
        }

        FillBeforePIDQAHistos(track);

        if (rejectTracks(track)) {
          return;
        }

        // For Pions:
        if (selPions(track) == cfgSelPi /* && selPionsInnerParam(track) == cfgSelPiInnerParam */) {
          N_Pi++;
          pt_Pi = track.pt();
          moments(pt_Pi, &Q1_Pi, &Q2_Pi, &Q3_Pi, &Q4_Pi);
          FillPionQAHistos(track);
        }

        // For Kaons:
        if (selKaons(track) == cfgSelKa /* && selKaonsInnerParam(track) == cfgSelKaInnerParam */) {
          N_Ka++;
          pt_Ka = track.pt();
          moments(pt_Ka, &Q1_Ka, &Q2_Ka, &Q3_Ka, &Q4_Ka);
          FillKaonQAHistos(track);
        }

        // For Protons:
        if (selProtons(track) == cfgSelPr /* && selProtonsInnerParam(track) == cfgSelPrInnerParam */) {
          N_Pr++;
          pt_Pr = track.pt();
          moments(pt_Pr, &Q1_Pr, &Q2_Pr, &Q3_Pr, &Q4_Pr);
          FillProtonQAHistos(track);
        }

      } else if constexpr (RecoFlag) {
        if (track.has_mcParticle() && track.mcParticle().isPhysicalPrimary()) {
          if (std::abs(track.eta()) < 0.8) {
            Nch++;
            pt_ch = track.pt();
            moments(pt_ch, &Q1_ch, &Q2_ch, &Q3_ch, &Q4_ch);
            FillChargedQAHistos(track);
          }
          FillBeforePIDQAHistos(track);

          if (rejectTracks(track))
            return;

          if (selPions(track) == cfgSelPi /* || selPionsInnerParam(track) == cfgSelPiInnerParam */) {
            hist.fill(HIST("QA/Reco/Pion/h_allPt"), track.pt());
            if (std::abs(track.mcParticle().pdgCode()) == 211) {
              N_Pi++;
              pt_Pi = track.pt();
              moments(pt_Pi, &Q1_Pi, &Q2_Pi, &Q3_Pi, &Q4_Pi);
              FillPionQAHistos(track);
            }
          }

          if (selKaons(track) == cfgSelKa /* || selKaonsInnerParam(track) == cfgSelKaInnerParam */) {
            hist.fill(HIST("QA/Reco/Kaon/h_allPt"), track.pt());
            if (std::abs(track.mcParticle().pdgCode()) == 321) {
              N_Ka++;
              pt_Ka = track.pt();
              moments(pt_Ka, &Q1_Ka, &Q2_Ka, &Q3_Ka, &Q4_Ka);
              FillKaonQAHistos(track);
            }
          }

          if (selProtons(track) == cfgSelPr /* || selProtonsInnerParam(track) == cfgSelPrInnerParam */) {
            hist.fill(HIST("QA/Reco/Proton/h_allPt"), track.pt());
            if (std::abs(track.mcParticle().pdgCode()) == 2212) {
              N_Pr++;
              pt_Pr = track.pt();
              moments(pt_Pr, &Q1_Pr, &Q2_Pr, &Q3_Pr, &Q4_Pr);
              FillProtonQAHistos(track);
            }
          }
        }
      }
    }

    N_FT0C = col.multFT0C();
    // Cent_FT0C = col.centFT0C();
    NTPC = col.multNTracksHasTPC();

    FillAfterQAHistos(col);
    hist.fill(HIST("QA/after/h2_NTPC_Nch"), NTPC, Nch);

    static constexpr std::string_view dire[] = {"Analysis/Charged/", "Analysis/Pion/", "Analysis/Kaon/", "Analysis/Proton/"};

    hist.fill(HIST(dire[0]) + HIST("h_Mult"), Nch);
    hist.fill(HIST(dire[1]) + HIST("h_Mult"), N_Pi);
    hist.fill(HIST(dire[2]) + HIST("h_Mult"), N_Ka);
    hist.fill(HIST(dire[3]) + HIST("h_Mult"), N_Pr);

    if (N_FT0C > 0) {
      parts(Q1_ch, Q2_ch, Q3_ch, Q4_ch, Nch, &mean_Q1_Ch, &twopart_Ch, &threepart_Ch, &fourpart_Ch);
      if (Nch > 1) {
        hist.fill(HIST(dire[0]) + HIST("p_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(Nch)))));
        hist.fill(HIST(dire[0]) + HIST("h_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(Nch)))), N_FT0C);
        if (mean_Q1_Ch != 0) {
          hist.fill(HIST(dire[0]) + HIST("h_mean_Q1"), mean_Q1_Ch);
          hist.fill(HIST(dire[0]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Ch);
          hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Ch, N_FT0C);
        }
        if (twopart_Ch != 0) {
          hist.fill(HIST(dire[0]) + HIST("p_twopart_Mult_var"), NTPC, twopart_Ch);
          hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Ch, N_FT0C);
        }
      }
      if (Nch > 2) {
        if (mean_Q1_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Ch, N_FT0C);
        if (twopart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Ch, N_FT0C);
        if (threepart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Ch, N_FT0C);
      }

      if (Nch > 3) {
        if (mean_Q1_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Ch, N_FT0C);
        if (twopart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Ch, N_FT0C);
        if (threepart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Ch, N_FT0C);
        if (fourpart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Ch, N_FT0C);
      }

      parts(Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi, N_Pi, &mean_Q1_Pi, &twopart_Pi, &threepart_Pi, &fourpart_Pi);
      if (N_Pi > 1) {
        hist.fill(HIST(dire[1]) + HIST("p_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Pi)))));
        hist.fill(HIST(dire[1]) + HIST("h_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Pi)))), N_FT0C);
        if (mean_Q1_Pi != 0) {
          hist.fill(HIST(dire[1]) + HIST("h_mean_Q1"), mean_Q1_Pi);
          hist.fill(HIST(dire[1]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Pi);
          hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Pi, N_FT0C);
        }
        if (twopart_Pi != 0) {
          hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Pi, N_FT0C);
        }
      }

      if (N_Pi > 2) {
        if (mean_Q1_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Pi, N_FT0C);
        if (twopart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Pi, N_FT0C);
        if (threepart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Pi, N_FT0C);
      }

      if (N_Pi > 3) {
        if (mean_Q1_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Pi, N_FT0C);
        if (twopart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Pi, N_FT0C);
        if (threepart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Pi, N_FT0C);
        if (fourpart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Pi, N_FT0C);
      }

      parts(Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka, N_Ka, &mean_Q1_Ka, &twopart_Ka, &threepart_Ka, &fourpart_Ka);
      if (N_Ka > 1) {
        hist.fill(HIST(dire[2]) + HIST("p_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Ka)))));
        hist.fill(HIST(dire[2]) + HIST("h_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Ka)))), N_FT0C);
        if (mean_Q1_Ka != 0) {
          hist.fill(HIST(dire[2]) + HIST("h_mean_Q1"), mean_Q1_Ka);
          hist.fill(HIST(dire[2]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Ka);
          hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Ka, N_FT0C);
        }
        if (twopart_Ka != 0) {
          hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Ka, N_FT0C);
        }
      }
      if (N_Ka > 2) {
        if (mean_Q1_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Ka, N_FT0C);
        if (twopart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Ka, N_FT0C);
        if (threepart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Ka, N_FT0C);
      }

      if (N_Ka > 3) {
        if (mean_Q1_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Ka, N_FT0C);
        if (twopart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Ka, N_FT0C);
        if (threepart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Ka, N_FT0C);
        if (fourpart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Ka, N_FT0C);
      }

      parts(Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr, N_Pr, &mean_Q1_Pr, &twopart_Pr, &threepart_Pr, &fourpart_Pr);
      if (N_Pr > 1) {
        hist.fill(HIST(dire[3]) + HIST("p_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Pr)))));
        hist.fill(HIST(dire[3]) + HIST("h_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Pr)))), N_FT0C);
        if (mean_Q1_Pr != 0) {
          hist.fill(HIST(dire[3]) + HIST("h_mean_Q1"), mean_Q1_Pr);
          hist.fill(HIST(dire[3]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Pr);
          hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Pr, N_FT0C);
        }
        if (twopart_Pr != 0) {
          hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Pr, N_FT0C);
        }
      }

      if (N_Pr > 2) {
        if (mean_Q1_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Pr, N_FT0C);
        if (twopart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Pr, N_FT0C);
        if (threepart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Pr, N_FT0C);
      }

      if (N_Pr > 3) {
        if (mean_Q1_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Pr, N_FT0C);
        if (twopart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Pr, N_FT0C);
        if (threepart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Pr, N_FT0C);
        if (fourpart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Pr, N_FT0C);
      }
    }
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

  void process_MCRecoRun3(MyMCCollisions::iterator const& col, aod::McCollisions const&, MyMCTracks const& tracks, aod::McParticles const&)
  {
    // Before Collision and Track Cuts:
    FillBeforeQAHistos(col, tracks);

    // After Collision and Track Cuts:
    if (selRun3Col(col)) {
      if (!col.has_mcCollision()) {
        return;
      }
      FillHistos<false, true>(col, tracks);
    }
  }
  PROCESS_SWITCH(meanPtFlucId, process_MCRecoRun3, "process MC Reconstructed Run-3", true);

  void process_MCGen(soa::Join<aod::McCollisions, aod::MultsExtraMC>::iterator const& mccol, aod::McParticles const& McParticles)
  {
    int N_Pi = 0, N_Ka = 0, N_Pr = 0;
    int Nch = 0, NTPC = 0, N_FT0C = 0;
    double pt_ch = 0, Q1_ch = 0, Q2_ch = 0, Q3_ch = 0, Q4_ch = 0;
    double pt_Pi = 0, Q1_Pi = 0, Q2_Pi = 0, Q3_Pi = 0, Q4_Pi = 0;
    double pt_Pr = 0, Q1_Pr = 0, Q2_Pr = 0, Q3_Pr = 0, Q4_Pr = 0;
    double pt_Ka = 0, Q1_Ka = 0, Q2_Ka = 0, Q3_Ka = 0, Q4_Ka = 0;
    double mean_Q1_Ch = 0, mean_Q1_Pi = 0, mean_Q1_Ka = 0, mean_Q1_Pr = 0;
    double twopart_Ch = 0, twopart_Pi = 0, twopart_Ka = 0, twopart_Pr = 0;
    double threepart_Ch = 0, threepart_Pi = 0, threepart_Ka = 0, threepart_Pr = 0;
    double fourpart_Ch = 0, fourpart_Pi = 0, fourpart_Ka = 0, fourpart_Pr = 0;

    if (abs(mccol.posZ()) > cfgCutPosZ)
      return;

    for (auto& mcParticle : McParticles) {
      if (!mcParticle.isPhysicalPrimary())
        continue;

      auto charge = 0.;
      auto* p = pdg->GetParticle(mcParticle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 1e-3) {
        continue; // reject neutral particles in counters
      }

      if (mcParticle.pt() > cfgCutPtMin && mcParticle.pt() < cfgCutPtMax && abs(mcParticle.y()) < cfgCutRap) {
        Nch++;
        pt_ch = mcParticle.pt();
        moments(pt_ch, &Q1_ch, &Q2_ch, &Q3_ch, &Q4_ch);
        hist.fill(HIST("Gen/Charged/h_Pt"), mcParticle.pt());

        if (std::abs(mcParticle.pdgCode()) == 211) {
          N_Pi++;
          pt_Pi = mcParticle.pt();
          moments(pt_Pi, &Q1_Pi, &Q2_Pi, &Q3_Pi, &Q4_Pi);
          hist.fill(HIST("Gen/Pion/h_Pt"), mcParticle.pt());
        }

        if (std::abs(mcParticle.pdgCode()) == 321) {
          N_Ka++;
          pt_Ka = mcParticle.pt();
          moments(pt_Ka, &Q1_Ka, &Q2_Ka, &Q3_Ka, &Q4_Ka);
          hist.fill(HIST("Gen/Kaon/h_Pt"), mcParticle.pt());
        }

        if (std::abs(mcParticle.pdgCode()) == 2212) {
          N_Pr++;
          pt_Pr = mcParticle.pt();
          moments(pt_Pr, &Q1_Pr, &Q2_Pr, &Q3_Pr, &Q4_Pr);
          hist.fill(HIST("Gen/Proton/h_Pt"), mcParticle.pt());
        }
      }
    }
    NTPC = mccol.multMCNParticlesEta08();
    N_FT0C = mccol.multMCFT0C();
    hist.fill(HIST("Gen/Counts"), 2);
    hist.fill(HIST("Gen/vtxZ"), mccol.posZ());
    hist.fill(HIST("Gen/NTPC"), NTPC);
    hist.fill(HIST("Gen/NFT0C"), N_FT0C);

    static constexpr std::string_view dire[] = {"Gen/Charged/", "Gen/Pion/", "Gen/Kaon/", "Gen/Proton/"};

    hist.fill(HIST(dire[0]) + HIST("h_Mult"), Nch);
    hist.fill(HIST(dire[1]) + HIST("h_Mult"), N_Pi);
    hist.fill(HIST(dire[2]) + HIST("h_Mult"), N_Ka);
    hist.fill(HIST(dire[3]) + HIST("h_Mult"), N_Pr);

    if (N_FT0C > 0) {
      parts(Q1_ch, Q2_ch, Q3_ch, Q4_ch, Nch, &mean_Q1_Ch, &twopart_Ch, &threepart_Ch, &fourpart_Ch);
      if (Nch > 1) {
        hist.fill(HIST(dire[0]) + HIST("p_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(Nch)))));
        hist.fill(HIST(dire[0]) + HIST("h_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(Nch)))), N_FT0C);
        if (mean_Q1_Ch != 0) {
          hist.fill(HIST(dire[0]) + HIST("h_mean_Q1"), mean_Q1_Ch);
          hist.fill(HIST(dire[0]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Ch);
          hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Ch, N_FT0C);
        }
        if (twopart_Ch != 0) {
          hist.fill(HIST(dire[0]) + HIST("p_twopart_Mult_var"), NTPC, twopart_Ch);
          hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Ch, N_FT0C);
        }
      }
      if (Nch > 2) {
        if (mean_Q1_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Ch, N_FT0C);
        if (twopart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Ch, N_FT0C);
        if (threepart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Ch, N_FT0C);
      }

      if (Nch > 3) {
        if (mean_Q1_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Ch, N_FT0C);
        if (twopart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Ch, N_FT0C);
        if (threepart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Ch, N_FT0C);
        if (fourpart_Ch != 0)
          hist.fill(HIST(dire[0]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Ch, N_FT0C);
      }

      parts(Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi, N_Pi, &mean_Q1_Pi, &twopart_Pi, &threepart_Pi, &fourpart_Pi);
      if (N_Pi > 1) {
        hist.fill(HIST(dire[1]) + HIST("p_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Pi)))));
        hist.fill(HIST(dire[1]) + HIST("h_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Pi)))), N_FT0C);
        if (mean_Q1_Pi != 0) {
          hist.fill(HIST(dire[1]) + HIST("h_mean_Q1"), mean_Q1_Pi);
          hist.fill(HIST(dire[1]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Pi);
          hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Pi, N_FT0C);
        }
        if (twopart_Pi != 0) {
          hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Pi, N_FT0C);
        }
      }

      if (N_Pi > 2) {
        if (mean_Q1_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Pi, N_FT0C);
        if (twopart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Pi, N_FT0C);
        if (threepart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Pi, N_FT0C);
      }

      if (N_Pi > 3) {
        if (mean_Q1_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Pi, N_FT0C);
        if (twopart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Pi, N_FT0C);
        if (threepart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Pi, N_FT0C);
        if (fourpart_Pi != 0)
          hist.fill(HIST(dire[1]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Pi, N_FT0C);
      }

      parts(Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka, N_Ka, &mean_Q1_Ka, &twopart_Ka, &threepart_Ka, &fourpart_Ka);
      if (N_Ka > 1) {
        hist.fill(HIST(dire[2]) + HIST("p_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Ka)))));
        hist.fill(HIST(dire[2]) + HIST("h_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Ka)))), N_FT0C);
        if (mean_Q1_Ka != 0) {
          hist.fill(HIST(dire[2]) + HIST("h_mean_Q1"), mean_Q1_Ka);
          hist.fill(HIST(dire[2]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Ka);
          hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Ka, N_FT0C);
        }
        if (twopart_Ka != 0) {
          hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Ka, N_FT0C);
        }
      }
      if (N_Ka > 2) {
        if (mean_Q1_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Ka, N_FT0C);
        if (twopart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Ka, N_FT0C);
        if (threepart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Ka, N_FT0C);
      }

      if (N_Ka > 3) {
        if (mean_Q1_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Ka, N_FT0C);
        if (twopart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Ka, N_FT0C);
        if (threepart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Ka, N_FT0C);
        if (fourpart_Ka != 0)
          hist.fill(HIST(dire[2]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Ka, N_FT0C);
      }

      parts(Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr, N_Pr, &mean_Q1_Pr, &twopart_Pr, &threepart_Pr, &fourpart_Pr);
      if (N_Pr > 1) {
        hist.fill(HIST(dire[3]) + HIST("p_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Pr)))));
        hist.fill(HIST(dire[3]) + HIST("h_CheckNCH"), NTPC, (1 / std::sqrt(1 - (1 / static_cast<double>(N_Pr)))), N_FT0C);
        if (mean_Q1_Pr != 0) {
          hist.fill(HIST(dire[3]) + HIST("h_mean_Q1"), mean_Q1_Pr);
          hist.fill(HIST(dire[3]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Pr);
          hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Pr, N_FT0C);
        }
        if (twopart_Pr != 0) {
          hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Pr, N_FT0C);
        }
      }

      if (N_Pr > 2) {
        if (mean_Q1_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Pr, N_FT0C);
        if (twopart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Pr, N_FT0C);
        if (threepart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Pr, N_FT0C);
      }

      if (N_Pr > 3) {
        if (mean_Q1_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Pr, N_FT0C);
        if (twopart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Pr, N_FT0C);
        if (threepart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Pr, N_FT0C);
        if (fourpart_Pr != 0)
          hist.fill(HIST(dire[3]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Pr, N_FT0C);
      }
    }
  }

  PROCESS_SWITCH(meanPtFlucId, process_MCGen, "process MC Generated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<meanPtFlucId>(cfgc)};
}
