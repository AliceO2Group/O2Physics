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
/// \brief Calculate EbyE <pt> fluctuations with moments method.
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
  Configurable<int> nPartBins{"nPartBins", 500, ""};
  Configurable<int> nCentBins{"nCentBins", 101, ""};
  Configurable<int> nEtaBins{"nEtaBins", 100, ""};
  Configurable<float> ptMax{"ptMax", 2.0, "maximum pT"};
  Configurable<float> ptMin{"ptMin", 0.15, "minimum pT"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> rapCut{"rapCut", 0.5, "Rapidity Cut"};
  Configurable<float> dcaXYCut{"dcaXYCut", 0.12, "DCAxy cut"};
  Configurable<float> dcaZCut{"dcaZCut", 1.0, "DCAz cut"};
  Configurable<float> posZCut{"posZCut", 10.0, "cut for vertex Z"};
  Configurable<float> nSigCuttpcEl{"nSigCuttpcEl", 1.5, "TPC nSigma Electron veto cut"};
  Configurable<float> nSigCutEl{"nSigCutEl", 1.5, "TOF nSigma Electron veto cut"};
  Configurable<float> nSigCut1{"nSigCut1", 1.0, "nSigma cut (1)"};
  Configurable<float> nSigCut2{"nSigCut2", 2.0, "nSigma cut (2)"};
  Configurable<float> nSigCut3{"nSigCut3", 3.0, "nSigma cut (3)"};
  Configurable<float> nSigCut4{"nSigCut4", 4.0, "nSigma cut (4)"};
  Configurable<float> nSigCut5{"nSigCut5", 5.0, "nSigma cut (5)"};
  Configurable<float> nSigCut15{"nSigCut15", 1.5, "nSigma cut (1.5)"};
  Configurable<float> nSigCut25{"nSigCut25", 2.5, "nSigma cut (2.5)"};
  Configurable<float> piP1{"piP1", 0.65, "pion p (1)"};
  Configurable<float> piP2{"piP2", 0.70, "pion p (2)"};
  Configurable<float> piP3{"piP3", 1.40, "pion p (3)"};
  Configurable<float> piP4{"piP4", 1.70, "pion p (4)"};
  Configurable<float> kaP1{"kaP1", 0.20, "min kaon p (1)"};
  Configurable<float> kaP2{"kaP2", 0.5, "kaon p (2)"};
  Configurable<float> kaP3{"kaP3", 0.55, "kaon p (3)"};
  Configurable<float> kaP4{"kaP4", 0.60, "kaon p (4)"};
  Configurable<float> kaP5{"kaP5", 0.65, "kaon p (5)"};
  Configurable<float> kaP6{"kaP6", 1.10, "kaon p (6)"};
  Configurable<float> kaP7{"kaP7", 1.28, "kaon p (7)"};
  Configurable<float> kaP8{"kaP8", 1.50, "kaon p (8)"};
  Configurable<float> prP1{"prP1", 0.40, "min proton p (1)"};
  Configurable<float> prP2{"prP2", 0.95, "proton p (2)"};
  Configurable<float> prP3{"prP3", 1.00, "proton p (3)"};
  Configurable<float> prP4{"prP4", 1.05, "proton p (4)"};
  Configurable<float> prP5{"prP5", 1.13, "proton p (5)"};
  Configurable<float> prP6{"prP6", 1.18, "proton p (6)"};
  ConfigurableAxis multTPCBins{"multTPCBins", {150, 0, 150}, "TPC Multiplicity bins"};
  ConfigurableAxis multFT0MBins{"multFT0MBins", {150, 0, 1500}, "Forward Multiplicity bins"};
  ConfigurableAxis dcaXYBins{"dcaXYBins", {100, -0.15, 0.15}, "dcaXY bins"};
  ConfigurableAxis dcaZBins{"dcaZBins", {100, -1.2, 1.2}, "dcaZ bins"};

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                                aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                                aod::pidTOFbeta>;
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Cs>;
  using MyMCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Cs, aod::McCollisionLabels>;
  using MyMCTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                               aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                               aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                               aod::pidTOFbeta, aod::McTrackLabels>;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    const AxisSpec axisEvents{5, 0, 5, "Counts"};
    const AxisSpec axisEta{nEtaBins, -1., +1., "#eta"};
    const AxisSpec axisY{nEtaBins, -1., +1., "Rapidity"};
    const AxisSpec axisPt{nPtBins, 0., 3., "p_{T} (GeV/c)"};
    const AxisSpec axisP{nPtBins, 0., 3., "p (GeV/c)"};
    const AxisSpec axisPart{nPartBins, 0., 5., " "};
    const AxisSpec axisMeanPt{100, 0., 3., "M(p_{T}) (GeV/c)"};
    const AxisSpec axisMult{100, 0, 100, "N_{ch}"};
    const AxisSpec axisMultTPC{multTPCBins, "N_{TPC} "};
    const AxisSpec axisMultFT0M{multFT0MBins, "N_{FT0M}"};
    const AxisSpec axisCentFT0M{nCentBins, 0, 101, "FT0M (%)"};
    const AxisSpec axisVtxZ{80, -20., 20., "V_{Z} (cm)"};
    const AxisSpec axisDCAz{dcaZBins, "DCA_{Z} (cm)"};
    const AxisSpec axisDCAxy{dcaXYBins, "DCA_{XY} (cm)"};
    const AxisSpec axisTPCNsigma{500, -5., 5., "n #sigma_{TPC}"};
    const AxisSpec axisTOFNsigma{500, -5., 5., "n #sigma_{TOF}"};
    const AxisSpec axisTPCSignal{180, 20., 200., "#frac{dE}{dx}"};
    const AxisSpec axisTOFSignal{100, 0.2, 1.2, "TOF #beta"};
    const AxisSpec axisChi2{40, 0., 40., "Chi2"};
    const AxisSpec axisCrossedTPC{300, 0, 300, "Crossed TPC"};

    HistogramConfigSpec QnHist({HistType::kTHnSparseD, {axisMultTPC, axisPart, axisCentFT0M}});
    HistogramConfigSpec TOFnSigmaHist({HistType::kTH2D, {axisP, axisTOFNsigma}});
    HistogramConfigSpec TOFSignalHist({HistType::kTH2D, {axisP, axisTOFSignal}});
    HistogramConfigSpec TPCnSigmaHist({HistType::kTH2D, {axisP, axisTPCNsigma}});
    HistogramConfigSpec TPCSignalHist({HistType::kTH2D, {axisP, axisTPCSignal}});
    HistogramConfigSpec TPCTOFHist({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});

    // QA Plots:
    hist.add("QA/before/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("QA/before/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});
    hist.add("QA/before/h_TPCChi2perCluster", "TPC #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/before/h_ITSChi2perCluster", "ITS #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/before/h_crossedTPC", "Crossed TPC", kTH1D, {axisCrossedTPC});
    hist.add("QA/before/h_Pt", "p_{T}", kTH1D, {axisPt});
    hist.add("QA/before/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/before/h2_Pt_Eta", "p_{T} vs #eta ", kTH2D, {{axisEta}, {axisPt}});
    hist.add("QA/before/h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/before/h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/before/h_NTPC", "N_{TPC}", kTH1D, {axisMultTPC});
    hist.add("QA/before/h_NFT0M", "FT0M Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/before/h_Cent", "FT0M (%)", kTH1D, {axisCentFT0M});
    hist.add("QA/before/h2_NTPC_Cent", "N_{TPC} vs FT0M(%)", kTH2D, {{axisCentFT0M}, {axisMultTPC}});
    hist.add("QA/before/h2_NTPC_NFT0M", "N_{TPC} vs N_{FT0M}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});
    hist.add("QA/before/h2_TPCSignal", "TPC Signal", TPCSignalHist);
    hist.add("QA/before/h2_TOFSignal", "TOF Signal", TOFSignalHist);

    hist.addClone("QA/before/", "QA/after/");

    hist.add("QA/after/p_NTPC_NFT0M", "N_{TPC} vs N_{FT0M} (Profile)", kTProfile, {{axisMultFT0M}});
    hist.add("QA/after/p_NTPC_Cent", "N_{TPC} vs FT0M(%) (Profile)", kTProfile, {{axisCentFT0M}});
    hist.add("QA/after/h2_NTPC_Nch", "N_{ch} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMult}});

    hist.add("QA/Pion/h_Pt", "p_{T} (TPC & TPC+TOF)", kTH1D, {axisPt});
    hist.add("QA/Pion/h_rap", "y (TPC & TPC+TOF)", kTH1D, {axisY});
    hist.add("QA/Pion/h_Eta", "Pseudorapidity (TPC & TPC+TOF)", kTH1D, {axisY});
    hist.add("QA/Pion/h2_Pt_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt}});
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
    hist.add("QA/Pion/h2_TPCSignal", "TPC Signal Pions", TPCSignalHist);
    hist.add("QA/Pion/h2_TOFSignal", "TOF Signal Pions", TOFSignalHist);
    hist.add("QA/Pion/h2_ExpTPCSignal", "Expected TPC Signal Pions", TPCSignalHist);

    hist.addClone("QA/Pion/", "QA/Kaon/");
    hist.addClone("QA/Pion/", "QA/Proton/");

    // Analysis Plots:
    hist.add("Analysis/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Analysis/Charged/h_mean_Q1", " <p_{T}> ", kTH1D, {axisMeanPt});
    hist.add("Analysis/Charged/p_mean_Q1_Mult_var", " <p_{T}> ", kTProfile, {axisMultTPC});
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

    hist.add("QA/Reco/Charged/h_Pt", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Reco/Charged/h_Eta", "Pseudorapidity ", kTH1D, {axisEta});
    hist.add("QA/Reco/Pion/h_Pt", "p_{T} ", kTH1D, {axisPt});
    hist.add("QA/Reco/Pion/h_Eta", "Pseudorapidity ", kTH1D, {axisEta});
    hist.add("QA/Reco/Pion/h_y", "Rapidity ", kTH1D, {axisY});
    hist.add("QA/Reco/Pion/h2_tofSignal", "TOF Signal ", TOFSignalHist);
    hist.add("QA/Reco/Pion/h2_tpcSignal", "TPC Signal #frac{dE}{dx}", TPCSignalHist);

    hist.addClone("QA/Reco/Pion/", "QA/Reco/Kaon/");
    hist.addClone("QA/Reco/Pion/", "QA/Reco/Proton/");

    hist.add("Gen/Counts", "Counts", kTH1D, {axisEvents});
    hist.add("Gen/vtxZ", "Vertex Z ", kTH1D, {axisVtxZ});
    hist.add("Gen/NTPC", "Mid rapidity Multiplicity", kTH1D, {axisMultTPC});
    hist.add("Gen/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Gen/Charged/h_Pt", "p_{T} ", kTH1D, {axisPt});
    hist.add("Gen/Charged/h_mean_pt", " <p_{T}> ", kTH1D, {axisMeanPt});
    hist.add("Gen/Charged/h_mean_Q1_Mult_var", " <p_{T}> vs N_{ch} ", kTProfile, {axisMultTPC});
    hist.add("Gen/Charged/h_twopart_Mult_var", "Twopart vs N_{ch} ", kTProfile, {axisMultTPC});

    hist.addClone("Gen/Charged/", "Gen/Pion/");
    hist.addClone("Gen/Charged/", "Gen/Kaon/");
    hist.addClone("Gen/Charged/", "Gen/Proton/");
  }

  template <typename T>
  bool selRun3Col(T const& col)
  {
    if (std::abs(col.posZ()) > posZCut)
      return false;

    if (!col.sel8())
      return false;

    return true;
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrackWoPtEta())
      return false;

    if (track.pt() < ptMin)
      return false;

    if (track.pt() > ptMax)
      return false;

    if (track.sign() == 0)
      return false;

    if (std::abs(track.dcaZ()) > dcaZCut)
      return false;

    if (std::abs(track.dcaXY()) > dcaXYCut)
      return false;

    return true;
  }

  template <typename T>
  bool selPions(T const& track)
  {
    if (((!track.hasTOF()) && std::abs(track.tpcNSigmaEl()) > nSigCuttpcEl &&
         ((std::abs(track.tpcNSigmaPi()) < nSigCut3 && track.p() <= piP1) || (std::abs(track.tpcNSigmaPi()) < nSigCut2 && track.p() > piP1 && track.p() <= piP2))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaPi()) < nSigCut4 && std::abs(track.tofNSigmaEl()) > nSigCutEl &&
         ((std::abs(track.tofNSigmaPi()) < nSigCut3 && track.p() <= piP3) || (std::abs(track.tofNSigmaPi()) < nSigCut25 && track.p() > piP3 && track.p() <= piP4) || (std::abs(track.tofNSigmaPi()) < nSigCut2 && track.p() > piP4)))) {
      if (abs(track.rapidity(massPi)) < rapCut)
        return true;
    }

    return false;
  }

  template <typename T>
  bool selKaons(T const& track)
  {
    if (((!track.hasTOF()) && std::abs(track.tpcNSigmaEl()) > nSigCuttpcEl &&
         ((std::abs(track.tpcNSigmaKa()) < nSigCut3 && track.pt() > kaP1 && track.p() <= kaP2) || (std::abs(track.tpcNSigmaKa()) < nSigCut25 && track.p() > kaP2 && track.p() <= kaP3) || (std::abs(track.tpcNSigmaKa()) < nSigCut2 && track.p() > kaP3 && track.p() <= kaP4) || (std::abs(track.tpcNSigmaKa()) < nSigCut15 && track.p() > kaP4 && track.p() <= kaP5))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaKa()) < nSigCut4 && std::abs(track.tofNSigmaEl()) > nSigCutEl &&
         ((std::abs(track.tofNSigmaKa()) < nSigCut3 && track.pt() > kaP1 && track.p() <= kaP6) || (std::abs(track.tofNSigmaKa()) < nSigCut2 && track.p() > kaP6 && track.p() <= kaP7) || (std::abs(track.tofNSigmaKa()) < nSigCut15 && track.p() > kaP7 && track.p() <= kaP8) || (std::abs(track.tofNSigmaKa()) < nSigCut1 && track.p() > kaP8)))) {
      if (abs(track.rapidity(massKa)) < rapCut)
        return true;
    }

    return false;
  }

  template <typename T>
  bool selProtons(T const& track)
  {
    if (((!track.hasTOF()) && std::abs(track.tpcNSigmaEl()) > nSigCuttpcEl &&
         ((std::abs(track.tpcNSigmaPr()) < nSigCut3 && track.pt() > prP1 && track.p() <= prP2) || (std::abs(track.tpcNSigmaPr()) < nSigCut25 && track.p() > prP2 && track.p() <= prP3) || (std::abs(track.tpcNSigmaPr()) < nSigCut2 && track.p() > prP3 && track.p() <= prP4) || (std::abs(track.tpcNSigmaPr()) < nSigCut15 && track.p() > prP4 && track.p() <= prP5) || (std::abs(track.tpcNSigmaPr()) < nSigCut1 && track.p() > prP5 && track.p() <= prP6))) ||
        (track.hasTOF() && std::abs(track.tpcNSigmaPr()) < nSigCut4 && std::abs(track.tofNSigmaEl()) > nSigCutEl && std::abs(track.tofNSigmaPr()) < nSigCut3 && track.pt() > prP1)) {
      if (abs(track.rapidity(massPr)) < rapCut)
        return true;
    }

    return false;
  }

  void moments(double pt, double* Q1, double* Q2, double* Q3, double* Q4)
  {
    *Q1 += pt;
    *Q2 += pt * pt;
    *Q3 += pt * pt * pt;
    *Q4 += pt * pt * pt * pt;
  }

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

  template <bool DataFlag, bool RecoFlag, typename T, typename U>
  void FillHistos(T const& col, U const& tracks)
  {
    int N_Pi = 0, N_Ka = 0, N_Pr = 0;
    int Nch = 0, NTPC = 0, N_FT0C = 0;
    double Cent_FT0C = 0;
    double pt_ch = 0, Q1_ch = 0, Q2_ch = 0, Q3_ch = 0, Q4_ch = 0;
    double pt_Pi = 0, Q1_Pi = 0, Q2_Pi = 0, Q3_Pi = 0, Q4_Pi = 0;
    double pt_Pr = 0, Q1_Pr = 0, Q2_Pr = 0, Q3_Pr = 0, Q4_Pr = 0;
    double pt_Ka = 0, Q1_Ka = 0, Q2_Ka = 0, Q3_Ka = 0, Q4_Ka = 0;
    double mean_Q1_Ch, mean_Q1_Pi, mean_Q1_Ka, mean_Q1_Pr;
    double twopart_Ch, twopart_Pi, twopart_Ka, twopart_Pr;
    double threepart_Ch, threepart_Pi, threepart_Ka, threepart_Pr;
    double fourpart_Ch, fourpart_Pi, fourpart_Ka, fourpart_Pr;

    for (auto& track : tracks) {
      if (!selTrack(track))
        continue;

      if constexpr (DataFlag) {
        if (std::abs(track.eta()) < 0.8) {
          Nch++;
          pt_ch = track.pt();
          moments(pt_ch, &Q1_ch, &Q2_ch, &Q3_ch, &Q4_ch);

          hist.fill(HIST("QA/after/h_Eta"), track.eta());
          hist.fill(HIST("QA/after/h_Pt"), track.pt());
          hist.fill(HIST("QA/after/h2_Pt_Eta"), track.eta(), track.pt());
          hist.fill(HIST("QA/after/h2_DcaXY"), track.pt(), track.dcaXY());
          hist.fill(HIST("QA/after/h2_DcaZ"), track.pt(), track.dcaZ());

          hist.fill(HIST("QA/after/h_TPCChi2perCluster"), track.tpcChi2NCl());
          hist.fill(HIST("QA/after/h_ITSChi2perCluster"), track.itsChi2NCl());
          hist.fill(HIST("QA/after/h_crossedTPC"), track.tpcNClsCrossedRows());
        }
        hist.fill(HIST("QA/before/h2_TOFSignal"), track.p(), track.beta());
        hist.fill(HIST("QA/before/h2_TPCSignal"), track.p(), track.tpcSignal());

        hist.fill(HIST("QA/Pion/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPi());
        hist.fill(HIST("QA/Pion/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPi());
        hist.fill(HIST("QA/Pion/before/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        hist.fill(HIST("QA/Proton/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaPr());
        hist.fill(HIST("QA/Proton/before/h2_TOFNsigma"), track.p(), track.tofNSigmaPr());
        hist.fill(HIST("QA/Proton/before/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        hist.fill(HIST("QA/Kaon/before/h2_TPCNsigma"), track.p(), track.tpcNSigmaKa());
        hist.fill(HIST("QA/Kaon/before/h2_TOFNsigma"), track.p(), track.tofNSigmaKa());
        hist.fill(HIST("QA/Kaon/before/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());

        // For Pions:
        if (selPions(track)) {
          N_Pi++;
          pt_Pi = track.pt();
          moments(pt_Pi, &Q1_Pi, &Q2_Pi, &Q3_Pi, &Q4_Pi);
          hist.fill(HIST("QA/Pion/h_Pt"), track.pt());
          hist.fill(HIST("QA/Pion/h_Eta"), track.eta());
          hist.fill(HIST("QA/Pion/h_rap"), track.rapidity(massPi));
          hist.fill(HIST("QA/Pion/h2_Pt_rap"), track.rapidity(massPi), track.pt());
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
          hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
          hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
        }

        // For Kaons:
        if (selKaons(track)) {
          N_Ka++;
          pt_Ka = track.pt();
          moments(pt_Ka, &Q1_Ka, &Q2_Ka, &Q3_Ka, &Q4_Ka);
          hist.fill(HIST("QA/Kaon/h_Pt"), track.pt());
          hist.fill(HIST("QA/Kaon/h_Eta"), track.eta());
          hist.fill(HIST("QA/Kaon/h_rap"), track.rapidity(massKa));
          hist.fill(HIST("QA/Kaon/h2_Pt_rap"), track.rapidity(massKa), track.pt());
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
          hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
          hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
        }

        // For Protons:
        if (selProtons(track)) {
          N_Pr++;
          pt_Pr = track.pt();
          moments(pt_Pr, &Q1_Pr, &Q2_Pr, &Q3_Pr, &Q4_Pr);
          hist.fill(HIST("QA/Proton/h_Pt"), track.pt());
          hist.fill(HIST("QA/Proton/h_Eta"), track.eta());
          hist.fill(HIST("QA/Proton/h_rap"), track.rapidity(massPr));
          hist.fill(HIST("QA/Proton/h2_Pt_rap"), track.rapidity(massPr), track.pt());
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
          hist.fill(HIST("QA/after/h2_TPCSignal"), track.p(), track.tpcSignal());
          hist.fill(HIST("QA/after/h2_TOFSignal"), track.p(), track.beta());
        }
      } else if constexpr (RecoFlag) {
        if (track.has_mcParticle() && track.mcParticle().isPhysicalPrimary()) {
          if (std::abs(track.eta()) < 0.8) {
            Nch++;
            pt_ch = track.pt();
            moments(pt_ch, &Q1_ch, &Q2_ch, &Q3_ch, &Q4_ch);
            hist.fill(HIST("QA/Reco/Charged/h_Pt"), track.pt());
            hist.fill(HIST("QA/Reco/Charged/h_Eta"), track.eta());
          }

          if (std::abs(track.mcParticle().pdgCode()) == 211 && selPions(track)) {
            N_Pi++;
            pt_Pi = track.pt();
            moments(pt_Pi, &Q1_Pi, &Q2_Pi, &Q3_Pi, &Q4_Pi);
            hist.fill(HIST("QA/Reco/Pion/h_Pt"), track.pt());
            hist.fill(HIST("QA/Reco/Pion/h_Eta"), track.eta());
            hist.fill(HIST("QA/Reco/Pion/h_y"), track.rapidity(massPi));
            hist.fill(HIST("QA/Reco/Pion/h2_tpcSignal"), track.p(), track.tpcSignal());
            hist.fill(HIST("QA/Reco/Pion/h2_tofSignal"), track.p(), track.beta());
          }

          if (std::abs(track.mcParticle().pdgCode()) == 321 && selKaons(track)) {
            N_Ka++;
            pt_Ka = track.pt();
            moments(pt_Ka, &Q1_Ka, &Q2_Ka, &Q3_Ka, &Q4_Ka);
            hist.fill(HIST("QA/Reco/Kaon/h_Pt"), track.pt());
            hist.fill(HIST("QA/Reco/Kaon/h_Eta"), track.eta());
            hist.fill(HIST("QA/Reco/Kaon/h_y"), track.rapidity(massKa));
            hist.fill(HIST("QA/Reco/Kaon/h2_tpcSignal"), track.p(), track.tpcSignal());
            hist.fill(HIST("QA/Reco/Kaon/h2_tofSignal"), track.p(), track.beta());
          }

          if (std::abs(track.mcParticle().pdgCode()) == 2212 && selProtons(track)) {
            N_Pr++;
            pt_Pr = track.pt();
            moments(pt_Pr, &Q1_Pr, &Q2_Pr, &Q3_Pr, &Q4_Pr);
            hist.fill(HIST("QA/Reco/Proton/h_Pt"), track.pt());
            hist.fill(HIST("QA/Reco/Proton/h_Eta"), track.eta());
            hist.fill(HIST("QA/Reco/Proton/h_y"), track.rapidity(massPr));
            hist.fill(HIST("QA/Reco/Proton/h2_tpcSignal"), track.p(), track.tpcSignal());
            hist.fill(HIST("QA/Reco/Proton/h2_tofSignal"), track.p(), track.beta());
          }
        }
      }
    }

    N_FT0C = col.multFT0C();
    Cent_FT0C = col.centFT0C();
    NTPC = col.multNTracksHasTPC();

    hist.fill(HIST("QA/after/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/after/h_Counts"), 2);
    hist.fill(HIST("QA/after/h_NTPC"), NTPC);
    hist.fill(HIST("QA/after/h_Cent"), Cent_FT0C);
    hist.fill(HIST("QA/after/h_NFT0M"), N_FT0C);
    hist.fill(HIST("QA/after/h2_NTPC_NFT0M"), N_FT0C, NTPC);
    hist.fill(HIST("QA/after/h2_NTPC_Cent"), Cent_FT0C, NTPC);
    hist.fill(HIST("QA/after/p_NTPC_Cent"), Cent_FT0C, NTPC);
    hist.fill(HIST("QA/after/p_NTPC_NFT0M"), N_FT0C, NTPC);
    hist.fill(HIST("QA/after/h2_NTPC_Nch"), NTPC, Nch);

    static constexpr std::string_view dire[] = {"Analysis/Charged/", "Analysis/Pion/", "Analysis/Kaon/", "Analysis/Proton/"};

    hist.fill(HIST(dire[0]) + HIST("h_Mult"), Nch);
    hist.fill(HIST(dire[1]) + HIST("h_Mult"), N_Pi);
    hist.fill(HIST(dire[2]) + HIST("h_Mult"), N_Ka);
    hist.fill(HIST(dire[3]) + HIST("h_Mult"), N_Pr);

    parts(Q1_ch, Q2_ch, Q3_ch, Q4_ch, Nch, &mean_Q1_Ch, &twopart_Ch, &threepart_Ch, &fourpart_Ch);
    if (Nch > 1) {
      if (mean_Q1_Ch != 0) {
        hist.fill(HIST(dire[0]) + HIST("h_mean_Q1"), mean_Q1_Ch);
        hist.fill(HIST(dire[0]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Ch);
        hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Ch, Cent_FT0C);
      }
      if (twopart_Ch != 0) {
        hist.fill(HIST(dire[0]) + HIST("p_twopart_Mult_var"), NTPC, twopart_Ch);
        hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Ch, Cent_FT0C);
      }
    }
    if (Nch > 2) {
      if (mean_Q1_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Ch, Cent_FT0C);
      if (twopart_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Ch, Cent_FT0C);
      if (threepart_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Ch, Cent_FT0C);
    }

    if (Nch > 3) {
      if (mean_Q1_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Ch, Cent_FT0C);
      if (twopart_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Ch, Cent_FT0C);
      if (threepart_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Ch, Cent_FT0C);
      if (fourpart_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Ch, Cent_FT0C);
    }

    parts(Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi, N_Pi, &mean_Q1_Pi, &twopart_Pi, &threepart_Pi, &fourpart_Pi);
    if (N_Pi > 1) {
      if (mean_Q1_Pi != 0) {
        hist.fill(HIST(dire[1]) + HIST("h_mean_Q1"), mean_Q1_Pi);
        hist.fill(HIST(dire[1]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Pi);
        hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Pi, Cent_FT0C);
      }
      if (twopart_Pi != 0) {
        hist.fill(HIST(dire[1]) + HIST("p_twopart_Mult_var"), NTPC, twopart_Pi);
        hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Pi, Cent_FT0C);
      }
    }

    if (N_Pi > 2) {
      if (mean_Q1_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Pi, Cent_FT0C);
      if (twopart_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Pi, Cent_FT0C);
      if (threepart_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Pi, Cent_FT0C);
    }

    if (N_Pi > 3) {
      if (mean_Q1_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Pi, Cent_FT0C);
      if (twopart_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Pi, Cent_FT0C);
      if (threepart_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Pi, Cent_FT0C);
      if (fourpart_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Pi, Cent_FT0C);
    }

    parts(Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka, N_Ka, &mean_Q1_Ka, &twopart_Ka, &threepart_Ka, &fourpart_Ka);
    if (N_Ka > 1) {
      if (mean_Q1_Ka != 0) {
        hist.fill(HIST(dire[2]) + HIST("h_mean_Q1"), mean_Q1_Ka);
        hist.fill(HIST(dire[2]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Ka);
        hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Ka, Cent_FT0C);
      }
      if (twopart_Ka != 0) {
        hist.fill(HIST(dire[2]) + HIST("p_twopart_Mult_var"), NTPC, twopart_Ka);
        hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Ka, Cent_FT0C);
      }
    }
    if (N_Ka > 2) {
      if (mean_Q1_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Ka, Cent_FT0C);
      if (twopart_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Ka, Cent_FT0C);
      if (threepart_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Ka, Cent_FT0C);
    }

    if (N_Ka > 3) {
      if (mean_Q1_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Ka, Cent_FT0C);
      if (twopart_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Ka, Cent_FT0C);
      if (threepart_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Ka, Cent_FT0C);
      if (fourpart_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Ka, Cent_FT0C);
    }

    parts(Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr, N_Pr, &mean_Q1_Pr, &twopart_Pr, &threepart_Pr, &fourpart_Pr);
    if (N_Pr > 1) {
      if (mean_Q1_Pr != 0) {
        hist.fill(HIST(dire[3]) + HIST("h_mean_Q1"), mean_Q1_Pr);
        hist.fill(HIST(dire[3]) + HIST("p_mean_Q1_Mult_var"), NTPC, mean_Q1_Pr);
        hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Pr, Cent_FT0C);
      }
      if (twopart_Pr != 0) {
        hist.fill(HIST(dire[3]) + HIST("p_twopart_Mult_var"), NTPC, twopart_Pr);
        hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Pr, Cent_FT0C);
      }
    }

    if (N_Pr > 2) {
      if (mean_Q1_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_skew"), NTPC, mean_Q1_Pr, Cent_FT0C);
      if (twopart_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_skew"), NTPC, twopart_Pr, Cent_FT0C);
      if (threepart_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_threepart_Mult_skew"), NTPC, threepart_Pr, Cent_FT0C);
    }

    if (N_Pr > 3) {
      if (mean_Q1_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_kurto"), NTPC, mean_Q1_Pr, Cent_FT0C);
      if (twopart_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_kurto"), NTPC, twopart_Pr, Cent_FT0C);
      if (threepart_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_threepart_Mult_kurto"), NTPC, threepart_Pr, Cent_FT0C);
      if (fourpart_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_fourpart_Mult_kurto"), NTPC, fourpart_Pr, Cent_FT0C);
    }
  }

  void process_Run3(MyCollisions::iterator const& col, MyAllTracks const& tracks)
  {
    // Before Collision and Track Cuts:
    for (auto& myTrack : tracks) {
      hist.fill(HIST("QA/before/h_Eta"), myTrack.eta());
      hist.fill(HIST("QA/before/h_Pt"), myTrack.pt());
      hist.fill(HIST("QA/before/h2_Pt_Eta"), myTrack.eta(), myTrack.pt());
      hist.fill(HIST("QA/before/h_TPCChi2perCluster"), myTrack.tpcChi2NCl());
      hist.fill(HIST("QA/before/h_ITSChi2perCluster"), myTrack.itsChi2NCl());
      hist.fill(HIST("QA/before/h_crossedTPC"), myTrack.tpcNClsCrossedRows());
      hist.fill(HIST("QA/before/h2_DcaXY"), myTrack.pt(), myTrack.dcaXY());
      hist.fill(HIST("QA/before/h2_DcaZ"), myTrack.pt(), myTrack.dcaZ());
    }
    hist.fill(HIST("QA/before/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/before/h_Counts"), 2);

    hist.fill(HIST("QA/before/h_NTPC"), col.multTPC());
    hist.fill(HIST("QA/before/h_Cent"), col.centFT0C());
    hist.fill(HIST("QA/before/h_NFT0M"), col.multFT0M());
    hist.fill(HIST("QA/before/h2_NTPC_NFT0M"), col.multFT0M(), col.multTPC());
    hist.fill(HIST("QA/before/h2_NTPC_Cent"), col.centFT0C(), col.multTPC());

    // After Collision and Track Cuts:
    if (selRun3Col(col)) {
      FillHistos<true, false>(col, tracks);
    }
  }
  PROCESS_SWITCH(meanPtFlucId, process_Run3, "Process for Run3", false);

  void process_MCRecoRun3(MyMCCollisions::iterator const& col, aod::McCollisions const&, MyMCTracks const& tracks, aod::McParticles const&)
  {
    if (selRun3Col(col)) {
      FillHistos<false, true>(col, tracks);
    }
  }
  PROCESS_SWITCH(meanPtFlucId, process_MCRecoRun3, "process MC Reconstructed Run-3", true);

  void process_MCGen(soa::Join<aod::McCollisions, aod::MultsExtraMC>::iterator const& mccol, aod::McParticles const& McParticles)
  {
    int N_Pi = 0, N_Ka = 0, N_Pr = 0;
    int Nch = 0, NTPC = 0;
    double pt_ch = 0, Q1_ch = 0, Q2_ch = 0, Q3_ch = 0, Q4_ch = 0;
    double pt_Pi = 0, Q1_Pi = 0, Q2_Pi = 0, Q3_Pi = 0, Q4_Pi = 0;
    double pt_Pr = 0, Q1_Pr = 0, Q2_Pr = 0, Q3_Pr = 0, Q4_Pr = 0;
    double pt_Ka = 0, Q1_Ka = 0, Q2_Ka = 0, Q3_Ka = 0, Q4_Ka = 0;
    double mean_Q1_Ch, mean_Q1_Pi, mean_Q1_Ka, mean_Q1_Pr;
    double twopart_Ch, twopart_Pi, twopart_Ka, twopart_Pr;
    double threepart_Ch, threepart_Pi, threepart_Ka, threepart_Pr;
    double fourpart_Ch, fourpart_Pi, fourpart_Ka, fourpart_Pr;

    if (abs(mccol.posZ()) > posZCut)
      return;

    for (auto& mcParticle : McParticles) {
      if (mcParticle.isPhysicalPrimary() && mcParticle.pt() > ptMin && mcParticle.pt() < ptMax && abs(mcParticle.y()) < rapCut) {
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
    hist.fill(HIST("Gen/Counts"), 2);
    hist.fill(HIST("Gen/vtxZ"), mccol.posZ());
    hist.fill(HIST("Gen/NTPC"), NTPC);

    static constexpr std::string_view dire[] = {"Gen/Charged/", "Gen/Pion/", "Gen/Kaon/", "Gen/Proton/"};

    hist.fill(HIST(dire[0]) + HIST("h_Mult"), Nch);
    hist.fill(HIST(dire[1]) + HIST("h_Mult"), N_Pi);
    hist.fill(HIST(dire[2]) + HIST("h_Mult"), N_Ka);
    hist.fill(HIST(dire[3]) + HIST("h_Mult"), N_Pr);

    parts(Q1_ch, Q2_ch, Q3_ch, Q4_ch, Nch, &mean_Q1_Ch, &twopart_Ch, &threepart_Ch, &fourpart_Ch);
    if (Nch > 1) {
      if (mean_Q1_Ch != 0) {
        hist.fill(HIST(dire[0]) + HIST("h_mean_pt"), mean_Q1_Ch);
        hist.fill(HIST(dire[0]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Ch);
      }
      if (twopart_Ch != 0)
        hist.fill(HIST(dire[0]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Ch);
    }

    parts(Q1_Pi, Q2_Pi, Q3_Pi, Q4_Pi, N_Pi, &mean_Q1_Pi, &twopart_Pi, &threepart_Pi, &fourpart_Pi);
    if (N_Pi > 1) {
      if (mean_Q1_Pi != 0) {
        hist.fill(HIST(dire[1]) + HIST("h_mean_pt"), mean_Q1_Pi);
        hist.fill(HIST(dire[1]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Pi);
      }
      if (twopart_Pi != 0)
        hist.fill(HIST(dire[1]) + HIST("h_twopart_Mult_var"), Nch, twopart_Pi);
    }

    parts(Q1_Ka, Q2_Ka, Q3_Ka, Q4_Ka, N_Ka, &mean_Q1_Ka, &twopart_Ka, &threepart_Ka, &fourpart_Ka);
    if (N_Ka > 1) {
      if (mean_Q1_Ka != 0) {
        hist.fill(HIST(dire[2]) + HIST("h_mean_pt"), mean_Q1_Ka);
        hist.fill(HIST(dire[2]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Ka);
      }
      if (twopart_Ka != 0)
        hist.fill(HIST(dire[2]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Ka);
    }

    parts(Q1_Pr, Q2_Pr, Q3_Pr, Q4_Pr, N_Pr, &mean_Q1_Pr, &twopart_Pr, &threepart_Pr, &fourpart_Pr);
    if (N_Pr > 1) {
      if (mean_Q1_Pr != 0) {
        hist.fill(HIST(dire[3]) + HIST("h_mean_pt"), mean_Q1_Pr);
        hist.fill(HIST(dire[3]) + HIST("h_mean_Q1_Mult_var"), NTPC, mean_Q1_Pr);
      }
      if (twopart_Pr != 0)
        hist.fill(HIST(dire[3]) + HIST("h_twopart_Mult_var"), NTPC, twopart_Pr);
    }
  }
  PROCESS_SWITCH(meanPtFlucId, process_MCGen, "process MC Generated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<meanPtFlucId>(cfgc)};
}
