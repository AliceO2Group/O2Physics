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
/// \brief This task provides the parameters required to calculate the observable
///           v0(pT) along with its statistical uncertainity using subsampling technique.
/// \author Anna Binoy (anna.binoy@niser.ac.in)

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

#include "TDatabasePDG.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
double massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
double massPr = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

struct Diff_pT_fluct_PID {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::vector<std::vector<std::shared_ptr<TProfile2D>>> Subsample;
  TRandom3* fRndm = new TRandom3(0);

  Configurable<int> nPtBins_qa{"nBinsPt_qa", 280, "N bins in pT histo qualitative analysis"};

  Configurable<int> nPtBins{"nBinsPt", 14, "N bins in pT histo"};
  Configurable<int> nEtaBins{"nEtaBins", 100, ""};

  Configurable<float> ptMax{"ptMax", 3.0, "maximum pT"};
  Configurable<float> ptMin{"ptMin", 0.2, "minimum pT"};

  Configurable<float> etaMin{"etaMin", 0.4, "Eta min"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> rapCut{"rapCut", 0.5, "Rapidity Cut"};

  Configurable<float> cfgCutPtUpperTPC{"cfgCutPtUpperTPC", 0.6f, "Upper pT cut for PID using TPC only"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgnSigmaCut{"cfgnSigmaCut", 2.0f, "PID nSigma cut"};

  Configurable<float> cfgnSigmaCut_TPC_pi{"cfgnSigmaCut_TPC_pi", 2.0f, "PID nSigma cut for TPC for pion"};
  Configurable<float> cfgnSigmaCut_TOF_pi{"cfgnSigmaCut_TOF_pi", 3.0f, "PID nSigma cut for TOF for pion"};
  Configurable<float> cfgnSigmaCut_TPC_ka{"cfgnSigmaCut_TPC_ka", 2.0f, "PID nSigma cut for TPC for kaon"};
  Configurable<float> cfgnSigmaCut_TOF_ka{"cfgnSigmaCut_TOF_ka", 3.0f, "PID nSigma cut for TOF for kaon"};
  Configurable<float> cfgnSigmaCut_TPC_pr{"cfgnSigmaCut_TPC_pr", 2.0f, "PID nSigma cut for TPC for proton"};
  Configurable<float> cfgnSigmaCut_TOF_pr{"cfgnSigmaCut_TOF_pr", 3.0f, "PID nSigma cut for TOF for proton"};

  // QualityCuts

  Configurable<float> dcaXYCut{"dcaXYCut", 0.2, "DCAxy cut"};
  Configurable<float> dcaZCut{"dcaZCut", 2.0, "DCAz cut"};
  Configurable<float> posZCut{"posZCut", 10.0, "cut for vertex Z"};

  Configurable<float> TPCNCrossedRowsCut{"TPCNCrossedRowsCut", 2.5, "n_TPC crossed rows Cut"};
  Configurable<float> chi2TPCperClstrCut{"chi2TPCperClstrCut", 4, "Chi2 TPC  per Cluster Cut"};
  Configurable<float> chi2ITSperClstrCut{"chi2ITSperClstrCut", 36, "Chi2 ITS per Cluster Cut"};

  ConfigurableAxis dcaXYBins{"dcaXYBins", {100, -0.15, 0.15}, "dcaXY bins"};
  ConfigurableAxis dcaZBins{"dcaZBins", {100, -1.2, 1.2}, "dcaZ bins"};

  O2_DEFINE_CONFIGURABLE(cfgUse22sEventCut, bool, true, "Use 22s event cut on mult correlations")

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // for the sub-sampling
  Configurable<int> cfgNSubsample{"cfgNSubsample", 10, "Number of subsamples"};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  // This is an example of a convenient declaration of "using"
  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                                aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                                aod::pidTOFbeta>;
  using MyRun2Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>;
  using MyRun3Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::CentFT0Ms>;

  void init(InitContext const&)
  {

    const AxisSpec axisEvents{5, 0, 5, "Counts"};
    const AxisSpec axisEta{nEtaBins, -1., +1., "#eta"};
    const AxisSpec axisY{nEtaBins, -1., +1., "Rapidity"};
    const AxisSpec axisPt{nPtBins, 0.2, 3., "p_{T} (GeV/c)"};
    const AxisSpec axisPt_qa{nPtBins_qa, 0.2, 3., "p_{T} (GeV/c)"};
    const AxisSpec axisP{nPtBins, 0.2, 3., "p (GeV/c)"};
    const AxisSpec axisCent{100, 0., 100, ""};
    const AxisSpec axis1Bin{1, 0., 1, ""};

    const AxisSpec axisNumberOfHadronEtaLess0{3000, 0, 3000, "Number of proton eta less than 0"};
    const AxisSpec axisNumberOfProtonEtaLess0{3000, 0, 3000, "Number of proton eta less than 0"};
    const AxisSpec axisNumberOfPionEtaLess0{3000, 0, 3000, "Number of pion eta less than 0"};
    const AxisSpec axisNumberOfKaonEtaLess0{3000, 0, 3000, "Number of kaon eta less than 0"};

    const AxisSpec axisVtxZ{80, -20., 20., "V_{Z} (cm)"};
    const AxisSpec axisDCAz{dcaZBins, "DCA_{Z} (cm)"};
    const AxisSpec axisDCAxy{dcaXYBins, "DCA_{XY} (cm)"};
    const AxisSpec axisTPCNsigma{500, -5., 5., "n #sigma_{TPC}"};
    const AxisSpec axisTOFNsigma{500, -5., 5., "n #sigma_{TOF}"};
    const AxisSpec axisTPCTOFNsigma{800, -8., 8., "n #sigma_{TOF+TPC}"};
    const AxisSpec axisTPCSignal{720, 20., 200., "#frac{dE}{dx}"};
    const AxisSpec axisTOFSignal{400, 0.2, 1.2, "TOF #beta"};
    const AxisSpec axisChi2{50, 0., 50., "Chi2"};
    const AxisSpec axisCrossedTPC{500, 0, 500, "Crossed TPC"};

    HistogramConfigSpec TOFnSigmaHist({HistType::kTH2D, {axisPt_qa, axisTOFNsigma}});
    HistogramConfigSpec TOFSignalHist({HistType::kTH2D, {axisPt_qa, axisTOFSignal}});
    HistogramConfigSpec TPCnSigmaHist({HistType::kTH2D, {axisPt_qa, axisTPCNsigma}});
    HistogramConfigSpec TPCSignalHist({HistType::kTH2D, {axisPt_qa, axisTPCSignal}});
    HistogramConfigSpec TPCTOFHist({HistType::kTH2D, {axisTPCNsigma, axisTOFNsigma}});
    HistogramConfigSpec TPCTOFnSigmaHist({HistType::kTH2D, {axisPt_qa, axisTPCTOFNsigma}});

    histos.add("fA_hadron", "", kTProfile2D, {axisCent, axisPt});
    histos.add("fA_pion", "", kTProfile2D, {axisCent, axisPt});
    histos.add("fA_kaon", "", kTProfile2D, {axisCent, axisPt});
    histos.add("fA_proton", "", kTProfile2D, {axisCent, axisPt});

    histos.add("fB1_hadron", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fB1_pion", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fB1_kaon", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fB1_proton", "", kTProfile2D, {axisCent, axis1Bin});

    histos.add("fB2_hadron", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fB2_pion", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fB2_kaon", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fB2_proton", "", kTProfile2D, {axisCent, axis1Bin});

    histos.add("fC_hadron", "", kTProfile2D, {axisCent, axisPt});
    histos.add("fC_pion", "", kTProfile2D, {axisCent, axisPt});
    histos.add("fC_kaon", "", kTProfile2D, {axisCent, axisPt});
    histos.add("fC_proton", "", kTProfile2D, {axisCent, axisPt});

    histos.add("fD_hadron", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fD_pion", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fD_kaon", "", kTProfile2D, {axisCent, axis1Bin});
    histos.add("fD_proton", "", kTProfile2D, {axisCent, axis1Bin});

    // QA Plots:
    histos.add("QA/before/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});
    histos.add("QA/before/h_TPCChi2perCluster", "TPC #Chi^{2}/Cluster", kTH1D, {axisChi2});
    histos.add("QA/before/h_ITSChi2perCluster", "ITS #Chi^{2}/Cluster", kTH1D, {axisChi2});
    histos.add("QA/before/h_crossedTPC", "Crossed TPC", kTH1D, {axisCrossedTPC});
    histos.add("QA/before/h2_DcaZ", "DCA_{Z}", kTH1D, {axisDCAz});
    histos.add("QA/before/h2_DcaXY", "DCA_{XY}", kTH1D, {axisDCAxy});
    histos.add("QA/before/h2_TPCSignal", "TPC Signal", TPCSignalHist);
    histos.add("QA/before/h2_TOFSignal", "TOF Signal", TOFSignalHist);

    histos.addClone("QA/before/", "QA/after/");

    histos.add("QA/Pion/h_Pt", "p_{T} (TPC & TPC+TOF)", kTH1D, {axisPt_qa});
    histos.add("QA/Pion/h_rap", "y (TPC & TPC+TOF)", kTH1D, {axisY});
    histos.add("QA/Pion/h2_Pt_rap", "p_{T} vs y", kTH2D, {{axisY}, {axisPt_qa}});
    histos.add("QA/Pion/h2_DcaZ", "DCA_{z}", kTH2D, {{axisPt_qa}, {axisDCAz}});
    histos.add("QA/Pion/h2_DcaXY", "DCA_{xy}", kTH2D, {{axisPt_qa}, {axisDCAxy}});

    histos.add("QA/Pion/before/h2_TPCNsigma", "n #sigma_{TPC}", TPCnSigmaHist);
    histos.add("QA/Pion/before/h2_TOFNsigma", "n #sigma_{TOF}", TOFnSigmaHist);
    histos.add("QA/Pion/before/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", TPCTOFHist);
    histos.add("QA/Pion/before/h2_TpcTofNsigma1", "n #sigma_{TPC+TOF}", TPCTOFnSigmaHist);

    histos.add("QA/Pion/h2_TPCNsigma", "n #sigma_{TPC}", TPCnSigmaHist);
    histos.add("QA/Pion/h2_TOFNsigma", "n #sigma_{TOF}", TOFnSigmaHist);
    histos.add("QA/Pion/h2_TpcTofNsigma", "n #sigma_{TPC} vs n #sigma_{TOF}", TPCTOFHist);
    histos.add("QA/Pion/h2_TpcTofNsigma1", "n #sigma_{TPC+TOF}", TPCTOFnSigmaHist);

    histos.add("QA/Pion/h2_TPCSignal", "TPC Signal vs pT", TPCSignalHist);
    histos.add("QA/Pion/h2_TOFSignal", "TOF Signal vs pT", TOFSignalHist);
    histos.add("QA/Pion/h2_ExpTPCSignal", "Expected TPC Signal vs pT", TPCSignalHist);

    histos.addClone("QA/Pion/", "QA/Kaon/");
    histos.addClone("QA/Pion/", "QA/Proton/");

    // Define Subsamples
    Subsample.resize(cfgNSubsample);
    for (int i = 0; i < cfgNSubsample; i++) {
      Subsample[i].resize(20);
    }
    for (int i = 0; i < cfgNSubsample; i++) {
      Subsample[i][0] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fA_hadron", i), "", {HistType::kTProfile2D, {axisCent, axisPt}}));
      Subsample[i][1] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fA_pion", i), "", {HistType::kTProfile2D, {axisCent, axisPt}}));
      Subsample[i][2] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fA_kaon", i), "", {HistType::kTProfile2D, {axisCent, axisPt}}));
      Subsample[i][3] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fA_proton", i), "", {HistType::kTProfile2D, {axisCent, axisPt}}));

      Subsample[i][4] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fB1_hadron", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][5] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fB1_pion", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][6] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fB1_kaon", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][7] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fB1_proton", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));

      Subsample[i][8] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fB2_hadron", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][9] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fB2_pion", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][10] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fB2_kaon", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][11] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fB2_proton", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));

      Subsample[i][12] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fC_hadron", i), "", {HistType::kTProfile2D, {axisCent, axisPt}}));
      Subsample[i][13] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fC_pion", i), "", {HistType::kTProfile2D, {axisCent, axisPt}}));
      Subsample[i][14] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fC_kaon", i), "", {HistType::kTProfile2D, {axisCent, axisPt}}));
      Subsample[i][15] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fC_proton", i), "", {HistType::kTProfile2D, {axisCent, axisPt}}));

      Subsample[i][16] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fD_hadron", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][17] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fD_pion", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][18] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fD_kaon", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
      Subsample[i][19] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/fD_proton", i), "", {HistType::kTProfile2D, {axisCent, axis1Bin}}));
    }
  }

  template <typename T>
  bool selRun2Col(T const& col)
  {
    if (std::abs(col.posZ()) > posZCut)
      return false;

    if (!col.sel7())
      return false;

    if (!col.alias_bit(kINT7))
      return false;

    return true;
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

    // kinematic cuts

    if (track.pt() < ptMin)
      return false;

    if (track.pt() > ptMax)
      return false;

    if (std::abs(track.eta()) > etaCut)
      return false;

    if (std::abs(track.dcaZ()) > dcaZCut)
      return false;

    if (std::abs(track.dcaXY()) > dcaXYCut)
      return false;

    if (track.tpcChi2NCl() > chi2TPCperClstrCut)
      return false;

    if (track.itsChi2NCl() > chi2ITSperClstrCut)
      return false;

    if (track.tpcNClsCrossedRows() > TPCNCrossedRowsCut)
      return false;

    if (!track.isGlobalTrack())
      return false;

    return true;
  }

  template <typename T>
  bool selPions(T const& track)
  {
    const float combNSigmaPi = std::sqrt(pow(track.tpcNSigmaPi(), 2.0) + pow(track.tofNSigmaPi(), 2.0));
    const float combNSigmaKa = std::sqrt(pow(track.tpcNSigmaKa(), 2.0) + pow(track.tofNSigmaKa(), 2.0));
    const float combNSigmaPr = std::sqrt(pow(track.tpcNSigmaPr(), 2.0) + pow(track.tofNSigmaPr(), 2.0));

    if (track.pt() <= cfgCutPtUpperTPC) {
      Int_t flag = 0;
      if (track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cfgnSigmaCut_TPC_pi)
        flag += 1;
      if (track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cfgnSigmaCut_TPC_pi)
        flag += 1;
      if (track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cfgnSigmaCut_TPC_pi)
        flag += 1;
      if (flag > 1)
        return false;
      else if (flag == 1 && std::abs(track.tpcNSigmaPi()) < cfgnSigmaCut_TPC_pi)
        return true;

    } else if (track.pt() > cfgCutPtUpperTPC) {
      Int_t flag = 0;
      if (track.hasTOF() && track.hasTPC() && combNSigmaPi < cfgnSigmaCut_TOF_pi)
        flag += 1;
      if (track.hasTOF() && track.hasTPC() && combNSigmaKa < cfgnSigmaCut_TOF_pi)
        flag += 1;
      if (track.hasTOF() && track.hasTPC() && combNSigmaPr < cfgnSigmaCut_TOF_pi)
        flag += 1;
      if (flag > 1)
        return false;
      else if (flag == 1 && combNSigmaPi < cfgnSigmaCut_TOF_pi)
        return true;
    }

    return false;
  }

  template <typename T>
  bool selKaons(T const& track)
  {
    const float combNSigmaPi = std::sqrt(pow(track.tpcNSigmaPi(), 2.0) + pow(track.tofNSigmaPi(), 2.0));
    const float combNSigmaKa = std::sqrt(pow(track.tpcNSigmaKa(), 2.0) + pow(track.tofNSigmaKa(), 2.0));
    const float combNSigmaPr = std::sqrt(pow(track.tpcNSigmaPr(), 2.0) + pow(track.tofNSigmaPr(), 2.0));

    if (track.pt() <= cfgCutPtUpperTPC) {
      Int_t flag = 0;
      if (track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cfgnSigmaCut_TPC_ka)
        flag += 1;
      if (track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cfgnSigmaCut_TPC_ka)
        flag += 1;
      if (track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cfgnSigmaCut_TPC_ka)
        flag += 1;
      if (flag > 1)
        return false;
      else if (flag == 1 && std::abs(track.tpcNSigmaKa()) < cfgnSigmaCut_TPC_ka)
        return true;

    } else if (track.pt() > cfgCutPtUpperTPC) {
      Int_t flag = 0;
      if (track.hasTOF() && track.hasTPC() && combNSigmaPi < cfgnSigmaCut_TOF_ka)
        flag += 1;
      if (track.hasTOF() && track.hasTPC() && combNSigmaKa < cfgnSigmaCut_TOF_ka)
        flag += 1;
      if (track.hasTOF() && track.hasTPC() && combNSigmaPr < cfgnSigmaCut_TOF_ka)
        flag += 1;
      if (flag > 1)
        return false;
      else if (flag == 1 && combNSigmaKa < cfgnSigmaCut_TOF_ka)
        return true;
    }

    return false;
  }

  template <typename T>
  bool selProtons(T const& track)
  {
    const float combNSigmaPi = std::sqrt(pow(track.tpcNSigmaPi(), 2.0) + pow(track.tofNSigmaPi(), 2.0));
    const float combNSigmaKa = std::sqrt(pow(track.tpcNSigmaKa(), 2.0) + pow(track.tofNSigmaKa(), 2.0));
    const float combNSigmaPr = std::sqrt(pow(track.tpcNSigmaPr(), 2.0) + pow(track.tofNSigmaPr(), 2.0));

    // if (abs(track.rapidity(massPr)) < 0.5)
    //   return true;
    if (track.pt() <= cfgCutPtUpperTPC && track.pt() > 0.4) {
      Int_t flag = 0;
      if (track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cfgnSigmaCut_TPC_pr)
        flag += 1;
      if (track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cfgnSigmaCut_TPC_pr)
        flag += 1;
      if (track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cfgnSigmaCut_TPC_pr)
        flag += 1;
      if (flag > 1)
        return false;
      else if (flag == 1 && std::abs(track.tpcNSigmaPr()) < cfgnSigmaCut_TPC_pr)
        return true;

    } else if (track.pt() > cfgCutPtUpperTPC) {
      Int_t flag = 0;
      if (track.hasTOF() && track.hasTPC() && combNSigmaPi < cfgnSigmaCut_TOF_pr)
        flag += 1;
      if (track.hasTOF() && track.hasTPC() && combNSigmaKa < cfgnSigmaCut_TOF_pr)
        flag += 1;
      if (track.hasTOF() && track.hasTPC() && combNSigmaPr < cfgnSigmaCut_TOF_pr)
        flag += 1;
      if (flag > 1)
        return false;
      else if (flag == 1 && combNSigmaPr < cfgnSigmaCut_TOF_pr)
        return true;
    }

    return false;
  }

  void process(MyRun3Collisions::iterator const& col, MyAllTracks const& tracks)
  {
    double Cent_FT0M = 0;

    if (selRun3Col(col)) {

      Cent_FT0M = col.centFT0M();

      Double_t pT_bin[14] = {0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9};

      Double_t N_Pi_eta_less_0 = 0;
      Double_t N_Ka_eta_less_0 = 0;
      Double_t N_Pr_eta_less_0 = 0;
      Double_t Nch_eta_less_0 = 0;

      Double_t pT_sum_etaLess0 = 0;
      Double_t pT_sum_etaGreaterEtamin = 0;
      Double_t N_sum_etaGreaterEtamin = 0;

      Double_t pt_Ch = 0, pt_Pi = 0, pt_Ka = 0, pt_Pr = 0;

      Double_t fA_hadron[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Double_t fA_pion[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Double_t fA_kaon[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Double_t fA_proton[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      Double_t fB1_hadron = 0;
      Double_t fB1_pion = 0;
      Double_t fB1_kaon = 0;
      Double_t fB1_proton = 0;

      Double_t fB2_hadron = 0;
      Double_t fB2_pion = 0;
      Double_t fB2_kaon = 0;
      Double_t fB2_proton = 0;

      Double_t fC_hadron[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Double_t fC_pion[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Double_t fC_kaon[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Double_t fC_proton[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      Double_t fD_hadron = 0;
      Double_t fD_pion = 0;
      Double_t fD_kaon = 0;
      Double_t fD_proton = 0;

      // normal creation of a histogram
      TH1D* fPt_profile = new TH1D("fPt_profile", "fPt_profile", 14, 0.2, 3);
      TH1D* fPt_profile_pion = new TH1D("fPt_profile_pion", "fPt_profile_pion", 14, 0.2, 3);
      TH1D* fPt_profile_kaon = new TH1D("fPt_profile_kaon", "fPt_profile_kaon", 14, 0.2, 3);
      TH1D* fPt_profile_proton = new TH1D("fPt_profile_proton", "fPt_profile_proton", 14, 0.2, 3);

      for (auto& track : tracks) {

        histos.fill(HIST("QA/before/h2_DcaXY"), track.dcaXY());
        histos.fill(HIST("QA/before/h2_DcaZ"), track.dcaZ());
        histos.fill(HIST("QA/before/h_VtxZ"), col.posZ());

        histos.fill(HIST("QA/before/h_TPCChi2perCluster"), track.tpcChi2NCl());
        histos.fill(HIST("QA/before/h_ITSChi2perCluster"), track.itsChi2NCl());
        histos.fill(HIST("QA/before/h_crossedTPC"), track.tpcNClsCrossedRows());

        if (selTrack(track))
          continue;

        histos.fill(HIST("QA/after/h2_DcaXY"), track.dcaXY());
        histos.fill(HIST("QA/after/h2_DcaZ"), track.dcaZ());

        histos.fill(HIST("QA/after/h_TPCChi2perCluster"), track.tpcChi2NCl());
        histos.fill(HIST("QA/after/h_ITSChi2perCluster"), track.itsChi2NCl());
        histos.fill(HIST("QA/after/h_crossedTPC"), track.tpcNClsCrossedRows());

        histos.fill(HIST("QA/before/h2_TOFSignal"), track.pt(), track.beta());
        histos.fill(HIST("QA/before/h2_TPCSignal"), track.pt(), track.tpcSignal());

        histos.fill(HIST("QA/Pion/before/h2_TPCNsigma"), track.pt(), track.tpcNSigmaPi());
        histos.fill(HIST("QA/Pion/before/h2_TOFNsigma"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("QA/Pion/before/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());

        histos.fill(HIST("QA/Proton/before/h2_TPCNsigma"), track.pt(), track.tpcNSigmaPr());
        histos.fill(HIST("QA/Proton/before/h2_TOFNsigma"), track.pt(), track.tofNSigmaPr());
        histos.fill(HIST("QA/Proton/before/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());

        histos.fill(HIST("QA/Kaon/before/h2_TPCNsigma"), track.pt(), track.tpcNSigmaKa());
        histos.fill(HIST("QA/Kaon/before/h2_TOFNsigma"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QA/Kaon/before/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());

        const float combNSigmaPi1 = std::sqrt(pow(track.tpcNSigmaPi(), 2.0) + pow(track.tofNSigmaPi(), 2.0));
        const float combNSigmaKa1 = std::sqrt(pow(track.tpcNSigmaKa(), 2.0) + pow(track.tofNSigmaKa(), 2.0));
        const float combNSigmaPr1 = std::sqrt(pow(track.tpcNSigmaPr(), 2.0) + pow(track.tofNSigmaPr(), 2.0));

        histos.fill(HIST("QA/Pion/before/h2_TpcTofNsigma1"), track.pt(), combNSigmaPi1);
        histos.fill(HIST("QA/Kaon/before/h2_TpcTofNsigma1"), track.pt(), combNSigmaKa1);
        histos.fill(HIST("QA/Proton/before/h2_TpcTofNsigma1"), track.pt(), combNSigmaPr1);

        if (track.eta() < 0) {
          Nch_eta_less_0++;
          pt_Ch = track.pt();

          pT_sum_etaLess0 += pt_Ch;
          fPt_profile->Fill(pt_Ch);

          // For Pions:
          if (selPions(track)) {
            N_Pi_eta_less_0++;
            pt_Pi = track.pt();

            fPt_profile_pion->Fill(pt_Pi);

            // QA
            histos.fill(HIST("QA/Pion/h_Pt"), track.pt());
            histos.fill(HIST("QA/Pion/h_rap"), track.rapidity(massPi));
            histos.fill(HIST("QA/Pion/h2_Pt_rap"), track.rapidity(massPi), track.pt());
            histos.fill(HIST("QA/Pion/h2_DcaXY"), track.pt(), track.dcaXY());
            histos.fill(HIST("QA/Pion/h2_DcaZ"), track.pt(), track.dcaZ());

            histos.fill(HIST("QA/Pion/h2_TPCNsigma"), track.pt(), track.tpcNSigmaPi());
            histos.fill(HIST("QA/Pion/h2_TOFNsigma"), track.pt(), track.tofNSigmaPi());
            histos.fill(HIST("QA/Pion/h2_TpcTofNsigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
            histos.fill(HIST("QA/Pion/h2_TOFSignal"), track.pt(), track.beta());
            histos.fill(HIST("QA/Pion/h2_TPCSignal"), track.pt(), track.tpcSignal());
            histos.fill(HIST("QA/Pion/h2_ExpTPCSignal"), track.pt(), track.tpcExpSignalPi(track.tpcSignal()));

            histos.fill(HIST("QA/after/h2_TOFSignal"), track.pt(), track.beta());
            histos.fill(HIST("QA/after/h2_TPCSignal"), track.pt(), track.tpcSignal());

            const float combNSigmaPi2 = std::sqrt(pow(track.tpcNSigmaPi(), 2.0) + pow(track.tofNSigmaPi(), 2.0));

            histos.fill(HIST("QA/Pion/h2_TpcTofNsigma1"), track.pt(), combNSigmaPi2);
          }

          // For Kaons:
          if (selKaons(track)) {
            N_Ka_eta_less_0++;
            pt_Ka = track.pt();

            fPt_profile_kaon->Fill(pt_Ka);

            // QA

            histos.fill(HIST("QA/Kaon/h_Pt"), track.pt());
            histos.fill(HIST("QA/Kaon/h_rap"), track.rapidity(massKa));
            histos.fill(HIST("QA/Kaon/h2_Pt_rap"), track.rapidity(massKa), track.pt());
            histos.fill(HIST("QA/Kaon/h2_DcaXY"), track.pt(), track.dcaXY());
            histos.fill(HIST("QA/Kaon/h2_DcaZ"), track.pt(), track.dcaZ());

            histos.fill(HIST("QA/Kaon/h2_TPCNsigma"), track.pt(), track.tpcNSigmaKa());
            histos.fill(HIST("QA/Kaon/h2_TOFNsigma"), track.pt(), track.tofNSigmaKa());
            histos.fill(HIST("QA/Kaon/h2_TpcTofNsigma"), track.tpcNSigmaKa(), track.tofNSigmaKa());
            histos.fill(HIST("QA/Kaon/h2_TOFSignal"), track.pt(), track.beta());
            histos.fill(HIST("QA/Kaon/h2_TPCSignal"), track.pt(), track.tpcSignal());
            histos.fill(HIST("QA/Kaon/h2_ExpTPCSignal"), track.pt(), track.tpcExpSignalKa(track.tpcSignal()));
            histos.fill(HIST("QA/after/h2_TOFSignal"), track.pt(), track.beta());
            histos.fill(HIST("QA/after/h2_TPCSignal"), track.pt(), track.tpcSignal());

            const float combNSigmaKa2 = std::sqrt(pow(track.tpcNSigmaKa(), 2.0) + pow(track.tofNSigmaKa(), 2.0));

            histos.fill(HIST("QA/Kaon/h2_TpcTofNsigma1"), track.pt(), combNSigmaKa2);
          }

          // For Protons:
          if (selProtons(track)) {
            N_Pr_eta_less_0++;
            pt_Pr = track.pt();

            fPt_profile_proton->Fill(pt_Pr);

            // QA

            histos.fill(HIST("QA/Proton/h_Pt"), track.pt());
            histos.fill(HIST("QA/Proton/h_rap"), track.rapidity(massPr));
            histos.fill(HIST("QA/Proton/h2_Pt_rap"), track.rapidity(massPr), track.pt());
            histos.fill(HIST("QA/Proton/h2_DcaZ"), track.pt(), track.dcaZ());
            histos.fill(HIST("QA/Proton/h2_DcaXY"), track.pt(), track.dcaXY());

            histos.fill(HIST("QA/Proton/h2_TPCNsigma"), track.pt(), track.tpcNSigmaPr());
            histos.fill(HIST("QA/Proton/h2_TOFNsigma"), track.pt(), track.tofNSigmaPr());
            histos.fill(HIST("QA/Proton/h2_TpcTofNsigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
            histos.fill(HIST("QA/Proton/h2_TPCSignal"), track.pt(), track.tpcSignal());
            histos.fill(HIST("QA/Proton/h2_TOFSignal"), track.pt(), track.beta());
            histos.fill(HIST("QA/Proton/h2_ExpTPCSignal"), track.pt(), track.tpcExpSignalPr(track.tpcSignal()));
            histos.fill(HIST("QA/after/h2_TPCSignal"), track.pt(), track.tpcSignal());
            histos.fill(HIST("QA/after/h2_TOFSignal"), track.pt(), track.beta());

            const float combNSigmaPr2 = std::sqrt(pow(track.tpcNSigmaPr(), 2.0) + pow(track.tofNSigmaPr(), 2.0));

            histos.fill(HIST("QA/Proton/h2_TpcTofNsigma1"), track.pt(), combNSigmaPr2);
          }

        } else if (track.eta() > etaMin) {
          pT_sum_etaGreaterEtamin += pt_Ch;
          N_sum_etaGreaterEtamin++;
        }
      }

      // selecting subsample and filling profiles
      float l_Random = fRndm->Rndm();
      int SampleIndex = static_cast<int>(cfgNSubsample * l_Random);

      // B1, B2, and D Calculation for hadrons
      if (N_sum_etaGreaterEtamin != 0 && Nch_eta_less_0 != 0) {
        fB1_hadron = pT_sum_etaLess0 / Nch_eta_less_0;
        histos.fill(HIST("fB1_hadron"), Cent_FT0M, 0.5, fB1_hadron);
        Subsample[SampleIndex][4]->Fill(Cent_FT0M, 0.5, fB1_hadron);

        fB2_hadron = pT_sum_etaGreaterEtamin / N_sum_etaGreaterEtamin;
        histos.fill(HIST("fB2_hadron"), Cent_FT0M, 0.5, fB2_hadron);
        Subsample[SampleIndex][8]->Fill(Cent_FT0M, 0.5, fB2_hadron);

        fD_hadron = fB1_hadron * fB2_hadron;
        histos.fill(HIST("fD_hadron"), Cent_FT0M, 0.5, fD_hadron);
        Subsample[SampleIndex][16]->Fill(Cent_FT0M, 0.5, fD_hadron);
      }

      // B1, B2, and D Calculation for pions
      if (N_sum_etaGreaterEtamin != 0 && Nch_eta_less_0 != 0 && N_Pi_eta_less_0 != 0) {
        fB1_pion = pT_sum_etaLess0 / Nch_eta_less_0;
        histos.fill(HIST("fB1_pion"), Cent_FT0M, 0.5, fB1_pion);
        Subsample[SampleIndex][5]->Fill(Cent_FT0M, 0.5, fB1_pion);

        fB2_pion = pT_sum_etaGreaterEtamin / N_sum_etaGreaterEtamin;
        histos.fill(HIST("fB2_pion"), Cent_FT0M, 0.5, fB2_pion);
        Subsample[SampleIndex][9]->Fill(Cent_FT0M, 0.5, fB2_pion);

        fD_pion = fB1_pion * fB2_pion;
        histos.fill(HIST("fD_pion"), Cent_FT0M, 0.5, fD_pion);
        Subsample[SampleIndex][17]->Fill(Cent_FT0M, 0.5, fD_pion);
      }

      // B1, B2, and D Calculation for kaons
      if (N_sum_etaGreaterEtamin != 0 && Nch_eta_less_0 != 0 && N_Ka_eta_less_0 != 0) {
        fB1_kaon = pT_sum_etaLess0 / Nch_eta_less_0;
        histos.fill(HIST("fB1_kaon"), Cent_FT0M, 0.5, fB1_kaon);
        Subsample[SampleIndex][6]->Fill(Cent_FT0M, 0.5, fB1_kaon);

        fB2_kaon = pT_sum_etaGreaterEtamin / N_sum_etaGreaterEtamin;
        histos.fill(HIST("fB2_kaon"), Cent_FT0M, 0.5, fB2_kaon);
        Subsample[SampleIndex][10]->Fill(Cent_FT0M, 0.5, fB2_kaon);

        fD_kaon = fB1_kaon * fB2_kaon;
        histos.fill(HIST("fD_kaon"), Cent_FT0M, 0.5, fD_kaon);
        Subsample[SampleIndex][18]->Fill(Cent_FT0M, 0.5, fD_kaon);
      }

      // B1, B2, and D Calculation for protons
      if (N_sum_etaGreaterEtamin != 0 && Nch_eta_less_0 != 0 && N_Pr_eta_less_0 != 0) {
        fB1_proton = pT_sum_etaLess0 / Nch_eta_less_0;
        histos.fill(HIST("fB1_proton"), Cent_FT0M, 0.5, fB1_proton);
        Subsample[SampleIndex][7]->Fill(Cent_FT0M, 0.5, fB1_proton);

        fB2_proton = pT_sum_etaGreaterEtamin / N_sum_etaGreaterEtamin;
        histos.fill(HIST("fB2_proton"), Cent_FT0M, 0.5, fB2_proton);
        Subsample[SampleIndex][11]->Fill(Cent_FT0M, 0.5, fB2_proton);

        fD_proton = fB1_proton * fB2_proton;
        histos.fill(HIST("fD_proton"), Cent_FT0M, 0.5, fD_proton);
        Subsample[SampleIndex][19]->Fill(Cent_FT0M, 0.5, fD_proton);
      }

      for (int i = 0; i < 14; i++) {
        // A_hadrone Calculation
        if (Nch_eta_less_0 != 0) {
          fA_hadron[i] = fPt_profile->GetBinContent(i + 1) / Nch_eta_less_0;
          histos.fill(HIST("fA_hadron"), Cent_FT0M, pT_bin[i], fA_hadron[i]);
          Subsample[SampleIndex][0]->Fill(Cent_FT0M, pT_bin[i], fA_hadron[i]);
        }

        // A_pion Calculation
        if (N_Pi_eta_less_0 != 0) {
          fA_pion[i] = fPt_profile_pion->GetBinContent(i + 1) / N_Pi_eta_less_0;
          histos.fill(HIST("fA_pion"), Cent_FT0M, pT_bin[i], fA_pion[i]);
          Subsample[SampleIndex][1]->Fill(Cent_FT0M, pT_bin[i], fA_pion[i]);
        }

        // A_kaon Calculation
        if (N_Ka_eta_less_0 != 0) {
          fA_kaon[i] = fPt_profile_kaon->GetBinContent(i + 1) / N_Ka_eta_less_0;
          histos.fill(HIST("fA_kaon"), Cent_FT0M, pT_bin[i], fA_kaon[i]);
          Subsample[SampleIndex][2]->Fill(Cent_FT0M, pT_bin[i], fA_kaon[i]);
        }

        // A_proton Calculation
        if (N_Pr_eta_less_0 != 0) {
          fA_proton[i] = fPt_profile_proton->GetBinContent(i + 1) / N_Pr_eta_less_0;
          histos.fill(HIST("fA_proton"), Cent_FT0M, pT_bin[i], fA_proton[i]);
          Subsample[SampleIndex][3]->Fill(Cent_FT0M, pT_bin[i], fA_proton[i]);
        }

        // C_hadron Calculation
        if (Nch_eta_less_0 != 0 && N_sum_etaGreaterEtamin != 0) {
          fC_hadron[i] = fA_hadron[i] * fB2_hadron;
          histos.fill(HIST("fC_hadron"), Cent_FT0M, pT_bin[i], fC_hadron[i]);
          Subsample[SampleIndex][12]->Fill(Cent_FT0M, pT_bin[i], fC_hadron[i]);
        }

        // C_pion Calculation
        if (N_Pi_eta_less_0 != 0 && N_sum_etaGreaterEtamin != 0) {
          fC_pion[i] = fA_pion[i] * fB2_pion;
          histos.fill(HIST("fC_pion"), Cent_FT0M, pT_bin[i], fC_pion[i]);
          Subsample[SampleIndex][13]->Fill(Cent_FT0M, pT_bin[i], fC_pion[i]);
        }

        // A_kaon Calculation
        if (N_Ka_eta_less_0 != 0 && N_sum_etaGreaterEtamin != 0) {
          fC_kaon[i] = fA_kaon[i] * fB2_kaon;
          histos.fill(HIST("fC_kaon"), Cent_FT0M, pT_bin[i], fC_kaon[i]);
          Subsample[SampleIndex][14]->Fill(Cent_FT0M, pT_bin[i], fC_kaon[i]);
        }

        // A_proton Calculation
        if (N_Pr_eta_less_0 != 0 && N_sum_etaGreaterEtamin != 0) {
          fC_proton[i] = fA_proton[i] * fB2_proton;
          histos.fill(HIST("fC_proton"), Cent_FT0M, pT_bin[i], fC_proton[i]);
          Subsample[SampleIndex][15]->Fill(Cent_FT0M, pT_bin[i], fC_proton[i]);
        }
      }

      fPt_profile->Delete();
      fPt_profile_pion->Delete();
      fPt_profile_kaon->Delete();
      fPt_profile_proton->Delete();
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Diff_pT_fluct_PID>(cfgc)};
}
