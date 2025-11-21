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
/// \brief This task is for Normalized Factorial Moments Analysis: PRC 85, 044914 (2012), nucl-ex:1411.6083
/// \author Salman Malik
/// \author Balwan Singh

#include <iostream>
#include <array>
#include <TH1F.h>
#include "TRandom.h"
// O2 includes
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "Framework/ASoAHelpers.h"
#include <TPDGCode.h>
#include <unordered_set> 
#include "TDatabasePDG.h"
using std::array;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

TH1D* tmpFqErr[6][5][52];

struct FactorialMoments {

  Configurable<Float_t> confEta{"centralEta", 0.9, "eta limit for tracks"};
  Configurable<Int_t> confNumPt{"numPt", 5, "number of pT bins"};
  Configurable<Float_t> confPtMin{"ptMin", 0.2f, "lower pT cut"};
  Configurable<Float_t> confDCAxy{"dcaXY", 2.4f, "DCA xy cut"};
  Configurable<Float_t> confDCAz{"dcaZ", 2.0f, "DCA z cut"};
  Configurable<Float_t> confMinTPCCls{"minTPCCls", 70.0f, "minimum number of TPC clusters"};
  Configurable<std::vector<Int_t>> confCentCut{"centLimits", {0, 5}, "centrality min and max"};
  Configurable<std::vector<Float_t>> confVertex{"vertexXYZ", {0.3f, 0.4f, 10.0f}, "vertex cuts"};
  Configurable<std::vector<Float_t>> confPtBins{"ptCuts", {0.2f, 2.0f}, "pT cuts"};
  Configurable<bool> IsApplySameBunchPileup{"IsApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> IsApplyGoodZvtxFT0vsPV{"IsApplyGoodZvtxFT0vsPV", true, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> IsApplyVertexITSTPC{"IsApplyVertexITSTPC", true, "Enable VertexITSTPC cut"};
  Configurable<bool> IsApplyVertexTOFmatched{"IsApplyVertexTOFmatched", true, "Enable VertexTOFmatched cut"};
  Configurable<bool> IsApplyVertexTRDmatched{"IsApplyVertexTRDmatched", true, "Enable VertexTRDmatched cut"};
  Configurable<bool> IsApplyExtraCorrCut{"IsApplyExtraCorrCut", false, "Enable extra NPVtracks vs FTOC correlation cut"};
  Configurable<bool> IsApplyExtraPhiCut{"IsApplyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<bool> includeGlobalTracks{"includeGlobalTracks", false, "Enable Global Tracks"};
  Configurable<bool> includeTPCTracks{"includeTPCTracks", false, "TPC Tracks"};
  Configurable<bool> includeITSTracks{"includeITSTracks", false, "ITS Tracks"};
  Configurable<int> confSamplesize{"samplesize", 100, "Sample size"};
  Configurable<bool> confUseMC{"useMC", false, "Use MC information"};
    Configurable<int> reduceOutput{"reduce-output", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};
  Filter filterTracks = (nabs(aod::track::eta) < confEta) && (aod::track::pt >= confPtMin) && (nabs(aod::track::dcaXY) < confDCAxy) && (nabs(aod::track::dcaZ) < confDCAz);
  Filter filterCollisions = (nabs(aod::collision::posZ) < confVertex.value[2]) && (nabs(aod::collision::posX) < confVertex.value[0]) && (nabs(aod::collision::posY) < confVertex.value[1]);
  Service<o2::framework::O2DatabasePDG> pdg;

  // Histograms
  HistogramRegistry histos1{
    "histos1",
    {
     {"hRecoPtBefore", "Reco pT before cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 20.0}}}},
    {"hGenPtBefore", "Gen pT before cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 20.0}}}},
    {"hRecoPtAfter", "Reco pT after cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 20.0}}}},
    {"hGenPtAfter", "Gen pT after cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 20.0}}}},
    {"hRecoEtaBefore", "Reco #eta before cuts;#eta;Counts", {HistType::kTH1F, {{200, -2.0, 2.0}}}},
    {"mCentMCImpactFull", "Impact Parameter (MC All Events);Impact Parameter b [fm];Counts", {HistType::kTH1F, {{200, 0.0, 50.0}}}},
    {"mCentMCImpact005", "Impact Parameter (MC 0–5%);Impact Parameter b [fm];Counts", {HistType::kTH1F, {{200, 0.0, 20.0}}}},
    {"hGenEtaBefore", "Gen #eta before cuts;#eta;Counts", {HistType::kTH1F, {{200, -2.0, 2.0}}}},
    {"hRecoEtaAfter", "Reco #eta after cuts;#eta;Counts", {HistType::kTH1F, {{200, -2.0, 2.0}}}},
    {"hGenEtaAfter", "Gen #eta after cuts;#eta;Counts", {HistType::kTH1F, {{200, -2.0, 2.0}}}},
    {"hRecoPhiBefore", "Reco #phi before cuts;#phi;Counts", {HistType::kTH1F, {{100, 0, 6.3}}}},
    {"hGenPhiBefore", "Gen #phi before cuts;#phi;Counts", {HistType::kTH1F, {{100, 0, 6.3}}}},
    {"hRecoPhiAfter", "Reco #phi after cuts;#phi;Counts", {HistType::kTH1F, {{100, 0, 6.3}}}},
    {"hGenCharge", "Gen particle charge;Charge;Counts", {HistType::kTH1F, {{13, -6.5, 6.5}}}},
    {"hGenPhiAfter", "Gen #phi after cuts;#phi;Counts", {HistType::kTH1F, {{100, 0, 6.3}}}},
      {"mChargeDist", "MC particle charge;Charge;Counts", {HistType::kTH1F, {{9, -4.5, 4.5}}}},
      {"mDiffPhi", "Δ#phi between selected MC particles;Δ#phi;Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
          {"mPrimariesPerEvent", "Primary MC particles per event;N_{primaries};Counts", {HistType::kTH1I, {{20000, 0, 20000}}}},
      {"mDiffPt", "Δp_{T} between selected MC particles;Δp_{T} (GeV/c);Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mDiffEta", "Δ#eta between selected MC particles;Δ#eta;Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mADiffPhi", "Δ#phi between selected MC particles;Δ#phi;Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mADiffPt", "Δp_{T} between selected MC particles;Δp_{T} (GeV/c);Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mADiffEta", "Δ#eta between selected MC particles;Δ#eta;Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
    },
  };

  HistogramRegistry histos{
    "histos",
    {
      {"mCollID", "collisionID", {HistType::kTH1I, {{1000, -10000, 10000}}}},
      {"mCentFT0M", "centFT0C_0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mCentFV0A", "centFV0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mCentFT0A", "centFT0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mCentFT0C", "centFT0C", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mVertexX", "vertexX", {HistType::kTH1F, {{100, -20, 20}}}},
      {"mVertexY", "vertexY", {HistType::kTH1F, {{100, -20, 20}}}},
      {"mVertexZ", "vertexZ", {HistType::kTH1F, {{100, -20, 20}}}},
      {"mEta", "#eta", {HistType::kTH1F, {{1000, -2, 2}}}},
      {"mPt", "#pt", {HistType::kTH1F, {{1000, -0.01, 50}}}},
      {"mPhi", "#phi", {HistType::kTH1F, {{100, 0, TMath::TwoPi()}}}},
      {"mEvents", "events", {HistType::kTH1D, {{5, -0.5, 4.5}}}},
      {"mNFindableClsTPC", "findable TPC clusters;findable clusters", {HistType::kTH1F, {{100, 0, 200}}}},
      {"mNClsTPC", "number of clusters TPC; nClusters TPC", {HistType::kTH1F, {{100, 0, 200}}}},
      {"mNClsITS", "number of clusters ITS; nClusters ITS", {HistType::kTH1F, {{100, 0, 10}}}},
      {"mChi2TPC", "chi2 TPC", {HistType::kTH1F, {{100, 0, 10}}}},
      {"mChi2ITS", "chi2 ITS", {HistType::kTH1F, {{100, 0, 10}}}},
      {"mChi2TRD", "chi2 TRD", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mDCAxy", "DCA xy", {HistType::kTH1F, {{500, -0.8, 0.8}}}},
      {"mDCAx", "DCA z", {HistType::kTH1F, {{500, -2.0, 2.0}}}},
      {"mDCAxyPt", "DCA xy vs #pt;#pt;DCAxy", {HistType::kTH2F, {{100, 0, 20}, {500, -0.5, 0.5}}}},
      {"mDCAxyPtbcut", "DCA xy vs #pt;#pt;DCAxycut", {HistType::kTH2F, {{100, 0, 20}, {500, -0.5, 0.5}}}},
      {"mDCAzPtbcut", "DCA z vs #pt;#pt;DCAzcut", {HistType::kTH2F, {{100, 0, 20}, {100, -2.0, 2.0}}}},
      {"mDCAzPt", "DCA z vs #pt;#pt;DCAz", {HistType::kTH2F, {{100, 0, 20}, {100, -2.0, 2.0}}}},
      {"mNSharedClsTPC", "shared clusters in TPC", {HistType::kTH1F, {{100, 0, 10}}}},
      {"mCrossedRowsTPC", "crossedrows in TPC", {HistType::kTH1F, {{100, 0, 200}}}},
      {"mNFinClsminusCRows", "findable cluster #minus crossed rows (TPC)", {HistType::kTH1F, {{100, 0, 200}}}},
      {"mNFractionShClsTPC", "fraction of shared clusters in TPC", {HistType::kTH1F, {{100, 0, 2}}}},
      {"mSharedClsvsPt", "shared cluster vs #pt", {HistType::kTH2F, {{100, 0, 50}, {100, 0, 10}}}},
      {"mSharedClsProbvsPt", "shared clusters ration vs #pt;#pt;sharedcls/ncrows", {HistType::kTH2F, {{100, 0, 50}, {100, 0, 5}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    true};
  static const Int_t nBins = 52;
  Int_t countSamples = 0;
  Int_t testc1 = 0, testc2 = 0, testc3 = 0;
  array<Int_t, nBins> binningM;
  array<Int_t, 5> countTracks{0, 0, 0, 0, 0};
  array<array<array<Double_t, nBins>, 5>, 6> fqEvent;
  array<array<array<Double_t, nBins>, 5>, 6> fqEventSampled;
  array<array<Double_t, nBins>, 5> binConEvent;
  array<array<array<Double_t, nBins>, 5>, 6> binConEventSampled;
  array<array<array<Double_t, nBins>, 5>, 6> errorFq = {{{{{0, 0, 0, 0, 0}}}}};
  std::vector<std::shared_ptr<TH2>> mHistArrReset;
  std::vector<std::shared_ptr<TH1>> mHistArrQA;
  std::vector<std::shared_ptr<TH1>> mFqBinFinal;
  std::vector<std::shared_ptr<TH1>> mBinConFinal;
  std::vector<std::shared_ptr<TH1>> mFqBinFinalSampled;
  std::vector<std::shared_ptr<TH1>> mBinConFinalSampled;
  std::vector<std::shared_ptr<TH1>> mFqError;

  // max number of bins restricted to 5
  static constexpr array<std::string_view, 5>
    mbinNames{"bin1/", "bin2/", "bin3/", "bin4/", "bin5/"};
  void init(o2::framework::InitContext&)
  {
    // NOTE: check to make number of pt and the vector consistent
    if (confNumPt != static_cast<int>(confPtBins.value.size()) / 2) {
      for (int i = confNumPt; i < static_cast<int>(confPtBins.value.size() / 2); i++) {
        confPtBins.value[2 * i] = 0;
        confPtBins.value[2 * i + 1] = 0;
      }
    }
    AxisSpec axisPt[5] = {{100, -0.01, 3 * confPtBins.value[1], ""}, {100, -0.01, 3 * confPtBins.value[3], ""}, {100, -0.01, 3 * confPtBins.value[5], ""}, {100, -0.01, 3 * confPtBins.value[7], ""}, {100, -0.01, 3 * confPtBins.value[9], ""}}; // pT axis
    auto mEventSelected = std::get<std::shared_ptr<TH1>>(histos.add("mEventSelected", "eventSelected", HistType::kTH1D, {{8, 0.5, 8.5}}));
    mEventSelected->GetXaxis()->SetBinLabel(1, "all");
    mEventSelected->GetXaxis()->SetBinLabel(2, "sel8");
    mEventSelected->GetXaxis()->SetBinLabel(3, "sameBunchPileup");
    mEventSelected->GetXaxis()->SetBinLabel(4, "goodZvtxFT0vsPV");
    mEventSelected->GetXaxis()->SetBinLabel(5, "vertexITSTPC");
    mEventSelected->GetXaxis()->SetBinLabel(6, "centrality");
    mEventSelected->GetXaxis()->SetBinLabel(7, "final");
    auto mTrackSelected = std::get<std::shared_ptr<TH1>>(histos.add(
      "mTrackSelected", "Track Selection Steps", HistType::kTH1D, {{5, 0.5, 5.5}}));

    mTrackSelected->GetXaxis()->SetBinLabel(1, "all");
    mTrackSelected->GetXaxis()->SetBinLabel(2, "charge");
    mTrackSelected->GetXaxis()->SetBinLabel(3, "PIDs");
    mTrackSelected->GetXaxis()->SetBinLabel(4, "Pids exclude pions");
    mTrackSelected->GetXaxis()->SetBinLabel(5, "Pids exclude pions + kions");
    if (confUseMC) {
      auto mMcTrackSelected = std::get<std::shared_ptr<TH1>>(histos.add("mMcTrackSelected", "mcTrackSelection", HistType::kTH1D, {{5, 0.5, 5.5}}));
    }

    for (auto iM = 0; iM < nBins; ++iM) {
      binningM[iM] = 2 * (iM + 2);
    }

    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mEta", iPt + 1), Form("#eta for bin %.2f-%.2f;#eta", confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{1000, -2, 2}})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mPt", iPt + 1), Form("pT for bin %.2f-%.2f;pT", confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {axisPt[iPt]})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mPhi", iPt + 1), Form("#phi for bin %.2f-%.2f;#phi", confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{1000, 0, 2 * TMath::Pi()}})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mMultiplicity", iPt + 1), Form("Multiplicity for bin %.2f-%.2f;Multiplicity", confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{1000, 0, 15000}})));
      for (auto iM = 0; iM < nBins; ++iM) {
        auto mHistsR = std::get<std::shared_ptr<TH2>>(histos.add(Form("bin%i/Reset/mEtaPhi%i", iPt + 1, iM), Form("#eta#phi_%i for bin %.2f-%.2f;#eta;#phi", iM, confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH2F, {{binningM[iM], -0.8, 0.8}, {binningM[iM], 0, 2 * TMath::Pi()}}));
        mHistArrReset.push_back(mHistsR);

        for (auto iq = 0; iq < 6; ++iq) {
          tmpFqErr[iq][iPt][iM] = new TH1D(Form("tmpFqErr%i%i%i", iq, iPt, iM), Form("tmpFqErr%i%i%i", iq, iPt, iM), 100, 0, 10);
        }
      }
      for (auto i = 0; i < 6; ++i) {
        auto mHistFq = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalFq%i_bin%i", i + 2, iPt + 1), Form("Final F_%i for bin %.2f-%.2f;M", i + 2, confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mFqBinFinal.push_back(mHistFq);
        auto mHistAv = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalAvBin%i_bin%i", i + 2, iPt + 1), Form("Final AvBin_%i for bin %.2f-%.2f;M", i + 2, confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mBinConFinal.push_back(mHistAv);

        auto mHistFqSampled = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalFq%iSampled_bin%i", i + 2, iPt + 1), Form("Final F_%i for bin %.2f-%.2f;M", i + 2, confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mFqBinFinalSampled.push_back(mHistFqSampled);
        auto mHistAvSampled = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalAvBin%iSampled_bin%i", i + 2, iPt + 1), Form("Final AvBin_%i for bin %.2f-%.2f;M", i + 2, confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mBinConFinalSampled.push_back(mHistAvSampled);

        auto mHistError = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFqError%i_bin%i", i + 2, iPt + 1), Form("Error for F_%i for bin %.2f-%.2f;M", i + 2, confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mFqError.push_back(mHistError);
      }
    }

    array<std::string, 3> charge{"pos", "neg", "all"};
    array<std::string, 2> detadphiRange{"", "close"};
    array<std::string, 4> ptOrd{"pt2>pt1_dphi<0", "pt2>pt1_dphi>0", "pt2<pt1_dphi<0", "pt2<pt1_dphi>0"};
    for (Int_t i = 0; i < confNumPt; ++i) {
      mHistArrQA.push_back(std::get<std::shared_ptr<TH2>>(histos.add(Form("mDetaDphiBef%i", i), Form("dEta dPhi for ptbin %d;#Delta#eta;#Delta#phi", i), HistType::kTH2F, {{35, -1.75, 1.75}, {50, -TMath::Pi(), TMath::Pi()}})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH2>>(histos.add(Form("mDetaDphiAf%i", i), Form("dEta dPhi for ptbin %d;#Delta#eta;#Delta#phi", i), HistType::kTH2F, {{50, -1.75, 1.75}, {50, -TMath::Pi(), TMath::Pi()}})));

      for (Int_t ch = 0; ch < 3; ++ch) {
        std::string histName = Form("mDetaDphiCh_%.2f-%.2f_%s%s_%d", confPtBins.value[2 * i], confPtBins.value[2 * i + 1], charge[ch].c_str(), detadphiRange[0].c_str(), i);
        std::replace(histName.begin(), histName.end(), '.', '_');

        mHistArrQA.push_back(std::get<std::shared_ptr<TH2>>(histos.add(histName.c_str(), Form("dEta dPhi for ptbin %d %s;#Delta#eta;#Delta#phi", i, charge[ch].data()), HistType::kTH2F, {{35, -1.75, 1.75}, {50, -TMath::Pi(), TMath::Pi()}})));
        std::string histName2 = Form("mDetaDphiCh_%.2f-%.2f_%s%s_%d", confPtBins.value[2 * i], confPtBins.value[2 * i + 1], charge[ch].data(), detadphiRange[1].data(), i);
        std::replace(histName2.begin(), histName2.end(), '.', '_');

        mHistArrQA.push_back(std::get<std::shared_ptr<TH2>>(histos.add(histName2.c_str(), Form("dEta dPhi for ptbin %d %s %s;#Delta#eta;#Delta#phi", i, charge[ch].data(), detadphiRange[1].data()), HistType::kTH2F, {{50, -0.02, 0.02}, {50, -0.05, 0.05}})));

        for (Int_t j = 0; j < 4; ++j) {
          std::string ptOrdSafe = ptOrd[j].data();
          std::replace(ptOrdSafe.begin(), ptOrdSafe.end(), '>', '_');
          std::replace(ptOrdSafe.begin(), ptOrdSafe.end(), '<', '_');
          std::replace(ptOrdSafe.begin(), ptOrdSafe.end(), '.', '_');
          std::string histNamePtOrd = Form("mDetaDphiPtOrd_%.2f-%.2f_%s_%s%s_%d", confPtBins.value[2 * i], confPtBins.value[2 * i + 1], ptOrdSafe.c_str(), charge[ch].data(), detadphiRange[1].data(), i);
          std::replace(histNamePtOrd.begin(), histNamePtOrd.end(), '.', '_');
          histNamePtOrd += std::to_string(j); // Append unique index to avoid collisions
          mHistArrQA.push_back(std::get<std::shared_ptr<TH2>>(histos1.add(histNamePtOrd.c_str(), Form("dEta dPhi for ptbin %d %s %s %s;#Delta#eta;#Delta#phi", i, ptOrdSafe.c_str(), charge[ch].data(), detadphiRange[1].data()), HistType::kTH2F, {{50, -0.02, 0.02}, {50, -0.05, 0.05}})));

          std::replace(ptOrdSafe.begin(), ptOrdSafe.end(), '>', '_');
          std::replace(ptOrdSafe.begin(), ptOrdSafe.end(), '<', '_');
          std::replace(ptOrdSafe.begin(), ptOrdSafe.end(), '.', '_');

          std::replace(histNamePtOrd.begin(), histNamePtOrd.end(), '.', '_');
          histNamePtOrd += std::to_string(j); // Append unique index to avoid collisions
          mHistArrQA.push_back(std::get<std::shared_ptr<TH2>>(histos1.add(histNamePtOrd.c_str(), Form("dEta dPhi for ptbin %d %s %s %s;#Delta#eta;#Delta#phi", i, ptOrdSafe.c_str(), charge[ch].data(), detadphiRange[1].data()), HistType::kTH2F, {{50, -0.02, 0.02}, {50, -0.05, 0.05}})));
        }
      }
      std::string histNameChDiff = Form("mDetaDphiChDiff_%.2f-%.2f%s", confPtBins.value[2 * i], confPtBins.value[2 * i + 1], detadphiRange[0].data());
      std::replace(histNameChDiff.begin(), histNameChDiff.end(), '.', '_');

      histNameChDiff += "_" + std::to_string(i); // Append unique index to avoid collisions

      mHistArrQA.push_back(std::get<std::shared_ptr<TH2>>(histos1.add(histNameChDiff.c_str(), Form("dEta dPhi for ptbin %d different ch;#Delta#eta;#Delta#phi", i), HistType::kTH2F, {{35, -1.75, 1.75}, {50, -TMath::Pi(), TMath::Pi()}})));

      std::string histNameChDiffClose = Form("mDetaDphiChDiffClose_%.2f-%.2f%s", confPtBins.value[2 * i], confPtBins.value[2 * i + 1], detadphiRange[1].data());
      std::replace(histNameChDiffClose.begin(), histNameChDiffClose.end(), '.', '_');

      histNameChDiffClose += "_" + std::to_string(i); // Append unique index to avoid collisions

      mHistArrQA.push_back(std::get<std::shared_ptr<TH2>>(histos1.add(histNameChDiffClose.c_str(), Form("dEta dPhi for ptbin %d different ch close;#Delta#eta;#Delta#phi", i), HistType::kTH2F, {{50, -0.02, 0.02}, {50, -0.05, 0.05}})));
    }
  }
  int chargeFromPDG(int pdg)
{
  if (pdg == 0) return 0;
  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  TParticlePDG* particle = pdgDB->GetParticle(pdg);
  return particle ? int(particle->Charge() / 3.0) : 0;
}
template <typename T>
void checkpT(const T& track)
{
  //int gIndex = track.globalIndex();
  //if (usedIndices.find(gIndex) != usedIndices.end()) {
   // return; // Already used this track
  //}
  //usedIndices.insert(gIndex); // Mark this track as used

  for (auto iPt = 0; iPt < confNumPt; ++iPt) {
    if (track.pt() > confPtBins.value[2 * iPt] && track.pt() < confPtBins.value[2 * iPt + 1]) {
      float iphi = track.phi();
      iphi = gRandom->Gaus(iphi, TMath::TwoPi());

      if (iphi < 0) {
        iphi += TMath::TwoPi();
      } else if (iphi > TMath::TwoPi()) {
        iphi -= TMath::TwoPi();
      }

      mHistArrQA[iPt * 4]->Fill(track.eta());
      mHistArrQA[iPt * 4 + 1]->Fill(track.pt());
      mHistArrQA[iPt * 4 + 2]->Fill(track.phi());
      countTracks[iPt]++;

      for (auto iM = 0; iM < nBins; ++iM) {
        mHistArrReset[iPt * nBins + iM]->Fill(track.eta(), track.phi());
      }
    }
  }
}

  void calculateMoments(std::vector<std::shared_ptr<TH2>> hist)
  {
    Double_t binContent = 0;
    countSamples++;
    Bool_t compSample = kFALSE;
    if (countSamples == confSamplesize) {
      compSample = kTRUE;
      countSamples = 0;
    }
    // Calculate the normalized factorial moments
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      for (auto iM = 0; iM < nBins; ++iM) {
        binContent = 0;
        Double_t sumfqBin[6] = {0};

        for (auto iEta = 1; iEta <= hist[iPt * nBins + iM]->GetNbinsX(); ++iEta) {
          for (auto iPhi = 1; iPhi <= hist[iPt * nBins + iM]->GetNbinsY(); ++iPhi) {
            double binconVal = 0;
            binconVal = hist[iPt * nBins + iM]->GetBinContent(iEta, iPhi);
            binContent += binconVal;
            for (auto iq = 0; iq < 6; ++iq) {
              Double_t fqBin = 0;
              if (binconVal >= iq + 2) {
                fqBin = TMath::Factorial(binconVal) / (TMath::Factorial(binconVal - (iq + 2)));
              }
              if (std::isnan(fqBin)) {
                break;
              }
              sumfqBin[iq] += fqBin;
            }
          }
        }
        binConEvent[iPt][iM] = binContent / (TMath::Power(binningM[iM], 2));
        for (auto iq = 0; iq < 6; ++iq) {
          if (sumfqBin[iq] > 0) {
            fqEvent[iq][iPt][iM] = sumfqBin[iq] / (TMath::Power(binningM[iM], 2));
            fqEventSampled[iq][iPt][iM] += fqEvent[iq][iPt][iM];
          }
          binConEventSampled[iq][iPt][iM] += binConEvent[iPt][iM];
          mFqBinFinal[iPt * 6 + iq]->Fill(iM, fqEvent[iq][iPt][iM]);
          mBinConFinal[iPt * 6 + iq]->Fill(iM, binConEvent[iPt][iM]);
          if (compSample) {
            mBinConFinalSampled[iPt * 6 + iq]->Fill(iM, binConEventSampled[iq][iPt][iM] / confSamplesize);

            double tmp = (fqEventSampled[iq][iPt][iM] / (confSamplesize)) / (pow(binConEventSampled[iq][iPt][iM] / (confSamplesize), (iq + 2)));
            mFqBinFinalSampled[iPt * 6 + iq]->Fill(iM, tmp);
            tmpFqErr[iq][iPt][iM]->Fill(tmp);
            errorFq[iq][iPt][iM] += TMath::Power(fqEventSampled[iq][iPt][iM] / (confSamplesize), 2);
            mFqError[iPt * 6 + iq]->SetBinContent(iM + 1, 0);
            mFqError[iPt * 6 + iq]->Fill(iM, tmpFqErr[iq][iPt][iM]->GetStdDev());

            fqEventSampled[iq][iPt][iM] = 0;
            binConEventSampled[iq][iPt][iM] = 0;
          }
        } // end of loop over order of F
      } // end of loop over M bins
    } // end of loop over pT bins
  }

  // write a template that takes coll and tracks from processRun3 and processMCRec
  using TracksFMs = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>::iterator const& coll, TracksFMs const& tracks)
  {
    // selection of events
    if (!coll.sel8()) {
      return;
    }
    if (IsApplySameBunchPileup && !coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (IsApplyGoodZvtxFT0vsPV && !coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (IsApplyVertexITSTPC && !coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    if (coll.centFT0C() < confCentCut.value[0] || coll.centFT0C() > confCentCut.value[1]) {
      return;
    }
    histos.fill(HIST("mVertexX"), coll.posX());
    histos.fill(HIST("mVertexY"), coll.posY());
    histos.fill(HIST("mVertexZ"), coll.posZ());
    histos.fill(HIST("mCentFT0M"), coll.centFT0M());
    histos.fill(HIST("mCentFV0A"), coll.centFV0A());
    histos.fill(HIST("mCentFT0A"), coll.centFT0A());
    histos.fill(HIST("mCentFT0C"), coll.centFT0C());
    for (auto const& h : mHistArrReset) {
      h->Reset();
    }
    countTracks = {0, 0, 0, 0, 0};
    fqEvent = {{{{{0, 0, 0, 0, 0, 0}}}}};
    binConEvent = {{{0, 0, 0, 0, 0}}};
    for (auto const& track : tracks) {
      if (track.hasTPC())
      // if (track.hasITS())
      // if (track.isGlobalTrack())
      {
        histos.fill(HIST("mCollID"), track.collisionId());
        histos.fill(HIST("mEta"), track.eta());
        histos.fill(HIST("mPt"), track.pt());
        histos.fill(HIST("mPhi"), track.phi());
        histos.fill(HIST("mNFindableClsTPC"), track.tpcNClsFindable());
        histos.fill(HIST("mNClsTPC"), track.tpcNClsFound());
        histos.fill(HIST("mNClsITS"), track.itsNCls());
        histos.fill(HIST("mChi2TPC"), track.tpcChi2NCl());
        histos.fill(HIST("mChi2ITS"), track.itsChi2NCl());
        histos.fill(HIST("mChi2TRD"), track.trdChi2());
        histos.fill(HIST("mDCAxy"), track.dcaXY());
        histos.fill(HIST("mDCAx"), track.dcaZ());
        histos.fill(HIST("mDCAxyPt"), track.pt(), track.dcaXY());
        histos.fill(HIST("mDCAzPt"), track.pt(), track.dcaZ());
        histos.fill(HIST("mNSharedClsTPC"), track.tpcNClsShared());
        histos.fill(HIST("mCrossedRowsTPC"), track.tpcNClsCrossedRows());
        histos.fill(HIST("mNFinClsminusCRows"), track.tpcNClsFindableMinusCrossedRows());
        histos.fill(HIST("mNFractionShClsTPC"), track.tpcFractionSharedCls());
        histos.fill(HIST("mSharedClsvsPt"), track.pt(), track.tpcNClsShared());
        histos.fill(HIST("mSharedClsProbvsPt"), track.pt(), track.tpcFractionSharedCls() / track.tpcNClsCrossedRows());
        checkpT(track);
      }
    }
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      // if (countTracks[iPt] > 0)
      if (countTracks[iPt] > 0) {
        mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
      }
    }
    // Calculate the normalized factorial moments
    calculateMoments(mHistArrReset);
  }
  PROCESS_SWITCH(FactorialMoments, processRun3, "main process function", false);

  using CollisionCandidateMCRec = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  using TracksMc = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>>;
  void processMCRec(soa::Filtered<CollisionCandidateMCRec>::iterator const& coll, TracksMc const& colltracks, aod::McParticles const&, aod::McCollisions const&)
  {
    if (!coll.has_mcCollision()) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 0);
    if (!coll.sel8()) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 1);
    if (IsApplySameBunchPileup && !coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 2);
    if (IsApplyGoodZvtxFT0vsPV && !coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 3);
    if (IsApplyVertexITSTPC && !coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 4);
    if (coll.centFT0C() < confCentCut.value[0] || coll.centFT0C() > confCentCut.value[1]) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 5);
    histos.fill(HIST("mVertexX"), coll.posX());
    histos.fill(HIST("mVertexY"), coll.posY());
    histos.fill(HIST("mVertexZ"), coll.posZ());
    // histos.fill(HIST("mCentFT0M"), coll.centFT0M());
    // histos.fill(HIST("mCentFV0A"), coll.centFV0A());
    // histos.fill(HIST("mCentFT0A"), coll.centFT0A());
    histos.fill(HIST("mCentFT0C"), coll.centFT0C());
    for (auto const& h : mHistArrReset) {
      h->Reset();
    }
    countTracks = {0, 0, 0, 0, 0};
    fqEvent = {{{{{0, 0, 0, 0, 0, 0}}}}};
    binConEvent = {{{0, 0, 0, 0, 0}}};
    for (auto const& track : colltracks) {
      // if (track.hasITS())
      // if (track.hasTPC())
      // if (track.isGlobalTrack()) {
      histos.fill(HIST("mCollID"), track.collisionId());
      histos.fill(HIST("mEta"), track.eta());
      histos.fill(HIST("mPt"), track.pt());
      histos.fill(HIST("mPhi"), track.phi());
      histos.fill(HIST("mNFindableClsTPC"), track.tpcNClsFindable());
      histos.fill(HIST("mNClsTPC"), track.tpcNClsFound());
      histos.fill(HIST("mNClsITS"), track.itsNCls());
      histos.fill(HIST("mChi2TPC"), track.tpcChi2NCl());
      histos.fill(HIST("mChi2ITS"), track.itsChi2NCl());
      histos.fill(HIST("mChi2TRD"), track.trdChi2());
      histos.fill(HIST("mDCAxy"), track.dcaXY());
      histos.fill(HIST("mDCAx"), track.dcaZ());
      histos.fill(HIST("mDCAxyPt"), track.pt(), track.dcaXY());
      histos.fill(HIST("mDCAzPt"), track.pt(), track.dcaZ());
      histos.fill(HIST("mNSharedClsTPC"), track.tpcNClsShared());
      histos.fill(HIST("mCrossedRowsTPC"), track.tpcNClsCrossedRows());
      histos.fill(HIST("mNFinClsminusCRows"), track.tpcNClsFindableMinusCrossedRows());
      histos.fill(HIST("mNFractionShClsTPC"), track.tpcFractionSharedCls());
      histos.fill(HIST("mSharedClsvsPt"), track.pt(), track.tpcNClsShared());
      histos.fill(HIST("mSharedClsProbvsPt"), track.pt(), track.tpcFractionSharedCls() / track.tpcNClsCrossedRows());
      checkpT(track);
      //}
    }
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      // if (countTracks[iPt] > 0)countTracks = {0, 0, 0, 0, 0};
      if (countTracks[iPt] > 0) {
        mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
      }
    }
    histos.fill(HIST("mEventSelected"), 6);
    // Calculate the normalized factorial moments
    calculateMoments(mHistArrReset);
  }
  PROCESS_SWITCH(FactorialMoments, processMCRec, "main process function", false);
   
 void processMCgenr3(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
{
   
    // access MC truth information with mcCollision() and mcParticle() methods
     
    if (reduceOutput < 2) {
      LOGF(info, "MC. vtx-z = %f", mcCollision.posZ());
      LOGF(info, "First: %d | Length: %d", mcParticles.begin().index(), mcParticles.size());
    }
    int count = 0;
    countTracks = {0, 0, 0, 0, 0};
    fqEvent = {{{{{0, 0, 0, 0, 0, 0}}}}};
    binConEvent = {{{0, 0, 0, 0, 0}}};
    if (mcCollision.impactParameter() >= 3.5) {
    return;
  }
    for (auto& mcParticle : mcParticles) {
    int pdgCode = mcParticle.pdgCode();
auto* pd = TDatabasePDG::Instance()->GetParticle(pdgCode);
double charge = (pd != nullptr) ? pd->Charge() : 0.;
      if (mcParticle.isPhysicalPrimary() && std::abs(mcParticle.eta()) < 0.8 && std::abs(charge) >= 1e-6) {
         histos.fill(HIST("mEta"), mcParticle.eta());
      histos.fill(HIST("mPt"), mcParticle.pt());
      histos.fill(HIST("mPhi"), mcParticle.phi());
        count++;
        // Loop over mothers and daughters
        if (mcParticle.has_mothers()) {
          // Check first mother
          auto const& mother = mcParticle.mothers_first_as<aod::McParticles>();
          if (reduceOutput == 0) {
            LOGF(info, "First mother: %d has pdg code %d", mother.globalIndex(), mother.pdgCode());
          }
          // Loop over all mothers (needed for some MCs with junctions etc.)
          for (auto& m : mcParticle.mothers_as<aod::McParticles>()) {
            LOGF(debug, "M2 %d %d", mcParticle.globalIndex(), m.globalIndex());
          }
        }
        if (mcParticle.has_daughters()) {
          for (auto& d : mcParticle.daughters_as<aod::McParticles>()) {
            LOGF(debug, "D2 %d %d", mcParticle.globalIndex(), d.globalIndex());
          }
        }
            checkpT(mcParticle);
            }
                histos1.fill(HIST("mPrimariesPerEvent"), count); 
    }
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      // if (countTracks[iPt] > 0)
      if (countTracks[iPt] > 0) {
        mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
      }
    }
    if (reduceOutput < 2) {
      LOGF(info, "Primaries for this collision: %d", count);
    }
    calculateMoments(mHistArrReset);
  }
PROCESS_SWITCH(FactorialMoments, processMCgenr3, "main process function", false);
using EventSelection_Run2 = soa::Join<aod::EvSels, aod::Mults, aod::CentRun2V0Ms, aod::CentRun2SPDTrks>;
  using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  using CollisionRecSim_Run2 = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection_Run2>>::iterator;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;  
Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
 
void processMcRun2(CollisionRecSim_Run2 const& coll,
                   aod::BCs const& bcs,
                   TracksRecSim const& tracks,
                   aod::McParticles const& mcParticles,
                   aod::McCollisions const& mcCollisions,
                   BCsWithRun2Info const& bcsInfo)
{
  auto bc = coll.bc_as<BCsWithRun2Info>();
  if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted))) {
    return;  // Skip event if it doesn't pass cuts
  }

  if (coll.centRun2V0M() < confCentCut.value[0] || coll.centRun2V0M() > confCentCut.value[1]) {
    return;
  }

  /// Fill basic vertex and centrality histograms
  histos.fill(HIST("mVertexX"), coll.posX());
  histos.fill(HIST("mVertexY"), coll.posY());
  histos.fill(HIST("mVertexZ"), coll.posZ());
  histos.fill(HIST("mCentFT0M"), coll.centRun2V0M());

  /// Reset histograms and counters
  for (auto const& h : mHistArrReset) {
    h->Reset();
  }

  countTracks = {0, 0, 0, 0, 0};
  fqEvent = {{{{{0, 0, 0, 0, 0, 0}}}}};
  binConEvent = {{{0, 0, 0, 0, 0}}};

  /// Loop over reconstructed tracks
  for (auto const& track : tracks) {
    double recoCharge = (track.sign() != 0) ? track.sign() : 0.;
    if (std::abs(track.eta()) < 0.8 && track.isGlobalTrack() && std::abs(recoCharge) >= 1e-6) {
      histos.fill(HIST("mCollID"), track.collisionId());
      
      histos.fill(HIST("mNFindableClsTPC"), track.tpcNClsFindable());
      histos.fill(HIST("mNClsTPC"), track.tpcNClsFound());
      histos.fill(HIST("mNClsITS"), track.itsNCls());
      histos.fill(HIST("mChi2TPC"), track.tpcChi2NCl());
      histos.fill(HIST("mChi2ITS"), track.itsChi2NCl());
      histos.fill(HIST("mChi2TRD"), track.trdChi2());
      histos.fill(HIST("mNSharedClsTPC"), track.tpcNClsShared());
      histos.fill(HIST("mCrossedRowsTPC"), track.tpcNClsCrossedRows());
      histos.fill(HIST("mNFinClsminusCRows"), track.tpcNClsFindableMinusCrossedRows());
      histos.fill(HIST("mNFractionShClsTPC"), track.tpcFractionSharedCls());
      histos.fill(HIST("mSharedClsvsPt"), track.pt(), track.tpcNClsShared());
      histos.fill(HIST("mSharedClsProbvsPt"), track.pt(), track.tpcFractionSharedCls() / track.tpcNClsCrossedRows());
    
   checkpT(track);
    }
  }

  /// Access associated mcCollision from index
  auto mccollision = coll.mcCollision_as<aod::McCollisions>();

  /// Impact parameter can be used here if needed
  float impactpar = mccollision.impactParameter();
  // Optionally: fillMCCollision<false>(coll, mcParticles, impactpar);

  /// Now slice over selectedMCParticles using mcCollision
  auto mcParts = mcParticles.sliceBy(perMcCollision, coll.mcCollision().globalIndex());

for (auto const& mc : mcParts) {
  int pdgCode = mc.pdgCode();
auto* pd = TDatabasePDG::Instance()->GetParticle(pdgCode);
double charge = (pd != nullptr) ? pd->Charge() : 0.;
      if (mc.isPhysicalPrimary() && std::abs(mc.eta()) < 0.8 && std::abs(charge) >= 1e-6) {
      histos.fill(HIST("mEta"), mc.eta());
      histos.fill(HIST("mPt"), mc.pt());
      histos.fill(HIST("mPhi"), mc.phi());

  // Use charged MC particle here
  //checkpT(mc);
}
}

  /// Apply cuts on tracks per pT bin
  for (auto iPt = 0; iPt < confNumPt; ++iPt) {
    if (countTracks[iPt] > 0) {
      mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
    } else {
      return;
    }
  }

  /// Final calculation of normalized factorial moments
  calculateMoments(mHistArrReset);
}

PROCESS_SWITCH(FactorialMoments, processMcRun2, "process MC Run2", true);

void processMcGenRun2(aod::McCollision const& coll, aod::McParticles const& mcParticles)
{
  auto bc = coll.bc_as<BCsWithRun2Info>();
  if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted))) {
    return;
  }
histos1.fill(HIST("mCentMCImpactFull"), coll.impactParameter());
   //if (coll.impactParameter() >= 3.5) {
   // return;
  //}
 histos1.fill(HIST("mCentMCImpact005"), coll.impactParameter());
  histos.fill(HIST("mVertexX"), coll.posX());
  histos.fill(HIST("mVertexY"), coll.posY());
  histos.fill(HIST("mVertexZ"), coll.posZ());
 

  for (auto const& h : mHistArrReset) {
    h->Reset();
  }

  countTracks = {0, 0, 0, 0, 0};
  fqEvent = {{{{{0, 0, 0, 0, 0, 0}}}}};
  binConEvent = {{{0, 0, 0, 0, 0}}};
  //using McParticle = std::decay_t<decltype(*particles.begin())>;
  //std::vector<McParticle> selectedParticles;
 // if (counter > 1) return;
  for (auto const& track : mcParticles) {
    histos1.fill(HIST("hRecoPtBefore"), track.pt());    
      histos1.fill(HIST("hRecoEtaBefore"), track.eta());  
    histos1.fill(HIST("hRecoPhiBefore"), track.phi());        
    if (std::abs(track.eta()) < 0.8 && track.isPhysicalPrimary()) {    
    histos1.fill(HIST("hRecoPtAfter"), track.pt());
    histos1.fill(HIST("hRecoEtaAfter"), track.eta());
    histos1.fill(HIST("hRecoPhiAfter"), track.phi());   
    checkpT(track);
    }
    histos.fill(HIST("mTrackSelected"), 1);
 
  }
   
 //counter++;
  for (auto iPt = 0; iPt < confNumPt; ++iPt) {
    if (countTracks[iPt] > 0) {
      mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
    } else {
      return;
    }
  }

  calculateMoments(mHistArrReset);
}

PROCESS_SWITCH(FactorialMoments, processMcGenRun2, "process MC Run2", false);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>>::iterator const& coll, TracksFMs const& tracks)
  {
    if ((!coll.alias_bit(kINT7)) || (!coll.sel7())) {
      return;
    }
    if (coll.centRun2V0M() < confCentCut.value[0] || coll.centRun2V0M() > confCentCut.value[1]) {
      return;
    }
    histos.fill(HIST("mVertexX"), coll.posX());
    histos.fill(HIST("mVertexY"), coll.posY());
    histos.fill(HIST("mVertexZ"), coll.posZ());
    histos.fill(HIST("mCentFT0M"), coll.centRun2V0M());

    for (auto const& h : mHistArrReset) {
      h->Reset();
    }

    countTracks = {0, 0, 0, 0, 0};
    fqEvent = {{{{{0, 0, 0, 0, 0, 0}}}}};
    binConEvent = {{{0, 0, 0, 0, 0}}};

    for (auto const& track : tracks) {
      if ((track.pt() < confPtMin) || (!track.isGlobalTrack()) || (track.tpcNClsFindable() < confMinTPCCls)) {
        continue;
      }
      histos.fill(HIST("mCollID"), track.collisionId());
      histos.fill(HIST("mEta"), track.eta());
      histos.fill(HIST("mPt"), track.pt());
      histos.fill(HIST("mPhi"), track.phi());
      histos.fill(HIST("mNFindableClsTPC"), track.tpcNClsFindable());
      histos.fill(HIST("mNClsTPC"), track.tpcNClsFound());
      histos.fill(HIST("mNClsITS"), track.itsNCls());
      histos.fill(HIST("mChi2TPC"), track.tpcChi2NCl());
      histos.fill(HIST("mChi2ITS"), track.itsChi2NCl());
      histos.fill(HIST("mChi2TRD"), track.trdChi2());
      histos.fill(HIST("mDCAxy"), track.dcaXY());
      histos.fill(HIST("mDCAx"), track.dcaZ());
      histos.fill(HIST("mDCAxyPt"), track.pt(), track.dcaXY());
      histos.fill(HIST("mDCAzPt"), track.pt(), track.dcaZ());
      histos.fill(HIST("mNSharedClsTPC"), track.tpcNClsShared());
      histos.fill(HIST("mCrossedRowsTPC"), track.tpcNClsCrossedRows());
      histos.fill(HIST("mNFinClsminusCRows"), track.tpcNClsFindableMinusCrossedRows());
      histos.fill(HIST("mNFractionShClsTPC"), track.tpcFractionSharedCls());
      histos.fill(HIST("mSharedClsvsPt"), track.pt(), track.tpcNClsShared());
      histos.fill(HIST("mSharedClsProbvsPt"), track.pt(), track.tpcFractionSharedCls() / track.tpcNClsCrossedRows());
      checkpT(track);
    }
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      if (countTracks[iPt] > 0) {
        mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
      } else {
        return;
      }
    }
    // Calculate the normalized factorial moments
    calculateMoments(mHistArrReset);
  }
  PROCESS_SWITCH(FactorialMoments, processRun2, "for RUN2", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FactorialMoments>(cfgc),
  };
}
