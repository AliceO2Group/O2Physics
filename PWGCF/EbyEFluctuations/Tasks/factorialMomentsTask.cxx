// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental Organization
// or submit itself to any jurisdiction.
///

/// \file factorialMomentsTask.cxx
/// \brief This task is for Normalized Factorial Moments Analysis: PRC 85, 044914 (2012), nucl-ex:1411.6083
/// \author Salman Malik,  Balwan Singh

#include "TRandom.h"
#include <TH1F.h>

// O2 includes
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
TH1D* tmpFqErr[6][5][52];

struct FactorialMomentsTask {
  Configurable<bool> applyCheckPtForRec{"applyCheckPtForRec", false, "Apply checkpT for reconstructed tracks"};
  Configurable<bool> applyCheckPtForMC{"applyCheckPtForMC", true, "Apply checkpT for MC-generated tracks"};
  Configurable<float> centralEta{"centralEta", 0.9, "eta limit for tracks"};
  Configurable<int> numPt{"numPt", 5, "number of pT bins"};
  Configurable<float> ptMin{"ptMin", 0.2f, "lower pT cut"};
  Configurable<float> dcaXY{"dcaXY", 2.4f, "DCA xy cut"};
  Configurable<float> dcaZ{"dcaZ", 2.0f, "DCA z cut"};
  Configurable<float> mintPCCls{"mintPCCls", 70.0f, "minimum number of TPC clusters"};
  Configurable<std::vector<int>> centLimits{"centLimits", {0, 5}, "centrality min and max"};
  Configurable<std::vector<float>> vertexXYZ{"vertexXYZ", {0.3f, 0.4f, 10.0f}, "vertex cuts"};
  Configurable<std::vector<float>> ptCuts{"ptCuts", {0.2f, 2.0f}, "pT cuts"};
  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", true, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> isApplyVertexITSTPC{"isApplyVertexITSTPC", true, "Enable VertexITSTPC cut"};
  Configurable<bool> isApplyVertexTOFmatched{"isApplyVertexTOFmatched", true, "Enable VertexTOFmatched cut"};
  Configurable<bool> isApplyVertexTRDmatched{"isApplyVertexTRDmatched", true, "Enable VertexTRDmatched cut"};
  Configurable<bool> isApplyExtraCorrCut{"isApplyExtraCorrCut", false, "Enable extra NPVtracks vs FTOC correlation cut"};
  Configurable<bool> isApplyExtraPhiCut{"isApplyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<bool> includeGlobalTracks{"includeGlobalTracks", false, "Enable Global Tracks"};
  Configurable<bool> includeTPCTracks{"includeTPCTracks", false, "TPC Tracks"};
  Configurable<bool> includeITSTracks{"includeITSTracks", false, "ITS Tracks"};
  Configurable<int> samplesize{"samplesize", 100, "Sample size"};
  Configurable<bool> useMC{"useMC", false, "Use MC information"};
  Configurable<int> reduceOutput{"reduceOutput", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};
  Filter filterTracks = (nabs(aod::track::eta) < centralEta) && (aod::track::pt >= ptMin) && (nabs(aod::track::dcaXY) < dcaXY) && (nabs(aod::track::dcaZ) < dcaZ);
  Filter filterCollisions = (nabs(aod::collision::posZ) < vertexXYZ.value[2]) && (nabs(aod::collision::posX) < vertexXYZ.value[0]) && (nabs(aod::collision::posY) < vertexXYZ.value[1]);
  Service<o2::framework::O2DatabasePDG> pdg;
  // Histograms
  HistogramRegistry histos1{
    "histos1",
    {
      {"hRecoPtBefore", "Reco pT before cuts;pt (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 20.0}}}},
      {"hGenPtBefore", "Gen pT before cuts;pt (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 20.0}}}},
      {"hRecoPtAfter", "Reco pT after cuts;pt (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 20.0}}}},
      {"hGenPtAfter", "Gen pT after cuts;pt (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 20.0}}}},
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
      {"mDiffPhi", "dphi between selected MC particles;dphi;Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mPrimariesPerEvent", "Primary MC particles per event;N_{primaries};Counts", {HistType::kTH1I, {{20000, 0, 20000}}}},
      {"mDiffPt", "dpt between selected MC particles;dpt (GeV/c);Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mDiffEta", "deta between selected MC particles;deta;Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mADiffPhi", "dphi between selected MC particles;dphi;Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mADiffPt", "dpt between selected MC particles;dpt (GeV/c);Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
      {"mADiffEta", "deta between selected MC particles;deta;Counts", {HistType::kTH1F, {{2000, -0.01, 0.01}}}},
    },
  };

  HistogramRegistry histos{
    "histos",
    {
      {"mChargeBefore", "Charge before MC cuts;charge;entries", {HistType::kTH1F, {{7, -3.5, 3.5}}}},
      {"mChargeAfter", "Charge after MC cuts;charge;entries", {HistType::kTH1F, {{7, -3.5, 3.5}}}},
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
      {"mPhi", "#phi", {HistType::kTH1F, {{100, 0, o2::constants::math::TwoPI}}}},
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
  static const int nBins = 52;
  double kMinCharge = 1e-6;
  static const int nfqOrder = 6;
  int countSamples = 0;
  int testc1 = 0, testc2 = 0, testc3 = 0;
  std::array<int, nBins> binningM;
  std::array<int, 5> countTracks{0, 0, 0, 0, 0};
  std::array<std::array<std::array<double, nBins>, 5>, 6> fqEvent;
  std::array<std::array<std::array<double, nBins>, 5>, 6> fqEventSampled;
  std::array<std::array<double, nBins>, 5> binConEvent;
  std::array<std::array<std::array<double, nBins>, 5>, 6> binConEventSampled;
  std::array<std::array<std::array<double, nBins>, 5>, 6> errorFq = {{{{{0, 0, 0, 0, 0}}}}};
  std::vector<std::shared_ptr<TH2>> mHistArrReset;
  std::vector<std::shared_ptr<TH1>> mHistArrQA;
  std::vector<std::shared_ptr<TH1>> mFqBinFinal;
  std::vector<std::shared_ptr<TH1>> mBinConFinal;
  std::vector<std::shared_ptr<TH1>> mFqBinFinalSampled;
  std::vector<std::shared_ptr<TH1>> mBinConFinalSampled;
  std::vector<std::shared_ptr<TH1>> mFqError;
  // max number of bins restricted to 5
  static constexpr std::array<std::string_view, 5>
    mbinNames{"bin1/", "bin2/", "bin3/", "bin4/", "bin5/"};
  void init(o2::framework::InitContext&)
  {
    // NOTE: check to make number of pt and the vector consistent
    if (numPt != static_cast<int>(ptCuts.value.size()) / 2) {
      for (int i = numPt; i < static_cast<int>(ptCuts.value.size() / 2); i++) {
        ptCuts.value[2 * i] = 0;
        ptCuts.value[2 * i + 1] = 0;
      }
    }
    AxisSpec axisPt[5] = {{100, -0.01, 3 * ptCuts.value[1], ""}, {100, -0.01, 3 * ptCuts.value[3], ""}, {100, -0.01, 3 * ptCuts.value[5], ""}, {100, -0.01, 3 * ptCuts.value[7], ""}, {100, -0.01, 3 * ptCuts.value[9], ""}}; // pT axis
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
    if (useMC) {
      auto mMcTrackSelected = std::get<std::shared_ptr<TH1>>(histos.add("mMcTrackSelected", "mcTrackSelection", HistType::kTH1D, {{5, 0.5, 5.5}}));
    }
    for (int iM = 0; iM < nBins; ++iM) {
      binningM[iM] = 2 * (iM + 2);
    }
    for (int iPt = 0; iPt < numPt; ++iPt) {
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mEta", iPt + 1), Form("#eta for bin %.2f-%.2f;#eta", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{1000, -2, 2}})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mPt", iPt + 1), Form("pT for bin %.2f-%.2f;pT", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {axisPt[iPt]})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mPhi", iPt + 1), Form("#phi for bin %.2f-%.2f;#phi", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{1000, 0,   o2::constants::math::TwoPI}})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mMultiplicity", iPt + 1), Form("Multiplicity for bin %.2f-%.2f;Multiplicity", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{1000, 0, 15000}})));
      for (int iM = 0; iM < nBins; ++iM) {
        auto mHistsR = std::get<std::shared_ptr<TH2>>(histos.add(Form("bin%i/Reset/mEtaPhi%i", iPt + 1, iM), Form("#eta#phi_%i for bin %.2f-%.2f;#eta;#phi", iM, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH2F, {{binningM[iM], -0.8, 0.8}, {binningM[iM], 0,  o2::constants::math::TwoPI}}));
        mHistArrReset.push_back(mHistsR);
        for (int iq = 0; iq < nfqOrder; ++iq) {
          tmpFqErr[iq][iPt][iM] = new TH1D(Form("tmpFqErr%i%i%i", iq, iPt, iM), Form("tmpFqErr%i%i%i", iq, iPt, iM), 100, 0, 10);
        }
      }
      for (int i = 0; i < nfqOrder; ++i) {
        auto mHistFq = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalFq%i_bin%i", i + 2, iPt + 1), Form("Final F_%i for bin %.2f-%.2f;M", i + 2, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mFqBinFinal.push_back(mHistFq);
        auto mHistAv = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalAvBin%i_bin%i", i + 2, iPt + 1), Form("Final AvBin_%i for bin %.2f-%.2f;M", i + 2, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mBinConFinal.push_back(mHistAv);
        auto mHistFqSampled = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalFq%iSampled_bin%i", i + 2, iPt + 1), Form("Final F_%i for bin %.2f-%.2f;M", i + 2, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mFqBinFinalSampled.push_back(mHistFqSampled);
        auto mHistAvSampled = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalAvBin%iSampled_bin%i", i + 2, iPt + 1), Form("Final AvBin_%i for bin %.2f-%.2f;M", i + 2, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mBinConFinalSampled.push_back(mHistAvSampled);

        auto mHistError = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFqError%i_bin%i", i + 2, iPt + 1), Form("Error for F_%i for bin %.2f-%.2f;M", i + 2, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mFqError.push_back(mHistError);
      }
    }
  }
  template <typename T>
  void checkpT(const T& track)
  {
    for (int iPt = 0; iPt < numPt; ++iPt) {
      if (track.pt() > ptCuts.value[2 * iPt] && track.pt() < ptCuts.value[2 * iPt + 1]) {
        float iphi = track.phi();
        iphi = gRandom->Gaus(iphi, o2::constants::math::TwoPI);
        iphi = RecoDecay::constrainAngle(iphi);

        mHistArrQA[iPt * 4]->Fill(track.eta());
        mHistArrQA[iPt * 4 + 1]->Fill(track.pt());
        mHistArrQA[iPt * 4 + 2]->Fill(track.phi());
        countTracks[iPt]++;

        for (int iM = 0; iM < nBins; ++iM) {
          mHistArrReset[iPt * nBins + iM]->Fill(track.eta(), track.phi());
        }
      }
    }
  }

  void calculateMoments(std::vector<std::shared_ptr<TH2>> hist)
  {
    double binContent = 0;
    countSamples++;
    bool compSample = kFALSE;
    if (countSamples == samplesize) {
      compSample = kTRUE;
      countSamples = 0;
    }
    // Calculate the normalized factorial moments
    for (int iPt = 0; iPt < numPt; ++iPt) {
      for (int iM = 0; iM < nBins; ++iM) {
        binContent = 0;
        double sumfqBin[6] = {0};

        for (int iEta = 1; iEta <= hist[iPt * nBins + iM]->GetNbinsX(); ++iEta) {
          for (int iPhi = 1; iPhi <= hist[iPt * nBins + iM]->GetNbinsY(); ++iPhi) {
            double binconVal = 0;
            binconVal = hist[iPt * nBins + iM]->GetBinContent(iEta, iPhi);
            binContent += binconVal;
            for (int iq = 0; iq < nfqOrder; ++iq) {
              double fqBin = 0;
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
        binConEvent[iPt][iM] = binContent / (std::pow(binningM[iM], 2));
        for (int iq = 0; iq < nfqOrder; ++iq) {
          if (sumfqBin[iq] > 0) {
            fqEvent[iq][iPt][iM] = sumfqBin[iq] / (std::pow(binningM[iM], 2));
            fqEventSampled[iq][iPt][iM] += fqEvent[iq][iPt][iM];
          }
          binConEventSampled[iq][iPt][iM] += binConEvent[iPt][iM];
          mFqBinFinal[iPt * 6 + iq]->Fill(iM, fqEvent[iq][iPt][iM]);
          mBinConFinal[iPt * 6 + iq]->Fill(iM, binConEvent[iPt][iM]);
          if (compSample) {
            mBinConFinalSampled[iPt * 6 + iq]->Fill(iM, binConEventSampled[iq][iPt][iM] / samplesize);

            double tmp = (fqEventSampled[iq][iPt][iM] / (samplesize)) / (std::pow(binConEventSampled[iq][iPt][iM] / (samplesize), (iq + 2)));
            mFqBinFinalSampled[iPt * 6 + iq]->Fill(iM, tmp);
            tmpFqErr[iq][iPt][iM]->Fill(tmp);
            errorFq[iq][iPt][iM] += std::pow(fqEventSampled[iq][iPt][iM] / (samplesize), 2);
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
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>::iterator const& coll, TracksFMs const& tracks)
  {
    // selection of events
    if (!coll.sel8()) {
      return;
    }
    if (isApplySameBunchPileup && !coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (isApplyGoodZvtxFT0vsPV && !coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (isApplyVertexITSTPC && !coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    if (coll.centFT0C() < centLimits.value[0] || coll.centFT0C() > centLimits.value[1]) {
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
      if (track.hasTPC()) {
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
    for (int iPt = 0; iPt < numPt; ++iPt) {
      if (countTracks[iPt] > 0) {
        mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
      }
    }
    calculateMoments(mHistArrReset);
  }
  PROCESS_SWITCH(FactorialMomentsTask, processRun3, "main process function", false);
  using CollisionCandidateMCRec = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  using TracksMc = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>>;
  void processMCRec(soa::Filtered<CollisionCandidateMCRec>::iterator const& coll, TracksMc const& colltracks, aod::McParticles const& mcParticles, aod::McCollisions const&)
  {
    if (!coll.has_mcCollision()) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 0);
    if (!coll.sel8()) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 1);
    if (isApplySameBunchPileup && !coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 2);
    if (isApplyGoodZvtxFT0vsPV && !coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 3);
    if (isApplyVertexITSTPC && !coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 4);
    if (coll.centFT0C() < centLimits.value[0] || coll.centFT0C() > centLimits.value[1]) {
      return;
    }
    histos.fill(HIST("mEventSelected"), 5);
    histos.fill(HIST("mVertexX"), coll.posX());
    histos.fill(HIST("mVertexY"), coll.posY());
    histos.fill(HIST("mVertexZ"), coll.posZ());
    histos.fill(HIST("mCentFT0C"), coll.centFT0C());
    for (auto const& h : mHistArrReset) {
      h->Reset();
    }
    countTracks = {0, 0, 0, 0, 0};
    fqEvent = {{{{{0, 0, 0, 0, 0, 0}}}}};
    binConEvent = {{{0, 0, 0, 0, 0}}};
    for (auto const& track : colltracks) {
      // if (track.hasITS())
       //if (track.hasTPC())
      //if (track.isGlobalTrack()) {
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
    auto mcParts = mcParticles.sliceBy(perMcCollision, coll.mcCollision().globalIndex());
    for (auto const& mc : mcParts) {
      int pdgCode = mc.pdgCode();
      auto pdgInfo = pdg->GetParticle(pdgCode);
      if (!pdgInfo) {
        continue;
      }
      double charge = pdgInfo->Charge();
      double physCharge = charge / 3.0;
      histos.fill(HIST("mChargeBefore"), physCharge);
      if (mc.isPhysicalPrimary() && std::abs(mc.eta()) < centralEta && std::abs(physCharge) >= kMinCharge) {
        histos.fill(HIST("mChargeAfter"), physCharge);
        histos.fill(HIST("mEta"), mc.eta());
        histos.fill(HIST("mPt"), mc.pt());
        histos.fill(HIST("mPhi"), mc.phi());
        if (applyCheckPtForMC && !applyCheckPtForRec) {
          checkpT(mc);
        }
      }
    }
    for (auto iPt = 0; iPt < numPt; ++iPt) {
      // if (countTracks[iPt] > 0)countTracks = {0, 0, 0, 0, 0};
      if (countTracks[iPt] > 0) {
        mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
      }
    }
    histos.fill(HIST("mEventSelected"), 6);
    // Calculate the normalized factorial moments
    calculateMoments(mHistArrReset);
  }
  PROCESS_SWITCH(FactorialMomentsTask, processMCRec, "main process function", false);
  using EventSelectionrun2 = soa::Join<aod::EvSels, aod::Mults, aod::CentRun2V0Ms, aod::CentRun2SPDTrks>;
  using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  using CollisionRecSimRun2 = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelectionrun2>>::iterator;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;
 void processMcRun2(CollisionRecSimRun2 const& coll,
                   aod::BCs const& ,
                   TracksRecSim const& tracks,
                   aod::McParticles const& mcParticles,
                   aod::McCollisions const&,
                   BCsWithRun2Info const&  )
 {
    auto bc = coll.bc_as<BCsWithRun2Info>();
    if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted))) {
      return;
    }
    if (coll.centRun2V0M() < centLimits.value[0] || coll.centRun2V0M() > centLimits.value[1]) {
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
      double recoCharge = (track.sign() != 0) ? track.sign() : 0.;
      if (std::abs(track.eta()) < centralEta && track.isGlobalTrack() && std::abs(recoCharge) >= kMinCharge) {
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
        if (applyCheckPtForRec && !applyCheckPtForMC) {
          checkpT(track);
        }
      }
    }
    auto mcParts = mcParticles.sliceBy(perMcCollision, coll.mcCollision().globalIndex());
    for (auto const& mc : mcParts) {
      int pdgCode = mc.pdgCode();
      auto pdgInfo = pdg->GetParticle(pdgCode);
      if (!pdgInfo) {
        continue;
      }
      double charge = pdgInfo->Charge();
      double physCharge = charge / 3.0;
      histos.fill(HIST("mChargeBefore"), physCharge);
      if (mc.isPhysicalPrimary() && std::abs(mc.eta()) < centralEta && std::abs(physCharge) >= kMinCharge) {
        histos.fill(HIST("mChargeAfter"), physCharge);
        histos.fill(HIST("mEta"), mc.eta());
        histos.fill(HIST("mPt"), mc.pt());
        histos.fill(HIST("mPhi"), mc.phi());
        if (applyCheckPtForMC && !applyCheckPtForRec) {
          checkpT(mc);
        }
      }
    }

    for (int iPt = 0; iPt < numPt; ++iPt) {
      if (countTracks[iPt] > 0) {
        mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
      } else {
        return;
      }
    }

    calculateMoments(mHistArrReset);
  }

  PROCESS_SWITCH(FactorialMomentsTask, processMcRun2, "process MC Run2", true);
  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>>::iterator const& coll, TracksFMs const& tracks)
  {
    if ((!coll.alias_bit(kINT7)) || (!coll.sel7())) {
      return;
    }
    if (coll.centRun2V0M() < centLimits.value[0] || coll.centRun2V0M() > centLimits.value[1]) {
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
      if ((track.pt() < ptMin) || (!track.isGlobalTrack()) || (track.tpcNClsFindable() < mintPCCls)) {
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
    for (int iPt = 0; iPt < numPt; ++iPt) {
      if (countTracks[iPt] > 0) {
        mHistArrQA[iPt * 4 + 3]->Fill(countTracks[iPt]);
      } else {
        return;
      }
    }
    // Calculate the normalized factorial moments
    calculateMoments(mHistArrReset);
  }
  PROCESS_SWITCH(FactorialMomentsTask, processRun2, "for RUN2", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FactorialMomentsTask>(cfgc),
  };
}
