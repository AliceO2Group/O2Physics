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
/// \file FactorialMomentsTask.cxx
/// \brief This task is for Normalized Factorial Moments Analysis: Hwa and Yang, Phys. Rev. C 85, 044914 (2012)
/// \author Salman Malik
/// \author Balwan Singh

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TRandom.h>
#include <TString.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FactorialMomentsTask {
  Configurable<bool> useITS{"useITS", false, "Select tracks with ITS"};
  Configurable<bool> useTPC{"useTPC", false, "Select tracks with TPC"};
  Configurable<bool> useGlobal{"useGlobal", true, "Select global tracks"};
  Configurable<bool> applyCheckPtForRec{"applyCheckPtForRec", false, "Apply checkpT for reconstructed tracks"};
  Configurable<bool> applyCheckPtForMC{"applyCheckPtForMC", true, "Apply checkpT for MC-generated tracks"};
  Configurable<bool> smearPhi{"smearPhi", false, "Randomise track phi with a Gaussian of sigma=2pi before filling the eta-phi lattices"};
  Configurable<bool> evSelNoITSROFrameBorder{"evSelNoITSROFrameBorder", true, "ITSROFrame border event selection cut"};
  Configurable<bool> evSelNoTimeFrameBorder{"evSelNoTimeFrameBorder", true, "TimeFrame border event selection cut"};
  Configurable<float> centralEta{"centralEta", 0.9, "eta limit for tracks"};
  Configurable<int> numPt{"numPt", 5, "number of pT bins"};
  Configurable<float> ptMin{"ptMin", 0.2f, "lower pT cut"};
  Configurable<float> dcaXY{"dcaXY", 0.1f, "DCA xy cut"};
  Configurable<float> dcaZ{"dcaZ", 1.0f, "DCA z cut"};
  Configurable<float> maxTPCChi2NCl{"maxTPCChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> maxITSChi2NCl{"maxITSChi2NCl", 40.0f, "Maximum ITSchi2NCl"};
  Configurable<float> minTPCNCls{"minTPCNCls", 70.0f, "minimum number of TPC clusters"};
  Configurable<std::vector<int>> centLimits{"centLimits", {0, 5}, "centrality min and max"};
  Configurable<std::vector<float>> vertexXYZ{"vertexXYZ", {0.3f, 0.4f, 10.0f}, "vertex cuts"};
  Configurable<std::vector<float>> ptCuts{"ptCuts", {0.2f, 2.0f}, "pT cuts"};
  Configurable<bool> evSelNoSameBunchPileup{"evSelNoSameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> evSelGoodZvtxFT0vsPV{"evSelGoodZvtxFT0vsPV", true, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> evSelGoodITSLayersAll{"evSelGoodITSLayersAll", true, "Remove time interval with dead ITS zone"};
  Configurable<bool> evSelVertexITSTPC{"evSelVertexITSTPC", true, "Enable VertexITSTPC cut"};
  Configurable<bool> evSelVertexTOFmatched{"evSelVertexTOFmatched", true, "Enable VertexTOFmatched cut"};
  Configurable<bool> evSelVertexTRDmatched{"evSelVertexTRDmatched", true, "Enable VertexTRDmatched cut"};
  Configurable<bool> applyExtraCorrCut{"applyExtraCorrCut", false, "Enable extra NPVtracks vs FTOC correlation cut"};
  Configurable<bool> applyExtraPhiCut{"applyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<bool> includeGlobalTracks{"includeGlobalTracks", false, "Enable Global Tracks"};
  Configurable<bool> includeTPCTracks{"includeTPCTracks", false, "TPC Tracks"};
  Configurable<bool> includeITSTracks{"includeITSTracks", false, "ITS Tracks"};
  Configurable<int> nSlots{"nSlots", 1250, "Number of hash slots for merge-safe subsampling (slot = splitmix64(globalBC) % nSlots); must be identical for every job"};
  Configurable<bool> useMC{"useMC", false, "Use MC information"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "Require global track in filter"};
  Configurable<int> minITSClusters{"minITSClusters", 6, "Minimum Number of ITS cluster"};
  Configurable<int> minTPCClusters{"minTPCClusters", 80, "Minimum Number of TPC cluster"};
  Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 70, "Minimum Number of TPC crossed-rows"};
  Configurable<float> minTPCCrossedRowsOverFindableCls{"minTPCCrossedRowsOverFindableCls", 0.8, "Minimum ratio of crossed rows over findable clusters TPC"};
  Configurable<int> reduceOutput{"reduceOutput", 0, "Suppress info level output (0 = all output, 1 = per collision, 2 = none)"};
  Filter filterTracks = (nabs(aod::track::eta) < centralEta) && (aod::track::pt >= ptMin);
  Filter filterCollisions = (nabs(aod::collision::posZ) < vertexXYZ.value[2]) && (nabs(aod::collision::posX) < vertexXYZ.value[0]) && (nabs(aod::collision::posY) < vertexXYZ.value[1]);
  Service<o2::framework::O2DatabasePDG> pdg{};
  // Histograms
  HistogramRegistry histos{
    "histos",
    {

      {"mtpcsignalvspt", "tpcsignal vs #pt", {HistType::kTH2F, {{900, 0, 10}, {1400, 0, 1400}}}},
      {"mChargeBefore", "Charge before MC cuts;charge;entries", {HistType::kTH1F, {{7, -3.5, 3.5}}}},
      {"mChargeAfter", "Charge after MC cuts;charge;entries", {HistType::kTH1F, {{7, -3.5, 3.5}}}},
      {"mCollID", "collisionID", {HistType::kTH1I, {{1000, -10000, 10000}}}},
      {"mCentFT0M", "centFT0M / Run2 V0M", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mCentFV0A", "centFV0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mCentFT0A", "centFT0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mCentFT0C", "centFT0C", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mVertexX", "vertexX", {HistType::kTH1F, {{100, -20, 20}}}},
      {"mVertexY", "vertexY", {HistType::kTH1F, {{100, -20, 20}}}},
      {"mVertexZ", "vertexZ", {HistType::kTH1F, {{100, -20, 20}}}},
      {"mEta", "#eta", {HistType::kTH1F, {{1000, -2, 2}}}},
      {"mPt", "#pt", {HistType::kTH1F, {{1000, -0.01, 50}}}},
      {"mPhi", "#phi", {HistType::kTH1F, {{100, 0, o2::constants::math::TwoPI}}}},
      {"mNFindableClsTPC", "findable TPC clusters;findable clusters", {HistType::kTH1F, {{100, 0, 200}}}},
      {"mNClsTPC", "number of clusters TPC; nClusters TPC", {HistType::kTH1F, {{100, 0, 200}}}},
      {"mNClsITS", "number of clusters ITS; nClusters ITS", {HistType::kTH1F, {{100, 0, 10}}}},
      {"mChi2TPC", "chi2 TPC", {HistType::kTH1F, {{100, 0, 10}}}},
      {"mChi2ITS", "chi2 ITS", {HistType::kTH1F, {{100, 0, 10}}}},
      {"mChi2TRD", "chi2 TRD", {HistType::kTH1F, {{100, 0, 100}}}},
      {"mDCAxy", "DCA xy", {HistType::kTH1F, {{500, -0.8, 0.8}}}},
      {"mDCAz", "DCA z", {HistType::kTH1F, {{500, -2.0, 2.0}}}},
      {"mDCAxyPt", "DCA xy vs #pt;#pt;DCAxy", {HistType::kTH2F, {{100, 0, 20}, {500, -0.5, 0.5}}}},
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
  float collisionZ = 0.f;
  static constexpr double kDcaXY0 = 0.0105;
  static constexpr double kDcaXY1 = 0.035;
  static constexpr double kDcaXY2 = 1.1;
  static const int nBins = 52;
  static constexpr double kMinCharge = 1e-6;
  static const int nfqOrder = 6;
  static constexpr int kEvAll = 1;
  static constexpr int kEvSel8 = 2;
  static constexpr int kEvITSROFrame = 3;
  static constexpr int kEvTimeFrame = 4;
  static constexpr int kEvSameBunch = 5;
  static constexpr int kEvGoodITS = 6;
  static constexpr int kEvGoodZvtx = 7;
  static constexpr int kEvVertexITSTPC = 8;
  static constexpr int kEvCent = 9;
  static constexpr int kEvAccepted = 10;
  static constexpr int kNEventSteps = 10;
  std::array<int, nBins> binningM{};
  std::array<int, 5> countTracks{};
  std::array<std::array<std::array<double, nBins>, 5>, nfqOrder> fqEvent{};
  std::array<std::array<double, nBins>, 5> binConEvent{};
  std::vector<std::shared_ptr<TH2>> mHistArrReset;
  std::vector<std::shared_ptr<TH1>> mHistEtaQA;
  std::vector<std::shared_ptr<TH1>> mHistPtQA;
  std::vector<std::shared_ptr<TH1>> mHistPhiQA;
  std::vector<std::shared_ptr<TH1>> mHistMultQA;
  std::vector<std::shared_ptr<TH3>> mHistArrEff;
  std::vector<std::shared_ptr<TH1>> mFqBinFinal;
  std::vector<std::shared_ptr<TH1>> mBinConFinal;
  // merge-safe slot output: additive per-slot sums, hadd-exact
  std::shared_ptr<TH1> mSlotNev;
  std::vector<std::shared_ptr<TH2>> mSlotFq;
  std::vector<std::shared_ptr<TH2>> mSlotOcc;
  uint64_t mEventGlobalBC = 0;

  void init(o2::framework::InitContext&)
  {
    if (numPt.value < 1 || numPt.value > 5) {
      LOG(fatal) << "numPt must be between 1 and 5, got " << numPt.value;
    }
    if (static_cast<int>(ptCuts.value.size()) < 2 * numPt.value) {
      LOG(fatal) << "ptCuts needs " << 2 * numPt.value << " values for numPt=" << numPt.value
                 << ", got " << ptCuts.value.size();
    }
    if (nSlots.value < 1 || nSlots.value > 1000000) {
      LOG(fatal) << "nSlots must be in [1, 1000000], got " << nSlots.value;
    }
    if (centLimits.value.size() != 2 || vertexXYZ.value.size() != 3) {
      LOG(fatal) << "centLimits needs 2 and vertexXYZ needs 3 values";
    }
    if (applyCheckPtForRec.value == applyCheckPtForMC.value) {
      LOG(fatal) << "exactly one of applyCheckPtForRec / applyCheckPtForMC must be enabled";
    }
    TString cfg = Form("FMtask centLimits=%d,%d numPt=%d centralEta=%g ptMin=%g dcaXY=%g dcaZ=%g useMC=%d nfqOrder=%d smearPhi=%d nSlots=%d",
                       centLimits.value[0], centLimits.value[1], numPt.value, centralEta.value,
                       ptMin.value, dcaXY.value, dcaZ.value, useMC.value ? 1 : 0, nfqOrder, smearPhi.value ? 1 : 0, nSlots.value);
    cfg += Form(" vertexXYZ=%g,%g,%g", vertexXYZ.value[0], vertexXYZ.value[1], vertexXYZ.value[2]);
    cfg += " ptCuts=";
    for (int i = 0; i < static_cast<int>(ptCuts.value.size()); i++) {
      cfg += Form("%s%g", i > 0 ? "," : "", ptCuts.value[i]);
    }
    histos.add("metaConfig", cfg.Data(), HistType::kTH1D, {{1, 0.5, 1.5}});
    auto mEventSelected = std::get<std::shared_ptr<TH1>>(histos.add("mEventSelected", "eventSelected", HistType::kTH1D, {{kNEventSteps, 0.5, kNEventSteps + 0.5}}));
    const std::array<const char*, kNEventSteps> evLabels = {"all", "sel8", "kNoITSROFrameBorder", "kNoTimeFrameBorder", "kNoSameBunchPileup", "kIsGoodITSLayersAll", "kIsGoodZvtxFT0vsPV", "kIsVertexITSTPC", "centrality", "accepted"};
    for (int iStep = 0; iStep < kNEventSteps; ++iStep) {
      mEventSelected->GetXaxis()->SetBinLabel(iStep + 1, evLabels[iStep]);
    }
    // merge-safe slot output: additive per-slot sums, hadd-exact (see slotpool.C)
    const AxisSpec slotAxis{nSlots.value, -0.5, nSlots.value - 0.5, ""};
    const AxisSpec slotMAxis{nBins, -0.5, nBins - 0.5, ""};
    mSlotNev = std::get<std::shared_ptr<TH1>>(histos.add("slots/nEv", "events per slot;slot;events", HistType::kTH1F, {slotAxis}));
    // bin1 = merged jobs, bin2 = sum of their nSlots: mismatch witness (hadd keeps the first axis silently)
    auto slotParams = std::get<std::shared_ptr<TH1>>(histos.add("slots/params", "slot campaign parameters;jobs;entries", HistType::kTH1D, {{2, 0.5, 2.5}}));
    slotParams->SetBinContent(1, 1.0);
    slotParams->SetBinContent(2, static_cast<double>(nSlots.value));
    for (int iPt = 0; iPt < numPt; ++iPt) {
      mSlotOcc.push_back(std::get<std::shared_ptr<TH2>>(histos.add(Form("slots/sumOcc_pt%d", iPt + 1),
                                                                   Form("sum of occupancy per slot for bin %.2f-%.2f;slot;M", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]),
                                                                   HistType::kTH2F, {slotAxis, slotMAxis})));
      for (int iq = 0; iq < nfqOrder; ++iq) {
        mSlotFq.push_back(std::get<std::shared_ptr<TH2>>(histos.add(Form("slots/sumFq_pt%d_o%d", iPt + 1, iq + 2),
                                                                    Form("sum of f_%d per slot for bin %.2f-%.2f;slot;M", iq + 2, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]),
                                                                    HistType::kTH2F, {slotAxis, slotMAxis})));
      }
    }
    for (int iM = 0; iM < nBins; ++iM) {
      binningM[iM] = 2 * (iM + 2);
    }
    for (int iPt = 0; iPt < numPt; ++iPt) {
      const AxisSpec ptAxis{100, -0.01, 3 * ptCuts.value[2 * iPt + 1], ""};
      mHistArrEff.push_back(std::get<std::shared_ptr<TH3>>(histos.add(Form("bin%i/m3DVtxZetaPhi", iPt + 1), Form("#eta #phi #vtxz for bin %.2f-%.2f;vz;#eta;#phi", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH3F, {{20, -10, 10}, {16, -0.8, +0.8}, {100, 0., o2::constants::math::TwoPI}})));
      mHistEtaQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mEta", iPt + 1), Form("#eta for bin %.2f-%.2f;#eta", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{1000, -2, 2}})));
      mHistPtQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mPt", iPt + 1), Form("pT for bin %.2f-%.2f;pT", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {ptAxis})));
      mHistPhiQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mPhi", iPt + 1), Form("#phi for bin %.2f-%.2f;#phi", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{1000, 0, o2::constants::math::TwoPI}})));
      mHistMultQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mMultiplicity", iPt + 1), Form("Multiplicity for bin %.2f-%.2f;Multiplicity", ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{1000, 0, 15000}})));
      for (int iM = 0; iM < nBins; ++iM) {
        auto mHistsR = std::get<std::shared_ptr<TH2>>(histos.add(Form("bin%i/Reset/mEtaPhi%i", iPt + 1, iM), Form("#eta#phi_%i for bin %.2f-%.2f;#eta;#phi", iM, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH2F, {{binningM[iM], -0.8, 0.8}, {binningM[iM], 0, o2::constants::math::TwoPI}}));
        mHistArrReset.push_back(mHistsR);
      }
      for (int i = 0; i < nfqOrder; ++i) {
        auto mHistFq = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalFq%i_bin%i", i + 2, iPt + 1), Form("Final F_%i for bin %.2f-%.2f;M", i + 2, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mFqBinFinal.push_back(mHistFq);
        auto mHistAv = std::get<std::shared_ptr<TH1>>(histos.add(Form("mFinalAvBin%i_bin%i", i + 2, iPt + 1), Form("Final AvBin_%i for bin %.2f-%.2f;M", i + 2, ptCuts.value[2 * iPt], ptCuts.value[2 * iPt + 1]), HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}}));
        mBinConFinal.push_back(mHistAv);
      }
    }
  }
  template <typename T>
  void checkpT(const T& track)
  {
    for (auto iPt = 0; iPt < numPt; ++iPt) {
      if (track.pt() > ptCuts.value[2 * iPt] && track.pt() < ptCuts.value[2 * iPt + 1]) {
        float iphi = track.phi();
        if (smearPhi) {
          iphi = gRandom->Gaus(iphi, o2::constants::math::TwoPI);
          iphi = RecoDecay::constrainAngle(iphi, 0.);
        }
        mHistArrEff[iPt]->Fill(collisionZ, track.eta(), iphi);
        mHistEtaQA[iPt]->Fill(track.eta());
        mHistPtQA[iPt]->Fill(track.pt());
        mHistPhiQA[iPt]->Fill(iphi);
        countTracks[iPt]++;
        for (auto iM = 0; iM < nBins; ++iM) {
          mHistArrReset[iPt * nBins + iM]->Fill(track.eta(), iphi);
        }
      }
    }
  }
  // splitmix64: deterministic stateless hash of the global BC
  static uint64_t slotOf(uint64_t x)
  {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  }
  // additive fold of the current accepted event into its slot (never Fill with weights: no fSumw2, hadd-exact)
  void accumulateSlot()
  {
    const int slot = static_cast<int>(slotOf(mEventGlobalBC) % static_cast<uint64_t>(nSlots.value));
    mSlotNev->SetBinContent(slot + 1, mSlotNev->GetBinContent(slot + 1) + 1.);
    for (int iPt = 0; iPt < numPt; ++iPt) {
      auto* occ = mSlotOcc[iPt].get();
      for (int iM = 0; iM < nBins; ++iM) {
        occ->SetBinContent(slot + 1, iM + 1, occ->GetBinContent(slot + 1, iM + 1) + binConEvent[iPt][iM]);
        for (int iq = 0; iq < nfqOrder; ++iq) {
          auto* fq = mSlotFq[iPt * nfqOrder + iq].get();
          fq->SetBinContent(slot + 1, iM + 1, fq->GetBinContent(slot + 1, iM + 1) + fqEvent[iq][iPt][iM]);
        }
      }
    }
  }
  void calculateMoments(const std::vector<std::shared_ptr<TH2>>& hist)
  {
    // normalized factorial moments per accepted event, per pT bin and M
    for (int iPt = 0; iPt < numPt; ++iPt) {
      for (int iM = 0; iM < nBins; ++iM) {
        double binContent = 0;
        std::array<double, nfqOrder> sumfqBin{0.};
        const double m2 = static_cast<double>(binningM[iM]) * static_cast<double>(binningM[iM]);
        const auto& histM = hist[iPt * nBins + iM];

        for (int iEta = 1; iEta <= histM->GetNbinsX(); ++iEta) {
          for (int iPhi = 1; iPhi <= histM->GetNbinsY(); ++iPhi) {
            const double binconVal = histM->GetBinContent(iEta, iPhi);
            binContent += binconVal;
            for (int iq = 0; iq < nfqOrder; ++iq) {
              double fqBin = 0;
              if (binconVal >= iq + 2) {
                const int nOcc = static_cast<int>(binconVal);
                fqBin = 1.0;
                for (int k = 0; k < iq + 2; ++k) {
                  fqBin *= (nOcc - k);
                }
              }
              sumfqBin[iq] += fqBin;
            }
          }
        }
        binConEvent[iPt][iM] = binContent / m2;
        for (int iq = 0; iq < nfqOrder; ++iq) {
          fqEvent[iq][iPt][iM] = sumfqBin[iq] / m2;
          mFqBinFinal[iPt * nfqOrder + iq]->Fill(iM, fqEvent[iq][iPt][iM]);
          mBinConFinal[iPt * nfqOrder + iq]->Fill(iM, binConEvent[iPt][iM]);
        }
      }
    }
    accumulateSlot();
  }

  template <typename T>
  bool passDcaCut(const T& track) const
  {
    return std::fabs(track.dcaXY()) < (kDcaXY0 + kDcaXY1 / std::pow(track.pt(), kDcaXY2));
  }

  template <typename T>
  void fillTrackQA(const T& track)
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
    histos.fill(HIST("mDCAz"), track.dcaZ());
    histos.fill(HIST("mtpcsignalvspt"), track.pt(), track.tpcSignal());
    histos.fill(HIST("mDCAxyPt"), track.pt(), track.dcaXY());
    histos.fill(HIST("mDCAzPt"), track.pt(), track.dcaZ());
    histos.fill(HIST("mNSharedClsTPC"), track.tpcNClsShared());
    histos.fill(HIST("mCrossedRowsTPC"), track.tpcNClsCrossedRows());
    histos.fill(HIST("mNFinClsminusCRows"), track.tpcNClsFindableMinusCrossedRows());
    histos.fill(HIST("mNFractionShClsTPC"), track.tpcFractionSharedCls());
    histos.fill(HIST("mSharedClsvsPt"), track.pt(), track.tpcNClsShared());
    const double crossedRows = track.tpcNClsCrossedRows();
    histos.fill(HIST("mSharedClsProbvsPt"), track.pt(), crossedRows > 0. ? track.tpcFractionSharedCls() / crossedRows : 0.);
  }

  template <typename T>
  void fillMcQA(const T& mcParts)
  {
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
  }

  template <typename T>
  bool selectEventRun3(const T& coll)
  {
    histos.fill(HIST("mEventSelected"), kEvAll);
    if (!coll.sel8()) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvSel8);
    if (evSelNoITSROFrameBorder && !(coll.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvITSROFrame);
    if (evSelNoTimeFrameBorder && !(coll.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvTimeFrame);
    if (evSelNoSameBunchPileup && !coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvSameBunch);
    if (evSelGoodITSLayersAll && !(coll.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvGoodITS);
    if (evSelGoodZvtxFT0vsPV && !coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvGoodZvtx);
    if (coll.centFT0C() < centLimits.value[0] || coll.centFT0C() > centLimits.value[1]) {
      return false;
    }
    collisionZ = coll.posZ();
    histos.fill(HIST("mEventSelected"), kEvCent);
    histos.fill(HIST("mEventSelected"), kEvAccepted);
    return true;
  }

  template <typename T>
  bool selectEventMCRec(const T& coll)
  {
    histos.fill(HIST("mEventSelected"), kEvAll);
    if (!coll.sel8()) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvSel8);
    if (evSelNoSameBunchPileup && !coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvSameBunch);
    if (evSelGoodZvtxFT0vsPV && !coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvGoodZvtx);
    if (evSelVertexITSTPC && !coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("mEventSelected"), kEvVertexITSTPC);
    if (coll.centFT0C() < centLimits.value[0] || coll.centFT0C() > centLimits.value[1]) {
      return false;
    }
    collisionZ = coll.posZ();
    histos.fill(HIST("mEventSelected"), kEvCent);
    histos.fill(HIST("mEventSelected"), kEvAccepted);
    return true;
  }

  void beginEvent()
  {
    for (int iPt = 0; iPt < numPt; ++iPt) {
      if (countTracks[iPt] == 0) {
        continue;
      }
      for (int iM = 0; iM < nBins; ++iM) {
        mHistArrReset[iPt * nBins + iM]->Reset();
      }
    }
    countTracks.fill(0);
    fqEvent = {};
    binConEvent = {};
  }

  void endEvent()
  {
    for (int iPt = 0; iPt < numPt; ++iPt) {
      if (countTracks[iPt] > 0) {
        mHistMultQA[iPt]->Fill(countTracks[iPt]);
      }
    }
    calculateMoments(mHistArrReset);
  }

  using TracksFMs = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>::iterator const& coll, TracksFMs const& tracks, aod::BCs const&)
  {
    mEventGlobalBC = coll.bc_as<aod::BCs>().globalBC();
    if (!selectEventRun3(coll)) {
      return;
    }
    histos.fill(HIST("mVertexX"), coll.posX());
    histos.fill(HIST("mVertexY"), coll.posY());
    histos.fill(HIST("mVertexZ"), coll.posZ());
    histos.fill(HIST("mCentFT0M"), coll.centFT0M());
    histos.fill(HIST("mCentFV0A"), coll.centFV0A());
    histos.fill(HIST("mCentFT0A"), coll.centFT0A());
    histos.fill(HIST("mCentFT0C"), coll.centFT0C());
    beginEvent();
    for (auto const& track : tracks) {
      if (useITS && !track.hasITS()) {
        continue;
      }
      if (useTPC && !track.hasTPC()) {
        continue;
      }
      if (useGlobal && !track.isGlobalTrack()) {
        continue;
      }
      if (!passDcaCut(track)) {
        continue;
      }
      fillTrackQA(track);
      checkpT(track);
    }
    endEvent();
  }
  PROCESS_SWITCH(FactorialMomentsTask, processRun3, "main process function", false);
  using CollisionCandidateMCRec = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  using TracksMc = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>>;
  void processMCRec(soa::Filtered<CollisionCandidateMCRec>::iterator const& coll, TracksMc const& colltracks, aod::McParticles const& mcParticles, aod::McCollisions const&, aod::BCs const&)
  {
    mEventGlobalBC = coll.bc_as<aod::BCs>().globalBC();
    if (!coll.has_mcCollision()) {
      return;
    }
    if (!selectEventMCRec(coll)) {
      return;
    }
    histos.fill(HIST("mVertexX"), coll.posX());
    histos.fill(HIST("mVertexY"), coll.posY());
    histos.fill(HIST("mVertexZ"), coll.posZ());
    histos.fill(HIST("mCentFT0C"), coll.centFT0C());
    beginEvent();
    for (auto const& track : colltracks) {
      if (useITS && !track.hasITS()) {
        continue;
      }
      if (useTPC && !track.hasTPC()) {
        continue;
      }
      if (useGlobal && !track.isGlobalTrack()) {
        continue;
      }
      if (!passDcaCut(track)) {
        continue;
      }
      fillTrackQA(track);
      if (applyCheckPtForRec && !applyCheckPtForMC) {
        checkpT(track);
      }
    }
    auto mcParts = mcParticles.sliceBy(perMcCollision, coll.mcCollision().globalIndex());
    fillMcQA(mcParts);
    endEvent();
  }
  PROCESS_SWITCH(FactorialMomentsTask, processMCRec, "main process function", false);
  using EventSelectionrun2 = soa::Join<aod::EvSels, aod::Mults, aod::CentRun2V0Ms, aod::CentRun2SPDTrks>;
  using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  using CollisionRecSimRun2 = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelectionrun2>>::iterator;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;
  void processMcRun2(CollisionRecSimRun2 const& coll,
                     aod::BCs const&,
                     TracksRecSim const& tracks,
                     aod::McParticles const& mcParticles,
                     aod::McCollisions const&,
                     BCsWithRun2Info const&)
  {
    auto bc = coll.bc_as<BCsWithRun2Info>();
    mEventGlobalBC = bc.globalBC();
    if (!(static_cast<bool>(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))) {
      return;
    }
    if (coll.centRun2V0M() < centLimits.value[0] || coll.centRun2V0M() > centLimits.value[1]) {
      return;
    }
    collisionZ = coll.posZ();
    histos.fill(HIST("mVertexX"), coll.posX());
    histos.fill(HIST("mVertexY"), coll.posY());
    histos.fill(HIST("mVertexZ"), coll.posZ());
    histos.fill(HIST("mCentFT0M"), coll.centRun2V0M());
    beginEvent();
    for (auto const& track : tracks) {
      double recoCharge = (track.sign() != 0) ? track.sign() : 0.;
      if (std::abs(track.eta()) < centralEta && track.isGlobalTrack() && std::abs(recoCharge) >= kMinCharge) {
        fillTrackQA(track);
        if (applyCheckPtForRec && !applyCheckPtForMC) {
          checkpT(track);
        }
      }
    }
    auto mcParts = mcParticles.sliceBy(perMcCollision, coll.mcCollision().globalIndex());
    fillMcQA(mcParts);
    endEvent();
  }

  PROCESS_SWITCH(FactorialMomentsTask, processMcRun2, "process MC Run2", false);
  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>>::iterator const& coll, TracksFMs const& tracks, aod::BCs const&)
  {
    mEventGlobalBC = coll.bc_as<aod::BCs>().globalBC();
    if ((!coll.alias_bit(kINT7)) || (!coll.sel7())) {
      return;
    }
    if (coll.centRun2V0M() < centLimits.value[0] || coll.centRun2V0M() > centLimits.value[1]) {
      return;
    }
    collisionZ = coll.posZ();
    histos.fill(HIST("mVertexX"), coll.posX());
    histos.fill(HIST("mVertexY"), coll.posY());
    histos.fill(HIST("mVertexZ"), coll.posZ());
    histos.fill(HIST("mCentFT0M"), coll.centRun2V0M());
    beginEvent();
    for (auto const& track : tracks) {
      if ((track.pt() < ptMin) || (!track.isGlobalTrack()) || (track.tpcNClsFindable() < minTPCNCls)) {
        continue;
      }
      fillTrackQA(track);
      checkpT(track);
    }
    endEvent();
  }
  PROCESS_SWITCH(FactorialMomentsTask, processRun2, "for RUN2", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<FactorialMomentsTask>(context),
  };
}
