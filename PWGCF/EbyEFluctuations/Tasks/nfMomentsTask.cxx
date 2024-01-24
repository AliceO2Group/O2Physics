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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <iostream>
#include <TH1F.h>
#include <array>

using std::array;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(mEvMult, eventMult, Int_t);
DECLARE_SOA_COLUMN(mVertexZ, vertexZ, Float_t);
DECLARE_SOA_COLUMN(mCentrality, cent, Float_t);
DECLARE_SOA_COLUMN(mBinMult, binmult, Int_t[5]);
DECLARE_SOA_COLUMN(mAvbCon, avbcontent, Double_t[200]);
DECLARE_SOA_COLUMN(mValM, m, Int_t[200]);
DECLARE_SOA_COLUMN(mPtBin, ptbin, Int_t[200]);
DECLARE_SOA_COLUMN(mF2, f2, Double_t[200]);
DECLARE_SOA_COLUMN(mF3, f3, Double_t[200]);
DECLARE_SOA_COLUMN(mF4, f4, Double_t[200]);
DECLARE_SOA_COLUMN(mF5, f5, Double_t[200]);
DECLARE_SOA_COLUMN(mF6, f6, Double_t[200]);
DECLARE_SOA_COLUMN(mF7, f7, Double_t[200]);
} // namespace full

DECLARE_SOA_TABLE(NFMFULL, "AOD", "NFMFULL", o2::soa::Index<>, full::mEvMult, full::mVertexZ, full::mCentrality, full::mBinMult, full::mAvbCon, full::mValM, full::mPtBin, full::mF2, full::mF3, full::mF4, full::mF5, full::mF6, full::mF7);
} // namespace o2::aod

// Write information to the tree
struct FactorialMoments {
  Produces<o2::aod::NFMFULL> tableFacMoments;

  Configurable<Float_t> confEta{"centralEta", 0.8, "eta limit for tracks"};
  Configurable<Int_t> confNumPt{"numPt", 5, "number of pT bins"};
  Configurable<Float_t> confPtMin{"ptMin", 0.2f, "lower pT cut"};
  Configurable<Float_t> confDCAxy{"dcaXY", 2.4f, "DCA xy cut"};
  Configurable<Float_t> confDCAz{"dcaZ", 3.2f, "DCA z cut"};
  Configurable<Float_t> confMinTPCCls{"minTPCCls", 70.0f, "minimum number of TPC clusters"};
  Configurable<std::vector<Int_t>> confCentCut{"centLimits", {0, 5}, "centrality min and max"};
  Configurable<std::vector<Float_t>> confVertex{"vertexXYZ", {0.3f, 0.4f, 10.0f}, "vertex cuts"};
  Configurable<std::vector<Float_t>> confPtBins{"ptCuts", {0.4f, 1.0f, 0.4f, 0.6f, 0.8f, 1.0f, 0.4f, 2.0f, 0.4f, 5.0f}, "pT cuts"};

  Filter filterTracks = (nabs(aod::track::eta) < confEta) && (aod::track::pt > confPtMin) && (nabs(aod::track::dcaXY) < confDCAxy) && (nabs(aod::track::dcaZ) < confDCAz);
  Filter filterCollisions = (nabs(aod::collision::posZ) < confVertex.value[2]) && (nabs(aod::collision::posX) < confVertex.value[0]) && (nabs(aod::collision::posY) < confVertex.value[1]);

  // Histograms
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
      {"mDCAxy", "DCA xy", {HistType::kTH1F, {{100, -0.8, 0.8}}}},
      {"mDCAx", "DCA z", {HistType::kTH1F, {{100, -2.0, 2.0}}}},
      {"mDCAxyPt", "DCA xy vs #pt;#pt;DCAxy", {HistType::kTH2F, {{100, 0, 20}, {100, -0.5, 0.5}}}},
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

  array<Int_t, 40> binningM;
  array<Int_t, 5> countTracks{0, 0, 0, 0, 0};
  array<array<array<Double_t, 40>, 5>, 6> fqEvent;
  array<array<Double_t, 40>, 5> averagebinCon;
  array<array<Int_t, 40>, 5> valueM;
  array<array<Int_t, 40>, 5> mPt;
  std::vector<std::shared_ptr<TH2>> mHistArrReset;
  std::vector<std::shared_ptr<TH2>> mHistArrUReset;
  std::vector<std::shared_ptr<TH1>> mHistArrQA;
  // max number of bins restricted to 5
  static constexpr array<std::string_view, 5> mbinNames{"bin1/", "bin2/", "bin3/", "bin4/", "bin5/"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec axisPt[5] = {{100, -0.01, 3 * confPtBins.value[1], ""}, {100, -0.01, 3 * confPtBins.value[3], ""}, {100, -0.01, 3 * confPtBins.value[5], ""}, {100, -0.01, 3 * confPtBins.value[7], ""}, {100, -0.01, 3 * confPtBins.value[9], ""}}; // pT axis

    if (confNumPt != static_cast<int>(confPtBins.value.size()) / 2) {
      throw std::runtime_error("Number of pT bins in confNumPt does not match the number of pT cuts in config.json!");
    }
    for (auto iM = 0; iM < 40; ++iM) {
      binningM[iM] = 2 * (iM + 2);
    }
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mEta", iPt + 1), Form("#eta for bin %.2f-%.2f;#eta", confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{1000, -2, 2}})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mPt", iPt + 1), Form("pT for bin %.2f-%.2f;pT", confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {axisPt[iPt]})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mPhi", iPt + 1), Form("#phi for bin %.2f-%.2f;#phi", confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{1000, 0, 2 * TMath::Pi()}})));
      mHistArrQA.push_back(std::get<std::shared_ptr<TH1>>(histos.add(Form("bin%i/mMultiplicity", iPt + 1), Form("Multiplicity for bin %.2f-%.2f;Multiplicity", confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH1F, {{1000, 0, 8000}})));
      for (auto iM = 0; iM < 40; ++iM) {
        auto mHistsR = std::get<std::shared_ptr<TH2>>(histos.add(Form("bin%i/Reset/mEtaPhi%i", iPt + 1, iM), Form("#eta#phi_%i for bin %.2f-%.2f;#eta;#phi", iM, confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH2F, {{binningM[iM], -0.8, 0.8}, {binningM[iM], 0, 2 * TMath::Pi()}}));
        mHistArrReset.push_back(mHistsR);
        auto mHistsUR = std::get<std::shared_ptr<TH2>>(histos.add(Form("bin%i/UReset/mEtaPhi%i", iPt + 1, iM), Form("#eta#phi_%i for bin %.2f-%.2f;#eta;#phi", iM, confPtBins.value[2 * iPt], confPtBins.value[2 * iPt + 1]), HistType::kTH2F, {{binningM[iM], -0.8, 0.8}, {binningM[iM], 0, 2 * TMath::Pi()}}));
        mHistArrUReset.push_back(mHistsUR);
      }
    }
  }

  template <class T>
  void checkpT(const T& track)
  {
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      if (track.pt() > confPtBins.value[2 * iPt] && track.pt() < confPtBins.value[2 * iPt + 1]) {
        mHistArrQA[iPt * 4]->Fill(track.eta());
        mHistArrQA[iPt * 4 + 1]->Fill(track.pt());
        mHistArrQA[iPt * 4 + 2]->Fill(track.phi());
        countTracks[iPt]++;
        for (auto iM = 0; iM < 40; ++iM) {
          mHistArrReset[iPt * 40 + iM]->Fill(track.eta(), track.phi());
          mHistArrUReset[iPt * 40 + iM]->Fill(track.eta(), track.phi());
        }
      }
    }
  }

  void calculateMoments(std::vector<std::shared_ptr<TH2>> hist)
  {
    Double_t binContent = 0;
    // Calculate the normalized factorial moments
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      for (auto iM = 0; iM < 40; ++iM) {
        binContent = 0;
        Double_t sumfqBin[6] = {0};

        for (auto iEta = 1; iEta <= hist[iPt * 40 + iM]->GetNbinsX(); ++iEta) {
          for (auto iPhi = 1; iPhi <= hist[iPt * 40 + iM]->GetNbinsY(); ++iPhi) {
            double binconVal = 0;
            binconVal = hist[iPt * 40 + iM]->GetBinContent(iEta, iPhi);
            binContent += binconVal;
            for (auto iOrder = 0; iOrder < 6; ++iOrder) {
              Double_t fqBin = 0;
              if (binconVal >= iOrder + 2) {
                fqBin = TMath::Factorial(binconVal) / (TMath::Factorial(binconVal - (iOrder + 2)));
              }
              if (isnan(fqBin)) {
                break;
              }
              sumfqBin[iOrder] += fqBin;
            }
          }
        }
        averagebinCon[iPt][iM] = binContent / (TMath::Power(binningM[iM], 2));
        valueM[iPt][iM] = binningM[iM];
        mPt[iPt][iM] = iPt;
        for (auto iOrder = 0; iOrder < 6; ++iOrder) {
          if (sumfqBin[iOrder] > 0) {
            fqEvent[iOrder][iPt][iM] = sumfqBin[iOrder] / (TMath::Power(binningM[iM], 2));
          }
        }
      } // end of loop over M bins
    }   // end of loop over pT bins
  }

  using TracksFMs = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>::iterator const& coll, TracksFMs const& tracks)
  {
    // selection of events
    if (!coll.sel8()) {
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
    fqEvent = {0, 0, 0, 0, 0, 0};
    averagebinCon = {0, 0, 0, 0, 0};
    valueM = {0, 0, 0, 0, 0};
    mPt = {0, 0, 0, 0, 0};

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
      }
    }
    // Calculate the normalized factorial moments
    calculateMoments(mHistArrReset);
    // Fill the table
    Int_t binMult[5] = {0};
    Double_t avBinconArr[200] = {0};
    Int_t valueMArr[200] = {0};
    Int_t mPtArr[200] = {0};
    Double_t fqEventArr2[200] = {0};
    Double_t fqEventArr3[200] = {0};
    Double_t fqEventArr4[200] = {0};
    Double_t fqEventArr5[200] = {0};
    Double_t fqEventArr6[200] = {0};
    Double_t fqEventArr7[200] = {0};
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      binMult[iPt] = countTracks[iPt];
      for (auto iM = 0; iM < 40; ++iM) {
        for (auto iOrder = 0; iOrder < 6; ++iOrder) {
          if (isinf(fqEvent[iOrder][iPt][iM])) {
            fqEvent[iOrder][iPt][iM] = 0;
          }
        }
        avBinconArr[iPt * 40 + iM] = averagebinCon[iPt][iM];
        valueMArr[iPt * 40 + iM] = valueM[iPt][iM];
        mPtArr[iPt * 40 + iM] = mPt[iPt][iM];
        fqEventArr2[iPt * 40 + iM] = fqEvent[0][iPt][iM];
        fqEventArr3[iPt * 40 + iM] = fqEvent[1][iPt][iM];
        fqEventArr4[iPt * 40 + iM] = fqEvent[2][iPt][iM];
        fqEventArr5[iPt * 40 + iM] = fqEvent[3][iPt][iM];
        fqEventArr6[iPt * 40 + iM] = fqEvent[4][iPt][iM];
        fqEventArr7[iPt * 40 + iM] = fqEvent[5][iPt][iM];
      }
    }
    tableFacMoments(coll.multFT0C(), coll.posZ(), coll.centFT0C(), binMult, avBinconArr, valueMArr, mPtArr, fqEventArr2, fqEventArr3, fqEventArr4, fqEventArr5, fqEventArr6, fqEventArr7);
  }
  PROCESS_SWITCH(FactorialMoments, processRun3, "main process function", false);

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
    fqEvent = {0, 0, 0, 0, 0, 0};
    averagebinCon = {0, 0, 0, 0, 0};
    valueM = {0, 0, 0, 0, 0};
    mPt = {0, 0, 0, 0, 0};

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
      }
    }
    // Calculate the normalized factorial moments
    calculateMoments(mHistArrReset);
    // Fill the table
    Int_t binMult[5] = {0};
    Double_t avBinconArr[200] = {0};
    Int_t valueMArr[200] = {0};
    Int_t mPtArr[200] = {0};
    Double_t fqEventArr2[200] = {0};
    Double_t fqEventArr3[200] = {0};
    Double_t fqEventArr4[200] = {0};
    Double_t fqEventArr5[200] = {0};
    Double_t fqEventArr6[200] = {0};
    Double_t fqEventArr7[200] = {0};
    for (auto iPt = 0; iPt < confNumPt; ++iPt) {
      binMult[iPt] = countTracks[iPt];
      for (auto iM = 0; iM < 40; ++iM) {
        for (auto iOrder = 0; iOrder < 6; ++iOrder) {
          if (isinf(fqEvent[iOrder][iPt][iM])) {
            fqEvent[iOrder][iPt][iM] = 0;
          }
        }
        avBinconArr[iPt * 40 + iM] = averagebinCon[iPt][iM];
        valueMArr[iPt * 40 + iM] = valueM[iPt][iM];
        mPtArr[iPt * 40 + iM] = mPt[iPt][iM];
        fqEventArr2[iPt * 40 + iM] = fqEvent[0][iPt][iM];
        fqEventArr3[iPt * 40 + iM] = fqEvent[1][iPt][iM];
        fqEventArr4[iPt * 40 + iM] = fqEvent[2][iPt][iM];
        fqEventArr5[iPt * 40 + iM] = fqEvent[3][iPt][iM];
        fqEventArr6[iPt * 40 + iM] = fqEvent[4][iPt][iM];
        fqEventArr7[iPt * 40 + iM] = fqEvent[5][iPt][iM];
      }
    }
    tableFacMoments(coll.centRun2V0M(), coll.posZ(), coll.centRun2V0M(), binMult, avBinconArr, valueMArr, mPtArr, fqEventArr2, fqEventArr3, fqEventArr4, fqEventArr5, fqEventArr6, fqEventArr7);
  }
  PROCESS_SWITCH(FactorialMoments, processRun2, "for RUN2", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FactorialMoments>(cfgc),
  };
}
