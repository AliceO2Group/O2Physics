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
//
/// \brief Hadron-nuclei correlation analysis task
/// \author Francesca Ercolessi
/// \since 21 April 2024

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <TParameter.h>
#include <TH1F.h>
#include <vector>
#include <TVector2.h>
#include <TVector3.h>

#include "Framework/ASoA.h"
#include "MathUtils/Utils.h"
#include "Framework/DataTypes.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"

#include "Framework/StaticFor.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"
#include "PWGCF/Femto3D/Core/femto3dPairTask.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct hadronnucleicorrelation {

  Configurable<bool> doQA{"doQA", true, "save QA histograms"};
  Configurable<bool> isMC{"isMC", false, "is MC"};
  Configurable<bool> disable_pantip{"disable_pantip", false, "disable_pantip"};

  // Event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0, "|vertexZ| value limit"};

  // Track selection
  Configurable<int16_t> min_TPC_nClusters{"min_TPC_nClusters", 80, "minimum number of found TPC clusters"};
  Configurable<float> min_TPC_nCrossedRowsOverFindableCls{"min_TPC_nCrossedRowsOverFindableCls", 0.8, "n TPC Crossed Rows Over Findable Cls"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4.0f, "maximum TPC chi^2/Ncls"};
  Configurable<float> max_chi2_ITS{"max_chi2_ITS", 36.0f, "maximum ITS chi^2/Ncls"};
  Configurable<float> etacut{"etacut", 0.8f, "eta cut"};
  Configurable<float> max_dcaxy{"max_dcaxy", 0.14f, "Maximum DCAxy"};
  Configurable<float> max_dcaz{"max_dcaz", 0.1f, "Maximum DCAz"};
  Configurable<float> nsigmaTPC{"nsigmaTPC", 3.0f, "cut nsigma TPC"};
  Configurable<float> nsigmaTOF{"nsigmaTOF", 3.5f, "cut nsigma TOF"};

  // Mixing parameters
  Configurable<int> _vertexNbinsToMix{"vertexNbinsToMix", 10, "Number of vertexZ bins for the mixing"};
  Configurable<int> _multNsubBins{"multSubBins", 10, "number of sub-bins to perform the mixing within"};

  // pT/A bins
  Configurable<std::vector<double>> pTA{"pTA", {0.4f, 0.6f, 0.8f}, "p_{T}/A bins"};
  ConfigurableAxis AxisNSigma{"AxisNSigma", {50, -10.f, 10.f}, "n#sigma"};

  using FilteredCollisions = aod::SingleCollSels;
  using FilteredTracks = aod::SingleTrackSels;
  using FilteredTracksMC = soa::Join<aod::SingleTrackSels, aod::SingleTrkMCs>;

  HistogramRegistry registry{"registry"};
  HistogramRegistry QA{"QA"};

  typedef std::shared_ptr<soa::Filtered<FilteredTracks>::iterator> trkType;
  typedef std::shared_ptr<soa::Filtered<FilteredCollisions>::iterator> colType;

  // key: int64_t - value: vector of trkType objects
  std::map<int64_t, std::vector<trkType>> selectedtracks_p;
  std::map<int64_t, std::vector<trkType>> selectedtracks_antid;
  std::map<int64_t, std::vector<trkType>> selectedtracks_antip;

  // key: pair of an integer and a float - value: vector of colType objects
  // for each key I have a vector of collisions
  std::map<std::pair<int, float>, std::vector<colType>> mixbins_antidantip;
  std::map<std::pair<int, float>, std::vector<colType>> mixbins_pantip;

  std::vector<std::shared_ptr<TH3>> hEtaPhi_PrAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhi_PrAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hEtaPhi_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhi_AntiDeAntiPr_ME;

  int nBins;

  void init(o2::framework::InitContext&)
  {
    AxisSpec ptAxis = {pTA, "#it{p}_{T}/A of #bar{p} (GeV/c)"};
    AxisSpec etaAxis = {100, -1.5, 1.5, "#Delta#eta"};
    AxisSpec phiAxis = {157, -TMath::Pi() / 2, 1.5 * TMath::Pi(), "#Delta#phi"};
    AxisSpec pTAxis = {200, -10.f, 10.f, "p_{T} GeV/c"};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1I, {{3, 0.f, 3.f}}});
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "Selected");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "#bar{d}-#bar{p}");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "p-#bar{p}");

    if (!isMC) {
      registry.add("hDebug", "hDebug", {HistType::kTH1I, {{4, 0.f, 4.f}}});
      registry.get<TH1>(HIST("hDebug"))->GetXaxis()->SetBinLabel(1, "all");
      registry.get<TH1>(HIST("hDebug"))->GetXaxis()->SetBinLabel(2, "ev. with #bar{d}");
      registry.get<TH1>(HIST("hDebug"))->GetXaxis()->SetBinLabel(3, "ev. with #bar{p}");
      registry.get<TH1>(HIST("hDebug"))->GetXaxis()->SetBinLabel(4, "ev. with p");

      registry.add("hDebugdp", "hDebugdp", {HistType::kTH1I, {{6, 0.f, 6.f}}});
      registry.get<TH1>(HIST("hDebugdp"))->GetXaxis()->SetBinLabel(1, "N coll with #bar{d}");
      registry.get<TH1>(HIST("hDebugdp"))->GetXaxis()->SetBinLabel(2, "N mixing bins");
      registry.get<TH1>(HIST("hDebugdp"))->GetXaxis()->SetBinLabel(3, "N coll with #bar{d}");
      registry.get<TH1>(HIST("hDebugdp"))->GetXaxis()->SetBinLabel(4, "#bar{d}-#bar{p} pairs SE");
      registry.get<TH1>(HIST("hDebugdp"))->GetXaxis()->SetBinLabel(5, "#bar{d}-#bar{p} pairs ME");

      nBins = pTA.value.size() - 1;
      for (int i = 0; i < nBins; i++) {
        if (!disable_pantip) {
          auto htempSE_PrAntiPr = registry.add<TH3>(Form("hEtaPhi_PrAntiPr_SE_ptA%02.0f%02.0f", pTA.value.at(i) * 10, pTA.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T}/A pr <%.1f)", pTA.value.at(i), pTA.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptAxis}});
          auto htempME_PrAntiPr = registry.add<TH3>(Form("hEtaPhi_PrAntiPr_ME_ptA%02.0f%02.0f", pTA.value.at(i) * 10, pTA.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T}/A pr <%.1f)", pTA.value.at(i), pTA.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptAxis}});

          hEtaPhi_PrAntiPr_SE.push_back(std::move(htempSE_PrAntiPr));
          hEtaPhi_PrAntiPr_ME.push_back(std::move(htempME_PrAntiPr));
        }

        auto htempSE_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_AntiDeAntiPr_SE_ptA%02.0f%02.0f", pTA.value.at(i) * 10, pTA.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T}/A #bar{d} <%.1f)", pTA.value.at(i), pTA.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptAxis}});
        auto htempME_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_AntiDeAntiPr_ME_ptA%02.0f%02.0f", pTA.value.at(i) * 10, pTA.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T}/A #bar{d} <%.1f)", pTA.value.at(i), pTA.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptAxis}});
        hEtaPhi_AntiDeAntiPr_SE.push_back(std::move(htempSE_AntiDeAntiPr));
        hEtaPhi_AntiDeAntiPr_ME.push_back(std::move(htempME_AntiDeAntiPr));
      }
    }

    if (doQA) {
      // Track QA
      QA.add("QA/hVtxZ_trk", "#it{z}_{vtx}", {HistType::kTH1F, {{200, -20.f, 20.f, "#it{z}_{vtx} (cm)"}}});
      QA.add("QA/hTPCnClusters", "N TPC Clusters; N TPC Clusters", {HistType::kTH1F, {{200, 0.f, 200.f}}});
      QA.add("QA/hTPCchi2", "TPC chi2/Ncls; TPC chi2/Ncls", {HistType::kTH1F, {{200, 0.f, 20.f}}});
      QA.add("QA/hTPCcrossedRowsOverFindableCls", "TPC crossed Rows Over Findable Cls; TPC Crossed Rows Over Findable Cls", {HistType::kTH1F, {{100, 0.f, 2.f}}});
      QA.add("QA/hITSchi2", "ITS chi2/Ncls; ITS chi2/Ncls", {HistType::kTH1F, {{200, 0.f, 40.f}}});
      QA.add("QA/hDCAxy", "DCAxy", {HistType::kTH1F, {{80, -0.4f, 0.4f, "DCA xy (cm)"}}});
      QA.add("QA/hDCAz", "DCAz", {HistType::kTH1F, {{80, -0.4f, 0.4f, "DCA z (cm)"}}});
      QA.add("QA/hnSigmaTPCVsPt_Pr", "n#sigma TPC vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      QA.add("QA/hnSigmaTPCVsPt_De", "n#sigma TPC vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      QA.add("QA/hnSigmaTOFVsPt_Pr", "n#sigma TOF vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      QA.add("QA/hnSigmaTOFVsPt_De", "n#sigma TOF vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});

      if (!isMC) {
        QA.add("QA/hEtaPr", Form("#eta ditribution for p"), {HistType::kTH1F, {{200, -1.f, 1.f, "#eta"}}});
        QA.add("QA/hPhiPr", Form("#phi ditribution for p"), {HistType::kTH1F, {{100, 0.f, 2 * TMath::Pi(), "#phi"}}});
        QA.add("QA/hEtaAntiPr", Form("#eta ditribution for #bar{p}"), {HistType::kTH1F, {{200, -1.f, 1.f, "#eta"}}});
        QA.add("QA/hPhiAntiPr", Form("#phi ditribution for #bar{p}"), {HistType::kTH1F, {{100, 0.f, 2 * TMath::Pi(), "#phi"}}});
        QA.add("QA/hEtaDe", Form("#eta ditribution for d"), {HistType::kTH1F, {{200, -1.f, 1.f, "#eta"}}});
        QA.add("QA/hPhiDe", Form("#phi ditribution for d"), {HistType::kTH1F, {{100, 0.f, 2 * TMath::Pi(), "#phi"}}});
        QA.add("QA/hEtaAntiDe", Form("#eta ditribution for #bar{d}"), {HistType::kTH1F, {{200, -1.f, 1.f, "#eta"}}});
        QA.add("QA/hPhiAntiDe", Form("#phi ditribution for #bar{d}"), {HistType::kTH1F, {{100, 0.f, 2 * TMath::Pi(), "#phi"}}});

        QA.add("QA/hnSigmaTPCVsPt_Pr_AfterSel", "n#sigma TPC vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsPt_De_AfterSel", "n#sigma TPC vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});

        QA.add("QA/hnSigmaTPCVsPhi_Pr", Form("n#sigma TPC vs #phi p; #phi; n#sigma TPC"), {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi()}, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsPhi_De", Form("n#sigma TPC vs #phi d; #phi; n#sigma TPC"), {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi()}, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsPhi_AntiPr", Form("n#sigma TPC vs #phi #bar{p}; #phi; n#sigma TPC"), {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi()}, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsPhi_AntiDe", Form("n#sigma TPC vs #phi #bar{d}; #phi; n#sigma TPC"), {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi()}, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsEta_Pr", Form("n#sigma TPC vs #eta p; #eta; n#sigma TPC"), {HistType::kTH2F, {{100, -1.f, +1.f}, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsEta_De", Form("n#sigma TPC vs #eta d; #eta; n#sigma TPC"), {HistType::kTH2F, {{100, -1.f, +1.f}, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsEta_AntiPr", Form("n#sigma TPC vs #eta #bar{p}; #eta; n#sigma TPC"), {HistType::kTH2F, {{100, -1.f, +1.f}, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsEta_AntiDe", Form("n#sigma TPC vs #eta #bar{d}; #eta; n#sigma TPC"), {HistType::kTH2F, {{100, -1.f, +1.f}, AxisNSigma}});

        QA.add("QA/hnSigmaTOFVsPt_Pr_AfterSel", "n#sigma TOF vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPt_De_AfterSel", "n#sigma TOF vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPhi_De", Form("n#sigma TOF vs #phi d; #phi; n#sigma TOF"), {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi()}, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPhi_AntiDe", Form("n#sigma TOF vs #phi #bar{d}; #phi; n#sigma TOF"), {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi()}, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsEta_De", Form("n#sigma TOF vs #eta d; #eta; n#sigma TOF"), {HistType::kTH2F, {{100, -1.f, +1.f}, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsEta_AntiDe", Form("n#sigma TOF vs #eta #bar{d}; #eta; n#sigma TOF"), {HistType::kTH2F, {{100, -1.f, +1.f}, AxisNSigma}});
      }
    }

    if (isMC) {
      registry.add("hReco_EtaPhiPt_Proton", "Gen (anti)protons in reco collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., 2 * TMath::Pi(), "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hReco_EtaPhiPt_Deuteron", "Gen (anti)deuteron in reco collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., 2 * TMath::Pi(), "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hReco_PID_EtaPhiPt_Proton", "Gen (anti)protons + PID in reco collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., 2 * TMath::Pi(), "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hReco_PID_EtaPhiPt_Deuteron", "Gen (anti)deuteron + PID in reco collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., 2 * TMath::Pi(), "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hReco_EtaPhiPtMC_Proton", "Gen (anti)protons in reco collisions (MC info used)", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., 2 * TMath::Pi(), "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hReco_EtaPhiPtMC_Deuteron", "Gen (anti)deuteron in reco collisions (MC info used)", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., 2 * TMath::Pi(), "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});

      registry.add("hnSigmaTPCVsPt_Pr_MC", "n#sigma TPC vs p_{T} for p hypothesis true MC; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      registry.add("hnSigmaTPCVsPt_De_MC", "n#sigma TPC vs p_{T} for d hypothesis true MC; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      registry.add("hnSigmaTOFVsPt_Pr_MC", "n#sigma TOF vs p_{T} for p hypothesis true MC; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      registry.add("hnSigmaTOFVsPt_De_MC", "n#sigma TOF vs p_{T} for d hypothesis true MC; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
    }
  }

  // Filters
  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) <= cutzvertex;
  Filter trackFilter = o2::aod::singletrackselector::tpcNClsFound >= min_TPC_nClusters &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedTpcChi2NCl) <= max_chi2_TPC &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::rowsOverFindable>(o2::aod::singletrackselector::storedTpcCrossedRowsOverFindableCls) >= min_TPC_nCrossedRowsOverFindableCls &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedItsChi2NCl) <= max_chi2_ITS &&
                       nabs(o2::aod::singletrackselector::unPack<singletrackselector::binning::dca>(o2::aod::singletrackselector::storedDcaXY)) <= max_dcaxy &&
                       nabs(o2::aod::singletrackselector::unPack<singletrackselector::binning::dca>(o2::aod::singletrackselector::storedDcaZ)) <= max_dcaz &&
                       nabs(o2::aod::singletrackselector::eta) <= etacut;

  template <int ME, typename Type>
  void mixTracks(Type const& tracks1, Type const& tracks2, bool isDe)
  { // last value: 0 -- SE; 1 -- ME
    for (auto it1 : tracks1) {
      for (auto it2 : tracks2) {

        // Variables
        float deltaEta = it2->eta() - it1->eta();
        float deltaPhi = it2->phi() - it1->phi();

        while (deltaPhi < -TMath::Pi() / 2) {
          deltaPhi += 2 * TMath::Pi();
        }

        while (deltaPhi >= 3 * TMath::Pi() / 2) {
          deltaPhi -= 2 * TMath::Pi();
        }

        for (int k = 0; k < nBins; k++) {
          if (!isDe && !disable_pantip) {
            if (it1->pt() > pTA.value.at(k) && it1->pt() <= pTA.value.at(k + 1)) {
              if (ME) {
                hEtaPhi_PrAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
              } else {
                hEtaPhi_PrAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
              }
            }
          } else {
            if (it1->pt() > pTA.value.at(k) * 2 && it1->pt() <= pTA.value.at(k + 1) * 2) {
              if (ME) {
                hEtaPhi_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
              } else {
                hEtaPhi_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
              }
            }
          }
        }
      }
    }
  }

  void processData(soa::Filtered<FilteredCollisions> const& collisions, soa::Filtered<FilteredTracks> const& tracks)
  {
    for (auto track : tracks) {
      if (abs(track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().posZ()) > cutzvertex)
        continue;

      if (doQA) {
        QA.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        QA.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        QA.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        QA.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        QA.fill(HIST("QA/hDCAxy"), track.dcaXY());
        QA.fill(HIST("QA/hDCAz"), track.dcaZ());
        QA.fill(HIST("QA/hVtxZ_trk"), track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().posZ());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
      }

      bool isPr = false;
      bool isAntiPr = false;
      bool isDeTPCTOF = false;
      bool isAntiDeTPCTOF = false;

      if (TMath::Abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() > 0)
        isPr = true;
      if (TMath::Abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() < 0)
        isAntiPr = true;
      if (TMath::Abs(track.tpcNSigmaDe()) < nsigmaTPC && TMath::Abs(track.tofNSigmaDe()) < nsigmaTOF && track.sign() > 0)
        isDeTPCTOF = true;
      if (TMath::Abs(track.tpcNSigmaDe()) < nsigmaTPC && TMath::Abs(track.tofNSigmaDe()) < nsigmaTOF && track.sign() < 0)
        isAntiDeTPCTOF = true;

      // Deuterons
      if (isAntiDeTPCTOF) {
        selectedtracks_antid[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        if (doQA) {
          QA.fill(HIST("QA/hEtaAntiDe"), track.eta());
          QA.fill(HIST("QA/hPhiAntiDe"), track.phi());
          QA.fill(HIST("QA/hnSigmaTPCVsPhi_AntiDe"), track.phi(), track.tpcNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTPCVsEta_AntiDe"), track.eta(), track.tpcNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTOFVsPhi_AntiDe"), track.phi(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTOFVsEta_AntiDe"), track.eta(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
        }
      }
      if (isDeTPCTOF) {
        if (doQA) {
          QA.fill(HIST("QA/hEtaDe"), track.eta());
          QA.fill(HIST("QA/hPhiDe"), track.phi());
          QA.fill(HIST("QA/hnSigmaTPCVsPhi_De"), track.phi(), track.tpcNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTPCVsEta_De"), track.eta(), track.tpcNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTOFVsPhi_De"), track.phi(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTOFVsEta_De"), track.eta(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
        }
      }

      // Protons
      if (isPr) {
        selectedtracks_p[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        if (doQA) {
          QA.fill(HIST("QA/hEtaPr"), track.eta());
          QA.fill(HIST("QA/hPhiPr"), track.phi());
          QA.fill(HIST("QA/hnSigmaTPCVsPhi_Pr"), track.phi(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsEta_Pr"), track.eta(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
        }
      } else if (isAntiPr) {
        selectedtracks_antip[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        if (doQA) {
          QA.fill(HIST("QA/hEtaAntiPr"), track.eta());
          QA.fill(HIST("QA/hPhiAntiPr"), track.phi());
          QA.fill(HIST("QA/hnSigmaTPCVsPhi_AntiPr"), track.phi(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsEta_AntiPr"), track.eta(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
        }
      }
    }

    for (auto collision : collisions) {

      if (TMath::Abs(collision.posZ()) > cutzvertex)
        continue;

      registry.fill(HIST("hNEvents"), 0.5);
      registry.fill(HIST("hDebug"), 0.5);

      if (selectedtracks_p.find(collision.globalIndex()) != selectedtracks_p.end() &&
          selectedtracks_antip.find(collision.globalIndex()) != selectedtracks_antip.end()) {
        registry.fill(HIST("hNEvents"), 2.5);
      }

      if (selectedtracks_antid.find(collision.globalIndex()) != selectedtracks_antid.end() &&
          selectedtracks_antip.find(collision.globalIndex()) != selectedtracks_antip.end()) {
        registry.fill(HIST("hNEvents"), 1.5);
      }

      int vertexBinToMix = std::floor((collision.posZ() + cutzvertex) / (2 * cutzvertex / _vertexNbinsToMix));
      int centBinToMix = std::floor(collision.multPerc() / (100.0 / _multNsubBins));

      if (selectedtracks_p.find(collision.globalIndex()) != selectedtracks_p.end()) {
        registry.fill(HIST("hDebug"), 3.5);
        mixbins_pantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }

      if (selectedtracks_antip.find(collision.globalIndex()) != selectedtracks_antip.end()) {
        registry.fill(HIST("hDebug"), 2.5);
      }

      if (selectedtracks_antid.find(collision.globalIndex()) != selectedtracks_antid.end()) {
        registry.fill(HIST("hDebug"), 1.5);
        registry.fill(HIST("hDebugdp"), 0.5); // numero tot di collisioni nella mappa mixbins_antidantip
        mixbins_antidantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }
    }

    registry.get<TH1>(HIST("hDebugdp"))->SetBinContent(6, mixbins_antidantip.size());

    if (!disable_pantip) {
      if (!mixbins_pantip.empty()) {
        for (auto i = mixbins_pantip.begin(); i != mixbins_pantip.end(); i++) { // iterating over all vertex&mult bins

          int EvPerBin = (i->second).size(); // number of collisions in each vertex&mult bin

          for (int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

            auto col1 = (i->second)[indx1];

            if (selectedtracks_antip.find(col1->index()) != selectedtracks_antip.end()) {
              mixTracks<0>(selectedtracks_p[col1->index()], selectedtracks_antip[col1->index()], 0); // mixing SE
            }

            int indx3 = EvPerBin;
            if (indx1 < (EvPerBin - 11)) {
              indx3 = indx1 + 11;
            }

            for (int indx2 = indx1 + 1; indx2 < indx3; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

              auto col2 = (i->second)[indx2];

              if (col1 == col2) {
                continue;
              }

              if (selectedtracks_antip.find(col2->index()) != selectedtracks_antip.end()) {
                mixTracks<1>(selectedtracks_p[col1->index()], selectedtracks_antip[col2->index()], 0); // mixing ME
              }
            }
          }
        }
      }
    }

    if (!mixbins_antidantip.empty()) {

      for (auto i = mixbins_antidantip.begin(); i != mixbins_antidantip.end(); i++) { // iterating over all vertex&mult bins

        registry.fill(HIST("hDebugdp"), 1.5); // numero di keys (vertex&mult bins) nella mappa mixbins_antidantip

        std::vector<colType> value = i->second;
        int EvPerBin = value.size(); // number of collisions in each vertex&mult bin

        for (int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          registry.fill(HIST("hDebugdp"), 2.5);

          auto col1 = value[indx1];

          if (selectedtracks_antip.find(col1->index()) != selectedtracks_antip.end()) {
            registry.fill(HIST("hDebugdp"), 3.5);
            mixTracks<0>(selectedtracks_antid[col1->index()], selectedtracks_antip[col1->index()], 1); // mixing SE
          }

          int indx3 = EvPerBin;
          if (indx1 < (EvPerBin - 11)) {
            indx3 = indx1 + 11;
          }

          for (int indx2 = 0; indx2 < indx3; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = value[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracks_antip.find(col2->index()) != selectedtracks_antip.end()) {
              registry.fill(HIST("hDebugdp"), 4.5);
              mixTracks<1>(selectedtracks_antid[col1->index()], selectedtracks_antip[col2->index()], 1); // mixing ME
            }
          }
        }
      }
    }

    // clearing up
    for (auto i = selectedtracks_antid.begin(); i != selectedtracks_antid.end(); i++)
      (i->second).clear();
    selectedtracks_antid.clear();

    for (auto i = selectedtracks_antip.begin(); i != selectedtracks_antip.end(); i++)
      (i->second).clear();
    selectedtracks_antip.clear();

    for (auto i = selectedtracks_p.begin(); i != selectedtracks_p.end(); i++)
      (i->second).clear();
    selectedtracks_p.clear();

    for (auto i = mixbins_antidantip.begin(); i != mixbins_antidantip.end(); i++)
      (i->second).clear();
    mixbins_antidantip.clear();

    for (auto i = mixbins_pantip.begin(); i != mixbins_pantip.end(); i++)
      (i->second).clear();
    mixbins_pantip.clear();
  }
  PROCESS_SWITCH(hadronnucleicorrelation, processData, "processData", true);

  void processMC(soa::Filtered<FilteredCollisions>::iterator const& collision, soa::Filtered<FilteredTracksMC> const& tracks)
  {
    if (TMath::Abs(collision.posZ()) > cutzvertex)
      return;

    registry.fill(HIST("hNEvents"), 0.5);

    for (auto track : tracks) {

      if (doQA) {
        QA.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        QA.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        QA.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        QA.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        QA.fill(HIST("QA/hDCAxy"), track.dcaXY());
        QA.fill(HIST("QA/hDCAz"), track.dcaZ());
        QA.fill(HIST("QA/hVtxZ_trk"), track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().posZ());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
      }

      if (track.origin() != 0)
        continue;

      if (abs(track.pdgCode()) != 2212 && abs(track.pdgCode()) != 1000010020)
        continue;

      if (track.pdgCode() == 2212) {
        registry.fill(HIST("hReco_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt());
        registry.fill(HIST("hReco_EtaPhiPtMC_Proton"), track.eta_MC(), track.phi_MC(), track.pt_MC());

        if (TMath::Abs(track.tpcNSigmaPr()) < nsigmaTPC) {
          registry.fill(HIST("hReco_PID_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt());
        }
        registry.fill(HIST("hnSigmaTPCVsPt_Pr_MC"), track.pt(), track.tpcNSigmaPr());
        registry.fill(HIST("hnSigmaTOFVsPt_Pr_MC"), track.pt(), track.tofNSigmaPr());
      }
      if (track.pdgCode() == -2212) {
        registry.fill(HIST("hReco_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        registry.fill(HIST("hReco_EtaPhiPtMC_Proton"), track.eta_MC(), track.phi_MC(), track.pt_MC() * -1);

        if (TMath::Abs(track.tpcNSigmaPr()) < nsigmaTPC) {
          registry.fill(HIST("hReco_PID_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        }
        registry.fill(HIST("hnSigmaTPCVsPt_Pr_MC"), track.pt() * -1, track.tpcNSigmaPr());
        registry.fill(HIST("hnSigmaTOFVsPt_Pr_MC"), track.pt() * -1, track.tofNSigmaPr());
      }
      if (track.pdgCode() == 1000010020) {
        registry.fill(HIST("hReco_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt());
        registry.fill(HIST("hReco_EtaPhiPtMC_Deuteron"), track.eta_MC(), track.phi_MC(), track.pt_MC());

        if (TMath::Abs(track.tpcNSigmaDe()) < nsigmaTPC && TMath::Abs(track.tofNSigmaDe()) < nsigmaTOF) {
          registry.fill(HIST("hReco_PID_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt());
        }
        registry.fill(HIST("hnSigmaTPCVsPt_De_MC"), track.pt(), track.tpcNSigmaDe());
        registry.fill(HIST("hnSigmaTOFVsPt_De_MC"), track.pt(), track.tofNSigmaDe());
      }
      if (track.pdgCode() == -1000010020) {
        registry.fill(HIST("hReco_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt() * -1);
        registry.fill(HIST("hReco_EtaPhiPtMC_Deuteron"), track.eta_MC(), track.phi_MC(), track.pt_MC() * -1);

        if (TMath::Abs(track.tpcNSigmaDe()) < nsigmaTPC && TMath::Abs(track.tofNSigmaDe()) < nsigmaTOF) {
          registry.fill(HIST("hReco_PID_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt() * -1);
        }
        registry.fill(HIST("hnSigmaTPCVsPt_De_MC"), track.pt() * -1, track.tpcNSigmaDe());
        registry.fill(HIST("hnSigmaTOFVsPt_De_MC"), track.pt() * -1, track.tofNSigmaDe());
      }
    }
  }
  PROCESS_SWITCH(hadronnucleicorrelation, processMC, "processMC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<hadronnucleicorrelation>(cfgc)};
}