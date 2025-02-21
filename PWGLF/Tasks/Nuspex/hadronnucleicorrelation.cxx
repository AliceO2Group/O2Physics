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

#include <TParameter.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <utility>
#include <TVector2.h>
#include <TVector3.h>
#include "TGrid.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

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

  // PDG codes and masses used in this analysis
  static constexpr int pdgProton = 2212;
  static constexpr int pdgDeuteron = 1000010020;

  Configurable<bool> doQA{"doQA", true, "save QA histograms"};
  Configurable<bool> doMCQA{"doMCQA", false, "save MC QA histograms"};
  Configurable<bool> isMC{"isMC", false, "is MC"};
  Configurable<bool> isMCGen{"isMCGen", false, "is isMCGen"};
  Configurable<bool> isPrim{"isPrim", true, "is isPrim"};
  Configurable<bool> domatterGen{"domatterGen", true, "domatterGen"};
  Configurable<bool> mcCorrelation{"mcCorrelation", false, "true: build the correlation function only for SE"};
  Configurable<bool> docorrection{"docorrection", false, "do efficiency correction"};

  Configurable<std::string> fCorrectionPath{"fCorrectionPath", "", "Correction path to file"};
  Configurable<std::string> fCorrectionHisto{"fCorrectionHisto", "", "Correction histogram"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0, "|vertexZ| value limit"};

  // Track selection
  Configurable<double> par0{"par0", 0.004, "par 0"};
  Configurable<double> par1{"par1", 0.013, "par 1"};
  Configurable<int16_t> min_TPC_nClusters{"min_TPC_nClusters", 80, "minimum number of found TPC clusters"};
  Configurable<float> min_TPC_nCrossedRowsOverFindableCls{"min_TPC_nCrossedRowsOverFindableCls", 0.8, "n TPC Crossed Rows Over Findable Cls"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4.0f, "maximum TPC chi^2/Ncls"};
  Configurable<float> max_chi2_ITS{"max_chi2_ITS", 36.0f, "maximum ITS chi^2/Ncls"};
  Configurable<float> etacut{"etacut", 0.8f, "eta cut"};
  Configurable<float> max_dcaxy{"max_dcaxy", 0.14f, "Maximum DCAxy"};
  Configurable<float> max_dcaz{"max_dcaz", 0.1f, "Maximum DCAz"};
  Configurable<float> nsigmaTPC{"nsigmaTPC", 3.0f, "cut nsigma TPC"};
  Configurable<float> nsigmaElPr{"nsigmaElPr", 1.0f, "cut nsigma TPC El for protons"};
  Configurable<float> nsigmaElDe{"nsigmaElDe", 3.0f, "cut nsigma TPC El for protons"};
  Configurable<float> nsigmaTOF{"nsigmaTOF", 3.5f, "cut nsigma TOF"};
  Configurable<float> pTthrpr_TOF{"pTthrpr_TOF", 0.8f, "threshold pT proton to use TOF"};
  Configurable<float> pTthrpr_TPCEl{"pTthrpr_TPCEl", 1.0f, "threshold pT proton to use TPC El rejection"};
  Configurable<float> pTthrde_TOF{"pTthrde_TOF", 1.0f, "threshold pT deuteron to use TOF"};
  Configurable<float> pTthrde_TPCEl{"pTthrde_TPCEl", 1.0f, "threshold pT deuteron to use TPC El rejection"};
  Configurable<bool> rejectionEl{"rejectionEl", true, "use TPC El rejection"};
  Configurable<float> max_tpcSharedCls{"max_tpcSharedCls", 0.4, "maximum fraction of TPC shared clasters"};
  Configurable<int> min_itsNCls{"min_itsNCls", 0, "minimum allowed number of ITS clasters"};
  Configurable<int> maxmixcollsGen{"maxmixcollsGen", 100, "maxmixcollsGen"};

  // Mixing parameters
  Configurable<int> _vertexNbinsToMix{"vertexNbinsToMix", 10, "Number of vertexZ bins for the mixing"};
  Configurable<int> _multNsubBins{"multSubBins", 10, "number of sub-bins to perform the mixing within"};

  // pT/A bins
  Configurable<std::vector<double>> pTBins{"pTBins", {0.6f, 1.0f, 1.2f, 2.f}, "p_{T} bins"};

  ConfigurableAxis AxisNSigma{"AxisNSigma", {35, -7.f, 7.f}, "n#sigma"};

  using FilteredCollisions = soa::Filtered<aod::SingleCollSels>;
  using SimCollisions = aod::McCollisions;
  using SimParticles = aod::McParticles;
  using FilteredTracks = soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkExtras, aod::SinglePIDEls, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDHes>>;                      // new tables (v3)
  using FilteredTracksMC = soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkMCs, aod::SingleTrkExtras, aod::SinglePIDEls, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDHes>>; // new tables (v3)

  HistogramRegistry registry{"registry"};
  HistogramRegistry QA{"QA"};

  typedef std::shared_ptr<FilteredTracks::iterator> trkType;
  typedef std::shared_ptr<FilteredTracksMC::iterator> trkTypeMC;
  typedef std::shared_ptr<SimParticles::iterator> partTypeMC;
  typedef std::shared_ptr<FilteredCollisions::iterator> colType;
  typedef std::shared_ptr<SimCollisions::iterator> MCcolType;

  // key: int64_t - value: vector of trkType objects
  std::map<int64_t, std::vector<trkType>> selectedtracks_p;
  std::map<int64_t, std::vector<trkType>> selectedtracks_d;
  std::map<int64_t, std::vector<trkType>> selectedtracks_antid;
  std::map<int64_t, std::vector<trkType>> selectedtracks_antip;

  // key: int64_t - value: vector of trkType objects
  std::map<int64_t, std::vector<partTypeMC>> selectedparticlesMC_d;
  std::map<int64_t, std::vector<partTypeMC>> selectedparticlesMC_p;
  std::map<int64_t, std::vector<partTypeMC>> selectedparticlesMC_antid;
  std::map<int64_t, std::vector<partTypeMC>> selectedparticlesMC_antip;

  // key: int64_t - value: vector of trkType objects
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksMC_d;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksMC_p;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksMC_antid;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksMC_antip;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksPIDMC_d;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksPIDMC_p;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksPIDMC_antid;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksPIDMC_antip;

  // key: pair of an integer and a float - value: vector of colType objects
  // for each key I have a vector of collisions
  std::map<std::pair<int, float>, std::vector<colType>> mixbins_antidantip;
  std::map<std::pair<int, float>, std::vector<colType>> mixbinsPID_antidantip;
  std::map<float, std::vector<MCcolType>> mixbinsMC_antidantip;
  std::map<float, std::vector<MCcolType>> mixbinsMC_dp;

  std::vector<std::shared_ptr<TH3>> hEtaPhi_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhi_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhi_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhi_AntiDeAntiPr_ME;

  std::vector<std::shared_ptr<TH3>> hEtaPhiRec_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhiGen_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhiRec_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hEtaPhiGen_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiRec_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiGen_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiRec_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiGen_AntiDeAntiPr_ME;

  int nBinspT;
  TH2F* hEffpTEta_proton;
  TH2F* hEffpTEta_antiproton;
  TH2F* hEffpTEta_deuteron;
  TH2F* hEffpTEta_antideuteron;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    if (docorrection) {
      GetCorrection(ccdb, TString(fCorrectionPath), TString(fCorrectionHisto));
    } else {
      hEffpTEta_proton = nullptr;
      hEffpTEta_antiproton = nullptr;
      hEffpTEta_deuteron = nullptr;
      hEffpTEta_antideuteron = nullptr;
    }

    AxisSpec ptBinnedAxis = {pTBins, "#it{p}_{T} of #bar{p} (GeV/c)"};
    AxisSpec etaAxis = {100, -1., 1., "#eta"};
    AxisSpec phiAxis = {157, 0., 2 * TMath::Pi(), "#phi (rad)"};
    AxisSpec pTAxis = {200, -10.f, 10.f, "p_{T} GeV/c"};
    AxisSpec pTAxis_small = {100, -5.f, 5.f, "p_{T} GeV/c"};

    AxisSpec DeltaEtaAxis = {100, -1.5, 1.5, "#Delta#eta"};
    AxisSpec DeltaPhiAxis = {60, -TMath::Pi() / 2, 1.5 * TMath::Pi(), "#Delta#phi (rad)"};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1D, {{5, 0.f, 5.f}}});
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "Selected");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "events with #bar{d}-#bar{p}");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "events with d-p");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "events with #bar{d}");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "events with d");

    nBinspT = pTBins.value.size() - 1;

    if (mcCorrelation) {
      for (int i = 0; i < nBinspT; i++) {
        auto htempSERec_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhiRec_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                         Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        hEtaPhiRec_AntiDeAntiPr_SE.push_back(std::move(htempSERec_AntiDeAntiPr));
        auto htempMERec_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhiRec_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                         Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        hEtaPhiRec_AntiDeAntiPr_ME.push_back(std::move(htempMERec_AntiDeAntiPr));

        auto htempPIDSERec_AntiDeAntiPr = registry.add<TH3>(Form("hPIDEtaPhiRec_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                            Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        auto htempPIDSEGen_AntiDeAntiPr = registry.add<TH3>(Form("hPIDEtaPhiGen_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                            Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        hPIDEtaPhiRec_AntiDeAntiPr_SE.push_back(std::move(htempPIDSERec_AntiDeAntiPr));
        hPIDEtaPhiGen_AntiDeAntiPr_SE.push_back(std::move(htempPIDSEGen_AntiDeAntiPr));
        auto htempPIDMERec_AntiDeAntiPr = registry.add<TH3>(Form("hPIDEtaPhiRec_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                            Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        auto htempPIDMEGen_AntiDeAntiPr = registry.add<TH3>(Form("hPIDEtaPhiGen_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                            Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        hPIDEtaPhiRec_AntiDeAntiPr_ME.push_back(std::move(htempPIDMERec_AntiDeAntiPr));
        hPIDEtaPhiGen_AntiDeAntiPr_ME.push_back(std::move(htempPIDMEGen_AntiDeAntiPr));
      }
    }

    if (isMCGen || mcCorrelation) {
      for (int i = 0; i < nBinspT; i++) {
        auto htempSEGen_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhiGen_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                         Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        hEtaPhiGen_AntiDeAntiPr_SE.push_back(std::move(htempSEGen_AntiDeAntiPr));
        auto htempMEGen_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhiGen_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                         Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        hEtaPhiGen_AntiDeAntiPr_ME.push_back(std::move(htempMEGen_AntiDeAntiPr));
      }
    }

    if (!isMC) {
      for (int i = 0; i < nBinspT; i++) {
        auto htempSE_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        auto htempME_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        hEtaPhi_AntiDeAntiPr_SE.push_back(std::move(htempSE_AntiDeAntiPr));
        hEtaPhi_AntiDeAntiPr_ME.push_back(std::move(htempME_AntiDeAntiPr));

        auto hCorrtempSE_AntiDeAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        auto hCorrtempME_AntiDeAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
        hCorrEtaPhi_AntiDeAntiPr_SE.push_back(std::move(hCorrtempSE_AntiDeAntiPr));
        hCorrEtaPhi_AntiDeAntiPr_ME.push_back(std::move(hCorrtempME_AntiDeAntiPr));
      }
    }

    if (doQA) {
      // Track QA
      QA.add("QA/hVtxZ_trk", "#it{z}_{vtx}", {HistType::kTH1D, {{150, -15.f, 15.f, "#it{z}_{vtx} (cm)"}}});
      QA.add("QA/hTPCnClusters", "N TPC Clusters; N TPC Clusters", {HistType::kTH1D, {{200, 0.f, 200.f}}});
      QA.add("QA/hTPCchi2", "TPC chi2/Ncls; TPC chi2/Ncls", {HistType::kTH1D, {{100, 0.f, 10.f}}});
      QA.add("QA/hTPCcrossedRowsOverFindableCls", "TPC crossed Rows Over Findable Cls; TPC Crossed Rows Over Findable Cls", {HistType::kTH1D, {{100, 0.f, 2.f}}});
      QA.add("QA/hITSchi2", "ITS chi2/Ncls; ITS chi2/Ncls", {HistType::kTH1D, {{100, 0.f, 20.f}}});
      QA.add("QA/hDCAxy", "DCAxy", {HistType::kTH2D, {{200, -0.2f, 0.2f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      QA.add("QA/hDCAz", "DCAz", {HistType::kTH2D, {{200, -0.2f, 0.2f, "DCA z (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      QA.add("QA/TPCChi2VsPZ", "TPCChi2VsPZ", {HistType::kTH2D, {{100, 0.f, 10.f, "p_{TPC}/Z (GeV/c)"}, {120, 0.f, 6.f, "TPC Chi2"}}});
      QA.add("QA/hnSigmaTPCVsPt_El", "n#sigma TPC vs p_{T} for e hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {pTAxis, AxisNSigma}});
      QA.add("QA/hnSigmaTPCVsPt_Pr", "n#sigma TPC vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {pTAxis, AxisNSigma}});
      QA.add("QA/hnSigmaTPCVsPt_De", "n#sigma TPC vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {pTAxis, AxisNSigma}});
      QA.add("QA/hnSigmaTOFVsPt_Pr", "n#sigma TOF vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {pTAxis, AxisNSigma}});
      QA.add("QA/hnSigmaTOFVsPt_De", "n#sigma TOF vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {pTAxis, AxisNSigma}});

      if (!isMC) {
        QA.add("QA/hEtaPr", Form("#eta ditribution for p"), {HistType::kTH1F, {etaAxis}});
        QA.add("QA/hPhiPr", Form("#phi ditribution for p"), {HistType::kTH1F, {phiAxis}});
        QA.add("QA/hEtaAntiPr", Form("#eta ditribution for #bar{p}"), {HistType::kTH1F, {etaAxis}});
        QA.add("QA/hPhiAntiPr", Form("#phi ditribution for #bar{p}"), {HistType::kTH1F, {phiAxis}});
        QA.add("QA/hEtaDe", Form("#eta ditribution for d"), {HistType::kTH1F, {etaAxis}});
        QA.add("QA/hPhiDe", Form("#phi ditribution for d"), {HistType::kTH1F, {phiAxis}});
        QA.add("QA/hEtaAntiDe", Form("#eta ditribution for #bar{d}"), {HistType::kTH1F, {etaAxis}});
        QA.add("QA/hPhiAntiDe", Form("#phi ditribution for #bar{d}"), {HistType::kTH1F, {phiAxis}});

        QA.add("QA/hnSigmaTPCVsPt_Pr_AfterSel", "n#sigma TPC vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsPt_De_AfterSel", "n#sigma TPC vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPt_Pr_AfterSel", "n#sigma TOF vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPt_De_AfterSel", "n#sigma TOF vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {pTAxis, AxisNSigma}});
      }
    }

    if (isMC) {
      registry.add("hReco_EtaPhiPt_Proton", "Gen (anti)protons in reco collisions", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});
      registry.add("hReco_EtaPhiPt_Deuteron", "Gen (anti)deuteron in reco collisions", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});
      registry.add("hReco_PID_EtaPhiPt_Proton", "Gen (anti)protons + PID in reco collisions", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});
      registry.add("hReco_PID_EtaPhiPt_Deuteron", "Gen (anti)deuteron + PID in reco collisions", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});
      registry.add("hReco_EtaPhiPtMC_Proton", "Gen (anti)protons in reco collisions (MC info used)", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});
      registry.add("hReco_EtaPhiPtMC_Deuteron", "Gen (anti)deuteron in reco collisions (MC info used)", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});
      registry.add("hReco_Pt_Proton", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
      registry.add("hReco_Pt_Deuteron", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {pTAxis_small}});

      registry.add("hSec_EtaPhiPt_Proton", "Secondary (anti)protons", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});
      registry.add("hPrimSec_EtaPhiPt_Proton", "Primary + Secondary (anti)protons", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});

      registry.add("hnSigmaTPCVsPt_Pr_MC", "n#sigma TPC vs p_{T} for p hypothesis true MC; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      registry.add("hnSigmaTPCVsPt_De_MC", "n#sigma TPC vs p_{T} for d hypothesis true MC; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      registry.add("hnSigmaTOFVsPt_Pr_MC", "n#sigma TOF vs p_{T} for p hypothesis true MC; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
      registry.add("hnSigmaTOFVsPt_De_MC", "n#sigma TOF vs p_{T} for d hypothesis true MC; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});

      registry.add("hResPt_Proton", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResPt_Deuteron", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResPt_AntiProton", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResPt_AntiDeuteron", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});

      registry.add("hNumeratorPurity_Proton", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
      registry.add("hNumeratorPurity_Deuteron", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
      registry.add("hDenominatorPurity_Proton", " p(#bar{p}); p_{T} (GeV/c);(S + B)", {HistType::kTH1F, {pTAxis_small}});
      registry.add("hDenominatorPurity_Deuteron", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});

      if (doMCQA) {
        registry.add("hResEta_Proton", "; #eta(gen); #eta(reco) - #eta(gen)  ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResEta_Deuteron", "; #eta(gen); #eta(reco) - #eta(gen) ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResPhi_Proton", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi(), "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResPhi_Deuteron", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi(), "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResEta_AntiProton", "; #eta(gen); #eta(reco) - #eta(gen)  ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResEta_AntiDeuteron", "; #eta(gen); #eta(reco) - #eta(gen) ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResPhi_AntiProton", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi(), "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResPhi_AntiDeuteron", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi(), "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});

        registry.add("hNumeratorPurity_Proton_TPC", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Deuteron_TPC", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Proton_TPCTOF", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Deuteron_TPCTOF", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Proton_TPC_or_TOF", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Deuteron_TPC_or_TOF", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Proton_TPCEl_or_TOF", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Proton_TPCEl", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Deuteron_TPCEl", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hNumeratorPurity_Deuteron_TPCEl_or_TOF", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Proton_TPC", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Deuteron_TPC", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Proton_TPCTOF", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Deuteron_TPCTOF", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Proton_TPC_or_TOF", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Deuteron_TPC_or_TOF", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Proton_TPCEl", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Proton_TPCEl_or_TOF", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Deuteron_TPCEl", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hDenominatorPurity_Deuteron_TPCEl_or_TOF", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {pTAxis_small}});

        registry.add("hReco_Pt_Proton_TPC", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Deuteron_TPC", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Proton_TPCTOF", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Deuteron_TPCTOF", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Proton_TPC_or_TOF", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Deuteron_TPC_or_TOF", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Proton_TPCEl", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Proton_TPCEl_or_TOF", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Deuteron_TPCEl", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
        registry.add("hReco_Pt_Deuteron_TPCEl_or_TOF", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {pTAxis_small}});
      }
    }

    if (isMCGen) {
      registry.add("Generated/hNEventsMC", "hNEventsMC", {HistType::kTH1D, {{1, 0.f, 1.f}}});
      registry.get<TH1>(HIST("Generated/hNEventsMC"))->GetXaxis()->SetBinLabel(1, "All");
    }
  }

  // Filters
  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) <= cutzvertex;
  Filter trackFilter = o2::aod::singletrackselector::tpcNClsFound >= min_TPC_nClusters &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedTpcChi2NCl) <= max_chi2_TPC &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::rowsOverFindable>(o2::aod::singletrackselector::storedTpcCrossedRowsOverFindableCls) >= min_TPC_nCrossedRowsOverFindableCls &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedItsChi2NCl) <= max_chi2_ITS &&
                       nabs(o2::aod::singletrackselector::unPack<singletrackselector::binning::dca>(o2::aod::singletrackselector::storedDcaXY)) <= max_dcaxy &&
                       nabs(o2::aod::singletrackselector::unPack<singletrackselector::binning::dca>(o2::aod::singletrackselector::storedDcaXY)) <= max_dcaz &&
                       nabs(o2::aod::singletrackselector::eta) <= etacut;

  template <typename Type>
  bool IsProton(Type const& track, int sign)
  {
    bool isProton = false;

    if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC) {
      if (track.pt() < pTthrpr_TOF) {
        if (sign > 0) {
          if (track.sign() > 0) {
            isProton = true;
          } else if (track.sign() < 0) {
            isProton = false;
          }
        } else if (sign < 0) {
          if (track.sign() > 0) {
            isProton = false;
          } else if (track.sign() < 0) {
            isProton = true;
          }
        }
      } else if (rejectionEl && track.beta() < -100 && track.pt() < pTthrpr_TPCEl && track.tpcNSigmaEl() >= nsigmaElPr) {
        if (sign > 0) {
          if (track.sign() > 0) {
            isProton = true;
          } else if (track.sign() < 0) {
            isProton = false;
          }
        } else if (sign < 0) {
          if (track.sign() > 0) {
            isProton = false;
          } else if (track.sign() < 0) {
            isProton = true;
          }
        }
      } else if (std::abs(track.tofNSigmaPr()) < nsigmaTOF) {
        if (sign > 0) {
          if (track.sign() > 0) {
            isProton = true;
          } else if (track.sign() < 0) {
            isProton = false;
          }
        } else if (sign < 0) {
          if (track.sign() > 0) {
            isProton = false;
          } else if (track.sign() < 0) {
            isProton = true;
          }
        }
      }
    }
    return isProton;
  }

  template <typename Type>
  bool IsDeuteron(Type const& track, int sign)
  {
    bool isDeuteron = false;

    if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC) {
      if (track.pt() < pTthrde_TOF) {
        if (sign > 0) {
          if (track.sign() > 0) {
            isDeuteron = true;
          } else if (track.sign() < 0) {
            isDeuteron = false;
          }
        } else if (sign < 0) {
          if (track.sign() > 0) {
            isDeuteron = false;
          } else if (track.sign() < 0) {
            isDeuteron = true;
          }
        }
      } else if (rejectionEl && track.beta() < -100 && track.pt() < pTthrde_TPCEl && track.tpcNSigmaEl() >= nsigmaElDe) {
        if (sign > 0) {
          if (track.sign() > 0) {
            isDeuteron = true;
          } else if (track.sign() < 0) {
            isDeuteron = false;
          }
        } else if (sign < 0) {
          if (track.sign() > 0) {
            isDeuteron = false;
          } else if (track.sign() < 0) {
            isDeuteron = true;
          }
        }
      } else if (std::abs(track.tofNSigmaDe()) < nsigmaTOF) {
        if (sign > 0) {
          if (track.sign() > 0) {
            isDeuteron = true;
          } else if (track.sign() < 0) {
            isDeuteron = false;
          }
        } else if (sign < 0) {
          if (track.sign() > 0) {
            isDeuteron = false;
          } else if (track.sign() < 0) {
            isDeuteron = true;
          }
        }
      }
    }
    return isDeuteron;
  }

  template <typename T1>
  bool applyDCAcut(const T1& track)
  {
    bool passcut = true;
    // pt-dependent selection
    if (TMath::Abs(track.dcaXY()) > (par0 + par1 / track.pt()))
      passcut = false;
    if (TMath::Abs(track.dcaZ()) > (par0 + par1 / track.pt()))
      passcut = false;

    return passcut;
  }

  template <int ME, bool doMC, typename Type>
  void mixTracks(Type const& tracks1, Type const& tracks2, bool isMCPID)
  { // last value: 0 -- SE; 1 -- ME
    for (auto it1 : tracks1) {
      for (auto it2 : tracks2) {

        // Calculate Delta-eta Delta-phi (reco)
        float deltaEta = it2->eta() - it1->eta();
        float deltaPhi = it2->phi() - it1->phi();
        deltaPhi = getDeltaPhi(deltaPhi);

        // Calculate Delta-eta Delta-phi (gen)
        float deltaEtaGen = -999.;
        float deltaPhiGen = -999.;
        if constexpr (doMC) {
          deltaEtaGen = it2->eta_MC() - it1->eta_MC();
          deltaPhiGen = it2->phi_MC() - it1->phi_MC();
          deltaPhiGen = getDeltaPhi(deltaPhiGen);
        }

        float antipcorr = 1, antidcorr = 1;

        for (int k = 0; k < nBinspT; k++) {

          if (it1->pt() >= pTBins.value.at(k) && it1->pt() < pTBins.value.at(k + 1)) {

            if (docorrection) { // Apply corrections
              antipcorr = hEffpTEta_antiproton->Interpolate(it2->pt(), it2->eta());
              antidcorr = hEffpTEta_antideuteron->Interpolate(it1->pt(), it1->eta());
            }

            if (ME) {
              if constexpr (!doMC) { // Data
                hEtaPhi_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                hCorrEtaPhi_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt(), 1. / (antipcorr * antidcorr));
              } else { // MC
                if (isMCPID) {
                  hPIDEtaPhiRec_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                  hPIDEtaPhiGen_AntiDeAntiPr_ME[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                } else {
                  hEtaPhiRec_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                  hEtaPhiGen_AntiDeAntiPr_ME[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                }
              }
            } else {
              if constexpr (!doMC) { // Data
                hEtaPhi_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                hCorrEtaPhi_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt(), 1. / (antipcorr * antidcorr));
              } else { // MC
                if (isMCPID) {
                  hPIDEtaPhiRec_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                  hPIDEtaPhiGen_AntiDeAntiPr_SE[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                } else {
                  hEtaPhiRec_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                  hEtaPhiGen_AntiDeAntiPr_SE[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                }
              }
            } // SE
          } // pT condition
        } // nBinspT loop
      } // tracks 2
    } // tracks 1
  }

  template <int ME, typename Type>
  void mixMCParticles(Type const& particles1, Type const& particles2)
  {
    for (auto it1 : particles1) {
      for (auto it2 : particles2) {
        // Calculate Delta-eta Delta-phi (gen)
        float deltaEtaGen = it2->eta() - it1->eta();
        float deltaPhiGen = getDeltaPhi(it2->phi() - it1->phi());

        // Loop over pT bins
        for (int k = 0; k < nBinspT; k++) {
          if (it1->pt() >= pTBins.value.at(k) && it1->pt() < pTBins.value.at(k + 1)) {
            // Use correct histogram based on ME flag
            if constexpr (ME) {
              hEtaPhiGen_AntiDeAntiPr_ME[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
            } else {
              hEtaPhiGen_AntiDeAntiPr_SE[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
            }
          }
        }
      }
    }
  }

  float getDeltaPhi(float deltaPhi)
  {
    if (deltaPhi < -TMath::Pi() / 2) {
      return deltaPhi += 2 * TMath::Pi();
    } else if (deltaPhi >= 3 * TMath::Pi() / 2) {
      return deltaPhi -= 2 * TMath::Pi();
    }
    return deltaPhi;
  }

  void GetCorrection(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj, TString filepath, TString histname)
  {
    TList* l = ccdbObj->get<TList>(filepath.Data());
    if (!l) {
      LOGP(error, "Could not open corrections file {}", Form("%s", filepath.Data()));
      return;
    }
    hEffpTEta_proton = static_cast<TH2F*>(l->FindObject(Form("%s_proton", histname.Data())));
    if (!hEffpTEta_proton) {
      LOGP(error, "Could not open histogram {}", Form("%s_proton", histname.Data()));
      return;
    }
    hEffpTEta_antiproton = static_cast<TH2F*>(l->FindObject(Form("%s_antiproton", histname.Data())));
    if (!hEffpTEta_antiproton) {
      LOGP(error, "Could not open histogram {}", Form("%s_antiproton", histname.Data()));
      return;
    }
    hEffpTEta_deuteron = static_cast<TH2F*>(l->FindObject(Form("%s_deuteron", histname.Data())));
    if (!hEffpTEta_deuteron) {
      LOGP(error, "Could not open histogram {}", Form("%s_deuteron", histname.Data()));
      return;
    }
    hEffpTEta_antideuteron = static_cast<TH2F*>(l->FindObject(Form("%s_antideuteron", histname.Data())));
    if (!hEffpTEta_antideuteron) {
      LOGP(error, "Could not open histogram {}", Form("%s_antideuteron", histname.Data()));
      return;
    }
    LOGP(info, "Opened histogram {}", Form("%s_proton", histname.Data()));
    LOGP(info, "Opened histogram {}", Form("%s_antiproton", histname.Data()));
    LOGP(info, "Opened histogram {}", Form("%s_deuteron", histname.Data()));
    LOGP(info, "Opened histogram {}", Form("%s_antideuteron", histname.Data()));
  }

  void processData(FilteredCollisions const& collisions, FilteredTracks const& tracks)
  {
    for (auto track : tracks) {
      if (std::abs(track.template singleCollSel_as<FilteredCollisions>().posZ()) > cutzvertex)
        continue;

      if (track.tpcFractionSharedCls() > max_tpcSharedCls)
        continue;
      if (track.itsNCls() < min_itsNCls)
        continue;
      if (!applyDCAcut(track))
        continue;

      if (doQA) {
        QA.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        QA.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        QA.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        QA.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        QA.fill(HIST("QA/hDCAxy"), track.dcaXY(), track.pt());
        QA.fill(HIST("QA/hDCAz"), track.dcaZ(), track.pt());
        QA.fill(HIST("QA/TPCChi2VsPZ"), track.tpcInnerParam() / track.sign(), track.tpcChi2NCl());
        QA.fill(HIST("QA/hVtxZ_trk"), track.template singleCollSel_as<FilteredCollisions>().posZ());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_El"), track.pt() * track.sign(), track.tpcNSigmaEl());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
      }

      // Discard candidates outside pT of interest
      if (track.pt() > pTBins.value.at(nBinspT) || track.pt() < pTBins.value.at(0))
        continue;

      bool isPr = IsProton(track, +1);
      bool isAntiPr = IsProton(track, -1);
      bool isDe = IsDeuteron(track, +1);
      bool isAntiDe = IsDeuteron(track, -1);

      if (!isPr && !isAntiPr && !isDe && !isAntiDe)
        continue;

      if (isPr && isDe) {
        isDe = 0;
      }
      if (isAntiPr && isAntiDe) {
        isAntiDe = 0;
      }

      // Deuterons Fill & QA
      if (isAntiDe) {
        selectedtracks_antid[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        if (doQA) {
          QA.fill(HIST("QA/hEtaAntiDe"), track.eta());
          QA.fill(HIST("QA/hPhiAntiDe"), track.phi());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
        }
      }
      if (isDe) {
        selectedtracks_d[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        if (doQA) {
          QA.fill(HIST("QA/hEtaDe"), track.eta());
          QA.fill(HIST("QA/hPhiDe"), track.phi());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
        }
      }

      // Protons Fill & QA
      if (isPr) {
        selectedtracks_p[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        if (doQA) {
          QA.fill(HIST("QA/hEtaPr"), track.eta());
          QA.fill(HIST("QA/hPhiPr"), track.phi());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
        }
      } else if (isAntiPr) {
        selectedtracks_antip[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        if (doQA) {
          QA.fill(HIST("QA/hEtaAntiPr"), track.eta());
          QA.fill(HIST("QA/hPhiAntiPr"), track.phi());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
        }
      }
    }

    for (auto collision : collisions) {

      if (std::abs(collision.posZ()) > cutzvertex)
        continue;

      registry.fill(HIST("hNEvents"), 0.5);

      if (selectedtracks_antid.find(collision.globalIndex()) != selectedtracks_antid.end() &&
          selectedtracks_antip.find(collision.globalIndex()) != selectedtracks_antip.end()) {
        registry.fill(HIST("hNEvents"), 1.5);
      }

      if (selectedtracks_d.find(collision.globalIndex()) != selectedtracks_d.end() &&
          selectedtracks_p.find(collision.globalIndex()) != selectedtracks_p.end()) {
        registry.fill(HIST("hNEvents"), 2.5);
      }

      int vertexBinToMix = std::floor((collision.posZ() + cutzvertex) / (2 * cutzvertex / _vertexNbinsToMix));
      int centBinToMix = std::floor(collision.multPerc() / (100.0 / _multNsubBins));

      if (selectedtracks_antid.find(collision.globalIndex()) != selectedtracks_antid.end()) {
        registry.fill(HIST("hNEvents"), 3.5);

        mixbins_antidantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }
    }

    if (!mixbins_antidantip.empty()) {

      for (auto i = mixbins_antidantip.begin(); i != mixbins_antidantip.end(); i++) { // iterating over all vertex&mult bins

        std::vector<colType> value = i->second;
        int EvPerBin = value.size(); // number of collisions in each vertex&mult bin

        for (int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = value[indx1];

          if (selectedtracks_antip.find(col1->index()) != selectedtracks_antip.end()) {
            mixTracks<0, 0>(selectedtracks_antid[col1->index()], selectedtracks_antip[col1->index()], 0); // mixing SE
          }

          for (int indx2 = 0; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = value[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracks_antip.find(col2->index()) != selectedtracks_antip.end()) {
              mixTracks<1, 0>(selectedtracks_antid[col1->index()], selectedtracks_antip[col2->index()], 0); // mixing ME
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

    for (auto i = selectedtracks_d.begin(); i != selectedtracks_d.end(); i++)
      (i->second).clear();
    selectedtracks_d.clear();

    for (auto& pair : mixbins_antidantip) {
      pair.second.clear(); // Clear the vector associated with the key
    }
    mixbins_antidantip.clear(); // Then clear the map itself
  }
  PROCESS_SWITCH(hadronnucleicorrelation, processData, "processData", true);

  void processMC(FilteredCollisions const& collisions, FilteredTracksMC const& tracks)
  {
    for (auto track : tracks) {
      if (std::abs(track.template singleCollSel_as<FilteredCollisions>().posZ()) > cutzvertex)
        continue;

      if (track.tpcFractionSharedCls() > max_tpcSharedCls)
        continue;
      if (track.itsNCls() < min_itsNCls)
        continue;
      if (!applyDCAcut(track))
        continue;

      // Keep only protons and deuterons
      // if (std::abs(track.pdgCode()) != pdgProton && std::abs(track.pdgCode()) != pdgDeuteron)
      // continue;

      if (doQA) {
        QA.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        QA.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        QA.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        QA.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        QA.fill(HIST("QA/hDCAxy"), track.dcaXY(), track.pt());
        QA.fill(HIST("QA/hDCAz"), track.dcaZ(), track.pt());
        QA.fill(HIST("QA/hVtxZ_trk"), track.template singleCollSel_as<FilteredCollisions>().posZ());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_El"), track.pt() * track.sign(), track.tpcNSigmaEl());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
      }

      bool isPr = (IsProton(track, +1) && track.pdgCode() == pdgProton);
      bool isAntiPr = (IsProton(track, -1) && track.pdgCode() == -pdgProton);
      bool isDe = (IsDeuteron(track, +1) && track.pdgCode() == pdgDeuteron);
      bool isAntiDe = (IsDeuteron(track, -1) && track.pdgCode() == -pdgDeuteron);

      if (track.origin() == 1 || track.origin() == 0) { // primaries and secondaries
        if (isPr) {
          registry.fill(HIST("hPrimSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * +1);
          if (track.origin() == 1) { // secondaries
            registry.fill(HIST("hSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * +1);
          }
        }
        if (isAntiPr) {
          registry.fill(HIST("hPrimSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
          if (track.origin() == 1) {
            registry.fill(HIST("hSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
          }
        }
      }

      if (track.origin() != 0)
        continue;

      if (track.pdgCode() == pdgProton) {
        registry.fill(HIST("hReco_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt());
        registry.fill(HIST("hReco_EtaPhiPtMC_Proton"), track.eta_MC(), track.phi_MC(), track.pt_MC());
        registry.fill(HIST("hResPt_Proton"), track.pt_MC(), track.pt() - track.pt_MC());
        if (doMCQA) {
          registry.fill(HIST("hResEta_Proton"), track.eta_MC(), track.eta() - track.eta_MC());
          registry.fill(HIST("hResPhi_Proton"), track.phi_MC(), track.phi() - track.phi_MC());
        }
        if (isPr) {
          registry.fill(HIST("hReco_PID_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt());
        }
        registry.fill(HIST("hnSigmaTPCVsPt_Pr_MC"), track.pt(), track.tpcNSigmaPr());
        registry.fill(HIST("hnSigmaTOFVsPt_Pr_MC"), track.pt(), track.tofNSigmaPr());
      }
      if (track.pdgCode() == -pdgProton) {
        registry.fill(HIST("hReco_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        registry.fill(HIST("hReco_EtaPhiPtMC_Proton"), track.eta_MC(), track.phi_MC(), track.pt_MC() * -1);
        registry.fill(HIST("hResPt_AntiProton"), track.pt_MC(), track.pt() - track.pt_MC());
        if (doMCQA) {
          registry.fill(HIST("hResEta_AntiProton"), track.eta_MC(), track.eta() - track.eta_MC());
          registry.fill(HIST("hResPhi_AntiProton"), track.phi_MC(), track.phi() - track.phi_MC());
        }
        if (isAntiPr) {
          registry.fill(HIST("hReco_PID_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        }
        registry.fill(HIST("hnSigmaTPCVsPt_Pr_MC"), track.pt() * -1, track.tpcNSigmaPr());
        registry.fill(HIST("hnSigmaTOFVsPt_Pr_MC"), track.pt() * -1, track.tofNSigmaPr());
      }
      if (track.pdgCode() == pdgDeuteron) {
        registry.fill(HIST("hReco_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt());
        registry.fill(HIST("hReco_EtaPhiPtMC_Deuteron"), track.eta_MC(), track.phi_MC(), track.pt_MC());
        registry.fill(HIST("hResPt_Deuteron"), track.pt_MC(), track.pt() - track.pt_MC());
        if (doMCQA) {
          registry.fill(HIST("hResEta_Deuteron"), track.eta_MC(), track.eta() - track.eta_MC());
          registry.fill(HIST("hResPhi_Deuteron"), track.phi_MC(), track.phi() - track.phi_MC());
        }
        if (isDe) {
          registry.fill(HIST("hReco_PID_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt());
        }
        registry.fill(HIST("hnSigmaTPCVsPt_De_MC"), track.pt(), track.tpcNSigmaDe());
        registry.fill(HIST("hnSigmaTOFVsPt_De_MC"), track.pt(), track.tofNSigmaDe());
      }
      if (track.pdgCode() == -pdgDeuteron) {
        registry.fill(HIST("hReco_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt() * -1);
        registry.fill(HIST("hReco_EtaPhiPtMC_Deuteron"), track.eta_MC(), track.phi_MC(), track.pt_MC() * -1);
        registry.fill(HIST("hResPt_AntiDeuteron"), track.pt_MC(), track.pt() - track.pt_MC());
        if (doMCQA) {
          registry.fill(HIST("hResEta_AntiDeuteron"), track.eta_MC(), track.eta() - track.eta_MC());
          registry.fill(HIST("hResPhi_AntiDeuteron"), track.phi_MC(), track.phi() - track.phi_MC());
        }
        if (isAntiDe) {
          registry.fill(HIST("hReco_PID_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt() * -1);
        }
        registry.fill(HIST("hnSigmaTPCVsPt_De_MC"), track.pt() * -1, track.tpcNSigmaDe());
        registry.fill(HIST("hnSigmaTOFVsPt_De_MC"), track.pt() * -1, track.tofNSigmaDe());
      }

      // Purity
      // Numerators
      if (isPr) {
        registry.fill(HIST("hNumeratorPurity_Proton"), track.pt());
        registry.fill(HIST("hReco_Pt_Proton"), track.pt());
      }
      if (isAntiPr) {
        registry.fill(HIST("hNumeratorPurity_Proton"), track.pt() * -1);
        registry.fill(HIST("hReco_Pt_Proton"), track.pt() * -1);
      }
      if (isDe) {
        registry.fill(HIST("hNumeratorPurity_Deuteron"), track.pt());
        registry.fill(HIST("hReco_Pt_Deuteron"), track.pt());
      }
      if (isAntiDe) {
        registry.fill(HIST("hNumeratorPurity_Deuteron"), track.pt() * -1);
        registry.fill(HIST("hReco_Pt_Deuteron"), track.pt() * -1);
      }
      if (IsProton(track, +1))
        registry.fill(HIST("hDenominatorPurity_Proton"), track.pt());
      if (IsProton(track, -1))
        registry.fill(HIST("hDenominatorPurity_Proton"), track.pt() * -1);
      if (IsDeuteron(track, +1))
        registry.fill(HIST("hDenominatorPurity_Deuteron"), track.pt());
      if (IsDeuteron(track, -1))
        registry.fill(HIST("hDenominatorPurity_Deuteron"), track.pt() * -1);

      if (doMCQA) {
        // Proton
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.pdgCode() == pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPC"), track.pt());
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF &&
            track.pdgCode() == pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCTOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPCTOF"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.pdgCode() == pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPCEl_or_TOF"), track.pt());
        }

        // AntiProton
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.pdgCode() == -pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPC"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF &&
            track.pdgCode() == -pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCTOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPCTOF"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == -pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPC_or_TOF"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.pdgCode() == -pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == -pdgProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPCEl_or_TOF"), track.pt() * -1);
        }

        // Deuteron
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.pdgCode() == pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPC"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF &&
            track.pdgCode() == pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCTOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPCTOF"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.pdgCode() == pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl_or_TOF"), track.pt());
        }

        // AntiDeuteron
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.pdgCode() == -pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPC"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF &&
            track.pdgCode() == -pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCTOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPCTOF"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == -pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPC_or_TOF"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.pdgCode() == -pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == -pdgDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
        }

        // Denominators
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPC"), track.pt());
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF && track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPCTOF"), track.pt());
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPC_or_TOF"), track.pt());
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl_or_TOF"), track.pt());
        }

        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPC"), track.pt() * -1);
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF && track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPCTOF"), track.pt() * -1);
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPC_or_TOF"), track.pt() * -1);
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl_or_TOF"), track.pt() * -1);
        }

        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC"), track.pt());
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF && track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCTOF"), track.pt());
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl_or_TOF"), track.pt());

        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC"), track.pt() * -1);
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF && track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCTOF"), track.pt() * -1);
        if ((
              (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < -100) ||
              (track.beta() > -100 && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC_or_TOF"), track.pt() * -1);
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < -100) ||
             (track.beta() > -100 && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
      }

      if (!mcCorrelation) {
        continue;
      }

      if (isDe) {
        selectedtracksPIDMC_d[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));
      }
      if (track.pdgCode() == pdgDeuteron) {
        selectedtracksMC_d[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));
      }
      if (isAntiDe) {
        selectedtracksPIDMC_antid[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));
      }
      if (track.pdgCode() == -pdgDeuteron) {
        selectedtracksMC_antid[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));
      }
      if (isPr) {
        selectedtracksPIDMC_p[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));
      }
      if (track.pdgCode() == pdgProton) {
        selectedtracksMC_p[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));
      }
      if (isAntiPr) {
        selectedtracksPIDMC_antip[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));
      }
      if (track.pdgCode() == -pdgProton) {
        selectedtracksMC_antip[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));
      }
    } // track

    if (!mcCorrelation) {
      return;
    }

    for (auto collision : collisions) {
      if (std::abs(collision.posZ()) > cutzvertex)
        continue;
      registry.fill(HIST("hNEvents"), 0.5);

      int vertexBinToMix = std::floor((collision.posZ() + cutzvertex) / (2 * cutzvertex / _vertexNbinsToMix));
      int centBinToMix = std::floor(collision.multPerc() / (100.0 / _multNsubBins));

      if (selectedtracksMC_antid.find(collision.globalIndex()) != selectedtracksMC_antid.end()) {
        mixbins_antidantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }

      if (selectedtracksPIDMC_antid.find(collision.globalIndex()) != selectedtracksPIDMC_antid.end()) {
        mixbinsPID_antidantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }
    } // coll

    if (!mixbins_antidantip.empty()) {

      for (auto i = mixbins_antidantip.begin(); i != mixbins_antidantip.end(); i++) { // iterating over all vertex&mult bins

        std::vector<colType> value = i->second;
        int EvPerBin = value.size(); // number of collisions in each vertex&mult bin

        for (int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = value[indx1];

          if (selectedtracksMC_antip.find(col1->index()) != selectedtracksMC_antip.end()) {
            mixTracks<0, 1>(selectedtracksMC_antid[col1->index()], selectedtracksMC_antip[col1->index()], 0); // mixing SE
          }

          for (int indx2 = indx1 + 1; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracksMC_antip.find(col2->index()) != selectedtracksMC_antip.end()) {
              mixTracks<1, 1>(selectedtracksMC_antid[col1->index()], selectedtracksMC_antip[col2->index()], 0); // mixing ME
            }
          }
        }
      }
    }

    if (!mixbinsPID_antidantip.empty()) {

      for (auto i = mixbinsPID_antidantip.begin(); i != mixbinsPID_antidantip.end(); i++) { // iterating over all vertex&mult bins

        std::vector<colType> value = i->second;
        int EvPerBin = value.size(); // number of collisions in each vertex&mult bin

        for (int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = value[indx1];

          if (selectedtracksPIDMC_antip.find(col1->index()) != selectedtracksPIDMC_antip.end()) {
            mixTracks<0, 1>(selectedtracksPIDMC_antid[col1->index()], selectedtracksPIDMC_antip[col1->index()], 1); // mixing SE
          }

          for (int indx2 = indx1 + 1; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracksPIDMC_antip.find(col2->index()) != selectedtracksPIDMC_antip.end()) {
              mixTracks<1, 1>(selectedtracksPIDMC_antid[col1->index()], selectedtracksPIDMC_antip[col2->index()], 1); // mixing ME
            }
          }
        }
      }
    }

    // clearing up
    for (auto i = selectedtracksMC_antid.begin(); i != selectedtracksMC_antid.end(); i++)
      (i->second).clear();
    selectedtracksMC_antid.clear();

    for (auto i = selectedtracksMC_d.begin(); i != selectedtracksMC_d.end(); i++)
      (i->second).clear();
    selectedtracksMC_d.clear();

    for (auto i = selectedtracksMC_antip.begin(); i != selectedtracksMC_antip.end(); i++)
      (i->second).clear();
    selectedtracksMC_antip.clear();

    for (auto i = selectedtracksMC_p.begin(); i != selectedtracksMC_p.end(); i++)
      (i->second).clear();
    selectedtracksMC_p.clear();

    for (auto i = selectedtracksPIDMC_antid.begin(); i != selectedtracksPIDMC_antid.end(); i++)
      (i->second).clear();
    selectedtracksPIDMC_antid.clear();

    for (auto i = selectedtracksPIDMC_d.begin(); i != selectedtracksPIDMC_d.end(); i++)
      (i->second).clear();
    selectedtracksPIDMC_d.clear();

    for (auto i = selectedtracksPIDMC_antip.begin(); i != selectedtracksPIDMC_antip.end(); i++)
      (i->second).clear();
    selectedtracksPIDMC_antip.clear();

    for (auto i = selectedtracksPIDMC_p.begin(); i != selectedtracksPIDMC_p.end(); i++)
      (i->second).clear();
    selectedtracksPIDMC_p.clear();

    for (auto& pair : mixbinsPID_antidantip) {
      pair.second.clear(); // clear the vector associated with the key
    }
    mixbinsPID_antidantip.clear(); // clear the map

    for (auto& pair : mixbins_antidantip) {
      pair.second.clear(); // clear the vector associated with the key
    }
    mixbins_antidantip.clear(); // clear the map
  }
  PROCESS_SWITCH(hadronnucleicorrelation, processMC, "processMC", false);

  void processGen(SimCollisions const& mcCollisions,
                  SimParticles const& mcParticles)
  {
    for (auto particle : mcParticles) {

      if (isPrim && !particle.isPhysicalPrimary()) {
        continue;
      }

      if (std::abs(particle.eta()) > etacut) {
        continue;
      }

      if (particle.pdgCode() == pdgDeuteron) {
        selectedparticlesMC_d[particle.mcCollisionId()].push_back(std::make_shared<decltype(particle)>(particle));
      }
      if (particle.pdgCode() == -pdgDeuteron) {
        selectedparticlesMC_antid[particle.mcCollisionId()].push_back(std::make_shared<decltype(particle)>(particle));
      }
      if (particle.pdgCode() == pdgProton) {
        if (!particle.has_daughters()) {
          selectedparticlesMC_p[particle.mcCollisionId()].push_back(std::make_shared<decltype(particle)>(particle));
        } else {
          bool isp = false;

          for (auto& dau : particle.daughters_as<aod::McParticles>()) {
            if (dau.pdgCode() == pdgProton) {
              isp = true;
            }
          }

          if (isp) {
            selectedparticlesMC_p[particle.mcCollisionId()].push_back(std::make_shared<decltype(particle)>(particle));
          }
        }
      }
      if (particle.pdgCode() == -pdgProton) {
        if (!particle.has_daughters()) {
          selectedparticlesMC_antip[particle.mcCollisionId()].push_back(std::make_shared<decltype(particle)>(particle));
        } else {
          bool isantip = false;

          for (auto& dau : particle.daughters_as<aod::McParticles>()) {
            if (dau.pdgCode() == -pdgProton) {
              isantip = true;
            }
          }

          if (isantip) {
            selectedparticlesMC_antip[particle.mcCollisionId()].push_back(std::make_shared<decltype(particle)>(particle));
          }
        }
      }
    }

    for (auto collision1 : mcCollisions) { // loop on collisions

      registry.fill(HIST("Generated/hNEventsMC"), 0.5);

      // anti-d - anti-p correlation
      if (selectedparticlesMC_antid.find(collision1.globalIndex()) != selectedparticlesMC_antid.end()) {
        if (selectedparticlesMC_antip.find(collision1.globalIndex()) != selectedparticlesMC_antip.end()) {
          mixMCParticles<0>(selectedparticlesMC_antid[collision1.globalIndex()], selectedparticlesMC_antip[collision1.globalIndex()]); // mixing SE
        }

        int stop1 = 0;

        for (auto collision2 : mcCollisions) { // nested loop on collisions

          if (collision1.globalIndex() == collision2.globalIndex()) {
            continue;
          }

          if (stop1 > maxmixcollsGen) {
            break;
          }

          if (selectedparticlesMC_antip.find(collision2.globalIndex()) != selectedparticlesMC_antip.end()) {
            mixMCParticles<1>(selectedparticlesMC_antid[collision1.globalIndex()], selectedparticlesMC_antip[collision2.globalIndex()]); // mixing ME
          }

          stop1++;
        }
      }

      // d - p correlation
      if (domatterGen) {
        if (selectedparticlesMC_d.find(collision1.globalIndex()) != selectedparticlesMC_d.end()) {
          if (selectedparticlesMC_p.find(collision1.globalIndex()) != selectedparticlesMC_p.end()) {
            mixMCParticles<0>(selectedparticlesMC_d[collision1.globalIndex()], selectedparticlesMC_p[collision1.globalIndex()]); // mixing SE
          }

          int stop2 = 0;

          for (auto collision2 : mcCollisions) { // nested loop on collisions

            if (collision1.globalIndex() == collision2.globalIndex()) {
              continue;
            }

            if (stop2 > maxmixcollsGen) {
              break;
            }

            if (selectedparticlesMC_p.find(collision2.globalIndex()) != selectedparticlesMC_p.end()) {
              mixMCParticles<1>(selectedparticlesMC_d[collision1.globalIndex()], selectedparticlesMC_p[collision2.globalIndex()]); // mixing ME
            }

            stop2++;
          }
        }
      }
    }

    // clearing up
    for (auto i = selectedparticlesMC_antid.begin(); i != selectedparticlesMC_antid.end(); i++)
      (i->second).clear();
    selectedparticlesMC_antid.clear();

    for (auto i = selectedparticlesMC_d.begin(); i != selectedparticlesMC_d.end(); i++)
      (i->second).clear();
    selectedparticlesMC_d.clear();

    for (auto i = selectedparticlesMC_antip.begin(); i != selectedparticlesMC_antip.end(); i++)
      (i->second).clear();
    selectedparticlesMC_antip.clear();

    for (auto i = selectedparticlesMC_p.begin(); i != selectedparticlesMC_p.end(); i++)
      (i->second).clear();
    selectedparticlesMC_p.clear();
  }
  PROCESS_SWITCH(hadronnucleicorrelation, processGen, "processGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<hadronnucleicorrelation>(cfgc)};
}
