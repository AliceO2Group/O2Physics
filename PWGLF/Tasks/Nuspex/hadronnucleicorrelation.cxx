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
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <vector>
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
  Configurable<bool> isMC{"isMC", false, "is MC"};
  Configurable<bool> mcCorrelation{"mcCorrelation", false, "true: build the correlation function only for SE"};
  Configurable<bool> disable_pantip{"disable_pantip", false, "disable_pantip"};
  Configurable<bool> docorrection{"docorrection", false, "do efficiency correction"};

  Configurable<std::string> fCorrectionPath{"fCorrectionPath", "", "Correction path to file"};
  Configurable<std::string> fCorrectionHisto{"fCorrectionHisto", "", "Correction histogram"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

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
  Configurable<float> pTthrpr_TOF{"pTthrpr_TOF", 0.8f, "threshold pT proton to use TOF"};
  Configurable<float> debug_ptthrp{"debug_ptthrp", 2.0f, "threshold pT proton for phi/eta debug"};
  Configurable<float> debug_ptthrd{"debug_ptthrd", 2.0f, "threshold pT deuteron for phi/eta debug"};
  Configurable<float> max_tpcSharedCls{"max_tpcSharedCls", 0.4, "maximum fraction of TPC shared clasters"};
  Configurable<int> min_itsNCls{"min_itsNCls", 0, "minimum allowed number of ITS clasters"};
  Configurable<float> threta{"threta", 0.1, "threshold for debug DeltaEta"};
  Configurable<float> thrphi{"thrphi", 0.5, "threshold for debug DeltaPhi"};

  // Mixing parameters
  Configurable<int> _vertexNbinsToMix{"vertexNbinsToMix", 10, "Number of vertexZ bins for the mixing"};
  Configurable<int> _multNsubBins{"multSubBins", 10, "number of sub-bins to perform the mixing within"};

  // pT/A bins
  Configurable<std::vector<double>> pTBins{"pTBins", {0.4f, 0.6f, 0.8f}, "p_{T} bins"};
  Configurable<std::vector<double>> etaBins{"etaBins", {-0.8f, 0.f, 0.8f}, "#eta bins"};
  Configurable<std::vector<double>> phiBins{"phiBins", {0.f, TMath::Pi(), 2 * TMath::Pi()}, "#phi bins"};

  ConfigurableAxis AxisNSigma{"AxisNSigma", {50, -10.f, 10.f}, "n#sigma"};

  using FilteredCollisions = soa::Filtered<aod::SingleCollSels>;
  using FilteredTracks = soa::Filtered<aod::SingleTrackSels>;
  using FilteredTracksMC = soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkMCs>>;

  HistogramRegistry registry{"registry"};
  HistogramRegistry QA{"QA"};

  typedef std::shared_ptr<FilteredTracks::iterator> trkType;
  typedef std::shared_ptr<FilteredTracksMC::iterator> trkTypeMC;
  typedef std::shared_ptr<FilteredCollisions::iterator> colType;

  // key: int64_t - value: vector of trkType objects
  std::map<int64_t, std::vector<trkType>> selectedtracks_p;
  std::map<int64_t, std::vector<trkType>> selectedtracks_antid;
  std::map<int64_t, std::vector<trkType>> selectedtracks_antip;

  // key: int64_t - value: vector of trkType objects
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksMC_p;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksMC_antid;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksMC_antip;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksPIDMC_p;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksPIDMC_antid;
  std::map<int64_t, std::vector<trkTypeMC>> selectedtracksPIDMC_antip;

  // key: pair of an integer and a float - value: vector of colType objects
  // for each key I have a vector of collisions
  std::map<std::pair<int, float>, std::vector<colType>> mixbins_antidantip;
  std::map<std::pair<int, float>, std::vector<colType>> mixbins_pantip;
  std::map<std::pair<int, float>, std::vector<colType>> mixbinsPID_antidantip;
  std::map<std::pair<int, float>, std::vector<colType>> mixbinsPID_pantip;

  std::vector<std::shared_ptr<TH3>> hEtaPhi_PrAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhi_PrAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hEtaPhi_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhi_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhi_PrAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhi_PrAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhi_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhi_AntiDeAntiPr_ME;

  std::vector<std::shared_ptr<TH3>> hEtaPhiRec_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhiGen_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhiRec_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hEtaPhiGen_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hEtaPhiRec_PrAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhiGen_PrAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhiRec_PrAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hEtaPhiGen_PrAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiRec_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiGen_AntiDeAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiRec_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiGen_AntiDeAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiRec_PrAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiGen_PrAntiPr_SE;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiRec_PrAntiPr_ME;
  std::vector<std::shared_ptr<TH3>> hPIDEtaPhiGen_PrAntiPr_ME;

  int nBinspT, nBinseta, nBinsphi;
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
    AxisSpec etaBinnedAxis = {etaBins, "#eta"};
    AxisSpec phiBinnedAxis = {phiBins, "#phi"};
    AxisSpec etaAxis = {100, -1.5, 1.5, "#Delta#eta"};
    AxisSpec phiAxis = {60, -TMath::Pi() / 2, 1.5 * TMath::Pi(), "#Delta#phi"};
    AxisSpec pTAxis = {200, -10.f, 10.f, "p_{T} GeV/c"};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1I, {{3, 0.f, 3.f}}});
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "Selected");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "#bar{d}-#bar{p}");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "p-#bar{p}");

    nBinspT = pTBins.value.size() - 1;
    nBinseta = etaBins.value.size() - 1;
    nBinsphi = phiBins.value.size() - 1;

    if (mcCorrelation) {
      for (int i = 0; i < nBinspT; i++) {
        auto htempSERec_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhiRec_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                         Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempSEGen_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhiGen_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                         Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hEtaPhiRec_AntiDeAntiPr_SE.push_back(std::move(htempSERec_AntiDeAntiPr));
        hEtaPhiGen_AntiDeAntiPr_SE.push_back(std::move(htempSEGen_AntiDeAntiPr));
        auto htempMERec_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhiRec_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                         Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempMEGen_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhiGen_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                         Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hEtaPhiRec_AntiDeAntiPr_ME.push_back(std::move(htempMERec_AntiDeAntiPr));
        hEtaPhiGen_AntiDeAntiPr_ME.push_back(std::move(htempMEGen_AntiDeAntiPr));

        auto htempSERec_PrAntiPr = registry.add<TH3>(Form("hEtaPhiRec_PrAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                     Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} p <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempSEGen_PrAntiPr = registry.add<TH3>(Form("hEtaPhiGen_PrAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                     Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} p <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hEtaPhiRec_PrAntiPr_SE.push_back(std::move(htempSERec_PrAntiPr));
        hEtaPhiGen_PrAntiPr_SE.push_back(std::move(htempSEGen_PrAntiPr));
        auto htempMERec_PrAntiPr = registry.add<TH3>(Form("hEtaPhiRec_PrAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                     Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} p <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempMEGen_PrAntiPr = registry.add<TH3>(Form("hEtaPhiGen_PrAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                     Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} p <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hEtaPhiRec_PrAntiPr_ME.push_back(std::move(htempMERec_PrAntiPr));
        hEtaPhiGen_PrAntiPr_ME.push_back(std::move(htempMEGen_PrAntiPr));

        auto htempPIDSERec_AntiDeAntiPr = registry.add<TH3>(Form("hPIDEtaPhiRec_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                            Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempPIDSEGen_AntiDeAntiPr = registry.add<TH3>(Form("hPIDEtaPhiGen_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                            Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hPIDEtaPhiRec_AntiDeAntiPr_SE.push_back(std::move(htempPIDSERec_AntiDeAntiPr));
        hPIDEtaPhiGen_AntiDeAntiPr_SE.push_back(std::move(htempPIDSEGen_AntiDeAntiPr));
        auto htempPIDMERec_AntiDeAntiPr = registry.add<TH3>(Form("hPIDEtaPhiRec_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                            Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempPIDMEGen_AntiDeAntiPr = registry.add<TH3>(Form("hPIDEtaPhiGen_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                            Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hPIDEtaPhiRec_AntiDeAntiPr_ME.push_back(std::move(htempPIDMERec_AntiDeAntiPr));
        hPIDEtaPhiGen_AntiDeAntiPr_ME.push_back(std::move(htempPIDMEGen_AntiDeAntiPr));

        auto htempPIDSERec_PrAntiPr = registry.add<TH3>(Form("hPIDEtaPhiRec_PrAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                        Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} p <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempPIDSEGen_PrAntiPr = registry.add<TH3>(Form("hPIDEtaPhiGen_PrAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                        Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} p <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hPIDEtaPhiRec_PrAntiPr_SE.push_back(std::move(htempPIDSERec_PrAntiPr));
        hPIDEtaPhiGen_PrAntiPr_SE.push_back(std::move(htempPIDSEGen_PrAntiPr));
        auto htempPIDMERec_PrAntiPr = registry.add<TH3>(Form("hPIDEtaPhiRec_PrAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                        Form("Rec #Delta#eta#Delta#phi (%.1f<p_{T} p <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempPIDMEGen_PrAntiPr = registry.add<TH3>(Form("hPIDEtaPhiGen_PrAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10),
                                                        Form("Gen #Delta#eta#Delta#phi (%.1f<p_{T} p <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hPIDEtaPhiRec_PrAntiPr_ME.push_back(std::move(htempPIDMERec_PrAntiPr));
        hPIDEtaPhiGen_PrAntiPr_ME.push_back(std::move(htempPIDMEGen_PrAntiPr));
      }
    }

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

      for (int i = 0; i < nBinspT; i++) {
        if (!disable_pantip) {
          auto htempSE_PrAntiPr = registry.add<TH3>(Form("hEtaPhi_PrAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta#eta#Delta#phi (%.1f<p_{T} pr <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
          auto htempME_PrAntiPr = registry.add<TH3>(Form("hEtaPhi_PrAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta#eta#Delta#phi (%.1f<p_{T} pr <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});

          hEtaPhi_PrAntiPr_SE.push_back(std::move(htempSE_PrAntiPr));
          hEtaPhi_PrAntiPr_ME.push_back(std::move(htempME_PrAntiPr));
        }

        auto htempSE_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto htempME_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hEtaPhi_AntiDeAntiPr_SE.push_back(std::move(htempSE_AntiDeAntiPr));
        hEtaPhi_AntiDeAntiPr_ME.push_back(std::move(htempME_AntiDeAntiPr));

        if (!disable_pantip) {
          auto hCorrtempSE_PrAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_PrAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T} pr <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
          auto hCorrtempME_PrAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_PrAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T} pr <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});

          hCorrEtaPhi_PrAntiPr_SE.push_back(std::move(hCorrtempSE_PrAntiPr));
          hCorrEtaPhi_PrAntiPr_ME.push_back(std::move(hCorrtempME_PrAntiPr));
        }

        auto hCorrtempSE_AntiDeAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_AntiDeAntiPr_SE_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        auto hCorrtempME_AntiDeAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_AntiDeAntiPr_ME_pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T} #bar{d} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {etaAxis, phiAxis, ptBinnedAxis}});
        hCorrEtaPhi_AntiDeAntiPr_SE.push_back(std::move(hCorrtempSE_AntiDeAntiPr));
        hCorrEtaPhi_AntiDeAntiPr_ME.push_back(std::move(hCorrtempME_AntiDeAntiPr));
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
      QA.add("QA/Pt_Pr", "Selected protons; p_{T} (GeV/c)", {HistType::kTH1F, {{100, 0.f, 10.f, "p_{T} (GeV/c)"}}});
      QA.add("QA/Pt_AntiDe", "Selected antideuterons; p_{T} (GeV/c)", {HistType::kTH1F, {{100, 0.f, 10.f, "p_{T} (GeV/c)"}}});
      QA.add("QA/Pt_AntiPr", "Selected antiprotons; p_{T} (GeV/c)", {HistType::kTH1F, {{100, 0.f, 10.f, "p_{T} (GeV/c)"}}});

      if (!isMC) {
        QA.add("QA/hnSigmaTPCVsPt_Pr_Debug", "n#sigma TPC vs p_{T} for p hypothesis; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsPt_De_Debug", "n#sigma TPC vs p_{T} for d hypothesis; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTPCVsPt_APrDe_Debug", "n#sigma TPC vs p_{T} for d hypothesis; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPt_Pr_Debug", "n#sigma TOF vs p_{T} for p hypothesis; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPt_De_Debug", "n#sigma TOF vs p_{T} for d hypothesis; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPt_APrDe_Debug", "n#sigma TOF vs p_{T} for d hypothesis; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {pTAxis, AxisNSigma}});
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
        QA.add("QA/hnSigmaTOFVsPhi_Pr", Form("n#sigma TOF vs #phi p; #phi; n#sigma TOF"), {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi()}, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsPhi_AntiPr", Form("n#sigma TOF vs #phi #bar{p}; #phi; n#sigma TOF"), {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi()}, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsEta_Pr", Form("n#sigma TOF vs #eta p; #eta; n#sigma TOF"), {HistType::kTH2F, {{100, -1.f, +1.f}, AxisNSigma}});
        QA.add("QA/hnSigmaTOFVsEta_AntiPr", Form("n#sigma TOF vs #eta #bar{p}; #eta; n#sigma TOF"), {HistType::kTH2F, {{100, -1.f, +1.f}, AxisNSigma}});
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

      registry.add("hResPt_Proton", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResPt_Deuteron", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResEta_Proton", "; #eta(gen); #eta(reco) - #eta(gen)  ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
      registry.add("hResEta_Deuteron", "; #eta(gen); #eta(reco) - #eta(gen) ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
      registry.add("hResPhi_Proton", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi(), "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
      registry.add("hResPhi_Deuteron", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi(), "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
      registry.add("hResPt_AntiProton", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResPt_AntiDeuteron", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResEta_AntiProton", "; #eta(gen); #eta(reco) - #eta(gen)  ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
      registry.add("hResEta_AntiDeuteron", "; #eta(gen); #eta(reco) - #eta(gen) ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
      registry.add("hResPhi_AntiProton", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi(), "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
      registry.add("hResPhi_AntiDeuteron", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, 2 * TMath::Pi(), "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
    }
  }

  // Filters
  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) <= cutzvertex;
  Filter trackFilter = o2::aod::singletrackselector::tpcNClsFound >= min_TPC_nClusters &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedTpcChi2NCl) <= max_chi2_TPC &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::rowsOverFindable>(o2::aod::singletrackselector::storedTpcCrossedRowsOverFindableCls) >= min_TPC_nCrossedRowsOverFindableCls &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedItsChi2NCl) <= max_chi2_ITS &&
                       //  nabs(o2::aod::singletrackselector::unPack<singletrackselector::binning::dca>(o2::aod::singletrackselector::storedDcaXY_v1)) <= max_dcaxy && // For now no filtering on the DCAxy or DCAz (casting not supported)
                       //  nabs(o2::aod::singletrackselector::unPack<singletrackselector::binning::dca>(o2::aod::singletrackselector::storedDcaXY_v1)) <= max_dcaz && // For now no filtering on the DCAxy or DCAz (casting not supported)
                       nabs(o2::aod::singletrackselector::eta) <= etacut;

  template <int ME, bool doMC, typename Type>
  void mixTracks(Type const& tracks1, Type const& tracks2, bool isDe, bool isMCPID)
  { // last value: 0 -- SE; 1 -- ME
    for (auto it1 : tracks1) {
      for (auto it2 : tracks2) {

        // Variables
        float deltaEta = it2->eta() - it1->eta();
        float deltaPhi = it2->phi() - it1->phi();
        deltaPhi = getDeltaPhi(deltaPhi);

        float deltaEtaGen = -999.;
        float deltaPhiGen = -999.;
        if constexpr (doMC) {
          deltaEtaGen = it2->eta_MC() - it1->eta_MC();
          deltaPhiGen = it2->phi_MC() - it1->phi_MC();
          deltaPhiGen = getDeltaPhi(deltaPhiGen);
        }

        float pcorr = 1, antipcorr = 1, antidcorr = 1;

        for (int k = 0; k < nBinspT; k++) {

          if (!isDe && !disable_pantip) {
            if (it1->pt() > pTBins.value.at(k) && it1->pt() <= pTBins.value.at(k + 1)) {

              if (doQA) {
                QA.fill(HIST("QA/Pt_Pr"), it1->pt());
                QA.fill(HIST("QA/Pt_AntiPr"), it2->pt());
              }

              if (docorrection) {
                pcorr = hEffpTEta_proton->Interpolate(it1->pt(), it1->eta());
                antipcorr = hEffpTEta_antiproton->Interpolate(it2->pt(), it2->eta());
              }

              if (ME) {
                if constexpr (!doMC) {
                  hEtaPhi_PrAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                  hCorrEtaPhi_PrAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt(), 1. / (pcorr * antipcorr));
                }
                if constexpr (doMC) {
                  if (isMCPID) {
                    hPIDEtaPhiRec_PrAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                    hPIDEtaPhiGen_PrAntiPr_ME[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                  } else {
                    hEtaPhiRec_PrAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                    hEtaPhiGen_PrAntiPr_ME[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                  }
                }
              } else {
                if constexpr (!doMC) {
                  hEtaPhi_PrAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                  hCorrEtaPhi_PrAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt(), 1. / (pcorr * antipcorr));

                  if (doQA && abs(deltaEta) < threta && abs(deltaPhi) < thrphi) {
                    QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_Debug"), it1->pt(), it1->tpcNSigmaPr());
                    QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_Debug"), -1.f * it2->pt(), it2->tpcNSigmaPr());
                    QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr_Debug"), it1->pt(), it1->tofNSigmaPr());
                    QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr_Debug"), -1.f * it2->pt(), it2->tofNSigmaPr());
                  }
                }
                if constexpr (doMC) {
                  if (isMCPID) {
                    hPIDEtaPhiRec_PrAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                    hPIDEtaPhiGen_PrAntiPr_SE[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                  } else {
                    hEtaPhiRec_PrAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                    hEtaPhiGen_PrAntiPr_SE[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                  }
                }
              }
            }
          } else {
            if (it1->pt() > pTBins.value.at(k) && it1->pt() <= pTBins.value.at(k + 1)) {

              if (doQA) {
                QA.fill(HIST("QA/Pt_AntiDe"), it1->pt());
                QA.fill(HIST("QA/Pt_AntiPr"), it2->pt());
              }

              if (docorrection) {
                antipcorr = hEffpTEta_antiproton->Interpolate(it2->pt(), it2->eta());
                antidcorr = hEffpTEta_antideuteron->Interpolate(it1->pt(), it1->eta());
              }

              if (ME) {
                if constexpr (!doMC) {
                  hEtaPhi_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                  hCorrEtaPhi_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt(), 1. / (antipcorr * antidcorr));
                }
                if constexpr (doMC) {
                  if (isMCPID) {
                    hPIDEtaPhiRec_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                    hPIDEtaPhiGen_AntiDeAntiPr_ME[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                  } else {
                    hEtaPhiRec_AntiDeAntiPr_ME[k]->Fill(deltaEta, deltaPhi, it2->pt());
                    hEtaPhiGen_AntiDeAntiPr_ME[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                  }
                }
              } else {
                if constexpr (!doMC) {
                  hEtaPhi_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                  hCorrEtaPhi_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt(), 1. / (antipcorr * antidcorr));

                  if (doQA && abs(deltaEta) < threta && abs(deltaPhi) < thrphi) {
                    QA.fill(HIST("QA/hnSigmaTPCVsPt_De_Debug"), -1.f * it1->pt(), it1->tpcNSigmaDe());
                    QA.fill(HIST("QA/hnSigmaTPCVsPt_APrDe_Debug"), -1.f * it2->pt(), it2->tpcNSigmaPr());
                    QA.fill(HIST("QA/hnSigmaTOFVsPt_De_Debug"), -1.f * it1->pt(), it1->tofNSigmaDe());
                    QA.fill(HIST("QA/hnSigmaTOFVsPt_APrDe_Debug"), -1.f * it2->pt(), it2->tofNSigmaPr());
                  }
                }
                if constexpr (doMC) {
                  if (isMCPID) {
                    hPIDEtaPhiRec_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                    hPIDEtaPhiGen_AntiDeAntiPr_SE[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                  } else {
                    hEtaPhiRec_AntiDeAntiPr_SE[k]->Fill(deltaEta, deltaPhi, it2->pt());
                    hEtaPhiGen_AntiDeAntiPr_SE[k]->Fill(deltaEtaGen, deltaPhiGen, it2->pt());
                  }
                }
              } // SE
            }
          }
        } // nBinspT loop
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
      if (abs(track.template singleCollSel_as<FilteredCollisions>().posZ()) > cutzvertex)
        continue;
      if (abs(track.dcaXY()) > max_dcaxy || abs(track.dcaZ()) > max_dcaz) { // For now no filtering on the DCAxy or DCAz (casting not supported)
        continue;
      }
      if (doQA) {
        QA.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        QA.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        QA.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        QA.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        QA.fill(HIST("QA/hDCAxy"), track.dcaXY());
        QA.fill(HIST("QA/hDCAz"), track.dcaZ());
        QA.fill(HIST("QA/hVtxZ_trk"), track.template singleCollSel_as<FilteredCollisions>().posZ());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
      }

      if (track.pt() > pTBins.value.at(nBinspT) || track.pt() < pTBins.value.at(0))
        continue;

      // Additional track cuts
      if (track.tpcFractionSharedCls() > max_tpcSharedCls || track.itsNCls() < min_itsNCls)
        continue;

      bool isPr = false;
      bool isAntiPr = false;
      bool isDeTPCTOF = false;
      bool isAntiDeTPCTOF = false;

      if (TMath::Abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() > 0) {
        if (track.pt() < pTthrpr_TOF) {
          isPr = true;
        } else if (TMath::Abs(track.tofNSigmaPr()) < nsigmaTOF) {
          isPr = true;
        }
      }
      if (TMath::Abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() < 0) {
        if (track.pt() < pTthrpr_TOF) {
          isAntiPr = true;
        } else if (TMath::Abs(track.tofNSigmaPr()) < nsigmaTOF) {
          isAntiPr = true;
        }
      }
      if (TMath::Abs(track.tpcNSigmaDe()) < nsigmaTPC && TMath::Abs(track.tpcNSigmaPr()) >= nsigmaTPC) {
        if (TMath::Abs(track.tofNSigmaDe()) < nsigmaTOF && TMath::Abs(track.tofNSigmaPr()) >= nsigmaTOF) {
          if (track.sign() > 0) {
            isDeTPCTOF = true;
          } else if (track.sign() < 0) {
            isAntiDeTPCTOF = true;
          }
        }
      }

      if (!isPr && !isAntiPr && !isDeTPCTOF && !isAntiDeTPCTOF)
        continue;

      if (isPr && isDeTPCTOF) {
        isDeTPCTOF = 0;
      }
      if (isAntiPr && isAntiDeTPCTOF) {
        isAntiDeTPCTOF = 0;
      }

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
          QA.fill(HIST("QA/hnSigmaTOFVsPhi_Pr"), track.phi(), track.tofNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTOFVsEta_Pr"), track.eta(), track.tofNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
        }
      } else if (isAntiPr) {
        selectedtracks_antip[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        if (doQA) {
          QA.fill(HIST("QA/hEtaAntiPr"), track.eta());
          QA.fill(HIST("QA/hPhiAntiPr"), track.phi());
          QA.fill(HIST("QA/hnSigmaTPCVsPhi_AntiPr"), track.phi(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsEta_AntiPr"), track.eta(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTOFVsPhi_AntiPr"), track.phi(), track.tofNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTOFVsEta_AntiPr"), track.eta(), track.tofNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
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
              mixTracks<0, 0>(selectedtracks_p[col1->index()], selectedtracks_antip[col1->index()], 0, 0); // mixing SE
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
                mixTracks<1, 0>(selectedtracks_p[col1->index()], selectedtracks_antip[col2->index()], 0, 0); // mixing ME
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
            mixTracks<0, 0>(selectedtracks_antid[col1->index()], selectedtracks_antip[col1->index()], 1, 0); // mixing SE
          }

          for (int indx2 = 0; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = value[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracks_antip.find(col2->index()) != selectedtracks_antip.end()) {
              registry.fill(HIST("hDebugdp"), 4.5);
              mixTracks<1, 0>(selectedtracks_antid[col1->index()], selectedtracks_antip[col2->index()], 1, 0); // mixing ME
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

    for (auto& pair : mixbins_antidantip) {
      pair.second.clear(); // Clear the vector associated with the key
    }
    mixbins_antidantip.clear(); // Then clear the map itself

    for (auto& pair : mixbins_antidantip) {
      pair.second.clear(); // Clear the vector associated with the key
    }
    mixbins_antidantip.clear(); // Then clear the map itself

    for (auto& pair : mixbins_pantip) {
      pair.second.clear(); // Clear the vector associated with the key
    }
    mixbins_pantip.clear(); // Then clear the map itself
  }
  PROCESS_SWITCH(hadronnucleicorrelation, processData, "processData", true);

  void processMC(FilteredCollisions const& collisions, FilteredTracksMC const& tracks)
  {
    for (auto track : tracks) {
      if (abs(track.template singleCollSel_as<FilteredCollisions>().posZ()) > cutzvertex)
        continue;

      if (abs(track.dcaXY()) > max_dcaxy || abs(track.dcaZ()) > max_dcaz) { // For now no filtering on the DCAxy or DCAz (casting not supported)
        continue;
      }

      if (doQA) {
        QA.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        QA.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        QA.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        QA.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        QA.fill(HIST("QA/hDCAxy"), track.dcaXY());
        QA.fill(HIST("QA/hDCAz"), track.dcaZ());
        QA.fill(HIST("QA/hVtxZ_trk"), track.template singleCollSel_as<FilteredCollisions>().posZ());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
      }

      if (track.origin() != 0)
        continue;

      if (abs(track.pdgCode()) != pdgProton && abs(track.pdgCode()) != pdgDeuteron)
        continue;

      bool isPr = false;
      bool isAntiPr = false;
      bool isAntiDeTPCTOF = false;

      if (track.pdgCode() == pdgProton) {
        registry.fill(HIST("hReco_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt());
        registry.fill(HIST("hReco_EtaPhiPtMC_Proton"), track.eta_MC(), track.phi_MC(), track.pt_MC());
        registry.fill(HIST("hResPt_Proton"), track.pt_MC(), track.pt() - track.pt_MC());
        registry.fill(HIST("hResEta_Proton"), track.eta_MC(), track.eta() - track.eta_MC());
        registry.fill(HIST("hResPhi_Proton"), track.phi_MC(), track.phi() - track.phi_MC());

        if (TMath::Abs(track.tpcNSigmaPr()) < nsigmaTPC) {
          if (track.pt() < pTthrpr_TOF) {
            isPr = true;
            registry.fill(HIST("hReco_PID_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt());
          } else if (TMath::Abs(track.tofNSigmaPr()) < nsigmaTOF) {
            isPr = true;
            registry.fill(HIST("hReco_PID_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt());
          }
        }
        registry.fill(HIST("hnSigmaTPCVsPt_Pr_MC"), track.pt(), track.tpcNSigmaPr());
        registry.fill(HIST("hnSigmaTOFVsPt_Pr_MC"), track.pt(), track.tofNSigmaPr());
      }
      if (track.pdgCode() == -pdgProton) {
        registry.fill(HIST("hReco_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        registry.fill(HIST("hReco_EtaPhiPtMC_Proton"), track.eta_MC(), track.phi_MC(), track.pt_MC() * -1);
        registry.fill(HIST("hResPt_AntiProton"), track.pt_MC(), track.pt() - track.pt_MC());
        registry.fill(HIST("hResEta_AntiProton"), track.eta_MC(), track.eta() - track.eta_MC());
        registry.fill(HIST("hResPhi_AntiProton"), track.phi_MC(), track.phi() - track.phi_MC());

        if (TMath::Abs(track.tpcNSigmaPr()) < nsigmaTPC) {
          if (track.pt() < pTthrpr_TOF) {
            isAntiPr = true;
            registry.fill(HIST("hReco_PID_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
          } else if (TMath::Abs(track.tofNSigmaPr()) < nsigmaTOF) {
            isAntiPr = true;
            registry.fill(HIST("hReco_PID_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
          }
        }
        registry.fill(HIST("hnSigmaTPCVsPt_Pr_MC"), track.pt() * -1, track.tpcNSigmaPr());
        registry.fill(HIST("hnSigmaTOFVsPt_Pr_MC"), track.pt() * -1, track.tofNSigmaPr());
      }
      if (track.pdgCode() == pdgDeuteron) {
        registry.fill(HIST("hReco_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt());
        registry.fill(HIST("hReco_EtaPhiPtMC_Deuteron"), track.eta_MC(), track.phi_MC(), track.pt_MC());
        registry.fill(HIST("hResPt_Deuteron"), track.pt_MC(), track.pt() - track.pt_MC());
        registry.fill(HIST("hResEta_Deuteron"), track.eta_MC(), track.eta() - track.eta_MC());
        registry.fill(HIST("hResPhi_Deuteron"), track.phi_MC(), track.phi() - track.phi_MC());

        if (TMath::Abs(track.tpcNSigmaDe()) < nsigmaTPC && TMath::Abs(track.tpcNSigmaPr()) >= nsigmaTPC) {
          if (TMath::Abs(track.tofNSigmaDe()) < nsigmaTOF && TMath::Abs(track.tofNSigmaPr()) >= nsigmaTOF) {
            registry.fill(HIST("hReco_PID_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt());
          }
        }
        registry.fill(HIST("hnSigmaTPCVsPt_De_MC"), track.pt(), track.tpcNSigmaDe());
        registry.fill(HIST("hnSigmaTOFVsPt_De_MC"), track.pt(), track.tofNSigmaDe());
      }
      if (track.pdgCode() == -pdgDeuteron) {
        registry.fill(HIST("hReco_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt() * -1);
        registry.fill(HIST("hReco_EtaPhiPtMC_Deuteron"), track.eta_MC(), track.phi_MC(), track.pt_MC() * -1);
        registry.fill(HIST("hResPt_AntiDeuteron"), track.pt_MC(), track.pt() - track.pt_MC());
        registry.fill(HIST("hResEta_AntiDeuteron"), track.eta_MC(), track.eta() - track.eta_MC());
        registry.fill(HIST("hResPhi_AntiDeuteron"), track.phi_MC(), track.phi() - track.phi_MC());

        if (TMath::Abs(track.tpcNSigmaDe()) < nsigmaTPC && TMath::Abs(track.tpcNSigmaPr()) >= nsigmaTPC) {
          if (TMath::Abs(track.tofNSigmaDe()) < nsigmaTOF && TMath::Abs(track.tofNSigmaPr()) >= nsigmaTOF) {
            isAntiDeTPCTOF = true;
            registry.fill(HIST("hReco_PID_EtaPhiPt_Deuteron"), track.eta(), track.phi(), track.pt() * -1);
          }
        }
        registry.fill(HIST("hnSigmaTPCVsPt_De_MC"), track.pt() * -1, track.tpcNSigmaDe());
        registry.fill(HIST("hnSigmaTOFVsPt_De_MC"), track.pt() * -1, track.tofNSigmaDe());
      }

      if (!mcCorrelation) {
        continue;
      }

      if (isAntiDeTPCTOF) {
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
      if (TMath::Abs(collision.posZ()) > cutzvertex)
        continue;
      registry.fill(HIST("hNEvents"), 0.5);

      int vertexBinToMix = std::floor((collision.posZ() + cutzvertex) / (2 * cutzvertex / _vertexNbinsToMix));
      int centBinToMix = std::floor(collision.multPerc() / (100.0 / _multNsubBins));

      if (selectedtracksMC_p.find(collision.globalIndex()) != selectedtracksMC_p.end()) {
        mixbins_pantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }

      if (selectedtracksMC_antid.find(collision.globalIndex()) != selectedtracksMC_antid.end()) {
        mixbins_antidantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }

      if (selectedtracksPIDMC_p.find(collision.globalIndex()) != selectedtracksPIDMC_p.end()) {
        mixbinsPID_pantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }

      if (selectedtracksPIDMC_antid.find(collision.globalIndex()) != selectedtracksPIDMC_antid.end()) {
        mixbinsPID_antidantip[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
      }
    } // coll

    if (!mixbins_pantip.empty()) {

      for (auto i = mixbins_pantip.begin(); i != mixbins_pantip.end(); i++) { // iterating over all vertex&mult bins

        std::vector<colType> value = i->second;
        int EvPerBin = value.size(); // number of collisions in each vertex&mult bin

        for (int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = value[indx1];

          if (selectedtracksMC_antip.find(col1->index()) != selectedtracksMC_antip.end()) {
            mixTracks<0, 1 /*MC qa*/>(selectedtracksMC_p[col1->index()], selectedtracksMC_antip[col1->index()], 0, 0); // mixing SE
          }

          for (int indx2 = indx1 + 1; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracksMC_antip.find(col2->index()) != selectedtracksMC_antip.end()) {
              mixTracks<1, 1>(selectedtracksMC_p[col1->index()], selectedtracksMC_antip[col2->index()], 0, 0); // mixing ME
            }
          }
        } // event
      }
    }

    if (!mixbinsPID_pantip.empty()) {

      for (auto i = mixbinsPID_pantip.begin(); i != mixbinsPID_pantip.end(); i++) { // iterating over all vertex&mult bins

        std::vector<colType> value = i->second;
        int EvPerBin = value.size(); // number of collisions in each vertex&mult bin

        for (int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = value[indx1];

          if (selectedtracksPIDMC_antip.find(col1->index()) != selectedtracksPIDMC_antip.end()) {
            mixTracks<0, 1 /*MC qa*/>(selectedtracksPIDMC_p[col1->index()], selectedtracksPIDMC_antip[col1->index()], 0, 1); // mixing SE
          }

          for (int indx2 = indx1 + 1; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracksPIDMC_antip.find(col2->index()) != selectedtracksPIDMC_antip.end()) {
              mixTracks<1, 1>(selectedtracksPIDMC_p[col1->index()], selectedtracksPIDMC_antip[col2->index()], 0, 1); // mixing ME
            }
          }
        }
      }
    }

    if (!mixbins_antidantip.empty()) {

      for (auto i = mixbins_antidantip.begin(); i != mixbins_antidantip.end(); i++) { // iterating over all vertex&mult bins

        std::vector<colType> value = i->second;
        int EvPerBin = value.size(); // number of collisions in each vertex&mult bin

        for (int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = value[indx1];

          if (selectedtracksMC_antip.find(col1->index()) != selectedtracksMC_antip.end()) {
            mixTracks<0, 1 /*MC qa*/>(selectedtracksMC_antid[col1->index()], selectedtracksMC_antip[col1->index()], 1, 0); // mixing SE
          }

          for (int indx2 = indx1 + 1; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracksMC_antip.find(col2->index()) != selectedtracksMC_antip.end()) {
              mixTracks<1, 1>(selectedtracksMC_antid[col1->index()], selectedtracksMC_antip[col2->index()], 1, 0); // mixing ME
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
            mixTracks<0, 1 /*MC qa*/>(selectedtracksPIDMC_antid[col1->index()], selectedtracksPIDMC_antip[col1->index()], 1, 1); // mixing SE
          }

          for (int indx2 = indx1 + 1; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            if (col1 == col2) {
              continue;
            }

            if (selectedtracksPIDMC_antip.find(col2->index()) != selectedtracksPIDMC_antip.end()) {
              mixTracks<1, 1>(selectedtracksPIDMC_antid[col1->index()], selectedtracksPIDMC_antip[col2->index()], 1, 1); // mixing ME
            }
          }
        }
      }
    }

    // clearing up
    for (auto i = selectedtracksMC_antid.begin(); i != selectedtracksMC_antid.end(); i++)
      (i->second).clear();
    selectedtracksMC_antid.clear();

    for (auto i = selectedtracksMC_antip.begin(); i != selectedtracksMC_antip.end(); i++)
      (i->second).clear();
    selectedtracksMC_antip.clear();

    for (auto i = selectedtracksMC_p.begin(); i != selectedtracksMC_p.end(); i++)
      (i->second).clear();
    selectedtracksMC_p.clear();

    for (auto i = selectedtracksPIDMC_antid.begin(); i != selectedtracksPIDMC_antid.end(); i++)
      (i->second).clear();
    selectedtracksPIDMC_antid.clear();

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

    for (auto& pair : mixbinsPID_pantip) {
      pair.second.clear(); // clear the vector associated with the key
    }
    mixbinsPID_pantip.clear(); // clear the map

    for (auto& pair : mixbins_pantip) {
      pair.second.clear(); // clear the vector associated with the key
    }
    mixbins_pantip.clear(); // clear the map
  }
  PROCESS_SWITCH(hadronnucleicorrelation, processMC, "processMC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<hadronnucleicorrelation>(cfgc)};
}
