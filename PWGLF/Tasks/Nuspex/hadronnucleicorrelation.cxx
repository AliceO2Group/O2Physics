// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for drtails of the copyright holders.
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

#include "PWGCF/Femto3D/Core/femto3dPairTask.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/Track.h"

#include "TGrid.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TParameter.h>
#include <TVector2.h>
#include <TVector3.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum Modes {
  kDbarPbar = 0,
  kDP,
  kDbarP,
  kDPbar,
  kPbarP,
  kPbarPbar,
  kPP,
  kPPbar
};

struct HadronNucleiCorrelation {

  static constexpr int betahasTOFthr = -100;

  SliceCache cache;

  Configurable<int> mode{"mode", 0, "0: antid-antip, 1: d-p, 2: antid-p, 3: d-antip, 4: antip-p, 5: antip-antip, 6: p-p, 7: p-antip"};

  Configurable<bool> dorapidity{"dorapidity", false, "do rapidity dependent analysis"};
  Configurable<bool> doQA{"doQA", true, "save QA histograms"};
  Configurable<bool> doMCQA{"doMCQA", false, "save MC QA histograms"};
  Configurable<bool> isMC{"isMC", false, "is MC"};
  Configurable<bool> isMCGen{"isMCGen", false, "is isMCGen"};
  Configurable<bool> isPrim{"isPrim", true, "is isPrim"};
  Configurable<bool> docorrection{"docorrection", false, "do efficiency correction"};

  Configurable<std::string> fCorrectionPath{"fCorrectionPath", "", "Correction path to file"};
  Configurable<std::string> fCorrectionHisto{"fCorrectionHisto", "", "Correction histogram"};
  Configurable<std::string> cfgUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Event selection
  Configurable<float> cutzVertex{"cutzVertex", 10.0, "|vertexZ| value limit"};

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
  Configurable<float> nsigmaITSPr{"nsigmaITSPr", -2.0f, "cut nsigma ITS Pr"};
  Configurable<float> nsigmaITSDe{"nsigmaITSDe", -2.0f, "cut nsigma ITS De"};
  Configurable<bool> doITSPID{"doITSPID", true, "do ITS PID"};
  Configurable<float> pTthrpr_TOF{"pTthrpr_TOF", 0.8f, "threshold pT proton to use TOF"};
  Configurable<float> pTthrpr_TPCEl{"pTthrpr_TPCEl", 1.0f, "threshold pT proton to use TPC El rejection"};
  Configurable<float> pTthrde_TOF{"pTthrde_TOF", 1.0f, "threshold pT deuteron to use TOF"};
  Configurable<float> pTthrde_TPCEl{"pTthrde_TPCEl", 1.0f, "threshold pT deuteron to use TPC El rejection"};
  Configurable<bool> rejectionEl{"rejectionEl", true, "use TPC El rejection"};
  Configurable<float> max_tpcSharedCls{"max_tpcSharedCls", 0.4, "maximum fraction of TPC shared clasters"};
  Configurable<int> min_itsNCls{"min_itsNCls", 0, "minimum allowed number of ITS clasters"};
  Configurable<int> maxmixcollsGen{"maxmixcollsGen", 100, "maxmixcollsGen"};
  Configurable<float> radiusTPC{"radiusTPC", 1.2, "TPC radius to calculate phi_star for"};
  Configurable<float> dEta{"dEta", 0.01, "minimum allowed difference in eta between two tracks in a pair"};
  Configurable<float> dPhi{"dPhi", 0.01, "minimum allowed difference in phi_star between two tracks in a pair"};

  // Mixing parameters
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 50.0f, 100.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ColumnBinningPolicy<aod::singletrackselector::PosZ, aod::singletrackselector::Mult> colBinning{{confVtxBins, confMultBins}, true};
  ColumnBinningPolicy<aod::mccollision::PosZ, o2::aod::mult::MultMCNParticlesEta05> colBinningGen{{confVtxBins, confMultBins}, true};

  // pT/A bins
  Configurable<std::vector<double>> pTBins{"pTBins", {0.6f, 1.0f, 1.2f, 2.f}, "p_{T} bins"};

  ConfigurableAxis AxisNSigma{"AxisNSigma", {35, -7.f, 7.f}, "n#sigma"};
  ConfigurableAxis DeltaPhiAxis = {"DeltaPhiAxis", {46, -1 * o2::constants::math::PIHalf, 3 * o2::constants::math::PIHalf}, "#Delta#phi (rad)"};

  using FilteredCollisions = soa::Filtered<aod::SingleCollSels>;
  using SimCollisions = soa::Join<aod::McCollisions, aod::MultsExtraMC>;
  using SimParticles = aod::McParticles;
  using FilteredTracks = soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkExtras, aod::SinglePIDEls, aod::SinglePIDPrs, aod::SinglePIDDes>>;                      // new tables (v3)
  using FilteredTracksMC = soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkMCs, aod::SingleTrkExtras, aod::SinglePIDEls, aod::SinglePIDPrs, aod::SinglePIDDes>>; // new tables (v3)

  HistogramRegistry registry{"registry"};
  HistogramRegistry QA{"QA"};

  using trkType = const FilteredTracks::iterator*;
  // using trkTypeMC = const FilteredTracksMC::iterator*;
  // typedef std::shared_ptr<FilteredCollisions::iterator> colType;
  // typedef std::shared_ptr<SimCollisions::iterator> MCcolType;

  std::unique_ptr<o2::aod::singletrackselector::FemtoPair<trkType>> Pair = std::make_unique<o2::aod::singletrackselector::FemtoPair<trkType>>();
  // std::unique_ptr<o2::aod::singletrackselector::FemtoPair<trkTypeMC>> PairMC = std::make_unique<o2::aod::singletrackselector::FemtoPair<trkTypeMC>>();

  // Data histograms
  std::vector<std::shared_ptr<TH3>> hEtaPhi_SE;
  std::vector<std::shared_ptr<TH3>> hEtaPhi_ME;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhi_SE;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhi_ME;

  int nBinspT;
  TH2F* hEffpTEta_proton;
  TH2F* hEffpTEta_antiproton;
  TH2F* hEffpTEta_deuteron;
  TH2F* hEffpTEta_antideuteron;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(cfgUrl.value);
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
    AxisSpec phiAxis = {157, 0., o2::constants::math::TwoPI, "#phi (rad)"};
    AxisSpec pTAxis = {200, -10.f, 10.f, "p_{T} GeV/c"};
    AxisSpec pTAxis_small = {100, -5.f, 5.f, "p_{T} GeV/c"};

    AxisSpec DeltaEtaAxis = {300, -1.5, 1.5, "#Delta#eta"};
    AxisSpec DeltaRapAxis = {300, -1.5, 1.5, "#Delta y"};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1D, {{7, 0.f, 7.f}}});
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "Selected");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "Mixing");

    registry.add("hNtrig_total", "hNtrig_total", {HistType::kTH1D, {ptBinnedAxis}});

    nBinspT = pTBins.value.size() - 1;

    TString name = "AntiDeAntiPr";
    switch (mode) {
      case kDP:
        name = "DePr";
        break;
      case kDbarP:
        name = "AntiDePr";
        break;
      case kDPbar:
        name = "DeAntiPr";
        break;
      case kPbarP:
        name = "AntiPrPr";
        break;
      case kPbarPbar:
        name = "AntiPrAntiPr";
        break;
      case kPP:
        name = "PrPr";
        break;
      case kPPbar:
        name = "PrAntiPr";
        break;
    }

    if (!isMC) {
      for (int i = 0; i < nBinspT; i++) {

        if (dorapidity) {

          auto htempSE_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_%s_SE_pt%02.0f%02.0f", name.Data(), pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta y #Delta#phi (%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaRapAxis, DeltaPhiAxis, ptBinnedAxis}});
          auto htempME_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_%s_ME_pt%02.0f%02.0f", name.Data(), pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta y #Delta#phi (%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaRapAxis, DeltaPhiAxis, ptBinnedAxis}});
          hEtaPhi_SE.push_back(std::move(htempSE_AntiDeAntiPr));
          hEtaPhi_ME.push_back(std::move(htempME_AntiDeAntiPr));

          auto hCorrtempSE_AntiDeAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_%s_SE_pt%02.0f%02.0f", name.Data(), pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta y #Delta#phi (%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaRapAxis, DeltaPhiAxis, ptBinnedAxis}});
          auto hCorrtempME_AntiDeAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_%s_ME_pt%02.0f%02.0f", name.Data(), pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta y #Delta#phi (%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaRapAxis, DeltaPhiAxis, ptBinnedAxis}});
          hCorrEtaPhi_SE.push_back(std::move(hCorrtempSE_AntiDeAntiPr));
          hCorrEtaPhi_ME.push_back(std::move(hCorrtempME_AntiDeAntiPr));
        } else {

          auto htempSE_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_%s_SE_pt%02.0f%02.0f", name.Data(), pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta#eta#Delta#phi (%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
          auto htempME_AntiDeAntiPr = registry.add<TH3>(Form("hEtaPhi_%s_ME_pt%02.0f%02.0f", name.Data(), pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("Raw #Delta#eta#Delta#phi (%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
          hEtaPhi_SE.push_back(std::move(htempSE_AntiDeAntiPr));
          hEtaPhi_ME.push_back(std::move(htempME_AntiDeAntiPr));

          auto hCorrtempSE_AntiDeAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_%s_SE_pt%02.0f%02.0f", name.Data(), pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
          auto hCorrtempME_AntiDeAntiPr = registry.add<TH3>(Form("hCorrEtaPhi_%s_ME_pt%02.0f%02.0f", name.Data(), pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10), Form("#Delta#eta#Delta#phi (%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1)), {HistType::kTH3F, {DeltaEtaAxis, DeltaPhiAxis, ptBinnedAxis}});
          hCorrEtaPhi_SE.push_back(std::move(hCorrtempSE_AntiDeAntiPr));
          hCorrEtaPhi_ME.push_back(std::move(hCorrtempME_AntiDeAntiPr));
        }
      }
    }

    registry.add("hPrDCAxy", "DCAxy p", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
    registry.add("hAntiPrDCAxy", "DCAxy #bar{p}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
    registry.add("hDeDCAxy", "DCAxy d", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
    registry.add("hAntiDeDCAxy", "DCAxy #bar{d}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
    registry.add("hMult", "multiplicity", {HistType::kTH1D, {{200, 0.f, 200.f, "N_{ch}"}}});

    if (doQA) {
      // Track QA
      QA.add("QA/hVtxZ_trk", "#it{z}_{vtx}", {HistType::kTH1D, {{150, -15.f, 15.f, "#it{z}_{vtx} (cm)"}}});
      QA.add("QA/hTPCnClusters", "N TPC Clusters; N TPC Clusters", {HistType::kTH1D, {{200, 0.f, 200.f}}});
      QA.add("QA/hTPCSharedClusters", "N TPC Shared Clusters; N TPC SharedClusters", {HistType::kTH1D, {{100, 0.f, 1.f}}});
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
      QA.add("QA/hnSigmaITSVsPt_Pr", "n#sigma ITS vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma ITS", {HistType::kTH2D, {pTAxis, AxisNSigma}});
      QA.add("QA/hnSigmaITSVsPt_De", "n#sigma ITS vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma ITS", {HistType::kTH2D, {pTAxis, AxisNSigma}});
      QA.add("QA/hdEtadPhistar", ";dPhi*;dEta ", {HistType::kTH2D, {{101, -0.2, 0.2, "dPhi*"}, {101, -0.2, 0.2, "dEta"}}});

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
        QA.add("QA/hnSigmaITSVsPt_Pr_AfterSel", "n#sigma ITS vs p_{T} for p hypothesis (all tracks); p_{T} (GeV/c); n#sigma ITS", {HistType::kTH2D, {pTAxis, AxisNSigma}});
        QA.add("QA/hnSigmaITSVsPt_De_AfterSel", "n#sigma ITS vs p_{T} for d hypothesis (all tracks); p_{T} (GeV/c); n#sigma ITS", {HistType::kTH2D, {pTAxis, AxisNSigma}});
      }
    }

    if (isMC) {
      registry.add("hPrimPrDCAxy", "DCAxy p", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hPrimAntiPrDCAxy", "DCAxy #bar{p}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hPrimDeDCAxy", "DCAxy d", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hPrimAntiDeDCAxy", "DCAxy #bar{d}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hSecMatPrDCAxy", "DCAxy p", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hSecMatAntiPrDCAxy", "DCAxy #bar{p}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hSecMatDeDCAxy", "DCAxy d", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hSecMatAntiDeDCAxy", "DCAxy #bar{d}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hSecWeakPrDCAxy", "DCAxy p", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hSecWeakAntiPrDCAxy", "DCAxy #bar{p}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hSecWeakDeDCAxy", "DCAxy d", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hSecWeakAntiDeDCAxy", "DCAxy #bar{d}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});

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
        registry.add("hResPhi_Proton", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, o2::constants::math::TwoPI, "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResPhi_Deuteron", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, o2::constants::math::TwoPI, "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResEta_AntiProton", "; #eta(gen); #eta(reco) - #eta(gen)  ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResEta_AntiDeuteron", "; #eta(gen); #eta(reco) - #eta(gen) ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResPhi_AntiProton", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, o2::constants::math::TwoPI, "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResPhi_AntiDeuteron", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, o2::constants::math::TwoPI, "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});

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

      registry.add("hGen_EtaPhiPt_Proton", "Gen (anti)protons in gen collisions", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});
      registry.add("hGen_EtaPhiPt_Deuteron", "Gen (anti)deuteron in gen collisions", {HistType::kTH3F, {etaAxis, phiAxis, pTAxis_small}});

      registry.add("Generated/hQAProtons", "hQAProtons", {HistType::kTH1D, {{5, 0.f, 5.f}}});
      registry.get<TH1>(HIST("Generated/hQAProtons"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("Generated/hQAProtons"))->GetXaxis()->SetBinLabel(2, "PhysicalPrimary");
      registry.get<TH1>(HIST("Generated/hQAProtons"))->GetXaxis()->SetBinLabel(3, "|#eta|<0.8");
      registry.get<TH1>(HIST("Generated/hQAProtons"))->GetXaxis()->SetBinLabel(4, "no daughters");
      registry.get<TH1>(HIST("Generated/hQAProtons"))->GetXaxis()->SetBinLabel(5, "d daughter");

      registry.add("Generated/hQADeuterons", "hQADeuterons", {HistType::kTH1D, {{3, 0.f, 3.f}}});
      registry.get<TH1>(HIST("Generated/hQADeuterons"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("Generated/hQADeuterons"))->GetXaxis()->SetBinLabel(2, "PhysicalPrimary");
      registry.get<TH1>(HIST("Generated/hQADeuterons"))->GetXaxis()->SetBinLabel(3, "|#eta|<0.8");

      registry.add("Generated/hDeuteronsVsPt", "hDeuteronsVsPt;  p_{T} (GeV/c);", {HistType::kTH1D, {{100, 0.f, 10.f}}});
      registry.add("Generated/hAntiDeuteronsVsPt", "hAntiDeuteronsVsPt;  p_{T} (GeV/c);", {HistType::kTH1D, {{100, 0.f, 10.f}}});
    }
  }

  // Filters
  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) <= cutzVertex;
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
    bool isTPCPID = std::abs(track.tpcNSigmaPr()) < nsigmaTPC;
    bool isTOFPID = std::abs(track.tofNSigmaPr()) < nsigmaTOF;
    bool isTPCElRejection = rejectionEl && track.beta() < betahasTOFthr && track.pt() < pTthrpr_TPCEl && track.tpcNSigmaEl() >= nsigmaElPr;
    bool isITSPID = track.itsNSigmaPr() > nsigmaITSPr;

    if (isTPCPID) {
      if (track.pt() < pTthrpr_TOF) {
        if (!doITSPID || isITSPID) {
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
      } else if (isTPCElRejection) {
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
      } else if (isTOFPID) {
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
    bool isTPCPID = std::abs(track.tpcNSigmaDe()) < nsigmaTPC;
    bool isTOFPID = std::abs(track.tofNSigmaDe()) < nsigmaTOF;
    bool isTPCElRejection = rejectionEl && track.beta() < betahasTOFthr && track.pt() < pTthrde_TPCEl && track.tpcNSigmaEl() >= nsigmaElDe;
    bool isITSPID = track.itsNSigmaDe() > nsigmaITSDe;

    if (isTPCPID) {
      if (track.pt() < pTthrde_TOF) {
        if (!doITSPID || isITSPID) {
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
      } else if (isTPCElRejection) {
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
      } else if (isTOFPID) {
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
    if (std::abs(track.dcaXY()) > (par0 + par1 / track.pt()))
      passcut = false;
    if (std::abs(track.dcaZ()) > (par0 + par1 / track.pt()))
      passcut = false;

    return passcut;
  }

  template <typename T1>
  void fillHistograms(T1 const& part0, T1 const& part1, bool ME, bool isIdentical)
  {
    Pair->SetPair(&part0, &part1);
    Pair->SetIdentical(isIdentical);
    if (isIdentical && Pair->IsClosePair(dEta, dPhi, radiusTPC)) {
      QA.fill(HIST("QA/hdEtadPhistar"), Pair->GetPhiStarDiff(radiusTPC), Pair->GetEtaDiff());
      return;
    }

    float deltaEta = part0.eta() - part1.eta();
    float deltaPhi = part0.phi() - part1.phi();
    deltaPhi = RecoDecay::constrainAngle(deltaPhi, -1 * o2::constants::math::PIHalf);

    for (int k = 0; k < nBinspT; k++) {

      if (part0.pt() >= pTBins.value.at(k) && part0.pt() < pTBins.value.at(k + 1)) {

        float corr0 = 1, corr1 = 1;

        if (docorrection) { // Apply corrections
          switch (mode) {
            case 0:
              corr0 = hEffpTEta_antideuteron->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffpTEta_antiproton->Interpolate(part1.pt(), part1.eta());
              break;
            case 1:
              corr0 = hEffpTEta_deuteron->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffpTEta_proton->Interpolate(part1.pt(), part1.eta());
              break;
            case 2:
              corr0 = hEffpTEta_antideuteron->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffpTEta_proton->Interpolate(part1.pt(), part1.eta());
              break;
            case 3:
              corr0 = hEffpTEta_deuteron->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffpTEta_antiproton->Interpolate(part1.pt(), part1.eta());
              break;
            case 4:
              corr0 = hEffpTEta_antiproton->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffpTEta_proton->Interpolate(part1.pt(), part1.eta());
              break;
            case 5:
              corr0 = hEffpTEta_antiproton->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffpTEta_antiproton->Interpolate(part1.pt(), part1.eta());
              break;
            case 6:
              corr0 = hEffpTEta_proton->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffpTEta_proton->Interpolate(part1.pt(), part1.eta());
              break;
            case 7:
              corr0 = hEffpTEta_proton->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffpTEta_antiproton->Interpolate(part1.pt(), part1.eta());
              break;
          }
        }

        if (ME) {
          hEtaPhi_ME[k]->Fill(deltaEta, deltaPhi, part1.pt());
          hCorrEtaPhi_ME[k]->Fill(deltaEta, deltaPhi, part1.pt(), 1. / (corr0 * corr1));
        } else {
          hEtaPhi_SE[k]->Fill(deltaEta, deltaPhi, part1.pt());
          hCorrEtaPhi_SE[k]->Fill(deltaEta, deltaPhi, part1.pt(), 1. / (corr0 * corr1));
        } // SE
      } // pT condition
    } // nBinspT loop

    Pair->ResetPair();
  }

  template <typename T1>
  void fillHistogramsGen(T1 const& part0, T1 const& part1, bool ME)
  {

    float deltaEta = part0.eta() - part1.eta();
    float deltaPhi = part0.phi() - part1.phi();
    deltaPhi = RecoDecay::constrainAngle(deltaPhi, -1 * o2::constants::math::PIHalf);

    for (int k = 0; k < nBinspT; k++) {

      if (part0.pt() >= pTBins.value.at(k) && part0.pt() < pTBins.value.at(k + 1)) {

        if (ME) {
          hEtaPhi_ME[k]->Fill(deltaEta, deltaPhi, part1.pt());
          hCorrEtaPhi_ME[k]->Fill(deltaEta, deltaPhi, part1.pt());
        } else {
          hEtaPhi_SE[k]->Fill(deltaEta, deltaPhi, part1.pt());
          hCorrEtaPhi_SE[k]->Fill(deltaEta, deltaPhi, part1.pt());
        } // SE
      } // pT condition
    } // nBinspT loop
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

  void processSameEvent(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks)
  {

    if (std::abs(collision.posZ()) > cutzVertex)
      return;

    registry.fill(HIST("hNEvents"), 0.5);
    registry.fill(HIST("hMult"), collision.mult());

    for (const auto& track : tracks) {
      if (track.tpcFractionSharedCls() > max_tpcSharedCls)
        continue;
      if (track.itsNCls() < min_itsNCls)
        continue;

      if (IsProton(track, +1))
        registry.fill(HIST("hPrDCAxy"), track.dcaXY(), track.pt());
      if (IsProton(track, -1))
        registry.fill(HIST("hAntiPrDCAxy"), track.dcaXY(), track.pt());
      if (IsDeuteron(track, +1))
        registry.fill(HIST("hDeDCAxy"), track.dcaXY(), track.pt());
      if (IsDeuteron(track, -1))
        registry.fill(HIST("hAntiDeDCAxy"), track.dcaXY(), track.pt());

      if (!applyDCAcut(track))
        continue;

      if (doQA) {
        QA.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        QA.fill(HIST("QA/hTPCSharedClusters"), track.tpcFractionSharedCls());
        QA.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        QA.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        QA.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        QA.fill(HIST("QA/hDCAxy"), track.dcaXY(), track.pt());
        QA.fill(HIST("QA/hDCAz"), track.dcaZ(), track.pt());
        QA.fill(HIST("QA/TPCChi2VsPZ"), track.tpcInnerParam() / track.sign(), track.tpcChi2NCl());
        QA.fill(HIST("QA/hVtxZ_trk"), collision.posZ());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_El"), track.pt() * track.sign(), track.tpcNSigmaEl());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        QA.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
        QA.fill(HIST("QA/hnSigmaITSVsPt_Pr"), track.pt() * track.sign(), track.itsNSigmaPr());
        QA.fill(HIST("QA/hnSigmaITSVsPt_De"), track.pt() * track.sign(), track.itsNSigmaDe());

        if (IsProton(track, +1)) {
          QA.fill(HIST("QA/hEtaAntiPr"), track.eta());
          QA.fill(HIST("QA/hPhiAntiPr"), track.phi());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaITSVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.itsNSigmaPr());
        }
        if (IsProton(track, -1)) {
          QA.fill(HIST("QA/hEtaPr"), track.eta());
          QA.fill(HIST("QA/hPhiPr"), track.phi());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          QA.fill(HIST("QA/hnSigmaITSVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.itsNSigmaPr());
        }
        if (IsDeuteron(track, +1)) {
          QA.fill(HIST("QA/hEtaAntiDe"), track.eta());
          QA.fill(HIST("QA/hPhiAntiDe"), track.phi());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
          QA.fill(HIST("QA/hnSigmaITSVsPt_De_AfterSel"), track.pt() * track.sign(), track.itsNSigmaDe());
        }
        if (IsDeuteron(track, -1)) {
          QA.fill(HIST("QA/hEtaDe"), track.eta());
          QA.fill(HIST("QA/hPhiDe"), track.phi());
          QA.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          QA.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
          QA.fill(HIST("QA/hnSigmaITSVsPt_De_AfterSel"), track.pt() * track.sign(), track.itsNSigmaDe());
        }
      }
    }

    Pair->SetMagField1(collision.magField());
    Pair->SetMagField2(collision.magField());

    if (mode == kPbarPbar || mode == kPP) { // Identical particle combinations

      for (const auto& [part0, part1] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {

        if (part0.tpcFractionSharedCls() > max_tpcSharedCls)
          continue;
        if (part0.itsNCls() < min_itsNCls)
          continue;
        if (part1.tpcFractionSharedCls() > max_tpcSharedCls)
          continue;
        if (part1.itsNCls() < min_itsNCls)
          continue;

        if (!applyDCAcut(part0))
          continue;
        if (!applyDCAcut(part1))
          continue;

        // mode 6
        if (mode == kPP) {
          if (!IsProton(part0, +1))
            continue;
          if (!IsProton(part1, +1))
            continue;
        }
        // mode 5
        if (mode == kPbarPbar) {
          if (!IsProton(part0, -1))
            continue;
          if (!IsProton(part1, -1))
            continue;
        }

        fillHistograms(part0, part1, false, true);
      }

    } else {

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

        if (part0.tpcFractionSharedCls() > max_tpcSharedCls)
          continue;
        if (part0.itsNCls() < min_itsNCls)
          continue;
        if (part1.tpcFractionSharedCls() > max_tpcSharedCls)
          continue;
        if (part1.itsNCls() < min_itsNCls)
          continue;

        if (!applyDCAcut(part0))
          continue;
        if (!applyDCAcut(part1))
          continue;

        // modes 0,1,2,3,4,7
        if (mode == kDbarPbar) {
          if (!IsDeuteron(part0, -1))
            continue;
          if (!IsProton(part1, -1))
            continue;
        }
        if (mode == kDP) {
          if (!IsDeuteron(part0, +1))
            continue;
          if (!IsProton(part1, +1))
            continue;
        }
        if (mode == kDbarP) {
          if (!IsDeuteron(part0, -1))
            continue;
          if (!IsProton(part1, +1))
            continue;
        }
        if (mode == kDPbar) {
          if (!IsDeuteron(part0, +1))
            continue;
          if (!IsProton(part1, -1))
            continue;
        }
        if (mode == kPbarP) {
          if (!IsProton(part0, -1))
            continue;
          if (!IsProton(part1, +1))
            continue;
        }
        if (mode == kPPbar) {
          if (!IsProton(part0, +1))
            continue;
          if (!IsProton(part1, -1))
            continue;
        }

        fillHistograms(part0, part1, false, false);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processSameEvent, "processSameEvent", true);

  void processMixedEvent(FilteredCollisions const& collisions, FilteredTracks const& tracks)
  {

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, collisions, collisions)) {

      // LOGF(info, "Mixed event collisions: (%d, %d) zvtx (%.1f, %.1f) mult (%d, %d)", collision1.globalIndex(), collision2.globalIndex(), collision1.posZ(), collision2.posZ(), collision1.mult(), collision2.mult());

      auto groupPartsOne = tracks.sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision1.globalIndex(), cache);
      auto groupPartsTwo = tracks.sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      Pair->SetMagField1(magFieldTesla1);
      Pair->SetMagField2(magFieldTesla2);

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        if (part0.tpcFractionSharedCls() > max_tpcSharedCls)
          continue;
        if (part0.itsNCls() < min_itsNCls)
          continue;
        if (part1.tpcFractionSharedCls() > max_tpcSharedCls)
          continue;
        if (part1.itsNCls() < min_itsNCls)
          continue;

        if (!applyDCAcut(part0))
          continue;
        if (!applyDCAcut(part1))
          continue;

        //{"mode", 0, "0: antid-antip, 1: d-p, 2: antid-p, 3: d-antip, 4: antip-p, 5: antip-antip, 6: p-p, 7: p-antip"};
        if (mode == kDbarPbar) {
          if (!IsDeuteron(part0, -1))
            continue;
          if (!IsProton(part1, -1))
            continue;
        }
        if (mode == kDP) {
          if (!IsDeuteron(part0, +1))
            continue;
          if (!IsProton(part1, +1))
            continue;
        }
        if (mode == kDbarP) {
          if (!IsDeuteron(part0, -1))
            continue;
          if (!IsProton(part1, +1))
            continue;
        }
        if (mode == kDPbar) {
          if (!IsDeuteron(part0, +1))
            continue;
          if (!IsProton(part1, -1))
            continue;
        }
        if (mode == kPbarP) {
          if (!IsProton(part0, -1))
            continue;
          if (!IsProton(part1, +1))
            continue;
        }
        if (mode == kPbarPbar) {
          if (!IsProton(part0, -1))
            continue;
          if (!IsProton(part1, -1))
            continue;
        }
        if (mode == kPP) {
          if (!IsProton(part0, +1))
            continue;
          if (!IsProton(part1, +1))
            continue;
        }
        if (mode == kPPbar) {
          if (!IsProton(part0, +1))
            continue;
          if (!IsProton(part1, -1))
            continue;
        }

        bool isIdentical = false;
        if (mode == kPbarPbar || mode == kPP)
          isIdentical = true;

        fillHistograms(part0, part1, true, isIdentical);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processMixedEvent, "processMixedEvent", true);

  void processMC(FilteredCollisions const&, FilteredTracksMC const& tracks)
  {
    for (const auto& track : tracks) {
      if (std::abs(track.template singleCollSel_as<FilteredCollisions>().posZ()) > cutzVertex)
        continue;

      if (track.tpcFractionSharedCls() > max_tpcSharedCls)
        continue;
      if (track.itsNCls() < min_itsNCls)
        continue;

      if (IsProton(track, +1) && track.pdgCode() == PDG_t::kProton) {
        registry.fill(HIST("hPrDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 0)
          registry.fill(HIST("hPrimPrDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 1)
          registry.fill(HIST("hSecWeakPrDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 2)
          registry.fill(HIST("hSecMatPrDCAxy"), track.dcaXY(), track.pt());
      }
      if (IsProton(track, -1) && track.pdgCode() == -PDG_t::kProton) {
        registry.fill(HIST("hAntiPrDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 0)
          registry.fill(HIST("hPrimAntiPrDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 1)
          registry.fill(HIST("hSecWeakAntiPrDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 2)
          registry.fill(HIST("hSecMatAntiPrDCAxy"), track.dcaXY(), track.pt());
      }
      if (IsDeuteron(track, +1) && track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("hDeDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 0)
          registry.fill(HIST("hPrimDeDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 1)
          registry.fill(HIST("hSecWeakDeDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 2)
          registry.fill(HIST("hSecMatDeDCAxy"), track.dcaXY(), track.pt());
      }
      if (IsDeuteron(track, -1) && track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("hAntiDeDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 0)
          registry.fill(HIST("hPrimAntiDeDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 1)
          registry.fill(HIST("hSecWeakAntiDeDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == 2)
          registry.fill(HIST("hSecMatAntiDeDCAxy"), track.dcaXY(), track.pt());
      }

      if (!applyDCAcut(track))
        continue;

      // Keep only protons and deuterons
      // if (std::abs(track.pdgCode()) != PDG_t::kProton && std::abs(track.pdgCode()) != o2::constants::physics::Pdg::kDeuteron)
      // continue;

      if (doQA) {
        QA.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        QA.fill(HIST("QA/hTPCSharedClusters"), track.tpcFractionSharedCls());
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

      bool isPr = (IsProton(track, +1) && track.pdgCode() == PDG_t::kProton);
      bool isAntiPr = (IsProton(track, -1) && track.pdgCode() == -PDG_t::kProton);
      bool isDe = (IsDeuteron(track, +1) && track.pdgCode() == o2::constants::physics::Pdg::kDeuteron);
      bool isAntiDe = (IsDeuteron(track, -1) && track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron);

      if (isPr) {
        registry.fill(HIST("hPrimSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * +1);
        if (track.origin() == 1 || track.origin() == 2) { // secondaries
          registry.fill(HIST("hSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * +1);
        }
      }
      if (isAntiPr) {
        registry.fill(HIST("hPrimSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        if (track.origin() == 1 || track.origin() == 2) {
          registry.fill(HIST("hSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        }
      }

      if (track.origin() != 0)
        continue;

      if (track.pdgCode() == PDG_t::kProton) {
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
      if (track.pdgCode() == -PDG_t::kProton) {
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
      if (track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
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
      if (track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
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
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.pdgCode() == PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPC"), track.pt());
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF &&
            track.pdgCode() == PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCTOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPCTOF"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.pdgCode() == PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPCEl_or_TOF"), track.pt());
        }

        // AntiProton
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.pdgCode() == -PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPC"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF &&
            track.pdgCode() == -PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCTOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPCTOF"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == -PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPC_or_TOF"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.pdgCode() == -PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == -PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPCEl_or_TOF"), track.pt() * -1);
        }

        // Deuteron
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPC"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF &&
            track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCTOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPCTOF"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl_or_TOF"), track.pt());
        }

        // AntiDeuteron
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPC"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF &&
            track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCTOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPCTOF"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPC_or_TOF"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
        }

        // Denominators
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPC"), track.pt());
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF && track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPCTOF"), track.pt());
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPC_or_TOF"), track.pt());
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl_or_TOF"), track.pt());
        }

        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPC"), track.pt() * -1);
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF && track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPCTOF"), track.pt() * -1);
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Proton_TPC_or_TOF"), track.pt() * -1);
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl_or_TOF"), track.pt() * -1);
        }

        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC"), track.pt());
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF && track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCTOF"), track.pt());
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() > 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl_or_TOF"), track.pt());

        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC"), track.pt() * -1);
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF && track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCTOF"), track.pt() * -1);
        if ((
              (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < betahasTOFthr) ||
              (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC_or_TOF"), track.pt() * -1);
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < betahasTOFthr) ||
             (track.beta() > betahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() < 0)
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
      }
    } // track
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processMC, "processMC", false);

  void processSameEventGen(SimCollisions::iterator const& mcCollision, SimParticles const& mcParticles)
  {

    if (std::abs(mcCollision.posZ()) > cutzVertex)
      return;

    registry.fill(HIST("Generated/hNEventsMC"), 0.5);

    for (const auto& particle : mcParticles) {

      if (particle.pdgCode() == PDG_t::kProton) {
        registry.fill(HIST("Generated/hQAProtons"), 0.5);
      }
      if (particle.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("Generated/hQADeuterons"), 0.5);
      }

      if (isPrim && !particle.isPhysicalPrimary()) {
        continue;
      }
      if (particle.pdgCode() == PDG_t::kProton) {
        registry.fill(HIST("Generated/hQAProtons"), 1.5);
      }
      if (particle.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("Generated/hQADeuterons"), 1.5);
      }

      if (particle.pdgCode() == o2::constants::physics::Pdg::kDeuteron && std::abs(particle.y()) < 0.5) {
        registry.fill(HIST("Generated/hDeuteronsVsPt"), particle.pt());
      }
      if (particle.pdgCode() == -o2::constants::physics::Pdg::kDeuteron && std::abs(particle.y()) < 0.5) {
        registry.fill(HIST("Generated/hAntiDeuteronsVsPt"), particle.pt());
      }

      if (std::abs(particle.eta()) > etacut) {
        continue;
      }
      if (particle.pdgCode() == PDG_t::kProton) {
        registry.fill(HIST("Generated/hQAProtons"), 2.5);
      }
      if (particle.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("Generated/hQADeuterons"), 2.5);
      }

      if (particle.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("hGen_EtaPhiPt_Deuteron"), particle.eta(), particle.phi(), particle.pt());
      }
      if (particle.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("hGen_EtaPhiPt_Deuteron"), particle.eta(), particle.phi(), -1. * particle.pt());
      }
      if (particle.pdgCode() == PDG_t::kProton) {
        registry.fill(HIST("hGen_EtaPhiPt_Proton"), particle.eta(), particle.phi(), particle.pt());
      }
      if (particle.pdgCode() == -PDG_t::kProton) {
        registry.fill(HIST("hGen_EtaPhiPt_Proton"), particle.eta(), particle.phi(), -1. * particle.pt());
      }
    }

    if (mode == kPbarPbar || mode == kPP) { // Identical particle combinations

      for (const auto& [part0, part1] : combinations(CombinationsStrictlyUpperIndexPolicy(mcParticles, mcParticles))) {

        // mode 6
        if (mode == kPP) {
          if (part0.pdgCode() != PDG_t::kProton)
            continue;
          if (part1.pdgCode() != PDG_t::kProton)
            continue;
        }
        // mode 5
        if (mode == kPbarPbar) {
          if (part0.pdgCode() != -PDG_t::kProton)
            continue;
          if (part1.pdgCode() != -PDG_t::kProton)
            continue;
        }

        fillHistogramsGen(part0, part1, false);
      }

    } else {

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(mcParticles, mcParticles))) {

        if (mode == kDbarPbar) {
          if (part0.pdgCode() != -o2::constants::physics::Pdg::kDeuteron)
            continue;
          if (part1.pdgCode() != -PDG_t::kProton)
            continue;
        }
        if (mode == kDP) {
          if (part0.pdgCode() != o2::constants::physics::Pdg::kDeuteron)
            continue;
          if (part1.pdgCode() != PDG_t::kProton)
            continue;
        }
        if (mode == kDbarP) {
          if (part0.pdgCode() != -o2::constants::physics::Pdg::kDeuteron)
            continue;
          if (part1.pdgCode() != PDG_t::kProton)
            continue;
        }
        if (mode == kDPbar) {
          if (part0.pdgCode() != o2::constants::physics::Pdg::kDeuteron)
            continue;
          if (part1.pdgCode() != -PDG_t::kProton)
            continue;
        }
        if (mode == kPbarP) {
          if (part0.pdgCode() != -PDG_t::kProton)
            continue;
          if (part1.pdgCode() != PDG_t::kProton)
            continue;
        }
        if (mode == kPPbar) {
          if (part0.pdgCode() != PDG_t::kProton)
            continue;
          if (part1.pdgCode() != -PDG_t::kProton)
            continue;
        }

        fillHistogramsGen(part0, part1, false);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processSameEventGen, "processSameEventGen", false);

  void processMixedEventGen(SimCollisions const& mcCollisions, SimParticles const& mcParticles)
  {

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningGen, 5, -1, mcCollisions, mcCollisions)) {

      //LOGF(info, "Mixed event collisions: (%d, %d) zvtx (%.1f, %.1f) mult (%d, %d)", collision1.globalIndex(), collision2.globalIndex(), collision1.posZ(), collision2.posZ(), collision1.multMCNParticlesEta05(), collision2.multMCNParticlesEta05());

      auto groupPartsOne = mcParticles.sliceByCached(o2::aod::mcparticle::mcCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = mcParticles.sliceByCached(o2::aod::mcparticle::mcCollisionId, collision2.globalIndex(), cache);

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        if (mode == kDbarPbar) {
          if (part0.pdgCode() != -o2::constants::physics::Pdg::kDeuteron)
            continue;
          if (part1.pdgCode() != -PDG_t::kProton)
            continue;
        }
        if (mode == kDP) {
          if (part0.pdgCode() != o2::constants::physics::Pdg::kDeuteron)
            continue;
          if (part1.pdgCode() != PDG_t::kProton)
            continue;
        }
        if (mode == kDbarP) {
          if (part0.pdgCode() != -o2::constants::physics::Pdg::kDeuteron)
            continue;
          if (part1.pdgCode() != PDG_t::kProton)
            continue;
        }
        if (mode == kDPbar) {
          if (part0.pdgCode() != o2::constants::physics::Pdg::kDeuteron)
            continue;
          if (part1.pdgCode() != -PDG_t::kProton)
            continue;
        }
        if (mode == kPbarP) {
          if (part0.pdgCode() != -PDG_t::kProton)
            continue;
          if (part1.pdgCode() != PDG_t::kProton)
            continue;
        }
        if (mode == kPbarPbar) {
          if (part0.pdgCode() != -PDG_t::kProton)
            continue;
          if (part1.pdgCode() != -PDG_t::kProton)
            continue;
        }
        if (mode == kPP) {
          if (part0.pdgCode() != PDG_t::kProton)
            continue;
          if (part1.pdgCode() != PDG_t::kProton)
            continue;
        }
        if (mode == kPPbar) {
          if (part0.pdgCode() != PDG_t::kProton)
            continue;
          if (part1.pdgCode() != -PDG_t::kProton)
            continue;
        }

        fillHistogramsGen(part0, part1, true);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processMixedEventGen, "processMixedEventGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HadronNucleiCorrelation>(cfgc)};
}
