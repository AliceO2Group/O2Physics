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
/// \file hadronnucleicorrelation.cxx
/// \brief Hadron-nuclei correlation analysis task
/// \author Francesca Ercolessi
/// \since 21 April 2024

#include "PWGCF/Femto3D/Core/femto3dPairTask.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Common/Core/RecoDecay.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TString.h>

#include <chrono>
#include <cmath>
#include <cstdint>
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

enum Origin {
  kPrimary = 0,
  kWeakDecay,
  kMaterial
};

struct HadronNucleiCorrelation {

  static constexpr int BetahasTOFthr = -100;

  SliceCache cache;

  Configurable<int> mode{"mode", 0, "0: antid-antip, 1: d-p, 2: antid-p, 3: d-antip, 4: antip-p, 5: antip-antip, 6: p-p, 7: p-antip"};

  Configurable<bool> doRapidity{"doRapidity", false, "do rapidity dependent analysis"};
  Configurable<bool> doQA{"doQA", true, "save QA histograms"};
  Configurable<bool> doMCQA{"doMCQA", false, "save MC QA histograms"};
  Configurable<bool> isMC{"isMC", false, "is MC"};
  Configurable<bool> isMCGen{"isMCGen", false, "is isMCGen"};
  Configurable<bool> isPrim{"isPrim", true, "is isPrim"};
  Configurable<bool> doCorrection{"doCorrection", false, "do efficiency correction"};
  Configurable<bool> doQuadraticPID{"doQuadraticPID", false, "do PID with sum in quadrature of TOF and TPC"};

  Configurable<std::string> fCorrectionPath{"fCorrectionPath", "", "Correction path to file"};
  Configurable<std::string> fCorrectionHisto{"fCorrectionHisto", "", "Correction histogram"};
  Configurable<std::string> cfgUrl{"cfgUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Event selection
  Configurable<float> cutzVertex{"cutzVertex", 10.0, "|vertexZ| value limit"};
  Configurable<bool> removeSameBunchPileup{"removeSameBunchPileup", false, "remove Same Bunch Pileup"};

  // Track selection
  Configurable<bool> doClosePairRejection{"doClosePairRejection", false, "doClosePairRejection"};
  Configurable<double> dcaPar0{"dcaPar0", 0.004, "par 0"};
  Configurable<double> dcaPar1{"dcaPar1", 0.013, "par 1"};
  Configurable<bool> doDCAZ{"doDCAZ", true, "do DCA z cut"};
  Configurable<int16_t> minTPCnClusters{"minTPCnClusters", 80, "minimum number of found TPC clusters"};
  Configurable<float> minTPCnCrossedRowsOverFindableCls{"minTPCnCrossedRowsOverFindableCls", 0.8, "n TPC Crossed Rows Over Findable Cls"};
  Configurable<float> maxchi2TPC{"maxchi2TPC", 4.0f, "maximum TPC chi^2/Ncls"};
  Configurable<float> maxchi2ITS{"maxchi2ITS", 36.0f, "maximum ITS chi^2/Ncls"};
  Configurable<float> etaCut{"etaCut", 0.8f, "eta cut"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.14f, "Maximum DCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.1f, "Maximum DCAz"};
  Configurable<float> nsigmaTPC{"nsigmaTPC", 3.0f, "cut nsigma TPC"};
  Configurable<float> nsigmaElPr{"nsigmaElPr", 1.0f, "cut nsigma TPC El for protons"};
  Configurable<float> nsigmaElDe{"nsigmaElDe", 3.0f, "cut nsigma TPC El for protons"};
  Configurable<float> nsigmaTOF{"nsigmaTOF", 3.5f, "cut nsigma TOF"};
  Configurable<float> nsigmaITSPr{"nsigmaITSPr", -2.0f, "cut nsigma ITS Pr"};
  Configurable<float> nsigmaITSDe{"nsigmaITSDe", -2.0f, "cut nsigma ITS De"};
  Configurable<bool> doITSPID{"doITSPID", true, "do ITS PID"};
  Configurable<float> pTthrprTOF{"pTthrprTOF", 0.8f, "threshold pT proton to use TOF"};
  Configurable<float> pTthrprTPCEl{"pTthrprTPCEl", 1.0f, "threshold pT proton to use TPC El rejection"};
  Configurable<float> pTthrdeTOF{"pTthrdeTOF", 1.0f, "threshold pT deuteron to use TOF"};
  Configurable<float> pTthrdeTPCEl{"pTthrdeTPCEl", 1.0f, "threshold pT deuteron to use TPC El rejection"};
  Configurable<bool> rejectionEl{"rejectionEl", true, "use TPC El rejection"};
  Configurable<float> maxtpcSharedCls{"maxtpcSharedCls", 0.4, "maximum fraction of TPC shared clasters"};
  Configurable<int> minitsNCls{"minitsNCls", 0, "minimum allowed number of ITS clasters"};
  Configurable<int> maxmixcollsGen{"maxmixcollsGen", 100, "maxmixcollsGen"};
  Configurable<float> radiusTPC{"radiusTPC", 1.2, "TPC radius to calculate phi_star for"};
  Configurable<float> dEta{"dEta", 0.01, "minimum allowed difference in eta between two tracks in a pair"};
  Configurable<float> dPhi{"dPhi", 0.01, "minimum allowed difference in phi_star between two tracks in a pair"};
  Configurable<float> yRap{"yRap", 0.5, "rapidity"};

  // Mixing parameters
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 50.0f, 100.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ColumnBinningPolicy<aod::singletrackselector::PosZ, aod::singletrackselector::Mult> colBinning{{confVtxBins, confMultBins}, true};

  // pT/A bins
  Configurable<std::vector<double>> pTBins{"pTBins", {0.6f, 1.0f, 1.2f, 2.f}, "p_{T} bins"};

  ConfigurableAxis axisNSigma{"axisNSigma", {35, -7.f, 7.f}, "n#sigma"};
  ConfigurableAxis deltaPhiAxis = {"deltaPhiAxis", {46, -1 * o2::constants::math::PIHalf, 3 * o2::constants::math::PIHalf}, "#Delta#phi (rad)"};

  using FilteredCollisions = soa::Filtered<aod::SingleCollSels>;
  using FilteredCollisionsExtra = soa::Filtered<soa::Join<aod::SingleCollSels, aod::SingleCollExtras>>;
  using SimCollisions = soa::Filtered<aod::McCollisions>;
  using SimParticles = aod::McParticles;
  using FilteredTracks = soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkExtras, aod::SinglePIDEls, aod::SinglePIDPrs, aod::SinglePIDDes>>;                      // new tables (v3)
  using FilteredTracksMC = soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkMCs, aod::SingleTrkExtras, aod::SinglePIDEls, aod::SinglePIDPrs, aod::SinglePIDDes>>; // new tables (v3)

  HistogramRegistry registry{"registry"};
  HistogramRegistry registryQa{"registryQa"};

  using TrkType = const FilteredTracks::iterator*;
  // using TrkTypeMC = const FilteredTracksMC::iterator*;
  // typedef std::shared_ptr<FilteredCollisions::iterator> ColType;
  // typedef std::shared_ptr<SimCollisions::iterator> MCcolType;

  std::unique_ptr<o2::aod::singletrackselector::FemtoPair<TrkType>> pair = std::make_unique<o2::aod::singletrackselector::FemtoPair<TrkType>>();
  // std::unique_ptr<o2::aod::singletrackselector::FemtoPair<TrkTypeMC>> PairMC = std::make_unique<o2::aod::singletrackselector::FemtoPair<TrkTypeMC>>();

  // Data histograms
  std::vector<std::shared_ptr<TH3>> hEtaPhiSameEv;
  std::vector<std::shared_ptr<TH3>> hEtaPhiMixdEv;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhiSameEv;
  std::vector<std::shared_ptr<TH3>> hCorrEtaPhiMixdEv;

  int nBinspT = 0;
  TH2F* hEffPtEtaProton = nullptr;
  TH2F* hEffPtEtaAntiProton = nullptr;
  TH2F* hEffPtEtaDeuteron = nullptr;
  TH2F* hEffPtEtaAntiDeuteron = nullptr;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::ccdb::CcdbApi ccdbApi;

  Service<o2::framework::O2DatabasePDG> pdgDB{};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(cfgUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    if (doCorrection) {
      getCorrection(ccdb, TString(fCorrectionPath), TString(fCorrectionHisto));
    }

    const AxisSpec ptBinnedAxis = {pTBins, "#it{p}_{T} of #bar{p} (GeV/#it{c})"};
    const AxisSpec etaAxis = {100, -1., 1., "#eta"};
    const AxisSpec phiAxis = {157, 0., o2::constants::math::TwoPI, "#phi (rad)"};
    const AxisSpec ptAxis = {200, -10.f, 10.f, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec ptAxisSmall = {100, -5.f, 5.f, "#it{p}_{T} GeV/#it{c}"};

    const AxisSpec deltaEtaAxis = {300, -1.5, 1.5, "#Delta#eta"};
    const AxisSpec deltaRapAxis = {300, -1.5, 1.5, "#Delta y"};

    if (doprocessSameEvent || doprocessSameEventEvSel) {
      registry.add("hNEvents", "hNEvents", {HistType::kTH1D, {{7, 0.f, 7.f}}});
      registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "Selected");
      registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "Mixing");
    }

    // Not used, commented out for now
    // registry.add("hNtrig_total", "hNtrig_total", {HistType::kTH1D, {ptBinnedAxis}});

    nBinspT = pTBins.value.size() - 1;

    TString name = "Undefined";
    switch (mode) {
      case kDbarPbar: // 0
        name = "AntiDeAntiPr";
        break;
      case kDP: // 1
        name = "DePr";
        break;
      case kDbarP: // 2
        name = "AntiDePr";
        break;
      case kDPbar: // 3
        name = "DeAntiPr";
        break;
      case kPbarP: // 4
        name = "AntiPrPr";
        break;
      case kPbarPbar: // 5
        name = "AntiPrAntiPr";
        break;
      case kPP: // 6
        name = "PrPr";
        break;
      case kPPbar: // 7
        name = "PrAntiPr";
        break;
      default:
        LOG(fatal) << "Unhandled case " << mode;
    }

    if (!isMC) {
      for (int i = 0; i < nBinspT; i++) {
        const TString ptTag = Form("pt%02.0f%02.0f", pTBins.value.at(i) * 10, pTBins.value.at(i + 1) * 10);
        const TString ptInterval = Form("(%.1f<p_{T}^{assoc} <%.1f GeV/c)", pTBins.value.at(i), pTBins.value.at(i + 1));
        if (doRapidity) {
          hEtaPhiSameEv.push_back(registry.add<TH3>(Form("hEtaPhi_%s_SE_%s", name.Data(), ptTag.Data()), "Raw #Delta y #Delta#phi " + ptInterval, {HistType::kTH3F, {deltaRapAxis, deltaPhiAxis, ptBinnedAxis}}));
          hEtaPhiMixdEv.push_back(registry.add<TH3>(Form("hEtaPhi_%s_ME_%s", name.Data(), ptTag.Data()), "Raw #Delta y #Delta#phi " + ptInterval, {HistType::kTH3F, {deltaRapAxis, deltaPhiAxis, ptBinnedAxis}}));

          hCorrEtaPhiSameEv.push_back(registry.add<TH3>(Form("hCorrEtaPhi_%s_SE_%s", name.Data(), ptTag.Data()), "#Delta y #Delta#phi " + ptInterval, {HistType::kTH3F, {deltaRapAxis, deltaPhiAxis, ptBinnedAxis}}));
          hCorrEtaPhiMixdEv.push_back(registry.add<TH3>(Form("hCorrEtaPhi_%s_ME_%s", name.Data(), ptTag.Data()), "#Delta y #Delta#phi " + ptInterval, {HistType::kTH3F, {deltaRapAxis, deltaPhiAxis, ptBinnedAxis}}));
        } else {
          hEtaPhiSameEv.push_back(registry.add<TH3>(Form("hEtaPhi_%s_SE_%s", name.Data(), ptTag.Data()), "Raw #Delta#eta#Delta#phi " + ptInterval, {HistType::kTH3F, {deltaEtaAxis, deltaPhiAxis, ptBinnedAxis}}));
          hEtaPhiMixdEv.push_back(registry.add<TH3>(Form("hEtaPhi_%s_ME_%s", name.Data(), ptTag.Data()), "Raw #Delta#eta#Delta#phi " + ptInterval, {HistType::kTH3F, {deltaEtaAxis, deltaPhiAxis, ptBinnedAxis}}));

          hCorrEtaPhiSameEv.push_back(registry.add<TH3>(Form("hCorrEtaPhi_%s_SE_%s", name.Data(), ptTag.Data()), "#Delta#eta#Delta#phi " + ptInterval, {HistType::kTH3F, {deltaEtaAxis, deltaPhiAxis, ptBinnedAxis}}));
          hCorrEtaPhiMixdEv.push_back(registry.add<TH3>(Form("hCorrEtaPhi_%s_ME_%s", name.Data(), ptTag.Data()), "#Delta#eta#Delta#phi " + ptInterval, {HistType::kTH3F, {deltaEtaAxis, deltaPhiAxis, ptBinnedAxis}}));
        }
      }
    }

    if (doprocessSameEvent || doprocessSameEventEvSel || doprocessMC) {
      registry.add("hPrDCAxy", "DCAxy p", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hAntiPrDCAxy", "DCAxy #bar{p}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hDeDCAxy", "DCAxy d", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registry.add("hAntiDeDCAxy", "DCAxy #bar{d}", {HistType::kTH2D, {{600, -3.f, 3.f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
    }
    registry.add("hMult", "multiplicity", {HistType::kTH1D, {{200, 0.f, 200.f, "N_{ch}"}}});

    if (doQA) {
      // Track QA
      registryQa.add("QA/hVtxZ_trk", "#it{z}_{vtx}", {HistType::kTH1D, {{150, -15.f, 15.f, "#it{z}_{vtx} (cm)"}}});
      registryQa.add("QA/hTPCnClusters", "N TPC Clusters; N TPC Clusters", {HistType::kTH1D, {{200, 0.f, 200.f}}});
      registryQa.add("QA/hTPCSharedClusters", "N TPC Shared Clusters; N TPC SharedClusters", {HistType::kTH1D, {{100, 0.f, 1.f}}});
      registryQa.add("QA/hTPCchi2", "TPC chi2/Ncls; TPC chi2/Ncls", {HistType::kTH1D, {{100, 0.f, 10.f}}});
      registryQa.add("QA/hTPCcrossedRowsOverFindableCls", "TPC crossed Rows Over Findable Cls; TPC Crossed Rows Over Findable Cls", {HistType::kTH1D, {{100, 0.f, 2.f}}});
      registryQa.add("QA/hITSchi2", "ITS chi2/Ncls; ITS chi2/Ncls", {HistType::kTH1D, {{100, 0.f, 20.f}}});
      registryQa.add("QA/hDCAxy", "DCAxy", {HistType::kTH2D, {{200, -0.2f, 0.2f, "DCA xy (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registryQa.add("QA/hDCAz", "DCAz", {HistType::kTH2D, {{200, -0.2f, 0.2f, "DCA z (cm)"}, {100, 0.f, 10.f, "p_{T} GeV/c"}}});
      registryQa.add("QA/TPCChi2VsPZ", "TPCChi2VsPZ", {HistType::kTH2D, {{100, 0.f, 10.f, "p_{TPC}/Z (GeV/c)"}, {120, 0.f, 6.f, "TPC Chi2"}}});
      const AxisSpec tofNSigmaAxis = {axisNSigma, "n#sigma TOF"};
      const AxisSpec tpcNSigmaAxis = {axisNSigma, "n#sigma TPC"};
      const AxisSpec itsNSigmaAxis = {axisNSigma, "n#sigma ITS"};
      registryQa.add("QA/h2dTPCTOF_Pr", "n#sigma TPC vs n#sigma TOF", {HistType::kTH2D, {tpcNSigmaAxis, tofNSigmaAxis}});
      registryQa.add("QA/h2dTPCTOF_AntiPr", "n#sigma TPC vs n#sigma TOF", {HistType::kTH2D, {tpcNSigmaAxis, tofNSigmaAxis}});
      registryQa.add("QA/hnSigmaTPCVsPt_El", "n#sigma TPC vs p_{T} for e hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tpcNSigmaAxis}});
      registryQa.add("QA/hnSigmaTPCVsPt_Pr", "n#sigma TPC vs p_{T} for p hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tpcNSigmaAxis}});
      registryQa.add("QA/hnSigmaTPCVsPt_De", "n#sigma TPC vs p_{T} for d hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tpcNSigmaAxis}});
      registryQa.add("QA/hnSigmaTOFVsPt_Pr", "n#sigma TOF vs p_{T} for p hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tofNSigmaAxis}});
      registryQa.add("QA/hnSigmaTOFVsPt_De", "n#sigma TOF vs p_{T} for d hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tofNSigmaAxis}});
      registryQa.add("QA/hnSigmaITSVsPt_Pr", "n#sigma ITS vs p_{T} for p hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, itsNSigmaAxis}});
      registryQa.add("QA/hnSigmaITSVsPt_De", "n#sigma ITS vs p_{T} for d hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, itsNSigmaAxis}});
      registryQa.add("QA/hdEtadPhistar", ";dPhi*;dEta ", {HistType::kTH2D, {{101, -0.2, 0.2, "dPhi*"}, {101, -0.2, 0.2, "dEta"}}});

      if (!isMC) {
        registryQa.add("QA/hEtaPr", Form("#eta ditribution for p"), {HistType::kTH1F, {etaAxis}});
        registryQa.add("QA/hPhiPr", Form("#phi ditribution for p"), {HistType::kTH1F, {phiAxis}});
        registryQa.add("QA/hEtaAntiPr", Form("#eta ditribution for #bar{p}"), {HistType::kTH1F, {etaAxis}});
        registryQa.add("QA/hPhiAntiPr", Form("#phi ditribution for #bar{p}"), {HistType::kTH1F, {phiAxis}});
        registryQa.add("QA/hEtaDe", Form("#eta ditribution for d"), {HistType::kTH1F, {etaAxis}});
        registryQa.add("QA/hPhiDe", Form("#phi ditribution for d"), {HistType::kTH1F, {phiAxis}});
        registryQa.add("QA/hEtaAntiDe", Form("#eta ditribution for #bar{d}"), {HistType::kTH1F, {etaAxis}});
        registryQa.add("QA/hPhiAntiDe", Form("#phi ditribution for #bar{d}"), {HistType::kTH1F, {phiAxis}});

        registryQa.add("QA/hnSigmaTPCVsPt_Pr_AfterSel", "n#sigma TPC vs p_{T} for p hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tpcNSigmaAxis}});
        registryQa.add("QA/hnSigmaTPCVsPt_De_AfterSel", "n#sigma TPC vs p_{T} for d hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tpcNSigmaAxis}});
        registryQa.add("QA/hnSigmaTOFVsPt_Pr_AfterSel", "n#sigma TOF vs p_{T} for p hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tofNSigmaAxis}});
        registryQa.add("QA/hnSigmaTOFVsPt_De_AfterSel", "n#sigma TOF vs p_{T} for d hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, tofNSigmaAxis}});
        registryQa.add("QA/hnSigmaITSVsPt_Pr_AfterSel", "n#sigma ITS vs p_{T} for p hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, itsNSigmaAxis}});
        registryQa.add("QA/hnSigmaITSVsPt_De_AfterSel", "n#sigma ITS vs p_{T} for d hypothesis (all tracks)", {HistType::kTH2D, {ptAxis, itsNSigmaAxis}});
        registryQa.add("QA/h2dTPCTOF_Pr_AfterSel", "n#sigma TPC vs n#sigma TOF", {HistType::kTH2D, {tpcNSigmaAxis, tofNSigmaAxis}});
        registryQa.add("QA/h2dTPCTOF_AntiPr_AfterSel", "n#sigma TPC vs n#sigma TOF", {HistType::kTH2D, {tpcNSigmaAxis, tofNSigmaAxis}});
      }
    }

    if (isMC && doprocessMC) {
      const AxisSpec dcaAxis = {600, -3.f, 3.f, "DCA xy (cm)"};
      const AxisSpec mcPtAxis = {100, 0.f, 10.f, "#it{p}_{T} GeV/#it{c}"};
      registry.add("hPrimPrDCAxy", "DCAxy p", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hPrimAntiPrDCAxy", "DCAxy #bar{p}", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hPrimDeDCAxy", "DCAxy d", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hPrimAntiDeDCAxy", "DCAxy #bar{d}", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hSecMatPrDCAxy", "DCAxy p", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hSecMatAntiPrDCAxy", "DCAxy #bar{p}", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hSecMatDeDCAxy", "DCAxy d", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hSecMatAntiDeDCAxy", "DCAxy #bar{d}", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hSecWeakPrDCAxy", "DCAxy p", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hSecWeakAntiPrDCAxy", "DCAxy #bar{p}", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hSecWeakDeDCAxy", "DCAxy d", {HistType::kTH2D, {dcaAxis, mcPtAxis}});
      registry.add("hSecWeakAntiDeDCAxy", "DCAxy #bar{d}", {HistType::kTH2D, {dcaAxis, mcPtAxis}});

      registry.add("hReco_EtaPhiPt_Proton", "Gen (anti)protons in reco collisions", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});
      registry.add("hReco_EtaPhiPt_Deuteron", "Gen (anti)deuteron in reco collisions", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});
      registry.add("hReco_PID_EtaPhiPt_Proton", "Gen (anti)protons + PID in reco collisions", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});
      registry.add("hReco_PID_EtaPhiPt_Deuteron", "Gen (anti)deuteron + PID in reco collisions", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});
      registry.add("hReco_EtaPhiPtMC_Proton", "Gen (anti)protons in reco collisions (MC info used)", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});
      registry.add("hReco_EtaPhiPtMC_Deuteron", "Gen (anti)deuteron in reco collisions (MC info used)", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});
      registry.add("hReco_Pt_Proton", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
      registry.add("hReco_Pt_Deuteron", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});

      registry.add("hSec_EtaPhiPt_Proton", "Secondary (anti)protons", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});
      registry.add("hPrimSec_EtaPhiPt_Proton", "Primary + Secondary (anti)protons", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});

      registry.add("hnSigmaTPCVsPt_Pr_MC", "n#sigma TPC vs p_{T} for p hypothesis true MC; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {ptAxis, axisNSigma}});
      registry.add("hnSigmaTPCVsPt_De_MC", "n#sigma TPC vs p_{T} for d hypothesis true MC; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2F, {ptAxis, axisNSigma}});
      registry.add("hnSigmaTOFVsPt_Pr_MC", "n#sigma TOF vs p_{T} for p hypothesis true MC; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {ptAxis, axisNSigma}});
      registry.add("hnSigmaTOFVsPt_De_MC", "n#sigma TOF vs p_{T} for d hypothesis true MC; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2F, {ptAxis, axisNSigma}});

      registry.add("hResPt_Proton", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResPt_Deuteron", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResPt_AntiProton", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});
      registry.add("hResPt_AntiDeuteron", "; p_{T}(gen) [GeV/c]; p_{T}(reco) - p_{T}(gen) ", {HistType::kTH2F, {{100, 0.f, 10.f, "p_{T}(gen) GeV/c"}, {200, -1.f, 1.f, "p_{T}(reco) - p_{T}(gen) "}}});

      registry.add("hNumeratorPurity_Proton", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
      registry.add("hNumeratorPurity_Deuteron", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
      registry.add("hDenominatorPurity_Proton", " p(#bar{p}); p_{T} (GeV/c);(S + B)", {HistType::kTH1F, {ptAxisSmall}});
      registry.add("hDenominatorPurity_Deuteron", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});

      if (doMCQA) {

        registry.add("hResEta_Proton", "; #eta(gen); #eta(reco) - #eta(gen)  ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResEta_Deuteron", "; #eta(gen); #eta(reco) - #eta(gen) ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResPhi_Proton", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, o2::constants::math::TwoPI, "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResPhi_Deuteron", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, o2::constants::math::TwoPI, "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResEta_AntiProton", "; #eta(gen); #eta(reco) - #eta(gen)  ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResEta_AntiDeuteron", "; #eta(gen); #eta(reco) - #eta(gen) ", {HistType::kTH2F, {{100, -1.f, 1.f, "#eta(gen)"}, {200, -0.5f, 0.5f, "#eta(reco) - #eta(gen) "}}});
        registry.add("hResPhi_AntiProton", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, o2::constants::math::TwoPI, "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});
        registry.add("hResPhi_AntiDeuteron", "; #phi(gen); #phi(reco) - #phi(gen)", {HistType::kTH2F, {{100, 0.f, o2::constants::math::TwoPI, "#phi(gen)"}, {200, -0.5f, 0.5f, "#phi(reco) - #phi(gen)"}}});

        registry.add("hNumeratorPurity_Proton_TPC", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Deuteron_TPC", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Proton_TPCTOF", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Deuteron_TPCTOF", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Proton_TPC_or_TOF", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Deuteron_TPC_or_TOF", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Proton_TPCEl_or_TOF", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Proton_TPCEl", " p(#bar{p}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Deuteron_TPCEl", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hNumeratorPurity_Deuteron_TPCEl_or_TOF", " d(#bar{d}); p_{T} (GeV/c);S", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Proton_TPC", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Deuteron_TPC", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Proton_TPCTOF", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Deuteron_TPCTOF", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Proton_TPC_or_TOF", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Deuteron_TPC_or_TOF", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Proton_TPCEl", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Proton_TPCEl_or_TOF", " p(#bar{p}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Deuteron_TPCEl", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hDenominatorPurity_Deuteron_TPCEl_or_TOF", " d(#bar{d}); p_{T} (GeV/c); (S + B)", {HistType::kTH1F, {ptAxisSmall}});

        registry.add("hReco_Pt_Proton_TPC", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Deuteron_TPC", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Proton_TPCTOF", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Deuteron_TPCTOF", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Proton_TPC_or_TOF", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Deuteron_TPC_or_TOF", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Proton_TPCEl", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Proton_TPCEl_or_TOF", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Deuteron_TPCEl", "Reco (anti)deuterons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
        registry.add("hReco_Pt_Deuteron_TPCEl_or_TOF", "Reco (anti)protons in reco collisions", {HistType::kTH1F, {ptAxisSmall}});
      }
    }

    if (isMCGen) {
      registry.add("Generated/hNEventsMC", "hNEventsMC", {HistType::kTH1D, {{1, 0.f, 1.f}}});
      registry.get<TH1>(HIST("Generated/hNEventsMC"))->GetXaxis()->SetBinLabel(1, "All");

      registry.add("hGen_EtaPhiPt_Proton", "Gen (anti)protons in gen collisions", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});
      registry.add("hGen_EtaPhiPt_Deuteron", "Gen (anti)deuteron in gen collisions", {HistType::kTH3F, {etaAxis, phiAxis, ptAxisSmall}});

      auto hp = registry.add<TH1>("Generated/hQAProtons", "hQAProtons", {HistType::kTH1D, {{5, 0.f, 5.f}}});
      hp->GetXaxis()->SetBinLabel(1, "All");
      hp->GetXaxis()->SetBinLabel(2, "PhysicalPrimary");
      hp->GetXaxis()->SetBinLabel(3, "|#eta|<0.8");
      hp->GetXaxis()->SetBinLabel(4, "no daughters");
      hp->GetXaxis()->SetBinLabel(5, "d daughter");
      registry.addClone("Generated/hQAProtons", "Generated/hQAAntiProtons");

      registry.addClone("Generated/hQAProtons", "Generated/hQANeutrons");
      registry.addClone("Generated/hQAProtons", "Generated/hQAAntiNeutrons");

      auto hd = registry.add<TH1>("Generated/hQADeuterons", "hQADeuterons", {HistType::kTH1D, {{3, 0.f, 3.f}}});
      hd->GetXaxis()->SetBinLabel(1, "All");
      hd->GetXaxis()->SetBinLabel(2, "PhysicalPrimary");
      hd->GetXaxis()->SetBinLabel(3, "|#eta|<0.8");
      registry.addClone("Generated/hQADeuterons", "Generated/hQAAntiDeuterons");

      const AxisSpec ptAxisGen = {100, 0.f, 10.f, "#it{p}_{T} GeV/#it{c}"};
      registry.add("Generated/hDeuteronsVsPt", "hDeuteronsVsPt;", {HistType::kTH1D, {ptAxisGen}});
      registry.add("Generated/hAntiDeuteronsVsPt", "hAntiDeuteronsVsPt;", {HistType::kTH1D, {ptAxisGen}});
    }
  }

  // Filters
  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) <= cutzVertex;
  Filter trackFilter = o2::aod::singletrackselector::tpcNClsFound >= minTPCnClusters &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedTpcChi2NCl) <= maxchi2TPC &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::rowsOverFindable>(o2::aod::singletrackselector::storedTpcCrossedRowsOverFindableCls) >= minTPCnCrossedRowsOverFindableCls &&
                       o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedItsChi2NCl) <= maxchi2ITS &&
                       nabs(o2::aod::singletrackselector::unPack<singletrackselector::binning::dca>(o2::aod::singletrackselector::storedDcaXY)) <= maxDCAxy &&
                       nabs(o2::aod::singletrackselector::unPack<singletrackselector::binning::dca>(o2::aod::singletrackselector::storedDcaXY)) <= maxDCAz &&
                       nabs(o2::aod::singletrackselector::eta) <= etaCut;

  Filter simvertexFilter = nabs(o2::aod::mccollision::posZ) <= cutzVertex;

  template <typename Type>
  bool isProton(Type const& track, const int sign)
  {
    const bool isTPCPID = std::abs(track.tpcNSigmaPr()) < nsigmaTPC;
    const bool isTOFPID = std::abs(track.tofNSigmaPr()) < nsigmaTOF;
    const bool isTPCElRejection = rejectionEl && track.beta() < BetahasTOFthr && track.pt() < pTthrprTPCEl && track.tpcNSigmaEl() >= nsigmaElPr;
    const bool isITSPID = track.itsNSigmaPr() > nsigmaITSPr;

    const bool isQuadraticPID = std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()) < nsigmaTPC;

    // Check if the sign of the track matches the expected sign for protons or antiprotons
    const bool signCheck = (sign > 0 && track.sign() > 0) || (sign < 0 && track.sign() < 0);
    if (!doQuadraticPID) {
      if (isTPCPID) {
        if (track.pt() < pTthrprTOF) {
          if (!doITSPID || isITSPID) {
            return signCheck;
          }
        } else if (isTPCElRejection || isTOFPID) {
          return signCheck;
        }
      }
    } else {
      if (track.pt() < pTthrprTOF) {
        if (isTPCPID) {
          if (!doITSPID || isITSPID) {
            return signCheck;
          }
        }
      } else if (isQuadraticPID) {
        return signCheck;
      }
    }
    return false;
  }

  template <typename Type>
  bool isDeuteron(Type const& track, const int sign)
  {
    const bool isTPCPID = std::abs(track.tpcNSigmaDe()) < nsigmaTPC;
    const bool isTOFPID = std::abs(track.tofNSigmaDe()) < nsigmaTOF;
    const bool isTPCElRejection = rejectionEl && track.beta() < BetahasTOFthr && track.pt() < pTthrdeTPCEl && track.tpcNSigmaEl() >= nsigmaElDe;
    const bool isITSPID = track.itsNSigmaDe() > nsigmaITSDe;

    const bool isQuadraticPID = std::hypot(track.tpcNSigmaDe(), track.tofNSigmaDe()) < nsigmaTPC;

    // Check if the sign of the track matches the expected sign for deuterons or antideuterons
    const bool signCheck = (sign > 0 && track.sign() > 0) || (sign < 0 && track.sign() < 0);
    if (!doQuadraticPID) {
      if (isTPCPID) {
        if (track.pt() < pTthrdeTOF) {
          if (!doITSPID || isITSPID) {
            return signCheck;
          }
        } else if (isTPCElRejection || isTOFPID) {
          return signCheck;
        }
      }
    } else {
      if (track.pt() < pTthrprTOF) {
        if (isTPCPID) {
          if (!doITSPID || isITSPID) {
            return signCheck;
          }
        }
      } else if (isQuadraticPID) {
        return signCheck;
      }
    }
    return false;
  }

  template <typename T1>
  bool applyDCAcut(const T1& track)
  {
    // pt-dependent selection
    if (std::abs(track.dcaXY()) > (dcaPar0 + dcaPar1 / track.pt())) {
      return false;
    }
    if (doDCAZ && std::abs(track.dcaZ()) > (dcaPar0 + dcaPar1 / track.pt())) {
      return false;
    }
    return true;
  }

  template <typename T1>
  void fillHistograms(T1 const& part0, T1 const& part1, bool ME, bool isIdentical)
  {
    pair->SetPair(&part0, &part1);
    pair->SetIdentical(isIdentical);
    if (isIdentical && pair->IsClosePair(dEta, dPhi, radiusTPC)) {
      registryQa.fill(HIST("QA/hdEtadPhistar"), pair->GetPhiStarDiff(radiusTPC), pair->GetEtaDiff());
      return;
    }

    if (doClosePairRejection && pair->IsClosePair(dEta, dPhi, radiusTPC)) {
      registryQa.fill(HIST("QA/hdEtadPhistar"), pair->GetPhiStarDiff(radiusTPC), pair->GetEtaDiff());
      return;
    }

    float deltaEta = part0.eta() - part1.eta();
    float deltaPhi = part0.phi() - part1.phi();
    deltaPhi = RecoDecay::constrainAngle(deltaPhi, -1 * o2::constants::math::PIHalf);

    for (int k = 0; k < nBinspT; k++) {

      if (part0.pt() >= pTBins.value.at(k) && part0.pt() < pTBins.value.at(k + 1)) {

        float corr0 = 1, corr1 = 1;

        if (doCorrection) { // Apply corrections
          switch (mode) {
            case kDbarPbar:
              corr0 = hEffPtEtaAntiDeuteron->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffPtEtaAntiProton->Interpolate(part1.pt(), part1.eta());
              break;
            case kDP:
              corr0 = hEffPtEtaDeuteron->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffPtEtaProton->Interpolate(part1.pt(), part1.eta());
              break;
            case kDbarP:
              corr0 = hEffPtEtaAntiDeuteron->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffPtEtaProton->Interpolate(part1.pt(), part1.eta());
              break;
            case kDPbar:
              corr0 = hEffPtEtaDeuteron->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffPtEtaAntiProton->Interpolate(part1.pt(), part1.eta());
              break;
            case kPbarP:
              corr0 = hEffPtEtaAntiProton->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffPtEtaProton->Interpolate(part1.pt(), part1.eta());
              break;
            case kPbarPbar:
              corr0 = hEffPtEtaAntiProton->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffPtEtaAntiProton->Interpolate(part1.pt(), part1.eta());
              break;
            case kPP:
              corr0 = hEffPtEtaProton->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffPtEtaProton->Interpolate(part1.pt(), part1.eta());
              break;
            case kPPbar:
              corr0 = hEffPtEtaProton->Interpolate(part0.pt(), part0.eta());
              corr1 = hEffPtEtaAntiProton->Interpolate(part1.pt(), part1.eta());
              break;
            default:
              LOG(error) << "Unknown mode for efficiency correction: " << mode;
              break;
          }
        }

        if (ME) {
          hEtaPhiMixdEv[k]->Fill(deltaEta, deltaPhi, part1.pt());
          if (corr0 != 0 && corr1 != 0) {
            hCorrEtaPhiMixdEv[k]->Fill(deltaEta, deltaPhi, part1.pt(), 1. / (corr0 * corr1));
          }
        } else {
          hEtaPhiSameEv[k]->Fill(deltaEta, deltaPhi, part1.pt());
          if (corr0 != 0 && corr1 != 0) {
            hCorrEtaPhiSameEv[k]->Fill(deltaEta, deltaPhi, part1.pt(), 1. / (corr0 * corr1));
          }
        } // SE
      } // pT condition
    } // nBinspT loop

    pair->ResetPair();
  }

  template <typename T1>
  void fillHistogramsGen(T1 const& part0, T1 const& part1, const bool ME)
  {

    float deltaEta = part0.eta() - part1.eta();
    float deltaPhi = part0.phi() - part1.phi();
    deltaPhi = RecoDecay::constrainAngle(deltaPhi, -1 * o2::constants::math::PIHalf);

    for (int k = 0; k < nBinspT; k++) {

      if (part0.pt() >= pTBins.value.at(k) && part0.pt() < pTBins.value.at(k + 1)) {

        if (ME) {
          hEtaPhiMixdEv[k]->Fill(deltaEta, deltaPhi, part1.pt());
          hCorrEtaPhiMixdEv[k]->Fill(deltaEta, deltaPhi, part1.pt());
        } else {
          hEtaPhiSameEv[k]->Fill(deltaEta, deltaPhi, part1.pt());
          hCorrEtaPhiSameEv[k]->Fill(deltaEta, deltaPhi, part1.pt());
        } // SE
      } // pT condition
    } // nBinspT loop
  }

  void getCorrection(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj, const TString& filepath, const TString& histname)
  {
    auto* l = ccdbObj->get<TList>(filepath.Data());
    if (!l) {
      LOGP(error, "Could not open corrections file {}", Form("%s", filepath.Data()));
      return;
    }
    hEffPtEtaProton = dynamic_cast<TH2F*>(l->FindObject(Form("%s_proton", histname.Data())));
    if (!hEffPtEtaProton) {
      LOGP(error, "Could not open histogram {}", Form("%s_proton", histname.Data()));
      return;
    }
    hEffPtEtaAntiProton = dynamic_cast<TH2F*>(l->FindObject(Form("%s_antiproton", histname.Data())));
    if (!hEffPtEtaAntiProton) {
      LOGP(error, "Could not open histogram {}", Form("%s_antiproton", histname.Data()));
      return;
    }
    hEffPtEtaDeuteron = dynamic_cast<TH2F*>(l->FindObject(Form("%s_deuteron", histname.Data())));
    if (!hEffPtEtaDeuteron) {
      LOGP(error, "Could not open histogram {}", Form("%s_deuteron", histname.Data()));
      return;
    }
    hEffPtEtaAntiDeuteron = dynamic_cast<TH2F*>(l->FindObject(Form("%s_antideuteron", histname.Data())));
    if (!hEffPtEtaAntiDeuteron) {
      LOGP(error, "Could not open histogram {}", Form("%s_antideuteron", histname.Data()));
      return;
    }
    LOGP(info, "Opened histogram {}", Form("%s_proton", histname.Data()));
    LOGP(info, "Opened histogram {}", Form("%s_antiproton", histname.Data()));
    LOGP(info, "Opened histogram {}", Form("%s_deuteron", histname.Data()));
    LOGP(info, "Opened histogram {}", Form("%s_antideuteron", histname.Data()));
  }

  template <typename TParticles>
  float getMCMultiplicity(TParticles const& particles)
  {
    float nCharged = 0.;
    for (const auto& mcParticle : particles) {

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }

      if (std::abs(mcParticle.eta()) > 1.0f) {
        continue;
      }

      TParticlePDG* p = pdgDB->GetParticle(mcParticle.pdgCode());
      if (std::abs(p->Charge()) > 1E-3) {
        nCharged++;
      }
    }

    registry.fill(HIST("hMult"), nCharged);
    return nCharged;
  }

  void processSameEvent(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks)
  {

    registry.fill(HIST("hNEvents"), 0.5);
    registry.fill(HIST("hMult"), collision.mult());

    for (const auto& track : tracks) {

      if (track.tpcFractionSharedCls() > maxtpcSharedCls) {
        continue;
      }
      if (track.itsNCls() < minitsNCls) {
        continue;
      }

      if (isProton(track, +1)) {
        registry.fill(HIST("hPrDCAxy"), track.dcaXY(), track.pt());
      }
      if (isProton(track, -1)) {
        registry.fill(HIST("hAntiPrDCAxy"), track.dcaXY(), track.pt());
      }
      if (isDeuteron(track, +1)) {
        registry.fill(HIST("hDeDCAxy"), track.dcaXY(), track.pt());
      }
      if (isDeuteron(track, -1)) {
        registry.fill(HIST("hAntiDeDCAxy"), track.dcaXY(), track.pt());
      }

      if (!applyDCAcut(track)) {
        continue;
      }

      if (doQA) {
        registryQa.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        registryQa.fill(HIST("QA/hTPCSharedClusters"), track.tpcFractionSharedCls());
        registryQa.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        registryQa.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        registryQa.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        registryQa.fill(HIST("QA/hDCAxy"), track.dcaXY(), track.pt());
        registryQa.fill(HIST("QA/hDCAz"), track.dcaZ(), track.pt());
        registryQa.fill(HIST("QA/TPCChi2VsPZ"), track.tpcInnerParam() / track.sign(), track.tpcChi2NCl());
        registryQa.fill(HIST("QA/hVtxZ_trk"), collision.posZ());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_El"), track.pt() * track.sign(), track.tpcNSigmaEl());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        registryQa.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        registryQa.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
        registryQa.fill(HIST("QA/hnSigmaITSVsPt_Pr"), track.pt() * track.sign(), track.itsNSigmaPr());
        registryQa.fill(HIST("QA/hnSigmaITSVsPt_De"), track.pt() * track.sign(), track.itsNSigmaDe());
        registryQa.fill(HIST("QA/h2dTPCTOF_AntiPr"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        registryQa.fill(HIST("QA/h2dTPCTOF_Pr"), track.tpcNSigmaPr(), track.tofNSigmaPr());

        if (isProton(track, -1)) {
          registryQa.fill(HIST("QA/hEtaAntiPr"), track.eta());
          registryQa.fill(HIST("QA/hPhiAntiPr"), track.phi());
          registryQa.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
          registryQa.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          registryQa.fill(HIST("QA/hnSigmaITSVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.itsNSigmaPr());
          registryQa.fill(HIST("QA/h2dTPCTOF_AntiPr_AfterSel"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        }
        if (isProton(track, +1)) {
          registryQa.fill(HIST("QA/hEtaPr"), track.eta());
          registryQa.fill(HIST("QA/hPhiPr"), track.phi());
          registryQa.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
          registryQa.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          registryQa.fill(HIST("QA/hnSigmaITSVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.itsNSigmaPr());
          registryQa.fill(HIST("QA/h2dTPCTOF_Pr_AfterSel"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        }
        if (isDeuteron(track, -1)) {
          registryQa.fill(HIST("QA/hEtaAntiDe"), track.eta());
          registryQa.fill(HIST("QA/hPhiAntiDe"), track.phi());
          registryQa.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          registryQa.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
          registryQa.fill(HIST("QA/hnSigmaITSVsPt_De_AfterSel"), track.pt() * track.sign(), track.itsNSigmaDe());
        }
        if (isDeuteron(track, +1)) {
          registryQa.fill(HIST("QA/hEtaDe"), track.eta());
          registryQa.fill(HIST("QA/hPhiDe"), track.phi());
          registryQa.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          registryQa.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
          registryQa.fill(HIST("QA/hnSigmaITSVsPt_De_AfterSel"), track.pt() * track.sign(), track.itsNSigmaDe());
        }
      }
    }

    pair->SetMagField1(collision.magField());
    pair->SetMagField2(collision.magField());

    if (mode == kPbarPbar || mode == kPP) { // Identical particle combinations

      for (const auto& [part0, part1] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {

        if (part0.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part0.itsNCls() < minitsNCls) {
          continue;
        }
        if (part1.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part1.itsNCls() < minitsNCls) {
          continue;
        }

        if (!applyDCAcut(part0)) {
          continue;
        }
        if (!applyDCAcut(part1)) {
          continue;
        }

        // remove tracks outside pt bins
        if (part0.pt() < pTBins.value.at(0) || part0.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }
        if (part1.pt() < pTBins.value.at(0) || part1.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }

        // mode 6
        if (mode == kPP) {
          if (!isProton(part0, +1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        // mode 5
        if (mode == kPbarPbar) {
          if (!isProton(part0, -1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }

        fillHistograms(part0, part1, false, true);
      }

    } else {

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

        if (part0.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part0.itsNCls() < minitsNCls) {
          continue;
        }
        if (part1.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part1.itsNCls() < minitsNCls) {
          continue;
        }

        if (!applyDCAcut(part0)) {
          continue;
        }
        if (!applyDCAcut(part1)) {
          continue;
        }

        // remove tracks outside pt bins
        if (part0.pt() < pTBins.value.at(0) || part0.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }
        if (part1.pt() < pTBins.value.at(0) || part1.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }

        // modes 0,1,2,3,4,7
        if (mode == kDbarPbar) {
          if (!isDeuteron(part0, -1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kDP) {
          if (!isDeuteron(part0, +1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kDbarP) {
          if (!isDeuteron(part0, -1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kDPbar) {
          if (!isDeuteron(part0, +1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kPbarP) {
          if (!isProton(part0, -1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kPPbar) {
          if (!isProton(part0, +1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }

        fillHistograms(part0, part1, false, false);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processSameEvent, "processSameEvent", true);

  void processSameEventEvSel(FilteredCollisionsExtra::iterator const& collision, FilteredTracks const& tracks)
  {

    registry.fill(HIST("hNEvents"), 0.5);
    registry.fill(HIST("hMult"), collision.mult());

    for (const auto& track : tracks) {

      if (removeSameBunchPileup && !track.template singleCollSel_as<FilteredCollisionsExtra>().isNoSameBunchPileup()) {
        continue;
      }

      if (track.tpcFractionSharedCls() > maxtpcSharedCls) {
        continue;
      }
      if (track.itsNCls() < minitsNCls) {
        continue;
      }

      if (isProton(track, +1)) {
        registry.fill(HIST("hPrDCAxy"), track.dcaXY(), track.pt());
      }
      if (isProton(track, -1)) {
        registry.fill(HIST("hAntiPrDCAxy"), track.dcaXY(), track.pt());
      }
      if (isDeuteron(track, +1)) {
        registry.fill(HIST("hDeDCAxy"), track.dcaXY(), track.pt());
      }
      if (isDeuteron(track, -1)) {
        registry.fill(HIST("hAntiDeDCAxy"), track.dcaXY(), track.pt());
      }

      if (!applyDCAcut(track)) {
        continue;
      }

      if (doQA) {
        registryQa.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        registryQa.fill(HIST("QA/hTPCSharedClusters"), track.tpcFractionSharedCls());
        registryQa.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        registryQa.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        registryQa.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        registryQa.fill(HIST("QA/hDCAxy"), track.dcaXY(), track.pt());
        registryQa.fill(HIST("QA/hDCAz"), track.dcaZ(), track.pt());
        registryQa.fill(HIST("QA/TPCChi2VsPZ"), track.tpcInnerParam() / track.sign(), track.tpcChi2NCl());
        registryQa.fill(HIST("QA/hVtxZ_trk"), collision.posZ());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_El"), track.pt() * track.sign(), track.tpcNSigmaEl());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        registryQa.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        registryQa.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
        registryQa.fill(HIST("QA/hnSigmaITSVsPt_Pr"), track.pt() * track.sign(), track.itsNSigmaPr());
        registryQa.fill(HIST("QA/hnSigmaITSVsPt_De"), track.pt() * track.sign(), track.itsNSigmaDe());
        registryQa.fill(HIST("QA/h2dTPCTOF_AntiPr"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        registryQa.fill(HIST("QA/h2dTPCTOF_Pr"), track.tpcNSigmaPr(), track.tofNSigmaPr());

        if (isProton(track, -1)) {
          registryQa.fill(HIST("QA/hEtaAntiPr"), track.eta());
          registryQa.fill(HIST("QA/hPhiAntiPr"), track.phi());
          registryQa.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
          registryQa.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          registryQa.fill(HIST("QA/hnSigmaITSVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.itsNSigmaPr());
          registryQa.fill(HIST("QA/h2dTPCTOF_AntiPr_AfterSel"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        }
        if (isProton(track, +1)) {
          registryQa.fill(HIST("QA/hEtaPr"), track.eta());
          registryQa.fill(HIST("QA/hPhiPr"), track.phi());
          registryQa.fill(HIST("QA/hnSigmaTOFVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tofNSigmaPr());
          registryQa.fill(HIST("QA/hnSigmaTPCVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaPr());
          registryQa.fill(HIST("QA/hnSigmaITSVsPt_Pr_AfterSel"), track.pt() * track.sign(), track.itsNSigmaPr());
          registryQa.fill(HIST("QA/h2dTPCTOF_Pr_AfterSel"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        }
        if (isDeuteron(track, -1)) {
          registryQa.fill(HIST("QA/hEtaAntiDe"), track.eta());
          registryQa.fill(HIST("QA/hPhiAntiDe"), track.phi());
          registryQa.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          registryQa.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
          registryQa.fill(HIST("QA/hnSigmaITSVsPt_De_AfterSel"), track.pt() * track.sign(), track.itsNSigmaDe());
        }
        if (isDeuteron(track, +1)) {
          registryQa.fill(HIST("QA/hEtaDe"), track.eta());
          registryQa.fill(HIST("QA/hPhiDe"), track.phi());
          registryQa.fill(HIST("QA/hnSigmaTOFVsPt_De_AfterSel"), track.pt() * track.sign(), track.tofNSigmaDe());
          registryQa.fill(HIST("QA/hnSigmaTPCVsPt_De_AfterSel"), track.pt() * track.sign(), track.tpcNSigmaDe());
          registryQa.fill(HIST("QA/hnSigmaITSVsPt_De_AfterSel"), track.pt() * track.sign(), track.itsNSigmaDe());
        }
      }
    }

    pair->SetMagField1(collision.magField());
    pair->SetMagField2(collision.magField());

    if (mode == kPbarPbar || mode == kPP) { // Identical particle combinations

      for (const auto& [part0, part1] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {

        if (removeSameBunchPileup && !part0.template singleCollSel_as<FilteredCollisionsExtra>().isNoSameBunchPileup()) {
          continue;
        }

        if (part0.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part0.itsNCls() < minitsNCls) {
          continue;
        }
        if (part1.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part1.itsNCls() < minitsNCls) {
          continue;
        }

        if (!applyDCAcut(part0)) {
          continue;
        }
        if (!applyDCAcut(part1)) {
          continue;
        }

        // remove tracks outside pt bins
        if (part0.pt() < pTBins.value.at(0) || part0.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }
        if (part1.pt() < pTBins.value.at(0) || part1.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }

        // mode 6
        if (mode == kPP) {
          if (!isProton(part0, +1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        // mode 5
        if (mode == kPbarPbar) {
          if (!isProton(part0, -1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }

        fillHistograms(part0, part1, false, true);
      }

    } else {

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

        if (removeSameBunchPileup && !part0.template singleCollSel_as<FilteredCollisionsExtra>().isNoSameBunchPileup()) {
          continue;
        }

        if (part0.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part0.itsNCls() < minitsNCls) {
          continue;
        }
        if (part1.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part1.itsNCls() < minitsNCls) {
          continue;
        }

        if (!applyDCAcut(part0)) {
          continue;
        }
        if (!applyDCAcut(part1)) {
          continue;
        }

        // remove tracks outside pt bins
        if (part0.pt() < pTBins.value.at(0) || part0.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }
        if (part1.pt() < pTBins.value.at(0) || part1.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }

        // modes 0,1,2,3,4,7
        if (mode == kDbarPbar) {
          if (!isDeuteron(part0, -1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kDP) {
          if (!isDeuteron(part0, +1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kDbarP) {
          if (!isDeuteron(part0, -1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kDPbar) {
          if (!isDeuteron(part0, +1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kPbarP) {
          if (!isProton(part0, -1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kPPbar) {
          if (!isProton(part0, +1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }

        fillHistograms(part0, part1, false, false);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processSameEventEvSel, "processSameEventEvSel", false);

  void processMixedEvent(FilteredCollisions const& collisions, FilteredTracks const& tracks)
  {

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, collisions, collisions)) {

      // LOGF(info, "Mixed event collisions: (%d, %d) zvtx (%.1f, %.1f) mult (%d, %d)", collision1.globalIndex(), collision2.globalIndex(), collision1.posZ(), collision2.posZ(), collision1.mult(), collision2.mult());

      auto groupPartsOne = tracks.sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision1.globalIndex(), cache);
      auto groupPartsTwo = tracks.sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      const float limit = 1e-4;

      if (std::abs(magFieldTesla1 - magFieldTesla2) > limit) {
        continue;
      }

      pair->SetMagField1(magFieldTesla1);
      pair->SetMagField2(magFieldTesla2);

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        if (part0.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part0.itsNCls() < minitsNCls) {
          continue;
        }
        if (part1.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part1.itsNCls() < minitsNCls) {
          continue;
        }

        if (!applyDCAcut(part0)) {
          continue;
        }
        if (!applyDCAcut(part1)) {
          continue;
        }

        // remove tracks outside pt bins
        if (part0.pt() < pTBins.value.at(0) || part0.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }
        if (part1.pt() < pTBins.value.at(0) || part1.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }

        //{"mode", 0, "0: antid-antip, 1: d-p, 2: antid-p, 3: d-antip, 4: antip-p, 5: antip-antip, 6: p-p, 7: p-antip"};
        if (mode == kDbarPbar) {
          if (!isDeuteron(part0, -1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kDP) {
          if (!isDeuteron(part0, +1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kDbarP) {
          if (!isDeuteron(part0, -1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kDPbar) {
          if (!isDeuteron(part0, +1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kPbarP) {
          if (!isProton(part0, -1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kPbarPbar) {
          if (!isProton(part0, -1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kPP) {
          if (!isProton(part0, +1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kPPbar) {
          if (!isProton(part0, +1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }

        bool isIdentical = false;
        if (mode == kPbarPbar || mode == kPP) {
          isIdentical = true;
        }

        fillHistograms(part0, part1, true, isIdentical);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processMixedEvent, "processMixedEvent", true);

  void processMixedEventEvSel(FilteredCollisionsExtra const& collisions, FilteredTracks const& tracks)
  {

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, collisions, collisions)) {

      // LOGF(info, "Mixed event collisions: (%d, %d) zvtx (%.1f, %.1f) mult (%d, %d)", collision1.globalIndex(), collision2.globalIndex(), collision1.posZ(), collision2.posZ(), collision1.mult(), collision2.mult());

      auto groupPartsOne = tracks.sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision1.globalIndex(), cache);
      auto groupPartsTwo = tracks.sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      const float limit = 1e-4;

      if (std::abs(magFieldTesla1 - magFieldTesla2) > limit) {
        continue;
      }

      pair->SetMagField1(magFieldTesla1);
      pair->SetMagField2(magFieldTesla2);

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        if (removeSameBunchPileup && !part0.template singleCollSel_as<FilteredCollisionsExtra>().isNoSameBunchPileup()) {
          continue;
        }
        if (removeSameBunchPileup && !part1.template singleCollSel_as<FilteredCollisionsExtra>().isNoSameBunchPileup()) {
          continue;
        }

        if (part0.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part0.itsNCls() < minitsNCls) {
          continue;
        }
        if (part1.tpcFractionSharedCls() > maxtpcSharedCls) {
          continue;
        }
        if (part1.itsNCls() < minitsNCls) {
          continue;
        }

        if (!applyDCAcut(part0)) {
          continue;
        }
        if (!applyDCAcut(part1)) {
          continue;
        }

        // remove tracks outside pt bins
        if (part0.pt() < pTBins.value.at(0) || part0.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }
        if (part1.pt() < pTBins.value.at(0) || part1.pt() >= pTBins.value.at(nBinspT)) {
          continue;
        }

        //{"mode", 0, "0: antid-antip, 1: d-p, 2: antid-p, 3: d-antip, 4: antip-p, 5: antip-antip, 6: p-p, 7: p-antip"};
        if (mode == kDbarPbar) {
          if (!isDeuteron(part0, -1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kDP) {
          if (!isDeuteron(part0, +1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kDbarP) {
          if (!isDeuteron(part0, -1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kDPbar) {
          if (!isDeuteron(part0, +1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kPbarP) {
          if (!isProton(part0, -1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kPbarPbar) {
          if (!isProton(part0, -1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }
        if (mode == kPP) {
          if (!isProton(part0, +1)) {
            continue;
          }
          if (!isProton(part1, +1)) {
            continue;
          }
        }
        if (mode == kPPbar) {
          if (!isProton(part0, +1)) {
            continue;
          }
          if (!isProton(part1, -1)) {
            continue;
          }
        }

        bool isIdentical = false;
        if (mode == kPbarPbar || mode == kPP) {
          isIdentical = true;
        }

        fillHistograms(part0, part1, true, isIdentical);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processMixedEventEvSel, "processMixedEventEvSel", false);

  void processMC(FilteredCollisions const&, FilteredTracksMC const& tracks)
  {
    for (const auto& track : tracks) {
      if (std::abs(track.template singleCollSel_as<FilteredCollisions>().posZ()) > cutzVertex) {
        continue;
      }

      if (track.tpcFractionSharedCls() > maxtpcSharedCls) {
        continue;
      }
      if (track.itsNCls() < minitsNCls) {
        continue;
      }

      if (isProton(track, +1) && track.pdgCode() == PDG_t::kProton) {
        registry.fill(HIST("hPrDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == kPrimary) {
          registry.fill(HIST("hPrimPrDCAxy"), track.dcaXY(), track.pt());
        }
        if (track.origin() == kWeakDecay) {
          registry.fill(HIST("hSecWeakPrDCAxy"), track.dcaXY(), track.pt());
        }
        if (track.origin() == kMaterial) {
          registry.fill(HIST("hSecMatPrDCAxy"), track.dcaXY(), track.pt());
        }
      }
      if (isProton(track, -1) && track.pdgCode() == -PDG_t::kProton) {
        registry.fill(HIST("hAntiPrDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == kPrimary) {
          registry.fill(HIST("hPrimAntiPrDCAxy"), track.dcaXY(), track.pt());
        }
        if (track.origin() == kWeakDecay) {
          registry.fill(HIST("hSecWeakAntiPrDCAxy"), track.dcaXY(), track.pt());
        }
        if (track.origin() == kMaterial) {
          registry.fill(HIST("hSecMatAntiPrDCAxy"), track.dcaXY(), track.pt());
        }
      }
      if (isDeuteron(track, +1) && track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("hDeDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == kPrimary) {
          registry.fill(HIST("hPrimDeDCAxy"), track.dcaXY(), track.pt());
        }
        if (track.origin() == kWeakDecay) {
          registry.fill(HIST("hSecWeakDeDCAxy"), track.dcaXY(), track.pt());
        }
        if (track.origin() == kMaterial) {
          registry.fill(HIST("hSecMatDeDCAxy"), track.dcaXY(), track.pt());
        }
      }
      if (isDeuteron(track, -1) && track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
        registry.fill(HIST("hAntiDeDCAxy"), track.dcaXY(), track.pt());
        if (track.origin() == kPrimary) {
          registry.fill(HIST("hPrimAntiDeDCAxy"), track.dcaXY(), track.pt());
        }
        if (track.origin() == kWeakDecay) {
          registry.fill(HIST("hSecWeakAntiDeDCAxy"), track.dcaXY(), track.pt());
        }
        if (track.origin() == kMaterial) {
          registry.fill(HIST("hSecMatAntiDeDCAxy"), track.dcaXY(), track.pt());
        }
      }

      if (!applyDCAcut(track)) {
        continue;
      }

      // Keep only protons and deuterons
      // if (std::abs(track.pdgCode()) != PDG_t::kProton && std::abs(track.pdgCode()) != o2::constants::physics::Pdg::kDeuteron)
      // continue;

      if (doQA) {
        registryQa.fill(HIST("QA/hTPCnClusters"), track.tpcNClsFound());
        registryQa.fill(HIST("QA/hTPCSharedClusters"), track.tpcFractionSharedCls());
        registryQa.fill(HIST("QA/hTPCchi2"), track.tpcChi2NCl());
        registryQa.fill(HIST("QA/hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        registryQa.fill(HIST("QA/hITSchi2"), track.itsChi2NCl());
        registryQa.fill(HIST("QA/hDCAxy"), track.dcaXY(), track.pt());
        registryQa.fill(HIST("QA/hDCAz"), track.dcaZ(), track.pt());
        registryQa.fill(HIST("QA/hVtxZ_trk"), track.template singleCollSel_as<FilteredCollisions>().posZ());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_El"), track.pt() * track.sign(), track.tpcNSigmaEl());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_Pr"), track.pt() * track.sign(), track.tpcNSigmaPr());
        registryQa.fill(HIST("QA/hnSigmaTPCVsPt_De"), track.pt() * track.sign(), track.tpcNSigmaDe());
        registryQa.fill(HIST("QA/hnSigmaTOFVsPt_Pr"), track.pt() * track.sign(), track.tofNSigmaPr());
        registryQa.fill(HIST("QA/hnSigmaTOFVsPt_De"), track.pt() * track.sign(), track.tofNSigmaDe());
      }

      bool isPr = (isProton(track, +1) && track.pdgCode() == PDG_t::kProton);
      bool isAntiPr = (isProton(track, -1) && track.pdgCode() == -PDG_t::kProton);
      bool isDe = (isDeuteron(track, +1) && track.pdgCode() == o2::constants::physics::Pdg::kDeuteron);
      bool isAntiDe = (isDeuteron(track, -1) && track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron);

      if (isPr) {
        registry.fill(HIST("hPrimSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * +1);
        if (track.origin() == kWeakDecay || track.origin() == kMaterial) { // secondaries
          registry.fill(HIST("hSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * +1);
        }
      }
      if (isAntiPr) {
        registry.fill(HIST("hPrimSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        if (track.origin() == kWeakDecay || track.origin() == kMaterial) {
          registry.fill(HIST("hSec_EtaPhiPt_Proton"), track.eta(), track.phi(), track.pt() * -1);
        }
      }

      if (track.origin() != 0) {
        continue;
      }

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
      if (isProton(track, +1)) {
        registry.fill(HIST("hDenominatorPurity_Proton"), track.pt());
      }
      if (isProton(track, -1)) {
        registry.fill(HIST("hDenominatorPurity_Proton"), track.pt() * -1);
      }
      if (isDeuteron(track, +1)) {
        registry.fill(HIST("hDenominatorPurity_Deuteron"), track.pt());
      }
      if (isDeuteron(track, -1)) {
        registry.fill(HIST("hDenominatorPurity_Deuteron"), track.pt() * -1);
      }

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
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.pdgCode() == PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl"), track.pt());
          registry.fill(HIST("hReco_Pt_Proton_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
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
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.pdgCode() == -PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPC_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPC_or_TOF"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.pdgCode() == -PDG_t::kProton) {
          registry.fill(HIST("hNumeratorPurity_Proton_TPCEl"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Proton_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
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
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC_or_TOF"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.pdgCode() == o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl"), track.pt());
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
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
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPC_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPC_or_TOF"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.pdgCode() == -o2::constants::physics::Pdg::kDeuteron) {
          registry.fill(HIST("hNumeratorPurity_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
          registry.fill(HIST("hReco_Pt_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
        }

        // Denominators
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPC"), track.pt());
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCTOF"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl_or_TOF"), track.pt());
        }

        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPC"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCTOF"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPC_or_TOF"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaPr()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElPr && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaPr()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElPr && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaPr()) < nsigmaTPC && std::abs(track.tofNSigmaPr()) < nsigmaTOF)) &&
            track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Proton_TPCEl_or_TOF"), track.pt() * -1);
        }

        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCTOF"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC_or_TOF"), track.pt());
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl"), track.pt());
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() > 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl_or_TOF"), track.pt());
        }

        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCTOF"), track.pt() * -1);
        }
        if ((
              (std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.beta() < BetahasTOFthr) ||
              (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPC_or_TOF"), track.pt() * -1);
        }
        if (std::abs(track.tpcNSigmaDe()) < nsigmaTPC &&
            track.tpcNSigmaEl() >= nsigmaElDe && track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl"), track.pt() * -1);
        }
        if (((std::abs(track.tpcNSigmaDe()) < nsigmaTPC && track.tpcNSigmaEl() >= nsigmaElDe && track.beta() < BetahasTOFthr) ||
             (track.beta() > BetahasTOFthr && std::abs(track.tpcNSigmaDe()) < nsigmaTPC && std::abs(track.tofNSigmaDe()) < nsigmaTOF)) &&
            track.sign() < 0) {
          registry.fill(HIST("hDenominatorPurity_Deuteron_TPCEl_or_TOF"), track.pt() * -1);
        }
      }
    } // track
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processMC, "processMC", false);

  void processSameEventGen(SimCollisions::iterator const&, SimParticles const& mcParticles)
  {

    registry.fill(HIST("Generated/hNEventsMC"), 0.5);

    for (const auto& particle : mcParticles) {
      auto fillGeneratedQa = [this, &particle](const float binPosition) {
        switch (particle.pdgCode()) {
          case PDG_t::kProton:
            registry.fill(HIST("Generated/hQAProtons"), binPosition);
            return true;
          case -PDG_t::kProton:
            registry.fill(HIST("Generated/hQAAntiProtons"), binPosition);
            return true;
          case o2::constants::physics::Pdg::kDeuteron:
            registry.fill(HIST("Generated/hQADeuterons"), binPosition);
            return true;
          case -o2::constants::physics::Pdg::kDeuteron:
            registry.fill(HIST("Generated/hQAAntiDeuterons"), binPosition);
            return true;
          default:
            return false;
        }
      };
      if (!fillGeneratedQa(0.5)) {
        continue;
      }
      if (isPrim && !particle.isPhysicalPrimary()) {
        continue;
      }

      fillGeneratedQa(1.5);

      if (std::abs(particle.y()) < yRap) {
        switch (particle.pdgCode()) {
          case o2::constants::physics::Pdg::kDeuteron:
            registry.fill(HIST("Generated/hDeuteronsVsPt"), particle.pt());
            break;
          case -o2::constants::physics::Pdg::kDeuteron:
            registry.fill(HIST("Generated/hAntiDeuteronsVsPt"), particle.pt());
            break;
          default:
            break;
        }
      }

      if (std::abs(particle.eta()) > etaCut) {
        continue;
      }
      fillGeneratedQa(2.5);

      switch (particle.pdgCode()) {
        case PDG_t::kProton:
          registry.fill(HIST("hGen_EtaPhiPt_Proton"), particle.eta(), particle.phi(), particle.pt());
          break;
        case -PDG_t::kProton:
          registry.fill(HIST("hGen_EtaPhiPt_Proton"), particle.eta(), particle.phi(), -1. * particle.pt());
          break;
        case o2::constants::physics::Pdg::kDeuteron:
          registry.fill(HIST("hGen_EtaPhiPt_Deuteron"), particle.eta(), particle.phi(), particle.pt());
          break;
        case -o2::constants::physics::Pdg::kDeuteron:
          registry.fill(HIST("hGen_EtaPhiPt_Deuteron"), particle.eta(), particle.phi(), -1. * particle.pt());
          break;
        default:
          LOG(fatal) << "Unhandled PDG code, should not happen, check the code!" << particle.pdgCode();
          break;
      }
    }

    if (mode == kPbarPbar || mode == kPP) { // Identical particle combinations

      for (const auto& [part0, part1] : combinations(CombinationsStrictlyUpperIndexPolicy(mcParticles, mcParticles))) {

        if (isPrim && !part0.isPhysicalPrimary()) {
          continue;
        }
        if (isPrim && !part1.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(part0.eta()) > etaCut) {
          continue;
        }
        if (std::abs(part1.eta()) > etaCut) {
          continue;
        }

        // mode 6
        if (mode == kPP) {
          if (part0.pdgCode() != PDG_t::kProton) {
            continue;
          }
          if (part1.pdgCode() != PDG_t::kProton) {
            continue;
          }
        }
        // mode 5
        if (mode == kPbarPbar) {
          if (part0.pdgCode() != -PDG_t::kProton) {
            continue;
          }
          if (part1.pdgCode() != -PDG_t::kProton) {
            continue;
          }
        }

        fillHistogramsGen(part0, part1, false);
      }

    } else {

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(mcParticles, mcParticles))) {

        if (isPrim && !part0.isPhysicalPrimary()) {
          continue;
        }
        if (isPrim && !part1.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(part0.eta()) > etaCut) {
          continue;
        }
        if (std::abs(part1.eta()) > etaCut) {
          continue;
        }

        if (mode == kDbarPbar) {
          if (part0.pdgCode() != -o2::constants::physics::Pdg::kDeuteron) {
            continue;
          }
          if (part1.pdgCode() != -PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kDP) {
          if (part0.pdgCode() != o2::constants::physics::Pdg::kDeuteron) {
            continue;
          }
          if (part1.pdgCode() != PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kDbarP) {
          if (part0.pdgCode() != -o2::constants::physics::Pdg::kDeuteron) {
            continue;
          }
          if (part1.pdgCode() != PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kDPbar) {
          if (part0.pdgCode() != o2::constants::physics::Pdg::kDeuteron) {
            continue;
          }
          if (part1.pdgCode() != -PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kPbarP) {
          if (part0.pdgCode() != -PDG_t::kProton) {
            continue;
          }
          if (part1.pdgCode() != PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kPPbar) {
          if (part0.pdgCode() != PDG_t::kProton) {
            continue;
          }
          if (part1.pdgCode() != -PDG_t::kProton) {
            continue;
          }
        }

        fillHistogramsGen(part0, part1, false);
      }
    }
  }
  PROCESS_SWITCH(HadronNucleiCorrelation, processSameEventGen, "processSameEventGen", false);

  Preslice<SimParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;

  void processMixedEventGen(SimCollisions const& mcCollisions, SimParticles const& mcParticles)
  {

    auto getMultiplicity = [this, &mcParticles](SimCollisions::iterator const& collision) {
      auto particlesPerCol = mcParticles.sliceBy(perMcCollision, collision.globalIndex());
      auto multiplicity = getMCMultiplicity(particlesPerCol);
      return multiplicity;
    };

    using BinningTypeMC = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::mccollision::PosZ, decltype(getMultiplicity)>;
    BinningTypeMC colBinningGen{{getMultiplicity}, {confVtxBins, confMultBins}, true};

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningGen, 5, -1, mcCollisions, mcCollisions)) {

      auto groupPartsOne = mcParticles.sliceBy(perMcCollision, collision1.globalIndex());
      auto groupPartsTwo = mcParticles.sliceBy(perMcCollision, collision2.globalIndex());

      // LOGF(info, "Mixed event collisions: (%d, %d) zvtx (%.1f, %.1f) mult (%.1f, %.1f)", collision1.globalIndex(), collision2.globalIndex(), collision1.posZ(), collision2.posZ(), getMCMultiplicity(groupPartsOne), getMCMultiplicity(groupPartsTwo));

      for (const auto& [part0, part1] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        if (isPrim && !part0.isPhysicalPrimary()) {
          continue;
        }
        if (isPrim && !part1.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(part0.eta()) > etaCut) {
          continue;
        }
        if (std::abs(part1.eta()) > etaCut) {
          continue;
        }

        if (mode == kDbarPbar) {
          if (part0.pdgCode() != -o2::constants::physics::Pdg::kDeuteron) {
            continue;
          }
          if (part1.pdgCode() != -PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kDP) {
          if (part0.pdgCode() != o2::constants::physics::Pdg::kDeuteron) {
            continue;
          }
          if (part1.pdgCode() != PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kDbarP) {
          if (part0.pdgCode() != -o2::constants::physics::Pdg::kDeuteron) {
            continue;
          }
          if (part1.pdgCode() != PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kDPbar) {
          if (part0.pdgCode() != o2::constants::physics::Pdg::kDeuteron) {
            continue;
          }
          if (part1.pdgCode() != -PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kPbarP) {
          if (part0.pdgCode() != -PDG_t::kProton) {
            continue;
          }
          if (part1.pdgCode() != PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kPbarPbar) {
          if (part0.pdgCode() != -PDG_t::kProton) {
            continue;
          }
          if (part1.pdgCode() != -PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kPP) {
          if (part0.pdgCode() != PDG_t::kProton) {
            continue;
          }
          if (part1.pdgCode() != PDG_t::kProton) {
            continue;
          }
        }
        if (mode == kPPbar) {
          if (part0.pdgCode() != PDG_t::kProton) {
            continue;
          }
          if (part1.pdgCode() != -PDG_t::kProton) {
            continue;
          }
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
