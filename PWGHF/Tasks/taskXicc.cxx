// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskXicc.cxx
/// \brief Ξcc±± analysis task
/// \note Inspired from taskLc.cxx
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Jinjoo Seo <jin.joo.seo@cern.ch>, Inha University

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_xicc;
//using namespace o2::aod::hf_cand_prong3;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Ξcc±± analysis task
struct HfTaskXicc {
  HistogramRegistry registry{
    "registry",
    {{"hptcand", "#Xi^{++}_{cc}-candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hptprong0", "#Xi^{++}_{cc}-candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hptprong1", "#Xi^{++}_{cc}-candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}}}};

  Configurable<int> d_selectionFlagXicc{"d_selectionFlagXicc", 1, "Selection Flag for Xicc"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_xicc_topkpipi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hmass", "#Xi^{++}_{cc} candidates;inv. mass (p K #pi #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{400, 3.2, 4.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "#Xi^{++}_{cc} candidates;decay length (cm);entries", {HistType::kTH2F, {{500, 0., 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCA", "#Xi^{++}_{cc} candidates;chi2 PCA (cm);entries", {HistType::kTH2F, {{500, 0., 0.01}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "#Xi^{++}_{cc} candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "#Xi^{++}_{cc} candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "#Xi^{++}_{cc} candidates;proper lifetime (#Xi^{++}_{cc}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "#Xi^{++}_{cc} candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{1000, -1.0, 1.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("habsCPA", "#Xi^{++}_{cc} candidates;abs. cosine of pointing angle;entries", {HistType::kTH2F, {{200, 0.8, 1.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "#Xi^{++}_{cc} candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{250, -5., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hY", "#Xi^{++}_{cc} candidates;candidate rapidity;entries", {HistType::kTH2F, {{250, -5., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hselectionstatus", "#Xi^{++}_{cc} candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr0", "#Xi^{++}_{cc} candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{200, 0, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr1", "#Xi^{++}_{cc} candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{200, 0, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "#Xi^{++}_{cc} candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_xicc::isSelXiccToPKPiPi >= d_selectionFlagXicc);
  void process(soa::Filtered<soa::Join<aod::HfCandXicc, aod::HFSelXiccToPKPiPiCandidate>> const& candidates)
  //void process(aod::HfCandXicc const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::XiccToXicPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YXicc(candidate)) > cutYCandMax) {
        continue;
      }
      registry.fill(HIST("hmass"), InvMassXiccToXicPi(candidate), candidate.pt()); //FIXME need to consider the two mass hp
      registry.fill(HIST("hptcand"), candidate.pt());
      registry.fill(HIST("hptprong0"), candidate.ptProng0());
      registry.fill(HIST("hptprong1"), candidate.ptProng1());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hChi2PCA"), candidate.chi2PCA(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), std::abs(candidate.impactParameter0()), candidate.pt());
      registry.fill(HIST("hd0Prong1"), std::abs(candidate.impactParameter1()), candidate.pt());
      registry.fill(HIST("hCt"), CtXicc(candidate), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("habsCPA"), std::abs(candidate.cpa()), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hY"), YXicc(candidate), candidate.pt());
      registry.fill(HIST("hselectionstatus"), candidate.isSelXiccToPKPiPi(), candidate.pt());
      registry.fill(HIST("hImpParErr0"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr1"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
    }
  }
};

/// Fills MC histograms.
struct HfTaskXiccMc {
  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "#Xi^{++}_{cc} candidates (rec. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtRecBg", "#Xi^{++}_{cc} candidates (rec. unmatched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtGen", "#Xi^{++}_{cc} MC particles (matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtGenSig", "#Xi^{++}_{cc} candidates (rec. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hEtaRecSig", "#Xi^{++}_{cc} candidates (rec. matched);#it{#eta};entries", {HistType::kTH1F, {{250, -5., 5.}}}},
     {"hEtaRecBg", "#Xi^{++}_{cc} candidates (rec. unmatched);#it{#eta};entries", {HistType::kTH1F, {{250, -5., 5.}}}},
     {"hEtaGen", "#Xi^{++}_{cc} MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{250, -5., 5.}}}},
     {"hYRecSig", "#Xi^{++}_{cc} candidates (rec. matched);#it{y};entries", {HistType::kTH1F, {{250, -5., 5.}}}},
     {"hYRecBg", "#Xi^{++}_{cc} candidates (rec. unmatched);#it{y};entries", {HistType::kTH1F, {{250, -5., 5.}}}},
     {"hYGen", "#Xi^{++}_{cc} MC particles (matched);#it{y};entries", {HistType::kTH1F, {{250, -5., 5.}}}},
     {"hPtvsEtavsYGen", "#Xi^{++}_{cc} MC particles (matched);#it{p}_{T};#it{#eta};#it{y}", {HistType::kTH3F, {{100, 0., 10.0}, {250, -5., 5.}, {250, -5., 5.}}}}}};

  Configurable<int> d_selectionFlagXicc{"d_selectionFlagXicc", 1, "Selection Flag for Xicc"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_xicc_topkpipi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hmassSig", "#Xi^{++}_{cc} (rec. matched) candidates;inv. mass (p K #pi #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{400, 3.2, 4.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hmassBg", "#Xi^{++}_{cc} (rec. unmatched) candidates;inv. mass (p K #pi #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{400, 3.2, 4.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCASig", "#Xi^{++}_{cc} (rec. matched) candidates;chi2 PCA (cm);entries", {HistType::kTH2F, {{500, 0., 0.01}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCABg", "#Xi^{++}_{cc} (rec. unmatched) candidates;chi2 PCA (cm);entries", {HistType::kTH2F, {{500, 0., 0.01}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthSig", "#Xi^{++}_{cc} (rec. matched) candidates;decay length (cm);entries", {HistType::kTH2F, {{500, 0., 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthBg", "#Xi^{++}_{cc} (rec. unmatched) candidates;decay length (cm);entries", {HistType::kTH2F, {{500, 0., 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0Sig", "#Xi^{++}_{cc} (rec. matched) candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0Bg", "#Xi^{++}_{cc} (rec. unmatched) candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1Sig", "#Xi^{++}_{cc} (rec. matched) candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1Bg", "#Xi^{++}_{cc} (rec. unmatched) candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtSig", "#Xi^{++}_{cc} (rec. matched) candidates;proper lifetime (#Xi_{cc}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtBg", "#Xi^{++}_{cc} (rec. unmatched) candidates;proper lifetime (#Xi_{cc}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPASig", "#Xi^{++}_{cc} (rec. matched) candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{1000, -1.0, 1.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPABg", "#Xi^{++}_{cc} (rec. unmatched) candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{1000, -1.0, 1.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("habsCPASig", "#Xi^{++}_{cc} (rec. matched) candidates;abs. cosine of pointing angle;entries", {HistType::kTH2F, {{200, 0.8, 1.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("habsCPABg", "#Xi^{++}_{cc} (rec. unmatched) candidates;abs. cosine of pointing angle;entries", {HistType::kTH2F, {{200, 0.8, 1.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaSig", "#Xi^{++}_{cc} (rec. matched) candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{250, -5., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaBg", "#Xi^{++}_{cc} (rec. unmatched) candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{250, -5., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYSig", "#Xi^{++}_{cc} (rec. matched) candidates;candidate rapidity;entries", {HistType::kTH2F, {{250, -5., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYBg", "#Xi^{++}_{cc} (rec. unmatched) candidates;candidate rapidity;entries", {HistType::kTH2F, {{250, -5., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hselectionstatusSig", "#Xi^{++}_{cc} (rec. matched) candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hselectionstatusBg", "#Xi^{++}_{cc} (rec. unmatched) candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr0Sig", "#Xi^{++}_{cc} (rec. matched) candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{200, 0, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr0Bg", "#Xi^{++}_{cc} (rec. unmatched) candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{200, 0, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr1Sig", "#Xi^{++}_{cc} (rec. matched) candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{200, 0, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr1Bg", "#Xi^{++}_{cc} (rec. unmatched) candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{200, 0, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // resolutions
    registry.add("hXSecVtxPosDiff", "#Xi^{++}_{cc} (rec. matched) candidates;x-axis sec. vertex pos. reco - gen (cm);entries", {HistType::kTH2F, {{400, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYSecVtxPosDiff", "#Xi^{++}_{cc} (rec. matched) candidates;y-axis sec. vertex pos. reco - gen (cm);entries", {HistType::kTH2F, {{400, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hZSecVtxPosDiff", "#Xi^{++}_{cc} (rec. matched) candidates;z-axis sec. vertex pos. reco - gen (cm);entries", {HistType::kTH2F, {{400, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecGenDiff", "#Xi^{++}_{cc} (rec. matched) candidates;pt reco - gen;entries (GeV/#it{c}})", {HistType::kTH2F, {{400, -1.0, 1.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    //debug
    registry.add("hDebugMCmatching", "#Xi^{++}_{cc} (rec. matched) candidates;debug MC matching bitmap;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // Check Y dependence (To be removed)
    registry.add("hmassSigvsPtvsY", "#Xi^{++}_{cc} (rec. matched) candidates;inv. mass (p K #pi #pi) (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{400, 3.2, 4.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hmassBgvsPtvsY", "#Xi^{++}_{cc} (rec. unmatched) candidates;inv. mass (p K #pi #pi) (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{400, 3.2, 4.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hChi2PCASigvsPtvsY", "#Xi^{++}_{cc} (rec. matched) candidates;chi2 PCA (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{500, 0., 0.01}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hChi2PCABgvsPtvsY", "#Xi^{++}_{cc} (rec. unmatched) candidates;chi2 PCA (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{500, 0., 0.01}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hDecLengthSigvsPtvsY", "#Xi^{++}_{cc} (rec. matched) candidates;decay length (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{500, 0., 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hDecLengthBgvsPtvsY", "#Xi^{++}_{cc} (rec. unmatched) candidates;decay length (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{500, 0., 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hd0Prong0SigvsPtvsY", "#Xi^{++}_{cc} (rec. matched) candidates;prong 0 DCAxy to prim. vertex (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hd0Prong0BgvsPtvsY", "#Xi^{++}_{cc} (rec. unmatched) candidates;prong 0 DCAxy to prim. vertex (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hd0Prong1SigvsPtvsY", "#Xi^{++}_{cc} (rec. matched) candidates;prong 1 DCAxy to prim. vertex (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hd0Prong1BgvsPtvsY", "#Xi^{++}_{cc} (rec. unmatched) candidates;prong 1 DCAxy to prim. vertex (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{200, -0.02, 0.02}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hCtSigvsPtvsY", "#Xi^{++}_{cc} (rec. matched) candidates;proper lifetime (#Xi_{cc}) * #it{c} (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{100, 0., 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hCtBgvsPtvsY", "#Xi^{++}_{cc} (rec. unmatched) candidates;proper lifetime (#Xi_{cc}) * #it{c} (cm); #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{100, 0., 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hCPASigvsPtvsY", "#Xi^{++}_{cc} (rec. matched) candidates;cosine of pointing angle; #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{1000, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("hCPABgvsPtvsY", "#Xi^{++}_{cc} (rec. unmatched) candidates;cosine of pointing angle; #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{1000, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("habsCPASigvsPtvsY", "#Xi^{++}_{cc} (rec. matched) candidates;abs. cosine of pointing angle;  #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{200, 0.8, 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
    registry.add("habsCPABgvsPtvsY", "#Xi^{++}_{cc} (rec. unmatched) candidates;abs. cosine of pointing angle; #it{p}_{T} (GeV/#it{c}); #it{y}", {HistType::kTH3F, {{200, 0.8, 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {250, -5., 5.}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_xicc::isSelXiccToPKPiPi >= d_selectionFlagXicc);
  //void process(soa::Filtered<soa::Join<aod::HfCandXicc, aod::HFSelXiccToPKPiPiCandidate>> const& candidates)
  void process(soa::Filtered<soa::Join<aod::HfCandXicc, aod::HFSelXiccToPKPiPiCandidate, aod::HfCandXiccMCRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandXiccMCGen> const& particlesMC, aod::BigTracksMC const& tracks)
  {
    // MC rec.
    //Printf("MC Candidates: %d", candidates.size());
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::XiccToXicPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YXicc(candidate)) > cutYCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMCMatchRec()) == 1 << DecayType::XiccToXicPi) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.index1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandXiccMCGen>>(), 4422, true);
        auto particleXicc = particlesMC.iteratorAt(indexMother);
        auto particleXic = particlesMC.iteratorAt(particleXicc.daughter0Id());
        /*
        auto daughter1 = particlesMC.iteratorAt(particleXicc.daughter1());
        auto p0xic = particlesMC.iteratorAt(particleXic.daughter0());
        auto p1xic = particlesMC.iteratorAt(particleXic.daughter0()+1);
        auto p2xic = particlesMC.iteratorAt(particleXic.daughter1());
        LOGF(info, "mother pdg %d", particleXicc.pdgCode());
        LOGF(info, "Xic pdg %d", particleXic.pdgCode());
        LOGF(info, "Xic prong 0 pdg %d", p0xic.pdgCode());
        LOGF(info, "Xic prong 1 pdg %d", p1xic.pdgCode());
        LOGF(info, "Xic prong 2 pdg %d", p2xic.pdgCode());
        LOGF(info, "Pion pdg %d", daughter1.pdgCode());
*/
        registry.fill(HIST("hPtGenSig"), particleXicc.pt()); // gen. level pT
        registry.fill(HIST("hPtRecSig"), candidate.pt());    // rec. level pT
        registry.fill(HIST("hEtaRecSig"), candidate.eta());
        registry.fill(HIST("hYRecSig"), YXicc(candidate));
        registry.fill(HIST("hmassSig"), InvMassXiccToXicPi(candidate), candidate.pt()); //FIXME need to consider the two mass hp
        registry.fill(HIST("hDecLengthSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hChi2PCASig"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCPASig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hd0Prong0Sig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1Sig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hCtSig"), CtXicc(candidate), candidate.pt());
        registry.fill(HIST("hCPASig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("habsCPASig"), std::abs(candidate.cpa()), candidate.pt());
        registry.fill(HIST("hEtaSig"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hYSig"), YXicc(candidate), candidate.pt());
        registry.fill(HIST("hImpParErr0Sig"), candidate.errorImpactParameter0(), candidate.pt());
        registry.fill(HIST("hImpParErr1Sig"), candidate.errorImpactParameter1(), candidate.pt());
        registry.fill(HIST("hXSecVtxPosDiff"), candidate.xSecondaryVertex() - particleXic.vx(), candidate.pt());
        registry.fill(HIST("hYSecVtxPosDiff"), candidate.ySecondaryVertex() - particleXic.vy(), candidate.pt());
        registry.fill(HIST("hZSecVtxPosDiff"), candidate.zSecondaryVertex() - particleXic.vz(), candidate.pt());
        registry.fill(HIST("hPtRecGenDiff"), candidate.pt() - particleXicc.pt(), candidate.pt());
        // Check Y dependence (To be removed)
        registry.fill(HIST("hmassSigvsPtvsY"), InvMassXiccToXicPi(candidate), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hDecLengthSigvsPtvsY"), candidate.decayLength(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hChi2PCASigvsPtvsY"), candidate.chi2PCA(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hCPASigvsPtvsY"), candidate.cpa(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hd0Prong0SigvsPtvsY"), candidate.impactParameter0(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hd0Prong1SigvsPtvsY"), candidate.impactParameter1(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hCtSigvsPtvsY"), CtXicc(candidate), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hCPASigvsPtvsY"), candidate.cpa(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("habsCPASigvsPtvsY"), std::abs(candidate.cpa()), candidate.pt(), YXicc(candidate));
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta());
        registry.fill(HIST("hYRecBg"), YXicc(candidate));
        registry.fill(HIST("hmassBg"), InvMassXiccToXicPi(candidate), candidate.pt()); //FIXME need to consider the two mass hp
        registry.fill(HIST("hDecLengthBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hChi2PCABg"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCPABg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("habsCPABg"), std::abs(candidate.cpa()), candidate.pt());
        registry.fill(HIST("hd0Prong0Bg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1Bg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hCtBg"), CtXicc(candidate), candidate.pt());
        registry.fill(HIST("hCPABg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaBg"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hYBg"), YXicc(candidate), candidate.pt());
        registry.fill(HIST("hImpParErr0Bg"), candidate.errorImpactParameter0(), candidate.pt());
        registry.fill(HIST("hImpParErr1Bg"), candidate.errorImpactParameter1(), candidate.pt());
        registry.fill(HIST("hDebugMCmatching"), candidate.debugMCRec(), candidate.pt());
        // Check Y dependence (To be removed)
        registry.fill(HIST("hmassBgvsPtvsY"), InvMassXiccToXicPi(candidate), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hDecLengthBgvsPtvsY"), candidate.decayLength(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hChi2PCABgvsPtvsY"), candidate.chi2PCA(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hCPABgvsPtvsY"), candidate.cpa(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("habsCPABgvsPtvsY"), std::abs(candidate.cpa()), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hd0Prong0BgvsPtvsY"), candidate.impactParameter0(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hd0Prong1BgvsPtvsY"), candidate.impactParameter1(), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hCtBgvsPtvsY"), CtXicc(candidate), candidate.pt(), YXicc(candidate));
        registry.fill(HIST("hCPABgvsPtvsY"), candidate.cpa(), candidate.pt(), YXicc(candidate));
      }
    } // end of loop over reconstructed candidates
    // MC gen.
    //Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::XiccToXicPi) {
        if (cutYCandMax >= 0. && std::abs(RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > cutYCandMax) {
          continue;
        }
        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta());
        registry.fill(HIST("hYGen"), RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode())));
        registry.fill(HIST("hPtvsEtavsYGen"), particle.pt(), particle.eta(), RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode())));
      }
    } // end of loop of MC particles
  }   // end of process function
};    // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfTaskXicc>(cfgc)};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfTaskXiccMc>(cfgc));
  }
  return workflow;
}
