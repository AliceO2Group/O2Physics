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

/// \file taskLc.cxx
/// \brief Λc± → p± K∓ π± analysis task
/// \note Extended from taskD0
///
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, University and INFN SALERNO
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::analysis::hf_cuts_lc_topkpi;

#include "Framework/runDataProcessing.h"

/// Λc± → p± K∓ π± analysis task
struct TaskLc {
  HistogramRegistry registry{
    "registry",
    {{"Data/hMass", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/signal/hMassRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/prompt/hMassRecSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/nonprompt/hMassRecSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/background/hMassRecBkg", "3-prong candidates (unmatched);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecSigPrompt", "3-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/background/hPtRecBkg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/signal/hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/prompt/hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/nonprompt/hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/signal/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng0Sig", "3-prong candidates (matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/background/hPtRecProng0Bkg", "3-prong candidates (unmatched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng1Sig", "3-prong candidates (matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/background/hPtRecProng1Bkg", "3-prong candidates (unmatched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng2Sig", "3-prong candidates (matched);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/background/hPtRecProng2Bkg", "3-prong candidates (unmatched);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"Data/hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng0Sig", "3-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/background/hd0RecProng0Bkg", "3-prong candidates (unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng1Sig", "3-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/background/hd0RecProng1Bkg", "3-prong candidates (unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng2Sig", "3-prong candidates (matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/background/hd0RecProng2Bkg", "3-prong candidates (unmatched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hDecLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/signal/hDecLengthRecSig", "3-prong candidates (matched);decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/prompt/hDecLengthRecSigPrompt", "3-prong candidates (matched, prompt);decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/nonprompt/hDecLengthRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/background/hDecLengthRecBkg", "3-prong candidates (unmatched);decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"Data/hDecLengthxy", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/signal/hDecLengthxyRecSig", "3-prong candidates (matched);decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/prompt/hDecLengthxyRecSigPrompt", "3-prong candidates (matched, prompt);decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/nonprompt/hDecLengthxyRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/background/hDecLengthxyRecBkg", "3-prong candidates (unmatched);decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"Data/hCt", "3-prong candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/signal/hCtRecSig", "3-prong candidates (matched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/prompt/hCtRecSigPrompt", "3-prong candidates (matched, prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/nonprompt/hCtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/background/hCtRecBkg", "3-prong candidates (unmatched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"Data/hCPA", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/signal/hCPARecSig", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/prompt/hCPARecSigPrompt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/nonprompt/hCPARecSigNonPrompt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/background/hCPARecBkg", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"Data/hCPAxy", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/signal/hCPAxyRecSig", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/prompt/hCPAxyRecSigPrompt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/nonprompt/hCPAxyRecSigNonPrompt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/background/hCPAxyRecBkg", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"Data/hDca2", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/signal/hDca2RecSig", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/prompt/hDca2RecSigPrompt", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/nonprompt/hDca2RecSigNonPrompt", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/background/hDca2RecBkg", "3-prong candidates (unmatched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"Data/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/signal/hEtaRecSig", "3-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/prompt/hEtaRecSigPrompt", "3-prong candidates (matched, prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/nonprompt/hEtaRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/background/hEtaRecBkg", "3-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/signal/hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/prompt/hEtaGenPrompt", "MC particles (matched, prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/nonprompt/hEtaGenNonPrompt", "MC particles (matched, non-prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/signal/hYGen", "MC particles (matched);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/prompt/hYGenPrompt", "MC particles (matched, prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/nonprompt/hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"Data/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/signal/hPhiRecSig", "3-prong candidates (matched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/prompt/hPhiRecSigPrompt", "3-prong candidates (matched, prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/nonprompt/hPhiRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/background/hPhiRecBkg", "3-prong candidates (unmatched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/signal/hPhiGen", "MC particles (matched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/prompt/hPhiGenPrompt", "MC particles (matched, prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/nonprompt/hPhiGenNonPrompt", "MC particles (matched, non-prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}}}};

  Configurable<int> d_selectionFlagLc{"d_selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_lc_topkpi::pTBins_v}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_selcandidate_lc::isSelLcpKpi >= d_selectionFlagLc || aod::hf_selcandidate_lc::isSelLcpiKp >= d_selectionFlagLc);

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    /// mass candidate
    registry.add("Data/hMassVsPtVsMult", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {5000, 0., 10000.}}});
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hMassVsPtRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hMassVsPtRecSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hMassVsPtRecBkg", "3-prong candidates (unmatched);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// DCAxy to prim. vertex prongs
    registry.add("Data/hd0VsPtProng0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng0Sig", "3-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0VsPtRecProng0Bkg", "3-prong candidates (unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng1Sig", "3-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0VsPtRecProng1Bkg", "3-prong candidates (unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng2Sig", "3-prong candidates (matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0VsPtRecProng2Bkg", "3-prong candidates (unmatched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// decay length candidate
    registry.add("Data/hDecLengthVsPt", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthVsPtRecSig", "3-prong candidates (matched);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLengthVsPtRecSigPrompt", "3-prong candidates (matched, prompt);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLengthVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLengthVsPtRecBkg", "3-prong candidates (unmatched);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// decay length xy candidate
    registry.add("Data/hDecLengthxyVsPt", "3-prong candidates;decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthxyVsPtRecSig", "3-prong candidates (matched);decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLengthxyVsPtRecSigPrompt", "3-prong candidates (matched, prompt);decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLengthxyVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLengthxyVsPtRecBkg", "3-prong candidates (unmatched);decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// proper lifetime
    registry.add("Data/hCtVsPt", "3-prong candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCtVsPtRecSig", "3-prong candidates (matched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCtVsPtRecSigPrompt", "3-prong candidates (matched, prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCtVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCtVsPtRecBkg", "3-prong candidates (unmatched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// cosine of pointing angle
    registry.add("Data/hCPAVsPt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAVsPtRecSig", "3-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCPAVsPtRecSigPrompt", "3-prong candidates (matched, prompt);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCPAVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCPAVsPtRecBkg", "3-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// cosine of pointing angle xy
    registry.add("Data/hCPAxyVsPt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAxyVsPtRecSig", "3-prong candidates (matched);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCPAxyVsPtRecSigPrompt", "3-prong candidates (matched, prompt);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCPAxyVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCPAxyVsPtRecBkg", "3-prong candidates (unmatched);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// Chi 2 PCA to sec. vertex
    registry.add("Data/hDca2VsPt", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDca2VsPtRecSig", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDca2VsPtRecSigPrompt", "3-prong candidates (matched, prompt);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDca2VsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDca2VsPtRecBkg", "3-prong candidates (unmatched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// eta
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hEtaVsPtRecSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hEtaVsPtRecSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hEtaVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hEtaVsPtRecBkg", "3-prong candidates (unmatched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hEtaVsPtGenSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hEtaVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hEtaVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// y
    // registry.add("Data/hYVsPt", "3-prong candidates;candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // registry.add("MC/reconstructed/signal/hYVsPtRecSig", "3-prong candidates (matched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // registry.add("MC/reconstructed/prompt/hYVsPtRecSigPrompt", "3-prong candidates (matched, prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // registry.add("MC/reconstructed/nonprompt/hYVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // registry.add("MC/reconstructed/background/hYVsPtRecBkg", "3-prong candidates (unmatched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hYVsPtGenSig", "3-prong candidates (matched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hYVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hYVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// phi
    registry.add("Data/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPhiVsPtRecSig", "3-prong candidates (matched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hPhiVsPtRecSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hPhiVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hPhiVsPtRecBkg", "3-prong candidates (unmatched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hPhiVsPtGenSig", "3-prong candidates (matched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hPhiVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hPhiVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// selection status
    registry.add("hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // /// impact parameter error
    // registry.add("Data/hImpParErr", "3-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // /// decay length error
    // registry.add("Data/hDecLenErr", "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  // FIXME: Add ALICE 2/3 switch!
  // void process(aod::HfCandProng3 const& candidates)
  void process(const o2::aod::Collision& collision, const soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>> const& candidates)
  {
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > 4.0) {
          continue;
        }
        if (std::abs(track.dcaXY()) > 0.0025 || std::abs(track.dcaZ()) > 0.0025) {
          continue;
        }
        nTracks++;
      }
    }
    registry.fill(HIST("Data/hMultiplicity"), nTracks);

    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YLc(candidate)) > cutYCandMax) {
        continue;
      }
      if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
        registry.fill(HIST("Data/hMass"), InvMassLcpKpi(candidate));
        registry.fill(HIST("Data/hMassVsPtVsMult"), InvMassLcpKpi(candidate), candidate.pt(), nTracks);
        registry.fill(HIST("Data/hMassVsPt"), InvMassLcpKpi(candidate), candidate.pt());
      }
      if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
        registry.fill(HIST("Data/hMass"), InvMassLcpiKp(candidate));
        registry.fill(HIST("Data/hMassVsPtVsMult"), InvMassLcpiKp(candidate), candidate.pt(), nTracks);
        registry.fill(HIST("Data/hMassVsPt"), InvMassLcpiKp(candidate), candidate.pt());
      }
      registry.fill(HIST("Data/hPt"), candidate.pt());
      registry.fill(HIST("Data/hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("Data/hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("Data/hPtProng2"), candidate.ptProng2());
      registry.fill(HIST("Data/hd0Prong0"), candidate.impactParameter0());
      registry.fill(HIST("Data/hd0Prong1"), candidate.impactParameter1());
      registry.fill(HIST("Data/hd0Prong2"), candidate.impactParameter2());
      registry.fill(HIST("Data/hd0VsPtProng0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("Data/hd0VsPtProng1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("Data/hd0VsPtProng2"), candidate.impactParameter2(), candidate.pt());
      registry.fill(HIST("Data/hDecLength"), candidate.decayLength());
      registry.fill(HIST("Data/hDecLengthVsPt"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("Data/hDecLengthxy"), candidate.decayLengthXY());
      registry.fill(HIST("Data/hDecLengthxyVsPt"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("Data/hCt"), CtLc(candidate));
      registry.fill(HIST("Data/hCtVsPt"), CtLc(candidate), candidate.pt());
      registry.fill(HIST("Data/hCPA"), candidate.cpa());
      registry.fill(HIST("Data/hCPAVsPt"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("Data/hCPAxy"), candidate.cpaXY());
      registry.fill(HIST("Data/hCPAxyVsPt"), candidate.cpaXY(), candidate.pt());
      registry.fill(HIST("Data/hDca2"), candidate.chi2PCA());
      registry.fill(HIST("Data/hDca2VsPt"), candidate.chi2PCA(), candidate.pt());
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), candidate.pt());
      // registry.fill(HIST("Data/hY"), candidate.y());
      // registry.fill(HIST("Data/hYVsPt"), candidate.y(), candidate.pt());
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPhiVsPt"), candidate.phi(), candidate.pt());
      registry.fill(HIST("hSelectionStatus"), candidate.isSelLcpKpi(), candidate.pt());
      registry.fill(HIST("hSelectionStatus"), candidate.isSelLcpiKp(), candidate.pt());
      // registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      // registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      // registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter2(), candidate.pt());
      // registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
    }
  }

  /// Fills MC histograms.
  void processMC(soa::Filtered<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate, aod::HfCandProng3MCRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandProng3MCGen> const& particlesMC, aod::BigTracksMC const& /*tracks*/)
  {
    for (auto& candidate : candidates) {
      /// Select Lc
      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        continue;
      }
      /// rapidity selection
      if (cutYCandMax >= 0. && std::abs(YLc(candidate)) > cutYCandMax) {
        continue;
      }

      if (std::abs(candidate.flagMCMatchRec()) == 1 << DecayType::LcToPKPi) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.index0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandProng3MCGen>>(), pdg::Code::kLambdaCPlus, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);
        registry.fill(HIST("MC/generated/signal/hPtGenSig"), particleMother.pt()); // gen. level pT

        /// MC reconstructed signal
        if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
          registry.fill(HIST("MC/reconstructed/signal/hMassRecSig"), InvMassLcpKpi(candidate));
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecSig"), InvMassLcpKpi(candidate), candidate.pt());
        }
        if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
          registry.fill(HIST("MC/reconstructed/signal/hMassRecSig"), InvMassLcpiKp(candidate));
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecSig"), InvMassLcpiKp(candidate), candidate.pt());
        }
        registry.fill(HIST("MC/reconstructed/signal/hPtRecSig"), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng0Sig"), candidate.ptProng0());
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng1Sig"), candidate.ptProng1());
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng2Sig"), candidate.ptProng2());

        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng0Sig"), candidate.impactParameter0());
        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng1Sig"), candidate.impactParameter1());
        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng2Sig"), candidate.impactParameter2());
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng0Sig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng1Sig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng2Sig"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthRecSig"), candidate.decayLength());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthVsPtRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthxyRecSig"), candidate.decayLengthXY());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthxyVsPtRecSig"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hCtRecSig"), CtLc(candidate));
        registry.fill(HIST("MC/reconstructed/signal/hCtVsPtRecSig"), CtLc(candidate), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hCPARecSig"), candidate.cpa());
        registry.fill(HIST("MC/reconstructed/signal/hCPAVsPtRecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hCPAxyRecSig"), candidate.cpaXY());
        registry.fill(HIST("MC/reconstructed/signal/hCPAxyVsPtRecSig"), candidate.cpaXY(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hDca2RecSig"), candidate.chi2PCA());
        registry.fill(HIST("MC/reconstructed/signal/hDca2VsPtRecSig"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecSig"), candidate.eta());
        registry.fill(HIST("MC/reconstructed/signal/hEtaVsPtRecSig"), candidate.eta(), candidate.pt());
        // registry.fill(HIST("MC/reconstructed/signal/hYRecSig"), candidate.y());
        // registry.fill(HIST("MC/reconstructed/signal/hYVsPtRecSig"), candidate.y(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hPhiRecSig"), candidate.phi());
        registry.fill(HIST("MC/reconstructed/signal/hPhiVsPtRecSig"), candidate.phi(), candidate.pt());

        /// reconstructed signal prompt
        if (candidate.originMCRec() == RecoDecay::OriginType::Prompt) {
          if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecSigPrompt"), InvMassLcpKpi(candidate));
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecSigPrompt"), InvMassLcpKpi(candidate), candidate.pt());
          }
          if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecSigPrompt"), InvMassLcpiKp(candidate));
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecSigPrompt"), InvMassLcpiKp(candidate), candidate.pt());
          }
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecSigPrompt"), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng0SigPrompt"), candidate.ptProng0());
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng1SigPrompt"), candidate.ptProng1());
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng2SigPrompt"), candidate.ptProng2());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng0SigPrompt"), candidate.impactParameter0());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng1SigPrompt"), candidate.impactParameter1());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng2SigPrompt"), candidate.impactParameter2());
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng0SigPrompt"), candidate.impactParameter0(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng1SigPrompt"), candidate.impactParameter1(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng2SigPrompt"), candidate.impactParameter2(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthRecSigPrompt"), candidate.decayLength());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthVsPtRecSigPrompt"), candidate.decayLength(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthxyRecSigPrompt"), candidate.decayLengthXY());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthxyVsPtRecSigPrompt"), candidate.decayLengthXY(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hCtRecSigPrompt"), CtLc(candidate));
          registry.fill(HIST("MC/reconstructed/prompt/hCtVsPtRecSigPrompt"), CtLc(candidate), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hCPARecSigPrompt"), candidate.cpa());
          registry.fill(HIST("MC/reconstructed/prompt/hCPAVsPtRecSigPrompt"), candidate.cpa(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hCPAxyRecSigPrompt"), candidate.cpaXY());
          registry.fill(HIST("MC/reconstructed/prompt/hCPAxyVsPtRecSigPrompt"), candidate.cpaXY(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hDca2RecSigPrompt"), candidate.chi2PCA());
          registry.fill(HIST("MC/reconstructed/prompt/hDca2VsPtRecSigPrompt"), candidate.chi2PCA(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hEtaRecSigPrompt"), candidate.eta());
          registry.fill(HIST("MC/reconstructed/prompt/hEtaVsPtRecSigPrompt"), candidate.eta(), candidate.pt());
          // registry.fill(HIST("MC/reconstructed/prompt/hYRecSigPrompt"), candidate.y());
          // registry.fill(HIST("MC/reconstructed/prompt/hYVsPtRecSigPrompt"), candidate.y(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hPhiRecSigPrompt"), candidate.phi());
          registry.fill(HIST("MC/reconstructed/prompt/hPhiVsPtRecSigPrompt"), candidate.phi(), candidate.pt());
        } else {
          if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecSigNonPrompt"), InvMassLcpKpi(candidate));
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt"), InvMassLcpKpi(candidate), candidate.pt());
          }
          if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecSigNonPrompt"), InvMassLcpiKp(candidate));
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt"), InvMassLcpiKp(candidate), candidate.pt());
          }
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecSigNonPrompt"), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng0SigNonPrompt"), candidate.ptProng0());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng1SigNonPrompt"), candidate.ptProng1());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng2SigNonPrompt"), candidate.ptProng2());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng0SigNonPrompt"), candidate.impactParameter0());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng1SigNonPrompt"), candidate.impactParameter1());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng2SigNonPrompt"), candidate.impactParameter2());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng0SigNonPrompt"), candidate.impactParameter0(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng1SigNonPrompt"), candidate.impactParameter1(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng2SigNonPrompt"), candidate.impactParameter2(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthRecSigNonPrompt"), candidate.decayLength());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthVsPtRecSigNonPrompt"), candidate.decayLength(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthxyRecSigNonPrompt"), candidate.decayLengthXY());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthxyVsPtRecSigNonPrompt"), candidate.decayLengthXY(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCtRecSigNonPrompt"), CtLc(candidate));
          registry.fill(HIST("MC/reconstructed/nonprompt/hCtVsPtRecSigNonPrompt"), CtLc(candidate), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPARecSigNonPrompt"), candidate.cpa());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAVsPtRecSigNonPrompt"), candidate.cpa(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAxyRecSigNonPrompt"), candidate.cpaXY());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAxyVsPtRecSigNonPrompt"), candidate.cpaXY(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDca2RecSigNonPrompt"), candidate.chi2PCA());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDca2VsPtRecSigNonPrompt"), candidate.chi2PCA(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaRecSigNonPrompt"), candidate.eta());
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaVsPtRecSigNonPrompt"), candidate.eta(), candidate.pt());
          // registry.fill(HIST("MC/reconstructed/nonprompt/hYRecSigNonPrompt"), candidate.y());
          // registry.fill(HIST("MC/reconstructed/nonprompt/hYVsPtRecSigNonPrompt"), candidate.y(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPhiRecSigNonPrompt"), candidate.phi());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPhiVsPtRecSigNonPrompt"), candidate.phi(), candidate.pt());
        }

      } else {
        if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
          registry.fill(HIST("MC/reconstructed/background/hMassRecBkg"), InvMassLcpKpi(candidate));
          registry.fill(HIST("MC/reconstructed/background/hMassVsPtRecBkg"), InvMassLcpKpi(candidate), candidate.pt());
        }
        if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
          registry.fill(HIST("MC/reconstructed/background/hMassRecBkg"), InvMassLcpiKp(candidate));
          registry.fill(HIST("MC/reconstructed/background/hMassVsPtRecBkg"), InvMassLcpiKp(candidate), candidate.pt());
        }
        registry.fill(HIST("MC/reconstructed/background/hPtRecBkg"), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hPtRecProng0Bkg"), candidate.ptProng0());
        registry.fill(HIST("MC/reconstructed/background/hPtRecProng1Bkg"), candidate.ptProng1());
        registry.fill(HIST("MC/reconstructed/background/hPtRecProng2Bkg"), candidate.ptProng2());
        registry.fill(HIST("MC/reconstructed/background/hd0RecProng0Bkg"), candidate.impactParameter0());
        registry.fill(HIST("MC/reconstructed/background/hd0RecProng1Bkg"), candidate.impactParameter1());
        registry.fill(HIST("MC/reconstructed/background/hd0RecProng2Bkg"), candidate.impactParameter2());
        registry.fill(HIST("MC/reconstructed/background/hd0VsPtRecProng0Bkg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hd0VsPtRecProng1Bkg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hd0VsPtRecProng2Bkg"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hDecLengthRecBkg"), candidate.decayLength());
        registry.fill(HIST("MC/reconstructed/background/hDecLengthVsPtRecBkg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hDecLengthxyRecBkg"), candidate.decayLengthXY());
        registry.fill(HIST("MC/reconstructed/background/hDecLengthxyVsPtRecBkg"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hCtRecBkg"), CtLc(candidate));
        registry.fill(HIST("MC/reconstructed/background/hCtVsPtRecBkg"), CtLc(candidate), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hCPARecBkg"), candidate.cpa());
        registry.fill(HIST("MC/reconstructed/background/hCPAVsPtRecBkg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hCPAxyRecBkg"), candidate.cpaXY());
        registry.fill(HIST("MC/reconstructed/background/hCPAxyVsPtRecBkg"), candidate.cpaXY(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hDca2RecBkg"), candidate.chi2PCA());
        registry.fill(HIST("MC/reconstructed/background/hDca2VsPtRecBkg"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hEtaRecBkg"), candidate.eta());
        registry.fill(HIST("MC/reconstructed/background/hEtaVsPtRecBkg"), candidate.eta(), candidate.pt());
        // registry.fill(HIST("MC/reconstructed/background/hYRecBkg"), RecoDecay::y(array{candidate.px(), candidate.py(), candidate.pz()}, RecoDecay::getMassPDG(candidate.pdgCode())));
        // registry.fill(HIST("MC/reconstructed/background/hYVsPtRecBkg"), RecoDecay::y(array{candidate.px(), candidate.py(), candidate.pz()}, RecoDecay::getMassPDG(candidate.pdgCode())), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hPhiRecBkg"), candidate.phi());
        registry.fill(HIST("MC/reconstructed/background/hPhiVsPtRecBkg"), candidate.phi(), candidate.pt());
      }
    }
    // MC gen.
    // Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::LcToPKPi) {
        if (cutYCandMax >= 0. && std::abs(RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > cutYCandMax) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        registry.fill(HIST("MC/generated/signal/hPtGen"), ptGen);
        registry.fill(HIST("MC/generated/signal/hEtaGen"), particle.eta());
        registry.fill(HIST("MC/generated/signal/hYGen"), yGen);
        registry.fill(HIST("MC/generated/signal/hPhiGen"), particle.phi());
        registry.fill(HIST("MC/generated/signal/hEtaVsPtGenSig"), particle.eta(), ptGen);
        registry.fill(HIST("MC/generated/signal/hYVsPtGenSig"), yGen, ptGen);
        registry.fill(HIST("MC/generated/signal/hPhiVsPtGenSig"), particle.phi(), ptGen);

        if (particle.originMCGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/generated/prompt/hPtGenPrompt"), ptGen);
          registry.fill(HIST("MC/generated/prompt/hEtaGenPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/prompt/hYGenPrompt"), yGen);
          registry.fill(HIST("MC/generated/prompt/hPhiGenPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/prompt/hEtaVsPtGenSigPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/prompt/hYVsPtGenSigPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/prompt/hPhiVsPtGenSigPrompt"), particle.phi(), ptGen);
        }
        if (particle.originMCGen() == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("MC/generated/nonprompt/hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hEtaGenNonPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/nonprompt/hYGenNonPrompt"), yGen);
          registry.fill(HIST("MC/generated/nonprompt/hPhiGenNonPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/nonprompt/hEtaVsPtGenSigNonPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hYVsPtGenSigNonPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hPhiVsPtGenSigNonPrompt"), particle.phi(), ptGen);
        }
      }
    }
  }

  PROCESS_SWITCH(TaskLc, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskLc>(cfgc, TaskName{"hf-task-lc"})};
}
