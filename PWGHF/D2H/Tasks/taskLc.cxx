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
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, University and INFN SALERNO
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Λc± → p± K∓ π± analysis task
struct HfTaskLc {
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

  HistogramRegistry registry{
    "registry",
    {/// mass candidate
     {"Data/hMass", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/signal/hMassRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/prompt/hMassRecSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/nonprompt/hMassRecSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     /// pT
     {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecSigPrompt", "3-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/signal/hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/prompt/hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/nonprompt/hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/signal/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng0Sig", "3-prong candidates (matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng1Sig", "3-prong candidates (matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng2Sig", "3-prong candidates (matched);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     /// DCAxy to prim. vertex prongs
     {"Data/hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng0Sig", "3-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng1Sig", "3-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng2Sig", "3-prong candidates (matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     /// decay length candidate
     {"Data/hDecLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/signal/hDecLengthRecSig", "3-prong candidates (matched);decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/prompt/hDecLengthRecSigPrompt", "3-prong candidates (matched, prompt);decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/nonprompt/hDecLengthRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     /// decay length xy candidate
     {"Data/hDecLengthxy", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/signal/hDecLengthxyRecSig", "3-prong candidates (matched);decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/prompt/hDecLengthxyRecSigPrompt", "3-prong candidates (matched, prompt);decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/nonprompt/hDecLengthxyRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     /// proper lifetime
     {"Data/hCt", "3-prong candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/signal/hCtRecSig", "3-prong candidates (matched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/prompt/hCtRecSigPrompt", "3-prong candidates (matched, prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/nonprompt/hCtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     /// cosine of pointing angle
     {"Data/hCPA", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/signal/hCPARecSig", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/prompt/hCPARecSigPrompt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/nonprompt/hCPARecSigNonPrompt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     /// cosine of pointing angle xy
     {"Data/hCPAxy", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/signal/hCPAxyRecSig", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/prompt/hCPAxyRecSigPrompt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/nonprompt/hCPAxyRecSigNonPrompt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     /// Chi 2 PCA to sec. vertex
     {"Data/hDca2", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/signal/hDca2RecSig", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/prompt/hDca2RecSigPrompt", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/nonprompt/hDca2RecSigNonPrompt", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     /// eta
     {"Data/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/signal/hEtaRecSig", "3-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/prompt/hEtaRecSigPrompt", "3-prong candidates (matched, prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/nonprompt/hEtaRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/signal/hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/prompt/hEtaGenPrompt", "MC particles (matched, prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/nonprompt/hEtaGenNonPrompt", "MC particles (matched, non-prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/signal/hYGen", "MC particles (matched);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/prompt/hYGenPrompt", "MC particles (matched, prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/nonprompt/hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     /// phi
     {"Data/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/signal/hPhiRecSig", "3-prong candidates (matched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/prompt/hPhiRecSigPrompt", "3-prong candidates (matched, prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/nonprompt/hPhiRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/signal/hPhiGen", "MC particles (matched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/prompt/hPhiGenPrompt", "MC particles (matched, prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/nonprompt/hPhiGenNonPrompt", "MC particles (matched, non-prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    /// mass candidate
    registry.add("Data/hMassVsPtVsMult", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {5000, 0., 10000.}}});
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hMassVsPtRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hMassVsPtRecSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// DCAxy to prim. vertex prongs
    registry.add("Data/hd0VsPtProng0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng0Sig", "3-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng1Sig", "3-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng2Sig", "3-prong candidates (matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// decay length candidate
    registry.add("Data/hDecLengthVsPt", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthVsPtRecSig", "3-prong candidates (matched);decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLengthVsPtRecSigPrompt", "3-prong candidates (matched, prompt);decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLengthVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// decay length xy candidate
    registry.add("Data/hDecLengthxyVsPt", "3-prong candidates;decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthxyVsPtRecSig", "3-prong candidates (matched);decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLengthxyVsPtRecSigPrompt", "3-prong candidates (matched, prompt);decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLengthxyVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// proper lifetime
    registry.add("Data/hCtVsPt", "3-prong candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCtVsPtRecSig", "3-prong candidates (matched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCtVsPtRecSigPrompt", "3-prong candidates (matched, prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCtVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// cosine of pointing angle
    registry.add("Data/hCPAVsPt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAVsPtRecSig", "3-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCPAVsPtRecSigPrompt", "3-prong candidates (matched, prompt);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCPAVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// cosine of pointing angle xy
    registry.add("Data/hCPAxyVsPt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAxyVsPtRecSig", "3-prong candidates (matched);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCPAxyVsPtRecSigPrompt", "3-prong candidates (matched, prompt);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCPAxyVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// Chi 2 PCA to sec. vertex
    registry.add("Data/hDca2VsPt", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDca2VsPtRecSig", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDca2VsPtRecSigPrompt", "3-prong candidates (matched, prompt);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDca2VsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// eta
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hEtaVsPtRecSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hEtaVsPtRecSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hEtaVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hEtaVsPtGenSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hEtaVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hEtaVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// y
    registry.add("MC/generated/signal/hYVsPtGenSig", "3-prong candidates (matched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hYVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hYVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// phi
    registry.add("Data/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPhiVsPtRecSig", "3-prong candidates (matched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hPhiVsPtRecSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hPhiVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hPhiVsPtGenSig", "3-prong candidates (matched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hPhiVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hPhiVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// selection status
    registry.add("hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// impact parameter error
    registry.add("Data/hImpParErrProng0", "3-prong candidates;prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng1", "3-prong candidates;prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng2", "3-prong candidates;prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hImpParErrProng0Sig", "3-prong candidates (matched);prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hImpParErrProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hImpParErrProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hImpParErrProng1Sig", "3-prong candidates (matched);prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hImpParErrProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hImpParErrProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hImpParErrProng2Sig", "3-prong candidates (matched);prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hImpParErrProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hImpParErrProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// decay length error
    registry.add("Data/hDecLenErr", "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLenErrSig", "3-prong candidates (matched);decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLenErrSigPrompt", "3-prong candidates (matched, prompt);decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLenErrSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(aod::Collision const& collision,
               aod::TracksWDca const& tracks,
               soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>> const& candidates)
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

    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLc(candidate)) > yCandRecoMax) {
        continue;
      }
      auto pt = candidate.pt();
      if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
        registry.fill(HIST("Data/hMass"), hfHelper.invMassLcToPKPi(candidate));
        registry.fill(HIST("Data/hMassVsPtVsMult"), hfHelper.invMassLcToPKPi(candidate), pt, nTracks);
        registry.fill(HIST("Data/hMassVsPt"), hfHelper.invMassLcToPKPi(candidate), pt);
      }
      if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
        registry.fill(HIST("Data/hMass"), hfHelper.invMassLcToPiKP(candidate));
        registry.fill(HIST("Data/hMassVsPtVsMult"), hfHelper.invMassLcToPiKP(candidate), pt, nTracks);
        registry.fill(HIST("Data/hMassVsPt"), hfHelper.invMassLcToPiKP(candidate), pt);
      }
      registry.fill(HIST("Data/hPt"), pt);
      registry.fill(HIST("Data/hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("Data/hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("Data/hPtProng2"), candidate.ptProng2());
      registry.fill(HIST("Data/hd0Prong0"), candidate.impactParameter0());
      registry.fill(HIST("Data/hd0Prong1"), candidate.impactParameter1());
      registry.fill(HIST("Data/hd0Prong2"), candidate.impactParameter2());
      registry.fill(HIST("Data/hd0VsPtProng0"), candidate.impactParameter0(), pt);
      registry.fill(HIST("Data/hd0VsPtProng1"), candidate.impactParameter1(), pt);
      registry.fill(HIST("Data/hd0VsPtProng2"), candidate.impactParameter2(), pt);
      registry.fill(HIST("Data/hDecLength"), candidate.decayLength());
      registry.fill(HIST("Data/hDecLengthVsPt"), candidate.decayLength(), pt);
      registry.fill(HIST("Data/hDecLengthxy"), candidate.decayLengthXY());
      registry.fill(HIST("Data/hDecLengthxyVsPt"), candidate.decayLengthXY(), pt);
      registry.fill(HIST("Data/hCt"), hfHelper.ctLc(candidate));
      registry.fill(HIST("Data/hCtVsPt"), hfHelper.ctLc(candidate), pt);
      registry.fill(HIST("Data/hCPA"), candidate.cpa());
      registry.fill(HIST("Data/hCPAVsPt"), candidate.cpa(), pt);
      registry.fill(HIST("Data/hCPAxy"), candidate.cpaXY());
      registry.fill(HIST("Data/hCPAxyVsPt"), candidate.cpaXY(), pt);
      registry.fill(HIST("Data/hDca2"), candidate.chi2PCA());
      registry.fill(HIST("Data/hDca2VsPt"), candidate.chi2PCA(), pt);
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), pt);
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPhiVsPt"), candidate.phi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelLcToPKPi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelLcToPiKP(), pt);
      registry.fill(HIST("Data/hImpParErrProng0"), candidate.errorImpactParameter0(), pt);
      registry.fill(HIST("Data/hImpParErrProng1"), candidate.errorImpactParameter1(), pt);
      registry.fill(HIST("Data/hImpParErrProng2"), candidate.errorImpactParameter2(), pt);
      registry.fill(HIST("Data/hDecLenErr"), candidate.errorDecayLength(), pt);
    }
  }

  /// Fills MC histograms.
  void processMc(soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles,
                 aod::TracksWMc const&)
  {
    for (const auto& candidate : candidates) {
      /// Select Lc
      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      /// rapidity selection
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLc(candidate)) > yCandRecoMax) {
        continue;
      }

      if (std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        // Get the corresponding MC particle.
        auto mcParticleProng0 = candidate.prong0_as<aod::TracksWMc>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>();
        auto pdgCodeProng0 = std::abs(mcParticleProng0.pdgCode());
        auto indexMother = RecoDecay::getMother(mcParticles, mcParticleProng0, pdg::Code::kLambdaCPlus, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);
        registry.fill(HIST("MC/generated/signal/hPtGenSig"), particleMother.pt()); // gen. level pT
        auto pt = candidate.pt();
        /// MC reconstructed signal
        if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
          registry.fill(HIST("MC/reconstructed/signal/hMassRecSig"), hfHelper.invMassLcToPKPi(candidate));
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecSig"), hfHelper.invMassLcToPKPi(candidate), pt);
        }
        if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
          registry.fill(HIST("MC/reconstructed/signal/hMassRecSig"), hfHelper.invMassLcToPiKP(candidate));
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecSig"), hfHelper.invMassLcToPiKP(candidate), pt);
        }
        registry.fill(HIST("MC/reconstructed/signal/hPtRecSig"), pt);
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng0Sig"), candidate.ptProng0());
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng1Sig"), candidate.ptProng1());
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng2Sig"), candidate.ptProng2());

        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng0Sig"), candidate.impactParameter0());
        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng1Sig"), candidate.impactParameter1());
        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng2Sig"), candidate.impactParameter2());
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng0Sig"), candidate.impactParameter0(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng1Sig"), candidate.impactParameter1(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng2Sig"), candidate.impactParameter2(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthRecSig"), candidate.decayLength());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthVsPtRecSig"), candidate.decayLength(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthxyRecSig"), candidate.decayLengthXY());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthxyVsPtRecSig"), candidate.decayLengthXY(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hCtRecSig"), hfHelper.ctLc(candidate));
        registry.fill(HIST("MC/reconstructed/signal/hCtVsPtRecSig"), hfHelper.ctLc(candidate), pt);
        registry.fill(HIST("MC/reconstructed/signal/hCPARecSig"), candidate.cpa());
        registry.fill(HIST("MC/reconstructed/signal/hCPAVsPtRecSig"), candidate.cpa(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hCPAxyRecSig"), candidate.cpaXY());
        registry.fill(HIST("MC/reconstructed/signal/hCPAxyVsPtRecSig"), candidate.cpaXY(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hDca2RecSig"), candidate.chi2PCA());
        registry.fill(HIST("MC/reconstructed/signal/hDca2VsPtRecSig"), candidate.chi2PCA(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecSig"), candidate.eta());
        registry.fill(HIST("MC/reconstructed/signal/hEtaVsPtRecSig"), candidate.eta(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hPhiRecSig"), candidate.phi());
        registry.fill(HIST("MC/reconstructed/signal/hPhiVsPtRecSig"), candidate.phi(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hImpParErrProng0Sig"), candidate.errorImpactParameter0(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hImpParErrProng1Sig"), candidate.errorImpactParameter1(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hImpParErrProng2Sig"), candidate.errorImpactParameter2(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hDecLenErrSig"), candidate.errorDecayLength(), pt);

        /// reconstructed signal prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecSigPrompt"), hfHelper.invMassLcToPKPi(candidate));
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecSigPrompt"), hfHelper.invMassLcToPKPi(candidate), pt);
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecSigPrompt"), hfHelper.invMassLcToPiKP(candidate));
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecSigPrompt"), hfHelper.invMassLcToPiKP(candidate), pt);
          }
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecSigPrompt"), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng0SigPrompt"), candidate.ptProng0());
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng1SigPrompt"), candidate.ptProng1());
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng2SigPrompt"), candidate.ptProng2());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng0SigPrompt"), candidate.impactParameter0());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng1SigPrompt"), candidate.impactParameter1());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng2SigPrompt"), candidate.impactParameter2());
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng0SigPrompt"), candidate.impactParameter0(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng1SigPrompt"), candidate.impactParameter1(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng2SigPrompt"), candidate.impactParameter2(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthRecSigPrompt"), candidate.decayLength());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthVsPtRecSigPrompt"), candidate.decayLength(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthxyRecSigPrompt"), candidate.decayLengthXY());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthxyVsPtRecSigPrompt"), candidate.decayLengthXY(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hCtRecSigPrompt"), hfHelper.ctLc(candidate));
          registry.fill(HIST("MC/reconstructed/prompt/hCtVsPtRecSigPrompt"), hfHelper.ctLc(candidate), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hCPARecSigPrompt"), candidate.cpa());
          registry.fill(HIST("MC/reconstructed/prompt/hCPAVsPtRecSigPrompt"), candidate.cpa(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hCPAxyRecSigPrompt"), candidate.cpaXY());
          registry.fill(HIST("MC/reconstructed/prompt/hCPAxyVsPtRecSigPrompt"), candidate.cpaXY(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hDca2RecSigPrompt"), candidate.chi2PCA());
          registry.fill(HIST("MC/reconstructed/prompt/hDca2VsPtRecSigPrompt"), candidate.chi2PCA(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hEtaRecSigPrompt"), candidate.eta());
          registry.fill(HIST("MC/reconstructed/prompt/hEtaVsPtRecSigPrompt"), candidate.eta(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hPhiRecSigPrompt"), candidate.phi());
          registry.fill(HIST("MC/reconstructed/prompt/hPhiVsPtRecSigPrompt"), candidate.phi(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hImpParErrProng0SigPrompt"), candidate.errorImpactParameter0(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hImpParErrProng1SigPrompt"), candidate.errorImpactParameter1(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hImpParErrProng2SigPrompt"), candidate.errorImpactParameter2(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hDecLenErrSigPrompt"), candidate.errorDecayLength(), pt);
        } else {
          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecSigNonPrompt"), hfHelper.invMassLcToPKPi(candidate));
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt"), hfHelper.invMassLcToPKPi(candidate), pt);
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecSigNonPrompt"), hfHelper.invMassLcToPiKP(candidate));
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt"), hfHelper.invMassLcToPiKP(candidate), pt);
          }
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecSigNonPrompt"), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng0SigNonPrompt"), candidate.ptProng0());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng1SigNonPrompt"), candidate.ptProng1());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng2SigNonPrompt"), candidate.ptProng2());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng0SigNonPrompt"), candidate.impactParameter0());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng1SigNonPrompt"), candidate.impactParameter1());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng2SigNonPrompt"), candidate.impactParameter2());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng0SigNonPrompt"), candidate.impactParameter0(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng1SigNonPrompt"), candidate.impactParameter1(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng2SigNonPrompt"), candidate.impactParameter2(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthRecSigNonPrompt"), candidate.decayLength());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthVsPtRecSigNonPrompt"), candidate.decayLength(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthxyRecSigNonPrompt"), candidate.decayLengthXY());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthxyVsPtRecSigNonPrompt"), candidate.decayLengthXY(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hCtRecSigNonPrompt"), hfHelper.ctLc(candidate));
          registry.fill(HIST("MC/reconstructed/nonprompt/hCtVsPtRecSigNonPrompt"), hfHelper.ctLc(candidate), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPARecSigNonPrompt"), candidate.cpa());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAVsPtRecSigNonPrompt"), candidate.cpa(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAxyRecSigNonPrompt"), candidate.cpaXY());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAxyVsPtRecSigNonPrompt"), candidate.cpaXY(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDca2RecSigNonPrompt"), candidate.chi2PCA());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDca2VsPtRecSigNonPrompt"), candidate.chi2PCA(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaRecSigNonPrompt"), candidate.eta());
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaVsPtRecSigNonPrompt"), candidate.eta(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hPhiRecSigNonPrompt"), candidate.phi());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPhiVsPtRecSigNonPrompt"), candidate.phi(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hImpParErrProng0SigNonPrompt"), candidate.errorImpactParameter0(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hImpParErrProng1SigNonPrompt"), candidate.errorImpactParameter1(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hImpParErrProng2SigNonPrompt"), candidate.errorImpactParameter2(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLenErrSigNonPrompt"), candidate.errorDecayLength(), pt);
        }
      }
    }
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        auto yGen = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::analysis::pdg::MassLambdaCPlus);
        if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
          continue;
        }
        auto ptGen = particle.pt();
        registry.fill(HIST("MC/generated/signal/hPtGen"), ptGen);
        registry.fill(HIST("MC/generated/signal/hEtaGen"), particle.eta());
        registry.fill(HIST("MC/generated/signal/hYGen"), yGen);
        registry.fill(HIST("MC/generated/signal/hPhiGen"), particle.phi());
        registry.fill(HIST("MC/generated/signal/hEtaVsPtGenSig"), particle.eta(), ptGen);
        registry.fill(HIST("MC/generated/signal/hYVsPtGenSig"), yGen, ptGen);
        registry.fill(HIST("MC/generated/signal/hPhiVsPtGenSig"), particle.phi(), ptGen);

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/generated/prompt/hPtGenPrompt"), ptGen);
          registry.fill(HIST("MC/generated/prompt/hEtaGenPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/prompt/hYGenPrompt"), yGen);
          registry.fill(HIST("MC/generated/prompt/hPhiGenPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/prompt/hEtaVsPtGenSigPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/prompt/hYVsPtGenSigPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/prompt/hPhiVsPtGenSigPrompt"), particle.phi(), ptGen);
        }
        if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
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

  PROCESS_SWITCH(HfTaskLc, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskLc>(cfgc)};
}
