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
    {
     /// mass
     {"Data/hMass", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{400, 2.1, 2.5}}}},
     {"MC/reconstructed/signal/hMassRecCandSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{400, 2.1, 2.5}}}},
     {"MC/reconstructed/prompt/hMassRecCandSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{400, 2.1, 2.5}}}},
     {"MC/reconstructed/nonprompt/hMassRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{400, 2.1, 2.5}}}},
     {"MC/reconstructed/background/hMassRecCandBkg", "3-prong candidates (unmatched);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{400, 2.1, 2.5}}}},

     /// pT candidate
     {"Data/hPtCand", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecCandSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecCandSigPrompt", "3-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/background/hPtRecCandBkg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/signal/hPtGenCand", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/prompt/hPtGenCandPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/nonprompt/hPtGenCandNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/background/hPtGenCandSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},

     /// pT prongs
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

     /// multiplicity
     {"Data/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     
     /// DCAxy to prim. vertex prongs
     {"Data/hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/signal/hd0RecProng0Sig", "3-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/prompt/hd0RecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/background/hd0RecProng0Bkg", "3-prong candidates (unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"Data/hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/signal/hd0RecProng1Sig", "3-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/prompt/hd0RecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/background/hd0RecProng1Bkg", "3-prong candidates (unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"Data/hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/signal/hd0RecProng2Sig", "3-prong candidates (matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/prompt/hd0RecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},
     {"MC/reconstructed/background/hd0RecProng2Bkg", "3-prong candidates (unmatched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -0.2, 0.2}}}},

     /// decay length candidate
     {"Data/hDecLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/signal/hDecLengthRecCandSig", "3-prong candidates (matched);decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/prompt/hDecLengthRecCandSigPrompt", "3-prong candidates (matched, prompt);decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/nonprompt/hDecLengthRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/background/hDecLengthRecCandBkg", "3-prong candidates (unmatched);decay length (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},

     /// decay length xy candidate
     {"Data/hDecLengthxy", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/signal/hDecLengthxyRecCandSig", "3-prong candidates (matched);decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/prompt/hDecLengthxyRecCandSigPrompt", "3-prong candidates (matched, prompt);decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/nonprompt/hDecLengthxyRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}},
     {"MC/reconstructed/background/hDecLengthxyRecCandBkg", "3-prong candidates (unmatched);decay length xy (cm);entries", {HistType::kTH1F, {{200, 0., 2.}}}}, 

     /// proper lifetime
     {"Data/hCt", "3-prong candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/signal/hCtRecCandSig", "3-prong candidates (matched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/prompt/hCtRecCandSigPrompt", "3-prong candidates (matched, prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/nonprompt/hCtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/background/hCtRecCandBkg", "3-prong candidates (unmatched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},

     /// cosine of pointing angle candidate 
     {"Data/hCPA", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/signal/hCPARecCandSig", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/prompt/hCPARecCandSigPrompt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/nonprompt/hCPARecCandSigNonPrompt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/background/hCPARecCandBkg", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},

     /// cosine of pointing angle xy candidate
     {"Data/hCPAxy", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/signal/hCPAxyRecCandSig", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/prompt/hCPAxyRecCandSigPrompt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/nonprompt/hCPAxyRecCandSigNonPrompt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/background/hCPAxyRecCandBkg", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},

     /// Chi2 PCA to sec. vertex
     {"Data/hDca2", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/signal/hDca2RecCandSig", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/signal/hDca2RecCandSigPrompt", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/signal/hDca2RecCandSigNonPrompt", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     {"MC/reconstructed/background/hDca2RecCandBkg", "3-prong candidates (unmatched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{100, 0, 0.5}}}},
     
     /// eta
     {"Data/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/signal/hEtaRecCandSig", "3-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/prompt/hEtaRecCandSigPrompt", "3-prong candidates (matched, prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/nonprompt/hEtaRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/background/hEtaRecCandBkg", "3-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/signal/hEtaGenCand", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/prompt/hEtaGenCandPrompt", "MC particles (matched, prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/nonprompt/hEtaGenCandNonPrompt", "MC particles (matched, non-prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     
     /// y
    //  {"Data/hY", "3-prong candidates;#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
    //  {"MC/reconstructed/signal/hYRecCandSig", "3-prong candidates (matched);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
    //  {"MC/reconstructed/prompt/hYRecCandSigPrompt", "3-prong candidates (matched, prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
    //  {"MC/reconstructed/nonprompt/hYRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
    //  {"MC/reconstructed/background/hYRecCandBkg", "3-prong candidates (unmatched);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/signal/hYGenCand", "MC particles (matched);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/prompt/hYGenCandPrompt", "MC particles (matched, prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/nonprompt/hYGenCandNonPrompt", "MC particles (matched, non-prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},

     /// phi
     {"Data/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/signal/hPhiRecCandSig", "3-prong candidates (matched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/prompt/hPhiRecCandSigPrompt", "3-prong candidates (matched, prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/nonprompt/hPhiRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/background/hPhiRecCandBkg", "3-prong candidates (unmatched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/signal/hPhiGenCand", "MC particles (matched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/prompt/hPhiGenCandPrompt", "MC particles (matched, prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/nonprompt/hPhiGenCandNonPrompt", "MC particles (matched, non-prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}}}};

  Configurable<int> d_selectionFlagLc{"d_selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_lc_topkpi::pTBins_v}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_selcandidate_lc::isSelLcpKpi >= d_selectionFlagLc || aod::hf_selcandidate_lc::isSelLcpiKp >= d_selectionFlagLc);

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    /// mass candidate
    registry.add("Data/hMassVsPtVsMult", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {5000, 0., 10000.}}});
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hMassVsPtRecCandSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hMassVsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hMassVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hMassVsPtRecCandBkg", "3-prong candidates (unmatched);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// DCAxy to prim. vertex prongs
    registry.add("Data/hd0VsPtProng0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng0Sig", "3-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0VsPtRecProng0Bkg", "3-prong candidates (unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng1Sig", "3-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0VsPtRecProng1Bkg", "3-prong candidates (unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng2Sig", "3-prong candidates (matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0VsPtRecProng2Bkg", "3-prong candidates (unmatched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.2, 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// decay length candidate
    registry.add("Data/hDecLengthVsPt", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthVsPtRecCandSig", "3-prong candidates (matched);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLengthVsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLengthVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLengthVsPtRecCandBkg", "3-prong candidates (unmatched);decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// decay length xy candidate
    registry.add("Data/hDecLengthxyVsPt", "3-prong candidates;decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthxyVsPtRecCandSig", "3-prong candidates (matched);decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLengthxyVsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLengthxyVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLengthxyVsPtRecCandBkg", "3-prong candidates (unmatched);decay length xy(cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// proper lifetime
    registry.add("Data/hCtVsPt", "3-prong candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCtVsPtRecCandSig", "3-prong candidates (matched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCtVsPtRecCandSigPromot", "3-prong candidates (matched, prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonpromptv/hCtVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCtVsPtRecCandBkg", "3-prong candidates (unmatched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// cosine of pointing angle
    registry.add("Data/hCPAVsPt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAVsPtRecCandSig", "3-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCPAVsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonpromptv/hCPAVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCPAVsPtRecCandBkg", "3-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// cosine of pointing angle xy
    registry.add("Data/hCPAxyVsPt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAxyVsPtRecCandSig", "3-prong candidates (matched);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCPAxyVsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonpromptv/hCPAxyVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCPAxyVsPtRecCandBkg", "3-prong candidates (unmatched);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// Chi 2 PCA to sec. vertex
    registry.add("Data/hDca2VsPt", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDca2VsPtRecCandSig", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDca2VsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonpromptv/hDca2VsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDca2VsPtRecCandBkg", "3-prong candidates (unmatched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// eta
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hEtaVsPtRecCandSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hEtaVsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonpromptv/hEtaVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hEtaVsPtRecCandBkg", "3-prong candidates (unmatched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hEtaVsPtGenCandSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hEtaVsPtGenCandSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonpromptv/hEtaVsPtGenCandSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// y
    // registry.add("Data/hYVsPt", "3-prong candidates;candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // registry.add("MC/reconstructed/signal/hYVsPtRecCandSig", "3-prong candidates (matched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // registry.add("MC/reconstructed/prompt/hYVsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // registry.add("MC/reconstructed/nonprompt/hYVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // registry.add("MC/reconstructed/background/hYVsPtRecCandBkg", "3-prong candidates (unmatched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hYVsPtGenCandSig", "3-prong candidates (matched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hYVsPtGenCandSigPrompt", "3-prong candidates (matched, prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonpromptv/hYVsPtGenCandSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// phi
    registry.add("Data/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPhiVsPtRecCandSig", "3-prong candidates (matched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hPhiVsPtRecCandSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonpromptv/hPhiVsPtRecCandSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hPhiVsPtRecCandBkg", "3-prong candidates (unmatched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hPhiVsPtGenCandSig", "3-prong candidates (matched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hPhiVsPtGenCandSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonpromptv/hPhiVsPtGenCandSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});


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
      registry.fill(HIST("Data/hPtCand"), candidate.pt());
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
      registry.fill(HIST("Data/hPhi"), candidate.phi(), candidate.pt());
      registry.fill(HIST("Data/hPhiVsPt"), candidate.phi());
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
        registry.fill(HIST("MC/generated/hPtGenCandSig"), particleMother.pt()); // gen. level pT

        /// MC reconstructed signal
        if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
          registry.fill(HIST("MC/reconstructed/signal/hMassRecCandSig"), InvMassLcpKpi(candidate));
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecCandSig"), InvMassLcpKpi(candidate), candidate.pt());

        }
        if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
          registry.fill(HIST("MC/reconstructed/signal/hMassRecCandSig"), InvMassLcpiKp(candidate));
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecCandSig"), InvMassLcpiKp(candidate), candidate.pt());
        }
        registry.fill(HIST("MC/reconstructed/signal/hPtRecCandSig"), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng0Sig"), candidate.ptProng0()); 
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng1Sig"), candidate.ptProng1()); 
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng2Sig"), candidate.ptProng2());

        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng0Sig"), candidate.impactParameter0());
        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng1Sig"), candidate.impactParameter1());
        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng2Sig"), candidate.impactParameter2());
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng0Sig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng1Sig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng2Sig"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthRecCandSig"), candidate.decayLength());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthVsPtRecCandSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthxyRecCandSig"), candidate.decayLengthXY());
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthxyVsPtRecCandSig"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hCtRecCandSig"), CtLc(candidate));
        registry.fill(HIST("MC/reconstructed/signal/hCtVsPtRecCandSig"), CtLc(candidate), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hCPARecCandSig"), candidate.cpa());
        registry.fill(HIST("MC/reconstructed/signal/hCPAVsPtRecCandSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hCPAxyRecCandSig"), candidate.cpaXY());
        registry.fill(HIST("MC/reconstructed/signal/hCPAxyVsPtRecCandSig"), candidate.cpaXY(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hDca2RecCandSig"), candidate.chi2PCA());
        registry.fill(HIST("MC/reconstructed/signal/hDca2VsPtRecCandSig"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecCandSig"), candidate.eta());
        registry.fill(HIST("MC/reconstructed/signal/hEtaVsPtRecCandSig"), candidate.eta(), candidate.pt());
        // registry.fill(HIST("MC/reconstructed/signal/hYRecCandSig"), candidate.y());
        // registry.fill(HIST("MC/reconstructed/signal/hYVsPtRecCandSig"), candidate.y(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hPhiRecCandSig"), candidate.phi(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hPhiVsPtRecCandSig"), candidate.phi());



        /// reconstructed signal prompt
        if (candidate.originMCRec() == RecoDecay::OriginType::Prompt) {
          if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecCandSigPrompt"), InvMassLcpKpi(candidate));
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecCandSigPrompt"), InvMassLcpKpi(candidate), candidate.pt());

          }
          if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecCandSigPrompt"), InvMassLcpiKp(candidate));
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecCandSigPrompt"), InvMassLcpiKp(candidate), candidate.pt());
          }
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecCandSigPrompt"), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng0SigPrompt"), candidate.ptProng0());
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng1SigPrompt"), candidate.ptProng1());
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng2SigPrompt"), candidate.ptProng2());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng0SigPrompt"), candidate.impactParameter0());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng1SigPrompt"), candidate.impactParameter1());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng2SigPrompt"), candidate.impactParameter2());
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng0SigPrompt"), candidate.impactParameter0(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng1SigPrompt"), candidate.impactParameter1(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng2SigPrompt"), candidate.impactParameter2(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthRecCandSigPrompt"), candidate.decayLength());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthVsPtRecCandSigPrompt"), candidate.decayLength(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthxyRecCandSigPrompt"), candidate.decayLengthXY());
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthxyVsPtRecCandSigPrompt"), candidate.decayLengthXY(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hCtRecCandSigPrompt"), CtLc(candidate));
          registry.fill(HIST("MC/reconstructed/prompt/hCtVsPtRecCandSigPrompt"), CtLc(candidate), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hCPARecCandSigPrompt"), candidate.cpa());
          registry.fill(HIST("MC/reconstructed/prompt/hCPAVsPtRecCandSigPrompt"), candidate.cpa(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hCPAxyRecCandSigPrompt"), candidate.cpaXY());
          registry.fill(HIST("MC/reconstructed/prompt/hCPAxyVsPtRecCandSigPrompt"), candidate.cpaXY(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hDca2RecCandSigPrompt"), candidate.chi2PCA());
          registry.fill(HIST("MC/reconstructed/prompt/hDca2VsPtRecCandSigPrompt"), candidate.chi2PCA(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hEtaRecCandSigPrompt"), candidate.eta());
          registry.fill(HIST("MC/reconstructed/prompt/hEtaVsPtRecCandSigPrompt"), candidate.eta(), candidate.pt());
          // registry.fill(HIST("MC/reconstructed/prompt/hYRecCandSigPrompt"), candidate.y());
          // registry.fill(HIST("MC/reconstructed/prompt/hYVsPtRecCandSigPrompt"), candidate.y(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hPhiRecCandSigPrompt"), candidate.phi(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/prompt/hPhiVsPtRecCandSigPrompt"), candidate.phi());     
        } 
        /// reconstructed signal non-prompt
        if (candidate.originMCRec() == RecoDecay::OriginType::NonPrompt) {
          if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecCandSigNonPrompt"), InvMassLcpKpi(candidate));
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecCandSigNonPrompt"), InvMassLcpKpi(candidate), candidate.pt());

          }
          if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecCandSigNonPrompt"), InvMassLcpiKp(candidate));
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecCandSigNonPrompt"), InvMassLcpiKp(candidate), candidate.pt());
          }
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecCandSigNonPrompt"), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng0SigNonPrompt"), candidate.ptProng0());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng1SigNonPrompt"), candidate.ptProng1());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng2SigNonPrompt"), candidate.ptProng2());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng0SigNonPrompt"), candidate.impactParameter0());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng1SigNonPrompt"), candidate.impactParameter1());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng2SigNonPrompt"), candidate.impactParameter2());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng0SigNonPrompt"), candidate.impactParameter0(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng1SigNonPrompt"), candidate.impactParameter1(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng2SigNonPrompt"), candidate.impactParameter2(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthRecCandSigNonPrompt"), candidate.decayLength());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthVsPtRecCandSigNonPrompt"), candidate.decayLength(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthxyRecCandSigNonPrompt"), candidate.decayLengthXY());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthxyVsPtRecCandSigNonPrompt"), candidate.decayLengthXY(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCtRecCandSigNonPrompt"), CtLc(candidate));
          registry.fill(HIST("MC/reconstructed/nonprompt/hCtVsPtRecCandSigNonPrompt"), CtLc(candidate), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPARecCandSigNonPrompt"), candidate.cpa());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAVsPtRecCandSigNonPrompt"), candidate.cpa(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAxyRecCandSigNonPrompt"), candidate.cpaXY());
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAxyVsPtRecCandSigNonPrompt"), candidate.cpaXY(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDca2RecCandSigNonPrompt"), candidate.chi2PCA());
          registry.fill(HIST("MC/reconstructed/nonprompt/hDca2VsPtRecCandSigNonPrompt"), candidate.chi2PCA(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaRecCandSigNonPrompt"), candidate.eta());
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaVsPtRecCandSigNonPrompt"), candidate.eta(), candidate.pt());
          // registry.fill(HIST("MC/reconstructed/nonprompt/hYRecCandSigNonPrompt"), candidate.y());
          // registry.fill(HIST("MC/reconstructed/nonprompt/hYVsPtRecCandSigNonPrompt"), candidate.y(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPhiRecCandSigNonPrompt"), candidate.phi(), candidate.pt());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPhiVsPtRecCandSigNonPrompt"), candidate.phi());     
        }

      } else {
        if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
          registry.fill(HIST("MC/reconstructed/background/hMassRecCandBkg"), InvMassLcpKpi(candidate));
          registry.fill(HIST("MC/reconstructed/background/hMassVsPtRecCandBkg"), InvMassLcpKpi(candidate), candidate.pt());

        }
        if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
          registry.fill(HIST("MC/reconstructed/background/hMassRecCandBkg"), InvMassLcpiKp(candidate));
          registry.fill(HIST("MC/reconstructed/background/hMassVsPtRecCandBkg"), InvMassLcpiKp(candidate), candidate.pt());
        }
        registry.fill(HIST("MC/reconstructed/background/hPtRecCandBkg"), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hPtRecProng0Bkg"), candidate.ptProng0());
        registry.fill(HIST("MC/reconstructed/background/hPtRecProng1Bkg"), candidate.ptProng1());
        registry.fill(HIST("MC/reconstructed/background/hPtRecProng2Bkg"), candidate.ptProng2());
        registry.fill(HIST("MC/reconstructed/background/hd0RecProng0Bkg"), candidate.impactParameter0());
        registry.fill(HIST("MC/reconstructed/background/hd0RecProng1Bkg"), candidate.impactParameter1());
        registry.fill(HIST("MC/reconstructed/background/hd0RecProng2Bkg"), candidate.impactParameter2());
        registry.fill(HIST("MC/reconstructed/background/hd0VsPtRecProng0Bkg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hd0VsPtRecProng1Bkg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hd0VsPtRecProng2Bkg"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hDecLengthRecCandBkg"), candidate.decayLength());
        registry.fill(HIST("MC/reconstructed/background/hDecLengthVsPtRecCandBkg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hDecLengthxyRecCandBkg"), candidate.decayLengthXY());
        registry.fill(HIST("MC/reconstructed/background/hDecLengthxyVsPtRecCandBkg"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hCtRecCandBkg"), CtLc(candidate));
        registry.fill(HIST("MC/reconstructed/background/hCtVsPtRecCandBkg"), CtLc(candidate), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hCPARecCandBkg"), candidate.cpa());
        registry.fill(HIST("MC/reconstructed/background/hCPAVsPtRecCandBkg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hCPAxyRecCandBkg"), candidate.cpaXY());
        registry.fill(HIST("MC/reconstructed/background/hCPAxyVsPtRecCandBkg"), candidate.cpaXY(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hDca2RecCandBkg"), candidate.chi2PCA());
        registry.fill(HIST("MC/reconstructed/background/hDca2VsPtRecCandBkg"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hEtaRecCandBkg"), candidate.eta());
        registry.fill(HIST("MC/reconstructed/background/hEtaVsPtRecCandBkg"), candidate.eta(), candidate.pt());
        // registry.fill(HIST("MC/reconstructed/background/hYRecCandBkg"), RecoDecay::y(array{candidate.px(), candidate.py(), candidate.pz()}, RecoDecay::getMassPDG(candidate.pdgCode())));
        // registry.fill(HIST("MC/reconstructed/background/hYVsPtRecCandBkg"), RecoDecay::y(array{candidate.px(), candidate.py(), candidate.pz()}, RecoDecay::getMassPDG(candidate.pdgCode())), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hPhiRecCandBkg"), candidate.phi(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/background/hPhiVsPtRecCandBkg"), candidate.phi());  
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
        registry.fill(HIST("MC/generated/signal/hPtGenCand"), ptGen);
        registry.fill(HIST("MC/generated/signal/hEtaGenCand"), particle.eta());
        registry.fill(HIST("MC/generated/signal/hYGenCand"), yGen);
        registry.fill(HIST("MC/generated/signal/hPhiGenCand"), particle.phi());
        registry.fill(HIST("MC/generated/signal/hEtaVsPtGenCandSig"), particle.eta(), ptGen);
        registry.fill(HIST("MC/generated/signal/hYVsPtGenCandSig"), yGen, ptGen);
        registry.fill(HIST("MC/generated/signal/hPhiVsPtGenCandSig"), particle.phi(), ptGen);

        if (particle.originMCGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/generated/prompt/hPtGenCandPrompt"), ptGen);
          registry.fill(HIST("MC/generated/prompt/hEtaGenCandPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/prompt/hYGenCandPrompt"), yGen);
          registry.fill(HIST("MC/generated/prompt/hPhiGenCandPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/prompt/hEtaVsPtGenCandSigPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/prompt/hYVsPtGenCandSigPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/prompt/hPhiVsPtGenCandSigPrompt"), particle.phi(), ptGen);
        } 
        if (particle.originMCGen() == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("MC/generated/nonprompt/hPtGenCandNonPrompt"), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hEtaGenCandNonPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/nonprompt/hYGenCandNonPrompt"), yGen);
          registry.fill(HIST("MC/generated/nonprompt/hPhiGenCandNonPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/nonprompt/hEtaVsPtGenCandSigNonPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hYVsPtGenCandSigNonPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hPhiVsPtGenCandSigNonPrompt"), particle.phi(), ptGen);
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
