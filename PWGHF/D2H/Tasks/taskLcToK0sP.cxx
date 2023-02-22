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

/// \file taskLcToK0sP.cxx
/// \brief Lc -> K0S+p analysis task
///
/// \author Chiara Zampolli, <Chiara.Zampolli@cern.ch>, CERN
///         Paul Buehler, <paul.buehler@oeaw.ac.at>, Vienna
///
/// based on taskD0.cxx, taskLc.cxx

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_casc;
using namespace o2::framework::expressions;

/// LcToK0sp analysis task
struct HfTaskLcToK0sP {
  Configurable<int> selectionFlagLcToK0sP{"selectionFlagLcToK0sP", 1, "Selection Flag for Lc"};
  Configurable<int> selectionFlagLcbarToK0sP{"selectionFlagLcbarToK0sP", 1, "Selection Flag for Lcbar"};
  Configurable<double> etaCandMax{"etaCandMax", -1., "max. cand. pseudorapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcToK0sP || aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcbarToK0sP);

  HistogramRegistry registry{"registry"};

  void init(InitContext& context)
  {
    // data
    registry.add("hPtCand", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
    registry.add("hEtaCand", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{500, -2.0f, 2.0f}}});
    registry.add("hEtaCandVsPtCand", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {{500, -2.0f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPhiCand", "cascade candidates;candidate #it{#phi};entries", {HistType::kTH1F, {{100, 0.f, 6.3f}}});
    registry.add("hPhiCandVsPtCand", "cascade candidates;candidate #it{#phi};p_{T}", {HistType::kTH2F, {{100, 0.f, 6.3f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMass", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{600, 1.98f, 2.58f}}});
    registry.add("hMassVsPtCand", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{600, 1.98f, 2.58f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtBach", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
    registry.add("hPtBachVsPtCand", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtV0", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
    registry.add("hPtV0VsPtCand", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Bach", "cascade candidates;bachelor DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{500, -0.5f, 0.5f}}});
    registry.add("hd0BachVsPtCand", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{500, -0.5f, 0.5f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0V0", "cascade candidates;V0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{500, -0.5f, 0.5f}}});
    registry.add("hd0V0VsPtCand", "cascade candidates;V0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{500, -0.5f, 0.5f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0V0pos", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{1000, -5.0f, 5.0f}}});
    registry.add("hd0V0posVsPtCand", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0V0neg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{1000, -5.0f, 5.0f}}});
    registry.add("hd0V0negVsPtCand", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtV0pos", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
    registry.add("hPtV0posVsPtCand", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtV0neg", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
    registry.add("hPtV0negVsPtCand", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0CPA", "cascade candidates;v0 cosine of pointing angle;entries", {HistType::kTH1F, {{500, 0.98f, 1.0001f}}});
    registry.add("hV0CPAVsPtCand", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98f, 1.0001f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0Radius", "cascade candidates;v0 radius (cm);entries", {HistType::kTH1F, {{1000, 0.f, 40.f}}});
    registry.add("hV0RadiusVsPtCand", "cascade candidates;v0 radius (cm);p_{T}", {HistType::kTH2F, {{1000, 0.f, 40.f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0DCADaughters", "cascade candidates;v0 dca daughters (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.f}}});
    registry.add("hV0DCADaughtersVsPtCand", "cascade candidates;v0 dca daughters (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0MK0Short", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0.4f, 0.6f}}});
    registry.add("hV0MK0ShortVsPtCand", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 0.4f, 0.6f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0MLambda", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.0f, 1.2f}}});
    registry.add("hV0MLambdaVsPtCand", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 1.0f, 1.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0MAntiLambda", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.0f, 1.2f}}});
    registry.add("hV0MAntiLambdaVsPtCand", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 1.0f, 1.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0MGamma", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0.0f, 0.4f}}});
    registry.add("hV0MGammaVsPtCand", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 0.0f, 0.4f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPACand", "cascade candidates;cosine pointing angle;entries", {HistType::kTH1F, {{110, -1.1f, 1.1f}}});
    registry.add("hCPACandVsPtCand", "cascade candidates;cosine pointing angle;p_{T}", {HistType::kTH2F, {{110, -1.1f, 1.1f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyCand", "cascade candidates;cosine pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1f, 1.1f}}});
    registry.add("hCPAxyCandVsPtCand", "cascade candidates;cosine pointing angle xy;p_{T}", {HistType::kTH2F, {{110, -1.1f, 1.1f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthCand", "cascade candidates;decay length (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.0f}}});
    registry.add("hDecLengthCandVsPtCand", "cascade candidates;decay length (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYCand", "cascade candidates;decay length xy (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.0f}}});
    registry.add("hDecLengthXYCandVsPtCand", "cascade candidates;decay length xy (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtCand", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0.f, 0.2f}}});
    registry.add("hCtCandVsPtCand", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);p_{T}", {HistType::kTH2F, {{100, 0.f, 0.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hTPCNSigmaPrBach", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {{100, -6.f, 6.f}}});
    registry.add("hPBachVsTPCNSigmaPrBach", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {{100, 0.f, 10.0f}, {100, -6.f, 6.f}}});
    registry.add("hTOFNSigmaPrBach", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {{100, -6.f, 6.f}}});
    registry.add("hPBachVsTOFNSigmaPrBach", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {{100, 0.f, 10.0f}, {100, -6.f, 6.f}}});

    // add MC histograms
    if (context.mOptions.get<bool>("processMc")) {
      registry.add("MC/Rec/hPtCandRecSig", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtCandRecBg", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hEtaCandRecSig", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{500, -2.0f, 2.0f}}});
      registry.add("MC/Rec/hEtaCandVsPtCandRecSig", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {{500, -2.0f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hEtaCandRecBg", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{500, -2.0f, 2.0f}}});
      registry.add("MC/Rec/hEtaCandVsPtCandRecBg", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {{500, -2.0f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPhiCandRecSig", "cascade candidates;candidate #it{#phi};entries", {HistType::kTH1F, {{100, 0.f, 6.3f}}});
      registry.add("MC/Rec/hPhiCandVsPtCandRecSig", "cascade candidates;candidate #it{#phi};p_{T}", {HistType::kTH2F, {{100, 0.f, 6.3f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPhiCandRecBg", "cascade candidates;candidate #it{#phi};entries", {HistType::kTH1F, {{100, 0.f, 6.3f}}});
      registry.add("MC/Rec/hPhiCandVsPtCandRecBg", "cascade candidates;candidate #it{#phi};p_{T}", {HistType::kTH2F, {{100, 0.f, 6.3f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Gen/hPtCandGen", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Gen/hEtaCandGen", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{500, -2.0f, 2.0f}}});
      registry.add("MC/Gen/hEtaCandVsPtCandGen", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {{500, -2.0f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Gen/hPhiCandGen", "cascade candidates;candidate #it{#phi};entries", {HistType::kTH1F, {{100, 0.f, 6.3f}}});
      registry.add("MC/Gen/hPhiCandVsPtCandGen", "cascade candidates;candidate #it{#phi};p_{T}", {HistType::kTH2F, {{100, 0.f, 6.3f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hMassRecSig", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH1F, {{600, 1.98f, 2.58f}}});
      registry.add("MC/Rec/hMassVsPtCandRecSig", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{600, 1.98f, 2.58f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hMassRecBg", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH1F, {{600, 1.98f, 2.58f}}});
      registry.add("MC/Rec/hMassVsPtCandRecBg", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{600, 1.98f, 2.58f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPtBachRecSig", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtBachVsPtCandRecSig", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPtBachRecBg", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtBachVsPtCandRecBg", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPtV0RecSig", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtV0VsPtCandRecSig", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPtV0RecBg", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtV0VsPtCandRecBg", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hd0BachRecSig", "cascade candidates;bachelor DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{500, -0.5f, 0.5f}}});
      registry.add("MC/Rec/hd0BachVsPtCandRecSig", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{500, -0.5f, 0.5f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hd0BachRecBg", "cascade candidates;bachelor DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{500, -0.5f, 0.5f}}});
      registry.add("MC/Rec/hd0BachVsPtCandRecBg", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{500, -0.5f, 0.5f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hd0V0RecSig", "cascade candidates;V0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{500, -0.5f, 0.5f}}});
      registry.add("MC/Rec/hd0V0VsPtCandRecSig", "cascade candidates;V0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{500, -0.5f, 0.5f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hd0V0RecBg", "cascade candidates;V0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{500, -0.5f, 0.5f}}});
      registry.add("MC/Rec/hd0V0VsPtCandRecBg", "cascade candidates;V0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{500, -0.5f, 0.5f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hd0V0posRecSig", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{1000, -5.0f, 5.0f}}});
      registry.add("MC/Rec/hd0V0posVsPtCandRecSig", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hd0V0posRecBg", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{1000, -5.0f, 5.0f}}});
      registry.add("MC/Rec/hd0V0posVsPtCandRecBg", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hd0V0negRecSig", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{1000, -5.0f, 5.0f}}});
      registry.add("MC/Rec/hd0V0negVsPtCandRecSig", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hd0V0negRecBg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{1000, -5.0f, 5.0f}}});
      registry.add("MC/Rec/hd0V0negVsPtCandRecBg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPtV0posRecSig", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtV0posVsPtCandRecSig", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPtV0posRecBg", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtV0posVsPtCandRecBg", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPtV0negRecSig", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtV0negVsPtCandRecSig", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hPtV0negRecBg", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("MC/Rec/hPtV0negVsPtCandRecBg", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{300, 0.0f, 30.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0CPARecSig", "cascade candidates;v0 cosine of pointing angle;entries", {HistType::kTH1F, {{500, 0.98f, 1.0001f}}});
      registry.add("MC/Rec/hV0CPAVsPtCandRecSig", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98f, 1.0001f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0CPARecBg", "cascade candidates;v0 cosine of pointing angle;entries", {HistType::kTH1F, {{500, 0.98f, 1.0001f}}});
      registry.add("MC/Rec/hV0CPAVsPtCandRecBg", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98f, 1.0001f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0RadiusRecSig", "cascade candidates;v0 radius (cm);entries", {HistType::kTH1F, {{1000, 0.f, 40.f}}});
      registry.add("MC/Rec/hV0RadiusVsPtCandRecSig", "cascade candidates;v0 radius (cm);p_{T}", {HistType::kTH2F, {{1000, 0.f, 40.f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0RadiusRecBg", "cascade candidates;v0 radius (cm);entries", {HistType::kTH1F, {{1000, 0.f, 40.f}}});
      registry.add("MC/Rec/hV0RadiusVsPtCandRecBg", "cascade candidates;v0 radius (cm);p_{T}", {HistType::kTH2F, {{1000, 0.f, 40.f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0DCADaughtersRecSig", "cascade candidates;v0 dca daughters (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.f}}});
      registry.add("MC/Rec/hV0DCADaughtersVsPtCandRecSig", "cascade candidates;v0 dca daughters (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0DCADaughtersRecBg", "cascade candidates;v0 dca daughters (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.f}}});
      registry.add("MC/Rec/hV0DCADaughtersVsPtCandRecBg", "cascade candidates;v0 dca daughters (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0MK0ShortRecSig", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0.4f, 0.6f}}});
      registry.add("MC/Rec/hV0MK0ShortVsPtCandRecSig", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 0.4f, 0.6f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0MK0ShortRecBg", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0.4f, 0.6f}}});
      registry.add("MC/Rec/hV0MK0ShortVsPtCandRecBg", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 0.4f, 0.6f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0MLambdaRecSig", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.0f, 1.2f}}});
      registry.add("MC/Rec/hV0MLambdaVsPtCandRecSig", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 1.0f, 1.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0MLambdaRecBg", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.0f, 1.2f}}});
      registry.add("MC/Rec/hV0MLambdaVsPtCandRecBg", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 1.0f, 1.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0MAntiLambdaRecSig", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.0f, 1.2f}}});
      registry.add("MC/Rec/hV0MAntiLambdaVsPtCandRecSig", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 1.0f, 1.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0MAntiLambdaRecBg", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.0f, 1.2f}}});
      registry.add("MC/Rec/hV0MAntiLambdaVsPtCandRecBg", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 1.0f, 1.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0MGammaRecSig", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0.0f, 0.4f}}});
      registry.add("MC/Rec/hV0MGammaVsPtCandRecSig", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 0.0f, 0.4f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hV0MGammaRecBg", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0.0f, 0.4f}}});
      registry.add("MC/Rec/hV0MGammaVsPtCandRecBg", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{500, 0.0f, 0.4f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hCPACandRecSig", "cascade candidates;cosine pointing angle;entries", {HistType::kTH1F, {{110, -1.1f, 1.1f}}});
      registry.add("MC/Rec/hCPACandVsPtCandRecSig", "cascade candidates;cosine pointing angle;p_{T}", {HistType::kTH2F, {{110, -1.1f, 1.1f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hCPACandRecBg", "cascade candidates;cosine pointing angle;entries", {HistType::kTH1F, {{110, -1.1f, 1.1f}}});
      registry.add("MC/Rec/hCPACandVsPtCandRecBg", "cascade candidates;cosine pointing angle;p_{T}", {HistType::kTH2F, {{110, -1.1f, 1.1f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hCPAxyCandRecSig", "cascade candidates;cosine pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1f, 1.1f}}});
      registry.add("MC/Rec/hCPAxyCandVsPtCandRecSig", "cascade candidates;cosine pointing angle xy;p_{T}", {HistType::kTH2F, {{110, -1.1f, 1.1f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hCPAxyCandRecBg", "cascade candidates;cosine pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1f, 1.1f}}});
      registry.add("MC/Rec/hCPAxyCandVsPtCandRecBg", "cascade candidates;cosine pointing angle xy;p_{T}", {HistType::kTH2F, {{110, -1.1f, 1.1f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hDecLengthCandRecSig", "cascade candidates;decay length (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.0f}}});
      registry.add("MC/Rec/hDecLengthCandVsPtCandRecSig", "cascade candidates;decay length (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hDecLengthCandRecBg", "cascade candidates;decay length (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.0f}}});
      registry.add("MC/Rec/hDecLengthCandVsPtCandRecBg", "cascade candidates;decay length (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hDecLengthXYCandRecSig", "cascade candidates;decay length xy (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.0f}}});
      registry.add("MC/Rec/hDecLengthXYCandVsPtCandRecSig", "cascade candidates;decay length xy (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hDecLengthXYCandRecBg", "cascade candidates;decay length xy (cm);entries", {HistType::kTH1F, {{200, 0.f, 2.0f}}});
      registry.add("MC/Rec/hDecLengthXYCandVsPtCandRecBg", "cascade candidates;decay length xy (cm);p_{T}", {HistType::kTH2F, {{200, 0.f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hCtCandRecSig", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0.f, 0.2f}}});
      registry.add("MC/Rec/hCtCandVsPtCandRecSig", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);p_{T}", {HistType::kTH2F, {{100, 0.f, 0.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hCtCandRecBg", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0.f, 0.2f}}});
      registry.add("MC/Rec/hCtCandVsPtCandRecBg", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);p_{T}", {HistType::kTH2F, {{100, 0.f, 0.2f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/Rec/hTPCNSigmaPrBachRecSig", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {{100, -6.f, 6.f}}});
      registry.add("MC/Rec/hPBachVsTPCNSigmaPrBachRecSig", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {{100, 0.f, 10.0f}, {100, -6.f, 6.f}}});
      registry.add("MC/Rec/hTPCNSigmaPrBachRecBg", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {{100, -6.f, 6.f}}});
      registry.add("MC/Rec/hPBachVsTPCNSigmaPrBachRecBg", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {{100, 0.f, 10.0f}, {100, -6.f, 6.f}}});
      registry.add("MC/Rec/hTOFNSigmaPrBachRecSig", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {{100, -6.f, 6.f}}});
      registry.add("MC/Rec/hPBachVsTOFNSigmaPrBachRecSig", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {{100, 0.f, 10.0f}, {100, -6.f, 6.f}}});
      registry.add("MC/Rec/hTOFNSigmaPrBachRecBg", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {{100, -6.f, 6.f}}});
      registry.add("MC/Rec/hPBachVsTOFNSigmaPrBachRecBg", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {{100, 0.f, 10.0f}, {100, -6.f, 6.f}}});
    }
  }

  void
    process(soa::Filtered<soa::Join<aod::HfCandCascExt, aod::HfSelLcToK0sP>> const& candidates, aod::BigTracksPID const&)
  {
    // Printf("Candidates: %d", candidates.size());
    for (auto& candidate : candidates) {
      /*
      // no such selection for LcK0sp for now - it is the only cascade
      if (!(candidate.hfflag() & 1 << D0ToPiK)) {
        continue;
      }
      */

      if (etaCandMax >= 0. && std::abs(candidate.eta()) > etaCandMax) {
        // Printf("Candidate: eta rejection: %g", candidate.eta());
        continue;
      }

      auto ptCand = candidate.pt();
      registry.fill(HIST("hPtCand"), ptCand);
      registry.fill(HIST("hEtaCand"), candidate.eta());
      registry.fill(HIST("hEtaCandVsPtCand"), candidate.eta(), ptCand);
      registry.fill(HIST("hPhiCand"), candidate.phi());
      registry.fill(HIST("hPhiCandVsPtCand"), candidate.phi(), ptCand);
      registry.fill(HIST("hMass"), invMassLcToK0sP(candidate));
      registry.fill(HIST("hMassVsPtCand"), invMassLcToK0sP(candidate), ptCand);
      registry.fill(HIST("hPtBach"), candidate.ptProng0());
      registry.fill(HIST("hPtBachVsPtCand"), candidate.ptProng0(), ptCand);
      registry.fill(HIST("hPtV0"), candidate.ptProng1());
      registry.fill(HIST("hPtV0VsPtCand"), candidate.ptProng1(), ptCand);
      registry.fill(HIST("hd0Bach"), candidate.impactParameter0());
      registry.fill(HIST("hd0BachVsPtCand"), candidate.impactParameter0(), ptCand);
      registry.fill(HIST("hd0V0"), candidate.impactParameter1());
      registry.fill(HIST("hd0V0VsPtCand"), candidate.impactParameter1(), ptCand);
      registry.fill(HIST("hd0V0pos"), candidate.dcapostopv());
      registry.fill(HIST("hd0V0posVsPtCand"), candidate.dcapostopv(), ptCand);
      registry.fill(HIST("hd0V0neg"), candidate.dcanegtopv());
      registry.fill(HIST("hd0V0negVsPtCand"), candidate.dcanegtopv(), ptCand);
      registry.fill(HIST("hPtV0pos"), candidate.ptV0Pos());
      registry.fill(HIST("hPtV0posVsPtCand"), candidate.ptV0Pos(), ptCand);
      registry.fill(HIST("hPtV0neg"), candidate.ptV0Neg());
      registry.fill(HIST("hPtV0negVsPtCand"), candidate.ptV0Neg(), ptCand);
      registry.fill(HIST("hV0CPA"), candidate.v0cosPA());
      registry.fill(HIST("hV0CPAVsPtCand"), candidate.v0cosPA(), ptCand);
      registry.fill(HIST("hV0Radius"), candidate.v0radius());
      registry.fill(HIST("hV0RadiusVsPtCand"), candidate.v0radius(), ptCand);
      registry.fill(HIST("hV0DCADaughters"), candidate.dcaV0daughters());
      registry.fill(HIST("hV0DCADaughtersVsPtCand"), candidate.dcaV0daughters(), ptCand);
      registry.fill(HIST("hV0MK0Short"), candidate.mK0Short());
      registry.fill(HIST("hV0MK0ShortVsPtCand"), candidate.mK0Short(), ptCand);
      registry.fill(HIST("hV0MLambda"), candidate.mLambda());
      registry.fill(HIST("hV0MLambdaVsPtCand"), candidate.mLambda(), ptCand);
      registry.fill(HIST("hV0MAntiLambda"), candidate.mAntiLambda());
      registry.fill(HIST("hV0MAntiLambdaVsPtCand"), candidate.mAntiLambda(), ptCand);
      registry.fill(HIST("hV0MGamma"), candidate.mGamma());
      registry.fill(HIST("hV0MGammaVsPtCand"), candidate.mGamma(), ptCand);
      registry.fill(HIST("hCPACand"), candidate.cpa());
      registry.fill(HIST("hCPACandVsPtCand"), candidate.cpa(), ptCand);
      registry.fill(HIST("hCPAxyCand"), candidate.cpaXY());
      registry.fill(HIST("hCPAxyCandVsPtCand"), candidate.cpaXY(), ptCand);
      registry.fill(HIST("hDecLengthCand"), candidate.decayLength());
      registry.fill(HIST("hDecLengthCandVsPtCand"), candidate.decayLength(), ptCand);
      registry.fill(HIST("hDecLengthXYCand"), candidate.decayLengthXY());
      registry.fill(HIST("hDecLengthXYCandVsPtCand"), candidate.decayLengthXY(), ptCand);
      registry.fill(HIST("hCtCand"), o2::aod::hf_cand_3prong::ctLc(candidate));
      registry.fill(HIST("hCtCandVsPtCand"), o2::aod::hf_cand_3prong::ctLc(candidate), ptCand);
      const auto& bach = candidate.prong0_as<aod::BigTracksPID>(); // bachelor track
      registry.fill(HIST("hTPCNSigmaPrBach"), bach.tpcNSigmaPr());
      registry.fill(HIST("hPBachVsTPCNSigmaPrBach"), bach.p(), bach.tpcNSigmaPr());
      if (bach.hasTOF()) {
        registry.fill(HIST("hTOFNSigmaPrBach"), bach.tofNSigmaPr());
        registry.fill(HIST("hPBachVsTOFNSigmaPrBach"), bach.p(), bach.tofNSigmaPr());
      }
    }
  }

  void processMc(soa::Filtered<soa::Join<aod::HfCandCascExt, aod::HfSelLcToK0sP, aod::HfCandCascadeMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& particlesMC,
                 aod::BigTracksMC const& tracks, aod::BigTracksPID const&)
  {
    // MC rec.
    // Printf("MC Candidates: %d", candidates.size());
    for (auto& candidate : candidates) {
      if (etaCandMax >= 0. && std::abs(candidate.eta()) > etaCandMax) {
        // Printf("MC Rec.: eta rejection: %g", candidate.eta());
        continue;
      }
      const auto& bach = candidate.prong0_as<aod::BigTracksPID>(); // bachelor track
      auto ptCand = candidate.pt();
      if (std::abs(candidate.flagMcMatchRec()) == 1) {
        registry.fill(HIST("MC/Rec/hPtCandRecSig"), ptCand);
        registry.fill(HIST("MC/Rec/hEtaCandRecSig"), candidate.eta());
        registry.fill(HIST("MC/Rec/hEtaCandVsPtCandRecSig"), candidate.eta(), ptCand);
        registry.fill(HIST("MC/Rec/hPhiCandRecSig"), candidate.phi());
        registry.fill(HIST("MC/Rec/hPhiCandVsPtCandRecSig"), candidate.phi(), ptCand);
        registry.fill(HIST("MC/Rec/hMassRecSig"), invMassLcToK0sP(candidate));
        registry.fill(HIST("MC/Rec/hMassVsPtCandRecSig"), invMassLcToK0sP(candidate), ptCand);
        registry.fill(HIST("MC/Rec/hPtBachRecSig"), candidate.ptProng0());
        registry.fill(HIST("MC/Rec/hPtBachVsPtCandRecSig"), candidate.ptProng0(), ptCand);
        registry.fill(HIST("MC/Rec/hPtV0RecSig"), candidate.ptProng1());
        registry.fill(HIST("MC/Rec/hPtV0VsPtCandRecSig"), candidate.ptProng1(), ptCand);
        registry.fill(HIST("MC/Rec/hd0BachRecSig"), candidate.impactParameter0());
        registry.fill(HIST("MC/Rec/hd0BachVsPtCandRecSig"), candidate.impactParameter0(), ptCand);
        registry.fill(HIST("MC/Rec/hd0V0RecSig"), candidate.impactParameter1());
        registry.fill(HIST("MC/Rec/hd0V0VsPtCandRecSig"), candidate.impactParameter1(), ptCand);
        registry.fill(HIST("MC/Rec/hd0V0posRecSig"), candidate.dcapostopv());
        registry.fill(HIST("MC/Rec/hd0V0posVsPtCandRecSig"), candidate.dcapostopv(), ptCand);
        registry.fill(HIST("MC/Rec/hd0V0negRecSig"), candidate.dcanegtopv());
        registry.fill(HIST("MC/Rec/hd0V0negVsPtCandRecSig"), candidate.dcanegtopv(), ptCand);
        registry.fill(HIST("MC/Rec/hPtV0posRecSig"), candidate.ptV0Pos());
        registry.fill(HIST("MC/Rec/hPtV0posVsPtCandRecSig"), candidate.ptV0Pos(), ptCand);
        registry.fill(HIST("MC/Rec/hPtV0negRecSig"), candidate.ptV0Neg());
        registry.fill(HIST("MC/Rec/hPtV0negVsPtCandRecSig"), candidate.ptV0Neg(), ptCand);
        registry.fill(HIST("MC/Rec/hV0CPARecSig"), candidate.v0cosPA());
        registry.fill(HIST("MC/Rec/hV0CPAVsPtCandRecSig"), candidate.v0cosPA(), ptCand);
        registry.fill(HIST("MC/Rec/hV0RadiusRecSig"), candidate.v0radius());
        registry.fill(HIST("MC/Rec/hV0RadiusVsPtCandRecSig"), candidate.v0radius(), ptCand);
        registry.fill(HIST("MC/Rec/hV0DCADaughtersRecSig"), candidate.dcaV0daughters());
        registry.fill(HIST("MC/Rec/hV0DCADaughtersVsPtCandRecSig"), candidate.dcaV0daughters(), ptCand);
        registry.fill(HIST("MC/Rec/hV0MK0ShortRecSig"), candidate.mK0Short());
        registry.fill(HIST("MC/Rec/hV0MK0ShortVsPtCandRecSig"), candidate.mK0Short(), ptCand);
        registry.fill(HIST("MC/Rec/hV0MLambdaRecSig"), candidate.mLambda());
        registry.fill(HIST("MC/Rec/hV0MLambdaVsPtCandRecSig"), candidate.mLambda(), ptCand);
        registry.fill(HIST("MC/Rec/hV0MAntiLambdaRecSig"), candidate.mAntiLambda());
        registry.fill(HIST("MC/Rec/hV0MAntiLambdaVsPtCandRecSig"), candidate.mAntiLambda(), ptCand);
        registry.fill(HIST("MC/Rec/hV0MGammaRecSig"), candidate.mGamma());
        registry.fill(HIST("MC/Rec/hV0MGammaVsPtCandRecSig"), candidate.mGamma(), ptCand);
        registry.fill(HIST("MC/Rec/hCPACandRecSig"), candidate.cpa());
        registry.fill(HIST("MC/Rec/hCPACandVsPtCandRecSig"), candidate.cpa(), ptCand);
        registry.fill(HIST("MC/Rec/hCPAxyCandRecSig"), candidate.cpaXY());
        registry.fill(HIST("MC/Rec/hCPAxyCandVsPtCandRecSig"), candidate.cpaXY(), ptCand);
        registry.fill(HIST("MC/Rec/hDecLengthCandRecSig"), candidate.decayLength());
        registry.fill(HIST("MC/Rec/hDecLengthCandVsPtCandRecSig"), candidate.decayLength(), ptCand);
        registry.fill(HIST("MC/Rec/hDecLengthXYCandRecSig"), candidate.decayLengthXY());
        registry.fill(HIST("MC/Rec/hDecLengthXYCandVsPtCandRecSig"), candidate.decayLengthXY(), ptCand);
        registry.fill(HIST("MC/Rec/hCtCandRecSig"), o2::aod::hf_cand_3prong::ctLc(candidate));
        registry.fill(HIST("MC/Rec/hCtCandVsPtCandRecSig"), o2::aod::hf_cand_3prong::ctLc(candidate), ptCand);
        registry.fill(HIST("MC/Rec/hTPCNSigmaPrBachRecSig"), bach.tpcNSigmaPr());
        registry.fill(HIST("MC/Rec/hPBachVsTPCNSigmaPrBachRecSig"), bach.p(), bach.tpcNSigmaPr());
        if (bach.hasTOF()) {
          registry.fill(HIST("MC/Rec/hTOFNSigmaPrBachRecSig"), bach.tofNSigmaPr());
          registry.fill(HIST("MC/Rec/hPBachVsTOFNSigmaPrBachRecSig"), bach.p(), bach.tofNSigmaPr());
        }
      } else {
        registry.fill(HIST("MC/Rec/hPtCandRecBg"), ptCand);
        registry.fill(HIST("MC/Rec/hEtaCandRecBg"), candidate.eta());
        registry.fill(HIST("MC/Rec/hEtaCandVsPtCandRecBg"), candidate.eta(), ptCand);
        registry.fill(HIST("MC/Rec/hPhiCandRecBg"), candidate.phi());
        registry.fill(HIST("MC/Rec/hPhiCandVsPtCandRecBg"), candidate.phi(), ptCand);
        registry.fill(HIST("MC/Rec/hMassRecBg"), invMassLcToK0sP(candidate));
        registry.fill(HIST("MC/Rec/hMassVsPtCandRecBg"), invMassLcToK0sP(candidate), ptCand);
        registry.fill(HIST("MC/Rec/hPtBachRecBg"), candidate.ptProng0());
        registry.fill(HIST("MC/Rec/hPtBachVsPtCandRecBg"), candidate.ptProng0(), ptCand);
        registry.fill(HIST("MC/Rec/hPtV0RecBg"), candidate.ptProng1());
        registry.fill(HIST("MC/Rec/hPtV0VsPtCandRecBg"), candidate.ptProng1(), ptCand);
        registry.fill(HIST("MC/Rec/hd0BachRecBg"), candidate.impactParameter0());
        registry.fill(HIST("MC/Rec/hd0BachVsPtCandRecBg"), candidate.impactParameter0(), ptCand);
        registry.fill(HIST("MC/Rec/hd0V0RecBg"), candidate.impactParameter1());
        registry.fill(HIST("MC/Rec/hd0V0VsPtCandRecBg"), candidate.impactParameter1(), ptCand);
        registry.fill(HIST("MC/Rec/hd0V0posRecBg"), candidate.dcapostopv());
        registry.fill(HIST("MC/Rec/hd0V0posVsPtCandRecBg"), candidate.dcapostopv(), ptCand);
        registry.fill(HIST("MC/Rec/hd0V0negRecBg"), candidate.dcanegtopv());
        registry.fill(HIST("MC/Rec/hd0V0negVsPtCandRecBg"), candidate.dcanegtopv(), ptCand);
        registry.fill(HIST("MC/Rec/hPtV0posRecBg"), candidate.ptV0Pos());
        registry.fill(HIST("MC/Rec/hPtV0posVsPtCandRecBg"), candidate.ptV0Pos(), ptCand);
        registry.fill(HIST("MC/Rec/hPtV0negRecBg"), candidate.ptV0Neg());
        registry.fill(HIST("MC/Rec/hPtV0negVsPtCandRecBg"), candidate.ptV0Neg(), ptCand);
        registry.fill(HIST("MC/Rec/hV0CPARecBg"), candidate.v0cosPA());
        registry.fill(HIST("MC/Rec/hV0CPAVsPtCandRecBg"), candidate.v0cosPA(), ptCand);
        registry.fill(HIST("MC/Rec/hV0RadiusRecBg"), candidate.v0radius());
        registry.fill(HIST("MC/Rec/hV0RadiusVsPtCandRecBg"), candidate.v0radius(), ptCand);
        registry.fill(HIST("MC/Rec/hV0DCADaughtersRecBg"), candidate.dcaV0daughters());
        registry.fill(HIST("MC/Rec/hV0DCADaughtersVsPtCandRecBg"), candidate.dcaV0daughters(), ptCand);
        registry.fill(HIST("MC/Rec/hV0MK0ShortRecBg"), candidate.mK0Short());
        registry.fill(HIST("MC/Rec/hV0MK0ShortVsPtCandRecBg"), candidate.mK0Short(), ptCand);
        registry.fill(HIST("MC/Rec/hV0MLambdaRecBg"), candidate.mLambda());
        registry.fill(HIST("MC/Rec/hV0MLambdaVsPtCandRecBg"), candidate.mLambda(), ptCand);
        registry.fill(HIST("MC/Rec/hV0MAntiLambdaRecBg"), candidate.mAntiLambda());
        registry.fill(HIST("MC/Rec/hV0MAntiLambdaVsPtCandRecBg"), candidate.mAntiLambda(), ptCand);
        registry.fill(HIST("MC/Rec/hV0MGammaRecBg"), candidate.mGamma());
        registry.fill(HIST("MC/Rec/hV0MGammaVsPtCandRecBg"), candidate.mGamma(), ptCand);
        registry.fill(HIST("MC/Rec/hCPACandRecBg"), candidate.cpa());
        registry.fill(HIST("MC/Rec/hCPACandVsPtCandRecBg"), candidate.cpa(), ptCand);
        registry.fill(HIST("MC/Rec/hCPAxyCandRecBg"), candidate.cpaXY());
        registry.fill(HIST("MC/Rec/hCPAxyCandVsPtCandRecBg"), candidate.cpaXY(), ptCand);
        registry.fill(HIST("MC/Rec/hDecLengthCandRecBg"), candidate.decayLength());
        registry.fill(HIST("MC/Rec/hDecLengthCandVsPtCandRecBg"), candidate.decayLength(), ptCand);
        registry.fill(HIST("MC/Rec/hDecLengthXYCandRecBg"), candidate.decayLengthXY());
        registry.fill(HIST("MC/Rec/hDecLengthXYCandVsPtCandRecBg"), candidate.decayLengthXY(), ptCand);
        registry.fill(HIST("MC/Rec/hCtCandRecBg"), o2::aod::hf_cand_3prong::ctLc(candidate));
        registry.fill(HIST("MC/Rec/hCtCandVsPtCandRecBg"), o2::aod::hf_cand_3prong::ctLc(candidate), ptCand);
        registry.fill(HIST("MC/Rec/hTPCNSigmaPrBachRecBg"), bach.tpcNSigmaPr());
        registry.fill(HIST("MC/Rec/hPBachVsTPCNSigmaPrBachRecBg"), bach.p(), bach.tpcNSigmaPr());
        if (bach.hasTOF()) {
          registry.fill(HIST("MC/Rec/hTOFNSigmaPrBachRecBg"), bach.tofNSigmaPr());
          registry.fill(HIST("MC/Rec/hPBachVsTOFNSigmaPrBachRecBg"), bach.p(), bach.tofNSigmaPr());
        }
      }
    }
    // MC gen.
    // Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (etaCandMax >= 0. && std::abs(particle.eta()) > etaCandMax) {
        // Printf("MC Gen.: eta rejection: %g", particle.eta());
        continue;
      }
      auto ptCand = particle.pt();
      if (std::abs(particle.flagMcMatchGen()) == 1) {
        registry.fill(HIST("MC/Gen/hPtCandGen"), ptCand);
        registry.fill(HIST("MC/Gen/hEtaCandGen"), particle.eta());
        registry.fill(HIST("MC/Gen/hEtaCandVsPtCandGen"), particle.eta(), ptCand);
        registry.fill(HIST("MC/Gen/hPhiCandGen"), particle.phi());
        registry.fill(HIST("MC/Gen/hPhiCandVsPtCandGen"), particle.phi(), ptCand);
      }
    }
  }

  PROCESS_SWITCH(HfTaskLcToK0sP, processMc, "Process MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskLcToK0sP>(cfgc),
  };
}
