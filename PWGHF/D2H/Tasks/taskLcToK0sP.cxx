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
/// \note based on taskD0.cxx, taskLc.cxx

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// LcToK0sp analysis task
struct HfTaskLcToK0sP {
  Configurable<int> selectionFlagLcToK0sP{"selectionFlagLcToK0sP", 1, "Selection Flag for Lc"};
  Configurable<int> selectionFlagLcbarToK0sP{"selectionFlagLcbarToK0sP", 1, "Selection Flag for Lcbar"};
  Configurable<double> etaCandMax{"etaCandMax", -1., "max. cand. pseudorapidity"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};

  using TracksWPid = soa::Join<aod::TracksWExtra, aod::TracksPidPr>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcToK0sP || aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcbarToK0sP);

  HistogramRegistry registry{"registry"};

  void init(InitContext& context)
  {
    // axes
    AxisSpec const axisBinsPt = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisPt = {300, 0.0f, 30.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisEta = {500, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec const axisPhi = {100, 0.f, 6.3f, "#it{#phi}"};
    AxisSpec const axisMassCand = {600, 1.98f, 2.58f, "inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2})"};
    AxisSpec const axisd0 = {500, -0.5f, 0.5f, "DCAxy (cm)"};
    AxisSpec const axisd0V0Daughters = {1000, -5.0f, 5.0f, "DCAxy (cm)"};
    AxisSpec const axisV0CPA = {500, 0.98f, 1.0001f, "v0 cos pointing angle"};
    AxisSpec const axisV0Radius = {1000, 0.f, 40.f, "V0 radius (cm)"};
    AxisSpec const axisV0DCADaughters = {200, 0.f, 2.f, "DCA (cm)"};
    AxisSpec const axisMassK0Short = {500, 0.4f, 0.6f, "#it{m}(K_{S}^{0}) (GeV/#it{c}^{2})"};
    AxisSpec const axisMassLambda = {500, 1.0f, 1.2f, "#it{m}(#Lambda) (GeV/#it{c}^{2})"};
    AxisSpec const axisMassGamma = {500, 0.0f, 0.4f, "#it{m}(#gamma) (GeV/#it{c}^{2})"};
    AxisSpec const axisCPACand = {110, -1.1f, 1.1f, "candiate cos pointing angle"};
    AxisSpec const axisDecLength = {200, 0.f, 2.0f, "decay length (cm)"};
    AxisSpec const axisProperLifetime = {100, 0.f, 0.2f, "#it{c#tau} (cm)"};
    AxisSpec const axisProperLifetimeV0 = {1000, 0.f, 80.f, "#it{c#tau} (cm)"};
    AxisSpec const axisNSigma = {100, -6.f, 6.f, "n#it{#sigma}_{p}"};
    AxisSpec const axisPidP = {100, 0.f, 10.0f, "#it{p} (GeV/#it{c})"};
    // data
    registry.add("hPtCand", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("hEtaCand", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}});
    registry.add("hEtaCandVsPtCand", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {axisEta, axisBinsPt}});
    registry.add("hPhiCand", "cascade candidates;candidate #it{#phi};entries", {HistType::kTH1F, {axisPhi}});
    registry.add("hPhiCandVsPtCand", "cascade candidates;candidate #it{#phi};p_{T}", {HistType::kTH2F, {axisPhi, axisBinsPt}});
    registry.add("hMass", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassCand}});
    registry.add("hMassVsPtCand", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassCand, axisBinsPt}});
    registry.add("hPtBach", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("hPtBachVsPtCand", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
    registry.add("hPtV0", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("hPtV0VsPtCand", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
    registry.add("hd0Bach", "cascade candidates;bachelor DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0}});
    registry.add("hd0BachVsPtCand", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0, axisBinsPt}});
    registry.add("hd0V0", "cascade candidates;V0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0}});
    registry.add("hd0V0VsPtCand", "cascade candidates;V0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0, axisBinsPt}});
    registry.add("hd0V0pos", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0V0Daughters}});
    registry.add("hd0V0posVsPtCand", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0V0Daughters, axisBinsPt}});
    registry.add("hd0V0neg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0V0Daughters}});
    registry.add("hd0V0negVsPtCand", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0V0Daughters, axisBinsPt}});
    registry.add("hPtV0pos", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("hPtV0posVsPtCand", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
    registry.add("hPtV0neg", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("hPtV0negVsPtCand", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
    registry.add("hV0CPA", "cascade candidates;v0 cosine of pointing angle;entries", {HistType::kTH1F, {axisV0CPA}});
    registry.add("hV0CPAVsPtCand", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {axisV0CPA, axisBinsPt}});
    registry.add("hV0Radius", "cascade candidates;v0 radius (cm);entries", {HistType::kTH1F, {axisV0Radius}});
    registry.add("hV0RadiusVsPtCand", "cascade candidates;v0 radius (cm);p_{T}", {HistType::kTH2F, {axisV0Radius, axisBinsPt}});
    registry.add("hV0DCADaughters", "cascade candidates;v0 dca daughters (cm);entries", {HistType::kTH1F, {axisV0DCADaughters}});
    registry.add("hV0DCADaughtersVsPtCand", "cascade candidates;v0 dca daughters (cm);p_{T}", {HistType::kTH2F, {axisV0DCADaughters, axisBinsPt}});
    registry.add("hV0MK0Short", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassK0Short}});
    registry.add("hV0MK0ShortVsPtCand", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassK0Short, axisBinsPt}});
    registry.add("hV0MLambda", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassLambda}});
    registry.add("hV0MLambdaVsPtCand", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassLambda, axisBinsPt}});
    registry.add("hV0MAntiLambda", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassLambda}});
    registry.add("hV0MAntiLambdaVsPtCand", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassLambda, axisBinsPt}});
    registry.add("hV0MGamma", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassGamma}});
    registry.add("hV0MGammaVsPtCand", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassGamma, axisBinsPt}});
    registry.add("hCtV0K0Short", "cascade candidates;proper lifetime (V0) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetimeV0}});
    registry.add("hCtV0K0ShortVsPtCand", "cascade candidates;proper lifetime (V0) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetimeV0, axisBinsPt}});
    registry.add("hCtV0Lambda", "cascade candidates;proper lifetime (V0) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetimeV0}});
    registry.add("hCtV0LambdaVsPtCand", "cascade candidates;proper lifetime (V0) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetimeV0, axisBinsPt}});
    registry.add("hCPACand", "cascade candidates;cosine pointing angle;entries", {HistType::kTH1F, {axisCPACand}});
    registry.add("hCPACandVsPtCand", "cascade candidates;cosine pointing angle;p_{T}", {HistType::kTH2F, {axisCPACand, axisBinsPt}});
    registry.add("hCPAxyCand", "cascade candidates;cosine pointing angle xy;entries", {HistType::kTH1F, {axisCPACand}});
    registry.add("hCPAxyCandVsPtCand", "cascade candidates;cosine pointing angle xy;p_{T}", {HistType::kTH2F, {axisCPACand, axisBinsPt}});
    registry.add("hDecLengthCand", "cascade candidates;decay length (cm);entries", {HistType::kTH1F, {axisDecLength}});
    registry.add("hDecLengthCandVsPtCand", "cascade candidates;decay length (cm);p_{T}", {HistType::kTH2F, {axisDecLength, axisBinsPt}});
    registry.add("hDecLengthXYCand", "cascade candidates;decay length xy (cm);entries", {HistType::kTH1F, {axisDecLength}});
    registry.add("hDecLengthXYCandVsPtCand", "cascade candidates;decay length xy (cm);p_{T}", {HistType::kTH2F, {axisDecLength, axisBinsPt}});
    registry.add("hCtCand", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetime}});
    registry.add("hCtCandVsPtCand", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetime, axisBinsPt}});
    registry.add("hTPCNSigmaPrBach", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
    registry.add("hPBachVsTPCNSigmaPrBach", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
    registry.add("hTOFNSigmaPrBach", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
    registry.add("hPBachVsTOFNSigmaPrBach", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});

    // add MC histograms
    if (context.mOptions.get<bool>("processMc")) {
      registry.add("MC/Rec/hPtCandRecSig", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtCandRecSigPrompt", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtCandRecSigNonPrompt", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtCandRecBg", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hEtaCandRecSig", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}});
      registry.add("MC/Rec/hEtaCandVsPtCandRecSig", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {axisEta, axisBinsPt}});
      registry.add("MC/Rec/hEtaCandRecBg", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}});
      registry.add("MC/Rec/hEtaCandVsPtCandRecBg", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {axisEta, axisBinsPt}});
      registry.add("MC/Rec/hPhiCandRecSig", "cascade candidates;candidate #it{#phi};entries", {HistType::kTH1F, {axisPhi}});
      registry.add("MC/Rec/hPhiCandVsPtCandRecSig", "cascade candidates;candidate #it{#phi};p_{T}", {HistType::kTH2F, {axisPhi, axisBinsPt}});
      registry.add("MC/Rec/hPhiCandRecBg", "cascade candidates;candidate #it{#phi};entries", {HistType::kTH1F, {axisPhi}});
      registry.add("MC/Rec/hPhiCandVsPtCandRecBg", "cascade candidates;candidate #it{#phi};p_{T}", {HistType::kTH2F, {axisPhi, axisBinsPt}});
      registry.add("MC/Gen/hPtCandGen", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Gen/hPtCandGenPrompt", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Gen/hPtCandGenNonPrompt", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Gen/hEtaCandGen", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}});
      registry.add("MC/Gen/hEtaCandVsPtCandGen", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {axisEta, axisBinsPt}});
      registry.add("MC/Gen/hPhiCandGen", "cascade candidates;candidate #it{#phi};entries", {HistType::kTH1F, {axisPhi}});
      registry.add("MC/Gen/hPhiCandVsPtCandGen", "cascade candidates;candidate #it{#phi};p_{T}", {HistType::kTH2F, {axisPhi, axisBinsPt}});
      registry.add("MC/Rec/hMassRecSig", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH1F, {axisMassCand}});
      registry.add("MC/Rec/hMassVsPtCandRecSig", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassCand, axisBinsPt}});
      registry.add("MC/Rec/hMassRecBg", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH1F, {axisMassCand}});
      registry.add("MC/Rec/hMassVsPtCandRecBg", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassCand, axisBinsPt}});
      registry.add("MC/Rec/hPtBachRecSig", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtBachVsPtCandRecSig", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("MC/Rec/hPtBachRecBg", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtBachVsPtCandRecBg", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("MC/Rec/hPtV0RecSig", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtV0VsPtCandRecSig", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("MC/Rec/hPtV0RecBg", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtV0VsPtCandRecBg", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("MC/Rec/hd0BachRecSig", "cascade candidates;bachelor DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0}});
      registry.add("MC/Rec/hd0BachVsPtCandRecSig", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0, axisBinsPt}});
      registry.add("MC/Rec/hd0BachRecBg", "cascade candidates;bachelor DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0}});
      registry.add("MC/Rec/hd0BachVsPtCandRecBg", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0, axisBinsPt}});
      registry.add("MC/Rec/hd0V0RecSig", "cascade candidates;V0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0}});
      registry.add("MC/Rec/hd0V0VsPtCandRecSig", "cascade candidates;V0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0, axisBinsPt}});
      registry.add("MC/Rec/hd0V0RecBg", "cascade candidates;V0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0}});
      registry.add("MC/Rec/hd0V0VsPtCandRecBg", "cascade candidates;V0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0, axisBinsPt}});
      registry.add("MC/Rec/hd0V0posRecSig", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0V0Daughters}});
      registry.add("MC/Rec/hd0V0posVsPtCandRecSig", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0V0Daughters, axisBinsPt}});
      registry.add("MC/Rec/hd0V0posRecBg", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0V0Daughters}});
      registry.add("MC/Rec/hd0V0posVsPtCandRecBg", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0V0Daughters, axisBinsPt}});
      registry.add("MC/Rec/hd0V0negRecSig", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0V0Daughters}});
      registry.add("MC/Rec/hd0V0negVsPtCandRecSig", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0V0Daughters, axisBinsPt}});
      registry.add("MC/Rec/hd0V0negRecBg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0V0Daughters}});
      registry.add("MC/Rec/hd0V0negVsPtCandRecBg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0V0Daughters, axisBinsPt}});
      registry.add("MC/Rec/hPtV0posRecSig", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtV0posVsPtCandRecSig", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("MC/Rec/hPtV0posRecBg", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtV0posVsPtCandRecBg", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("MC/Rec/hPtV0negRecSig", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtV0negVsPtCandRecSig", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("MC/Rec/hPtV0negRecBg", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtV0negVsPtCandRecBg", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("MC/Rec/hV0CPARecSig", "cascade candidates;v0 cosine of pointing angle;entries", {HistType::kTH1F, {axisV0CPA}});
      registry.add("MC/Rec/hV0CPAVsPtCandRecSig", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {axisV0CPA, axisBinsPt}});
      registry.add("MC/Rec/hV0CPARecBg", "cascade candidates;v0 cosine of pointing angle;entries", {HistType::kTH1F, {axisV0CPA}});
      registry.add("MC/Rec/hV0CPAVsPtCandRecBg", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {axisV0CPA, axisBinsPt}});
      registry.add("MC/Rec/hV0RadiusRecSig", "cascade candidates;v0 radius (cm);entries", {HistType::kTH1F, {axisV0Radius}});
      registry.add("MC/Rec/hV0RadiusVsPtCandRecSig", "cascade candidates;v0 radius (cm);p_{T}", {HistType::kTH2F, {axisV0Radius, axisBinsPt}});
      registry.add("MC/Rec/hV0RadiusRecBg", "cascade candidates;v0 radius (cm);entries", {HistType::kTH1F, {axisV0Radius}});
      registry.add("MC/Rec/hV0RadiusVsPtCandRecBg", "cascade candidates;v0 radius (cm);p_{T}", {HistType::kTH2F, {axisV0Radius, axisBinsPt}});
      registry.add("MC/Rec/hV0DCADaughtersRecSig", "cascade candidates;v0 dca daughters (cm);entries", {HistType::kTH1F, {axisV0DCADaughters}});
      registry.add("MC/Rec/hV0DCADaughtersVsPtCandRecSig", "cascade candidates;v0 dca daughters (cm);p_{T}", {HistType::kTH2F, {axisV0DCADaughters, axisBinsPt}});
      registry.add("MC/Rec/hV0DCADaughtersRecBg", "cascade candidates;v0 dca daughters (cm);entries", {HistType::kTH1F, {axisV0DCADaughters}});
      registry.add("MC/Rec/hV0DCADaughtersVsPtCandRecBg", "cascade candidates;v0 dca daughters (cm);p_{T}", {HistType::kTH2F, {axisV0DCADaughters, axisBinsPt}});
      registry.add("MC/Rec/hV0MK0ShortRecSig", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassK0Short}});
      registry.add("MC/Rec/hV0MK0ShortVsPtCandRecSig", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassK0Short, axisBinsPt}});
      registry.add("MC/Rec/hV0MK0ShortRecBg", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassK0Short}});
      registry.add("MC/Rec/hV0MK0ShortVsPtCandRecBg", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassK0Short, axisBinsPt}});
      registry.add("MC/Rec/hV0MLambdaRecSig", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassLambda}});
      registry.add("MC/Rec/hV0MLambdaVsPtCandRecSig", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassLambda, axisBinsPt}});
      registry.add("MC/Rec/hV0MLambdaRecBg", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassLambda}});
      registry.add("MC/Rec/hV0MLambdaVsPtCandRecBg", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassLambda, axisBinsPt}});
      registry.add("MC/Rec/hV0MAntiLambdaRecSig", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassLambda}});
      registry.add("MC/Rec/hV0MAntiLambdaVsPtCandRecSig", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassLambda, axisBinsPt}});
      registry.add("MC/Rec/hV0MAntiLambdaRecBg", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassLambda}});
      registry.add("MC/Rec/hV0MAntiLambdaVsPtCandRecBg", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassLambda, axisBinsPt}});
      registry.add("MC/Rec/hV0MGammaRecSig", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassGamma}});
      registry.add("MC/Rec/hV0MGammaVsPtCandRecSig", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassGamma, axisBinsPt}});
      registry.add("MC/Rec/hV0MGammaRecBg", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassGamma}});
      registry.add("MC/Rec/hV0MGammaVsPtCandRecBg", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassGamma, axisBinsPt}});
      registry.add("MC/Rec/hCtV0K0ShortRecSig", "cascade candidates;proper lifetime (V0) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetimeV0}});
      registry.add("MC/Rec/hCtV0K0ShortVsPtCandRecSig", "cascade candidates;proper lifetime (V0) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetimeV0, axisBinsPt}});
      registry.add("MC/Rec/hCtV0K0ShortRecBg", "cascade candidates;proper lifetime (V0) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetimeV0}});
      registry.add("MC/Rec/hCtV0K0ShortVsPtCandRecBg", "cascade candidates;proper lifetime (V0) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetimeV0, axisBinsPt}});
      registry.add("MC/Rec/hCtV0LambdaRecSig", "cascade candidates;proper lifetime (V0) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetimeV0}});
      registry.add("MC/Rec/hCtV0LambdaVsPtCandRecSig", "cascade candidates;proper lifetime (V0) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetimeV0, axisBinsPt}});
      registry.add("MC/Rec/hCtV0LambdaRecBg", "cascade candidates;proper lifetime (V0) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetimeV0}});
      registry.add("MC/Rec/hCtV0LambdaVsPtCandRecBg", "cascade candidates;proper lifetime (V0) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetimeV0, axisBinsPt}});
      registry.add("MC/Rec/hCPACandRecSig", "cascade candidates;cosine pointing angle;entries", {HistType::kTH1F, {axisCPACand}});
      registry.add("MC/Rec/hCPACandVsPtCandRecSig", "cascade candidates;cosine pointing angle;p_{T}", {HistType::kTH2F, {axisCPACand, axisBinsPt}});
      registry.add("MC/Rec/hCPACandRecBg", "cascade candidates;cosine pointing angle;entries", {HistType::kTH1F, {axisCPACand}});
      registry.add("MC/Rec/hCPACandVsPtCandRecBg", "cascade candidates;cosine pointing angle;p_{T}", {HistType::kTH2F, {axisCPACand, axisBinsPt}});
      registry.add("MC/Rec/hCPAxyCandRecSig", "cascade candidates;cosine pointing angle xy;entries", {HistType::kTH1F, {axisCPACand}});
      registry.add("MC/Rec/hCPAxyCandVsPtCandRecSig", "cascade candidates;cosine pointing angle xy;p_{T}", {HistType::kTH2F, {axisCPACand, axisBinsPt}});
      registry.add("MC/Rec/hCPAxyCandRecBg", "cascade candidates;cosine pointing angle xy;entries", {HistType::kTH1F, {axisCPACand}});
      registry.add("MC/Rec/hCPAxyCandVsPtCandRecBg", "cascade candidates;cosine pointing angle xy;p_{T}", {HistType::kTH2F, {axisCPACand, axisBinsPt}});
      registry.add("MC/Rec/hDecLengthCandRecSig", "cascade candidates;decay length (cm);entries", {HistType::kTH1F, {axisDecLength}});
      registry.add("MC/Rec/hDecLengthCandVsPtCandRecSig", "cascade candidates;decay length (cm);p_{T}", {HistType::kTH2F, {axisDecLength, axisBinsPt}});
      registry.add("MC/Rec/hDecLengthCandRecBg", "cascade candidates;decay length (cm);entries", {HistType::kTH1F, {axisDecLength}});
      registry.add("MC/Rec/hDecLengthCandVsPtCandRecBg", "cascade candidates;decay length (cm);p_{T}", {HistType::kTH2F, {axisDecLength, axisBinsPt}});
      registry.add("MC/Rec/hDecLengthXYCandRecSig", "cascade candidates;decay length xy (cm);entries", {HistType::kTH1F, {axisDecLength}});
      registry.add("MC/Rec/hDecLengthXYCandVsPtCandRecSig", "cascade candidates;decay length xy (cm);p_{T}", {HistType::kTH2F, {axisDecLength, axisBinsPt}});
      registry.add("MC/Rec/hDecLengthXYCandRecBg", "cascade candidates;decay length xy (cm);entries", {HistType::kTH1F, {axisDecLength}});
      registry.add("MC/Rec/hDecLengthXYCandVsPtCandRecBg", "cascade candidates;decay length xy (cm);p_{T}", {HistType::kTH2F, {axisDecLength, axisBinsPt}});
      registry.add("MC/Rec/hCtCandRecSig", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetime}});
      registry.add("MC/Rec/hCtCandVsPtCandRecSig", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetime, axisBinsPt}});
      registry.add("MC/Rec/hCtCandRecBg", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetime}});
      registry.add("MC/Rec/hCtCandVsPtCandRecBg", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetime, axisBinsPt}});
      registry.add("MC/Rec/hTPCNSigmaPrBachRecSig", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/hPBachVsTPCNSigmaPrBachRecSig", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/hTPCNSigmaPrBachRecBg", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/hPBachVsTPCNSigmaPrBachRecBg", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/hTOFNSigmaPrBachRecSig", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/hPBachVsTOFNSigmaPrBachRecSig", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/hTOFNSigmaPrBachRecBg", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/hPBachVsTOFNSigmaPrBachRecBg", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});
    }
  }

  void process(soa::Filtered<soa::Join<aod::HfCandCascExt, aod::HfSelLcToK0sP>> const& candidates,
               TracksWPid const&)
  {
    for (const auto& candidate : candidates) {
      if (etaCandMax >= 0. && std::abs(candidate.eta()) > etaCandMax) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLc(candidate)) > yCandRecoMax) {
        continue;
      }

      auto ptCand = candidate.pt();
      auto eta = candidate.eta();
      auto phi = candidate.phi();
      auto invMassLcToK0sP = HfHelper::invMassLcToK0sP(candidate);
      auto ptProng0 = candidate.ptProng0();
      auto ptProng1 = candidate.ptProng1();
      auto impactParameter0 = candidate.impactParameter0();
      auto impactParameter1 = candidate.impactParameter1();
      auto dcaPosToPV = candidate.dcapostopv();
      auto dcaNegToPV = candidate.dcanegtopv();
      auto ptV0Pos = candidate.ptV0Pos();
      auto ptV0Neg = candidate.ptV0Neg();
      auto v0CosPA = candidate.v0cosPA();
      auto v0Radius = candidate.v0radius();
      auto dcaV0Daughters = candidate.dcaV0daughters();
      auto mK0Short = candidate.mK0Short();
      auto mLambda = candidate.mLambda();
      auto mAntiLambda = candidate.mAntiLambda();
      auto mGamma = candidate.mGamma();
      auto ctV0K0Short = HfHelper::ctV0K0s(candidate);
      auto ctV0Lambda = HfHelper::ctV0Lambda(candidate);
      auto cpa = candidate.cpa();
      auto cpaXY = candidate.cpaXY();
      auto decayLength = candidate.decayLength();
      auto decayLengthXY = candidate.decayLengthXY();
      auto ctLc = HfHelper::ctLc(candidate);

      registry.fill(HIST("hPtCand"), ptCand);
      registry.fill(HIST("hEtaCand"), eta);
      registry.fill(HIST("hEtaCandVsPtCand"), eta, ptCand);
      registry.fill(HIST("hPhiCand"), phi);
      registry.fill(HIST("hPhiCandVsPtCand"), phi, ptCand);
      registry.fill(HIST("hMass"), invMassLcToK0sP);
      registry.fill(HIST("hMassVsPtCand"), invMassLcToK0sP, ptCand);
      registry.fill(HIST("hPtBach"), ptProng0);
      registry.fill(HIST("hPtBachVsPtCand"), ptProng0, ptCand);
      registry.fill(HIST("hPtV0"), ptProng1);
      registry.fill(HIST("hPtV0VsPtCand"), ptProng1, ptCand);
      registry.fill(HIST("hd0Bach"), impactParameter0);
      registry.fill(HIST("hd0BachVsPtCand"), impactParameter0, ptCand);
      registry.fill(HIST("hd0V0"), impactParameter1);
      registry.fill(HIST("hd0V0VsPtCand"), impactParameter1, ptCand);
      registry.fill(HIST("hd0V0pos"), dcaPosToPV);
      registry.fill(HIST("hd0V0posVsPtCand"), dcaPosToPV, ptCand);
      registry.fill(HIST("hd0V0neg"), dcaNegToPV);
      registry.fill(HIST("hd0V0negVsPtCand"), dcaNegToPV, ptCand);
      registry.fill(HIST("hPtV0pos"), ptV0Pos);
      registry.fill(HIST("hPtV0posVsPtCand"), ptV0Pos, ptCand);
      registry.fill(HIST("hPtV0neg"), ptV0Neg);
      registry.fill(HIST("hPtV0negVsPtCand"), ptV0Neg, ptCand);
      registry.fill(HIST("hV0CPA"), v0CosPA);
      registry.fill(HIST("hV0CPAVsPtCand"), v0CosPA, ptCand);
      registry.fill(HIST("hV0Radius"), v0Radius);
      registry.fill(HIST("hV0RadiusVsPtCand"), v0Radius, ptCand);
      registry.fill(HIST("hV0DCADaughters"), dcaV0Daughters);
      registry.fill(HIST("hV0DCADaughtersVsPtCand"), dcaV0Daughters, ptCand);
      registry.fill(HIST("hV0MK0Short"), mK0Short);
      registry.fill(HIST("hV0MK0ShortVsPtCand"), mK0Short, ptCand);
      registry.fill(HIST("hV0MLambda"), mLambda);
      registry.fill(HIST("hV0MLambdaVsPtCand"), mLambda, ptCand);
      registry.fill(HIST("hV0MAntiLambda"), mAntiLambda);
      registry.fill(HIST("hV0MAntiLambdaVsPtCand"), mAntiLambda, ptCand);
      registry.fill(HIST("hV0MGamma"), mGamma);
      registry.fill(HIST("hV0MGammaVsPtCand"), mGamma, ptCand);
      registry.fill(HIST("hCtV0K0Short"), ctV0K0Short);
      registry.fill(HIST("hCtV0K0ShortVsPtCand"), ctV0K0Short, ptCand);
      registry.fill(HIST("hCtV0Lambda"), ctV0Lambda);
      registry.fill(HIST("hCtV0LambdaVsPtCand"), ctV0Lambda, ptCand);
      registry.fill(HIST("hCPACand"), cpa);
      registry.fill(HIST("hCPACandVsPtCand"), cpa, ptCand);
      registry.fill(HIST("hCPAxyCand"), cpaXY);
      registry.fill(HIST("hCPAxyCandVsPtCand"), cpaXY, ptCand);
      registry.fill(HIST("hDecLengthCand"), decayLength);
      registry.fill(HIST("hDecLengthCandVsPtCand"), decayLength, ptCand);
      registry.fill(HIST("hDecLengthXYCand"), decayLengthXY);
      registry.fill(HIST("hDecLengthXYCandVsPtCand"), decayLengthXY, ptCand);
      registry.fill(HIST("hCtCand"), ctLc);
      registry.fill(HIST("hCtCandVsPtCand"), ctLc, ptCand);

      const auto& bach = candidate.prong0_as<TracksWPid>(); // bachelor track
      auto tpcNSigmaPr = bach.tpcNSigmaPr();
      auto pBach = bach.p();
      registry.fill(HIST("hTPCNSigmaPrBach"), tpcNSigmaPr);
      registry.fill(HIST("hPBachVsTPCNSigmaPrBach"), pBach, tpcNSigmaPr);
      if (bach.hasTOF()) {
        auto tofNSigmaPr = bach.tofNSigmaPr();
        registry.fill(HIST("hTOFNSigmaPrBach"), tofNSigmaPr);
        registry.fill(HIST("hPBachVsTOFNSigmaPrBach"), pBach, tofNSigmaPr);
      }
    }
  }

  void processMc(soa::Filtered<soa::Join<aod::HfCandCascExt, aod::HfSelLcToK0sP, aod::HfCandCascadeMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& mcParticles,
                 aod::TracksWMc const&,
                 TracksWPid const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (etaCandMax >= 0. && std::abs(candidate.eta()) > etaCandMax) {
        continue;
      }

      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLc(candidate)) > yCandRecoMax) {
        continue;
      }

      auto ptCand = candidate.pt();
      auto eta = candidate.eta();
      auto phi = candidate.phi();
      auto invMassLcToK0sP = HfHelper::invMassLcToK0sP(candidate);
      auto ptProng0 = candidate.ptProng0();
      auto ptProng1 = candidate.ptProng1();
      auto impactParameter0 = candidate.impactParameter0();
      auto impactParameter1 = candidate.impactParameter1();
      auto dcaPosToPV = candidate.dcapostopv();
      auto dcaNegToPV = candidate.dcanegtopv();
      auto ptV0Pos = candidate.ptV0Pos();
      auto ptV0Neg = candidate.ptV0Neg();
      auto v0CosPA = candidate.v0cosPA();
      auto v0Radius = candidate.v0radius();
      auto dcaV0Daughters = candidate.dcaV0daughters();
      auto mK0Short = candidate.mK0Short();
      auto mLambda = candidate.mLambda();
      auto mAntiLambda = candidate.mAntiLambda();
      auto mGamma = candidate.mGamma();
      auto ctV0K0Short = HfHelper::ctV0K0s(candidate);
      auto ctV0Lambda = HfHelper::ctV0Lambda(candidate);
      auto cpa = candidate.cpa();
      auto cpaXY = candidate.cpaXY();
      auto decayLength = candidate.decayLength();
      auto decayLengthXY = candidate.decayLengthXY();
      auto ctLc = HfHelper::ctLc(candidate);

      const auto& bach = candidate.prong0_as<TracksWPid>(); // bachelor track
      auto tpcNSigmaPr = bach.tpcNSigmaPr();
      auto pBach = bach.p();

      if (std::abs(candidate.flagMcMatchRec()) == 1) {
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/Rec/hPtCandRecSigPrompt"), ptCand);
        } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("MC/Rec/hPtCandRecSigNonPrompt"), ptCand);
        }
        registry.fill(HIST("MC/Rec/hPtCandRecSig"), ptCand);
        registry.fill(HIST("MC/Rec/hEtaCandRecSig"), eta);
        registry.fill(HIST("MC/Rec/hEtaCandVsPtCandRecSig"), eta, ptCand);
        registry.fill(HIST("MC/Rec/hPhiCandRecSig"), phi);
        registry.fill(HIST("MC/Rec/hPhiCandVsPtCandRecSig"), phi, ptCand);
        registry.fill(HIST("MC/Rec/hMassRecSig"), invMassLcToK0sP);
        registry.fill(HIST("MC/Rec/hMassVsPtCandRecSig"), invMassLcToK0sP, ptCand);
        registry.fill(HIST("MC/Rec/hPtBachRecSig"), ptProng0);
        registry.fill(HIST("MC/Rec/hPtBachVsPtCandRecSig"), ptProng0, ptCand);
        registry.fill(HIST("MC/Rec/hPtV0RecSig"), ptProng1);
        registry.fill(HIST("MC/Rec/hPtV0VsPtCandRecSig"), ptProng1, ptCand);
        registry.fill(HIST("MC/Rec/hd0BachRecSig"), impactParameter0);
        registry.fill(HIST("MC/Rec/hd0BachVsPtCandRecSig"), impactParameter0, ptCand);
        registry.fill(HIST("MC/Rec/hd0V0RecSig"), impactParameter1);
        registry.fill(HIST("MC/Rec/hd0V0VsPtCandRecSig"), impactParameter1, ptCand);
        registry.fill(HIST("MC/Rec/hd0V0posRecSig"), dcaPosToPV);
        registry.fill(HIST("MC/Rec/hd0V0posVsPtCandRecSig"), dcaPosToPV, ptCand);
        registry.fill(HIST("MC/Rec/hd0V0negRecSig"), dcaNegToPV);
        registry.fill(HIST("MC/Rec/hd0V0negVsPtCandRecSig"), dcaNegToPV, ptCand);
        registry.fill(HIST("MC/Rec/hPtV0posRecSig"), ptV0Pos);
        registry.fill(HIST("MC/Rec/hPtV0posVsPtCandRecSig"), ptV0Pos, ptCand);
        registry.fill(HIST("MC/Rec/hPtV0negRecSig"), ptV0Neg);
        registry.fill(HIST("MC/Rec/hPtV0negVsPtCandRecSig"), ptV0Neg, ptCand);
        registry.fill(HIST("MC/Rec/hV0CPARecSig"), v0CosPA);
        registry.fill(HIST("MC/Rec/hV0CPAVsPtCandRecSig"), v0CosPA, ptCand);
        registry.fill(HIST("MC/Rec/hV0RadiusRecSig"), v0Radius);
        registry.fill(HIST("MC/Rec/hV0RadiusVsPtCandRecSig"), v0Radius, ptCand);
        registry.fill(HIST("MC/Rec/hV0DCADaughtersRecSig"), dcaV0Daughters);
        registry.fill(HIST("MC/Rec/hV0DCADaughtersVsPtCandRecSig"), dcaV0Daughters, ptCand);
        registry.fill(HIST("MC/Rec/hV0MK0ShortRecSig"), mK0Short);
        registry.fill(HIST("MC/Rec/hV0MK0ShortVsPtCandRecSig"), mK0Short, ptCand);
        registry.fill(HIST("MC/Rec/hV0MLambdaRecSig"), mLambda);
        registry.fill(HIST("MC/Rec/hV0MLambdaVsPtCandRecSig"), mLambda, ptCand);
        registry.fill(HIST("MC/Rec/hV0MAntiLambdaRecSig"), mAntiLambda);
        registry.fill(HIST("MC/Rec/hV0MAntiLambdaVsPtCandRecSig"), mAntiLambda, ptCand);
        registry.fill(HIST("MC/Rec/hV0MGammaRecSig"), mGamma);
        registry.fill(HIST("MC/Rec/hV0MGammaVsPtCandRecSig"), mGamma, ptCand);
        registry.fill(HIST("MC/Rec/hCtV0K0ShortRecSig"), ctV0K0Short);
        registry.fill(HIST("MC/Rec/hCtV0K0ShortVsPtCandRecSig"), ctV0K0Short, ptCand);
        registry.fill(HIST("MC/Rec/hCtV0LambdaRecSig"), ctV0Lambda);
        registry.fill(HIST("MC/Rec/hCtV0LambdaVsPtCandRecSig"), ctV0Lambda, ptCand);
        registry.fill(HIST("MC/Rec/hCPACandRecSig"), cpa);
        registry.fill(HIST("MC/Rec/hCPACandVsPtCandRecSig"), cpa, ptCand);
        registry.fill(HIST("MC/Rec/hCPAxyCandRecSig"), cpaXY);
        registry.fill(HIST("MC/Rec/hCPAxyCandVsPtCandRecSig"), cpaXY, ptCand);
        registry.fill(HIST("MC/Rec/hDecLengthCandRecSig"), decayLength);
        registry.fill(HIST("MC/Rec/hDecLengthCandVsPtCandRecSig"), decayLength, ptCand);
        registry.fill(HIST("MC/Rec/hDecLengthXYCandRecSig"), decayLengthXY);
        registry.fill(HIST("MC/Rec/hDecLengthXYCandVsPtCandRecSig"), decayLengthXY, ptCand);
        registry.fill(HIST("MC/Rec/hCtCandRecSig"), ctLc);
        registry.fill(HIST("MC/Rec/hCtCandVsPtCandRecSig"), ctLc, ptCand);
        registry.fill(HIST("MC/Rec/hTPCNSigmaPrBachRecSig"), tpcNSigmaPr);
        registry.fill(HIST("MC/Rec/hPBachVsTPCNSigmaPrBachRecSig"), pBach, tpcNSigmaPr);
        if (bach.hasTOF()) {
          auto tofNSigmaPr = bach.tofNSigmaPr();
          registry.fill(HIST("MC/Rec/hTOFNSigmaPrBachRecSig"), tofNSigmaPr);
          registry.fill(HIST("MC/Rec/hPBachVsTOFNSigmaPrBachRecSig"), pBach, tofNSigmaPr);
        }
      } else {
        registry.fill(HIST("MC/Rec/hPtCandRecBg"), ptCand);
        registry.fill(HIST("MC/Rec/hEtaCandRecBg"), eta);
        registry.fill(HIST("MC/Rec/hEtaCandVsPtCandRecBg"), eta, ptCand);
        registry.fill(HIST("MC/Rec/hPhiCandRecBg"), phi);
        registry.fill(HIST("MC/Rec/hPhiCandVsPtCandRecBg"), phi, ptCand);
        registry.fill(HIST("MC/Rec/hMassRecBg"), invMassLcToK0sP);
        registry.fill(HIST("MC/Rec/hMassVsPtCandRecBg"), invMassLcToK0sP, ptCand);
        registry.fill(HIST("MC/Rec/hPtBachRecBg"), ptProng0);
        registry.fill(HIST("MC/Rec/hPtBachVsPtCandRecBg"), ptProng0, ptCand);
        registry.fill(HIST("MC/Rec/hPtV0RecBg"), ptProng1);
        registry.fill(HIST("MC/Rec/hPtV0VsPtCandRecBg"), ptProng1, ptCand);
        registry.fill(HIST("MC/Rec/hd0BachRecBg"), impactParameter0);
        registry.fill(HIST("MC/Rec/hd0BachVsPtCandRecBg"), impactParameter0, ptCand);
        registry.fill(HIST("MC/Rec/hd0V0RecBg"), impactParameter1);
        registry.fill(HIST("MC/Rec/hd0V0VsPtCandRecBg"), impactParameter1, ptCand);
        registry.fill(HIST("MC/Rec/hd0V0posRecBg"), dcaPosToPV);
        registry.fill(HIST("MC/Rec/hd0V0posVsPtCandRecBg"), dcaPosToPV, ptCand);
        registry.fill(HIST("MC/Rec/hd0V0negRecBg"), dcaNegToPV);
        registry.fill(HIST("MC/Rec/hd0V0negVsPtCandRecBg"), dcaNegToPV, ptCand);
        registry.fill(HIST("MC/Rec/hPtV0posRecBg"), ptV0Pos);
        registry.fill(HIST("MC/Rec/hPtV0posVsPtCandRecBg"), ptV0Pos, ptCand);
        registry.fill(HIST("MC/Rec/hPtV0negRecBg"), ptV0Neg);
        registry.fill(HIST("MC/Rec/hPtV0negVsPtCandRecBg"), ptV0Neg, ptCand);
        registry.fill(HIST("MC/Rec/hV0CPARecBg"), v0CosPA);
        registry.fill(HIST("MC/Rec/hV0CPAVsPtCandRecBg"), v0CosPA, ptCand);
        registry.fill(HIST("MC/Rec/hV0RadiusRecBg"), v0Radius);
        registry.fill(HIST("MC/Rec/hV0RadiusVsPtCandRecBg"), v0Radius, ptCand);
        registry.fill(HIST("MC/Rec/hV0DCADaughtersRecBg"), dcaV0Daughters);
        registry.fill(HIST("MC/Rec/hV0DCADaughtersVsPtCandRecBg"), dcaV0Daughters, ptCand);
        registry.fill(HIST("MC/Rec/hV0MK0ShortRecBg"), mK0Short);
        registry.fill(HIST("MC/Rec/hV0MK0ShortVsPtCandRecBg"), mK0Short, ptCand);
        registry.fill(HIST("MC/Rec/hV0MLambdaRecBg"), mLambda);
        registry.fill(HIST("MC/Rec/hV0MLambdaVsPtCandRecBg"), mLambda, ptCand);
        registry.fill(HIST("MC/Rec/hV0MAntiLambdaRecBg"), mAntiLambda);
        registry.fill(HIST("MC/Rec/hV0MAntiLambdaVsPtCandRecBg"), mAntiLambda, ptCand);
        registry.fill(HIST("MC/Rec/hV0MGammaRecBg"), mGamma);
        registry.fill(HIST("MC/Rec/hV0MGammaVsPtCandRecBg"), mGamma, ptCand);
        registry.fill(HIST("MC/Rec/hCtV0K0ShortRecBg"), ctV0K0Short);
        registry.fill(HIST("MC/Rec/hCtV0K0ShortVsPtCandRecBg"), ctV0K0Short, ptCand);
        registry.fill(HIST("MC/Rec/hCtV0LambdaRecBg"), ctV0Lambda);
        registry.fill(HIST("MC/Rec/hCtV0LambdaVsPtCandRecBg"), ctV0Lambda, ptCand);
        registry.fill(HIST("MC/Rec/hCPACandRecBg"), cpa);
        registry.fill(HIST("MC/Rec/hCPACandVsPtCandRecBg"), cpa, ptCand);
        registry.fill(HIST("MC/Rec/hCPAxyCandRecBg"), cpaXY);
        registry.fill(HIST("MC/Rec/hCPAxyCandVsPtCandRecBg"), cpaXY, ptCand);
        registry.fill(HIST("MC/Rec/hDecLengthCandRecBg"), decayLength);
        registry.fill(HIST("MC/Rec/hDecLengthCandVsPtCandRecBg"), decayLength, ptCand);
        registry.fill(HIST("MC/Rec/hDecLengthXYCandRecBg"), decayLengthXY);
        registry.fill(HIST("MC/Rec/hDecLengthXYCandVsPtCandRecBg"), decayLengthXY, ptCand);
        registry.fill(HIST("MC/Rec/hCtCandRecBg"), ctLc);
        registry.fill(HIST("MC/Rec/hCtCandVsPtCandRecBg"), ctLc, ptCand);
        registry.fill(HIST("MC/Rec/hTPCNSigmaPrBachRecBg"), tpcNSigmaPr);
        registry.fill(HIST("MC/Rec/hPBachVsTPCNSigmaPrBachRecBg"), pBach, tpcNSigmaPr);
        if (bach.hasTOF()) {
          auto tofNSigmaPr = bach.tofNSigmaPr();
          registry.fill(HIST("MC/Rec/hTOFNSigmaPrBachRecBg"), tofNSigmaPr);
          registry.fill(HIST("MC/Rec/hPBachVsTOFNSigmaPrBachRecBg"), pBach, tofNSigmaPr);
        }
      }
    }
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (etaCandMax >= 0. && std::abs(particle.eta()) > etaCandMax) {
        continue;
      }

      if (std::abs(particle.flagMcMatchGen()) == 1) {

        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus);
        if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
          continue;
        }
        auto ptCand = particle.pt();
        auto eta = particle.eta();
        auto phi = particle.phi();
        registry.fill(HIST("MC/Gen/hPtCandGen"), ptCand);
        registry.fill(HIST("MC/Gen/hEtaCandGen"), eta);
        registry.fill(HIST("MC/Gen/hEtaCandVsPtCandGen"), eta, ptCand);
        registry.fill(HIST("MC/Gen/hPhiCandGen"), phi);
        registry.fill(HIST("MC/Gen/hPhiCandVsPtCandGen"), phi, ptCand);

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/Gen/hPtCandGenPrompt"), ptCand);
        } else if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("MC/Gen/hPtCandGenNonPrompt"), ptCand);
        }
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
