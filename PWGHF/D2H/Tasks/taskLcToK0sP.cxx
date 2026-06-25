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
///         Xufei Xue, <xufei.xue@cern.ch>, CUG
///
/// \note based on taskD0.cxx, taskLc.cxx

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;

/// LcToK0sp analysis task
struct HfTaskLcToK0sP {
  Configurable<bool> fillTHn{"fillTHn", true, "fill THnSparse"};
  Configurable<int> selectionFlagLcToK0sP{"selectionFlagLcToK0sP", 1, "Selection Flag for Lc"};
  Configurable<double> etaCandMax{"etaCandMax", -1., "max. cand. pseudorapidity"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};

  SliceCache cache;

  using FilteredCandLcToPK0SWSelFlag = soa::Filtered<soa::Join<aod::HfCandCascade, aod::HfSelLcToK0sP>>;
  using FilteredCandLcToPK0SWSelFlagAndMl = soa::Filtered<soa::Join<aod::HfCandCascade, aod::HfSelLcToK0sP, aod::HfMlLcToK0sP>>;
  using FilteredCandLcToPK0SWSelFlagAndMc = soa::Filtered<soa::Join<aod::HfCandCascade, aod::HfSelLcToK0sP, aod::HfCandCascadeMcRec>>;
  using FilteredCandLcToPK0SWSelFlagAndMcAndMl = soa::Filtered<soa::Join<aod::HfCandCascade, aod::HfSelLcToK0sP, aod::HfCandCascadeMcRec, aod::HfMlLcToK0sP>>;
  using McParticlesCascadeMatched = soa::Join<aod::McParticles, aod::HfCandCascadeMcGen>;

  using TracksWPid = soa::Join<aod::TracksWExtra, aod::TracksPidPr>;

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsMc = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsMcWithFT0C = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsMcWithFT0M = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Ms>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcToK0sP;

  Preslice<aod::HfCandCascade> lcToPK0SPerCollision = aod::hf_cand::collisionId;
  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mcparticle::mcCollisionId;

  // THnSparse configurable axes
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {300, 2.0, 2.6}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {360, 0., 36.}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreBkg{"thnConfigAxisBdtScoreBkg", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisBdtScorePrompt{"thnConfigAxisBdtScorePrompt", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreNonPrompt{"thnConfigAxisBdtScoreNonPrompt", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {3000, 0., 300.}, ""};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, ""};
  ConfigurableAxis thnConfigAxisCentrality{"thnConfigAxisCentrality", {100, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1., 1.}, "candidate rapidity"};

  HistogramRegistry registry{"registry"};

  // ML class indices
  enum MlClasses : int {
    MlClassBackground = 0,
    MlClassPrompt,
    MlClassNonPrompt,
    NumberOfMlClasses
  };

  static constexpr int DecayChannelLcToK0sP = 1;

  // Names of folders and suffixes for MC signal histograms
  constexpr static std::string_view SignalFolders[] = {"signal", "prompt", "nonprompt"};
  constexpr static std::string_view SignalSuffixes[] = {"", "Prompt", "NonPrompt"};

  enum SignalClass : int {
    Signal = 0,
    Prompt,
    NonPrompt
  };

  void init(InitContext&)
  {
    // Check that only one process function is enabled
    std::array<bool, 12> processFlags{
      doprocessDataStd, doprocessDataStdWithFT0C, doprocessDataStdWithFT0M,
      doprocessDataWithMl, doprocessDataWithMlWithFT0C, doprocessDataWithMlWithFT0M,
      doprocessMcStd, doprocessMcStdWithFT0C, doprocessMcStdWithFT0M,
      doprocessMcWithMl, doprocessMcWithMlWithFT0C, doprocessMcWithMlWithFT0M};
    if ((std::accumulate(processFlags.begin(), processFlags.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }

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
    AxisSpec const axisCPACand = {110, -1.1f, 1.1f, "candidate cos pointing angle"};
    AxisSpec const axisDecLength = {200, 0.f, 2.0f, "decay length (cm)"};
    AxisSpec const axisProperLifetime = {100, 0.f, 0.2f, "#it{c#tau} (cm)"};
    AxisSpec const axisProperLifetimeV0 = {1000, 0.f, 80.f, "#it{c#tau} (cm)"};
    AxisSpec const axisNSigma = {100, -6.f, 6.f, "n#it{#sigma}_{p}"};
    AxisSpec const axisPidP = {100, 0.f, 10.0f, "#it{p} (GeV/#it{c})"};
    // data - helper lambda: adds 1D histogram and its VsPtCand 2D variant
    auto addHistosData = [&](const std::string& name, const std::string& title, const AxisSpec& axis) {
      registry.add(("Data/h" + name).c_str(), ("cascade candidates;" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
      registry.add(("Data/h" + name + "VsPtCand").c_str(), ("cascade candidates;" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
    };

    addHistosData("PtCand", "candidate #it{p}_{T} (GeV/#it{c})", axisPt);
    addHistosData("EtaCand", "candidate #it{#eta}", axisEta);
    addHistosData("PhiCand", "candidate #it{#phi}", axisPhi);
    addHistosData("Mass", "inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2})", axisMassCand);
    addHistosData("PtBach", "bachelor #it{p}_{T} (GeV/#it{c})", axisPt);
    addHistosData("PtV0", "v0 #it{p}_{T} (GeV/#it{c})", axisPt);
    addHistosData("d0Bach", "bachelor DCAxy to prim. vertex (cm)", axisd0);
    addHistosData("d0V0", "V0 DCAxy to prim. vertex (cm)", axisd0);
    addHistosData("d0V0pos", "pos daugh v0 DCAxy to prim. vertex (cm)", axisd0V0Daughters);
    addHistosData("d0V0neg", "neg daugh v0 DCAxy to prim. vertex (cm)", axisd0V0Daughters);
    addHistosData("PtV0pos", "pos daugh v0 #it{p}_{T} (GeV/#it{c})", axisPt);
    addHistosData("PtV0neg", "neg daugh v0 #it{p}_{T} (GeV/#it{c})", axisPt);
    addHistosData("V0CPA", "v0 cosine of pointing angle", axisV0CPA);
    addHistosData("V0Radius", "V0 radius (cm)", axisV0Radius);
    addHistosData("V0DCADaughters", "v0 dca daughters (cm)", axisV0DCADaughters);
    addHistosData("V0MK0Short", "v0 mass K0s (GeV/#it{c}^{2})", axisMassK0Short);
    addHistosData("V0MLambda", "v0 mass #Lambda (GeV/#it{c}^{2})", axisMassLambda);
    addHistosData("V0MAntiLambda", "v0 mass Anti#Lambda (GeV/#it{c}^{2})", axisMassLambda);
    addHistosData("V0MGamma", "v0 mass #gamma (GeV/#it{c}^{2})", axisMassGamma);
    addHistosData("CtV0K0Short", "proper lifetime (V0) * #it{c} (cm)", axisProperLifetimeV0);
    addHistosData("CtV0Lambda", "proper lifetime (V0) * #it{c} (cm)", axisProperLifetimeV0);
    addHistosData("CPACand", "cosine pointing angle", axisCPACand);
    addHistosData("CPAxyCand", "cosine pointing angle xy", axisCPACand);
    addHistosData("DecLengthCand", "decay length (cm)", axisDecLength);
    addHistosData("DecLengthXYCand", "decay length xy (cm)", axisDecLength);
    addHistosData("CtCand", "proper lifetime (#Lambda_{c}) * #it{c} (cm)", axisProperLifetime);

    // PID histograms (non-standard axis pattern)
    registry.add("Data/hTPCNSigmaPrBach", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
    registry.add("Data/hPBachVsTPCNSigmaPrBach", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
    registry.add("Data/hTOFNSigmaPrBach", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
    registry.add("Data/hPBachVsTOFNSigmaPrBach", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});

    // add MC histograms
    bool const isMc = doprocessMcStd || doprocessMcStdWithFT0C || doprocessMcStdWithFT0M ||
                      doprocessMcWithMl || doprocessMcWithMlWithFT0C || doprocessMcWithMlWithFT0M;
    if (isMc) {
      // MC Rec helper: signal/prompt/nonprompt + background
      auto addHistosMcRec = [&](const std::string& name, const std::string& title, const AxisSpec& axis) {
        registry.add(("MC/Rec/signal/h" + name + "RecSig").c_str(), ("cascade candidates (matched);" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Rec/signal/h" + name + "VsPtCandRecSig").c_str(), ("cascade candidates (matched);" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
        registry.add(("MC/Rec/prompt/h" + name + "RecSigPrompt").c_str(), ("cascade candidates (matched, prompt);" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Rec/prompt/h" + name + "VsPtCandRecSigPrompt").c_str(), ("cascade candidates (matched, prompt);" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
        registry.add(("MC/Rec/nonprompt/h" + name + "RecSigNonPrompt").c_str(), ("cascade candidates (matched, non-prompt);" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Rec/nonprompt/h" + name + "VsPtCandRecSigNonPrompt").c_str(), ("cascade candidates (matched, non-prompt);" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
        registry.add(("MC/Rec/background/h" + name + "RecBg").c_str(), ("cascade candidates (unmatched);" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Rec/background/h" + name + "VsPtCandRecBg").c_str(), ("cascade candidates (unmatched);" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
      };

      // MC Gen helper: signal/prompt/nonprompt
      auto addHistosMcGen = [&](const std::string& name, const std::string& title, const AxisSpec& axis) {
        registry.add(("MC/Gen/signal/h" + name + "Gen").c_str(), ("MC particles (matched);" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Gen/signal/h" + name + "VsPtCandGen").c_str(), ("MC particles (matched);" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
        registry.add(("MC/Gen/prompt/h" + name + "GenPrompt").c_str(), ("MC particles (matched, prompt);" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Gen/prompt/h" + name + "VsPtCandGenPrompt").c_str(), ("MC particles (matched, prompt);" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
        registry.add(("MC/Gen/nonprompt/h" + name + "GenNonPrompt").c_str(), ("MC particles (matched, non-prompt);" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Gen/nonprompt/h" + name + "VsPtCandGenNonPrompt").c_str(), ("MC particles (matched, non-prompt);" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
      };

      // MC Rec: PtCand has extra signal/prompt/nonprompt + background variants
      registry.add("MC/Rec/signal/hPtCandRecSig", "cascade candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/prompt/hPtCandRecSigPrompt", "cascade candidates (matched, prompt);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/nonprompt/hPtCandRecSigNonPrompt", "cascade candidates (matched, non-prompt);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/background/hPtCandRecBg", "cascade candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});

      addHistosMcRec("EtaCand", "candidate #it{#eta}", axisEta);
      addHistosMcRec("PhiCand", "candidate #it{#phi}", axisPhi);
      addHistosMcRec("Mass", "inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2})", axisMassCand);
      addHistosMcRec("PtBach", "bachelor #it{p}_{T} (GeV/#it{c})", axisPt);
      addHistosMcRec("PtV0", "v0 #it{p}_{T} (GeV/#it{c})", axisPt);
      addHistosMcRec("d0Bach", "bachelor DCAxy to prim. vertex (cm)", axisd0);
      addHistosMcRec("d0V0", "V0 DCAxy to prim. vertex (cm)", axisd0);
      addHistosMcRec("d0V0pos", "pos daugh v0 DCAxy to prim. vertex (cm)", axisd0V0Daughters);
      addHistosMcRec("d0V0neg", "neg daugh v0 DCAxy to prim. vertex (cm)", axisd0V0Daughters);
      addHistosMcRec("PtV0pos", "pos daugh v0 #it{p}_{T} (GeV/#it{c})", axisPt);
      addHistosMcRec("PtV0neg", "neg daugh v0 #it{p}_{T} (GeV/#it{c})", axisPt);
      addHistosMcRec("V0CPA", "v0 cosine of pointing angle", axisV0CPA);
      addHistosMcRec("V0Radius", "V0 radius (cm)", axisV0Radius);
      addHistosMcRec("V0DCADaughters", "v0 dca daughters (cm)", axisV0DCADaughters);
      addHistosMcRec("V0MK0Short", "v0 mass K0s (GeV/#it{c}^{2})", axisMassK0Short);
      addHistosMcRec("V0MLambda", "v0 mass #Lambda (GeV/#it{c}^{2})", axisMassLambda);
      addHistosMcRec("V0MAntiLambda", "v0 mass Anti#Lambda (GeV/#it{c}^{2})", axisMassLambda);
      addHistosMcRec("V0MGamma", "v0 mass #gamma (GeV/#it{c}^{2})", axisMassGamma);
      addHistosMcRec("CtV0K0Short", "proper lifetime (V0) * #it{c} (cm)", axisProperLifetimeV0);
      addHistosMcRec("CtV0Lambda", "proper lifetime (V0) * #it{c} (cm)", axisProperLifetimeV0);
      addHistosMcRec("CPACand", "cosine pointing angle", axisCPACand);
      addHistosMcRec("CPAxyCand", "cosine pointing angle xy", axisCPACand);
      addHistosMcRec("DecLengthCand", "decay length (cm)", axisDecLength);
      addHistosMcRec("DecLengthXYCand", "decay length xy (cm)", axisDecLength);
      addHistosMcRec("CtCand", "proper lifetime (#Lambda_{c}) * #it{c} (cm)", axisProperLifetime);

      // MC Rec PID
      registry.add("MC/Rec/signal/hTPCNSigmaPrBachRecSig", "cascade candidates (matched);n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/signal/hPBachVsTPCNSigmaPrBachRecSig", "cascade candidates (matched);#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/signal/hTOFNSigmaPrBachRecSig", "cascade candidates (matched);n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/signal/hPBachVsTOFNSigmaPrBachRecSig", "cascade candidates (matched);#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/prompt/hTPCNSigmaPrBachRecSigPrompt", "cascade candidates (matched, prompt);n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/prompt/hPBachVsTPCNSigmaPrBachRecSigPrompt", "cascade candidates (matched, prompt);#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/prompt/hTOFNSigmaPrBachRecSigPrompt", "cascade candidates (matched, prompt);n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/prompt/hPBachVsTOFNSigmaPrBachRecSigPrompt", "cascade candidates (matched, prompt);#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/nonprompt/hTPCNSigmaPrBachRecSigNonPrompt", "cascade candidates (matched, non-prompt);n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/nonprompt/hPBachVsTPCNSigmaPrBachRecSigNonPrompt", "cascade candidates (matched, non-prompt);#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/nonprompt/hTOFNSigmaPrBachRecSigNonPrompt", "cascade candidates (matched, non-prompt);n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/nonprompt/hPBachVsTOFNSigmaPrBachRecSigNonPrompt", "cascade candidates (matched, non-prompt);#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/background/hTPCNSigmaPrBachRecBg", "cascade candidates (unmatched);n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/background/hPBachVsTPCNSigmaPrBachRecBg", "cascade candidates (unmatched);#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/background/hTOFNSigmaPrBachRecBg", "cascade candidates (unmatched);n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/background/hPBachVsTOFNSigmaPrBachRecBg", "cascade candidates (unmatched);#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});

      // MC Gen: PtCand with signal/prompt/nonprompt variants
      registry.add("MC/Gen/signal/hPtCandGen", "MC particles (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Gen/prompt/hPtCandGenPrompt", "MC particles (matched, prompt);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Gen/nonprompt/hPtCandGenNonPrompt", "MC particles (matched, non-prompt);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});

      addHistosMcGen("EtaCand", "candidate #it{#eta}", axisEta);
      addHistosMcGen("PhiCand", "candidate #it{#phi}", axisPhi);
    }

    // THnSparse for ML analysis
    if (fillTHn) {
      const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec thnAxisBdtScoreBkg{thnConfigAxisBdtScoreBkg, "BDT bkg score"};
      const AxisSpec thnAxisBdtScorePrompt{thnConfigAxisBdtScorePrompt, "BDT prompt score"};
      const AxisSpec thnAxisBdtScoreNonPrompt{thnConfigAxisBdtScoreNonPrompt, "BDT non-prompt score"};
      bool const isFT0C = doprocessDataStdWithFT0C || doprocessDataWithMlWithFT0C || doprocessMcStdWithFT0C || doprocessMcWithMlWithFT0C;
      bool const isFT0M = doprocessDataStdWithFT0M || doprocessDataWithMlWithFT0M || doprocessMcStdWithFT0M || doprocessMcWithMlWithFT0M;
      std::string centLabel = "centrality";
      if (isFT0C) {
        centLabel = "centrality (FT0C)";
      } else if (isFT0M) {
        centLabel = "centrality (FT0M)";
      }
      const AxisSpec thnAxisCentrality{thnConfigAxisCentrality, centLabel};
      const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of PV contributors"};
      const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
      const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "origin"};
      const AxisSpec thnAxisY{thnConfigAxisY, "rapidity"};

      // Data with ML: {mass, pt, centrality, bkg, prompt, non-prompt, numPvContr}
      if (doprocessDataWithMl || doprocessDataWithMlWithFT0C || doprocessDataWithMlWithFT0M) {
        registry.add("hnLcK0sPDataWithBdt", "THn for Lc->K0sP data with BDT", HistType::kTHnSparseF,
                     {thnAxisMass, thnAxisPt, thnAxisCentrality, thnAxisBdtScoreBkg, thnAxisBdtScorePrompt, thnAxisBdtScoreNonPrompt, thnAxisNumPvContr});
      }

      // MC Rec with ML: {mass, pt, centrality, bkg, prompt, non-prompt, numPvContr, origin}
      if (doprocessMcWithMl || doprocessMcWithMlWithFT0C || doprocessMcWithMlWithFT0M) {
        registry.add("hnLcK0sPRecMcWithBdt", "THn for Lc->K0sP MC rec with BDT", HistType::kTHnSparseF,
                     {thnAxisMass, thnAxisPt, thnAxisCentrality, thnAxisBdtScoreBkg, thnAxisBdtScorePrompt, thnAxisBdtScoreNonPrompt, thnAxisNumPvContr, thnAxisOrigin});
      }

      // MC Gen: {pt, centrality, rapidity, numPvContr, ptB, origin}
      if (isMc) {
        registry.add("hnLcK0sPGenMc", "THn for Lc->K0sP MC gen", HistType::kTHnSparseF,
                     {thnAxisPt, thnAxisCentrality, thnAxisY, thnAxisNumPvContr, thnAxisPtB, thnAxisOrigin});
      }
    }
  }

  /// Evaluate centrality/multiplicity percentile
  /// \param collision is collision
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll>
  float evaluateCentralityColl(const Coll& collision)
  {
    return o2::hf_centrality::getCentralityColl<Coll>(collision);
  }

  /// Bundled candidate variables with conversion constructor
  struct HfCandVars {
    double ptCand{}, eta{}, phi{}, invMassLcToK0sP{};
    double ptProng0{}, ptProng1{};
    double impactParameter0{}, impactParameter1{};
    double dcaPosToPV{}, dcaNegToPV{};
    double ptV0Pos{}, ptV0Neg{};
    double v0CosPA{}, v0Radius{}, dcaV0Daughters{};
    double mK0Short{}, mLambda{}, mAntiLambda{}, mGamma{};
    double ctV0K0Short{}, ctV0Lambda{};
    double cpa{}, cpaXY{};
    double decayLength{}, decayLengthXY{}, ctLc{};

    template <typename CandType>
    explicit HfCandVars(CandType const& candidate)
      : ptCand(candidate.pt()),
        eta(candidate.eta()),
        phi(candidate.phi()),
        invMassLcToK0sP(HfHelper::invMassLcToK0sP(candidate)),
        ptProng0(candidate.ptProng0()),
        ptProng1(candidate.ptProng1()),
        impactParameter0(candidate.impactParameter0()),
        impactParameter1(candidate.impactParameter1()),
        dcaPosToPV(candidate.dcapostopv()),
        dcaNegToPV(candidate.dcanegtopv()),
        ptV0Pos(candidate.ptV0Pos()),
        ptV0Neg(candidate.ptV0Neg()),
        v0CosPA(candidate.v0cosPA()),
        v0Radius(candidate.v0radius()),
        dcaV0Daughters(candidate.dcaV0daughters()),
        mK0Short(candidate.mK0Short()),
        mLambda(candidate.mLambda()),
        mAntiLambda(candidate.mAntiLambda()),
        mGamma(candidate.mGamma()),
        ctV0K0Short(HfHelper::ctV0K0s(candidate)),
        ctV0Lambda(HfHelper::ctV0Lambda(candidate)),
        cpa(candidate.cpa()),
        cpaXY(candidate.cpaXY()),
        decayLength(candidate.decayLength()),
        decayLengthXY(candidate.decayLengthXY()),
        ctLc(HfHelper::ctLc(candidate))
    {
    }
  };

  /// Helper function to fill candidate histograms
  /// \param candidate is the candidate
  template <typename CandType>
  void fillCandHistograms(CandType const& candidate)
  {
    HfCandVars cand(candidate);

    registry.fill(HIST("Data/hPtCand"), cand.ptCand);
    registry.fill(HIST("Data/hEtaCand"), cand.eta);
    registry.fill(HIST("Data/hEtaCandVsPtCand"), cand.eta, cand.ptCand);
    registry.fill(HIST("Data/hPhiCand"), cand.phi);
    registry.fill(HIST("Data/hPhiCandVsPtCand"), cand.phi, cand.ptCand);
    registry.fill(HIST("Data/hMass"), cand.invMassLcToK0sP);
    registry.fill(HIST("Data/hMassVsPtCand"), cand.invMassLcToK0sP, cand.ptCand);
    registry.fill(HIST("Data/hPtBach"), cand.ptProng0);
    registry.fill(HIST("Data/hPtBachVsPtCand"), cand.ptProng0, cand.ptCand);
    registry.fill(HIST("Data/hPtV0"), cand.ptProng1);
    registry.fill(HIST("Data/hPtV0VsPtCand"), cand.ptProng1, cand.ptCand);
    registry.fill(HIST("Data/hd0Bach"), cand.impactParameter0);
    registry.fill(HIST("Data/hd0BachVsPtCand"), cand.impactParameter0, cand.ptCand);
    registry.fill(HIST("Data/hd0V0"), cand.impactParameter1);
    registry.fill(HIST("Data/hd0V0VsPtCand"), cand.impactParameter1, cand.ptCand);
    registry.fill(HIST("Data/hd0V0pos"), cand.dcaPosToPV);
    registry.fill(HIST("Data/hd0V0posVsPtCand"), cand.dcaPosToPV, cand.ptCand);
    registry.fill(HIST("Data/hd0V0neg"), cand.dcaNegToPV);
    registry.fill(HIST("Data/hd0V0negVsPtCand"), cand.dcaNegToPV, cand.ptCand);
    registry.fill(HIST("Data/hPtV0pos"), cand.ptV0Pos);
    registry.fill(HIST("Data/hPtV0posVsPtCand"), cand.ptV0Pos, cand.ptCand);
    registry.fill(HIST("Data/hPtV0neg"), cand.ptV0Neg);
    registry.fill(HIST("Data/hPtV0negVsPtCand"), cand.ptV0Neg, cand.ptCand);
    registry.fill(HIST("Data/hV0CPA"), cand.v0CosPA);
    registry.fill(HIST("Data/hV0CPAVsPtCand"), cand.v0CosPA, cand.ptCand);
    registry.fill(HIST("Data/hV0Radius"), cand.v0Radius);
    registry.fill(HIST("Data/hV0RadiusVsPtCand"), cand.v0Radius, cand.ptCand);
    registry.fill(HIST("Data/hV0DCADaughters"), cand.dcaV0Daughters);
    registry.fill(HIST("Data/hV0DCADaughtersVsPtCand"), cand.dcaV0Daughters, cand.ptCand);
    registry.fill(HIST("Data/hV0MK0Short"), cand.mK0Short);
    registry.fill(HIST("Data/hV0MK0ShortVsPtCand"), cand.mK0Short, cand.ptCand);
    registry.fill(HIST("Data/hV0MLambda"), cand.mLambda);
    registry.fill(HIST("Data/hV0MLambdaVsPtCand"), cand.mLambda, cand.ptCand);
    registry.fill(HIST("Data/hV0MAntiLambda"), cand.mAntiLambda);
    registry.fill(HIST("Data/hV0MAntiLambdaVsPtCand"), cand.mAntiLambda, cand.ptCand);
    registry.fill(HIST("Data/hV0MGamma"), cand.mGamma);
    registry.fill(HIST("Data/hV0MGammaVsPtCand"), cand.mGamma, cand.ptCand);
    registry.fill(HIST("Data/hCtV0K0Short"), cand.ctV0K0Short);
    registry.fill(HIST("Data/hCtV0K0ShortVsPtCand"), cand.ctV0K0Short, cand.ptCand);
    registry.fill(HIST("Data/hCtV0Lambda"), cand.ctV0Lambda);
    registry.fill(HIST("Data/hCtV0LambdaVsPtCand"), cand.ctV0Lambda, cand.ptCand);
    registry.fill(HIST("Data/hCPACand"), cand.cpa);
    registry.fill(HIST("Data/hCPACandVsPtCand"), cand.cpa, cand.ptCand);
    registry.fill(HIST("Data/hCPAxyCand"), cand.cpaXY);
    registry.fill(HIST("Data/hCPAxyCandVsPtCand"), cand.cpaXY, cand.ptCand);
    registry.fill(HIST("Data/hDecLengthCand"), cand.decayLength);
    registry.fill(HIST("Data/hDecLengthCandVsPtCand"), cand.decayLength, cand.ptCand);
    registry.fill(HIST("Data/hDecLengthXYCand"), cand.decayLengthXY);
    registry.fill(HIST("Data/hDecLengthXYCandVsPtCand"), cand.decayLengthXY, cand.ptCand);
    registry.fill(HIST("Data/hCtCand"), cand.ctLc);
    registry.fill(HIST("Data/hCtCandVsPtCand"), cand.ctLc, cand.ptCand);
  }

  /// Data processing template
  /// \tparam FillMl switch to fill ML histograms
  template <bool FillMl, typename CollType, typename CandType>
  void runAnalysisPerCollisionData(CollType const& collisions, CandType const& candidates)
  {
    for (const auto& collision : collisions) {
      fillHistosData<FillMl>(collision, candidates);
    }
  }

  /// Helper function to fill data histograms
  /// \tparam FillMl switch to fill ML histograms
  template <bool FillMl, typename CollType, typename CandType>
  void fillHistosData(CollType const& collision, CandType const& candidates)
  {
    const auto thisCollId = collision.globalIndex();
    const auto& groupedCandidates = candidates.sliceBy(lcToPK0SPerCollision, thisCollId);

    for (const auto& candidate : groupedCandidates) {
      if (etaCandMax >= 0. && std::abs(candidate.eta()) > etaCandMax) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLc(candidate)) > yCandRecoMax) {
        continue;
      }
      if (candidate.isSelLcToK0sP() < selectionFlagLcToK0sP) {
        continue;
      }

      fillCandHistograms(candidate);

      // PID histograms
      const auto& bach = candidate.template prong0_as<TracksWPid>();
      auto tpcNSigmaPr = bach.tpcNSigmaPr();
      auto pBach = bach.p();
      registry.fill(HIST("Data/hTPCNSigmaPrBach"), tpcNSigmaPr);
      registry.fill(HIST("Data/hPBachVsTPCNSigmaPrBach"), pBach, tpcNSigmaPr);
      if (bach.hasTOF()) {
        registry.fill(HIST("Data/hTOFNSigmaPrBach"), bach.tofNSigmaPr());
        registry.fill(HIST("Data/hPBachVsTOFNSigmaPrBach"), pBach, bach.tofNSigmaPr());
      }

      if (fillTHn && FillMl) {
        if constexpr (CandType::template contains<aod::HfMlLcToK0sP>()) {
          if (candidate.mlProbLcToK0sP().size() == NumberOfMlClasses) {
            float ml0 = candidate.mlProbLcToK0sP()[MlClassBackground];
            float ml1 = candidate.mlProbLcToK0sP()[MlClassPrompt];
            float ml2 = candidate.mlProbLcToK0sP()[MlClassNonPrompt];
            if (ml0 >= 0.f && ml1 >= 0.f && ml2 >= 0.f) {
              float cent = evaluateCentralityColl(collision);
              int numPvContr = collision.numContrib();
              registry.get<THnSparse>(HIST("hnLcK0sPDataWithBdt"))->Fill(HfHelper::invMassLcToK0sP(candidate), candidate.pt(), cent, ml0, ml1, ml2, static_cast<double>(numPvContr));
            }
          }
        }
      }
    }
  }

  /// Fill MC reconstructed signal histograms
  /// \tparam SignalType signal class (Signal, Prompt, NonPrompt)
  template <int SignalType, typename CandidateType>
  void fillHistogramsRecSig(CandidateType const& candidate)
  {
    HfCandVars cand(candidate);

    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.eta);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaCandVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.eta, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.phi);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiCandVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.phi, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hMassRecSig") + HIST(SignalSuffixes[SignalType]), cand.invMassLcToK0sP);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hMassVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.invMassLcToK0sP, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtBachRecSig") + HIST(SignalSuffixes[SignalType]), cand.ptProng0);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtBachVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ptProng0, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtV0RecSig") + HIST(SignalSuffixes[SignalType]), cand.ptProng1);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtV0VsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ptProng1, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0BachRecSig") + HIST(SignalSuffixes[SignalType]), cand.impactParameter0);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0BachVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.impactParameter0, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0V0RecSig") + HIST(SignalSuffixes[SignalType]), cand.impactParameter1);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0V0VsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.impactParameter1, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0V0posRecSig") + HIST(SignalSuffixes[SignalType]), cand.dcaPosToPV);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0V0posVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.dcaPosToPV, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0V0negRecSig") + HIST(SignalSuffixes[SignalType]), cand.dcaNegToPV);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0V0negVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.dcaNegToPV, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtV0posRecSig") + HIST(SignalSuffixes[SignalType]), cand.ptV0Pos);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtV0posVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ptV0Pos, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtV0negRecSig") + HIST(SignalSuffixes[SignalType]), cand.ptV0Neg);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtV0negVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ptV0Neg, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0CPARecSig") + HIST(SignalSuffixes[SignalType]), cand.v0CosPA);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0CPAVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.v0CosPA, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0RadiusRecSig") + HIST(SignalSuffixes[SignalType]), cand.v0Radius);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0RadiusVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.v0Radius, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0DCADaughtersRecSig") + HIST(SignalSuffixes[SignalType]), cand.dcaV0Daughters);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0DCADaughtersVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.dcaV0Daughters, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0MK0ShortRecSig") + HIST(SignalSuffixes[SignalType]), cand.mK0Short);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0MK0ShortVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.mK0Short, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0MLambdaRecSig") + HIST(SignalSuffixes[SignalType]), cand.mLambda);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0MLambdaVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.mLambda, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0MAntiLambdaRecSig") + HIST(SignalSuffixes[SignalType]), cand.mAntiLambda);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0MAntiLambdaVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.mAntiLambda, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0MGammaRecSig") + HIST(SignalSuffixes[SignalType]), cand.mGamma);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hV0MGammaVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.mGamma, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCtV0K0ShortRecSig") + HIST(SignalSuffixes[SignalType]), cand.ctV0K0Short);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCtV0K0ShortVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ctV0K0Short, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCtV0LambdaRecSig") + HIST(SignalSuffixes[SignalType]), cand.ctV0Lambda);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCtV0LambdaVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ctV0Lambda, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCPACandRecSig") + HIST(SignalSuffixes[SignalType]), cand.cpa);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCPACandVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.cpa, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCPAxyCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.cpaXY);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCPAxyCandVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.cpaXY, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.decayLength);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthCandVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.decayLength, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthXYCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.decayLengthXY);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthXYCandVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.decayLengthXY, cand.ptCand);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ctLc);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCtCandVsPtCandRecSig") + HIST(SignalSuffixes[SignalType]), cand.ctLc, cand.ptCand);

    // PID
    const auto& bach = candidate.template prong0_as<TracksWPid>();
    auto tpcNSigmaPr = bach.tpcNSigmaPr();
    auto pBach = bach.p();
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hTPCNSigmaPrBachRecSig") + HIST(SignalSuffixes[SignalType]), tpcNSigmaPr);
    registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPBachVsTPCNSigmaPrBachRecSig") + HIST(SignalSuffixes[SignalType]), pBach, tpcNSigmaPr);
    if (bach.hasTOF()) {
      registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hTOFNSigmaPrBachRecSig") + HIST(SignalSuffixes[SignalType]), bach.tofNSigmaPr());
      registry.fill(HIST("MC/Rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPBachVsTOFNSigmaPrBachRecSig") + HIST(SignalSuffixes[SignalType]), pBach, bach.tofNSigmaPr());
    }
  }

  /// MC processing template
  /// \tparam FillMl switch to fill ML histograms
  template <bool FillMl, typename CollType, typename CandType, typename CandMcGen>
  void runAnalysisPerCollisionMc(CollType const& collisions,
                                 CandType const& candidates,
                                 CandMcGen const& mcParticles)
  {
    for (const auto& collision : collisions) {
      fillHistosMcRec<FillMl>(collision, candidates, mcParticles);
    }
    fillHistosMcGen(mcParticles, collisions);
  }

  /// Helper function to fill MC reconstructed histograms
  /// \tparam FillMl switch to fill ML histograms
  template <bool FillMl, typename CollType, typename CandType, typename CandMcGen>
  void fillHistosMcRec(CollType const& collision, CandType const& candidates, CandMcGen const&)
  {
    const auto thisCollId = collision.globalIndex();
    const auto& groupedCandidates = candidates.sliceBy(lcToPK0SPerCollision, thisCollId);

    for (const auto& candidate : groupedCandidates) {
      if (etaCandMax >= 0. && std::abs(candidate.eta()) > etaCandMax) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLc(candidate)) > yCandRecoMax) {
        continue;
      }

      if (std::abs(candidate.flagMcMatchRec()) == DecayChannelLcToK0sP) {
        fillCandHistograms(candidate);

        // MC reconstructed signal
        fillHistogramsRecSig<Signal>(candidate);

        // reconstructed signal prompt / nonprompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          fillHistogramsRecSig<Prompt>(candidate);
        } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
          fillHistogramsRecSig<NonPrompt>(candidate);
        }

        if (fillTHn && FillMl) {
          float mass = HfHelper::invMassLcToK0sP(candidate);
          float pt = candidate.pt();
          std::array<float, 3> mlScores{-1.f, -1.f, -1.f};
          if constexpr (CandType::template contains<aod::HfMlLcToK0sP>()) {
            if (candidate.mlProbLcToK0sP().size() == NumberOfMlClasses) {
              mlScores = {candidate.mlProbLcToK0sP()[MlClassBackground],
                          candidate.mlProbLcToK0sP()[MlClassPrompt],
                          candidate.mlProbLcToK0sP()[MlClassNonPrompt]};
            }
          }
          float cent = evaluateCentralityColl(collision);
          int numPvContr = collision.numContrib();
          int8_t origin = candidate.originMcRec();
          registry.get<THnSparse>(HIST("hnLcK0sPRecMcWithBdt"))->Fill(mass, pt, cent, mlScores[0], mlScores[1], mlScores[2], static_cast<double>(numPvContr), origin);
        }
      } else {
        // MC reconstructed background
        HfCandVars cand(candidate);

        registry.fill(HIST("MC/Rec/background/hPtCandRecBg"), cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hEtaCandRecBg"), cand.eta);
        registry.fill(HIST("MC/Rec/background/hEtaCandVsPtCandRecBg"), cand.eta, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hPhiCandRecBg"), cand.phi);
        registry.fill(HIST("MC/Rec/background/hPhiCandVsPtCandRecBg"), cand.phi, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hMassRecBg"), cand.invMassLcToK0sP);
        registry.fill(HIST("MC/Rec/background/hMassVsPtCandRecBg"), cand.invMassLcToK0sP, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hPtBachRecBg"), cand.ptProng0);
        registry.fill(HIST("MC/Rec/background/hPtBachVsPtCandRecBg"), cand.ptProng0, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hPtV0RecBg"), cand.ptProng1);
        registry.fill(HIST("MC/Rec/background/hPtV0VsPtCandRecBg"), cand.ptProng1, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hd0BachRecBg"), cand.impactParameter0);
        registry.fill(HIST("MC/Rec/background/hd0BachVsPtCandRecBg"), cand.impactParameter0, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hd0V0RecBg"), cand.impactParameter1);
        registry.fill(HIST("MC/Rec/background/hd0V0VsPtCandRecBg"), cand.impactParameter1, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hd0V0posRecBg"), cand.dcaPosToPV);
        registry.fill(HIST("MC/Rec/background/hd0V0posVsPtCandRecBg"), cand.dcaPosToPV, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hd0V0negRecBg"), cand.dcaNegToPV);
        registry.fill(HIST("MC/Rec/background/hd0V0negVsPtCandRecBg"), cand.dcaNegToPV, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hPtV0posRecBg"), cand.ptV0Pos);
        registry.fill(HIST("MC/Rec/background/hPtV0posVsPtCandRecBg"), cand.ptV0Pos, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hPtV0negRecBg"), cand.ptV0Neg);
        registry.fill(HIST("MC/Rec/background/hPtV0negVsPtCandRecBg"), cand.ptV0Neg, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hV0CPARecBg"), cand.v0CosPA);
        registry.fill(HIST("MC/Rec/background/hV0CPAVsPtCandRecBg"), cand.v0CosPA, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hV0RadiusRecBg"), cand.v0Radius);
        registry.fill(HIST("MC/Rec/background/hV0RadiusVsPtCandRecBg"), cand.v0Radius, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hV0DCADaughtersRecBg"), cand.dcaV0Daughters);
        registry.fill(HIST("MC/Rec/background/hV0DCADaughtersVsPtCandRecBg"), cand.dcaV0Daughters, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hV0MK0ShortRecBg"), cand.mK0Short);
        registry.fill(HIST("MC/Rec/background/hV0MK0ShortVsPtCandRecBg"), cand.mK0Short, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hV0MLambdaRecBg"), cand.mLambda);
        registry.fill(HIST("MC/Rec/background/hV0MLambdaVsPtCandRecBg"), cand.mLambda, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hV0MAntiLambdaRecBg"), cand.mAntiLambda);
        registry.fill(HIST("MC/Rec/background/hV0MAntiLambdaVsPtCandRecBg"), cand.mAntiLambda, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hV0MGammaRecBg"), cand.mGamma);
        registry.fill(HIST("MC/Rec/background/hV0MGammaVsPtCandRecBg"), cand.mGamma, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hCtV0K0ShortRecBg"), cand.ctV0K0Short);
        registry.fill(HIST("MC/Rec/background/hCtV0K0ShortVsPtCandRecBg"), cand.ctV0K0Short, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hCtV0LambdaRecBg"), cand.ctV0Lambda);
        registry.fill(HIST("MC/Rec/background/hCtV0LambdaVsPtCandRecBg"), cand.ctV0Lambda, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hCPACandRecBg"), cand.cpa);
        registry.fill(HIST("MC/Rec/background/hCPACandVsPtCandRecBg"), cand.cpa, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hCPAxyCandRecBg"), cand.cpaXY);
        registry.fill(HIST("MC/Rec/background/hCPAxyCandVsPtCandRecBg"), cand.cpaXY, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hDecLengthCandRecBg"), cand.decayLength);
        registry.fill(HIST("MC/Rec/background/hDecLengthCandVsPtCandRecBg"), cand.decayLength, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hDecLengthXYCandRecBg"), cand.decayLengthXY);
        registry.fill(HIST("MC/Rec/background/hDecLengthXYCandVsPtCandRecBg"), cand.decayLengthXY, cand.ptCand);
        registry.fill(HIST("MC/Rec/background/hCtCandRecBg"), cand.ctLc);
        registry.fill(HIST("MC/Rec/background/hCtCandVsPtCandRecBg"), cand.ctLc, cand.ptCand);

        // PID background
        const auto& bach = candidate.template prong0_as<TracksWPid>();
        auto tpcNSigmaPr = bach.tpcNSigmaPr();
        auto pBach = bach.p();
        registry.fill(HIST("MC/Rec/background/hTPCNSigmaPrBachRecBg"), tpcNSigmaPr);
        registry.fill(HIST("MC/Rec/background/hPBachVsTPCNSigmaPrBachRecBg"), pBach, tpcNSigmaPr);
        if (bach.hasTOF()) {
          registry.fill(HIST("MC/Rec/background/hTOFNSigmaPrBachRecBg"), bach.tofNSigmaPr());
          registry.fill(HIST("MC/Rec/background/hPBachVsTOFNSigmaPrBachRecBg"), pBach, bach.tofNSigmaPr());
        }
      }
    }
  }

  /// Fill MC generated signal histograms
  /// \tparam SignalType signal class (Signal, Prompt, NonPrompt)
  template <int SignalType, typename ParticleType>
  void fillHistogramsGenSig(ParticleType const& particle)
  {
    registry.fill(HIST("MC/Gen/") + HIST(SignalFolders[SignalType]) + HIST("/hPtCandGen") + HIST(SignalSuffixes[SignalType]), particle.pt());
    registry.fill(HIST("MC/Gen/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaCandGen") + HIST(SignalSuffixes[SignalType]), particle.eta());
    registry.fill(HIST("MC/Gen/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiCandGen") + HIST(SignalSuffixes[SignalType]), particle.phi());
  }

  /// Helper function to fill MC generated histograms
  template <typename CandMcGen, typename Coll>
  void fillHistosMcGen(CandMcGen const& mcParticles, Coll const& recoCollisions)
  {
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == DecayChannelLcToK0sP) {
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus);
        if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
          continue;
        }
        fillHistogramsGenSig<Signal>(particle);
        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          fillHistogramsGenSig<Prompt>(particle);
        } else if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          fillHistogramsGenSig<NonPrompt>(particle);
        }

        if (fillTHn) {
          int8_t origin = particle.originMcGen();
          float ptB = (origin == RecoDecay::OriginType::NonPrompt) ? mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt() : -1.;
          const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
          int numPvContr = 0;
          for (const auto& recCol : recoCollsPerMcColl) {
            numPvContr = recCol.numContrib() > numPvContr ? recCol.numContrib() : numPvContr;
          }
          float cent = o2::hf_centrality::getCentralityGenColl(recoCollsPerMcColl);
          registry.get<THnSparse>(HIST("hnLcK0sPGenMc"))->Fill(particle.pt(), cent, yGen, static_cast<double>(numPvContr), ptB, origin);
        }
      }
    }
  }

  void processDataStd(Collisions const& collisions,
                      FilteredCandLcToPK0SWSelFlag const& candidates,
                      TracksWPid const&)
  {
    runAnalysisPerCollisionData<false>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processDataStd, "Process Data", false);

  void processMcStd(CollisionsMc const& collisions,
                    FilteredCandLcToPK0SWSelFlagAndMc const& candidates,
                    McParticlesCascadeMatched const& mcParticles,
                    aod::TracksWMc const&,
                    aod::McCollisions const&,
                    TracksWPid const&)
  {
    runAnalysisPerCollisionMc<false>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcStd, "Process MC data", false);

  // Data with ML
  void processDataWithMl(Collisions const& collisions,
                         FilteredCandLcToPK0SWSelFlagAndMl const& candidates,
                         TracksWPid const&)
  {
    runAnalysisPerCollisionData<true>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processDataWithMl, "Process Data with ML", false);

  // MC Rec and Gen with ML
  void processMcWithMl(CollisionsMc const& collisions,
                       FilteredCandLcToPK0SWSelFlagAndMcAndMl const& candidates,
                       McParticlesCascadeMatched const& mcParticles,
                       aod::TracksWMc const&,
                       aod::McCollisions const&,
                       TracksWPid const&)
  {
    runAnalysisPerCollisionMc<true>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcWithMl, "Process MC with ML", false);

  // Data with FT0C centrality
  void processDataStdWithFT0C(CollisionsWithFT0C const& collisions,
                              FilteredCandLcToPK0SWSelFlag const& candidates,
                              TracksWPid const&)
  {
    runAnalysisPerCollisionData<false>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processDataStdWithFT0C, "Process Data with FT0C centrality", false);

  // Data with FT0M centrality
  void processDataStdWithFT0M(CollisionsWithFT0M const& collisions,
                              FilteredCandLcToPK0SWSelFlag const& candidates,
                              TracksWPid const&)
  {
    runAnalysisPerCollisionData<false>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processDataStdWithFT0M, "Process Data with FT0M centrality", false);

  // Data with ML + FT0C centrality
  void processDataWithMlWithFT0C(CollisionsWithFT0C const& collisions,
                                 FilteredCandLcToPK0SWSelFlagAndMl const& candidates,
                                 TracksWPid const&)
  {
    runAnalysisPerCollisionData<true>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processDataWithMlWithFT0C, "Process Data with ML and FT0C centrality", false);

  // Data with ML + FT0M centrality
  void processDataWithMlWithFT0M(CollisionsWithFT0M const& collisions,
                                 FilteredCandLcToPK0SWSelFlagAndMl const& candidates,
                                 TracksWPid const&)
  {
    runAnalysisPerCollisionData<true>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processDataWithMlWithFT0M, "Process Data with ML and FT0M centrality", false);

  // MC Std with FT0C centrality
  void processMcStdWithFT0C(CollisionsMcWithFT0C const& collisions,
                            FilteredCandLcToPK0SWSelFlagAndMc const& candidates,
                            McParticlesCascadeMatched const& mcParticles,
                            aod::TracksWMc const&,
                            aod::McCollisions const&,
                            TracksWPid const&)
  {
    runAnalysisPerCollisionMc<false>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcStdWithFT0C, "Process MC with FT0C centrality", false);

  // MC Std with FT0M centrality
  void processMcStdWithFT0M(CollisionsMcWithFT0M const& collisions,
                            FilteredCandLcToPK0SWSelFlagAndMc const& candidates,
                            McParticlesCascadeMatched const& mcParticles,
                            aod::TracksWMc const&,
                            aod::McCollisions const&,
                            TracksWPid const&)
  {
    runAnalysisPerCollisionMc<false>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcStdWithFT0M, "Process MC with FT0M centrality", false);

  // MC with ML + FT0C centrality
  void processMcWithMlWithFT0C(CollisionsMcWithFT0C const& collisions,
                               FilteredCandLcToPK0SWSelFlagAndMcAndMl const& candidates,
                               McParticlesCascadeMatched const& mcParticles,
                               aod::TracksWMc const&,
                               aod::McCollisions const&,
                               TracksWPid const&)
  {
    runAnalysisPerCollisionMc<true>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcWithMlWithFT0C, "Process MC with ML and FT0C centrality", false);

  // MC with ML + FT0M centrality
  void processMcWithMlWithFT0M(CollisionsMcWithFT0M const& collisions,
                               FilteredCandLcToPK0SWSelFlagAndMcAndMl const& candidates,
                               McParticlesCascadeMatched const& mcParticles,
                               aod::TracksWMc const&,
                               aod::McCollisions const&,
                               TracksWPid const&)
  {
    runAnalysisPerCollisionMc<true>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcWithMlWithFT0M, "Process MC with ML and FT0M centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskLcToK0sP>(cfgc)};
}
