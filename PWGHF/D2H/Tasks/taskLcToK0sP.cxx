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

#include <cstdlib>
#include <numeric>
#include <string>
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
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {2, -0.5, 1.5}, ""};
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

  static constexpr int KDecayChannelLcToK0sP = 1;

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
    AxisSpec const axisCPACand = {110, -1.1f, 1.1f, "candiate cos pointing angle"};
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
      // MC helper lambdas
      auto addHistosMcRec = [&](const std::string& name, const std::string& title, const AxisSpec& axis) {
        registry.add(("MC/Rec/h" + name + "RecSig").c_str(), ("cascade candidates;" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Rec/h" + name + "VsPtCandRecSig").c_str(), ("cascade candidates;" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
        registry.add(("MC/Rec/h" + name + "RecBg").c_str(), ("cascade candidates;" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Rec/h" + name + "VsPtCandRecBg").c_str(), ("cascade candidates;" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
      };
      auto addHistosMcGen = [&](const std::string& name, const std::string& title, const AxisSpec& axis) {
        registry.add(("MC/Gen/h" + name + "Gen").c_str(), ("cascade candidates;" + title + ";entries").c_str(), {HistType::kTH1F, {axis}});
        registry.add(("MC/Gen/h" + name + "VsPtCandGen").c_str(), ("cascade candidates;" + title + ";p_{T}").c_str(), {HistType::kTH2F, {axis, axisBinsPt}});
      };

      // MC Rec: PtCand has extra Prompt/NonPrompt variants
      registry.add("MC/Rec/hPtCandRecSig", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtCandRecSigPrompt", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtCandRecSigNonPrompt", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Rec/hPtCandRecBg", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});

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
      registry.add("MC/Rec/hTPCNSigmaPrBachRecSig", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/hPBachVsTPCNSigmaPrBachRecSig", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/hTPCNSigmaPrBachRecBg", "cascade candidates;n#it{#sigma}_{p} TPC;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/hPBachVsTPCNSigmaPrBachRecBg", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/hTOFNSigmaPrBachRecSig", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/hPBachVsTOFNSigmaPrBachRecSig", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});
      registry.add("MC/Rec/hTOFNSigmaPrBachRecBg", "cascade candidates;n#it{#sigma}_{p} TOF;entries", {HistType::kTH1F, {axisNSigma}});
      registry.add("MC/Rec/hPBachVsTOFNSigmaPrBachRecBg", "cascade candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigma}});

      // MC Gen
      registry.add("MC/Gen/hPtCandGen", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Gen/hPtCandGenPrompt", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("MC/Gen/hPtCandGenNonPrompt", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});

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
      const AxisSpec thnAxisCentrality{thnConfigAxisCentrality, "centrality (FT0C)"};
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
      if (doprocessMcWithMl || doprocessMcWithMlWithFT0C || doprocessMcWithMlWithFT0M) {
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

  /// Helper function to fill candidate histograms
  /// \param candidate is the candidate
  template <typename CandType>
  void fillCandHistograms(CandType const& candidate)
  {
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

    registry.fill(HIST("Data/hPtCand"), ptCand);
    registry.fill(HIST("Data/hEtaCand"), eta);
    registry.fill(HIST("Data/hEtaCandVsPtCand"), eta, ptCand);
    registry.fill(HIST("Data/hPhiCand"), phi);
    registry.fill(HIST("Data/hPhiCandVsPtCand"), phi, ptCand);
    registry.fill(HIST("Data/hMass"), invMassLcToK0sP);
    registry.fill(HIST("Data/hMassVsPtCand"), invMassLcToK0sP, ptCand);
    registry.fill(HIST("Data/hPtBach"), ptProng0);
    registry.fill(HIST("Data/hPtBachVsPtCand"), ptProng0, ptCand);
    registry.fill(HIST("Data/hPtV0"), ptProng1);
    registry.fill(HIST("Data/hPtV0VsPtCand"), ptProng1, ptCand);
    registry.fill(HIST("Data/hd0Bach"), impactParameter0);
    registry.fill(HIST("Data/hd0BachVsPtCand"), impactParameter0, ptCand);
    registry.fill(HIST("Data/hd0V0"), impactParameter1);
    registry.fill(HIST("Data/hd0V0VsPtCand"), impactParameter1, ptCand);
    registry.fill(HIST("Data/hd0V0pos"), dcaPosToPV);
    registry.fill(HIST("Data/hd0V0posVsPtCand"), dcaPosToPV, ptCand);
    registry.fill(HIST("Data/hd0V0neg"), dcaNegToPV);
    registry.fill(HIST("Data/hd0V0negVsPtCand"), dcaNegToPV, ptCand);
    registry.fill(HIST("Data/hPtV0pos"), ptV0Pos);
    registry.fill(HIST("Data/hPtV0posVsPtCand"), ptV0Pos, ptCand);
    registry.fill(HIST("Data/hPtV0neg"), ptV0Neg);
    registry.fill(HIST("Data/hPtV0negVsPtCand"), ptV0Neg, ptCand);
    registry.fill(HIST("Data/hV0CPA"), v0CosPA);
    registry.fill(HIST("Data/hV0CPAVsPtCand"), v0CosPA, ptCand);
    registry.fill(HIST("Data/hV0Radius"), v0Radius);
    registry.fill(HIST("Data/hV0RadiusVsPtCand"), v0Radius, ptCand);
    registry.fill(HIST("Data/hV0DCADaughters"), dcaV0Daughters);
    registry.fill(HIST("Data/hV0DCADaughtersVsPtCand"), dcaV0Daughters, ptCand);
    registry.fill(HIST("Data/hV0MK0Short"), mK0Short);
    registry.fill(HIST("Data/hV0MK0ShortVsPtCand"), mK0Short, ptCand);
    registry.fill(HIST("Data/hV0MLambda"), mLambda);
    registry.fill(HIST("Data/hV0MLambdaVsPtCand"), mLambda, ptCand);
    registry.fill(HIST("Data/hV0MAntiLambda"), mAntiLambda);
    registry.fill(HIST("Data/hV0MAntiLambdaVsPtCand"), mAntiLambda, ptCand);
    registry.fill(HIST("Data/hV0MGamma"), mGamma);
    registry.fill(HIST("Data/hV0MGammaVsPtCand"), mGamma, ptCand);
    registry.fill(HIST("Data/hCtV0K0Short"), ctV0K0Short);
    registry.fill(HIST("Data/hCtV0K0ShortVsPtCand"), ctV0K0Short, ptCand);
    registry.fill(HIST("Data/hCtV0Lambda"), ctV0Lambda);
    registry.fill(HIST("Data/hCtV0LambdaVsPtCand"), ctV0Lambda, ptCand);
    registry.fill(HIST("Data/hCPACand"), cpa);
    registry.fill(HIST("Data/hCPACandVsPtCand"), cpa, ptCand);
    registry.fill(HIST("Data/hCPAxyCand"), cpaXY);
    registry.fill(HIST("Data/hCPAxyCandVsPtCand"), cpaXY, ptCand);
    registry.fill(HIST("Data/hDecLengthCand"), decayLength);
    registry.fill(HIST("Data/hDecLengthCandVsPtCand"), decayLength, ptCand);
    registry.fill(HIST("Data/hDecLengthXYCand"), decayLengthXY);
    registry.fill(HIST("Data/hDecLengthXYCandVsPtCand"), decayLengthXY, ptCand);
    registry.fill(HIST("Data/hCtCand"), ctLc);
    registry.fill(HIST("Data/hCtCandVsPtCand"), ctLc, ptCand);
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

      if (std::abs(candidate.flagMcMatchRec()) == KDecayChannelLcToK0sP) {
        fillCandHistograms(candidate);

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
      }
    }
  }

  /// Helper function to fill MC generated histograms
  template <typename CandMcGen, typename Coll>
  void fillHistosMcGen(CandMcGen const& mcParticles, Coll const& recoCollisions)
  {
    for (const auto& particle : mcParticles) {
      // if (etaCandMax >= 0. && std::abs(particle.eta()) > etaCandMax) {
      //   continue;
      // }
      if (std::abs(particle.flagMcMatchGen()) == KDecayChannelLcToK0sP) {
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus);
        if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
          continue;
        }
        registry.fill(HIST("MC/Gen/hPtCandGen"), particle.pt());
        registry.fill(HIST("MC/Gen/hEtaCandGen"), particle.eta());
        registry.fill(HIST("MC/Gen/hPhiCandGen"), particle.phi());
        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/Gen/hPtCandGenPrompt"), particle.pt());
        } else if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("MC/Gen/hPtCandGenNonPrompt"), particle.pt());
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
                    soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& mcParticles,
                    aod::TracksWMc const&,
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
                       soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& mcParticles,
                       aod::TracksWMc const&,
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
                            soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& mcParticles,
                            aod::TracksWMc const&,
                            TracksWPid const&)
  {
    runAnalysisPerCollisionMc<false>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcStdWithFT0C, "Process MC with FT0C centrality", false);

  // MC Std with FT0M centrality
  void processMcStdWithFT0M(CollisionsMcWithFT0M const& collisions,
                            FilteredCandLcToPK0SWSelFlagAndMc const& candidates,
                            soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& mcParticles,
                            aod::TracksWMc const&,
                            TracksWPid const&)
  {
    runAnalysisPerCollisionMc<false>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcStdWithFT0M, "Process MC with FT0M centrality", false);

  // MC with ML + FT0C centrality
  void processMcWithMlWithFT0C(CollisionsMcWithFT0C const& collisions,
                               FilteredCandLcToPK0SWSelFlagAndMcAndMl const& candidates,
                               soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& mcParticles,
                               aod::TracksWMc const&,
                               TracksWPid const&)
  {
    runAnalysisPerCollisionMc<true>(collisions, candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLcToK0sP, processMcWithMlWithFT0C, "Process MC with ML and FT0C centrality", false);

  // MC with ML + FT0M centrality
  void processMcWithMlWithFT0M(CollisionsMcWithFT0M const& collisions,
                               FilteredCandLcToPK0SWSelFlagAndMcAndMl const& candidates,
                               soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& mcParticles,
                               aod::TracksWMc const&,
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
