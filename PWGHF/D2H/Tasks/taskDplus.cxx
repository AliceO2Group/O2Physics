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

/// \file taskDplus.cxx
/// \brief D± analysis task
/// \note Extended from taskD0
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Luca Aglietta <luca.aglietta@cern.ch>, University and INFN Torino

#include "CommonConstants/PhysicsConstants.h"
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

/// D± analysis task
struct HfTaskDplus {
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for DPlus"}; // 7 corresponds to topo+PID cuts
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  ConfigurableAxis axisMlScore0{"axisMlScore0", {100, 0., 1.}, "axis for ML output score 0"};
  ConfigurableAxis axisMlScore1{"axisMlScore1", {100, 0., 1.}, "axis for ML output score 1"};
  ConfigurableAxis axisMlScore2{"axisMlScore2", {100, 0., 1.}, "axis for ML output score 2"};

  HfHelper hfHelper;

  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandDplusDataWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandDplusMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using CandDplusMcRecoWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec, aod::HfMlDplusToPiKPi>>;
  using McParticles = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Filter filterDplusFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi))) != static_cast<uint8_t>(0);

  // data
  Partition<CandDplusData> selectedDPlusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<CandDplusDataWithMl> selectedDPlusCandidatesWithMl = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  // Matched MC
  Partition<CandDplusMcReco> recoDPlusCandidates = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<CandDplusMcRecoWithMl> recoDPlusCandidatesWithMl = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  // MC Bkg
  Partition<CandDplusMcReco> recoBkgCandidates = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<CandDplusMcRecoWithMl> recoBkgCandidatesWithMl = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  // Generated particles

  HistogramRegistry registry{
    "registry",
    {{"hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hCPARecSig", "3-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBg", "3-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "3-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaRecBg", "3-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}}}};

  void init(InitContext&)
  {
    std::array<bool, 4> doprocess{doprocessData, doprocessDataWithMl, doprocessMc, doprocessMcWithMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function should be enabled! Please check your configuration!");
    }
    auto vbins = (std::vector<double>)binsPt;
    AxisSpec ptbins = {vbins, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massbins = {600, 1.67, 2.27, "inv. mass (K#pi#pi) (GeV/#it{c}^{2})"};
    registry.add("hMass", "3-prong candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{350, 1.7, 2.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "3-prong candidates;proper lifetime (D^{#pm}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthXY", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecayLengthXY", "3-prong candidates;norm. decay length xy;entries", {HistType::kTH2F, {{80, 0., 80.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "3-prong candidates;cos. pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxy", "3-prong candidates;cos. pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterXY", "3-prong candidates;impact parameter xy (cm);entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMaxNormalisedDeltaIP", "3-prong candidates;norm. IP;entries", {HistType::kTH2F, {{200, -20., 20.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterProngSqSum", "3-prong candidates;squared sum of prong imp. par. (cm^{2});entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthError", "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthXYError", "3-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterError", "3-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigPrompt", "3-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecBg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenSig", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtVsYRecSig_RecoPID", "3-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigPromptRecoPID", "3-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigNonPromptRecoPID", "3-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigRecoTopol", "3-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigPromptRecoTopol", "3-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigNonPromptRecoTopol", "3-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSig_RecoSkim", "3-prong candidates (RecoSkim - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigPrompt_RecoSkim", "3-prong candidates (RecoSkim - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigNonPrompt_RecoSkim", "3-prong candidates (RecoSkim - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYGen", "MC particles (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    if (doprocessDataWithMl) {
      registry.add("hSparseMass", "THn for Dplus", HistType::kTHnSparseF, {massbins, ptbins, axisMlScore0, axisMlScore1, axisMlScore2});
    }
    if (doprocessMcWithMl) {
      registry.add("hSparseMassPrompt", "THn for Dplus Prompt", HistType::kTHnSparseF, {massbins, ptbins, axisMlScore0, axisMlScore1, axisMlScore2});
      registry.add("hSparseMassFD", "THn for Dplus FD", HistType::kTHnSparseF, {massbins, ptbins, axisMlScore0, axisMlScore1, axisMlScore2});
      registry.add("hSparseMassBkg", "THn for Dplus Bkg", HistType::kTHnSparseF, {massbins, ptbins, axisMlScore0, axisMlScore1, axisMlScore2});
    }
  }

  // Fill histograms of quantities for the reconstructed Dplus candidates
  /// \param candidate is candidate
  template <typename T1>
  void fillHisto(const T1& candidate)
  {
    float pt = candidate.pt();
    registry.fill(HIST("hMass"), hfHelper.invMassDplusToPiKPi(candidate), pt);
    registry.fill(HIST("hPt"), pt);
    registry.fill(HIST("hEta"), candidate.eta(), pt);
    registry.fill(HIST("hCt"), hfHelper.ctDplus(candidate), pt);
    registry.fill(HIST("hDecayLength"), candidate.decayLength(), pt);
    registry.fill(HIST("hDecayLengthXY"), candidate.decayLengthXY(), pt);
    registry.fill(HIST("hNormalisedDecayLengthXY"), candidate.decayLengthXYNormalised(), pt);
    registry.fill(HIST("hCPA"), candidate.cpa(), pt);
    registry.fill(HIST("hCPAxy"), candidate.cpaXY(), pt);
    registry.fill(HIST("hImpactParameterXY"), candidate.impactParameterXY(), pt);
    registry.fill(HIST("hMaxNormalisedDeltaIP"), candidate.maxNormalisedDeltaIP(), pt);
    registry.fill(HIST("hImpactParameterProngSqSum"), candidate.impactParameterProngSqSum(), pt);
    registry.fill(HIST("hDecayLengthError"), candidate.errorDecayLength(), pt);
    registry.fill(HIST("hDecayLengthXYError"), candidate.errorDecayLengthXY(), pt);
    registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter0(), pt);
    registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter1(), pt);
    registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter2(), pt);
    registry.fill(HIST("hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1"), candidate.ptProng1());
    registry.fill(HIST("hPtProng2"), candidate.ptProng2());
    registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), pt);
    registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), pt);
    registry.fill(HIST("hd0Prong2"), candidate.impactParameter2(), pt);
  }

  // Fill THnSparses for the ML analysis
  /// \param candidate is a particle candidate
  template <bool isMc, bool isMatched, typename T1>
  void fillSparseML(const T1& candidate)
  {
    std::vector<float> outputMl = {-999., -999., -999.};
    for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
      outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
    }
    if constexpr (isMc) {
      if constexpr (isMatched) {
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hSparseMassPrompt"), hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
        } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("hSparseMassFD"), hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
        }
      } else {
        registry.fill(HIST("hSparseMassBkg"), hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
      }
    } else {
      registry.fill(HIST("hSparseMass"), hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
    }
  }

  // Fill histograms of quantities for the reconstructed Dplus candidates with MC matching
  /// \param candidate is candidate
  template <bool isMatched, typename T1>
  void fillHistoMCRec(const T1& candidate)
  {
    if constexpr (isMatched) {
      auto ptRec = candidate.pt();
      auto yRec = hfHelper.yDplus(candidate);
      registry.fill(HIST("hPtVsYRecSig_RecoSkim"), ptRec, yRec);
      if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoTopol)) {
        registry.fill(HIST("hPtVsYRecSigRecoTopol"), ptRec, yRec);
      }
      if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoPID)) {
        registry.fill(HIST("hPtVsYRecSig_RecoPID"), ptRec, yRec);
      }
      if (candidate.isSelDplusToPiKPi() >= selectionFlagDplus) {
        registry.fill(HIST("hPtRecSig"), ptRec); // rec. level pT
      }
      if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hPtVsYRecSigPrompt_RecoSkim"), ptRec, yRec);
        if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoTopol)) {
          registry.fill(HIST("hPtVsYRecSigPromptRecoTopol"), ptRec, yRec);
        }
        if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoPID)) {
          registry.fill(HIST("hPtVsYRecSigPromptRecoPID"), ptRec, yRec);
        }
        if (candidate.isSelDplusToPiKPi() >= selectionFlagDplus) {
          registry.fill(HIST("hPtRecSigPrompt"), ptRec); // rec. level pT, prompt
        }
      } else {
        registry.fill(HIST("hPtVsYRecSigNonPrompt_RecoSkim"), ptRec, yRec);
        if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoTopol)) {
          registry.fill(HIST("hPtVsYRecSigNonPromptRecoTopol"), ptRec, yRec);
        }
        if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoPID)) {
          registry.fill(HIST("hPtVsYRecSigNonPromptRecoPID"), ptRec, yRec);
        }
        if (candidate.isSelDplusToPiKPi() >= selectionFlagDplus) {
          registry.fill(HIST("hPtRecSigNonPrompt"), ptRec); // rec. level pT, non-prompt
        }
      }
      registry.fill(HIST("hCPARecSig"), candidate.cpa());
      registry.fill(HIST("hEtaRecSig"), candidate.eta());
    } else {
      registry.fill(HIST("hPtRecBg"), candidate.pt());
      registry.fill(HIST("hCPARecBg"), candidate.cpa());
      registry.fill(HIST("hEtaRecBg"), candidate.eta());
    }
  }

  // Fill histograms of quantities for generated Dplus particles
  /// \param particle is a particle with MC information
  template <typename T2>
  void fillHistoMCGen(const T2& particle)
  {
    auto ptGen = particle.pt();
    auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus);
    registry.fill(HIST("hPtGen"), ptGen);
    registry.fill(HIST("hPtVsYGen"), ptGen, yGen);
    if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtGenPrompt"), ptGen);
      registry.fill(HIST("hPtVsYGenPrompt"), ptGen, yGen);
    } else {
      registry.fill(HIST("hPtGenNonPrompt"), ptGen);
      registry.fill(HIST("hPtVsYGenNonPrompt"), ptGen, yGen);
    }
    registry.fill(HIST("hEtaGen"), particle.eta());
  }

  // Run analysis for the reconstructed Dplus candidates from data
  /// \param candidates are reconstructed candidates
  template <bool fillMl, typename T1>
  void runDataAnalysis(const T1& /*candidates*/)
  {
    if constexpr (!fillMl) {
      for (const auto& candidate : selectedDPlusCandidates) {
        if ((yCandRecoMax >= 0. && std::abs(hfHelper.yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHisto(candidate);
      }
    } else {
      for (const auto& candidate : selectedDPlusCandidatesWithMl) {
        if ((yCandRecoMax >= 0. && std::abs(hfHelper.yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHisto(candidate);
        fillSparseML<false, false>(candidate);
      }
    }
  }
  // Run analysis for the reconstructed Dplus candidates with MC matching
  /// \param candidates are reconstructed candidates
  /// \param mcParticles are particles with MC information
  template <bool fillMl, typename T1, typename T2>
  void runMCAnalysis(const T1& /*recoCandidates*/, const T2& mcParticles)
  {
    // MC rec. w/o Ml
    if constexpr (!fillMl) {
      for (const auto& candidate : recoDPlusCandidates) {
        if ((yCandRecoMax >= 0. && std::abs(hfHelper.yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHisto(candidate);
        fillHistoMCRec<true>(candidate);
      }
      // Bkg
      for (const auto& candidate : recoBkgCandidates) {
        if ((yCandRecoMax >= 0. && std::abs(hfHelper.yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHistoMCRec<false>(candidate);
      }
    } else {
      for (const auto& candidate : recoDPlusCandidatesWithMl) {
        if ((yCandRecoMax >= 0. && std::abs(hfHelper.yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHisto(candidate);
        fillHistoMCRec<true>(candidate);
        fillSparseML<true, true>(candidate);
      }
      // Bkg
      for (const auto& candidate : recoBkgCandidatesWithMl) {
        if ((yCandRecoMax >= 0. && std::abs(hfHelper.yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHistoMCRec<false>(candidate);
        fillSparseML<true, false>(candidate);
      }
    }
    // MC gen.
    for (const auto& particle : mcParticles) {
      auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus);
      if ((yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) || (std::abs(particle.flagMcMatchGen()) != 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi)) {
        continue;
      }
      fillHistoMCGen(particle);
    }
  }

  // process functions
  void processData(CandDplusData const& candidates)
  {
    runDataAnalysis<false>(candidates);
  }
  PROCESS_SWITCH(HfTaskDplus, processData, "Process data w/o ML", true);

  void processDataWithMl(CandDplusDataWithMl const& candidates)
  {
    runDataAnalysis<true>(candidates);
  }
  PROCESS_SWITCH(HfTaskDplus, processDataWithMl, "Process data with ML", false);

  void processMc(CandDplusMcReco const& candidates,
                 McParticles const& mcParticles)
  {
    runMCAnalysis<false>(candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDplus, processMc, "Process MC w/o ML", false);

  void processMcWithMl(CandDplusMcRecoWithMl const& candidates,
                       McParticles const& mcParticles)
  {
    runMCAnalysis<true>(candidates, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDplus, processMcWithMl, "Process MC with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDplus>(cfgc)};
}
