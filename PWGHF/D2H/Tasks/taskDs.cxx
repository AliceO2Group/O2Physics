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

/// \file taskDs.cxx
/// \brief Ds± analysis task
/// \note Extended from taskD0 and taskDplus
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Universita and INFN Torino
/// \author Stefano Politanò <stefano.politano@cern.ch>, Politecnico & INFN Torino

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

/// Ds± analysis task
struct HfTaskDs {
  Configurable<int> decayChannel{"decayChannel", 1, "Switch between decay channels: 1 for Ds->PhiPi->KKpi, 2 for Ds->K0*K->KKPi"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Filter filterDsFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi))) != static_cast<uint8_t>(0);

  Partition<CandDsData> selectedDsToKKPiCand = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs;
  Partition<CandDsData> selectedDsToPiKKCand = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  Partition<CandDsMcReco> reconstructedCandSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && aod::hf_cand_3prong::flagMcDecayChanRec == decayChannel;
  Partition<CandDsMcReco> reconstructedCandBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi));

  HistogramRegistry registry{
    "registry",
    {{"hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hCPARecSig", "3-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBkg", "3-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "3-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaRecBkg", "3-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    AxisSpec ybins = {100, -5., 5, "#it{y}"};
    registry.add("hMass", "3-prong candidates;inv. mass (KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{400, 1.77, 2.17}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "3-prong candidates;proper lifetime (D_{s}^{#pm}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 100}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthXY", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecayLengthXY", "3-prong candidates;norm. decay length xy;entries", {HistType::kTH2F, {{80, 0., 80.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "3-prong candidates;cos. pointing angle;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxy", "3-prong candidates;cos. pointing angle xy;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterXY", "3-prong candidates;impact parameter xy (cm);entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMaxNormalisedDeltaIP", "3-prong candidates;norm. IP;entries", {HistType::kTH2F, {{200, -20., 20.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCos3PiK", "3-prong candidates;cos^{3} #theta'(K);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hAbsCos3PiK", "3-prong candidates;|cos^{3} #theta'(K)|;entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeltaMassPhi", "3-prong candidates;|M(KK) - M(#phi)| (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{100, 0., 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterProngSqSum", "3-prong candidates;squared sum of prong imp. par. (cm^{2});entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthError", "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthXYError", "3-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterError", "3-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "3-prong candidates;prong 0 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "3-prong candidates;prong 1 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2", "3-prong candidates;prong 2 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigPrompt", "3-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecBkg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenSig", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtVsYRecSigRecoPID", "3-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigPromptRecoPID", "3-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigNonPromptRecoPID", "3-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigRecoTopol", "3-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigPromptRecoTopol", "3-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigNonPromptRecoTopol", "3-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigRecoSkim", "3-prong candidates (RecoSkim - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigPromptRecoSkim", "3-prong candidates (RecoSkim - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigNonPromptRecoSkim", "3-prong candidates (RecoSkim - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGen", "MC particles (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
  }

  /// Fill histograms of quantities independent from the daugther-mass hypothesis
  /// \param candidate is candidate
  template <typename T1>
  void fillHisto(const T1& candidate)
  {
    auto pt = candidate.pt();
    registry.fill(HIST("hPt"), pt);
    registry.fill(HIST("hEta"), candidate.eta(), pt);
    registry.fill(HIST("hCt"), hfHelper.ctDs(candidate), pt);
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
    return;
  }

  /// Fill histograms of quantities for the KKPi daugther-mass hypothesis
  /// \param candidate is candidate
  template <typename T1>
  void fillHistoKKPi(const T1& candidate)
  {
    auto pt = candidate.pt();
    registry.fill(HIST("hMass"), hfHelper.invMassDsToKKPi(candidate), pt);
    registry.fill(HIST("hCos3PiK"), hfHelper.cos3PiKDsToKKPi(candidate), pt);
    registry.fill(HIST("hAbsCos3PiK"), std::abs(hfHelper.cos3PiKDsToKKPi(candidate)), pt);
    registry.fill(HIST("hDeltaMassPhi"), hfHelper.deltaMassPhiDsToKKPi(candidate), pt);
    return;
  }

  /// Fill histograms of quantities for the PiKK daugther-mass hypothesis
  /// \param candidate is candidate
  template <typename T1>
  void fillHistoPiKK(const T1& candidate)
  {
    auto pt = candidate.pt();
    registry.fill(HIST("hMass"), hfHelper.invMassDsToPiKK(candidate), pt);
    registry.fill(HIST("hCos3PiK"), hfHelper.cos3PiKDsToPiKK(candidate), pt);
    registry.fill(HIST("hAbsCos3PiK"), std::abs(hfHelper.cos3PiKDsToPiKK(candidate)), pt);
    registry.fill(HIST("hDeltaMassPhi"), hfHelper.deltaMassPhiDsToPiKK(candidate), pt);
    return;
  }

  /// Fill MC histograms at reconstruction level
  /// \param candidate is candidate
  /// \param flag is the selection flag, obtained from either isSelDsToKKPi() or isSelDsToPiKK()
  template <typename T1>
  void fillHistoMCRec(const T1& candidate, int flag)
  {
    auto pt = candidate.pt(); // rec. level pT
    auto y = hfHelper.yDs(candidate);

    registry.fill(HIST("hPtRecSig"), pt);
    registry.fill(HIST("hCPARecSig"), candidate.cpa());
    registry.fill(HIST("hEtaRecSig"), candidate.eta());
    registry.fill(HIST("hPtVsYRecSigRecoSkim"), pt, y);
    if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
      registry.fill(HIST("hPtVsYRecSigRecoTopol"), pt, y);
    }
    if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
      registry.fill(HIST("hPtVsYRecSigRecoPID"), pt, y);
    }

    // prompt
    if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtRecSigPrompt"), pt);
      registry.fill(HIST("hPtVsYRecSigPromptRecoSkim"), pt, y);
      if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
        registry.fill(HIST("hPtVsYRecSigPromptRecoTopol"), pt, y);
      }
      if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
        registry.fill(HIST("hPtVsYRecSigPromptRecoPID"), pt, y);
      }
    }

    // non-prompt
    if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtRecSigNonPrompt"), pt);
      registry.fill(HIST("hPtVsYRecSigNonPromptRecoSkim"), pt, y);
      if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
        registry.fill(HIST("hPtVsYRecSigNonPromptRecoTopol"), pt, y);
      }
      if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
        registry.fill(HIST("hPtVsYRecSigNonPromptRecoPID"), pt, y);
      }
    }

    return;
  }

  void process(CandDsData const& candidates)
  {
    for (const auto& candidate : selectedDsToKKPiCand) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillHisto(candidate);
      fillHistoKKPi(candidate);
    }

    for (const auto& candidate : selectedDsToPiKKCand) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillHisto(candidate);
      fillHistoPiKK(candidate);
    }
  }

  void processMc(CandDsMcReco const& candidates,
                 CandDsMcGen const& mcParticles,
                 aod::TracksWMc const&)
  {
    // MC rec.
    for (const auto& candidate : reconstructedCandSig) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
        continue;
      }

      auto prong0McPart = candidate.prong0_as<aod::TracksWMc>().mcParticle_as<CandDsMcGen>();
      auto indexMother = RecoDecay::getMother(mcParticles, prong0McPart, pdg::Code::kDS, true);
      auto particleMother = mcParticles.iteratorAt(indexMother);
      registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT

      // KKPi
      if (candidate.isCandidateSwapped() == 0) { // 0 corresponds to KKPi
        fillHistoMCRec(candidate, candidate.isSelDsToKKPi());
      }

      // PiKK
      if (candidate.isCandidateSwapped() == 1) { // 1 corresponds to PiKK
        fillHistoMCRec(candidate, candidate.isSelDsToPiKK());
      }
    }

    for (const auto& candidate : reconstructedCandBkg) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
        continue;
      }

      registry.fill(HIST("hPtRecBkg"), candidate.pt());
      registry.fill(HIST("hCPARecBkg"), candidate.cpa());
      registry.fill(HIST("hEtaRecBkg"), candidate.eta());
    }
    // TODO: add histograms for reflections

    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        if (particle.flagMcDecayChanGen() != decayChannel) {
          continue;
        }
        auto pt = particle.pt();
        auto y = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::analysis::pdg::MassDS);
        if (yCandGenMax >= 0. && std::abs(y) > yCandGenMax) {
          continue;
        }

        registry.fill(HIST("hPtGen"), pt);
        registry.fill(HIST("hPtVsYGen"), pt, y);
        registry.fill(HIST("hEtaGen"), particle.eta());

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hPtGenPrompt"), pt);
          registry.fill(HIST("hPtVsYGenPrompt"), pt, y);
        }

        if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("hPtGenNonPrompt"), pt);
          registry.fill(HIST("hPtVsYGenNonPrompt"), pt, y);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDs, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDs>(cfgc)};
}
