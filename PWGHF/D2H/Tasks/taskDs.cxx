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
/// \author Stefano Politanò <stefano.politano@cern.ch>, Politecnico & INFN Torino

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;

/// Ds± analysis task
struct HfTaskDs {
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits"};

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>> selectedDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>> recoFlagDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  HistogramRegistry registry{
    "registry",
    {{"hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}}},
     {"hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}}},
     {"hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}}},
     {"hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}}},
     {"hCPARecSig", "3-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{100, -1., 1.}}}},
     {"hCPARecBg", "3-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{100, -1., 1.}}}},
     {"hEtaRecSig", "3-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaRecBg", "3-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}}}};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    AxisSpec ybins = {100, -5., 5, "#it{y}"};
    registry.add("hMass", "3-prong candidates;inv. mass (K K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{350, 1.7, 2.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "3-prong candidates;proper lifetime (D_{s}^{#pm}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 100}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthXY", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecayLengthXY", "3-prong candidates;norm. decay length xy;entries", {HistType::kTH2F, {{80, 0., 80.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "3-prong candidates;cos. pointing angle;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxy", "3-prong candidates;cos. pointing angle xy;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterXY", "3-prong candidates;impact parameter xy (cm);entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMaxNormalisedDeltaIP", "3-prong candidates;norm. IP;entries", {HistType::kTH2F, {{200, -20., 20.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
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
    registry.add("hPtRecBg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
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

  void process(soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi> const& candidates)
  {
    for (auto& candidate : selectedDsCandidates) {
      // not possible in Filter since expressions do not support binary operators
      if (!(candidate.hfflag() & 1 << DecayType::DsToKKPi)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yDs(candidate)) > yCandMax) {
        continue;
      }
      registry.fill(HIST("hMass"), invMassDsToKKPi(candidate), candidate.pt());
      registry.fill(HIST("hPt"), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hCt"), ctDs(candidate), candidate.pt());
      registry.fill(HIST("hDecayLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecayLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hNormalisedDecayLengthXY"), candidate.decayLengthXYNormalised(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hCPAxy"), candidate.cpaXY(), candidate.pt());
      registry.fill(HIST("hImpactParameterXY"), candidate.impactParameterXY(), candidate.pt());
      registry.fill(HIST("hMaxNormalisedDeltaIP"), candidate.maxNormalisedDeltaIP(), candidate.pt());
      registry.fill(HIST("hImpactParameterProngSqSum"), candidate.impactParameterProngSqSum(), candidate.pt());
      registry.fill(HIST("hDecayLengthError"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecayLengthXYError"), candidate.errorDecayLengthXY(), candidate.pt());
      registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter2(), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hPtProng2"), candidate.ptProng2());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0Prong2"), candidate.impactParameter2(), candidate.pt());
    }
  }

  void processMc(soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particlesMC, aod::BigTracksMC const&)
  {
    // MC rec.
    for (auto& candidate : recoFlagDsCandidates) {
      // not possible in Filter since expressions do not support binary operators
      if (!(candidate.hfflag() & 1 << DecayType::DsToKKPi)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yDs(candidate)) > yCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMcMatchRec()) == 1 << DecayType::DsToKKPi) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.prong0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>(), pdg::Code::kDS, true);
        auto particleMother = particlesMC.iteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT
        auto ptRec = candidate.pt();
        auto yRec = yDs(candidate);
        auto DsToKKPi = candidate.isSelDsToKKPi();
        auto DsToPiKK = candidate.isSelDsToPiKK();
        registry.fill(HIST("hPtvsYRecSigRecoSkim"), ptRec, yRec);
        if (TESTBIT(DsToKKPi, aod::SelectionStep::RecoTopol) || TESTBIT(DsToPiKK, aod::SelectionStep::RecoTopol)) {
          registry.fill(HIST("hPtvsYRecSigRecoTopol"), ptRec, yRec);
        }
        if (TESTBIT(DsToKKPi, aod::SelectionStep::RecoPID) || TESTBIT(DsToPiKK, aod::SelectionStep::RecoPID)) {
          registry.fill(HIST("hPtvsYRecSigRecoPID"), ptRec, yRec);
        }
        registry.fill(HIST("hPtRecSig"), ptRec); // rec. level pT
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hPtvsYRecSigPromptRecoSkim"), ptRec, yRec);
          if (TESTBIT(DsToKKPi, aod::SelectionStep::RecoTopol) || TESTBIT(DsToPiKK, aod::SelectionStep::RecoTopol)) {
            registry.fill(HIST("hPtvsYRecSigPromptRecoTopol"), ptRec, yRec);
          }
          if (TESTBIT(DsToKKPi, aod::SelectionStep::RecoPID) || TESTBIT(DsToPiKK, aod::SelectionStep::RecoPID)) {
            registry.fill(HIST("hPtvsYRecSigPromptRecoPID"), ptRec, yRec);
          }
          registry.fill(HIST("hPtRecSigPrompt"), ptRec); // rec. level pT, prompt
        } else {                                         // FD
          registry.fill(HIST("hPtvsYRecSigNonPromptRecoSkim"), ptRec, yRec);
          if (TESTBIT(DsToKKPi, aod::SelectionStep::RecoTopol) || TESTBIT(DsToPiKK, aod::SelectionStep::RecoTopol)) {
            registry.fill(HIST("hPtvsYRecSigNonPromptRecoTopol"), ptRec, yRec);
          }
          if (TESTBIT(DsToKKPi, aod::SelectionStep::RecoPID) || TESTBIT(DsToPiKK, aod::SelectionStep::RecoPID)) {
            registry.fill(HIST("hPtvsYRecSigNonPromptRecoPID"), ptRec, yRec);
          }
          registry.fill(HIST("hPtRecSigNonPrompt"), ptRec); // rec. level pT, non-prompt
        }
        registry.fill(HIST("hCPARecSig"), candidate.cpa());
        registry.fill(HIST("hEtaRecSig"), candidate.eta());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa());
        registry.fill(HIST("hEtaRecBg"), candidate.eta());
      }
    }
    // MC gen.
    // Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::DsToKKPi) {
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        if (yCandMax >= 0. && std::abs(yGen) > yCandMax) {
          continue;
        }
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
    }
  }
  PROCESS_SWITCH(HfTaskDs, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDs>(cfgc)};
}
