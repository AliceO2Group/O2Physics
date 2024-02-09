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

/// Ds± analysis task
struct HfTaskDs {
  Configurable<int> decayChannelDs{"decayChannelDs", 1, "Switch between decay channels: 1 for Ds->PhiPi->KKpi, 2 for Ds->K0*K->KKPi"};
  Configurable<int> decayChannelDplusMc{"decayChannelDplusMc", 1, "Switch between decay channels for MC Dplus: 1 for Dplus->PhiPi->KKpi, 2 for Dplus->K0*K->KKPi"};
  Configurable<bool> FillDplusMc{"fillDplusMc", false, "Switch to fill Dplus MC information"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  int offsetDplusDecayChannel = aod::hf_cand_3prong::DecayChannelDToKKPi::DplusToPhiPi - aod::hf_cand_3prong::DecayChannelDToKKPi::DsToPhiPi; // Offset between Dplus and Ds to use the same decay channel. See aod::hf_cand_3prong::DecayChannelDToKKPi

  Filter filterDsFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi))) != static_cast<uint8_t>(0);

  Partition<CandDsData> selectedDsToKKPiCand = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs;
  Partition<CandDsData> selectedDsToPiKKCand = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  Partition<CandDsMcReco> reconstructedCandSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && (aod::hf_cand_3prong::flagMcDecayChanRec == decayChannelDs || (FillDplusMc && aod::hf_cand_3prong::flagMcDecayChanRec == (decayChannelDplusMc + offsetDplusDecayChannel)));
  Partition<CandDsMcReco> reconstructedCandBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi));

  HistogramRegistry registry{
    "registry",
    {{"hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hCPARecSigDs", "3-prong Ds candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecSigDplus", "3-prong Dplus candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBkg", "3-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSigDs", "3-prong Ds candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaRecSigDplus", "3-prong Dplus candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaRecBkg", "3-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaGenDs", "MC Ds particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaGenDplus", "MC Dplus particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    AxisSpec ybins = {100, -5., 5, "#it{y}"};
    registry.add("hMass", "3-prong candidates;inv. mass (KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{600, 1.67, 2.27}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
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
    registry.add("hMassKK", "3-prong candidates;M(KK) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{100, o2::constants::physics::MassPhi - 0.05, o2::constants::physics::MassPhi + 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterProngSqSum", "3-prong candidates;squared sum of prong imp. par. (cm^{2});entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthError", "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthXYError", "3-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterError", "3-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "3-prong candidates;prong 0 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "3-prong candidates;prong 1 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2", "3-prong candidates;prong 2 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigDs", "3-prong candidates (Ds matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigDsPrompt", "3-prong candidates (Ds matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigDsNonPrompt", "3-prong candidates (Ds matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigDplus", "3-prong candidates (Dplus matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigDplusPrompt", "3-prong candidates (Dplus matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigDplusNonPrompt", "3-prong candidates (Dplus matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecBkg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenDs", "MC Ds particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenDplus", "MC Dplus particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenSigDs", "MC Ds particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenSigDplus", "MC Dplus particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenDsPrompt", "MC Ds particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenDsNonPrompt", "MC Ds particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenDplusPrompt", "MC Dplus particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenDplusNonPrompt", "MC Dplus particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtVsYRecSigDsRecoPID", "3-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDsPromptRecoPID", "3-prong candidates (RecoPID - Ds matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDsNonPromptRecoPID", "3-prong candidates (RecoPID - Ds matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDsRecoTopol", "3-prong candidates (RecoTopol - Ds matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDsPromptRecoTopol", "3-prong candidates (RecoTopol - Ds matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDsNonPromptRecoTopol", "3-prong candidates (RecoTopol - Ds matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDsRecoSkim", "3-prong candidates (RecoSkim - Ds matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDsPromptRecoSkim", "3-prong candidates (RecoSkim - Ds matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDsNonPromptRecoSkim", "3-prong candidates (RecoSkim - Ds matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusRecoPID", "3-prong candidates (RecoPID - Dplus matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusPromptRecoPID", "3-prong candidates (RecoPID - Dplus matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusNonPromptRecoPID", "3-prong candidates (RecoPID - Dplus matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusRecoTopol", "3-prong candidates (RecoTopol - Dplus matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusPromptRecoTopol", "3-prong candidates (RecoTopol - Dplus matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusNonPromptRecoTopol", "3-prong candidates (RecoTopol - Dplus matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusRecoSkim", "3-prong candidates (RecoSkim - Dplus matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusPromptRecoSkim", "3-prong candidates (RecoSkim - Dplus matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYRecSigDplusNonPromptRecoSkim", "3-prong candidates (RecoSkim - Dplus matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGenDs", "MC Ds particles (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGenDsPrompt", "MC Ds particles (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGenDsNonPrompt", "MC Ds particles (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGenDplus", "MC Dplus particles (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGenDplusPrompt", "MC Dplus particles (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
    registry.add("hPtVsYGenDplusNonPrompt", "MC Dplus particles (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {ybins}}});
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
    registry.fill(HIST("hMassKK"), hfHelper.massKKPairDsToKKPi(candidate), pt);
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
    registry.fill(HIST("hMassKK"), hfHelper.massKKPairDsToPiKK(candidate), pt);
    return;
  }

  /// Fill MC histograms at reconstruction level
  /// \param candidate is candidate
  /// \param flag is the selection flag, obtained from either isSelDsToKKPi() or isSelDsToPiKK()
  template <typename T1>
  void fillHistoMCRec(const T1& candidate, int flag, bool isDplus)
  {
    auto pt = candidate.pt(); // rec. level pT
    // Dplus
    if (isDplus) {
      auto y = hfHelper.yDplus(candidate);

      registry.fill(HIST("hPtRecSigDplus"), pt);
      registry.fill(HIST("hCPARecSigDplus"), candidate.cpa());
      registry.fill(HIST("hEtaRecSigDplus"), candidate.eta());
      registry.fill(HIST("hPtVsYRecSigDplusRecoSkim"), pt, y);
      if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
        registry.fill(HIST("hPtVsYRecSigDplusRecoTopol"), pt, y);
      }
      if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
        registry.fill(HIST("hPtVsYRecSigDplusRecoPID"), pt, y);
      }

      // prompt
      if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hPtRecSigDplusPrompt"), pt);
        registry.fill(HIST("hPtVsYRecSigDplusPromptRecoSkim"), pt, y);
        if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
          registry.fill(HIST("hPtVsYRecSigDplusPromptRecoTopol"), pt, y);
        }
        if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
          registry.fill(HIST("hPtVsYRecSigDplusPromptRecoPID"), pt, y);
        }
      }

      // non-prompt
      if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
        registry.fill(HIST("hPtRecSigDplusNonPrompt"), pt);
        registry.fill(HIST("hPtVsYRecSigDplusNonPromptRecoSkim"), pt, y);
        if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
          registry.fill(HIST("hPtVsYRecSigDplusNonPromptRecoTopol"), pt, y);
        }
        if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
          registry.fill(HIST("hPtVsYRecSigDplusNonPromptRecoPID"), pt, y);
        }
      }

      return;
    }

    // Ds
    auto y = hfHelper.yDs(candidate);

    registry.fill(HIST("hPtRecSigDs"), pt);
    registry.fill(HIST("hCPARecSigDs"), candidate.cpa());
    registry.fill(HIST("hEtaRecSigDs"), candidate.eta());
    registry.fill(HIST("hPtVsYRecSigDsRecoSkim"), pt, y);
    if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
      registry.fill(HIST("hPtVsYRecSigDsRecoTopol"), pt, y);
    }
    if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
      registry.fill(HIST("hPtVsYRecSigDsRecoPID"), pt, y);
    }

    // prompt
    if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtRecSigDsPrompt"), pt);
      registry.fill(HIST("hPtVsYRecSigDsPromptRecoSkim"), pt, y);
      if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
        registry.fill(HIST("hPtVsYRecSigDsPromptRecoTopol"), pt, y);
      }
      if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
        registry.fill(HIST("hPtVsYRecSigDsPromptRecoPID"), pt, y);
      }
    }

    // non-prompt
    if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtRecSigDsNonPrompt"), pt);
      registry.fill(HIST("hPtVsYRecSigDsNonPromptRecoSkim"), pt, y);
      if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
        registry.fill(HIST("hPtVsYRecSigDsNonPromptRecoTopol"), pt, y);
      }
      if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
        registry.fill(HIST("hPtVsYRecSigDsNonPromptRecoPID"), pt, y);
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
      auto prong0McPart = candidate.prong0_as<aod::TracksWMc>().mcParticle_as<CandDsMcGen>();
      auto indexMother = RecoDecay::getMother(mcParticles, prong0McPart, o2::constants::physics::Pdg::kDS, true);
      if (indexMother != -1) {
        if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
          continue;
        }
        auto particleMother = mcParticles.iteratorAt(indexMother);
        registry.fill(HIST("hPtGenSigDs"), particleMother.pt()); // gen. level pT
        // KKPi
        if (candidate.isCandidateSwapped() == 0) { // 0 corresponds to KKPi
          fillHistoMCRec(candidate, candidate.isSelDsToKKPi(), false);
        }

        // PiKK
        if (candidate.isCandidateSwapped() == 1) { // 1 corresponds to PiKK
          fillHistoMCRec(candidate, candidate.isSelDsToPiKK(), false);
        }
      } else if (FillDplusMc) {
        if (yCandRecoMax >= 0. && std::abs(hfHelper.yDplus(candidate)) > yCandRecoMax) {
          continue;
        }
        indexMother = RecoDecay::getMother(mcParticles, prong0McPart, o2::constants::physics::Pdg::kDPlus, true);
        auto particleMother = mcParticles.iteratorAt(indexMother);
        registry.fill(HIST("hPtGenSigDplus"), particleMother.pt()); // gen. level pT
        // KKPi
        if (candidate.isCandidateSwapped() == 0) { // 0 corresponds to KKPi
          fillHistoMCRec(candidate, candidate.isSelDsToKKPi(), true);
        }

        // PiKK
        if (candidate.isCandidateSwapped() == 1) { // 1 corresponds to PiKK
          fillHistoMCRec(candidate, candidate.isSelDsToPiKK(), true);
        }
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
        if (particle.flagMcDecayChanGen() == decayChannelDs || (FillDplusMc && particle.flagMcDecayChanGen() == (decayChannelDplusMc + offsetDplusDecayChannel))) {
          auto pt = particle.pt();
          auto y = 0;
          if (particle.flagMcDecayChanGen() == decayChannelDs) {
            y = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassDS);
          } else if (FillDplusMc) {
            y = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassDPlus);
          }
          if (yCandGenMax >= 0. && std::abs(y) > yCandGenMax) {
            continue;
          }

          if (particle.flagMcDecayChanGen() == decayChannelDs) {
            registry.fill(HIST("hPtGenDs"), pt);
            registry.fill(HIST("hPtVsYGenDs"), pt, y);
            registry.fill(HIST("hEtaGenDs"), particle.eta());
            if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hPtGenDsPrompt"), pt);
              registry.fill(HIST("hPtVsYGenDsPrompt"), pt, y);
            }
            if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hPtGenDsNonPrompt"), pt);
              registry.fill(HIST("hPtVsYGenDsNonPrompt"), pt, y);
            }
          } else if (FillDplusMc) {
            registry.fill(HIST("hPtGenDplus"), pt);
            registry.fill(HIST("hPtVsYGenDplus"), pt, y);
            registry.fill(HIST("hEtaGenDplus"), particle.eta());
            if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hPtGenDplusPrompt"), pt);
              registry.fill(HIST("hPtVsYGenDplusPrompt"), pt, y);
            }
            if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hPtGenDplusNonPrompt"), pt);
              registry.fill(HIST("hPtVsYGenDplusNonPrompt"), pt, y);
            }
          }
        } else {
          continue;
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
