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

/// \file taskDstarToD0Pi.cxx
/// \file D* analysis task

/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Dstar analysis task
struct HfTaskDstarToD0Pi {
  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 & Pi"};
  Configurable<double> yCandDstarRecoMax{"yCandDstarRecoMax", 0.8, "max. candidate Dstar rapidity"};
  Configurable<double> yCandDstarGenMax{"yCandDstarGenMax", 0.5, "max. rapidity of Generator level Particle"};
  Configurable<bool> selectionFlagHfD0ToPiK{"selectionFlagHfD0ToPiK", true, "Selection Flag for HF flagged candidates"};
  Configurable<std::vector<double>> ptBins{"ptBins", std::vector<double>{hf_cuts_dstar_to_d0_pi::vecBinsPt}, "pT bin limits for Dstar"};

  using CandDstarWSelFlag = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi>;
  /// @brief specially for MC data
  // full reconstructed Dstar candidate
  using CandDstarWSelFlagMcRec = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfCandDstarMcRec>;
  using CandDstarMcGen = soa::Join<aod::McParticles, aod::HfCandDstarMcGen>;

  Partition<CandDstarWSelFlag> rowsSelectedCandDstar = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  Partition<CandDstarWSelFlagMcRec> rowsSelectedCandDstarMcRec = aod::hf_sel_candidate_dstar::isRecoD0Flag == selectionFlagHfD0ToPiK;

  ConfigurableAxis binningImpactParam{"binningImpactParam", {1000, 0.1, -0.1}, " Bins of Impact Parameter"};
  ConfigurableAxis binningDecayLength{"binningDecayLength", {1000, 0.0, 0.7}, "Bins of Decay Length"};
  ConfigurableAxis binningNormDecayLength{"binningNormDecayLength", {1000, 0.0, 40.0}, "Bins of Normalised Decay Length"};

  HistogramRegistry registry{
    "registry",
    {{"QA/hPtDstar", "Dstar Candidates; Dstar candidate #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtD0", "D0 Candiades; D0 Candidate #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtProng0D0", "Prong0 of D0 Candidates; Prong0 #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtProng1D0", "Prong1 of D0 Candidates; Prong1 #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtProng0D0Bar", "Prong0 of D0Bar Candidates; Prong0 #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtProng1D0Bar", "Prong1 of D0Bar Candidates; Prong1 #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtSoftPi", "Soft Pi ; Soft Pi #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, 0., 1.}}}}}};

  void init(InitContext&)
  {
    auto vecPtBins = (std::vector<double>)ptBins;

    AxisSpec axisImpactParam = {binningImpactParam, "impact parameter (cm)"};
    AxisSpec axisDecayLength = {binningDecayLength, " decay length (cm)"};
    AxisSpec axisNormDecayLength = {binningNormDecayLength, "normalised decay length (cm)"};

    registry.add("Yield/hDeltaInvMassDstar2D", "#Delta #it{M}_{inv} D* Candidate; inv. mass ((#pi #pi k) - (#pi k)) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0.13, 0.16}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}}, true);
    registry.add("Yield/hDeltaInvMassDstar1D", "#Delta #it{M}_{inv} D* Candidate; inv. mass ((#pi #pi k) - (#pi k)) (GeV/#it{c}^{2}); entries", {HistType::kTH1F, {{100, 0.13, 0.16}}}, true);
    registry.add("Yield/hInvMassDstar", "#Delta #it{M}_{inv} D* Candidate; inv. mass (#pi #pi k) (GeV/#it{c}^{2}); entries", {HistType::kTH1F, {{500, 0., 5.0}}}, true);
    registry.add("Yield/hInvMassD0", "#it{M}_{inv}D^{0} candidate;#it{M}_{inv} D^{0} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{500, 0., 5.0}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}}, true);
    // only QA
    registry.add("QA/hEtaDstar", "D* Candidate; D* Candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hCtD0", "D0 Candidate; D0 Candidate's proper life time #it{c}#tau (cm) ; entries ", {HistType::kTH2F, {{1000, -0.1, 14.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthD0", "D0 Candidate; D0 Candidate's decay length (cm); entries", {HistType::kTH2F, {axisDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthXYD0", "D0 Candidate; D0 Candidate's decay length xy (cm);  entries", {HistType::kTH2F, {axisDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthNormalisedD0", "D0 Candidates;Do Candidate's normalised decay length (cm); entries", {HistType::kTH2F, {axisNormDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthXYNormalisedD0", "D0 candidate; D0 Candidate's normalised decay length xy (cm); entries", {HistType::kTH2F, {axisNormDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hCPAD0", "D0 Candidates; D0 Candidate's cosine pointing angle; entries", {HistType::kTH2F, {{110, -1., 1.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hCPAxyD0", "D0 candidates; D0 Candidate's cosine of pointing angle xy; entries", {HistType::kTH2F, {{110, -1., 1.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hImpactParameterXYD0", "D0 Candidates; D0 Candidate's reconstructed impact parameter xy (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDeltaIPMaxNormalisedD0", "D0 Candidate; D0 Candidate's Norm. Delta IP; entries", {HistType::kTH2F, {{1000, -20., 20.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hSqSumProngsImpactParameterD0", "D0 Candidates; Sqr Sum of Impact params of D0 Prongs; entries ", {HistType::kTH2F, {{1000, 0., 0.25}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthErrorD0", "D0 Candidates; D0 Candidate's Decay Length Error (cm); entries", {HistType::kTH2F, {axisDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthXYErrorD0", "D0 Candidates; D0 Candidate's Decay Length Error XY (cm); entries", {HistType::kTH2F, {axisDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hImpactParameterError", "D0 Prongs; Impact param error of different D0 Prongs (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hd0Prong0", "Prong0; DCAxy of Prong0 (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hd0Prong1", "Prong1; DCAxy of Prong1 (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hd0ProngSoftPi", "ProngSoftPi; DCAxy of Prong Soft Pi (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // MC Matching at Reconstruction Level Successful
    registry.add("QA/hCPASkimD0RecSig", "MC Matched Skimed D* Candidates at Reconstruction Level; cosine of pointing angle", {HistType::kTH1F, {{110, -1.1, 1.1}}});
    registry.add("QA/hEtaSkimD0RecSig", "MC Matched Skimed D* Candidates at Reconstruction Level; #it{#eta} of D0 Prong", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("QA/hEtaSkimDstarRecSig", "MC Matched Skimed D* Candidates at Reconstruction Level; #it{#eta} of D* Candidate", {HistType::kTH1F, {{100, -2., 2.}}});
    // pt vs y
    registry.add("QA/hPtSkimDstarGenSig", "MC Matched Skimed D* Reconstructed Candidates at Generator Level; #it{p}_{T} of D* at Generator Level (GeV/#it{c})", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hPtVsYSkimDstarRecSig", "MC Matched Skimed D* Candidates at Reconstruction Level; #it{p}_{T} of D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoTopolDstarRecSig", "MC Matched RecoTopol D* Candidates at Reconstruction Level; #it{p}_{T} of D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoPidDstarRecSig", "MC Matched RecoPid D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtFullRecoDstarRecSig", "MC Matched FullReco D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c})", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // Only Prompt RecSig
    registry.add("QA/hPtVsYSkimPromptDstarRecSig", "MC Matched Skimed Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}; #it{y})", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoTopolPromptDstarRecSig", "MC Matched RecoTopol Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoPidPromptDstarRecSig", "MC Matched RecoPid Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtFullRecoPromptDstarRecSig", "MC Matched FullReco Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c})", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // Only Non-Prompt RecSig
    registry.add("QA/hPtVsYSkimNonPromptDstarRecSig", "MC Matched Skimed Non-Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoTopolNonPromptDstarRecSig", "MC Matched RecoTopol Non-Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoPidNonPromptDstarRecSig", "MC Matched RecoPid Non-Prompt D*  Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYFullRecoNonPromptDstarRecSig", "MC Matched FullReco Non-Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c})", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // MC Matching UnSuccessful
    registry.add("QA/hCPASkimD0RecBg", "MC UnMatched Skimmed D0 Candidates at Reconstruction Level; cosine of pointing angle", {HistType::kTH1F, {{110, -1., 1.}}});
    registry.add("QA/hEtaSkimD0RecBg", "MC UnMatched Skimmed D0 Candidates at Reconstruction Level; #it{#eta} of D0 Prong", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("QA/hEtaSkimDstarRecBg", "MC UnMatched Skimmed D* Candidates at Reconstruction Level; #it{#eta} of D* Candidate", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("QA/hPtSkimD0RecBg", "MC UnMatched Skimmed D0 Candidates at Reconstruction Level;  #it{p}_{T} of  D0", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hPtSkimDstarRecBg", "MC UnMatched Skimmed D* Candidates at Reconstruction Level; #it{p}_{T} of  D*", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // MC Matching at Generator level Successful
    registry.add("QA/hEtaDstarGen", "MC Matched D* Candidates at Generator Level; #it{#eta}", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("QA/hPtDstarGen", "MC Matched D* Candidates at Generator Level; #it{p}_{T} of  D*", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hPtVsYDstarGen", "MC Matched D* Candidates at Generator Level; #it{p}_{T} of  D*; #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    // Prompt Gen
    registry.add("QA/hPtPromptDstarGen", "MC Matched Prompt D* Candidates at Generator Level; #it{p}_{T} of  D*", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hPtVsYPromptDstarGen", "MC Matched Prompt D* Candidates at Generator Level; #it{p}_{T} of  D*; #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    // Non Prmpt Gen
    registry.add("QA/hPtNonPromptDstarGen", "MC Matched Non-Prompt D* Candidates at Generator Level; #it{p}_{T} of  D*", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hPtVsYNonPromptDstarGen", "MC Matched Non-Prompt D* Candidates at Generator Level; #it{p}_{T} of  D*; #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
  }

  void process(CandDstarWSelFlag const&)
  {
    for (const auto& candDstar : rowsSelectedCandDstar) {
      auto yDstar = candDstar.y(constants::physics::MassDStar);
      if (yCandDstarRecoMax >= 0. && std::abs(yDstar) > yCandDstarRecoMax) {
        continue;
      }
      registry.fill(HIST("QA/hPtDstar"), candDstar.pt());
      registry.fill(HIST("QA/hPtD0"), candDstar.ptD0());
      registry.fill(HIST("QA/hPtSoftPi"), candDstar.ptSoftPi());
      registry.fill(HIST("QA/hEtaDstar"), candDstar.eta(), candDstar.pt());
      registry.fill(HIST("QA/hCtD0"), candDstar.ctD0(), candDstar.pt());
      registry.fill(HIST("QA/hDecayLengthD0"), candDstar.decayLengthD0(), candDstar.pt());
      registry.fill(HIST("QA/hDecayLengthXYD0"), candDstar.decayLengthXYD0(), candDstar.pt());
      registry.fill(HIST("QA/hDecayLengthNormalisedD0"), candDstar.decayLengthNormalisedD0(), candDstar.pt());
      registry.fill(HIST("QA/hDecayLengthXYNormalisedD0"), candDstar.decayLengthXYNormalisedD0(), candDstar.pt());
      registry.fill(HIST("QA/hCPAD0"), candDstar.cpaD0(), candDstar.pt());
      registry.fill(HIST("QA/hCPAxyD0"), candDstar.cpaXYD0(), candDstar.pt());
      registry.fill(HIST("QA/hImpactParameterXYD0"), candDstar.impactParameterXYD0(), candDstar.pt());
      registry.fill(HIST("QA/hDeltaIPMaxNormalisedD0"), candDstar.deltaIPNormalisedMaxD0(), candDstar.pt());
      registry.fill(HIST("QA/hSqSumProngsImpactParameterD0"), candDstar.impactParameterProngSqSumD0(), candDstar.pt());
      registry.fill(HIST("QA/hDecayLengthErrorD0"), candDstar.errorDecayLengthD0(), candDstar.pt());
      registry.fill(HIST("QA/hDecayLengthXYErrorD0"), candDstar.errorDecayLengthXYD0(), candDstar.pt());
      registry.fill(HIST("QA/hImpactParameterError"), candDstar.errorImpactParameter0(), candDstar.pt());
      registry.fill(HIST("QA/hImpactParameterError"), candDstar.errorImpactParameter1(), candDstar.pt());
      registry.fill(HIST("QA/hImpactParameterError"), candDstar.errorImpParamSoftPi(), candDstar.pt());
      registry.fill(HIST("QA/hd0Prong0"), candDstar.impactParameter0(), candDstar.pt());
      registry.fill(HIST("QA/hd0Prong1"), candDstar.impactParameter1(), candDstar.pt());
      registry.fill(HIST("QA/hd0ProngSoftPi"), candDstar.impParamSoftPi(), candDstar.pt());

      auto invDstar = candDstar.invMassDstar();
      auto invAntiDstar = candDstar.invMassAntiDstar();
      auto invD0 = candDstar.invMassD0();
      auto invD0Bar = candDstar.invMassD0Bar();

      auto signDstar = candDstar.signSoftPi();
      if (signDstar > 0) {
        registry.fill(HIST("Yield/hDeltaInvMassDstar2D"), (invDstar - invD0), candDstar.pt());
        registry.fill(HIST("Yield/hInvMassD0"), invD0, candDstar.ptD0());
        registry.fill(HIST("Yield/hDeltaInvMassDstar1D"), (invDstar - invD0));
        registry.fill(HIST("Yield/hInvMassDstar"), invDstar);
        // filling pt of two pronges of D0
        registry.fill(HIST("QA/hPtProng0D0"), candDstar.ptProng0());
        registry.fill(HIST("QA/hPtProng1D0"), candDstar.ptProng1());
      } else if (signDstar < 0) {
        registry.fill(HIST("Yield/hDeltaInvMassDstar2D"), (invAntiDstar - invD0Bar), candDstar.pt());
        registry.fill(HIST("Yield/hInvMassD0"), invD0Bar, candDstar.ptD0());
        registry.fill(HIST("Yield/hDeltaInvMassDstar1D"), (invAntiDstar - invD0Bar));
        registry.fill(HIST("Yield/hInvMassDstar"), invAntiDstar);
        // filling pt of two pronges of D0Bar
        registry.fill(HIST("QA/hPtProng0D0Bar"), candDstar.ptProng0());
        registry.fill(HIST("QA/hPtProng1D0Bar"), candDstar.ptProng1());
      }
    }
  }

  void processMC(CandDstarWSelFlagMcRec const&,
                 CandDstarMcGen const& rowsMcPartilces,
                 aod::TracksWMc const&)
  {

    int8_t signDstar = 0;
    // MC at Reconstruction level
    for (const auto& candDstarMcRec : rowsSelectedCandDstarMcRec) {
      auto ptDstarRecSig = candDstarMcRec.pt();
      auto yDstarRecSig = candDstarMcRec.y(constants::physics::MassDStar);
      if (yCandDstarRecoMax >= 0. && std::abs(yDstarRecSig) > yCandDstarRecoMax) {
        continue;
      }
      if (TESTBIT(std::abs(candDstarMcRec.flagMcMatchRec()), aod::hf_cand_dstar::DecayType::DstarToD0Pi)) { // if MC matching is successful at Reconstruction Level
        // get MC Mother particle
        auto indexMother = RecoDecay::getMother(rowsMcPartilces, candDstarMcRec.prong0_as<aod::TracksWMc>().mcParticle_as<CandDstarMcGen>(), o2::constants::physics::Pdg::kDStar, true, &signDstar, 2);
        auto particleMother = rowsMcPartilces.rawIteratorAt(indexMother); // What is difference between rawIterator() or iteratorAt() methods?
        registry.fill(HIST("QA/hPtDstarGenSig"), particleMother.pt());    // generator level pt

        registry.fill(HIST("QA/hPtVsYSkimDstarRecSig"), ptDstarRecSig, yDstarRecSig); // Skimed at level of trackIndexSkimCreator
        if (candDstarMcRec.isRecoTopol()) {                                           // if Topological selection are passed
          registry.fill(HIST("QA/hPtVsYRecoTopolDstarRecSig"), ptDstarRecSig, yDstarRecSig);
        }
        if (candDstarMcRec.isRecoPid()) { // if PID selection is passed
          registry.fill(HIST("QA/hPtVsYRecoPidDstarRecSig"), ptDstarRecSig, yDstarRecSig);
        }
        if (candDstarMcRec.isRecoCand()) { // if all selection passed
          registry.fill(HIST("QA/hPtFullRecoDstarRecSig"), ptDstarRecSig);
        }
        registry.fill(HIST("QA/hCPASkimD0RecSig"), candDstarMcRec.cpaD0());
        registry.fill(HIST("QA/hEtaSkimD0RecSig"), candDstarMcRec.etaD0());
        registry.fill(HIST("QA/hEtaSkimDstarRecSig"), candDstarMcRec.eta());

        // only prompt signal at reconstruction level
        if (candDstarMcRec.originMcRec() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("QA/hPtVsYSkimPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig); // Skimed at level of trackIndexSkimCreator
          if (candDstarMcRec.isRecoTopol()) {                                                 // if Topological selection are passed
            registry.fill(HIST("QA/hPtVsYRecoTopolPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          }
          if (candDstarMcRec.isRecoPid()) { // if PID selection is passed
            registry.fill(HIST("QA/hPtVsYRecoPidPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          }
          if (candDstarMcRec.isRecoCand()) { // if all selection passed
            registry.fill(HIST("QA/hPtFullRecoPromptDstarRecSig"), ptDstarRecSig);
          }

        } else if (candDstarMcRec.originMcRec() == RecoDecay::OriginType::NonPrompt) { // only non-prompt signal at reconstruction level
          registry.fill(HIST("QA/hPtVsYSkimNonPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          if (candDstarMcRec.isRecoTopol()) { // if Topological selection are passed
            registry.fill(HIST("QA/hPtVsYRecoTopolNonPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          }
          if (candDstarMcRec.isRecoPid()) {
            registry.fill(HIST("QA/PtVsYRecoPidNonPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          }
          if (candDstarMcRec.isRecoCand()) {
            registry.fill(HIST("QA/PtFullRecoPromptDstarRecSig"), ptDstarRecSig);
          }
        }
      } else { // MC Unmatched (Baground at Reconstruction Level)
        registry.fill(HIST("QA/hCPASkimD0RecBg"), candDstarMcRec.cpaD0());
        registry.fill(HIST("QA/hEtaSkimD0RecBg"), candDstarMcRec.etaD0());
        registry.fill(HIST("QA/hEtaSkimDstarRecBg"), candDstarMcRec.eta());
        registry.fill(HIST("QA/hPtSkimD0RecBg"), candDstarMcRec.ptD0());
        registry.fill(HIST("QA/hPtSkimDstarRecBg"), candDstarMcRec.pt());
      }
    }

    // MC at Generator Level
    for (const auto& mcParticle : rowsMcPartilces) {
      if (TESTBIT(std::abs(mcParticle.flagMcMatchGen()), aod::hf_cand_dstar::DecayType::DstarToD0Pi)) { // MC Matching is successful at Generator Level
        auto ptGen = mcParticle.pt();
        // auto yGen = mcParticle.y(); // Can we use this definition?
        auto yGen = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassDStar);
        if (yCandDstarGenMax >= 0. && std::abs(yGen) > yCandDstarGenMax) {
          continue;
        }
        registry.fill(HIST("QA/hEtaDstarGen"), mcParticle.eta());
        registry.fill(HIST("QA/hPtDstarGen"), ptGen);
        registry.fill(HIST("QA/hPtVsYDstarGen"), ptGen, yGen);
        // only promt Dstar candidate at Generator level
        if (mcParticle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("QA/hPtPromptDstarGen"), ptGen);
          registry.fill(HIST("QA/hPtVsYPromptDstarGen"), ptGen, yGen);
        } else if (mcParticle.originMcGen() == RecoDecay::OriginType::NonPrompt) { // only non-prompt Dstar candidate at Generator level
          registry.fill(HIST("QA/hPtNonPromptDstarGen"), ptGen);
          registry.fill(HIST("QA/hPtVsYNonPromptDstarGen"), ptGen, yGen);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDstarToD0Pi, processMC, "Process MC Data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDstarToD0Pi>(cfgc)};
}
