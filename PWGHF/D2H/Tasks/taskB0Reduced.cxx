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

/// \file taskB0Reduced.cxx
/// \brief B0 → D- π+ → (π- K+ π-) π+ analysis task
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/RecoDecay.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// B0 analysis task
struct HfTaskB0Reduced {
  Configurable<int> selectionFlagB0{"selectionFlagB0", 1, "Selection Flag for B0"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_b0_to_d_pi::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_b0::isSelB0ToDPi >= selectionFlagB0);

  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "B0 candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hPtProng1", "B0 candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hPtCand", "B0 candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}}}};

  void init(InitContext&)
  {
    registry.add("hMass", "B^{0} candidates;inv. mass D^{#minus}#pi^{#plus} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{300, 4.5, 6.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "B^{0} candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXY", "B^{0} candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "B^{0} candidates;prong 0 (D^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "B^{0} candidates;prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "B^{0} candidates;B^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "B^{0} candidates;B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidity", "B^{0} candidates;B^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "B^{0} candidates;B^{0} candidate impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "B^{0} candidates;B^{0} candidate decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "B^{0} candidates;B^{0} candidate decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hIPProd", "B^{0} candidates;B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hInvMassD", "B^{0} candidates;prong0, D^{#minus} inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{350, 1.7, 2.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hEtaGen", "MC particles (generated);B^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGen", "MC particles (generated);B^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0Gen", "MC particles (generated);prong 0 (D^{#minus}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1Gen", "MC particles (generated);prong 1 (#pi^{-}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYProng0Gen", "MC particles (generated);prong 0 (D^{#minus}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYProng1Gen", "MC particles (generated);prong 1 (#pi^{-}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaProng0Gen", "MC particles (generated);prong 0 (B^{0}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaProng1Gen", "MC particles (generated);prong 1 (#pi^{-}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecSig", "B^{0} candidates (matched);B^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "B^{0} candidates (unmatched);B^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyRecSig", "B^{0} candidates (matched);B^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyRecBg", "B^{0} candidates (unmatched);B^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPADRecSig", "B^{0} candidates (matched);prong 0 (D^{#minus}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPADRecBg", "B^{0} candidates (unmatched);prong 0 (D^{#minus}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "B^{0} candidates (matched);B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "B^{0} candidates (unmatched);B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidityRecSig", "B^{0} candidates (matched);B^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidityRecBg", "B^{0} candidates (unmatched);B^{0} candidate #it{#y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtProng0RecSig", "B^{0} candidates (matched);prong 0 (D^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "B^{0} candidates (matched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "B^{0} candidates (unmatched);prong 0 (D^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "B^{0} candidates (unmatched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecSig", "B^{0} candidates (matched);inv. mass D^{#minus}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.00}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecBg", "B^{0} candidates (unmatched);inv. mass D^{#minus}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "B^{0} candidates (matched);prong 0 (D^{#minus}}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "B^{0} candidates (matched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "B^{0} candidates (unmatched);prong 0 (D^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "B^{0} candidates (unmatched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecSig", "B^{0} candidates (matched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYRecSig", "B^{0} candidates (matched);B^{0} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecBg", "B^{0} candidates (unmatched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYRecBg", "B^{0} candidates (unmatched);B^{0} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthDRecSig", "B^{0} candidates (matched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthDRecBg", "B^{0} candidates (unmatched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthNormRecSig", "B^{0} candidates (matched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthNormRecBg", "B^{0} candidates (unmatched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParProdB0RecSig", "B^{0} candidates (matched);B^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.01, 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParProdB0RecBg", "B^{0} candidates (unmatched);B^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.01, 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hChi2PCARecSig", "B^{0} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCARecBg", "B^{0} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtRecSig", "B0 candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtRecBg", "B0 candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenSig", "B0 candidates (gen+rec);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 10.}}});
    registry.add("hPtGen", "MC particles (generated);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenWithRapidityBelowHalf", "MC particles (generated - |#it{y}^{gen}|<0.5);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
  }

  /// Selection of B0 daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of B0 prong
  /// \param ptProng is the pT of B0 prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  void process(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfSelB0ToDPi>> const& candidates,
               aod::HfRed3Prongs const&)
  {
    for (const auto& candidate : candidates) {
      if (!TESTBIT(candidate.hfflag(), hf_cand_b0::DecayType::B0ToDPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }

      auto ptCandB0 = candidate.pt();
      auto candD = candidate.prong0_as<aod::HfRed3Prongs>();

      registry.fill(HIST("hMass"), hfHelper.invMassB0ToDPi(candidate), ptCandB0);
      registry.fill(HIST("hPtCand"), ptCandB0);
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hIPProd"), candidate.impactParameterProduct(), ptCandB0);
      registry.fill(HIST("hDecLength"), candidate.decayLength(), ptCandB0);
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), ptCandB0);
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), ptCandB0);
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), ptCandB0);
      registry.fill(HIST("hCPA"), candidate.cpa(), ptCandB0);
      registry.fill(HIST("hEta"), candidate.eta(), ptCandB0);
      registry.fill(HIST("hRapidity"), hfHelper.yB0(candidate), ptCandB0);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), ptCandB0);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), ptCandB0);
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), ptCandB0);
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), ptCandB0);
      registry.fill(HIST("hInvMassD"), candD.invMass(), ptCandB0);
    } // candidate loop
  }   // process

  /// B0 MC analysis and fill histograms
  void processMc(soa::Join<aod::HfRedCandB0, aod::HfMcRecRedB0s> const& candidates,
                 aod::HfMcGenRedB0s const& mcParticles,
                 aod::HfRed3Prongs const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (!TESTBIT(candidate.hfflag(), hf_cand_b0::DecayType::B0ToDPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }

      auto ptCandB0 = candidate.pt();
      auto candD = candidate.prong0_as<aod::HfRed3Prongs>();
      std::array<float, 3> posPv{candidate.posX(), candidate.posY(), candidate.posZ()};
      std::array<float, 3> posSvD{candD.xSecondaryVertex(), candD.ySecondaryVertex(), candD.zSecondaryVertex()};
      std::array<float, 3> momD{candD.px(), candD.py(), candD.pz()};
      auto cospD = RecoDecay::cpa(posPv, posSvD, momD);
      auto decLenD = RecoDecay::distance(posPv, posSvD);

      if (TESTBIT(std::abs(candidate.flagMcMatchRec()), hf_cand_b0::DecayType::B0ToDPi)) {
        registry.fill(HIST("hPtGenSig"), candidate.ptMother());
        registry.fill(HIST("hPtRecSig"), ptCandB0);
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), ptCandB0);
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpaXY(), ptCandB0);
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), ptCandB0);
        registry.fill(HIST("hRapidityRecSig"), hfHelper.yB0(candidate), ptCandB0);
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), ptCandB0);
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), ptCandB0);
        registry.fill(HIST("hMassRecSig"), hfHelper.invMassB0ToDPi(candidate), ptCandB0);
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), ptCandB0);
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), ptCandB0);
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), ptCandB0);
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), ptCandB0);
        registry.fill(HIST("hImpParProdB0RecSig"), candidate.impactParameterProduct(), ptCandB0);
        registry.fill(HIST("hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), ptCandB0);
        registry.fill(HIST("hCPADRecSig"), cospD, ptCandB0);
        registry.fill(HIST("hDecLengthDRecSig"), decLenD, ptCandB0);
        registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), ptCandB0);
      } else {
        registry.fill(HIST("hPtRecBg"), ptCandB0);
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), ptCandB0);
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpaXY(), ptCandB0);
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), ptCandB0);
        registry.fill(HIST("hRapidityRecBg"), hfHelper.yB0(candidate), ptCandB0);
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), ptCandB0);
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), ptCandB0);
        registry.fill(HIST("hMassRecBg"), hfHelper.invMassB0ToDPi(candidate), ptCandB0);
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), ptCandB0);
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), ptCandB0);
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), ptCandB0);
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), ptCandB0);
        registry.fill(HIST("hImpParProdB0RecBg"), candidate.impactParameterProduct(), ptCandB0);
        registry.fill(HIST("hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), ptCandB0);
        registry.fill(HIST("hCPADRecBg"), cospD, ptCandB0);
        registry.fill(HIST("hDecLengthDRecBg"), decLenD, ptCandB0);
        registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), ptCandB0);
      }
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      auto ptParticle = particle.ptTrack();
      auto yParticle = particle.yTrack();
      auto etaParticle = particle.etaTrack();
      if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
        continue;
      }

      std::array<float, 2> ptProngs = {particle.ptProng0(), particle.ptProng1()};
      std::array<float, 2> yProngs = {particle.yProng0(), particle.yProng1()};
      std::array<float, 2> etaProngs = {particle.etaProng0(), particle.etaProng1()};

      registry.fill(HIST("hPtProng0Gen"), ptProngs[0], ptParticle);
      registry.fill(HIST("hPtProng1Gen"), ptProngs[1], ptParticle);
      registry.fill(HIST("hYProng0Gen"), yProngs[0], ptParticle);
      registry.fill(HIST("hYProng1Gen"), yProngs[1], ptParticle);
      registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], ptParticle);
      registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], ptParticle);

      registry.fill(HIST("hPtGen"), ptParticle);
      registry.fill(HIST("hYGen"), yParticle, ptParticle);
      registry.fill(HIST("hEtaGen"), etaParticle, ptParticle);

      // generated B0 with |y|<0.5
      if (std::abs(yParticle) < 0.5) {
        registry.fill(HIST("hPtGenWithRapidityBelowHalf"), ptParticle);
      }

      // generated B0 with daughters in geometrical acceptance
      if (isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1])) {
        registry.fill(HIST("hPtGenWithProngsInAcceptance"), ptParticle);
        registry.fill(HIST("hYGenWithProngsInAcceptance"), yParticle, ptParticle);
        registry.fill(HIST("hEtaGenWithProngsInAcceptance"), etaParticle, ptParticle);
      }
    } // gen
  }   // process
  PROCESS_SWITCH(HfTaskB0Reduced, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskB0Reduced>(cfgc)};
}
