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

/// \file taskBplusReduced.cxx
/// \brief B+ → D0bar π+ → (π+ K-) π+ analysis task
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Sezione di Bari

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

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

// string definitions, used for histogram axis labels
const TString stringPt = "#it{p}_{T} (GeV/#it{c})";
const TString bPlusCandTitle = "B+ candidates;";
const TString entries = "entries";
const TString bPlusCandMatch = "B+ candidates (matched);";
const TString bPlusCandUnmatch = "B+ candidates (unmatched);";
const TString mcParticleMatched = "MC particles (matched);";

/// B+ analysis task
struct HfTaskBplusReduced {
  Configurable<int> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for Bplus"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bplus_to_d0_pi::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBplus);

  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", bPlusCandTitle + "prong 0 #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hPtProng1", bPlusCandTitle + "prong 1 #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hPtCand", bPlusCandTitle + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hCentrality", "centrality;centrality percentile;" + entries, {HistType::kTH1F, {{100, 0., 100.}}}},
     {"hPtRecSig", bPlusCandMatch + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtRecBg", bPlusCandUnmatch + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtGenSig", bPlusCandMatch + "candidate #it{p}_{T}^{gen.} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtGen", mcParticleMatched + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}}}}};

  void init(InitContext&)
  {
    const AxisSpec axisMass{150, 4.5, 6.0};
    const AxisSpec axisCPA{120, -1.1, 1.1};
    const AxisSpec axisCPAFiner{300, 0.85, 1.0};
    const AxisSpec axisPtProng{100, 0., 10.};
    const AxisSpec axisD0Prong{200, -0.05, 0.05};
    const AxisSpec axisImpParProd{200, -0.001, 0.001};
    const AxisSpec axisDecLength{100, 0., 0.5};
    const AxisSpec axisNormDecLength{40, 0., 20};
    const AxisSpec axisEta{100, -2., 2.};
    const AxisSpec axisRapidity{100, -2., 2.};
    const AxisSpec axisPtB{(std::vector<double>)binsPt, "#it{p}_{T}^{B^{+}} (GeV/#it{c})"};
    // histograms process
    registry.add("hMass", bPlusCandTitle + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMass, axisPtB}});
    registry.add("hDecLength", bPlusCandTitle + "decay length (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
    registry.add("hDecLengthXY", bPlusCandTitle + "decay length xy (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
    registry.add("hd0Prong0", bPlusCandTitle + "prong 0 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
    registry.add("hd0Prong1", bPlusCandTitle + "prong 1 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
    registry.add("hCPA", bPlusCandTitle + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCPA, axisPtB}});
    registry.add("hCPAxy", bPlusCandTitle + "candidate cosine of pointing angle xy;" + stringPt, {HistType::kTH2F, {axisCPA, axisPtB}});
    registry.add("hEta", bPlusCandTitle + "candidate #it{#eta};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hRapidity", bPlusCandTitle + "candidate #it{y};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hImpParErr", bPlusCandTitle + "candidate impact parameter error (cm);" + stringPt, {HistType::kTH2F, {{100, -1., 1.}, axisPtB}});
    registry.add("hDecLenErr", bPlusCandTitle + "candidate decay length error (cm);" + stringPt, {HistType::kTH2F, {{100, 0., 1.}, axisPtB}});
    registry.add("hDecLenXYErr", bPlusCandTitle + "candidate decay length xy error (cm);" + stringPt, {HistType::kTH2F, {{100, 0., 1.}, axisPtB}});
    registry.add("hd0d0", bPlusCandTitle + "candidate product of DCAxy to prim. vertex (cm^{2});" + stringPt, {HistType::kTH2F, {axisImpParProd, axisPtB}});
    registry.add("hInvMassD0", bPlusCandTitle + "prong0, D0 inv. mass (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {{500, 1.4, 2.4}, axisPtB}});
    registry.add("hCPAFinerBinning", bPlusCandTitle + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCPAFiner, axisPtB}});
    registry.add("hCPAxyFinerBinning", bPlusCandTitle + "candidate cosine of pointing angle xy;" + stringPt, {HistType::kTH2F, {axisCPAFiner, axisPtB}});
    // histograms processMC - Gen Level
    registry.add("hEtaGen", mcParticleMatched + "candidate #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hYGen", mcParticleMatched + "candidate #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hPtProng0Gen", mcParticleMatched + "prong 0 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
    registry.add("hPtProng1Gen", mcParticleMatched + "prong 1 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
    registry.add("hYProng0Gen", mcParticleMatched + "prong 0 #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hYProng1Gen", mcParticleMatched + "prong 1 #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hEtaProng0Gen", mcParticleMatched + "prong 0 #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hEtaProng1Gen", mcParticleMatched + "prong 1 #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hPtProngsVsPtBGen", mcParticleMatched + "prong 0 #it{p}_{T}^{gen} (GeV/#it{c});prong 1 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH3F, {axisPtProng, axisPtProng, axisPtB}});
    registry.add("hYProngsVsBplusGen", mcParticleMatched + "prong 0 #it{y}^{gen};prong 1 #it{y}^{gen};" + stringPt, {HistType::kTH3F, {axisRapidity, axisRapidity, axisPtB}});
    registry.add("hEtaProngsVsBplusGen", mcParticleMatched + "prong 0 #it{#eta}^{gen};prong 1 #it{#eta}^{gen};" + stringPt, {HistType::kTH3F, {axisEta, axisEta, axisPtB}});
    registry.add("hPtGenWithRapidityBelowHalf", "MC particles (generated - |#it{y}^{gen}|<0.5);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    // histograms processMC - Reco Level
    registry.add("hCPARecSig", bPlusCandMatch + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCPA, axisPtB}});
    registry.add("hCPARecBg", bPlusCandUnmatch + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCPA, axisPtB}});
    registry.add("hCPAxyRecSig", bPlusCandMatch + "candidate CPAxy;" + stringPt, {HistType::kTH2F, {axisCPA, axisPtB}});
    registry.add("hCPAxyRecBg", bPlusCandUnmatch + "candidate CPAxy;" + stringPt, {HistType::kTH2F, {axisCPA, axisPtB}});
    registry.add("hEtaRecSig", bPlusCandMatch + "candidate #it{#eta};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hEtaRecBg", bPlusCandUnmatch + "candidate #it{#eta};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hRapidityRecSig", bPlusCandMatch + "candidate #it{y};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hRapidityRecBg", bPlusCandUnmatch + "candidate #it{#y};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hPtProng0RecSig", bPlusCandMatch + "prong 0 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
    registry.add("hPtProng1RecSig", bPlusCandMatch + "prong 1 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
    registry.add("hPtProngsVsBplusRecBg", bPlusCandUnmatch + "prong 0 #it{p}_{T} (GeV/#it{c});prong 1 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH3F, {axisPtProng, axisPtProng, axisPtB}});
    registry.add("hPtProng0RecBg", bPlusCandUnmatch + "prong 0 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
    registry.add("hPtProng1RecBg", bPlusCandUnmatch + "prong 1 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
    registry.add("hMassRecSig", bPlusCandMatch + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMass, axisPtB}});
    registry.add("hMassRecBg", bPlusCandUnmatch + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMass, axisPtB}});
    registry.add("hd0Prong0RecSig", bPlusCandMatch + "prong 0 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
    registry.add("hd0Prong1RecSig", bPlusCandMatch + "prong 1 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
    registry.add("hd0Prong0RecBg", bPlusCandUnmatch + "prong 0 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
    registry.add("hd0Prong1RecBg", bPlusCandUnmatch + "prong 1 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
    registry.add("hDecLengthRecSig", bPlusCandMatch + "candidate decay length (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
    registry.add("hDecLengthXYRecSig", bPlusCandMatch + "candidate decay length xy (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
    registry.add("hDecLengthRecBg", bPlusCandUnmatch + "candidate decay length (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
    registry.add("hDecLengthXYRecBg", bPlusCandUnmatch + "candidate decay length xy(cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
    registry.add("hDecLengthNormRecSig", bPlusCandMatch + "candidate normalized decay length (cm);" + stringPt, {HistType::kTH2F, {axisNormDecLength, axisPtB}});
    registry.add("hDecLengthNormRecBg", bPlusCandUnmatch + "candidate normalized decay length (cm);" + stringPt, {HistType::kTH2F, {axisNormDecLength, axisPtB}});
    registry.add("hd0d0RecSig", bPlusCandMatch + "product of DCAxy to prim. vertex (cm^{2});" + stringPt, {HistType::kTH2F, {axisImpParProd, axisPtB}});
    registry.add("hd0d0RecBg", bPlusCandUnmatch + "product of DCAxy to prim. vertex (cm^{2});" + stringPt, {HistType::kTH2F, {axisImpParProd, axisPtB}});
    // MC histograms with finer binning
    registry.add("hCPAFinerBinningRecSig", bPlusCandMatch + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCPAFiner, axisPtB}});
    registry.add("hCPAFinerBinningRecBg", bPlusCandMatch + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCPAFiner, axisPtB}});
    registry.add("hCPAxyFinerBinningRecSig", bPlusCandMatch + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCPAFiner, axisPtB}});
    registry.add("hCPAxyFinerBinningRecBg", bPlusCandMatch + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCPAFiner, axisPtB}});
    // histograms prong0(D0) - Reco Level
    registry.add("hCPAD0RecSig", bPlusCandMatch + "prong0 (D^{0}) cosine of pointing angle;#it{p}_{T}(D0) (GeV/#it{c})" + entries, {HistType::kTH2F, {{220, 0., 1.1}, {120, 0., 60.}}});
    registry.add("hCPAD0RecBg", bPlusCandUnmatch + "prong0 (D^{0}) cosine of pointing angle;#it{p}_{T}(D0) (GeV/#it{c})" + entries, {HistType::kTH2F, {{220, 0., 1.1}, {120, 0., 60.}}});
    registry.add("hDecLengthD0RecSig", bPlusCandMatch + "prong0 D^{0} decay length (cm);#it{p}_{T}(D0) (GeV/#it{c})" + entries, {HistType::kTH2F, {{100, 0., 0.5}, {120, 0., 60.}}});
    registry.add("hDecLengthD0RecBg", bPlusCandUnmatch + "prong0 D^{0} candidate decay length (cm);#it{p}_{T}(D0) (GeV/#it{c})" + entries, {HistType::kTH2F, {{100, 0., 0.5}, {120, 0., 60.}}});
  }

  /// Selection of B+ daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of B+ prong
  /// \param ptProng is the pT of B+ prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  void process(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfSelBplusToD0Pi>> const& candidates,
               aod::HfRed2Prongs const&,
               aod::HfRedTracks const&)
  {
    for (const auto& candidate : candidates) {
      if (!TESTBIT(candidate.hfflag(), hf_cand_bplus::DecayType::BplusToD0Pi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBplus(candidate)) > yCandRecoMax) {
        continue;
      }

      auto candD0 = candidate.prong0_as<HfRed2Prongs>();
      auto candPi = candidate.prong1_as<aod::HfRedTracks>();
      auto ptCandBplus = candidate.pt();
      auto invMassD0 = (candPi.signed1Pt() < 0) ? candD0.invMassD0() : candD0.invMassD0Bar();

      registry.fill(HIST("hMass"), hfHelper.invMassBplusToD0Pi(candidate), ptCandBplus);
      registry.fill(HIST("hPtCand"), ptCandBplus);
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), ptCandBplus);
      registry.fill(HIST("hDecLength"), candidate.decayLength(), ptCandBplus);
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), ptCandBplus);
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), ptCandBplus);
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), ptCandBplus);
      registry.fill(HIST("hCPA"), candidate.cpa(), ptCandBplus);
      registry.fill(HIST("hCPAxy"), candidate.cpaXY(), ptCandBplus);
      registry.fill(HIST("hEta"), candidate.eta(), ptCandBplus);
      registry.fill(HIST("hRapidity"), hfHelper.yBplus(candidate), ptCandBplus);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), ptCandBplus);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), ptCandBplus);
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), ptCandBplus);
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), ptCandBplus);
      registry.fill(HIST("hInvMassD0"), invMassD0, ptCandBplus);
      registry.fill(HIST("hCPAFinerBinning"), candidate.cpa(), ptCandBplus);
      registry.fill(HIST("hCPAxyFinerBinning"), candidate.cpaXY(), ptCandBplus);
    } // candidate loop
  }   // process

  /// B+ MC analysis and fill histograms
  void processMc(soa::Join<aod::HfRedCandBplus, aod::HfMcRecRedBps> const& candidates,
                 aod::HfMcGenRedBps const& mcParticles,
                 aod::HfRed2Prongs const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (!TESTBIT(candidate.hfflag(), hf_cand_bplus::DecayType::BplusToD0Pi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBplus(candidate)) > yCandRecoMax) {
        continue;
      }

      auto ptCandBplus = candidate.pt();
      auto candD0 = candidate.prong0_as<aod::HfRed2Prongs>();
      std::array<float, 3> posPv{candidate.posX(), candidate.posY(), candidate.posZ()};
      std::array<float, 3> posSvD{candD0.xSecondaryVertex(), candD0.ySecondaryVertex(), candD0.zSecondaryVertex()};
      std::array<float, 3> momD{candD0.px(), candD0.py(), candD0.pz()};
      auto cospD0 = RecoDecay::cpa(posPv, posSvD, momD);
      auto decLenD0 = RecoDecay::distance(posPv, posSvD);

      if (TESTBIT(std::abs(candidate.flagMcMatchRec()), hf_cand_bplus::DecayType::BplusToD0Pi)) {

        registry.fill(HIST("hPtGenSig"), candidate.ptMother());
        registry.fill(HIST("hPtRecSig"), ptCandBplus);
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), ptCandBplus);
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpaXY(), ptCandBplus);
        registry.fill(HIST("hCPAFinerBinningRecSig"), candidate.cpa(), ptCandBplus);
        registry.fill(HIST("hCPAxyFinerBinningRecSig"), candidate.cpaXY(), ptCandBplus);
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), ptCandBplus);
        registry.fill(HIST("hRapidityRecSig"), hfHelper.yBplus(candidate), ptCandBplus);
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), ptCandBplus);
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), ptCandBplus);
        registry.fill(HIST("hMassRecSig"), hfHelper.invMassBplusToD0Pi(candidate), ptCandBplus);
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), ptCandBplus);
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), ptCandBplus);
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), ptCandBplus);
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), ptCandBplus);
        registry.fill(HIST("hd0d0RecSig"), candidate.impactParameterProduct(), ptCandBplus);
        registry.fill(HIST("hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), ptCandBplus);
        registry.fill(HIST("hCPAD0RecSig"), cospD0, candidate.ptProng0());
        registry.fill(HIST("hDecLengthD0RecSig"), decLenD0, candidate.ptProng0());
      } else {
        registry.fill(HIST("hPtRecBg"), ptCandBplus);
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), ptCandBplus);
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpaXY(), ptCandBplus);
        registry.fill(HIST("hCPAFinerBinningRecBg"), candidate.cpa(), ptCandBplus);
        registry.fill(HIST("hCPAxyFinerBinningRecBg"), candidate.cpaXY(), ptCandBplus);
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), ptCandBplus);
        registry.fill(HIST("hRapidityRecBg"), hfHelper.yBplus(candidate), ptCandBplus);
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), ptCandBplus);
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), ptCandBplus);
        registry.fill(HIST("hMassRecBg"), hfHelper.invMassBplusToD0Pi(candidate), ptCandBplus);
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), ptCandBplus);
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), ptCandBplus);
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), ptCandBplus);
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), ptCandBplus);
        registry.fill(HIST("hd0d0RecBg"), candidate.impactParameterProduct(), ptCandBplus);
        registry.fill(HIST("hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), ptCandBplus);
        registry.fill(HIST("hCPAD0RecBg"), cospD0, candidate.ptProng0());
        registry.fill(HIST("hDecLengthD0RecBg"), decLenD0, candidate.ptProng0());
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
      registry.fill(HIST("hPtProngsVsPtBGen"), ptProngs[0], ptProngs[1], ptParticle);
      registry.fill(HIST("hYProng0Gen"), yProngs[0], ptParticle);
      registry.fill(HIST("hYProng1Gen"), yProngs[1], ptParticle);
      registry.fill(HIST("hYProngsVsBplusGen"), yProngs[0], yProngs[1], ptParticle);
      registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], ptParticle);
      registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], ptParticle);
      registry.fill(HIST("hEtaProngsVsBplusGen"), etaProngs[0], etaProngs[1], ptParticle);

      registry.fill(HIST("hPtGen"), ptParticle);
      registry.fill(HIST("hYGen"), yParticle, ptParticle);
      registry.fill(HIST("hEtaGen"), etaParticle, ptParticle);

      // generated B+ with |y|<0.5
      if (std::abs(yParticle) < 0.5) {
        registry.fill(HIST("hPtGenWithRapidityBelowHalf"), ptParticle);
      }

      // generated B+ with daughters in geometrical acceptance
      if (isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1])) {
        registry.fill(HIST("hPtGenWithProngsInAcceptance"), ptParticle);
        registry.fill(HIST("hYGenWithProngsInAcceptance"), yParticle, ptParticle);
        registry.fill(HIST("hEtaGenWithProngsInAcceptance"), etaParticle, ptParticle);
      }
    } // gen
  }   // process
  PROCESS_SWITCH(HfTaskBplusReduced, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBplusReduced>(cfgc)};
}
