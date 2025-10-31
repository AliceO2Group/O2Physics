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

/// \file taskBplus.cxx
/// \brief B± analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Sezione di Bari
/// \author Deepa Thomas <deepa.thomas@cern.ch>, UT Austin

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TString.h>

#include <Rtypes.h>

#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_decay::hf_cand_beauty;

// string definitions, used for histogram axis labels
const TString stringPt = "#it{p}_{T} (GeV/#it{c})";
const TString bPlusCandTitle = "B+ candidates;";
const TString entries = "entries";
const TString bPlusCandMatch = "B+ candidates (matched);";
const TString bPlusCandUnmatch = "B+ candidates (unmatched);";
const TString mcParticleMatched = "MC particles (matched);";

/// B± analysis task
struct HfTaskBplus {
  Configurable<int> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for B+"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bplus_to_d0_pi::vecBinsPt}, "pT bin limits"};

  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg;
  HfHelper hfHelper;

  Partition<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>> selectedBPlusCandidates = aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBplus;
  Partition<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>> selectedBPlusCandidatesMC = aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBplus;

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
    // histograms MC: Gen Level
    registry.add("hEtaGen", mcParticleMatched + "candidate #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hYGen", mcParticleMatched + "candidate #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hPtProng0Gen", mcParticleMatched + "prong 0 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
    registry.add("hPtProng1Gen", mcParticleMatched + "prong 1 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
    registry.add("hYProng0Gen", mcParticleMatched + "prong 0 #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hYProng1Gen", mcParticleMatched + "prong 1 #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
    registry.add("hEtaProng0Gen", mcParticleMatched + "prong 0 #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hEtaProng1Gen", mcParticleMatched + "prong 1 #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
    registry.add("hPtGenWithRapidityBelowHalf", "MC particles (generated - |#it{y}^{gen}|<0.5);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    // histograms MC: Reco Level
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

  void process(aod::Collisions const&,
               soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi> const&,
               soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&,
               aod::Tracks const&)
  {

    for (const auto& candidate : selectedBPlusCandidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      auto ptCandBplus = candidate.pt();
      auto candD0 = candidate.prong0_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>();
      auto candPi = candidate.prong1();

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
      if (candPi.sign() > 0) {
        registry.fill(HIST("hInvMassD0"), hfHelper.invMassD0barToKPi(candD0), ptCandBplus);
      } else {
        registry.fill(HIST("hInvMassD0"), hfHelper.invMassD0ToPiK(candD0), ptCandBplus);
      }
      registry.fill(HIST("hCPAFinerBinning"), candidate.cpa(), ptCandBplus);
      registry.fill(HIST("hCPAxyFinerBinning"), candidate.cpaXY(), ptCandBplus);
    } // candidate loop
  } // process

  void processMc(soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec> const&,
                 soa::Join<aod::McParticles, aod::HfCandBplusMcGen> const& mcParticles,
                 aod::TracksWMc const&,
                 aod::HfCand2Prong const&)
  {
    // MC rec
    for (const auto& candidate : selectedBPlusCandidatesMC) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      auto ptCandBplus = candidate.pt();
      // auto candD0 = candidate.prong0_as<aod::HfCand2Prong>();
      if (std::abs(candidate.flagMcMatchRec()) == DecayChannelMain::BplusToD0Pi) {

        auto indexMother = RecoDecay::getMother(mcParticles, candidate.prong1_as<aod::TracksWMc>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandBplusMcGen>>(), o2::constants::physics::Pdg::kBPlus, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt());
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
      }
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == DecayChannelMain::BplusToD0Pi) {

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), o2::constants::physics::MassBPlus);
        if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
          continue;
        }

        float ptProngs[2], yProngs[2], etaProngs[2];
        int counter = 0;
        for (const auto& daught : particle.daughters_as<soa::Join<aod::McParticles, aod::HfCandBplusMcGen>>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }

        registry.fill(HIST("hPtProng0Gen"), ptProngs[0], ptParticle);
        registry.fill(HIST("hPtProng1Gen"), ptProngs[1], ptParticle);
        registry.fill(HIST("hYProng0Gen"), yProngs[0], ptParticle);
        registry.fill(HIST("hYProng1Gen"), yProngs[1], ptParticle);
        registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], ptParticle);
        registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], ptParticle);

        //  if (yCandMax >= 0. && (std::abs(yProngs[0]) > yCandMax || std::abs(yProngs[1]) > yCandMax))
        //    continue;

        registry.fill(HIST("hPtGen"), ptParticle);
        registry.fill(HIST("hYGen"), yParticle, ptParticle);
        registry.fill(HIST("hEtaGen"), particle.eta(), ptParticle);
        // generated B0 with |y|<0.5
        if (std::abs(yParticle) < 0.5) {
          registry.fill(HIST("hPtGenWithRapidityBelowHalf"), ptParticle);
        }

        // reject B+ daughters that are not in geometrical acceptance
        if (!isProngInAcceptance(etaProngs[0], ptProngs[0]) || !isProngInAcceptance(etaProngs[1], ptProngs[1])) {
          continue;
        }
        registry.fill(HIST("hPtGenWithProngsInAcceptance"), ptParticle);
        registry.fill(HIST("hYGenWithProngsInAcceptance"), yParticle, ptParticle);
        registry.fill(HIST("hEtaGenWithProngsInAcceptance"), particle.eta(), ptParticle);
      }
    } // gen
  } // processMc

  PROCESS_SWITCH(HfTaskBplus, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBplus>(cfgc)};
}
