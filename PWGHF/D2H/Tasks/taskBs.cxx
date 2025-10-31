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

/// \file taskBs.cxx
/// \brief Bs → Ds π+ → (K- K+ π-) π+ analysis task
/// \note adapted from taskB0.cxx
///
/// \author Phil Stahlhut <phil.lennart.stahlhut@cern.ch>

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

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

#include <TH3.h>
#include <TString.h>

#include <Rtypes.h>

#include <array>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_decay::hf_cand_beauty;

/// Bs analysis task
struct HfTaskBs {
  Configurable<int> selectionFlagBs{"selectionFlagBs", 1, "Selection Flag for Bs"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bs_to_ds_pi::vecBinsPt}, "pT bin limits"};
  // MC checks
  Configurable<bool> checkDecayTypeMc{"checkDecayTypeMc", false, "Flag to enable DecayType histogram"};

  Service<o2::framework::O2DatabasePDG> pdg;
  HfHelper hfHelper;

  using TracksWithSel = soa::Join<aod::Tracks, aod::TrackSelection>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_bs::isSelBsToDsPi >= selectionFlagBs);

  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "B^{0}_{s} candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 40.}}}},
     {"hPtProng1", "B^{0}_{s} candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 16.}}}},
     {"hPtCand", "B^{0}_{s} candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 40.}}}}}};

  void init(InitContext const&)
  {
    static const AxisSpec axisMassBs = {300, 4.5, 6.0, "inv. mass (GeV/#it{c}^{2})"};
    static const AxisSpec axisPt = {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"};

    registry.add("hEta", "B^{0}_{s} candidates;B^{0}_{s} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hRapidity", "B^{0}_{s} candidates;B^{0}_{s} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hCPA", "B^{0}_{s} candidates;B^{0}_{s} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
    registry.add("hMass", "B^{0}_{s} candidates;inv. mass D^{#mp}_{s} #pi^{#pm} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisMassBs, axisPt}});
    registry.add("hDecLength", "B^{0}_{s} candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, axisPt}});
    registry.add("hDecLenErr", "B^{0}_{s} candidates;B^{0}_{s} candidate decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, axisPt}});
    registry.add("hDecLengthXY", "B^{0}_{s} candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, axisPt}});
    registry.add("hDecLenXYErr", "B^{0}_{s} candidates;B^{0}_{s} candidate decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, axisPt}});
    registry.add("hd0Prong0", "B^{0}_{s} candidates;prong 0 (D^{#pm}_{s}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong1", "B^{0}_{s} candidates;prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, axisPt}});
    registry.add("hImpParErr", "B^{0}_{s} candidates;B^{0}_{s} candidate impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, axisPt}});
    registry.add("hIPProd", "B^{0}_{s} candidates;B^{0}_{s} candidate impact parameter product;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, axisPt}});
    registry.add("hInvMassDs", "B^{0}_{s} candidates;prong 0 (D^{#pm}_{s}) inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0, 5}, axisPt}});

    registry.add("hPtGenSig", "B^{0}_{s} candidates (gen+rec);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtRecSig", "B^{0}_{s} candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtRecBg", "B^{0}_{s} candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hEtaRecSig", "B^{0}_{s} candidates (matched);B^{0}_{s} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hEtaRecBg", "B^{0}_{s} candidates (unmatched);B^{0}_{s} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hRapidityRecSig", "B^{0}_{s} candidates (matched);B^{0}_{s} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hRapidityRecBg", "B^{0}_{s} candidates (unmatched);B^{0}_{s} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hCPARecSig", "B^{0}_{s} candidates (matched);B^{0}_{s} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPARecBg", "B^{0}_{s} candidates (unmatched);B^{0}_{s} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPAxyRecSig", "B^{0}_{s} candidates (matched);B^{0}_{s} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPAxyRecBg", "B^{0}_{s} candidates (unmatched);B^{0}_{s} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hMassRecSig", "B^{0}_{s} candidates (matched);inv. mass D^{#mp}_{s} #pi^{#pm} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.00}, axisPt}});
    registry.add("hMassRecBg", "B^{0}_{s} candidates (unmatched);inv. mass D^{#mp}_{s} #pi^{#pm} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.0}, axisPt}});
    registry.add("hDecLengthRecSig", "B^{0}_{s} candidates (matched);B^{0}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthRecBg", "B^{0}_{s} candidates (unmatched);B^{0}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthXYRecSig", "B^{0}_{s} candidates (matched);B^{0}_{s} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthXYRecBg", "B^{0}_{s} candidates (unmatched);B^{0}_{s} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthNormRecSig", "B^{0}_{s} candidates (matched);B^{0}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthNormRecBg", "B^{0}_{s} candidates (unmatched);B^{0}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hImpParProdBsRecSig", "B^{0}_{s} candidates (matched);B^{0}_{s} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.01, 0.01}, axisPt}});
    registry.add("hImpParProdBsRecBg", "B^{0}_{s} candidates (unmatched);B^{0}_{s} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.01, 0.01}, axisPt}});
    registry.add("hPtProng0RecSig", "B^{0}_{s} candidates (matched);prong 0 (D^{#pm}_{s}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hPtProng0RecBg", "B^{0}_{s} candidates (unmatched);prong 0 (D^{#pm}_{s}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hPtProng1RecSig", "B^{0}_{s} candidates (matched);prong 1 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hPtProng1RecBg", "B^{0}_{s} candidates (unmatched);prong 1 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hd0Prong0RecSig", "B^{0}_{s} candidates (matched);prong 0 (D^{#pm}_{s}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong0RecBg", "B^{0}_{s} candidates (unmatched);prong 0 (D^{#pm}_{s}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong1RecSig", "B^{0}_{s} candidates (matched);prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong1RecBg", "B^{0}_{s} candidates (unmatched);prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
    registry.add("hCPADsRecSig", "B^{0}_{s} candidates (matched);prong 0 (D^{#pm}_{s}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPADsRecBg", "B^{0}_{s} candidates (unmatched);prong 0 (D^{#pm}_{s}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hDecLengthDsRecSig", "B^{0}_{s} candidates (matched);D^{#pm}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthDsRecBg", "B^{0}_{s} candidates (unmatched);D^{#pm}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hChi2PCARecSig", "B^{0}_{s} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, axisPt}});
    registry.add("hChi2PCARecBg", "B^{0}_{s} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, axisPt}});

    registry.add("hPtProng0Gen", "MC particles (generated);prong 0 (D_{s}^{#pm}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hPtProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hEtaProng0Gen", "MC particles (generated);prong 0 (D^{#pm}_{s}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
    registry.add("hEtaProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
    registry.add("hYProng0Gen", "MC particles (generated);prong 0 (D^{#pm}_{s}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
    registry.add("hYProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
    registry.add("hPtGen", "MC particles (generated);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hEtaGen", "MC particles (generated);B^{0}_{s} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hYGen", "MC particles (generated);B^{0}_{s} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hPtGenWithRapidityBelowHalf", "MC particles (generated - |#it{y}^{gen}|<0.5);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0}_{s} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0}_{s} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});

    if (checkDecayTypeMc) {
      constexpr uint8_t NBinsDecayTypeMc = hf_cand_bs::DecayTypeMc::NDecayTypeMc + 1;
      TString labels[NBinsDecayTypeMc];
      labels[hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi] = "B^{0}_{s} #rightarrow (D^{#mp}_{s} #rightarrow K^{#minus} K^{#plus} #pi^{#mp}) #pi^{#pm}";
      labels[hf_cand_bs::DecayTypeMc::B0ToDsPiToPhiPiPiToKKPiPi] = "B^{0} #rightarrow (D^{#pm}_{s} #rightarrow K^{#minus} K^{#plus} #pi^{#pm}) #pi^{#mp}";
      labels[hf_cand_bs::DecayTypeMc::PartlyRecoDecay] = "Partly reconstructed decay channel";
      labels[hf_cand_bs::DecayTypeMc::NDecayTypeMc] = "Other decays";
      static const AxisSpec axisDecayType = {NBinsDecayTypeMc, 0.5, NBinsDecayTypeMc + 0.5, ""};
      registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassBs, axisPt}});
      for (uint8_t iBin = 0; iBin < NBinsDecayTypeMc; ++iBin) {
        registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
      }
    }
  }

  /// Selection of Bs daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of Bs prong
  /// \param ptProng is the pT of Bs prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  void process(soa::Filtered<soa::Join<aod::HfCandBs, aod::HfSelBsToDsPi>> const& candidates,
               soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi> const&,
               TracksWithSel const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }

      auto candDs = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>();

      auto ptCandBs = candidate.pt();

      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hPtCand"), ptCandBs);
      registry.fill(HIST("hEta"), candidate.eta(), ptCandBs);
      registry.fill(HIST("hRapidity"), hfHelper.yBs(candidate), ptCandBs);
      registry.fill(HIST("hCPA"), candidate.cpa(), ptCandBs);
      registry.fill(HIST("hMass"), hfHelper.invMassBsToDsPi(candidate), ptCandBs);
      registry.fill(HIST("hDecLength"), candidate.decayLength(), ptCandBs);
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), ptCandBs);
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), ptCandBs);
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), ptCandBs);
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), ptCandBs);
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), ptCandBs);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), ptCandBs);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), ptCandBs);
      registry.fill(HIST("hIPProd"), candidate.impactParameterProduct(), ptCandBs);
      registry.fill(HIST("hInvMassDs"), hfHelper.invMassDsToKKPi(candDs), ptCandBs);
    } // candidate loop
  } // process

  /// Bs MC analysis and fill histograms
  void processMc(soa::Filtered<soa::Join<aod::HfCandBs, aod::HfSelBsToDsPi, aod::HfCandBsMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandBsMcGen> const& mcParticles,
                 aod::TracksWMc const&,
                 soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec> const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }

      auto ptCandBs = candidate.pt();
      auto candDs = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec>>();
      auto invMassCandBs = hfHelper.invMassBsToDsPi(candidate);
      auto flagMcMatchRecBs = std::abs(candidate.flagMcMatchRec());

      if (flagMcMatchRecBs == DecayChannelMain::BsToDsPi) {
        auto indexMother = RecoDecay::getMother(mcParticles, candidate.prong1_as<aod::TracksWMc>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandBsMcGen>>(), o2::constants::physics::Pdg::kBS, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);

        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), ptCandBs);
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), ptCandBs);
        registry.fill(HIST("hRapidityRecSig"), hfHelper.yBs(candidate), ptCandBs);
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), ptCandBs);
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpaXY(), ptCandBs);
        registry.fill(HIST("hMassRecSig"), hfHelper.invMassBsToDsPi(candidate), ptCandBs);
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), ptCandBs);
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), ptCandBs);
        registry.fill(HIST("hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), ptCandBs);
        registry.fill(HIST("hImpParProdBsRecSig"), candidate.impactParameterProduct(), ptCandBs);
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), ptCandBs);
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), ptCandBs);
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), ptCandBs);
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), ptCandBs);
        registry.fill(HIST("hCPADsRecSig"), candDs.cpa(), ptCandBs);
        registry.fill(HIST("hDecLengthDsRecSig"), candDs.decayLength(), ptCandBs);
        registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), ptCandBs);

        if (checkDecayTypeMc) {
          registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi, invMassCandBs, ptCandBs);
        }
      } else {
        registry.fill(HIST("hPtRecBg"), ptCandBs);
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), ptCandBs);
        registry.fill(HIST("hRapidityRecBg"), hfHelper.yBs(candidate), ptCandBs);
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), ptCandBs);
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpaXY(), ptCandBs);
        registry.fill(HIST("hMassRecBg"), hfHelper.invMassBsToDsPi(candidate), ptCandBs);
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), ptCandBs);
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), ptCandBs);
        registry.fill(HIST("hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), ptCandBs);
        registry.fill(HIST("hImpParProdBsRecBg"), candidate.impactParameterProduct(), ptCandBs);
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), ptCandBs);
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), ptCandBs);
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), ptCandBs);
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), ptCandBs);
        registry.fill(HIST("hCPADsRecBg"), candDs.cpa(), ptCandBs);
        registry.fill(HIST("hDecLengthDsRecBg"), candDs.decayLength(), ptCandBs);
        registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), ptCandBs);

        if (checkDecayTypeMc) {
          if (flagMcMatchRecBs == DecayChannelMain::B0ToDsPi) { // B0(bar) → Ds± π∓ → (K- K+ π±) π∓
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bs::DecayTypeMc::B0ToDsPiToPhiPiPiToKKPiPi, invMassCandBs, ptCandBs);
          } else if (flagMcMatchRecBs == hf_cand_bs::DecayTypeMc::PartlyRecoDecay) { // FIXME, Partly reconstructed decay channel
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bs::DecayTypeMc::PartlyRecoDecay, invMassCandBs, ptCandBs);
          } else {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bs::DecayTypeMc::NDecayTypeMc, invMassCandBs, ptCandBs);
          }
        }
      }
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == DecayChannelMain::BsToDsPi) {

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), o2::constants::physics::MassBS);
        if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
          continue;
        }

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }

        registry.fill(HIST("hPtProng0Gen"), ptProngs[0], ptParticle);
        registry.fill(HIST("hPtProng1Gen"), ptProngs[1], ptParticle);
        registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], ptParticle);
        registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], ptParticle);
        registry.fill(HIST("hYProng0Gen"), yProngs[0], ptParticle);
        registry.fill(HIST("hYProng1Gen"), yProngs[1], ptParticle);
        registry.fill(HIST("hPtGen"), ptParticle);
        registry.fill(HIST("hYGen"), yParticle, ptParticle);
        registry.fill(HIST("hEtaGen"), particle.eta(), ptParticle);

        // generated Bs with |y|<0.5
        if (std::abs(yParticle) < 0.5) {
          registry.fill(HIST("hPtGenWithRapidityBelowHalf"), ptParticle);
        }

        // reject Bs daughters that are not in geometrical acceptance
        if (!isProngInAcceptance(etaProngs[0], ptProngs[0]) || !isProngInAcceptance(etaProngs[1], ptProngs[1])) {
          continue;
        }
        registry.fill(HIST("hPtGenWithProngsInAcceptance"), ptParticle);
        registry.fill(HIST("hEtaGenWithProngsInAcceptance"), particle.eta(), ptParticle);
        registry.fill(HIST("hYGenWithProngsInAcceptance"), yParticle, ptParticle);
      }
    } // gen
  } // process
  PROCESS_SWITCH(HfTaskBs, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBs>(cfgc)};
}
