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

/// \file taskB0.cxx
/// \brief B0 → D- π+ → (π- K+ π-) π+ analysis task
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

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

/// B0 analysis task
struct HfTaskB0 {
  Configurable<int> selectionFlagB0{"selectionFlagB0", 1, "Selection Flag for B0"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_b0_to_d_pi::vecBinsPt}, "pT bin limits"};
  // MC checks
  Configurable<bool> checkDecayTypeMc{"checkDecayTypeMc", false, "Flag to enable DecayType histogram"};

  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg{};

  using TracksWithSel = soa::Join<aod::Tracks, aod::TrackSelection>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_b0::isSelB0ToDPi >= selectionFlagB0);

  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "B0 candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hPtProng1", "B0 candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hPtCand", "B0 candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}}}};

  void init(InitContext const&)
  {
    static const AxisSpec axisMassB0 = {300, 4.5, 6.0, "inv. mass (GeV/#it{c}^{2})"};
    static const AxisSpec axisPt = {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"};

    registry.add("hMass", "B^{0} candidates;inv. mass D^{#minus}#pi^{#plus} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisMassB0, axisPt}});
    registry.add("hDecLength", "B^{0} candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, axisPt}});
    registry.add("hDecLengthXY", "B^{0} candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, axisPt}});
    registry.add("hd0Prong0", "B^{0} candidates;prong 0 (D^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong1", "B^{0} candidates;prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, axisPt}});
    registry.add("hCPA", "B^{0} candidates;B^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
    registry.add("hEta", "B^{0} candidates;B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hRapidity", "B^{0} candidates;B^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hImpParErr", "B^{0} candidates;B^{0} candidate impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, axisPt}});
    registry.add("hDecLenErr", "B^{0} candidates;B^{0} candidate decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, axisPt}});
    registry.add("hDecLenXYErr", "B^{0} candidates;B^{0} candidate decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, axisPt}});
    registry.add("hIPProd", "B^{0} candidates;B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, axisPt}});
    registry.add("hInvMassD", "B^{0} candidates;prong0, D^{#minus} inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0, 5}, axisPt}});

    registry.add("hEtaGen", "MC particles (generated);B^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hYGen", "MC particles (generated);B^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hPtProng0Gen", "MC particles (generated);prong 0 (D^{#minus}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hPtProng1Gen", "MC particles (generated);prong 1 (#pi^{-}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hYProng0Gen", "MC particles (generated);prong 0 (D^{#minus}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
    registry.add("hYProng1Gen", "MC particles (generated);prong 1 (#pi^{-}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
    registry.add("hEtaProng0Gen", "MC particles (generated);prong 0 (B^{0}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
    registry.add("hEtaProng1Gen", "MC particles (generated);prong 1 (#pi^{-}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
    registry.add("hCPARecSig", "B^{0} candidates (matched);B^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPARecBg", "B^{0} candidates (unmatched);B^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPAxyRecSig", "B^{0} candidates (matched);B^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPAxyRecBg", "B^{0} candidates (unmatched);B^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPADRecSig", "B^{0} candidates (matched);prong 0 (D^{#minus}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hCPADRecBg", "B^{0} candidates (unmatched);prong 0 (D^{#minus}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, axisPt}});
    registry.add("hEtaRecSig", "B^{0} candidates (matched);B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hEtaRecBg", "B^{0} candidates (unmatched);B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hRapidityRecSig", "B^{0} candidates (matched);B^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hRapidityRecBg", "B^{0} candidates (unmatched);B^{0} candidate #it{#y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});

    registry.add("hPtProng0RecSig", "B^{0} candidates (matched);prong 0 (D^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hPtProng1RecSig", "B^{0} candidates (matched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hPtProng0RecBg", "B^{0} candidates (unmatched);prong 0 (D^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hPtProng1RecBg", "B^{0} candidates (unmatched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, axisPt}});
    registry.add("hMassRecSig", "B^{0} candidates (matched);inv. mass D^{#minus}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.00}, axisPt}});
    registry.add("hMassRecBg", "B^{0} candidates (unmatched);inv. mass D^{#minus}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.0}, axisPt}});
    registry.add("hd0Prong0RecSig", "B^{0} candidates (matched);prong 0 (D^{#minus}}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong1RecSig", "B^{0} candidates (matched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong0RecBg", "B^{0} candidates (unmatched);prong 0 (D^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong1RecBg", "B^{0} candidates (unmatched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
    registry.add("hDecLengthRecSig", "B^{0} candidates (matched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthXYRecSig", "B^{0} candidates (matched);B^{0} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthRecBg", "B^{0} candidates (unmatched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthXYRecBg", "B^{0} candidates (unmatched);B^{0} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthDRecSig", "B^{0} candidates (matched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthDRecBg", "B^{0} candidates (unmatched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthNormRecSig", "B^{0} candidates (matched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hDecLengthNormRecBg", "B^{0} candidates (unmatched);B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
    registry.add("hImpParProdB0RecSig", "B^{0} candidates (matched);B^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.01, 0.01}, axisPt}});
    registry.add("hImpParProdB0RecBg", "B^{0} candidates (unmatched);B^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.01, 0.01}, axisPt}});

    registry.add("hChi2PCARecSig", "B^{0} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, axisPt}});
    registry.add("hChi2PCARecBg", "B^{0} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, axisPt}});

    registry.add("hPtRecSig", "B0 candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtRecBg", "B0 candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenSig", "B0 candidates (gen+rec);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGen", "MC particles (generated);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenWithRapidityBelowHalf", "MC particles (generated - |#it{y}^{gen}|<0.5);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
    registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});

    if (checkDecayTypeMc) {
      constexpr uint8_t NBinsDecayTypeMc = hf_cand_b0::DecayTypeMc::NDecayTypeMc; // FIXME
      TString labels[NBinsDecayTypeMc];
      labels[hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi] = "B^{0} #rightarrow (D^{#minus} #rightarrow #pi^{#minus} K^{#plus} #pi^{#minus}) #pi^{#plus}";
      labels[hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi] = "B^{0} #rightarrow (D^{#minus}_{s} #rightarrow K^{#minus} K^{#plus} #pi^{#minus}) #pi^{#plus}";
      labels[hf_cand_b0::DecayTypeMc::PartlyRecoDecay] = "Partly reconstructed decay channel";
      labels[hf_cand_b0::DecayTypeMc::OtherDecay] = "Other decays";
      static const AxisSpec axisDecayType = {NBinsDecayTypeMc, 0.5, NBinsDecayTypeMc + 0.5, ""};
      registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassB0, axisPt}});
      for (uint8_t iBin = 0; iBin < NBinsDecayTypeMc; ++iBin) {
        registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
      }
    }
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

  void process(soa::Filtered<soa::Join<aod::HfCandB0, aod::HfSelB0ToDPi>> const& candidates,
               soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const&,
               TracksWithSel const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }

      auto candD = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>();
      // auto candPi = candidate.prong1_as<TracksWithSel>();

      auto ptCandB0 = candidate.pt();

      registry.fill(HIST("hMass"), HfHelper::invMassB0ToDPi(candidate), ptCandB0);
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
      registry.fill(HIST("hRapidity"), HfHelper::yB0(candidate), ptCandB0);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), ptCandB0);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), ptCandB0);
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), ptCandB0);
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), ptCandB0);
      registry.fill(HIST("hInvMassD"), HfHelper::invMassDplusToPiKPi(candD), ptCandB0);
    } // candidate loop
  } // process

  /// B0 MC analysis and fill histograms
  void processMc(soa::Filtered<soa::Join<aod::HfCandB0, aod::HfSelB0ToDPi, aod::HfCandB0McRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandB0McGen> const& mcParticles,
                 aod::TracksWMc const&,
                 soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec> const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }

      auto ptCandB0 = candidate.pt();
      auto candD = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec>>();
      auto invMassCandB0 = HfHelper::invMassB0ToDPi(candidate);
      auto flagMcMatchRecB0 = std::abs(candidate.flagMcMatchRec());

      if (flagMcMatchRecB0 == DecayChannelMain::B0ToDminusPi) {
        auto indexMother = RecoDecay::getMother(mcParticles, candidate.prong1_as<aod::TracksWMc>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandB0McGen>>(), o2::constants::physics::Pdg::kB0, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);

        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), ptCandB0);
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), ptCandB0);
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpaXY(), ptCandB0);
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), ptCandB0);
        registry.fill(HIST("hRapidityRecSig"), HfHelper::yB0(candidate), ptCandB0);
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), ptCandB0);
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), ptCandB0);
        registry.fill(HIST("hMassRecSig"), HfHelper::invMassB0ToDPi(candidate), ptCandB0);
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), ptCandB0);
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), ptCandB0);
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), ptCandB0);
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), ptCandB0);
        registry.fill(HIST("hImpParProdB0RecSig"), candidate.impactParameterProduct(), ptCandB0);
        registry.fill(HIST("hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), ptCandB0);
        registry.fill(HIST("hCPADRecSig"), candD.cpa(), ptCandB0);
        registry.fill(HIST("hDecLengthDRecSig"), candD.decayLength(), ptCandB0);
        registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), ptCandB0);

        if (checkDecayTypeMc) {
          registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi, invMassCandB0, ptCandB0);
        }
      } else {
        registry.fill(HIST("hPtRecBg"), ptCandB0);
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), ptCandB0);
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpaXY(), ptCandB0);
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), ptCandB0);
        registry.fill(HIST("hRapidityRecBg"), HfHelper::yB0(candidate), ptCandB0);
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), ptCandB0);
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), ptCandB0);
        registry.fill(HIST("hMassRecBg"), HfHelper::invMassB0ToDPi(candidate), ptCandB0);
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), ptCandB0);
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), ptCandB0);
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), ptCandB0);
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), ptCandB0);
        registry.fill(HIST("hImpParProdB0RecBg"), candidate.impactParameterProduct(), ptCandB0);
        registry.fill(HIST("hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), ptCandB0);
        registry.fill(HIST("hCPADRecBg"), candD.cpa(), ptCandB0);
        registry.fill(HIST("hDecLengthDRecBg"), candD.decayLength(), ptCandB0);
        registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), ptCandB0);

        if (checkDecayTypeMc) {
          if (flagMcMatchRecB0 == DecayChannelMain::B0ToDsPi) { // B0 → Ds- π+ → (K- K+ π-) π+
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi, invMassCandB0, ptCandB0);
          } else if (flagMcMatchRecB0 == hf_cand_b0::DecayTypeMc::PartlyRecoDecay) { // FIXME, Partly reconstructed decay channel
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::PartlyRecoDecay, invMassCandB0, ptCandB0);
          } else {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::OtherDecay, invMassCandB0, ptCandB0);
          }
        }
      }
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_beauty::DecayChannelMain::B0ToDminusPi) {

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), o2::constants::physics::MassB0);
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
        registry.fill(HIST("hYProng0Gen"), yProngs[0], ptParticle);
        registry.fill(HIST("hYProng1Gen"), yProngs[1], ptParticle);
        registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], ptParticle);
        registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], ptParticle);

        registry.fill(HIST("hPtGen"), ptParticle);
        registry.fill(HIST("hYGen"), yParticle, ptParticle);
        registry.fill(HIST("hEtaGen"), particle.eta(), ptParticle);

        // generated B0 with |y|<0.5
        if (std::abs(yParticle) < 0.5) {
          registry.fill(HIST("hPtGenWithRapidityBelowHalf"), ptParticle);
        }

        // reject B0 daughters that are not in geometrical acceptance
        if (!isProngInAcceptance(etaProngs[0], ptProngs[0]) || !isProngInAcceptance(etaProngs[1], ptProngs[1])) {
          continue;
        }
        registry.fill(HIST("hPtGenWithProngsInAcceptance"), ptParticle);
        registry.fill(HIST("hYGenWithProngsInAcceptance"), yParticle, ptParticle);
        registry.fill(HIST("hEtaGenWithProngsInAcceptance"), particle.eta(), ptParticle);
      }
    } // gen
  } // process
  PROCESS_SWITCH(HfTaskB0, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskB0>(cfgc)};
}
