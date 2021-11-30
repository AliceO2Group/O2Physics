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
/// \brief Bs analysis task
///
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/Core/HFSelectorCuts.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/DataModel/Centrality.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_cand_bs;
using namespace o2::analysis::hf_cuts_bs_todspi;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Λb0 analysis task
struct HfTaskBs {
  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "Bs candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hPtProng1", "Bs candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hPtCand", "Bs candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hCentrality", "centrality;centrality percentile;entries", {HistType::kTH1F, {{100, 0., 100.}}}}}};

  Configurable<int> selectionFlagBs{"selectionFlagBs", 1, "Selection Flag for Bs"};
  Configurable<double> cutYCandMax{"cutYCandMax", 1.44, "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_bs_todspi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hMass", "B_{s} candidates;inv. mass D_{s}^{#plus}#pi^{#minus} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c}); centrality", {HistType::kTH3F, {{500, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, 0., 100.}}});
    registry.add("hDecLength", "B_{s} candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXY", "B_{s} candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "B_{s} candidates;prong 0 (D_{s}^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "B_{s} candidates;prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "B_{s} candidates;B_{s} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "B_{s} candidates;B_{s} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidity", "B_{s} candidates;B_{s} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "B_{s} candidates;B_{s} candidate impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "B_{s} candidates;B_{s} candidate decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "B_{s} candidates;B_{s} candidate decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hIPProd", "B_{s} candidates;B_{s} candidate impact parameter product;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hInvMassDs", "B_{s} candidates;prong0, D_{s}^{+} inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0, 5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_bs::isSelBsToDsPi >= selectionFlagBs);

  void process(soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, soa::Filtered<soa::Join<aod::HfCandBs, aod::HFSelBsToDsPiCandidate>> const& candidates, soa::Join<aod::HfCandProng3, aod::HFSelDsCandidate>, aod::BigTracks)
  {
    float centrality = collision.centV0M();
    registry.fill(HIST("hCentrality"), centrality);

    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << hf_cand_bs::DecayType::BsToDsPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YBs(candidate)) > cutYCandMax) {
        continue;
      }

      auto candDs = candidate.index0_as<soa::Join<aod::HfCandProng3, aod::HFSelDsCandidate>>();
      auto candPi = candidate.index1_as<aod::BigTracks>();

      registry.fill(HIST("hMass"), InvMassBsToDsPi(candidate), candidate.pt(), centrality);
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hIPProd"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hRapidity"), YBs(candidate), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
      if (candPi.sign() < 0) {
        registry.fill(HIST("hInvMassDs"), InvMassDsKKpi(candDs), candidate.pt());
      }
    } // candidate loop
  }   // process
};    // struct

/// Bs MC analysis and fill histograms
struct HfTaskBsMc {
  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "Bs candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtRecBg", "Bs candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtGenSig", "Bs candidates (matched);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 10.}}}},
     {"hPtGen", "MC particles (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}}}};

  Configurable<int> selectionFlagBs{"selectionFlagBs", 1, "Selection Flag for Bs"};
  Configurable<double> cutYCandMax{"cutYCandMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_bs_todspi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hEtaGen", "MC particles (matched);B_{s} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGen", "MC particles (matched);B_{s} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0Gen", "MC particles (matched);prong 0 (D_{s}^{+}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYProng0Gen", "MC particles (matched);prong 0 (D_{s}^{+}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaProng0Gen", "MC particles (matched);prong 0 (B_{s}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecSig", "B_{s} candidates (matched);B_{s} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "B_{s} candidates (unmatched);B_{s} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyRecSig", "B_{s} candidates (matched);B_{s} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyRecBg", "B_{s} candidates (unmatched);B_{s} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPADsRecSig", "B_{s} candidates (matched);prong 0 (D_{s}^{+}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPADsRecBg", "B_{s} candidates (unmatched);prong 0 (D_{s}^{+}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "B_{s} candidates (matched);B_{s} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "B_{s} candidates (unmatched);B_{s} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidityRecSig", "B_{s} candidates (matched);B_{s} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidityRecBg", "B_{s} candidates (unmatched);B_{s} candidate #it{#y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtProng0RecSig", "B_{s} candidates (matched);prong 0 (D_{s}^{+}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "B_{s} candidates (matched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "B_{s} candidates (unmatched);prong 0 (D_{s}^{+}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "B_{s} candidates (unmatched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecSig", "B_{s} candidates (matched);inv. mass D_{s}^{+}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.00}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecBg", "B_{s} candidates (unmatched);inv. mass D_{s}^{+}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.0}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "B_{s} candidates (matched);prong 0 (D_{s}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "B_{s} candidates (matched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "B_{s} candidates (unmatched);prong 0 (D_{s}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "B_{s} candidates (unmatched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecSig", "B_{s} candidates (matched);B_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYRecSig", "B_{s} candidates (matched);B_{s} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecBg", "B_{s} candidates (unmatched);B_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYRecBg", "B_{s} candidates (unmatched);B_{s} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthDsRecSig", "B_{s} candidates (matched);B_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthDsRecBg", "B_{s} candidates (unmatched);B_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthNormRecSig", "B_{s} candidates (matched);B_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthNormRecBg", "B_{s} candidates (unmatched);B_{s} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParProdBsRecSig", "B_{s} candidates (matched);B_{s} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParProdBsRecBg", "B_{s} candidates (unmatched);B_{s} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hChi2PCARecSig", "B_{s} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCARecBg", "B_{s} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hThetaStarRecSig", "B_{s} candidates (matched);B_{s} #cos(#theta^{*});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hThetaStarRecBg", "B_{s} candidates (unmatched);B_{s} #cos(#theta^{*});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_bs::isSelBsToDsPi >= selectionFlagBs);

  void process(soa::Filtered<soa::Join<aod::HfCandBs, aod::HFSelBsToDsPiCandidate, aod::HfCandBsMCRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandBsMCGen> const& particlesMC, aod::BigTracksMC const& tracks, aod::HfCandProng3 const&)
  {
    // MC rec
    for (auto& candidate : candidates) {
      // Printf("(Panos) MC candidate: pT: %lf",candidate.pt());
      if (!(candidate.hfflag() & 1 << hf_cand_bs::DecayType::BsToDsPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YBs(candidate)) > cutYCandMax) {
        continue;
      }
      auto candDs = candidate.index0_as<aod::HfCandProng3>();
      if (std::abs(candidate.flagMCMatchRec()) == 1 << hf_cand_bs::DecayType::BsToDsPi) {

        auto indexMother = RecoDecay::getMother(particlesMC, candidate.index1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandBsMCGen>>(), pdg::Code::kBs, true);
        auto particleMother = particlesMC.iteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hRapidityRecSig"), YBs(candidate), candidate.pt());
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hMassRecSig"), InvMassBsToDsPi(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hImpParProdBsRecSig"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), candidate.pt());
        registry.fill(HIST("hCPADsRecSig"), candDs.cpa(), candidate.pt());
        registry.fill(HIST("hDecLengthDsRecSig"), candDs.decayLength(), candidate.pt());
        registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), candidate.pt());
        // registry.fill(HIST("hThetaStarRecSig"), candidate.cosThetaStar(), candidate.pt());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hRapidityRecBg"), YBs(candidate), candidate.pt());
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hMassRecBg"), InvMassBsToDsPi(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hImpParProdBsRecBg"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), candidate.pt());
        registry.fill(HIST("hCPADsRecBg"), candDs.cpa(), candidate.pt());
        registry.fill(HIST("hDecLengthDsRecBg"), candDs.decayLength(), candidate.pt());
        registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), candidate.pt());
        // registry.fill(HIST("hThetaStarRecBg"), candidate.cosThetaStar(), candidate.pt());
      }
    } // rec

    // MC gen. level
    // Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << hf_cand_bs::DecayType::BsToDsPi) {

        auto yParticle = RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(pdg::Code::kBs));
        if (cutYCandMax >= 0. && std::abs(yParticle) > cutYCandMax) {
          continue;
        }

        float ptProngs[2], yProngs[2], etaProngs[2];
        for (int iD = particle.daughter0Id(), counter = 0; iD <= particle.daughter1Id(); ++iD, counter++) {
          auto daught = particlesMC.iteratorAt(iD);
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::Y(array{daught.px(), daught.py(), daught.pz()}, RecoDecay::getMassPDG(daught.pdgCode()));
        }

        registry.fill(HIST("hPtProng0Gen"), ptProngs[0], particle.pt());
        registry.fill(HIST("hPtProng1Gen"), ptProngs[1], particle.pt());
        registry.fill(HIST("hYProng0Gen"), yProngs[0], particle.pt());
        registry.fill(HIST("hYProng1Gen"), yProngs[1], particle.pt());
        registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], particle.pt());
        registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], particle.pt());

        //  if (cutYCandMax >= 0. && (std::abs(yProngs[0]) > cutYCandMax || std::abs(yProngs[1]) > cutYCandMax))
        //    continue;

        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hYGen"), yParticle, particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta(), particle.pt());
      }
    } // gen
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  const bool doMC = cfgc.options().get<bool>("doMC");
  workflow.push_back(adaptAnalysisTask<HfTaskBs>(cfgc));
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfTaskBsMc>(cfgc));
  }
  return workflow;
}
