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

/// \file taskLb.cxx
/// \brief Λb0 analysis task
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
using namespace o2::aod::hf_cand_lb;
using namespace o2::analysis::hf_cuts_lb_tolcpi;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Λb0 analysis task
struct HfTaskLb {
  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "Lb candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hPtProng1", "Lb candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hPtCand", "Lb candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hCentrality", "centrality;centrality percentile;entries", {HistType::kTH1F, {{100, 0., 100.}}}}}};

  Configurable<int> selectionFlagLb{"selectionFlagLb", 1, "Selection Flag for Lb"};
  Configurable<double> cutYCandMax{"cutYCandMax", 1.44, "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_lb_tolcpi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    registry.add("hMass", "#Lambda_{b}^{0} candidates;inv. mass #Lambda_{c}^{#plus}#pi^{#minus} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c}); centrality", {HistType::kTH3F, {{500, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}, {100, 0., 100.}}});
    registry.add("hDecLength", "#Lambda_{b}^{0} candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXY", "#Lambda_{b}^{0} candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "#Lambda_{b}^{0} candidates;prong 0 (#Lambda_{c}^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "#Lambda_{b}^{0} candidates;prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidity", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hIPProd", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate impact parameter product;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hInvMassLc", "#Lambda_{b}^{0} candidates;prong0, #Lambda_{c}^{+} inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0, 5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_lb::isSelLbToLcPi >= selectionFlagLb);

  void process(soa::Join<aod::Collisions, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<soa::Join<aod::HfCandLb, aod::HFSelLbToLcPiCandidate>> const& candidates, soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate> const&, aod::BigTracks const&)
  {
    float centrality = collision.centRun2V0M();
    registry.fill(HIST("hCentrality"), centrality);

    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << hf_cand_lb::DecayType::LbToLcPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YLb(candidate)) > cutYCandMax) {
        continue;
      }

      auto candLc = candidate.index0_as<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>>();
      auto candPi = candidate.index1_as<aod::BigTracks>();

      registry.fill(HIST("hMass"), InvMassLbToLcPi(candidate), candidate.pt(), centrality);
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
      registry.fill(HIST("hRapidity"), YLb(candidate), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
      if (candPi.sign() < 0) {
        registry.fill(HIST("hInvMassLc"), InvMassLcpKpi(candLc), candidate.pt());
      }
    } // candidate loop
  }   // process
};    // struct

/// Lb MC analysis and fill histograms
struct HfTaskLbMc {
  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "Lb candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtRecBg", "Lb candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtGenSig", "Lb candidates (matched);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 10.}}}},
     {"hPtGen", "MC particles (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}}}};

  Configurable<int> selectionFlagLb{"selectionFlagLb", 1, "Selection Flag for Lb"};
  Configurable<double> cutYCandMax{"cutYCandMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_lb_tolcpi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    registry.add("hEtaGen", "MC particles (matched);#Lambda_{b}^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGen", "MC particles (matched);#Lambda_{b}^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0Gen", "MC particles (matched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYProng0Gen", "MC particles (matched);prong 0 (#Lambda_{c}^{+}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaProng0Gen", "MC particles (matched);prong 0 (#Lambda_{b}^{0}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPALcRecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPALcRecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidityRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidityRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate #it{#y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtProng0RecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "#Lambda_{b}^{0} candidates (matched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecSig", "#Lambda_{b}^{0} candidates (matched);inv. mass #Lambda_{c}^{+}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.00}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecBg", "#Lambda_{b}^{0} candidates (unmatched);inv. mass #Lambda_{c}^{+}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.0}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "#Lambda_{b}^{0} candidates (matched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthLcRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthNormRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthNormRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParProdLbRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParProdLbRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hChi2PCARecSig", "#Lambda_{b}^{0} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCARecBg", "#Lambda_{b}^{0} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hThetaStarRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} #cos(#theta^{*});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hThetaStarRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} #cos(#theta^{*});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_lb::isSelLbToLcPi >= selectionFlagLb);

  void process(soa::Filtered<soa::Join<aod::HfCandLb, aod::HFSelLbToLcPiCandidate, aod::HfCandLbMCRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandLbMCGen> const& particlesMC, aod::BigTracksMC const& tracks, aod::HfCandProng3 const&)
  {
    //MC rec
    for (auto& candidate : candidates) {
      //Printf("(Panos) MC candidate: pT: %lf",candidate.pt());
      if (!(candidate.hfflag() & 1 << hf_cand_lb::DecayType::LbToLcPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YLb(candidate)) > cutYCandMax) {
        continue;
      }
      auto candLc = candidate.index0_as<aod::HfCandProng3>();
      if (std::abs(candidate.flagMCMatchRec()) == 1 << hf_cand_lb::DecayType::LbToLcPi) {

        auto indexMother = RecoDecay::getMother(candidate.index1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandLbMCGen>>(), pdg::Code::kLambdaB0, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hRapidityRecSig"), YLb(candidate), candidate.pt());
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hMassRecSig"), InvMassLbToLcPi(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hImpParProdLbRecSig"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), candidate.pt());
        registry.fill(HIST("hCPALcRecSig"), candLc.cpa(), candidate.pt());
        registry.fill(HIST("hDecLengthLcRecSig"), candLc.decayLength(), candidate.pt());
        registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), candidate.pt());
        //registry.fill(HIST("hThetaStarRecSig"), candidate.cosThetaStar(), candidate.pt());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hRapidityRecBg"), YLb(candidate), candidate.pt());
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hMassRecBg"), InvMassLbToLcPi(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hImpParProdLbRecBg"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), candidate.pt());
        registry.fill(HIST("hCPALcRecBg"), candLc.cpa(), candidate.pt());
        registry.fill(HIST("hDecLengthLcRecBg"), candLc.decayLength(), candidate.pt());
        registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), candidate.pt());
        //registry.fill(HIST("hThetaStarRecBg"), candidate.cosThetaStar(), candidate.pt());
      }
    } // rec

    // MC gen. level
    //Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << hf_cand_lb::DecayType::LbToLcPi) {

        auto yParticle = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(pdg::Code::kLambdaB0));
        if (cutYCandMax >= 0. && std::abs(yParticle) > cutYCandMax) {
          continue;
        }

        float ptProngs[2], yProngs[2], etaProngs[2];
        int counter = 0;
        for (auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(array{daught.px(), daught.py(), daught.pz()}, RecoDecay::getMassPDG(daught.pdgCode()));
          counter++;
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
    } //gen
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  const bool doMC = cfgc.options().get<bool>("doMC");
  workflow.push_back(adaptAnalysisTask<HfTaskLb>(cfgc));
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfTaskLbMc>(cfgc));
  }
  return workflow;
}
