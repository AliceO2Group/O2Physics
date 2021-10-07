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

/// \file taskChic.cxx
/// \brief Chi_c1(1P) analysis task. Adapted from X 
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Alessandro De Falco <alessandro.de.falco@ca.infn.it>, Cagliari University

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/Core/HFSelectorCuts.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::analysis::hf_cuts_chic_tojpsigamma;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_cand_chic;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_prong2;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// chi_c1(1P) analysis task
struct TaskChic {
  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}}}};

  Configurable<int> d_selectionFlagChic{"d_selectionFlagChic", 1, "Selection Flag for Chic"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_chic_tojpsigamma::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    registry.add("hMass", "2-prong candidates;inv. mass (J/#psi #gamma) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hdeclength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hdeclengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_chic::isSelChicToJpsiGamma >= d_selectionFlagChic);

  void process(soa::Filtered<soa::Join<aod::HfCandChic, aod::HFSelChicToJpsiGammaCandidate>> const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << ChicToJpsiGamma)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YChic(candidate)) > cutYCandMax) {
        continue;
      }

      registry.fill(HIST("hMass"), InvMassChicToJpsiGamma(candidate), candidate.pt());
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hdeclength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hdeclengthxy"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    } // candidate loop
  }   // process
};    // struct

struct TaskChicMC {
  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "2-prong candidates (rec. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtRecBg", "2-prong candidates (rec. unmatched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtGen", "2-prong candidates (gen. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtGenSig", "2-prong candidates (rec. matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}}}};

  Configurable<int> d_selectionFlagChic{"d_selectionFlagChic", 1, "Selection Flag for Chic"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_chic_tojpsigamma::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    registry.add("hCPARecSig", "2-prong candidates (rec. matched);cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "2-prong candidates (rec. unmatched);cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "2-prong candidates (rec. matched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "2-prong candidates (rec. unmatched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaGen", "2-prong candidates (gen. matched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtProng0RecSig", "2-prong candidates (rec. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "2-prong candidates (rec. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "2-prong candidates (rec. unmatched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "2-prong candidates (rec. unmatched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng0", "2-prong candidates (gen. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng1", "2-prong candidates (gen. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hMassRecSig", "2-prong candidates (rec. matched);inv. mass (J/#psi #gamma) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecBg", "2-prong candidates (rec. unmatched);inv. mass (J/#psi #gamma) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "2-prong candidates (rec. matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "2-prong candidates (rec. matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "2-prong candidates (rec. unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "2-prong candidates (rec. unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeclengthRecSig", "2-prong candidates (rec. matched);decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeclengthRecBg", "2-prong candidates (rec. unmatched);decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hChi2PCASig", "2-prong candidates (rec. matched);chi2 PCA (cm);entries", {HistType::kTH2F, {{1000, 0., 0.0001}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCABg", "2-prong candidates (rec. unmatched);chi2 PCA (cm);entries", {HistType::kTH2F, {{1000, 0., 0.0001}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtSig", "2-prong candidates (rec. matched);proper lifetime chi_c * #it{c} (cm);entries", {HistType::kTH2F, {{400, 0., 0.001}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtBg", "2-prong candidates (rec. unmatched);proper lifetime chi_c * #it{c} (cm);entries", {HistType::kTH2F, {{400, 0., 0.001}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYSig", "2-prong candidates (rec. matched);candidate rapidity;entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYBg", "2-prong candidates (rec. unmatched);candidate rapidity;entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_chic::isSelChicToJpsiGamma >= d_selectionFlagChic);

  void process(soa::Filtered<soa::Join<aod::HfCandChic, aod::HFSelChicToJpsiGammaCandidate, aod::HfCandChicMCRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandChicMCGen> const& particlesMC, aod::BigTracksMC const& tracks)
  {
    // MC rec.
    //Printf("MC Candidates: %d", candidates.size());
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << ChicToJpsiGamma)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YChic(candidate)) > cutYCandMax) {
        continue;
      }
      if (candidate.flagMCMatchRec() == 1 << ChicToJpsiGamma) {
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.index1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandChicMCGen>>(), 20443, true);
        auto particleMother = particlesMC.iteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());

        registry.fill(HIST("hDeclengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hMassRecSig"), InvMassChicToJpsiGamma(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hChi2PCASig"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCtSig"), CtChic(candidate), candidate.pt());
        registry.fill(HIST("hYSig"), YChic(candidate), candidate.pt());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());

        registry.fill(HIST("hDeclengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hMassRecBg"), InvMassChicToJpsiGamma(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hChi2PCABg"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCtBg"), CtChic(candidate), candidate.pt());
        registry.fill(HIST("hYBg"), YChic(candidate), candidate.pt());
      }
    } // rec
    // MC gen.
    //Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (particle.flagMCMatchGen() == 1 << ChicToJpsiGamma) {
	auto mchic = RecoDecay::getMassPDG(20443); // chi_c1(1p)
        if (cutYCandMax >= 0. && std::abs(RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, mchic)) > cutYCandMax) {
          // Printf("MC Gen.: Y rejection: %g", RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, 3.87168));
          continue;
        }
        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta(), particle.pt());

        // properties of gen matched chic, to get a first look at some cuts
        float ptProngs[3];
        int counter = 0;
        for (int iD = particle.daughter0Id(); iD <= particle.daughter1Id(); ++iD) {  //adf check this
          ptProngs[counter] = particlesMC.iteratorAt(iD).pt();
          counter++;
        }
        registry.fill(HIST("hPtGenProng0"), ptProngs[0], particle.pt());
        registry.fill(HIST("hPtGenProng1"), ptProngs[1], particle.pt());
      }
    } // gen
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TaskChic>(cfgc, TaskName{"hf-task-chic"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<TaskChicMC>(cfgc, TaskName{"hf-task-chic-mc"}));
  }
  return workflow;
}
