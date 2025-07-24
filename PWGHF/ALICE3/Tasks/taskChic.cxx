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

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// chi_c1(1P) analysis task
struct HfTaskChic {
  Configurable<int> selectionFlagChic{"selectionFlagChic", 1, "Selection Flag for Chic"};
  Configurable<double> yCandMax{"yCandMax", 1., "max. cand. rapidity"};
  Configurable<bool> modeChicToJpsiToMuMuGamma{"modeChicToJpsiToMuMuGamma", true, "Perform Jpsi to mu+mu- analysis"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_chic_to_jpsi_gamma::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_chic::isSelChicToJpsiToEEGamma >= selectionFlagChic || aod::hf_sel_candidate_chic::isSelChicToJpsiToMuMuGamma >= selectionFlagChic);

  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}}}};

  void init(InitContext&)
  {
    registry.add("hMass", "2-prong candidates;inv. mass (J/#psi #gamma) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeltaMass", "2-prong candidates;inv. mass (J/#psi #gamma) - inv. mass (J/#psi) + mass^{PDG} (J/#psi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXY", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    //    registry.add("hEGamma", "Photon energy", {HistType::kTH1F, {{200, 0., 10.}}});
  }

  void process(soa::Filtered<soa::Join<aod::HfCandChic, aod::HfSelChicToJpsiGamma>> const& candidates)
  {
    int decayMode = modeChicToJpsiToMuMuGamma ? hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma : hf_cand_chic::DecayType::ChicToJpsiToEEGamma;
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << decayMode)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(hfHelper.yChic(candidate)) > yCandMax) {
        continue;
      }

      registry.fill(HIST("hMass"), hfHelper.invMassChicToJpsiGamma(candidate), candidate.pt());
      registry.fill(HIST("hDeltaMass"), hfHelper.invMassChicToJpsiGamma(candidate) - candidate.jpsiToMuMuMass() + o2::constants::physics::MassJPsi, candidate.pt());
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
      //      registry.fill(HIST("hEGamma"), candidate.prong1().e());
    } // candidate loop
  }   // process
};    // struct

struct HfTaskChicMc {
  Configurable<int> selectionFlagChic{"selectionFlagChic", 1, "Selection Flag for Chic"};
  Configurable<double> yCandMax{"yCandMax", 1., "max. cand. rapidity"};
  Configurable<bool> modeChicToJpsiToMuMuGamma{"modeChicToJpsiToMuMuGamma", true, "Perform Jpsi to mu+mu- analysis"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_chic_to_jpsi_gamma::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_chic::isSelChicToJpsiToEEGamma >= selectionFlagChic || aod::hf_sel_candidate_chic::isSelChicToJpsiToMuMuGamma >= selectionFlagChic);

  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "2-prong candidates (rec. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtRecBg", "2-prong candidates (rec. unmatched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtGen", "2-prong candidates (gen. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtGenSig", "2-prong candidates (rec. matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}}}};

  void init(InitContext&)
  {
    registry.add("hCPARecSig", "2-prong candidates (rec. matched);cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "2-prong candidates (rec. unmatched);cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "2-prong candidates (rec. matched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "2-prong candidates (rec. unmatched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaGen", "2-prong candidates (gen. matched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtProng0RecSig", "2-prong candidates (rec. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "2-prong candidates (rec. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "2-prong candidates (rec. unmatched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "2-prong candidates (rec. unmatched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0Gen", "2-prong candidates (gen. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1Gen", "2-prong candidates (gen. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hDeltaMassRecSig", "2-prong candidates (rec. matched);inv. mass (J/#psi #gamma) - inv. mass (J/#psi) + mass^{PDG} (J/#psi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeltaMassRecBg", "2-prong candidates (rec. unmatched);inv. mass (J/#psi #gamma) - inv. mass (J/#psi) + mass^{PDG} (J/#psi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecSig", "2-prong candidates (rec. matched);inv. mass (J/#psi #gamma) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecBg", "2-prong candidates (rec. unmatched);inv. mass (J/#psi #gamma) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "2-prong candidates (rec. matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "2-prong candidates (rec. matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "2-prong candidates (rec. unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "2-prong candidates (rec. unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecSig", "2-prong candidates (rec. matched);decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecBg", "2-prong candidates (rec. unmatched);decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hChi2PCARecSig", "2-prong candidates (rec. matched);chi2 PCA (J/#psi)(cm);entries", {HistType::kTH2F, {{5000, 0., 5e-6}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCARecBg", "2-prong candidates (rec. unmatched);chi2 PCA (J/#psi)(cm);entries", {HistType::kTH2F, {{5000, 0., 5e-6}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtRecSig", "2-prong candidates (rec. matched);proper lifetime chi_c * #it{c} (cm);entries", {HistType::kTH2F, {{400, 0., 0.001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtRecBg", "2-prong candidates (rec. unmatched);proper lifetime chi_c * #it{c} (cm);entries", {HistType::kTH2F, {{400, 0., 0.001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYRecSig", "2-prong candidates (rec. matched);candidate rapidity;entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYRecBg", "2-prong candidates (rec. unmatched);candidate rapidity;entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(soa::Filtered<soa::Join<aod::HfCandChic, aod::HfSelChicToJpsiGamma, aod::HfCandChicMcRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandChicMcGen> const& mcParticles,
               aod::TracksWMc const&)
  {
    // MC rec.
    int decayMode = modeChicToJpsiToMuMuGamma ? hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma : hf_cand_chic::DecayType::ChicToJpsiToEEGamma;
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << decayMode)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(hfHelper.yChic(candidate)) > yCandMax) {
        continue;
      }
      if (candidate.flagMcMatchRec() == 1 << decayMode) {
        // FIXME the access to the MC particle gen not yet functional
        // int indexMother = RecoDecay::getMother(mcParticles, mcParticles.rawIteratorAt(candidate.prong1().mcParticle_as<aod::McParticles>().globalIndex()), 20443);
        // auto particleMother = mcParticles.rawIteratorAt(indexMother);
        // registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDeltaMassRecSig"), hfHelper.invMassChicToJpsiGamma(candidate) - candidate.jpsiToMuMuMass() + o2::constants::physics::MassJPsi), candidate.pt();
        registry.fill(HIST("hMassRecSig"), hfHelper.invMassChicToJpsiGamma(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCtRecSig"), hfHelper.ctChic(candidate), candidate.pt());
        registry.fill(HIST("hYRecSig"), hfHelper.yChic(candidate), candidate.pt());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDeltaMassRecBg"), hfHelper.invMassChicToJpsiGamma(candidate) - candidate.jpsiToMuMuMass() + o2::constants::physics::MassJPsi), candidate.pt();
        registry.fill(HIST("hMassRecBg"), hfHelper.invMassChicToJpsiGamma(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCtRecBg"), hfHelper.ctChic(candidate), candidate.pt());
        registry.fill(HIST("hYRecBg"), hfHelper.yChic(candidate), candidate.pt());
      }
    } // rec
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (particle.flagMcMatchGen() == 1 << decayMode) {
        auto mchic = o2::constants::physics::MassChiC1; // chi_c1(1p)
        if (yCandMax >= 0. && std::abs(RecoDecay::y(particle.pVector(), mchic)) > yCandMax) {
          continue;
        }

        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta(), particle.pt());

        // properties of gen matched chic, to get a first look at some cuts
        float ptProngs[3];
        int counter = 0;
        for (const auto& dau : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = dau.pt();
          counter++;
        }

        registry.fill(HIST("hPtProng0Gen"), ptProngs[0], particle.pt());
        registry.fill(HIST("hPtProng1Gen"), ptProngs[1], particle.pt());
      }
    } // gen
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfTaskChic>(cfgc)};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfTaskChicMc>(cfgc));
  }
  return workflow;
}
