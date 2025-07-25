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

/// \file taskJpsi.cxx
/// \brief Jpsi analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// jpsitoee analysis task
struct HfTaskJpsi {
  Configurable<int> selectionFlagJpsi{"selectionFlagJpsi", 0, "Selection Flag for Jpsi"};
  Configurable<bool> modeJpsiToMuMu{"modeJpsiToMuMu", false, "Perform Jpsi to mu+mu- analysis"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<bool> selectedTof{"selectedTof", false, "select TOF for Jpsi"};
  Configurable<bool> selectedRich{"selectedRich", false, "select RICH for Jpsi"};
  Configurable<bool> selectedTofRich{"selectedTofRich", false, "select TOF and RICH for Jpsi"};
  Configurable<bool> selectedMid{"selectedMid", false, "select MID for Jpsi to mu+mu-"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_jpsi_to_e_e::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_jpsi::isSelJpsiToEETopol >= selectionFlagJpsi || aod::hf_sel_candidate_jpsi::isSelJpsiToMuMuTopol >= selectionFlagJpsi);

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 20.}}}},
     {"hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 20.}}}},
     {"hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 20.}}}}}};

  void init(InitContext&)
  {
    if (modeJpsiToMuMu) {
      registry.add("hMass", "2-prong candidates;inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    } else {
      registry.add("hMass", "2-prong candidates;inv. mass (e^{#plus} e^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    }
    registry.add("hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{1000, 0.5, 1.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelJpsi>> const& candidates)
  {
    int decayMode = modeJpsiToMuMu ? aod::hf_cand_2prong::DecayType::JpsiToMuMu : aod::hf_cand_2prong::DecayType::JpsiToEE;

    for (const auto& candidate : candidates) {

      if (!(candidate.hfflag() & 1 << decayMode)) {
        continue;
      }
      if (selectionFlagJpsi > 0) {
        if (modeJpsiToMuMu) {
          if (candidate.isSelJpsiToMuMuTopol() <= 0) {
            continue;
          }
          if (selectedMid && candidate.isSelJpsiToMuMuMid() <= 0) {
            continue;
          }
        } else {
          if (candidate.isSelJpsiToEETopol() <= 0) {
            continue;
          }
          if (selectedTof && candidate.isSelJpsiToEETof() <= 0) {
            continue;
          }
          if (selectedRich && candidate.isSelJpsiToEERich() <= 0) {
            continue;
          }
          if (selectedTofRich && candidate.isSelJpsiToEETofRich() <= 0) {
            continue;
          }
        }
      }
      if (yCandMax >= 0. && std::abs(hfHelper.yJpsi(candidate)) > yCandMax) {
        continue;
      }

      if (modeJpsiToMuMu) {
        registry.fill(HIST("hMass"), hfHelper.invMassJpsiToMuMu(candidate), candidate.pt());
      } else {
        registry.fill(HIST("hMass"), hfHelper.invMassJpsiToEE(candidate), candidate.pt());
      }
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthxy"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    }
  }
};

/// Fills MC histograms.
struct HfTaskJpsiMc {
  Configurable<int> selectionFlagJpsi{"selectionFlagJpsi", 1, "Selection Flag for Jpsi"};
  Configurable<bool> modeJpsiToMuMu{"modeJpsiToMuMu", false, "Perform Jpsi to mu+mu- analysis"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<bool> selectedTof{"selectedTof", false, "select TOF for Jpsi"};
  Configurable<bool> selectedRich{"selectedRich", false, "select RICH for Jpsi"};
  Configurable<bool> selectedTofRich{"selectedTofRich", false, "select TOF and RICH for Jpsi"};
  Configurable<bool> selectedMid{"selectedMid", false, "select MID for Jpsi to mu+mu-"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_jpsi_to_e_e::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  using McParticlesHf = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_jpsi::isSelJpsiToEETopol >= selectionFlagJpsi || aod::hf_sel_candidate_jpsi::isSelJpsiToMuMuTopol >= selectionFlagJpsi);

  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "2-prong candidates (rec. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 20.}}}},
     {"hPtRecBg", "2-prong candidates (rec. unmatched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 20.}}}},
     {"hPtGen", "2-prong candidates (gen. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 20.}}}},
     {"hPtGenSig", "2-prong candidates (rec. matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 20.}}}},
     {"hCPARecSig", "2-prong candidates (rec. matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBg", "2-prong candidates (rec. unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "2-prong candidates (rec. matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaRecBg", "2-prong candidates (rec. unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaGen", "2-prong candidates (gen. matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}}}};

  void init(InitContext&)
  {
    if (modeJpsiToMuMu) {
      registry.add("hMassSig", "2-prong candidates (rec matched);inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hMassBg", "2-prong candidates (rec unmatched);inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    } else {
      registry.add("hMassSig", "2-prong candidates (rec matched);inv. mass (e^{#plus} e^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hMassBg", "2-prong candidates (rec unmatched);inv. mass (e^{#plus} e^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    }
    registry.add("hDecLengthSig", "2-prong candidates (rec matched);decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthBg", "2-prong candidates (rec unmatched);decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYSig", "2-prong candidates (rec matched);decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxyBg", "2-prong candidates (rec unmatched);decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0Sig", "2-prong candidates (rec matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0Bg", "2-prong candidates (rec unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1Sig", "2-prong candidates (rec matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1Bg", "2-prong candidates (rec unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0Sig", "2-prong candidates (rec matched);product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0Bg", "2-prong candidates (rec unmatched);product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCASig", "2-prong candidates (rec. matched);chi2 PCA (cm);entries", {HistType::kTH2F, {{1000, 0., 0.0001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCABg", "2-prong candidates (rec. unmatched);chi2 PCA (cm);entries", {HistType::kTH2F, {{1000, 0., 0.0001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtSig", "2-prong candidates (rec. matched);proper lifetime X(3872) * #it{c} (cm);entries", {HistType::kTH2F, {{400, 0., 0.001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtBg", "2-prong candidates (rec. unmatched);proper lifetime X(3872) * #it{c} (cm);entries", {HistType::kTH2F, {{400, 0., 0.001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGenSig", "2-prong candidates (rec. matched);candidate rapidity;entries", {HistType::kTH2F, {{10, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}^{gen}_{T} (GeV/#it{c})"}}});
    registry.add("hYSig", "2-prong candidates (rec. matched);candidate rapidity;entries", {HistType::kTH2F, {{10, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYBg", "2-prong candidates (rec. unmatched);candidate rapidity;entries", {HistType::kTH2F, {{10, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGen", "2-prong MC particles (gen. matched);candidate rapidity;entries", {HistType::kTH2F, {{10, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng0", "2-prong candidates (gen. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng1", "2-prong candidates (gen. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelJpsi, aod::HfCand2ProngMcRec>> const& candidates,
               McParticlesHf const& mcParticles,
               aod::TracksWMc const&)
  {
    // MC rec.
    int decayMode = modeJpsiToMuMu ? aod::hf_cand_2prong::DecayType::JpsiToMuMu : aod::hf_cand_2prong::DecayType::JpsiToEE;

    for (const auto& candidate : candidates) {

      if (!(candidate.hfflag() & 1 << decayMode)) {
        continue;
      }
      if (selectionFlagJpsi > 0) {
        if (modeJpsiToMuMu) {
          if (candidate.isSelJpsiToMuMuTopol() <= 0) {
            continue;
          }
          if (selectedMid && candidate.isSelJpsiToMuMuMid() <= 0) {
            continue;
          }
        } else {
          if (candidate.isSelJpsiToEETopol() <= 0) {
            continue;
          }
          if (selectedTof && candidate.isSelJpsiToEETof() <= 0) {
            continue;
          }
          if (selectedRich && candidate.isSelJpsiToEERich() <= 0) {
            continue;
          }
          if (selectedTofRich && candidate.isSelJpsiToEETofRich() <= 0) {
            continue;
          }
        }
      }

      if (yCandMax >= 0. && std::abs(hfHelper.yJpsi(candidate)) > yCandMax) {
        continue;
      }
      if (candidate.flagMcMatchRec() == 1 << decayMode) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(mcParticles, candidate.prong0_as<aod::TracksWMc>().mcParticle_as<McParticlesHf>(), o2::constants::physics::Pdg::kJPsi, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT
        registry.fill(HIST("hPtRecSig"), candidate.pt());      // rec. level pT
        registry.fill(HIST("hCPARecSig"), candidate.cpa());
        registry.fill(HIST("hEtaRecSig"), candidate.eta());
        if (modeJpsiToMuMu) {
          registry.fill(HIST("hMassSig"), hfHelper.invMassJpsiToMuMu(candidate), candidate.pt());
        } else {
          registry.fill(HIST("hMassSig"), hfHelper.invMassJpsiToEE(candidate), candidate.pt());
        }
        registry.fill(HIST("hDecLengthSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYSig"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hd0Prong0Sig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1Sig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hd0d0Sig"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hChi2PCASig"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCtSig"), hfHelper.ctJpsi(candidate), candidate.pt());
        registry.fill(HIST("hYSig"), hfHelper.yJpsi(candidate), candidate.pt());
        registry.fill(HIST("hYGenSig"), RecoDecay::y(particleMother.pVector(), o2::constants::physics::MassJPsi), particleMother.pt());

      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa());
        registry.fill(HIST("hEtaRecBg"), candidate.eta());
        if (modeJpsiToMuMu) {
          registry.fill(HIST("hMassBg"), hfHelper.invMassJpsiToMuMu(candidate), candidate.pt());
        } else {
          registry.fill(HIST("hMassBg"), hfHelper.invMassJpsiToEE(candidate), candidate.pt());
        }
        registry.fill(HIST("hDecLengthBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthxyBg"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hd0Prong0Bg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1Bg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hd0d0Bg"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hChi2PCABg"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCtBg"), hfHelper.ctJpsi(candidate), candidate.pt());
        registry.fill(HIST("hYBg"), hfHelper.yJpsi(candidate), candidate.pt());
      }
    }
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (particle.flagMcMatchGen() == 1 << decayMode) {
        if (yCandMax >= 0. && std::abs(RecoDecay::y(particle.pVector(), o2::constants::physics::MassJPsi)) > yCandMax) {
          continue;
        }
        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta());
        registry.fill(HIST("hYGen"), RecoDecay::y(particle.pVector(), o2::constants::physics::MassJPsi), particle.pt());
        // registry.fill(HIST("hPtGenProng0"), particle.daughter0_as<McParticlesHf>().pt(), particle.pt());
        // registry.fill(HIST("hPtGenProng1"), particle.daughter1_as<McParticlesHf>().pt(), particle.pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfTaskJpsi>(cfgc)};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfTaskJpsiMc>(cfgc));
  }
  return workflow;
}
