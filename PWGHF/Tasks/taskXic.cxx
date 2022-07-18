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

/// \file taskXic.cxx
/// \brief Ξc± analysis task
/// \note Inspired from taskLc.cxx
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, University and INFN PADOVA
/// \author Anton Alkin <anton.alkin@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::analysis::hf_cuts_xic_topkpi;

#include "Framework/runDataProcessing.h"

/// Ξc± analysis task
struct TaskXic {
  HistogramRegistry registry{
    "registry",
    {
      {"hPtCand", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"hPtRecBg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}} ///
    }};

  Configurable<int> d_selectionFlagXic{"d_selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_xic_topkpi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hMass", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "3-prong candidates;selection status;;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "3-prong candidates;impact parameter error (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "3-prong candidates;decay length error (cm);;entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hMassSig", "Invariant mass (matched);m (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassBg", "Invariant mass (unmatched);m (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaGen", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hDecLengthRecSig", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2RecSig", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtRecSig", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecSig", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecSig", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng2RecSig", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hDecLengthRecBg", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2RecBg", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtRecBg", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng2RecBg", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Partition<soa::Join<aod::HfCandProng3, aod::HFSelXicToPKPiCandidate>> selectedXicCandidates = aod::hf_selcandidate_xic::isSelXicToPKPi >= d_selectionFlagXic || aod::hf_selcandidate_xic::isSelXicToPiKP >= d_selectionFlagXic;

  void process(soa::Join<aod::HfCandProng3, aod::HFSelXicToPKPiCandidate>& candidates)
  {
    for (auto& candidate : selectedXicCandidates) {
      if (!(candidate.hfflag() & 1 << DecayType::XicToPKPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YXic(candidate)) > cutYCandMax) {
        continue;
      }

      if (candidate.isSelXicToPKPi() >= d_selectionFlagXic) {
        registry.fill(HIST("hMass"), InvMassXicToPKPi(candidate), candidate.pt());
      }
      if (candidate.isSelXicToPiKP() >= d_selectionFlagXic) {
        registry.fill(HIST("hMass"), InvMassXicToPiKP(candidate), candidate.pt());
      }

      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0(), candidate.pt());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1(), candidate.pt());
      registry.fill(HIST("hPtProng2"), candidate.ptProng2(), candidate.pt());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0Prong2"), candidate.impactParameter2(), candidate.pt());
      registry.fill(HIST("hCt"), CtXic(candidate), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hSelectionStatus"), candidate.isSelXicToPKPi(), candidate.pt());
      registry.fill(HIST("hSelectionStatus"), candidate.isSelXicToPiKP(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter2(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.chi2PCA(), candidate.pt());
    }
  }

  Partition<soa::Join<aod::HfCandProng3, aod::HFSelXicToPKPiCandidate, aod::HfCandProng3MCRec>> recoFlag3Prong = (aod::hf_selcandidate_xic::isSelXicToPKPi >= d_selectionFlagXic || aod::hf_selcandidate_xic::isSelXicToPiKP >= d_selectionFlagXic);

  void processMC(soa::Join<aod::HfCandProng3, aod::HFSelXicToPKPiCandidate, aod::HfCandProng3MCRec> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandProng3MCGen> const& particlesMC, aod::BigTracksMC const&)
  {
    // MC rec.
    for (auto& candidate : recoFlag3Prong) {
      if (!(candidate.hfflag() & 1 << DecayType::XicToPKPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YXic(candidate)) > cutYCandMax) {
        continue;
      }

      auto mass = 0.;
      auto massC = 0.;

      if (candidate.isSelXicToPKPi() >= d_selectionFlagXic) {
        mass = InvMassXicToPKPi(candidate);
      }
      if (candidate.isSelXicToPiKP() >= d_selectionFlagXic) {
        massC = InvMassXicToPiKP(candidate);
      }

      if (std::abs(candidate.flagMCMatchRec()) == 1 << DecayType::XicToPKPi) {
        // Signal
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.index0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandProng3MCGen>>(), pdg::Code::kXiCPlus, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);

        registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT
        registry.fill(HIST("hPtRecSig"), candidate.pt());      // rec. level pT

        if (mass != 0.) {
          registry.fill(HIST("hMassSig"), mass, candidate.pt());
        }
        if (massC != 0.) {
          registry.fill(HIST("hMassSig"), massC, candidate.pt());
        }

        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hPtProng2RecSig"), candidate.ptProng2(), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hd0Prong2RecSig"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("hCtRecSig"), CtXic(candidate), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());
      } else {
        // Background
        registry.fill(HIST("hPtRecBg"), candidate.pt());

        if (mass != 0.) {
          registry.fill(HIST("hMassBg"), mass, candidate.pt());
        }
        if (massC != 0.) {
          registry.fill(HIST("hMassBg"), massC, candidate.pt());
        }

        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hPtProng2RecBg"), candidate.ptProng2(), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hd0Prong2RecBg"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("hCtRecBg"), CtXic(candidate), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());
      }
    }
    // MC gen.
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::XicToPKPi) {
        if (cutYCandMax >= 0. && std::abs(RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > cutYCandMax) {
          continue;
        }
        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta(), particle.pt());
      }
    }
  }

  PROCESS_SWITCH(TaskXic, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskXic>(cfgc, TaskName{"hf-task-xic"})};
}
