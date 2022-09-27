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

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/Core/HFSelectorCuts.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_bplus;
using namespace o2::analysis::hf_cuts_bplus_tod0pi;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

// string definitions, used for histogram axis labels
const TString stringPt = "#it{p}_{T} (GeV/#it{c})";
const TString BPlusCandTitle = "B+ candidates;";
const TString entries = "entries";
const TString BPlusCandMatch = "B+ candidates (matched);";
const TString BPlusCandUnmatch = "B+ candidates (unmatched);";
const TString mcParticleMatched = "MC particles (matched);";

/// B± analysis task
struct HfTaskBplus {
  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", BPlusCandTitle + "prong 0 #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hPtProng1", BPlusCandTitle + "prong 1 #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hPtCand", BPlusCandTitle + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hCentrality", "centrality;centrality percentile;" + entries, {HistType::kTH1F, {{100, 0., 100.}}}},
     {"hPtRecSig", BPlusCandMatch + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtRecBg", BPlusCandUnmatch + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtGenSig", BPlusCandMatch + "candidate #it{p}_{T}^{gen.} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 10.}}}},
     {"hPtGen", mcParticleMatched + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}}}}};

  Configurable<int> selectionFlagBPlus{"selectionFlagBPlus", 1, "Selection Flag for B+"};
  Configurable<double> cutYCandMax{"cutYCandMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_bplus_tod0pi::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    // registry.add("hMass", BPlusCandTitle + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});centrality", {HistType::kTH3F, {{500, 0., 10.}, {(std::vector<double>)bins, "#it{p}_{T} (GeV/#it{c})"}, {100, 0., 100.}}});
    registry.add("hMass", BPlusCandTitle + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {{500, 0., 10.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLength", BPlusCandTitle + "decay length (cm);" + entries, {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLengthXY", BPlusCandTitle + "decay length xy (cm);" + entries, {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0Prong0", BPlusCandTitle + "prong 0 DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0Prong1", BPlusCandTitle + "prong 1 DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hCPA", BPlusCandTitle + "candidate cosine of pointing angle;" + entries, {HistType::kTH2F, {{110, -1.0, 1.1}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hEta", BPlusCandTitle + "candidate #it{#eta};" + entries, {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hRapidity", BPlusCandTitle + "candidate #it{y};" + entries, {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hImpParErr", BPlusCandTitle + "candidate impact parameter error (cm);" + entries, {HistType::kTH2F, {{100, -1., 1.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLenErr", BPlusCandTitle + "candidate decay length error (cm);" + entries, {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLenXYErr", BPlusCandTitle + "candidate decay length xy error (cm);" + entries, {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0d0", BPlusCandTitle + "candidate product of DCAxy to prim. vertex (cm^{2});" + entries, {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hInvMassD0", BPlusCandTitle + "prong0, D0 inv. mass (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {{500, 0, 5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hEtaGen", mcParticleMatched + "candidate #it{#eta}^{gen};" + entries, {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hYGen", mcParticleMatched + "candidate #it{y}^{gen};" + entries, {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hPtProng0Gen", mcParticleMatched + "prong 0 #it{p}_{T}^{gen} (GeV/#it{c});" + entries, {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hPtProng1Gen", mcParticleMatched + "prong 1 #it{p}_{T}^{gen} (GeV/#it{c});" + entries, {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hYProng0Gen", mcParticleMatched + "prong 0 #it{y}^{gen};" + entries, {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hYProng1Gen", mcParticleMatched + "prong 1 #it{y}^{gen};" + entries, {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hEtaProng0Gen", mcParticleMatched + "prong 0 #it{#eta}^{gen};" + entries, {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hEtaProng1Gen", mcParticleMatched + "prong 1 #it{#eta}^{gen};" + entries, {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hCPARecSig", BPlusCandMatch + "candidate cosine of pointing angle;" + entries, {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hCPARecBg", BPlusCandUnmatch + "candidate cosine of pointing angle;" + entries, {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hCPAxyRecSig", BPlusCandMatch + "candidate CPAxy;" + entries, {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hCPAxyRecBg", BPlusCandUnmatch + "candidate CPAxy;" + entries, {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hEtaRecSig", BPlusCandMatch + "candidate #it{#eta};" + entries, {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hEtaRecBg", BPlusCandUnmatch + "candidate #it{#eta};" + entries, {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hRapidityRecSig", BPlusCandMatch + "candidate #it{y};" + entries, {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hRapidityRecBg", BPlusCandUnmatch + "candidate #it{#y};" + entries, {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hPtProng0RecSig", BPlusCandMatch + "prong 0 #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hPtProng1RecSig", BPlusCandMatch + "prong 1 #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hPtProng0RecBg", BPlusCandUnmatch + "prong 0 #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hPtProng1RecBg", BPlusCandUnmatch + "prong 1 #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hMassRecSig", BPlusCandMatch + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {{300, 4.0, 7.00}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hMassRecBg", BPlusCandUnmatch + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {{300, 4.0, 7.0}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0Prong0RecSig", BPlusCandMatch + "prong 0 DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0Prong1RecSig", BPlusCandMatch + "prong 1 DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0Prong0RecBg", BPlusCandUnmatch + "prong 0 DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0Prong1RecBg", BPlusCandUnmatch + "prong 1 DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLengthRecSig", BPlusCandMatch + "candidate decay length (cm);" + entries, {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLengthXYRecSig", BPlusCandMatch + "candidate decay length xy (cm);" + entries, {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLengthRecBg", BPlusCandUnmatch + "candidate decay length (cm);" + entries, {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLengthXYRecBg", BPlusCandUnmatch + "candidate decay length xy(cm);" + entries, {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLengthNormRecSig", BPlusCandMatch + "candidate normalized decay length (cm);" + entries, {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hDecLengthNormRecBg", BPlusCandUnmatch + "candidate normalized decay length (cm);" + entries, {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0d0RecSig", BPlusCandMatch + "product of DCAxy to prim. vertex (cm^{2});" + entries, {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
    registry.add("hd0d0RecBg", BPlusCandUnmatch + "product of DCAxy to prim. vertex (cm^{2});" + entries, {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)bins, stringPt.Data()}}});
  }

  Partition<soa::Join<aod::HfCandBPlus, aod::HFSelBPlusToD0PiCandidate>> selectedBPlusCandidates = aod::hf_selcandidate_bplus::isSelBPlusToD0Pi >= selectionFlagBPlus;

  void process(aod::Collisions const& collision, soa::Join<aod::HfCandBPlus, aod::HFSelBPlusToD0PiCandidate> const&, soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate> const&, aod::BigTracks const&)
  {

    for (auto& candidate : selectedBPlusCandidates) {
      if (!(candidate.hfflag() & 1 << hf_cand_bplus::DecayType::BPlusToD0Pi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YBPlus(candidate)) > cutYCandMax) {
        continue;
      }

      auto candD0 = candidate.index0_as<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>>();
      auto candPi = candidate.index1_as<aod::BigTracks>();

      registry.fill(HIST("hMass"), InvMassBPlus(candidate), candidate.pt());
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hRapidity"), YBPlus(candidate), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
      if (candPi.sign() > 0) {
        registry.fill(HIST("hInvMassD0"), InvMassD0bar(candD0), candidate.pt());
      } else {
        registry.fill(HIST("hInvMassD0"), InvMassD0(candD0), candidate.pt());
      }
    } // candidate loop
  }   // process

  Partition<soa::Join<aod::HfCandBPlus, aod::HFSelBPlusToD0PiCandidate, aod::HfCandBPMCRec>> selectedBPlusCandidatesMC = aod::hf_selcandidate_bplus::isSelBPlusToD0Pi >= selectionFlagBPlus;

  void processMC(soa::Join<aod::HfCandBPlus, aod::HFSelBPlusToD0PiCandidate, aod::HfCandBPMCRec> const&,
                 soa::Join<aod::McParticles, aod::HfCandBPMCGen> const& particlesMC, aod::BigTracksMC const& tracks, aod::HfCandProng2 const&)
  {
    // MC rec
    for (auto& candidate : selectedBPlusCandidatesMC) {
      if (!(candidate.hfflag() & 1 << hf_cand_bplus::DecayType::BPlusToD0Pi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YBPlus(candidate)) > cutYCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMCMatchRec()) == 1 << hf_cand_bplus::DecayType::BPlusToD0Pi) {

        auto indexMother = RecoDecay::getMother(particlesMC, candidate.index1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandBPMCGen>>(), pdg::Code::kBPlus, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hRapidityRecSig"), YBPlus(candidate), candidate.pt());
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hMassRecSig"), InvMassBPlus(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hd0d0RecSig"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), candidate.pt());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hRapidityRecBg"), YBPlus(candidate), candidate.pt());
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hMassRecBg"), InvMassBPlus(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hd0d0RecBg"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), candidate.pt());
      }
    } // rec

    // MC gen. level
    // Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << hf_cand_bplus::DecayType::BPlusToD0Pi) {

        auto yParticle = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(pdg::Code::kBPlus));
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
    } // gen
  }   // processMC

  PROCESS_SWITCH(HfTaskBplus, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBplus>(cfgc)};
}
