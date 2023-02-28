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

/// \file correlatorDsHadrons.cxx
/// \brief Ds-Hadrons correlator task - data-like, MC-reco and MC-Gen analyses
/// \author Grazia Luparello <grazia.luparello@cern.ch>
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::aod::hf_correlation_ds_hadron;
using namespace o2::analysis::hf_cuts_ds_to_k_k_pi;
using namespace o2::constants::math;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

/// definition of variables for Ds hadron pairs (in data-like, MC-reco and MC-kine tasks)
const int nBinsPtMassAndEfficiency = o2::analysis::hf_cuts_ds_to_k_k_pi::nBinsPt;
const double efficiencyDmesonDefault[nBinsPtMassAndEfficiency] = {};
auto vecEfficiencyDmeson = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + nBinsPtMassAndEfficiency};

// histogram axes definition
AxisSpec axisMassD = {300, 1.7, 2.2, ""};
AxisSpec axisPhi = {128, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf, ""};
AxisSpec axisEta = {100, -2., 2., ""};
AxisSpec axisY = {100, -2., 2., ""};
AxisSpec axisPtD = {180, 0., 36., ""};

using MCParticlesPlus3Prong = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

/// Ds-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDsHadrons {
  Produces<aod::DsHadronPair> entryDsHadronPair;
  Produces<aod::DsHadronRecoInfo> entryDsHadronRecoInfo;

  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<double> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<double> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<double> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<double> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<double> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<double> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<double> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for Ds meson"};

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>> selectedDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>> recoFlagDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "Ds,Hadron candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng0", "Ds,Hadron candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng1", "Ds,Hadron candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng2", "Ds,Hadron candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hSelectionStatus", "Ds,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hEta", "Ds,Hadron candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}}},
     {"hPhi", "Ds,Hadron candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {axisPhi}}},
     {"hY", "Ds,Hadron candidates;candidate #it{#y};entries", {HistType::kTH1F, {axisY}}},
     {"hPtCandMCRec", "Ds,Hadron candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng0MCRec", "Ds,Hadron candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng1MCRec", "Ds,Hadron candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng2MCRec", "Ds,Hadron candidates - MC reco;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hSelectionStatusMCRec", "Ds,Hadron candidates - MC reco;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hEtaMCRec", "Ds,Hadron candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}}},
     {"hPhiMCRec", "Ds,Hadron candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {axisPhi}}},
     {"hYMCRec", "Ds,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisY}}},
     {"hMCEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMCGen", "Ds,Hadron particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}}},
     {"hEtaMCGen", "Ds,Hadron particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {axisEta}}},
     {"hPhiMCGen", "Ds,Hadron particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {axisPhi}}},
     {"hYMCGen", "Ds,Hadron candidates - MC gen;candidate #it{#y};entries", {HistType::kTH1F, {axisY}}},
     {"hcountDsHadronPerEvent", "Ds,Hadron particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}}}};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMassDs2D", "Ds candidates;inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDsData", "Ds candidates;inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassDsMCRec", "Ds candidates;inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassDsMCRecSig", "Ds signal candidates - MC reco;inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDsMCRecBkg", "Ds background candidates - MC reco;inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcountDstriggersMCGen", "Ds trigger particles - MC gen;;N of trigger Ds", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  /// Ds-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi> const& candidates)
  {
    // number of tracks
    if (selectedDsCandidates.size() > 0) {
      int nTracks = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
            continue;
          }
          nTracks++;
        }
      }

      auto selectedDsCandidatesGrouped = selectedDsCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex());

      for (const auto& candidate1 : selectedDsCandidatesGrouped) {
        if (yCandMax >= 0. && std::abs(yDs(candidate1)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
          continue;
        }
        if (candidate1.pt() > ptTrackMax) {
          continue;
        }
        // check decay channel flag for candidate1
        if (!(candidate1.hfflag() & 1 << DecayType::DsToKKPi)) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
        }

        // fill invariant mass plots and generic info from all Ds candidates
        registry.fill(HIST("hMassDs2D"), invMassDsToKKPi(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassDsData"), invMassDsToKKPi(candidate1), efficiencyWeight);
        registry.fill(HIST("hPtCand"), candidate1.pt());
        registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
        registry.fill(HIST("hPtProng2"), candidate1.ptProng2());
        registry.fill(HIST("hEta"), candidate1.eta());
        registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate1.phi(), -o2::constants::math::PIHalf));
        registry.fill(HIST("hY"), yDs(candidate1));
        registry.fill(HIST("hSelectionStatus"), candidate1.isSelDsToKKPi());

        // Ds-Hadron correlation dedicated section
        // if the candidate is a Ds, search for Hadrons and evaluate correlations
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (track.pt() < ptTrackMin) {
            continue;
          }
          // Remove secondary tracks
          if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
            continue;
          }
          // Removing Ds daughters by checking track indices
          if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex()) || (candidate1.prong2Id() == track.globalIndex())) {
            continue;
          }
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                            track.eta() - candidate1.eta(),
                            candidate1.pt(),
                            track.pt());
          entryDsHadronRecoInfo(invMassDsToKKPi(candidate1), false);
        } // Hadron Tracks loop
      }   // end outer Ds loop
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processData, "Process data", true);

  /// Ds-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec> const& candidates)
  {
    if (recoFlagDsCandidates.size() > 0) {
      int nTracks = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
            continue;
          }
          nTracks++;
        }
      }
      registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
      if (nTracks < multMin || nTracks > multMax) {
        return;
      }
      registry.fill(HIST("hMultiplicity"), nTracks);

      auto selectedDsCandidatesGroupedMC = recoFlagDsCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex());
      // MC reco level
      bool flagDsSignal = false;
      for (const auto& candidate1 : selectedDsCandidatesGroupedMC) {
        // check decay channel flag for candidate1
        if (!(candidate1.hfflag() & 1 << DecayType::DsToKKPi)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(yDs(candidate1)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
          continue;
        }
        if (candidate1.pt() >= ptTrackMax) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
        }

        if (std::abs(candidate1.flagMcMatchRec()) == 1 << DecayType::DsToKKPi) {
          // fill per-candidate distributions from Ds true candidates
          registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
          registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
          registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
          registry.fill(HIST("hPtProng2MCRec"), candidate1.ptProng2());
          registry.fill(HIST("hEtaMCRec"), candidate1.eta());
          registry.fill(HIST("hPhiMCRec"), RecoDecay::constrainAngle(candidate1.phi(), -o2::constants::math::PIHalf));
          registry.fill(HIST("hYMCRec"), yDs(candidate1));
          registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelDsToKKPi());
        }
        // fill invariant mass plots from Ds signal and background candidates
        registry.fill(HIST("hMassDsMCRec"), invMassDsToKKPi(candidate1), efficiencyWeight);
        if (candidate1.flagMcMatchRec() == 1 << DecayType::DsToKKPi) { // also matched as Ds
          registry.fill(HIST("hMassDsMCRecSig"), invMassDsToKKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassDsMCRecBkg"), invMassDsToKKPi(candidate1), candidate1.pt(), efficiencyWeight);
        }

        // Ds-Hadron correlation dedicated section
        // if the candidate is selected as Ds, search for Hadron and evaluate correlations
        flagDsSignal = candidate1.flagMcMatchRec() == 1 << DecayType::DsToKKPi;
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (track.pt() < ptTrackMin) {
            continue;
          }
          // Removing Ds daughters by checking track indices
          if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex()) || (candidate1.prong2Id() == track.globalIndex())) {
            continue;
          }
          if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
            continue; // Remove secondary tracks
          }
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                            track.eta() - candidate1.eta(),
                            candidate1.pt(),
                            track.pt());
          entryDsHadronRecoInfo(invMassDsToKKPi(candidate1), flagDsSignal);
        } // end inner loop (Tracks)

      } // end outer Ds loop
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRec, "Process MC Reco mode", false);

  /// Ds-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::McCollision const&, MCParticlesPlus3Prong const& particlesMc)
  {
    int counterDsHadron = 0;
    registry.fill(HIST("hMCEvtCount"), 0);
    // MC gen level
    for (auto const& particle1 : particlesMc) {
      // check if the particle is Ds  (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != pdg::Code::kDS) {
        continue;
      }
      double yD = RecoDecay::y(array{particle1.px(), particle1.py(), particle1.pz()}, RecoDecay::getMassPDG(particle1.pdgCode()));
      if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandMCGen"), particle1.pt());
      registry.fill(HIST("hEtaMCGen"), particle1.eta());
      registry.fill(HIST("hPhiMCGen"), RecoDecay::constrainAngle(particle1.phi(), -o2::constants::math::PIHalf));
      registry.fill(HIST("hYMCGen"), yD);
      counterDsHadron++;
      // Ds Hadron correlation dedicated section
      // if it's a Ds particle, search for Hadron and evaluate correlations
      if (std::abs(particle1.pdgCode()) != pdg::Code::kDS) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hcountDstriggersMCGen"), 0, particle1.pt()); // to count trigger Ds for normalisation

      for (auto const& particle2 : particlesMc) {
        if (std::abs(particle2.eta()) > etaTrackMax) {
          continue;
        }
        if (particle2.pt() < ptTrackMin) {
          continue;
        }

        if ((std::abs(particle2.pdgCode()) != kElectron) && (std::abs(particle2.pdgCode()) != kMuonMinus) && (std::abs(particle2.pdgCode()) != kPiPlus) && (std::abs(particle2.pdgCode()) != kKPlus) && (std::abs(particle2.pdgCode()) != kProton)) {
          continue;
        }

        entryDsHadronPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                          particle2.eta() - particle1.eta(),
                          particle1.pt(),
                          particle2.pt());

      } // end inner loop
    }   // end outer loop
    registry.fill(HIST("hcountDsHadronPerEvent"), counterDsHadron);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDsHadrons>(cfgc)};
}
