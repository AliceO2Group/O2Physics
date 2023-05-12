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

/// \file correlatorDplusHadrons.cxx
/// \author Shyam Kumar <shyam.kumar@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::aod::hf_correlation_dplus_hadron;
using namespace o2::analysis::hf_cuts_dplus_to_pi_k_pi;
using namespace o2::constants::math;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies

double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PI / 2.);
}

/// definition of variables for Dplus hadron pairs (in data-like, MC-reco and MC-kine tasks)
const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_dplus_to_pi_k_pi::nBinsPt;
const double efficiencyDmesonDefault[npTBinsMassAndEfficiency] = {};
auto efficiencyDmeson_v = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + npTBinsMassAndEfficiency};

// histogram binning definition
const int massAxisBins = 350;
const double massAxisMin = 1.7;
const double massAxisMax = 2.05;
const int phiAxisBins = 32;
const double phiAxisMin = -o2::constants::math::PI / 2.;
const double phiAxisMax = 3. * o2::constants::math::PI / 2.;
const int yAxisBins = 100;
const double yAxisMin = -2.;
const double yAxisMax = 2.;
const int ptDAxisBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;

using MCParticlesPlus3Prong = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

/// Dplus-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDplusHadrons {
  SliceCache cache;
  Preslice<aod::HfCand3Prong> perCol = aod::hf_cand::collisionId;

  Produces<aod::DplusHadronPair> entryDplusHadronPair;
  Produces<aod::DplusHadronRecoInfo> entryDplusHadronRecoInfo;

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Selection Flag for Dplus"};
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
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{efficiencyDmeson_v}, "Efficiency values for Dplus meson"};

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>> selectedDPlusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>> recoFlagDPlusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi > 0;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "Dplus,Hadron candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "Dplus,Hadron candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "Dplus,Hadron candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng2", "Dplus,Hadron candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "Dplus,Hadron candidates;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEta", "Dplus,Hadron candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "Dplus,Hadron candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "Dplus,Hadron candidates;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPtCandMCRec", "Dplus,Hadron candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0MCRec", "Dplus,Hadron candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1MCRec", "Dplus,Hadron candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng2MCRec", "Dplus,Hadron candidates - MC reco;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusMCRec", "Dplus,Hadron candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEtaMCRec", "Dplus,Hadron candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCRec", "Dplus,Hadron candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCRec", "Dplus,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMCEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMCGen", "Dplus,Hadron particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaMCGen", "Dplus,Hadron particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCGen", "Dplus,Hadron particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCGen", "Dplus,Hadron candidates - MC gen;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hcountDplusHadronPerEvent", "Dplus,Hadron particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}}}};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMassDplus_2D", "Dplus candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDplusData", "Dplus candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassDplusMCRec", "Dplus candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassDplusMCRecSig", "Dplus signal candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDplusMCRecBkg", "Dplus background candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcountDplustriggersMCGen", "Dplus trigger particles - MC gen;;N of trigger Dplus", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  /// Dplus-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const& candidates)
  {
    if (selectedDPlusCandidates.size() > 0) {
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

      auto selectedDPlusCandidatesGrouped = selectedDPlusCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

      for (auto& candidate1 : selectedDPlusCandidatesGrouped) {
        if (yCandMax >= 0. && std::abs(yDplus(candidate1)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
          continue;
        }
        if (candidate1.pt() > ptTrackMax) {
          continue;
        }
        // check decay channel flag for candidate1
        if (!(candidate1.hfflag() & 1 << DecayType::DplusToPiKPi)) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
        }
        // fill invariant mass plots and generic info from all Dplus candidates
        registry.fill(HIST("hMassDplus_2D"), invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassDplusData"), invMassDplusToPiKPi(candidate1), efficiencyWeight);
        registry.fill(HIST("hPtCand"), candidate1.pt());
        registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
        registry.fill(HIST("hPtProng2"), candidate1.ptProng2());
        registry.fill(HIST("hEta"), candidate1.eta());
        registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate1.phi(), -o2::constants::math::PI / 2.));
        registry.fill(HIST("hY"), yDplus(candidate1));
        registry.fill(HIST("hSelectionStatus"), candidate1.isSelDplusToPiKPi());

        // Dplus-Hadron correlation dedicated section
        // if the candidate is a Dplus, search for Hadrons and evaluate correlations
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (track.pt() < ptTrackMin) {
            continue;
          }
          if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
            continue; // Remove secondary tracks
          }
          // Removing Dplus daughters by checking track indices
          if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex()) || (candidate1.prong2Id() == track.globalIndex())) {
            continue;
          }
          entryDplusHadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                               track.eta() - candidate1.eta(),
                               candidate1.pt(),
                               track.pt());
          entryDplusHadronRecoInfo(invMassDplusToPiKPi(candidate1), 0);
        } // Hadron Tracks loop
      }   // end outer Dplus loop
    }
  }

  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processData, "Process data", false);

  /// Dplus-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec> const& candidates)
  {
    if (recoFlagDPlusCandidates.size() > 0) {
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

      auto selectedDPlusCandidatesGroupedMC = recoFlagDPlusCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      // MC reco level
      bool flagDplusSignal = false;
      for (auto& candidate1 : selectedDPlusCandidatesGroupedMC) {
        // check decay channel flag for candidate1
        if (!(candidate1.hfflag() & 1 << DecayType::DplusToPiKPi)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(yDplus(candidate1)) > yCandMax) {
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

        if (std::abs(candidate1.flagMcMatchRec()) == 1 << DecayType::DplusToPiKPi) {
          // fill per-candidate distributions from Dplus true candidates
          registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
          registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
          registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
          registry.fill(HIST("hPtProng2MCRec"), candidate1.ptProng2());
          registry.fill(HIST("hEtaMCRec"), candidate1.eta());
          registry.fill(HIST("hPhiMCRec"), RecoDecay::constrainAngle(candidate1.phi(), -o2::constants::math::PI / 2.));
          registry.fill(HIST("hYMCRec"), yDplus(candidate1));
          registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelDplusToPiKPi());
        }
        // fill invariant mass plots from Dplus signal and background candidates
        registry.fill(HIST("hMassDplusMCRec"), invMassDplusToPiKPi(candidate1), efficiencyWeight);
        if (candidate1.flagMcMatchRec() == 1 << DecayType::DplusToPiKPi) { // also matched as Dplus
          registry.fill(HIST("hMassDplusMCRecSig"), invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassDplusMCRecBkg"), invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        }

        // Dplus-Hadron correlation dedicated section
        // if the candidate is selected as Dplus, search for Hadron and evaluate correlations
        flagDplusSignal = candidate1.flagMcMatchRec() == 1 << DecayType::DplusToPiKPi;
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (track.pt() < ptTrackMin) {
            continue;
          }
          // Removing Dplus daughters by checking track indices
          if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex()) || (candidate1.prong2Id() == track.globalIndex())) {
            continue;
          }
          if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
            continue; // Remove secondary tracks
          }
          entryDplusHadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                               track.eta() - candidate1.eta(),
                               candidate1.pt(),
                               track.pt());
          entryDplusHadronRecoInfo(invMassDplusToPiKPi(candidate1), static_cast<int>(flagDplusSignal));
        } // end inner loop (Tracks)

      } // end outer Dplus loop
    }
  }

  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcRec, "Process MC Reco mode", true);
  /// Dplus-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::McCollision const& mccollision, MCParticlesPlus3Prong const& particlesMC)
  {
    int counterDplusHadron = 0;
    registry.fill(HIST("hMCEvtCount"), 0);
    // MC gen level
    for (auto& particle1 : particlesMC) {
      // check if the particle is Dplus  (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != pdg::Code::kDPlus) {
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
      registry.fill(HIST("hPhiMCGen"), RecoDecay::constrainAngle(particle1.phi(), -o2::constants::math::PI / 2.));
      registry.fill(HIST("hYMCGen"), yD);
      counterDplusHadron++;
      // Dplus Hadron correlation dedicated section
      // if it's a Dplus particle, search for Hadron and evaluate correlations
      if (std::abs(particle1.pdgCode()) != pdg::Code::kDPlus) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hcountDplustriggersMCGen"), 0, particle1.pt()); // to count trigger Dplus for normalisation)
      for (auto& particle2 : particlesMC) {
        if (std::abs(particle2.eta()) > etaTrackMax) {
          continue;
        }
        if (particle2.pt() < ptTrackMin) {
          continue;
        }

        if ((std::abs(particle2.pdgCode()) != 11) && (std::abs(particle2.pdgCode()) != 13) && (std::abs(particle2.pdgCode()) != 211) && (std::abs(particle2.pdgCode()) != 321) && (std::abs(particle2.pdgCode()) != 2212)) {
          continue;
        }
        entryDplusHadronPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                             particle2.eta() - particle1.eta(),
                             particle1.pt(),
                             particle2.pt());

      } // end inner loop
    }   // end outer loop
    registry.fill(HIST("hcountDplusHadronPerEvent"), counterDplusHadron);
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDplusHadrons>(cfgc)};
}
