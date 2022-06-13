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

/// \file HfCorrelatorDplusDminus.cxx
/// \brief Dplus-Dminus correlator task - data-like, MC-reco and MC-kine analyses. For ULS and LS pairs
///
/// \author Fabio Colamaria <fabio.colamaria@ba.infn.it>, INFN Bari

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_correlation_ddbar;
using namespace o2::analysis::hf_cuts_dplus_topikpi;
using namespace o2::constants::math;

#include "Framework/runDataProcessing.h"

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PI / 2.);
}

/// definition of variables for DplusDminus pairs vs eta acceptance studies (hDDbarVsEtaCut, in data-like, MC-reco and MC-kine tasks)
const double maxEtaCut = 5.;
const double ptThresholdForMaxEtaCut = 10.;
const double incrementEtaCut = 0.1;
const double incrementPtThreshold = 0.5;
const double epsilon = 1E-5;

const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_dplus_topikpi::npTBins;
const double efficiencyDmesonDefault[npTBinsMassAndEfficiency] = {};
auto efficiencyDmeson_v = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + npTBinsMassAndEfficiency};

// histogram binning definition
const int massAxisBins = 120;
const double massAxisMin = 1.5848;
const double massAxisMax = 2.1848;
const int phiAxisBins = 32;
const double phiAxisMin = 0.;
const double phiAxisMax = 2. * o2::constants::math::PI;
const int yAxisBins = 100;
const double yAxisMin = -5.;
const double yAxisMax = 5.;
const int ptDAxisBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;

using MCParticlesPlus2Prong = soa::Join<aod::McParticles, aod::HfCandProng2MCGen>;
using MCParticlesPlus3Prong = soa::Join<aod::McParticles, aod::HfCandProng3MCGen>;

struct HfCorrelatorDplusDminus {
  Produces<aod::DDbarPair> entryDplusDminusPair;
  Produces<aod::DDbarRecoInfo> entryDplusDminusRecoInfo;

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassDplus for trigger normalisation (S*0.955), and hMass2DCorrelationPairs (in final task) for 2D-sideband-subtraction purposes
    {{"hPtCand", "Dplus,Dminus candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "Dplus,Dminus candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "Dplus,Dminus candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng2", "Dplus,Dminus candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "Dplus,Dminus candidates;selection status;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
     {"hEta", "Dplus,Dminus candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "Dplus,Dminus candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "Dplus,Dminus candidates;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hDDbarVsEtaCut", "Dplus,Dminus pairs vs #eta cut;#eta_{max};candidates #it{p}_{T} threshold (GeV/#it{c});entries", {HistType::kTH2F, {{(int)(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {(int)(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}},
     {"hPtCandMCRec", "Dplus,Dminus candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0MCRec", "Dplus,Dminus candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1MCRec", "Dplus,Dminus candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng2MCRec", "Dplus,Dminus candidates - MC reco;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusMCRec", "Dplus,Dminus candidates - MC reco;selection status;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
     {"hEtaMCRec", "Dplus,Dminus candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCRec", "Dplus,Dminus candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCRec", "Dplus,Dminus candidates - MC reco;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMCEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMCGen", "Dplus,Dminus particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaMCGen", "Dplus,Dminus particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCGen", "Dplus,Dminus particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCGen", "Dplus,Dminus candidates - MC gen;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hCountDplusDminusPerEvent", "Dplus,Dminus particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hDDbarVsDaughterEtaCut", "Dplus,Dminus pairs vs #eta cut on D daughters;#eta_{max};candidates #it{p}_{T} threshold (GeV/#it{c});entries", {HistType::kTH2F, {{(int)(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {(int)(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}},
     {"hCountCCbarPerEvent", "c,cbar particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountCCbarPerEventBeforeEtaCut", "c,cbar particles - MC gen;Number per event pre #eta cut;entries", {HistType::kTH1F, {{20, 0., 20.}}}}}};

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Selection Flag for Dplus,Dminus"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<double> cutPtCandMin{"cutPtCandMin", -1., "min. cand. pT"};
  Configurable<std::vector<double>> bins{"ptBinsForMassAndEfficiency", std::vector<double>{o2::analysis::hf_cuts_dplus_topikpi::pTBins_v}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{efficiencyDmeson_v}, "Efficiency values for Dplus meson"};
  Configurable<int> flagApplyEfficiency{"efficiencyFlagD", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hMass", "Dplus,Dminus candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDplus", "Dplus,Dminus candidates;inv. mass Dplus only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDminus", "Dplus,Dminus candidates;inv. mass Dminus only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDplusMCRecSig", "Dplus signal candidates - MC reco;inv. mass D+ only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDminusMCRecSig", "Dminus signal candidates - MC reco;inv. mass D- only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDplusMCRecBkg", "Dplus background candidates - MC reco;inv. mass D+ only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDminusMCRecBkg", "Dminus background candidates - MC reco;inv. mass D- only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCountDplustriggersMCGen", "Dplus trigger particles - MC gen;;N of trigger D0", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCountCtriggersMCGen", "c trigger particles - MC gen;;N of trigger c quark", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus);

  /// Dplus-Dminus correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCandProng3, aod::HFSelDplusToPiKPiCandidate>> const& candidates, aod::BigTracks const& bigtracks)
  {
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < -4.0 || track.eta() > 4.0) {
          continue;
        }
        if (abs(track.dcaXY()) > 0.0025 || abs(track.dcaZ()) > 0.0025) {
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

    for (auto& candidate1 : candidates) {
      if (cutYCandMax >= 0. && std::abs(YDPlus(candidate1)) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && candidate1.pt() < cutPtCandMin) {
        continue;
      }
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << DecayType::DPlusToPiKPi)) { // probably dummy since already selected? not sure...
        continue;
      }

      double efficiencyWeight = 1.;
      if (flagApplyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }

      int outerParticleSign = 1; // Dplus
      auto outerSecondTrack = candidate1.index1_as<aod::BigTracks>();
      if (outerSecondTrack.sign() == 1) {
        outerParticleSign = -1; // Dminus (second daughter track is positive)
      }

      // fill invariant mass plots and generic info from all Dplus/Dminus candidates
      if (outerParticleSign == 1) {
        registry.fill(HIST("hMass"), InvMassDPlus(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassDplus"), InvMassDPlus(candidate1), candidate1.pt(), efficiencyWeight);
      } else {
        registry.fill(HIST("hMass"), InvMassDPlus(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassDminus"), InvMassDPlus(candidate1), candidate1.pt(), efficiencyWeight);
      }
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hPtProng2"), candidate1.ptProng2());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), YDPlus(candidate1));
      registry.fill(HIST("hSelectionStatus"), candidate1.isSelDplusToPiKPi());

      // D-Dbar correlation dedicated section
      // if the candidate is a Dplus, search for Dminus and evaluate correlations
      if (outerParticleSign != 1) {
        continue;
      }
      for (auto& candidate2 : candidates) {
        // check decay channel flag for candidate2
        if (!(candidate2.hfflag() & 1 << DecayType::DPlusToPiKPi)) { // probably dummy since already selected? not sure...
          continue;
        }
        auto innerSecondTrack = candidate2.index1_as<aod::BigTracks>();
        if (innerSecondTrack.sign() != 1) { // keep only Dminus (with second daughter track positive)
          continue;
        }
        if (cutYCandMax >= 0. && std::abs(YDPlus(candidate2)) > cutYCandMax) {
          continue;
        }
        if (cutPtCandMin >= 0. && candidate2.pt() < cutPtCandMin) {
          continue;
        }
        entryDplusDminusPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                             candidate2.eta() - candidate1.eta(),
                             candidate1.pt(),
                             candidate2.pt());
        entryDplusDminusRecoInfo(InvMassDPlus(candidate1),
                                 InvMassDPlus(candidate2),
                                 0);
        double etaCut = 0.;
        double ptCut = 0.;
        do { // fill pairs vs etaCut plot
          ptCut = 0.;
          etaCut += incrementEtaCut;
          do { // fill pairs vs etaCut plot
            if (std::abs(candidate1.eta()) < etaCut && std::abs(candidate2.eta()) < etaCut && candidate1.pt() > ptCut && candidate2.pt() > ptCut) {
              registry.fill(HIST("hDDbarVsEtaCut"), etaCut - epsilon, ptCut + epsilon);
            }
            ptCut += incrementPtThreshold;
          } while (ptCut < ptThresholdForMaxEtaCut - epsilon);
        } while (etaCut < maxEtaCut - epsilon);
      } // end inner loop (Dminus)
    }   // end outer loop (Dplus)
  }

  PROCESS_SWITCH(HfCorrelatorDplusDminus, processData, "Process data", false);

  /// Dplus-Dminus correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCandProng3, aod::HFSelDplusToPiKPiCandidate, aod::HfCandProng3MCRec>> const& candidates, aod::BigTracks const& bigtracks)
  {
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < -4.0 || track.eta() > 4.0) {
          continue;
        }
        if (abs(track.dcaXY()) > 0.0025 || abs(track.dcaZ()) > 0.0025) {
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

    // MC reco level
    bool flagDplusSignal = false;
    bool flagDminusSignal = false;
    for (auto& candidate1 : candidates) {
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << DecayType::DPlusToPiKPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YDPlus(candidate1)) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && candidate1.pt() < cutPtCandMin) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (flagApplyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }

      int outerParticleSign = 1; // Dplus
      auto outerSecondTrack = candidate1.index1_as<aod::BigTracks>();
      if (outerSecondTrack.sign() == 1) {
        outerParticleSign = -1; // Dminus (second daughter track is positive)
      }
      if (std::abs(candidate1.flagMCMatchRec()) == 1 << DecayType::DPlusToPiKPi) {
        // fill invariant mass plots and per-candidate distributions from Dplus/Dminus signal candidates
        if (outerParticleSign == 1) { // reco and matched as Dplus
          registry.fill(HIST("hMassDplusMCRecSig"), InvMassDPlus(candidate1), candidate1.pt(), efficiencyWeight);
        } else { // reco and matched as Dminus
          registry.fill(HIST("hMassDminusMCRecSig"), InvMassDPlus(candidate1), candidate1.pt(), efficiencyWeight);
        }
        registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
        registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
        registry.fill(HIST("hPtProng2MCRec"), candidate1.ptProng2());
        registry.fill(HIST("hEtaMCRec"), candidate1.eta());
        registry.fill(HIST("hPhiMCRec"), candidate1.phi());
        registry.fill(HIST("hYMCRec"), YDPlus(candidate1));
        registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelDplusToPiKPi());
      } else {
        // fill invariant mass plots from Dplus/Dminus background candidates
        if (outerParticleSign == 1) { // reco as Dplus
          registry.fill(HIST("hMassDplusMCRecBkg"), InvMassDPlus(candidate1), candidate1.pt(), efficiencyWeight);
        } else { // matched as Dminus
          registry.fill(HIST("hMassDminusMCRecBkg"), InvMassDPlus(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }

      // D-Dbar correlation dedicated section
      if (outerParticleSign == -1) {
        continue; // reject Dminus in outer loop
      }
      flagDplusSignal = std::abs(candidate1.flagMCMatchRec()) == 1 << DecayType::DPlusToPiKPi; // flagDplusSignal 'true' if candidate1 matched to Dplus
      for (auto& candidate2 : candidates) {
        if (!(candidate2.hfflag() & 1 << DecayType::DPlusToPiKPi)) { // check decay channel flag for candidate2
          continue;
        }
        auto innerSecondTrack = candidate2.index1_as<aod::BigTracks>();
        if (innerSecondTrack.sign() != 1) { // keep only Dminus (with second daughter track positive)
          continue;
        }
        flagDminusSignal = std::abs(candidate2.flagMCMatchRec()) == 1 << DecayType::DPlusToPiKPi; // flagDminusSignal 'true' if candidate2 matched to Dminus
        if (cutYCandMax >= 0. && std::abs(YDPlus(candidate2)) > cutYCandMax) {
          continue;
        }
        if (cutPtCandMin >= 0. && candidate2.pt() < cutPtCandMin) {
          continue;
        }
        // choice of options (Dplus/Dminus signal/bkg)
        int pairSignalStatus = 0; // 0 = bkg/bkg, 1 = bkg/ref, 2 = bkg/sig, 3 = ref/bkg, 4 = ref/ref, 5 = ref/sig, 6 = sig/bkg, 7 = sig/ref, 8 = sig/sig. Of course only 0,2,6,8 are relevant for D+D-
        if (flagDplusSignal) {
          pairSignalStatus += 6;
        }
        if (flagDminusSignal) {
          pairSignalStatus += 2;
        }
        entryDplusDminusPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                             candidate2.eta() - candidate1.eta(),
                             candidate1.pt(),
                             candidate2.pt());
        entryDplusDminusRecoInfo(InvMassDPlus(candidate1),
                                 InvMassDPlus(candidate2),
                                 pairSignalStatus);
        double etaCut = 0.;
        double ptCut = 0.;
        do { // fill pairs vs etaCut plot
          ptCut = 0.;
          etaCut += incrementEtaCut;
          do { // fill pairs vs etaCut plot
            if (std::abs(candidate1.eta()) < etaCut && std::abs(candidate2.eta()) < etaCut && candidate1.pt() > ptCut && candidate2.pt() > ptCut) {
              registry.fill(HIST("hDDbarVsEtaCut"), etaCut - epsilon, ptCut + epsilon);
            }
            ptCut += incrementPtThreshold;
          } while (ptCut < ptThresholdForMaxEtaCut - epsilon);
        } while (etaCut < maxEtaCut - epsilon);
      } // end inner loop (Dbars)

    } // end outer loop
  }

  PROCESS_SWITCH(HfCorrelatorDplusDminus, processMcRec, "Process MC Reco mode", true);

  /// Dplus-Dminus correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::McCollision const& mccollision, MCParticlesPlus3Prong const& particlesMC)
  {
    int counterDplusDminus = 0;
    registry.fill(HIST("hMCEvtCount"), 0);
    // MC gen level
    for (auto& particle1 : particlesMC) {
      // check if the particle is Dplus or Dminus (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != pdg::Code::kDPlus) {
        continue;
      }
      double yD = RecoDecay::Y(array{particle1.px(), particle1.py(), particle1.pz()}, RecoDecay::getMassPDG(particle1.pdgCode()));
      if (cutYCandMax >= 0. && std::abs(yD) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && particle1.pt() < cutPtCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandMCGen"), particle1.pt());
      registry.fill(HIST("hEtaMCGen"), particle1.eta());
      registry.fill(HIST("hPhiMCGen"), particle1.phi());
      registry.fill(HIST("hYMCGen"), yD);
      counterDplusDminus++;

      // D-Dbar correlation dedicated section
      // if it's a Dplus particle, search for Dminus and evaluate correlations
      if (particle1.pdgCode() != pdg::Code::kDPlus) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hCountDplustriggersMCGen"), 0, particle1.pt()); // to count trigger Dplus (for normalisation)
      for (auto& particle2 : particlesMC) {
        if (particle2.pdgCode() != -pdg::Code::kDPlus) { // check that inner particle is a Dminus
          continue;
        }
        if (cutYCandMax >= 0. && std::abs(RecoDecay::Y(array{particle2.px(), particle2.py(), particle2.pz()}, RecoDecay::getMassPDG(particle2.pdgCode()))) > cutYCandMax) {
          continue;
        }
        if (cutPtCandMin >= 0. && particle2.pt() < cutPtCandMin) {
          continue;
        }
        entryDplusDminusPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                             particle2.eta() - particle1.eta(),
                             particle1.pt(),
                             particle2.pt());
        entryDplusDminusRecoInfo(1.869,
                                 1.869,
                                 8); // Dummy
        double etaCut = 0.;
        double ptCut = 0.;

        // fill pairs vs etaCut plot
        bool rightDecayChannels = false;
        if ((std::abs(particle1.flagMCMatchGen()) == 1 << DecayType::DPlusToPiKPi) && (std::abs(particle2.flagMCMatchGen()) == 1 << DecayType::DPlusToPiKPi)) {
          rightDecayChannels = true;
        }
        do {
          ptCut = 0.;
          etaCut += incrementEtaCut;
          do { // fill pairs vs etaCut plot
            if (std::abs(particle1.eta()) < etaCut && std::abs(particle2.eta()) < etaCut && particle1.pt() > ptCut && particle2.pt() > ptCut) {
              registry.fill(HIST("hDDbarVsEtaCut"), etaCut - epsilon, ptCut + epsilon);
            }
            if (rightDecayChannels) { //fill with D and Dbar daughter particls acceptance checks
              bool candidate1DauInAcc = true;
              bool candidate2DauInAcc = true;
              for (auto& dau : particle1.daughters_as<MCParticlesPlus3Prong>()) {
                if (std::abs(dau.eta()) > etaCut) {
                  candidate1DauInAcc = false;
                  break;
                }
              }
              for (auto& dau : particle2.daughters_as<MCParticlesPlus3Prong>()) {
                if (std::abs(dau.eta()) > etaCut) {
                  candidate2DauInAcc = false;
                  break;
                }
              }
              if (candidate1DauInAcc && candidate2DauInAcc && particle1.pt() > ptCut && particle2.pt() > ptCut) {
                registry.fill(HIST("hDDbarVsDaughterEtaCut"), etaCut - epsilon, ptCut + epsilon);
              }
            }
            ptCut += incrementPtThreshold;
          } while (ptCut < ptThresholdForMaxEtaCut - epsilon);
        } while (etaCut < maxEtaCut - epsilon);
      } // end inner loop
    }   // end outer loop
    registry.fill(HIST("hCountDplusDminusPerEvent"), counterDplusDminus);
  }

  PROCESS_SWITCH(HfCorrelatorDplusDminus, processMcGen, "Process MC Gen mode", false);

  /// c-cbar correlator table builder - for MC gen-level analysis
  void processccbar(aod::McCollision const& mccollision, MCParticlesPlus2Prong const& particlesMC)
  {
    registry.fill(HIST("hMCEvtCount"), 0);
    int counterCCbar = 0, counterCCbarBeforeEtasel = 0;

    // loop over particles at MC gen level
    for (auto& particle1 : particlesMC) {
      if (std::abs(particle1.pdgCode()) != PDG_t::kCharm) { // search c or cbar particles
        continue;
      }
      int partMothPDG = particle1.mothers_as<MCParticlesPlus2Prong>().front().pdgCode();
      //check whether mothers of quark c/cbar are still '4'/'-4' particles - in that case the c/cbar quark comes from its own fragmentation, skip it
      if (partMothPDG == particle1.pdgCode()) {
        continue;
      }
      counterCCbarBeforeEtasel++; // count c or cbar (before kinematic selection)
      double yC = RecoDecay::Y(array{particle1.px(), particle1.py(), particle1.pz()}, RecoDecay::getMassPDG(particle1.pdgCode()));
      if (cutYCandMax >= 0. && std::abs(yC) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && particle1.pt() < cutPtCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandMCGen"), particle1.pt());
      registry.fill(HIST("hEtaMCGen"), particle1.eta());
      registry.fill(HIST("hPhiMCGen"), particle1.phi());
      registry.fill(HIST("hYMCGen"), yC);
      counterCCbar++; // count if c or cbar don't come from themselves during fragmentation (after kinematic selection)

      // c-cbar correlation dedicated section
      // if it's c, search for cbar and evaluate correlations.
      if (particle1.pdgCode() != PDG_t::kCharm) {
        continue;
      }
      registry.fill(HIST("hCountCtriggersMCGen"), 0, particle1.pt()); // to count trigger c quark (for normalisation)

      for (auto& particle2 : particlesMC) {
        if (particle2.pdgCode() != PDG_t::kCharmBar) {
          continue;
        }
        if (cutYCandMax >= 0. && std::abs(RecoDecay::Y(array{particle2.px(), particle2.py(), particle2.pz()}, RecoDecay::getMassPDG(particle2.pdgCode()))) > cutYCandMax) {
          continue;
        }
        if (cutPtCandMin >= 0. && particle2.pt() < cutPtCandMin) {
          continue;
        }
        //check whether mothers of quark cbar (from associated loop) are still '-4' particles - in that case the cbar quark comes from its own fragmentation, skip it
        if (particle2.mothers_as<MCParticlesPlus2Prong>().front().pdgCode() == PDG_t::kCharmBar) {
          continue;
        }
        entryDplusDminusPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                             particle2.eta() - particle1.eta(),
                             particle1.pt(),
                             particle2.pt());
        entryDplusDminusRecoInfo(1.869,
                                 1.869,
                                 8); // Dummy
      }                              // end inner loop
    }                                // end outer loop
    registry.fill(HIST("hCountCCbarPerEvent"), counterCCbar);
    registry.fill(HIST("hCountCCbarPerEventBeforeEtaCut"), counterCCbarBeforeEtasel);
  }

  PROCESS_SWITCH(HfCorrelatorDplusDminus, processccbar, "Process ccbar pairs", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDplusDminus>(cfgc)};
}
