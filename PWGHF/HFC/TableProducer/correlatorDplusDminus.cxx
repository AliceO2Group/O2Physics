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

/// \file correlatorDplusDminus.cxx
/// \brief Dplus-Dminus correlator task - data-like, MC-reco and MC-kine analyses. For ULS and LS pairs
///
/// \author Fabio Colamaria <fabio.colamaria@ba.infn.it>, INFN Bari

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PIHalf);
}

/// definition of variables for DplusDminus pairs vs eta acceptance studies (hDDbarVsEtaCut, in data-like, MC-reco and MC-kine tasks)
const double maxEtaCut = 5.;
const double ptThresholdForMaxEtaCut = 10.;
const double incrementEtaCut = 0.1;
const double incrementPtThreshold = 0.5;
const double epsilon = 1E-5;

const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_dplus_to_pi_k_pi::NBinsPt;
const double efficiencyDmesonDefault[npTBinsMassAndEfficiency] = {};
const auto efficiencyDmesonV = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + npTBinsMassAndEfficiency};

// histogram binning definition
const int massAxisBins = 120;
const double massAxisMin = 1.5848;
const double massAxisMax = 2.1848;
const int phiAxisBins = 32;
const double phiAxisMin = 0.;
const double phiAxisMax = o2::constants::math::TwoPI;
const int yAxisBins = 100;
const double yAxisMin = -5.;
const double yAxisMax = 5.;
const int ptDAxisBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;

using McParticlesPlus2Prong = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
using McParticlesPlus3Prong = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

struct HfCorrelatorDplusDminus {
  Produces<aod::DDbarPair> entryDplusDminusPair;
  Produces<aod::DDbarRecoInfo> entryDplusDminusRecoInfo;

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Selection Flag for Dplus,Dminus"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<double> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
  Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{efficiencyDmesonV}, "Efficiency values for Dplus meson"};

  HfHelper hfHelper;
  SliceCache cache;

  Preslice<aod::HfCand3Prong> perCol = aod::hf_cand::collisionId;

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>> selectedDPlusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>> selectedDPlusCandidatesMC = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

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
     {"hDDbarVsEtaCut", "Dplus,Dminus pairs vs #eta cut;#eta_{max};candidates #it{p}_{T} threshold (GeV/#it{c});entries", {HistType::kTH2F, {{static_cast<int>(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {static_cast<int>(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}},
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
     {"hDDbarVsDaughterEtaCut", "Dplus,Dminus pairs vs #eta cut on D daughters;#eta_{max};candidates #it{p}_{T} threshold (GeV/#it{c});entries", {HistType::kTH2F, {{static_cast<int>(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {static_cast<int>(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}},
     {"hCountCCbarPerEvent", "c,cbar particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountCCbarPerEventBeforeEtaCut", "c,cbar particles - MC gen;Number per event pre #eta cut;entries", {HistType::kTH1F, {{20, 0., 20.}}}}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
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

  /// Dplus-Dminus correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(aod::Collision const& collision,
                   aod::TracksWDca const& tracks,
                   soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const&)
  {
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < -4.0 || track.eta() > 4.0) {
          continue;
        }
        if (std::abs(track.dcaXY()) > 0.0025 || std::abs(track.dcaZ()) > 0.0025) {
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

    for (const auto& candidate1 : selectedDPlusCandidatesGrouped) {
      if (yCandMax >= 0. && std::abs(hfHelper.yDplus(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }
      // check decay channel flag for candidate1
      if ((candidate1.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi) == 0) { // probably dummy since already selected? not sure...
        continue;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
      }

      int outerParticleSign = 1; // Dplus
      auto outerSecondTrack = candidate1.prong1_as<aod::TracksWDca>();
      if (outerSecondTrack.sign() == 1) {
        outerParticleSign = -1; // Dminus (second daughter track is positive)
      }

      // fill invariant mass plots and generic info from all Dplus/Dminus candidates
      if (outerParticleSign == 1) {
        registry.fill(HIST("hMass"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassDplus"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
      } else {
        registry.fill(HIST("hMass"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassDminus"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
      }
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hPtProng2"), candidate1.ptProng2());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), hfHelper.yDplus(candidate1));
      registry.fill(HIST("hSelectionStatus"), candidate1.isSelDplusToPiKPi());

      // D-Dbar correlation dedicated section
      // if the candidate is a Dplus, search for Dminus and evaluate correlations
      if (outerParticleSign != 1) {
        continue;
      }
      for (const auto& candidate2 : selectedDPlusCandidatesGrouped) {
        // check decay channel flag for candidate2
        if ((candidate2.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi) == 0) { // probably dummy since already selected? not sure...
          continue;
        }
        auto innerSecondTrack = candidate2.prong1_as<aod::TracksWDca>();
        if (innerSecondTrack.sign() != 1) { // keep only Dminus (with second daughter track positive)
          continue;
        }
        if (yCandMax >= 0. && std::abs(hfHelper.yDplus(candidate2)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate2.pt() < ptCandMin) {
          continue;
        }
        entryDplusDminusPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                             candidate2.eta() - candidate1.eta(),
                             candidate1.pt(),
                             candidate2.pt());
        entryDplusDminusRecoInfo(hfHelper.invMassDplusToPiKPi(candidate1),
                                 hfHelper.invMassDplusToPiKPi(candidate2),
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
    } // end outer loop (Dplus)
  }
  PROCESS_SWITCH(HfCorrelatorDplusDminus, processData, "Process data", false);

  /// Dplus-Dminus correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(aod::Collision const& collision,
                    aod::TracksWDca const& tracks,
                    soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec> const&)
  {
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < -4.0 || track.eta() > 4.0) {
          continue;
        }
        if (std::abs(track.dcaXY()) > 0.0025 || std::abs(track.dcaZ()) > 0.0025) {
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

    auto selectedDPlusCandidatesGroupedMC = selectedDPlusCandidatesMC->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

    // MC reco level
    bool flagDplusSignal = false;
    bool flagDminusSignal = false;
    for (const auto& candidate1 : selectedDPlusCandidatesGroupedMC) {
      // check decay channel flag for candidate1
      if ((candidate1.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi) == 0) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(hfHelper.yDplus(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
      }

      int outerParticleSign = 1; // Dplus
      auto outerSecondTrack = candidate1.prong1_as<aod::TracksWDca>();
      if (outerSecondTrack.sign() == 1) {
        outerParticleSign = -1; // Dminus (second daughter track is positive)
      }
      if (std::abs(candidate1.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
        // fill invariant mass plots and per-candidate distributions from Dplus/Dminus signal candidates
        if (outerParticleSign == 1) { // reco and matched as Dplus
          registry.fill(HIST("hMassDplusMCRecSig"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else { // reco and matched as Dminus
          registry.fill(HIST("hMassDminusMCRecSig"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        }
        registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
        registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
        registry.fill(HIST("hPtProng2MCRec"), candidate1.ptProng2());
        registry.fill(HIST("hEtaMCRec"), candidate1.eta());
        registry.fill(HIST("hPhiMCRec"), candidate1.phi());
        registry.fill(HIST("hYMCRec"), hfHelper.yDplus(candidate1));
        registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelDplusToPiKPi());
      } else {
        // fill invariant mass plots from Dplus/Dminus background candidates
        if (outerParticleSign == 1) { // reco as Dplus
          registry.fill(HIST("hMassDplusMCRecBkg"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else { // matched as Dminus
          registry.fill(HIST("hMassDminusMCRecBkg"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }

      // D-Dbar correlation dedicated section
      if (outerParticleSign == -1) {
        continue; // reject Dminus in outer loop
      }
      flagDplusSignal = std::abs(candidate1.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi; // flagDplusSignal 'true' if candidate1 matched to Dplus
      for (const auto& candidate2 : selectedDPlusCandidatesGroupedMC) {
        if ((candidate2.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi) == 0) { // check decay channel flag for candidate2
          continue;
        }
        auto innerSecondTrack = candidate2.prong1_as<aod::TracksWDca>();
        if (innerSecondTrack.sign() != 1) { // keep only Dminus (with second daughter track positive)
          continue;
        }
        flagDminusSignal = std::abs(candidate2.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi; // flagDminusSignal 'true' if candidate2 matched to Dminus
        if (yCandMax >= 0. && std::abs(hfHelper.yDplus(candidate2)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate2.pt() < ptCandMin) {
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
        entryDplusDminusRecoInfo(hfHelper.invMassDplusToPiKPi(candidate1),
                                 hfHelper.invMassDplusToPiKPi(candidate2),
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
  void processMcGen(aod::McCollision const&,
                    McParticlesPlus3Prong const& mcParticles)
  {
    int counterDplusDminus = 0;
    registry.fill(HIST("hMCEvtCount"), 0);
    // MC gen level
    for (const auto& particle1 : mcParticles) {
      // check if the particle is Dplus or Dminus (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != Pdg::kDPlus) {
        continue;
      }
      double const yD = RecoDecay::y(particle1.pVector(), MassDPlus);
      if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandMCGen"), particle1.pt());
      registry.fill(HIST("hEtaMCGen"), particle1.eta());
      registry.fill(HIST("hPhiMCGen"), particle1.phi());
      registry.fill(HIST("hYMCGen"), yD);
      counterDplusDminus++;

      // D-Dbar correlation dedicated section
      // if it's a Dplus particle, search for Dminus and evaluate correlations
      if (particle1.pdgCode() != Pdg::kDPlus) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hCountDplustriggersMCGen"), 0, particle1.pt()); // to count trigger Dplus (for normalisation)
      for (const auto& particle2 : mcParticles) {
        if (particle2.pdgCode() != -Pdg::kDPlus) { // check that inner particle is a Dminus
          continue;
        }
        if (yCandMax >= 0. && std::abs(RecoDecay::y(particle2.pVector(), MassDPlus)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle2.pt() < ptCandMin) {
          continue;
        }
        entryDplusDminusPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                             particle2.eta() - particle1.eta(),
                             particle1.pt(),
                             particle2.pt());
        entryDplusDminusRecoInfo(MassDPlus,
                                 MassDPlus,
                                 8); // Dummy
        double etaCut = 0.;
        double ptCut = 0.;

        // fill pairs vs etaCut plot
        bool rightDecayChannels = false;
        if ((std::abs(particle1.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) && (std::abs(particle2.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) {
          rightDecayChannels = true;
        }
        do {
          ptCut = 0.;
          etaCut += incrementEtaCut;
          do { // fill pairs vs etaCut plot
            if (std::abs(particle1.eta()) < etaCut && std::abs(particle2.eta()) < etaCut && particle1.pt() > ptCut && particle2.pt() > ptCut) {
              registry.fill(HIST("hDDbarVsEtaCut"), etaCut - epsilon, ptCut + epsilon);
            }
            if (rightDecayChannels) { // fill with D and Dbar daughter particls acceptance checks
              bool candidate1DauInAcc = true;
              bool candidate2DauInAcc = true;
              for (const auto& dau : particle1.daughters_as<McParticlesPlus3Prong>()) {
                if (std::abs(dau.eta()) > etaCut) {
                  candidate1DauInAcc = false;
                  break;
                }
              }
              for (const auto& dau : particle2.daughters_as<McParticlesPlus3Prong>()) {
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
    } // end outer loop
    registry.fill(HIST("hCountDplusDminusPerEvent"), counterDplusDminus);
  }
  PROCESS_SWITCH(HfCorrelatorDplusDminus, processMcGen, "Process MC Gen mode", false);

  /// c-cbar correlator table builder - for MC gen-level analysis
  void processCCbar(aod::McCollision const&,
                    McParticlesPlus2Prong const& mcParticles)
  {
    registry.fill(HIST("hMCEvtCount"), 0);
    int counterCCbar = 0, counterCCbarBeforeEtasel = 0;

    // loop over particles at MC gen level
    for (const auto& particle1 : mcParticles) {
      if (std::abs(particle1.pdgCode()) != PDG_t::kCharm) { // search c or cbar particles
        continue;
      }
      int const partMothPDG = particle1.mothers_as<McParticlesPlus2Prong>().front().pdgCode();
      // check whether mothers of quark c/cbar are still '4'/'-4' particles - in that case the c/cbar quark comes from its own fragmentation, skip it
      if (partMothPDG == particle1.pdgCode()) {
        continue;
      }
      counterCCbarBeforeEtasel++; // count c or cbar (before kinematic selection)
      double const yC = RecoDecay::y(particle1.pVector(), MassCharm);
      if (yCandMax >= 0. && std::abs(yC) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
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

      for (const auto& particle2 : mcParticles) {
        if (particle2.pdgCode() != PDG_t::kCharmBar) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(RecoDecay::y(particle2.pVector(), MassCharmBar)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle2.pt() < ptCandMin) {
          continue;
        }
        // check whether mothers of quark cbar (from associated loop) are still '-4' particles - in that case the cbar quark comes from its own fragmentation, skip it
        if (particle2.mothers_as<McParticlesPlus2Prong>().front().pdgCode() == PDG_t::kCharmBar) {
          continue;
        }
        entryDplusDminusPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                             particle2.eta() - particle1.eta(),
                             particle1.pt(),
                             particle2.pt());
        entryDplusDminusRecoInfo(1.869,
                                 1.869,
                                 8); // Dummy
      } // end inner loop
    } // end outer loop
    registry.fill(HIST("hCountCCbarPerEvent"), counterCCbar);
    registry.fill(HIST("hCountCCbarPerEventBeforeEtaCut"), counterCCbarBeforeEtasel);
  }
  PROCESS_SWITCH(HfCorrelatorDplusDminus, processCCbar, "Process c-cbar pairs", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDplusDminus>(cfgc)};
}
