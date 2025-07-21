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

/// \file alice3-correlatorDDbar.cxx
/// \brief D0-D0bar correlator task - data-like and MC-reco analysis performance on ALICE 3 detector.
///
/// \author Fabio Colamaria <fabio.colamaria@ba.infn.it>, INFN Bari

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

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

/// definition of variables for D0D0bar pairs vs eta acceptance studies (hDDbarVsEtaCut, in data-like, MC-reco and MC-kine tasks)
const double maxEtaCut = 5.;
const double ptThresholdForMaxEtaCut = 10.;
const double incrementEtaCut = 0.1;
const double incrementPtThreshold = 0.5;
const double epsilon = 1E-5;

const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_d0_to_pi_k::NBinsPt;
const double efficiencyDmesonDefault[npTBinsMassAndEfficiency] = {};
auto efficiencyDmeson_v = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + npTBinsMassAndEfficiency};

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

struct alice3correlatorddbar {
  SliceCache cache;
  Preslice<aod::Alice3D0Meson> perCol = aod::a3D0meson::collisionId;
  Produces<aod::DDbarPair> entryD0D0barPair;
  Produces<aod::DDbarRecoInfo> entryD0D0barRecoInfo;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<float> yCandMax{"yCandMax", 999., "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{efficiencyDmeson_v}, "Efficiency values for D0 meson"};

  // HfHelper hfHelper; //not needed for now

  Partition<soa::Join<aod::Alice3D0Meson, aod::Alice3D0Sel>> selectedCandidates = aod::a3D0meson::y > -yCandMax&& aod::a3D0meson::y<yCandMax && aod::a3D0meson::pt> ptCandMin && (aod::a3D0Selection::isSelD0 >= selectionFlagD0 || aod::a3D0Selection::isSelD0bar >= selectionFlagD0bar);
  Partition<soa::Join<aod::Alice3D0Meson, aod::Alice3D0Sel, aod::Alice3D0MCTruth>> selectedCandidatesMC = aod::a3D0meson::y > -yCandMax&& aod::a3D0meson::y<yCandMax && aod::a3D0meson::pt> ptCandMin && (aod::a3D0Selection::isSelD0 >= selectionFlagD0 || aod::a3D0Selection::isSelD0bar >= selectionFlagD0bar);

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassD0 for trigger normalisation (S*0.955), and hMass2DCorrelationPairs (in final task) for 2D-sideband-subtraction purposes
    {{"hPtCand", "D0,D0bar candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "D0,D0bar candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "D0,D0bar candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "D0,D0bar candidates;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEta", "D0,D0bar candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "D0,D0bar candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "D0,D0bar candidates;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hDDbarVsEtaCut", "D0,D0bar pairs vs #eta cut;#eta_{max};candidates #it{p}_{T} threshold (GeV/#it{c});entries", {HistType::kTH2F, {{static_cast<int>(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {static_cast<int>(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}},
     {"hPtCandMCRec", "D0,D0bar candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0MCRec", "D0,D0bar candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1MCRec", "D0,D0bar candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusMCRec", "D0,D0bar candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEtaMCRec", "D0,D0bar candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCRec", "D0,D0bar candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCRec", "D0,D0bar candidates - MC reco;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMass", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0", "D0,D0bar candidates;inv. mass D0 only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0bar", "D0,D0bar candidates;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0MCRecSig", "D0 signal candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0MCRecRefl", "D0 reflection candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0MCRecBkg", "D0 background candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecSig", "D0bar signal candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecRefl", "D0bar reflection candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecBkg", "D0bar background candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMass_NoEff", "D0,D0bar candidates (wo efficiency);inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0_NoEff", "D0,D0bar candidates (wo efficiency);inv. mass D0 only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0bar_NoEff", "D0,D0bar candidates (wo efficiency);inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0MCRecSig_NoEff", "D0 signal candidates - MC reco (wo efficiency);inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0MCRecRefl_NoEff", "D0 reflection candidates - MC reco (wo efficiency);inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0MCRecBkg_NoEff", "D0 background candidates - MC reco (wo efficiency);inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecSig_NoEff", "D0bar signal candidates - MC reco (wo efficiency);inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecRefl_NoEff", "D0bar reflection candidates - MC reco (wo efficiency);inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecBkg_NoEff", "D0bar background candidates - MC reco (wo efficiency);inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  /// D0-D0bar correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(aod::Collision const& collision,
                   soa::Join<aod::Alice3D0Meson, aod::Alice3D0Sel> const&)
  {
    auto selectedCandidatesGrouped = selectedCandidates->sliceByCached(aod::a3D0meson::collisionId, collision.globalIndex(), cache);

    for (const auto& candidate1 : selectedCandidatesGrouped) { // loop over reconstructed and selected D0 and D0bar (together, to fill mass plots first)
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
      }

      // fill invariant mass plots and generic info from all D0/D0bar candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), candidate1.m(), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassD0"), candidate1.m(), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMass_NoEff"), candidate1.m(), candidate1.pt());
        registry.fill(HIST("hMassD0_NoEff"), candidate1.m(), candidate1.pt());
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), candidate1.m(), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassD0bar"), candidate1.m(), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMass_NoEff"), candidate1.m(), candidate1.pt());
        registry.fill(HIST("hMassD0bar_NoEff"), candidate1.m(), candidate1.pt());
      }
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), candidate1.y());
      registry.fill(HIST("hSelectionStatus"), candidate1.isSelD0() + (candidate1.isSelD0bar() * 2));

      // D-Dbar correlation dedicated section

      // if the candidate is a D0, search for D0bar and evaluate correlations
      if (candidate1.isSelD0() < selectionFlagD0) {
        continue;
      }
      for (const auto& candidate2 : selectedCandidatesGrouped) {
        if (candidate2.isSelD0bar() < selectionFlagD0bar) { // keep only D0bar candidates passing the selection
          continue;
        }
        // excluding trigger self-correlations (possible in case of both mass hypotheses accepted)
        if (candidate1.mRowIndex == candidate2.mRowIndex) { // this by definition should never happen, since each candidate is either D0 or D0bar
          continue;
        }
        if ((candidate1.pt() - candidate2.pt()) < 1e-5 && (candidate1.eta() - candidate2.eta()) < 1e-5 && (candidate1.phi() - candidate2.phi()) < 1e-5) { // revised, temporary condition to avoid self-correlations (the best would be check the daughterIDs, but we don't store them at the moment)
          continue;
        }
        entryD0D0barPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                         candidate2.eta() - candidate1.eta(),
                         candidate1.pt(),
                         candidate2.pt());
        entryD0D0barRecoInfo(candidate1.m(), // mD0
                             candidate2.m(), // mD0bar
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
        // note: candidates selected as both D0 and D0bar are used, and considered in both situation (but not auto-correlated): reflections could play a relevant role.
        // another, more restrictive, option, could be to consider only candidates selected with a single option (D0 xor D0bar)

      } // end inner loop (Dbars)

    } // end outer loop
  }

  PROCESS_SWITCH(alice3correlatorddbar, processData, "Process data", true);

  /// D0-D0bar correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(aod::Collision const& collision,
                    soa::Join<aod::Alice3D0Meson, aod::Alice3D0Sel, aod::Alice3D0MCTruth> const&)
  {
    auto selectedCandidatesGroupedMC = selectedCandidatesMC->sliceByCached(aod::a3D0meson::collisionId, collision.globalIndex(), cache);

    // MC reco level
    bool flagD0Signal = false;
    bool flagD0Reflection = false;
    bool flagD0barSignal = false;
    bool flagD0barReflection = false;
    for (const auto& candidate1 : selectedCandidatesGroupedMC) {
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
      }

      if (candidate1.mcTruthInfo()) { // 1 or 2, i.e. true D0 or D0bar
        // fill per-candidate distributions from D0/D0bar true candidates
        registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
        registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
        registry.fill(HIST("hEtaMCRec"), candidate1.eta());
        registry.fill(HIST("hPhiMCRec"), candidate1.phi());
        registry.fill(HIST("hYMCRec"), candidate1.y());
        registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelD0() + (candidate1.isSelD0bar() * 2));
      }
      // fill invariant mass plots from D0/D0bar signal and background candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {                                                 // only reco as D0
        if (candidate1.mcTruthInfo() == 1) {                                                         // also matched as D0
          registry.fill(HIST("hMassD0MCRecSig"), candidate1.m(), candidate1.pt(), efficiencyWeight); // here m is univoque, since a given candidate passes the selection with only a single mass option
          registry.fill(HIST("hMassD0MCRecSig_NoEff"), candidate1.m(), candidate1.pt());
        } else if (candidate1.mcTruthInfo() == 2) {
          registry.fill(HIST("hMassD0MCRecRefl"), candidate1.m(), candidate1.pt(), efficiencyWeight);
          registry.fill(HIST("hMassD0MCRecRefl_NoEff"), candidate1.m(), candidate1.pt());
        } else {
          registry.fill(HIST("hMassD0MCRecBkg"), candidate1.m(), candidate1.pt(), efficiencyWeight);
          registry.fill(HIST("hMassD0MCRecBkg_NoEff"), candidate1.m(), candidate1.pt());
        }
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {                                              // only reco as D0bar
        if (candidate1.mcTruthInfo() == 2) {                                                            // also matched as D0bar
          registry.fill(HIST("hMassD0barMCRecSig"), candidate1.m(), candidate1.pt(), efficiencyWeight); // here m is univoque, since a given candidate passes the selection with only a single mass option
          registry.fill(HIST("hMassD0barMCRecSig_NoEff"), candidate1.m(), candidate1.pt());
        } else if (candidate1.mcTruthInfo() == 1) {
          registry.fill(HIST("hMassD0barMCRecRefl"), candidate1.m(), candidate1.pt(), efficiencyWeight);
          registry.fill(HIST("hMassD0barMCRecRefl_NoEff"), candidate1.m(), candidate1.pt());
        } else {
          registry.fill(HIST("hMassD0barMCRecBkg"), candidate1.m(), candidate1.pt(), efficiencyWeight);
          registry.fill(HIST("hMassD0barMCRecBkg_NoEff"), candidate1.m(), candidate1.pt());
        }
      }

      // D-Dbar correlation dedicated section
      // if the candidate is selected ad D0, search for D0bar and evaluate correlations
      if (candidate1.isSelD0() < selectionFlagD0) { // discard candidates not selected as D0 in outer loop
        continue;
      }

      flagD0Signal = candidate1.mcTruthInfo() == 1;     // flagD0Signal 'true' if candidate1 matched to D0 (particle)
      flagD0Reflection = candidate1.mcTruthInfo() == 2; // flagD0Reflection 'true' if candidate1, selected as D0 (particle), is matched to D0bar (antiparticle)

      for (const auto& candidate2 : selectedCandidatesGroupedMC) {
        if (candidate2.isSelD0bar() < selectionFlagD0bar) { // discard candidates not selected as D0bar in inner loop
          continue;
        }
        flagD0barSignal = candidate2.mcTruthInfo() == 2;     // flagD0barSignal 'true' if candidate2 matched to D0bar (antiparticle)
        flagD0barReflection = candidate2.mcTruthInfo() == 1; // flagD0barReflection 'true' if candidate2, selected as D0bar (antiparticle), is matched to D0 (particle)

        // Excluding trigger self-correlations (possible in case of both mass hypotheses of the same real particle, reconstructed as candidate1 for D0 and candidate2 for D0bar)
        if (candidate1.mRowIndex == candidate2.mRowIndex) { // this by definition should never happen, since each candidate is either D0 or D0bar
          continue;
        }
        if ((candidate1.pt() - candidate2.pt()) < 1e-5 && (candidate1.eta() - candidate2.eta()) < 1e-5 && (candidate1.phi() - candidate2.phi()) < 1e-5) { // revised, temporary condition to avoid self-correlations (the best would be check the daughterIDs, but we don't store them at the moment)
          continue;
        }
        // choice of options (D0/D0bar signal/bkg)
        int pairSignalStatus = 0; // 0 = bkg/bkg, 1 = bkg/ref, 2 = bkg/sig, 3 = ref/bkg, 4 = ref/ref, 5 = ref/sig, 6 = sig/bkg, 7 = sig/ref, 8 = sig/sig
        if (flagD0Signal) {
          pairSignalStatus += 6;
        }
        if (flagD0Reflection) {
          pairSignalStatus += 3;
        }
        if (flagD0barSignal) {
          pairSignalStatus += 2;
        }
        if (flagD0barReflection) {
          pairSignalStatus += 1;
        }
        entryD0D0barPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                         candidate2.eta() - candidate1.eta(),
                         candidate1.pt(),
                         candidate2.pt());
        entryD0D0barRecoInfo(candidate1.m(), // mD0
                             candidate2.m(), // mD0bar
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

  PROCESS_SWITCH(alice3correlatorddbar, processMcRec, "Process MC Reco mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3correlatorddbar>(cfgc)};
}
