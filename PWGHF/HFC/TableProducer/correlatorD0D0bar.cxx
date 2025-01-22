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

/// \file correlatorD0D0bar.cxx
/// \brief D0-D0bar correlator task - data-like, MC-reco and MC-kine analyses. For ULS and LS pairs
///
/// \author Fabio Colamaria <fabio.colamaria@ba.infn.it>, INFN Bari

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

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

const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_d0_to_pi_k::nBinsPt;
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

using MCParticlesPlus = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

struct HfCorrelatorD0D0bar {
  SliceCache cache;
  Preslice<aod::HfCand2Prong> perCol = aod::hf_cand::collisionId;
  Produces<aod::DDbarPair> entryD0D0barPair;
  Produces<aod::DDbarRecoInfo> entryD0D0barRecoInfo;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<double> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
  Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{efficiencyDmeson_v}, "Efficiency values for D0 meson"};

  HfHelper hfHelper;

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> selectedD0candidatesMC = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

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
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hDDbarVsEtaCut", "D0,D0bar pairs vs #eta cut;#eta_{max};candidates #it{p}_{T} threshold (GeV/#it{c});entries", {HistType::kTH2F, {{static_cast<int>(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {static_cast<int>(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}},
     {"hPtCandMCRec", "D0,D0bar candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0MCRec", "D0,D0bar candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1MCRec", "D0,D0bar candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusMCRec", "D0,D0bar candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEtaMCRec", "D0,D0bar candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCRec", "D0,D0bar candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCRec", "D0,D0bar candidates - MC reco;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMCEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMCGen", "D0,D0bar particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaMCGen", "D0,D0bar particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCGen", "D0,D0bar particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCGen", "D0,D0bar candidates - MC gen;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hCountD0D0barPerEvent", "D0,D0bar particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hDDbarVsDaughterEtaCut", "D0,D0bar pairs vs #eta cut on D daughters;#eta_{max};candidates #it{p}_{T} threshold (GeV/#it{c});entries", {HistType::kTH2F, {{static_cast<int>(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {static_cast<int>(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}},
     {"hCountCCbarPerEvent", "c,cbar particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountCCbarPerEventBeforeEtaCut", "c,cbar particles - MC gen;Number per event pre #eta cut;entries", {HistType::kTH1F, {{20, 0., 20.}}}}}};

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
    registry.add("hCountD0triggersMCGen", "D0 trigger particles - MC gen;;N of trigger D0", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCountCtriggersMCGen", "c trigger particles - MC gen;;N of trigger c quark", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  /// D0-D0bar correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(aod::Collision const& collision,
                   aod::TracksWDca const& tracks,
                   soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&)
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

    auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

    for (const auto& candidate1 : selectedD0CandidatesGrouped) {
      if (yCandMax >= 0. && std::abs(hfHelper.yD0(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
      }

      // fill invariant mass plots and generic info from all D0/D0bar candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassD0"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassD0bar"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
      }
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), hfHelper.yD0(candidate1));
      registry.fill(HIST("hSelectionStatus"), candidate1.isSelD0bar() + (candidate1.isSelD0() * 2));

      // D-Dbar correlation dedicated section
      // if the candidate is a D0, search for D0bar and evaluate correlations
      if (candidate1.isSelD0() < selectionFlagD0) {
        continue;
      }
      for (const auto& candidate2 : selectedD0CandidatesGrouped) {
        if (!(candidate2.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) { // check decay channel flag for candidate2
          continue;
        }
        if (candidate2.isSelD0bar() < selectionFlagD0bar) { // keep only D0bar candidates passing the selection
          continue;
        }
        // kinematic selection on D0bar candidates
        if (yCandMax >= 0. && std::abs(hfHelper.yD0(candidate2)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate2.pt() < ptCandMin) {
          continue;
        }
        // excluding trigger self-correlations (possible in case of both mass hypotheses accepted)
        if (candidate1.mRowIndex == candidate2.mRowIndex) {
          continue;
        }
        entryD0D0barPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                         candidate2.eta() - candidate1.eta(),
                         candidate1.pt(),
                         candidate2.pt());
        entryD0D0barRecoInfo(hfHelper.invMassD0ToPiK(candidate1),
                             hfHelper.invMassD0barToKPi(candidate2),
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

  PROCESS_SWITCH(HfCorrelatorD0D0bar, processData, "Process data", false);

  /// D0-D0bar correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(aod::Collision const& collision,
                    aod::TracksWDca const& tracks,
                    soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const&)
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

    auto selectedD0CandidatesGroupedMC = selectedD0candidatesMC->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

    // MC reco level
    bool flagD0Signal = false;
    bool flagD0Reflection = false;
    bool flagD0barSignal = false;
    bool flagD0barReflection = false;
    for (const auto& candidate1 : selectedD0CandidatesGroupedMC) {
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(hfHelper.yD0(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
      }

      if (std::abs(candidate1.flagMcMatchRec()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {
        // fill per-candidate distributions from D0/D0bar true candidates
        registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
        registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
        registry.fill(HIST("hEtaMCRec"), candidate1.eta());
        registry.fill(HIST("hPhiMCRec"), candidate1.phi());
        registry.fill(HIST("hYMCRec"), hfHelper.yD0(candidate1));
        registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelD0bar() + (candidate1.isSelD0() * 2));
      }
      // fill invariant mass plots from D0/D0bar signal and background candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {                  // only reco as D0
        if (candidate1.flagMcMatchRec() == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) { // also matched as D0
          registry.fill(HIST("hMassD0MCRecSig"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        } else if (candidate1.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassD0MCRecRefl"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0MCRecBkg"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {               // only reco as D0bar
        if (candidate1.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) { // also matched as D0bar
          registry.fill(HIST("hMassD0barMCRecSig"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else if (candidate1.flagMcMatchRec() == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {
          registry.fill(HIST("hMassD0barMCRecRefl"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0barMCRecBkg"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }

      // D-Dbar correlation dedicated section
      // if the candidate is selected ad D0, search for D0bar and evaluate correlations
      if (candidate1.isSelD0() < selectionFlagD0) { // discard candidates not selected as D0 in outer loop
        continue;
      }
      flagD0Signal = candidate1.flagMcMatchRec() == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;        // flagD0Signal 'true' if candidate1 matched to D0 (particle)
      flagD0Reflection = candidate1.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK); // flagD0Reflection 'true' if candidate1, selected as D0 (particle), is matched to D0bar (antiparticle)
      for (const auto& candidate2 : selectedD0CandidatesGroupedMC) {
        if (!(candidate2.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) { // check decay channel flag for candidate2
          continue;
        }
        if (candidate2.isSelD0bar() < selectionFlagD0bar) { // discard candidates not selected as D0bar in inner loop
          continue;
        }
        flagD0barSignal = candidate2.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK);  // flagD0barSignal 'true' if candidate2 matched to D0bar (antiparticle)
        flagD0barReflection = candidate2.flagMcMatchRec() == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK; // flagD0barReflection 'true' if candidate2, selected as D0bar (antiparticle), is matched to D0 (particle)
        if (yCandMax >= 0. && std::abs(hfHelper.yD0(candidate2)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate2.pt() < ptCandMin) {
          continue;
        }
        // Excluding trigger self-correlations (possible in case of both mass hypotheses accepted)
        if (candidate1.mRowIndex == candidate2.mRowIndex) {
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
        entryD0D0barRecoInfo(hfHelper.invMassD0ToPiK(candidate1),
                             hfHelper.invMassD0barToKPi(candidate2),
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

  PROCESS_SWITCH(HfCorrelatorD0D0bar, processMcRec, "Process MC Reco mode", true);

  /// D0-D0bar correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::McCollision const&,
                    MCParticlesPlus const& mcParticles)
  {
    int counterD0D0bar = 0;
    registry.fill(HIST("hMCEvtCount"), 0);
    // MC gen level
    for (const auto& particle1 : mcParticles) {
      // check if the particle is D0 or D0bar (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != Pdg::kD0) {
        continue;
      }
      double yD = RecoDecay::y(particle1.pVector(), MassD0);
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
      counterD0D0bar++;

      // D-Dbar correlation dedicated section
      // if it's a D0 particle, search for D0bar and evaluate correlations
      if (particle1.pdgCode() != Pdg::kD0) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hCountD0triggersMCGen"), 0, particle1.pt()); // to count trigger D0 (for normalisation)
      for (const auto& particle2 : mcParticles) {
        if (particle2.pdgCode() != Pdg::kD0Bar) { // check that inner particle is D0bar
          continue;
        }
        if (yCandMax >= 0. && std::abs(RecoDecay::y(particle2.pVector(), MassD0Bar)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle2.pt() < ptCandMin) {
          continue;
        }
        entryD0D0barPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                         particle2.eta() - particle1.eta(),
                         particle1.pt(),
                         particle2.pt());
        entryD0D0barRecoInfo(1.864,
                             1.864,
                             8); // dummy information
        double etaCut = 0.;
        double ptCut = 0.;

        // fill pairs vs etaCut plot
        bool rightDecayChannels = false;
        if ((std::abs(particle1.flagMcMatchGen()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) && (std::abs(particle2.flagMcMatchGen()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          rightDecayChannels = true;
        }
        do {
          ptCut = 0.;
          etaCut += incrementEtaCut;
          do {                                                                                                                                  // fill pairs vs etaCut plot
            if (std::abs(particle1.eta()) < etaCut && std::abs(particle2.eta()) < etaCut && particle1.pt() > ptCut && particle2.pt() > ptCut) { // fill with D and Dbar acceptance checks
              registry.fill(HIST("hDDbarVsEtaCut"), etaCut - epsilon, ptCut + epsilon);
            }
            if (rightDecayChannels) { // fill with D and Dbar daughter particls acceptance checks
              bool candidate1DauInAcc = true;
              bool candidate2DauInAcc = true;
              for (const auto& dau : particle1.daughters_as<MCParticlesPlus>()) {
                if (std::abs(dau.eta()) > etaCut) {
                  candidate1DauInAcc = false;
                  break;
                }
              }
              for (const auto& dau : particle2.daughters_as<MCParticlesPlus>()) {
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
    registry.fill(HIST("hCountD0D0barPerEvent"), counterD0D0bar);
  }

  PROCESS_SWITCH(HfCorrelatorD0D0bar, processMcGen, "Process MC Gen mode", false);

  /// c-cbar correlator table builder - for MC gen-level analysis
  void processCCbar(aod::McCollision const&,
                    MCParticlesPlus const& mcParticles)
  {
    registry.fill(HIST("hMCEvtCount"), 0);
    int counterCCbar = 0, counterCCbarBeforeEtasel = 0;

    // loop over particles at MC gen level
    for (const auto& particle1 : mcParticles) {
      if (std::abs(particle1.pdgCode()) != PDG_t::kCharm) { // search c or cbar particles
        continue;
      }
      int partMothPDG = particle1.mothers_as<MCParticlesPlus>().front().pdgCode();
      // check whether mothers of quark c/cbar are still '4'/'-4' particles - in that case the c/cbar quark comes from its own fragmentation, skip it
      if (partMothPDG == particle1.pdgCode()) {
        continue;
      }
      counterCCbarBeforeEtasel++; // count c or cbar (before kinematic selection)
      double yC = RecoDecay::y(particle1.pVector(), MassCharm);
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
        if (particle2.pdgCode() != PDG_t::kCharmBar) { // check that inner particle is a cbar
          continue;
        }
        if (yCandMax >= 0. && std::abs(RecoDecay::y(particle2.pVector(), MassCharmBar)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle2.pt() < ptCandMin) {
          continue;
        }
        // check whether mothers of quark cbar (from associated loop) are still '-4' particles - in that case the cbar quark comes from its own fragmentation, skip it
        if (particle2.mothers_as<MCParticlesPlus>().front().pdgCode() == PDG_t::kCharmBar) {
          continue;
        }
        entryD0D0barPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                         particle2.eta() - particle1.eta(),
                         particle1.pt(),
                         particle2.pt());
        entryD0D0barRecoInfo(1.864,
                             1.864,
                             8); // dummy information
      }                          // end inner loop
    }                            // end outer loop
    registry.fill(HIST("hCountCCbarPerEvent"), counterCCbar);
    registry.fill(HIST("hCountCCbarPerEventBeforeEtaCut"), counterCCbarBeforeEtasel);
  }

  PROCESS_SWITCH(HfCorrelatorD0D0bar, processCCbar, "Process c-cbar pairs", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorD0D0bar>(cfgc)};
}
