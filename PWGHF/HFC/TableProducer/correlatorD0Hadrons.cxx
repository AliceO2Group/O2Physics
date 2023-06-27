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

/// \file correlatorD0Hadrons.cxx
/// \brief D0-Hadron correlator task - data-like, MC-reco and MC-kine analyses.
///
/// \author Samrangy Sadhu <samrangy.sadhu@cern.ch>, INFN Bari
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>, IIT Indore

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::aod::hf_correlation_d0_hadron;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;
using namespace o2::constants::math;

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PIHalf);
}

const int nPtBinsMassAndEfficiency = o2::analysis::hf_cuts_d0_to_pi_k::nBinsPt;
const double efficiencyDmesonDefault[nPtBinsMassAndEfficiency] = {};
auto vecEfficiencyDmeson = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + nPtBinsMassAndEfficiency};

// histogram binning definition
const int massAxisNBins = 120;
const double massAxisMin = 1.5848;
const double massAxisMax = 2.1848;
const int phiAxisNBins = 32;
const double phiAxisMin = 0.;
const double phiAxisMax = o2::constants::math::TwoPI;
const int yAxisNBins = 100;
const double yAxisMin = -5.;
const double yAxisMax = 5.;
const int ptDAxisNBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;
const double massD0 = RecoDecay::getMassPDG(pdg::Code::kD0);
const double softPiMass = 0.14543; // pion mass + Q-value of the D*->D0pi decay
auto massPi = RecoDecay::getMassPDG(kPiPlus);
auto massK = RecoDecay::getMassPDG(kKPlus);

struct HfCorrelatorD0Hadrons {
  SliceCache cache;
  Preslice<aod::HfCand2Prong> perCol = aod::hf_cand::collisionId;

  Produces<aod::DHadronPair> entryD0HadronPair;
  Produces<aod::DHadronRecoInfo> entryD0HadronRecoInfo;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<double> etaTrackMax{"etaTrackMax", 4., "max. eta of tracks"};
  Configurable<double> dcaXYTrackMax{"dcaXYTrackMax", 0.0025, "max. DCAxy of tracks"};
  Configurable<double> dcaZTrackMax{"dcaZTrackMax", 0.0025, "max. DCAz of tracks"};
  Configurable<double> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
  Configurable<double> ptTrackMin{"ptTrackMin", -1., "min. track pT"};
  Configurable<std::vector<double>> bins{"ptBinsForMassAndEfficiency", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for D0 meson"};
  Configurable<int> applyEfficiency{"efficiencyFlagD", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<double> ptSoftPionMax{"ptSoftPionMax", 3 * 800. * pow(10., -6.), "max. pT cut for soft pion identification"};

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> selectedD0candidatesMC = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassD0 for trigger normalisation (S*0.955), and hMass2DCorrelationPairs (in final task) for 2D-sideband-subtraction purposes
    {{"hPtCand", "D0,D0bar candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "D0,D0bar candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "D0,D0bar candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "D0,D0bar candidates;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEta", "D0,D0bar candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "D0,D0bar candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisNBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "D0,D0bar candidates;candidate #it{y};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hPtCandRec", "D0,D0bar candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0Rec", "D0,D0bar candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1Rec", "D0,D0bar candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusRec", "D0,D0bar candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEtaRec", "D0,D0bar candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hPhiRec", "D0,D0bar candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisNBins, phiAxisMin, phiAxisMax}}}},
     {"hYRec", "D0,D0bar candidates - MC reco;candidate #it{y};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hEvtCountGen", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandGen", "D0,D0bar particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaGen", "D0,D0bar particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hPhiGen", "D0,D0bar particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisNBins, phiAxisMin, phiAxisMax}}}},
     {"hYGen", "D0,D0bar candidates - MC gen;candidate #it{y};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hTrackCounter", "soft pion counter -  Data", {HistType::kTH1F, {{5, 0., 5.}}}},
     {"hTrackCounterRec", "soft pion counter - MC rec", {HistType::kTH1F, {{5, 0., 5.}}}},
     {"hTrackCounterGen", "soft pion counter - MC gen", {HistType::kTH1F, {{5, 0., 5.}}}}}};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hMass", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMass1D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
    registry.add("hMassD01D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
    registry.add("hMassD0bar1D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
    // mass histogram for D0 signal only
    registry.add("hMassD0RecSig", "D0 signal candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0 Reflection candidates only
    registry.add("hMassD0RecRef", "D0 reflection candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0 background candidates only
    registry.add("hMassD0RecBg", "D0 background candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0bar signal only
    registry.add("hMassD0barRecSig", "D0bar signal candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0bar Reflection candidates only
    registry.add("hMassD0barRecRef", "D0bar reflection candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0bar background candidates only
    registry.add("hMassD0barRecBg", "D0bar background candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCountD0TriggersGen", "D0 trigger particles - MC gen;;N of trigger D0", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // =============================================================================  Process starts for Data ==================================================================================
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// D0-h correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates)
  {
    // protection against empty tables to be sliced
    if (selectedD0Candidates.size() == 0) {
      return;
    }
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

    auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

    for (auto const& candidate1 : selectedD0CandidatesGrouped) {
      if (yCandMax >= 0. && std::abs(yD0(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }

      // ========================== Define parameters for soft pion removal ================================
      auto ePiK = RecoDecay::e(candidate1.pVectorProng0(), massPi) + RecoDecay::e(candidate1.pVectorProng1(), massK);
      auto eKPi = RecoDecay::e(candidate1.pVectorProng0(), massK) + RecoDecay::e(candidate1.pVectorProng1(), massPi);

      // ========================== trigger efficiency ================================
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }
      // ========================== Fill mass histo  ================================
      if (candidate1.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMass1D"), invMassD0ToPiK(candidate1), efficiencyWeight);
        registry.fill(HIST("hMassD01D"), invMassD0ToPiK(candidate1), efficiencyWeight);
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMass1D"), invMassD0barToKPi(candidate1), efficiencyWeight);
        registry.fill(HIST("hMassD0bar1D"), invMassD0barToKPi(candidate1), efficiencyWeight);
      }
      // ========================== Fill general histos ================================
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), yD0(candidate1));
      registry.fill(HIST("hSelectionStatus"), candidate1.isSelD0bar() + (candidate1.isSelD0() * 2));

      // ================================================================================= D-h correlation dedicated section =====================================================

      // ========================== track loop starts here ================================
      for (const auto& track : tracks) {
        registry.fill(HIST("hTrackCounter"), 1); // fill total no. of tracks
        // Remove D0 daughters by checking track indices
        if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex())) {
          continue;
        }
        if (std::abs(track.dcaXY()) >= 1. || std::abs(track.dcaZ()) >= 1.)
          continue; // Remove secondary tracks

        registry.fill(HIST("hTrackCounter"), 2); // fill no. of tracks before soft pion removal

        // ===== soft pion removal ===================================================
        double invMassDstar1 = 0., invMassDstar2 = 0.;
        bool isSoftpiD0 = false, isSoftpiD0bar = false;
        auto pSum2 = RecoDecay::p2(candidate1.px() + track.px(), candidate1.py() + track.py(), candidate1.pz() + track.pz());
        auto ePion = track.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (candidate1.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - invMassD0ToPiK(candidate1)) - softPiMass) < ptSoftPionMax) {
            isSoftpiD0 = true;
            continue;
          }
        }

        if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - invMassD0barToKPi(candidate1)) - softPiMass) < ptSoftPionMax) {
            isSoftpiD0bar = true;
            continue;
          }
        }
        registry.fill(HIST("hTrackCounter"), 3); // fill no. of tracks after soft pion removal

        int signalStatus = 0;
        if ((candidate1.isSelD0() >= selectionFlagD0) && (isSoftpiD0 == false)) {
          signalStatus += 1;
        }
        if ((candidate1.isSelD0bar() >= selectionFlagD0bar) && (isSoftpiD0bar == false)) {
          signalStatus += 2;
        }

        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                          track.eta() - candidate1.eta(),
                          candidate1.pt(),
                          track.pt());
        entryD0HadronRecoInfo(invMassD0ToPiK(candidate1), invMassD0barToKPi(candidate1), signalStatus);

      } // end inner loop (tracks)

    } // end outer loop
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processData, "Process data", false);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // =============================================================================  Process starts for MCRec ==================================================================================
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void processMcRec(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const& candidates)
  {
    // protection against empty tables to be sliced
    if (selectedD0candidatesMC.size() == 0) {
      return;
    }
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

    auto selectedD0CandidatesGroupedMC = selectedD0candidatesMC->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    // MC reco level
    bool flagD0 = false;
    bool flagD0bar = false;

    for (auto const& candidate1 : selectedD0CandidatesGroupedMC) {
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yD0(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }

      if (std::abs(candidate1.flagMcMatchRec()) == 1 << DecayType::D0ToPiK) {
        // fill per-candidate distributions from D0/D0bar true candidates
        registry.fill(HIST("hPtCandRec"), candidate1.pt());
        registry.fill(HIST("hPtProng0Rec"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1Rec"), candidate1.ptProng1());
        registry.fill(HIST("hEtaRec"), candidate1.eta());
        registry.fill(HIST("hPhiRec"), candidate1.phi());
        registry.fill(HIST("hYRec"), yD0(candidate1));
        registry.fill(HIST("hSelectionStatusRec"), candidate1.isSelD0bar() + (candidate1.isSelD0() * 2));
      }
      // fill invariant mass plots from D0/D0bar signal and background candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {                  // only reco as D0
        if (candidate1.flagMcMatchRec() == 1 << DecayType::D0ToPiK) { // also matched as D0
          registry.fill(HIST("hMassD0RecSig"), invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        } else if (candidate1.flagMcMatchRec() == -(1 << DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassD0RecRef"), invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0RecBg"), invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {               // only reco as D0bar
        if (candidate1.flagMcMatchRec() == -(1 << DecayType::D0ToPiK)) { // also matched as D0bar
          registry.fill(HIST("hMassD0barRecSig"), invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else if (candidate1.flagMcMatchRec() == 1 << DecayType::D0ToPiK) {
          registry.fill(HIST("hMassD0barRecRef"), invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0barRecBg"), invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }

      // ========================== Define parameters for soft pion removal ================================
      auto ePiK = RecoDecay::e(candidate1.pVectorProng0(), massPi) + RecoDecay::e(candidate1.pVectorProng1(), massK);
      auto eKPi = RecoDecay::e(candidate1.pVectorProng0(), massK) + RecoDecay::e(candidate1.pVectorProng1(), massPi);

      // D0-h correlation dedicated section

      flagD0 = candidate1.flagMcMatchRec() == (1 << DecayType::D0ToPiK);     // flagD0Signal 'true' if candidate1 matched to D0 (particle)
      flagD0bar = candidate1.flagMcMatchRec() == -(1 << DecayType::D0ToPiK); // flagD0Reflection 'true' if candidate1, selected as D0 (particle), is matched to D0bar (antiparticle)

      // ========== track loop starts here ========================

      for (const auto& track : tracks) {
        registry.fill(HIST("hTrackCounterRec"), 1); // fill total no. of tracks
        if (std::abs(track.eta()) > etaTrackMax) {
          continue;
        }
        if (track.pt() < ptTrackMin) {
          continue;
        }
        // Removing D0 daughters by checking track indices
        if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex())) {
          continue;
        }
        if (std::abs(track.dcaXY()) >= 1. || std::abs(track.dcaZ()) >= 1.) {
          continue; // Remove secondary tracks
        }
        registry.fill(HIST("hTrackCounterRec"), 2); // fill no. of tracks before soft pion removal

        // ===== soft pion removal ===================================================
        double invMassDstar1 = 0, invMassDstar2 = 0;
        bool isSoftpiD0 = false, isSoftpiD0bar = false;
        auto pSum2 = RecoDecay::p2(candidate1.px() + track.px(), candidate1.py() + track.py(), candidate1.pz() + track.pz());
        auto ePion = track.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (candidate1.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - invMassD0ToPiK(candidate1)) - softPiMass) < ptSoftPionMax) {
            isSoftpiD0 = true;
            continue;
          }
        }

        if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - invMassD0barToKPi(candidate1)) - softPiMass) < ptSoftPionMax) {
            isSoftpiD0bar = true;
            continue;
          }
        }

        registry.fill(HIST("hTrackCounterRec"), 3); // fill no. of tracks after soft pion removal

        int signalStatus = 0;
        if ((flagD0 == true) && (candidate1.isSelD0() >= selectionFlagD0) && (isSoftpiD0 == false)) {
          signalStatus += 1;
        } // signal case D0
        if ((flagD0bar == true) && (candidate1.isSelD0() >= selectionFlagD0) && (isSoftpiD0 == false)) {
          signalStatus += 2;
        } // reflection case D0
        if ((flagD0 == false) && (flagD0bar == false) && (candidate1.isSelD0() >= selectionFlagD0) && (isSoftpiD0 == false)) {
          signalStatus += 4;
        } // background case D0

        if ((flagD0bar == true) && (candidate1.isSelD0bar() >= selectionFlagD0bar) && (isSoftpiD0bar == false)) {
          signalStatus += 8;
        } // signal case D0bar
        if ((flagD0 == true) && (candidate1.isSelD0bar() >= selectionFlagD0bar) && (isSoftpiD0bar == false)) {
          signalStatus += 16;
        } // reflection case D0bar
        if ((flagD0 == false) && (flagD0bar == false) && (candidate1.isSelD0bar() >= selectionFlagD0bar) && (isSoftpiD0bar == false)) {
          signalStatus += 32;
        } // background case D0bar

        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                          track.eta() - candidate1.eta(),
                          candidate1.pt(),
                          track.pt());
        entryD0HadronRecoInfo(invMassD0ToPiK(candidate1), invMassD0barToKPi(candidate1), signalStatus);
      } // end inner loop (Tracks)

    } // end of outer loop (D0)
  }

  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcRec, "Process MC Reco mode", true);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // =============================================================================  Process starts for MCGen ==================================================================================
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void processMcGen(aod::McCollision const& mccollision, soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& particlesMC)
  {

    registry.fill(HIST("hEvtCountGen"), 0);
    // MC gen level
    for (auto const& particle1 : particlesMC) {
      // check if the particle is D0 or D0bar (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != pdg::Code::kD0) {
        continue;
      }
      double yD = RecoDecay::y(array{particle1.px(), particle1.py(), particle1.pz()}, RecoDecay::getMassPDG(particle1.pdgCode()));
      if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandGen"), particle1.pt());
      registry.fill(HIST("hEtaGen"), particle1.eta());
      registry.fill(HIST("hPhiGen"), particle1.phi());
      registry.fill(HIST("hYGen"), yD);

      // D-h correlation dedicated section

      if (std::abs(particle1.pdgCode()) != pdg::Code::kD0) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hCountD0TriggersGen"), 0, particle1.pt()); // to count trigger D0 (for normalisation)

      for (auto const& particle2 : particlesMC) {
        registry.fill(HIST("hTrackCounterGen"), 1); // total no. of tracks
        if (std::abs(particle2.eta()) > etaTrackMax) {
          continue;
        }
        if (particle2.pt() < ptTrackMin) {
          continue;
        }
        if ((std::abs(particle2.pdgCode()) != kElectron) && (std::abs(particle2.pdgCode()) != kMuonMinus) && (std::abs(particle2.pdgCode()) != kPiPlus) && (std::abs(particle2.pdgCode()) != kKPlus) && (std::abs(particle2.pdgCode()) != kProton)) {
          continue;
        }

        // ==============================soft pion removal================================
        registry.fill(HIST("hTrackCounterGen"), 2); // fill before soft pi removal
        // method used: indexMother = -1 by default if the mother doesn't match with given PID of the mother. We find mother of pion if it is D* and mother of D0 if it is D*. If they are both positive and they both match each other, then it is detected as a soft pion

        auto indexMotherPi = RecoDecay::getMother(particlesMC, particle2, pdg::Code::kDStar, true, nullptr, 1); // last arguement 1 is written to consider immediate decay mother only
        auto indexMotherD0 = RecoDecay::getMother(particlesMC, particle1, pdg::Code::kDStar, true, nullptr, 1);
        if (std::abs(particle2.pdgCode()) == kPiPlus && indexMotherPi >= 0 && indexMotherD0 >= 0 && indexMotherPi == indexMotherD0)
          continue;

        registry.fill(HIST("hTrackCounterGen"), 3); // fill after soft pion removal
        entryD0HadronPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                          particle2.eta() - particle1.eta(),
                          particle1.pt(),
                          particle2.pt());
        entryD0HadronRecoInfo(massD0,
                              massD0,
                              0); // dummy info
      }                           // end inner loop

    } // end outer loop
  }

  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorD0Hadrons>(cfgc)};
}
