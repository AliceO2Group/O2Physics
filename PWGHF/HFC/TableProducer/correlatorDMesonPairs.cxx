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

/// \file correlatorDMesonPairs.cxx
/// \brief D0(bar) and DPlus(Minus) correlator task - data-like, MC-reco and MC-kine analyses.
///
/// \author Fabio Colamaria <fabio.colamaria@ba.infn.it>, INFN Bari
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DMesonPairsTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PIHalf);
}

namespace
{
enum CandidateTypeSel {
  SelectedD = 0, // This particle is selected as a D
  SelectedDbar,  // This particle is selected as a Dbar
  TrueD,         // This particle is a true D
  TrueDbar       // This particle is a true Dbar
};
} // namespace

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

struct HfCorrelatorDMesonPairs {
  SliceCache cache;
  Preslice<aod::HfCand2Prong> perCol2Prong = aod::hf_cand::collisionId;
  Preslice<aod::HfCand3Prong> perCol3Prong = aod::hf_cand::collisionId;

  Produces<aod::D0Pair> entryD0Pair;
  Produces<aod::D0PairRecoInfo> entryD0PairRecoInfo;
  Produces<aod::DplusPair> entryDplusPair;
  Produces<aod::DplusPairRecoInfo> entryDplusPairRecoInfo;

  Configurable<bool> verbose{"verbose", false, "Enable debug mode"};
  // Selection Flags
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagDplus{"selectionFlagDPlus", 1, "Selection Flag for DPlus"};
  // Kinematics
  Configurable<float> etaCandMax{"etaCandMax", 4.0, "max. eta accepted"};
  Configurable<float> etaCandMin{"etaCandMin", -4.0, "min. eta accepted"};
  Configurable<float> dcaXYMax{"dcaXYMax", 0.0025, "max. dcaXY accepted"};
  Configurable<float> dcaZMax{"dcaZMax", 0.0025, "max. dcaZ accepted"};

  Configurable<float> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
  Configurable<float> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<float> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots"};

  HfHelper hfHelper;

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> selectedD0CandidatesMc = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>> selectedDPlusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>> selectedDPlusCandidatesMc = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  // HistoTypes
  HistogramConfigSpec hTH1Pt{HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}};
  HistogramConfigSpec hTH1Y{HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}};
  HistogramConfigSpec hTH1Phi{HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}};
  HistogramConfigSpec hTH1Matched{HistType::kTH1F, {{3, -1.5, 1.5}}};
  HistogramConfigSpec hTH1Origin{HistType::kTH1F, {{3, -0.5, 2.5}}};
  HistogramConfigSpec hTH1Multiplicity{HistType::kTH1F, {{10000, 0., 10000.}}};

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassD0 for trigger normalisation (S*0.955), and hMass2DCorrelationPairs (in final task) for 2D-sideband-subtraction purposes
    {{"hPtCand", "D Meson pair candidates;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng0", "D Meson pair candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng1", "D Meson pair candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEta", "D Meson pair candidates;candidate #it{#eta};entries", hTH1Y},
     {"hPhi", "D Meson pair candidates;candidate #it{#varphi};entries", hTH1Phi},
     {"hY", "D Meson pair candidates;candidate #it{y};entries", hTH1Y},
     // Mc Reco
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", hTH1Multiplicity},
     {"hMultiplicity", "multiplicity;multiplicity;entries", hTH1Multiplicity},
     {"hPtCandMcRec", "D Meson pair candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng0McRec", "D Meson pair candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng1McRec", "D Meson pair candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEtaMcRec", "D Meson pair candidates - MC reco;candidate #it{#eta};entries", hTH1Y},
     {"hPhiMcRec", "D Meson pair candidates - MC reco;candidate #it{#varphi};entries", hTH1Phi},
     {"hYMcRec", "D Meson pair candidates - MC reco;candidate #it{y};entries", hTH1Y},
     {"hMatchedMcRec", "D Meson pair candidates - MC reco;MC Matched;entries", hTH1Matched},
     {"hOriginMcRec", "D Meson pair candidates - MC reco;prompt vs. non-prompt;entries", hTH1Origin},
     // Mc Gen
     {"hMcEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMcGen", "D Meson pair particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEtaMcGen", "D Meson pair particles - MC gen;particle #it{#eta};entries", hTH1Y},
     {"hPhiMcGen", "D Meson pair particles - MC gen;particle #it{#varphi};entries", hTH1Phi},
     {"hYMcGen", "D Meson pair candidates - MC gen;candidate #it{y};entries", hTH1Y},
     {"hMatchedMcGen", "D Meson pair candidates - MC gen;MC Matched;entries", hTH1Matched},
     {"hOriginMcGen", "D Meson pair candidates - MC gen;prompt vs. non-prompt;entries", hTH1Origin}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    constexpr int kNBinsSelStatus = 6;
    std::string labels[kNBinsSelStatus];
    labels[0] = "total # of Selected pairs";
    labels[1] = "# of Selected D + Dbar";
    labels[2] = "# of Selected D";
    labels[3] = "# of Selected Dbar";
    labels[4] = "# of True D";
    labels[5] = "# of True Dbar";
    AxisSpec axisSelStatus = {kNBinsSelStatus, 0.5, kNBinsSelStatus + 0.5, ""};
    registry.add("hSelectionStatus", "D Meson candidates;selection status;entries", HistType::kTH1F, {axisSelStatus});
    for (int iBin = 0; iBin < kNBinsSelStatus; iBin++) {
      registry.get<TH1>(HIST("hSelectionStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    registry.add("hMass", "D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  // Performs an analysis on multiplicity and fills histograms
  template <typename T, typename U>
  void analyseMultiplicity(const T& collision, const U& tracks)
  {
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < etaCandMin || track.eta() > etaCandMax) {
          continue;
        }
        if (std::abs(track.dcaXY()) > dcaXYMax || std::abs(track.dcaZ()) > dcaZMax) {
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
  }

  // Returns false if the candidate does not pass cuts on decay type, y max, and pt min. Used for data and MC reco.
  template <typename T, bool isD0>
  bool kinematicCuts(const T& candidate)
  {
    // check decay channel flag for candidate
    if constexpr (isD0) {
      if (!(TESTBIT(candidate.hfflag(), o2::aod::hf_cand_2prong::DecayType::D0ToPiK))) {
        return false;
      }
      if (yCandMax >= 0. && std::abs(candidate.y(o2::analysis::pdg::MassD0)) > yCandMax) {
        return false;
      }
    } else {
      if (!(TESTBIT(candidate.hfflag(), o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi))) {
        return false;
      }
      if (yCandMax >= 0. && std::abs(candidate.y(o2::analysis::pdg::MassDPlus)) > yCandMax) {
        return false;
      }
    }
    if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
      return false;
    }
    return true;
  }

  // Returns false if the candidate does not pass cuts on pdgCode, y max, and pt min. Used for MC gen.
  template <typename T>
  bool kinematicCutsGen(const T& particle)
  {
    // check if the particle is D or Dbar (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
    if (std::abs(particle.pdgCode()) != pdg::Code::kDPlus && std::abs(particle.pdgCode()) != pdg::Code::kD0) {
      return false;
    }
    if (yCandMax >= 0. && std::abs(particle.y()) > yCandMax) {
      return false;
    }
    if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
      return false;
    }
    return true;
  }

  // Fills histograms with basic kinematic info.
  template <typename T>
  void fillInfoHists(const T& candidate, bool const& isReco, bool const& isD0)
  {
    if (isReco) {
      registry.fill(HIST("hPtCandMcRec"), candidate.pt());
      registry.fill(HIST("hPtProng0McRec"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1McRec"), candidate.ptProng1());
      registry.fill(HIST("hEtaMcRec"), candidate.eta());
      registry.fill(HIST("hPhiMcRec"), candidate.phi());
      if (isD0) {
        registry.fill(HIST("hYMcRec"), candidate.y(o2::analysis::pdg::MassD0));
      } else {
        registry.fill(HIST("hYMcRec"), candidate.y(o2::analysis::pdg::MassDPlus));
      }
    } else {
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hEta"), candidate.eta());
      registry.fill(HIST("hPhi"), candidate.phi());
      if (isD0) {
        registry.fill(HIST("hY"), candidate.y(o2::analysis::pdg::MassD0));
      } else {
        registry.fill(HIST("hY"), candidate.y(o2::analysis::pdg::MassDPlus));
      }
    }
  }

  // Sets bits to select candidate type for D0.
  // SelectedD and SelectedDbar bits look at whether the candidate passed the selection flags, and TrueD and TrueDbar check if it is matched with Mc.
  template <typename T, bool isRecoMc>
  uint8_t assignCandidateTypeD0(const T& candidate)
  {
    uint8_t candidateType(0);
    if (candidate.isSelD0() >= selectionFlagD0) {
      SETBIT(candidateType, SelectedD);
      registry.fill(HIST("hSelectionStatus"), 2);
      registry.fill(HIST("hSelectionStatus"), 3);
    }
    if (candidate.isSelD0bar() >= selectionFlagD0bar) {
      SETBIT(candidateType, SelectedDbar);
      registry.fill(HIST("hSelectionStatus"), 2);
      registry.fill(HIST("hSelectionStatus"), 4);
    }
    if constexpr (isRecoMc) {
      if (candidate.flagMcMatchRec() == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) { // matched as D0
        SETBIT(candidateType, TrueD);
        registry.fill(HIST("hSelectionStatus"), 5);
      }
      if (candidate.flagMcMatchRec() == -(1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) { // matched as D0bar
        SETBIT(candidateType, TrueDbar);
        registry.fill(HIST("hSelectionStatus"), 6);
      }
    }
    return candidateType;
  }

  // Sets bits to select candidate type for D plus
  template <typename T, bool isRecoMc>
  uint8_t assignCandidateTypeDPlus(const T& candidate, int const& particleSign)
  {
    uint8_t candidateType(0);
    if (particleSign == 1) {
      SETBIT(candidateType, SelectedD);
      registry.fill(HIST("hSelectionStatus"), 2);
      registry.fill(HIST("hSelectionStatus"), 3);
    } else {
      SETBIT(candidateType, SelectedDbar);
      registry.fill(HIST("hSelectionStatus"), 2);
      registry.fill(HIST("hSelectionStatus"), 4);
    }
    if constexpr (isRecoMc) {
      if (std::abs(candidate.flagMcMatchRec()) == 1 << o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi) { // matched as DPlus
        SETBIT(candidateType, TrueD);
        registry.fill(HIST("hSelectionStatus"), 5);
      } else { // matched as D0bar
        SETBIT(candidateType, TrueDbar);
        registry.fill(HIST("hSelectionStatus"), 6);
      }
    }
    return candidateType;
  }

  // Sets bits to select candidate type at generator level
  template <typename T>
  uint8_t assignCandidateTypeGen(const T& candidate)
  {
    uint8_t candidateType(0);
    if (candidate.pdgCode() == pdg::Code::kDPlus || candidate.pdgCode() == pdg::Code::kD0) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
      SETBIT(candidateType, SelectedD);
      SETBIT(candidateType, TrueD);
    } else if (candidate.pdgCode() == -pdg::Code::kDPlus || candidate.pdgCode() == -pdg::Code::kD0) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
      SETBIT(candidateType, SelectedDbar);
      SETBIT(candidateType, TrueDbar);
    }
    return candidateType;
  }

  // Common code to analyse D0's and D+'s at Gen level.
  template <typename T>
  void analyseMcGen(const T& mcParticles)
  {
    registry.fill(HIST("hMcEvtCount"), 0);
    for (const auto& particle1 : mcParticles) {
      // check if the particle is D0, D0bar, DPlus or DMinus (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      auto pdgCode = std::abs(particle1.pdgCode());
      if (pdgCode != pdg::Code::kD0 && pdgCode != pdg::Code::kDPlus) {
        continue;
      }
      auto massD = pdgCode == pdg::Code::kD0 ? o2::analysis::pdg::MassD0 : o2::analysis::pdg::MassDPlus;
      double yD = RecoDecay::y(std::array{particle1.px(), particle1.py(), particle1.pz()}, massD);
      if (!kinematicCutsGen(particle1)) {
        continue;
      }

      registry.fill(HIST("hPtCandMcGen"), particle1.pt());
      registry.fill(HIST("hEtaMcGen"), particle1.eta());
      registry.fill(HIST("hPhiMcGen"), particle1.phi());
      registry.fill(HIST("hYMcGen"), yD);

      auto candidateType1 = assignCandidateTypeGen(particle1); // Candidate sign attribution

      // check if it's prompt or non-prompt
      int8_t originGen1 = particle1.originMcGen();
      registry.fill(HIST("hOriginMcGen"), originGen1);
      // check if it's MC matched
      int8_t matchedGen1 = particle1.flagMcMatchGen();
      registry.fill(HIST("hMatchedMcGen"), matchedGen1);

      for (const auto& particle2 : mcParticles) {
        // Candidate sign attribution.
        auto candidateType2 = assignCandidateTypeGen(particle2);
        if (!kinematicCutsGen(particle2)) {
          continue;
        }

        // check if it's prompt or non-prompt
        int8_t originGen2 = particle2.originMcGen();
        // check if it's MC matched
        int8_t matchedGen2 = particle2.flagMcMatchGen();

        // If both particles are D0's, fill D0Pair table
        if (std::abs(particle1.pdgCode()) == pdg::Code::kD0 && std::abs(particle2.pdgCode()) == pdg::Code::kD0) {
          entryD0Pair(getDeltaPhi(particle2.phi(), particle1.phi()),
                      particle2.eta() - particle1.eta(),
                      particle1.pt(),
                      particle2.pt(),
                      particle1.y(),
                      particle2.y(),
                      o2::analysis::pdg::MassD0,
                      o2::analysis::pdg::MassD0,
                      candidateType1,
                      candidateType2,
                      2);
          entryD0PairRecoInfo(originGen1,
                              originGen2,
                              matchedGen1,
                              matchedGen2);
          // If both particles are DPlus, fill DplusPair table
        } else if (std::abs(particle1.pdgCode()) == pdg::Code::kDPlus && std::abs(particle2.pdgCode()) == pdg::Code::kDPlus) {
          entryDplusPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                         particle2.eta() - particle1.eta(),
                         particle1.pt(),
                         particle2.pt(),
                         particle1.y(),
                         particle2.y(),
                         o2::analysis::pdg::MassDPlus,
                         o2::analysis::pdg::MassDPlus,
                         candidateType1,
                         candidateType2,
                         2);
          entryDplusPairRecoInfo(originGen1,
                                 originGen2,
                                 matchedGen1,
                                 matchedGen2);
        }
      } // end inner loop
    }   // end outer loop
  }

  /// D0(bar)-D0(bar) correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processDataD0(aod::Collision const& collision,
                     aod::TracksWDca const& tracks,
                     soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&)
  {
    // protection against empty tables to be sliced
    if (selectedD0Candidates.size() <= 1) {
      return;
    }
    analyseMultiplicity(collision, tracks);
    auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    for (const auto& candidate1 : selectedD0CandidatesGrouped) {
      fillInfoHists(candidate1, false, true);
      if (!kinematicCuts<decltype(candidate1), true>(candidate1)) {
        continue;
      }

      registry.fill(HIST("hMass"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt());
      auto candidateType1 = assignCandidateTypeD0<decltype(candidate1), false>(candidate1); // Candidate type attribution.

      for (const auto& candidate2 : selectedD0CandidatesGrouped) {
        if (!kinematicCuts<decltype(candidate2), true>(candidate2)) {
          continue;
        }
        // avoid double counting
        if (candidate1.mRowIndex >= candidate2.mRowIndex) {
          continue;
        }
        auto candidateType2 = assignCandidateTypeD0<decltype(candidate2), false>(candidate2); // Candidate type attribution
        registry.fill(HIST("hSelectionStatus"), 1);
        // fill tables
        entryD0Pair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                    candidate2.eta() - candidate1.eta(),
                    candidate1.pt(),
                    candidate2.pt(),
                    hfHelper.yD0(candidate1),
                    hfHelper.yD0(candidate2),
                    hfHelper.invMassD0ToPiK(candidate1),
                    hfHelper.invMassD0barToKPi(candidate2),
                    candidateType1,
                    candidateType2,
                    0);
      } // end inner loop (Cand2)
    }   // end outer loop (Cand1)
  }
  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processDataD0, "Process data D0", true);

  /// D0(bar)-D0(bar) correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRecD0(aod::Collision const& collision,
                      aod::TracksWDca const& tracks,
                      soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const&)
  {
    // protection against empty tables to be sliced
    if (selectedD0CandidatesMc.size() <= 1) {
      return;
    }
    analyseMultiplicity(collision, tracks);
    auto selectedD0CandidatesGroupedMc = selectedD0CandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    for (const auto& candidate1 : selectedD0CandidatesGroupedMc) {
      if (!kinematicCuts<decltype(candidate1), true>(candidate1)) {
        continue;
      }
      if (std::abs(candidate1.flagMcMatchRec()) == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) {
        fillInfoHists(candidate1, true, true);
      }

      registry.fill(HIST("hMass"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt());
      auto candidateType1 = assignCandidateTypeD0<decltype(candidate1), true>(candidate1); // Candidate type attribution

      int8_t origin1 = 0, matchedRec1 = 0;
      if (!(TESTBIT(candidateType1, TrueD) && TESTBIT(candidateType1, TrueDbar))) { // if our event is not bkg
        // check if it's prompt or non-prompt
        origin1 = candidate1.originMcRec();
        registry.fill(HIST("hOriginMcRec"), origin1);
        // check if it's MC matched
        matchedRec1 = candidate1.flagMcMatchRec();
        registry.fill(HIST("hMatchedMcRec"), matchedRec1);
      }

      for (const auto& candidate2 : selectedD0CandidatesGroupedMc) {
        if (!kinematicCuts<decltype(candidate2), true>(candidate2)) {
          continue;
        }
        // avoid double counting
        if (candidate1.mRowIndex >= candidate2.mRowIndex) {
          continue;
        }

        auto candidateType2 = assignCandidateTypeD0<decltype(candidate2), true>(candidate2); // Candidate type attribution
        int8_t origin2 = 0, matchedRec2 = 0;
        if (!(TESTBIT(candidateType2, TrueD) && TESTBIT(candidateType2, TrueDbar))) { // if our event is not bkg
          // check if it's prompt or non-prompt
          origin2 = candidate2.originMcRec();
          // check if it's MC matched
          matchedRec2 = candidate2.flagMcMatchRec();
        }
        registry.fill(HIST("hSelectionStatus"), 1);
        // fill tables
        entryD0Pair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                    candidate2.eta() - candidate1.eta(),
                    candidate1.pt(),
                    candidate2.pt(),
                    hfHelper.yD0(candidate1),
                    hfHelper.yD0(candidate2),
                    hfHelper.invMassD0ToPiK(candidate1),
                    hfHelper.invMassD0barToKPi(candidate2),
                    candidateType1,
                    candidateType2,
                    1);
        entryD0PairRecoInfo(origin1,
                            origin2,
                            matchedRec1,
                            matchedRec2);
      } // end inner loop
    }   // end outer loop
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcRecD0, "Process D0 Mc Reco mode", false);

  /// D0(bar)-D0(bar) correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGenD0(aod::McCollision const&,
                      McParticlesPlus2Prong const& mcParticles)
  {
    analyseMcGen(mcParticles);
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcGenD0, "Process D0 Mc Gen mode", false);

  /// Dplus(minus)-Dplus(minus) correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processDataDPlus(aod::Collision const& collision,
                        aod::TracksWDca const& tracks,
                        soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const&)
  {
    // protection against empty tables to be sliced
    if (selectedDPlusCandidates.size() <= 1) {
      return;
    }
    analyseMultiplicity(collision, tracks);
    auto selectedDPlusCandidatesGrouped = selectedDPlusCandidates->sliceByCached(o2::aod::hf_cand::collisionId, collision.globalIndex(), cache);
    for (const auto& candidate1 : selectedDPlusCandidatesGrouped) {
      if (!kinematicCuts<decltype(candidate1), false>(candidate1)) {
        continue;
      }

      int outerParticleSign = 1; // Dplus
      auto outerSecondTrack = candidate1.prong1();
      if (outerSecondTrack.sign() == 1) {
        outerParticleSign = -1; // Dminus (second daughter track is positive)
      }

      auto candidateType1 = assignCandidateTypeDPlus<decltype(candidate1), false>(candidate1, outerParticleSign);
      fillInfoHists(candidate1, false, false);
      registry.fill(HIST("hMass"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt());

      for (const auto& candidate2 : selectedDPlusCandidatesGrouped) {
        if (!kinematicCuts<decltype(candidate2), false>(candidate2)) {
          continue;
        }
        // avoid double counting
        if (candidate1.mRowIndex >= candidate2.mRowIndex) {
          continue;
        }

        int innerParticleSign = 1; // Dplus
        auto innerSecondTrack = candidate2.prong1();
        if (innerSecondTrack.sign() == 1) {
          innerParticleSign = -1; // Dminus (second daughter track is positive)
        }
        auto candidateType2 = assignCandidateTypeDPlus<decltype(candidate2), false>(candidate2, innerParticleSign);
        registry.fill(HIST("hSelectionStatus"), 1);
        // fill tables
        entryDplusPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                       candidate2.eta() - candidate1.eta(),
                       candidate1.pt(),
                       candidate2.pt(),
                       hfHelper.yDplus(candidate1),
                       hfHelper.yDplus(candidate2),
                       hfHelper.invMassDplusToPiKPi(candidate1),
                       hfHelper.invMassDplusToPiKPi(candidate2),
                       candidateType1,
                       candidateType2,
                       0);
      } // end inner loop (cand2)
    }   // end outer loop (cand1)
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processDataDPlus, "Process Data DPlus", false);

  /// Dplus(minus)-Dplus(minus) correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRecDPlus(aod::Collision const& collision,
                         aod::TracksWDca const& tracks,
                         soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec> const&)
  {
    // protection against empty tables to be sliced
    if (selectedDPlusCandidatesMc.size() <= 1) {
      return;
    }
    analyseMultiplicity(collision, tracks);
    auto selectedDPlusCandidatesGroupedMc = selectedDPlusCandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    for (const auto& candidate1 : selectedDPlusCandidatesGroupedMc) {
      if (!kinematicCuts<decltype(candidate1), false>(candidate1)) {
        continue;
      }

      int outerParticleSign = 1; // Dplus
      auto outerSecondTrack = candidate1.prong1();
      if (outerSecondTrack.sign() == 1) {
        outerParticleSign = -1; // Dminus (second daughter track is positive)
      }

      auto candidateType1 = assignCandidateTypeDPlus<decltype(candidate1), true>(candidate1, outerParticleSign);
      registry.fill(HIST("hMass"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt());

      int8_t origin1 = 0, matchedRec1 = 0;
      if (!(TESTBIT(candidateType1, TrueD) && TESTBIT(candidateType1, TrueDbar))) { // if our event is not bkg
        // check if it's prompt or non-prompt
        origin1 = candidate1.originMcRec();
        registry.fill(HIST("hOriginMcRec"), origin1);
        // check if it's MC matched
        matchedRec1 = candidate1.flagMcMatchRec();
        registry.fill(HIST("hMatchedMcRec"), matchedRec1);
      }

      for (const auto& candidate2 : selectedDPlusCandidatesGroupedMc) {
        if (!kinematicCuts<decltype(candidate2), false>(candidate2)) {
          continue;
        }
        // avoid double counting
        if (candidate1.mRowIndex >= candidate2.mRowIndex) {
          continue;
        }

        int innerParticleSign = 1; // Dplus
        auto innerSecondTrack = candidate2.prong1();
        if (innerSecondTrack.sign() == 1) {
          innerParticleSign = -1; // Dminus (second daughter track is positive)
        }

        uint candidateType2 = assignCandidateTypeDPlus<decltype(candidate2), true>(candidate2, innerParticleSign);
        int8_t origin2 = 0, matchedRec2 = 0;
        if (!(TESTBIT(candidateType2, TrueD) && TESTBIT(candidateType2, TrueDbar))) { // if our event is not bkg
          // check if it's prompt or non-prompt
          origin2 = candidate1.originMcRec();
          // check if it's MC matched
          matchedRec2 = candidate1.flagMcMatchRec();
        }
        registry.fill(HIST("hSelectionStatus"), 1);
        // fill tables
        entryDplusPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                       candidate2.eta() - candidate1.eta(),
                       candidate1.pt(),
                       candidate2.pt(),
                       hfHelper.yDplus(candidate1),
                       hfHelper.yDplus(candidate2),
                       hfHelper.invMassDplusToPiKPi(candidate1),
                       hfHelper.invMassDplusToPiKPi(candidate2),
                       candidateType1,
                       candidateType2,
                       1);
        entryDplusPairRecoInfo(origin1,
                               origin2,
                               matchedRec1,
                               matchedRec2);
      } // end inner loop (cand2)
    }   // end outer loop (cand1)
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcRecDPlus, "Process DPlus Mc Reco", false);

  /// Dplus(minus)-Dplus(minus) correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGenDPlus(aod::McCollision const&,
                         McParticlesPlus3Prong const& mcParticles)
  {
    analyseMcGen(mcParticles);
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcGenDPlus, "Process DPlus Mc Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDMesonPairs>(cfgc)};
}
