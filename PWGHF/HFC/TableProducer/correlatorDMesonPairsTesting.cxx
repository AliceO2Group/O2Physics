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
/// \brief D0(bar) correlator task - data-like, MC-reco and MC-kine analyses.
/// \note Temporary code for tests of D0-D0 correlations
///
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DMesonPairsTablesTesting.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
enum CandidateType {
  SelectedD = 0, // This particle is selected as a D
  SelectedDbar,  // This particle is selected as a Dbar
  TrueD,         // This particle is a true D
  TrueDbar       // This particle is a true Dbar
};

enum PairTypeOfSelMassSel {
  DD = 0,   // This is a D0-D0 pair
  DbarDbar, // This is a D0bar-D0bar pair
  DDbar,    // This is a D0-D0bar pair
  DbarD     // This is a D0bar-D0 pair
};
} // namespace

using McParticlesPlus2Prong = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

struct HfCorrelatorDMesonPairsTesting {
  SliceCache cache;
  Preslice<aod::HfCand2Prong> perCol2Prong = aod::hf_cand::collisionId;

  Produces<aod::D0PairTesting> entryD0Pair;
  Produces<aod::D0PairMcInfoTesting> entryD0PairMcInfo;
  Produces<aod::D0PairMcGenTesting> entryD0PairMcGen;
  Produces<aod::D0PairMcGenInfoTesting> entryD0PairMcGenInfo;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "maxmum |y| of D0 candidates"};
  Configurable<float> ptCandMin{"ptCandMin", -1., "minimum pT of D0 candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots"};

  HfHelper hfHelper;

  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPiExt, aod::TracksPidKaExt>;

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> selectedD0CandidatesMc = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  HistogramConfigSpec hTH1Pt{HistType::kTH1F, {{180, 0., 36.}}};
  HistogramConfigSpec hTH1Y{HistType::kTH1F, {{100, -5., 5.}}};
  HistogramConfigSpec hTH1Phi{HistType::kTH1F, {{32, 0., o2::constants::math::TwoPI}}};
  HistogramConfigSpec hTH2Pid{HistType::kTH2F, {{500, 0., 10.}, {400, -20., 20.}}};

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "D meson candidates;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtCandAfterCut", "D meson candidates after pT cut;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng0", "D meson candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng1", "D meson candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEta", "D meson candidates;candidate #it{#eta};entries", hTH1Y},
     {"hPhi", "D meson candidates;candidate #it{#varphi};entries", hTH1Phi},
     {"hY", "D meson candidates;candidate #it{y};entries", hTH1Y},
     // MC Gen plots
     {"hPtCandMcGen", "D meson candidates MC Gen;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtCandAfterCutMcGen", "D meson candidates after pT cut;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEtaMcGen", "D meson candidates MC Gen;candidate #it{#eta};entries", hTH1Y},
     {"hPhiMcGen", "D meson candidates MC Gen;candidate #it{#varphi};entries", hTH1Phi},
     // PID plots ----- Not definitively here
     {"PID/hTofNSigmaPi", "(TOFsignal-time#pi)/tofSigPid;p[GeV/c];(TOFsignal-time#pi)/tofSigPid", hTH2Pid},
     {"PID/hTofNSigmaKa", "(TOFsignal-timeK)/tofSigPid;p[GeV/c];(TOFsignal-timeK)/tofSigPid", hTH2Pid},
     {"PID/hTpcNSigmaPi", "(TPCsignal-time#pi)/tpcSigPid;p[GeV/c];(TPCsignal-time#pi)/tpcSigPid", hTH2Pid},
     {"PID/hTpcNSigmaKa", "(TPCsignal-timeK)/tpcSigPid;p[GeV/c];(TPCsignal-timeK)/tpcSigPid", hTH2Pid},
     {"PID/hTpcTofNSigmaPi", "(TPC+TOFsignal-time#pi)/tpcTofSigPid;p[GeV/#it{c}];(TPC+TOFsignal-time#pi)/tpcTofSigPid", hTH2Pid},
     {"PID/hTpcTofNSigmaKa", "(TPC+TOFsignal-timeK)/tpcTofSigPid;p[GeV/c];(TPC+TOFsignal-timeK)/tpcTofSigPid", hTH2Pid}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    constexpr int kNBinsSelStatus = 25;
    std::string labels[kNBinsSelStatus];

    labels[0] = "total # of Selected pairs";
    // Cand1 analysis
    labels[1] = "total # of Selected Cands 1";
    labels[2] = "# of Selected D Cand 1 ONLY";
    labels[3] = "# of Selected Dbar Cand 1 ONLY";
    labels[4] = "# of Selected Simultaneous D + Dbar Cand 1";
    labels[5] = "# of True D Cand 1";
    labels[6] = "# of True Dbar Cand 1";
    // Cand2 analysis
    labels[7] = "total # of Selected Cands 2";
    labels[8] = "# of Selected D Cand 2 ONLY";
    labels[9] = "# of Selected Dbar Cand 2 ONLY";
    labels[10] = "# of Selected Simultaneous D + Dbar Cand 2";
    labels[11] = "# of True D Cand 2";
    labels[12] = "# of True Dbar Cand 2";
    // Pair analysis
    labels[13] = "# of D+D Pairs";
    labels[14] = "# of Dbar+Dbar Pairs";
    labels[15] = "# of D+Dbar Pairs";
    labels[16] = "# of Dbar+D Pairs";
    labels[17] = "# of D+D ONLY Pairs";
    labels[18] = "# of Dbar+Dbar ONLY Pairs";
    labels[19] = "# of D+Dbar ONLY Pairs";
    labels[20] = "# of Dbar+D ONLY Pairs";
    // True pair analysis
    labels[21] = "# of True D+D Pairs";
    labels[22] = "# of True Dbar+Dbar Pairs";
    labels[23] = "# of True D+Dbar Pairs";
    labels[24] = "# of True Dbar+D Pairs";

    AxisSpec axisSelStatus = {kNBinsSelStatus, 0.5, kNBinsSelStatus + 0.5, ""};
    registry.add("hSelectionStatus", "D Meson candidates;selection status;entries", HistType::kTH1F, {axisSelStatus});
    registry.add("hSelectionStatusMcGen", "D Meson candidates MC Gen;selection status;entries", HistType::kTH1F, {axisSelStatus});

    for (int iBin = 0; iBin < kNBinsSelStatus; iBin++) {
      registry.get<TH1>(HIST("hSelectionStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      registry.get<TH1>(HIST("hSelectionStatusMcGen"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    constexpr int kNBinsMatching = 8;
    std::string labelsMatching[kNBinsMatching];
    // Cand1 analysis
    labelsMatching[0] = "total # of Cand 1";
    labelsMatching[1] = "# of matched D Cand 1";
    labelsMatching[2] = "# of matched Dbar Cand 1";
    labelsMatching[3] = "# of unmatched Cand 1";
    // Cand2 analysis
    labelsMatching[4] = "total # of Cand 2";
    labelsMatching[5] = "# of matched D Cand 2";
    labelsMatching[6] = "# of matched Dbar Cand 2";
    labelsMatching[7] = "# of unmatched Cand 2";

    AxisSpec axisMatching = {kNBinsMatching, 0.5, kNBinsMatching + 0.5, ""};
    registry.add("hMatchingMcRec", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisMatching});
    registry.add("hMatchingMcGen", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisMatching});

    for (int iBin = 0; iBin < kNBinsMatching; iBin++) {
      registry.get<TH1>(HIST("hMatchingMcRec"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMatching[iBin].data());
      registry.get<TH1>(HIST("hMatchingMcGen"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMatching[iBin].data());
    }

    constexpr int kNBinsSinglePart = 6;
    std::string labelsSinglePart[kNBinsSinglePart];
    // Candidate analysis
    labelsSinglePart[0] = "total # of Candidates";
    labelsSinglePart[1] = "# of selected D";
    labelsSinglePart[2] = "# of selected Dbar";
    labelsSinglePart[3] = "# of selected D and Dbar";
    labelsSinglePart[4] = "# of true D";
    labelsSinglePart[5] = "# of true Dbar";

    AxisSpec axisSinglePart = {kNBinsSinglePart, 0.5, kNBinsSinglePart + 0.5, ""};
    registry.add("hStatusSinglePart", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisSinglePart});
    registry.add("hStatusSinglePartMcGen", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisSinglePart});

    for (int iBin = 0; iBin < kNBinsSinglePart; iBin++) {
      registry.get<TH1>(HIST("hStatusSinglePart"))->GetXaxis()->SetBinLabel(iBin + 1, labelsSinglePart[iBin].data());
      registry.get<TH1>(HIST("hStatusSinglePartMcGen"))->GetXaxis()->SetBinLabel(iBin + 1, labelsSinglePart[iBin].data());
    }

    AxisSpec axisInputD0 = {200, -0.5, 199.5};
    registry.add("hInputCheckD0", "Check on input D0 meson candidates/event", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0bar", "Check on input D0bar meson candidates/event", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0AndD0bar", "Check on input D0 & D0bar meson candidates/event", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0OrD0bar", "Check on input D0 | D0bar meson candidates/event", {HistType::kTH1F, {axisInputD0}});
    // MC Gen
    registry.add("hInputCheckD0McGen", "Check on input D0 meson candidates/event MC Gen", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0barMcGen", "Check on input D0bar meson candidates/event MC Gen", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0AndD0barMcGen", "Check on input D0 & D0bar meson candidates/event MC Gen", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0OrD0barMcGen", "Check on input D0 | D0bar meson candidates/event MC Gen", {HistType::kTH1F, {axisInputD0}});

    registry.add("hMass", "D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  /// Sets bits to select candidate type for D0
  /// SelectedD and SelectedDbar bits look at whether the candidate passed the selection flags.
  /// \param candidate is candidate
  /// \return bitmap with type of candidate
  template <typename T>
  uint8_t assignCandidateTypeD0(const T& candidate)
  {
    uint8_t candidateType(0);
    if (candidate.isSelD0() >= selectionFlagD0) {
      SETBIT(candidateType, SelectedD);
    }
    if (candidate.isSelD0bar() >= selectionFlagD0bar) {
      SETBIT(candidateType, SelectedDbar);
    }
    return candidateType;
  }

  /// Sets bits to select true candidate type for D0
  /// SelectedD and SelectedDbar bits look at whether the candidate passed the selection flags.
  /// \param candidate is candidate
  /// \return bitmap with true type of candidate
  template <typename T>
  uint8_t assignCandidateTypeD0True(const T& candidate)
  {
    uint8_t candidateType(0);
    if (candidate.flagMcMatchRec() == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) { // matched as D0
      SETBIT(candidateType, TrueD);
    }
    if (candidate.flagMcMatchRec() == -(1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) { // matched as D0bar
      SETBIT(candidateType, TrueDbar);
    }
    return candidateType;
  }

  /// Sets bits to select candidate type at generator level
  /// \param particle is particle
  /// \return bitmap with type of gen particle - they are selected as True
  template <typename T>
  uint8_t assignParticleTypeD0Gen(const T& particle)
  {
    uint8_t particleType(0);
    if (particle.pdgCode() == Pdg::kD0) { // just checking the particle PDG, not the decay channel
      SETBIT(particleType, TrueD);
    }
    if (particle.pdgCode() == (-Pdg::kD0)) { // just checking the particle PDG, not the decay channel
      SETBIT(particleType, TrueDbar);
    }
    return particleType;
  }

  /// Check if the candidate is a D
  /// \param candidate is candidate
  /// \return true if is a D
  bool isD(const uint8_t candidateType)
  {
    return (TESTBIT(candidateType, SelectedD));
  }

  /// Check if the candidate is a Dbar
  /// \param candidate is candidate
  /// \return true if is a Dbar
  bool isDbar(const uint8_t candidateType)
  {
    return (TESTBIT(candidateType, SelectedDbar));
  }

  /// Check if the candidate is a true D
  /// \param candidate is candidate
  /// \return true if is a true D
  bool isTrueD(const uint8_t candidateType)
  {
    return (TESTBIT(candidateType, TrueD));
  }

  /// Check if the candidate is a true Dbar
  /// \param candidate is candidate
  /// \return true if is a true Dbar
  bool isTrueDbar(const uint8_t candidateType)
  {
    return (TESTBIT(candidateType, TrueDbar));
  }

  /// Fill PID related plots for D0 and D0bar
  /// \param candidate is candidate
  template <typename T>
  void AnalysePid(const T& candidate)
  {
    auto prong0 = candidate.template prong0_as<TracksWPid>();
    auto prong1 = candidate.template prong1_as<TracksWPid>();
    if (candidate.isSelD0() >= selectionFlagD0) {
      registry.fill(HIST("PID/hTofNSigmaPi"), candidate.ptProng0(), prong0.tofNSigmaPi());
      registry.fill(HIST("PID/hTofNSigmaKa"), candidate.ptProng1(), prong1.tofNSigmaKa());
      registry.fill(HIST("PID/hTpcNSigmaPi"), candidate.ptProng0(), prong0.tpcNSigmaPi());
      registry.fill(HIST("PID/hTpcNSigmaKa"), candidate.ptProng1(), prong1.tpcNSigmaKa());
      registry.fill(HIST("PID/hTpcTofNSigmaPi"), candidate.ptProng0(), prong0.tpcTofNSigmaPi());
      registry.fill(HIST("PID/hTpcTofNSigmaKa"), candidate.ptProng1(), prong1.tpcTofNSigmaKa());
    }
    if (candidate.isSelD0bar() >= selectionFlagD0bar) {
      registry.fill(HIST("PID/hTofNSigmaPi"), candidate.ptProng1(), prong1.tofNSigmaPi());
      registry.fill(HIST("PID/hTofNSigmaKa"), candidate.ptProng0(), prong0.tofNSigmaKa());
      registry.fill(HIST("PID/hTpcNSigmaPi"), candidate.ptProng1(), prong1.tpcNSigmaPi());
      registry.fill(HIST("PID/hTpcNSigmaKa"), candidate.ptProng0(), prong0.tpcNSigmaKa());
      registry.fill(HIST("PID/hTpcTofNSigmaPi"), candidate.ptProng1(), prong1.tpcTofNSigmaPi());
      registry.fill(HIST("PID/hTpcTofNSigmaKa"), candidate.ptProng0(), prong0.tpcTofNSigmaKa());
    }
  }

  /// Fill counters for D0 and D0bar
  /// \param selectedD0Candidates contains all D0 candidates
  template <typename T>
  void GetCountersPerEvent(const T& selectedD0Candidates)
  {
    int nDevent = 0, nDbarevent = 0, nDDbarevent = 0, nDorDbarevent = 0;
    for (const auto& candidate : selectedD0Candidates) {
      // Get counters per event
      auto candidateType1 = assignCandidateTypeD0(candidate); // Candidate type attribution
      registry.fill(HIST("hPtCand"), candidate.pt());
      if (abs(hfHelper.yD0(candidate)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hEta"), candidate.eta());
      registry.fill(HIST("hPhi"), candidate.phi());
      registry.fill(HIST("hY"), candidate.y(MassD0));
      registry.fill(HIST("hPtCandAfterCut"), candidate.pt());
      registry.fill(HIST("hMass"), hfHelper.invMassD0ToPiK(candidate), candidate.pt());

      bool isDCand1 = isD(candidateType1);
      bool isDbarCand1 = isDbar(candidateType1);
      if (isDCand1) {
        nDevent++;
      }
      if (isDbarCand1) {
        nDbarevent++;
      }
      if (isDCand1 && isDbarCand1) {
        nDDbarevent++;
      }
      if (isDCand1 || isDbarCand1) {
        nDorDbarevent++;
      }

      // Fill selection status single particle
      registry.fill(HIST("hStatusSinglePart"), 1);
      if (isDCand1 && !isDbarCand1) {
        registry.fill(HIST("hStatusSinglePart"), 2);
      } else if (isDbarCand1 && !isDCand1) {
        registry.fill(HIST("hStatusSinglePart"), 3);
      } else if (isDCand1 && isDbarCand1) {
        registry.fill(HIST("hStatusSinglePart"), 4);
      }
    }
    if (nDevent > 0) {
      registry.fill(HIST("hInputCheckD0"), nDevent);
    }
    if (nDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0bar"), nDbarevent);
    }
    if (nDDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0AndD0bar"), nDDbarevent);
    }
    if (nDorDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0OrD0bar"), nDorDbarevent);
    }
  }

  /// Fill selection status histogram
  /// \param candidate1 is the first candidate of the pair
  /// \param candidate2 is the second candidate of the pair
  template <typename T>
  void fillEntry(const T& candidate1, const T& candidate2, const bool& isDCand1, const bool& isDbarCand1,
                 const bool& isDCand2, const bool& isDbarCand2)
  {

    /// Fill information on the D candidates
    registry.fill(HIST("hSelectionStatus"), 2);
    if (isDCand1 && !isDbarCand1) {
      registry.fill(HIST("hSelectionStatus"), 3);
    } else if (isDbarCand1 && !isDCand1) {
      registry.fill(HIST("hSelectionStatus"), 4);
    } else if (isDCand1 && isDbarCand1) {
      registry.fill(HIST("hSelectionStatus"), 5);
    }

    registry.fill(HIST("hSelectionStatus"), 8);
    if (isDCand2 && !isDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 9);
    } else if (isDbarCand2 && !isDCand2) {
      registry.fill(HIST("hSelectionStatus"), 10);
    } else if (isDCand2 && isDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 11);
    }

    /// Collect information on the D pairs
    uint8_t pairType(0);
    registry.fill(HIST("hSelectionStatus"), 1);
    float yCand1 = hfHelper.yD0(candidate1);
    float yCand2 = hfHelper.yD0(candidate2);
    float massDCand1 = hfHelper.invMassD0ToPiK(candidate1);
    float massDbarCand1 = hfHelper.invMassD0barToKPi(candidate1);
    float massDCand2 = hfHelper.invMassD0ToPiK(candidate2);
    float massDbarCand2 = hfHelper.invMassD0barToKPi(candidate1);
    if (isDCand1 && isDCand2) {
      SETBIT(pairType, DD);
      registry.fill(HIST("hSelectionStatus"), 14);
      if ((!isDbarCand1) && (!isDbarCand2)) {
        registry.fill(HIST("hSelectionStatus"), 18);
      }
    }
    if (isDbarCand1 && isDbarCand2) {
      SETBIT(pairType, DbarDbar);
      registry.fill(HIST("hSelectionStatus"), 15);
      if ((!isDCand1) && (!isDCand2)) {
        registry.fill(HIST("hSelectionStatus"), 19);
      }
    }
    if (isDCand1 && isDbarCand2) {
      SETBIT(pairType, DDbar);
      registry.fill(HIST("hSelectionStatus"), 16);
      if ((!isDbarCand1) && (!isDCand2)) {
        registry.fill(HIST("hSelectionStatus"), 20);
      }
    }
    if (isDbarCand1 && isDCand2) {
      SETBIT(pairType, DbarD);
      registry.fill(HIST("hSelectionStatus"), 17);
      if ((!isDCand1) && (!isDbarCand2)) {
        registry.fill(HIST("hSelectionStatus"), 21);
      }
    }

    entryD0Pair(candidate1.pt(), candidate2.pt(), yCand1, yCand2, massDCand1, massDbarCand1, massDCand2, massDbarCand2, pairType);
  }

  /// D0(bar)-D0(bar) correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(aod::Collision const& collision,
                   soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates, TracksWPid const&)
  {
    for (const auto& candidate : candidates) {
      AnalysePid(candidate);
    }
    auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    GetCountersPerEvent(selectedD0CandidatesGrouped);
    // protection against empty tables to be sliced
    if (selectedD0Candidates.size() <= 1) {
      return;
    }
    for (const auto& candidate1 : selectedD0CandidatesGrouped) {
      if (abs(hfHelper.yD0(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }

      auto candidateType1 = assignCandidateTypeD0(candidate1); // Candidate type attribution
      bool isDCand1 = isD(candidateType1);
      bool isDbarCand1 = isDbar(candidateType1);

      for (auto candidate2 = candidate1 + 1; candidate2 != selectedD0CandidatesGrouped.end(); ++candidate2) {
        if (abs(hfHelper.yD0(candidate2)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate2.pt() < ptCandMin) {
          continue;
        }
        auto candidateType2 = assignCandidateTypeD0(candidate2); // Candidate type attribution

        bool isDCand2 = isD(candidateType2);
        bool isDbarCand2 = isDbar(candidateType2);

        fillEntry(candidate1, candidate2, isDCand1, isDbarCand1, isDCand2, isDbarCand2);
      } // end inner loop (Cand2)
    }   // end outer loop (Cand1)
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairsTesting, processData, "Process data mode", true);

  void processMcRec(aod::Collision const& collision, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const& candidates, TracksWPid const&)
  {
    for (const auto& candidate : candidates) {
      AnalysePid(candidate);
    }
    auto selectedD0CandidatesGroupedMc = selectedD0CandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    GetCountersPerEvent(selectedD0CandidatesGroupedMc);
    // protection against empty tables to be sliced
    if (selectedD0CandidatesMc.size() <= 1) {
      return;
    }
    for (const auto& candidate1 : selectedD0CandidatesGroupedMc) {
      if (abs(hfHelper.yD0(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }

      auto candidateType1 = assignCandidateTypeD0(candidate1); // Candidate type attribution

      bool isDCand1 = isD(candidateType1);
      bool isDbarCand1 = isDbar(candidateType1);

      // assign true type
      auto candidateTypeTrue1 = assignCandidateTypeD0True(candidate1);
      bool isTrueDCand1 = isTrueD(candidateTypeTrue1);
      bool isTrueDbarCand1 = isTrueDbar(candidateTypeTrue1);

      int8_t matchedRec1 = candidate1.flagMcMatchRec();
      int8_t originRec1 = candidate1.originMcRec();

      if (isTrueDCand1) {
        registry.fill(HIST("hStatusSinglePart"), 5);
      } else if (isTrueDbarCand1) {
        registry.fill(HIST("hStatusSinglePart"), 6);
      }

      for (auto candidate2 = candidate1 + 1; candidate2 != selectedD0CandidatesGroupedMc.end(); ++candidate2) {
        if (abs(hfHelper.yD0(candidate2)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate2.pt() < ptCandMin) {
          continue;
        }
        auto candidateType2 = assignCandidateTypeD0(candidate2); // Candidate type attribution

        bool isDCand2 = isD(candidateType2);
        bool isDbarCand2 = isDbar(candidateType2);
        // assign true type
        auto candidateTypeTrue2 = assignCandidateTypeD0True(candidate2);
        bool isTrueDCand2 = isTrueD(candidateTypeTrue2);
        bool isTrueDbarCand2 = isTrueDbar(candidateTypeTrue2);

        int8_t matchedRec2 = candidate2.flagMcMatchRec();
        int8_t originRec2 = candidate2.originMcRec();

        // Fill hMatchingMcRec - Cand 1
        registry.fill(HIST("hMatchingMcRec"), 1);
        if (matchedRec1 == 1) {
          registry.fill(HIST("hMatchingMcRec"), 2);
        } else if (matchedRec1 == -1) {
          registry.fill(HIST("hMatchingMcRec"), 3);
        } else if (matchedRec1 == 0) {
          registry.fill(HIST("hMatchingMcRec"), 4);
        }
        // Fill hMatchingMcRec - Cand 2
        registry.fill(HIST("hMatchingMcRec"), 5);
        if (matchedRec2 == 1) {
          registry.fill(HIST("hMatchingMcRec"), 6);
        } else if (matchedRec2 == -1) {
          registry.fill(HIST("hMatchingMcRec"), 7);
        } else if (matchedRec2 == 0) {
          registry.fill(HIST("hMatchingMcRec"), 8);
        }
        fillEntry(candidate1, candidate2, isDCand1, isDbarCand1, isDCand2, isDbarCand2);
        entryD0PairMcInfo(originRec1, originRec2, matchedRec1, matchedRec2);

        if (isTrueDCand1) {
          registry.fill(HIST("hSelectionStatus"), 6);
        } else if (isTrueDbarCand1) {
          registry.fill(HIST("hSelectionStatus"), 7);
        }
        if (isTrueDCand2) {
          registry.fill(HIST("hSelectionStatus"), 12);
        } else if (isTrueDbarCand2) {
          registry.fill(HIST("hSelectionStatus"), 13);
        }
        if (isTrueDCand1 && isTrueDCand2) {
          registry.fill(HIST("hSelectionStatus"), 22);
        } else if (isTrueDbarCand1 && isTrueDbarCand2) {
          registry.fill(HIST("hSelectionStatus"), 23);
        } else if (isTrueDCand1 && isTrueDbarCand2) {
          registry.fill(HIST("hSelectionStatus"), 24);
        } else if (isTrueDbarCand1 && isTrueDCand2) {
          registry.fill(HIST("hSelectionStatus"), 25);
        }
      } // end inner loop (Cand2)
    }   // end outer loop (Cand1)
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairsTesting, processMcRec, "Process Mc reco mode", false);

  void processMcGen(aod::McCollision const&, McParticlesPlus2Prong const& mcParticles)
  {
    // Get counters per event
    int nDevent = 0, nDbarevent = 0, nDDbarevent = 0, nDorDbarevent = 0;
    for (const auto& particle : mcParticles) {
      // check if the particle is D0 or D0bar
      if (std::abs(particle.pdgCode()) != Pdg::kD0) {
        continue;
      }
      if (abs(particle.y()) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
        continue;
      }
      auto particleType = assignParticleTypeD0Gen(particle); // Candidate type attribution
      bool isDParticle = isTrueD(particleType);
      bool isDbarParticle = isTrueDbar(particleType);
      if (isDParticle) {
        nDevent++;
      }
      if (isDbarParticle) {
        nDbarevent++;
      }
      if (isDParticle && isDbarParticle) {
        nDDbarevent++;
      }
      if (isDParticle || isDbarParticle) {
        nDorDbarevent++;
      }
    }
    if (nDevent > 0) {
      registry.fill(HIST("hInputCheckD0McGen"), nDevent);
    }
    if (nDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0barMcGen"), nDbarevent);
    }
    if (nDDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0AndD0barMcGen"), nDDbarevent);
    }
    if (nDorDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0OrD0barMcGen"), nDorDbarevent);
    }

    auto massD = MassD0;
    auto massDbar = MassD0Bar;

    for (const auto& particle1 : mcParticles) {
      // check if the particle is D0 or D0bar
      if (std::abs(particle1.pdgCode()) != Pdg::kD0) {
        continue;
      }
      registry.fill(HIST("hPtCandMcGen"), particle1.pt());

      if (abs(particle1.y()) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hEtaMcGen"), particle1.eta());
      registry.fill(HIST("hPhiMcGen"), particle1.phi());
      registry.fill(HIST("hPtCandAfterCutMcGen"), particle1.pt());

      auto particleType1 = assignParticleTypeD0Gen(particle1); // Candidate sign attribution
      bool isDParticle1 = isTrueD(particleType1);
      bool isDbarParticle1 = isTrueDbar(particleType1);

      // check if it's MC matched
      int8_t matchedGen1 = particle1.flagMcMatchGen();
      // check origin
      int8_t originGen1 = particle1.originMcGen();

      // Fill selection status single particle
      registry.fill(HIST("hStatusSinglePartMcGen"), 1);
      if (isDParticle1 && !isDbarParticle1) {
        registry.fill(HIST("hStatusSinglePartMcGen"), 2);
      } else if (isDbarParticle1 && !isDParticle1) {
        registry.fill(HIST("hStatusSinglePartMcGen"), 3);
      } else if (isDParticle1 && isDbarParticle1) {
        registry.fill(HIST("hStatusSinglePartMcGen"), 4);
      }

      for (auto particle2 = particle1 + 1; particle2 != mcParticles.end(); ++particle2) {
        // check if the particle is D0 or D0bar
        if (std::abs(particle2.pdgCode()) != Pdg::kD0) {
          continue;
        }
        if (abs(particle2.y()) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle2.pt() < ptCandMin) {
          continue;
        }
        // Candidate sign attribution.
        auto particleType2 = assignParticleTypeD0Gen(particle2);
        bool isDParticle2 = isTrueD(particleType2);
        bool isDbarParticle2 = isTrueDbar(particleType2);

        // check if it's MC matched
        int8_t matchedGen2 = particle2.flagMcMatchGen();
        // check origin
        int8_t originGen2 = particle2.originMcGen();

        // Fill hMatchingMcGen - Cand 1
        registry.fill(HIST("hMatchingMcGen"), 1);
        if (matchedGen1 == 1) {
          registry.fill(HIST("hMatchingMcGen"), 2);
        } else if (matchedGen1 == -1) {
          registry.fill(HIST("hMatchingMcGen"), 3);
        } else if (matchedGen1 == 0) {
          registry.fill(HIST("hMatchingMcGen"), 4);
        }
        // Fill hMatchingMcRec - Cand 2
        registry.fill(HIST("hMatchingMcGen"), 5);
        if (matchedGen2 == 1) {
          registry.fill(HIST("hMatchingMcGen"), 6);
        } else if (matchedGen2 == -1) {
          registry.fill(HIST("hMatchingMcGen"), 7);
        } else if (matchedGen2 == 0) {
          registry.fill(HIST("hMatchingMcGen"), 8);
        }

        // Fill particle1's Selection Status
        registry.fill(HIST("hSelectionStatusMcGen"), 2);
        if (isDParticle1 && !isDbarParticle1) {
          registry.fill(HIST("hSelectionStatusMcGen"), 6);
        } else if (isDbarParticle1 && !isDParticle1) {
          registry.fill(HIST("hSelectionStatusMcGen"), 7);
        }

        // Fill particle2's Selection Status
        registry.fill(HIST("hSelectionStatusMcGen"), 8);
        if (isDParticle2 && !isDbarParticle2) {
          registry.fill(HIST("hSelectionStatusMcGen"), 12);
        } else if (isDbarParticle2 && !isDParticle2) {
          registry.fill(HIST("hSelectionStatusMcGen"), 13);
        }

        // Fill pair Selection Status
        uint8_t pairType(0);
        registry.fill(HIST("hSelectionStatusMcGen"), 1);

        if (isDParticle1 && isDParticle2) {
          SETBIT(pairType, DD);
          registry.fill(HIST("hSelectionStatusMcGen"), 22);
        }
        if (isDbarParticle1 && isDbarParticle2) {
          SETBIT(pairType, DbarDbar);
          registry.fill(HIST("hSelectionStatusMcGen"), 23);
        }
        if (isDParticle1 && isDbarParticle2) {
          SETBIT(pairType, DDbar);
          registry.fill(HIST("hSelectionStatusMcGen"), 24);
        }
        if (isDbarParticle1 && isDParticle2) {
          SETBIT(pairType, DbarD);
          registry.fill(HIST("hSelectionStatusMcGen"), 25);
        }

        // Fill pair Selection Status
        entryD0PairMcGen(particle1.pt(), particle2.pt(), particle1.y(), particle2.y(), massD, massDbar, massD, massDbar, pairType);
        entryD0PairMcGenInfo(originGen1, originGen2, matchedGen1, matchedGen2);

      } // end inner loop
    }   // end outer loop
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairsTesting, processMcGen, "Process D0 Mc Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDMesonPairsTesting>(cfgc)};
}
