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
/// \brief D+-Hadrons correlator task - data-like, MC-reco and MC-Gen analyses
/// \author Shyam Kumar <shyam.kumar@cern.ch>

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies

double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

/// definition of variables for Dplus hadron pairs (in data-like, MC-reco and MC-kine tasks)
const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_dplus_to_pi_k_pi::NBinsPt;
const std::vector<double> efficiencyDmeson(npTBinsMassAndEfficiency + 1);

// definition of ME variables
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
using BinningTypeMcGen = ColumnBinningPolicy<aod::mccollision::PosZ, o2::aod::mult::MultMCFT0A>;

// Code to select a Dmeson in a collision
struct HfCorrelatorDplusHadronsDplusSelection {
  Produces<aod::DmesonSelection> dplusSel;

  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<bool> doSelDplusCollision{"doSelDplusCollision", true, "Select collisions with at least one D+"};
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for D+"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  HfHelper hfHelper;
  SliceCache cache;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CandidatesDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandidatesDplusMcRec = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using CandDplusMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  // filter on selection of Dplus meson and decay channel Dplus->KPiPi
  Filter dplusFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi)) != static_cast<uint8_t>(0)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  void processDplusSelectionData(SelCollisions::iterator const& collision,
                                 CandidatesDplusData const& candidates)
  {
    bool isSelColl = true;
    bool isDplusFound = true;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    if (doSelDplusCollision) {
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yDplus(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          isDplusFound = false;
          continue;
        }
        isDplusFound = true;
        break;
      }
    }
    if (useSel8) {
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
    }
    isSelColl = isDplusFound && isSel8 && isNosameBunchPileUp;
    dplusSel(isSelColl);
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadronsDplusSelection, processDplusSelectionData, "Process Dplus Selection Data", false);

  void processDplusSelectionMcRec(SelCollisions::iterator const& collision,
                                  CandidatesDplusMcRec const& candidates)
  {
    bool isSelColl = true;
    bool isDplusFound = false;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    if (doSelDplusCollision) {
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yDplus(candidate)) >= yCandMax || candidate.pt() <= ptCandMin) {
          continue;
        }
        isDplusFound = true;
        break;
      }
    }
    if (useSel8) {
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
    }
    isSelColl = isDplusFound && isSel8 && isNosameBunchPileUp;
    dplusSel(isSelColl);
  }

  PROCESS_SWITCH(HfCorrelatorDplusHadronsDplusSelection, processDplusSelectionMcRec, "Process Dplus Selection MCRec", false);

  void processDplusSelectionMcGen(aod::McCollision const&,
                                  CandDplusMcGen const& mcParticles)
  {
    bool isDplusFound = false;
    for (const auto& particle1 : mcParticles) {
      if (std::abs(particle1.pdgCode()) != Pdg::kDPlus) {
        continue;
      }
      double const yD = RecoDecay::y(particle1.pVector(), MassDPlus);
      if (std::abs(yD) >= yCandMax || particle1.pt() <= ptCandMin) {
        continue;
      }
      isDplusFound = true;
      break;
    }
    dplusSel(isDplusFound);
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadronsDplusSelection, processDplusSelectionMcGen, "Process Dplus Selection MCGen", false);
};

/// Dplus-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDplusHadrons {
  Produces<aod::DplusHadronPair> entryDplusHadronPair;
  Produces<aod::DplusHadronRecoInfo> entryDplusHadronRecoInfo;
  Produces<aod::DplusHadronMlInfo> entryDplusHadronMlInfo;
  Produces<aod::DplusRecoInfo> entryDplusCandRecoInfo;
  Produces<aod::DplusHadronGenInfo> entryDplusHadronGenInfo;
  Produces<aod::DplusGenInfo> entryDplusCandGenInfo;
  Produces<aod::TrkRecInfoDplus> entryTrackRecoInfo;
  Produces<aod::Dplus> entryDplus;
  Produces<aod::Hadron> entryHadron;
  static constexpr std::size_t NDaughters{3u};

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for Dplus"}; // 7 corresponds to topo+PID cuts
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<bool> applyEfficiency{"applyEfficiency", true, "Flag for applying D-meson efficiency weights"};
  Configurable<bool> removeDaughters{"removeDaughters", true, "Flag for removing D-meson daughters from correlations"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 100., "max. track pT"};
  Configurable<float> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<float> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{o2::analysis::hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<float>> efficiencyD{"efficiencyD", {1., 1., 1., 1., 1., 1.}, "efficiency values for D+ meson"};

  HfHelper hfHelper;
  SliceCache cache;

  // Event Mixing for the Data Mode
  using SelCollisionsWithDplus = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::DmesonSelection>>;
  using SelCollisionsWithDplusMc = soa::Filtered<soa::Join<aod::McCollisions, aod::DmesonSelection, aod::MultsExtraMC>>; // collisionFilter applied
  using CandidatesDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  // Event Mixing for the MCRec Mode
  using CandidatesDplusMcRec = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using CandDplusMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>; // flagDplusFilter applied
  // Event Mixing for the MCGen Mode
  using McCollisionsSel = soa::Filtered<soa::Join<aod::McCollisions, aod::DmesonSelection>>;
  using McParticlesSel = soa::Filtered<aod::McParticles>;
  // Tracks used in Data and MC
  using TracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;                           // trackFilter applied
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels>>; // trackFilter applied

  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  // filter on selection of Dplus meson and decay channel Dplus->KPiPi
  Filter dplusFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi)) != static_cast<uint8_t>(0)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (nabs(aod::track::pt) > ptTrackMin) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);
  // Filter particlesFilter = nabs(aod::mcparticle::pdgCode) == 411 || ((aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsMultiplicityMc{"binsMultiplicityMc", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"}; // In MCGen multiplicity is defined by counting tracks
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 6000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsMassD{"binsMassD", {200, 1.7, 2.10}, "inv. mass (#pi^{+}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
  BinningType corrBinning{{binsZVtx, binsMultiplicity}, true};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMassD = {binsMassD, "inv. mass (pi^{+}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
    AxisSpec const axisEta = {binsEta, "#it{#eta}"};
    AxisSpec const axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec const axisMultiplicity = {binsMultiplicity, "Multiplicity"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec const axisPosZ = {binsZVtx, "PosZ"};
    AxisSpec const axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec const axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec const axisStatus = {15, 0.5, 15.5, "Selection status"};
    AxisSpec const axisRapidity = {100, -2, 2, "Rapidity"};

    registry.add("hPtCand", "Dplus,Hadron candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng0", "Dplus,Hadron candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng1", "Dplus,Hadron candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng2", "Dplus,Hadron candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}});
    registry.add("hSelectionStatus", "Dplus,Hadron candidates;selection status;entries", {HistType::kTH1F, {{10, -0.5, 9.5}}});
    registry.add("hEta", "Dplus,Hadron candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}});
    registry.add("hPhi", "Dplus,Hadron candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {axisPhi}});
    registry.add("hY", "Dplus,Hadron candidates;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hcountDplusHadronPerEvent", "Dplus,Hadron particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}});
    registry.add("hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hMultFT0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}});
    registry.add("hDplusBin", "Dplus selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}});
    registry.add("hTracksBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}});
    registry.add("hMassDplus_2D", "Dplus candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassDplusData", "Dplus candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{axisMassD}}});
    registry.add("hDplusPoolBin", "D+ candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
    // Histograms for MC Reco analysis
    registry.add("hSelectionStatusMCRec", "Dplus,Hadron candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}});
    registry.add("hMCEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}});
    registry.add("hPtProng0MCRec", "Dplus,Hadron candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng1MCRec", "Dplus,Hadron candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng2MCRec", "Dplus,Hadron candidates - MC reco;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}});
    registry.add("hMassDplusMcRec", "D+ candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{axisMassD}}});
    registry.add("hMassDplusVsPtMcRec", "D+ signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassDplusMcRecSig", "D+ signal candidates - MC reco;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassDplusMcRecBkg", "D+ background candidates - MC reco;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hPtCandMcRecSig", "D+,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcRecSigPrompt", "D+,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcRecSigNonPrompt", "D+,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcRecBkg", "D+,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hEtaMcRecSig", "D+,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcRecSig", "D+,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hYMCRecSig", "D+,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hEtaMcRecBkg", "D+,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcRecBkg", "D+,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hYMCRecBkg", "Dplus,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hFakeTracksMcRec", "Fake tracks - MC Rec", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hPtParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtD}}});
    registry.add("hPtPrimaryParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtD}}});
    registry.add("hPtVsMultiplicityMcRecPrompt", "Multiplicity FT0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultFT0M}}});
    registry.add("hPtVsMultiplicityMcRecNonPrompt", "Multiplicity FT0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultFT0M}}});
    // Histograms for MC Gen analysis
    registry.add("hcountDplustriggersMCGen", "D+ trigger particles - MC gen;;N of trigger Dplus", {HistType::kTH2F, {{1, -0.5, 0.5}, {axisPtD}}});
    registry.add("hPtCandMCGen", "Dplus,Hadron particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtD}});
    registry.add("hYMCGen", "Dplus,Hadron candidates - MC gen;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hPtCandMcGenPrompt", "D+,Hadron particles - MC Gen Prompt", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcGenNonPrompt", "D+,Hadron particles - MC Gen Non Prompt", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hEtaMcGen", "D+,Hadron particles - MC Gen", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcGen", "D+,Hadron particles - MC Gen", {HistType::kTH1F, {axisPhi}});
    registry.add("hMultFT0AMcGen", "D+,Hadron multiplicity FT0A - MC Gen", {HistType::kTH1F, {axisMultiplicity}});
    corrBinning = {{binsZVtx, binsMultiplicity}, true};
  }

  /// Dplus-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(SelCollisionsWithDplus::iterator const& collision,
                   TracksData const& tracks,
                   CandidatesDplusData const& candidates, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int gCollisionId = collision.globalIndex();
    int64_t timeStamp = bc.timestamp();

    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    if (candidates.size() > 0) {
      int nTracks = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) >= etaTrackMax || std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
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

      int cntDplus = 0;
      std::vector<float> outputMl = {-1., -1., -1.};
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yDplus(candidate)) >= yCandMax || candidate.pt() <= ptCandMin || candidate.pt() >= ptCandMax) {
          continue;
        }
        int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, candidate.pt());
        double efficiencyWeightD = 1.;
        if (applyEfficiency) {
          efficiencyWeightD = 1. / efficiencyD->at(effBinD);
        }
        // fill invariant mass plots and generic info from all Dplus candidates
        registry.fill(HIST("hMassDplus_2D"), hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), efficiencyWeightD);
        registry.fill(HIST("hMassDplusData"), hfHelper.invMassDplusToPiKPi(candidate), efficiencyWeightD);
        registry.fill(HIST("hPtCand"), candidate.pt());
        registry.fill(HIST("hPtProng0"), candidate.ptProng0());
        registry.fill(HIST("hPtProng1"), candidate.ptProng1());
        registry.fill(HIST("hPtProng2"), candidate.ptProng2());
        registry.fill(HIST("hEta"), candidate.eta());
        registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
        registry.fill(HIST("hY"), hfHelper.yDplus(candidate));
        registry.fill(HIST("hSelectionStatus"), candidate.isSelDplusToPiKPi());
        registry.fill(HIST("hDplusBin"), poolBin);
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
        }
        entryDplusCandRecoInfo(hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]); // 0: BkgBDTScore, 1:PromptBDTScore, 2: FDScore
        entryDplus(candidate.phi(), candidate.eta(), candidate.pt(), hfHelper.invMassDplusToPiKPi(candidate), poolBin, gCollisionId, timeStamp);

        // Dplus-Hadron correlation dedicated section
        // if the candidate is a Dplus, search for Hadrons and evaluate correlations
        for (const auto& track : tracks) {
          // apply track selection
          if (!track.isGlobalTrackWoDCA()) {
            continue;
          }
          // Removing Dplus daughters by checking track indices
          if (removeDaughters) {
            if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
              continue;
            }
          }
          entryDplusHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                               track.eta() - candidate.eta(),
                               candidate.pt(),
                               track.pt(), poolBin);
          entryDplusHadronRecoInfo(hfHelper.invMassDplusToPiKPi(candidate), false);
          entryDplusHadronGenInfo(false, false, 0);
          entryDplusHadronMlInfo(outputMl[0], outputMl[1], outputMl[2]);
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
          if (cntDplus == 0) {
            entryHadron(track.phi(), track.eta(), track.pt(), poolBin, gCollisionId, timeStamp);
            registry.fill(HIST("hTracksBin"), poolBin);
          }
        } // Hadron Tracks loop
        cntDplus++;
      } // end outer Dplus loop
      registry.fill(HIST("hZvtx"), collision.posZ());
      registry.fill(HIST("hMultFT0M"), collision.multFT0M());
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processData, "Process data", false);

  /// Dplus-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(SelCollisionsWithDplus::iterator const& collision,
                    TracksWithMc const& tracks,
                    CandidatesDplusMcRec const& candidates,
                    aod::McParticles const& mcParticles)
  {
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    if (candidates.size() > 0) {
      int nTracks = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) >= etaTrackMax || std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
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

      float const multiplicityFT0M = collision.multFT0M();

      // MC reco level
      bool isDplusPrompt = false;
      bool isDplusNonPrompt = false;
      bool isDplusSignal = false;
      for (const auto& candidate : candidates) {
        // rapidity and pT selections
        if (std::abs(hfHelper.yDplus(candidate)) >= yCandMax || candidate.pt() <= ptCandMin || candidate.pt() >= ptCandMax) {
          continue;
        }
        // efficiency weight determination
        int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, candidate.pt());
        double efficiencyWeightD = 1.;
        if (applyEfficiency) {
          efficiencyWeightD = 1. / efficiencyD->at(effBinD);
        }
        // Dplus flag
        isDplusSignal = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi;
        // prompt and non-prompt division
        isDplusPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        isDplusNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;

        std::vector<float> outputMl = {-1., -1., -1.};

        // fill invariant mass plots from Dplus signal and background candidates
        registry.fill(HIST("hMassDplusMcRec"), hfHelper.invMassDplusToPiKPi(candidate), efficiencyWeightD);
        registry.fill(HIST("hDplusBin"), poolBin);

        if (isDplusSignal) {
          // fill per-candidate distributions from Dplus true candidates
          registry.fill(HIST("hPtProng0MCRec"), candidate.ptProng0());
          registry.fill(HIST("hPtProng1MCRec"), candidate.ptProng1());
          registry.fill(HIST("hPtProng2MCRec"), candidate.ptProng2());
          registry.fill(HIST("hMassDplusVsPtMcRec"), hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusMCRec"), candidate.isSelDplusToPiKPi());
          registry.fill(HIST("hPtCandMcRecSig"), candidate.pt());
          registry.fill(HIST("hEtaMcRecSig"), candidate.eta());
          registry.fill(HIST("hYMCRecSig"), hfHelper.yDplus(candidate));
          registry.fill(HIST("hPhiMcRecSig"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));

          // prompt and non-prompt division
          if (isDplusPrompt) {
            registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), multiplicityFT0M);
          } else if (isDplusNonPrompt) {
            registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), multiplicityFT0M);
          }
          // Storing ML scores for signal reco candidates and charm/beauty origin
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
          }
          registry.fill(HIST("hMassDplusMcRecSig"), hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), efficiencyWeightD);
          entryDplusCandRecoInfo(hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
          entryDplusCandGenInfo(isDplusPrompt);
        } else {
          registry.fill(HIST("hPtCandMcRecBkg"), candidate.pt());
          registry.fill(HIST("hEtaMcRecBkg"), candidate.eta());
          registry.fill(HIST("hPhiMcRecBkg"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
          registry.fill(HIST("hMassDplusMcRecBkg"), hfHelper.invMassDplusToPiKPi(candidate), candidate.pt(), efficiencyWeightD);
        }

        // Dplus-Hadron correlation dedicated section
        // if the candidate is selected as Dplus, search for Hadron and evaluate correlations
        for (const auto& track : tracks) {
          bool isPhysicalPrimary = false;
          int trackOrigin = -1;
          // apply track selection
          if (!track.isGlobalTrackWoDCA()) {
            continue;
          }
          // Removing Dplus daughters by checking track indices
          if (removeDaughters) {
            if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
              continue;
            }
          }
          entryDplusHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                               track.eta() - candidate.eta(),
                               candidate.pt(),
                               track.pt(), poolBin);
          entryDplusHadronRecoInfo(hfHelper.invMassDplusToPiKPi(candidate), isDplusSignal);
          entryDplusHadronMlInfo(outputMl[0], outputMl[1], outputMl[2]);
          if (track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            isPhysicalPrimary = mcParticle.isPhysicalPrimary();
            trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
            entryDplusHadronGenInfo(isDplusPrompt, isPhysicalPrimary, trackOrigin);
          } else {
            entryDplusHadronGenInfo(isDplusPrompt, false, 0);
            registry.fill(HIST("hFakeTracksMcRec"), track.pt());
          }
          // for secondary particle fraction estimation
          registry.fill(HIST("hPtParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
          if (isPhysicalPrimary) {
            registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
          }
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
        } // end inner loop (Tracks)

      } // end outer Dplus loop
      registry.fill(HIST("hZvtx"), collision.posZ());
      registry.fill(HIST("hMultFT0M"), collision.multFT0M());
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcRec, "Process MC Reco mode", true);

  /// Dplus-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(SelCollisionsWithDplusMc::iterator const& mcCollision,
                    CandDplusMcGen const& mcParticles)
  {
    int counterDplusHadron = 0;
    registry.fill(HIST("hMCEvtCount"), 0);

    BinningTypeMcGen const corrBinningMcGen{{binsZVtx, binsMultiplicityMc}, true};
    int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), mcCollision.multMCFT0A()));
    registry.fill(HIST("hMultFT0AMcGen"), mcCollision.multMCFT0A());

    bool isDplusPrompt = false;
    bool isDplusNonPrompt = false;
    // MC gen level
    for (const auto& particle1 : mcParticles) {
      // check if the particle is Dplus  (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != Pdg::kDPlus) {
        continue;
      }
      if (std::abs(particle1.flagMcMatchGen()) != hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
        continue;
      }
      double const yD = RecoDecay::y(particle1.pVector(), MassDPlus);
      if (std::abs(yD) >= yCandMax || particle1.pt() <= ptCandMin) {
        continue;
      }
      registry.fill(HIST("hDplusBin"), poolBin);
      registry.fill(HIST("hPtCandMCGen"), particle1.pt());
      registry.fill(HIST("hEtaMcGen"), particle1.eta());
      registry.fill(HIST("hPhiMcGen"), RecoDecay::constrainAngle(particle1.phi(), -PIHalf));
      registry.fill(HIST("hYMCGen"), yD);

      // prompt and non-prompt division
      isDplusPrompt = particle1.originMcGen() == RecoDecay::OriginType::Prompt;
      isDplusNonPrompt = particle1.originMcGen() == RecoDecay::OriginType::NonPrompt;
      if (isDplusPrompt) {
        registry.fill(HIST("hPtCandMcGenPrompt"), particle1.pt());
      } else if (isDplusNonPrompt) {
        registry.fill(HIST("hPtCandMcGenNonPrompt"), particle1.pt());
      }

      // prompt and non-prompt division
      std::vector<int> listDaughters{};
      std::array<int, NDaughters> const arrDaughDplusPDG = {+kPiPlus, -kKPlus, kPiPlus};
      std::array<int, NDaughters> prongsId{};
      listDaughters.clear();
      RecoDecay::getDaughters(particle1, &listDaughters, arrDaughDplusPDG, 2);
      int counterDaughters = 0;
      if (listDaughters.size() == NDaughters) {
        for (const auto& dauIdx : listDaughters) {
          auto daughI = mcParticles.rawIteratorAt(dauIdx - mcParticles.offset());
          counterDaughters += 1;
          prongsId[counterDaughters - 1] = daughI.globalIndex();
        }
      }
      counterDplusHadron++;
      // Dplus Hadron correlation dedicated section
      // if it's a Dplus particle, search for Hadron and evaluate correlations
      registry.fill(HIST("hcountDplustriggersMCGen"), 0, particle1.pt()); // to count trigger Dplus for normalisation)
      for (const auto& particleAssoc : mcParticles) {
        if (std::abs(particleAssoc.eta()) > etaTrackMax || particleAssoc.pt() < ptTrackMin || particleAssoc.pt() > ptTrackMax) {
          continue;
        }
        if (removeDaughters) {
          if (particleAssoc.globalIndex() == prongsId[0] || particleAssoc.globalIndex() == prongsId[1] || particleAssoc.globalIndex() == prongsId[2]) {
            continue;
          }
        }
        if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
          continue;
        }
        if (!particleAssoc.isPhysicalPrimary()) {
          continue;
        }

        int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
        registry.fill(HIST("hPtParticleAssocMcGen"), particleAssoc.pt());
        entryDplusHadronPair(getDeltaPhi(particleAssoc.phi(), particle1.phi()),
                             particleAssoc.eta() - particle1.eta(),
                             particle1.pt(),
                             particleAssoc.pt(),
                             poolBin);
        entryDplusHadronRecoInfo(MassDPlus, true);
        entryDplusHadronGenInfo(isDplusPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
      } // end associated loop
    } // end trigger
    registry.fill(HIST("hcountDplusHadronPerEvent"), counterDplusHadron);
    registry.fill(HIST("hZvtx"), mcCollision.posZ());
    // registry.fill(HIST("hMultiplicity"), getTracksSize(mcCollision));
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcGen, "Process MC Gen mode", false);

  void processDataMixedEvent(SelCollisionsWithDplus const& collisions,
                             CandidatesDplusData const& candidates,
                             TracksData const& tracks)
  {
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelCollisionsWithDplus, CandidatesDplusData, TracksData, BinningType> const pairData{corrBinning, 5, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      // LOGF(info, "Mixed event collisions: Index = (%d, %d), tracks Size: (%d, %d), Z Vertex: (%f, %f), Pool Bin: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), c1.posZ(), c2.posZ(), corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M())),corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()))); // For debug
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      for (const auto& [trigDplus, assocParticle] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!assocParticle.isGlobalTrackWoDCA() || std::abs(hfHelper.yDplus(trigDplus)) >= yCandMax) {
          continue;
        }
        entryDplusHadronPair(getDeltaPhi(trigDplus.phi(), assocParticle.phi()), trigDplus.eta() - assocParticle.eta(), trigDplus.pt(), assocParticle.pt(), poolBin);
        entryDplusHadronRecoInfo(hfHelper.invMassDplusToPiKPi(trigDplus), 0);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processDataMixedEvent, "Process Mixed Event Data", false);

  void processMcRecMixedEvent(SelCollisionsWithDplus const& collisions,
                              CandidatesDplusMcRec const& candidates,
                              TracksWithMc const& tracks,
                              aod::McParticles const& mcParticles)
  {
    BinningType const corrBinning{{binsZVtx, binsMultiplicityMc}, true};
    for (const auto& candidate : candidates) {
      if (std::abs(hfHelper.yDplus(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      // Dplus flag
      bool const isDplusSignal = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi;
      // prompt and non-prompt division
      bool const isDplusPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
      bool const isDplusNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
      if (isDplusSignal) {
        if (isDplusPrompt) {
          registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
          registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), 0);
        } else if (isDplusNonPrompt) {
          registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
          registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), 0);
        }
      } else {
        registry.fill(HIST("hPtCandMcRecBkg"), candidate.pt());
        registry.fill(HIST("hEtaMcRecBkg"), candidate.eta());
        registry.fill(HIST("hPhiMcRecBkg"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
      }
    }
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelCollisionsWithDplus, CandidatesDplusMcRec, TracksWithMc, BinningType> const pairMcRec{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      int const poolBinDplus = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hMultFT0M"), c1.multFT0M());
      registry.fill(HIST("hZVtx"), c1.posZ());
      registry.fill(HIST("hTracksPoolBin"), poolBin);     // note that the selections here are not yet applied
      registry.fill(HIST("hDplusPoolBin"), poolBinDplus); // note that the selections here are not yet applied
      for (const auto& [candidate, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (std::abs(hfHelper.yDplus(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
          continue;
        }
        if (!pAssoc.isGlobalTrackWoDCA()) {
          continue;
        }
        std::vector<float> outputMl = {-1., -1., -1.};
        bool isPhysicalPrimary = false;
        int trackOrigin = -1;
        bool isDplusSignal = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi;
        // prompt and non-prompt division
        bool isDplusPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        if (pAssoc.has_mcParticle()) {
          auto mcParticle = pAssoc.template mcParticle_as<aod::McParticles>();
          isPhysicalPrimary = mcParticle.isPhysicalPrimary();
          trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
        } else {
          registry.fill(HIST("hFakeTracksMcRec"), pAssoc.pt());
        }
        entryDplusHadronPair(getDeltaPhi(pAssoc.phi(), candidate.phi()),
                             pAssoc.eta() - candidate.eta(),
                             candidate.pt(),
                             pAssoc.pt(),
                             poolBin);
        entryDplusHadronRecoInfo(hfHelper.invMassDplusToPiKPi(candidate), isDplusSignal);
        entryDplusHadronGenInfo(isDplusPrompt, isPhysicalPrimary, trackOrigin);
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
        }
        entryDplusHadronMlInfo(outputMl[0], outputMl[1], outputMl[2]);
        entryTrackRecoInfo(pAssoc.dcaXY(), pAssoc.dcaZ(), pAssoc.tpcNClsCrossedRows());
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcRecMixedEvent, "Process Mixed Event MCRec", false);

  void processMcGenMixedEvent(SelCollisionsWithDplusMc const& collisions,
                              CandDplusMcGen const& mcParticles)
  {
    BinningTypeMcGen const corrBinningMcGen{{binsZVtx, binsMultiplicityMc}, true};
    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<SelCollisionsWithDplusMc, CandDplusMcGen, CandDplusMcGen, BinningTypeMcGen> const pairMcGen{corrBinningMcGen, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      int poolBin = corrBinningMcGen.getBin(std::make_tuple(c1.posZ(), c1.multMCFT0A()));
      for (const auto& [candidate, particleAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (std::abs(candidate.pdgCode()) != Pdg::kDPlus) {
          continue;
        }
        double const yD = RecoDecay::y(candidate.pVector(), MassDPlus);
        if (std::abs(yD) > yCandGenMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
          continue;
        }
        if (std::abs(particleAssoc.eta()) > etaTrackMax || particleAssoc.pt() < ptTrackMin || particleAssoc.pt() > ptTrackMax) {
          continue;
        }
        if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
          continue;
        }
        if (!particleAssoc.isPhysicalPrimary()) {
          continue;
        }
        int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
        bool isDplusPrompt = candidate.originMcGen() == RecoDecay::OriginType::Prompt;
        entryDplusHadronPair(getDeltaPhi(particleAssoc.phi(), candidate.phi()),
                             particleAssoc.eta() - candidate.eta(),
                             candidate.pt(),
                             particleAssoc.pt(),
                             poolBin);
        entryDplusHadronRecoInfo(MassDPlus, true);
        entryDplusHadronGenInfo(isDplusPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcGenMixedEvent, "Process Mixed Event MCGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDplusHadronsDplusSelection>(cfgc),
                      adaptAnalysisTask<HfCorrelatorDplusHadrons>(cfgc)};
}
