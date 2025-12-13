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

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/Utils/utilsCorrelations.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
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
#include <Framework/DataTypes.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::hf_correlations;

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiHadron, double phiD)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

// Types
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
using BinningTypeMcGen = ColumnBinningPolicy<aod::mccollision::PosZ, o2::aod::mult::MultMCFT0A>;
using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::DmesonSelection, aod::CentFT0Ms>>;
using SelectedTracks = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;
using SelectedCandidatesData = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
using SelectedCandidatesDataMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
using SelectedTracksMcRec = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels>>;
using SelectedCandidatesMcRec = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;
using SelectedCandidatesMcRecMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0, aod::HfCand2ProngMcRec>>;
using SelectedCollisionsMcGen = soa::Filtered<soa::Join<aod::McCollisions, aod::DmesonSelection, aod::MultsExtraMC>>;
using SelectedParticlesMcGen = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

// Code to select collisions with at least one D0
struct HfCorrelatorD0HadronsSelection {
  Produces<aod::DmesonSelection> collisionsWithSelD0;

  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<float> yCandMax{"yCandMax", 4.0, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
  Configurable<float> centMin{"centMin", 0., "Minimum Centrality"};
  Configurable<float> centMax{"centMax", 100., "Maximum Centrality"};
  Configurable<bool> useCentrality{"useCentrality", false, "Flag for centrality dependent analyses"};

  SliceCache cache;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;

  Preslice<aod::HfCand2Prong> perCol = aod::hf_cand::collisionId;

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> selectedD0candidatesMc = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  void processD0SelectionData(SelCollisions::iterator const& collision,
                              soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&)
  {
    bool isSelColl = true;
    bool isD0Found = true;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    bool isCentInRange = false;
    if (selectedD0Candidates.size() > 0) {
      auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

      for (const auto& candidate : selectedD0CandidatesGrouped) {
        // check decay channel flag for candidate
        if (!TESTBIT(candidate.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          isD0Found = false;
          continue;
        }
        if (std::abs(HfHelper::yD0(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          isD0Found = false;
          continue;
        }
        isD0Found = true;
        break;
      }
    }
    float cent = 0.;
    if (useCentrality) {
      cent = collision.centFT0M();
    }
    if (useSel8) {
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
    }
    isCentInRange = (cent >= centMin && cent <= centMax);
    isSelColl = isD0Found && isSel8 && isNosameBunchPileUp && isCentInRange;
    collisionsWithSelD0(isSelColl);
  }
  PROCESS_SWITCH(HfCorrelatorD0HadronsSelection, processD0SelectionData, "Process D0 Selection Data", false);

  void processD0SelectionMcRec(SelCollisions::iterator const& collision,
                               soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const&)
  {
    bool isSelColl = true;
    bool isD0Found = true;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    if (selectedD0candidatesMc.size() > 0) {
      auto selectedD0CandidatesGroupedMc = selectedD0candidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (const auto& candidate : selectedD0CandidatesGroupedMc) {
        // check decay channel flag for candidate
        if (!TESTBIT(candidate.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          isD0Found = false;
          continue;
        }
        if (std::abs(HfHelper::yD0(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          isD0Found = false;
          continue;
        }
        isD0Found = true;
        break;
      }
    }
    if (useSel8) {
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
    }
    isSelColl = isD0Found && isSel8 && isNosameBunchPileUp;
    collisionsWithSelD0(isSelColl);
  }
  PROCESS_SWITCH(HfCorrelatorD0HadronsSelection, processD0SelectionMcRec, "Process D0 Selection MCRec", true);

  void processD0SelectionMcGen(aod::McCollision const&,
                               aod::McParticles const& mcParticles)
  {
    bool isD0Found = false;
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != Pdg::kD0) {
        continue;
      }
      double const yD = RecoDecay::y(particle.pVector(), MassD0);
      if (std::abs(yD) > yCandMax || particle.pt() < ptCandMin) {
        continue;
      }
      isD0Found = true;
      break;
    }
    collisionsWithSelD0(isD0Found);
  }
  PROCESS_SWITCH(HfCorrelatorD0HadronsSelection, processD0SelectionMcGen, "Process D0 Selection MCGen", false);
};

/// D0-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorD0Hadrons {
  Produces<aod::D0HadronPair> entryD0HadronPair;
  Produces<aod::D0HadronRecoInfo> entryD0HadronRecoInfo;
  Produces<aod::D0HadronGenInfo> entryD0HadronGenInfo;
  Produces<aod::D0HadronMlInfo> entryD0HadronMlInfo;
  Produces<aod::D0CandRecoInfo> entryD0CandRecoInfo;
  Produces<aod::D0CandGenInfo> entryD0CandGenInfo;
  Produces<aod::D0TrackRecoInfo> entryTrackRecoInfo;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCAxy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCAz of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 99., "max. track pT"};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", {1., 1., 1., 1., 1., 1.}, "Efficiency values for D0 meson"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<float> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<float> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<float> ptSoftPionMax{"ptSoftPionMax", 3.f * 800.f * std::pow(10.f, -6.f), "max. pT cut for soft pion identification"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<bool> correlateD0WithLeadingParticle{"correlateD0WithLeadingParticle", false, "Switch for correlation of D0 mesons with leading particle only"};
  Configurable<bool> storeAutoCorrelationFlag{"storeAutoCorrelationFlag", false, "Store flag that indicates if the track is paired to its D-meson mother instead of skipping it"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<bool> useCentrality{"useCentrality", false, "Flag for centrality dependent analyses"};

  int leadingIndex = 0;
  double massD0{0.};
  double massPi{0.};
  double massK{0.};
  double softPiMass = 0.14543; // pion mass + Q-value of the D*->D0pi decay

  SliceCache cache;

  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter trackFilter = requireGlobalTrackWoDCAInFilter() && (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);
  Filter d0Filter = (aod::hf_sel_candidate_d0::isSelD0 >= 1) || (aod::hf_sel_candidate_d0::isSelD0bar >= 1);
  Filter collisionFilterGen = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter particlesFilter = nabs(aod::mcparticle::pdgCode) == static_cast<int>(Pdg::kD0) || ((aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);

  Preslice<aod::HfCand2Prong> perCol = aod::hf_cand::collisionId;
  Preslice<aod::Tracks> perCollisionID = aod::track::collisionId;
  Preslice<aod::McParticles> perTrueCollision = o2::aod::mcparticle::mcCollisionId;

  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 10000.0f}, "event multiplicity pools (FT0M)"};
  ConfigurableAxis multPoolBinsMcGen{"multPoolBinsMcGen", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"}; // In MCGen multiplicity is defined by counting tracks
  ConfigurableAxis binsMassD{"binsMassD", {200, 1.3848, 2.3848}, "inv. mass (#pi K) (GeV/#it{c}^{2});entries"};
  ConfigurableAxis binsEta{"binsEta", {100, -5., 5.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {10000, 0., 10000.}, "Multiplicity"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {10000, 0., 10000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsCentFt0m{"binsCentFt0m", {100, 0., 100.}, "Centrality percentile (FT0M)"};

  BinningType corrBinning{{zPoolBins, multPoolBins}, true};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    massD0 = MassD0;
    massPi = MassPiPlus;
    massK = MassKPlus;

    AxisSpec axisMassD = {binsMassD, "inv. mass (#pi K) (GeV/#it{c}^{2})"};
    AxisSpec const axisEta = {binsEta, "#it{#eta}"};
    AxisSpec const axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec const axisRapidity = {100, -5., 5., "Rapidity"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec const axisMultiplicity = {binsMultiplicity, "Multiplicity"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec const axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec const axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec const axisStatus = {4, -0.5, 3.5, "Selection status"};
    AxisSpec const axisSignalStatus = {200, 0., 200., "Signal status"};
    AxisSpec axisEvtCount = {1, -0.5, 0.5};
    AxisSpec const axisTrkCount = {5, 0., 5.};
    AxisSpec axisBdtScoreBkg = {100, 0., 1., "Bdt score background"};
    AxisSpec axisBdtScorePrompt = {100, 0., 1., "Bdt score prompt"};
    AxisSpec axisBdtScoreNonPrompt = {100, 0., 1., "Bdt score Nonprompt"};
    AxisSpec axisOrigin = {10, 0., 10., "Candidate origin"};
    AxisSpec axisCent = {binsCentFt0m, "Centrality"};

    // Histograms for Data
    registry.add("hPtCand", "D0, D0bar candidates", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng0", "D0, D0bar candidates prong 0", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng1", "D0, D0bar candidates prong 1", {HistType::kTH1F, {axisPtD}});
    registry.add("hSelectionStatus", "D0, D0bar candidates selection status", {HistType::kTH1F, {axisStatus}});
    registry.add("hEta", "D0,D0bar candidates", {HistType::kTH1F, {axisEta}});
    registry.add("hPhi", "D0,D0bar candidates", {HistType::kTH1F, {axisPhi}});
    registry.add("hY", "D0,D0bar candidates", {HistType::kTH1F, {axisRapidity}});
    registry.add("hCentFT0M", "Centrality FT0M;centrality;entries", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("hMass", "D0, D0bar candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMass1D", "D0, D0bar candidates mass", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassD01D", "D0 candidates mass", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassD0bar1D", "D0bar candidates mass", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassD0VsPtVsCent", "D0 candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCent}}});
    registry.add("hMLScoresVsMassVsPtVsEtaVsOrigin", "D0, D0bar candidates BkgVspromptVsNonPromptVsMassVsPtVsEtaVsOrigin", {HistType::kTHnSparseD, {{axisBdtScoreBkg}, {axisBdtScorePrompt}, {axisBdtScoreNonPrompt}, {axisMassD}, {axisPtD}, {axisEta}, {axisOrigin}}});
    // Histograms for MC Reco
    registry.add("hPtCandRec", "D0, D0bar candidates - MC reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng0Rec", "D0, D0bar candidates prong 0 - MC reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtProng1Rec", "D0, D0bar candidates prong 1 - MC reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hSelectionStatusRec", "D0, D0bar candidates selection status - MC reco", {HistType::kTH1F, {axisStatus}});
    registry.add("hSignalStatusMERec", "Signal Status - MC reco ME", {HistType::kTH1F, {axisSignalStatus}});
    registry.add("hEtaRec", "D0,D0bar candidates - MC reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiRec", "D0,D0bar candidates - MC reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hYRec", "D0,D0bar candidates - MC reco", {HistType::kTH1F, {axisRapidity}});
    registry.add("hMassD0RecSig", "D0 signal candidates massVsPt - MC reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassD0RecRef", "D0 reflection candidates massVsPt - MC reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassD0RecBg", "D0 background candidates massVsPt - MC reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassD0barRecSig", "D0bar signal candidates massVsPt - MC reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassD0barRecRef", "D0bar reflection candidates massVsPt - MC reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassD0barRecBg", "D0bar background candidates massVsPt - MC reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hPtCandRecSigPrompt", "D0,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtVsMLScoresVsEtaRecSigPrompt", "Prompt D0-D0bar signal candidates MLVsPtVsEta - MC reco", {HistType::kTHnSparseD, {{axisBdtScoreBkg}, {axisBdtScorePrompt}, {axisBdtScoreNonPrompt}, {axisPtD}, {axisEta}}});
    registry.add("hPtCandRecSigNonPrompt", "D0,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtVsMLScoresVsEtaRecSigNonPrompt", "NonPrompt D0-D0bar signal candidates MLVsPtVsEta - MC reco", {HistType::kTHnSparseD, {{axisBdtScoreBkg}, {axisBdtScorePrompt}, {axisBdtScoreNonPrompt}, {axisPtD}, {axisEta}}});
    registry.add("hPtVsMultiplicityRecPrompt", "Multiplicity FT0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultFT0M}}});
    registry.add("hPtVsMultiplicityRecNonPrompt", "Multiplicity FT0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultFT0M}}});
    registry.add("hPtParticleAssocVsCandRec", "Associated Particle - MC reco", {HistType::kTH2F, {{axisPtHadron}, {axisPtD}}});
    registry.add("hPtPrimaryParticleAssocVsCandRec", "Associated Particle - MC reco", {HistType::kTH2F, {{axisPtHadron}, {axisPtD}}});
    // Histograms for MC Gen
    registry.add("hEvtCountGen", "Event counter - MC gen", {HistType::kTH1F, {axisEvtCount}});
    registry.add("hPtCandGen", "D0, D0bar candidates - MC gen", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandGenPrompt", "D0, D0bar candidates - MC gen prompt", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtVsEtaCandGenSigPrompt", "D0,Hadron candidates PtvsEta - MC Gen prompt", {HistType::kTH2F, {{axisPtD}, {axisEta}}});
    registry.add("hPtCandGenNonPrompt", "D0, D0bar candidates - MC gen non prompt", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtVsEtaCandGenSigNonPrompt", "D0,Hadron candidates PtvsEta - MC Gen non-prompt", {HistType::kTH2F, {{axisPtD}, {axisEta}}});
    registry.add("hEtaGen", "D0,D0bar candidates - MC gen", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiGen", "D0,D0bar candidates - MC gen", {HistType::kTH1F, {axisPhi}});
    registry.add("hYGen", "D0,D0bar candidates - MC gen", {HistType::kTH1F, {axisRapidity}});
    registry.add("hCountD0TriggersGen", "D0 trigger particles - MC gen;;N of trigger D0", {HistType::kTH2F, {{axisEvtCount}, {axisPtD}}});
    // Common histograms
    registry.add("hTrackCounter", "Track counter", {HistType::kTH1F, {axisTrkCount}});
    registry.get<TH1>(HIST("hTrackCounter"))->GetXaxis()->SetBinLabel(1, "all");
    registry.get<TH1>(HIST("hTrackCounter"))->GetXaxis()->SetBinLabel(2, "before softpi");
    registry.get<TH1>(HIST("hTrackCounter"))->GetXaxis()->SetBinLabel(3, "after softpi");
    registry.get<TH1>(HIST("hTrackCounter"))->GetXaxis()->SetBinLabel(4, "with leading particles");
    registry.get<TH1>(HIST("hTrackCounter"))->GetXaxis()->SetBinLabel(5, "fake tracks");
    registry.add("hZvtx", "z vertex", {HistType::kTH1F, {axisPosZ}});
    registry.add("hMultFT0M", "Multiplicity FT0M", {HistType::kTH1F, {axisMultFT0M}});
    registry.add("hCollisionPoolBin", "collision pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hD0PoolBin", "D0 selected in pool Bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
  }

  // =======  Process starts for Data, Same event ============

  /// D0-h correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(SelectedCollisions::iterator const& collision,
                   SelectedTracks const& tracks,
                   SelectedCandidatesDataMl const& candidates)
  {
    // find leading particle
    if (correlateD0WithLeadingParticle) {
      leadingIndex = findLeadingParticle(tracks, etaTrackMax.value);
    }
    float cent = 0.;
    if (useCentrality) {
      cent = collision.centFT0M();
    }

    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());

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
    std::vector<float> outputMlD0 = {-1., -1., -1.};
    std::vector<float> outputMlD0bar = {-1., -1., -1.};

    for (const auto& candidate : candidates) {
      if (std::abs(HfHelper::yD0(candidate)) >= yCandMax || candidate.pt() <= ptCandMin || candidate.pt() >= ptTrackMax) {
        continue;
      }
      // check decay channel flag for candidate
      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }

      // ========================== Define parameters for soft pion removal ================================
      auto ePiK = RecoDecay::e(candidate.pVectorProng0(), massPi) + RecoDecay::e(candidate.pVectorProng1(), massK);
      auto eKPi = RecoDecay::e(candidate.pVectorProng0(), massK) + RecoDecay::e(candidate.pVectorProng1(), massPi);

      // ========================== trigger efficiency ================================
      double efficiencyWeight = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()));
      }

      // Invariant mass of D0 and D0bar
      const auto invMassD0 = HfHelper::invMassD0ToPiK(candidate);
      const auto invMassD0bar = HfHelper::invMassD0barToKPi(candidate);

      // ========================== Fill mass histo  ================================
      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), invMassD0, candidate.pt(), efficiencyWeight);
        registry.fill(HIST("hMass1D"), invMassD0, efficiencyWeight);
        registry.fill(HIST("hMassD01D"), invMassD0, efficiencyWeight);
        registry.fill(HIST("hMassD0VsPtVsCent"), invMassD0, candidate.pt(), cent, efficiencyWeight);
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMlD0[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
        }
        registry.fill(HIST("hMLScoresVsMassVsPtVsEtaVsOrigin"), outputMlD0[0], outputMlD0[1], outputMlD0[2], invMassD0, candidate.pt(), candidate.eta(), (candidate.isSelD0bar() != 0) ? o2::aod::hf_correlation_d0_hadron::D0D0barBoth : o2::aod::hf_correlation_d0_hadron::D0Only);
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), invMassD0bar, candidate.pt(), efficiencyWeight);
        registry.fill(HIST("hMass1D"), invMassD0bar, efficiencyWeight);
        registry.fill(HIST("hMassD0bar1D"), invMassD0bar, efficiencyWeight);
        registry.fill(HIST("hMassD0VsPtVsCent"), invMassD0bar, candidate.pt(), cent, efficiencyWeight);
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMlD0bar[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
        }
        registry.fill(HIST("hMLScoresVsMassVsPtVsEtaVsOrigin"), outputMlD0bar[0], outputMlD0bar[1], outputMlD0bar[2], invMassD0bar, candidate.pt(), candidate.eta(), (candidate.isSelD0() != 0) ? o2::aod::hf_correlation_d0_hadron::D0D0barBoth : o2::aod::hf_correlation_d0_hadron::D0barOnly);
      }
      entryD0CandRecoInfo(invMassD0, invMassD0bar, candidate.pt(), outputMlD0[0], outputMlD0[2], outputMlD0bar[0], outputMlD0bar[2]);

      // ========================== Fill general histos ================================
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hEta"), candidate.eta());
      registry.fill(HIST("hPhi"), candidate.phi());
      registry.fill(HIST("hY"), HfHelper::yD0(candidate));
      registry.fill(HIST("hSelectionStatus"), candidate.isSelD0bar() + (candidate.isSelD0() * 2));
      registry.fill(HIST("hD0PoolBin"), poolBin);

      // ============ D-h correlation dedicated section ==================================

      // ========================== track loop starts here ================================
      for (const auto& track : tracks) {
        registry.fill(HIST("hTrackCounter"), 0); // fill total no. of tracks
        // Remove D0 daughters by checking track indices
        bool correlationStatus = false;
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex())) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }
        registry.fill(HIST("hTrackCounter"), 1); // fill no. of tracks before soft pion removal

        // ========== soft pion removal ===================================================
        double invMassDstar1 = 0., invMassDstar2 = 0.;
        auto pSum2 = RecoDecay::p2(candidate.pVector(), track.pVector());
        auto ePion = track.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (candidate.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - invMassD0) - softPiMass) < ptSoftPionMax) {
            continue;
          }
        }

        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - invMassD0bar) - softPiMass) < ptSoftPionMax) {
            continue;
          }
        }
        registry.fill(HIST("hTrackCounter"), 2); // fill no. of tracks after soft pion removal

        int signalStatus = 0;
        if (candidate.isSelD0() >= selectionFlagD0) {
          signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0Only;
        }
        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0barOnly;
        }

        if (correlateD0WithLeadingParticle) {
          if (track.globalIndex() != leadingIndex) {
            continue;
          }
          registry.fill(HIST("hTrackCounter"), 3); // fill no. of tracks  have leading particle
        }
        entryD0HadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                          track.eta() - candidate.eta(),
                          candidate.pt(),
                          track.pt(),
                          poolBin,
                          correlationStatus,
                          cent);
        entryD0HadronRecoInfo(invMassD0, invMassD0bar, signalStatus);
        entryD0HadronGenInfo(false, false, 0);
        entryD0HadronMlInfo(outputMlD0[0], outputMlD0[1], outputMlD0[2], outputMlD0bar[0], outputMlD0bar[1], outputMlD0bar[2]);
        entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
        registry.fill(HIST("hCentFT0M"), cent);

      } // end inner loop (tracks)

    } // end outer loop
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processData, "Process data", false);

  // ================  Process starts for MCRec, same event ========================

  void processMcRec(SelectedCollisions::iterator const& collision,
                    SelectedTracksMcRec const& tracks,
                    SelectedCandidatesMcRecMl const& candidates,
                    aod::McParticles const& mcParticles)
  {
    // find leading particle
    if (correlateD0WithLeadingParticle) {
      leadingIndex = findLeadingParticle(tracks, etaTrackMax.value);
    }
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());

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

    // MC reco level
    bool flagD0 = false;
    bool flagD0bar = false;
    std::vector<float> outputMlD0 = {-1., -1., -1.};
    std::vector<float> outputMlD0bar = {-1., -1., -1.};

    for (const auto& candidate : candidates) {
      bool isD0Prompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
      bool isD0NonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
      // check decay channel flag for candidate
      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (std::abs(HfHelper::yD0(candidate)) >= yCandMax || candidate.pt() <= ptCandMin || candidate.pt() >= ptTrackMax) {
        continue;
      }

      registry.fill(HIST("hD0PoolBin"), poolBin);

      double efficiencyWeight = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()));
      }

      const auto invMassD0 = HfHelper::invMassD0ToPiK(candidate);
      const auto invMassD0bar = HfHelper::invMassD0barToKPi(candidate);

      if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
        // fill per-candidate distributions from D0/D0bar true candidates
        registry.fill(HIST("hPtCandRec"), candidate.pt());
        registry.fill(HIST("hPtProng0Rec"), candidate.ptProng0());
        registry.fill(HIST("hPtProng1Rec"), candidate.ptProng1());
        registry.fill(HIST("hEtaRec"), candidate.eta());
        registry.fill(HIST("hPhiRec"), candidate.phi());
        registry.fill(HIST("hYRec"), HfHelper::yD0(candidate));
        registry.fill(HIST("hSelectionStatusRec"), candidate.isSelD0bar() + (candidate.isSelD0() * 2));
      }
      // fill invariant mass plots from D0/D0bar signal and background candidates
      if (candidate.isSelD0() >= selectionFlagD0) {                                                  // only reco as D0
        if (candidate.flagMcMatchRec() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) { // also matched as D0
          registry.fill(HIST("hMassD0RecSig"), invMassD0, candidate.pt(), efficiencyWeight);
          if (isD0Prompt) {
            registry.fill(HIST("hPtCandRecSigPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityRecPrompt"), candidate.pt(), collision.multFT0M());
            registry.fill(HIST("hPtVsMLScoresVsEtaRecSigPrompt"), outputMlD0[0], outputMlD0[1], outputMlD0[2], candidate.pt(), candidate.eta());
          } else if (isD0NonPrompt) {
            registry.fill(HIST("hPtCandRecSigNonPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityRecNonPrompt"), candidate.pt(), collision.multFT0M());
            registry.fill(HIST("hPtVsMLScoresVsEtaRecSigNonPrompt"), outputMlD0[0], outputMlD0[1], outputMlD0[2], candidate.pt(), candidate.eta());
          }
        } else if (candidate.flagMcMatchRec() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          registry.fill(HIST("hMassD0RecRef"), invMassD0, candidate.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0RecBg"), invMassD0, candidate.pt(), efficiencyWeight);
        }
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMlD0[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
        }
        registry.fill(HIST("hMLScoresVsMassVsPtVsEtaVsOrigin"), outputMlD0[0], outputMlD0[1], outputMlD0[2], invMassD0, candidate.pt(), candidate.eta(), isD0Prompt);
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {                                             // only reco as D0bar
        if (candidate.flagMcMatchRec() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) { // also matched as D0bar
          registry.fill(HIST("hMassD0barRecSig"), invMassD0bar, candidate.pt(), efficiencyWeight);
          if (isD0Prompt) {
            registry.fill(HIST("hPtCandRecSigPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityRecPrompt"), candidate.pt(), collision.multFT0M());
            registry.fill(HIST("hPtVsMLScoresVsEtaRecSigPrompt"), outputMlD0bar[0], outputMlD0bar[1], outputMlD0bar[2], candidate.pt(), candidate.eta());
          } else if (isD0NonPrompt) {
            registry.fill(HIST("hPtCandRecSigNonPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityRecNonPrompt"), candidate.pt(), collision.multFT0M());
            registry.fill(HIST("hPtVsMLScoresVsEtaRecSigNonPrompt"), outputMlD0bar[0], outputMlD0bar[1], outputMlD0bar[2], candidate.pt(), candidate.eta());
          }
        } else if (candidate.flagMcMatchRec() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          registry.fill(HIST("hMassD0barRecRef"), invMassD0bar, candidate.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0barRecBg"), invMassD0bar, candidate.pt(), efficiencyWeight);
        }
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMlD0bar[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
        }
        registry.fill(HIST("hMLScoresVsMassVsPtVsEtaVsOrigin"), outputMlD0bar[0], outputMlD0bar[1], outputMlD0bar[2], invMassD0bar, candidate.pt(), candidate.eta(), isD0Prompt);
      }
      entryD0CandRecoInfo(invMassD0, invMassD0bar, candidate.pt(), outputMlD0[0], outputMlD0[2], outputMlD0bar[0], outputMlD0bar[2]);
      entryD0CandGenInfo(isD0Prompt);

      // ===================== Define parameters for soft pion removal ========================
      auto ePiK = RecoDecay::e(candidate.pVectorProng0(), massPi) + RecoDecay::e(candidate.pVectorProng1(), massK);
      auto eKPi = RecoDecay::e(candidate.pVectorProng0(), massK) + RecoDecay::e(candidate.pVectorProng1(), massPi);

      // ============== D-h correlation dedicated section ====================================

      flagD0 = candidate.flagMcMatchRec() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;     // flagD0Signal 'true' if candidate matched to D0 (particle)
      flagD0bar = candidate.flagMcMatchRec() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK; // flagD0Reflection 'true' if candidate, selected as D0 (particle), is matched to D0bar (antiparticle)
      float cent = 100.0;                                                                                 // Centrality Placeholder: will be updated later
      // ========== track loop starts here ========================

      for (const auto& track : tracks) {
        registry.fill(HIST("hTrackCounter"), 0); // fill total no. of tracks
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }
        // Removing D0 daughters by checking track indices
        bool correlationStatus = false;
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex())) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }
        registry.fill(HIST("hTrackCounter"), 1); // fill no. of tracks before soft pion removal

        bool isPhysicalPrimary = false;
        // ===== soft pion removal ===================================================
        double invMassDstar1 = 0, invMassDstar2 = 0;
        bool isSoftPiD0 = false, isSoftPiD0bar = false;
        auto pSum2 = RecoDecay::p2(candidate.pVector(), track.pVector());
        auto ePion = track.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (candidate.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - invMassD0) - softPiMass) < ptSoftPionMax) {
            continue;
          }
        }

        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - invMassD0bar) - softPiMass) < ptSoftPionMax) {
            continue;
          }
        }

        registry.fill(HIST("hTrackCounter"), 2); // fill no. of tracks after soft pion removal

        if (correlateD0WithLeadingParticle) {
          if (track.globalIndex() != leadingIndex) {
            continue;
          }
          registry.fill(HIST("hTrackCounter"), 3); // fill no. of tracks  have leading particle
        }

        int signalStatus = 0;
        if (flagD0 && (candidate.isSelD0() >= selectionFlagD0) && !isSoftPiD0) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Sig);
        } // signal case D0
        if (flagD0bar && (candidate.isSelD0() >= selectionFlagD0) && !isSoftPiD0) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Ref);
        } // reflection case D0
        if (!flagD0 && !flagD0bar && (candidate.isSelD0() >= selectionFlagD0) && !isSoftPiD0) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Bg);
        } // background case D0

        if (flagD0bar && (candidate.isSelD0bar() >= selectionFlagD0bar) && !isSoftPiD0bar) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barSig);
        } // signal case D0bar
        if (flagD0 && (candidate.isSelD0bar() >= selectionFlagD0bar) && !isSoftPiD0bar) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barRef);
        } // reflection case D0bar
        if (!flagD0 && !flagD0bar && (candidate.isSelD0bar() >= selectionFlagD0bar) && !isSoftPiD0bar) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barBg);
        } // background case D0bar

        entryD0HadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                          track.eta() - candidate.eta(),
                          candidate.pt(),
                          track.pt(),
                          poolBin,
                          correlationStatus,
                          cent);
        entryD0HadronRecoInfo(invMassD0, invMassD0bar, signalStatus);
        entryD0HadronMlInfo(outputMlD0[0], outputMlD0[1], outputMlD0[2], outputMlD0bar[0], outputMlD0bar[1], outputMlD0bar[2]);
        if (track.has_mcParticle()) {
          auto mcParticle = track.template mcParticle_as<aod::McParticles>();
          isPhysicalPrimary = mcParticle.isPhysicalPrimary();
          auto trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
          entryD0HadronGenInfo(isD0Prompt, isPhysicalPrimary, trackOrigin);
        } else {
          entryD0HadronGenInfo(isD0Prompt, isPhysicalPrimary, 0);
          registry.fill(HIST("hTrackCounter"), 4); // fill no. of fake tracks
        }
        // for secondary particle fraction estimation
        registry.fill(HIST("hPtParticleAssocVsCandRec"), track.pt(), candidate.pt());
        if (isPhysicalPrimary) {
          registry.fill(HIST("hPtPrimaryParticleAssocVsCandRec"), track.pt(), candidate.pt());
        }
        entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
      } // end inner loop (Tracks)
    } // end of outer loop (D0)
  }

  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcRec, "Process MC Reco mode", true);

  // =================  Process starts for MCGen, same event ===================

  void processMcGen(SelectedCollisionsMcGen::iterator const& mcCollision,
                    SelectedParticlesMcGen const& mcParticles)
  {
    BinningTypeMcGen const corrBinningMcGen{{zPoolBins, multPoolBinsMcGen}, true};
    int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), mcCollision.multMCFT0A()));
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hEvtCountGen"), 0);
    // MC gen level
    // find leading particle
    if (correlateD0WithLeadingParticle) {
      leadingIndex = findLeadingParticleMcGen(mcParticles, etaTrackMax.value, ptTrackMin.value);
    }
    bool isD0Prompt = false;
    bool isD0NonPrompt = false;
    int trackOrigin = -1;
    float cent = 100.; // Centrality Placeholder: will be updated later

    for (const auto& particleTrigg : mcParticles) {
      if (std::abs(particleTrigg.pdgCode()) != Pdg::kD0) {
        continue;
      }
      if (std::abs(particleTrigg.flagMcMatchGen()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
        double const yD = RecoDecay::y(particleTrigg.pVector(), MassD0);
        if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particleTrigg.pt() < ptCandMin) {
          continue;
        }

        registry.fill(HIST("hD0PoolBin"), poolBin);
        registry.fill(HIST("hPtCandGen"), particleTrigg.pt());
        registry.fill(HIST("hEtaGen"), particleTrigg.eta());
        registry.fill(HIST("hPhiGen"), particleTrigg.phi());
        registry.fill(HIST("hYGen"), yD);
        registry.fill(HIST("hCountD0TriggersGen"), 0, particleTrigg.pt()); // to count trigger D0 (for normalisation)

        isD0Prompt = particleTrigg.originMcGen() == RecoDecay::OriginType::Prompt;
        isD0NonPrompt = particleTrigg.originMcGen() == RecoDecay::OriginType::NonPrompt;

        // prompt and non-prompt division
        if (isD0Prompt) {
          registry.fill(HIST("hPtCandGenPrompt"), particleTrigg.pt());
          registry.fill(HIST("hPtVsEtaCandGenSigPrompt"), particleTrigg.pt(), particleTrigg.eta());
        } else if (isD0NonPrompt) {
          registry.fill(HIST("hPtCandGenNonPrompt"), particleTrigg.pt());
          registry.fill(HIST("hPtVsEtaCandGenSigNonPrompt"), particleTrigg.pt(), particleTrigg.eta());
        }

        // =============== D-h correlation dedicated section =====================
        for (const auto& particleAssoc : mcParticles) {
          registry.fill(HIST("hTrackCounter"), 0); // total no. of tracks
          if (std::abs(particleAssoc.eta()) > etaTrackMax) {
            continue;
          }
          if (particleAssoc.pt() < ptTrackMin) {
            continue;
          }
          if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
            continue;
          }
          if (!particleAssoc.isPhysicalPrimary()) {
            continue;
          }
          // ==============================soft pion removal================================
          registry.fill(HIST("hTrackCounter"), 1); // fill before soft pi removal
          // method used: indexMother = -1 by default if the mother doesn't match with given PID of the mother. We find mother of pion if it is D* and mother of D0 if it is D*. If they are both positive and they both match each other, then it is detected as a soft pion

          auto indexMotherPi = RecoDecay::getMother(mcParticles, particleAssoc, Pdg::kDStar, true, nullptr, 1); // last arguement 1 is written to consider immediate decay mother only
          auto indexMotherD0 = RecoDecay::getMother(mcParticles, particleTrigg, Pdg::kDStar, true, nullptr, 1);
          bool correlationStatus = false;
          if (std::abs(particleAssoc.pdgCode()) == kPiPlus && indexMotherPi >= 0 && indexMotherD0 >= 0 && indexMotherPi == indexMotherD0) {
            if (!storeAutoCorrelationFlag) {
              continue;
            }
            correlationStatus = true;
          }

          registry.fill(HIST("hTrackCounter"), 2); // fill after soft pion removal

          if (correlateD0WithLeadingParticle) {
            if (particleAssoc.globalIndex() != leadingIndex) {
              continue;
            }
            registry.fill(HIST("hTrackCounter"), 3); // fill no. of tracks  have leading particle
          }
          trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
          entryD0HadronPair(getDeltaPhi(particleAssoc.phi(), particleTrigg.phi()),
                            particleAssoc.eta() - particleTrigg.eta(),
                            particleTrigg.pt(),
                            particleAssoc.pt(),
                            poolBin,
                            correlationStatus,
                            cent);
          entryD0HadronRecoInfo(massD0, massD0, 0); // dummy info
          entryD0HadronGenInfo(isD0Prompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
        } // end inner loop (Tracks)
      }
    } // end outer loop (D0)
  }

  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcGen, "Process MC Gen mode", false);

  // ====================== Implement Event mixing on Data ===================================

  void processDataMixedEvent(SelectedCollisions const& collisions,
                             SelectedCandidatesDataMl const& candidates,
                             SelectedTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      registry.fill(HIST("hMultFT0M"), collision.multFT0M());
      registry.fill(HIST("hZvtx"), collision.posZ());
    }

    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelectedCollisions, SelectedCandidatesDataMl, SelectedTracks, BinningType> const pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      if (tracks1.size() == 0) {
        continue;
      }
      // LOGF(info, "Mixed event collisions: Index = (%d, %d), tracks Size: (%d, %d), Z Vertex: (%f, %f), Pool Bin: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), c1.posZ(), c2.posZ(), corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M())),corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()))); // For debug
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      int const poolBinD0 = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hTracksPoolBin"), poolBin);
      registry.fill(HIST("hD0PoolBin"), poolBinD0);
      for (const auto& [candidate, particleAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (std::abs(HfHelper::yD0(candidate)) >= yCandMax || candidate.pt() < ptCandMin) {
          continue;
        }

        // soft pion removal, signal status 1,3 for D0 and 2,3 for D0bar (SoftPi removed), signal status 11,13 for D0  and 12,13 for D0bar (only SoftPi)
        const auto invMassD0 = HfHelper::invMassD0ToPiK(candidate);
        const auto invMassD0bar = HfHelper::invMassD0barToKPi(candidate);
        auto ePiK = RecoDecay::e(candidate.pVectorProng0(), massPi) + RecoDecay::e(candidate.pVectorProng1(), massK);
        auto eKPi = RecoDecay::e(candidate.pVectorProng0(), massK) + RecoDecay::e(candidate.pVectorProng1(), massPi);
        double invMassDstar1 = 0., invMassDstar2 = 0.;
        bool isSoftPiD0 = false, isSoftPiD0bar = false;
        auto pSum2 = RecoDecay::p2(candidate.pVector(), particleAssoc.pVector());
        auto ePion = particleAssoc.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);
        std::vector<float> outputMlD0 = {-1., -1., -1.};
        std::vector<float> outputMlD0bar = {-1., -1., -1.};
        float cent = 100.; // Centrality Placeholder: will be updated later

        if (candidate.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - invMassD0) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0 = true;
          }
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMlD0[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
          }
        }
        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - invMassD0bar) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0bar = true;
          }
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMlD0bar[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
          }
        }

        int signalStatus = 0;
        if (candidate.isSelD0() >= selectionFlagD0) {
          if (!isSoftPiD0) {
            signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0Only;
          } else {
            signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0OnlySoftPi;
          }
        }
        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          if (!isSoftPiD0bar) {
            signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0barOnly;
          } else {
            signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0barOnlySoftPi;
          }
        }
        bool correlationStatus = false;
        entryD0HadronPair(getDeltaPhi(candidate.phi(), particleAssoc.phi()), candidate.eta() - particleAssoc.eta(), candidate.pt(), particleAssoc.pt(), poolBin, correlationStatus, cent);
        entryD0HadronRecoInfo(invMassD0, invMassD0bar, signalStatus);
        entryD0HadronGenInfo(false, false, 0);
        entryD0HadronMlInfo(outputMlD0[0], outputMlD0[1], outputMlD0[2], outputMlD0bar[0], outputMlD0bar[1], outputMlD0bar[2]);
        entryTrackRecoInfo(particleAssoc.dcaXY(), particleAssoc.dcaZ(), particleAssoc.tpcNClsCrossedRows());
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processDataMixedEvent, "Process data mixed event", false);

  // ====================== Implement Event mixing on McRec ===================================

  void processMcRecMixedEvent(SelectedCollisions const& collisions,
                              SelectedCandidatesMcRecMl const& candidates,
                              SelectedTracksMcRec const& tracks,
                              aod::McParticles const& mcParticles)
  {
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelectedCollisions, SelectedCandidatesMcRecMl, SelectedTracksMcRec, BinningType> const pairMcRec{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    bool isD0Prompt = false;
    bool flagD0 = false;
    bool flagD0bar = false;
    bool isPhysicalPrimary = false;
    int trackOrigin = 0;
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      int const poolBinD0 = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hTracksPoolBin"), poolBin);
      registry.fill(HIST("hD0PoolBin"), poolBinD0);
      registry.fill(HIST("hMultFT0M"), c1.multFT0M());
      registry.fill(HIST("hZvtx"), c1.posZ());
      float cent = 100.; // Centrality Placeholder: will be updated later

      for (const auto& [candidate, particleAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (std::abs(HfHelper::yD0(candidate)) >= yCandMax || candidate.pt() < ptCandMin) {
          continue;
        }
        if (!particleAssoc.isGlobalTrackWoDCA()) {
          continue;
        }
        std::vector<float> outputMlD0 = {-1., -1., -1.};
        std::vector<float> outputMlD0bar = {-1., -1., -1.};
        isD0Prompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        if (particleAssoc.has_mcParticle()) {
          auto mcParticle = particleAssoc.template mcParticle_as<aod::McParticles>();
          isPhysicalPrimary = mcParticle.isPhysicalPrimary();
          trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
        }
        if (candidate.isSelD0() >= selectionFlagD0) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMlD0[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
          }
        } else if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMlD0bar[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
          }
        }

        // soft pion removal
        const auto invMassD0 = HfHelper::invMassD0ToPiK(candidate);
        const auto invMassD0bar = HfHelper::invMassD0barToKPi(candidate);
        auto ePiK = RecoDecay::e(candidate.pVectorProng0(), massPi) + RecoDecay::e(candidate.pVectorProng1(), massK);
        auto eKPi = RecoDecay::e(candidate.pVectorProng0(), massK) + RecoDecay::e(candidate.pVectorProng1(), massPi);
        double invMassDstar1 = 0., invMassDstar2 = 0.;
        bool isSoftPiD0 = false, isSoftPiD0bar = false;
        auto pSum2 = RecoDecay::p2(candidate.pVector(), particleAssoc.pVector());
        auto ePion = particleAssoc.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (candidate.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - invMassD0) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0 = true;
          }
        }

        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - invMassD0bar) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0bar = true;
          }
        }

        flagD0 = candidate.flagMcMatchRec() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;     // flagD0Signal 'true' if candidate matched to D0 (particle)
        flagD0bar = candidate.flagMcMatchRec() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK; // flagD0Reflection 'true' if candidate, selected as D0 (particle), is matched to D0bar (antiparticle)
        int signalStatus = 0;

        if (flagD0 && (candidate.isSelD0() >= selectionFlagD0)) {
          if (!isSoftPiD0) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Sig); //  signalStatus += 1;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi); // signalStatus += 64;
          }
        } // signal case D0

        if (flagD0bar && (candidate.isSelD0() >= selectionFlagD0)) {
          if (!isSoftPiD0) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Ref); //   signalStatus += 2;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi); // signalStatus += 64;
          }
        } // reflection case D0

        if (!flagD0 && !flagD0bar && (candidate.isSelD0() >= selectionFlagD0)) {
          if (!isSoftPiD0) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Bg); //  signalStatus += 4;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi);
          }
        } // background case D0

        if (flagD0bar && (candidate.isSelD0bar() >= selectionFlagD0bar)) {
          if (!isSoftPiD0bar) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barSig); //  signalStatus += 8;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi);
          }
        } // signal case D0bar

        if (flagD0 && (candidate.isSelD0bar() >= selectionFlagD0bar)) {
          if (!isSoftPiD0bar) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barRef); // signalStatus += 16;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi);
          }
        } // reflection case D0bar

        if (!flagD0 && !flagD0bar && (candidate.isSelD0bar() >= selectionFlagD0bar)) {
          if (!isSoftPiD0bar) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barBg); //   signalStatus += 32;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi);
          }
        } // background case D0bar
        registry.fill(HIST("hSignalStatusMERec"), signalStatus);
        bool correlationStatus = false;
        entryD0HadronPair(getDeltaPhi(candidate.phi(), particleAssoc.phi()), candidate.eta() - particleAssoc.eta(), candidate.pt(), particleAssoc.pt(), poolBin, correlationStatus, cent);
        entryD0HadronRecoInfo(invMassD0, invMassD0bar, signalStatus);
        entryD0HadronGenInfo(isD0Prompt, isPhysicalPrimary, trackOrigin);
        entryD0HadronMlInfo(outputMlD0[0], outputMlD0[1], outputMlD0[2], outputMlD0bar[0], outputMlD0bar[1], outputMlD0bar[2]);
        entryTrackRecoInfo(particleAssoc.dcaXY(), particleAssoc.dcaZ(), particleAssoc.tpcNClsCrossedRows());
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcRecMixedEvent, "Process Mixed Event MCRec", false);

  // ====================== Implement Event mixing on McGen ===================================

  void processMcGenMixedEvent(SelectedCollisionsMcGen const& collisions,
                              SelectedParticlesMcGen const& mcParticles)
  {
    BinningTypeMcGen const corrBinningMcGen{{zPoolBins, multPoolBinsMcGen}, true};
    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<SelectedCollisionsMcGen, SelectedParticlesMcGen, SelectedParticlesMcGen, BinningTypeMcGen> const pairMcGen{corrBinningMcGen, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      int poolBin = corrBinningMcGen.getBin(std::make_tuple(c1.posZ(), c1.multMCFT0A()));
      for (const auto& [particleTrigg, particleAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (std::abs(particleTrigg.pdgCode()) != Pdg::kD0) {
          continue;
        }
        if (std::abs(particleTrigg.flagMcMatchGen()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          double const yD = RecoDecay::y(particleTrigg.pVector(), MassD0);
          if (std::abs(yD) >= yCandMax || particleTrigg.pt() <= ptCandMin || std::abs(particleAssoc.eta()) >= etaTrackMax || particleAssoc.pt() <= ptTrackMin) {
            continue;
          }
          if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
            continue;
          }
          if (!particleAssoc.isPhysicalPrimary()) {
            continue;
          }

          // ==============================soft pion removal================================
          // method used: indexMother = -1 by default if the mother doesn't match with given PID of the mother. We find mother of pion if it is D* and mother of D0 if it is D*. If they are both positive and they both match each other, then it is detected as a soft pion

          auto indexMotherPi = RecoDecay::getMother(mcParticles, particleAssoc, Pdg::kDStar, true, nullptr, 1); // last arguement 1 is written to consider immediate decay mother only
          auto indexMotherD0 = RecoDecay::getMother(mcParticles, particleTrigg, Pdg::kDStar, true, nullptr, 1);
          if (std::abs(particleAssoc.pdgCode()) == kPiPlus && indexMotherPi >= 0 && indexMotherD0 >= 0 && indexMotherPi == indexMotherD0) {
            continue;
          }
          float cent = 100.; // Centrality Placeholder: will be updated later
          bool correlationStatus = false;
          int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
          bool isD0Prompt = particleTrigg.originMcGen() == RecoDecay::OriginType::Prompt;
          entryD0HadronPair(getDeltaPhi(particleAssoc.phi(), particleTrigg.phi()), particleAssoc.eta() - particleTrigg.eta(), particleTrigg.pt(), particleAssoc.pt(), poolBin, correlationStatus, cent);
          entryD0HadronRecoInfo(massD0, massD0, 0); // dummy info
          entryD0HadronGenInfo(isD0Prompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcGenMixedEvent, "Process Mixed Event MCGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCorrelatorD0HadronsSelection>(cfgc),
    adaptAnalysisTask<HfCorrelatorD0Hadrons>(cfgc)};
}
