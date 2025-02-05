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

/// \file correlatorDsHadrons.cxx
/// \brief Ds-Hadrons correlator task - data-like, MC-reco and MC-Gen analyses
/// \author Grazia Luparello <grazia.luparello@cern.ch>
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
double getDeltaPhi(double phiHadron, double phiD)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -PIHalf);
}

// binning type
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
using BinningTypeMcGen = ColumnBinningPolicy<aod::mccollision::PosZ, o2::aod::mult::MultMCFT0A>;

/// Code to select collisions with at least one Ds meson
struct HfCorrelatorDsHadronsSelCollision {
  Produces<aod::DmesonSelection> collisionsWithSelDs;

  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<bool> doSelDsCollision{"doSelDsCollision", true, "Select collisions with at least one Ds"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  HfHelper hfHelper;
  SliceCache cache;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;

  Filter dsFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);

  /// Code to select collisions with at least one Ds meson - for real data and data-like analysis
  void processDsSelCollisionsData(SelCollisions::iterator const& collision,
                                  CandDsData const& candidates)
  {
    bool isSelColl = true;
    bool isDsFound = true;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    if (doSelDsCollision) {
      isDsFound = false; // if candidate table is empty for-loop is not performed
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          isDsFound = false;
          continue;
        }
        isDsFound = true;
        break;
      }
    }
    if (useSel8) {
      isSel8 = false;
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = false;
      isNosameBunchPileUp = static_cast<bool>(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup));
    }
    isSelColl = isDsFound && isSel8 && isNosameBunchPileUp;
    collisionsWithSelDs(isSelColl);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsSelCollision, processDsSelCollisionsData, "Process Ds Collision Selection Data", true);

  /// Code to select collisions with at least one Ds meson - for MC-level analysis
  void processDsSelCollisionsMc(SelCollisions::iterator const& collision,
                                CandDsMcReco const& candidates)
  {
    bool isSelColl = true;
    bool isDsFound = true;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    if (doSelDsCollision) { // to enable only for the MC reco part
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          isDsFound = false;
          continue;
        }
        isDsFound = true;
        break;
      }
    }
    if (useSel8) {
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = static_cast<bool>(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup));
    }
    isSelColl = isDsFound && isSel8 && isNosameBunchPileUp;
    collisionsWithSelDs(isSelColl);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsSelCollision, processDsSelCollisionsMc, "Process Ds Collision Selection MCRec", false);
};

/// Ds-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDsHadrons {
  Produces<aod::DsHadronPair> entryDsHadronPair;
  Produces<aod::DsHadronRecoInfo> entryDsHadronRecoInfo;
  Produces<aod::DsHadronGenInfo> entryDsHadronGenInfo;
  Produces<aod::DsHadronMlInfo> entryDsHadronMlInfo;
  Produces<aod::DsCandRecoInfo> entryDsCandRecoInfo;
  Produces<aod::DsCandGenInfo> entryDsCandGenInfo;
  Produces<aod::TrackRecoInfo> entryTrackRecoInfo;
  Produces<aod::HfcRedCollisions> collReduced;
  Produces<aod::DsCandReduceds> candReduced;
  Produces<aod::AssocTrackReds> assocTrackReduced;

  Configurable<bool> fillHistoData{"fillHistoData", true, "Flag for filling histograms in data processes"};
  Configurable<bool> fillHistoMcRec{"fillHistoMcRec", true, "Flag for filling histograms in MC Rec processes"};
  Configurable<bool> fillHistoMcGen{"fillHistoMcGen", true, "Flag for filling histograms in MC Gen processes"};
  Configurable<bool> removeCollWSplitVtx{"removeCollWSplitVtx", false, "Flag for rejecting the splitted collisions"};
  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection (used only in MC processes)"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing (used only in MC processes)"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds (avoid the case of flag = 0, no outputMlScore)"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<int> decayChannel{"decayChannel", 1, "Decay channels: 1 for Ds->PhiPi->KKpi, 2 for Ds->K0*K->KKPi"};
  Configurable<bool> applyEfficiency{"applyEfficiency", true, "Flag for applying D-meson efficiency weights"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 2., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 2., "max. DCA_z of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", {1., 1., 1., 1., 1., 1.}, "efficiency values for Ds meson"};

  int hfcReducedCollisionIndex = 0;

  HfHelper hfHelper;
  SliceCache cache;

  using SelCollisionsWithDs = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::DmesonSelection>>; // collisionFilter applied
  // using SelCollisionsWithDsWithMc = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::DmesonSelection, aod::McCollisionLabels>>; // collisionFilter applied
  using SelCollisionsMc = soa::Join<aod::McCollisions, aod::MultsExtraMC>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;                           // flagDsFilter applied
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi, aod::HfCand3ProngMcRec>>; // flagDsFilter applied
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;                                                         // flagDsFilter applied
  using MyTracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;                           // trackFilter applied
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels>>;   // trackFilter applied

  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter flagDsFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  Preslice<CandDsData> candsDsPerCollision = aod::hf_cand::collisionId;
  Preslice<MyTracksData> trackIndicesPerCollision = aod::track::collisionId;
  Preslice<CandDsMcGen> perCollisionCandMc = o2::aod::mcparticle::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels, aod::McCollisionLabels>> collPerCollMc = o2::aod::mccollisionlabel::mcCollisionId;

  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 900., 1800., 6000.}, "event multiplicity pools (FT0M)"};
  ConfigurableAxis binsMassD{"binsMassD", {200, 1.7, 2.25}, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "#it{#eta}"};

  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {200, 0., 800.}, "Multiplicity"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 6000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMassD = {binsMassD, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
    AxisSpec axisEta = {binsEta, "#it{#eta}"};
    AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec axisMultiplicity = {binsMultiplicity, "Multiplicity"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec axisStatus = {15, 0.5, 15.5, "Selection status"};

    // Histograms for data analysis
    if (fillHistoData) {
      registry.add("hPtCand", "Ds candidates pt", {HistType::kTH1F, {axisPtD}});
      registry.add("hSelectionStatusDsToKKPi", "Ds candidates selection", {HistType::kTH1F, {axisStatus}});
      registry.add("hSelectionStatusDsToPiKK", "Ds candidates selection", {HistType::kTH1F, {axisStatus}});
      registry.add("hCountSelectionStatusDsToKKPiAndToPiKK", "Ds candidates selection", {HistType::kTH1F, {{1, -0.5, 0.5}}});
      registry.add("hEta", "Ds candidates eta", {HistType::kTH1F, {axisEta}});
      registry.add("hEtaVsPtCand", "Ds candidates etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtD}}});
      registry.add("hEtaVsPtPartAssoc", "Particles associated etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtHadron}}});
      registry.add("hPhi", "Ds candidates phi", {HistType::kTH1F, {axisPhi}});
      registry.add("hPhiVsPtCand", "Ds candidates phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtD}}});
      registry.add("hPhiVsPtPartAssoc", "Particles associated phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtHadron}}});
      registry.add("hMultiplicity", "Multiplicity", {HistType::kTH1F, {axisMultiplicity}});
      registry.add("hMultFT0M", "Multiplicity FT0M", {HistType::kTH1F, {axisMultFT0M}});
      registry.add("hZVtx", "z vertex", {HistType::kTH1F, {axisPosZ}});
      registry.add("hMassDsVsPt", "Ds candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
      registry.add("hMassDsData", "Ds candidates mass", {HistType::kTH1F, {axisMassD}});
      registry.add("hCollisionPoolBin", "Ds candidates collision pool bin", {HistType::kTH1F, {axisPoolBin}});
      registry.add("hDsPoolBin", "Ds candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
      registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
    }
    // Histograms for MC Reco analysis
    if (fillHistoMcRec) {
      registry.add("hPtCandMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcRecSigPrompt", "Ds,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcRecSigNonPrompt", "Ds,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}});
      registry.add("hSelectionStatusDsToKKPiMcRec", "Ds candidates selection", {HistType::kTH1F, {axisStatus}});
      registry.add("hSelectionStatusDsToPiKKMcRec", "Ds candidates selection", {HistType::kTH1F, {axisStatus}});
      registry.add("hCountSelectionStatusDsToKKPiAndToPiKKMcRec", "Ds candidates selection", {HistType::kTH1F, {{1, -0.5, 0.5}}});
      registry.add("hPtParticleAssocMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtD}}});
      registry.add("hPtPrimaryParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtD}}});
      registry.add("hEtaMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
      registry.add("hPhiMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
      registry.add("hEtaMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
      registry.add("hPhiMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
      registry.add("hPtVsMultiplicityMcRecPrompt", "Multiplicity FT0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultFT0M}}});
      registry.add("hPtVsMultiplicityMcRecNonPrompt", "Multiplicity FT0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultFT0M}}});
      registry.add("hMassDsMcRec", "Ds candidates", {HistType::kTH1F, {axisMassD}});
      registry.add("hMassDsVsPtMcRec", "Ds signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
      registry.add("hMassDsMcRecSig", "Ds signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
      registry.add("hMassDsMcRecBkg", "Ds background candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
      registry.add("hFakeTracksMcRec", "Fake tracks - MC Rec", {HistType::kTH1F, {axisPtHadron}});
    }
    // Histograms for MC Gen analysis
    if (fillHistoMcGen) {
      registry.add("hPtCandMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcGenPrompt", "Ds,Hadron particles - MC Gen Prompt", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcGenNonPrompt", "Ds,Hadron particles - MC Gen Non Prompt", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hEtaMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisEta}});
      registry.add("hPhiMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisPhi}});
      registry.add("hMultFT0AMcGen", "Ds,Hadron multiplicity FT0A - MC Gen", {HistType::kTH1F, {axisMultiplicity}});
    }
  }

  /// Fill histograms of quantities independent from the daugther-mass hypothesis for data
  /// \param candidate is candidate
  template <typename T1>
  void fillHisto(const T1& candidate)
  {
    registry.fill(HIST("hPtCand"), candidate.pt());
    registry.fill(HIST("hEta"), candidate.eta());
    registry.fill(HIST("hEtaVsPtCand"), candidate.eta(), candidate.pt());
    registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
    registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf), candidate.pt());
  }

  /// Fill histograms of quantities for the KKPi daugther-mass hypothesis for data
  /// \param candidate is candidate
  /// \param efficiencyWeight is the efficiency correction
  template <typename T1>
  void fillHistoKKPi(const T1& candidate, double efficiencyWeight)
  {
    registry.fill(HIST("hMassDsVsPt"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
    registry.fill(HIST("hMassDsData"), hfHelper.invMassDsToKKPi(candidate), efficiencyWeight);
    registry.fill(HIST("hSelectionStatusDsToKKPi"), candidate.isSelDsToKKPi());
  }

  /// Fill histograms of quantities for the PiKK daugther-mass hypothesis for data
  /// \param candidate is candidate
  /// \param efficiencyWeight is the efficiency correction
  template <typename T1>
  void fillHistoPiKK(const T1& candidate, double efficiencyWeight)
  {
    registry.fill(HIST("hMassDsVsPt"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
    registry.fill(HIST("hMassDsData"), hfHelper.invMassDsToPiKK(candidate), efficiencyWeight);
    registry.fill(HIST("hSelectionStatusDsToPiKK"), candidate.isSelDsToPiKK());
  }

  /// Fill histograms of quantities for the Ds signal for MC reco-level
  /// \param candidate is candidate
  /// \param multiplicityFT0M is the multiplicity
  template <typename T1>
  void fillHistoMcRecSig(const T1& candidate, float multiplicityFT0M)
  {
    registry.fill(HIST("hPtCandMcRecSig"), candidate.pt());
    registry.fill(HIST("hEtaMcRecSig"), candidate.eta());
    registry.fill(HIST("hPhiMcRecSig"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));

    // prompt and non-prompt division
    if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
      registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), multiplicityFT0M);
    } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
      registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), multiplicityFT0M);
    }
  }

  /// Fill histograms of quantities for the Ds backgroung for MC reco-level
  /// \param candidate is candidate
  template <typename T1>
  void fillHistoMcRecBkg(const T1& candidate)
  {
    registry.fill(HIST("hPtCandMcRecBkg"), candidate.pt());
    registry.fill(HIST("hEtaMcRecBkg"), candidate.eta());
    registry.fill(HIST("hPhiMcRecBkg"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
  }

  /// Fill histograms of quantities for the Ds signal for MC reco-level
  /// \param particle is particle, Ds
  template <typename T1>
  void fillMcGenHisto(const T1& particle)
  {
    registry.fill(HIST("hPtCandMcGen"), particle.pt());
    registry.fill(HIST("hEtaMcGen"), particle.eta());
    registry.fill(HIST("hPhiMcGen"), RecoDecay::constrainAngle(particle.phi(), -PIHalf));

    // prompt and non-prompt division
    if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtCandMcGenPrompt"), particle.pt());
    } else if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtCandMcGenNonPrompt"), particle.pt());
    }
  }

  /// Ds-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(SelCollisionsWithDs::iterator const& collision,
                   CandDsData const& candidates,
                   MyTracksData const& tracks)
  {
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    registry.fill(HIST("hZVtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    int nTracks = tracks.size();
    registry.fill(HIST("hMultiplicity"), nTracks);

    // Ds fill histograms and Ds-Hadron correlation for DsToKKPi
    for (const auto& candidate : candidates) {
      if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      double efficiencyWeightD = 1.;
      if (applyEfficiency) {
        efficiencyWeightD = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()));
      }
      std::vector<float> outputMl = {-1., -1., -1.};
      fillHisto(candidate);
      if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
        fillHistoKKPi(candidate, efficiencyWeightD);
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
        }
        entryDsCandRecoInfo(hfHelper.invMassDsToKKPi(candidate), candidate.pt(), outputMl[0], outputMl[2]);
      } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
        fillHistoPiKK(candidate, efficiencyWeightD);
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
        }
        entryDsCandRecoInfo(hfHelper.invMassDsToPiKK(candidate), candidate.pt(), outputMl[0], outputMl[2]);
      }
      if (candidate.isSelDsToKKPi() >= selectionFlagDs && candidate.isSelDsToPiKK() >= selectionFlagDs) {
        registry.fill(HIST("hCountSelectionStatusDsToKKPiAndToPiKK"), 0.);
      }

      // Ds-Hadron correlation dedicated section
      for (const auto& track : tracks) {
        // Removing Ds daughters by checking track indices
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
          continue;
        }
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }

        registry.fill(HIST("hEtaVsPtPartAssoc"), track.eta(), candidate.pt());
        registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(track.phi(), -PIHalf), candidate.pt());
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(candidate), false, false);
          // entryDsHadronGenInfo(false, false, 0);
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), false, false);
          // entryDsHadronGenInfo(false, false, 0);
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
        }
      } // end track loop
    }   // end candidate loop
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processData, "Process data", true);

  /// Ds-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(SelCollisionsWithDs::iterator const& collision,
                    CandDsMcReco const& candidates,
                    TracksWithMc const& tracks,
                    aod::McParticles const& mcParticles)
  {
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    registry.fill(HIST("hZVtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    registry.fill(HIST("hCollisionPoolBin"), poolBin);

    // MC reco level
    bool isDsPrompt = false;
    bool isDsSignal = false;
    bool isCorrectInvMassHypo = false;
    bool isDecayChan = false;
    bool isAlreadyFilledEvent = false;
    float multiplicityFT0M = collision.multFT0M();
    for (const auto& candidate : candidates) {
      // prompt and non-prompt division
      isDsPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
      // Ds Signal
      isDsSignal = std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi;
      isDecayChan = candidate.flagMcDecayChanRec() == decayChannel;

      if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }

      auto prong0McPart = candidate.template prong0_as<TracksWithMc>().template mcParticle_as<aod::McParticles>();
      isCorrectInvMassHypo = ((std::abs(prong0McPart.pdgCode()) == kKPlus) && (candidate.isSelDsToKKPi() >= selectionFlagDs)) || ((std::abs(prong0McPart.pdgCode()) == kPiPlus) && (candidate.isSelDsToPiKK() >= selectionFlagDs));

      double efficiencyWeightD = 1.;
      if (applyEfficiency) {
        efficiencyWeightD = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()));
      }

      std::vector<float> outputMl = {-1., -1., -1.};

      if (isDsSignal && isDecayChan && isCorrectInvMassHypo) {
        fillHistoMcRecSig(candidate, multiplicityFT0M);
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
          }
          registry.fill(HIST("hMassDsMcRec"), hfHelper.invMassDsToKKPi(candidate), efficiencyWeightD);
          registry.fill(HIST("hMassDsMcRecSig"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hMassDsVsPtMcRec"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusDsToKKPiMcRec"), candidate.isSelDsToKKPi());
          entryDsCandRecoInfo(hfHelper.invMassDsToKKPi(candidate), candidate.pt(), outputMl[0], outputMl[2]);
          entryDsCandGenInfo(isDsPrompt);
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
          }
          registry.fill(HIST("hMassDsMcRec"), hfHelper.invMassDsToPiKK(candidate), efficiencyWeightD);
          registry.fill(HIST("hMassDsMcRecSig"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hMassDsVsPtMcRec"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusDsToPiKKMcRec"), candidate.isSelDsToPiKK());
          entryDsCandRecoInfo(hfHelper.invMassDsToPiKK(candidate), candidate.pt(), outputMl[0], outputMl[2]);
          entryDsCandGenInfo(isDsPrompt);
        }
      } else {
        fillHistoMcRecBkg(candidate);
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
          }
          registry.fill(HIST("hMassDsMcRec"), hfHelper.invMassDsToKKPi(candidate), efficiencyWeightD);
          registry.fill(HIST("hMassDsMcRecBkg"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hMassDsVsPtMcRec"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusDsToKKPi"), candidate.isSelDsToKKPi());
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
          }
          registry.fill(HIST("hMassDsMcRec"), hfHelper.invMassDsToPiKK(candidate), efficiencyWeightD);
          registry.fill(HIST("hMassDsMcRecBkg"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hMassDsVsPtMcRec"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusDsToPiKK"), candidate.isSelDsToPiKK());
        }
      }

      // Ds-Hadron correlation dedicated section
      // if the candidate is selected as Ds, search for Hadron and evaluate correlations
      for (const auto& track : tracks) {
        // Removing Ds daughters by checking track indices
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
          continue;
        }
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }
        bool isPhysicalPrimary = false;
        int trackOrigin = -1;
        // DsToKKPi and DsToPiKK division
        if (isCorrectInvMassHypo && candidate.isSelDsToKKPi() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(candidate), isDsSignal, isDecayChan);
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          if (track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            isPhysicalPrimary = mcParticle.isPhysicalPrimary();
            trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
            entryDsHadronGenInfo(isDsPrompt, isPhysicalPrimary, trackOrigin);
          } else {
            entryDsHadronGenInfo(isDsPrompt, isPhysicalPrimary, 0);
            registry.fill(HIST("hFakeTracksMcRec"), track.pt());
          }
          // for secondary particle fraction estimation
          if (!isAlreadyFilledEvent) {
            registry.fill(HIST("hPtParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
            if (isPhysicalPrimary) {
              registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
            }
          }
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
        } else if (isCorrectInvMassHypo && candidate.isSelDsToPiKK() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), isDsSignal, isDecayChan);
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          if (track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            isPhysicalPrimary = mcParticle.isPhysicalPrimary();
            trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
            entryDsHadronGenInfo(isDsPrompt, isPhysicalPrimary, trackOrigin);
          } else {
            entryDsHadronGenInfo(isDsPrompt, false, 0);
            registry.fill(HIST("hFakeTracksMcRec"), track.pt());
          }
          // for secondary particle fraction estimation
          if (!isAlreadyFilledEvent) {
            registry.fill(HIST("hPtParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
            if (isPhysicalPrimary) {
              registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
            }
          }
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
        }
      } // end track loop
      isAlreadyFilledEvent = true;
    } // end candidate loop
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRec, "Process MC Reco mode", false);

  /// Ds-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(SelCollisionsMc const& mcCollisions,
                    soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels, aod::McCollisionLabels> const& collisions,
                    CandDsMcGen const& mcParticles)
  {
    BinningTypeMcGen corrBinningMcGen{{zPoolBins, multPoolBins}, true};

    for (const auto& mcCollision : mcCollisions) {

      // auto mcCollision = collision.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();

      int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), mcCollision.multMCFT0A()));
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      registry.fill(HIST("hMultFT0AMcGen"), mcCollision.multMCFT0A());

      const auto groupedMcParticles = mcParticles.sliceBy(perCollisionCandMc, mcCollision.globalIndex());
      const auto groupedCollisions = collisions.sliceBy(collPerCollMc, mcCollision.globalIndex());

      if (groupedCollisions.size() < 1) { // Skipping MC events that have no reconstructed collisions
        continue;
      }
      if (groupedCollisions.size() > 1 && removeCollWSplitVtx) { // Skipping MC events that have more than one reconstructed collision
        continue;
      }

      /// loop over reconstructed collisions
      for (const auto& collision : groupedCollisions) {

        // reco collision selection
        if (useSel8 && !collision.sel8()) {
          continue;
        }
        if (std::abs(collision.posZ()) > 10.) {
          continue;
        }
        if (selNoSameBunchPileUpColl && !(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
          continue;
        }
        if (!collision.has_mcCollision()) {
          registry.fill(HIST("hFakeCollision"), 0.);
          continue;
        }

        bool isDsPrompt = false;
        bool isDecayChan = false;
        int trackOrigin = -1;

        // MC gen level
        for (const auto& particle : groupedMcParticles) {
          // check if the particle is Ds
          if ((std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) && (particle.flagMcDecayChanGen() == decayChannel)) {
            double yD = RecoDecay::y(particle.pVector(), MassDS);
            if (std::abs(yD) > yCandGenMax || particle.pt() < ptCandMin || particle.pt() > ptCandMax) {
              continue;
            }
            fillMcGenHisto(particle);
            // prompt and non-prompt division
            isDsPrompt = particle.originMcGen() == RecoDecay::OriginType::Prompt;
            isDecayChan = particle.flagMcDecayChanGen() == decayChannel;
            std::vector<int> listDaughters{};
            std::array<int, 3> arrDaughDsPDG = {+kKPlus, -kKPlus, kPiPlus};
            std::array<int, 3> prongsId;
            listDaughters.clear();
            RecoDecay::getDaughters(particle, &listDaughters, arrDaughDsPDG, 2);
            int counterDaughters = 0;
            if (listDaughters.size() == 3) {
              for (const auto& dauIdx : listDaughters) {
                // auto daughI = mcParticles.rawIteratorAt(dauIdx - mcParticles.offset());
                auto daughI = groupedMcParticles.rawIteratorAt(dauIdx - groupedMcParticles.offset());
                counterDaughters += 1;
                prongsId[counterDaughters - 1] = daughI.globalIndex();
              }
            }
            // Ds Hadron correlation dedicated section
            for (const auto& particleAssoc : groupedMcParticles) {
              if (std::abs(particleAssoc.eta()) > etaTrackMax || particleAssoc.pt() < ptTrackMin || particleAssoc.pt() > ptTrackMax) {
                continue;
              }
              if (particleAssoc.globalIndex() == prongsId[0] || particleAssoc.globalIndex() == prongsId[1] || particleAssoc.globalIndex() == prongsId[2]) {
                continue;
              }
              if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
                continue;
              }
              if (!particleAssoc.isPhysicalPrimary()) {
                continue;
              }
              // trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
              trackOrigin = RecoDecay::getCharmHadronOrigin(groupedMcParticles, particleAssoc, true);
              registry.fill(HIST("hPtParticleAssocMcGen"), particleAssoc.pt());
              entryDsHadronPair(getDeltaPhi(particleAssoc.phi(), particle.phi()),
                                particleAssoc.eta() - particle.eta(),
                                particle.pt(),
                                particleAssoc.pt(),
                                poolBin);
              entryDsHadronRecoInfo(MassDS, true, isDecayChan);
              entryDsHadronGenInfo(isDsPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
            }
          } // end loop generated particles
        } // end loop generated Ds
      } // end loop reconstructed collision
    }   // end loop generated collision
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcGen, "Process MC Gen mode", false);

  void processDerivedDataDs(SelCollisionsWithDs const& collisions,
                            CandDsData const& candidates,
                            MyTracksData const& tracks)
  {

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDsThisColl = candidates.sliceBy(candsDsPerCollision, thisCollId);
      auto tracksThisColl = tracks.sliceBy(trackIndicesPerCollision, thisCollId);

      // Ds fill histograms and Ds candidates information stored
      for (const auto& candidate : candsDsThisColl) {
        // candidate selected
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          candReduced(hfcReducedCollisionIndex, candidate.phi(), candidate.eta(), candidate.pt(), hfHelper.invMassDsToKKPi(candidate));
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          candReduced(hfcReducedCollisionIndex, candidate.phi(), candidate.eta(), candidate.pt(), hfHelper.invMassDsToPiKK(candidate));
        }
      }

      // tracks information
      for (const auto& track : tracksThisColl) {
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }
        assocTrackReduced(hfcReducedCollisionIndex, track.phi(), track.eta(), track.pt());
      }

      collReduced(collision.multFT0M(), collision.posZ());
      hfcReducedCollisionIndex++;
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processDerivedDataDs, "Process derived data Ds", false);

  void processDerivedDataDsLastIndex(SelCollisionsWithDs const& collisions,
                                     CandDsData const& candidates,
                                     MyTracksData const& tracks)
  {

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDsThisColl = candidates.sliceBy(candsDsPerCollision, thisCollId);
      auto tracksThisColl = tracks.sliceBy(trackIndicesPerCollision, thisCollId);

      int indexHfcReducedCollision = collReduced.lastIndex() + 1;

      // Ds fill histograms and Ds candidates information stored
      for (const auto& candidate : candsDsThisColl) {
        // candidate selected
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          candReduced(indexHfcReducedCollision, candidate.phi(), candidate.eta(), candidate.pt(), hfHelper.invMassDsToKKPi(candidate));
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          candReduced(indexHfcReducedCollision, candidate.phi(), candidate.eta(), candidate.pt(), hfHelper.invMassDsToPiKK(candidate));
        }
      }

      // tracks information
      for (const auto& track : tracksThisColl) {
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }
        assocTrackReduced(indexHfcReducedCollision, track.phi(), track.eta(), track.pt());
      }

      collReduced(collision.multFT0M(), collision.posZ());
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processDerivedDataDsLastIndex, "Process derived data Ds w lastIndex", false);

  // Event Mixing
  void processDataME(SelCollisionsWithDs const& collisions,
                     CandDsData const& candidates,
                     MyTracksData const& tracks)
  {
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    for (const auto& collision : collisions) {
      registry.fill(HIST("hMultFT0M"), collision.multFT0M());
      registry.fill(HIST("hZVtx"), collision.posZ());
    }

    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelCollisionsWithDs, CandDsData, MyTracksData, BinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      if (tracks1.size() == 0) {
        continue;
      }

      // LOGF(info, "Mixed event collisions: Index = (%d, %d), tracks Size: (%d, %d), Z Vertex: (%f, %f), Pool Bin: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), c1.posZ(), c2.posZ(), corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M())), corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M())));
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      int poolBinDs = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hTracksPoolBin"), poolBin); // note that the selections here are not yet applied
      registry.fill(HIST("hDsPoolBin"), poolBinDs);   // note that the selections here are not yet applied
      for (const auto& [cand, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!(cand.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) {
          continue;
        }
        if (std::abs(hfHelper.yDs(cand)) > yCandMax || cand.pt() < ptCandMin || cand.pt() > ptCandMax) {
          continue;
        }
        if (!pAssoc.isGlobalTrackWoDCA()) {
          continue;
        }
        std::vector<float> outputMl = {-1., -1., -1.};
        // DsToKKPi and DsToPiKK division
        if (cand.isSelDsToKKPi() >= selectionFlagDs) {
          // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d), KKPi", cand.index(), pAssoc.index(), c1.index(), c2.index(), cand.collision().index(), pAssoc.collision().index());
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), cand.phi()),
                            pAssoc.eta() - cand.eta(),
                            cand.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(cand), false, false);
          // entryDsHadronGenInfo(false, false, 0);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = cand.mlProbDsToKKPi()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          entryTrackRecoInfo(pAssoc.dcaXY(), pAssoc.dcaZ(), pAssoc.tpcNClsCrossedRows());
        } else if (cand.isSelDsToPiKK() >= selectionFlagDs) {
          // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d), PiKK", cand.index(), pAssoc.index(), c1.index(), c2.index(), cand.collision().index(), pAssoc.collision().index());
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), cand.phi()),
                            pAssoc.eta() - cand.eta(),
                            cand.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(cand), false, false);
          // entryDsHadronGenInfo(false, false, 0);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = cand.mlProbDsToPiKK()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          entryTrackRecoInfo(pAssoc.dcaXY(), pAssoc.dcaZ(), pAssoc.tpcNClsCrossedRows());
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processDataME, "Process Mixed Event Data", false);

  void processMcRecME(SelCollisionsWithDs const& collisions,
                      CandDsMcReco const& candidates,
                      TracksWithMc const& tracks,
                      aod::McParticles const& mcParticles)
  {
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    for (const auto& candidate : candidates) {
      if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        // DsToKKPi and DsToPiKK division
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          fillHistoMcRecSig(candidate, 0.);
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          fillHistoMcRecSig(candidate, 0.);
        }
      } else {
        fillHistoMcRecBkg(candidate);
      }
    }
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelCollisionsWithDs, CandDsMcReco, TracksWithMc, BinningType> pairMcRec{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    bool isDsPrompt = false;
    bool isDsSignal = false;
    bool isDecayChan = false;
    bool isPhysicalPrimary = false;
    int trackOrigin = 0;
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      int poolBinDs = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hMultFT0M"), c1.multFT0M());
      registry.fill(HIST("hZVtx"), c1.posZ());
      registry.fill(HIST("hTracksPoolBin"), poolBin); // note that the selections here are not yet applied
      registry.fill(HIST("hDsPoolBin"), poolBinDs);   // note that the selections here are not yet applied
      for (const auto& [candidate, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
          continue;
        }
        if (!pAssoc.isGlobalTrackWoDCA()) {
          continue;
        }
        std::vector<float> outputMl = {-1., -1., -1.};
        // prompt and non-prompt division
        isDsPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        // Ds Signal
        isDsSignal = std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi;
        isDecayChan = candidate.flagMcDecayChanRec() == decayChannel;
        if (pAssoc.has_mcParticle()) {
          auto mcParticle = pAssoc.template mcParticle_as<aod::McParticles>();
          isPhysicalPrimary = mcParticle.isPhysicalPrimary();
          trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
        } else {
          registry.fill(HIST("hFakeTracksMcRec"), pAssoc.pt());
        }
        // DsToKKPi and DsToPiKK division
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), candidate.phi()),
                            pAssoc.eta() - candidate.eta(),
                            candidate.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(candidate), isDsSignal, isDecayChan);
          entryDsHadronGenInfo(isDsPrompt, isPhysicalPrimary, trackOrigin);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          entryTrackRecoInfo(pAssoc.dcaXY(), pAssoc.dcaZ(), pAssoc.tpcNClsCrossedRows());
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), candidate.phi()),
                            pAssoc.eta() - candidate.eta(),
                            candidate.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), isDsSignal, isDecayChan);
          entryDsHadronGenInfo(isDsPrompt, isPhysicalPrimary, trackOrigin);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          entryTrackRecoInfo(pAssoc.dcaXY(), pAssoc.dcaZ(), pAssoc.tpcNClsCrossedRows());
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRecME, "Process Mixed Event MC Rec", false);

  void processMcGenME(SelCollisionsMc const& collisions,
                      CandDsMcGen const& mcParticles)
  {
    BinningTypeMcGen corrBinningMcGen{{zPoolBins, multPoolBins}, true};
    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<SelCollisionsMc, CandDsMcGen, CandDsMcGen, BinningTypeMcGen> pairMcGen{corrBinningMcGen, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      int poolBin = corrBinningMcGen.getBin(std::make_tuple(c1.posZ(), c1.multMCFT0A()));
      for (const auto& [candidate, particleAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if ((std::abs(candidate.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) && (candidate.flagMcDecayChanGen() == decayChannel)) {
          double yD = RecoDecay::y(candidate.pVector(), MassDS);
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
          bool isDsPrompt = candidate.originMcGen() == RecoDecay::OriginType::Prompt;
          entryDsHadronPair(getDeltaPhi(particleAssoc.phi(), candidate.phi()),
                            particleAssoc.eta() - candidate.eta(),
                            candidate.pt(),
                            particleAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(MassDS, true, true);
          entryDsHadronGenInfo(isDsPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcGenME, "Process Mixed Event MC Gen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDsHadronsSelCollision>(cfgc),
                      adaptAnalysisTask<HfCorrelatorDsHadrons>(cfgc)};
}
