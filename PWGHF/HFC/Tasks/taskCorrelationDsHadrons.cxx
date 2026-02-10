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

/// \file taskCorrelationDsHadrons.cxx
/// \brief Ds-Hadrons azimuthal correlations analysis task - data-like, MC-reco and MC-Gen analyses
/// \author Grazia Luparello <grazia.luparello@cern.ch>
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#include "PWGHF/Core/CentralityEstimation.h"
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
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
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
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TPDGCode.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::analysis::hf_correlations;

enum ResonantChannel : int8_t {
  PhiPi = 1,
  Kstar0K = 2
};

namespace
{
std::unordered_map<int8_t, int8_t> channelsResonant = {{{ResonantChannel::PhiPi, hf_decay::hf_cand_3prong::DecayChannelResonant::DsToPhiPi}, {ResonantChannel::Kstar0K, hf_decay::hf_cand_3prong::DecayChannelResonant::DsToKstar0K}}};
}

/// Ds-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfTaskCorrelationDsHadrons {
  Configurable<bool> fillHistoData{"fillHistoData", true, "Flag for filling histograms in data processes"};
  Configurable<bool> fillHistoMcRec{"fillHistoMcRec", true, "Flag for filling histograms in MC Rec processes"};
  Configurable<bool> fillHistoMcGen{"fillHistoMcGen", true, "Flag for filling histograms in MC Gen processes"};
  Configurable<bool> fillHistoMcEff{"fillHistoMcEff", true, "Flag for filling histograms in efficiency processes"};
  Configurable<bool> applyEfficiency{"applyEfficiency", false, "Flag for applying efficiency weights"};
  Configurable<bool> useSel8ForEff{"useSel8ForEff", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<bool> removeCollWSplitVtx{"removeCollWSplitVtx", false, "Flag for rejecting the splitted collisions"};
  Configurable<bool> loadAccXEffFromCCDB{"loadAccXEffFromCCDB", false, "Flag for loading efficiency distributions from CCDB"};
  Configurable<bool> separateTrackOrigins{"separateTrackOrigins", false, "Flag to enable separation of track origins (from c or b)"};
  Configurable<bool> useHighDimHistoForEff{"useHighDimHistoForEff", false, "Flag to create/use higher dimension histograms in the efficiency processes/correction"};
  Configurable<bool> applyDeltaPhiCorrEff{"applyDeltaPhiCorrEff", false, "Flag to use higher dimension histograms with delta phi correction in the efficiency correction"};
  Configurable<bool> doLSpair{"doLSpair", false, "Flag to evaluate correlations for like-sign pairs"};
  Configurable<bool> doULSpair{"doULSpair", false, "Flag to evaluate correlations for unlike-sign pairs"};
  Configurable<bool> pidTrkApplied{"pidTrkApplied", false, "Apply PID selection for associated tracks"};
  Configurable<bool> forceTOF{"forceTOF", false, "force the TOF signal for the PID"};
  // Configurable<bool> doMcCollisionCheck{"doMcCollisionCheck", false, "Flag for applying the collision check and selection based on MC collision info"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds (avoid the case of flag = 0, no outputMlScore)"};
  Configurable<int> nTpcCrossedRaws{"nTpcCrossedRaws", 70, "Number of crossed TPC Rows"};
  // Configurable<int> eventGeneratorType{"eventGeneratorType", -1, "If positive, enable event selection using subGeneratorId information. The value indicates which events to keep (0 = MB, 4 = charm triggered, 5 = beauty triggered)"};
  Configurable<int> decayChannel{"decayChannel", 1, "Resonant decay channels: 1 for Ds->PhiPi->KKpi, 2 for Ds->K0*K->KKPi"};
  Configurable<float> cutCollPosZMc{"cutCollPosZMc", 10., "max z-vertex position for collision acceptance"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT (used in eff. process only)"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand pT (used in eff. process only)"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT (used in eff. process only)"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT (used in eff. process only)"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity (used in eff. process only)"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen. cand. rapidity (used in eff. process only)"};
  Configurable<float> ptDaughterMin{"ptDaughterMin", 0.1, "min. daughter pT (used in eff. process only)"};
  Configurable<float> tofPIDThreshold{"tofPIDThreshold", 0.75, "minimum pT after which TOF PID is applicable"};
  Configurable<float> centralityMin{"centralityMin", 0., "min. centrality"};
  Configurable<float> centralityMax{"centralityMax", 100., "max. centrality"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Proton, o2::track::PID::Pion, o2::track::PID::Kaon}, "Trk sel: Particles species for PID, proton, pion, kaon"};
  Configurable<std::vector<float>> pidTPCMax{"pidTPCMax", std::vector<float>{3., 0., 0.}, "maximum nSigma TPC"};
  Configurable<std::vector<float>> pidTOFMax{"pidTOFMax", std::vector<float>{3., 0., 0.}, "maximum nSigma TOF"};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle efficiency"};
  Configurable<std::vector<double>> mlOutputPromptMin{"mlOutputPromptMin", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for prompt"};
  Configurable<std::vector<double>> mlOutputPromptMax{"mlOutputPromptMax", {1.0, 1.0, 1.0, 1.0}, "Machine learning scores for prompt"};
  Configurable<std::vector<double>> mlOutputBkg{"mlOutputBkg", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for bkg"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for D-meson efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyHad{"binsPtEfficiencyHad", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for associated particle efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", {1., 1., 1., 1., 1., 1.}, "efficiency values for Ds meson"};
  Configurable<std::vector<double>> efficiencyHad{"efficiencyHad", {1., 1., 1., 1., 1., 1.}, "efficiency values for associated particles"};
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", {1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", {1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", {1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", {1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", {2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", {2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320}, "Outer values of right sideband vs pT"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> associatedEffCcdbPath{"associatedEffCcdbPath", "", "CCDB path for associated efficiency"};
  Configurable<std::string> promptEffCcdbPath{"promptEffCcdbPath", "", "CCDB path for trigger efficiency"};
  Configurable<std::string> fdEffCcdbPath{"fdEffCcdbPath", "", "CCDB path for trigger efficiency"};
  Configurable<int64_t> timestampCcdb{"timestampCcdb", -1, "timestamp of the efficiency files used to query in CCDB"};

  std::shared_ptr<TH1> hEfficiencyD = nullptr;
  std::shared_ptr<TH1> hEfficiencyAssociated = nullptr;
  std::shared_ptr<TH2> hEfficiencyDMult = nullptr;
  std::shared_ptr<TH2> hEfficiencyAssociatedMult = nullptr;
  std::shared_ptr<TH3> hEfficiencyAssociatedDeltaPhiCorr = nullptr;

  static constexpr float Epsilon{1.e-8};

  SliceCache cache;

  Service<ccdb::BasicCCDBManager> ccdb{};

  using DsHadronPair = soa::Filtered<soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo>>;
  using DsHadronPairFull = soa::Filtered<soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronGenInfo>>;
  using DsHadronPairWithMl = soa::Filtered<soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronMlInfo, aod::TrackRecoInfo>>;
  using DsHadronPairFullWithMl = soa::Filtered<soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronGenInfo, aod::DsHadronMlInfo, aod::TrackRecoInfo>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi, aod::HfCand3ProngMcRec>>;                                                                                                                 // flagDsFilter applied
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;                                                                                                                                                                         // flagDsFilter applied
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>; // trackFilter applied

  Filter flagDsFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);
  Filter trackPairFilter = (nabs(aod::hf_correlation_ds_hadron::signedPtHadron) > ptTrackMin) && (nabs(aod::hf_correlation_ds_hadron::signedPtHadron) < ptTrackMax);

  Preslice<CandDsMcReco> perCollisionCand = o2::aod::hf_cand::collisionId;
  Preslice<CandDsMcGen> perCollisionCandMc = o2::aod::mcparticle::mcCollisionId;
  Preslice<TracksWithMc> perCollision = o2::aod::track::collisionId;
  Preslice<o2::aod::McParticles> perCollisionMc = o2::aod::mcparticle::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels, aod::McCollisionLabels>> collPerCollMc = o2::aod::mccollisionlabel::mcCollisionId;

  ConfigurableAxis binsMassD{"binsMassD", {200, 1.7, 2.25}, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {100, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 8000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsNumPvContr{"binsNumPvContr", {100, 0., 1000.}, "number PV contributors"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMassD = {binsMassD, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
    AxisSpec axisDetlaEta = {binsEta, "#it{#eta}^{Hadron}-#it{#eta}^{D}"};
    AxisSpec axisEta = {binsEta, "#it{#eta}"};
    AxisSpec axisDetlaPhi = {binsPhi, "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad)"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T}^{D} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T}^{Hadron} (GeV/#it{c})"};
    AxisSpec axisPoolBin = {binsPoolBin, "poolBin"};
    AxisSpec const axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec const axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec axisNumPvContr = {binsNumPvContr, "Num PV contributors"};
    AxisSpec axisDsPrompt = {2, -0.5, 1.5, "Prompt Ds"};

    // Histograms for data analysis
    registry.add("hBdtScorePrompt", "Ds BDT prompt score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hBdtScoreBkg", "Ds BDT bkg score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hMassDsVsPt", "Ds candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    if (doLSpair) {
      registry.add("hCorrel2DVsPtSignalRegionLS", "Ds-h correlations signal region - LS pairs", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSidebandLeftLS", "Ds-h correlations sideband left region - LS pairs", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSidebandRightLS", "Ds-h correlations sideband right region - LS pairs", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});

      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionLS"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandLeftLS"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandRightLS"))->Sumw2();
    }
    if (doULSpair) {
      registry.add("hCorrel2DVsPtSignalRegionULS", "Ds-h correlations signal region - ULS pairs", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSidebandLeftULS", "Ds-h correlations sideband left region - ULS pairs", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSidebandRightULS", "Ds-h correlations sideband right region - ULS pairs", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});

      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionULS"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandLeftULS"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandRightULS"))->Sumw2();
    }
    if (fillHistoData) {
      registry.add("hDeltaEtaPtIntSignalRegion", "Ds-h deltaEta signal region", {HistType::kTH1F, {axisDetlaEta}});
      registry.add("hDeltaPhiPtIntSignalRegion", "Ds-h deltaPhi signal region", {HistType::kTH1F, {axisDetlaPhi}});
      registry.add("hCorrel2DVsPtSignalRegion", "Ds-h correlations signal region", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hDeltaEtaPtIntSidebandLeft", "Ds-h deltaEta sideband left region", {HistType::kTH1F, {axisDetlaEta}});
      registry.add("hDeltaPhiPtIntSidebandLeft", "Ds-h deltaPhi sideband left region", {HistType::kTH1F, {axisDetlaPhi}});
      registry.add("hDeltaEtaPtIntSidebandRight", "Ds-h deltaEta sideband right region", {HistType::kTH1F, {axisDetlaEta}});
      registry.add("hDeltaPhiPtIntSidebandRight", "Ds-h deltaPhi sideband right region", {HistType::kTH1F, {axisDetlaPhi}});
      registry.add("hCorrel2DVsPtSidebandLeft", "Ds-h correlations sideband left region", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSidebandRight", "Ds-h correlations sideband right region", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});

      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandLeft"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandRight"))->Sumw2();
    }
    // Histograms for MC Reco analysis
    if (fillHistoMcRec) {
      registry.add("hMassPromptDsVsPt", "Ds prompt candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
      registry.add("hMassNonPromptDsVsPt", "Ds non prompt candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
      registry.add("hDeltaEtaPtIntSignalRegionMcRec", "Ds-h deltaEta signal region MC reco", {HistType::kTH1F, {axisDetlaEta}});
      registry.add("hDeltaPhiPtIntSignalRegionMcRec", "Ds-h deltaPhi signal region MC reco", {HistType::kTH1F, {axisDetlaPhi}});
      registry.add("hCorrel2DVsPtSignalRegionMcRec", "Ds-h correlations signal region Prompt-NonPrompt MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisDsPrompt}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSignalRegionPromptDsPromptHadronMcRec", "Ds-h correlations signal region Prompt-NonPrompt MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSignalRegionNonPromptDsNonPromptHadronMcRec", "Ds-h correlations signal region Prompt-NonPrompt MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtPhysicalPrimaryMcRec", "Ds-h correlations signal region (only true primary particles) MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisDsPrompt}, {axisPoolBin}}});
      registry.add("hDeltaEtaPtIntSidebandLeftMcRec", "Ds-h deltaEta sideband left region MC reco", {HistType::kTH1F, {axisDetlaEta}});
      registry.add("hDeltaPhiPtIntSidebandLeftMcRec", "Ds-h deltaPhi sideband left region MC reco", {HistType::kTH1F, {axisDetlaPhi}});
      registry.add("hCorrel2DVsPtSidebandLeftMcRec", "Ds-h correlations sideband left region MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hDeltaEtaPtIntSidebandRightMcRec", "Ds-h deltaEta sideband right region MC reco", {HistType::kTH1F, {axisDetlaEta}});
      registry.add("hDeltaPhiPtIntSidebandRightMcRec", "Ds-h deltaPhi sideband right region MC reco", {HistType::kTH1F, {axisDetlaPhi}});
      registry.add("hCorrel2DVsPtSidebandRightMcRec", "Ds-h correlations sideband right region MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      if (doLSpair) {
        registry.add("hCorrel2DVsPtPhysicalPrimaryMcRecLS", "Ds-h correlations signal region (only true primary particles) MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisDsPrompt}, {axisPoolBin}}});
      }
      if (doULSpair) {
        registry.add("hCorrel2DVsPtPhysicalPrimaryMcRecULS", "Ds-h correlations signal region (only true primary particles) MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisDsPrompt}, {axisPoolBin}}});
      }
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandLeftMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandRightMcRec"))->Sumw2();
    }
    // Histograms for MC Gen analysis
    if (fillHistoMcGen) {
      registry.add("hDeltaEtaPtIntMcGen", "Ds-h deltaEta MC Gen", {HistType::kTH1F, {axisDetlaEta}});
      registry.add("hDeltaPhiPtIntMcGen", "Ds-h deltaPhi MC Gen", {HistType::kTH1F, {axisDetlaPhi}});
      registry.add("hCorrel2DVsPtMcGen", "Ds-h correlations MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtMcGenPrompt", "Ds-h correlations Prompt MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtMcGenPromptDsPromptHadron", "Ds-h correlations prompt Ds prompt h MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtMcGenNonPromptDsNonPromptHadron", "Ds-h correlations non prompt Ds non prompt h MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtMcGenNonPrompt", "Ds-h correlations NonPrompt MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});

      if (doLSpair) {
        registry.add("hCorrel2DVsPtMcGenPromptLS", "Ds-h correlations Prompt MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtMcGenNonPromptLS", "Ds-h correlations NonPrompt MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      }
      if (doULSpair) {
        registry.add("hCorrel2DVsPtMcGenPromptULS", "Ds-h correlations Prompt MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtMcGenNonPromptULS", "Ds-h correlations NonPrompt MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      }
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenPrompt"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenNonPrompt"))->Sumw2();
    }
    // Histograms for efficiencies
    if (fillHistoMcEff) {
      registry.add("hFakeCollision", "Fake collision counter", {HistType::kTH1F, {{1, -0.5, 0.5, "n fake coll"}}});
      registry.add("hFakeTracks", "Fake tracks counter", {HistType::kTH1F, {{1, -0.5, 0.5, "n fake tracks"}}});
      registry.add("hNumPvContrib", "Num PV contributors", {HistType::kTH1F, {axisNumPvContr}});
      registry.add("hPtPrmPromptPartMcGen", "Primary prompt particles - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmNonPromptPartMcGen", "Primary non-prompt particles - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmPromptPartMcRec", "Primary prompt particles - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmNonPromptPartMcRec", "Primary non-prompt particles - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtMcParticleAssocSpecieMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtParticleAssocMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtCandMcGenDaughterInAcc", "Ds,Hadron particles non prompt - MC Gen", {HistType::kTH1F, {axisPtD}});

      if (useHighDimHistoForEff) {
        registry.add("hPtCandMcRecPrompt", "Ds prompt candidates", {HistType::kTH2F, {{axisPtD}, {axisNumPvContr}}});
        registry.add("hPtCandMcRecNonPrompt", "Ds non prompt candidates pt", {HistType::kTH2F, {{axisPtD}, {axisNumPvContr}}});
        registry.add("hPtCandMcGenPrompt", "Ds,Hadron particles prompt - MC Gen", {HistType::kTH2F, {{axisPtD}, {axisNumPvContr}}});
        registry.add("hPtCandMcGenNonPrompt", "Ds,Hadron particles non prompt - MC Gen", {HistType::kTH2F, {{axisPtD}, {axisNumPvContr}}});
        registry.add("hPtParticleAssocSpecieMcRec", "Associated Particle - MC Rec", {HistType::kTHnSparseF, {axisPtHadron}});
        registry.add("hPtPrmPionMcRec", "Primary pions - MC Rec", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmKaonMcRec", "Primary kaons - MC Rec", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmProtonMcRec", "Primary protons - MC Rec", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmElectronMcRec", "Primary electrons - MC Rec", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmMuonMcRec", "Primary muons - MC Rec", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmPionMcGen", "Primary pions - MC Gen", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmKaonMcGen", "Primary kaons - MC Gen", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmProtonMcGen", "Primary protons - MC Gen", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmElectronMcGen", "Primary electrons - MC Gen", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
        registry.add("hPtPrmMuonMcGen", "Primary muons - MC Gen", {HistType::kTHnSparseF, {{axisPtHadron}, {axisEta}, {axisPosZ}, {axisNumPvContr}}});
      } else {
        registry.add("hPtCandMcRecPrompt", "Ds prompt candidates pt", {HistType::kTH1F, {axisPtD}});
        registry.add("hPtCandMcRecNonPrompt", "Ds non prompt candidates pt", {HistType::kTH1F, {axisPtD}});
        registry.add("hPtCandMcGenPrompt", "Ds,Hadron particles prompt - MC Gen", {HistType::kTH1F, {axisPtD}});
        registry.add("hPtCandMcGenNonPrompt", "Ds,Hadron particles non prompt - MC Gen", {HistType::kTH1F, {axisPtD}});
        registry.add("hPtParticleAssocSpecieMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmPionMcRec", "Primary pions - MC Rec", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmKaonMcRec", "Primary kaons - MC Rec", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmProtonMcRec", "Primary protons - MC Rec", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmElectronMcRec", "Primary electrons - MC Rec", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmMuonMcRec", "Primary muons - MC Rec", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmPionMcGen", "Primary pions - MC Gen", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmKaonMcGen", "Primary kaons - MC Gen", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmProtonMcGen", "Primary protons - MC Gen", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmElectronMcGen", "Primary electrons - MC Gen", {HistType::kTH1F, {axisPtHadron}});
        registry.add("hPtPrmMuonMcGen", "Primary muons - MC Gen", {HistType::kTH1F, {axisPtHadron}});
      }
    }

    // Loading efficiency histograms from CCDB
    if (applyEfficiency && loadAccXEffFromCCDB) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

      if (useHighDimHistoForEff) {
        hEfficiencyDMult = std::shared_ptr<TH2>(ccdb->getForTimeStamp<TH2F>(promptEffCcdbPath, timestampCcdb));
        if (hEfficiencyDMult == nullptr) {
          LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", promptEffCcdbPath.value.c_str());
        }
        LOGF(info, "Loaded trigger efficiency (prompt D) histogram from %s", promptEffCcdbPath.value.c_str());

        if (applyDeltaPhiCorrEff) {
          hEfficiencyAssociatedDeltaPhiCorr = std::shared_ptr<TH3>(ccdb->getForTimeStamp<TH3F>(associatedEffCcdbPath, timestampCcdb));
          if (hEfficiencyAssociatedDeltaPhiCorr == nullptr) {
            LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", associatedEffCcdbPath.value.c_str());
          }
          LOGF(info, "Loaded associated efficiency histogram from %s", associatedEffCcdbPath.value.c_str());
        } else {
          hEfficiencyAssociatedMult = std::shared_ptr<TH2>(ccdb->getForTimeStamp<TH2F>(associatedEffCcdbPath, timestampCcdb));
          if (hEfficiencyAssociatedMult == nullptr) {
            LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", associatedEffCcdbPath.value.c_str());
          }
          LOGF(info, "Loaded associated efficiency histogram from %s", associatedEffCcdbPath.value.c_str());
        }
      } else {
        hEfficiencyD = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1F>(promptEffCcdbPath, timestampCcdb));
        if (hEfficiencyD == nullptr) {
          LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", promptEffCcdbPath.value.c_str());
        }
        LOGF(info, "Loaded trigger efficiency (prompt D) histogram from %s", promptEffCcdbPath.value.c_str());

        hEfficiencyAssociated = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1F>(associatedEffCcdbPath, timestampCcdb));
        if (hEfficiencyAssociated == nullptr) {
          LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", associatedEffCcdbPath.value.c_str());
        }
        LOGF(info, "Loaded associated efficiency histogram from %s", associatedEffCcdbPath.value.c_str());
      }
    }
  }

  enum class EfficiencyMode {
    DsOnly = 1,
    DsHadronPair = 2
  };

  bool isSelectedCandidate(const int ptBinD, const float bdtScorePrompt, const float bdtScoreBkg)
  {

    return (ptBinD != -1 && bdtScorePrompt >= mlOutputPromptMin->at(ptBinD) && bdtScorePrompt <= mlOutputPromptMax->at(ptBinD) && bdtScoreBkg <= mlOutputBkg->at(ptBinD));
  }

  double getEfficiencyWeight(float ptD, std::optional<int> multPvContrib = std::nullopt, std::optional<float> ptAssoc = std::nullopt, std::optional<float> deltaPhi = std::nullopt, EfficiencyMode mode = EfficiencyMode::DsOnly)
  {
    if (!applyEfficiency) {
      return 1.;
    }

    double weight = 1.;

    switch (mode) {
      case EfficiencyMode::DsOnly:
        if (loadAccXEffFromCCDB) {
          if (useHighDimHistoForEff) {
            if (hEfficiencyDMult->GetBinContent(hEfficiencyDMult->FindBin(ptD, static_cast<double>(*multPvContrib))) <= Epsilon) {
              LOG(fatal) << "A bin content in Ds-meson efficiency histogram is zero!";
            }
            weight = 1. / hEfficiencyDMult->GetBinContent(hEfficiencyDMult->FindBin(ptD, static_cast<double>(*multPvContrib)));
          } else {
            if (hEfficiencyD->GetBinContent(hEfficiencyD->FindBin(ptD)) <= Epsilon) {
              LOG(fatal) << "A bin content in Ds-meson efficiency histogram is zero!";
            }
            weight = 1. / hEfficiencyD->GetBinContent(hEfficiencyD->FindBin(ptD));
          }
        } else {
          if (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) <= Epsilon) {
            LOG(fatal) << "A bin content in Ds-meson efficiency vector is zero!";
          }
          weight = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD));
        }
        break;
      case EfficiencyMode::DsHadronPair:
        if (loadAccXEffFromCCDB) {
          if (ptAssoc && (hEfficiencyAssociated || hEfficiencyAssociatedMult || hEfficiencyAssociatedDeltaPhiCorr)) {
            if (useHighDimHistoForEff) {
              if (applyDeltaPhiCorrEff) {
                if (hEfficiencyAssociatedDeltaPhiCorr->GetBinContent(hEfficiencyAssociatedDeltaPhiCorr->FindBin(*ptAssoc, ptD, static_cast<double>(*deltaPhi))) <= Epsilon) {
                  LOG(fatal) << "A bin content in associated particle efficiency histogram is zero!";
                }
                weight = 1. / (hEfficiencyDMult->GetBinContent(hEfficiencyDMult->FindBin(ptD, static_cast<double>(*multPvContrib))) * hEfficiencyAssociatedDeltaPhiCorr->GetBinContent(hEfficiencyAssociatedDeltaPhiCorr->FindBin(*ptAssoc, ptD, static_cast<double>(*deltaPhi))));
              } else {
                if (hEfficiencyAssociatedMult->GetBinContent(hEfficiencyAssociatedMult->FindBin(*ptAssoc, static_cast<double>(*multPvContrib))) <= Epsilon) {
                  LOG(fatal) << "A bin content in associated particle efficiency histogram is zero!";
                }
                weight = 1. / (hEfficiencyDMult->GetBinContent(hEfficiencyDMult->FindBin(ptD, static_cast<double>(*multPvContrib))) * hEfficiencyAssociatedMult->GetBinContent(hEfficiencyAssociatedMult->FindBin(*ptAssoc, static_cast<double>(*multPvContrib))));
              }
            } else {
              if (hEfficiencyAssociated->GetBinContent(hEfficiencyAssociated->FindBin(*ptAssoc)) <= Epsilon) {
                LOG(fatal) << "A bin content in associated particle efficiency histogram is zero!";
              }
              weight = 1. / (hEfficiencyD->GetBinContent(hEfficiencyD->FindBin(ptD)) * hEfficiencyAssociated->GetBinContent(hEfficiencyAssociated->FindBin(*ptAssoc)));
            }
          }
        } else {
          if (ptAssoc) {
            if (efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, *ptAssoc)) <= Epsilon) {
              LOG(fatal) << "A bin content in associated particle efficiency vector is zero!";
            }
            weight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, *ptAssoc)));
          }
        }
        break;
      default:
        LOG(fatal) << "Unknown efficiency mode!";
        break;
    }

    return weight;
  }

  /// Check event selections for collision and fill the collision table
  /// \param collision is the collision
  template <typename Coll>
  int getCentrality(Coll const& collision)
  {
    int cent{-1};
    if (centEstimator == CentralityEstimator::FT0A) {
      cent = collision.centFT0A();
    } else if (centEstimator == CentralityEstimator::FT0C) {
      cent = collision.centFT0C();
    } else if (centEstimator == CentralityEstimator::FT0M) {
      cent = collision.centFT0M();
    } else if (centEstimator == CentralityEstimator::FV0A) {
      cent = collision.centFV0A();
    } else {
      LOG(fatal) << "Centrality estimator not recognized for collision selection";
      std::abort();
    }
    return cent;
  }

  void processData(DsHadronPairWithMl const& pairEntries,
                   aod::DsCandRecoInfo const& candidates)
  {
    for (const auto& candidate : candidates) {
      float const massD = candidate.mD();
      float const ptD = candidate.signedPtD();
      float const bdtScorePrompt = candidate.mlScorePrompt();
      float const bdtScoreBkg = candidate.mlScoreBkg();
      int multPvContrib = candidate.numPvContrib();
      int const ptBinD = o2::analysis::findBin(binsPtD, std::abs(ptD));

      if (!isSelectedCandidate(ptBinD, bdtScorePrompt, bdtScoreBkg)) {
        continue;
      }

      double efficiencyWeightD = 1.;
      if (useHighDimHistoForEff) {
        efficiencyWeightD = getEfficiencyWeight(std::abs(ptD), multPvContrib);
      } else {
        efficiencyWeightD = getEfficiencyWeight(std::abs(ptD));
      }

      registry.fill(HIST("hMassDsVsPt"), massD, std::abs(ptD), efficiencyWeightD);
      registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
      registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.signedPtD();
      float const ptHadron = pairEntry.signedPtHadron();
      float const massD = pairEntry.mD();
      float const bdtScorePrompt = pairEntry.mlScorePrompt();
      float const bdtScoreBkg = pairEntry.mlScoreBkg();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int multPvContrib = pairEntry.numPvContrib();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int const poolBin = pairEntry.poolBin();
      int const ptBinD = o2::analysis::findBin(binsPtD, std::abs(ptD));

      if (!isSelectedCandidate(ptBinD, bdtScorePrompt, bdtScoreBkg)) {
        continue;
      }

      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (useHighDimHistoForEff) {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), multPvContrib, std::abs(ptHadron), deltaPhi, EfficiencyMode::DsHadronPair);
      } else {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), std::nullopt, std::abs(ptHadron), std::nullopt, EfficiencyMode::DsHadronPair);
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSignalRegionLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSignalRegionULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
        }
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandLeftLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandLeftULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSidebandLeft"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSidebandLeft"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSidebandLeft"), deltaPhi, efficiencyWeight);
        }
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandRightLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandRightULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSidebandRight"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSidebandRight"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSidebandRight"), deltaPhi, efficiencyWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processData, "Process data", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(DsHadronPairFullWithMl const& pairEntries,
                    soa::Join<aod::DsCandRecoInfo, aod::DsCandGenInfo> const& candidates)
  {
    for (const auto& candidate : candidates) {
      float const massD = candidate.mD();
      float const ptD = candidate.signedPtD();
      float const bdtScorePrompt = candidate.mlScorePrompt();
      float const bdtScoreBkg = candidate.mlScoreBkg();
      int const ptBinD = o2::analysis::findBin(binsPtD, std::abs(ptD));
      int multPvContrib = candidate.numPvContrib();
      bool const isDsPrompt = candidate.isPrompt();

      if (!isSelectedCandidate(ptBinD, bdtScorePrompt, bdtScoreBkg)) {
        continue;
      }

      double efficiencyWeightD = 1.;
      if (useHighDimHistoForEff) {
        efficiencyWeightD = getEfficiencyWeight(std::abs(ptD), multPvContrib);
      } else {
        efficiencyWeightD = getEfficiencyWeight(std::abs(ptD));
      }

      if (isDsPrompt) {
        registry.fill(HIST("hMassPromptDsVsPt"), massD, std::abs(ptD), efficiencyWeightD);
        registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
        registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
      } else {
        registry.fill(HIST("hMassNonPromptDsVsPt"), massD, std::abs(ptD), efficiencyWeightD);
        registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
        registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
      }
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.signedPtD();
      float const ptHadron = pairEntry.signedPtHadron();
      float const massD = pairEntry.mD();
      float const bdtScorePrompt = pairEntry.mlScorePrompt();
      float const bdtScoreBkg = pairEntry.mlScoreBkg();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int multPvContrib = pairEntry.numPvContrib();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int const poolBin = pairEntry.poolBin();
      int const statusDsPrompt = static_cast<int>(pairEntry.isPrompt());
      int const statusPromptHadron = pairEntry.trackOrigin();
      int const ptBinD = o2::analysis::findBin(binsPtD, std::abs(ptD));
      bool const isPhysicalPrimary = pairEntry.isPhysicalPrimary();

      if (!isSelectedCandidate(ptBinD, bdtScorePrompt, bdtScoreBkg)) {
        continue;
      }

      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (useHighDimHistoForEff) {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), multPvContrib, std::abs(ptHadron), deltaPhi, EfficiencyMode::DsHadronPair);
      } else {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), std::nullopt, std::abs(ptHadron), std::nullopt, EfficiencyMode::DsHadronPair);
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        // prompt and non-prompt division
        if (pairEntry.isSignal() && pairEntry.isDecayChan()) {
          registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
          registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), statusDsPrompt, poolBin, efficiencyWeight);
          if (isPhysicalPrimary) {
            if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
              registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRecLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), statusDsPrompt, poolBin, efficiencyWeight);
            } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
              registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRecULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), statusDsPrompt, poolBin, efficiencyWeight);
            } else { // default case
              registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), statusDsPrompt, poolBin, efficiencyWeight);
            }
            if (statusDsPrompt == 1 && statusPromptHadron == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionPromptDsPromptHadronMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
            } else if (statusDsPrompt == 0 && statusPromptHadron == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionNonPromptDsNonPromptHadronMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
            }
          }
        }
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeftMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeftMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeftMcRec"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRightMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandRightMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRightMcRec"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcRec, "Process MC Reco mode", false);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(DsHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float const deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.signedPtD();
      float const ptHadron = pairEntry.signedPtHadron();
      int const poolBin = pairEntry.poolBin();
      int const statusPromptHadron = pairEntry.trackOrigin();
      bool const isDsPrompt = pairEntry.isPrompt();

      registry.fill(HIST("hCorrel2DVsPtMcGen"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
      registry.fill(HIST("hDeltaEtaPtIntMcGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMcGen"), deltaPhi);
      if (isDsPrompt) {
        registry.fill(HIST("hCorrel2DVsPtMcGenPrompt"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
        if (doULSpair) {
          registry.fill(HIST("hCorrel2DVsPtMcGenPromptULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
        }
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) {
          registry.fill(HIST("hCorrel2DVsPtMcGenPromptLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
        }
        if (statusPromptHadron == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hCorrel2DVsPtMcGenPromptDsPromptHadron"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
        }
      } else {
        registry.fill(HIST("hCorrel2DVsPtMcGenNonPrompt"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
        if (doULSpair) {
          registry.fill(HIST("hCorrel2DVsPtMcGenNonPromptULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
        }
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) {
          registry.fill(HIST("hCorrel2DVsPtMcGenNonPromptLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
        }
        if (statusPromptHadron == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("hCorrel2DVsPtMcGenNonPromptDsNonPromptHadron"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcGen, "Process MC Gen mode", false);

  void processDataME(DsHadronPairWithMl const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.signedPtD();
      float const ptHadron = pairEntry.signedPtHadron();
      float const massD = pairEntry.mD();
      float const bdtScorePrompt = pairEntry.mlScorePrompt();
      float const bdtScoreBkg = pairEntry.mlScoreBkg();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int multPvContrib = pairEntry.numPvContrib();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int const poolBin = pairEntry.poolBin();
      int const ptBinD = o2::analysis::findBin(binsPtD, std::abs(ptD));

      if (!isSelectedCandidate(ptBinD, bdtScorePrompt, bdtScoreBkg)) {
        continue;
      }

      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (useHighDimHistoForEff) {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), multPvContrib, std::abs(ptHadron), deltaPhi, EfficiencyMode::DsHadronPair);
      } else {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), std::nullopt, std::abs(ptHadron), std::nullopt, EfficiencyMode::DsHadronPair);
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSignalRegionLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSignalRegionULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
        }
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandLeftLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandLeftULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSidebandLeft"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSidebandLeft"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSidebandLeft"), deltaPhi, efficiencyWeight);
        }
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandRightLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandRightULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSidebandRight"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSidebandRight"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSidebandRight"), deltaPhi, efficiencyWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processDataME, "Process data ME", false);

  void processDerivedDataME(DsHadronPair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.signedPtD();
      float const ptHadron = pairEntry.signedPtHadron();
      float const massD = pairEntry.mD();
      int multPvContrib = pairEntry.numPvContrib();
      int const poolBin = pairEntry.poolBin();
      int const ptBinD = o2::analysis::findBin(binsPtD, std::abs(ptD));

      double efficiencyWeight = 1.;
      if (useHighDimHistoForEff) {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), multPvContrib, std::abs(ptHadron), deltaPhi, EfficiencyMode::DsHadronPair);
      } else {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), std::nullopt, std::abs(ptHadron), std::nullopt, EfficiencyMode::DsHadronPair);
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSignalRegionLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSignalRegionULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
        }
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandLeftLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandLeftULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSidebandLeft"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSidebandLeft"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSidebandLeft"), deltaPhi, efficiencyWeight);
        }
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandRightLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
          registry.fill(HIST("hCorrel2DVsPtSidebandRightULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        } else if (fillHistoData) { // default case
          registry.fill(HIST("hCorrel2DVsPtSidebandRight"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
          registry.fill(HIST("hDeltaEtaPtIntSidebandRight"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSidebandRight"), deltaPhi, efficiencyWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processDerivedDataME, "Process derived data ME", false);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRecME(DsHadronPairFullWithMl const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.signedPtD();
      float const ptHadron = pairEntry.signedPtHadron();
      float const massD = pairEntry.mD();
      float const bdtScorePrompt = pairEntry.mlScorePrompt();
      float const bdtScoreBkg = pairEntry.mlScoreBkg();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int multPvContrib = pairEntry.numPvContrib();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int const poolBin = pairEntry.poolBin();
      int const statusDsPrompt = static_cast<int>(pairEntry.isPrompt());
      int const statusPromptHadron = pairEntry.trackOrigin();
      int const ptBinD = o2::analysis::findBin(binsPtD, std::abs(ptD));
      bool const isPhysicalPrimary = pairEntry.isPhysicalPrimary();

      if (!isSelectedCandidate(ptBinD, bdtScorePrompt, bdtScoreBkg)) {
        continue;
      }

      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (useHighDimHistoForEff) {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), multPvContrib, std::abs(ptHadron), deltaPhi, EfficiencyMode::DsHadronPair);
      } else {
        efficiencyWeight = getEfficiencyWeight(std::abs(ptD), std::nullopt, std::abs(ptHadron), std::nullopt, EfficiencyMode::DsHadronPair);
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        // prompt and non-prompt division
        if (pairEntry.isSignal() && pairEntry.isDecayChan()) {
          registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
          registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), statusDsPrompt, poolBin, efficiencyWeight);
          if (isPhysicalPrimary) {
            if (doLSpair && ((ptD > 0. && ptHadron > 0.) || (ptD < 0. && ptHadron < 0.))) { // like-sign pairs
              registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRecLS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), statusDsPrompt, poolBin, efficiencyWeight);
            } else if (doULSpair && ((ptD > 0. && ptHadron < 0.) || (ptD < 0. && ptHadron > 0.))) { // unlike-sign pairs
              registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRecULS"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), statusDsPrompt, poolBin, efficiencyWeight);
            } else { // default case
              registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), statusDsPrompt, poolBin, efficiencyWeight);
            }
            if (statusDsPrompt == 1 && statusPromptHadron == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionPromptDsPromptHadronMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
            } else if (statusDsPrompt == 0 && statusPromptHadron == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionNonPromptDsNonPromptHadronMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
            }
          }
        }
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeftMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeftMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeftMcRec"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRightMcRec"), deltaPhi, deltaEta, std::abs(ptD), std::abs(ptHadron), poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandRightMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRightMcRec"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcRecME, "Process MC Reco ME", false);

  /// Ds-Hadron correlation - for calculating candidate reconstruction efficiency using MC reco-level analysis
  void processMcCandEfficiency(soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels, aod::McCollisionLabels> const& collisions,
                               soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                               CandDsMcGen const& mcParticles,
                               CandDsMcReco const& candidates,
                               aod::TracksWMc const&)
  {
    /// loop over generated collisions
    for (const auto& mcCollision : mcCollisions) {

      const auto groupedCollisions = collisions.sliceBy(collPerCollMc, mcCollision.globalIndex());
      const auto groupedMcParticles = mcParticles.sliceBy(perCollisionCandMc, mcCollision.globalIndex());

      if (groupedCollisions.size() < 1) { // Skipping MC events that have no reconstructed collisions
        continue;
      }
      if (groupedCollisions.size() > 1 && removeCollWSplitVtx) { // Skipping MC events that have more than one reconstructed collision
        continue;
      }

      /// loop over reconstructed collisions
      for (const auto& collision : groupedCollisions) {

        // reco collision selection
        if (useSel8ForEff && !collision.sel8()) {
          continue;
        }
        if (std::abs(collision.posZ()) > cutCollPosZMc) {
          continue;
        }
        if (selNoSameBunchPileUpColl && !(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
          continue;
        }
        if (!collision.has_mcCollision()) {
          registry.fill(HIST("hFakeCollision"), 0.);
          continue;
        }

        registry.fill(HIST("hNumPvContrib"), collision.numContrib());

        const auto groupedCandidates = candidates.sliceBy(perCollisionCand, collision.globalIndex());

        // generated candidate loop
        for (const auto& mcParticle : groupedMcParticles) {
          if ((std::abs(mcParticle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) && (mcParticle.flagMcDecayChanGen() == channelsResonant[decayChannel])) {
            auto yDs = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassDS);
            if (std::abs(yDs) <= yCandGenMax) {
              if (mcParticle.originMcGen() == RecoDecay::OriginType::Prompt) {
                if (useHighDimHistoForEff) {
                  registry.fill(HIST("hPtCandMcGenPrompt"), mcParticle.pt(), collision.numContrib());
                } else {
                  registry.fill(HIST("hPtCandMcGenPrompt"), mcParticle.pt());
                }
              }
              if (mcParticle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
                if (useHighDimHistoForEff) {
                  registry.fill(HIST("hPtCandMcGenNonPrompt"), mcParticle.pt(), collision.numContrib());
                } else {
                  registry.fill(HIST("hPtCandMcGenNonPrompt"), mcParticle.pt());
                }
              }

              bool isDaughterInAcceptance = true;
              auto daughters = mcParticle.template daughters_as<CandDsMcGen>();
              for (const auto& daughter : daughters) {
                if (daughter.pt() < ptDaughterMin || std::abs(daughter.eta()) > etaTrackMax) {
                  isDaughterInAcceptance = false;
                }
              }
              if (isDaughterInAcceptance) {
                registry.fill(HIST("hPtCandMcGenDaughterInAcc"), mcParticle.pt());
              }
            }
          }
        } // end loop candidate gen

        // reconstructed candidate loop
        for (const auto& candidate : groupedCandidates) {
          if (candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
            continue;
          }
          std::vector<float> outputMl = {-1., -1., -1.};
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
            for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
              outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
            }
          } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
            for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
              outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
            }
          }
          if (outputMl[0] < mlOutputPromptMin->at(o2::analysis::findBin(binsPtD, candidate.pt())) || outputMl[0] > mlOutputPromptMax->at(o2::analysis::findBin(binsPtD, candidate.pt())) || outputMl[2] > mlOutputBkg->at(o2::analysis::findBin(binsPtD, candidate.pt()))) {
            continue;
          }

          if ((std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) && (candidate.flagMcDecayChanRec() == channelsResonant[decayChannel])) {
            auto prong0McPart = candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<CandDsMcGen>();
            // DsToKKPi and DsToPiKK division
            if (((std::abs(prong0McPart.pdgCode()) == kKPlus) && (candidate.isSelDsToKKPi() >= selectionFlagDs)) || ((std::abs(prong0McPart.pdgCode()) == kPiPlus) && (candidate.isSelDsToPiKK() >= selectionFlagDs))) {
              if (std::abs(HfHelper::yDs(candidate)) <= yCandMax) {
                if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
                  if (useHighDimHistoForEff) {
                    registry.fill(HIST("hPtCandMcRecPrompt"), candidate.pt(), collision.numContrib());
                  } else {
                    registry.fill(HIST("hPtCandMcRecPrompt"), candidate.pt());
                  }
                }
                if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
                  if (useHighDimHistoForEff) {
                    registry.fill(HIST("hPtCandMcRecNonPrompt"), candidate.pt(), collision.numContrib());
                  } else {
                    registry.fill(HIST("hPtCandMcRecNonPrompt"), candidate.pt());
                  }
                }
              }
            }
          }
        } // end loop candidate reco
      } // end loop reconstructed collision
    } // end loop generated collisions
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcCandEfficiency, "Process MC for calculating candidate reconstruction efficiency", false);

  void processMcCandEfficiencyWoColl(soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels> const&,
                                     soa::Join<aod::McCollisions, aod::MultsExtraMC> const&,
                                     CandDsMcGen const& mcParticles,
                                     CandDsMcReco const& candidates,
                                     aod::TracksWMc const&)
  {
    /// Gen loop
    for (const auto& mcParticle : mcParticles) {
      // generated candidates
      if ((std::abs(mcParticle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) && (mcParticle.flagMcDecayChanGen() == channelsResonant[decayChannel])) {
        auto yDs = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassDS);
        if (std::abs(yDs) <= yCandGenMax) {
          if (mcParticle.originMcGen() == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("hPtCandMcGenPrompt"), mcParticle.pt());
          }
          if (mcParticle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("hPtCandMcGenNonPrompt"), mcParticle.pt());
          }
        }
        bool isDaughterInAcceptance = true;
        auto daughters = mcParticle.template daughters_as<CandDsMcGen>();
        for (const auto& daughter : daughters) {
          if (daughter.pt() < ptDaughterMin || std::abs(daughter.eta()) > etaTrackMax) {
            isDaughterInAcceptance = false;
          }
        }
        if (isDaughterInAcceptance) {
          registry.fill(HIST("hPtCandMcGenDaughterInAcc"), mcParticle.pt());
        }
      }
    }

    // recontructed candidates loop
    for (const auto& candidate : candidates) {
      if (candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      std::vector<float> outputMl = {-1., -1., -1.};
      if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
        }
      } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
        }
      }
      if (outputMl[0] < mlOutputPromptMin->at(o2::analysis::findBin(binsPtD, candidate.pt())) || outputMl[0] > mlOutputPromptMax->at(o2::analysis::findBin(binsPtD, candidate.pt())) || outputMl[2] > mlOutputBkg->at(o2::analysis::findBin(binsPtD, candidate.pt()))) {
        continue;
      }
      auto collision = candidate.template collision_as<soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels>>();
      if (selNoSameBunchPileUpColl && !(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      if ((std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) && (candidate.flagMcDecayChanRec() == channelsResonant[decayChannel])) {
        auto prong0McPart = candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<CandDsMcGen>();
        // DsToKKPi and DsToPiKK division
        if (((std::abs(prong0McPart.pdgCode()) == kKPlus) && (candidate.isSelDsToKKPi() >= selectionFlagDs)) || ((std::abs(prong0McPart.pdgCode()) == kPiPlus) && (candidate.isSelDsToPiKK() >= selectionFlagDs))) {
          if (std::abs(HfHelper::yDs(candidate)) <= yCandMax) {
            if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hPtCandMcRecPrompt"), candidate.pt());
            }
            if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hPtCandMcRecNonPrompt"), candidate.pt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcCandEfficiencyWoColl, "Process MC for calculating candidate reconstruction efficiency", false);

  /// Ds-Hadron correlation - for calculating associated particle tracking efficiency using MC reco-level analysis
  void processMcTrackEfficiency(soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::McCollisionLabels> const& collisions,
                                soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                aod::McParticles const& mcParticles,
                                TracksWithMc const& tracksData)
  {

    /// loop over generated collisions
    for (const auto& mcCollision : mcCollisions) {

      const auto groupedCollisions = collisions.sliceBy(collPerCollMc, mcCollision.globalIndex());
      const auto groupedMcParticles = mcParticles.sliceBy(perCollisionMc, mcCollision.globalIndex());

      if (groupedCollisions.size() < 1) { // Skipping MC events that have no reconstructed collisions
        continue;
      }
      if (groupedCollisions.size() > 1 && removeCollWSplitVtx) { // Skipping MC events that have more than one reconstructed collision
        continue;
      }

      /// loop over reconstructed collisions
      for (const auto& collision : groupedCollisions) {

        // reco collision selection
        if (useSel8ForEff && !collision.sel8()) {
          continue;
        }
        if (std::abs(collision.posZ()) > cutCollPosZMc) {
          continue;
        }
        if (selNoSameBunchPileUpColl && !(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
          continue;
        }
        int centrality = getCentrality(collision);
        if (centrality < centralityMin || centrality > centralityMax) {
          continue;
        }
        if (!collision.has_mcCollision()) {
          registry.fill(HIST("hFakeCollision"), 0.);
          continue;
        }

        const auto groupedTracks = tracksData.sliceBy(perCollision, collision.globalIndex());

        // generated track loop
        for (const auto& mcParticle : groupedMcParticles) {
          if (mcParticle.isPhysicalPrimary() && ((std::abs(mcParticle.pdgCode()) == kElectron) || (std::abs(mcParticle.pdgCode()) == kMuonMinus) || (std::abs(mcParticle.pdgCode()) == kPiPlus) || (std::abs(mcParticle.pdgCode()) == kKPlus) || (std::abs(mcParticle.pdgCode()) == kProton))) {
            if (mcParticle.pt() > ptTrackMin && mcParticle.pt() < ptTrackMax) {
              if (std::abs(mcParticle.eta()) < etaTrackMax) {
                if (useHighDimHistoForEff) {
                  registry.fill(HIST("hPtParticleAssocMcGen"), mcParticle.pt(), mcParticle.eta(), collision.posZ(), collision.numContrib());
                  if (std::abs(mcParticle.pdgCode()) == kPiPlus) {
                    registry.fill(HIST("hPtPrmPionMcGen"), mcParticle.pt(), mcParticle.eta(), collision.posZ(), collision.numContrib());
                  } else if (std::abs(mcParticle.pdgCode()) == kKPlus) {
                    registry.fill(HIST("hPtPrmKaonMcGen"), mcParticle.pt(), mcParticle.eta(), collision.posZ(), collision.numContrib());
                  } else if (std::abs(mcParticle.pdgCode()) == kProton) {
                    registry.fill(HIST("hPtPrmProtonMcGen"), mcParticle.pt(), mcParticle.eta(), collision.posZ(), collision.numContrib());
                  } else if (std::abs(mcParticle.pdgCode()) == kElectron) {
                    registry.fill(HIST("hPtPrmElectronMcGen"), mcParticle.pt(), mcParticle.eta(), collision.posZ(), collision.numContrib());
                  } else if (std::abs(mcParticle.pdgCode()) == kMuonMinus) {
                    registry.fill(HIST("hPtPrmMuonMcGen"), mcParticle.pt(), mcParticle.eta(), collision.posZ(), collision.numContrib());
                  }
                } else {
                  registry.fill(HIST("hPtParticleAssocMcGen"), mcParticle.pt());
                  if (std::abs(mcParticle.pdgCode()) == kPiPlus) {
                    registry.fill(HIST("hPtPrmPionMcGen"), mcParticle.pt());
                  } else if (std::abs(mcParticle.pdgCode()) == kKPlus) {
                    registry.fill(HIST("hPtPrmKaonMcGen"), mcParticle.pt());
                  } else if (std::abs(mcParticle.pdgCode()) == kProton) {
                    registry.fill(HIST("hPtPrmProtonMcGen"), mcParticle.pt());
                  } else if (std::abs(mcParticle.pdgCode()) == kElectron) {
                    registry.fill(HIST("hPtPrmElectronMcGen"), mcParticle.pt());
                  } else if (std::abs(mcParticle.pdgCode()) == kMuonMinus) {
                    registry.fill(HIST("hPtPrmMuonMcGen"), mcParticle.pt());
                  }
                }
                if (separateTrackOrigins) {
                  int const trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
                  if (trackOrigin == RecoDecay::OriginType::Prompt) { // charm orgin
                    registry.fill(HIST("hPtPrmPromptPartMcGen"), mcParticle.pt());
                  } else if (trackOrigin == RecoDecay::OriginType::NonPrompt) { // beauty origin
                    registry.fill(HIST("hPtPrmNonPromptPartMcGen"), mcParticle.pt());
                  }
                }
              }
            }
          }
        }

        // reconstructed track loop
        for (const auto& track : groupedTracks) {
          if (!track.isGlobalTrackWoDCA() || track.tpcNClsCrossedRows() < nTpcCrossedRaws) {
            continue;
          }
          if (track.has_mcParticle()) {
            if (pidTrkApplied) {
              if (!passPIDSelection(track, trkPIDspecies, pidTPCMax, pidTOFMax, tofPIDThreshold, forceTOF)) {
                continue;
              }
            }
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            if (mcParticle.isPhysicalPrimary()) {
              registry.fill(HIST("hPtParticleAssocMcRec"), track.pt());
              if ((std::abs(mcParticle.pdgCode()) == kElectron) || (std::abs(mcParticle.pdgCode()) == kMuonMinus) || (std::abs(mcParticle.pdgCode()) == kPiPlus) || (std::abs(mcParticle.pdgCode()) == kKPlus) || (std::abs(mcParticle.pdgCode()) == kProton)) {
                // check the pt spectra of mcParticle
                registry.fill(HIST("hPtMcParticleAssocSpecieMcRec"), mcParticle.pt());
                if (useHighDimHistoForEff) {
                  registry.fill(HIST("hPtParticleAssocSpecieMcRec"), track.pt());
                  if (std::abs(mcParticle.pdgCode()) == kPiPlus) {
                    registry.fill(HIST("hPtPrmPionMcRec"), track.pt(), track.eta(), collision.posZ(), collision.numContrib());
                  } else if (std::abs(mcParticle.pdgCode()) == kKPlus) {
                    registry.fill(HIST("hPtPrmKaonMcRec"), track.pt(), track.eta(), collision.posZ(), collision.numContrib());
                  } else if (std::abs(mcParticle.pdgCode()) == kProton) {
                    registry.fill(HIST("hPtPrmProtonMcRec"), track.pt(), track.eta(), collision.posZ(), collision.numContrib());
                  } else if (std::abs(mcParticle.pdgCode()) == kElectron) {
                    registry.fill(HIST("hPtPrmElectronMcRec"), track.pt(), track.eta(), collision.posZ(), collision.numContrib());
                  } else if (std::abs(mcParticle.pdgCode()) == kMuonMinus) {
                    registry.fill(HIST("hPtPrmMuonMcRec"), track.pt(), track.eta(), collision.posZ(), collision.numContrib());
                  }
                } else {
                  registry.fill(HIST("hPtParticleAssocSpecieMcRec"), track.pt());
                  if (std::abs(mcParticle.pdgCode()) == kPiPlus) {
                    registry.fill(HIST("hPtPrmPionMcRec"), track.pt());
                  } else if (std::abs(mcParticle.pdgCode()) == kKPlus) {
                    registry.fill(HIST("hPtPrmKaonMcRec"), track.pt());
                  } else if (std::abs(mcParticle.pdgCode()) == kProton) {
                    registry.fill(HIST("hPtPrmProtonMcRec"), track.pt());
                  } else if (std::abs(mcParticle.pdgCode()) == kElectron) {
                    registry.fill(HIST("hPtPrmElectronMcRec"), track.pt());
                  } else if (std::abs(mcParticle.pdgCode()) == kMuonMinus) {
                    registry.fill(HIST("hPtPrmMuonMcRec"), track.pt());
                  }
                }
                // check track origin
                if (separateTrackOrigins) {
                  int const trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
                  if (trackOrigin == RecoDecay::OriginType::Prompt) { // charm orgin
                    registry.fill(HIST("hPtPrmPromptPartMcRec"), track.pt());
                  } else if (trackOrigin == RecoDecay::OriginType::NonPrompt) { // beauty origin
                    registry.fill(HIST("hPtPrmNonPromptPartMcRec"), track.pt());
                  }
                }
              }
            }
          } else {
            // fake track
            registry.fill(HIST("hFakeTracks"), 0.);
          }
        }

      } // end loop reconstructed collision
    } // end loop generated collisions
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcTrackEfficiency, "Process MC for calculating associated particle tracking efficiency", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDsHadrons>(cfgc)};
}
