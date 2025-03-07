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

#include <memory>
#include <string>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

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
  // Configurable<bool> doMcCollisionCheck{"doMcCollisionCheck", false, "Flag for applying the collision check and selection based on MC collision info"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds (avoid the case of flag = 0, no outputMlScore)"};
  Configurable<int> nTpcCrossedRaws{"nTpcCrossedRaws", 70, "Number of crossed TPC Rows"};
  // Configurable<int> eventGeneratorType{"eventGeneratorType", -1, "If positive, enable event selection using subGeneratorId information. The value indicates which events to keep (0 = MB, 4 = charm triggered, 5 = beauty triggered)"};
  Configurable<int> decayChannel{"decayChannel", 1, "Decay channels: 1 for Ds->PhiPi->KKpi, 2 for Ds->K0*K->KKPi"};
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
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle efficiency"};
  Configurable<std::vector<double>> mlOutputPrompt{"mlOutputPrompt", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for prompt"};
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

  std::shared_ptr<TH1> mEfficiencyD = nullptr;
  std::shared_ptr<TH1> mEfficiencyAssociated = nullptr;

  enum CandidateStep {
    kCandidateStepMcGenDsToKKPi = 0,
    kCandidateStepMcCandInAcceptance,
    kCandidateStepMcDaughtersInAcceptance,
    kCandidateStepMcReco,
    kCandidateStepMcRecoInAcceptance,
    kCandidateNSteps
  };

  enum AssocTrackStep { kAssocTrackStepMcGen = 0,
                        kAssocTrackStepMcGenInAcceptance,
                        kAssocTrackStepRecoAll,
                        kAssocTrackStepRecoMcMatch,
                        kAssocTrackStepRecoPrimaries,
                        kAssocTrackStepRecoSpecies,
                        kAssocTrackNSteps };

  HfHelper hfHelper;
  SliceCache cache;

  Service<ccdb::BasicCCDBManager> ccdb;

  using DsHadronPair = soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo>;
  using DsHadronPairFull = soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronGenInfo>;
  using DsHadronPairWithMl = soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronMlInfo, aod::TrackRecoInfo>;
  using DsHadronPairFullWithMl = soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronGenInfo, aod::DsHadronMlInfo, aod::TrackRecoInfo>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi, aod::HfCand3ProngMcRec>>; // flagDsFilter applied
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;                                                         // flagDsFilter applied
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels>>;   // trackFilter applied

  Filter flagDsFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

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
    AxisSpec axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec axisDsPrompt = {2, -0.5, 1.5, "Prompt Ds"};

    // Histograms for data analysis
    registry.add("hBdtScorePrompt", "Ds BDT prompt score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hBdtScoreBkg", "Ds BDT bkg score", {HistType::kTH1F, {axisBdtScore}});
    if (fillHistoData) {
      registry.add("hMassDsVsPt", "Ds candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
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

      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenPrompt"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenNonPrompt"))->Sumw2();
    }
    // Histograms for efficiencies
    if (fillHistoMcEff) {
      registry.add("hPtCandMcRecPrompt", "Ds prompt candidates pt", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcRecNonPrompt", "Ds non prompt candidates pt", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcGenPrompt", "Ds,Hadron particles prompt - MC Gen", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcGenNonPrompt", "Ds,Hadron particles non prompt - MC Gen", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtCandMcGenDaughterInAcc", "Ds,Hadron particles non prompt - MC Gen", {HistType::kTH1F, {axisPtD}});
      registry.add("hPtParticleAssocMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtParticleAssocSpecieMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtMcParticleAssocSpecieMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmPionMcRec", "Primary pions - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmKaonMcRec", "Primary kaons - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmProtonMcRec", "Primary protons - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmElectronMcRec", "Primary electrons - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmMuonMcRec", "Primary muons - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmPromptPartMcRec", "Primary prompt particles - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmNonPromptPartMcRec", "Primary non-prompt particles - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmPionMcGen", "Primary pions - MC Gen", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmKaonMcGen", "Primary kaons - MC Gen", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmProtonMcGen", "Primary protons - MC Gen", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmElectronMcGen", "Primary electrons - MC Gen", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmMuonMcGen", "Primary muons - MC Gen", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmPromptPartMcGen", "Primary prompt particles - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hPtPrmNonPromptPartMcGen", "Primary non-prompt particles - MC Rec", {HistType::kTH1F, {axisPtHadron}});
      registry.add("hFakeCollision", "Fake collision counter", {HistType::kTH1F, {{1, -0.5, 0.5, "n fake coll"}}});
      registry.add("hFakeTracks", "Fake tracks counter", {HistType::kTH1F, {{1, -0.5, 0.5, "n fake tracks"}}});

      auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPtD, axisMultFT0M, {RecoDecay::OriginType::NonPrompt + 1, +RecoDecay::OriginType::None - 0.5, +RecoDecay::OriginType::NonPrompt + 0.5}}, kCandidateNSteps});
      hCandidates->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
      hCandidates->GetAxis(1)->SetTitle("multiplicity");
      hCandidates->GetAxis(2)->SetTitle("Charm hadron origin");
      auto hAssocTracks = registry.add<StepTHn>("hAssocTracks", "Associated tracks at different steps", {HistType::kStepTHnF, {axisEta, axisPtHadron, axisMultFT0M, axisPosZ}, kAssocTrackNSteps});
      hAssocTracks->GetAxis(0)->SetTitle("#eta");
      hAssocTracks->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
      hAssocTracks->GetAxis(2)->SetTitle("multiplicity");
      hAssocTracks->GetAxis(3)->SetTitle("pos z");
    }

    // Loading efficiency histograms from CCDB
    if (applyEfficiency && loadAccXEffFromCCDB) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

      mEfficiencyD = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1F>(promptEffCcdbPath, timestampCcdb));
      if (mEfficiencyD == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", promptEffCcdbPath.value.c_str());
      }
      LOGF(info, "Loaded trigger efficiency (prompt D) histogram from %s", promptEffCcdbPath.value.c_str());

      mEfficiencyAssociated = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1F>(associatedEffCcdbPath, timestampCcdb));
      if (mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", associatedEffCcdbPath.value.c_str());
      }
      LOGF(info, "Loaded associated efficiency histogram from %s", associatedEffCcdbPath.value.c_str());
    }
  }

  void processData(DsHadronPairWithMl const& pairEntries,
                   aod::DsCandRecoInfo const& candidates)
  {
    for (const auto& candidate : candidates) {
      float massD = candidate.mD();
      float ptD = candidate.ptD();
      float bdtScorePrompt = candidate.mlScorePrompt();
      float bdtScoreBkg = candidate.mlScoreBkg();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      double efficiencyWeightD = 1.;
      if (applyEfficiency) {
        efficiencyWeightD = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD));
        if (loadAccXEffFromCCDB) {
          efficiencyWeightD = 1. / mEfficiencyD->GetBinContent(mEfficiencyD->FindBin(ptD));
        }
      }
      registry.fill(HIST("hMassDsVsPt"), massD, ptD, efficiencyWeightD);
      registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
      registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      float massD = pairEntry.mD();
      float bdtScorePrompt = pairEntry.mlScorePrompt();
      float bdtScoreBkg = pairEntry.mlScoreBkg();
      float trackDcaXY = pairEntry.trackDcaXY();
      float trackDcaZ = pairEntry.trackDcaZ();
      int trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int poolBin = pairEntry.poolBin();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
        if (loadAccXEffFromCCDB) {
          efficiencyWeight = 1. / (mEfficiencyD->GetBinContent(mEfficiencyD->FindBin(ptD)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
        }
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeft"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeft"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeft"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRight"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandRight"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRight"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processData, "Process data", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(DsHadronPairFullWithMl const& pairEntries,
                    soa::Join<aod::DsCandRecoInfo, aod::DsCandGenInfo> const& candidates)
  {
    for (const auto& candidate : candidates) {
      float massD = candidate.mD();
      float ptD = candidate.ptD();
      float bdtScorePrompt = candidate.mlScorePrompt();
      float bdtScoreBkg = candidate.mlScoreBkg();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);
      bool isDsPrompt = candidate.isPrompt();

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      double efficiencyWeightD = 1.;
      if (applyEfficiency) {
        efficiencyWeightD = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD));
        if (loadAccXEffFromCCDB) {
          efficiencyWeightD = 1. / mEfficiencyD->GetBinContent(mEfficiencyD->FindBin(ptD));
        }
      }
      if (isDsPrompt) {
        registry.fill(HIST("hMassPromptDsVsPt"), massD, ptD, efficiencyWeightD);
        registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
        registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
      } else {
        registry.fill(HIST("hMassNonPromptDsVsPt"), massD, ptD, efficiencyWeightD);
        registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
        registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
      }
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      float massD = pairEntry.mD();
      float bdtScorePrompt = pairEntry.mlScorePrompt();
      float bdtScoreBkg = pairEntry.mlScoreBkg();
      float trackDcaXY = pairEntry.trackDcaXY();
      float trackDcaZ = pairEntry.trackDcaZ();
      int trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int poolBin = pairEntry.poolBin();
      int statusDsPrompt = static_cast<int>(pairEntry.isPrompt());
      int statusPromptHadron = pairEntry.trackOrigin();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);
      bool isPhysicalPrimary = pairEntry.isPhysicalPrimary();

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
        if (loadAccXEffFromCCDB) {
          efficiencyWeight = 1. / (mEfficiencyD->GetBinContent(mEfficiencyD->FindBin(ptD)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
        }
      }
      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        // prompt and non-prompt division
        if (pairEntry.isSignal() && pairEntry.isDecayChan()) {
          registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
          registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, ptD, ptHadron, statusDsPrompt, poolBin, efficiencyWeight);
          if (isPhysicalPrimary) {
            registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"), deltaPhi, deltaEta, ptD, ptHadron, statusDsPrompt, poolBin, efficiencyWeight);
            if (statusDsPrompt == 1 && statusPromptHadron == 1) {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionPromptDsPromptHadronMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
            } else if (statusDsPrompt == 0 && statusPromptHadron == 2) {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionNonPromptDsNonPromptHadronMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
            }
          }
        }
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeftMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeftMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeftMcRec"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRightMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
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
      float deltaPhi = pairEntry.deltaPhi();
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      int statusPromptHadron = pairEntry.trackOrigin();
      bool isDsPrompt = pairEntry.isPrompt();

      registry.fill(HIST("hCorrel2DVsPtMcGen"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
      registry.fill(HIST("hDeltaEtaPtIntMcGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMcGen"), deltaPhi);
      if (isDsPrompt) {
        registry.fill(HIST("hCorrel2DVsPtMcGenPrompt"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        if (statusPromptHadron == 1) {
          registry.fill(HIST("hCorrel2DVsPtMcGenPromptDsPromptHadron"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        }
      } else {
        registry.fill(HIST("hCorrel2DVsPtMcGenNonPrompt"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        if (statusPromptHadron == 2) {
          registry.fill(HIST("hCorrel2DVsPtMcGenNonPromptDsNonPromptHadron"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
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
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      float massD = pairEntry.mD();
      float bdtScorePrompt = pairEntry.mlScorePrompt();
      float bdtScoreBkg = pairEntry.mlScoreBkg();
      float trackDcaXY = pairEntry.trackDcaXY();
      float trackDcaZ = pairEntry.trackDcaZ();
      int trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int poolBin = pairEntry.poolBin();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
        if (loadAccXEffFromCCDB) {
          efficiencyWeight = 1. / (mEfficiencyD->GetBinContent(mEfficiencyD->FindBin(ptD)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
        }
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeft"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeft"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeft"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRight"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandRight"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRight"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processDataME, "Process data ME", false);

  void processDerivedDataME(DsHadronPair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      float massD = pairEntry.mD();
      int poolBin = pairEntry.poolBin();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);

      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
        if (loadAccXEffFromCCDB) {
          efficiencyWeight = 1. / (mEfficiencyD->GetBinContent(mEfficiencyD->FindBin(ptD)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
        }
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeft"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeft"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeft"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRight"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandRight"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRight"), deltaPhi, efficiencyWeight);
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
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      float massD = pairEntry.mD();
      float bdtScorePrompt = pairEntry.mlScorePrompt();
      float bdtScoreBkg = pairEntry.mlScoreBkg();
      float trackDcaXY = pairEntry.trackDcaXY();
      float trackDcaZ = pairEntry.trackDcaZ();
      int trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int poolBin = pairEntry.poolBin();
      int statusDsPrompt = static_cast<int>(pairEntry.isPrompt());
      int statusPromptHadron = pairEntry.trackOrigin();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);
      bool isPhysicalPrimary = pairEntry.isPhysicalPrimary();

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
        if (loadAccXEffFromCCDB) {
          efficiencyWeight = 1. / (mEfficiencyD->GetBinContent(mEfficiencyD->FindBin(ptD)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
        }
      }
      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        // prompt and non-prompt division
        if (pairEntry.isSignal() && pairEntry.isDecayChan()) {
          registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
          registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
          registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, ptD, ptHadron, statusDsPrompt, poolBin, efficiencyWeight);
          if (isPhysicalPrimary) {
            registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"), deltaPhi, deltaEta, ptD, ptHadron, statusDsPrompt, poolBin, efficiencyWeight);
            if (statusDsPrompt == 1 && statusPromptHadron == 1) {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionPromptDsPromptHadronMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
            } else if (statusDsPrompt == 0 && statusPromptHadron == 2) {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionNonPromptDsNonPromptHadronMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
            }
          }
        }
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeftMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeftMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeftMcRec"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRightMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
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
    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));

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

        float multiplicityReco = collision.multFT0M();
        float multiplicityGen = mcCollision.multMCFT0A() + mcCollision.multMCFT0C(); // multFT0M = multFt0A + multFT0C

        const auto groupedCandidates = candidates.sliceBy(perCollisionCand, collision.globalIndex());

        // generated candidate loop
        for (const auto& mcParticle : groupedMcParticles) {
          if ((std::abs(mcParticle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) && (mcParticle.flagMcDecayChanGen() == decayChannel)) {
            hCandidates->Fill(kCandidateStepMcGenDsToKKPi, mcParticle.pt(), multiplicityGen, mcParticle.originMcGen());
            auto yDs = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassDS);
            if (std::abs(yDs) <= yCandGenMax) {
              hCandidates->Fill(kCandidateStepMcCandInAcceptance, mcParticle.pt(), multiplicityGen, mcParticle.originMcGen());
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
              hCandidates->Fill(kCandidateStepMcDaughtersInAcceptance, mcParticle.pt(), multiplicityGen, mcParticle.originMcGen());
              registry.fill(HIST("hPtCandMcGenDaughterInAcc"), mcParticle.pt());
            }
          }
        }

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
          if (outputMl[0] < mlOutputPrompt->at(o2::analysis::findBin(binsPtD, candidate.pt())) || outputMl[2] > mlOutputBkg->at(o2::analysis::findBin(binsPtD, candidate.pt()))) {
            continue;
          }

          if ((std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) && (candidate.flagMcDecayChanRec() == decayChannel)) {
            auto prong0McPart = candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<CandDsMcGen>();
            // DsToKKPi and DsToPiKK division
            if (((std::abs(prong0McPart.pdgCode()) == kKPlus) && (candidate.isSelDsToKKPi() >= selectionFlagDs)) || ((std::abs(prong0McPart.pdgCode()) == kPiPlus) && (candidate.isSelDsToPiKK() >= selectionFlagDs))) {
              hCandidates->Fill(kCandidateStepMcReco, candidate.pt(), multiplicityReco, candidate.originMcRec());
              if (std::abs(hfHelper.yDs(candidate)) <= yCandMax) {
                hCandidates->Fill(kCandidateStepMcRecoInAcceptance, candidate.pt(), multiplicityReco, candidate.originMcRec());
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
    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));

    /// Gen loop
    float multiplicity = -1.;
    for (const auto& mcParticle : mcParticles) {
      // generated candidates
      if ((std::abs(mcParticle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) && (mcParticle.flagMcDecayChanGen() == decayChannel)) {
        auto mcCollision = mcParticle.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();
        multiplicity = mcCollision.multMCFT0A() + mcCollision.multMCFT0C(); // multFT0M = multFt0A + multFT0C
        hCandidates->Fill(kCandidateStepMcGenDsToKKPi, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
        auto yDs = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassDS);
        if (std::abs(yDs) <= yCandGenMax) {
          hCandidates->Fill(kCandidateStepMcCandInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
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
          hCandidates->Fill(kCandidateStepMcDaughtersInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
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
      if (outputMl[0] < mlOutputPrompt->at(o2::analysis::findBin(binsPtD, candidate.pt())) || outputMl[2] > mlOutputBkg->at(o2::analysis::findBin(binsPtD, candidate.pt()))) {
        continue;
      }
      auto collision = candidate.template collision_as<soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels>>();
      if (selNoSameBunchPileUpColl && !(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      multiplicity = collision.multFT0M();
      if ((std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) && (candidate.flagMcDecayChanRec() == decayChannel)) {
        auto prong0McPart = candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<CandDsMcGen>();
        // DsToKKPi and DsToPiKK division
        if (((std::abs(prong0McPart.pdgCode()) == kKPlus) && (candidate.isSelDsToKKPi() >= selectionFlagDs)) || ((std::abs(prong0McPart.pdgCode()) == kPiPlus) && (candidate.isSelDsToPiKK() >= selectionFlagDs))) {
          hCandidates->Fill(kCandidateStepMcReco, candidate.pt(), multiplicity, candidate.originMcRec());
          if (std::abs(hfHelper.yDs(candidate)) <= yCandMax) {
            hCandidates->Fill(kCandidateStepMcRecoInAcceptance, candidate.pt(), multiplicity, candidate.originMcRec());
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
  void processMcTrackEfficiency(soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels, aod::McCollisionLabels> const& collisions,
                                soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                aod::McParticles const& mcParticles,
                                TracksWithMc const& tracksData)
  {
    auto hAssocTracks = registry.get<StepTHn>(HIST("hAssocTracks"));

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
        if (!collision.has_mcCollision()) {
          registry.fill(HIST("hFakeCollision"), 0.);
          continue;
        }

        float multiplicityReco = collision.multFT0M();
        float posZReco = collision.posZ();
        float multiplicityGen = mcCollision.multMCFT0A() + mcCollision.multMCFT0C(); // multFT0M = multFt0A + multFT0C
        float posZGen = mcCollision.posZ();

        const auto groupedTracks = tracksData.sliceBy(perCollision, collision.globalIndex());

        // generated track loop
        for (const auto& mcParticle : groupedMcParticles) {
          if (mcParticle.isPhysicalPrimary() && ((std::abs(mcParticle.pdgCode()) == kElectron) || (std::abs(mcParticle.pdgCode()) == kMuonMinus) || (std::abs(mcParticle.pdgCode()) == kPiPlus) || (std::abs(mcParticle.pdgCode()) == kKPlus) || (std::abs(mcParticle.pdgCode()) == kProton))) {
            if (mcParticle.pt() > ptTrackMin && mcParticle.pt() < ptTrackMax) {
              hAssocTracks->Fill(kAssocTrackStepMcGen, mcParticle.eta(), mcParticle.pt(), multiplicityGen, posZGen);
              if (std::abs(mcParticle.eta()) < etaTrackMax) {
                hAssocTracks->Fill(kAssocTrackStepMcGenInAcceptance, mcParticle.eta(), mcParticle.pt(), multiplicityGen, posZGen);
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
                if (separateTrackOrigins) {
                  int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
                  if (trackOrigin == 1) { // charm orgin
                    registry.fill(HIST("hPtPrmPromptPartMcGen"), mcParticle.pt());
                  } else if (trackOrigin == 2) { // beauty origin
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
            hAssocTracks->Fill(kAssocTrackStepRecoMcMatch, track.eta(), track.pt(), multiplicityReco, posZReco);
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            if (mcParticle.isPhysicalPrimary()) {
              hAssocTracks->Fill(kAssocTrackStepRecoPrimaries, track.eta(), track.pt(), multiplicityReco, posZReco);
              registry.fill(HIST("hPtParticleAssocMcRec"), track.pt());
              if ((std::abs(mcParticle.pdgCode()) == kElectron) || (std::abs(mcParticle.pdgCode()) == kMuonMinus) || (std::abs(mcParticle.pdgCode()) == kPiPlus) || (std::abs(mcParticle.pdgCode()) == kKPlus) || (std::abs(mcParticle.pdgCode()) == kProton)) {
                hAssocTracks->Fill(kAssocTrackStepRecoSpecies, track.eta(), track.pt(), multiplicityReco, posZReco);
                registry.fill(HIST("hPtParticleAssocSpecieMcRec"), track.pt());
                // check the pt spectra of mcParticle
                registry.fill(HIST("hPtMcParticleAssocSpecieMcRec"), mcParticle.pt());
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
                // check track origin
                if (separateTrackOrigins) {
                  int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
                  if (trackOrigin == 1) { // charm orgin
                    registry.fill(HIST("hPtPrmPromptPartMcRec"), track.pt());
                  } else if (trackOrigin == 2) { // beauty origin
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
