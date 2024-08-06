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

/// \file taskDs.cxx
/// \brief Ds± analysis task
/// \note Extended from taskD0 and taskDplus
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Universita and INFN Torino
/// \author Stefano Politanò <stefano.politano@cern.ch>, Politecnico & INFN Torino
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Universita and INFN Torino

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum FinalState { KKPi = 0,
                  PiKK };

enum DataType { Data = 0,
                McDsPrompt,
                McDsNonPrompt,
                McDplusPrompt,
                McDplusNonPrompt,
                McDplusBkg,
                McBkg,
                kDataTypes };

enum SpeciesAndDecay { DsToKKPi = 0,
                       DplusToKKPi,
                       DplusToPiKPi,
                       kSpeciesAndDecay };

/// Ds± analysis task
struct HfTaskDs {

  Configurable<int> decayChannel{"decayChannel", 1, "Switch between decay channels: 1 for Ds/Dplus->PhiPi->KKpi, 2 for Ds/Dplus->K0*K->KKPi"};
  Configurable<bool> fillDplusMc{"fillDplusMc", true, "Switch to fill Dplus MC information"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2, 3}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<double> minMassDsSignal{"minMassDsSignal", 1.934, "min mass for Ds signal"};
  Configurable<double> maxMassDsSignal{"maxMassDsSignal", 1.994, "max mass for Ds signal"};
  Configurable<double> minMassDplusSignal{"minMassDplusSignal", 1.866, "min mass for Dplus signal"};
  Configurable<double> maxMassDplusSignal{"maxMassDplusSignal", 1.906, "max mass for Dplus signal"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f}, "axis for pT"};
  ConfigurableAxis axisNPvContributors{"axisNPvContributors", {200, -0.5f, 199.5f}, "axis for NPvContributors"};
  ConfigurableAxis axisMlScore0{"axisMlScore0", {100, 0., 1.}, "axis for ML output score 0"};
  ConfigurableAxis axisMlScore1{"axisMlScore1", {100, 0., 1.}, "axis for ML output score 1"};
  ConfigurableAxis axisMlScore2{"axisMlScore2", {100, 0., 1.}, "axis for ML output score 2"};

  HfHelper hfHelper;

  using CentralityEstimator = o2::hf_centrality::CentralityEstimator;
  using TH1_ptr = std::shared_ptr<TH1>;
  using TH2_ptr = std::shared_ptr<TH2>;
  using THnSparse_ptr = std::shared_ptr<THnSparse>;
  using histTypes = std::variant<TH1_ptr, TH2_ptr, THnSparse_ptr>;

  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsWithNTracksPV = soa::Join<aod::Collisions, aod::EvSels, aod::CentNTPVs>;

  using CollisionsMc = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using CollisionsMcWithFT0C = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsMcWithFT0M = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsMcWithNTracksPV = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentNTPVs>;

  PresliceUnsorted<CollisionsMc> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsMcWithFT0C> colPerMcCollisionWithFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsMcWithFT0M> colPerMcCollisionWithFT0M = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsMcWithNTracksPV> colPerMcCollisionWithNTracksPV = aod::mccollisionlabel::mcCollisionId;
  SliceCache cache;

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsDataWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcRecoWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec, aod::HfMlDsToKKPi>>;
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Preslice<CandDsData> candDsDataPerCollision = aod::hf_cand::collisionId;
  Preslice<CandDsDataWithMl> candDsDataWithMlPerCollision = aod::hf_cand::collisionId;
  Preslice<CandDsMcReco> candDsMcRecoPerCollision = aod::hf_cand::collisionId;
  Preslice<CandDsMcRecoWithMl> candDsMcRecoWithMlPerCollision = aod::hf_cand::collisionId;

  int offsetDplusDecayChannel = aod::hf_cand_3prong::DecayChannelDToKKPi::DplusToPhiPi - aod::hf_cand_3prong::DecayChannelDToKKPi::DsToPhiPi; // Offset between Dplus and Ds to use the same decay channel. See aod::hf_cand_3prong::DecayChannelDToKKPi

  Filter filterDsFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi))) != static_cast<uint8_t>(0);

  // Data
  Partition<CandDsData> selectedDsToKKPiCandData = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs;
  Partition<CandDsDataWithMl> selectedDsToKKPiCandWithMlData = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs;
  Partition<CandDsData> selectedDsToPiKKCandData = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Partition<CandDsDataWithMl> selectedDsToPiKKCandWithMlData = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  // MC
  Partition<CandDsMcReco> selectedDsToKKPiCandMc = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs;
  Partition<CandDsMcRecoWithMl> selectedDsToKKPiCandWithMlMc = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs;
  Partition<CandDsMcReco> selectedDsToPiKKCandMc = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Partition<CandDsMcRecoWithMl> selectedDsToPiKKCandWithMlMc = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  // Matched MC, no ML
  Partition<CandDsMcReco> reconstructedCandDsSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && aod::hf_cand_3prong::flagMcDecayChanRec == decayChannel;
  Partition<CandDsMcReco> reconstructedCandDplusSig = fillDplusMc && nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && aod::hf_cand_3prong::flagMcDecayChanRec == (decayChannel + offsetDplusDecayChannel);
  Partition<CandDsMcReco> reconstructedCandDplusBkg = fillDplusMc && nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi));
  Partition<CandDsMcReco> reconstructedCandBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi));

  // Matched MC, with ML
  Partition<CandDsMcRecoWithMl> reconstructedCandDsSigWithMl = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && aod::hf_cand_3prong::flagMcDecayChanRec == decayChannel;
  Partition<CandDsMcRecoWithMl> reconstructedCandDplusSigWithMl = fillDplusMc && nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && aod::hf_cand_3prong::flagMcDecayChanRec == (decayChannel + offsetDplusDecayChannel);
  Partition<CandDsMcRecoWithMl> reconstructedCandDplusBkgWithMl = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi)) && aod::hf_cand_3prong::flagMcDecayChanRec == decayChannel;
  Partition<CandDsMcRecoWithMl> reconstructedCandBkgWithMl = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi));

  HistogramRegistry registry{"registry", {}};

  std::array<std::string, DataType::kDataTypes> folders = {"Data/", "MC/Ds/Prompt/", "MC/Ds/NonPrompt/", "MC/Dplus/Prompt/", "MC/Dplus/NonPrompt/", "MC/Dplus/Bkg/", "MC/Bkg/"};

  std::unordered_map<std::string, histTypes> dataHistograms = {};
  std::unordered_map<std::string, histTypes> mcDsPromptHistograms = {};
  std::unordered_map<std::string, histTypes> mcDsNonPromptHistograms = {};
  std::unordered_map<std::string, histTypes> mcDplusPromptHistograms = {};
  std::unordered_map<std::string, histTypes> mcDplusNonPromptHistograms = {};
  std::unordered_map<std::string, histTypes> mcDplusBkgHistograms = {};
  std::unordered_map<std::string, histTypes> mcBkgHistograms = {};

  std::array<std::unordered_map<std::string, histTypes>, DataType::kDataTypes> histosPtr = {dataHistograms, mcDsPromptHistograms, mcDsNonPromptHistograms, mcDplusPromptHistograms, mcDplusNonPromptHistograms, mcDplusBkgHistograms, mcBkgHistograms};

  void init(InitContext&)
  {
    std::array<int, 16> processes = {doprocessDataWithCentFT0C, doprocessDataWithCentFT0M, doprocessDataWithCentNTracksPV, doprocessData, doprocessDataWithMlAndCentFT0C, doprocessDataWithMlAndCentFT0M, doprocessDataWithMlAndCentNTracksPV, doprocessDataWithMl, doprocessMcWithCentFT0C, doprocessMcWithCentFT0M, doprocessMcWithCentNTracksPV, doprocessMc, doprocessMcWithMlAndCentFT0C, doprocessMcWithMlAndCentFT0M, doprocessMcWithMlAndCentNTracksPV, doprocessMcWithMl};

    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    AxisSpec ptbins{axisPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ybins = {100, -5., 5, "#it{y}"};
    AxisSpec massbins = {600, 1.67, 2.27, "inv. mass (KK#pi) (GeV/#it{c}^{2})"};
    AxisSpec centralitybins = {100, 0., 100., "Centrality"};

    histosPtr[DataType::Data]["hNPvContribAll"] = registry.add<TH2>((folders[DataType::Data] + "hNPvContribAll").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, {100, 0., 100}});

    for (auto i = 0; i < DataType::kDataTypes; ++i) {
      if (doprocessDataWithCentFT0C || doprocessDataWithCentFT0M || doprocessDataWithCentNTracksPV ||
          doprocessMcWithCentFT0C || doprocessMcWithCentFT0M || doprocessMcWithCentNTracksPV) {
        if (i == DataType::Data) // If data do not fill PV contributors in sparse
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins});
        else
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins, axisNPvContributors});

        histosPtr[i]["hNPvContribCands"] = registry.add<TH2>((folders[i] + "hNPvContribCands").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
        histosPtr[i]["hNPvContribCandsInSignalRegionDs"] = registry.add<TH2>((folders[i] + "hNPvContribCandsInSignalRegionDs").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
        histosPtr[i]["hNPvContribCandsInSignalRegionDplus"] = registry.add<TH2>((folders[i] + "hNPvContribCandsInSignalRegionDplus").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
      } else if (doprocessDataWithMlAndCentFT0C || doprocessDataWithMlAndCentFT0M || doprocessDataWithMlAndCentNTracksPV ||
                 doprocessMcWithMlAndCentFT0C || doprocessMcWithMlAndCentFT0M || doprocessMcWithMlAndCentNTracksPV) {
        if (i == DataType::Data) // If data do not fill PV contributors in sparse
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins, axisMlScore0, axisMlScore1, axisMlScore2});
        else
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins, axisMlScore0, axisMlScore1, axisMlScore2, axisNPvContributors});
        histosPtr[i]["hNPvContribCands"] = registry.add<TH2>((folders[i] + "hNPvContribCands").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
        histosPtr[i]["hNPvContribCandsInSignalRegionDs"] = registry.add<TH2>((folders[i] + "hNPvContribCandsInSignalRegionDs").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
        histosPtr[i]["hNPvContribCandsInSignalRegionDplus"] = registry.add<TH2>((folders[i] + "hNPvContribCandsInSignalRegionDplus").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
      } else if (doprocessData || doprocessMc) {
        if (i == DataType::Data) // If data do not fill PV contributors in sparse
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins});
        else
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, axisNPvContributors});
      } else if (doprocessDataWithMl || doprocessMcWithMl) {
        if (i == DataType::Data) // If data do not fill PV contributors in sparse
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, axisMlScore0, axisMlScore1, axisMlScore2});
        else
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, axisMlScore0, axisMlScore1, axisMlScore2, axisNPvContributors});
      }
      histosPtr[i]["hPt"] = registry.add<TH1>((folders[i] + "hPt").c_str(), "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      histosPtr[i]["hPtProng0"] = registry.add<TH1>((folders[i] + "hPtProng0").c_str(), "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      histosPtr[i]["hPtProng1"] = registry.add<TH1>((folders[i] + "hPtProng1").c_str(), "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      histosPtr[i]["hPtProng2"] = registry.add<TH1>((folders[i] + "hPtProng2").c_str(), "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      histosPtr[i]["hEta"] = registry.add<TH2>((folders[i] + "hEta").c_str(), "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, ptbins}});
      histosPtr[i]["hCt"] = registry.add<TH2>((folders[i] + "hCt").c_str(), "3-prong candidates;proper lifetime (D_{s}^{#pm}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 100}, ptbins}});
      histosPtr[i]["hDecayLength"] = registry.add<TH2>((folders[i] + "hDecayLength").c_str(), "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, ptbins}});
      histosPtr[i]["hDecayLengthXY"] = registry.add<TH2>((folders[i] + "hDecayLengthXY").c_str(), "3-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, ptbins}});
      histosPtr[i]["hNormalisedDecayLengthXY"] = registry.add<TH2>((folders[i] + "hNormalisedDecayLengthXY").c_str(), "3-prong candidates;norm. decay length xy;entries", {HistType::kTH2F, {{80, 0., 80.}, ptbins}});
      histosPtr[i]["hCPA"] = registry.add<TH2>((folders[i] + "hCPA").c_str(), "3-prong candidates;cos. pointing angle;entries", {HistType::kTH2F, {{100, -1., 1.}, ptbins}});
      histosPtr[i]["hCPAxy"] = registry.add<TH2>((folders[i] + "hCPAxy").c_str(), "3-prong candidates;cos. pointing angle xy;entries", {HistType::kTH2F, {{100, -1., 1.}, ptbins}});
      histosPtr[i]["hImpactParameterXY"] = registry.add<TH2>((folders[i] + "hImpactParameterXY").c_str(), "3-prong candidates;impact parameter xy (cm);entries", {HistType::kTH2F, {{200, -1., 1.}, ptbins}});
      histosPtr[i]["hMaxNormalisedDeltaIP"] = registry.add<TH2>((folders[i] + "hMaxNormalisedDeltaIP").c_str(), "3-prong candidates;norm. IP;entries", {HistType::kTH2F, {{200, -20., 20.}, ptbins}});
      histosPtr[i]["hCos3PiK"] = registry.add<TH2>((folders[i] + "hCos3PiK").c_str(), "3-prong candidates;cos^{3} #theta'(K);entries", {HistType::kTH2F, {{100, -1., 1.}, ptbins}});
      histosPtr[i]["hAbsCos3PiK"] = registry.add<TH2>((folders[i] + "hAbsCos3PiK").c_str(), "3-prong candidates;|cos^{3} #theta'(K)|;entries", {HistType::kTH2F, {{100, 0., 1.}, ptbins}});
      histosPtr[i]["hDeltaMassPhi"] = registry.add<TH2>((folders[i] + "hDeltaMassPhi").c_str(), "3-prong candidates;|M(KK) - M(#phi)| (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{100, 0., 0.1}, ptbins}});
      histosPtr[i]["hMassKK"] = registry.add<TH2>((folders[i] + "hMassKK").c_str(), "3-prong candidates;M(KK) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{100, o2::constants::physics::MassPhi - 0.05, o2::constants::physics::MassPhi + 0.05}, ptbins}});
      histosPtr[i]["hImpactParameterProngSqSum"] = registry.add<TH2>((folders[i] + "hImpactParameterProngSqSum").c_str(), "3-prong candidates;squared sum of prong imp. par. (cm^{2});entries", {HistType::kTH2F, {{100, 0., 1.}, ptbins}});
      histosPtr[i]["hDecayLengthError"] = registry.add<TH2>((folders[i] + "hDecayLengthError").c_str(), "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, ptbins}});
      histosPtr[i]["hDecayLengthXYError"] = registry.add<TH2>((folders[i] + "hDecayLengthXYError").c_str(), "3-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, ptbins}});
      histosPtr[i]["hImpactParameterError"] = registry.add<TH2>((folders[i] + "hImpactParameterError").c_str(), "3-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, ptbins}});
      histosPtr[i]["hd0Prong0"] = registry.add<TH2>((folders[i] + "hd0Prong0").c_str(), "3-prong candidates;prong 0 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, ptbins}});
      histosPtr[i]["hd0Prong1"] = registry.add<TH2>((folders[i] + "hd0Prong1").c_str(), "3-prong candidates;prong 1 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, ptbins}});
      histosPtr[i]["hd0Prong2"] = registry.add<TH2>((folders[i] + "hd0Prong2").c_str(), "3-prong candidates;prong 2 DCA to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, ptbins}});
    }

    if (doprocessMcWithCentFT0C || doprocessMcWithCentFT0M || doprocessMcWithCentNTracksPV ||
        doprocessMcWithMlAndCentFT0C || doprocessMcWithMlAndCentFT0M || doprocessMcWithMlAndCentNTracksPV ||
        doprocessMc || doprocessMcWithMl) { // processing MC

      for (auto i = 0; i < DataType::kDataTypes; ++i) {
        if (i == DataType::McDsPrompt || i == DataType::McDsNonPrompt || i == DataType::McDplusPrompt || i == DataType::McDplusNonPrompt || i == DataType::McDplusBkg) {

          histosPtr[i]["hEtaGen"] = registry.add<TH1>((folders[i] + "hEtaGen").c_str(), "3-prong candidates (matched);#eta;entries", {HistType::kTH1F, {{100, -2., 2.}}});
          histosPtr[i]["hEtaRecSig"] = registry.add<TH1>((folders[i] + "hEtaRecSig").c_str(), "3-prong candidates (matched);#eta;entries", {HistType::kTH1F, {{100, -2., 2.}}});
          histosPtr[i]["hCPARecSig"] = registry.add<TH1>((folders[i] + "hCPARecSig").c_str(), "3-prong candidates (matched);cos. pointing angle;entries", {HistType::kTH1F, {{100, -1., 1.}}});
          histosPtr[i]["hPtRecSig"] = registry.add<TH1>((folders[i] + "hPtRecSig").c_str(), "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {ptbins}});
          histosPtr[i]["hPtGenSig"] = registry.add<TH1>((folders[i] + "hPtGenSig").c_str(), "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {ptbins}});
          histosPtr[i]["hPtGen"] = registry.add<TH1>((folders[i] + "hPtGen").c_str(), "MC particles (unmatched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {ptbins}});
          histosPtr[i]["hPtVsYRecSigRecoPID"] = registry.add<TH2>((folders[i] + "hPtVsYRecSigRecoPID").c_str(), "3-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
          histosPtr[i]["hPtVsYRecSigRecoTopol"] = registry.add<TH2>((folders[i] + "hPtVsYRecSigRecoTopol").c_str(), "3-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
          histosPtr[i]["hPtVsYRecSigRecoSkim"] = registry.add<TH2>((folders[i] + "hPtVsYRecSigRecoSkim").c_str(), "3-prong candidates (RecoSkim - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
          histosPtr[i]["hPtYNPvContribGen"] = registry.add<THnSparse>((folders[i] + "hPtYNPvContribGen").c_str(), "Thn for generated candidates", {HistType::kTHnSparseF, {ptbins, {ybins}, axisNPvContributors}});
        }
      }
    }
  }

  /// Evaluate centrality/multiplicity percentile
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <CentralityEstimator centDetector, typename Coll>
  int evaluateCentralityColl(const Coll& collision)
  {
    if constexpr (centDetector == CentralityEstimator::FT0C)
      return collision.centFT0C();
    else if constexpr (centDetector == CentralityEstimator::FT0M)
      return collision.centFT0M();
    else if constexpr (centDetector == CentralityEstimator::NTracksPV)
      return collision.centNTPV();
  }

  /// Evaluate centrality/multiplicity percentile
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision associated to the candidate
  template <CentralityEstimator centDetector, typename Coll, typename T1>
  int evaluateCentralityCand(const T1& candidate)
  {
    if constexpr (centDetector == CentralityEstimator::FT0C)
      return evaluateCentralityColl<CentralityEstimator::FT0C>(candidate.template collision_as<Coll>());
    else if constexpr (centDetector == CentralityEstimator::FT0M)
      return evaluateCentralityColl<CentralityEstimator::FT0M>(candidate.template collision_as<Coll>());
    else if constexpr (centDetector == CentralityEstimator::NTracksPV)
      return evaluateCentralityColl<CentralityEstimator::NTracksPV>(candidate.template collision_as<Coll>());
  }

  /// Fill histograms of quantities independent from the daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <CentralityEstimator centDetector, typename Coll, typename T1>
  void fillHisto(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    std::get<TH1_ptr>(histosPtr[dataType]["hPt"])->Fill(pt);
    std::get<TH1_ptr>(histosPtr[dataType]["hPtProng0"])->Fill(candidate.ptProng0());
    std::get<TH1_ptr>(histosPtr[dataType]["hPtProng1"])->Fill(candidate.ptProng1());
    std::get<TH1_ptr>(histosPtr[dataType]["hPtProng2"])->Fill(candidate.ptProng2());
    std::get<TH2_ptr>(histosPtr[dataType]["hEta"])->Fill(candidate.eta(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hCt"])->Fill(hfHelper.ctDs(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hDecayLength"])->Fill(candidate.decayLength(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hDecayLengthXY"])->Fill(candidate.decayLengthXY(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hNormalisedDecayLengthXY"])->Fill(candidate.decayLengthXYNormalised(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hCPA"])->Fill(candidate.cpa(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hCPAxy"])->Fill(candidate.cpaXY(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hImpactParameterXY"])->Fill(candidate.impactParameterXY(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hMaxNormalisedDeltaIP"])->Fill(candidate.maxNormalisedDeltaIP(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hImpactParameterProngSqSum"])->Fill(candidate.impactParameterProngSqSum(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hDecayLengthError"])->Fill(candidate.errorDecayLength(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hDecayLengthXYError"])->Fill(candidate.errorDecayLengthXY(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hImpactParameterError"])->Fill(candidate.errorImpactParameter0(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hImpactParameterError"])->Fill(candidate.errorImpactParameter1(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hImpactParameterError"])->Fill(candidate.errorImpactParameter2(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hd0Prong0"])->Fill(candidate.impactParameter0(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hd0Prong1"])->Fill(candidate.impactParameter1(), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hd0Prong2"])->Fill(candidate.impactParameter2(), pt);

    return;
  }

  /// Fill histograms of quantities for the KKPi daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <CentralityEstimator centDetector, bool useMl, typename Coll, typename T1>
  void fillHistoKKPi(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    if (dataType == DataType::Data) { // If data do not fill PV contributors in sparse
      if constexpr (useMl) {
        std::vector<float> outputMl = {-999., -999., -999.};
        for (unsigned int iclass = 0; iclass < classMl->size() && candidate.mlProbDsToKKPi().size() != 0; iclass++) { // TODO: add checks for classMl size
          outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
        }
        if constexpr (centDetector != CentralityEstimator::None) {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, evaluateCentralityCand<centDetector, Coll>(candidate), outputMl[0], outputMl[1], outputMl[2]);
        } else {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, outputMl[0], outputMl[1], outputMl[2]);
        }
      } else {
        if constexpr (centDetector != CentralityEstimator::None) {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, evaluateCentralityCand<centDetector, Coll>(candidate));
        } else {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt);
        }
      }
    } else {
      if constexpr (useMl) {
        std::vector<float> outputMl = {-999., -999., -999.};
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) { // TODO: add checks for classMl size
        if (candidate.mlProbDsToKKPi().size() == 0)
            continue;
          outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
        }
        if constexpr (centDetector != CentralityEstimator::None) {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, evaluateCentralityCand<centDetector, Coll>(candidate), outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib());
        } else {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib());
        }
      } else {
        if constexpr (centDetector != CentralityEstimator::None) {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, evaluateCentralityCand<centDetector, Coll>(candidate), candidate.template collision_as<Coll>().numContrib());
        } else {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, candidate.template collision_as<Coll>().numContrib());
        }
      }
    }

    std::get<TH2_ptr>(histosPtr[dataType]["hCos3PiK"])->Fill(hfHelper.cos3PiKDsToKKPi(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hAbsCos3PiK"])->Fill(hfHelper.absCos3PiKDsToKKPi(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hDeltaMassPhi"])->Fill(hfHelper.deltaMassPhiDsToKKPi(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hMassKK"])->Fill(hfHelper.massKKPairDsToKKPi(candidate), pt);

    return;
  }

  /// Fill histograms of quantities for the PiKK daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <CentralityEstimator centDetector, bool useMl, typename Coll, typename T1>
  void fillHistoPiKK(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();

    if (dataType == DataType::Data) { // If data do not fill PV contributors in sparse
      if constexpr (useMl) {
        std::vector<float> outputMl = {-999., -999., -999.};
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) { // TODO: add checks for classMl size
          if (candidate.mlProbDsToPiKK().size() == 0)
            continue;
          outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
        }
        if constexpr (centDetector != CentralityEstimator::None) {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, evaluateCentralityCand<centDetector, Coll>(candidate), outputMl[0], outputMl[1], outputMl[2]);
        } else {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, outputMl[0], outputMl[1], outputMl[2]);
        }
      } else {
        if constexpr (centDetector != CentralityEstimator::None) {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, evaluateCentralityCand<centDetector, Coll>(candidate));
        } else {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt);
        }
      }
    } else {
      if constexpr (useMl) {
        std::vector<float> outputMl = {-999., -999., -999.};
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) { // TODO: add checks for classMl size
          outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
        }
        if constexpr (centDetector != CentralityEstimator::None) {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, evaluateCentralityCand<centDetector, Coll>(candidate), outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib());
        } else {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib());
        }
      } else {
        if constexpr (centDetector != CentralityEstimator::None) {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, evaluateCentralityCand<centDetector, Coll>(candidate), candidate.template collision_as<Coll>().numContrib());
        } else {
          std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, candidate.template collision_as<Coll>().numContrib());
        }
      }
    }

    std::get<TH2_ptr>(histosPtr[dataType]["hCos3PiK"])->Fill(hfHelper.cos3PiKDsToPiKK(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hAbsCos3PiK"])->Fill(hfHelper.absCos3PiKDsToPiKK(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hDeltaMassPhi"])->Fill(hfHelper.deltaMassPhiDsToPiKK(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hMassKK"])->Fill(hfHelper.massKKPairDsToPiKK(candidate), pt);

    return;
  }

  /// Fill MC histograms at reconstruction level
  /// \param candidate is candidate
  /// \param mcParticles are particles with MC information
  /// \param whichSpeciesDecay defines which histogram to fill
  template <CentralityEstimator centDetector, bool useMl, typename Coll, typename T1>
  void fillHistoMCRec(const T1& candidate, const CandDsMcGen& mcParticles, SpeciesAndDecay whichSpeciesDecay)
  {
    auto indexMother = RecoDecay::getMother(mcParticles,
                                            candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<CandDsMcGen>(),
                                            whichSpeciesDecay == SpeciesAndDecay::DsToKKPi ? o2::constants::physics::Pdg::kDS : o2::constants::physics::Pdg::kDPlus, true);

    if (indexMother != -1) {
      if (yCandRecoMax >= 0. && std::abs(whichSpeciesDecay == SpeciesAndDecay::DsToKKPi ? hfHelper.yDs(candidate) : hfHelper.yDplus(candidate)) > yCandRecoMax) {
        return;
      }

      auto particleMother = mcParticles.iteratorAt(indexMother);

      int flag = candidate.isCandidateSwapped() ? candidate.isSelDsToPiKK() : candidate.isSelDsToKKPi(); // 0 corresponds to KKPi, 1 to PiKK

      auto pt = candidate.pt(); // rec. level pT

      // Ds
      if (whichSpeciesDecay == SpeciesAndDecay::DsToKKPi) {

        double y = hfHelper.yDs(candidate);

        // prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
            fillHisto<centDetector, Coll>(candidate, DataType::McDsPrompt);
            fillHistoKKPi<centDetector, useMl, Coll>(candidate, DataType::McDsPrompt);
          }
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
            fillHisto<centDetector, Coll>(candidate, DataType::McDsPrompt);
            fillHistoPiKK<centDetector, useMl, Coll>(candidate, DataType::McDsPrompt);
          }

          std::get<TH1_ptr>(histosPtr[DataType::McDsPrompt]["hPtRecSig"])->Fill(pt);
          std::get<TH1_ptr>(histosPtr[DataType::McDsPrompt]["hPtGenSig"])->Fill(particleMother.pt()); // gen. level pT
          std::get<TH2_ptr>(histosPtr[DataType::McDsPrompt]["hPtVsYRecSigRecoSkim"])->Fill(pt, y);
          std::get<TH1_ptr>(histosPtr[DataType::McDsPrompt]["hCPARecSig"])->Fill(candidate.cpa());
          std::get<TH1_ptr>(histosPtr[DataType::McDsPrompt]["hEtaRecSig"])->Fill(candidate.eta());
          if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
            std::get<TH2_ptr>(histosPtr[DataType::McDsPrompt]["hPtVsYRecSigRecoTopol"])->Fill(pt, y);
          }
          if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
            std::get<TH2_ptr>(histosPtr[DataType::McDsPrompt]["hPtVsYRecSigRecoPID"])->Fill(pt, y);
          }
        }

        // non-prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
            fillHisto<centDetector, Coll>(candidate, DataType::McDsNonPrompt);
            fillHistoKKPi<centDetector, useMl, Coll>(candidate, DataType::McDsNonPrompt);
          }
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
            fillHisto<centDetector, Coll>(candidate, DataType::McDsNonPrompt);
            fillHistoPiKK<centDetector, useMl, Coll>(candidate, DataType::McDsNonPrompt);
          }

          std::get<TH1_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtRecSig"])->Fill(pt);
          std::get<TH1_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtGenSig"])->Fill(particleMother.pt()); // gen. level pT
          std::get<TH2_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtVsYRecSigRecoSkim"])->Fill(pt, y);
          std::get<TH1_ptr>(histosPtr[DataType::McDsNonPrompt]["hCPARecSig"])->Fill(candidate.cpa());
          std::get<TH1_ptr>(histosPtr[DataType::McDsNonPrompt]["hEtaRecSig"])->Fill(candidate.eta());
          if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
            std::get<TH2_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtVsYRecSigRecoTopol"])->Fill(pt, y);
          }
          if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
            std::get<TH2_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtVsYRecSigRecoPID"])->Fill(pt, y);
          }
        }
        return;
      } // end Ds

      // D+→ K± K∓ π±
      if (whichSpeciesDecay == SpeciesAndDecay::DplusToKKPi) {
        double y = hfHelper.yDplus(candidate);

        // prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
            fillHisto<centDetector, Coll>(candidate, DataType::McDplusPrompt);
            fillHistoKKPi<centDetector, useMl, Coll>(candidate, DataType::McDplusPrompt);
          }
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
            fillHisto<centDetector, Coll>(candidate, DataType::McDplusPrompt);
            fillHistoPiKK<centDetector, useMl, Coll>(candidate, DataType::McDplusPrompt);
          }

          std::get<TH1_ptr>(histosPtr[DataType::McDplusPrompt]["hPtRecSig"])->Fill(pt);
          std::get<TH1_ptr>(histosPtr[DataType::McDplusPrompt]["hPtGenSig"])->Fill(particleMother.pt()); // gen. level pT
          std::get<TH2_ptr>(histosPtr[DataType::McDplusPrompt]["hPtVsYRecSigRecoSkim"])->Fill(pt, y);
          std::get<TH1_ptr>(histosPtr[DataType::McDplusPrompt]["hCPARecSig"])->Fill(candidate.cpa());
          std::get<TH1_ptr>(histosPtr[DataType::McDplusPrompt]["hEtaRecSig"])->Fill(candidate.eta());
          if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
            std::get<TH2_ptr>(histosPtr[DataType::McDplusPrompt]["hPtVsYRecSigRecoTopol"])->Fill(pt, y);
          }
          if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
            std::get<TH2_ptr>(histosPtr[DataType::McDplusPrompt]["hPtVsYRecSigRecoPID"])->Fill(pt, y);
          }
        }

        // non-prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
            fillHisto<centDetector, Coll>(candidate, DataType::McDplusNonPrompt);
            fillHistoKKPi<centDetector, useMl, Coll>(candidate, DataType::McDplusNonPrompt);
          }
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
            fillHisto<centDetector, Coll>(candidate, DataType::McDplusNonPrompt);
            fillHistoPiKK<centDetector, useMl, Coll>(candidate, DataType::McDplusNonPrompt);
          }

          std::get<TH1_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtRecSig"])->Fill(pt);
          std::get<TH1_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtGenSig"])->Fill(particleMother.pt()); // gen. level pT
          std::get<TH2_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtVsYRecSigRecoSkim"])->Fill(pt, y);
          std::get<TH1_ptr>(histosPtr[DataType::McDplusNonPrompt]["hCPARecSig"])->Fill(candidate.cpa());
          std::get<TH1_ptr>(histosPtr[DataType::McDplusNonPrompt]["hEtaRecSig"])->Fill(candidate.eta());
          if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
            std::get<TH2_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtVsYRecSigRecoTopol"])->Fill(pt, y);
          }
          if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
            std::get<TH2_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtVsYRecSigRecoPID"])->Fill(pt, y);
          }
        }

        return;
      } // end D+→ K± K∓ π±

      // D+→ π± K∓ π±
      if (whichSpeciesDecay == SpeciesAndDecay::DplusToPiKPi) {
        double y = hfHelper.yDplus(candidate);

        // Fill whether it is prompt or non-prompt
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
          fillHisto<centDetector, Coll>(candidate, DataType::McDplusBkg);
          fillHistoKKPi<centDetector, useMl, Coll>(candidate, DataType::McDplusBkg);
        }
        if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
          fillHisto<centDetector, Coll>(candidate, DataType::McDplusBkg);
          fillHistoPiKK<centDetector, useMl, Coll>(candidate, DataType::McDplusBkg);
        }

        std::get<TH1_ptr>(histosPtr[DataType::McDplusBkg]["hPtRecSig"])->Fill(pt);
        std::get<TH1_ptr>(histosPtr[DataType::McDplusBkg]["hPtGenSig"])->Fill(particleMother.pt()); // gen. level pT
        std::get<TH2_ptr>(histosPtr[DataType::McDplusBkg]["hPtVsYRecSigRecoSkim"])->Fill(pt, y);
        std::get<TH1_ptr>(histosPtr[DataType::McDplusBkg]["hCPARecSig"])->Fill(candidate.cpa());
        std::get<TH1_ptr>(histosPtr[DataType::McDplusBkg]["hEtaRecSig"])->Fill(candidate.eta());
        if (TESTBIT(flag, aod::SelectionStep::RecoTopol)) {
          std::get<TH2_ptr>(histosPtr[DataType::McDplusBkg]["hPtVsYRecSigRecoTopol"])->Fill(pt, y);
        }
        if (TESTBIT(flag, aod::SelectionStep::RecoPID)) {
          std::get<TH2_ptr>(histosPtr[DataType::McDplusBkg]["hPtVsYRecSigRecoPID"])->Fill(pt, y);
        }

        return;
      } // end D+→ π± K∓ π±
    }
  }

  template <FinalState decayChannel, CentralityEstimator centDetector, bool useMl, typename Coll, typename CandsDs>
  void runDataAnalysis(CandsDs const& candidates)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
        continue;
      }
      if constexpr (decayChannel == FinalState::KKPi) { // KKPi
        fillHisto<centDetector, Coll>(candidate, DataType::Data);
        fillHistoKKPi<centDetector, useMl, Coll>(candidate, DataType::Data);
      } else if constexpr (decayChannel == FinalState::PiKK) { // PiKK
        fillHisto<centDetector, Coll>(candidate, DataType::Data);
        fillHistoPiKK<centDetector, useMl, Coll>(candidate, DataType::Data);
      }
    }
  }

  template <CentralityEstimator centDetector, bool useMl, typename Coll, typename CandsDs>
  void runMcAnalysis(CandsDs const& /*candidates*/,
                     CandDsMcGen const& mcParticles,
                     Coll const& recoCollisions)
  {
    // MC rec.
    // Ds± → K± K∓ π±
    if constexpr (useMl) {
      // Ds
      for (const auto& candidate : reconstructedCandDsSigWithMl)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl, Coll>(candidate, mcParticles, SpeciesAndDecay::DsToKKPi);
      // D+→ K± K∓ π±
      for (const auto& candidate : reconstructedCandDplusSigWithMl)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl, Coll>(candidate, mcParticles, SpeciesAndDecay::DplusToKKPi);
      // D+→ π± K∓ π±
      for (const auto& candidate : reconstructedCandDplusBkgWithMl)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl, Coll>(candidate, mcParticles, SpeciesAndDecay::DplusToPiKPi);
      // Bkg
      for (const auto& candidate : reconstructedCandBkgWithMl) {
        if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
          continue;
        }

        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHisto<centDetector, Coll>(candidate, DataType::McBkg);
      }
    } else {
      // Ds
      for (const auto& candidate : reconstructedCandDsSig)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl, Coll>(candidate, mcParticles, SpeciesAndDecay::DsToKKPi);

      // D+→ K± K∓ π±
      for (const auto& candidate : reconstructedCandDplusSig)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl, Coll>(candidate, mcParticles, SpeciesAndDecay::DplusToKKPi);

      // D+→ π± K∓ π±
      for (const auto& candidate : reconstructedCandDplusBkg)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl, Coll>(candidate, mcParticles, SpeciesAndDecay::DplusToPiKPi);

      // Bkg
      for (const auto& candidate : reconstructedCandBkg) {
        if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
          continue;
        }

        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs) {
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
            fillHisto<centDetector, Coll>(candidate, DataType::McBkg);
            fillHistoKKPi<centDetector, useMl, Coll>(candidate, DataType::McDsPrompt);
          }
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
            fillHisto<centDetector, Coll>(candidate, DataType::McBkg);
            fillHistoPiKK<centDetector, useMl, Coll>(candidate, DataType::McDsPrompt);
          }
        }
      }
    }

    // TODO: add histograms for reflections

    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        if (particle.flagMcDecayChanGen() == decayChannel || (fillDplusMc && particle.flagMcDecayChanGen() == (decayChannel + offsetDplusDecayChannel))) {
          auto pt = particle.pt();
          double y{0.f};

          unsigned maxNumContrib = 0; // Search for reco. collisions of the same MC collision

          if constexpr (centDetector == CentralityEstimator::None) {
            const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
            for (const auto& recCol : recoCollsPerMcColl) {
              maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
            }
          } else if constexpr (centDetector == CentralityEstimator::FT0C) {
            const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollisionWithFT0C, particle.mcCollision().globalIndex());
            for (const auto& recCol : recoCollsPerMcColl) {
              maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
            }
          } else if constexpr (centDetector == CentralityEstimator::FT0M) {
            const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollisionWithFT0M, particle.mcCollision().globalIndex());
            for (const auto& recCol : recoCollsPerMcColl) {
              maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
            }
          } else if constexpr (centDetector == CentralityEstimator::NTracksPV) {
            const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollisionWithNTracksPV, particle.mcCollision().globalIndex());
            for (const auto& recCol : recoCollsPerMcColl) {
              maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
            }
          }

          if (particle.flagMcDecayChanGen() == decayChannel) {
            y = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDS);
            if (yCandGenMax >= 0. && std::abs(y) > yCandGenMax) {
              continue;
            }

            if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
              std::get<TH1_ptr>(histosPtr[DataType::McDsPrompt]["hPtGen"])->Fill(pt); // gen. level pT
              std::get<TH1_ptr>(histosPtr[DataType::McDsPrompt]["hEtaGen"])->Fill(particle.eta());
              std::get<THnSparse_ptr>(histosPtr[DataType::McDsPrompt]["hPtYNPvContribGen"])->Fill(pt, y, maxNumContrib); // gen. level pT
            }
            if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              std::get<TH1_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtGen"])->Fill(pt);                                    // gen. level pT
              std::get<TH1_ptr>(histosPtr[DataType::McDsNonPrompt]["hEtaGen"])->Fill(particle.eta());                       // gen. level pT
              std::get<THnSparse_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtYNPvContribGen"])->Fill(pt, y, maxNumContrib); // gen. level pT
            }
          } else if (fillDplusMc) {
            y = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus);
            if (yCandGenMax >= 0. && std::abs(y) > yCandGenMax) {
              continue;
            }
            if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
              std::get<TH1_ptr>(histosPtr[DataType::McDplusPrompt]["hPtGen"])->Fill(pt); // gen. level pT
              std::get<TH1_ptr>(histosPtr[DataType::McDplusPrompt]["hEtaGen"])->Fill(particle.eta());
              std::get<THnSparse_ptr>(histosPtr[DataType::McDplusPrompt]["hPtYNPvContribGen"])->Fill(pt, y, maxNumContrib); // gen. level pT
            }
            if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              std::get<TH1_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtGen"])->Fill(pt); // gen. level pT
              std::get<TH1_ptr>(histosPtr[DataType::McDplusNonPrompt]["hEtaGen"])->Fill(particle.eta());
              std::get<THnSparse_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtYNPvContribGen"])->Fill(pt, y, maxNumContrib); // gen. level pT
            }
          }
        }
      }
    }
  }

  /// Checks wheter the candidate is in the signal region of either the Ds or D+ decay
  /// \param candidate is the candidate
  /// \param isDs is true if we check for the Ds signal region, false for the D+ signal region
  /// \return true if the candidate is in the signal region, false otherwise
  template <typename CandDs>
  bool checkCandInSignalRegion(const CandDs& candidate, bool isDs)
  {
    bool isKKPi = candidate.isSelDsToKKPi() >= selectionFlagDs;
    float invMass = isKKPi ? hfHelper.invMassDsToKKPi(candidate) : hfHelper.invMassDsToPiKK(candidate);
    if (isDs && (invMass < minMassDsSignal || invMass > minMassDsSignal))
      return false;
    if (!isDs && (invMass < minMassDplusSignal || invMass > maxMassDplusSignal))
      return false;
    return true;
  }

  template <typename CandDs>
  bool isDsPrompt(const CandDs& candidate)
  {
    return abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && candidate.flagMcDecayChanRec() == decayChannel && candidate.originMcRec() == RecoDecay::OriginType::Prompt;
  }

  template <typename CandDs>
  bool isDplusPrompt(const CandDs& candidate)
  {
    return abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && candidate.flagMcDecayChanRec() == decayChannel + offsetDplusDecayChannel && candidate.originMcRec() == RecoDecay::OriginType::Prompt;
  }

  template <typename CandDs>
  bool isDsNonPrompt(const CandDs& candidate)
  {
    return abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && candidate.flagMcDecayChanRec() == decayChannel && candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
  }

  template <typename CandDs>
  bool isDplusNonPrompt(const CandDs& candidate)
  {
    return abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && candidate.flagMcDecayChanRec() == decayChannel + offsetDplusDecayChannel && candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
  }

  template <CentralityEstimator centDetector, bool doMc, bool useMl, typename Coll, typename CandsDs>
  void fillNPvContribHisto(const Coll& collisions, const CandsDs& candidates)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      std::array<int, DataType::kDataTypes> nCandsPerType{0};
      std::array<int, DataType::kDataTypes> nCandsInSignalRegionDsPerType{0};
      std::array<int, DataType::kDataTypes> nCandsInSignalRegionDplusPerType{0};
      if constexpr (doMc) {
        if constexpr (useMl) {
          auto groupedDsCandidates = candidates.sliceBy(candDsMcRecoWithMlPerCollision, thisCollId);
          for (const auto& candidate : groupedDsCandidates) {
            if (candidate.isSelDsToKKPi() < selectionFlagDs && candidate.isSelDsToPiKK() < selectionFlagDs)
              continue;

            if (isDsPrompt(candidate)) {
              ++nCandsPerType[DataType::McDsPrompt];
              nCandsInSignalRegionDsPerType[DataType::McDsPrompt] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McDsPrompt] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            } else if (isDsNonPrompt(candidate)) {
              ++nCandsPerType[DataType::McDsNonPrompt];
              nCandsInSignalRegionDsPerType[DataType::McDsNonPrompt] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McDsNonPrompt] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            } else if (isDplusPrompt(candidate)) {
              ++nCandsPerType[DataType::McDplusPrompt];
              nCandsInSignalRegionDsPerType[DataType::McDplusPrompt] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McDplusPrompt] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            } else if (isDplusNonPrompt(candidate)) {
              ++nCandsPerType[DataType::McDplusNonPrompt];
              nCandsInSignalRegionDsPerType[DataType::McDplusNonPrompt] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McDplusNonPrompt] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            } else {
              ++nCandsPerType[DataType::McBkg];
              nCandsInSignalRegionDsPerType[DataType::McBkg] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McBkg] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            }
          }
        } else {
          auto groupedDsCandidates = candidates.sliceBy(candDsMcRecoPerCollision, thisCollId);
          for (const auto& candidate : groupedDsCandidates) {
            if (candidate.isSelDsToKKPi() < selectionFlagDs && candidate.isSelDsToPiKK() < selectionFlagDs)
              continue;

            if (isDsPrompt(candidate)) {
              ++nCandsPerType[DataType::McDsPrompt];
              nCandsInSignalRegionDsPerType[DataType::McDsPrompt] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McDsPrompt] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            } else if (isDsNonPrompt(candidate)) {
              ++nCandsPerType[DataType::McDsNonPrompt];
              nCandsInSignalRegionDsPerType[DataType::McDsNonPrompt] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McDsNonPrompt] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            } else if (isDplusPrompt(candidate)) {
              ++nCandsPerType[DataType::McDplusPrompt];
              nCandsInSignalRegionDsPerType[DataType::McDplusPrompt] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McDplusPrompt] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            } else if (isDplusNonPrompt(candidate)) {
              ++nCandsPerType[DataType::McDplusNonPrompt];
              nCandsInSignalRegionDsPerType[DataType::McDplusNonPrompt] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McDplusNonPrompt] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            } else {
              ++nCandsPerType[DataType::McBkg];
              nCandsInSignalRegionDsPerType[DataType::McBkg] += checkCandInSignalRegion(candidate, true) ? 1 : 0;
              nCandsInSignalRegionDplusPerType[DataType::McBkg] += checkCandInSignalRegion(candidate, false) ? 1 : 0;
            }
          }
        }
        nCandsPerType[DataType::Data] = nCandsPerType[DataType::McDsPrompt] + nCandsPerType[DataType::McDsNonPrompt] + nCandsPerType[DataType::McDplusPrompt] + nCandsPerType[DataType::McDplusNonPrompt] + nCandsPerType[DataType::McBkg];
        nCandsInSignalRegionDsPerType[DataType::Data] = nCandsInSignalRegionDsPerType[DataType::McDsPrompt] + nCandsInSignalRegionDsPerType[DataType::McDsNonPrompt] + nCandsInSignalRegionDsPerType[DataType::McDplusPrompt] + nCandsInSignalRegionDsPerType[DataType::McDplusNonPrompt] + nCandsInSignalRegionDsPerType[DataType::McBkg];
        nCandsInSignalRegionDplusPerType[DataType::Data] = nCandsInSignalRegionDplusPerType[DataType::McDsPrompt] + nCandsInSignalRegionDplusPerType[DataType::McDsNonPrompt] + nCandsInSignalRegionDplusPerType[DataType::McDplusPrompt] + nCandsInSignalRegionDplusPerType[DataType::McDplusNonPrompt] + nCandsInSignalRegionDplusPerType[DataType::McBkg];
      } else {
        if constexpr (useMl) {
          auto groupedDsCandidates = candidates.sliceBy(candDsDataWithMlPerCollision, thisCollId);
          for (const auto& candidate : groupedDsCandidates) {
            if (candidate.isSelDsToKKPi() < selectionFlagDs && candidate.isSelDsToPiKK() < selectionFlagDs)
              continue;

            ++nCandsPerType[DataType::Data];
            if (checkCandInSignalRegion(candidate, true))
              ++nCandsInSignalRegionDsPerType[DataType::Data];
            if (checkCandInSignalRegion(candidate, false))
              ++nCandsInSignalRegionDplusPerType[DataType::Data];
          }
        } else {
          auto groupedDsCandidates = candidates.sliceBy(candDsDataPerCollision, thisCollId);
          for (const auto& candidate : groupedDsCandidates) {
            if (candidate.isSelDsToKKPi() < selectionFlagDs && candidate.isSelDsToPiKK() < selectionFlagDs)
              continue;

            ++nCandsPerType[DataType::Data];
            if (checkCandInSignalRegion(candidate, true))
              ++nCandsInSignalRegionDsPerType[DataType::Data];
            if (checkCandInSignalRegion(candidate, false))
              ++nCandsInSignalRegionDplusPerType[DataType::Data];
          }
        }
      }
      std::get<TH2_ptr>(histosPtr[DataType::Data]["hNPvContribAll"])->Fill(numPvContributors, evaluateCentralityColl<centDetector>(collision));
      for (int i = 0; i < DataType::kDataTypes; i++) {
        if (nCandsPerType[i])
          std::get<TH2_ptr>(histosPtr[i]["hNPvContribCands"])->Fill(numPvContributors, evaluateCentralityColl<centDetector>(collision));
        if (nCandsInSignalRegionDsPerType[i])
          std::get<TH2_ptr>(histosPtr[i]["hNPvContribCandsInSignalRegionDs"])->Fill(numPvContributors, evaluateCentralityColl<centDetector>(collision));
        if (nCandsInSignalRegionDplusPerType[i])
          std::get<TH2_ptr>(histosPtr[i]["hNPvContribCandsInSignalRegionDplus"])->Fill(numPvContributors, evaluateCentralityColl<centDetector>(collision));
      }
    }
  }

  void processDataWithCentFT0C(CollisionsWithFT0C const& collisions,
                               CandDsData const& candsDs,
                               aod::Tracks const&)
  {
    fillNPvContribHisto<CentralityEstimator::FT0C, false /*doMC*/, false /*useML*/, CollisionsWithFT0C>(collisions, candsDs);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0C, false, CollisionsWithFT0C>(selectedDsToKKPiCandData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0C, false, CollisionsWithFT0C>(selectedDsToPiKKCandData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentFT0C, "Process data w/o ML information on Ds, with information on centrality from FT0C", false);

  void processDataWithCentFT0M(CollisionsWithFT0M const& collisions,
                               CandDsData const& candsDs,
                               aod::Tracks const&)
  {
    fillNPvContribHisto<CentralityEstimator::FT0M, false /*doMC*/, false /*useML*/, CollisionsWithFT0M>(collisions, candsDs);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0M, false, CollisionsWithFT0M>(selectedDsToKKPiCandData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0M, false, CollisionsWithFT0M>(selectedDsToPiKKCandData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentFT0M, "Process data w/o ML information on Ds, with information on centrality from FT0M", false);

  void processDataWithCentNTracksPV(CollisionsWithNTracksPV const& collisions,
                                    CandDsData const& candsDs,
                                    aod::Tracks const&)
  {
    fillNPvContribHisto<CentralityEstimator::NTracksPV, false /*doMC*/, false /*useML*/, CollisionsWithNTracksPV>(collisions, candsDs);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::NTracksPV, false, CollisionsWithNTracksPV>(selectedDsToKKPiCandData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::NTracksPV, false, CollisionsWithNTracksPV>(selectedDsToPiKKCandData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentNTracksPV, "Process data w/o ML information on Ds, with information on centrality from NTracksPV", false);

  void processData(aod::Collisions const&,
                   CandDsData const&,
                   aod::Tracks const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::None, false, aod::Collisions>(selectedDsToKKPiCandData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::None, false, aod::Collisions>(selectedDsToPiKKCandData);
  }
  PROCESS_SWITCH(HfTaskDs, processData, "Process data w/o ML information on Ds, w/o information on centrality", true);

  void processDataWithMlAndCentFT0C(CollisionsWithFT0C const& collisions,
                                    CandDsDataWithMl const& candsDs,
                                    aod::Tracks const&)
  {
    fillNPvContribHisto<CentralityEstimator::FT0C, false /*doMC*/, true /*useML*/, CollisionsWithFT0C>(collisions, candsDs);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0C, true, CollisionsWithFT0C>(selectedDsToKKPiCandWithMlData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0C, true, CollisionsWithFT0C>(selectedDsToPiKKCandWithMlData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentFT0C, "Process data with ML information on Ds, with information on centrality from FT0C", false);

  void processDataWithMlAndCentFT0M(CollisionsWithFT0M const& collisions,
                                    CandDsDataWithMl const& candsDs,
                                    aod::Tracks const&)
  {
    fillNPvContribHisto<CentralityEstimator::FT0M, false /*doMC*/, true /*useML*/, CollisionsWithFT0M>(collisions, candsDs);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0M, true, CollisionsWithFT0M>(selectedDsToKKPiCandWithMlData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0M, true, CollisionsWithFT0M>(selectedDsToPiKKCandWithMlData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentFT0M, "Process data with ML information on Ds, with information on centrality from FT0M", false);

  void processDataWithMlAndCentNTracksPV(CollisionsWithNTracksPV const& collisions,
                                         CandDsDataWithMl const& candsDs,
                                         aod::Tracks const&)
  {
    fillNPvContribHisto<CentralityEstimator::NTracksPV, false /*doMC*/, true /*useML*/, CollisionsWithNTracksPV>(collisions, candsDs);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::NTracksPV, true, CollisionsWithNTracksPV>(selectedDsToKKPiCandWithMlData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::NTracksPV, true, CollisionsWithNTracksPV>(selectedDsToPiKKCandWithMlData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentNTracksPV, "Process data with ML information on Ds, with information on centrality", false);

  void processDataWithMl(soa::Join<aod::Collisions, aod::EvSels> const&,
                         CandDsDataWithMl const&,
                         aod::Tracks const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::None, true, soa::Join<aod::Collisions, aod::EvSels>>(selectedDsToKKPiCandWithMlData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::None, true, soa::Join<aod::Collisions, aod::EvSels>>(selectedDsToPiKKCandWithMlData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMl, "Process data with ML information on Ds, w/o information on centrality", false);

  void processMcWithCentFT0C(CollisionsMcWithFT0C const& recoCollisions,
                             CandDsMcReco const& candsDs,
                             CandDsMcGen const& mcParticles,
                             aod::McCollisions const&,
                             aod::TracksWMc const&)
  {
    fillNPvContribHisto<CentralityEstimator::FT0C, true /*doMC*/, false /*useML*/, CollisionsMcWithFT0C>(recoCollisions, candsDs);
    runMcAnalysis<CentralityEstimator::FT0C, false, CollisionsMcWithFT0C>(candsDs, mcParticles, recoCollisions);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0C, false, CollisionsMcWithFT0C>(selectedDsToKKPiCandMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0C, false, CollisionsMcWithFT0C>(selectedDsToPiKKCandMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentFT0C, "Process MC w/o ML information on Ds, with information on centrality from FT0C", false);

  void processMcWithCentFT0M(CollisionsMcWithFT0M const& recoCollisions,
                             CandDsMcReco const& candsDs,
                             CandDsMcGen const& mcParticles,
                             aod::McCollisions const&,
                             aod::TracksWMc const&)
  {
    fillNPvContribHisto<CentralityEstimator::FT0M, true /*doMC*/, false /*useML*/, CollisionsMcWithFT0M>(recoCollisions, candsDs);
    runMcAnalysis<CentralityEstimator::FT0M, false, CollisionsMcWithFT0M>(candsDs, mcParticles, recoCollisions);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0M, false, CollisionsMcWithFT0M>(selectedDsToKKPiCandMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0M, false, CollisionsMcWithFT0M>(selectedDsToPiKKCandMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentFT0M, "Process MC w/o ML information on Ds, with information on centrality from FT0M", false);

  void processMcWithCentNTracksPV(CollisionsMcWithNTracksPV const& recoCollisions,
                                  CandDsMcReco const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    fillNPvContribHisto<CentralityEstimator::NTracksPV, true /*doMC*/, false /*useML*/, CollisionsMcWithNTracksPV>(recoCollisions, candsDs);
    runMcAnalysis<CentralityEstimator::NTracksPV, false, CollisionsMcWithNTracksPV>(candsDs, mcParticles, recoCollisions);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::NTracksPV, false, CollisionsMcWithNTracksPV>(selectedDsToKKPiCandMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::NTracksPV, false, CollisionsMcWithNTracksPV>(selectedDsToPiKKCandMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentNTracksPV, "Process MC w/o ML information on Ds, with information on centrality from NTracksPV", false);

  void processMc(CollisionsMc const& recoCollisions,
                 CandDsMcReco const& candsDs,
                 CandDsMcGen const& mcParticles,
                 aod::McCollisions const&,
                 aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::None, false, CollisionsMc>(candsDs, mcParticles, recoCollisions);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::None, false, CollisionsMc>(selectedDsToKKPiCandMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::None, false, CollisionsMc>(selectedDsToPiKKCandMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMc, "Process MC w/o ML information on Ds, w/o information on centrality", false);

  void processMcWithMlAndCentFT0C(CollisionsMcWithFT0C const& recoCollisions,
                                  CandDsMcRecoWithMl const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    fillNPvContribHisto<CentralityEstimator::FT0C, true /*doMC*/, true /*useML*/, CollisionsMcWithFT0C>(recoCollisions, candsDs);
    runMcAnalysis<CentralityEstimator::FT0C, true, CollisionsMcWithFT0C>(candsDs, mcParticles, recoCollisions);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0C, true, CollisionsMcWithFT0C>(selectedDsToKKPiCandWithMlMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0C, true, CollisionsMcWithFT0C>(selectedDsToPiKKCandWithMlMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentFT0C, "Process MC with ML information on Ds, with information on centrality from FT0C", false);

  void processMcWithMlAndCentFT0M(CollisionsMcWithFT0M const& recoCollisions,
                                  CandDsMcRecoWithMl const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    fillNPvContribHisto<CentralityEstimator::FT0M, true /*doMC*/, true /*useML*/, CollisionsMcWithFT0M>(recoCollisions, candsDs);
    runMcAnalysis<CentralityEstimator::FT0M, true, CollisionsMcWithFT0M>(candsDs, mcParticles, recoCollisions);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0M, true, CollisionsMcWithFT0M>(selectedDsToKKPiCandWithMlMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0M, true, CollisionsMcWithFT0M>(selectedDsToPiKKCandWithMlMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentFT0M, "Process MC with ML information on Ds, with information on centrality from FT0M", false);

  void processMcWithMlAndCentNTracksPV(CollisionsMcWithNTracksPV const& recoCollisions,
                                       CandDsMcRecoWithMl const& candsDs,
                                       CandDsMcGen const& mcParticles,
                                       aod::McCollisions const&,
                                       aod::TracksWMc const&)
  {
    fillNPvContribHisto<CentralityEstimator::NTracksPV, true /*doMC*/, true /*useML*/, CollisionsMcWithNTracksPV>(recoCollisions, candsDs);
    runMcAnalysis<CentralityEstimator::NTracksPV, true, CollisionsMcWithNTracksPV>(candsDs, mcParticles, recoCollisions);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::NTracksPV, true, CollisionsMcWithNTracksPV>(selectedDsToKKPiCandWithMlMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::NTracksPV, true, CollisionsMcWithNTracksPV>(selectedDsToPiKKCandWithMlMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentNTracksPV, "Process MC with ML information on Ds, with information on centrality from NTracksPV", false);

  void processMcWithMl(CollisionsMc const& recoCollisions,
                       CandDsMcRecoWithMl const& candsDs,
                       CandDsMcGen const& mcParticles,
                       aod::McCollisions const&,
                       aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::None, true, CollisionsMc>(candsDs, mcParticles, recoCollisions);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::None, true, CollisionsMc>(selectedDsToKKPiCandWithMlMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::None, true, CollisionsMc>(selectedDsToPiKKCandWithMlMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMl, "Process MC with ML information on Ds, w/o information on centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDs>(cfgc)};
}
