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

template <typename T>
concept hasDsMlInfo = requires(T candidate)
{
  candidate.mlProbDsToKKPi();
  candidate.mlProbDsToPiKK();
};

/// Ds± analysis task
struct HfTaskDs {

  Configurable<int> decayChannel{"decayChannel", 1, "Switch between decay channels: 1 for Ds/Dplus->PhiPi->KKpi, 2 for Ds/Dplus->K0*K->KKPi"};
  Configurable<bool> fillDplusMc{"fillDplusMc", true, "Switch to fill Dplus MC information"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2, 3}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> massDsSignalMin{"massDsSignalMin", 1.934, "min mass for Ds signal"};
  Configurable<float> massDsSignalMax{"massDsSignalMax", 1.994, "max mass for Ds signal"};
  Configurable<float> massDplusSignalMin{"massDplusSignalMin", 1.866, "min mass for Dplus signal"};
  Configurable<float> massDplusSignalMax{"massDplusSignalMax", 1.906, "max mass for Dplus signal"};

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
  template <typename CandDs>
  using MemberFunctionPointer = bool (HfTaskDs::*)(const CandDs&);

  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsWithNTracksPV = soa::Join<aod::Collisions, aod::EvSels, aod::CentNTPVs>;

  using CollisionsMc = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using CollisionsMcWithFT0C = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsMcWithFT0M = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsMcWithNTracksPV = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentNTPVs>;

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsDataWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcRecoWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec, aod::HfMlDsToKKPi>>;
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Preslice<CandDsData> candDsDataPerCollision = aod::hf_cand::collisionId;
  Preslice<CandDsDataWithMl> candDsDataWithMlPerCollision = aod::hf_cand::collisionId;
  Preslice<CandDsMcReco> candDsMcRecoPerCollision = aod::hf_cand::collisionId;
  Preslice<CandDsMcRecoWithMl> candDsMcRecoWithMlPerCollision = aod::hf_cand::collisionId;

  PresliceUnsorted<CollisionsMc> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsMcWithFT0C> colPerMcCollisionWithFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsMcWithFT0M> colPerMcCollisionWithFT0M = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsMcWithNTracksPV> colPerMcCollisionWithNTracksPV = aod::mccollisionlabel::mcCollisionId;
  SliceCache cache;

  int offsetDplusDecayChannel = aod::hf_cand_3prong::DecayChannelDToKKPi::DplusToPhiPi - aod::hf_cand_3prong::DecayChannelDToKKPi::DsToPhiPi; // Offset between Dplus and Ds to use the same decay channel. See aod::hf_cand_3prong::DecayChannelDToKKPi

  Filter filterDsFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi))) != static_cast<uint8_t>(0);

  HistogramRegistry registry{"registry", {}};

  std::array<std::string, DataType::kDataTypes> folders = {"Data/", "MC/Ds/Prompt/", "MC/Ds/NonPrompt/", "MC/Dplus/Prompt/", "MC/Dplus/NonPrompt/", "MC/Dplus/Bkg/", "MC/Bkg/"};

  std::unordered_map<std::string, histTypes> dataHistograms = {};
  std::unordered_map<std::string, histTypes> mcDsPromptHistograms = {};
  std::unordered_map<std::string, histTypes> mcDsNonPromptHistograms = {};
  std::unordered_map<std::string, histTypes> mcDplusPromptHistograms = {};
  std::unordered_map<std::string, histTypes> mcDplusNonPromptHistograms = {};
  std::unordered_map<std::string, histTypes> mcDplusBkgHistograms = {};
  std::unordered_map<std::string, histTypes> mcBkgHistograms = {};

  std::map<CentralityEstimator, std::variant<PresliceUnsorted<CollisionsMc>, PresliceUnsorted<CollisionsMcWithFT0C>, PresliceUnsorted<CollisionsMcWithFT0M>, PresliceUnsorted<CollisionsMcWithNTracksPV>>> colPerMcCollisionMap{
    {CentralityEstimator::None, colPerMcCollision},
    {CentralityEstimator::FT0C, colPerMcCollisionWithFT0C},
    {CentralityEstimator::FT0M, colPerMcCollisionWithFT0M},
    {CentralityEstimator::NTracksPV, colPerMcCollisionWithNTracksPV}};

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
      if (doprocessDataWithCentFT0C || doprocessDataWithCentFT0M || doprocessDataWithCentNTracksPV || doprocessData || doprocessMcWithCentFT0C || doprocessMcWithCentFT0M || doprocessMcWithCentNTracksPV || doprocessMc) {
        if (i == DataType::Data) { // If data do not fill PV contributors in sparse
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins});
        } else {
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins, axisNPvContributors});
        }
      } else if (doprocessDataWithMlAndCentFT0C || doprocessDataWithMlAndCentFT0M || doprocessDataWithMlAndCentNTracksPV || doprocessDataWithMl || doprocessMcWithMlAndCentFT0C || doprocessMcWithMlAndCentFT0M || doprocessMcWithMlAndCentNTracksPV || doprocessMcWithMl) {
        if (i == DataType::Data) { // If data do not fill PV contributors in sparse
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins, axisMlScore0, axisMlScore1, axisMlScore2});
        } else {
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins, axisMlScore0, axisMlScore1, axisMlScore2, axisNPvContributors});
        }
      }
      histosPtr[i]["hNPvContribCands"] = registry.add<TH2>((folders[i] + "hNPvContribCands").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
      histosPtr[i]["hNPvContribCandsInSignalRegionDs"] = registry.add<TH2>((folders[i] + "hNPvContribCandsInSignalRegionDs").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
      histosPtr[i]["hNPvContribCandsInSignalRegionDplus"] = registry.add<TH2>((folders[i] + "hNPvContribCandsInSignalRegionDplus").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, centralitybins});
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
          histosPtr[i]["hPtGen"] = registry.add<TH1>((folders[i] + "hPtGen").c_str(), "MC particles (unmatched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {ptbins}});
          histosPtr[i]["hPtVsYRecoPID"] = registry.add<TH2>((folders[i] + "hPtVsYRecoPID").c_str(), "3-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
          histosPtr[i]["hPtVsYRecoTopol"] = registry.add<TH2>((folders[i] + "hPtVsYRecoTopol").c_str(), "3-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
          histosPtr[i]["hPtVsYRecoSkim"] = registry.add<TH2>((folders[i] + "hPtVsYRecoSkim").c_str(), "3-prong candidates (RecoSkim - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
          histosPtr[i]["hPtYNPvContribGen"] = registry.add<THnSparse>((folders[i] + "hPtYNPvContribGen").c_str(), "Thn for generated candidates", {HistType::kTHnSparseF, {ptbins, {ybins}, axisNPvContributors}});
        }
      }
    }
  }

  template <typename CandDs>
  bool isDsPrompt(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && candidate.flagMcDecayChanRec() == decayChannel && candidate.originMcRec() == RecoDecay::OriginType::Prompt;
  }

  template <typename CandDs>
  bool isDplusPrompt(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && candidate.flagMcDecayChanRec() == decayChannel + offsetDplusDecayChannel && candidate.originMcRec() == RecoDecay::OriginType::Prompt;
  }

  template <typename CandDs>
  bool isDsNonPrompt(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && candidate.flagMcDecayChanRec() == decayChannel && candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
  }

  template <typename CandDs>
  bool isDplusNonPrompt(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && candidate.flagMcDecayChanRec() == decayChannel + offsetDplusDecayChannel && candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
  }

  template <typename CandDs>
  bool isDplusBkg(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi));
  }

  /// Checks whether the candidate is in the signal region of either the Ds or D+ decay
  /// \param candidate is the candidate
  /// \param isDs is true if we check for the Ds signal region, false for the D+ signal region
  /// \return true if the candidate is in the signal region, false otherwise
  template <typename CandDs>
  bool isCandInSignalRegion(const CandDs& candidate, bool isDs)
  {
    bool isKKPi = candidate.isSelDsToKKPi() >= selectionFlagDs;
    float invMass = isKKPi ? hfHelper.invMassDsToKKPi(candidate) : hfHelper.invMassDsToPiKK(candidate);
    if (isDs && (invMass < massDsSignalMin || invMass > massDsSignalMax)) {
      return false;
    }
    if (!isDs && (invMass < massDplusSignalMin || invMass > massDplusSignalMax)) {
      return false;
    }
    return true;
  }

  /// Evaluate centrality/multiplicity percentile (centrality estimator is automatically selected based on the used table)
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll>
  float evaluateCentralityColl(const Coll& collision)
  {
    return o2::hf_centrality::getCentralityColl<Coll>(collision);
  }

  /// Evaluate centrality/multiplicity percentile
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision associated to the candidate
  template <typename Coll, typename T1>
  float evaluateCentralityCand(const T1& candidate)
  {
    return evaluateCentralityColl(candidate.template collision_as<Coll>());
  }

  /// Fill histograms of quantities independent from the daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <typename T1>
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

  /// Fill mass sparse if ML information is present
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  /// \param finalState is either KKPi or PiKK, as defined in FinalState enum
  template <typename Coll, hasDsMlInfo Cand>
  void fillSparse(const Cand& candidate, DataType dataType, FinalState finalState)
  {
    auto mass = finalState == FinalState::KKPi ? hfHelper.invMassDsToKKPi(candidate) : hfHelper.invMassDsToPiKK(candidate);
    auto pt = candidate.pt();
    auto mlScore = finalState == FinalState::KKPi ? candidate.mlProbDsToKKPi() : candidate.mlProbDsToPiKK();

    std::vector<float> outputMl = {-999., -999., -999.};
    for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) { // TODO: add checks for classMl size
      if (mlScore.size() == 0) {
        continue;
      }
      outputMl[iclass] = mlScore[classMl->at(iclass)];
    }

    if (dataType == DataType::Data) { // If data do not fill PV contributors in sparse
      std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), outputMl[0], outputMl[1], outputMl[2]);
      return;
    }
    std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib());

    return;
  }

  /// Fill mass sparse if ML information is not present
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  /// \param finalState is either KKPi or PiKK, as defined in FinalState enum
  template <typename Coll, typename Cand>
  void fillSparse(const Cand& candidate, DataType dataType, FinalState finalState)
  {
    auto mass = finalState == FinalState::KKPi ? hfHelper.invMassDsToKKPi(candidate) : hfHelper.invMassDsToPiKK(candidate);
    auto pt = candidate.pt();

    if (dataType == DataType::Data) { // If data do not fill PV contributors in sparse
      std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate));
      return;
    }
    std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), candidate.template collision_as<Coll>().numContrib());

    return;
  }

  /// Fill histograms of quantities for the KKPi daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <typename Coll, typename T1>
  void fillHistoKKPi(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    fillSparse<Coll>(candidate, dataType, FinalState::KKPi);

    std::get<TH2_ptr>(histosPtr[dataType]["hCos3PiK"])->Fill(hfHelper.cos3PiKDsToKKPi(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hAbsCos3PiK"])->Fill(hfHelper.absCos3PiKDsToKKPi(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hDeltaMassPhi"])->Fill(hfHelper.deltaMassPhiDsToKKPi(candidate), pt);
    std::get<TH2_ptr>(histosPtr[dataType]["hMassKK"])->Fill(hfHelper.massKKPairDsToKKPi(candidate), pt);

    return;
  }

  /// Fill histograms of quantities for the PiKK daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <typename Coll, typename T1>
  void fillHistoPiKK(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    fillSparse<Coll>(candidate, dataType, FinalState::PiKK);

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
  template <typename Coll, typename T1>
  void fillHistoMCRec(const T1& candidate, const CandDsMcGen& mcParticles, DataType dataType)
  {

    SpeciesAndDecay whichSpeciesDecay = SpeciesAndDecay::DsToKKPi;
    if (dataType == DataType::McDplusPrompt || dataType == DataType::McDplusNonPrompt) {
      whichSpeciesDecay = SpeciesAndDecay::DplusToKKPi;
    } else if (dataType == DataType::McDplusBkg) {
      whichSpeciesDecay = SpeciesAndDecay::DplusToPiKPi;
    }

    auto indexMother = RecoDecay::getMother(mcParticles,
                                            candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<CandDsMcGen>(),
                                            whichSpeciesDecay == SpeciesAndDecay::DsToKKPi ? o2::constants::physics::Pdg::kDS : o2::constants::physics::Pdg::kDPlus, true);

    if (indexMother != -1) {
      auto yCand = whichSpeciesDecay == SpeciesAndDecay::DsToKKPi ? hfHelper.yDs(candidate) : hfHelper.yDplus(candidate);
      if (yCandRecoMax >= 0. && std::abs(yCand) > yCandRecoMax) {
        return;
      }

      auto pt = candidate.pt(); // rec. level pT

      if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
        fillHisto(candidate, dataType);
        fillHistoKKPi<Coll>(candidate, dataType);

        if (TESTBIT(candidate.isSelDsToKKPi(), aod::SelectionStep::RecoSkims)) {
          std::get<TH2_ptr>(histosPtr[dataType]["hPtVsYRecoSkim"])->Fill(pt, yCand);
        }
        if (TESTBIT(candidate.isSelDsToKKPi(), aod::SelectionStep::RecoTopol)) {
          std::get<TH2_ptr>(histosPtr[dataType]["hPtVsYRecoTopol"])->Fill(pt, yCand);
        }
        if (TESTBIT(candidate.isSelDsToKKPi(), aod::SelectionStep::RecoPID)) {
          std::get<TH2_ptr>(histosPtr[dataType]["hPtVsYRecoPID"])->Fill(pt, yCand);
        }
      }
      if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
        fillHisto(candidate, dataType);
        fillHistoPiKK<Coll>(candidate, dataType);

        if (TESTBIT(candidate.isSelDsToPiKK(), aod::SelectionStep::RecoSkims)) {
          std::get<TH2_ptr>(histosPtr[dataType]["hPtVsYRecoSkim"])->Fill(pt, yCand);
        }
        if (TESTBIT(candidate.isSelDsToPiKK(), aod::SelectionStep::RecoTopol)) {
          std::get<TH2_ptr>(histosPtr[dataType]["hPtVsYRecoTopol"])->Fill(pt, yCand);
        }
        if (TESTBIT(candidate.isSelDsToPiKK(), aod::SelectionStep::RecoPID)) {
          std::get<TH2_ptr>(histosPtr[dataType]["hPtVsYRecoPID"])->Fill(pt, yCand);
        }
      }
    }
    return;
  }

  template <typename Coll, typename CandDs>
  void runDataAnalysisPerCandidate(CandDs const& candidate)
  {
    if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
      return;
    }

    if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
      fillHisto(candidate, DataType::Data);
      fillHistoKKPi<Coll>(candidate, DataType::Data);
    }
    if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
      fillHisto(candidate, DataType::Data);
      fillHistoPiKK<Coll>(candidate, DataType::Data);
    }
  }

  template <typename Coll, typename CandDs>
  void runMcAnalysisPerCandidate(CandDs const& candidate,
                                 CandDsMcGen const& mcParticles)
  {
    // MC rec.
    std::array<MemberFunctionPointer<CandDs>, 5> isOfType = {// Contains the functions to check if the candidate is of a certain type
                                                             &HfTaskDs::isDsPrompt<CandDs>,
                                                             &HfTaskDs::isDsNonPrompt<CandDs>,
                                                             &HfTaskDs::isDplusPrompt<CandDs>,
                                                             &HfTaskDs::isDplusNonPrompt<CandDs>,
                                                             &HfTaskDs::isDplusBkg<CandDs>};

    bool isBkg = true;
    for (int i = DataType::McDsPrompt; i <= DataType::McDplusBkg; i++) { // Check what type of MC signal candidate it is, and fill the corresponding histograms
      if ((this->*isOfType[i - DataType::McDsPrompt])(candidate)) {
        isBkg = false;
        fillHistoMCRec<Coll>(candidate, mcParticles, static_cast<DataType>(i));
        break;
      }
    }
    if (isBkg) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
        return;
      }

      if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs) {
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
          fillHisto(candidate, DataType::McBkg);
          fillHistoKKPi<Coll>(candidate, DataType::McBkg);
        }
        if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
          fillHisto(candidate, DataType::McBkg);
          fillHistoPiKK<Coll>(candidate, DataType::McBkg);
        }
      }
    }

    // TODO: add histograms for reflections
  }

  template <CentralityEstimator centDetector, typename Coll>
  void fillMcGenHistos(CandDsMcGen const& mcParticles,
                       Coll const& recoCollisions)
  {
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        if (particle.flagMcDecayChanGen() == decayChannel || (fillDplusMc && particle.flagMcDecayChanGen() == (decayChannel + offsetDplusDecayChannel))) {
          auto pt = particle.pt();
          double y{0.f};

          unsigned maxNumContrib = 0;
          const auto& recoCollsPerMcColl = recoCollisions.sliceBy(std::get<PresliceUnsorted<Coll>>(colPerMcCollisionMap.at(centDetector)), particle.mcCollision().globalIndex());
          for (const auto& recCol : recoCollisions) {
            maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
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

  template <typename Coll>
  void fillNPvContribHisto(const Coll& collision,
                           std::array<int, DataType::kDataTypes>& nCandsPerType,
                           std::array<int, DataType::kDataTypes>& nCandsInSignalRegionDsPerType,
                           std::array<int, DataType::kDataTypes>& nCandsInSignalRegionDplusPerType)
  {
    int numPvContributors = collision.numContrib();
    float centrality = evaluateCentralityColl(collision);
    std::get<TH2_ptr>(histosPtr[DataType::Data]["hNPvContribAll"])->Fill(numPvContributors, centrality);
    for (int i = 0; i < DataType::kDataTypes; i++) {
      if (nCandsPerType[i]) {
        std::get<TH2_ptr>(histosPtr[i]["hNPvContribCands"])->Fill(numPvContributors, centrality);
      }
      if (nCandsInSignalRegionDsPerType[i]) {
        std::get<TH2_ptr>(histosPtr[i]["hNPvContribCandsInSignalRegionDs"])->Fill(numPvContributors, centrality);
      }
      if (nCandsInSignalRegionDplusPerType[i]) {
        std::get<TH2_ptr>(histosPtr[i]["hNPvContribCandsInSignalRegionDplus"])->Fill(numPvContributors, centrality);
      }
    }
  }

  template <typename Coll, typename CandsDs>
  void runDataAnalysisPerCollision(const Coll& collisions, const CandsDs& candsDs, Preslice<CandsDs> candDsPerCollision)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      std::array<int, DataType::kDataTypes> nCandsPerType{0};
      std::array<int, DataType::kDataTypes> nCandsInSignalRegionDsPerType{0};
      std::array<int, DataType::kDataTypes> nCandsInSignalRegionDplusPerType{0};

      auto groupedDsCandidates = candsDs.sliceBy(candDsPerCollision, thisCollId);
      for (const auto& candidate : groupedDsCandidates) {
        if (candidate.isSelDsToKKPi() < selectionFlagDs && candidate.isSelDsToPiKK() < selectionFlagDs) {
          continue;
        }
        runDataAnalysisPerCandidate<Coll>(candidate);

        ++nCandsPerType[DataType::Data];
        if (isCandInSignalRegion(candidate, true)) {
          ++nCandsInSignalRegionDsPerType[DataType::Data];
        }
        if (isCandInSignalRegion(candidate, false)) {
          ++nCandsInSignalRegionDplusPerType[DataType::Data];
        }
      }
      fillNPvContribHisto(collision, nCandsPerType, nCandsInSignalRegionDsPerType, nCandsInSignalRegionDplusPerType);
    }
  }

  template <CentralityEstimator centDetector, typename Coll, typename CandsDs, typename CandDsMcGen>
  void runMcAnalysisPerCollision(const Coll& collisions,
                                 const CandsDs& candsDs,
                                 const CandDsMcGen& mcParticles,
                                 Preslice<CandsDs> candDsPerCollision)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      std::array<int, DataType::kDataTypes> nCandsPerType{0};
      std::array<int, DataType::kDataTypes> nCandsInSignalRegionDsPerType{0};
      std::array<int, DataType::kDataTypes> nCandsInSignalRegionDplusPerType{0};

      auto groupedDsCandidates = candsDs.sliceBy(candDsPerCollision, thisCollId);
      for (const auto& candidate : groupedDsCandidates) {
        if (candidate.isSelDsToKKPi() < selectionFlagDs && candidate.isSelDsToPiKK() < selectionFlagDs) {
          continue;
        }
        runDataAnalysisPerCandidate<Coll>(candidate);
        runMcAnalysisPerCandidate<Coll>(candidate, mcParticles);

        // Increase the number of candidates of the corresponding type to fill the NPvContrib histos
        std::array<MemberFunctionPointer<typename CandsDs::iterator>, 4> isOfType = {// Contains the functions to check if the candidate is of a certain type
                                                                                     &HfTaskDs::isDsPrompt<typename CandsDs::iterator>,
                                                                                     &HfTaskDs::isDsNonPrompt<typename CandsDs::iterator>,
                                                                                     &HfTaskDs::isDplusPrompt<typename CandsDs::iterator>,
                                                                                     &HfTaskDs::isDplusNonPrompt<typename CandsDs::iterator>};
        bool isBkg = true;
        for (int i = DataType::McDsPrompt; i <= DataType::McDplusNonPrompt; i++) { // Check what type of MC signal candidate it is, and fill the corresponding arrays
          if ((this->*isOfType[i - DataType::McDsPrompt])(candidate)) {
            isBkg = false;
            ++nCandsPerType[i];
            if (isCandInSignalRegion(candidate, true)) {
              ++nCandsInSignalRegionDsPerType[i];
            }
            if (isCandInSignalRegion(candidate, false)) {
              ++nCandsInSignalRegionDplusPerType[i];
            }
            break;
          }
        }
        if (isBkg) {
          ++nCandsPerType[DataType::McBkg];
          if (isCandInSignalRegion(candidate, true)) {
            ++nCandsInSignalRegionDsPerType[DataType::McBkg];
          }
          if (isCandInSignalRegion(candidate, false)) {
            ++nCandsInSignalRegionDplusPerType[DataType::McBkg];
          }
        }

        nCandsPerType[DataType::Data] = nCandsPerType[DataType::McDsPrompt] + nCandsPerType[DataType::McDsNonPrompt] + nCandsPerType[DataType::McDplusPrompt] + nCandsPerType[DataType::McDplusNonPrompt] + nCandsPerType[DataType::McBkg];

        nCandsInSignalRegionDsPerType[DataType::Data] = nCandsInSignalRegionDsPerType[DataType::McDsPrompt] + nCandsInSignalRegionDsPerType[DataType::McDsNonPrompt] + nCandsInSignalRegionDsPerType[DataType::McDplusPrompt] + nCandsInSignalRegionDsPerType[DataType::McDplusNonPrompt] + nCandsInSignalRegionDsPerType[DataType::McBkg];

        nCandsInSignalRegionDplusPerType[DataType::Data] = nCandsInSignalRegionDplusPerType[DataType::McDsPrompt] + nCandsInSignalRegionDplusPerType[DataType::McDsNonPrompt] + nCandsInSignalRegionDplusPerType[DataType::McDplusPrompt] + nCandsInSignalRegionDplusPerType[DataType::McDplusNonPrompt] + nCandsInSignalRegionDplusPerType[DataType::McBkg];
      }
      fillNPvContribHisto(collision, nCandsPerType, nCandsInSignalRegionDsPerType, nCandsInSignalRegionDplusPerType);
    }
    fillMcGenHistos<centDetector>(mcParticles, collisions);
  }

  void processDataWithCentFT0C(CollisionsWithFT0C const& collisions,
                               CandDsData const& candsDs,
                               aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs, candDsDataPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentFT0C, "Process data w/o ML information on Ds, with information on centrality from FT0C", false);

  void processDataWithCentFT0M(CollisionsWithFT0M const& collisions,
                               CandDsData const& candsDs,
                               aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs, candDsDataPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentFT0M, "Process data w/o ML information on Ds, with information on centrality from FT0M", false);

  void processDataWithCentNTracksPV(CollisionsWithNTracksPV const& collisions,
                                    CandDsData const& candsDs,
                                    aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs, candDsDataPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentNTracksPV, "Process data w/o ML information on Ds, with information on centrality from NTracksPV", false);

  void processData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   CandDsData const& candsDs,
                   aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs, candDsDataPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processData, "Process data w/o ML information on Ds, w/o information on centrality", true);

  void processDataWithMlAndCentFT0C(CollisionsWithFT0C const& collisions,
                                    CandDsDataWithMl const& candsDs,
                                    aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs, candDsDataWithMlPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentFT0C, "Process data with ML information on Ds, with information on centrality from FT0C", false);

  void processDataWithMlAndCentFT0M(CollisionsWithFT0M const& collisions,
                                    CandDsDataWithMl const& candsDs,
                                    aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs, candDsDataWithMlPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentFT0M, "Process data with ML information on Ds, with information on centrality from FT0M", false);

  void processDataWithMlAndCentNTracksPV(CollisionsWithNTracksPV const& collisions,
                                         CandDsDataWithMl const& candsDs,
                                         aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs, candDsDataWithMlPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentNTracksPV, "Process data with ML information on Ds, with information on centrality", false);

  void processDataWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                         CandDsDataWithMl const& candsDs,
                         aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs, candDsDataWithMlPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMl, "Process data with ML information on Ds, w/o information on centrality", false);

  void processMcWithCentFT0C(CollisionsMcWithFT0C const& collisions,
                             CandDsMcReco const& candsDs,
                             CandDsMcGen const& mcParticles,
                             aod::McCollisions const&,
                             aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CentralityEstimator::FT0C, CollisionsMcWithFT0C>(collisions, candsDs, mcParticles, candDsMcRecoPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentFT0C, "Process MC w/o ML information on Ds, with information on centrality from FT0C", false);

  void processMcWithCentFT0M(CollisionsMcWithFT0M const& collisions,
                             CandDsMcReco const& candsDs,
                             CandDsMcGen const& mcParticles,
                             aod::McCollisions const&,
                             aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CentralityEstimator::FT0M, CollisionsMcWithFT0M>(collisions, candsDs, mcParticles, candDsMcRecoPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentFT0M, "Process MC w/o ML information on Ds, with information on centrality from FT0M", false);

  void processMcWithCentNTracksPV(CollisionsMcWithNTracksPV const& collisions,
                                  CandDsMcReco const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CentralityEstimator::NTracksPV, CollisionsMcWithNTracksPV>(collisions, candsDs, mcParticles, candDsMcRecoPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentNTracksPV, "Process MC w/o ML information on Ds, with information on centrality from NTracksPV", false);

  void processMc(CollisionsMc const& collisions,
                 CandDsMcReco const& candsDs,
                 CandDsMcGen const& mcParticles,
                 aod::McCollisions const&,
                 aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CentralityEstimator::None, CollisionsMc>(collisions, candsDs, mcParticles, candDsMcRecoPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processMc, "Process MC w/o ML information on Ds, w/o information on centrality", false);

  void processMcWithMlAndCentFT0C(CollisionsMcWithFT0C const& collisions,
                                  CandDsMcRecoWithMl const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CentralityEstimator::FT0C, CollisionsMcWithFT0C>(collisions, candsDs, mcParticles, candDsMcRecoWithMlPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentFT0C, "Process MC with ML information on Ds, with information on centrality from FT0C", false);

  void processMcWithMlAndCentFT0M(CollisionsMcWithFT0M const& collisions,
                                  CandDsMcRecoWithMl const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CentralityEstimator::FT0M, CollisionsMcWithFT0M>(collisions, candsDs, mcParticles, candDsMcRecoWithMlPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentFT0M, "Process MC with ML information on Ds, with information on centrality from FT0M", false);

  void processMcWithMlAndCentNTracksPV(CollisionsMcWithNTracksPV const& collisions,
                                       CandDsMcRecoWithMl const& candsDs,
                                       CandDsMcGen const& mcParticles,
                                       aod::McCollisions const&,
                                       aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CentralityEstimator::NTracksPV, CollisionsMcWithNTracksPV>(collisions, candsDs, mcParticles, candDsMcRecoWithMlPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentNTracksPV, "Process MC with ML information on Ds, with information on centrality from NTracksPV", false);

  void processMcWithMl(CollisionsMc const& collisions,
                       CandDsMcRecoWithMl const& candsDs,
                       CandDsMcGen const& mcParticles,
                       aod::McCollisions const&,
                       aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CentralityEstimator::None, CollisionsMc>(collisions, candsDs, mcParticles, candDsMcRecoWithMlPerCollision);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMl, "Process MC with ML information on Ds, w/o information on centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDs>(cfgc)};
}
