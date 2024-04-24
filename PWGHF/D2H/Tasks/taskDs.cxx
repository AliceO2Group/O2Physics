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
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f}, "axis for pT"};
  ConfigurableAxis axisMlScore0{"axisMlScore0", {100, 0., 1.}, "axis for ML output score 0"};
  ConfigurableAxis axisMlScore1{"axisMlScore1", {100, 0., 1.}, "axis for ML output score 1"};
  ConfigurableAxis axisMlScore2{"axisMlScore2", {100, 0., 1.}, "axis for ML output score 2"};

  HfHelper hfHelper;

  using CentralityEstimator = o2::aod::hf_collision_centrality::CentralityEstimator;
  using TH1_ptr = std::shared_ptr<TH1>;
  using TH2_ptr = std::shared_ptr<TH2>;
  using THnSparse_ptr = std::shared_ptr<THnSparse>;
  using histTypes = std::variant<TH1_ptr, TH2_ptr, THnSparse_ptr>;

  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsWithNTracksPV = soa::Join<aod::Collisions, aod::EvSels, aod::CentNTPVs>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsDataWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcRecoWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec, aod::HfMlDsToKKPi>>;
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

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

    for (auto i = 0; i < DataType::kDataTypes; ++i) {
      if (doprocessDataWithCentFT0C || doprocessDataWithCentFT0M || doprocessDataWithCentNTracksPV ||
          doprocessMcWithCentFT0C || doprocessMcWithCentFT0M || doprocessMcWithCentNTracksPV) {
        histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins});
      } else if (doprocessDataWithMlAndCentFT0C || doprocessDataWithMlAndCentFT0M || doprocessDataWithMlAndCentNTracksPV ||
                 doprocessMcWithMlAndCentFT0C || doprocessMcWithMlAndCentFT0M || doprocessMcWithMlAndCentNTracksPV) {
        histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, centralitybins, axisMlScore0, axisMlScore1, axisMlScore2});
      } else if (doprocessData || doprocessMc) {
        histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins});
      } else if (doprocessDataWithMl || doprocessMcWithMl) {
        histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, {massbins, ptbins, axisMlScore0, axisMlScore1, axisMlScore2});
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
          histosPtr[i]["hPtVsYGen"] = registry.add<TH2>((folders[i] + "hPtVsYGen").c_str(), "MC particles (unmatched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
        }
      }
    }
  }

  /// Evaluate multiplicity
  /// \param candidate is candidate
  /// \return multiplicity of the collision associated to the candidate
  template <CentralityEstimator centDetector, typename T1>
  int evaluateCentrality(const T1& candidate)
  {
    if constexpr (centDetector == CentralityEstimator::FT0C)
      return candidate.template collision_as<CollisionsWithFT0C>().centFT0C();
    else if constexpr (centDetector == CentralityEstimator::FT0M)
      return candidate.template collision_as<CollisionsWithFT0M>().centFT0M();
    else if constexpr (centDetector == CentralityEstimator::NTracksPV)
      return candidate.template collision_as<CollisionsWithNTracksPV>().centNTPV();
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

  /// Fill histograms of quantities for the KKPi daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <CentralityEstimator centDetector, bool useMl, typename T1>
  void fillHistoKKPi(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    if constexpr (useMl) {
      std::vector<float> outputMl = {-999., -999., -999.};
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) { // TODO: add checks for classMl size
        outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
      }
      if constexpr (centDetector != CentralityEstimator::None) {
        std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, evaluateCentrality<centDetector>(candidate), outputMl[0], outputMl[1], outputMl[2]);
      } else {
        std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, outputMl[0], outputMl[1], outputMl[2]);
      }
    } else {
      if constexpr (centDetector != CentralityEstimator::None) {
        std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt, evaluateCentrality<centDetector>(candidate));
      } else {
        std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToKKPi(candidate), pt);
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
  template <CentralityEstimator centDetector, bool useMl, typename T1>
  void fillHistoPiKK(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    if constexpr (useMl) {
      std::vector<float> outputMl = {-999., -999., -999.};
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) { // TODO: add checks for classMl size
        outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
      }
      if constexpr (centDetector != CentralityEstimator::None) {
        std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, evaluateCentrality<centDetector>(candidate), outputMl[0], outputMl[1], outputMl[2]);
      } else {
        std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, outputMl[0], outputMl[1], outputMl[2]);
      }
    } else {
      if constexpr (centDetector != CentralityEstimator::None) {
        std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt, evaluateCentrality<centDetector>(candidate));
      } else {
        std::get<THnSparse_ptr>(histosPtr[dataType]["hSparseMass"])->Fill(hfHelper.invMassDsToPiKK(candidate), pt);
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
  template <CentralityEstimator centDetector, bool useMl, typename T1>
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

        auto y = hfHelper.yDs(candidate);

        // prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          fillHisto(candidate, DataType::McDsPrompt);
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) // KKPi
            fillHistoKKPi<centDetector, useMl>(candidate, DataType::McDsPrompt);
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) // PiKK
            fillHistoPiKK<centDetector, useMl>(candidate, DataType::McDsPrompt);

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
          fillHisto(candidate, DataType::McDsNonPrompt);
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) // KKPi
            fillHistoKKPi<centDetector, useMl>(candidate, DataType::McDsNonPrompt);
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) // PiKK
            fillHistoPiKK<centDetector, useMl>(candidate, DataType::McDsNonPrompt);

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
        auto y = hfHelper.yDplus(candidate);

        // prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          fillHisto(candidate, DataType::McDplusPrompt);
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) // KKPi
            fillHistoKKPi<centDetector, useMl>(candidate, DataType::McDplusPrompt);
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) // PiKK
            fillHistoPiKK<centDetector, useMl>(candidate, DataType::McDplusPrompt);

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
          fillHisto(candidate, DataType::McDplusNonPrompt);
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) // KKPi
            fillHistoKKPi<centDetector, useMl>(candidate, DataType::McDplusNonPrompt);
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) // PiKK
            fillHistoPiKK<centDetector, useMl>(candidate, DataType::McDplusNonPrompt);

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
        auto y = hfHelper.yDplus(candidate);

        // Fill whether it is prompt or non-prompt
        fillHisto(candidate, DataType::McDplusBkg);

        if (candidate.isSelDsToKKPi() >= selectionFlagDs) // KKPi
          fillHistoKKPi<centDetector, useMl>(candidate, DataType::McDplusBkg);
        if (candidate.isSelDsToPiKK() >= selectionFlagDs) // PiKK
          fillHistoPiKK<centDetector, useMl>(candidate, DataType::McDplusBkg);

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

  template <FinalState decayChannel, CentralityEstimator centDetector, bool useMl, typename CandsDs>
  void runDataAnalysis(CandsDs const& candidates)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillHisto(candidate, DataType::Data);
      if constexpr (decayChannel == FinalState::KKPi) { // KKPi
        fillHistoKKPi<centDetector, useMl>(candidate, DataType::Data);
      } else if constexpr (decayChannel == FinalState::PiKK) { // PiKK
        fillHistoPiKK<centDetector, useMl>(candidate, DataType::Data);
      }
    }
  }

  template <CentralityEstimator centDetector, bool useMl, typename CandsDs>
  void runMcAnalysis(CandsDs const& /*candidates*/,
                     CandDsMcGen const& mcParticles)
  {
    // MC rec.
    // Ds± → K± K∓ π±
    if constexpr (useMl) {
      // Ds
      for (const auto& candidate : reconstructedCandDsSigWithMl)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl>(candidate, mcParticles, SpeciesAndDecay::DsToKKPi);
      // D+→ K± K∓ π±
      for (const auto& candidate : reconstructedCandDplusSigWithMl)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl>(candidate, mcParticles, SpeciesAndDecay::DplusToKKPi);
      // D+→ π± K∓ π±
      for (const auto& candidate : reconstructedCandDplusBkgWithMl)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl>(candidate, mcParticles, SpeciesAndDecay::DplusToPiKPi);
      // Bkg
      for (const auto& candidate : reconstructedCandBkgWithMl) {
        if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
          continue;
        }

        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHisto(candidate, DataType::McBkg);
      }
    } else {
      // Ds
      for (const auto& candidate : reconstructedCandDsSig)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl>(candidate, mcParticles, SpeciesAndDecay::DsToKKPi);

      // D+→ K± K∓ π±
      for (const auto& candidate : reconstructedCandDplusSig)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl>(candidate, mcParticles, SpeciesAndDecay::DplusToKKPi);

      // D+→ π± K∓ π±
      for (const auto& candidate : reconstructedCandDplusBkg)
        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs)
          fillHistoMCRec<centDetector, useMl>(candidate, mcParticles, SpeciesAndDecay::DplusToPiKPi);

      // Bkg
      for (const auto& candidate : reconstructedCandBkg) {
        if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
          continue;
        }

        if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs) {
          fillHisto(candidate, DataType::McBkg);
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) // KKPi
            fillHistoKKPi<centDetector, useMl>(candidate, DataType::McDsPrompt);
          if (candidate.isSelDsToPiKK() >= selectionFlagDs) // PiKK
            fillHistoPiKK<centDetector, useMl>(candidate, DataType::McDsPrompt);
        }
      }
    }

    // TODO: add histograms for reflections

    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        if (particle.flagMcDecayChanGen() == decayChannel || (fillDplusMc && particle.flagMcDecayChanGen() == (decayChannel + offsetDplusDecayChannel))) {
          auto pt = particle.pt();
          auto y = 0;

          if (particle.flagMcDecayChanGen() == decayChannel) {
            y = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDS);
            if (yCandGenMax >= 0. && std::abs(y) > yCandGenMax) {
              continue;
            }
            if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
              std::get<TH1_ptr>(histosPtr[DataType::McDsPrompt]["hPtGen"])->Fill(pt);              // gen. level pT
              std::get<TH2_ptr>(histosPtr[DataType::McDsPrompt]["hPtVsYGen"])->Fill(pt, y);        // gen. level pT
              std::get<TH1_ptr>(histosPtr[DataType::McDsPrompt]["hEtaGen"])->Fill(particle.eta()); // gen. level pT
            }
            if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              std::get<TH1_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtGen"])->Fill(pt);              // gen. level pT
              std::get<TH2_ptr>(histosPtr[DataType::McDsNonPrompt]["hPtVsYGen"])->Fill(pt, y);        // gen. level pT
              std::get<TH1_ptr>(histosPtr[DataType::McDsNonPrompt]["hEtaGen"])->Fill(particle.eta()); // gen. level pT
            }
          } else if (fillDplusMc) {
            y = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus);
            if (yCandGenMax >= 0. && std::abs(y) > yCandGenMax) {
              continue;
            }
            if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
              std::get<TH1_ptr>(histosPtr[DataType::McDplusPrompt]["hPtGen"])->Fill(pt);       // gen. level pT
              std::get<TH2_ptr>(histosPtr[DataType::McDplusPrompt]["hPtVsYGen"])->Fill(pt, y); // gen. level pT
              std::get<TH1_ptr>(histosPtr[DataType::McDplusPrompt]["hEtaGen"])->Fill(particle.eta());
            }
            if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              std::get<TH1_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtGen"])->Fill(pt); // gen. level pT
              std::get<TH2_ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtVsYGen"])->Fill(pt, y);
              std::get<TH1_ptr>(histosPtr[DataType::McDplusNonPrompt]["hEtaGen"])->Fill(particle.eta());
            }
          }
        }
      }
    }
  }

  void processDataWithCentFT0C(CollisionsWithFT0C const&,
                               CandDsData const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0C, false>(selectedDsToKKPiCandData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0C, false>(selectedDsToPiKKCandData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentFT0C, "Process data w/o ML information on Ds, with information on centrality from FT0C", false);

  void processDataWithCentFT0M(CollisionsWithFT0M const&,
                               CandDsData const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0M, false>(selectedDsToKKPiCandData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0M, false>(selectedDsToPiKKCandData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentFT0M, "Process data w/o ML information on Ds, with information on centrality from FT0M", false);

  void processDataWithCentNTracksPV(CollisionsWithNTracksPV const&,
                                    CandDsData const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::NTracksPV, false>(selectedDsToKKPiCandData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::NTracksPV, false>(selectedDsToPiKKCandData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentNTracksPV, "Process data w/o ML information on Ds, with information on centrality from NTracksPV", false);

  void processData(aod::Collisions const&,
                   CandDsData const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::None, false>(selectedDsToKKPiCandData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::None, false>(selectedDsToPiKKCandData);
  }
  PROCESS_SWITCH(HfTaskDs, processData, "Process data w/o ML information on Ds, w/o information on centrality", true);

  void processDataWithMlAndCentFT0C(CollisionsWithFT0C const&,
                                    CandDsDataWithMl const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0C, true>(selectedDsToKKPiCandWithMlData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0C, true>(selectedDsToPiKKCandWithMlData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentFT0C, "Process data with ML information on Ds, with information on centrality from FT0C", false);

  void processDataWithMlAndCentFT0M(CollisionsWithFT0M const&,
                                    CandDsDataWithMl const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0M, true>(selectedDsToKKPiCandWithMlData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0M, true>(selectedDsToPiKKCandWithMlData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentFT0M, "Process data with ML information on Ds, with information on centrality from FT0M", false);

  void processDataWithMlAndCentNTracksPV(CollisionsWithNTracksPV const&,
                                         CandDsDataWithMl const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::NTracksPV, true>(selectedDsToKKPiCandWithMlData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::NTracksPV, true>(selectedDsToPiKKCandWithMlData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentNTracksPV, "Process data with ML information on Ds, with information on centrality", false);

  void processDataWithMl(aod::Collisions const&,
                         CandDsDataWithMl const&)
  {
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::None, true>(selectedDsToKKPiCandWithMlData);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::None, true>(selectedDsToPiKKCandWithMlData);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMl, "Process data with ML information on Ds, w/o information on centrality", false);

  void processMcWithCentFT0C(CollisionsWithFT0C const&,
                             CandDsMcReco const& candidates,
                             CandDsMcGen const& mcParticles,
                             aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::FT0C, false>(candidates, mcParticles);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0C, false>(selectedDsToKKPiCandMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0C, false>(selectedDsToPiKKCandMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentFT0C, "Process MC w/o ML information on Ds, with information on centrality from FT0C", false);

  void processMcWithCentFT0M(CollisionsWithFT0M const&,
                             CandDsMcReco const& candidates,
                             CandDsMcGen const& mcParticles,
                             aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::FT0M, false>(candidates, mcParticles);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0M, false>(selectedDsToKKPiCandMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0M, false>(selectedDsToPiKKCandMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentFT0M, "Process MC w/o ML information on Ds, with information on centrality from FT0M", false);

  void processMcWithCentNTracksPV(CollisionsWithNTracksPV const&,
                                  CandDsMcReco const& candidates,
                                  CandDsMcGen const& mcParticles,
                                  aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::NTracksPV, false>(candidates, mcParticles);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::NTracksPV, false>(selectedDsToKKPiCandMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::NTracksPV, false>(selectedDsToPiKKCandMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentNTracksPV, "Process MC w/o ML information on Ds, with information on centrality from NTracksPV", false);

  void processMc(aod::Collisions const&,
                 CandDsMcReco const& candidates,
                 CandDsMcGen const& mcParticles,
                 aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::None, false>(candidates, mcParticles);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::None, false>(selectedDsToKKPiCandMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::None, false>(selectedDsToPiKKCandMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMc, "Process MC w/o ML information on Ds, w/o information on centrality", false);

  void processMcWithMlAndCentFT0C(CollisionsWithFT0C const&,
                                  CandDsMcRecoWithMl const& candidates,
                                  CandDsMcGen const& mcParticles,
                                  aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::FT0C, true>(candidates, mcParticles);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0C, true>(selectedDsToKKPiCandWithMlMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0C, true>(selectedDsToPiKKCandWithMlMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentFT0C, "Process MC with ML information on Ds, with information on centrality from FT0C", false);

  void processMcWithMlAndCentFT0M(CollisionsWithFT0M const&,
                                  CandDsMcRecoWithMl const& candidates,
                                  CandDsMcGen const& mcParticles,
                                  aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::FT0M, true>(candidates, mcParticles);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::FT0M, true>(selectedDsToKKPiCandWithMlMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::FT0M, true>(selectedDsToPiKKCandWithMlMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentFT0M, "Process MC with ML information on Ds, with information on centrality from FT0M", false);

  void processMcWithMlAndCentNTracksPV(CollisionsWithNTracksPV const&,
                                       CandDsMcRecoWithMl const& candidates,
                                       CandDsMcGen const& mcParticles,
                                       aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::NTracksPV, true>(candidates, mcParticles);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::NTracksPV, true>(selectedDsToKKPiCandWithMlMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::NTracksPV, true>(selectedDsToPiKKCandWithMlMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentNTracksPV, "Process MC with ML information on Ds, with information on centrality from NTracksPV", false);

  void processMcWithMl(aod::Collisions const&,
                       CandDsMcRecoWithMl const& candidates,
                       CandDsMcGen const& mcParticles,
                       aod::TracksWMc const&)
  {
    runMcAnalysis<CentralityEstimator::None, true>(candidates, mcParticles);
    runDataAnalysis<FinalState::KKPi, CentralityEstimator::None, true>(selectedDsToKKPiCandWithMlMc);
    runDataAnalysis<FinalState::PiKK, CentralityEstimator::None, true>(selectedDsToPiKKCandWithMlMc);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMl, "Process MC with ML information on Ds, w/o information on centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDs>(cfgc)};
}
