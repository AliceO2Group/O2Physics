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

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/MetadataHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
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
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TProfile.h>

#include <Rtypes.h>

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

enum FinalState { KKPi = 0,
                  PiKK };

enum DataType { Data = 0,
                McDsPrompt,
                McDsNonPrompt,
                McDplusPrompt,
                McDplusNonPrompt,
                McDplusBkg,
                McLcBkg,
                McBkg,
                kDataTypes };

enum Mother : int8_t {
  Ds,
  Dplus
};

enum ResonantChannel : int8_t {
  PhiPi = 1,
  Kstar0K = 2
};

namespace o2::aod
{
namespace hf_cand_ds_mini
{
DECLARE_SOA_COLUMN(Pt, pt, float);                                           //! Transverse momentum of D-meson candidate (GeV/c)
DECLARE_SOA_COLUMN(M, m, float);                                             //! Invariant mass of D-meson candidates (GeV/c)
DECLARE_SOA_COLUMN(Centrality, centrality, float);                           //! Centrality of collision
DECLARE_SOA_COLUMN(ImpactParameter, impactParameter, float);                 //! Impact parameter of D-meson candidate
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                         //! Decay length of D-meson candidate
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                     //! Transverse decay length of D-meson candidate
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);     //! Normalised decay length of D-meson candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float); //! Normalised transverse decay length of D-meson candidate
DECLARE_SOA_COLUMN(FlagMc, flagMc, int);                                     //! MC flag (according to DataType enum)
} // namespace hf_cand_ds_mini

DECLARE_SOA_TABLE(HfCandDsMinis, "AOD", "HFCANDDSMINI", //! Table with few Ds properties
                  hf_cand_ds_mini::M,
                  hf_cand_ds_mini::Pt,
                  hf_cand_ds_mini::Centrality);

DECLARE_SOA_TABLE(HfCandDsDlMinis, "AOD", "HFCANDDSDLMINI", //! Table with decay length Ds properties
                  hf_cand_ds_mini::DecayLength,
                  hf_cand_ds_mini::DecayLengthXY,
                  hf_cand_ds_mini::DecayLengthNormalised,
                  hf_cand_ds_mini::DecayLengthXYNormalised);

DECLARE_SOA_TABLE(HfCandDsD0Minis, "AOD", "HFCANDDSD0MINI", //! Table with impact parameter (d0)
                  hf_cand_ds_mini::ImpactParameter);

DECLARE_SOA_TABLE(HfCandDsMcMinis, "AOD", "HFCANDDSMCMINI", //! Table with MC decay type check
                  hf_cand_ds_mini::FlagMc);
} // namespace o2::aod

static std::unordered_map<int8_t, std::unordered_map<int8_t, int8_t>> channelsResonant = {{{Mother::Ds, {{ResonantChannel::PhiPi, hf_decay::hf_cand_3prong::DecayChannelResonant::DsToPhiPi}, {ResonantChannel::Kstar0K, hf_decay::hf_cand_3prong::DecayChannelResonant::DsToKstar0K}}},
                                                                                           {Mother::Dplus, {{ResonantChannel::PhiPi, hf_decay::hf_cand_3prong::DecayChannelResonant::DplusToPhiPi}, {ResonantChannel::Kstar0K, hf_decay::hf_cand_3prong::DecayChannelResonant::DplusToKstar0K}}}}};

template <typename T>
concept HasDsMlInfo = requires(T candidate) {
  candidate.mlProbDsToKKPi();
  candidate.mlProbDsToPiKK();
};

/// Ds± analysis task
struct HfTaskDs {
  Produces<aod::HfCandDsMinis> hfCandDsMinis;
  Produces<aod::HfCandDsDlMinis> hfCandDsDlMinis;
  Produces<aod::HfCandDsD0Minis> hfCandDsD0Minis;
  Produces<aod::HfCandDsMcMinis> hfCandDsMcMinis;

  Configurable<int> decayChannel{"decayChannel", 1, "Switch between resonant decay channels: 1 for Ds/Dplus->PhiPi->KKpi, 2 for Ds/Dplus->K0*K->KKPi"};
  Configurable<bool> fillDplusMc{"fillDplusMc", true, "Switch to fill Dplus MC information"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2, 3}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> massDsSignalMin{"massDsSignalMin", 1.934, "min mass for Ds signal"};
  Configurable<float> massDsSignalMax{"massDsSignalMax", 1.994, "max mass for Ds signal"};
  Configurable<float> massDplusSignalMin{"massDplusSignalMin", 1.866, "min mass for Dplus signal"};
  Configurable<float> massDplusSignalMax{"massDplusSignalMax", 1.906, "max mass for Dplus signal"};
  Configurable<bool> fillPercentiles{"fillPercentiles", true, "Wheter to fill multiplicity axis with percentiles or raw information"};
  Configurable<bool> storeOccupancy{"storeOccupancy", false, "Flag to store occupancy information"};
  Configurable<int> occEstimator{"occEstimator", 0, "Occupancy estimation (None: 0, ITS: 1, FT0C: 2)"};
  Configurable<bool> fillMcBkgHistos{"fillMcBkgHistos", false, "Flag to fill and store histograms for MC background"};
  struct : ConfigurableGroup {
    Configurable<bool> produceMiniTrees{"produceMiniTrees", false, "Flag to produce mini trees"};
    Configurable<bool> extendWithDecayLength{"extendWithDecayLength", false, "Flag to extend trees with decay length information"};
    Configurable<bool> extendWithImpactParameter{"extendWithImpactParameter", false, "Flag to extend trees with impact parameter information"};
  } miniTrees;

  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "The CCDB endpoint url address"};
    Configurable<std::string> ccdbPath{"ccdbPath", "Centrality/Calibration", "The CCDB path for centrality/multiplicity information"};
    Configurable<std::string> reconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};
  } ccdbConfig;

  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  using TH1Ptr = std::shared_ptr<TH1>;
  using TH2Ptr = std::shared_ptr<TH2>;
  using THnSparsePtr = std::shared_ptr<THnSparse>;
  using HistTypes = std::variant<TH1Ptr, TH2Ptr, THnSparsePtr>;
  template <typename CandDs>
  using MemberFunctionPointer = bool (HfTaskDs::*)(const CandDs&);

  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::CentFT0Ms>;
  using CollisionsWithNTracksPV = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentNTPVs>;

  using CollisionsMc = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using CollisionsMcWithFT0C = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::FT0Mults, aod::CentFT0Cs>;
  using CollisionsMcWithFT0M = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::FT0Mults, aod::CentFT0Ms>;
  using CollisionsMcWithNTracksPV = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::PVMults, aod::CentNTPVs>;

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsDataWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcRecoWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec, aod::HfMlDsToKKPi>>;
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Filter filterDsFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi))) != static_cast<uint8_t>(0);

  Preslice<aod::HfCand3Prong> candDsPerCollision = aod::hf_cand::collisionId;
  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f}, "axis for pT"};
  ConfigurableAxis axisPtBHad{"axisPtBHad", {50, 0., 100}, "axis for pt of B hadron decayed into D candidate"};
  ConfigurableAxis axisFlagBHad{"axisFlagBHad", {5, 0, 5}, "axis for B hadron mother flag"};
  ConfigurableAxis axisNPvContributors{"axisNPvContributors", {200, -0.5f, 199.5f}, "axis for NPvContributors"};
  ConfigurableAxis axisMlScore0{"axisMlScore0", {100, 0., 1.}, "axis for ML output score 0"};
  ConfigurableAxis axisMlScore1{"axisMlScore1", {100, 0., 1.}, "axis for ML output score 1"};
  ConfigurableAxis axisMlScore2{"axisMlScore2", {100, 0., 1.}, "axis for ML output score 2"};
  ConfigurableAxis axisCentrality{"axisCentrality", {100, 0, 100}, "axis for centrality/multiplicity"};
  ConfigurableAxis axisOccupancy{"axisOccupancy", {14, 0., 14000.}, "axis for occupancy"};

  int mRunNumber{0};
  bool lCalibLoaded{};
  TList* lCalibObjects{};
  TProfile* hVtxZFT0A{};
  TProfile* hVtxZFT0C{};
  TProfile* hVtxZNTracks{};

  HistogramRegistry registry{"registry", {}};

  std::array<std::string, DataType::kDataTypes> folders = {"Data/", "MC/Ds/Prompt/", "MC/Ds/NonPrompt/", "MC/Dplus/Prompt/", "MC/Dplus/NonPrompt/", "MC/Dplus/Bkg/", "MC/Lc/", "MC/Bkg/"};

  std::unordered_map<std::string, HistTypes> dataHistograms;
  std::unordered_map<std::string, HistTypes> mcDsPromptHistograms;
  std::unordered_map<std::string, HistTypes> mcDsNonPromptHistograms;
  std::unordered_map<std::string, HistTypes> mcDplusPromptHistograms;
  std::unordered_map<std::string, HistTypes> mcDplusNonPromptHistograms;
  std::unordered_map<std::string, HistTypes> mcDplusBkgHistograms;
  std::unordered_map<std::string, HistTypes> mcLcBkgHistograms;
  std::unordered_map<std::string, HistTypes> mcBkgHistograms;

  std::array<std::unordered_map<std::string, HistTypes>, DataType::kDataTypes> histosPtr = {dataHistograms, mcDsPromptHistograms, mcDsNonPromptHistograms, mcDplusPromptHistograms, mcDplusNonPromptHistograms, mcDplusBkgHistograms, mcLcBkgHistograms, mcBkgHistograms};

  void init(InitContext&)
  {
    std::array<int, 16> processes = {doprocessDataWithCentFT0C, doprocessDataWithCentFT0M, doprocessDataWithCentNTracksPV, doprocessData, doprocessDataWithMlAndCentFT0C, doprocessDataWithMlAndCentFT0M, doprocessDataWithMlAndCentNTracksPV, doprocessDataWithMl, doprocessMcWithCentFT0C, doprocessMcWithCentFT0M, doprocessMcWithCentNTracksPV, doprocessMc, doprocessMcWithMlAndCentFT0C, doprocessMcWithMlAndCentFT0M, doprocessMcWithMlAndCentNTracksPV, doprocessMcWithMl};

    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    if (decayChannel != ResonantChannel::PhiPi && decayChannel != ResonantChannel::Kstar0K) {
      LOGP(fatal, "Invalid value of decayChannel");
    }

    AxisSpec const ptbins{axisPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const ptBHad{axisPtBHad, "#it{p}_{T}(B) (GeV/#it{c})"};
    AxisSpec const flagBHad{axisFlagBHad, "B Hadron flag"};
    AxisSpec ybins = {100, -5., 5, "#it{y}"};
    AxisSpec const massbins = {600, 1.67, 2.27, "inv. mass (KK#pi) (GeV/#it{c}^{2})"};
    AxisSpec const centralitybins = {axisCentrality, "Centrality"};
    AxisSpec const npvcontributorsbins = {axisNPvContributors, "NPvContributors"};
    AxisSpec const mlscore0bins = {axisMlScore0, "Score 0"};
    AxisSpec const mlscore1bins = {axisMlScore1, "Score 1"};
    AxisSpec const mlscore2bins = {axisMlScore2, "Score 2"};
    AxisSpec const occupancybins = {axisOccupancy, "Occupancy"};

    histosPtr[DataType::Data]["hNPvContribAll"] = registry.add<TH2>((folders[DataType::Data] + "hNPvContribAll").c_str(), "3-prong candidates;NPvContributors;Centrality;Entries", HistType::kTH2F, {axisNPvContributors, {100, 0., 100}});

    std::vector<AxisSpec> axes = {massbins, ptbins, centralitybins};
    std::vector<AxisSpec> axesMl = {massbins, ptbins, centralitybins, mlscore0bins, mlscore1bins, mlscore2bins};
    std::vector<AxisSpec> axesFdWithNpv = {massbins, ptbins, centralitybins, npvcontributorsbins, ptBHad, flagBHad};
    std::vector<AxisSpec> axesFdWithNpvMl = {massbins, ptbins, centralitybins, mlscore0bins, mlscore1bins, mlscore2bins, npvcontributorsbins, ptBHad, flagBHad};
    std::vector<AxisSpec> axesWithNpv = {massbins, ptbins, centralitybins, npvcontributorsbins};
    std::vector<AxisSpec> axesWithNpvMl = {massbins, ptbins, centralitybins, mlscore0bins, mlscore1bins, mlscore2bins, npvcontributorsbins};
    std::vector<AxisSpec> axesGenPrompt = {ptbins, ybins, npvcontributorsbins, centralitybins};
    std::vector<AxisSpec> axesGenFd = {ptbins, ybins, npvcontributorsbins, centralitybins, ptBHad, flagBHad};

    if (storeOccupancy) {
      axes.insert(axes.end(), {occupancybins});
      axesMl.insert(axesMl.end(), {occupancybins});
      axesFdWithNpv.insert(axesFdWithNpv.end(), {occupancybins});
      axesFdWithNpvMl.insert(axesFdWithNpvMl.end(), {occupancybins});
      axesWithNpv.insert(axesWithNpv.end(), {occupancybins});
      axesWithNpvMl.insert(axesWithNpvMl.end(), {occupancybins});
      axesGenPrompt.insert(axesGenPrompt.end(), {occupancybins});
      axesGenFd.insert(axesGenFd.end(), {occupancybins});
    }

    for (auto i = 0; i < DataType::kDataTypes; ++i) {
      if (doprocessDataWithCentFT0C || doprocessDataWithCentFT0M || doprocessDataWithCentNTracksPV || doprocessData || doprocessMcWithCentFT0C || doprocessMcWithCentFT0M || doprocessMcWithCentNTracksPV || doprocessMc) {
        if (i == DataType::Data) { // If data do not fill PV contributors in sparse
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, axes);
        } else if (i == DataType::McDsNonPrompt || i == DataType::McDplusNonPrompt) {
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, axesFdWithNpv);
        } else {
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, axesWithNpv);
        }
      } else if (doprocessDataWithMlAndCentFT0C || doprocessDataWithMlAndCentFT0M || doprocessDataWithMlAndCentNTracksPV || doprocessDataWithMl || doprocessMcWithMlAndCentFT0C || doprocessMcWithMlAndCentFT0M || doprocessMcWithMlAndCentNTracksPV || doprocessMcWithMl) {
        if (i == DataType::McBkg && !fillMcBkgHistos) {
          continue;
        }

        if (i == DataType::Data) { // If data do not fill PV contributors in sparse
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, axesMl);
        } else if (i == DataType::McDsNonPrompt || i == DataType::McDplusNonPrompt) {
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, axesFdWithNpvMl);
        } else {
          histosPtr[i]["hSparseMass"] = registry.add<THnSparse>((folders[i] + "hSparseMass").c_str(), "THn for Ds", HistType::kTHnSparseF, axesWithNpvMl);
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
        if (i == DataType::McDsPrompt || i == DataType::McDsNonPrompt || i == DataType::McDplusPrompt || i == DataType::McDplusNonPrompt || i == DataType::McDplusBkg || i == DataType::McLcBkg) {
          histosPtr[i]["hEtaGen"] = registry.add<TH1>((folders[i] + "hEtaGen").c_str(), "3-prong candidates (matched);#eta;entries", {HistType::kTH1F, {{100, -2., 2.}}});
          histosPtr[i]["hPtGen"] = registry.add<TH1>((folders[i] + "hPtGen").c_str(), "MC particles (unmatched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {ptbins}});
          histosPtr[i]["hPtVsYRecoPID"] = registry.add<TH2>((folders[i] + "hPtVsYRecoPID").c_str(), "3-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
          histosPtr[i]["hPtVsYRecoTopol"] = registry.add<TH2>((folders[i] + "hPtVsYRecoTopol").c_str(), "3-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
          histosPtr[i]["hPtVsYRecoSkim"] = registry.add<TH2>((folders[i] + "hPtVsYRecoSkim").c_str(), "3-prong candidates (RecoSkim - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {ptbins, {ybins}}});
        }
        if (i == DataType::McDsPrompt || i == DataType::McDplusPrompt) {
          histosPtr[i]["hSparseGen"] = registry.add<THnSparse>((folders[i] + "hSparseGen").c_str(), "Thn for generated prompt candidates", HistType::kTHnSparseF, axesGenPrompt);
        }
        if (i == DataType::McDsNonPrompt || i == DataType::McDplusNonPrompt) {
          histosPtr[i]["hSparseGen"] = registry.add<THnSparse>((folders[i] + "hSparseGen").c_str(), "Thn for generated nonprompt candidates", HistType::kTHnSparseF, axesGenFd);
        }
      }
    }
  }

  template <typename CandDs>
  bool isDsPrompt(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) && candidate.flagMcDecayChanRec() == channelsResonant[Mother::Ds][decayChannel] && candidate.originMcRec() == RecoDecay::OriginType::Prompt;
  }

  template <typename CandDs>
  bool isDplusPrompt(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKK) && candidate.flagMcDecayChanRec() == channelsResonant[Mother::Dplus][decayChannel] && candidate.originMcRec() == RecoDecay::OriginType::Prompt;
  }

  template <typename CandDs>
  bool isDsNonPrompt(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) && candidate.flagMcDecayChanRec() == channelsResonant[Mother::Ds][decayChannel] && candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
  }

  template <typename CandDs>
  bool isDplusNonPrompt(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKK) && candidate.flagMcDecayChanRec() == channelsResonant[Mother::Dplus][decayChannel] && candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
  }

  template <typename CandDs>
  bool isDplusBkg(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi);
  }

  template <typename CandDs>
  bool isLcBkg(const CandDs& candidate)
  {
    return std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi);
  }

  /// Checks whether the candidate is in the signal region of either the Ds or D+ decay
  /// \param candidate is the candidate
  /// \param isDs is true if we check for the Ds signal region, false for the D+ signal region
  /// \return true if the candidate is in the signal region, false otherwise
  template <typename CandDs>
  bool isCandInSignalRegion(const CandDs& candidate, bool isDs)
  {
    bool const isKKPi = candidate.isSelDsToKKPi() >= selectionFlagDs;
    float const invMass = isKKPi ? HfHelper::invMassDsToKKPi(candidate) : HfHelper::invMassDsToPiKK(candidate);
    if (isDs && (invMass < massDsSignalMin || invMass > massDsSignalMax)) {
      return false;
    }
    if (!isDs && (invMass < massDplusSignalMin || invMass > massDplusSignalMax)) {
      return false;
    }
    return true;
  }

  /// Evaluate centrality/multiplicity percentile using FT0M estimator
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <o2::hf_centrality::HasFT0MCent Coll>
  float getZEqMultColl(const Coll& collision, uint8_t nProngsContributorsPV)
  {
    auto multFT0A = collision.multFT0A() - nProngsContributorsPV;
    auto multFT0C = collision.multFT0C() - nProngsContributorsPV;
    float const multZeqFT0A = hVtxZFT0A->Interpolate(0.0) * multFT0A / hVtxZFT0A->Interpolate(collision.posZ());
    float const multZeqFT0C = hVtxZFT0C->Interpolate(0.0) * multFT0C / hVtxZFT0C->Interpolate(collision.posZ());
    return multZeqFT0A + multZeqFT0C;
  }

  /// Evaluate centrality/multiplicity percentile using NTracksPV estimator
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <o2::hf_centrality::HasNTracksPvCent Coll>
  float getZEqMultColl(const Coll& collision, uint8_t nProngsContributorsPV)
  {
    auto multNTracksPV = collision.multNTracksPV() - nProngsContributorsPV;
    float const multZeqNTracksPV = hVtxZNTracks->Interpolate(0.0) * multNTracksPV / hVtxZNTracks->Interpolate(collision.posZ());
    return multZeqNTracksPV;
  }

  /// Default case if no centrality/multiplicity estimator is provided
  /// \param candidate is candidate
  /// \return dummy value for centrality/multiplicity percentile of the collision
  template <typename Coll>
  float getZEqMultColl(const Coll&, uint8_t)
  {
    return -1.f;
  }

  /// Evaluate centrality/multiplicity percentile (centrality estimator is automatically selected based on the used table)
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll, typename CandDs>
  float evaluateCentralityColl(const Coll& collision, const CandDs& candidate)
  {
    if (fillPercentiles) {
      return o2::hf_centrality::getCentralityColl<Coll>(collision);
    }
    return getZEqMultColl<Coll>(collision, candidate.nProngsContributorsPV());
  }

  /// Evaluate centrality/multiplicity percentile (centrality estimator is automatically selected based on the used table)
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll>
  float evaluateCentralityColl(const Coll& collision)
  {
    if (fillPercentiles) {
      return o2::hf_centrality::getCentralityColl<Coll>(collision);
    }
    return getZEqMultColl<Coll>(collision, 0);
  }

  /// Evaluate centrality/multiplicity percentile
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision associated to the candidate
  template <typename Coll, typename T1>
  float evaluateCentralityCand(const T1& candidate)
  {
    return evaluateCentralityColl(candidate.template collision_as<Coll>(), candidate);
  }

  /// Fill histograms of quantities independent from the daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <typename T1>
  void fillHisto(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    std::get<TH1Ptr>(histosPtr[dataType]["hPt"])->Fill(pt);
    std::get<TH1Ptr>(histosPtr[dataType]["hPtProng0"])->Fill(candidate.ptProng0());
    std::get<TH1Ptr>(histosPtr[dataType]["hPtProng1"])->Fill(candidate.ptProng1());
    std::get<TH1Ptr>(histosPtr[dataType]["hPtProng2"])->Fill(candidate.ptProng2());
    std::get<TH2Ptr>(histosPtr[dataType]["hEta"])->Fill(candidate.eta(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hCt"])->Fill(HfHelper::ctDs(candidate), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hDecayLength"])->Fill(candidate.decayLength(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hDecayLengthXY"])->Fill(candidate.decayLengthXY(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hNormalisedDecayLengthXY"])->Fill(candidate.decayLengthXYNormalised(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hCPA"])->Fill(candidate.cpa(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hCPAxy"])->Fill(candidate.cpaXY(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hImpactParameterXY"])->Fill(candidate.impactParameterXY(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hMaxNormalisedDeltaIP"])->Fill(candidate.maxNormalisedDeltaIP(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hImpactParameterProngSqSum"])->Fill(candidate.impactParameterProngSqSum(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hDecayLengthError"])->Fill(candidate.errorDecayLength(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hDecayLengthXYError"])->Fill(candidate.errorDecayLengthXY(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hImpactParameterError"])->Fill(candidate.errorImpactParameter0(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hImpactParameterError"])->Fill(candidate.errorImpactParameter1(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hImpactParameterError"])->Fill(candidate.errorImpactParameter2(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hd0Prong0"])->Fill(candidate.impactParameter0(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hd0Prong1"])->Fill(candidate.impactParameter1(), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hd0Prong2"])->Fill(candidate.impactParameter2(), pt);
  }

  /// Fill mass sparse if ML information is present
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  /// \param finalState is either KKPi or PiKK, as defined in FinalState enum
  template <bool IsMc, typename Coll, HasDsMlInfo Cand>
  void fillSparse(const Cand& candidate, DataType dataType, FinalState finalState)
  {
    auto mass = finalState == FinalState::KKPi ? HfHelper::invMassDsToKKPi(candidate) : HfHelper::invMassDsToPiKK(candidate);
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
      if (storeOccupancy) {
        std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), outputMl[0], outputMl[1], outputMl[2], o2::hf_occupancy::getOccupancyColl(candidate.template collision_as<Coll>(), occEstimator));
        return;
      }
      std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), outputMl[0], outputMl[1], outputMl[2]);
      return;
    }
    if constexpr (IsMc) {
      if (dataType == DataType::McDsNonPrompt || dataType == DataType::McDplusNonPrompt) {
        if (storeOccupancy) {
          std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib(), candidate.ptBhadMotherPart(), getBHadMotherFlag(candidate.pdgBhadMotherPart()), o2::hf_occupancy::getOccupancyColl(candidate.template collision_as<Coll>(), occEstimator));
          return;
        }
        std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib(), candidate.ptBhadMotherPart(), getBHadMotherFlag(candidate.pdgBhadMotherPart()));
        return;
      }
      if (storeOccupancy) {
        std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib(), o2::hf_occupancy::getOccupancyColl(candidate.template collision_as<Coll>(), occEstimator));
        return;
      }
      std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), outputMl[0], outputMl[1], outputMl[2], candidate.template collision_as<Coll>().numContrib());
      return;
    }
  }

  /// Fill mass sparse if ML information is not present
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  /// \param finalState is either KKPi or PiKK, as defined in FinalState enum
  template <bool IsMc, typename Coll, typename Cand>
  void fillSparse(const Cand& candidate, DataType dataType, FinalState finalState)
  {
    auto mass = finalState == FinalState::KKPi ? HfHelper::invMassDsToKKPi(candidate) : HfHelper::invMassDsToPiKK(candidate);
    auto pt = candidate.pt();

    if (dataType == DataType::Data) { // If data do not fill PV contributors in sparse
      if (storeOccupancy) {
        std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), o2::hf_occupancy::getOccupancyColl(candidate.template collision_as<Coll>(), occEstimator));
        return;
      }
      std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate));
      return;
    }
    if constexpr (IsMc) {
      if (dataType == DataType::McDsNonPrompt || dataType == DataType::McDplusNonPrompt) {
        if (storeOccupancy) {
          std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), candidate.template collision_as<Coll>().numContrib(), candidate.ptBhadMotherPart(), getBHadMotherFlag(candidate.pdgBhadMotherPart()), o2::hf_occupancy::getOccupancyColl(candidate.template collision_as<Coll>(), occEstimator));
          return;
        }
        std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), candidate.template collision_as<Coll>().numContrib(), candidate.ptBhadMotherPart(), getBHadMotherFlag(candidate.pdgBhadMotherPart()));
        return;
      }
      if (storeOccupancy) {
        std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), candidate.template collision_as<Coll>().numContrib(), o2::hf_occupancy::getOccupancyColl(candidate.template collision_as<Coll>(), occEstimator));
        return;
      }
      std::get<THnSparsePtr>(histosPtr[dataType]["hSparseMass"])->Fill(mass, pt, evaluateCentralityCand<Coll>(candidate), candidate.template collision_as<Coll>().numContrib());
      return;
    }
  }

  template <bool IsMc, typename Coll, typename Cand>
  void fillMiniTrees(const Cand& candidate, DataType dataType, FinalState finalState)
  {
    auto mass = finalState == FinalState::KKPi ? HfHelper::invMassDsToKKPi(candidate) : HfHelper::invMassDsToPiKK(candidate);
    auto pt = candidate.pt();

    hfCandDsMinis(mass, pt, evaluateCentralityCand<Coll>(candidate));
    if (miniTrees.extendWithDecayLength) {
      hfCandDsDlMinis(candidate.decayLength(), candidate.decayLengthXY(), candidate.decayLengthNormalised(), candidate.decayLengthXYNormalised());
    }
    if (miniTrees.extendWithImpactParameter) {
      hfCandDsD0Minis(candidate.impactParameterXY());
    }
    if constexpr (IsMc) {
      hfCandDsMcMinis(dataType);
    }
  }

  /// Fill histograms of quantities for the KKPi daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <bool IsMc, typename Coll, typename T1>
  void fillHistoKKPi(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    fillSparse<IsMc, Coll>(candidate, dataType, FinalState::KKPi);

    std::get<TH2Ptr>(histosPtr[dataType]["hCos3PiK"])->Fill(HfHelper::cos3PiKDsToKKPi(candidate), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hAbsCos3PiK"])->Fill(HfHelper::absCos3PiKDsToKKPi(candidate), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hDeltaMassPhi"])->Fill(HfHelper::deltaMassPhiDsToKKPi(candidate), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hMassKK"])->Fill(HfHelper::massKKPairDsToKKPi(candidate), pt);
  }

  /// Fill histograms of quantities for the PiKK daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param dataType is data class, as defined in DataType enum
  template <bool IsMc, typename Coll, typename T1>
  void fillHistoPiKK(const T1& candidate, DataType dataType)
  {
    auto pt = candidate.pt();
    fillSparse<IsMc, Coll>(candidate, dataType, FinalState::PiKK);

    std::get<TH2Ptr>(histosPtr[dataType]["hCos3PiK"])->Fill(HfHelper::cos3PiKDsToPiKK(candidate), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hAbsCos3PiK"])->Fill(HfHelper::absCos3PiKDsToPiKK(candidate), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hDeltaMassPhi"])->Fill(HfHelper::deltaMassPhiDsToPiKK(candidate), pt);
    std::get<TH2Ptr>(histosPtr[dataType]["hMassKK"])->Fill(HfHelper::massKKPairDsToPiKK(candidate), pt);
  }

  /// Fill MC histograms at reconstruction level
  /// \param candidate is candidate
  /// \param mcParticles are particles with MC information
  /// \param whichSpeciesDecay defines which histogram to fill
  template <typename Coll, typename T1>
  void fillHistoMCRec(const T1& candidate, const CandDsMcGen& mcParticles, DataType dataType)
  {

    int id = o2::constants::physics::Pdg::kDS;
    if (dataType == DataType::McDplusPrompt || dataType == DataType::McDplusNonPrompt || dataType == DataType::McDplusBkg) {
      id = o2::constants::physics::Pdg::kDPlus;
    } else if (dataType == DataType::McLcBkg) {
      id = o2::constants::physics::Pdg::kLambdaCPlus;
    }

    auto indexMother = RecoDecay::getMother(mcParticles,
                                            candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<CandDsMcGen>(),
                                            id, true);

    if (indexMother != -1) {

      auto pt = candidate.pt(); // rec. level pT

      if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
        auto yCand = candidate.y(HfHelper::invMassDsToKKPi(candidate));
        if (yCandRecoMax >= 0. && std::abs(yCand) > yCandRecoMax) {
          return;
        }
        fillHisto(candidate, dataType);
        fillHistoKKPi<true, Coll>(candidate, dataType);
        if (miniTrees.produceMiniTrees) {
          fillMiniTrees<true, Coll>(candidate, dataType, FinalState::KKPi);
        }
        if (TESTBIT(candidate.isSelDsToKKPi(), aod::SelectionStep::RecoSkims)) {
          std::get<TH2Ptr>(histosPtr[dataType]["hPtVsYRecoSkim"])->Fill(pt, yCand);
        }
        if (TESTBIT(candidate.isSelDsToKKPi(), aod::SelectionStep::RecoTopol)) {
          std::get<TH2Ptr>(histosPtr[dataType]["hPtVsYRecoTopol"])->Fill(pt, yCand);
        }
        if (TESTBIT(candidate.isSelDsToKKPi(), aod::SelectionStep::RecoPID)) {
          std::get<TH2Ptr>(histosPtr[dataType]["hPtVsYRecoPID"])->Fill(pt, yCand);
        }
      }
      if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
        auto yCand = candidate.y(HfHelper::invMassDsToPiKK(candidate));
        if (yCandRecoMax >= 0. && std::abs(yCand) > yCandRecoMax) {
          return;
        }
        fillHisto(candidate, dataType);
        fillHistoPiKK<true, Coll>(candidate, dataType);
        if (miniTrees.produceMiniTrees) {
          fillMiniTrees<true, Coll>(candidate, dataType, FinalState::PiKK);
        }

        if (TESTBIT(candidate.isSelDsToPiKK(), aod::SelectionStep::RecoSkims)) {
          std::get<TH2Ptr>(histosPtr[dataType]["hPtVsYRecoSkim"])->Fill(pt, yCand);
        }
        if (TESTBIT(candidate.isSelDsToPiKK(), aod::SelectionStep::RecoTopol)) {
          std::get<TH2Ptr>(histosPtr[dataType]["hPtVsYRecoTopol"])->Fill(pt, yCand);
        }
        if (TESTBIT(candidate.isSelDsToPiKK(), aod::SelectionStep::RecoPID)) {
          std::get<TH2Ptr>(histosPtr[dataType]["hPtVsYRecoPID"])->Fill(pt, yCand);
        }
      }
    }
  }

  template <typename Coll, typename CandDs>
  void runDataAnalysisPerCandidate(CandDs const& candidate)
  {
    if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
      if (yCandRecoMax >= 0. && std::abs(candidate.y(HfHelper::invMassDsToKKPi(candidate))) > yCandRecoMax) {
        return;
      }
      fillHisto(candidate, DataType::Data);
      fillHistoKKPi<false, Coll>(candidate, DataType::Data);
      if (miniTrees.produceMiniTrees) {
        fillMiniTrees<true, Coll>(candidate, DataType::Data, FinalState::KKPi);
      }
    }
    if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
      if (yCandRecoMax >= 0. && std::abs(candidate.y(HfHelper::invMassDsToPiKK(candidate))) > yCandRecoMax) {
        return;
      }
      fillHisto(candidate, DataType::Data);
      fillHistoPiKK<false, Coll>(candidate, DataType::Data);
      if (miniTrees.produceMiniTrees) {
        fillMiniTrees<true, Coll>(candidate, DataType::Data, FinalState::PiKK);
      }
    }
  }

  template <typename Coll, typename CandDs>
  void runMcAnalysisPerCandidate(CandDs const& candidate,
                                 CandDsMcGen const& mcParticles)
  {
    // MC rec.
    std::array<MemberFunctionPointer<CandDs>, 6> isOfType = {// Contains the functions to check if the candidate is of a certain type
                                                             &HfTaskDs::isDsPrompt<CandDs>,
                                                             &HfTaskDs::isDsNonPrompt<CandDs>,
                                                             &HfTaskDs::isDplusPrompt<CandDs>,
                                                             &HfTaskDs::isDplusNonPrompt<CandDs>,
                                                             &HfTaskDs::isDplusBkg<CandDs>,
                                                             &HfTaskDs::isLcBkg<CandDs>};

    bool isBkg = true;
    for (int i = DataType::McDsPrompt; i <= DataType::McLcBkg; i++) { // Check what type of MC signal candidate it is, and fill the corresponding histograms
      if ((this->*isOfType[i - DataType::McDsPrompt])(candidate)) {
        isBkg = false;
        fillHistoMCRec<Coll>(candidate, mcParticles, static_cast<DataType>(i));
        break;
      }
    }
    if (isBkg && fillMcBkgHistos) {

      if (candidate.isSelDsToKKPi() >= selectionFlagDs || candidate.isSelDsToPiKK() >= selectionFlagDs) {
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) { // KKPi
          if (yCandRecoMax >= 0. && std::abs(candidate.y(HfHelper::invMassDsToKKPi(candidate))) > yCandRecoMax) {
            return;
          }
          fillHisto(candidate, DataType::McBkg);
          fillHistoKKPi<true, Coll>(candidate, DataType::McBkg);
        }
        if (candidate.isSelDsToPiKK() >= selectionFlagDs) { // PiKK
          if (yCandRecoMax >= 0. && std::abs(candidate.y(HfHelper::invMassDsToPiKK(candidate))) > yCandRecoMax) {
            return;
          }
          fillHisto(candidate, DataType::McBkg);
          fillHistoPiKK<true, Coll>(candidate, DataType::McBkg);
        }
      }
    }

    // TODO: add histograms for reflections
  }

  template <typename Coll>
  void fillMcGenHistosSparse(CandDsMcGen const& mcParticles,
                             Coll const& recoCollisions)
  {
    // MC gen.
    for (const auto& particle : mcParticles) {

      if (std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK || std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKK) {
        const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
        if (particle.flagMcDecayChanGen() == channelsResonant[Mother::Ds][decayChannel] || (fillDplusMc && particle.flagMcDecayChanGen() == channelsResonant[Mother::Dplus][decayChannel])) {
          auto pt = particle.pt();
          double y{0.f};

          unsigned maxNumContrib = 0;
          for (const auto& recCol : recoCollsPerMcColl) {
            maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
          }
          float const cent = o2::hf_centrality::getCentralityGenColl(recoCollsPerMcColl);
          float occ{-1.};
          if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
            occ = o2::hf_occupancy::getOccupancyGenColl(recoCollsPerMcColl, occEstimator);
          }

          if (particle.flagMcDecayChanGen() == channelsResonant[Mother::Ds][decayChannel]) {
            y = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDS);
            if (yCandGenMax >= 0. && std::abs(y) > yCandGenMax) {
              continue;
            }

            if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
              std::get<TH1Ptr>(histosPtr[DataType::McDsPrompt]["hPtGen"])->Fill(pt); // gen. level pT
              std::get<TH1Ptr>(histosPtr[DataType::McDsPrompt]["hEtaGen"])->Fill(particle.eta());
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                std::get<THnSparsePtr>(histosPtr[DataType::McDsPrompt]["hSparseGen"])->Fill(pt, y, maxNumContrib, cent, occ);
              } else {
                std::get<THnSparsePtr>(histosPtr[DataType::McDsPrompt]["hSparseGen"])->Fill(pt, y, maxNumContrib, cent);
              }
            }
            if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              std::get<TH1Ptr>(histosPtr[DataType::McDsNonPrompt]["hPtGen"])->Fill(pt); // gen. level pT
              std::get<TH1Ptr>(histosPtr[DataType::McDsNonPrompt]["hEtaGen"])->Fill(particle.eta());
              auto bHadMother = mcParticles.rawIteratorAt(particle.idxBhadMotherPart() - mcParticles.offset());
              int const flagGenB = getBHadMotherFlag(bHadMother.pdgCode());
              float const ptGenB = bHadMother.pt();
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                std::get<THnSparsePtr>(histosPtr[DataType::McDsNonPrompt]["hSparseGen"])->Fill(pt, y, maxNumContrib, cent, occ, ptGenB, flagGenB);
              } else {
                std::get<THnSparsePtr>(histosPtr[DataType::McDsNonPrompt]["hSparseGen"])->Fill(pt, y, maxNumContrib, cent, ptGenB, flagGenB);
              }
            }
          } else if (fillDplusMc) {
            y = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus);
            if (yCandGenMax >= 0. && std::abs(y) > yCandGenMax) {
              continue;
            }
            if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
              std::get<TH1Ptr>(histosPtr[DataType::McDplusPrompt]["hPtGen"])->Fill(pt); // gen. level pT
              std::get<TH1Ptr>(histosPtr[DataType::McDplusPrompt]["hEtaGen"])->Fill(particle.eta());
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                std::get<THnSparsePtr>(histosPtr[DataType::McDplusPrompt]["hSparseGen"])->Fill(pt, y, maxNumContrib, cent, occ);
              } else {
                std::get<THnSparsePtr>(histosPtr[DataType::McDplusPrompt]["hSparseGen"])->Fill(pt, y, maxNumContrib, cent);
              }
            }
            if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              std::get<TH1Ptr>(histosPtr[DataType::McDplusNonPrompt]["hPtGen"])->Fill(pt); // gen. level pT
              std::get<TH1Ptr>(histosPtr[DataType::McDplusNonPrompt]["hEtaGen"])->Fill(particle.eta());
              auto bHadMother = mcParticles.rawIteratorAt(particle.idxBhadMotherPart() - mcParticles.offset());
              int const flagGenB = getBHadMotherFlag(bHadMother.pdgCode());
              float const ptGenB = bHadMother.pt();
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                std::get<THnSparsePtr>(histosPtr[DataType::McDplusNonPrompt]["hSparseGen"])->Fill(pt, y, maxNumContrib, cent, occ, ptGenB, flagGenB);
              } else {
                std::get<THnSparsePtr>(histosPtr[DataType::McDplusNonPrompt]["hSparseGen"])->Fill(pt, y, maxNumContrib, cent, ptGenB, flagGenB);
              }
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
    int const numPvContributors = collision.numContrib();
    float const centrality = evaluateCentralityColl(collision);
    std::get<TH2Ptr>(histosPtr[DataType::Data]["hNPvContribAll"])->Fill(numPvContributors, centrality);
    for (int i = 0; i < DataType::kDataTypes; i++) {
      if (i == DataType::McBkg && !fillMcBkgHistos) {
        continue;
      }
      if (nCandsPerType[i]) {
        std::get<TH2Ptr>(histosPtr[i]["hNPvContribCands"])->Fill(numPvContributors, centrality);
      }
      if (nCandsInSignalRegionDsPerType[i]) {
        std::get<TH2Ptr>(histosPtr[i]["hNPvContribCandsInSignalRegionDs"])->Fill(numPvContributors, centrality);
      }
      if (nCandsInSignalRegionDplusPerType[i]) {
        std::get<TH2Ptr>(histosPtr[i]["hNPvContribCandsInSignalRegionDplus"])->Fill(numPvContributors, centrality);
      }
    }
  }

  template <typename Coll, typename CandsDs>
  void runDataAnalysisPerCollision(const Coll& collisions, const CandsDs& candsDs)
  {
    for (const auto& collision : collisions) {
      /* check the previous run number */
      const auto& bc = collision.bc();
      if (bc.runNumber() != mRunNumber) {
        mRunNumber = bc.runNumber(); // mark this run as at least tried
        if (ccdbConfig.reconstructionPass.value.empty()) {
          lCalibObjects = ccdb->getForRun<TList>(ccdbConfig.ccdbPath, mRunNumber);
        } else if (ccdbConfig.reconstructionPass.value == "metadata") {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
          LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
          lCalibObjects = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, mRunNumber, metadata);
        } else {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = ccdbConfig.reconstructionPass.value;
          LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", ccdbConfig.reconstructionPass.value);
          lCalibObjects = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, mRunNumber, metadata);
        }

        if (lCalibObjects) {
          LOG(info) << "CCDB objects loaded successfully";
          hVtxZFT0A = dynamic_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFT0A"));
          hVtxZFT0C = dynamic_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFT0C"));
          hVtxZNTracks = dynamic_cast<TProfile*>(lCalibObjects->FindObject("hVtxZNTracksPV"));
          lCalibLoaded = true;
          // Capture error
          if (!hVtxZFT0A || !hVtxZFT0C || !hVtxZNTracks) {
            LOGF(error, "Problem loading CCDB objects! Please check");
            lCalibLoaded = false;
          }
        } else {
          LOGF(error, "Problem loading CCDB object! Please check");
          lCalibLoaded = false;
        }
      }

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

  template <typename Coll, typename CandsDs, typename CandDsMcGen>
  void runMcAnalysisPerCollision(const Coll& collisions,
                                 const CandsDs& candsDs,
                                 const CandDsMcGen& mcParticles)
  {
    for (const auto& collision : collisions) {
      /* check the previous run number */
      const auto& bc = collision.bc();
      if (bc.runNumber() != mRunNumber) {
        mRunNumber = bc.runNumber(); // mark this run as at least tried
        if (ccdbConfig.reconstructionPass.value.empty()) {
          lCalibObjects = ccdb->getForRun<TList>(ccdbConfig.ccdbPath, mRunNumber);
        } else if (ccdbConfig.reconstructionPass.value == "metadata") {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
          LOGF(info, "Loading CCDB for reconstruction pass (from metadata): %s", metadataInfo.get("RecoPassName"));
          lCalibObjects = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, mRunNumber, metadata);
        } else {
          std::map<std::string, std::string> metadata;
          metadata["RecoPassName"] = ccdbConfig.reconstructionPass.value;
          LOGF(info, "Loading CCDB for reconstruction pass (from provided argument): %s", ccdbConfig.reconstructionPass.value);
          lCalibObjects = ccdb->getSpecificForRun<TList>(ccdbConfig.ccdbPath, mRunNumber, metadata);
        }

        if (lCalibObjects) {
          LOG(info) << "CCDB objects loaded successfully";
          hVtxZFT0A = dynamic_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFT0A"));
          hVtxZFT0C = dynamic_cast<TProfile*>(lCalibObjects->FindObject("hVtxZFT0C"));
          hVtxZNTracks = dynamic_cast<TProfile*>(lCalibObjects->FindObject("hVtxZNTracksPV"));
          lCalibLoaded = true;
          // Capture error
          if (!hVtxZFT0A || !hVtxZFT0C || !hVtxZNTracks) {
            LOGF(error, "Problem loading CCDB objects! Please check");
            lCalibLoaded = false;
          }
        } else {
          LOGF(error, "Problem loading CCDB object! Please check");
          lCalibLoaded = false;
        }
      }

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
    fillMcGenHistosSparse(mcParticles, collisions);
  }

  void processDataWithCentFT0C(CollisionsWithFT0C const& collisions,
                               CandDsData const& candsDs,
                               aod::BCs const&,
                               aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentFT0C, "Process data w/o ML information on Ds, with information on centrality from FT0C", false);

  void processDataWithCentFT0M(CollisionsWithFT0M const& collisions,
                               CandDsData const& candsDs,
                               aod::BCs const&,
                               aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentFT0M, "Process data w/o ML information on Ds, with information on centrality from FT0M", false);

  void processDataWithCentNTracksPV(CollisionsWithNTracksPV const& collisions,
                                    CandDsData const& candsDs,
                                    aod::BCs const&,
                                    aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithCentNTracksPV, "Process data w/o ML information on Ds, with information on centrality from NTracksPV", false);

  void processData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   CandDsData const& candsDs,
                   aod::BCs const&,
                   aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs);
  }
  PROCESS_SWITCH(HfTaskDs, processData, "Process data w/o ML information on Ds, w/o information on centrality", true);

  void processDataWithMlAndCentFT0C(CollisionsWithFT0C const& collisions,
                                    CandDsDataWithMl const& candsDs,
                                    aod::BCs const&,
                                    aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentFT0C, "Process data with ML information on Ds, with information on centrality from FT0C", false);

  void processDataWithMlAndCentFT0M(CollisionsWithFT0M const& collisions,
                                    CandDsDataWithMl const& candsDs,
                                    aod::BCs const&,
                                    aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentFT0M, "Process data with ML information on Ds, with information on centrality from FT0M", false);

  void processDataWithMlAndCentNTracksPV(CollisionsWithNTracksPV const& collisions,
                                         CandDsDataWithMl const& candsDs,
                                         aod::BCs const&,
                                         aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMlAndCentNTracksPV, "Process data with ML information on Ds, with information on centrality", false);

  void processDataWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                         CandDsDataWithMl const& candsDs,
                         aod::BCs const&,
                         aod::Tracks const&)
  {
    runDataAnalysisPerCollision(collisions, candsDs);
  }
  PROCESS_SWITCH(HfTaskDs, processDataWithMl, "Process data with ML information on Ds, w/o information on centrality", false);

  void processMcWithCentFT0C(CollisionsMcWithFT0C const& collisions,
                             CandDsMcReco const& candsDs,
                             CandDsMcGen const& mcParticles,
                             aod::BCs const&,
                             aod::McCollisions const&,
                             aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CollisionsMcWithFT0C>(collisions, candsDs, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentFT0C, "Process MC w/o ML information on Ds, with information on centrality from FT0C", false);

  void processMcWithCentFT0M(CollisionsMcWithFT0M const& collisions,
                             CandDsMcReco const& candsDs,
                             CandDsMcGen const& mcParticles,
                             aod::BCs const&,
                             aod::McCollisions const&,
                             aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CollisionsMcWithFT0M>(collisions, candsDs, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentFT0M, "Process MC w/o ML information on Ds, with information on centrality from FT0M", false);

  void processMcWithCentNTracksPV(CollisionsMcWithNTracksPV const& collisions,
                                  CandDsMcReco const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::BCs const&,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CollisionsMcWithNTracksPV>(collisions, candsDs, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithCentNTracksPV, "Process MC w/o ML information on Ds, with information on centrality from NTracksPV", false);

  void processMc(CollisionsMc const& collisions,
                 CandDsMcReco const& candsDs,
                 CandDsMcGen const& mcParticles,
                 aod::BCs const&,
                 aod::McCollisions const&,
                 aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CollisionsMc>(collisions, candsDs, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDs, processMc, "Process MC w/o ML information on Ds, w/o information on centrality", false);

  void processMcWithMlAndCentFT0C(CollisionsMcWithFT0C const& collisions,
                                  CandDsMcRecoWithMl const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::BCs const&,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CollisionsMcWithFT0C>(collisions, candsDs, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentFT0C, "Process MC with ML information on Ds, with information on centrality from FT0C", false);

  void processMcWithMlAndCentFT0M(CollisionsMcWithFT0M const& collisions,
                                  CandDsMcRecoWithMl const& candsDs,
                                  CandDsMcGen const& mcParticles,
                                  aod::BCs const&,
                                  aod::McCollisions const&,
                                  aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CollisionsMcWithFT0M>(collisions, candsDs, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentFT0M, "Process MC with ML information on Ds, with information on centrality from FT0M", false);

  void processMcWithMlAndCentNTracksPV(CollisionsMcWithNTracksPV const& collisions,
                                       CandDsMcRecoWithMl const& candsDs,
                                       CandDsMcGen const& mcParticles,
                                       aod::BCs const&,
                                       aod::McCollisions const&,
                                       aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CollisionsMcWithNTracksPV>(collisions, candsDs, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMlAndCentNTracksPV, "Process MC with ML information on Ds, with information on centrality from NTracksPV", false);

  void processMcWithMl(CollisionsMc const& collisions,
                       CandDsMcRecoWithMl const& candsDs,
                       CandDsMcGen const& mcParticles,
                       aod::BCs const&,
                       aod::McCollisions const&,
                       aod::TracksWMc const&)
  {
    runMcAnalysisPerCollision<CollisionsMc>(collisions, candsDs, mcParticles);
  }
  PROCESS_SWITCH(HfTaskDs, processMcWithMl, "Process MC with ML information on Ds, w/o information on centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<HfTaskDs>(cfgc)};
}
