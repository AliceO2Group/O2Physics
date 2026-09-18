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
//
// Step-by-step event selection efficiency for MC particle-level jets.
// Two ordering variants exposed as separate process functions:
//   - CollRecoFirst: reco collision required first, BC bits read from the reco-coll EvSel bitmask.
//   - BcBitsFirst:   BC bits read from the MC truth BC, then truth-side selections, then reco at the end.
//
/// \author Joonsuk Bae <joonsuk.bae@cern.ch>
/// \author Wooseok Ham <wooseok.ham@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsITSMFT/DPLAlpideParam.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>

#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetCrossSectionEfficiency {

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "selTVX | selMC | selMCFull | sel8 | sel8Full | sel8FullPbPb"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "reject min-bias gap events from hybrid MB+JJ MC productions"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0f, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0f, "maximum centrality"};
  Configurable<bool> checkCentFT0M{"checkCentFT0M", false, "0: centFT0C as default, 1: use centFT0M estimator"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9f, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9f, "maximum eta acceptance for tracks"};
  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4f, "resolution parameter for histograms without radius"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0f, "minimum jet pT"};
  Configurable<double> jetPtMax{"jetPtMax", 200.0, "set jet pT bin max"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.5f, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.5f, "maximum jet pseudorapidity"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0f, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMinMCP{"leadingConstituentPtMinMCP", -99.0f, "minimum pT selection on MCP jet constituent"};
  Configurable<float> leadingConstituentPtMaxMCP{"leadingConstituentPtMaxMCP", 9999.0f, "maximum pT selection on MCP jet constituent"};

  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> acceptSplitCollisions{"acceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0f, "maximum jet pT in units of pTHat for outlier rejection at MC particle level"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0f, "minimum pTHat (drops events with pTHat below this)"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0f, "exponent in the back-calculation of pTHat from the event weight"};
  Configurable<float> simPtRef{"simPtRef", 10.0f, "reference pT for the back-calculation of pTHat from the event weight"};
  Configurable<float> ptHardCalcMethodSwitch{"ptHardCalcMethodSwitch", 999.0f, "sentinel and threshold for the stored-ptHard branch (set <=0 to force weight-derived)"};
  Configurable<bool> applyRCT{"applyRCT", true, "apply RCT_pass step in the cascade (false: force the RCT step to pass)"};
  Configurable<std::string> rctSelectionsLabel{"rctSelectionsLabel", "CBT_hadronPID", "RCT selection preset name (see RCTSelectionFlags)"};

  // The reconstructed occupancy bits cannot be evaluated before hasColl.
  // These configurables define the truth-only approximation used by sel8FullPbPb.
  Configurable<std::string> truthRofCcdbUrl{"truthRofCcdbUrl", "http://alice-ccdb.cern.ch", "CCDB URL used to retrieve ITS ROF parameters for sel8FullPbPb truth selections"};
  Configurable<int> truthRofOffsetInBC{"truthRofOffsetInBC", -1, "ITS ROF bias in BC for sel8FullPbPb truth selections; -1 retrieves DPLAlpideParam from CCDB"};
  Configurable<int> truthRofLengthInBC{"truthRofLengthInBC", -1, "ITS ROF length in BC for sel8FullPbPb truth selections; -1 retrieves DPLAlpideParam from CCDB"};
  Configurable<float> truthTimeRangeStandardMinUs{"truthTimeRangeStandardMinUs", -4.0f, "truth lower time range in us; follows kNoCollInTimeRangeStandard"};
  Configurable<float> truthTimeRangeStandardMaxUs{"truthTimeRangeStandardMaxUs", 2.0f, "truth upper time range in us; follows kNoCollInTimeRangeStandard"};
  Configurable<float> truthTimeRangeNarrowUs{"truthTimeRangeNarrowUs", 0.25f, "truth narrow time range in us; follows kNoCollInTimeRangeStandard"};
  Configurable<float> truthFT0CActivityMinForTrackProxy{"truthFT0CActivityMinForTrackProxy", 0.0f, "minimum multFT0C truth activity used as an ITS-track-presence proxy"};
  Configurable<float> truthFT0CActivityThresholdTimeRange{"truthFT0CActivityThresholdTimeRange", 8000.0f, "multFT0C truth proxy threshold for NoCollInTimeRangeStandard; reco uses FT0C amplitude, so validate or tune this value"};
  Configurable<float> truthFT0CActivityThresholdROF{"truthFT0CActivityThresholdROF", 5000.0f, "multFT0C truth proxy threshold for NoCollInRofStandard; reco uses FT0C amplitude, so validate or tune this value"};
  Configurable<float> truthRofCloseVzMax{"truthRofCloseVzMax", 0.3f, "maximum truth |delta z| in cm for the NoCollInRofStandard track-presence proxy"};

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  uint64_t rctMask = 0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int cachedTruthRofRun = std::numeric_limits<int>::min();
  int64_t cachedTruthRofOffsetInBC = -1;
  int64_t cachedTruthRofLengthInBC = -1;

  bool applyTFB = true;
  bool applyROFB = true;
  bool applySBP = true;
  bool applyNoCollInTimeRangeStandard = false;
  bool applyNoCollInRofStandard = false;
  bool isSel8FullPbPb = false;

  std::vector<std::string> collRecoFirstLabels;
  std::vector<std::string> bcBitsFirstLabels;

  enum AcceptSplitCollisionsOptions {
    NonSplitOnly = 0,
    SplitOkCheckAnyAssocColl,      // 1
    SplitOkCheckFirstAssocCollOnly // 2
  };

  enum BinPbPbTruthSelectionOnly {
    PbPbTruthSelectionOnlyInel = 1,
    PbPbTruthSelectionOnlyRct,
    PbPbTruthSelectionOnlyTvx,
    PbPbTruthSelectionOnlyNoTimeFrameBorder,
    PbPbTruthSelectionOnlyNoItsRofBorder,
    PbPbTruthSelectionOnlySelection,
    PbPbTruthSelectionOnlyHasCollision,
    PbPbTruthSelectionOnlyVertexZ,
    PbPbTruthSelectionOnlyNoSplit,
    PbPbTruthSelectionOnlyNBins = PbPbTruthSelectionOnlyNoSplit
  };

  static constexpr float ConfigSwitchLow = -98.0f;
  static constexpr float ConfigSwitchHigh = 9998.0f;
  static constexpr float BrokenPtHardSentinel = 1.0f;
  static constexpr int MinITSClustersForOccupancy = 5;

  struct TruthPbPbSelections {
    bool valid = false;
    bool noCollInTimeRangeStandard = false;
    bool noCollInRofStandard = false;
  };

  using JetCollisionsMCDWithParent = soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs>;
  using CollisionsWithEvSels = soa::Join<aod::Collisions, aod::EvSels>;
  using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;

  Partition<FullTracksIU> pvTracks = ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Preslice<FullTracksIU> pvTracksPerCollision = aod::track::collisionId;

  enum EventSelectionPreset {
    PresetSelTvx = 0,
    PresetSelMc,
    PresetSelMcFull,
    PresetSel8,
    PresetSel8Full,
    PresetSel8FullPbPb,
    PresetInvalid
  };

  static EventSelectionPreset getEventSelectionPreset(const std::string& preset)
  {
    if (preset == "selTVX") {
      return PresetSelTvx;
    }
    if (preset == "selMC") {
      return PresetSelMc;
    }
    if (preset == "selMCFull") {
      return PresetSelMcFull;
    }
    if (preset == "sel8") {
      return PresetSel8;
    }
    if (preset == "sel8Full") {
      return PresetSel8Full;
    }
    if (preset == "sel8FullPbPb") {
      return PresetSel8FullPbPb;
    }
    return PresetInvalid;
  }

  Preslice<aod::JetMcCollisions> mcCollsPerBC = aod::jmccollision::bcId;

  void init(InitContext&)
  {
    if (!(acceptSplitCollisions == NonSplitOnly || acceptSplitCollisions == SplitOkCheckAnyAssocColl || acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly)) {
      LOGF(fatal, "Configurable acceptSplitCollisions has wrong input value; stopping workflow");
    }

    rctChecker.init(static_cast<std::string>(rctSelectionsLabel));
    rctMask = rctChecker.value();

    switch (getEventSelectionPreset(static_cast<std::string>(eventSelections))) {
      case PresetSelTvx:
        applyTFB = false;
        applyROFB = false;
        applySBP = false;
        break;
      case PresetSelMc:
        applyTFB = true;
        applyROFB = false;
        applySBP = false;
        break;
      case PresetSelMcFull:
        applyTFB = true;
        applyROFB = false;
        applySBP = true;
        break;
      case PresetSel8:
        applyTFB = true;
        applyROFB = true;
        applySBP = false;
        break;
      case PresetSel8Full:
        applyTFB = true;
        applyROFB = true;
        applySBP = true;
        break;
      case PresetSel8FullPbPb:
        applyTFB = true;
        applyROFB = true;
        applySBP = false;
        applyNoCollInTimeRangeStandard = true;
        applyNoCollInRofStandard = true;
        isSel8FullPbPb = true;
        break;
      default:
        LOGF(fatal, "Configurable eventSelections=%s not supported; use selTVX, selMC, selMCFull, sel8, sel8Full, or sel8FullPbPb", static_cast<std::string>(eventSelections).c_str());
        break;
    }

    if (doprocessCrossSectionEfficiencyBcBitsFirst && isSel8FullPbPb &&
        (truthRofOffsetInBC < 0 || truthRofLengthInBC <= 0)) {
      ccdb->setURL(static_cast<std::string>(truthRofCcdbUrl));
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
    }

    // Preserve the upstream binning and order for existing presets.
    // sel8FullPbPb replaces the SBP stage with the two PbPb occupancy selections.
    collRecoFirstLabels = {"INEL", "+RCT_pass", "+hasRecoColl", "+|zReco|<10", "+noSplit", "+kTVX", "+kNoTFB", "+kNoITSROFB"};
    bcBitsFirstLabels = {"INEL", "+RCT_pass", "+kTVX(truth)", "+kNoTFB(truth)", "+kNoITSROFB(truth)"};
    if (isSel8FullPbPb) {
      collRecoFirstLabels.emplace_back("+kNoCollInTimeRangeStandard");
      collRecoFirstLabels.emplace_back("+kNoCollInRofStandard");
      bcBitsFirstLabels.emplace_back("+kNoCollInTimeRangeStandard(truth)");
      bcBitsFirstLabels.emplace_back("+kNoCollInRofStandard(truth)");
      bcBitsFirstLabels.emplace_back("+hasColl");
      bcBitsFirstLabels.emplace_back("+|zReco|<10");
      bcBitsFirstLabels.emplace_back("+noSplit");
    } else {
      collRecoFirstLabels.emplace_back("+kNoSBP");
      bcBitsFirstLabels.emplace_back("+kNoSBP(truth)");
      bcBitsFirstLabels.emplace_back("+hasColl");
      bcBitsFirstLabels.emplace_back("+|zReco|<10");
      bcBitsFirstLabels.emplace_back("+noSplit");
    }

    auto setAxisLabels = [](TAxis* axis, const std::vector<std::string>& labels) {
      for (size_t i = 0; i < labels.size(); ++i) {
        axis->SetBinLabel(static_cast<int>(i + 1), labels[i].c_str());
      }
    };

    AxisSpec jetPtAxis = {200, 0., jetPtMax, "#it{p}_{T} (GeV/#it{c})"};

    if (doprocessCrossSectionEfficiency) {
      AxisSpec axisSelectionCollRecoFirst = {static_cast<int>(collRecoFirstLabels.size()), 0.5, static_cast<double>(collRecoFirstLabels.size()) + 0.5, "event selection (CollRecoFirst)"};
      registry.add("h2_jet_pt_part_eventselection_collRecoFirst",
                   "part jet pT vs event selection (CollRecoFirst);#it{p}_{T,jet}^{part} (GeV/#it{c});event selection;counts",
                   {HistType::kTH2F, {jetPtAxis, axisSelectionCollRecoFirst}});
      setAxisLabels(registry.get<TH2>(HIST("h2_jet_pt_part_eventselection_collRecoFirst"))->GetYaxis(), collRecoFirstLabels);

      registry.add("h_mccollisions_eventselection_collRecoFirst",
                   "number of mc events vs event selection (CollRecoFirst);event selection;entries",
                   {HistType::kTH1F, {axisSelectionCollRecoFirst}});
      setAxisLabels(registry.get<TH1>(HIST("h_mccollisions_eventselection_collRecoFirst"))->GetXaxis(), collRecoFirstLabels);
    }

    if (doprocessCrossSectionEfficiencyBcBitsFirst) {
      AxisSpec axisSelectionBcBitsFirst = {static_cast<int>(bcBitsFirstLabels.size()), 0.5, static_cast<double>(bcBitsFirstLabels.size()) + 0.5, "event selection (BcBitsFirst)"};
      registry.add("h2_jet_pt_part_eventselection_bcBitsFirst",
                   "part jet pT vs event selection (BcBitsFirst);#it{p}_{T,jet}^{part} (GeV/#it{c});event selection;counts",
                   {HistType::kTH2F, {jetPtAxis, axisSelectionBcBitsFirst}});
      setAxisLabels(registry.get<TH2>(HIST("h2_jet_pt_part_eventselection_bcBitsFirst"))->GetYaxis(), bcBitsFirstLabels);

      registry.add("h_mccollisions_eventselection_bcBitsFirst",
                   "number of mc events vs event selection (BcBitsFirst);event selection;entries",
                   {HistType::kTH1F, {axisSelectionBcBitsFirst}});
      setAxisLabels(registry.get<TH1>(HIST("h_mccollisions_eventselection_bcBitsFirst"))->GetXaxis(), bcBitsFirstLabels);

      if (isSel8FullPbPb) {
        AxisSpec axisPbPbTruthSelectionOnly = {PbPbTruthSelectionOnlyNBins, 0.5, static_cast<double>(PbPbTruthSelectionOnlyNBins) + 0.5, "event selection (PbPb truth-only)"};
        std::vector<std::string> timeRangeOnlyLabels = {"INEL", "+RCT_pass", "+kTVX(truth)", "+kNoTFB(truth)", "+kNoITSROFB(truth)", "+kNoCollInTimeRangeStandard(truth)", "+hasColl", "+|zReco|<10", "+noSplit"};
        std::vector<std::string> rofOnlyLabels = {"INEL", "+RCT_pass", "+kTVX(truth)", "+kNoTFB(truth)", "+kNoITSROFB(truth)", "+kNoCollInRofStandard(truth)", "+hasColl", "+|zReco|<10", "+noSplit"};

        registry.add("h2_jet_pt_part_eventselection_bcBitsFirst_timeRangeOnly_truth",
                     "part jet pT vs TimeRange-only truth selection;#it{p}_{T,jet}^{part} (GeV/#it{c});event selection;weighted counts",
                     {HistType::kTH2F, {jetPtAxis, axisPbPbTruthSelectionOnly}});
        registry.add("h_mccollisions_eventselection_bcBitsFirst_timeRangeOnly_truth",
                     "weighted MC events vs TimeRange-only truth selection;event selection;weighted events",
                     {HistType::kTH1F, {axisPbPbTruthSelectionOnly}});
        setAxisLabels(registry.get<TH2>(HIST("h2_jet_pt_part_eventselection_bcBitsFirst_timeRangeOnly_truth"))->GetYaxis(), timeRangeOnlyLabels);
        setAxisLabels(registry.get<TH1>(HIST("h_mccollisions_eventselection_bcBitsFirst_timeRangeOnly_truth"))->GetXaxis(), timeRangeOnlyLabels);

        registry.add("h2_jet_pt_part_eventselection_bcBitsFirst_rofOnly_truth",
                     "part jet pT vs ROF-only truth selection;#it{p}_{T,jet}^{part} (GeV/#it{c});event selection;weighted counts",
                     {HistType::kTH2F, {jetPtAxis, axisPbPbTruthSelectionOnly}});
        registry.add("h_mccollisions_eventselection_bcBitsFirst_rofOnly_truth",
                     "weighted MC events vs ROF-only truth selection;event selection;weighted events",
                     {HistType::kTH1F, {axisPbPbTruthSelectionOnly}});
        setAxisLabels(registry.get<TH2>(HIST("h2_jet_pt_part_eventselection_bcBitsFirst_rofOnly_truth"))->GetYaxis(), rofOnlyLabels);
        setAxisLabels(registry.get<TH1>(HIST("h_mccollisions_eventselection_bcBitsFirst_rofOnly_truth"))->GetXaxis(), rofOnlyLabels);
      }
    }

    if (doprocessTruthMultFT0CvsRecoAmplitude) {
      AxisSpec truthMultFT0CAxis = {500, -0.5, 499.5, "truth multFT0C"};
      AxisSpec recoFT0CAmplitudeAxis = {1000, 0., 20000., "A_{FT0C}^{reco}"};
      registry.add("h2_mccollision_mult_ft0c_found_ft0_sum_amp_c",
                   "truth multFT0C vs reconstructed A_{FT0C};truth multFT0C;A_{FT0C}^{reco};weighted counts",
                   {HistType::kTH2F, {truthMultFT0CAxis, recoFT0CAmplitudeAxis}});
    }

    if (doprocessTruthMultFT0CvsRecoITS567) {
      AxisSpec truthMultFT0CAxis = {500, -0.5, 499.5, "truth multFT0C"};
      AxisSpec recoNITS567Axis = {1000, -0.5, 999.5, "N_{ITS567}^{reco}"};
      AxisSpec truthMultFT0CPositiveAxis = {2, -0.5, 1.5, "truth multFT0C > 0"};
      AxisSpec recoNITS567PositiveAxis = {2, -0.5, 1.5, "N_{ITS567}^{reco} > 0"};
      registry.add("h2_mccollision_mult_ft0c_reco_n_its567",
                   "truth multFT0C vs reconstructed N_{ITS567};truth multFT0C;N_{ITS567}^{reco};weighted counts",
                   {HistType::kTH2F, {truthMultFT0CAxis, recoNITS567Axis}});
      registry.add("h2_mccollision_mult_ft0c_positive_vs_reco_its567_positive",
                   "truth multFT0C activity vs reconstructed ITS567 activity;truth multFT0C > 0;N_{ITS567}^{reco} > 0;weighted counts",
                   {HistType::kTH2F, {truthMultFT0CPositiveAxis, recoNITS567PositiveAxis}});
      auto binaryActivity = registry.get<TH2>(HIST("h2_mccollision_mult_ft0c_positive_vs_reco_its567_positive"));
      binaryActivity->GetXaxis()->SetBinLabel(1, "multFT0C = 0");
      binaryActivity->GetXaxis()->SetBinLabel(2, "multFT0C > 0");
      binaryActivity->GetYaxis()->SetBinLabel(1, "N_{ITS567} = 0");
      binaryActivity->GetYaxis()->SetBinLabel(2, "N_{ITS567} > 0");
    }
  }

  template <typename TMcCollision>
  float computePtHat(TMcCollision const& mccollision)
  {
    float ptHardFromMc = ptHardCalcMethodSwitch;
    float storedPtHard = mccollision.ptHard();
    if (storedPtHard > BrokenPtHardSentinel && storedPtHard < ptHardCalcMethodSwitch) {
      ptHardFromMc = storedPtHard;
    }
    float weight = mccollision.weight();
    return ptHardFromMc < ptHardCalcMethodSwitch
             ? ptHardFromMc
             : simPtRef / std::pow(weight, 1.0f / pTHatExponent);
  }

  template <typename TBC>
  bool configureTruthRofParameters(TBC const& truthBC)
  {
    if (truthRofOffsetInBC >= 0 && truthRofLengthInBC > 0) {
      cachedTruthRofOffsetInBC = truthRofOffsetInBC;
      cachedTruthRofLengthInBC = truthRofLengthInBC;
      return true;
    }

    if (cachedTruthRofRun == truthBC.runNumber() && cachedTruthRofLengthInBC > 0) {
      return true;
    }

    auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", truthBC.timestamp());
    if (alppar == nullptr) {
      LOGF(fatal, "Could not retrieve ITS/Config/AlpideParam for sel8FullPbPb truth selections (run %d, timestamp %" PRIu64 ")", truthBC.runNumber(), static_cast<uint64_t>(truthBC.timestamp()));
      return false;
    }

    cachedTruthRofRun = truthBC.runNumber();
    cachedTruthRofOffsetInBC = truthRofOffsetInBC >= 0 ? truthRofOffsetInBC : alppar->roFrameBiasInBC;
    cachedTruthRofLengthInBC = truthRofLengthInBC > 0 ? truthRofLengthInBC : alppar->roFrameLengthInBC;
    if (cachedTruthRofLengthInBC <= 0) {
      LOGF(fatal, "Invalid ITS ROF length %" PRId64 " BC for sel8FullPbPb truth selections", static_cast<int64_t>(cachedTruthRofLengthInBC));
      return false;
    }
    LOGF(info, "sel8FullPbPb truth selections use ITS ROF offset %" PRId64 " and length %" PRId64 " BC for run %d", static_cast<int64_t>(cachedTruthRofOffsetInBC), static_cast<int64_t>(cachedTruthRofLengthInBC), cachedTruthRofRun);
    return true;
  }

  int64_t truthRofId(uint64_t globalBC) const
  {
    // Match EventSelectionModule: use the ITS ROF bias and length from DPLAlpideParam, with one orbit added to avoid a negative numerator.
    return (static_cast<int64_t>(globalBC) + o2::constants::lhc::LHCMaxBunches - cachedTruthRofOffsetInBC) / cachedTruthRofLengthInBC;
  }

  template <typename TMcCollision>
  void requireValidTruthFT0CActivity(TMcCollision const& mccollision) const
  {
    if (mccollision.multFT0C() < 0.0f) {
      LOGF(fatal, "sel8FullPbPb truth selection requires valid MC multFT0C, but multFT0C < 0 was found. Enable MC multiplicity information in the derived-data production.");
    }
  }

  template <typename TMcCollision>
  float truthFT0CActivity(TMcCollision const& mccollision) const
  {
    requireValidTruthFT0CActivity(mccollision);
    return mccollision.multFT0C();
  }

  template <typename TMcCollision, typename TAllMcCollisions>
  TruthPbPbSelections evaluateTruthPbPbSelections(TMcCollision const& mccollision, TAllMcCollisions const& allMcCollisions)
  {
    TruthPbPbSelections selections;
    requireValidTruthFT0CActivity(mccollision);
    auto truthBC = mccollision.template bc_as<aod::JBCs>();
    const uint64_t currentGlobalBC = truthBC.globalBC();
    if (currentGlobalBC == std::numeric_limits<uint64_t>::max() || !configureTruthRofParameters(truthBC)) {
      return selections;
    }

    const int64_t currentRof = truthRofId(currentGlobalBC);
    bool hasNarrowActivity = false;
    bool hasHighActivityInTimeRange = false;
    bool hasHighActivityInSameRof = false;
    bool hasCloseVzActivityInSameRof = false;

    // EventSelectionModule uses reconstructed foundGlobalBC, FT0C digit amplitudes and ITS layer 5-7 tracks.
    // Those are unavailable for hasColl=false MC collisions.
    // This truth-only approximation uses nominal MC BCs and multFT0C; its thresholds are intentionally configurable.
    // The scan is limited to the current input chunk because JMcCollisions has no truth-side found-BC or time-frame association.
    for (auto const& otherMcCollision : allMcCollisions) {
      if (otherMcCollision.globalIndex() == mccollision.globalIndex()) {
        continue;
      }

      auto otherTruthBC = otherMcCollision.template bc_as<aod::JBCs>();
      const uint64_t otherGlobalBC = otherTruthBC.globalBC();
      if (otherTruthBC.runNumber() != truthBC.runNumber() || otherGlobalBC == std::numeric_limits<uint64_t>::max()) {
        continue;
      }

      const int64_t deltaGlobalBC = static_cast<int64_t>(otherGlobalBC) - static_cast<int64_t>(currentGlobalBC);
      const float deltaTimeUs = static_cast<float>(deltaGlobalBC) * o2::constants::lhc::LHCBunchSpacingNS / 1000.0f;
      const bool inNarrowWindow = std::abs(deltaTimeUs) < truthTimeRangeNarrowUs;
      const bool inStandardTimeWindow = deltaTimeUs > truthTimeRangeStandardMinUs && deltaTimeUs < truthTimeRangeStandardMaxUs;
      const bool isSameRof = truthRofId(otherGlobalBC) == currentRof;
      if (!inNarrowWindow && !inStandardTimeWindow && !isSameRof) {
        continue;
      }

      const float otherActivity = truthFT0CActivity(otherMcCollision);
      const bool hasTrackPresenceProxy = otherActivity > truthFT0CActivityMinForTrackProxy;

      if (inNarrowWindow && hasTrackPresenceProxy) {
        hasNarrowActivity = true;
      }
      if (inStandardTimeWindow && otherActivity > truthFT0CActivityThresholdTimeRange) {
        hasHighActivityInTimeRange = true;
      }

      if (isSameRof) {
        if (otherActivity > truthFT0CActivityThresholdROF) {
          hasHighActivityInSameRof = true;
        }
        // The reco bit uses ITS layer 5-7 tracks and reconstructed PV z.
        // The truth approximation uses multFT0C and MC collision z instead.
        if (hasTrackPresenceProxy && std::abs(otherMcCollision.posZ() - mccollision.posZ()) < truthRofCloseVzMax) {
          hasCloseVzActivityInSameRof = true;
        }
      }
    }

    selections.valid = true;
    selections.noCollInTimeRangeStandard = !hasNarrowActivity && !hasHighActivityInTimeRange;
    selections.noCollInRofStandard = !hasHighActivityInSameRof && !hasCloseVzActivityInSameRof;
    return selections;
  }

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(TJets const& jet)
  {
    if (jetAreaFractionMin > ConfigSwitchLow) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentMinPt = (leadingConstituentPtMinMCP > ConfigSwitchLow);
    bool checkConstituentMaxPt = (leadingConstituentPtMaxMCP < ConfigSwitchHigh);
    bool checkConstituentPt = checkConstituentMinPt || checkConstituentMaxPt;

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<TTracks>()) {
        double constituentPt = constituent.pt();

        if (checkConstituentMinPt && constituentPt >= leadingConstituentPtMinMCP) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && constituentPt > leadingConstituentPtMaxMCP) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }
    return true;
  }

  void processTruthMultFT0CvsRecoAmplitude(aod::JetMcCollisions::iterator const& mccollision,
                                           soa::SmallGroups<JetCollisionsMCDWithParent> const& collisions,
                                           CollisionsWithEvSels const&,
                                           aod::FT0s const& ft0s)
  {
    if (skipMBGapEvents && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (mccollision.multFT0C() < 0.0f) {
      return;
    }
    // A single associated JCollision gives a one-to-one MC-collision to reconstructed-collision association.
    // The parent index restores the original EvSels relation used by EventSelectionModule.
    if (collisions.size() != 1) {
      return;
    }
    auto const& collision = collisions.begin();
    auto originalCollision = collision.template collision_as<CollisionsWithEvSels>();
    if (!originalCollision.has_foundFT0()) {
      return;
    }
    auto foundFT0 = ft0s.rawIteratorAt(originalCollision.foundFT0Id());
    registry.fill(HIST("h2_mccollision_mult_ft0c_found_ft0_sum_amp_c"), mccollision.multFT0C(), foundFT0.sumAmpC(), mccollision.weight());
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processTruthMultFT0CvsRecoAmplitude,
                 "truth multFT0C vs reconstructed foundFT0 A_FT0C for one-to-one MC/reco collision associations", false);

  void processTruthMultFT0CvsRecoITS567(aod::JetMcCollisions::iterator const& mccollision,
                                        soa::SmallGroups<JetCollisionsMCDWithParent> const& collisions,
                                        CollisionsWithEvSels const&,
                                        FullTracksIU const&)
  {
    if (skipMBGapEvents && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (mccollision.multFT0C() < 0.0f) {
      return;
    }
    // This matches EventSelectionModule::vTracksITS567perColl: only PV-contributor tracks with at least five ITS clusters are counted.
    if (collisions.size() != 1) {
      return;
    }
    auto const& collision = collisions.begin();
    auto originalCollision = collision.template collision_as<CollisionsWithEvSels>();
    auto const& collisionPvTracks = pvTracks.sliceBy(pvTracksPerCollision, originalCollision.globalIndex());
    int nITS567 = 0;
    for (const auto& track : collisionPvTracks) {
      if (track.itsNCls() >= MinITSClustersForOccupancy) {
        ++nITS567;
      }
    }

    const float weight = mccollision.weight();
    registry.fill(HIST("h2_mccollision_mult_ft0c_reco_n_its567"), mccollision.multFT0C(), nITS567, weight);
    registry.fill(HIST("h2_mccollision_mult_ft0c_positive_vs_reco_its567_positive"),
                  mccollision.multFT0C() > 0.0f, nITS567 > 0, weight);
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processTruthMultFT0CvsRecoITS567,
                 "truth multFT0C activity vs reconstructed PV-contributor ITS567 activity for one-to-one MC/reco collision associations", false);

  void processCrossSectionEfficiency(aod::JetMcCollisions::iterator const& mccollision,
                                     soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                     soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                     aod::JetParticles const&,
                                     aod::JBCs const&)
  {
    if (skipMBGapEvents && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }

    bool hasRecoColl = (collisions.size() >= 1);
    bool passesZvtxCutReco = false;
    if (hasRecoColl) {
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
        auto const& col = collisions.begin();
        passesZvtxCutReco = (std::abs(col.posZ()) <= vertexZCut);
      } else {
        for (auto const& col : collisions) {
          if (std::abs(col.posZ()) <= vertexZCut) {
            passesZvtxCutReco = true;
            break;
          }
        }
      }
    }
    bool noSplitPass = (acceptSplitCollisions == NonSplitOnly) ? (collisions.size() == 1) : true;

    bool passesTVX = false;
    bool passesNoTFB = false;
    bool passesNoITSROFB = false;
    bool passesNoSBP = false;
    bool passesNoCollInTimeRangeStandard = false;
    bool passesNoCollInRofStandard = false;
    if (hasRecoColl) {
      auto const& col = collisions.begin();
      auto evSel = col.eventSel();
      passesTVX = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selTVX)) != 0u;
      passesNoTFB = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder)) != 0u;
      passesNoITSROFB = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder)) != 0u;
      passesNoSBP = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selNoSameBunchPileup)) != 0u;
      passesNoCollInTimeRangeStandard = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selNoCollInTimeRangeStandard)) != 0u;
      passesNoCollInRofStandard = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selNoCollInRofStandard)) != 0u;
    }

    bool passesRct = applyRCT ? (mccollision.bc_as<aod::JBCs>().rct_raw() & rctMask) == 0 : true;
    std::vector<bool> pass = {true, passesRct, hasRecoColl, passesZvtxCutReco,
                              hasRecoColl && noSplitPass, passesTVX,
                              applyTFB ? passesNoTFB : true,
                              applyROFB ? passesNoITSROFB : true};
    if (applyNoCollInTimeRangeStandard) {
      pass.push_back(passesNoCollInTimeRangeStandard);
    }
    if (applyNoCollInRofStandard) {
      pass.push_back(passesNoCollInRofStandard);
    }
    if (!isSel8FullPbPb) {
      pass.push_back(applySBP ? passesNoSBP : true);
    }

    float weight = mccollision.weight();
    int sMax = 0;
    for (size_t s = 0; s < pass.size(); ++s) {
      if (!pass[s]) {
        break;
      }
      registry.fill(HIST("h_mccollisions_eventselection_collRecoFirst"), static_cast<double>(s + 1), weight);
      sMax = static_cast<int>(s + 1);
    }
    if (sMax == 0) {
      return;
    }

    float pTHat = computePtHat(mccollision);
    if (pTHat < pTHatAbsoluteMin) {
      return;
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) ||
          jet.pt() < jetPtMin || jet.pt() > pTHatMaxMCP * pTHat ||
          !isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      for (int s = 1; s <= sMax; ++s) {
        registry.fill(HIST("h2_jet_pt_part_eventselection_collRecoFirst"), jet.pt(), static_cast<double>(s), weight);
      }
    }
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processCrossSectionEfficiency,
                 "Cascade efficiency with collision-reco required first, BC bits read from the reco-coll EvSels bitmask", true);

  void processCrossSectionEfficiencyBcBitsFirst(aod::JetMcCollisions::iterator const& mccollision,
                                                soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                                soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                                aod::JetParticles const&,
                                                aod::JBCs const&,
                                                aod::JetMcCollisions const& allMcCollisions)
  {
    if (skipMBGapEvents && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }

    auto truthBC = mccollision.bc_as<aod::JBCs>();
    bool passesTVXTruth = truthBC.selection_bit(aod::evsel::kIsTriggerTVX);
    bool passesNoTFBTruth = truthBC.selection_bit(aod::evsel::kNoTimeFrameBorder);
    bool passesNoITSROFBTruth = truthBC.selection_bit(aod::evsel::kNoITSROFrameBorder);

    bool truthNoSBP = true;
    if (applySBP) {
      auto sameBC = allMcCollisions.sliceBy(mcCollsPerBC, mccollision.bcId());
      truthNoSBP = (sameBC.size() == 1);
    }

    bool hasRecoColl = (collisions.size() >= 1);
    bool passesZvtxCutReco = false;
    if (hasRecoColl) {
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
        auto const& col = collisions.begin();
        passesZvtxCutReco = (std::abs(col.posZ()) <= vertexZCut);
      } else {
        for (auto const& col : collisions) {
          if (std::abs(col.posZ()) <= vertexZCut) {
            passesZvtxCutReco = true;
            break;
          }
        }
      }
    }
    bool noSplitPass = (acceptSplitCollisions == NonSplitOnly) ? (collisions.size() == 1) : true;

    TruthPbPbSelections truthPbPbSelections;
    if (isSel8FullPbPb) {
      truthPbPbSelections = evaluateTruthPbPbSelections(mccollision, allMcCollisions);
    }

    bool passesRct = applyRCT ? (truthBC.rct_raw() & rctMask) == 0 : true;
    std::vector<bool> pass = {true, passesRct, passesTVXTruth,
                              applyTFB ? passesNoTFBTruth : true,
                              applyROFB ? passesNoITSROFBTruth : true};
    if (applyNoCollInTimeRangeStandard) {
      pass.push_back(truthPbPbSelections.valid && truthPbPbSelections.noCollInTimeRangeStandard);
    }
    if (applyNoCollInRofStandard) {
      pass.push_back(truthPbPbSelections.valid && truthPbPbSelections.noCollInRofStandard);
    }
    if (!isSel8FullPbPb) {
      pass.push_back(applySBP ? truthNoSBP : true);
    }
    pass.push_back(hasRecoColl);
    pass.push_back(passesZvtxCutReco);
    pass.push_back(hasRecoColl && noSplitPass);

    float weight = mccollision.weight();

    int sMaxTimeRangeOnly = 0;
    int sMaxRofOnly = 0;
    if (isSel8FullPbPb) {
      std::vector<bool> timeRangeOnlyPass = {true, passesRct, passesTVXTruth,
                                             passesNoTFBTruth, passesNoITSROFBTruth,
                                             truthPbPbSelections.valid && truthPbPbSelections.noCollInTimeRangeStandard,
                                             hasRecoColl, passesZvtxCutReco, hasRecoColl && noSplitPass};
      std::vector<bool> rofOnlyPass = {true, passesRct, passesTVXTruth,
                                       passesNoTFBTruth, passesNoITSROFBTruth,
                                       truthPbPbSelections.valid && truthPbPbSelections.noCollInRofStandard,
                                       hasRecoColl, passesZvtxCutReco, hasRecoColl && noSplitPass};
      for (size_t s = 0; s < timeRangeOnlyPass.size(); ++s) {
        if (!timeRangeOnlyPass[s]) {
          break;
        }
        registry.fill(HIST("h_mccollisions_eventselection_bcBitsFirst_timeRangeOnly_truth"), static_cast<double>(s + 1), weight);
        sMaxTimeRangeOnly = static_cast<int>(s + 1);
      }
      for (size_t s = 0; s < rofOnlyPass.size(); ++s) {
        if (!rofOnlyPass[s]) {
          break;
        }
        registry.fill(HIST("h_mccollisions_eventselection_bcBitsFirst_rofOnly_truth"), static_cast<double>(s + 1), weight);
        sMaxRofOnly = static_cast<int>(s + 1);
      }
    }

    int sMax = 0;
    for (size_t s = 0; s < pass.size(); ++s) {
      if (!pass[s]) {
        break;
      }
      registry.fill(HIST("h_mccollisions_eventselection_bcBitsFirst"), static_cast<double>(s + 1), weight);
      sMax = static_cast<int>(s + 1);
    }
    if (sMax == 0) {
      return;
    }

    float pTHat = computePtHat(mccollision);
    if (pTHat < pTHatAbsoluteMin) {
      return;
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) ||
          jet.pt() < jetPtMin || jet.pt() > pTHatMaxMCP * pTHat ||
          !isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      for (int s = 1; s <= sMax; ++s) {
        registry.fill(HIST("h2_jet_pt_part_eventselection_bcBitsFirst"), jet.pt(), static_cast<double>(s), weight);
      }
      if (isSel8FullPbPb) {
        for (int s = 1; s <= sMaxTimeRangeOnly; ++s) {
          registry.fill(HIST("h2_jet_pt_part_eventselection_bcBitsFirst_timeRangeOnly_truth"), jet.pt(), static_cast<double>(s), weight);
        }
        for (int s = 1; s <= sMaxRofOnly; ++s) {
          registry.fill(HIST("h2_jet_pt_part_eventselection_bcBitsFirst_rofOnly_truth"), jet.pt(), static_cast<double>(s), weight);
        }
      }
    }
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processCrossSectionEfficiencyBcBitsFirst,
                 "Cascade efficiency with BC bits evaluated on MC truth BC before collision reco", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetCrossSectionEfficiency>(cfgc)};
}
