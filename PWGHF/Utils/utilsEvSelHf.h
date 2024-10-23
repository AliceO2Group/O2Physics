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

/// \file utilsEvSelHf.h
/// \brief Event selection utilities for HF analyses
/// \author Mattia Faggin <mfaggin@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, CERN

#ifndef PWGHF_UTILS_UTILSEVSELHF_H_
#define PWGHF_UTILS_UTILSEVSELHF_H_

#include <memory> // std::shared_ptr
#include <string> // std::string

#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/OutputObjHeader.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"
#include "PWGHF/Core/CentralityEstimation.h"

namespace o2::hf_evsel
{
// event rejection types
enum EventRejection {
  None = 0,
  SoftwareTrigger,
  Centrality,
  Trigger,
  TvxTrigger,
  TimeFrameBorderCut,
  ItsRofBorderCut,
  IsGoodZvtxFT0vsPV,
  NoSameBunchPileup,
  NumTracksInTimeRange,
  NContrib,
  Chi2,
  PositionZ,
  NEventRejection
};

o2::framework::AxisSpec axisEvents = {EventRejection::NEventRejection, -0.5f, +EventRejection::NEventRejection - 0.5f, ""};

/// \brief Function to put labels on monitoring histogram
/// \param hRejection monitoring histogram
/// \param softwareTriggerLabel bin label for software trigger rejection
template <typename Histo>
void setEventRejectionLabels(Histo& hRejection, std::string softwareTriggerLabel = "")
{
  // Puts labels on the collision monitoring histogram.
  hRejection->GetXaxis()->SetBinLabel(EventRejection::None + 1, "All");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::SoftwareTrigger + 1, softwareTriggerLabel.data());
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Centrality + 1, "Centrality");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Trigger + 1, "Trigger");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::TvxTrigger + 1, "TVX Trigger");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::TimeFrameBorderCut + 1, "TF border");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::ItsRofBorderCut + 1, "ITS ROF border");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::IsGoodZvtxFT0vsPV + 1, "PV #it{z} consistency FT0 timing");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::NoSameBunchPileup + 1, "No same-bunch pile-up"); // POTENTIALLY BAD FOR BEAUTY ANALYSES
  hRejection->GetXaxis()->SetBinLabel(EventRejection::NumTracksInTimeRange + 1, "Occupancy");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::NContrib + 1, "# of PV contributors");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::Chi2 + 1, "PV #it{#chi}^{2}");
  hRejection->GetXaxis()->SetBinLabel(EventRejection::PositionZ + 1, "PV #it{z}");
}

struct HfEventSelection : o2::framework::ConfigurableGroup {
  std::string prefix = "hfEvSel"; // JSON group name
  // event selection parameters (in chronological order of application)
  o2::framework::Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality"};
  o2::framework::Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality"};
  o2::framework::Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "Apply the sel8 event selection"};
  o2::framework::Configurable<int> triggerClass{"triggerClass", -1, "Trigger class different from sel8 (e.g. kINT7 for Run2) used only if useSel8Trigger is false"};
  o2::framework::Configurable<bool> useTvxTrigger{"useTvxTrigger", true, "Apply TVX trigger sel"};
  o2::framework::Configurable<bool> useTimeFrameBorderCut{"useTimeFrameBorderCut", true, "Apply TF border cut"};
  o2::framework::Configurable<bool> useItsRofBorderCut{"useItsRofBorderCut", true, "Apply ITS ROF border cut"};
  o2::framework::Configurable<bool> useIsGoodZvtxFT0vsPV{"useIsGoodZvtxFT0vsPV", false, "Check consistency between PVz from central barrel with that from FT0 timing"};
  o2::framework::Configurable<bool> useNoSameBunchPileup{"useNoSameBunchPileup", false, "Exclude collisions in bunches with more than 1 reco. PV"}; // POTENTIALLY BAD FOR BEAUTY ANALYSES
  o2::framework::Configurable<bool> useNumTracksInTimeRange{"useNumTracksInTimeRange", false, "Apply occupancy selection (num. ITS tracks with at least 5 clusters in +-100us from current collision)"};
  o2::framework::Configurable<int> numTracksInTimeRangeMin{"numTracksInTimeRangeMin", 0, "Minimum occupancy"};
  o2::framework::Configurable<int> numTracksInTimeRangeMax{"numTracksInTimeRangeMax", 1000000, "Maximum occupancy"};
  o2::framework::Configurable<int> nPvContributorsMin{"nPvContributorsMin", 0, "Minimum number of PV contributors"};
  o2::framework::Configurable<float> chi2PvMax{"chi2PvMax", -1.f, "Maximum PV chi2"};
  o2::framework::Configurable<float> zPvPosMin{"zPvPosMin", -10.f, "Minimum PV posZ (cm)"};
  o2::framework::Configurable<float> zPvPosMax{"zPvPosMax", 10.f, "Maximum PV posZ (cm)"};
  o2::framework::Configurable<std::string> softwareTrigger{"softwareTrigger", "", "Label of software trigger. Multiple triggers can be selected dividing them by a comma. Set None if you want bcs that are not selected by any trigger"};
  o2::framework::Configurable<uint64_t> bcMarginForSoftwareTrigger{"bcMarginForSoftwareTrigger", 100, "Number of BCs of margin for software triggers"};
  o2::framework::Configurable<std::string> ccdbPathSoftwareTrigger{"ccdbPathSoftwareTrigger", "Users/m/mpuccio/EventFiltering/OTS/", "ccdb path for ZORRO objects"};

  // histogram names
  static constexpr char nameHistCollisions[] = "hCollisions";
  static constexpr char nameHistSelCollisionsCent[] = "hSelCollisionsCent";
  static constexpr char nameHistPosZBeforeEvSel[] = "hPosZBeforeEvSel";
  static constexpr char nameHistPosZAfterEvSel[] = "hPosZAfterEvSel";
  static constexpr char nameHistPosXAfterEvSel[] = "hPosXAfterEvSel";
  static constexpr char nameHistPosYAfterEvSel[] = "hPosYAfterEvSel";
  static constexpr char nameHistNumPvContributorsAfterSel[] = "hNumPvContributorsAfterSel";

  std::shared_ptr<TH1> hCollisions, hSelCollisionsCent, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel;

  // util to retrieve trigger mask in case of software triggers
  Zorro zorro;
  o2::framework::OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  int currentRun{-1};

  /// \brief Adds collision monitoring histograms in the histogram registry.
  /// \param registry reference to the histogram registry
  void addHistograms(o2::framework::HistogramRegistry& registry)
  {
    hCollisions = registry.add<TH1>(nameHistCollisions, "HF event counter;;# of accepted collisions", {o2::framework::HistType::kTH1D, {axisEvents}});
    hSelCollisionsCent = registry.add<TH1>(nameHistSelCollisionsCent, "HF event counter;T0M;# of accepted collisions", {o2::framework::HistType::kTH1D, {{100, 0., 100.}}});
    hPosZBeforeEvSel = registry.add<TH1>(nameHistPosZBeforeEvSel, "all events;#it{z}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{400, -20., 20.}}});
    hPosZAfterEvSel = registry.add<TH1>(nameHistPosZAfterEvSel, "selected events;#it{z}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{400, -20., 20.}}});
    hPosXAfterEvSel = registry.add<TH1>(nameHistPosXAfterEvSel, "selected events;#it{x}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{200, -0.5, 0.5}}});
    hPosYAfterEvSel = registry.add<TH1>(nameHistPosYAfterEvSel, "selected events;#it{y}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{200, -0.5, 0.5}}});
    hNumPvContributorsAfterSel = registry.add<TH1>(nameHistNumPvContributorsAfterSel, "selected events;#it{y}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{500, -0.5, 499.5}}});
    setEventRejectionLabels(hCollisions, softwareTrigger);

    // we initialise the summary object
    if (softwareTrigger.value != "") {
      zorroSummary.setObject(zorro.getZorroSummary());
    }
  }

  /// \brief Applies event selection.
  /// \tparam useEvSel use information from the EvSel table
  /// \tparam centEstimator centrality estimator
  /// \param collision collision to test against the selection criteria
  /// \param centrality collision centrality variable to be set in this function
  /// \param ccdb ccdb service needed to retrieve the needed info for zorro
  /// \param registry reference to the histogram registry needed for zorro
  /// \return bitmask with the event selection criteria not satisfied by the collision
  template <bool useEvSel, o2::hf_centrality::CentralityEstimator centEstimator, typename BCs, typename Coll>
  uint16_t getHfCollisionRejectionMask(const Coll& collision, float& centrality, o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, o2::framework::HistogramRegistry& registry)
  {
    uint16_t rejectionMask{0}; // 16 bits, in case new ev. selections will be added

    if constexpr (centEstimator != o2::hf_centrality::CentralityEstimator::None) {
      if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0A) {
        centrality = collision.centFT0A();
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0C) {
        centrality = collision.centFT0C();
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
        centrality = collision.centFT0M();
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FV0A) {
        centrality = collision.centFV0A();
      } else {
        LOGP(fatal, "Unsupported centrality estimator!");
      }
      if (centrality < centralityMin || centrality > centralityMax) {
        SETBIT(rejectionMask, EventRejection::Centrality);
      }
    }

    if constexpr (useEvSel) {
      /// trigger condition
      if ((useSel8Trigger && !collision.sel8()) || (!useSel8Trigger && triggerClass > -1 && !collision.alias_bit(triggerClass))) {
        SETBIT(rejectionMask, EventRejection::Trigger);
      }
      /// TVX trigger selection
      if (useTvxTrigger && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        SETBIT(rejectionMask, EventRejection::TvxTrigger);
      }
      /// time frame border cut
      if (useTimeFrameBorderCut && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        SETBIT(rejectionMask, EventRejection::TimeFrameBorderCut);
      }
      /// ITS rof border cut
      if (useItsRofBorderCut && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        SETBIT(rejectionMask, EventRejection::ItsRofBorderCut);
      }
      /// PVz consistency tracking - FT0 timing
      if (useIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        SETBIT(rejectionMask, EventRejection::IsGoodZvtxFT0vsPV);
      }
      /// remove collisions in bunches with more than 1 reco collision
      /// POTENTIALLY BAD FOR BEAUTY ANALYSES
      if (useNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        SETBIT(rejectionMask, EventRejection::NoSameBunchPileup);
      }
      /// occupancy estimator (ITS tracks with at least 5 clusters in +-10us from current collision)
      if (useNumTracksInTimeRange) {
        const int numTracksInTimeRange = collision.trackOccupancyInTimeRange();
        if (numTracksInTimeRange < numTracksInTimeRangeMin || numTracksInTimeRange > numTracksInTimeRangeMax) {
          SETBIT(rejectionMask, EventRejection::NumTracksInTimeRange);
        }
      }
    }

    /// number of PV contributors
    if (collision.numContrib() < nPvContributorsMin) {
      SETBIT(rejectionMask, EventRejection::NContrib);
    }

    /// max PV chi2
    if (chi2PvMax > 0. && collision.chi2() > chi2PvMax) {
      SETBIT(rejectionMask, EventRejection::Chi2);
    }

    /// primary vertex z
    if (collision.posZ() < zPvPosMin || collision.posZ() > zPvPosMax) {
      SETBIT(rejectionMask, EventRejection::PositionZ);
    }

    if (softwareTrigger.value != "") {
      // we might have to update it from CCDB
      auto bc = collision.template bc_as<BCs>();

      int runNumber = bc.runNumber();
      if (runNumber != currentRun) { // We might need to update Zorro from CCDB if the run number changes
        zorro.setCCDBpath(ccdbPathSoftwareTrigger);
        zorro.setBCtolerance(bcMarginForSoftwareTrigger);
        zorro.initCCDB(ccdb.service, runNumber, bc.timestamp(), softwareTrigger.value);
        currentRun = runNumber;
      }
      zorro.populateHistRegistry(registry, runNumber);

      if (softwareTrigger.value != "None") {
        if (!zorro.isSelected(bc.globalBC(), bcMarginForSoftwareTrigger)) { /// Just let Zorro do the accounting
          SETBIT(rejectionMask, EventRejection::SoftwareTrigger);
        }
      } else {
        if (!zorro.isNotSelectedByAny(bc.globalBC(), bcMarginForSoftwareTrigger)) { /// Just let Zorro do the accounting of not selected BCs
          SETBIT(rejectionMask, EventRejection::SoftwareTrigger);
        }
      }
    }

    return rejectionMask;
  }

  /// \brief Fills histograms for monitoring event selections satisfied by the collision.
  /// \param collision analysed collision
  /// \param rejectionMask bitmask storing the info about which ev. selections are not satisfied by the collision
  template <typename Coll>
  void fillHistograms(Coll const& collision, const uint16_t rejectionMask, float& centrality)
  {
    hCollisions->Fill(EventRejection::None);
    const float posZ = collision.posZ();
    hPosZBeforeEvSel->Fill(posZ);

    for (std::size_t reason = 1; reason < EventRejection::NEventRejection; reason++) {
      if (TESTBIT(rejectionMask, reason)) {
        return;
      }
      hCollisions->Fill(reason);
    }

    hPosXAfterEvSel->Fill(collision.posX());
    hPosYAfterEvSel->Fill(collision.posY());
    hPosZAfterEvSel->Fill(posZ);
    hNumPvContributorsAfterSel->Fill(collision.numContrib());
    hSelCollisionsCent->Fill(centrality);
  }
};

struct HfEventSelectionMc {
  // event selection parameters (in chronological order of application)
  bool useSel8Trigger{false};       // Apply the Sel8 selection
  bool useTvxTrigger{false};        // Apply the TVX trigger
  bool useTimeFrameBorderCut{true}; // Apply TF border cut
  bool useItsRofBorderCut{false};   // Apply the ITS RO frame border cut
  float zPvPosMin{-1000.f};         // Minimum PV posZ (cm)
  float zPvPosMax{1000.f};          // Maximum PV posZ (cm)
  float centralityMin{0.f};         // Minimum centrality
  float centralityMax{100.f};       // Maximum centrality

  // histogram names
  static constexpr char nameHistParticles[] = "hParticles";
  std::shared_ptr<TH1> hParticles;

  /// \brief Adds collision monitoring histograms in the histogram registry.
  /// \param registry reference to the histogram registry
  void addHistograms(o2::framework::HistogramRegistry& registry)
  {
    hParticles = registry.add<TH1>(nameHistParticles, "HF particle counter;;# of accepted particles", {o2::framework::HistType::kTH1D, {axisEvents}});
    // Puts labels on the collision monitoring histogram.
    setEventRejectionLabels(hParticles);
  }

  void configureFromDevice(const o2::framework::DeviceSpec& device)
  {
    for (const auto& option : device.options) {
      if (option.name.compare("hfEvSel.useSel8Trigger") == 0) {
        useSel8Trigger = option.defaultValue.get<bool>();
      } else if (option.name.compare("hfEvSel.useTvxTrigger") == 0) {
        useTvxTrigger = option.defaultValue.get<bool>();
      } else if (option.name.compare("hfEvSel.useTimeFrameBorderCut") == 0) {
        useTimeFrameBorderCut = option.defaultValue.get<bool>();
      } else if (option.name.compare("hfEvSel.useItsRofBorderCut") == 0) {
        useItsRofBorderCut = option.defaultValue.get<bool>();
      } else if (option.name.compare("hfEvSel.zPvPosMin") == 0) {
        zPvPosMin = option.defaultValue.get<float>();
      } else if (option.name.compare("hfEvSel.zPvPosMax") == 0) {
        zPvPosMax = option.defaultValue.get<float>();
      } else if (option.name.compare("hfEvSel.centralityMin") == 0) {
        centralityMin = option.defaultValue.get<float>();
      } else if (option.name.compare("hfEvSel.centralityMax") == 0) {
        centralityMax = option.defaultValue.get<float>();
      }
    }
  }

  /// \brief Function to apply event selections to generated MC collisions
  /// \param mcCollision MC collision to test against the selection criteria
  /// \param collSlice collection of reconstructed collisions
  /// \param centrality centrality variable to be set in this function
  /// \return a bitmask with the event selections not satisfied by the analysed collision
  template <typename TBc, o2::hf_centrality::CentralityEstimator centEstimator, typename CCs, typename TMcColl>
  uint16_t getHfMcCollisionRejectionMask(TMcColl const& mcCollision, CCs const& collSlice, float& centrality)
  {
    uint16_t rejectionMask{0};
    float zPv = mcCollision.posZ();
    auto bc = mcCollision.template bc_as<TBc>();

    if constexpr (centEstimator != o2::hf_centrality::CentralityEstimator::None) {
      float multiplicity{0.f};
      for (const auto& collision : collSlice) {
        float collCent{0.f};
        float collMult{0.f};
        if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0A) {
          collCent = collision.centFT0A();
        } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0C) {
          collCent = collision.centFT0C();
        } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
          collCent = collision.centFT0M();
        } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FV0A) {
          collCent = collision.centFV0A();
        } else {
          LOGP(fatal, "Unsupported centrality estimator!");
        }
        collMult = collision.numContrib();
        if (collMult > multiplicity) {
          centrality = collCent;
          multiplicity = collMult;
        }
      }
      /// centrality selection
      if (centrality < centralityMin || centrality > centralityMax) {
        SETBIT(rejectionMask, EventRejection::Centrality);
      }
    }
    /// Sel8 trigger selection
    if (useSel8Trigger && (!bc.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
      SETBIT(rejectionMask, EventRejection::Trigger);
    }
    /// TVX trigger selection
    if (useTvxTrigger && !bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      SETBIT(rejectionMask, EventRejection::TvxTrigger);
    }
    /// time frame border cut
    if (useTimeFrameBorderCut && !bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      SETBIT(rejectionMask, EventRejection::TimeFrameBorderCut);
    }
    /// ITS RO frame border cut
    if (useItsRofBorderCut && !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      SETBIT(rejectionMask, EventRejection::ItsRofBorderCut);
    }
    /// primary vertex z
    if (zPv < zPvPosMin || zPv > zPvPosMax) {
      SETBIT(rejectionMask, EventRejection::PositionZ);
    }

    return rejectionMask;
  }

  /// \brief Fills histogram for monitoring event selections satisfied by the collision.
  /// \param rejectionMask bitmask storing the info about which ev. selections are not satisfied by the collision
  void fillHistograms(const uint16_t rejectionMask)
  {
    hParticles->Fill(EventRejection::None);

    for (std::size_t reason = 1; reason < EventRejection::NEventRejection; reason++) {
      if (TESTBIT(rejectionMask, reason)) {
        return;
      }
      hParticles->Fill(reason);
    }
  }
};
} // namespace o2::hf_evsel

#endif // PWGHF_UTILS_UTILSEVSELHF_H_
