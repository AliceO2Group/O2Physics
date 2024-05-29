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

#ifndef PWGHF_UTILS_UTILSEVSELHF_H_
#define PWGHF_UTILS_UTILSEVSELHF_H_

#include <memory> // std::shared_ptr
#include <string> // std::string

#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include "PWGHF/Core/CentralityEstimation.h"

namespace o2::hf_evsel
{
// event rejection types
enum EventRejection {
  None = 0,
  Centrality,
  Trigger,
  TimeFrameBorderCut,
  IsGoodZvtxFT0vsPV,
  NoSameBunchPileup,
  NumTracksInTimeRange,
  NContrib,
  Chi2,
  PositionZ,
  NEventRejection
};

struct HfEventSelection : o2::framework::ConfigurableGroup {
  std::string prefix = "hfEvSel"; // JSON group name
  // event selection parameters (in chronological order of application)
  o2::framework::Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality"};
  o2::framework::Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality"};
  o2::framework::Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "Apply the sel8 event selection"};
  o2::framework::Configurable<int> triggerClass{"triggerClass", -1, "Trigger class different from sel8 (e.g. kINT7 for Run2) used only if useSel8Trigger is false"};
  o2::framework::Configurable<bool> useTimeFrameBorderCut{"useTimeFrameBorderCut", true, "Apply TF border cut"};
  o2::framework::Configurable<bool> useIsGoodZvtxFT0vsPV{"useIsGoodZvtxFT0vsPV", false, "Check consistency between PVz from central barrel with that from FT0 timing"};
  o2::framework::Configurable<bool> useNoSameBunchPileup{"useNoSameBunchPileup", false, "Exclude collisions in bunches with more than 1 reco. PV"}; // POTENTIALLY BAD FOR BEAUTY ANALYSES
  o2::framework::Configurable<bool> useNumTracksInTimeRange{"useNumTracksInTimeRange", false, "Apply occupancy selection (num. ITS tracks with at least 5 clusters in +-100us from current collision)"};
  o2::framework::Configurable<int> numTracksInTimeRangeMin{"numTracksInTimeRangeMin", 0, "Minimum occupancy"};
  o2::framework::Configurable<int> numTracksInTimeRangeMax{"numTracksInTimeRangeMax", 1000000, "Maximum occupancy"};
  o2::framework::Configurable<int> nPvContributorsMin{"nPvContributorsMin", 0, "Minimum number of PV contributors"};
  o2::framework::Configurable<float> chi2PvMax{"chi2PvMax", -1.f, "Maximum PV chi2"};
  o2::framework::Configurable<float> zPvPosMin{"zPvPosMin", -10.f, "Minimum PV posZ (cm)"};
  o2::framework::Configurable<float> zPvPosMax{"zPvPosMax", 10.f, "Maximum PV posZ (cm)"};

  // histogram names
  static constexpr char nameHistCollisions[] = "hCollisions";
  static constexpr char nameHistPosZBeforeEvSel[] = "hPosZBeforeEvSel";
  static constexpr char nameHistPosZAfterEvSel[] = "hPosZAfterEvSel";
  static constexpr char nameHistPosXAfterEvSel[] = "hPosXAfterEvSel";
  static constexpr char nameHistPosYAfterEvSel[] = "hPosYAfterEvSel";
  static constexpr char nameHistNumPvContributorsAfterSel[] = "hNumPvContributorsAfterSel";

  std::shared_ptr<TH1> hCollisions, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel;

  /// \brief Adds collision monitoring histograms in the histogram registry.
  /// \param registry reference to the histogram registry
  void addHistograms(o2::framework::HistogramRegistry& registry)
  {
    o2::framework::AxisSpec axisEvents = {EventRejection::NEventRejection, -0.5f, +EventRejection::NEventRejection - 0.5f, ""};
    hCollisions = registry.add<TH1>(nameHistCollisions, "HF event counter;;entries", {o2::framework::HistType::kTH1D, {axisEvents}});
    hPosZBeforeEvSel = registry.add<TH1>(nameHistPosZBeforeEvSel, "all events;#it{z}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{400, -20., 20.}}});
    hPosZAfterEvSel = registry.add<TH1>(nameHistPosZAfterEvSel, "selected events;#it{z}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{400, -20., 20.}}});
    hPosXAfterEvSel = registry.add<TH1>(nameHistPosXAfterEvSel, "selected events;#it{x}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{200, -0.5, 0.5}}});
    hPosYAfterEvSel = registry.add<TH1>(nameHistPosYAfterEvSel, "selected events;#it{y}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{200, -0.5, 0.5}}});
    hNumPvContributorsAfterSel = registry.add<TH1>(nameHistNumPvContributorsAfterSel, "selected events;#it{y}_{prim. vtx.} (cm);entries", {o2::framework::HistType::kTH1D, {{500, -0.5, 499.5}}});
    // Puts labels on the collision monitoring histogram.
    hCollisions->SetTitle("HF event counter;;# of accepted collisions");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::None + 1, "All");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::Centrality + 1, "Centrality");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::Trigger + 1, "Trigger");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::TimeFrameBorderCut + 1, "TF border");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::IsGoodZvtxFT0vsPV + 1, "PV #it{z} consistency FT0 timing");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::NoSameBunchPileup + 1, "No same-bunch pile-up"); // POTENTIALLY BAD FOR BEAUTY ANALYSES
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::NumTracksInTimeRange + 1, "Occupancy");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::NContrib + 1, "# of PV contributors");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::Chi2 + 1, "PV #it{#chi}^{2}");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::PositionZ + 1, "PV #it{z}");
  }

  /// \brief Applies event selection.
  /// \tparam useEvSel use information from the EvSel table
  /// \tparam centEstimator centrality estimator
  /// \param collision collision to test against the selection criteria
  /// \param centrality collision centrality variable to be set in this function
  /// \return bitmask with the event selection criteria not satisfied by the collision
  template <bool useEvSel, o2::hf_centrality::CentralityEstimator centEstimator, typename Coll>
  uint16_t getHfCollisionRejectionMask(const Coll& collision, float& centrality)
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
      /// time frame border cut
      if (useTimeFrameBorderCut && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        SETBIT(rejectionMask, EventRejection::TimeFrameBorderCut);
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

    return rejectionMask;
  }

  /// \brief Fills histograms for monitoring event selections satisfied by the collision.
  /// \param collision analysed collision
  /// \param rejectionMask bitmask storing the info about which ev. selections are not satisfied by the collision
  template <typename Coll>
  void fillHistograms(Coll const& collision, const uint16_t rejectionMask)
  {
    hCollisions->Fill(EventRejection::None);
    const float posZ = collision.posZ();
    hPosZBeforeEvSel->Fill(posZ);

    for (size_t reason = 1; reason < EventRejection::NEventRejection; reason++) {
      if (TESTBIT(rejectionMask, reason)) {
        return;
      }
      hCollisions->Fill(reason);
    }

    hPosXAfterEvSel->Fill(collision.posX());
    hPosYAfterEvSel->Fill(collision.posY());
    hPosZAfterEvSel->Fill(posZ);
    hNumPvContributorsAfterSel->Fill(collision.numContrib());
  }
};
} // namespace o2::hf_evsel

#endif // PWGHF_UTILS_UTILSEVSELHF_H_
