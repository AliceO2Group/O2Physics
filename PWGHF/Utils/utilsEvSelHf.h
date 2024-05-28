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

#include <string>

#include "Framework/Configurable.h"
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

o2::framework::AxisSpec axisEvents = {EventRejection::NEventRejection, -0.5f, +EventRejection::NEventRejection - 0.5f, ""};

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

  /// \brief Function to put labels on collision monitoring histogram
  /// \param hCollisions is the histogram
  template <typename Histo>
  void setLabelHistoEvSel(Histo& hCollisions)
  {
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

  /// \brief Function to apply event selections in HF analyses
  /// \tparam applyEvSel use information from the EvSel table
  /// \tparam centEstimator centrality estimator
  /// \param collision collision that has to satisfy the selection criteria
  /// \param centrality collision centrality to be set in this function
  /// \return bitmask with the event selections not satisfied by the collision
  template <bool applyEvSel, o2::hf_centrality::CentralityEstimator centEstimator, typename Coll>
  uint16_t getHfCollisionRejectionMask(const Coll& collision, float& centrality)
  {
    uint16_t statusCollision{0}; // 16 bits, in case new ev. selections will be added

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
        SETBIT(statusCollision, EventRejection::Centrality);
      }
    }

    if constexpr (applyEvSel) {
      /// trigger condition
      if ((useSel8Trigger && !collision.sel8()) || (!useSel8Trigger && triggerClass > -1 && !collision.alias_bit(triggerClass))) {
        SETBIT(statusCollision, EventRejection::Trigger);
      }
      /// time frame border cut
      if (useTimeFrameBorderCut && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        SETBIT(statusCollision, EventRejection::TimeFrameBorderCut);
      }
      /// PVz consistency tracking - FT0 timing
      if (useIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        SETBIT(statusCollision, EventRejection::IsGoodZvtxFT0vsPV);
      }
      /// remove collisions in bunches with more than 1 reco collision
      /// POTENTIALLY BAD FOR BEAUTY ANALYSES
      if (useNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        SETBIT(statusCollision, EventRejection::NoSameBunchPileup);
      }
      /// occupancy estimator (ITS tracks with at least 5 clusters in +-10us from current collision)
      if (useNumTracksInTimeRange) {
        const int numTracksInTimeRange = collision.trackOccupancyInTimeRange();
        if (numTracksInTimeRange < numTracksInTimeRangeMin || numTracksInTimeRange > numTracksInTimeRangeMax) {
          SETBIT(statusCollision, EventRejection::NumTracksInTimeRange);
        }
      }
    }

    /// number of PV contributors
    if (collision.numContrib() < nPvContributorsMin) {
      SETBIT(statusCollision, EventRejection::NContrib);
    }

    /// max PV chi2
    if (chi2PvMax > 0. && collision.chi2() > chi2PvMax) {
      SETBIT(statusCollision, EventRejection::Chi2);
    }

    /// primary vertex z
    if (collision.posZ() < zPvPosMin || collision.posZ() > zPvPosMax) {
      SETBIT(statusCollision, EventRejection::PositionZ);
    }

    return statusCollision;
  }

  /// \brief function to monitor the event selection satisfied by collisions used for HF analyses
  /// \param collision is the analysed collision
  /// \param rejectionMask is the bitmask storing the info about which ev. selections are not satisfied by the collision
  /// \param hCollisions is a histogram to keep track of the satisfied event selections
  /// \param hPosZBeforeEvSel is a histogram for the PV position Z for all analysed collisions
  /// \param hPosZAfterEvSel is a histogram for the PV position Z only for collisions satisfying the event selections
  /// \param hPosXAfterEvSel is a histogram for the PV position X only for collisions satisfying the event selections
  /// \param hPosYAfterEvSel is a histogram for the PV position Y only for collisions satisfying the event selections
  /// \param hNumContributors is a histogram for the number of PV contributors only for collisions satisfying the event selections
  template <typename Coll, typename Hist>
  void monitorCollision(Coll const& collision, const uint16_t rejectionMask, Hist& hCollisions, Hist& hPosZBeforeEvSel, Hist& hPosZAfterEvSel, Hist& hPosXAfterEvSel, Hist& hPosYAfterEvSel, Hist& hNumContributors)
  {
    hCollisions->Fill(EventRejection::None); // all collisions
    const float posZ = collision.posZ();
    hPosZBeforeEvSel->Fill(posZ);

    for (size_t reason = 1; reason < EventRejection::NEventRejection; reason++) {
      if (TESTBIT(rejectionMask, reason)) {
        return;
      }
      hCollisions->Fill(reason); // Centrality ok
    }

    hPosXAfterEvSel->Fill(collision.posX());
    hPosYAfterEvSel->Fill(collision.posY());
    hPosZAfterEvSel->Fill(posZ);
    hNumContributors->Fill(collision.numContrib());
  }
};
} // namespace o2::hf_evsel

#endif // PWGHF_UTILS_UTILSEVSELHF_H_
