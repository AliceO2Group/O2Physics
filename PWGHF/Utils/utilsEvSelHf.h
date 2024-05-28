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
/// \brief Utility set the event selections for HF analyses
/// \author Mattia Faggin <mfaggin@cern.ch>, CERN

#ifndef PWGHF_UTILS_UTILSEVSELHF_H_
#define PWGHF_UTILS_UTILSEVSELHF_H_

#include "Framework/Configurable.h"
#include "Framework/HistogramSpec.h"

namespace o2::framework::hf_evsel
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

o2::framework::AxisSpec axisEvents = {EventRejection::NEventRejection, -0.5f, static_cast<float>(EventRejection::NEventRejection) - 0.5f, ""};

/// @brief Function to put labels on collision monitoring histogram
/// \param hCollisions is the histogram
template <typename Histo>
void setLabelHistoEvSel(Histo& hCollisions)
{
  hCollisions->SetTitle("HF event counter;;accepted collisions");
  hCollisions->GetXaxis()->SetBinLabel(EventRejection::None + 1, "All collisions");
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
/// \tparam centEstimator select the centrality estimator
/// \param collision collision that has to satisfy the selection criteria
/// \param centrality collision centrality to be initialised in this function
/// \param centralityMin minimum centrality accepted
/// \param centralityMax maximum centrality accepted
/// \param useSel8Trigger switch to activate the sel8() event selection
/// \param triggerClass trigger class different from sel8 (e.g. kINT7 for Run2) used only if useSel8Trigger is false
/// \param useTimeFrameBorderCut switch to activate the time frame border cut
/// \param zPvPosMin minimum primary-vertex z
/// \param zPvPosMax maximum primary-vertex z
/// \param nPvContributorsMin minimum number of PV contributors
/// \param chi2PvMax maximum PV chi2
/// \return a bitmask with the event selections not satisfied by the analysed collision
template <bool applyEvSel, o2::aod::hf_collision_centrality::CentralityEstimator centEstimator, typename Coll>
uint16_t getHfCollisionRejectionMask(const Coll& collision, float& centrality, float centralityMin, float centralityMax,
                                     bool useSel8Trigger,
                                     int triggerClass,
                                     bool useTimeFrameBorderCut,
                                     float zPvPosMin, float zPvPosMax,
                                     int nPvContributorsMin,
                                     float chi2PvMax,
                                     bool useIsGoodZvtxFT0vsPV, bool useNoSameBunchPileup,
                                     bool useNumTracksInTimeRange, int numTracksInTimeRangeMin, int numTracksInTimeRangeMax)
{

  uint16_t statusCollision{0}; // 16 bits, in case new ev. selections will be added

  if constexpr (centEstimator != o2::aod::hf_collision_centrality::CentralityEstimator::None) {
    if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FT0A) {
      centrality = collision.centFT0A();
    } else if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FT0C) {
      centrality = collision.centFT0C();
    } else if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FT0M) {
      centrality = collision.centFT0M();
    } else if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FV0A) {
      centrality = collision.centFV0A();
    } else {
      LOGP(fatal, "Centrality estimator different from FT0A, FT0C, FT0M, and FV0A, fix it!");
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

  /// primary vertex z
  if (collision.posZ() < zPvPosMin || collision.posZ() > zPvPosMax) {
    SETBIT(statusCollision, EventRejection::PositionZ);
  }

  /// number of PV contributors
  if (collision.numContrib() < nPvContributorsMin) {
    SETBIT(statusCollision, EventRejection::NContrib);
  }

  /// max PV chi2
  if (chi2PvMax > 0. && collision.chi2() > chi2PvMax) {
    SETBIT(statusCollision, EventRejection::Chi2);
  }

  return statusCollision;
}

/// @brief function to monitor the event selection satisfied by collisions used for HF analyses
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

  /// centrality
  if (TESTBIT(rejectionMask, EventRejection::Centrality)) {
    return;
  }
  hCollisions->Fill(EventRejection::Centrality); // Centrality ok

  /// sel8()
  if (TESTBIT(rejectionMask, EventRejection::Trigger)) {
    return;
  }
  hCollisions->Fill(EventRejection::Trigger); // Centrality + sel8 ok

  /// Time Frame border cut
  if (TESTBIT(rejectionMask, EventRejection::TimeFrameBorderCut)) {
    return;
  }
  hCollisions->Fill(EventRejection::TimeFrameBorderCut); // Centrality + sel8 + TF border ok

  /// PVz consistency tracking - FT0 timing
  if (TESTBIT(rejectionMask, EventRejection::IsGoodZvtxFT0vsPV)) {
    return;
  }
  hCollisions->Fill(EventRejection::IsGoodZvtxFT0vsPV); // Centrality + sel8 + TF border + PVz FTO ok

  /// same bunch pile-up
  if (TESTBIT(rejectionMask, EventRejection::NoSameBunchPileup)) {
    return;
  }
  hCollisions->Fill(EventRejection::NoSameBunchPileup); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up ok

  /// occupancy
  if (TESTBIT(rejectionMask, EventRejection::NumTracksInTimeRange)) {
    return;
  }
  hCollisions->Fill(EventRejection::NumTracksInTimeRange); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up + occupancy ok

  /// PV contributors
  if (TESTBIT(rejectionMask, EventRejection::NContrib)) {
    return;
  }
  hCollisions->Fill(EventRejection::NContrib); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up + occupancy + PV contr ok

  /// PV chi2
  if (TESTBIT(rejectionMask, EventRejection::Chi2)) {
    return;
  }
  hCollisions->Fill(EventRejection::Chi2); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up + occupancy + PV contr + chi2 ok

  /// PV position Z
  if (TESTBIT(rejectionMask, EventRejection::PositionZ)) {
    return;
  }
  hCollisions->Fill(EventRejection::PositionZ); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up + occupancy + PV contr + chi2 + posZ ok

  hPosXAfterEvSel->Fill(collision.posX());
  hPosYAfterEvSel->Fill(collision.posY());
  hPosZAfterEvSel->Fill(posZ);
  hNumContributors->Fill(collision.numContrib());
}

struct HfEvSel : o2::framework::ConfigurableGroup {
  std::string prefix = "hfEvSel"; // JSON group name
  // event selection parameters (in chronological order of application)
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality"};
  Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "apply the sel8 event selection"};
  Configurable<int> triggerClass{"triggerClass", -1, "trigger class different from sel8 (e.g. kINT7 for Run2) used only if useSel8Trigger is false"};
  Configurable<bool> useTimeFrameBorderCut{"useTimeFrameBorderCut", true, "apply TF border cut"};
  Configurable<bool> useIsGoodZvtxFT0vsPV{"useIsGoodZvtxFT0vsPV", false, "check consistency between PVz from central barrel with that from FT0 timing"};
  Configurable<bool> useNoSameBunchPileup{"useNoSameBunchPileup", false, "exclude collisions in bunches with more than 1 reco. PV"}; // POTENTIALLY BAD FOR BEAUTY ANALYSES
  Configurable<bool> useNumTracksInTimeRange{"useNumTracksInTimeRange", false, "apply occupancy selection (num. ITS tracks with at least 5 clusters in +-100us from current collision)"};
  Configurable<int> numTracksInTimeRangeMin{"numTracksInTimeRangeMin", 0, "min. value for occupancy selection"};
  Configurable<int> numTracksInTimeRangeMax{"numTracksInTimeRangeMax", 1000000, "max. value for occupancy selection"};
  Configurable<int> nPvContributorsMin{"nPvContributorsMin", 0, "minimum number of PV contributors"};
  Configurable<float> chi2PvMax{"chi2PvMax", -1.f, "maximum PV chi2"};
  Configurable<float> zPvPosMin{"zPvPosMin", -10.f, "min. PV posZ (cm)"};
  Configurable<float> zPvPosMax{"zPvPosMax", 10.f, "max. PV posZ (cm)"};

  /// @brief Function to put labels on collision monitoring histogram
  /// \param hCollisions is the histogram
  template <typename Histo>
  void setLabelHistoEvSel(Histo& hCollisions)
  {
    hCollisions->SetTitle("HF event counter;;accepted collisions");
    hCollisions->GetXaxis()->SetBinLabel(EventRejection::None + 1, "All collisions");
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
  /// \tparam centEstimator select the centrality estimator
  /// \param collision collision that has to satisfy the selection criteria
  /// \param centrality collision centrality to be initialised in this function
  /// \return a bitmask with the event selections not satisfied by the analysed collision
  template <bool applyEvSel, o2::aod::hf_collision_centrality::CentralityEstimator centEstimator, typename Coll>
  uint16_t getHfCollisionRejectionMask(const Coll& collision, float& centrality)
  {

    uint16_t statusCollision{0}; // 16 bits, in case new ev. selections will be added

    if constexpr (centEstimator != o2::aod::hf_collision_centrality::CentralityEstimator::None) {
      if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FT0A) {
        centrality = collision.centFT0A();
      } else if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FT0C) {
        centrality = collision.centFT0C();
      } else if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FT0M) {
        centrality = collision.centFT0M();
      } else if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FV0A) {
        centrality = collision.centFV0A();
      } else {
        LOGP(fatal, "Centrality estimator different from FT0A, FT0C, FT0M, and FV0A, fix it!");
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

  /// @brief function to monitor the event selection satisfied by collisions used for HF analyses
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

    /// centrality
    if (TESTBIT(rejectionMask, EventRejection::Centrality)) {
      return;
    }
    hCollisions->Fill(EventRejection::Centrality); // Centrality ok

    /// sel8()
    if (TESTBIT(rejectionMask, EventRejection::Trigger)) {
      return;
    }
    hCollisions->Fill(EventRejection::Trigger); // Centrality + sel8 ok

    /// Time Frame border cut
    if (TESTBIT(rejectionMask, EventRejection::TimeFrameBorderCut)) {
      return;
    }
    hCollisions->Fill(EventRejection::TimeFrameBorderCut); // Centrality + sel8 + TF border ok

    /// PVz consistency tracking - FT0 timing
    if (TESTBIT(rejectionMask, EventRejection::IsGoodZvtxFT0vsPV)) {
      return;
    }
    hCollisions->Fill(EventRejection::IsGoodZvtxFT0vsPV); // Centrality + sel8 + TF border + PVz FTO ok

    /// same bunch pile-up
    if (TESTBIT(rejectionMask, EventRejection::NoSameBunchPileup)) {
      return;
    }
    hCollisions->Fill(EventRejection::NoSameBunchPileup); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up ok

    /// occupancy
    if (TESTBIT(rejectionMask, EventRejection::NumTracksInTimeRange)) {
      return;
    }
    hCollisions->Fill(EventRejection::NumTracksInTimeRange); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up + occupancy ok

    /// PV contributors
    if (TESTBIT(rejectionMask, EventRejection::NContrib)) {
      return;
    }
    hCollisions->Fill(EventRejection::NContrib); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up + occupancy + PV contr ok

    /// PV chi2
    if (TESTBIT(rejectionMask, EventRejection::Chi2)) {
      return;
    }
    hCollisions->Fill(EventRejection::Chi2); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up + occupancy + PV contr + chi2 ok

    /// PV position Z
    if (TESTBIT(rejectionMask, EventRejection::PositionZ)) {
      return;
    }
    hCollisions->Fill(EventRejection::PositionZ); // Centrality + sel8 + TF border + PVz FTO + sam bunch pile-up + occupancy + PV contr + chi2 + posZ ok

    hPosXAfterEvSel->Fill(collision.posX());
    hPosYAfterEvSel->Fill(collision.posY());
    hPosZAfterEvSel->Fill(posZ);
    hNumContributors->Fill(collision.numContrib());
  }
};
} // namespace o2::framework::hf_evsel

#endif // PWGHF_UTILS_UTILSEVSELHF_H_
