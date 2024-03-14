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

// event rejection types
enum EventRejection {
  Trigger = 0,
  TimeFrameBorderCut,
  PositionX,
  PositionY,
  PositionZ,
  NContrib,
  Chi2,
  Centrality,
  NEventRejection
};

enum ValuesEvSel : int {
  All = 0,
  Cent,
  CentSel8,
  CentSel8PosZ,
  CentSel8PosZTFBorder,
  NEvSel
};

/// @brief Function to put labels on collision monitoring histogram
/// \param hCollisions is the histogram
template <typename Histo>
void setLabelHistoEvSel(Histo& hCollisions)
{
  hCollisions->GetXaxis()->SetBinLabel(ValuesEvSel::All + 1, "All collisions");
  hCollisions->GetXaxis()->SetBinLabel(ValuesEvSel::Cent + 1, "Centrality ok");
  hCollisions->GetXaxis()->SetBinLabel(ValuesEvSel::CentSel8 + 1, "Centrality + sel8 ok");
  hCollisions->GetXaxis()->SetBinLabel(ValuesEvSel::CentSel8PosZ + 1, "Centrality + sel8 + posZ ok");
  hCollisions->GetXaxis()->SetBinLabel(ValuesEvSel::CentSel8PosZTFBorder + 1, "Centrality + sel8 + posZ + TF border ok");
}

/// \brief Function to apply event selections in HF analyses
/// \param collision collision that has to satisfy the selection criteria
/// \param useSel8Trigger switch to activate the sel8() event selection
/// \param zPvPosMax maximum primary-vertex z
/// \param useTimeFrameBorderCut switch to activate the time frame border cut
/// \return a bitmask with the event selections not satisfied by the analysed collision
template <o2::aod::hf_collision_centrality::CentralityEstimator centEstimator, typename Coll>
uint16_t getHfCollisionRejectionMask(const Coll& collision, float centralityMin, float centralityMax, bool useSel8Trigger, float maxPvPosZ, bool useTimeFrameBorderCut)
{

  uint16_t statusCollision = 0; // 16 bits, in case new ev. selections will be added
  float centrality = -1.;

  if constexpr (centEstimator != o2::aod::hf_collision_centrality::CentralityEstimator::None) {
    if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FT0C) {
      centrality = collision.centFT0C();
    } else if constexpr (centEstimator == o2::aod::hf_collision_centrality::CentralityEstimator::FT0M) {
      centrality = collision.centFT0M();
    } else {
      LOGP(fatal, "Centrality estimator different from FT0C and FT0M, fix it!");
    }
    if (centrality < centralityMin || centrality > centralityMax) {
      SETBIT(statusCollision, EventRejection::Centrality);
    }
  }

  /// sel8() condition
  if (useSel8Trigger && !collision.sel8()) {
    SETBIT(statusCollision, EventRejection::Trigger);
  }

  /// primary vertex z
  if (std::fabs(collision.posZ()) > maxPvPosZ) {
    SETBIT(statusCollision, EventRejection::PositionZ);
  }

  /// time frame border cut
  if (useTimeFrameBorderCut && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    SETBIT(statusCollision, EventRejection::TimeFrameBorderCut);
  }

  /// TODO: add other selections, to extend it to the trackIndexSkimCreator

  return statusCollision;
}

/// @brief function to monitor the event selection satisfied by collisions used for HF analyses
/// \param collision is the analysed collision
/// \param rejectionMask is the bitmask storing the info about which ev. selections are not satisfied by the collision
/// \param hCollisions is a histogram to keep track of the satisfied event selections
/// \param hPosZBeforeEvSel is PV position Z for all analysed collisions
/// \param hPosZAfterEvSel is PV position Z only for collisions satisfying the event selections
template <typename Coll, typename Hist>
void monitorCollision(Coll const& collision, const uint16_t rejectionMask, Hist& hCollisions, Hist& hPosZBeforeEvSel, Hist& hPosZAfterEvSel)
{

  hCollisions->Fill(ValuesEvSel::All); // all collisions
  const float posZ = collision.posZ();
  hPosZBeforeEvSel->Fill(posZ);

  /// centrality
  if (TESTBIT(rejectionMask, EventRejection::Centrality)) {
    return;
  }
  hCollisions->Fill(ValuesEvSel::Cent); // Centrality ok

  /// sel8()
  if (TESTBIT(rejectionMask, EventRejection::Trigger)) {
    return;
  }
  hCollisions->Fill(ValuesEvSel::CentSel8); // Centrality + sel8 ok

  /// PV position Z
  if (TESTBIT(rejectionMask, EventRejection::PositionZ)) {
    return;
  }
  hCollisions->Fill(ValuesEvSel::CentSel8PosZ); // Centrality + sel8 + posZ ok

  /// Time Frame border cut
  if (TESTBIT(rejectionMask, EventRejection::TimeFrameBorderCut)) {
    return;
  }
  hCollisions->Fill(ValuesEvSel::CentSel8PosZTFBorder); // Centrality + sel8 + posZ + TF border ok
  hPosZAfterEvSel->Fill(posZ);
}

#endif // PWGHF_UTILS_UTILSEVSELHF_H_
