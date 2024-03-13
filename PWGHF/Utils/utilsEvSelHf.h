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

using namespace o2::aod::hf_collision_centrality;

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

/// \brief Function to apply event selections in HF analyses
/// \param collision collision that has to satisfy the selection criteria
/// \param useSel8Trigger switch to activate the sel8() event selection
/// \param zPvPosMax maximum primary-vertex z
/// \param useTimeFrameBorderCut switch to activate the time frame border cut
/// \return true if collision satisfies all criteria, false otherwise
template <int centEstimator = 0, typename Coll>
bool isHfCollisionSelected(const Coll& collision, bool useCentrality, std::array<float, 2> centralityLimits, bool useSel8Trigger, float maxPvPosZ, bool useTimeFrameBorderCut, float& centrality)
{

  uint16_t statusCollision = 0; // 16 bits, in case new ev. selections will be added

  if constexpr (centEstimator != CentralityEstimator::None) {
    if constexpr (centEstimator == CentralityEstimator::FT0C) {
      centrality = collision.centFT0C();
    } else if constexpr (centEstimator == CentralityEstimator::FT0M) {
      centrality = collision.centFT0M();
    } else {
      LOGP(fatal, "Centrality estimator different from FT0C and FT0M, fix it!");
    }
    if (centrality < centralityLimits.at(0) || centrality > centralityLimits.at(1)) {
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
#endif // PWGHF_UTILS_UTILSEVSELHF_H_
