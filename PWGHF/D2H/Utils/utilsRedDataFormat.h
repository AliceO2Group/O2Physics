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

/// \file utilsRedDataFormat.h
/// \brief Event selection utilities for reduced data format analyses
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin

#ifndef PWGHF_D2H_UTILS_UTILSREDDATAFORMAT_H_
#define PWGHF_D2H_UTILS_UTILSREDDATAFORMAT_H_

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

namespace o2::hf_evsel
{
/// Helper function to count collisions at different event selection stages
/// \tparam useEvSel use information from the EvSel table
/// \tparam centEstimator centrality estimator
/// \param collision collision to test against the selection criteria
template <bool useEvSel, o2::hf_centrality::CentralityEstimator centEstimator, typename Coll>
void checkEvSel(Coll const& collision, o2::hf_evsel::HfEventSelection& hfEvSel, int& zvtxColl, int& sel8Coll, int& zvtxAndSel8Coll, int& allSelColl)
{
  float centrality{-1.f};
  const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<useEvSel, o2::hf_centrality::CentralityEstimator::None>(collision, centrality);
  if (!TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::Trigger))
    sel8Coll++;
  if (!TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::PositionZ))
    zvtxColl++;
  if (!TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::PositionZ) && !TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::Trigger))
    zvtxAndSel8Coll++;
  if (rejectionMask == 0)
    allSelColl++;
}
} // namespace o2::hf_evsel

#endif // PWGHF_D2H_UTILS_UTILSREDDATAFORMAT_H_
