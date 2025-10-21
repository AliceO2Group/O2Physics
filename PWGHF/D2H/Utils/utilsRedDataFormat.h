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
/// \brief Utilities for reduced data format analyses
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin

#ifndef PWGHF_D2H_UTILS_UTILSREDDATAFORMAT_H_
#define PWGHF_D2H_UTILS_UTILSREDDATAFORMAT_H_

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/HistogramRegistry.h>

#include <Rtypes.h>

#include <cmath>

namespace o2::hf_evsel
{
/// Helper function to count collisions at different event selection stages
/// \tparam useEvSel use information from the EvSel table
/// \tparam centEstimator centrality estimator
/// \param collision collision to test against the selection criteria
template <bool UseEvSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename BCs, typename Coll>
void checkEvSel(Coll const& collision, o2::hf_evsel::HfEventSelection& hfEvSel, int& zvtxColl, int& sel8Coll, int& zvtxAndSel8Coll, int& zvtxAndSel8CollAndSoftTrig, int& allSelColl, o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, o2::framework::HistogramRegistry& registry)
{
  float centrality{-1.f};
  const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<UseEvSel, o2::hf_centrality::CentralityEstimator::None, BCs>(collision, centrality, ccdb, registry);
  if (!TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::Trigger)) {
    sel8Coll++;
  }
  if (!TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::PositionZ)) {
    zvtxColl++;
  }
  if (!TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::PositionZ) && !TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::Trigger)) {
    zvtxAndSel8Coll++;
  }
  if (!TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::PositionZ) && !TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::Trigger) && !TESTBIT(rejectionMask, o2::hf_evsel::EventRejection::SoftwareTrigger)) {
    zvtxAndSel8CollAndSoftTrig++;
  }
  if (rejectionMask == 0) {
    allSelColl++;
  }
}
} // namespace o2::hf_evsel

namespace o2::pid_tpc_tof_utils
{
/// Helper function to retrive PID information of bachelor pion from b-hadron decay
/// \param prong1 pion track from reduced data format, soa::Join<HfRedTracks, HfRedTracksPid>
template <typename T1>
float getTpcTofNSigmaPi1(const T1& prong1)
{
  float const defaultNSigma = -999.f; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

  bool hasTpc = prong1.hasTPC();
  bool hasTof = prong1.hasTOF();

  if (hasTpc && hasTof) {
    float tpcNSigma = prong1.tpcNSigmaPi();
    float tofNSigma = prong1.tofNSigmaPi();
    return std::sqrt(.5f * tpcNSigma * tpcNSigma + .5f * tofNSigma * tofNSigma);
  }
  if (hasTpc) {
    return std::abs(prong1.tpcNSigmaPi());
  }
  if (hasTof) {
    return std::abs(prong1.tofNSigmaPi());
  }
  return defaultNSigma;
}

/// Helper function to retrive PID information of bachelor pion from b-hadron decay
/// \param prongSoftPi soft pion track
template <typename T1>
float getTpcTofNSigmaSoftPi(const T1& prongSoftPi)
{
  float const defaultNSigma = -999.f; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

  bool hasTpc = prongSoftPi.hasTPC();
  bool hasTof = prongSoftPi.hasTOF();

  if (hasTpc && hasTof) {
    float tpcNSigma = prongSoftPi.tpcNSigmaPiSoftPi();
    float tofNSigma = prongSoftPi.tofNSigmaPiSoftPi();
    return std::sqrt(.5f * tpcNSigma * tpcNSigma + .5f * tofNSigma * tofNSigma);
  }
  if (hasTpc) {
    return std::abs(prongSoftPi.tpcNSigmaPiSoftPi());
  }
  if (hasTof) {
    return std::abs(prongSoftPi.tofNSigmaPiSoftPi());
  }
  return defaultNSigma;
}

/// Helper function to retrive PID information of bachelor kaon from b-hadron decay
/// \param prong1 kaon track from reduced data format, aod::HfRedBachProng0Tracks
/// \return the combined TPC and TOF n-sigma for kaon
template <typename T1>
float getTpcTofNSigmaKa1(const T1& prong1)
{
  float const defaultNSigma = -999.f; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

  bool hasTpc = prong1.hasTPC();
  bool hasTof = prong1.hasTOF();

  if (hasTpc && hasTof) {
    float tpcNSigma = prong1.tpcNSigmaKa();
    float tofNSigma = prong1.tofNSigmaKa();
    return std::sqrt(.5f * tpcNSigma * tpcNSigma + .5f * tofNSigma * tofNSigma);
  }
  if (hasTpc) {
    return std::abs(prong1.tpcNSigmaKa());
  }
  if (hasTof) {
    return std::abs(prong1.tofNSigmaKa());
  }
  return defaultNSigma;
}
} // namespace o2::pid_tpc_tof_utils

#endif // PWGHF_D2H_UTILS_UTILSREDDATAFORMAT_H_
