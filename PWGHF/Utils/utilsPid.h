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

/// \file CandidateReconstructionTables.h
/// \brief Definitions of tables produced by candidate reconstruction workflows
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, CERN

#ifndef PWGHF_UTILS_PID_H_
#define PWGHF_UTILS_PID_H_

namespace o2::aod
{

namespace pid_tpc_tof_utils
{

/// Function to combine TPC and TOF NSigma
/// \param tpcNSigma is the (binned) NSigma separation in TPC (if tiny = true)
/// \param tofNSigma is the (binned) NSigma separation in TOF (if tiny = true)
/// \return combined NSigma of TPC and TOF
template <typename T1>
T1 combineNSigma(T1 tpcNSigma, T1 tofNSigma)
{
  static constexpr float defaultNSigmaTolerance = .1f;
  static constexpr float defaultNSigma = -999.f + defaultNSigmaTolerance; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

  if ((tpcNSigma > defaultNSigma) && (tofNSigma > defaultNSigma)) { // TPC and TOF
    return std::sqrt(.5f * (tpcNSigma * tpcNSigma + tofNSigma * tofNSigma));
  }
  if (tpcNSigma > defaultNSigma) { // only TPC
    return std::abs(tpcNSigma);
  }
  if (tofNSigma > defaultNSigma) { // only TOF
    return std::abs(tofNSigma);
  }
  return tofNSigma; // no TPC nor TOF
}

} // namespace pid_tpc_tof_utils

} // namespace o2::aod

#endif // PWGHF_UTILS_PID_H_
