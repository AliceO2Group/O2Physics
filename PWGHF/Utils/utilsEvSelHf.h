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

/// \brief Function to apply event selections in HF analyses
/// \param collision collision that has to satisfy the selection criteria
/// \param useSel8Trigger switch to activate the sel8() event selection
/// \param zPvPosMax maximum primary-vertex z
/// \param useTimeFrameBorderCut switch to activate the time frame border cut
/// \return true if collision satisfies all criteria, false otherwise
template <typename Coll>
bool isHfCollisionSelected(const Coll& collision, bool useSel8Trigger, float maxPvPosZ, bool useTimeFrameBorderCut)
{

  /// sel8() condition
  if (useSel8Trigger && !collision.sel8()) {
    return false;
  }

  /// primary vertex z
  if (std::fabs(collision.posZ()) > maxPvPosZ) {
    return false;
  }

  /// time frame border cut
  if (useTimeFrameBorderCut && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    return false;
  }

  /// all conditions satisfied
  return true;
}
#endif // PWGHF_UTILS_UTILSEVSELHF_H_
