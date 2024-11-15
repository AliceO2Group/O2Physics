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

/// \file utilsDerivedData.h
/// \brief Utilities for derived-data creators
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_UTILS_UTILSDERIVEDDATA_H_
#define PWGHF_UTILS_UTILSDERIVEDDATA_H_

#include <Framework/Configurable.h>

// Macro to store nSigma for prong _id_ with PID hypothesis _hyp_ in an array
#define GET_N_SIGMA_PRONG(_array_, _candidate_, _id_, _hyp_) \
  _array_[0] = _candidate_.nSigTpc##_hyp_##_id_();           \
  _array_[1] = _candidate_.nSigTof##_hyp_##_id_();           \
  _array_[2] = _candidate_.tpcTofNSigma##_hyp_##_id_();

namespace o2::analysis::hf_derived
{
/// Reserve space in the filled table for all entries in the source table.
/// \param cursor  cursor of the filled table
/// \param enabled  switch for filling the table
/// \param size  size of the source table
template <typename T>
void reserveTable(T& cursor, const o2::framework::Configurable<bool>& enabled, const uint64_t size)
{
  if (enabled.value) {
    cursor.reserve(size);
  }
};
} // namespace o2::analysis::hf_derived

#endif // PWGHF_UTILS_UTILSDERIVEDDATA_H_
