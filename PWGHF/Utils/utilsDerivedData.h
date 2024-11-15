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

#include <Framework/Core/include/Framework/Configurable.h>

namespace o2::analysis::hf_derived
{
  template <typename T>
  void reserveTable(T& cursor, const o2::framework::Configurable<bool>& enabled, const uint64_t size)
  {
    if (enabled.value) {
      cursor.reserve(size);
    }
  };
} // namespace o2::analysis::hf_derived

#endif // PWGHF_UTILS_UTILSDERIVEDDATA_H_
