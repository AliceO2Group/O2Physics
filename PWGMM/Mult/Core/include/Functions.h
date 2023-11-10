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

#ifndef PWGMM_FUNCTIONS_H
#define PWGMM_FUNCTIONS_H
#include <cmath>
#include "Framework/AnalysisDataModel.h"

namespace pwgmm::mult {
using namespace o2;

// helper function to determine if collision/mccollison type contains centrality
template <typename T>
static constexpr bool hasCent()
{
  if constexpr (!soa::is_soa_join_v<T>) {
    return false;
  } else if constexpr (T::template contains<aod::HepMCHeavyIons>()) {
    return true;
  } else {
    return false;
  }
}
}

#endif // PWGMM_FUNCTIONS_H
