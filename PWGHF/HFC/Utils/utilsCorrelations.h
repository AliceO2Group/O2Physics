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

/// \file utilsCorrelations.h
/// \brief Utilities for HFC analyses

#ifndef UTILS_CORRELATIONS_H
#define UTILS_CORRELATIONS_H

#include <TString.h>
#include <cmath>

enum Region { Toward,
              Away,
              Transverse };

Region getRegion(double deltaPhi)
{
  if (std::abs(deltaPhi) < o2::constants::math::PI / 3.) {
    return Toward;
  } else if (deltaPhi > 2. * o2::constants::math::PI / 3. && deltaPhi < 4. * o2::constants::math::PI / 3.) {
    return Away;
  } else {
    return Transverse;
  }
}

#endif // UTILS_CORRELATIONS_H