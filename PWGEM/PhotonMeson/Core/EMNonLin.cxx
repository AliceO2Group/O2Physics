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

/// \file EMNonLin.cxx
/// \brief Source file of NonLin class for photons.
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "EMNonLin.h"

#include <algorithm>
#include <cmath>

using namespace o2::pwgem::nonlin;

float EMNonLin::getCorrectionFactor(float var, const Context& ctx)
{
  if (!ctx.params || var == 0.f) [[unlikely]] {
    return 1.f;
  }

  int maxIter = std::min(ctx.nIter, MaxIter - 1);
  float scale = 1.f;  // cumulative scale
  float refVal = var; // reference value updated each iteration

  for (int i = 0; i <= maxIter; ++i) {
    if (refVal == 0.f) {
      break;
    }
    const auto& p = ctx.params[i];

    // evaluate pol1 for each parameter at this centrality
    float a = p.a0 + p.a1 * ctx.cent;
    float b = p.b0 + p.b1 * ctx.cent;
    float c = p.c0 + p.c1 * ctx.cent;

    // guard against c <= 0 which would make pow(x, -c) diverge
    if (c <= 0.f) {
      continue;
    }

    float iterScale = a + b * std::pow(refVal, -c);
    scale *= iterScale;
    refVal = var * scale; // next iteration uses scaled original input
  }

  return scale;
}
