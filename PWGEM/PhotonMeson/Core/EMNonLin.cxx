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

using namespace o2::pwgem::nonlin;

float EMNonLin::getCorrectionFactor(float x, const Context& ctx)
{
  if (!ctx.params || x == 0.f) [[unlikely]] {
    return x;
  }

  float val = x;
  // safety measure
  int maxIter = std::min(ctx.nIter, MaxIter - 1);

  for (int i = 0; i <= maxIter; ++i) {
    if (val == 0.f) {
      break;
    }

    float inv = 1.f / val;
    val *= (1.f + ctx.params[i].par0 * inv +
            ctx.params[i].par1 * inv * inv) /
           (1.f + ctx.params[i].par2 * inv);
  }

  return val;
}

const EMNonLin::NonLinParams* EMNonLin::resolveParams(PhotonType type, float cent)
{
  int centBin = static_cast<int>(cent / 10.f);
  if (centBin < 0)
    centBin = 0;
  if (centBin >= CentBins)
    centBin = CentBins - 1;

  return &kNonLinTable[static_cast<int>(type)][centBin][0];
}

const EMNonLin::NonLinParams* EMNonLin::resolveParamsMC(PhotonType type, float cent)
{
  int centBin = static_cast<int>(cent / 10.f);
  if (centBin < 0)
    centBin = 0;
  if (centBin >= CentBins)
    centBin = CentBins - 1;

  return &kNonLinTableMC[static_cast<int>(type)][centBin][0];
}
