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

///
/// \file   EventPlaneHelper.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Helper methods to calculate information needed for the event plane
///         calculations with FIT.
///

#include "Common/Core/EventPlaneHelper.h"

#include <FT0Base/Geometry.h>
#include <FV0Base/Geometry.h>

#include <TComplex.h>
#include <TH2.h>
#include <TMath.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <memory>
#include <vector>

double EventPlaneHelper::GetPhiFV0(int chno, o2::fv0::Geometry* fv0geom)
{
  /* Calculate the azimuthal angle in FV0 for the channel number 'chno'. The offset
  on the A-side is taken into account here. */

  // Set the offset of the left and right side of FV0.
  float offsetX = 0.;
  float offsetY = 0.;

  const int cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27,
                             32, 40, 33, 41, 34, 42, 35, 43};
  bool isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);

  if (isChnoInLeft) {
    offsetX = mOffsetFV0leftX;
    offsetY = mOffsetFV0leftY;
  } else {
    offsetX = mOffsetFV0rightX;
    offsetY = mOffsetFV0rightY;
  }

  auto chPos = fv0geom->getReadoutCenter(chno);

  return TMath::ATan2(chPos.y + offsetY, chPos.x + offsetX);
}

double EventPlaneHelper::GetPhiFT0(int chno, o2::ft0::Geometry ft0geom)
{
  /* Calculate the azimuthal angle in FT0 for the channel number 'chno'. The offset
    of FT0-A is taken into account if chno is between 0 and 95. */

  float offsetX = 0.;
  float offsetY = 0.; // No offset for FT0-C (default case).

  if (chno < 96) { // Channel in FT0-A, non-zero offset must be applied. // LOKI: make general.
    offsetX = mOffsetFT0AX;
    offsetY = mOffsetFT0AY;
  }

  ft0geom.calculateChannelCenter();
  auto chPos = ft0geom.getChannelCenter(chno);
  /// printf("Channel id: %d X: %.3f Y: %.3f\n", chno, chPos.X(), chPos.Y());

  return TMath::ATan2(chPos.Y() + offsetY, chPos.X() + offsetX);
}

void EventPlaneHelper::SumQvectors(int det, int chno, float ampl, int nmod, TComplex& Qvec, float& sum, o2::ft0::Geometry ft0geom, o2::fv0::Geometry* fv0geom)
{
  /* Calculate the complex Q-vector for the provided detector and channel number,
    before adding it to the total Q-vector given as argument. */
  double phi = -999.;

  switch (det) {
    case 0:                           // FT0. Note: the channel number for FT0-C should
      phi = GetPhiFT0(chno, ft0geom); // already be given in the right range.
      break;
    case 1: // FV0.
      phi = GetPhiFV0(chno, fv0geom);
      break;
    default:
      printf("'int det' value does not correspond to any accepted case.\n");
      break;
  }

  /// printf("Phi: %.3f\n", phi);

  if (phi < -900) {
    printf("Error on phi. Skip\n");
    return;
  } // TODO: ensure proper safety check.
  Qvec += TComplex(ampl * TMath::Cos(phi * nmod), ampl * TMath::Sin(phi * nmod));
  sum += ampl;
}

int EventPlaneHelper::GetCentBin(float cent)
{
  const float centClasses[] = {0., 5., 10., 20., 30., 40., 50., 60., 80.};

  for (int i = 0; i < 9; i++) {
    if (cent >= centClasses[i]) {
      continue;
    } else {
      return i - 1;
    }
  }

  // We went through all centrality edges without returning --> The measured percentile
  // is larger than the final class we consider.
  return -1;
}

void EventPlaneHelper::DoCorrections(float& qx, float& qy,
                                     const std::vector<float>& corrections)
{
  // Recentering of the Qx-Qy distribution to (0,0).
  qx -= corrections[0];
  qy -= corrections[1];

  // Twisting of the Qx-Qy distribution.
  qx = (qx - corrections[3] * qy) / (1.0 - corrections[3] * corrections[2]);
  qy = (qy - corrections[2] * qx) / (1.0 - corrections[3] * corrections[2]);

  // Rescaling of the Qx-Qy into a circle.
  if (std::fabs(corrections[4]) > 1e-8) {
    qx /= corrections[4];
  }
  if (std::fabs(corrections[5]) > 1e-8) {
    qy /= corrections[5];
  }
}

void EventPlaneHelper::DoRecenter(float& qx, float& qy, float x0, float y0)
{
  qx -= x0;
  qy -= y0;
}

void EventPlaneHelper::DoTwist(float& qx, float& qy, float lp, float lm)
{
  qx = (qx - lm * qy) / (1.0 - lm * lp);
  qy = (qy - lp * qx) / (1.0 - lm * lp);
}

void EventPlaneHelper::DoRescale(float& qx, float& qy, float ap, float am)
{
  if (std::fabs(ap) > 1e-8)
    qx /= ap;
  if (std::fabs(am) > 1e-8)
    qy /= am;
}

void EventPlaneHelper::GetCorrRecentering(const std::shared_ptr<TH2> histQ, float& meanX, float& meanY)
{
  meanX = histQ->GetMean(1);
  meanY = histQ->GetMean(2);
}

void EventPlaneHelper::GetCorrWidth(const std::shared_ptr<TH2> histQ, float& stdX, float& stdY)
{
  stdX = histQ->GetStdDev(1);
  stdY = histQ->GetStdDev(2);
}

void EventPlaneHelper::GetCorrTwistRecale(const std::shared_ptr<TH2> histQ,
                                          float& aPlus, float& aMinus,
                                          float& lambdaPlus, float& lambdaMinus)
{
  // Get first information from the provided TH2D.
  float rho = histQ->GetCorrelationFactor();
  float sigmax = histQ->GetStdDev(1);
  float sigmay = histQ->GetStdDev(2);

  // Combine them in the "b" parameter LOKI: define it better in the comment.
  float b = rho * sigmax * sigmay * TMath::Sqrt(2.0 * (sigmax * sigmax + sigmay * sigmay - 2.0 * sigmax * sigmay * TMath::Sqrt(1.0 - rho * rho)) / ((sigmax * sigmax - sigmay * sigmay) * (sigmax * sigmax - sigmay * sigmay) + 4.0 * (sigmax * sigmay * rho) * (sigmax * sigmay * rho)));

  // Calculate finally the correction constants.
  aPlus = TMath::Sqrt(2. * TMath::Power(sigmax, 2.) - TMath::Power(b, 2.));
  aMinus = TMath::Sqrt(2. * TMath::Power(sigmay, 2.) - TMath::Power(b, 2.));
  lambdaPlus = b / aPlus;
  lambdaMinus = b / aMinus;
}

float EventPlaneHelper::GetEventPlane(const float qx, const float qy, int nmode)
{
  return (1. / nmode) * (TMath::ATan2(qy, qx));
}

float EventPlaneHelper::GetResolution(const float RefA, const float RefB, int nmode)
{
  return std::cos((RefA - RefB) * nmode);
}
