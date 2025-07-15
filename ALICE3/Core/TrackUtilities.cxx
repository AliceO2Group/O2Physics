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
/// \file   TrackUtilities.cxx
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
/// \brief  Set of utilities for the ALICE3 track handling
/// \since  May 21, 2025
///

#include "TrackUtilities.h"

#include <vector>

void o2::upgrade::convertTLorentzVectorToO2Track(const int charge,
                                                 const TLorentzVector particle,
                                                 const std::vector<double> productionVertex,
                                                 o2::track::TrackParCov& o2track)
{
  std::array<float, 5> params;
  std::array<float, 15> covm = {0.};
  float s, c, x;
  o2::math_utils::sincos(static_cast<float>(particle.Phi()), s, c);
  o2::math_utils::rotateZInv(static_cast<float>(productionVertex[0]), static_cast<float>(productionVertex[1]), x, params[0], s, c);
  params[1] = static_cast<float>(productionVertex[2]);
  params[2] = 0.; // since alpha = phi
  const auto theta = 2. * std::atan(std::exp(-particle.PseudoRapidity()));
  params[3] = 1. / std::tan(theta);
  params[4] = charge / particle.Pt();

  // Initialize TrackParCov in-place
  new (&o2track)(o2::track::TrackParCov)(x, particle.Phi(), params, covm);
}
