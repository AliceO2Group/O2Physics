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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN

#ifndef O2_ANALYSIS_DERIVEDEXAMPLETABLES_H
#define O2_ANALYSIS_DERIVEDEXAMPLETABLES_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
DECLARE_SOA_TABLE(DerivedCollisions, "AOD", "DERCOLLS", o2::soa::Index<>,
                  o2::aod::collision::PosZ);
using DerivedCollision = DerivedCollisions::iterator;

namespace exampleTrackSpace
{
DECLARE_SOA_INDEX_COLUMN(DerivedCollision, derivedCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace exampleTrackSpace
DECLARE_SOA_TABLE(DerivedTracks, "AOD", "DERTRACKS", o2::soa::Index<>, exampleTrackSpace::DerivedCollisionId,
                  exampleTrackSpace::Pt, exampleTrackSpace::Eta, exampleTrackSpace::Phi);
using DerivedTrack = DerivedTracks::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_DERIVEDEXAMPLETABLES_H
