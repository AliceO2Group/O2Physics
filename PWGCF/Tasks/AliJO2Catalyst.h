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
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \since Sep 2022
#ifndef JFLUC_CATALYST_H
#define JFLUC_CATALYST_H

namespace o2::aod
{
namespace particleTrack
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(WeightNUA, weightNUA, float);
DECLARE_SOA_COLUMN(WeightEff, weightEff, float);
} // namespace particleTrack

namespace collisionData
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(CBin, cbin, Int_t);
} // namespace collisionData

DECLARE_SOA_TABLE(ParticleTrack, "AOD", "PARTICLETRACK",
                  particleTrack::CollisionId,
                  particleTrack::Pt, particleTrack::Eta, particleTrack::Phi,
                  particleTrack::WeightNUA, particleTrack::WeightEff);
DECLARE_SOA_TABLE(CollisionData, "AOD", "COLLISIONDATA",
                  collisionData::CollisionId,
                  collisionData::Cent,
                  collisionData::CBin);
} // namespace o2::aod

static float jflucCentBins[] = {0.0f, 1.0f, 2.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f};

#endif
