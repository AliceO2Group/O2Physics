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

#ifndef COMMON_DATAMODEL_CALOCLUSTERS_H_
#define COMMON_DATAMODEL_CALOCLUSTERS_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace calocluster
{
// Columns to store momenta of "photons"
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //! collisionID used as index for matched clusters
DECLARE_SOA_INDEX_COLUMN(BC, bc);               //! BC index
DECLARE_SOA_COLUMN(Px, px, float);              //! momenta components
DECLARE_SOA_COLUMN(Py, py, float);              //!
DECLARE_SOA_COLUMN(Pz, pz, float);              //!
DECLARE_SOA_COLUMN(E, e, float);                //!

// Columns to store cluster PID parameters
DECLARE_SOA_COLUMN(Module, mod, uint8_t);                //! module/supermodule number
DECLARE_SOA_COLUMN(Ncell, ncell, uint8_t);               //! cluster multiplicity
DECLARE_SOA_COLUMN(X, x, float);                         //! cluster local coordinates
DECLARE_SOA_COLUMN(Z, z, float);                         //!
DECLARE_SOA_COLUMN(GlobalX, globalx, float);             //! cluster global coordinates
DECLARE_SOA_COLUMN(GlobalY, globaly, float);             //! cluster global coordinates
DECLARE_SOA_COLUMN(GlobalZ, globalz, float);             //! cluster global coordinates
DECLARE_SOA_COLUMN(Time, time, float);                   //! cluster time (seconds)
DECLARE_SOA_COLUMN(NLM, nlm, uint8_t);                   //! number of local maxima
DECLARE_SOA_COLUMN(M02, m02, float);                     //! longer dispersion axis
DECLARE_SOA_COLUMN(M20, m20, float);                     //! shorter dispersion axis
DECLARE_SOA_COLUMN(TrackDist, trackdist, float);         //! distance to closest track
DECLARE_SOA_COLUMN(TrackIndex, trackIndex, uint8_t);     //! index of closest track
DECLARE_SOA_COLUMN(FiredTrigger, firedTrigger, uint8_t); //! Matched with trigger tile
DECLARE_SOA_COLUMN(DistBad, distBad, float);             //! distance to closest bad channel

} // namespace calocluster

DECLARE_SOA_TABLE(CaloClusters, "AOD", "CALOCLUSTERS", //!
                  o2::soa::Index<>, calocluster::CollisionId,
                  calocluster::Px, calocluster::Py, calocluster::Pz, calocluster::E,
                  calocluster::Module, calocluster::Ncell,
                  calocluster::X, calocluster::Z,
                  calocluster::GlobalX, calocluster::GlobalY, calocluster::GlobalZ,
                  calocluster::Time, calocluster::NLM, calocluster::M02, calocluster::M20,
                  calocluster::TrackDist, calocluster::TrackIndex, calocluster::FiredTrigger, calocluster::DistBad);

// table of ambiguous clusters that could not be matched to a collision
DECLARE_SOA_TABLE(CaloAmbiguousClusters, "AOD", "CALOAMBCLUS", //!
                  o2::soa::Index<>, calocluster::BCId,
                  calocluster::Px, calocluster::Py, calocluster::Pz, calocluster::E,
                  calocluster::Module, calocluster::Ncell,
                  calocluster::X, calocluster::Z,
                  calocluster::GlobalX, calocluster::GlobalY, calocluster::GlobalZ,
                  calocluster::Time, calocluster::NLM, calocluster::M02, calocluster::M20,
                  calocluster::TrackDist, calocluster::TrackIndex, calocluster::FiredTrigger, calocluster::DistBad);

using CaloCluster = CaloClusters::iterator;
using CaloAMBCluster = CaloAmbiguousClusters::iterator;

} // namespace o2::aod

#endif // COMMON_DATAMODEL_CALOCLUSTERS_H_
