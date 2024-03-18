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
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \since Sep 2022
#ifndef JFLUC_CATALYST_H
#define JFLUC_CATALYST_H

namespace o2::aod
{
namespace jcollision
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Collision centrality or multiplicity
} // namespace collisionData

DECLARE_SOA_TABLE(JCollisions, "AOD", "JCOLLISION",
                  o2::soa::Index<>,
				  bc::RunNumber,
				  collision::PosZ,
				  jcollision::Multiplicity);
using JCollision = JCollisions::iterator;

namespace jtrack
{
DECLARE_SOA_INDEX_COLUMN(JCollision, jcollision);  //! collision ID
DECLARE_SOA_COLUMN(Pt, pt, float);               //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);             //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);             //! Phi
DECLARE_SOA_COLUMN(Sign, sign, int8_t);             //! Phi
} // namespace particleTrack

DECLARE_SOA_TABLE(JTracks, "AOD", "JTRACK",
				  o2::soa::Index<>,
                  jtrack::JCollisionId,
                  jtrack::Pt, jtrack::Eta, jtrack::Phi, jtrack::Sign);
using JTrack = JTracks::iterator;
} // namespace o2::aod

#endif
