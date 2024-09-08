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
#ifndef PWGCF_JCORRAN_DATAMODEL_JCATALYST_H_
#define PWGCF_JCORRAN_DATAMODEL_JCATALYST_H_

namespace o2::aod
{
namespace jcollision
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Collision centrality or multiplicity
} // namespace jcollision

DECLARE_SOA_TABLE(JCollisions, "AOD", "JCOLLISION", //! Reduced collision table
                  o2::soa::Index<>,
                  bc::RunNumber,
                  collision::PosZ,
                  jcollision::Multiplicity);
using JCollision = JCollisions::iterator;

namespace jtrack
{
DECLARE_SOA_INDEX_COLUMN(JCollision, jcollision); //! collision ID
DECLARE_SOA_COLUMN(Pt, pt, float);                //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);              //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);              //! Phi
DECLARE_SOA_COLUMN(Sign, sign, int8_t);           //! Phi
} // namespace jtrack

DECLARE_SOA_TABLE(JTracks, "AOD", "JTRACK", //! Reduced track table
                  o2::soa::Index<>,
                  jtrack::JCollisionId,
                  jtrack::Pt, jtrack::Eta, jtrack::Phi, jtrack::Sign);
using JTrack = JTracks::iterator;

namespace jweight
{
DECLARE_SOA_COLUMN(WeightNUA, weightNUA, float); //! Non-uniform acceptance weight
DECLARE_SOA_COLUMN(WeightEff, weightEff, float); //! Non-uniform efficiency weight
} // namespace jweight
DECLARE_SOA_TABLE(JWeights, "AOD", "JWEIGHT", jweight::WeightNUA, jweight::WeightEff); //! JFluc table for weights

namespace j2prongweight
{
DECLARE_SOA_COLUMN(WeightNUA, weightNUA, float); //! Non-uniform acceptance weight
DECLARE_SOA_COLUMN(WeightEff, weightEff, float); //! Non-uniform efficiency weight
} // namespace j2prongweight
DECLARE_SOA_TABLE(J2ProngWeights, "AOD", "J2PRONGWEIGHT", j2prongweight::WeightNUA, j2prongweight::WeightEff); //! JFluc table for weights, associated with 2Prong particles

} // namespace o2::aod

#endif // PWGCF_JCORRAN_DATAMODEL_JCATALYST_H_
