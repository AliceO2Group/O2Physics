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
#ifndef O2_ANALYSIS_CFDERIVED_H
#define O2_ANALYSIS_CFDERIVED_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Centrality.h"

namespace o2::aod
{
namespace cfmccollision
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Centrality/multiplicity value
} // namespace cfmccollision
DECLARE_SOA_TABLE(CFMcCollisions, "AOD", "CFMCCOLLISION", //! Reduced MC collision table
                  o2::soa::Index<>,
                  mccollision::PosZ, cfmccollision::Multiplicity);
using CFMcCollision = CFMcCollisions::iterator;

namespace cfmcparticle
{
DECLARE_SOA_INDEX_COLUMN(CFMcCollision, cfMcCollision); //! Index to reduced MC collision
DECLARE_SOA_COLUMN(Pt, pt, float);                      //! pT (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                    //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                    //! Phi angle
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                 //! Sign (positive, negative)
const uint8_t kReconstructed = 0x80;
} // namespace cfmcparticle
DECLARE_SOA_TABLE(CFMcParticles, "AOD", "CFMCPARTICLE", //! Reduced MC particle table
                  o2::soa::Index<>,
                  cfmcparticle::CFMcCollisionId,
                  cfmcparticle::Pt, cfmcparticle::Eta, cfmcparticle::Phi, cfmcparticle::Sign,
                  mcparticle::PdgCode, mcparticle::Flags,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);
using CFMcParticle = CFMcParticles::iterator;

namespace cfcollision
{
DECLARE_SOA_INDEX_COLUMN(CFMcCollision, cfMcCollision); //! Index to reduced MC collision
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);  //! Centrality/multiplicity value
} // namespace cfcollision
DECLARE_SOA_TABLE(CFCollisions, "AOD", "CFCOLLISION", //! Reduced collision table
                  o2::soa::Index<>,
                  bc::RunNumber, collision::PosZ,
                  cfcollision::Multiplicity, timestamp::Timestamp);
DECLARE_SOA_TABLE(CFCollLabels, "AOD", "CFCOLLLABEL", //! Labels for reduced collision table
                  cfcollision::CFMcCollisionId);
using CFCollision = CFCollisions::iterator;
using CFCollLabel = CFCollLabels::iterator;
using CFCollisionsWithLabel = soa::Join<CFCollisions, CFCollLabels>;
using CFCollisionWithLabel = CFCollisionsWithLabel::iterator;

namespace cftrack
{
DECLARE_SOA_INDEX_COLUMN(CFCollision, cfCollision);   //! Index to collision
DECLARE_SOA_INDEX_COLUMN(CFMcParticle, cfMCParticle); //! Index to MC particle
DECLARE_SOA_COLUMN(Pt, pt, float);                    //! pT (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                  //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                  //! Phi angle
DECLARE_SOA_COLUMN(Sign, sign, int8_t);               //! Sign (positive, negative)
} // namespace cftrack
DECLARE_SOA_TABLE(CFTracks, "AOD", "CFTRACK", //! Reduced track table
                  o2::soa::Index<>,
                  cftrack::CFCollisionId,
                  cftrack::Pt, cftrack::Eta, cftrack::Phi,
                  cftrack::Sign, track::TrackType);
DECLARE_SOA_TABLE(CFTrackLabels, "AOD", "CFTRACKLABEL", //! Labels for reduced track table
                  cftrack::CFMcParticleId);
using CFTrack = CFTracks::iterator;
using CFTrackLabel = CFTrackLabels::iterator;
using CFTracksWithLabel = soa::Join<CFTracks, CFTrackLabels>;
using CFTrackWithLabel = CFTracksWithLabel::iterator;

//------transient CF-filter to CF-2prong-filter
DECLARE_SOA_TABLE(CFCollRefs, "AOD", "CFCOLLREF", o2::soa::Index<>, track::CollisionId); //! Transient cf collision index table

using CFCollRef = CFCollRefs::iterator;

namespace cftrackref
{
DECLARE_SOA_INDEX_COLUMN(Track, track);
} // namespace cftrackref
DECLARE_SOA_TABLE(CFTrackRefs, "AOD", "CFTRACKREF", o2::soa::Index<>, track::CollisionId, cftrackref::TrackId); //! Transient cf track index table

using CFTrackRef = CFTrackRefs::iterator;
//------

namespace cf2prongtrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(CFTrackProng0, cfTrackProng0, int, CFTracks, "_0"); //! Index to prong 1 CFTrack
DECLARE_SOA_INDEX_COLUMN_FULL(CFTrackProng1, cfTrackProng1, int, CFTracks, "_1"); //! Index to prong 2 CFTrack
// type
DECLARE_SOA_COLUMN(Pt, pt, float);       //! pT (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);     //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);     //! Phi angle
DECLARE_SOA_COLUMN(Mask, mask, uint8_t); //! Selection mask
const uint8_t kD0ToPiK = 0x1;
} // namespace cf2prongtrack
DECLARE_SOA_TABLE(CF2ProngTracks, "AOD", "CF2PRONGTRACK", //! Reduced track table
                  o2::soa::Index<>,
                  cftrack::CFCollisionId,
                  cf2prongtrack::CFTrackProng0Id,
                  cf2prongtrack::CFTrackProng1Id,
                  cf2prongtrack::Pt, cf2prongtrack::Eta, cf2prongtrack::Phi, cf2prongtrack::Mask);
using CF2ProngTrack = CF2ProngTracks::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_CFDERIVED_H
