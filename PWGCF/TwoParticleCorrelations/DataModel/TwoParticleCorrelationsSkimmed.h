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

/// \file TwoParticleCorrelationsSkimmed.h
/// \brief The skimmed datamodel
/// \author victor.gonzalez.sebastian@gmail.com

#ifndef PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_TWOPARTICLECORRELATIONSSKIMMED_H_
#define PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_TWOPARTICLECORRELATIONSSKIMMED_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>
#include <vector>

namespace o2::aod
{
/* we have to change from int to bool when bool columns work properly */
namespace cfskim
{
DECLARE_SOA_COLUMN(CentMult, centMult, std::vector<float>); //! The centrality/multiplicity pecentile
DECLARE_SOA_COLUMN(SelFlags, selFlags, uint64_t);           //! The skimming flags for collision selection
DECLARE_SOA_COLUMN(McSelFlags, mcSelFlags, uint64_t);       //! The skimming flags for generated collision selection
} // namespace cfskim
DECLARE_SOA_TABLE(CFCollisions, "AOD", "CFCOLLISION", //! Accepted reconstructed collisions/events filtered table
                  o2::soa::Index<>,
                  collision::PosZ,
                  bc::RunNumber,
                  timestamp::Timestamp,
                  cfskim::SelFlags,
                  cfskim::CentMult);
using CFCollision = CFCollisions::iterator;
DECLARE_SOA_TABLE(CFCollMasks, "AOD", "CFCOLLMASK", //! collision/event minimal filtered table
                  cfskim::SelFlags,
                  cfskim::CentMult);
using CFCollMask = CFCollMasks::iterator;
DECLARE_SOA_TABLE(CFMCCollisions, "AOD", "CFMCCOLLISION", //! Accepted generated collisions/events filtered table
                  o2::soa::Index<>,
                  mccollision::PosZ,
                  cfskim::McSelFlags,
                  cfskim::CentMult); //! TODO: To check if to have this columng does make sense for generated collisions
using CFMCCollision = CFMCCollisions::iterator;
DECLARE_SOA_TABLE(CFMCCollMasks, "AOD", "CFMCCOLLMASK", //! Generated collision/event minimal filtered table
                  cfskim::McSelFlags,
                  cfskim::CentMult); //! TODO: To check if to have this columng does make sense for generated collisions
using CFMCCollMask = CFMCCollMasks::iterator;
namespace cfskim
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);           //! Reconstructed collision/event
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);       //! Generated collision/event
DECLARE_SOA_COLUMN(TrackFlags, trackFlags, uint64_t);     //! The skimming flags for track selection, B0 track/particle positive charge, B1 track/particle negative charge
DECLARE_SOA_COLUMN(McTrackFlags, mcTrackFlags, uint64_t); //! The skimming flags for particle selection, B0 track/particle positive charge, B1 track/particle negative charge
DECLARE_SOA_COLUMN(PidFlags, pidFlags, uint64_t);         //! The PID skimming flags for track selection
DECLARE_SOA_COLUMN(Pt, pt, float);                        //! The track transverse momentum
DECLARE_SOA_COLUMN(Eta, eta, float);                      //! The track pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                      //! The track azimuthal angle
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign,                    //! Charge: positive: 1, negative: -1
                           [](uint64_t mask) -> int8_t { return ((mask & 0x1L) == 0x1L) ? 1 : ((mask & 0x2L) == 0x2L) ? -1
                                                                                                                      : 0; });

} // namespace cfskim
DECLARE_SOA_TABLE(CFTracks, "AOD", "CFTRACK", //! The reconstructed tracks filtered table
                  o2::soa::Index<>,
                  cfskim::CollisionId,
                  cfskim::TrackFlags,
                  cfskim::Pt,
                  cfskim::Eta,
                  cfskim::Phi,
                  cfskim::Sign<cfskim::TrackFlags>);
using CFTrack = CFTracks::iterator;
DECLARE_SOA_TABLE(CFTrackMasks, "AOD", "CFTRACKMASK", //! The reconstructed tracks filtered table
                  cfskim::TrackFlags);
using CFTrackMask = CFTrackMasks::iterator;
DECLARE_SOA_TABLE(CFTrackPIDs, "AOD", "CFTRACKPID", //! The reconstructed tracks PID filtered table
                  cfskim::PidFlags);
using CFTrackPID = CFTrackPIDs::iterator;
DECLARE_SOA_TABLE(CFMCParticles, "AOD", "CFMCPARICLE", //! The generated particles filtered table
                  o2::soa::Index<>,
                  cfskim::McCollisionId,
                  cfskim::McTrackFlags,
                  cfskim::Pt,
                  cfskim::Eta,
                  cfskim::Phi);
using CFMCParticle = CFMCParticles::iterator;
DECLARE_SOA_TABLE(CFMCPartMask, "AOD", "CFMCPARTMASK", //! The generated particles filtered table
                  cfskim::McTrackFlags);
using CFMCParticle = CFMCParticles::iterator;
namespace cfskim
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
}
DECLARE_SOA_TABLE(CFMCCollisionLbls, "AOD", "CFMCCOLLISONLBL",
                  cfskim::McCollisionId,
                  mccollisionlabel::McMask);
DECLARE_SOA_TABLE(CFMCTrackLabels, "AOD", "CFMCTRACKLABEL",
                  cfskim::McParticleId,
                  mctracklabel::McMask);
} // namespace o2::aod

#endif // PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_TWOPARTICLECORRELATIONSSKIMMED_H_
