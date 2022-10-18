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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "PWGCF/TwoParticleCorrelations/Core/SkimmingConfigurableCuts.h"
#include "PWGCF/TwoParticleCorrelations/Core/EventSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/TrackSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/PIDSelectionFilterAndAnalysis.h"
#include "Framework/runDataProcessing.h"

/* TODO: this will go to its proper header at some point */
#ifndef O2_ANALYSIS_cfskimMED_H
#define O2_ANALYSIS_cfskimMED_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2
{
namespace aod
{
/* we have to change from int to bool when bool columns work properly */
namespace cfskim
{
DECLARE_SOA_COLUMN(CFCollisionCentMult, centmult, std::vector<float>); //! The centrality/multiplicity pecentile
DECLARE_SOA_COLUMN(CFCollisionFlags, selflags, uint64_t);              //! The skimming flags for collision selection
} // namespace cfskim
DECLARE_SOA_TABLE(CFCollisions, "AOD", "CFCOLLISION", //! Accepted reconstructed collisions/events filtered table
                  o2::soa::Index<>,
                  collision::PosZ,
                  bc::RunNumber,
                  timestamp::Timestamp,
                  cfskim::CFCollisionFlags,
                  cfskim::CFCollisionCentMult);
using CFCollision = CFCollisions::iterator;
DECLARE_SOA_TABLE(CFMCCollisions, "AOD", "CFMCCOLLISION", //! Accepted generated collisions/events filtered table
                  o2::soa::Index<>,
                  mccollision::PosZ,
                  cfskim::CFCollisionFlags,
                  cfskim::CFCollisionCentMult); //! To check if to have this columng does make sense for generated collisions
using CFMCCollision = CFMCCollisions::iterator;
namespace cfskim
{
DECLARE_SOA_INDEX_COLUMN(CFCollision, cfcollision);     //! Reconstructed collision/event
DECLARE_SOA_INDEX_COLUMN(CFMCCollision, cfmccollision); //! Generated collision/event
DECLARE_SOA_COLUMN(CFTrackFlags, trackflags, uint64_t); //! The skimming flags for track selection, B0 track/particle positive charge, B1 track/particle negative charge
DECLARE_SOA_COLUMN(CFPidFlags, pidflags, uint64_t);     //! The PID skimming flags for track selection
DECLARE_SOA_COLUMN(Pt, pt, float);                      //! The track transverse momentum
DECLARE_SOA_COLUMN(Eta, eta, float);                    //! The track pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                    //! The track azimuthal angle
} // namespace cfskim
DECLARE_SOA_TABLE(CFTracks, "AOD", "CFTRACK", //! The reconstructed tracks filtered table
                  cfskim::CFCollisionId,
                  cfskim::CFTrackFlags,
                  cfskim::Pt,
                  cfskim::Eta,
                  cfskim::Phi);
using CFTrack = CFTracks::iterator;
DECLARE_SOA_TABLE(CFTrackPIDs, "AOD", "CFTRACKPID", //! The reconstructed tracks PID filtered table
                  cfskim::CFPidFlags);
using CFTrackPID = CFTrackPIDs::iterator;
DECLARE_SOA_TABLE(CFMCParticles, "AOD", "CFMCPARICLE", //! The generated particles filtered table
                  o2::soa::Index<>,
                  cfskim::CFMCCollisionId,
                  cfskim::CFTrackFlags,
                  cfskim::Pt,
                  cfskim::Eta,
                  cfskim::Phi);
using CFMCParticle = CFMCParticles::iterator;
namespace cfskim
{
DECLARE_SOA_INDEX_COLUMN(CFMCParticle, mcparticle);
}
DECLARE_SOA_TABLE(CFMCCollisionLbls, "AOD", "CFMCCOLLISONLBL",
                  cfskim::CFMCCollisionId,
                  mccollisionlabel::McMask);
DECLARE_SOA_TABLE(CFMCTrackLabels, "AOD", "CFMCTRACKLABEL",
                  cfskim::CFMCParticleId,
                  mctracklabel::McMask);
} // namespace aod
} // namespace o2

#endif // O2_ANALYSIS_cfskimMED_H
