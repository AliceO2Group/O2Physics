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
#ifndef O2_ANALYSIS_TWOPSKIMMED_H
#define O2_ANALYSIS_TWOPSKIMMED_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2
{
namespace aod
{
/* we have to change from int to bool when bool columns work properly */
namespace twopskim
{
DECLARE_SOA_COLUMN(TwoPSkimmedCollisionCentMult, centmult, float); //! The centrality/multiplicity pecentile
DECLARE_SOA_COLUMN(TwoPSkimmedCollisionFlags, selflags, uint64_t); //! The skimming flags for collision selection
} // namespace twopskim
DECLARE_SOA_TABLE(TwoPSkimmedCollisions, "AOD", "TWOPSKMDCOLL", //! Accepted reconstructed collisions/events filtered table
                  o2::soa::Index<>,
                  collision::PosZ,
                  twopskim::TwoPSkimmedCollisionCentMult,
                  twopskim::TwoPSkimmedCollisionFlags);
using TwoPSkimmedCollision = TwoPSkimmedCollisions::iterator;
DECLARE_SOA_TABLE(TwoPSkimmedGenCollisions, "AOD", "TWOPSKMDGENCOLL", //! Accepted generated collisions/events filtered table
                  o2::soa::Index<>,
                  mccollision::PosZ,
                  twopskim::TwoPSkimmedCollisionCentMult,
                  twopskim::TwoPSkimmedCollisionFlags);
using TwoPSkimmedGenCollision = TwoPSkimmedGenCollisions::iterator;
namespace twopskim
{
DECLARE_SOA_INDEX_COLUMN(TwoPSkimmedCollision, event);           //! Reconstructed collision/event
DECLARE_SOA_INDEX_COLUMN(TwoPSkimmedGenCollision, mcevent);      //! Generated collision/event
DECLARE_SOA_COLUMN(TwoPSkimmedTrackFlags, trackflags, uint64_t); //! The skimming flags for track selection
DECLARE_SOA_COLUMN(sPt, spt, float);                             //! The track charge signed transverse momentum
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,                               //! The track transverse momentum
                           [](float signedpt) -> float {
                             return std::abs(signedpt);
                           });
DECLARE_SOA_COLUMN(Eta, eta, float);        //! The track pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);        //! The track azimuthal angle
DECLARE_SOA_COLUMN(Charge, charge, int8_t); //! Charge: in units of electric charge
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign,      //! The track transverse momentum
                           [](float signedpt) -> short {
                             return (signedpt >= 0) ? 1 : -1;
                           });
} // namespace twopskim
DECLARE_SOA_TABLE(TwoPSkimmedTracks, "AOD", "TWOPSKMDTRKS", //! The reconstructed tracks filtered table
                  twopskim::TwoPSkimmedCollisionId,
                  twopskim::TwoPSkimmedTrackFlags,
                  twopskim::sPt,
                  twopskim::Eta,
                  twopskim::Phi,
                  twopskim::Pt<twopskim::sPt>,
                  twopskim::Sign<twopskim::sPt>);
DECLARE_SOA_TABLE(TwoPSkimmedParticles, "AOD", "TWOPSKMDPARTS", //! The generated particles filtered table
                  twopskim::TwoPSkimmedGenCollisionId,
                  twopskim::TwoPSkimmedTrackFlags,
                  twopskim::sPt,
                  twopskim::Eta,
                  twopskim::Phi,
                  twopskim::Pt<twopskim::sPt>,
                  twopskim::Sign<twopskim::sPt>);
} // namespace aod
} // namespace o2

#endif // O2_ANALYSIS_TWOPSKIMMED_H
