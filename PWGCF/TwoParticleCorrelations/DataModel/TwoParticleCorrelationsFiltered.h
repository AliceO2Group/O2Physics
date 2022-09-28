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
#ifndef O2_ANALYSIS_TWOPARTFILTERED_H
#define O2_ANALYSIS_TWOPARTFILTERED_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "TwoParticleCorrelationsSkimmed.h"

namespace o2
{
namespace aod
{
/* we have to change from int to bool when bool columns work properly */
namespace twopfilter
{
DECLARE_SOA_COLUMN(TwoPCollisionCentMult, centmult, float); //! The centrality/multiplicity pecentile
} // namespace twopfilter
DECLARE_SOA_TABLE(TwoPAcceptedCollisions, "AOD", "TWOPACCCOLL", //! Accepted reconstructed collisions/events filtered table
                  o2::soa::Index<>,
                  collision::PosZ,
                  twopfilter::TwoPCollisionCentMult);
using TowPAcceptedCollision = TwoPAcceptedCollisions::iterator;
DECLARE_SOA_TABLE(TwoPAcceptedGenCollisions, "AOD", "TWOPACCGENCOLL", //! Accepted generated collisions/events filtered table
                  o2::soa::Index<>,
                  mccollision::PosZ,
                  twopfilter::TwoPCollisionCentMult);
using TwoPAcceptedGenCollision = TwoPAcceptedGenCollisions::iterator;
namespace twopfilter
{
DECLARE_SOA_INDEX_COLUMN(TwoPAcceptedCollision, collision);            //! Reconstructed collision/event
DECLARE_SOA_INDEX_COLUMN(TwoPAcceptedGenCollision, mccollision);       //! Generated collision/event
DECLARE_SOA_COLUMN(TwoPTrackacceptedAs, twoptrackacceptedas, uint8_t); //! Track accepted as type 0..255, even positive or particle, odd negative or antiparticle
} // namespace twopfilter
DECLARE_SOA_TABLE(TwoPFilteredTracks, "AOD", "FILTEREDTRKS", //! The reconstructed tracks filtered table
                  twopfilter::TwoPAcceptedCollisionId,
                  twopfilter::TwoPTrackacceptedAs,
                  cfskim::Pt,
                  cfskim::Eta,
                  cfskim::Phi);
DECLARE_SOA_TABLE(TwoPFilteredParticles, "AOD", "FILTEREDGENTRKS", //! The generated particles filtered table
                  twopfilter::TwoPAcceptedGenCollisionId,
                  twopfilter::TwoPTrackacceptedAs,
                  cfskim::Pt,
                  cfskim::Eta,
                  cfskim::Phi);
} // namespace aod
} // namespace o2

#endif // O2_ANALYSIS_TWOPARTFILTERED_H
