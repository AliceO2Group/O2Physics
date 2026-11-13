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

/// \file TwoParticleCorrelationsFiltered.h
/// \brief The filtered datamodel
/// \author victor.gonzalez.sebastian@gmail.com

#ifndef PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_TWOPARTICLECORRELATIONSFILTERED_H_
#define PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_TWOPARTICLECORRELATIONSFILTERED_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
/* we have to change from int to bool when bool columns work properly */
namespace twopfilter
{
DECLARE_SOA_COLUMN(Accepted, accepted, uint8_t);   //! If the collision/event has been accepted or not
DECLARE_SOA_COLUMN(CentMult, centMult, float);     //! The centrality/multiplicity pecentile
DECLARE_SOA_COLUMN(McCentMult, mcCentMult, float); //! The centrality/multiplicity pecentile
} // namespace twopfilter
DECLARE_SOA_TABLE(TwoPAcceptedCollisions, "AOD", "TWOPACCCOLL", //! Accepted reconstructed collisions/events filtered table
                  twopfilter::CentMult,
                  twopfilter::Accepted);
using TowPAcceptedCollision = TwoPAcceptedCollisions::iterator;
DECLARE_SOA_TABLE(TwoPAcceptedGenCollisions, "AOD", "TWOPACCGENCOLL", //! Accepted generated collisions/events filtered table
                  twopfilter::McCentMult,
                  twopfilter::Accepted);
using TwoPAcceptedGenCollision = TwoPAcceptedGenCollisions::iterator;
namespace twopfilter
{
DECLARE_SOA_COLUMN(AcceptedAs, acceptedAs, uint8_t);     //! Track accepted as type 0..255, even positive, odd negative, < 0 not accepted
DECLARE_SOA_COLUMN(McAcceptedAs, mcAcceptedAs, uint8_t); //! Particle accepted as type 0..255, even positive, odd negative, < 0 not accepted
} // namespace twopfilter
DECLARE_SOA_TABLE(TwoPFilteredTracks, "AOD", "FILTEREDTRKS", //! The reconstructed tracks filtered table
                  twopfilter::AcceptedAs);
DECLARE_SOA_TABLE(TwoPFilteredParticles, "AOD", "FILTEREDGENTRKS", //! The generated particles filtered table
                  twopfilter::McAcceptedAs);
} // namespace o2::aod

#endif // PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_TWOPARTICLECORRELATIONSFILTERED_H_
