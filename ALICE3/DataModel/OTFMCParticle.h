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

///
/// \file   OTFMCParticle.h
/// \author Jesper Karlsson Gumprecht
/// \since  16/12/2025
/// \brief  Redefinition of the mcparticles table specifically for the fast sim
///

#ifndef ALICE3_DATAMODEL_OTFMCPARTICLE_H_
#define ALICE3_DATAMODEL_OTFMCPARTICLE_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{

namespace otfmcparticle
{
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(IsAlive, isAlive, bool);
DECLARE_SOA_COLUMN(IsPrimary, isPrimary, bool);

DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother0, mother0, int, "McPartsWithDau_Mother0");       //! Track index of the first mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother1, mother1, int, "McPartsWithDau_Mother1");       //! Track index of the last mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter0, daughter0, int, "McPartsWithDau_Daughter0"); //! Track index of the first daugther
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter1, daughter1, int, "McPartsWithDau_Daughter1"); //! Track index of the last daugther
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);                                     //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters);                                 //! Daughter tracks (possibly empty) slice. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
} // namespace otfmcparticle

DECLARE_SOA_TABLE_FULL(McPartWithDaus, "McPartWithDaus", "AOD", "MCPARTSWITHDAU",
                       o2::soa::Index<>,
                       mcparticle::McCollisionId,
                       mcparticle::PdgCode,
                       mcparticle::StatusCode,
                       mcparticle::Flags,
                       otfmcparticle::MothersIds,
                       otfmcparticle::DaughtersIdSlice,
                       mcparticle::Weight,
                       mcparticle::Px,
                       mcparticle::Py,
                       mcparticle::Pz,
                       mcparticle::E,
                       mcparticle::Vx,
                       mcparticle::Vy,
                       mcparticle::Vz,
                       mcparticle::Vt,
                       otfmcparticle::Phi,
                       otfmcparticle::Eta,
                       otfmcparticle::Pt,
                       otfmcparticle::P,
                       otfmcparticle::Y,
                       otfmcparticle::IsAlive,
                       otfmcparticle::IsPrimary,
                       mcparticle::PVector<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::GetHepMCStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using McPartWithDau = McPartWithDaus::iterator;

namespace otfmctracklable
{
DECLARE_SOA_INDEX_COLUMN(McPartWithDau, mcPartWithDau); //! MC particle
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);           //! Bit mask to indicate detector mismatches (bit ON means mismatch). Bit 0-6: mismatch at ITS layer. Bit 12: ITSAB tracklet mismatch. Bit 13: ITS-TPC mismatch. Bit 14: isNoise == True (global track), Bit 15: isFake == True (global track)
} // namespace otfmctracklable

DECLARE_SOA_TABLE(McTrackWithDauLabels, "AOD", "MCTRACKWithDAULABEL", //! Table joined to the track table containing the MC index
                  otfmctracklable::McPartWithDauId, otfmctracklable::McMask);

using McTrackWithDauLabel = McTrackWithDauLabels::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFMCPARTICLE_H_
