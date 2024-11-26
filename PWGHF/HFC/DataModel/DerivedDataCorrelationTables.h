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

/// \file DerivedDataCorrelations.h
/// \brief Tables for producing derived data for correlation analysis
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#ifndef PWGHF_HFC_DATAMODEL_DERIVEDDATACORRELATIONTABLES_H_
#define PWGHF_HFC_DATAMODEL_DERIVEDDATACORRELATIONTABLES_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace hf_collisions_reduced
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Event multiplicity
DECLARE_SOA_COLUMN(PosZ, posZ, float);                 //! Primary vertex z position

} // namespace hf_collisions_reduced

DECLARE_SOA_TABLE(HfcRedCollisions, "AOD", "HFCREDCOLLISION", //! Table with collision info
                  soa::Index<>,
                  aod::hf_collisions_reduced::Multiplicity,
                  aod::hf_collisions_reduced::PosZ);

using HfcRedCollision = HfcRedCollisions::iterator;

// DECLARE_SOA_TABLE(HfCandColCounts, "AOD", "HFCANDCOLCOUNT", //! Table with number of collisions which contain at least one candidate
//                   aod::hf_collisions_reduced::OriginalCollisionCount);

namespace hf_candidate_reduced
{
DECLARE_SOA_INDEX_COLUMN(HfcRedCollision, hfcRedCollision); //! ReducedCollision index
DECLARE_SOA_COLUMN(PhiCand, phiCand, float);              //! Phi of the candidate
DECLARE_SOA_COLUMN(EtaCand, etaCand, float);              //! Eta of the candidate
DECLARE_SOA_COLUMN(PtCand, ptCand, float);                //! Pt of the candidate
DECLARE_SOA_COLUMN(InvMassDs, invMassDs, float);          //! Invariant mass of Ds candidate
} // namespace hf_candidate_reduced
DECLARE_SOA_TABLE(DsCandReduceds, "AOD", "DSCANDREDUCED", //! Table with Ds candidate info (rectangular selection)
                  soa::Index<>,
                  aod::hf_candidate_reduced::HfcRedCollisionId,
                  aod::hf_candidate_reduced::PhiCand,
                  aod::hf_candidate_reduced::EtaCand,
                  aod::hf_candidate_reduced::PtCand,
                  aod::hf_candidate_reduced::InvMassDs);

namespace hf_assoc_track_reduced
{
DECLARE_SOA_COLUMN(TrackId, trackId, int);               //! Original track index
DECLARE_SOA_COLUMN(EtaAssocTrack, etaAssocTrack, float); //! Eta of the track
DECLARE_SOA_COLUMN(PhiAssocTrack, phiAssocTrack, float); //! Phi of the track
DECLARE_SOA_COLUMN(PtAssocTrack, ptAssocTrack, float);   //! Pt of the track
} // namespace hf_assoc_track_reduced
DECLARE_SOA_TABLE(AssocTrackReds, "AOD", "ASSOCTRACKRED", //! Table with associated track info
                  soa::Index<>,
                  aod::hf_candidate_reduced::HfcRedCollisionId,
                  aod::hf_assoc_track_reduced::PhiAssocTrack,
                  aod::hf_assoc_track_reduced::EtaAssocTrack,
                  aod::hf_assoc_track_reduced::PtAssocTrack)
} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_DERIVEDDATACORRELATIONTABLES_H_
