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

/// \file DerivedDataCorrelationTables.h
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
DECLARE_SOA_COLUMN(Prong0Id, prong0Id, int);                //! Prong 0 index
DECLARE_SOA_COLUMN(Prong1Id, prong1Id, int);                //! Prong 1 index
DECLARE_SOA_COLUMN(Prong2Id, prong2Id, int);                //! Prong2 index
DECLARE_SOA_COLUMN(PhiCand, phiCand, float);                //! Phi of the candidate
DECLARE_SOA_COLUMN(EtaCand, etaCand, float);                //! Eta of the candidate
DECLARE_SOA_COLUMN(PtCand, ptCand, float);                  //! Pt of the candidate
DECLARE_SOA_COLUMN(InvMassDs, invMassDs, float);            //! Invariant mass of Ds candidate
DECLARE_SOA_COLUMN(BdtScorePrompt, bdtScorePrompt, float);  //! BDT output score for prompt hypothesis
DECLARE_SOA_COLUMN(BdtScoreBkg, bdtScoreBkg, float);        //! BDT output score for backgronud hypothesis
} // namespace hf_candidate_reduced
DECLARE_SOA_TABLE(DsCandReduceds, "AOD", "DSCANDREDUCED", //! Table with Ds candidate info
                  soa::Index<>,
                  aod::hf_candidate_reduced::HfcRedCollisionId,
                  aod::hf_candidate_reduced::PhiCand,
                  aod::hf_candidate_reduced::EtaCand,
                  aod::hf_candidate_reduced::PtCand,
                  aod::hf_candidate_reduced::InvMassDs,
                  aod::hf_candidate_reduced::Prong0Id,
                  aod::hf_candidate_reduced::Prong1Id,
                  aod::hf_candidate_reduced::Prong2Id);

DECLARE_SOA_TABLE(DsCandSelInfos, "AOD", "DSCANDSELINFO", //! Table with Ds candidate selection info
                  soa::Index<>,
                  aod::hf_candidate_reduced::HfcRedCollisionId,
                  aod::hf_candidate_reduced::BdtScorePrompt,
                  aod::hf_candidate_reduced::BdtScoreBkg);

namespace hf_assoc_track_reduced
{
DECLARE_SOA_COLUMN(OriginTrackId, originTrackId, int);     //! Original track index
DECLARE_SOA_COLUMN(NTpcCrossedRows, nTpcCrossedRows, int); //! Number of crossed TPC Rows
DECLARE_SOA_COLUMN(ItsClusterMap, itsClusterMap, int);     //! ITS cluster map, one bit per a layer, starting from the innermost
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, int);                 //! Number of ITS clusters
DECLARE_SOA_COLUMN(EtaAssocTrack, etaAssocTrack, float);   //! Eta of the track
DECLARE_SOA_COLUMN(PhiAssocTrack, phiAssocTrack, float);   //! Phi of the track
DECLARE_SOA_COLUMN(PtAssocTrack, ptAssocTrack, float);     //! Pt of the track
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                   //! Impact parameter in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                     //! Impact parameter in Z of the track to the primary vertex
} // namespace hf_assoc_track_reduced
DECLARE_SOA_TABLE(AssocTrackReds, "AOD", "ASSOCTRACKRED", //! Table with associated track info
                  soa::Index<>,
                  aod::hf_candidate_reduced::HfcRedCollisionId,
                  aod::hf_assoc_track_reduced::OriginTrackId,
                  aod::hf_assoc_track_reduced::PhiAssocTrack,
                  aod::hf_assoc_track_reduced::EtaAssocTrack,
                  aod::hf_assoc_track_reduced::PtAssocTrack);

DECLARE_SOA_TABLE(AssocTrackSels, "AOD", "ASSOCTRACKSEL", //! Table with associated track info
                  soa::Index<>,
                  aod::hf_candidate_reduced::HfcRedCollisionId,
                  aod::hf_assoc_track_reduced::NTpcCrossedRows,
                  aod::hf_assoc_track_reduced::ItsClusterMap,
                  aod::hf_assoc_track_reduced::ItsNCls,
                  aod::hf_assoc_track_reduced::DcaXY,
                  aod::hf_assoc_track_reduced::DcaZ)
} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_DERIVEDDATACORRELATIONTABLES_H_
