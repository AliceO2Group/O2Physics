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

#include <Framework/ASoA.h>

namespace o2::aod
{
namespace hf_collisions_reduced
{
DECLARE_SOA_COLUMN(NumPvContrib, numPvContrib, int);   //! Event multiplicity from PV contributors
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Event multiplicity
DECLARE_SOA_COLUMN(Centrality, centrality, float);     //! Event centrality
DECLARE_SOA_COLUMN(PosZ, posZ, float);                 //! Primary vertex z position

} // namespace hf_collisions_reduced

DECLARE_SOA_TABLE(HfcRedCollisions, "AOD", "HFCREDCOLLISION", //! Table with collision info
                  soa::Index<>,
                  aod::hf_collisions_reduced::Multiplicity,
                  aod::hf_collisions_reduced::NumPvContrib,
                  aod::hf_collisions_reduced::PosZ);

DECLARE_SOA_TABLE(HfcRedCorrColls, "AOD", "HFCREDCORRCOLL", //! Table with collision info
                  soa::Index<>,
                  aod::hf_collisions_reduced::Multiplicity,
                  aod::hf_collisions_reduced::NumPvContrib,
                  aod::hf_collisions_reduced::Centrality,
                  aod::hf_collisions_reduced::PosZ);

using HfcRedCollision = HfcRedCollisions::iterator;

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
DECLARE_SOA_COLUMN(BdtScoreBkg, bdtScoreBkg, float);        //! BDT output score for background hypothesis
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
DECLARE_SOA_COLUMN(NSigmaTpc, nSigmaTpc, float);           //! Number of sigma TPC
DECLARE_SOA_COLUMN(NSigmaTof, nSigmaTpc, float);           //! Number of sigma TOF
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
                  aod::hf_assoc_track_reduced::DcaZ);

DECLARE_SOA_TABLE(AssocTrackPids, "AOD", "ASSOCTRACKPID", //! Table with associated track pid info
                  soa::Index<>,
                  aod::hf_assoc_track_reduced::NSigmaTpc,
                  aod::hf_assoc_track_reduced::NSigmaTof);

// definition of columns and tables for Charm-Hadron and Hadron-Hadron correlation pairs
namespace hf_correl_charm_had_reduced
{
// Correlation columns
DECLARE_SOA_INDEX_COLUMN(HfcRedCorrColl, hfcRedCorrColl); //! ReducedCollision index
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);            //! DeltaPhi between charm hadron and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);            //! DeltaEta between charm hadron and Hadrons
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);                //! Pool Bin for the MixedEvent
// General trigger particle columns
DECLARE_SOA_COLUMN(PhiTrig, phiTrig, float); //! Phi of the trigger candidate
DECLARE_SOA_COLUMN(EtaTrig, etaTrig, float); //! Eta of the trigger candidate
DECLARE_SOA_COLUMN(PtTrig, ptTrig, float);   //! Pt of the trigger candidate
// Charm trigger particle selection columns
DECLARE_SOA_COLUMN(InvMassTrig, invMassTrig, float);     //! Invariant mass of Charm trigger candidate
DECLARE_SOA_COLUMN(BdtScore0Trig, bdtScore0Trig, float); //! First BDT output score
DECLARE_SOA_COLUMN(BdtScore1Trig, bdtScore1Trig, float); //! Second BDT output score
// Hadron trigger particle selection columns
DECLARE_SOA_COLUMN(NTpcCrossedRowsTrig, nTpcCrossedRowsTrig, int); //! Number of crossed TPC Rows
DECLARE_SOA_COLUMN(ItsClsMapTrig, itsClsMapTrig, int);             //! ITS cluster map, one bit per a layer, starting from the innermost
DECLARE_SOA_COLUMN(ItsNClsTrig, itsNClsTrig, int);                 //! Number of ITS clusters
DECLARE_SOA_COLUMN(DcaXYTrig, dcaXYTrig, float);                   //! Impact parameter in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(DcaZTrig, dcaZTrig, float);                     //! Impact parameter in Z of the track to the primary vertex
// General associated particle columns
DECLARE_SOA_COLUMN(EtaAssoc, etaAssoc, float); //! Eta of the associated candidate
DECLARE_SOA_COLUMN(PhiAssoc, phiAssoc, float); //! Phi of the associated candidate
DECLARE_SOA_COLUMN(PtAssoc, ptAssoc, float);   //! Pt of the associated candidate
// Hadron associated particle selection columns
DECLARE_SOA_COLUMN(NTpcCrossedRowsAssoc, nTpcCrossedRowsAssoc, int); //! Number of crossed TPC Rows
DECLARE_SOA_COLUMN(ItsClsMapAssoc, itsClsMapAssoc, int);             //! ITS cluster map, one bit per a layer, starting from the innermost
DECLARE_SOA_COLUMN(ItsNClsAssoc, itsNClsAssoc, int);                 //! Number of ITS clusters
DECLARE_SOA_COLUMN(DcaXYAssoc, dcaXYAssoc, float);                   //! Impact parameter in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(DcaZAssoc, dcaZAssoc, float);                     //! Impact parameter in Z of the track to the primary vertex
} // namespace hf_correl_charm_had_reduced

DECLARE_SOA_TABLE(HfcRedTrigBases, "AOD", "HFCREDTRIGBASE", //! Table with trigger candidate base info
                  soa::Index<>,
                  aod::hf_correl_charm_had_reduced::PhiTrig,
                  aod::hf_correl_charm_had_reduced::EtaTrig);

DECLARE_SOA_TABLE(HfcRedTrigCharms, "AOD", "HFCREDTRIGCHARM", //! Table with Same Event Charm-Hadron pairs information
                  aod::hf_correl_charm_had_reduced::HfcRedCorrCollId,
                  aod::hf_correl_charm_had_reduced::PtTrig,
                  aod::hf_correl_charm_had_reduced::InvMassTrig,
                  aod::hf_correl_charm_had_reduced::BdtScore0Trig,
                  aod::hf_correl_charm_had_reduced::BdtScore1Trig);

DECLARE_SOA_TABLE(HfcRedTrigTracks, "AOD", "HFCREDTRIGTRACK", //! Table with Same Event Charm-Hadron pairs information
                  aod::hf_correl_charm_had_reduced::HfcRedCorrCollId,
                  aod::hf_correl_charm_had_reduced::PtTrig,
                  aod::hf_correl_charm_had_reduced::NTpcCrossedRowsTrig,
                  aod::hf_correl_charm_had_reduced::ItsClsMapTrig,
                  aod::hf_correl_charm_had_reduced::ItsNClsTrig,
                  aod::hf_correl_charm_had_reduced::DcaXYTrig,
                  aod::hf_correl_charm_had_reduced::DcaZTrig);

namespace hf_correl_charm_had_reduced
{
DECLARE_SOA_INDEX_COLUMN(HfcRedTrigCharm, hfcRedTrigCharm); //! Same Event pair index
DECLARE_SOA_INDEX_COLUMN(HfcRedTrigTrack, hfcRedTrigTrack); //! Same Event pair index
} // namespace hf_correl_charm_had_reduced

DECLARE_SOA_TABLE(HfcRedSEChBases, "AOD", "HFCREDSECHBASE", //! Table with Same Event Trig-Assoc pairs
                  aod::hf_correl_charm_had_reduced::HfcRedCorrCollId,
                  aod::hf_correl_charm_had_reduced::HfcRedTrigCharmId,
                  aod::hf_correl_charm_had_reduced::PtAssoc,
                  aod::hf_correl_charm_had_reduced::DeltaEta,
                  aod::hf_correl_charm_had_reduced::DeltaPhi);

DECLARE_SOA_TABLE(HfcRedSEHadBases, "AOD", "HFCREDSEHADBASE", //! Table with Same Event Trig-Assoc pairs
                  aod::hf_correl_charm_had_reduced::HfcRedCorrCollId,
                  aod::hf_correl_charm_had_reduced::HfcRedTrigTrackId,
                  aod::hf_correl_charm_had_reduced::PtAssoc,
                  aod::hf_correl_charm_had_reduced::DeltaEta,
                  aod::hf_correl_charm_had_reduced::DeltaPhi);

DECLARE_SOA_TABLE(HfcRedAssBases, "AOD", "HFCREDASSBASE", //! Table with associated candidate base info
                  soa::Index<>,
                  aod::hf_correl_charm_had_reduced::HfcRedCorrCollId,
                  aod::hf_correl_charm_had_reduced::PhiAssoc,
                  aod::hf_correl_charm_had_reduced::EtaAssoc,
                  aod::hf_correl_charm_had_reduced::PtAssoc);

DECLARE_SOA_TABLE(HfcRedAssTracks, "AOD", "HFCREDASSTRACK", //! Table with Same Event Track Selections information
                  aod::hf_correl_charm_had_reduced::NTpcCrossedRowsAssoc,
                  aod::hf_correl_charm_had_reduced::ItsClsMapAssoc,
                  aod::hf_correl_charm_had_reduced::ItsNClsAssoc,
                  aod::hf_correl_charm_had_reduced::DcaXYAssoc,
                  aod::hf_correl_charm_had_reduced::DcaZAssoc);

DECLARE_SOA_TABLE(HfcRedSEChHads, "AOD", "HFCREDSECHHAD", //! Correlation pairs information Same Event
                  aod::hf_correl_charm_had_reduced::PoolBin,
                  aod::hf_correl_charm_had_reduced::PtTrig,
                  aod::hf_correl_charm_had_reduced::PtAssoc,
                  aod::hf_correl_charm_had_reduced::DeltaEta,
                  aod::hf_correl_charm_had_reduced::DeltaPhi,
                  aod::hf_correl_charm_had_reduced::InvMassTrig,
                  aod::hf_correl_charm_had_reduced::BdtScore0Trig,
                  aod::hf_correl_charm_had_reduced::BdtScore1Trig,
                  aod::hf_correl_charm_had_reduced::NTpcCrossedRowsAssoc,
                  aod::hf_correl_charm_had_reduced::ItsClsMapAssoc,
                  aod::hf_correl_charm_had_reduced::ItsNClsAssoc,
                  aod::hf_correl_charm_had_reduced::DcaXYAssoc,
                  aod::hf_correl_charm_had_reduced::DcaZAssoc,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(HfcRedMEChHads, "AOD", "HFCREDMECHHAD", //! Correlation pairs information Same Event
                  aod::hf_correl_charm_had_reduced::PoolBin,
                  aod::hf_correl_charm_had_reduced::PtTrig,
                  aod::hf_correl_charm_had_reduced::PtAssoc,
                  aod::hf_correl_charm_had_reduced::DeltaEta,
                  aod::hf_correl_charm_had_reduced::DeltaPhi,
                  aod::hf_correl_charm_had_reduced::InvMassTrig,
                  aod::hf_correl_charm_had_reduced::BdtScore0Trig,
                  aod::hf_correl_charm_had_reduced::BdtScore1Trig,
                  aod::hf_correl_charm_had_reduced::NTpcCrossedRowsAssoc,
                  aod::hf_correl_charm_had_reduced::ItsClsMapAssoc,
                  aod::hf_correl_charm_had_reduced::ItsNClsAssoc,
                  aod::hf_correl_charm_had_reduced::DcaXYAssoc,
                  aod::hf_correl_charm_had_reduced::DcaZAssoc,
                  soa::Marker<2>);

DECLARE_SOA_TABLE(HfcRedSEHadHads, "AOD", "HFCREDSEHADHAD", //! Correlation pairs information Same Event
                  aod::hf_correl_charm_had_reduced::PoolBin,
                  aod::hf_correl_charm_had_reduced::PtTrig,
                  aod::hf_correl_charm_had_reduced::PtAssoc,
                  aod::hf_correl_charm_had_reduced::DeltaEta,
                  aod::hf_correl_charm_had_reduced::DeltaPhi,
                  aod::hf_correl_charm_had_reduced::NTpcCrossedRowsTrig,
                  aod::hf_correl_charm_had_reduced::ItsClsMapTrig,
                  aod::hf_correl_charm_had_reduced::ItsNClsTrig,
                  aod::hf_correl_charm_had_reduced::DcaXYTrig,
                  aod::hf_correl_charm_had_reduced::DcaZTrig,
                  aod::hf_correl_charm_had_reduced::NTpcCrossedRowsAssoc,
                  aod::hf_correl_charm_had_reduced::ItsClsMapAssoc,
                  aod::hf_correl_charm_had_reduced::ItsNClsAssoc,
                  aod::hf_correl_charm_had_reduced::DcaXYAssoc,
                  aod::hf_correl_charm_had_reduced::DcaZAssoc,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(HfcRedMEHadHads, "AOD", "HFCREDMEHADHAD", //! Correlation pairs information Same Event
                  aod::hf_correl_charm_had_reduced::PoolBin,
                  aod::hf_correl_charm_had_reduced::PtTrig,
                  aod::hf_correl_charm_had_reduced::PtAssoc,
                  aod::hf_correl_charm_had_reduced::DeltaEta,
                  aod::hf_correl_charm_had_reduced::DeltaPhi,
                  aod::hf_correl_charm_had_reduced::NTpcCrossedRowsTrig,
                  aod::hf_correl_charm_had_reduced::ItsClsMapTrig,
                  aod::hf_correl_charm_had_reduced::ItsNClsTrig,
                  aod::hf_correl_charm_had_reduced::DcaXYTrig,
                  aod::hf_correl_charm_had_reduced::DcaZTrig,
                  aod::hf_correl_charm_had_reduced::NTpcCrossedRowsAssoc,
                  aod::hf_correl_charm_had_reduced::ItsClsMapAssoc,
                  aod::hf_correl_charm_had_reduced::ItsNClsAssoc,
                  aod::hf_correl_charm_had_reduced::DcaXYAssoc,
                  aod::hf_correl_charm_had_reduced::DcaZAssoc,
                  soa::Marker<2>);

DECLARE_SOA_TABLE(HfcRedCollInfos, "AOD", "HFCREDCOLLINFO", //! Table with collision info
                  aod::hf_collisions_reduced::Multiplicity,
                  aod::hf_collisions_reduced::NumPvContrib,
                  aod::hf_collisions_reduced::Centrality);
} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_DERIVEDDATACORRELATIONTABLES_H_
