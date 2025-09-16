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
using HfcRedCorrColl = HfcRedCorrColls::iterator;

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

namespace hf_correlation_trigger_reduced
{
DECLARE_SOA_INDEX_COLUMN(HfcRedCorrColl, hfcRedCorrColl);           //! ReducedCollision index
DECLARE_SOA_COLUMN(PhiTrig, phiTrig, float);                        //! Phi of the trigger candidate
DECLARE_SOA_COLUMN(EtaTrig, etaTrig, float);                        //! Eta of the trigger candidate
DECLARE_SOA_COLUMN(PtTrig, ptTrig, float);                          //! Pt of the trigger candidate
DECLARE_SOA_COLUMN(InvMassTrig, invMassTrig, float);                //! Invariant mass of Charm trigger candidate
DECLARE_SOA_COLUMN(BdtScore0Trig, bdtScore0Trig, float);            //! First BDT output score
DECLARE_SOA_COLUMN(BdtScore1Trig, bdtScore1Trig, float);            //! Second BDT output score
DECLARE_SOA_COLUMN(NTpcCrossedRowsTrig, nTpcCrossedRowsTrig, int);  //! Number of crossed TPC Rows
DECLARE_SOA_COLUMN(ItsClusterMapTrig, itsClusterMapTrig, int);      //! ITS cluster map, one bit per a layer, starting from the innermost
DECLARE_SOA_COLUMN(ItsNClsTrig, itsNClsTrig, int);                  //! Number of ITS clusters
DECLARE_SOA_COLUMN(EtaTrigTrack, etaTrigTrack, float);              //! Eta of the track
DECLARE_SOA_COLUMN(PhiTrigTrack, phiTrigTrack, float);              //! Phi of the track
DECLARE_SOA_COLUMN(PtTrigTrack, ptTrigTrack, float);                //! Pt of the track
DECLARE_SOA_COLUMN(DcaXYTrig, dcaXYTrig, float);                    //! Impact parameter in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(DcaZTrig, dcaZTrig, float);                      //! Impact parameter in Z of the track to the primary vertex
} // namespace hf_correlation_trigger_reduced

DECLARE_SOA_TABLE(HfcRedTrigs, "AOD", "HFCREDTRIG", //! Table with charm hadron candidate info
                  soa::Index<>,
                  aod::hf_correlation_trigger_reduced::HfcRedCorrCollId,
                  aod::hf_correlation_trigger_reduced::PhiTrig,
                  aod::hf_correlation_trigger_reduced::EtaTrig,
                  aod::hf_correlation_trigger_reduced::PtTrig);

DECLARE_SOA_TABLE(HfcRedTrigCharms, "AOD", "HFCREDTRIGCHARM", //! Table with Same Event Charm-Hadron pairs information
                  aod::hf_correlation_trigger_reduced::InvMassTrig,
                  aod::hf_correlation_trigger_reduced::BdtScore0Trig,
                  aod::hf_correlation_trigger_reduced::BdtScore1Trig);

DECLARE_SOA_TABLE(HfcRedTrigHads, "AOD", "HFCREDTRIGHAD", //! Table with Same Event Charm-Hadron pairs information
                  aod::hf_correlation_trigger_reduced::NTpcCrossedRowsTrig,
                  aod::hf_correlation_trigger_reduced::ItsClusterMapTrig,
                  aod::hf_correlation_trigger_reduced::ItsNClsTrig,
                  aod::hf_correlation_trigger_reduced::DcaXYTrig,
                  aod::hf_correlation_trigger_reduced::DcaZTrig);

namespace hf_assoc_track_reduced
{
DECLARE_SOA_COLUMN(OriginTrackId, originTrackId, int);      //! Original track index
DECLARE_SOA_COLUMN(NTpcCrossedRows, nTpcCrossedRows, int);  //! Number of crossed TPC Rows
DECLARE_SOA_COLUMN(ItsClusterMap, itsClusterMap, int);      //! ITS cluster map, one bit per a layer, starting from the innermost
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, int);                  //! Number of ITS clusters
DECLARE_SOA_COLUMN(EtaAssocTrack, etaAssocTrack, float);    //! Eta of the track
DECLARE_SOA_COLUMN(PhiAssocTrack, phiAssocTrack, float);    //! Phi of the track
DECLARE_SOA_COLUMN(PtAssocTrack, ptAssocTrack, float);      //! Pt of the track
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                    //! Impact parameter in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                      //! Impact parameter in Z of the track to the primary vertex
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

DECLARE_SOA_TABLE(HfcRedTrkAssocs, "AOD", "HFCREDTRKASSOC", //! Table with associated track info
                  soa::Index<>,
                  aod::hf_correlation_trigger_reduced::HfcRedCorrCollId,
                  aod::hf_assoc_track_reduced::PhiAssocTrack,
                  aod::hf_assoc_track_reduced::EtaAssocTrack,
                  aod::hf_assoc_track_reduced::PtAssocTrack);

DECLARE_SOA_TABLE(HfcRedTrkSels, "AOD", "HFCREDSETRKSEL", //! Table with Same Event Track Selections information
                  aod::hf_assoc_track_reduced::NTpcCrossedRows,
                  aod::hf_assoc_track_reduced::ItsClusterMap,
                  aod::hf_assoc_track_reduced::ItsNCls,
                  aod::hf_assoc_track_reduced::DcaXY,
                  aod::hf_assoc_track_reduced::DcaZ);

// definition of columns and tables for Charm-Hadron and Hadron-Hadron correlation pairs
namespace hf_correlation_charm_hadron_reduced
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);                                    //! DeltaPhi between charm hadron and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);                                    //! DeltaEta between charm hadron and Hadrons
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);                                        //! Pool Bin for the MixedEvent
} // namespace hf_correlation_charm_hadron_reduced

DECLARE_SOA_TABLE(HfcRedSEPairs, "AOD", "HFCREDSEPAIR", //! Table with Same Event Trig-Assoc pairs
                  aod::hf_correlation_trigger_reduced::HfcRedCorrCollId,
                  aod::hf_correlation_trigger_reduced::PtTrig,
                  aod::hf_assoc_track_reduced::PtAssocTrack,
                  aod::hf_correlation_charm_hadron_reduced::DeltaEta,
                  aod::hf_correlation_charm_hadron_reduced::DeltaPhi);

DECLARE_SOA_TABLE(HfcRedCorrPair, "AOD", "HFCREDCORRPAIR", //! Correlation pairs information
                  aod::hf_correlation_charm_hadron_reduced::PoolBin,
                  aod::hf_correlation_trigger_reduced::PtTrig,
                  aod::hf_assoc_track_reduced::PtAssocTrack,
                  aod::hf_correlation_charm_hadron_reduced::DeltaEta,
                  aod::hf_correlation_charm_hadron_reduced::DeltaPhi);

} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_DERIVEDDATACORRELATIONTABLES_H_
