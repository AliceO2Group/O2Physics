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

/// \file ReducedDataModel.h
/// \brief Header file with definition of methods and tables
//  used to fold (unfold) track and primary vertex information by writing (reading) AO2Ds
/// \note
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Bari

#ifndef PWGHF_D2H_DATAMODEL_REDUCEDDATAMODEL_H_
#define PWGHF_D2H_DATAMODEL_REDUCEDDATAMODEL_H_

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

namespace o2
{
namespace aod
{
namespace hf_reduced_collision
{
DECLARE_SOA_COLUMN(Bz, bz, float); //! Magnetic field in z-direction
// keep track of the number of studied events (for normalization purposes)
DECLARE_SOA_COLUMN(OriginalCollisionCount, originalCollisionCount, int); //! Size of COLLISION table processed
} // namespace hf_reduced_collision

DECLARE_SOA_TABLE(HfRedCollisions, "AOD", "HFREDCOLLISION", //! Table with collision for reduced workflow
                  soa::Index<>,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::CovXX,
                  collision::CovXY,
                  collision::CovYY,
                  collision::CovXZ,
                  collision::CovYZ,
                  collision::CovZZ,
                  hf_reduced_collision::Bz);

using HfRedCollision = HfRedCollisions::iterator;

DECLARE_SOA_TABLE(HfOrigColCounts, "AOD", "HFORIGCOLCOUNT", //! Table with original number of collisions
                  hf_reduced_collision::OriginalCollisionCount);

namespace hf_track_par_cov
{
// CAREFUL: the getters names shall be the same as the ones of the getTrackParCov method in Common/Core/trackUtilities.h
DECLARE_SOA_COLUMN(CYY, cYY, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CZY, cZY, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CZZ, cZZ, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpY, cSnpY, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpZ, cSnpZ, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpSnp, cSnpSnp, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(CTglY, cTglY, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CTglZ, cTglZ, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CTglSnp, cTglSnp, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(CTglTgl, cTglTgl, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtY, c1PtY, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtZ, c1PtZ, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtSnp, c1PtSnp, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtTgl, c1PtTgl, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1Pt21Pt2, c1Pt21Pt2, float); //! Covariance matrix
} // namespace hf_track_par_cov

// general columns
#define HFTRACKPAR_COLUMNS \
  aod::track::X,           \
    aod::track::Alpha,     \
    aod::track::Y,         \
    aod::track::Z,         \
    aod::track::Snp,       \
    aod::track::Tgl,       \
    aod::track::Signed1Pt

#define HFTRACKPARCOV_COLUMNS  \
  hf_track_par_cov::CYY,       \
    hf_track_par_cov::CZY,     \
    hf_track_par_cov::CZZ,     \
    hf_track_par_cov::CSnpY,   \
    hf_track_par_cov::CSnpZ,   \
    hf_track_par_cov::CSnpSnp, \
    hf_track_par_cov::CTglY,   \
    hf_track_par_cov::CTglZ,   \
    hf_track_par_cov::CTglSnp, \
    hf_track_par_cov::CTglTgl, \
    hf_track_par_cov::C1PtY,   \
    hf_track_par_cov::C1PtZ,   \
    hf_track_par_cov::C1PtSnp, \
    hf_track_par_cov::C1PtTgl, \
    hf_track_par_cov::C1Pt21Pt2

namespace hf_track_index_reduced
{
DECLARE_SOA_INDEX_COLUMN(HfRedCollision, hfRedCollision); //! ReducedCollision index
DECLARE_SOA_COLUMN(TrackId, trackId, int);                //! Original track index
DECLARE_SOA_COLUMN(Prong0Id, prong0Id, int);              //! Original track index
DECLARE_SOA_COLUMN(Prong1Id, prong1Id, int);              //! Original track index
DECLARE_SOA_COLUMN(Prong2Id, prong2Id, int);              //! Original track index
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);                 //! Flag to check if track has a TPC match
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);                 //! Flag to check if track has a TOF match
} // namespace hf_track_index_reduced

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRedTrackBases, "AOD", "HFREDTRACKBASE", //! Table with track information for reduced workflow
                  soa::Index<>,
                  hf_track_index_reduced::TrackId,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  aod::track::Px<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Py<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>);

DECLARE_SOA_TABLE(HfRedTracksCov, "AOD", "HFREDTRACKCOV", //! Table with track covariance information for reduced workflow
                  soa::Index<>,
                  HFTRACKPARCOV_COLUMNS);

// table with all attributes needed to call statusTpcAndTof() in the selector task
DECLARE_SOA_TABLE(HfRedTracksPid, "AOD", "HFREDTRACKPID", //! Table with PID track information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_index_reduced::HasTPC,
                  hf_track_index_reduced::HasTOF,
                  pidtpc::TPCNSigmaPi,
                  pidtof::TOFNSigmaPi);

DECLARE_SOA_EXTENDED_TABLE_USER(HfRedTracksExt, HfRedTrackBases, "HFREDTRACKEXT", //! Track parameters at collision vertex
                                aod::track::Pt);

using HfRedTracks = HfRedTracksExt;

namespace hf_charm_cand_reduced
{
DECLARE_SOA_COLUMN(InvMass, invMass, float);           //! Invariant mass of 2prong candidate in GeV/c2
DECLARE_SOA_COLUMN(InvMassD0, invMassD0, float);       //! Invariant mass of 2prong candidate in GeV/c2
DECLARE_SOA_COLUMN(InvMassD0Bar, invMassD0Bar, float); //! Invariant mass of 2prong candidate in GeV/c2
} // namespace hf_charm_cand_reduced

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRed2Prongs, "AOD", "HFRED2PRONG", //! Table with 2prong candidate information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_charm_cand_reduced::InvMassD0, hf_charm_cand_reduced::InvMassD0Bar,
                  aod::track::Px<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Py<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>);

DECLARE_SOA_TABLE(HfRed2ProngsCov, "AOD", "HFRED2PRONGSCOV", //! Table with 2prong candidate covariance for reduced workflow
                  o2::soa::Index<>,
                  HFTRACKPARCOV_COLUMNS,
                  o2::soa::Marker<1>);

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRed3Prongs, "AOD", "HFRED3PRONG", //! Table with 3prong candidate information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id, hf_track_index_reduced::Prong2Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_charm_cand_reduced::InvMass,
                  aod::track::Px<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Py<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>);

DECLARE_SOA_TABLE(HfRed3ProngsCov, "AOD", "HFRED3PRONGSCOV", //! Table with 3prong candidate covariance for reduced workflow
                  o2::soa::Index<>,
                  HFTRACKPARCOV_COLUMNS,
                  o2::soa::Marker<2>);

// Beauty candidates prongs
namespace hf_cand_b0_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed3Prongs, "_0");    //! Prong0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrackBases, "_1"); //! Prong1 index
} // namespace hf_cand_b0_reduced

DECLARE_SOA_TABLE(HfRedB0Prongs, "AOD", "HFREDB0PRONG",
                  hf_cand_b0_reduced::Prong0Id, hf_cand_b0_reduced::Prong1Id);

using HfRedCandB0 = soa::Join<HfCandB0Ext, HfRedB0Prongs>;

namespace hf_cand_bplus_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed2Prongs, "_0");    //! Prong0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrackBases, "_1"); //! Prong1 index
} // namespace hf_cand_bplus_reduced

DECLARE_SOA_TABLE(HfRedBplusProngs, "AOD", "HFREDBPPRONG",
                  hf_cand_bplus_reduced::Prong0Id, hf_cand_bplus_reduced::Prong1Id);

using HfRedCandBplus = soa::Join<HfCandBplusExt, HfRedBplusProngs>;

namespace hf_b0_mc
{
// MC Rec
DECLARE_SOA_COLUMN(PtMother, ptMother, float); //! Transverse momentum of the mother in GeV/c
// MC Gen
DECLARE_SOA_COLUMN(PtTrack, ptTrack, float);     //! Transverse momentum of the track in GeV/c
DECLARE_SOA_COLUMN(YTrack, yTrack, float);       //! Rapidity of the track
DECLARE_SOA_COLUMN(EtaTrack, etaTrack, float);   //! Pseudorapidity of the track
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);   //! Transverse momentum of the track's prong0 in GeV/c
DECLARE_SOA_COLUMN(YProng0, yProng0, float);     //! Rapidity of the track's prong0
DECLARE_SOA_COLUMN(EtaProng0, etaProng0, float); //! Pseudorapidity of the track's prong0
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);   //! Transverse momentum of the track's prong1 in GeV/c
DECLARE_SOA_COLUMN(YProng1, yProng1, float);     //! Rapidity of the track's prong1
DECLARE_SOA_COLUMN(EtaProng1, etaProng1, float); //! Pseudorapidity of the track's prong1
} // namespace hf_b0_mc

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedDpPis, "AOD", "HFMCRECREDDPPI", //! Table with reconstructed MC information on DPi(<-B0) pairs for reduced workflow
                  hf_cand_b0_reduced::Prong0Id,
                  hf_cand_b0_reduced::Prong1Id,
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::DebugMcRec,
                  hf_b0_mc::PtMother);

// Table with same size as HFCANDB0
DECLARE_SOA_TABLE(HfMcRecRedB0s, "AOD", "HFMCRECREDB0", //! Reconstruction-level MC information on B0 candidates for reduced workflow
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::DebugMcRec,
                  hf_b0_mc::PtMother);

DECLARE_SOA_TABLE(HfMcGenRedB0s, "AOD", "HFMCGENREDB0", //! Generation-level MC information on B0 candidates for reduced workflow
                  hf_cand_b0::FlagMcMatchGen,
                  hf_b0_mc::PtTrack,
                  hf_b0_mc::YTrack,
                  hf_b0_mc::EtaTrack,
                  hf_b0_mc::PtProng0,
                  hf_b0_mc::YProng0,
                  hf_b0_mc::EtaProng0,
                  hf_b0_mc::PtProng1,
                  hf_b0_mc::YProng1,
                  hf_b0_mc::EtaProng1);

// store all configurables values used in the first part of the workflow
// so we can use them in the B0 part
namespace hf_cand_b0_config
{
DECLARE_SOA_COLUMN(MySelectionFlagD, mySelectionFlagD, int8_t); //! Flag to filter selected D+ mesons
DECLARE_SOA_COLUMN(MyInvMassWindowDPi, myInvMassWindowDPi, float); //! Half-width of the B0 invariant-mass window in GeV/c2
} // namespace hf_cand_b0_config

DECLARE_SOA_TABLE(HfCandB0Configs, "AOD", "HFCANDB0CONFIG", //! Table with configurables information for reduced workflow
                  hf_cand_b0_config::MySelectionFlagD,
                  hf_cand_b0_config::MyInvMassWindowDPi);

namespace hf_bplus_mc
{
// MC Rec
DECLARE_SOA_COLUMN(PtMother, ptMother, float); //! Transverse momentum of the mother in GeV/c
// MC Gen
DECLARE_SOA_COLUMN(PtTrack, ptTrack, float);     //! Transverse momentum of the track in GeV/c
DECLARE_SOA_COLUMN(YTrack, yTrack, float);       //! Rapidity of the track
DECLARE_SOA_COLUMN(EtaTrack, etaTrack, float);   //! Pseudorapidity of the track
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);   //! Transverse momentum of the track's prong0 in GeV/c
DECLARE_SOA_COLUMN(YProng0, yProng0, float);     //! Rapidity of the track's prong0
DECLARE_SOA_COLUMN(EtaProng0, etaProng0, float); //! Pseudorapidity of the track's prong0
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);   //! Transverse momentum of the track's prong1 in GeV/c
DECLARE_SOA_COLUMN(YProng1, yProng1, float);     //! Rapidity of the track's prong1
DECLARE_SOA_COLUMN(EtaProng1, etaProng1, float); //! Pseudorapidity of the track's prong1
} // namespace hf_bplus_mc

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedD0Pis, "AOD", "HFMCRECREDD0PI", //! Table with reconstructed MC information on D0Pi(<-B+) pairs for reduced workflow
                  hf_cand_bplus_reduced::Prong0Id,
                  hf_cand_bplus_reduced::Prong1Id,
                  hf_cand_bplus::FlagMcMatchRec,
                  hf_bplus_mc::PtMother);

// Table with same size as HFCANDBPLUS
DECLARE_SOA_TABLE(HfMcRecRedBps, "AOD", "HFMCRECREDBP", //! Reconstruction-level MC information on B+ candidates for reduced workflow
                  hf_cand_bplus::FlagMcMatchRec,
                  hf_bplus_mc::PtMother);

DECLARE_SOA_TABLE(HfMcGenRedBps, "AOD", "HFMCGENREDBP", //! Generation-level MC information on B+ candidates for reduced workflow
                  hf_cand_bplus::FlagMcMatchGen,
                  hf_bplus_mc::PtTrack,
                  hf_bplus_mc::YTrack,
                  hf_bplus_mc::EtaTrack,
                  hf_bplus_mc::PtProng0,
                  hf_bplus_mc::YProng0,
                  hf_bplus_mc::EtaProng0,
                  hf_bplus_mc::PtProng1,
                  hf_bplus_mc::YProng1,
                  hf_bplus_mc::EtaProng1);

// store all configurables values used in the first part of the workflow
// so we can use them in the Bplus part
namespace hf_cand_bplus_config
{
DECLARE_SOA_COLUMN(MySelectionFlagD0, mySelectionFlagD0, int8_t);       //! Flag to filter selected D0 mesons
DECLARE_SOA_COLUMN(MySelectionFlagD0bar, mySelectionFlagD0bar, int8_t); //! Flag to filter selected D0 mesons
DECLARE_SOA_COLUMN(MyInvMassWindowD0Pi, myInvMassWindowD0Pi, float);    //! Half-width of the Bplus invariant-mass window in GeV/c2
} // namespace hf_cand_bplus_config

DECLARE_SOA_TABLE(HfCandBpConfigs, "AOD", "HFCANDBPCONFIG", //! Table with configurables information for reduced workflow
                  hf_cand_bplus_config::MySelectionFlagD0,
                  hf_cand_bplus_config::MySelectionFlagD0bar,
                  hf_cand_bplus_config::MyInvMassWindowD0Pi);
} // namespace aod

namespace soa
{
DECLARE_EQUIVALENT_FOR_INDEX(aod::HfCand2ProngBase, aod::HfRed2Prongs);
DECLARE_EQUIVALENT_FOR_INDEX(aod::HfCand3ProngBase, aod::HfRed3Prongs);
DECLARE_EQUIVALENT_FOR_INDEX(aod::StoredTracks, aod::HfRedTrackBases);
DECLARE_EQUIVALENT_FOR_INDEX(aod::Collisions, aod::HfRedCollisions);
} // namespace soa
} // namespace o2
#endif // PWGHF_D2H_DATAMODEL_REDUCEDDATAMODEL_H_
