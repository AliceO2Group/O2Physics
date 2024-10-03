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
DECLARE_SOA_COLUMN(HfCollisionRejectionMap, hfCollisionRejectionMap, uint16_t); //! Bitmask with failed selection criteria
// keep track of the number of studied events (for normalization purposes)
DECLARE_SOA_COLUMN(OriginalCollisionCount, originalCollisionCount, int); //! Size of COLLISION table processed
DECLARE_SOA_COLUMN(ZvtxSelectedCollisionCount, zvtxSelectedCollisionCount, int);                                                 //! Number of COLLISIONS with |zvtx| < zvtxMax
DECLARE_SOA_COLUMN(TriggerSelectedCollisionCount, triggerSelectedCollisionCount, int);                                           //! Number of COLLISIONS with sel8
DECLARE_SOA_COLUMN(ZvtxAndTriggerSelectedCollisionCount, zvtxAndTriggerSelectedCollisionCount, int);                             //! Number of COLLISIONS with |zvtx| < zvtxMax and sel8
DECLARE_SOA_COLUMN(ZvtxAndTriggerAndSoftTriggerSelectedCollisionCount, zvtxAndTriggerAndSoftTriggerSelectedCollisionCount, int); //! Number of COLLISIONS with |zvtx| < zvtxMax, sel8, and selected by the software trigger
DECLARE_SOA_COLUMN(AllSelectionsCollisionCount, allSelectionsCollisionCount, int);                                               //! Number of COLLISIONS that passed all selections
} // namespace hf_reduced_collision

DECLARE_SOA_TABLE(HfRedCollisions, "AOD", "HFREDCOLLISION", //! Table with collision for reduced workflow
                  o2::soa::Index<>,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::NumContrib,
                  hf_reduced_collision::HfCollisionRejectionMap,
                  hf_reduced_collision::Bz,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfRedCollExtras, "AOD", "HFREDCOLLEXTRA", //! Table with collision extras for reduced workflow
                  collision::CovXX,
                  collision::CovXY,
                  collision::CovYY,
                  collision::CovXZ,
                  collision::CovYZ,
                  collision::CovZZ);

using HfRedCollision = HfRedCollisions::iterator;

DECLARE_SOA_TABLE(HfOrigColCounts, "AOD", "HFORIGCOLCOUNT", //! Table with original number of collisions
                  hf_reduced_collision::OriginalCollisionCount,
                  hf_reduced_collision::ZvtxSelectedCollisionCount,
                  hf_reduced_collision::TriggerSelectedCollisionCount,
                  hf_reduced_collision::ZvtxAndTriggerSelectedCollisionCount,
                  hf_reduced_collision::ZvtxAndTriggerAndSoftTriggerSelectedCollisionCount,
                  hf_reduced_collision::AllSelectionsCollisionCount);

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
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>,
                  aod::track::PVector<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha, aod::track::Tgl>);

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
DECLARE_SOA_COLUMN(InvMass, invMass, float);                                     //! Invariant mass of 2prong candidate in GeV/c2
DECLARE_SOA_COLUMN(InvMassD0, invMassD0, float);                                 //! Invariant mass of 2prong candidate in GeV/c2
DECLARE_SOA_COLUMN(InvMassD0Bar, invMassD0Bar, float);                           //! Invariant mass of 2prong candidate in GeV/c2
DECLARE_SOA_COLUMN(MlScoreBkgMassHypo0, mlScoreBkgMassHypo0, float);             //! ML score for background class (mass hypothesis 0)
DECLARE_SOA_COLUMN(MlScorePromptMassHypo0, mlScorePromptMassHypo0, float);       //! ML score for prompt class (mass hypothesis 0)
DECLARE_SOA_COLUMN(MlScoreNonpromptMassHypo0, mlScoreNonpromptMassHypo0, float); //! ML score for non-prompt class (mass hypothesis 0)
DECLARE_SOA_COLUMN(MlScoreBkgMassHypo1, mlScoreBkgMassHypo1, float);             //! ML score for background class (mass hypothesis 1)
DECLARE_SOA_COLUMN(MlScorePromptMassHypo1, mlScorePromptMassHypo1, float);       //! ML score for prompt class (mass hypothesis 1)
DECLARE_SOA_COLUMN(MlScoreNonpromptMassHypo1, mlScoreNonpromptMassHypo1, float); //! ML score for non-prompt class (mass hypothesis 1)
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
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>,
                  aod::track::PVector<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha, aod::track::Tgl>);

DECLARE_SOA_TABLE(HfRed2ProngsCov, "AOD", "HFRED2PRONGSCOV", //! Table with 2prong candidate covariance for reduced workflow
                  o2::soa::Index<>,
                  HFTRACKPARCOV_COLUMNS,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfRed2ProngsMl, "AOD", "HFRED2PRONGML", //! Table with 2prong candidate ML scores
                  hf_charm_cand_reduced::MlScoreBkgMassHypo0,
                  hf_charm_cand_reduced::MlScorePromptMassHypo0,
                  hf_charm_cand_reduced::MlScoreNonpromptMassHypo0,
                  hf_charm_cand_reduced::MlScoreBkgMassHypo1,
                  hf_charm_cand_reduced::MlScorePromptMassHypo1,
                  hf_charm_cand_reduced::MlScoreNonpromptMassHypo1);

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
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>,
                  aod::track::PVector<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha, aod::track::Tgl>);

DECLARE_SOA_TABLE(HfRed3ProngsCov, "AOD", "HFRED3PRONGSCOV", //! Table with 3prong candidate covariance for reduced workflow
                  o2::soa::Index<>,
                  HFTRACKPARCOV_COLUMNS,
                  o2::soa::Marker<2>);

DECLARE_SOA_TABLE(HfRed3ProngsMl, "AOD", "HFRED3PRONGML", //! Table with 3prong candidate ML scores
                  hf_charm_cand_reduced::MlScoreBkgMassHypo0,
                  hf_charm_cand_reduced::MlScorePromptMassHypo0,
                  hf_charm_cand_reduced::MlScoreNonpromptMassHypo0);

// Beauty candidates prongs
namespace hf_cand_b0_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed3Prongs, "_0");    //! Prong0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrackBases, "_1"); //! Prong1 index
DECLARE_SOA_COLUMN(Prong0MlScoreBkg, prong0MlScoreBkg, float);             //! Bkg ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScorePrompt, prong0MlScorePrompt, float);       //! Prompt ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScoreNonprompt, prong0MlScoreNonprompt, float); //! Nonprompt ML score of the D daughter
} // namespace hf_cand_b0_reduced

DECLARE_SOA_TABLE(HfRedB0Prongs, "AOD", "HFREDB0PRONG", //! Table with B0 daughter indices
                  hf_cand_b0_reduced::Prong0Id, hf_cand_b0_reduced::Prong1Id);

DECLARE_SOA_TABLE(HfRedB0DpMls, "AOD", "HFREDB0DPML", //! Table with ML scores for the D+ daughter
                  hf_cand_b0_reduced::Prong0MlScoreBkg,
                  hf_cand_b0_reduced::Prong0MlScorePrompt,
                  hf_cand_b0_reduced::Prong0MlScoreNonprompt,
                  o2::soa::Marker<1>);

using HfRedCandB0 = soa::Join<HfCandB0Ext, HfRedB0Prongs>;

namespace hf_cand_bplus_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed2Prongs, "_0");    //! Prong0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrackBases, "_1"); //! Prong1 index
DECLARE_SOA_COLUMN(Prong0MlScoreBkg, prong0MlScoreBkg, float);             //! Bkg ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScorePrompt, prong0MlScorePrompt, float);       //! Prompt ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScoreNonprompt, prong0MlScoreNonprompt, float); //! Nonprompt ML score of the D daughter
} // namespace hf_cand_bplus_reduced

DECLARE_SOA_TABLE(HfRedBplusProngs, "AOD", "HFREDBPPRONG",
                  hf_cand_bplus_reduced::Prong0Id, hf_cand_bplus_reduced::Prong1Id);

DECLARE_SOA_TABLE(HfRedBplusD0Mls, "AOD", "HFREDBPLUSD0ML", //! Table with ML scores for the D+ daughter
                  hf_cand_bplus_reduced::Prong0MlScoreBkg,
                  hf_cand_bplus_reduced::Prong0MlScorePrompt,
                  hf_cand_bplus_reduced::Prong0MlScoreNonprompt,
                  o2::soa::Marker<1>);

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

DECLARE_SOA_COLUMN(PdgCodeBeautyMother, pdgCodeBeautyMother, int); //! Pdg code of beauty mother
DECLARE_SOA_COLUMN(PdgCodeCharmMother, pdgCodeCharmMother, int);   //! Pdg code of charm mother
DECLARE_SOA_COLUMN(PdgCodeProng0, pdgCodeProng0, int);             //! Pdg code of prong0
DECLARE_SOA_COLUMN(PdgCodeProng1, pdgCodeProng1, int);             //! Pdg code of prong1
DECLARE_SOA_COLUMN(PdgCodeProng2, pdgCodeProng2, int);             //! Pdg code of prong2
DECLARE_SOA_COLUMN(PdgCodeProng3, pdgCodeProng3, int);             //! Pdg code of prong3
} // namespace hf_b0_mc

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedDpPis, "AOD", "HFMCRECREDDPPI", //! Table with reconstructed MC information on DPi(<-B0) pairs for reduced workflow
                  hf_cand_b0_reduced::Prong0Id,
                  hf_cand_b0_reduced::Prong1Id,
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::DebugMcRec,
                  hf_b0_mc::PtMother);

// try with extended table ?
// DECLARE_SOA_EXTENDED_TABLE_USER(ExTable, Tracks, "EXTABLE",
DECLARE_SOA_TABLE(HfMcCheckDpPis, "AOD", "HFMCCHECKDPPI", //! Table with reconstructed MC information on DPi(<-B0) pairs for MC checks in reduced workflow
                  hf_b0_mc::PdgCodeBeautyMother,
                  hf_b0_mc::PdgCodeCharmMother,
                  hf_b0_mc::PdgCodeProng0,
                  hf_b0_mc::PdgCodeProng1,
                  hf_b0_mc::PdgCodeProng2,
                  hf_b0_mc::PdgCodeProng3,
                  o2::soa::Marker<1>);

// Table with same size as HFCANDB0
DECLARE_SOA_TABLE(HfMcRecRedB0s, "AOD", "HFMCRECREDB0", //! Reconstruction-level MC information on B0 candidates for reduced workflow
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::DebugMcRec,
                  hf_b0_mc::PtMother);

DECLARE_SOA_TABLE(HfMcCheckB0s, "AOD", "HFMCCHECKB0", //! Table with reconstructed MC information on B0 candidates for MC checks in reduced workflow
                  hf_b0_mc::PdgCodeBeautyMother,
                  hf_b0_mc::PdgCodeCharmMother,
                  hf_b0_mc::PdgCodeProng0,
                  hf_b0_mc::PdgCodeProng1,
                  hf_b0_mc::PdgCodeProng2,
                  hf_b0_mc::PdgCodeProng3,
                  o2::soa::Marker<2>);

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
DECLARE_SOA_COLUMN(MySelectionFlagD, mySelectionFlagD, int8_t);    //! Flag to filter selected D+ mesons
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

DECLARE_SOA_COLUMN(PdgCodeBeautyMother, pdgCodeBeautyMother, int); //! Pdg code of beauty mother
DECLARE_SOA_COLUMN(PdgCodeProng0, pdgCodeProng0, int);             //! Pdg code of prong0
DECLARE_SOA_COLUMN(PdgCodeProng1, pdgCodeProng1, int);             //! Pdg code of prong1
DECLARE_SOA_COLUMN(PdgCodeProng2, pdgCodeProng2, int);             //! Pdg code of prong2
} // namespace hf_bplus_mc

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedD0Pis, "AOD", "HFMCRECREDD0PI", //! Table with reconstructed MC information on D0Pi(<-B+) pairs for reduced workflow
                  hf_cand_bplus_reduced::Prong0Id,
                  hf_cand_bplus_reduced::Prong1Id,
                  hf_cand_bplus::FlagMcMatchRec,
                  hf_cand_bplus::DebugMcRec,
                  hf_bplus_mc::PtMother);

// DECLARE_SOA_EXTENDED_TABLE_USER(ExTable, Tracks, "EXTABLE",
DECLARE_SOA_TABLE(HfMcCheckD0Pis, "AOD", "HFMCCHECKD0PI", //! Table with reconstructed MC information on D0Pi(<-B0) pairs for MC checks in reduced workflow
                  hf_bplus_mc::PdgCodeBeautyMother,
                  hf_bplus_mc::PdgCodeProng0,
                  hf_bplus_mc::PdgCodeProng1,
                  hf_bplus_mc::PdgCodeProng2,
                  o2::soa::Marker<1>);

// Table with same size as HFCANDBPLUS
DECLARE_SOA_TABLE(HfMcRecRedBps, "AOD", "HFMCRECREDBP", //! Reconstruction-level MC information on B+ candidates for reduced workflow
                  hf_cand_bplus::FlagMcMatchRec,
                  hf_cand_bplus::DebugMcRec,
                  hf_bplus_mc::PtMother);

DECLARE_SOA_TABLE(HfMcCheckBps, "AOD", "HFMCCHECKBP", //! Table with reconstructed MC information on B+ candidates for MC checks in reduced workflow
                  hf_bplus_mc::PdgCodeBeautyMother,
                  hf_bplus_mc::PdgCodeProng0,
                  hf_bplus_mc::PdgCodeProng1,
                  hf_bplus_mc::PdgCodeProng2,
                  o2::soa::Marker<2>);

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

// Charm resonances analysis
namespace hf_reso_cand_reduced
{
DECLARE_SOA_COLUMN(InvMass, invMass, float);             //! Invariant mass in GeV/c2
DECLARE_SOA_COLUMN(InvMassProng0, invMassProng0, float); //! Invariant Mass of D daughter in GeV/c
DECLARE_SOA_COLUMN(InvMassProng1, invMassProng1, float); //! Invariant Mass of V0 daughter in GeV/c
DECLARE_SOA_COLUMN(MlScoreBkgProng0, mlScoreBkgProng0, float);             //! Bkg ML score of the D daughter
DECLARE_SOA_COLUMN(MlScorePromptProng0, mlScorePromptProng0, float);       //! Prompt ML score of the D daughter
DECLARE_SOA_COLUMN(MlScoreNonpromptProng0, mlScoreNonpromptProng0, float); //! Nonprompt ML score of the D daughter

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float pxProng0, float pxProng1, float pyProng0, float pyProng1) -> float { return RecoDecay::pt((1.f * pxProng0 + 1.f * pxProng1), (1.f * pyProng0 + 1.f * pyProng1)); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //!
                           [](float pxProng0, float pyProng0) -> float { return RecoDecay::pt(pxProng0, pyProng0); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float pxProng1, float pyProng1) -> float { return RecoDecay::pt(pxProng1, pyProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStarDs1, cosThetaStarDs1, //! costhetastar under Ds1 hypothesis
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float invMass) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{o2::constants::physics::MassDStar, o2::constants::physics::MassK0}, invMass, 1); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStarDs2Star, cosThetaStarDs2Star, //! costhetastar under Ds2Star hypothesis
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float invMass) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassK0}, invMass, 1); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStarXiC3055, cosThetaStarXiC3055, //! costhetastar under XiC3055 hypothesis
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float invMass) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassLambda0}, invMass, 1); });
} // namespace hf_reso_cand_reduced

namespace hf_reso_3_prong
{
DECLARE_SOA_COLUMN(DType, dType, int8_t); //! Integer with selected D candidate type: 1 = Dplus, -1 = Dminus, 2 = DstarPlus, -2 = DstarMinus

DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //!
                           [](float pxProng0, float pxProng1, float pxProng2) -> float { return 1.f * pxProng0 + 1.f * pxProng1 + 1.f * pxProng2; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pyProng0, float pyProng1, float pyProng2) -> float { return 1.f * pyProng0 + 1.f * pyProng1 + 1.f * pyProng2; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pzProng0, float pzProng1, float pzProng2) -> float { return 1.f * pzProng0 + 1.f * pzProng1 + 1.f * pzProng2; });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float pxProng0, float pxProng1, float pxProng2, float pyProng0, float pyProng1, float pyProng2) -> float { return RecoDecay::pt((1.f * pxProng0 + 1.f * pxProng1 + 1.f * pxProng2), (1.f * pyProng0 + 1.f * pyProng1 + 1.f * pyProng2)); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassDplus, invMassDplus,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus, constants::physics::MassPiPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassDstar, invMassDstar,
                           [](float pxSoftPi, float pySoftPi, float pzSoftPi, float pxProng0, float pyProng0, float pzProng0, float pxProng1, float pyProng1, float pzProng1)
                             -> float { return RecoDecay::m(std::array{std::array{pxSoftPi, pySoftPi, pzSoftPi}, std::array{pxProng0, pyProng0, pzProng0}, std::array{pxProng1, pyProng1, pzProng1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassPiPlus, constants::physics::MassKPlus}) - RecoDecay::m(std::array{std::array{pxProng0, pyProng0, pzProng0}, std::array{pxProng1, pyProng1, pzProng1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassAntiDstar, invMassAntiDstar,
                           [](float pxSoftPi, float pySoftPi, float pzSoftPi, float pxProng0, float pyProng0, float pzProng0, float pxProng1, float pyProng1, float pzProng1)
                             -> float { return RecoDecay::m(std::array{std::array{pxSoftPi, pySoftPi, pzSoftPi}, std::array{pxProng0, pyProng0, pzProng0}, std::array{pxProng1, pyProng1, pzProng1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus, constants::physics::MassPiPlus}) - RecoDecay::m(std::array{std::array{pxProng0, pyProng0, pzProng0}, std::array{pxProng1, pyProng1, pzProng1}}, std::array{constants::physics::MassKPlus, constants::physics::MassPiPlus}); });
} // namespace hf_reso_3_prong

namespace hf_reso_v0
{
DECLARE_SOA_COLUMN(Cpa, cpa, float);         //! Cosine of Pointing Angle of V0 candidate
DECLARE_SOA_COLUMN(Dca, dca, float);         //! DCA of V0 candidate
DECLARE_SOA_COLUMN(Radius, radius, float);   //! Radius of V0 candidate
DECLARE_SOA_COLUMN(V0Type, v0Type, uint8_t); //! Bitmap with mass hypothesis of the V0
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,           //!
                           [](float pxProng0, float pxProng1) -> float { return 1.f * pxProng0 + 1.f * pxProng1; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pyProng0, float pyProng1) -> float { return 1.f * pyProng0 + 1.f * pyProng1; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pzProng0, float pzProng1) -> float { return 1.f * pzProng0 + 1.f * pzProng1; });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float pxProng0, float pxProng1, float pyProng0, float pyProng1) -> float { return RecoDecay::pt((1.f * pxProng0 + 1.f * pxProng1), (1.f * pyProng0 + 1.f * pyProng1)); });
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0Radius, //! V0 decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassLambda, invMassLambda, //! mass under lambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassAntiLambda, invMassAntiLambda, //! mass under antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton}); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassK0s, invMassK0s, //! mass under K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged}); });
} // namespace hf_reso_v0

namespace hf_reso_track
{
DECLARE_SOA_COLUMN(Px, px, float);                   //! x-component of momentum
DECLARE_SOA_COLUMN(Py, py, float);                   //! y-component of momentum
DECLARE_SOA_COLUMN(Pz, pz, float);                   //! z-component of momentum
DECLARE_SOA_COLUMN(Sign, sign, uint8_t);             //! charge sign
DECLARE_SOA_COLUMN(NSigmaTpcPi, nSigmaTpcPi, float); //! TPC Nsigma for pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcKa, nSigmaTpcKa, float); //! TPC Nsigma for kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcPr, nSigmaTpcPr, float); //! TPC Nsigma for proton hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPi, nSigmaTofPi, float); //! TOF Nsigma for pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTofKa, nSigmaTofKa, float); //! TOF Nsigma for kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPr, nSigmaTofPr, float); //! TOF Nsigma for proton hypothesis
DECLARE_SOA_COLUMN(HasTof, hasTof, bool);            //! flag for presence of TOF
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,                   //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //!
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array<float, 3>{px, py, pz}); });

} // namespace hf_reso_track

DECLARE_SOA_TABLE(HfRedVzeros, "AOD", "HFREDVZERO", //! Table with V0 candidate information for resonances reduced workflow
                  o2::soa::Index<>,
                  // Indices
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  // Static
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_reso_v0::Cpa,
                  hf_reso_v0::Dca,
                  hf_reso_v0::V0Type,
                  // Dynamic
                  hf_reso_v0::Px<hf_cand::PxProng0, hf_cand::PxProng1>,
                  hf_reso_v0::Py<hf_cand::PyProng0, hf_cand::PyProng1>,
                  hf_reso_v0::Pz<hf_cand::PzProng0, hf_cand::PzProng1>,
                  hf_reso_v0::InvMassK0s<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_v0::InvMassLambda<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_v0::InvMassAntiLambda<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_v0::V0Radius<hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,
                  hf_reso_v0::Pt<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PyProng0, hf_cand::PyProng1>);

DECLARE_SOA_TABLE(HfRedTrkNoParams, "AOD", "HFREDTRKNOPARAM", //! Table with tracks without track parameters for resonances reduced workflow
                  o2::soa::Index<>,
                  // Indices
                  hf_track_index_reduced::HfRedCollisionId,
                  // Static
                  hf_reso_track::Px,
                  hf_reso_track::Py,
                  hf_reso_track::Pz,
                  hf_reso_track::Sign,
                  hf_reso_track::NSigmaTpcPi,
                  hf_reso_track::NSigmaTpcKa,
                  hf_reso_track::NSigmaTpcPr,
                  hf_reso_track::NSigmaTofPi,
                  hf_reso_track::NSigmaTofKa,
                  hf_reso_track::NSigmaTofPr,
                  hf_reso_track::HasTof,
                  // Dynamic
                  hf_reso_track::Pt<hf_reso_track::Px, hf_reso_track::Py>,
                  hf_reso_track::Eta<hf_reso_track::Px, hf_reso_track::Py, hf_reso_track::Pz>,
                  hf_reso_track::Phi<hf_reso_track::Px, hf_reso_track::Py>);

DECLARE_SOA_TABLE(HfRed3PrNoTrks, "AOD", "HFRED3PRNOTRK", //! Table with 3 prong candidate information for resonances reduced workflow
                  o2::soa::Index<>,
                  // Indices
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id, hf_track_index_reduced::Prong2Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  // Static
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_reso_3_prong::DType,
                  // Dynamic
                  hf_reso_3_prong::Px<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PxProng2>,
                  hf_reso_3_prong::Py<hf_cand::PyProng0, hf_cand::PyProng1, hf_cand::PyProng2>,
                  hf_reso_3_prong::Pz<hf_cand::PzProng0, hf_cand::PzProng1, hf_cand::PzProng2>,
                  hf_reso_3_prong::InvMassDplus<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_reso_3_prong::InvMassDstar<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_3_prong::InvMassAntiDstar<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_3_prong::Pt<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PxProng2, hf_cand::PyProng0, hf_cand::PyProng1, hf_cand::PyProng2>);

DECLARE_SOA_TABLE(HfCandCharmReso, "AOD", "HFCANDCHARMRESO", //! Table with Resonance candidate information for resonances reduced workflow
                  o2::soa::Index<>,
                  // Indices
                  hf_track_index_reduced::HfRedCollisionId,
                  // Static
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_reso_cand_reduced::InvMass,
                  hf_reso_cand_reduced::InvMassProng0,
                  hf_reso_cand_reduced::InvMassProng1,
                  hf_reso_v0::Cpa,
                  hf_reso_v0::Dca,
                  hf_reso_v0::Radius,
                  // Dynamic
                  hf_reso_cand_reduced::Pt<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PyProng0, hf_cand::PyProng1>,
                  hf_reso_cand_reduced::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_reso_cand_reduced::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_cand_reduced::CosThetaStarDs1<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_reso_cand_reduced::InvMass>,
                  hf_reso_cand_reduced::CosThetaStarDs2Star<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_reso_cand_reduced::InvMass>,
                  hf_reso_cand_reduced::CosThetaStarXiC3055<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_reso_cand_reduced::InvMass>);

DECLARE_SOA_TABLE(HfCharmResoMLs, "AOD", "HFCHARMRESOML", //! Table with ML scores for the D daughter
                  hf_reso_cand_reduced::MlScoreBkgProng0,
                  hf_reso_cand_reduced::MlScorePromptProng0,
                  hf_reso_cand_reduced::MlScoreNonpromptProng0,
                  o2::soa::Marker<1>);
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
