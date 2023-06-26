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

#ifndef PWGHF_D2H_DATAMODEL_REDUCEDDATAMODEL_H_
#define PWGHF_D2H_DATAMODEL_REDUCEDDATAMODEL_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

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

DECLARE_SOA_TABLE(HfReducedCollisions, "AOD", "HFREDCOLLISION", //! Table with collision for reduced workflow
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

using HfReducedCollision = HfReducedCollisions::iterator;

DECLARE_SOA_TABLE(HfOriginalCollisionsCounter, "AOD", "HFCOLCOUNTER", //! Table with original number of collisions
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

DECLARE_SOA_COLUMN(Px, px, float); //! Momentum in x-direction in GeV/c
DECLARE_SOA_COLUMN(Py, py, float); //! Momentum in y-direction in GeV/c
DECLARE_SOA_COLUMN(Pz, pz, float); //! Momentum in z-direction in GeV/c
} // namespace hf_track_par_cov

// general columns
#define HFTRACKPARCOV_COLUMNS    \
  aod::track::X,                 \
    aod::track::Alpha,           \
    aod::track::Y,               \
    aod::track::Z,               \
    aod::track::Snp,             \
    aod::track::Tgl,             \
    aod::track::Signed1Pt,       \
    hf_track_par_cov::CYY,       \
    hf_track_par_cov::CZY,       \
    hf_track_par_cov::CZZ,       \
    hf_track_par_cov::CSnpY,     \
    hf_track_par_cov::CSnpZ,     \
    hf_track_par_cov::CSnpSnp,   \
    hf_track_par_cov::CTglY,     \
    hf_track_par_cov::CTglZ,     \
    hf_track_par_cov::CTglSnp,   \
    hf_track_par_cov::CTglTgl,   \
    hf_track_par_cov::C1PtY,     \
    hf_track_par_cov::C1PtZ,     \
    hf_track_par_cov::C1PtSnp,   \
    hf_track_par_cov::C1PtTgl,   \
    hf_track_par_cov::C1Pt21Pt2, \
    hf_track_par_cov::Px,        \
    hf_track_par_cov::Py,        \
    hf_track_par_cov::Pz

namespace hf_track_index_reduced
{
DECLARE_SOA_INDEX_COLUMN(HfReducedCollision, hfReducedCollision); //! ReducedCollision index
DECLARE_SOA_INDEX_COLUMN(Track, track);                           //! Track index
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);                         //! Flag to check if track has a TPC match
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);                         //! Flag to check if track has a TOF match
} // namespace hf_track_index_reduced

DECLARE_SOA_TABLE(HfTracksReduced, "AOD", "HFTRACKRED", //! Table with track information for reduced workflow
                  soa::Index<>,
                  hf_track_index_reduced::TrackId,
                  hf_track_index_reduced::HfReducedCollisionId,
                  HFTRACKPARCOV_COLUMNS);

namespace hf_track_pid_reduced
{
DECLARE_SOA_COLUMN(Pt, pt, float); //! Transverse momentum of the track in GeV/c
} // namespace hf_track_pid_reduced

// table with all attributes needed to call getStatusTrackPIDTpcAndTof() in the selector task
DECLARE_SOA_TABLE(HfTracksPidReduced, "AOD", "HFTRACKPIDRED", //! Table with PID track information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_index_reduced::HfReducedCollisionId,
                  hf_track_pid_reduced::Pt,
                  hf_track_index_reduced::HasTPC,
                  hf_track_index_reduced::HasTOF,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCNSigmaMu,
                  pidtpc::TPCNSigmaPi,
                  pidtpc::TPCNSigmaKa,
                  pidtpc::TPCNSigmaPr,
                  pidtof::TOFNSigmaEl,
                  pidtof::TOFNSigmaMu,
                  pidtof::TOFNSigmaPi,
                  pidtof::TOFNSigmaKa,
                  pidtof::TOFNSigmaPr);

namespace hf_cand_3prong_reduced
{
DECLARE_SOA_COLUMN(CPA, cpa, float);                 //! Cosinus pointing angle
DECLARE_SOA_COLUMN(DecayLength, decayLength, float); //! Decay length in cm
DECLARE_SOA_COLUMN(InvMass, invMass, float);         //! Invariant mass of 3prong candidate in GeV/c2

template <typename T>
auto invMassDplusToPiKPi(const T& pVec0, const T& pVec1, const T& pVec2)
{
  return RecoDecay::m(std::array{pVec0, pVec1, pVec2},
                      std::array{RecoDecay::getMassPDG(kPiPlus),
                                 RecoDecay::getMassPDG(kKPlus),
                                 RecoDecay::getMassPDG(kPiPlus)});
}
} // namespace hf_cand_3prong_reduced

DECLARE_SOA_TABLE(HfCand3ProngReduced, "AOD", "HFCAND3PRONGRED", //! Table with 3prong candidate information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id,
                  hf_track_index_reduced::HfReducedCollisionId,
                  HFTRACKPARCOV_COLUMNS,
                  hf_cand_3prong_reduced::CPA,
                  hf_cand_3prong_reduced::DecayLength,
                  hf_cand_3prong_reduced::InvMass);

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
DECLARE_SOA_TABLE(HfDPiMcRecReduced, "AOD", "HFDPIMCRECRED", //! Table with reconstructed MC information on DPi(<-B0) pairs for reduced workflow
                  hf_cand_b0::Prong0Id,
                  hf_track_index::Prong1Id,
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::OriginMcRec,
                  hf_cand_b0::DebugMcRec,
                  hf_b0_mc::PtMother);

// Table with same size as HFCANDB0
DECLARE_SOA_TABLE(HfB0McRecReduced, "AOD", "HFB0MCRECRED", //! Reconstruction-level MC information on B0 candidates for reduced workflow
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::OriginMcRec,
                  hf_cand_b0::DebugMcRec,
                  hf_b0_mc::PtMother);

DECLARE_SOA_TABLE(HfB0McGenReduced, "AOD", "HFB0MCGENRED", //! Generation-level MC information on B0 candidates for reduced workflow
                  hf_cand_b0::FlagMcMatchGen,
                  hf_cand_b0::OriginMcGen,
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
} // namespace hf_cand_b0_config

DECLARE_SOA_TABLE(HfCandB0Config, "AOD", "HFCANDB0CONFIG", //! Table with configurables information for reduced workflow
                  hf_cand_b0_config::MySelectionFlagD);

} // namespace aod

namespace soa
{
DECLARE_EQUIVALENT_FOR_INDEX(aod::HfCand3ProngBase, aod::HfCand3ProngReduced);
DECLARE_EQUIVALENT_FOR_INDEX(aod::StoredTracks, aod::HfTracksReduced);
DECLARE_EQUIVALENT_FOR_INDEX(aod::StoredTracks, aod::HfTracksPidReduced);
} // namespace soa
} // namespace o2
#endif // PWGHF_D2H_DATAMODEL_REDUCEDDATAMODEL_H_
