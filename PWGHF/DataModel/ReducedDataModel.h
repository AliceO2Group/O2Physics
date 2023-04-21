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

#ifndef PWGHF_DATAMODEL_REDUCED_DATA_MODEL_H_
#define PWGHF_DATAMODEL_REDUCED_DATA_MODEL_H_

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponse.h"
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
DECLARE_SOA_COLUMN(OriginalAODId, originalAODId, int);     //! User-defined index of COLLISION table processed
DECLARE_SOA_COLUMN(OriginalAODSize, originalAODSize, int); //! Size of COLLISION table processed
} // namespace hf_reduced_collision

DECLARE_SOA_TABLE(HfReducedCollisions, "AOD", "HFREDCOLLISION", //! Table with collision for reduced workflow
                  track::CollisionId,                           // FIXME : soa::Index<> does not work event with DECLARE_EQUIVALENT_FOR_INDEX(aod::Collisions, aod::HfReducedCollisions)
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::CovXX,
                  collision::CovXY,
                  collision::CovYY,
                  collision::CovXZ,
                  collision::CovYZ,
                  collision::CovZZ,
                  hf_reduced_collision::Bz,
                  hf_reduced_collision::OriginalAODId,
                  hf_reduced_collision::OriginalAODSize);

namespace hf_track_par_cov
{
// CAREFUL: the getters names shall be the same as the ones of the getTrackParCov method in Common/Core/trackUtilities.h
DECLARE_SOA_COLUMN(X, x, float);            //! X of track evaluation
DECLARE_SOA_COLUMN(Alpha, alpha, float);    //! Track frame angle
DECLARE_SOA_COLUMN(Y, y, float);            //! Y of track evaluation
DECLARE_SOA_COLUMN(Z, z, float);            //! Z of track evaluation
DECLARE_SOA_COLUMN(Snp, snp, float);        //! sin(phi)
DECLARE_SOA_COLUMN(Tgl, tgl, float);        //! tan(lambda)
DECLARE_SOA_COLUMN(Q2Pt, signed1Pt, float); //! (sign of charge)/Pt in c/GeV

DECLARE_SOA_COLUMN(SigY2, cYY, float);          //! Covariance matrix
DECLARE_SOA_COLUMN(SigZY, cZY, float);          //! Covariance matrix
DECLARE_SOA_COLUMN(SigZ2, cZZ, float);          //! Covariance matrix
DECLARE_SOA_COLUMN(SigSnpY, cSnpY, float);      //! Covariance matrix
DECLARE_SOA_COLUMN(SigSnpZ, cSnpZ, float);      //! Covariance matrix
DECLARE_SOA_COLUMN(SigSnp2, cSnpSnp, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(SigTglY, cTglY, float);      //! Covariance matrix
DECLARE_SOA_COLUMN(SigTglZ, cTglZ, float);      //! Covariance matrix
DECLARE_SOA_COLUMN(SigTglSnp, cTglSnp, float);  //! Covariance matrix
DECLARE_SOA_COLUMN(SigTgl2, cTglTgl, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(SigQ2PtY, c1PtY, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(SigQ2PtZ, c1PtZ, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(SigQ2PtSnp, c1PtSnp, float); //! Covariance matrix
DECLARE_SOA_COLUMN(SigQ2PtTgl, c1PtTgl, float); //! Covariance matrix
DECLARE_SOA_COLUMN(SigQ2Pt2, c1Pt21Pt2, float); //! Covariance matrix

DECLARE_SOA_COLUMN(Px, px, float); //! Momentum in x-direction in GeV/c
DECLARE_SOA_COLUMN(Py, py, float); //! Momentum in y-direction in GeV/c
DECLARE_SOA_COLUMN(Pz, pz, float); //! Momentum in z-direction in GeV/c
} // namespace hf_track_par_cov

// general columns
#define HFTRACKPARCOV_COLUMNS     \
  hf_track_par_cov::X,            \
    hf_track_par_cov::Alpha,      \
    hf_track_par_cov::Y,          \
    hf_track_par_cov::Z,          \
    hf_track_par_cov::Snp,        \
    hf_track_par_cov::Tgl,        \
    hf_track_par_cov::Q2Pt,       \
    hf_track_par_cov::SigY2,      \
    hf_track_par_cov::SigZY,      \
    hf_track_par_cov::SigZ2,      \
    hf_track_par_cov::SigSnpY,    \
    hf_track_par_cov::SigSnpZ,    \
    hf_track_par_cov::SigSnp2,    \
    hf_track_par_cov::SigTglY,    \
    hf_track_par_cov::SigTglZ,    \
    hf_track_par_cov::SigTglSnp,  \
    hf_track_par_cov::SigTgl2,    \
    hf_track_par_cov::SigQ2PtY,   \
    hf_track_par_cov::SigQ2PtZ,   \
    hf_track_par_cov::SigQ2PtSnp, \
    hf_track_par_cov::SigQ2PtTgl, \
    hf_track_par_cov::SigQ2Pt2,   \
    hf_track_par_cov::Px,         \
    hf_track_par_cov::Py,         \
    hf_track_par_cov::Pz

namespace hf_reduced_track_index
{
DECLARE_SOA_INDEX_COLUMN(Track, track);   //! Track index
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool); //! Flag to check if track has a TPC match
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool); //! Flag to check if track has a TOF match
} // namespace hf_reduced_track_index

DECLARE_SOA_TABLE(HfReducedTracksWithSel, "AOD", "HFREDTRACKWSEL", //! Table with track information for reduced workflow
                  soa::Index<>,
                  hf_reduced_track_index::TrackId,
                  track::CollisionId,
                  HFTRACKPARCOV_COLUMNS);

namespace hf_track_pid_with_sel
{
DECLARE_SOA_COLUMN(Pt, pt, float); //! Transverse momentum of the track in GeV/c
} // namespace hf_track_pid_with_sel

// table with all attributes needed to call getStatusTrackPIDTpcAndTof() in the selector task
DECLARE_SOA_TABLE(HfReducedTracksPIDWithSel, "AOD", "HFREDTRACKPID", //! Table with PID track information for reduced workflow
                  o2::soa::Index<>,
                  track::CollisionId,
                  hf_track_pid_with_sel::Pt,
                  hf_reduced_track_index::HasTPC,
                  hf_reduced_track_index::HasTOF,
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

namespace hf_reduced_cand_3prong
{
DECLARE_SOA_COLUMN(CPA, cpa, float);                 //! Cosinus pointing angle
DECLARE_SOA_COLUMN(DecayLength, decayLength, float); //! Decay length in cm
DECLARE_SOA_COLUMN(InvMass, invMass, float);         //! Invariant mass of 3prong candidate in GeV/c2

template <typename T>
auto invMassDplusToPiKPi(const T& pVec0, const T& pVec1, const T& pVec2)
{
  return RecoDecay::m(array{array{pVec0}, array{pVec1}, array{pVec2}},
                      array{RecoDecay::getMassPDG(kPiPlus),
                            RecoDecay::getMassPDG(kKPlus),
                            RecoDecay::getMassPDG(kPiPlus)});
}
} // namespace hf_reduced_cand_3prong

DECLARE_SOA_TABLE(HfReducedCand3Prong, "AOD", "HFREDCAND3PRONG", //! Table with 3prong candidate information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id,
                  track::CollisionId,
                  HFTRACKPARCOV_COLUMNS,
                  hf_reduced_cand_3prong::CPA,
                  hf_reduced_cand_3prong::DecayLength,
                  hf_reduced_cand_3prong::InvMass);

namespace hf_b0_mc
{
// MC Rec
DECLARE_SOA_COLUMN(PtMother, ptMother, float); //!
// MC Gen
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(YProng0, yProng0, float);
DECLARE_SOA_COLUMN(EtaProng0, etaProng0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(YProng1, yProng1, float);
DECLARE_SOA_COLUMN(EtaProng1, etaProng1, float);
} // namespace hf_b0_mc

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfReducedDPiMcRec, "AOD", "HFREDDPIMCREC", //! Table with reconstructed MC information on DPi(<-B0) pairs for reduced workflow
                  hf_cand_b0::Prong0Id,
                  hf_track_index::Prong1Id,
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::OriginMcRec,
                  hf_cand_b0::DebugMcRec,
                  hf_b0_mc::PtMother);

// Table with same size as HFCANDB0
DECLARE_SOA_TABLE(HfReducedB0McRec, "AOD", "HFREDB0MCREC", //! Reconstruction-level MC information on B0 candidates for reduced workflow
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::OriginMcRec,
                  hf_cand_b0::DebugMcRec,
                  hf_b0_mc::PtMother);

DECLARE_SOA_TABLE(HfReducedB0McGen, "AOD", "HFREDB0MCGEN", //! Generation-level MC information on B0 candidates for reduced workflow
                  hf_cand_b0::FlagMcMatchGen,
                  hf_cand_b0::OriginMcGen,
                  hf_b0_mc::Pt,
                  hf_b0_mc::Y,
                  hf_b0_mc::Eta,
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
DECLARE_SOA_COLUMN(MySelectionFlagD, mySelectionFlagD, int); //! Flag to filter selected D+ mesons
} // namespace hf_cand_b0_config

DECLARE_SOA_TABLE(HfCandB0Config, "AOD", "HFCANDB0CONFIG", //! Table with configurables information for reduced workflow
                  hf_cand_b0_config::MySelectionFlagD);

} // namespace aod

namespace soa
{
DECLARE_EQUIVALENT_FOR_INDEX(aod::HfCand3ProngBase, aod::HfReducedCand3Prong);
DECLARE_EQUIVALENT_FOR_INDEX(aod::StoredTracks, aod::HfReducedTracksWithSel);
DECLARE_EQUIVALENT_FOR_INDEX(aod::StoredTracks, aod::HfReducedTracksPIDWithSel);
} // namespace soa
} // namespace o2

#endif // PWGHF_DATAMODEL_REDUCED_DATA_MODEL_H_