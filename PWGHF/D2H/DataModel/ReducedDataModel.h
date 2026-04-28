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
/// \note used to fold (unfold) track and primary vertex information by writing (reading) AO2Ds
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Bari
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Luca Aglietta <luca.aglietta@cern.ch>, Università degli Studi di Torino (UniTO)
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi di Torino (UniTO)

#ifndef PWGHF_D2H_DATAMODEL_REDUCEDDATAMODEL_H_
#define PWGHF_D2H_DATAMODEL_REDUCEDDATAMODEL_H_

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsPid.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>

namespace o2
{
namespace aod
{
namespace hf_reduced_collision
{
DECLARE_SOA_COLUMN(Bz, bz, float);                                                                            //! Magnetic field in z-direction
DECLARE_SOA_COLUMN(HfCollisionRejectionMap, hfCollisionRejectionMap, o2::hf_evsel::HfCollisionRejectionMask); //! Bitmask with failed selection criteria
// keep track of the number of studied events (for normalization purposes)
DECLARE_SOA_COLUMN(OriginalCollisionCount, originalCollisionCount, int);                                                         //! Size of COLLISION table processed
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

DECLARE_SOA_TABLE(HfRedCollCents, "AOD", "HFREDCOLLCENT", //! Table with collision centrality for reduced workflow
                  cent::CentFT0C,
                  cent::CentFT0M,
                  evsel::NumTracksInTimeRange,
                  evsel::SumAmpFT0CInTimeRange);

DECLARE_SOA_TABLE(HfRedQvectors, "AOD", "HFREDQVECTOR", //! Table with collision centrality for reduced workflow
                  qvec::QvecFT0CRe,
                  qvec::QvecFT0CIm,
                  qvec::SumAmplFT0C,
                  qvec::QvecFT0ARe,
                  qvec::QvecFT0AIm,
                  qvec::SumAmplFT0A,
                  qvec::QvecFT0MRe,
                  qvec::QvecFT0MIm,
                  qvec::SumAmplFT0M,
                  qvec::QvecTPCposRe,
                  qvec::QvecTPCposIm,
                  qvec::NTrkTPCpos,
                  qvec::QvecTPCnegRe,
                  qvec::QvecTPCnegIm,
                  qvec::NTrkTPCneg,
                  qvec::QvecTPCallRe,
                  qvec::QvecTPCallIm,
                  qvec::NTrkTPCall);

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
} // namespace hf_track_index_reduced

namespace hf_track_vars_reduced
{
// CAREFUL: the getters names shall be the same as the ones of the getTrackParCov method in Common/Core/trackUtilities.h
DECLARE_SOA_COLUMN(Px, px, float);                                               //! x-component of momentum
DECLARE_SOA_COLUMN(Py, py, float);                                               //! y-component of momentum
DECLARE_SOA_COLUMN(Pz, pz, float);                                               //! z-component of momentum
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                          //! charge sign
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);                                        //! Flag to check if track has a TPC match
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);                                        //! Flag to check if track has a TOF match
DECLARE_SOA_COLUMN(HasTPCProng0, hasTPCProng0, bool);                            //! Flag to check if prong0 has a TPC match
DECLARE_SOA_COLUMN(HasTOFProng0, hasTOFProng0, bool);                            //! Flag to check if prong0 has a TOF match
DECLARE_SOA_COLUMN(HasTPCProng1, hasTPCProng1, bool);                            //! Flag to check if prong1 has a TPC match
DECLARE_SOA_COLUMN(HasTOFProng1, hasTOFProng1, bool);                            //! Flag to check if prong1 has a TOF match
DECLARE_SOA_COLUMN(HasTPCProng2, hasTPCProng2, bool);                            //! Flag to check if prong2 has a TPC match
DECLARE_SOA_COLUMN(HasTOFProng2, hasTOFProng2, bool);                            //! Flag to check if prong2 has a TOF match
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, int);                                       //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, int);                 //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float);                               //! TPC chi2
DECLARE_SOA_COLUMN(ItsChi2NCl, itsChi2NCl, float);                               //! ITS chi2
DECLARE_SOA_COLUMN(ItsNClsProngMin, itsNClsProngMin, int);                       //! minimum value of number of ITS clusters for the decay daughter tracks
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsProngMin, tpcNClsCrossedRowsProngMin, int); //! minimum value of number of TPC crossed rows for the decay daughter tracks
DECLARE_SOA_COLUMN(TpcChi2NClProngMax, tpcChi2NClProngMax, float);               //! maximum value of TPC chi2 for the decay daughter tracks
DECLARE_SOA_COLUMN(PtProngMin, ptProngMin, float);                               //! minimum value of transverse momentum for the decay daughter tracks
DECLARE_SOA_COLUMN(AbsEtaProngMin, absEtaProngMin, float);                       //! minimum value of absolute pseudorapidity for the decay daughter tracks

// dynamic columns
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! transverse momentum
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! azimuthal angle
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! pseudorapidity
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array<float, 3>{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //!
                           [](float pxProng0, float pyProng0) -> float { return RecoDecay::pt(pxProng0, pyProng0); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float pxProng1, float pyProng1) -> float { return RecoDecay::pt(pxProng1, pyProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProng2, ptProng2, //!
                           [](float pxProng2, float pyProng2) -> float { return RecoDecay::pt(pxProng2, pyProng2); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaProng0, etaProng0, //!
                           [](float pxProng0, float pyProng0, float pzProng0) -> float { return RecoDecay::eta(std::array<float, 3>{pxProng0, pyProng0, pzProng0}); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaProng1, etaProng1, //!
                           [](float pxProng1, float pyProng1, float pzProng1) -> float { return RecoDecay::eta(std::array<float, 3>{pxProng1, pyProng1, pzProng1}); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaProng2, etaProng2, //!
                           [](float pxProng2, float pyProng2, float pzProng2) -> float { return RecoDecay::eta(std::array<float, 3>{pxProng2, pyProng2, pzProng2}); });
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector, //! 3-momentum vector
                           [](float px, float py, float pz) -> std::array<float, 3> { return {px, py, pz}; });
} // namespace hf_track_vars_reduced

namespace hf_b_to_jpsi_track_vars_reduced
{
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! transverse momentum
                           [](float signed1Pt) -> float { return std::abs(signed1Pt) <= o2::constants::math::Almost0 ? o2::constants::math::VeryBig : 1.f / std::abs(signed1Pt); });
} // namespace hf_b_to_jpsi_track_vars_reduced

namespace hf_track_pid_reduced
{
DECLARE_SOA_COLUMN(TPCNSigmaPiProng0, tpcNSigmaPiProng0, float); //! NsigmaTPCPi for prong0, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TPCNSigmaPiProng1, tpcNSigmaPiProng1, float); //! NsigmaTPCPi for prong1, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TPCNSigmaPiProng2, tpcNSigmaPiProng2, float); //! NsigmaTPCPi for prong2, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TPCNSigmaKaProng0, tpcNSigmaKaProng0, float); //! NsigmaTPCKa for prong0, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TPCNSigmaKaProng1, tpcNSigmaKaProng1, float); //! NsigmaTPCKa for prong1, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TPCNSigmaKaProng2, tpcNSigmaKaProng2, float); //! NsigmaTPCKa for prong2, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TPCNSigmaPrProng0, tpcNSigmaPrProng0, float); //! NsigmaTPCPr for prong0, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TPCNSigmaPrProng1, tpcNSigmaPrProng1, float); //! NsigmaTPCPr for prong1, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TPCNSigmaPrProng2, tpcNSigmaPrProng2, float); //! NsigmaTPCPr for prong2, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaPiProng0, tofNSigmaPiProng0, float); //! NsigmaTOFPi for prong0, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaPiProng1, tofNSigmaPiProng1, float); //! NsigmaTOFPi for prong1, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaPiProng2, tofNSigmaPiProng2, float); //! NsigmaTOFPi for prong2, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaKaProng0, tofNSigmaKaProng0, float); //! NsigmaTOFKa for prong0, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaKaProng1, tofNSigmaKaProng1, float); //! NsigmaTOFKa for prong1, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaKaProng2, tofNSigmaKaProng2, float); //! NsigmaTOFKa for prong2, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaPrProng0, tofNSigmaPrProng0, float); //! NsigmaTOFPr for prong0, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaPrProng1, tofNSigmaPrProng1, float); //! NsigmaTOFPr for prong1, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaPrProng2, tofNSigmaPrProng2, float); //! NsigmaTOFPr for prong2, o2-linter: disable=name/o2-column (written to disk)
// dynamic columns
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPi, tpcTofNSigmaPi, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPi, float tofNSigmaPi) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi, tofNSigmaPi); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaKa, tpcTofNSigmaKa, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPi, float tofNSigmaPi) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi, tofNSigmaPi); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPr, tpcTofNSigmaPr, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPi, float tofNSigmaPi) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi, tofNSigmaPi); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPiProng0, tpcTofNSigmaPiProng0, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPi, float tofNSigmaPi) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi, tofNSigmaPi); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPiProng1, tpcTofNSigmaPiProng1, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPi, float tofNSigmaPi) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi, tofNSigmaPi); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPiProng2, tpcTofNSigmaPiProng2, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPi, float tofNSigmaPi) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi, tofNSigmaPi); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaKaProng0, tpcTofNSigmaKaProng0, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaKa, float tofNSigmaKa) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa, tofNSigmaKa); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaKaProng1, tpcTofNSigmaKaProng1, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaKa, float tofNSigmaKa) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa, tofNSigmaKa); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaKaProng2, tpcTofNSigmaKaProng2, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaKa, float tofNSigmaKa) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa, tofNSigmaKa); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPrProng0, tpcTofNSigmaPrProng0, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPr, float tofNSigmaPr) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr, tofNSigmaPr); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPrProng1, tpcTofNSigmaPrProng1, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPr, float tofNSigmaPr) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr, tofNSigmaPr); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPrProng2, tpcTofNSigmaPrProng2, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPr, float tofNSigmaPr) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr, tofNSigmaPr); });
} // namespace hf_track_pid_reduced

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRedTrackBases, "AOD", "HFREDTRACKBASE", //! Table with track information for reduced workflow
                  soa::Index<>,
                  hf_track_index_reduced::TrackId,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  hf_track_vars_reduced::ItsNCls,
                  hf_track_vars_reduced::TpcNClsCrossedRows,
                  hf_track_vars_reduced::TpcChi2NCl,
                  aod::track::Px<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Py<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>,
                  aod::track::PVector<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha, aod::track::Tgl>,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfRedTracksCov, "AOD", "HFREDTRACKCOV", //! Table with track covariance information for reduced workflow
                  soa::Index<>,
                  HFTRACKPARCOV_COLUMNS,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfRedTracksMom, "AOD", "HFREDTRACKMOM", //! Table with track momentum information for reduced workflow
                  soa::Index<>,
                  hf_track_vars_reduced::Px,
                  hf_track_vars_reduced::Py,
                  hf_track_vars_reduced::Pz,
                  hf_track_vars_reduced::Sign);

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRedBach0Bases, "AOD", "HFREDBACH0BASE", //! Table with track information for reduced workflow
                  soa::Index<>,
                  hf_track_index_reduced::TrackId,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  hf_b_to_jpsi_track_vars_reduced::Pt<aod::track::Signed1Pt>,
                  hf_track_vars_reduced::ItsNCls,
                  hf_track_vars_reduced::TpcNClsCrossedRows,
                  hf_track_vars_reduced::TpcChi2NCl,
                  hf_track_vars_reduced::ItsChi2NCl,
                  hf_track_vars_reduced::HasTPC,
                  hf_track_vars_reduced::HasTOF,
                  pidtpc::TPCNSigmaPi,
                  pidtof::TOFNSigmaPi,
                  pidtpc::TPCNSigmaKa,
                  pidtof::TOFNSigmaKa,
                  pidtpc::TPCNSigmaPr,
                  pidtof::TOFNSigmaPr,
                  hf_track_pid_reduced::TPCTOFNSigmaPi<pidtpc::TPCNSigmaPi, pidtof::TOFNSigmaPi>,
                  hf_track_pid_reduced::TPCTOFNSigmaKa<pidtpc::TPCNSigmaKa, pidtof::TOFNSigmaKa>,
                  hf_track_pid_reduced::TPCTOFNSigmaPr<pidtpc::TPCNSigmaPr, pidtof::TOFNSigmaPr>,
                  aod::track::Px<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Py<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>,
                  aod::track::PVector<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha, aod::track::Tgl>);

DECLARE_SOA_TABLE(HfRedBach0Cov, "AOD", "HFREDBACH0COV", //! Table with track covariance information for reduced workflow
                  soa::Index<>,
                  HFTRACKPARCOV_COLUMNS);

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRedBach1Bases, "AOD", "HFREDBACH1BASE", //! Table with track information for reduced workflow
                  soa::Index<>,
                  hf_track_index_reduced::TrackId,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  hf_b_to_jpsi_track_vars_reduced::Pt<aod::track::Signed1Pt>,
                  hf_track_vars_reduced::ItsNCls,
                  hf_track_vars_reduced::TpcNClsCrossedRows,
                  hf_track_vars_reduced::TpcChi2NCl,
                  hf_track_vars_reduced::ItsChi2NCl,
                  hf_track_vars_reduced::HasTPC,
                  hf_track_vars_reduced::HasTOF,
                  pidtpc::TPCNSigmaPi,
                  pidtof::TOFNSigmaPi,
                  pidtpc::TPCNSigmaKa,
                  pidtof::TOFNSigmaKa,
                  pidtpc::TPCNSigmaPr,
                  pidtof::TOFNSigmaPr,
                  hf_track_pid_reduced::TPCTOFNSigmaPi<pidtpc::TPCNSigmaPi, pidtof::TOFNSigmaPi>,
                  hf_track_pid_reduced::TPCTOFNSigmaKa<pidtpc::TPCNSigmaKa, pidtof::TOFNSigmaKa>,
                  hf_track_pid_reduced::TPCTOFNSigmaPr<pidtpc::TPCNSigmaPr, pidtof::TOFNSigmaPr>,
                  aod::track::Px<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Py<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>,
                  aod::track::PVector<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha, aod::track::Tgl>);

DECLARE_SOA_TABLE(HfRedBach1Cov, "AOD", "HFREDBACH1COV", //! Table with track covariance information for reduced workflow
                  soa::Index<>,
                  HFTRACKPARCOV_COLUMNS);

// table with all attributes needed to call statusTpcAndTof() in the selector task
DECLARE_SOA_TABLE(HfRedTracksPid, "AOD", "HFREDTRACKPID", //! Table with PID track information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_vars_reduced::HasTPC,
                  hf_track_vars_reduced::HasTOF,
                  pidtpc::TPCNSigmaPi,
                  pidtof::TOFNSigmaPi,
                  hf_track_pid_reduced::TPCTOFNSigmaPi<pidtpc::TPCNSigmaPi, pidtof::TOFNSigmaPi>);

DECLARE_SOA_EXTENDED_TABLE_USER(HfRedTracksExt, HfRedTrackBases, "HFREDTRACKEXT", //! Track parameters at collision vertex
                                aod::track::Pt);
DECLARE_SOA_EXTENDED_TABLE_USER(HfRedBach0Ext, HfRedBach0Bases, "HFREDBACH0EXT", //! Track parameters at collision vertex
                                aod::track::Pt);
DECLARE_SOA_EXTENDED_TABLE_USER(HfRedBach1Ext, HfRedBach1Bases, "HFREDBACH1EXT", //! Track parameters at collision vertex
                                aod::track::Pt);

using HfRedTracks = HfRedTracksExt;
using HfRedBach0Tracks = HfRedBach0Bases;
using HfRedBach1Tracks = HfRedBach1Bases;

namespace hf_charm_cand_reduced
{
DECLARE_SOA_COLUMN(InvMassHypo0, invMassHypo0, float);                           //! Invariant mass of candidate in GeV/c2 (mass hypothesis 0)
DECLARE_SOA_COLUMN(InvMassHypo1, invMassHypo1, float);                           //! Invariant mass of candidate in GeV/c2 (mass hypothesis 1)
DECLARE_SOA_COLUMN(MlScoreBkgMassHypo0, mlScoreBkgMassHypo0, float);             //! ML score for background class (mass hypothesis 0)
DECLARE_SOA_COLUMN(MlScorePromptMassHypo0, mlScorePromptMassHypo0, float);       //! ML score for prompt class (mass hypothesis 0)
DECLARE_SOA_COLUMN(MlScoreNonpromptMassHypo0, mlScoreNonpromptMassHypo0, float); //! ML score for non-prompt class (mass hypothesis 0)
DECLARE_SOA_COLUMN(MlScoreBkgMassHypo1, mlScoreBkgMassHypo1, float);             //! ML score for background class (mass hypothesis 1)
DECLARE_SOA_COLUMN(MlScorePromptMassHypo1, mlScorePromptMassHypo1, float);       //! ML score for prompt class (mass hypothesis 1)
DECLARE_SOA_COLUMN(MlScoreNonpromptMassHypo1, mlScoreNonpromptMassHypo1, float); //! ML score for non-prompt class (mass hypothesis 1)
} // namespace hf_charm_cand_reduced

namespace hf_jpsi_cand_reduced
{
DECLARE_SOA_COLUMN(ProngPosId, prongPosId, int);             //! Original track index
DECLARE_SOA_COLUMN(ProngNegId, prongNegId, int);             //! Original track index
DECLARE_SOA_COLUMN(HfRedCollisionId, hfRedCollisionId, int); //! Collision index
DECLARE_SOA_COLUMN(M, m, float);                             //! Invariant mass of candidate in GeV/c2

DECLARE_SOA_COLUMN(ItsNClsDauPos, itsNClsDauPos, int);                       //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsDauPos, tpcNClsCrossedRowsDauPos, int); //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TpcChi2NClDauPos, tpcChi2NClDauPos, float);               //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(ItsChi2NClDauPos, itsChi2NClDauPos, float);               //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(ItsNClsDauNeg, itsNClsDauNeg, int);                       //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsDauNeg, tpcNClsCrossedRowsDauNeg, int); //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TpcChi2NClDauNeg, tpcChi2NClDauNeg, float);               //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(ItsChi2NClDauNeg, itsChi2NClDauNeg, float);               //! ITS chi2 / Number of clusters

DECLARE_SOA_COLUMN(XDauPos, xDauPos, float);                 //! x
DECLARE_SOA_COLUMN(XDauNeg, xDauNeg, float);                 //! x
DECLARE_SOA_COLUMN(YDauPos, yDauPos, float);                 //! y
DECLARE_SOA_COLUMN(YDauNeg, yDauNeg, float);                 //! y
DECLARE_SOA_COLUMN(ZDauPos, zDauPos, float);                 //! z
DECLARE_SOA_COLUMN(ZDauNeg, zDauNeg, float);                 //! z
DECLARE_SOA_COLUMN(AlphaDauPos, alphaDauPos, float);         //! alpha of the J/Psi positive decay daughter
DECLARE_SOA_COLUMN(AlphaDauNeg, alphaDauNeg, float);         //! alpha of the J/Psi negative decay daughter
DECLARE_SOA_COLUMN(SnpDauPos, snpDauPos, float);             //! snp of the J/Psi positive decay daughter
DECLARE_SOA_COLUMN(SnpDauNeg, snpDauNeg, float);             //! snp of the J/Psi negative decay daughter
DECLARE_SOA_COLUMN(TglDauPos, tglDauPos, float);             //! tgl of the J/Psi positive decay daughter
DECLARE_SOA_COLUMN(TglDauNeg, tglDauNeg, float);             //! tgl of the J/Psi negative decay daughter
DECLARE_SOA_COLUMN(Signed1PtDauPos, signed1PtDauPos, float); //! signed1Pt of the J/Psi positive decay daughter
DECLARE_SOA_COLUMN(Signed1PtDauNeg, signed1PtDauNeg, float); //! signed1Pt of the J/Psi negative decay daughter

DECLARE_SOA_DYNAMIC_COLUMN(PxDauPos, pxDauPos, //! Momentum in x-direction in GeV/c
                           [](float signed1Pt, float snp, float alpha) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             // FIXME: GCC & clang should optimize to sincosf
                             float cs = cosf(alpha), sn = sinf(alpha);
                             auto r = std::sqrt((1.f - snp) * (1.f + snp));
                             return pt * (r * cs - snp * sn);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(PyDauPos, pyDauPos, //! Momentum in y-direction in GeV/c
                           [](float signed1Pt, float snp, float alpha) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             // FIXME: GCC & clang should optimize to sincosf
                             float cs = cosf(alpha), sn = sinf(alpha);
                             auto r = std::sqrt((1.f - snp) * (1.f + snp));
                             return pt * (snp * cs + r * sn);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(PzDauPos, pzDauPos, //! Momentum in z-direction in GeV/c
                           [](float signed1Pt, float tgl) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             return pt * tgl;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(PxDauNeg, pxDauNeg, //! Momentum in x-direction in GeV/c
                           [](float signed1Pt, float snp, float alpha) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             // FIXME: GCC & clang should optimize to sincosf
                             float cs = cosf(alpha), sn = sinf(alpha);
                             auto r = std::sqrt((1.f - snp) * (1.f + snp));
                             return pt * (r * cs - snp * sn);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(PyDauNeg, pyDauNeg, //! Momentum in y-direction in GeV/c
                           [](float signed1Pt, float snp, float alpha) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             // FIXME: GCC & clang should optimize to sincosf
                             float cs = cosf(alpha), sn = sinf(alpha);
                             auto r = std::sqrt((1.f - snp) * (1.f + snp));
                             return pt * (snp * cs + r * sn);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(PzDauNeg, pzDauNeg, //! Momentum in z-direction in GeV/c
                           [](float signed1Pt, float tgl) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             return pt * tgl;
                           });

// Covariance matrix of the J/Psi positive decay daughter
DECLARE_SOA_COLUMN(CYYDauPos, cYYDauPos, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CZYDauPos, cZYDauPos, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CZZDauPos, cZZDauPos, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpYDauPos, cSnpYDauPos, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpZDauPos, cSnpZDauPos, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpSnpDauPos, cSnpSnpDauPos, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(CTglYDauPos, cTglYDauPos, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CTglZDauPos, cTglZDauPos, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CTglSnpDauPos, cTglSnpDauPos, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(CTglTglDauPos, cTglTglDauPos, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtYDauPos, c1PtYDauPos, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtZDauPos, c1PtZDauPos, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtSnpDauPos, c1PtSnpDauPos, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtTglDauPos, c1PtTglDauPos, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1Pt21Pt2DauPos, c1Pt21Pt2DauPos, float); //! Covariance matrix

// Covariance matrix of the J/Psi negative decay daughter
DECLARE_SOA_COLUMN(CYYDauNeg, cYYDauNeg, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CZYDauNeg, cZYDauNeg, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CZZDauNeg, cZZDauNeg, float);             //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpYDauNeg, cSnpYDauNeg, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpZDauNeg, cSnpZDauNeg, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CSnpSnpDauNeg, cSnpSnpDauNeg, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(CTglYDauNeg, cTglYDauNeg, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CTglZDauNeg, cTglZDauNeg, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(CTglSnpDauNeg, cTglSnpDauNeg, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(CTglTglDauNeg, cTglTglDauNeg, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtYDauNeg, c1PtYDauNeg, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtZDauNeg, c1PtZDauNeg, float);         //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtSnpDauNeg, c1PtSnpDauNeg, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1PtTglDauNeg, c1PtTglDauNeg, float);     //! Covariance matrix
DECLARE_SOA_COLUMN(C1Pt21Pt2DauNeg, c1Pt21Pt2DauNeg, float); //! Covariance matrix
} // namespace hf_jpsi_cand_reduced

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRed2Prongs, "AOD", "HFRED2PRONG", //! Table with 2prong candidate information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_charm_cand_reduced::InvMassHypo0, hf_charm_cand_reduced::InvMassHypo1,
                  hf_track_vars_reduced::PtProngMin, hf_track_vars_reduced::AbsEtaProngMin,
                  hf_track_vars_reduced::ItsNClsProngMin, hf_track_vars_reduced::TpcNClsCrossedRowsProngMin, hf_track_vars_reduced::TpcChi2NClProngMax,
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
DECLARE_SOA_TABLE(HfRedSoftPiBases, "AOD", "HFREDSOFTPIBASE", //! Table with track information for reduced workflow
                  soa::Index<>,
                  hf_track_index_reduced::TrackId,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  hf_track_vars_reduced::ItsNCls,
                  hf_track_vars_reduced::TpcNClsCrossedRows,
                  hf_track_vars_reduced::TpcChi2NCl,
                  aod::track::Px<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Py<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>,
                  aod::track::PVector<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha, aod::track::Tgl>);

DECLARE_SOA_TABLE(HfRedSoftPiCov, "AOD", "HFREDSOFTPICOV", //! Table with track covariance information for reduced workflow
                  soa::Index<>,
                  HFTRACKPARCOV_COLUMNS,
                  o2::soa::Marker<2>);

DECLARE_SOA_TABLE(HfRedSoftPiPid, "AOD", "HFREDSOFTPIPID",
                  soa::Index<>,
                  hf_cand_dstar::TPCNSigmaPiSoftPi,
                  hf_cand_dstar::TOFNSigmaPiSoftPi,
                  hf_track_vars_reduced::HasTOF,
                  hf_track_vars_reduced::HasTPC,
                  hf_cand_dstar::TPCTOFNSigmaPiSoftPi<hf_cand_dstar::TPCNSigmaPiSoftPi, hf_cand_dstar::TOFNSigmaPiSoftPi>)

namespace hf_track_index_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(SoftPi, softPi, int, HfRedSoftPiBases, ""); //! ReducedCollision index
}; // namespace hf_track_index_reduced

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRed3Prongs, "AOD", "HFRED3PRONG", //! Table with 3prong candidate information for reduced workflow
                  o2::soa::Index<>,
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id, hf_track_index_reduced::Prong2Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  HFTRACKPAR_COLUMNS,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_charm_cand_reduced::InvMassHypo0, hf_charm_cand_reduced::InvMassHypo1,
                  hf_track_vars_reduced::PtProngMin, hf_track_vars_reduced::AbsEtaProngMin,
                  hf_track_vars_reduced::ItsNClsProngMin, hf_track_vars_reduced::TpcNClsCrossedRowsProngMin, hf_track_vars_reduced::TpcChi2NClProngMax,
                  aod::track::Px<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Py<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha>,
                  aod::track::Pz<aod::track::Signed1Pt, track::Tgl>,
                  aod::track::PVector<aod::track::Signed1Pt, aod::track::Snp, aod::track::Alpha, aod::track::Tgl>);

DECLARE_SOA_TABLE(HfRed3ProngsCov, "AOD", "HFRED3PRONGSCOV", //! Table with 3prong candidate covariance for reduced workflow
                  o2::soa::Index<>,
                  HFTRACKPARCOV_COLUMNS,
                  o2::soa::Marker<2>);

DECLARE_SOA_TABLE(HfRed3ProngsMl_000, "AOD", "HFRED3PRONGML", //! Table with 3prong candidate ML scores
                  hf_charm_cand_reduced::MlScoreBkgMassHypo0,
                  hf_charm_cand_reduced::MlScorePromptMassHypo0,
                  hf_charm_cand_reduced::MlScoreNonpromptMassHypo0);

DECLARE_SOA_TABLE_VERSIONED(HfRed3ProngsMl_001, "AOD", "HFRED3PRONGML", 1, //! Table with 3prong candidate ML scores (format for 2 mass hypotheses needed for Ds and Lc)
                            hf_charm_cand_reduced::MlScoreBkgMassHypo0,
                            hf_charm_cand_reduced::MlScorePromptMassHypo0,
                            hf_charm_cand_reduced::MlScoreNonpromptMassHypo0,
                            hf_charm_cand_reduced::MlScoreBkgMassHypo1,
                            hf_charm_cand_reduced::MlScorePromptMassHypo1,
                            hf_charm_cand_reduced::MlScoreNonpromptMassHypo1,
                            o2::soa::Marker<1>);

using HfRed3ProngsMl = HfRed3ProngsMl_001;

DECLARE_SOA_TABLE(HfRedMomDDaugs, "AOD", "HFREDMOMDDAUGS", //! Table with 2prong candidate ML scores
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::PxProng2,
                  hf_cand::PyProng2,
                  hf_cand::PzProng2);

// CAREFUL: need to follow convention [Name = Description + 's'] in DECLARE_SOA_TABLE(Name, "AOD", Description)
// to call DECLARE_SOA_INDEX_COLUMN_FULL later on
DECLARE_SOA_TABLE(HfRedJpsis, "AOD", "HFREDJPSI", //! Table with J/Psi candidate information for reduced workflow
                  o2::soa::Index<>,
                  hf_jpsi_cand_reduced::ProngPosId,
                  hf_jpsi_cand_reduced::ProngNegId,
                  hf_track_index_reduced::HfRedCollisionId,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_jpsi_cand_reduced::M,
                  hf_jpsi_cand_reduced::ItsNClsDauPos,
                  hf_jpsi_cand_reduced::TpcNClsCrossedRowsDauPos,
                  hf_jpsi_cand_reduced::TpcChi2NClDauPos,
                  hf_jpsi_cand_reduced::ItsChi2NClDauPos,
                  hf_jpsi_cand_reduced::ItsNClsDauNeg,
                  hf_jpsi_cand_reduced::TpcNClsCrossedRowsDauNeg,
                  hf_jpsi_cand_reduced::TpcChi2NClDauNeg,
                  hf_jpsi_cand_reduced::ItsChi2NClDauNeg,
                  hf_jpsi_cand_reduced::XDauPos, hf_jpsi_cand_reduced::XDauNeg,
                  hf_jpsi_cand_reduced::YDauPos, hf_jpsi_cand_reduced::YDauNeg,
                  hf_jpsi_cand_reduced::ZDauPos, hf_jpsi_cand_reduced::ZDauNeg,
                  hf_jpsi_cand_reduced::AlphaDauPos, hf_jpsi_cand_reduced::AlphaDauNeg,
                  hf_jpsi_cand_reduced::SnpDauPos, hf_jpsi_cand_reduced::SnpDauNeg,
                  hf_jpsi_cand_reduced::TglDauPos, hf_jpsi_cand_reduced::TglDauNeg,
                  hf_jpsi_cand_reduced::Signed1PtDauPos, hf_jpsi_cand_reduced::Signed1PtDauNeg,
                  hf_jpsi_cand_reduced::PxDauPos<hf_jpsi_cand_reduced::Signed1PtDauPos, hf_jpsi_cand_reduced::SnpDauPos, hf_jpsi_cand_reduced::AlphaDauPos>,
                  hf_jpsi_cand_reduced::PxDauNeg<hf_jpsi_cand_reduced::Signed1PtDauNeg, hf_jpsi_cand_reduced::SnpDauNeg, hf_jpsi_cand_reduced::AlphaDauNeg>,
                  hf_jpsi_cand_reduced::PyDauPos<hf_jpsi_cand_reduced::Signed1PtDauPos, hf_jpsi_cand_reduced::SnpDauPos, hf_jpsi_cand_reduced::AlphaDauPos>,
                  hf_jpsi_cand_reduced::PyDauNeg<hf_jpsi_cand_reduced::Signed1PtDauNeg, hf_jpsi_cand_reduced::SnpDauNeg, hf_jpsi_cand_reduced::AlphaDauNeg>,
                  hf_jpsi_cand_reduced::PzDauPos<hf_jpsi_cand_reduced::Signed1PtDauPos, hf_jpsi_cand_reduced::TglDauPos>,
                  hf_jpsi_cand_reduced::PzDauNeg<hf_jpsi_cand_reduced::Signed1PtDauNeg, hf_jpsi_cand_reduced::TglDauNeg>);

DECLARE_SOA_TABLE(HfRedJpsiCov, "AOD", "HFREDJPSICOV", //! Table with J/Psi candidate covariance for reduced workflow
                  o2::soa::Index<>,
                  hf_jpsi_cand_reduced::CYYDauPos, hf_jpsi_cand_reduced::CYYDauNeg,
                  hf_jpsi_cand_reduced::CZYDauPos, hf_jpsi_cand_reduced::CZYDauNeg,
                  hf_jpsi_cand_reduced::CZZDauPos, hf_jpsi_cand_reduced::CZZDauNeg,
                  hf_jpsi_cand_reduced::CSnpYDauPos, hf_jpsi_cand_reduced::CSnpYDauNeg,
                  hf_jpsi_cand_reduced::CSnpZDauPos, hf_jpsi_cand_reduced::CSnpZDauNeg,
                  hf_jpsi_cand_reduced::CSnpSnpDauPos, hf_jpsi_cand_reduced::CSnpSnpDauNeg,
                  hf_jpsi_cand_reduced::CTglYDauPos, hf_jpsi_cand_reduced::CTglYDauNeg,
                  hf_jpsi_cand_reduced::CTglZDauPos, hf_jpsi_cand_reduced::CTglZDauNeg,
                  hf_jpsi_cand_reduced::CTglSnpDauPos, hf_jpsi_cand_reduced::CTglSnpDauNeg,
                  hf_jpsi_cand_reduced::CTglTglDauPos, hf_jpsi_cand_reduced::CTglTglDauNeg,
                  hf_jpsi_cand_reduced::C1PtYDauPos, hf_jpsi_cand_reduced::C1PtYDauNeg,
                  hf_jpsi_cand_reduced::C1PtZDauPos, hf_jpsi_cand_reduced::C1PtZDauNeg,
                  hf_jpsi_cand_reduced::C1PtSnpDauPos, hf_jpsi_cand_reduced::C1PtSnpDauNeg,
                  hf_jpsi_cand_reduced::C1PtTglDauPos, hf_jpsi_cand_reduced::C1PtTglDauNeg,
                  hf_jpsi_cand_reduced::C1Pt21Pt2DauPos, hf_jpsi_cand_reduced::C1Pt21Pt2DauNeg);

DECLARE_SOA_TABLE(HfRedPidDau0s_000, "AOD", "HFREDPIDDAU0", //!
                  hf_track_pid_reduced::TPCNSigmaPiProng0,
                  hf_track_pid_reduced::TOFNSigmaPiProng0,
                  hf_track_pid_reduced::TPCNSigmaKaProng0,
                  hf_track_pid_reduced::TOFNSigmaKaProng0,
                  hf_track_vars_reduced::HasTOFProng0,
                  hf_track_vars_reduced::HasTPCProng0,
                  hf_track_pid_reduced::TPCTOFNSigmaPiProng0<hf_track_pid_reduced::TPCNSigmaPiProng0, hf_track_pid_reduced::TOFNSigmaPiProng0>,
                  hf_track_pid_reduced::TPCTOFNSigmaKaProng0<hf_track_pid_reduced::TPCNSigmaKaProng0, hf_track_pid_reduced::TOFNSigmaKaProng0>);

DECLARE_SOA_TABLE(HfRedPidDau1s_000, "AOD", "HFREDPIDDAU1", //!
                  hf_track_pid_reduced::TPCNSigmaPiProng1,
                  hf_track_pid_reduced::TOFNSigmaPiProng1,
                  hf_track_pid_reduced::TPCNSigmaKaProng1,
                  hf_track_pid_reduced::TOFNSigmaKaProng1,
                  hf_track_vars_reduced::HasTOFProng1,
                  hf_track_vars_reduced::HasTPCProng1,
                  hf_track_pid_reduced::TPCTOFNSigmaPiProng1<hf_track_pid_reduced::TPCNSigmaPiProng1, hf_track_pid_reduced::TOFNSigmaPiProng1>,
                  hf_track_pid_reduced::TPCTOFNSigmaKaProng1<hf_track_pid_reduced::TPCNSigmaKaProng1, hf_track_pid_reduced::TOFNSigmaKaProng1>);

DECLARE_SOA_TABLE(HfRedPidDau2s_000, "AOD", "HFREDPIDDAU2", //!
                  hf_track_pid_reduced::TPCNSigmaPiProng2,
                  hf_track_pid_reduced::TOFNSigmaPiProng2,
                  hf_track_pid_reduced::TPCNSigmaKaProng2,
                  hf_track_pid_reduced::TOFNSigmaKaProng2,
                  hf_track_vars_reduced::HasTOFProng2,
                  hf_track_vars_reduced::HasTPCProng2,
                  hf_track_pid_reduced::TPCTOFNSigmaPiProng2<hf_track_pid_reduced::TPCNSigmaPiProng2, hf_track_pid_reduced::TOFNSigmaPiProng2>,
                  hf_track_pid_reduced::TPCTOFNSigmaKaProng2<hf_track_pid_reduced::TPCNSigmaKaProng2, hf_track_pid_reduced::TOFNSigmaKaProng2>);

DECLARE_SOA_TABLE_VERSIONED(HfRedPidDau0s_001, "AOD", "HFREDPIDDAU0", 1, //!
                            hf_track_pid_reduced::TPCNSigmaPiProng0,
                            hf_track_pid_reduced::TOFNSigmaPiProng0,
                            hf_track_pid_reduced::TPCNSigmaKaProng0,
                            hf_track_pid_reduced::TOFNSigmaKaProng0,
                            hf_track_pid_reduced::TPCNSigmaPrProng0,
                            hf_track_pid_reduced::TOFNSigmaPrProng0,
                            hf_track_vars_reduced::HasTOFProng0,
                            hf_track_vars_reduced::HasTPCProng0,
                            hf_track_pid_reduced::TPCTOFNSigmaPiProng0<hf_track_pid_reduced::TPCNSigmaPiProng0, hf_track_pid_reduced::TOFNSigmaPiProng0>,
                            hf_track_pid_reduced::TPCTOFNSigmaKaProng0<hf_track_pid_reduced::TPCNSigmaKaProng0, hf_track_pid_reduced::TOFNSigmaKaProng0>,
                            hf_track_pid_reduced::TPCTOFNSigmaPrProng0<hf_track_pid_reduced::TPCNSigmaPrProng0, hf_track_pid_reduced::TOFNSigmaPrProng0>,
                            o2::soa::Marker<1>);

DECLARE_SOA_TABLE_VERSIONED(HfRedPidDau1s_001, "AOD", "HFREDPIDDAU1", 1, //!
                            hf_track_pid_reduced::TPCNSigmaPiProng1,
                            hf_track_pid_reduced::TOFNSigmaPiProng1,
                            hf_track_pid_reduced::TPCNSigmaKaProng1,
                            hf_track_pid_reduced::TOFNSigmaKaProng1,
                            hf_track_pid_reduced::TPCNSigmaPrProng1,
                            hf_track_pid_reduced::TOFNSigmaPrProng1,
                            hf_track_vars_reduced::HasTOFProng1,
                            hf_track_vars_reduced::HasTPCProng1,
                            hf_track_pid_reduced::TPCTOFNSigmaPiProng1<hf_track_pid_reduced::TPCNSigmaPiProng1, hf_track_pid_reduced::TOFNSigmaPiProng1>,
                            hf_track_pid_reduced::TPCTOFNSigmaKaProng1<hf_track_pid_reduced::TPCNSigmaKaProng1, hf_track_pid_reduced::TOFNSigmaKaProng1>,
                            hf_track_pid_reduced::TPCTOFNSigmaPrProng1<hf_track_pid_reduced::TPCNSigmaPrProng1, hf_track_pid_reduced::TOFNSigmaPrProng1>,
                            o2::soa::Marker<1>);

DECLARE_SOA_TABLE_VERSIONED(HfRedPidDau2s_001, "AOD", "HFREDPIDDAU2", 1, //!
                            hf_track_pid_reduced::TPCNSigmaPiProng2,
                            hf_track_pid_reduced::TOFNSigmaPiProng2,
                            hf_track_pid_reduced::TPCNSigmaKaProng2,
                            hf_track_pid_reduced::TOFNSigmaKaProng2,
                            hf_track_pid_reduced::TPCNSigmaPrProng2,
                            hf_track_pid_reduced::TOFNSigmaPrProng2,
                            hf_track_vars_reduced::HasTOFProng2,
                            hf_track_vars_reduced::HasTPCProng2,
                            hf_track_pid_reduced::TPCTOFNSigmaPiProng2<hf_track_pid_reduced::TPCNSigmaPiProng2, hf_track_pid_reduced::TOFNSigmaPiProng2>,
                            hf_track_pid_reduced::TPCTOFNSigmaKaProng2<hf_track_pid_reduced::TPCNSigmaKaProng2, hf_track_pid_reduced::TOFNSigmaKaProng2>,
                            hf_track_pid_reduced::TPCTOFNSigmaPrProng2<hf_track_pid_reduced::TPCNSigmaPrProng2, hf_track_pid_reduced::TOFNSigmaPrProng2>,
                            o2::soa::Marker<1>);

using HfRedPidDau0s = HfRedPidDau0s_001;
using HfRedPidDau1s = HfRedPidDau1s_001;
using HfRedPidDau2s = HfRedPidDau2s_001;

using HfRedPidDau0 = HfRedPidDau0s::iterator;
using HfRedPidDau1 = HfRedPidDau1s::iterator;
using HfRedPidDau2 = HfRedPidDau2s::iterator;

// Beauty candidates prongs
namespace hf_cand_b0_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed3Prongs, "_0");               //! Prong0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrackBases, "_1");            //! Prong1 index
DECLARE_SOA_INDEX_COLUMN_FULL(ProngD0, prongD0, int, HfRed2Prongs, "_0");             //! ProngD0 index
DECLARE_SOA_INDEX_COLUMN_FULL(ProngBachPi, prongBachPi, int, HfRedTrackBases, "_1");  //! ProngBachPi index
DECLARE_SOA_INDEX_COLUMN_FULL(ProngSoftPi, prongSoftPi, int, HfRedSoftPiBases, "_2"); //! ProngSoftPi index
DECLARE_SOA_COLUMN(Prong0MlScoreBkg, prong0MlScoreBkg, float);                        //! Bkg ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScorePrompt, prong0MlScorePrompt, float);                  //! Prompt ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScoreNonprompt, prong0MlScoreNonprompt, float);            //! Nonprompt ML score of the D daughter
} // namespace hf_cand_b0_reduced

DECLARE_SOA_TABLE(HfRedB0Prongs, "AOD", "HFREDB0PRONG", //! Table with B0 daughter indices
                  hf_cand_b0_reduced::Prong0Id, hf_cand_b0_reduced::Prong1Id);

DECLARE_SOA_TABLE(HfRedB0ProngDStars, "AOD", "HFREDB0PRONGDST", //! Table with B0 daughter indices
                  hf_cand_b0_reduced::ProngD0Id, hf_cand_b0_reduced::ProngBachPiId, hf_cand_b0_reduced::ProngSoftPiId);

DECLARE_SOA_TABLE(HfRedB0DpMls, "AOD", "HFREDB0DPML", //! Table with ML scores for the D+ daughter
                  hf_cand_b0_reduced::Prong0MlScoreBkg,
                  hf_cand_b0_reduced::Prong0MlScorePrompt,
                  hf_cand_b0_reduced::Prong0MlScoreNonprompt,
                  o2::soa::Marker<1>);

using HfRedCandB0 = soa::Join<HfCandB0Ext, HfRedB0Prongs>;
using HfRedCandB0DStar = soa::Join<HfCandB0DStExt, HfRedB0ProngDStars>;

namespace hf_cand_bplus_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed2Prongs, "_0");    //! Prong0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrackBases, "_1"); //! Prong1 index
DECLARE_SOA_INDEX_COLUMN_FULL(Jpsi, jpsi, int, HfRedJpsis, "_0");          //! J/Psi index
DECLARE_SOA_INDEX_COLUMN_FULL(BachKa, bachKa, int, HfRedBach0Bases, "_0"); //! J/Psi index
DECLARE_SOA_COLUMN(Prong0MlScoreBkg, prong0MlScoreBkg, float);             //! Bkg ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScorePrompt, prong0MlScorePrompt, float);       //! Prompt ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScoreNonprompt, prong0MlScoreNonprompt, float); //! Nonprompt ML score of the D daughter
} // namespace hf_cand_bplus_reduced

DECLARE_SOA_TABLE(HfRedBplusProngs, "AOD", "HFREDBPPRONG",
                  hf_cand_bplus_reduced::Prong0Id, hf_cand_bplus_reduced::Prong1Id);

DECLARE_SOA_TABLE(HfRedBplus2JpsiDaus, "AOD", "HFREDBP2JPSIDAU",
                  hf_cand_bplus_reduced::JpsiId, hf_cand_bplus_reduced::BachKaId);

DECLARE_SOA_TABLE(HfRedBplusD0Mls, "AOD", "HFREDBPLUSD0ML", //! Table with ML scores for the D0 daughter
                  hf_cand_bplus_reduced::Prong0MlScoreBkg,
                  hf_cand_bplus_reduced::Prong0MlScorePrompt,
                  hf_cand_bplus_reduced::Prong0MlScoreNonprompt,
                  o2::soa::Marker<1>);

using HfRedCandBplus = soa::Join<HfCandBplusExt, HfRedBplusProngs>;
using HfRedCandBplusToJpsiK = soa::Join<HfCandBpJPExt, HfRedBplus2JpsiDaus>;

namespace hf_cand_bs_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed3Prongs, "_0");          //! Prong0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrackBases, "_1");       //! Prong1 index
DECLARE_SOA_INDEX_COLUMN_FULL(Jpsi, jpsi, int, HfRedJpsis, "_0");                //! J/Psi index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0Phi, prong0Phi, int, HfRedBach0Bases, "_0"); //! J/Psi index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1Phi, prong1Phi, int, HfRedBach1Bases, "_0"); //! J/Psi index
DECLARE_SOA_COLUMN(Prong0MlScoreBkg, prong0MlScoreBkg, float);                   //! Bkg ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScorePrompt, prong0MlScorePrompt, float);             //! Prompt ML score of the D daughter
DECLARE_SOA_COLUMN(Prong0MlScoreNonprompt, prong0MlScoreNonprompt, float);       //! Nonprompt ML score of the D daughter
} // namespace hf_cand_bs_reduced

DECLARE_SOA_TABLE(HfRedBsProngs, "AOD", "HFREDBSPRONG", //! Table with Bs daughter indices
                  hf_cand_bs_reduced::Prong0Id, hf_cand_bs_reduced::Prong1Id);

DECLARE_SOA_TABLE(HfRedBs2JpsiDaus, "AOD", "HFREDBS2JPSIDAU",
                  hf_cand_bs_reduced::JpsiId, hf_cand_bs_reduced::Prong0PhiId, hf_cand_bs_reduced::Prong1PhiId);

DECLARE_SOA_TABLE(HfRedBsDsMls, "AOD", "HFREDBSDSML", //! Table with ML scores for the Ds daughter
                  hf_cand_bs_reduced::Prong0MlScoreBkg,
                  hf_cand_bs_reduced::Prong0MlScorePrompt,
                  hf_cand_bs_reduced::Prong0MlScoreNonprompt,
                  o2::soa::Marker<1>);

using HfRedCandBs = soa::Join<HfCandBsExt, HfRedBsProngs>;
using HfRedCandBsToJpsiPhi = soa::Join<HfCandBsJPExt, HfRedBs2JpsiDaus>;

namespace hf_cand_lb_reduced
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed3Prongs, "_0");    //! Prong0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrackBases, "_1"); //! Prong1 index
DECLARE_SOA_COLUMN(Prong0MlScoreBkg, prong0MlScoreBkg, float);             //! Bkg ML score of the Lc daughter
DECLARE_SOA_COLUMN(Prong0MlScorePrompt, prong0MlScorePrompt, float);       //! Prompt ML score of the Lc daughter
DECLARE_SOA_COLUMN(Prong0MlScoreNonprompt, prong0MlScoreNonprompt, float); //! Nonprompt ML score of the Lc daughter
} // namespace hf_cand_lb_reduced

DECLARE_SOA_TABLE(HfRedLbProngs, "AOD", "HFREDLBPRONG", //! Table with Lb daughter indices
                  hf_cand_lb_reduced::Prong0Id, hf_cand_lb_reduced::Prong1Id);

DECLARE_SOA_TABLE(HfRedLbLcMls, "AOD", "HFREDLBLCML", //! Table with ML scores for the Lc daughter
                  hf_cand_lb_reduced::Prong0MlScoreBkg,
                  hf_cand_lb_reduced::Prong0MlScorePrompt,
                  hf_cand_lb_reduced::Prong0MlScoreNonprompt,
                  o2::soa::Marker<1>);

using HfRedCandLb = soa::Join<HfCandLbExt, HfRedLbProngs>;

namespace hf_cand_mc_flag
{
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t); //! reconstruction level
}

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
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
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

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedDStarPis, "AOD", "HFMCRECREDDSTPI", //! Table with reconstructed MC information on DStarPi pairs for reduced workflow
                  hf_cand_b0_reduced::ProngD0Id,
                  hf_cand_b0_reduced::ProngBachPiId,
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_b0_mc::PtMother);

// Table with same size as HFCANDB0
DECLARE_SOA_TABLE(HfMcRecRedB0s, "AOD", "HFMCRECREDB0", //! Reconstruction-level MC information on B0 candidates for reduced workflow
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
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
                  hf_cand_mc_flag::FlagMcMatchGen,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_b0_mc::PtTrack,
                  hf_b0_mc::YTrack,
                  hf_b0_mc::EtaTrack,
                  hf_b0_mc::PtProng0,
                  hf_b0_mc::YProng0,
                  hf_b0_mc::EtaProng0,
                  hf_b0_mc::PtProng1,
                  hf_b0_mc::YProng1,
                  hf_b0_mc::EtaProng1,
                  hf_reduced_collision::HfCollisionRejectionMap,
                  cent::CentFT0C,
                  cent::CentFT0M);

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
DECLARE_SOA_COLUMN(PdgCodeCharmMother, pdgCodeCharmMother, int);   //! Pdg code of charm mother
DECLARE_SOA_COLUMN(PdgCodeProng0, pdgCodeProng0, int);             //! Pdg code of prong0
DECLARE_SOA_COLUMN(PdgCodeProng1, pdgCodeProng1, int);             //! Pdg code of prong1
DECLARE_SOA_COLUMN(PdgCodeProng2, pdgCodeProng2, int);             //! Pdg code of prong2
} // namespace hf_bplus_mc

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedD0Pis, "AOD", "HFMCRECREDD0PI", //! Table with reconstructed MC information on D0Pi(<-B+) pairs for reduced workflow
                  hf_cand_bplus_reduced::Prong0Id,
                  hf_cand_bplus_reduced::Prong1Id,
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_bplus_mc::PtMother);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedJPKs, "AOD", "HFMCRECREDJPK", //! Table with reconstructed MC information on J/PsiK(<-B+) pairs for reduced workflow
                  hf_cand_bplus_reduced::JpsiId,
                  hf_cand_bplus_reduced::BachKaId,
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_bplus_mc::PtMother);

// DECLARE_SOA_EXTENDED_TABLE_USER(ExTable, Tracks, "EXTABLE",
DECLARE_SOA_TABLE(HfMcCheckD0Pis, "AOD", "HFMCCHECKD0PI", //! Table with reconstructed MC information on D0Pi(<-B0) pairs for MC checks in reduced workflow
                  hf_bplus_mc::PdgCodeBeautyMother,
                  hf_bplus_mc::PdgCodeCharmMother,
                  hf_bplus_mc::PdgCodeProng0,
                  hf_bplus_mc::PdgCodeProng1,
                  hf_bplus_mc::PdgCodeProng2,
                  o2::soa::Marker<1>);

// Table with same size as HFCANDBPLUS
DECLARE_SOA_TABLE(HfMcRecRedBps, "AOD", "HFMCRECREDBP", //! Reconstruction-level MC information on B+ candidates for reduced workflow
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_bplus_mc::PtMother);

DECLARE_SOA_TABLE(HfMcCheckBps, "AOD", "HFMCCHECKBP", //! Table with reconstructed MC information on B+ candidates for MC checks in reduced workflow
                  hf_bplus_mc::PdgCodeBeautyMother,
                  hf_bplus_mc::PdgCodeCharmMother,
                  hf_bplus_mc::PdgCodeProng0,
                  hf_bplus_mc::PdgCodeProng1,
                  hf_bplus_mc::PdgCodeProng2,
                  o2::soa::Marker<2>);

DECLARE_SOA_TABLE(HfMcGenRedBps, "AOD", "HFMCGENREDBP", //! Generation-level MC information on B+ candidates for reduced workflow
                  hf_cand_mc_flag::FlagMcMatchGen,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_bplus_mc::PtTrack,
                  hf_bplus_mc::YTrack,
                  hf_bplus_mc::EtaTrack,
                  hf_bplus_mc::PtProng0,
                  hf_bplus_mc::YProng0,
                  hf_bplus_mc::EtaProng0,
                  hf_bplus_mc::PtProng1,
                  hf_bplus_mc::YProng1,
                  hf_bplus_mc::EtaProng1,
                  hf_reduced_collision::HfCollisionRejectionMap,
                  cent::CentFT0C,
                  cent::CentFT0M);

// store all configurables values used in the first part of the workflow
// so we can use them in the Bplus part
namespace hf_cand_bplus_config
{
DECLARE_SOA_COLUMN(MySelectionFlagD0, mySelectionFlagD0, int8_t);       //! Flag to filter selected D0 mesons
DECLARE_SOA_COLUMN(MySelectionFlagD0bar, mySelectionFlagD0bar, int8_t); //! Flag to filter selected D0 mesons
DECLARE_SOA_COLUMN(MyInvMassWindowD0Pi, myInvMassWindowD0Pi, float);    //! Half-width of the Bplus invariant-mass window in GeV/c2
DECLARE_SOA_COLUMN(MyInvMassWindowJpsiK, myInvMassWindowJpsiK, float);  //! Half-width of the Bplus invariant-mass window in GeV/c2
} // namespace hf_cand_bplus_config

DECLARE_SOA_TABLE(HfCandBpConfigs, "AOD", "HFCANDBPCONFIG", //! Table with configurables information for reduced workflow
                  hf_cand_bplus_config::MySelectionFlagD0,
                  hf_cand_bplus_config::MySelectionFlagD0bar,
                  hf_cand_bplus_config::MyInvMassWindowD0Pi);

DECLARE_SOA_TABLE(HfCfgBpToJpsi, "AOD", "HFCFGBPTOJPSI", //! Table with configurables information for reduced workflow
                  hf_cand_bplus_config::MyInvMassWindowJpsiK);

namespace hf_bs_mc
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
} // namespace hf_bs_mc

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedDsPis, "AOD", "HFMCRECREDDSPI", //! Table with reconstructed MC information on DsPi(<-Bs) pairs for reduced workflow
                  hf_cand_bs_reduced::Prong0Id,
                  hf_cand_bs_reduced::Prong1Id,
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_bs_mc::PtMother);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedJPPhis, "AOD", "HFMCRECREDJPPHI", //! Table with reconstructed MC information on DsPi(<-Bs) pairs for reduced workflow
                  hf_cand_bs_reduced::JpsiId,
                  hf_cand_bs_reduced::Prong0PhiId,
                  hf_cand_bs_reduced::Prong1PhiId,
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_bs_mc::PtMother);

// try with extended table ?
// DECLARE_SOA_EXTENDED_TABLE_USER(ExTable, Tracks, "EXTABLE",
DECLARE_SOA_TABLE(HfMcCheckDsPis, "AOD", "HFMCCHECKDSPI", //! Table with reconstructed MC information on DsPi(<-Bs) pairs for MC checks in reduced workflow
                  hf_bs_mc::PdgCodeBeautyMother,
                  hf_bs_mc::PdgCodeCharmMother,
                  hf_bs_mc::PdgCodeProng0,
                  hf_bs_mc::PdgCodeProng1,
                  hf_bs_mc::PdgCodeProng2,
                  hf_bs_mc::PdgCodeProng3,
                  o2::soa::Marker<1>);

// Table with same size as HFCANDBS
DECLARE_SOA_TABLE(HfMcRecRedBss, "AOD", "HFMCRECREDBS", //! Reconstruction-level MC information on Bs candidates for reduced workflow
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_bs_mc::PtMother);

DECLARE_SOA_TABLE(HfMcCheckBss, "AOD", "HFMCCHECKBS", //! Table with reconstructed MC information on Bs candidates for MC checks in reduced workflow
                  hf_bs_mc::PdgCodeBeautyMother,
                  hf_bs_mc::PdgCodeCharmMother,
                  hf_bs_mc::PdgCodeProng0,
                  hf_bs_mc::PdgCodeProng1,
                  hf_bs_mc::PdgCodeProng2,
                  hf_bs_mc::PdgCodeProng3,
                  o2::soa::Marker<2>);

DECLARE_SOA_TABLE(HfMcGenRedBss, "AOD", "HFMCGENREDBS", //! Generation-level MC information on Bs candidates for reduced workflow
                  hf_cand_mc_flag::FlagMcMatchGen,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_bs_mc::PtTrack,
                  hf_bs_mc::YTrack,
                  hf_bs_mc::EtaTrack,
                  hf_bs_mc::PtProng0,
                  hf_bs_mc::YProng0,
                  hf_bs_mc::EtaProng0,
                  hf_bs_mc::PtProng1,
                  hf_bs_mc::YProng1,
                  hf_bs_mc::EtaProng1,
                  hf_reduced_collision::HfCollisionRejectionMap,
                  cent::CentFT0C,
                  cent::CentFT0M);

// store all configurables values used in the first part of the workflow
// so we can use them in the Bs part
namespace hf_cand_bs_config
{
DECLARE_SOA_COLUMN(MySelectionFlagD, mySelectionFlagD, int8_t);            //! Flag to filter selected Ds mesons
DECLARE_SOA_COLUMN(MyInvMassWindowDPi, myInvMassWindowDPi, float);         //! Half-width of the Bs invariant-mass window in GeV/c2
DECLARE_SOA_COLUMN(MyInvMassWindowJpsiPhi, myInvMassWindowJpsiPhi, float); //! Half-width of the Bs invariant-mass window in GeV/c2
} // namespace hf_cand_bs_config

DECLARE_SOA_TABLE(HfCandBsConfigs, "AOD", "HFCANDBSCONFIG", //! Table with configurables information for reduced workflow
                  hf_cand_bs_config::MySelectionFlagD,
                  hf_cand_bs_config::MyInvMassWindowDPi);

DECLARE_SOA_TABLE(HfCfgBsToJpsis, "AOD", "HFCFGBSTOJPSI", //! Table with configurables information for reduced workflow
                  hf_cand_bs_config::MyInvMassWindowJpsiPhi);
namespace hf_lb_mc
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
} // namespace hf_lb_mc

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfMcRecRedLcPis, "AOD", "HFMCRECREDLCPI", //! Table with reconstructed MC information on LcPi(<-Lb) pairs for reduced workflow
                  hf_cand_lb_reduced::Prong0Id,
                  hf_cand_lb_reduced::Prong1Id,
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_lb_mc::PtMother);

DECLARE_SOA_TABLE(HfMcCheckLcPis, "AOD", "HFMCCHECKLCPI", //! Table with reconstructed MC information on LcPi(<-Lb) pairs for MC checks in reduced workflow
                  hf_lb_mc::PdgCodeBeautyMother,
                  hf_lb_mc::PdgCodeCharmMother,
                  hf_lb_mc::PdgCodeProng0,
                  hf_lb_mc::PdgCodeProng1,
                  hf_lb_mc::PdgCodeProng2,
                  hf_lb_mc::PdgCodeProng3,
                  o2::soa::Marker<1>);

// Table with same size as HFCANDLc
DECLARE_SOA_TABLE(HfMcRecRedLbs, "AOD", "HFMCRECREDLB", //! Reconstruction-level MC information on Lb candidates for reduced workflow
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagWrongCollision,
                  hf_cand_mc_flag::DebugMcRec,
                  hf_lb_mc::PtMother);

DECLARE_SOA_TABLE(HfMcCheckLbs, "AOD", "HFMCCHECKLB", //! Table with reconstructed MC information on Lb candidates for MC checks in reduced workflow
                  hf_lb_mc::PdgCodeBeautyMother,
                  hf_lb_mc::PdgCodeCharmMother,
                  hf_lb_mc::PdgCodeProng0,
                  hf_lb_mc::PdgCodeProng1,
                  hf_lb_mc::PdgCodeProng2,
                  hf_lb_mc::PdgCodeProng3,
                  o2::soa::Marker<2>);

DECLARE_SOA_TABLE(HfMcGenRedLbs, "AOD", "HFMCGENREDLB", //! Generation-level MC information on Lb candidates for reduced workflow
                  hf_cand_mc_flag::FlagMcMatchGen,
                  hf_lb_mc::PtTrack,
                  hf_lb_mc::YTrack,
                  hf_lb_mc::EtaTrack,
                  hf_lb_mc::PtProng0,
                  hf_lb_mc::YProng0,
                  hf_lb_mc::EtaProng0,
                  hf_lb_mc::PtProng1,
                  hf_lb_mc::YProng1,
                  hf_lb_mc::EtaProng1,
                  hf_reduced_collision::HfCollisionRejectionMap,
                  cent::CentFT0C,
                  cent::CentFT0M);

// store all configurables values used in the first part of the workflow
// so we can use them in the B0 part
namespace hf_cand_lb_config
{
DECLARE_SOA_COLUMN(MySelectionFlagLc, mySelectionFlagLc, int8_t);    //! Flag to filter selected Lc baryons
DECLARE_SOA_COLUMN(MyInvMassWindowLcPi, myInvMassWindowLcPi, float); //! Half-width of the Lb invariant-mass window in GeV/c2
} // namespace hf_cand_lb_config

DECLARE_SOA_TABLE(HfCandLbConfigs, "AOD", "HFCANDLBCONFIG", //! Table with configurables information for reduced workflow
                  hf_cand_lb_config::MySelectionFlagLc,
                  hf_cand_lb_config::MyInvMassWindowLcPi);

// Charm resonances analysis
namespace hf_reso_3_prong
{
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                      //! Integer with selected D candidate sign
DECLARE_SOA_COLUMN(ItsNClsSoftPi, itsNClsSoftPi, int);                       //! minimum value of number of ITS clusters for the decay daughter tracks
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsSoftPi, tpcNClsCrossedRowsSoftPi, int); //! minimum value of number of TPC crossed rows for the decay daughter tracks
DECLARE_SOA_COLUMN(TpcChi2NClSoftPi, tpcChi2NClSoftPi, float);               //! maximum value of TPC chi2 for the decay daughter tracks
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,                                           //!
                           [](float pxProng0, float pxProng1, float pxProng2) -> float { return 1.f * pxProng0 + 1.f * pxProng1 + 1.f * pxProng2; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pyProng0, float pyProng1, float pyProng2) -> float { return 1.f * pyProng0 + 1.f * pyProng1 + 1.f * pyProng2; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pzProng0, float pzProng1, float pzProng2) -> float { return 1.f * pzProng0 + 1.f * pzProng1 + 1.f * pzProng2; });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float pxProng0, float pxProng1, float pxProng2, float pyProng0, float pyProng1, float pyProng2) -> float { return RecoDecay::pt((1.f * pxProng0 + 1.f * pxProng1 + 1.f * pxProng2), (1.f * pyProng0 + 1.f * pyProng1 + 1.f * pyProng2)); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassDplus, invMassDplus,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus, constants::physics::MassPiPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2) -> std::array<float, 3> { return std::array{px0 + px1 + px2, py0 + py1 + py2, pz0 + pz1 + pz2}; });
} // namespace hf_reso_3_prong

namespace hf_reso_2_prong
{
DECLARE_SOA_COLUMN(SelFlagD0, selFlagD0, uint8_t); //! Integer with D0 selection flag: 1 = selected as D0, 2 = selected as D0bar, 3 = selected as D0 and D0bar
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,                 //!
                           [](float pxProng0, float pxProng1) -> float { return 1.f * pxProng0 + 1.f * pxProng1; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pyProng0, float pyProng1) -> float { return 1.f * pyProng0 + 1.f * pyProng1; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pzProng0, float pzProng1) -> float { return 1.f * pzProng0 + 1.f * pzProng1; });
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector,
                           [](float pxProng0, float pyProng0, float pzProng0, float pxProng1, float pyProng1, float pzProng1) -> std::array<float, 3> { return std::array{pxProng0 + pxProng1, pyProng0 + pyProng1, pzProng0 + pzProng1}; });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float pxProng0, float pxProng1, float pyProng0, float pyProng1) -> float { return RecoDecay::pt((1.f * pxProng0 + 1.f * pxProng1), (1.f * pyProng0 + 1.f * pyProng1)); });
} // namespace hf_reso_2_prong

namespace hf_reso_v0
{
DECLARE_SOA_COLUMN(Cpa, cpa, float);         //! Cosine of Pointing Angle of V0 candidate
DECLARE_SOA_COLUMN(Dca, dca, float);         //! DCA of V0 candidate
DECLARE_SOA_COLUMN(Radius, radius, float);   //! Radius of V0 candidate
DECLARE_SOA_COLUMN(V0Type, v0Type, uint8_t); //! Bitmap with mass hypothesis of the V0

DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //!
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
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector,
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> std::array<float, 3> { return std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}; });
} // namespace hf_reso_v0

DECLARE_SOA_TABLE(HfRedVzeros, "AOD", "HFREDVZERO", //! Table with V0 candidate information for resonances reduced workflow
                  o2::soa::Index<>,
                  // Indices
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  // Static
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_reso_v0::Cpa, hf_reso_v0::Dca,
                  hf_track_vars_reduced::ItsNClsProngMin, hf_track_vars_reduced::TpcNClsCrossedRowsProngMin, hf_track_vars_reduced::TpcChi2NClProngMax,
                  hf_reso_v0::V0Type,
                  // Dynamic
                  hf_reso_v0::Px<hf_cand::PxProng0, hf_cand::PxProng1>,
                  hf_reso_v0::Py<hf_cand::PyProng0, hf_cand::PyProng1>,
                  hf_reso_v0::Pz<hf_cand::PzProng0, hf_cand::PzProng1>,
                  hf_track_vars_reduced::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_track_vars_reduced::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_track_vars_reduced::EtaProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_track_vars_reduced::EtaProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_v0::InvMassK0s<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_v0::InvMassLambda<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_v0::InvMassAntiLambda<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_v0::V0Radius<hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,
                  hf_reso_v0::Pt<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PyProng0, hf_cand::PyProng1>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_v0::PVector<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>);

DECLARE_SOA_TABLE(HfRedTrkNoParams, "AOD", "HFREDTRKNOPARAM", //! Table with tracks without track parameters for resonances reduced workflow
                  o2::soa::Index<>,
                  // Indices
                  hf_track_index_reduced::TrackId,
                  hf_track_index_reduced::HfRedCollisionId,
                  // Static
                  hf_track_vars_reduced::Px,
                  hf_track_vars_reduced::Py,
                  hf_track_vars_reduced::Pz,
                  hf_track_vars_reduced::Sign,
                  pidtpc::TPCNSigmaPi,
                  pidtpc::TPCNSigmaKa,
                  pidtpc::TPCNSigmaPr,
                  pidtof::TOFNSigmaPi,
                  pidtof::TOFNSigmaKa,
                  pidtof::TOFNSigmaPr,
                  hf_track_vars_reduced::HasTOF,
                  hf_track_vars_reduced::HasTPC,
                  hf_track_vars_reduced::ItsNCls,
                  hf_track_vars_reduced::TpcNClsCrossedRows,
                  hf_track_vars_reduced::TpcChi2NCl,
                  // Dynamic
                  hf_track_vars_reduced::Pt<hf_track_vars_reduced::Px, hf_track_vars_reduced::Py>,
                  hf_track_vars_reduced::Eta<hf_track_vars_reduced::Px, hf_track_vars_reduced::Py, hf_track_vars_reduced::Pz>,
                  hf_track_vars_reduced::Phi<hf_track_vars_reduced::Px, hf_track_vars_reduced::Py>,
                  hf_track_pid_reduced::TPCTOFNSigmaPi<pidtpc::TPCNSigmaPi, pidtof::TOFNSigmaPi>,
                  hf_track_pid_reduced::TPCTOFNSigmaKa<pidtpc::TPCNSigmaKa, pidtof::TOFNSigmaKa>,
                  hf_track_pid_reduced::TPCTOFNSigmaPr<pidtpc::TPCNSigmaPr, pidtof::TOFNSigmaPr>,
                  hf_track_vars_reduced::PVector<hf_track_vars_reduced::Px, hf_track_vars_reduced::Py, hf_track_vars_reduced::Pz>);

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
                  hf_track_vars_reduced::ItsNClsProngMin, hf_track_vars_reduced::TpcNClsCrossedRowsProngMin, hf_track_vars_reduced::TpcChi2NClProngMax,
                  hf_reso_3_prong::Sign,
                  // Dynamic
                  hf_reso_3_prong::Px<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PxProng2>,
                  hf_reso_3_prong::Py<hf_cand::PyProng0, hf_cand::PyProng1, hf_cand::PyProng2>,
                  hf_reso_3_prong::Pz<hf_cand::PzProng0, hf_cand::PzProng1, hf_cand::PzProng2>,
                  hf_track_vars_reduced::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_track_vars_reduced::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_track_vars_reduced::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_track_vars_reduced::EtaProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_track_vars_reduced::EtaProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_track_vars_reduced::EtaProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_reso_3_prong::InvMassDplus<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_reso_3_prong::Pt<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PxProng2, hf_cand::PyProng0, hf_cand::PyProng1, hf_cand::PyProng2>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_reso_3_prong::PVector<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>);

DECLARE_SOA_TABLE(HfRed2PrNoTrks, "AOD", "HFRED2PRNOTRK", //! Table with 2 prong candidate information for resonances reduced workflow
                  o2::soa::Index<>,
                  // Indices
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  // Static
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_track_vars_reduced::ItsNClsProngMin, hf_track_vars_reduced::TpcNClsCrossedRowsProngMin, hf_track_vars_reduced::TpcChi2NClProngMax,
                  hf_reso_2_prong::SelFlagD0,
                  // Dynamic
                  hf_reso_2_prong::Px<hf_cand::PxProng0, hf_cand::PxProng1>,
                  hf_reso_2_prong::Py<hf_cand::PyProng0, hf_cand::PyProng1>,
                  hf_reso_2_prong::Pz<hf_cand::PzProng0, hf_cand::PzProng1>,
                  hf_track_vars_reduced::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_track_vars_reduced::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_track_vars_reduced::EtaProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_track_vars_reduced::EtaProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_2_prong::PVector<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_2_prong::Pt<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PyProng0, hf_cand::PyProng1>,
                  // InvMasses
                  hf_cand_dstar::InvMassD0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassD0Bar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>);

DECLARE_SOA_TABLE(HfRedDstarNoTrks, "AOD", "HFREDDSTARNOTRK", //! Table with 3 prong candidate information for resonances reduced workflow
                  o2::soa::Index<>,
                  // Indices
                  hf_track_index_reduced::Prong0Id, hf_track_index_reduced::Prong1Id, hf_track_index_reduced::Prong2Id,
                  hf_track_index_reduced::HfRedCollisionId,
                  // Static
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_track_vars_reduced::ItsNClsProngMin, hf_track_vars_reduced::TpcNClsCrossedRowsProngMin, hf_track_vars_reduced::TpcChi2NClProngMax,
                  hf_reso_3_prong::ItsNClsSoftPi, hf_reso_3_prong::TpcNClsCrossedRowsSoftPi, hf_reso_3_prong::TpcChi2NClSoftPi,
                  hf_reso_3_prong::Sign,
                  // Dynamic
                  hf_reso_3_prong::Px<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PxProng2>,
                  hf_reso_3_prong::Py<hf_cand::PyProng0, hf_cand::PyProng1, hf_cand::PyProng2>,
                  hf_reso_3_prong::Pz<hf_cand::PzProng0, hf_cand::PzProng1, hf_cand::PzProng2>,
                  hf_track_vars_reduced::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_track_vars_reduced::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_track_vars_reduced::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_track_vars_reduced::EtaProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_track_vars_reduced::EtaProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_track_vars_reduced::EtaProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_dstar::InvMassDstar<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassAntiDstar<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassD0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassD0Bar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_reso_3_prong::Pt<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PxProng2, hf_cand::PyProng0, hf_cand::PyProng1, hf_cand::PyProng2>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_reso_3_prong::PVector<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>);

namespace hf_reso_cand_reduced
{
DECLARE_SOA_COLUMN(InvMass, invMass, float);             //! Invariant mass in GeV/c2
DECLARE_SOA_COLUMN(InvMassProng0, invMassProng0, float); //! Invariant Mass of D daughter in GeV/c
DECLARE_SOA_COLUMN(InvMassProng1, invMassProng1, float); //! Invariant Mass of V0/Tr daughter in GeV/c
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                  //! Sign of the Resonance candidate
DECLARE_SOA_COLUMN(IsWrongSign, isWrongSign, int8_t);    //! Flag for wrong sign of the Resonance candidate, 1 = wrong sign, 0 = right sign

DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);     // flag for resonance decay channel classification reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchRecD, flagMcMatchRecD, int8_t);   // flag for D meson bachelor decay channel classification reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchChanD, flagMcMatchChanD, int8_t); // flag for D meson resonant channel classification reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);     // flag for decay channel classification generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, uint16_t);           // debug flag for mis-association at reconstruction level
DECLARE_SOA_COLUMN(Origin, origin, int8_t);                     // Flag for origin of MC particle 1=promt, 2=FD
DECLARE_SOA_COLUMN(SignD0, signD0, int8_t);                     // Sign of the D0 in the channels with D* -> D0 pi, needed in case of non-matched D*
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                        // Pt at generation level in GeV/c
DECLARE_SOA_COLUMN(InvMassGen, invMassGen, float);              //! Invariant mass at generation level in GeV/c2
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,                              //!
                           [](float pxProng0, float pxProng1, float pyProng0, float pyProng1) -> float { return RecoDecay::pt((1.f * pxProng0 + 1.f * pxProng1), (1.f * pyProng0 + 1.f * pyProng1)); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //!
                           [](float pxProng0, float pyProng0) -> float { return RecoDecay::pt(pxProng0, pyProng0); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float pxProng1, float pyProng1) -> float { return RecoDecay::pt(pxProng1, pyProng1); });
} // namespace hf_reso_cand_reduced

namespace hf_reso_3pr_v0
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed3PrNoTrks, "_0"); //! Prong0 index (D daughter)
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedVzeros, "_1");    //! Prong1 index (V0 daughter)
} // namespace hf_reso_3pr_v0
namespace hf_reso_dstar_v0
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRedDstarNoTrks, "_0"); //! Prong0 index (D daughter)
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedVzeros, "_1");      //! Prong1 index (V0 daughter)
} // namespace hf_reso_dstar_v0
namespace hf_reso_2pr_v0
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed2PrNoTrks, "_0"); //! Prong0 index (D daughter)
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedVzeros, "_1");    //! Prong1 index (V0 daughter)
} // namespace hf_reso_2pr_v0
namespace hf_reso_3pr_trk
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed3PrNoTrks, "_0");   //! Prong0 index (D daughter)
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrkNoParams, "_1"); //! Prong1 index (Track daughter)
} // namespace hf_reso_3pr_trk
namespace hf_reso_dstar_trk
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRedDstarNoTrks, "_0"); //! Prong0 index (D daughter)
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrkNoParams, "_1"); //! Prong1 index (Track daughter)
} // namespace hf_reso_dstar_trk
namespace hf_reso_2pr_trk
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfRed2PrNoTrks, "_0");   //! Prong0 index (D daughter)
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, HfRedTrkNoParams, "_1"); //! Prong1 index (Track daughter)
} // namespace hf_reso_2pr_trk

DECLARE_SOA_TABLE(HfCandCharmReso, "AOD", "HFCANDCHARMRESO", //! Table with Resonance candidate information for resonances reduced workflow
                  o2::soa::Index<>,
                  // Static
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_reso_cand_reduced::InvMass,
                  hf_reso_cand_reduced::InvMassProng0,
                  hf_reso_cand_reduced::InvMassProng1,
                  hf_reso_cand_reduced::Sign,
                  hf_reso_cand_reduced::IsWrongSign,
                  // Dynamic
                  hf_reso_cand_reduced::Pt<hf_cand::PxProng0, hf_cand::PxProng1, hf_cand::PyProng0, hf_cand::PyProng1>,
                  hf_reso_cand_reduced::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_reso_cand_reduced::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_reso_v0::Px<hf_cand::PxProng0, hf_cand::PxProng1>,
                  hf_reso_v0::Py<hf_cand::PyProng0, hf_cand::PyProng1>,
                  hf_reso_v0::Pz<hf_cand::PzProng0, hf_cand::PzProng1>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>);

DECLARE_SOA_TABLE(Hf3PrV0Ids, "AOD", "HF3PRV0ID",
                  hf_track_index_reduced::HfRedCollisionId,
                  hf_reso_3pr_v0::Prong0Id,
                  hf_reso_3pr_v0::Prong1Id);
DECLARE_SOA_TABLE(HfDstarV0Ids, "AOD", "HFDSTARV0ID",
                  hf_track_index_reduced::HfRedCollisionId,
                  hf_reso_dstar_v0::Prong0Id,
                  hf_reso_dstar_v0::Prong1Id);
DECLARE_SOA_TABLE(Hf2PrV0Ids, "AOD", "HF2PRV0ID",
                  hf_track_index_reduced::HfRedCollisionId,
                  hf_reso_2pr_v0::Prong0Id,
                  hf_reso_2pr_v0::Prong1Id);
DECLARE_SOA_TABLE(Hf3PrTrkIds, "AOD", "HF3PRTRKID",
                  hf_track_index_reduced::HfRedCollisionId,
                  hf_reso_3pr_trk::Prong0Id,
                  hf_reso_3pr_trk::Prong1Id);
DECLARE_SOA_TABLE(HfDstarTrkIds, "AOD", "HFDSTARTRKID",
                  hf_track_index_reduced::HfRedCollisionId,
                  hf_reso_dstar_trk::Prong0Id,
                  hf_reso_dstar_trk::Prong1Id);
DECLARE_SOA_TABLE(Hf2PrTrkIds, "AOD", "HF2PRTRKID",
                  hf_track_index_reduced::HfRedCollisionId,
                  hf_reso_2pr_trk::Prong0Id,
                  hf_reso_2pr_trk::Prong1Id);

// Tables for MC Resonance analysis
// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(Hf3PrV0McRec, "AOD", "HF3PRV0MCREC",
                  hf_reso_3pr_v0::Prong0Id,
                  hf_reso_3pr_v0::Prong1Id,
                  hf_reso_cand_reduced::FlagMcMatchRec,
                  hf_reso_cand_reduced::FlagMcMatchRecD,
                  hf_reso_cand_reduced::FlagMcMatchChanD,
                  hf_reso_cand_reduced::DebugMcRec,
                  hf_reso_cand_reduced::Origin,
                  hf_reso_cand_reduced::PtGen,
                  hf_reso_cand_reduced::InvMassGen,
                  hf_cand_mc_flag::NTracksDecayed,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfDstarV0McRec, "AOD", "HFDSTARV0MCREC",
                  hf_reso_dstar_v0::Prong0Id,
                  hf_reso_dstar_v0::Prong1Id,
                  hf_reso_cand_reduced::FlagMcMatchRec,
                  hf_reso_cand_reduced::FlagMcMatchRecD,
                  hf_reso_cand_reduced::FlagMcMatchChanD,
                  hf_reso_cand_reduced::DebugMcRec,
                  hf_reso_cand_reduced::Origin,
                  hf_reso_cand_reduced::PtGen,
                  hf_reso_cand_reduced::InvMassGen,
                  hf_cand_mc_flag::NTracksDecayed,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(Hf2PrV0McRec, "AOD", "HF2PRV0MCREC",
                  hf_reso_2pr_v0::Prong0Id,
                  hf_reso_2pr_v0::Prong1Id,
                  hf_reso_cand_reduced::FlagMcMatchRec,
                  hf_reso_cand_reduced::FlagMcMatchRecD,
                  hf_reso_cand_reduced::FlagMcMatchChanD,
                  hf_reso_cand_reduced::DebugMcRec,
                  hf_reso_cand_reduced::Origin,
                  hf_reso_cand_reduced::PtGen,
                  hf_reso_cand_reduced::InvMassGen,
                  hf_cand_mc_flag::NTracksDecayed,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(Hf3PrTrkMcRec, "AOD", "HF3PRTRKMCREC",
                  hf_reso_3pr_trk::Prong0Id,
                  hf_reso_3pr_trk::Prong1Id,
                  hf_reso_cand_reduced::FlagMcMatchRec,
                  hf_reso_cand_reduced::FlagMcMatchRecD,
                  hf_reso_cand_reduced::FlagMcMatchChanD,
                  hf_reso_cand_reduced::DebugMcRec,
                  hf_reso_cand_reduced::Origin,
                  hf_reso_cand_reduced::PtGen,
                  hf_reso_cand_reduced::InvMassGen,
                  hf_cand_mc_flag::NTracksDecayed,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfDstarTrkMcRec, "AOD", "HFDSTARTRKMCREC",
                  hf_reso_dstar_trk::Prong0Id,
                  hf_reso_dstar_trk::Prong1Id,
                  hf_reso_cand_reduced::FlagMcMatchRec,
                  hf_reso_cand_reduced::FlagMcMatchRecD,
                  hf_reso_cand_reduced::FlagMcMatchChanD,
                  hf_reso_cand_reduced::DebugMcRec,
                  hf_reso_cand_reduced::Origin,
                  hf_reso_cand_reduced::PtGen,
                  hf_reso_cand_reduced::InvMassGen,
                  hf_cand_mc_flag::NTracksDecayed,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(Hf2PrTrkMcRec, "AOD", "HF2PRTRKMCREC",
                  hf_reso_2pr_trk::Prong0Id,
                  hf_reso_2pr_trk::Prong1Id,
                  hf_reso_cand_reduced::FlagMcMatchRec,
                  hf_reso_cand_reduced::FlagMcMatchRecD,
                  hf_reso_cand_reduced::FlagMcMatchChanD,
                  hf_reso_cand_reduced::DebugMcRec,
                  hf_reso_cand_reduced::Origin,
                  hf_reso_cand_reduced::PtGen,
                  hf_reso_cand_reduced::InvMassGen,
                  hf_cand_mc_flag::NTracksDecayed,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfMcGenRedResos, "AOD", "HFMCGENREDRESO", //! Generation-level MC information on Ds-Resonances candidates for reduced workflow
                  hf_cand_mc_flag::FlagMcMatchGen,
                  hf_reso_cand_reduced::Origin,
                  hf_b0_mc::PtTrack,
                  hf_b0_mc::YTrack,
                  hf_b0_mc::EtaTrack,
                  hf_b0_mc::PtProng0,
                  hf_b0_mc::YProng0,
                  hf_b0_mc::EtaProng0,
                  hf_b0_mc::PtProng1,
                  hf_b0_mc::YProng1,
                  hf_b0_mc::EtaProng1,
                  hf_reso_cand_reduced::InvMassGen,
                  hf_reduced_collision::HfCollisionRejectionMap,
                  o2::soa::Marker<1>);

// Table with same size as HfCandCharmReso
DECLARE_SOA_TABLE(HfMcRecRedResos, "AOD", "HFMCRECREDRESO", //! Reconstruction-level MC information on Ds-Resonances candidates for reduced workflow
                  hf_reso_cand_reduced::FlagMcMatchRec,
                  hf_reso_cand_reduced::FlagMcMatchRecD,
                  hf_reso_cand_reduced::FlagMcMatchChanD,
                  hf_reso_cand_reduced::DebugMcRec,
                  hf_reso_cand_reduced::Origin,
                  hf_reso_cand_reduced::PtGen,
                  hf_reso_cand_reduced::InvMassGen,
                  hf_cand_mc_flag::NTracksDecayed,
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
