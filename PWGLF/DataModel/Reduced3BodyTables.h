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

/// \file Reduced3BodyTables.h
/// \brief Definitions of tables for reduced data of 3body decayed hypertriton analysis
/// \author Carolina Reetz <c.reetz@cern.ch>
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#ifndef PWGLF_DATAMODEL_REDUCED3BODYTABLES_H_
#define PWGLF_DATAMODEL_REDUCED3BODYTABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"

namespace o2::aod
{

DECLARE_SOA_TABLE(RedCollisions, "AOD", "REDCOLLISION", //! reduced collision table (same structure as the original collision table)
                  o2::soa::Index<>,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::CovXX, collision::CovXY, collision::CovYY, collision::CovXZ, collision::CovYZ, collision::CovZZ,
                  collision::Flags, collision::Chi2, collision::NumContrib,
                  collision::CollisionTime, collision::CollisionTimeRes,
                  bc::RunNumber);

DECLARE_SOA_TABLE(RedPVMults, "AOD", "REDPVMULT", //! Multiplicity from the PV contributors, joinable with reducedCollisions
                  mult::MultNTracksPV);

DECLARE_SOA_TABLE(RedCentFT0Cs, "AOD", "REDCENTFT0C", //! Reduced Run 3 FT0C centrality table, joinable with reducedCollisions
                  cent::CentFT0C);

namespace reducedtracks3body
{
// track parameter definition
DECLARE_SOA_INDEX_COLUMN_FULL(Collision, collision, int, RedCollisions, ""); //! Collision index

// track PID definition
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float); //! Nsigma separation with the TPC detector for proton
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float); //! Nsigma separation with the TPC detector for deuteron
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float); //! Nsigma separation with the TPC detector for pion
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float); //! Nsigma separation with the TOF detector for deuteron (recalculated)

} // namespace reducedtracks3body

DECLARE_SOA_TABLE_FULL(StoredRedIUTracks, "RedIUTracks", "AOD", "REDIUTRACK", //! On disk version of the track parameters at inner most update (e.g. ITS) as it comes from the tracking
                       o2::soa::Index<>, reducedtracks3body::CollisionId,
                       track::X, track::Alpha,
                       track::Y, track::Z, track::Snp, track::Tgl,
                       track::Signed1Pt,
                       // cov matrix
                       track::SigmaY, track::SigmaZ, track::SigmaSnp, track::SigmaTgl, track::Sigma1Pt,
                       track::RhoZY, track::RhoSnpY, track::RhoSnpZ, track::RhoTglY, track::RhoTglZ,
                       track::RhoTglSnp, track::Rho1PtY, track::Rho1PtZ, track::Rho1PtSnp, track::Rho1PtTgl,
                       // tracks extra
                       track::TPCInnerParam, track::Flags, track::ITSClusterSizes,
                       track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                       track::TRDPattern, track::TPCChi2NCl, track::TOFChi2,
                       track::TPCSignal, track::TOFExpMom,
                       // TPC PID
                       reducedtracks3body::TPCNSigmaPr, reducedtracks3body::TPCNSigmaPi, reducedtracks3body::TPCNSigmaDe,
                       reducedtracks3body::TOFNSigmaDe,

                       // ----------- dynmaic columns ------------
                       // tracks IU
                       track::Px<track::Signed1Pt, track::Snp, track::Alpha>,
                       track::Py<track::Signed1Pt, track::Snp, track::Alpha>,
                       track::Pz<track::Signed1Pt, track::Tgl>,
                       track::Rapidity<track::Signed1Pt, track::Tgl>,
                       track::Sign<track::Signed1Pt>,
                       // tracks extra
                       track::PIDForTracking<track::Flags>,
                       track::IsPVContributor<track::Flags>,
                       track::HasITS<track::v001::DetectorMap>,
                       track::HasTPC<track::v001::DetectorMap>,
                       track::HasTOF<track::v001::DetectorMap>,
                       track::HasTRD<track::v001::DetectorMap>,
                       track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                       track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                       track::v001::ITSClsSizeInLayer<track::ITSClusterSizes>,
                       track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>);

DECLARE_SOA_EXTENDED_TABLE_USER(RedIUTracks, StoredRedIUTracks, "REDIUTRACKEXT", //! Track parameters at inner most update (e.g. ITS) as it comes from the tracking
                                track::Pt,
                                track::P,
                                track::Eta,
                                track::Phi,
                                // cov matrix
                                track::CYY,
                                track::CZY,
                                track::CZZ,
                                track::CSnpY,
                                track::CSnpZ,
                                track::CSnpSnp,
                                track::CTglY,
                                track::CTglZ,
                                track::CTglSnp,
                                track::CTglTgl,
                                track::C1PtY,
                                track::C1PtZ,
                                track::C1PtSnp,
                                track::C1PtTgl,
                                track::C1Pt21Pt2,
                                // tracks extra
                                track::v001::DetectorMap);

namespace reduceddecay3body
{
DECLARE_SOA_INDEX_COLUMN_FULL(Track0, track0, int, RedIUTracks, "_0");       //! Track 0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Track1, track1, int, RedIUTracks, "_1");       //! Track 1 index
DECLARE_SOA_INDEX_COLUMN_FULL(Track2, track2, int, RedIUTracks, "_2");       //! Track 2 index
DECLARE_SOA_INDEX_COLUMN_FULL(Collision, collision, int, RedCollisions, ""); //! Collision index
DECLARE_SOA_COLUMN(RadiusKF, radiusKF, float);                               //! phi of momentum of mother particle calculated by KF
DECLARE_SOA_COLUMN(PhiKF, phiKF, float);                                     //! SV radius in x-y plane calculated by KF
DECLARE_SOA_COLUMN(PosZKF, poszKF, float);                                   //! z position of SV calculated by KF
DECLARE_SOA_COLUMN(RadiusDCA, radiusDCA, float);                             //! phi of momentum of mother particle calculated by dcaFitter
DECLARE_SOA_COLUMN(PhiDCA, phiDCA, float);                                   //! SV radius in x-y plane calculated by dcaFitter
DECLARE_SOA_COLUMN(PosZDCA, poszDCA, float);                                 //! z position of SV calculated by dcaFitter
DECLARE_SOA_COLUMN(TrackedClSize, trackedClSize, float);                     //! average ITS cluster size (if tracked)
} // namespace reduceddecay3body

DECLARE_SOA_TABLE(RedDecay3Bodys, "AOD", "REDDECAY3BODY", //! reduced 3-body decay table
                  o2::soa::Index<>, reduceddecay3body::CollisionId, reduceddecay3body::Track0Id, reduceddecay3body::Track1Id, reduceddecay3body::Track2Id);

DECLARE_SOA_TABLE(Red3BodyInfo, "AOD", "RED3BODYINFO", //! joinable with RedDecay3Bodys
                  reduceddecay3body::RadiusKF, reduceddecay3body::PhiKF, reduceddecay3body::PosZKF,
                  reduceddecay3body::RadiusDCA, reduceddecay3body::PhiDCA, reduceddecay3body::PosZDCA,
                  reduceddecay3body::TrackedClSize);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_REDUCED3BODYTABLES_H_
