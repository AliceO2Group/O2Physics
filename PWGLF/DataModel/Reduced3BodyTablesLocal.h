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

#ifndef PWGLF_DATAMODEL_REDUCED3BODYTABLESLOCAL_H_
#define PWGLF_DATAMODEL_REDUCED3BODYTABLESLOCAL_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "PWGLF/DataModel/Reduced3BodyTables.h"

namespace o2::aod
{

DECLARE_SOA_TABLE(RedLocCollisions, "AOD", "REDCOLLISIONLOC", //! reduced collision table (same structure as the original collision table)
                  o2::soa::Index<>,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::CovXX, collision::CovXY, collision::CovYY, collision::CovXZ, collision::CovYZ, collision::CovZZ,
                  collision::Flags, collision::Chi2, collision::NumContrib,
                  collision::CollisionTime, collision::CollisionTimeRes,
                  bc::RunNumber);

DECLARE_SOA_TABLE(RedLocPVMults, "AOD", "REDLOCPVMULT", //! Multiplicity from the PV contributors, joinable with reducedCollisions
                  mult::MultNTracksPV);

DECLARE_SOA_TABLE(RedLocCentFT0Cs, "AOD", "REDLOCCENTFT0C", //! Reduced Run 3 FT0C centrality table, joinable with reducedCollisions
                  cent::CentFT0C);

DECLARE_SOA_TABLE_FULL(StoredRedLocIUTracks, "RedLocIUTracks", "AOD", "REDLOCIUTRACK", //! On disk version of the track parameters at inner most update (e.g. ITS) as it comes from the tracking
                       o2::soa::Index<>, reducedtracks3body::CollisionId,
                       reducedtracks3body::X, reducedtracks3body::Alpha,
                       reducedtracks3body::Y, reducedtracks3body::Z, reducedtracks3body::Snp, reducedtracks3body::Tgl,
                       reducedtracks3body::Signed1Pt,
                       // cov matrix
                       reducedtracks3body::SigmaY, reducedtracks3body::SigmaZ, reducedtracks3body::SigmaSnp, reducedtracks3body::SigmaTgl, reducedtracks3body::Sigma1Pt,
                       reducedtracks3body::RhoZY, reducedtracks3body::RhoSnpY, reducedtracks3body::RhoSnpZ, reducedtracks3body::RhoTglY, reducedtracks3body::RhoTglZ,
                       reducedtracks3body::RhoTglSnp, reducedtracks3body::Rho1PtY, reducedtracks3body::Rho1PtZ, reducedtracks3body::Rho1PtSnp, reducedtracks3body::Rho1PtTgl,
                       // tracks extra
                       reducedtracks3body::TPCInnerParam, reducedtracks3body::Flags, reducedtracks3body::ITSClusterSizes,
                       reducedtracks3body::TPCNClsFindable, reducedtracks3body::TPCNClsFindableMinusFound, reducedtracks3body::TPCNClsFindableMinusCrossedRows,
                       reducedtracks3body::TRDPattern, reducedtracks3body::TPCChi2NCl, reducedtracks3body::TOFChi2,
                       reducedtracks3body::TPCSignal, reducedtracks3body::TOFExpMom,
                       // TPC PID
                       reducedtracks3body::TPCNSigmaPr, reducedtracks3body::TPCNSigmaPi, reducedtracks3body::TPCNSigmaDe,
                       reducedtracks3body::TOFNSigmaDe,

                       // ----------- dynmaic columns ------------
                       // tracks IU
                       reducedtracks3body::Px<reducedtracks3body::Signed1Pt, reducedtracks3body::Snp, reducedtracks3body::Alpha>,
                       reducedtracks3body::Py<reducedtracks3body::Signed1Pt, reducedtracks3body::Snp, reducedtracks3body::Alpha>,
                       reducedtracks3body::Pz<reducedtracks3body::Signed1Pt, reducedtracks3body::Tgl>,
                       reducedtracks3body::Rapidity<reducedtracks3body::Signed1Pt, reducedtracks3body::Tgl>,
                       reducedtracks3body::Sign<reducedtracks3body::Signed1Pt>,
                       // tracks extra
                       reducedtracks3body::PIDForTracking<reducedtracks3body::Flags>,
                       reducedtracks3body::IsPVContributor<reducedtracks3body::Flags>,
                       reducedtracks3body::HasTPC<reducedtracks3body::DetectorMap>,
                       reducedtracks3body::HasTOF<reducedtracks3body::DetectorMap>,
                       reducedtracks3body::TPCNClsFound<reducedtracks3body::TPCNClsFindable, reducedtracks3body::TPCNClsFindableMinusFound>,
                       reducedtracks3body::TPCNClsCrossedRows<reducedtracks3body::TPCNClsFindable, reducedtracks3body::TPCNClsFindableMinusCrossedRows>,
                       reducedtracks3body::ITSClsSizeInLayer<reducedtracks3body::ITSClusterSizes>,
                       reducedtracks3body::TPCCrossedRowsOverFindableCls<reducedtracks3body::TPCNClsFindable, reducedtracks3body::TPCNClsFindableMinusCrossedRows>);

DECLARE_SOA_EXTENDED_TABLE_USER(RedLocIUTracks, StoredRedLocIUTracks, "REDLOCIUTRACKEXT", //! Track parameters at inner most update (e.g. ITS) as it comes from the tracking
                                reducedtracks3body::Pt,
                                reducedtracks3body::P,
                                reducedtracks3body::Eta,
                                reducedtracks3body::Phi,
                                // cov matrix
                                reducedtracks3body::CYY,
                                reducedtracks3body::CZY,
                                reducedtracks3body::CZZ,
                                reducedtracks3body::CSnpY,
                                reducedtracks3body::CSnpZ,
                                reducedtracks3body::CSnpSnp,
                                reducedtracks3body::CTglY,
                                reducedtracks3body::CTglZ,
                                reducedtracks3body::CTglSnp,
                                reducedtracks3body::CTglTgl,
                                reducedtracks3body::C1PtY,
                                reducedtracks3body::C1PtZ,
                                reducedtracks3body::C1PtSnp,
                                reducedtracks3body::C1PtTgl,
                                reducedtracks3body::C1Pt21Pt2,
                                // tracks extra
                                reducedtracks3body::DetectorMap);

namespace reduceddecay3bodylocal
{
DECLARE_SOA_COLUMN(Phi, phi, float);       //! decay3body radius
DECLARE_SOA_COLUMN(Radius, radius, float); //! decay3body phi
DECLARE_SOA_COLUMN(PosZ, posz, float);     //! decay3body z position
} // namespace reduceddecay3bodylocal

DECLARE_SOA_TABLE(RedLocDecay3Bodys, "AOD", "REDLOCDECAY3BODY", //! reduced 3-body decay table
                  o2::soa::Index<>, reduceddecay3body::CollisionId, reduceddecay3body::Track0Id, reduceddecay3body::Track1Id, reduceddecay3body::Track2Id);

DECLARE_SOA_TABLE(RedLoc3BodyInfo, "AOD", "REDLOC3BODYINFO", //! joinable with RedLocDecay3Bodys
                  reduceddecay3bodylocal::Radius, reduceddecay3bodylocal::Phi, reduceddecay3bodylocal::PosZ);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_REDUCED3BODYTABLESLOCAL_H_
