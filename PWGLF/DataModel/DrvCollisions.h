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
#ifndef PWGLF_DATAMODEL_DRVCOLLISIONS_H_
#define PWGLF_DATAMODEL_DRVCOLLISIONS_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace drvcollision
{
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);    //! Run number
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t); //! Bunch crossing number (globally unique in this run)

DECLARE_SOA_INDEX_COLUMN(BC, bc);        //! Most probably BC to where this collision has occured
DECLARE_SOA_COLUMN(PosX, posX, float);   //! X Vertex position in cm
DECLARE_SOA_COLUMN(PosY, posY, float);   //! Y Vertex position in cm
DECLARE_SOA_COLUMN(PosZ, posZ, float);   //! Z Vertex position in cm
DECLARE_SOA_COLUMN(CovXX, covXX, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovXY, covXY, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovXZ, covXZ, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovYY, covYY, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovYZ, covYZ, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovZZ, covZZ, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(MultNTracksPVeta1, multNTracksPVeta1, int);
DECLARE_SOA_COLUMN(IsInelGt0, isInelGt0, bool);
DECLARE_SOA_COLUMN(IsInelGt1, isInelGt1, bool);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);                    //! Run 2: see CollisionFlagsRun2 | Run 3: see Vertex::Flags
DECLARE_SOA_COLUMN(Chi2, chi2, float);                         //! Chi2 of vertex fit
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint16_t);          //! Number of tracks used for the vertex
DECLARE_SOA_COLUMN(CollisionTime, collisionTime, float);       //! Collision time in ns relative to BC stored in bc()
DECLARE_SOA_COLUMN(CollisionTimeRes, collisionTimeRes, float); //! Resolution of collision time

DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int16_t);
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, int16_t);
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, uint32_t);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(CollisionId, collisionId, int32_t);
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float);
} // namespace drvcollision
DECLARE_SOA_TABLE(DrvTracks, "AOD", "DRVTRACKS",
                  o2::soa::Index<>,
                  drvcollision::CollisionId,
                  drvcollision::Pt,
                  drvcollision::Eta,
                  drvcollision::Phi,
                  drvcollision::Sign,
                  drvcollision::Signed1Pt,
                  drvcollision::Px,
                  drvcollision::Py,
                  drvcollision::Pz,
                  drvcollision::TPCNClsCrossedRows,
                  drvcollision::HasITS,
                  drvcollision::TPCNClsShared,
                  drvcollision::ITSClusterMap,
                  drvcollision::DcaXY,
                  drvcollision::TPCNSigmaPi,
                  drvcollision::TPCNSigmaPr);

DECLARE_SOA_TABLE(DrvCollisions, "AOD", "DRVCOLLISION", //! Time and vertex information of collision
                  o2::soa::Index<>, drvcollision::BCId,
                  drvcollision::PosX, drvcollision::PosY, drvcollision::PosZ, drvcollision::MultNTracksPVeta1,
                  drvcollision::CovXX, drvcollision::CovXY, drvcollision::CovXZ, drvcollision::CovYY, drvcollision::CovYZ, drvcollision::CovZZ,
                  drvcollision::CentFT0C, drvcollision::CentFT0M, drvcollision::IsInelGt0, drvcollision::IsInelGt1, drvcollision::Flags, drvcollision::Chi2, drvcollision::NumContrib,
                  drvcollision::CollisionTime, drvcollision::CollisionTimeRes);

DECLARE_SOA_TABLE(BCandTime, "AOD", "BCANDTIME",
                  o2::soa::Index<>,
                  drvcollision::Timestamp,
                  drvcollision::RunNumber, drvcollision::GlobalBC);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_DRVCOLLISIONS_H_
