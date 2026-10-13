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
///
/// \file LambdaSpinPolarizationDerived.h
///
/// \brief Derived tables of Lambda candidates with loose selections, input of lambdaSpinPolarization
/// \author Subhadeep Roy <subhadeep.roy@cern.ch>

#ifndef PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LAMBDASPINPOLARIZATIONDERIVED_H_
#define PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LAMBDASPINPOLARIZATIONDERIVED_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace lambdaloosecollision
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Mult, mult, float);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, uint64_t);
} // namespace lambdaloosecollision

DECLARE_SOA_TABLE(LambdaLooseCollisions, "AOD", "LAMBDALOOSECOLS",
                  o2::soa::Index<>,
                  lambdaloosecollision::Cent,
                  lambdaloosecollision::Mult,
                  lambdaloosecollision::PosX,
                  lambdaloosecollision::PosY,
                  lambdaloosecollision::PosZ,
                  lambdaloosecollision::TimeStamp);
using LambdaLooseCollision = LambdaLooseCollisions::iterator;

namespace lambdaloosetrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaLooseCollision, lambdaLooseCollision);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(MassOther, massOther, float);
DECLARE_SOA_COLUMN(PrPx, prPx, float);
DECLARE_SOA_COLUMN(PrPy, prPy, float);
DECLARE_SOA_COLUMN(PrPz, prPz, float);
DECLARE_SOA_COLUMN(PosTrackId, posTrackId, int32_t);
DECLARE_SOA_COLUMN(NegTrackId, negTrackId, int32_t);
DECLARE_SOA_COLUMN(V0Type, v0Type, int8_t);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DcaDau, dcaDau, float);
DECLARE_SOA_COLUMN(DcaV0ToPV, dcaV0ToPV, float);
DECLARE_SOA_COLUMN(CTau, cTau, float);
DECLARE_SOA_COLUMN(MK0Short, mK0Short, float);
DECLARE_SOA_COLUMN(PosNSigmaPr, posNSigmaPr, float);
DECLARE_SOA_COLUMN(PosNSigmaPi, posNSigmaPi, float);
DECLARE_SOA_COLUMN(NegNSigmaPr, negNSigmaPr, float);
DECLARE_SOA_COLUMN(NegNSigmaPi, negNSigmaPi, float);
DECLARE_SOA_COLUMN(PosTpcCrossedRows, posTpcCrossedRows, int16_t);
DECLARE_SOA_COLUMN(NegTpcCrossedRows, negTpcCrossedRows, int16_t);
} // namespace lambdaloosetrack

DECLARE_SOA_TABLE(LambdaLooseTracks, "AOD", "LAMBDALOOSETRKS",
                  o2::soa::Index<>,
                  lambdaloosetrack::LambdaLooseCollisionId,
                  lambdaloosetrack::Px,
                  lambdaloosetrack::Py,
                  lambdaloosetrack::Pz,
                  lambdaloosetrack::Pt,
                  lambdaloosetrack::Eta,
                  lambdaloosetrack::Phi,
                  lambdaloosetrack::Rap,
                  lambdaloosetrack::Mass,
                  lambdaloosetrack::MassOther,
                  lambdaloosetrack::PrPx,
                  lambdaloosetrack::PrPy,
                  lambdaloosetrack::PrPz,
                  lambdaloosetrack::PosTrackId,
                  lambdaloosetrack::NegTrackId,
                  lambdaloosetrack::V0Type,
                  lambdaloosetrack::CosPA,
                  lambdaloosetrack::DcaDau,
                  lambdaloosetrack::DcaV0ToPV,
                  lambdaloosetrack::CTau,
                  lambdaloosetrack::MK0Short,
                  lambdaloosetrack::PosNSigmaPr,
                  lambdaloosetrack::PosNSigmaPi,
                  lambdaloosetrack::NegNSigmaPr,
                  lambdaloosetrack::NegNSigmaPi,
                  lambdaloosetrack::PosTpcCrossedRows,
                  lambdaloosetrack::NegTpcCrossedRows);
using LambdaLooseTrack = LambdaLooseTracks::iterator;
} // namespace o2::aod

#endif // PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LAMBDASPINPOLARIZATIONDERIVED_H_
