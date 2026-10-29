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
/// \author Sourav Kundu <sourav.kundu@cern.ch>

#ifndef PWGLF_DATAMODEL_REDUCEDDOUBLEPHITABLES_H_
#define PWGLF_DATAMODEL_REDUCEDDOUBLEPHITABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace redphievent
{
DECLARE_SOA_COLUMN(NumPos, numPos, int); //! Number of positive Kaon
DECLARE_SOA_COLUMN(NumNeg, numNeg, int); //! Number of negative Kaon
DECLARE_SOA_COLUMN(Centrality, centrality, float); //! Number of negative Kaon
} // namespace redphievent
DECLARE_SOA_TABLE(RedPhiEvents, "AOD", "REDPHIEVENT",
                  o2::soa::Index<>,
                  bc::GlobalBC,
                  bc::RunNumber,
                  timestamp::Timestamp,
                  collision::PosZ,
                  collision::NumContrib,
                  redphievent::NumPos,
                  redphievent::NumNeg,
                  redphievent::Centrality);
using RedPhiEvent = RedPhiEvents::iterator;

namespace phitrack
{
DECLARE_SOA_INDEX_COLUMN(RedPhiEvent, redPhiEvent);
DECLARE_SOA_COLUMN(PhiPx, phiPx, float);             //! Phi Px
DECLARE_SOA_COLUMN(PhiPy, phiPy, float);             //! Phi Py
DECLARE_SOA_COLUMN(PhiPz, phiPz, float);             //! Phi Pz
DECLARE_SOA_COLUMN(Phid1Px, phid1Px, float);         //! Phi d1 Px
DECLARE_SOA_COLUMN(Phid1Py, phid1Py, float);         //! Phi d1 Py
DECLARE_SOA_COLUMN(Phid1Pz, phid1Pz, float);         //! Phi d1 Pz
DECLARE_SOA_COLUMN(Phid2Px, phid2Px, float);         //! Phi d2 Px
DECLARE_SOA_COLUMN(Phid2Py, phid2Py, float);         //! Phi d2 Py
DECLARE_SOA_COLUMN(Phid2Pz, phid2Pz, float);         //! Phi d2 Pz
DECLARE_SOA_COLUMN(PhiMass, phiMass, float);         //! Phi Mass
DECLARE_SOA_COLUMN(Phid1Index, phid1Index, int64_t); //! Phi d1 index
DECLARE_SOA_COLUMN(Phid2Index, phid2Index, int64_t); //! Phi d2 index
DECLARE_SOA_COLUMN(Phid1Charge, phid1Charge, float); //! Phi d1 charge
DECLARE_SOA_COLUMN(Phid2Charge, phid2Charge, float); //! Phi d1 charge
DECLARE_SOA_COLUMN(Phid1TPC, phid1TPC, float);       //! TPC nsigma d1
DECLARE_SOA_COLUMN(Phid2TPC, phid2TPC, float);       //! TPC nsigma d2
DECLARE_SOA_COLUMN(Phid1TOFHit, phid1TOFHit, int);   //! TOF hit d1
DECLARE_SOA_COLUMN(Phid2TOFHit, phid2TOFHit, int);   //! TOF hit d2
DECLARE_SOA_COLUMN(Phid1TOF, phid1TOF, float);       //! TOF nsigma d1
DECLARE_SOA_COLUMN(Phid2TOF, phid2TOF, float);       //! TOF nsigma d2
} // namespace phitrack
DECLARE_SOA_TABLE(PhiTracks, "AOD", "PHITRACK",
                  o2::soa::Index<>,
                  phitrack::RedPhiEventId,
                  phitrack::PhiPx,
                  phitrack::PhiPy,
                  phitrack::PhiPz,
                  phitrack::Phid1Px,
                  phitrack::Phid1Py,
                  phitrack::Phid1Pz,
                  phitrack::Phid2Px,
                  phitrack::Phid2Py,
                  phitrack::Phid2Pz,
                  phitrack::PhiMass,
                  phitrack::Phid1Index,
                  phitrack::Phid2Index,
                  phitrack::Phid1Charge,
                  phitrack::Phid2Charge,
                  phitrack::Phid1TPC,
                  phitrack::Phid2TPC,
                  phitrack::Phid1TOFHit,
                  phitrack::Phid2TOFHit,
                  phitrack::Phid1TOF,
                  phitrack::Phid2TOF);

using PhiTrack = PhiTracks::iterator;

namespace phiphipair
{
DECLARE_SOA_INDEX_COLUMN(RedPhiEvent, redPhiEvent);

DECLARE_SOA_COLUMN(Phi1Row, phi1Row, int64_t);
DECLARE_SOA_COLUMN(Phi2Row, phi2Row, int64_t);

DECLARE_SOA_COLUMN(PairMass, pairMass, float);
DECLARE_SOA_COLUMN(PairPx, pairPx, float);
DECLARE_SOA_COLUMN(PairPy, pairPy, float);
DECLARE_SOA_COLUMN(PairPz, pairPz, float);

DECLARE_SOA_COLUMN(Phi1Mass, phi1Mass, float);
DECLARE_SOA_COLUMN(Phi1Px, phi1Px, float);
DECLARE_SOA_COLUMN(Phi1Py, phi1Py, float);
DECLARE_SOA_COLUMN(Phi1Pz, phi1Pz, float);

DECLARE_SOA_COLUMN(Phi2Mass, phi2Mass, float);
DECLARE_SOA_COLUMN(Phi2Px, phi2Px, float);
DECLARE_SOA_COLUMN(Phi2Py, phi2Py, float);
DECLARE_SOA_COLUMN(Phi2Pz, phi2Pz, float);

// Kaon 1: positive daughter of phi 1
DECLARE_SOA_COLUMN(K1Index, k1Index, int64_t);
DECLARE_SOA_COLUMN(K1Charge, k1Charge, int8_t);
DECLARE_SOA_COLUMN(K1DcaXY, k1DcaXY, float);
DECLARE_SOA_COLUMN(K1DcaZ, k1DcaZ, float);
DECLARE_SOA_COLUMN(K1DcaXYSig, k1DcaXYSig, float);
DECLARE_SOA_COLUMN(K1DcaZSig, k1DcaZSig, float);
DECLARE_SOA_COLUMN(K1DcaChi2, k1DcaChi2, float);
DECLARE_SOA_COLUMN(K1Px, k1Px, float);
DECLARE_SOA_COLUMN(K1Py, k1Py, float);
DECLARE_SOA_COLUMN(K1Pz, k1Pz, float);
DECLARE_SOA_COLUMN(K1TOFHit, k1TOFHit, int8_t);
DECLARE_SOA_COLUMN(K1TPC, k1TPC, float);
DECLARE_SOA_COLUMN(K1TOF, k1TOF, float);
DECLARE_SOA_COLUMN(K1ITS, k1ITS, float);

// Kaon 2: negative daughter of phi 1
DECLARE_SOA_COLUMN(K2Index, k2Index, int64_t);
DECLARE_SOA_COLUMN(K2Charge, k2Charge, int8_t);
DECLARE_SOA_COLUMN(K2DcaXY, k2DcaXY, float);
DECLARE_SOA_COLUMN(K2DcaZ, k2DcaZ, float);
DECLARE_SOA_COLUMN(K2DcaXYSig, k2DcaXYSig, float);
DECLARE_SOA_COLUMN(K2DcaZSig, k2DcaZSig, float);
DECLARE_SOA_COLUMN(K2DcaChi2, k2DcaChi2, float);
DECLARE_SOA_COLUMN(K2Px, k2Px, float);
DECLARE_SOA_COLUMN(K2Py, k2Py, float);
DECLARE_SOA_COLUMN(K2Pz, k2Pz, float);
DECLARE_SOA_COLUMN(K2TOFHit, k2TOFHit, int8_t);
DECLARE_SOA_COLUMN(K2TPC, k2TPC, float);
DECLARE_SOA_COLUMN(K2TOF, k2TOF, float);
DECLARE_SOA_COLUMN(K2ITS, k2ITS, float);

// Kaon 3: positive daughter of phi 2
DECLARE_SOA_COLUMN(K3Index, k3Index, int64_t);
DECLARE_SOA_COLUMN(K3Charge, k3Charge, int8_t);
DECLARE_SOA_COLUMN(K3DcaXY, k3DcaXY, float);
DECLARE_SOA_COLUMN(K3DcaZ, k3DcaZ, float);
DECLARE_SOA_COLUMN(K3DcaXYSig, k3DcaXYSig, float);
DECLARE_SOA_COLUMN(K3DcaZSig, k3DcaZSig, float);
DECLARE_SOA_COLUMN(K3DcaChi2, k3DcaChi2, float);
DECLARE_SOA_COLUMN(K3Px, k3Px, float);
DECLARE_SOA_COLUMN(K3Py, k3Py, float);
DECLARE_SOA_COLUMN(K3Pz, k3Pz, float);
DECLARE_SOA_COLUMN(K3TOFHit, k3TOFHit, int8_t);
DECLARE_SOA_COLUMN(K3TPC, k3TPC, float);
DECLARE_SOA_COLUMN(K3TOF, k3TOF, float);
DECLARE_SOA_COLUMN(K3ITS, k3ITS, float);

// Kaon 4: negative daughter of phi 2
DECLARE_SOA_COLUMN(K4Index, k4Index, int64_t);
DECLARE_SOA_COLUMN(K4Charge, k4Charge, int8_t);
DECLARE_SOA_COLUMN(K4DcaXY, k4DcaXY, float);
DECLARE_SOA_COLUMN(K4DcaZ, k4DcaZ, float);
DECLARE_SOA_COLUMN(K4DcaXYSig, k4DcaXYSig, float);
DECLARE_SOA_COLUMN(K4DcaZSig, k4DcaZSig, float);
DECLARE_SOA_COLUMN(K4DcaChi2, k4DcaChi2, float);
DECLARE_SOA_COLUMN(K4Px, k4Px, float);
DECLARE_SOA_COLUMN(K4Py, k4Py, float);
DECLARE_SOA_COLUMN(K4Pz, k4Pz, float);
DECLARE_SOA_COLUMN(K4TOFHit, k4TOFHit, int8_t);
DECLARE_SOA_COLUMN(K4TPC, k4TPC, float);
DECLARE_SOA_COLUMN(K4TOF, k4TOF, float);
DECLARE_SOA_COLUMN(K4ITS, k4ITS, float);

// Four-track vertex fit
DECLARE_SOA_COLUMN(FitStatus, fitStatus, int8_t);
DECLARE_SOA_COLUMN(FitChi2, fitChi2, float);
DECLARE_SOA_COLUMN(FitNdf, fitNdf, int8_t);
DECLARE_SOA_COLUMN(FitChi2Ndf, fitChi2Ndf, float);

// Primary vertex and fitted common vertex
DECLARE_SOA_COLUMN(PvX, pvX, float);
DECLARE_SOA_COLUMN(PvY, pvY, float);
DECLARE_SOA_COLUMN(PvZ, pvZ, float);

DECLARE_SOA_COLUMN(VtxX, vtxX, float);
DECLARE_SOA_COLUMN(VtxY, vtxY, float);
DECLARE_SOA_COLUMN(VtxZ, vtxZ, float);

DECLARE_SOA_COLUMN(DeltaVtxX, deltaVtxX, float);
DECLARE_SOA_COLUMN(DeltaVtxY, deltaVtxY, float);
DECLARE_SOA_COLUMN(DeltaVtxZ, deltaVtxZ, float);

DECLARE_SOA_COLUMN(VertexLxy, vertexLxy, float);
DECLARE_SOA_COLUMN(VertexL3D, vertexL3D, float);
DECLARE_SOA_COLUMN(VertexLxyErr, vertexLxyErr, float);
DECLARE_SOA_COLUMN(VertexL3DErr, vertexL3DErr, float);
DECLARE_SOA_COLUMN(VertexLxySig, vertexLxySig, float);
DECLARE_SOA_COLUMN(VertexL3DSig, vertexL3DSig, float);

// Resolution-weighted DCA quantities
DECLARE_SOA_COLUMN(NValidDca, nValidDca, int8_t);
DECLARE_SOA_COLUMN(SumDcaXYSig2, sumDcaXYSig2, float);
DECLARE_SOA_COLUMN(SumDcaZSig2, sumDcaZSig2, float);
DECLARE_SOA_COLUMN(SumDcaChi2, sumDcaChi2, float);
DECLARE_SOA_COLUMN(RmsDcaSig, rmsDcaSig, float);
DECLARE_SOA_COLUMN(MaxDcaChi2, maxDcaChi2, float);
DECLARE_SOA_COLUMN(MaxAbsDcaXYSig, maxAbsDcaXYSig, float);
DECLARE_SOA_COLUMN(MaxAbsDcaZSig, maxAbsDcaZSig, float);
} // namespace phiphipair

// clang-format off
DECLARE_SOA_TABLE(PhiPhiPairs, "AOD", "PHIPHIPAIR",
                  o2::soa::Index<>,
                  phiphipair::RedPhiEventId,
                  phiphipair::Phi1Row,
                  phiphipair::Phi2Row,

                  phiphipair::PairMass,
                  phiphipair::PairPx,
                  phiphipair::PairPy,
                  phiphipair::PairPz,

                  phiphipair::Phi1Mass,
                  phiphipair::Phi1Px,
                  phiphipair::Phi1Py,
                  phiphipair::Phi1Pz,

                  phiphipair::Phi2Mass,
                  phiphipair::Phi2Px,
                  phiphipair::Phi2Py,
                  phiphipair::Phi2Pz,

                  phiphipair::K1Index,
                  phiphipair::K1Charge,
                  phiphipair::K1DcaXY,
                  phiphipair::K1DcaZ,
                  phiphipair::K1DcaXYSig,
                  phiphipair::K1DcaZSig,
                  phiphipair::K1DcaChi2,
                  phiphipair::K1Px,
                  phiphipair::K1Py,
                  phiphipair::K1Pz,
                  phiphipair::K1TOFHit,
                  phiphipair::K1TPC,
                  phiphipair::K1TOF,
                  phiphipair::K1ITS,

                  phiphipair::K2Index,
                  phiphipair::K2Charge,
                  phiphipair::K2DcaXY,
                  phiphipair::K2DcaZ,
                  phiphipair::K2DcaXYSig,
                  phiphipair::K2DcaZSig,
                  phiphipair::K2DcaChi2,
                  phiphipair::K2Px,
                  phiphipair::K2Py,
                  phiphipair::K2Pz,
                  phiphipair::K2TOFHit,
                  phiphipair::K2TPC,
                  phiphipair::K2TOF,
                  phiphipair::K2ITS,

                  phiphipair::K3Index,
                  phiphipair::K3Charge,
                  phiphipair::K3DcaXY,
                  phiphipair::K3DcaZ,
                  phiphipair::K3DcaXYSig,
                  phiphipair::K3DcaZSig,
                  phiphipair::K3DcaChi2,
                  phiphipair::K3Px,
                  phiphipair::K3Py,
                  phiphipair::K3Pz,
                  phiphipair::K3TOFHit,
                  phiphipair::K3TPC,
                  phiphipair::K3TOF,
                  phiphipair::K3ITS,

                  phiphipair::K4Index,
                  phiphipair::K4Charge,
                  phiphipair::K4DcaXY,
                  phiphipair::K4DcaZ,
                  phiphipair::K4DcaXYSig,
                  phiphipair::K4DcaZSig,
                  phiphipair::K4DcaChi2,
                  phiphipair::K4Px,
                  phiphipair::K4Py,
                  phiphipair::K4Pz,
                  phiphipair::K4TOFHit,
                  phiphipair::K4TPC,
                  phiphipair::K4TOF,
                  phiphipair::K4ITS,

                  phiphipair::FitStatus,
                  phiphipair::FitChi2,
                  phiphipair::FitNdf,
                  phiphipair::FitChi2Ndf,

                  phiphipair::PvX,
                  phiphipair::PvY,
                  phiphipair::PvZ,
                  phiphipair::VtxX,
                  phiphipair::VtxY,
                  phiphipair::VtxZ,
                  phiphipair::DeltaVtxX,
                  phiphipair::DeltaVtxY,
                  phiphipair::DeltaVtxZ,

                  phiphipair::VertexLxy,
                  phiphipair::VertexL3D,
                  phiphipair::VertexLxyErr,
                  phiphipair::VertexL3DErr,
                  phiphipair::VertexLxySig,
                  phiphipair::VertexL3DSig,

                  phiphipair::NValidDca,
                  phiphipair::SumDcaXYSig2,
                  phiphipair::SumDcaZSig2,
                  phiphipair::SumDcaChi2,
                  phiphipair::RmsDcaSig,
                  phiphipair::MaxDcaChi2,
                  phiphipair::MaxAbsDcaXYSig,
                  phiphipair::MaxAbsDcaZSig);
// clang-format on

using PhiPhiPair = PhiPhiPairs::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_REDUCEDDOUBLEPHITABLES_H_
