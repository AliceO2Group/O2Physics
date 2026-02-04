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
/// \file LFAntinCexTables.h
/// \brief Slim tables for nucleiAntineutronCex
/// \author Fabiola Lugo
///

#ifndef PWGLF_DATAMODEL_LFANTINCEXTABLES_H_
#define PWGLF_DATAMODEL_LFANTINCEXTABLES_H_

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace antin_cex
{
// Metadata
DECLARE_SOA_COLUMN(IsCex, isCex, bool);            // 1=CEX (from antin), 0=BG
DECLARE_SOA_COLUMN(MotherPdg, motherPdg, int32_t); // mother PDG
DECLARE_SOA_COLUMN(ColId, colId, int32_t);         // mcCollisionId
DECLARE_SOA_COLUMN(PId, pId, int32_t);             // proton MC id
DECLARE_SOA_COLUMN(AntipId, antipId, int32_t);     // antiproton MC id

// MC (pair)
DECLARE_SOA_COLUMN(McPairP, mcPairP, float);
DECLARE_SOA_COLUMN(McPairPt, mcPairPt, float);
DECLARE_SOA_COLUMN(McPairPz, mcPairPz, float);
DECLARE_SOA_COLUMN(McDplane, mcDplane, float);
DECLARE_SOA_COLUMN(McAngleDeg, mcAngleDeg, float);
DECLARE_SOA_COLUMN(McVtxX, mcVtxX, float);
DECLARE_SOA_COLUMN(McVtxY, mcVtxY, float);
DECLARE_SOA_COLUMN(McVtxZ, mcVtxZ, float);

// Tracks (pair, fitter)
DECLARE_SOA_COLUMN(TrkPairP, trkPairP, float);
DECLARE_SOA_COLUMN(TrkPairPt, trkPairPt, float);
DECLARE_SOA_COLUMN(TrkPairPz, trkPairPz, float);
DECLARE_SOA_COLUMN(TrkAngleDeg, trkAngleDeg, float);
DECLARE_SOA_COLUMN(TrkVtxfitDcaPair, trkVtxfitDcaPair, float);
DECLARE_SOA_COLUMN(TrkVtxfitR, trkVtxfitR, float);
DECLARE_SOA_COLUMN(TrkVtxfitDistToPv, trkVtxfitDistToPv, float);
DECLARE_SOA_COLUMN(TrkVtxfitSecVtxX, trkVtxfitSecVtxX, float);
DECLARE_SOA_COLUMN(TrkVtxfitSecVtxY, trkVtxfitSecVtxY, float);
DECLARE_SOA_COLUMN(TrkVtxfitSecVtxZ, trkVtxfitSecVtxZ, float);

// Fit quality (fit − MC)
DECLARE_SOA_COLUMN(VtxfitChi2, vtxfitChi2, float);
DECLARE_SOA_COLUMN(VtxfitStatus, vtxfitStatus, int32_t);
DECLARE_SOA_COLUMN(NCand, nCand, int32_t);
DECLARE_SOA_COLUMN(VtxfitDX, vtxfitDX, float);
DECLARE_SOA_COLUMN(VtxfitDY, vtxfitDY, float);
DECLARE_SOA_COLUMN(VtxfitDZ, vtxfitDZ, float);
DECLARE_SOA_COLUMN(VtxfitD3D, vtxfitD3D, float);

// Proton track
DECLARE_SOA_COLUMN(PTrkP, pTrkP, float);
DECLARE_SOA_COLUMN(PTrkPx, pTrkPx, float);
DECLARE_SOA_COLUMN(PTrkPy, pTrkPy, float);
DECLARE_SOA_COLUMN(PTrkPz, pTrkPz, float);
DECLARE_SOA_COLUMN(PTrkEta, pTrkEta, float);
DECLARE_SOA_COLUMN(PTrkTpcSignal, pTrkTpcSignal, float);
DECLARE_SOA_COLUMN(PTrkNClsIts, pTrkNClsIts, int16_t);

// Antiproton track
DECLARE_SOA_COLUMN(AntipTrkP, antipTrkP, float);
DECLARE_SOA_COLUMN(AntipTrkPx, antipTrkPx, float);
DECLARE_SOA_COLUMN(AntipTrkPy, antipTrkPy, float);
DECLARE_SOA_COLUMN(AntipTrkPz, antipTrkPz, float);
DECLARE_SOA_COLUMN(AntipTrkEta, antipTrkEta, float);
DECLARE_SOA_COLUMN(AntipTrkTpcSignal, antipTrkTpcSignal, float);
DECLARE_SOA_COLUMN(AntipTrkNClsIts, antipTrkNClsIts, int16_t);

// Cuts Mask
DECLARE_SOA_COLUMN(SelMask, selMask, uint32_t);

DECLARE_SOA_COLUMN(PairPointingAngleDeg, pairPointingAngleDeg, float);
DECLARE_SOA_COLUMN(PairPBalance, pairPBalance, float);
DECLARE_SOA_COLUMN(PairPtBalance, pairPtBalance, float);
DECLARE_SOA_COLUMN(PairQ, pairQ, float);

DECLARE_SOA_COLUMN(DPairP, dPairP, float);
DECLARE_SOA_COLUMN(DPairPt, dPairPt, float);
DECLARE_SOA_COLUMN(DPairPz, dPairPz, float);
DECLARE_SOA_COLUMN(DOpenAngle, dOpenAngle, float);

DECLARE_SOA_COLUMN(SVNearestLayerId, svNearestLayerId, int16_t);
DECLARE_SOA_COLUMN(SVDeltaRToLayer, svDeltaRToLayer, float);

DECLARE_SOA_COLUMN(PTrkItsHitMap, pTrkItsHitMap, uint16_t);
DECLARE_SOA_COLUMN(APTrkItsHitMap, apTrkItsHitMap, uint16_t);
DECLARE_SOA_COLUMN(PLayersOk, pLayersOk, int8_t);
DECLARE_SOA_COLUMN(APLayersOk, apLayersOk, int8_t);

DECLARE_SOA_COLUMN(PVtxZ, pVtxZ, float);

// Proton ITS PID
DECLARE_SOA_COLUMN(PTrkItsNSigmaPr, pTrkItsNSigmaPr, float);
DECLARE_SOA_COLUMN(PTrkItsPidValid, pTrkItsPidValid, int8_t);
DECLARE_SOA_COLUMN(PTrkTgl, pTrkTgl, float);

// Antiproton ITS PID
DECLARE_SOA_COLUMN(AntipTrkItsNSigmaPr, antipTrkItsNSigmaPr, float);
DECLARE_SOA_COLUMN(AntipTrkItsPidValid, antipTrkItsPidValid, int8_t);
DECLARE_SOA_COLUMN(AntipTrkTgl, antipTrkTgl, float);
} // namespace antin_cex

// Table
DECLARE_SOA_TABLE(AntinCexPairs, "AOD", "ANTINCEX",
                  antin_cex::IsCex,
                  antin_cex::MotherPdg, antin_cex::ColId, antin_cex::PId, antin_cex::AntipId,
                  antin_cex::McPairP, antin_cex::McPairPt, antin_cex::McPairPz,
                  antin_cex::McDplane, antin_cex::McAngleDeg, antin_cex::McVtxX, antin_cex::McVtxY, antin_cex::McVtxZ,
                  antin_cex::TrkPairP, antin_cex::TrkPairPt, antin_cex::TrkPairPz, antin_cex::TrkAngleDeg,
                  antin_cex::TrkVtxfitDcaPair, antin_cex::TrkVtxfitR, antin_cex::TrkVtxfitDistToPv,
                  antin_cex::TrkVtxfitSecVtxX, antin_cex::TrkVtxfitSecVtxY, antin_cex::TrkVtxfitSecVtxZ,
                  antin_cex::VtxfitChi2, antin_cex::VtxfitStatus, antin_cex::NCand,
                  antin_cex::VtxfitDX, antin_cex::VtxfitDY, antin_cex::VtxfitDZ, antin_cex::VtxfitD3D,
                  antin_cex::PTrkP, antin_cex::PTrkPx, antin_cex::PTrkPy, antin_cex::PTrkPz, antin_cex::PTrkEta, antin_cex::PTrkTpcSignal, antin_cex::PTrkNClsIts,
                  antin_cex::AntipTrkP, antin_cex::AntipTrkPx, antin_cex::AntipTrkPy, antin_cex::AntipTrkPz, antin_cex::AntipTrkEta, antin_cex::AntipTrkTpcSignal, antin_cex::AntipTrkNClsIts,
                  antin_cex::SelMask,
                  antin_cex::PairPointingAngleDeg, antin_cex::PairPBalance, antin_cex::PairPtBalance, antin_cex::PairQ,
                  antin_cex::DPairP, antin_cex::DPairPt, antin_cex::DPairPz, antin_cex::DOpenAngle,
                  antin_cex::SVNearestLayerId, antin_cex::SVDeltaRToLayer,
                  antin_cex::PTrkItsHitMap, antin_cex::APTrkItsHitMap, antin_cex::PLayersOk, antin_cex::APLayersOk,
                  antin_cex::PVtxZ,
                  antin_cex::PTrkItsNSigmaPr, antin_cex::PTrkItsPidValid, antin_cex::PTrkTgl,
                  antin_cex::AntipTrkItsNSigmaPr, antin_cex::AntipTrkItsPidValid, antin_cex::AntipTrkTgl);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFANTINCEXTABLES_H_
