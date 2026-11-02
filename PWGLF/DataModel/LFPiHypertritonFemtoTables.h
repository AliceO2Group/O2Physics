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
/// \file LFPiHypertritonFemtoTables.h
/// \brief Slim tables for pi-hypertriton femto pairs
///

#ifndef PWGLF_DATAMODEL_LFPIHYPERTRITONFEMTOTABLES_H_
#define PWGLF_DATAMODEL_LFPIHYPERTRITONFEMTOTABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace pihypertritonfemto
{
constexpr uint8_t ClosePairAngular = 1 << 0;
constexpr uint8_t ClosePairDistance = 1 << 1;
constexpr uint8_t ClosePairAngularUnavailable = 1 << 2;
constexpr uint8_t ClosePairDistanceUnavailable = 1 << 3;
DECLARE_SOA_COLUMN(IsMixed, isMixed, bool);
// Offline CPR bits, OR-ed across same-sign daughters: 1=angular rejection, 2=distance rejection,
// 4=angular unavailable, 8=distance unavailable. Zero means all applicable checks passed.
// Opposite-sign comparisons require no CPR. All flagged pairs remain in the table.
DECLARE_SOA_COLUMN(IsClosePairRejected, isClosePairRejected, uint8_t);
DECLARE_SOA_COLUMN(MixingDepth, mixingDepth, int);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(TrackOccupancy, trackOccupancy, int);
DECLARE_SOA_COLUMN(Ft0cOccupancy, ft0cOccupancy, float);
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float);
DECLARE_SOA_COLUMN(XPrimVtx, xPrimVtx, float);
DECLARE_SOA_COLUMN(YPrimVtx, yPrimVtx, float);
DECLARE_SOA_COLUMN(ZPrimVtx, zPrimVtx, float);
DECLARE_SOA_COLUMN(HypIsMatter, hypIsMatter, bool);
DECLARE_SOA_COLUMN(HypPtHe3, hypPtHe3, float);
DECLARE_SOA_COLUMN(HypEtaHe3, hypEtaHe3, float);
DECLARE_SOA_COLUMN(HypPhiHe3, hypPhiHe3, float);
DECLARE_SOA_COLUMN(HypPtPi, hypPtPi, float);
DECLARE_SOA_COLUMN(HypEtaPi, hypEtaPi, float);
DECLARE_SOA_COLUMN(HypPhiPi, hypPhiPi, float);
DECLARE_SOA_COLUMN(HypDcaV0Daug, hypDcaV0Daug, float);
DECLARE_SOA_COLUMN(HypDcaHe, hypDcaHe, float);
DECLARE_SOA_COLUMN(HypDcaPi, hypDcaPi, float);
DECLARE_SOA_COLUMN(HypNSigmaHe, hypNSigmaHe, float);
DECLARE_SOA_COLUMN(HypNTPCCrossedRowsHe, hypNTPCCrossedRowsHe, uint8_t);
DECLARE_SOA_COLUMN(HypNTPCCrossedRowsPi, hypNTPCCrossedRowsPi, uint8_t);
DECLARE_SOA_COLUMN(HypTpcMomHe, hypTpcMomHe, float);
DECLARE_SOA_COLUMN(HypTpcMomPi, hypTpcMomPi, float);
DECLARE_SOA_COLUMN(HypTpcSignalHe, hypTpcSignalHe, uint16_t);
DECLARE_SOA_COLUMN(HypTpcSignalPi, hypTpcSignalPi, uint16_t);
DECLARE_SOA_COLUMN(HypItsClusterSizesHe, hypItsClusterSizesHe, uint32_t);
DECLARE_SOA_COLUMN(HypItsClusterSizesPi, hypItsClusterSizesPi, uint32_t);
DECLARE_SOA_COLUMN(HypXDecVtx, hypXDecVtx, float);
DECLARE_SOA_COLUMN(HypYDecVtx, hypYDecVtx, float);
DECLARE_SOA_COLUMN(HypZDecVtx, hypZDecVtx, float);
DECLARE_SOA_COLUMN(HadPt, hadPt, float);
DECLARE_SOA_COLUMN(HadEta, hadEta, float);
DECLARE_SOA_COLUMN(HadPhi, hadPhi, float);
DECLARE_SOA_COLUMN(HadSign, hadSign, int8_t);
DECLARE_SOA_COLUMN(HadDcaXY, hadDcaXY, float);
DECLARE_SOA_COLUMN(HadTpcNClsCrossedRows, hadTpcNClsCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(HadTpcNClsPID, hadTpcNClsPID, uint8_t);
DECLARE_SOA_COLUMN(HadTpcChi2NCl, hadTpcChi2NCl, float);
DECLARE_SOA_COLUMN(HadItsClusterSizes, hadItsClusterSizes, uint32_t);
DECLARE_SOA_COLUMN(HadItsChi2NCl, hadItsChi2NCl, float);
DECLARE_SOA_COLUMN(HadHasTOF, hadHasTOF, bool);
DECLARE_SOA_COLUMN(HadTpcNSigmaPi, hadTpcNSigmaPi, float);
DECLARE_SOA_COLUMN(HadTofNSigmaPi, hadTofNSigmaPi, float);
DECLARE_SOA_COLUMN(HypGenPt, hypGenPt, float);
DECLARE_SOA_COLUMN(HypGenEta, hypGenEta, float);
DECLARE_SOA_COLUMN(HypGenPhi, hypGenPhi, float);
DECLARE_SOA_COLUMN(HypGenXDecVtx, hypGenXDecVtx, float);
DECLARE_SOA_COLUMN(HypGenYDecVtx, hypGenYDecVtx, float);
DECLARE_SOA_COLUMN(HypGenZDecVtx, hypGenZDecVtx, float);
DECLARE_SOA_COLUMN(HypIsReco, hypIsReco, bool);
DECLARE_SOA_COLUMN(HypIsSignal, hypIsSignal, bool);
DECLARE_SOA_COLUMN(HypIsRecoMCCollision, hypIsRecoMCCollision, bool);
DECLARE_SOA_COLUMN(HypIsSurvEvSel, hypIsSurvEvSel, bool);
DECLARE_SOA_COLUMN(HypIsTwoBodyDecay, hypIsTwoBodyDecay, bool);
DECLARE_SOA_COLUMN(HypStatusCode, hypStatusCode, int16_t);
DECLARE_SOA_COLUMN(HeGenPt, heGenPt, float);
DECLARE_SOA_COLUMN(HeIsPhysicalPrimary, heIsPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(DecayPiIsPhysicalPrimary, decayPiIsPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(HadRecoPt, hadRecoPt, float);
DECLARE_SOA_COLUMN(HadRecoEta, hadRecoEta, float);
DECLARE_SOA_COLUMN(HadRecoPhi, hadRecoPhi, float);
DECLARE_SOA_COLUMN(HadGenPt, hadGenPt, float);
DECLARE_SOA_COLUMN(HadGenEta, hadGenEta, float);
DECLARE_SOA_COLUMN(HadGenPhi, hadGenPhi, float);
DECLARE_SOA_COLUMN(HadIsPhysicalPrimary, hadIsPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(HadProcess, hadProcess, int16_t);
DECLARE_SOA_COLUMN(SameMCCollision, sameMCCollision, bool);
DECLARE_SOA_COLUMN(MatchesHypRecoMCCollision, matchesHypRecoMCCollision, bool);
DECLARE_SOA_COLUMN(MatchesPairRecoMCCollision, matchesPairRecoMCCollision, bool);
DECLARE_SOA_COLUMN(IsTruthSelfCorrelation, isTruthSelfCorrelation, bool);
DECLARE_SOA_COLUMN(IsTruePrimaryHadHyperPair, isTruePrimaryHadHyperPair, bool);
DECLARE_SOA_COLUMN(HadIsTruePion, hadIsTruePion, bool); // Truth species, independent of SE/ME and primary status.
} // namespace pihypertritonfemto

DECLARE_SOA_TABLE(PiHypertritonFemtoTable, "AOD", "PIHYPFEMTO",
                  pihypertritonfemto::IsMixed,
                  pihypertritonfemto::IsClosePairRejected,
                  pihypertritonfemto::MixingDepth,
                  pihypertritonfemto::PosZ,
                  pihypertritonfemto::CentFT0C,
                  pihypertritonfemto::TrackOccupancy,
                  pihypertritonfemto::Ft0cOccupancy,
                  pihypertritonfemto::MultFT0C,
                  pihypertritonfemto::XPrimVtx,
                  pihypertritonfemto::YPrimVtx,
                  pihypertritonfemto::ZPrimVtx,
                  pihypertritonfemto::HypIsMatter,
                  pihypertritonfemto::HypPtHe3,
                  pihypertritonfemto::HypEtaHe3,
                  pihypertritonfemto::HypPhiHe3,
                  pihypertritonfemto::HypPtPi,
                  pihypertritonfemto::HypEtaPi,
                  pihypertritonfemto::HypPhiPi,
                  pihypertritonfemto::HypDcaV0Daug,
                  pihypertritonfemto::HypDcaHe,
                  pihypertritonfemto::HypDcaPi,
                  pihypertritonfemto::HypNSigmaHe,
                  pihypertritonfemto::HypNTPCCrossedRowsHe,
                  pihypertritonfemto::HypNTPCCrossedRowsPi,
                  pihypertritonfemto::HypTpcMomHe,
                  pihypertritonfemto::HypTpcMomPi,
                  pihypertritonfemto::HypTpcSignalHe,
                  pihypertritonfemto::HypTpcSignalPi,
                  pihypertritonfemto::HypItsClusterSizesHe,
                  pihypertritonfemto::HypItsClusterSizesPi,
                  pihypertritonfemto::HypXDecVtx,
                  pihypertritonfemto::HypYDecVtx,
                  pihypertritonfemto::HypZDecVtx,
                  pihypertritonfemto::HadPt,
                  pihypertritonfemto::HadEta,
                  pihypertritonfemto::HadPhi,
                  pihypertritonfemto::HadSign,
                  pihypertritonfemto::HadDcaXY,
                  pihypertritonfemto::HadTpcNClsCrossedRows,
                  pihypertritonfemto::HadTpcNClsPID,
                  pihypertritonfemto::HadTpcChi2NCl,
                  pihypertritonfemto::HadItsClusterSizes,
                  pihypertritonfemto::HadItsChi2NCl,
                  pihypertritonfemto::HadHasTOF,
                  pihypertritonfemto::HadTpcNSigmaPi,
                  pihypertritonfemto::HadTofNSigmaPi);

DECLARE_SOA_TABLE(PiHypertritonFemtoTableMC, "AOD", "PIHYPFEMTOMC",
                  pihypertritonfemto::IsMixed,
                  pihypertritonfemto::IsClosePairRejected,
                  pihypertritonfemto::MixingDepth,
                  pihypertritonfemto::PosZ,
                  pihypertritonfemto::CentFT0C,
                  pihypertritonfemto::TrackOccupancy,
                  pihypertritonfemto::Ft0cOccupancy,
                  pihypertritonfemto::MultFT0C,
                  pihypertritonfemto::XPrimVtx,
                  pihypertritonfemto::YPrimVtx,
                  pihypertritonfemto::ZPrimVtx,
                  pihypertritonfemto::HypIsMatter,
                  pihypertritonfemto::HypPtHe3,
                  pihypertritonfemto::HypEtaHe3,
                  pihypertritonfemto::HypPhiHe3,
                  pihypertritonfemto::HypPtPi,
                  pihypertritonfemto::HypEtaPi,
                  pihypertritonfemto::HypPhiPi,
                  pihypertritonfemto::HypDcaV0Daug,
                  pihypertritonfemto::HypDcaHe,
                  pihypertritonfemto::HypDcaPi,
                  pihypertritonfemto::HypNSigmaHe,
                  pihypertritonfemto::HypNTPCCrossedRowsHe,
                  pihypertritonfemto::HypNTPCCrossedRowsPi,
                  pihypertritonfemto::HypTpcMomHe,
                  pihypertritonfemto::HypTpcMomPi,
                  pihypertritonfemto::HypTpcSignalHe,
                  pihypertritonfemto::HypTpcSignalPi,
                  pihypertritonfemto::HypItsClusterSizesHe,
                  pihypertritonfemto::HypItsClusterSizesPi,
                  pihypertritonfemto::HypXDecVtx,
                  pihypertritonfemto::HypYDecVtx,
                  pihypertritonfemto::HypZDecVtx,
                  pihypertritonfemto::HadPt,
                  pihypertritonfemto::HadEta,
                  pihypertritonfemto::HadPhi,
                  pihypertritonfemto::HadSign,
                  pihypertritonfemto::HadDcaXY,
                  pihypertritonfemto::HadTpcNClsCrossedRows,
                  pihypertritonfemto::HadTpcNClsPID,
                  pihypertritonfemto::HadTpcChi2NCl,
                  pihypertritonfemto::HadItsClusterSizes,
                  pihypertritonfemto::HadItsChi2NCl,
                  pihypertritonfemto::HadHasTOF,
                  pihypertritonfemto::HadTpcNSigmaPi,
                  pihypertritonfemto::HadTofNSigmaPi,
                  pihypertritonfemto::HypGenPt,
                  pihypertritonfemto::HypGenEta,
                  pihypertritonfemto::HypGenPhi,
                  pihypertritonfemto::HypGenXDecVtx,
                  pihypertritonfemto::HypGenYDecVtx,
                  pihypertritonfemto::HypGenZDecVtx,
                  pihypertritonfemto::HypIsReco,
                  pihypertritonfemto::HypIsSignal,
                  pihypertritonfemto::HypIsRecoMCCollision,
                  pihypertritonfemto::HypIsSurvEvSel,
                  pihypertritonfemto::HypIsTwoBodyDecay,
                  pihypertritonfemto::HypStatusCode,
                  pihypertritonfemto::HeGenPt,
                  pihypertritonfemto::HeIsPhysicalPrimary,
                  pihypertritonfemto::DecayPiIsPhysicalPrimary,
                  pihypertritonfemto::HadRecoPt,
                  pihypertritonfemto::HadRecoEta,
                  pihypertritonfemto::HadRecoPhi,
                  pihypertritonfemto::HadGenPt,
                  pihypertritonfemto::HadGenEta,
                  pihypertritonfemto::HadGenPhi,
                  pihypertritonfemto::HadIsPhysicalPrimary,
                  pihypertritonfemto::HadProcess,
                  pihypertritonfemto::SameMCCollision,
                  pihypertritonfemto::MatchesHypRecoMCCollision,
                  pihypertritonfemto::MatchesPairRecoMCCollision,
                  pihypertritonfemto::IsTruthSelfCorrelation,
                  pihypertritonfemto::IsTruePrimaryHadHyperPair,
                  pihypertritonfemto::HadIsTruePion);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFPIHYPERTRITONFEMTOTABLES_H_
