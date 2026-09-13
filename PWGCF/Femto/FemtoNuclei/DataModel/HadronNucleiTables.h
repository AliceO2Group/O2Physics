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
/// \file HadronNucleiTables.h
/// \brief Slim tables for piNuclei
/// \author CMY
/// \date 2025-04-10

#ifndef PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_HADRONNUCLEITABLES_H_
#define PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_HADRONNUCLEITABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace hadron_nuclei_tables
{

DECLARE_SOA_COLUMN(PtNu, ptNu, float);
DECLARE_SOA_COLUMN(EtaNu, etaNu, float);
DECLARE_SOA_COLUMN(PhiNu, phiNu, float);
DECLARE_SOA_COLUMN(PtHyp, ptHyp, float);
DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);
DECLARE_SOA_COLUMN(EtaHyp, etaHyp, float);
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);
DECLARE_SOA_COLUMN(PhiHyp, phiHyp, float);
DECLARE_SOA_COLUMN(PtHad, ptHad, float);
DECLARE_SOA_COLUMN(EtaHad, etaHad, float);
DECLARE_SOA_COLUMN(PhiHad, phiHad, float);

DECLARE_SOA_COLUMN(DcaxyNu, dcaxyNu, float);
DECLARE_SOA_COLUMN(DcazNu, dcazNu, float);
DECLARE_SOA_COLUMN(DcaxyHad, dcaxyHad, float);
DECLARE_SOA_COLUMN(DcazHad, dcazHad, float);
DECLARE_SOA_COLUMN(DcaPair, dcaPair, float);

DECLARE_SOA_COLUMN(SignalTPCNu, signalTPCNu, float);
DECLARE_SOA_COLUMN(InnerParamTPCNu, innerParamTPCNu, float);
DECLARE_SOA_COLUMN(SignalTPCHad, signalTPCHad, float);
DECLARE_SOA_COLUMN(InnerParamTPCHad, innerParamTPCHad, float);
DECLARE_SOA_COLUMN(NClsTPCNu, nClsTPCNu, uint8_t);
DECLARE_SOA_COLUMN(NClsTPCHad, nClsTPCHad, uint8_t);
DECLARE_SOA_COLUMN(NCrossedRowsTPCNu, nCrossedRowsTPCNu, uint8_t);
DECLARE_SOA_COLUMN(NCrossedRowsTPCHad, nCrossedRowsTPCHad, uint8_t);
DECLARE_SOA_COLUMN(NSigmaTPCNu, nSigmaTPCNu, float);
DECLARE_SOA_COLUMN(NSigmaTPCHad, nSigmaTPCHad, float);
DECLARE_SOA_COLUMN(NSigmaTOFNu, nSigmaTOFNu, float);
DECLARE_SOA_COLUMN(NSigmaITSNu, nSigmaITSNu, float);
DECLARE_SOA_COLUMN(NSigmaTOFHad, nSigmaTOFHad, float);
DECLARE_SOA_COLUMN(NSigmaITSHad, nSigmaITSHad, float);
DECLARE_SOA_COLUMN(NSigmaTPCHadPi, nSigmaTPCHadPi, float);
DECLARE_SOA_COLUMN(NSigmaTPCHadKa, nSigmaTPCHadKa, float);
DECLARE_SOA_COLUMN(NSigmaTPCHadPr, nSigmaTPCHadPr, float);
DECLARE_SOA_COLUMN(NSigmaTOFHadPi, nSigmaTOFHadPi, float);
DECLARE_SOA_COLUMN(NSigmaTOFHadKa, nSigmaTOFHadKa, float);
DECLARE_SOA_COLUMN(NSigmaTOFHadPr, nSigmaTOFHadPr, float);
DECLARE_SOA_COLUMN(Chi2TPCNu, chi2TPCNu, float);
DECLARE_SOA_COLUMN(Chi2TPCHad, chi2TPCHad, float);
DECLARE_SOA_COLUMN(MassTOFNu, massTOFNu, float);
DECLARE_SOA_COLUMN(MassTOFHad, massTOFHad, float);
DECLARE_SOA_COLUMN(PidTrkNu, pidTrkNu, uint32_t);
DECLARE_SOA_COLUMN(PidTrkHad, pidTrkHad, uint32_t);
DECLARE_SOA_COLUMN(TrackIDHad, trackIDHad, int);
DECLARE_SOA_COLUMN(TrackIDNu, trackIDNu, int);

DECLARE_SOA_COLUMN(ItsClusterSizeNu, itsClusterSizeNu, uint32_t);
DECLARE_SOA_COLUMN(ItsClusterSizeHad, itsClusterSizeHad, uint32_t);

DECLARE_SOA_COLUMN(SharedClustersNu, sharedClustersNu, uint8_t);
DECLARE_SOA_COLUMN(SharedClustersHad, sharedClustersHad, uint8_t);

DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);

// Reconstructed-MC pair information. The signed generated pT follows the
// convention used by PtNu/PtHad: particles are positive and antiparticles
// are negative.
DECLARE_SOA_COLUMN(PtNuMC, ptNuMC, float);
DECLARE_SOA_COLUMN(EtaNuMC, etaNuMC, float);
DECLARE_SOA_COLUMN(PhiNuMC, phiNuMC, float);
DECLARE_SOA_COLUMN(PtHadMC, ptHadMC, float);
DECLARE_SOA_COLUMN(EtaHadMC, etaHadMC, float);
DECLARE_SOA_COLUMN(PhiHadMC, phiHadMC, float);
DECLARE_SOA_COLUMN(KstarMC, kstarMC, float);
DECLARE_SOA_COLUMN(PdgCodeNuMC, pdgCodeNuMC, int32_t);
DECLARE_SOA_COLUMN(PdgCodeHadMC, pdgCodeHadMC, int32_t);
DECLARE_SOA_COLUMN(IsPhysicalPrimaryNuMC, isPhysicalPrimaryNuMC, bool);
DECLARE_SOA_COLUMN(IsPhysicalPrimaryHadMC, isPhysicalPrimaryHadMC, bool);
DECLARE_SOA_COLUMN(SameMCCollision, sameMCCollision, bool);
DECLARE_SOA_COLUMN(MatchesRecoMCCollision, matchesRecoMCCollision, bool);

DECLARE_SOA_COLUMN(IsBkgUS, isBkgUS, bool);
DECLARE_SOA_COLUMN(IsBkgEM, isBkgEM, bool);

DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(ZVertex, zVertex, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, uint16_t);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(MultiplicityFT0C, multiplicityFT0C, float);

} // namespace hadron_nuclei_tables

DECLARE_SOA_TABLE(HadronNucleiTable, "AOD", "HADNUCLEITABLE",
                  hadron_nuclei_tables::PtNu,
                  hadron_nuclei_tables::EtaNu,
                  hadron_nuclei_tables::PhiNu,
                  hadron_nuclei_tables::PtHad,
                  hadron_nuclei_tables::EtaHad,
                  hadron_nuclei_tables::PhiHad,
                  hadron_nuclei_tables::DcaxyNu,
                  hadron_nuclei_tables::DcazNu,
                  hadron_nuclei_tables::DcaxyHad,
                  hadron_nuclei_tables::DcazHad,
                  hadron_nuclei_tables::DcaPair,
                  hadron_nuclei_tables::SignalTPCNu,
                  hadron_nuclei_tables::InnerParamTPCNu,
                  hadron_nuclei_tables::SignalTPCHad,
                  hadron_nuclei_tables::InnerParamTPCHad,
                  hadron_nuclei_tables::NClsTPCNu,
                  hadron_nuclei_tables::NClsTPCHad,
                  hadron_nuclei_tables::NCrossedRowsTPCNu,
                  hadron_nuclei_tables::NCrossedRowsTPCHad,
                  hadron_nuclei_tables::NSigmaTPCNu,
                  hadron_nuclei_tables::NSigmaTPCHadPi,
                  hadron_nuclei_tables::NSigmaTPCHadKa,
                  hadron_nuclei_tables::NSigmaTPCHadPr,
                  hadron_nuclei_tables::NSigmaTOFHadPi,
                  hadron_nuclei_tables::NSigmaTOFHadKa,
                  hadron_nuclei_tables::NSigmaTOFHadPr,
                  hadron_nuclei_tables::Chi2TPCNu,
                  hadron_nuclei_tables::Chi2TPCHad,
                  hadron_nuclei_tables::MassTOFNu,
                  hadron_nuclei_tables::MassTOFHad,
                  hadron_nuclei_tables::PidTrkNu,
                  hadron_nuclei_tables::PidTrkHad,
                  hadron_nuclei_tables::ItsClusterSizeNu,
                  hadron_nuclei_tables::ItsClusterSizeHad,
                  hadron_nuclei_tables::SharedClustersNu,
                  hadron_nuclei_tables::SharedClustersHad,
                  hadron_nuclei_tables::DeltaEta,
                  hadron_nuclei_tables::DeltaPhi,
                  hadron_nuclei_tables::NSigmaTPCHad,
                  hadron_nuclei_tables::NSigmaTOFNu,
                  hadron_nuclei_tables::NSigmaITSNu,
                  hadron_nuclei_tables::NSigmaTOFHad,
                  hadron_nuclei_tables::NSigmaITSHad)
DECLARE_SOA_TABLE(HadronNucleiTableMC, "AOD", "HADNUCLEIMC",
                  hadron_nuclei_tables::PtNuMC,
                  hadron_nuclei_tables::EtaNuMC,
                  hadron_nuclei_tables::PhiNuMC,
                  hadron_nuclei_tables::PtHadMC,
                  hadron_nuclei_tables::EtaHadMC,
                  hadron_nuclei_tables::PhiHadMC,
                  hadron_nuclei_tables::KstarMC,
                  hadron_nuclei_tables::PdgCodeNuMC,
                  hadron_nuclei_tables::PdgCodeHadMC,
                  hadron_nuclei_tables::IsPhysicalPrimaryNuMC,
                  hadron_nuclei_tables::IsPhysicalPrimaryHadMC,
                  hadron_nuclei_tables::SameMCCollision,
                  hadron_nuclei_tables::MatchesRecoMCCollision)
DECLARE_SOA_TABLE(HadronNucleiMult, "AOD", "HADNUCLEIMULT",
                  hadron_nuclei_tables::CollisionId,
                  hadron_nuclei_tables::ZVertex,
                  hadron_nuclei_tables::Multiplicity,
                  hadron_nuclei_tables::CentFT0C,
                  hadron_nuclei_tables::MultiplicityFT0C)

// Reduced hadron-hypertriton pair output. Reconstructable quantities such as
// k*, invariant mass and pointing variables are intentionally derived offline.
namespace hadronhyperpair
{
DECLARE_SOA_COLUMN(IsMixed, isMixed, bool);
DECLARE_SOA_COLUMN(MixingDepth, mixingDepth, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(TrackOccupancy, trackOccupancy, int);
DECLARE_SOA_COLUMN(Ft0cOccupancy, ft0cOccupancy, float);
DECLARE_SOA_COLUMN(NContributors, nContributors, uint16_t);
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
DECLARE_SOA_COLUMN(HypTofMass, hypTofMass, float);
DECLARE_SOA_COLUMN(HypNTPCCrossedRowsHe, hypNTPCCrossedRowsHe, uint8_t);
DECLARE_SOA_COLUMN(HypNTPCCrossedRowsPi, hypNTPCCrossedRowsPi, uint8_t);
DECLARE_SOA_COLUMN(HypTpcMomHe, hypTpcMomHe, float);
DECLARE_SOA_COLUMN(HypTpcMomPi, hypTpcMomPi, float);
DECLARE_SOA_COLUMN(HypTpcSignalHe, hypTpcSignalHe, uint16_t);
DECLARE_SOA_COLUMN(HypTpcSignalPi, hypTpcSignalPi, uint16_t);
DECLARE_SOA_COLUMN(HypTpcChi2He, hypTpcChi2He, float);
DECLARE_SOA_COLUMN(HypItsChi2He, hypItsChi2He, float);
DECLARE_SOA_COLUMN(HypItsChi2Pi, hypItsChi2Pi, float);
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
DECLARE_SOA_COLUMN(HadTpcNSigmaKa, hadTpcNSigmaKa, float);
DECLARE_SOA_COLUMN(HadTpcNSigmaPr, hadTpcNSigmaPr, float);
DECLARE_SOA_COLUMN(HadTofNSigmaPi, hadTofNSigmaPi, float);
DECLARE_SOA_COLUMN(HadTofNSigmaKa, hadTofNSigmaKa, float);
DECLARE_SOA_COLUMN(HadTofNSigmaPr, hadTofNSigmaPr, float);
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
} // namespace hadronhyperpair

DECLARE_SOA_TABLE(HadronHyperTable, "AOD", "HADHYPERTABLE",
                  hadronhyperpair::IsMixed,
                  hadronhyperpair::MixingDepth,
                  hadronhyperpair::RunNumber,
                  hadronhyperpair::PosZ,
                  hadronhyperpair::CentFT0A,
                  hadronhyperpair::CentFT0C,
                  hadronhyperpair::CentFT0M,
                  hadronhyperpair::TrackOccupancy,
                  hadronhyperpair::Ft0cOccupancy,
                  hadronhyperpair::NContributors,
                  hadronhyperpair::MultFT0C,
                  hadronhyperpair::XPrimVtx,
                  hadronhyperpair::YPrimVtx,
                  hadronhyperpair::ZPrimVtx,
                  hadronhyperpair::HypIsMatter,
                  hadronhyperpair::HypPtHe3,
                  hadronhyperpair::HypEtaHe3,
                  hadronhyperpair::HypPhiHe3,
                  hadronhyperpair::HypPtPi,
                  hadronhyperpair::HypEtaPi,
                  hadronhyperpair::HypPhiPi,
                  hadronhyperpair::HypDcaV0Daug,
                  hadronhyperpair::HypDcaHe,
                  hadronhyperpair::HypDcaPi,
                  hadronhyperpair::HypNSigmaHe,
                  hadronhyperpair::HypTofMass,
                  hadronhyperpair::HypNTPCCrossedRowsHe,
                  hadronhyperpair::HypNTPCCrossedRowsPi,
                  hadronhyperpair::HypTpcMomHe,
                  hadronhyperpair::HypTpcMomPi,
                  hadronhyperpair::HypTpcSignalHe,
                  hadronhyperpair::HypTpcSignalPi,
                  hadronhyperpair::HypTpcChi2He,
                  hadronhyperpair::HypItsChi2He,
                  hadronhyperpair::HypItsChi2Pi,
                  hadronhyperpair::HypItsClusterSizesHe,
                  hadronhyperpair::HypItsClusterSizesPi,
                  hadronhyperpair::HypXDecVtx,
                  hadronhyperpair::HypYDecVtx,
                  hadronhyperpair::HypZDecVtx,
                  hadronhyperpair::HadPt,
                  hadronhyperpair::HadEta,
                  hadronhyperpair::HadPhi,
                  hadronhyperpair::HadSign,
                  hadronhyperpair::HadDcaXY,
                  hadronhyperpair::HadTpcNClsCrossedRows,
                  hadronhyperpair::HadTpcNClsPID,
                  hadronhyperpair::HadTpcChi2NCl,
                  hadronhyperpair::HadItsClusterSizes,
                  hadronhyperpair::HadItsChi2NCl,
                  hadronhyperpair::HadHasTOF,
                  hadronhyperpair::HadTpcNSigmaPi,
                  hadronhyperpair::HadTpcNSigmaKa,
                  hadronhyperpair::HadTpcNSigmaPr,
                  hadronhyperpair::HadTofNSigmaPi,
                  hadronhyperpair::HadTofNSigmaKa,
                  hadronhyperpair::HadTofNSigmaPr);

DECLARE_SOA_TABLE(HadronHyperTableMC, "AOD", "HADHYPERMC",
                  hadronhyperpair::IsMixed,
                  hadronhyperpair::MixingDepth,
                  hadronhyperpair::RunNumber,
                  hadronhyperpair::PosZ,
                  hadronhyperpair::CentFT0A,
                  hadronhyperpair::CentFT0C,
                  hadronhyperpair::CentFT0M,
                  hadronhyperpair::TrackOccupancy,
                  hadronhyperpair::Ft0cOccupancy,
                  hadronhyperpair::NContributors,
                  hadronhyperpair::MultFT0C,
                  hadronhyperpair::XPrimVtx,
                  hadronhyperpair::YPrimVtx,
                  hadronhyperpair::ZPrimVtx,
                  hadronhyperpair::HypIsMatter,
                  hadronhyperpair::HypPtHe3,
                  hadronhyperpair::HypEtaHe3,
                  hadronhyperpair::HypPhiHe3,
                  hadronhyperpair::HypPtPi,
                  hadronhyperpair::HypEtaPi,
                  hadronhyperpair::HypPhiPi,
                  hadronhyperpair::HypDcaV0Daug,
                  hadronhyperpair::HypDcaHe,
                  hadronhyperpair::HypDcaPi,
                  hadronhyperpair::HypNSigmaHe,
                  hadronhyperpair::HypTofMass,
                  hadronhyperpair::HypNTPCCrossedRowsHe,
                  hadronhyperpair::HypNTPCCrossedRowsPi,
                  hadronhyperpair::HypTpcMomHe,
                  hadronhyperpair::HypTpcMomPi,
                  hadronhyperpair::HypTpcSignalHe,
                  hadronhyperpair::HypTpcSignalPi,
                  hadronhyperpair::HypTpcChi2He,
                  hadronhyperpair::HypItsChi2He,
                  hadronhyperpair::HypItsChi2Pi,
                  hadronhyperpair::HypItsClusterSizesHe,
                  hadronhyperpair::HypItsClusterSizesPi,
                  hadronhyperpair::HypXDecVtx,
                  hadronhyperpair::HypYDecVtx,
                  hadronhyperpair::HypZDecVtx,
                  hadronhyperpair::HadPt,
                  hadronhyperpair::HadEta,
                  hadronhyperpair::HadPhi,
                  hadronhyperpair::HadSign,
                  hadronhyperpair::HadDcaXY,
                  hadronhyperpair::HadTpcNClsCrossedRows,
                  hadronhyperpair::HadTpcNClsPID,
                  hadronhyperpair::HadTpcChi2NCl,
                  hadronhyperpair::HadItsClusterSizes,
                  hadronhyperpair::HadItsChi2NCl,
                  hadronhyperpair::HadHasTOF,
                  hadronhyperpair::HadTpcNSigmaPi,
                  hadronhyperpair::HadTpcNSigmaKa,
                  hadronhyperpair::HadTpcNSigmaPr,
                  hadronhyperpair::HadTofNSigmaPi,
                  hadronhyperpair::HadTofNSigmaKa,
                  hadronhyperpair::HadTofNSigmaPr,
                  hadronhyperpair::HypGenPt,
                  hadronhyperpair::HypGenEta,
                  hadronhyperpair::HypGenPhi,
                  hadronhyperpair::HypGenXDecVtx,
                  hadronhyperpair::HypGenYDecVtx,
                  hadronhyperpair::HypGenZDecVtx,
                  hadronhyperpair::HypIsReco,
                  hadronhyperpair::HypIsSignal,
                  hadronhyperpair::HypIsRecoMCCollision,
                  hadronhyperpair::HypIsSurvEvSel,
                  hadronhyperpair::HypIsTwoBodyDecay,
                  hadronhyperpair::HypStatusCode,
                  hadronhyperpair::HeGenPt,
                  hadronhyperpair::HeIsPhysicalPrimary,
                  hadronhyperpair::DecayPiIsPhysicalPrimary,
                  hadronhyperpair::HadRecoPt,
                  hadronhyperpair::HadRecoEta,
                  hadronhyperpair::HadRecoPhi,
                  hadronhyperpair::HadGenPt,
                  hadronhyperpair::HadGenEta,
                  hadronhyperpair::HadGenPhi,
                  hadronhyperpair::HadIsPhysicalPrimary,
                  hadronhyperpair::HadProcess,
                  hadronhyperpair::SameMCCollision,
                  hadronhyperpair::MatchesHypRecoMCCollision,
                  hadronhyperpair::MatchesPairRecoMCCollision,
                  hadronhyperpair::IsTruthSelfCorrelation,
                  hadronhyperpair::IsTruePrimaryHadHyperPair);

} // namespace o2::aod

#endif // PWGCF_FEMTO_FEMTONUCLEI_DATAMODEL_HADRONNUCLEITABLES_H_
