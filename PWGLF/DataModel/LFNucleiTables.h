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
/// \file LFNucleiTables.h
///
/// \author Rutuparna Rath <rutuparna.rath@cern.ch> and Giovanni Malfattore <giovanni.malfattore@cern.ch>
///

#include "Common/CCDB/EventSelectionParams.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGLF_DATAMODEL_LFNUCLEITABLES_H_
#define PWGLF_DATAMODEL_LFNUCLEITABLES_H_

namespace o2::aod
{
namespace fullEvent
{ // Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(CentFV0M, centFV0M, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_DYNAMIC_COLUMN(Selection_Bit, selection_bit, //! Dummy
                           [](o2::aod::evsel::EventSelectionFlags /*v*/) -> bool { return true; });
} // namespace fullEvent
DECLARE_SOA_TABLE(LfNuclEvents, "AOD", "LFNUCLEvent",
                  o2::soa::Index<>,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  fullEvent::CentFV0M,
                  fullEvent::CentFT0M,
                  fullEvent::IsEventReject,
                  fullEvent::RunNumber,
                  fullEvent::Selection_Bit<>);
using LfNuclEvent = LfNuclEvents::iterator;

namespace full
{
DECLARE_SOA_INDEX_COLUMN(LfNuclEvent, lfNuclEvent);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * cosh(eta); });
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Sign, sign, int16_t);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity,
                           [](float pt, float eta, float mass) -> float {
                             const auto p = pt * cosh(eta);
                             const auto pz = pt * sinh(eta);
                             const auto energy = sqrt(p * p + mass * mass);
                             return 0.5f * log((energy + pz) / (energy - pz));
                           });
// ITS
DECLARE_SOA_COLUMN(ITSClusterSizes, itsClusterSizes, uint32_t); //! ITS cluster sizes per layer
// TPC
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float);
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float);
DECLARE_SOA_COLUMN(TPCNSigmaTr, tpcNSigmaTr, float);
DECLARE_SOA_COLUMN(TPCNSigmaHe, tpcNSigmaHe, float);
DECLARE_SOA_COLUMN(TPCNSigmaAl, tpcNSigmaAl, float);
// TOF
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNSigmaPi, float);
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float);
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float);
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float);
DECLARE_SOA_COLUMN(TOFNSigmaTr, tofNSigmaTr, float);
DECLARE_SOA_COLUMN(TOFNSigmaHe, tofNSigmaHe, float);
DECLARE_SOA_COLUMN(TOFNSigmaAl, tofNSigmaAl, float);
DECLARE_SOA_COLUMN(TPCExpSignalDiffPr, tpcExpSignalDiffPr, float);
DECLARE_SOA_COLUMN(TPCExpSignalDiffDe, tpcExpSignalDiffDe, float);
DECLARE_SOA_COLUMN(TPCExpSignalDiffHe, tpcExpSignalDiffHe, float);
DECLARE_SOA_COLUMN(TOFExpSignalDiffPr, tofExpSignalDiffPr, float);
DECLARE_SOA_COLUMN(TOFExpSignalDiffDe, tofExpSignalDiffDe, float);
DECLARE_SOA_COLUMN(TOFExpSignalDiffHe, tofExpSignalDiffHe, float);
DECLARE_SOA_COLUMN(IsEvTimeTOF, isEvTimeTOF, bool);
DECLARE_SOA_COLUMN(IsEvTimeT0AC, isEvTimeT0AC, bool);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);
DECLARE_SOA_COLUMN(TOFExpMom, tofExpMom, float);
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);
DECLARE_SOA_COLUMN(Beta, beta, float);
// TPC and ITS QA
DECLARE_SOA_COLUMN(PIDForTracking, pidForTracking, uint8_t);
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, int16_t);
DECLARE_SOA_COLUMN(TPCChi2Ncl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);
DECLARE_SOA_COLUMN(TpcPassed, tpcPassed, bool);
DECLARE_SOA_COLUMN(ItsPassed, itsPassed, bool);
DECLARE_SOA_COLUMN(FakeHitsFlag, fakeHitsFlag, bool);

// For MC
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
DECLARE_SOA_COLUMN(GetProcess, getProcess, int);

} // namespace full
namespace dummy
{
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPi, tpcNSigmaPi,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa, tpcNSigmaKa,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaTr, tpcNSigmaTr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaAl, tpcNSigmaAl,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPi, tofNSigmaPi,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaKa, tofNSigmaKa,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaTr, tofNSigmaTr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaAl, tofNSigmaAl,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffPr, tpcExpSignalDiffPr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffDe, tpcExpSignalDiffDe,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCExpSignalDiffHe, tpcExpSignalDiffHe,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffPr, tofExpSignalDiffPr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffDe, tofExpSignalDiffDe,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpSignalDiffHe, tofExpSignalDiffHe,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TOFExpMom, tofExpMom,
                           [](bool /*b*/) -> float { return 0.f; });
} // namespace dummy

/*
namespace fullMC
{
DECLARE_SOA_INDEX_COLUMN(LfNuclEvent, lfCandNucleusFullEvent);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
}
*/

DECLARE_SOA_TABLE(LfCandNucleus, "AOD", "LFNUCL",
                  o2::soa::Index<>,
                  full::LfNuclEventId,
                  full::DcaXY, full::DcaZ,
                  full::TPCNSigmaDe, full::TPCNSigmaHe,
                  full::TOFNSigmaDe, full::TOFNSigmaHe,
                  full::IsEvTimeTOF,
                  full::IsEvTimeT0AC,
                  full::HasTOF,
                  full::HasTRD,
                  full::TPCInnerParam,
                  full::Beta,
                  full::PIDForTracking,
                  full::TPCSignal,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Sign,
                  full::ITSNCls,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  full::TPCChi2Ncl,
                  full::ITSChi2NCl,
                  track::ITSClusterMap,
                  full::IsPVContributor,
                  full::P<full::Pt, full::Eta>,
                  full::Rapidity<full::Pt, full::Eta>,
                  full::ITSClusterSizes,
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>);
DECLARE_SOA_TABLE_VERSIONED(LfCandNucleusDummy, "AOD", "LFNUCL", 1,
                            o2::soa::Index<>,
                            full::LfNuclEventId,
                            full::DcaXY, full::DcaZ,
                            full::TPCNSigmaDe, full::TPCNSigmaHe,
                            full::TOFNSigmaDe, full::TOFNSigmaHe,
                            full::IsEvTimeTOF,
                            full::IsEvTimeT0AC,
                            full::HasTOF,
                            full::HasTRD,
                            full::TPCInnerParam,
                            full::Beta,
                            full::PIDForTracking,
                            full::TPCSignal,
                            full::Pt,
                            full::Eta,
                            full::Phi,
                            full::Sign,
                            full::ITSNCls,
                            track::TPCNClsFindable,
                            track::TPCNClsFindableMinusFound,
                            track::TPCNClsFindableMinusCrossedRows,
                            full::TPCChi2Ncl,
                            full::ITSChi2NCl,
                            track::ITSClusterMap,
                            full::IsPVContributor,
                            full::P<full::Pt, full::Eta>,
                            full::ITSClusterSizes,
                            dummy::TPCNSigmaPi<full::HasTOF>, dummy::TPCNSigmaKa<full::HasTOF>, dummy::TPCNSigmaPr<full::HasTOF>,
                            dummy::TPCNSigmaTr<full::HasTOF>, dummy::TPCNSigmaAl<full::HasTOF>,
                            dummy::TOFNSigmaPi<full::HasTOF>, dummy::TOFNSigmaKa<full::HasTOF>, dummy::TOFNSigmaPr<full::HasTOF>,
                            dummy::TOFNSigmaTr<full::HasTOF>, dummy::TOFNSigmaAl<full::HasTOF>,
                            dummy::TPCExpSignalDiffPr<full::HasTOF>, dummy::TPCExpSignalDiffDe<full::HasTOF>, dummy::TPCExpSignalDiffHe<full::HasTOF>,
                            dummy::TOFExpSignalDiffPr<full::HasTOF>, dummy::TOFExpSignalDiffDe<full::HasTOF>, dummy::TOFExpSignalDiffHe<full::HasTOF>,
                            dummy::TOFExpMom<full::HasTOF>,
                            full::Rapidity<full::Pt, full::Eta>,
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>);

DECLARE_SOA_TABLE(LfCandNucleusExtra, "AOD", "LFNUCLEXTRA",
                  full::TPCNSigmaPi, full::TPCNSigmaKa, full::TPCNSigmaPr,
                  full::TPCNSigmaTr, full::TPCNSigmaAl,
                  full::TOFNSigmaPi, full::TOFNSigmaKa, full::TOFNSigmaPr,
                  full::TOFNSigmaTr, full::TOFNSigmaAl,
                  full::TPCExpSignalDiffPr, full::TPCExpSignalDiffDe, full::TPCExpSignalDiffHe,
                  full::TOFExpSignalDiffPr, full::TOFExpSignalDiffDe, full::TOFExpSignalDiffHe,
                  full::TOFExpMom);

DECLARE_SOA_TABLE(LfCandNucleusMC, "AOD", "LFNUCLMC",
                  mcparticle::PdgCode,
                  full::IsPhysicalPrimary,
                  full::ProducedByGenerator,
                  full::GetProcess,
                  full::ItsPassed,
                  full::TpcPassed,
                  mcparticle::Px,
                  mcparticle::Py,
                  mcparticle::Pz,
                  full::FakeHitsFlag);

using LfCandNucleusFull = soa::Join<LfCandNucleus, LfCandNucleusExtra>;

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFNUCLEITABLES_H_
