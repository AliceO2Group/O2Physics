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

#ifndef PWGLF_DATAMODEL_LFNUCLEITABLES_H_
#define PWGLF_DATAMODEL_LFNUCLEITABLES_H_

#include "Common/CCDB/EventSelectionParams.h"

#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace fullEvent
{ // Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(CentFV0M, centFV0M, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(MultNTracksPVeta1, multNTracksPVeta1, int);
DECLARE_SOA_DYNAMIC_COLUMN(IsInelGt0, isInelGt0, // is INEL > 0
                           [](int multPveta1) -> bool { return multPveta1 > 0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsInelGt1, isInelGt1, // is INEL > 1
                           [](int multPveta1) -> bool { return multPveta1 > 1; });
DECLARE_SOA_DYNAMIC_COLUMN(Selection_bit, selection_bit, // o2-linter: disable=name/o2-column (temporary fix)
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
                  fullEvent::MultNTracksPVeta1,
                  fullEvent::IsInelGt0<fullEvent::MultNTracksPVeta1>,
                  fullEvent::IsInelGt1<fullEvent::MultNTracksPVeta1>,
                  fullEvent::IsEventReject,
                  fullEvent::RunNumber,
                  fullEvent::Selection_bit<>);
using LfNuclEvent = LfNuclEvents::iterator;

namespace full
{
DECLARE_SOA_INDEX_COLUMN(LfNuclEvent, lfNuclEvent);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Sign, sign, int16_t);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity,
                           [](float pt, float eta, float mass) -> float {
                             const auto p = pt * std::cosh(eta);
                             const auto pz = pt * std::sinh(eta);
                             const auto energy = std::sqrt(p * p + mass * mass);
                             return 0.5f * std::log((energy + pz) / (energy - pz));
                           });
// ITS
DECLARE_SOA_COLUMN(ItsClusterSizes, itsClusterSizes, uint32_t); //! ITS cluster sizes per layer
// TPC
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TpcNSigmaPr, tpcNSigmaPr, float);
DECLARE_SOA_COLUMN(TpcNSigmaDe, tpcNSigmaDe, float);
DECLARE_SOA_COLUMN(TpcNSigmaTr, tpcNSigmaTr, float);
DECLARE_SOA_COLUMN(TpcNSigmaHe, tpcNSigmaHe, float);
DECLARE_SOA_COLUMN(TpcNSigmaAl, tpcNSigmaAl, float);
// TOF
DECLARE_SOA_COLUMN(TofNSigmaPi, tofNSigmaPi, float);
DECLARE_SOA_COLUMN(TofNSigmaKa, tofNSigmaKa, float);
DECLARE_SOA_COLUMN(TofNSigmaPr, tofNSigmaPr, float);
DECLARE_SOA_COLUMN(TofNSigmaDe, tofNSigmaDe, float);
DECLARE_SOA_COLUMN(TofNSigmaTr, tofNSigmaTr, float);
DECLARE_SOA_COLUMN(TofNSigmaHe, tofNSigmaHe, float);
DECLARE_SOA_COLUMN(TofNSigmaAl, tofNSigmaAl, float);
DECLARE_SOA_COLUMN(TpcExpSignalDiffPr, tpcExpSignalDiffPr, float);
DECLARE_SOA_COLUMN(TpcExpSignalDiffDe, tpcExpSignalDiffDe, float);
DECLARE_SOA_COLUMN(TpcExpSignalDiffHe, tpcExpSignalDiffHe, float);
DECLARE_SOA_COLUMN(TofExpSignalDiffPr, tofExpSignalDiffPr, float);
DECLARE_SOA_COLUMN(TofExpSignalDiffDe, tofExpSignalDiffDe, float);
DECLARE_SOA_COLUMN(TofExpSignalDiffHe, tofExpSignalDiffHe, float);
DECLARE_SOA_COLUMN(IsEvTimeTOF, isEvTimeTOF, bool);
DECLARE_SOA_COLUMN(IsEvTimeT0AC, isEvTimeT0AC, bool);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(TpcInnerParam, tpcInnerParam, float);
DECLARE_SOA_COLUMN(TofExpMom, tofExpMom, float);
DECLARE_SOA_COLUMN(TpcSignal, tpcSignal, float);
DECLARE_SOA_COLUMN(Beta, beta, float);
// TPC and ITS QA
DECLARE_SOA_COLUMN(PidForTracking, pidForTracking, uint8_t);
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, int16_t);
DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(ItsChi2NCl, itsChi2NCl, float);
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
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaPi, tpcNSigmaPi,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaKa, tpcNSigmaKa,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaPr, tpcNSigmaPr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaTr, tpcNSigmaTr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaAl, tpcNSigmaAl,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaPi, tofNSigmaPi,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaKa, tofNSigmaKa,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaPr, tofNSigmaPr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaTr, tofNSigmaTr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaAl, tofNSigmaAl,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcExpSignalDiffPr, tpcExpSignalDiffPr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcExpSignalDiffDe, tpcExpSignalDiffDe,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TpcExpSignalDiffHe, tpcExpSignalDiffHe,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofExpSignalDiffPr, tofExpSignalDiffPr,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofExpSignalDiffDe, tofExpSignalDiffDe,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofExpSignalDiffHe, tofExpSignalDiffHe,
                           [](bool /*b*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(TofExpMom, tofExpMom,
                           [](bool /*b*/) -> float { return 0.f; });
} // namespace dummy

DECLARE_SOA_TABLE(LfCandNucleus, "AOD", "LFNUCL",
                  o2::soa::Index<>,
                  full::LfNuclEventId,
                  full::DcaXY, full::DcaZ,
                  full::TpcNSigmaDe, full::TpcNSigmaHe,
                  full::TofNSigmaDe, full::TofNSigmaHe,
                  full::IsEvTimeTOF,
                  full::IsEvTimeT0AC,
                  full::HasTOF,
                  full::HasTRD,
                  full::TpcInnerParam,
                  full::Beta,
                  full::PidForTracking,
                  full::TpcSignal,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Sign,
                  full::ItsNCls,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  full::TpcChi2NCl,
                  full::ItsChi2NCl,
                  track::ITSClusterMap,
                  full::IsPVContributor,
                  full::P<full::Pt, full::Eta>,
                  full::Rapidity<full::Pt, full::Eta>,
                  full::ItsClusterSizes,
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>);
DECLARE_SOA_TABLE_VERSIONED(LfCandNucleusDummy, "AOD", "LFNUCL", 1,
                            o2::soa::Index<>,
                            full::LfNuclEventId,
                            full::DcaXY, full::DcaZ,
                            full::TpcNSigmaDe, full::TpcNSigmaHe,
                            full::TofNSigmaDe, full::TofNSigmaHe,
                            full::IsEvTimeTOF,
                            full::IsEvTimeT0AC,
                            full::HasTOF,
                            full::HasTRD,
                            full::TpcInnerParam,
                            full::Beta,
                            full::PidForTracking,
                            full::TpcSignal,
                            full::Pt,
                            full::Eta,
                            full::Phi,
                            full::Sign,
                            full::ItsNCls,
                            track::TPCNClsFindable,
                            track::TPCNClsFindableMinusFound,
                            track::TPCNClsFindableMinusCrossedRows,
                            full::TpcChi2NCl,
                            full::ItsChi2NCl,
                            track::ITSClusterMap,
                            full::IsPVContributor,
                            full::P<full::Pt, full::Eta>,
                            full::ItsClusterSizes,
                            dummy::TpcNSigmaPi<full::HasTOF>, dummy::TpcNSigmaKa<full::HasTOF>, dummy::TpcNSigmaPr<full::HasTOF>,
                            dummy::TpcNSigmaTr<full::HasTOF>, dummy::TpcNSigmaAl<full::HasTOF>,
                            dummy::TofNSigmaPi<full::HasTOF>, dummy::TofNSigmaKa<full::HasTOF>, dummy::TofNSigmaPr<full::HasTOF>,
                            dummy::TofNSigmaTr<full::HasTOF>, dummy::TofNSigmaAl<full::HasTOF>,
                            dummy::TpcExpSignalDiffPr<full::HasTOF>, dummy::TpcExpSignalDiffDe<full::HasTOF>, dummy::TpcExpSignalDiffHe<full::HasTOF>,
                            dummy::TofExpSignalDiffPr<full::HasTOF>, dummy::TofExpSignalDiffDe<full::HasTOF>, dummy::TofExpSignalDiffHe<full::HasTOF>,
                            dummy::TofExpMom<full::HasTOF>,
                            full::Rapidity<full::Pt, full::Eta>,
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>);

DECLARE_SOA_TABLE(LfCandNucleusExtra, "AOD", "LFNUCLEXTRA",
                  full::TpcNSigmaPi, full::TpcNSigmaKa, full::TpcNSigmaPr,
                  full::TpcNSigmaTr, full::TpcNSigmaAl,
                  full::TofNSigmaPi, full::TofNSigmaKa, full::TofNSigmaPr,
                  full::TofNSigmaTr, full::TofNSigmaAl,
                  full::TpcExpSignalDiffPr, full::TpcExpSignalDiffDe, full::TpcExpSignalDiffHe,
                  full::TofExpSignalDiffPr, full::TofExpSignalDiffDe, full::TofExpSignalDiffHe,
                  full::TofExpMom);

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
