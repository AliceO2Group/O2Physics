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

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFNUCLEITABLES_H_
#define PWGLF_DATAMODEL_LFNUCLEITABLES_H_

using namespace o2;

namespace o2::aod
{
namespace fullEvent
{                                 // Events
DECLARE_SOA_INDEX_COLUMN(BC, bc); //! Most probably BC to where this collision has occurred
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(MultFV0M, multFV0M, float);
} // namespace fullEvent
DECLARE_SOA_TABLE(LfCandNucleusEvents, "AOD", "LFNUCLEvent",
                  o2::soa::Index<>,
                  fullEvent::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  fullEvent::MultFV0M,
                  fullEvent::IsEventReject,
                  fullEvent::RunNumber);
using LfCandNucleusEvent = LfCandNucleusEvents::iterator;

namespace full
{
DECLARE_SOA_INDEX_COLUMN(LfCandNucleusEvent, lfCandNucleusFullEvent);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,
                           [](float px, float py) -> float { return sqrt(px * px + py * py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p,
                           [](float px, float py, float pz) -> float { return sqrt(px * px + py * py + pz * pz); });
DECLARE_SOA_COLUMN(Sign, sign, short);
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta,
                           [](float px, float py, float pz) -> float {
                             const auto p = sqrt(px * px + py * py + pz * pz);
                             return 0.5f * log((p + pz) / (p - pz));
                           });
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity,
                           [](float px, float py, float pz, float mass) -> float {
                             const auto energy = sqrt(px * px + py * py + pz * pz + mass * mass);
                             return 0.5f * log((energy + pz) / (energy - pz));
                           });
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
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, int16_t);
DECLARE_SOA_COLUMN(TPCChi2Ncl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);
// For MC
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
} // namespace full

/*
namespace fullMC
{
DECLARE_SOA_INDEX_COLUMN(LfCandNucleusEvent, lfCandNucleusFullEvent);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
}
*/

DECLARE_SOA_TABLE(LfCandNucleus, "AOD", "LFNUCL",
                  o2::soa::Index<>,
                  full::LfCandNucleusEventId,
                  full::DcaXY, full::DcaZ,
                  full::TPCNSigmaDe, full::TPCNSigmaHe,
                  full::TOFNSigmaDe, full::TOFNSigmaHe,
                  full::IsEvTimeTOF,
                  full::IsEvTimeT0AC,
                  full::HasTOF,
                  full::HasTRD,
                  full::TPCInnerParam,
                  full::Beta,
                  full::TPCSignal,
                  full::Px, full::Py, full::Pz,
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
                  full::Pt<full::Px, full::Py>,
                  full::P<full::Px, full::Py, full::Pz>,
                  full::Eta<full::Px, full::Py, full::Pz>,
                  full::Rapidity<full::Px, full::Py, full::Pz>,
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
                  full::ProducedByGenerator);

using LfCandNucleusFull = soa::Join<LfCandNucleus, LfCandNucleusExtra>;

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFNUCLEITABLES_H_
