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

#ifndef O2_ANALYSIS_LFNUCLEITABLES_H_
#define O2_ANALYSIS_LFNUCLEITABLES_H_

using namespace o2;

namespace o2::aod
{
namespace fullEvent
{                                 // Events
DECLARE_SOA_INDEX_COLUMN(BC, bc); //! Most probably BC to where this collision has occured
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(Vz, vz, float);
DECLARE_SOA_COLUMN(V0M, v0m, int);
} // namespace fullEvent
DECLARE_SOA_TABLE(LfCandNucleusFullEvents, "AOD", "LFNUCLEvent",
                  o2::soa::Index<>,
                  fullEvent::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  fullEvent::V0M,
                  fullEvent::IsEventReject,
                  fullEvent::RunNumber);
using LfCandNucleusFullEvent = LfCandNucleusFullEvents::iterator;

namespace full
{
DECLARE_SOA_INDEX_COLUMN(LfCandNucleusFullEvent, lfCandNucleusFullEvent);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(NSigTPCPi, nsigTPCPi, float);
DECLARE_SOA_COLUMN(NSigTPCKa, nsigTPCKa, float);
DECLARE_SOA_COLUMN(NSigTPCPr, nsigTPCPr, float);
DECLARE_SOA_COLUMN(NSigTPCDe, nsigTPCD, float);
DECLARE_SOA_COLUMN(NSigTPC3He, nsigTPC3He, float);
DECLARE_SOA_COLUMN(NSigTOFPi, nsigTOFPi, float);
DECLARE_SOA_COLUMN(NSigTOFKa, nsigTOFKa, float);
DECLARE_SOA_COLUMN(NSigTOFPr, nsigTOFPr, float);
DECLARE_SOA_COLUMN(NSigTOFDe, nsigTOFD, float);
DECLARE_SOA_COLUMN(NSigTOF3He, nsigTOF3He, float);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(DCAxy, dcaxy, float);
DECLARE_SOA_COLUMN(DCAz, dcaz, float);
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);
DECLARE_SOA_COLUMN(Beta, beta, float);
// TPC and ITS QA
DECLARE_SOA_COLUMN(NcrTPC, ncrTPC, int);
DECLARE_SOA_COLUMN(RTPC, rTPC, float);
DECLARE_SOA_COLUMN(Chi2TPC, chi2TPC, float);
DECLARE_SOA_COLUMN(Chi2ITS, chi2ITS, float);
} // namespace full

/*
namespace fullMC
{
DECLARE_SOA_INDEX_COLUMN(LfCandNucleusFullEvent, lfCandNucleusFullEvent);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
}
*/

DECLARE_SOA_TABLE(LfCandNucleusFull, "AOD", "LFNUCL",
                  o2::soa::Index<>,
                  full::LfCandNucleusFullEventId,
                  full::DCAxy,
                  full::DCAz,
                  full::NSigTPCPi, full::NSigTPCKa, full::NSigTPCPr,
                  full::NSigTPCDe, full::NSigTPC3He,
                  full::NSigTOFPi, full::NSigTOFKa, full::NSigTOFPr,
                  full::NSigTOFDe, full::NSigTOF3He,
                  full::HasTOF,
                  full::TPCInnerParam,
                  full::TPCSignal,
                  full::Beta,
                  full::Px,
                  full::Py,
                  full::Pz,
                  full::Pt,
                  full::P,
                  full::Eta,
                  full::Phi,
                  full::Sign,
                  full::NcrTPC,
                  full::RTPC,
                  full::Chi2TPC,
                  full::Chi2ITS);
DECLARE_SOA_TABLE(LfCandNucleusMC, "AOD", "LFNUCLMC",
                  mcparticle::PdgCode);

} // namespace o2::aod
#endif
