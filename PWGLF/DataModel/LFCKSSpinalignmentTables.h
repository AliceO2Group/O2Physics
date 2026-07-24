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

/// \file LFCKSSpinalignmentTables.h
/// \brief DataModel for reduced K0s and pion tables for charged KStar analysis
///
/// \author Prottay Das <prottay.das@cern.ch>

#ifndef PWGLF_DATAMODEL_LFCKSSPINALIGNMENTTABLES_H_
#define PWGLF_DATAMODEL_LFCKSSPINALIGNMENTTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace kshortpionevent
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Posz, posz, float);
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);
DECLARE_SOA_COLUMN(PsiFT0A, psiFT0A, float);
DECLARE_SOA_COLUMN(PsiTPC, psiTPC, float);
} // namespace kshortpionevent

DECLARE_SOA_TABLE(KShortpionEvents, "AOD", "KSHORTPIONEVENT",
                  o2::soa::Index<>,
                  kshortpionevent::Cent,
                  kshortpionevent::Posz,
                  kshortpionevent::PsiFT0C,
                  kshortpionevent::PsiFT0A,
                  kshortpionevent::PsiTPC);

using KShortpionEvent = KShortpionEvents::iterator;

// ------------------------------------------------------------
// K0s candidate table
// ------------------------------------------------------------
namespace kshorttrack
{
DECLARE_SOA_INDEX_COLUMN(KShortpionEvent, kshortpionevent);

DECLARE_SOA_COLUMN(V0Cospa, v0Cospa, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(DcaPositive, dcaPositive, float);
DECLARE_SOA_COLUMN(DcaNegative, dcaNegative, float);
DECLARE_SOA_COLUMN(DcaBetweenDaughter, dcaBetweenDaughter, float);
// DECLARE_SOA_COLUMN(V0Lifetime, v0Lifetime, float);

DECLARE_SOA_COLUMN(KShortPx, kShortPx, float);
DECLARE_SOA_COLUMN(KShortPy, kShortPy, float);
DECLARE_SOA_COLUMN(KShortPz, kShortPz, float);
DECLARE_SOA_COLUMN(KShortMass, kShortMass, float);

DECLARE_SOA_COLUMN(PionIndex1, pionIndex1, int64_t);
DECLARE_SOA_COLUMN(PionIndex2, pionIndex2, int64_t);
} // namespace kshorttrack

DECLARE_SOA_TABLE(KShortTracks, "AOD", "KSHORTTRACK",
                  o2::soa::Index<>,
                  kshorttrack::KShortpionEventId,
                  kshorttrack::V0Cospa,
                  kshorttrack::V0Radius,
                  kshorttrack::DcaPositive,
                  kshorttrack::DcaNegative,
                  kshorttrack::DcaBetweenDaughter,
                  // kshorttrack::V0Lifetime,
                  kshorttrack::KShortPx,
                  kshorttrack::KShortPy,
                  kshorttrack::KShortPz,
                  kshorttrack::KShortMass,
                  kshorttrack::PionIndex1,
                  kshorttrack::PionIndex2);

using KShortTrack = KShortTracks::iterator;

// ------------------------------------------------------------
// Bachelor pion table
// ------------------------------------------------------------
namespace piontrack
{
DECLARE_SOA_INDEX_COLUMN(KShortpionEvent, kshortpionevent);

DECLARE_SOA_COLUMN(PionBachPx, pionBachPx, float);
DECLARE_SOA_COLUMN(PionBachPy, pionBachPy, float);
DECLARE_SOA_COLUMN(PionBachPz, pionBachPz, float);
DECLARE_SOA_COLUMN(Charge, charge, int8_t);
DECLARE_SOA_COLUMN(PionBachIndex, pionBachIndex, int64_t);
DECLARE_SOA_COLUMN(PionPidMask, pionPidMask, uint8_t);
} // namespace piontrack

DECLARE_SOA_TABLE(PionTracks, "AOD", "PIONTRACK",
                  o2::soa::Index<>,
                  piontrack::KShortpionEventId,
                  piontrack::PionBachPx,
                  piontrack::PionBachPy,
                  piontrack::PionBachPz,
                  piontrack::Charge,
                  piontrack::PionBachIndex,
                  piontrack::PionPidMask);

using PionTrack = PionTracks::iterator;

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFCKSSPINALIGNMENTTABLES_H_
