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
/// \brief DataModel for Charged KStar spin alignment
///
/// \author Prottay Das <prottay.das@cern.ch>

#ifndef PWGLF_DATAMODEL_LFCKSSPINALIGNMENTTABLES_H_
#define PWGLF_DATAMODEL_LFCKSSPINALIGNMENTTABLES_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace kshortpionevent
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Posz, posz, float);
// DECLARE_SOA_COLUMN(CollIndex, collIndex, float);
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);
DECLARE_SOA_COLUMN(PsiFT0A, psiFT0A, float);
DECLARE_SOA_COLUMN(PsiTPC, psiTPC, float);
} // namespace kshortpionevent
DECLARE_SOA_TABLE(KShortpionEvents, "AOD", "KSHORTPIONEVENT",
                  o2::soa::Index<>,
                  kshortpionevent::Cent,
                  kshortpionevent::Posz,
                  // kshortpionevent::CollIndex,
                  kshortpionevent::PsiFT0C,
                  kshortpionevent::PsiFT0A,
                  kshortpionevent::PsiTPC)
using KShortpionEvent = KShortpionEvents::iterator;

namespace kshortpionpair
{
DECLARE_SOA_INDEX_COLUMN(KShortpionEvent, kshortpionevent);
DECLARE_SOA_COLUMN(V0Cospa, v0Cospa, float);                       //! V0 Cospa
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);                     //! V0 Radius
DECLARE_SOA_COLUMN(DcaPositive, dcaPositive, float);               //! DCA Positive
DECLARE_SOA_COLUMN(DcaNegative, dcaNegative, float);               //! DCA Negative
DECLARE_SOA_COLUMN(DcaBetweenDaughter, dcaBetweenDaughter, float); //! DCA between daughters
DECLARE_SOA_COLUMN(V0Lifetime, v0Lifetime, float);                 //! KShort lifetime
DECLARE_SOA_COLUMN(KShortPx, kShortPx, float);                     //! KShort Px
DECLARE_SOA_COLUMN(KShortPy, kShortPy, float);                     //! KShort Py
DECLARE_SOA_COLUMN(KShortPz, kShortPz, float);                     //! KShort Pz
DECLARE_SOA_COLUMN(KShortMass, kShortMass, float);                 //! KShort Mass
DECLARE_SOA_COLUMN(PionBachPx, pionBachPx, float);                 //! Bachelor Pion Px
DECLARE_SOA_COLUMN(PionBachPy, pionBachPy, float);                 //! Bachelor Pion Py
DECLARE_SOA_COLUMN(PionBachPz, pionBachPz, float);                 //! Bachelor Pion Pz
/*DECLARE_SOA_COLUMN(PionBachTPC, pionBachTPC, float);               //! Bachelor Pion nsigmatpc
DECLARE_SOA_COLUMN(PionBachTOFHit, pionBachTOFHit, int);           //! Bachelor Pion tof hit availability
DECLARE_SOA_COLUMN(PionBachTOF, pionBachTOF, float);               //! Bachelor Pion nsigmatof*/
DECLARE_SOA_COLUMN(PionBachIndex, pionBachIndex, int); //! Bachelor Pion index
DECLARE_SOA_COLUMN(PionIndex1, pionIndex1, int);       //! Daughter Pion index1
DECLARE_SOA_COLUMN(PionIndex2, pionIndex2, int);       //! Daughter Pion index2
DECLARE_SOA_COLUMN(PionPidMask, pionPidMask, uint8_t); //! bitmask for PID selections
} // namespace kshortpionpair
DECLARE_SOA_TABLE(KShortTracks, "AOD", "KSHORTTRACK",
                  o2::soa::Index<>,
                  kshortpionpair::KShortpionEventId,
                  kshortpionpair::V0Cospa,
                  kshortpionpair::V0Radius,
                  kshortpionpair::DcaPositive,
                  kshortpionpair::DcaNegative,
                  kshortpionpair::DcaBetweenDaughter,
                  kshortpionpair::V0Lifetime,
                  // kshortpionpair::Armenteros,
                  kshortpionpair::KShortPx,
                  kshortpionpair::KShortPy,
                  kshortpionpair::KShortPz,
                  kshortpionpair::KShortMass,
                  kshortpionpair::PionIndex1,
                  kshortpionpair::PionIndex2);

using KShortTrack = KShortTracks::iterator;

DECLARE_SOA_TABLE(PionTracks, "AOD", "PIONTRACK",
                  o2::soa::Index<>,
                  kshortpionpair::KShortpionEventId,
                  kshortpionpair::PionBachPx,
                  kshortpionpair::PionBachPy,
                  kshortpionpair::PionBachPz,
                  // kshortpionpair::PionBachSign,
                  // kshortpionpair::PionBachTPC,
                  // kshortpionpair::PionBachTOFHit,
                  // kshortpionpair::PionBachTOF,
                  kshortpionpair::PionBachIndex,
                  kshortpionpair::PionPidMask);

using PionTrack = PionTracks::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFCKSSPINALIGNMENTTABLES_H_
