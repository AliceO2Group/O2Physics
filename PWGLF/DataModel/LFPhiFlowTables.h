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

/// \file LFPhiFlowTables.h
/// \brief DataModel for Phi flow
///
/// \author Prottay Das <prottay.das@cern.ch>

#ifndef PWGLF_DATAMODEL_LFPHIFLOWTABLES_H_
#define PWGLF_DATAMODEL_LFPHIFLOWTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace kaonevent
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Posz, posz, float);
DECLARE_SOA_COLUMN(QxA, qxA, float);
DECLARE_SOA_COLUMN(QxC, qxC, float);
DECLARE_SOA_COLUMN(QyA, qyA, float);
DECLARE_SOA_COLUMN(QyC, qyC, float);
} // namespace kaonevent
DECLARE_SOA_TABLE(KaonEvents, "AOD", "KAONEVENT",
                  o2::soa::Index<>,
                  kaonevent::Cent,
                  kaonevent::Posz,
                  kaonevent::QxA,
                  kaonevent::QxC,
                  kaonevent::QyA,
                  kaonevent::QyC)
using KaonEvent = KaonEvents::iterator;

namespace kaonpair
{
DECLARE_SOA_INDEX_COLUMN(KaonEvent, kaonevent);
DECLARE_SOA_COLUMN(Px, px, float);                      //! Bachelor Kaon Px
DECLARE_SOA_COLUMN(Py, py, float);                      //! Bachelor Kaon Py
DECLARE_SOA_COLUMN(Pz, pz, float);                      //! Bachelor Kaon Pz
DECLARE_SOA_COLUMN(Charge, charge, int8_t);             //! Charge
DECLARE_SOA_COLUMN(KaonIndex, kaonIndex, int64_t);      //! Daughter Kaon index1
DECLARE_SOA_COLUMN(KaonPidMask, kaonPidMask, uint16_t); //! bitmask for PID selections
} // namespace kaonpair

DECLARE_SOA_TABLE(KaonTracks, "AOD", "KAONTRACK",
                  o2::soa::Index<>,
                  kaonpair::KaonEventId,
                  kaonpair::Px,
                  kaonpair::Py,
                  kaonpair::Pz,
                  kaonpair::Charge,
                  kaonpair::KaonIndex,
                  kaonpair::KaonPidMask);

using KaonTrack = KaonTracks::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFPHIFLOWTABLES_H_
