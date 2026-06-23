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
namespace kaonkaonevent
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Posz, posz, float);
DECLARE_SOA_COLUMN(QxA, qxA, float);
DECLARE_SOA_COLUMN(QxC, qxC, float);
DECLARE_SOA_COLUMN(QyA, qyA, float);
DECLARE_SOA_COLUMN(QyC, qyC, float);
} // namespace kaonkaonevent
DECLARE_SOA_TABLE(KaonkaonEvents, "AOD", "KAONKAONEVENT",
                  o2::soa::Index<>,
                  kaonkaonevent::Cent,
                  kaonkaonevent::Posz,
                  kaonkaonevent::QxA,
                  kaonkaonevent::QxC,
                  kaonkaonevent::QyA,
                  kaonkaonevent::QyC)
using KaonkaonEvent = KaonkaonEvents::iterator;

namespace kaonkaonpair
{
DECLARE_SOA_INDEX_COLUMN(KaonkaonEvent, kaonkaonevent);
DECLARE_SOA_COLUMN(D1Px, d1Px, float);                  //! Bachelor Kaon Px
DECLARE_SOA_COLUMN(D1Py, d1Py, float);                  //! Bachelor Kaon Py
DECLARE_SOA_COLUMN(D1Pz, d1Pz, float);                  //! Bachelor Kaon Pz
DECLARE_SOA_COLUMN(D2Px, d2Px, float);                  //! Bachelor Kaon Px
DECLARE_SOA_COLUMN(D2Py, d2Py, float);                  //! Bachelor Kaon Py
DECLARE_SOA_COLUMN(D2Pz, d2Pz, float);                  //! Bachelor Kaon Pz
DECLARE_SOA_COLUMN(PhiM, phiM, float);                  //! Phi Mass
DECLARE_SOA_COLUMN(KaonIndex1, kaonIndex1, int64_t);    //! Daughter Kaon index1
DECLARE_SOA_COLUMN(KaonIndex2, kaonIndex2, int64_t);    //! Daughter Kaon index2
DECLARE_SOA_COLUMN(KaonPidMask, kaonPidMask, uint16_t); //! bitmask for PID selections
} // namespace kaonkaonpair

DECLARE_SOA_TABLE(KaonTracks, "AOD", "KAONTRACK",
                  o2::soa::Index<>,
                  kaonkaonpair::KaonkaonEventId,
                  kaonkaonpair::D1Px,
                  kaonkaonpair::D1Py,
                  kaonkaonpair::D1Pz,
                  kaonkaonpair::D2Px,
                  kaonkaonpair::D2Py,
                  kaonkaonpair::D2Pz,
                  kaonkaonpair::PhiM,
                  kaonkaonpair::KaonIndex1,
                  kaonkaonpair::KaonIndex2,
                  kaonkaonpair::KaonPidMask);

using KaonTrack = KaonTracks::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFPHIFLOWTABLES_H_
