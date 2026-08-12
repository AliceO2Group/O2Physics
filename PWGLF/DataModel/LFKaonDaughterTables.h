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

/// \file LFKaonDaughterTables.h
/// \brief DataModel for shared kaon/pion daughter derived data (K*(892), phi)
///
/// \author sarjeeta.gami@cern.ch

#ifndef PWGLF_DATAMODEL_LFKAONDAUGHTERTABLES_H_
#define PWGLF_DATAMODEL_LFKAONDAUGHTERTABLES_H_

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace kadaughterevent
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Posz, posz, float);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);
DECLARE_SOA_COLUMN(QFT0C, qFT0C, float);
DECLARE_SOA_COLUMN(PsiZDCC, psiZDCC, float); // reserved for a future phi/SP consumer; 0 for kstar
DECLARE_SOA_COLUMN(BcId, bcId, int);         // mirrors collision.bcId(), used to dedup mixed pairs sharing a BC
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(Timestamp, timestamp, int64_t); // for CCDB EP-weight fetch
} // namespace kadaughterevent
DECLARE_SOA_TABLE(KaDaughterEvents, "AOD", "KADEVENT", o2::soa::Index<>,
                  kadaughterevent::Cent, kadaughterevent::Posz, kadaughterevent::Occupancy,
                  kadaughterevent::PsiFT0C, kadaughterevent::QFT0C, kadaughterevent::PsiZDCC,
                  kadaughterevent::BcId, kadaughterevent::RunNumber, kadaughterevent::Timestamp);
using KaDaughterEvent = KaDaughterEvents::iterator;

namespace kadaughtertrack
{
DECLARE_SOA_INDEX_COLUMN(KaDaughterEvent, kaDaughterEvent);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(TrackId, trackId, int64_t); // original AOD track global index, for bookkeeping only
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TofNSigmaKa, tofNSigmaKa, float);
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TofNSigmaPi, tofNSigmaPi, float);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(Beta, beta, float);
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return std::hypot(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::constrainAngle(std::atan2(py, px)); });
} // namespace kadaughtertrack
DECLARE_SOA_TABLE(KaDaughterTracks, "AOD", "KADTRACK", o2::soa::Index<>,
                  kadaughtertrack::KaDaughterEventId, kadaughtertrack::Px, kadaughtertrack::Py,
                  kadaughtertrack::Pz, kadaughtertrack::Sign, kadaughtertrack::TrackId,
                  kadaughtertrack::DcaXY, kadaughtertrack::DcaZ,
                  kadaughtertrack::TpcNSigmaKa, kadaughtertrack::TofNSigmaKa,
                  kadaughtertrack::TpcNSigmaPi, kadaughtertrack::TofNSigmaPi,
                  kadaughtertrack::HasTOF, kadaughtertrack::Beta,
                  kadaughtertrack::Pt<kadaughtertrack::Px, kadaughtertrack::Py>,
                  kadaughtertrack::Phi<kadaughtertrack::Px, kadaughtertrack::Py>);
using KaDaughterTrack = KaDaughterTracks::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFKAONDAUGHTERTABLES_H_
