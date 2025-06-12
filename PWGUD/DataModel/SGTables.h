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

#ifndef PWGUD_DATAMODEL_SGTABLES_H_
#define PWGUD_DATAMODEL_SGTABLES_H_

#include <vector>
#include <cmath>
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/DataTypes.h"
#include "MathUtils/Utils.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

namespace o2::aod
{
namespace sgevent
{
DECLARE_SOA_COLUMN(Run, run, int32_t);
DECLARE_SOA_COLUMN(Flag, flag, int);
DECLARE_SOA_COLUMN(GS, gs, int);
DECLARE_SOA_COLUMN(ZNA, zna, float);
DECLARE_SOA_COLUMN(ZNC, znc, float);
DECLARE_SOA_COLUMN(Ntr, ntr, int);
DECLARE_SOA_COLUMN(Occ, occ, int);
DECLARE_SOA_COLUMN(Ir, ir, float);
} // namespace sgevent
DECLARE_SOA_TABLE(SGEvents, "AOD", "SGEVENT", // o2::soa::Index<>,
                  sgevent::Run, sgevent::Flag, sgevent::GS, sgevent::ZNA, sgevent::ZNC, sgevent::Ntr, sgevent::Occ, sgevent::Ir);
// sgevent::Run, sgevent::Flag);
using SGEvent = SGEvents::iterator;
namespace sgtrack
{
DECLARE_SOA_INDEX_COLUMN(SGEvent, sgEvent);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(TPCpi, tpcpi, float);
DECLARE_SOA_COLUMN(TPCka, tpcka, float);
DECLARE_SOA_COLUMN(TPCpr, tpcpr, float);
DECLARE_SOA_COLUMN(TPCel, tpcel, float);
DECLARE_SOA_COLUMN(TOFpi, tofpi, float);
DECLARE_SOA_COLUMN(TOFka, tofka, float);
DECLARE_SOA_COLUMN(TOFpr, tofpr, float);
DECLARE_SOA_COLUMN(TOFel, tofel, float);
DECLARE_SOA_COLUMN(TPCmu, tpcmu, float);
DECLARE_SOA_COLUMN(TOFmu, tofmu, float);
DECLARE_SOA_COLUMN(TPCde, tpcde, float);
DECLARE_SOA_COLUMN(TOFde, tofde, float);
} // namespace sgtrack
DECLARE_SOA_TABLE(SGTracks, "AOD", "SGTRACK",
                  o2::soa::Index<>, sgtrack::SGEventId,
                  sgtrack::Pt, sgtrack::Eta, sgtrack::Phi, sgtrack::Sign, sgtrack::TPCpi, sgtrack::TPCka, sgtrack::TPCpr, sgtrack::TPCel, sgtrack::TOFpi, sgtrack::TOFka, sgtrack::TOFpr, sgtrack::TOFel, sgtrack::TPCmu, sgtrack::TOFmu, sgtrack::TPCde, sgtrack::TOFde);
using SGTrack = SGTracks::iterator;
} // namespace o2::aod

#endif // PWGUD_DATAMODEL_SGTABLES_H_
