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

/// \file ZDCCalTables.h
///
/// author: prottay das 07/09/2024
/// email: prottay.das@cern.ch

#ifndef PWGLF_DATAMODEL_ZDCCALTABLES_H_
#define PWGLF_DATAMODEL_ZDCCALTABLES_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"

#include <cmath>

namespace o2::aod
{
namespace zdccaltable
{
DECLARE_SOA_COLUMN(TriggerEventZDC, triggerEventZDC, bool);
DECLARE_SOA_COLUMN(TriggerEventRunNo, triggerEventRunNo, int);
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Vx, vx, float);
DECLARE_SOA_COLUMN(Vy, vy, float);
DECLARE_SOA_COLUMN(Vz, vz, float);
DECLARE_SOA_COLUMN(ZnaC, znaC, float);
DECLARE_SOA_COLUMN(ZncC, zncC, float);
DECLARE_SOA_COLUMN(ZnaE0, znaE0, float);
DECLARE_SOA_COLUMN(ZnaE1, znaE1, float);
DECLARE_SOA_COLUMN(ZnaE2, znaE2, float);
DECLARE_SOA_COLUMN(ZnaE3, znaE3, float);
DECLARE_SOA_COLUMN(ZncE0, zncE0, float);
DECLARE_SOA_COLUMN(ZncE1, zncE1, float);
DECLARE_SOA_COLUMN(ZncE2, zncE2, float);
DECLARE_SOA_COLUMN(ZncE3, zncE3, float);
} // namespace zdccaltable
DECLARE_SOA_TABLE(ZDCCalTables, "AOD", "ZDCCALTABLE",
                  zdccaltable::TriggerEventZDC,
                  zdccaltable::TriggerEventRunNo,
                  zdccaltable::Cent,
                  zdccaltable::Vx,
                  zdccaltable::Vy,
                  zdccaltable::Vz,
                  zdccaltable::ZnaC,
                  zdccaltable::ZncC,
                  zdccaltable::ZnaE0,
                  zdccaltable::ZnaE1,
                  zdccaltable::ZnaE2,
                  zdccaltable::ZnaE3,
                  zdccaltable::ZncE0,
                  zdccaltable::ZncE1,
                  zdccaltable::ZncE2,
                  zdccaltable::ZncE3);
using ZDCCalTable = ZDCCalTables::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_ZDCCALTABLES_H_
