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
DECLARE_SOA_COLUMN(QxA, qxA, float);
DECLARE_SOA_COLUMN(QxC, qxC, float);
DECLARE_SOA_COLUMN(QyA, qyA, float);
DECLARE_SOA_COLUMN(QyC, qyC, float);
} // namespace zdccaltable
DECLARE_SOA_TABLE(ZDCCalTables, "AOD", "ZDCCALTABLE",
                  zdccaltable::TriggerEventZDC,
                  zdccaltable::TriggerEventRunNo,
                  zdccaltable::Cent,
                  zdccaltable::Vx,
                  zdccaltable::Vy,
                  zdccaltable::Vz,
                  zdccaltable::QxA,
                  zdccaltable::QxC,
                  zdccaltable::QyA,
                  zdccaltable::QyC);
using ZDCCalTable = ZDCCalTables::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_ZDCCALTABLES_H_
