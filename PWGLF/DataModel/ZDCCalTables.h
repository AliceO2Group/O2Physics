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

#include <Framework/AnalysisDataModel.h>

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

// Extra optional linked table.
// This table does NOT duplicate cent, vx, vy, vz, run number, trigger, etc.
// It only stores the ZDC energies and links back to ZDCCalTables.
namespace zdcenergytable
{
DECLARE_SOA_INDEX_COLUMN(ZDCCalTable, zdcCalTable);

DECLARE_SOA_COLUMN(ZNACommon, znaCommon, float);
DECLARE_SOA_COLUMN(ZNCCommon, zncCommon, float);

DECLARE_SOA_COLUMN(ZNA0, zna0, float);
DECLARE_SOA_COLUMN(ZNA1, zna1, float);
DECLARE_SOA_COLUMN(ZNA2, zna2, float);
DECLARE_SOA_COLUMN(ZNA3, zna3, float);

DECLARE_SOA_COLUMN(ZNC0, znc0, float);
DECLARE_SOA_COLUMN(ZNC1, znc1, float);
DECLARE_SOA_COLUMN(ZNC2, znc2, float);
DECLARE_SOA_COLUMN(ZNC3, znc3, float);
} // namespace zdcenergytable

DECLARE_SOA_TABLE(ZDCEnergyTables, "AOD", "ZDCENERGY",
                  zdcenergytable::ZDCCalTableId,
                  zdcenergytable::ZNACommon,
                  zdcenergytable::ZNCCommon,
                  zdcenergytable::ZNA0,
                  zdcenergytable::ZNA1,
                  zdcenergytable::ZNA2,
                  zdcenergytable::ZNA3,
                  zdcenergytable::ZNC0,
                  zdcenergytable::ZNC1,
                  zdcenergytable::ZNC2,
                  zdcenergytable::ZNC3);

using ZDCEnergyTable = ZDCEnergyTables::iterator;

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_ZDCCALTABLES_H_
