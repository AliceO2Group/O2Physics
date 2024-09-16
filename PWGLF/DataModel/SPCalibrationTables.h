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

/// \file SPCalibrationTables.h
///
/// author: prottay das 07/09/2024
/// email: prottay.das@cern.ch

#ifndef PWGLF_DATAMODEL_SPCALIBRATIONTABLES_H_
#define PWGLF_DATAMODEL_SPCALIBRATIONTABLES_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace spcalibrationtable
{
DECLARE_SOA_COLUMN(TriggerEvent, triggerevent, bool);
DECLARE_SOA_COLUMN(TriggerEventRunNo, triggereventrunno, int);
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Vz, vz, float);
DECLARE_SOA_COLUMN(ZNAEN1, znaen1, float);
DECLARE_SOA_COLUMN(ZNAEN2, znaen2, float);
DECLARE_SOA_COLUMN(ZNAEN3, znaen3, float);
DECLARE_SOA_COLUMN(ZNAEN4, znaen4, float);
DECLARE_SOA_COLUMN(ZNCEN1, zncen1, float);
DECLARE_SOA_COLUMN(ZNCEN2, zncen2, float);
DECLARE_SOA_COLUMN(ZNCEN3, zncen3, float);
DECLARE_SOA_COLUMN(ZNCEN4, zncen4, float);
DECLARE_SOA_COLUMN(QXZDCA, qxZDCA, float);
DECLARE_SOA_COLUMN(QXZDCC, qxZDCC, float);
DECLARE_SOA_COLUMN(QYZDCA, qyZDCA, float);
DECLARE_SOA_COLUMN(QYZDCC, qyZDCC, float);
DECLARE_SOA_COLUMN(PsiZDCC, psiZDCC, float);
DECLARE_SOA_COLUMN(PsiZDCA, psiZDCA, float);
} // namespace spcalibrationtable
DECLARE_SOA_TABLE(SPCalibrationTables, "AOD", "SPCALCOLS",
                  spcalibrationtable::TriggerEvent,
                  spcalibrationtable::TriggerEventRunNo,
                  spcalibrationtable::Cent,
                  spcalibrationtable::Vz,
                  spcalibrationtable::ZNAEN1,
                  spcalibrationtable::ZNAEN2,
                  spcalibrationtable::ZNAEN3,
                  spcalibrationtable::ZNAEN4,
                  spcalibrationtable::ZNCEN1,
                  spcalibrationtable::ZNCEN2,
                  spcalibrationtable::ZNCEN3,
                  spcalibrationtable::ZNCEN4,
                  spcalibrationtable::QXZDCA,
                  spcalibrationtable::QXZDCC,
                  spcalibrationtable::QYZDCA,
                  spcalibrationtable::QYZDCC,
                  spcalibrationtable::PsiZDCC,
                  spcalibrationtable::PsiZDCA);
using SPCalibrationTable = SPCalibrationTables::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_SPCALIBRATIONTABLES_H_
