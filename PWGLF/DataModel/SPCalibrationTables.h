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

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"

#include <cmath>

namespace o2::aod
{
namespace spcalibrationtable
{
DECLARE_SOA_COLUMN(TriggerEventSP, triggereventsp, bool);
DECLARE_SOA_COLUMN(TriggerEventRunNo, triggereventrunno, int);
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Vx, vx, float);
DECLARE_SOA_COLUMN(Vy, vy, float);
DECLARE_SOA_COLUMN(Vz, vz, float);
DECLARE_SOA_COLUMN(ZNAC, znaC, float);
DECLARE_SOA_COLUMN(ZNCC, zncC, float);
DECLARE_SOA_COLUMN(ZNAE0, znaE0, float);
DECLARE_SOA_COLUMN(ZNAE1, znaE1, float);
DECLARE_SOA_COLUMN(ZNAE2, znaE2, float);
DECLARE_SOA_COLUMN(ZNAE3, znaE3, float);
DECLARE_SOA_COLUMN(ZNCE0, zncE0, float);
DECLARE_SOA_COLUMN(ZNCE1, zncE1, float);
DECLARE_SOA_COLUMN(ZNCE2, zncE2, float);
DECLARE_SOA_COLUMN(ZNCE3, zncE3, float);
DECLARE_SOA_COLUMN(QXZDCA, qxZDCA, float);
DECLARE_SOA_COLUMN(QXZDCC, qxZDCC, float);
DECLARE_SOA_COLUMN(QYZDCA, qyZDCA, float);
DECLARE_SOA_COLUMN(QYZDCC, qyZDCC, float);
DECLARE_SOA_COLUMN(PsiZDCC, psiZDCC, float);
DECLARE_SOA_COLUMN(PsiZDCA, psiZDCA, float);
} // namespace spcalibrationtable
DECLARE_SOA_TABLE(SPCalibrationTables, "AOD", "SPCALCOLS",
                  spcalibrationtable::TriggerEventSP,
                  spcalibrationtable::TriggerEventRunNo,
                  spcalibrationtable::Cent,
                  spcalibrationtable::Vx,
                  spcalibrationtable::Vy,
                  spcalibrationtable::Vz,
                  spcalibrationtable::ZNAC,
                  spcalibrationtable::ZNCC,
                  spcalibrationtable::ZNAE0,
                  spcalibrationtable::ZNAE1,
                  spcalibrationtable::ZNAE2,
                  spcalibrationtable::ZNAE3,
                  spcalibrationtable::ZNCE0,
                  spcalibrationtable::ZNCE1,
                  spcalibrationtable::ZNCE2,
                  spcalibrationtable::ZNCE3,
                  spcalibrationtable::QXZDCA,
                  spcalibrationtable::QXZDCC,
                  spcalibrationtable::QYZDCA,
                  spcalibrationtable::QYZDCC,
                  spcalibrationtable::PsiZDCC,
                  spcalibrationtable::PsiZDCA);
using SPCalibrationTable = SPCalibrationTables::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_SPCALIBRATIONTABLES_H_
