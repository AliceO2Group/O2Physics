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

/// \file SPCalibrationTablesFT0C.h
///
/// author: prottay das 11/03/2025
/// email: prottay.das@cern.ch

#ifndef PWGLF_DATAMODEL_SPCALIBRATIONTABLESFT0C_H_
#define PWGLF_DATAMODEL_SPCALIBRATIONTABLESFT0C_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace spcalibrationtableFT0C
{
DECLARE_SOA_COLUMN(TriggerEventSPFT0C, triggereventspft0c, bool);
DECLARE_SOA_COLUMN(TriggerEventRunNoFT0C, triggereventrunnoft0c, int);
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(QXFT0C, qxFT0C, float);
DECLARE_SOA_COLUMN(QYFT0C, qyFT0C, float);
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);
} // namespace spcalibrationtableFT0C
DECLARE_SOA_TABLE(SPCalibrationTablesFT0C, "AOD", "SPCAL",
                  spcalibrationtableFT0C::TriggerEventSPFT0C,
                  spcalibrationtableFT0C::TriggerEventRunNoFT0C,
                  spcalibrationtableFT0C::Cent,
                  spcalibrationtableFT0C::QXFT0C,
                  spcalibrationtableFT0C::QYFT0C,
                  spcalibrationtableFT0C::PsiFT0C);
using SPCalibrationTableFT0C = SPCalibrationTablesFT0C::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_SPCALIBRATIONTABLESFT0C_H_
