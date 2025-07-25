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

/// \file SPTableZDC.h

#ifndef PWGCF_DATAMODEL_SPTABLEZDC_H_
#define PWGCF_DATAMODEL_SPTABLEZDC_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace sptablezdc
{
DECLARE_SOA_COLUMN(Runnumber, runnumber, int);
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Vx, vx, float);
DECLARE_SOA_COLUMN(Vy, vy, float);
DECLARE_SOA_COLUMN(Vz, vz, float);
DECLARE_SOA_COLUMN(QXA, qxA, float);
DECLARE_SOA_COLUMN(QYA, qyA, float);
DECLARE_SOA_COLUMN(QXC, qxC, float);
DECLARE_SOA_COLUMN(QYC, qyC, float);
DECLARE_SOA_COLUMN(IsSelected, isSelected, bool);
DECLARE_SOA_COLUMN(Iteration, iteration, int);
DECLARE_SOA_COLUMN(Step, step, int);

} // namespace sptablezdc

DECLARE_SOA_TABLE(SPTableZDC, "AOD", "SPZDC",
                  sptablezdc::Runnumber,
                  sptablezdc::Cent,
                  sptablezdc::Vx,
                  sptablezdc::Vy,
                  sptablezdc::Vz,
                  sptablezdc::QXA,
                  sptablezdc::QYA,
                  sptablezdc::QXC,
                  sptablezdc::QYC,
                  sptablezdc::IsSelected,
                  sptablezdc::Iteration,
                  sptablezdc::Step);
} // namespace o2::aod
#endif // PWGCF_DATAMODEL_SPTABLEZDC_H_
