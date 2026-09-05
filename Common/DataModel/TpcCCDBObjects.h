// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TpcCCDBObjects.h
/// \brief Declarative CCDB columns for TPC calibration objects.
///
/// Unlike the geometry/material family in GloCCDBObjects.h, the drift velocity
/// genuinely varies within a run, so the table keeps the default uniformity — one
/// object per distinct timestamp — rather than collapsing per run.
///
/// Usage:
/// \code
///   using BCsWithVDrift = soa::Join<aod::BCsWithTimestamps, aod::TpcCalibCCDBObjects>;
///   vdriftManager.update(bc.vdriftTgl());
/// \endcode

#ifndef COMMON_DATAMODEL_TPCCCDBOBJECTS_H_
#define COMMON_DATAMODEL_TPCCCDBOBJECTS_H_

#include <DataFormatsTPC/VDriftCorrFact.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace ccdbTpc
{
DECLARE_SOA_CCDB_COLUMN(VDriftTgl, vdriftTgl, o2::tpc::VDriftCorrFact, "TPC/Calib/VDriftTgl"); //!
} // namespace ccdbTpc

DECLARE_SOA_TIMESTAMPED_TABLE(TpcCalibCCDBObjects, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "TPCCALIBCCDB", //!
                              ccdbTpc::VDriftTgl);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_TPCCCDBOBJECTS_H_
