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

/// \file GloCCDBObjects.h
/// \brief Declarative CCDB columns for the most-commonly needed GLO/ objects.
///
/// Provides a single shared timestamped table (`aod::GloCCDBObjects`) that
/// tasks can join with `aod::BCsWithTimestamps` to obtain the four GLO/
/// calibration objects without using `Service<BasicCCDBManager>`.
///
/// Usage:
/// \code
///   #include "Common/DataModel/GloCCDBObjects.h"
///   // ...
///   using BCsWithCCDB = soa::Join<aod::BCsWithTimestamps, aod::GloCCDBObjects>;
///   void process(BCsWithCCDB const& bcs) {
///     const auto& grpmag  = bcs.begin().grpMagField();
///     const auto& meanvtx = bcs.begin().meanVertex();
///     const auto& grpecs  = bcs.begin().grpECS();
///     const auto& grplhcif= bcs.begin().grpLHCIF();
///   }
/// \endcode
///
/// If you need only a subset of the four objects, declare your own
/// `DECLARE_SOA_TIMESTAMPED_TABLE` with the relevant subset of columns from
/// the `o2::aod::ccdbGlo` namespace rather than joining `aod::GloCCDBObjects`.
///
/// Note: MatLayerCylSet is intentionally omitted — it requires
/// `MatLayerCylSet::rectifyPtrFromFile()` after deserialisation, which the
/// CCDB column mechanism does not perform.

#ifndef COMMON_DATAMODEL_GLOCCDBOBJECTS_H_
#define COMMON_DATAMODEL_GLOCCDBOBJECTS_H_

#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPECSObject.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace ccdbGlo
{
DECLARE_SOA_CCDB_COLUMN(GRPMagField, grpMagField, o2::parameters::GRPMagField, "GLO/Config/GRPMagField");   //!
DECLARE_SOA_CCDB_COLUMN(MeanVertex, meanVertex, o2::dataformats::MeanVertexObject, "GLO/Calib/MeanVertex"); //!
DECLARE_SOA_CCDB_COLUMN(GRPECSObject, grpECS, o2::parameters::GRPECSObject, "GLO/Config/GRPECS");           //!
DECLARE_SOA_CCDB_COLUMN(GRPLHCIFData, grpLHCIF, o2::parameters::GRPLHCIFData, "GLO/Config/GRPLHCIF");       //!
} // namespace ccdbGlo

/// Full table — join with aod::BCsWithTimestamps to obtain all four objects.
DECLARE_SOA_TIMESTAMPED_TABLE(GloCCDBObjects, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "GLOCCDBOBJ", //!
                              ccdbGlo::GRPMagField, ccdbGlo::MeanVertex, ccdbGlo::GRPECSObject, ccdbGlo::GRPLHCIFData);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_GLOCCDBOBJECTS_H_
